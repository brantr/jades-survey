.. _jwst_image_stack:

=========================================
``make_jwst_stack.py``: JWST Image Stacking and Weight Normalisation
=========================================

This document provides a detailed walk‑through of the reference implementation
``make_stack_jwst.py`` (the script shown below).  The script is designed to combine
multiple *JWST* (James Webb Space Telescope) exposures – typically the *SCI*
(science), *ERR* (error), *WHT* (weight) and optional *NIM* (mask) extensions –
into a single, well‑behaved FITS file.  It also offers a number of command line
options for handling instrument‑specific quirks (e.g. MIRI vs. NIRCam) and for
pre‑processing the *NIM* layer to mask diffraction spikes.

The documentation is written in reStructuredText (RST) format so it can be
hosted directly on a ReadTheDocs site.

--------------------------------------------------------------------
Table of contents
--------------------------------------------------------------------

.. contents::
   :local:
   :depth: 2

--------------------------------------------------------------------
1.  High‑level purpose
--------------------------------------------------------------------

The script performs three high‑level tasks:

1. **Parse a list of input FITS files** (one per line) supplied by the user.
2. **Optionally pre‑process the NIM mask** to identify diffraction spikes that
   should be ignored during stacking.
3. **Iteratively accumulate** the science data, its variance, the weight map,
   exposure time and (optionally) the NIM mask, applying instrument‑specific
   corrections such as inverse‑variance weighting for NIRCam data.

The result is a single FITS file containing the combined ``SCI``, ``ERR``,
``EXP``, ``WHT`` and (if requested) ``NIM`` extensions, together with the
original primary header.

--------------------------------------------------------------------
2.  Dependencies
--------------------------------------------------------------------

The script relies on the following third‑party Python packages:

- :mod:`numpy` – array handling
- :mod:`astropy.io.fits` – FITS I/O
- :mod:`scipy.ndimage` – binary morphology (`binary_dilation`,
  `binary_erosion`)
- :mod:`skimage.measure` – connected‑component labeling
- :mod:`skimage.segmentation` – label expansion
- :mod:`photutils.segmentation` – source catalog generation
- :mod:`sep` – (imported but not used in the current version)

These packages are all available from *PyPI* and can be installed with:

.. code-block:: bash

   pip install numpy astropy scipy scikit-image photutils sep

--------------------------------------------------------------------
3.  Command‑line interface
--------------------------------------------------------------------

The entry point is the :func:`create_parser` function which builds an
:class:`argparse.ArgumentParser`.  The supported options are summarised in the
table below.

+----------------------+----------------------+--------------------------------------+
| Short / Long flag    | Type / Action        | Description                          |
+======================+======================+======================================+
| ``-i`` / ``--input_list`` | ``str`` (default: ``input_list.txt``) | Text file containing one FITS path per line. |
+----------------------+----------------------+--------------------------------------+
| ``-o`` / ``--output`` | ``str`` (default: ``out.fits``) | Name of the stacked output file. |
+----------------------+----------------------+--------------------------------------+
| ``--renorm-err``     | ``float`` (optional) | Divide the final error image by this factor (useful for MIRI). |
+----------------------+----------------------+--------------------------------------+
| ``-e`` / ``--exp``   | ``store_true`` (default: ``False``) | Use the ``WHT`` extension as a fake ``EXP`` when set. |
+----------------------+----------------------+--------------------------------------+
| ``-m`` / ``--multimodule`` | ``store_true`` (default: ``False``) | Force the ``MODULE`` keyword in all HDUs to the value ``MULTIPLE``. |
+----------------------+----------------------+--------------------------------------+
| ``--preprocess-nim``| ``store_false`` (default: ``True``) | Disable the NIM pre‑processing step. |
+----------------------+----------------------+--------------------------------------+
| ``--nim``            | ``store_false`` (default: ``True``) | Do **not** propagate the NIM extension to the output. |
+----------------------+----------------------+--------------------------------------+
| ``--inv-var``        | ``store_false`` (default: ``True``) | Disable the new inverse‑variance weighting scheme. |
+----------------------+----------------------+--------------------------------------+
| ``-v`` / ``--verbose`` | ``store_true`` (default: ``False``) | Print progress information. |
+----------------------+----------------------+--------------------------------------+
| ``-b`` / ``--bit``   | ``store_true`` (default: ``False``) | Cast floating‑point images to 32‑bit; NIM to 16‑bit. |
+----------------------+----------------------+--------------------------------------+

The parser is invoked inside :func:`main` with ``parser.parse_args()`` and the
resulting namespace is stored in ``args``.

--------------------------------------------------------------------
4.  Core functions
--------------------------------------------------------------------

The script is deliberately short, with most logic residing in three helper
functions and the ``main`` routine.

--------------------------------------------------------------------
4.1 ``create_parser()``
--------------------------------------------------------------------

*Purpose*: Build and return an :class:`argparse.ArgumentParser` configured with
the options described above.

*Key points*:

- ``dest`` arguments map the command line flag to an attribute of the
  ``args`` namespace.
- ``action='store_true'`` / ``store_false`` toggle Boolean flags.
- ``default`` values provide sensible fall‑backs when the user omits a flag.

--------------------------------------------------------------------
4.2 ``preprocess_nim(fl_input_list)``
--------------------------------------------------------------------

*Purpose*: Scan all input files and create a Boolean mask (``nim_flag``) that
identifies pixels which **should not** be masked by the NIM extension during
the final stack.

*Algorithm*:

1. Allocate three temporary arrays of the same shape as the first ``NIM`` layer:
   ``nim_plus`` (accumulates good data), ``nim_two`` (counts diffraction‑spike
   pixels), and ``nim_flag`` (the final Boolean mask).
2. Loop over every file in ``fl_input_list``:
   - Load the ``NIM`` and ``WHT`` extensions.
   - Increment ``nim_plus`` where both ``NIM`` and ``WHT`` are positive
     (valid data).
   - Increment ``nim_two`` where ``NIM == -2`` (diffraction spike) **and**
     ``WHT`` is positive.
3. After the loop, set ``nim_flag`` to ``True`` for pixels where there is **no**
   good data (``nim_plus == 0``) but at least one spike pixel was found
   (``nim_two > 0``).

*Result*: ``nim_flag`` is a ``bool`` array, the same size as the input images,
where ``True`` marks regions that should be ignored when applying the NIM mask.

--------------------------------------------------------------------
4.3 ``get_instrument(header_img)``
--------------------------------------------------------------------

*Purpose*: Determine which instrument (``HST``, ``NIRCam`` or ``MIRI``) produced
the input data, based on header keywords.

*Logic*:

- Prefer the ``TELESCOP`` keyword.  If it equals ``'HST'`` the function returns
  ``'HST'``.
- Otherwise look for ``INSTRUME`` – the value is returned unchanged.
- If neither keyword is present, fall back to ``'MIRI'`` (the most common case
  for JWST NIRCam/MIRI data).

The function also prints a diagnostic message when ``verbose`` is enabled.

--------------------------------------------------------------------
4.4 ``main()``
--------------------------------------------------------------------

The ``main`` function orchestrates the complete workflow:

1. **Parse arguments** and optionally print them (``--verbose``).
2. **Read the input list** from ``args.input_list``; each line is stripped of
   trailing new‑lines.
3. **Pre‑process NIM** if ``args.prenim`` is ``True`` – the Boolean mask
   ``nim_flag`` is obtained from :func:`preprocess_nim`.
4. **Iterate over all input files** and accumulate data:

   - On the **first file** the primary header and all relevant extension
     headers are cached.  Empty output arrays are allocated:
     ``sci_out``, ``err_out``, ``wht_out``, ``exp_out`` and optionally
     ``nim_out``.
   - For each file the script loads the ``WHT`` extension, and (if the data
     are from NIRCam) also the ``NIM`` extension.
   - **NIM handling**:
     * If ``args.nim`` is ``True`` and the instrument is NIRCam, the script
       determines which pixels to mask based on ``nim_flag``.
     * When ``args.prenim`` is active the code performs a sophisticated
       spike‑detection routine (see Section 5.2).  Pixels belonging to a
       diffraction spike that are *partially* covered by good data are
       completely masked (``nim == -4``).
     * The weight of all masked pixels is set to zero.
   - **Inverse‑variance weighting** (enabled with ``--inv-var``):
     * If the primary header contains the keyword ``VAR_P_BG`` (background
       Poisson variance), the script adds this variance term to the weight
       map in quadrature:

       ``wht = 1 / (1/wht + var_p_bkg)``

   - **Accumulate** the extensions:
     * ``wht_out``  ← ``wht_out`` + ``wht``
     * ``sci_out``  ← ``sci_out`` + ``wht * SCI``
     * ``err_out``  ← ``err_out`` + ``(wht * ERR)^2``
     * ``exp_out``  ← ``exp_out`` + ``EXP`` (or ``WHT`` if ``--exp`` is set)
     * ``nim_out``  ← ``nim_out`` + ``NIM`` (only for NIRCam when ``args.nim``)
5. **Normalise** the stacked images:
   - For all pixels where ``wht_out > 0``:
     * ``SCI = SCI / wht_out``
     * ``ERR = sqrt( ERR / wht_out**2 )``
   - If ``--exp`` was requested, ``exp_out`` is recomputed from the primary
     header keyword ``EFFEXPTM``.
6. **Optional post‑processing**:
   - Cast to 32‑bit floats if ``--bit`` is set.
   - Apply the user‑supplied ``--renorm-err`` scaling to the error image.
   - Overwrite the ``MODULE`` keyword with ``MULTIPLE`` if ``--multimodule``.
7. **Write the output**:
   - Create a new :class:`astropy.io.fits.HDUList` consisting of the primary
     HDU and the stacked extensions.
   - Write the file to ``args.output`` with ``overwrite=True``.

The script ends with the conventional ``if __name__ == "__main__":`` guard
that calls ``main()``.

--------------------------------------------------------------------
5.  Detailed processing steps
--------------------------------------------------------------------

Below we expand on two of the more intricate parts of the pipeline.

--------------------------------------------------------------------
5.1  NIM pre‑processing (``preprocess_nim``)
--------------------------------------------------------------------

The NIM extension encodes per‑pixel data quality flags used by JWST
pipelines:

- ``>0`` – good data
- ``-2`` – diffraction spike (masked)
- ``-4`` – saturated or otherwise defective pixel

The pre‑processing routine determines *where* spikes lie **without** any
overlapping good data from other exposures.  Such isolated spikes are then
treated as *masked* for the entire stack (they would otherwise introduce
artificially low weight).

The algorithm is a simple two‑pass accumulator:

1. **First pass** – count the number of good pixels (`nim_plus`) and the
   number of spike pixels (`nim_two`) for each location across *all* images.
2. **Second pass** – any location where ``nim_plus == 0`` but ``nim_two > 0``
   is flagged as ``True`` in ``nim_flag``.

This Boolean mask is later used to decide whether a pixel belonging to a
diffraction spike should be kept or zeroed out.

--------------------------------------------------------------------
5.2  Spike detection and masking (inside the main loop)
--------------------------------------------------------------------

When ``args.prenim`` is active **and** the data are from NIRCam, the script
executes a more refined spike handling routine:

1. **Identify spike pixels** (`nim == -2`).
2. **Create a binary image** where spike pixels are 1 and everything else is 0.
3. **Erode** the binary image by 10 pixels (``binary_erosion``).  This removes
   thin connections that could cause distinct spikes to be merged.
4. **Label** the eroded spikes using ``skimage.measure.label``.
5. **Expand** the labelled regions by 11 pixels (``expand_labels``) to restore
   the original spike size while keeping them separated.
6. **Generate a source catalog** with ``photutils.SourceCatalog`` on the
   original weight image, using the labelled spikes as segmentation.
7. **Loop over each labelled spike**:
   - Extract the pixel indices belonging to the spike.
   - Check the pre‑computed ``nim_flag`` for those indices.
   - If **any** pixel of the spike is *not* flagged (i.e. underlying good data
     exists), the entire spike is set to ``nim == -4`` – effectively masking it
     completely.

The routine also writes two diagnostic FITS files for each exposure:

- ``<filename>.spikes.fits`` – the binary spike mask before erosion.
- ``<filename>.spike_labels.fits`` – the labelled spike map after expansion.

These files are useful for debugging and for visual inspection of the
masking process.

--------------------------------------------------------------------
6.  Example usage
--------------------------------------------------------------------

Assume you have a text file ``my_images.txt`` containing the absolute paths of
five NIRCam mosaics you wish to combine.

.. code-block:: bash

   # Basic stacking, keep NIM, verbose output
   python stack_jwst.py -i my_images.txt -o stacked.fits -v

   # Disable NIM propagation and force 32‑bit output
   python stack_jwst.py -i my_images.txt -o stacked.fits --nim -b

   # Use the inverse‑variance weighting scheme and renormalise the error
   python stack_jwst.py -i my_images.txt -o stacked.fits --renorm-err 1.5

   # Treat the WHT extension as a fake EXP map (useful for certain pipelines)
   python stack_jwst.py -i my_images.txt -o stacked.fits --exp

   # Override the MODULE keyword for all HDUs (required by some downstream
   # processing tools)
   python stack_jwst.py -i my_images.txt -o stacked.fits -m

--------------------------------------------------------------------
7.  Configuration notes
--------------------------------------------------------------------

- **Input list format** – one absolute or relative path per line, no
  surrounding whitespace.  Blank lines are ignored because they are stripped
  before processing.
- **Instrument detection** – if a file lacks both ``TELESCOP`` and ``INSTRUME``,
  the script assumes ``MIRI``.  This behaviour can be customised by editing
  :func:`get_instrument`.
- **Weight handling** – the script expects the ``WHT`` extension to contain the
  *inverse variance* (i.e. ``1/σ²``) of the science data.  When ``--inv-var`` is
  enabled, a background Poisson variance term (``VAR_P_BG``) is added in
  quadrature.
- **Memory usage** – the script loads each image *twice* (once for the NIM
  pre‑processing, once for the main accumulation).  For very large mosaics
  you may need to increase available RAM or modify the code to stream the
  data in smaller blocks.

--------------------------------------------------------------------
8.  Extending the script
--------------------------------------------------------------------

The modular structure makes it straightforward to add new features:

* **Additional instruments** – extend ``get_instrument`` and add any
  instrument‑specific weighting or mask handling in the main loop.
* **Alternative weighting schemes** – replace the block under
  ``#new weighting scheme`` with a custom implementation.
* **Parallel processing** – the per‑file accumulation can be parallelised
  with ``concurrent.futures`` or ``dask`` if the I/O bandwidth permits.
* **More diagnostics** – write out intermediate stacked products (e.g.
  ``sci_partial.fits``) for each iteration to track convergence.

--------------------------------------------------------------------
9.  Full source listing
--------------------------------------------------------------------

For reference, the complete script is reproduced below with inline comments
highlighting the major sections.

.. literalinclude:: stack_jwst.py
   :language: python
   :linenos:

--------------------------------------------------------------------
10.  License
--------------------------------------------------------------------

The original author has not supplied a license header.  If you intend to
redistribute or modify the script, please contact the maintainer or apply an
appropriate open‑source license (e.g. MIT, BSD) consistent with the
dependencies used.

--------------------------------------------------------------------
11.  Bibliography
--------------------------------------------------------------------

* JWST Data Management System – Calibration Reference Data System.
* Astropy Collaboration, *Astropy: A community Python package for astronomy*,
  2013, *A&A*, 558, A33.
* van der Walt et al., *scikit-image: Image processing in Python*, 2014,
  *PeerJ*, 2:e453.
* Bradley et al., *Photutils: Photometry tools for astronomy*, 2022,
  *Astropy* package.
