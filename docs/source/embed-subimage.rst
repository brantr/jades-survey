.. _fits_mosaic_script:

=========================================================
``embed_subimage.py`` Mosaic a Sub‑Image into a Large FITS Image (Python Script)
=========================================================

This document explains the purpose, design and operation of the
``mosaic_fits.py`` script.  The script is a small command‑line utility that
takes a *sub‑image* (a cut‑out of a larger astronomical image stored in a
FITS file) and inserts it into a *full‑size* image – either an existing
mosaic or a newly created empty image defined by a header file.

The script is written for **Python 3**, relies on :mod:`numpy` and
:mod:`astropy.io.fits`, and is intended to be run from the command line.
Below you will find a line‑by‑line description of the code, the
command‑line interface, the core functions and the overall workflow.

-----------------------------------------------------------------
Table of contents
-----------------------------------------------------------------

.. contents::
   :depth: 2
   :local:

-----------------------------------------------------------------
1.  High‑level purpose
-----------------------------------------------------------------

The script performs three related tasks:

1. **Read the headers** of a *sub‑image* and a *full‑size* image (or a
   header file that describes the full‑size image).

2. **Determine the pixel region** (slice) of the full‑size image that
   corresponds to the sub‑image using World Coordinate System (WCS)
   information (``CRPIX`` and ``CRVAL`` keywords).

3. **Create or open a full‑size FITS file** and copy the data from the
   sub‑image into the appropriate region, optionally adding the data
   (i.e. summing) rather than overwriting it.

The script can be used to build large mosaics from many small cut‑outs,
or to insert a new cut‑out into an existing mosaic while preserving
existing data.

-----------------------------------------------------------------
2.  Dependencies
-----------------------------------------------------------------

* ``numpy`` – for numerical arrays and simple arithmetic.
* ``astropy.io.fits`` – for reading and writing FITS files and handling
  headers.
* Standard library modules: :mod:`os`, :mod:`argparse`.

-----------------------------------------------------------------
3.  Command‑line interface
-----------------------------------------------------------------

The script uses :class:`argparse.ArgumentParser` to define the following
options:

.. code-block:: bash

   --sub_image   PATH   (default: small_image.fits)
       Full path to the FITS file that contains the sub‑image.

   --full_image  PATH   (default: giant_image.fits)
       Destination file for the full‑size mosaic.  If ``--full_header`` is
       supplied, this file **must not already exist** (the script never
       overwrites an existing file for safety).

   --sub_header  PATH   (default: None)
       Optional header file for the sub‑image.  If omitted, the header
       is taken directly from the ``SCI`` extension of ``--sub_image``.

   --full_header PATH   (default: hlf_v2.0.1_30mas_cropped.hdr)
       Header file describing the full‑size image.  If omitted the
       script expects an existing FITS file whose ``SCI`` extension
       already contains the required header.

   --sci                (store_true)
       Only process the ``SCI`` extension.  By default all scientific
       extensions are processed.

   --nim                (store_false, default=True)
       Propagate the ``NIM`` extension.  ``--nim`` flips the default to
       *False*.

   --telescope  STRING
       Manually set the ``TELESCOP`` keyword in the output header.

   --instrument STRING
       Manually set the ``INSTRUME`` keyword in the output header.

The script is executed as:

.. code-block:: bash

   python mosaic_fits.py \
       --sub_image  my_cutout.fits \
       --full_image my_mosaic.fits \
       --full_header full_header.hdr \
       --telescope JWST \
       --instrument NIRCam

-----------------------------------------------------------------
4.  Core functions
-----------------------------------------------------------------

The implementation is split into three reusable functions and a ``main``
section that orchestrates the workflow.

-----------------------------------------------------------------
4.1 ``find_slices(sub_header, full_header)``

``find_slices`` determines the integer pixel offsets that map the sub‑image
onto the full‑size image.  It works with *WCS* information stored in the
FITS headers.

**Signature**

.. code-block:: python

   def find_slices(sub_header, full_header):
       """
       Parameters
       ----------
       sub_header : fits.Header or dict
           Header of the sub‑image (must contain CD, CRVAL, CRPIX, NAXIS*).

       full_header : fits.Header or dict
           Header of the full‑size image.

       Returns
       -------
       xstart, xend, ystart, yend, flag_transpose
           Integer start/end indices for the X and Y axes (Python uses
           ``[y, x]`` order for 2‑D arrays) and a boolean indicating whether
           the image must be transposed before insertion.
       """
       ...

**What it does**

1. **Sanity checks** – prints a few WCS keywords (``PC*`` and ``CD*``)
   from both headers for debugging; asserts that the scale and rotation
   keywords (``CD1_1``, ``CD2_2``, ``CD1_2``, ``CD2_1``) match within a
   relative tolerance of ``1e‑6``.
2. **Size of the sub‑image** – reads ``NAXIS1`` (X size) and ``NAXIS2`` (Y
   size) from ``sub_header``.
3. **Pixel offset** – computes the offset between reference pixels
   (``CRPIX``) of the two images:

   .. math::

      \text{xstart} = \text{CRPIX1}_\text{full} - \text{CRPIX1}_\text{sub}
      \\
      \text{ystart} = \text{CRPIX2}_\text{full} - \text{CRPIX2}_\text{sub}

   The offsets must be integer values; the function asserts that the
   modulus with 1.0 is zero.
4. **End indices** – adds the sub‑image size to the start indices:

   .. code-block:: python

      xend = xstart + sub_naxis1
      yend = ystart + sub_naxis2

5. **Diagnostics** – prints header sizes and computed offsets.
6. **Return** – ``(xstart, xend, ystart, yend, flag_transpose)`` where
   ``flag_transpose`` is always ``False`` in the current version (the
   logic for detecting a transpose based on the ``PC`` matrix is
   commented out).

The function does **not** modify any data; it merely provides the slice
coordinates needed by later steps.

-----------------------------------------------------------------
4.2 ``empty_image(hdr, subh, extensions, fill, primary_hdr=None,
                  types=None, telescope=None, instrument=None)``

Creates an empty FITS HDU list with the dimensions of the full‑size image
and fills each requested extension with a constant value.

**Parameters**

* ``hdr`` – Header of the full‑size image (provides ``NAXIS1``,
  ``NAXIS2``, ``CRPIX*`` etc.).
* ``subh`` – Header of the sub‑image (used to copy WCS keywords into the
  new HDUs so that each extension inherits the correct coordinate system).
* ``extensions`` – List of extension names to create (e.g. ``["SCI",
  "ERR", "EXP", "WHT", "NIM"]``).
* ``fill`` – List of constant values with which each extension is initialised.
* ``primary_hdr`` – Optional primary header to store in the 0‑th HDU.
* ``types`` – Optional list of NumPy dtypes (``np.float32``, ``np.int16``,
  …) that control the datatype of each extension.
* ``telescope`` / ``instrument`` – Optional strings that, if supplied,
  overwrite the corresponding header keywords in **all** HDUs.

**What it does**

1. Determines the image size ``(NAXIS2, NAXIS1)`` from ``hdr``.
2. Makes a copy of ``subh`` and overwrites its ``NAXIS*`` and ``CRPIX*``
   values with those of the full‑size image – this ensures that each new
   HDU carries the correct WCS.
3. Builds an :class:`astropy.io.fits.HDUList`:
   * If a ``primary_hdr`` is supplied, a ``PrimaryHDU`` with that header
     is added as the first HDU; otherwise the list starts empty.
4. For each extension name and fill value:
   * Creates a NumPy array of the full size filled with the constant.
   * Casts it to the requested dtype (if ``types`` is given).
   * Appends an :class:`astropy.io.fits.ImageHDU` with the prepared data
     and the *full‑size* header created in step 2.
   * Sets the ``EXTNAME`` keyword to the supplied extension name.
5. If ``telescope`` or ``instrument`` were provided, the corresponding
   keywords are added/overwritten in **every** HDU.
6. Returns the populated ``HDUList`` ready to be written to disk.

-----------------------------------------------------------------
4.3 ``combine(outvals, invals_data, xstart, xend, ystart, yend,
             flag_transpose)``

Copies (or adds) the sub‑image data into the appropriate region of a full‑size
extension.

**Parameters**

* ``outvals`` – An :class:`astropy.io.fits.ImageHDU` (the destination
  extension of the full‑size image).
* ``invals_data`` – 2‑D NumPy array containing the sub‑image data that
  should be inserted.
* ``xstart``, ``xend``, ``ystart``, ``yend`` – Integer slice limits
  returned by :func:`find_slices`.
* ``flag_transpose`` – If ``True`` the sub‑image data is transposed
  before insertion (not used in the current script).

**Implementation**

.. code-block:: python

   if not flag_transpose:
       outvals.data[ystart:yend, xstart:xend] += invals_data
   else:
       outvals.data[ystart:yend, xstart:xend] += invals_data.T

The function **adds** the sub‑image values to any existing data
(``+=``) rather than overwriting them.  This behaviour is convenient when
building mosaics from many overlapping cut‑outs.  The function also
prints a few debugging statements about array shapes and slice limits.

-----------------------------------------------------------------
5.  The ``if __name__ == "__main__":`` block – overall workflow
-----------------------------------------------------------------

The script follows these high‑level steps when executed directly:

1. **Parse arguments** – using the parser defined in section 3.
2. **Determine which extensions** to process and the fill values based
   on ``--sci`` and ``--nim`` flags.
3. **Read the full‑size header**:
   * If ``--full_header`` points to a FITS file, the header of the
     ``PRIMARY`` or ``SCI`` extension is used (the script checks the
     ``TELESCOP`` keyword to decide).
   * If ``--full_header`` is a plain ``.hdr`` file, it is read with
     :meth:`fits.Header.fromtextfile`.
   * The script asserts that the destination ``--full_image`` does **not**
     already exist when a new header is supplied.
4. **Read the sub‑image header** – either from the supplied header file or
   directly from the ``SCI`` extension of ``--sub_image``.
5. **Calculate slice indices** – ``find_slices`` returns the region of the
   full image that corresponds to the sub‑image.
6. **Open the sub‑image FITS file** (``fits.open``) inside a ``with`` block.
7. **Create the full‑size image**:
   * The primary header of the sub‑image is copied to the new image.
   * ``empty_image`` builds an HDU list with the appropriate extensions,
     fill values, data types, and optional telescope/instrument keywords.
8. **Loop over each extension** (``SCI``, ``ERR``, …):
   * Optionally cast the ``NIM`` extension to ``int16``.
   * Compute *clipping* indices to handle cases where the sub‑image would
     extend beyond the borders of the full‑size mosaic.  The calculation
     uses ``np.max`` and ``np.min`` to ensure that only the overlapping
     region is copied.
   * Call ``combine`` with the appropriate slice limits.
9. **Propagate additional header keywords** – each HDU in the final image
   receives any missing keywords from the sub‑image header via
   ``hdu.header.extend(subh, unique=True)``.
10. **Write the result** – ``full_image.writeto(args.full_image)`` creates
    the new FITS file on disk.

-----------------------------------------------------------------
6.  Example usage scenarios
-----------------------------------------------------------------

**A.  Build a new empty mosaic and insert a single cut‑out**

.. code-block:: bash

   python mosaic_fits.py \
       --sub_image   cutout1.fits \
       --full_image  mosaic.fits \
       --full_header full_template.hdr \
       --telescope JWST \
       --instrument NIRCam

*The script creates ``mosaic.fits`` from the template header and inserts
the data from ``cutout1.fits``.*

**B.  Add a cut‑out to an existing mosaic (no new header)**

Assume ``existing_mosaic.fits`` already exists and contains a ``SCI``
extension.

.. code-block:: bash

   python mosaic_fits.py \
       --sub_image   cutout2.fits \
       --full_image  existing_mosaic.fits \
       --full_header ""      # empty string forces the “existing image” path

Because ``--full_header`` is omitted the script expects ``existing_mosaic.fits`` to
already be present; it will read the header from its ``SCI`` extension,
create a temporary in‑memory copy, add the new data and overwrite the file.

**C.  Only propagate the science data (no error, weight, …)**

.. code-block:: bash

   python mosaic_fits.py \
       --sub_image   cutout3.fits \
       --full_image  big_mosaic.fits \
       --full_header big_template.hdr \
       --sci

The ``--sci`` flag restricts processing to the ``SCI`` extension only.

-----------------------------------------------------------------
7.  Extending the script
-----------------------------------------------------------------

The current implementation is intentionally simple, but several
enhancements are straightforward:

* **Support for rotation or different pixel scales** – presently the
  script aborts if the CD matrix differs between the two headers.  A more
  general WCS reprojection (e.g. using :mod:`reproject`) could be added.
* **Weighted combination** – replace the ``+=`` in ``combine`` with an
  average or a variance‑weighted sum.
* **Parallel processing** – when inserting many cut‑outs, the loop over
  extensions could be parallelised with :mod:`concurrent.futures`.
* **Better header handling** – propagate additional keywords such as
  ``DATE‑OBS`` or ``EXPTIME`` while preserving provenance.

-----------------------------------------------------------------
8.  Summary
-----------------------------------------------------------------

* ``find_slices`` calculates where a sub‑image belongs inside a larger
  mosaic using WCS reference pixels.
* ``empty_image`` builds a blank FITS file of the correct size, filling
  each requested extension with a constant value and copying essential
  header information.
* ``combine`` inserts (or adds) the sub‑image data into the appropriate
  slice of the full‑size image, handling optional transposition.
* The ``main`` block ties everything together: argument parsing, header
  handling, slice computation, image creation, data insertion, header
  merging, and final FITS output.

The script provides a lightweight, command‑line driven way to assemble
large astronomical mosaics from many smaller cut‑outs while preserving
the scientific extensions required by downstream analysis pipelines.

-----------------------------------------------------------------
9.  Full source listing
-----------------------------------------------------------------

For completeness, the full source code is reproduced below (identical to
the original file).  The documentation above refers to the line numbers
and sections of this listing.

.. literalinclude:: mosaic_fits.py
   :language: python
   :linenos:
   :encoding: utf-8
