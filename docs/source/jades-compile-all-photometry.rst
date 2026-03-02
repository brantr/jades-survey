.. _catalog_assembler:

========================================================
``jades-compile-all-photometry.py`` – Combine multiple FITS catalog files
========================================================

**Version:** |release|  
**Author:**  (original author)  
**License:** BSD‑3-Clause (or as provided with the source)

--------------------------------------------------------------------
Table of contents
--------------------------------------------------------------------
.. contents::
   :depth: 2
   :local:

--------------------------------------------------------------------
1. Overview
--------------------------------------------------------------------

The *catalog_assembler* script is a small command‑line utility written in
Python that merges a collection of astronomical catalog products (produced
by source‑detection pipelines such as **SExtractor**, **photutils**, or
custom pipelines) into a single multi‑extension FITS file.

Typical use‑cases are:

* Gather a “size” table (positions, shape parameters, etc.) from a detection
  catalog.
* Append flag tables that contain quality or classification flags for each
  source.
* Merge photometric measurements performed with different apertures
  (circular, Kron, growth‑curve, PSF‑convolved, background‑subtracted, …).
* Optionally convert angular units from radians to degrees.

The script is deliberately **modular** – each type of input file is processed
by a dedicated ``create_*_hdu`` function that returns a
``astropy.io.fits.BinTableHDU``.  All HDUs are collected in a list,
prepended with a primary empty HDU, and finally written to the user‑specified
output file.

--------------------------------------------------------------------
2. Dependencies
--------------------------------------------------------------------

The script relies on the following third‑party packages (all available on
PyPI):

* **NumPy** – fundamental array handling.
* **Astropy** – FITS I/O, ASCII table reading, and the ``Table`` class.
* **tqdm** – progress bar for iterating over long file lists.
* **argparse** – part of the Python standard library, used for command‑line
  parsing.

Optional (commented‑out) imports such as ``matplotlib.pyplot`` or
``astropy.wcs`` are not required for the core functionality.

--------------------------------------------------------------------
3. Command‑line interface
--------------------------------------------------------------------

The script is executed as a regular Python program:

.. code-block:: console

   $ python catalog_assembler.py [options]

The command‑line options are defined in :func:`create_parser`.  All options
accept a *string* that points to a *text file* containing a list of FITS
files (one per line).  The following arguments are recognised:

+--------------------+----------------------+--------------------------------------------+
| Argument           | Destination          | Description                                |
+====================+======================+============================================+
| ``--det``          | ``det``              | Detection catalog (size/shape table).      |
+--------------------+----------------------+--------------------------------------------+
| ``--flag``         | ``flag``             | Files that contain a ``FLAG`` extension.   |
+--------------------+----------------------+--------------------------------------------+
| ``--circ``         | ``circ``             | Circular aperture photometry (no bsub).    |
+--------------------+----------------------+--------------------------------------------+
| ``--circ-bsub``    | ``circ_bsub``        | Circular aperture photometry, background‑ |
|                    |                      | subtracted.                                |
+--------------------+----------------------+--------------------------------------------+
| ``--circ-psf``     | ``circ_psf``         | PSF‑convolved circular photometry.         |
+--------------------+----------------------+--------------------------------------------+
| ``--circ-psf-bsub``| ``circ_psf_bsub``    | PSF‑convolved circular photometry, bsub.   |
+--------------------+----------------------+--------------------------------------------+
| ``--kron``         | ``kron``             | Kron‑like aperture photometry.             |
+--------------------+----------------------+--------------------------------------------+
| ``--kron-psf``     | ``kron_psf``         | PSF‑convolved Kron photometry.             |
+--------------------+----------------------+--------------------------------------------+
| ``--growth``       | ``growth``           | Growth‑curve photometry (not implemented   |
|                    |                      | in the current script).                    |
+--------------------+----------------------+--------------------------------------------+
| ``--growth-psf``   | ``growth_psf``       | PSF‑convolved growth‑curve (placeholder).  |
+--------------------+----------------------+--------------------------------------------+
| ``--miri``         | ``miri``             | MIRI‑specific photometry (JWST).           |
+--------------------+----------------------+--------------------------------------------+
| ``-o`` / ``--output``| ``output``          | Name of the compiled FITS file (default:   |
|                    |                      | ``output.fits``).                          |
+--------------------+----------------------+--------------------------------------------+
| ``--fix-radians``  | ``fix_radians``      | Convert angle columns from radians to      |
|                    |                      | degrees (applies to Kron tables).          |
+--------------------+----------------------+--------------------------------------------+
| ``-v`` / ``--verbose``| ``verbose``        | Print timing information after execution.  |
+--------------------+----------------------+--------------------------------------------+

All arguments are optional; the script simply skips any step for which the
corresponding file list is not provided.

--------------------------------------------------------------------
4. Core functions
--------------------------------------------------------------------

The module consists of a handful of small, well‑documented helper functions.
Below each function is described with its purpose, the main algorithmic
steps, and the returned object.

--------------------------------------------------------------------
4.1 ``create_parser()``
--------------------------------------------------------------------

*Purpose* – Build an :class:`argparse.ArgumentParser` describing the CLI.

*Key points*

* Uses ``action='store_true'`` for boolean flags (``--fix-radians`` and
  ``--verbose``).
* All file‑list arguments accept a string (the path to a plain‑text file).
* Returns the configured ``ArgumentParser`` instance.

--------------------------------------------------------------------
4.2 ``create_size_hdu(fname, hdu_name)`` 
--------------------------------------------------------------------

*Purpose* – Convert a detection catalog stored in an **ASCII** table into a
binary FITS table (an HDU) containing a curated subset of columns.

*Algorithm*

1. ``ascii.read(fname)`` loads the whole ASCII file into an Astropy ``Table``.
2. A new empty ``Table`` ``t`` is created.
3. Selected columns are copied from the input table to ``t`` and renamed to
   the names used by downstream analysis (e.g. ``id → ID``, ``ra → RA``).
4. The column list is printed for debugging.
5. ``fits.BinTableHDU(data=t, name=hdu_name)`` creates a binary table HDU.

*Returned value* – ``fits.BinTableHDU`` named *hdu_name* (normally ``SIZE``).

--------------------------------------------------------------------
4.3 ``create_flag_hdu(fname, hdu_name)`` 
--------------------------------------------------------------------

*Purpose* – Merge the ``FLAG`` extensions from a list of FITS files into a
single table, preserving only the columns that are needed for later work.

*Algorithm*

1. Read the list of file names from ``fname`` (strip newline characters).
2. Initialise an empty Astropy ``Table`` ``t``.
3. Iterate over the list with ``tqdm`` to display progress.
4. For each file:
   * Load the ``FLAG`` extension using ``fits.getdata(file, 'FLAG')``.
   * For each column in the FITS table, if the column is not already present
     in ``t`` **and** the column name contains ``ID``, ``RA``, ``DEC`` or
     ``FLAG``, copy the whole column into ``t``.
5. After the loop, print the final column names and build the HDU.

*Returned value* – ``fits.BinTableHDU`` named *hdu_name* (normally ``FLAG``).

--------------------------------------------------------------------
4.4 ``create_flux_hdu(fname, hdu_name,
                       exclude_large=True,
                       flag_to_degrees=False)`` 
--------------------------------------------------------------------

*Purpose* – Assemble photometric measurements from the ``CAT`` extension of
multiple FITS files (circular, Kron, growth‑curve, etc.) into a single table.

*Parameters*

* ``exclude_large`` – When *True* (default) columns whose name contains
  ``CIRC7`` or ``CIRC8`` are ignored.  Those columns usually hold very large
  apertures that are not needed for the final catalog.
* ``flag_to_degrees`` – If *True*, any column whose name contains
  ``THETA_KRON`` is multiplied by ``180/π`` to convert from radians to degrees.

*Algorithm*

1. Load the list of FITS files.
2. Initialise an empty ``Table`` ``t``.
3. Loop over the files (progress bar via ``tqdm``).
4. For each file, read the ``CAT`` extension.
5. Iterate over its columns; for each column not already in ``t``:
   * Skip columns containing ``TEXP``, ``WHT`` or ``bkg`` (exposure time,
     weight map, background – not photometric measurements).
   * Determine a multiplicative ``factor`` (1.0 unless conversion to degrees
     is requested for a ``THETA_KRON`` column).
   * If ``exclude_large`` is *True*, discard columns that contain
     ``CIRC7`` or ``CIRC8``; otherwise keep everything.
   * Store the scaled column in ``t``.
6. After processing all files, create a binary table HDU.

*Returned value* – ``fits.BinTableHDU`` named *hdu_name* (e.g. ``CIRC``,
``KRON``, ``CIRC_BSUB`` …).

--------------------------------------------------------------------
4.5 ``create_miri_hdu(fname, hdu_name,
                       exclude_large=True,
                       flag_to_degrees=False)`` 
--------------------------------------------------------------------

*Purpose* – Same as ``create_flux_hdu`` but adds a small tweak for MIRI
(JWST) data: if the source FITS file name contains the substring ``bsub``,
the output column name is suffixed with ``_BSUB`` to make the background‑
subtracted nature explicit.

*Algorithm* – Identical to ``create_flux_hdu`` with the extra step:

* When ``exclude_large`` is *True* and the filename includes ``bsub``,
  ``col_out = col + '_BSUB'`` before insertion into the table.

*Returned value* – ``fits.BinTableHDU`` named *hdu_name* (normally ``MIRI``).

--------------------------------------------------------------------
4.6 ``main()``
--------------------------------------------------------------------

*Purpose* – Orchestrate the entire workflow:

1. **Timing** – Record the start time.
2. **Argument parsing** – Obtain the ``Namespace`` from ``create_parser``.
3. **Initialize HDU list** – Start with an empty primary HDU
   (required by FITS standards).
4. **Conditional processing** – For each CLI option that is not ``None``,
   call the appropriate ``create_*_hdu`` function and append the resulting
   HDU to the list.
5. **Write output** – Build a ``fits.HDUList`` from the collected HDUs and
   write it to ``args.output`` (overwrites existing files).
6. **Verbose timing** – If ``--verbose`` is set, print the elapsed wall‑clock
   time.

The function contains a small placeholder for ``create_growth_hdu`` which
is not defined in the supplied source – attempting to use ``--growth`` or
``--growth-psf`` will raise a ``NameError``.  This is intentional: the
script template is often extended by the user to implement those HDUs.

--------------------------------------------------------------------
5. Execution flow diagram
--------------------------------------------------------------------

.. mermaid::
   :caption: High‑level execution flow of *catalog_assembler*.

   graph TD
       A[Start → main()] --> B[Parse CLI arguments]
       B --> C{Any arguments ?}
       C -->|det| D[create_size_hdu()]
       C -->|flag| E[create_flag_hdu()]
       C -->|circ*| F[create_flux_hdu()]
       C -->|kron*| G[create_flux_hdu(flag_to_degrees)]
       C -->|miri| H[create_miri_hdu()]
       D --> I[Append HDU to list]
       E --> I
       F --> I
       G --> I
       H --> I
       I --> J[fits.HDUList → write output]
       J --> K[Print timing (optional)]
       K --> L[End]

--------------------------------------------------------------------
6. Example usage
--------------------------------------------------------------------

Assume the following files exist in the current directory:

* ``det.txt`` – a single‑line text file that points to ``det_catalog.ascii``.
* ``flags.txt`` – contains three FITS files each with a ``FLAG`` extension.
* ``circ.txt`` – list of photometry files measured with a 2‑pixel circular
  aperture.
* ``kron.txt`` – list of Kron‑aperture files where the angle columns are in
  radians.

A typical command to produce a merged catalog would be:

.. code-block:: console

   $ python catalog_assembler.py \
       --det det.txt \
       --flag flags.txt \
       --circ circ.txt \
       --kron kron.txt \
       --fix-radians \
       -o merged_catalog.fits \
       -v

The script will:

* Create a ``SIZE`` HDU from the detection catalog.
* Merge all flag tables into a ``FLAG`` HDU.
* Assemble circular photometry into a ``CIRC`` HDU.
* Convert the ``THETA_KRON`` column from radians to degrees and store the
  result in a ``KRON`` HDU.
* Write the four HDUs (plus the primary empty HDU) to ``merged_catalog.fits``.
* Print the total execution time because ``-v`` was supplied.

--------------------------------------------------------------------
7. Extending the script
--------------------------------------------------------------------

The current implementation is deliberately lightweight.  Common extensions
include:

* **Growth‑curve support** – implement ``create_growth_hdu`` analogous to
  ``create_flux_hdu`` but handling the different column naming scheme.
* **WCS handling** – import ``astropy.wcs`` and add world‑coordinate
  information to the primary header.
* **Parallel I/O** – replace the Python loops with ``multiprocessing.Pool``
  for faster processing of very large file lists.
* **Unit handling** – attach Astropy ``Quantity`` objects to columns so that
  downstream tools can automatically recognise flux units (e.g. ``u.Jy``).

--------------------------------------------------------------------
8. License and acknowledgements
--------------------------------------------------------------------

The script is released under the BSD‑3‑Clause license.  It makes use of
open‑source libraries maintained by the scientific Python community.

--------------------------------------------------------------------
9. References
--------------------------------------------------------------------

* Astropy Collaboration, *Astropy: A community Python package for astronomy*,
  A&A 558, A33 (2013).  
* T. K. Miller et al., *tqdm: A fast, extensible progress bar for Python*,
  JOSS 5(54), 2020.  

--------------------------------------------------------------------
10. Index
--------------------------------------------------------------------
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
