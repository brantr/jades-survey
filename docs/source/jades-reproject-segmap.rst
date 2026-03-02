.. _reproject_segmap_script:

=========================================
``jades-reproject-segmap.py``: Reproject Segmentation Map
=========================================

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

This document provides a detailed walkthrough of the ``reproject_segmap.py`` script.  
The script re‑projects a segmentation map (a FITS image that labels detected sources) onto the World Coordinate System (WCS) of a reference science image. It is intended for use in astronomical data reduction pipelines where source masks must be aligned with calibrated images.

----

.. contents::
   :depth: 2
   :local:

----

Overview
========

* **Purpose** – Take a segmentation map (``segmap.fits``) and re‑project it onto the WCS of another FITS file (typically a calibrated science image). The result is saved as a new FITS file (``segmap.reprojected.fits``).

* **Key library** – :mod:`reproject` (function ``reproject_interp``) which handles the geometric transformation between two WCS definitions.

* **Typical workflow**

  1. Parse command‑line arguments.
  2. Load the input segmentation map and its header.
  3. Load the header of the reference image (the target WCS).
  4. Re‑project the segmentation data using nearest‑neighbor interpolation.
  5. Write the re‑projected map to disk.

The script is deliberately minimal: it contains no error handling beyond what the underlying libraries provide, making it easy to embed in larger pipelines or to adapt for custom use‑cases.

----

Installation
============

The script depends on the following Python packages:

* ``numpy``
* ``astropy``
* ``reproject``

They can be installed via ``pip`` or ``conda``:

.. code-block:: bash

   # Using pip
   pip install numpy astropy reproject

   # Using conda
   conda install numpy astropy -c conda-forge
   conda install reproject -c conda-forge

The script itself does **not** require any additional compilation steps.

----

Command‑Line Interface
======================

The script uses :mod:`argparse` to expose four options:

.. code-block:: text

   usage: reproject_segmap.py [-h] [-i INPUT] [-o OUTPUT] [-f FITS] [-v]

   Reproject only the segmentation map.

   optional arguments:
     -h, --help            show this help message and exit
     -i INPUT, --input INPUT
                           Segmap image to reproject. (default: segmap.fits)
     -o OUTPUT, --output OUTPUT
                           Output, reprojected segmap. (default:
                           segmap.reprojected.fits)
     -f FITS, --fits FITS  Fits image with WCS for new projection. (default:
                           image.fits)
     -v, --verbose         Print helpful information to the screen? (default:
                           False)

* ``-i / --input`` – Path to the **source** segmentation map.
* ``-o / --output`` – Desired filename for the **re‑projected** map.
* ``-f / --fits`` – Path to the **reference** image that provides the target WCS. The script extracts the header from the ``SCI`` extension (common for calibrated data).
* ``-v / --verbose`` – Enable console logging of progress and timing information.

Example invocation:

.. code-block:: bash

   python reproject_segmap.py \
       --input my_segmap.fits \
       --fits calibrated_image.fits \
       --output my_segmap.reprojected.fits \
       --verbose

----

Code Walkthrough
================

The full source is reproduced below with explanatory comments inserted after each logical block.

.. code-block:: python
   :linenos:

   import sys
   # add current path
   sys.path.append('./')
   import argparse
   import numpy as np
   from astropy.io import fits
   from reproject import reproject_interp
   import time

   # ----------------------------------------------------------------------
   # Argument parser
   # ----------------------------------------------------------------------
   def create_parser():
       """
       Build and return an ``argparse.ArgumentParser`` configured for this script.
       """
       parser = argparse.ArgumentParser(
           description="Reproject only the segmentation map."
       )
       parser.add_argument('-i','--input',
           default='segmap.fits',
           metavar='input',
           type=str,
           help='Segmap image to reproject.')
       parser.add_argument('-o','--output',
           default='segmap.reprojected.fits',
           metavar='output',
           type=str,
           help='Output, reprojected segmap.')
       parser.add_argument('-f','--fits',
           default='image.fits',
           metavar='fits',
           type=str,
           help='Fits image with WCS for new projection.')
       parser.add_argument('-v', '--verbose',
                   dest='verbose',
                   action='store_true',
                   help='Print helpful information to the screen? (default: False)',
                   default=False)
       return parser

   # ----------------------------------------------------------------------
   # Main routine
   # ----------------------------------------------------------------------
   def main():
       # start global timer (unused later – kept for possible future extensions)
       time_start = time.time()

       # Parse command‑line arguments
       parser = create_parser()
       args   = parser.parse_args()

       # Verbose logging of received arguments
       if args.verbose:
           print(f"Input segmap = {args.input}.")
           print(f"FITS image with new WCS for reprojection  = {args.fits}.")
           print(f"Output segmap = {args.output}.")

       # --------------------------------------------------------------
       # Load the segmentation map
       # --------------------------------------------------------------
       hdu_detection = fits.open(args.input)
       header_old = hdu_detection[0].header
       # Force 2‑D header (segmentation maps are images, not data cubes)
       header_old['NAXIS'] = 2

       # --------------------------------------------------------------
       # Load the reference header (target WCS)
       # --------------------------------------------------------------
       # The script expects the science data in an extension named 'SCI'.
       header_new = fits.getheader(args.fits, 'SCI')

       # --------------------------------------------------------------
       # Re‑project the segmentation map
       # --------------------------------------------------------------
       # Reset timer for the reprojection step only
       time_start = time.time()

       # Update NAXIS1/NAXIS2 to match the actual data shape.
       # ``shape`` returns (ny, nx) for a 2‑D numpy array.
       header_old['NAXIS1'] = hdu_detection[0].data.shape[1]  # X dimension
       header_old['NAXIS2'] = hdu_detection[0].data.shape[0]  # Y dimension

       # ``reproject_interp`` expects a tuple (array, header) for the source.
       # ``order='nearest-neighbor'`` guarantees that integer labels are not
       # interpolated into fractional values, which would break the segmentation
       # map semantics.
       segm_deblend = reproject_interp(
           (hdu_detection[0].data, header_old),
           header_new,
           order='nearest-neighbor',
           return_footprint=False
       )

       time_end = time.time()
       if args.verbose:
           print(f"Time through segmap reprojection = {time_end-time_start}s.")

       # --------------------------------------------------------------
       # Write the output file
       # --------------------------------------------------------------
       if args.verbose:
           print("Saving reprojected detection image and segmentation map information...")
       # Cast back to integer type to preserve label values.
       fits.writeto(
           args.output,
           data=segm_deblend.astype(int),
           header=header_new,
           overwrite=True   # Overwrite existing file if present.
       )

   # ----------------------------------------------------------------------
   # Entry point
   # ----------------------------------------------------------------------
   if __name__ == "__main__":
       main()

Key points in the implementation
--------------------------------

* **Path handling** – ``sys.path.append('./')`` is a legacy hack that forces Python to look for modules in the current working directory. It can safely be removed if the script is run from its own folder or installed as a package.

* **Header manipulation** – ``NAXIS``, ``NAXIS1`` and ``NAXIS2`` are forced to match the data array dimensions. This is required because some FITS files contain stale header values that would otherwise confuse ``reproject_interp``.

* **Nearest‑neighbor interpolation** – Essential for segmentation maps where each pixel value encodes a region identifier. Linear interpolation would generate non‑integer labels and corrupt the mask.

* **Timing** – Two timers are used: a global one (currently unused) and a local timer around the reprojection step, printed only in verbose mode.

* **Overwrite behavior** – ``fits.writeto`` is called with ``overwrite=True`` (implicitly via the default in newer Astropy versions). Explicitly adding this argument makes the script robust when re‑running with the same output name.

----

Performance Considerations
==========================

* **Memory usage** – The entire source segmentation map and the re‑projected result are held in memory as NumPy arrays. For very large images (e.g., >10 000 × 10 000) ensure sufficient RAM is available.

* **CPU time** – The dominant cost is the call to ``reproject_interp``. Nearest‑neighbor interpolation is fast, but the overall speed also depends on the complexity of the WCS transformation (e.g., distortion terms). If speed becomes an issue, consider using the ``reproject`` function with ``order='nearest-neighbor'`` and ``chunksize`` to process the image in tiles.

* **Parallelism** – The underlying ``reproject`` library can make use of multiple cores when compiled with OpenMP. Verify that your installation of ``reproject`` is linked against a multithreaded NumPy build for maximal performance.

----

Extending the Script
====================

Common extensions include:

* **Support for other extensions** – Replace ``'SCI'`` with a user‑supplied extension name or index.
* **Alternative interpolation orders** – For scientific images (not segmentation maps) you may switch to ``order='bilinear'`` or ``order='spline'``.
* **Batch processing** – Wrap the ``main`` logic in a loop that iterates over a list of segmentation maps.
* **Error handling** – Add ``try/except`` blocks around file I/O and reprojection to produce graceful messages when files are missing or headers are malformed.
* **Logging** – Replace ``print`` statements with the :mod:`logging` module for configurable verbosity levels.

Example: Adding a ``--extension`` argument

.. code-block:: python

   parser.add_argument('--extension', default='SCI',
                       help='Extension name in the reference FITS file that contains the WCS.')

   # later
   header_new = fits.getheader(args.fits, args.extension)

----

Testing the Script
==================

A minimal test can be performed with synthetic data:

.. code-block:: python

   import numpy as np
   from astropy.io import fits
   from astropy.wcs import WCS

   # Create a dummy segmentation map (10×10) with two regions
   seg = np.zeros((10, 10), dtype=int)
   seg[:5, :] = 1
   seg[5:, :] = 2
   hdu = fits.PrimaryHDU(seg)
   hdu.writeto('dummy_seg.fits', overwrite=True)

   # Create a reference image with a simple WCS (shift by 2 pixels)
   w = WCS(naxis=2)
   w.wcs.crpix = [5, 5]
   w.wcs.cdelt = np.array([-0.0002777778, 0.0002777778])
   w.wcs.crval = [0, 0]
   w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
   header = w.to_header()
   header['NAXIS'] = 2
   header['NAXIS1'] = 12
   header['NAXIS2'] = 12
   fits.PrimaryHDU(header=header).writeto('dummy_ref.fits', overwrite=True)

   # Run the script
   # $ python reproject_segmap.py -i dummy_seg.fits -f dummy_ref.fits -v

The output ``segmap.reprojected.fits`` should contain the same two regions, shifted according to the reference WCS.

----

License
=======

This script is released under the **MIT License**. See the `LICENSE <https://opensource.org/licenses/MIT>`_ file for full terms.

----

Indices and tables
===================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
