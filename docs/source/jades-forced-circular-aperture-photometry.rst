.. _photometry_script:

========================================
``jades-forced-circular-aperture-photometry.py``: Circular Aperture Photometry Script
========================================

|Version| |License|

**Author:** *Original author of the script*  
**Date:** 2025‑08‑12  

This document describes the behaviour of the Python script that performs
circular aperture photometry on JWST/HST images.  It is intended for
inclusion in a *Read the Docs* website; therefore the file is written in
reStructuredText (``.rst``) format.

--------------------------------------------------------------------
Table of contents
--------------------------------------------------------------------

.. contents::
   :local:
   :depth: 2


--------------------------------------------------------------------
1. Overview
--------------------------------------------------------------------

The script is a command‑line tool that reads a science image (``SCI``),
its associated error, weight and exposure‑time extensions (or separate
files), a source catalog, and optionally a segmentation map.  For each
source in the catalog it:

* converts pixel coordinates from world coordinates (RA/Dec) using the
  image WCS,
* builds a circular aperture of a user‑specified radius,
* measures the weighted exposure time, the weight, the raw flux,
* applies an aperture correction derived from a supplied PSF,
* optionally performs a local background subtraction using an annulus,
* estimates the total uncertainty (photon + background + RMS), and
* writes all results to a FITS binary table.

The script is deliberately instrument‑agnostic: it works with HST,
JWST/NIRCam and JWST/MIRI data by inspecting the FITS header.

--------------------------------------------------------------------
2. Module imports
--------------------------------------------------------------------

The script imports the following third‑party libraries:

* ``numpy`` – numerical arrays.
* ``astropy.io.fits`` – reading and writing FITS files.
* ``astropy.io.ascii`` – reading the input catalog (ASCII table).
* ``astropy.table.Table`` – constructing the output table.
* ``astropy.wcs.WCS`` – World Coordinate System handling.
* ``argparse`` – parsing command‑line arguments.
* ``tqdm`` – progress bars for loops.
* ``photutils`` – aperture definitions and photometry utilities.
* ``time`` – simple timing of the whole program.

No modifications are made to these imports; they are used throughout
the script.

--------------------------------------------------------------------
3. Functions
--------------------------------------------------------------------

Each function is documented separately below.  The order follows the
definition order in the source file.

--------------------------------------------------------------------
3.1 ``create_parser()``
--------------------------------------------------------------------

**Purpose**

Creates and returns an ``argparse.ArgumentParser`` instance that defines
all command‑line options accepted by the script.

**Behaviour**

* The parser description is ``"Flags and options from user."``.
* The following arguments are added (default values shown in *italics*):

  - ``-a / --aper`` – aperture identifier (e.g. ``CIRC1``). Default:
    ``'CIRC1'``.
  - ``-i / --input`` – primary science FITS image. Default:
    ``'input.fits'``.
  - ``--segmap`` – segmentation map FITS image (used for background
    subtraction). Default: ``'segmap.fits'``.
  - ``-c / --catalog`` – source catalog (ASCII). Default:
    ``'cat.txt'``.
  - ``-e / --error`` – error image FITS file. Default:
    ``'err.fits'``.
  - ``--wht`` – optional weight image FITS file. Default: ``None``.
  - ``--wht-ext`` – extension name for the weight image (if a multi‑
    extension FITS). Default: ``None``.
  - ``--texp`` – exposure‑time image FITS file. Default: ``None``.
  - ``--exp-ext`` – extension name for the exposure‑time image.
    Default: ``None``.
  - ``-o / --output`` – name of the output FITS catalog. Default:
    ``'output.fits'``.
  - ``--psf`` – PSF image used for aperture correction. Default:
    ``'psf.fits'``.
  - ``--band`` – photometric band identifier (e.g. ``f444w``). Default:
    ``'f444w'``.
  - ``--no-aper-corr`` – a *store_false* flag that disables aperture
    correction.  When present ``args.aper_corr`` becomes ``False``.
    Default value is ``True``.
  - ``--bsub`` – a *store_true* flag that activates background
    subtraction using the segmentation map. Default ``False``.
  - ``-v / --verbose`` – a *store_true* flag that prints the total
    execution time at the end. Default ``False``.

* The function returns the fully configured ``ArgumentParser`` object.

--------------------------------------------------------------------
3.2 ``get_instrument(header_img, instrument=None)``
--------------------------------------------------------------------

**Purpose**

Determine which instrument (``HST``, ``NIRCAM`` or ``MIRI``) produced the
image based on its FITS header.

**Behaviour**

* If the optional ``instrument`` argument is supplied, it is used to set
  ``flag_inst`` directly (with a fallback to ``'HST'`` for any unknown
  value).
* When ``instrument`` is ``None`` the function inspects the header:

  - If the keyword ``TELESCOP`` exists and its value is ``'HST'``,
    ``flag_inst`` is set to ``'HST'`` and a message is printed.
  - Otherwise the function tries to read ``INSTRUME`` (which for JWST
    images will be ``'NIRCAM'`` or ``'MIRI'``) and uses that value.
  - If neither keyword is present, it defaults to ``'MIRI'``.

* The resulting string ``flag_inst`` is returned.

--------------------------------------------------------------------
3.3 ``get_image(args, ext='SCI', flag_inst='JWST')``
--------------------------------------------------------------------

**Purpose**

Read a specific data extension (science, error, weight, exposure, or RMS)
from the FITS files supplied via the command‑line arguments.

**Parameters**

* ``args`` – the namespace returned by ``parser.parse_args()``.
* ``ext`` – one of ``'SCI'``, ``'ERR'``, ``'WHT'``, ``'EXP'`` or ``'RMS'``.
* ``flag_inst`` – instrument identifier (as returned by
  ``get_instrument``).

**Behaviour**

The function follows a series of conditional branches:

* ``ext == 'SCI'``  

  - For HST images the science data are read from the primary HDU:
    ``fits.getdata(args.input)``.
  - For JWST images the data are taken from the ``'SCI'`` extension of
    ``args.input``.

* ``ext == 'ERR'``  

  - For HST the error image is read from ``args.err`` (primary HDU).
  - For JWST the error is taken from the ``'ERR'`` extension of
    ``args.input`` (the same file as the science image).

* ``ext == 'WHT'``  

  - If ``args.wht`` is ``None`` the weight is read from the ``'WHT'``
    extension of ``args.input``.
  - Otherwise ``fits.getdata(args.wht, args.wht_ext)`` is used.

* ``ext == 'EXP'``  

  - If ``args.texp`` is ``None`` the exposure is read from the ``'EXP'``
    extension of ``args.input``.
  - Otherwise ``fits.getdata(args.texp, args.exp_ext)`` is used.

* ``ext == 'RMS'``  

  - The RMS image is read from ``args.err`` (the error file).

Finally the data are cast to ``np.float32`` and returned.

--------------------------------------------------------------------
3.4 ``compute_flux_to_nJy(header_flux, flag_inst)``
--------------------------------------------------------------------

**Purpose**

Calculate the multiplicative factor that converts the image units to
nano‑Jansky (nJy).  The function also returns the pixel scale and a WCS
object derived from the header.

**Behaviour**

1. A ``WCS`` object is created from ``header_flux``.
2. The pixel scale (arcseconds per pixel) is obtained from the WCS
   transformation matrix using ``wcs.pixel_scale_matrix`` implicitly via
   ``np.sqrt(np.sum(wcs.pixel_scale_matrix**2, axis=0))`` (the script
   does not compute it explicitly; the returned ``pixel_scale`` is the
   ``abs`` value of the CD matrix diagonal).
3. The conversion factor ``flux_to_nJy`` is set according to the
   instrument:

   * **HST** – the script leaves the factor at its default value
     (``1.0``) because HST images are already in calibrated flux units.
   * **JWST (NIRCam or MIRI)** – the factor is obtained from the header
     keyword ``PHOTFNU`` (or a similar keyword).  In the script the
     factor is stored directly in ``flux_to_nJy``; however the actual
     conversion from MJy sr⁻¹ to nJy is performed later using the
     exposure time.

4. The function appends a header keyword ``('FNUTONJY', flux_to_nJy,
   'Multiply to convert MJySR to nJy.')`` when the script writes the
   output catalog.

5. The function returns a three‑element tuple:

   ``(wcs, flux_to_nJy, pixel_scale)`` where

   * ``wcs`` – an ``astropy.wcs.WCS`` instance,
   * ``flux_to_nJy`` – float conversion factor,
   * ``pixel_scale`` – pixel size in arcseconds.

--------------------------------------------------------------------
3.5 ``get_aperture_radius_and_correction(args, flag_inst)``
--------------------------------------------------------------------

*The actual function name in the source is
``get_aperture_radius_and_correction``; the documentation follows that
exact name.*

**Purpose**

Derive the aperture radius (in pixels) corresponding to the user‑chosen
aperture identifier and compute the aperture correction using the PSF
image.

**Behaviour**

1. The PSF image is loaded with ``fits.getdata(args.psf)``.
2. The script extracts the central part of the PSF to estimate the
   encircled energy for the requested aperture.  The exact algorithm
   is:

   * The PSF image dimensions are used to locate the centre.
   * ``photutils`` creates a circular aperture with the radius
     ``aperture_radius`` that corresponds to the identifier
     (e.g. ``CIRC1`` → 1 pixel, ``CIRC0`` → variable size, etc.).
   * The sum of the PSF inside this aperture is compared with the total
     PSF flux (sum of the whole image) to obtain the correction factor.

3. If the user disables aperture correction via ``--no-aper-corr``,
   the function still returns the computed correction but the calling
   code will ignore it.

4. The function returns a tuple ``(r_pixel, aperture_correction)``:

   * ``r_pixel`` – aperture radius in pixels (float).
   * ``aperture_correction`` – multiplicative factor that converts a
     measured flux to a total flux (float).

--------------------------------------------------------------------
3.6 ``compute_rms_uncertainties(args, mask, apertures, r_pixel,
flag_inst, flux_to_nJy)``
--------------------------------------------------------------------

**Purpose**

Estimate the RMS uncertainty for a given aperture size when the user
asks for the special ``CIRC0`` aperture (the smallest aperture, for which
the RMS image does not contain a pre‑computed value).

**Behaviour**

* For apertures other than ``CIRC0`` the function is essentially a thin
  wrapper around reading the RMS image and applying the aperture
  correction.
* When ``args.aper`` equals ``CIRC0`` the function:

  1. Calls ``load_rms_uncertainties`` (or performs its own interpolation)
     to obtain the per‑aperture RMS values.
  2. Multiplies the RMS by the aperture correction.
  3. Returns the resulting array.

* The function is **not** used in the final program flow; instead the
  script calls ``load_rms_uncertainties`` directly.  The code path that
  would invoke this function is commented out.

--------------------------------------------------------------------
3.7 ``load_rms_uncertainties(args, r_pixel, naps)``
--------------------------------------------------------------------

**Purpose**

Read the RMS (or error) values stored in the error image (``ERR`` or
``RMS``) and provide an uncertainty estimate for each source.  For the
standard apertures (``CIRC1``–``CIRC8``) the uncertainties are read
directly from the RMS image; for ``CIRC0`` the function interpolates.

**Parameters**

* ``args`` – parsed command‑line arguments.
* ``r_pixel`` – aperture radius in pixels (float).
* ``naps`` – total number of sources (length of the catalog).

**Behaviour**

1. The function calls ``load_rms_uncertainties`` (the name in the source
   is ``load_rms_uncertainties``) which:

   * Opens the RMS/ERROR FITS file via ``get_image`` with ``ext='ERR'``
     (or ``'RMS'`` when appropriate).
   * For each source it extracts the RMS value inside the aperture
     using ``ApertureStats``.
   * If the aperture identifier is ``CIRC0`` the function performs a
     linear interpolation of the RMS as a function of aperture radius.
     The implementation uses the source radius ``r_pixel`` to find the
     two nearest pre‑computed RMS tables and interpolates between them.

2. The resulting 1‑D ``numpy`` array ``drms`` (length ``naps``) is
   returned.  No further scaling is performed here; the calling code
   applies the aperture correction later.

--------------------------------------------------------------------
3.8 ``load_rms_uncertainties(args, r_pixel, naps)``
--------------------------------------------------------------------

**Purpose**

Read a pre‑computed RMS uncertainty table from the error FITS file (or
interpolate for ``CIRC0``) and return a per‑source array of uncertainties
in nJy.

**Behaviour**

* Calls ``load_rms_uncertainties`` (the function defined in the source
  code) which:

  1. Opens the RMS/ERROR image using ``get_image(args, ext='ERR')``.
  2. For each source it extracts the RMS value inside the user‑defined
     aperture using ``ApertureStats``.
  3. If the aperture identifier is not ``CIRC0`` the values are taken
     directly from the image; otherwise the function interpolates as
     described in the previous section.

* Returns a ``numpy`` array of length ``naps`` containing the RMS
  uncertainties (in nJy) **before** aperture correction.

--------------------------------------------------------------------
3.9 ``main()``
--------------------------------------------------------------------

**Purpose**

Coordinate the whole workflow: parse arguments, read data, perform the
photometry, compute uncertainties, write the output catalog, and report
timing information if requested.

**Step‑by‑step behaviour**

1. **Argument parsing** – ``args = parser.parse_args()``.
2. **Timer start** – ``time_global_start = time.time()``.
3. **Instrument detection** – ``flag_inst = get_instrument(header)``.
4. **Aperture definition**

   * Reads the PSF image and calls
     ``get_aperture_radius_and_correction`` to obtain ``r_pixel`` and
     ``aperture_correction``.
   * Prints the chosen aperture radius in pixels and arcseconds.

5. **Source catalog handling**

   * Loads the ASCII catalog with ``ascii.read``.
   * Extracts world coordinates (RA, DEC) and source identifiers.
   * Converts RA/DEC to pixel coordinates using the WCS of the science
     image.

6. **Weight (WHT) measurement**

   * Reads the weight image (``WHT``) via ``get_image``.
   * Uses ``photutils.ApertureStats`` to compute the mean weight inside
     each aperture.
   * Stores the result in the output table column ``{BAND}_WHT``.

7. **Exposure time (EXP) measurement**

   * Reads the exposure‑time image (``EXP``) via ``get_image``.
   * Loops over each aperture, applying ``ApertureStats`` to obtain the
     mean exposure time per source.
   * Saves the array as ``{BAND}_TEXP`` in the final table.

8. **RMS uncertainty handling**

   * Calls ``load_rms_uncertainties`` (or the commented out
     ``compute_rms_uncertainties``) to obtain ``drms`` – an array of
     per‑source RMS uncertainties in nJy.

9. **Flux measurement**

   * Reads the science image (``SCI``) and, for JWST, the error image
     (``ERR``) via ``get_image``.
   * Scales both arrays by the conversion factor ``flux_to_nJy`` to
     obtain fluxes in nJy.
   * Extends the mask to flag any NaN or infinite values.
   * Performs aperture photometry with
     ``photutils.aperture_photometry``.  For JWST the ``error=`` keyword
     supplies the per‑pixel error image; for HST it is omitted.

10. **Optional local background subtraction** (activated with
    ``--bsub``)

    * Constructs a thin annulus (inner radius 1.5″, outer radius 1.55″)
      around each source.
    * Masks all segmentation regions, then temporarily unmasks the
      current source to avoid self‑contamination.
    * Measures the *median* background level inside the annulus using
      ``ApertureStats`` and multiplies by the aperture area
      (π r²) to obtain a background flux in nJy.
    * Stores the background values in the table column
      ``{BAND}_{aper}_bkg`` (if background subtraction is requested).

11. **Flux conversion to instrumental units**

    * Computes ``nJy_to_electrons = texp_source / flux_to_nJy``.
    * For JWST instruments the result is further divided by the header
      keyword ``PHOTMJSR`` to obtain electrons per nJy.
    * Raw fluxes ``F`` (nJy) are multiplied by ``nJy_to_electrons`` to
      obtain counts ``F_in_counts``.
    * RMS uncertainties are also converted to counts.

12. **Total uncertainty**

    * Calculates the quadrature sum of the RMS noise (in counts) and
      Poisson noise (approximated by ``|F_in_counts|``):

      ``total_uncertainty_counts = sqrt(rms_counts² + |F_counts|)``

    * Converts back to nJy by dividing by ``nJy_to_electrons``.
    * Applies the aperture correction factor to both flux and
      uncertainty.

13. **Table assembly**

    The script builds an ``astropy.table.Table`` ``t`` with the
    following columns (all stored as ``float32`` unless otherwise noted):

    * ``ID``, ``RA``, ``DEC`` – copied from the input catalog.
    * ``{BAND}_WHT`` – mean weight inside each aperture.
    * ``{BAND}_TEXP`` – exposure time per source (seconds).
    * ``{BAND}_{aper}`` – aperture‑corrected flux in nJy.
    * ``{BAND}_{aper}_bkg`` – background flux (if ``--bsub`` is used).
    * ``{BAND}_{aper}_e`` – total uncertainty (nJy) after aperture
      correction and background subtraction.
    * ``{BAND}_{aper}_ei`` – *single‑pixel* error estimate:
      * for HST this is the RMS value from the error image,
      * for JWST it is the error column returned by
        ``aperture_photometry`` (``aperture_sum_err``).

14. **FITS output**

    * Generates a primary HDU (empty) and a binary table HDU named
      ``'CAT'`` containing the table ``t``.
    * Constructs a header from the WCS with additional keywords that
      record the input files, the chosen band, aperture, pixel scale,
      conversion factor, and whether an aperture correction was applied.
    * Writes the HDU list to the file indicated by ``args.output``,
      overwriting any existing file.

15. **Timing and verbose output**

    * If ``args.verbose`` is ``True`` the total execution time is printed.

--------------------------------------------------------------------
3.10 ``if __name__ == "__main__":``
--------------------------------------------------------------------

When the script is executed as a stand‑alone program, the ``main()``
function is called, launching the full photometry pipeline described
above.

--------------------------------------------------------------------
4. Command‑line usage example
--------------------------------------------------------------------

Assuming the script is saved as ``aperture_photometry.py`` and is
available on the ``PATH``:

.. code-block:: bash

   python aperture_photometry.py \
       --input   jwst_nircam_f444w_sci.fits \
       --error   jwst_nircam_f444w_err.fits \
       --wht     jwst_nircam_f444w_wht.fits \
       --texp    jwst_nircam_f444w_exp.fits \
       --psf     jwst_nircam_f444w_psf.fits \
       --catalog source_list.txt \
       --band    f444w \
       --aper    CIRC2 \
       --bsub \
       --output  f444w_circ2_catalog.fits \
       -v

The above command will:

* use the ``CIRC2`` aperture,
* perform a local background subtraction,
* write verbose timing information, and
* store the resulting catalog in ``f444w_circ2_catalog.fits``.

--------------------------------------------------------------------
5. Frequently asked questions
--------------------------------------------------------------------

**Q:** *Why is an aperture correction needed?*  
**A:** The PSF of a space‑based telescope spreads flux over many pixels.
A small circular aperture captures only a fraction of the total source
light.  The script measures the encircled energy of the supplied PSF
within the chosen aperture and computes ``aperture_correction`` = 1 /
(encircled‑energy).  Multiplying the measured flux by this factor
produces an estimate of the total source flux.

**Q:** *What does the ``--no-aper-corr`` flag do?*  
**A:** It sets ``args.aper_corr`` to ``False``.  The script will still
measure the raw aperture sum, but it will *not* multiply by the
aperture correction factor before writing the output columns.

**Q:** *Can I use this script with images that have different extension
names?*  
**A:** Yes.  Use the ``--wht-ext`` and ``--exp-ext`` arguments to point
to the correct HDU names when the weight or exposure images are stored
in non‑standard extensions.

**Q:** *How is the total uncertainty computed?*  
**A:**  

``total_uncertainty = sqrt( rms_uncertainty² + |flux_counts| )``  

where ``rms_uncertainty`` is the RMS noise inside the aperture (converted
to counts) and ``|flux_counts|`` approximates the Poisson noise
(``σ = sqrt(N)``) assuming the measured flux in electrons follows a
Poisson distribution.

--------------------------------------------------------------------
6. References
--------------------------------------------------------------------

* Astropy: https://www.astropy.org/
* Photutils: https://photutils.readthedocs.io/
* JWST Calibration Documentation: https://jwst-pipeline.readthedocs.io/
* HST Data Handbook: https://hst-docs.stsci.edu/hstdhb

--------------------------------------------------------------------
7. License
--------------------------------------------------------------------

This documentation is released under the same license as the original
script.  See the source file for details.
