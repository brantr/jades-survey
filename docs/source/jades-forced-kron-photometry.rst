.. _kron_photometry_script:

=========================================================
``jades-forced-kron-photometry.py``: Kron Photometry Script – Detailed Documentation
=========================================================

This document provides a line‑by‑line explanation of the Python script used to perform
aperture photometry on astronomical images (e.g. JWST/NIRCam, MIRI, HST).  
The script reads a FITS image, an input catalog, optional weight, exposure‑time and
segmentation maps, builds elliptical (Kron) apertures, applies aperture corrections,
computes uncertainties and writes the results to a new FITS table.

The documentation is written in reStructuredText (RST) format so it can be included
directly in a ReadTheDocs site.

-----------------------------------------------------------------
Table of Contents
-----------------------------------------------------------------

.. contents::
   :depth: 2
   :local:

-----------------------------------------------------------------
Imports
-----------------------------------------------------------------

The script imports a small set of scientific Python packages:

``python
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.wcs import WCS
import argparse
from tqdm import tqdm
from photutils.aperture import (EllipticalAperture, EllipticalAnnulus,
                               CircularAperture, aperture_photometry,
                               ApertureStats)
import time
```

* **numpy** – array handling.
* **astropy.io.fits** – reading/writing FITS files.
* **astropy.io.ascii** – reading the input catalog (plain‑text).
* **astropy.table.Table** – building the output table.
* **astropy.wcs.WCS** – world‑coordinate system handling.
* **argparse** – command‑line argument parsing.
* **tqdm** – progress bars for long loops.
* **photutils.aperture** – classes for elliptical/circular apertures and
  photometry utilities.
* **time** – simple wall‑clock timing.

-----------------------------------------------------------------
Utility Functions
-----------------------------------------------------------------

The script is organised into a set of helper functions that are called from
``main()``.  Each function is documented below.

.. _create_parser:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``create_parser()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Creates an :class:`argparse.ArgumentParser` that defines all command‑line
options accepted by the program.

**Behaviour**

* Returns an ``argparse.ArgumentParser`` instance.
* The parser defines the following options (the most important ones are
  highlighted):

  ``-a / --aper`` – name of the aperture used for Kron photometry (default
  ``KRON``).

  ``-i / --input`` – path to the primary science FITS image (default
  ``input.fits``).

  ``-c / --catalog`` – path to the ASCII source catalog (default
  ``catalog``).

  ``--wht`` – optional weight image; if omitted the script will read the
  weight extension from the primary input image.

  ``--texp`` – optional exposure‑time image; if omitted the program will read
  the ``EXP`` extension from the primary input image.

  ``--err`` – path to the RMS error image (default ``err.fits``).

  ``--wht`` – optional separate weight map.

  ``--texp`` – optional separate exposure‑time map.

  ``--segmap`` – optional segmentation map used to mask neighbouring sources.

  ``--bbox_*`` – bounding‑box columns are read from the catalog when a
  segmentation map is supplied.

  ``--bsub`` – flag that activates background subtraction using an elliptical
  annulus.

  ``--of`` – multiplier for the outer annulus when ``--bsub`` is set.

  ``--of`` – *outer factor* (default = 1.0) that controls the size of the
  background annulus.

  ``--quick`` – a shortcut flag used during development; when set the script
  skips many expensive steps (weight, exposure, RMS loading, etc.).

  ``--verbose`` – prints timing information.

* Boolean flags are implemented with ``action='store_true'`` where appropriate
  (e.g. ``--bsub``).

* The parser is returned to the caller; it is later used in ``main()`` to
  obtain the ``Namespace`` object that holds all user‑provided values.

-----------------------------------------------------------------

.. _get_instrument:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``get_instrument(args, flag_inst=None)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Selects the correct FITS extension (``SCI``, ``ERR``, ``WHT``, ``EXP``) based on
the instrument that produced the data (JWST NIRCam, MIRI or HST).

**Behaviour**

* The function is *not* explicitly defined in the script; instead, the logic for
  picking the appropriate extension lives inside :func:`get_image`.  The
  ``flag_inst`` argument passed to ``get_image`` tells the routine whether the
  image follows the JWST convention (separate ``SCI``/``ERR`` extensions) or the
  HST convention (single extension only).

* ``flag_inst`` is a string such as ``'NIRCAM'``, ``'MIRI'`` or ``'HST'``.
  It is obtained by calling :func:`get_instrument` (see below) and then passed
  to :func:`get_image` to decide which HDU name to request.

-----------------------------------------------------------------
.. _get_image:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``get_image(args, ext='SCI', flag_inst='NIRCAM')``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Read a specific image extension (science, error, weight, or exposure) from the
input FITS files, handling the differences between HST and JWST data formats.

**Parameters**

* ``args`` – the ``argparse.Namespace`` produced by ``create_parser``.
* ``ext`` – a string that selects which image to read:

  * ``'SCI'`` – the science data (default).
  * ``'ERR'`` – the error image (used for JWST; ignored for HST).
  * ``'WHT'`` – the weight map (or the same file if ``--wht`` is not supplied).
  * ``'EXP'`` – the exposure‑time map.

* ``flag_inst`` – instrument identifier (``'NIRCAM'``, ``'MIRI'``, ``'HST'``).

**Behaviour**

1. Determines which FITS file contains the requested extension:
   * For ``'WHT'`` and ``'EXP'`` the function checks whether the user supplied a
     separate file via ``--wht`` / ``--texp``; otherwise it falls back to the
     primary input image.
   * For ``'SCI'``, ``'ERR'``, and ``'RMS'`` the function reads the appropriate
     HDU name (``SCI`` for JWST, ``ERR`` for JWST error, ``RMS`` for the RMS
     catalog).

2. Uses :func:`astropy.io.fits.getdata` (or ``fits.open`` internally) to
   retrieve the image as a NumPy array.

3. Returns the image array unchanged (no scaling is applied here – scaling
   to physical units is performed later in ``main``).

-----------------------------------------------------------------
.. _compute_flux_to_nJy:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``compute_flux_to_nJy(args, flag_inst, header_flux)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Computes the multiplicative factor that converts the pixel values in the
science image (usually MJy sr⁻¹ for JWST) to nano‑Jansky (nJy).  The conversion
depends on the instrument and on header keywords that describe the photometric
calibration.

**Behaviour**

* Calls :func:`astropy.wcs.WCS` on the science image to obtain the pixel scale
  (arcsec pixel⁻¹) – used later for area calculations.
* For JWST instruments (NIRCam or MIRI) the conversion factor is taken from
  ``header_flux['PHOTMJSR']`` and the ``PHOTMJSR`` keyword; for HST the factor
  is derived directly from the header (the script references a non‑existent
  ``jpyn`` module – in practice this part of the code would raise an error unless
  the missing import is added; however the documentation reflects the *intended*
  behaviour).

* Returns a single floating‑point number ``flux_to_nJy`` that must be multiplied
  with the raw science data to obtain fluxes in nJy.

-----------------------------------------------------------------
.. _get_ellipsoidal_aperture_correction:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``get_ellipsoidal_aperture_correction(psf, apertures, args)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Derives an aperture‑correction factor for each source by measuring the fraction
of a point‑spread function (PSF) that falls inside the user‑defined aperture.
The script implements this by convolving the PSF with each elliptical aperture
and comparing the measured flux to the total PSF flux.

**Behaviour**

* Loads the PSF image (not shown in the posted code – the function expects a
  variable ``psf`` that must be supplied by the caller).
* For every source it creates an :class:`EllipticalAperture` using the same
  centre and semi‑major/minor axes as the Kron aperture.
* Calls :class:`photutils.aperture.ApertureStats` on the PSF image inside each
  aperture and records the ``sum`` value.
* The aperture‑correction factor is the ratio

  ``correction = total_PSF_flux / measured_PSF_flux``

  for each source.
* Returns a ``numpy.ndarray`` of correction factors with length equal to the
  number of sources.

*Note*: In the provided script the function is named
``get_ellipsoidal_aperture_correction`` but the actual implementation is
``get_ellipsoidal_aperture_correction`` – the description above follows the
actual code.

-----------------------------------------------------------------
.. _load_rms_uncertainties:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``load_rms_uncertainties(args, aperture_equiv_radius, nobj)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Estimate the root‑mean‑square (RMS) uncertainty for each source based on a
pre‑computed RMS image (generated by a separate pipeline).  The uncertainty
depends on the *linear* equivalent radius of the elliptical aperture.

**Parameters**

* ``args`` – the parsed command‑line arguments.
* ``aperture_area`` – *not* the full ellipse area but the linear equivalent
  radius (``sqrt(area)``) used to query the RMS table.
* ``nobj`` – total number of objects in the catalog.

**Behaviour**

1. Reads the RMS catalog (a FITS table) that contains columns ``AREA`` and
   ``RMS`` – these give the measured RMS for a range of aperture sizes.
2. Interpolates the RMS value for each source’s linear size using
   :func:`numpy.interp`.
3. Returns a ``numpy.ndarray`` (length ``nobj``) with the RMS noise expressed
   in the same units as the science image (later scaled to nJy).

-----------------------------------------------------------------
.. _main:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
``main()``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**

Orchestrates the whole photometric measurement pipeline:

* parses command‑line arguments,
* reads the science image, error/weight/exposure maps,
* builds Kron (elliptical) apertures from the input catalog,
* applies optional background subtraction,
* computes aperture corrections using a PSF,
* estimates uncertainties (RMS + Poisson),
* writes the final catalog to a FITS binary table.

**Detailed Flow**

1. **Argument parsing**

   ``args = create_parser().parse_args()``

2. **Timing start**

   ``time_start_global = time.time()``

3. **Read input catalog**

   ``cat = ascii.read(args.catalog)`` (the script expects columns such as
   ``id``, ``bbox_xmin`` … ``bbox_ymax``, ``bbox_*`` and Kron parameters).

4. **Instrument identification**

   ``instrument = get_instrument(args, None)`` – determines whether the data are
   from JWST/NIRCam, MIRI, or HST and stores the string in ``flag_inst``.

5. **Kron aperture preparation**

   * Compute semi‑major axis ``A`` and semi‑minor axis ``B`` from the catalog’s
     Kron radius using the user‑supplied Kron factor.
   * Apply a *minimum* aperture size (``args.kron_unit``) – if the calculated
     ``A`` or ``B`` fall below the minimum, they are reset to that value.
   * Build a list of :class:`EllipticalAperture` objects for each source:

     ``apertures = [EllipticalAperture(pos, A[i], B[i], theta[i]) ...]``

   * For weight and exposure maps a *circular* aperture of 5 pixel radius is
     created (``CircularAperture``) because those maps are typically sampled at
     lower resolution.

6. **Background subtraction (optional)** – if ``--bsub`` is set the script
   constructs an :class:`EllipticalAnnulus` for each source that defines the
   background region (inner and outer elliptical radii are derived from the
   ``--of`` factor).

7. **Weight map (WHT) and exposure map (EXP) handling**

   * ``data_wht = get_image(args, ext='WHT', flag_inst=instrument)``  
     The weight image is masked where its value is zero.
   * For each source the mean weight inside the *circular* 5‑pixel aperture is
     measured using :class:`ApertureStats`.  The result is stored in ``wht``.
   * The same procedure is repeated for the exposure‑time map, giving ``texp``.

8. **Segmentation map masking (optional)**

   If a segmentation map is supplied, the script temporarily *un‑masks* the
   pixels belonging to the current source while measuring its flux.  This is
   required because the weight map may have masked the source itself.

9. **Science image and error image**

   * ``data_flux = get_image(args, ext='SCI', flag_inst=instrument)``  
     The science data are multiplied by the conversion factor ``flux_to_nJy``
     obtained from :func:`compute_flux_to_nJy`.
   * For JWST instruments an error image (``ERR``) is also read and scaled
     to nJy.

10. **Photometry loop**

    For each source (loop over ``nobj`` with ``tqdm``):

    * If a segmentation map is present, the source’s own pixels are
      temporarily un‑masked.
    * ``ApertureStats`` (or ``aperture_photometry`` when no segmap is used)
      measures the summed flux inside the elliptical Kron aperture,
      the median pixel value, and the associated error.
    * The script corrects for masked pixels by estimating the missing flux
      from the median value and the *difference* between the theoretical
      aperture area and the actually measured area.
    * When ``--bsub`` is active, a background value (median of the annulus)
      is computed and later subtracted from the source flux.

11. **Aperture correction**

    * The PSF‑based correction factors are loaded with
      ``aperture_correction = get_ellipsoidal_aperture_correction(...)``.
    * If the user disabled corrections with ``--no-aper-corr`` the
      correction array is set to 1.

12. **Uncertainty estimation**

    * The RMS noise for the appropriate aperture size is read from the RMS
      catalog via :func:`load_rms_uncertainties`.
    * The conversion from nJy to electrons is derived from the exposure time
      and the photometric calibration keyword ``PHOTMJSR`` (for JWST).
    * Total uncertainty combines sky RMS and Poisson noise:

      ``total = sqrt( (sky_rms)^2 + |flux_in_counts| )``

    * The final uncertainty is converted back to nJy and multiplied by the
      aperture‑correction factor (if applicable).

13. **Populate output table**

    The script creates an :class:`astropy.table.Table` named ``t`` with the
    following columns (example for the *F200W* filter):

    * ``ID`` – source identifier (taken from the input catalog).
    * ``RA``, ``DEC`` – world coordinates (written via the WCS header).
    * ``F200W_KRON`` – aperture‑corrected flux in nJy.
    * ``F200W_KRON_e`` – total 1‑σ uncertainty (nJy) *including* the aperture
      correction.
    * ``F200W_KRON_ei`` – per‑pixel (or error‑image) uncertainty.
    * ``F200W_WHT`` – mean weight value inside a 5‑pixel circular aperture.
    * ``F200W_TEXP`` – exposure time (seconds) at the source position.
    * Optional columns ``*_bkg`` (background) are added when ``--bsub`` is used.

14. **Write FITS output**

    * The script builds a primary HDU (empty) and a binary table HDU named
      ``CAT`` that stores the table together with a header generated from the
      WCS and a number of provenance keywords (e.g. ``SCI``, ``WHT``, ``EXP``,
      ``RMS``, ``BAND``, ``APER`` etc.).
    * The resulting file is written to the path supplied by ``--output``
      (default ``output.fits``) with ``overwrite=True``.

15. **Timing report**

    If ``--verbose`` is set, the script prints the elapsed wall‑clock time for
    each major step and the total runtime.

-----------------------------------------------------------------
Function Reference
-----------------------------------------------------------------

Below each function is presented with its signature, a concise description of
its inputs and outputs, and any side effects.

.. function:: create_parser()

   Returns
   -------
   ``argparse.ArgumentParser``
       Fully populated parser ready for ``parse_args()``.

.. function:: get_instrument(args, flag_inst=None)

   *Not explicitly defined in the posted script.*  
   In the original version this function would inspect the FITS header
   (keywords such as ``INSTRUME``) to decide whether the data are from
   ``NIRCAM``, ``MIRI`` or ``HST``.  In the current script the logic is folded
   into ``main()`` where ``flag_inst`` is set manually.

.. function:: get_image(args, ext='SCI', flag_inst='NIRCAM')

   **Parameters**

   * ``args`` – the parsed command‑line arguments.
   * ``ext`` – one of ``'SCI'``, ``'ERR'``, ``'WHT'``, ``'EXP'``.
   * ``flag_inst`` – instrument identifier.

   **Returns**

   * ``numpy.ndarray`` containing the requested image data.

   **Details**

   * Chooses the correct FITS file based on ``args.wht``, ``args.texp``,
     ``args.err`` and the primary ``args.input``.
   * For HST images the weight and exposure extensions are not separate,
     so ``ext='WHT'`` and ``ext='EXP'`` fallback to the primary image.

.. function:: compute_flux_to_nJy(args, flag_inst, header_flux)

   **Purpose**

   Compute the multiplicative factor that converts the image units (normally
   MJy sr⁻¹ for JWST or the native units for HST) to nano‑Jansky.

   **Algorithm**

   1. Builds a WCS object from the science image header.
   2. Determines the pixel scale (arcsec pixel⁻¹) from the WCS.
   3. For JWST instruments reads the ``PHOTMJSR`` keyword (MJy sr⁻¹ per count) and
      combines it with the exposure time and ``PHOTFLAM``‑type conversion
      to produce ``flux_to_nJy``.
   4. Returns the scalar ``flux_to_nJy`` and the WCS object.

.. function:: get_ellipsoidal_aperture_correction(psf, apertures, args)

   **Purpose**

   Derive an aperture‑correction factor for each source by measuring how much of a
   known PSF falls inside the elliptical Kron aperture.

   **Steps**

   1. For each aperture, creates an :class:`EllipticalAperture` on the PSF image.
   2. Uses :class:`photutils.aperture.ApertureStats` to obtain the summed PSF
      flux inside the aperture.
   3. The correction factor is ``total_PSF_flux / measured_PSF_flux``.
   4. Returns a ``numpy.ndarray`` of correction factors (length = number of sources).

   *If the ``--no-aper-corr`` flag is set the function is bypassed and the
   correction array is set to 1.*

.. function:: load_rms_uncertainties(args, aperture_equiv_radius, nobj)

   **Purpose**

   Interpolate the RMS noise for a given aperture size from a pre‑computed RMS
   catalog (the script expects a FITS table where each row contains an ``AREA``
   column and a corresponding ``RMS`` column).

   **Parameters**

   * ``aperture_equiv_radius`` – ``sqrt(area)`` of the elliptical aperture
     (in pixels).
   * ``nobj`` – number of sources.

   **Return**

   * ``numpy.ndarray`` of RMS values for each source (nJy after scaling).

.. function:: main()

   See the *Detailed Flow* section above for a step‑by‑step description.
   This is the only function that is executed when the module is run as a
   script (``if __name__ == '__main__': main()``).

-----------------------------------------------------------------
Example Usage
-----------------------------------------------------------------

```bash
# Basic run on JWST NIRCam data
python krontable.py \
    --input jwst_image.fits \
    --err jwst_err.fits \
    --wht jwst_weight.fits \
    --texp jwst_exptime.fits \
    --catalog sources.txt \
    --kron 1.0 \
    --bsub \
    --of 1.5 \
    --output photometry.fits \
    --verbose
```

*The command reads the science image, associated weight/exp maps, computes
Kron apertures using a factor of 1.0, applies background subtraction with an
outer factor of 1.5, and writes the final catalog to ``photometry.fits``.*

-----------------------------------------------------------------
Glossary of Key Variables
-----------------------------------------------------------------

* ``args`` – Namespace with all user‑provided options.
* ``cat`` – Input source catalog (ASCII).
* ``instrument`` – String identifier of the instrument.
* ``flag_inst`` – Same as ``instrument``; used throughout the script.
* ``nobj`` – Number of entries in the catalog.
* ``A``, ``B`` – Semi‑major and semi‑minor axes of the elliptical Kron
  aperture (pixels).
* ``theta`` – Position angle of each ellipse (radians).
* ``apertures`` – List of :class:`EllipticalAperture` objects.
* ``wht`` – Mean weight value inside a 5‑pixel circular aperture.
* ``texp`` – Exposure time (seconds) at the source position.
* ``flux_to_nJy`` – Scalar conversion factor from image units to nJy.
* ``aperture_correction`` – Array of per‑source aperture‑correction factors.
* ``total_uncertainty`` – Final 1‑σ error (nJy) for each source.
* ``t`` – Output :class:`astropy.table.Table`.

-----------------------------------------------------------------
Conclusion
-----------------------------------------------------------------

The script provides a complete pipeline for measuring fluxes of extended
sources using Kron apertures, with optional background subtraction and PSF‑based
aperture corrections.  It is designed to work with data from JWST (NIRCam,
MIRI) and HST, handling the subtle differences in FITS structure between the
two observatories.  By separating the heavy‑lifting steps into dedicated
functions (:func:`create_parser`, :func:`get_image`, :func:`compute_flux_to_nJy`,
:func:`get_ellipsoidal_aperture_correction`, :func:`load_rms_uncertainties`),
the code remains modular and can be extended to support additional instruments
or calibration schemes.

-----------------------------------------------------------------
