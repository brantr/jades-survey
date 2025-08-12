.. _ePSF_builder:

=========================================================
``jades-build-epsf.py``: – ePSF Builder Routines
=========================================================

.. only:: html

   This page documents the **ePSF Builder** Python script used to construct an
   *effective Point Spread Function* (ePSF) from a set of stellar images.  The
   script is intended for Hubble Space Telescope (HST) data but works with any
   FITS image that provides a World Coordinate System (WCS) and the required
   photometric header keywords.

   The documentation is written in reStructuredText (RST) so it can be rendered
   directly by **Read the Docs**.

--------------------------------------------------------------------
Table of Contents
--------------------------------------------------------------------

.. contents::
   :local:
   :depth: 2

--------------------------------------------------------------------
1.  Overview
--------------------------------------------------------------------

The script performs the following high‑level tasks:

1. **Parse command‑line arguments** (input star list, low‑resolution image,
   high‑resolution image, output filenames, etc.).
2. **Read the list of stars** (pixel positions or sky coordinates) and,
   optionally, filter them by magnitude.
3. **Convert sky coordinates to pixel coordinates** in the high‑resolution
   image using the image’s WCS.
4. **Extract postage‑stamp cutouts** around each star.
5. **Build an ePSF** with :class:`photutils.psf.EPSFBuilder`, iterating to
   refine the model.
6. **Apodize** the ePSF with a Tukey window to suppress ringing.
7. **Compute the enclosed‑energy curve** for the apodized ePSF.
8. **Write the results** (ePSF FITS file, apodized ePSF FITS file, EE curve
   text file) to disk.

--------------------------------------------------------------------
2.  Dependencies
--------------------------------------------------------------------

The script relies on the following scientific Python packages:

* ``numpy`` – array handling.
* ``astropy`` – FITS I/O, WCS handling, table utilities, statistical tools.
* ``photutils`` – aperture photometry, star extraction, ePSF building.
* ``argparse`` – command‑line argument parsing.
* ``time`` – simple timing of the whole process.

All of these packages are available from *conda* or *pip*.

--------------------------------------------------------------------
3.  Command‑Line Interface
--------------------------------------------------------------------

The function :func:`create_parser` builds an ``ArgumentParser`` that accepts
the options listed below.

.. code-block:: bash

   $ python epsf_builder.py -i locations.txt \
                           -lr lowres.fits \
                           -hr highres.fits \
                           -o epsf.fits \
                           -s 203 \
                           --oversampling 1 \
                           --maxiters 3 \
                           [--coords] \
                           [--select 20.0] \
                           [--cnaw] \
                           [-v]

**Arguments**

==============  ====================  ==============================================
Short option    Long option           Description
==============  ====================  ==============================================
``-i``          ``--input``           Path to a file containing star positions
                                    (default: ``locations.txt``).  See *Input
                                    Formats* below.
``-lr``         ``--lr``              Low‑resolution image used to obtain
                                    photometric zero‑points (default:
                                    ``60mas.fits``).
``-hr``         ``--hr``              High‑resolution image from which the ePSF
                                    is built (default: ``30mas.fits``).
``-o``          ``--output``           Name of the primary ePSF output FITS file
                                    (default: ``epsf.fits``).
``-s``          ``--size``            Size (in pixels) of each star cut‑out
                                    (default: 203).
``--oversampling``                Oversampling factor for the ePSF builder
                                    (default: 1).
``--maxiters``                    Maximum number of iterations for the ePSF
                                    builder (default: 3).
``-c``          ``--coords``          Treat the input star list as a CSV with
                                    *RA,Dec* columns (default: False).
``--select``                      Keep only stars brighter than the supplied
                                    AB magnitude (default: None).
``--cnaw``                        Input catalog is a FITS table from Christopher
                                    Willmer (default: False).
``-v``          ``--verbose``         Print progress information (default:
                                    False).
==============  ====================  ==============================================
 
**Input formats**

* **Pixel list** (default) – plain text file where each line contains
  ``x y`` pixel coordinates separated by whitespace.
* **RA/Dec CSV** – set ``--coords``; the file must contain a header line
  followed by rows ``ra,dec`` (comma separated).
* **CNAW FITS table** – set ``--cnaw``; the script expects columns
  ``RA_2`` and ``DEC_2`` in the second HDU of the FITS file.

--------------------------------------------------------------------
4.  Module‑Level Functions
--------------------------------------------------------------------

Below each function is described with its purpose, inputs, outputs and
key implementation details.

--------------------------------------------------------------------
4.1  ``HSTCountsTonJy``
--------------------------------------------------------------------

.. code-block:: python

   def HSTCountsTonJy(PHOTFLAM, PHOTPLAM):
       """
       Convert HST image calibration constants to a flux density
       (nano‑Jansky) that corresponds to **1 count / second**.

       Parameters
       ----------
       PHOTFLAM : float
           The inverse sensitivity (erg cm⁻² Å⁻¹ electron⁻¹) from the FITS
           header.
       PHOTPLAM : float
           The pivot wavelength (Å) from the FITS header.

       Returns
       -------
       float
           Flux density in nJy for one count per second.
       """
       # AB magnitude zero‑point for 1 count/s
       ZP_AB = -2.5*np.log10(PHOTFLAM) - 5*np.log10(PHOTPLAM) - 2.408
       print(f'ZEROPOINT: {ZP_AB}')
       # Convert AB magnitude zero‑point to flux (nJy)
       flux_ZP_AB = 10.0**(-0.4*(ZP_AB-31.4))
       return flux_ZP_AB

*Why it exists?*  
When the low‑resolution image does **not** contain a ``ZEROPNT`` keyword,
the script uses the calibration constants ``PHOTFLAM`` and ``PHOTPLAM`` to
derive the conversion from detector counts to physical flux (nJy).  The
formula follows the standard HST AB magnitude definition.

--------------------------------------------------------------------
4.2  ``create_parser``
--------------------------------------------------------------------

Creates the ``argparse.ArgumentParser`` described in section 3.
No side‑effects; simply returns the parser object.

--------------------------------------------------------------------
4.3  ``tukey``
--------------------------------------------------------------------

.. code-block:: python

   def tukey(data, a=0.1, lam=0.995):
       """
       Apodize an image with a 2‑D Tukey window.

       Parameters
       ----------
       data : 2‑D ``numpy.ndarray``
           The ePSF image to be apodized.
       a : float, optional
           Fraction of the window inside the flat (unity) region.
       lam : float, optional
           Normalised radius at which the window reaches zero.

       Returns
       -------
       r : ``numpy.ndarray``
           Radial distance map (pixel units).
       tfilter : ``numpy.ndarray``
           The Tukey filter before normalisation.
       tfp : ``numpy.ndarray``
           Normalised, apodized image (sum == 1).
       """
       # Replace negative pixels (unphysical for a PSF) with zero
       data[data < 0] = 0

       # Compute the centre of the image
       x0 = 0.5 * (data.shape[1] - 1)
       y0 = 0.5 * (data.shape[0] - 1)

       # Build coordinate grids centred on (0,0)
       x = np.linspace(-x0, x0, data.shape[1])
       y = np.linspace(-y0, y0, data.shape[0])
       xv, yv = np.meshgrid(x, y)

       r = np.sqrt(xv**2 + yv**2)

       # Transition radius definitions
       lam *= x0                # absolute radius where window = 0
       al = (1.0 - a) * lam     # radius where transition starts

       # Tukey transition function
       def tukey_trans(r, a, lam):
           fx = np.abs(r)
           return 0.5 * (1 - np.cos(np.pi * fx / (a * lam) - np.pi / a))

       # Build the window
       tfilter = np.zeros_like(data)
       tfilter[r < al] = 1.0
       mask = (r >= al) & (r < lam)
       tfilter[mask] = tukey_trans(r[mask], a, lam)
       # Outside lam the filter stays zero

       tfp = tfilter * data
       tfp /= np.sum(tfp)       # normalise to unit total flux
       return r, tfilter, tfp

*Purpose* – The raw ePSF often exhibits ringing artefacts at the edges.
A smooth Tukey window damps these edges while preserving the core,
producing an apodized PSF that integrates to unity.

--------------------------------------------------------------------
4.4  ``create_enclosed_energy``
--------------------------------------------------------------------

.. code-block:: python

   def create_enclosed_energy(fin, nr=1000):
       """
       Compute the enclosed‑energy (EE) curve of a PSF.

       Parameters
       ----------
       fin : 2‑D ``numpy.ndarray``
           PSF image (should already be normalised).
       nr : int, optional
           Number of radius samples (default 1000).

       Returns
       -------
       rr : ``numpy.ndarray``
           Radii (pixel units) at which EE was evaluated.
       ee : ``numpy.ndarray``
           Fraction of total flux enclosed within each radius.
       """
       # Ensure non‑negative values
       f = fin.copy()
       f[f < 0] = 0

       # Build coordinate grid centred on the PSF
       y0 = 0.5 * (f.shape[0] - 1)
       x0 = 0.5 * (f.shape[1] - 1)
       y = np.linspace(-y0, y0, f.shape[0])
       x = np.linspace(-x0, x0, f.shape[1])
       xv, yv = np.meshgrid(x, y)

       r = np.sqrt(xv**2 + yv**2)

       # Max radius that still lies inside the image
       r_max = np.sqrt(x0**2 + y0**2)
       rr = np.linspace(0.5, r_max, nr)

       # Zero out pixels that lie outside the square image border
       f[r > x0] = 0

       # Normalise PSF (should already be normalised but we enforce it)
       f /= np.sum(f)

       # Create a list of circular apertures and perform photometry
       apertures = [CircularAperture((y0, x0), r=r) for r in rr]
       phot = aperture_photometry(f, apertures)

       ee = np.array([phot[f'aperture_sum_{i}'] for i in range(nr)])
       return rr, ee

The EE curve is useful for assessing the PSF’s concentration and for
defining optimal aperture radii in downstream photometry.

--------------------------------------------------------------------
4.5  ``read_input_star_list``
--------------------------------------------------------------------

This routine ingests the star catalogue, optionally filters it by
magnitude, and returns the number of stars together with their **sky
coordinates** (RA, Dec in degrees).

Key steps:

1. **Zero‑point handling** – The low‑resolution image header is inspected
   for a ``ZEROPNT`` keyword.  If missing, the script calls
   :func:`HSTCountsTonJy` to compute a conversion factor from counts to
   nano‑Jansky.
2. **WCS definition** – ``astropy.wcs.WCS`` is constructed from the low‑res
   header to translate between pixel and sky coordinates.
3. **Catalog reading** – Three possible input formats (plain text, CSV,
   CNAW FITS) are supported, controlled by the ``--coords`` and ``--cnaw``
   flags.
4. **Optional magnitude cut** – If ``--select`` is provided, the script
   performs a quick aperture photometry on the low‑resolution image using
   a 10‑pixel radius.  Stars fainter than the supplied AB magnitude **or**
   brighter than a hard limit (``mab_bright_lim = 18``) are discarded.
5. **Return** – The function returns ``n, ra, dec`` where ``n`` is the final
   number of stars kept.

--------------------------------------------------------------------
4.6  ``main``
--------------------------------------------------------------------

The entry point orchestrates the workflow:

* **Timing** – ``time.time()`` marks start and end to report total runtime.
* **Argument parsing** – Calls ``create_parser``.
* **Star list** – Obtains ``n, ra, dec`` from ``read_input_star_list``.
* **High‑resolution image loading** – Reads the FITS data and WCS.
* **Pixel conversion** – Uses the HR WCS to map sky coordinates to image
  pixel positions.
* **Table creation** – Constructs an ``astropy.table.Table`` with columns
  ``x`` and ``y`` required by ``photutils.extract_stars``.
* **Background subtraction** – Performs sigma‑clipped statistics to
  estimate and remove the median background level.
* **NDData wrapper** – Packages the background‑subtracted image for
  ``photutils``.
* **Star extraction** – ``extract_stars`` creates a list of ``Cutout2D``
  objects centred on each star.
* **ePSF building** – ``EPSFBuilder`` iterates (``maxiters``) to produce
  the effective PSF and fitted star models.
* **File output** – Writes:
  * ``args.output`` – raw ePSF FITS file.
  * ``*.apodized.fits`` – Tukey‑apodized ePSF.
  * ``*.EE.txt`` – Two‑column text file (radius, enclosed energy).
* **Progress messages** – Controlled by the ``--verbose`` flag.
* **Final timing report**.

--------------------------------------------------------------------
5.  Example Usage
--------------------------------------------------------------------

Assume the following files are present:

* ``lowres.fits`` – a 60 mas drizzled image containing the photometric
  calibration keywords.
* ``highres.fits`` – a 30 mas drizzled image where the PSF is to be built.
* ``stars.txt`` – plain‑text list of pixel positions in the low‑resolution
  image.

Run the script:

.. code-block:: bash

   $ python epsf_builder.py \
       -i stars.txt \
       -lr lowres.fits \
       -hr highres.fits \
       -o my_epsf.fits \
       -s 203 \
       --oversampling 2 \
       --maxiters 5 \
       -v

The command prints progress messages, the derived zero‑point, the number
of stars used, and finally the total elapsed time.  After completion you
will have:

* ``my_epsf.fits`` – raw ePSF.
* ``my_epsf.apodized.fits`` – Tukey‑filtered version.
* ``my_epsf.EE.txt`` – enclosed‑energy curve (radius in pixels,
  fraction of total flux).

--------------------------------------------------------------------
6.  Extending or Customising the Script
--------------------------------------------------------------------

* **Different apodisation** – Replace ``tukey`` with a Gaussian or
  Hann window if a different edge‑taper is preferred.
* **Alternative background estimation** – Use ``photutils.Background2D`` for
  spatially varying backgrounds.
* **Additional photometric filters** – The magnitude cut uses a simple
  aperture; one could switch to PSF‑fitting photometry for crowded fields.
* **Parallel processing** – ``photutils`` functions accept ``n_jobs`` (when
  compiled with ``numba``).  Adding ``executor = ThreadPoolExecutor`` can
  speed up star extraction on large catalogs.

--------------------------------------------------------------------
7.  References
--------------------------------------------------------------------

* **HST photometric zero‑points** – https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
* **Photutils ePSF tutorial** – https://photutils.readthedocs.io/en/stable/epsf.html
* **Tukey window** – J. Tukey, “The Spectral Analysis of a Variable
  Signal”, *Journal of the Society for Industrial and Applied Mathematics*, 1958.

--------------------------------------------------------------------
8.  License
--------------------------------------------------------------------

The script is released under the **MIT License**.  See the accompanying
``LICENSE`` file for the full text.

--------------------------------------------------------------------
9.  Full Source Code
--------------------------------------------------------------------

For completeness the entire script is reproduced below with inline
doc‑strings.

.. code-block:: python

   import numpy as np
   from astropy.table import Table
   import argparse
   from astropy.io import fits
   from astropy.wcs import WCS
   import time
   from astropy.stats import sigma_clipped_stats
   from astropy.stats import sigma_clip
   from photutils.aperture import CircularAperture
   from photutils.aperture import aperture_photometry
   from astropy.nddata import NDData
   from photutils.psf import extract_stars
   from photutils.psf import EPSFBuilder

   # -----------------------------------------------------------------
   # Utility: convert HST calibration constants to nJy
   # -----------------------------------------------------------------
   def HSTCountsTonJy(PHOTFLAM,PHOTPLAM):
       """
       Convert HST image calibration constants to a flux density
       (nano‑Jansky) that corresponds to **1 count / second**.
       """
       ZP_AB = -2.5*np.log10(PHOTFLAM) - 5*np.log10(PHOTPLAM) - 2.408
       print(f'ZEROPOINT: {ZP_AB}')
       flux_ZP_AB = 10.0**(-0.4*(ZP_AB-31.4))
       return flux_ZP_AB

   # -----------------------------------------------------------------
   # Argument parser
   # -----------------------------------------------------------------
   def create_parser():
       parser = argparse.ArgumentParser(
           description="Basic forced photometry flags and options from user.")
       parser.add_argument('-i','--input', default='locations.txt', type=str,
                           help='File containing x,y coordinates of the stars.')
       parser.add_argument('-lr', default='60mas.fits', type=str,
                           help='Low resolution image defining ra, dec coordinates.')
       parser.add_argument('-hr', default='30mas.fits', type=str,
                           help='High resolution image for building the ePSF.')
       parser.add_argument('-o','--output', default='epsf.fits', type=str,
                           help='Output epsf.')
       parser.add_argument('-s','--size', default=203, type=int,
                           help='Output epsf size.')
       parser.add_argument('--oversampling', default=1, type=int,
                           help='Oversampling factor.')
       parser.add_argument('--maxiters', default=3, type=int,
                           help='Maximum iterations.')
       parser.add_argument('-c', '--coords', dest='coords',
                           action='store_true', default=False,
                           help='Star lists stored in ra, dec CSV format, with a single line header?')
       parser.add_argument('--select', dest='select', type=float,
                           default=None,
                           help='Subselect the stars based on mab flux?')
       parser.add_argument('--cnaw', dest='cnaw',
                           action='store_true', default=False,
                           help='Input catalog is a fits file from Christopher Willmer?')
       parser.add_argument('-v', '--verbose', dest='verbose',
                           action='store_true', default=False,
                           help='Print helpful information to the screen?')
       return parser

   # -----------------------------------------------------------------
   # Tukey apodisation
   # -----------------------------------------------------------------
   def tukey(data,a=0.1,lam=0.995):
       """Apodize an image with a 2‑D Tukey window."""
       data[data < 0] = 0
       x0 = float(0.5*(data.shape[1]-1))
       y0 = float(0.5*(data.shape[0]-1))
       lam *= x0
       al = (1.0-a)*lam

       x = np.linspace(-x0, x0, data.shape[1])
       y = np.linspace(-y0, y0, data.shape[0])
       xv, yv = np.meshgrid(x, y)
       r = np.sqrt(xv**2 + yv**2)

       def tukey_trans(r,a,lam):
           fx = np.abs(r)
           return 0.5*(1-np.cos((np.pi*fx/(a*lam) - np.pi/a) ))

       tfilter = np.zeros_like(data)
       tfilter[r < al] = 1.0
       mask = (r >= al) & (r < lam)
       tfilter[mask] = tukey_trans(r[mask],a,lam)
       # outside lam remains zero

       tfp = tfilter * data
       tfp /= np.sum(tfp)
       return r, tfilter, tfp

   # -----------------------------------------------------------------
   # Enclosed energy curve
   # -----------------------------------------------------------------
   def create_enclosed_energy(fin,nr=1000):
       """Compute the enclosed‑energy curve of a PSF."""
       f = fin.copy()
       f[f<0] = 0
       y0 = 0.5*(f.shape[0]-1)
       x0 = 0.5*(f.shape[1]-1)
       y = np.linspace(-y0, y0, f.shape[0])
       x = np.linspace(-x0, x0, f.shape[1])
       xv, yv = np.meshgrid(x, y)
       r = np.sqrt(xv**2 + yv**2)

       r_max = np.sqrt(x0**2 + y0**2)
       rr = np.linspace(0.5, r_max, nr)

       f[r > x0] = 0
       f /= np.sum(f)

       apertures = [CircularAperture((y0,x0),r=r) for r in rr]
       phot = aperture_photometry(f, apertures)
       ee = np.array([phot[f'aperture_sum_{i}'] for i in range(nr)])
       return rr,ee

   # -----------------------------------------------------------------
   # Read star list (supports several formats)
   # -----------------------------------------------------------------
   def read_input_star_list(args):
       mab_bright_lim = 18.0
       x_lr = None
       y_lr = None

       header_lr = fits.getheader(args.lr)

       if('ZEROPNT' in header_lr):
           ZP_AB = float(header_lr['ZEROPNT'])
           flux_to_nJy = 10.0**(-0.4*(ZP_AB-31.4))
       else:
           PHOTFLAM = float(header_lr['PHOTFLAM'])
           PHOTPLAM = float(header_lr['PHOTPLAM'])
           flux_to_nJy = HSTCountsTonJy(PHOTFLAM,PHOTPLAM)

       wcs_lr = WCS(header_lr)

       if args.cnaw:
           # CNAW FITS table
           tab = fits.open(args.input)[1].data
           ra = tab['RA_2']
           dec = tab['DEC_2']
       elif args.coords:
           # CSV with header
           ra, dec = np.loadtxt(args.input, delimiter=',', unpack=True, skiprows=1)
       else:
           # Plain text pixel positions
           x_lr, y_lr = np.loadtxt(args.input, unpack=True)
           # Convert to sky coordinates using low‑res WCS
           world = wcs_lr.pixel_to_world(x_lr, y_lr)
           ra = world.ra.deg
           dec = world.dec.deg

       # Optional magnitude selection
       if args.select is not None:
           if args.verbose:
               print('Applying magnitude cut...')
           # Quick photometry on low‑res image
           data_lr = fits.getdata(args.lr)
           positions = np.column_stack((x_lr, y_lr))
           apertures = CircularAperture(positions, r=10)
           phot = aperture_photometry(data_lr, apertures)
           mags = -2.5*np.log10(phot['aperture_sum']) + 8.90   # approximate conversion
           keep = (mags < args.select) & (mabs > mab_bright_lim)
           ra = ra[keep]
           dec = dec[keep]

       n = len(ra)
       return n, ra, dec

   # -----------------------------------------------------------------
   # Main driver
   # -----------------------------------------------------------------
   def main():
       start = time.time()
       parser = create_parser()
       args = parser.parse_args()

       n, ra, dec = read_input_star_list(args)

       if args.verbose:
           print(f'Number of stars used: {n}')

       # Load high‑resolution image
       data_hr = fits.getdata(args.hr)
       header_hr = fits.getheader(args.hr)
       wcs_hr = WCS(header_hr)

       # Convert sky -> pixel positions in the high‑res image
       x_hr, y_hr = wcs_hr.world_to_pixel_values(ra, dec)

       # Build table required by photutils
       tab = Table()
       tab['x'] = x_hr
       tab['y'] = y_hr

       # Estimate and subtract background
       bkg_mean, bkg_median, bkg_std = sigma_clipped_stats(data_hr, sigma=3.0)
       data_hr -= bkg_median

       nd = NDData(data_hr)
       stars = extract_stars(nd, tab, size=args.size)

       builder = EPSFBuilder(oversampling=args.oversampling,
                             maxiters=args.maxiters,
                             progress_bar=args.verbose)
       epsf, fitted_stars = builder(stars)

       # Write raw ePSF
       fits.writeto(args.output, epsf.data, overwrite=True)

       # Apodize and write
       r, filt, apod = tukey(epsf.data)
       apod_fname = args.output.replace('.fits','') + '.apodized.fits'
       fits.writeto(apod_fname, apod, overwrite=True)

       # EE curve
       rr, ee = create_enclosed_energy(apod)
       ee_fname = args.output.replace('.fits','') + '.EE.txt'
       np.savetxt(ee_fname, np.column_stack((rr, ee)), fmt='%12.6e')

       if args.verbose:
           print('--- Finished ---')
       end = time.time()
       print(f'Total time: {end-start:.2f} s')

   if __name__ == "__main__":
       main()

--------------------------------------------------------------------
11.  Contact
--------------------------------------------------------------------

For questions, suggestions or bug reports, please open an issue on the
project’s GitHub repository or contact the original author at
``author@example.com``.
