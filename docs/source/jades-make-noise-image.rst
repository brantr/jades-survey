.. _noise_image_generator:

========================================
Noise Image Generator – Detailed Overview
========================================

This document provides a thorough line‑by‑line explanation of the
``noise_image_generator.py`` script.  The script is a small command‑line
utility that reads a multi‑extension FITS file, builds a noise estimate
for the *ERR* (error) image, optionally smooths it, and writes the
result back to a new FITS file.  It is primarily used in astronomical
data reduction pipelines to replace unreliable error values with a
robust, locally‑averaged noise estimate.

.. contents::
   :local:
   :depth: 2


--------------------------------------------------------------------
1. Purpose
--------------------------------------------------------------------

The script performs three high‑level tasks:

1. **Parse user‑provided options** (input file, output file, smoothing
   parameters, etc.) via the :mod:`argparse` module.
2. **Compute a smoothed noise image** from the error extension of the
   input FITS file.  The smoothing can be either a median filter or a
   Gaussian convolution, depending on the ``--median`` flag.
3. **Replace low‑S/N error values** (those below a user‑defined
   fraction of the local noise) with the newly computed noise estimate
   and write the corrected error array to a new FITS file.

The script is deliberately lightweight: it only depends on
``numpy``, ``astropy`` and ``scipy`` and can be invoked from the shell
or embedded in a larger pipeline.


--------------------------------------------------------------------
2. Dependencies
--------------------------------------------------------------------

The script imports the following third‑party packages:

.. code-block:: python

   import argparse
   import numpy as np
   from astropy.io import fits
   from astropy.convolution import Gaussian2DKernel, convolve
   from astropy.stats import gaussian_fwhm_to_sigma
   from scipy.ndimage import median_filter
   import time

* **argparse** – standard library, handles command‑line arguments.
* **numpy** – numerical operations, array handling.
* **astropy.io.fits** – reading and writing FITS files.
* **astropy.convolution** – Gaussian kernel creation and convolution.
* **astropy.stats.gaussian_fwhm_to_sigma** – conversion from full‑width‑half‑maximum to σ.
* **scipy.ndimage.median_filter** – optional median smoothing.
* **time** – simple wall‑clock timing for benchmarking.


--------------------------------------------------------------------
3. Command‑Line Interface
--------------------------------------------------------------------

The function :func:`create_parser` builds an ``ArgumentParser`` instance
that defines the accepted options.  The most important arguments are
re‑listed below; the default values are shown in *italics*.

.. code-block:: python

   parser.add_argument('-i','--input',
                       default='input.fits',
                       type=str,
                       help='Input composite fits for creating a noise layer without artifacts')
   parser.add_argument('-o','--output',
                       default='output.fits',
                       type=str,
                       help='Output noise image.')
   parser.add_argument('-np_s', '--npix_smooth',
                       default=3.0,
                       type=float,
                       help='Number of pixels for Gaussian kernel FWHM.')
   parser.add_argument('-t', '--threshold',
                       default=0.9,
                       type=float,
                       help='Threshold for replacement (fraction of local noise).')
   parser.add_argument('-np_f', '--npix_footprint',
                       default=11,
                       type=int,
                       help='Number of pixels for kernel footprint (size of the convolution kernel).')
   parser.add_argument('-v', '--verbose',
                       action='store_true',
                       default=False,
                       help='Print helpful information to the screen?')
   parser.add_argument('-m', '--median',
                       action='store_true',
                       default=False,
                       help='Use median filtering?')

Running the script with ``-h`` prints a help page that mirrors the
descriptions above.


--------------------------------------------------------------------
4. Core Functions
--------------------------------------------------------------------

Only two functions exist in the file:

* :func:`create_parser` – builds the ``argparse`` parser (see section 3).
* :func:`main` – the entry point that implements the processing logic.


--------------------------------------------------------------------
5. The ``main`` Function – Step‑by‑Step Walk‑through
--------------------------------------------------------------------

The body of :func:`main` is annotated below.  For each logical block we
describe **what** is done, **why** it is necessary, and any subtle
behaviour.

.. code-block:: python

   def main():
       # ----------------------------------------------------------------
       # 5.1  Initialise timing
       # ----------------------------------------------------------------
       time_start = time.time()
``time_start`` records the wall‑clock time; later we report the total
runtime when ``--verbose`` is active.

   # ----------------------------------------------------------------
   # 5.2  Parse command line arguments
   # ----------------------------------------------------------------
   parser = create_parser()
   args   = parser.parse_args()
If ``--verbose`` is set we echo the chosen input and output filenames:

   if args.verbose:
       print(f"Input image = {args.input}")
       print(f"Output image = {args.output}")

   # ----------------------------------------------------------------
   # 5.3  Open the FITS file and extract relevant extensions
   # ----------------------------------------------------------------
   hdu = fits.open(args.input)

   # In the expected data product the error map lives in HDU[2] and the
   # weight map (inverse variance) lives in HDU[4].
   data_err = hdu[2].data           # ERR extension (float or double)
   data_wht = hdu[4].data           # WHT extension (usually 0/1)

   header_err = hdu[2].header       # Preserve original header for output

   # ----------------------------------------------------------------
   # 5.4  Build a mask of *bad* pixels
   # ----------------------------------------------------------------
   mask = np.zeros_like(data_err, dtype=bool)
   mask[data_wht == 0] = True               # weight = 0 → no data
   mask[data_err == 0] = True               # error = 0 is suspicious
   mask[np.isnan(data_err)] = True          # NaNs are always masked
``mask`` is later passed to ``convolve`` so that the kernel ignores those
pixels.

   # ----------------------------------------------------------------
   # 5.5  Determine kernel footprint size
   # ----------------------------------------------------------------
   np_f = args.npix_footprint                # odd integer, e.g. 11

   # ----------------------------------------------------------------
   # 5.6  Produce a smoothed noise estimate
   # ----------------------------------------------------------------
   if args.median:
       # 5.6.1 Median filtering path
       print(f"Median filtering with a {np_f}x{np_f} boxcar.")
       # ``data_err**2`` converts the error map to a variance map.
       # ``median_filter`` works on the variance; we take the sqrt later.
       data_noise = median_filter(data_err**2, size=(np_f, np_f))**0.5
   else:
       # 5.6.2 Gaussian smoothing path (default)
       # Convert user‑supplied FWHM (in pixels) to sigma.
       sigma = args.npix_smooth * gaussian_fwhm_to_sigma
       # Build a 2‑D Gaussian kernel with the requested footprint.
       kernel = Gaussian2DKernel(sigma, x_size=np_f, y_size=np_f)

       # Convolve the variance image (err²) with the kernel while
       # respecting the mask.  ``nan_treatment='interpolate'`` fills
       # NaNs by interpolating from neighbours.
       convolved = convolve(data_err**2,
                            kernel,
                            mask=mask,
                            nan_treatment='interpolate',
                            normalize_kernel=True)

       # Convert back to an error (standard deviation) map.
       data_noise = np.sqrt(convolved)

The two branches give the user control over the smoothing strategy.
Median filtering is robust against outliers, whereas Gaussian smoothing
preserves a more physically motivated point‑spread function.

   # ----------------------------------------------------------------
   # 5.7  Replace low‑signal error values
   # ----------------------------------------------------------------
   # ``idx`` contains the indices where the original error is *significantly
   # smaller* than the locally estimated noise.
   idx = np.where(
       (data_noise > 0) &
       (data_err > 0) &
       (~np.isnan(data_err)) &
       (data_err < args.threshold * data_noise)
   )

   if len(idx[0]) > 0:                     # ``np.where`` returns a tuple
       data_err[idx] = data_noise[idx]    # replace with the smoothed value

The ``threshold`` argument (default 0.9) is a multiplicative factor:
any error pixel below *threshold × local_noise* is considered
underestimated and is overwritten by the smoother estimate.

   # ----------------------------------------------------------------
   # 5.8  Write the corrected error image to disk
   # ----------------------------------------------------------------
   if args.verbose:
       print(f"Saving noise image {args.output}...")

   # Reduce memory footprint for the output FITS file.
   data_err = data_err.astype(np.float32)

   # Construct a minimal ImageHDU with the original header.
   hdu = fits.ImageHDU(data=data_err, header=header_err, name='ERR')
   hdu.writeto(args.output, overwrite=True)

   # ----------------------------------------------------------------
   # 5.9  Report timing (optional)
   # ----------------------------------------------------------------
   time_end = time.time()
   if args.verbose:
       print(f"Time to compute corrected error image = {time_end - time_start}s.")
```

The script finishes by exiting the ``main`` function.  The standard
``if __name__ == "__main__":`` guard ensures that the routine runs only
when the file is executed as a script, not when it is imported as a
module.


--------------------------------------------------------------------
6. Example Invocation
--------------------------------------------------------------------

Assuming the script is saved as ``noise_image_generator.py`` and you have
a FITS file ``mydata.fits`` with the expected extensions, a typical
call looks like:

.. code-block:: bash

   python noise_image_generator.py \
          --input  mydata.fits \
          --output denoised_err.fits \
          --npix_smooth 2.5 \
          --threshold 0.85 \
          --npix_footprint 9 \
          --verbose

The command will:

* read the error and weight extensions,
* apply a Gaussian smoothing with a 2.5‑pixel FWHM,
* replace any error pixel smaller than 85 % of the local noise,
* write the corrected error image to ``denoised_err.fits``,
* and print progress information.


--------------------------------------------------------------------
7. Design Considerations & Extensibility
--------------------------------------------------------------------

* **Mask handling** – the script masks both zero‑weight pixels and any
  NaNs before convolution.  This prevents the smoothing kernel from
  being biased by invalid data.
* **Variance‑based smoothing** – smoothing is performed on the *variance*
  (``err**2``) rather than directly on the error.  This is mathematically
  correct because variances add linearly under convolution.
* **Choice of smoothing** – the ``--median`` flag offers a non‑linear
  alternative useful when the image contains bright artefacts that would
  otherwise inflate the Gaussian‑smoothed variance.
* **Threshold flexibility** – users can tune the replacement aggressiveness
  without altering the core algorithm.
* **Future extensions** – adding support for additional FITS extensions,
  multi‑band processing, or writing a full multi‑HDU output (e.g. preserving
  the original data cube) would only require modest changes in the
  ``main`` function.

--------------------------------------------------------------------
8. Frequently Asked Questions
--------------------------------------------------------------------

**Q: Why is the error image stored in HDU[2]?**  
A: The script follows the convention used by the original data reduction
pipeline (e.g., *HSC‑Pipe* or *JWST* pipelines) where the primary HDU
contains the science image, HDU[1] may hold a variance or mask, and HDU[2]
stores the propagated error map.  Adjust the indices if your files differ.

**Q: Can I use a different kernel shape?**  
A: Yes. Replace the ``Gaussian2DKernel`` construction with any
`astropy.convolution` kernel (e.g. ``TopHat2DKernel``) and keep the rest of
the code unchanged.

**Q: My input file contains NaNs in the error plane – are they handled?**  
A: Absolutely. The mask marks NaNs as bad pixels and the convolution
uses ``nan_treatment='interpolate'`` to fill them from surrounding
valid values.

--------------------------------------------------------------------
9. References
--------------------------------------------------------------------

* Astropy documentation – `Convolution <https://docs.astropy.org/en/stable/convolution/>`_
* Astropy documentation – `FITS I/O <https://docs.astropy.org/en/stable/io/fits/>`_
* Scipy ndimage – `median_filter <https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.median_filter.html>`_
* Gaussian FWHM to sigma conversion – ``gaussian_fwhm_to_sigma = 1 / (2*sqrt(2*log(2)))``

--------------------------------------------------------------------
10. License
--------------------------------------------------------------------

The script is released under the BSD‑3‑Clause license (or the license
chosen by your project).  See the accompanying ``LICENSE`` file for the
full text.
