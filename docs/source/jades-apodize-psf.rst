=========================================================
epsf_apodizer – Apodize an Empirical PSF and Compute EE Curve
=========================================================

**Version:** |release|  
**Author:**  (original author)  
**License:** MIT (or as specified in the source)

.. contents:: Table of Contents
   :depth: 2
   :local:

--------------------------------------------------------------------
Overview
--------------------------------------------------------------------

``epsf_apodizer`` is a small command‑line utility written in pure Python
that

* reads an empirical point‑spread function (ePSF) stored in a FITS file,
* applies a *Tukey* (cosine‑tapered) apodisation window to the PSF,
* optionally normalises the resulting PSF,
* computes the **enclosed‑energy (EE) growth curve** as a function of radius,
* writes the apodised PSF back to a FITS file and the EE curve to plain
  text files, and
* can convert the radius axis from pixels to arcseconds using either the
  WCS information in the header or a user‑supplied pixel scale.

The script is intended for astronomers who need a smooth, finite‑support
PSF for forced photometry, image convolution or PSF modelling.

--------------------------------------------------------------------
Dependencies
--------------------------------------------------------------------

The script relies on the scientific Python ecosystem:

* ``numpy`` – array handling.
* ``astropy`` – FITS I/O, WCS handling and basic statistics.
* ``photutils`` – aperture photometry utilities.
* ``argparse`` – command‑line parsing (standard library).

All of these packages are available from ``pip`` or ``conda`.

--------------------------------------------------------------------
Command‑line Interface
--------------------------------------------------------------------

The parser is defined in :func:`create_parser`.  The following options are
available:

+----------------------+----------------------+---------------------------------------+
| Option               | Default              | Meaning                               |
+======================+======================+=======================================+
| ``-i`` / ``--input`` | ``epsf.fits``        | Path to the input FITS image that     |
|                      |                      | contains the empirical PSF.           |
+----------------------+----------------------+---------------------------------------+
| ``--nr``             | ``1000``             | Number of radial bins used when       |
|                      |                      | constructing the EE curve.            |
+----------------------+----------------------+---------------------------------------+
| ``-o`` / ``--output``| ``epsf.apodized.fits``| Filename for the apodised PSF output. |
+----------------------+----------------------+---------------------------------------+
| ``--normalize``      | ``False``            | If given, the apodised PSF is         |
|                      |                      | normalised to unit total flux.        |
+----------------------+----------------------+---------------------------------------+
| ``--psf_pixel_scale``| ``None``             | Floating‑point arcsec/pixel value to   |
|                      |                      | use when converting the EE radius to   |
|                      |                      | arcseconds.                            |
+----------------------+----------------------+---------------------------------------+
| ``--wcs_pixel_scale``| ``False``            | If set, the script extracts the pixel  |
|                      |                      | scale from the FITS WCS header and    |
|                      |                      | writes an EE curve in arcseconds.      |
+----------------------+----------------------+---------------------------------------+
| ``-v`` / ``--verbose``| ``False``           | Print additional progress information.|
+----------------------+----------------------+---------------------------------------+

Typical invocation::

    $ python epsf_apodizer.py -i my_epsf.fits -o my_epsf.apod.fits --normalize \
        --nr 2000 --wcs_pixel_scale -v

--------------------------------------------------------------------
Module Structure
--------------------------------------------------------------------

The source file is divided into three logical parts:

1. **Argument parsing** – :func:`create_parser`.
2. **Utility functions** – :func:`tukey` and :func:`create_enclosed_energy`.
3. **Main driver** – :func:`main`.

Below each component is explained in detail.

--------------------------------------------------------------------
1. Argument parser – ``create_parser()``
--------------------------------------------------------------------

.. autofunction:: create_parser
   :noindex:

The function builds an :class:`argparse.ArgumentParser` with the options
listed in the table above.  The ``store_true`` actions (`--normalize`,
`--wcs_pixel_scale`, `--verbose`) set the corresponding attribute to
``True`` when the flag appears on the command line.

--------------------------------------------------------------------
2. Tukey apodisation – ``tukey(data, a=0.1, lam=0.995, flag_norm=True)``
--------------------------------------------------------------------

The **Tukey window** is a cosine‑tapered top‑hat that smoothly reduces a
function to zero at a chosen radius.  In the context of an ePSF it
prevents sharp edges that would otherwise cause ringing artefacts in
convolutions.

Parameters
~~~~~~~~~~

* ``data`` – 2‑D ``numpy`` array containing the raw PSF.
* ``a`` – Fraction of the transition region (default ``0.1``).  Smaller
  values give a narrower cosine roll‑off.
* ``lam`` – Outer radius of the window as a *fraction* of the half‑size
  of the image (default ``0.995``).  ``lam`` is multiplied by the
  distance from the centre to obtain the absolute radius.
* ``flag_norm`` – If ``True`` the returned PSF is normalised to unit sum.

Algorithm
~~~~~~~~~

1. **Negative pixels are set to zero** – PSFs are strictly non‑negative.
2. Compute the centre coordinates ``(x0, y0)`` (half the array size).
3. Build two 2‑D coordinate grids ``x`` and ``y`` that span
   ``[-x0, +x0]`` and ``[-y0, +y0]`` respectively.
4. Compute the radial distance ``r = sqrt(x**2 + y**2)``.
5. Define the inner radius ``al = (1‑a) * lam`` – inside this region the
   filter equals 1.
6. For ``al <= r < lam`` the filter follows the Tukey transition function

   .. math::

      T(r) = \frac12\Bigl[1-\cos\bigl(\pi\frac{r}{a\,\lambda}
          -\frac{\pi}{a}\bigr)\Bigr].

7. Outside ``lam`` the filter is zero.
8. Multiply the original image by the filter to obtain the apodised PSF.
9. If ``flag_norm`` is ``True`` divide by the total sum so that the PSF
   integrates to one.

Return values
~~~~~~~~~~~~~

* ``r`` – 2‑D array of radial distances (pixels) – useful for diagnostics.
* ``tfilter`` – The pure Tukey window (same shape as ``data``).
* ``tfp`` – The apodised (and possibly normalised) PSF.

--------------------------------------------------------------------
3. Enclosed energy – ``create_enclosed_energy(fin, nr=1000, flag_norm=True)``
--------------------------------------------------------------------

The **enclosed‑energy (EE) curve** quantifies the fraction of total PSF
flux contained within a circular aperture of radius *r*.  It is often
referred to as the *growth curve*.

Parameters
~~~~~~~~~~

* ``fin`` – 2‑D array (the apodised PSF) from which to compute EE.
* ``nr`` – Number of radius samples (default ``1000``).
* ``flag_norm`` – Normalise the PSF before measuring EE (default ``True``).

Procedure
~~~~~~~~~

1. Clip negative values to zero.
2. Construct the same ``x``/``y`` coordinate grids as in :func:`tukey`.
3. Compute the radial map ``rf``.
4. Determine ``r_max`` – the distance from the centre to the furthest
   corner of the image.
5. Build a linearly spaced radius array ``rr`` ranging from
   ``0.5`` to ``r_max`` with ``nr`` points.
6. Zero out any pixels whose radius exceeds the half‑width of the image
   (they fall outside the valid region).
7. Normalise the image if required.
8. Create a list of :class:`photutils.aperture.CircularAperture`
   objects centred on the PSF with radii given by ``rr``.
9. Run :func:`photutils.aperture_photometry` to obtain the summed flux
   inside each aperture.
10. Store the results in ``ee`` – an array of the same length as ``rr``.

Return values
~~~~~~~~~~~~~

* ``rr`` – Radii (pixels) at which EE was evaluated.
* ``ee`` – Fractional enclosed energy (0 ≤ EE ≤ 1).

--------------------------------------------------------------------
4. Main driver – ``main()``
--------------------------------------------------------------------

The ``main`` function orchestrates the workflow:

1. **Timer start** – ``time.time()`` for simple performance reporting.
2. **Parse arguments** via :func:`create_parser`.
3. **Load the FITS image**:

   * If the file contains multiple HDUs, the script assumes a science
     extension named ``SCI``; otherwise it reads the primary HDU.
   * Data are cast to ``float32`` for memory efficiency.

4. **Apodise the PSF** using :func:`tukey`, honouring the ``--normalize``
   flag.
5. **Write the apodised PSF** to the output FITS file (overwrites if it
   already exists).
6. **Compute the EE curve** with :func:`create_enclosed_energy`.
7. **Save the EE curve** (pixel units) to a ``.EE.txt`` file – only the
   portion where ``EE < 1`` is written to avoid a trailing line with a
   value of exactly 1.
8. **Optional WCS conversion** – if ``--wcs_pixel_scale`` is set:

   * A :class:`astropy.wcs.WCS` object is built from the FITS header.
   * The pixel scale (arcsec / pixel) is derived from the determinant of
     the ``pixel_scale_matrix``.
   * The radius column is multiplied by this scale and saved to
     ``.EE.arcsec.wcs_pixel_scale.txt``.

9. **Optional user‑supplied pixel scale** – if ``--psf_pixel_scale`` is
   supplied, the same conversion is performed using that value and the
   output file name ends with ``.EE.arcsec.psf_pixel_scale.txt``.
10. **Timer end** – print total elapsed time.

The script is guarded by the usual ``if __name__ == "__main__":`` block,
so it can also be imported as a module without executing the workflow.

--------------------------------------------------------------------
5. Example Use‑Cases
--------------------------------------------------------------------

**a) Simple apodisation and normalisation**

.. code-block:: console

    $ python epsf_apodizer.py -i raw_epsf.fits -o raw_epsf.apod.fits --normalize

Result:

* ``raw_epsf.apod.fits`` – PSF whose sum of pixel values equals 1.
* ``raw_epsf.EE.txt`` – growth curve in pixel units.

**b) Growth curve in arcseconds using the FITS WCS**

.. code-block:: console

    $ python epsf_apodizer.py -i raw_epsf.fits -o raw_epsf.apod.fits \
        --wcs_pixel_scale --nr 1500 -v

Additional output:

* ``raw_epsf.EE.arcsec.wcs_pixel_scale.txt`` – radius (arcsec) vs. EE.

**c) Providing a custom pixel scale**

If the FITS header lacks a valid WCS (e.g. simulated data), you can
specify the scale directly:

.. code-block:: console

    $ python epsf_apodizer.py -i simulated_epsf.fits \
        --psf_pixel_scale 0.04 -o simulated_epsf.apod.fits

The file ``simulated_epsf.EE.arcsec.psf_pixel_scale.txt`` contains the
curve in arcseconds using a scale of 0.04″ pixel⁻¹.

--------------------------------------------------------------------
6. Extending the Code
--------------------------------------------------------------------

The module is deliberately lightweight, making it easy to adapt:

* **Different window functions** – replace :func:`tukey` with a Gaussian
  or a custom radial mask.
* **Alternative photometry** – use `photutils`'s ``RadialProfile`` or
  integrate the 2‑D array directly with ``np.cumsum`` for speed.
* **Batch processing** – import ``epsf_apodizer`` in a larger pipeline
  and call :func:`tukey`/``create_enclosed_energy`` programmatically.

--------------------------------------------------------------------
7. API Reference
--------------------------------------------------------------------

.. automodule:: epsf_apodizer
   :members:
   :undoc-members:
   :show-inheritance:

--------------------------------------------------------------------
8. License
--------------------------------------------------------------------

The code is released under the MIT License.  See the ``LICENSE`` file in
the repository for the full text.

--------------------------------------------------------------------
9. Changelog
--------------------------------------------------------------------

* **v1.0.0** – Initial public release.
* Add more detailed docstrings and RST documentation.

--------------------------------------------------------------------
10. Contact & Contributions
--------------------------------------------------------------------

For bug reports, feature requests or contributions, please open an issue
or pull request on the project's GitHub repository.
