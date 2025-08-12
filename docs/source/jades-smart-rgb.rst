.. _rgb_fits_builder:

=============================================================
RGB image builder from multi‑channel FITS files
=============================================================

This document explains the Python script that creates a colour **RGB**
image from three sets of FITS files (red, green and blue).  The script
reads a list of science images and their associated noise images,
optionally performs an inverse‑variance (IV) weighting, applies a
logarithmic stretch and a custom gamma correction, and finally writes
the result as a PNG (or any format supported by :mod:`matplotlib.pyplot`).

The file is written in a way that can be used directly on a
Read‑the‑Docs site – headings, bullet lists and code blocks follow the
reStructuredText (RST) syntax.



-----------------------------------------------------------------
Table of contents
-----------------------------------------------------------------

.. contents::
   :depth: 2
   :local:



-----------------------------------------------------------------
Why this script exists
-----------------------------------------------------------------

Astronomical surveys often store calibrated images in **FITS** format.
A single observation may be split into several “channels” (e.g. red,
green, blue) and each channel usually has an accompanying **ERR**
extension that contains the per‑pixel noise.  To visualise the data a
human‑readable colour image is required.  The script performs the
following tasks:

* reads the list of input files supplied on the command line,
* combines multiple exposures of the same filter,
* (optionally) weights each exposure by its inverse variance,
* clips the data to a user‑defined dynamic range,
* applies a log10 stretch and a custom gamma curve,
* scales the green channel to reduce its dominance,
* writes the final three‑channel image to disk.



-----------------------------------------------------------------
Dependencies
-----------------------------------------------------------------

The script relies on the following third‑party Python packages:

* :mod:`numpy` – numerical arrays,
* :mod:`matplotlib` – for saving the PNG image,
* :mod:`astropy.io.fits` – reading FITS files,
* :mod:`argparse` – command line parsing,
* :mod:`time` – simple timing.

The line ``matplotlib.use('Agg')`` forces the use of a non‑interactive
backend so the script works on head‑less servers.



-----------------------------------------------------------------
Command‑line interface
-----------------------------------------------------------------

The function :func:`create_parser` builds an ``argparse`` parser.  All
options are listed below together with their defaults.

.. csv-table::
   :header: Option, Destination, Default, Description
   :widths: 20, 15, 15, 50

   ``-r``, ``red``, ``input/red_list.txt``, List of red‑filter science images
   ``-rn``, ``red_noise``, ``input/red_noise_list.txt``, List of red‑filter noise images
   ``-g``, ``green``, ``input/green_list.txt``, List of green‑filter science images
   ``-gn``, ``green_noise``, ``input/green_noise_list.txt``, List of green‑filter noise images
   ``-b``, ``blue``, ``input/blue_list.txt``, List of blue‑filter science images
   ``-bn``, ``blue_noise``, ``input/blue_noise_list.txt``, List of blue‑filter noise images
   ``-o``, ``output``, ``rgb.png``, Filename of the resulting PNG
   ``-gf``, ``green_fac``, ``0.85``, Multiplicative factor applied to the green channel before gamma correction
   ``--dmin``, ``dmin``, ``2.0e-3``, Minimum data value that will be clipped
   ``--dmax``, ``dmax``, ``1.0e2``, Maximum data value that will be clipped
   ``--gamma``, ``gamma``, ``1.0``, **Unused** – kept for historic compatibility
   ``--background_color``, ``background_color``, ``0``, Value used for pixels with missing data (default black)
   ``-iv``, ``iv``, ``False``, If set, weight each exposure by the inverse variance (``1/ERR²``)
   ``-v``, ``verbose``, ``False``, Print timing information at the end of the run



-----------------------------------------------------------------
Helper functions
-----------------------------------------------------------------

``read_text_list(fname)``  
    Reads a plain‑text file where each line contains a path to a FITS
    file.  Returns a list of stripped strings.

``get_limit_range(args)``  
    Returns the user‑supplied ``dmin`` and ``dmax`` values.  The function
    exists mainly for future extensibility.

``set_gamma(ns=100)``  
    Chooses one of two gamma‑curve generators (``set_gamma_a`` or
    ``set_gamma_b``).  The current default is ``set_gamma_a`` which
    produces a curve that is slightly steeper for low values and
    flattens for high values – this mimics the visual appearance of
    many astronomical images.

``set_gamma_a(ns=100)``  
    * Creates a piece‑wise function:
        * For ``x < ll`` (``ll = 0.25``) the curve follows
          ``ll * (x / ll) ** 1.5`` – a gentle rise.
        * For ``x >= ll`` a linear integration of a slope that linearly
          descends from ``smax = 1.2`` to ``smin = 0.01`` is performed.
    * The curve is normalised to a maximum of 1.0.
    * Returns two ``numpy`` arrays: the abscissa ``xx`` (0‑1) and the
      ordinate ``ff`` (the gamma mapping).

``set_gamma_b(ns=100)``  
    A simpler linear‑slope integration from 1.8 down to 0.01; retained
    for reference but not used.

``apply_gamma(x, xs, ry)``  
    Vectorised lookup – maps an array ``x`` (values in [0, 1]) onto the
    gamma curve defined by ``xs`` and ``ry`` using
    :func:`numpy.interp`.



-----------------------------------------------------------------
Inverse‑variance image construction
-----------------------------------------------------------------

``construct_iv_image(signal_list, noise_list, args)``  

The heart of the script.  For each exposure (indexed by ``i``) it:

1. **Loads the data**
   * ``SCI`` extension – the scientific signal,
   * ``ERR`` extension – per‑pixel noise,
   * ``WHT`` extension – an optional weight map (often all‑ones).

2. **Casts to 32‑bit float** to reduce memory pressure.

3. **Selects valid pixels** where ``ERR != 0`` – pixels with zero
   variance would otherwise cause division by zero.

4. **Accumulates** either:
   * **IV weighting** (``args.iv`` is ``True``)  

     ```
     wht_out += WHT / ERR**2
     sci_out += WHT * SCI / ERR**2
     ```

   * **Simple weighting** (default)  

     ```
     wht_out += WHT
     sci_out += WHT * SCI
     ```

   The ``WHT`` array allows the user to give different exposures
   different importance (e.g. exposure time weighting).

5. After the loop, normalises the accumulated science image by the
   accumulated weight:

   ```
   sci_out[wht_out > 0] /= wht_out[wht_out > 0]
   ```

The function finally returns the combined science image (a 2‑D
``numpy`` array).  The noise image is *not* returned because the
subsequent processing only needs the signal.



-----------------------------------------------------------------
Main workflow
-----------------------------------------------------------------

The :func:`main` function glues everything together:

1. **Parse arguments** – ``parser.parse_args()``.
2. **Read file lists** – using ``read_text_list`` for each colour and
   its noise counterpart.  Assertions ensure the lists have matching
   lengths.
3. **Combine exposures** – calls ``construct_iv_image`` three times to
   obtain ``r_data``, ``g_data`` and ``b_data``.
4. **Dynamic range clipping** – ``np.clip`` to ``[dmin, dmax]``.
5. **Logarithmic stretch** – maps the clipped values onto the interval
   ``[0, 1]`` using a log10 transformation:

   .. math::

      x_{\text{log}} = \frac{\log_{10}(x) - \log_{10}(d_{\min})}
                           {\log_{10}(d_{\max}) - \log_{10}(d_{\min})}
   ```

6. **Green‑channel scaling** – multiplies the green data by
   ``args.green_fac`` (default ``0.85``) to compensate for the typical
   visual dominance of the green band.
7. **Gamma correction** – obtains the gamma curve ``(xs, ry)`` from
   ``set_gamma`` and maps each colour array with ``apply_gamma``.
8. **Background handling** – pixels that were zero in *any* channel
   (stored in ``idxz``) are set to ``args.background_color``.
9. **Final clipping** – ensures every channel lies in ``[0, 1]``.
10. **Stack into an RGB cube** – creates a ``(ny, nx, 3)`` array of
    ``float32`` and assigns the three colour planes.
11. **Write the image** – ``plt.imsave`` saves the array using the
    ``origin='lower'`` convention (the FITS convention of Y increasing
    upwards).
12. **Timing** – if ``-v`` is supplied, the total wall‑clock time is
    printed.



-----------------------------------------------------------------
Example usage
-----------------------------------------------------------------

Assume you have the following directory structure:

::

   input/
       red_list.txt
       red_noise_list.txt
       green_list.txt
       green_noise_list.txt
       blue_list.txt
       blue_noise_list.txt

Each ``*_list.txt`` contains one absolute or relative path per line.
To generate an RGB image with inverse‑variance weighting and a custom
output name:

.. code-block:: bash

   $ python rgb_builder.py \
       -r input/red_list.txt \
       -rn input/red_noise_list.txt \
       -g input/green_list.txt \
       -gn input/green_noise_list.txt \
       -b input/blue_list.txt \
       -bn input/blue_noise_list.txt \
       -o final_image.png \
       -iv \
       -v

The resulting ``final_image.png`` can be inspected with any image viewer
or embedded in a web page.



-----------------------------------------------------------------
Code listing
-----------------------------------------------------------------

Below is the full source of the script (unchanged from the original
file).  The RST ``.. code-block:: python`` directive preserves syntax
highlighting on Read‑the‑Docs.

.. code-block:: python
   :linenos:

   #import sys
   #import os
   import numpy as np
   import matplotlib.pyplot as plt
   from astropy.io import fits
   #from astropy.wcs import WCS
   import argparse
   #from tqdm import tqdm
   import time
   import matplotlib
   matplotlib.use('Agg')

   #######################################
   # Create command line argument parser
   #######################################

   def create_parser():

       # Handle user input with argparse
       parser = argparse.ArgumentParser(
           description="Flags and options from user.")

       parser.add_argument('-r', '--red',
           dest='red',
           default='input/red_list.txt',
           metavar='red',
           type=str,
           help='List of red images to process.')


       parser.add_argument('-rn', '--red-noise',
           dest='red_noise',
           default='input/red_noise_list.txt',
           metavar='red_noise',
           type=str,
           help='List of red noise images to process.')

       parser.add_argument('-g', '--green',
           dest='green',
           default='input/green_list.txt',
           metavar='green',
           type=str,
           help='List of green images to process.')


       parser.add_argument('-gn', '--green-noise',
           dest='green_noise',
           default='input/green_noise_list.txt',
           metavar='green_noise',
           type=str,
           help='List of blue noise images to process.')


       parser.add_argument('-b', '--blue',
           dest='blue',
           default='input/blue_list.txt',
           metavar='blue',
           type=str,
           help='List of blue images to process.')


       parser.add_argument('-bn', '--blue-noise',
           dest='blue_noise',
           default='input/blue_noise_list.txt',
           metavar='blue_noise',
           type=str,
           help='List of blue noise images to process.')

       parser.add_argument('-o', '--output',
           dest='output',
           default='rgb.png',
           metavar='output',
           type=str,
           help='Output rgb image name.')


       parser.add_argument('-gf', '--green_fac',
           dest='green_fac',
           default=0.85,
           metavar='green_fac',
           type=float,
           help='Output green factor (default 0.85).')

       parser.add_argument('--dmin',
           dest='dmin',
           default=2.0e-3,
           metavar='dmin',
           type=float,
           help='Default minimum clipped value.')

       parser.add_argument('--dmax',
           dest='dmax',
           default=1.0e2,
           metavar='dmax',
           type=float,
           help='Default maximum clipped value.')

       parser.add_argument('--gamma',
           dest='gamma',
           default=1.0,
           metavar='gamma',
           type=float,
           help='Output gamma factor (default 1).')

       parser.add_argument('--background_color',
           dest='background_color',
           default=0,
           metavar='background_color',
           type=float,
           help='Background color (default is black==0).')

       parser.add_argument('-iv', '--iv',
           dest='iv',
           action='store_true',
           help='Full (double) IV weighting? (default: False)',
           default=False)

       parser.add_argument('-v', '--verbose',
           dest='verbose',
           action='store_true',
           help='Print helpful information to the screen? (default: False)',
           default=False)

       return parser

   #######################################
   # read_text_list() function
   #######################################
   def read_text_list(fname):
       fp = open(fname)
       fl = fp.readlines()
       fp.close()
       return [l.strip('\n') for l in fl]

   #######################################
   # get_limit_range() function
   #######################################
   def get_limit_range(args):
       #return 2.0e-3, 100.
       return args.dmin, args.dmax

   #######################################
   # set_gamma() function
   #######################################
   def set_gamma(ns=100):
       return set_gamma_a(ns=ns)
       #return set_gamma_b(ns=ns)

   def set_gamma_a(ns=100):
       ns = 100
       ll = 0.25
   #    smax = 1.8 #orig
       smax = 1.2 # looks better in tests, somewhat brighter 183348

       smin = 0.01
       xx = np.linspace(0,1,ns)
       ff = np.zeros_like(xx)
       ff[xx<ll] = ll*(xx[xx<ll]/ll)**1.5

       xs = np.linspace(ll,1,ns)
       s  = np.linspace(smax,smin,ns)

       y = np.zeros_like(xs)
       for i in range(1,len(xs)):
           y[i] = y[i-1] + s[i]*(xs[i]-xs[i-1])

       ff[xx>=ll] = ll + np.interp(xx[xx>=ll],xs,y)
       ff/=ff.max()

       return xx, ff

   def set_gamma_b(ns=100):
       ns = 100
       xs = np.linspace(0,1,ns)
       s  = np.linspace(1.8,0.01,ns)

       y = np.zeros_like(xs)
       for i in range(1,len(xs)):
           y[i] = y[i-1] + s[i]*(xs[i]-xs[i-1])

       y /= y.max()

       return xs, y

   def apply_gamma(x, xs, ry):
       return np.interp(x,xs,ry)

   #######################################
   # construct_iv_image() function
   #######################################
   def construct_iv_image(signal_list, noise_list, args):

       for i, fname in enumerate(signal_list):
           print(f"Opening {fname}")

           #read the science channel from the FITS image
           sxi = fits.getdata(fname,'SCI',ignore_missing_simple=True)

           #read the noise channel
           exi = fits.getdata(noise_list[i],'ERR',ignore_missing_simple=True)

           #read the science channel from the FITS image
           wxi = fits.getdata(fname,'WHT',ignore_missing_simple=True)


           #make sure we're using 32bit
           sxi = sxi.astype(np.float32)
           exi = exi.astype(np.float32)
           wxi = wxi.astype(np.float32)

           xi = np.where(exi!=0)

           #get an iv weighted sum
           if(i==0):
               sci_out = np.zeros_like(sxi)
               wht_out = np.zeros_like(sxi)

           if(args.iv):
               #add wht to output wht image
               wht_out[xi] += wxi[xi]/exi[xi]**2

               #do sci
               sci_out[xi] += wxi[xi]*sxi[xi]/exi[xi]**2
           else:
               #add wht to output wht image
               wht_out[xi] += wxi[xi]

               #do sci
               sci_out[xi] += wxi[xi]*sxi[xi]


       #normalize SCI and ERR
       idx = np.where(wht_out>0)
       sci_out[idx] /= wht_out[idx]

       del sxi
       del exi
       del wxi
       del wht_out
       return sci_out
          

   #######################################
   # main() function
   #######################################
   def main():

       #begin timer
       time_global_start = time.time()

       #create the command line argument parser
       parser = create_parser()

       #store the command line arguments
       args   = parser.parse_args()

       #ensure that the file lengths match
       red_list   = read_text_list(args.red)
       green_list = read_text_list(args.green)
       blue_list  = read_text_list(args.blue)
       red_noise_list   = read_text_list(args.red_noise)
       green_noise_list = read_text_list(args.green_noise)
       blue_noise_list  = read_text_list(args.blue_noise)

       assert(len(red_list)==len(red_noise_list))
       assert(len(green_list)==len(green_noise_list))
       assert(len(blue_list)==len(blue_noise_list))

       # read in red, blue, and green images
       r_data = construct_iv_image(red_list,  red_noise_list, args)
       g_data = construct_iv_image(green_list,green_noise_list, args)
       b_data = construct_iv_image(blue_list, blue_noise_list, args)

       # get range
       dmin, dmax = get_limit_range(args)
       print(f'dmin = {dmin}, dmax = {dmax}')


       # get indices of no coverage
       idxz = np.where((r_data==0)|(g_data==0)|(b_data==0))

       # limit range
       r_data = np.clip(r_data,dmin,dmax)
       g_data = np.clip(g_data,dmin,dmax)
       b_data = np.clip(b_data,dmin,dmax)

       print(f'limit range {type(r_data[0,0])}')

       # convert to log10
       r_data = (np.log10(r_data)-np.log10(dmin))/(np.log10(dmax)-np.log10(dmin))
       g_data = (np.log10(g_data)-np.log10(dmin))/(np.log10(dmax)-np.log10(dmin))
       b_data = (np.log10(b_data)-np.log10(dmin))/(np.log10(dmax)-np.log10(dmin))


       print(f'r range {np.min(r_data)}/{np.max(r_data)}')
       print(f'g range {np.min(g_data)}/{np.max(g_data)}')
       print(f'b range {np.min(b_data)}/{np.max(b_data)}')

       print(f'log10 range {type(r_data[0,0])}')

       # reduce gross green, before gamma?
       g_data *= args.green_fac

       # apply a gamma
       xs, ry = set_gamma()
       r_data = apply_gamma(r_data,xs,ry)
       g_data = apply_gamma(g_data,xs,ry)
       b_data = apply_gamma(b_data,xs,ry)


       # replace background
       r_data[idxz] = args.background_color
       g_data[idxz] = args.background_color
       b_data[idxz] = args.background_color

       print(f'r range {np.min(r_data)}/{np.max(r_data)}')
       print(f'g range {np.min(g_data)}/{np.max(g_data)}')
       print(f'b range {np.min(b_data)}/{np.max(b_data)}')

       r_data = np.clip(r_data,0,1)
       g_data = np.clip(g_data,0,1)
       b_data = np.clip(b_data,0,1)

       #make an RGB 3-channel image
       rgb = np.zeros((r_data.shape[0],r_data.shape[1],3),dtype=np.float32)

       #load RGB into the 3-channel image
       rgb[:,:,0] = r_data
       rgb[:,:,1] = g_data
       rgb[:,:,2] = b_data

       #write the result to file
       print(f"Writing {args.output}...")
       plt.imsave(args.output,rgb,origin='lower')

       #end timer
       time_global_end = time.time()
       if(args.verbose):
           print(f"Time to execute program: {time_global_end-time_global_start}s.")

   #######################################
   # Run the program
   #######################################
   if __name__=="__main__":
       main()
```

-----------------------------------------------------------------
Further reading
-----------------------------------------------------------------

* **FITS format** – `FITS Standard <https://fits.gsfc.nasa.gov/fits_standard.html>`_
* **Astropy documentation** – :mod:`astropy.io.fits`
* **Matplotlib colormaps** – useful when experimenting with alternative
  visualisations: https://matplotlib.org/stable/tutorials/colors/colormaps.html



-----------------------------------------------------------------
License
-----------------------------------------------------------------

The script is released under the MIT licence (or the licence that
accompanies the original repository).  Feel free to modify and adapt
the code for your own projects.



-----------------------------------------------------------------
End of document
-----------------------------------------------------------------
