======================================
``jades-psf-match-lite-split.py``: Common PSF Image Generator
======================================
**In plain English – what this script actually does**

1. **What it’s for**  
   Astronomers take pictures of the sky with space‑telescopes (HST, JWST NIRCam, JWST MIRI, etc.).  
   Each picture has a different amount of “blur” (the point‑spread‑function or PSF) that comes from the optics of the instrument.  
   If you want to compare two images – for example a Hubble image and a JWST image – you first need to make their blurs the same.  
   This program reads an image, smooths it so that its blur matches a *target* blur you supply, and (for JWST data) also produces a matching error‑map.

2. **How you run it**  
   You call the script from the command line and tell it:  

   * the science image to be processed (`‑i …`)  
   * the image that contains the current PSF (`‑c …`)  
   * the image that contains the PSF you want to end up with (`‑t …`)  
   * a segmentation‑mask that flags where objects are (`‑m …`) – the mask is used to avoid noisy background when estimating the noise level  
   * optional arguments like a scaling factor for the noise, whether to use the GPU, whether to apply a window in Fourier space, etc.  

   The script prints a short summary of the options when you ask for verbosity.

3. **Step‑by‑step what the code does**

   | Step | What the code does (in everyday language) |
   |------|-------------------------------------------|
   | **0** | Starts a timer so it can tell you how long it took. |
   | **1** | Reads the segmentation mask (`SEGMENTATION.fits`) and the two PSF files (the “old” PSF and the “new” PSF you want). |
   | **2** | Looks at the header of the science image to figure out which instrument made it (HST, JWST‑NIRCam, JWST‑MIRI). |
   | **3** | Loads the science image (the actual picture) into a NumPy array. |
   | **4** | Because the algorithm works best on images with an odd number of pixels, it pads the image by one pixel if needed. |
   | **5** | **Builds two big “filter” images** – the old PSF and the target PSF – by putting each PSF in the centre of a zero‑filled array that has the same size as the science image. |
   | **6** | **Estimates the noise level** in the image (the variance of the background).  This is used later to regularise the de‑convolution so that it does not amplify noise. |
   | **7** | **Performs the PSF‑matching**: it works in Fourier space (the “frequency” version of the image).  Using a Wiener‑like filter it figures out how to multiply the Fourier transform of the image so that, after transforming back, the image now has the *target* blur.  The heavy‑lifting can be done on the GPU if you asked for it. |
   | **8** | (Optional) **Applies a smooth “window”** in Fourier space (a cosine‑bell function).  This reduces ringing artefacts that sometimes appear when you do large‑scale Fourier transforms. |
   | **9** | **If the data come from JWST**, the script repeats steps 5‑8 for the associated error (ERR) extension, but this time it treats the error image as a variance map and finally returns an error image that matches the newly‑blurred science image. |
   | **10**| **Because very large images can be memory‑intensive**, the script cuts the image into four roughly equal quadrants, adds a buffer around each piece (large enough to contain the PSF), runs the whole PSF‑matching process on each piece, and then stitches the pieces back together. |
   | **11**| **Writes the result** to a new FITS file.  For HST images it just replaces the primary data array; for JWST it writes out the usual extensions (`SCI`, `ERR`, `EXP`, `WHT`, and optionally `NIM`). |
   | **12**| Prints how long the whole operation took. |

4. **Key ideas behind the math (but still in layman terms)**  

   * **Fourier transforms** – turning the image into a “frequency” picture lets you manipulate the blur much more easily, because a blur (the PSF) becomes a simple multiplication in Fourier space.  
   * **Wiener‑like filter** – the script builds a filter that balances two things: (i) the known shape of the old PSF, (ii) the desired shape of the new PSF, and (iii) the amount of random noise in the image.  This prevents the filter from wildly amplifying noise while still giving the right amount of extra smoothing.  
   * **Noise regularisation** – before the filter is built the script measures the typical background variance (how much the pixel values fluctuate in empty sky) and scales it by a user‑supplied factor (`‑s`).  This “regularisation” makes the algorithm stable.  
   * **Windowing** – the optional cosine‑bell window gently tapers the high‑frequency edges of the Fourier transform, which reduces ringing artefacts that sometimes appear after an inverse transform.

5. **Why it matters**  

   * After you run the script, the science image and its associated error map now share a *common* PSF – a blur that is the same everywhere and matches the target you supplied.  
   * This is essential when you want to combine data from different telescopes, or when you want to compare measurements (e.g., surface‑brightness profiles) across images that originally had different resolutions.  

6. **Bottom line**  

   The program is a fairly sophisticated “image‑blurring” tool for astronomical data.  
   It reads a picture, figures out how blurry it currently is, computes how to make it look like a different, user‑specified blur, does the convolution (optionally on a GPU and in manageable chunks), updates the error map accordingly, and writes out a new FITS file that can be used in downstream scientific analysis.
