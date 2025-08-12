.. _bithash_tool:


=========================================
``jades-create-program-bit-hash.py``: Bithash Image Generator
=========================================

A small command‑line utility that builds a *bithash* FITS image from a list of
individual exposure files.  
Each input image contributes a single bit (or a set of bits) that encodes the
observing program(s) that produced the data.  The resulting ``BITHASH`` image
can be used as a quick look‑up table to see which programs cover any given
pixel in a mosaic.

The script was written for the **JADES** (JWST Advanced Deep Extragalactic Survey)
programs but works for any collection of FITS images that contain a
``WHT`` (weight) extension indicating valid data regions.


Why a bithash?
--------------

* **Compact provenance** – a 32‑bit integer can store up to 32 distinct
  programs.  By setting the *n*‑th bit for a program we obtain a single image
  that tells us *exactly* which programs contributed to each pixel.
* **Fast masking** – downstream tools can simply test a pixel with a bitwise
  ``AND`` operation to know whether a particular program contributed.
* **Human readable** – the script also writes a ``DECODER`` table that maps
  file names to the bit they occupy, making the hash self‑documenting.



Table of bit assignments
------------------------

The comment block at the top of the script defines the mapping from bit
positions (0‑27) to JWST program identifiers.  For example:

+--------+------------------------------------------+
| Bit    | Program (example)                        |
+========+==========================================+
| 0      | 1180 Deep (Obs 7‑18)                     |
+--------+------------------------------------------+
| 1      | 1180 Medium/HST (Obs 25‑30 + 136)        |
+--------+------------------------------------------+
| …      | …                                        |
+--------+------------------------------------------+
| 27     | PRIMER (Obs 1837)                        |
+--------+------------------------------------------+

You may edit this list or add new bits – just keep the number of bits ≤ 31
(the script stores the result in a signed 32‑bit integer).



Installation
------------

The script depends only on a few widely‑used scientific Python packages:

* ``numpy``
* ``astropy``
* ``tqdm`` (optional – provides the progress bar)
* ``argparse`` (standard library)

Install them with ``pip`` if they are not already present:

.. code-block:: bash

    pip install numpy astropy tqdm


Running the tool
----------------

The program is invoked from the command line:

.. code-block:: bash

    python bithash.py -i input.txt -o output.fits [-v]

* ``-i / --input`` – text file that lists **one FITS filename per line**
  followed by a space and the **bit number** that should be set for that file.
  Example ``input.txt``:

  .. code-block:: text

      jw011800_deep.fits   0
      jw011800_medium.fits 1
      jw011800_medium23.fits 2
      jw011810_hst.fits   3
      ...

* ``-o / --output`` – name of the FITS file that will contain the bithash
  image and the decoder table.
* ``-v / --verbose`` – print additional progress information.



How the code works
------------------

Below is a line‑by‑line walk‑through of the most important sections.  The
original source is reproduced in a ``code-block`` for reference.

.. code-block:: python

    # ----------------------------------------------------------------------
    # Imports
    # ----------------------------------------------------------------------
    import numpy as np
    from astropy.io import fits
    from astropy.table import Table
    import argparse
    from tqdm import tqdm
    import time

The script only needs ``numpy`` for array handling, ``astropy`` for FITS I/O,
``argparse`` for command‑line parsing, ``tqdm`` for a nice progress bar, and
the standard ``time`` module to measure execution time.

----------------------------------------------------------------------
Parser creation
----------------------------------------------------------------------
The ``create_parser`` function builds the ``argparse.ArgumentParser`` object.
It defines three options (``-i``, ``-o`` and ``-v``) and returns the parser.

----------------------------------------------------------------------
Main workflow
----------------------------------------------------------------------
The ``main`` function orchestrates everything:

1. **Timing start** – ``time.time()`` is stored so the total runtime can be
   printed when ``--verbose`` is used.

2. **Parse arguments** – ``args = parser.parse_args()``.

3. **Read the input list** – the file given with ``-i`` is read line‑by‑line.
   Each line is split on whitespace; the first token is the filename,
   the second token is converted to an integer bit number.

   .. code-block:: python

        fp_input_list = open(args.input)
        fl_input_list = fp_input_list.readlines()
        fp_input_list.close()

        n_images = len(fl_input_list)

        fnames = []                     # list of filenames
        bits   = np.zeros(n_images, dtype=np.int32)   # corresponding bits

        for i in range(n_images):
            fnames.append(fl_input_list[i].split(' ')[0])
            bits[i] = int(fl_input_list[i].split(' ')[1])
            if args.verbose:
                print(f'fnames[{i}] {fnames[i]} bits[{i}] {bits[i]}')

4. **Iterate over *unique* bits** – many images may share the same bit (e.g.
   several exposures from the same program).  ``np.unique(bits)`` yields the
   set of distinct bit numbers.  The outer loop (wrapped by ``tqdm``) processes
   each bit **once**.

   .. code-block:: python

        unique_bits = np.unique(bits)
        print(f'unique_bits {unique_bits}')

        k = 0
        for bit in tqdm(unique_bits):

5. **Select images belonging to the current bit** – ``np.where(bits==bit)``.
   ``idx_image`` contains the indices of all files that should contribute to
   the current bit.

6. **Load the weight maps** – for each selected image the script reads the
   ``WHT`` extension (weight map) from the FITS file.  The weight map is a
   2‑D array where positive values indicate *valid* data.

   .. code-block:: python

            wht = fits.getdata(fnames[idx_image[j]], 'WHT')
            wht = wht.astype(np.float32)

7. **Create a binary mask for the current bit** – the first image for the
   bit creates ``current_bit_image`` (zeros of the same shape as the weight
   map).  Wherever ``wht > 0`` the mask is set to ``1``.

   .. code-block:: python

            if j == 0:
                current_bit_image = np.zeros_like(wht, dtype=np.int64)

            idx = np.where(wht > 0)
            current_bit_image[idx] = 1

8. **Initialize the output image** – on the very first iteration
   (``k == 0`` and ``j == 0``) the script also creates ``bit_image``, a zero
   array that will hold the final hash.  It also grabs the primary header and
   the ``SCI`` header from the first FITS file – these are copied to the
   output file unchanged.

   .. code-block:: python

            if (k == 0) & (j == 0):
                bit_image = np.zeros_like(wht, dtype=np.int32)
                header = fits.getheader(fnames[idx_image[j]], 'SCI')
                header_pri = fits.getheader(fnames[idx_image[j]])

9. **Shift the binary mask by the appropriate bit** – after all images that
   share the same bit have been processed, ``current_bit_image`` contains a
   mask of ``1`` where *any* of those images have coverage.  The mask is then
   left‑shifted by ``bit`` positions (``<< bit``) so that the ``bit``‑th bit
   of the integer is set.

   .. code-block:: python

            bit_image += (current_bit_image << bit)

   The addition works because each pixel of ``bit_image`` is an integer; the
   shifted mask adds the appropriate power‑of‑two to the pixel value.

10. **Book‑keeping** – ``k`` is incremented, the temporary ``current_bit_image``
    is deleted to free memory, and the outer loop proceeds to the next unique
    bit.

11. **Write the decoder table** – a small ``astropy.table.Table`` maps each
    filename to its assigned bit.  This table is stored as a binary table HDU
    called ``DECODER`` in the final FITS file.

    .. code-block:: python

        t = Table([fnames, bits], names=['image', 'bit'])

12. **Create the output FITS file** – three HDUs are assembled:

    * Primary HDU (only the header, no image data)
    * Image HDU named ``BITHASH`` containing the final integer array
    * Binary table HDU named ``DECODER``

    .. code-block:: python

        pri_hdu = fits.PrimaryHDU(header=header_pri)
        img_hdu = fits.ImageHDU(data=bit_image, header=header, name='BITHASH')
        tbl_hdu = fits.BinTableHDU(data=t, name='DECODER')
        hdu_out = fits.HDUList([pri_hdu, img_hdu, tbl_hdu])
        hdu_out.writeto(args.output, overwrite=True)

13. **Timing end** – if ``--verbose`` is set the total runtime is printed.

    .. code-block:: python

        time_global_end = time.time()
        if args.verbose:
            print(f"Time to execute program: {time_global_end - time_global_start}s.")

----------------------------------------------------------------------
Key implementation details
----------------------------------------------------------------------

* **Bitwise logic** – The core operation is ``current_bit_image << bit``.
  ``<<`` shifts the binary representation of the integer left by ``bit``
  positions, which is equivalent to multiplying by ``2**bit``.  Adding the
  shifted mask to ``bit_image`` therefore sets the appropriate bit(s) for
  each pixel.

* **Data type choices** – ``bit_image`` is stored as ``np.int32``.  This
  allows up to 31 usable bits (the sign bit is reserved).  If you need more
  programs you can change the dtype to ``np.int64`` and adjust the header
  accordingly.

* **Weight‑map based coverage** – The script assumes that a positive value in
  the ``WHT`` extension means *the pixel is covered*.  This is the standard
  convention for JWST calibrated products.

* **Memory handling** – Only one mask (``current_bit_image``) is kept in
  memory at a time; after each bit is processed the mask is deleted.  This
  makes the script feasible even for large mosaics (e.g. 10 000 × 10 000
  pixels).

----------------------------------------------------------------------
Example usage
-------------

Suppose you have three calibrated JWST exposures:

* ``jw011800_deep.fits``  – belongs to bit **0**
* ``jw011800_medium.fits`` – belongs to bit **1**
* ``jw011810_hst.fits`` – belongs to bit **3**

Create an ``input.txt`` file:

.. code-block:: text

    jw011800_deep.fits   0
    jw011800_medium.fits 1
    jw011810_hst.fits    3

Run the tool:

.. code-block:: bash

    python bithash.py -i input.txt -o jw_bithash.fits -v

The resulting ``jw_bithash.fits`` will contain:

* ``BITHASH`` image – each pixel value is a sum of powers of two, e.g.
  ``5`` (binary ``0101``) means the pixel is covered by programs with bits
  0 and 2.
* ``DECODER`` table – mapping of the three filenames to bits 0, 1 and 3.

You can now query the hash with ``astropy`` or ``numpy``:

.. code-block:: python

    from astropy.io import fits
    hdul = fits.open('jw_bithash.fits')
    bithash = hdul['BITHASH'].data

    # Pixels covered by program bit 1?
    mask_bit1 = (bithash & (1 << 1)) != 0

----------------------------------------------------------------------
Limitations & future extensions
-------------------------------

* **Maximum of 31 bits** – limited by the signed 32‑bit integer storage.
  Switching to ``np.uint64`` would raise the limit to 63 bits.
* **Assumes a ``WHT`` extension** – If your data use a different convention
  (e.g. a mask called ``MASK``) you need to edit the ``fits.getdata`` call.
* **No handling of overlapping bit assignments** – The script simply adds
  bits; if the same bit is assigned to two different files the result is
  still correct (the bit stays set).  However, duplicate entries in the
  decoder table are not filtered out.
* **No WCS propagation** – The WCS from the first image is copied verbatim.
  All input images should be on the same pixel grid; otherwise the hash
  will be meaningless.

----------------------------------------------------------------------
Reference
---------

* Astropy documentation – https://docs.astropy.org/
* JWST data reduction pipeline – https://jwst-pipeline.readthedocs.io/
* `tqdm – Fast, extensible progress bar`_.

.. _tqdm – Fast, extensible progress bar: https://tqdm.github.io/

----------------------------------------------------------------------
License
-------

The script is released under the BSD 3‑Clause license (the same license
as Astropy).  Feel free to modify it for your own projects.

----------------------------------------------------------------------
``bithash.py`` – full source
----------------------------

.. code-block:: python

    #!/usr/bin/env python
    # -*- coding: utf-8 -*-

    import numpy as np
    from astropy.io import fits
    from astropy.table import Table
    import argparse
    from tqdm import tqdm
    import time

    '''
    Bit 0: 1180 Deep: Obs 7-18 (jw011800_deep, jw011800_deep23)
    Bit 1 : 1180 Medium/HST: Obs 25-30 + 136 (jw011800_medium, jw011800_medium_redo)
    Bit 2 : 1180 Medium/MIRI: Obs 19-24, 219-223 (jw011800_medium23, jw011800_medium_obs022, jw011800_medium_obs219, jw011800_medium_obs220+222, jw011800_medium_obs223)
    Bit 3: 1181 Medium HST+MIRI Obs 1-7 (jw011810_hst, jw011810_miri)
    Bit 4: 1181 Medium JWST Obs 9-11+98 (jw011810_jwst)
    Bit 5: 1210
    Bit 6: 1286
    Bit 7: 1287
    Bit 8: 3215 (JOF)
    Bit 9: 4540 (JADES GRISM)
    Bit 10: 5997 (OASIS)
    Bit 11: 6434 (SAPPHIRES)
    Bit 12: 6541 (Egami DDT GOODS-S)
    ---- end JADES Programs
    Bit 13: 1176 (Windhorst GTO)
    Bit 14: 1283 (EU MIRI UDF)
    Bit 15: 1895 (FRESCO)
    Bit 16: 1963 (JEMS)
    Bit 17: 2079 (NGDEEP)
    Bit 18: 2198 (Barrufet GOODS-S)
    Bit 19: 2514 (PANORAMIC)
    Bit 20: 2516 (Hodge GOODS-S)
    Bit 21: 2674 (Arrabal Haro GOODS-N)
    Bit 22: 3577 (CONGRESS GOODS-N)
    Bit 23: 3990 (Morishita GOODS-S)
    Bit 24: 6511 (MIRI/EC)
    --- end GOODS-S+N programs
    Bit 25: 1345 (CEERS)
    Bit 26: 1727 (COSMOS Web)
    Bit 27: 1837 (PRIMER)
    '''

    # ----------------------------------------------------------------------
    # Argument parser
    # ----------------------------------------------------------------------
    def create_parser():
        """Create the command‑line argument parser."""
        parser = argparse.ArgumentParser(
            description="Create a bithash image encoding contributing programs."
        )
        parser.add_argument(
            "-i",
            "--input",
            dest="input",
            default="input.txt",
            metavar="input",
            type=str,
            help="List of input filenames (column 1) and the corresponding bit (column 2).",
        )
        parser.add_argument(
            "-o",
            "--output",
            dest="output",
            default="output.fits",
            metavar="output",
            type=str,
            help="Output bithash FITS image filename.",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            dest="verbose",
            action="store_true",
            help="Print helpful information to the screen? (default: False)",
            default=False,
        )
        return parser

    # ----------------------------------------------------------------------
    # Main routine
    # ----------------------------------------------------------------------
    def main():
        """Entry point of the program."""
        # Timing
        time_global_start = time.time()

        # Parse arguments
        parser = create_parser()
        args = parser.parse_args()

        # ------------------------------------------------------------------
        # Load input list
        # ------------------------------------------------------------------
        with open(args.input) as fp:
            fl_input_list = fp.readlines()

        n_images = len(fl_input_list)

        fnames = []
        bits = np.zeros(n_images, dtype=np.int32)

        for i in range(n_images):
            parts = fl_input_list[i].split()
            fnames.append(parts[0])
            bits[i] = int(parts[1])
            if args.verbose:
                print(f"fnames[{i}] {fnames[i]} bits[{i}] {bits[i]}")

        # ------------------------------------------------------------------
        # Build the bithash
        # ------------------------------------------------------------------
        unique_bits = np.unique(bits)
        print(f"unique_bits {unique_bits}")

        k = 0
        for bit in tqdm(unique_bits):
            idx_image = np.where(bits == bit)[0]

            for j, img_idx in enumerate(idx_image):
                if args.verbose:
                    print(f"Processing {fnames[img_idx]}")

                wht = fits.getdata(fnames[img_idx], "WHT").astype(np.float32)

                if j == 0:
                    current_bit_image = np.zeros_like(wht, dtype=np.int64)

                # set mask where weight > 0
                idx = np.where(wht > 0)
                current_bit_image[idx] = 1

                # Initialise the final image on the very first iteration
                if (k == 0) and (j == 0):
                    bit_image = np.zeros_like(wht, dtype=np.int32)
                    header = fits.getheader(fnames[img_idx], "SCI")
                    header_pri = fits.getheader(fnames[img_idx])

            # Add the shifted mask to the global bithash
            bit_image += (current_bit_image << bit)
            k += 1
            del current_bit_image

        # ------------------------------------------------------------------
        # Write decoder table
        # ------------------------------------------------------------------
        t = Table([fnames, bits], names=["image", "bit"])

        # ------------------------------------------------------------------
        # Write output FITS file
        # ------------------------------------------------------------------
        pri_hdu = fits.PrimaryHDU(header=header_pri)
        img_hdu = fits.ImageHDU(data=bit_image, header=header, name="BITHASH")
        tbl_hdu = fits.BinTableHDU(data=t, name="DECODER")
        hdu_out = fits.HDUList([pri_hdu, img_hdu, tbl_hdu])
        hdu_out.writeto(args.output, overwrite=True)

        # Timing end
        time_global_end = time.time()
        if args.verbose:
            print(f"Time to execute program: {time_global_end - time_global_start}s.")

    # ----------------------------------------------------------------------
    # Run when executed as a script
    # ----------------------------------------------------------------------
    if __name__ == "__main__":
        main()
