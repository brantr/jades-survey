.. _flagging_tool:

=============================================
``jades-flag-bad-pixels-hdu.py``: Photometric Catalog Flagging Utility (Python)
=============================================

.. image:: https://img.shields.io/badge/python-3.8%2B-blue.svg
   :target: https://www.python.org/
.. image:: https://img.shields.io/badge/astropy-5.x-green.svg
   :target: https://www.astropy.org/
.. image:: https://img.shields.io/badge/license-MIT-brightgreen.svg
   :target: https://opensource.org/licenses/MIT

-----------------------------------------------------------------
A command‑line script that reads a list of FITS images, inspects the
pixels that belong to each object in an input photometric catalogue and
stores per‑filter “bad‑pixel” flags in a new FITS binary table.

The document below explains the source line‑by‑line, the expected
behaviour, and how to invoke the program from a terminal or from a
Python environment.  It is written in *reStructuredText* (RST) so it can
be rendered directly on a `Read the Docs <https://readthedocs.org/>`_
site.

-----------------------------------------------------------------

.. contents::
   :depth: 2
   :local:

-----------------------------------------------------------------
1.  Overview
-----------------------------------------------------------------

The script is intended for the following workflow:

* A **photometric catalogue** (ASCII table) contains a list of detected
  sources together with their bounding boxes (`bbox_xmin`,
  `bbox_xmax`, `bbox_ymin`, `bbox_ymax`).

* One or more **science images** (and optionally associated weight or
  error images) are supplied in a plain‑text list file.

* For each image the script determines the filter name (either from the
  file name or from the FITS header) and creates a new column in a
  `astropy.table.Table` called ``<FILTER>_FLAG``.

* Inside each object's bounding box the script counts:

  - NaN pixels in the science image,
  - Pixels equal to zero in the weight image (if supplied),
  - (Optionally) extreme or NaN values in the error image.

* The maximum number of “bad” pixels found for each object across all
  images of the same filter is written into the corresponding flag
  column.

* Finally the table is saved as a FITS binary table (extension name
  ``FLAG``) together with an empty primary HDU.

-----------------------------------------------------------------
2.  External Dependencies
-----------------------------------------------------------------

The script imports the following third‑party packages:

* **numpy** – fast array operations
* **astropy** – FITS I/O (`astropy.io.fits`), ASCII table reading
  (`astropy.io.ascii`) and the high‑level table class (`astropy.table.Table`)
* **tqdm** – progress‑bar for the inner loop (optional but highly
  recommended)
* **argparse** – standard library for command‑line parsing
* **time**, **os**, **sys** – standard library utilities

All of them are pure‑Python and available on *PyPI*:

.. code-block:: bash

   pip install numpy astropy tqdm

-----------------------------------------------------------------
3.  Command‑Line Interface
-----------------------------------------------------------------

The helper function :func:`create_parser` builds an ``ArgumentParser`` that
exposes the following options.  The table below mirrors the definition
in the source code.

+----------------------+----------------------+-----------------------------------+
| Argument             | Default              | Meaning                           |
+======================+======================+===================================+
| ``-i`` / ``--input_list`` | ``input_list.txt`` | Text file with one FITS image path per line |
+----------------------+----------------------+-----------------------------------+
| ``-w`` / ``--wht_list``   | ``None``          | Optional list of weight images   |
+----------------------+----------------------+-----------------------------------+
| ``--sci-ext``            | ``SCI``           | Extension name for the science   |
|                          |                  | image when the weight and science |
|                          |                  | are stored in the same file       |
+----------------------+----------------------+-----------------------------------+
| ``--wht-ext``            | ``PRIMARY``       | Extension name for the weight image |
+----------------------+----------------------+-----------------------------------+
| ``--err-ext``            | ``ERR``           | Extension name for the error image |
+----------------------+----------------------+-----------------------------------+
| ``-e`` / ``--err_list``   | ``None``          | Optional list of error images    |
+----------------------+----------------------+-----------------------------------+
| ``-c`` / ``--cat``        | ``phot_cat_uncorrected.fits`` | Input photometric catalogue (ASCII) |
+----------------------+----------------------+-----------------------------------+
| ``--nlim``               | ``None``          | Truncate the catalogue to the first *n* rows |
+----------------------+----------------------+-----------------------------------+
| ``-o`` / ``--output``     | ``phot_cat_corrected.fits`` | Output FITS file containing the flag table |
+----------------------+----------------------+-----------------------------------+
| ``-v`` / ``--verbose``   | ``False``         | Print diagnostic information     |
+----------------------+----------------------+-----------------------------------+

Typical invocation:

.. code-block:: bash

   python flagging_tool.py \
       -i images.txt \
       -w wht.txt \
       -e err.txt \
       -c phot_cat_uncorrected.fits \
       -o phot_cat_corrected.fits \
       -v

-----------------------------------------------------------------
4.  Core Functions
-----------------------------------------------------------------

The script contains only two public functions:

* :func:`create_parser`
* :func:`main`

Both are explained in the following subsections.

-----------------------------------------------------------------
4.1  ``create_parser()``
-----------------------------------------------------------------

```python
def create_parser():
    parser = argparse.ArgumentParser(
        description="Detection flags and options from user.")
    # ... add_argument calls ...
    return parser
```

* Builds an ``ArgumentParser`` with a description.
* Each ``add_argument`` mirrors one row of the table above.
* ``action='store_true'`` for ``--verbose`` creates a boolean flag.
* The function returns the ready‑to‑use parser; the actual parsing is
  performed later in :func:`main`.

-----------------------------------------------------------------
4.2  ``main()``
-----------------------------------------------------------------

The *heart* of the script.  Its logical flow is visualised below
(see also the inline comments in the source).

````text
┌─────────────────────┐
│ Parse CLI arguments │
└─────────┬───────────┘
          │
          ▼
   ┌───────────────┐
   │ Open list(s)  │
   └───────┬───────┘
           │
   ┌───────▼───────┐
   │ Read catalogue│
   └───────┬───────┘
           │
   ┌───────▼───────┐
   │ Initialise  │
   │ output Table │
   └───────┬───────┘
           │
   ┌───────▼───────────────────────┐
   │ Loop over every FITS image    │
   │   – read data, header          │
   │   – infer FILTER name          │
   │   – (optional) read weight/err │
   │   – create column <FILTER>_FLAG│
   │   – Loop over every object    │
   │       * clip bbox to image    │
   │       * count NaNs in SCI     │
   │       * count zeros in WHT    │
   │       * (optional) bad err    │
   │       * store max bad count   │
   └───────────────┬───────────────┘
                   ▼
   ┌─────────────────────────────────────┐
   │ Write output FITS (Primary + BinTable)│
   └─────────────────────────────────────┘
````

Below we walk through the most important blocks.

-----------------------------------------------------------------
4.2.1  Argument handling & timing
-----------------------------------------------------------------

```python
parser = create_parser()
args   = parser.parse_args()
t_start = time.time()
```

* ``args`` holds all CLI values.
* ``t_start`` is used later to report total execution time.

If ``--verbose`` is set, the script prints the supplied file names and
other options – useful for debugging.

-----------------------------------------------------------------
4.2.2  Loading the input lists
-----------------------------------------------------------------

```python
fp_input_list = open(args.input_list)
fl_input_list = [l.strip('\n') for l in fp_input_list.readlines()]
fp_input_list.close()
```

* ``fl_input_list`` becomes a Python ``list`` of image file paths.
* Similar code loads the optional weight and error lists, and sets
  ``flag_wht`` / ``flag_err`` booleans that gate later processing.

-----------------------------------------------------------------
4.2.3  Reading the catalogue
-----------------------------------------------------------------

```python
cat = ascii.read(args.cat)
if args.nlim is not None:
    cat = cat[:args.nlim]
nobj = len(cat)
```

* ``ascii.read`` automatically detects the delimiter (CSV, space‑delimited,
  etc.) and returns an ``astropy.table.Table``.
* The catalogue must contain at least the columns:
  ``id``, ``ra``, ``dec``, ``bbox_xmin``, ``bbox_xmax``,
  ``bbox_ymin``, ``bbox_ymax``.
* ``nobj`` is the number of sources to be examined.

-----------------------------------------------------------------
4.2.4  Preparing the output table
-----------------------------------------------------------------

```python
t = Table()
t['ID']  = cat['id'].copy()
t['RA']  = cat['ra'].copy()
t['DEC'] = cat['dec'].copy()
```

Only the identifier and sky coordinates are copied; flag columns are
added later inside the image‑loop.

The bounding boxes are *expanded* by a 10‑pixel margin (both sides) to
provide a safety buffer:

```python
bbox_xmin = cat['bbox_xmin'].copy() - 10
bbox_xmax = cat['bbox_xmax'].copy() + 11   # +10 + 1
bbox_ymin = cat['bbox_ymin'].copy() - 10
bbox_ymax = cat['bbox_ymax'].copy() + 11
```

-----------------------------------------------------------------
4.2.5  Main image loop
-----------------------------------------------------------------

```python
for i in range(n_images):
    fname_image = fl_input_list[i].strip("\n")
    print(f'Image: {fname_image}')
```

* ``n_images`` = length of the input list.
* ``fname_image`` is the path of the current science FITS file.

**Weight‑image handling**

If a weight list is supplied, ``fname_wht`` is taken from the same index.
If the weight file name is identical to the science file name the
script expects a *separate* extension (``WHT``) inside the same FITS
container; otherwise it reads the extension given by ``--wht-ext``.

**Header inspection & filter name extraction**

```python
header_file = fits.getheader(fname_image, 'PRIMARY')
if '560' in fname_image:
    FILTER = 'F560W'
elif '770' in fname_image:
    FILTER = 'F770W'
...
elif 'TELESCOP' in header_file:
    # HST case, JWST case, etc.
```

The logic first looks for known numeric substrings in the filename
(e.g., “560” → MIRI filter ``F560W``).  If none match, it falls back to
the ``TELESCOP`` keyword:

* For HST it prefers the ``FILTER`` keyword, otherwise ``FILTER1`` or
  ``FILTER2``.
* For other telescopes it strips blanks from the header ``FILTER``
  value and, if a ``MODULE`` keyword containing ``A`` or ``B`` exists,
  appends its first character (e.g., ``F770WA``).

The final string is used to create a column name:

```python
flag = FILTER + '_FLAG'      # e.g. "F770W_FLAG"
t[flag] = np.zeros(nobj)     # initialise with zeros
```

**Weight / error image loading (optional)**

```python
if flag_wht:
    data_wht = fits.getdata(fname_wht, args.wht_ext).astype(np.float32)
if flag_err:
    data_err = fits.getdata(fname_err, args.err_ext).astype(np.float32)
```

All image data are cast to ``float32`` to save memory.

-----------------------------------------------------------------
4.2.6  Per‑object loop (inner loop)
-----------------------------------------------------------------

The inner loop iterates over every catalogue source:

```python
for j in tqdm(range(nobj)):
    # Clip the bounding box to the image dimensions
    xmin = np.max([bbox_xmin[j], 0])
    xmax = np.min([bbox_xmax[j], shape[1]])
    ymin = np.max([bbox_ymin[j], 0])
    ymax = np.min([bbox_ymax[j], shape[0]])
```

* ``shape`` = ``data_image.shape`` (``(ny, nx)``).
* The clipping guarantees that array slicing never exceeds the image
  borders.

**Counting bad pixels**

*Science image* – NaN values:

```python
di   = data_image[ymin:ymax, xmin:xmax].copy()
idx  = np.where(np.isnan(di.ravel()))[0]
nbad = len(idx)
```

*Weight image* – zeros (or values below ``1e‑10``) indicate a pixel
that was masked during the reduction:

```python
if flag_wht:
    di   = data_wht[ymin:ymax, xmin:xmax].copy()
    idx  = np.where(np.abs(di) < 1.0e-10)[0]
    nbad = len(idx)
```

*Error image* – currently disabled (`if False:`).  The original author
left the code as a placeholder for future use.

**Storing the result**

The script keeps a temporary array ``flag_tmp`` that holds the *maximum*
number of bad pixels encountered for each object across all images of
the same filter:

```python
flag_tmp[j] = np.max([flag_tmp[j], nbad])
```

After the inner loop finishes, ``flag_tmp`` is copied into the table
column:

```python
t[flag] = flag_tmp.copy()
```

Thus each flag column contains, for a given filter, the worst‑case
bad‑pixel count among all images of that filter.

-----------------------------------------------------------------
4.2.7  Writing the output FITS file
-----------------------------------------------------------------

```python
hdu_pri = fits.PrimaryHDU()
hdu_tab = fits.BinTableHDU(data=t, name='FLAG')
hdul = fits.HDUList([hdu_pri, hdu_tab])
hdul.writeto(args.output, overwrite=True)
```

* An empty primary HDU is required by the FITS standard.
* ``t`` (the `astropy.table.Table` built above) becomes a binary table
  HDU named ``FLAG``.
* ``overwrite=True`` lets the script replace an existing file.

-----------------------------------------------------------------
4.2.8  Timing and final message
-----------------------------------------------------------------

```python
t_end = time.time()
print(f"Time to set flags = {t_end-t_start} seconds.")
```

-----------------------------------------------------------------
5.  Execution Flow Summary (Pseudo‑code)
-----------------------------------------------------------------

```python
def main():
    args = parse_cli()
    start = now()
    sci_list = read_list(args.input_list)
    wht_list = read_list(args.wht_list) if args.wht_list else []
    err_list = read_list(args.err_list) if args.err_list else []

    cat = read_ascii(args.cat)
    if args.nlim: cat = cat[:args.nlim]

    out = Table()
    out['ID'], out['RA'], out['DEC'] = cat['id'], cat['ra'], cat['dec']

    # expand bounding boxes
    xbmin, xbmax = cat['bbox_xmin']-10, cat['bbox_xmax']+11
    ybmin, ybmax = cat['bbox_ymin']-10, cat['bbox_ymax']+11

    for img_path, wht_path, err_path in zip(sci_list, wht_list, err_list):
        data_sci = fits.getdata(img_path, args.sci_ext).astype(np.float32)
        header   = fits.getheader(img_path, 'PRIMARY')
        FILTER   = infer_filter(img_path, header)
        flag_col = f'{FILTER}_FLAG'
        out[flag_col] = np.zeros(len(cat))

        if wht_path:
            data_wht = fits.getdata(wht_path, args.wht_ext).astype(np.float32)
        if err_path:
            data_err = fits.getdata(err_path, args.err_ext).astype(np.float32)

        for i, (xmin, xmax, ymin, ymax) in enumerate(zip(xbmin, xbmax, ybmin, ybmax)):
            # clip to image size
            xmin, xmax = clip(xmin, xmax, data_sci.shape[1])
            ymin, ymax = clip(ymin, ymax, data_sci.shape[0])

            bad = count_nan(data_sci[ymin:ymax, xmin:xmax])
            if wht_path:   bad = max(bad, count_zeros(data_wht[ymin:ymax, xmin:xmax]))
            # if err_path:   bad = max(bad, count_error_pixels(...))

            out[flag_col][i] = max(out[flag_col][i], bad)

    write_fits(args.output, out)
    print(f'Done in {now()-start:.1f}s')
```

-----------------------------------------------------------------
6.  Usage Examples
-----------------------------------------------------------------

**Basic run (science images only)**

.. code-block:: bash

   python flagging_tool.py -i sci_images.txt -c mycat.txt -o flagged.fits -v

**Including weight images**

.. code-block:: bash

   python flagging_tool.py \
       -i sci_images.txt \
       -w wht_images.txt \
       --wht-ext WHT \
       -c mycat.txt \
       -o flagged.fits

**Limiting the catalogue to the first 500 objects**

.. code-block:: bash

   python flagging_tool.py -i sci.txt -c mycat.txt --nlim 500 -o out.fits

-----------------------------------------------------------------
7.  Remarks & Potential Improvements
-----------------------------------------------------------------

* **Error‑image handling** – the code block for ``flag_err`` is wrapped
  in ``if False``; activating it would require a robust definition of
  “bad” error values (e.g., ``np.isinf`` or thresholds).

* **Memory usage** – for very large images the ``copy()`` of the
  sub‑array inside the inner loop can be avoided by working directly on
  the slice view (e.g., ``sub = data_image[ymin:ymax, xmin:xmax]``).

* **Parallelism** – the per‑object loop is a perfect candidate for
  ``multiprocessing`` or ``numba`` JIT compilation, especially when
  dealing with thousands of sources.

* **Filter inference** – the current heuristic works for JWST MIRI
  filters and HST but may fail for other instruments.  A more generic
  approach would read the ``FILTER`` keyword first and only fall back
  to filename parsing if the keyword is missing.

* **Unit tests** – adding a small suite that creates synthetic FITS
  files and a mock catalogue would guarantee that future refactoring
  does not break the flag‑counting logic.

-----------------------------------------------------------------
8.  License
-----------------------------------------------------------------

The script is released under the **MIT License**.  See the accompanying
``LICENSE`` file for the full text.

-----------------------------------------------------------------
9.  Author & Contact
-----------------------------------------------------------------

*Original author*: *[Your Name]*  
*Maintainer*: *[Your Email]*  

Feel free to open an issue on the project's GitHub repository for
questions, bug reports, or feature requests.

