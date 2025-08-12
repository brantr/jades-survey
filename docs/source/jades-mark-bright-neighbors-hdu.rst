.. _bright-neighbor-flagger:

Bright‑Neighbor Flagging Utility
================================

This document explains the purpose, functionality and inner workings of the
``bright_neighbor_flagger.py`` script.  The script is a small command‑line
utility used in astronomical image processing pipelines to flag objects that
have a **bright neighbour** in a deblended catalog.  The resulting flag is
written to a FITS binary table that can be ingested by downstream analysis
steps (e.g. source classification, photometric redshift estimation, …).

The description below follows the structure of the source file and provides
detailed commentary on each block of code, the external libraries it relies
on and the algorithmic decisions that were made.

.. contents::
   :local:
   :depth: 2


--------------------------------------------------------------------
1.  High‑level purpose
--------------------------------------------------------------------

* **Input**  
  * A *catalog* (ASCII table) produced by a source‑detection program (e.g.
    `SExtractor` or `sep`).  The catalog must contain at least the columns
    ``id``, ``ra``, ``dec`` and a photometric measurement named
    ``kron_flux`` (the script uses this as a proxy for source brightness).
  * A *segmentation map* (FITS image) where each pixel value is the integer
    identifier of the object that occupies that pixel.  The segmentation map
    is assumed to be the *deblended* version – i.e. each source has been
    separated from its neighbours.

* **Output**  
  * A FITS binary table with three columns: ``ID``, ``RA``, ``DEC`` and a
    newly added column ``FLAG_BN``.  The flag encodes the presence of a bright
    neighbour:

    ===============  ==============================
    Flag value       Meaning
    ===============  ==============================
    ``0``            No bright neighbour found
    ``1``            At least one neighbour with flux > 2 × the target
    ``2``            At least one neighbour with flux > 10 × the target
    ===============  ==============================

The script is intentionally lightweight – it does **not** re‑run any
photometry; it only inspects the existing catalog and segmentation map to
derive the flag.

--------------------------------------------------------------------
2.  Dependencies
--------------------------------------------------------------------

The script imports the following third‑party packages:

* :mod:`tqdm` – progress bar for the main loop.
* :mod:`argparse` – command‑line argument parsing.
* :mod:`numpy` – numerical array handling.
* :mod:`sep` – (Source Extraction and Photometry) – imported but not used in the
  current version; retained for compatibility with other pipeline scripts.
* :mod:`astropy.io.fits` – reading/writing FITS files.
* :mod:`astropy.io.ascii` – reading the input ASCII catalog.
* :mod:`astropy.table.Table` – convenient container for the output catalog.
* :mod:`time` – simple wall‑clock timing.
* :mod:`scipy.stats.mode` – imported but not used (historical artifact).

All of these packages are part of the standard scientific Python stack and are
available from *pip* or *conda*.

--------------------------------------------------------------------
3.  Command‑line interface
--------------------------------------------------------------------

The function :func:`create_parser` builds an ``ArgumentParser`` instance:

.. code-block:: python

    def create_parser():
        parser = argparse.ArgumentParser(
            description="Set parent/daughter flags and bright neighbor flags.")
        parser.add_argument('-c', '--cat',
                            default='det_v10_final.cat',
                            metavar='cat',
                            type=str,
                            help='Input catalog.')
        parser.add_argument('--output',
                            default='output.fits',
                            metavar='output',
                            type=str,
                            help='Output catalog.')
        parser.add_argument('-s', '--segmap',
                            default='segmap.fits',
                            metavar='segmap',
                            type=str,
                            help='Deblended catalog segmap.')
        parser.add_argument('-v', '--verbose',
                            dest='verbose',
                            action='store_true',
                            help='Print helpful information to the screen? '
                                 '(default: False)',
                            default=False)
        return parser

* ``-c / --cat`` – path to the ASCII source catalog (default:
  ``det_v10_final.cat``).
* ``--output`` – name of the FITS file that will contain the flag table
  (default: ``output.fits``).
* ``-s / --segmap`` – path to the deblended segmentation map (default:
  ``segmap.fits``).
* ``-v / --verbose`` – if supplied, the script prints the names of the input
  files and a few status messages.

Running the script with ``-h`` shows the help text generated by ``argparse``.

--------------------------------------------------------------------
4.  Core algorithm (``main`` function)
--------------------------------------------------------------------

The bulk of the work happens inside :func:`main`.  The function can be split
into logical stages that are described below.

--------------------------------------------------------------------
4.1  Initialisation & data loading
--------------------------------------------------------------------

.. code-block:: python

    parser = create_parser()
    args = parser.parse_args()

    t_start = time.time()

    if args.verbose:
        print(f"Catalog {args.cat}")
        print(f"Segmap  {args.segmap}")

    cat = ascii.read(args.cat)                # Astropy Table
    segmap = fits.getdata(args.segmap)        # 2‑D NumPy array
    segmap = segmap.astype(np.int32)          # ensure integer IDs

* ``cat`` now holds the source catalog; column access is performed via the
  table’s field names (e.g. ``cat['id']``).
* ``segmap`` is a 2‑D image where each pixel contains the *object ID* of the
  source that occupies that pixel.  Converting to ``int32`` guarantees that
  later NumPy comparisons are fast and memory‑efficient.

--------------------------------------------------------------------
4.2  Defining a region of interest around each object
--------------------------------------------------------------------

Every object in the catalog already has a bounding box (``bbox_xmin``,
``bbox_xmax``, ``bbox_ymin``, ``bbox_ymax``) supplied by the detection
software.  The script expands this box by a 10‑pixel margin on each side so
that a modest surrounding area is examined.  This margin is hard‑coded and
provides a safety buffer when the segmentation map contains small
inaccuracies at the edges of objects.

.. code-block:: python

    bbox_xmin = cat['bbox_xmin'].copy() - 10
    bbox_xmax = cat['bbox_xmax'].copy() + 11
    bbox_ymin = cat['bbox_ymin'].copy() - 10
    bbox_ymax = cat['bbox_ymax'].copy() + 11

The ``+11`` (instead of ``+10``) compensates for Python’s half‑open slice
semantics (``[start:stop]`` does **not** include ``stop``).

--------------------------------------------------------------------
4.3  Preparing the output table
--------------------------------------------------------------------

Only a subset of the original catalog columns is retained in the final FITS
file.  The script creates a fresh :class:`astropy.table.Table` and copies the
relevant fields:

.. code-block:: python

    t = Table()
    t['ID']  = cat['id'].copy()
    t['RA']  = cat['ra'].copy()
    t['DEC'] = cat['dec'].copy()

The ``FLAG_BN`` column will be added later once the neighbour check is
completed.

--------------------------------------------------------------------
4.4  Loop over every object – neighbour detection
--------------------------------------------------------------------

The most important part of the script is the ``for`` loop that iterates over
all objects (``nobj = len(cat)``).  ``tqdm`` supplies a progress bar that
updates in the terminal, which is very handy for large catalogs.

The loop proceeds as follows for each object ``i``:

1. **Extract the expanded bounding box** and clip it to the image limits.
   This prevents ``IndexError`` when an object lies close to the edge of the
   frame.

   .. code-block:: python

        xmin = max(bbox_xmin[i], 0)
        ymin = max(bbox_ymin[i], 0)
        xmax = min(bbox_xmax[i], segmap.shape[1])
        ymax = min(bbox_ymax[i], segmap.shape[0])

2. **Slice the segmentation map** to obtain ``ds_cut``, a small cut‑out that
   contains the current object *and* any other IDs that overlap the region.

   .. code-block:: python

        ds_cut = segmap[ymin:ymax, xmin:xmax].copy()

3. **Identify pixels that do *not* belong to the current object**.  ``np.where``
   returns the indices where the condition is true; we then take the values
   at those indices to obtain the set of *other* IDs present in the cut‑out.

   .. code-block:: python

        idx_cut = np.where(ds_cut != id_obj)          # id_obj = cat['id'][i]
        du = np.unique(ds_cut[idx_cut])               # unique neighbour IDs
        du = du[du > 0]                                # discard background (0)

4. **Initialize the flag** for this object to ``0`` (no bright neighbour).

   .. code-block:: python

        flag_bn[i] = 0

5. **Compare fluxes** of the neighbours with the flux of the current object.
   The script uses the column ``kron_flux`` (the default in many SExtractor
   runs) and applies two thresholds:

   * ``> 2 ×`` → flag becomes ``1``.
   * ``> 10 ×`` → flag becomes ``2`` (overwrites a possible previous ``1``).

   The inner ``for`` loop iterates over each neighbour ID ``idc`` and
   retrieves its catalog row via ``np.where(cat['id'] == idc)``.  Because
   ``cat['id']`` is guaranteed to be unique, the result is a one‑element
   array, and indexing ``cat[FLUX_CHECK][idx]`` yields the neighbour’s flux.

   .. code-block:: python

        FLUX_CHECK = 'kron_flux'
        if len(du) > 0:
            for j, idc in enumerate(du):
                idx = np.where(cat['id'] == idc)
                if cat[FLUX_CHECK][idx] > 2 * cat[FLUX_CHECK][i]:
                    flag_bn[i] = 1
                if cat[FLUX_CHECK][idx] > 10 * cat[FLUX_CHECK][i]:
                    flag_bn[i] = 2

   Note that the flag can be upgraded from ``1`` to ``2`` if a neighbour
   satisfies the stricter criterion.

The result of the loop is a NumPy array ``flag_bn`` of length ``nobj`` that
holds the neighbour flag for every source.

--------------------------------------------------------------------
4.5  Assemble final table and write to FITS
--------------------------------------------------------------------

After the loop finishes, the flag column is appended to the output table
``t`` and the table is written as a binary FITS extension named ``FLAG``.

.. code-block:: python

    t['FLAG_BN'] = flag_bn.copy()

    hdu_pri = fits.PrimaryHDU()
    hdu_tab = fits.BinTableHDU(data=t, name='FLAG')
    hdul = fits.HDUList([hdu_pri, hdu_tab])
    hdul.writeto(args.output, overwrite=True)

The primary HDU is empty (just a placeholder) because the useful data lives
in the binary table extension.

--------------------------------------------------------------------
4.6  Timing and termination
--------------------------------------------------------------------

The script reports how long the flag‑setting stage took:

.. code-block:: python

    t_end = time.time()
    print(f"Time to set flags = {t_end - t_start} seconds.")
```

--------------------------------------------------------------------
5.  Usage example
--------------------------------------------------------------------

Assuming the script is saved as ``bright_neighbor_flagger.py`` and that the
current working directory contains the files ``det_v10_final.cat`` and
``segmap.fits``:

.. code-block:: bash

    $ python bright_neighbor_flagger.py \
          -c det_v10_final.cat \
          -s segmap.fits \
          --output flagged_sources.fits \
          -v

The ``-v`` switch will echo the input file names.  After a few seconds (or
minutes for very large catalogs) a file ``flagged_sources.fits`` will be
created.  Its content can be inspected with ``fitslist`` or with Astropy:

.. code-block:: python

    >>> from astropy.io import fits
    >>> hdul = fits.open('flagged_sources.fits')
    >>> hdul['FLAG'].data[:5]   # first five rows
    (ID, RA, DEC, FLAG_BN) ...

--------------------------------------------------------------------
6.  Extending or modifying the script
--------------------------------------------------------------------

* **Different photometric proxy** – If the catalog stores a different flux
  column (e.g. ``aperture_flux``), change the ``FLUX_CHECK`` variable near the
  top of the neighbour loop.

* **Additional flag criteria** – Insert extra ``if`` statements inside the
  neighbour loop; remember to respect the hierarchy (e.g. a flag of ``3`` can
  be introduced for a new condition).

* **Parallelisation** – The per‑object loop is embarrassingly parallel.  One
  could replace the ``for i in tqdm(range(nobj))`` construct with a
  ``multiprocessing.Pool`` map, taking care to pass a *read‑only* copy of the
  segmentation map to each worker.

* **Removing unused imports** – ``sep`` and ``scipy.stats.mode`` are not used
  in the current version; they can be safely removed to speed up import time.

--------------------------------------------------------------------
7.  Summary
--------------------------------------------------------------------

The ``bright_neighbor_flagger`` script provides a lightweight, pure‑Python
method for annotating a source catalog with a flag that signals the presence
of significantly brighter neighbours within a modest surrounding region.  It
leverages Astropy for I/O, NumPy for fast array manipulation, and tqdm for a
user‑friendly progress bar.  The resulting FITS binary table can be fed into
any downstream astronomical analysis pipeline that needs to treat bright‑
neighbour contamination specially.

--------------------------------------------------------------------
8.  Full source listing (for reference)
--------------------------------------------------------------------

.. literalinclude:: bright_neighbor_flagger.py
   :language: python
   :linenos:
