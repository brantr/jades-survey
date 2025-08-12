========================================
extract_layers – FITS science‑extension copier
========================================

.. image:: https://img.shields.io/badge/python-3.8%2B-blue.svg
   :target: https://www.python.org/
   :alt: Python version

.. image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

--------------------------------------------------------------------
**Purpose**

The module provides a tiny command‑line utility that copies the *SCI* (science)
extension of a FITS image to a new FITS file while preserving its header.  
It is useful when you only need the calibrated image data from a multi‑extension
FITS file (for example, a raw image that also contains variance, mask, or
instrumental extensions) and want to work with a single‑extension file.

--------------------------------------------------------------------
**Table of contents**

.. contents::
   :local:
   :depth: 2

--------------------------------------------------------------------
1.  Overview
-----------

The script consists of three logical parts:

* **Imports** – loads the required third‑party libraries.
* **`extract_layers` function** – reads the *SCI* extension from ``input_fits`` and
  writes it (data + header) to ``output_fits``.
* **Command‑line interface** – parses ``sys.argv`` to allow the user to specify
  input and output filenames and then calls ``extract_layers``.

--------------------------------------------------------------------
2.  Dependencies
----------------

* **NumPy** – imported as ``np`` but not used in the current version.  It is
  retained for possible future numeric processing.
* **Astropy** – the :mod:`astropy.io.fits` package provides robust FITS I/O.
* **sys** – used to read command‑line arguments.

The module can be installed with::

   pip install numpy astropy

--------------------------------------------------------------------
3.  Code walk‑through
---------------------

Below is the full source code with line‑by‑line commentary.

.. code-block:: python
   :linenos:

   import numpy as np
   import sys
   from astropy.io import fits

   # ------------------------------------------------------------------
   # extract the header from a fits image
   # ------------------------------------------------------------------
   def extract_layers(input_fits, output_fits):
       """
       Copy the *SCI* extension of ``input_fits`` to ``output_fits``.

       Parameters
       ----------
       input_fits : str
           Path to the source FITS file that contains a ``'SCI'`` HDU.
       output_fits : str
           Destination path for the new single‑extension FITS file.

       The function performs three steps:
       1. **Read the header** of the ``'SCI'`` HDU using :func:`fits.getheader`.
       2. **Read the data array** of the same HDU using :func:`fits.getdata`.
       3. **Write a new FITS file** that contains only this data and header
          with :func:`fits.writeto`.  ``overwrite=True`` allows the target file
          to be replaced if it already exists.
       """
       # SCIENCE
       header = fits.getheader(input_fits, 'SCI')
       data   = fits.getdata(input_fits, 'SCI')
       # The original code built the name from ``input_fits``; the caller now
       # passes the exact output name, so we simply reuse it.
       fname = output_fits
       fits.writeto(fname, data=data, header=header, overwrite=True)

   # ------------------------------------------------------------------
   # Command‑line entry point
   # ------------------------------------------------------------------
   def main():
       """
       Entry point for the ``python -m extract_layers`` style invocation.

       Usage
       -----
       ``python extract_layers.py [input_fits] [output_fits]``

       * If *input_fits* is omitted, the script defaults to ``image.fits``.
       * If *output_fits* is omitted, the script derives it by replacing the
         ``.fits`` suffix of *input_fits* with ``_SCI.fits``.
       """
       # Default input filename
       input_fits = 'image.fits'

       # ``sys.argv`` contains the script name at index 0.
       # Any extra arguments are interpreted as filenames.
       if len(sys.argv) > 1:
           input_fits = sys.argv[1]

       # Derive a default output name: ``myfile.fits`` → ``myfile_SCI.fits``
       output_fits = input_fits.replace('.fits', '_SCI.fits')
       if len(sys.argv) > 2:
           output_fits = sys.argv[2]

       # Perform the extraction
       extract_layers(input_fits, output_fits)

   # ------------------------------------------------------------------
   # ``python extract_layers.py`` will execute ``main()``
   # ------------------------------------------------------------------
   if __name__ == "__main__":
       main()

--------------------------------------------------------------------
4.  Function reference
----------------------

``extract_layers(input_fits, output_fits)``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Purpose** – Isolate the science data of a multi‑extension FITS file.
* **Parameters**
  * ``input_fits`` (``str``) – Path to the source file.
  * ``output_fits`` (``str``) – Destination path; existing files are overwritten.
* **Returns** – ``None`` (the side‑effect is the creation of ``output_fits``).
* **Raises**
  * :class:`OSError` if the file cannot be opened.
  * :class:`KeyError` if the ``'SCI'`` extension does not exist.
  * :class:`ValueError` if the header or data are malformed.

--------------------------------------------------------------------
5.  Command‑line interface
--------------------------

The script can be invoked directly from a terminal:

.. code-block:: bash

   # Use the default names (input = image.fits, output = image_SCI.fits)
   $ python extract_layers.py

   # Provide only an input file; output is auto‑derived
   $ python extract_layers.py my_observation.fits
   # → creates ``my_observation_SCI.fits``

   # Provide both input and output explicitly
   $ python extract_layers.py raw.fits calibrated_science.fits

**Explanation of the logic**

* ``sys.argv[0]`` – script name (ignored).
* ``sys.argv[1]`` – optional input FITS file.
* ``sys.argv[2]`` – optional output FITS file.
* If the user omits arguments, sensible defaults are used.

--------------------------------------------------------------------
6.  Example workflow
--------------------

Suppose you have a calibrated multi‑extension FITS file ``survey.fits`` that
contains the following HDUs:

* ``PRIMARY`` – empty header.
* ``SCI`` – the 2‑D image you want to analyse.
* ``ERR`` – variance array.
* ``DQ`` – data‑quality mask.

You only need the science image for a downstream pipeline that expects a
single‑extension FITS. Run:

.. code-block:: bash

   $ python extract_layers.py survey.fits
   # creates ``survey_SCI.fits``

You can now open the new file with any library that reads simple FITS images:

.. code-block:: python

   from astropy.io import fits
   data, header = fits.getdata('survey_SCI.fits', header=True)
   print(data.shape)
   print(header['OBJECT'])

--------------------------------------------------------------------
7.  Limitations & possible extensions
-------------------------------------

* **Hard‑coded extension name** – The function only works with an HDU named
  ``'SCI'``.  A future version could accept the extension name as an argument.
* **Unused import** – ``numpy`` is imported but not used.  Remove it or
  replace it with a small sanity check (e.g., ``np.isnan`` on the data).
* **Error handling** – The current script aborts with a traceback if the file
  does not exist or the ``SCI`` extension is missing.  Wrapping the I/O calls
  in ``try/except`` blocks would produce friendlier messages.
* **Logging** – Adding the :mod:`logging` module would allow users to enable
  verbose output without cluttering ``stdout``.

--------------------------------------------------------------------
8.  License
-----------

The code is released under the **MIT License**.  See the accompanying
``LICENSE`` file for the full text.

--------------------------------------------------------------------
9.  Contributing
----------------

Contributions are welcome!  Please fork the repository, make your changes,
and submit a pull request.  Ensure that:

* New functionality is documented in this ``.rst`` file.
* Unit tests are added under the ``tests/`` directory.
* The code follows the project's ``PEP 8`` style guide.

--------------------------------------------------------------------
10.  API reference (generated by Sphinx)
---------------------------------------

.. automodule:: extract_layers
   :members:
   :undoc-members:
   :show-inheritance:

--------------------------------------------------------------------
*Generated on 2025‑08‑12 by the documentation script.*
