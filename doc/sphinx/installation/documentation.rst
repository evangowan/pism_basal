.. include:: ../global.txt

.. _sec-install-documentation:

Rebuilding PISM documentation
=============================

You might want to rebuild the documentation from source, as PISM and its documentation
evolve together. These tools are required:

.. list-table::
   :header-rows: 1

   * - Tool
     - Comment
   * - Sphinx_
     - needed to rebuild this Manual
   * - sphinxcontrib.bibtex_
     - needed to rebuild this Manual and the documentation of PISM's Python bindings
       (below)
   * - LaTeX_
     - needed to rebuild the PDF version of this Manual, the |pism-browser|, and the
       documentation of PISM's Python bindings
   * - ``dvipng``
     - needed to rebuild the documentation of PISM's Python bindings
   * - Latexmk_
     - needed to rebuild the PDF version of this Manual
   * - doxygen_
     - required to rebuild the |pism-browser|
   * - graphviz_
     - required to rebuild the |pism-browser|

Note that if you install Sphinx using MacPorts_, you will install a version that
corresponds to your Python version, and its executables will have names with suffixes
corresponding to this version, e.g. ``sphinx-build-2.7`` rather than ``sphinx-build`` for
Python 2.7. You will want to set up aliases so that the standard names work as well. To do
this, run

.. code-block:: none

   sudo port select sphinx py27-sphinx

(replacing ``py27-sphinx`` with ``py26-sphinx`` for Python 2.6, etc.) If you opt not to do
this, you can tell CMake the name of your Sphinx executable using

.. code-block:: none

   cmake -DSPHINX_EXECUTABLE=sphinx-build-2.7 ...

for example.


Manual
------

To rebuild this manual, change to the PISM build directory and run

.. code-block:: bash

   make manual_html # for the HTML version of the manual
   make manual_pdf  # for the PDF version of the manual

The main page for this manual is then in ``doc/sphinx/html/index.html`` inside your build
directory.

The PDF manual will be in ``doc/sphinx/pism_manual.pdf`` in your build directory.

Source Code Browser
-------------------

To rebuild the |pism-browser|, change to the PISM build directory and run

.. code-block:: bash

   make browser

The main page for the documentation is then in ``doc/browser/html/index.html`` inside
your build directory.

Documentation for PISM???s Python bindings and inversion tools
------------------------------------------------------------

These bindings make scripting and interactive use of PISM possible, but many PISM users
will not need them. Installing them is required to use PISM for inversion of surface
velocities for basal shear stress and ice hardness. Building their documentation is
strongly-recommended before use.

To build the documentation of PISM's Python bindings, change to the PISM build directory
and run

.. code-block:: bash

   make pismpython_docs

The main page for the documentation is then in ``doc/pismpython/html/index.html`` inside
your build directory. The documentation build can take some time while it builds a large
number of small images from LaTeX formulas.

Re-building documentation without PISM's prerequisites
------------------------------------------------------

To build documentation on a system without PISM???s prerequisite libraries (such as MPI and
PETSc), assuming that PISM sources are in ``~/pism-stable``, do the following:

.. code-block:: bash

   cd ~/pism-stable
   mkdir doc-build # create a build directory
   cd doc-build
   cmake ../doc

then commands "``make manual_html``", "``make manual_pdf``" and others (see above) will
work as expected.
