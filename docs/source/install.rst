Installation
============

The permanent package allows you to solve the permanent of a given
matrix using the optimal algorithm for your matrix dimensions. You can
either use the pre-defined parameters or fine tune them to your machine.

Installing from PyPI
--------------------

Simply run:

.. code:: bash

   pip install qc-permanent

Installing manually
-------------------

1. Install Python on your machine. Depending on your operating system,
   the instructions may vary.

2. Install gcc on your machine. Depending on your operating system, the
   instructions may vary.

3. Create and activate a virtual environment for this project named
   ``permanents``. One way to do this is with pip.

   .. code:: bash

      pip install virtualenv
      virtualenv permanents

4. Activate the virtual environment.

   .. code:: bash

      source permanents/bin/activate

5. Install ``qc-permanent``.

   .. code:: bash

      pip install .

   Optionally, install dependencies for building documentation
   (``doc``), running the tuning algorithm (``tune``), and/or running
   the tests (``test``) by specifying them in square brackets:

   .. code:: bash

      pip install `.[doc,tune,test]`

If you want to generate a machine-specific tuning header, preface the
``pip`` command with the corresponding environment variable like so:
``bash   PERMANENT_TUNE=ON pip install '.[tune]'``

This compiles the code with machine specific tuning for algorithm
swapping. *Note that machine specific tuning will run a series of tests.
This will take anywhere from 10 minutes to 1 hour depending on your
system.*

Using the C++ library
---------------------

The C++ library can be used by including the (CMake project for
``matrix-permanent``)[/CMakeLists.txt] in your own CMake project. The
(Makefile)[/Makefile] also acts as a convenience wrapper around the
CMake build for quickly compiling the C++ library.