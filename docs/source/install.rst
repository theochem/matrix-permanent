Installation
============
The permanent package allows you to solve the permanent of a given matrix using the optimal algorithm for your matrix dimensions. You can either use the pre-defined parameters or fine tune them to your machine.

Setting up your environment
---------------------------

#. Install Python on your machine. Depending on your operating system, the instructions may vary.

#. Install gcc on your machine. Depending on your operating system, the instructions may vary.

#. Create and activate a virtual environment for this project named ``permanents``. One way to do this is with pip.

	.. code-block:: bash

		pip install virtualenv
		virtualenv permanents

#. Activate the virtual environment.

	.. code-block:: bash

		source ~/permanents/bin/activate

#. Install Sphinx and other dependencies.

	.. code-block:: bash

		pip install Sphinx sphinx-rtd-theme sphinx-copybutton

#. Install Python dependencies.

	.. code-block:: bash

		pip install numpy pandas scikit-learn

#. (Optional) Install Pyest if you wish to run tests.

	.. code-block:: bash

		pip install pytest


Now that you have your environment set up and activated you are ready to compile the source code into an executable. Here you have two options - compile the code as is with the pre-defined parameters for algorithm swapping, **or** compile the code with machine specific tuning for algorithm swapping. *Note that machine specific tuning will run a series of tests. This will take anywhere from 10 minutes to 1 hour depending on your system.*

Option 1: Use given parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#. Compile the permanent code.

	.. code-block:: bash

		make BUILD_NATIVE=1

	**Note: if using M1 architecture, or want a portable build, simply run the following.**

	.. code-block:: bash

		make

#. (Optional) Run tests on the algorithms.

	.. code-block:: bash

		make test

#. Compile the website.

	.. code-block:: bash

		cd docs && make html

#. Load the website.

	.. code-block:: bash

		<browser> build/html/index.html

Option 2: Tune parameters
~~~~~~~~~~~~~~~~~~~~~~~~~

#. Compile the permanent code with the ``tuning`` flag.

	.. code-block:: bash

		make RUN_TUNING=1

**Note: it will take some time to run the tuning tests on your machine.**

#. (Optional) Run tests on the algorithms.

	.. code-block:: bash

		make test

#. Compile the website.

	.. code-block:: bash

		cd docs && make html

#. Load the website.

	.. code-block:: bash

		<browser> build/html/index.html

Notes about the ``Makefile``
----------------------------

The Makefile in this project is used to compile C and Python libraries and includes rules for installation, testing, and cleaning. Here's a breakdown of its sections:

#. Variables:

* ``CXX``, ``AR``, ``PYTHON``: Define compiler, archiver, and Python executable.
* ``CXXFLAGS``: Compiler flags including C++ version, warnings, debugging, optimization, and platform-specific options.

#. Conditional Compilation:

* ``ifeq ($(shell uname -s),Darwin)``: Additional flags for macOS.
* ``ifneq ($(BUILD_NATIVE),)``: Optimization flags if building for native architecture.
* ``ifneq ($(RUN_TUNING),)``: Flag for runtime tuning.
* ``ifeq ($(PREFIX),)``: Default installation prefix.

#. Targets:

* ``all``, ``c``, ``python``: Phony targets for building all, C, or Python libraries.
* ``install``: Installs C libraries and headers
* ``test``: Runs tests using pytest.
* ``clean``: Removes generated files.

#. File generation:

* ``compile_flags.txt``: Generates compilation flags for clangd.
* ``src/tuning.h``: Generates tuning parameters header file.

#. Compilation Rules:

* ``permanent/permanent.so``: Compiles Python extension module.
* ``src/libpermanent.o``: Compiles object code.
* ``libpermanent.a, libpermanent.so``: Compiles static and shared C libraries respectively.
