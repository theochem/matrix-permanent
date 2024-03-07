[![Python 3](http://img.shields.io/badge/python-3-blue.svg)](https://docs.python.org/3/)
[![gcc](https://img.shields.io/badge/-C++-blue?logo=cplusplus)](https://gcc.gnu.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/theochem/matrix-permanent/actions/workflows/pull_request.yml)
[![GNU GPLv3](https://img.shields.io/badge/license-%20%20GNU%20GPLv3%20-green?style=plastic)](https://www.gnu.org/licenses/gpl-3.0.en.html)


# Permanent

The permanent of a (square) matrix, like the determinant is a polynomial in the entries of the matrix. Unlike the determinant, the signatures of the permutations are not taken into account making the permanent much more difficult to compute because decomposition methods cannot be used.

The permanent commonly appears in problems related to quantum mechanics, and the most common brute-force combinatorial method has time complexity $\mathcal{O}(N!N)$, thus it is useful to look for more efficient algorithms. The two algorithms considered to be the fastest are one by Ryser (based on the inclusion-exclusion principle), and one by Glynn (based on invariant theory).

This library aims to solve the need for an efficient library that solves the permenent of a given matrix.

# Algorithms

`permanent.opt()`

Compute the permanent of a matrix using the best algorithm for the shape of the given matrix.

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

---

`permanent.combinatoric()`

Compute the permanent of a matrix combinatorically.

**Formula:**
```math
\text{per}(A) = \sum_{\sigma \in P(N,M)}{\prod_{i=1}^M{a_{i,{\sigma(i)}}}}
```

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

---

`permanent.glynn()`

**Formula:**

```math
\text{per}(A) = \frac{1}{2^{N-1}} \cdot \sum_{\delta}{
    \left(\sum_{k=1}^N{\delta_k}\right){\prod_{j=1}^N{\sum_{i=1}^N{\delta_i a_{i,j}}}}}
```

**Additional Information:**
The original formula has been generalized here to work with $M$-by-$N$ rectangular permanents with $M \leq N$ by use of the following identity (shown here for $M \geq N$):

```math
\text{per}\left(\begin{matrix}a_{1,1} & \cdots & a_{1,N} \\ \vdots & \ddots & \vdots \\ a_{M,1} & \cdots & a_{M,N}\end{matrix}\right) = \frac{1}{(M - N + 1)!} \cdot \text{per}\left(\begin{matrix}a_{1,1} & \cdots & a_{1,N} & 1_{1,N+1} & \cdots & 1_{1,M} \\ \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\ a_{M,1} & \cdots & a_{M,N} & 1_{M,N+1} & \cdots & 1_{M,M}\end{matrix}\right)
```

This can be neatly fit into the original formula by extending the inner sums over $\delta$ from $[1,M]$ to $[1,N]$:

```math
\text{per}(A) = \frac{1}{2^{N-1}} \cdot \frac{1}{(N - M + 1)!}\cdot \sum_{\delta}{
        \left(\sum_{k=1}^N{\delta_k}\right)
        \prod_{j=1}^N{\left(
            \sum_{i=1}^M{\delta_i a_{i,j}} + \sum_{i=M+1}^N{\delta_i}
        \right)}
    }
```

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

---

`permanent.ryser()`

**Formula:**

```math
\text{per}(A) = \sum_{k=0}^{M-1}{
        {(-1)}^k
        \binom{N - M + k}{k}
        \sum_{\sigma \in P(N,M-k)}{
            \prod_{i=1}^M{
                \sum_{j=1}^{M-k}{a_{i,{\sigma(j)}}}
            }
        }
    }
```

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

# Installation

The permanent package allows you to solve the permanent of a given matrix using the optimal algorithm for your matrix dimensions. You can either use the pre-defined parameters or fine tune them to your machine.

## Setting up your environment

1. Install Python on your machine. Depending on your operating system, the instructions may vary.

2. Install gcc on your machine. Depending on your operating system, the instructions may vary.

3. Create and activate a virtual environment for this project named `permanents`. One way to do this is with pip.

   ```bash
   pip install virtualenv
   virtualenv permanents
   ```

4. Activate the virtual environment.

   ```bash
   source permanents/bin/activate
   ```

5. Install Sphinx and other dependencies.

   ```bash
   pip install sphinx sphinx-rtd-theme sphinx-copybutton
   ```

6. Install Python dependencies.

   ```bash
   pip install numpy pandas scikit-learn
   ```

7. (Optional) Install Pytest if you wish to run tests.

   ```bash
   pip install pytest
   ```

Now that you have your environment set up and activated you are ready to compile the source code into an executable. Here you have two options - compile the code as is with the pre-defined parameters for algorithm swapping, **or** compile the code with machine specific tuning for algorithm swapping. _Note that machine specific tuning will run a series of tests. This will take anywhere from 10 minutes to 1 hour depending on your system._

## Option 1: Use given parameters

1. Compile the permanent code (natively for your CPU architecture).

   ```bash
   make BUILD_NATIVE=1
   ```

   **Note: if using M1 architecture, or want a portable build, simply run the following.**

   ```bash
   make
   ```

2. (Optional) Run tests on the algorithms.

   ```bash
   make test
   ```

3. Compile the website.

   ```bash
   cd docs && make html
   ```

4. Load the website.

   ```bash
   open build/html/index.html
   ```

## Option 2: Tune parameters

1. Compile the permanent code with the `tuning` flag.

   ```bash
   make RUN_TUNING=1
   ```

   **Note: it will take some time to run the tuning tests on your machine.**

2. (Optional) Run tests on the algorithms.

   ```bash
   make test
   ```

3. Compile the website.

   ```bash
   cd docs && make html
   ```

4. Load the website using your web browser.

   ```bash
   <browser> build/html/index.html
   ```

## Notes about the `Makefile`

The Makefile in this project is used to compile C and Python libraries and includes rules for installation, testing, and cleaning. Here's a breakdown of its sections:

1. Variables:

- `CXX`, `AR`, `PYTHON`: Define compiler, archiver, and Python executable.
- `CXXFLAGS`: Compiler flags including C++ version, warnings, debugging, optimization, and platform-specific options.

2. Conditional Compilation:

- `ifeq ($(shell uname -s),Darwin)`: Additional flags for macOS.
- `ifneq ($(BUILD_NATIVE),)`: Optimization flags if building for native architecture.
- `ifneq ($(RUN_TUNING),)`: Flag for runtime tuning.
- `ifeq ($(PREFIX),)`: Default installation prefix.

3. Targets:

- `all`, `c`, `python`: Phony targets for building all, C, or Python libraries.
- `install`: Installs C libraries and headers
- `test`: Runs tests using pytest.
- `clean`: Removes generated files.

4. File generation:

- `compile_flags.txt`: Generates compilation flags for clangd.
- `src/tuning.h`: Generates tuning parameters header file.

5. Compilation Rules:

- `permanent/permanent.so`: Compiles Python extension module.
- `src/libpermanent.o`: Compiles object code.
- `libpermanent.a`, `libpermanent.so`: Compiles static and shared C libraries respectively.

# License

This code is distributed under the GNU General Public License version 3 (GPLv3).
See <http://www.gnu.org/licenses/> for more information.

# Dependencies

The following programs/libraries are required to compile this package:

- [Python](http://python.org/) (≥3.6)
- [gcc](https://gcc.gnu.org/) (≥11.4)
- [Sphinx](https://pypi.org/project/Sphinx/) (≥7.2)
- [sphinx-rtd-theme](https://pypi.org/project/sphinx-rtd-theme/) (≥2.0)
- [NumPy](http://numpy.org/) (≥1.13)
- [pandas](https://pypi.org/project/pandas/) (≥2.2)
- [scikit-learn](https://pypi.org/project/scikit-learn/) (≥1.4)
- [Pytest](http://docs.pytest.org/en/latest/) (optional: to run tests)
