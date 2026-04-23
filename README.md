[![This project supports Python 3.10+](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org/downloads)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/theochem/matrix-permanent/actions/workflows/pull_request.yml)
[![GNU GPLv3](https://img.shields.io/badge/license-%20%20GNU%20GPLv3%20-green?style=plastic)](https://www.gnu.org/licenses/gpl-3.0.en.html)

# Permanent

The permanent of a (square) matrix, like the determinant is a polynomial in the entries of the
matrix. Unlike the determinant, the signatures of the permutations are not taken into account making
the permanent much more difficult to compute because decomposition methods cannot be used.

The permanent commonly appears in problems related to quantum mechanics, and the most common
brute-force combinatorial method has time complexity $\mathcal{O}(N!N)$, thus it is useful to look
for more efficient algorithms. The two algorithms considered to be the fastest are one by Ryser
(based on the inclusion-exclusion principle), and one by Glynn (based on invariant theory).

This library aims to solve the need for an efficient library that solves the permanent of a given
matrix.

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

$$
\text{per}(A) = \sum_{\sigma \in P(N,M)}{\prod_{i=1}^M{a_{i,{\sigma(i)}}}}
$$

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

---

`permanent.glynn()`

**Formula:**

$$
\text{per}(A) = \frac{1}{2^{N-1}} \cdot \sum_{\delta \in \left[\delta_1 = 1,~ \delta_2 \dots \delta_N=\pm1\right]}{
    \left(\sum_{k=1}^N{\delta_k}\right){\prod_{j=1}^N{\sum_{i=1}^N{\delta_i a_{i,j}}}}}
$$

**Additional Information:**
The original formula has been generalized here to work with $M$ by $N$ rectangular permanents with
$M \leq N$ by use of the following identity (shown here for $M \geq N$):

$$
\begin{aligned}
\text{per}\left(\begin{matrix}a_{1,1} & \cdots & a_{1,N} \\ \vdots & \ddots & \vdots \\ a_{M,1} & \cdots & a_{M,N}\end{matrix}\right) = \frac{1}{(M - N + 1)!} \cdot \text{per}\left(\begin{matrix}a_{1,1} & \cdots & a_{1,N} & 1_{1,N+1} & \cdots & 1_{1,M} \\ \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\ a_{M,1} & \cdots & a_{M,N} & 1_{M,N+1} & \cdots & 1_{M,M}\end{matrix}\right)
\end{aligned}
$$

This can be neatly fit into the original formula by extending the inner sums over $\delta$ from $[1,M]$ to $[1,N]$:

$$
\text{per}(A) = \frac{1}{2^{N-1}} \cdot \frac{1}{(N - M + 1)!}\cdot \sum_{\delta \in \left[\delta_1 = 1,~ \delta_2 \dots \delta_N=\pm1\right]}{
        \left(\sum_{k=1}^N{\delta_k}\right)
        \prod_{j=1}^N{\left(
            \sum_{i=1}^M{\delta_i a_{i,j}} + \sum_{i=M+1}^N{\delta_i}
        \right)}
    }
$$

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

---

`permanent.ryser()`

**Formula:**

$$
\text{per}(A) = \sum_{k=0}^{M-1}{
        {(-1)}^k
        \binom{N - M + k}{k}
        \sum_{\sigma \in P(N,M-k)}{
            \prod_{i=1}^M{
                \sum_{j=1}^{M-k}{a_{i,{\sigma(j)}}}
            }
        }
    }
$$

**Parameters:**

- `matrix`: `np.ndarray(M, N, dtype=(np.double|np.complex))`

**Returns:**

- `permanent`: `(np.double|np.complex)` - Permanent of matrix.

# Installation

The `permanent`  package allows you to solve the permanent of a given matrix using the
**optimal algorithm** for your matrix dimensions.

## Installing from PyPI

Simply run:
```bash
pip install qc-permanent
```

This will install the package with pre-set parameters with a good performance for most cases.
Advanced users can also **compile the code locally** and fine tune it for their specific
architecture. They can either use the pre-defined parameters or fine tune them to their machine.

## Installing manually

1. Install Python on your machine. Depending on your operating system, the instructions may vary.

2. Install gcc on your machine. Depending on your operating system, the instructions may vary.

3. Create and activate a virtual environment for this project named `permanents`. One way to do this
is with pip.

   ```bash
   pip install virtualenv
   virtualenv permanents
   ```

4. Activate the virtual environment.

   ```bash
   source permanents/bin/activate
   ```

5. Install `qc-permanent`.

   ```bash
   pip install .
   ```

   Optionally, install dependencies for building documentation (`doc`), running the tuning
   algorithm (`tune`), and/or running the tests (`test`) by specifying them in square brackets:
   ```bash
   pip install '.[doc,tune,test]'
   ```

  If you want to generate a machine-specific tuning header, preface the `pip` command with the
  corresponding environment variable like so:
   ```bash
  PERMANENT_TUNE=ON pip install '.[tune]'
   ```

  This compiles the code with machine specific tuning for algorithm swapping.
  _Note that machine specific tuning will run a series of tests.
  This will take anywhere from 10 minutes to 1 hour depending on your system._

## Using the C++ library

  The C++ library can be used by including the
  [CMake project for `matrix-permanent`](/CMakeLists.txt) in your own CMake project.
  The [Makefile](/Makefile) also acts as a convenience wrapper around the CMake
  build for quickly compiling the C++ library.

# License

This code is distributed under the GNU General Public License version 3 (GPLv3).
See <https://www.gnu.org/licenses/> for more information.
