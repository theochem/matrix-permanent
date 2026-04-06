[![Python 3](http://img.shields.io/badge/python-3-blue.svg)](https://docs.python.org/3/)
[![gcc](https://img.shields.io/badge/-C++-blue?logo=cplusplus)](https://gcc.gnu.org/)
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

Compute the permanent of a matrix using an automatically selected algorithm. The library uses a polynomial logistic regression model (degree 4) trained on benchmarks to predict whether Ryser's or Glynn's algorithm will be faster for the given matrix dimensions.

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
The original formula has been generalized here to work with $M$-by-$N$ rectangular permanents with
$M \leq N$ by use of the following identity (shown here for $M \geq N$):

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

The `permanent`  package allows you to solve the permanent of a given matrix using the
**optimal algorithm** for your matrix dimensions.

## Quick Start (Recommended)

The easiest way to install and use the `permanent` package:

### From PyPI
```bash
pip install qc-permanent
```

### For Development
```bash
# Clone the repository and install with test dependencies
# Note: This project uses pyproject.toml, editable mode requires pip >= 21.3
pip install ".[test]"

# For editable mode (requires pip >= 21.3 with PEP 660 support):
# pip install --editable ".[test]"

# Run tests
pytest tests/
```

This will install the package with pre-set parameters that work well for most cases.

## Advanced Installation

For users who want to compile from source with custom optimizations or machine-specific tuning.

### Prerequisites

- **Python** ≥ 3.9
- **CMake** ≥ 3.23 (required for building from source)
- **C++ Compiler**: gcc ≥ 11.4 or equivalent
- **make**: System build tool (not a Python package)

### Installing Prerequisites

#### Ubuntu/Debian
```bash
# Install build tools and newer CMake
sudo apt-get update
sudo apt-get install build-essential

# Install CMake 3.23+ (Ubuntu's default may be older)
# Option 1: Via Kitware's APT repository
sudo apt install ca-certificates gpg wget

# If the kitware-archive-keyring package has not been installed previously, manually obtain a copy of our signing key:
test -f /usr/share/doc/kitware-archive-keyring/copyright ||
wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | gpg --dearmor - | sudo tee /usr/share/keyrings/kitware-archive-keyring.gpg >/dev/null

# Add kitware's repository to your sources list and update.
# For Ubuntu Focal Fossa (20.04):
echo 'deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ focal main' | sudo tee /etc/apt/sources.list.d/kitware.list >/dev/null

# If the kitware-archive-keyring package has not been installed previously, remove the manually obtained signed key to make room for the package:
test -f /usr/share/doc/kitware-archive-keyring/copyright ||
sudo rm /usr/share/keyrings/kitware-archive-keyring.gpg

# Install the kitware-archive-keyring package to ensure that your keyring stays up to date as we rotate our keys:
sudo apt install kitware-archive-keyring

# Finally we can install the cmake package
sudo apt install cmake
```

# Option 2: Via official CMake website
```bash
# Navigate to https://cmake.org/download/ and download version 3.23+ for Linux x86_64
# For example, for CMake 3.29.0:
wget https://cmake.org/files/v3.29/cmake-3.29.0-linux-x86_64.tar.gz

# Extract the archive
tar -xzf cmake-3.29.0-linux-x86_64.tar.gz

# Move to /opt or preferred location
sudo mv cmake-3.29.0-linux-x86_64 /opt/cmake

# Add to PATH
export PATH=/opt/cmake/bin:$PATH
# To make permanent, add the above line to ~/.bashrc
```

#### macOS
```bash
# Install Xcode Command Line Tools (includes make)
xcode-select --install

# Install CMake via Homebrew
brew install cmake
```

#### Conda (All Platforms)
```bash
# Install within your conda environment
conda install -c conda-forge make cmake compilers
```

### Setting up your environment

1. Create and activate a virtual environment:

   ```bash
   python -m venv permanents
   source permanents/bin/activate  
   ```

2. Install Python dependencies:

   ```bash
   pip install numpy pandas scikit-learn pytest
   ```

3. (Optional) For documentation:

   ```bash
   pip install sphinx sphinx-rtd-theme sphinx-copybutton
   ```

### Build Options

#### Option 1: Standard Build
```bash
# Basic build with default optimizations
make
```

#### Option 2: Native CPU Optimizations
```bash
# Optimize for your specific CPU architecture
make BUILD_NATIVE=1
```

**Note:** Use standard build for M1 Macs or if you need a portable build.

#### Option 3: Machine-Specific Tuning
```bash
# Run extensive benchmarks to find optimal algorithm thresholds for your machine
make PERMANENT_TUNE=1
```

**Important Notes:**
- Tuning will take 10-60 minutes depending on your system
- Creates `include/permanent/tuning.h` with custom parameters
- Creates `build/tuning.csv` with benchmark data

#### Verify Tuning Success
```bash
# Check if tuning generated custom parameters
diff include/permanent/tuning.h include/permanent/tuning.default.h

# If files differ, tuning succeeded!
# View generated parameters
cat include/permanent/tuning.h
```

### Running Tests
```bash
# After building
pytest tests/

# Or using make
make test
```

### Building Documentation
```bash
cd docs && make html
# Open build/html/index.html in your browser
```

### Troubleshooting

#### CMake version too old
If you get "CMake 3.23 or higher is required", see the Prerequisites section above for installation instructions.

#### make: command not found
`make` is a system tool, not a Python package. Install it using your system's package manager (see Prerequisites).

#### Build cache issues
```bash
# Clean and rebuild
make clean
# or
rm -rf build/
make
```

## Notes about the `Makefile`

The Makefile in this project is used to compile C and Python libraries and includes rules for
installation, testing, and cleaning. Here's a breakdown of its sections:

1. Variables:

- `CXX`, `AR`, `PYTHON`: Define compiler, archiver, and Python executable.
- `CXXFLAGS`: Compiler flags including C++ version, warnings, debugging, optimization, and
platform-specific options.

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
