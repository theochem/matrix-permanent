Permanent
=========

.. toctree::
   :maxdepth: 2
   :caption: Contents

   install
   api

About
-----

The permanent of a (square) matrix, like the determinant is a polynomial
in the entries of the matrix. Unlike the determinant, the signatures of
the permutations are not taken into account making the permanent much
more difficult to compute because decomposition methods cannot be used.

The permanent commonly appears in problems related to quantum mechanics,
and the most common brute-force combinatorial method has time complexity
:math:`\mathcal{O}(N!N)`, thus it is useful to look for more efficient
algorithms. The two algorithms considered to be the fastest are one by
Ryser (based on the inclusion-exclusion principle), and one by Glynn
(based on invariant theory).

This library aims to solve the need for an efficient library that solves
the permenent of a given matrix.

Algorithms
----------

``permanent.opt()``

Compute the permanent of a matrix using the best algorithm for the shape
of the given matrix.

**Parameters:**

-  ``matrix``: ``np.ndarray(M, N, dtype=(np.double|np.complex))``

**Returns:**

-  ``permanent``: ``(np.double|np.complex)`` - Permanent of matrix.

--------------

``permanent.combinatoric()``

Compute the permanent of a matrix combinatorically.

**Formula:**

.. math::

   \text{per}(A) = \sum_{\sigma \in P(N,M)}{\prod_{i=1}^M{a_{i,{\sigma(i)}}}}

**Parameters:**

-  ``matrix``: ``np.ndarray(M, N, dtype=(np.double|np.complex))``

**Returns:**

-  ``permanent``: ``(np.double|np.complex)`` - Permanent of matrix.

--------------

``permanent.glynn()``

**Formula:**

.. math::

   \text{per}(A) = \frac{1}{2^{N-1}} \cdot \sum_{\delta}{
       \left(\sum_{k=1}^N{\delta_k}\right){\prod_{j=1}^N{\sum_{i=1}^N{\delta_i a_{i,j}}}}}

**Additional Information:** The original formula has been generalized
here to work with :math:`M`-by-:math:`N` rectangular permanents with
:math:`M \leq N` by use of the following identity (shown here for
:math:`M \geq N`):

.. math::

   \text{per}\left(\begin{matrix}a_{1,1} & \cdots & a_{1,N} \\ \vdots & \ddots & \vdots \\ a_{M,1} & \cdots & a_{M,N}\end{matrix}\right) = \frac{1}{(M - N + 1)!} \cdot \text{per}\left(\begin{matrix}a_{1,1} & \cdots & a_{1,N} & 1_{1,N+1} & \cdots & 1_{1,M} \\ \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\ a_{M,1} & \cdots & a_{M,N} & 1_{M,N+1} & \cdots & 1_{M,M}\end{matrix}\right)

This can be neatly fit into the original formula by extending the inner
sums over :math:`\delta` from :math:`[1,M]` to :math:`[1,N]`:

.. math::

   \text{per}(A) = \frac{1}{2^{N-1}} \cdot \frac{1}{(N - M + 1)!}\cdot \sum_{\delta}{
           \left(\sum_{k=1}^N{\delta_k}\right)
           \prod_{j=1}^N{\left(
               \sum_{i=1}^M{\delta_i a_{i,j}} + \sum_{i=M+1}^N{\delta_i}
           \right)}
       }

**Parameters:**

-  ``matrix``: ``np.ndarray(M, N, dtype=(np.double|np.complex))``

**Returns:**

-  ``permanent``: ``(np.double|np.complex)`` - Permanent of matrix.

--------------

``permanent.ryser()``

**Formula:**

.. math::

   \text{per}(A) = \sum_{k=0}^{M-1}{
           {(-1)}^k
           \binom{N - M + k}{k}
           \sum_{\sigma \in P(N,M-k)}{
               \prod_{i=1}^M{
                   \sum_{j=1}^{M-k}{a_{i,{\sigma(j)}}}
               }
           }
       }

**Parameters:**

-  ``matrix``: ``np.ndarray(M, N, dtype=(np.double|np.complex))``

**Returns:**

-  ``permanent``: ``(np.double|np.complex)`` - Permanent of matrix.

License
-------

This code is distributed under the GNU General Public License version 3
(GPLv3). See http://www.gnu.org/licenses/ for more information.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
