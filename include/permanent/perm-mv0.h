#if !defined HAVE_PERM_MV0_H__
#define HAVE_PERM_MV0_H__
// This file is part of the FXT library.
// Copyright (C) 2010, 2011, 2012, 2014, 2019, 2023 Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the main directory.

#include <cstdlib>

#include "swap.h"

// define to update d[0] with each step:
#define PERM_MV0_UPDATE_D0  // default is on

// whether to use arrays instead of pointers:
#define PERM_MV0_FIXARRAYS  // small speedup; default off

class perm_mv0
// Inverse permutations corresponding to falling factorial numbers.
// CAT algorithm based on mixed radix Gray code
//   for the factorial number system (falling base).
{
 public:
#if !defined PERM_MV0_FIXARRAYS
  std::size_t *d_;  // mixed radix digits with radix = [n-1, n-2, n-3, ..., 2]
  std::size_t *x_;  // permutation
#else
  std::size_t d_[64];
  std::size_t x_[64];
#endif
  std::size_t ect_;  // counter for easy case
  std::size_t n_;    // permutations of n elements

  perm_mv0(const perm_mv0 &) = delete;
  perm_mv0 &operator=(const perm_mv0 &) = delete;

 public:
  explicit perm_mv0(std::size_t n)
  // Must have n>=2
  {
    n_ = n;
#if !defined PERM_MV0_FIXARRAYS
    d_ = new std::size_t[n_];
    x_ = new std::size_t[n_];
#endif
    d_[n - 1] = 1;  // sentinel (must be nonzero)
    first();
  }

  ~perm_mv0()
  {
#if !defined PERM_MV0_FIXARRAYS
    delete[] d_;
    delete[] x_;
#endif
  }

  const std::size_t *data() const { return x_; }

  void first()
  {
    for (std::size_t k = 0; k < n_; ++k) x_[k] = k;
    for (std::size_t k = 0; k < n_ - 1; ++k) d_[k] = 0;
    ect_ = 0;
  }

  bool next()
  {
    if (++ect_ < n_)  // easy case
    {
#if defined PERM_MV0_UPDATE_D0
      d_[0] = ect_;
#endif
      swap2(x_[ect_], x_[ect_ - 1]);
      return true;
    } else {
      ect_ = 0;
#if defined PERM_MV0_UPDATE_D0
      d_[0] = ect_;
#endif

      std::size_t j = 1;
      std::size_t m1 = n_ - 2;  // nine in falling factorial base
      while (d_[j] == m1)       // find digit to increment
      {
        d_[j] = 0;
        --m1;
        ++j;
      }

      if (j == n_ - 1) return false;  // current permutation is last

      const std::size_t dj = d_[j];
      d_[j] = dj + 1;

      // element at d[j] moves one position to the right:
      swap2(x_[dj], x_[dj + 1]);

      {  // move n-j elements to end:
        std::size_t s = n_ - j, d = n_;
        do {
          --s;
          --d;
          x_[d] = x_[s];
        } while (s);
      }

      // fill in 0,1,2,..,j-1 at start:
      for (std::size_t k = 0; k < j; ++k) x_[k] = k;

      return true;
    }
  }
};
// -------------------------

// #undef PERM_MV0_UPDATE_D0
// #undef PERM_MV0_FIXARRAYS

#endif  // !defined HAVE_PERM_MV0_H__
