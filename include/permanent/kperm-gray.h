#if !defined HAVE_KPERM_GRAY_H__
#define HAVE_KPERM_GRAY_H__
// This file is part of the FXT library.
// Copyright (C) 2010, 2012, 2013, 2014, 2019, 2023 Joerg Arndt
// License: GNU General Public License version 3 or later,
// see the file COPYING.txt in the main directory.

#include <cstdlib>

#include "swap.h"

class kperm_gray
// Gray code for k-permutations of n elements.
// Same as: k-prefixes of permutations of n elements.
// Same as: arrangements of k out of n elements.
// CAT algorithm based on mixed radix Gray code
//   for the factorial number system (falling base).
{
 protected:
  std::size_t d_[64];      // mixed radix digits with radix = [n-1, n-2, ..., 2]
  std::size_t i_[64];      // directions
  std::size_t ix_[64];     // permutation (inverse perms in Trotter's order)
  std::size_t x_[64];      // inverse permutation (perms in Trotter's order)
  std::size_t n_;          // total of n elements
  std::size_t k_;          // prefix length: permutations of k elements
  std::size_t sw1_, sw2_;  // indices of elements swapped most recently

  kperm_gray(const kperm_gray &) = delete;
  kperm_gray &operator=(const kperm_gray &) = delete;

 public:
  explicit kperm_gray(std::size_t n)
  {
    n_ = n;
    k_ = n;
    // d_ = new std::size_t[n_];
    d_[n - 1] = -1UL;  // sentinel
    // i_ = new std::size_t[n_];
    // x_ = new std::size_t[n_];
    // ix_ = new std::size_t[n_];
    i_[n_ - 1] = 0;
    first(n_);
  }

  ~kperm_gray()
  {
    // delete [] i_;
    // delete [] d_;
    // delete [] x_;
    // delete [] ix_;
  }

  const std::size_t *data() const { return ix_; }
  const std::size_t *invdata() const { return x_; }
  void get_swap(std::size_t &s1, std::size_t &s2) const
  {
    s1 = sw1_;
    s2 = sw2_;
  }

  const std::size_t *mr_digits() const { return d_; }

  void first(std::size_t k)
  {
    k_ = k;
    for (std::size_t j = 0; j < n_ - 1; ++j) d_[j] = 0;
    for (std::size_t j = 0; j < n_ - 1; ++j) i_[j] = +1;
    for (std::size_t j = 0; j < n_; ++j) x_[j] = ix_[j] = j;
    sw1_ = n_ - 1;
    sw2_ = n_ - 2;
  }

  //    void last(std::size_t k)
  //    {
  //        first(k);  // k_ is set with call
  //        d_[n_-2] = 1;
  //        swap2(x_[n_-1], x_[n_-2]);
  //        swap2(ix_[n_-1], ix_[n_-2]);
  //        for (std::size_t j=0; j<n_-2; ++j)  i_[j] = -1UL;
  //        i_[n_-2] = +1;
  //    }

 private:
  void swap(std::size_t j, std::size_t im)
  {
    const std::size_t x1 = j;        // element j
    const std::size_t i1 = ix_[x1];  // position of j
    const std::size_t i2 = i1 + im;  // neighbor
    const std::size_t x2 = x_[i2];   // position of neighbor
    x_[i1] = x2;
    x_[i2] = x1;  // swap2(x_[i1], x_[i2]);
    ix_[x1] = i2;
    ix_[x2] = i1;  // swap2(ix_[x1], ix_[x2]);
    sw1_ = i1;
    sw2_ = i2;
  }

 public:
  bool next()
  {
    std::size_t j = 0;
    std::size_t m1 = n_ - 1;  // nine in falling factorial base
    std::size_t ij;
    while ((ij = i_[j])) {
      std::size_t im = i_[j];
      std::size_t dj = d_[j] + im;
      if (dj > m1)  // =^= if ( (dj>m1) || ((long)dj<0) )
      {
        i_[j] = -ij;
      } else {
        d_[j] = dj;
        swap(j, im);
        return true;
      }

      --m1;
      ++j;
      if (j >= k_) return false;
    }
    return false;
  }

  //    bool prev()
  //    {
  //        std::size_t j = 0;
  //        std::size_t m1 = n_ - 1;
  //        std::size_t ij;
  //        while ( (ij=i_[j]) )
  //        {
  //            std::size_t im = -i_[j];
  //            std::size_t dj = d_[j] + im;
  //            if ( dj>m1 )  // =^= if ( (dj>m1) || ((long)dj<0) )
  //            {
  //                i_[j] = -ij;
  //            }
  //            else
  //            {
  //                d_[j] = dj;
  //                swap(j, im);
  //                return true;
  //            }
  //
  //            --m1;
  //            ++j;
  //            if ( j>=k_ )  return false;
  //        }
  //        return false;
  //    }
};
// -------------------------

#endif  // !defined HAVE_KPERM_GRAY_H__
