// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Hilbert.h"

#include "OmpInterface.h"

#include <cmath>
#include <cstring>

namespace gmshfem::reorder
{


  template< unsigned int T_dim, unsigned int T_N >
  void Hilbert< T_dim, T_N >::computeGrayCode()
  {
    uint8_t mask = (unsigned int)~0 ^ ((unsigned int)~0 << T_dim);
    /* Generate the Gray code sequence:
   *    For 2D : 00 01 11 10
   *
   *          *------*
   *                 |
   *                 |
   *          *------*
   *
   *    For 3D : 000 001 011 010 110 111 101 100
   *
   *            *------*
   *            |     /
   *          *------*
   *            *------*
   *                  /
   *          *------*
   */
    uint8_t graycode[T_N];
    for(auto i = 0U; i < T_N; ++i) {
      graycode[i] = i ^ (i >> 1);
    }

    /* Generate the permuted Gray code sequence:
   *  start point: e
   */
    uint8_t transgraycode[T_N][T_dim][T_N];
    for(auto e = 0U; e < T_N; ++e) {
      for(auto d = 0U; d < T_dim; ++d) {
        for(auto i = 0U; i < T_N; ++i) {
          const unsigned int k = i ^ e;
          // Calculate the permuted Gray code by rotate k.
          transgraycode[e][d][i] =
            (k >> (d + 1)) | ((k << (T_dim - d - 1)) & mask);
        }
      }
    }

    /* Generate the inverse Gray code sequence:
   */
    uint8_t invgraycode[T_N];
    invgraycode[0] = 0;
    for(auto i = 1U; i < T_N; ++i) {
      invgraycode[i] = i;
      for(auto j = 1U; j < T_dim; ++j) {
        invgraycode[i] ^= i >> j;
      }
    }

    for(auto i = 0U; i < T_N; ++i) {
      for(auto j = 0U; j < T_dim; ++j) {
        for(auto k = 0U; k < T_N; ++k) {
          _transinvgraycode[i][j][k] = invgraycode[transgraycode[i][j][k]];
        }
      }
    }

    /* Generate the start point sequence:
   */
    uint8_t e[T_N];
    e[0] = 0;
    for(auto i = 1U; i < T_N; ++i) {
      e[i] = graycode[2 * ((i - 1) / 2)];
    }

    /* Generate the direction sequence:
   */
    uint8_t d[T_N];
    d[0] = 0;
    for(auto i = 1U; i < T_N; ++i) {
      if(i < T_N / 2)
        d[i] = (i % 2 == 0) ? graycode[i - 1] % T_dim : graycode[i] % T_dim;
      else
        d[i] = d[T_N - i - 1];
    }

    uint8_t nextE[T_N][T_dim][T_N];
    for(auto i = 0U; i < T_N; ++i) {
      for(auto j = 0U; j < T_dim; ++j) {
        for(auto k = 0U; k < T_N; ++k) {
          uint8_t l = e[_transinvgraycode[i][j][k]];
          nextE[i][j][k] =
            i ^ (((l << (j + 1)) & mask) | (l >> (T_dim - j - 1)));
        }
      }
    }

    uint8_t nextD[T_N][T_dim][T_N];
    for(auto i = 0U; i < T_N; ++i) {
      for(auto j = 0U; j < T_dim; ++j) {
        for(auto k = 0U; k < T_N; ++k) {
          nextD[i][j][k] = (j + d[_transinvgraycode[i][j][k]] + 1) % T_dim;
        }
      }
    }

    // ed = ooxxx in 3D and oxx in 2D with x the digits for e and o the digits for d
    for(auto i = 0U; i < T_N * T_dim; ++i) {
      for(auto j = 0U; j < T_N; ++j) {
        uint8_t e = i & ((unsigned int)~0 ^ ((unsigned int)~0 << T_dim));
        uint8_t d = i >> T_dim;
        if(d >= T_dim) continue;
        uint8_t nE = nextE[e][d][j];
        uint8_t nD = nextD[e][d][j];
        _nextED[i][j] = nD << T_dim | nE;
        _transinvgraycodeED[i][j] = _transinvgraycode[e][d][j];
      }
    }
  }

  template< unsigned int T_dim, unsigned int T_N >
  uint64_t Hilbert< T_dim, T_N >::getHilbertIndex(const uint32_t *const element, const int begBit) const noexcept
  {
    alignas(64) uint64_t h = 0;
    alignas(64) uint8_t ed = 0;
    for(int8_t m = begBit; m >= 0; m--) {
      alignas(64) uint8_t l = 0;
      for(int8_t n = T_dim - 1; n >= 0; n--) {
        l ^= ((element[n] >> m) & 0x1) << n;
      }
      h <<= T_dim;
      h ^= _transinvgraycodeED[ed][l];
      ed = _nextED[ed][l];
    }
    return h;
  }

  template< unsigned int T_dim, unsigned int T_N >
  uint64_t Hilbert< T_dim, T_N >::getMortonIndex(const uint32_t *const element, const int begBit) const noexcept
  {
    uint64_t h = 0;
    for(int8_t m = begBit; m >= 0; m--) {
      for(int8_t n = T_dim - 1; n >= 0; n--) {
        h |= ((element[n] >> m) & 0x1) << n;
      }
      h <<= T_dim;
    }
    return h;
  }

  template< unsigned int T_dim, unsigned int T_N >
  Hilbert< T_dim, T_N >::Hilbert() :
    _transinvgraycode(), _nextED(), _transinvgraycodeED()
  {
    computeGrayCode();
  }

  template< unsigned int T_dim, unsigned int T_N >
  Hilbert< T_dim, T_N >::~Hilbert()
  {
  }

  template< unsigned int T_dim, unsigned int T_N >
  void
  Hilbert< T_dim, T_N >::apply(std::vector< SortedEntity< T_dim > > &elements, const float min[T_dim], const float max[T_dim])
  {
    const unsigned long long nbrElements = elements.size();
    const unsigned int nbrBits = floor(64.0 / T_dim);
    const unsigned int twoPowNbrBits = 0x1 << (nbrBits - 1);

    alignas(64) float scale[T_dim];
    alignas(64) float scaleMin[T_dim];

    for(auto i = 0U; i < T_dim; ++i) {
      scale[i] = twoPowNbrBits / (max[i] - min[i]);
      scaleMin[i] = -min[i] * scale[i];
    }

    const uint8_t shift = 64 - 11;
    const unsigned int nbrBuckets = 2048;
    const unsigned int numThreads = omp::getNumThreads();

    SortedEntity< T_dim > *elementsCopy = new SortedEntity< T_dim >[nbrElements];
    uint64_t *histogramFull = new uint64_t[nbrBuckets * (numThreads + 1) + 1];
    uint64_t *histogram = &histogramFull[numThreads * nbrBuckets];

#pragma omp parallel num_threads(omp::getNumThreads())
    {
      const unsigned int myThreadID = omp::getThreadNum();
      uint64_t *histogramLocal = &histogramFull[myThreadID * nbrBuckets];
      memset(histogramLocal, 0, nbrBuckets * sizeof(uint64_t));

#pragma omp single
      {
        memset(histogram, 0, (nbrBuckets + 1) * sizeof(uint64_t));
      }

#pragma omp for
      for(auto i = 0ULL; i < nbrElements; ++i) {
        alignas(64) uint32_t u32[T_dim];
        for(auto j = 0U; j < T_dim; ++j) {
#ifdef FP_FAST_FMAF
          u32[j] = fma(elements[i].f[j], scale[j], scaleMin[j]);
#else
          u32[j] = elements[i].f[j] * scale[j] + scaleMin[j];
#endif
        }
        elements[i].ptrB = elements[i].ptrA;
        elements[i].h = getHilbertIndex(u32, nbrBits - 1);
        histogramLocal[elements[i].h >> shift]++;
      }

#pragma omp for
      for(auto i = 0U; i < nbrBuckets; ++i) {
        uint64_t sum = 0;
        for(auto j = 0U; j < numThreads + 1; ++j) {
          uint64_t tsum = histogramFull[j * nbrBuckets + i] + sum;
          histogramFull[j * nbrBuckets + i] = sum;
          sum = tsum;
        }
      }

#pragma omp single
      {
        uint64_t sum = 0;
        for(auto i = 0U; i < nbrBuckets + 1; ++i) {
          uint64_t tsum = sum + histogram[i];
          histogram[i] = sum;
          sum = tsum;
        }
      }

      for(auto i = 0U; i < nbrBuckets; ++i) {
        histogramLocal[i] += histogram[i];
      }

#pragma omp for
      for(auto i = 0ULL; i < nbrElements; ++i) {
        elementsCopy[histogramLocal[elements[i].h >> shift]++] = elements[i];
      }

#pragma omp for schedule(static, 1)
      for(auto i = 0U; i < nbrBuckets; ++i) {
        const unsigned int nbrElementsInBucket =
          histogram[i + 1] - histogram[i];
        const uint32_t mask = UINT32_C(0XFF);
        const unsigned int offset = histogram[i];
        SortedEntity< T_dim > *elementsLocal = &elements[offset];
        SortedEntity< T_dim > *elementsCopyLocal = &elementsCopy[offset];

        memset(histogramLocal, 0, nbrBuckets * sizeof(uint64_t));

        for(auto j = 0U; j < nbrElementsInBucket; ++j) {
          for(auto k = 0; k < 8; ++k) {
            histogramLocal[((elementsCopyLocal[j].h >> k * 8) & mask) + 256 * k]++;
          }
        }

        alignas(64) uint64_t sum[8] = {0};
        for(auto j = 0; j < 256; ++j) {
          for(auto k = 0; k < 8; ++k) {
            uint64_t tsum = sum[k] + histogramLocal[256 * k + j];
            histogramLocal[256 * k + j] = sum[k];
            sum[k] = tsum;
          }
        }

        for(auto k = 0; k < 7; ++k) {
          for(auto j = 0U; j < nbrElementsInBucket; ++j) {
            elementsLocal[histogramLocal[((elementsCopyLocal[j].h >> k * 8) & mask) + 256 * k]++] = elementsCopyLocal[j];
          }
          std::swap(elementsLocal, elementsCopyLocal);
        }
      }
    }

    delete[] elementsCopy;
    delete[] histogramFull;
  }


  template class Hilbert< 2, 4 >;
  template class Hilbert< 3, 8 >;


} // namespace gmshfem::reorder
