// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_HILBERT
#define H_GMSHFEM_HILBERT

#include <cstdint>
#include <vector>

namespace gmshfem::reorder
{


  template< unsigned int T_dim >
  struct SortedEntity {
    union {
      struct {
        void *ptrB;
        float f[T_dim];
      };
      struct {
        void *ptrA;
        uint64_t h;
      };
    };
  };

  template< unsigned int T_dim, unsigned int T_N >
  class Hilbert
  {
   private:
    uint8_t _transinvgraycode[T_N][T_dim][T_N];
    alignas(64) uint8_t _nextED[T_N * (T_dim + T_dim % 2)][T_N];
    alignas(64) uint8_t _transinvgraycodeED[T_N * T_dim][T_N];

    void computeGrayCode();

   public:
    Hilbert();
    ~Hilbert();

    void apply(std::vector< SortedEntity< T_dim > > &elements, const float min[T_dim], const float max[T_dim]);

    uint64_t getHilbertIndex(const uint32_t *const element, const int begBit) const noexcept;
    uint64_t getMortonIndex(const uint32_t *const element, const int begBit) const noexcept;
  };


} // namespace gmshfem::reorder

#endif // H_GMSHFEM_HILBERT
