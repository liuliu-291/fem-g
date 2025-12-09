// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_MEMORY
#define H_GMSHFEM_MEMORY

#include <string>

namespace gmshfem::common
{
  enum class MessageType;

  template< MessageType T_Type >
  class Message;
} // namespace gmshfem::common

namespace gmshfem::common
{


  class Memory
  {
   private:
    unsigned long long _byte;

   public:
    Memory();
    Memory(const unsigned long long byte);
    Memory(const Memory &other);
    Memory(Memory &&other);
    ~Memory();

    Memory &operator=(const unsigned long long byte);
    Memory &operator=(const Memory &other);
    Memory &operator=(Memory &&other);

    void reset();
    std::string memory() const;
    bool isNull() const;

    unsigned long long byte() const;

    unsigned long long kilo() const;
    unsigned long long mega() const;
    unsigned long long giga() const;
    unsigned long long tera() const;

    unsigned long long kibi() const;
    unsigned long long mebi() const;
    unsigned long long gibi() const;
    unsigned long long tebi() const;

    Memory &operator+=(const Memory &other);
    Memory &operator-=(const Memory &other);

    Memory operator+(const Memory &other) const;
    Memory operator-(const Memory &other) const;

    bool operator==(const Memory &other) const;
    bool operator!=(const Memory &other) const;
    bool operator<(const Memory &other) const;
    bool operator>(const Memory &other) const;
    bool operator<=(const Memory &other) const;
    bool operator>=(const Memory &other) const;
  };

  template< MessageType T_Type >
  Message< T_Type > &operator<<(Message< T_Type > &message, const Memory &memory);


} // namespace gmshfem::common

#endif // H_GMSHFEM_MEMORY
