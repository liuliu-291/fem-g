// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include "Memory.h"

#include "Message.h"
#include "instantiate.h"

#include <utility>

namespace gmshfem::common
{


  Memory::Memory() :
    _byte(0)
  {
  }

  Memory::Memory(const unsigned long long byte) :
    _byte(byte)
  {
  }

  Memory::Memory(const Memory &other) :
    _byte(other._byte)
  {
  }

  Memory::Memory(Memory &&other) :
    _byte(std::move(other._byte))
  {
  }

  Memory::~Memory()
  {
  }

  Memory &Memory::operator=(const unsigned long long byte)
  {
    _byte = byte;
    return *this;
  }

  Memory &Memory::operator=(const Memory &other)
  {
    _byte = other._byte;
    return *this;
  }

  Memory &Memory::operator=(Memory &&other)
  {
    _byte = std::move(other._byte);
    return *this;
  }

  void Memory::reset()
  {
    _byte = 0;
  }

  std::string Memory::memory() const
  {
    double copy = _byte;
    unsigned int unit = 0;
    while(copy >= 1024.) {
      copy /= 1024;
      unit++;
    }

    std::string str = std::to_string(copy);
    unsigned int pos = 0;
    while(str[pos] != '.') ++pos;
    if(copy < 10.) {
      str.erase(pos + 3);
    }
    else if(copy >= 10. && copy < 100.) {
      str.erase(pos + 2);
    }
    else {
      str.erase(pos);
    }
    str += " ";

    switch(unit) {
    case 0:
      str += "[B]";
      break;
    case 1:
      str += "[KiB]";
      break;
    case 2:
      str += "[MiB]";
      break;
    case 3:
      str += "[GiB]";
      break;
    case 4:
      str += "[TiB]";
      break;
    case 5:
      str += "[PiB]";
      break;
    case 6:
      str += "[EiB]";
      break;
    default:
      str += "[??]";
      break;
    }

    return str;
  }

  bool Memory::isNull() const
  {
    return _byte == 0;
  }

  unsigned long long Memory::byte() const
  {
    return _byte;
  }

  unsigned long long Memory::kilo() const
  {
    return _byte / 1000ULL;
  }

  unsigned long long Memory::mega() const
  {
    return _byte / (1000ULL * 1000ULL);
  }

  unsigned long long Memory::giga() const
  {
    return _byte / (1000ULL * 1000ULL * 1000ULL);
  }

  unsigned long long Memory::tera() const
  {
    return _byte / (1000ULL * 1000ULL * 1000ULL * 1000ULL);
  }

  unsigned long long Memory::kibi() const
  {
    return _byte / 1024ULL;
  }

  unsigned long long Memory::mebi() const
  {
    return _byte / (1024ULL * 1024ULL);
  }

  unsigned long long Memory::gibi() const
  {
    return _byte / (1024ULL * 1024ULL * 1024ULL);
  }

  unsigned long long Memory::tebi() const
  {
    return _byte / (1024ULL * 1024ULL * 1024ULL * 1024ULL);
  }

  Memory &Memory::operator+=(const Memory &other)
  {
    _byte += other._byte;
    return *this;
  }

  Memory &Memory::operator-=(const Memory &other)
  {
    _byte -= other._byte;
    return *this;
  }

  Memory Memory::operator+(const Memory &other) const
  {
    Memory memory(_byte + other._byte);
    return memory;
  }

  Memory Memory::operator-(const Memory &other) const
  {
    Memory memory(_byte - other._byte);
    return memory;
  }

  bool Memory::operator==(const Memory &other) const
  {
    return _byte == other._byte;
  }

  bool Memory::operator!=(const Memory &other) const
  {
    return _byte != other._byte;
  }

  bool Memory::operator<(const Memory &other) const
  {
    return _byte < other._byte;
  }

  bool Memory::operator>(const Memory &other) const
  {
    return _byte > other._byte;
  }

  bool Memory::operator<=(const Memory &other) const
  {
    return _byte <= other._byte;
  }

  bool Memory::operator>=(const Memory &other) const
  {
    return _byte >= other._byte;
  }

  template< MessageType T_Type >
  Message< T_Type > &operator<<(Message< T_Type > &message, const Memory &memory)
  {
    return (message << memory.memory());
  }

  INSTANTIATE_OPLL(Message< TEMPLATE_PARAM_1 > &, , 0, 5, MessageType, TEMPLATE_ARGS(MessageType::None, MessageType::Error, MessageType::Warning, MessageType::Info, MessageType::Debug), TEMPLATE_PARAMS(Message< TEMPLATE_PARAM_1 > &, const Memory &))


} // namespace gmshfem::common
