// GmshFEM - Copyright (C) 2019-2022, A. Royer, E. Béchet, C. Geuzaine, Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#ifndef H_GMSHFEM_HTMLTREE
#define H_GMSHFEM_HTMLTREE

#include "Color.h"
#include "Memory.h"
#include "NodeObject.h"
#include "io.h"

#include <fstream>
#include <vector>

namespace gmshfem::common
{


  struct HTMLNode {
    std::string name;
    std::string scalar;
    common::Memory localMemory;
    Degree degree;
    function::NodeType nodeType;
    bool constant;
    bool swappable;
    std::vector< HTMLNode > leaves;
  };

  class HTMLTree
  {
   private:
    std::ofstream _file;
    std::string _name;
    HTMLNode _head;
    common::Memory _peakMemory;

   public:
    HTMLTree();
    HTMLTree(const std::string &name, const std::string &path = "");
    ~HTMLTree();

    bool open(const std::string &name, const std::string &path = "");
    bool isOpen() const;
    void close();
    void write();

    void setPeakMemory(const common::Memory &peakMemory);

    HTMLNode *getHead();
  };


} // namespace gmshfem::common

#endif // H_GMSHFEM_HTMLTREE
