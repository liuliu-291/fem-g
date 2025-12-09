#include <gmshfem/GmshFem.h>
#include <gmshfem/Domain.h>
#include <gmshfem/FieldInterface.h>

int main(int argc, char **argv)
{
  gmshfem::common::GmshFem gmshFem(argc, argv);
  
  bool derivative = false;
  gmshFem.userDefinedParameter(derivative, "derivative");
  std::string element = "line"; // "point", "line", "triangle", "quadrangle", "tetrahedron", "hexahedron", "prism", "pyramid"
  gmshFem.userDefinedParameter(element, "element");
  unsigned int order = 1;
  gmshFem.userDefinedParameter(order, "order");
  std::string form = ""; // "form0", "form1", "form2", "form3"
  gmshFem.userDefinedParameter(form, "form");
  std::string basisFunction = ""; // "isoParametric", "hierarchical", "p_isoParametric", "p_hierarchical", "d_isoParametric", "d_isoParametric", "dp_isoParametric", "dp_isoParametric"
  gmshFem.userDefinedParameter(basisFunction, "bf");
  double lc = 0.1;
  gmshFem.userDefinedParameter(lc, "lc");
  unsigned int orientation = 0;
  gmshFem.userDefinedParameter(orientation, "orientation");
  
  gmshfem::domain::Domain empty;
  
  if(form == "form0") {
    if(basisFunction == "isoParametric") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form0 > field("field", empty, gmshfem::functionSpaceH1::Lagrange, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form0 > field("field", empty, gmshfem::functionSpaceH1::HierarchicalH1, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else {
      gmshfem::msg::error << "Unknown basis function type '" << basisFunction << "' for differential form '" << form << "'" << gmshfem::msg::endl;
    }
  }
  else if(form == "form1") {
    if(basisFunction == "hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::HierarchicalHCurl, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "p_isoParametric") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::P_Lagrange, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "p_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::P_HierarchicalH1, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "d_isoParametric") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::D_Lagrange, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "d_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form1 > field("field", empty, gmshfem::functionSpaceHCurl::D_HierarchicalH1, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else {
      gmshfem::msg::error << "Unknown basis function type '" << basisFunction << "' for differential form '" << form << "'" << gmshfem::msg::endl;
    }
  }
  else if(form == "form2") {
    if(basisFunction == "p_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form2 > field("field", empty, gmshfem::functionSpaceHDiv::P_HierarchicalHCurl, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "d_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form2 > field("field", empty, gmshfem::functionSpaceHDiv::D_HierarchicalHCurl, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "dp_isoParametric") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form2 > field("field", empty, gmshfem::functionSpaceHDiv::DP_Lagrange, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "dp_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form2 > field("field", empty, gmshfem::functionSpaceHDiv::DP_HierarchicalH1, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else {
      gmshfem::msg::error << "Unknown basis function type '" << basisFunction << "' for differential form '" << form << "'" << gmshfem::msg::endl;
    }
  }
  else if(form == "form3") {
    if(basisFunction == "p_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form3 > field("field", empty, gmshfem::functionSpaceH0::P_Constant, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form3 > field("field", empty, gmshfem::functionSpaceH0::Constant, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else if(basisFunction == "dp_hierarchical") {
      gmshfem::field::Field< double, gmshfem::field::Form::Form3 > field("field", empty, gmshfem::functionSpaceH0::DP_HierarchicalHCurl, order);
      field.getFunctionSpace()->draw(derivative, element, orientation, lc);
    }
    else {
      gmshfem::msg::error << "Unknown basis function type '" << basisFunction << "' for differential form '" << form << "'" << gmshfem::msg::endl;
    }
  }
  else {
    gmshfem::msg::error << "Unknown differential form '" << form << "'" << gmshfem::msg::endl;
  }

  return 0;
}
