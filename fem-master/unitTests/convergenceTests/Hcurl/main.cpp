#include "CSVio.h"
#include "Formulation.h"
#include "Function.h"
#include "GmshFem.h"
#include "Message.h"
#include "PGFPlotsio.h"
#include "ColorMap.h"
#include "gmsh.h"
#include <numeric>
using namespace gmshfem::domain;
using namespace gmshfem::common;
using namespace gmshfem;
using namespace gmshfem::problem;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;

//***************************************
// LINEAR REGRESSION
//***************************************
void linRegress(const std::vector<double> &x, const std::vector<double> &y,
                double &coefficientOfDetermination, double &slope) {
  int n = x.size();
  double avgX = std::accumulate(x.begin(), x.end(), 0.0) / n;
  double avgY = std::accumulate(y.begin(), y.end(), 0.0) / n;
  double sxy = 0.0;
  double sxx = 0.0;
  double syy = 0.0;
  for (int i = 0; i < n; ++i) {
    sxy += (x[i] - avgX) * (y[i] - avgY);
    sxx += (x[i] - avgX) * (x[i] - avgX);
    syy += (y[i] - avgY) * (y[i] - avgY);
  }
  slope = sxy / sxx;
  coefficientOfDetermination = sxy * sxy / (sxx * syy);
}

//***************************************
// GEOMETRY
//***************************************
void meshGeo(Domain &omega, const double lc, const int geometry) {
  gmsh::model::add("geometry");
  gmsh::model::geo::addPoint(0., 0., 0., lc, 1);
  gmsh::model::geo::addPoint(1., 0., 0., lc, 2);
  gmsh::model::geo::addLine(1, 2, 1);
  int tagOmega = 1;
  int dim = 1;
  if (geometry != 0) {
    gmsh::model::geo::addPoint(1., 1., 0., lc, 3);
    gmsh::model::geo::addPoint(0., 1., 0., lc, 4);
    gmsh::model::geo::addLine(2, 3, 2);
    gmsh::model::geo::addLine(3, 4, 3);
    gmsh::model::geo::addLine(4, 1, 4);
    gmsh::model::geo::addCurveLoop({1, 2, 3, 4}, 1);
    gmsh::model::geo::addPlaneSurface({1}, 1);
    if (geometry == 2 || geometry == 4) {
      gmsh::model::geo::mesh::setRecombine(2, 1);
    }
    if (geometry < 3) {
      dim = 2;
    } else {
      std::vector<std::pair<int, int>> outDimTags;
      switch (geometry) {
      case 3:
        gmsh::model::geo::extrude({std::make_pair(2, tagOmega)}, 0., 0., 1.,
                                  outDimTags);
        break;
      case 4:
      case 5:
        gmsh::model::geo::extrude({std::make_pair(2, tagOmega)}, 0., 0., 1.,
                                  outDimTags, std::vector<int>(1, 1 / lc),
                                  std::vector<double>(), true);
        break;
      default:
        break;
      }
      tagOmega = outDimTags[1].second;
      dim = 3;
    }
  }
  gmsh::model::addPhysicalGroup(dim, {tagOmega}, 1);
  gmsh::model::setPhysicalName(dim, 1, "omega");
  gmsh::model::geo::synchronize();
  gmsh::model::mesh::generate();
  omega = Domain(dim, 1);
}

int main(int argc, char **argv) {
  GmshFem gmshFem(argc, argv);
  Domain omega;
  std::vector<std::string> elementType{
      "line", "triangle", "quadrangle", "tetrahedron", "hexahedron", "prism"};
  bool withRefine = true;
  gmshFem.userDefinedParameter(withRefine, "withRefine");
  int FEorderMax = 3;
  gmshFem.userDefinedParameter(FEorderMax, "FEorderMax");
  int numberOfPoints = 2;
  gmshFem.userDefinedParameter(numberOfPoints, "numberOfPoints");
  const int maxNumOfGeo = 6;
  std::vector<double> slopes(maxNumOfGeo * (FEorderMax+1));
  std::vector<double> coefficientsOfDetermination(maxNumOfGeo * (FEorderMax+1));
  ColorMap colormap(ColorMapName::Neutral, FEorderMax);
  int it = 0;
  for (int iGeometry = 0; iGeometry < maxNumOfGeo; iGeometry++) {
    VectorFunction<double> f;
    std::string gauss = "Gauss30";
    switch (iGeometry) {
    case 0:
      f = vector<double>(sin(10 * x<double>()), 0., 0.);
      break;
    case 1:
    case 2:
      f = vector<double>(sin(10 * x<double>()), sin(10 * y<double>()), 0.);
      break;
    case 3:
    case 4:
    case 5:
      gauss = "Gauss8";
      f = vector<double>(sin(10 * x<double>()), sin(10 * y<double>()),
                         sin(10 * z<double>()));
      break;
    }
    PGFPlotsio latexPlot("HcurlerrorPlotFor" + elementType[iGeometry]);
    for (int iFEorder = 0; iFEorder <= FEorderMax; iFEorder++) {
      double lc = 0.8;
      meshGeo(omega, lc, iGeometry);
      std::vector<double> l2RelativeError(numberOfPoints, 0);
      std::vector<double> meshRefinementLevel(numberOfPoints, 0);
      for (int i = 1; i <= numberOfPoints; i++) {
        if (withRefine) {
          lc = lc / 2.;
          gmsh::model::mesh::refine();
        } else {
          lc = 1. / (i * 2);
          meshGeo(omega, lc, iGeometry);
        }
        Formulation<double> formulation("L2-Projection");
        Field<double, Form::Form1> v(
            "v", omega, FunctionSpaceTypeForm1::HierarchicalHCurl, iFEorder);

        formulation.integral(dof(v), tf(v), omega, gauss);
        formulation.integral(-f, tf(v), omega, gauss);

        formulation.pre();
        formulation.assemble();
        formulation.solve();
        double numx = integrate(pow(abs(xComp(f) - xComp(v)), 2), omega, gauss);
        double denx = integrate(pow(abs(xComp(f)), 2), omega, gauss);
        double error = numx / denx;
        switch (iGeometry) {
        case 1:
        case 2:
          error += integrate(pow(abs(yComp(f) - yComp(v)), 2), omega, gauss) /
                   integrate(pow(abs(yComp(f)), 2), omega, gauss);
          break;
        case 3:
        case 4:
        case 5:
          error += integrate(pow(abs(yComp(f) - yComp(v)), 2), omega, gauss) /
                       integrate(pow(abs(yComp(f)), 2), omega, gauss) +
                   integrate(pow(abs(zComp(f) - zComp(v)), 2), omega, gauss) /
                       integrate(pow(abs(zComp(f)), 2), omega, gauss);
          break;
        }
          l2RelativeError[i-1] = log10(sqrt(error));
        meshRefinementLevel[i-1] = log10(lc);
      }
      Curve<double> curveForOneFEorder;
      double slopeForOneFEorder;
      double coefficientOfDeterminationForOneFEorder;
      linRegress(meshRefinementLevel, l2RelativeError,
                 coefficientOfDeterminationForOneFEorder, slopeForOneFEorder);
      curveForOneFEorder.setName("FEorder/slope : " + std::to_string(iFEorder) +
                                 "/" + std::to_string(slopeForOneFEorder));
      curveForOneFEorder.setData(meshRefinementLevel, l2RelativeError);
      curveForOneFEorder.setMark("*");
      curveForOneFEorder.setColor(colormap(iFEorder));
      curveForOneFEorder.setLineStyle("dashed");
      latexPlot.addCurve(curveForOneFEorder);
      slopes[it] = slopeForOneFEorder;
      coefficientsOfDetermination[it] = coefficientOfDeterminationForOneFEorder;
      it++;
    }
    latexPlot.xLabel("log(h)");
    latexPlot.yLabel("log(Error)");
    latexPlot.legendPos("outer north east");
    latexPlot.printLegend(true);
    latexPlot.grid(true);
    latexPlot.write(true);
    latexPlot.close();
  }
  msg::info << "********************************" << msg::endl;
  msg::info << "********************************" << msg::endl;
  it = 0;
  msg::info << " In this code, the developed conforming hierarchical "
               "finiteelement approximations for the space Hcurl is  validated."
            << msg::endl;
  msg::info << " To meet this objective, the L2  projection of a known "
               "function is performed on a domain and the relative  error "
               "between the known function and the FE solution is computed for "
               "different mesh densities and FE discretization orders"
            << msg::endl;
  for (int iGeometry = 0; iGeometry < maxNumOfGeo; iGeometry++) {
    msg::info << "Element type : " << elementType[iGeometry] << msg::endl;
    msg::info << "********************************" << msg::endl;
    for (int iFEorder = 0; iFEorder <= FEorderMax; iFEorder++) {
      msg::info << "FE order : " << iFEorder << msg::endl;
      msg::info << "Convergence rate (slope) : " << slopes[it] << msg::endl;
      msg::info << "coefficient of determination: "
                << coefficientsOfDetermination[it] << msg::endl;
      it++;
    }
  }
return 0;
}
