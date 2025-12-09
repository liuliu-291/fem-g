// GmshFEM - Copyright (C) 2019-2021 Université de Liège
//
// See the LICENSE.txt file for license information. Please report all
// issues on https://gitlab.onelab.info/gmsh/fem/issues

#include <ColorMap.h>
#include <Formulation.h>
#include <Function.h>
#include <Post.h>
#include <PGFPlotsio.h>

#include "convergences.h"
#include "geo.h"

#include <numeric>

using gmshfem::equation::dof;
using gmshfem::equation::dt_dof;
using gmshfem::equation::dt2_dof;
using gmshfem::equation::tf;
using gmshfem::function::operator-;

void convergences(gmshfem::common::GmshFem &gmshfem, unsigned int &numTest)
{
  gmshfem::msg::info << ++numTest << ") Convergence" << gmshfem::msg::endl;

  try {
    gmshfem::msg::info << "Convergences H1:" << gmshfem::msg::endl;
    convergencesH1();
  } catch(...) { throw; }
  
  // TODO : finish it
}

static void meshGeo(gmshfem::domain::Domain &omega, const double lc, const int geometry) {
  gmsh::model::add("geometry");
  int p0 = gmsh::model::geo::addPoint(0., 0., 0., lc);
  int p1 = gmsh::model::geo::addPoint(1., 0., 0., lc);
  int tagOmega = gmsh::model::geo::addLine(p0, p1);
  int dim = 1;
  if (geometry != 0) {
    int p2 = gmsh::model::geo::addPoint(1., 1., 0., lc);
    int p3 = gmsh::model::geo::addPoint(0., 1., 0., lc);
    int l1 = gmsh::model::geo::addLine(p1, p2);
    int l2 = gmsh::model::geo::addLine(p2, p3);
    int l3 = gmsh::model::geo::addLine(p3, p0);
    int cl = gmsh::model::geo::addCurveLoop({tagOmega, l1, l2, l3});
    tagOmega = gmsh::model::geo::addPlaneSurface({cl});
    dim = 2;
    if(geometry >= 3) {
      std::vector<std::pair<int, int>> outDimTags;
      switch (geometry) {
      case 3:
        gmsh::model::geo::extrude({std::make_pair(2, tagOmega)}, 0., 0., 0.5, outDimTags, std::vector<int>(1, 0.5 / lc), std::vector<double>(), false);
        break;
      case 4:
      case 5:
        gmsh::model::geo::extrude({std::make_pair(2, tagOmega)}, 0., 0., 0.5, outDimTags, std::vector<int>(1, 0.5 / lc), std::vector<double>(), true);
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
  omega = gmshfem::domain::Domain(dim, tagOmega);
}

static void linRegress(const std::vector<double> &x, const std::vector<double> &y, double &coefficientOfDetermination, double &slope) {
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

void convergencesH1()
{
  gmshfem::domain::Domain omega;
  std::vector< std::string >  elementType{"line", "triangle", "quadrangle", "tetrahedron", "hexahedron", "prism"};
  int FEorderMax = 3;
  int numberOfPoints = 4;
  const int maxNumOfGeo = 3; // 3D elements are too costly for the CI
  std::vector< double > slopes(maxNumOfGeo * FEorderMax);
  std::vector< double > coefficientsOfDetermination(maxNumOfGeo * FEorderMax);
  gmshfem::common::ColorMap colormap(gmshfem::common::ColorMapName::Neutral, FEorderMax);
  int it = 0;
  for(int iGeometry = 0; iGeometry < maxNumOfGeo; iGeometry++) {
    gmshfem::function::ScalarFunction<double> f;
    switch (iGeometry) {
    case 0:
      f = gmshfem::function::sin(10. * gmshfem::function::x<double>());
      break;
    case 1:
    case 2:
      f = gmshfem::function::sin(10. * gmshfem::function::x<double>()) + gmshfem::function::sin(10. * gmshfem::function::y<double>());
      break;
    case 3:
    case 4:
    case 5:
      f = gmshfem::function::sin(10. * gmshfem::function::x<double>()) + gmshfem::function::sin(10. * gmshfem::function::y<double>()) + gmshfem::function::sin(10. * gmshfem::function::z<double>());
      break;
    }
    gmshfem::common::PGFPlotsio latexPlot("H1errorPlotFor" + elementType[iGeometry]);
    for(int iFEorder = 1; iFEorder <= FEorderMax; iFEorder++) {
      std::string gauss = "Gauss" + std::to_string(5*iFEorder);
      double lc = 0.5;
      meshGeo(omega, lc, iGeometry);
      std::vector< double > l2RelativeError(numberOfPoints, 0);
      std::vector< double > meshRefinementLevel(numberOfPoints, 0);
      for (int i = 1; i <= numberOfPoints; i++) {
        lc /= 2.;
        gmsh::model::mesh::refine();
        {
          gmshfem::problem::Formulation< double > formulation("L2-Projection");
          gmshfem::field::Field< double, gmshfem::field::Form::Form0 > v("v", omega, gmshfem::field::FunctionSpaceTypeForm0::HierarchicalH1, iFEorder);

          formulation.integral(dof(v), tf(v), omega, gauss);
          formulation.integral(-f, tf(v), omega, gauss);

          formulation.pre();
          formulation.assemble();
          formulation.solve();

          double num = gmshfem::post::integrate(pow(abs(f - v), 2), omega, gauss);
          double den = gmshfem::post::integrate(pow(abs(f), 2), omega, gauss);
          l2RelativeError[i - 1] = log10(sqrt(num / den));
        }
        meshRefinementLevel[i - 1] = log10(lc);
      }
      gmshfem::common::Curve< double > curveForOneFEorder;
      double slopeForOneFEorder;
      double coefficientOfDeterminationForOneFEorder;
      linRegress(meshRefinementLevel, l2RelativeError, coefficientOfDeterminationForOneFEorder, slopeForOneFEorder);
      curveForOneFEorder.setName("FEorder/slope : " + std::to_string(iFEorder) + "/" + std::to_string(slopeForOneFEorder));
      curveForOneFEorder.setData(meshRefinementLevel, l2RelativeError);
      curveForOneFEorder.setMark("*");
      curveForOneFEorder.setColor(colormap(iFEorder));
      curveForOneFEorder.setLineStyle("dashed");
      latexPlot.addCurve(curveForOneFEorder);
      slopes[it] = slopeForOneFEorder;
      coefficientsOfDetermination[it] = coefficientOfDeterminationForOneFEorder;
      it++;
      gmsh::model::remove();
    }
    latexPlot.xLabel("log(h)");
    latexPlot.yLabel("log(Error)");
    latexPlot.legendPos("outer north east");
    latexPlot.printLegend(true);
    latexPlot.grid(true);
    latexPlot.write(true);
    latexPlot.close();
  }

  it = 0;
  for(int iGeometry = 0; iGeometry < maxNumOfGeo; iGeometry++) {
    gmshfem::msg::info << "Element type : " << elementType[iGeometry] << gmshfem::msg::endl;
    gmshfem::msg::info << "********************************" << gmshfem::msg::endl;
    for(int iFEorder = 1; iFEorder <= FEorderMax; iFEorder++) {
      gmshfem::msg::info << "FE order : " << iFEorder << gmshfem::msg::endl;
      gmshfem::msg::info << "Convergence rate (slope) : " << slopes[it] << gmshfem::msg::endl;
      gmshfem::msg::info << "coefficient of determination: " << coefficientsOfDetermination[it] << gmshfem::msg::endl;
      if(slopes[it] < ((double)iFEorder+1.) - 0.5 || slopes[it] > ((double)iFEorder+1.) + 0.5) {
        throw gmshfem::common::Exception("Error in " + std::string(__FILE__) + ", line " + std::to_string(__LINE__) + " -> function: convergencesH1");
      }
      it++;
    }
  }
}
