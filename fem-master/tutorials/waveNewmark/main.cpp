#include <gmsh.h>
#include <gmshfem/Formulation.h>
#include <gmshfem/Function.h>
#include <gmshfem/GmshFem.h>
#include <gmshfem/Matrix.h>
#include <gmshfem/Message.h>
#include <gmshfem/Vector.h>
#include <petsc.h>

using namespace gmshfem;
using namespace gmshfem::common;
using namespace gmshfem::problem;
using namespace gmshfem::domain;
using namespace gmshfem::field;
using namespace gmshfem::function;
using namespace gmshfem::post;
using namespace gmshfem::equation;
using namespace gmshfem::algebra;

namespace gmodel = gmsh::model;
namespace factory = gmsh::model::geo;

#define OMEGA_TAG 1
#define DOMEGA_TAG 2
#define EMITTER_TAG 3

void mesh()
{
  double h = 5e-3;
  int p0 = factory::addPoint(0., 0., 0., h);
  int p1 = factory::addPoint(1., 0., 0., h);
  int p2 = factory::addPoint(1., 1., 0., h);
  int p3 = factory::addPoint(0., 1., 0., h);
  int pe = factory::addPoint(0.5, 0.5, 0., h);
  int l0 = factory::addLine(p0, p1);
  int l1 = factory::addLine(p1, p2);
  int l2 = factory::addLine(p2, p3);
  int l3 = factory::addLine(p3, p0);

  int cl = factory::addCurveLoop({l0, l1, l2, l3});
  int s1 = factory::addPlaneSurface({cl});

  factory::synchronize();
  gmodel::mesh::embed(0, {pe}, 2, s1);

  gmodel::addPhysicalGroup(2, {s1}, OMEGA_TAG);
  gmodel::setPhysicalName(2, OMEGA_TAG, "omega");
  gmodel::addPhysicalGroup(1, {l0, l1, l2, l3}, DOMEGA_TAG);
  gmodel::setPhysicalName(1, DOMEGA_TAG, "domega");
  gmodel::addPhysicalGroup(0, {pe}, EMITTER_TAG);
  gmodel::setPhysicalName(0, EMITTER_TAG, "emitter");

  factory::synchronize();
  gmodel::mesh::generate();
  gmsh::write("mesh.msh");
}

void Vec2Std(Vec v, std::vector< double > &s)
{
  PetscInt N;
  VecGetSize(v, &N);
  PetscScalar *array;
  VecGetArray(v, &array);
  s.resize(N);
  for(int i = 0; i < N; i++) s[i] = std::real(array[i]);
  VecRestoreArray(v, &array);
}

double excitation(double time)
{
  double offset = 1.0;
  double t = time - offset;
  double tau = 0.125;
  double T = t * t / tau / tau;
  return (1. - T) * std::exp(-T / 2.);
}

int main(int argc, char **argv)
{
  GmshFem gmshFem(argc, argv);
  mesh();

  Domain omega(2, OMEGA_TAG);
  Domain domega(1, DOMEGA_TAG);
  Domain emitter(0, EMITTER_TAG);

  int order = 1;
  std::string gauss = "Gauss2";

  /* Helmholtz */
  Field< double, Form::Form0 > p("p", omega | domega, FunctionSpaceTypeForm0::HierarchicalH1, order);


  Formulation< double > formulation("wave");
  /* LHS */
  formulation.integral(dt2_dof(p), tf(p), omega, gauss);
  formulation.integral(grad(dof(p)), grad(tf(p)), omega, gauss);
  formulation.integral(dt_dof(p), tf(p), domega, gauss);
  /* RHS */
  formulation.integral(-1., tf(p), emitter, "Gauss0");

  formulation.pre();
  formulation.assemble(true);

  MatrixCRS< double > mass, damping, stiffness;
  formulation.getMass(mass);
  formulation.getDamping(damping);
  formulation.getStiffness(stiffness);

  Vector< double > source;
  formulation.getRHS(source);

  Mat mK = stiffness.getPetsc();
  Mat mC = damping.getPetsc();
  Mat mM = mass.getPetsc();

  double newmark_beta = 0.25;
  double newmark_gamma = 0.5;
  double dt = 1e-2;

  Mat mA, mIter1, mIter2;
  /* Make global system matrix */
  MatDuplicate(mM, MAT_COPY_VALUES, &mA);
  MatAXPY(mA, newmark_gamma * dt, mC, SAME_NONZERO_PATTERN);
  MatAXPY(mA, newmark_beta * dt * dt, mK, SAME_NONZERO_PATTERN);

  /* Make x_n update matrix */
  MatDuplicate(mM, MAT_COPY_VALUES, &mIter1);
  MatScale(mIter1, 2.0);
  MatAXPY(mIter1, -(1 - 2.0 * newmark_gamma) * dt, mC, SAME_NONZERO_PATTERN);
  MatAXPY(mIter1, -(0.5 + newmark_gamma - 2.0 * newmark_beta) * dt * dt, mK, SAME_NONZERO_PATTERN);

  /* Make x_{n-1} update matrix */
  MatDuplicate(mM, MAT_COPY_VALUES, &mIter2);
  MatScale(mIter2, -1.0);
  MatAXPY(mIter2, -(newmark_gamma - 1.0) * dt, mC, SAME_NONZERO_PATTERN);
  MatAXPY(mIter2, -(0.5 - newmark_gamma + newmark_beta) * dt * dt, mK, SAME_NONZERO_PATTERN);

  Vec vx_prev, vx_curr, vx_next;
  Vec vb_prev, vb_curr, vb_next;
  Vec vRhs, vSource;
  vSource = source.getPetsc();
  VecDuplicate(vSource, &vRhs); //veccopy
  VecDuplicate(vSource, &vx_prev);
  VecZeroEntries(vx_prev);
  VecDuplicate(vSource, &vx_curr);
  VecZeroEntries(vx_curr);
  VecDuplicate(vSource, &vx_next);
  VecZeroEntries(vx_next);

  KSP ksp;
  KSPCreate(PETSC_COMM_SELF, &ksp);
  KSPSetReusePreconditioner(ksp, PETSC_TRUE);
  KSPSetOperators(ksp, mA, mA);
  KSPSetType(ksp, "preonly");

  PC pc;
  KSPGetPC(ksp, &pc);
  PCSetType(pc, PCCHOLESKY);
  PCFactorSetMatSolverType(pc, "mumps");

  double current_time = 0.;
  int num_timesteps = 100;

  VecDuplicate(vSource, &vb_prev);
  VecScale(vb_prev, excitation(current_time - dt));
  VecDuplicate(vSource, &vb_curr);
  VecScale(vb_curr, excitation(current_time));
  VecDuplicate(vSource, &vb_next);
  VecScale(vb_next, excitation(current_time + dt));

  PostproMap< double, Degree::Degree0 > map("p");
  for(int n = 0; n < num_timesteps; n++) {
    msg::info << "** Time step " << n << msg::endl;

    MatMult(mIter1, vx_curr, vRhs);
    MatMultAdd(mIter2, vx_prev, vRhs, vRhs);
    VecAXPY(vRhs, dt * dt * newmark_beta, vb_next);
    VecAXPY(vRhs, dt * dt * (0.5 + newmark_gamma - 2.0 * newmark_beta), vb_curr);
    VecAXPY(vRhs, dt * dt * (0.5 - newmark_gamma + newmark_beta), vb_prev);

    // solve for x^{n+1}
    KSPSolve(ksp, vRhs, vx_next);

    // store solution back in the formulation
    std::vector< double > sol;
    Vec2Std(vx_next, sol);
    formulation.setSolutionIntoFields(sol);

    VecCopy(vx_curr, vx_prev);
    VecCopy(vx_next, vx_curr);

    current_time += dt;
    VecCopy(vb_curr, vb_prev);
    VecCopy(vb_next, vb_curr);
    VecCopy(vSource, vb_next);
    VecScale(vb_next, excitation(current_time + dt));

    if(n % 10 == 0) {
      map.append(p, omega | domega, n, n * dt);
    }
  }
  map.write("pos");
}
