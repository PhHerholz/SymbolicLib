#pragma warning(push, 0)
#include <Eigen/Dense>
#include <igl/read_triangle_mesh.h>
#pragma warning(pop)

#include "../../src/scalar/Symbolic.hpp"
#include "../../src/compiler/ComputeUnit.hpp"
#include "../../src/support/Timer.hpp"
#include "../../src/scalar/SymbolicDifferentiation.hpp"
#include "stretch_shear_forces.hpp"

#include "../dataPath.hpp"

int main(int argc, char* argv[])
{

  using namespace Sym;
  using namespace std;
  typedef float RealT;

  string meshFile = filePath() + "/cloth3.off";
  cout << "read " << meshFile << endl;

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;

  if (!igl::read_triangle_mesh(meshFile, V, F))
  {
    std::cout << "could not read file\n";
    return 0;
  }

  Eigen::MatrixXd X;
  X.resizeLike(V);
  X.setRandom();
  V += X;

  X.setRandom();
  X = 0.01 * X + V;

  StretchShear<double> ss;
  ss.init(800., 200., V, F);
  Eigen::VectorXd grad;
  vector<Eigen::Triplet<double>> tripH;

  Timer t;
  ss.compute_grad_hessian(X, grad, tripH, false);
  t.printTime("build values");

  Eigen::SparseMatrix<RealT> H(X.size(), X.size());
  H.setFromTriplets(tripH.begin(), tripH.end());
  t.printTime("build matrix");

  StretchShear<Symbolic> ssS;
  ssS.init(800., 200., V, F);

  Eigen::Matrix<Symbolic, -1, -1> XS;
  XS.resizeLike(X);
  setVariables(XS, 0);

  Timer tx;

  // direct compilation
  Eigen::Matrix<Symbolic, -1, 1> gradS;
  vector<Eigen::Triplet<Symbolic>> tripHS;
  ssS.compute_grad_hessian(XS, gradS, tripHS, false);
  Eigen::SparseMatrix<Symbolic> HS(X.size(), X.size());
  HS.setFromTriplets(tripHS.begin(), tripHS.end());

  //  ComputeUnit<RealT> unit(Device(UseCuda(), ThreadsPerBlock(64)), XS, HS);
  ComputeUnit<RealT> unit(Device(VecWidth(8), NumThreads(1)), XS, HS);
  unit.compile();

  Eigen::Matrix<RealT, -1, -1> Xf = X.cast<RealT>();
  unit.execute(Xf);

  tx.reset();
  for (int i = 0; i < 100; ++i)
    unit.execute(Xf);
  tx.printTime("100 runs");

  vector<RealT> H2(HS.nonZeros());
  vector<RealT> grad2(gradS.size());

  unit.getResults(H2);
  cout << "residual: " << squaredError(H, H2) << endl;
  ;
  unit.close();

  // automatic differentiation
  auto en = Symbolic(ADD, ssS.energy(XS));

  tx.reset();
  auto tripHS_ad = hessianSparse(en, XS);
  tx.printTime("construct Hessian (AD)");

  Eigen::SparseMatrix<Symbolic> HS_ad;
  toSparseMatrix(HS_ad, tripHS_ad, XS.size(), true);

  vector<RealT> H3(HS_ad.nonZeros());
  ComputeUnit<RealT> unitAD(Device(VecWidth(8), NumThreads(1)), XS, HS_ad);
  unitAD.compile();
  unitAD.execute(Xf);

  Timer t2;
  for (int i = 0; i < 100; ++i)
    unitAD.execute(Xf);
  t2.printTime("100 runs (AD)");

  unitAD.getResults(H3);
  cout << "residual (AD): " << squaredError(H, H3) << endl;

  return 0;
}
