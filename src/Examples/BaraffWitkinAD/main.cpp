#pragma warning(push, 0)  
#include <Eigen/Dense>
#include <igl/read_triangle_mesh.h>
#pragma warning(pop)  

#include "Symbolic.hpp"
#include "ComputeUnit.hpp"
#include "Timer.hpp"
#include "SymbolicDifferentiation.hpp"
#include "stretch_shear_forces.hpp"

int main(int argc, char *argv[]) {

    using namespace Sym;
    using namespace std;
    typedef double RealT;
    
    string meshFile = "../../data/cloth1.off";
    cout << "read " << meshFile << endl;
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    
    if(!igl::read_triangle_mesh(meshFile, V, F))  {
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
    
    Eigen::SparseMatrix<double> H(X.size(), X.size());
    H.setFromTriplets(tripH.begin(), tripH.end());
    t.printTime("build matrix");

    StretchShear<Symbolic> ssS;
    ssS.init(800., 200., V, F);
     
    Eigen::Matrix<Symbolic, -1, -1> XS;
    XS.resizeLike(X);
    setVariables(XS, 0);
   
    // direct compilation
    Eigen::Matrix<Symbolic, -1, 1> gradS;
    vector<Eigen::Triplet<Symbolic>> tripHS;
    ssS.compute_grad_hessian(XS, gradS, tripHS);
    Eigen::SparseMatrix<Symbolic> HS(X.size(), X.size());
    HS.setFromTriplets(tripHS.begin(), tripHS.end());
 
  //  ComputeUnit<RealT> unit(Device(UseCuda(), ThreadsPerBlock(64)), XS, HS);
    ComputeUnit<RealT> unit(Device(VecWidth(4), NumThreads(8)), XS, HS);
    unit.compile();

    Eigen::Matrix<RealT, -1, -1> Xf = X.cast<RealT>();
    unit.execute(Xf);
    
    Timer tx;
    for(int i = 0; i < 100; ++i) unit.execute(Xf);
    tx.printTime("100 runs");

    vector<RealT> H2(HS.nonZeros());
    vector<double> grad2(gradS.size());

    unit.getResults(H2);
    cout << "residual: " << squaredError(H, H2) << endl;;
    unit.close();

    // automatic differentiation
    auto en = Symbolic(ADD, ssS.energy(XS));
 
    tx.reset();
    auto tripHS_ad = hessianSparse(en, XS);
    tx.printTime("construct Hessian (AD)");
    
    Eigen::SparseMatrix<Symbolic> HS_ad;
    toSparseMatrix(HS_ad, tripHS_ad, XS.size(), true);
    
    vector<double> H3(HS_ad.nonZeros());
    ComputeUnit<double> unitAD(Device(VecWidth(4), NumThreads(8), DecompositionThreshold(6)), XS, HS_ad);
    unitAD.compile();
    unitAD.execute(X);
    
    Timer t2;
    for(int i = 0; i < 100; ++i) unitAD.execute(X);
    t2.printTime("100 runs (AD)");

    unitAD.getResults(H3);
    
    cout << "residual (AD): " << squaredError(H, H3) << endl;

    return 0;
}
