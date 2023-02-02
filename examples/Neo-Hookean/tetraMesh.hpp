#include <vector>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
class tetraMesh
{
private:
    int num_nodes;
    int num_elements;
    Eigen::MatrixXd node_positions;                       // record the position for each node
    Eigen::MatrixXi element_indices;                      // record the index of node for each element
    std::vector<Eigen::Matrix<double, 3, 3>> dm_inverses; // store the D_m inverse for each tetrahedral
    std::vector<float> volumes;                           // store the volumes of each element
    void computeVolume(int ind);                          // compute the volume of a specific tetrahedron
    void computeDmInverse(int ind);                       // compute the Dm inverse of a specific tetrahedron
    double u;                                             // the u factor for material
    double lambda;                                        // the lambda factor for material
    std::vector<std::vector<int>> node_to_element_list;   // record for each node, the element list
    void allocateGlobalHessianSpace();                    // allocate the space for global hessian
    void addLocalHessianToGlobalHessian(int ind);         // assemble a local hessian and add it to global hessian
    std::vector<std::set<int>> connectivity;              // records the connectivity
public:
    Eigen::SparseMatrix<double> globalHessian; // store the global hessian of the matrix
    tetraMesh() : num_nodes(0), num_elements(0){};
    tetraMesh(const std::string &input);
    Eigen::Matrix<double, 12, 12> computeLocalHessian(int ind, double coeff);         // compute the 12 by 12 local hessian, here we only consider the neo hookean energy hessian
    void addSmallMovements();                                                         // add small movement to each node, we don't care about the inversion or anything
    void computeGlobalHessian();                                                      // compute the global hessian
    void computeSPMVOneRow(int ind, const Eigen::VectorXd &in, Eigen::VectorXd &out); // compute one row(actually 3 rows) of spmv, the row is dependent on ind
    // void computeNodeToElement();                                               // build the connectivity list
};