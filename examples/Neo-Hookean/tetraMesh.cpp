#include "tetraMesh.hpp"
#include <iostream>
#include <fstream>
#include <stdlib.h> /* srand, rand */
// #include "AutoFlipSVD.hpp"
// #include "utils.hpp"
#include <set>
using namespace Eigen;
using namespace std;
#define my_print(x) std::cout << #x ": " << x << std::endl
tetraMesh::tetraMesh(const std::string& input)
{
    // this loading function is strictly designed
    // for .msh format with version 4.1
    // any different version may not work with this code
    // also, so far let's only consider the tetrahedron which has 4 nodes per element
    this->num_nodes = 1;
    this->num_elements = 1;
    ifstream infile(input);
    string info;
    float version_number;
    int file_type, data_size;
    int num_entities, min_node_index, max_node_index;
    int num_dim, place_holder1, param;
    string eol; // for end of the line

    infile >> info;                                     // this gives us begining of mesh format
    infile >> version_number >> file_type >> data_size; // this gives us mesh file version info
    infile >> info;                                     // this gives us end of mesh format

    infile >> info;                                                                  // this gives us beginning of nodes, nodes are vertices
    infile >> num_entities >> (this->num_nodes) >> min_node_index >> max_node_index; // this gives us node info
    this->node_positions.resize(3, this->num_nodes);                                 // resize our position matrix
    infile >> num_dim >> place_holder1 >> param >> (this->num_nodes);                // we will ignore this line for now
    assert(num_dim == 3);                                                            // assert it is a 3d mesh
    for (unsigned int i = 0; i < this->num_nodes; i++)                               // get the node vertex indices
        infile >> info;                                                              // ignore those infos because it's literally just 1 to N
    double x, y, z;                                                                  // placeholder to record the positions
    for (unsigned int i = 0; i < this->num_nodes; i++)                               // for loop
    {                                                                                // for loop
        infile >> x >> y >> z;                                                       // read positions
        this->node_positions.col(i) = Vector3d(x, y, z);                             // store positions
    }                                                                                // end of for loop
    infile >> info;                                                                  // this will give us end of nodes

    infile >> info;                                                                     // this will give us beginning of elements
    infile >> num_entities >> (this->num_elements) >> min_node_index >> max_node_index; // this gives us elements info
    infile >> num_dim >> place_holder1 >> param >> (this->num_elements);                // ignore this line for now
    assert(num_dim == 3);                                                               // assert this is a 3d mesh
    assert(param == 4);                                                                 // assert this is a tetrahedron
    this->element_indices.resize(param, this->num_elements);                            // resize our indices matrix
    this->volumes.resize(this->num_elements);                                           // resize storage for volumes
    this->dm_inverses.resize(this->num_elements);                                       // resize storage for dm_inverse
    this->node_to_element_list.resize(this->num_nodes);                                 // resize the storage for element to node list
    int index, n1, n2, n3, n4;                                                          // for recording the indices of each element
    for (unsigned int i = 0; i < this->num_elements; i++)                               // for loop
    {                                                                                   // for loop
        infile >> index >> n1 >> n2 >> n3 >> n4;                                        // read indices
        this->element_indices.col(i) = Vector4i(n1 - 1, n2 - 1, n3 - 1, n4 - 1);        // store the tetrahedron's indices
        this->node_to_element_list[n1 - 1].push_back(i);                                // store the correspondence
        this->node_to_element_list[n2 - 1].push_back(i);                                // store the correspondence
        this->node_to_element_list[n3 - 1].push_back(i);                                // store the correspondence
        this->node_to_element_list[n4 - 1].push_back(i);                                // store the correspondence
    }                                                                                   // end for loop
    infile >> info;                                                                     // this will give us end of elements

    cout << "Tet mesh read, number of nodes: " << this->num_nodes << ", number of elements: " << this->num_elements << endl;

    // here we set k and v to set for the material coefficients
    double k = 0.25;
    double v = 0.3;
    this->u = k / (2 * (1 + v));
    this->lambda = k * v / (1 + v) / (1 - 2 * v);

    // // compute the volume and dm inverse here
    // tbb::parallel_for(
    //     size_t(0), size_t(this->num_elements), [&](size_t i)
    //     { this->computeVolume(i);  
    //     this->computeDmInverse(i); });
    // this->addSmallMovements(); // add small displacements to each node
}

// this function computes the volume of a tetrahedron
// the input is the index of the tetrahedron
void tetraMesh::computeVolume(int i)
{
    // first we get the positions of the 4 nodes
    const Vector4i tetra_indices = this->element_indices.col(i);
    const Vector3d p0 = this->node_positions.col(tetra_indices[0]);
    const Vector3d p1 = this->node_positions.col(tetra_indices[1]);
    const Vector3d p2 = this->node_positions.col(tetra_indices[2]);
    const Vector3d p3 = this->node_positions.col(tetra_indices[3]);
    const Vector3d l01 = p1 - p0; // line from p0 to p1
    const Vector3d l02 = p2 - p0; // line from p0 to p2
    const Vector3d l03 = p3 - p0; // line from p0 to p3
    MatrixXd result_matrix(3, 3);
    result_matrix << l01, l02, l03;
    double volume = result_matrix.determinant() / 6;
    assert(volume >= 0);
    this->volumes[i] = volume;
}

// this function computes the D_m inverse in neo-hookean energy
// the input is the index of the tetrahedron
void tetraMesh::computeDmInverse(int i)
{
    const Vector4i tetra_indices = this->element_indices.col(i);
    const Vector3d p0 = this->node_positions.col(tetra_indices[0]);
    const Vector3d p1 = this->node_positions.col(tetra_indices[1]);
    const Vector3d p2 = this->node_positions.col(tetra_indices[2]);
    const Vector3d p3 = this->node_positions.col(tetra_indices[3]);
    MatrixXd Dm(3, 3);
    Dm.col(0) = p1 - p0;
    Dm.col(1) = p2 - p0;
    Dm.col(2) = p3 - p0;
    this->dm_inverses[i] = Dm.inverse();
}

// // this function adds a small displacement to all nodes
// // used so that there's some stretching energy
// void tetraMesh::addSmallMovements()
// {
//     tbb::parallel_for(size_t(0), size_t(this->num_nodes), [&](size_t i)
//         {
//             // add random number to positions
//             this->node_positions(0, i) += (double(rand()) / (double(RAND_MAX) + 1.0) - 0.5) / 1000.0;
//             this->node_positions(1, i) += (double(rand()) / (double(RAND_MAX) + 1.0) - 0.5) / 1000.0;
//             this->node_positions(2, i) += (double(rand()) / (double(RAND_MAX) + 1.0) - 0.5) / 1000.0; });
// }

// Matrix<double, 12, 12> tetraMesh::computeLocalHessian(int ind, double coeff)
// {
//     Eigen::Matrix<double, 12, 12> hessian;
//     // this hessian computation only considers 3d tetrahedron
//     const Vector4i tetra_indices = this->element_indices.col(ind);
//     const Matrix<double, 3, 3> A = this->dm_inverses[ind];
//     const Vector3d p0 = this->node_positions.col(tetra_indices[0]);
//     const Vector3d p1 = this->node_positions.col(tetra_indices[1]);
//     const Vector3d p2 = this->node_positions.col(tetra_indices[2]);
//     const Vector3d p3 = this->node_positions.col(tetra_indices[3]);

//     Matrix<double, 3, 3> F;
//     IPC::AutoFlipSVD<Matrix<double, 3, 3>> svd;

//     // here we redo the svd and compute the values of u, v, and sigma
//     Matrix<double, 3, 3> Xt;
//     Xt.col(0) = p1 - p0;
//     Xt.col(1) = p2 - p0;
//     Xt.col(2) = p3 - p0;
//     // my_print(A);
//     // my_print(Xt);
//     F = Xt * A;
//     svd.compute(F, ComputeFullU | ComputeFullV); // values are stored in svd
//     // my_print(F);
//     // my_print(svd.singularValues());

//     // compute derivatives
//     Eigen::Matrix<double, 9, 9> wdP_div_dF;
//     const double w = coeff * (this->volumes[ind]);
//     compute_dP_div_dF(svd, this->u, this->lambda, wdP_div_dF, w, true);
//     // my_print(wdP_div_dF);
//     // my_print(w);

//     Eigen::Matrix<double, 12, 9> wdP_div_dx;
//     dF_div_dx_mult<9>(wdP_div_dF.transpose(), A, wdP_div_dx, false);
//     dF_div_dx_mult<12>(wdP_div_dx.transpose(), A, hessian, true);
//     double square = (hessian * hessian.transpose()).norm();
//     double diff = (hessian - hessian.transpose()).norm();
//     // my_print(square);
//     // my_print(diff);
//     assert(diff <= 0.00001);

//     return hessian;
// }

// void tetraMesh::allocateGlobalHessianSpace()
// {
//     this->connectivity.resize(num_nodes);                                                                                      // records for each node, the connectivity
//     for (unsigned int i = 0; i < this->num_elements; i++)                                                                      // for loop
//     {                                                                                                                          // for loop
//         Vector4i indices = this->element_indices.col(i);                                                                       // get the indices
//         vector<int> sorted_node_indices({indices(0), indices(1), indices(2), indices(3)});                                     // sort the indice
//         sort(sorted_node_indices.begin(), sorted_node_indices.end());                                                          // sort the indices
//         for (int j = 0; j < 3; j++)                                                                                            // begin for loop to store indices
//         {                                                                                                                      // begin for loop to store indices
//             this->connectivity[sorted_node_indices[j]].insert(sorted_node_indices.begin() + j + 1, sorted_node_indices.end()); // put the rest of indices in the connectivity list
//         }                                                                                                                      // end of for loop
//     }
//     // here we allocate the space for global hessian
//     vector<Triplet<double>> global_hessian_allocation_triplets;          // the triplet for allocating the hessian matrix;
//     global_hessian_allocation_triplets.reserve(this->num_elements * 16); // allocate some space first
//     for (unsigned int i = 0; i < this->num_nodes; i++)
//     {
//         // first allocate the diagonal elements, which is 9 by 9
//         for (unsigned int m = 0; m < 3; m++)
//         {
//             for (unsigned int n = 0; n < 3; n++)
//             {
//                 global_hessian_allocation_triplets.push_back(Triplet<double>(i * 3 + m, i * 3 + n, 1.0));
//             }
//         }
//         // now allocate the rest connecting elements
//         // for each connecting element, it will allocate 2 3 by 3 block, one for row one for column
//         for (int neighbor_node_index : this->connectivity[i]) // neighbor_node_index is guaranteed to be larger than i
//         {
//             for (unsigned int m = 0; m < 3; m++)
//             {
//                 for (unsigned int n = 0; n < 3; n++)
//                 {
//                     global_hessian_allocation_triplets.push_back(Triplet<double>(i * 3 + m, neighbor_node_index * 3 + n, 1.0));
//                     global_hessian_allocation_triplets.push_back(Triplet<double>(neighbor_node_index * 3 + n, i * 3 + m, 1.0));
//                 }
//             }
//         }
//     }
//     this->globalHessian.resize(this->num_nodes * 3, this->num_nodes * 3);                                                      // set the size of global hessian
//     this->globalHessian.setFromTriplets(global_hessian_allocation_triplets.begin(), global_hessian_allocation_triplets.end()); // set from triplets
// }

// void tetraMesh::addLocalHessianToGlobalHessian(int ind)
// {
//     Matrix<double, 12, 12> local_hessian = this->computeLocalHessian(ind, 0.1);
//     Vector4i indices = this->element_indices.col(ind);
//     for (unsigned int m = 0; m < 4; m++)
//     {
//         int first_ind = indices(m) * 3;
//         for (unsigned int n = 0; n < 4; n++)
//         {
//             int second_ind = indices(n) * 3;
//             Matrix<double, 3, 3> block = local_hessian.block<3, 3>(m * 3, n * 3);
//             // add this block to global hessian
//             for (unsigned int i = 0; i < 3; i++)
//             {
//                 for (unsigned int j = 0; j < 3; j++)
//                 {
//                     this->globalHessian.coeffRef(first_ind + i, second_ind + j) += block(i, j);
//                 }
//             }
//         }
//     }
// }

// void tetraMesh::computeGlobalHessian()
// {
//     if (this->globalHessian.rows() == 0)
//     {
//         // we have not allocated anything for global hessian, assemble it once now
//         this->allocateGlobalHessianSpace();
//     }

//     // set zero to the global hessian

//     for (unsigned int i = 0; i < this->globalHessian.nonZeros(); i++)
//     {
//         this->globalHessian.valuePtr()[i] = 0;
//     }

//     for (unsigned int i = 0; i < this->num_elements; i++)
//     {
//         this->addLocalHessianToGlobalHessian(i);
//     }
//     // here we check if hessian is nonzero and symmetric
//     my_print((this->globalHessian).norm());
//     SparseMatrix<double> hessian_transpose = (this->globalHessian).transpose();
//     my_print(((this->globalHessian) - hessian_transpose).norm());
// }

// void tetraMesh::computeSPMVOneRow(int ind, const Eigen::VectorXd &in, Eigen::VectorXd &out)
// {
//     for (int i = 0; i < this->node_to_element_list[ind].size(); i++)
//     {
//         // get the element that contains this node
//         int element_ind = this->node_to_element_list[ind][i];
//         Vector4i node_indices = this->element_indices.col(element_ind);
//         Matrix<double, 12, 12> local_hessian = this->computeLocalHessian(element_ind, 0.1);
//         int node_position_in_element = 0;
//         for (int j = 0; j < 4; j++)
//         {
//             if (ind == node_indices(j))
//             {
//                 node_position_in_element = j;
//                 break;
//             }
//         }
//         // extract the 3 by 12 matrix
//         // multiply each of the 3 by 3 block
//         // my_print(out.block(ind * 3, 0, 3, 1));
//         for (int j = 0; j < 4; j++)
//         {
//             out.block(ind * 3, 0, 3, 1) += local_hessian.block(node_position_in_element * 3, j * 3, 3, 3) * in.block(node_indices(j) * 3, 0, 3, 1);
//         }
//         // my_print(value_place_holder);
//         // my_print(out.block(ind * 3, 0, 3, 1));
//         // break;
//     }
// }