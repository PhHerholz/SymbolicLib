#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <igl/read_triangle_mesh.h>
#include <iostream>
#include <fstream>
#include <filesystem>

#include "../../src/scalar/Symbolic.hpp"
#include "../../src/scalar/SymbolicUtilities.hpp"
#include "../../src/compiler/ComputeUnit.hpp"
#include "../../src/support/Timer.hpp"
#include "../../src/support/Utilities.hpp"

#include "../dataPath.hpp"

using namespace std;
struct SymIndex
{
  int index;
  int pattern;
};

inline int convertAccessPatternToUniqueID(const std::vector<int>& pos, int size) {
  // find the unique orientation id
  std::vector<int> positions;
  std::vector<int> positionsAfterErasure;
  for (int i = 0; i < size; i++) {
    positions.push_back(i);
  }
  for (int i = 0; i < pos.size(); i++) {
    int p = pos[i];
    for (int j = 0; j < positions.size(); j++) {
      if (positions[j] == p) {
        positionsAfterErasure.push_back(j);
        positions.erase(positions.begin() + j);
        break;
      }
    }
  }
  // now we just neeed to compute the unique id from those positions
  int id = 0;
  for (int i = 0; i < pos.size(); i++) {
    int p = positionsAfterErasure[i];
    int length = 1;
    for (int j = 0; j < pos.size() - i - 1; j++) {
      length *= size - j - 1 - i;
    }
    if (length <= 0) {
      length = 1;
    }
    id += length * p;
  }
  return id;
}


template <class TNum>
std::array<long, 2> buildCotan(Eigen::Matrix < TNum, -1, -1 > & V, Eigen::MatrixXi& F, Eigen::SparseMatrix<TNum>& L, Eigen::SparseMatrix<TNum>& M)
{
  std::array<long, 2> ret;

  std::vector<Eigen::Triplet<TNum>> triplets;
  triplets.reserve(F.rows() * 12);

  // buildMass(V, F, M);

  for (int i = 0; i < F.rows(); ++i)
  {
    Eigen::Matrix<TNum, 3, 1> e0 = V.row(F(i, 0)) - V.row(F(i, 1));
    Eigen::Matrix<TNum, 3, 1> e1 = V.row(F(i, 0)) - V.row(F(i, 2));
    const auto a = e0.cross(e1).norm();

    for (int j = 0; j < 3; ++j)
    {
      const int j1 = F(i, (j + 1) % 3);
      const int j2 = F(i, (j + 2) % 3);

      Eigen::Matrix<TNum, 3, 1> e0 = V.row(F(i, j)) - V.row(j1);
      Eigen::Matrix<TNum, 3, 1> e1 = V.row(F(i, j)) - V.row(j2);

      TNum w = e0.dot(e1) / a;

      triplets.emplace_back(j1, j1, -w); // a node value is updated
      triplets.emplace_back(j2, j2, -w); // a node value is updated
      triplets.emplace_back(j1, j2, w); // an edge value is updated
      triplets.emplace_back(j2, j1, w); // an edge value
    }
  }

  L.resize(V.rows(), V.rows());
  L.setFromTriplets(triplets.begin(), triplets.end());

  return ret;
}

int main(int argc, char* argv[])
{

  using namespace Sym;
  std::string meshFile = filePath() + "/d200k.off";
  ;

  Eigen::MatrixXi F;
  Eigen::MatrixXd V;

  igl::read_triangle_mesh(meshFile, V, F);

  // initializaiton of our stuffs
  std::vector<std::vector<SymIndex>> FACE;
  std::vector<std::vector<SymIndex>> VERTEX(V.rows());
  std::map<int, std::map<int, std::map<int, int>>> edge_to_face_map; // for each starting vertex, for each ending vertex, what is the triangle id and its pattern
  FACE.reserve(F.rows());
  for (int i = 0; i < F.rows(); i++) {
    std::vector<SymIndex> row;
    for (int j = 0; j < 3; j++) {
      row.push_back({F(i, j), -1});
    }
    FACE.push_back(row);
  }
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = j + 1; k < 3; k++) {
        std::vector<int> access_pattern = {j, k};
        int v0 = F(i, j);
        int v1 = F(i, k);
        if (v0 > v1) {
          int v2 = v0;
          v0 = v1;
          v1 = v2;
          access_pattern = {k, j};
        }
        if (edge_to_face_map.find(v0) == edge_to_face_map.end()) {
          edge_to_face_map[v0] = {};
        }
        if (edge_to_face_map[v0].find(v1) == edge_to_face_map[v0].end()) {
          edge_to_face_map[v0][v1] = {};
        }
        if (edge_to_face_map[v0][v1].find(i) == edge_to_face_map[v0][v1].end()) {
          int ap = convertAccessPatternToUniqueID(access_pattern, 3);
          edge_to_face_map[v0][v1][i] = ap;
        }
      }
    }
  }

  // initialize VERTEX2FACE
  std::vector<std::vector<SymIndex>> VERTEX2FACE(V.rows());
  for (int i = 0; i < F.rows(); i++) {
    for (int j = 0; j < 3; j++) {
      VERTEX2FACE[F(i, j)].push_back({i, j});
    }
  }
  // initialize EDGE and EDGE2FACE
  std::vector<std::vector<SymIndex>> EDGE;
  std::vector<std::vector<SymIndex>> EDGE2FACE;
  for (auto const& [v0, v1_map] : edge_to_face_map) {
    for (auto const& [v1, face_map] : v1_map) {
      EDGE.push_back({{v0, -1}, {v1, -1}});
      std::vector<SymIndex> face_and_patterns;
      for (auto const& [face_id, access_pattern] : face_map) {
        face_and_patterns.push_back({face_id, access_pattern});
      }
      EDGE2FACE.push_back(face_and_patterns);
    }
  }
  // std::cout << "Finished initialization\n";

  ;
  Eigen::MatrixXd VERTEX_POSITIONS = V;
  // end initialization

  Eigen::SparseMatrix<double> L, L2, M, M2;
  Timer t;
  buildCotan(V, F, L, M);
  t.printTime("buildCotan direct");


  // technically this is also initialization
  Eigen::SparseMatrix<double> OUTPUT_MATRIX;
  std::vector<std::vector<int>> VARYING_INPUT;
  OUTPUT_MATRIX.resize(V.rows(), V.rows());
  OUTPUT_MATRIX = L;

// here we have the output storage for this varying function0
  std::vector<std::vector<double>> VARYING_STORAGE_ACCESS_0 = std::vector<std::vector<double>>(FACE.size(), std::vector<double>(6));
// here we have the output storage for this varying function1
  std::vector<std::vector<double>> VARYING_STORAGE_ACCESS_1 = std::vector<std::vector<double>>(FACE.size(), std::vector<double>(3));
// here we have the output storage for this varying function2
  std::vector<std::vector<double>> VARYING_STORAGE_ACCESS_2 = std::vector<std::vector<double>>(FACE.size(), std::vector<double>(3));
// here we have the output storage for fixed function0
  std::vector<double> FIXED_STORAGE_ACCESS_0 = std::vector<double>(FACE.size());
  const std::vector<std::vector<int>> VARYING_INPUT_0 = {{2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1}, {1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0}, {2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1}, {0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2}, {1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0}, {0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2}};
  const std::vector<std::vector<int>> VARYING_INPUT_1 = {{2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1}, {2, 0, 2, 1, 2, 0, 2, 1, 2, 0, 2, 1}, {1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0}};
  const std::vector<std::vector<int>> VARYING_INPUT_2 = {{1, 2, 1, 0, 1, 2, 1, 0, 1, 2, 1, 0}, {0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2}, {0, 1, 0, 2, 0, 1, 0, 2, 0, 1, 0, 2}};


  t.reset();
  // ============================================================================================================

// VARYING FUNCTION 0
// THIS IS A FUNCTION OVER FACE

// here we have the actual for loop
  for (int FOR_INDEX_0 = 0; FOR_INDEX_0 < FACE.size(); FOR_INDEX_0++) {
    const std::vector<SymIndex>& INPUT_LEVEL_0 = FACE[FOR_INDEX_0];
    for (int pattern = 0; pattern < VARYING_INPUT_0.size(); pattern++) {
      VARYING_STORAGE_ACCESS_0[FOR_INDEX_0][pattern] = (((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][0]].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][1]].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][2]].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][3]].index, 0)))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][4]].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][5]].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][6]].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][7]].index, 2)))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][8]].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][9]].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][10]].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_0[pattern][11]].index, 1)))));
    }
  }

// VARYING FUNCTION 1
// THIS IS A FUNCTION OVER FACE
// here we have the actual for loop
  for (int FOR_INDEX_0 = 0; FOR_INDEX_0 < FACE.size(); FOR_INDEX_0++) {
    const std::vector<SymIndex>& INPUT_LEVEL_0 = FACE[FOR_INDEX_0];
    for (int pattern = 0; pattern < VARYING_INPUT_1.size(); pattern++) {
      VARYING_STORAGE_ACCESS_1[FOR_INDEX_0][pattern] = (((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][0]].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][1]].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][2]].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][3]].index, 0)))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][4]].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][5]].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][6]].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][7]].index, 2)))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][8]].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][9]].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][10]].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_1[pattern][11]].index, 1)))));
    }
  }

// VARYING FUNCTION 2
// THIS IS A FUNCTION OVER FACE
// here we have the actual for loop
  for (int FOR_INDEX_0 = 0; FOR_INDEX_0 < FACE.size(); FOR_INDEX_0++) {
    const std::vector<SymIndex>& INPUT_LEVEL_0 = FACE[FOR_INDEX_0];
    for (int pattern = 0; pattern < VARYING_INPUT_2.size(); pattern++) {
      VARYING_STORAGE_ACCESS_2[FOR_INDEX_0][pattern] = (((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][0]].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][1]].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][2]].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][3]].index, 0)))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][4]].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][5]].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][6]].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][7]].index, 2)))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][8]].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][9]].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][10]].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[VARYING_INPUT_2[pattern][11]].index, 1)))));
    }
  }

// FIXED FUNCTION 0
// THIS IS A FUNCTION OVER FACE
// here we have the actual for loop
  for (int FOR_INDEX_0 = 0; FOR_INDEX_0 < FACE.size(); FOR_INDEX_0++) {
    const std::vector<SymIndex>& INPUT_LEVEL_0 = FACE[FOR_INDEX_0];
    FIXED_STORAGE_ACCESS_0[FOR_INDEX_0] = sqrt((((((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 1)))) + (-1.000000 * ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 0)))))) * (((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 1)))) + (-1.000000 * ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 0))))))) + (((((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 2)))) + (-1.000000 * ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 1)))))) * (((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 1))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 2)))) + (-1.000000 * ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 1) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 1))))))) + (((-1.000000 * ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 2))))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 0))))) * ((-1.000000 * ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 0))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 2))))) + ((VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 2) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[1].index, 2))) * (VERTEX_POSITIONS(INPUT_LEVEL_0[0].index, 0) + (-1.000000 * VERTEX_POSITIONS(INPUT_LEVEL_0[2].index, 0)))))))));
  }

// OUTPUT FUNCTION 0
// THIS IS A FUNCTION OVER EDGE2FACE
// THIS FUNCTION REFERENCES FROM STORAGE EDGE
// THIS FUNCTION REFERENCES TO STORAGE FACE
  for (int FOR_INDEX_0 = 0; FOR_INDEX_0 < EDGE.size(); FOR_INDEX_0 ++) {
    double result = 0.0;
    const std::vector<SymIndex>& INPUT_LEVEL_0 = EDGE2FACE[FOR_INDEX_0];
    double INTERMEDIATE_0 = 0.0;

    for (int FOR_INDEX_1 = 0; FOR_INDEX_1 < INPUT_LEVEL_0.size(); FOR_INDEX_1++) {
      int pattern = INPUT_LEVEL_0[FOR_INDEX_1].pattern;
      double in1 = VARYING_STORAGE_ACCESS_0[INPUT_LEVEL_0[FOR_INDEX_1].index][pattern];
      double in2 = FIXED_STORAGE_ACCESS_0[INPUT_LEVEL_0[FOR_INDEX_1].index];
      INTERMEDIATE_0 += ( in1 / in2);

    }
    ;
    OUTPUT_MATRIX.coeffRef(EDGE[FOR_INDEX_0][0].index, EDGE[FOR_INDEX_0][1].index) = INTERMEDIATE_0;
  }

// OUTPUT FUNCTION 1
// THIS IS A FUNCTION OVER VERTEX2FACE
// THIS FUNCTION REFERENCES FROM STORAGE VERTEX
// THIS FUNCTION REFERENCES TO STORAGE FACE
  for (int FOR_INDEX_0 = 0; FOR_INDEX_0 < VERTEX.size(); FOR_INDEX_0 ++) {
    double result = 0.0;
    const std::vector<SymIndex>& INPUT_LEVEL_0 = VERTEX2FACE[FOR_INDEX_0];
    double INTERMEDIATE_0 = 0.0;

    for (int FOR_INDEX_1 = 0; FOR_INDEX_1 < INPUT_LEVEL_0.size(); FOR_INDEX_1++) {
      int pattern = INPUT_LEVEL_0[FOR_INDEX_1].pattern;
      double in0 = VARYING_STORAGE_ACCESS_1[INPUT_LEVEL_0[FOR_INDEX_1].index][pattern];
      double in1 = FIXED_STORAGE_ACCESS_0[INPUT_LEVEL_0[FOR_INDEX_1].index];
      double in2 = VARYING_STORAGE_ACCESS_2[INPUT_LEVEL_0[FOR_INDEX_1].index][pattern];
      INTERMEDIATE_0 += ((-1.000000 * (in0 / in1)) + (-1.000000 * ( in2 / in1)));

    }
    ;
    OUTPUT_MATRIX.coeffRef(FOR_INDEX_0, FOR_INDEX_0) = INTERMEDIATE_0;
  }
  t.printTime("Output");
  // ============================================================================================================
  // // double diff = 0;
  std::cout << L.norm() << std::endl;
  std::cout << "diff: " << (L - OUTPUT_MATRIX).norm() << std::endl;


  M2 = M;
  L2 = L;
  fill_n(L2.valuePtr(), L2.nonZeros(), .0);
  fill_n(M2.valuePtr(), M2.nonZeros(), .0);

  auto Vs = Sym::makeSymbolic(V, 0);
  Eigen::SparseMatrix<Sym::Symbolic> Ls, Ms;

  t.reset();
  buildCotan(Vs, F, Ls, Ls);
  t.printTime("build cotan");

  Sym::ComputeUnit<double> unit(Device(VecWidth(4), NumThreads(2)), Vs, Ls);
  t.printTime("build");

  unit.compile();
  t.printTime("compile");
  for (int i = 0; i < 100; i++) {
    unit.execute(V);
  }
  t.printTime("execute");
  unit.getResults(L2);

  std::cout << "residual: " << squaredError(L, L2) << endl;
  printDifferences(L, L2, 1e-8);

  return 0;
}
