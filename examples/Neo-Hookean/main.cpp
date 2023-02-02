#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <omp.h>

using namespace std;
int main() {
#pragma omp parallel
    {
        // if (omp_get_thread_num() == 0) {
        std::cout << "Num threads: " << omp_get_num_threads() << endl;
        // }
    }
    std::cout << "Compile\n";
}