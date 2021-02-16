//
// Created by Bean Juice on 16/12/2020.
//

#ifndef SCHICTOOLS_REPRODUCIBILITY_H
#define SCHICTOOLS_REPRODUCIBILITY_H

#include <utility>
#include <Eigen/Dense>
#include <vector>
using namespace std;
using namespace Eigen;
double pearsoncoeff(MatrixXd x, MatrixXd y, double size);

MatrixXd zscore_prop(MatrixXd a, int axis);

vector<double>
pairwise_distance(vector<MatrixXd> &all_strata, string similarity_method,
                  bool print_time = false, double sigma = .5,
                  unsigned window_size = 10);

#endif //SCHICTOOLS_REPRODUCIBILITY_H
