//
// Created by Bean Juice on 08/01/2021.
//

#ifndef SCHICTOOLS_PROCESSING_UTILS_H
#define SCHICTOOLS_PROCESSING_UTILS_H
#include <cmath>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
MatrixXd matrix_operation(MatrixXd &mat,string operation);
MatrixXd convolution(MatrixXd &mat,int kernel_shape =3);
#endif //SCHICTOOLS_PROCESSING_UTILS_H
