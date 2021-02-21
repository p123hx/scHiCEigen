//
// Created by Bean Juice on 08/01/2021.
//

#include "processing_utils.h"

using namespace std;
using namespace Eigen;

MatrixXd matrix_operation(MatrixXd &mat, string operation) {
    MatrixXd tmp;

    transform(operation.begin(), operation.end(), operation.begin(), ::tolower);
    if (operation == "convolution") {
        tmp = convolution(mat, 3);
    } else {
        throw "Not implemented";
    }


    return tmp;
}

MatrixXd convolution(MatrixXd &mat, int kernel_shape) {
//credit: https://github.com/chaowang15/fast-image-convolution-cpp/blob/master/src/convolution.cpp

    MatrixXd conv = MatrixXd::Ones(kernel_shape, kernel_shape) / pow
            (kernel_shape, 2.0);
    int h_dst = mat.rows(), w_dst = mat.cols();
    MatrixXd outMat = MatrixXd::Zero(h_dst, w_dst);

    double temp;
    int k, l;
    int low_k, high_k, low_l, high_l;

    int h_src = h_dst, w_src = w_dst;
    for (int i = 0; i < h_dst; ++i) {
        low_k = max(0, i - int((kernel_shape - 1.0) / 2.0));
        high_k = min(h_src - 1, i + int(kernel_shape / 2.0));
        for (int j = 0; j < w_dst; j++) {
            low_l = std::max(0, j - int((kernel_shape - 1.0) / 2.0));
            high_l = std::min(w_src - 1, j + int(kernel_shape / 2.0));
            temp = 0.0;
            for (k = low_k; k <= high_k; ++k) {
                for (l = low_l; l <= high_l; ++l) {
                    temp += mat(k, l) *
                            conv((i - k + int(kernel_shape / 2.0)),
                                 (j - l + int(kernel_shape / 2.0)));
                }
            }
            outMat(i, j) = temp;
        }
    }
    return outMat;

}