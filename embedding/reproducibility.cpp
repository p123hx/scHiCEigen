
//
// Created by Hongxi on 16/12/2020.
//
#include <utility>
#include <fstream>
#include "reproducibility.h"
#include <vector>
#include <string>
#include <iostream>
#include <chrono>
#include <ctime>
#include "omp.h"
#include <mkl.h>
#include <thread>
#include <cctype>
#include <algorithm>
#include <cfloat>
#include <cmath>
#define EIGEN_USE_MKL_ALL
#include <Eigen/Core>
#define MAXD 1e12
using namespace std;
using namespace Eigen;
using namespace std::chrono;

MatrixXd euc_pdist_square(MatrixXd &x, int row, int col, double sigma) {
    vector<double> tmp;
    tmp.reserve(0.5*row*(row-1));
//NOt include #pragma omp parallel for
    {
        for (int i = 0; i < row-1; i++) {
            for (int j = i + 1; j < row; j++) {
                tmp.push_back((x.row(i) - x.row(j)).norm());
            }
        }
    }
    MatrixXd ans = MatrixXd::Zero(row, row);
//#pragma omp parallel for

    for (int i = 0; i < row; i++)
        for (int j = i+1; j < row; j++) {
            int ij=0.5*(2*row-1-i)*i+j-i-1;
            double tmpv = sqrt(2.0-2.0*exp(-sigma*tmp[ij]));
            ans(j, i) = tmpv;
            ans(i, j) = tmpv;
        }
//    cout << "\n ans: \n" << ans << endl;
    return ans;
}

vector<int> z_pos(MatrixXd & s1, MatrixXd & s2) {
    vector<int> z;
//#pragma omp parallel for
    for (int i = 0; i < s1.size(); i++) {
        if (!(s1(i) == 0 && s2(i) == 0)) z.push_back(i);
    }
    sort(z.begin(), z.end());
    return z;
}

MatrixXd zero_delete(MatrixXd& s, vector<int>& z) {
    VectorXd ans(z.size());
//#pragma omp parallel for
        for (int i = 0; i < z.size(); i++) {
            ans(i) = s(z[i]);
        }
    return ans;
}

vector<double>
pairwise_distance(vector<MatrixXd> &all_strata, string similarity_method, bool
print_time, double sigma, unsigned window_size
) {
    //#TODO Change thread
//    mkl_set_num_threads(56);
//    setNbThreads(56);

    transform(similarity_method.begin(), similarity_method.end(),
              similarity_method.begin(), ::tolower);
    chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;

    int n_cells = all_strata[0].rows(), n_bins = all_strata[0].cols();
    MatrixXd distance_mat(n_cells, n_cells);
    int n_strata = all_strata.size();
    t0 = high_resolution_clock::now();
    if (similarity_method == "hicrep") {
        int score_col = 0;
        MatrixXd weighted_std = MatrixXd::Zero(n_cells, n_bins);
//#pragma omp parallel for reduction(+:score_col)

        for (int i = 0; i < n_strata; i++) {
            int icol = all_strata[i].cols();
            auto mean = all_strata[i].rowwise().mean();
            all_strata[i].colwise() -= mean;
            MatrixXd std = (all_strata[i].array().square().rowwise().sum() /
                            icol).sqrt();
            weighted_std.col(i) = sqrt(n_bins - i) * std;
            score_col += icol;
        }


        t1 = high_resolution_clock::now();
        MatrixXd scores(n_cells, score_col);
        scores << all_strata[0], all_strata[1], all_strata[2],
                all_strata[3], all_strata[4], all_strata[5],
                all_strata[6], all_strata[7], all_strata[8],
                all_strata[9];

        t2 = high_resolution_clock::now();
        MatrixXd tmp1(n_cells, n_cells), tmp2(n_cells, n_cells);
        tmp1.noalias() = (scores * scores.transpose());
        tmp2.noalias() = (weighted_std * (weighted_std.transpose()));
        t3 = high_resolution_clock::now();
//#pragma omp parallel for schedule(guided) collapse(2)
            for (int i = 0; i < n_cells; i++)
                for (int j = 0; j < n_cells; j++) {
                    if (tmp2(i, j) == 0) {
                        distance_mat(i, j) = 0.0;
                        continue;
                    }
                    double tmpdiv = tmp1(i, j) / (tmp2(i, j));
                    if (tmpdiv >= 1) distance_mat(i, j) = 0.0;
                    else if (tmpdiv <= -1) distance_mat(i, j) = 2;
                    else distance_mat(i, j) = sqrt(2 - 2.0 * tmpdiv);
                }
        t4 = high_resolution_clock::now();


    } else if (similarity_method == "inner_product" or
               similarity_method == "innerproduct") {
        int score_col = 0;
        for (int i = 0; i < n_strata; i++) {
            int icol = all_strata[i].cols();
            auto mean = all_strata[i].rowwise().mean();
            all_strata[i].colwise() -= mean;
            MatrixXd std = (all_strata[i].array().square().rowwise().sum() /
                            icol).sqrt();
//            for (int j = 0; j < all_strata[i].rows(); j++) {
//                all_strata[i].row(j).array() /= std(j);
//            }

            for (int j = 0; j < all_strata[i].rows(); j++)
                for (int k = 0; k < all_strata[i].cols(); k++) {
                    all_strata[i](j, k) /= std(j);
                    if (isnan(all_strata[i](j, k))) all_strata[i](j, k) = 0.0;
                }
            score_col += icol;
        }
        MatrixXd scores(n_cells, score_col);
        scores << all_strata[0], all_strata[1], all_strata[2],
                all_strata[3], all_strata[4], all_strata[5],
                all_strata[6], all_strata[7], all_strata[8],
                all_strata[9];
        t1 = high_resolution_clock::now();
        MatrixXd inner; inner.noalias() = (scores * scores.transpose()) / score_col;
        t2 = high_resolution_clock::now();
//#pragma omp parallel for schedule(guided) collapse(2)
            for (int i = 0; i < n_cells; i++)
                for (int j = 0; j < n_cells; j++) {
                    double tmp = inner(i, j);
                    if (tmp >= 1)distance_mat(i, j) = 0.0;
                    else if (tmp <= -1) distance_mat(i, j) = 2.0;
                    else distance_mat(i, j) = sqrt(2 - 2.0 * tmp);
                }
        t4 = high_resolution_clock::now();
    } else if (similarity_method == "selfish") {
        int n_windows = n_bins / window_size;
        MatrixXd all_windows = MatrixXd::Zero(n_cells, n_windows);
#pragma omp parallel for
            for (int i = 0; i < n_strata; i++)
                for (int j = 0; j < n_windows; j++) {
                    all_windows.col(j) += all_strata[i].block(0, j * window_size,
                                                              n_cells,
                                                              window_size - i).rowwise
                            ().sum();
                }
        t1 = high_resolution_clock::now();
        int f_col = n_windows * (n_windows - 1) / 2;
        MatrixXd fingerprints = MatrixXd::Zero(n_cells, f_col);
#pragma omp parallel for
            for (int i = 0; i < n_windows; i++)
                for (int j = 0; j < n_windows - i - 1; j++) {
                    int k = (int) ((2 * n_windows - i - 1) * i * 0.5 + j);
                    for (int z = 0; z < n_cells; z++) {
                        fingerprints(z, k) = double(
                                all_windows(z, i) > all_windows(z, j));
                    }
                }
        distance_mat = euc_pdist_square(fingerprints, n_cells, f_col, sigma);
        t4 = high_resolution_clock::now();
    }
    else if (similarity_method == "old_hicrep") {
        MatrixXd similarity = MatrixXd::Zero(n_cells, n_cells);
        for (int i = 0; i < n_cells; i++)
            for (int j = i + 1; j < n_cells; j++) {
                RowVectorXd corrs(n_strata), weights(n_strata);
                for (int k = 0; k < n_strata; k++) {
                    MatrixXd stratum = all_strata[k];
                    MatrixXd s1 = stratum.row(i), s2 = stratum.row(j);
                    double std1 = sqrt(
                            (s1.array() - s1.mean()).square().sum() / s1.size());
                    double std2 = sqrt(
                            (s2.array() - s2.mean()).square().sum() / s2.size());
                    if (std1 == 0 || std2 == 0) {
                        weights(k) = 0;
                        corrs(k) = 0;
                    } else {
                        vector<int> zP = z_pos(s1, s2);
                        s1 = zero_delete(s1, zP);
                        s2 = zero_delete(s2, zP);

                        s1.array() -= s1.mean();
                        s2.array() -=s2.mean();
                        double std1_n = sqrt(
                                s1.array().square().sum() / s1.size());
                        double std2_n = sqrt(
                                s2.array().square().sum() / s2.size());
                        weights(k) = s1.size() * std1_n * std2_n;
                        //cout<<"weight: "<<weights(k)<<endl;
                        double tmp_nume = (s2.adjoint()*s1).value(), tmp_d = (s1.size() *
                                                                             std1_n *
                                                                             std2_n);
                        if (tmp_d == 0 && tmp_nume == 0) corrs(k) = 1.0;
                        else if (tmp_d == 0 && tmp_nume != 0) corrs(k) = MAXD;
                        else corrs(k) = tmp_nume / tmp_d;
                    }
                }
//                cout << "corrs: \n" << corrs << endl;
//                cout<<"weights: \n"<<weights<<endl;
                double tmp_v, tmp_n = corrs.dot(weights), tmp_d = weights.sum();
                if (tmp_n == 0 && tmp_d == 0) tmp_v = 1.0;
                else if (tmp_n != 0 && tmp_d == 0) tmp_v = MAXD;
                else tmp_v = tmp_n / tmp_d;
                double s = sqrt(2 - 2 * tmp_v);
                similarity(i, j) = s;
                similarity(j, i) = s;
            }
        distance_mat = similarity;
        t4 = high_resolution_clock::now();
    } else {
        throw "Method {0} not supported. Only \"inner_product\", \"HiCRep\", \"old_hicrep\" and \"Selfish\".";
    }

    std::chrono::duration<double, std::milli> duration1 = (t1 - t0);
    std::chrono::duration<double, std::milli> duration2 = (t4 - t1);
//    std::chrono::duration<double, std::milli> duration2 = (t2 - t1);
//    std::chrono::duration<double, std::milli> duration3 = (t3 - t2);
//    std::chrono::duration<double, std::milli> duration4 = (t4 - t3);
    std::chrono::duration<double, std::milli> duration_total = (t4 - t0);
    double tout = duration_total.count(), tout1 = duration1.count(), tout2 = duration2.count();
//    cout << "Time 1:" << duration1.count() << endl
//         << "Time 2:" << duration2.count() << endl
//         << "Time 3:" << duration3.count() << endl
//         << "Time 4:" << duration4.count() << endl
//    cout << "total: " << tout << endl;
//    cout << distance_mat << endl;
    vector<double> tv = {(double) tout, (double) tout1, (double) tout2};
    return tv;
}

