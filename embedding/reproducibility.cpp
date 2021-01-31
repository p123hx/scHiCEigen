
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
#include <Eigen/Core>
#include <thread>
#include <cctype>
#include <algorithm>
#include <cfloat>

using namespace std;
using namespace Eigen;
using namespace std::chrono;

vector<int> z_pos(MatrixXd s1, MatrixXd s2) {
    vector<int> z;
#pragma omp parallel for
    {
        for (int i = 0; i < s1.size(); i++) {
            if (!(s1(i) == 0 && s2(i) == 0)) z.push_back(i);
        }
    }
    sort(z.begin(), z.end());
    return z;
}

MatrixXd zero_delete(MatrixXd s, vector<int> z) {
    VectorXd ans(z.size());
#pragma omp parallel for
    {
        for (int i = 0; i < z.size(); i++) {
            ans(i) = s(z[i]);
        }
    }
    return ans;
}

double
pairwise_distance(vector<MatrixXd> all_strata, string similarity_method, bool
print_time, double sigma, unsigned window_size
) {
    //#TODO Change thread
    setNbThreads(16);
    transform(similarity_method.begin(), similarity_method.end(),
              similarity_method.begin(), ::tolower);
    chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;
    MatrixXd zscores;

    int n_cells = all_strata[0].rows(), n_bins = all_strata[0].cols();
    MatrixXd distance_mat(n_cells, n_cells);
    int n_strata = all_strata.size();
    t0 = high_resolution_clock::now();
    if (similarity_method == "hicrep") {
        int score_col = 0;
        MatrixXd weighted_std = MatrixXd::Zero(n_cells, n_bins);
#pragma omp parallel for reduction(+:score_col)

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
        int j;
#pragma omp parallel for
        {
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
        }
        t4 = high_resolution_clock::now();


    } else if (similarity_method == "inner_product" or
               similarity_method == "innerproduct") {
        MatrixXd similarity = MatrixXd::Ones(n_cells, n_cells);

    } else if (similarity_method == "old_hicrep") {
        MatrixXd similarity = MatrixXd::Zero(n_cells, n_cells);
        for (int i = 0; i < n_cells; i++)
            for (int j = i + 1; j < n_cells; j++) {
                VectorXd corrs(n_strata), weights(n_strata);
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
                        double std1_n = sqrt(
                                (s1.array() - s1.mean()).square().sum() / s1.size());
                        double std2_n = sqrt(
                                (s2.array() - s2.mean()).square().sum() / s2.size());
                        weights(k) = s1.size() * std1_n * std2_n;
                        double tmp_nume= (s1 * s2.transpose())(0), tmp_d =(s1.size() *std1_n *
                                                                           std2_n);
                        if(tmp_d==0 && tmp_nume==0) corrs(k)=1.0;
                        else if(tmp_d==0 &&tmp_nume!=0) corrs(k)=1e12;
                        else  corrs(k)=tmp_nume/tmp_d;

                    }
                }
//                cout << "corrs: \n" << corrs << endl;
//                cout<<"weights: \n"<<weights<<endl;
                cout<<(corrs * weights.transpose())(0)<<" / "<<weights.sum() <<endl;
                double tmp_v, tmp_n = (corrs * weights.transpose())(0),tmp_d=weights.sum();
                if(tmp_n ==0 && tmp_d==0) tmp_v=1.0;
                else if(tmp_n!=0 && tmp_d==0) tmp_v = 1e12;
                else tmp_v=tmp_n/tmp_d;
                double s = sqrt(2 - 2 * tmp_v);
                similarity(i, j) = s;
                similarity(j, i) = s;
            }
        distance_mat = similarity;
        t4 = high_resolution_clock::now();
    } else {
        throw "Method {0} not supported. Only \"inner_product\", \"HiCRep\", \"old_hicrep\" and \"Selfish\".";
    }

//    std::chrono::duration<double, std::milli> duration1 = (t1 - t0);
//    std::chrono::duration<double, std::milli> duration2 = (t2 - t1);
//    std::chrono::duration<double, std::milli> duration3 = (t3 - t2);
//    std::chrono::duration<double, std::milli> duration4 = (t4 - t3);
    std::chrono::duration<double, std::milli> duration_total = (t4 - t0);
    double tout = duration_total.count();
//    cout << "Time 1:" << duration1.count() << endl
//         << "Time 2:" << duration2.count() << endl
//         << "Time 3:" << duration3.count() << endl
//         << "Time 4:" << duration4.count() << endl
    cout << "total: " << tout << endl;
    cout << distance_mat << endl;
    return tout;
}


