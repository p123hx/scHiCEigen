//
// Created by Bean Juice on 16/12/2020.
//

#ifndef SCHICTOOLS_LOAD_HIC_FILE_H
#define SCHICTOOLS_LOAD_HIC_FILE_H
#include <vector>
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

pair<set<string>, map<string, int> >
get_chromosome_lengths(const string &ref_str, const string &chromosomes,
                       int res);

pair<set<string>, map<string, int> >
get_chromosome_lengths(map<string, int> ref_map, string chromosomes,
                       int res);

pair<MatrixXd, vector<MatrixXd> > load_HiC(string file, map<string,
        int> genome_length,
                                          string
                                          format = "",
                                          int custom_format =
                                          0,
                                          int header = 0,
                                          string chromosome = "",
                                          int resolution = 10000,
                                          bool resolution_adjust = true,
                                          double map_filter = 0.,
                                          bool sparse = false,
                                          bool gzip = false,
                                          int keep_n_strata = 0,
                                          string
                                          operations = "");

#endif //SCHICTOOLS_LOAD_HIC_FILE_H

