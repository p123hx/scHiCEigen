//
// Created by Bean Juice on 16/12/2020.
//

#ifndef SCHICTOOLS_CONTACTMAPS_H
#define SCHICTOOLS_CONTACTMAPS_H

#include <string>
#include "load_hic_file.h"
#include <vector>
#include<set>
#include <map>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;
class scHiCs {
public:
    scHiCs(vector<string> list_of_files, string reference_genome, int resolution, int
    kernel_shape,int max_distance, bool
           adjust_resolution = true,  string chromosomes = "all", string
           format = "customized", int keep_n_strata = 10, bool store_full_map = false,
           string operations="", int header = 0, int customized_format =
    0,
           double map_filter = 0., bool gzip = false,
           bool parallelize = false, int n_processes = 0,bool sparse = false

    );


    map<string, vector<MatrixXd>> get_strata();
private:
    vector<string> files, operations;
    string reference_genome, format;
    set<string> chromosomes;
    map<string, vector<MatrixXd>> strata;
    map<string, int> chromosome_lengths;
    int
            resolution, keep_n_strata, header, customized_format, n_processes,
            num_of_cells;
    bool
            adjust_resolution, sparse, store_full_map, gzip, parallelize;

};

#endif //SCHICTOOLS_CONTACTMAPS_H
