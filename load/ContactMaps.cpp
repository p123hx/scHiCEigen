//
// Created by Bean Juice on 16/12/2020.
//

#include <string>
#include "ContactMaps.h"
#include "load_hic_file.h"
#include <vector>
#include <set>
#include <map>
#include <tuple>
using namespace std;

scHiCs::scHiCs(vector<string> &list_of_files, string reference_genome, int resolution,int
kernel_shape,int max_distance,
               bool adjust_resolution, string chromosomes, string format,
               int keep_n_strata, bool store_full_map, string operations,
               int header, int customized_format, double map_filter, bool gzip,
               bool parallelize, int n_processes, bool sparse) {
    {
        this->resolution = resolution;
        tie(this->chromosomes,this->chromosome_lengths) = get_chromosome_lengths(reference_genome,
                                                               chromosomes,
                                                               resolution);
        this->num_of_cells = list_of_files.size();
        this->sparse = sparse;
        this->keep_n_strata = keep_n_strata;
        this->files = list_of_files;

        if (keep_n_strata) {
            for (string ch : this->chromosomes) {

                for (int i = 0; i < keep_n_strata; i++) {
                    this->strata[ch].push_back(MatrixXd::Zero(this->num_of_cells,
                                                              this->chromosome_lengths[ch] -
                                                              i));
                }

            }
        }
//full_maps is not implemented

        cout << "Loading HiC data...\n";
        if (parallelize) {
            throw "Not implemented yet";
        }
        else {
            this->idx = 0;
            for (string file : this->files) {
                cout<<"loading: "<<file<<endl;
                for (string ch : this->chromosomes) {
                    size_t index = ch.find("ch");
                    if (index != string::npos && ch.find("chr") == string::npos) {
                        ch.replace(index, 2, "chr");
                    }
                    MatrixXd mat;
                    vector<MatrixXd> strata_local;
                    tie(mat, strata_local) = load_HiC(
                            file, this->chromosome_lengths,
                            format, customized_format,
                            header, ch, resolution,
                            adjust_resolution,
                            map_filter, sparse, gzip,
                            keep_n_strata, operations = operations
                    );

                    if(keep_n_strata){
                        int strata_idx = 0;
                        for(MatrixXd stratum : strata_local){
                            this->strata[ch][strata_idx].row(idx)=stratum.transpose();
                            strata_idx++;
//                            cout<<stratum<<endl;
                        }
                    }
                }
                idx++;
            }
        }
    }
}

void scHiCs::load100(vector<string> &list_of_files, string reference_genome,
                     int resolution, int kernel_shape, int max_distance,
                     bool adjust_resolution, string chromosomes, string format,
                     int keep_n_strata, bool store_full_map, string operations,
                     int header, int customized_format, double map_filter, bool gzip,
                     bool parallelize, int n_processes, bool sparse) {
    this->files = list_of_files;
    this->num_of_cells += list_of_files.size();
    if (keep_n_strata) {
        for (string ch : this->chromosomes) {

            for (int i = 0; i < keep_n_strata; i++) {

                this->strata[ch][i].conservativeResize(this->num_of_cells,
                                                       this->chromosome_lengths[ch] -
                                                       i);
            }

        }
    }
    cout << this->num_of_cells<<": Loading HiC data...\n";
    if (parallelize) {
        throw "Not implemented yet";
    }
    else {
        for (string file : this->files) {
            cout<<"loading: "<<file<<endl;
            for (string ch : this->chromosomes) {
                size_t index = ch.find("ch");
                if (index != string::npos && ch.find("chr") == string::npos) {
                    ch.replace(index, 2, "chr");
                }
                MatrixXd mat;
                vector<MatrixXd> strata_local;
                tie(mat, strata_local) = load_HiC(
                        file, this->chromosome_lengths,
                        format, customized_format,
                        header, ch, resolution,
                        adjust_resolution,
                        map_filter, sparse, gzip,
                        keep_n_strata, operations = operations
                );

                if(keep_n_strata){
                    int strata_idx = 0;
                    for(MatrixXd stratum : strata_local){
                        this->strata[ch][strata_idx].row(idx)=stratum.transpose();
                        strata_idx++;
//                            cout<<stratum<<endl;
                    }
                }
            }
            idx++;
        }
    }

}
map<string, vector<MatrixXd>> scHiCs::get_strata() {
    return this->strata;
}