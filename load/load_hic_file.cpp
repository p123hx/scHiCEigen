//
// Created by Bean Juice on 16/12/2020.
//

#include "load_hic_file.h"
#include <cmath>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include "processing_utils.h"
#include <Eigen/Dense>
#include <vector>
#include <iterator>
using namespace std;
using namespace Eigen;

int _res_change(int length, int res) {
    return ceil((double) length / (double) res);
}

struct Flg {
    vector<int> p1s, p2s;
    vector<double> vs;
};

Flg file_line_generator(const string& file,
                        vector<int> format =
                        vector<int>(), const string&
chrom = "", int header = 0, int resolution = 1,
                        bool resolution_adjust = true, double
                        mapping_filter = 0., bool gzip = false) {
    if (gzip) throw "gzip opening not implemented yet";
    ifstream fin;
//    cout << file << endl;
    fin.open(file);
    Flg ans;
    if (fin.is_open()) {
//        cout << "opening " << file << endl;
        if (header) {
            for (int tmp = 0; tmp < header; tmp++) {
                string dum;
                getline(fin, dum);
            }
        }
        string line;

        while (getline(fin, line)) {
            istringstream iss(line);
            vector<std::string> lst((istream_iterator<string>(iss)),
                                    istream_iterator<string>());
            if (lst.size() != 4 && lst.size() != 5 && lst.size() != 6) {
                throw "Wrong custom format!";
            }
            if (format[0] != 0 && format[2] != 0) {
                string c1 = lst[format[0] - 1],
                        c2 = lst[format[2] - 1];
                if ((c1 != chrom && "chr" + c1 != chrom) || (c2 != chrom && "chr" + c2 !=
                                                                            chrom))
                    continue;

            }
            if (format.size() == 6) {
                double q1 = stod(lst[format[4] - 1]), q2 = stod(lst[format[5] - 1]);
                if (q1 < mapping_filter or q2 < mapping_filter) continue;
            }

            int p1 = stoi(lst[format[1] - 1]), p2 = stoi(lst[format[3] - 1]);
            if (resolution_adjust) {
                p1 = p1 / resolution;
                p2 = p2 / resolution;
            }
            double v;
            if (format.size() == 4 || format.size() == 6) v = 1.0;
            else if (format.size() == 5) {
                v = stod(lst[format[4] - 1]);
            } else {
                throw "Wrong custom format!";
            }
            ans.p1s.push_back(p1);
            ans.p2s.push_back(p2);
            ans.vs.push_back(v);
        }
        fin.close();

    } else {
        throw "File opening fail";
    }
    return ans;

}

pair<set<string>, map<string, int>>
get_chromosome_lengths(const string &ref_str, const string &chromosomes, int res = 1) {
//path is in the same folder

    map<string, int> length;
    ifstream fin;

    string filename = "../load/reference_genome/" + ref_str;
    try {
        fin.open(filename);
        if (!(fin.is_open() && fin.good())) {
            cout << "not open" << endl;
            throw "not open";
        }
    } catch (const ifstream::failure &e) {
        try {
            fin.open("../load/reference_genome\\" + ref_str);
        } catch (const ifstream::failure &e) {
            cout << "Exception opening/reading file";
        }
    }
    string name, count;
    set<string> chroms;
    while (fin >> name >> count) {
        int res_change = _res_change(stoi(count), res);
        length.insert(pair<string, int>(name, res_change));
        chroms.insert(name);
    }
    if (chromosomes == "except X") {
        length.erase("chrX");
        chroms.erase("chrX");
    } else if (chromosomes == "except Y") {
        length.erase("chrY");
        chroms.erase("chrY");
    } else if (chromosomes == "except XY") {
        length.erase("chrX");
        length.erase("chrY");
        chroms.erase("chrX");
        chroms.erase("chrY");
    }

    return pair<set<string>, map<string, int >>(chroms, length);
//should we set up threshold for bad reference genome type?
}


//ref_map: char int or string int
pair<set<string>, map<string, int>>
get_chromosome_lengths(const map<string, int>& ref_map, const string& chromosomes,
                       const set<string>& chrom_set, int res = 1) {
    map<string, int> length;
    set<string> chroms;
    for (auto & it : ref_map) {
        length.insert(pair<string, int>(it.first, _res_change(it.second, res)));
        chroms.insert(it.first);
    }
    if (chromosomes == "except X") length.erase("chrX");
    else if (chromosomes == "except Y") length.erase("chrY");
    else if (chromosomes == "except XY") {
        length.erase("chrX");
        length.erase("chrY");
    }
    map<string, int> length2;

    if (!chrom_set.empty()) {
        for (const string& elm : chrom_set) {
            length2[elm] = length[elm];
        }
        return make_pair(chrom_set, length2);
    }
    return make_pair(chroms, length2);
}
pair<MatrixXd, vector<MatrixXd>> load_HiC(string file, map<string,
        int> &genome_length,
                                          string
                                          format ,
                                          int custom_format,
                                          int header ,
                                          string chromosome ,
                                          int resolution ,
                                          bool resolution_adjust ,
                                          double map_filter ,
                                          bool sparse ,
                                          bool gzip ,
                                          int keep_n_strata ,
                                          string
                                          operations
){
    vector<MatrixXd> strata = vector<MatrixXd>();
    int size = genome_length[chromosome];
    transform(format.begin(), format.end(), format.begin(), ::tolower);
    Flg gen;
    if (format == "shortest_score") {
        gen = file_line_generator(file, vector<int>{1, 2, 3, 4, 5}, chromosome, header,
                                  resolution,
                                  resolution_adjust, 0, gzip);
    } else {
        throw "Not implemented yet";
    }
    MatrixXd mat = MatrixXd::Zero(size, size);
    int gen_size = gen.vs.size();
    if (gen.vs.empty()) throw "empty vs";

    int p1_m = *min_element(gen.p1s.begin(), gen.p1s.end()), p2_m = *min_element(gen.p2s
                                                                                         .begin(),
                                                                                 gen.p2s.end());

    for (int count = 0; count < gen_size; count++) {
        int p1 = gen.p1s[count] - p1_m;
        int p2 = gen.p2s[count] - p2_m;
        double val = gen.vs[count];
        mat(p1, p2) += val;
        if (p1 != p2) mat(p2, p1) += val;
    }

    MatrixXd mat2 = matrix_operation(mat, operations);

    if (keep_n_strata) {
        int matSize = mat2.rows();
        for (int i = 0; i < keep_n_strata; i++) {
            strata.push_back(mat2.block(i,0,matSize-i,matSize-i).diagonal());
        }
    } else {
        throw "Not implemented yet";
    }
    if (sparse) throw "Not implemented yet";
    return make_pair(mat2, strata);

}


