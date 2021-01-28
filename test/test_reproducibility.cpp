#include "xtensor/xio.hpp"
#include "../load/ContactMaps.h"
#include "../embedding/reproducibility.h"
#include <Eigen/Dense>

using namespace std;using namespace Eigen;

double fastHicP(vector<MatrixXd> all_strata) {

    double times
            =
            pairwise_distance(all_strata, "hicrep");

    return times;
}
int main(){
    vector<string> fileLst{"../test/data/cell_03","../test/data/cell_01","../test/data/cell_02"};

    string operation = "convolution";
    scHiCs y = scHiCs(fileLst, "mm9", 500000, 3, 4000000, true, "except Y",
                      "shortest_score",
                      10, true,
                      operation);
    vector<string> chrs{"chr1", "chr2", "chrX", "chr3", "chr4", "chr5", "chr6", "chr7",
                        "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                        "chr15", "chr16", "chr17", "chr18", "chr19"}; //since "except Y"
    double tsum = 0.0;
    for (string s:chrs) {
        cout << "\n" << s << ":\n";
        vector<MatrixXd> chr = y.get_strata()[s];
        tsum += (fastHicP(chr));
    }
    cout << "time1 + time2 fastHiCrep 100 cells: " << tsum << " in milliseconds\n";
}