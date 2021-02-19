# scHiCTools

### Summary
A C++ version used in benchmarking https://github.com/liu-bioinfo-lab/scHiCTools.git
mainly to compare the running time of reproducibility process.

### Installation
**Required C++ Package**
- C++11
- Eigen
- MKL
- OpenMP

## Usage
To reproduce the table of Magano 1000 cells:
Linux cmd: 
Go to cmake loaction: cd scHiCEigen/cmake-build-debug/
Source mkl location to configure if needed: for example
source /opt/intel/mkl/2021.1/setvars.sh

cmake -D CMAKE_CXX_COMPILER=/opt/gcc/8.4/modbin/g++ CMAKE_BUILD_TYPE=Release ..
make
./scHiCEigen

Output will be in scHiCTest/cmake-build-debug/ files named $#ofcell$out.txt as well as print out


### Citation
Xinjun Li, Fan Feng, Wai Yan Leung and Jie Liu. "scHiCTools: a computational toolbox for analyzing single cell Hi-C data."
bioRxiv (2019): 769513.

### References
A. R. Ardakany, F. Ay, and S. Lonardi.  Selfish: Discovery of differential chromatininteractions via a self-similarity measure.bioRxiv, 2019.

N. C. Durand, J. T. Robinson, M. S. Shamim, I. Machol, J. P. Mesirov, E. S. Lander, and E. Lieberman Aiden. "Juicebox provides a visualization system for Hi-C contact maps with unlimited zoom." Cell Systems 3(1), 2016.

J. Liu, D. Lin, G. Yardimci, and W. S. Noble. Unsupervised embedding of single-cellHi-C data.Bioinformatics, 34:96–104, 2018.

T. Nagano, Y. Lubling, C. Várnai, C. Dudley, W. Leung, Y. Baran, N. M. Cohen,S.  Wingett,  P.  Fraser,  and  A.  Tanay.    Cell-cycle  dynamics  of  chromosomal organization at single-cell resolution.Nature, 547:61–67, 2017.

B. Wang, A. Pourshafeie, M. Zitnik, J. Zhu, C. D. Bustamante, S. Batzoglou, andJ.  Leskovec.   Network  enhancement  as  a  general  method  to  denoise  weighted biological networks.Nature Communications, 9(1):3108, 2018.

T. Yang, F. Zhang, G. G. Y. mcı, F. Song, R. C. Hardison, W. S. Noble, F. Yue, andQ. Li. HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient.Genome Research, 27(11):1939–1949, 2017.

G.  G.  Yardımcı,   H.  Ozadam,   M.  E.  Sauria,   O.  Ursu,   K. K.  Yan,   T.  Yang,A. Chakraborty, A. Kaul, B. R. Lajoie, F. Song, et al. Measuring the reproducibilityand quality of hi-c data.Genome Biology, 20(1):57, 2019.