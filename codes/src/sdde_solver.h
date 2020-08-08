#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <complex>
#include <random>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <assert.h>
#include <algorithm>
#include <functional>
//#include <fftw3.h>
#include <iomanip>
#include <Eigen/Dense>
#include <omp.h>

using namespace Eigen;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define SQUARE(X) ((X) * (X))
#define REAL 0
#define IMAG 1

typedef std::vector<double> StateVec;
typedef std::vector<double> dim1;
typedef std::vector<std::vector<double>> dim2;
typedef std::vector<std::vector<long unsigned int> > dim2I;
// typedef Matrix<unsigned long int, Dynamic, Dynamic> Mdim2I;

class SDDE
{
private:
    const int N;
    dim1 pt_t_intp, pt_t_ar;
    dim2 pt_x_intp, pt_x_ar;
    dim2 pt_x_Var;
    dim1 pt_w;
    Eigen::MatrixXd Cij;
    Eigen::MatrixXd Dij;
    std::vector<std::vector<int>> adjlistC;
    int len_t_intp;
    int seed;
    double PARk;
    double K_over_N;
    double PARtau;
    double maxdelay;
    long unsigned int MaxIter;
    int nstart;
    double tinitial;
    double tfinal;
    double dt_intp; // dt of interpolated time
    double dtmax;
    double dtmin;
    double dt0;
    double AbsTol;
    double RelTol;
    double PARD;
    double R_ave;

public:
    SDDE(int iN) : N(iN)
    {
    }
    virtual ~SDDE() {}

    void set_history(const std::vector<double> &);
    void set_matrices(
        const Eigen::MatrixXd &, 
        const Eigen::MatrixXd &,
        const std::vector<std::vector<int>>& );
    void set_params(
        double itinitial,
        double itfinal,
        double iPARk,
        double imaxdelay,
        double idtmax);

    inline double mean(const std::vector<double> &vec, const int id)
    {
        /*average the vector from element "id" to end of the vector */
        return accumulate(vec.begin() + id, vec.end(), 0.0) / (vec.size() - id);
    }

    void set_initial_frequencies(const std::vector<double> &);
    void integrate();
    dim1 KuramotoK1(const double, dim2I &, const int);
    dim1 KuramotoK2(const double, dim2I &, const int, const double, const dim1 &);
    dim1 KuramotoK3(const double, dim2I &, const int, const double, const dim1 &);
    dim1 KuramotoK4(const double, dim2I &, const int, const double, const dim1 &);
    double interp_x(double, long unsigned int &, int);
    double hermite_x(double , double , double , double ,
                     double , double , double );
    double interpolate(const std::vector<double> &X,
                       const std::vector<double> &Y,
                       const double xnew, int &n0,
                       const std::string kind);

    void mean_std(const dim1 &, const int, double &, double &);
    dim1 get_times(const std::string kind = "interpolated");
    dim1 linspace(double a, double b, int n);

    dim1 interpolate_order_parameter(
        const std::vector<int> &nodes,
        const double ti,
        const double tf,
        const double dt_intp,
        const std::string kind);

    double order_parameter(const dim1 &x);
    dim1 order_parameter_array(const std::vector<int>& nodes);

    Eigen::MatrixXd correlation(const dim1 &x);
    Eigen::MatrixXd get_correlation();
};

dim2 kuramoto_correlation(const dim1 &x);
