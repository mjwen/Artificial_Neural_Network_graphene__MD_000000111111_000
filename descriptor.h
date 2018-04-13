#ifndef DESCRIPTOR_H_
#define DESCRIPTOR_H_

#include <cmath>
#include <cstring>
#include <vector>
#include <iostream>
#include "helper.h"

#define MY_PI 3.1415926535897932

// Symmetry functions taken from:

typedef double (*CutoffFunction)(double r, double rcut);
typedef double (*dCutoffFunction)(double r, double rcut);


class Descriptor
{
  public:

		std::vector<char*> name;    // name of each descriptor
		std::vector<int> starting_index;  // starting index of each descriptor
                                      // in generalized coords
		std::vector<double**> params;     // params of each descriptor
		std::vector<int> num_param_sets;  // number of parameter sets of each descriptor
		std::vector<int> num_params;      // size of parameters of each descriptor
    bool has_three_body;

    bool center_and_normalize;        // whether to center and normalize the data
    std::vector<double> features_mean;
    std::vector<double> features_std;


    // distinct values of parameters
    std::vector<double> g4_distinct_zeta;
    std::vector<double> g4_distinct_lambda;
    std::vector<double> g4_distinct_eta;
    // index of each parameter in distinct values
    std::vector<int> g4_lookup_zeta;
    std::vector<int> g4_lookup_lambda;
    std::vector<int> g4_lookup_eta;


    Descriptor();
		~Descriptor();

		// initialization helper
		void set_cutfunc(char* name);
		void add_descriptor(char* name, double** values, int row, int col);
		void set_center_and_normalize(bool do_center_and_normalize, int size,
        double* means, double* stds);

    int get_num_descriptors();

		// symmetry functions
    void sym_g1(double r, double rcut, double &phi);
    void sym_g2(double eta, double Rs, double r, double rcut, double &phi);
    void sym_g3(double kappa, double r, double rcut, double &phi);
    void sym_g4(double zeta, double lambda, double eta,
        const double* r, const double* rcut, double &phi);
    void sym_g5(double zeta, double lambda, double eta,
        const double* r, const double* rcut, double &phi);

    void sym_d_g1(double r, double rcut, double &phi, double &dphi);
    void sym_d_g2(const double eta, const double Rs, const double r, const double rcut,
        const double fcij, const double dfcij, double &phi,
        double &dphi);
    void sym_d_g3(double kappa, double r, double rcut, double &phi, double &dphi);
    void sym_d_g4(const double zeta, const double lambda, const double eta,
        const double* const r, const double* const rcut,
        const double fcij, const double fcik, const double fcjk, const double dfcij, const double dfcik, const double dfcjk,
        double &phi, double* const dphi);

    void sym_d_g4_2(const double* const r, const double* const rcut,
       const double fcprod,  const double* dfcprod_dr,
       const double costerm, const double* dcosterm_dr,
       const double eterm, const double* determ_dr,
       double &phi, double* const dphi);

    void sym_d_g5(double zeta, double lambda, double eta,
        const double* r, const double* rcut, double &phi, double* const dphi);

    void create_g4_lookup();

    void precompute_g4(const double rijmag, const double rikmag, const double rjkmag,
        const double rijsq, const double riksq, const double rjksq,
        const int n_lambda, const int n_zeta, const int n_eta,
        double** const costerm, double*** const dcosterm_dr,
        double* const eterm, double** const determ_dr);


//TODO delete; for debug purpose
    void echo_input() {
      std::cout<<"====================================="<<std::endl;
      for (size_t i=0; i<name.size(); i++) {
        int rows = num_param_sets.at(i);
        int cols = num_params.at(i);
        std::cout<<"name: "<<name.at(i)<<", rows: "<<rows<<", cols: "<<cols<<std::endl;
        for (int m=0; m<rows; m++) {
          for (int n=0; n<cols; n++) {
            std::cout<<params.at(i)[m][n]<< " ";
          }
          std::cout<<std::endl;
        }
        std::cout<<std::endl;
      }

      // centering and normalization
      std::cout<<"centering and normalizing params"<<std::endl;
      std::cout<<"means:"<<std::endl;
      for (size_t i=0; i<features_mean.size(); i++) {
        std::cout<<features_mean.at(i)<< std::endl;
      }
      std::cout<<"stds:"<<std::endl;
      for (size_t i=0; i<features_std.size(); i++) {
        std::cout<<features_std.at(i)<< std::endl;
      }

    }


	private:
		CutoffFunction cutoff;
		dCutoffFunction d_cutoff;
};


// cutoffs
inline double cut_cos(double r, double rcut) {
	if (r < rcut)
		return 0.5 * (cos(MY_PI*r/rcut) + 1);
	else
		return 0.0;
}

inline double d_cut_cos(double r, double rcut) {
	if (r < rcut)
		return -0.5*MY_PI/rcut * sin(MY_PI*r/rcut);
	else
		return 0.0;
}


//TODO correct it
inline double cut_exp(double r, double rcut) {
	if (r < rcut)
		return 1;
	else
		return 0.0;
}

inline double d_cut_exp(double r, double rcut) {
	if (r < rcut)
		return 0.0;
	else
		return 0.0;
}



inline double pow2(double x) {
  return x*x;
}

inline double pow4(double x) {
  double x2 = x*x;
  return x2*x2;
}

inline double pow8(double x) {
  double x2 = x*x;
  double x4 = x2*x2;
  return x4*x4;
}

inline double pow16(double x) {
  double x2 = x*x;
  double x4 = x2*x2;
  double x8 = x4*x4;
  return x8*x8;
}


inline double fast_pow(double base, int n) {

  double power = 0.0;
  switch(n) {
    case 1:
      power = base;
      break;
    case 2:
      power = pow2(base);
      break;
    case 4:
      power = pow4(base);
      break;
    case 8:
      power = pow8(base);
      break;
    case 16:
      power = pow16(base);
      break;
    default:
      power = std::pow(base, n);
      std::cerr <<"Warning: KIM potential, `fast_pow` does not support n = "
        << n << ". Using `std::pow`, which may be slow."<<std::endl;
  }

  return power;
}


// check wheter a double is a whole number
inline bool check_whole(double x) {

  if (std::ceil(x) == x) {
    return true;
  }
  else {
    return false;
  }

}


void add_distinct_value(double v, std::vector<double>& v_vec, double eps=1e-10);
int find_index(double v, std::vector<double>& v_vec, double eps=1e-10);


#endif // DESCRIPTOR_H_


