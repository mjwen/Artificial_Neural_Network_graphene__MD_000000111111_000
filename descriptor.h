#ifndef DESCRIPTOR_H_
#define DESCRIPTOR_H_

#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include "helper.h"

#define MY_PI 3.1415926535897932


typedef double (*CutoffFunction)(double r, double rcut);
typedef double (*dCutoffFunction)(double r, double rcut);


class Descriptor
{
  public:
    Descriptor();
		~Descriptor();

		// initialization helper
		void set_cutfunc(char* name);
		void add_descriptor(char* name, double** values, int row, int col);

		// symmetry functions
    double sym_g2(double r, double rcut, double eta, double Rs);
    double sym_d_g2(double r, double rcut, double eta, double Rs);
    double sym_g3(double r, double rcut, double kappa);
    double sym_d_g3(double r, double rcut, double kappa);

	private:
		CutoffFunction cutoff;
		dCutoffFunction d_cutoff;
		std::vector<std::string> desc_name;  // name of each descriptor
		std::vector<double**> desc_param;  // params of each descriptor
		std::vector<int> num_desc;   // number of parameter sets of each descriptor
		std::vector<int> num_param;  // number of parameters of each descriptor
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



#endif // DESCRIPTOR_H_


