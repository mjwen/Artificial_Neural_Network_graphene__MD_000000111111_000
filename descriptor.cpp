#include "helper.h"
#include "descriptor.h"
#include <iostream>


// nothing to do at this moment
Descriptor::Descriptor(){}

Descriptor::~Descriptor() {
	for (size_t i=0; i<desc_params.size(); i++) {
		Deallocate2DArray(desc_params.at(i));
	}
}


void Descriptor::set_cutfunc(char* name)
{
	if (strcmp(name, "cos") == 0) {
		cutoff = &cut_cos;
		d_cutoff = &d_cut_cos;
	}
	else if (strcmp(name, "exp") == 0) {
		cutoff = &cut_exp;
		d_cutoff = &d_cut_exp;
	}

//TODO delete
std::cout<<"flag using cutfun "<<name<<std::endl;

}

void Descriptor::add_descriptor(char* name, double** values, int row, int col)
{
	double ** params = 0;
	AllocateAndInitialize2DArray(params, row, col);
	for (int i=0; i<row; i++) {
		for (int j=0; j<col; j++) {
			params[i][j] = values[i][j];
		}
	}
	desc_name.push_back(name);
	desc_params.push_back(params);
	num_param_sets.push_back(row);
	num_params.push_back(col);

//TODO delete
  std::cout<<"flag add_descriptor"<<std::endl;
  print_mat(desc_params.back(), num_param_sets.back(), num_params.back());


}


double Descriptor::sym_g2(double r, double rcut, double eta, double Rs) {
	return exp(-eta*pow((r-Rs),2)) * cutoff(r, rcut);
}

double Descriptor::sym_d_g2(double r, double rcut, double eta, double Rs) {
	return exp(-eta*pow((r-Rs),2))*d_cutoff(r, rcut) - 2*eta*(r-Rs)*sym_g2(r, rcut, eta, Rs);
}


double Descriptor::sym_g3(double r, double rcut, double kappa) {
	return cos(kappa*r)*cutoff(r, rcut);
}
double Descriptor::sym_d_g3(double r, double rcut, double kappa) {
	return cos(kappa*r)*d_cutoff(r, rcut) - kappa*sin(kappa*r)*cutoff(r, rcut);
}



int main(){}
