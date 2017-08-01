#include "helper.h"
#include "descriptor.h"
#include <iostream>


Descriptor::Descriptor(){
    has_three_body = false;
  }

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

  // set t
  if (strcmp(name, "g4") == 0 || strcmp(name, "g5") ==0 ) {
    has_three_body = true;
  }
}

int Descriptor::get_num_descriptors() {
  int N = 0;
  for (size_t i=0; i<num_param_sets.size(); i++) {
    N += num_param_sets.at(i);
  }
  return N;
}


//*****************************************************************************
// Symmetry functions: Jorg Behler, J. Chem. Phys. 134, 074106, 2011.
//*****************************************************************************

void Descriptor::sym_g1(double r, double rcut, double &phi) {
  phi = cutoff(r, rcut);
}

void Descriptor::sym_d_g1(double r, double rcut, double &phi, double &dphi) {
  phi = cutoff(r, rcut);
  dphi = d_cutoff(r, rcut);
}

void Descriptor::sym_g2(double eta, double Rs, double r, double rcut, double &phi) {
  phi = exp(-eta*(r-Rs)*(r-Rs)) * cutoff(r, rcut);
}

void Descriptor::sym_d_g2(double eta, double Rs, double r, double rcut,
    double &phi, double &dphi)
{
  double eterm = exp(-eta*(r-Rs)*(r-Rs));
  double determ = -2*eta*(r-Rs)*eterm;
  double fc = cutoff(r, rcut);
  double dfc = d_cutoff(r, rcut);

  phi = eterm*fc;
  dphi = determ*fc + eterm*dfc;
}

void Descriptor::sym_g3(double kappa, double r, double rcut, double &phi) {
	phi = cos(kappa*r) * cutoff(r, rcut);
}

void Descriptor::sym_d_g3(double kappa, double r, double rcut, double &phi,
    double &dphi)
{
  double costerm = cos(kappa*r);
  double dcosterm = -kappa*sin(kappa*r);
  double fc = cutoff(r, rcut);
  double dfc = d_cutoff(r, rcut);

	phi = costerm*fc;
	dphi = dcosterm*fc + costerm*dfc;
}

void Descriptor::sym_g4(double zeta, double lambda, double eta,
    const double* r, const double* rcut, double &phi)
{
  double rij = r[0];
  double rik = r[1];
  double rjk = r[2];
  double rcutij = rcut[0];
  double rcutik = rcut[1];
  double rcutjk = rcut[2];
  double rijsq = rij*rij;
  double riksq = rik*rik;
  double rjksq = rjk*rjk;

  // i is the apex atom
  double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
  double costerm = pow(1+lambda*cos_ijk, zeta);
  double eterm = exp(-eta*(rijsq + riksq + rjksq));

  phi = pow(2, 1-zeta) * costerm * eterm * cutoff(rij, rcutij)
    *cutoff(rik, rcutik) * cutoff(rjk, rcutjk);
}

void Descriptor::sym_d_g4(double zeta, double lambda, double eta,
    const double* r, const double* rcut, double &phi, double* const dphi)
{
  double rij = r[0];
  double rik = r[1];
  double rjk = r[2];
  double rcutij = rcut[0];
  double rcutik = rcut[1];
  double rcutjk = rcut[2];
  double rijsq = rij*rij;
  double riksq = rik*rik;
  double rjksq = rjk*rjk;


  // cosine term, i is the apex atom
  double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
  double costerm = pow(1+lambda*cos_ijk, zeta);
  double dcos_dij = (rijsq - riksq + rjksq)/(2*rijsq*rik);
  double dcos_dik = (riksq - rijsq + rjksq)/(2*rij*riksq);
  double dcos_djk = -rjk/(rij*rik);
  double dcosterm_dcos = zeta * pow(1+lambda*cos_ijk, zeta-1) * lambda;
  double dcosterm_dij = dcosterm_dcos * dcos_dij;
  double dcosterm_dik = dcosterm_dcos * dcos_dik;
  double dcosterm_djk = dcosterm_dcos * dcos_djk;

  // exponential term
  double eterm = exp(-eta*(rijsq + riksq + rjksq));
  double determ_dij = -2*eterm*eta*rij;
  double determ_dik = -2*eterm*eta*rik;
  double determ_djk = -2*eterm*eta*rjk;

  // power 2 term
  double p2 = pow(2, 1-zeta);

  // cutoff
  double fcij = cutoff(rij, rcutij);
  double fcik = cutoff(rik, rcutik);
  double fcjk = cutoff(rjk, rcutjk);
  double fcprod = fcij*fcik*fcjk;
  double dfcprod_dij = d_cutoff(rij, rcutij)*fcik*fcjk;
  double dfcprod_dik = d_cutoff(rik, rcutik)*fcij*fcjk;
  double dfcprod_djk = d_cutoff(rjk, rcutjk)*fcij*fcik;

  // phi
  phi =  p2 * costerm * eterm * fcprod;
  // dphi_dij
  dphi[0] = p2 * (dcosterm_dij*eterm*fcprod + costerm*determ_dij*fcprod
      + costerm*eterm*dfcprod_dij);
  // dphi_dik
  dphi[1] = p2 * (dcosterm_dik*eterm*fcprod + costerm*determ_dik*fcprod
      + costerm*eterm*dfcprod_dik);
  // dphi_djk
  dphi[2] = p2 * (dcosterm_djk*eterm*fcprod + costerm*determ_djk*fcprod
      + costerm*eterm*dfcprod_djk);
}

void Descriptor::sym_g5(double zeta, double lambda, double eta,
    const double* r, const double* rcut, double &phi)
{
  double rij = r[0];
  double rik = r[1];
  double rjk = r[2];
  double rcutij = rcut[0];
  double rcutik = rcut[1];
  double rijsq = rij*rij;
  double riksq = rik*rik;
  double rjksq = rjk*rjk;

  // i is the apex atom
  double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
  double costerm = pow(1+lambda*cos_ijk, zeta);
  double eterm = exp(-eta*(rijsq + riksq));

  phi = pow(2, 1-zeta)*costerm*eterm*cutoff(rij, rcutij)*cutoff(rik, rcutik);
}

void Descriptor::sym_d_g5(double zeta, double lambda, double eta,
    const double* r, const double* rcut, double &phi, double* const dphi)
{
  double rij = r[0];
  double rik = r[1];
  double rjk = r[2];
  double rcutij = rcut[0];
  double rcutik = rcut[1];
  double rijsq = rij*rij;
  double riksq = rik*rik;
  double rjksq = rjk*rjk;

  // cosine term, i is the apex atom
  double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
  double costerm = pow(1+lambda*cos_ijk, zeta);
  double dcos_dij = (rijsq - riksq + rjksq)/(2*rijsq*rik);
  double dcos_dik = (riksq - rijsq + rjksq)/(2*rij*riksq);
  double dcos_djk = -rjk/(rij*rik);
  double dcosterm_dcos = zeta * pow(1+lambda*cos_ijk, zeta-1) * lambda;
  double dcosterm_dij = dcosterm_dcos * dcos_dij;
  double dcosterm_dik = dcosterm_dcos * dcos_dik;
  double dcosterm_djk = dcosterm_dcos * dcos_djk;

  // exponential term
  double eterm = exp(-eta*(rijsq + riksq));
  double determ_dij = -2*eterm*eta*rij;
  double determ_dik = -2*eterm*eta*rik;

  // power 2 term
  double p2 = pow(2, 1-zeta);

  // cutoff
  double fcij = cutoff(rij, rcutij);
  double fcik = cutoff(rik, rcutik);
  double fcprod = fcij*fcik;
  double dfcprod_dij = d_cutoff(rij, rcutij)*fcik;
  double dfcprod_dik = d_cutoff(rik, rcutik)*fcij;

  // phi
  phi =  p2 * costerm * eterm * fcprod;
  // dphi_dij
  dphi[0] = p2 * (dcosterm_dij*eterm*fcprod + costerm*determ_dij*fcprod
      + costerm*eterm*dfcprod_dij);
  // dphi_dik
  dphi[1] = p2 * (dcosterm_dik*eterm*fcprod + costerm*determ_dik*fcprod
      + costerm*eterm*dfcprod_dik);
  // dphi_djk
  dphi[2] = p2*dcosterm_djk*eterm*fcprod;
}



