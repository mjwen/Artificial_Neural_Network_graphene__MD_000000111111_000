#include "helper.h"
#include "descriptor.h"
#include <iostream>


Descriptor::Descriptor(){
    has_three_body = false;
  }

Descriptor::~Descriptor() {
	for (size_t i=0; i<params.size(); i++) {
		Deallocate2DArray(params.at(i));
    delete [] name.at(i);
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

  int SIZE = 8;
  char* nm;
  nm = new char[SIZE];
  strcpy(nm, name);


  int index = 0;
  for (size_t i=0; i<num_param_sets.size(); i++) {
    index += num_param_sets[i];
  }

	this->name.push_back(nm);
	this->params.push_back(params);
	num_param_sets.push_back(row);
	num_params.push_back(col);
	starting_index.push_back(index);

  // set t
  if (strcmp(name, "g4") == 0 || strcmp(name, "g5") ==0 ) {
    has_three_body = true;
  }
}

void Descriptor::set_center_and_normalize(bool do_center_and_normalize, int size,
    double* means, double* stds) {
  center_and_normalize = do_center_and_normalize;
  for (int i=0; i<size; i++) {
    features_mean.push_back(means[i]);
    features_std.push_back(stds[i]);
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

  if (rij > rcutij || rik > rcutik || rjk > rcutjk) {
    phi = 0.0;
  }
  else {
    // i is the apex atom
    double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);

    double costerm;
    double base = 1+lambda*cos_ijk;
    if (base <= 0) { // prevent numerical unstability (when lambd=-1 and cos_ijk=1)
      costerm = 0;
    }
    else {
      costerm = pow(base, zeta);
    }

    double eterm = exp(-eta*(rijsq + riksq + rjksq));

    phi = pow(2, 1-zeta) * costerm * eterm * cutoff(rij, rcutij)
      *cutoff(rik, rcutik) * cutoff(rjk, rcutjk);
  }
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

  if (rij > rcutij || rik > rcutik || rjk > rcutjk) {
    phi = 0.0;
    dphi[0] = 0.0;
    dphi[1] = 0.0;
    dphi[2] = 0.0;
  }
  else {
    // cosine term, i is the apex atom
    double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
    double dcos_dij = (rijsq - riksq + rjksq)/(2*rijsq*rik);
    double dcos_dik = (riksq - rijsq + rjksq)/(2*rij*riksq);
    double dcos_djk = -rjk/(rij*rik);

    double costerm;
    double dcosterm_dcos;
    double base = 1+lambda*cos_ijk;
    if (base <= 0) { // prevent numerical unstability (when lambd=-1 and cos_ijk=1)
      costerm = 0.0;
      dcosterm_dcos = 0.0;
    }
    else {
      double power = fast_pow(base, (int)zeta);
      double power_minus1 = power/base;
      costerm = power;
      dcosterm_dcos = zeta * power_minus1 * lambda;
    }

    double dcosterm_dij = dcosterm_dcos * dcos_dij;
    double dcosterm_dik = dcosterm_dcos * dcos_dik;
    double dcosterm_djk = dcosterm_dcos * dcos_djk;

    // exponential term
    double eterm = exp(-eta*(rijsq + riksq + rjksq));
    double determ_dij = -2*eterm*eta*rij;
    double determ_dik = -2*eterm*eta*rik;
    double determ_djk = -2*eterm*eta*rjk;

    // power 2 term
    //double p2 = std::pow(2, (int)(1-zeta));
    int tmp = 1 << (int)zeta;  // compute 2^(zeta)
    double p2 = 2. / tmp;  // compute 2^(1-zeta)

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

  if (rij > rcutij || rik > rcutik) {
    phi = 0.0;
  }
  else {
    // i is the apex atom
    double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);

    double costerm;
    double base = 1+lambda*cos_ijk;
    if (base <= 0) { // prevent numerical unstability (when lambd=-1 and cos_ijk=1)
      costerm = 0;
    }
    else {
      costerm = pow(base, zeta);
    }

    double eterm = exp(-eta*(rijsq + riksq));

    phi = pow(2, 1-zeta)*costerm*eterm*cutoff(rij, rcutij)*cutoff(rik, rcutik);
  }
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

  if (rij > rcutij || rik > rcutik) {
    phi = 0.0;
    dphi[0] = 0.0;
    dphi[1] = 0.0;
    dphi[2] = 0.0;
  }
  else {
    // cosine term, i is the apex atom
    double cos_ijk = (rijsq + riksq - rjksq)/(2*rij*rik);
    double dcos_dij = (rijsq - riksq + rjksq)/(2*rijsq*rik);
    double dcos_dik = (riksq - rijsq + rjksq)/(2*rij*riksq);
    double dcos_djk = -rjk/(rij*rik);

    double costerm;
    double dcosterm_dcos;
    double base = 1+lambda*cos_ijk;
    if (base <= 0) { // prevent numerical unstability (when lambd=-1 and cos_ijk=1)
      costerm = 0.0;
      dcosterm_dcos = 0.0;
    }
    else {
      costerm = pow(base, zeta);
      dcosterm_dcos = zeta * pow(base, zeta-1) * lambda;
    }
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
}



