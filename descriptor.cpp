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

void Descriptor::sym_d_g2(const double eta, const double Rs, const double r, const double rcut,
    const double fc, const double dfc, double &phi, double &dphi)
{
  double eterm = exp(-eta*(r-Rs)*(r-Rs));
  double determ = -2*eta*(r-Rs)*eterm;
  //double fc = cutoff(r, rcut);
  //double dfc = d_cutoff(r, rcut);

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

void Descriptor::sym_d_g4(const double zeta, const double lambda, const double eta,
    const double* const r, const double* const rcut,
    const double fcij,  const double fcik, const double fcjk,
    const double dfcij, const double dfcik, const double dfcjk,
    double &phi, double* const dphi)
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
    //double fcij = cutoff(rij, rcutij);
    //double fcik = cutoff(rik, rcutik);
    //double fcjk = cutoff(rjk, rcutjk);
    double fcprod = fcij*fcik*fcjk;
    //double dfcprod_dij = d_cutoff(rij, rcutij)*fcik*fcjk;
    //double dfcprod_dik = d_cutoff(rik, rcutik)*fcij*fcjk;
    //double dfcprod_djk = d_cutoff(rjk, rcutjk)*fcij*fcik;
    double dfcprod_dij = dfcij*fcik*fcjk;
    double dfcprod_dik = dfcik*fcij*fcjk;
    double dfcprod_djk = dfcjk*fcij*fcik;


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






void Descriptor::sym_d_g4_2(const double* const r, const double* const rcut,
    const double fcprod,  const double* dfcprod_dr,
    const double costerm, const double* dcosterm_dr,
    const double eterm, const double* determ_dr,
    double &phi, double* const dphi)
{
  // 0->rij, 1->rik, 2->rjk
  if (r[0] > rcut[0] || r[1] > rcut[1] || r[2] > rcut[2]) {
    phi = 0.0;
    dphi[0] = 0.0;
    dphi[1] = 0.0;
    dphi[2] = 0.0;
  }
  else {
    // phi
    phi =  costerm * eterm * fcprod;
    // dphi_dij
    dphi[0] = dcosterm_dr[0]*eterm*fcprod + costerm*determ_dr[0]*fcprod
        + costerm*eterm*dfcprod_dr[0];
    // dphi_dik
    dphi[1] = dcosterm_dr[1]*eterm*fcprod + costerm*determ_dr[1]*fcprod
        + costerm*eterm*dfcprod_dr[1];
    // dphi_djk
    dphi[2] = dcosterm_dr[2]*eterm*fcprod + costerm*determ_dr[2]*fcprod
        + costerm*eterm*dfcprod_dr[2];
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



// precompute the terms associated with parameters of g4
void Descriptor::precompute_g4(
    const double rijmag, const double rikmag, const double rjkmag,
    const double rijsq, const double riksq, const double rjksq,
    const int n_lambda, const int n_zeta, const int n_eta,
    double**const costerm, double*** const dcosterm_dr,
    double*const eterm, double** const determ_dr)
{

  // cosine term, all terms associated with lambda and zeta, including 2^(1-zeta)
  // i is the apex atom
  double cos_ijk = (rijsq + riksq - rjksq)/(2*rijmag*rikmag);
  double dcos_dij = (rijsq - riksq + rjksq)/(2*rijsq*rikmag);
  double dcos_dik = (riksq - rijsq + rjksq)/(2*rijmag*riksq);
  double dcos_djk = -rjkmag/(rijmag*rikmag);


  for (int ilam=0; ilam<n_lambda; ilam++) {

    double lambda = g4_distinct_lambda[ilam];
    double base = 1 +  lambda * cos_ijk;

    for (int izeta=0; izeta<n_zeta; izeta++) {

      double zeta = g4_distinct_zeta[izeta];

      int tmp = 1 << (int)zeta;  // 2^(zeta)
      double p2 = 2. / tmp;  // 2^(1-zeta)



      if (base <= 0) { // prevent numerical unstability (when lambd=-1 and cos_ijk=1)
        costerm[ilam][izeta] = 0.0;
        dcosterm_dr[ilam][izeta][0] = 0.0;
        dcosterm_dr[ilam][izeta][1] = 0.0;
        dcosterm_dr[ilam][izeta][2] = 0.0;
      }
      else {
        double power = p2 * fast_pow(base, (int) zeta);
        double dcosterm_dcos = zeta * power/base * lambda;
        costerm[ilam][izeta] = power;
        dcosterm_dr[ilam][izeta][0] = dcosterm_dcos * dcos_dij;
        dcosterm_dr[ilam][izeta][1] = dcosterm_dcos * dcos_dik;
        dcosterm_dr[ilam][izeta][2] = dcosterm_dcos * dcos_djk;
      }
    }
  }


  // exponential term
  double rsq = rijsq + riksq + rjksq;
  for (int ieta=0; ieta<n_eta; ieta++) {
    double eta = g4_distinct_eta[ieta];
    eterm[ieta] = exp(-eta * rsq);
    double factor = -2*eterm[ieta]*eta;
    determ_dr[ieta][0] = factor*rijmag;
    determ_dr[ieta][1] = factor*rikmag;
    determ_dr[ieta][2] = factor*rjkmag;
  }

}




// precompute the distinct of parameter values of g4, and find the
// index of each parameter in the distinct values array
void Descriptor::create_g4_lookup() {

  double eps;

  // find distinct values of zeta, lambda, and eta
  eps = 1e-10;
  for (size_t p=0; p<this->name.size(); p++) {

    if (strcmp(this->name[p], "g4") == 0) {

      for(int q=0; q<this->num_param_sets[p]; q++) {
        double zeta = this->params[p][q][0];
        double lambda = this->params[p][q][1];
        double eta = this->params[p][q][2];

        // check wheter zeta is whole number. `fast_power` we implemented only supports integers
        if (check_whole(zeta) == false) {
          std::cerr<<"Error in KIM Potential: this model only supports integer `zeta` in `g4`."<<std::endl;
          exit(1);
        }

        add_distinct_value(zeta, g4_distinct_zeta, eps);
        add_distinct_value(lambda, g4_distinct_lambda, eps);
        add_distinct_value(eta, g4_distinct_eta, eps);
      }
    }
  }

  // find index of each parameter in distinct values g4_distinct_zeta, g4_distinct_lambda, and g4_distinct_eta
  eps = 1e-10;
  for (size_t p=0; p<this->name.size(); p++) {

    if (strcmp(this->name[p], "g4") == 0) {

      for(int q=0; q<this->num_param_sets[p]; q++) {
        double zeta = this->params[p][q][0];
        double lambda = this->params[p][q][1];
        double eta = this->params[p][q][2];

        int idx;
        idx = find_index(zeta, g4_distinct_zeta, eps);
        g4_lookup_zeta.push_back(idx);

        idx = find_index(lambda, g4_distinct_lambda, eps);
        g4_lookup_lambda.push_back(idx);

        idx = find_index(eta, g4_distinct_eta, eps);
        g4_lookup_eta.push_back(idx);

      }
    }
  }

  //@DEBUG
/*  print_vector<double> ("g4_distinct_zeta", g4_distinct_zeta);
  print_vector<double> ("g4_distinct_lambda", g4_distinct_lambda);
  print_vector<double> ("g4_distinct_eta", g4_distinct_eta);
  print_vector<int> ("g4_lookup_zeta", g4_lookup_zeta);
  print_vector<int> ("g4_lookup_lambda", g4_lookup_lambda);
  print_vector<int> ("g4_lookup_eta", g4_lookup_eta);
*/

}


// add v to v_vec if v is not in v_vec
// two double considered the same if their abs difference is smaller than eps
void add_distinct_value(double v, std::vector<double>& v_vec, double eps) {

  bool already_in = false;
  for (size_t i=0; i<v_vec.size(); i++) {
    if (std::abs(v - v_vec[i]) < eps) {
      already_in = true;
      break;
    }
  }
  if (already_in == false) {
    v_vec.push_back(v);
  }

}

// find the index of v in v_vec
// two double considered the same if their abs difference is smaller than eps
int find_index(double v, std::vector<double>& v_vec, double eps) {
  int idx = -1;

  for (size_t i=0; i<v_vec.size(); i++) {
    if (std::abs(v - v_vec[i]) < eps) {
      idx = i;
      break;
    }
  }
  if (idx == -1) {
    std::cerr<<"KIM model, cannot find v = " <<v<< " int v_vec."<< std::endl;
    exit(1);
  }

  return idx;
}


