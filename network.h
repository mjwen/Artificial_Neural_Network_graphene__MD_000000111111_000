#ifndef NETWORK_H_
#define NETWORK_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include "helper.h"

using namespace Eigen;

typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;
typedef void(*ActivationFunction) (double x, double &y);


class NeuralNetwork
{
  public:
    NeuralNetwork();
    ~NeuralNetwork();

    void set_nn_structure(int size_input, int num_layers, int* num_perceptrons);
    void set_activation(char* name);
    void add_weight_bias(double** weight, double* bias, int layer);
    double forward(const double ** zeta, double* atom_energy);
    void backward(const double *** dzetadr);



//TODO maybe delete,  for debug purpose delete
    void echo_input() {
      std::cout<<"==================================="<<std::endl;
      std::cout<<"Input data for class NeuralNetwork"<<std::endl;
      std::cout<<"Ninput_: "<<Ninput_<<std::endl;
      std::cout<<"Nlayers_: "<<Nlayers_<<std::endl;
      std::cout<<"Nperceptrons_: ";
      for (size_t i=0; i<Nperceptrons_.size(); i++) {
        std::cout<< Nperceptrons_.at(i) <<" ";
      }
      std::cout<<std::endl;

      std::cout<<"weights and biases:"<<std::endl;
      for (size_t i=0; i<weights_.size(); i++) {
        std::cout<<"w_"<<i<<std::endl<<weights_.at(i)<<std::endl;
        std::cout<<"b_"<<i<<std::endl<<biases_.at(i)<<std::endl;
      }



    }



  private:
    int Ninput_;         // size of input layer
    int Nlayers_;        // number of layers, including output
    std::vector<int> Nperceptrons_;  // number of perceptrons in each layer
    ActivationFunction act_;
    std::vector<RowMatrixXd> weights_;
    std::vector<RowVectorXd> biases_;
    RowMatrixXd dEdzeta_;



};


inline void tanh(double x, double &y) {
  y = (1 - exp(-2*x)) / (1 + exp(-2*x));
}


inline void relu(double x, double &y) {
//TODO
}



#endif // NETWORK_H_
