#ifndef NETWORK_H_
#define NETWORK_H_

#include <cmath>
#include <vector>
#include <iostream>
#include <Eigen/Core>
#include "helper.h"

using namespace Eigen;

// typedef function pointer
typedef Matrix<double, Dynamic, Dynamic, RowMajor> RowMatrixXd;

typedef  RowMatrixXd(*ActivationFunction) (RowMatrixXd const& x);
typedef  RowMatrixXd(*ActivationFunctionDerivative) (RowMatrixXd const& x);


class NeuralNetwork
{
  public:
    NeuralNetwork();
    ~NeuralNetwork();

    void set_nn_structure(int input_size, int num_layers, int* layer_sizes);
    void set_activation(char* name);
    void add_weight_bias(double** weight, double* bias, int layer);
    void forward(double * zeta, const int rows, const int cols);
    void backward(double* dEdzeta, const int rows, const int cols);

    double reduce_sum_output() {
      return activOutputLayer_.sum();
    }


//TODO maybe delete,  for debug purpose delete
    void echo_input() {
      std::cout<<"==================================="<<std::endl;
      std::cout<<"Input data for class NeuralNetwork"<<std::endl;
      std::cout<<"inputSize_: "<<inputSize_<<std::endl;
      std::cout<<"Nlayers_: "<<Nlayers_<<std::endl;
      std::cout<<"Nperceptrons_: ";
      for (size_t i=0; i<layerSizes_.size(); i++) {
        std::cout<< layerSizes_.at(i) <<" ";
      }
      std::cout<<std::endl;

      std::cout<<"weights and biases:"<<std::endl;
      for (size_t i=0; i<weights_.size(); i++) {
        std::cout<<"w_"<<i<<std::endl<<weights_.at(i)<<std::endl;
        std::cout<<"b_"<<i<<std::endl<<biases_.at(i)<<std::endl;
      }



    }



  private:
    int inputSize_;         // size of input layer
    int Nlayers_;           // number of layers, including output, excluding input
    std::vector<int> layerSizes_;  // number of perceptrons in each layer
    ActivationFunction activFunc_;
    ActivationFunctionDerivative activFuncDeriv_;
    std::vector<RowMatrixXd> weights_;
    std::vector<RowVectorXd> biases_;
    std::vector<RowMatrixXd> preactiv_;
    RowMatrixXd activOutputLayer_;
    RowVectorXd delta_;
    RowMatrixXd dEdzeta_;



};


// activation fucntion and derivatives
RowMatrixXd relu(RowMatrixXd const& x);
RowMatrixXd relu_derivative(RowMatrixXd const& x);
RowMatrixXd tanh(RowMatrixXd const& x);
RowMatrixXd tanh_derivative(RowMatrixXd const& x);
RowMatrixXd sigmoid(RowMatrixXd const& x);
RowMatrixXd sigmoid_derivative(RowMatrixXd const& x);


#endif // NETWORK_H_
