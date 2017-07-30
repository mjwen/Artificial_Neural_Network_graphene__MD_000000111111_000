#include "network.h"


// Nothing to do at this moment
NeuralNetwork::NeuralNetwork(){}

NeuralNetwork::~NeuralNetwork(){}

void NeuralNetwork::set_nn_structure(int size_input, int num_layers,
    int* layer_sizes)
{
  inputSize_ = size_input;
  Nlayers_ = num_layers;
  for (int i=0; i<Nlayers_; i++) {
    layerSizes_.push_back(layer_sizes[i]);
  }

  weights_.resize(Nlayers_);
  biases_.resize(Nlayers_);
  preactiv_.resize(Nlayers_);
}

void NeuralNetwork::set_activation(char* name) {
  if (strcmp(name, "relu") == 0) {
    activFunc_ = &relu;
    activFuncDeriv_ = &relu_derivative;
  }
  else if (strcmp(name, "tanh") == 0) {
    activFunc_ = &tanh;
    activFuncDeriv_ = &tanh_derivative;
  }
}

void NeuralNetwork::add_weight_bias(double** weight, double* bias, int layer)
{
  int rows;
  int cols;

  if (layer == 0) {
    rows = inputSize_;
    cols = layerSizes_[layer];
  }
  else {
    rows = layerSizes_[layer-1];
    cols = layerSizes_[layer];
  }

  // copy data
  RowMatrixXd w(rows, cols);
  RowVectorXd b(cols);
  for (int i=0; i<rows; i++) {
    for (int j=0; j<cols; j++) {
      w(i,j) = weight[i][j];
    }
  }
  for (int j=0; j<cols; j++) {
    b(j) = bias[j];
  }

  // store in vector
  weights_[layer] = w;
  biases_[layer] = b;

}

void NeuralNetwork::forward(double * zeta, const int rows, const int cols)
{
  // map raw C++ data into Matrix data
  Map<RowMatrixXd> activation(zeta, rows, cols);

  for (int i=0; i<Nlayers_; i++) {
    preactiv_[i] = (activation * weights_[i]).rowwise() + biases_[i];
    if (i == Nlayers_ - 1) {  // output layer (no activation function applied)
      activOutputLayer_ = preactiv_[i];
    }
    else {
      activation = activFunc_(preactiv_[i]);
    }
  }

}

void NeuralNetwork::backward(double* dEdzeta, const int rows, const int cols)
{
  // map raw C++ data into Matrix data
  Map<RowMatrixXd> maped_dEdzeta(dEdzeta, rows, cols);


  // our `cost' is the sum of activations at output layer, and no activation
  // is employed in the output layer
  int extent1 = preactiv_[Nlayers_-1].rows();
  int extent2  = preactiv_[Nlayers_-1].cols();
  delta_ = RowMatrixXd::Constant(extent1, extent2, 1.0);

  for (int i = Nlayers_ - 2; i>=0; i--) {
    // eval() is used to prevent aliasing since delta_ is both lvalue and rvalue.
    delta_ = ( ( delta_ * weights_[i+1].transpose() )
        .cwiseProduct( activFuncDeriv_(preactiv_[i]) ) ).eval() ;
  }

  // derivative of cost (energy here) w.r.t to generalized coords zeta
  maped_dEdzeta = delta_ * weights_[0].transpose();
}














RowMatrixXd relu(RowMatrixXd const& x)
{
  return x.cwiseMax(0.0);
}

RowMatrixXd relu_derivative(RowMatrixXd const& x)
{
//TODO Eigenic way should exist
  RowMatrixXd deriv(x.rows(), x.cols());
  for (int i=0; i<x.rows(); i++) {
    for (int j=0; j<x.cols(); j++) {
      if (x(i,j) < 0.) {
        deriv(i,j) = 0.;
      } else {
        deriv(i,j) = 1.;
      }
    }
  }
  return deriv;
}

RowMatrixXd tanh(RowMatrixXd const& x)
{
  return (x.array().tanh()).matrix();
}

RowMatrixXd tanh_derivative(RowMatrixXd const& x)
{
  return (1.0 - x.array().tanh().square()).matrix();
}

RowMatrixXd sigmoid(RowMatrixXd const& x)
{
  return (1.0 / ( 1.0 + (-x).array().exp() )).matrix();
}

//TODO determine which of the following is fastest
RowMatrixXd sigmoid_derivative(RowMatrixXd const& x)
{
  RowMatrixXd s = sigmoid(x);
  return ( s.array() * (1.0 - s.array()) ).matrix();

  // return ( sigmoid(x).array() * (1.0 - sigmoid(x).array()) ).matrix();
  // The above seems fine, the compiler will optimize it such that sigmoid will only
  // be called once?


  //return ((1.0/(1.0 + (-x).array().exp())) * (1-1.0/(1.0+(-x).array().exp()))).matrix();
}




