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
  if (strcmp(name, "sigmoid") == 0) {
    activFunc_ = &sigmoid;
    activFuncDeriv_ = &sigmoid_derivative;
  }
  else if (strcmp(name, "tanh") == 0) {
    activFunc_ = &tanh;
    activFuncDeriv_ = &tanh_derivative;
  }
  else if (strcmp(name, "relu") == 0) {
    activFunc_ = &relu;
    activFuncDeriv_ = &relu_derivative;
  }
  else if (strcmp(name, "elu") == 0) {
    activFunc_ = &elu;
    activFuncDeriv_ = &elu_derivative;
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
  RowMatrixXd act;

  // map raw C++ data into Matrix data
  // see: https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
  Map<RowMatrixXd> activation(zeta, rows, cols);

  for (int i=0; i<Nlayers_; i++) {
    preactiv_[i] = (activation * weights_[i]).rowwise() + biases_[i];
    if (i == Nlayers_ - 1) {  // output layer (no activation function applied)
      activOutputLayer_ = preactiv_[i];
    }
    else {
      act = activFunc_(preactiv_[i]);
      // cannot assign activFunc_(...) directly to activation.
      // Changing the mapped matrix `activation' does not invoke memory reallocation
      new (&activation) Map<RowMatrixXd> (act.data(), act.rows(), act.cols());
    }
  }
}

void NeuralNetwork::backward()
{
  // our cost (energy E) is the sum of activations at output layer, and no activation
  // function is employed in the output layer
  int rows = preactiv_[Nlayers_-1].rows();
  int cols  = preactiv_[Nlayers_-1].cols();

  // error at output layer
  RowMatrixXd delta_ = RowMatrixXd::Constant(rows, cols, 1.0);

  for (int i = Nlayers_ - 2; i>=0; i--) {
    // eval() is used to prevent aliasing since delta_ is both lvalue and rvalue.
    delta_ =  ( delta_ * weights_[i+1].transpose() ).eval()
        .cwiseProduct( activFuncDeriv_(preactiv_[i]) )  ;
  }

  // derivative of cost (energy E) w.r.t to input (generalized coords)
  gradInput_ = delta_ * weights_[0].transpose();
}




//*****************************************************************************
// activation functions and derivatives
//*****************************************************************************

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

RowMatrixXd elu(RowMatrixXd const& x)
{
  double alpha = 1.0;
  // the following is invalid for large alpha, e.g. alpha=10
  return x.cwiseMax((alpha*x.array().exp() - alpha).matrix());
}

RowMatrixXd elu_derivative(RowMatrixXd const& x)
{
  //TODO Eigenic way should exist
  double alpha = 1.0;
  RowMatrixXd deriv(x.rows(), x.cols());
  for (int i=0; i<x.rows(); i++) {
    for (int j=0; j<x.cols(); j++) {
      if (x(i,j) < 0.) {
        deriv(i,j) = alpha*exp(x(i,j));
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

RowMatrixXd sigmoid_derivative(RowMatrixXd const& x)
{
  RowMatrixXd s = sigmoid(x);
  return ( s.array() * (1.0 - s.array()) ).matrix();
}



