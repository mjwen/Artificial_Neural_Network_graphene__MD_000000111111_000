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
  keep_prob_.resize(Nlayers_);
  keep_prob_binary_.resize(Nlayers_);
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

void NeuralNetwork::set_keep_prob(double* keep_prob) {
  for (int i=0; i<Nlayers_; i++) {
    keep_prob_[i] = keep_prob[i];
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
  act = activation;

  for (int i=0; i<Nlayers_; i++) {

    // apply dropout
    if (static_cast<int>(keep_prob_[i]) != 1) {
      act = dropout_(act, i);  // no aliasing will occur for act
    }

    preactiv_[i] = (act * weights_[i]).rowwise() + biases_[i];

    if (i == Nlayers_ - 1) {  // output layer (no activation function applied)
      activOutputLayer_ = preactiv_[i];
    }
    else {
      act = activFunc_(preactiv_[i]);
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
  RowMatrixXd delta = RowMatrixXd::Constant(rows, cols, 1.0);

  for (int i = Nlayers_ - 2; i>=0; i--) {
    // eval() is used to prevent aliasing since delta is both lvalue and rvalue.
    delta = ( delta * weights_[i+1].transpose() ).eval()
        .cwiseProduct(activFuncDeriv_(preactiv_[i]));

    // apply dropout
    if (static_cast<int>(keep_prob_[i+1]) != 1) {
      delta = delta.cwiseProduct(keep_prob_binary_[i+1]) / keep_prob_[i+1];  // no aliasing will occur
    }
  }

  // derivative of cost (energy E) w.r.t to input (generalized coords)
  if (static_cast<int>(keep_prob_[0]) != 1) {
    gradInput_ = (delta * weights_[0].transpose()).cwiseProduct(keep_prob_binary_[0])/keep_prob_[0];
  }
  else {
    gradInput_ = delta * weights_[0].transpose();
  }

}


// dropout
RowMatrixXd NeuralNetwork::dropout_(RowMatrixXd const& x, int layer)
{

  RowMatrixXd y;
  double keep_prob = keep_prob_[layer];

  if (static_cast<int>(keep_prob) != 1) {
    // uniform [-1, 1]
    RowMatrixXd random = RowMatrixXd::Random(x.rows(), x.cols());
    // uniform [keep_prob, 1+keep_prob]
    random = (random/2.).array() + (0.5 + keep_prob);

    keep_prob_binary_[layer] = random.array().floor();

    //TODO delete  This is for debug pourpose, should be used together with `openkim-fit/tests/test_ann_force_dropout.py`
    bool debug = false;
    if (debug) {
      keep_prob_binary_[layer] = RowMatrixXd::Ones(x.rows(), x.cols());
      keep_prob_binary_[layer](0,0) = 0.;
      keep_prob_binary_[layer](0,5) = 0.;
      keep_prob_binary_[layer](1,2) = 0.;
      keep_prob_binary_[layer](1,7) = 0.;
      keep_prob_binary_[layer](3,2) = 0.;
      keep_prob_binary_[layer](3,5) = 0.;
    }

    y = (x/keep_prob).cwiseProduct(keep_prob_binary_[layer]);
  }
  else {
    y = x;
  }

  return y;
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



