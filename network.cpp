#include "network.h"


// Nothing to do at this moment
NeuralNetwork::NeuralNetwork(){}

NeuralNetwork::~NeuralNetwork(){}



void NeuralNetwork::set_nn_structure(int size_input, int num_layers,
    int* num_perceptrons)
{
  Ninput_ = size_input;
  Nlayers_ = num_layers;
  for (int i=0; i<Nlayers_; i++) {
    Nperceptrons_.push_back(num_perceptrons[i]);
  }
}

void NeuralNetwork::set_activation(char* name) {
  if (strcmp(name, "relu") == 0) {
    act_ = &relu;
  }
  else if (strcmp(name, "tanh") == 0) {
    act_ = &relu;
  }
}

void NeuralNetwork::add_weight_bias(double** weight, double* bias, int layer)
{

  int row;
  int col;

  if (layer == 0) {
    row = Ninput_;
    col = Nperceptrons_[layer];
  }
  else {
    row = Nperceptrons_[layer-1];
    col = Nperceptrons_[layer];
  }

  // copy data
  RowMatrixXd w(row, col);
  RowVectorXd b(col);
  for (int i=0; i<row; i++) {
    for (int j=0; j<col; j++) {
      w(i,j) = weight[i][j];
    }
  }
  for (int j=0; j<col; j++) {
    b(j) = bias[j];
  }

  // store in vector
  weights_.push_back(w);
  biases_.push_back(b);

}




//    double forward(const double ** zeta, double* atom_energy);
//    void backward(const double *** dzetadr);


