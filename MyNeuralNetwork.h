//////////////////////////////////////////////////////////////////////////
////Dartmouth CS70.01 starter code provided by Bo Zhu
////Linear algebra and vector math library
//////////////////////////////////////////////////////////////////////////

#ifndef __MyNeuralNetwork_h__
#define __MyNeuralNetwork_h__
#include "MyMatrix.h"
#include <utility>
#include <fstream>
#include <algorithm>

using namespace std;

class Appetizer
{
public:
    /// constructor
    Appetizer()=default;

    double function_approximation(double (*unknown_function)(double)) {
        // Function Approximation
        Matrix train_data(1000,3);
        Matrix train_targets(1000,1);
        Matrix θ(3,1);
        
        
        for (int i=0;i<1000;i++) {
            double x=(double)(rand()%20000-10000)/(double)10000;
            double y=(*unknown_function)(x);
            // set up the A matrix
            train_data(i, 0)=1;
            train_data(i, 1)=pow(x,2);
            train_data(i, 2)=pow(x,3);
            // set up y matrix
            train_targets(i, 0)=y;
        }
        
        // find θ matrix
        Matrix Q(train_data.m, train_data.n); Matrix R(train_data.n, train_data.n);
        train_data.Gram_Schmidt(Q, R);
        θ = train_data.QRFactorization(Q, R, train_targets);
        
        // evaluate approximator
        Matrix test_data(200,3);
        Matrix test_targets(200,1);
        
        for (int i=0;i<200;i++) {
            double x=(double)(rand()%20000-10000)/(double)10000;
            double y=(*unknown_function)(x);
            // set up the A matrix
            test_data(i, 0)=1;
            test_data(i, 1)=pow(x,2);
            test_data(i, 2)=pow(x,3);
            // set up y matrix
            test_targets(i, 0)=y;

        }
        Matrix yhat(test_data.m, 1);
        yhat = test_data*θ;
        
        // calculate the mean squared error (MSE)
        double MES = 0;
        for (int i=0; i<yhat.m; i++) {
            MES += pow((test_targets(i,0)-yhat(i,0)),2);
        }
        MES = MES/yhat.m;
        return MES;
    }

    double function_approximation1(double (*unknown_function)(double)) {
        // Function Approximation
        Matrix train_data(1000,4);
        Matrix train_targets(1000,1);
        Matrix θ(4,1);
        
        
        for (int i=0;i<1000;i++) {
            double x=(double)(rand()%20000-10000)/(double)10000;
            double y=(*unknown_function)(x);
            // set up the A matrix
            train_data(i, 0)=1;
            train_data(i, 1)=pow(x,2);
            train_data(i, 2)=pow(x,3);
            train_data(i, 3)=pow(x,4);

            // set up y matrix
            train_targets(i, 0)=y;
        }
        
        // find θ matrix
        Matrix Q(train_data.m, train_data.n); Matrix R(train_data.n, train_data.n);
        
        train_data.Gram_Schmidt(Q, R);
        θ = train_data.QRFactorization(Q, R, train_targets);
        
        // evaluate approximator
        Matrix test_data(200,4);
        Matrix test_targets(200,1);
        
        for (int i=0;i<200;i++) {
            double x=(double)(rand()%20000-10000)/(double)10000;
            double y=(*unknown_function)(x);
            // set up the A matrix
            test_data(i, 0)=1;
            test_data(i, 1)=pow(x,2);
            test_data(i, 2)=pow(x,3);
            test_data(i, 3)=pow(x,4);
            // set up y matrix
            test_targets(i, 0)=y;

        }
        Matrix yhat(test_data.m, 1);

        yhat = test_data*θ;
        
        // calculate the mean squared error (MSE)
        double MES = 0;
        for (int i=0; i<yhat.m; i++) {
            MES += pow((test_targets(i,0)-yhat(i,0)),2);
        }
        MES = MES/yhat.m;
        return MES;
    }

    void openFiles(string data_path, Matrix& train_data, vector<int>& train_label, Matrix& test_data, vector<int>& test_labels) {

        int tempint;
        ifstream file(data_path+"/train_data.txt");
        if (!file.is_open()) {
            std::cout<<"file does not exist"<<std::endl; }
        for (int i=0; i< 1000; i++) {
            for (int j=0; j<28*28+1; j++) {
                if (j == 784) {
                    train_data(i,j) = 1;
                }
                else{
                    if(file >> tempint) {
                        train_data(i,j)=((double)(tempint)/255.0);
                    }
                }
            }
        }
        file.close();

        file = ifstream(data_path+"/train_labels.txt");
        if (!file.is_open()) {
            std::cout<<"file does not exist"<<std::endl;}
        for(int i=0;i<1000;i++) {
            if(file >> tempint) {
                train_label.push_back(tempint);
            }
        }
        file.close();

        file = ifstream(data_path+"/test_data.txt");
        if (!file.is_open()) {
            std::cout<<"file does not exist"<<std::endl; }
        for(int i=0;i<200;i++) {
            for(int j=0;j<28*28+1;j++) {
                if (j==784) {
                    test_data(i,j) = 1;
                }
                else {
                    if(file >> tempint) {
                        test_data(i,j)=((double)tempint/255.0);
                    }
                }
            }
        }
        file.close();

        
        file = ifstream(data_path+"/test_labels.txt");
        if (!file.is_open()) {
            std::cout<<"file does not exist"<<std::endl; }
        for(int i=0;i<200;i++) {
            if(file >> tempint) {
                test_labels.push_back(tempint);
            }
        }
        file.close();
        
    }

    Matrix training_QR(Matrix& train_data, vector<int>& train_labels) {
        
        Matrix θ(train_data.n, 1);
        Matrix Q(train_data.m, train_data.n); Matrix R(train_data.n, train_data.n);
        
        train_data.Gram_Schmidt(Q, R);
        
        Matrix T(785, 10);
        // for every #(1-9)
        
        for (int k=0; k<10; k++) {
            Matrix actual_labels((int)train_data.m, 1);

            // make a y made up of 1 and -1
            for (int m=0; m<train_labels.size(); m++) {
                if (train_labels[m]==k) {
                    actual_labels(m,0)=1;
                }
                else {
                    actual_labels(m,0)=-1;
                }
            }
            
            // find Matrix T
            θ = train_data.QRFactorization(Q, R, actual_labels);
            
            for (int i=0; i<T.m; i++) {
                T(i, k) = θ(i, 0);
            }
        }

        return T;
    }
    
    Matrix training_LeastSquares(Matrix& train_data, vector<int>& train_labels) {
        
        Matrix θ(train_data.n, 1);
                
        Matrix T(785, 10);
        // for every #(1-9)
        
        for (int k=0; k<10; k++) {
            Matrix actual_labels((int)train_data.m, 1);

            // make a y made up of 1 and -1
            for (int m=0; m<train_labels.size(); m++) {
                if (train_labels[m]==k) {
                    actual_labels(m,0)=1;
                }
                else {
                    actual_labels(m,0)=-1;
                }
            }
            
            // find Matrix T
            train_data.LeastSquares(actual_labels, θ);
            
            for (int i=0; i<T.m; i++) {
                T(i, k) = θ(i, 0);
            }
        }

        return T;
    }

    vector<int> image_Classification(Matrix& test_data, Matrix& Theta) {
        vector<int> approximations;
        Matrix yhat(test_data.m, Theta.n);
        
        yhat = test_data*Theta;
        
        // calculate which image has maximum
        for (int i=0; i<yhat.m; i++) {
            int max_j = -1;
            double max_val = yhat(i,0);
            
            for (int j=0; j<yhat.n; j++) {
                if (max_val <= yhat(i,j)) {
                    max_j = j;
                    max_val = yhat(i,j);
                }
            }
            approximations.push_back(max_j);
        }
        
        return approximations;
    }

    void calculate_Accuracy(vector<int>& test_labels, vector<int> approximations){
        double correct = 0;
        assert(test_labels.size() == approximations.size());
            
        for (int i=0; i<test_labels.size(); i++) {
            if (test_labels[i]== approximations[i]) {
                correct ++;
            }
        }
         cout<<"Accuracy: "<<correct<<" out of "<<test_labels.size()<<" correct, "<<((correct/test_labels.size())*100)<<"%"<<endl;
    }
};

////Task 1: linear layer
class LinearLayer
{
private:
    Matrix stored_input;     //// Here we should store the input matrix A for Backward
public:
    Matrix weight;
    Matrix weight_grad;     //// record the gradient of the weight

    ////linear layer constructor
    explicit LinearLayer(const int _m=1,const int _n=1):
            weight(Matrix(_m,_n)),weight_grad(Matrix(_m,_n))
    {
        /* _m is the input hidden size and _n is the output hidden size
         * "Kaiming initialization" is important for neural network to converge. The NN will not converge without it!
         */
        for (auto &item : this->weight.data) {
            item=(double)(rand()%20000-10000)/(double)10000*sqrt(6.0/(double)(this->weight.m)); //// Kaiming initialization
        }
    }

    Matrix Forward(Matrix& input)
    {
        /* input.m is batch size and input.n is the #features.
         * 1) Store theinput in stored_data for Backward.
         * 2) Return input * weight.
         */

        assert(input.n==this->weight.m);
        this->stored_input = input;
        Matrix output = input * this->weight;
        return output;
    }

    ////BE CAREFUL! THIS IS THE MOST CONFUSING FUNCTION. YOU SHOULD READ THE OVERLEAF CAREFULLY BEFORE DIVING INTO THIS!
    Matrix Backward(Matrix& output_grad)
    {
        /* output_grad.m is batch size and output_grad.n is the # output features (this->weight.n)
         * 1) Calculate the gradient of the output (the result of the Forward method) w.r.t. the **weight** and store the product of the gradient and output_grad in weight_grad
         * 2) Calculate the gradient of the output (the result of the Forward method) w.r.t. the **input** and return the product of the gradient and output_grad
         */

        assert(output_grad.n == this->weight.n);
        this->weight_grad = this->stored_input.Transpose() * output_grad;
        return output_grad * this->weight.Transpose();
    }
};

////Task 2: non-linear activation
class ReLU
{
private:
    Matrix stored_input;   //// Here we should store the input matrix A for Backward
public:

    bool dessert = false;
    
    ////ReLU layer constructor
    ReLU(bool extra) {
        dessert = extra;
    }

    Matrix Forward(const Matrix& input)
    {
        /*
         *  input_data.m is batch size and input.n is the #features.
         *  This method returns the relu result for each element.
         *  TODO: 1) Go though each element in input and perform relu=max(0,x)
         *  TODO: 2) Store the input in this->stored_data for Backward.
         */
        
        Matrix output_data(input.m, input.n);
        // store the input for backwards
        this->stored_input = input;

        // for each element in input
        for (int i=0; i<input.m; i++) {
            for (int j=0; j<input.n; j++) {
                // perform relu
                if (!dessert) {
                    output_data(i,j) = max((double)0.0, (double)input(i,j));
                }
                // perform leaky relu
                else {
                    if (input(i,j) > 0) {
                        output_data(i,j) = input(i,j);
                    }
                    else {
                        output_data(i,j) = (double)((0.001)*input(i,j));
                    }
                }
            }
        }
        // return the output data
        return output_data;
    }

    Matrix Backward(const Matrix& output_grad)
    {
        /*  grad(relu)=1 if relu(x)=x
         *  grad(relu)=0 if relu(x)=0
         *  TODO: returns the gradient of the input data
         *  ATTENTION: Do not forget to multiply the grad(relu) with the output_grad
         */
        Matrix gradient(output_grad.m, output_grad.n);


        for (int i=0; i<output_grad.m; i++) {
            for (int j=0; j<output_grad.n; j++) {
                if (!dessert) {
                    if (stored_input(i,j) > 0) {
                        gradient(i,j) = output_grad(i,j);
                    }
                    else {
                        gradient(i,j) = 0;
                    }
                }
                // perform leaky relu derivative
                else {
                    if (stored_input(i,j) > 0) {
                        gradient(i,j) = output_grad(i,j);
                    }
                    else {
                        gradient(i,j) = (double)((0.001)*output_grad(i,j));
                    }
                }
            }
        }
        return gradient;
    }
};

////Task 3: Loss function
class MSELoss
{
private:
    Matrix stored_data;
public:

    ////cross entropy loss constructor
    MSELoss()= default;

    ////return the mse loss mean(y_j-y_pred_i)^2
    double Forward(Matrix& pred,const Matrix& truth)
    {
        /*  TODO: 1) return MSE(X_pred, X_truth) = ||X_pred-X_truth||^2 / n
         *  TODO: 2) store the difference in this->stored_data for Backward.
         */

        
        stored_data.Resize(pred.m, pred.n);
        assert(pred.n==truth.n && pred.m==truth.m && pred.n==stored_data.n &&pred.m==stored_data.m);
        
        double MSE = 0;
        
        // for every element in matrix
        for (int i=0; i<pred.m; i++) {
            for (int j=0; j<pred.n; j++) {
                // store the difference between the predicted and the truth
                stored_data(i,j) = (double) (pred(i,j) - truth(i,j));
                // calculate  Σ||X_pred-X_truth||^2
                MSE += (double) pow(pred(i,j) - truth(i,j), 2);
            }
        }
        // divide Σ(X_pred-X_truth)^2 by the number of samples
        MSE = (double) (MSE/pred.m);
        
        // return the MSE
        return MSE;
    }

    ////return the gradient of the input data
    Matrix Backward()
    {
        /* TODO 1) return the gradient of the MSE loss: grad(MSE) = 2(X_pred-X_truth) / n
         */
        Matrix gradient(stored_data.m, stored_data.n);
        
        for (int i=0; i<stored_data.m; i++){
            for (int j=0; j<stored_data.n; j++) {
                gradient(i,j) = (double)stored_data(i,j)*((double)2/(double)stored_data.m);
            }
        }
        return gradient;
    }
};

////Task 4: Network architecture
class Network
{
public:
    int n_layers=0;
    vector<LinearLayer> linear_layers;
    vector<ReLU> activation_layers;
    
    bool dessert1 = false;
    bool dessert2 = false;

    ////MNISTNetwork constructor
    Network(bool extra1, bool extra2, const vector<pair<int, int>>& feature_sizes)
    {
        assert(feature_sizes.size()!=0);
        for (int i=0;i<feature_sizes.size()-1;i++) {assert(feature_sizes[i].second==feature_sizes[i+1].first);}
        
        /*  TODO: 1) Initialize the array for the linear layers with the feature size specified in the vector feature_sizes.
         *                      In each pair (in_size, out_size), the in_size is the feature size of the previous layer and the out_size is the feature size of the output (that goes to the next layer)
         *                      In the linear layer, the weight should have the shape (in_size, out_size).
         *  TODO: 2) Initialize the array for the non-linear layers, the number of which should be feature_size.size()-1.
         *
         *  For example, if feature_size={{256, 128}, {128, 64}, {64, 32}},
       *                              then there are three linear layers whose weights are with shapes (256, 128), (128, 64), (64, 32),
       *                             and there are two non-linear layers.
         *  Attention: There should be one non-linear layer between two linear layers
         *                       However, if there is only one linear layer, there should be no non-linear layer at all.
         *                          The output feature size of the linear layer i should always equal to the input feature size of the linear layer i+1.
       */

        
        this->dessert1 = extra1;
        this->dessert2 = extra2;
        
        if (!extra1) {
            for (int k=0; k<feature_sizes.size(); k++) {
                if (k!= feature_sizes.size()-1) {
                    if (!extra2) {
                        activation_layers.push_back(ReLU(false));
                    }
                    else {
                        activation_layers.push_back(ReLU(true));
                    }
                }
                
                LinearLayer layer(feature_sizes[k].first, feature_sizes[k].second);
                linear_layers.push_back(layer);
            }
        }
        
        // new network architecture (1 linear) --> (2 non-linear) --> (1 linear)...
        else {
            for (int k=0; k<feature_sizes.size(); k++) {
                if (k!= feature_sizes.size()-1) {
                    if (!extra2) {
                        activation_layers.push_back(ReLU(false));
                        activation_layers.push_back(ReLU(false));
                    }
                    else {
                        activation_layers.push_back(ReLU(true));
                        activation_layers.push_back(ReLU(true));
                    }
                }
                
                LinearLayer layer(feature_sizes[k].first, feature_sizes[k].second);
                linear_layers.push_back(layer);
            }
        }

    }

    Matrix Forward(Matrix& input)
    {
        /* Propagate the input from the first layer to the last layer (before the loss function) by going through the forward functions of all the layers in linear_layers and activation_layers
         * For example, for a network with k linear layers and k-1 activation layers, the data flow is:
         * linear[0] -> activation[0] -> linear[1] ->activation[1] -> ... -> linear[k-2] -> activation[k-2] -> linear[k-1]
         * TODO: 1) propagate the input data throught the network.
         */


        Matrix curr_input = input;
        
        if (!dessert1) {
            for (int i=0; i<activation_layers.size(); i++) {
                curr_input = linear_layers[i].Forward(curr_input);
                curr_input = activation_layers[i].Forward(curr_input);
            }
        }
        // new network architecture (1 linear) --> (2 non-linear) --> (1 linear)...
        else {
            int k=0;
            for (int i=0; i<linear_layers.size()-1; i++) {
                curr_input = linear_layers[i].Forward(curr_input);

                curr_input = activation_layers[k].Forward(curr_input);
                k++;
                curr_input = activation_layers[k].Forward(curr_input);
                k++;
            }
        }
        
        curr_input = linear_layers[linear_layers.size()-1].Forward(curr_input);

        return curr_input;

    }

    ////return the gradient of the input data
    Matrix Backward(const Matrix& output_grad)
    {
        /* Propagate the gradient from the last layer to the first layer by going through the backward functions of all the layers in linear_layers and activation_layers
         * TODO: 1) propagate the gradient of the output (the one we got from the Forward method) back throught the network.
         * Notice: We should use the chain rule for the backward.
         * Notice: The order is opposite to the forward.
         */

        Matrix curr_input = output_grad;
        curr_input = linear_layers[linear_layers.size()-1].Backward(curr_input);
            
        if (!dessert1) {
            for (int i=(int)linear_layers.size()-2; i>=0; i--) {
                curr_input = activation_layers[i].Backward(curr_input);
                curr_input = linear_layers[i].Backward(curr_input);
           }
        }
        else {
            // new network architecture (1 linear) --> (2 non-linear) --> (1 linear)...
            int k = (int)(activation_layers.size()-1);
            for (int i=(int)linear_layers.size()-2; i>=0; i--) {
                 curr_input = activation_layers[k].Backward(curr_input);
                 k--;
                 curr_input = activation_layers[k].Backward(curr_input);
                 k--;
                 curr_input = linear_layers[i].Backward(curr_input);
            }
        }
       return curr_input;
    }
};

////Task 5: Matrix slicing
Matrix Matrix_Slice(Matrix& A, const int start, const int end)
{
    /*  We need to slice the matrix for batch stochastic gradient decent
     *  TODO: 1) Return a matrix with rows of the input A from row 'start' to row 'end-1'.
     */
    Matrix sliced;
    sliced = A.Block(start, 0, end-start, A.n);
    return sliced;
}

////Task 6: Regression
class Regressor
{
public:

    Network net;
    MSELoss loss_function=MSELoss();
    Matrix train_data;
    Matrix train_targets;
    Matrix test_data;
    Matrix test_targets;
    double learning_rate=1e-3;
    int max_epoch=200;
    int batch_size=32;

    ////Regressor constructor
    Regressor(bool extra1, bool extra2, vector<pair<int, int>> feature_sizes, double (*unknown_function)(const double)):
            net(Network(extra1, extra2, feature_sizes)),
            train_data(Matrix(1000,1)), train_targets(Matrix(1000,1)),
            test_data(Matrix(200,1)), test_targets(Matrix(200,1))
    {
        if (extra2) {
            learning_rate = 1e-4;
        }
        for (int i=0;i<1000;i++) {
            double x=(double)(rand()%20000-10000)/(double)10000;
            double y=unknown_function(x);
            train_data(i, 0)=x;
            train_targets(i, 0)=y;
        }

        for (int i=0;i<200;i++) {
            double x=(double)(rand()%20000-10000)/(double)10000;
            double y=unknown_function(x);
            test_data(i, 0)=x;
            test_targets(i, 0)=y;
        }
    }

    //// Here we train the network using gradient descent
    double Train_One_Epoch()
    {
        double loss=0;
        int n_loop=this->train_data.m/this->batch_size;
        for (int i=1;i<n_loop;i++){
            Matrix batch_data=Matrix_Slice(this->train_data, (i-1)*this->batch_size, i*this->batch_size);
            Matrix batch_targets=Matrix_Slice(this->train_targets, (i-1)*this->batch_size, i*this->batch_size);
            /*  Forward the data to the network.
             *  Forward the result to the loss function.
             *  Backward.
             *  Update the weights with weight gradients.
             *  Do not forget the learning rate!
             */
            Matrix pred = this->net.Forward(batch_data);
            loss += this->loss_function.Forward(pred, batch_targets);
            Matrix pred_grad = this->loss_function.Backward();
            net.Backward(pred_grad);   //// we do not need the gradient for train_data,but just the parameters.
            for (auto& item : this->net.linear_layers) {
                item.weight -= item.weight_grad*this->learning_rate;
            }
        }
        return loss/(double)n_loop;
    }

    double Test()
    {
        Matrix pred=this->net.Forward(this->test_data);
        double loss=this->loss_function.Forward(pred,this->test_targets);
        return loss;
    }

    void Train()
    {
        for (int i=0;i<this->max_epoch;i++) {
            double train_loss=Train_One_Epoch();
            double test_loss=Test();
            std::cout<<"Epoch: "<<(i+1)<<"/"<<this->max_epoch<<" | Train loss: "<<train_loss<<" | Test loss: "<<test_loss<<std::endl;
        }
    }
};


Matrix One_Hot_Encode(vector<int> labels, int classes=10)
{
    /*  Make the labels one-hot.
     *  For example, if there are 5 classes {0, 1, 2, 3, 4} then
     *  [0, 2, 4] -> [[1, 0, 0, 0, 0],
     *               [0, 0, 1, 0, 0],
     *               [0, 0, 0, 0, 1]]
     */
    Matrix output((int)labels.size(), classes);
    for (int i=0; i<labels.size(); i++) {
        output(i, labels[i]) = 1;
    }

    return output;
}


class Classifier
{
public:

    Network net;
    MSELoss loss_function=MSELoss();
    Matrix train_data;  //// The shape should be (m=n_samples,n=28^2)
    vector<int> train_labels;
    Matrix test_data;
    vector<int> test_labels;
    double learning_rate=1e-3;
    int max_epoch=200;
    int batch_size=32;
    

    // DESSERT POINT 1: increase matrix size to 28*28+1, to include the bias colum
    ////Classifier constructor
    Classifier(bool extra1, bool extra2, int num, const string& data_dict,const vector<pair<int, int>>& feature_sizes):
            net(Network(extra1, extra2, feature_sizes)),
            train_data(Matrix(1000,num)),test_data(Matrix(200,num))
    {
        
        if (extra2) {
            learning_rate = 1e-4;
        }
        
        int tempint;
        ifstream file(data_dict+"/train_data.txt");
        if (!file.is_open()) std::cout<<"file does not exist"<<std::endl;
        // // DESSERT POINT 1
        if (num != 28*28+1) {
            for(int i=0;i<1000*28*28;i++) if(file >> tempint) {train_data.data[i]=(double)tempint/255.0;}}
        else {
            for (int i=0; i< 1000; i++) {
                for (int j=0; j<28*28+1; j++) {
                    if (j == 784) {train_data(i,j) = 1;}
                    else {
                        if(file >> tempint) {
                            train_data(i,j)=((double)(tempint)/255.0);}}}}}
        file.close();

        file = ifstream(data_dict+"/train_labels.txt");
        if (!file.is_open()) std::cout<<"file does not exist"<<std::endl;
        for(int i=0;i<1000;i++) if(file >> tempint) train_labels.push_back(tempint);
        file.close();

        file = ifstream(data_dict+"/test_data.txt");
        if (!file.is_open()) std::cout<<"file does not exist"<<std::endl;
        // // DESSERT POINT 1:
        if (num != 28*28+1) {
            for(int i=0;i<200*28*28;i++) if(file >> tempint) {test_data.data[i]=(double)tempint/255.0;}}
        else {
            for(int i=0;i<200;i++) {
                for(int j=0;j<28*28+1;j++) {
                    if (j==784) {test_data(i,j) = 1;}
                    else{
                        if(file >> tempint) {
                            test_data(i,j)=((double)tempint/255.0);}}}}}
        file.close();

        file = ifstream(data_dict+"/test_labels.txt");
        if (!file.is_open()) std::cout<<"file does not exist"<<std::endl;
        for(int i=0;i<200;i++) if(file >> tempint) test_labels.push_back(tempint);
        file.close();
    }

    double Train_One_Epoch()
    {
        double loss=0;
        int n_loop=this->train_data.m/this->batch_size;
        for (int i=1;i<n_loop;i++){
            Matrix batch_data=Matrix_Slice(this->train_data, (i-1)*this->batch_size, i*this->batch_size);
            auto start=this->train_labels.begin()+(i-1)*this->batch_size;
            vector<int> batch_labels = vector<int>(start, start+this->batch_size);
           
            /*  Forward the data to the network.
             *  Forward the result to the loss function.
             *  Backward.
             *  Update the weights with weight gradients.
             *  Do not forget the learning rate!
             */
            Matrix pred = this->net.Forward(batch_data);
            
            Matrix truth;
            truth = One_Hot_Encode(batch_labels, 10);
            
            loss += this->loss_function.Forward(pred, truth);
            Matrix pred_grad = this->loss_function.Backward();
            net.Backward(pred_grad);
            for (auto& item : this->net.linear_layers) {
                item.weight -= item.weight_grad*this->learning_rate;
            }
        }
        return loss/(double)n_loop;
    }


    double Test()
    {
        Matrix score=this->net.Forward(this->test_data);    //// the class with max score is our predicted label
        double accuracy=0;
        for (int i=0;i<score.m;i++) {
            int max_index=0;
            for (int j=0;j<score.n;j++) { if (score(i,j)>score(i,max_index)) {max_index=j;} }
            if (max_index==test_labels[i]) {accuracy+=1;}
        }
        cout<<accuracy<<", "<<score.m<<endl;
        return accuracy/(double)score.m;

    }

    void Train()
    {
        for (int i=0;i<this->max_epoch;i++) {
            double loss=Train_One_Epoch();
            double accuracy=Test();
            std::cout<<"Epoch: "<<(i+1)<<"/"<<this->max_epoch<<" | Train loss: "<<loss<<" | Test Accuracy: "<<accuracy<<std::endl;
        }
    }
};

#endif
