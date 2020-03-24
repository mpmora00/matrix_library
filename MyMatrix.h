//////////////////////////////////////////////////////////////////////////
////Dartmouth CS70.01 starter code provided by Bo Zhu
////Linear algebra and vector math library
//////////////////////////////////////////////////////////////////////////

#ifndef __MyMatrix_h__
#define __MyMatrix_h__
#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <chrono>

using namespace std;

class Matrix
{
public:
    int m;                        ////number of rows
    int n;                        ////number of columns
    vector<double> data;        ////element values, we use double for the data type

    ////matrix constructor
    Matrix(const int _m=1,const int _n=1)
    {
        Resize(_m,_n);
    }

    void Resize(const int _m,const int _n)
    {
        m=_m;
        n=_n;
        data.resize(m*n);
        for(int i=0;i<m*n;i++){
            data[i]=0.;
        }
    }

    ////A=B
    void operator = (const Matrix& B)
    {
        Resize(B.m,B.n);

        for(int i=0;i<(int)data.size();i++){
            data[i]=B.data[i];
        }
    }

    ////A={1.,4.,2.,...}, assigning a std::vector to A. A should initialized beforehand
    void operator = (const vector<double>& input_data)
    {
        assert(input_data.size()<=data.size());

        for(int i=0;i<(int)input_data.size();i++){
            data[i]=input_data[i];
        }
    }

    ////return whether A==B
    bool operator == (const Matrix& B)
    {
        assert(m == B.m&&n == B.n);
        for (int i = 0; i < (int)data.size(); i++)
        {
            if (fabs(data[i] - B.data[i]) > 1e-6) return false;        }
        return true;
    }

    ////return -A
    Matrix operator - ()
    {
        Matrix C(m,n);
        for(int i=0;i<(int)data.size();i++){
            C.data[i]=-data[i];
        }
        return C;
    }

    ////random access of a matrix element
    double& operator() (int i,int j)
    {
        assert(i>=0 &&i<m&& j>=0 && j<n);
        return data[i*n+j];
    }

    const double& operator() (int i,int j) const
    {
        assert(i>=0&&i<m&&j>=0&&j<n);
        return data[i*n+j];
    }

    ////display matrix in terminal
    friend ostream & operator << (ostream &out,const Matrix &mtx)
    {
        for(int i=0;i<mtx.m;i++){
            for(int j=0;j<mtx.n;j++){
                out<<mtx(i,j)<<", ";
            }
            out<<std::endl;
        }
        return out;
    }

    //////////////////////////////////////////////////////////////////////////
    ////overloaded operators

    ////matrix-matrix additions
    ////return C = A + B
    ////Notice: I use A to refer to the object itself in all my comments,
    ////if you want to self-access in the C++ code, you should use (*this),
    ////e.g., return (*this); means returning the object itself
    ////the comment A+=B; means (*this)+=B; in the code
    
    Matrix operator + (const Matrix& B)
    {
        assert(m==B.m&&n==B.n);

        Matrix C(m,n);
        for(int i=0;i<(int)data.size();i++){
            C.data[i]=data[i]+B.data[i];
        }
        return C;
    }

    ////A+=B
    void operator += (const Matrix& B)
    {
        assert(m==B.m&&n==B.n);

        for(int i=0;i<(int)data.size();i++){
            data[i]+=B.data[i];
        }
    }

    //////////////////////////////////////////////////////////////////////////

    ////Task 1: Mimic the "+" and "+=" operators,
    ////implement four new operators: "-", "-=", matrix-scalar multiplications "*" and "*="
    
    ////return A-B
    Matrix operator - (const Matrix& B)
    {
        
        assert(m==B.m && n==B.n);

        Matrix C(m,n);
        for (int i=0; i<(int)data.size(); i++){
            C.data[i]=data[i]-B.data[i];
        }
        return C;
    }

    ////A=A-B
    void operator -= (const Matrix& B)
    {
        assert(m==B.m && n ==B.n);
        
        for (int i=0; i<(int)data.size(); i++){
            data[i] = data[i] - B.data[i];
        }
    }

    ////return A*s, with s as a scalar
    Matrix operator * (const double& s)
    {
        Matrix C(m,n);
        for (int i=0; i<(int)data.size(); i++){
            C.data[i] = data[i] * s;
        }
        return C;
    }

    ////A=A*s, with s as a scalar
    void operator *= (const double& s)
    {
        for (int i=0; i<(int)data.size(); i++){
            data[i] = data[i] * s;
        }
    }

    ////Task 2: matrix-matrix multiplication
    ////Hints: there are four steps:
    ////1, check compatibility by an assert;
    ////2, allocate a matrix C with proper size;
    ////3, calculate each element in C by a (left)row-(right)column multiplication
    //// when accessing an element (i,j) in the object itself, use (*this)(i,j)
    ////4, return c
    
    //ORIGINAL
//    Matrix operator * (const Matrix& B)
//    {
//        assert(n==B.m);
//
//        Matrix C(m,B.n);
//        for (int i=0; i<(int)m; i++){
//            for (int j=0; j<(int)B.n; j++){
//              for (int k=0; k<(int)n; k++)
//                    C(i,j) += (*this)(i,k) * B(k,j);
//            }
//        }
//        return C;
//    }
    
    // EDITED VERSION
    Matrix operator * (const Matrix& B)
    {
        assert(n==B.m);

        Matrix C(m,B.n);
        for (int i=0; i<(int)m; i++){
            for (int k=0; k<(int)n; k++)
                for (int j=0; j<(int)B.n; j++){
                    C(i,j) += (*this)(i,k) * B(k,j);
            }
        }
        return C;
    }
    
    ////Task 3: identity, transpose(), block

    ////return an identity matrix
    Matrix Identity() {
        assert(m==n); // make sure it's a square
        Matrix C(m,n);
        for (int i=0; i<(int)m; i++){
            for (int j=0; j<(int)n; j++){
                if (i==j) {
                    C(i,j)=1;
                }
                else {
                    C(i,j)=0;
                }
            }
        }
        return C;
    }

    ////return A^T
    Matrix Transpose() {
        Matrix C(n,m);
        for (int i=0; i<(int)m; i++){
            for (int j=0; j<(int)n; j++){
                C(j,i) = (*this)(i,j);
            }
        }
        return C;
    }

    ////return a submatrix block A_ijab,
    ////with i,j as the starting element and a,b as the block size
    Matrix Block(const int i,const int j,const int a,const int b)
    {
        Matrix C(a,b);
        for (int c=0; c<(int)a; c++){
            for (int d=0; d<(int)b; d++){
                C(c,d) = (*this)(c+i,d+j);
            }
        }
        return C;
    }

    ////Task 4: implement a function or a set of functions that were not specified in class
    
    // Find the inverse of a 2x2 matrix
    Matrix Inverse() {
        assert(m==2 && n==2);
        Matrix C(m,n);
       
        // check if the determinant of the matrix is 0
        double determinant = (data[0]*data[3])-(data[1]*data[2]);
        assert(determinant!=0);
        
        // if it's not zero, then the matrix has an inverse
        for (int i=0; i<(int)m; i++){
            for (int j=0; j<(int)n; j++){
                if (i!=j) { // if it's not in the diagonal
                    C(i,j)= -(*this)(i,j);  // simply change it's sign
                }
                else {  // if it's part of the diagonal
                    // flip the numbers
                    if(i+1>=m || j+1>=n) {
                        C(i,j)=(*this)(n-j-1,m-i-1);
                    }
                    else {
                        C(i,j)=(*this)(i+1,j+1);
                    }
                }
            }
        }
        // return the Matrix x 1/det
        return C * (double)(1/determinant);
    }
    
    //////////////////////////////////////////////////////////////////////////
    ////Assignment 2: Solve Linear Equations using Gaussian Elimination
    ////return the solution of linear equations Ax=b
    
    Matrix Solve(const Matrix& b)
    {
        // Solving will only work if:
                // A is a square matrix
                // b has the same amount of rows as A
        assert (m==n && m==b.m);
        
        // create a matrix with format
        //  [ a1,1  a1,2  a1,3  | b1 ]
        //  [ a2,1  a2,2  a2,3  | b2 ]
        //  [ a3,1  a3,2  a3,3  | b3 ]
    
        Matrix C(m, n+1);
        for (int i=0; i < m; i++) {
            for (int j=0; j < n; j++) {
                C(i,j) = (*this)(i,j);
            }
            // b has 1 column so it's column[0]
            C(i, n) = b(i, 0);
        }
        
    ////Forward Elimination: transform A into an upper triangular matrix
        
    //For each row:
    //1. Search for the maximum element in current column
    //2. Swap maximum row with current row
    //3. Make all of the element 0 below current row in current column
        
        // Only go until n (ignore b values)
        for (int i=0; i < n; i++) {
            int maxRow = i;             // current max row
            
            // search for the row with the maximum element in this column
            // (starting one idx after the row we are currently in)
            for (int k = i+1; k < n; k++) {
                if (abs(C(k,i)) > abs(C(maxRow,i))) {  // if this element is bigger than the maximum element we had found, save it
                    maxRow = k;
                }
            }
            
            // swap max row with current column
            // until n+1 (because B is counted in the C matrix, there is an extra column)
            if (maxRow != i) {      // if the max row is not the current row
                for (int k = i; k < n+1; k++) {
                    swap(C(maxRow, k),C(i, k));
                }
            }
            
            // Make all of the element 0 below current row in current column
            // for every row
            for (int k = i+1; k < m; k++) {
                
                // number needed to cancel the current column
                double eliminVar = - C(k, i)/C(i,i);
                
                // for every column
                // until n+1 (because B is counted in the C matrix)
                for (int j = i; j < n+1; j++) {
                    C(k,j) += eliminVar*C(i,j);
                }
            }
        }
        
    ////Backward Substitution: solve unknowns in a reverse order
        
        // create the resulting matrix (row size == originsl row size, but only 1 column)
        Matrix result(m, 1);
        
        // starting from the last row and moving upwards
        for (int i = n - 1; i>=0; i--) {
            
            // find what the variable in last (unchecked) row equals
            // since there is only 1 column, result's column number will always be 0
            result(i, 0) = C(i,n)/C(i,i);
            
            // update the value of the last found variable on all rows above
            for (int k = i - 1; k>=0; k--) {
                C(k,n) -= C(k,i)*result(i,0);
            }
        }
        return result;
    }
    
    ///// EXTRA CREDIT: LU Decompsition
    void LUDecomp(Matrix& L, Matrix& U)
    {
        // start the upper triangle matrix
        for (int i=0; i < m; i++) {
            for (int j=0; j < n; j++) {
                U(i,j) = (*this)(i,j);
            }
        }
    
        // Create an identity matrix
        L = (*this).Identity();

        for (int i=0; i < n; i++) {
            // no need for row swapping step
            
            // Make all of the element 0 below current row in current column
            // for every row
            for (int k = i+1; k < m; k++) {
                
                // number needed to cancel the current column
                double eliminVar = - U(k, i)/U(i,i);
                L(k, i) = - eliminVar;
                
                // for every column
                for (int j = i; j < n; j++) {
                    U(k,j) += eliminVar*U(i,j);
                }
            }
        }
    }
    
    // COMBO ASSIGNMENT HOMEWORK
    // least square linear solver Ax = b for an over-determined problem
    // solving the normal equation AT Ax = AT b
    void LeastSquares(const Matrix& b, Matrix& x) {
    
        assert (m==b.m);
        
        // build the A transpoise matrix
        Matrix A_transpose;
        A_transpose = (*this).Transpose();
        
        // Use the gaussian elimination method we had written before to solve for AT*A x = AT*b
        x.Resize(n, b.n);
        x = (A_transpose*(*this)).Solve(A_transpose*b);
    }
    
    // Calculate associated residual norm squared |Ax − b|^2
    double residualNorm(const Matrix& x, const Matrix& b) {
        
        assert(m==b.m && x.n==b.n);
        
        Matrix residual(m, 1);
        residual = ((*this)*x)-b;
        
        double residual_norm = 0;
        for (int i=0; i<residual.m; i++) {
            residual_norm+= residual(i, 0)*residual(i, 0);
        }
        
        return residual_norm;
    }

    // create a matrix with random variables in it from (0-max)
    // max variable is not inclusive
    void randomMatrix(int max) {
        srand(5);
        for (int i=0; i<m; i++) {
            for (int j=0; j<n; j++) {
                (*this)(i,j) = ((double) rand()/ (double) RAND_MAX) - max;
            }
        }
    }
    
    //Solving the Lewis Richardson
    void richardsonAlgorithm(const Matrix& b, Matrix& result, const int maxInterations)
    {
        assert(m == b.m);
    
        Matrix x(n, 1);
    
        // make x1 = 0;
        for (int i = 0; i < n; i++ ) {
            x(i,0) = 0;}
          
        // build the necessary matrices for the calculations
        Matrix AT; AT = (*this).Transpose();
        Matrix I(x.m, x.m); I = I.Identity();

        
        Matrix ATA(AT.m, n); ATA = AT*(*this);
        Matrix AT_b(AT.m, b.n); AT_b= AT*b;
        Matrix xhat; (*this).LeastSquares(b, xhat);
        
        // find the sum of squares of all the elements of the matrix
        double u = 0;
        for (int i=0; i<(*this).m;i++) {
            for (int j=0; j<(*this).n;j++) {
                u += (*this)(i,j)*(*this)(i,j);
            }
        }
        u = 1/u;
        
        // open a new file
        string file_name="OutputValues.txt";
        ofstream myfile(file_name);
        // check if the file is open
        if (!myfile) {
            cout<<"Cannot open file "<<file_name<<endl;
            return;
        }
        
        // fill in x and following ones
        for (int k = 0; k < maxInterations; k++) {
            
            // the actual equation
            x = ((I-(ATA*u))*x) + (AT_b*(u));
        
            // calculate ||x-xhat||
            double res = 0;
            for (int i=0; i<x.m; i++){
                for (int j=0; j<x.n; j++) {
                    res += abs(x(i,j) - xhat(i, j));
                }
            }
            // add values into the text file
            myfile<<res<<endl;
        }
        
        // if we went over the maximum iterations placed, then return that vertex
        result.Resize(x.m, x.n);
        result = x;
    }
    
    ////Backward Substitution: solve unknowns in a reverse order
    Matrix BackwardSubstitution(Matrix& b) {
    
        assert (m==n && m==b.m);
                   
        // create the resulting matrix (row size == originsl row size, but only 1 column)
        Matrix result(m, 1);
               
        // starting from the last row and moving upwards
        for (int i = n - 1; i>=0; i--) {
                   
            // find what the variable in last (unchecked) row equals
            // since there is only 1 column, result's column number will always be 0
            result(i, 0) = b(i,0)/(*this)(i,i);
                   
            // update the value of the last found variable on all rows above
            for (int k = i - 1; k>=0; k--) {
                b(k,0) -= (*this)(k,i)*result(i,0);
            }
        }
        return result;
    }
    
    // COMBO ASSIGNMENT 2 HOMEWORK
    //QR Factorization Algorithm
    void Gram_Schmidt(Matrix& Q, Matrix& R) {
        
        // breaking A into separate column vectors
        vector<Matrix> a_column;
        for (int k=0; k<n; k++) {
            a_column.push_back((*this).Block(0, k, m, 1));
        }
        
        R.Resize(n,n);
        // breaking A into separate column vectors
        vector<Matrix> Q_columns;
        for (int i=0; i<n; i++) {
            Q_columns.push_back(a_column[i]); // Q1 = a1
            for (int j=0; j<i; j++) {
                Matrix tranpose_a;
                tranpose_a = a_column[i].Transpose();
                
                // orthogonalize
                Q_columns[i] -= Q_columns[j]*(tranpose_a*Q_columns[j])(0,0);
                R(j,i) = (tranpose_a*Q_columns[j])(0,0);
            }
            Matrix zero(Q_columns[i].m);
            
            // find magniture of Q[i]
            double magnitude = 0;
            for (int size=0; size<Q_columns[i].m; size++) {
                magnitude +=pow(Q_columns[i](size,0), 2);
            }
            
            // make sure bectore is linearly independent
            if (magnitude == 0) {cout<<"Not linearly independent"<<endl;}
            assert(magnitude!=0);
            
            R(i,i) = sqrt(magnitude);
            Q_columns[i] = Q_columns[i]*(1/sqrt(magnitude));// find qhat
            
            // make the Q Matrix
            for (int i=0; i<Q_columns.size(); i++) {
                Matrix current = Q_columns[i];
                for (int j=0; j<current.m; j++)
                    Q(j, i) = current(j, 0);
            }
        }
    }
    
    Matrix QRFactorization(Matrix& Q, Matrix& R, Matrix& b) {
        // Solve
        
        Matrix y(Q.n, Q.m);
        y = Q.Transpose()*b;
        return R.BackwardSubstitution(y);
    }
    
};

#endif
