//////////////////////////////////////////////////////////////////////////
////Dartmouth CS70.01 starter code provided by Bo Zhu
////Linear algebra and vector math library
//////////////////////////////////////////////////////////////////////////

#ifndef __MySparseMatrix_h__
#define __MySparseMatrix_h__

#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <tuple>
#include <algorithm>
#include "MyMatrix.h"

using namespace std;

class SparseMatrix
{
public:
    int m;                        ////number of rows
    int n;                        ////number of columns

    ////HW3 Task 0: memory storage for a sparse matrix
    
    vector<int> row;
    vector<int> col;
    vector<double> val;
    

    ////matrix constructor
    SparseMatrix(const int _m=1,const int _n=1)
    {
        Resize(_m,_n);
    }

    void Resize(const int _m,const int _n)
    {
        m=_m;
        n=_n;
    }

    ////A=B
    void operator = (const SparseMatrix& B)
    {
        Resize(B.m,B.n);
        
        ////HW3 Task 0: copy the data from another sparse matrix
        
        row=B.row;
        col=B.col;
        val=B.val;
        
    }

    ////display matrix in terminal
    friend ostream & operator << (ostream &out,const SparseMatrix &mtx)
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
    ////HW3 Task 1: initialize the three arrays using an array of tuples
    void operator = (vector<tuple<int,int,double> >& elements)
    {
        row.clear();
        col.clear();
        val.clear();
        //// sort the elements by 1) row index i; 2) col index j
//        sort(elements.begin(),elements.end(),
//            [](const tuple<int,int,double>& a,const tuple<int,int,double>& b)->bool
//            {return get<0>(a)<= get<1>b)||(get<0>(a)==get<0>(b)&&get<1>(a)<=std::get<1>(b));});

        ////Using the information in elements to initialize row, col, val
        ////You may uncomment the following code to start or implement your own
        
        int r = -1;      // to compare r-1 value with current r value
        
        for(int p=0;p<elements.size();p++){
            int i=get<0>(elements[p]);        ////access row index i in tuple
            int j=get<1>(elements[p]);        ////access column index j in tuple
            double v=get<2>(elements[p]);    ////access value in tuple

            
            col.push_back(j);       // add column index to array
            val.push_back(v);       // add value to array
                       
            // if it's the first element or a different row index
            if (r != i) {
                // add the number of elements we have traversed to in the row array
                row.push_back(p);
                
                // change the number we are comparing i value to
                r = i;
            }
                        
        }
        
        row.push_back((int) elements.size());

    }
    
    ////HW3, Task 2: random access of a matrix element
    ////notice that we are not using a reference in this case
    double operator() (int i,int j) const
    {
        assert(i>=0&&i<m&&j>=0&&j<n);
        
        ////random access element (i,j) by using row,col,val
        int rowStart = row[i];      //where the row starts
        int nextRow = row[i+1];     //where the next row starts
        
        //Check all the column indices in row i
        for (int k = rowStart; k < nextRow; k++) {
            
            // if j is column index of a non-zero value, return it
            if (col[k] == j) {return val[k]; }
            
            // if we passed j, then it must be zero
            else if (col[k] > j) {return 0;}
        }
            
        return 0;
    }

    ////HW3, Task 3: sparse matrix-vector multiplication
    ////implement sparse prod=Ax, assuming x is a vector (with n==1)
    Matrix operator * (const Matrix& x)
    {
        assert(x.n==1);
        Matrix prod(x.m,1);
        
        // for every  row
        for (int i = 0; i < m; i++) {
            // for all of the non-zero columns in this row
            for( int j = row[i]; j < row[i+1]; j++){
                prod(i, 0) += val[j] * x(col[j], 0);
                
            }
        }

        return prod;
    }

    /// EXTRA CREDIT: Jacobi's Iteration
    
    Matrix Solve(const Matrix& b, const int tolerance, const int maxInterations)
    {
        assert(n == b.m && m == n);
    
        Matrix x(m, 1);   // x initial guess
        Matrix x1(m, 1);  // x prev
        
        // make k1 = 0;
        for (int i = 0; i < n; i++ ) {
            x(i,0) = 0;
        }
                
        // fill in k2 and following ones
        for (int k = 0; k < maxInterations; k++) {
            x1 = x;
            double e = 0;   // the error value

            for (int i = 0; i < n; i++) {
                double Aijxj = 0;   // finding the Aii*xj case
                for (int j=0; j < n; j++) {
                    if (j != i) {
                        Aijxj += (*this)(i,j) * x1(j,0);
                    }
                }
                
                // subtract x(k+1) - x(k)
                x(i,0) = (b(i,0) - Aijxj) / (*this)(i,i);
                // the error value
                e += abs(x(i,0) - x1(i,0));
            }
            // if the error value is smaller than the tolerance
            if (e < tolerance) {
                // return the vertex found
                return x;
            }
        }
        // if we went over the maximum iterations placed, then return that vertex
        return x;
    }
    
};

#endif
