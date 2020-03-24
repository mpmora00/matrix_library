//////////////////////////////////////////////////////////////////////////
////Dartmouth CS70.01 starter code provided by Bo Zhu
////Linear algebra and vector math library
//////////////////////////////////////////////////////////////////////////

#ifndef __MyGraph_h__
#define __MyGraph_h__
#include <utility>
#include <tuple>
#include <iostream>
#include "MyMatrix.h"
#include "MySparseMatrix.h"

class Graph
{
public:
    vector<double> node_values;                ////values on nodes
    vector<double> edge_values;                ////values on edges
    vector<pair<int,int> > edges;        ////edges connecting nodes

    void Add_Node(const double& value)
    {
        node_values.push_back(value);
    }

    void Add_Edge(const int i,const int j,const double value=1.)
    {
        edges.push_back(pair<int,int>(i,j));
        edge_values.push_back(value);
    }

    ////display graph in terminal
    friend ostream & operator << (ostream &out,const Graph &graph)
    {
        cout<<"graph node values: "<<graph.node_values.size()<<endl;
        for(int i=0;i<(int)graph.node_values.size();i++){
            cout<<"["<<i<<", "<<graph.node_values[i]<<"] ";
        }
        cout<<endl;

        cout<<"graph edge values: "<<graph.edge_values.size()<<endl;
        for(int i=0;i<(int)graph.edge_values.size();i++){
            cout<<"["<<i<<", "<<graph.edge_values[i]<<"] ";
        }
        cout<<endl;

        cout<<"graph edges: "<<graph.edges.size()<<endl;
        for(int i=0;i<(int)graph.edges.size();i++){
            cout<<"["<<graph.edges[i].first<<", "<<graph.edges[i].second<<"] ";
        }
        cout<<endl;

        return out;
    }

    //////////////////////////////////////////////////////////////////////////

    ////HW4 Task 0: build incidence matrix
    void Incidence_Matrix(/*result*/Matrix& inc_m)
    {
        
        // change the incidence matrix size to fit the number of nodes and edges on the graph
        inc_m.Resize((int) node_values.size(), (int) edges.size());
        
        // the number of edges define which column to place number in
        for (int i=0; i < (int) edges.size(); i++) {
            
            // add a 1 on the row value of where the edge is headed towards
            // and a -1 for the row value where the edge is coming from
            inc_m(edges[i].first, i) = - 1;
            inc_m(edges[i].second, i) = 1;
        }

    }

    ////HW4 Task 1: build adjancency matrix
    // for a directed graphs
    void Adjacency_Matrix(/*result*/Matrix& adj_m)
    {

        // change the adjencency matrix size to fit the number of nodes on the graph
        adj_m.Resize((int) node_values.size(), (int) node_values.size());
                
        // for every edge add an entry on the matrix
        for (int i=0; i < (int) edges.size(); i++) {
            // add a 1 for edges coming out of this row value, and a -1 for edges coming in
            adj_m(edges[i].first, edges[i].second) = 1;
            adj_m(edges[i].second, edges[i].first) = -1;
        }
        
    }

    ////HW4 Task 3: build the negative Laplacian matrix
    void Laplacian_Matrix(/*result*/Matrix& lap_m)
    {
        
        // change the laplacian matrix size to fit the number of nodes on the graph
        lap_m.Resize((int) node_values.size(), (int) node_values.size());
        
        // for every edge add and entry on the matrix
        for (int i=0; i < (int) edges.size(); i++) {
            lap_m(edges[i].first, edges[i].second) = 1;
            lap_m(edges[i].second, edges[i].first) = 1;
        }
        
        // for each row
        for (int i=0; i < lap_m.m; i++) {
            // count the number of neighbors this node has
            int sum = 0;
            for (int j=0; j < lap_m.n; j ++) {
                if (lap_m(i, j) == 1) {
                    sum +=1;
                }
            }
            // on the diagonal, write the -# of neighbors this node has
            lap_m(i,i) = -sum;
        }
    }

    ////HW4 Task 4: calculate the Dirichlet energy
    double Dirichlet_Energy(const vector<double>& v)
    {
        double de=(double)0;
        
        // create an incidence matrix and transpose it
        Matrix inc_m;
        Incidence_Matrix(inc_m);
        
        Matrix inc_transpose;
        inc_transpose = inc_m.Transpose();
        
        Matrix p(inc_transpose.n,1);
        for(int i=0;i<v.size();i++){
            p(i,0)=v[i];
        }
        
        Matrix q(inc_transpose.m, 1);
        q = inc_transpose*p;
        for (int i=0;i<q.m;i++){
            de += q(i,0)*q(i,0);
        }
        

        return de;
    }

    ////HW4 Task 5: smooth the node values on the graph by iteratively applying the Laplacian matrix
    void Smooth_Node_Values(const double dt,const int iter_num)
    {
        ////copy node_values to local variables
        int m=(int)node_values.size();
        
        Matrix v(m,1);
        for(int i=0;i<m;i++){
            v(i,0)=node_values[i];
        }
        Matrix v2=v;

        ////smooth v
        
        
        for (int k = 0; k < iter_num; k++) {
            // build the laplacian matrix
            Matrix lap_m;
            Laplacian_Matrix(lap_m);
            v = v2 + (lap_m*v2)*dt;

            cout<<Dirichlet_Energy(v.data)<<endl;
            v2 = v;
        }

        ////copy local variables to node_values
        for(int i=0;i<m;i++){
            node_values[i]=v(i,0);
        }
    }
    
    ///// EXTRA CREDIT: Implementing with a Sparse Matrix
    void Incidence_Matrix(/*result*/SparseMatrix& inc_mtx)
       {
           
           // change the incidence matrix size to fit the number of nodes and edges on the graph
           inc_mtx.Resize((int) node_values.size(), (int) edges.size());
           
           // create the elements that will go into the sparse matrix
           vector<tuple<int,int,double> > elements;
           
           // the number of edges define which column to place number in
           for (int i=0; i < (int) edges.size(); i++) {
               
               // add a 1 on the row value of where the edge is headed towards
               // and a -1 for the row value where the edge is coming from
               
               elements.push_back(make_tuple<int,int,double>((int) edges[i].first,(int) i, -1));
               elements.push_back(make_tuple<int,int,double>((int) edges[i].second,(int) i, 1));
            }
           
           // sort the elements since my computer doesn't allow the other sort function
           sort(elements.begin(), elements.end());
           
           inc_mtx=elements;
        }
    
    // for a directed graphs
    void Adjacency_Matrix(/*result*/SparseMatrix& adj_mtx)
    {

        // change the adjencency matrix size to fit the number of nodes on the graph
        adj_mtx.Resize((int) node_values.size(), (int) node_values.size());
        
        // create the elements that will go into the sparse matrix
                  vector<tuple<int,int,double> > elements;
                
        // for every edge add an entry on the matrix
        for (int i=0; i < (int) edges.size(); i++) {
            
            // add a 1 for edges coming out of this row value, and a -1 for edges coming in
            elements.push_back(make_tuple<int,int,double>((int) edges[i].first, (int) edges[i].second, 1));
            elements.push_back(make_tuple<int,int,double>((int) edges[i].second, (int) edges[i].first,-1));
        }
        
        // sort the elements since my computer doesn't allow the other sort function
        sort(elements.begin(), elements.end());
        
        adj_mtx=elements;
    }
    
    void Laplacian_Matrix(/*result*/SparseMatrix& lap_mtx)
    {
        // change the laplacian matrix size to fit the number of nodes on the graph
        lap_mtx.Resize((int) node_values.size(), (int) node_values.size());

        // create the elements that will go into the sparse matrix
        vector<tuple<int,int,double> > elements;

        // for every edge add and entry on the matrix
        for (int i=0; i < (int) edges.size(); i++) {
            elements.push_back(make_tuple<int,int,double>((int) edges[i].first, (int) edges[i].second, 1));
            elements.push_back(make_tuple<int,int,double>((int)edges[i].second, (int) edges[i].first, 1));
        }

        sort(elements.begin(), elements.end());

        lap_mtx=elements;

        // for each row
        for (int i=0; i < lap_mtx.m; i++) {
            // count the number of neighbors this node has
            int sum = 0;
            for (int j=0; j < lap_mtx.n; j ++) {
                if (lap_mtx(i, j) == 1) {
                    sum +=1;
                }
            }
            // on the diagonal, write the -# of neighbors this node has
            elements.push_back(make_tuple<int,int,double>((int) i, (int) i, (double) -sum));
        }

        sort(elements.begin(), elements.end());

        lap_mtx=elements;
    }
};

#endif
