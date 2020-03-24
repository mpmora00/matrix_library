#define HW_ID 7
#include <iostream>
#include "MyMatrix.h"
#if HW_ID >= 3
#include "MySparseMatrix.h"        ////uncomment this for HW3
#endif
#if HW_ID >= 4
#include "MyGraph.h"            ////uncomment this for HW4
#endif
#if HW_ID >= 7
#include "MyNeuralNetwork.h"            ////uncomment this for Final Project
#endif


void Test_HW1()
{
#if HW_ID == 1

    ////test the sample code
    Matrix m1(3,3),m2(3,3);
    
    m1={1.,2.,3.,4.,5.,6.,7.,8.,9.};
    m2={1.,0.,0.,0.,2.,0.,0.,0.,3.};
    cout<<"m1:\n"<<m1<<endl;
    cout<<"m2:\n"<<m2<<endl;

    cout<<"m1+m2:\n"<<(m1+m2)<<endl;

    Matrix m3=m1;
    cout<<"m3=m1, m3:\n"<<m3<<endl;

    //////////////////////////////////////////////////////////////////////////
    ////start to test your implementation
    ////Notice the code will not compile before you implement your corresponding operator
    ////so uncomment them one by one if you want to test each function separately
    
    ////test subtractions
    Matrix m4=-m3;
    cout<<"m4=-m3, m4:\n"<<m4<<endl;
    
    m4=m1-m3;
    cout<<"m4=m1-m3,m4:\n"<<m4<<endl;

    m4-=m1;
    cout<<"m4-=m1, m4:\n"<<m4<<endl;

    ////test matrix-scalar products
    double s=2;
    Matrix m5=m4*s;
    cout<<"m5=m4*s, m5:\n"<<m5<<endl;
    m5*=s;
    cout<<"m5*=s, m5:\n"<<m5<<endl;

    ////check matrix-matrix multiplication
    Matrix v1(3,1);    ////column vector
    v1={1.,2.,3.};
    cout<<"column vector v1:\n"<<v1<<endl;

    Matrix v2(1,3);    ////row vector
    v2={-3.,-2.,-1.};
    cout<<"row vector v2:\n"<<v2<<endl;

    Matrix v3=v1*v2;
    cout<<"v3=v1*v2, dimension: ["<<v3.m<<", "<<v3.n<<"]"<<endl;
    cout<<"v3 values:\n"<<v3<<endl;

    Matrix v4=v2*v1;
    cout<<"v4=v2*v1, dimension: ["<<v4.m<<", "<<v4.n<<"]"<<endl;
    cout<<"v4 values:\n"<<v4<<endl;

    ////test identity, transpose, and block
    Matrix m6(3,3);
    cout<<"m6:\n"<<m6.Identity()<<endl;

    Matrix m7(4,2);
    cout<<"m7.Transpose():\n"<<m7.Transpose()<<endl;

    cout<<"m2.Block(0,0,2,2):\n"<<m2.Block(0,0,2,2)<<std::endl;
#endif
}


 void Test_HW2()
{
        #if HW_ID == 2
        ////test Gaussian Elimination
        Matrix A1(3, 3);
        Matrix b1(3, 1);
        Matrix x1(3, 1);
        A1 = { 1.,1.,1.,2.,2.,5.,4.,6.,8. };
        b1 = { 1.,2.,3. };
        x1 = A1.Solve(b1);
        if ((A1*x1) == b1) cout << "Solve A1x=b1:\n " << x1 << endl;
        else cout << "Wrong Answer for A1x=b1" << endl;
        Matrix A2(3, 3);
        Matrix b2(3, 1);
        A2 = { 5.,1.,3.,4.,5.,3.,1.,5.,2. };
        b2 = { 3.,6.,-1. };
        Matrix x2(3, 1);
        x2 = A2.Solve(b2);
        if (A2*x2 == b2) cout << "Solve A2x=b2:\n " << x2 << endl;
        else cout << "Wrong Answer for A2x=b2" << endl;
        Matrix A3(5, 5);
        Matrix b3(5, 1);
        A3 = { 2.,4.,5.,3.,2.,
            4.,8.,3.,4.,3.,
            3.,3.,2.,7.,2.,
            1.,2.,2.,1.,3.,
            3.,4.,2.,5.,7. };
        b3 = { 7.,-4.,-15.,14.,16. };
        Matrix x3(5, 1);
        x3 = A3.Solve(b3);
        if (A3*x3 == b3) cout << "Solve A3x=b3:\n " << x3 << endl;
        else cout << "Wrong Answer for A3x=b3" << endl;
        Matrix A4(1000, 1000);
        Matrix b4(1000, 1);
        for (int i = 1; i < A4.m -1; i++) {
            A4(i, i) = 2.;
            A4(i, i + 1) = -1.;
            A4(i, i - 1) = -1.;
        }
        A4(0, 0) = A4(A4.m-1, A4.n-1) = 2.;
        A4(0, 1) = A4(A4.m - 1, A4.n - 2) = -1.;
        for (int i = 0; i < b4.m; i++) {
            b4(i, 0) = (double)i/(double)(b4.m*b4.m); }
        Matrix x4(1000, 1);
        x4 = A4.Solve(b4);
        if (A4*x4 == b4) {
            cout << "Solve A4x=b4:\n " << x4 << endl;}
        else {
            cout << "Wrong Answer for A4x=b4" << endl; }

    ////test LU Decomposition
    //// EXTRA CREDIT
    Matrix c1(5, 5);
    c1 = { 2.,4.,5.,3.,2.,
           5.,4.,3.,4.,1.,
           3.,3.,2.,7.,2.,
           1.,2.,2.,1.,3.,
           3.,4.,2.,5.,7. };
    Matrix Lc1(5, 5);
    Matrix Uc1(5, 5);
    c1.LUDecomp(Lc1, Uc1);
    cout << "U c1:\n" << Uc1 << endl;
    cout << "L c1:\n" << Lc1 << endl;
    cout << "A=LU proof:\n" << Lc1*Uc1 <<endl;
    Matrix c2(3, 3);
    c2 = { 1.,2.,3.,2.,3.,4.,6.,4.,1. };
    Matrix Lc2(3, 3);
    Matrix Uc2(3, 3);
    c2.LUDecomp(Lc2, Uc2);
    cout << "U c2:\n" << Uc2 << endl;
    cout << "L c2:\n" << Lc2 << endl;
    cout << "A=LU proof:\n" << Lc2*Uc2 <<endl;
    Matrix c3(3, 3);
    c3 = { 5.,1.,3.,4.,5.,3.,1.,5.,2. };
    Matrix Lc3(3, 3);
    Matrix Uc3(3, 3);
    c3.LUDecomp(Lc3, Uc3);
    cout << "U c3:\n" << Uc3 << endl;
    cout << "L c3:\n" << Lc3 << endl;
    cout << "A=LU proof:\n" << Lc3*Uc3 <<endl;
#endif
}

void Test_HW3()
{
    #if HW_ID == 3
        std::cout<<"Test sparse matrix"<<std::endl;
        SparseMatrix mtx(5,5);
        vector<tuple<int,int,double> > elements;
        elements.push_back(make_tuple<int,int,double>(0,0,7));
        elements.push_back(make_tuple<int,int,double>(0,1,5));
        elements.push_back(make_tuple<int,int,double>(1,0,1));
        elements.push_back(make_tuple<int,int,double>(1,2,3));
        elements.push_back(make_tuple<int,int,double>(2,3,5));
        elements.push_back(make_tuple<int,int,double>(2,4,4));
        elements.push_back(make_tuple<int,int,double>(3,3,1));
        elements.push_back(make_tuple<int,int,double>(4,1,7));
        elements.push_back(make_tuple<int,int,double>(4,4,3));
        mtx=elements;

        cout<<"sparse matrix:\n"<<mtx<<endl;

        Matrix v(5,1);v={1,2,3,4,5};
        Matrix prod(5,1);
        prod=mtx*v;
        cout<<"sparse matrix-vector multiplication:\n";
        cout<<prod<<endl;
    
    /// EXTRA CREDIT: Jacobi's Iteration
    SparseMatrix mtx2(3,3);
    vector<tuple<int,int,double> > elements2;
    elements2.push_back(make_tuple<int,int,double>(0,0,2));
    elements2.push_back(make_tuple<int,int,double>(0,1,-1));
    elements2.push_back(make_tuple<int,int,double>(1,0,-1));
    elements2.push_back(make_tuple<int,int,double>(1,1,4));
    elements2.push_back(make_tuple<int,int,double>(1,2,-1));
    elements2.push_back(make_tuple<int,int,double>(2,1,-1));
    elements2.push_back(make_tuple<int,int,double>(2,2,2));
    
    mtx2=elements2;
    Matrix b(3,1);
    b={0,4,4};
    cout<<"sparse matrix 2:\n"<<mtx2<<endl;
    cout<<mtx2.Solve(b, 1e-4, 1000);
#endif
}

void Test_HW4()
{
#if HW_ID ==4
    std::cout<<"Test graph matrix"<<std::endl;
    Graph g;
    g.Add_Node(0.);
    g.Add_Node(1.);
    g.Add_Node(2.);
    g.Add_Node(3.);
    g.Add_Node(4.);
    g.Add_Node(5.);

    g.Add_Edge(0,1);
    g.Add_Edge(1,2);
    g.Add_Edge(1,3);
    g.Add_Edge(2,3);
    g.Add_Edge(2,4);
    g.Add_Edge(3,4);
    g.Add_Edge(4,5);

    Matrix adj_m; g.Adjacency_Matrix(adj_m);
    Matrix inc_m; g.Incidence_Matrix(inc_m);
    Matrix lap_m; g.Laplacian_Matrix(lap_m);
    double energy=g.Dirichlet_Energy(g.node_values);

    cout<<g<<endl;
    cout<<"Adjacency matrix\n"<<adj_m<<endl;
    cout<<"Incidency matrix\n"<<inc_m<<endl;
    cout<<"Laplacian matrix\n"<<lap_m<<endl;
    cout << "Dirichlet energy before smoothing: "<<energy<<endl;

    g.Smooth_Node_Values(.1,10);
    energy=g.Dirichlet_Energy(g.node_values);
    cout<<"Dirichlet energy after smoothing: "<<energy<<endl;
    
    
    //// EXTRA CREDIT: Implement a Sparse Matrix
    SparseMatrix inc_mtx; g. Incidence_Matrix(inc_mtx);
    SparseMatrix adj_mtx; g.Adjacency_Matrix(adj_mtx);
    SparseMatrix lap_mtx; g.Laplacian_Matrix(lap_mtx);

    cout<<"\nExtra Credit"<<endl;
    cout<<"Adjacency Sparse matrix\n"<<adj_mtx<<endl;
    cout<<"Incidency Sparse matrix\n"<<inc_mtx<<endl;
    cout<<"Laplacian Sparse matrix\n"<<lap_mtx<<endl;
    
#endif
}

// COMBO ASSIGMENT SOLUTION
void Test_HW5()
{
#if HW_ID == 5
    // Practice testing
//    Matrix A(3, 2);
//    Matrix b(3, 1);
//    Matrix x;
//    A = { 2, 0, -1, 1, 0, 2};
//    b = {1,0,-1};
//    A.LeastSquares(b, x);
//    cout<<"Least square approximate solution matrix:\n"<<x<<endl;
//    cout<<"Residual Norm:\n"<<A.residualNorm(x,b)<<endl;
    
    // VMLS 12.10 problem solution
    Matrix A_random(30,10); A_random.randomMatrix(10);
    Matrix b_random(30,1); b_random.randomMatrix(10);
    Matrix x_random;
    
    A_random.LeastSquares(b_random, x_random);
    cout<<"\nLeast square approximate solution matrix:\n"<<x_random<<endl;
    
    cout<<"Residual Norm x:\n"<<A_random.residualNorm(x_random,b_random)<<endl;
    
    Matrix d1(x_random.m, x_random.n); d1.randomMatrix(10);
    cout<<"Residual Norm d1:\n"<<A_random.residualNorm(x_random+d1,b_random)<<endl;

    Matrix d2(x_random.m, x_random.n); d2.randomMatrix(10);
    cout<<"Residual Norm d2:\n"<<A_random.residualNorm(x_random+d2,b_random)<<endl;

    Matrix d3(x_random.m, x_random.n); d3.randomMatrix(10);
    cout<<"Residual Norm d3:\n"<<A_random.residualNorm(x_random+d3,b_random)<<endl;
    
    // VMLS 12.13 problem solution
    Matrix A1_random(20,10); A1_random.randomMatrix(10);
    Matrix b1_random(20,1); b1_random.randomMatrix(10);
    Matrix x1_random;
    A1_random.richardsonAlgorithm(b1_random, x1_random, 500);
    
    cout<<"\nRichardson Algorithm:\n"<<x1_random<<endl;
    
#endif
}

// COMBO ASSIGMENT 2 SOLUTION
void Test_HW6()
{
#if HW_ID == 6
    
    // Practice Testing Gram-Schmidt Algorithm
//    Matrix A(4,3);
//    A = {-1,-1,1,1,3,3,-1,-1,5,1,3,7};
//    Matrix b(4,1); b.randomMatrix(10);
//    Matrix x1;
//    Matrix x2;
//
//    cout<<"A:\n"<<A<<endl;
//    Matrix R(3,3);
//    Matrix Q(4,3);
//    A.GramSchmidt(Q, R);
//    cout<<"Q:\n"<<Q<<endl;
//    cout<<"R:\n"<<R<<endl;
//
//    A.LeastSquares(b, x1);
//    A.QRFactorization(b, x2);
//    cout<<"\nLeast square solution:\n"<<x1<<endl;
//    cout<<"QR Factorization solution:\n"<<x2<<endl;
//
//    Matrix QT(Q.n, Q.m); QT = Q.Transpose();
//
//    // Show that Axˆ = QQTb
//    if (A*x2 == (Q*QT)*b) {
//        cout<<"True, Axˆ= QQTb"<<endl;
//    }
//    else {
//        cout<<"False, Axˆ=/ QQTb"<<endl;
//    }

    // Real Solution
    Matrix A(20,10); A.randomMatrix(100);
    Matrix b(20,1); b.randomMatrix(100);
    Matrix x1;
    Matrix x2;

    Matrix R(10,10);
    Matrix Q(20,10);

    A.Gram_Schmidt(Q, R);
    x2 = A.QRFactorization(Q, R, b);
    A.LeastSquares(b, x1);
        
        cout<<"\nLeast square solution:\n"<<x1<<endl;
        cout<<"QR Factorization solution:\n"<<x2<<endl;

        Matrix QT(Q.n, Q.m); QT = Q.Transpose();
        // Show that Axˆ = QQTb
        if (A*x2 == (Q*QT)*b) {
            cout<<"True, Axˆ= QQTb\n"<<endl;
        }
        else {
            cout<<"False, Axˆ=/ QQTb\n"<<endl;
    }
    
    // Polish the matrix operations and optimize their performance
    // times vary based on each calculation,
    // FOR AVERAGED TIMES CHECK SUBMISSION DOCUMENT ***
    
    using namespace std::chrono;

    vector<int> num;
    num = {10, 100, 1000};
    
    // original version
    for (int i=0; i < num.size(); i++) {
        Matrix F(num[i],num[i]); F.randomMatrix(10);
        Matrix M(num[i],num[i]); M.randomMatrix(10);
        Matrix C(F.m,M.n);
        
        //ORIGINAL VERSION
        // Get starting timepoint
        auto start = high_resolution_clock::now();

        for (int i=0; i<(int)F.m; i++){
            for (int j=0; j<(int)M.n; j++){
              for (int k=0; k<(int)F.n; k++)
                    C(i,j) += F(i,k) * M(k,j);
            }
        }
        // Get ending timepoint
        auto stop = high_resolution_clock::now();
        
        // Get duration
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by original "<<num[i]<<"x"<<num[i]<<": "<<duration.count() << " microseconds" << endl;
        
        //EDITED VERSION
        // Get starting timepoint
        start = high_resolution_clock::now();

        for (int i=0; i<(int)F.m; i++){
            for (int k=0; k<(int)F.n; k++)
                for (int j=0; j<(int)M.n; j++){
                    C(i,j) += F(i,k) * M(k,j);
            }
        }
        // Get ending timepoint
        stop = high_resolution_clock::now();
        
        // Get duration
        duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken by edited "<<num[i]<<"x"<<num[i]<<": "<<duration.count() << " microseconds" << endl;
    }
    
#endif
}

// using the function y = 1+x^2+x^3
double data_function(const double x)  {
    return 1+pow(x,2)+pow(x,3);}

// using the function y = 1+x^2+x^3+x^4
double second_function(double x) {
    return 1+pow(x,2)+pow(x,3)+pow(x,4);
}

void Test_Final()
{
#if HW_ID == 7
    
    // change boolean according to the part of the final project you want to run
    bool appetizer = true;
    bool performance = false;

    bool entree = false;
    
    // DESSERT POINT 1: add a 1 to include bias
    bool dessert1 = false;
    // DESSERT POINT 2: new network architecture
    bool dessert2 = false;
    // DESSERT POINT 3: new reLU function
    bool dessert3 = false;

    
    // change the data path to find files
    string data_path = "/Users/mpmora/Desktop/COSC 70.01/Coding/dartmouth-cs70-starter/MyMatrixLibrary/MNIST_Sub";
    
    if (appetizer) {
        using namespace std::chrono;
           
        Appetizer app;
        // compare the accuracy of different feature models by trying different sets of basis functions.
        double MES1 = app.function_approximation(&data_function);
        cout<<"MES 1:\n"<<MES1<<endl;
        double MES1_b= app.function_approximation(&data_function);
        cout<<MES1_b<<endl;
        
        double MES2 = app.function_approximation1(&second_function);
        cout<<"\nMES 2:\n"<<MES2<<endl;
        double MES2_b = app.function_approximation1(&second_function);
        cout<<MES2_b<<endl;
        
        Matrix train_data(1000,28*28+1);  //// The shape should be (m=n_samples,n=28^2)
        vector<int> train_labels;
        Matrix test_data(200,28*28+1);
        vector<int> test_labels;
        Matrix T(785, 10);
        
        app.openFiles(data_path, train_data, train_labels, test_data, test_labels);
        
        Matrix noise(train_data.m, train_data.n);
        noise.randomMatrix(0.001);
        train_data = train_data+noise;
        
        Matrix noise2 (test_data.m, test_data.n);
        noise2 = noise.Block(0,0, test_data.m, test_data.n);
        test_data = test_data+noise2;
        
        vector<int> approximations;

        // Get starting timepoint
        auto start = high_resolution_clock::now();
        
        if (!performance) {
            T = app.training_QR(train_data, train_labels);
        }
        else {
            T = app.training_LeastSquares(train_data, train_labels);
        }
        
        approximations = app.image_Classification(test_data, T);
        app.calculate_Accuracy(test_labels, approximations);
        
        // Get ending timepoint
        auto stop = high_resolution_clock::now();
               
        // Get duration
        auto duration = duration_cast<microseconds>(stop - start);
        cout << "Time taken: "<<duration.count() << " microseconds" << endl;
               
    }
    
    if (entree) {
        bool architecture = false;
        bool relu = false;
        int bias = 28*28;
        
        if (dessert1) {
            bias = 28*28+1;
        }
        if (dessert2) {
            architecture = true;
        }
        if (dessert3) {
            relu = true;
        }
        
        vector<pair<int, int>> regressor_feature_sizes={{1, 16}, {16, 16}, {16, 16}, {16, 1}};
        Regressor reg(architecture, relu, regressor_feature_sizes,&data_function);

        reg.Train();

        vector<pair<int, int>> classifier_feature_sizes;
        classifier_feature_sizes={{bias, 256}, {256, 256}, {256, 256}, {256, 10}};

        Classifier cls(architecture, relu, bias, data_path, classifier_feature_sizes);

        cls.Train();
    }
#endif
}

int main()
{
    std::cout<<"Hello CS70!"<<std::endl;

    Test_HW1();
    Test_HW2();
    Test_HW3();
    Test_HW4();
    Test_HW5();
    Test_HW6();
    Test_Final();

    system("PAUSE");
}
