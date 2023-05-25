#include <inmost.h>
#include <iostream>
#include <chrono>

using namespace INMOST;
using namespace std;

int main(int argc, char *argv[]) {
    Solver::Initialize(&argc, &argv);
    // Load matrix
    Sparse::Matrix A;
    A.Load("/Users/deti/Desktop/3_курс/6_сем/СВТ/Task_1/2/A_tet1.mtx");
    // Load vector
    Sparse::Vector rhs;
    rhs.Load("/Users/deti/Desktop/3_курс/6_сем/СВТ/Task_1/2/rhs_tet1.mtx");
    // Solution
    Sparse::Vector sol = rhs;
    // Linear solver   
    Solver S(Solver::INNER_ILU2);
    // Set parameters
    S.SetParameter("drop_tolerance", "0.1");
    S.SetParameter("relative_tolerance", "1e-7");
    
    // ILU(2)
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    S.SetMatrix(A);
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    cout << "ILU(2): " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " ms\n";

    // BiCGStab
    begin = chrono::steady_clock::now();
    bool res = S.Solve(rhs, sol);
    end = chrono::steady_clock::now();
    cout << "BiCGStab: " << chrono::duration_cast<chrono::milliseconds>(end - begin).count() << " ms\n";
    
    // Check
    if(!res){
	    cout << "Error: " << S.GetReason() << endl;
	}
	
    cout << "Number of iterations: " << S.Iterations() << endl;
    cout << "Residual: " << S.Residual() << endl;
    Solver::Finalize();
    return 0;
}