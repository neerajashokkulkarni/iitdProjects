#include <iostream>

#include "simplexSolve.h"

using namespace assignment1_MAL704;

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        cout << "usage: <program-name> <path-to-problem-file>" << endl;
        return 1;
    }
    CSimplex problem;
    const bool success = problem.loadProblemFromFile(argv[1]);
    if (success == true)
    {
        problem.set_InitialBFSMethod(USE_BIG_M_METHOD);
        problem.setNumberOfIterations(100);
        problem.solve();
        std::vector<double> X = problem.getOptimumValuesOfOriginalVariables();
        cout<<endl<<"Solution:";
        for(size_t i=0; i< X.size(); i++)
        {
            cout<<" "<<X[i];
        }
        cout<<endl<<"Optimum Objective Function Value: "<< problem.getOptimumValueOfObjectiveFunction();
        cout<<endl;
        return 0;
    }
    else
    {
        return 1;
    }

}