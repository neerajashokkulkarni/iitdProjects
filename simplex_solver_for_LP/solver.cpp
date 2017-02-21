#include "simplexSolve.h"
#include <iostream>

using namespace assignment1_MAL704;



int main(int argc, char* argv[])
{
	CSimplex problem;
	bool status=problem.loadProblemFromFile("./problems/problem0.prob");
	if (status==true)
	{
		problem.set_InitialBFSMethod(USE_BIG_M_METHOD);
		problem.setNumberOfIterations(100);
		problem.solve();
		std::vector<double> X = problem.getOptimumValuesOfOriginalVariables();
		cout<<endl<<"Solution:";
		for( int i=0; i< X.size(); i++)
		{
			cout<<" "<<X[i];
		}
		cout<<endl<<"Optimum Objective Function Value: "<< problem.getOptimumValueOfObjectiveFunction();
		cout<<endl;
	}
	else
		return 1;
	return 0;
}
