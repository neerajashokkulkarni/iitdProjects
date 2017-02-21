#ifndef __SIMPLEX_SOLVE_H__
#define __SIMPLEX_SOLVE_H__

#include<Eigen/Geometry>
#include <iostream>
#include <vector>


namespace assignment1_MAL704
{
    /*using*/
    using std::cout;
    using std::cerr;
    using std::endl;
    using std::vector;
    using Eigen::MatrixXd;
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    class CSimplex;
    enum initialBFSMethodEnum{USE_TWO_PHASE_METHOD, USE_BIG_M_METHOD};
    class CSimplex
    {
        private:
            unsigned int MAX_ITERATIONS;
            enum {STATUS_EMPTY, STATUS_PROBLEM_READ, STATUS_PROBLEM_IN_STD_FORM, STATUS_PHASE_ONE_PROBLEM_CONSTRUCTED, STATUS_PROBLEM_SOLVED} status;
            enum {PROBLEM_MINIMIZATION, PROBLEM_CONVERTED_TO_MINIMIZATION} problemType;
            /* A note about PROBLEM_CONVERTED_TO_MINIMIZATION
             * ---------------------------------------------
             * Note that there is no PROBLEM_MAXIMIZATION listed in problemType enum. Our solve() method basically 
             * solves strictly a minimization problem. But what if the input problem is originially a maximization 
             * problem? We convert it to a minimization problem by simply negating all the objective function
             * coefficients. As we proceed and find the optimum solution, we check whether the problemType is 
             * PROBLEM_CONVERTED_TO_MINIMIZATION, if yes, we simply return the negative value of the optimum 
             * objective function value. As simple as that! That's how we 'solve' a maximization problem.*/
             
            unsigned int nConstraints, nVariables;
            VectorXd vec_coefficientsObj; //Objective function coefficients
            MatrixXd mat_A; //constraints matrix A
            VectorXd vec_b; //constraints RHS matrix b
            typedef enum {SIGN_EQUAL, SIGN_GREATERTHANOREQUAL, SIGN_LESSERTHANOREQUAL} ConstrainSign;
            std::vector<ConstrainSign> constraints_signs;
            initialBFSMethodEnum initialBFSMethod;
            
            /*'current-iteration' specific data*/
            vector<int> current_BasisIndices, current_NonbasisIndices; //the basis indices (Bi) and non-basis indices (Ni)
            MatrixXd current_BasisMatrix, current_InverseBasisMatrix, current_NonbasisMatrix; //the basis matrix B, its inverse and non-Basis matrix (N)
            VectorXd current_XB;//The basic solution X_B
            VectorXd current_XN;//The non-basic solution X_N
            VectorXd current_CB;//The basic solution X_B
            VectorXd current_CN;//The basic solution X_B
            VectorXd current_ys; //Ys
            VectorXd current_SN; //Pricing
            VectorXd current_d; //used for computing leaving index while updating Basis and also to update XB at the end of every iteration.
            double current_Xq_plus; //scalar update for newer XB, used while updating current_XB at the end of every iteration.
            
            /*Associated Equivalent Problems*/
            CSimplex *standardFormLPP;
            CSimplex *phaseOneLPP;
    
    
        protected:
            void allocateFundamentalDataStructures();
            void allocateIterationDataStructures();
            void conservativeResizeOfIterationDataStructures();
            void convertToStandardForm();
            void computeInitialBasis();
            void assertFeasiblityOfCurrentBasicSolution() const;
            void computeBasisMatrix();
            void computeCurrrent_CB_and_CN();
            
            int chooseIndexOfEnteringVariable();
            int chooseIndexOfLeavingVariable(int const& enteringIndex);
            
            void addPositiveVariable(unsigned int constraintID); //adds slack variable
            void addNegativeVariable(unsigned int constraintID); //adds surplus variable
            
            /*helper functions*/
            void clearAll();
            bool isProblemInStandardForm() const;
            double evaluateObjectiveFunctionValatCurrentPoint() const;
            
        public:
            CSimplex() //constructor
            {
                this->status = STATUS_EMPTY;
                problemType = PROBLEM_MINIMIZATION;
                this->nConstraints = this->nVariables = 0;
                standardFormLPP = NULL;
                phaseOneLPP = NULL;
                current_Xq_plus = 0;
                MAX_ITERATIONS = 100; //default
                initialBFSMethod = USE_BIG_M_METHOD; //default
            }
            void copyFrom(CSimplex const* other)
            {
                this->problemType       = other->problemType;
                this->nVariables        = other->nVariables;
                this->nConstraints      = other->nConstraints;
                this->allocateFundamentalDataStructures();
                this->vec_coefficientsObj = other->vec_coefficientsObj;
                this->mat_A             = other->mat_A;
                this->vec_b             = other->vec_b;
                this->constraints_signs = other->constraints_signs;
                this->standardFormLPP   = other->standardFormLPP;
                this->status            = other->status;
            }
            bool loadProblemFromFile(const char filename[]);
            void set_InitialBFSMethod(int in_method);
            void setNumberOfIterations(int in_n_iterations);
            void solve();
            std::vector<double> getOptimumValuesOfOriginalVariables() const;
            double getOptimumValueOfObjectiveFunction() const;

    };
    
}





#endif //__SIMPLEX_SOLVE_H__
