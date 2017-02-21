#include "simplexSolve.h"
#include "utils.hpp"
#include <iostream>
#include <float.h>
#include <fstream>
#include <string>

bool FLAG_VERBOSE = false;

namespace assignment1_MAL704
{
    
void CSimplex::allocateFundamentalDataStructures()
/*Once we have read nVariables and nConstraints from the problem file,
 *we allocate the data-structures with the appropriate sizes.*/
{
    assert(this->nVariables > 0);
    assert(this->nConstraints > 0);
    this->vec_coefficientsObj = VectorXd(this->nVariables);
    this->mat_A = MatrixXd(this->nConstraints, this->nVariables);
    this->vec_b = VectorXd(this->nConstraints);
    this->constraints_signs.resize(this->nConstraints);
}

void CSimplex::allocateIterationDataStructures()
{
    assert(this->nVariables > 0);
    assert(this->nConstraints > 0);

    /*allocate all iteration-specific data structures to right sizes*/
    current_BasisIndices.resize(nConstraints); //the basis indices
    current_NonbasisIndices.resize(nVariables-nConstraints);//the non-basis indices N
    current_BasisMatrix.resize(nConstraints,nConstraints); //the basis matrix
    current_InverseBasisMatrix.resize(nConstraints,nConstraints); //the inverse basis matrix
    current_NonbasisMatrix.resize(nConstraints,nVariables-nConstraints); //the non-basis matrix
    current_XB.resize(nConstraints);//The basic solution
    current_XN.resize(nVariables-nConstraints);//The non-basic solution
    current_CB.resize(nConstraints);//The basic-indexed coefficients
    current_CN.resize(nVariables-nConstraints);//The non-basic-indexed coefficients
    current_ys.resize(nConstraints);
    current_SN.resize(nConstraints);
}


bool CSimplex::loadProblemFromFile(const char filename[])
{
    /*Part I: Read the problem in our member variables*/
    std::ifstream probfile;
    probfile.open(filename);
    if(!probfile.is_open()){
        cerr<<endl<<"Couldn't open input problem file "<<filename<<" !";
        return false;
    }
    
    if(FLAG_VERBOSE) cout<<endl<<"Loading LP Problem from "<<filename;
    char *token;
    std::string line;
    unsigned int stuffReadFromFile_counter=0;
    //~ vector<std::string> tokens;
    std::string mdata_string("[METADATA]");
    std::string objfn_string("[OBJECTIVE]"); 
    std::string constr_string("[CONSTRAINTS]");
    while(probfile.good())
    {
        getline(probfile,line);
        if (line.compare(mdata_string) == 0)
        /*Read Matadata (nVariables and nConstraints)*/
        {
            stuffReadFromFile_counter++;
            assert(probfile.good());
            getline(probfile,line);
            token = strtok((char*)line.c_str()," ,");
            assert(std::string(token).compare("vars")==0);
            token = strtok(NULL, " ,");
            this->nVariables=atoi(token);
            assert(this->nVariables > 0);
            assert(probfile.good());
            getline(probfile,line);
            token = strtok((char*)line.c_str()," ,");
            assert(std::string(token).compare("constraints")==0);
            token = strtok(NULL, " ,");
            assert(token != NULL);
            this->nConstraints=atoi(token);
            assert(this->nConstraints > 0);
            //now that we know both nVars and nConstraints...
            this->allocateFundamentalDataStructures();
        }
        
        if (line.compare(objfn_string) == 0)
        /*Read Objective function*/
        {
            stuffReadFromFile_counter++;
            assert(probfile.good());
            getline(probfile,line);
            token = strtok((char*)line.c_str()," ,");
            if(std::string(token).compare("minimize")==0)
            {
                this->problemType = PROBLEM_MINIMIZATION;
            }
            else if(std::string(token).compare("maximize")==0)
            {
                this->problemType = PROBLEM_CONVERTED_TO_MINIMIZATION;
                //refer to the 'PROBLEM_CONVERTED_TO_MINIMIZATION' note in simplexSolve.h.
            }
            else assert(!"Invalid keyword (problem type) in the input problem file!");
            for(int i=0; i<this->nVariables; i++)
            {
                token = strtok(NULL, " ,");
                assert(token != NULL);
                this->vec_coefficientsObj[i] = atof(token);
                if (this->problemType == PROBLEM_CONVERTED_TO_MINIMIZATION)
                //refer to the 'PROBLEM_CONVERTED_TO_MINIMIZATION' note in simplexSolve.h.
                {
                    this->vec_coefficientsObj[i] *= -1.0;
                }
            }
        }
        if (line.compare(constr_string) == 0)
        {
            stuffReadFromFile_counter++;
            for(int i=0; i< this->nConstraints; i++)
            //each ith constraint
            {
                assert(probfile.good());
                getline(probfile,line);
                token = strtok((char*)line.c_str()," ,");
                for(int j=0; j < this->nVariables; j++)
                {
                    assert(token != NULL);
                    this->mat_A(i,j) = atof(token);
                    token = strtok(NULL, " ,");
                }
                assert(token != NULL);
                if(std::string(token).compare("=")==0)
                {
                    this->constraints_signs[i] = SIGN_EQUAL;
                }
                else if(std::string(token).compare(">=")==0)
                {
                    this->constraints_signs[i] = SIGN_GREATERTHANOREQUAL;
                }
                else if(std::string(token).compare("<=")==0)
                {
                    this->constraints_signs[i] = SIGN_LESSERTHANOREQUAL;
                }
                else assert(!"invalid constraint sign in the input file");
                // and Now read 'b' (the RHS of the current i'th constraint)
                token = strtok(NULL, " ,");
                assert(token != NULL);
                this->vec_b(i) = atof(token);
            }
        }
    }
    assert(stuffReadFromFile_counter == 3); //3 things = metadata+objective+constraints
    probfile.close();
        
    if(FLAG_VERBOSE) cout<<endl<<"C = "<<endl<< this->vec_coefficientsObj;
    if(FLAG_VERBOSE) cout<<endl<<"A = "<<endl<< this->mat_A;
    if(FLAG_VERBOSE) cout<<endl<<"b = "<<endl<< this->vec_b;
    
    
    if(FLAG_VERBOSE) cout<<endl<<"assert(A has full Row Rank)";
    //~ assert(mat_A.rank() == mat_A.rows());
    this->status = STATUS_PROBLEM_READ;
    if(FLAG_VERBOSE) cout<<endl<<"loaded LPP.";
    return true;
}


void CSimplex::solve()
/*Our algorithm is based upon:  Procedure 13.1, Page 370, section 13.3, Chapter 13 ('Linear Programming: The Simplex Method')
 *                              from the book 'Numerical Optimization' by Nocedal & Wright.*/
{
    //assert our assumptions first
    assert(MAX_ITERATIONS > 0);
    
    /*convert To Standard Form Linear Programming Problem*/
    //if the original problem is already in standard form then we just
    //make 'this->standardFormPP' point to 'this' itself. If the current
    //problem is not in standard form, then we convert it to standard form
    //and make 'this->standardFormPP' point to that new standard-form problem.
    //In any case, once we return from convertToStandardForm(), the 
    //this->standardFormPP should be used when we proceed further to solve.
    //Note that we have avoided any recursive structure (calling this->
    //standardFormLPP->solve() from this->solve()), to avoid confusion.
    convertToStandardForm();
    assert(this->status == STATUS_PROBLEM_IN_STD_FORM);


    /*The main optimization-iterations*/
    this->standardFormLPP->allocateIterationDataStructures();
    int enteringNonbasisIndex=-1, leavingBasisIndex=-1;
    int enteringVariableIndex=-1, leavingVariableIndex=-1;
    for(unsigned int i=0; i< MAX_ITERATIONS ; ++i)
    {
        if(FLAG_VERBOSE) cout<<endl<<"----------------Iteration "<<i<<"---------------------------";
        if(i==0)
        {
            //compute the initial Basis Indices.
            standardFormLPP->computeInitialBasis();
            this->standardFormLPP->conservativeResizeOfIterationDataStructures();
            //the conservative resize of iteration-data-structures is needed 
            //after Initial-basis computation as there 
            //exists a possibility of on-the-fly addition of artificial
            //variables during the initial-basis-computation procedure.
        }
        
        /*Step 0: According to basis indices, compute the basis matrix and other relevant quantities.*/
        standardFormLPP->computeBasisMatrix();
        if(FLAG_VERBOSE) {cout<<endl<<"BasisIndices :\n"; printVector(standardFormLPP->current_BasisIndices);}
        if(FLAG_VERBOSE) {cout<<endl<<"NonbasisIndices :\n"; printVector(standardFormLPP->current_NonbasisIndices);}
        if(FLAG_VERBOSE) cout<<endl<<"Basis :\n"<<standardFormLPP->current_BasisMatrix;
        if(FLAG_VERBOSE) cout<<endl<<"BasisInverse :\n"<<standardFormLPP->current_InverseBasisMatrix;
        if(FLAG_VERBOSE) cout<<endl<<"NonBasis :\n"<<standardFormLPP->current_NonbasisMatrix;
        //secondly compute initial Basic Feasible solution
        standardFormLPP->current_XB = standardFormLPP->current_InverseBasisMatrix*standardFormLPP->vec_b;
        if(FLAG_VERBOSE) cout<<endl<<"XB :\n"<<standardFormLPP->current_XB;
        if(FLAG_VERBOSE) cout<<endl<<"XN :\n"<<standardFormLPP->current_XN;
        standardFormLPP->assertFeasiblityOfCurrentBasicSolution();
        //thirdly compute all Ys accordingly 
        standardFormLPP->computeCurrrent_CB_and_CN();
        standardFormLPP->current_ys = standardFormLPP->current_InverseBasisMatrix.transpose()*standardFormLPP->current_CB;
        if(FLAG_VERBOSE) cout<<endl<<"current_ys:\n"<<standardFormLPP->current_ys;
        
        
        /*Step 1: with the updated BFS solve for Sn*/
        standardFormLPP->current_SN = standardFormLPP->current_CN-standardFormLPP->current_NonbasisMatrix.transpose()*standardFormLPP->current_ys;
        if(FLAG_VERBOSE) cout<<endl<<"Sn :\n"<<standardFormLPP->current_SN;
        
        /*Step 2: entering variable-- given the new SN, see if we have reached an 
         * optimum. If no, choose an entering variable. */
        enteringNonbasisIndex = standardFormLPP->chooseIndexOfEnteringVariable();
        if( enteringNonbasisIndex == -1)
        {
            cout<<endl<<"------------------------------------------------------";
            cout<<endl<<endl<<"Reached an OPTIMUM solution in iteration "<<i<<" !!";
            cout<<endl<<"Object Value at opt solution: "<<standardFormLPP->evaluateObjectiveFunctionValatCurrentPoint();
            cout<<endl<<"------------------------------------------------------"<<endl;
            this->status = STATUS_PROBLEM_SOLVED;
            standardFormLPP->status = STATUS_PROBLEM_SOLVED;
            return;
        }
        enteringVariableIndex = standardFormLPP->current_NonbasisIndices[enteringNonbasisIndex];
        if(FLAG_VERBOSE) cout<<endl<<"enteringVariableIndex = "<<enteringVariableIndex;
        
        /*Step 3: leaving variable-- choose the leaving variable index and in the 
         *process detect unboundedness. */
        leavingBasisIndex = standardFormLPP->chooseIndexOfLeavingVariable(enteringVariableIndex);
        if( leavingBasisIndex == -1)
        {
            cout<<endl<<"------------------------------------------------------";
            cout<<endl<<endl<<"Detected UNBOUNDEDNESS at iteration "<<i<<" !!";
            cout<<endl<<"------------------------------------------------------";
            return;
        }
        leavingVariableIndex = standardFormLPP->current_BasisIndices[leavingBasisIndex];
        if(FLAG_VERBOSE) cout<<endl<<"leavingVariableIndex  = "<<leavingVariableIndex;
        
        /*step 4: Now accordingly update the basis, and update the relevant members.*/
        standardFormLPP->current_BasisIndices[leavingBasisIndex] = enteringVariableIndex;
        standardFormLPP->current_NonbasisIndices[enteringNonbasisIndex] = leavingVariableIndex;
        standardFormLPP->current_XB = standardFormLPP->current_XB - standardFormLPP->current_d*standardFormLPP->current_Xq_plus;
        standardFormLPP->current_XN.setZero(standardFormLPP->nVariables-standardFormLPP->nConstraints);
        //~ standardFormLPP->current_XN[enteringNonbasisIndex] = standardFormLPP->current_Xq_plus;
    }
    
    
    
}


int CSimplex::chooseIndexOfLeavingVariable(int const& enteringIndex)
{
    assert(enteringIndex >= 0);
    int out_leavingBasisIndex=-1;
    /*Solve B*d = Aq for d*/
    this->current_d = current_InverseBasisMatrix*singleColumnOf(standardFormLPP->mat_A, enteringIndex);
    if (FLAG_VERBOSE) cout<<endl<<"d=\n"<<this->current_d;
    
    /*Find the index that minimizes the XB/d ratio*/
    int minRatio_index = -1;
    double minRatio_value=DBL_MAX,ratio=0.0;
    unsigned int nNegativeElements=0;
    for(int i=0; i<this->current_d.size(); i++){
        if(current_d(i) <= 0)
        {
            nNegativeElements++; //required for the unboundedness check below
        }
        else //consider only positive quotients
        {
            assert(current_d(i) > 0);
            ratio = current_XB(i)/current_d(i);
            if (ratio < minRatio_value)
            {
                minRatio_value = ratio;
                minRatio_index = i;
            }
        }
    }
    if (nNegativeElements == current_d.size()) //all elements were negative
    /*Unbounded Problem detected*/
    {
        if (FLAG_VERBOSE) cout<<endl<<endl<<"Unbounded problem."<<endl;
        return -1;
    }
    assert(minRatio_index >= 0); //otherwise the problem is unbounded already.
    out_leavingBasisIndex = minRatio_index;
    this->current_Xq_plus = minRatio_value;
    return  out_leavingBasisIndex;
}


void CSimplex::addPositiveVariable(unsigned int constraintID)
/*Adds a "positive" variable to the specified constraint. Note that,
 *by "positive" we mean that its coefficient in the specified constraint 
 *will be +1. 
 *It should be noted that, in fact, any variable we add is always and always
 * going to be non-negative (>=0), because of the non-negativity constraint
 * on each variable of a standard LPP.*/
{
    this->nVariables++;
    
    //update mat_A 
    this->mat_A.conservativeResize(this->nConstraints, this->nVariables);
    //conservativeResize() in Eigen library resizes with data preservation
    for(int i=0; i< this->mat_A.rows(); i++)
    {
        if( i==constraintID)
        {
            mat_A(i,nVariables-1) = 1;
        }
        else
        {
            mat_A(i,nVariables-1) = 0;
        }
    }
    
    // update coefficients
    this->vec_coefficientsObj.conservativeResize(nVariables);
    this->vec_coefficientsObj(nVariables-1) = 0;
}


void CSimplex::addNegativeVariable(unsigned int constraintID)
/*Adds a "negative" variable to the specified constraint. Note that,
 *by "negative" we mean that its coefficient in the specified constraint 
 *will be -1.
 *It should be noted that, in fact, any variable we add is always and always
 * going to be non-negative (>=0), because of the non-negativity constraint
 * on each variable of a standard LPP.*/
{
    this->nVariables++;
    
    //update mat_A 
    this->mat_A.conservativeResize(this->nConstraints, this->nVariables);
    //conservativeResize() in Eigen library resizes with data preservation
    for(int i=0; i< this->mat_A.rows(); i++)
    {
        if( i==constraintID)
        {
            mat_A(i,nVariables-1) = -1;
        }
        else
        {
            mat_A(i,nVariables-1) = 0;
        }
    }
    
    // update coefficients
    this->vec_coefficientsObj.conservativeResize(nVariables);
    this->vec_coefficientsObj(nVariables-1) = 0;
}


int CSimplex::chooseIndexOfEnteringVariable()
/*Checks values of the member vector current_SN and accordingly chooses
 *which index should enter the current basis.
 *It returns index (-1) if no such eligible entering variable is possible. 
 *That also means we have reached an optimum point already.*/
{
    unsigned int i=0;
    int outNonbasisIndex=-1;
    int minNegCost=0;
    for(i=0; i< this->current_SN.size(); i++)
    {
        if(current_SN[i] < 0 && current_SN[i]<minNegCost){
            minNegCost=current_SN[i];
            outNonbasisIndex=i;
        }
    }
    
    return outNonbasisIndex;
}


void CSimplex::computeCurrrent_CB_and_CN()
/*Given the basis indices Bi[], we equate CB to corresponding objective-function coefficients, in the same order.*/
{
    //compute current_CB
    assert(current_CB.size() == nConstraints);
    assert(current_BasisIndices.size() == nConstraints);
    for(unsigned int i=0;i<current_BasisIndices.size(); ++i){
        current_CB[i] = vec_coefficientsObj[current_BasisIndices[i]]; 
    }
    if(FLAG_VERBOSE) cout<<endl<<"current_CB:\n"<<current_CB;
    //compute current_CN
    assert(current_CN.size() == nVariables-nConstraints);
    assert(current_NonbasisIndices.size() == nVariables-nConstraints);
    for(unsigned int i=0;i<current_NonbasisIndices.size(); ++i){
        current_CN[i] = vec_coefficientsObj[current_NonbasisIndices[i]]; 
    }
    if(FLAG_VERBOSE) cout<<endl<<"current_CN:\n"<<current_CN;
}


void CSimplex::computeBasisMatrix()
/*Given the basis indices Bi[], we equate B to corresponding columns of A, in the same order.*/
{
    assert(current_BasisMatrix.rows() == nConstraints);
    assert(current_BasisMatrix.cols() == nConstraints);
    assert(current_NonbasisMatrix.rows() == nConstraints);
    assert(current_NonbasisMatrix.cols() == nVariables-nConstraints);
    
    //build basis matrix
    for(unsigned int i=0;i< current_BasisIndices.size(); ++i)
    {
        /*copy Bi[i]th column of A into ith column of B*/
        for(unsigned int j=0; j< mat_A.rows(); ++j)
        {
            current_BasisMatrix(j,i) = mat_A(j,current_BasisIndices[i]);
        }
    }
    //build non-basis matrix
    for(unsigned int i=0;i< current_NonbasisIndices.size(); ++i)
    {
        /*copy Bi[i]th column of A into ith column of B*/
        for(unsigned int j=0; j< mat_A.rows(); ++j)
        {
            current_NonbasisMatrix(j,i) = mat_A(j,current_NonbasisIndices[i]);
        }
    }   
    //build inverse of basis matrix
    current_InverseBasisMatrix = current_BasisMatrix.inverse();
}



bool CSimplex::isProblemInStandardForm() const
{
    /*Sign of every constraint should be equality('=')*/
    for(int i=0; i< constraints_signs.size(); i++)
    {
        if(constraints_signs[i] != SIGN_EQUAL)
            return false;
    }
    return true;
}

void CSimplex::clearAll()
{
    this->nVariables = this->nConstraints = 0;
    //now, accordingly zero-sized-allocate all data structures
    this->vec_coefficientsObj.resize(0);
    this->mat_A.resize(0, 0);
    this->vec_b.resize(0);
    this->constraints_signs.resize(0);
    
    /*allocate all iteration-specific data structures to right sizes*/
    current_BasisIndices.resize(0); //the basis indices
    current_NonbasisIndices.resize(0);//the non-basis indices N
    current_BasisMatrix.resize(0,0); //the basis matrix
    current_InverseBasisMatrix.resize(0,0); //the inverse basis matrix
    current_NonbasisMatrix.resize(0,0); //the non-basis matrix
    current_XB.resize(0);//The basic solution
    current_XN.resize(0);//The non-basic solution
    current_CB.resize(0);//The basic-indexed coefficients
    current_CN.resize(0);//The non-basic-indexed coefficients
    current_ys.resize(0);
    current_SN.resize(0);
    //~ this->status = STATUS_EMPTY;
}

void CSimplex::convertToStandardForm()
/*This function converts the given problem (stored in 'this'), into a standard form Linear
 * Programming Problem. The given problem may have its constraints in any form (>= or <=
 * inequalities or equality), so this function produces another equivalent problem with all 
 * equality constraints by suitably appending additional slack/surplus variables to individual 
 * constraints.
 * Writes: this->standardFormLPP.*/
{
    
    if(status == STATUS_PROBLEM_IN_STD_FORM || this->standardFormLPP != NULL){
        if(FLAG_VERBOSE) cout<<endl<<"Already CONVERTED to standard form.";
        return;
    }
    
    if (isProblemInStandardForm())
    {
        if(FLAG_VERBOSE) cout<<endl<<"Problem is already in standard form!";
        this->status = STATUS_PROBLEM_IN_STD_FORM;
        /*if the current problem is itself in standard form, then we
         *just set this->standardFormLPP = this. That is, a self-loop.*/
        standardFormLPP = this; //self-loop
        return;
    }
    else
    {
        assert(standardFormLPP == NULL);
        standardFormLPP = new CSimplex;
        standardFormLPP->copyFrom(this);
        standardFormLPP->standardFormLPP = standardFormLPP; //self-loop.
        for(int i=0; i< standardFormLPP->nConstraints; i++)
        {
            switch(standardFormLPP->constraints_signs[i])
            {
                case SIGN_EQUAL:
                break;
                    
                case SIGN_LESSERTHANOREQUAL:
                /*Add a slack variable to current constraint if necessary*/
                 standardFormLPP->addPositiveVariable(i);
                break;
                
                case SIGN_GREATERTHANOREQUAL:
                /*Add a surplus variable to current constraint if necessary*/
                 standardFormLPP->addNegativeVariable(i);
                break;
                
                default: assert(0);//impossible.
            }
        }
    }
    if(FLAG_VERBOSE) cout<<endl<<"C = "<<endl<< standardFormLPP->vec_coefficientsObj;
    if(FLAG_VERBOSE) cout<<endl<<"A = "<<endl<< standardFormLPP->mat_A;
    if(FLAG_VERBOSE) cout<<endl<<"b = "<<endl<< standardFormLPP->vec_b;
    
    /*Verify that our standard assumption of m<n holds on this standard problem.*/
    assert(standardFormLPP->mat_A.rows() < standardFormLPP->mat_A.cols());
    this->status = this->standardFormLPP->status = STATUS_PROBLEM_IN_STD_FORM;
    this->standardFormLPP->standardFormLPP = standardFormLPP;//self-loop
}


void CSimplex::computeInitialBasis()
/*This method searches for all columns of the Identity matrix in the constraint matrix 'A' to 
 *do the first estimate of the Basic Feasible Solution(BFS). In this process, if such columns
 *are not found in matrix A, it also may add additional 'artificial' variables, and hence append
 *additional respective missing Identity columns to A.
 *Writes: (i) First Basis indices: Indices corresponding to Indentity columns.
 *       (ii) may add additional missing Identity columns to A.
 *       (iii) And also, current_NonbasisIndices[].
 * */
{
    assert(status == STATUS_PROBLEM_IN_STD_FORM);
    assert(this->current_BasisIndices.size() == nConstraints);
    assert(this->current_NonbasisIndices.size() == nVariables-nConstraints);
    unsigned int i,j,k,nb=0;
    for (i=0; i< this->current_BasisIndices.size(); i++)
    {
        if(FLAG_VERBOSE) cout<<endl<<"-searching for "<<i<<"th column of Identity Matrix.";
        for(j=0; j<mat_A.cols(); ++j)
        {
            if(FLAG_VERBOSE) cout<<endl<<"\t-matching against column "<<j<<".";
            /*match column contents*/
            for (k=0; k< mat_A.rows(); ++k)
            {
                if( (k == i && mat_A(k,j) != 1) ||
                    (k != i && mat_A(k,j) != 0) )//mismatch
                {
                    if(FLAG_VERBOSE) cout<<"\t..failed.";
                    break;
                }
            }
            if (k== mat_A.rows()){
                if(FLAG_VERBOSE) cout<<"\t..success";
                current_BasisIndices[i] = j;
                break;
            }
        }
        
        if(j== mat_A.cols())
        {
            if(FLAG_VERBOSE) cout<<endl<<"-Tried all columns: We need to introduce artificial variables";
            if (this->initialBFSMethod == USE_BIG_M_METHOD)
            {
                /*add a "positive" artificial variable with very high coefficient value.*/
                if(FLAG_VERBOSE) cout<<endl<<"-Adding a Big-M artificial variable.";
                this->addPositiveVariable(i);
                this->vec_coefficientsObj(this->nVariables-1) = DBL_MAX;
                current_BasisIndices[i] = this->nVariables-1;
            }
            else
            {
                assert(!"ASSERT: Two Phase Method is not supported by the source yet.");
            }
        }
    }
    
    /*now initialize non-basis indices*/
    for(i=0,j=0,k=0; j< current_BasisIndices.size(); j++,i++)
    {
        while(i!= current_BasisIndices[j])
        {
            current_NonbasisIndices[k++] = i;
            i++;
        }
    }
    assert(k == nVariables-nConstraints);
}


void CSimplex::assertFeasiblityOfCurrentBasicSolution() const
{
    for(unsigned int i=0; i< current_XB.size(); i++)
    {
        assert(current_XB[i] >= 0);
    }
    
    if(FLAG_VERBOSE) cout<<endl<<"-Current Basic Solution Feasibility Check: OK";   
}


double CSimplex::evaluateObjectiveFunctionValatCurrentPoint() const
{
    double objValue=0.0;
    for(unsigned int i=0; i<current_BasisIndices.size(); i++){
        objValue += current_XB[i]*vec_coefficientsObj[current_BasisIndices[i]];
    }
    /*contribution by the non-basic part is zero since current_XN is always zero
     *(hence we don't even store or keep track of it).
    */
    if (this->problemType == PROBLEM_CONVERTED_TO_MINIMIZATION){
    //refer to the 'PROBLEM_CONVERTED_TO_MINIMIZATION' note in simplexSolve.h.
        return (-1.0)*objValue;
    }
    else{
        return objValue;
    }
}


void CSimplex::set_InitialBFSMethod(int in_method)
{
    assert(in_method == USE_BIG_M_METHOD || in_method == USE_TWO_PHASE_METHOD);
    
    if(in_method == USE_TWO_PHASE_METHOD)
    {
        assert(!"ASSERT: Two Phase Method is not supported by the source yet.");
    }
    
    this->initialBFSMethod = initialBFSMethodEnum(in_method);
    
}

void CSimplex::conservativeResizeOfIterationDataStructures()
{
    assert(this->nVariables > 0);
    assert(this->nConstraints > 0);

    /*allocate all iteration-specific data structures to right sizes*/
    current_BasisIndices.resize(nConstraints); //the basis indices
    current_NonbasisIndices.resize(nVariables-nConstraints);//the non-basis indices N
    current_BasisMatrix.conservativeResize(nConstraints,nConstraints); //the basis matrix
    current_InverseBasisMatrix.conservativeResize(nConstraints,nConstraints); //the inverse basis matrix
    current_NonbasisMatrix.conservativeResize(nConstraints,nVariables-nConstraints); //the non-basis matrix
    current_XB.conservativeResize(nConstraints);//The basic solution
    current_XN.conservativeResize(nVariables-nConstraints);//The non-basic solution
    current_CB.conservativeResize(nConstraints);//The basic-indexed coefficients
    current_CN.conservativeResize(nVariables-nConstraints);//The non-basic-indexed coefficients
    current_ys.conservativeResize(nConstraints);
    current_SN.conservativeResize(nConstraints);
}

std::vector<double> CSimplex::getOptimumValuesOfOriginalVariables() const
{
    assert(this->status == STATUS_PROBLEM_SOLVED);
    std::vector<double> outVec;
    outVec.resize(this->nVariables); //this->nVariables is the original number of variables
    for(unsigned int i=0; i < standardFormLPP->current_BasisIndices.size(); i++)
    {
        if(standardFormLPP->current_BasisIndices[i] < this->nVariables) //only original variables
        {
            outVec[standardFormLPP->current_BasisIndices[i]] = standardFormLPP->current_XB(i);
        }
    }
    for(unsigned int i=0; i < standardFormLPP->current_NonbasisIndices.size(); i++)
    {
        if(standardFormLPP->current_NonbasisIndices[i] < this->nVariables)//only original variables
        {
            outVec[standardFormLPP->current_NonbasisIndices[i]] = standardFormLPP->current_XN(i);//is going to be zero as it is
        }
    }
    return outVec;
}

void CSimplex::setNumberOfIterations(int in_n_iterations)
{
    assert(in_n_iterations> 0);
    this->MAX_ITERATIONS = in_n_iterations;
}


double CSimplex::getOptimumValueOfObjectiveFunction() const
{
    assert(this->status == STATUS_PROBLEM_SOLVED);
    return standardFormLPP->evaluateObjectiveFunctionValatCurrentPoint();
}


}//namespace assignment1_MAL704 over

