#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include<Eigen/Geometry>
#include <iostream>
#include <vector>


namespace assignment1_MAL704
{
template <typename mytype> void printVector(std::vector<mytype> vec_in);
VectorXd singleColumnOf(MatrixXd A, unsigned int index);

using Eigen::VectorXd;
using Eigen::MatrixXd;
	
template <typename mytype>
void printVector(std::vector<mytype> vec_in)
{
	typename std::vector< mytype>::iterator i;
	for(i=vec_in.begin(); i!=vec_in.end(); i++)
	{
		std::cout<<"\t"<<*i;
	}
}

VectorXd singleColumnOf(MatrixXd A, unsigned int index)
{
	VectorXd outColumn(A.rows());
	for(unsigned int i=0; i<A.rows(); i++)
	{
		outColumn(i) = A(i,index);
	}
	return outColumn;
}
}

#endif //__UTILS_HPP__
