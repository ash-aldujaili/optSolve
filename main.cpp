/*
 * main.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: ash
 */
#include"Dictionary.hpp"
#include<iostream>

int main(int argc, char** argv)
{
	// construct the dictionary
	/*
	Dictionary dict(argv[1]);
	// view it
	print(dict);
	// solve the dict:
	if (dict.pivotToFinalDict(false))
	{
		std::cout.precision(5);
		std::cout<<dict.getObjVal()<<"\n";
		std::cout<<dict.getNumPivots()<<"\n";
	}
	else
		std::cout<<"UNBOUNDED\n";*/
	// Anohter approach for initialization:

	// get the feasible dictionary
	// construct the dictionary
	//Dictionary dict2(argv[1]);
	//dict2.solve();

	Dictionary dict(argv[1]);
	dict.solveILP();
	//dict2.toFeasibleDict(true);
	/*
	if (dict2.pivotToFinalDict(false))
		{
			std::cout.precision(5);
			std::cout<<dict2.getObjVal()<<"\n";
			//std::cout<<dict2.getNumPivots()<<"\n";
		}
		else
			std::cout<<"UNBOUNDED\n";
			*/
	// pivot
	//dict.pivot();
	//std::cout<<dict.getObjVal()<<"\n";
	//print(dict2);



}


