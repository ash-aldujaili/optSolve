/*
 * Dictionary.hpp
 *
 *  Created on: Nov 1, 2014
 *      Author: Abdullah Al-Dujaili
 */

#ifndef DICTIONARY_HPP_
#define DICTIONARY_HPP_

#include<vector>
#include<string>

class Dictionary
{

	public:
		// constructor
		Dictionary(const std::string fileName);
		// copy constructor * currently it supports copying for dulaizing puroposse only
		Dictionary(const Dictionary& dict);
		// destructor
		~Dictionary();
		// performs a single pivot step
		int pivot(const bool isVerbose = true);
		// get the final dictionary
		bool pivotToFinalDict(const bool isVerbose = false);
		// return the current objective value
		float getObjVal() const { return z0; };
		// initialize the dictionary to get a feasible dictionary
		// return true if the a feasible dictionary found, false otherwize
		bool toFeasibleDict(bool isVerbose = false, bool isObjChng = true);
		int getNumPivots() const {return numPivots;}
		// test if the dictionary is final
		bool isFinal() const;
		// return UNBOUNDED, INFEASIBLE or the achieved objective value:
		void solve();
		// return UNBOUNDED, INFEASIBLE or the achieved objective value with integer constraints
		void solveILP();
		// print dictionary
		friend void print(const Dictionary &dict);
	private:
		int m , n; 	// number of variables, constraints
		int numPivots; // number of pivots to the final
		std::vector<int> basicVar; 	// list of basic indices [m integers]
		std::vector<int> nBasicVar; // list of non-basic indices
		float* b;   // m floating point numbers
		float* a; 	// A matrix
		float** aRow; // rows of A matrix
		float* c; // last m objective coefficients
		float z0; // 1st objective coefficient
		// change the vector b to +1
		void setB();
		// swap A b c with their dual according to the dict
		void dualize(const Dictionary& dict, bool isObjChng = true, bool isAllocated = true);
		// dualize by self
		void dualize();


};


#endif /* DICTIONARY_HPP_ */
