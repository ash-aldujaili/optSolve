/*
 * Dictionary.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: Abdullah Al-Dujaili
 */

#include"Dictionary.hpp"
#include<fstream>
#include<iostream>
#include<string>
#include<cmath>
#include<cassert>
#include<cfloat>
// constructor
Dictionary::Dictionary(const std::string fileName)
{
	std::ifstream inputFile(fileName.c_str());
	// updating info:
	inputFile >> m;
	inputFile >> n;
	// updating basic indices
	basicVar.reserve(m);
	for(int i = 0; i < m; i++)
	{
		int val;
		inputFile >> val;
		basicVar.push_back(val);
	}
	// update nonbasic indices:
	nBasicVar.reserve(n);
	for(int i = 0; i < n; i++)
	{
		int val;
		inputFile >> val;
		nBasicVar.push_back(val);
	}
	// b matrix:
	b = new float[m];
	for(int i = 0; i < m; i++)
	{
		inputFile >> b[i];
	}
	// A matrix
	a = new float[m*n];
	for(int i = 0; i < m * n; i++)
	{
		inputFile >> a[i];
	}
	// set A row pointers
	aRow = new float*[m];
	for(int i = 0; i < m; i++)
		aRow[i]=a + n * i;
	// coefficient value:
	inputFile >> z0;
	c = new float[n];
	for(int i = 0; i < n; i++)
	{
		inputFile >> c[i];
	}
	// numPivots
	numPivots = -1;
}
// destructor
Dictionary::~Dictionary()
{
	delete [] a;
	delete [] aRow;
	delete [] b;
	delete [] c;
}
// copy constructor * currently it supports copying for dulaizing puroposse only
Dictionary::Dictionary(const Dictionary& dict)
{
	dualize(dict,true, false);
}

void Dictionary::dualize(const Dictionary& dict, bool isObjChng, bool isAllocated)
{
	m = dict.n;
	n = dict.m;
	numPivots = -1;
	// vectors
	if (isAllocated) // this is because this function is used as a copy constructor and as dualizing function
	{
		delete [] a;
		delete [] aRow;
		delete [] b;
	}
	a = new float[n*m];
	aRow = new float*[m];
	b = new float[m];
	// basic , nonbasic variables
	basicVar=dict.nBasicVar;//.reserve(m); 	// list of basic indices [m integers]
	nBasicVar=dict.basicVar;//.reserve(n); // list of non-basic indices
	/*
	for(int i = 0; i < n; i++)
			nBasicVar.push_back((dict.basicVar[i] > m)? dict.basicVar[i]-m : dict.basicVar[i] + n) ;
	for(int i = 0; i < m; i++)
			basicVar.push_back((dict.nBasicVar[i] > m)? dict.nBasicVar[i]-m : dict.nBasicVar[i] + n) ;
*/
	for(int i = 0; i < m; i++)
	{
		b[i] = -dict.c[i];
		aRow[i]= a + i * n;
	}

	for(int j = 0; j < m; j++)
		for(int i = 0; i < n; i++)
			aRow[j][i]=-dict.a[i * m + j];

	// dual C ?
	if (isObjChng)
	{
		if (isAllocated)
			delete [] c;
		c = new float[n];
		z0 = -dict.z0;
		for(int i = 0; i < n; i++)
			c[i] = -dict.b[i];
	}

}


// dualize on self:
void Dictionary::dualize()
{
	int tempN = n;
	n = m;
	m = tempN;
	numPivots = -1;
	// vectors
	float *aNew = new float[n*m];
	float **aRowNew = new float*[m];
	float *bNew = new float[m];
	// basic , nonbasic variables
	std::vector<int> tempVar =nBasicVar;//.reserve(m); 	// list of basic indices [m integers]
	nBasicVar=basicVar;//.reserve(n); // list of non-basic indices
	basicVar = tempVar;

	for(int i = 0; i < m; i++)
	{
		bNew[i] = -c[i];
		aRowNew[i]= aNew + i * n;
	}

	for(int j = 0; j < m; j++)
		for(int i = 0; i < n; i++)
			aRowNew[j][i]=-a[i * m + j];

	// dual C ?
	float *cNew = new float[n];
	z0 = -z0;
	for(int i = 0; i < n; i++)
		cNew[i] = -b[i];

	delete [] a;
	delete [] b;
	delete [] c;
	delete [] aRow;

	a = aNew;
	b = bNew;
	c = cNew;
	aRow = aRowNew;

}
// initialize the dictionary to get a feasible region:
bool Dictionary::toFeasibleDict(bool isVerbose, bool isObjChng)
{
	// create the dual dictionary:
	Dictionary dualDict(*this);
	// change the objective function to negative
	if (isObjChng)
		dualDict.setB();
	if (isVerbose)
	{
		std::cout<<"Dual dict\n";
		print(dualDict);
	}
	// pivot to a final objective value:
	if(dualDict.pivotToFinalDict())
	{
		//
		if (isVerbose)
		{
			std::cout<<"Final Dual DICT:\n";
			print(dualDict);
		}
		// before dualizing see what the actual object function nonbasic var are to substitute
		// the current nonbasic var into it
		std::vector<int> objVar = nBasicVar;
		float *newC = new float[n];
		for (int i = 0; i < n; i++)
			newC[i] = 0.0;
		// dual the dual feasible dict to get the feasible dict (with no change on the objective val c)
		dualize(dualDict, false);


		// rectify c according to the current nonbasic var
		// go through the first
		for (int i = 0; i < n; i++)
		{
			bool isNBasic = false;
			for (int k = 0;(!isNBasic && k < n); k++)
			{
				if (nBasicVar[k] == objVar[i])
				{
					isNBasic = true;
					newC[k] += c[i];
				}
			}
			if (!isNBasic)
				for(int j = 0; j < m; j++)
				{
					if (basicVar[j] == objVar[i])
					{
						for(int k = 0; k < n; k++)
							newC[k] += c[i] * aRow[j][k];
						z0 += c[i] * b[j];
					}
					else
						continue;
				}
		}
		// set the final c:
		for(int i = 0; i < n; i++)
			c[i] = newC[i];
		// free newC;
		// new primal:
		if (isVerbose)
		{
			std::cout<<"New Primal\n";
			print(*this);
		}
		delete [] newC;
		// write the
		return true;
	}
	else
		return false;

}
// set the b vector to 1
void Dictionary::setB()
{
	for(int i = 0; i < m; i++)
		b[i]=  +1.0;
}
// performs pivotint till it reach a final or unbounded dictionary
// returns true if the dictionary is final, false if unbounded
bool Dictionary::pivotToFinalDict(const bool isVerbose)
{
	int indicator;
	do
	{
		indicator = pivot(false);
		if (isVerbose) print(*this);
		numPivots++;
		if (indicator < 0)
		{
			if (isVerbose) std::cout<<"After "<<numPivots<<", the dictionary is unbounded !\n";
			return false;
		}
	}while(indicator == 0);

	return true;
}
/* performs a pivoting step and modifies the dictionary:
 * @return: 0 if the dictionary is not final, +1 if it is final, -1 if its unbounded
 */
int Dictionary::pivot(const bool isVerbose)
{
	// bounded
	bool isBounded = false;
	bool isFinal = true;
	// decide upon the entering variable (Bald's Rule)
	int enterVarId = m + n + 1;
	int enterId=0;
	for(int i = 0; i < n; i++)
	{
		if (c[i] > 1e-12) // a variabel is a possible entering var if it has positive coefficient
		{
			isFinal = false;
			if (nBasicVar[i] < enterVarId)
			{
				enterVarId = nBasicVar[i];
				enterId = i;
			}
		}
	}

	if (isFinal)
	{
		if (isVerbose) std::cout<<"Final Dictionary";
		return +1;
	}

	assert(enterVarId > 0);
	// decide the leaving variable:
	int leaveVarId = m + n + 1;
	int leaveId = 0;
	float val = FLT_MAX;
	//
	for(int j = 0; j < m; j++)
	{
		float tempVal = 0.0;
		bool skip = false;
		//if (b[j] <= FLT_EPSILON)
		//	tempVal = 0.0;
		if (-a[j * n + enterId] > 1e-5)
			tempVal = (-1.0 * b[j] / a[j * n + enterId]);
		else skip = true;
					//std::cout<<"leave var " << b[j] << " " << a[j*n+ enterId] << "tempval: " <<tempVal<<" \n";
		if (!skip)
			if ( val >= tempVal)
			{
				if (fabs(val - tempVal) <= 1e-5 )
				{
					//leaveVarId = (fabs(val -tempVal) <= FLT_EPSILON )? std::min(basicVar[j], leaveVarId) :basicVar[j]  ;
					if (leaveVarId > basicVar[j])
					{
						leaveVarId = basicVar[j];
						leaveId = j;
						isBounded = true;
						val = tempVal;
					}
				}
				else
				{
					leaveVarId = basicVar[j]  ;
					leaveId = j;
					isBounded = true;
					val = tempVal;
				}
			}
	}
	if (!isBounded)
	{
		if (isVerbose) std::cout<<"UNBOUNDED\n";
		return -1;
	}
	// calculate the objective function:
	val = z0 + val * c[enterId];


	// update the dictionary
	// 1. Update the row that corresponds to the leaving var
	float enterVarCof = -a[leaveId * n + enterId];
	b[leaveId] = b[leaveId] / enterVarCof;
	for(int i = 0; i < n; i++)
	{
		if (i == enterId)
			aRow[leaveId][i]= -1.0 / enterVarCof;
		else
			aRow[leaveId][i]= aRow[leaveId][i]/ enterVarCof;
	}
	// 2. update the rest of the rows
	for(int j = 0; j < m; j++)
	{
		if (j == leaveId) // the leaving var
		{
			continue;
		}
		enterVarCof = a[j * n + enterId];
		b[j] = b[j] + enterVarCof * b[leaveId];
		for(int i = 0; i < n; i++)
		{
			if(i == enterId)
				aRow[j][i] = enterVarCof * aRow[leaveId][i];
			else
				aRow[j][i] += enterVarCof * aRow[leaveId][i];
		}

	}
	// update z0
	enterVarCof = c[enterId];
	z0 = z0 + enterVarCof * b[leaveId];
	// update c
	for(int i = 0; i < n; i++)
	{
		if(i == enterId)
			c[i] = enterVarCof * aRow[leaveId][i];
		else
			c[i] += enterVarCof * aRow[leaveId][i];
	}


	// change the positions of the varialbes:
	basicVar[leaveId] = enterVarId;
	nBasicVar[enterId] = leaveVarId;
	// printing things:
	if (isVerbose)
	{
		std::cout.precision(5);
		std::cout << enterVarId << "\n" << leaveVarId << "\n" << val << "\n";
	}

	return 0;
}
// for viewing purposes:
void print(const Dictionary &dict)
{
	std::cout<<"************************************************************\n";
	// upper part
	for(int i = 0; i < dict.m; i++)
	{
		std::cout << "x" << dict.basicVar[i] << " : " << dict.b[i];
		for(int j = 0; j < dict.n; j++)
			std::cout << "|\t" << dict.aRow[i][j] << "*x" <<  dict.nBasicVar[j];
		std::cout << "\n";
	}
	// lower part:
	std::cout << "Z" << " : " << dict.z0;
	for(int i = 0; i < dict.n; i++)
		std::cout << "|\t" << dict.c[i] << "*x" <<  dict.nBasicVar[i];
	std::cout << "\n";
	std::cout<<"************************************************************\n";
}
// return UNBOUNDED, INFEASIBLE or the achieved objective value:
void Dictionary::solve()
{
	// check for feasibility:
	if(!this->toFeasibleDict())
	{
		std::cout<<"INFEASIBLE\n";
		return;
	}
	if (this->pivotToFinalDict(false))
	{
		std::cout.precision(5);
		std::cout<<this->getObjVal()<<"\n";
	}
	else
		std::cout<<"UNBOUNDED\n";

	return;
}
//
// return UNBOUNDED, INFEASIBLE or the achieved objective value with integer constraint
// using cutting plane method:
void Dictionary::solveILP()
{

	std::vector<int> nonIntVar;
	bool isFirstDict = true;

	// do a first iteration first
	// check for unbounded, or infeasibility
	if(!this->toFeasibleDict(false,isFirstDict))
	{
		std::cout<<"INFEASIBLE\n";
		return;
	}
	if (this->pivotToFinalDict(false))
	{
		// clear the vector integer
		nonIntVar.clear();
		// collect non integer variables indices
		for(int i = 0; i < m; i++)
			if(fabs(b[i]-roundf(b[i])) > 1e-3)
			{
				//std::cout<<fabs(b[i]-(long) b[i])<<"\n";
				nonIntVar.push_back(i);
			}
		// check if there is any
		if (nonIntVar.size() == 0)
		{
			std::cout.precision(5);
			std::cout<<this->getObjVal()<<"\n";
			return;
		}
	}
	else
	{
		std::cout<<"UNBOUNDED\n";
		return;
	}
	//print(*this);
	// cutting plane:
	do
	{
		// add another constraint (cutting planes)
		int mPlus = nonIntVar.size();
		// construct a new dictionary
		float *bNew = new float[m+mPlus];
		float *aNew = new float[(m+mPlus)*n];
		float **aRowNew = new float*[(m+mPlus)];

		// copy:
		for (int i = 0; i < m*n; i++)
			aNew[i]=a[i];
		for (int i = 0; i < m; i++)
			bNew[i]=b[i];
		for (int i = 0; i < m+mPlus; i++)
			aRowNew[i]= aNew + (i * n);
		// fill the new constraints:
		for (int i = 0; i < mPlus; i++)
		{
			// b part
			bNew[i + m]=-(bNew[nonIntVar[i]] - std::floor(bNew[nonIntVar[i]]));
			// new constraint basic variable
			basicVar.push_back(i + m + n + 1);
			// a matrix
			for (int j = 0; j < n; j++)
				aRowNew[i+m][j]= (-aRowNew[nonIntVar[i]][j] - std::floor(-aRowNew[nonIntVar[i]][j]));
		}

		// reassign the arrays:
		delete [] a;
		delete [] aRow;
		delete [] b;

		a = aNew;
		b = bNew;
		aRow = aRowNew;
		m = m + mPlus;
		//print(*this);


		// solve the new problem with dual since the primal is infeasible
		// check for feasibility:
		if(!this->toFeasibleDict())
		{
			std::cout<<"INFEASIBLE\n";
			return;
		}
		if (this->pivotToFinalDict(false))
		{
			//clean b
			for(int i = 0; i < m; i ++)
				if(fabs(b[i]) < 1e-5)
					b[i] = 0.f;
			//clean c
			for(int i = 0; i < n; i ++)
				if(fabs(c[i]) < 1e-5)
					c[i] = 0.f;
			//print(*this);
			// clear the vector integer
			nonIntVar.clear();
			// collect non integer variables indices
			for(int i = 0; i < m; i++)
				if(fabs(b[i]-roundf(b[i])) > 1e-2)
				{
					//std::cout<<fabs(b[i]-roundf(b[i]))<<" "<<b[i]<<" "<<(roundf(b[i]))<<"\n";
					nonIntVar.push_back(i);
				}
			// check if there is any
			if (nonIntVar.size() == 0)
			{
				std::cout.precision(5);
				std::cout<<this->getObjVal()<<"\n";
				return;
			}

			//std::cout.precision(5);
			//std::cout<<this->getObjVal()<<"\n";
		}
		else
			std::cout<<"UNBOUNDED\n";
		/*
		// dualize:
		dualize();
		//print(*this);
		// pivot to the final dict
		if (this->pivotToFinalDict(false))
		{
			//print(*this);
			// dualize to the final primal
			dualize();
			//print(*this);
			// clear the vector integer
			nonIntVar.clear();
			// collect non integer variables indices
			for(int i = 0; i < m; i++)
				if(fabs(b[i]-roundf(b[i])) > 1e-5)
				{
					//std::cout<<fabs(b[i]-roundf(b[i]))<<" "<<b[i]<<" "<<((int) b[i])<<"\n";
					nonIntVar.push_back(i);
				}
			// check if there is any
			if (nonIntVar.size() == 0)
			{
				std::cout.precision(5);
				std::cout<<this->getObjVal()<<"\n";
				return;
			}
		}
		else
		{
			std::cout<<"UNBOUNDED\n";
			return;
		}
		*/
	}while(true);

}
// test if the dictionary is final
bool Dictionary::isFinal() const
{
	bool isFinal = true;
	for(int i = 0; i < n; i++)
	{
		if (c[i] > 1e-12)
		{
			isFinal = false;
			break;
		}
	}

	return isFinal;
}


