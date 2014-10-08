#ifndef _BAdjust_H__
#define _BAdjust_H__

#include <mex.h>
#include <string>

#include <ceres/ceres.h>
#include <Eigen/Core>

#include "ceresBA_datastruct.h"

#define nullptr NULL

template <class T1, class T2> 
class BAdjust
{
	public:
		BAdjust(string);  // Load from file
		BAdjust(int, const mxArray **, mxClassID, mxClassID);   // Parse mex input  arguments
		BAdjust(BAdata<T1,T2>*);   // Set pre-initialized datastruct directly
		~BAdjust(){ delete(da); da=nullptr; };
		BAdata<T1,T2>* getDataPtr() {return da;}; // return pointer to datastruct, allow ownership of data to other objects
		void WriteToMex(int*, mxArray **); // Write to mex output
		void WriteToFile(string); // Write to file
		void NormalizeData(); // Normalize data to zero median and unit variance
		void UnnormalizeData(); // Revert normalization
		void EvalResidual();
		void CallSolver();
    
	private:
		void CheckErrors(); // check data for inconsistencies
    
		//vector<T2> RotatePoint(vector<T2>, vector<T2>); // Put in extra file
    //void SetSolverOptions(ceres::Solver::Options*);
		
		
		BAdata<T1,T2> *da;
		mxClassID T1mexclass;
		mxClassID T2mexclass;
		valarray<T2> median3d;
		T2 stdev3d;
};
		
		
#include "BAdjust_parser.hxx"
#include "BAdjust_ceresBA.hxx"




#endif // _BAdjust_H__
