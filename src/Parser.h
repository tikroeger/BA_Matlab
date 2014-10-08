#ifndef _Parser_H__
#define _Parser_H__

//#include <mex.h>
#include <string>

//#include <ceres/ceres.h>
//#include <Eigen/Core>

#include "ceresBA_datastruct.h"

#define nullptr NULL

template <class T1, class T2> 
class Parser
{
	public:
    Parser() {};  // do-nothing constructor
		Parser(string);  // Load from file

		Parser(BAdata<T1,T2>*);   // Set pre-initialized datastruct directly
		~Parser(){ delete(da); da=nullptr; };
		BAdata<T1,T2>* getDataPtr() {return da;}; // return pointer to datastruct, allow ownership of data to other objects

		void WriteToFile(string); // Write to file
		void NormalizeData(); // Normalize data to zero median and unit variance
		void UnnormalizeData(); // Revert normalization
    
	protected:
		bool CheckErrors(); // check data for inconsistencies
    
		
		
		BAdata<T1,T2> *da;
};
		
		
#include "Parser.hxx"




#endif // _Parser_H__
