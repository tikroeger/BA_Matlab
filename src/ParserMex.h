#ifndef _ParserMex_H__
#define _ParserMex_H__

#include <mex.h>
#include <string>

//#include <ceres/ceres.h>
//#include <Eigen/Core>

#include "Parser.h"


#include "ceresBA_datastruct.h"

#define nullptr NULL

template <class T1, class T2> 
class ParserMex : public Parser<T1,T2>
{
  public:
    ParserMex(int, const mxArray **, mxClassID, mxClassID);   // Parse mex input  arguments
    void WriteToMex(int*, mxArray **);                        // Write to mex output
    
  private:
    mxClassID T1mexclass;
    mxClassID T2mexclass;
};
    
    
#include "ParserMex.hxx"


#endif // _ParserMex__
