#include "mex.h"
#include"ParserMex.h"
#include"CeresBA.h"
#include"ceresBAmex.h"

void __mexFunction__(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	// Verify Matlab T1, T2 datatypes
	long int fieldnum = mxGetFieldNumber(prhs[0], "formatversion");  
	mxArray* tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
	mxClassID T2mexclass = mxGetClassID(tmp);
	fieldnum = mxGetFieldNumber(prhs[0], "BAopt");
	tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
	fieldnum = mxGetFieldNumber(tmp, "cores");
	mxArray *tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
	mxClassID T1mexclass = mxGetClassID(tmp2);
	
	
	if (T1mexclass == mxUINT64_CLASS &&  T2mexclass == mxDOUBLE_CLASS)
	{
	    ParserMex<uint64_t, double> parser(nrhs, prhs, T1mexclass, T2mexclass);

      parser.NormalizeData();
	    
      CeresBA<uint64_t, double> ceresba(parser.getDataPtr());
      
      ceresba.CallSolver();
            
      parser.UnnormalizeData();
      
      parser.WriteToMex(&nlhs, plhs);
	}
}



void __at_exit__()
{}
