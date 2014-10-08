

#include"Parser.h"
#include"CeresBA.h"


using namespace std;


int main(int argc, char **argv)
{
  Parser<uint64_t, double> parser(argv[1]);
  
  parser.NormalizeData();

  CeresBA<uint64_t, double> ceresba(parser.getDataPtr());

  ceresba.CallSolver();

  parser.UnnormalizeData();

  parser.WriteToFile(argv[2]);

  return 0;
}
