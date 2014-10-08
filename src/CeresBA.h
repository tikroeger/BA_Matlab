#ifndef _CeresBA_H__
#define _CeresBA_H__

#include <ceres/ceres.h>
#include <Eigen/Core>


#include "ceresBA_datastruct.h"
#include "residual.h"

template <class T1, class T2> 
class CeresBA
{
	public:
		CeresBA(BAdata<T1,T2>*);

    void EvalResidualTest();  // test residual functors
    void EvalResidualTest_SM(); // test smoothing residual functors
    void CallSolver();
    
  private:
    void SetSolverOptions(ceres::Solver::Options*);
    void AddPointCameraResiduals(ceres::Problem*);
    void AddCameraDerivativeResiduals(ceres::Problem*);
    void AddPointPlaneResiduals(ceres::Problem*);
    void AddCameraCameraResiduals(ceres::Problem*);
    void AddCameraPlaneResiduals(ceres::Problem*);
    void AddLineVPResiduals(ceres::Problem*);
    
    void OptionSwitchFixParameterAndAddResidual_PointCamera(ceres::Problem*, ceres::LossFunction*, int, int, int,
                                                T2*, T2* , T2* , T2* , T2* , T2* , vector<bool>, residual::RepError*);
    
    
    void OptionSwitchFixParameterAndAddResidual_VPCamera(ceres::Problem*, ceres::LossFunction*, int, int, int,
                                                T2*, T2* , T2* , T2* , T2* , vector<bool>, residual::RepVPError*);    
    

    
    ceres::LossFunction* InitializeLossFunction(vector<T2>);
    
    BAdata<T1,T2>* da;

};
		
		
#include "CeresBA.hxx"




#endif // _CeresBA_H__
