#include<iostream>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm> 
#include <iterator>
#include <sstream>
#include <fstream>
#include <iomanip>    
#include <numeric>
#include <valarray>



#include "CeresBA.h"
#include "ceresBA_datastruct.h"
#include "resutil.h"
#include "residual.h"

#include <ceres/ceres.h>
#include <Eigen/Core>

#define nullptr NULL


using namespace std;

template <class T1, class T2>
CeresBA<T1, T2>::CeresBA(BAdata<T1,T2>* da_in){
  da = da_in;
};


template <class T1, class T2>
void CeresBA<T1, T2>::EvalResidualTest() 
{
  T2 res[2];
  T1 pointno = 0;
  for (T1 pointno_reprojno=7 ;  pointno_reprojno  < 10;  ++pointno_reprojno)
  {

      cout <<  "*** Point no: " << pointno_reprojno << ": " << endl;
      //T1 pointno_reprojno = 0;
      
      double camweight, ptweight;  
      int64_t noview = (int64_t) da->PtReprojView[pointno][pointno_reprojno];
      int64_t noview_tr, noview_ori;

      int64_t camno = 0;  
      while (noview - (int64_t)da->Cams_noviews[camno]  > 0 )
      {
        noview -= da->Cams_noviews[camno];
        ++camno;
      }
      noview --; // indexing starts at 0 in C++ 
      
      if (da->Cams_fixOrientation[camno].size()==1) // when only one set of parameters present
        camweight = da->Cams_Weight[camno][0];
      else
        camweight = da->Cams_Weight[camno][noview];
      
      if (da->PtWeight.size()==1) // when only one set of parameters present
        ptweight = da->PtWeight[0];
      else
        ptweight = da->PtWeight[pointno];
        
      vector<bool> paramoptions(6);
      if (da->PtFixPosition.size()==1) // when only one set of parameters present
        paramoptions[0] = da->PtFixPosition[0];
      else
        paramoptions[0] = da->PtFixPosition[pointno]; 
      
      paramoptions[1] = da->Cams_fixInternals[camno][0]; // fc
      paramoptions[2] = da->Cams_fixInternals[camno][1]; // cc
      paramoptions[3] = da->Cams_fixInternals[camno][2]; // kc
      
      if (da->Cams_fixOrientation[camno].size()==1) // when only one set of parameters present
        noview_ori = 0;
      else
        noview_ori = noview;    

      paramoptions[4] = da->Cams_fixOrientation[camno][noview_ori]; 

      if (da->Cams_fixTranslation[camno].size()==1) // when only one set of parameters present
        noview_tr = 0;
      else
        noview_tr = noview;    
      
      paramoptions[5] = da->Cams_fixTranslation[camno][noview_tr]; 

      residual::RepError reperror(da->PtReprojPos[pointno][pointno_reprojno], 
                                  da->Cams_view_orien[camno][noview_ori], da->Cams_view_trans[camno][noview_tr],
                                  da->Pt3d[pointno], 
                                  da->Cams_fc[camno], da->Cams_cc[camno], da->Cams_kc[camno],                                   
                                  ptweight, camweight,
                                  paramoptions);
      reperror.operator()<T2>(&(da->Cams_view_orien[camno][noview_ori][0]), 
                                &(da->Cams_view_trans[camno][noview_tr][0]),
                                &(da->Pt3d[pointno][0]), 
                                &(da->Cams_fc[camno][0]), 
                                &(da->Cams_cc[camno][0]),
                                &(da->Cams_kc[camno][0]),
                                &(res[0]));                
      
      cout <<  "TEST, Point residual (AutoDiff): " << res[0] << " " << res[1] << endl;        

 
  }
}

template <class T1, class T2>
void CeresBA<T1, T2>::SetSolverOptions(ceres::Solver::Options* options)
{

  options->max_num_iterations = uint32_t(da->BAopt_maxiter);
  options->minimizer_progress_to_stdout = true;
  options->num_threads = uint32_t(da->BAopt_nocores);
  options->eta = double(1e-2);
  options->max_solver_time_in_seconds = double(da->BAopt_timeout_sec);

  options->use_nonmonotonic_steps = 0;

  ceres::StringToLinearSolverType("sparse_schur", &options->linear_solver_type); //dense_schur, iterative_schur, sparse_normal_cholesky, dense_qr, dense_normal_cholesky and cgnr
  ceres::StringToPreconditionerType("jacobi", &options->preconditioner_type); // identity, jacobi, schur_jacobi, cluster_jacobi, cluster_tridiagonal
  //ceres::StringToSparseLinearAlgebraLibraryType("suite_sparse", &options->sparse_linear_algebra_library_type); // cx_sparse
  ceres::StringToSparseLinearAlgebraLibraryType("cx_sparse", &options->sparse_linear_algebra_library_type);      // suitesparse
  options->num_linear_solver_threads = uint32_t(da->BAopt_nocores);
  
  options->solver_log = "";
  options->gradient_tolerance = 1e-16;
  options->function_tolerance = 1e-16;

}


template <class T1, class T2>
void CeresBA<T1, T2>::AddLineVPResiduals(ceres::Problem* problem)
{ 
  for (T1 camno=0 ;  camno < da->Cams_noviews.size()  ;  ++camno)
  //for (T1 camno=0 ;  camno < 1  ;  ++camno)    
  {
    if (da->Cams_VP_ID[camno].size() != 0) // VP Constraints are used 
    {
      for (T1 noview=0 ;  noview < da->Cams_noviews[camno]  ;  ++noview)
      //for (T1 noview=0 ;  noview < 1  ;  ++noview)        
      {
        for (T1 vpno = 0; vpno < da->Cams_VP_ID[camno][noview].size(); ++vpno)
        //for (T1 vpno = 0; vpno < 1; ++vpno)
        {
          T1 nolines = da->Cams_VPlines[camno][noview][vpno].size()/4;
          T1 vpid = da->Cams_VP_ID[camno][noview][vpno] - 1;
          T1 noview_ori;
          
//cout << "vpid : " <<  vpid << " nolines " << nolines << endl;
          double camweight, vpweight; 
          if (da->Cams_fixOrientation[camno].size()==1) // when only one value was set for all views in camera
            camweight = da->Cams_Weight[camno][0];
          else
            camweight = da->Cams_Weight[camno][noview];
          
          vpweight = da->VPConstr_VPw[vpid];
                
          vector<bool> paramoptions(6);
            
          paramoptions[0] = da->VPConstr_fixVP[vpid];
          
          paramoptions[1] = da->Cams_fixInternals[camno][0]; // fc
          paramoptions[2] = da->Cams_fixInternals[camno][1]; // cc
          paramoptions[3] = da->Cams_fixInternals[camno][2]; // kc

          if (da->Cams_fixOrientation[camno].size()==1) // when only one vector was set for all views in camera, i.e. shared orientation parameters
            noview_ori = 0;
          else
            noview_ori = noview;    
          
          paramoptions[4] = da->Cams_fixOrientation[camno][noview_ori]; 
          
          paramoptions[5] = 1;      // translation is irrelevant for VP fitting, set to constant. This is necessary for branching in OptionSwitchFixParameterAndAddResidual_PointCamera() to supply the correct functor

          if ( !(paramoptions[0] & paramoptions[1] & paramoptions[2] & paramoptions[3] & paramoptions[4])  )
          {
            T2* ptr1=NULL; T2* ptr3=NULL; T2* ptr4=NULL; T2* ptr5=NULL; T2* ptr6=NULL;

            ptr1 = &(da->Cams_view_orien[camno][noview_ori][0]);
            //ptr2 = NULL;
            ptr3 = &(da->VPConstr_VP[vpid][0]);
            ptr4 = &(da->Cams_fc[camno][0]);
            ptr5 = &(da->Cams_cc[camno][0]);
            ptr6 = &(da->Cams_kc[camno][0]);
            
            for (T1 nl = 0; nl < nolines; ++nl)
            //for (T1 nl = 0; nl < 1; ++nl)
            {
              int fcsize=da->Cams_fc[camno].size();
              int ccsize=da->Cams_cc[camno].size();
              int kcsize=da->Cams_kc[camno].size();
              
              vector<T2> midpoint(2), endpoint_A(2), endpoint_B(2);
              
              
              endpoint_A[0] = da->Cams_VPlines[camno][noview][vpno][0*nolines + nl];
              endpoint_A[1] = da->Cams_VPlines[camno][noview][vpno][1*nolines + nl];

              endpoint_B[0] = da->Cams_VPlines[camno][noview][vpno][2*nolines + nl];
              endpoint_B[1] = da->Cams_VPlines[camno][noview][vpno][3*nolines + nl];

              midpoint[0] = (endpoint_A[0] + endpoint_B[0]) / 2;
              midpoint[1] = (endpoint_A[1] + endpoint_B[1]) / 2;
              
              // Create functor
              residual::RepVPError* rpvp =  new residual::RepVPError(endpoint_A, midpoint, 
                                    da->Cams_view_orien[camno][noview_ori],
                                    da->VPConstr_VP[vpid], 
                                    da->Cams_fc[camno], da->Cams_cc[camno], da->Cams_kc[camno],                                   
                                    vpweight, camweight, paramoptions);
      
              // Choose loss function
              ceres::LossFunction* loss_function = InitializeLossFunction(da->BAopt_VPLoss); 
              
              // Add functor to Ceres
              OptionSwitchFixParameterAndAddResidual_VPCamera(problem, loss_function, fcsize, ccsize, kcsize,
                                                  ptr1, ptr3, ptr4, ptr5, ptr6, paramoptions, rpvp);
              
              
              // Same for other line endpoint
              
              rpvp =  new residual::RepVPError(endpoint_B, midpoint, 
                                    da->Cams_view_orien[camno][noview_ori],
                                    da->VPConstr_VP[vpid], 
                                    da->Cams_fc[camno], da->Cams_cc[camno], da->Cams_kc[camno],                                   
                                    vpweight, camweight, paramoptions);
              
              loss_function = InitializeLossFunction(da->BAopt_VPLoss); 

              OptionSwitchFixParameterAndAddResidual_VPCamera(problem, loss_function, fcsize, ccsize, kcsize,
                                                  ptr1, ptr3, ptr4, ptr5, ptr6, paramoptions, rpvp);
              
              
                          
            }
            
          }
            
        }
      }
    }
  }
}

template <class T1, class T2>
void CeresBA<T1, T2>::AddCameraDerivativeResiduals(ceres::Problem* problem)
{ 
  for (T1 camno=0 ;  camno < da->Cams_noviews.size()  ;  ++camno)
  {
    if (da->Cams_smootherM[camno].size()>0) // Smoother Matrix is set
    {
      for (T1 noview=0 ;  noview < da->Cams_noviews[camno]  ;  ++noview)
      {
          vector<T2> sm_val;
          vector<valarray<T2>> camTr;
          vector<valarray<T2>> camOr;
          sm_val.resize(0);

          T1 el_cent;
          T1 noviews_sm;
          if (da->Cams_smootherM[camno].size()>1) // switch if full smoother Matrix is available or same vector for all cameras
          {   // sm_val has length da->Cams_noviews[camno], element noview-1 is center element, remove zero padding
            noviews_sm = noview;
            el_cent = noview;  // element index of this noview camera in sm_val
          }
          else
          { // sm_val has odd length <= da->Cams_noviews[camno], middle element is element noview, remove zero padding
            noviews_sm  = 0;
            el_cent = T1(da->Cams_smootherM[camno][noviews_sm].size()/2);   // element index of this noview camera in sm_val
          }
          
          // Crop vector:  find left bound (lb) and right bound (rb) of elements which are non-zero, and within valid range of existing cameras
          int32_t lb=0;
          while ( lb < int32_t(el_cent) && (lb - int32_t(el_cent) + int32_t(noview) < 0 || abs(da->Cams_smootherM[camno][noviews_sm][lb])<1e-10)) // lb is smaller_equal then el_center and sm_val(lb) zero or lb invalid camera
            lb++;
          int32_t rb=da->Cams_smootherM[camno][noviews_sm].size()-1;
          while ( rb > int32_t(el_cent) && (rb - int32_t(el_cent) + int32_t(noview) >= int32_t(da->Cams_noviews[camno]) || abs(da->Cams_smootherM[camno][noviews_sm][rb])<1e-10)) // rb is larger_equal then el_center and sm_val(rb) zero or rb invalid camera
            rb--;

          T2 normabs=0, normfull=0;
          
          for (T1 el=lb ;  el <= T1(rb);  ++el)
          {
            sm_val.push_back(  da->Cams_smootherM[camno][noviews_sm][el]);
            normabs += abs(da->Cams_smootherM[camno][noviews_sm][el]);
            normfull += da->Cams_smootherM[camno][noviews_sm][el];
          }
          
          int nocamssm = (rb-lb+1); // number of cameras used in smoothing;
          vector<bool> paramoptions(2*nocamssm); // save information whether tr,or are fixed or free for all cameras used in smoothing

          camTr.resize(nocamssm);
          camOr.resize(nocamssm);
          T1 camcnt=0;
          T1 nofree=0; // count how many free and R's and t's we have
          for (T1 el= T1(lb) ;  el <= T1(rb);  ++el)
          {               
            T1 idx   = el - el_cent + noview; // index into camera array
            
            camTr[camcnt].resize(3);
            camOr[camcnt].resize(3);
            
            if (da->Cams_fixOrientation[camno].size()==1) // when only one set of parameters present
              paramoptions[camcnt         ] =   da->Cams_fixOrientation[camno][0]; 
            else
              paramoptions[camcnt         ] =   da->Cams_fixOrientation[camno][idx];
              
            if (da->Cams_fixTranslation[camno].size()==1) // when only one set of parameters present            
              paramoptions[camcnt+nocamssm] =   da->Cams_fixTranslation[camno][0];
            else
              paramoptions[camcnt+nocamssm] =   da->Cams_fixTranslation[camno][idx];
            
            nofree += !paramoptions[camcnt         ];
            nofree += !paramoptions[camcnt+nocamssm];
            
            if (da->Cams_fixOrientation[camno].size()==1) // when only one set of parameters present
            { 
              camOr[camcnt][0] = da->Cams_view_orien[camno][0][0]; 
              camOr[camcnt][1] = da->Cams_view_orien[camno][0][1]; 
              camOr[camcnt][2] = da->Cams_view_orien[camno][0][2]; 
            }
            else
            { 
              camOr[camcnt][0] = da->Cams_view_orien[camno][idx][0]; 
              camOr[camcnt][1] = da->Cams_view_orien[camno][idx][1]; 
              camOr[camcnt][2] = da->Cams_view_orien[camno][idx][2]; 
            } 
            
            if (da->Cams_fixTranslation[camno].size()==1) // when only one set of parameters present  
            { 
              camTr[camcnt][0] = da->Cams_view_trans[camno][0][0];  
              camTr[camcnt][1] = da->Cams_view_trans[camno][0][1]; 
              camTr[camcnt][2] = da->Cams_view_trans[camno][0][2]; 
            }
            else
            { 
              camTr[camcnt][0] = da->Cams_view_trans[camno][idx][0]; 
              camTr[camcnt][1] = da->Cams_view_trans[camno][idx][1]; 
              camTr[camcnt][2] = da->Cams_view_trans[camno][idx][2]; 
            }
            camcnt++;
          }   
          
          
          // create functor
          double camweight; 
          if (da->Cams_fixOrientation[camno].size()==1) // when only one value was set for all views in camera
            camweight = da->Cams_Weight[camno][0];
          else
            camweight = da->Cams_Weight[camno][noview];          
          residual::CamPathDerError* rp =  new residual::CamPathDerError(camOr, camTr, camweight, da->stdev3d, sm_val, paramoptions);
        
          // Add functor to Ceres
          ceres::DynamicAutoDiffCostFunction<residual::CamPathDerError, 4>* cost_function;
          cost_function = new ceres::DynamicAutoDiffCostFunction<residual::CamPathDerError, 4>(rp);  // 4: stride parameter, only affects speed,  -> see ceres modeling documentation
          
          vector<double*> vars(nofree);
          
          int32_t vars_idx=0;

          camcnt=0; 
          for (T1 el= T1(lb) ;  el <= T1(rb);  ++el)
          {               
            T1 idx   = el - el_cent + noview; // index into camera array

            if (paramoptions[camcnt         ]==0)  // If orientation is free variable, add this variable
            {
              if (da->Cams_fixOrientation[camno].size()==1) // when only one set of parameters present
                vars[vars_idx] = (&(da->Cams_view_orien[camno][0][0]));
              else
                vars[vars_idx] = (&(da->Cams_view_orien[camno][idx][0]));
              vars_idx++;
              cost_function->AddParameterBlock(3); // Add free orientation variable
            }
            camcnt++;
          }

          camcnt=0; 
          for (T1 el= T1(lb) ;  el <= T1(rb);  ++el)
          {               
            T1 idx   = el - el_cent + noview; // index into camera array
            
            if (paramoptions[4+camcnt+nocamssm]==0)  // If translation is free variable, add this variable
            {
              if (da->Cams_fixTranslation[camno].size()==1) // when only one set of parameters present  
                vars[vars_idx] = (&(da->Cams_view_trans[camno][0][0]));
              else
                vars[vars_idx] = (&(da->Cams_view_trans[camno][idx][0]));
              vars_idx++;
              cost_function->AddParameterBlock(3); // Add free translation variable
            }
            camcnt++; 
          }

          cost_function->SetNumResiduals(3);
          
          ceres::LossFunction* loss_function = InitializeLossFunction(da->BAopt_DerivLoss);               
          problem->AddResidualBlock(cost_function, loss_function, vars);
      }
    }
  }
}

  
template <class T1, class T2>
void CeresBA<T1, T2>::AddPointCameraResiduals(ceres::Problem* problem)
{ 
    for (T1 pointno=0 ;  pointno < da->Pt3d.size()  ;  ++pointno) //  
    {
      for (T1 pointno_reprojno=0 ;  pointno_reprojno <  da->PtReprojView[pointno].size()  ;  ++pointno_reprojno) 
      {
        int64_t noview = (T1) da->PtReprojView[pointno][pointno_reprojno];
        int64_t noview_tr, noview_ori;
                
        // Get camera number and view number in camera
        int64_t camno = 0;  
        while (noview - (int64_t)da->Cams_noviews[camno]  > 0 )
        {
          noview -= da->Cams_noviews[camno];
          ++camno;
        }
        noview --; // indexing starts at 0 in C++ 

        double camweight, ptweight; 
        if (da->Cams_fixOrientation[camno].size()==1) // when only one value was set for all views in camera
          camweight = da->Cams_Weight[camno][0];
        else
          camweight = da->Cams_Weight[camno][noview];
        
        if (da->PtWeight.size()==1) // when only one value was set for all points
          ptweight = da->PtWeight[0];
        else
          ptweight = da->PtWeight[pointno];
               
        vector<bool> paramoptions(6);
        if (da->PtFixPosition.size()==1) // when only one value was set for all points
          paramoptions[0] = da->PtFixPosition[0];
        else
          paramoptions[0] = da->PtFixPosition[pointno]; 
        
        paramoptions[1] = da->Cams_fixInternals[camno][0]; // fc
        paramoptions[2] = da->Cams_fixInternals[camno][1]; // cc
        paramoptions[3] = da->Cams_fixInternals[camno][2]; // kc

        if (da->Cams_fixOrientation[camno].size()==1) // when only one vector was set for all views in camera, i.e. shared orientation parameters
          noview_ori = 0;
        else
          noview_ori = noview;    
        paramoptions[4] = da->Cams_fixOrientation[camno][noview_ori]; 

        if (da->Cams_fixTranslation[camno].size()==1) // when only one vector was set for all views in camera, i.e. shared translation parameters
          noview_tr = 0;
        else
          noview_tr = noview;    
        paramoptions[5] = da->Cams_fixTranslation[camno][noview_tr];      
  
          
        // if at least one free parameter
        if ( !(paramoptions[0] & paramoptions[1] & paramoptions[2] & paramoptions[3] & paramoptions[4] & paramoptions[5])) 
        {
          // parameters+dims, 1:orientation (3d), 2:translation (3d), 3: 3dpoint (3d), 
          //                  4: foc (1/2-dim), 5: cc (0/2-dim), 6: kc (0/1/2/4/5-dim)
          // For parameters which are static, the dim value is set to 0, the passed pointer is NULL
          // The first parameter for the AutoDiffCostFunction is the dimensionality of the residual (here: 2)
          T2* ptr1=NULL; T2* ptr2=NULL; T2* ptr3=NULL; T2* ptr4=NULL; T2* ptr5=NULL; T2* ptr6=NULL;
          //int dim1=0; int dim2=0; int dim3=0; int dim4=0; int dim5=0; int dim6=0;
          //if (paramoptions[4]==0) {
          ptr1 = &(da->Cams_view_orien[camno][noview_ori][0]);
          ptr2 = &(da->Cams_view_trans[camno][noview_tr][0]);
          ptr3 = &(da->Pt3d[pointno][0]);
          ptr4 = &(da->Cams_fc[camno][0]);
          ptr5 = &(da->Cams_cc[camno][0]);
          ptr6 = &(da->Cams_kc[camno][0]);
                
          // Create functor
          residual::RepError* rp =  new residual::RepError(da->PtReprojPos[pointno][pointno_reprojno], 
                                da->Cams_view_orien[camno][noview_ori], da->Cams_view_trans[camno][noview_tr],
                                da->Pt3d[pointno], 
                                da->Cams_fc[camno], da->Cams_cc[camno], da->Cams_kc[camno],                                   
                                ptweight, camweight, paramoptions);
  

          // Choose loss function
          ceres::LossFunction* loss_function = InitializeLossFunction(da->BAopt_PTLoss); 

          int fcsize=da->Cams_fc[camno].size();
          int ccsize=da->Cams_cc[camno].size();
          int kcsize=da->Cams_kc[camno].size();
          
          // Add functor to Ceres
          OptionSwitchFixParameterAndAddResidual_PointCamera(problem, loss_function, fcsize, ccsize, kcsize,
                                              ptr1, ptr2, ptr3, ptr4, ptr5, ptr6, paramoptions, rp);                  
        }
      }
    }
}


template <class T1, class T2>
ceres::LossFunction* CeresBA<T1, T2>::InitializeLossFunction(vector<T2> lossparam)
{
  ceres::LossFunction* loss_function=NULL;
  if      (lossparam[0]==1) loss_function = new ceres::HuberLoss(lossparam[1]);
  else if (lossparam[0]==2) loss_function = new ceres::SoftLOneLoss(lossparam[1]);
  else if (lossparam[0]==3) loss_function = new ceres::CauchyLoss(lossparam[1]);
  else if (lossparam[0]==4) loss_function = new ceres::ArctanLoss(lossparam[1]);

  return loss_function;
}

template <class T1, class T2>
void CeresBA<T1, T2>::AddCameraCameraResiduals(ceres::Problem* problem)
{ 
    for (T1 reprono =0;  reprono <  da->Cams_reproj_view.size() ;  ++reprono) //  
    {
        int64_t noview = (int64_t) da->Cams_reproj_view[reprono][1];
        int64_t noview_repro = (int64_t) da->Cams_reproj_view[reprono][0];
        
        int64_t noview_tr, noview_ori;

        // Get camera number and view number in observing camera
        int64_t camno = 0;  
        while (noview - (int64_t)da->Cams_noviews[camno]  > 0 )
        {
          noview -= da->Cams_noviews[camno];
          ++camno;
        }
        noview --; // indexing starts at 0 in C++ 
        
        // Get camera number and view number in observed camera
        int64_t camno_repro = 0;  
        while (noview_repro - (int64_t)da->Cams_noviews[camno_repro]  > 0 )
        {
          noview_repro -= da->Cams_noviews[camno_repro];
          ++camno_repro;
        }
        noview_repro --; // indexing starts at 0 in C++ 
        // Generate options vector for residual functor, collect weights for view and point
        double camweight, camweight_repro;          
        vector<bool> paramoptions(6);
        if (da->Cams_fixTranslation[camno_repro].size()==1) // when only one value was set for all points
          paramoptions[0] = da->Cams_fixTranslation[camno_repro][0];
        else
          paramoptions[0] = da->Cams_fixTranslation[camno_repro][noview_repro]; 
        
        paramoptions[1] = da->Cams_fixInternals[camno][0]; // fc
        paramoptions[2] = da->Cams_fixInternals[camno][1]; // cc
        paramoptions[3] = da->Cams_fixInternals[camno][2]; // kc

        if (da->Cams_fixOrientation[camno].size()==1) // when only one value was set for all views in camera
          camweight = da->Cams_Weight[camno][0];
        else
          camweight = da->Cams_Weight[camno][noview];
        
        camweight_repro = da->Cams_reproj_weight[reprono];

        if (da->Cams_fixOrientation[camno].size()==1) // when only one vector was set for all views in camera, i.e. shared orientation parameters
          noview_ori = 0;
        else
          noview_ori = noview;    
        paramoptions[4] = da->Cams_fixOrientation[camno][noview_ori]; 

        if (da->Cams_fixTranslation[camno].size()==1) // when only one vector was set for all views in camera, i.e. shared translation parameters
          noview_tr = 0;
        else
          noview_tr = noview;    
        paramoptions[5] = da->Cams_fixTranslation[camno][noview_tr]; 
       
        // if at least one free parameter
        if ( !(paramoptions[0] & paramoptions[1] & paramoptions[2] & paramoptions[3] & paramoptions[4] & paramoptions[5])) 
        {

          // parameters+dims, 1:orientation (3d), 2:translation (3d), 3: 3dpoint (3d), 
          //                  4: foc (1/2-dim), 5: cc (0/2-dim), 6: kc (0/1/2/4/5-dim)
          // For parameters which are static, the dim value is set to 0, the passed pointer is NULL
          // The first parameter for the AutoDiffCostFunction is the dimensionality of the residual (here: 2)
          T2* ptr1=NULL; T2* ptr2=NULL; T2* ptr3=NULL; T2* ptr4=NULL; T2* ptr5=NULL; T2* ptr6=NULL;
        
          ptr1 = &(da->Cams_view_orien[camno][noview_ori][0]);
          ptr2 = &(da->Cams_view_trans[camno][noview_tr][0]);
          if (da->Cams_view_trans[camno_repro].size()==1)
              ptr3 = &(da->Cams_view_trans[camno_repro][0][0]);
          else
              ptr3 = &(da->Cams_view_trans[camno_repro][noview_repro][0]);
          ptr4 = &(da->Cams_fc[camno][0]);
          ptr5 = &(da->Cams_cc[camno][0]);
          ptr6 = &(da->Cams_kc[camno][0]);
          residual::RepError* rp;
          if (da->Cams_view_trans[camno_repro].size()==1)
          {
              rp =  new residual::RepError(da->Cams_reproj_pos[reprono], 
                          da->Cams_view_orien[camno][noview_ori], da->Cams_view_trans[camno][noview_tr],
                          da->Cams_view_trans[camno_repro][0], 
                          da->Cams_fc[camno], da->Cams_cc[camno], da->Cams_kc[camno],                                   
                          camweight_repro, camweight, paramoptions);
          }
          else
          {
              rp =  new residual::RepError(da->Cams_reproj_pos[reprono], 
                          da->Cams_view_orien[camno][noview_ori], da->Cams_view_trans[camno][noview_tr],
                          da->Cams_view_trans[camno_repro][noview_repro], 
                          da->Cams_fc[camno], da->Cams_cc[camno], da->Cams_kc[camno],                                   
                          camweight_repro, camweight, paramoptions);
          }

          // Choose loss function
          ceres::LossFunction* loss_function = InitializeLossFunction(da->BAopt_PTLoss);           

          int fcsize=da->Cams_fc[camno].size();
          int ccsize=da->Cams_cc[camno].size();
          int kcsize=da->Cams_kc[camno].size();
          
          OptionSwitchFixParameterAndAddResidual_PointCamera(problem, loss_function, fcsize, ccsize, kcsize,
                                              ptr1, ptr2, ptr3, ptr4, ptr5, ptr6, paramoptions, rp);
        }
    }
}


template <class T1, class T2>
void CeresBA<T1, T2>::OptionSwitchFixParameterAndAddResidual_PointCamera(ceres::Problem* problem, ceres::LossFunction* loss_function, int fcsize, int ccsize, int kcsize,
                                             T2* ptr1, T2* ptr2, T2* ptr3, T2* ptr4, T2* ptr5, T2* ptr6, vector<bool> paramoptions, residual::RepError* rp)
{
	/* Create a functor instantiation for every possible camera, point parameter combination. 
   *We could also use DynamicAutoDiffCostFunction, but declaring them statically should be faster in runtume,  (although much slower in compile-time)*/
ceres::CostFunction* cost_function;

if (paramoptions[4]==0) {             // ori is free parameter
    if (paramoptions[5]==0) {                      // tr is free parameter
        if (paramoptions[0]==0) {                    // pt is free parameter
            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr3); 
                  }
                }}
                  
        } else {  // pt is fixed

            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2); 
                  }
                }}
       }
    } else {             // tr is fixed
        if (paramoptions[0]==0) {                    // pt is free parameter
            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3); 
                  }
                }}
                  
        } else {  // pt is fixed

            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1); 
                  }
                }}
       }
    }
}
else            // ori is fixed parameter
{
    if (paramoptions[5]==0) {                      // tr is free parameter
        if (paramoptions[0]==0) {                    // pt is free parameter
            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr3); 
                  }
                }}
                  
        } else {  // pt is fixed

            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr2, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr2); 
                  }
                }}
       }
    } else {             // tr is fixed
        if (paramoptions[0]==0) {                    // pt is free parameter
            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3); 
                  }
                }}
                  
        } else {  // pt is fixed

            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr6);}
// *** remove the following option, 0 free parameters are not allowed                    
//                  else                                // kc is not modeled or static   
//                  { if (fcsize==2) // focal length fx,fy
//                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2>(rp);
//                    else                             // focal length fxy
//                      cost_function = new ceres::AutoDiffCostFunction<residual::RepError, 2>(rp);                      
//                    problem->AddResidualBlock(cost_function, loss_function); 
//                  } 
                }}
       }
    }
}
}

template <class T1, class T2>
void CeresBA<T1, T2>::OptionSwitchFixParameterAndAddResidual_VPCamera(ceres::Problem* problem, ceres::LossFunction* loss_function, int fcsize, int ccsize, int kcsize,
                                             T2* ptr1, T2* ptr3, T2* ptr4, T2* ptr5, T2* ptr6, vector<bool> paramoptions, residual::RepVPError* rp)
{
  /* Create a functor instantiation for every possible camera, point parameter combination. 
   *We could also use DynamicAutoDiffCostFunction, but declaring them statically should be faster in runtume,  (although much slower in compile-time)*/
ceres::CostFunction* cost_function;

if (paramoptions[4]==0) {             // ori is free parameter
        if (paramoptions[0]==0) {                    // vp is free parameter
            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr3); 
                  }
                }}
                  
        } else {  // vp is fixed

            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr1); 
                  }
                }}
       }
}
else            // ori is fixed parameter
{
        if (paramoptions[0]==0) {                    // vp is free parameter
            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr3, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 3>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr3); 
                  }
                }}
                  
        } else {  // vp is fixed


            if (paramoptions[1]==0) {                 // fc is free parameter
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr4, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr4, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr4, ptr6);}
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr4); 
                  }
                }
            }else{                       // fc is fixed
                if (paramoptions[2]==0 && ccsize > 0) // cc is specified and free parameter
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr5, ptr6);  }
                  else                                // kc is not modeled or static
                  { if (fcsize==2) // focal length fx,fy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2>(rp);
                    else                             // focal length fxy
                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2>(rp);                      
                    problem->AddResidualBlock(cost_function, loss_function, ptr5); 
                  }
                }
                else                                           // cc is not modelled or set as static
                { if (paramoptions[3]==0 && kcsize > 0)  // kc is specified and free parameter
                  { if (fcsize==2) // focal length fx,fy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 5>(rp); break;}
                    else                             // focal length fxy
                      switch (kcsize) {
                      case 1: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 1>(rp); break;
                      case 2: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 2>(rp); break;                        
                      case 4: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 4>(rp); break;                       
                      case 5: cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2, 5>(rp); break;}
                    problem->AddResidualBlock(cost_function, loss_function, ptr6);}
// *** remove the following option, 0 free parameters are not allowed                    
//                  else                                // kc is not modeled or static   
//                  { if (fcsize==2) // focal length fx,fy
//                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2>(rp);
//                    else                             // focal length fxy
//                      cost_function = new ceres::AutoDiffCostFunction<residual::RepVPError, 2>(rp);                      
//                    problem->AddResidualBlock(cost_function, loss_function); 
//                  } 
                }}
       }
    }
}


template <class T1, class T2>
void CeresBA<T1, T2>::AddPointPlaneResiduals(ceres::Problem* problem)
{ 
    for (T1 pointno=0 ;  pointno <  da->Pt3d.size() ;  ++pointno) //  
    {
        //cout << pointno ;
        double planeweight, ptweight;  
        
        // Get camera number and view number in camera
                             
        ptweight = da->PtWeight[pointno];
          
        for (T1 planeno_cnt=0 ;  planeno_cnt <  da->PtOnPlane[pointno].size() ;  ++planeno_cnt) //  
        {
            T1 planeno = da->PtOnPlane[pointno][planeno_cnt] - 1;
            //cout << " " << planeno  ;
            
            planeweight = da->PlaneConstr_PlanesW[planeno];

            vector<bool> paramoptions(2);
          //cout << da->PtFixPosition[pointno] << da->PlaneConstr_fixPlane[planeno] << endl;
            
            paramoptions[0] = da->PtFixPosition[pointno]; // fixPoint
            paramoptions[1] = da->PlaneConstr_fixPlane[planeno]; // fixPlane
            
            // if at least one free parameter
            if ( !(paramoptions[0] & paramoptions[1])) 
            {
              residual::PlaneError* plerr = new  residual::PlaneError(da->PlaneConstr_Planes[planeno], da->Pt3d[pointno],
                                                  planeweight, ptweight, paramoptions);
                                                  //camweight*stdev3d, ptweight*stdev3d);
              
              ceres::LossFunction* loss_function = InitializeLossFunction(da->BAopt_PlaneLoss);                         
            
              ceres::CostFunction* cost_function;
            
              T2 *ptr1 = &(da->PlaneConstr_Planes[planeno][0]);
              T2 *ptr2 = &(da->Pt3d[pointno][0]);
              
              
              if (paramoptions[0]==0 && paramoptions[1]==0) // point and plane are free
              {
                cost_function = new ceres::AutoDiffCostFunction<residual::PlaneError, 1, 3, 3>(plerr);
                problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2);
              }
              else
              {
                cost_function = new ceres::AutoDiffCostFunction<residual::PlaneError, 1, 3>(plerr);
                if (paramoptions[0]==0) // plane is free, point is fixed
                  problem->AddResidualBlock(cost_function, loss_function, ptr1);
                else  // point is free, plane is fixed
                  problem->AddResidualBlock(cost_function, loss_function, ptr2);          
              }
            }
        
        }  
    }
}



template <class T1, class T2>
void CeresBA<T1, T2>::AddCameraPlaneResiduals(ceres::Problem* problem)
{ 
    
    for (T1 camno=0 ;  camno <  da->Cams_noviews.size() ;  ++camno) //  
    {
      for (T1 planeno_cnt=0 ;  planeno_cnt <  da->Cams_OnPlane[camno].size() ;  ++planeno_cnt) // over all planes that this camera lies on
      {
          T1 planeno = da->Cams_OnPlane[camno][planeno_cnt];
          if (planeno !=0)
          {            
            planeno--; // indices in c++ start from zero;
          
          
            T2 planeweight = da->PlaneConstr_PlanesW[planeno];
            
            for (T1 noview=0 ;  noview < da->Cams_noviews[camno] ;  ++noview) //  
            {
              T2 camweight;
              if (da->Cams_Weight[camno].size() == 1)
                camweight = da->Cams_Weight[camno][0];
              else
                camweight = da->Cams_Weight[camno][noview];
              
                
              vector<bool> paramoptions(2);
              
              if (da->Cams_fixTranslation[camno].size()==1 || da->Cams_view_trans[camno].size()==1)
                paramoptions[0] = da->Cams_fixTranslation[camno][0]; // fixCameraPos
              else
                paramoptions[0] = da->Cams_fixTranslation[camno][noview]; // fixCameraPos
                  
              paramoptions[1] = da->PlaneConstr_fixPlane[planeno]; // fixPlane

              // if at least one free parameter
              if ( !(paramoptions[0] & paramoptions[1])) 
              {
                T2 *ptr1 = &(da->PlaneConstr_Planes[planeno][0]);
                T2 *ptr2;

                residual::PlaneError* plerr;
                if (da->Cams_view_trans[camno].size()==1)
                {
                  plerr = new  residual::PlaneError(da->PlaneConstr_Planes[planeno], da->Cams_view_trans[camno][0],
                                                    planeweight, camweight, paramoptions);
                  ptr2 = &(da->Cams_view_trans[camno][0][0]);                                                
                }
                else
                {
                  plerr = new  residual::PlaneError(da->PlaneConstr_Planes[planeno], da->Cams_view_trans[camno][noview],
                                                    planeweight, camweight, paramoptions);
                  ptr2 = &(da->Cams_view_trans[camno][noview][0]);                                                                                                
                }
                
                ceres::LossFunction* loss_function = InitializeLossFunction(da->BAopt_PlaneLoss);                     
              
                ceres::CostFunction* cost_function;
                
                if (paramoptions[0]==0 && paramoptions[1]==0)
                {
                  cost_function = new ceres::AutoDiffCostFunction<residual::PlaneError, 1, 3, 3>(plerr);
                  problem->AddResidualBlock(cost_function, loss_function, ptr1, ptr2);
                }
                else
                {
                  cost_function = new ceres::AutoDiffCostFunction<residual::PlaneError, 1, 3>(plerr);
                  if (paramoptions[0]==0)
                    problem->AddResidualBlock(cost_function, loss_function, ptr1);
                  else
                    problem->AddResidualBlock(cost_function, loss_function, ptr2);          
                }
              }
              
            }
          }
      }
    }
}


template <class T1, class T2>
void CeresBA<T1, T2>::CallSolver()  
{
  ceres::Problem problem; 
  
  AddPointCameraResiduals(&problem);
  AddCameraDerivativeResiduals(&problem);
  AddPointPlaneResiduals(&problem);
  AddCameraCameraResiduals(&problem);
  AddCameraPlaneResiduals(&problem);
  AddLineVPResiduals(&problem);

  //EvalResidualTest();
  
  ceres::Solver::Options options;
  SetSolverOptions(&options);
  ceres::Solver::Summary summary;
  Solve(options, &problem, &summary);

  std::cout << summary.FullReport() << "\n";

  da->Residual=summary.final_cost;

}







