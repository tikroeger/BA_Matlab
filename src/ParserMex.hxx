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


#include "ParserMex.h"
#include "ceresBA_datastruct.h"


template <class T1, class T2>
ParserMex<T1, T2>::ParserMex(int nrhs, const mxArray **prhs, mxClassID T1mexClass_in, mxClassID T2mexClass_in) {
  T1mexclass = T1mexClass_in;
  T2mexclass = T2mexClass_in;
  
  const char * mexfieldnames_root[] = {"formatversion", "BAopt", "planeconstraint", "VPconstraint", "cams", "cam_reproj", "points"};
  const char * mexfieldnames_BAopt[] = {"cores", "timeout_sec", "maxiter", "PTLoss", "PlaneLoss", "DerivLoss", "VPLoss"};
  const char * mexfieldnames_planeconstraints[] = {"plane", "planew", "fixPlane"};
  const char * mexfieldnames_VPconstraints[] = {"vp", "vpw", "fixVP"};
  const char * mexfieldnames_cams[] = {"noviews", "fc", "cc", "kc", "views_orien", "views_trans", "viewids", "fixInternals", "fixOrientation", "fixTranslation", "OnPlane", "camweight", "smootherM", "vplines"};
  const char * mexfieldnames_cam_reproj[] = {"view", "pos", "weight"};
  const char * mexfieldnames_points[] = {"pt3d", "OnPlane", "fixPosition", "pointw", "reproj_view", "reproj_pos"};
  
  this->da = new BAdata<T1,T2>();
  T1* dataT1;
  T2* dataT2;
  bool* dataBool;
  mwIndex Cellidx;
  T1 novps, nolines;
  
  
  long int fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[0]);  
  mxArray *tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  mxArray *tmp3, *tmp4, *tmp5;
  dataT2  = (T2*)mxGetData(tmp);
  this->da->Formatversion = dataT2[0];

  // Parse BAopt
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[1]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[0]);
  mxArray *tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT1  = (T1*)mxGetData(tmp2);
  this->da->BAopt_nocores = dataT1[0];
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT1  = (T1*) mxGetData(tmp2);
  this->da->BAopt_timeout_sec = dataT1[0];

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT1  = (T1*) mxGetData(tmp2);
  this->da->BAopt_maxiter = dataT1[0];
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[3]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  this->da->BAopt_PTLoss.resize(2);
  this->da->BAopt_PTLoss[0] = dataT2[0];
  this->da->BAopt_PTLoss[1] = dataT2[1];  
    
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[4]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  this->da->BAopt_PlaneLoss.resize(2);
  this->da->BAopt_PlaneLoss[0] = dataT2[0];
  this->da->BAopt_PlaneLoss[1] = dataT2[1];  

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[5]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  this->da->BAopt_DerivLoss.resize(2);
  this->da->BAopt_DerivLoss[0] = dataT2[0];
  this->da->BAopt_DerivLoss[1] = dataT2[1];  
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[6]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  this->da->BAopt_VPLoss.resize(2);
  this->da->BAopt_VPLoss[0] = dataT2[0];
  this->da->BAopt_VPLoss[1] = dataT2[1];
  
  
  // Parse plane constraints
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[2]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_planeconstraints[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  uint64_t numelem = mxGetNumberOfElements(tmp2);
  copy(dataT2, dataT2+numelem, back_inserter(this->da->PlaneConstr_PlanesW));

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_planeconstraints[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataBool  = (bool*) mxGetData(tmp2);
  this->da->PlaneConstr_fixPlane.resize(numelem);
  for (T1 i = 0; i<numelem; ++i)
    this->da->PlaneConstr_fixPlane[i] = dataBool[i]; 

  this->da->PlaneConstr_Planes.resize(numelem);
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_planeconstraints[0]);  
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);  
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<numelem; ++i)
  {
    this->da->PlaneConstr_Planes[i].resize(3);
    this->da->PlaneConstr_Planes[i][0] = dataT2[i*3];
    this->da->PlaneConstr_Planes[i][1] = dataT2[i*3+1];
    this->da->PlaneConstr_Planes[i][2] = dataT2[i*3+2];
  }

  // Parse VP constraints
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[3]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_VPconstraints[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  numelem = mxGetNumberOfElements(tmp2);
  copy(dataT2, dataT2+numelem, back_inserter(this->da->VPConstr_VPw));

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_VPconstraints[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataBool  = (bool*) mxGetData(tmp2);
  this->da->VPConstr_fixVP.resize(numelem);
  for (T1 i = 0; i<numelem; ++i)
    this->da->VPConstr_fixVP[i] = dataBool[i]; 

  this->da->VPConstr_VP.resize(numelem);
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_VPconstraints[0]);  
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);  
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<numelem; ++i)
  {
    this->da->VPConstr_VP[i].resize(3);
    this->da->VPConstr_VP[i][0] = dataT2[i*3];
    this->da->VPConstr_VP[i][1] = dataT2[i*3+1];
    this->da->VPConstr_VP[i][2] = dataT2[i*3+2];
  }
  

  
 
  // Parse Cameras
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[4]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum); 
  
  T1 nocams = mxGetNumberOfElements(tmp);
  this->da->Cams_fc.resize(nocams);
  this->da->Cams_cc.resize(nocams);
  this->da->Cams_kc.resize(nocams);
  this->da->Cams_fixInternals.resize(nocams);
  this->da->Cams_viewids.resize(nocams);
  this->da->Cams_view_orien.resize(nocams);
  this->da->Cams_view_trans.resize(nocams);
  this->da->Cams_fixOrientation.resize(nocams);
  this->da->Cams_fixTranslation.resize(nocams);
  this->da->Cams_OnPlane.resize(nocams);
  this->da->Cams_Weight.resize(nocams);
  this->da->Cams_smootherM.resize(nocams);
  this->da->Cams_VP_ID.resize(nocams);
  this->da->Cams_VPlines.resize(nocams);

  
  for (T1 i = 0; i<nocams; ++i)
  {
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[0]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT1  = (T1*) mxGetData(tmp2);
    this->da->Cams_noviews.push_back(dataT1[0]);

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[1]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(this->da->Cams_fc[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[2]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(this->da->Cams_cc[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[3]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(this->da->Cams_kc[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[4]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    this->da->Cams_view_orien[i].resize(numelem/3);
    for (T1 j = 0; j<this->da->Cams_view_orien[i].size(); ++j)
    {
      this->da->Cams_view_orien[i][j].resize(3);
      this->da->Cams_view_orien[i][j][0] = dataT2[j*3];
      this->da->Cams_view_orien[i][j][1] = dataT2[j*3+1];
      this->da->Cams_view_orien[i][j][2] = dataT2[j*3+2];
    }      


    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[5]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    this->da->Cams_view_trans[i].resize(numelem/3);
    for (T1 j = 0; j<this->da->Cams_view_trans[i].size(); ++j)
    {
      this->da->Cams_view_trans[i][j].resize(3);
      this->da->Cams_view_trans[i][j][0] = dataT2[j*3];
      this->da->Cams_view_trans[i][j][1] = dataT2[j*3+1];
      this->da->Cams_view_trans[i][j][2] = dataT2[j*3+2];
    }        

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[6]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT1  = (T1*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT1, dataT1+numelem, back_inserter(this->da->Cams_viewids[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[7]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    this->da->Cams_fixInternals[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      this->da->Cams_fixInternals[i][j] = dataBool[j]; 

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[8]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    this->da->Cams_fixOrientation[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      this->da->Cams_fixOrientation[i][j] = dataBool[j]; 

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[9]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    this->da->Cams_fixTranslation[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      this->da->Cams_fixTranslation[i][j] = dataBool[j]; 

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[10]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT1  = (T1*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT1, dataT1+numelem, back_inserter(this->da->Cams_OnPlane[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[11]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(this->da->Cams_Weight[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[12]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    if (numelem == 0) // no smoothing
      this->da->Cams_smootherM[i].resize(0);
    else if (numelem == this->da->Cams_noviews[i] * this->da->Cams_noviews[i]) // full nxn smoother matrix
    {
        this->da->Cams_smootherM[i].resize(this->da->Cams_noviews[i]);
        for (T1 j = 0; j<this->da->Cams_noviews[i]; ++j)
        {
          this->da->Cams_smootherM[i][j].resize(this->da->Cams_noviews[i]);
        }
        for (T1 j = 0; j<this->da->Cams_noviews[i]; ++j)
          for (T1 k = 0; k<this->da->Cams_noviews[i]; ++k)
            this->da->Cams_smootherM[i][k][j] = dataT2[j + k*this->da->Cams_noviews[i]];
    }
    else // same smoothing vector for all data points
    {
      this->da->Cams_smootherM[i].resize(1);
      this->da->Cams_smootherM[i][0].resize(numelem);
      for (T1 j = 0; j<numelem; ++j)
        this->da->Cams_smootherM[i][0][j] = dataT2[j];
    }

    // parse VP constraints
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[13]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    
    if (mxGetNumberOfElements(tmp2)!=this->da->Cams_noviews[i] )
    {
        this->da->Cams_VP_ID[i].resize(0);
        this->da->Cams_VPlines[i].resize(0);
    }      
    else
    {
      this->da->Cams_VP_ID[i].resize(this->da->Cams_noviews[i]);
      this->da->Cams_VPlines[i].resize(this->da->Cams_noviews[i]);
      for (T1 j = 0; j<this->da->Cams_noviews[i]; ++j)
      {
    
        fieldnum = mxGetFieldNumber(tmp2, "l");
        tmp3 = mxGetFieldByNumber(tmp2, j, fieldnum);
          
        novps = mxGetM(tmp3);
        this->da->Cams_VP_ID[i][j].resize(novps);
        this->da->Cams_VPlines[i][j].resize(novps);
        
        for (T1 k = 0; k < novps; ++k)
        {
          
          mwSize subs[2];
          subs[0]=k; subs[1]=0;
          Cellidx = mxCalcSingleSubscript(tmp3, 2, subs);
          tmp4 = mxGetCell(tmp3, Cellidx);
          dataT1 = (T1*) mxGetData(tmp4);
          
          this->da->Cams_VP_ID[i][j][k] = dataT1[0] ;
          
          subs[0]=k; subs[1]=1;
          Cellidx = mxCalcSingleSubscript(tmp3, 2, subs);
          tmp4 = mxGetCell(tmp3, Cellidx);
          dataT2 = (T2*) mxGetData(tmp4);
          
          nolines = (T1) mxGetM(tmp4);
          this->da->Cams_VPlines[i][j][k].resize(4*nolines);
          
          for (T1 ki = 0; ki < 4*nolines; ++ki)
            this->da->Cams_VPlines[i][j][k][ki] = dataT2[ki];

        }
          
      }      
    }

      

  }

  // Parse mutual camera visibility
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[5]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum); 
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cam_reproj[0]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  T1 nomutualvis = mxGetNumberOfElements(tmp2) / 2;
  
  this->da->Cams_reproj_view.resize(nomutualvis);
  this->da->Cams_reproj_pos.resize(nomutualvis);

  dataT1  = (T1*) mxGetPr(tmp2);
  for (T1 i = 0; i<nomutualvis; ++i)
  {
    this->da->Cams_reproj_view[i].resize(2);
    this->da->Cams_reproj_view[i][0] = dataT1[i+0];
    this->da->Cams_reproj_view[i][1] = dataT1[i+nomutualvis];
  }

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cam_reproj[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<nomutualvis; ++i)
  {
    this->da->Cams_reproj_pos[i].resize(2);
    this->da->Cams_reproj_pos[i][0] = dataT2[i+0]; 
    this->da->Cams_reproj_pos[i][1] = dataT2[i+nomutualvis];
  }
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cam_reproj[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  copy(dataT2, dataT2+nomutualvis, back_inserter(this->da->Cams_reproj_weight));

  // Parse points
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[6]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum); 
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[0]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  T1 nopts = mxGetNumberOfElements(tmp2) / 3;
  this->da->Pt3d.resize(nopts);
  this->da->PtReprojView.resize(nopts);
  this->da->PtReprojPos.resize(nopts);
  this->da->PtOnPlane.resize(nopts);  
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<nopts ; ++i) 
  {
    this->da->Pt3d[i].resize(3);
    this->da->Pt3d[i][0] = dataT2[i*3];
    this->da->Pt3d[i][1] = dataT2[i*3+1];
    this->da->Pt3d[i][2] = dataT2[i*3+2];
  }        

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataBool = (bool*) mxGetPr(tmp2);
  this->da->PtFixPosition.resize(nopts);
  for (T1 j = 0; j<nopts; ++j)
    this->da->PtFixPosition[j] = dataBool[j];

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[3]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2 = (T2*) mxGetPr(tmp2);
  copy(dataT2, dataT2+nopts, back_inserter(this->da->PtWeight));

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[4]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[5]);
  tmp3 = mxGetFieldByNumber(tmp, 0, fieldnum);
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[1]);
  tmp5 = mxGetFieldByNumber(tmp, 0, fieldnum);
  
  for (T1 i = 0; i<nopts ; ++i) 
  {

    tmp4 = mxGetCell(tmp2, i);
    dataT1 = (T1*) mxGetPr(tmp4);
    numelem = mxGetNumberOfElements(tmp4);
    copy(dataT1, dataT1+numelem, back_inserter(this->da->PtReprojView[i]));

    tmp4 = mxGetCell(tmp3, i);
    dataT2 = (T2*) mxGetPr(tmp4);
    numelem = mxGetNumberOfElements(tmp4);
    this->da->PtReprojPos[i].resize(numelem/2);
    for (T1 j = 0; j<this->da->PtReprojPos[i].size(); ++j) 
      copy(dataT2+j*2, dataT2+(j+1)*2, back_inserter(this->da->PtReprojPos[i][j]));
  
    tmp4 = mxGetCell(tmp5, i);
    if (tmp4!=nullptr) 
    {
      numelem = mxGetNumberOfElements(tmp4);
      dataT1 = (T1*) mxGetPr(tmp4);
      copy(dataT1, dataT1+numelem, back_inserter(this->da->PtOnPlane[i]));
    }
  }
  
  this->da->Residual=-1;
  
  this->CheckErrors();
}




template <class T1, class T2>
void ParserMex<T1, T2>::WriteToMex(int* nlhs, mxArray **plhs) {
  const char * mexfieldnames_root[] = {"formatversion", "BAopt", "planeconstraint", "VPconstraint", "cams", "cam_reproj", "points"};
  const char * mexfieldnames_BAopt[] = {"cores", "timeout_sec", "maxiter", "PTLoss", "PlaneLoss", "DerivLoss", "VPLoss"};
  const char * mexfieldnames_planeconstraints[] = {"plane", "planew", "fixPlane"};
  const char * mexfieldnames_VPconstraints[] = {"vp", "vpw", "fixVP"};  
  const char * mexfieldnames_cams[] = {"noviews", "fc", "cc", "kc", "views_orien", "views_trans", "viewids", "fixInternals", "fixOrientation", "fixTranslation", "OnPlane", "camweight", "smootherM", "vplines"};
  const char * mexfieldnames_cam_reproj[] = {"view", "pos", "weight"};
  const char * mexfieldnames_points[] = {"pt3d", "OnPlane", "fixPosition", "pointw", "reproj_view", "reproj_pos"};
  const char * mexfieldnames_vplines[] = {"l"};
  
  mxArray *tmp, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6;
  mwIndex Cellidx;
  T1 novps, nolines;

  T1* dataT1;
  T2* dataT2;
  
  plhs[0] = mxCreateStructMatrix(1, 1, 7, mexfieldnames_root);
  (*nlhs)++;
  
  tmp = mxCreateNumericMatrix(1, 1, T2mexclass, mxREAL);
  dataT2 = (T2*)mxGetData(tmp);
  dataT2[0] = this->da->Formatversion;
  mxSetField(plhs[0], 0, mexfieldnames_root[0], tmp);
       
  // Write BAopt
  tmp2 = mxCreateStructMatrix(1, 1, 7, mexfieldnames_BAopt);

      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = this->da->BAopt_nocores;
      mxSetField(tmp2, 0, mexfieldnames_BAopt[0], tmp);
      
      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = this->da->BAopt_timeout_sec;
      mxSetField(tmp2, 0, mexfieldnames_BAopt[1], tmp);
    
      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = this->da->BAopt_maxiter;
      mxSetField(tmp2, 0, mexfieldnames_BAopt[2], tmp);

      tmp = mxCreateNumericMatrix(1, 2, T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      dataT2[0] = this->da->BAopt_PTLoss[0];
      dataT2[1] = this->da->BAopt_PTLoss[1];      
      mxSetField(tmp2, 0, mexfieldnames_BAopt[3], tmp);

      tmp = mxCreateNumericMatrix(1, 2, T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      dataT2[0] = this->da->BAopt_PlaneLoss[0];
      dataT2[1] = this->da->BAopt_PlaneLoss[1];      
      mxSetField(tmp2, 0, mexfieldnames_BAopt[4], tmp);

      tmp = mxCreateNumericMatrix(1, 2, T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      dataT2[0] = this->da->BAopt_DerivLoss[0];
      dataT2[1] = this->da->BAopt_DerivLoss[1];      
      mxSetField(tmp2, 0, mexfieldnames_BAopt[5], tmp);

      tmp = mxCreateNumericMatrix(1, 2, T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      dataT2[0] = this->da->BAopt_VPLoss[0];
      dataT2[1] = this->da->BAopt_VPLoss[1];      
      mxSetField(tmp2, 0, mexfieldnames_BAopt[6], tmp);
      
  mxSetField(plhs[0], 0, mexfieldnames_root[1], tmp2);
  
  // Write plane constraints
  tmp2 = mxCreateStructMatrix(1, 1, 3, mexfieldnames_planeconstraints);
  
  uint8_t *uint8data;
    if (this->da->PlaneConstr_Planes.size()==0) // no plane entered
    {
      tmp = mxCreateNumericMatrix(0, this->da->PlaneConstr_Planes.size(), T2mexclass, mxREAL);
      mxSetField(tmp2, 0, mexfieldnames_planeconstraints[0], tmp); 
      tmp = mxCreateNumericMatrix(0, this->da->PlaneConstr_Planes.size(), T2mexclass, mxREAL);
      mxSetField(tmp2, 0, mexfieldnames_planeconstraints[1], tmp); 
      
      tmp = mxCreateNumericMatrix(0, this->da->PlaneConstr_Planes.size(), mxUINT8_CLASS, mxREAL);
      mxSetField(tmp2, 0, mexfieldnames_planeconstraints[2], tmp);
    }
    else
    {
      tmp = mxCreateNumericMatrix(3, this->da->PlaneConstr_Planes.size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetPr(tmp);
      for (T1 i = 0; i<this->da->PlaneConstr_Planes.size() ; ++i) {
        dataT2[3*i+0] = this->da->PlaneConstr_Planes[i][0];
        dataT2[3*i+1] = this->da->PlaneConstr_Planes[i][1];
        dataT2[3*i+2] = this->da->PlaneConstr_Planes[i][2];
      }
      mxSetField(tmp2, 0, mexfieldnames_planeconstraints[0], tmp); 
      
      tmp = mxCreateNumericMatrix(1, this->da->PlaneConstr_Planes.size(), T2mexclass, mxREAL);
      copy(this->da->PlaneConstr_PlanesW.begin(), this->da->PlaneConstr_PlanesW.end(), mxGetPr(tmp));   
      mxSetField(tmp2, 0, mexfieldnames_planeconstraints[1], tmp); 
      
      tmp = mxCreateNumericMatrix(1, this->da->PlaneConstr_Planes.size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(this->da->PlaneConstr_fixPlane.begin(), this->da->PlaneConstr_fixPlane.end(), uint8data); 
      mxSetField(tmp2, 0, mexfieldnames_planeconstraints[2], tmp);
    }


    
  mxSetField(plhs[0], 0, mexfieldnames_root[2], tmp2); 
  
  // Write VP lane constraints
  tmp2 = mxCreateStructMatrix(1, 1, 3, mexfieldnames_VPconstraints);

  if (this->da->VPConstr_VP.size()==0) // no plane entered
    {
      tmp = mxCreateNumericMatrix(0, this->da->VPConstr_VP.size(), T2mexclass, mxREAL);
      mxSetField(tmp2, 0, mexfieldnames_VPconstraints[0], tmp); 
      tmp = mxCreateNumericMatrix(0, this->da->VPConstr_VP.size(), T2mexclass, mxREAL);
      mxSetField(tmp2, 0, mexfieldnames_VPconstraints[1], tmp); 
      
      tmp = mxCreateNumericMatrix(0, this->da->VPConstr_VP.size(), mxUINT8_CLASS, mxREAL);
      mxSetField(tmp2, 0, mexfieldnames_VPconstraints[2], tmp);
    }
    else
    {
      tmp = mxCreateNumericMatrix(3, this->da->VPConstr_VP.size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetPr(tmp);
      for (T1 i = 0; i<this->da->VPConstr_VP.size() ; ++i) {
        dataT2[3*i+0] = this->da->VPConstr_VP[i][0];
        dataT2[3*i+1] = this->da->VPConstr_VP[i][1];
        dataT2[3*i+2] = this->da->VPConstr_VP[i][2];
      }
      mxSetField(tmp2, 0, mexfieldnames_VPconstraints[0], tmp); 
      
      tmp = mxCreateNumericMatrix(1, this->da->VPConstr_VP.size(), T2mexclass, mxREAL);
      copy(this->da->VPConstr_VPw.begin(), this->da->VPConstr_VPw.end(), mxGetPr(tmp));   
      mxSetField(tmp2, 0, mexfieldnames_VPconstraints[1], tmp); 
      
      tmp = mxCreateNumericMatrix(1, this->da->VPConstr_VP.size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);

      copy(this->da->VPConstr_fixVP.begin(), this->da->VPConstr_fixVP.end(), uint8data); 
      mxSetField(tmp2, 0, mexfieldnames_VPconstraints[2], tmp);
    }
    
  mxSetField(plhs[0], 0, mexfieldnames_root[3], tmp2); 
  

  // Write Cameras
  T1 nocams = this->da->Cams_noviews.size();
  tmp2 = mxCreateStructMatrix(1, nocams , 14, mexfieldnames_cams);
  for (T1 i = 0; i<nocams; ++i) {
      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = this->da->Cams_noviews[i];
      mxSetField(tmp2, i, mexfieldnames_cams[0], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_fc[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(this->da->Cams_fc[i].begin(), this->da->Cams_fc[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[1], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_cc[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(this->da->Cams_cc[i].begin(), this->da->Cams_cc[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[2], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_kc[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(this->da->Cams_kc[i].begin(), this->da->Cams_kc[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[3], tmp);

      tmp = mxCreateNumericMatrix(3, this->da->Cams_view_orien[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      for (T1 j=0; j< this->da->Cams_view_orien[i].size(); ++j)
      {
        dataT2[j*3+0] = this->da->Cams_view_orien[i][j][0];
        dataT2[j*3+1] = this->da->Cams_view_orien[i][j][1];
        dataT2[j*3+2] = this->da->Cams_view_orien[i][j][2];
      }
      mxSetField(tmp2, i, mexfieldnames_cams[4], tmp);

      tmp = mxCreateNumericMatrix(3, this->da->Cams_view_trans[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      for (T1 j=0; j< this->da->Cams_view_trans[i].size(); ++j)
      {
        dataT2[j*3+0] = this->da->Cams_view_trans[i][j][0];
        dataT2[j*3+1] = this->da->Cams_view_trans[i][j][1];
        dataT2[j*3+2] = this->da->Cams_view_trans[i][j][2];
      }
      mxSetField(tmp2, i, mexfieldnames_cams[5], tmp);
      
      tmp = mxCreateNumericMatrix(1, this->da->Cams_viewids[i].size(), T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      copy(this->da->Cams_viewids[i].begin(), this->da->Cams_viewids[i].end(), dataT1); 
      mxSetField(tmp2, i, mexfieldnames_cams[6], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_fixInternals[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(this->da->Cams_fixInternals[i].begin(), this->da->Cams_fixInternals[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[7], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_fixOrientation[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(this->da->Cams_fixOrientation[i].begin(), this->da->Cams_fixOrientation[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[8], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_fixTranslation[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(this->da->Cams_fixTranslation[i].begin(), this->da->Cams_fixTranslation[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[9], tmp);

      tmp = mxCreateNumericMatrix(1, this->da->Cams_OnPlane[i].size(), T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      copy(this->da->Cams_OnPlane[i].begin(), this->da->Cams_OnPlane[i].end(), dataT1); 
      mxSetField(tmp2, i, mexfieldnames_cams[10], tmp);
      
      tmp = mxCreateNumericMatrix(1, this->da->Cams_Weight[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(this->da->Cams_Weight[i].begin(), this->da->Cams_Weight[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[11], tmp);
      
      if (this->da->Cams_smootherM[i].size()==1)
      {
          tmp = mxCreateNumericMatrix(this->da->Cams_smootherM[i][0].size(), 1, T2mexclass, mxREAL);
          dataT2 = (T2*)mxGetData(tmp);
          copy(this->da->Cams_smootherM[i][0].begin(), this->da->Cams_smootherM[i][0].end(), dataT2); 
      }
      else if (this->da->Cams_smootherM[i].size()==this->da->Cams_noviews[i])
      {
        tmp = mxCreateNumericMatrix(this->da->Cams_noviews[i], this->da->Cams_noviews[i], T2mexclass, mxREAL);
        dataT2 = (T2*)mxGetData(tmp);
        for (T1 j=0; j< this->da->Cams_noviews[i]; ++j)
          for (T1 k=0; k< this->da->Cams_noviews[i]; ++k)
            dataT2[k + j*this->da->Cams_noviews[i]] = this->da->Cams_smootherM[i][j][k];
      }
      else
      { 
        tmp = mxCreateNumericMatrix(0, 0, T2mexclass, mxREAL);
      }
      mxSetField(tmp2, i, mexfieldnames_cams[12], tmp);
      
      // set VP constraints per camera
      if (this->da->Cams_VP_ID[i].size()==0)
        tmp = mxCreateNumericMatrix(0, 0, T2mexclass, mxREAL);
      else
      {
        tmp = mxCreateStructMatrix(1, this->da->Cams_noviews[i] , 1, mexfieldnames_vplines);
        
        for (T1 j = 0; j<this->da->Cams_noviews[i]; ++j)
        {
          novps = this->da->Cams_VP_ID[i][j].size();
          tmp3 = mxCreateCellMatrix(novps, 2);
          
          for (T1 k = 0; k < novps; ++k)
          {
            mwSize subs[2];
            subs[0]=k; subs[1]=0;
            Cellidx = mxCalcSingleSubscript(tmp3, 2, subs);
            
            tmp5 = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
            dataT1 = (T1*)mxGetPr(tmp5);
            dataT1[0] = this->da->Cams_VP_ID[i][j][k];

            mxSetCell(tmp3, Cellidx, tmp5);
            
            
            subs[0]=k; subs[1]=1;
            Cellidx = mxCalcSingleSubscript(tmp3, 2, subs);
            nolines = this->da->Cams_VPlines[i][j][k].size()/4;
            tmp5 = mxCreateNumericMatrix(nolines, 4, T2mexclass, mxREAL);
            dataT2 = (T2*)mxGetPr(tmp5);
            for (T1 ki = 0; ki < 4*nolines; ++ki)
              dataT2[ki] = this->da->Cams_VPlines[i][j][k][ki];
            mxSetCell(tmp3, Cellidx, tmp5);            
          }
          
          mxSetField(tmp, j, mexfieldnames_vplines[0], tmp3);
        }
      }
      mxSetField(tmp2, i, mexfieldnames_cams[13], tmp);
      
  }
  mxSetField(plhs[0], 0, mexfieldnames_root[4], tmp2);
  
  // Write mutual camera visibility
  tmp2 = mxCreateStructMatrix(1, 1, 3, mexfieldnames_cam_reproj);
  T1 nomutualvis = this->da->Cams_reproj_view.size();

  if (nomutualvis==0)
  {
    tmp = mxCreateNumericMatrix(0, 0, T1mexclass, mxREAL);
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[0], tmp); 
    tmp = mxCreateNumericMatrix(0, 0,  T2mexclass, mxREAL);
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[1], tmp); 
    tmp = mxCreateNumericMatrix(0, 0, T2mexclass, mxREAL);
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[2], tmp); 
  }
  else
  {
    tmp = mxCreateNumericMatrix(nomutualvis, 2, T1mexclass, mxREAL);
    dataT1 = (T1*)mxGetPr(tmp);
    for (T1 i = 0; i< nomutualvis  ; ++i) {
      dataT1[i            ] = this->da->Cams_reproj_view[i][0];
      dataT1[i+nomutualvis] = this->da->Cams_reproj_view[i][1];
    }
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[0], tmp); 

    tmp = mxCreateNumericMatrix(nomutualvis, 2,  T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    for (T1 i = 0; i< nomutualvis  ; ++i) {
      dataT2[i            ] = this->da->Cams_reproj_pos[i][0];
      dataT2[i+nomutualvis] = this->da->Cams_reproj_pos[i][1];
    }
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[1], tmp); 

    tmp = mxCreateNumericMatrix(1, nomutualvis, T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    copy(this->da->Cams_reproj_weight.begin(), this->da->Cams_reproj_weight.end(), dataT2); 
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[2], tmp); 
  }   
  mxSetField(plhs[0], 0, mexfieldnames_root[5], tmp2);

  // Write points
  tmp2 = mxCreateStructMatrix(1, 1, 6, mexfieldnames_points);
  T1 nopts = this->da->Pt3d.size();

    tmp = mxCreateNumericMatrix(1, nopts, mxUINT8_CLASS, mxREAL);
    uint8data= (uint8_t*)mxGetPr(tmp);
    copy(this->da->PtFixPosition.begin(), this->da->PtFixPosition.end(), uint8data); 
    mxSetField(tmp2, 0, mexfieldnames_points[2], tmp); 

    tmp = mxCreateNumericMatrix(1, nopts, T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    copy(this->da->PtWeight.begin(), this->da->PtWeight.end(), dataT2); 
    mxSetField(tmp2, 0, mexfieldnames_points[3], tmp); 

    tmp = mxCreateNumericMatrix(3, nopts, T2mexclass, mxREAL);
    tmp3 = mxCreateCellMatrix(1, nopts);
    tmp4 = mxCreateCellMatrix(1, nopts);
    tmp5 = mxCreateCellMatrix(1, nopts);

    for (T1 i = 0; i< nopts; ++i)
    {
      dataT2 = (T2*)mxGetPr(tmp);
      dataT2[i*3+0] = this->da->Pt3d[i][0];
      dataT2[i*3+1] = this->da->Pt3d[i][1];
      dataT2[i*3+2] = this->da->Pt3d[i][2];

      if (this->da->PtOnPlane[i].size() > 0)
      {
  tmp6 = mxCreateNumericMatrix(1, this->da->PtOnPlane[i].size(), T1mexclass, mxREAL);
  dataT1 = (T1*)mxGetPr(tmp6);
  copy(this->da->PtOnPlane[i].begin(), this->da->PtOnPlane[i].end(), dataT1); 
  mxSetCell(tmp3, i, tmp6);
      }
      
      tmp6 = mxCreateNumericMatrix(1, this->da->PtReprojView[i].size(), T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetPr(tmp6);
      copy(this->da->PtReprojView[i].begin(), this->da->PtReprojView[i].end(), dataT1); 
      mxSetCell(tmp4, i, tmp6);

      tmp6 = mxCreateNumericMatrix(2, this->da->PtReprojPos[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetPr(tmp6);
      for (T1 j=0; j< this->da->PtReprojPos[i].size(); ++j)
  copy(this->da->PtReprojPos[i][j].begin(), this->da->PtReprojPos[i][j].end(), dataT2+j*2);   
      mxSetCell(tmp5, i, tmp6);

    }
    mxSetField(tmp2, 0, mexfieldnames_points[0], tmp);
    mxSetField(tmp2, 0, mexfieldnames_points[1], tmp3);
    mxSetField(tmp2, 0, mexfieldnames_points[4], tmp4); 
    mxSetField(tmp2, 0, mexfieldnames_points[5], tmp5);
  
  mxSetField(plhs[0], 0, mexfieldnames_root[6], tmp2); 
  

  // return residual
  plhs[1] = mxCreateNumericMatrix(1, 1, T2mexclass, mxREAL);
  (*nlhs)++;
  
  dataT2 = (T2*)mxGetData(plhs[1]);
  dataT2[0] = this->da->Residual;
  

};

