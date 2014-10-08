  
#ifndef _Residual_H__
#define _Residual_H__


#include <valarray>
#include <iostream>
#include <numeric>

//#include "ceres/rotation.h"

// Residual function objects are declared and given all initial parameters at construction. 
// The selection of fixed or free variables is passed also at construction time.
// Correspondingly, the operator() is overloaded to accept the required free parameters in each call.
//
// If the functor is used in DynamicAutoDiffCostFunction and not AutoDiffCostFunction 
// the number of free variables is unknown at compile time. A pointer to an array of pointers 
// to parameters is passed. This is used in the functor CamPathDerError.

namespace residual 
{
 
// VP LINE IMAGE BASED ENDPOINT ERROR FUNCTOR
struct RepVPError {
  RepVPError(vector<double> endpoint_, vector<double> midpoint_,
           valarray<double> ori_, 
           valarray<double> vp_, 
           vector<double> fc_, vector<double> cc_,  vector<double> kc_, 
           double vpweight_, double camweight_, vector<bool> paramoptions_)
              
  // pass all camera + vp parameters as if they were static. 
  // Last valarray determines which arguments should be considered static. 
  // Non-static arguments are ignored and accepted as free parameters to overloaded operator().
  {   
    endpoint_s_2D = endpoint_;
    midpoint_s_2D = midpoint_;
    endpoint_s_3D.resize(3);
    midpoint_s_3D.resize(3);
    linedir_org.resize(2);
    
    linedir_org[0] = midpoint_s_2D[0] - endpoint_s_2D[0];
    linedir_org[1] = midpoint_s_2D[1] - endpoint_s_2D[1];
    double norm = sqrt(linedir_org[0]*linedir_org[0] + linedir_org[1]*linedir_org[1]);
    linedir_org[0] /= norm;
    linedir_org[1] /= norm;
    
    vp_s = vp_;
    fc_s = fc_;
    cc_s = cc_;
    kc_s = kc_;
    ori_s = ori_;
    paramoptions = paramoptions_;
    vpweight = vpweight_;
    camweight = camweight_;
    
    if (cc_s.size()==2)
      modelpp = 1;
    else
      modelpp = 0;
      
    radmodel = kc_s.size();
    fcnoelem = fc_s.size();
    
    callopt=0;
    callopt += (paramoptions[4] << 5);
    callopt += (paramoptions[5] << 4);
    callopt += (paramoptions[0] << 3);
    callopt += (paramoptions[1] << 2);
    callopt += (paramoptions[2] << 1);
    callopt += (paramoptions[3] << 0);
    

    //note: paramoptions[5] (for translation) is also passed, but is always 1 for this VP functor since translation is irrelevant for VPs
    
    if ( paramoptions[1] & paramoptions[2] & paramoptions[3] ) // if fc,cc,kc are fixed, precompute radial distortion for line mid/endpoint
    {
      // 1/2. bring 2D line endpoints to normalized image coordinates, and undistort points
      // 3. compute original 2D line midpoint, project to normalized image coordinates, undistort
      double res[3];
      resutil::ImgPt_to_NormCoordSphere(&(fc_s[0]), &(cc_s[0]), &(kc_s[0]), fcnoelem, modelpp, radmodel, &(midpoint_s_2D[0]), res);
      midpoint_s_3D[0] = res[0]; 
      midpoint_s_3D[1] = res[1];
      midpoint_s_3D[2] = res[2];
      
      resutil::ImgPt_to_NormCoordSphere(&(fc_s[0]), &(cc_s[0]), &(kc_s[0]), fcnoelem, modelpp, radmodel, &(endpoint_s_2D[0]), res);
      endpoint_s_3D[0] = res[0]; 
      endpoint_s_3D[1] = res[1];
      endpoint_s_3D[2] = res[2];
    }

  }
   template <typename T>
   bool evalresidual(const T* const ori,
                     const void *,
                   const T* const vp,
                   const T* const fc,
                   const T* const cc,
                   const T* const kc, 
                   T* residual) const {
                     
                    T midpoint[3], endpoint[3], planevec[3], vp_rot[3], p[2], linedir_proj[2], norm, ang;     
                    T ori_tmp[3]; ori_tmp[0] = T(ori_s[0]); ori_tmp[1] = T(ori_s[1]); ori_tmp[2] = T(ori_s[2]);
                    T kc_tmp[radmodel]; for (int i = 0; i < radmodel; ++i) kc_tmp[i] = T(kc_s[i]);                    
                    T fc_tmp[fcnoelem]; for (int i = 0; i < fcnoelem; ++i) fc_tmp[i] = T(fc_s[i]);                    
                    T cc_tmp[2]; if (modelpp) { cc_tmp[0] = T(cc_s[0]); cc_tmp[1] = T(cc_s[1]); }
                    T vp_tmp[3]; vp_tmp[0] = T(vp_s[0]); vp_tmp[1] = T(vp_s[1]); vp_tmp[2] = T(vp_s[2]);
                    T midpoint_tmp[2]; midpoint_tmp[0] = T(midpoint_s_2D[0]); midpoint_tmp[1] = T(midpoint_s_2D[1]);
                    T endpoint_tmp[2]; endpoint_tmp[0] = T(endpoint_s_2D[0]); endpoint_tmp[1] = T(endpoint_s_2D[1]);
                    
                    // 0. Apply rotation to 3D vanishing direction  
                    if (paramoptions[0]) { // vp is fixed
                      if (paramoptions[4]) { // orientation is fixed
                        
                        resutil::AngleAxisRotatePoint(ori_tmp, vp_tmp, vp_rot); }
                      else {
                        resutil::AngleAxisRotatePoint(ori,     vp_tmp, vp_rot); }  //
                    } else {        // vp is free
                      if (paramoptions[4]) { // orientation is fixed
                        resutil::AngleAxisRotatePoint(ori_tmp, vp, vp_rot); }
                      else {
                        resutil::AngleAxisRotatePoint(ori,     vp, vp_rot); }  
                    }
                    
                    
                    // 1/2. bring 2D line endpoints to normalized image coordinates, and undistort points
                    // 3. compute original 2D line midpoint, project to normalized image coordinates, undistort

                    if (paramoptions[1]) { // fc is fixed
                      if (paramoptions[2]) { // cc is fixed
                        if (paramoptions[3]) {  // kc is fixed
                          midpoint[0] = T(midpoint_s_3D[0]); midpoint[1] = T(midpoint_s_3D[1]); midpoint[2] = T(midpoint_s_3D[2]);
                          endpoint[0] = T(endpoint_s_3D[0]); endpoint[1] = T(endpoint_s_3D[1]); endpoint[2] = T(endpoint_s_3D[2]); }
                        else {                  // kc is free
                          resutil::ImgPt_to_NormCoordSphere(fc_tmp, cc_tmp, kc, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc_tmp, cc_tmp, kc, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                        
                      } else {            // cc is free
                        if (paramoptions[3]) {  // kc is fixed
                          resutil::ImgPt_to_NormCoordSphere(fc_tmp, cc, kc_tmp, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc_tmp, cc, kc_tmp, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                        else {                  // kc is free
                          resutil::ImgPt_to_NormCoordSphere(fc_tmp, cc, kc, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc_tmp, cc, kc, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                    } } else {               // fc is free
                      if (paramoptions[2]) { // cc is fixed
                        if (paramoptions[3]) {  // kc is fixed
                          resutil::ImgPt_to_NormCoordSphere(fc, cc_tmp, kc_tmp, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc, cc_tmp, kc_tmp, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                        else {                  // kc is free
                          resutil::ImgPt_to_NormCoordSphere(fc, cc_tmp, kc, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc, cc_tmp, kc, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                        
                      } else {                 // cc is free  
                        if (paramoptions[3]) {  // kc is fixed
                          resutil::ImgPt_to_NormCoordSphere(fc, cc, kc_tmp, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc, cc, kc_tmp, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                        else {                  // kc is free
//cout << "Recompute Line mid/endpoints in image normalized coordinates"<< endl;      
                          resutil::ImgPt_to_NormCoordSphere(fc, cc, kc, fcnoelem, modelpp, radmodel, midpoint_tmp, midpoint);
                          resutil::ImgPt_to_NormCoordSphere(fc, cc, kc, fcnoelem, modelpp, radmodel, endpoint_tmp, endpoint); } 
                    } }

                     // 4. use line midpoint in normalized image coordinates to compute greater circle (VP plane through origin)
                     resutil::CrossProduct(vp_rot, midpoint, planevec);

                     resutil::VecNormalizeL2_D3(planevec);
                     
                      // 5). project line endpoint (in n. i. c.) to greater circle
                     T proj_error = resutil::DotProduct(planevec, endpoint); // 3D distance  of endpoint to interpretation plane (because both vectors are unit vectors)
                     
                     endpoint[0] = endpoint[0] - proj_error*planevec[0];
                     endpoint[1] = endpoint[1] - proj_error*planevec[1];
                     endpoint[2] = endpoint[2] - proj_error*planevec[2];
                     
                     resutil::VecNormalizeL2_D3(endpoint);

                          
                    // 6) transform projected point on greater circle back to 2D image coordinates, distort!  
                    endpoint[0] = endpoint[0] / endpoint[2];
                    endpoint[1] = endpoint[1] / endpoint[2];
                    if (paramoptions[3]==0) // kc is free parameter        
                      resutil::RadDist(kc, radmodel, endpoint, p);
                    else    // ori is fixed
                      resutil::RadDist(kc_tmp, radmodel, endpoint, p);
    
                    // Project into camera
                    if (paramoptions[1]==0)  // fc is free parameter      
                    {
                      if (fcnoelem==2)
                      {
                        p[0] = p[0]*fc[0];
                        p[1] = p[1]*fc[1];
                      }
                      else
                      {
                        p[0] = p[0]*fc[0];
                        p[1] = p[1]*fc[0];
                      }
                    }
                    else      // fc is fixed
                    { 
                      if (fcnoelem==2)
                      {
                        p[0] = p[0]*T(fc_s[0]);
                        p[1] = p[1]*T(fc_s[1]);
                      }
                      else
                      {
                        p[0] = p[0]*T(fc_s[0]);
                        p[1] = p[1]*T(fc_s[0]);
                      }      
                    }
                    
                    // Principal point
                    if (modelpp)
                    {
                      //if (cc!=NULL)  // cc is free parameter
                      if (paramoptions[2]==0)  // cc is free parameter
                      {
                          p[0] += cc[0];
                          p[1] += cc[1];
                      }
                      else      // cc is fixed
                      {
                          p[0] += T(cc_s[0]);
                          p[1] += T(cc_s[1]);
                      }
                    }
                    
                    // 7) Transformed endpoints and original 2D endpoints should overlap. Project transformed endpoints on original 2D image line.
                    linedir_proj[0] = midpoint_s_2D[0] - p[0];
                    linedir_proj[1] = midpoint_s_2D[1] - p[1];
                    norm = sqrt(linedir_proj[0]*linedir_proj[0] + linedir_proj[1]*linedir_proj[1]);
                    
                    
                    ang = (T(linedir_org[0])*linedir_proj[0] +  T(linedir_org[1])*linedir_proj[1]) / norm; // cosine between vectors in image space
                    
                    endpoint[0] =  (norm * ang) * T(linedir_org[0]); // vector along observed 2D line to right-angle closest point for projected line endpoint
                    endpoint[1] =  (norm * ang) * T(linedir_org[1]);

                    // 8) Projection error vector of both endpoints is residual 
                    
                    residual[0] = linedir_proj[0] - endpoint[0];
                    residual[1] = linedir_proj[1] - endpoint[1];
                    
                    //cout <<  residual[0] << " " <<  residual[1] << endl;
                    
                    residual[0] *= T(camweight*vpweight);
                    residual[1] *= T(camweight*vpweight);                    
                     
                    return true;
                   }
                   
    // 6 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3,const T* const p4,const T* const p5,const T* const p6,T* residual) const {  bool ret;
    ret=this->evalresidual<T>(p1            ,p2            ,p3            ,p4            ,p5            ,p6            , residual);
    return ret; };

                   
    // 5 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3,const T* const p4,const T* const p5, T* residual) const {  bool ret;
    switch( callopt ) {
    case  1: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,p4            ,p5            ,     NULL     , residual); break;
    case  2: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,p4            ,     NULL     ,p5            , residual); break; 
    case  4: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,p4            ,p5            , residual); break; 
    case  8: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,p4            ,p5            , residual); break; 
    case 16: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,p4            ,p5            , residual); break; 
    case 32: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,p4            ,p5            , residual); break; 
    }
    return ret; };

    // 4 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3,const T* const p4, T* residual) const {  bool ret;
    switch( callopt ) {
    case  3: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,p4            ,     NULL     ,     NULL     , residual); break;
    case  5: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,p4            ,     NULL     , residual); break;
    case  6: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,     NULL     ,p4            , residual); break;    
    case  9: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,p4            ,     NULL     , residual); break; 
    case 10: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,     NULL     ,p4            , residual); break; 
    case 12: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,p3            ,p4            , residual); break; 
    case 17: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,p4            ,     NULL     , residual); break; 
    case 18: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,     NULL     ,p4            , residual); break; 
    case 20: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,p3            ,p4            , residual); break; 
    case 24: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,p3            ,p4            , residual); break; 
    case 33: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,p4            ,     NULL     , residual); break;  
    case 34: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,     NULL     ,p4            , residual); break; 
    case 36: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,p3            ,p4            , residual); break;
    case 40: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,p3            ,p4            , residual); break; 
    case 48: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,p3            ,p4            , residual); break; 
    }
    return ret; };

    // 3 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3, T* residual) const {  bool ret;
    //cout << "HERE3"  << callopt << endl;
    switch( callopt ) {
    case  7: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,     NULL     ,     NULL     , residual); break;
    case 11: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,     NULL     ,     NULL     , residual); break; 
    case 13: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,p3            ,     NULL     , residual); break; 
    case 14: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,     NULL     ,p3            , residual); break; 
    case 19: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,     NULL     ,     NULL     , residual); break; 
    case 21: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,p3            ,     NULL     , residual); break; 
    case 22: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,     NULL     ,p3            , residual); break; 
    case 25: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,p3            ,     NULL     , residual); break; 
    case 26: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,     NULL     ,p3            , residual); break; 
    case 28: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,p2            ,p3            , residual); break; 
    case 35: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,     NULL     ,     NULL     , residual); break; 
    case 37: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,p3            ,     NULL     , residual); break; 
    case 38: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,     NULL     ,p3            , residual); break; 
    case 41: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,p3            ,     NULL     , residual); break; 
    case 42: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,     NULL     ,p3            , residual); break; 
    case 44: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,p2            ,p3            , residual); break; 
    case 49: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,p3            ,     NULL     , residual); break; 
    case 50: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,     NULL     ,p3            , residual); break; 
    case 52: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,p2            ,p3            , residual); break; 
    case 56: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,p2            ,p3            , residual); break; 
    }
    return ret; };

    // 2 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2, T* residual) const {  bool ret;
    switch( callopt ) {
    case 15: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 23: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 27: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,     NULL     ,     NULL     , residual); break; 
    case 29: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,p2            ,     NULL     , residual); break; 
    case 30: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,     NULL     ,p2            , residual); break; 
    case 39: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 43: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,     NULL     ,     NULL     , residual); break; 
    case 45: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,p2            ,     NULL     , residual); break; 
    case 46: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,     NULL     ,p2            , residual); break; 
    case 51: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,     NULL     ,     NULL     , residual); break; 
    case 53: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,p2            ,     NULL     , residual); break; 
    case 54: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,     NULL     ,p2            , residual); break; 
    case 57: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,p2            ,     NULL     , residual); break; 
    case 58: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,     NULL     ,p2            , residual); break; 
    case 60: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,     NULL     ,p1            ,p2            , residual); break; 
    }
    return ret; };

    // 1 Parameter + residual
    template <typename T>  bool operator()(const T* const p1, T* residual) const {  bool ret;
    switch( callopt ) {
    case 31: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 47: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 55: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 59: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,     NULL     ,     NULL     , residual); break;  
    case 61: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,     NULL     ,p1            ,     NULL     , residual); break; 
    case 62: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,     NULL     ,     NULL     ,p1            , residual); break; 
    }
    return ret; };
    
    
    
    
    
    
                   
  private: 
   bool modelpp; // model principal point
   int radmodel; // which rad. dist model to use
   int fcnoelem; // number of elements for focal length, |{fx,fy}|=2 or |{fxy}|=1

   double camweight;
   double vpweight;
   
   vector<double> endpoint_s_2D; // line observations midpoint in image, static, used for computing line reprojection error
   vector<double> midpoint_s_2D; // line observations endpoint in image, static, used for computing VP interpretation plane through 2D line midpoint
   vector<double> linedir_org; // 2D direction from midpoint to endpoint
   vector<double> endpoint_s_3D; // same as above, but in normalized image coordinates, mapped onto unit sphere
   vector<double> midpoint_s_3D; // same as above, but in normalized image coordinates, mapped onto unit sphere
   
   valarray<double> vp_s; // 3d vp direction, static
   vector<double> fc_s; // focal length, 1 or 2-dim, static
   vector<double> cc_s; // principal point, static
   vector<double> kc_s; // 1,2,4,5 radial distortion, Model after D.C. Brown, et al.'Decentering Distortion of lenses', 1966, static
   valarray<double> ori_s; // orientation in rodrigues format, static
   vector<bool> paramoptions; // fix parameters: elem 0: vp, 1: fc, 2: cc, 3:kc, 4:ori
   
   uint32_t callopt; // [0:63] switch for choosing which reproj. function to call
    // 5-bit = 1/0:  ori is fixed/ free
    // 4-bit = 1/0:  tr is fixed/ free
    // 3-bit = 1/0:  pt is fixed/ free
    // 2-bit = 1/0:  fc is fixed/ free
    // 1-bit = 1/0:  cc is fixed/ free
    // 0-bit = 1/0:  kk is fixed/ free

};

  
// NON-SMOOTHING POINT REPROJECTION FUNCTOR
struct RepError {
  RepError(vector<double> obs_, 
           valarray<double> ori_, valarray<double> tr_, 
           valarray<double> pt_, 
           vector<double> fc_, vector<double> cc_,  vector<double> kc_, 
           double ptweight_, double camweight_, vector<bool> paramoptions_)
              
  // pass all camera + point parameters as if they were static. 
  // Last valarray determines which arguments should be considered static. 
  // Non-static arguments are ignored and accepted as free parameters to overloaded operator().
  {   
    obs_s = obs_;
    pt_s = pt_;
    fc_s = fc_;
    cc_s = cc_;
    kc_s = kc_;
    ori_s = ori_;
    tr_s = tr_;
    paramoptions = paramoptions_;
    ptweight = ptweight_;
    camweight = camweight_;
    
    if (cc_s.size()==2)
      modelpp = 1;
    else
      modelpp = 0;
      
    radmodel = kc_s.size();
    fcnoelem = fc_s.size();
    
    callopt=0;
    callopt += (paramoptions[4] << 5);
    callopt += (paramoptions[5] << 4);
    callopt += (paramoptions[0] << 3);
    callopt += (paramoptions[1] << 2);
    callopt += (paramoptions[2] << 1);
    callopt += (paramoptions[3] << 0);

  }
  
   template <typename T>
   bool evalresidual(const T* const ori, 
                   const T* const tr,
                   const T* const pt,
                   const T* const fc,
                   const T* const cc,
                   const T* const kc, 
                   T* residual) const {
     
     // Rotate and shift point according to R,t of camera
     T p[3], tmp[3];
     T kc_tmp[radmodel]; for (int i = 0; i < radmodel; ++i) kc_tmp[i] = T(kc_s[i]);
     
     //p[0] = tr[0];
     //p[1] = tr[1];
     //p[2] = tr[2];
     
    // tr is given as in world ref. frame
    //if (tr!=NULL)  // tr is free parameter
    if (paramoptions[0]==0)  // pt is free parameter            
    {
        if (paramoptions[5]==0)  // tr is free parameter          
        {
          tmp[0] = pt[0] - tr[0]; 
          tmp[1] = pt[1] - tr[1]; 
          tmp[2] = pt[2] - tr[2];
        }
        else        // tr is fixed
        {
          tmp[0] = pt[0] - T(tr_s[0]); 
          tmp[1] = pt[1] - T(tr_s[1]); 
          tmp[2] = pt[2] - T(tr_s[2]);
        }
    }
    else                   // pt is fixed parameter
    {
        if (paramoptions[5]==0)  // tr is free parameter          
        {
          tmp[0] = T(pt_s[0]) - tr[0]; 
          tmp[1] = T(pt_s[1]) - tr[1]; 
          tmp[2] = T(pt_s[2]) - tr[2];
        }
        else        // tr is fixed
        {
          tmp[0] = T(pt_s[0]) - T(tr_s[0]); 
          tmp[1] = T(pt_s[1]) - T(tr_s[1]); 
          tmp[2] = T(pt_s[2]) - T(tr_s[2]);
        }
    }
            
    //if (ori!=NULL)  // ori is free parameter
    if (paramoptions[4]==0) // ori is free parameter        
      resutil::AngleAxisRotatePoint(ori, tmp, p);
    else    // ori is fixed
    {
      T ori_tmp[3]; ori_tmp[0] = T(ori_s[0]); ori_tmp[1] = T(ori_s[1]); ori_tmp[2] = T(ori_s[2]);
      resutil::AngleAxisRotatePoint(ori_tmp, tmp, p);
      //resutil::AngleAxisRotatePoint((T*)&(ori_s[0]), tmp, p);
    }

    // Radial Distortion
    tmp[0] = p[0] / p[2];
    tmp[1] = p[1] / p[2]; 
    if (paramoptions[3]==0) // kc is free parameter        
      resutil::RadDist(kc, radmodel, tmp, p);
    else    // kc is fixed
      resutil::RadDist(kc_tmp, radmodel, tmp, p);
      //resutil::RadDist((T*)&(kc_s[0]), radmodel, tmp, p);
     
    // Project into camera
    //if (fc!=NULL)  // fc is free parameter
    if (paramoptions[1]==0)  // fc is free parameter      
    {
      //cout << "FC is free" << endl;
      if (fcnoelem==2)
      {
        p[0] = p[0]*fc[0];
        p[1] = p[1]*fc[1];
      }
      else
      {
        p[0] = p[0]*fc[0];
        p[1] = p[1]*fc[0];
      }
    }
    else      // fc is fixed
    { 
      if (fcnoelem==2)
      {
        p[0] = p[0]*T(fc_s[0]);
        p[1] = p[1]*T(fc_s[1]);
      }
      else
      {
        p[0] = p[0]*T(fc_s[0]);
        p[1] = p[1]*T(fc_s[0]);
      }      
    }
     
    // Principal point
    if (modelpp)
    {
      //if (cc!=NULL)  // cc is free parameter
      if (paramoptions[2]==0)  // cc is free parameter
      {
          p[0] += cc[0];
          p[1] += cc[1];
      }
      else      // cc is fixed
      {
          p[0] += T(cc_s[0]);
          p[1] += T(cc_s[1]);
      }
    }
     
     residual[0] = p[0] - T(obs_s[0]);
     residual[1] = p[1] - T(obs_s[1]);
     
     residual[0] *= T(camweight*ptweight);
     residual[1] *= T(camweight*ptweight);
     
     
     return true;
   };
  
   
    // 6 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3,const T* const p4,const T* const p5,const T* const p6,T* residual) const {  bool ret;
    ret=this->evalresidual<T>(p1            ,p2            ,p3            ,p4            ,p5            ,p6            , residual);
    return ret; };

    // 5 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3,const T* const p4,const T* const p5, T* residual) const {  bool ret;
    switch( callopt ) {
    case  1: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,p4            ,p5            ,     NULL     , residual); break;
    case  2: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,p4            ,     NULL     ,p5            , residual); break; 
    case  4: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,p4            ,p5            , residual); break; 
    case  8: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,p4            ,p5            , residual); break; 
    case 16: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,p4            ,p5            , residual); break; 
    case 32: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,p4            ,p5            , residual); break; 
    }
    return ret; };

    // 4 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3,const T* const p4, T* residual) const {  bool ret;
    switch( callopt ) {
    case  3: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,p4            ,     NULL     ,     NULL     , residual); break;
    case  5: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,p4            ,     NULL     , residual); break;
    case  6: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,     NULL     ,p4            , residual); break;    
    case  9: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,p4            ,     NULL     , residual); break; 
    case 10: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,     NULL     ,p4            , residual); break; 
    case 12: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,p3            ,p4            , residual); break; 
    case 17: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,p4            ,     NULL     , residual); break; 
    case 18: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,     NULL     ,p4            , residual); break; 
    case 20: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,p3            ,p4            , residual); break; 
    case 24: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,p3            ,p4            , residual); break; 
    case 33: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,p4            ,     NULL     , residual); break;  
    case 34: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,     NULL     ,p4            , residual); break; 
    case 36: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,p3            ,p4            , residual); break;
    case 40: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,p3            ,p4            , residual); break; 
    case 48: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,p3            ,p4            , residual); break; 
    }
    return ret; };

    // 3 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2,const T* const p3, T* residual) const {  bool ret;
    //cout << "HERE3"  << callopt << endl;
    switch( callopt ) {
    case  7: ret=this->evalresidual<T>(p1             ,p2            ,p3            ,     NULL     ,     NULL     ,     NULL     , residual); break;
    case 11: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,p3            ,     NULL     ,     NULL     , residual); break; 
    case 13: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,p3            ,     NULL     , residual); break; 
    case 14: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,     NULL     ,p3            , residual); break; 
    case 19: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,p3            ,     NULL     ,     NULL     , residual); break; 
    case 21: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,p3            ,     NULL     , residual); break; 
    case 22: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,     NULL     ,p3            , residual); break; 
    case 25: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,p3            ,     NULL     , residual); break; 
    case 26: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,     NULL     ,p3            , residual); break; 
    case 28: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,p2            ,p3            , residual); break; 
    case 35: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,p3            ,     NULL     ,     NULL     , residual); break; 
    case 37: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,p3            ,     NULL     , residual); break; 
    case 38: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,     NULL     ,p3            , residual); break; 
    case 41: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,p3            ,     NULL     , residual); break; 
    case 42: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,     NULL     ,p3            , residual); break; 
    case 44: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,p2            ,p3            , residual); break; 
    case 49: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,p3            ,     NULL     , residual); break; 
    case 50: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,     NULL     ,p3            , residual); break; 
    case 52: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,p2            ,p3            , residual); break; 
    case 56: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,p2            ,p3            , residual); break; 
    }
    return ret; };

    // 2 Parameters + residual
    template <typename T>  bool operator()(const T* const p1,const T* const p2, T* residual) const {  bool ret;
    switch( callopt ) {
    case 15: ret=this->evalresidual<T>(p1             ,p2            ,     NULL     ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 23: ret=this->evalresidual<T>(p1             ,     NULL     ,p2            ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 27: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,p2            ,     NULL     ,     NULL     , residual); break; 
    case 29: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,p2            ,     NULL     , residual); break; 
    case 30: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,     NULL     ,p2            , residual); break; 
    case 39: ret=this->evalresidual<T>(     NULL      ,p1            ,p2            ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 43: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,p2            ,     NULL     ,     NULL     , residual); break; 
    case 45: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,p2            ,     NULL     , residual); break; 
    case 46: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,     NULL     ,p2            , residual); break; 
    case 51: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,p2            ,     NULL     ,     NULL     , residual); break; 
    case 53: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,p2            ,     NULL     , residual); break; 
    case 54: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,     NULL     ,p2            , residual); break; 
    case 57: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,p2            ,     NULL     , residual); break; 
    case 58: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,     NULL     ,p2            , residual); break; 
    case 60: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,     NULL     ,p1            ,p2            , residual); break; 
    }
    return ret; };

    // 1 Parameter + residual
    template <typename T>  bool operator()(const T* const p1, T* residual) const {  bool ret;
    switch( callopt ) {
    case 31: ret=this->evalresidual<T>(p1             ,     NULL     ,     NULL     ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 47: ret=this->evalresidual<T>(     NULL      ,p1            ,     NULL     ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 55: ret=this->evalresidual<T>(     NULL      ,     NULL     ,p1            ,     NULL     ,     NULL     ,     NULL     , residual); break; 
    case 59: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,p1            ,     NULL     ,     NULL     , residual); break;  
    case 61: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,     NULL     ,p1            ,     NULL     , residual); break; 
    case 62: ret=this->evalresidual<T>(     NULL      ,     NULL     ,     NULL     ,     NULL     ,     NULL     ,p1            , residual); break; 
    }
    return ret; };

    // 0 Parameters + residual, is not allowed by ceres

  private: 
   bool modelpp; // model principal point
   int radmodel; // which rad. dist model to use
   int fcnoelem; // number of elements for focal length, |{fx,fy}|=2 or |{fxy}|=1

   double camweight;
   double ptweight;
   
   vector<double> obs_s; // observation in image, static
   valarray<double> pt_s; // 3d point, static
   vector<double> fc_s; // focal length, 1 or 2-dim, static
   vector<double> cc_s; // principal point, static
   vector<double> kc_s; // 1,2,4,5 radial distortion, Model after D.C. Brown, et al.'Decentering Distortion of lenses', 1966, static
   valarray<double> ori_s; // orientation in rodrigues format, static
   valarray<double> tr_s; // camera center in world, static
   vector<bool> paramoptions; // fix parameters: elem 0: pt3d, 1: fc, 2: cc, 3:kc, 4:ori, 5:tr
   
   uint32_t callopt; // [0:63] switch for choosing which reproj. function to call
    // 5-bit = 1/0:  ori is fixed/ free
    // 4-bit = 1/0:  tr is fixed/ free
    // 3-bit = 1/0:  pt is fixed/ free
    // 2-bit = 1/0:  fc is fixed/ free
    // 1-bit = 1/0:  cc is fixed/ free
    // 0-bit = 1/0:  kk is fixed/ free

};


// FRAME-TO-FRAME POSITINAL DERIVATIVE FUNCTOR
struct CamPathDerError {
  CamPathDerError(vector<valarray<double>> camOr_, vector<valarray<double>> camTr_, 
                    double camweight_, double stdev3d_,  vector<double> sm_val_, 
                    vector<bool> paramoptions_)
              
  // pass all camera parameters as if they were static. 
  // Last valarray determines which arguments should be considered static. 
  // Non-static arguments are ignored and accepted as free parameters to overloaded operator().
  {   
    ori_s = camOr_;
    tr_s = camTr_;
    
    paramoptions = paramoptions_;
    camweight = camweight_;
    stdev3d = stdev3d_;
    sm_val = sm_val_;
  }
   
   
   
  template<typename T>
  bool operator()(T const* const* vars, T* residuals) const 
  { 
    int32_t vars_idx = 0;    
    
    T cp[3] = {}; // camera position
    T ct[3] = {}; // camera target
    T cu[3] = {}; // camera upvector
    
    T RotM[9] = {};
    
    for (uint32_t i=0; i < sm_val.size(); ++i) // over all cameras that should be used in smoothing
    {
      // if camera has fixed orientation, use static value, other use passed variable
      if (paramoptions[i]) // or is fixed
      {
        T ori_tmp[3]; ori_tmp[0] = T(ori_s[i][0]); ori_tmp[1] = T(ori_s[i][1]); ori_tmp[2] = T(ori_s[i][2]);
        resutil::AngleAxisToRMatrix(ori_tmp, RotM);  
      
      }
      else   // or is free
      {
        resutil::AngleAxisToRMatrix(&(vars[vars_idx][0]), RotM);
        vars_idx++;
      }
    
      ct[0] += T(sm_val[i])*RotM[6];
      ct[1] += T(sm_val[i])*RotM[7];
      ct[2] += T(sm_val[i])*RotM[8];

      cu[0] += T(sm_val[i])*RotM[3];
      cu[1] += T(sm_val[i])*RotM[4];
      cu[2] += T(sm_val[i])*RotM[5];
    }
         
    T pos[3] = {};
    //T last[3] = {};
    for (uint32_t i=0; i < sm_val.size(); ++i) // over all cameras that should be used in smoothing
    {
      
      // if camera has fixed orientation, use static value, other use passed variable
      if (paramoptions[sm_val.size() + i]) // tr fixed
      {
        pos[0] = T(tr_s[i][0]);
        pos[1] = T(tr_s[i][1]);
        pos[2] = T(tr_s[i][2]);
      }
      else   // or is free
      {
        pos[0] = vars[vars_idx][0];
        pos[1] = vars[vars_idx][1];
        pos[2] = vars[vars_idx][2];
        vars_idx++;
      }

      cp[0] += T(sm_val[i])*pos[0];
      cp[1] += T(sm_val[i])*pos[1];
      cp[2] += T(sm_val[i])*pos[2];
    }
    
    residuals[0] = T(camweight) * cp[0];
    residuals[1] = T(camweight) * cp[1];
    residuals[2] = T(camweight) * cp[2];
    
    return true;    
  };
  
  
  private: 
   double camweight;
   double stdev3d; // downscaling of whole dataset, use to scale unit vectors for directional smoothing
   
   vector<valarray<double>> ori_s; // orientation in rodrigues format, static, for all cameras considered in this smoothing
   vector<valarray<double>> tr_s; // camera center in world, static, for all cameras considered in this smoothing
   vector<double> sm_val;  // Smoothing values from smoothing spline matrix, 1xN
   vector<bool> paramoptions; // fix parameters: 1:(N-1) ori, (N):(2*N-1) tr
};


struct PlaneError {
  PlaneError(valarray<double> pl_, valarray<double> pt_, 
             double planeweight_, double pointweight_, vector<bool> paramoptions_)
  {
      pl_s.resize(3); pl_s[0] = pl_[0]; pl_s[1] = pl_[1]; pl_s[2] = pl_[2];
      pt_s.resize(3); pt_s[0] = pt_[0]; pt_s[1] = pt_[1]; pt_s[2] = pt_[2];
      
      planeweight = planeweight_;
      pointweight = pointweight_;
      
      // 0: fix nothing, 1: fix point, 2: fix plane, 3: fix point and fix plane
      paramoptions = 0;
      paramoptions += paramoptions_[0] << 0; 
      paramoptions += paramoptions_[1] << 1; 
      
}

  // 2 parameters + residual
  template <typename T>  bool operator()(const T* const p1,const T* const p2, T* residual) const {  bool ret;
  ret=this->evalresidual<T>(p1            ,p2            ,residual);
  return ret; };

  // 1 Parameter + residual
  template <typename T>  bool operator()(const T* const p1, T* residual) const {  bool ret;
  switch( paramoptions ) {
  case 1: ret=this->evalresidual<T>(p1             ,     NULL     , residual); break; 
  case 2: ret=this->evalresidual<T>(     NULL      ,p1            , residual); break; 
  }
  return ret; };  
  
  // 0 Parameters + residual
  template <typename T>  bool operator()(T* residual) const {  bool ret;
  ret=this->evalresidual<T>(     NULL      ,     NULL     , residual); 
  return ret; };  
    
  template <typename T>
  bool evalresidual(const T* const pl,
                    const T* const pt,
                    T* residuals) const {
                    
    
    T p[3];
    
    if (paramoptions==1 || paramoptions==3) // point is fixed
    {
      p[0] = T(pt_s[0]);
      p[1] = T(pt_s[1]);
      p[2] = T(pt_s[2]);
    }
    else     // point is free parameter
    {
      p[0] = pt[0];
      p[1] = pt[1];
      p[2] = pt[2];
    }
    
    T plane[3];
    T norm;

    if (paramoptions==2 || paramoptions==3)   // plane is fixed
    {
      norm = sqrt(T(pl_s[0])*T(pl_s[0]) + T(pl_s[1])*T(pl_s[1]) + T(pl_s[2])*T(pl_s[2]));
      plane[0] = T(pl_s[0])/norm;
      plane[1] = T(pl_s[1])/norm;
      plane[2] = T(pl_s[2])/norm;
    }
    else   // plane is free parameter
    {
      norm = sqrt(pl[0]*pl[0] + pl[1]*pl[1] + pl[2]*pl[2]);
      plane[0] = pl[0]/norm;
      plane[1] = pl[1]/norm;
      plane[2] = pl[2]/norm;
    }


    residuals[0] = abs((plane[0]*p[0] + plane[1]*p[1] + plane[2]*p[2] ) - norm);

    residuals[0] *= T(planeweight)*T(pointweight);
    
    return true;
  }

  valarray<double> pl_s; // plane vector , static
  valarray<double> pt_s; // point ,  static
  double planeweight;
  double pointweight;
  int paramoptions; // 0: fix nothing, 1: fix point, 2: fix plane, 3: fix point and fix plane
}; 

}

#endif // def Residual