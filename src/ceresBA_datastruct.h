#ifndef _CeresBA_DataStruct__
#define _CeresBA_DataStruct__

#include <vector>
#include <valarray>

using namespace std;

template<class T1, class T2> 
struct BAdata // set T1 as uint32/uint64 and T2 as double/float
{
	T2 Formatversion;  // Version of DataStruct format
	T1 BAopt_nocores;  //  Number of cores for Ceres Solver
	T1 BAopt_timeout_sec;  // Timeout in seconds for Ceres Solver  
	T1 BAopt_maxiter;  // Max. iterations for Ceres Solver
	vector<T2> BAopt_PTLoss;  // 2-elem vector: [0]: Loss function, [1] parameter for loss function, 
                                                // [0]==0 Trivial loss, 
                                                // [0]==1 HuberLoss,  [0]==2 SoftLOneLoss, [0]==3 CauchyLoss, [0]==4 ArctanLoss,
  vector<T2> BAopt_PlaneLoss;  // 2-elem vector, same options as for PTLoss
  vector<T2> BAopt_DerivLoss;  // 2-elem vector, same options as for PTLoss
  vector<T2> BAopt_VPLoss;  // 2-elem vector, same options as for PTLoss

  valarray<T2> median3d; // for normalization of data, is computed automatically by Parser.NormalizeData()
  T2 stdev3d;            // for normalization of data, is computed automatically by Parser.NormalizeData()
	  
	vector<valarray<T2>> PlaneConstr_Planes; // 3xN Planes, for point and camera constraints, plane is 3-dim parametrization for plane n: [n(1)/norm(n) n(2)/norm(n) n(3)/norm(n) -norm(n)]
	vector<T2> PlaneConstr_PlanesW; // 1xN Planes, plane weights
	vector<bool> PlaneConstr_fixPlane;  // Fix groundplane or keep it as free variable

  vector<valarray<T2>> VPConstr_VP; // 3xN VPs as 3D direction 
  vector<T2> VPConstr_VPw; // 1xN VP weights
  vector<bool> VPConstr_fixVP;  // Fix VP or keep it as free variable

	vector<T1> Cams_noviews; // number of views per camera
	vector<vector<T2>> Cams_fc; // focal length, [fx,fy] or [fxy], 
	vector<vector<T2>> Cams_cc; // principal point, [cx,cy] or [] (=implicitly [0,0])
	vector<vector<T2>> Cams_kc; // radial distortion, noelem=0: No radial distortion, 
				    //                    noelem>0: Model after D.C. Brown, et al.'Decentering Distortion of lenses', 1966, 
				    //                              Order: 2nd, 4th symm. rad., 1st, 2nd tangent., 6th symm. rad dist.
				    //                              Simplified model: use only components 1, or 1+2, or 1+2+3+4, or 1+2+3+4+5 	
				    // Note: Skew is assumed to be zero
	vector<vector<bool>> Cams_fixInternals; // Fix focal length, principal point, radial distortion parameters
	vector<vector<T1>> Cams_viewids; // continuous enumeration of all views over all cameras, (starting with 1)
	vector<vector<valarray<T2>>> Cams_view_orien; // 3xN indep. orientation, 3x1 && noviews>1: shared orientation across all views
	vector<vector<valarray<T2>>> Cams_view_trans; // 3xN indep. translation, 3x1 && noviews>1: shared translation across all views,  supply t in world ref.frame (i.e. cam.center in world = t)
  vector<vector<vector<T2>>> Cams_smootherM; // smoother matrix; 
        // n x n matrix, n=noviews:   Full smoother Matrix for given p, Hastie&Tibshirani 2009, eq 5.17, PASS WITH SMOOTHING VECTORS AS COLUMN VECTORS ! (i.e. TRANSPOSED smoother matrix)
        // n x 1 column vector,n odd: Same smoothing for all views, n odd and n <= noviews , middle element will be smoothing center, vector will be renormalized to 1, high n destroys matrix sparsity in Jacobian ! PASS AS COLUMN VECTOR!
        // 0x0:                       No smoothing 
  
  vector<vector<vector<T1>>>         Cams_VP_ID;   // For every camera, and every view, list all visible VPs  1xN. Leave empty if unused.
  vector<vector<vector<vector<T2>>>> Cams_VPlines; // For all N visible VP in Cams_VPlines list endpoints of M consistent 2D lines (4xM). . Leave empty if unused.
  
	
	vector<vector<bool>> Cams_fixOrientation; // 1xN or 1x1, Fix camera orientation or keep it as free variable(s) 
	vector<vector<bool>> Cams_fixTranslation; // 1xN or 1x1, Fix camera translation or keep it as free variable(s)
	vector<vector<T1>> Cams_OnPlane; // Camera (with all views) lies on plane id, set to 0 if camera is not on plane, possibly multiple independent planes
	vector<vector<T2>> Cams_Weight; // Camera weight, 1xN or 1x1 for uniform weight

	vector<vector<T1>> Cams_reproj_view; //  [View_Idx1, View_Idx2] x K: Link camera center of view View_Idx1 to reprojection of camera center into view View_Idx2
	vector<vector<T2>> Cams_reproj_pos; // [X,Y] x K:  Position of reprojected camera center of view View_Idx1 in view view View_Idx2
	vector<T2> Cams_reproj_weight; // 1 x K:  Weight of reprojection in BA

	vector<valarray<T2>> Pt3d; // 3xM: 3D points, 
	vector<bool> PtFixPosition; // 1xM: Fix 3d position of point. If empty assume no point is fixed.
	vector<T2> PtWeight;  // Weights for point-wise reprojection error. If empty assume uniform weight of 1 
	vector<vector<T1>> PtOnPlane; // 1xM: Cell array, associate point with plane (possibly multiple planes)
	vector<vector<T1>> PtReprojView; // 1xM Vector of 1xO_i view indices, indicates view into which point is projected
	vector<vector<vector<T2>>> PtReprojPos;  // 1xM Vector of 2xO_i [x,y] position, indicates reproj. position in view, 
  
  T2 Residual;
};



#endif // _CeresBA_DataStruct__
