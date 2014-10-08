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


#include "BAdjust.h"
#include "ceresBA_datastruct.h"
#include "rotutil.h"

#define nullptr NULL


using namespace std;


template <class T1, class T2>
BAdjust<T1, T2>::BAdjust(string filename) {
	cout << "Load data from file: " << filename << endl;

	vector<T2> tmp3d(3);
	vector<T2> tmp2d_T2(2);
	vector<T1> tmp2d_T1(2);
	//T1 tmpT1;
	T1 no3Dpoints, nocams, nocamreproj, noplanes;
	//T2 tmpT2;  
	da = new BAdata<T1,T2>();

	ifstream infile(filename);
	string line;

	getline(infile, line); stringstream lineStream(line);
	lineStream >> da->Formatversion;
	
	getline(infile, line); lineStream.str(line);  lineStream.clear();
	lineStream >> da->BAopt_nocores >> da->BAopt_timeout_sec >> da->BAopt_maxiter;

	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> no3Dpoints;
	//vector<T1> da->noreproj((istream_iterator<T1>(lineStream)), istream_iterator<T1>());

	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> nocams;
	copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->Cams_noviews));
	da->Cams_fc.resize(nocams);
	da->Cams_kc.resize(nocams);
	da->Cams_cc.resize(nocams);
	da->Cams_fixInternals.resize(nocams);
	da->Cams_viewids.resize(nocams);
	da->Cams_view_orien.resize(nocams);
	da->Cams_fixOrientation.resize(nocams);	
	da->Cams_view_trans.resize(nocams);
	da->Cams_fixTranslation.resize(nocams);
	da->Cams_OnPlane.resize(nocams);
	da->Cams_Weight.resize(nocams);
	da->Cams_camcenter.resize(nocams);
	
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> nocamreproj;
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> noplanes;

	
	/*getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> da->GlobConstr_fixGroundplane;
	//if (! lineStream.eof()) {
	  //lineStream >> tmp3d[0] >> tmp3d[1] >> tmp3d[2];
	  //da->GlobConstr_groundplane = tmp3d;
	  copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->GlobConstr_groundplane));
	//}*/
	da->PlaneConstr_Planes.resize(noplanes);
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	for (T1 i = 0; i<noplanes; ++i) {
	  valarray<T2> tmp3dar(3);
	  lineStream >> tmp3dar[0] >> tmp3dar[1] >> tmp3dar[2];
	  da->PlaneConstr_Planes[i] = tmp3dar;
	}
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->PlaneConstr_PlanesW));
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->PlaneConstr_fixPlane));
	
	
	for (T1 i = 0; i<nocams; ++i) {
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_fc[i]));
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_cc[i]));
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_kc[i]));	    
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->Cams_fixInternals[i]));
	    
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    vector<T2> tmp_trans;
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(tmp_trans)); 
	    da->Cams_view_orien[i].resize(tmp_trans.size()/3);
	    for (T1 j=0; j<da->Cams_view_orien[i].size; ++j)
	      copy(tmp_trans.begin()+j*3, tmp_trans.begin()+(j+1)*3, back_inserter(da->Cams_view_orien[i][j])); 
    
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->Cams_fixOrientation[i]));

	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    vector<T2> tmp_oriens;
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(tmp_oriens)); 
	    da->Cams_view_trans[i].resize(tmp_oriens.size()/3);
	    for (T1 j=0; j<da->Cams_view_orien[i].size; ++j)
	      copy(tmp_oriens.begin()+j*3, tmp_oriens.begin()+(j+1)*3, back_inserter(da->Cams_view_orien[i][j])); 

	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->Cams_fixTranslation[i]));
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->Cams_viewids[i])); 
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->Cams_OnPlane[i]));
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_Weight[i])); 
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->Cams_camcenter[i])); 
	    
	    //copy(da->Cams_viewids[i].begin(), da->Cams_viewids[i].end(), ostream_iterator<T1>(cout, " "));
	    //cout << endl;
	}

	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	da->Cams_reproj_view.resize(nocamreproj);
	for (T1 i = 0; i<nocamreproj; ++i) {
	  lineStream >> tmp2d_T1[0] >> tmp2d_T1[1];
	  da->Cams_reproj_view[i] = tmp2d_T1;
	}
	getline(infile, line); lineStream.str(line);  lineStream.clear();		
	da->Cams_reproj_pos.resize(nocamreproj);
	for (T1 i = 0; i<nocamreproj; ++i) {
	  lineStream >> tmp2d_T2[0] >> tmp2d_T2[1];
	  da->Cams_reproj_pos[i] = tmp2d_T2;
	}
	getline(infile, line); lineStream.str(line);  lineStream.clear();		
	copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_reproj_weight)); 
	
	da->Pt3d.resize(no3Dpoints);
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	for (T1 i = 0; i<no3Dpoints; ++i) {
	  lineStream >> tmp3d[0] >> tmp3d[1] >> tmp3d[2];
	  da->Pt3d[i] = tmp3d;
	  //opy(da->Pt3d[i].begin(), da->Pt3d[i].end(), ostream_iterator<T2>(cout, " "));
	  //cout << endl;
	}
		
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->PtFixPosition)); 
	    
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->PtWeight)); 	

	getline(infile, line); lineStream.str(line);  lineStream.clear();
	T1 nopointsonplane;
	da->PtOnPlane.resize(no3Dpoints);
	lineStream >> nopointsonplane;
	for (T1 i = 0; i<nopointsonplane; ++i) {
	  getline(infile, line); lineStream.str(line);  lineStream.clear();	
	  int ptid;
	  lineStream >> ptid;
	  copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->PtOnPlane[ptid-1])); 
	}
	

	da->PtReprojView.resize(no3Dpoints);
	da->PtReprojPos.resize(no3Dpoints);
	for (T1 i = 0; i<no3Dpoints; ++i) {
	  getline(infile, line); lineStream.str(line);  lineStream.clear();
	  copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->PtReprojView[i])); 
	  getline(infile, line); lineStream.str(line);  lineStream.clear();
	  da->PtReprojPos[i].resize(da->PtReprojView[i].size());
	  for (T1 j=0; j < da->PtReprojView[i].size(); ++j)
	  {
	    lineStream >> tmp2d_T2[0] >> tmp2d_T2[1];
	    da->PtReprojPos[i][j] = tmp2d_T2;
	  }  
	  //copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->PtReprojPos[i])); 
	}
	
	infile.close();
	
	/*cout << "Finished loading " << nocams << " cameras with ";	
	for (int i = 0; i<nocams; ++i) cout << da->Cams_noviews[i] << "+";
	cout << "=" <<  accumulate(da->Cams_noviews.begin(),da->Cams_noviews.end(),0);
	cout << " views and " << no3Dpoints << " 3D points in ";
	cout << accumulate(da->PtReprojView.begin(),da->PtReprojView.end(), 0, [](int i, vector<T1> d) { return i+d.size(); });
	cout << " reprojections" << endl; */
	
	CheckErrors();
};


template <class T1, class T2>
BAdjust<T1, T2>::BAdjust(int nrhs, const mxArray **prhs, mxClassID T1mexClass_in, mxClassID T2mexClass_in) {
  T1mexclass = T1mexClass_in;
  T2mexclass = T2mexClass_in;
  
  const char * mexfieldnames_root[] = {"formatversion", "BAopt", "planeconstraint", "cams", "cam_reproj", "points"};
  const char * mexfieldnames_BAopt[] = {"cores", "timeout_sec", "maxiter"};
  const char * mexfieldnames_planeconstraints[] = {"plane", "planew", "fixPlane"};
  const char * mexfieldnames_cams[] = {"noviews", "fc", "cc", "kc", "views_orien", "views_trans", "viewids", "fixInternals", "fixOrientation", "fixTranslation", "OnPlane", "camweight", "camcenter"};
  const char * mexfieldnames_cam_reproj[] = {"view", "pos", "weight"};
  const char * mexfieldnames_points[] = {"pt3d", "OnPlane", "fixPosition", "pointw", "reproj_view", "reproj_pos"};
  
  da = new BAdata<T1,T2>();
  T1* dataT1;
  T2* dataT2;
  bool* dataBool;
  
  long int fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[0]);  
  mxArray *tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  mxArray *tmp3, *tmp4, *tmp5;
  dataT2  = (T2*)mxGetData(tmp);
  da->Formatversion = dataT2[0];

  // Parse BAopt
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[1]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[0]);
  mxArray *tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT1  = (T1*)mxGetData(tmp2);
  da->BAopt_nocores = dataT1[0];
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT1  = (T1*) mxGetData(tmp2);
  da->BAopt_timeout_sec = dataT1[0];

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_BAopt[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT1  = (T1*) mxGetData(tmp2);
  da->BAopt_maxiter = dataT1[0];
  
  // Parse plane constraints
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[2]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum);
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_planeconstraints[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  long int numelem = mxGetNumberOfElements(tmp2);
  copy(dataT2, dataT2+numelem, back_inserter(da->PlaneConstr_PlanesW));

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_planeconstraints[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataBool  = (bool*) mxGetData(tmp2);
  da->PlaneConstr_fixPlane.resize(numelem);
  for (T1 i = 0; i<numelem; ++i)
    da->PlaneConstr_fixPlane[i] = dataBool[i];  // use this to avoid gcc 4.4 operator= overload error
  //copy(dataBool, dataBool+numelem, back_inserter(da->PlaneConstr_fixPlane));

  da->PlaneConstr_Planes.resize(numelem);
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_planeconstraints[0]);  
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);  
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<numelem; ++i)
  {
    //da->PlaneConstr_Planes[i] = valarray<T2> {dataT2[i*3],dataT2[i*3+1],dataT2[i*3+2]};
    da->PlaneConstr_Planes[i].resize(3);
    da->PlaneConstr_Planes[i][0] = dataT2[i*3];
    da->PlaneConstr_Planes[i][1] = dataT2[i*3+1];
    da->PlaneConstr_Planes[i][2] = dataT2[i*3+2];
  }
  //cout << da->PlaneConstr_Planes[0][0] << da->PlaneConstr_Planes[0][1] << da->PlaneConstr_Planes[0][2] << endl;
  //copy(da->GlobConstr_groundplane.begin(), da->GlobConstr_groundplane.end(), ostream_iterator<T2>(cout, " "));
 
  // Parse Cameras
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[3]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum); 
  
  T1 nocams = mxGetNumberOfElements(tmp);
  da->Cams_fc.resize(nocams);
  da->Cams_cc.resize(nocams);
  da->Cams_kc.resize(nocams);
  da->Cams_fixInternals.resize(nocams);
  da->Cams_viewids.resize(nocams);
  da->Cams_view_orien.resize(nocams);
  da->Cams_view_trans.resize(nocams);
  da->Cams_fixOrientation.resize(nocams);
  da->Cams_fixTranslation.resize(nocams);
  da->Cams_OnPlane.resize(nocams);
  da->Cams_Weight.resize(nocams);
  da->Cams_camcenter.resize(nocams);
  
  for (T1 i = 0; i<nocams; ++i)
  {
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[0]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT1  = (T1*) mxGetData(tmp2);
    da->Cams_noviews.push_back(dataT1[0]);

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[1]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(da->Cams_fc[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[2]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(da->Cams_cc[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[3]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(da->Cams_kc[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[4]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    da->Cams_view_orien[i].resize(numelem/3);
    for (T1 j = 0; j<da->Cams_view_orien[i].size(); ++j)
    {
      da->Cams_view_orien[i][j].resize(3);
      da->Cams_view_orien[i][j][0] = dataT2[j*3];
      da->Cams_view_orien[i][j][1] = dataT2[j*3+1];
      da->Cams_view_orien[i][j][2] = dataT2[j*3+2];
    }      
    //{
      //da->Cams_view_orien[i][j] = valarray<T2> {dataT2[j*3],dataT2[j*3+1],dataT2[j*3+2]};
      //dataT2[j*3] = 
      //copy(dataT2+j*3, dataT2+(j+1)*3, back_inserter(da->Cams_view_orien[i][j]));
    //}

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[5]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    da->Cams_view_trans[i].resize(numelem/3);
    for (T1 j = 0; j<da->Cams_view_trans[i].size(); ++j)
    {
      da->Cams_view_trans[i][j].resize(3);
      da->Cams_view_trans[i][j][0] = dataT2[j*3];
      da->Cams_view_trans[i][j][1] = dataT2[j*3+1];
      da->Cams_view_trans[i][j][2] = dataT2[j*3+2];
    }        
      //da->Cams_view_trans[i][j] = valarray<T2> {dataT2[j*3],dataT2[j*3+1],dataT2[j*3+2]};      
      //copy(dataT2+j*3, dataT2+(j+1)*3, back_inserter(da->Cams_view_trans[i][j]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[6]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT1  = (T1*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT1, dataT1+numelem, back_inserter(da->Cams_viewids[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[7]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    da->Cams_fixInternals[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      da->Cams_fixInternals[i][j] = dataBool[j];  // use this to avoid gcc 4.4 operator= overload error
    //copy(dataBool, dataBool+numelem, back_inserter(da->Cams_fixInternals[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[8]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    da->Cams_fixOrientation[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      da->Cams_fixOrientation[i][j] = dataBool[j];  // use this to avoid gcc 4.4 operator= overload error
    //copy(dataBool, dataBool+numelem, back_inserter(da->Cams_fixOrientation[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[9]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    da->Cams_fixTranslation[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      da->Cams_fixTranslation[i][j] = dataBool[j];  // use this to avoid gcc 4.4 operator= overload error
    //copy(dataBool, dataBool+numelem, back_inserter(da->Cams_fixTranslation[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[10]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT1  = (T1*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT1, dataT1+numelem, back_inserter(da->Cams_OnPlane[i]));

    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[11]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataT2  = (T2*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    copy(dataT2, dataT2+numelem, back_inserter(da->Cams_Weight[i]));
    
    fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cams[12]);
    tmp2 = mxGetFieldByNumber(tmp, i, fieldnum);
    dataBool = (bool*) mxGetData(tmp2);
    numelem = mxGetNumberOfElements(tmp2);
    da->Cams_camcenter[i].resize(numelem);
    for (T1 j = 0; j<numelem; ++j)
      da->Cams_camcenter[i][j] = dataBool[j];  // use this to avoid gcc 4.4 operator= overload error
    //copy(dataBool, dataBool+numelem, back_inserter(da->Cams_camcenter[i]));

  }

  // Parse mutual camera visibility
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[4]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum); 
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cam_reproj[0]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  T1 nomutualvis = mxGetNumberOfElements(tmp2) / 2;
  
  da->Cams_reproj_view.resize(nomutualvis);
  da->Cams_reproj_pos.resize(nomutualvis);

  dataT1  = (T1*) mxGetPr(tmp2);
  for (T1 i = 0; i<nomutualvis; ++i)
  {
    copy(dataT1, dataT1+2, back_inserter(da->Cams_reproj_view[i]));
    dataT1 = dataT1 + 2;
  }

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cam_reproj[1]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<nomutualvis; ++i)
  {
    copy(dataT2, dataT2+2, back_inserter(da->Cams_reproj_pos[i]));
    dataT2 = dataT2 + 2;
  }
  
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_cam_reproj[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2  = (T2*) mxGetData(tmp2);
  copy(dataT2, dataT2+nomutualvis, back_inserter(da->Cams_reproj_weight));

  // Parse points
  fieldnum = mxGetFieldNumber(prhs[0], mexfieldnames_root[5]);
  tmp = mxGetFieldByNumber(prhs[0], 0, fieldnum); 
  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[0]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  T1 nopts = mxGetNumberOfElements(tmp2) / 3;
  da->Pt3d.resize(nopts);
  da->PtReprojView.resize(nopts);
  da->PtReprojPos.resize(nopts);
  da->PtOnPlane.resize(nopts);  
  dataT2  = (T2*) mxGetPr(tmp2);
  for (T1 i = 0; i<nopts ; ++i) 
  {
    da->Pt3d[i].resize(3);
    da->Pt3d[i][0] = dataT2[i*3];
    da->Pt3d[i][1] = dataT2[i*3+1];
    da->Pt3d[i][2] = dataT2[i*3+2];
  }        
  //{
    //da->Pt3d[i] = valarray<T2> {dataT2[i*3],dataT2[i*3+1],dataT2[i*3+2]};    
    //copy(dataT2, dataT2+3, back_inserter(da->Pt3d[i]));
    //dataT2 = dataT2 + 3;
  //}
  

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[2]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataBool = (bool*) mxGetPr(tmp2);
  da->PtFixPosition.resize(nopts);
  for (T1 j = 0; j<nopts; ++j)
    da->PtFixPosition[j] = dataBool[j];  // use this to avoid gcc 4.4 operator= overload error
  //copy(dataBool, dataBool+nopts, back_inserter(da->PtFixPosition));

  fieldnum = mxGetFieldNumber(tmp, mexfieldnames_points[3]);
  tmp2 = mxGetFieldByNumber(tmp, 0, fieldnum);
  dataT2 = (T2*) mxGetPr(tmp2);
  copy(dataT2, dataT2+nopts, back_inserter(da->PtWeight));

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
    copy(dataT1, dataT1+numelem, back_inserter(da->PtReprojView[i]));

    tmp4 = mxGetCell(tmp3, i);
    dataT2 = (T2*) mxGetPr(tmp4);
    numelem = mxGetNumberOfElements(tmp4);
    da->PtReprojPos[i].resize(numelem/2);
    for (T1 j = 0; j<da->PtReprojPos[i].size(); ++j) 
      copy(dataT2+j*2, dataT2+(j+1)*2, back_inserter(da->PtReprojPos[i][j]));
	
    tmp4 = mxGetCell(tmp5, i);
    if (tmp4!=nullptr) 
    {
      numelem = mxGetNumberOfElements(tmp4);
      dataT1 = (T1*) mxGetPr(tmp4);
      copy(dataT1, dataT1+numelem, back_inserter(da->PtOnPlane[i]));
/*      
      for (T1 l = 0; l<da->PtOnPlane[i].size(); ++l) 
	cout << da->PtOnPlane[i][l] << endl;*/
    }
  }

  CheckErrors();
}

template <class T1, class T2>
void BAdjust<T1, T2>::WriteToMex(int* nlhs, mxArray **plhs) {
  const char * mexfieldnames_root[] = {"formatversion", "BAopt", "planeconstraint", "cams", "cam_reproj", "points"};
  const char * mexfieldnames_BAopt[] = {"cores", "timeout_sec", "maxiter"};
  const char * mexfieldnames_planeconstraints[] = {"plane", "planew", "fixPlane"};
  const char * mexfieldnames_cams[] = {"noviews", "fc", "cc", "kc", "views_orien", "views_trans", "viewids", "fixInternals", "fixOrientation", "fixTranslation", "OnPlane", "camweight", "camcenter"};
  const char * mexfieldnames_cam_reproj[] = {"view", "pos", "weight"};
  const char * mexfieldnames_points[] = {"pt3d", "OnPlane", "fixPosition", "pointw", "reproj_view", "reproj_pos"};

  
  mxArray *tmp, *tmp2, *tmp3, *tmp4, *tmp5, *tmp6;
    
  T1* dataT1;
  T2* dataT2;
  
  plhs[0] = mxCreateStructMatrix(1, 1, 6, mexfieldnames_root);
  (*nlhs)++;
  
  tmp = mxCreateNumericMatrix(1, 1, T2mexclass, mxREAL);
  dataT2 = (T2*)mxGetData(tmp);
  dataT2[0] = da->Formatversion;
  mxSetField(plhs[0], 0, mexfieldnames_root[0], tmp);
	     
  // Write BAopt
  tmp2 = mxCreateStructMatrix(1, 1, 3, mexfieldnames_BAopt);

      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = da->BAopt_nocores;
      mxSetField(tmp2, 0, mexfieldnames_BAopt[0], tmp);
      
      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = da->BAopt_timeout_sec;
      mxSetField(tmp2, 0, mexfieldnames_BAopt[1], tmp);
		
      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = da->BAopt_maxiter;
      mxSetField(tmp2, 0, mexfieldnames_BAopt[2], tmp);
      
  mxSetField(plhs[0], 0, mexfieldnames_root[1], tmp2);
  
  // Write plane constraints
  tmp2 = mxCreateStructMatrix(1, 1, 3, mexfieldnames_planeconstraints);
  
    tmp = mxCreateNumericMatrix(3, da->PlaneConstr_Planes.size(), T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    for (T1 i = 0; i<da->PlaneConstr_Planes.size() ; ++i) {
      dataT2[3*i+0] = da->PlaneConstr_Planes[i][0];
      dataT2[3*i+1] = da->PlaneConstr_Planes[i][1];
      dataT2[3*i+2] = da->PlaneConstr_Planes[i][2];
    }
    mxSetField(tmp2, 0, mexfieldnames_planeconstraints[0], tmp); 

    tmp = mxCreateNumericMatrix(1, da->PlaneConstr_Planes.size(), T2mexclass, mxREAL);
    copy(da->PlaneConstr_PlanesW.begin(), da->PlaneConstr_PlanesW.end(), mxGetPr(tmp));   
    mxSetField(tmp2, 0, mexfieldnames_planeconstraints[1], tmp); 
    
    
    tmp = mxCreateNumericMatrix(1, da->PlaneConstr_Planes.size(), mxUINT8_CLASS, mxREAL);
    uint8_t *uint8data = (uint8_t*)mxGetPr(tmp);
    //for (T1 i = 0; i<da->PlaneConstr_Planes.size() ; ++i)
      //uint8data[i] = da->PlaneConstr_fixPlane[i];
    copy(da->PlaneConstr_fixPlane.begin(), da->PlaneConstr_fixPlane.end(), uint8data); 
    mxSetField(tmp2, 0, mexfieldnames_planeconstraints[2], tmp);
    
  mxSetField(plhs[0], 0, mexfieldnames_root[2], tmp2); 

  // Write Cameras
  T1 nocams = da->Cams_noviews.size();
  tmp2 = mxCreateStructMatrix(1, nocams , 13, mexfieldnames_cams);
  for (T1 i = 0; i<nocams; ++i) {
      tmp = mxCreateNumericMatrix(1, 1, T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      dataT1[0] = da->Cams_noviews[i];
      mxSetField(tmp2, i, mexfieldnames_cams[0], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_fc[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(da->Cams_fc[i].begin(), da->Cams_fc[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[1], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_cc[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(da->Cams_cc[i].begin(), da->Cams_cc[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[2], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_kc[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(da->Cams_kc[i].begin(), da->Cams_kc[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[3], tmp);

      tmp = mxCreateNumericMatrix(3, da->Cams_view_orien[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      for (T1 j=0; j< da->Cams_view_orien[i].size(); ++j)
      {
	dataT2[j*3+0] = da->Cams_view_orien[i][j][0];
	dataT2[j*3+1] = da->Cams_view_orien[i][j][1];
	dataT2[j*3+2] = da->Cams_view_orien[i][j][2];
      }
	//copy(da->Cams_view_orien[i][j].begin(), da->Cams_view_orien[i][j].end(), dataT2+j*3); 	
      mxSetField(tmp2, i, mexfieldnames_cams[4], tmp);

      tmp = mxCreateNumericMatrix(3, da->Cams_view_trans[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      for (T1 j=0; j< da->Cams_view_trans[i].size(); ++j)
      {
	dataT2[j*3+0] = da->Cams_view_trans[i][j][0];
	dataT2[j*3+1] = da->Cams_view_trans[i][j][1];
	dataT2[j*3+2] = da->Cams_view_trans[i][j][2];
      }
	//copy(da->Cams_view_trans[i][j].begin(), da->Cams_view_trans[i][j].end(), dataT2+j*3); 	
      mxSetField(tmp2, i, mexfieldnames_cams[5], tmp);
      
      tmp = mxCreateNumericMatrix(1, da->Cams_viewids[i].size(), T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      copy(da->Cams_viewids[i].begin(), da->Cams_viewids[i].end(), dataT1); 
      mxSetField(tmp2, i, mexfieldnames_cams[6], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_fixInternals[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(da->Cams_fixInternals[i].begin(), da->Cams_fixInternals[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[7], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_fixOrientation[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(da->Cams_fixOrientation[i].begin(), da->Cams_fixOrientation[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[8], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_fixTranslation[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(da->Cams_fixTranslation[i].begin(), da->Cams_fixTranslation[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[9], tmp);

      tmp = mxCreateNumericMatrix(1, da->Cams_OnPlane[i].size(), T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetData(tmp);
      copy(da->Cams_OnPlane[i].begin(), da->Cams_OnPlane[i].end(), dataT1); 
      mxSetField(tmp2, i, mexfieldnames_cams[10], tmp);
      
      tmp = mxCreateNumericMatrix(1, da->Cams_Weight[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetData(tmp);
      copy(da->Cams_Weight[i].begin(), da->Cams_Weight[i].end(), dataT2); 
      mxSetField(tmp2, i, mexfieldnames_cams[11], tmp);
      
      tmp = mxCreateNumericMatrix(1, da->Cams_camcenter[i].size(), mxUINT8_CLASS, mxREAL);
      uint8data = (uint8_t*)mxGetPr(tmp);
      copy(da->Cams_camcenter[i].begin(), da->Cams_camcenter[i].end(), uint8data); 
      mxSetField(tmp2, i, mexfieldnames_cams[12], tmp);
      
  }
  mxSetField(plhs[0], 0, mexfieldnames_root[3], tmp2);
  
  // Write mutual camera visibility
  tmp2 = mxCreateStructMatrix(1, 1, 3, mexfieldnames_cam_reproj);
  T1 nomutualvis = da->Cams_reproj_view.size();

  tmp = mxCreateNumericMatrix(2, nomutualvis, T1mexclass, mxREAL);
    dataT1 = (T1*)mxGetPr(tmp);
    for (T1 i = 0; i< nomutualvis  ; ++i) {
      copy(da->Cams_reproj_view[i].begin(), da->Cams_reproj_view[i].end(), dataT1);       
      dataT1 = dataT1 + 2;      
    }
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[0], tmp); 

    tmp = mxCreateNumericMatrix(2, nomutualvis , T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    for (T1 i = 0; i< nomutualvis  ; ++i) {
      copy(da->Cams_reproj_pos[i].begin(), da->Cams_reproj_pos[i].end(), dataT2);       
      dataT2 = dataT2 + 2;      
    }
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[1], tmp); 

    tmp = mxCreateNumericMatrix(1, nomutualvis, T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    copy(da->Cams_reproj_weight.begin(), da->Cams_reproj_weight.end(), dataT2); 
    mxSetField(tmp2, 0, mexfieldnames_cam_reproj[2], tmp); 
      
  mxSetField(plhs[0], 0, mexfieldnames_root[4], tmp2);

  // Write points
  tmp2 = mxCreateStructMatrix(1, 1, 6, mexfieldnames_points);
  T1 nopts = da->Pt3d.size();

    tmp = mxCreateNumericMatrix(1, nopts, mxUINT8_CLASS, mxREAL);
    uint8data= (uint8_t*)mxGetPr(tmp);
    copy(da->PtFixPosition.begin(), da->PtFixPosition.end(), uint8data); 
    mxSetField(tmp2, 0, mexfieldnames_points[2], tmp); 

    tmp = mxCreateNumericMatrix(1, nopts, T2mexclass, mxREAL);
    dataT2 = (T2*)mxGetPr(tmp);
    copy(da->PtWeight.begin(), da->PtWeight.end(), dataT2); 
    mxSetField(tmp2, 0, mexfieldnames_points[3], tmp); 

    tmp = mxCreateNumericMatrix(3, nopts, T2mexclass, mxREAL);
    tmp3 = mxCreateCellMatrix(1, nopts);
    tmp4 = mxCreateCellMatrix(1, nopts);
    tmp5 = mxCreateCellMatrix(1, nopts);

    for (T1 i = 0; i< nopts; ++i)
    {
      dataT2 = (T2*)mxGetPr(tmp);
      dataT2[i*3+0] = da->Pt3d[i][0];
      dataT2[i*3+1] = da->Pt3d[i][1];
      dataT2[i*3+2] = da->Pt3d[i][2];
      //copy(da->Pt3d[i].begin(), da->Pt3d[i].end(), dataT2+i*3); 

      if (da->PtOnPlane[i].size() > 0)
      {
	tmp6 = mxCreateNumericMatrix(1, da->PtOnPlane[i].size(), T1mexclass, mxREAL);
	dataT1 = (T1*)mxGetPr(tmp6);
	copy(da->PtOnPlane[i].begin(), da->PtOnPlane[i].end(), dataT1); 
	mxSetCell(tmp3, i, tmp6);
      }
      
      tmp6 = mxCreateNumericMatrix(1, da->PtReprojView[i].size(), T1mexclass, mxREAL);
      dataT1 = (T1*)mxGetPr(tmp6);
      copy(da->PtReprojView[i].begin(), da->PtReprojView[i].end(), dataT1); 
      mxSetCell(tmp4, i, tmp6);

      tmp6 = mxCreateNumericMatrix(2, da->PtReprojPos[i].size(), T2mexclass, mxREAL);
      dataT2 = (T2*)mxGetPr(tmp6);
      for (T1 j=0; j< da->PtReprojPos[i].size(); ++j)
	copy(da->PtReprojPos[i][j].begin(), da->PtReprojPos[i][j].end(), dataT2+j*2); 	
      mxSetCell(tmp5, i, tmp6);

    }
    mxSetField(tmp2, 0, mexfieldnames_points[0], tmp);
    mxSetField(tmp2, 0, mexfieldnames_points[1], tmp3);
    mxSetField(tmp2, 0, mexfieldnames_points[4], tmp4); 
    mxSetField(tmp2, 0, mexfieldnames_points[5], tmp5);
  
  mxSetField(plhs[0], 0, mexfieldnames_root[5], tmp2); 
};


template <class T1, class T2>
BAdjust<T1, T2>::BAdjust(BAdata<T1,T2>* da_in){
  da = da_in;
  
  CheckErrors();
};


template <class T1, class T2>
void BAdjust<T1, T2>::WriteToFile(string filename)
{
	T1 nocams = da->Cams_noviews.size();
	T1 no3Dpoints = da->Pt3d.size();
	T1 noreprojcams = da->Cams_reproj_view.size();
	T1 noplanes = da->PlaneConstr_Planes.size();	
	
  	cout << "Save data to file: " << filename << endl;

	ofstream outfile(filename, ofstream::out);
	outfile << setprecision(20);
	
	outfile << da->Formatversion <<  endl;
	outfile << da->BAopt_nocores <<  " " << da->BAopt_timeout_sec << " " << da->BAopt_maxiter << endl;

	outfile << no3Dpoints << " ";
	for (T1 i=0; i<no3Dpoints; ++i)
	  outfile << da->PtReprojView[i].size() << " ";
	outfile << endl;
	
	outfile << nocams << " ";
	copy(da->Cams_noviews.begin(), da->Cams_noviews.end(), ostream_iterator<T1>(outfile, " "));  outfile << endl;
	
	outfile << noreprojcams << endl;
	outfile << noplanes << endl;	
	
	for (T1 i=0; i<noplanes; ++i)
	  outfile << da->PlaneConstr_Planes[i][0] << " " << da->PlaneConstr_Planes[i][1] << " " << da->PlaneConstr_Planes[i][2] << " ";
	outfile << endl;	
	copy(da->PlaneConstr_PlanesW.begin(), da->PlaneConstr_PlanesW.end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;	
	copy(da->PlaneConstr_fixPlane.begin(), da->PlaneConstr_fixPlane.end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;	
	
 	//outfile << da->GlobConstr_fixGroundplane << " ";
	//copy(da->GlobConstr_groundplane.begin(), da->GlobConstr_groundplane.end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
	
	for (T1 i=0; i<nocams; ++i) {
	  	copy(da->Cams_fc[i].begin(), da->Cams_fc[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
		copy(da->Cams_cc[i].begin(), da->Cams_cc[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;				
		copy(da->Cams_kc[i].begin(), da->Cams_kc[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
		copy(da->Cams_fixInternals[i].begin(), da->Cams_fixInternals[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;		
		
		for (T1 j=0; j<da->Cams_view_orien[i].size(); ++j)
		  outfile << da->Cams_view_orien[i][j][0] << " " << da->Cams_view_orien[i][j][1] << " " << da->Cams_view_orien[i][j][2] << " ";
		outfile << endl;
		
		//copy(da->Cams_view_orien[i].begin(), da->Cams_view_orien[i].end(), ostream_iterator<T2>(outfile, " "));  
		//outfile << endl;
		
		copy(da->Cams_fixOrientation[i].begin(), da->Cams_fixOrientation[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;
		
		for (T1 j=0; j<da->Cams_view_trans[i].size(); ++j)
		  outfile << da->Cams_view_trans[i][j][0] << " " << da->Cams_view_trans[i][j][1] << " " << da->Cams_view_trans[i][j][2] << " ";
		outfile << endl;

		//copy(da->Cams_view_trans[i].begin(), da->Cams_view_trans[i].end(), ostream_iterator<T2>(outfile, " "));  
		//outfile << endl;
		
		copy(da->Cams_fixTranslation[i].begin(), da->Cams_fixTranslation[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;		
		copy(da->Cams_viewids[i].begin(), da->Cams_viewids[i].end(), ostream_iterator<T1>(outfile, " "));  outfile << endl;
		copy(da->Cams_OnPlane[i].begin(), da->Cams_OnPlane[i].end(), ostream_iterator<T1>(outfile, " "));  outfile << endl;		
		copy(da->Cams_Weight[i].begin(), da->Cams_Weight[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
		copy(da->Cams_camcenter[i].begin(), da->Cams_camcenter[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;				
	}
	
	for (T1 i=0; i<noreprojcams; ++i)
	    copy(da->Cams_reproj_view[i].begin(), da->Cams_reproj_view[i].end(), ostream_iterator<T1>(outfile, " "));
	outfile << endl;
	for (T1 i=0; i<noreprojcams; ++i)
	    copy(da->Cams_reproj_pos[i].begin(), da->Cams_reproj_pos[i].end(), ostream_iterator<T2>(outfile, " "));
	outfile << endl;
	copy(da->Cams_reproj_weight.begin(), da->Cams_reproj_weight.end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
		
	for (T1 i=0; i<no3Dpoints; ++i)
	  outfile << da->Pt3d[i][0] << " " << da->Pt3d[i][1] << " " << da->Pt3d[i][2] << " ";
	outfile << endl;

	copy(da->PtFixPosition.begin(), da->PtFixPosition.end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;	
	copy(da->PtWeight.begin(), da->PtWeight.end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;	
	
	T1 nopointsonplane=0;
	for (T1 i=0; i<no3Dpoints; ++i)
	  if (da->PtOnPlane[i].size() > 0)
	  //{
	    nopointsonplane++;
	   // cout << i << endl;
	 // }
	  	//cout << nopointsonplane << endl;
		
		
	outfile << nopointsonplane << endl;  
	for (T1 i=0; i<no3Dpoints; ++i)
	{
	  if (da->PtOnPlane[i].size() > 0)
	  {
	    outfile << (i+1) << " ";
	    copy(da->PtOnPlane[i].begin(), da->PtOnPlane[i].end(), ostream_iterator<T1>(outfile, " "));  
	    outfile << endl;
	  }
	}
	
	for (T1 i=0; i<no3Dpoints; ++i) {
	  copy(da->PtReprojView[i].begin(), da->PtReprojView[i].end(), ostream_iterator<T1>(outfile, " "));  outfile << endl;
	  for (T1 j=0; j< da->PtReprojPos[i].size(); ++j)
	    outfile << da->PtReprojPos[i][j][0] << " " << da->PtReprojPos[i][j][1] << " ";
	  outfile << endl;
	  //copy(da->PtReprojPos[i].begin(), da->PtReprojPos[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl; 	 
	}

	outfile.close(); 
};


template <class T1, class T2>
void BAdjust<T1, T2>::CheckErrors()
{
  
  // Check cameras
  // Check if a camera has shared translation over all views, and camera center set to -R'*t instead of t
  T1 nocams = da->Cams_noviews.size();
    
  for(T1 i=0; i <nocams ; ++i)
  {
    if (da->Cams_view_trans[i].size() == 1 &&  da->Cams_noviews[i] > 1 && da->Cams_camcenter[i][0]==1)
      cout << "Warning: Camera " << i+1 << " has shared translation over multiple views, but camera center is not set in world frame. Set camcenter to 0." << endl;
  }
};




template <class T1, class T2>
void BAdjust<T1, T2>::NormalizeData()
{
  vector<vector<T2>> varvect(3); 
  bool ct;
  
  for (T1 i=0; i< da->Cams_noviews.size(); ++i)
  //for (T1 i=0; i< 1; ++i)
  {
    for (T1 j=0; j< da->Cams_noviews[i]; ++j)
    //for (T1 j=0; j< 1; ++j)
    {
      valarray<T2> tr(T2(0),3);
      valarray<T2> ori(T2(0),3);

      if (da->Cams_view_trans[i].size()==1) // when only one set of parameters
	tr = da->Cams_view_trans[i][0];
      else
	tr = da->Cams_view_trans[i][j];
	
      if (da->Cams_view_orien[i].size()==1) // when only one set of parameters
	ori = da->Cams_view_orien[i][0];
      else
	ori = da->Cams_view_orien[i][j];
	
      if (da->Cams_camcenter[i].size()==1) // when only one set of parameters
	ct = da->Cams_camcenter[i][0];
      else
	ct = da->Cams_camcenter[i][j];  
    
      //cout << ori[0] << " " << ori[1] << " " << ori[2] << endl;
      //cout << tr[0] << " " << tr[1] << " " << tr[2] << endl;      
      
      if (ct==1)
      {
	tr *= -1;
	ori *= -1;
	valarray<T2> res(T2(0),3);
	rotutil::RotatePoint(tr, ori, &res);
	tr = res;
      }
	
      varvect[0].push_back(tr[0]);
      varvect[1].push_back(tr[1]);
      varvect[2].push_back(tr[2]);
      //varvect.push_back(move(tr)); // expliclty use std::move() to invoke move semantics
    }
  }
  
  for (T1 i=0; i< da->Pt3d.size(); ++i)
  {
      varvect[0].push_back(da->Pt3d[i][0]);
      varvect[1].push_back(da->Pt3d[i][1]);
      varvect[2].push_back(da->Pt3d[i][2]);
  }
  
  median3d.resize(3);
  
  // Compute median
  for (T1 i=0; i< 3; ++i)
  {
    vector<T2> tmp = varvect[i];
    sort(tmp.begin(), tmp.end());
    median3d[i] = tmp[int(tmp.size()/2)];
  }
  
  //cout <<median3d[0] <<  " " <<median3d[1] <<  " " <<median3d[2] << endl;
  
  // Compute Variance of L1 distances of 
  vector<T2> tmp(varvect[0].size());
  for (T1 i=0; i< varvect[0].size(); ++i)
    tmp[i] = abs((varvect[0][i] - median3d[0])) +  abs((varvect[1][i] - median3d[1])) + abs((varvect[2][i] - median3d[2])); // L1 norm
  
  //cout << tmp[0] << tmp[1] << tmp[2] << endl; 
  T2 accum = inner_product( tmp.begin(), tmp.end(), tmp.begin(), 0 );
  stdev3d = sqrt(accum / (tmp.size() - 1));
  
  //T2 accum = 0.0;
  //for_each (tmp.begin(), tmp.end(), [&](const double d) { accum += d * d; });
  //stdev3d[0] = sqrt(accum / (tmp.size()-1));
  //cout << stdev3d << endl;

  // move cameras and points to zero median and unit isotropic stdev.
  for (T1 i=0; i< da->Pt3d.size(); ++i)
  {
    da->Pt3d[i] -= median3d;
    da->Pt3d[i] /= stdev3d;
  }
  
  for (T1 i=0; i< da->Cams_noviews.size(); ++i)
  {
    for (T1 j=0; j< da->Cams_view_trans[i].size(); ++j)
    {
      bool ct;
      if (da->Cams_camcenter[i].size()==1) // when only one set of parameters
        ct = da->Cams_camcenter[i][0];
      else
        ct = da->Cams_camcenter[i][j];  
      
      if (ct==0)
      {
        da->Cams_view_trans[i][j] -= median3d;
        da->Cams_view_trans[i][j] /= stdev3d;
      }
      else
      {
        valarray<T2> ori(T2(0),3);
        if (da->Cams_view_orien[i].size()==1) // when only one set of parameters
          ori = da->Cams_view_orien[i][0];
        else
          ori = da->Cams_view_orien[i][j];
            
        valarray<T2> tr = da->Cams_view_trans[i][j];
        
        // get camera center in world coordinates
        tr *= -1;
        ori *= -1; // invert rotation
        valarray<T2> res(T2(0),3);
        rotutil::RotatePoint(tr, ori, &res);
        tr = res;

        tr -= median3d;
        tr /= stdev3d;
        
        // convert back to camera centered translation
        tr *= -1;
        ori *= -1;  // undo inversion of rotation
        res.resize(3);
        rotutil::RotatePoint(tr, ori, &res);
        da->Cams_view_trans[i][j] = res;
      }
      
    }
  }
  
  // Normalize planes
  for (T1 i=0; i< da->PlaneConstr_Planes.size(); ++i)
  {
    T2 planenorms = sqrt(da->PlaneConstr_Planes[i][0]*da->PlaneConstr_Planes[i][0] +
			  da->PlaneConstr_Planes[i][1]*da->PlaneConstr_Planes[i][1] +
			  da->PlaneConstr_Planes[i][2]*da->PlaneConstr_Planes[i][2]);

    da->PlaneConstr_Planes[i][0] /= planenorms;
    da->PlaneConstr_Planes[i][1] /= planenorms;
    da->PlaneConstr_Planes[i][2] /= planenorms;

    
    planenorms -= rotutil::DotProduct(da->PlaneConstr_Planes[i], median3d);
    planenorms /= stdev3d;
    
    da->PlaneConstr_Planes[i][0] *= planenorms;
    da->PlaneConstr_Planes[i][1] *= planenorms;
    da->PlaneConstr_Planes[i][2] *= planenorms;
  }
      

  
};

template <class T1, class T2>
void BAdjust<T1, T2>::UnnormalizeData()
{
  for (T1 i=0; i< da->Pt3d.size(); ++i)
  {
    da->Pt3d[i] *= stdev3d;
    da->Pt3d[i] += median3d;
  }
  
  for (T1 i=0; i< da->Cams_noviews.size(); ++i)
  {
    for (T1 j=0; j< da->Cams_view_trans[i].size(); ++j)
    {
      bool ct;
      if (da->Cams_camcenter[i].size()==1) // when only one set of parameters
        ct = da->Cams_camcenter[i][0];
      else
        ct = da->Cams_camcenter[i][j];  
      
      if (ct==0)
      {
        da->Cams_view_trans[i][j] *= stdev3d;
        da->Cams_view_trans[i][j] += median3d;
      }
      else
      {
        valarray<T2> ori(T2(0),3);
        if (da->Cams_view_orien[i].size()==1) // when only one set of parameters
          ori = da->Cams_view_orien[i][0];
        else
          ori = da->Cams_view_orien[i][j];
      
          valarray<T2> tr = da->Cams_view_trans[i][j];
          
          // get camera center in world coordinates
          tr *= -1;
          ori *= -1; // invert rotation
          valarray<T2> res(T2(0),3);
          rotutil::RotatePoint(tr, ori, &res);
          tr = res;

          tr *= stdev3d;
          tr += median3d;

          // convert back to camera centered translation
          tr *= -1;
          ori *= -1;  // undo inversion of rotation
          res.resize(3);
          rotutil::RotatePoint(tr, ori, &res);
          da->Cams_view_trans[i][j] = res;
      }
    }
  }
  
  // Unnormalize planes
  for (T1 i=0; i< da->PlaneConstr_Planes.size(); ++i)
  {
    T2 planenorms = sqrt(da->PlaneConstr_Planes[i][0]*da->PlaneConstr_Planes[i][0] +
			  da->PlaneConstr_Planes[i][1]*da->PlaneConstr_Planes[i][1] +
			  da->PlaneConstr_Planes[i][2]*da->PlaneConstr_Planes[i][2]);

    da->PlaneConstr_Planes[i][0] /= planenorms;
    da->PlaneConstr_Planes[i][1] /= planenorms;
    da->PlaneConstr_Planes[i][2] /= planenorms;

    planenorms *= stdev3d;    
    planenorms += rotutil::DotProduct(da->PlaneConstr_Planes[i], median3d);
    
    da->PlaneConstr_Planes[i][0] *= planenorms;
    da->PlaneConstr_Planes[i][1] *= planenorms;
    da->PlaneConstr_Planes[i][2] *= planenorms;
  }  
}