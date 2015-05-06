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


#include "Parser.h"
#include "ceresBA_datastruct.h"
#include "resutil.h"

#define nullptr NULL


using namespace std;


template <class T1, class T2>
Parser<T1, T2>::Parser(string filename) {
	cout << "Load data from file: " << filename << endl;

	vector<T2> tmp3d(3);
  valarray<T2> tmp3d_valarr(3);
	vector<T2> tmp2d_T2(2);
	vector<T1> tmp2d_T1(2);
  
	T1 no3Dpoints, nocams, nocamreproj, noplanes, novps, vpid, nolines, cnt;

  da = new BAdata<T1,T2>();

	ifstream infile(filename);
	string line;

	getline(infile, line); stringstream lineStream(line);
	lineStream >> da->Formatversion;
	
	getline(infile, line); lineStream.str(line);  lineStream.clear();
	lineStream >> da->BAopt_nocores >> da->BAopt_timeout_sec >> da->BAopt_maxiter;

  da->BAopt_PTLoss.resize(2);
  da->BAopt_PlaneLoss.resize(2);
  da->BAopt_DerivLoss.resize(2);
  da->BAopt_VPLoss.resize(2);
  
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  lineStream >> da->BAopt_PTLoss[0] >> da->BAopt_PTLoss[1];
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  lineStream >> da->BAopt_PlaneLoss[0] >> da->BAopt_PlaneLoss[1];
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  lineStream >> da->BAopt_DerivLoss[0] >> da->BAopt_DerivLoss[1];
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  lineStream >> da->BAopt_VPLoss[0] >> da->BAopt_VPLoss[1];

  
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
  da->Cams_smootherM.resize(nocams);
  da->Cams_VP_ID.resize(nocams);
  da->Cams_VPlines.resize(nocams);
	
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> nocamreproj;
  
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	lineStream >> noplanes;
	da->PlaneConstr_Planes.resize(noplanes);
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	for (T1 i = 0; i<noplanes; ++i) {
	  valarray<T2> tmp3dar(3);
	  lineStream >> tmp3dar[0] >> tmp3dar[1] >> tmp3dar[2];
	  da->PlaneConstr_Planes[i].resize(3);
    da->PlaneConstr_Planes[i][0] = tmp3dar[0];
    da->PlaneConstr_Planes[i][1] = tmp3dar[1];
    da->PlaneConstr_Planes[i][2] = tmp3dar[2];
	}
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->PlaneConstr_PlanesW));
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->PlaneConstr_fixPlane));

  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  lineStream >> novps;
  da->VPConstr_VP.resize(novps);
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  for (T1 i = 0; i<novps; ++i) {
    valarray<T2> tmp3dar(3);
    lineStream >> tmp3dar[0] >> tmp3dar[1] >> tmp3dar[2];
    da->VPConstr_VP[i].resize(3);
    da->VPConstr_VP[i][0] = tmp3dar[0];
    da->VPConstr_VP[i][1] = tmp3dar[1];
    da->VPConstr_VP[i][2] = tmp3dar[2];
  }
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->VPConstr_VPw));
  getline(infile, line); lineStream.str(line);  lineStream.clear(); 
  copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->VPConstr_fixVP));

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
	    vector<T2> tmp_oriens;
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(tmp_oriens)); 
	    da->Cams_view_orien[i].resize(tmp_oriens.size()/3);
	    for (T1 j=0; j < da->Cams_view_orien[i].size() ; ++j)
      {
        da->Cams_view_orien[i][j].resize(3);
        da->Cams_view_orien[i][j][0] = tmp_oriens[j*3+0];
        da->Cams_view_orien[i][j][1] = tmp_oriens[j*3+1];
        da->Cams_view_orien[i][j][2] = tmp_oriens[j*3+2];
	      //copy(tmp_oriens.begin()+j*3, tmp_oriens.begin()+(j+1)*3, back_inserter(da->Cams_view_orien[i][j])); 
      }
    
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->Cams_fixOrientation[i]));

	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    vector<T2> tmp_trans;
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(tmp_trans)); 
	    da->Cams_view_trans[i].resize(tmp_trans.size()/3);
	    for (T1 j=0; j<da->Cams_view_trans[i].size() ; ++j)
      {
        da->Cams_view_trans[i][j].resize(3);
        da->Cams_view_trans[i][j][0] = tmp_trans[j*3+0];
        da->Cams_view_trans[i][j][1] = tmp_trans[j*3+1];
        da->Cams_view_trans[i][j][2] = tmp_trans[j*3+2];        
	      //copy(tmp_oriens.begin()+j*3, tmp_oriens.begin()+(j+1)*3, back_inserter(da->Cams_view_trans[i][j])); 
      }

	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<bool>(lineStream), istream_iterator<bool>(), back_inserter(da->Cams_fixTranslation[i]));
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->Cams_viewids[i])); 
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T1>(lineStream), istream_iterator<T1>(), back_inserter(da->Cams_OnPlane[i]));
	    getline(infile, line); lineStream.str(line);  lineStream.clear();	
	    copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_Weight[i])); 
	    
      getline(infile, line); lineStream.str(line);  lineStream.clear();       
      vector<T2> vechelp;
      vechelp.resize(0);
      copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(vechelp));
      
      if (vechelp.size()==0) // no smoothing 
        da->Cams_smootherM[i].resize(0);
      else if (vechelp.size()== da->Cams_noviews[i]*da->Cams_noviews[i]) // full smoother matrix
      {
        da->Cams_smootherM[i].resize(da->Cams_noviews[i]);
        for (T1 j = 0; j<da->Cams_noviews[i]; ++j)
          da->Cams_smootherM[i][j].resize(da->Cams_noviews[i]);
        
        for (T1 j = 0; j<da->Cams_noviews[i]; ++j)
          for (T1 k = 0; k<da->Cams_noviews[i]; ++k)
            da->Cams_smootherM[i][k][j] = vechelp[j + k*da->Cams_noviews[i]];
      }
      else // same smoothing vector for all data points
      {
        da->Cams_smootherM[i].resize(1);
        da->Cams_smootherM[i][0].resize(vechelp.size());
        for (T1 j = 0; j<vechelp.size(); ++j)
          da->Cams_smootherM[i][0][j] = vechelp[j];
      }

      getline(infile, line); lineStream.str(line);  lineStream.clear();      
      vechelp.resize(0);
      copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(vechelp));
      if (vechelp.size()==0) // no VP constraints
      {
        da->Cams_VP_ID[i].resize(0);
        da->Cams_VPlines[i].resize(0);
      }
      else
      {
        da->Cams_VP_ID[i].resize(da->Cams_noviews[i]);
        da->Cams_VPlines[i].resize(da->Cams_noviews[i]);
        cnt = 0;        
        for (T1 j = 0; j<da->Cams_noviews[i]; ++j)
        {
          novps = (T1)vechelp[cnt]; cnt++;
          //lineStream >> novps;
          da->Cams_VP_ID[i][j].resize(novps);
          da->Cams_VPlines[i][j].resize(novps);
          for (T1 k = 0; k < novps; ++k)
          {
            vpid = (T1)vechelp[cnt]; cnt++;
            nolines = (T1)vechelp[cnt]; cnt++;
            
            //lineStream >> vpid >> nolines;
            da->Cams_VP_ID[i][j][k] = vpid;
            da->Cams_VPlines[i][j][k].resize(4*nolines);
            
            for (T1 ki = 0; ki < 4*nolines; ++ki)
            {
              da->Cams_VPlines[i][j][k][ki] = vechelp[cnt]; cnt++;;
            }
          }
        }
      }
          
	}

	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	da->Cams_reproj_view.resize(nocamreproj);
	for (T1 i = 0; i<nocamreproj; ++i) {
	  lineStream >> tmp2d_T1[0] >> tmp2d_T1[1];
    da->Cams_reproj_view[i].resize(2);
    da->Cams_reproj_view[i][0] = tmp2d_T1[0];
    da->Cams_reproj_view[i][1] = tmp2d_T1[1];
	}
	getline(infile, line); lineStream.str(line);  lineStream.clear();		
	da->Cams_reproj_pos.resize(nocamreproj);
	for (T1 i = 0; i<nocamreproj; ++i) {
	  lineStream >> tmp2d_T2[0] >> tmp2d_T2[1];
    da->Cams_reproj_pos[i].resize(2);
    da->Cams_reproj_pos[i][0] = tmp2d_T2[0];
    da->Cams_reproj_pos[i][1] = tmp2d_T2[1];
	}
	getline(infile, line); lineStream.str(line);  lineStream.clear();		
	copy(istream_iterator<T2>(lineStream), istream_iterator<T2>(), back_inserter(da->Cams_reproj_weight)); 
	
	da->Pt3d.resize(no3Dpoints);
	getline(infile, line); lineStream.str(line);  lineStream.clear();	
	for (T1 i = 0; i<no3Dpoints; ++i) {
	  lineStream >> tmp3d_valarr[0] >> tmp3d_valarr[1] >> tmp3d_valarr[2];
	  da->Pt3d[i].resize(3);
    da->Pt3d[i][0] = tmp3d_valarr[0];
    da->Pt3d[i][1] = tmp3d_valarr[1];
    da->Pt3d[i][2] = tmp3d_valarr[2];
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
	}
	
	
	da->Residual=-1;
    
	infile.close();
	
	CheckErrors();
};

template <class T1, class T2>
Parser<T1, T2>::Parser(BAdata<T1,T2>* da_in){
  da = da_in;
  CheckErrors();
};


template <class T1, class T2>
void Parser<T1, T2>::WriteToFile(string filename)
{
	T1 nocams = da->Cams_noviews.size();
	T1 no3Dpoints = da->Pt3d.size();
	T1 noreprojcams = da->Cams_reproj_view.size();
	T1 noplanes = da->PlaneConstr_Planes.size();	
  T1 novps = da->VPConstr_VP.size();  
	
  cout << "Save data to file: " << filename << endl;

	ofstream outfile(filename, ofstream::out);
	outfile << setprecision(20);
	
	outfile << da->Formatversion <<  endl;
	outfile << da->BAopt_nocores <<  " " << da->BAopt_timeout_sec << " " << da->BAopt_maxiter << endl;

  outfile << da->BAopt_PTLoss[0] << " " << da->BAopt_PTLoss[1] << endl;
  outfile << da->BAopt_PlaneLoss[0] << " " << da->BAopt_PlaneLoss[1] << endl;
  outfile << da->BAopt_DerivLoss[0] << " " << da->BAopt_DerivLoss[1] << endl;
  outfile << da->BAopt_VPLoss[0] << " " << da->BAopt_VPLoss[1] << endl;
  
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

  outfile << novps << endl;  
  
  for (T1 i=0; i<novps; ++i)
    outfile << da->VPConstr_VP[i][0] << " " << da->VPConstr_VP[i][1] << " " << da->VPConstr_VP[i][2] << " ";
  outfile << endl;  
  copy(da->VPConstr_VPw.begin(), da->VPConstr_VPw.end(), ostream_iterator<T2>(outfile, " "));  outfile << endl; 
  copy(da->VPConstr_fixVP.begin(), da->VPConstr_fixVP.end(), ostream_iterator<bool>(outfile, " "));  outfile << endl; 

  
  
	for (T1 i=0; i<nocams; ++i) {
    copy(da->Cams_fc[i].begin(), da->Cams_fc[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
		copy(da->Cams_cc[i].begin(), da->Cams_cc[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;				
		copy(da->Cams_kc[i].begin(), da->Cams_kc[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
		copy(da->Cams_fixInternals[i].begin(), da->Cams_fixInternals[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;		
		
		for (T1 j=0; j<da->Cams_view_orien[i].size(); ++j)
		  outfile << da->Cams_view_orien[i][j][0] << " " << da->Cams_view_orien[i][j][1] << " " << da->Cams_view_orien[i][j][2] << " ";
		outfile << endl;
		
		copy(da->Cams_fixOrientation[i].begin(), da->Cams_fixOrientation[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;
		
		for (T1 j=0; j<da->Cams_view_trans[i].size(); ++j)
		  outfile << da->Cams_view_trans[i][j][0] << " " << da->Cams_view_trans[i][j][1] << " " << da->Cams_view_trans[i][j][2] << " ";
		outfile << endl;
		
		copy(da->Cams_fixTranslation[i].begin(), da->Cams_fixTranslation[i].end(), ostream_iterator<bool>(outfile, " "));  outfile << endl;		
		copy(da->Cams_viewids[i].begin(), da->Cams_viewids[i].end(), ostream_iterator<T1>(outfile, " "));  outfile << endl;
		copy(da->Cams_OnPlane[i].begin(), da->Cams_OnPlane[i].end(), ostream_iterator<T1>(outfile, " "));  outfile << endl;		
		copy(da->Cams_Weight[i].begin(), da->Cams_Weight[i].end(), ostream_iterator<T2>(outfile, " "));  outfile << endl;
    
    for (T1 j=0; j< da->Cams_smootherM[i].size(); ++j)
      copy(da->Cams_smootherM[i][j].begin(), da->Cams_smootherM[i][j].end(), ostream_iterator<T2>(outfile, " "));
    outfile << endl;
    
    if (da->Cams_VP_ID[i].size() > 0)
    {
      for (T1 j=0; j< da->Cams_noviews[i] ; ++j)
      {
        outfile << da->Cams_VP_ID[i][j].size() << " ";
        for (T1 k=0; k< da->Cams_VP_ID[i][j].size() ; ++k)
        {
          outfile << da->Cams_VP_ID[i][j][k] << " ";
          outfile << da->Cams_VPlines[i][j][k].size()/4 << " ";
          for (T1 ki=0; ki < da->Cams_VPlines[i][j][k].size() ; ++ki)
            outfile << da->Cams_VPlines[i][j][k][ki] << " ";
        }
      }
      outfile << endl;  
    }
    else
    {
      outfile << endl;  
    }

	}
	
	for (T1 i=0; i<noreprojcams; ++i)
	    outfile << da->Cams_reproj_view[i][0] << " " << da->Cams_reproj_view[i][1] << " ";
	outfile << endl;
  
	for (T1 i=0; i<noreprojcams; ++i)
      outfile << da->Cams_reproj_pos[i][0]  << " " << da->Cams_reproj_pos[i][1]  << " ";    

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
	    nopointsonplane++;
		
  
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
	}
	
	outfile << da->Residual << endl;
	
	outfile.close(); 
};


template <class T1, class T2>
bool Parser<T1, T2>::CheckErrors()
{
  bool ret = 1; // data struct is consistent
  
  T1 nopts = da->Pt3d.size();
  T1 noplanes = da->PlaneConstr_Planes.size();
  T1 nocams = da->Cams_noviews.size();
  
  if (nopts==0 || nocams==0)
  {
    ret *= 0;
    cout << "Warning: No 3d points or cameras supplied" << endl;
  }
  
  if (da->PlaneConstr_PlanesW.size() != noplanes || da->PlaneConstr_fixPlane.size() != noplanes)
  {
    ret *= 0;
    cout << "Warning: Incorrect number of plane options/weights" << endl;
  }
  
  for (T1 i=0; i < nocams ; ++i)
  {
    if (da->Cams_noviews[i] == 0)
    {
      ret *= 0;
      cout << "Warning: Zero views in camera " << i << endl;    
    }
    if (da->Cams_fc[i].size() < 1 || da->Cams_fc[i].size() > 2)
    {
      ret *= 0;
      cout << "Warning: Incorrect number of values for focal length in camera " <<  i << endl;    
    }
    if (  !(da->Cams_cc[i].size() == 0  || da->Cams_cc[i].size() == 2)  )
    {
      ret *= 0;
      cout << "Warning: Incorrect number of values for principal point in camera " <<  i << endl;    
    }
    if (  !(da->Cams_kc[i].size() == 0  || da->Cams_kc[i].size() == 1 || da->Cams_kc[i].size() == 2 || da->Cams_kc[i].size() == 4 || da->Cams_kc[i].size() == 5)  )
    {
      ret *= 0;
      cout << "Warning: Incorrect number of values for radial distortion in camera " <<  i << endl;    
    }    
    if ( da->Cams_fixInternals[i].size() != 3)
    {
      ret *= 0;
      cout << "Warning: Incorrect number of values for fixed parameters in camera " <<  i << endl;    
    }    
    if ( da->Cams_fixInternals[i].size() != 3)
    {
      ret *= 0;
      cout << "Warning: Incorrect number of values for fixed parameters in camera " <<  i << endl;    
    }    
    if ( !(da->Cams_view_orien[i].size() == 1 ||  da->Cams_view_orien[i].size() == da->Cams_noviews[i]))
    {
      ret *= 0;
      cout << "Warning: Number of views in camera " <<  i << " and number of orientation parameters mismatch" << endl;    
    }    
    if ( !(da->Cams_view_trans[i].size() == 1 ||  da->Cams_view_trans[i].size() == da->Cams_noviews[i]))
    {
      ret *= 0;
      cout << "Warning: Number of views in camera " <<  i << " and number of translation parameters mismatch" << endl;    
    }    
    for (T1 k=0; k < da->Cams_view_orien[i].size() ; ++k)
    {
      if (da->Cams_view_orien[i][k].size() != 3)
      {
        ret *= 0;
        cout << "Warning: Orientation vector in camera " <<  i << " is too long" << endl;    
      }
    }
    for (T1 k=0; k < da->Cams_view_trans[i].size() ; ++k)
    {
      if (da->Cams_view_trans[i][k].size() != 3)
      {
        ret *= 0;
        cout << "Warning: Translation vector in camera " <<  i << " is too long" << endl;    
      }
    }
    
    if ( !(da->Cams_Weight[i].size() == 1 ||  da->Cams_Weight[i].size() == da->Cams_noviews[i]))
    {
      ret *= 0;
      cout << "Warning: Number of views in camera " <<  i << " and number of camera weights mismatch" << endl;    
    }            
  }

for (T1 k=0; k < nopts ; ++k)
{
    if (da->Pt3d[k].size() != 3)
    {
      ret *= 0;
      cout << "Warning: 3D Point vector is not of correct size" << endl;    
    }
    
    if (da->PtReprojView[k].size() != da->PtReprojPos[k].size())
    {
      ret *= 0;
      cout << "Warning: 3D Point reprojections (view, pos) have different length (" << k << " " << da->PtReprojView[k].size() << ", " << da->PtReprojPos[k].size()  << ")" << endl;    
    }
}

if ((da->Pt3d.size() !=  da->PtReprojView.size()) || (da->Pt3d.size() !=  da->PtReprojPos.size()) )
{
  ret *= 0;
  cout << "Warning: Incorrect number of reprojections for 3d point" << endl;    
}
  
if (da->PtOnPlane.size() !=  da->Pt3d.size() &&  da->PtOnPlane.size() != 0)
{
  ret *= 0;
  cout << "Warning: Incorrect number of plane ids for points" << endl;    
}
  
if (da->PtWeight.size() !=  da->Pt3d.size() &&  da->PtWeight.size() != 1)
{
  ret *= 0;
  cout << "Warning: Incorrect number of weights for cameras" << endl;    
}

if (da->PtFixPosition.size() !=  da->Pt3d.size() &&  da->PtFixPosition.size() != 1)
{
  ret *= 0;
  cout << "Warning: Incorrect number of fixed cameras" << endl;    
}
  return ret;
};




template <class T1, class T2>
void Parser<T1, T2>::NormalizeData()
{
  vector<vector<T2>> varvect(3); 
  
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
	
      varvect[0].push_back(tr[0]);
      varvect[1].push_back(tr[1]);
      varvect[2].push_back(tr[2]);
    }
  }
  
  for (T1 i=0; i< da->Pt3d.size(); ++i)
  {
      varvect[0].push_back(da->Pt3d[i][0]);
      varvect[1].push_back(da->Pt3d[i][1]);
      varvect[2].push_back(da->Pt3d[i][2]);
  }
  
  da->median3d.resize(3);
  
  // Compute median
  for (T1 i=0; i< 3; ++i)
  {
    vector<T2> tmp = varvect[i];
    sort(tmp.begin(), tmp.end());
    da->median3d[i] = tmp[int(tmp.size()/2)];
  }
  
//cout << da->median3d[0] << " " << da->median3d[1] << " " << da->median3d[2] << endl;

  // Compute Variance of L1 distances of 
//   vector<T2> tmp(varvect[0].size());
//   for (T1 i=0; i< varvect[0].size(); ++i)
//   {
//     tmp[i] = abs((varvect[0][i] - da->median3d[0])) +  abs((varvect[1][i] - da->median3d[1])) + abs((varvect[2][i] - da->median3d[2])); // L1 norm
//     tmp[i] = tmp[i]*tmp[i];
//     cout << "T " << tmp[i] << endl;
//   }

  //cout << tmp[0] << tmp[1] << tmp[2] << endl; 
  //T2 accum = inner_product( tmp.begin(), tmp.end(), tmp.begin(), 0 );
  //da->stdev3d = sqrt(accum / (tmp.size() - 1));


  T2 di = 0;
  for (T1 i=0; i< varvect[0].size(); ++i)
  {
    di +=  sqrt((varvect[0][i] - da->median3d[0]) * (varvect[0][i] - da->median3d[0]) +  
                (varvect[1][i] - da->median3d[1]) * (varvect[1][i] - da->median3d[1]) +
                (varvect[2][i] - da->median3d[2]) * (varvect[2][i] - da->median3d[2]));
  }

  da->stdev3d = sqrt(   di / (varvect[0].size()-1) );
  
  // move cameras and points to zero median and unit isotropic stdev.
  for (T1 i=0; i< da->Pt3d.size(); ++i)
  {
    da->Pt3d[i] -= da->median3d;
    da->Pt3d[i] /= da->stdev3d;
  }
  
  for (T1 i=0; i< da->Cams_noviews.size(); ++i)
  {
    for (T1 j=0; j< da->Cams_view_trans[i].size(); ++j)
    {
        da->Cams_view_trans[i][j] -= da->median3d;
        da->Cams_view_trans[i][j] /= da->stdev3d;
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

    
    planenorms -= resutil::DotProduct(da->PlaneConstr_Planes[i], da->median3d);
    planenorms /= da->stdev3d;
    
    da->PlaneConstr_Planes[i][0] *= planenorms;
    da->PlaneConstr_Planes[i][1] *= planenorms;
    da->PlaneConstr_Planes[i][2] *= planenorms;
  }
};

template <class T1, class T2>
void Parser<T1, T2>::UnnormalizeData()
{
  for (T1 i=0; i< da->Pt3d.size(); ++i)
  {
    da->Pt3d[i] *= da->stdev3d;
    da->Pt3d[i] += da->median3d;
  }
  
  for (T1 i=0; i< da->Cams_noviews.size(); ++i)
  {
    for (T1 j=0; j< da->Cams_view_trans[i].size(); ++j)
    {
        da->Cams_view_trans[i][j] *= da->stdev3d;
        da->Cams_view_trans[i][j] += da->median3d;
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

    planenorms *= da->stdev3d;    
    planenorms += resutil::DotProduct(da->PlaneConstr_Planes[i], da->median3d);
    
    da->PlaneConstr_Planes[i][0] *= planenorms;
    da->PlaneConstr_Planes[i][1] *= planenorms;
    da->PlaneConstr_Planes[i][2] *= planenorms;
  }  
  
  // Clean Up:
  // Scale VP directions back to unit norm (small deviations possible after optimization)
  for (T1 i=0; i<da->VPConstr_VP.size(); i++)
  {
    T2 norm;
    norm = sqrt(da->VPConstr_VP[i][0]*da->VPConstr_VP[i][0] + da->VPConstr_VP[i][1]*da->VPConstr_VP[i][1] + da->VPConstr_VP[i][2]*da->VPConstr_VP[i][2]);
    da->VPConstr_VP[i][0] /= norm;
    da->VPConstr_VP[i][1] /= norm;
    da->VPConstr_VP[i][2] /= norm;
  }
  
}
