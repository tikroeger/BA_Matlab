
#ifndef _ResUtil_H__
#define _ResUtil_H__


#include <valarray>
#include <iostream>

namespace resutil 
{

template<typename T>
void VecNormalizeL2_D3(T x[3])
{
  T n = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  
  x[0] /= n; 
  x[1] /= n;
  x[2] /= n;
};


template<class T2>
valarray<T2> CrossProduct(valarray<T2> x, valarray<T2> y) {
  valarray<T2> x_cross_y(3,0);
  x_cross_y[0] = x[1] * y[2] - x[2] * y[1];
  x_cross_y[1] = x[2] * y[0] - x[0] * y[2];
  x_cross_y[2] = x[0] * y[1] - x[1] * y[0];
  return (x_cross_y);
};


template<class T2>
T2 DotProduct(valarray<T2> x, valarray<T2> y) {
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
};


template<typename T>
void RadDist(const T* kc, const int radmodel, const T pt[3], T res[3]) {
    if (radmodel>0)
    {
      
      //cout << radmodel << " " << kc[0]<< endl;
        
        
        T r2 = T(pt[0]*pt[0] + pt[1] * pt[1]);
        T rc= T(1.0) + kc[0]*r2;
        T dx[2]; dx[0] = T(0.0); dx[1] = T(0.0);

        if (radmodel>1) 
        {
          rc += kc[1]*r2*r2;
          
          if(radmodel==5)
            rc += kc[4]*r2*r2*r2;
          
          if (radmodel>=3)
          {
            // tangential distortion
            dx[0] = T(2.0)*kc[2]*pt[0]*pt[1] + kc[3]*(r2 + T(2.0)*pt[0]*pt[0]);
            dx[1] = T(2.0)*kc[3]*pt[0]*pt[1] + kc[2]*(r2 + T(2.0)*pt[1]*pt[1]);
          }
            
        }
        
        res[0] = pt[0]*rc + dx[0];
        res[1] = pt[1]*rc + dx[1];
    }  
    else
    {
      res[0] = pt[0];
      res[1] = pt[1];
    }
};

template<typename T>
void RadUndist(const T* kc, const int radmodel, const T pt[3], T res[3]) {
  
  
    res[0] = pt[0];
    res[1] = pt[1];

    if (radmodel>0)
    {
      T r2, rc;
      T dx[2];
      for (int i = 0; i < 20; ++i) // iterative solution
      {
        dx[0]=T(0.0);
        dx[1]=T(0.0);

        r2 = T(res[0]*res[0] + res[1] * res[1]);
        rc = T(1.0) + kc[0]*r2;
          
        if (radmodel>1) 
        {
          rc += kc[1]*r2*r2;
          
          if (radmodel==5) 
            rc += kc[4]*r2*r2*r2;
            
          if (radmodel>=3)
          {
            // tangential distortion
            dx[0] = T(2.0)*kc[2]*res[0]*res[1] + kc[3]*(r2 + T(2.0)*res[0]*res[0]);
            dx[1] = T(2.0)*kc[3]*res[0]*res[1] + kc[2]*(r2 + T(2.0)*res[1]*res[1]);
          }
        }
        
        res[0] = (pt[0] - dx[0]) / rc;
        res[1] = (pt[1] - dx[1]) / rc;
      }
    }  

};

template <class T2>
void RotatePoint(valarray<T2> pt,valarray<T2> ori, valarray<T2> *res) {
  valarray<T2> w(T2(0),3);
  T2 sintheta;
  T2 costheta;
  
     
  const T2 theta2 = DotProduct(ori,ori);
  //const T2 theta2 = (ori * ori).sum;   // C++11, not supported in gcc 4.4
  if (theta2 > 0.0) {
    // Away from zero, use the rodriguez formula
    //
    //   result = pt costheta +
    //            (w x pt) * sintheta +
    //            w (w . pt) (1 - costheta)
    //
    // We want to be careful to only evaluate the square root if the
    // norm of the angle_axis vector is greater than zero. Otherwise
    // we get a division by zero.
    //
    const T2 theta = sqrt(theta2);
    
    w[0] = ori[0] / theta;
    w[1] = ori[1] / theta;
    w[2] = ori[2] / theta;
    costheta = cos(theta);
    sintheta = sin(theta);
    valarray<T2> w_cross_pt = CrossProduct(w, pt);
    T2 w_dot_pt = DotProduct(w,pt);
    for (int i = 0; i < 3; ++i) 
    {
      //cout << pt[i] * costheta << " " << w_cross_pt[i] * sintheta << " " << w[i] * (T2(1.0) - costheta) * w_dot_pt << " " << pt[i] * costheta + w_cross_pt[i] * sintheta + w[i] * (T2(1.0) - costheta) * w_dot_pt;      
      (*res)[i] = pt[i] * costheta +
          w_cross_pt[i] * sintheta +
          w[i] * (T2(1.0) - costheta) * w_dot_pt;
	//  cout << " " << res[i] << endl;
    }
  //cout << res[0] << res[1] << res[2] << endl;    
  } else {
    // Near zero, the first order Taylor approximation of the rotation
    // matrix R corresponding to a vector w and angle w is
    //
    //   R = I + hat(w) * sin(theta)
    //
    // But sintheta ~ theta and theta * w = angle_axis, which gives us
    //
    //  R = I + hat(w)
    //
    // and actually performing multiplication with the point pt, gives us
    // R * pt = pt + w x pt.
    //
    // Switching to the Taylor expansion at zero helps avoid all sorts
    // of numerical nastiness.
    valarray<T2> w_cross_pt = CrossProduct(ori, pt);
    for (int i = 0; i < 3; ++i) {
      (*res)[i] = pt[i] + w_cross_pt[i];
    }
  }

  //cout << (*res)[0] << (*res)[1] << (*res)[2] << " end" << endl;
  //return res;
}

template<typename T> inline
void AngleAxisToRMatrix(const T aaxis[3], T result[9]) {

    T alpha=sqrt(aaxis[0]*aaxis[0] + aaxis[1]*aaxis[1] + aaxis[2]*aaxis[2]); // angle in radians
    
    result[0] = T(1);  result[1] = T(0); result[2] = T(0);
    result[3] = T(0);  result[4] = T(1); result[5] = T(0);
    result[6] = T(0);  result[7] = T(0); result[8] = T(1);

    /*result[0] = 0;  result[1] = 0; result[2] = 0;
    result[3] = 0;  result[4] = 0; result[5] = 0;
    result[6] = 0;  result[7] = 0; result[8] = 0;*/
    
    if (abs(alpha)> T(1e-10))
    {
      T c = sin(alpha) / alpha;
      T d = (T(1)-cos(alpha)) / (alpha*alpha);

//cout << "Alpha: " << alpha << " " << " c: " << c<< " d: " << d << endl;
      
      result[0] += + T(0)       + d*(-aaxis[2]*aaxis[2] - aaxis[1]*aaxis[1] );
      result[1] += - c*aaxis[2] + d*( aaxis[1]*aaxis[0] );
      result[2] += + c*aaxis[1] + d*( aaxis[2]*aaxis[0] );
      result[3] += + c*aaxis[2] + d*( aaxis[0]*aaxis[1] );
      result[4] += + T(0)       + d*(-aaxis[2]*aaxis[2] - aaxis[0]*aaxis[0]);
      result[5] += - c*aaxis[0] + d*( aaxis[1]*aaxis[2] );
      result[6] += - c*aaxis[1] + d*( aaxis[0]*aaxis[2] );
      result[7] += + c*aaxis[0] + d*( aaxis[1]*aaxis[2] );
      result[8] += + T(0)       + d*(-aaxis[1]*aaxis[1] - aaxis[0]*aaxis[0]);
      
/*cout << result[0] << " " << result[1] << " " <<   result[2] << endl;
cout << result[3] << " " << result[4] << " " <<   result[5] << endl;
cout << result[6] << " " << result[7] << " " <<   result[8] << endl;*/      
    }
}


// Copied from ceressolver/include/ceres/rotation.h

// xy = x cross y;
template<typename T> inline
void CrossProduct(const T x[3], const T y[3], T x_cross_y[3]) {
  x_cross_y[0] = x[1] * y[2] - x[2] * y[1];
  x_cross_y[1] = x[2] * y[0] - x[0] * y[2];
  x_cross_y[2] = x[0] * y[1] - x[1] * y[0];
}

template<typename T> inline
T DotProduct(const T x[3], const T y[3]) {
  return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}

template<typename T> inline
void AngleAxisRotatePoint(const T angle_axis[3], const T pt[3], T result[3]) {
  T w[3];
  T sintheta;
  T costheta;

  const T theta2 = DotProduct(angle_axis, angle_axis);
  if (theta2 > 0.0) {
    // Away from zero, use the rodriguez formula
    //
    //   result = pt costheta +
    //            (w x pt) * sintheta +
    //            w (w . pt) (1 - costheta)
    //
    // We want to be careful to only evaluate the square root if the
    // norm of the angle_axis vector is greater than zero. Otherwise
    // we get a division by zero.
    //
    const T theta = sqrt(theta2);
    w[0] = angle_axis[0] / theta;
    w[1] = angle_axis[1] / theta;
    w[2] = angle_axis[2] / theta;
    costheta = cos(theta);
    sintheta = sin(theta);
    T w_cross_pt[3];
    CrossProduct(w, pt, w_cross_pt);
    T w_dot_pt = DotProduct(w, pt);
    for (int i = 0; i < 3; ++i) {
      result[i] = pt[i] * costheta +
          w_cross_pt[i] * sintheta +
          w[i] * (T(1.0) - costheta) * w_dot_pt;
    }
  } else {
    // Near zero, the first order Taylor approximation of the rotation
    // matrix R corresponding to a vector w and angle w is
    //
    //   R = I + hat(w) * sin(theta)
    //
    // But sintheta ~ theta and theta * w = angle_axis, which gives us
    //
    //  R = I + hat(w)
    //
    // and actually performing multiplication with the point pt, gives us
    // R * pt = pt + w x pt.
    //
    // Switching to the Taylor expansion at zero helps avoid all sorts
    // of numerical nastiness.
    T w_cross_pt[3];
    CrossProduct(angle_axis, pt, w_cross_pt);
    for (int i = 0; i < 3; ++i) {
      result[i] = pt[i] + w_cross_pt[i];
    }
  }
}

template<typename T> inline
void AngleAxisRotatePoint_Backward(const T angle_axis[3], const T pt[3], T result[3]) {
  T w[3];
  T sintheta;
  T costheta;

  const T theta2 = DotProduct(angle_axis, angle_axis);
  if (theta2 > 0.0) {
    // Away from zero, use the rodriguez formula
    //
    //   result = pt costheta +
    //            (w x pt) * sintheta +
    //            w (w . pt) (1 - costheta)
    //
    // We want to be careful to only evaluate the square root if the
    // norm of the angle_axis vector is greater than zero. Otherwise
    // we get a division by zero.
    //
    const T theta = sqrt(theta2);
    w[0] = -angle_axis[0] / theta;
    w[1] = -angle_axis[1] / theta;
    w[2] = -angle_axis[2] / theta;
    costheta = cos(theta);
    sintheta = sin(theta);
    T w_cross_pt[3];
    CrossProduct(w, pt, w_cross_pt);
    T w_dot_pt = DotProduct(w, pt);
    for (int i = 0; i < 3; ++i) {
      result[i] = pt[i] * costheta +
          w_cross_pt[i] * sintheta +
          w[i] * (T(1.0) - costheta) * w_dot_pt;
    }
  } else {
    // Near zero, the first order Taylor approximation of the rotation
    // matrix R corresponding to a vector w and angle w is
    //
    //   R = I + hat(w) * sin(theta)
    //
    // But sintheta ~ theta and theta * w = angle_axis, which gives us
    //
    //  R = I + hat(w)
    //
    // and actually performing multiplication with the point pt, gives us
    // R * pt = pt + w x pt.
    //
    // Switching to the Taylor expansion at zero helps avoid all sorts
    // of numerical nastiness.
    T w_cross_pt[3];
    CrossProduct(angle_axis, pt, w_cross_pt);
    for (int i = 0; i < 3; ++i) {
      result[i] = pt[i] - w_cross_pt[i];
    }
  }
}

template<typename T>
void ImgPt_to_NormCoordSphere(const T* fc, const T* cc, const T* kc, const int fcnoelem, const int  modelpp, const int radmodel,   const T pt[2], T res[3])
{
 
      T tmp[2];
      tmp[0] = pt[0];
      tmp[1] = pt[1];

        
      if (modelpp)
      {
        tmp[0] = (tmp[0] - cc[0]);
        tmp[1] = (tmp[1] - cc[1]);
      }

      tmp[0] = tmp[0] / fc[0];
      
      if (fcnoelem==2)
        tmp[1] = tmp[1] / fc[1];
      else
        tmp[1] = tmp[1] / fc[0];
      
      //cout << "A: " << tmp[0] << " " << tmp[1] << endl;
      
      
      RadUndist(&(kc[0]), radmodel, tmp, res);
      res[2]=T(1.0);
      
      
      //cout << "B: " << res[0] << " " << res[1] << " " << res[2] << endl;
      //cout << "C: " << radmodel << " " << kc[0] << " " << kc[1] << endl;
      
      

      VecNormalizeL2_D3(res);
}



}


#endif // def Rotutil