BA_Matlab.
=========

Bundle Adjustment for Matlab


*** OVERVIEW

This code package constaints Bundle Adjustment software with support for multiple residual terms

1) Std. 3D->2D point reprojection error.

The internal camera calibration can be set fixed or as part of the optimization. 
Multiple sets of views with shared calibration can be defined. 
3D points and cameras can (independently) be set as fixed or part of the optimization.

2) 3D point plane fitting error.

Every point can be assigned to [0, inf] planes. 
The 3D point-plane projection error is used as a optimization residual.
The plane normals can be set fixed or as part of the optimization.

3) Camera location plane fitting error.

Same as 2) but for camera location.
This is useful if a camera plane / ground plane is known.
This is used in: [1].

4) Camera-Camera visibility

Use mutual camera-camera visibility as 3D (camera location) to 2D 
(observed camera location) reprojection constraint. Otherwise similar to 1).
Ref: [2].

5) Vanishing point reprojection error

For given line segments in image, and associated VPs (3D direction unit vectors)
use VP reprojection in image as residual. Every VP can be set fixed or as 
part of the optimization.
Ref: [3].

6) Camera smoothing for sequences
A smoother matrix (Reinsch Form, See [4]) is used to penalize non-smooth
camera paths.
This is used in: [1].

If you use this software please cite [1].


[1] T. Kroeger, Video Registration to SfM Models, ECCV 2014 <br />
[2] T. Thormählen, et al., Exploiting Mutual Camera Visibility in Multi-Camera Motion Estimation, ISVC 2009 <br />
[3] J.-P. Tardif, Non-Iterative Approach for Fast and Accurate Vanishing Point Detection, ICCV 2009. Fig. 4.<br />
[4] T. Hastie, et al. The Elements of Statistical Learning, Data Mining, Inference, and Prediction. Springer, 2nd, 2009 edn. (2009), Page 154 f.<br />

The software package makes use of the Ceres Solver: http://ceres-solver.org/



*** INSTALLATION / SETUP

1) Download and install the Ceres Solver: http://ceres-solver.org/

Make sure to activate shared library support.
Last tested with ceres v1.8.0.rc2.

2) Adapt: src/CMakeLists.txt and src/mex_link.sh.

Change $MATLAB_INCLUDE_DIR in CMakeLists.txt to you local Matlab include directory.
Change $CERESLIBSPATH in CMakeLists.txt to you Ceres install directory if it is not already visible from matlab.
In mex_link.sh change path to your local matlab binary.

3) Compile with 'cmake' + 'make'.

Test binary 'BAdjustbin' with 'make test'.
Test mex-binary 'BAdjustMex' by opening MATLAB and running BATestSuite/run_BA_TestSuite.m

4) Currently no detailed documentation is available.

From matlab you can call the library directly with BAdjustMax(data), with 'data' being the input and output structure for a BA problem. 
See examples in BATestSuite/run_BA_TestSuite.m of how the datastructure is created and filled.
Please make sure all data types, dimensions are correct in the input struct.
Input verification is still very basic.

You can also call the libray over the shell by ./BAdjustbin inputfile outputfile.
Please see BATestSuite/func_readwrite_BA_problem.m for example code for writing and parsing of the text files.

If you need to call the library from C directly, have a look at src/ceresBA_datastruct.h for documentation of the data structure, and src/ceresBAbin.cpp of an example of how the solver is called.

If you use this software please cite [1].


