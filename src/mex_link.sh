#!/bin/sh


SOURCEPATH=$1
STATIC_LIB_NAME=$2
CERES_LIB_PATH=$3
CERES_LIBS=$4


/usr/local/MATLAB/MATLAB_Production_Server/R2013a/bin/matlab -nodesktop -nojvm -nosplash -r "mex CFLAGS='-Wall -fpic -std=c++11'  ${SOURCEPATH}/mex_stub.cpp -L./ -l$STATIC_LIB_NAME $CERES_LIB_PATH $CERES_LIBS -output $STATIC_LIB_NAME;  exit"
#/usr/pack/matlab-8.3r2014a-fg/bin/matlab -nodesktop -nojvm -nosplash -r "mex CFLAGS='-Wall -fpic -std=c++11'  ${SOURCEPATH}/mex_stub.cpp -L./ -l$STATIC_LIB_NAME $CERES_LIB_PATH $CERES_LIBS -output $STATIC_LIB_NAME;  exit"
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Change to path of your local MATLAB installation 










#mex CFLAGS='-Wall -fpic -std=c++11'  /home/kroegert/local/Code/BA_Matlab/src/mex_stub.cpp  -L./ -lBAdjustMex -L/scratch_net/zinc/kroegert/usrlocal/lib -lglog -lceres -output BAdjustMex

