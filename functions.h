#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "FFTW_Mat_R_2D.h"

//calculate the characteristic length, and this functions is only valid when Dx = Dy
double length_R(FFTW_Mat_R_2D *psi, double dx);
//calculate the average of a real Matrix
double cal_ave(FFTW_Mat_R_2D *mat);
//calculate the max absolute value of a real Matrix
double cal_max_abs(FFTW_Mat_R_2D *mat);


#endif  // FUNCTIONS_H_




