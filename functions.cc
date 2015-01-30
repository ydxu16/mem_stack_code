#include "functions.h"
#include "FFTW_Mat_R_2D.h"
#include <math.h>

//calculate the characteristic length, and this functions is only valid when Dx = Dy
double length_R(FFTW_Mat_R_2D *psi, double dx){
	int i,j;
	int right, up;
	int sum = 0;
	double rc;
	
	for (i = 0; i < psi->nrow_; i++) {
		for (j = 0; j < psi->ncol_; j++) {
			right = i + 1;
			up = j + 1;			
			if (i == (psi->nrow_ - 1)) {
				right = 0;
			}			
			if (j == (psi->ncol_ - 1)) {
				up = 0;
			}
			//sum up the edge numbers
			if ( (psi->val_[i][j]*psi->val_[right][j]) < 0.0) {
				sum = sum + 1;
			}
			if ( (psi->val_[i][j]*psi->val_[i][up]) < 0.0) {
				sum = sum + 1;
			}
		}
	}
	
	rc = dx * psi->nrow_ * psi->ncol_ / (double)sum; 
	
	return rc;
	
}
//calculate the average of a real Matrix
double cal_ave(FFTW_Mat_R_2D *mat){
	int i,j;
	double ave = 0.0;
	
	for (i = 0; i < mat->nrow_; i++) {
		for (j = 0; j < mat->ncol_; j++) {
			ave = ave + mat->val_[i][j];			
		}
	}
	
	ave = ave / (mat->nrow_*mat->ncol_);	
	return ave;
}
//calculate the max absolute value of a real Matrix
double cal_max_abs(FFTW_Mat_R_2D *mat){
	int i,j;
	double max_abs = 0.0;
	
	for (i = 0; i < mat->nrow_; i++) {
		for (j = 0; j < mat->ncol_; j++) {
			if (max_abs < fabs(mat->val_[i][j])) {
				max_abs = fabs(mat->val_[i][j]);
			}
		}
	}
	
	return max_abs;
}






