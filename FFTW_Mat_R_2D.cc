//implement the 2D real Matrix class

#include <stdlib.h>
#include <fftw3.h>

#include "FFTW_Mat_R_2D.h"

//constructor
FFTW_Mat_R_2D::FFTW_Mat_R_2D(int nr, int nc)
: nrow_(nr), ncol_(nc){
	val_ = (double**)malloc(nrow_ * sizeof(double*)); //assign of pointers to rows
	val_s_ = (double *)fftw_malloc(ncol_ * nrow_ * sizeof(double)); //array of matrix values
	for (int i = 0; i < nrow_; i++) {
		val_[i] = val_s_ + i*ncol_; // set pointers to rows of data
	}	
}

//destructor
FFTW_Mat_R_2D::~FFTW_Mat_R_2D() {
	free(val_);
	fftw_free(val_s_);  	
}
