//Implement the 2D complex Matrix class

#include <stdlib.h>
#include <fftw3.h>

#include "FFTW_Mat_C_2D.h"

//constructor
FFTW_Mat_C_2D::FFTW_Mat_C_2D(int nr, int nc)
:nrow_(nr), ncol_(nc){
	val_ = (fftw_complex**)malloc(nrow_ * sizeof(fftw_complex*)); //assign of pointers to rows
	val_s_ = (fftw_complex*)fftw_malloc(ncol_ * nrow_ *sizeof(fftw_complex)); //array of matrix values
	for (int i = 0; i < nrow_; i++)// set pointers to rows of data
		val_[i] = val_s_ + i*ncol_;
}

//destructor
FFTW_Mat_C_2D::~FFTW_Mat_C_2D(){
	free(val_);
	fftw_free(val_s_);  
}



