//Implement a 3D real Matrix class

#include <stdlib.h>
#include <fftw3.h>

#include "FFTW_Mat_R_3D.h"

//constructor
FFTW_Mat_R_3D::FFTW_Mat_R_3D(int n1, int n2, int n3)
: n1_(n1), n2_(n2), n3_(n3){
	val_s_ = (double*)fftw_malloc(n1_*n2_*n3_*sizeof(double)); //array of matrix values
	val_c_ = (double**)malloc(n1_*n2_*sizeof(double*)); //array of pointers to columns
	val_ = (double***)malloc(n1_*sizeof(double**)); //array of pointers to rows
	for (int i = 0; i < n1_; i++) {
		for (int j = 0; j < n2_; j++) {
			int index = i * n2_ + j;
			val_c_[index] = val_s_ + index*n3_;
		}
		val_[i] = val_c_ + i*n2_;
	}
	
}

//destructor
FFTW_Mat_R_3D::~FFTW_Mat_R_3D(){
	fftw_free(val_s_);
	free(val_c_);
	free(val_);
}



