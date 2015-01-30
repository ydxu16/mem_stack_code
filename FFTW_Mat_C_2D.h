#ifndef FFTW_MAT_C_2D_H_
#define FFTW_MAT_C_2D_H_

#include <fftw3.h>

//define a 2D complex Matrix class
class FFTW_Mat_C_2D {
public:
	FFTW_Mat_C_2D(int nr, int nc); //constructor
	~FFTW_Mat_C_2D();  //destructor
	fftw_complex* val_s_; // space for 2-D array of values
	fftw_complex** val_;  // space for pointers to rows
	const int nrow_, ncol_; // number of rows/columns
};

#endif  // FFTW_MAT_C_2D_H_




