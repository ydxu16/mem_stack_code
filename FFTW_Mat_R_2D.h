#ifndef FFTW_MAT_R_2D_H_
#define FFTW_MAT_R_2D_H_

//define a 2D real Matrix class
class FFTW_Mat_R_2D {
public:
	FFTW_Mat_R_2D(int nr, int nc); //constructor
	~FFTW_Mat_R_2D();  //destructor
	double* val_s_; // space for 2-D array of values
	double** val_;  // space for pointers to rows
	const int nrow_, ncol_; // number of rows/columns
};

#endif  // FFTW_MAT_R_2D_H_




