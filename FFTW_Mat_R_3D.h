#ifndef FFTW_MAT_R_3D_H_
#define FFTW_MAT_R_3D_H_

//define a 3D real Matrix class
class FFTW_Mat_R_3D {
public:
	FFTW_Mat_R_3D(int n1, int n2, int n3); //constructor
	~FFTW_Mat_R_3D();  //destructor
	double* val_s_; // space for 3-D array of values
	double** val_c_; // space for pointers to columns
	double*** val_; //space for pointers to rows
	const int n1_, n2_, n3_; // number of each dimension
};

#endif  // FFTW_MAT_R_3D_H_




