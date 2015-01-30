#ifndef OneLayer_H_
#define OneLayer_H_


#include <fftw3.h>

#include "FFTW_Mat_C_2D.h"
#include "FFTW_Mat_R_2D.h"
#include "FFTW_Mat_R_3D.h"

//define a OneLayer class
class OneLayer {
public:

	
	//Constructor
	OneLayer(int Nx, int Ny, double Dx, double Dy,
			 double a, double b, double w, double lambda, 
			 double eta_M, double eta_S, double gamma, double M, int i, double thick_sol1);
	//Destructor
	~OneLayer(); 
	
	int init_psi_Gauss(double ave, double sigma); //initialize psi by Gaussian with average and sigma	
	int init_psi_Sin(int nwave, double delta); //initialize psi by sin shape curved interface in middle of y, 
											//between bulk phases (psi = 1 and -1)
											//nwave: wavenumber, delta: amplitude in units of grid size
	int cal_mu_boundary(FFTW_Mat_R_2D *psi_s);     //Calculate chemical potential from psi and psi_s
    int cal_mu_middle(FFTW_Mat_R_2D *psi_up, FFTW_Mat_R_2D *psi_down);
	int cal_W_hat();                      //Calculate driving force in Fourier space from psi and mu
	int cal_uM();  //from W_hat, W_hat_s, T_ii_hat and T_ij_hat
    int init_psi_circle(double, double, double);
    int init_psi_2circle(double, double, double, double, double);
	int cal_rhs();                        //Calculate the right hand side of governing equation  
				                          //in Fourier space from psi, uM and mu, 
	                                      //and store the result into rhs_hat
	int Nx(){ return Nx_; };              //return system size in x 
	int Ny(){ return Ny_; };              //return system size in y
	double Dx(){ return Dx_; };           //return step size in x
	double Dy(){ return Dy_; };           //return step size in y
	double eta_M(){ return eta_M_; };     //return membrane's viscosity
	double eta_S(){ return eta_S_; };     //return solvent's viscosity
	
	//DFT functions
	int psi_to_psi_hat(); //transform psi from real to Fourier space
	int psi_hat_to_psi(); //transform psi from Fourier to real space
	
	//Matrix
	FFTW_Mat_R_2D *psi_, *mu_, *uMx_, *uMy_; 
	FFTW_Mat_C_2D *psi_hat_, *Wx_hat_, *Wy_hat_, *rhs_hat_;
	const int index;
	const double thick_sol;
    FFTW_Mat_C_2D *mu_hat_, *uMx_hat_, *uMy_hat_, *psi_uMx_hat_, *psi_uMy_hat_;
//the following are added to test the code
		fftw_plan p_out_psi_, p_out_uMx_, p_out_uMy_;	
    //Matrix
//	FFTW_Mat_R_2D *Wx_, *Wy_, *psi_uMx_, *psi_uMy_; 
//	FFTW_Mat_R_3D *T_ii_hat_, *T_ij_hat_;
//	FFTW_Mat_C_2D *mu_hat_, *uMx_hat_, *uMy_hat_, *psi_uMx_hat_, *psi_uMy_hat_;
//the above are added to test the code	
	
	int cal_uM_hat(FFTW_Mat_C_2D *Wx_hat_s,   //Calculate uMx and uMy in Fourier space
				   FFTW_Mat_C_2D *Wy_hat_s);  //from W_hat, W_hat_s, T_ii_hat and T_ij_hat
	int cal_T_hat();            //and eta_M_s and eta_S_s from the other layer
    FFTW_Mat_R_2D *Wx_, *Wy_, *psi_uMx_, *psi_uMy_;

private:	                                               	                                           
	int cal_W();                              //Calculate driving force in real space from psi and mu
	//int cal_uM_hat(FFTW_Mat_C_2D *Wx_hat_s,   //Calculate uMx and uMy in Fourier space
	//			   FFTW_Mat_C_2D *Wy_hat_s);  //from W_hat, W_hat_s, T_ii_hat and T_ij_hat
	                                                                  
	//DFT functions
	int W_to_W_hat();           //transform driving force from real to Fourier space
	int uM_hat_to_uM();         //transform uM from Fourier to real space
	int psi_uM_to_psi_uM_hat(); //transform psi_uM from real to Fourier space
	int mu_to_mu_hat();         //transform the chemical potential from real to Fourier space
	
	//Matrix
	//FFTW_Mat_R_2D *Wx_, *Wy_, *psi_uMx_, *psi_uMy_;
	FFTW_Mat_R_3D *T_ii_hat_, *T_ij_hat_;
	//FFTW_Mat_C_2D *mu_hat_, *uMx_hat_, *uMy_hat_, *psi_uMx_hat_, *psi_uMy_hat_;
	
	//Parameters
	const int Nx_, Ny_;          // system size
	const double Dx_, Dy_;       // space step
	const double eta_M_, eta_S_; // (eta_M_, eta_S_): membrane and solvent viscosity
	const double a_, b_, w_;     // (a,b,w): coefficients for free energy density 
	const double lambda_;        // lambda: coupling strength
	const double gamma_;         // gamma_: friction coefficient
	const double M_;             // M: mobility
	//FFTW plans
	fftw_plan p_in_psi_, p_in_Wx_, p_in_Wy_, p_in_psi_uMx_, p_in_psi_uMy_, p_in_mu_;
	
};

#endif  // OneLayer_H_



