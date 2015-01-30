 //implement the OneLayer class

#include "OneLayer.h"
#include "FFTW_Mat_C_2D.h"
#include "FFTW_Mat_R_2D.h"
#include "FFTW_Mat_R_3D.h"

#include <fftw3.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
using namespace std;
#define PI 3.1415926

//Constructor
OneLayer::OneLayer(int Nx, int Ny, double Dx, double Dy,
				   double a, double b, double w, double lambda, 
				   double eta_M, double eta_S, 
				   double gamma, double M, int i, double thick_sol1)
: Nx_(Nx), Ny_(Ny), Dx_(Dx), Dy_(Dy), eta_M_(eta_M), eta_S_(eta_S),
a_(a), b_(b), w_(w), lambda_(lambda), gamma_(gamma), M_(M), index(i), thick_sol(thick_sol1){
	
	//create real 2D matrix
	psi_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	mu_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	uMx_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	uMy_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	Wx_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	Wy_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	psi_uMx_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	psi_uMy_ = new FFTW_Mat_R_2D(Nx_, Ny_);
	
	//create complex 2D matrix
	psi_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	mu_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	uMx_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	uMy_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	Wx_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	Wy_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	psi_uMx_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	psi_uMy_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	rhs_hat_ = new FFTW_Mat_C_2D(Nx_, Ny_/2 + 1);
	
	//create real 3D matrix
	T_ii_hat_ = new FFTW_Mat_R_3D(Nx_, Ny_/2 + 1, 4);
	T_ij_hat_ = new FFTW_Mat_R_3D(Nx_, Ny_/2 + 1, 4);
	
	//create plans for forward DFT
	p_in_psi_ = fftw_plan_dft_r2c_2d(Nx_, Ny_, psi_->val_s_, psi_hat_->val_s_, FFTW_PATIENT);
	p_in_mu_  = fftw_plan_dft_r2c_2d(Nx_, Ny_, mu_->val_s_, mu_hat_->val_s_, FFTW_PATIENT);
	p_in_Wx_  = fftw_plan_dft_r2c_2d(Nx_, Ny_, Wx_->val_s_, Wx_hat_->val_s_, FFTW_PATIENT);
	p_in_Wy_  = fftw_plan_dft_r2c_2d(Nx_, Ny_, Wy_->val_s_, Wy_hat_->val_s_, FFTW_PATIENT);
	p_in_psi_uMx_ = fftw_plan_dft_r2c_2d(Nx_, Ny_, psi_uMx_->val_s_, psi_uMx_hat_->val_s_, FFTW_PATIENT);
	p_in_psi_uMy_ = fftw_plan_dft_r2c_2d(Nx_, Ny_, psi_uMy_->val_s_, psi_uMy_hat_->val_s_, FFTW_PATIENT);
	//create plans for backward DFT
	p_out_psi_ = fftw_plan_dft_c2r_2d(Nx_, Ny_, psi_hat_->val_s_, psi_->val_s_, FFTW_PATIENT);
	p_out_uMx_ = fftw_plan_dft_c2r_2d(Nx_, Ny_, uMx_hat_->val_s_, uMx_->val_s_, FFTW_PATIENT);
	p_out_uMy_ = fftw_plan_dft_c2r_2d(Nx_, Ny_, uMy_hat_->val_s_, uMy_->val_s_, FFTW_PATIENT);	
	
	//initialize T_ii_hat and T_ij_hat
	//cal_T_hat(eta_M_s, eta_S_s);
}

//Destructor
OneLayer::~OneLayer(){
	//destroy the plan
	fftw_destroy_plan(p_in_psi_);
	fftw_destroy_plan(p_in_mu_);
	fftw_destroy_plan(p_in_Wx_);
	fftw_destroy_plan(p_in_Wy_);
	fftw_destroy_plan(p_in_psi_uMx_);
	fftw_destroy_plan(p_in_psi_uMy_);
	
	fftw_destroy_plan(p_out_psi_);
	fftw_destroy_plan(p_out_uMx_);
	fftw_destroy_plan(p_out_uMy_);
	
	//delete the matrix
	delete psi_;
	delete mu_;
	delete uMx_;
	delete uMy_;
	delete Wx_;
	delete Wy_;
	delete psi_uMx_;
	delete psi_uMy_;
	
	delete psi_hat_;
	delete mu_hat_;
	delete uMx_hat_;
	delete uMy_hat_;
	delete Wx_hat_;
	delete Wy_hat_;
	delete psi_uMx_hat_;
	delete psi_uMy_hat_;
	delete rhs_hat_;
	
	delete T_ii_hat_;
	delete T_ij_hat_;
}

//Initialize psi by Gaussian with mean and sigma
int OneLayer::init_psi_Gauss(double ave, double sigma){
	double sum = 0.0;
	double x1, x2, w, y1, y2;
	int iset = 0;
	//srand(time(NULL));  // initialize random seed
	
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
			if (iset == 0) {
				do {
					x1 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
					x2 = 2.0 * (double)rand()/(double)RAND_MAX - 1.0;
					w = x1*x1 + x2*x2;
				} while ( w >= 1.0 );
				
				w = sqrt(-2.0*log(w)/w);
				y1 = x1 * w;
				y2 = x2 * w;				
				psi_->val_[i][j] = y2;
				iset = 1;
			}
			else {
				psi_->val_[i][j] = y1;
				iset = 0;
			}	
			sum = sum + psi_->val_[i][j];
		}
	}	
	
	double average = sum / (Nx_ * Ny_);
	
	//rescale
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
			psi_->val_[i][j] = ave + (psi_->val_[i][j] - average) * sigma; 
		}
	}
	return 0;
}

//initialize psi by sin shape curved interface in middle of y, 
//between bulk phases (psi = 1 and -1)
//nwave: wavenumber, delta: amplitude in units of grid size
int OneLayer::init_psi_Sin(int nwave, double delta){
	double cut; 
	int i, j;
	
	for (i = 0; i < psi_->nrow_; i++) {		
		cut = (Ny_-1)/2.0 + delta*sin(2.0*PI*nwave*(double)i/Nx_);
		for (j = 0; j < psi_->ncol_; j++) {			
			if (j <= cut) {
				psi_->val_[i][j] = -1.0;
			}
			else {
				psi_->val_[i][j] = 1.0;
			}			
		}
	}
	
	return 0;
}

int OneLayer::init_psi_circle(double x, double y, double r){
    int i, j;
    for (i = 0; i < psi_->nrow_; i++){
        for (j = 0; j < psi_->ncol_; j++){
            if ( ((i-x)*(i-x) + (j-y)*(j-y)) <= r*r )
                psi_->val_[i][j] = 1;
            else psi_->val_[i][j] = -1;
        }
    }
}

int OneLayer::init_psi_2circle(double x1, double y1, double x2, double y2, double r){
    int i, j;
    for (i = 0; i < psi_->nrow_; i++){
        for (j = 0; j < psi_->ncol_; j++){
            if ( (((i-x1)*(i-x1) + (j-y1)*(j-y1)) <= r*r) || (((i-x2)*(i-x2) + (j-y2)*(j-y2)) <= r*r))
                psi_->val_[i][j] = 1;
            else psi_->val_[i][j] = -1;
        }
    }
}


//Calculate chemical potential from psi and psi_s
int OneLayer::cal_mu_boundary(FFTW_Mat_R_2D *psi_s){
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
			int left = i - 1;
			int right = i + 1;
			int down = j - 1;
			int up = j + 1;
			
			if (i == 0) { left = Nx_ - 1; }
			if (i == (Nx_ - 1)) { right = 0; }
			if (j == 0) { down = Ny_ - 1; }
			if (j == (Ny_ - 1)) { up = 0; }
			
			mu_->val_[i][j] = -0.5*w_*w_*( (psi_->val_[right][j] + psi_->val_[left][j] - 2.0*psi_->val_[i][j])/(Dx_*Dx_) 
										  + (psi_->val_[i][up] + psi_->val_[i][down] - 2.0*psi_->val_[i][j])/(Dy_*Dy_) ) 
			- a_*psi_->val_[i][j] + b_*pow(psi_->val_[i][j],3) + 2.0*lambda_*(psi_->val_[i][j] - psi_s->val_[i][j]);		
			//debug:output
            
	//		cout<<"The difference is "<<(psi_->val_[i][j] - psi_s->val_[i][j])<<endl;
		}
	}
	return 0;
}

int OneLayer::cal_mu_middle(FFTW_Mat_R_2D *psi_up, FFTW_Mat_R_2D *psi_down){
    for (int i = 0; i < Nx_; i++) {
        for (int j = 0; j < Ny_; j++) {
            int left = i - 1;
            int right = i + 1;
            int down = j - 1;
            int up = j + 1;
            
            if (i == 0) { left = Nx_ - 1; }
            if (i == (Nx_ - 1)) { right = 0; }
            if (j == 0) { down = Ny_ - 1; }
            if (j == (Ny_ - 1)) { up = 0; }
            
            mu_->val_[i][j] = -0.5*w_*w_*( (psi_->val_[right][j] + psi_->val_[left][j] - 2.0*psi_->val_[i][j])/(Dx_*Dx_)
                                          + (psi_->val_[i][up] + psi_->val_[i][down] - 2.0*psi_->val_[i][j])/(Dy_*Dy_) )
            - a_*psi_->val_[i][j] + b_*pow(psi_->val_[i][j],3) + 2.0*lambda_*(psi_->val_[i][j] - psi_up->val_[i][j]) + 2.0 * lambda_*(psi_->val_[i][j] - psi_down->val_[i][j] );
        }
    }
    return 0;
}




//Calculate driving force in Fourier space from psi and mu
int OneLayer::cal_W_hat(){
	cal_W(); //calcualte drving force in real space at first
	W_to_W_hat(); //then do the DFT
// debug
/*
char file_name[50] = "Record W_hat";
 FILE* fp_Wx = fopen(file_name, "w");	
 fprintf(fp_Wx, "Wx_hat is"); 
	for(int i = 0; i < Nx_; i++){ 
		for(int j = 0; j < Ny_/2 + 1; j++){ 
		 fprintf(fp_Wx, "index_i is %d index_j is %d Wx_hat_real is %lf Wx_hat_ima is %lf Wy_hat_real is %lf Wy_hat_ima is %lf", i,j ,Wx_hat_->val_[i][j][0],  Wx_hat_->val_[i][j][1],  Wy_hat_->val_[i][j][0],  Wy_hat_->val_[i][j][1]);
		}
	}
fclose(fp_Wx);
*/
	return 0;
}

//Calculate uMx and uMy in real space from W_hat, W_hat_s, T_ii_hat and T_ij_hat
int OneLayer::cal_uM(){
	//cal_uM_hat(Wx_hat_s, Wy_hat_s); //calculate uM in Fourier space at first
	//uM_hat_to_uM(); //then do the inverse DFT
	return 0;
}

//Calculate the right hand side of governing equation in Fourier space from psi, uM and mu,
//and store the result into rhs_hat
int OneLayer::cal_rhs(){
	//calculate psi_uMx and psi_uMy in real space at first
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
            //debug
            //cout<< "psi_uMx_hat_->val_[i][j][0]" << psi_uMx_->val_[i][j] <<endl;
   ///         cout<< "psi_uMx_hat_->val_[i][j][0]" << psi_uMx_->val_[i][j]<<endl;
		
            // Restore these later
            psi_uMx_->val_[i][j] = psi_->val_[i][j] * uMx_->val_[i][j];
            psi_uMy_->val_[i][j] = psi_->val_[i][j] * uMy_->val_[i][j];
            //psi_uMx_->val_[i][j] = 0;
            //psi_uMy_->val_[i][j] = 0;
            
            
            //debug
           // printf("psi is %f \n", psi_->val_[i][j]);
           // printf("uMx is %f \n", uMx_->val_[i][j]);
            //debug
            //printf("%f",psi_uMy_->val_[i][j]); printf("\n");
		}
	}
    //debug
    //printf("check psi_uMx_ value %f", psi_uMx_->val_[5][5]);
    
	//then calculate the DFT of mu and psi_uM
	mu_to_mu_hat();
    
    //debug
    //printf("psi_uM_hat is %f  \n",psi_uMx_hat_->val_[5][5][0]);
    /*for (int i = 0; i < 64; i++){
        for (int j = 0; j < 64; j++){
            printf("i is %d, j is %d, psi_uMx_is  %f",i, j, psi_uMx_->val_[i][j]);
        }
        printf("\n");
    }*/
  //  printf("end of this step");
  //   printf("psi_uM_hat is %f  \n",psi_uMx_hat_->val_[5][5][0]);
	psi_uM_to_psi_uM_hat();
   //debug
   // printf("psi_uM_hat is %f  \n",psi_uMx_hat_->val_[5][5][0]);
    /*for (int i = 0; i < 64; i++){
        for (int j = 0; j < 64; j++){
            printf("i is %d, j is %d, psi_uMx_hat_ is  %f",i, j, psi_uMx_hat_->val_[i][j]);
        }
        printf("\n");
    }*/
    //finally calculate the rhs
	for (int i = 0; i < Nx_; i++) {
		int ii = i;
		if (ii > Nx_/2) {
			ii = ii - Nx_;
		}		
		double q1 = 2.0*PI*ii / (Nx_*Dx_);
		
		for (int j = 0; j < (Ny_/2 + 1); j++) {
			double q2 = 2.0*PI*j / (Ny_*Dy_);
			double q = sqrt(q1*q1 + q2*q2);
            //debug
           // printf("psi_uMx_hat_->val_[i][j][0] is %f \n",psi_uMx_hat_->val_[i][j][0]);
          //  printf("mu_hat_->val_[i][j][1] %f \n", mu_hat_->val_[i][j][1]);
            
            // restore later
			double x = q1 * psi_uMx_hat_->val_[i][j][0] + q2 * psi_uMy_hat_->val_[i][j][0];
			double y = q1 * psi_uMx_hat_->val_[i][j][1] + q2 * psi_uMy_hat_->val_[i][j][1];
    
     //double x = 0;
            //double y = 0;
       //debug
         //   printf("mu is %f \n", x);
            
			rhs_hat_->val_[i][j][0] = y - M_*q*q*mu_hat_->val_[i][j][0];
			rhs_hat_->val_[i][j][1] = -x - M_*q*q*mu_hat_->val_[i][j][1];
//printf("rhs_hat is %lf \n", rhs_hat_->val_[i][j][0]);
            
		}
        //debug
     //   printf("rhs hat is %f \n",rhs_hat_->val_[1][1][0]);
	}
	
	return 0;
}

//Calculate T_ii_hat_ and T_ij_hat_ from its own properties 
//and eta_M_s and eta_S_s from the other layer
int OneLayer::cal_T_hat(){
	for (int i = 0; i < Nx_; i++) {
		int ii = i;
		if (ii > Nx_/2) {
			ii = ii - Nx_;
		}		
		double qi = 2.0*PI*ii / (Nx_*Dx_);
		for (int j = 0; j < (Ny_/2 + 1); j++) {
            if (i!= 0 || j != 0){
                double qj = 2.0*PI*j / (Ny_*Dy_);
                double q = sqrt(qi*qi + qj*qj);
          //      double a1 = eta_M_*q*q - eta_S_;
           //     double a3 = eta_S_;
                //double factorii = -a1/(a3*a3-a1*a1);
                //double factorij = a3/(a3*a3 - a1*a1);
        //debug:old situation     
        double a1 = eta_M_*q*q + eta_S_*q; double a2 = eta_S_*q/tanh(2*q); double a3 = eta_S_*q / sinh(2*q); 
        if (q == 0) { a2 = 0; a3 = 0;} 
		double factorii = (a2 + a1) / ((a1 + a2)*(a1 + a2) - a3*a3); 
		double factorij = a3/ ((a1 + a2)*(a1 + a2) - a3*a3);	
		T_ii_hat_->val_[i][j][0] = factorii * (qj * qj)/ (q * q);
                T_ii_hat_->val_[i][j][1] = factorii * (-qi * qj)/ (q * q);
                T_ii_hat_->val_[i][j][2] = T_ii_hat_->val_[i][j][1];
                T_ii_hat_->val_[i][j][3] = factorii * (qi * qi)/ (q * q);
			
                T_ij_hat_->val_[i][j][0] = factorij * (qj * qj)/ (q * q);
                T_ij_hat_->val_[i][j][1] = factorij * (-qi * qj)/ (q * q);
                T_ij_hat_->val_[i][j][2] = T_ij_hat_->val_[i][j][1];
                T_ij_hat_->val_[i][j][3] = factorij * (qi * qi)/ (q * q);
            }
        }
	}
	
	//the following is some modification when i=j=0 
	T_ii_hat_->val_[0][0][0] = 0.0;
	T_ii_hat_->val_[0][0][1] = 0.0;
	T_ii_hat_->val_[0][0][2] = 0.0;
	T_ii_hat_->val_[0][0][3] = 0.0;
	
	T_ij_hat_->val_[0][0][0] = 0.0;
	T_ij_hat_->val_[0][0][1] = 0.0;
	T_ij_hat_->val_[0][0][2] = 0.0;
	T_ij_hat_->val_[0][0][3] = 0.0;
//debug
	
 char record_filename[50] = "Value_Oseen_Tensor";
    FILE* fp_test = fopen(record_filename, "w");
    
for (int i = 0; i < Nx_; i++){ 
		for (int j = 0 ; j < (Ny_ / 2 +1); j ++){ 
fprintf(fp_test, "%lf \n", T_ii_hat_->val_[i][j][0]);  	
}
}   	
        fclose(fp_test);
	return 0;
}

//Calculate driving force in real space from psi and mu
int OneLayer::cal_W(){	
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
			int left = i - 1;
			int right = i + 1;
			int down = j - 1;
			int up = j + 1;
			if (i == 0) {
				left = Nx_ - 1;
			}
			if (i == (Nx_ - 1)) {
				right = 0;
			}
			if (j == 0) {
				down =  Ny_ - 1;
			}
			if (j == (Ny_ - 1)) {
				up = 0;
			}
			
			Wx_->val_[i][j] = mu_->val_[i][j] * (psi_->val_[right][j]-psi_->val_[left][j]) / (2.0*Dx_);
			Wy_->val_[i][j] = mu_->val_[i][j] * (psi_->val_[i][up]-psi_->val_[i][down]) / (2.0*Dy_);
		}
	}
	return 0;
}

//Calculate uMx and uMy in Fourier space from W_hat, W_hat_s, T_ii_hat and T_ij_hat
int OneLayer::cal_uM_hat(FFTW_Mat_C_2D *Wx_hat_s, FFTW_Mat_C_2D *Wy_hat_s){	

		for (int i = 0; i < Nx_; i++) {		
			for (int j = 0; j < (Ny_/2 + 1); j++) {		
				int ii = i;
				if (ii > Nx_/2) {
					ii = ii - Nx_;
				}	
				double qi = 2.0*PI*ii / (Nx_*Dx_);
				double qj = 2.0*PI*j / (Ny_*Dy_);
				double q = sqrt(qi*qi + qj*qj);
				
					uMx_hat_->val_[i][j][0] = T_ii_hat_->val_[i][j][0]*Wx_hat_->val_[i][j][0]
											+T_ii_hat_->val_[i][j][1]*Wy_hat_->val_[i][j][0]
											+T_ij_hat_->val_[i][j][0]*Wx_hat_s->val_[i][j][0]
											+T_ij_hat_->val_[i][j][1]*Wy_hat_s->val_[i][j][0];
					uMx_hat_->val_[i][j][1] = T_ii_hat_->val_[i][j][0]*Wx_hat_->val_[i][j][1]
											+T_ii_hat_->val_[i][j][1]*Wy_hat_->val_[i][j][1]
											+T_ij_hat_->val_[i][j][0]*Wx_hat_s->val_[i][j][1]
											+T_ij_hat_->val_[i][j][1]*Wy_hat_s->val_[i][j][1];

					uMy_hat_->val_[i][j][0] = T_ii_hat_->val_[i][j][2]*Wx_hat_->val_[i][j][0]
											+T_ii_hat_->val_[i][j][3]*Wy_hat_->val_[i][j][0]
											+T_ij_hat_->val_[i][j][2]*Wx_hat_s->val_[i][j][0]
											+T_ij_hat_->val_[i][j][3]*Wy_hat_s->val_[i][j][0];
										
					uMy_hat_->val_[i][j][1] = T_ii_hat_->val_[i][j][2]*Wx_hat_->val_[i][j][1]
											+T_ii_hat_->val_[i][j][3]*Wy_hat_->val_[i][j][1]
											+T_ij_hat_->val_[i][j][2]*Wx_hat_s->val_[i][j][1]
											+T_ij_hat_->val_[i][j][3]*Wy_hat_s->val_[i][j][1];				
				}
			}
		
		return 0;
	}


//Define DFT functions:
//Transform psi from real to Fourier space:
int OneLayer::psi_to_psi_hat(){
	fftw_execute(p_in_psi_);
	return 0;
}
//Transform driving force from real to Fourier space
int OneLayer::W_to_W_hat(){
	fftw_execute(p_in_Wx_);
	fftw_execute(p_in_Wy_);
	return 0;
}
//Transform psi_uM from real to Fourier space
int OneLayer::psi_uM_to_psi_uM_hat(){
	fftw_execute(p_in_psi_uMx_);
	fftw_execute(p_in_psi_uMy_);
	return 0;
}
//Transform the chemical potential from real to Fourier space
int OneLayer::mu_to_mu_hat(){
	fftw_execute(p_in_mu_);
	return 0;
}
//Transform psi from Fourier to real space
int OneLayer::psi_hat_to_psi(){
	fftw_execute(p_out_psi_);
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
			psi_->val_[i][j] = psi_->val_[i][j] / (Nx_ * Ny_);
		}
	}
	return 0;
}
//Transform uM from Fourier to real space
int OneLayer::uM_hat_to_uM(){
	fftw_execute(p_out_uMx_);
	fftw_execute(p_out_uMy_);
	for (int i = 0; i < Nx_; i++) {
		for (int j = 0; j < Ny_; j++) {
			uMx_->val_[i][j] = uMx_->val_[i][j] / (Nx_ * Ny_);
			uMy_->val_[i][j] = uMy_->val_[i][j] / (Nx_ * Ny_);			
		}
	}
	return 0;
}



