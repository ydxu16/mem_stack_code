#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#include <iostream>
#include <ctime>
#include "OneLayer.h"
#include "FFTW_Mat_C_2D.h"
#include "FFTW_Mat_R_2D.h"
#include "FFTW_Mat_R_3D.h"
#include "functions.h"
#include <vector>
#include <unistd.h>
#include <omp.h>

using namespace std;

#define PI 3.1415926


int Complex_Mutiply(const fftw_complex x, const fftw_complex y, fftw_complex result) {
	double a = x[0]*y[0] - x[1]*y[1];
	double b = x[1]*y[0] + x[0]*y[1];
	
	result[0] = a;
	result[1] = b;
	
	return 0;
}

class N_membrane{
    
public:
    N_membrane(double dt_, int Nstep_, int num_layers_,int Nx_, int Ny_, double Dx_,double Dy_,double a_, double b_,double w_,double M_,double eta_M_,double eta_S_, double lamda_,double gamma_,double thick_sol_, double ave_psi);
    ~N_membrane();

    const double dt; //time step size
    const int Nstep; //total time steps
    const int num_layers; // number of layers
    int Nprint; //every Nprint steps, print the result
    int Nanalysis; //every Nanalysis steps, do analysis and print out
    int Nx;
    int Ny;
    double Dx;
    double Dy;
    double a;
    double b;
    double w;
    double M;
    double eta_M;
    double eta_S;
    double lambda;
    double gamma;
    double thick_sol;
    
    double** c1;
    double** c2;
    double** c3;
 
    vector<OneLayer*> vector_layer;
    int cal_uM_x_hat_real();
    int cal_uM_x_hat_ima();
    int cal_uM_y_hat_real();
    int cal_uM_y_hat_ima();
    double* reverse_tridiag(int, double*, double, double, double); 
    int cal_mu();
    int cal_force();
    int cal_rhs();

    };

N_membrane::N_membrane(double dt_, int Nstep_, int num_layers_,int Nx_, int Ny_, double Dx_,double Dy_,double a_, double b_,double w_,double M_,double eta_M_,double eta_S_, double lamda_,double gamma_,double thick_sol_, double ave_psi):dt(dt_), Nstep(Nstep_), num_layers(num_layers_)
    {
        
        // dt = dt_; //time step size Nstep = Nstep_; //total time steps
        //num_layers = num_layers_; // number of layers
        Nprint = Nstep_/100; //every Nprint steps, print the result
        Nanalysis = Nstep_/20; //every Nanalysis steps, do analysis and print out
        
        Nx = Nx_; Ny = Ny_; Dx = Dx_;  Dy = Dy_; a = a_; b = b_; w = w_;
        
        M = M_; eta_M = eta_M_; eta_S = eta_S_;
        
        lambda = lamda_; gamma = gamma_; thick_sol = thick_sol_;
        c1 = new double* [Nx];
        c2 = new double* [Nx];
        c3 = new double* [Nx];
        for (int i = 0; i < Nx; i++){
            c1[i] = new double[Ny];
            c2[i] = new double[Ny];
            c3[i] = new double[Ny];
        }
        
        for (int i = 0; i < Nx; i++){
            int ii = i;
            if (ii > Nx/2) {
                ii = ii - Nx;
            }
            double q1 = 2.0*PI*ii / (Nx*Dx);
            for (int j = 0; j < Ny; j++){
                double q2 = 2.0*PI*j / (Ny*Dy);
                double q = sqrt(q1*q1 + q2*q2);
                c1[i][j] = eta_M*q*q + 2* eta_S/ thick_sol;
                c2[i][j] = eta_M*q*q + eta_S / thick_sol;
                c3[i][j] = -eta_S / thick_sol;

            }
        }

        
        //create the object
        for (int i = 0; i < num_layers; i++){
            OneLayer* layer_i = new OneLayer(Nx, Ny, Dx, Dy, a, b, w, lambda, eta_M, eta_S, gamma, M, i, thick_sol);
            //if (i == 0 || i ==2){
             //   layer_i->init_psi_circle(32, 32, 10);
                srand(time(NULL));
                //layer_i->init_psi_Sin(3, 10);
		layer_i->init_psi_Gauss(ave_psi, 0.05);
                sleep(2);
                vector_layer.push_back( layer_i );
            printf("average value is %lf \n", cal_ave(layer_i->psi_));
            //}
            //else{
             //   layer_i->init_psi_2circle(32, 43, 32, 21, 5);
             //   for (int p = 0; p < Nx; p++){
              //      for (int q = 0; q < Ny; q++){
              //          printf(" %d", (int)layer_i->psi_->val_[p][q]);
               //     }
               //     printf("\n");
               // }
                //vector_layer.push_back( layer_i );
           // }
            //delete layer_i;
        }
    
    }


int N_membrane::cal_mu(){
    for(vector<OneLayer*>:: size_type i = 0; i < vector_layer.size(); i++){
        if (vector_layer[i]->index == 0)
            {vector_layer[i]->cal_mu_boundary(vector_layer[i+1]-> psi_);}
        else if ((1 <= vector_layer[i]->index) &&(vector_layer[i]->index < num_layers-1) )
            {vector_layer[i]->cal_mu_middle(vector_layer[i-1]-> psi_, vector_layer[i+1]-> psi_);}
        else if (vector_layer[i]->index == num_layers-1)
            {vector_layer[i]->cal_mu_boundary(vector_layer[i-1]-> psi_);}
                    //calculate chemical potential
    }

  return 0;
}

int N_membrane::cal_force(){
    
    for(vector<OneLayer*>::size_type i = 0; i < vector_layer.size(); i++){
        /*//debug:output
        char filename_W[50] = "Wy_0";
        FILE* fp_W = fopen(filename_W, "w");
        for (int i = 0; i < Nx; i++){
            for (int j = 0; j < Ny; j++){
                fprintf(fp_W, "%lf \n", vector_layer[0]->Wy_->val_[i][j]);
            }
        }
        fclose(fp_W);
        */
        vector_layer[i]->cal_W_hat();
        //calculate driving force in Fourier space
       
    }
    return 0;
    
   
}

int N_membrane::cal_rhs(){
    for (vector<OneLayer*>::size_type i = 0; i < vector_layer.size(); i++){
        vector_layer[i]->cal_rhs();
        //calculate rhs in Fourier space
    }
    return 0;
}

double* N_membrane::reverse_tridiag(int dim, double* right_eq, double c1, double c2, double c3){ 
	
            // temporary array to reverse the matrix
            double* c_prime = new double[dim];
            double* d_prime = new double[dim];

            
            for (int i = 0; i < dim - 1; i++){
                if (i == 0)
                {c_prime[i] = c3 / c2;}
  
                else
                {c_prime[i] = c3 / (c1 - c_prime[i-1]* c3);}
                
            }
            
        
            for (int i = 0; i < dim; i++){
                if (i == 0)
                {d_prime[i] = right_eq[i] / c2;}
                else if(i < dim - 1)
                {d_prime[i] = (right_eq[i] - c3 * d_prime[i-1]) / (c1 - c_prime[i-1]* c3);}
                else if(i == dim - 1)
                {d_prime[i] = (right_eq[i] - c3 * d_prime[i-1]) / (c2 - c_prime[i-1]* c3);}
                
            }
	   double* ans = new double[dim];

            for (int i = dim - 1; i >= 0; i--){
                if (i == dim - 1)
                {ans[i] = d_prime[i];}
                else
                {ans[i] = d_prime[i] - c_prime[i] * ans[i+1];}
            }
            
            delete c_prime;
            delete d_prime;
		return ans;
}

int N_membrane::cal_uM_x_hat_real(){
    
    for (int index_i = 0; index_i < Nx; index_i++){
        
        int index_ii = index_i;
        if (index_ii > Nx/2) {
            index_ii = index_ii - Nx;
        }
        
        double qi = 2.0*PI*index_ii / (Nx*Dx);
        
        for (int index_j = 0; index_j <(Ny/2 + 1); index_j++){
            double qj = 2.0*PI*index_j / (Ny*Dy);
            double q = sqrt(qi * qi + qj* qj);
            
            
            
            //for each i, j calculate layer_i_u_hat_mx[i][j]
            double* right_eq = new double[num_layers];
            
            //?????? needed to be further considered
            if( (index_ii == 0) && (index_j == 0) ){   //discuss the case that q_i = q_j = 0;
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = 0;
                }
            }
            //??????????????????????????
            
            else{
                vector<OneLayer*>::size_type j = 0;
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = (1 - qi*qi/(q*q)) *vector_layer[j+i]->Wx_hat_->val_[index_i][index_j][0] + (- qi*qj/(q*q)) * vector_layer[j+i]->Wy_hat_->val_[index_i][index_j][0];
                }
            }
            
 	   double* ans = reverse_tridiag(num_layers, right_eq, c1[index_i][index_j], c2[index_i][index_j], c3[index_i][index_j]); 
            
            vector<OneLayer*>::size_type j = 0;
            for (int i = num_layers - 1; i >= 0; i--){
                {vector_layer[j+i]->uMx_hat_->val_[index_i][index_j][0] = ans[i];}
                vector_layer[j+i]->uMx_hat_->val_[0][0][0] = 0;
            }
            delete ans;  
            delete right_eq;
        }
    }
    
    
    return 0;
}









int N_membrane::cal_uM_x_hat_ima(){
    
    for (int index_i = 0; index_i < Nx; index_i++){
        
        int index_ii = index_i;
        if (index_ii > Nx/2) {
            index_ii = index_ii - Nx;
        }
        
        double qi = 2.0*PI*index_ii / (Nx*Dx);
        
        for (int index_j = 0; index_j < (Ny/2 + 1); index_j++){
            double qj = 2.0*PI*index_j / (Ny*Dy);
            double q = sqrt(qi * qi + qj* qj);
            
            
            
            //for each i, j calculate layer_i_u_hat_mx[i][j]
            double* right_eq = new double[num_layers];
            
            //?????? needed to be further considered
            if( (index_ii == 0) && (index_j == 0) ){   //discuss the case that q_i = q_j = 0;
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = 0;
                }
            }
            //??????????????????????????
           
            else{
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = (1 - qi*qi/(q*q)) *vector_layer[i]->Wx_hat_->val_[index_i][index_j][1] + (- qi*qj/(q*q)) * vector_layer[i]->Wy_hat_->val_[index_i][index_j][1];
                }
            }
            
            // temporary array to reverse the matrix
            double* c_prime = new double[num_layers];
            double* d_prime = new double[num_layers];
            
            for (int i = 0; i < num_layers - 1; i++){
                if (i == 0)
                {c_prime[i] = c3[index_i][index_j] / c2[index_i][index_j];}
                else
                {c_prime[i] = c3[index_i][index_j] / (c1[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                
            }
            
            
            for (int i = 0; i < num_layers; i++){
                if (i == 0)
                {d_prime[i] = right_eq[i] / c2[index_i][index_j];}
                else if(i < num_layers - 1)
                {d_prime[i] = (right_eq[i] - c3[index_i][index_j] * d_prime[i-1]) / (c1[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                else if(i == num_layers - 1)
                {d_prime[i] = (right_eq[i] - c3[index_i][index_j] * d_prime[i-1]) / (c2[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                
            }
            
            for (int i = num_layers - 1; i >= 0; i--){
                if (i == num_layers - 1)
                {vector_layer[i]->uMx_hat_->val_[index_i][index_j][1] = d_prime[i];}
                else
                {vector_layer[i]->uMx_hat_->val_[index_i][index_j][1] = d_prime[i] - c_prime[i] * vector_layer[i+1]->uMx_hat_->val_[index_i][index_j][1];}
                vector_layer[i]->uMx_hat_->val_[0][0][1] = 0;
            }
            
            delete c_prime;
            delete d_prime;
            delete right_eq;
        }
             
    }
           
    return 0;
}







int N_membrane::cal_uM_y_hat_real(){
    
    for (int index_i = 0; index_i < Nx; index_i++){
        
        int index_ii = index_i;
        if (index_ii > Nx/2) {
            index_ii = index_ii - Nx;
        }
        
        double qi = 2.0*PI*index_ii / (Nx*Dx);
        
        for (int index_j = 0; index_j < (Ny/2 + 1); index_j++){
            double qj = 2.0*PI*index_j / (Ny*Dy);
            double q = sqrt(qi * qi + qj* qj);
            
            
            
            //for each i, j calculate layer_i_u_hat_mx[i][j]
            double* right_eq = new double[num_layers];
            
            //?????? needed to be further considered
            if( (index_ii == 0) && (index_j == 0) ){   //discuss the case that q_i = q_j = 0;
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = 0;
                }
            }
            //??????????????????????????
            
            else{
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = (1 - qj*qj/(q*q)) *vector_layer[i]->Wy_hat_->val_[index_i][index_j][0] + (- qi*qj/(q*q)) * vector_layer[i]->Wx_hat_->val_[index_i][index_j][0];
                }
            }
            
            // temporary array to reverse the matrix
            double* c_prime = new double[num_layers];
            double* d_prime = new double[num_layers];
            
            for (int i = 0; i < num_layers - 1; i++){
                if (i == 0)
                {c_prime[i] = c3[index_i][index_j] / c2[index_i][index_j];}
                else
                {c_prime[i] = c3[index_i][index_j] / (c1[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                
            }
            
            
            for (int i = 0; i < num_layers; i++){
                if (i == 0)
                {d_prime[i] = right_eq[i] / c2[index_i][index_j];}
                else if(i < num_layers - 1)
                {d_prime[i] = (right_eq[i] - c3[index_i][index_j] * d_prime[i-1]) / (c1[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                else if(i == num_layers - 1)
                {d_prime[i] = (right_eq[i] - c3[index_i][index_j] * d_prime[i-1]) / (c2[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                
            }
            
            for (int i = num_layers - 1; i >= 0; i--){
                if (i == num_layers - 1)
                {vector_layer[i]->uMy_hat_->val_[index_i][index_j][0] = d_prime[i];}
                else
                {vector_layer[i]->uMy_hat_->val_[index_i][index_j][0] = d_prime[i] - c_prime[i] * vector_layer[i+1]->uMy_hat_->val_[index_i][index_j][0];}
                vector_layer[i]->uMy_hat_->val_[0][0][0] = 0;
            }
            
            delete c_prime;
            delete d_prime;
            delete right_eq;
            
        }
    }
    return 0;
}










int N_membrane::cal_uM_y_hat_ima(){
    
    for (int index_i = 0; index_i < Nx; index_i++){
        
        int index_ii = index_i;
        if (index_ii > Nx/2) {
            index_ii = index_ii - Nx;
        }
        
        double qi = 2.0*PI*index_ii / (Nx*Dx);
        
        for (int index_j = 0; index_j < (Ny/2 + 1); index_j++){
            double qj = 2.0*PI*index_j / (Ny*Dy);
            double q = sqrt(qi * qi + qj* qj);
            
            //for each i, j calculate layer_i_u_hat_mx[i][j]
            double* right_eq = new double[num_layers];
            
            //?????? needed to be further considered
            if( (index_ii == 0) && (index_j == 0) ){   //discuss the case that q_i = q_j = 0;
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = 0;
                }
            }
            //??????????????????????????
            else{
                for(int i = 0; i < num_layers; i++){
                    right_eq[i] = (1 - qj*qj/(q*q)) *vector_layer[i]->Wy_hat_->val_[index_i][index_j][1] + (- qi*qj/(q*q)) * vector_layer[i]->Wx_hat_->val_[index_i][index_j][1];
                }
            }
            
            // temporary array to reverse the matrix
            double* c_prime = new double[num_layers];
            double* d_prime = new double[num_layers];
            
            for (int i = 0; i < num_layers - 1; i++){
                if (i == 0)
                {c_prime[i] = c3[index_i][index_j] / c2[index_i][index_j];}
                else
                {c_prime[i] = c3[index_i][index_j] / (c1[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                
            }
            
            for (int i = 0; i < num_layers; i++){
                if (i == 0)
                {d_prime[i] = right_eq[i] / c2[index_i][index_j];}
                else if(i < num_layers - 1)
                {d_prime[i] = (right_eq[i] - c3[index_i][index_j] * d_prime[i-1]) / (c1[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                else if(i == num_layers - 1)
                {d_prime[i] = (right_eq[i] - c3[index_i][index_j] * d_prime[i-1]) / (c2[index_i][index_j] - c_prime[i-1]* c3[index_i][index_j]);}
                
            }
            
            for (int i = num_layers - 1; i >= 0; i--){
                if (i == num_layers - 1)
                {vector_layer[i]->uMy_hat_->val_[index_i][index_j][1] = d_prime[i];}
                else
                {vector_layer[i]->uMy_hat_->val_[index_i][index_j][1] = d_prime[i] - c_prime[i] * vector_layer[i+1]->uMy_hat_->val_[index_i][index_j][1];}
                vector_layer[i]->uMy_hat_->val_[0][0][1] = 0;
            }
            
            delete c_prime;
            delete d_prime;
            delete right_eq;
        }
    }
   
   
    return 0;
}
    


int main(int argc, char **argv) {
	//check for the correct number of input arguments
	/*if (argc != 5) {
		fprintf(stderr, "There should be 3 arguments, dt (double), Nstep (integer) and Number of Layers(int) \n");
		fprintf(stderr, "Usage: %s dt Nstep \n", argv[0]);
		return 1;
	}*/
	
    double dt_ = atof(argv[1]);
    int Nstep_ = atoi(argv[2]);
    int num_layers_ = atoi(argv[3]);
    double lambda = atof(argv[4]);
	    cout << lambda<<" is lambda"<<endl;
    int Nx = atoi(argv[5]);
    int Ny = atoi(argv[6]);
    int Nanalysis = Nstep_ / 20;
    double Dx = 0.5;
    double Dy = 0.5;
	    int a = 1;
	    int b = 1;
	    int w = 1;
	    int M = 1;
     double eta_M_ = atof(argv[7]);
     double eta_S_ = atof(argv[8]);
     int N_print_gap = atof(argv[9]);	
     int Nprint = Nstep_ / N_print_gap;
     double ave_psi = atof(argv[10]);  
     cout<< dt_<< " " <<Nstep_<< " "<<num_layers_<< " "<<lambda<< " "<<Nx<<" "<<Ny<<" "<<eta_M_<<" "<<eta_S_<<" "<<N_print_gap<<" "<<ave_psi<<" "<< endl;

//double lambda = 0.07;
	    double thick_sol_ = 1;
	    N_membrane *system = new N_membrane(dt_, Nstep_, num_layers_, Nx,Ny, Dx,Dy,a,b,w,M,eta_M_,eta_S_,lambda,0,thick_sol_, ave_psi);
	    /*N_membrane(dt_, Nstep_, num_layers_, Nx_, Ny_, Dx_, Dy_, a_, b_, w_, M_, eta_M_, eta_S_,
		       lamda_, gamma_, thick_sol_):dt(dt_), Nstep(Nstep_), num_layers(num_layers_)*/
		
		
	    //vector<char*> name_psi;
	    //vector<char*> name_uM_x;
	    //vector<char*> name_uM_y;
		char filename1[60], filename2[60], filename3[60];
	    
	    FILE* fp_psi;
	    FILE* fp_uM_x;
	    FILE* fp_uM_y;
    FILE* fp_length;
    FILE* fp_record;
   /* vector<FILE*> vector_file_psi;
    vector<FILE*> vector_file_uM_x;
    vector<FILE*> vector_file_uM_y;
*/
    
	/*for (int i = 0; i < num_layers; i++){
        OneLayer* layer_i = new OneLayer(Nx, Ny, Dx, Dy, a, b, w, lambda, eta_M, eta_S, eta_M_s, eta_S_s, gamma, M, i, thick_sol);
		layer_i->init_psi_Gauss(0.0, 0.05);
        vector_layer.push_back( layer_i ) ;
	}*/
	
    
    //debug
    //cout<< system->vector_layer[1]->psi_->val_[10][10]<< "initialize" <<endl;
    
    char filename4[50], filename5[50];
    
    for(vector<OneLayer*>::size_type p = 0; p < system->vector_layer.size(); p++){
        
        sprintf(filename4, "./data/membrane_%u_file_length", p);
        sprintf(filename5, "./data/membrane_%u_file_record", p);
        
        if ( (fp_length = fopen(filename4,"a")) == NULL) {
            printf("can't open file\n");
            exit(0);
        }
        
        if ( (fp_record = fopen(filename5,"a")) == NULL) {
            printf("can't open file\n");
            exit(0);
        }
        fprintf(fp_record, "ave(psi) max_abs(psi) ave(uMx) max_abs(uMx) ave(uMy) max_abs(uMy)\n");
        fclose(fp_length);
        fclose(fp_record);

    }
    //time
    //clock_t tStart = clock();
//debug
/*
double* test = new double[5];
test[0] = 0.1; test[1] = 0.2; test[2] = 0.3; test[3] = 0.4;test[4] = 0.5;
double* ans = system->reverse_tridiag(5, test, 1, 2, 3);
for (int i = 0; i < 5; i++){ 
cout<<"answer is "<< ans[i]<<endl;
}*/

	clock_t tStart = clock();
		for (int n = 0; n <= system->Nstep; n++) {
            
            
            /*
            if (n == 0){
                char filename_psi[50] = "psi_at_0";
                FILE* fp_psi = fopen(filename_psi, "w");
                for (int i = 0; i < Nx; i++){
                    for (int j = 0; j < Ny; j++){
                        
                        fprintf(fp_psi, "%lf \n", system->vector_layer[0]->psi_->val_[i][j]);
                    }
                }
                fclose(fp_psi);
            } 
*/
            
            
            
            system->cal_mu();
  /*          //debug:output
             if (n == 0){ 
		char filename_mu[50] = "mu_at_0";
		FILE* fp_mu = fopen(filename_mu, "w");
		for (int i = 0; i < Nx; i++){ 
			for (int j = 0; j < Ny; j++){ 
			
		fprintf(fp_mu, "%lf \n", system->vector_layer[0]->mu_->val_[i][j]);
}
}
		fclose(fp_mu);
		} 
            //debug
          //  cout<< system->vector_layer[1]->psi_->val_[10][10]<<"psi is "<<endl;
    */        
            system->cal_force();
            //debug
            //cout<< system->vector_layer[1]->Wx_hat_->val_[10][10]<<"Wx is"<<endl;
 
            //for testing
           // tStart = clock();
           

            omp_set_num_threads(4);
            #pragma omp parallel num_threads(4)
            {
                int ID = omp_get_thread_num();
    
                if(ID == 0){
                // these funtions calculate the SYSTEM(all the membrane)'s velocity field
                    system->cal_uM_x_hat_real();
                //cout<< system->vector_layer[1]->uMx_hat_->val_[10][10][0]<<"uMx is"<<endl;
            }
            else if (ID == 1){
                system->cal_uM_x_hat_ima();
            }
            else if (ID == 2){
                system->cal_uM_y_hat_real();
            }
            else { system->cal_uM_y_hat_ima();}
            }


 /* 
  
  ////using oseen tensor to calculate the velocity
	system->vector_layer[0]->cal_T_hat();
	system->vector_layer[1]->cal_T_hat();
	system->vector_layer[0]->cal_W_hat();
	system->vector_layer[1]->cal_W_hat();
	system->vector_layer[0]->cal_uM_hat(system->vector_layer[1]->Wx_hat_, system->vector_layer[1]->Wy_hat_);
	system->vector_layer[1]->cal_uM_hat(system->vector_layer[0]->Wx_hat_, system->vector_layer[0]->Wy_hat_);
   
            //debug:output
            if (n == 0){
                char filename_uM[50] = "uM_hat_at_0_real";
                FILE* fp_uM = fopen(filename_uM, "w");
                for (int i = 0; i < Nx; i++){
                    for (int j = 0; j < Ny; j++){
                        
                        fprintf(fp_uM, "%lf \n", system->vector_layer[0]->uMx_hat_->val_[i][j][0]);
                    }
                }
                fclose(fp_uM);
            }
            
            if (n == 0){
                char filename_force[50] = "Wx_hat_at_0_real";
                FILE* fp_force = fopen(filename_force, "w");
                for (int i = 0; i < Nx; i++){
                    for (int j = 0; j < Ny; j++){
                        
                        fprintf(fp_force, "%lf \n", system->vector_layer[0]->Wx_hat_->val_[i][j][0]);
                    }
                }
                fclose(fp_force);
            } 

   */         
            
 //////////////////////////////////////////////////////////////////
 //These commands will generate a file recording the driving force W_hat, the velocity uM_hat using different methods.
/* 
if(n == 0){   
char record_filename1[50] = "./uMx_hat_real_part_using_reverse_matrix_at_0";
   FILE* fp_test = fopen(record_filename1, "w");
   cout << "Nx is "<< Nx;
	for (int i = 0; i < Nx; i++){ 
	for (int j = 0; j < (Ny/2 + 1 ); j++){ 
   	fprintf(fp_test,"%lf \n", system->vector_layer[0]->uMx_hat_->val_[i][j][0]);
	}
	}
	fclose(fp_test);

char record_filename2[50] = "./uMy_hat_real_part_using_reverse_matrix_at_0";
     fp_test = fopen(record_filename2, "w");
         for (int i = 0; i < Nx; i++){
         for (int j = 0; j < (Ny/2 + 1 ); j++){
         fprintf(fp_test,"%lf \n", system->vector_layer[0]->uMy_hat_->val_[i][j][0]);
         }
         }
         fclose(fp_test);

char record_filename3[50] = "./uMx_hat_ima_part_using_reverse_matrix_at_0";
     fp_test = fopen(record_filename3, "w");
         for (int i = 0; i < Nx; i++){
         for (int j = 0; j < (Ny/2 + 1 ); j++){
         fprintf(fp_test,"%lf \n", system->vector_layer[0]->uMx_hat_->val_[i][j][1]);
         }
         }
         fclose(fp_test);

char record_filename4[50] = "./uMy_hat_ima_part_using_reverse_matrix_at_0";
     fp_test = fopen(record_filename4, "w");
         for (int i = 0; i < Nx; i++){
         for (int j = 0; j < (Ny/2 + 1 ); j++){
         fprintf(fp_test,"%lf \n", system->vector_layer[0]->uMy_hat_->val_[i][j][1]);
         }
         }
         fclose(fp_test);

} 
*/
// cout<<"Time for cal_velocity(matriex inverse) taken"<< (double)(clock() - tStart)<<endl;
            //for testing
           // tStart = clock();
            
            //system->vector<OneLayer*>::size_type j = 0;
            for (int i = 0; i < system->num_layers; i++){
                fftw_execute(system->vector_layer[i]->p_out_uMx_);
                fftw_execute(system->vector_layer[i]->p_out_uMy_);
                for (int p = 0; p < Nx; p++) {
                    for (int q = 0; q < Ny; q++) {
                        system->vector_layer[i]->uMx_->val_[p][q] = system->vector_layer[i]->uMx_->val_[p][q] / (Nx * Ny);
                        system->vector_layer[i]->uMy_->val_[p][q] = system->vector_layer[i]->uMy_->val_[p][q] / (Nx * Ny);
                    }
                }
            }

            system->cal_rhs();
          //  cout<<system->vector_layer[1]->rhs_hat_->val_[10][10][1]<<"rhs is "<<endl;

			//print the result if every Nprint steps
			if ((n % Nprint) == 0) {
                for ( vector<OneLayer*>::size_type p = 0; p < system->vector_layer.size(); p++) {
              
                    sprintf(filename1, "./data/file_membrane_%u_psi", p);
                    sprintf(filename2, "./data/file_membrane_%u_uM_x", p);
                    sprintf(filename3, "./data/file_membrane_%u_uM_y", p);
			    //create files
                    if ((fp_psi = fopen(filename1,"a")) == NULL) {
                        printf("cannot open file\n");
                        exit(0);
                    }
                    
                    if ((fp_uM_x = fopen(filename2,"a")) == NULL) {
                        printf("cannot open file\n");
                        exit(0);
                    }
				
                    if ((fp_uM_y = fopen(filename3,"a")) == NULL) {
                        printf("cannot open file\n");
                        exit(0);
                    }
                    //print out the time step
                    fprintf(fp_psi, "t, %d\n",n);
		    fprintf(fp_uM_x, "t, %d\n",n);
  		    fprintf(fp_uM_y, "t, %d\n",n);
                    //output data
                    for (int i = 0; i < Nx; i++) {
                        for (int j = 0; j < Ny; j++) {
                            fprintf(fp_psi, "%lf\n", system->vector_layer[p]->psi_->val_[i][j]);
                            fprintf(fp_uM_x, "%lf\n", system->vector_layer[p]->uMx_->val_[i][j]);
                            fprintf(fp_uM_y, "%lf\n", system->vector_layer[p]->uMy_->val_[i][j]);
                            }
                        }
			
                    double time = (double)(clock()-tStart)/(double)CLOCKS_PER_SEC;
 		    cout << "Time for each step is "<< time<<endl;
                    //close the file
                    fclose(fp_psi);
                    fclose(fp_uM_x);
                    fclose(fp_uM_y);
                }
            }
            
            
                
                if ((n % Nanalysis) == 0) {
                    for(vector<OneLayer*>::size_type p = 0; p < system->vector_layer.size(); p++){
                        
                        sprintf(filename4, "./data/membrane_%u_file_length", p);
                        sprintf(filename5, "./data/membrane_%u_file_record", p);
                        
                        if ( (fp_length = fopen(filename4,"w")) == NULL) {
                            printf("can't open file\n");
                            exit(0);
                        }
                        
                        if ( (fp_record = fopen(filename5,"w")) == NULL) {
                            printf("can't open file\n");
                            exit(0);
                        }
                        fseek(fp_length,0,2);
                        fseek(fp_record,0,2);
                        fprintf(fp_length,"%lf\n", length_R(system->vector_layer[p]->psi_, Dx));
                        fprintf(fp_record,"%lf %lf %lf %lf %lf %lf\n",
                                cal_ave(system->vector_layer[p]->psi_), cal_max_abs(system->vector_layer[p]->psi_));
                              //  cal_ave(system->vector_layer[p]->uMx_), cal_max_abs(system->vector_layer[p]->uMx_),
                              //  cal_ave(system->vector_layer[p]->uMy_), cal_max_abs(system->vector_layer[p]->uMy_));
                    
                    printf("%lf %lf \n",
                           cal_ave(system->vector_layer[p]->psi_), cal_max_abs(system->vector_layer[p]->psi_));
                    // cal_ave(system->vector_layer[p]->uMx_), cal_max_abs(system->vector_layer[p]->uMx_),
                    //       cal_ave(system->vector_layer[p]->uMy_), cal_max_abs(system->vector_layer[p]->uMy_));
                    
                    fclose(fp_length);
                    fclose(fp_record);
                        
                    }

                }
            
    
 
       
		//calculate psi_hat
            for (vector<OneLayer*>::size_type i = 0; i < system->vector_layer.size(); i++) {
                system->vector_layer[i]->psi_to_psi_hat();
            }
           
            
		//update psi_hat by Euler Forward method
            for (vector<OneLayer*>::size_type i = 0; i < system->vector_layer.size(); i++ ){
                for (int p = 0; p < system->vector_layer[i]->psi_hat_->nrow_; p++) {
                        for (int q = 0; q < system->vector_layer[i]->psi_hat_->ncol_; q++) {
                            system->vector_layer[i]->psi_hat_->val_[p][q][0] += dt_*system->vector_layer[i]->rhs_hat_->val_[p][q][0];
                            system->vector_layer[i]->psi_hat_->val_[p][q][1] += dt_*system->vector_layer[i]->rhs_hat_->val_[p][q][1];
                        }
                }
            }
/*
            //debug:output
            if (n == 0){
                char filename_rhs[50] = "rhs_0";
                FILE* fp_rhs = fopen(filename_rhs, "w");
                for (int i = 0; i < Nx; i++){
                    for (int j = 0; j < Ny; j++){
                        
                        fprintf(fp_rhs, "%lf \n", system->vector_layer[0]->rhs_hat_->val_[i][j][0]);
                                }
                                }
                                fclose(fp_rhs);
                                }
*/
            
            for(vector<OneLayer*>::size_type i = 0; i < system->vector_layer.size(); i++){
                    //calculate psi by inverse DFT
                    system->vector_layer[i]->psi_hat_to_psi();
            }
	
            //time
       //     double time = (double)(clock()-tStart);
       //     cout<<"Time for fftw the rest is"<< time<<endl;
            
        
	}
			
	//delete layer1;
	
	return 0;
}



    
    
    

    
    

