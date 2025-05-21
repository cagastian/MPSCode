#ifndef WESSEL_ULA_LAPACK
#define WESSEL_ULA_LAPACK

#include <ula/matrixtypes.h>
#include <complex>
#include <string>


extern "C" void dsaupd_(int* ido, char* bmat, int* n, char* which, int* nev,
			double* tol, double* resid, int* ncv, double** v,
			int* ldv, int* iparam, int* ipntr, double*  workd,
			double* workl, int* lworkl, int *info);
extern "C" void dseupd_(bool* rvec, char* howmny, bool* select, double** d,
			double** z, int* ldz, double* sigma, char* bmat, int* n,
			char* which, int* nev, double* tol, double* resid,
			int* ncv, double** v, int* ldv, int* iparam, int* ipntr,
			double* workd, double* workl, int* lworkl, int* ierr );
//v, resid, iparam, ipntr, workl, workd : are vectors
//SSAUPD computes the quantities needed to construct the desired eigenvalues and eigenvectors(if requested)
//SSEUPD extracts the desired eigenvalues and eigenvectors: rvec=true eigenvectors are computed
// HOWMNY  Specifies how many Ritz vectors are wanted and the form of Zthe matrix of Ritz vectors. See remark 1 below. = 'A': compute NEV Ritz vectors;  = 'S': compute some of the Ritz vectors, specified by the logical array SELECT.

namespace ula{
  
 void  TriDiag_MV(int dim, double* X, double* Y){
  //Problema PARTICULAR 
  //     The matrix used is the 2 dimensional discrete Laplacian on unit
  //     square with zero Dirichlet boundary condition. 
   //Tridiagnoal matrix-Vector multiplication
    double* X_ = new double [dim];
    double* Y_ = new double [dim];
    for(int j=0; j<dim; j++){
      X_[j]=X[j];
      Y_[j]=Y[j];
    }
    Real DD = 4.0;
    Real DL = -1.0;
    Real DU = - 1.0;
    Y_[0] = DD*X_[0] + DU*X_[1];
    for(int j=1; j<dim-1; j++){
      Y_[j] = DL*X_[j-1] + DD*X_[j] + DU*X_[j+1]; 
    }
    Y_[dim-1] = DL*X_[dim-2] + DD*X_[dim-1];
    
    for(int j=0; j<dim; j++){
      X[j]=X_[j];
      Y[j]=Y_[j];
    }
    delete[] X_;
    delete[] Y_;
  }
  
  void svpv(int dim, double scalar,double* v, double* y){
    //svpv: Scalar_Vector_plus_Vector y=y+scalar*v
    double* v_ = new double [dim];
    double* y_ = new double [dim];
    for(int j=0; j<dim; j++){
      v_[j]=v[j];
      y_[j]=y[j];
    }
    for(int j=0; j<dim; j++)
      y_[j] += scalar*v_[j]; 
    
    for(int j=0; j<dim; j++)
      y[j]=y_[j];
    
    delete[] v_;
    delete[] y_;
  }
  double norm (int dim, double* v){
    double norm=0.;
    for(int i=0; i<dim; i++)
      norm+= v[i]*v[i];
    return sqrt(norm);
  }
  
  RealVector Matrix_Vec(int nx, RealVector &X){
    //     Compute the matrix vector multiplication Y=T*X
    //     where T is a nx by nx tridiagonal matrix with DD on the 
    //     diagonal, DL on the subdiagonal, and DU on the superdiagonal.
    int dim=nx*nx;
    RealVector X_(dim);
    RealVector Y(dim);
    X_ = X;
    Real DD = 4.0;
    Real DL = -1.0;
    Real DU = - 1.0;
    Y(0) = DD*X(0) + DU*X(1);
    for(int j=1; j<dim-1; j++){
      Y(j) = DL*X(j-1) + DD*X(j) + DU*X(j+1); 
    }
    Y(dim-1) = DL*X(dim-2) + DD*X(dim-1);
    
    return Y;
  }
  
  
   

  void arnoldi_sim(int num_ev, int basis, RealVector& Evals){
    //basis: basis elements
    //num_ev   : elements (values or vectors) to be computed 
    //%--------------------------------------%
    //           | Perform matrix vector multiplication |
    //           |              y <--- OP*x             |
    //           | The user should supply his/her own   |
    //           | matrix vector multiplication routine |
    //           | here that takes workd(ipntr(1)) as   |
    //           | the input, and return the result to  |
    //           | workd(ipntr(2)).                     |
    //           %--------------------------------------%
    
    //     | MAXN:   Maximum dimension of the A allowed.
    //     | MAXNEV: Maximum NEV allowed.                 
    //     | MAXNCV: Maximum NCV allowed
    //     | NEV <= MAXNEV
    //     | NEV + 1 <= NCV <= MAXNCV
    
    //options for which: SM, LA, SA, LI, SI
    //'LA' - largest (algebraic) eigenvalues.
    //'SA' - smallest (algebraic) eigenvalues.
    //'LM' - largest (in magnitude) eigenvalues.
    //'SM' - smallest (in magnitude) eigenvalues.
    //'BE' - half from each end of the spectrum.
    
    //________________________________________________
    int maxn=basis, maxncv=num_ev, maxnev= maxncv/2;          //|
    int ncv = maxncv, nev = ncv/2 , ldv = basis;// ldz = basis;  //|
    maxn=256; maxnev=10; maxncv=25; ldv=maxn;
    int nx=10; nev=4; ncv=20;
    int n = nx*nx;//nev*nev;
    //________________________________________________           
    if(n>maxn){
      std::cout<<"Error: wrong dimension nx*nx is grater than maxn"<<std::endl;
    }
    if(nev>maxnev){
      std::cout<<"Error: wrong dimension nev is grater maxnev "<<std::endl;
    }
    if(ncv>maxncv){
      std::cout<<"Error: wrong dimension ncv is grater maxncv "<<std::endl;
      std::cout<<"basis="<<basis<<" maxncv="<<maxncv<<std::endl;
    }
    //--------ssaupd--------
    int ido = 0;
    char bmat = 'I'; //El sistema A*x=lambdda*x, con B=Matriz Identidad
  
    char which[] = {'L','M'}; 
    double zero = 0.0e0;
    double tol = zero;
    
    double* resid = new double[maxn];
    double** v = new double* [ldv];
    for(int i=0; i<ldv; i++)
      v[i] = new double [ncv];
    int* iparam = new int [11];
    int ishfrt = 1;  //0; ???
    int maxitr = 300;
    int model = 1;
    iparam[0] = ishfrt;
    iparam[2] = maxitr;
    iparam[6] = model;
    int* ipntr = new int [11];
    double* workd = new double [3*maxn];
    double* workl = new double [maxncv*(maxncv+8)];
    int lworkl = ncv*(ncv+8);
    int info=0;
    //----------------------
    //--------sseupd--------
    bool rvec;
    char howmny = 'A';
    bool* select = new bool [maxncv];
    //double** z = new double* [ldz];
    //for(int i=0; i<ldz; i++)
    //  z[i] = new double [maxncv];
    double **d = new double *[maxncv];
    for(int i=0; i<maxncv; i++)
      d[i]=new double [2];
    double sigma;
    int ierr;
    int nconv;
    double* ax = new double [maxn];
    //double* p = new double [maxn];
    //----------------------
    
    for( int i=0; i<3*maxn; i++){
      workd[i]=1.;
    }
    
    do{
      dsaupd_(&ido, &bmat, &basis, which, &nev, &tol, resid, &ncv, v, &ldv,
	      iparam, ipntr, workd, workl, &lworkl, &info);
      

      
      //Particular problem to solve
      TriDiag_MV(maxn,workd+ipntr[0]-1,workd+ipntr[1]-1); 
      
      
    }while(ido== -1 || ido == 1);
    if(info<0)
      std::cout<<" Error with ssaupd, info = "<<info<<std::endl;
    else{
      rvec=true;
      
      dseupd_(&rvec, &howmny, select, d, v, &ldv, &sigma, &bmat, &n,
	      which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr,
	      workd, workl, &lworkl, &ierr );
      //Eigenvalues are returned in the first column |
      //of the two dimensional array D and the       |
      //corresponding eigenvectors are returned in   |
      //the first NCONV (=IPARAM[4]) columns of the  |
      //two dimensional array V if requested.        |
      if (ierr != 0)
	std::cout<<"Error with sseupd, info = "<<ierr<<std::endl;
      else{
	//Calc of the resudual norm-------------
	//A*v=a*v Problem to solve
	nconv = iparam[4];
	
	for(int j=0; j<nconv; j++){
	  //ax=A*v_cal
	  
	  TriDiag_MV(maxn, v[0]+j, ax);


	  //ax=ax-a*v_cal
	  svpv(maxn, d[j][0], v[0]+j, ax);
	  //||ax||
	  d[j][1]=norm(maxn,ax);
	  //residual
	  d[j][1]/=abs(d[j][0]);
	}
	
	for(int j=0; j<nconv; j++){
	  std::cout<<j<<" "<<d[j][1]<<std::endl;
	}	  
      }
    }
    for(int i=0; i<num_ev; i++)
      Evals(i)=d[i][0];
  }
  
 
}





#endif
