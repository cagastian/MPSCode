#ifndef WESSEL_ULA_LAPACK
#define WESSEL_ULA_LAPACK

#include <ula/matrixtypes.h>
#include <complex>

namespace ula {
  
  void _message( int info, int lwork ) {
    if( info > 0 ) 
      std::cerr << "Warning: diagonalization failed to converge." << std::endl;
    else 
      if( info < 0 ) 
	std::cerr << "Warning: illegal argument value in diagonalization." << std::endl;
#ifdef BL_DEBUG
      else std::cerr << "Optimal cache size is lwork = " << lwork << "." << std::endl;
#endif
    
  }
}

/*
  ---------------------------------------------------
  
  DIAGONALIZATION OF REAL NON-SYMMETRIC MATRIX
  
  H: n-by-n input matrix
  E: eigenvalues
  
  ---------------------------------------------------
*/

extern "C" void dgeev_(char* jobvl, char* jobvr, 
		       int* n, double* a, int* lda, 
		       double* er, double* ei, 
		       double* vl, int* ldvl, double* vr, int* ldvr,
                       double* work, int* lwork, int* info);

extern "C" void dgeevx_(char* balanc, 
                        char* jobvl, char* jobvr, 
                        char* sense,
		        int* n, double* a, int* lda, 
		        double* er, double* ei, 
		        double* vl, int* ldvl, double* vr, int* ldvr,
                        int* ilo, int* ihi, double* scale, double* abnrm, double* rconde, double* rcondv,
                        double* work, int* lwork, 
                        int* iwork,
                        int* info);

namespace ula {
  
  void ns_diag( RealMatrix& H, ComplexVector& E) {
    RealMatrix H_=H;
    int calc_evl=0;
    int calc_evr=0;
    char jobvl = ( calc_evl ? 'V' : 'N' );  // calculate left  eigenvectors?
    char jobvr = ( calc_evr ? 'V' : 'N' );  // calculate right eigenvectors?
    int n = H.size1();                      // squareness is not checked
    int lda = n;                            // leading dimension;
    double* er = new double[n];             // real part of eigenvalues
    double* ei = new double[n];             // imag part of eigenvalues
    int ldvl=1;
    double* vl = new double[ldvl];
    int ldvr=1;
    double* vr = new double[ldvr];
    int lwork=-1;                          //  LWORK >= max(1,3*N)                                    
    double* work_ = new double[1];
    int info;
    dgeev_( &jobvl, &jobvr, 
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            work_, &lwork, &info );  
    _message( info, lwork );
    lwork = (int)work_[0];           
    double* work = new double[lwork];
    dgeev_( &jobvl, &jobvr, 
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            work, &lwork, &info ); 
    _message( info, lwork );
    
    for (int i=0;i<n;++i)
      E(i)=std::complex<double>(er[i],ei[i]);
    delete[] er;
    delete[] ei;
    delete[] vl;
    delete[] vr;
    delete[] work;
    delete[] work_;
  }

  void ns_diag( RealMatrix& H, ComplexVector& E, ComplexMatrix& V) {
    RealMatrix H_=H;
    int calc_evl=0;
    int calc_evr=1;
    char jobvl = ( calc_evl ? 'V' : 'N' );  // calculate left  eigenvectors?
    char jobvr = ( calc_evr ? 'V' : 'N' );  // calculate right eigenvectors?
    int n = H.size1();                      // squareness is not checked
    int lda = n;                            // leading dimension;
    double* er = new double[n];             // real part of eigenvalues
    double* ei = new double[n];             // imag part of eigenvalues
    int ldvl=1;
    double* vl = new double[ldvl*ldvl];
    int ldvr=n;
    double* vr = new double[ldvr*ldvr];
    int lwork=-1;                         //  if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N. (5N was used)
    double* work_ = new double[1];
    int info;
    dgeev_( &jobvl, &jobvr, 
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            work_, &lwork, &info );  
    _message( info, lwork );
    lwork = (int)work_[0];           
    double* work = new double[lwork];
    dgeev_( &jobvl, &jobvr, 
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            work, &lwork, &info ); 
    _message( info, lwork );
    for (int i=0;i<n;++i) {
      if (ei[i]) {
        E(i)=std::complex<double>(er[i],ei[i]);
        E(i+1)=std::complex<double>(er[i+1],ei[i+1]);
        for (int j=0;j<n;++j) {
          V(j,i)=std::complex<double>(vr[n*i+j],vr[n*(i+1)+j]);
          V(j,i+1)=std::complex<double>(vr[n*i+j],-vr[n*(i+1)+j]);
        }
        ++i;
      }
      else {
        E(i)=std::complex<double>(er[i],ei[i]);
        for (int j=0;j<n;++j)
          V(j,i)=vr[n*i+j];
      }
    }   
    delete[] er;
    delete[] ei;
    delete[] vl;
    delete[] vr;
    delete[] work;
    delete[] work_;
  }

  void ns_diagx( RealMatrix& H, ComplexVector& E, ComplexMatrix& V) {
    RealMatrix H_=H;
    int calc_evl=1;
    int calc_evr=1;
    char balanc= 'B';
    char jobvl = ( calc_evl ? 'V' : 'N' );  // calculate left  eigenvectors?
    char jobvr = ( calc_evr ? 'V' : 'N' );  // calculate right eigenvectors?
    char sense = 'B';
    int n = H.size1();                      // squareness is not checked
    int lda = n;                            // leading dimension;
    double* er = new double[n];             // real part of eigenvalues
    double* ei = new double[n];             // imag part of eigenvalues
    int ldvl=n;
    double* vl = new double[ldvl*ldvl];
    int ldvr=n;
    double* vr = new double[ldvr*ldvr];
    int ilo;
    int ihi;
    double* scale=new double[n];
    double abnrm;
    double* rconde=new double[n];
    double* rcondv=new double[n];
    int lwork=-1;                           //n*(n+6) was used
    double* work_ = new double[1];
    int* iwork = new int[2*n-2];
    int info;
    dgeevx_( &balanc,
            &jobvl, &jobvr, 
            &sense,
	    &n, &(H_(0,0)), &lda, 
	    er, ei,
	    vl, &ldvl, vr, &ldvr,
            &ilo, &ihi , scale, &abnrm, rconde, rcondv,
	     work_, &lwork, 
	     iwork,
	     &info );  
    _message( info, lwork );
    lwork = (int)work_[0];           
    double* work = new double[lwork];
    dgeevx_( &balanc,
            &jobvl, &jobvr, 
	     &sense,
	     &n, &(H_(0,0)), &lda, 
	     er, ei,
	     vl, &ldvl, vr, &ldvr,
	     &ilo, &ihi , scale, &abnrm, rconde, rcondv,
	     work, &lwork, 
	     iwork,
	     &info );  
    _message( info, lwork );
    
    for (int i=0;i<n;++i) {
      if (ei[i]) {
        E(i)=std::complex<double>(er[i],ei[i]);
        E(i+1)=std::complex<double>(er[i+1],ei[i+1]);
        for (int j=0;j<n;++j) {
          V(j,i)=std::complex<double>(vr[n*i+j],vr[n*(i+1)+j]);
          V(j,i+1)=std::complex<double>(vr[n*i+j],-vr[n*(i+1)+j]);
        }
        ++i;
      }
      else {
        E(i)=std::complex<double>(er[i],ei[i]);
        for (int j=0;j<n;++j)
          V(j,i)=vr[n*i+j];
      }
    }   
    delete[] er;
    delete[] ei;
    delete[] vl;
    delete[] vr;
    delete[] work;
    delete[] work_;
  }
}

/*
  ---------------------------------------------------
  
  SINGULAR VALUE DECOMPOSITION OF A COMPLEX MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

extern "C" void zgesvd_(char *jobu, char *jobvt, int *m, int *n, 
			std::complex<double> *A, int *lda, double *S, 
			std::complex<double> *U, int *ldu, 
			std::complex<double> *vt, int *ldvt, 
			std::complex<double> *work, int *lwork, 
			double *rwork, int *info);

namespace ula {

  void svd( ComplexMatrix& A_, ComplexMatrix& U, RealVector& D, 
	    ComplexMatrix& V){
    char jobu ='S';
    char jobvt='S';
    int m=A_.size1();
    int n=A_.size2();
    ComplexMatrix A(n,m);
    A=A_;
    int lda=m;
    int ldu=U.size1();
    int d;
    int ldvt=V.size1();
    if(m<n)d=m;
    else d=n;
    std::complex<double>* work1 = new std::complex<double>[1];
    int lwork=-1;
    double* rwork = new double[5*d];
    int info;
    zgesvd_(&jobu, &jobvt, &m, &n, 
	    &(A(0,0)), &lda, &(D(0)), 
	    &(U(0,0)),  &ldu, 
	    &(V(0,0)), &ldvt,
	    work1, &lwork,
	    rwork, &info);
    lwork=(int)real(work1[0]);
    delete[] work1;
    std::complex<double>* work = new std::complex<double>[lwork];
    zgesvd_(&jobu, &jobvt, &m, &n, 
	    &(A(0,0)), &lda, &(D(0)), 
	    &(U(0,0)),  &ldu, 
	    &(V(0,0)), &ldvt,
	    work, &lwork,
	    rwork, &info);
    if(info!=0)std::cout<<"INFO="<<info<<"\n";
    delete[] rwork;
    delete[] work;
  }
}


/*
  ---------------------------------------------------
  
  DIAGONALIZATION OF REAL SYMMETRIC MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

extern "C" void dsyev_(char* jobz, char* uplo, 
                       int* n, double* a, int* lda, 
                       double* w, 
                       double* work, int* lwork,
                       int* info);

namespace ula {
  
  void _diag( RealMatrix& H, RealVector& E, bool calc_ev ) {
    char jobz = ( calc_ev ? 'V' : 'N' );  // calculate eigenvectors?
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    int lda = n;                          // leading dimension;
    int lwork = 3*n-1;
    double* work = new double[lwork];
    int info;
    dsyev_( &jobz, &uplo, 
	    &n, &(H(0,0)), &lda, 
	    &(E(0)), 
	    work, &lwork, 
	    &info);  
    _message( info, lwork );
    delete[] work;
  }
  
  
  inline void diag( RealMatrix& H, RealVector& E ) {
    _diag( H, E, false );
  }
  
  inline void diag( RealMatrix& H, RealVector& E, RealMatrix& V ) {
    V=H;
    _diag( V, E, true );
  }
}

/*
  ---------------------------------------------------
  
  INVERSE OF REAL MATRIX
  
  A : n-by-n input matrix
  A_: A^{-1}
  
  ----------------------------------------------------
*/

extern "C" void dgesv_(int *n, int *nrhs, 
		       double *a, int *lda,
		       int *ipiv,
		       double *b, int *ldb, 
		       int* info );

namespace ula {
  
  void invr(RealMatrix& A, RealMatrix& A_) {
    int n = A.size1();
    RealMatrix B(n,n);
    int lda = n;
    int ldb = n;
    int info;
    int nrhs=n;
    int* ipiv = new int[n];
    A_=A;
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++)
	if(i==j)B(i,i)=1;
	else B(i,j)=0;
    dgesv_( &n, &nrhs ,
	    &(A_(0,0)), &lda, 
	    ipiv,
	    &(B(0,0)), &ldb, 	    
	    &info );
    if(info==0)A_=B;
    else std::cout<<"Info="<<info<<std::endl;    
  }
}
/*
  ---------------------------------------------------
  
  INVERSE OF COMPLEX MATRIX
  
  A : n-by-n input matrix
  A_: A^{-1}
  
  ----------------------------------------------------
*/

extern "C" void zgesv_(int *n, int *nrhs, 
		       std::complex<double> *a, int *lda,
		       int *ipiv,
		       std::complex<double> *b, int *ldb, 
		       int* info );

namespace ula {
  
  void invr(ComplexMatrix& A, ComplexMatrix& A_) {
    int n = A.size1();
    ComplexMatrix B(n,n);
    int lda = n;
    int ldb = n;
    int info;
    int nrhs=n;
    int* ipiv = new int[n];
    A_=A;
    for(int i=0;i<n;i++)
      B(i,i)=1; 
    zgesv_( &n, &nrhs ,
	    &(A_(0,0)), &lda, 
	    ipiv,
	    &(B(0,0)), &ldb, 	    
	    &info );
    if(info==0)A_=B;
    else {
      std::cout<<"Info="<<info<<std::endl;    
      exit(1);
    }
  }
}

/*
  ---------------------------------------------------
  
  FORWARD FORIER TRANSFORM
  
  Ft: Initial vector
  Fw: FFT[Ft]

  ---------------------------------------------------
*/
/*
extern "C" void dffti_(int *n, double *wsave);

extern "C" void dfftf_(int *n, double *r, double *wsave);

extern "C" void zffti_(int *n, double *wsave);

extern "C" void zfftf_(int *n, std::complex<double> *r, double *wsave);

namespace ula {
  
  void FFT(RealVector& Ft,RealVector& Fw){
    int n=Ft.size();
    double* wsave = new double[2*n+15];
    Fw=Ft;
    dffti_(&n,wsave);
    dfftf_(&n,&(Fw(0)),wsave);
    delete[] wsave;
  }

  void FFT(ComplexVector& Ft,ComplexVector& Fw){
    int n=Ft.size();
    double* wsave = new double[4*n+15];
    Fw=Ft;
    zffti_(&n,wsave);
    zfftf_(&n,&(Fw(0)),wsave);
    delete[] wsave;
  }

  void FFT(ComplexMatrix& T,ComplexMatrix& W,bool d){
    int r=T.size1(),c=T.size2(),i,j;
    if(d){
      ComplexVector v(c),v_(c); 
      for(i=0;i<r;i++){
	for(j=0;j<c;j++)
	  v(j)=T(i,j);
	FFT(v,v_);
	for(j=0;j<c;j++)
	  W(i,j)=v_(j);	
      }
    }else{
      ComplexVector v(r),v_(r); 
      for(i=0;i<c;i++){
	for(j=0;j<r;j++)
	  v(j)=T(j,i);
	FFT(v,v_);
	for(j=0;j<c;j++)
	  W(j,i)=v_(j);	
      }
    }
  }
}
*/
/*
  ---------------------------------------------------
  
  INVERSE FORIER TRANSFORM
  
  Ft: Initial vector
  Fw: FFT[Ft]

  ---------------------------------------------------
*/
/*
extern "C" void dffti_(int *n, double *wsave);

extern "C" void dfftb_(int *n, double *r, double *wsave);

extern "C" void zffti_(int *n, double *wsave);

extern "C" void zfftb_(int *n, std::complex<double> *r, double *wsave);

namespace ula {
  
  void IFT(RealVector& Ft,RealVector& Fw){
    int n=Ft.size();
    double* wsave = new double[2*n+15];
    Fw=Ft;
    dffti_(&n,wsave);
    dfftb_(&n,&(Fw(0)),wsave);
    delete[] wsave;
  }

  void IFT(ComplexVector& Ft,ComplexVector& Fw){
    int n=Ft.size();
    double* wsave = new double[4*n+15];
    Fw=Ft;
    zffti_(&n,wsave);
    zfftb_(&n,&(Fw(0)),wsave);
    delete[] wsave;
  }

  void IFT(ComplexMatrix& T,ComplexMatrix& W,bool d){
    int r=T.size1(),c=T.size2(),i,j;
    if(d){
      ComplexVector v(c),v_(c); 
      for(i=0;i<r;i++){
	for(j=0;j<c;j++)
	  v(j)=T(i,j);
	IFT(v,v_);
	for(j=0;j<c;j++)
	  W(i,j)=v_(j)*(1./c);	
      }
    }else{
      ComplexVector v(r),v_(r); 
      for(i=0;i<c;i++){
	for(j=0;j<r;j++)
	  v(j)=T(j,i);
	IFT(v,v_);
	for(j=0;j<c;j++)
	  W(j,i)=v_(j)*(1./r);	
      }
    }
  }
}
*/
/*
  ---------------------------------------------------
  
  LOWEST EIGENVECTOR OF A HERMITIAN MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

extern "C" void zheevx_(char* jobz,char* range,char* uplo,int* n,
			std::complex<double>* a,int* lda,double* vl,
			double* vu,int *il,int* iu,double* abstol,
			int* m,double* w,std::complex<double>* z,int* ldz,
			std::complex<double>* work,int* lwork,
			double* rwork,int* iwork,int* ifail,int* info);

extern "C" void dsyevx_(char* jobz,char* range,char* uplo,int* n,
			double* a,int* lda,double* vl,
			double* vu,int *il,int* iu,double* abstol,
			int* m,double* w,double* z,int* ldz,
			double* work,int* lwork,
			int* iwork,int* ifail,int* info);

namespace ula {
  
  void lowest(ComplexMatrix& H,Real &E,ComplexMatrix& V) {
    char jobz ='V';  
    char range= 'I';    
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    double vl;
    double vu;
    int il=1;
    int iu=1;
    double abstol=1e-9;
    int m=1;
    double* w=new double[n];
    int ldz=n;       
    int lda = n;                          // leading dimension
    int lwork = 2*n;  //before 2*n-1 changed for actualization 25.09.2009          
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[7*n];
    int* iwork=new int[5*n];
    int* ifail=new int[n];    
    int info;
    zheevx_( &jobz,&range, &uplo,&n, &(H(0,0)), &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work,&lwork, rwork, iwork, ifail,&info);
    E=w[0];
    delete[] work;
    delete[] rwork;
    delete[] w;
    delete[] iwork;
    delete[] ifail;
  }
  
 void first(ComplexMatrix& H,Real &E,ComplexMatrix& V) {
    char jobz ='V';  
    char range= 'I';    
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    double vl;
    double vu;
    int il=2;
    int iu=2;
    double abstol=1e-9;
    int m=1;
    double* w=new double[n];
    int ldz=n;       
    int lda = n;                          // leading dimension
    int lwork = 2*n;  //before 2*n-1 changed for actualization 25.09.2009          
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[7*n];
    int* iwork=new int[5*n];
    int* ifail=new int[n];    
    int info;
    zheevx_( &jobz,&range, &uplo,&n, &(H(0,0)), &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work,&lwork, rwork, iwork, ifail,&info);
    E=w[0];
    delete[] work;
    delete[] rwork;
    delete[] w;
    delete[] iwork;
    delete[] ifail;
  }

  
  void lowest(RealMatrix& H,Real &E,RealMatrix& V) {
    char jobz ='V'; 
    char range= 'I';   
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    double vl;
    double vu;
    int il=1;
    int iu=1;
    double abstol=1e-9;
    int m=1;
    double* w=new double[n];
    int ldz=n;       
    int lda = n;                          // leading dimension
    int lwork = -1;           
    double* work_= new double[1];
    int* iwork=new int[5*n];
    int* ifail=new int[n];   
    int info;
    dsyevx_( &jobz,&range, &uplo,&n, &(H(0,0)),
	     &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work_,&lwork, iwork, ifail,&info);
    lwork = (int)work_[0];           
    double* work = new double[lwork];
    dsyevx_( &jobz,&range, &uplo,&n, &(H(0,0)),
	     &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work,&lwork, iwork, ifail,&info);   
    E=w[0];
    delete[] work;
    delete[] work_;
    delete[] w;
    delete[] iwork;
    delete[] ifail;
  }

   void first(RealMatrix& H,Real &E,RealMatrix& V) {
    char jobz ='V'; 
    char range= 'I';   
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    double vl;
    double vu;
    int il=2;
    int iu=2;
    double abstol=1e-9;
    int m=1;
    double* w=new double[n];
    int ldz=n;       
    int lda = n;                          // leading dimension
    int lwork = -1;           
    double* work_= new double[1];
    int* iwork=new int[5*n];
    int* ifail=new int[n];   
    int info;
    dsyevx_( &jobz,&range, &uplo,&n, &(H(0,0)),
	     &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work_,&lwork, iwork, ifail,&info);
    lwork = (int)work_[0];           
    double* work = new double[lwork];
    dsyevx_( &jobz,&range, &uplo,&n, &(H(0,0)),
	     &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work,&lwork, iwork, ifail,&info);   
    E=w[0];
    delete[] work;
    delete[] work_;
    delete[] w;
    delete[] iwork;
    delete[] ifail;
  }
}

/*
  ---------------------------------------------------
  
  LARGEST EIGENVECTOR OF A HERMITIAN MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

     extern "C" void zheevx_(char* jobz,char* range,char* uplo,int* n,
			     std::complex<double>* a,int* lda,double* vl,
			     double* vu,int *il,int* iu,double* abstol,
			     int* m,double* w,std::complex<double>* z,int* ldz,
			     std::complex<double>* work,int* lwork,
			     double* rwork,int* iwork,int* ifail,int* info);

namespace ula {
  
  void largest(ComplexMatrix& H_,Real &E,ComplexMatrix& V) {
    char jobz ='V';  
    char range= 'I';    
    char uplo = 'U';                      // upper triangle is stored
    int n = H_.size1();                    // squareness is not checked    
    ComplexMatrix H(n,n);
    H=H_;
    double vl;
    double vu;
    int il=n;
    int iu=n;
    double abstol=1e-9;
    int m=1;
    double* w=new double[n];
    int ldz=n;       
    int lda = n;                          // leading dimension
    int lwork = 2*n-1;            
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[7*n];
    int* iwork=new int[5*n];
    int* ifail=new int[n];    
    int info;
    zheevx_( &jobz,&range, &uplo,&n, &(H(0,0)), &lda,&vl,&vu,&il,&iu,&abstol,
	     &m,w,&(V(0,0)),&ldz,work,&lwork, rwork, iwork, ifail,&info);
    E=w[0];
    delete[] work;
    delete[] rwork;
    delete[] w;
    delete[] iwork;
    delete[] ifail;
  }
}

/*
  ---------------------------------------------------
  
  DIAGONALIZATION OF A HERMITIAN MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  E: eigenvalues
  V: eigenvectors
  
  ---------------------------------------------------
*/

extern "C" void zheev_( char* jobz, char* uplo, 
                        int* n, std::complex<double>* a, int* lda, 
                        double* w,
                        std::complex<double>* work, int* lwork, double* rwork,
                        int* info);

namespace ula {
  
  void _cdiag(ComplexMatrix& H, RealVector& E, bool calc_ev ) {
    char jobz = ( calc_ev? 'V' : 'N' );   // calculate eigenvectors?
    char uplo = 'U';                      // upper triangle is stored
    int n = H.size1();                    // squareness is not checked
    int lda = n;                          // leading dimension
    int lwork = 3*n-1;                    // check optimal size in lwork after call
    int lrw = 3*n-2;                      // this size is hardwired
    std::complex<double>* work = new std::complex<double>[lwork];
    double* rwork = new double[lrw];
    int info;
    zheev_( &jobz, &uplo, 
	    &n, &(H(0,0)), &lda, 
	    &(E(0)), 
	    work, &lwork, rwork, 
	    &info);
    _message( info, lwork );
    delete[] work;
    delete[] rwork;
  }
  
  void diag(ComplexMatrix& H, RealVector& E) {
    _cdiag( H, E, false );
  }
  
  void diag(ComplexMatrix& H, RealVector& E, ComplexMatrix& V) {
    V=H;
    _cdiag( V, E, true );
  }
}



/*
  ---------------------------------------------------
  
  DETERMINANT OF A HERMITIAN MATRIX
  
  H: n-by-n input matrix, with upper triangle stored
  
  ----------------------------------------------------
*/

extern "C" void zgetrf_(int *m, int *n, 
                        std::complex<double> *a, int *lda, 
                        int* ipiv, 
                        int* info );

namespace ula {
  
  std::complex<double> determinant( ComplexMatrix& H ) {
    int n = H.size1();      // squareness is not checked!!!
    int lda = n;
    int* ipiv = new int[n];
    int info;
    zgetrf_( &n, &n, 
	     &(H(0,0)), &lda, 
	     ipiv, 
	     &info );
    int sign = 1;
    for( int i = 1; i <= n; i++ )
      if( ipiv[i-1] != i ) sign = -sign;
    std::complex<double> det = 1.0*sign;
    for( int i = 0; i < n; i++ ) det *= H(i,i*n);
    delete[] ipiv;
    return det;
  }
}
 


/*
  ---------------------------------------------------
  
    SORT ROUTINE FOR REAL EIGENVALUES
    
    ----------------------------------------------------
  */

namespace ula {
  
  template <typename valuetype>
    void eigensort(RealVector& d, boost::numeric::ublas::matrix<valuetype, boost::numeric::ublas::column_major, std::vector<valuetype > >& v) {
    int k,j,i;
    int n=v.size1();
    for (i=0;i<n-1;++i) {
      double p=d(k=i);
      for (j=i+1;j<n;++j)
	if (d(j)>=p)
	  p=d(k=j);
      if (k!=i) {
	d(k)=d(i);
	d(i)=p;
	for (j=0;j<n;++j) {
	  valuetype v_=v(j,i);
	  v(j,i)=v(j,k);
	  v(j,k)=v_;
	};
      }
    };
  }
}

#endif

