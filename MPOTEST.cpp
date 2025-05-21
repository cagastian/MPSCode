#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <iomanip>
#include <iostream>
#include <lapack.h>
#include <matrixtypes.h>


#include<iostream>
#include<fstream>

using std::cout;
using std::cin;
using std::endl;
using std::flush;
using std::streambuf;
using std::ofstream;

using namespace ula;

#define PI 3.14159265358979323846

int m=16,M=20,L=16, w=2;

#include <MPSLib.h>
//#include <MPOlib.h>

//############################
//#### Previous Functions ####
//############################

void A_mat(ComplexMatrix& A){
	for(int S=0;S<L;S++){
		if(S==1) A(5*M,S*M)=1.;
		else if(S==3) A(10*M,S*M)=1.; 
		else if(S==2) A(15*M,S*M)=1.;
		else if(S==4) A(8*M,S*M)=1.;
		else if(S==5) A(13*M,S*M)=1.;
		else if(S==6) A(7*M,S*M)=1.;
		//if(S==5){ A(3*M,S*M)=0.5; A(5*M,S*M)=0.5;} 
		else A(0,S*M)=1.; //No particles on other sites 
	}
	Correct(A,0); 
	Correct(A,1);
}

void A_short(ComplexMatrix& A){
	
	for(int S=0;S<L;S++){
		if( S==1){ A(3*M,S*M)=1.; }
		else if( S==2){ A(7*M,S*M)=1.; }
		else if( S==3){ A(2*M,S*M)=1.; }
		else{A(0*M,S*M)=1.;}
	}
	Correct(A,0);
	Correct(A,1);
}

void A_short2(ComplexMatrix& A){
	
	for(int mm=0;mm<m;mm++){
		for(int S=0;S<L;S++){
			if(S==mm){ 
			if(mm==0){A(1*M,S*M)=1.;}
			else{A(mm*M,S*M)=1.;}
			  }
		}
	}
	Correct(A,0);
	Correct(A,1);
}

void A_short3(ComplexMatrix& A){
	
		for(int S=0;S<L;S++){
			A(14*M,S*M)=1.;
			  }
	
	Correct(A,0);
	Correct(A,1);
}

void ACalc(ComplexMatrix& A){
  int S;
  for(S=0;S<L;S++)
    //    A(0,S*M)=1.;
    { 
      A(0,S*M)  =sqrt(1.-2.*Gauss(S,10,1.5));
      A(M,S*M)  =sqrt(2.*Gauss(S,10,1.5))*exp(I*(double)S*PI/2.);
      A(2*M,S*M)=sqrt(2.*Gauss(S,10,1.5))*exp(-I*(double)S*3.*PI/2.);
    }
    /*  if(S==20)
    A(2*M,S*M)=1.;
  else
  A(0,S*M)=1.; */
  Correct(A,0);
  Correct(A,1); 
}

void Atest(ComplexMatrix& A){
	int rows=A.size1(),cols=A.size2();
	for(int i=0;i<rows;i++){
		for(int j=0;j<cols;j++){
			A(i,j)=1.;
		}
	}
}




//#########################################################################################################################################################################################################################
//#### Testing MPO Functions ##############################################################################################################################################################################################
//#########################################################################################################################################################################################################################

ComplexMatrix dot_prod(ComplexMatrix& A,ComplexMatrix& B){
    int A_r=A.size1() , A_c=A.size2(), B_c=B.size2();  
    // A_columns == B_rows  
    Complex temp;
    //https://stackoverflow.com/questions/9936132/why-does-the-order-of-the-loops-affect-performance-when-iterating-over-a-2d-arra
    //https://stackoverflow.com/questions/20467117/for-matrix-operation-why-is-ikj-faster-than-ijk

    ComplexMatrix C(A_r,B_c);
            
    for(int i=0;i<A_r;i++){
        for(int k=0;k<A_c;k++){
            temp = A(i,k);
            for(int j=0;j<B_c;j++){
                C(i,j) += temp*B(k,j);
            }        
        }
    
    }
    return C;
}

void CheckingNonZero(ComplexMatrix A){
	Real Win=0.;
	int Ar=A.size1(), Ac=A.size2();
	
	for (int Sr=0;Sr<Ar;Sr++){
		  	for (int Sc=0;Sc<Ac;Sc++){
		  		Win = real(A(Sr,Sc));
		  		if(Win != 0){
		  			cout<<"A("<<Sr<<","<<Sc<<")="<<Win<<endl;
		  		}
		  	}
		  }
}

void CheckingNonZero(ComplexVector A){
	Real Win=0.;
	int Ar=A.size();
	
	for (int Sr=0;Sr<Ar;Sr++){
		  		Win = real(A(Sr));
		  		if(Win != 0){
		  			cout<<"A("<<Sr<<")="<<Win<<endl;
		  		}
	}
}



//-------------------------Operators------------------------//

//Number operator MPO
ComplexMatrix Number_MPO(){
	int w=2;
	//m is the internal degrees of freedom
	ComplexMatrix N(w*m,w*m);
	for(int alpha=0;alpha<w;alpha++){
		for(int beta=0;beta<w;beta++){
			SubComplexMatrix Op(N,Range(alpha*m, (alpha+1)*m),Range(beta*m, (beta+1)*m));
			if(alpha==beta){
				Op = ident_op();
				}
			else if(alpha==1 && beta==0){
				Op = Number_Op2();
				}
			else{
				continue;
				}
		}
	}
	return N;
}

//----------------- Functions MPO------------------//

ComplexMatrix MPO_correct(ComplexMatrix &H,int k){
	//k is either 0 (first site) or L (last site)
	//H(w*m,w*m)
	int w=H.size1()/m;
	ComplexMatrix H_end(w*m,w*m);
	ComplexMatrix w_auxiliar(w*m,w*m); // Might be necessary to initialize with zeros.
	
	if(k==0){
		SubComplexMatrix identity(w_auxiliar,Range((w-1)*m, w*m),Range((w-1)*m, w*m));
		identity = ident_op();
		H_end = prod(w_auxiliar,H);
		
	}
	else if(k==L-1){
		
		SubComplexMatrix identity(w_auxiliar,Range(0, m),Range(0, m));
		identity = ident_op();
		H_end = prod(H,w_auxiliar);

	}
	else{
	cout<<"be aware, k is not a lattice end."<<endl;
	H_end=H;
	}
	
	//TESTING
	//cout<<"MPO for #"<<k<<" site"<<endl;
	//cout<<H<<endl;
	
	return H_end;
}


ComplexMatrix flatten_MPO(ComplexMatrix &H,int k){
	//H is the MatrixOp, dimensions (m*w,m*w)
	int w = H.size1()/m;
	ComplexMatrix H_mid(m*w,m*w);
	ComplexMatrix H_flat(m,w*w*m);
	
	H_mid=H;
	
	if(k==0){
		H_mid=MPO_correct(H,k);
	}
	else if(k==L-1){
		H_mid=MPO_correct(H,k);
	}
	else{
		H_mid=H;//Optimizable
	}
	
	
	
	//alpha= beta_s ; beta = beta_{s+1}
	for(int alpha=0;alpha<w;alpha++){
		for(int beta=0;beta<w;beta++){
			SubComplexMatrix Op(H_mid,Range(alpha*m, (alpha+1)*m),Range(beta*m, (beta+1)*m));
			SubComplexMatrix Op_flat(H_flat,Range(0,m),Range((alpha*w+beta)*m, (alpha*w+beta+1)*m));
			Op_flat = Op;
		
		}
	}
	
			//TESTING
			//cout<<"MPO Flatten for #"<<k<<" site"<<endl;
			//cout<<H_flat<<endl;
	
	return H_flat;
}

//Function to compute single site MPO product  with tensor
ComplexMatrix Single_site_opten_MPO(ComplexMatrix &O,ComplexMatrix &A, int k){
		//A(m*M,L*M) matrix orthonormalized	O operator (m,m)	k site
	int A_r = A.size1();
	SubComplexMatrix Ak(A,Range(0, A_r),Range(k*M, (k+1)*M));
	ComplexMatrix Ak_(A_r, M);//
	Ak_=Ak;
	//TESTING
	//cout<<"Ak; k="<<k<<endl;
	//cout<<Ak_<<endl;
	
	Ak_=opten(O,Ak_);//HERE'S THE ERROR 
	
	//TESTING
	//cout<<"Operator"<<endl;
	//cout<<O<<endl;
	//cout<<"After Opten k="<<k<<endl;
	//cout<<Ak_<<endl<<endl<<endl;
	
	return Ak_;
}

ComplexMatrix opten_MPO(ComplexMatrix &H_flat, ComplexMatrix& A, int k){//No flaw
	//H_flat is the MatrixOp Flattened to (m, m w^2)
	//We have in mind that A.size = (mM,M)
	ComplexMatrix ON(m*M*w*w,M);
	ComplexMatrix OA(m*M,M), Op(m,m);
	
	ComplexMatrix Op_squared(m,m); 
	
	//TESTING
	//cout<<"H flat (inside Opten MPO)"<<endl;
	//cout<<H_flat<<endl<<endl;
	
	for(int gamma=0;gamma<w*w;gamma++){
		SubComplexMatrix OA(ON,Range(gamma*m*M,(gamma+1)*m*M),Range(0,M));
		SubComplexMatrix Op(H_flat,Range(0,m),Range(gamma*m,(gamma+1)*m));
		
		Op_squared=Op;

			
		
		OA = Single_site_opten_MPO(Op_squared,A,k);
		
		if(k==5){
			//TESTING (Good)
			//cout<<"Op #"<<gamma<<endl;
			//CheckingNonZero(Op_squared);	
			
			
			}

	}
		//TESTING  (Good)
		//cout<<"ON k="<<k<<endl;
		//CheckingNonZero(ON);
	return ON;
}

//Trace over left and right (d) MPO transfer matrix
ComplexMatrix Trmu_MPO(ComplexMatrix &supE,bool d){
  
	ComplexMatrix a_MPO(w*M,w*M);
	ComplexMatrix a(M,M);
	ComplexMatrix temp_supE(M*M,M*M);
	for(int gamma=0;gamma<w;gamma++){
	   for(int gamma_=0;gamma_<w;gamma_++){
	   	SubComplexMatrix loc_supE(supE,Range(gamma_*M*M,(gamma_+1)*M*M),Range(gamma*M*M,(gamma+1)*M*M));
	   	temp_supE=loc_supE;
	   	
	   	SubComplexMatrix a(a_MPO,Range(gamma_*M,(gamma_+1)*M),Range(gamma*M,(gamma+1)*M));
		a=Trmu(temp_supE,d);
				}
		}
  return a_MPO;
}

//Transfer matrix with the MPO implementation
ComplexMatrix TrnsfrMtrx_MPO(ComplexMatrix &H, ComplexMatrix& A, int k){
	int alpha,beta,gamma,Ns,alpha_,beta_,gamma_;
	int w = H.size1()/m; //MPO bond dimension
		  	
  	SubComplexMatrix Ak_(A,Range(0,m*M),Range(k*M,(k+1)*M));
  	//First flatten the MPO in rows
  	ComplexMatrix H_flat(m,w*w*m);
  	H_flat=flatten_MPO(H,k);
  	
		  	//TESTING
		  	//cout<<"MPO Flatten for #"<<k<<" site"<<endl;
		  	//cout<<H_flat<<endl;
  	
  	//Apply the MPO over the state at k-th site 
  	ComplexMatrix ON(m*M*w*w,M);
  	ON=opten_MPO(H_flat,A,k);
			
			//TESTING
		 //cout<<"Checking for ON (k="<<k<<")"<<endl;
		 //CheckingNonZero(ON);		
  	
	
	ComplexMatrix OA(m*M,M);
	ComplexMatrix Ak(m,M*M);
  	ComplexMatrix Bk(m,M*M);
  	
  	ComplexMatrix E_(M*M,M*M);
  	ComplexMatrix E(M*M,M*M);
  	
  	ComplexMatrix loc_regE(M*M,M*M);
  	ComplexMatrix loc_supE_(M*M,M*M);
  	ComplexMatrix loc_supE(M*M,M*M);
  	
  	ComplexMatrix supE_(w*w*M*M,M*M);
  	ComplexMatrix regE(w*M*M,M*M);
  	ComplexMatrix supE(w*M*M,w*M*M);	
	
	for(gamma=0;gamma<w*w;gamma++){
		SubComplexMatrix OA(ON,Range(gamma*m*M,(gamma+1)*m*M),Range(0,M));
		SubComplexMatrix loc_supE_(supE_,Range(gamma*M*M,(gamma+1)*M*M),Range(0,M*M));
	  for(Ns=0;Ns<m;Ns++)
	    for(alpha=0;alpha<M;alpha++)
	      for(beta=0;beta<M;beta++){
		Bk(Ns,alpha*M+beta)=OA(Ns*M+alpha,beta);
		Ak(Ns,alpha*M+beta)=Ak_(Ns*M+alpha,beta);
	      }
	  noalias(E_)=prod(herm(Ak),Bk);
	  
	  for(alpha=0;alpha<M;alpha++){
	    for(beta=0;beta<M;beta++){
	      for(alpha_=0;alpha_<M;alpha_++){
		for(beta_=0;beta_<M;beta_++){
		  E(alpha*M+alpha_,beta*M+beta_)=E_(alpha*M+beta,alpha_*M+beta_);
		  			}
		 	 	}
		 	 }
		  }

	  loc_supE_=E;
	  
	}
	
	for(gamma=0;gamma<w;gamma++){
	   SubComplexMatrix regE(supE_,Range(gamma*w*M*M,(gamma+1)*w*M*M),Range(0,M*M));
	   ComplexMatrix temp_regE = regE;
	   for(gamma_=0;gamma_<w;gamma_++){
		SubComplexMatrix loc_regE(temp_regE,Range(gamma_*M*M,(gamma_+1)*M*M),Range(0,M*M));
		SubComplexMatrix loc_supE(supE,Range(gamma_*M*M,(gamma_+1)*M*M),Range(gamma*M*M,(gamma+1)*M*M));
		loc_supE=loc_regE;
		
		}
  	}
  	
	return supE;
}


//original (forward)
Real MPO_measurement_forward(ComplexMatrix& A, ComplexMatrix& H, int L) {
    // Compute the expectation value of an MPO using transfer matrices
    Real measure;
    //int Ar = A.size1(), Ac = A.size2();
    
    // Copy MPS tensors
    ComplexMatrix Abra = A;
    ComplexMatrix Aket = A;

    // Initialize vectors and matrices-
    ComplexVector E(w * M * M), Ef(w * M * M);
    ComplexMatrix supEE_1(w * M, w * M), supEE_2(w * M * M, w * M * M), EE_3(w * M * M, w * M * M);
    ComplexMatrix supE(w * M * M, w);

    // Compute the first transfer matrix incorporating the MPO
    supEE_2 = TrnsfrMtrx_MPO(H, Aket, 0);
    supEE_1 = Trmu_MPO(supEE_2, 0);	  

	 
    // Extract relevant elements into vector E
    for (int gamma = 0; gamma < w; gamma++) {
        for (int gamma_ = 0; gamma_ < w; gamma_++) {
            SubComplexMatrix EE_1(supEE_1, Range(gamma * M, (gamma + 1) * M), Range(gamma_ * M, (gamma_ + 1) * M));
            
            for (int alpha = 0; alpha < M; alpha++) {
                for (int alpha_ = 0; alpha_ < M; alpha_++) {
                    int idx = gamma * M * M + alpha * M + alpha_;
                    E(idx) = EE_1(alpha, alpha_);
                }
            }
        }
    }
	
    // Iterate through lattice sites to contract the transfer matrices
    //CheckingNonZero(Aket);
    for (int ss = 1; ss < L-1; ss++) {
    		
            EE_3 = TrnsfrMtrx_MPO(H, Aket, ss);
            E = prod(trans(EE_3),E);//Qué pasa si no transpongo 
            	
            	//TESTING
		cout<<"Checking for TrnsfrMtrx_ k="<<ss<<endl;
		CheckingNonZero(EE_3);	      
          
    }

    // Final site operation
    	
	supEE_2 = TrnsfrMtrx_MPO(H, Aket, L-1);
	supEE_1 = Trmu_MPO(supEE_2, 1);
	
	cout<<"Checking for TrnsfrMtrx_ k=L"<<endl;
	CheckingNonZero(supEE_2);	 
	
	for (int gamma = 0; gamma < w; gamma++) {
		for (int alpha = 0; alpha < M; alpha++) {
			for (int alpha_ = 0; alpha_ < M; alpha_++) {
				int idx = gamma * M * M + alpha * M + alpha_;
				Ef(idx) = supEE_1(gamma * M + alpha, gamma * M + alpha_);
                    }
                }
            }
        
	
    // Compute final measurement value as real inner product
   
    measure = real(inner_prod(E, Ef));

    return measure;
}

//backwards
Real MPO_measurement_backwards(ComplexMatrix& A, ComplexMatrix& H, int L) { //THIS ONE DOES IT FROM RIGHT TO LEFT.
    // Compute the expectation value of an MPO using transfer matrices
    Real measure;
    
    // Copy MPS tensors
    ComplexMatrix Abra = A;
    ComplexMatrix Aket = A;

    // Initialize vectors and matrices-
    ComplexVector E(w * M * M), Ef(w * M * M);
    ComplexMatrix supEE_1(w * M, w * M), supEE_2(w * M * M, w * M * M), EE_3(w * M * M, w * M * M);
    ComplexMatrix supE(w * M * M, w);

    // Compute the LAST transfer matrix incorporating the MPO
    supEE_2 = TrnsfrMtrx_MPO(H, Aket, L-1);
    supEE_1 = Trmu_MPO(supEE_2, 1);

		  //TESTING
		   cout<<"Checking for TrnsfrMtrx_ k=L"<<endl;
		   CheckingNonZero(supEE_2);	
		   
    // Extract relevant elements into vector E
    for (int gamma = 0; gamma < w; gamma++) {
        for (int gamma_ = 0; gamma_ < w; gamma_++) {
            SubComplexMatrix EE_1(supEE_1, Range(gamma * M, (gamma + 1) * M), Range(gamma_ * M, (gamma_ + 1) * M));
            
            for (int alpha = 0; alpha < M; alpha++) {
                for (int alpha_ = 0; alpha_ < M; alpha_++) {
                    int idx = gamma * M * M + alpha * M + alpha_;
                    E(idx) = EE_1(alpha, alpha_);
                }
            }
        }
    }
	
    // Iterate through lattice sites to contract the transfer matrices
    //CheckingNonZero(Aket);
    for (int ss = L-2; ss > 0; ss--) {
    		
            EE_3 = TrnsfrMtrx_MPO(H, Aket, ss);
            E = prod(E,EE_3);
            	
            	//TESTING
		cout<<"Checking for TrnsfrMtrx_ k="<<ss<<endl;
		CheckingNonZero(E);	      
          
    }

    // First site operation
    	
	supEE_2 = TrnsfrMtrx_MPO(H, Aket, 0);
	supEE_1 = Trmu_MPO(supEE_2, 0);
	
	cout<<"Checking for TrnsfrMtrx_ k=0"<<endl;
	CheckingNonZero(supEE_2);	 
	
	for (int gamma = 0; gamma < w; gamma++) {
		for (int alpha = 0; alpha < M; alpha++) {
			for (int alpha_ = 0; alpha_ < M; alpha_++) {
				int idx = gamma * M * M + alpha * M + alpha_;
				Ef(idx) = supEE_1(gamma * M + alpha, gamma * M + alpha_);
                    }
                }
            }
         
    // Compute final measurement value as real inner product
   
    measure = real(inner_prod(Ef,E));

    return measure;
}


//MPO_measurement_slow
Real MPO_measurement(ComplexMatrix& A, ComplexMatrix& H, int L) { 
    // Compute the expectation value of an MPO using transfer matrices
    Real measure;
    
    // Copy MPS tensors
    ComplexMatrix Abra = A;
    ComplexMatrix Aket = A;

    // Initialize vectors and matrices-
    ComplexVector E(w * M * M), Ef(w * M * M);
    ComplexMatrix supEE_1(w * M, w * M), supEE_2(w * M * M, w * M * M), EE_3(w * M * M, w * M * M);
    ComplexMatrix supE(w * M * M, w);

    // Compute the LAST transfer matrix incorporating the MPO
    
        cout<<"----------------------------"<<endl;
	cout<<"PRODUCT OF TRANSFER MATRICES (Last to First)"<<endl;
	cout<<"----------------------------"<<endl;
    
    
	ComplexMatrix Total_E(w * M * M, w * M * M), Temporal_E(w * M * M, w * M * M) ;
	for(int omega=L-1;omega>0;omega--){
		if(omega==L-1){
		Total_E = TrnsfrMtrx_MPO(H, Aket, omega);
		
		}else{
		Temporal_E = TrnsfrMtrx_MPO(H, Aket, omega);
			cout<<"Checking for TrnfrMTX k="<<omega<<endl;
		        CheckingNonZero(Temporal_E);	
		
		Total_E = prod(Total_E,Temporal_E);
		
			
		}
		//TESTING
			cout<<"Checking for Total_E k="<<omega<<endl;
		    	CheckingNonZero(Total_E);
		    	cout<<endl;
	}	 
	cout<<"----------------------------"<<endl;
	cout<<"            END             "<<endl;
	cout<<"----------------------------"<<endl;	  
	//Let's calculate the traces
	supEE_1 = Trmu_MPO(Total_E, 0);
	supEE_2 = Trmu_MPO(Total_E, 1);	   
		   
		   
    // Extract relevant elements into vector E
    for (int gamma = 0; gamma < w; gamma++) {
        for (int gamma_ = 0; gamma_ < w; gamma_++) {
            SubComplexMatrix EE_1(supEE_1, Range(gamma * M, (gamma + 1) * M), Range(gamma_ * M, (gamma_ + 1) * M));
            
            for (int alpha = 0; alpha < M; alpha++) {
                for (int alpha_ = 0; alpha_ < M; alpha_++) {
                    int idx = gamma * M * M + alpha * M + alpha_;
                    E(idx) = EE_1(alpha, alpha_);
                }
            }
        }
    }
	 
	
	for (int gamma = 0; gamma < w; gamma++) {
		for (int alpha = 0; alpha < M; alpha++) {
			for (int alpha_ = 0; alpha_ < M; alpha_++) {
				int idx = gamma * M * M + alpha * M + alpha_;
				Ef(idx) = supEE_2(gamma * M + alpha, gamma * M + alpha_);
                    }
                }
            }
         
         //This one is currently working 
         cout<<"--------------------------------------------------------"<<endl;
        cout<<"Vectors after traces"<<endl;
        cout<<"--------------------------------------------------------"<<endl;
        
        cout<<"Trace Right"<<endl;
        CheckingNonZero(Ef);
        
        cout<<"Trace left"<<endl;
        CheckingNonZero(E);
         
         
        cout<<"--------------------------------------------------------"<<endl;
        cout<<"Checking for the product of Transfer Matrices "<<endl;
        cout<<"--------------------------------------------------------"<<endl;
        
       	
    // Compute final measurement value as real inner product
   
    measure = real(inner_prod(Ef,E));

    return measure;
}








//#######################################################################################################################################################################################
//####         Main       ###############################################################################################################################################################
//#######################################################################################################################################################################################

int main(){

ofstream out("kokito.txt");
streambuf* oldCout = cout.rdbuf(); 
cout.rdbuf(out.rdbuf()); 


FILE *dskw;
char archive[300];
snprintf(archive, sizeof(archive), "output.txt");
dskw=fopen(archive,"w+");

ComplexMatrix A(M*m,M*L);

int A_r = A.size1();
ComplexMatrix Ak__(A_r, M);
ComplexMatrix Ak_(A_r, M);

//RACalc(A);
//ACalc(A);
A_short2(A);


//int k= 2; //Center of orthogonality
ComplexMatrix Op(m,m);

Op = Number_MPO();

Real result;
result=MPO_measurement(A,Op,L);
cout.rdbuf(oldCout); 


fflush(dskw);
cout<<result<<endl;

return 0;
}



