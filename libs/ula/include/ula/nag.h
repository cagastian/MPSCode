#ifndef WESSEL_ULA_NAG
#define WESSEL_ULA_NAG

//equation solver
extern "C"  void c05nbf_(void (*) (int&, double*, double*, int&) ,
                    int*, double*, double*, double* , double*, int*, int*);

void c05nbf(void (*fnc) (int&, double*, double*, int&),
                    int  n, double* x, double* fvec, double tol , double* wa, int lwa, int& ifail) { 
  c05nbf_(fnc,&n,x,fvec,&tol,wa,&lwa,&ifail); 
}


//simplex minimalization
extern "C" void e04ccf_(int*, double*, double*, double*, int*,
                        double*, double*, double*, double*, double*, double*, 
                        void(*) (int&, double*, double&),
                        void(*) (double&, double&, double*, int&, int&, int&),
                        int*, int*);             
                        
void e04ccf(int n, double* x, double& f, double tol, int iw, 
            double *w1, double* w2, double* w3, double* w4, double* w5, double* w6, 
            void(*fnc) (int&, double*, double&), 
            void(*monit) (double&, double&, double*, int&, int&, int&), 
            int maxcall, int& ifail) {
  e04ccf_(&n,x,&f,&tol,&iw,w1,w2,w3,w4,w5,w6,fnc,monit,&maxcall,&ifail);
}      

#endif
