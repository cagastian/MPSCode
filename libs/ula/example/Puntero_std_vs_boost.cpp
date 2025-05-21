#include <iostream>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
//hacen parte de ula
#include <lapack.h>
#include <matrixtypes.h> //tipos de datos

using namespace std;
using namespace ula;

int main(){

  double* work;
  int N =5; 
  //hacer un arreglo a partir de uno puntero usamos el comando new

  work = new double [N];

  for(int i=0; i<N; i++){
    work[i]=i;
    cout<<work[i]<<endl;
  }

  //borrar un arreglo dinámico
  delete[] work; //la memora previamente ocupada se encuentra disponible

  //cout<<work[3]<<endl;//información NO CONFIABLE

  //-------a lo boost------------

  RealVector arreglo(4);
  double* p;
  p=&(arreglo(0)); //apunta al elemento arreglo(0)
  //p=arreglo; // No funciona para boost

  //------a lo std-----------

  double vector[4];
  double* pe;
  pe=vector; //apunta al elemento vector[0]
  
  return 0;;
}
