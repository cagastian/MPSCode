#include <iostream>
#include <ula/matrixtypes.h>
#include <ula/lapack.h>

int main () 
{

 ula::RealMatrix m(4,4);
 ula::ComplexVector e(4);
 ula::ComplexMatrix v(4,4);

 m(0,0)= 0.35; m(0,1)= 0.45; m(0,2)=-0.14; m(0,3)=-0.17;
 m(1,0)= 0.09; m(1,1)= 0.07; m(1,2)=-0.54; m(1,3)= 0.35;
 m(2,0)=-0.44; m(2,1)=-0.33; m(2,2)=-0.03; m(2,3)= 0.17;
 m(3,0)= 0.25; m(3,1)=-0.32; m(3,2)=-0.13; m(3,3)= 0.11;

 ula::ns_diag(m,e,v);

 std::cout << m << std::endl;
 std::cout << e << std::endl;
 std::cout << v << std::endl;

 for (int j=0;j<4;++j) {
   ula::ComplexVector t(4);
   for (int i=0;i<4;++i)
     t(i)=v(i,j);
   ula::ComplexVector res(4);
   res=prod(m,t)-e(j)*t;
   std::cout << res << std::endl;
 }

 return 0;
 
}
