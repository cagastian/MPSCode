#include <lapack.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <matrixtypes.h>

int main () 
{

 ula::ComplexMatrix m(3,3);
 ula::ComplexMatrix v(3,3);
 ula::RealVector e(3);

 for (unsigned int i = 0; i < m.size1 (); ++ i) 
  for (unsigned int j = 0; j < m.size2 (); ++ j) 
    m (i, j) = i + j;

 std::cout << m+m << std::endl;

 ula::diag(m,e,v);
 ula::eigensort(e,v);

 std::cout << m << std::endl;
 std::cout << e << std::endl;
 std::cout << v << std::endl;

 return 0;
 
}
