//make test_LAPACK1
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
//hacen parte de ula
#include <lapack.h>
#include <matrixtypes.h> //tipos de datos

using namespace std;
using namespace ula;

int main () 
{

 ComplexMatrix m(3,3);
 ComplexMatrix v(3,3);
 RealVector e(3);

 for (unsigned int i = 0; i < m.size1 (); ++ i) 
  for (unsigned int j = 0; j < m.size2 (); ++ j) 
    m (i, j) = i + j;

 cout << m+m << endl;

 diag(m,e,v);
 eigensort(e,v);

 cout << m << endl;
 cout << e << endl;
 cout << v << endl;

 return 0;
 
}
