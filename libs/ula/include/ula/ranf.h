#ifndef WESSEL_ULA_RANF
#define WESSEL_ULA_RANF

#include <cmath>
#include <boost/config.hpp>
#include <boost/random/lagged_fibonacci.hpp>

using namespace std;

typedef boost::lagged_fibonacci607 Gen;

Gen mygen;

inline void set_seed(long s)
{
 mygen.seed(static_cast<boost::uint32_t>(s));
}

inline double ranf()
{
return mygen();
}

inline  int ranf(int a,int b)
{ // ranf is in [a,b] and integer
return  int(a+(b-a+1)*ranf());
}

inline long ranf(long a, long b)
{ // ranf is in [a,b] and long integer
return long(a+(b-a+1)*ranf());
}

inline double ranf(double a, double b)
{ // ranf is in (a,b)
return a+(b-a)*ranf();
}

#endif
