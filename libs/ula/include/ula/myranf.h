#ifndef WESSEL_ULA_MYRANF
#define WESSEL_ULA_MYRANF

long _myranf_iseed_ = 42l;
long _myranf_state_  = _myranf_iseed_;

inline void set_myseed(long l) {
  _myranf_iseed_ = l;
  _myranf_state_  = _myranf_iseed_;
}

inline long get_myseed() {
  return _myranf_iseed_;
}

inline double myranf() {
  const long ia=16807l,ic=2147483647l,iq=127773l,ir=2836l;
  long  il,ih,it;
  double rc;
  ih = _myranf_state_/iq;
  il = _myranf_state_%iq;
  it = ia*il-ir*ih;
  if (it > 0)
    _myranf_state_ = it;
  else
   _myranf_state_ = ic+it;
  rc = ic;
  return _myranf_state_/rc;
}

inline  int myranf(const  int& a, const int& b) {
// myranf is in [a,b] and integer
  return  int(a+(b-a+1)*myranf());
}
 
inline long myranf(const long& a, const long& b) {
// myranf is in [a,b] and integer
  return long(a+(b-a+1)*myranf());
}
 
inline double myranf(const double& a, const double& b) {
// myranf is in (a,b)
  return a+(b-a)*myranf();
}

#endif
