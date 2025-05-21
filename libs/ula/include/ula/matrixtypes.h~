#ifndef WESSEL_MATRIXTYPES_
#define WESSEL_MATRIXTYPES_

#include <complex>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
 

namespace ula {
  
  typedef double Real;
  typedef std::complex<double> Complex;
  typedef boost::numeric::ublas::vector<double, std::vector<double> > RealVector;
  typedef boost::numeric::ublas::vector<std::complex<double>, std::vector<std::complex<double> > > ComplexVector;
    
  typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major, std::vector<double > >  RealMatrix;
  typedef boost::numeric::ublas::matrix<std::complex<double>, boost::numeric::ublas::column_major, std::vector<std::complex<double> > >  ComplexMatrix;
  typedef boost::numeric::ublas::matrix<double> SparseRealMatrix;
  typedef boost::numeric::ublas::matrix<std::complex<double> > SparseComplexMatrix;
  
  typedef  boost::numeric::ublas::vector_range<ComplexVector> SubComplexVector;
 
  typedef  boost::numeric::ublas::matrix_range<ComplexMatrix> SubComplexMatrix;
  typedef  boost::numeric::ublas::matrix_range<RealMatrix> SubRealMatrix;
  typedef  boost::numeric::ublas::matrix_row<ComplexMatrix> ComplexMatrixRow;
  typedef  boost::numeric::ublas::matrix_column<ComplexMatrix> ComplexMatrixColumn; 
  typedef  boost::numeric::ublas::matrix_row<RealMatrix> RealMatrixRow;
  typedef  boost::numeric::ublas::matrix_column<RealMatrix> RealMatrixColumn; 
  
}

#endif

