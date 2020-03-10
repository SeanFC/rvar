#ifndef lin_maths_h 
#define lin_maths_h 

#include <boost/math/special_functions/bessel.hpp> 
#include "lin_maths_types.h"

using base_type = float;

using Mat_d = matrix<double>;
using Vec_d = vec<double>;
using Symmetric_d = symmetric<double>;
using Square_d = square<double>;
using Diagonal_d = diagonal<double>;
using Block_d = block<double>;
using Block_Diag_d = block_diag<double>;
using Kroneker_d = kroneker<double>;

using Mat_f = matrix<float>;
using Vec_f = vec<float>;
using Symmetric_f = symmetric<float>;
using Square_f = square<float>;
using Diagonal_f = diagonal<float>;
using Block_f = block<float>;
using Block_Diag_f = block_diag<float>;
using Kroneker_f = kroneker<float>;

Vec_d lin_gemv(double alpha, const Mat_d& A, const Vec_d& x, double beta, const Vec_d& y);
Vec_f lin_gemv(double alpha, const Mat_f& A, const Vec_f& x, double beta, const Vec_f& y);

/////////////////////////////////////////////////////////////////////////////////
//Generic Matricies//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
Square_d get_bessel_correlation_matrix_d(int matrix_size, double length_scale);
Square_f get_bessel_correlation_matrix_f(int matrix_size, double length_scale);
Square_d get_bessel_correlaton_matrix_d(vector<loc> state_location_grid, double state_grid_scale);
Square_f get_bessel_correlaton_matrix_f(vector<loc> state_location_grid, double state_grid_scale);
//square<double> get_laplacian_correlation_matrix_inv(int matrix_size, double length_scale);

Square_d Ident_d(int size); 
Mat_d Ident_d(int r, int c);

Square_f Ident_f(int size);
Mat_f Ident_f(int r, int c); 

/////////////////////////////////////////////////////////////////////////////////
//Utility functions//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
int mod (int a, int b);
vector<double> temp_averages(vector<double>);
double sum(vector<double>);
double average(vector<double>);
vector<double> spread(double, int);

double scaled_bessel_correlation(double angle_between_points, double length_scale, double radius);

#endif
