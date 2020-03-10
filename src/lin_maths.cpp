#include "lin_maths.h"
#include "lin_maths_types.cpp"

template<>
matrix<double> operator*(const matrix<double>& lhs, const matrix<double>& rhs) {
    matrix<double> output = matrix<double>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    cblas_dgemm(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
            lhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M, 
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N, 
            rhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint K, 
            1,                //OPENBLAS_CONST double alpha,
            lhs.get_all(),    //OPENBLAS_CONST double *A
            lhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            rhs.get_all(),    //OPENBLAS_CONST double *B,
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint ldb,
            0,                //OPENBLAS_CONST double beta, 
            output.get_all(), //double *C, 
            rhs.get_dim(Dimension::column)    //OPENBLAS_CONST blasint ldc,
            );
    return output;
}

template<>
matrix<float> operator*(const matrix<float>& lhs, const matrix<float>& rhs) {
    matrix<float> output = matrix<float>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    cblas_sgemm(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
            lhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M, 
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N, 
            rhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint K, 
            1,                //OPENBLAS_CONST double alpha,
            lhs.get_all(),    //OPENBLAS_CONST float *A
            lhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            rhs.get_all(),    //OPENBLAS_CONST float *B,
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint ldb,
            0,                //OPENBLAS_CONST double beta, 
            output.get_all(), //float *C, 
            rhs.get_dim(Dimension::column)    //OPENBLAS_CONST blasint ldc,
            );
    return output;
}

template<> 
matrix<double> operator*(const mat<double>& lhs, const mat<double>& rhs) {
    matrix<double> output = matrix<double>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    double* lhs_data = lhs.get_as_full_array();
    double* rhs_data = rhs.get_as_full_array();

    //output.set(lhs.get(0,0)*rhs.get(0,0));
    cblas_dgemm(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
            lhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M, 
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N, 
            rhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint K, 
            1,                //OPENBLAS_CONST double alpha,
            lhs_data,    //OPENBLAS_CONST double *A
            lhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            rhs_data,    //OPENBLAS_CONST double *B,
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint ldb,
            0,                //OPENBLAS_CONST double beta, 
            output.get_all(), //double *C, 
            rhs.get_dim(Dimension::column)    //OPENBLAS_CONST blasint ldc,
            );
    delete [] lhs_data;
    delete [] rhs_data;

    return output;
}

template<>
matrix<float> operator*(const mat<float>& lhs, const mat<float>& rhs) {
    matrix<float> output = matrix<float>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    float* lhs_data = lhs.get_as_full_array();
    float* rhs_data = rhs.get_as_full_array();
    
    cblas_sgemm(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
            lhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M, 
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N, 
            rhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint K, 
            1,                //OPENBLAS_CONST double alpha,
            lhs_data,    //OPENBLAS_CONST float *A
            lhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            rhs_data,    //OPENBLAS_CONST float *B,
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint ldb,
            0,                //OPENBLAS_CONST double beta, 
            output.get_all(), //float *C, 
            rhs.get_dim(Dimension::column)    //OPENBLAS_CONST blasint ldc,
            );

    delete [] lhs_data;
    delete [] rhs_data;

    return output;
}

template<> 
matrix<float> operator*(const matrix<float>& lhs, const mat<float>& rhs) {
    matrix<float> output = matrix<float>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    float* rhs_data = rhs.get_as_full_array();
    cblas_sgemm(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
            lhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M, 
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N, 
            rhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint K, 
            1,                //OPENBLAS_CONST double alpha,
            lhs.get_all(),    //OPENBLAS_CONST float *A
            lhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            rhs_data,    //OPENBLAS_CONST float *B,
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint ldb,
            0,                //OPENBLAS_CONST double beta, 
            output.get_all(), //float *C, 
            rhs.get_dim(Dimension::column)    //OPENBLAS_CONST blasint ldc,
            );
    
    delete [] rhs_data;
    return output;
}

template<>
matrix<float> operator*(const mat<float>& lhs, const matrix<float>& rhs) {
    matrix<float> output = matrix<float>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    float* lhs_data = lhs.get_as_full_array();
    cblas_sgemm(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA, 
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, 
            lhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M, 
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N, 
            rhs.get_dim(Dimension::row),   //OPENBLAS_CONST blasint K, 
            1,                //OPENBLAS_CONST double alpha,
            lhs_data,    //OPENBLAS_CONST float *A
            lhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            rhs.get_all(),    //OPENBLAS_CONST float *B,
            rhs.get_dim(Dimension::column),   //OPENBLAS_CONST blasint ldb,
            0,                //OPENBLAS_CONST double beta, 
            output.get_all(), //float *C, 
            rhs.get_dim(Dimension::column)    //OPENBLAS_CONST blasint ldc,
            );
    
    delete [] lhs_data;
    return output;
}



template <>
double vec<double>::L2_norm() const {
    return cblas_dnrm2(get_size(), this->get_all(), 1);
}


template <>
float vec<float>::L2_norm() const {
    return cblas_snrm2(get_size(), this->get_all(), 1);
}

//TODO: This needs to be rewritten much like the inplace method
template <>
symmetric<float> symmetric<float>::eigen_square_root() {
    square<float> cur_copy = square<float>(this->get_size());
    cur_copy.set(*this);

    vec<float> evals = vec<float>(this->get_size());
    vec<float> evals_im = vec<float>(this->get_size());
    square<float> evecs = square<float>(this->get_size());

    int spectral_info = LAPACKE_sgeev( 
            LAPACK_ROW_MAJOR,//int matrix_layout, 
            'N',//char jobvl, Don't compute left eigenvectors, although it's just the right evecs transpose so might be worth it
            'V',//char jobvr, 
            this->get_size(),//lapack_int n, 
            cur_copy.get_all(),//float* a, 
            this->get_size(),//lapack_int lda, 
            evals.get_all(),//float* wr, 
            evals_im.get_all(),//float* wi, 
            NULL,//float* vl, 
            evecs.get_size(),//lapack_int ldvl, 
            evecs.get_all(),//float* vr, 
            evecs.get_size()//lapack_int ldvr 
            );
    
    if(spectral_info)
        throw "Spectral solver error code:"+ to_string(spectral_info) + "\n";

    //This step requires positive evals
    diagonal<float> D_h = diagonal<float>(evals.get_size());
    for(int i=0; i<D_h.get_dim(Dimension::row); i++) {
        if(evals.get(i) < 0)  {
            //if(abs(evals.get(i)) < 1e-3)
            if (evals.get(i) < this->tol_value*10) printf("%i Negative eigen %f\n", i, evals.get(i));
            D_h.set(i,i,0);
            //else
            //    throw "Negative eval " + to_string(evals.get(i)) + " at position " + to_string(i) + " was just too small to ignore and hence a real square root of the matrix given doesn't exist\n";
        } else {
            D_h.set(i,i,sqrt(evals.get(i)));//sqrt(abs(evals.get(i))));
        }
    }
    
    //If the evecs are orthogonal then we should only need transpose here instead of inverse but this isn't always the case for some reason?
    //You need evecs.inverse() in here or else all the variables get jumbled around. It should still be technically correct but confusing

    symmetric<float> output = symmetric(evecs.get_dim(Dimension::row));
    output.set(evecs * D_h * evecs.inverse());
    return output;
}

// Note: This original used LAPACKE_ssyev but this has a bug maybe? and can't create the correct eigenvectors 
template <>
void symmetric<float>::eigen_square_root_inplace() {
    vec<float> evals = vec<float>(this->get_size());

    char jobz = 'V';
    char uplo = 'U';
    int n = get_size();
    int LDA = n;
    int lwork = -1; 
    float wkopt;
    float *work;
    int spectral_info;
    
    // Find the size and create the work array
    LAPACK_ssyev(
            &jobz, &uplo,
            &n,
            this->get_all(), &LDA,
            evals.get_all(),
            &wkopt, &lwork,
            &spectral_info
            );
    lwork = (int)wkopt;
    work = new float[lwork]; //TODO: Delete work

    // Find the eigenvalues and eigenvectors
    LAPACK_ssyev(
            &jobz, &uplo,
            &n,
            this->get_all(), &LDA,
            evals.get_all(),
            work, &lwork,
            &spectral_info);
    delete [] work;
    // Here *this rows are populated with the eigenvectors

    if(spectral_info) {
        printf("ERROR: Spectral solver error code: %s\n", to_string(spectral_info).c_str());
        throw std::runtime_error("ERROR: Spectral solver");
    }

    //This step requires positive evals
    diagonal<float> D_h = diagonal<float>(evals.get_size());

    for(int i=0; i<D_h.get_dim_r(); i++) 
        D_h.set(i, i, sqrt(abs(evals.get(i))));
    
    square<float> trans = square<float>(get_size()); trans.set(*this); trans.transpose(); // We have to set up this transpose because *this is currently no symmetric so Transpose(*this) wont be corret

    // Set this matrix to the output
    this->set(trans*D_h*(*this));
}

//////////////////////////////////////////////////////////////////////////////////
// Instantiations ////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//Double instansiation
template matrix<double> operator*<double>(const mat<double>& lhs, const diagonal<double>& rhs); 
template matrix<double> operator*<double>(const matrix<double>& lhs, const diagonal<double>& rhs);
template matrix<double> operator*<double>(const diagonal<double>& lhs, const mat<double>& rhs);
template matrix<double> operator*<double>(const diagonal<double>& lhs, const matrix<double>& rhs);

template vec<double> operator*<double>(const kroneker<double>&, const vec<double>&);

template vec<double> operator*<double>(const block_diag<double>& lhs, const vec<double>& rhs); 
template block_diag<double> operator*<double>(const block_diag<double>& lhs, const diagonal<double>& rhs);
template double diag_square_over_space(const vec<double>&, const diagonal<double>&);

template matrix<double> operator*<double>(const mat<double>& lhs, const mat<double>& rhs);
template matrix<double> operator*<double>(const matrix<double>& lhs, const mat<double>& rhs);
template matrix<double> operator*<double>(const mat<double>& lhs, const matrix<double>& rhs);
template matrix<double> operator*<double>(const matrix<double>& lhs, const matrix<double>& rhs);

template vec<double> operator*(const matrix<double>& lhs, const vec<double>& rhs);
template vec<double> piecewise_times<double>(const vec<double>& lhs, const vec<double>& rhs);

template vec<double> operator*<double>(const diagonal<double>& lhs, vec<double> rhs);
template diagonal<double> operator*<double>(const diagonal<double> &lhs, const diagonal<double>& rhs);
template matrix<double> operator*(const block_diag<double>& lhs, const matrix<double>& rhs); 

template vec<double> square_over_space_main_diag<double>(const kroneker<double>& weight, const block_diag<double>& to_square);

//Float instansiation
template matrix<float> operator*<float>(const mat<float>& lhs, const diagonal<float>& rhs); 
template matrix<float> operator*<float>(const matrix<float>& lhs, const diagonal<float>& rhs);
template matrix<float> operator*<float>(const diagonal<float>& lhs, const mat<float>& rhs);
template matrix<float> operator*<float>(const diagonal<float>& lhs, const matrix<float>& rhs);

template vec<float> operator*<float>(const kroneker<float>&, const vec<float>&);

template vec<float> operator*<float>(const block<float>& lhs, const vec<float>& rhs); 
template block<float> operator*<float>(const block<float>& lhs, const diagonal<float>& rhs);
template matrix<float> operator*(const block<float>& lhs, const matrix<float>& rhs); 

template vec<float> operator*<float>(const block_diag<float>& lhs, const vec<float>& rhs); 
template block_diag<float> operator*<float>(const block_diag<float>& lhs, const diagonal<float>& rhs);
template matrix<float> operator*(const block_diag<float>& lhs, const matrix<float>& rhs); 

template float diag_square_over_space(const vec<float>&, const diagonal<float>&);

template matrix<float> operator*<float>(const mat<float>& lhs, const mat<float>& rhs);
template matrix<float> operator*<float>(const matrix<float>& lhs, const mat<float>& rhs);
template matrix<float> operator*<float>(const mat<float>& lhs, const matrix<float>& rhs);
template matrix<float> operator*<float>(const matrix<float>& lhs, const matrix<float>& rhs);

template vec<float> operator*(const matrix<float>& lhs, const vec<float>& rhs);
template vec<float> piecewise_times<float>(const vec<float>& lhs, const vec<float>& rhs);

template vec<float> operator*<float>(const diagonal<float>& lhs, vec<float> rhs);
template diagonal<float> operator*<float>(const diagonal<float>& lhs, const diagonal<float>& rhs);

template vec<float> square_over_space_main_diag<float>(const kroneker<float>& weight, const block_diag<float>& to_square);
template matrix<float> square_over_space<float>(const kroneker<float>& weight, const block<float>& to_square);
template matrix<float> operator*(const block<float>& lhs, const kroneker<float>& rhs);
template matrix<float> operator*(const kroneker<float>& lhs, const block<float>& rhs);
template matrix<float> operator*(const kroneker<float>& lhs, const block_diag<float>& rhs);

template class matrix<double>;
template class vec<double>;
template class symmetric<double>;
template class square<double>;
template class diagonal<double>;
template class block<double>;
template class block_diag<double>;
template class kroneker<double>;

template class matrix<float>;
template class vec<float>;
template class symmetric<float>;
template class square<float>;
template class diagonal<float>;
template class block<float>;
template class block_diag<float>;
template class kroneker<float>;

///////////////////////////////////////////////////////////////////////////////////////////////
//Utility functions////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
int mod (int a, int b) {
   if(b < 0) //you can check for b == 0 separately and do what you want
     return mod(a, -b);   
   int ret = a % b;
   if(ret < 0)
     ret+=b;
   return ret;
}

//Find the sum of all the numbers in a vector
double sum(vector<double> v) {
    double running_total = 0;

    int size = v.size();

    for(int i=0; i<size; i++) {
        running_total += v[i];
    }

    return running_total;
}

//Find the average of all the numbers in a vector
double average(const vector<double> &v) {
    double running_total = 0;

    int size = v.size();

    for(int i=0; i<size; i++) {
        running_total += v[i];
    }

    return running_total/size; 
}

//Take the total and spread it over the amount of segements given by over_amount
vector<double> spread(double total, int over_amount) {
    vector<double> spread_segments = vector<double>();
    double var_to_set = total/over_amount;

    for(int i=0; i<over_amount; i++) {
        spread_segments.push_back(var_to_set);
    }

    return spread_segments;
}

double scaled_bessel_correlation(double angle_between_points, double length_scale, double radius) {
    double k = radius/length_scale;
    double scaled_angle = abs(sin(angle_between_points/2));

    if(scaled_angle == 0)
        return 1.0;

    double bessel = boost::math::cyl_bessel_k(1, k * scaled_angle);
    return k * scaled_angle * bessel;
}

///////////////////////////////////////////////////////////////////////////////////
//Generic Matricies////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//TODO:These matricies are circulant, could save a lot of memory by specifying this?
square<double> get_bessel_correlation_matrix_d(int matrix_size, double length_scale) {
    square<double> bessel_correlation = square<double>(matrix_size);
    double radius = 1; //TODO:This is sort of just set through the length scale 
    
    for(int i=0; i<matrix_size; i++) {
        for(int j=0; j<matrix_size; j++) {
            //TODO: This doesn't account for different days the each month
            double angle_between_points = abs(i-j);
            if(angle_between_points > matrix_size/2.0)
                angle_between_points = matrix_size - angle_between_points;
            angle_between_points *= (2*M_PI)/matrix_size;

            bessel_correlation.set(i, j, scaled_bessel_correlation(angle_between_points, length_scale, radius));
        }

    }
    return bessel_correlation;
}

square<float> get_bessel_correlation_matrix_f(int matrix_size, double length_scale) {
    square<float> bessel_correlation = square<float>(matrix_size);
    double radius = 1; //TODO:This is sort of just set through the length scale 
    
    for(int i=0; i<matrix_size; i++) {
        for(int j=0; j<matrix_size; j++) {
            double angle_between_points = abs(i-j);
            if(angle_between_points > matrix_size/2.0)
                angle_between_points = matrix_size - angle_between_points;
            angle_between_points *= (2*M_PI)/matrix_size;

            bessel_correlation.set(i, j, scaled_bessel_correlation(angle_between_points, length_scale, radius));
        }

    }
    return bessel_correlation;
}

square<double> get_bessel_correlaton_matrix_d(vector<loc> location_vector, double length_scale) {
    double radius = location_vector[0].get_earth_radius_km(); //Earth radius in KM
    int matrix_size = location_vector.size();

    square<double> bessel_correlation = square<double>(matrix_size); //This can be symmetric

    for(int i=0; i<matrix_size; i++) {
        for(int j=0; j<i; j++) {
            double scale = scaled_bessel_correlation(location_vector[i].angle_distance_between(location_vector[j]), length_scale, radius);
            bessel_correlation.set(i, j, scale);
            bessel_correlation.set(j, i, scale);
        }
        bessel_correlation.set(i, i, 1.0);
    }

    return bessel_correlation;
}

square<float> get_bessel_correlaton_matrix_f(vector<loc> location_vector, double length_scale) {
    double radius = location_vector[0].get_earth_radius_km(); //Earth radius in KM
    int matrix_size = location_vector.size();

    square<float> bessel_correlation = square<float>(matrix_size); //This can be symmetric

    for(int i=0; i<matrix_size; i++) {
        for(int j=0; j<i; j++) {
            double scale = scaled_bessel_correlation(location_vector[i].angle_distance_between(location_vector[j]), length_scale, radius);
            bessel_correlation.set(i, j, scale);
            bessel_correlation.set(j, i, scale);
        }
        bessel_correlation.set(i, i, 1.0);
    }

    return bessel_correlation;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//Extra functions///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
Vec_d lin_gemv(double alpha, const Mat_d& A, const Vec_d& x, double beta, const Vec_d& y)  {
    Vec_d output = y;

    cblas_dgemv(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order,
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
            A.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M,
            A.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N,
            alpha,
            A.get_all(),
            A.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            x.get_all(),
            1,
            beta,
            output.get_all(),
            1 
            );

    return output;
}

Vec_f lin_gemv(double alpha, const Mat_f& A, const Vec_f& x, double beta, const Vec_f& y)  {
    Vec_f output = y;

    cblas_sgemv(
            CblasRowMajor,    //OPENBLAS_CONST enum CBLAS_ORDER Order,
            CblasNoTrans,     //OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
            A.get_dim(Dimension::row),   //OPENBLAS_CONST blasint M,
            A.get_dim(Dimension::column),   //OPENBLAS_CONST blasint N,
            alpha,
            A.get_all(),
            A.get_dim(Dimension::column),   //OPENBLAS_CONST blasint lda, //How big is the leading dimension (so how long is a row, the amount of columns)
            x.get_all(),
            1,
            beta,
            output.get_all(),
            1 
            );

    return output;
}

Square_d Ident_d(int size) { return Square_d(size).eyes(); }
Mat_d Ident_d(int r, int c) { return Mat_d(r,c).eyes(); }

Square_f Ident_f(int size) { return Square_f(size).eyes(); }
Mat_f Ident_f(int r, int c) { return Mat_f(r,c).eyes(); }
