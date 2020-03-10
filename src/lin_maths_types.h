#ifndef lin_maths_types_h
#define lin_maths_types_h 

#include <vector>   //To interface with regular vectors
#include <iostream> //For cerr stream
#include <iomanip>  //For setw
#include <limits>   //For numeric limits and epsilon
#include <cmath>    //For floor

//For BLAS methods
#include "cblas.h"
#include "lapacke.h"

//For location based matricies #TODO:Specific matricies should be in their files
#include "location.h"

using namespace std;

enum class Dimension { row, column };

///////////////////////////////////////////////////////////////////////
// Matrix Types ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

// Mat interface /////////////////////////////////////////////////////
template <class T> class mat {
    protected:
        //The size of the matrix
        unsigned int dim_r, dim_c;

        const T tol_value = sqrt(std::numeric_limits<T>::epsilon());

    public:
        virtual ~mat() {} //So this can be virtualised out

        virtual T get(int, int) const = 0 ;
        virtual void set(int, int, const T) = 0;

        //Give a copy of the data as a whole array
        virtual T* get_as_full_array() const ;

        virtual void transpose() = 0;
        //virtual mat<T>& Tra();

        virtual mat<T>& operator=(const mat<T>&);

        virtual mat<T>& operator+=(const mat<T>&);
        virtual mat<T>& operator-=(const mat<T>&);
        virtual mat<T>& operator/=(const T);
        virtual mat<T>& operator*=(const T);

        virtual bool operator==(const mat<T>&) const;
        
        virtual mat<T>& eyes();
        virtual mat<T>& eyes(mat<T>&); 
        virtual mat<T>& eye_square(mat<T>&); //Unsafe method

         //Dimension getters and setters
        virtual int smallest_dim_size() const;
        virtual int get_dim(Dimension d) const { return (d==Dimension::row) ? this->dim_r : this->dim_c; };
        virtual int get_dim_r() const { return get_dim(Dimension::row); };
        virtual int get_dim_c() const { return get_dim(Dimension::column); };
        
        virtual void set(T); //Set all values to this one
        virtual void set(int, int, const mat<T>&); 
        virtual void set(const mat<T> &vals) { set(0, 0, vals); }; 

        virtual void print(ostream& where_to_output = std::cout, bool space=true) const;
        virtual void print_dims() const { printf("%i %i\n", get_dim(Dimension::row), get_dim(Dimension::column)); };

        virtual T get_tolerance() const { return tol_value; };

    protected:
        virtual void set_dim(Dimension d, int i) { if(d==Dimension::row) this->dim_r=i; else this->dim_c=i; };
        virtual void set_dim_r(int i) { return set_dim(Dimension::row, i); };
        virtual void set_dim_c(int i) { return set_dim(Dimension::column, i); };

};

// Matrix /////////////////////////////////////////////////////////////
template <class T> class matrix : public mat<T> {
    protected:
        //The values of the matrix on disk
        T *values = NULL;

    public:
        //Constructors and destructors
        matrix<T>(int, int);
        matrix<T>(int, int, const T vals[]); 
        matrix<T>(const mat<T>&);

        //The rule of three
        matrix<T>(const matrix<T>& M);
        virtual matrix<T>& operator=(const matrix<T>&); 
        virtual ~matrix<T>();

        virtual void transpose();
        
        //Low level access
        virtual void generate_data_holder();// { this->data = new Full_Numeric_Array<T>(this->dim_r, this->dim_c); };

        virtual matrix<T>& operator+=(const matrix<T>& sq) { this->mat<T>::operator+=(sq); return *this; }
        virtual matrix<T>& operator-=(const matrix<T>& sq) { this->mat<T>::operator-=(sq); return *this; }
        virtual matrix<T>& operator*=(const T dub) { this->mat<T>::operator*=(dub); return *this; }
        virtual matrix<T>& operator/=(const T dub) { this->mat<T>::operator/=(dub); return *this; }
        //virtual matrix<T>& operator*=(const mat<T>& sm);

        //Standard matrixhematrixical operations
        //virtual matrix<double> truncate_SVD();
        virtual void multiply_row(int i, double to_mult_by);
        virtual void bracket_diag(const mat<T>&);

        //Getters and setters 
        virtual T get(int, int) const;
        virtual matrix<T> get(int, int, int, int) const;

        virtual void set(T v) { mat<T>::set(v); } //pass through 
        virtual void set(int, int, const T);
        virtual void set(int i, int j, const mat<T>& m) { mat<T>::set(i, j, m); } //pass through
        virtual void set(const mat<T>& m) { mat<T>::set(m); } //pass through

        virtual T max() const;
        virtual T min() const;

        virtual T* get_all() const { return this->values; } //Give the underlying data structure here

        // For fiddling around with 0
        virtual bool check_zero() const;
        virtual matrix<T>& boil_off_smalls();

        //Printing methods
        virtual ostream& operator<<(ostream& stream) { 
            for(unsigned int i=0; i<this->dim_r; i++) {
                for(unsigned int j=0; j<this->dim_c; j++) {
                    stream << get(i,j);
                }
            }
            return stream;
        };
};

// Vector //////////////////////////////////////////////////////////////////////////////////////////
template <class T> class vec : public matrix<T> {
    public:
        vec<T>(vector<T> in);
        vec<T>(int r) : matrix<T>(r, 1) {};
        vec<T>(int i, const T vals[]) : matrix<T>(i, 1, vals) {};
        vec<T>(const matrix<T> &to_clone) : matrix<T>(to_clone) {};

        //The rule of three
        vec<T>(const vec<T> &to_clone) : matrix<T>(to_clone) {};
        virtual vec<T>& operator=(const vec<T>& rhs) { matrix<T>::operator=(rhs); return *this; };
        virtual ~vec<T>() {};

        //TODO: Transposing written pretty poorly, should be more like /peartor overloads?
        virtual void transpose();
        //virtual vec<T> Tra() const { vec<T> output = vec<T>(*this); output.transpose(); return output; }

        vec<T> piecewise_square(); 
        virtual vec<T> piecewise_root(); 
        virtual double average() const;
        virtual T sum() const;
        virtual T min() const;
        virtual T max() const;
    
        virtual vec<T>& piecewise_times_in_place(const vec<T>&); 
        virtual vec<T>& operator+=(const vec<T>&); 
        virtual vec<T>& operator-=(const vec<T>&); 
        virtual vec<T>& operator/=(const T);
        virtual vec<T>& operator*=(const T dub) { this->matrix<T>::operator*=(dub); return *this; }

        virtual T get(int) const;
        virtual void add(int, const vec<T>&);
        
        //virtual int get_dim(Dimension d) const { return matrix<T>::get_dim(d); } //This function shouldn't be needed? We should just inheret the matrix method?
        virtual int get_size() const; 
        virtual Dimension get_prominent_dim() const;

        virtual vec<T> get_sub_vec(int, int) const;
        virtual vector<T> as_vector();
        
        virtual T L2_norm() const;

        virtual void set(T val) { matrix<T>::set(val); }
        virtual void set(int i, T val) { matrix<T>::set(i, 0, val); }
        virtual void set(int, const vec<T>&); //This can be underwritten by the mat version
        virtual void set(const vec<T>& vals) { set(0, vals); }; 
};

// Diagonal //////////////////////////////////////////////////////////////////////////////////////////
template <class T> class diagonal : public mat<T> {
    protected:
        //The values of the matrix on disk
        T *values;

        unsigned int access_dim;

    public:
        diagonal<T>(int);
    
        //The rule of three
        diagonal<T>(const diagonal<T> &to_clone);
        virtual diagonal<T>& operator=(const diagonal<T>& rhs);
        virtual ~diagonal<T>();

        virtual void transpose() {};

        virtual diagonal<T>& operator+=(const diagonal<T>& d) { this->mat<T>::operator+=(d); return *this; } 
        virtual diagonal<T>& operator-=(const diagonal<T>& d) { this->mat<T>::operator-=(d); return *this; } 
        virtual diagonal<T>& operator/=(const T dub) { this->mat<T>::operator/=(dub); return *this; }
        virtual diagonal<T>& operator*=(const T dub) { this->mat<T>::operator*=(dub); return *this; }
        //virtual diagonal<T>& operator*=(const diagonal<T>& sm);
                     
        virtual void generate_data_holder();

        virtual T* get_data() const { return this->values; };

        virtual int smallest_dim_size() const { return this->get_dim(Dimension::row); };
        virtual int get_access_dim() const { return access_dim; };


        //TODO:I think these aren't being called corretly or someting is messing around with the diag_data, possibly deleting it
        virtual T get(int i) const { return get(i,i); }
        virtual T get(int i, int j) const; //{ return get_data()->get(i,j); } 
        virtual void set(int i, T val) { return set(i,i, val); }
        virtual void set(int i, int j, T val);// { return get_data()->set(i, j, val); }

        virtual diagonal<T> get_sub_diag(int start_i, int end_i) const { 
            diagonal<T> output = diagonal<T>(end_i-start_i);

            for(int i=start_i; i<end_i; i++)
                output.set(i-start_i, get(i));

            return output; 
        }

        T* get_as_full_array() const;
};

// Square //////////////////////////////////////////////////////////////////////////////////////////
template <class T> class LU_dual;

template <class T> class square : public matrix<T> {
    public:
        square<T>(int i) : matrix<T>(i,i) {}
        square<T>(const mat<T>& s) : matrix<T>(s) {} //TODO:This doesn't guarrente any squareness at all!

        //The rule of three
        square<T>(const square<T> &to_clone) : matrix<T>(to_clone) {};
        virtual square<T>& operator=(const square<T>& rhs) { matrix<T>::operator=(rhs); return *this; };
        virtual ~square<T>() {};

        virtual square<T>& operator+=(const square<T>& sq) { this->matrix<T>::operator+=(sq); return *this; }
        virtual square<T>& operator-=(const square<T>& sq) { this->matrix<T>::operator-=(sq); return *this; }
        virtual square<T>& operator*=(const T dub) { this->matrix<T>::operator*=(dub); return *this; }
        virtual square<T>& operator/=(const T dub) { this->matrix<T>::operator/=(dub); return *this; }
        //virtual square<T>& operator*=(const square<T>& sm) { this->matrix<T>::operator*=(sm); return *this; }

        virtual LU_dual<T> get_LU_decom(); 
        virtual square<T> grab_main_diagonoal_mat();
        virtual vec<T> eigen_values();
        virtual square<T> inverse();
        virtual double get_condition_number();

        virtual void set(int i, int j, const T v) { matrix<T>::set(i, j, v); } //pass through
        virtual void set(int i, int j, const mat<T>& m) { matrix<T>::set(i, j, m); } //pass through
        virtual void set(const mat<T>& m) { matrix<T>::set(m); } //pass through

        virtual int smallest_dim_size() const { return this->get_dim(Dimension::row); }
        virtual int get_size() const { return this->get_dim(Dimension::row); }
        //virtual int get_dim(Dimension d) const {return mat<T>::get_dim(d); } //TODO:Why did I have to make this, it should just come from mat!
};

template <class T> class LU_dual {
    public:
        square<T> lower;
        square<T> upper;

        LU_dual(int size) : lower(size), upper(size) {}
};

// Symmetric ///////////////////////////////////////////////////////////////////////////////////////
template <class T> class symmetric : public square<T> {
    public:
        symmetric<T>(int s) : square<T>(s) {};
        symmetric<T>(const mat<T> &s) : square<T>(s) {}; //TODO:Doesn't guarentee symmetry

        //The rule of three
        symmetric<T>(const symmetric<T> &to_clone) : square<T>(to_clone) {};
        virtual symmetric<T>& operator=(const symmetric<T>& rhs) { square<T>::operator=(rhs); return *this; };
        virtual ~symmetric<T>() {};


        //TODO: Write custom operators here since we don't need to apply everything to every value with a symetric matrix
        virtual symmetric<T>& operator+=(const symmetric<T>& sq) { this->square<T>::operator+=(sq); return *this; }
        virtual symmetric<T>& operator-=(const symmetric<T>& sq) { this->square<T>::operator-=(sq); return *this; }
        virtual symmetric<T>& operator*=(const T dub) { this->square<T>::operator*=(dub); return *this; }
        virtual symmetric<T>& operator/=(const T dub) { this->square<T>::operator/=(dub); return *this; }
        //virtual symmetric<T>& operator*=(const symmetric<T>& sm) { this->square<T>::operator*=(sm); return *this; }
        
        virtual void transpose() { }; //It's symmetric, so transpose does nothing
        //virtual symmetric<T> Tra() const { symmetric<T> output = symmetric<T>(*this); output.transpose(); return output; }

        virtual square<T> cholesky_square_root();
        virtual symmetric<T> eigen_square_root();
        virtual void eigen_square_root_inplace();
        virtual square<T> circulant_eigen_square_root();
        //virtual square<T> symmetric_square_root();
        
        virtual void set(int i, int j, const T v) { square<T>::set(i, j, v); } //pass through
        virtual void set(int i, int j, const mat<T>& m) { square<T>::set(i, j, m); } //pass through
        virtual void set(const mat<T>& m) { square<T>::set(m); } //pass through

        //TODO: This is allowed with a mat sometimes and only symmetric<T> sometimes and sometimes not allowed with either, what would be the best way to code this so that it doesn't break at compile time
        //Maybe give only one index for symetric ? Still doesn't solve the mat problem
        //virtual void set(int, int,  const mat&); 
};

// Block matrix /////////////////////////////////////////////////////////////
template <class T> class block : public mat<T> {
    protected:
        //The values of the matrix on disk
        vector<matrix<T>*> blocks; 

        //TODO:Only need int here, not double
        matrix<double> block_use;
    
        //The sub matrix dimensions
        unsigned int block_dim_r, block_dim_c, block_amount_r, block_amount_c;

    public:
        //Constructors and destructors
        block<T>(int, int, int, int);

        //The rule of three
        block<T>(const block<T>& M);
        virtual block<T>& operator=(const block<T>&); 
        virtual ~block<T>();

        virtual void transpose(); 

        virtual void generate_data_holder();
        //virtual vector<matrix<T>*> get_data() const { return this->blocks; }; 

        virtual int get_block_dim(Dimension d) const { return (d==Dimension::row) ? this->block_dim_r: this->block_dim_c; };
        virtual int get_block_amount(Dimension d) const { return (d==Dimension::row) ? this->block_amount_r: this->block_amount_c; };
        virtual int get_block_amount_held() const { return blocks.size(); };

        virtual matrix<T>& get_block(int i) const { return *(blocks[i]); };

        //void set_block(int block_i, int i, int j, const T val);
        virtual void set_block(int block_i_r, int block_j_c, const mat<T> &block);
        virtual T get(int, int) const;
        virtual void set(int, int, const T);

        virtual T* get_as_full_array() const;

        virtual int get_block_use(int i, int j) const { return block_use.get(i,j); };
        virtual matrix<double> get_full_block_use() const { return block_use; };
};

template <class T> class block_diag : public mat<T> {
    protected:
        vector<matrix<T>*> blocks; //Note that one would prefer to have a vector of mat<T>, that way each block could be any type of matrix. Howver, vector doesn't allow for such behavior (I think because it needs to know the size of the elements it holds)

        unsigned int block_dim_r, block_dim_c, block_amount;

    public:
        block_diag(int block_r, int block_c, int block_amount);

        //The rule of three
        block_diag<T>(const block_diag<T> &to_clone);
        virtual block_diag<T>& operator=(const block_diag<T>& rhs);
        virtual ~block_diag();

        //TODO: make opeator*(mat<T>), can be made much more efficient with all the blocks stuff

        virtual void transpose(); 

        virtual void generate_data_holder();
        virtual vector<matrix<T>*> get_data() const { return this->blocks; }; 

        virtual int get_block_dim(Dimension d) const { return (d==Dimension::row) ? this->block_dim_r: this->block_dim_c; };
        virtual matrix<T>& get_block(int i) const { return *(blocks[i]); };

        void set_block(int block_i, int i, int j, const T val);
        virtual void set_block(int block_i, const mat<T> &block) { 
            if((unsigned int)block_i<block_amount) 
                *blocks[block_i] = block; //TODO:Check we're not writing to the heap or seomthing
            else 
                printf("Error: block_diag, writing to incorrect block\n");
        }; 

        virtual T get(int, int) const;
        virtual void set(int, int, const T);

        virtual void set(T v) { mat<T>::set(v); }; //TODO:Can be made more efficient by overriding 
        virtual void set(int i, int j, const mat<T>& M) { mat<T>::set(i, j, M); }; //TODO:Can be made more efficient by overriding
        virtual void set(const mat<T> &vals) { this->set(0, 0, vals); }; 

        virtual T* get_as_full_array() const;

        int get_block_amount() const { return block_amount; }
};

// Kronker product matricies //////////////////////////////////////////////////////////////////////
template <class T> class kroneker: public mat<T> {
    protected:
        //The values of the matrix on disk
        matrix<T> left_matrix, right_matrix;

    public:
        kroneker(int left_r, int left_c, int right_r, int right_c);

        //The rule of three
        kroneker<T>(const kroneker<T> &to_clone);
        virtual kroneker<T>& operator=(const kroneker<T>& rhs);
        virtual ~kroneker<T>() {}; //Everything will delete itself

        virtual void transpose();

        virtual kroneker<T>& operator+=(const kroneker<T>& d) { this->mat<T>::operator+=(d); return *this; } 
        virtual kroneker<T>& operator-=(const kroneker<T>& d) { this->mat<T>::operator-=(d); return *this; } 
        virtual kroneker<T>& operator/=(const T dub) { this->mat<T>::operator/=(dub); return *this; }
        virtual kroneker<T>& operator*=(const T dub) { this->mat<T>::operator*=(dub); return *this; }
                     
        virtual void redo_dims();

        virtual T get(int i, int j) const; 
        virtual void set(int i, int j, T val);

        virtual void set_left_matrix(matrix<T> M) { left_matrix = M; redo_dims(); }
        virtual void set_right_matrix(matrix<T> M) { right_matrix = M; redo_dims(); }

        virtual matrix<T> get_left_matrix() const { return left_matrix; }
        virtual matrix<T> get_right_matrix() const { return right_matrix; }

        virtual matrix<T>& get_left_matrix_ref() { return left_matrix; }
        virtual matrix<T>& get_right_matrix_ref() { return right_matrix; }

};

//////////////////////////////////////////////////////////////////////////////////////////////////
// Matrix Operations /////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////

// Matrix operations ////////////////////////////////////////////////////////////////////////////// 
//TODO: Many of these can be simplified with gemm
template<class T> matrix<T> operator+(const mat<T> &lhs, const mat<T>& rhs) { matrix<T> output = matrix<T>(lhs); output += rhs; return output; }
template<class T> matrix<T> operator-(const mat<T> &lhs, const mat<T>& rhs) { matrix<T> output = matrix<T>(lhs); output -= rhs; return output; }
template<class T> matrix<T> operator/(const mat<T> &lhs, const T rhs) { matrix<T> output = matrix<T>(lhs); output /= rhs; return output; }
template<class T> matrix<T> operator*(const mat<T> &lhs, const T rhs) { matrix<T> output = matrix<T>(lhs); output *= rhs; return output; }
template<class T> matrix<T> operator*(const T lhs, const mat<T> &rhs) { matrix<T> output = matrix<T>(rhs); output *= lhs; return output; }
template<class T> matrix<T> operator*(const mat<T>& lhs, const mat<T>& rhs);

template<class T> matrix<T> operator*(const matrix<T>& lhs, const mat<T>& rhs);
template<class T> matrix<T> operator*(const mat<T>& lhs, const matrix<T>& rhs);
template<class T> matrix<T> operator*(const matrix<T>& lhs, const matrix<T>& rhs);

template<class T> matrix<T> Transpose(const matrix<T>& M) { matrix<T> output = matrix<T>(M); output.transpose(); return output; }

// Vector operations ////////////////////////////////////////////////////////////////////////////// 
template<class T> vec<T> piecewise_times(const vec<T>& lhs, const vec<T>& rhs);

//TODO: Simplify many of these with axpy
template<class T> vec<T> operator+(vec<T> lhs, const vec<T>& rhs) { lhs += rhs; return lhs;}
template<class T> vec<T> operator-(vec<T> lhs, const vec<T>& rhs) { lhs -= rhs; return lhs; }
template<class T> vec<T> operator*(vec<T> lhs, const T rhs) { lhs *= rhs; return lhs; }
template<class T> vec<T> Transpose(const vec<T>& M) { vec<T> output = vec<T>(M); output.transpose(); return output; }

template<class T> vec<T> mat_vec_mult(const mat<T>& lhs, const vec<T>& rhs);
template<class T> vec<T> operator*(const mat<T>& lhs, const vec<T>& rhs) { return mat_vec_mult(lhs, rhs); }
template<class T> vec<T> operator*(const matrix<T>& lhs, const vec<T>& rhs);

// Diagonal operations ////////////////////////////////////////////////////////////////////////////// 
template<class T> vec<T> operator*(const diagonal<T>& lhs, vec<T> rhs);
template<class T> diagonal<T> operator*(const diagonal<T>& lhs, const diagonal<T>& rhs);

template<class T> matrix<T> diag_mult(const mat<T>& lhs, const diagonal<T>& rhs);
template<class T> matrix<T> diag_mult(const diagonal<T>& rhs, const mat<T>& lhs);
template<class T> matrix<T> operator*(const mat<T>& lhs, const diagonal<T>& rhs); 
template<class T> matrix<T> operator*(const matrix<T>& lhs, const diagonal<T>& rhs);
template<class T> matrix<T> operator*(const diagonal<T>& lhs, const mat<T>& rhs);
template<class T> matrix<T> operator*(const diagonal<T>& lhs, const matrix<T>& rhs);

template<class T> diagonal<T> Transpose(const diagonal<T>& M) { diagonal<T> output = diagonal<T>(M); output.transpose(); return output; }

template<class T> T diag_square_over_space(const vec<T>&, const diagonal<T>&);

// Square operations ////////////////////////////////////////////////////////////////////////////// 
// In these instances the copies are done for you with the first grument
template<class T> square<T> operator+(square<T> lhs, const square<T>& rhs) { lhs += rhs; return lhs; }
template<class T> square<T> operator-(square<T> lhs, const square<T>& rhs) { lhs -= rhs; return lhs; }
template<class T> square<T> operator/(square<T> lhs, const T rhs) { lhs /= rhs; return lhs; }
template<class T> square<T> operator*(const square<T> &lhs, const T rhs) { square<T> output = square<T>(lhs); output *= rhs; return output; }
template<class T> square<T> operator*(const T lhs, const square<T> &rhs) { square<T> output = square<T>(rhs); output *= lhs; return output; }

template<class T> square<T> Transpose(const square<T>& M) { square<T> output = square<T>(M); output.transpose(); return output; }

// Symmetric operations ////////////////////////////////////////////////////////////////////////////// 
template<class T> symmetric<T> operator+(symmetric<T> lhs, const symmetric<T>& rhs) { lhs += rhs; return lhs; }
template<class T> symmetric<T> operator-(symmetric<T> lhs, const symmetric<T>& rhs) { lhs -= rhs; return lhs; }
template<class T> symmetric<T> operator*(symmetric<T> lhs, const T rhs) { lhs *= rhs; return lhs; }
template<class T> symmetric<T> operator/(symmetric<T> lhs, const T rhs) { lhs /= rhs; return lhs; }
template<class T> symmetric<T> Transpose(const symmetric<T>& M) { symmetric<T> output = symmetric<T>(M); output.transpose(); return output; }

// Block operations /////////////////////////////////////////////////////////////////////////
template<class T> vec<T> operator*(const block<T>& lhs, const vec<T>& rhs); 
template<class T> matrix<T> operator*(const block<T>& lhs, const matrix<T>& rhs); 
template<class T> block<T> operator*(const block<T>& lhs, const diagonal<T>& rhs);
template<class T> block<T> Transpose(const block<T>& M) { block<T> output = block<T>(M); output.transpose(); return output; }


// Block Diagonal operations /////////////////////////////////////////////////////////////////////////
template<class T> vec<T> operator*(const block_diag<T>& lhs, const vec<T>& rhs); 
template<class T> matrix<T> operator*(const block_diag<T>& lhs, const matrix<T>& rhs); 

template<class T> block_diag<T> operator*(const block_diag<T>& lhs, const diagonal<T>& rhs);
template<class T> block_diag<T> Transpose(const block_diag<T>& M) { block_diag<T> output = block_diag<T>(M); output.transpose(); return output; }

// Kroneker Operations ///////////////////////////////////////////////////////////////////////////////
template<class T> vec<T> operator*(const kroneker<T>& lhs, const vec<T>& rhs);

template<class T> vec<T> square_over_space_main_diag(const kroneker<T>& weight, const block_diag<T>& to_square);

template<class T> matrix<T> operator*(const block<T>& lhs, const kroneker<T>& rhs);
template<class T> matrix<T> operator*(const kroneker<T>& lhs, const block<T>& rhs );
template<class T> matrix<T> square_over_space(const kroneker<T>& weight, const block<T>& to_square);

template<class T> matrix<T> operator*(const kroneker<T>& lhs, const block_diag<T>& rhs);

#endif
