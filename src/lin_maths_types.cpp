#include "lin_maths_types.h"

//Note that putting the implemenation in the header file means that it will take longer to compile everything since each individual unit will need to compile mat<whatever> however if we don't do this the template has to be compiled for every type of class ever.
//The solution is to let each complilation unit compile what they want or to compile specific examples such as double and int (which we've done here)
//Read more here http://stackoverflow.com/questions/8752837/undefined-reference-to-template-class-constructor#8752879

//////////////////////////////////////////////////////////////////////////////////////////
// Generic mat interface /////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//Set the vector to the values of another given vector
template <class T>
void mat<T>::set(int start_position_r, int start_position_c, const mat<T> &setting_from) {
    //Get the size of this mat<T> and the one that we're setting from
    int amount_to_set_r = setting_from.get_dim(Dimension::row);
    int amount_to_set_c = setting_from.get_dim(Dimension::column);
    
    //Check that we're not setting more values than this matrix has space for, if we are just set as many as we have space for
    if(start_position_r+amount_to_set_r > get_dim_r())
        amount_to_set_r = get_dim_r() - start_position_r; 
    if(start_position_c+amount_to_set_c > get_dim_c())
        amount_to_set_c = get_dim_c() - start_position_c; 

    //Set the values of this matrix to the ones of the other matrix 
    for(int i=0; i<amount_to_set_r; i++)
        for(int j=0; j<amount_to_set_c; j++)
            set(start_position_r+i, start_position_c+j, setting_from.get(i,j));
}

template <class T>
T* mat<T>::get_as_full_array() const { 
    T* output = new T[this->dim_r*this->dim_c]();

    for(unsigned int i=0; i<this->dim_r; i++) 
        for(unsigned int j=0; j<this->dim_c; j++) 
            output[i*this->dim_c + j] = get(i,j);

    return output; 
}

template <class T>
mat<T>& mat<T>::eye_square(mat<T>& v) {
    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            set(i,j, get(i,j) * v.get(i, 0) * v.get(j, 0));
        }
    }

    return *this;
}

template <class T>
mat<T>& mat<T>::operator=(const mat<T>& rhs) {
    if(this != &rhs) {
        this->dim_r = rhs.get_dim(Dimension::row);
        this->dim_c = rhs.get_dim(Dimension::column);

        for(unsigned int i=0; i< this->dim_r; i++) {
            for(unsigned int j=0; j< this->dim_c; j++) {
                set(i,j, rhs.get(i,j));
            }
        }
    }

    return *this;
}


template <class T>
mat<T>& mat<T>::operator+=(const mat<T>& rhs) {
    for(int i=0; i<get_dim_r(); i++) {
        for(int j=0; j<get_dim_c(); j++) {
            set(i,j, get(i,j) + rhs.get(i,j));
        }
    }

    return *this;
}

template <class T>
mat<T>& mat<T>::operator-=(const mat<T>& rhs) {
    for(int i=0; i<get_dim_r(); i++) {
        for(int j=0; j<get_dim_c(); j++) {
            set(i,j, get(i,j) - rhs.get(i,j));
        }
    }

    return *this;
}

template <class T>
bool mat<T>::operator==(const mat<T>& o) const {
    if (o.get_dim(Dimension::row) != get_dim(Dimension::row))  return false;
    if (o.get_dim(Dimension::column) != get_dim(Dimension::column))  return false;

    for(int i=0; i<get_dim_r(); i++) 
        for(int j=0; j<get_dim_c(); j++) 
            if (!(get(i,j) == o.get(i,j))) return false;

    return false;
}

template <class T>
void mat<T>::set(T val) { 
    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            set(i,j, val);
        }
    }
}

//Find the dimension that has the smallest size and return said size
template <class T>
int mat<T>::smallest_dim_size() const {
    if(get_dim_r() > get_dim_c())
        return get_dim_c();

    return get_dim_r();
}

////////////////////////////////////////////////////////////////////////////////////////
// Specific matricies //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////

// Matrix ///////////////////////////////////////////////////
template <class T> 
matrix<T>::matrix(int r, int c) {
    this->dim_r = r< 1 ? 1 : r;
    this->dim_c = c< 1 ? 1 : c;

    generate_data_holder();
}

//TODO:Shouldn't this be a soft copy rather and just copying over the values?
template <class T>
matrix<T>::matrix(int r, int c, const T vals[]) : matrix(r,c) {
    for(int i=0; i<r; i++) {
        for(int j=0; j<c; j++) {
            set(i, j, vals[i*c + j]);
        }
    }
}

template <class T>
matrix<T>::matrix(const mat<T>& to_clone) {
    this->dim_r = to_clone.get_dim_r();
    this->dim_c = to_clone.get_dim_c();

    generate_data_holder();

    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            set(i,j, to_clone.get(i,j));
        }
    }
}

//We can't call the mat<T> constructor but we do the exact same thing
template <class T>
matrix<T>::matrix(const matrix<T>& to_clone) {
    this->dim_r = to_clone.get_dim_r();
    this->dim_c = to_clone.get_dim_c();

    generate_data_holder();

    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            set(i,j, to_clone.get(i,j));
        }
    }
}

//TODO:Need to reassign memory
template <class T>
matrix<T>& matrix<T>::operator=(const matrix<T>& rhs) {
    this->dim_r = rhs.get_dim_r();
    this->dim_c = rhs.get_dim_c();

    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            set(i,j, rhs.get(i,j));
        }
    }

    return *this;
}

template <class T>
matrix<T>::~matrix() {
    if(this->values) delete [] this->values;
}

//Make the thing that holds all the values
template <class T> 
void matrix<T>::generate_data_holder() {
    try {
        this->values = new T[this->dim_r*this->dim_c]();
    } catch (std::bad_alloc& ba) {
        std::cerr << "bad_alloc caught when getting memory for matrix: " << ba.what() << '\n';
        std::cerr << "sizes " << this->dim_r << " and " << this->dim_c << " for " << this->dim_r*this->dim_c << '\n';
    }
}

template <class T> 
T matrix<T>::get(int i, int j) const {
    if((unsigned int)i< this->dim_r || (unsigned int)j < this->dim_c)
        return this->values[(unsigned int)i*this->dim_c + (unsigned int)j];

    printf("Getting outside of range at %i %i for dimensions %i %i\n", i, j, this->dim_r, this->dim_c);
    return 0;
}


template <class T> 
void matrix<T>::set(int i, int j, T val) {
    if((unsigned int)i>= this->dim_r || (unsigned int)j >= this->dim_c) {
        printf("Setting outside of bounds\n");
        return;
    }

    this->values[(unsigned int)i*this->dim_c + (unsigned int)j] = val;
}

template <class T>
matrix<T> matrix<T>::get(int row_min, int col_min, int row_max, int col_max) const {
    matrix<T> output = matrix<T>(row_max - row_min, col_max - col_min);

    for(int i=0; i<output.get_dim(Dimension::row); i++) 
        for(int j=0; j<output.get_dim(Dimension::column); j++) 
            output.set(i,j, get(i+row_min, j+col_min));

    return output;
}

// Diagonal ///////////////////////////////////////////////////
template <class T> 
diagonal<T>::diagonal(int d) {
    this->dim_r = d<1 ? 1: d;
    this->dim_c = d<1 ? 1: d;

    this->access_dim = d<1 ? 1 : d;
    
    generate_data_holder();
}

template <class T> 
diagonal<T>::diagonal(const diagonal<T> &to_clone) {
    this->dim_r = to_clone.get_dim_r();
    this->dim_c = to_clone.get_dim_c();

    this->access_dim = to_clone.get_access_dim();
    
    generate_data_holder();

    for(unsigned int i=0; i<this->access_dim; i++)
        set(i, to_clone.get(i));

}

template <class T> 
diagonal<T>& diagonal<T>::operator=(const diagonal<T>& rhs) {
    if((unsigned int)rhs.get_dim_r() != this->dim_r) {
        printf("Error: Wrong sizes in diagonal copy assignment operator\n");
    }

    this->dim_r = rhs.get_dim_r();
    this->dim_c = rhs.get_dim_c();

    this->access_dim = rhs.get_access_dim();
    
    for(unsigned int i=0; i<this->access_dim; i++)
        set(i, rhs.get(i));
    
    return *this;
}

template <class T>
diagonal<T>::~diagonal() {
    if(this->values) delete [] this->values;
}

template <class T> 
void diagonal<T>::generate_data_holder() {
    try {
        this->values = new T[access_dim]();
    } catch (std::bad_alloc& ba) {
        std::cerr << "bad_alloc caught when getting memory for matrix: " << ba.what() << '\n';
    }

}

template <class T> 
T diagonal<T>::get(int r, int c) const {
    if(r==c)
        return this->values[(unsigned int)r];
    else
        return 0;
}

template <class T> 
void diagonal<T>::set(int r, int c, const T val) {
    if(r==c)
        this->values[(unsigned int)r] = val;
}

//Really slow functions, have to make a whole new array and populate it
//TODO:Implement this using Full_Numeric_Array
template <class T> 
T* diagonal<T>::get_as_full_array() const {
    printf("Diagonal: Realising full matrix\n");
    T *output = new T[access_dim*access_dim]();

    for(unsigned int i=0; i<access_dim; i++) {
        for(unsigned int j=0; j<access_dim; j++) {
            if(i==j)
                output[i*access_dim + j] = get(i,j);
            else
                output[i*access_dim + j] = 0;
        }
    }

    return output;
}

////////////////////////////////////////////////////////////////////////
// Block ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class T> 
block<T>::block(int block_dim_r, int block_dim_c, int block_amount_r, int block_amount_c) : block_use(block_amount_r, block_amount_c) {
    block_use.set(-1);
    this->block_dim_r = block_dim_r > 1 ? block_dim_r : 1;
    this->block_dim_c = block_dim_c > 1 ? block_dim_c : 1;
    this->block_amount_r = block_amount_r > 1 ? block_amount_r : 1;
    this->block_amount_c = block_amount_c > 1 ? block_amount_c : 1;
    this->dim_r = this->block_dim_r*block_amount_r;
    this->dim_c = this->block_dim_c*block_amount_c;

    generate_data_holder(); 
}

template <class T> 
block<T>::block(const block<T> &to_clone) : block_use(to_clone.get_block_amount(Dimension::row), to_clone.get_block_amount(Dimension::column)) {
    this->dim_r = to_clone.get_dim_r();
    this->dim_c = to_clone.get_dim_c();

    this->block_dim_r = to_clone.get_block_dim(Dimension::row);
    this->block_dim_c = to_clone.get_block_dim(Dimension::column);
    this->block_amount_r = to_clone.get_block_amount(Dimension::row);
    this->block_amount_c = to_clone.get_block_amount(Dimension::column);

    generate_data_holder();

    //Do copy assignment on every block    
    for(int b=0; b<to_clone.get_block_amount_held(); b++)
        blocks.push_back(new matrix<T>(to_clone.get_block(b)));

    block_use = to_clone.get_full_block_use();
}

//TODO:This needs to check that rhs is the correct size since this isn't a constructor
template <class T> 
block<T>& block<T>::operator=(const block<T>& rhs) {
    if((unsigned int)rhs.get_dim_r() != this->dim_r || (unsigned int)rhs.get_dim_c() != this->dim_c) {
        printf("Error: Wrong sizes in block copy operator\n");
    }

    this->dim_r = rhs.get_dim_r();
    this->dim_c = rhs.get_dim_c();

    this->block_dim_c = rhs.get_block_dim(Dimension::column);
    this->block_dim_r = rhs.get_block_dim(Dimension::row);
    this->block_amount_r = rhs.get_block_amount(Dimension::row);
    this->block_amount_c = rhs.get_block_amount(Dimension::column);

    generate_data_holder();

    //Do copy assignment on every block    
    for(int b=0; b<rhs.get_block_amount_held(); b++)
        blocks.push_back(new matrix<T>(rhs.get_block(b)));

    block_use = rhs.get_full_block_use();

    return *this;
}

template <class T> 
block<T>::~block() {
    for( auto &it : blocks )
        delete it;

    blocks.clear(); //Not sure this really does anything useful, salts the earth more than anything
}

template <class T> 
void block<T>::generate_data_holder() {
    //It's difficult to make a vector from a class that doesn't have a default constructor, instead we just make an empty vector and fill it as needs be
    blocks = vector<matrix<T>*>();
}

template <class T> 
T block<T>::get(int i, int j) const {
    //TODO:Error checking
    int block_num = block_use.get((int)(i/this->block_dim_r), (int)(j/this->block_dim_c));

    if(block_num >= 0)
        return blocks[block_num]->get(i%this->block_dim_r, j%this->block_dim_c);
    else
        return 0;
}

template <class T> 
void block<T>::set(int i, int j, const T val) {
    //TODO:Error checking
    int block_num = block_use.get((int)(i/this->block_dim_r), (int)(j/this->block_dim_c));

    if(block_num >= 0) {
        blocks[block_num]->set(i%this->block_dim_r, j%this->block_dim_c, val);
    } else {
        blocks.push_back(new matrix<T>(this->block_dim_r, this->block_dim_c));
        block_use.set((int)(i/this->block_dim_r), (int)(j/this->block_dim_c), blocks.size() - 1);
        blocks[blocks.size() - 1]->set(i%this->block_dim_r, j%this->block_dim_c, val);
    }
}

template <class T> 
void block<T>::transpose() { 
    block_use.transpose();

    for(auto &b : blocks) 
        b->transpose(); 

    int tmp = this->dim_c; 
    this->dim_c = this->dim_r; 
    this->dim_r = tmp; 

    tmp = this->block_dim_c; 
    this->block_dim_c = this->block_dim_r; 
    this->block_dim_r = tmp; 

    tmp = this->block_amount_c; 
    this->block_amount_c = this->block_amount_r; 
    this->block_amount_r = tmp; 
}

template <class T> 
void block<T>::set_block(int block_r, int block_c, const mat<T> &block) { 
    //TODO:Error checking
    int block_num = block_use.get(block_r, block_c);

    if(block_num >= 0) {
        (*blocks[block_num]) = block;
    } else {
        blocks.push_back(new matrix<T>(this->block_dim_r, this->block_dim_c));
        block_use.set(block_r, block_c, blocks.size() - 1);
        (*blocks[blocks.size() - 1]) = block;
    }
}

template <class T> 
T* block<T>::get_as_full_array() const {
    printf("Block: Realising full matrix\n");
    T *output = new T[this->block_dim_r*this->block_dim_c*block_amount_c*block_amount_r]();

    //TODO:Ugly for loops, is there some way to go through the blocks to do this?
    for(unsigned int b_r=0; b_r<this->block_amount_r; b_r++) {
        for(unsigned int b_c=0; b_c<this->block_amount_c; b_c++) {
            int block_num = block_use.get(b_r, b_c);

            if(block_num >=0) {
                for(unsigned int i=0; i<this->block_dim_r; i++) {
                    for(unsigned int j=0; j<this->block_dim_c; j++) {
                        //     block row                           block column            row in block    col in block
                        output[b_r*this->block_dim_r*this->dim_c + b_c*this->block_dim_c + i*this->dim_c + j           ] = blocks[block_num]->get(i,j);
                    }
                }
            }
        }
    }

    return output;
}


// Block Diagonal //////////////////////////////////////////////////////
template <class T> 
block_diag<T>::block_diag(int block_dim_r, int block_dim_c, int block_amount) {
    this->block_dim_r = block_dim_r > 1 ? block_dim_r : 1;
    this->block_dim_c = block_dim_c > 1 ? block_dim_c : 1;
    this->block_amount = block_amount > 1 ? block_amount : 1;
    this->dim_r = this->block_dim_r*block_amount;
    this->dim_c = this->block_dim_c*block_amount;

    generate_data_holder(); 
}

template <class T> 
block_diag<T>::block_diag(const block_diag<T> &to_clone) {
    this->dim_r = to_clone.get_dim_r();
    this->dim_c = to_clone.get_dim_c();

    this->block_dim_c = to_clone.get_block_dim(Dimension::column);
    this->block_dim_r = to_clone.get_block_dim(Dimension::row);
    this->block_amount= to_clone.get_block_amount();

    generate_data_holder();
    
    //Do copy assignment on every block    
    for(int b=0; b<this->get_block_amount(); b++)
        *blocks[b] = to_clone.get_block(b);
}

template <class T> 
block_diag<T>& block_diag<T>::operator=(const block_diag<T>& rhs) {
    if((unsigned int)rhs.get_dim_r() != this->dim_r || (unsigned int)rhs.get_dim_c() != this->dim_c) {
        printf("Error: Wrong sizes in block copy operator\n");
    }

    this->dim_r = rhs.get_dim_r();
    this->dim_c = rhs.get_dim_c();

    this->block_dim_c = rhs.get_block_dim(Dimension::column);
    this->block_dim_r = rhs.get_block_dim(Dimension::row);
    this->block_amount= rhs.get_block_amount();

    //Do copy assignment on every block    
    for(int b=0; b<this->get_block_amount(); b++)
        *blocks[b] = rhs.get_block(b);

    return *this;
}


template <class T> 
block_diag<T>::~block_diag() {
    for( auto &it : blocks )
        delete it;

    blocks.clear(); //Not sure this really does anythin useful, salts the earth more than anything
}

template <class T> 
void block_diag<T>::generate_data_holder() {
    //It's difficult to make a vector from a class that doesn't have a default constructor, instead we just make an empty vector and fill it
    blocks = vector<matrix<T>*>();

    for(unsigned int i=0; i< block_amount; i++)
        blocks.push_back(new matrix<T>(this->block_dim_r, this->block_dim_c));
}

template <class T> 
T block_diag<T>::get(int i, int j) const {
    if((int)(i/this->block_dim_r) == (int)(j/this->block_dim_c))
        return blocks[(int)(i/this->block_dim_r)]->get(i%this->block_dim_r, j%this->block_dim_c);
    else
        return 0;
}

template <class T> 
void block_diag<T>::set(int i, int j, const T val) {
    //TODO:Error checking
    set_block((int)(i/this->block_dim_r), i%this->block_dim_r, j%this->block_dim_c, val);
}

template <class T> 
void block_diag<T>::set_block(int block_i, int i, int j, const T val) {
    blocks[block_i]->set(i, j, val);
}

template <class T> 
void block_diag<T>::transpose() { 
    for(auto &b : blocks) 
        b->transpose(); 

    int tmp = this->dim_c; 
    this->dim_c = this->dim_r; 
    this->dim_r = tmp; 

    tmp = this->block_dim_c; 
    this->block_dim_c = this->block_dim_r; 
    this->block_dim_r = tmp; 
}

template <class T> 
T* block_diag<T>::get_as_full_array() const {
    printf("Block Diag: Realising full matrix\n");
    T *output = new T[this->block_dim_r*this->block_dim_c*block_amount*block_amount]();

    for(unsigned int b=0; b<block_amount; b++) {
        for(unsigned int i=0; i<this->block_dim_r; i++) {
            for(unsigned int j=0; j<this->block_dim_c; j++) {
                //     Which row is the block on        Row of the block       row in block    col in block
                output[b*this->block_dim_r*this->dim_c + b*this->block_dim_c + i*this->dim_c + j           ] = blocks[b]->get(i,j);
            }
        }
    }

    return output;
}

// Kroneker ////////////////////////////////////////////////////////////////////////
template <class T> 
kroneker<T>::kroneker(int left_r, int left_c, int right_r, int right_c) : left_matrix(left_r, left_c), right_matrix(right_r, right_c) {
    redo_dims();
}

//The rule of three
template <class T> 
kroneker<T>::kroneker(const kroneker<T> &to_clone): left_matrix(to_clone.get_left_matrix()), right_matrix(to_clone.get_right_matrix()) {
    redo_dims();
}

template <class T> 
kroneker<T>& kroneker<T>::operator=(const kroneker<T>& rhs) {
    left_matrix = rhs.get_left_matrix();
    right_matrix = rhs.get_right_matrix();

    redo_dims();

    return *this;
}

template <class T> 
void kroneker<T>::transpose() {
    left_matrix.transpose();
    right_matrix.transpose();

    redo_dims();
}

template <class T> 
void kroneker<T>::redo_dims() {
    this->dim_r = left_matrix.get_dim(Dimension::row)*right_matrix.get_dim(Dimension::row);
    this->dim_c = left_matrix.get_dim(Dimension::column)*right_matrix.get_dim(Dimension::column);
}

template <class T> 
T kroneker<T>::get(int i, int j) const {
    return left_matrix.get((int)i/right_matrix.get_dim(Dimension::row), (int)j/right_matrix.get_dim(Dimension::column)) * right_matrix.get(i%right_matrix.get_dim(Dimension::row), j%right_matrix.get_dim(Dimension::column));
}

template <class T> 
void kroneker<T>::set(int i, int j, T val) {
    printf("Error: Kroneker matrix, what are you setting here?");
}

template<class T> 
vec<T> operator*(const kroneker<T>& lhs, const vec<T>& rhs) {
    vec<T> output = vec<T>(rhs.get_size());

    matrix<T> C_l = lhs.get_left_matrix(); 
    matrix<T> C_r = lhs.get_right_matrix(); 

    int C_l_dim_r = C_l.get_dim(Dimension::row);
    int C_l_dim_c = C_l.get_dim(Dimension::column);

    int C_r_dim_c = C_r.get_dim(Dimension::column);

    //Split the vector into blocks
    vector<vec<T>> cur_subs = vector<vec<T>>(); 
    for(int i=0; i<C_l_dim_c; i++) 
        cur_subs.push_back(rhs.get_sub_vec(i*C_r_dim_c, (i+1)*C_r_dim_c));
   
    //Do all the left matrix work 
    for(int i=0; i<C_l_dim_c; i++) { //Every column block
        for(int j=0; j<C_l_dim_r; j++) { //Every row block
            output.add(j*C_r_dim_c, C_l.get(i, j) *cur_subs[i]); 
            //TODO:The profiler says that there is a ~matrix<T> being called here but I don't understand where from
        }
    }
    
    //Do the right matrix work
    for(int i=0; i<C_l_dim_r; i++) 
        output.set(i*C_r_dim_c, C_r*cur_subs[i]);

    return output;
}



////////////////////////////////////////////////////////////////////////
//mat///////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class T>
mat<T>& mat<T>::eyes() {
    int smallest_dim = this->smallest_dim_size();
    set(T());

    for(int i=0; i<smallest_dim; i++) {
        set(i,i,1.0);
    }
    return *this;
}

template <class T>
mat<T>& mat<T>::eyes(mat<T>& diag) {
    int smallest_dim = this->smallest_dim_size();
    set(T());

    for(int i=0; i<smallest_dim; i++) {
        set(i,i,diag.get(i,0));
    }

    return *this;
}

//I would need to write a lot of stuff to make SVD work for int matricies so I'm just going to leave it out for now since I don't want to call it anyway
/*template <>
mat<double> mat<int>::truncate_SVD() {
    printf("Don't be callin this!\n");
    return mat<double>();
}

template <>
mat<double> mat<double>::truncate_SVD() {
    mat<double> cur_copy = mat<double>(this->get_dim(Dimension::row), this->get_dim(Dimension::column));
    cur_copy.set(*this);

    int s_dim = this->get_dim(Dimension::row) < this->get_dim(Dimension::column) ? this->get_dim(Dimension::row) : this->get_dim(Dimension::column);
    vec<double> S = vec<double>(s_dim);
    mat<double> U = square<double>(this->get_dim(Dimension::row));
    mat<double> VT = square<double>(this->get_dim(Dimension::column));

    int svd_info = LAPACKE_dgesdd( 
            LAPACK_ROW_MAJOR, //int matrix_layout
            'A', //char jobz, 
            this->get_dim(Dimension::row),            //lapack_int m,
            this->get_dim(Dimension::column),            //lapack_int n, 
            cur_copy.get_all(),//double* a, 
            this->get_dim(Dimension::column), //lapack_int lda, 
            S.get_all(),//double* s,
            U.get_all(),//double* u, 
            U.get_dim(Dimension::column), //lapack_int ldu, 
            VT.get_all(),//double* vt,
            VT.get_dim(Dimension::column)//lapack_int ldvt 
            );

    if(svd_info)
        printf("SVD error code:%i\n", svd_info);

    for(int i=0; i<S.get_size(); i++) 
        if(S.get(i) < 0.001)  S.set(i, 0);

    mat<double>(this->get_dim(Dimension::row), this->get_dim(Dimension::column)).eyes(S).print();
    U.print();
    VT.print();
    return U*mat<double>(this->get_dim(Dimension::row), this->get_dim(Dimension::column)).eyes(S)*VT;
}*/

template <class T>
void matrix<T>::multiply_row(int i, double to_mult_by) {
    for(int j=0; j<this->get_dim_c(); j++) {
        set(i, j, get(i,j)*to_mult_by);
    }
}

template <class T>
void matrix<T>::bracket_diag(const mat<T>& main_diag) {
    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            set(i,j, main_diag.get(j,j)*get(i,j)*main_diag.get(i,i));
        }
    }
}

template <class T>
mat<T>& mat<T>::operator/=(const T rhs) {
    for(int i=0; i<get_dim_r(); i++) {
        for(int j=0; j<get_dim_c(); j++) {
            set(i,j, get(i,j)/rhs);
        }
    }

    return *this;
}

template <class T>
mat<T>& mat<T>::operator*=(const T rhs) {
    for(int i=0; i<get_dim_r(); i++) {
        for(int j=0; j<get_dim_c(); j++) {
            set(i,j, get(i,j)*rhs);
        }
    }

    return *this;
}



template <class T>
void matrix<T>::transpose() {
    matrix<T> old = matrix<T>(*this);

    int tmp_dim = this->dim_c;
    this->dim_c = this->dim_r;
    this->dim_r = tmp_dim;

    for(int i=0; i<old.get_dim_r(); i++)
        for(int j=0; j<old.get_dim_c(); j++)
            this->set(j,i, old.get(i,j));
}


template <class T>
T matrix<T>::max() const {
    T maximum = T();

    for(int i=0; i<this->get_dim(Dimension::row); i++) {
        for(int j=0; j<this->get_dim(Dimension::column); j++) {
            if(get(i,j) > maximum) {
                maximum = get(i,j);
            }
        }
    }

    return maximum;
}

template <class T>
T matrix<T>::min() const {
    T minimum = T();

    for(int i=0; i<this->get_dim(Dimension::row); i++) {
        for(int j=0; j<this->get_dim(Dimension::column); j++) {
            if(get(i,j) < minimum) {
                minimum = get(i,j);
            }
        }
    }

    return minimum;
}

//Print out the whole matrix
template <class T>
void mat<T>::print(std::ostream& where_to_output, bool space) const {
    for(int i=0; i<get_dim(Dimension::row); i++) {
        for(int j=0; j<get_dim(Dimension::column); j++) {
            if(space)
                where_to_output << std::setw(9) << std::right << std::setprecision(2) << std::scientific << get(i,j);
            else  {
                if(j!=0)
                    where_to_output << " ";
                where_to_output << get(i,j);
            }
        }
        where_to_output << "\n";
    }

    where_to_output << "\n";
}

template <class T>
bool matrix<T>::check_zero() const {
    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
           if(this->get(i,j) != this->get(i,j) or abs(this->get(i,j)) > this->tol_value)
               return false;
        }
    }

    return true;
}

template <class T>
matrix<T>& matrix<T>::boil_off_smalls() {
    for(int i=0; i<this->get_dim_r(); i++) {
        for(int j=0; j<this->get_dim_c(); j++) {
            if(abs(this->get(i,j)) < this->tol_value)
                this->set(i,j, 0);
        }
    }

    return *this;
}



/////////////////////////////////////////////////////////////////////////////////////
//vec////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//Copy initialise vector
template <class T>
vec<T>::vec(vector<T> in) : matrix<T>(in.size(), 1) {
    for(int i=0; i<get_size(); i++) {
        set(i, in[i]);
    }
}

template <class T>
vec<T>& vec<T>::piecewise_times_in_place(const vec<T>& rhs) {
    for(int i=0; i<get_size(); i++)
        set(i, get(i) * rhs.get(i));

    return *this;
}

template <class T>
vec<T>& vec<T>::operator+=(const vec<T>& rhs) {
    for(int i=0; i<get_size(); i++) {
        set(i, get(i) + rhs.get(i));
    }

    return *this;
}

template <class T>
vec<T>& vec<T>::operator-=(const vec<T>& rhs) {
    for(int i=0; i<this->get_dim_r(); i++) {
        set(i, get(i) - rhs.get(i));
    }

    return *this;
}

template <class T>
vec<T>& vec<T>::operator/=(const T rhs) {
    if(rhs==0) 
        printf("ERROR: Div by 0\n");

    for(int i=0; i<get_size(); i++)
        set(i, get(i)/rhs);

    return *this;
}

template <class T>
void vec<T>::transpose(){
    int tmp = this->get_dim_r();
    this->set_dim_r(this->get_dim_c());
    this->set_dim_c(tmp);
}

/*template <class T>
vec<T>* vec<T>::Tr() const {
    vec<T> *output = new vec<T>(*this);

    output->transpose();

    return output;
}

template <class T>
vec<T> vec<T>::Tra() const {
    vec<T> output = vec<T>(*this);

    output.transpose();

    return output;
}*/


//////////////////////////////////////////////////////////////////////////////////
// Vector ////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
template <class T>
vec<T> vec<T>::piecewise_square() {
    vec<T> output = *this;
    
    int size = get_size();

    for(int i=0; i<size; i++) 
        output.set(i, pow(get(i), 2));
    

    return output;
}

template <class T>
vec<T> vec<T>::piecewise_root() {
    vec<T> output = *this;
    
    int size = get_size();

    for(int i=0; i<size; i++) {
        output.set(i, pow(get(i), 0.5));
    }

    return output;
}

template <class T>
double vec<T>::average() const {
    double running_total = 0;

    int size = get_size();

    for(int i=0; i<size; i++) {
        running_total += get(i);
    }

    return running_total/size; 
}

template <class T>
T vec<T>::sum() const {
    T output = 0;

    int size = get_size();

    for(int i=0; i<size; i++) {
        output += get(i);
    }

    return output;
}


template <class T>
T vec<T>::min() const {
    T minimum = get(0);

    int size = get_size();

    for(int i=1; i<size; i++) {
        if(get(i) < minimum)
            minimum = get(i);
    }

    return minimum;
}

template <class T>
T vec<T>::max() const {
    T maximum = get(0);

    int size = get_size();

    for(int i=1; i<size; i++) {
        if(get(i) > maximum)
            maximum = get(i);
    }

    return maximum;
}

//TODO:Doesn't work for row vecs?
template <class T>
T vec<T>::get(int i) const {
    return (get_prominent_dim() == Dimension::row) ? matrix<T>::get(i, 0) : matrix<T>::get(0, i);
}


template <class T>
void vec<T>::add(int start, const vec<T>& to_add) {
    int add_size = to_add.get_size();
    for(int i=0; i<add_size; i++)
        this->values[(unsigned int)i+(unsigned int)start] += to_add.get(i);
}

template <class T>
int vec<T>::get_size() const { 
    return (this->get_dim(Dimension::row) > this->get_dim(Dimension::column)) ? this->get_dim(Dimension::row) : this->get_dim(Dimension::column); 
}

template <class T>
Dimension vec<T>::get_prominent_dim() const {
    return (this->get_dim(Dimension::row) > this->get_dim(Dimension::column)) ? Dimension::row : Dimension::column; 
}

//Set the vector to the values of another given vector
template <class T>
void vec<T>::set(int start_position, const vec<T>& setting_from) {
    //Get the size of this vector and the one that we're setting from
    int amount_to_set = setting_from.get_size();
    int this_dim = this->get_size();
    
    //Check that we're not setting more values than this vector has space for, if we are just set as many as we have space for
    if(start_position+amount_to_set > this_dim) {
        amount_to_set = this_dim - start_position; 
    }

    //Set the values of this vector to the ones of the other vector 
    for(int i=0; i<amount_to_set; i++) {
        set(start_position + i, setting_from.get(i));
    }
}

template <class T>
vec<T> vec<T>::get_sub_vec(int start, int finish) const {
    if(start > finish || (unsigned int)finish > this->dim_r) 
        printf("Error: vec, accessing %i %i for dimension %i\n", start, finish, this->dim_r);

    vec<T> output = vec<T>(finish-start);

    for(int i=start; i<finish; i++) 
        output.set(i-start, get(i));

    return output;
}

template <class T>
vector<T> vec<T>::as_vector() {
    vector<T> output = vector<T>();
    int s = this->get_size();

    for(int i=0; i<s; i++) {
        output.push_back(get(i));
    }

    return output;
}

template <class T>
T vec<T>::L2_norm() const {
    int output=0;
    for(int i=0; i<get_size(); i++)
        output += pow(get(i),2);
    return sqrt(output);
}

template<class T> 
vec<T> mat_vec_mult(const mat<T>& lhs, const vec<T>& rhs) {
    if(rhs.get_dim(Dimension::row) == 1) //If rhs is row matrix
        printf("Right times by row vector is invalid operator\n");
    vec<T> output = vec<T>(rhs.get_size()); 
    T* lhs_data =  lhs.get_as_full_array();

    cblas_dgemv(
            CblasRowMajor,  //const CBLAS_LAYOUT layout,
            CblasNoTrans,   //const CBLAS_TRANSPOSE TransA,
            lhs.get_dim(Dimension::row), //const int M, 
            lhs.get_dim(Dimension::column), //const int N,
            1,              //const double alpha, 
            lhs_data,  //const double *A, 
            lhs.get_dim(Dimension::column), //const int lda,
            rhs.get_all(),  //const double *X, 
            1,              //const int incX, 
            0,              //const double beta,
            output.get_all(),  //double *Y, 
            1               //const int incY
            );

    delete [] lhs_data;

    return output;
}

template<> 
vec<float> mat_vec_mult(const mat<float>& lhs, const vec<float>& rhs) {
    vec<float> output = vec<float>(lhs.get_dim_r()); 
    float* lhs_data = lhs.get_as_full_array();

    cblas_sgemv(
            CblasRowMajor, //const CBLAS_LAYOUT layout,
            CblasNoTrans, //const CBLAS_TRANSPOSE TransA,
            lhs.get_dim_r(), //const int M, 
            lhs.get_dim_c(), //const int N,
            1, //const double alpha, 
            lhs_data, //const double *A, 
            lhs.get_dim_c(), //const int lda,
            rhs.get_all(), //const double *X, 
            1, //const int incX, 
            0, //const double beta,
            output.get_all(), //double *Y, 
            1 //const int incY
            );

    delete [] lhs_data;

    return output;
}

template<class T>
vec<T> operator*(const matrix<T>& lhs, const vec<T>& rhs) {
    if(rhs.get_dim(Dimension::row) == 1) //If rhs is row matrix
        printf("Right times by row vector is invalid operator\n");
    vec<double> output = vec<double>(lhs.get_dim(Dimension::row)); 

    cblas_dgemv(
            CblasRowMajor,  //const CBLAS_LAYOUT layout,
            CblasNoTrans,   //const CBLAS_TRANSPOSE TransA,
            lhs.get_dim(Dimension::row), //const int M, 
            lhs.get_dim(Dimension::column), //const int N,
            1,              //const double alpha, 
            lhs.get_all(),  //const double *A, 
            lhs.get_dim(Dimension::column), //const int lda,
            rhs.get_all(),  //const double *X, 
            1,              //const int incX, 
            0,              //const double beta,
            output.get_all(),  //double *Y, 
            1               //const int incY
            );

    return output;
}

template<>
vec<float> operator*(const matrix<float>& lhs, const vec<float>& rhs) {
    if((rhs.get_dim(Dimension::row) == 1) and !(rhs.get_dim(Dimension::column) == 1)) //If rhs is row matrix
        printf("Right times by row vector is invalid operator for lots of matricies?\n");
    vec<float> output = vec<float>(lhs.get_dim(Dimension::row)); 

    cblas_sgemv(
            CblasRowMajor,  //const CBLAS_LAYOUT layout,
            CblasNoTrans,   //const CBLAS_TRANSPOSE TransA,
            lhs.get_dim(Dimension::row), //const int M, 
            lhs.get_dim(Dimension::column), //const int N,
            1,              //const double alpha, 
            lhs.get_all(),  //const double *A, 
            lhs.get_dim(Dimension::column), //const int lda,
            rhs.get_all(),  //const double *X, 
            1,              //const int incX, 
            0,              //const double beta,
            output.get_all(),  //double *Y, 
            1               //const int incY
            );

    return output;
}



///////////////////////////////////////////////////////////////////////////
//square///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <class T>
LU_dual<T> square<T>::get_LU_decom() {
    LU_dual<T> output = LU_dual<T>(get_size());
    output.lower.eyes();

    for(int r=0; r<get_size(); r++) {
        for(int i=r; i<get_size(); i++) {

            //Go down all the previous columns and calculate the dot product sum
            double running_sum =0;
            for(int j=0; j<r; j++) {
                running_sum += output.lower.get(r, j)*output.upper.get(j,i);
            }
            output.upper.set(r,i, this->get(r,i) - running_sum);
        } 

        for(int i=r+1; i<get_size(); i++) {
            double running_sum =0;
            for(int j=0; j<r; j++) {
                running_sum += output.lower.get(i, j)*output.upper.get(j,r);
            }
            if(output.upper.get(r,r) !=0) 
                output.lower.set(i,r, (this->get(i,r) - running_sum)/output.upper.get(r,r));
        }
    }

    return output;
}

template <class T> 
square<T> square<T>::grab_main_diagonoal_mat() {
    square<T> output = square<T>(this->get_size());
    for(int i=0; i<output.get_size(); i++) {
        output.set(i,i, this->get(i,i));        
    }

    return output;
}

template <> 
vec<float> square<float>::eigen_values() {
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
        printf("Spectral solver error code:%i\n", spectral_info);

    return evals;
}

template <class T> 
vec<T> square<T>::eigen_values() {
    square<T> cur_copy = square<T>(this->get_size());
    cur_copy.set(*this);

    vec<double> evals = vec<double>(this->get_size());
    vec<double> evals_im = vec<double>(this->get_size());
    square<double> evecs = square<double>(this->get_size());
    int spectral_info = LAPACKE_dgeev( 
            LAPACK_ROW_MAJOR,//int matrix_layout, 
            'N',//char jobvl, Don't compute left eigenvectors, although it's just the right evecs transpose so might be worth it
            'V',//char jobvr, 
            this->get_size(),//lapack_int n, 
            cur_copy.get_all(),//double* a, 
            this->get_size(),//lapack_int lda, 
            evals.get_all(),//double* wr, 
            evals_im.get_all(),//double* wi, 
            NULL,//double* vl, 
            evecs.get_size(),//lapack_int ldvl, 
            evecs.get_all(),//double* vr, 
            evecs.get_size()//lapack_int ldvr 
            );

     if(spectral_info)
        printf("Spectral solver error code:%i\n", spectral_info);

    return evals;
}

template <>
square<float> square<float>::inverse() {
    //There are problems with the types of template clashing here, dgetri can only accept doubles and I want to handle the int case as well
    square<float> output = square<float>(this->get_size());

    for(int i=0; i<output.get_size(); i++) 
        for(int j=0; j<output.get_size(); j++) 
            output.set(i, j, this->get(i,j));

    vec<int> pivot_indicies = vec<int>(get_size());
    int sgetrf_info = LAPACKE_sgetrf(
            LAPACK_ROW_MAJOR, //int matrix_layout, 
            get_size(),
            get_size(),
            output.get_all(),
            get_size(),
            pivot_indicies.get_all()
            );

    if(sgetrf_info) 
        printf("Problem with LU decomposition in inverse matrix method, error code %i\n", sgetrf_info);

    int sgetri_info = LAPACKE_sgetri( 
            LAPACK_ROW_MAJOR, //int matrix_layout, 
            get_size(),
            output.get_all(), 
            get_size(),
            pivot_indicies.get_all()
            );

    if(sgetri_info) 
        printf("Problem with inverse matrix method, error code %i\n", sgetri_info);

    return output;
}


//TODO:Write this function using https://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c#3525136
template <class T> 
square<T> square<T>::inverse() {
    //There are problems with the types of template clashing here, dgetri can only accept doubles and I want to handle the int case as well
    square<T> output = square<T>(this->get_size());
    for(int i=0; i<output.get_size(); i++) 
        for(int j=0; j<output.get_size(); j++) 
            output.set(i, j, this->get(i,j));

    vec<int> pivot_indicies = vec<int>(get_size());
    int dgetrf_info = LAPACKE_dgetrf(
            LAPACK_ROW_MAJOR, //int matrix_layout, 
            get_size(),
            get_size(),
            output.get_all(),
            get_size(),
            pivot_indicies.get_all()
            );

    if(dgetrf_info) 
        printf("Problem with LU decomposition in inverse matrix method, error code %i\n", dgetrf_info);

    int dgetri_info = LAPACKE_dgetri( 
            LAPACK_ROW_MAJOR, //int matrix_layout, 
            get_size(),
            output.get_all(), 
            get_size(),
            pivot_indicies.get_all()
            );

    if(dgetri_info) 
        printf("Problem with inverse matrix method, error code %i\n", dgetri_info);

    return output;
}

template <class T>
double square<T>::get_condition_number() {
    square<double> cur_copy = square<double>(this->get_size());
    for(int i=0; i<cur_copy.get_size(); i++) 
        for(int j=0; j<cur_copy.get_size(); j++) 
            cur_copy.set(i, j, this->get(i,j));

    int m;
    vec<int> supports = vec<int>(2*this->get_size());
    vec<double> evals = vec<double>(this->get_size());
    square<double> evecs = square<double>(this->get_size());

    //TODO:Change condition number method to only compute the eigen values
    int spectral_info = LAPACKE_dsyevr( 
            LAPACK_ROW_MAJOR, //int matrix_layout, 
            'V',//char jobz, 
            'A', //char range, 
            'L', //char uplo,
            this->get_size(), //lapack_int n, 
            cur_copy.get_all(), //double* a, 
            this->get_size(), //lapack_int lda, 
            0, //double vl,
            0, //double vu, 
            0,//lapack_int il, 
            0,//lapack_int iu,
            LAPACKE_dlamch('S'), //double abstol, 
            &m, //lapack_int* m, 
            evals.get_all(), //double* w, 
            evecs.get_all(), //double* z,
            evecs.get_size(), //lapack_int ldz, 
            supports.get_all()//lapack_int* isuppz 
            );

    if(spectral_info)
        printf("Spectral solver error code:%i\n", spectral_info);

    return evals.max()/evals.min();
}

template <class T> 
matrix<double> solve_lin_eq(LU_dual<T> &LU, const vec<double> &rhs, vec<double> &result) {
    //Solve the LU system
    vec<double> lower_sol = vec<double>(rhs.get_size());
    double temp_row_sol = 0;

    //Solve lower system
    for(int r=0; r<rhs.get_size(); r++) {
        temp_row_sol = rhs.get(r);

        for(int c=0; c<r; c++) {
            temp_row_sol -= LU.lower.get(r, c)*lower_sol.get(c);
        }

        lower_sol.set(r, temp_row_sol/LU.lower.get(r,r));
    }

    //Solve upper system
    for(int r=rhs.get_size()-1; r>-1; r--) {
        temp_row_sol=lower_sol.get(r);

        for(int c=rhs.get_size()-1; c>r; c--) {
            temp_row_sol -= LU.upper.get(r, c)*result.get(c);
        }

        result.set(r, temp_row_sol/LU.upper.get(r,r));
    }

    return result;
     
}

///////////////////////////////////////////////////////////////////////////
//symmetric////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//Find the cholesky factor of the matrix (find L s.t. A = L * L^T where A is this matrix and L is lower diagonal)
template <>
square<float> symmetric<float>::cholesky_square_root() {
    square<float> output = square<float>(this->get_size());
    output.set(*this);
    for(int i=0; i<get_dim_r(); i++)
        for(int j=i; j<get_dim_c(); j++)
            if(i!=j)
                output.set(i, j, 0);

    //Look here for how to understand row major and lda (note we use row major) https://stackoverflow.com/questions/34698550/understanding-lapack-row-major-and-lapack-col-major-with-lda
    int info = LAPACKE_spotrf(
            LAPACK_ROW_MAJOR, 
            'L', 
            this->get_size(), 
            output.get_all(), 
            this->get_size()
            );

    if(info)
        throw "Cholesky factorisation error, lapacke info code: "+to_string(info) +"\n";

    return output;
}

//Find the cholesky factor of the matrix (find L s.t. A = L * L^T where A is this matrix and L is lower diagonal)
template <class T>
square<T> symmetric<T>::cholesky_square_root() {
    square<T> output = square<T>(this->get_size());
    output.set(*this);
    for(int i=0; i<this->get_dim_r(); i++)
        for(int j=i; j<this->get_dim_c(); j++)
            if(i!=j)
                output.set(i, j, 0);

    //Look here for how to understand row major and lda (note we use row major) https://stackoverflow.com/questions/34698550/understanding-lapack-row-major-and-lapack-col-major-with-lda
    int info = LAPACKE_dpotrf(
            LAPACK_ROW_MAJOR, 
            'L', 
            this->get_size(), 
            output.get_all(), 
            this->get_size()
            );

    if(info)
        throw "Cholesky factorisation error, lapacke info code: "+to_string(info) +"\n";

    return output;
}


//Find the eigenvalue decomposition square root of the matrix (find U s.t. A = U * U^T where A is this matrix)
//TODO:This isn't optimised in the doubles case, look at the float specialisation
template <class T>
symmetric<T> symmetric<T>::eigen_square_root() {
    square<double> cur_copy = square<double>(this->get_size());
    cur_copy.set(*this);

    vec<double> evals = vec<double>(this->get_size());
    vec<double> evals_im = vec<double>(this->get_size());
    square<double> evecs = square<double>(this->get_size());

    int spectral_info = LAPACKE_dgeev( 
            LAPACK_ROW_MAJOR,//int matrix_layout, 
            'N',//char jobvl, Don't compute left eigenvectors, although it's just the right evecs transpose so might be worth it
            'V',//char jobvr, 
            this->get_size(),//lapack_int n, 
            cur_copy.get_all(),//double* a, 
            this->get_size(),//lapack_int lda, 
            evals.get_all(),//double* wr, 
            evals_im.get_all(),//double* wi, 
            NULL,//double* vl, 
            evecs.get_size(),//lapack_int ldvl, 
            evecs.get_all(),//double* vr, 
            evecs.get_size()//lapack_int ldvr 
            );
    
    if(spectral_info)
        throw "Spectral solver error code:"+ to_string(spectral_info) + "\n";

    //This step requires positive evals
    square<double> D_h = square<double>(evals.get_size());
    for(int i=0; i<D_h.get_size(); i++) {
        if(evals.get(i) < 0)  {
            //if(abs(evals.get(i)) < 1e-3)
            if (evals.get(i) < -this->tol_value*10) printf("%i Negative eigen %f\n", i, evals.get(i));
            D_h.set(i,i,0);
            //else
            //    throw "Negative eval " + to_string(evals.get(i)) + " at position " + to_string(i) + " was just too small to ignore and hence a real square root of the matrix given doesn't exist\n";
        } else {
            D_h.set(i,i,sqrt(evals.get(i)));//sqrt(abs(evals.get(i))));
        }
    }
    
    //If the evecs are orthogonal then we should only need transpose here instead of inverse but this isn't always the case for some reason?
    //You need evecs.inverse() in here or else all the variables get jumbled around. It should still be technically correct but confusing
    symmetric<T> output = symmetric(evecs.get_dim(Dimension::row));
    output.set(evecs * D_h * evecs.inverse());
    return output;
}


template <class T>
void symmetric<T>::eigen_square_root_inplace() {

    vec<double> evals = vec<double>(this->get_size());
    vec<double> evals_im = vec<double>(this->get_size());
    square<double> evecs = square<double>(this->get_size());

    int spectral_info = LAPACKE_dgeev( 
            LAPACK_ROW_MAJOR,//int matrix_layout, 
            'N',//char jobvl, Don't compute left eigenvectors, although it's just the right evecs transpose so might be worth it
            'V',//char jobvr, 
            this->get_size(),//lapack_int n, 
            this->get_all(),//double* a, 
            this->get_size(),//lapack_int lda, 
            evals.get_all(),//double* wr, 
            evals_im.get_all(),//double* wi, 
            NULL,//double* vl, 
            evecs.get_size(),//lapack_int ldvl, 
            evecs.get_all(),//double* vr, 
            evecs.get_size()//lapack_int ldvr 
            );
    
    if(spectral_info)
        throw "Spectral solver error code:"+ to_string(spectral_info) + "\n";

    //This step requires positive evals
    square<double> D_h = square<double>(evals.get_size());
    for(int i=0; i<D_h.get_size(); i++) {
        if(evals.get(i) < 0)  {
            //if(abs(evals.get(i)) < 1e-3)
            if (evals.get(i) < -this->tol_value*10) printf("%i Negative eigen %f\n", i, evals.get(i));
            D_h.set(i,i,0);
            //else
            //    throw "Negative eval " + to_string(evals.get(i)) + " at position " + to_string(i) + " was just too small to ignore and hence a real square root of the matrix given doesn't exist\n";
        } else {
            D_h.set(i,i,sqrt(evals.get(i)));//sqrt(abs(evals.get(i))));
        }
    }
    
    //If the evecs are orthogonal then we should only need transpose here instead of inverse but this isn't always the case for some reason?
    //You need evecs.inverse() in here or else all the variables get jumbled around. It should still be technically correct but confusing
    this->set(0,0, (*this) * D_h * (*this).inverse()); 
}

template <class T>
square<T> symmetric<T>::circulant_eigen_square_root() {
    double N = this->get_dim(Dimension::row);
    square<T> output = square<T>(N);
    printf("Implement me!\n");
    /*vec<T> c = vec<T>(N);

    for(int i=0; i<N; i++) 
        c.set(i, this->get(i,0));

    c.Tra().print();

    
    complex<double> im = complex<double>(0,1);
    for(double j=0; j<N; j++)  {
        complex<double> sum =0;
        for(double k=0; k<N; k++)
            sum += c.get(N-k)*pow(exp(2 * M_PI * im * j/N),k);
        printf("sum %f %f\n", sum.real(), sum.imag());
    }*/


    /*vec<complex> evecs = vec<complex>(this->get_dim(Dimension::row));
    vec<complex> roots_of_unity = vec<complex>(this->get_dim(Dimension::row));

    for(int i=0; i<this->get_dim(Dimension::row); i++) 
        roots_of_unity.set(i, );

    for(int i=0; i<this->get_dim(Dimension::row); i++) 
        evecs.set(i, );*/

    return output;
}


//This method doesn't work very and will often diverge for useful tolerances
/*template <class T>
square<T> symmetric<T>::symmetric_square_root() {
    square<T> output = Ident(this->get_size());
    output.set(*this);
    square<T> update = Ident(this->get_size());
    double tol = 1e-3;
    bool converge = 1;
    for(int i=0; i<20; i++) {
        update = 0.5*(output + (*this) * output.inverse());

        square<T> change= output - update;
        change.print();
        for(int j=0; j < change.get_dim(Dimension::row); j++) {
            for(int k=0; k < change.get_dim(Dimension::column); k++) {
                if(abs(change.get(j,k)) > tol)
                    converge = 0;
            }
        }
        if(converge)
            break;
        else
            converge = 1;


        output = update;
    }


    (*this - output*output).print();
    return output;// * evecs.Tra();
}*/

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Operations //////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

//The default gemm implementation, it's terrible, don't use this if at all possible 
//TODO: This wont work if *this and rhs don't have the same dimension
template<class T> 
matrix<T> operator*(const mat<T>& lhs, const mat<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    T* holding_row = new T[lhs.get_dim(Dimension::column)];
        
    for(int i=0; i<lhs.get_dim(Dimension::row); i++) {
        for(int j=0; j<rhs.get_dim(Dimension::column); j++)
            holding_row[j] = lhs.get(i, j);

        for(int j=0; j<rhs.get_dim(Dimension::column); j++) {
            T tmp = T();

                for(int k=0; k<lhs.get_dim(Dimension::column); k++)
                    tmp += holding_row[k]  * lhs.get(k, j);

            output.set(i, j, tmp);
        }
    }

    delete [] holding_row;

    return output;
}

template<class T> 
matrix<T> operator*(const matrix<T>& lhs, const mat<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    T* rhs_data = rhs.get_as_full_array();
    cblas_dgemm(
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

template<class T>
matrix<T> operator*(const mat<T>& lhs, const matrix<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    T* lhs_data = lhs.get_as_full_array();
    cblas_dgemm(
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


template<class T>
matrix<T> operator*(const matrix<T>& lhs, const matrix<T>& rhs) {
    printf("Implement me!\n");
    return matrix<T>(1,1);
}


//Diagonal////////////////////////////////////////////////////////////////////////
template<class T> 
diagonal<T> operator*(const diagonal<T>& lhs, const diagonal<T>& rhs) {
    diagonal<T> output = diagonal<T>(lhs.get_dim(Dimension::row));

    for(int i=0; i<rhs.get_dim(Dimension::row); i++)
        output.set(i, lhs.get(i) * rhs.get(i));
    
    return output;
}

template<class T> 
vec<T> operator*(const diagonal<T>& lhs, vec<T> rhs) {
    for(int i=0; i<rhs.get_size(); i++)
        rhs.set(i, lhs.get(i) * rhs.get(i));
    
    return rhs;
}

template<class T> 
vec<T> piecewise_times(const vec<T>& lhs, const vec<T>& rhs) {
    vec<T> output = vec<T>(lhs.get_size());

    for(int i=0; i<output.get_size(); i++)
        output.set(i, lhs.get(i)*rhs.get(i));

    return output;
}

template<class T> 
matrix<T> diag_mult(const mat<T>& lhs, const diagonal<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    for(int i=0; i<lhs.get_dim(Dimension::row); i++) {
        for(int j=0; j<rhs.get_dim(Dimension::column); j++) {
            output.set(i, j, lhs.get(i,j)*rhs.get(j));
        }
    }

    return output;
}

template<class T> 
matrix<T> diag_mult(const diagonal<T>& lhs, const mat<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    for(int i=0; i<lhs.get_dim(Dimension::row); i++) {
        for(int j=0; j<rhs.get_dim(Dimension::column); j++) {
            output.set(i, j, lhs.get(i)*rhs.get(i,j));
        }
    }

    return output;
}

template<class T> T diag_square_over_space(const vec<T>& to_square, const diagonal<T>& weight) {
    T output = T();

    for(int i=0; i<to_square.get_size(); i++) 
        output += pow(to_square.get(i), 2) * weight.get(i);

    return output;

}

// Block operators ////////////////////////////////////////////////////////////////////////////////
template<class T> vec<T> operator*(const block<T>& lhs, const vec<T>& rhs) {
    vec<T> output = vec<T>(lhs.get_dim(Dimension::row));
    int block_dim_c = lhs.get_block_dim(Dimension::column);

    vec<T> cur_sub_vec = vec<T>(lhs.get_block_dim(Dimension::row));
    int block_num = -1;

    for(int i=0; i<lhs.get_block_amount(Dimension::row); i++) {
        cur_sub_vec.set(0);

        for(int j=0; j<lhs.get_block_amount(Dimension::column); j++) {
            block_num = lhs.get_block_use(i,j);
            if(block_num >= 0)
                cur_sub_vec += lhs.get_block(block_num) * rhs.get_sub_vec(j*block_dim_c, (j+1)*block_dim_c);
        }
        output.set(i*lhs.get_block_dim(Dimension::row), cur_sub_vec);
    }

    return output;
}

template<class T> matrix<T> operator*(const block<T>& lhs, const matrix<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    int left_block_r = lhs.get_block_dim(Dimension::row);
    int left_block_c = lhs.get_block_dim(Dimension::column);
    int left_amount_r = lhs.get_block_amount(Dimension::row);
    int left_amount_c = lhs.get_block_amount(Dimension::column);

    int right_block_r = left_block_c;
    int right_block_c = rhs.get_dim(Dimension::column);

    int block_num = -1;

    matrix<T> cur_mat = matrix<T>(left_block_r, right_block_c);

    for(int i=0; i<left_amount_r; i++) {
        cur_mat.set(0);

        for(int k=0; k<left_amount_c; k++) {
            block_num = lhs.get_block_use(i, k);

            if(block_num >= 0) {
                cur_mat += lhs.get_block(block_num)*rhs.get(k*right_block_r, 0, (k+1)*right_block_r, right_block_c);
            }
        }

        output.set(i*left_block_r, 0, cur_mat);
    }

    return output;
}

template<class T> block<T> operator*(const block<T>& lhs, const diagonal<T>& rhs) { 
    block<T> output = block<T>(lhs.get_block_dim(Dimension::row), lhs.get_block_dim(Dimension::column), lhs.get_block_amount(Dimension::row), lhs.get_block_amount(Dimension::column));
    int block_dim_c = lhs.get_block_dim(Dimension::column);

    int block_num = -1;

    for(int i=0; i<lhs.get_block_amount(Dimension::row); i++) {
        for(int j=0; j<lhs.get_block_amount(Dimension::column); j++) {
            block_num = lhs.get_block_use(i,j);
            if(block_num >= 0) 
                output.set_block(i,j, lhs.get_block(block_num) * rhs.get_sub_diag(j*block_dim_c, (j+1)*block_dim_c));
        }
    }

    return output;
}

// Block diag operators ////////////////////////////////////////////////////////////////////////////
template<class T> vec<T> operator*(const block_diag<T>& lhs, const vec<T>& rhs) { 
    vec<T> output = vec<T>(lhs.get_dim(Dimension::row));

    for(int i=0; i<lhs.get_block_amount(); i++) {
        output.set(lhs.get_block_dim(Dimension::row)*i, lhs.get_block(i) * rhs.get_sub_vec(i*lhs.get_block_dim(Dimension::column), (i+1)*lhs.get_block_dim(Dimension::column)));
    }

    return output;
}

template<class T> block_diag<T> operator*(const block_diag<T>& lhs, const diagonal<T>& rhs) { 
    block_diag<T> output = block_diag<T>(lhs.get_block_dim(Dimension::row), lhs.get_block_dim(Dimension::column), lhs.get_block_amount());

    for(int i=0; i<lhs.get_block_amount(); i++)
        output.set_block(i, lhs.get_block(i) * rhs.get_sub_diag(i*lhs.get_block_dim(Dimension::column), (i+1)*lhs.get_block_dim(Dimension::column)));
    
    return output;
}

template<class T> matrix<T> operator*(const block_diag<T>& lhs, const matrix<T>& rhs) { 
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));

    int rhs_sub_block_row_dim = lhs.get_block_dim(Dimension::column);
    int rhs_block_row_dim = (int)((double)rhs.get_dim(Dimension::row) / (double)rhs_sub_block_row_dim);

    //Here we only split rhs into lots of different 'long' matrices, there is no natrual split to use, might as well just give BLAS as big a matrix as possible 
    //int rhs_block_col_dim = 1;
    int rhs_sub_block_col_dim = rhs.get_dim(Dimension::column);

    for(int block_row=0; block_row<rhs_block_row_dim; block_row++) {
        output.set(
                block_row*lhs.get_block_dim(Dimension::row), 
                0,
                lhs.get_block(block_row)*rhs.get(
                    block_row*rhs_sub_block_row_dim, 
                    0,
                    (block_row+1)*rhs_sub_block_row_dim, 
                    rhs_sub_block_col_dim
                    )
                );
    }

    return output;
}

template<class T> matrix<T> operator*(const mat<T>& lhs, const diagonal<T>& rhs) { return diag_mult(lhs, rhs); }
template<class T> matrix<T> operator*(const matrix<T>& lhs, const diagonal<T>& rhs) { return diag_mult(lhs, rhs); }
template<class T> matrix<T> operator*(const diagonal<T>& lhs, const mat<T>& rhs) { return diag_mult(lhs, rhs); }
template<class T> matrix<T> operator*(const diagonal<T>& lhs, const matrix<T>& rhs) { return diag_mult(lhs, rhs); }

//Main diagonal of to_square * weight * to_square^T
template<class T> vec<T> square_over_space_main_diag(const kroneker<T>& weight, const block_diag<T>& to_square) {
    vec<T> output = vec<T>(to_square.get_dim(Dimension::row));
    matrix<T> cur = matrix<T>(to_square.get_block_dim(Dimension::row), to_square.get_block_dim(Dimension::row));

    for(int i=0; i<to_square.get_block_amount(); i++) {
        cur = weight.get_left_matrix().get(i,i) * to_square.get_block(i) * weight.get_right_matrix() * Transpose(to_square.get_block(i));
        
        for(int j=0; j<cur.get_dim(Dimension::row); j++)
            output.set(i* cur.get_dim(Dimension::row) + j, cur.get(j,j));
    }

    return output;
}

//Does lhs * rhs_1 x rhs_2
template<class T> matrix<T> operator*(const block<T>& lhs, const kroneker<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    matrix<T> rhs_left_matrix = rhs.get_left_matrix();
    matrix<T> rhs_right_matrix = rhs.get_right_matrix();

    matrix<T> cur_mat = matrix<T>(lhs.get_block_dim(Dimension::row), rhs_right_matrix.get_dim(Dimension::column));
    for(int i=0; i<lhs.get_block_amount(Dimension::row); i++) {
        for(int j=0; j<lhs.get_block_amount(Dimension::column); j++) {
            cur_mat.set(0);

            for(int k=0; k<lhs.get_block_amount(Dimension::column); k++) {
                int block_num = lhs.get_block_use(i,k);
                if(block_num >= 0) 
                    cur_mat += lhs.get_block(block_num)*rhs_left_matrix.get(k,j);
            }

            output.set(i*lhs.get_block_dim(Dimension::row),j*rhs_right_matrix.get_dim(Dimension::column), cur_mat * rhs_right_matrix);
        }
    }

    return output;
}

//Does lhs_1 x lhs_2 * rhs
template<class T> matrix<T> operator*(const kroneker<T>& lhs, const block<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));
    matrix<T> lhs_left_matrix = lhs.get_left_matrix();
    matrix<T> lhs_right_matrix = lhs.get_right_matrix();

    matrix<T> cur_mat = matrix<T>(lhs_right_matrix.get_dim(Dimension::row), rhs.get_block_dim(Dimension::column));
    for(int i=0; i<lhs_left_matrix.get_dim(Dimension::row); i++) {
        for(int j=0; j<rhs.get_block_amount(Dimension::column); j++) {
            cur_mat.set(0);

            for(int k=0; k<rhs.get_block_amount(Dimension::row); k++) {
                int block_num = rhs.get_block_use(k, j);
                if(block_num >= 0) 
                    cur_mat += lhs_left_matrix.get(i, k) * rhs.get_block(block_num);
            }

            output.set(i*lhs_right_matrix.get_dim(Dimension::row) , j*rhs.get_block_dim(Dimension::column), lhs_right_matrix * cur_mat );
        }
    }

    return output;
}

//Performs to_square * weight * to_square^T = B * C_1 x C_2 * B^T
template<class T> matrix<T> square_over_space(const kroneker<T>& weight, const block<T>& to_square) {
    block<T> left_part_output = block<T>(to_square.get_block_dim(Dimension::row), weight.get_right_matrix().get_dim(Dimension::column), to_square.get_block_amount(Dimension::row), weight.get_left_matrix().get_dim(Dimension::column));

    matrix<T> output = matrix<T>(to_square.get_dim(Dimension::row), to_square.get_dim(Dimension::row));

    //Does B * C_1 x C_2
    matrix<T> cur_mat = matrix<T>(to_square.get_block_dim(Dimension::row), weight.get_right_matrix().get_dim(Dimension::column));
    for(int i=0; i<to_square.get_block_amount(Dimension::row); i++) {
        for(int j=0; j<to_square.get_block_amount(Dimension::column); j++) {
            cur_mat.set(0);

            for(int k=0; k<to_square.get_block_amount(Dimension::column); k++) {
                int block_num = to_square.get_block_use(i,k);
                if(block_num >= 0) 
                    cur_mat += to_square.get_block(block_num)*weight.get_left_matrix().get(k,j);
            }

            left_part_output.set_block(i,j, cur_mat * weight.get_right_matrix());
        }
    }

    //Does the rhs part of (B * C_1 x C_2) * B^T
    matrix<T> cur_mat_square = matrix<T>(to_square.get_block_dim(Dimension::row), to_square.get_block_dim(Dimension::row));
    for(int i=0; i<to_square.get_block_amount(Dimension::row); i++) {
        for(int j=0; j<to_square.get_block_amount(Dimension::row); j++) {
            cur_mat_square.set(0);

            for(int k=0; k<to_square.get_block_amount(Dimension::column); k++) {
                int block_num_left = left_part_output.get_block_use(i, k);
                int block_num_right = to_square.get_block_use(j, k);
                
                if(block_num_left >= 0 and block_num_right >= 0)  {
                    cur_mat_square += left_part_output.get_block(block_num_left)*Transpose(to_square.get_block(block_num_right));
                }
            }

            output.set(i*to_square.get_block_dim(Dimension::row),j*to_square.get_block_dim(Dimension::row), cur_mat_square);
        }

    }

    return output;
}

template<class T> matrix<T> operator*(const kroneker<T>& lhs, const block_diag<T>& rhs) {
    matrix<T> output = matrix<T>(lhs.get_dim(Dimension::row), rhs.get_dim(Dimension::column));

    int block_dim_r = lhs.get_right_matrix().get_dim(Dimension::row);
    int block_dim_c = rhs.get_block_dim(Dimension::column);
    int block_amount_r = lhs.get_left_matrix().get_dim(Dimension::row);
    int block_amount_c = lhs.get_left_matrix().get_dim(Dimension::column);

    matrix<T> lhs_left_matrix = lhs.get_left_matrix();
    matrix<T> lhs_right_matrix = lhs.get_right_matrix();

    matrix<T> tmp_matrix = matrix<T>(block_dim_r, block_dim_c);

    for(int j=0; j<block_amount_c; j++) {
         tmp_matrix = lhs_right_matrix * rhs.get_block(j);

         for(int i=0; i<block_amount_r; i++) {
            output.set(
                    i*block_dim_r, 
                    j*block_dim_c, 
                    lhs_left_matrix.get(i,j) * tmp_matrix
                    );
        }
    }
    
    return output;
}


