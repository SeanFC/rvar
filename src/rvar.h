#ifndef rvar_h 
#define rvar_h 

#include <cmath>
#include <iostream>
#include <fstream>
#include <thread>
#include <limits>

#include "lin_maths.h"              // For mathematical operators
#include "observation_function.h"   // For general observation function defintions
#include "lbfgs.h"                  // For L-BFGS minimser

// Minimiser choices
enum {
    NEDLER_MEAD_MINIMISE, //TODO:Nedler mean currently not supported
    LBFGS_MINIMISE
};

// Define the settings we're going to use
#define MINIMISATION_METHOD LBFGS_MINIMISE //Which minimiser to use

// Which verbose options to use 
#define MINIMISE_VERBOSE 1 
#define SINGULAR_MINIMISE_VERBOSE 0
#define TIMING_VERBOSE 0

// Struct to hold an experiment. Where the experiment is in climate space (x and w) and other data about the current position 
struct Point_Diagnostics {
    double cost = 0;
    Vec_f w, x, gradient;
    
    Point_Diagnostics(int size) : w(size), x(size), gradient(size) {  }
};

// Conditioned Cost Minimiser
// Used to solve a conditioned 3D-Var data assimilation problem 
class ccm {
    protected:
        //Size of the problem
        int state_size, obs_size;

        //Observation error covariance 
        //Square_f obs_weight;
        Diagonal_f obs_weight;

        //Background error covariance
        Symmetric_f M;

        //Observation, background and curren state
        Vec_f observation, x_b;

        //Observation fucntion
        sto *observer = NULL;

        //The current state of the experiment 
        Point_Diagnostics current_result;

    public:
        //Constructor and descructor
        //The largest matrix here should be M 
        ccm(int ss, int os) : state_size(ss), obs_size(os), obs_weight(os), M(ss), observation(os), x_b(ss), current_result(ss) {}
        /*virtual ~ccm() { 
            //TODO: Can't delete observer for some reason
            //if(observer != NULL)
            //    delete observer;
        }*/

        //Cost functions
        double conditioned_cost_function(Vec_f);
        double conditioned_cost_function(Vec_f, Vec_f);

        //Caluclate the gradient along with the cost
        double calculate_cost_and_grad(Vec_f w, double gradient[]);
        void calculate_quick_gradient(double gradient[], Vec_f, Vec_f, Block_f);

        //Minimise a givne cost function
        virtual Point_Diagnostics minimise_cost_function();
    
        //Convert from conditioned state to unconditioend
        Vec_f scaling_conversion(Vec_f w) { 
            return lin_gemv(1.0, M, w, 1.0, x_b);
        }

        //Update the current result of the experiemtn
        void update_current(); 
               
        //Calculate the hessian at a speicifc state 
        Square_f hessian(Vec_f w);
    
        //Accessors
        int get_state_size() { return state_size; }
        int get_obs_size() { return obs_size; }

        Point_Diagnostics& get_current_result() { return current_result; }
};

////////////////////////////////////////////////////////////////////////
//Static functions//////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//Find the cost given a ccm class, for use when using the minimser
double stat_conditioned_cost_call(double* x, void*);
//Find the cost and gradien tgiven a ccm class, for use when using the minimser
double stat_conditioned_cost_lbfgs(void *user_class, const double x[], double gradient[], int n, double step);

//Print the table header for the L-BFGS progress
void print_progress_table_header();
//Print the progress of the L-BFGS minimisation 
int stat_lbfgs_progress(void *user_class, const double x[], const double grad[], const double cur_cost, const double xnorm, const double gnorm, const double step, int n, int k, int ls);
#endif
