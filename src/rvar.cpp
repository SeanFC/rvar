#include "rvar.h"

/////////////////////////////////////////////////////////////////////////////
// Conditioned climate minimiser ////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
//Find the observation increment and call the rest of the cost function
double ccm::conditioned_cost_function(Vec_f w) {
    return conditioned_cost_function(w, observation - observer->state_to_observation(scaling_conversion(w)));
}

//Find the current cost from the current w and observaiton increment
double ccm::conditioned_cost_function(Vec_f w, Vec_f cur_obs_increment) {
    //Start timing 
    struct timespec start, end;
    if(TIMING_VERBOSE) clock_gettime(CLOCK_MONOTONIC, &start);
    
    float bac_sec = pow(w.L2_norm(), 2); 
    double obs_sec = diag_square_over_space(cur_obs_increment, obs_weight);

    //Print the obs and bac costs if in verbose mode
    if(SINGULAR_MINIMISE_VERBOSE) {
        int dec_places = 7;
        cout << fixed;
        cout << "Total:" << setw(dec_places*2+1) << setprecision(dec_places) << bac_sec+obs_sec << " Bac:" << setw(dec_places*2 +1) << setprecision(dec_places) << bac_sec << " Obs:" << setw(dec_places*2+1) << setprecision(dec_places) << obs_sec << endl;
    }

    //Print how long the cost took to compute if in verbose timing
    if(TIMING_VERBOSE) {
        clock_gettime(CLOCK_MONOTONIC, &end);
        double elapsed = (end.tv_sec - start.tv_sec);
        elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("In cost %f\n", elapsed);
    }
   
    //Return the final cost, the 0.5 here means that the gradient can cancel out the half 
    return 0.5*(obs_sec + bac_sec);
}

//TODO: The cost function and the minimisers could be reconfigured to use less repeating code
//Given the current set up of ccm, find the state which minimses the cost function
Point_Diagnostics ccm::minimise_cost_function() {
    //Get the current cost, current_result.w should be 0 as we start at the background
    current_result.cost = conditioned_cost_function(current_result.w); 
    
    //Chooce which minimisation method to use
    switch(MINIMISATION_METHOD) {
        //TODO:Nedler mean currently not supported
        case NEDLER_MEAD_MINIMISE:
            {
                //Nedler Mead Simplex Method
                /*
                double var_lim = 1;
                int convg_itter = 2;
                int max_itter = 100000;// *state_size;
                int used_itter = 0;
                int restarts = 0;
                int fflag = 0;

                double *step = new double[state_size];
                for(int i=0; i<state_size; i++) 
                    step[i] = 1;//TODO:What should initial step size be here?
    

                nelmin(
                        &stat_conditioned_cost_call, //Function to minimise
                        state_size, //Dimensions of the function 
                        current_result.w.get_all(), //The output
                        &(current_result.cost), //The minimum value of the function
                        var_lim,    //The terminating limit for the variance of function values 
                        step, //The size and shape of the initial simplex, relative values should reflect units 
                        convg_itter, //How often the convergence check is carried out
                        max_itter, //Maximum number of funciton evaluations
                        &used_itter, //Number of function evaluations used
                        &restarts, //Number of restarts
                        &fflag, //Fault error indicator
                        this);

                delete[] step;
                if(MINIMISE_VERBOSE) {
                    printf(
                            "Minimisation complete \n"
                            "    %i iterations performed\n"
                            "    %i fault\n",
                            used_itter, fflag
                          );
                }*/
                break;
            }
        case LBFGS_MINIMISE:
        default:
            {
                if(MINIMISE_VERBOSE) printf("Starting cost %f\n", current_result.cost);

                //Set the terminating parameters
                lbfgs_parameter_t params_to_use;
                lbfgs_parameter_init(&params_to_use);
                params_to_use.xtol = std::numeric_limits<float>::epsilon(); //Set the precision of lbfgs to be the machine precision
                
                //Limit the number of itterations, useful for debugging
                //params_to_use.max_iterations = 2;

                //Print headers for the progression table 
                if(MINIMISE_VERBOSE) print_progress_table_header();

                //Create and set the working conditioned state
                Vec_d w_d = Vec_d(current_result.w.get_size());    

                for(int i=0; i<current_result.w.get_size(); i++)
                    w_d.set(i, current_result.w.get(i));

                //Perform the L-BFGS method to find the minimum conditioned state 
                int lbfgs_output_code = lbfgs(
                        state_size,                     //Number of variables
                        w_d.get_all(),                  //Array of variables (also the result)
                        &(current_result.cost),         //Final lowest cost
                        stat_conditioned_cost_lbfgs,    //Evaluation function
                        stat_lbfgs_progress,            //Call back function with the progress of the minimisation, NULL for not inuse
                        this,                           //Pass through parameter
                        &params_to_use                  //Settings for the minimisation, NULL for defaults
                        );

                //Set the result to output to the one found by L-BFGS
                for(int i=0; i<current_result.w.get_size(); i++)
                    current_result.w.set(i, w_d.get(i));

                //Print out the result code (and any errors if found) if in verbose mode 
                if(MINIMISE_VERBOSE) {
                    printf("Result Code: %i ", lbfgs_output_code);
                    switch(lbfgs_output_code) {
                        case LBFGSERR_ROUNDING_ERROR:
                            printf("A rounding error occurred; alternatively, no line-search step satisfies the sufficient decrease and curvature conditions.\n"); break;
                        case LBFGSERR_MAXIMUMLINESEARCH:
                            printf("The line-search routine reaches the maximum number of evaluations.\n"); break;
                        case 0:
                            printf("Minimisation finished as planned\n"); break;
                        default:
                            printf("Unknown result code\n");
                    }
                }
            }
    }
    
    //Update the current result based on the minimisation output
    update_current();
    return current_result;
}

//Set the current result based on the state
void ccm::update_current() {
    Vec_d grad = Vec_d(current_result.gradient.get_size());
    for(int i=0; i<grad.get_size(); i++) grad.set(i, current_result.gradient.get(i)); //TODO:Might not actually need to set this?

    //Find the current cost, unconditioned state and analysis error
    calculate_cost_and_grad(current_result.w, grad.get_all());
    current_result.x = scaling_conversion(current_result.w);

    for(int i=0; i<grad.get_size(); i++) current_result.gradient.set(i, grad.get(i)); 
}

//Given a conditioned state, find the cost and gradient of that state
double ccm::calculate_cost_and_grad(Vec_f cur_w, double gradient[]){ 
    //Precalute the terms common to both the cost and gradient functions 
    Vec_f x = scaling_conversion(cur_w);
    Vec_f cur_obs_increment = observation - observer->state_to_observation(x);
    Block_f H_T = Transpose(observer->block_jacobian(x)); 

    //Find the cost and gradient
    double cost = conditioned_cost_function(cur_w, cur_obs_increment);
    calculate_quick_gradient(gradient, cur_w, cur_obs_increment, H_T); 

    return cost;
}

//Given a the conditioned state and some precalculated terms, calculate the graident and put it in the container
//TODO: Pass const Vec_f &cur_w etc instead
//TODO: Change end of minimisation stuff 
void ccm::calculate_quick_gradient(double gradient[], Vec_f cur_w, Vec_f cur_obs_increment, Block_f H_T) {
    //Calculate gradient
    Vec_f grad = cur_w - M * (H_T * (obs_weight * cur_obs_increment));//piecewise_times(obs_weight, cur_obs_increment); 
    
    //Set gradient to container 
    for(int i=0; i<state_size; i++)
        gradient[i] = grad.get(i);  
}

//Find the hessian at a particular state
Square_f ccm::hessian(Vec_f w) {
    Mat_f Ha      = observer->jacobian(scaling_conversion(w));

    return Square_f(Ident_f(M.get_dim(Dimension::row)) + M * Transpose(Ha) * obs_weight * Ha * M);
}

/////////////////////////////////////////////////////////////////////////////////
//Static functions///////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//A funciton to find the cost and gradient of a given conditioned_cost_minimiser class
//Can be used if you need to send the cost and gradient function to another abitrary but don't have a way to send the ccm class along with it
//TODO:I don't know how to just use member functions like you can do with thread. What makes thread special?
double stat_conditioned_cost_call(double x[], void * user_class) {
    ccm *self = static_cast<ccm*>(user_class);

    Vec_f w = Vec_f(self->get_state_size());
    for(int i=0; i<self->get_state_size(); i++) w.set(i, x[i]);

    return self->conditioned_cost_function(w);
}

//A funciton to find the cost and gradient of a given conditioned_cost_minimiser class for use with the L-BFGS function solver
double stat_conditioned_cost_lbfgs(void *user_class, const double x[], double gradient[], int n, double step) {
    ccm *self = static_cast<ccm*>(user_class);
    Vec_f w = Vec_f(self->get_state_size());
    for(int i=0; i<self->get_state_size(); i++) w.set(i, x[i]);
    return self->calculate_cost_and_grad(w, gradient);
} 

//Prints the headings for the progress function when using the L-BFGS minimiser
void print_progress_table_header() {
    std::cout << std::left << std::setw(5) << "N"; 
    std::cout << std::left << std::setw(15) << "Cost";
    std::cout << std::left << std::setw(15) << "Xnorm" ;
    std::cout << std::left << std::setw(15) << "Gnorm" ;
    std::cout << std::left << std::setw(15) << "Step" ;
    std::cout << std::left << std::setw(3) << "Evals" <<"\n";
}

//Prints the progress of the L-BFGS minimser 
int stat_lbfgs_progress(void *user_class, const double x[], const double grad[], const double cur_cost, const double xnorm, const double gnorm, const double step, int n, int k, int ls) {
    if(MINIMISE_VERBOSE) {
        if(mod(k, 20) == 0 and k!=0)
            print_progress_table_header();

        std::cout << std::left << std::setw(5) << k;
        std::cout << std::left << std::setw(15) << cur_cost ;
        std::cout << std::left << std::setw(15) << xnorm ;
        std::cout << std::left << std::setw(15) << gnorm ;
        std::cout << std::left << std::setw(15) << step ;
        std::cout << std::left << std::setw(3) << ls <<"\n";

        if(MINIMISE_VERBOSE >1) {
            printf("    x          grad\n");
            for(int i=0; i<n; i++) printf("%i   %f    %f\n", i, x[i], grad[i]);
        } 
    }
 
    return 0;
}

