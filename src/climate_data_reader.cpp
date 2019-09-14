#include "climate_data_reader.h"

////////////////////////////////////////////////////////////////////////////////////////////////
//Seasonal_correlated_minimiser/////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
// Sets up the x_b, y vectors and B, R matricies so they are non-dimensionalised and ready for minimisation 
Seasonal_correlated_minimiser::Seasonal_correlated_minimiser(
        vector<OBS_accessor::OBS_obs_reading*> &obs, 
        vector<Seasonal_sd_reading*> &bac, 
        vector<Wang_P_Budyko*> &s2o, 
        vector<Wang_P_Budyko*> &full_grid_sto
        ) : 
    ccm(bac.size()*bac[0]->get_size(), obs.size()*obs[0]->get_size()), 
    full_grid_stos(full_grid_sto), 
    view_size(bac.size()*obs[0]->get_size()),
    s2o(s2o),
    bac_sd_eye(state_size),
    C(bac.size(), bac.size(), bac[0]->get_size(), bac[0]->get_size()),
    R(obs_size) {

    //Non-dimensionalise the background in place (This should be done in place as the background can be large and having two of them can be quite awkward)
    Seasonal_space_reading avg_sds = Seasonal_space_reading();
    for(auto &it : bac) {
        avg_sds += it->get_sd();
    }
    avg_sds /= 2*bac.size();

    Seasonal_space_reading new_sd = Seasonal_space_reading();

    // If we're performig a test run and set several variables to specifc values
    bool single_site_test = false;

    Vec_f bac_err_scaling = Vec_f(13);
    if(single_site_test) {
        bac_err_scaling.set(1.5);
        bac_err_scaling.set(7, 2.5);
        bac_err_scaling.set(6, 2.5);
        bac_err_scaling.set(8, 2.5);
    } else {
        bac_err_scaling.set(1);
    }
        
    for(auto &it : bac) {
        new_sd = Seasonal_space_reading((avg_sds + it->get_sd()/(float)2)*bac_err_scaling);
        for(int i=1; i<13; i++)
            new_sd.set(i, new_sd.get(i)*bac_err_scaling.get(i));
        it->set_sd(new_sd);
        MI_alpha_calc::nondimensionalise_state_with_sd(*it);
    }

    //Non-dimensionalise the observations in place 
    for(auto &it : obs) {

        if(single_site_test) {
            (*it).get_sd().set_MAP(0);
            (*it).get_sd().set_MAT(0);

            (*it).set_MTWA(30);
            (*it).set_MTCO(-15);
            (*it).get_sd().set_MTWA(2.29);
            (*it).get_sd().set_MTCO(2.34);
        }

        MI_alpha_calc::nondimensionalise_obs_with_sd(*it);
    }
    
    // Caluclate dimensions of the various spaces 
    amount_of_backgrounds = bac.size();
    amount_of_observations = obs.size();
    inner_state_dim = bac[0]->get_size();
    inner_obs_dim = obs[0]->get_size();

    bac_sd_eye = Diagonal_f(state_size);
    state_location_grid = vector<loc>();

    //Set the background and save all locations of backgrounds
    for(int i=0; i<amount_of_backgrounds; i++) {
        x_b.set(i*inner_state_dim, *bac[i]); 
        state_location_grid.push_back(full_grid_stos[i]->get_position());
    }

    //Set up the sd diagonal matrix 
    for(int i=0; i<amount_of_backgrounds; i++) {
        for(int j=0; j<inner_state_dim; j++) {
            //Set the sd matrix
            if(bac[i]->get_sd().get(j) == 0){
                bac_sd_eye.set(i*inner_state_dim + j, 0); //TODO: If this is sd then this 0 should be inf?
            } else {
                bac_sd_eye.set(i*inner_state_dim + j, bac[i]->get_sd().get(j));
            }
        }
    }

    double cur_sd_inv = 0;
    for(int i=0; i<amount_of_observations; i++) {
        observation.set(i*inner_obs_dim, *obs[i]);
        for(int j=0; j<obs[i]->get_size(); j++) {
            
            //TODO: This should be based on NAN, not 0           
            if(obs[i]->get_sd().get(j) == 0 or obs[i]->get_sd().get(j) != obs[i]->get_sd().get(j)){
                cur_sd_inv = 0;
            } else {
                cur_sd_inv = 1.0/obs[i]->get_sd().get(j);
            }

            obs_weight.set(i*inner_obs_dim + j, cur_sd_inv*cur_sd_inv);

            if(cur_sd_inv == 0) {
                R.set(i*inner_obs_dim + j, 1e10);//1/R.get_tolerance());
            } else {
                R.set(i*inner_obs_dim + j, 1/(cur_sd_inv*cur_sd_inv));
            }                
        }
    }

}

// Delete both the observation functions
Seasonal_correlated_minimiser::~Seasonal_correlated_minimiser() {
    if(observer) delete observer;
    if(no_CO2_observer) delete no_CO2_observer;
}

//Define both the observation error covariance and background error covariance matrices, as well as the observer
void Seasonal_correlated_minimiser::set_scales(double rh_scale, double temp_scale, double grid_scale, int number_of_threads) {

    //Figure out the covariance related to the state grid 
    Square_f state_distance_grid = get_bessel_correlaton_matrix_f(state_location_grid, grid_scale);//.boil_off_smalls();

    //Figure out the covariance related to the temporal grid 
    Square_f mini_covariance = Ident_f(inner_state_dim);
    mini_covariance.set(1,1, get_bessel_correlation_matrix_f(12, temp_scale));

    //Symmetric_f(state_distance_grid).eigen_square_root_inplace(); return;  
    C.set_left_matrix(state_distance_grid);
    C.set_right_matrix(mini_covariance);
    
    M = C;

    //Add in background variance
    M.bracket_diag(bac_sd_eye);

    printf("Square rooting matrix B start\n");

    //Square root the B matrix
    M.eigen_square_root_inplace();

    printf("Square rooting matrix B finished\n");
   
    //Create observation functions 
    observer = new Wang_P_Budyko_grid(s2o, state_location_grid, inner_state_dim, inner_obs_dim, grid_scale, number_of_threads);
    no_CO2_observer = new Wang_P_Budyko_grid(full_grid_stos, state_location_grid, inner_state_dim, inner_obs_dim, grid_scale, number_of_threads);
}

//Find the analysis error associated with a particular state
void Seasonal_correlated_minimiser::calculate_uncertainty_shifted(const Vec_f &x, const Block_Diag_f &G_a, const Block_Diag_f &G_b, Vec_f& ana_err, Vec_f& bac_err) {
    //Precalculate several reused terms
    Block_f Hb    = observer->block_jacobian(x_b);
    Block_f Hb_T  = Transpose(Hb); 
    Block_f Ha    = observer->block_jacobian(x);
    Block_f Ha_T  = Transpose(Ha);

    Vec_f e_i = Vec_f(view_size);//G.get_dim_c());

    Block_f M_1 = Hb * bac_sd_eye; //Hb scaled by the bac uncertainty
    Block_f M_1_T = Transpose(M_1);
    Block_Diag_f M_2 = G_a * bac_sd_eye; //G scaled by the bac uncertainty
    Block_Diag_f M_2_T = Transpose(M_2);
    
    Square_f obs_space_inverse = Square_f(M_1 * (C * M_1_T) + R).inverse();
    
    Vec_f v_1 = Vec_f(state_size);
    Vec_f v_2 = Vec_f(obs_size);
    Vec_f v_3 = Vec_f(obs_size);
    Vec_f bac_err_at_analysis = square_over_space_main_diag(C, M_2);

    // The vector method tries to make the compulation using only vectors, in theory only vector sized objects need to be held in RAM. This might not work fully in practice depending on the lin_maths module
    bool vector_method = false; 
    if(vector_method) {
        for(int i=0; i< view_size; i++) {
            if(i!=0) e_i.set(i-1, 0);
            e_i.set(i, 1);

            v_2 = M_1 * (C * (M_2_T * e_i));
            v_3 = obs_space_inverse * v_2;  //TODO: Use trsv from BLAS instead of full matrix inversion

            ana_err.set(i, bac_err_at_analysis.get(i) - (M_2 * (C * (M_1_T * v_3))).get(i));
        } 
    } 
    // Compute the covariance using a matrix based method. Depending on how lin_maths is setup this could be more viable than the vector method
    else {
        Mat_f A_improve = M_2 * (C * (M_1_T * (obs_space_inverse * (M_1 * (C * M_2_T)))));

        for(int i=0; i<view_size; i++)
            ana_err.set(i, bac_err_at_analysis.get(i) - A_improve.get(i,i)); 
    }

    bac_err = square_over_space_main_diag(C, G_b*bac_sd_eye);
}

// Parse the current state of the minimiser
Seasonal_correlated_minimiser::os_map Seasonal_correlated_minimiser::parse_current_result() { 
    //The output for this function
    os_map output_map = os_map(state_size, obs_size);

    //Get the current result that we'll transform into the output_map
    Point_Diagnostics result = get_current_result(); 

    //The non dim result for state and obs
    Vec_f obs_analysis =        no_CO2_observer->state_to_observation(result.x); 
    Vec_f obs_background =      no_CO2_observer->state_to_observation(x_b);

    Vec_f obs_analysis_sds =    Vec_f(view_size);
    Vec_f obs_background_sds =  Vec_f(view_size);  
    Vec_f attrib =              Vec_f(view_size);  
    
    calculate_uncertainty_shifted(
            result.x, 
            no_CO2_observer->block_diag_jacobian(result.x), 
            no_CO2_observer->block_diag_jacobian(x_b), 
            obs_analysis_sds, 
            obs_background_sds
            );

    //TODO:Re implment considering we no longer actually figure these out
    //for(int i=0; i<state_sds.get_size(); i++) 
    //    state_sds.set(i, sqrt(result.analysis_covariance_main.get(i)));
   
    for(int i=0; i<view_size; i++) {
        //Square root the analysis variance to get sds
        obs_analysis_sds.set(i, sqrt(obs_analysis_sds.get(i)));

        //Square root background variance to get sds
        obs_background_sds.set(i, sqrt(obs_background_sds.get(i)));

        //This is for if GDD5 is 0. In that case it's not possible to know the uncertainty 
        if(i%7 == 6 and obs_analysis_sds.get(i) == 0) {
            obs_analysis_sds.set(i, 0);
        }

        if(obs_background_sds.get(i) != 0)
            attrib.set(i, 100 * (obs_analysis_sds.get(i) / obs_background_sds.get(i)));
    }

    Vec_f state_sds = Vec_f(result.x.get_size());

    //Add the obs and state result to the output maps after dimensionalising them
    for(int i=0; i<result.x.get_size()/inner_state_dim; i++) {
        //Redimensionalise and save state and state sd
        Seasonal_sd_reading current_state = Seasonal_sd_reading(
                result.x.get_sub_vec(i*inner_state_dim, (i+1)*inner_state_dim),
                state_sds.get_sub_vec(i*inner_state_dim, (i+1)*inner_state_dim)
                ); 
        output_map.states.push_back(MI_alpha_calc::dimensionalise_state_with_sd(current_state)); 
        
        //Redimensionalise and save obs and obs sd
        am_sd_reading current_obs = am_sd_reading(
                obs_analysis.get_sub_vec(i*inner_obs_dim, (i+1)*inner_obs_dim),
                obs_analysis_sds.get_sub_vec(i*inner_obs_dim, (i+1)*inner_obs_dim)
                );
        

        MI_alpha_calc::dimensionalise_obs_with_sd(current_obs);
        output_map.obs.push_back(current_obs);

        
        output_map.attrib.push_back(am_space_reading(attrib.get_sub_vec(i*inner_obs_dim, (i+1)*inner_obs_dim)));
    }

    //Set all the relevent outputs
    output_map.x = result.x;
    output_map.w = result.w;
    output_map.cost = result.cost;
    output_map.gradient = result.gradient; 

    // Calculate extra diagnostics
    bool calculate_resolution = false;
    bool calculate_cond = false;

    if(calculate_resolution) {
        Block_f Hb = observer->block_jacobian(x_b);
        
        Mat_f B = M * M;
        Mat_f K = B * Transpose(Hb) * Square_f(Hb * B * Transpose(Hb) + R).inverse(); //TODO: Should this be R or obs_weight?
        output_map.resolution.set(M.inverse() * K * Hb * M);
    }

    //Mat_f Hb = observer->jacobian(x_b);
    if(calculate_cond) {
        Vec_f hess_eigen = hessian(result.w).eigen_values();
        output_map.cond_number = hess_eigen.get(0)/hess_eigen.get(inner_obs_dim-1);
    }

    return output_map;
}

// Get the hessian for a given point in state space
Mat_f Seasonal_correlated_minimiser::get_hessian(Vec_f point_of_hess) {
    //I + B^(T/2) H^T R^(-1) H B^(1/2)
    //B^(1/2) == M
    //R^(1/2) == obs_weight
    //I == Mat_f(state_size, state_size).eyes()
    
    Mat_f obs_jac = observer->jacobian(point_of_hess);
    return Ident_f(state_size) + Transpose(M) * Transpose(obs_jac) * obs_weight * obs_jac * M;
}

// Get the unconditioned condition number at the start of a run
double Seasonal_correlated_minimiser::get_condition_number() {
    return Square_f(get_hessian(x_b)).get_condition_number();
}

// Get the resolution matrix for the current run 
// TODO: Unimplemented
Mat_f Seasonal_correlated_minimiser::get_resolution_matrix() {
    //Mat_f H = observer->jacobian(x_b);
    //Mat_f T = M.Tra() * H.Tra() * obs_weight * H * M;
    //Mat_f K = Square_f(Mat_f(state_size, state_size).eyes() + T).inverse(); 
    //Wang_P_no_CO2_grid wpn2 = Wang_P_no_CO2_grid(full_grid_stos, inner_state_dim, inner_obs_dim);
    //Mat_f H_c = wpn2.jacobian(x_b);
    ////return H_c * M * K * M.Tra() * H.Tra() * obs_weight * H;// * H_c.inverse();
    //set_gain(Vec_f(state_size));
    //return get_gain() * H * M;

    return Mat_f(1,1);
}

//////////////////////////////////////////////////////////////////////////////
//pmip seasonal assimilator///////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Constructor for the assimilator
// Set up all the data and observation operators for a given experiment
pmip_seasonal_assimilator::pmip_seasonal_assimilator(
        Single_PMIP_accessor &back,
        OBS_accessor &o,
        CRU_accessor &mod,
        int grid_de,
        double lat_max,
        double lat_min,
        double lon_max,
        double lon_min,
        int number_of_threads,
        int time_period_id,
        bool use_CO2,
        bool exclude_PSM_obs
        ) : background(back), obs(o), modern(mod), grid_degrade(grid_de), lat_max(lat_max), lat_min(lat_min), lon_max(lon_max), lon_min(lon_min), number_of_threads(number_of_threads), use_CO2(use_CO2) {
    
    // Get all the observations
    obs_readings = obs.get_area(lat_max, lat_min, lon_max, lon_min, time_period_id);

    // Exclude any observations that have the CO2_corrected flag. This is because these observations may have correlated errors between variables
    if(exclude_PSM_obs) {
        for(vector<OBS_accessor::OBS_obs_reading*>::iterator it= obs_readings.begin(); it != obs_readings.end();) {
            if((*it)->get_CO2_corrected()) {
                it = obs_readings.erase(it);
            } else {
                ++it;
            }
        }
    }

    if(obs_readings.size() < 1)
        printf("ERROR: There aren't any observations so no assimilation can be performed\n");

    // Set up the modern and background data 
    all_back_readings = background.get_area(lat_max, lat_min, lon_max, lon_min, grid_degrade);
    back_reading_position_ids = vector<int>(); 
    back_readings = vector<Seasonal_sd_reading*>(); 
    modern_readings = modern.get_area(lat_max, lat_min, lon_max, lon_min, grid_degrade);
    stos = vector<Wang_P_Budyko*>(); 
    state_grid_stos = vector<Wang_P_Budyko*>(); 
    
    double distance_cutoff = 0.1;

    // These need to only be available to the background where there are modern observations
    for(unsigned int i=0; i<modern_readings.size(); i++) {
        for(unsigned int j=0; j<all_back_readings.size(); j++) {
             if(all_back_readings[j]->get_loc().equal_flat(modern_readings[i]->get_loc()) && all_back_readings[j]->get_ice_percentage() == 0) {

                for(unsigned int k=0; k<obs_readings.size(); k++) {
                    double distance = all_back_readings[j]->get_loc().angle_distance_between(obs_readings[k]->get_loc());

                    if(scaled_bessel_correlation(distance, 400.0, EARTH_RADIUS_KM) > distance_cutoff) {
                        //We never use CO2 correction here. This is just for the final transform.
                        state_grid_stos.push_back(
                                new Wang_P_Budyko(
                                    *(modern_readings[i]), 
                                    *(all_back_readings[j]), 
                                    get_orbital(time_period_id), 
                                    modern_readings[i]->get_loc(), 
                                    false)
                                );

                        back_readings.push_back(all_back_readings[j]->get_seasonal_reading_part());
                        back_reading_position_ids.push_back(all_back_readings[j]->get_position_id());

                        break;
                    }
                }

                break;
            } 
        }
    }

    printf("Bac R:%lu, Sta:%lu, Bac:%lu, Obs:%lu\n", all_back_readings.size(), modern_readings.size(), back_readings.size(), obs_readings.size());


    for(unsigned int i=0; i<obs_readings.size(); i++) {
        int closest_mod_index = 0;
        loc current_obs_loc = obs_readings[i]->get_loc(); 
        double smallest_mod_distance = 9000000;
        for(unsigned int j=0; j<modern_readings.size(); j++) {
            double current_distance = modern_readings[j]->get_loc().distance_between(current_obs_loc); 
            if(current_distance < smallest_mod_distance) {
                closest_mod_index = j;
                smallest_mod_distance = current_distance;
            }
        }

        //TODO: Something needs to delete these Wang_P_Budykos    
        stos.push_back(new Wang_P_Budyko(
                    *(modern_readings[closest_mod_index]), 
                    *(modern_readings[closest_mod_index]), 
                    get_orbital(time_period_id), 
                    obs_readings[i]->get_loc(),
                    use_CO2 && !(obs_readings[i]->get_CO2_corrected())
                    )); //Note that we can't use current_obs_loc here since it will go out of scope after this
    }

    scm = new Seasonal_correlated_minimiser(obs_readings, back_readings, stos, state_grid_stos);
}

// Delete all the readings information on the deconstructor 
pmip_seasonal_assimilator::~pmip_seasonal_assimilator() {
    for(auto &it : obs_readings) delete it;
    for(auto &it : all_back_readings) delete it;
    for(auto &it : back_readings) delete it;
    for(auto &it : state_grid_stos) delete it;
    for(auto &it : modern_readings) delete it;
    for(auto &it : stos) delete it;

    delete scm;
}

// Find the condition number for some particular scales at the start of a run
double pmip_seasonal_assimilator::get_condition_number(double rh_scale, double temp_scale, double grid_scale) {
    scm->set_scales(rh_scale, temp_scale, grid_scale, number_of_threads);
    return scm->get_condition_number();
}

// Find the resolution matrix for some particular scales at the start of a run
Mat_f pmip_seasonal_assimilator::get_resolution_matrix(double rh_scale, double temp_scale, double grid_scale) {
    scm->set_scales(rh_scale, temp_scale, grid_scale, number_of_threads);
    return scm->get_resolution_matrix();
}

// Find the analysis for some particular scales at the start of a run
Seasonal_correlated_minimiser::os_map pmip_seasonal_assimilator::static_analysis(double rh_scale, double temp_scale, double grid_scale) {
    //Set the relevant grid scales
    scm->set_scales(rh_scale, temp_scale, grid_scale, number_of_threads);
    scm->update_current();
    return scm->parse_current_result();
}

// Find the analysis for a run where the cost function has been fully minimised
Seasonal_correlated_minimiser::os_map pmip_seasonal_assimilator::find_analysis(double rh_scale, double temp_scale, double grid_scale) {
    //Set the relevant grid scales
    scm->set_scales(rh_scale, temp_scale, grid_scale, number_of_threads);

    //Minimise the cost function and find the analysis
    scm->minimise_cost_function();
    return scm->parse_current_result();
}

// Save an experiment to the database and save extra information to a file
void pmip_seasonal_assimilator::save_analysis(ANA_experiment &analysis, Seasonal_correlated_minimiser::os_map computed_ana) {
    printf("Saving to database\n");
    //Save the important data in the database
    for(unsigned int i=0; i<computed_ana.states.size(); i++) {
        analysis.add_seasonal_record(
                back_reading_position_ids[i], 
                computed_ana.states[i], 
                computed_ana.obs[i], 
                computed_ana.obs[i].get_sd(), 
                computed_ana.attrib[i]
                );
    }
   
    // Whether to save the diagnostic information to a file
    // TODO: Move this to not be selected at compile time
    bool save_file = false;  

    if(save_file) {
        printf("Saving to file\n");
        ofstream outfile("res/output" + to_string(analysis.get_experiment_id()) + ".rvar");
        
        outfile << computed_ana.cost << "\n\n";
        outfile << computed_ana.cond_number << "\n\n";
        //computed_ana.x.print(outfile, false);
        //computed_ana.w.print(outfile, false);
        //computed_ana.gradient.print(outfile, false);
        //scm->get_hessian(computed_ana.x).print(outfile, false);
        //computed_ana.H_a.print(outfile, false);
        computed_ana.resolution.print(outfile, false);
        outfile.close();
    }

    //Delete readings vectors
}
