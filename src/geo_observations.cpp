#include "geo_observations.h"

////////////////////////////////////////////////////////////////////////////
//MI_alpha_calc/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Does the simple averaging part of the sto function
//TODO: Update this with the improved average function (seen in process_nondim_averages)
am_space_reading& MI_alpha_calc::process_averages(am_space_reading &output_obs, const Vec_f &temps, float MAP) {

    output_obs.set_MTCO(temps.min());
    output_obs.set_MTWA(temps.max());

    float gdd5 = 0;
    float mat = 0;
    for(int i=0; i<temps.get_size(); i++) {
        mat += temps.get(i);//*days_in_month[i];
        if(temps.get(i) > 5) 
            gdd5 += (temps.get(i) - 5)*days_in_month[i]; 
    }

    output_obs.set_MAT(mat/12.0);///number_of_days_in_year);
    output_obs.set_GDD5(gdd5);
    output_obs.set_MAP(MAP>0 ? MAP : 0);

    return output_obs;
}

// Calculate the non-dimensional averages
am_space_reading& MI_alpha_calc::process_nondim_averages(am_space_reading &output_obs, const Vec_f &temps, float MAP) {
    output_obs.set_MTCO(temps.min());
    output_obs.set_MTWA(temps.max());

    float running_avg = 0;
    float gdd5 = 0;
    float temp_cutoff = 5/temp_scaler;
    
    for(int i=0; i<temps.get_size(); i++) {
        running_avg += temps.get(i)*days_in_month[i];

        if(temps.get(i) > temp_cutoff) 
            gdd5 += (temps.get(i) - temp_cutoff)*days_in_month[i]; 
    }

    output_obs.set_MAT(running_avg/number_of_days_in_year);
    output_obs.set_GDD5(gdd5/number_of_days_in_year);
    output_obs.set_MAP(MAP>0 ? MAP : 0);

    return output_obs;
}

//Calculates the obs without doing all the CO2 correction stuff
am_space_reading MI_alpha_calc::sto_no_CO2(const Seasonal_Modern_reading &seasonal_climate) {
    am_space_reading output_obs = am_space_reading(); 

    //Do the averages
    process_nondim_averages(output_obs, seasonal_climate.get_temp(), seasonal_climate.get_MAP());

    //Get the straight MI
    output_obs.set_moisture(Moisture_Corrector::true_mi(seasonal_climate.get_MAP(), seasonal_climate.get_temp(), seasonal_climate.get_sunf(), orb, location));

    //Get the alpha from the MI
    output_obs.set_alpha(budyko_relationship(output_obs.get_moisture()));
    return output_obs;
}

///////////////////////////////////////////////////////////////////////////////////////
//Dimensional and Non Dim functions////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
am_space_reading& MI_alpha_calc::dimensionalise_obs(am_space_reading& output_obs) {
    output_obs.set_MAT(temp_scaler*output_obs.get_MAT());
    output_obs.set_GDD5(number_of_days_in_year * temp_scaler *output_obs.get_GDD5());
    output_obs.set_MTCO(temp_scaler*output_obs.get_MTCO());
    output_obs.set_MTWA(temp_scaler*output_obs.get_MTWA());
    
    output_obs.set_MAP(Wang_Dimension::dim_precip(output_obs.get_MAP()));

    output_obs.set_moisture(output_obs.get_moisture());
    output_obs.set_alpha(output_obs.get_alpha());

    return output_obs;
}

am_sd_reading& MI_alpha_calc::dimensionalise_obs_with_sd(am_sd_reading& output_obs) {
    //Set all the sd values. Note that we need the original values so we set them after 
    output_obs.get_sd().set_MAT(temp_scaler*output_obs.get_sd().get_MAT());
    output_obs.get_sd().set_GDD5(number_of_days_in_year*temp_scaler*output_obs.get_sd().get_GDD5());
    output_obs.get_sd().set_MTCO(temp_scaler*output_obs.get_sd().get_MTCO());
    output_obs.get_sd().set_MTWA(temp_scaler*output_obs.get_sd().get_MTWA());
    
    output_obs.get_sd().set_MAP(Wang_Dimension::dim_precip_derivative(output_obs.get_MAP()) * output_obs.get_sd().get_MAP());

    output_obs.get_sd().set_moisture(output_obs.get_sd().get_moisture());
    output_obs.get_sd().set_alpha(output_obs.get_sd().get_alpha());

    //Set the regular values using a sub function
    output_obs.set(dimensionalise_obs(output_obs));

    return output_obs; 
}

am_space_reading& MI_alpha_calc::nondimensionalise_obs(am_space_reading& output_obs) {
    output_obs.set_MAT(output_obs.get_MAT()/temp_scaler);
    output_obs.set_GDD5(output_obs.get_GDD5()/temp_scaler/number_of_days_in_year);
    output_obs.set_MTCO(output_obs.get_MTCO()/temp_scaler);
    output_obs.set_MTWA(output_obs.get_MTWA()/temp_scaler);

    output_obs.set_MAP(Wang_Dimension::nondim_precip(output_obs.get_MAP()));

    output_obs.set_moisture(output_obs.get_moisture());
    output_obs.set_alpha(output_obs.get_alpha());


    return output_obs;
}

am_sd_reading& MI_alpha_calc::nondimensionalise_obs_with_sd(am_sd_reading& output_obs) {

    //Set all the sd values. Note that we need the original values so we set them after 
    output_obs.get_sd().set_MAT(output_obs.get_sd().get_MAT()/temp_scaler);
    output_obs.get_sd().set_GDD5(output_obs.get_sd().get_GDD5()/temp_scaler/number_of_days_in_year);
    output_obs.get_sd().set_MTCO(output_obs.get_sd().get_MTCO()/temp_scaler);
    output_obs.get_sd().set_MTWA(output_obs.get_sd().get_MTWA()/temp_scaler);

    output_obs.get_sd().set_moisture(output_obs.get_sd().get_moisture());
    output_obs.get_sd().set_alpha(output_obs.get_sd().get_alpha());

    output_obs.get_sd().set_MAP(Wang_Dimension::nondim_precip_derivative(output_obs.get_MAP()) * output_obs.get_sd().get_MAP());

    //Set the regular values using a sub function
    output_obs.set(nondimensionalise_obs(output_obs));

    return output_obs; 
}

Seasonal_space_reading& MI_alpha_calc::dimensionalise_state(Seasonal_space_reading& output_state) {
    Vec_f real_tp = Vec_f(output_state.get_temp().get_size());
    for(int i=0; i<real_tp.get_size(); i++) {
        real_tp.set(i, temp_scaler*output_state.get_temp().get(i));
    }

    output_state.set_temp(real_tp);
    output_state.set_MAP(Wang_Dimension::dim_precip(output_state.get_MAP()));
    return output_state;
}

Seasonal_sd_reading& MI_alpha_calc::dimensionalise_state_with_sd(Seasonal_sd_reading& output_state) {
    //Set all the sd values. Note that we need the original values so we set them after 
    //Vec_f fake_rh_sd = Vec_f(output_state.get_RH().get_size());
    Vec_f fake_tp_sd = Vec_f(output_state.get_temp().get_size());

    for(int i=0; i<fake_tp_sd.get_size(); i++) {
        //fake_rh_sd.set(i, output_state.get_sd().get_RH().get(i)/abs(output_state.get_RH().get(i)+1));
        fake_tp_sd.set(i, output_state.get_sd().get_temp().get(i) * temp_scaler);
    }

    output_state.get_sd().set_temp(fake_tp_sd);
    
    output_state.get_sd().set_MAP(output_state.get_sd().get_MAP() * Wang_Dimension::dim_precip_derivative(output_state.get_MAP()));

    //Set the regular values using a sub function
    output_state.set(dimensionalise_state(output_state));

    return output_state;
}

Seasonal_space_reading& MI_alpha_calc::nondimensionalise_state(Seasonal_space_reading& output_state) {
    //Vec_f fake_rh = Vec_f(output_state.get_RH().get_size());
    Vec_f fake_tp = Vec_f(output_state.get_temp().get_size());

    for(int i=0; i<fake_tp.get_size(); i++) {
        fake_tp.set(i, Wang_Dimension::nondim_temp(output_state.get_temp().get(i))); 
    }

    output_state.set_temp(fake_tp);
    output_state.set_MAP(Wang_Dimension::nondim_precip(output_state.get_MAP()));

    return output_state;
}

Seasonal_sd_reading& MI_alpha_calc::nondimensionalise_state_with_sd(Seasonal_sd_reading& output_state) {
    //Set all the sd values. Note that we need the original values here so we change them after
    Vec_f fake_tp_sd = Vec_f(output_state.get_temp().get_size());

    for(int i=0; i<fake_tp_sd.get_size(); i++) {
        fake_tp_sd.set(i, output_state.get_sd().get_temp().get(i) /temp_scaler); //TODO: Use Wang_Dimension function?
    }

    output_state.get_sd().set_temp(fake_tp_sd);
    output_state.get_sd().set_MAP(output_state.get_sd().get_MAP() * Wang_Dimension::nondim_precip_derivative(output_state.get_MAP()));

    //Set the regular values after using a sub function
    nondimensionalise_state(output_state);

    return output_state;
}

////////////////////////////////////////////////////////////////////////////
//Wang_P_Budyko/////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// The full observation function
am_space_reading Wang_P_Budyko::PB_obs_func(const Seasonal_space_reading &input_state) {
    am_space_reading output_obs = am_space_reading(); 

    //Do all the averaging part of the sto
    process_nondim_averages(output_obs, input_state.get_temp(), input_state.get_MAP());

    //TODO:Use the correct orbitals here
    //The wang correction which we invert to uncorrect moisture index
    output_obs.set_moisture(
            mc.uncorrect_mi(
                input_state.get_MAP(), 
                input_state.get_temp())
            );

    if(output_obs.get_moisture() != output_obs.get_moisture())
        printf("Weird moisture stuff\n");

    //Budyko Relationship
    output_obs.set_alpha(budyko_relationship(output_obs.get_moisture()));

    return output_obs;
}

// The derivative of the alpha MI calculation, with the CO2 correction applied
am_space_reading Wang_P_Budyko::PB_obs_func_derivative(const Seasonal_space_reading &point, const am_space_reading &reference_point, int index) {
    am_space_reading output_obs = am_space_reading(); 
    if(index == 0) {
        output_obs.set_MAP(1);
        //output_obs.set_moisture(reference_point.get_moisture());
        if(point.get_MAP() < 1)
            output_obs.set_moisture(reference_point.get_moisture());
        else 
            output_obs.set_moisture(reference_point.get_moisture()/point.get_MAP());

    } else {
        output_obs.set_MAP(0); 
        Vec_f temps = point.get_temp();
        if(temps.max() == temps.get(index-1))
            output_obs.set_MTWA(1);
        if(temps.min() == temps.get(index-1))
            output_obs.set_MTCO(1);
        
        output_obs.set_MAT((float)days_in_month[index-1]/ (float)number_of_days_in_year);
        if(temps.get(index-1) > 5/temp_scaler) 
            output_obs.set_GDD5(output_obs.get_MAT()); //TODO:Double check this with the thesis implementation

        float x = temps.get(index-1);
        float h = sqrt(eps) * x;
        volatile float xph = x + h;
        float dx = xph - x;

        temps.set(index-1, xph);

        //float old_mi_wrt_T = (mc.uncorrect_mi(point.get_MAP(), temps)- reference_point.get_moisture())/dx;
        output_obs.set_moisture((mc.uncorrect_mi(point.get_MAP(), temps) - reference_point.get_moisture())/dx);

        //New derivative/////
        //output_obs.set_moisture(mc.uncorrect_mi_wrt_T(reference_point.get_moisture(), point.get_MAP(), temps.get(index), index-1));
        //printf("This should be much closer to 0: %f\n", output_obs.get_moisture() - old_mi_wrt_T);
    }
    output_obs.set_alpha(budyko_relationship_derivative(reference_point.get_moisture()) * output_obs.get_moisture());

    return output_obs;
}

// The derivative of the full observation function, without CO2 correction applied 
am_space_reading Wang_P_Budyko::PB_obs_func_no_CO2(const Seasonal_space_reading &input_state) {
    Seasonal_Modern_reading full_read = Seasonal_Modern_reading(); 
    full_read.set_MAP(input_state.get_MAP());
    full_read.set_temp(input_state.get_temp());
    full_read.set_RH(sta_model_climate.get_RH());
    full_read.set_sunf(sta_model_climate.get_sunf());

    return sto_no_CO2(full_read);
}

// The derivative of the alpha MI function, without CO2 correction applied 
am_space_reading Wang_P_Budyko::PB_obs_func_no_CO2_derivative(const Seasonal_space_reading &point, const am_space_reading &reference_point, int index) {
    am_space_reading output_obs = am_space_reading(); 
    if(index == 0) {
        output_obs.set_MAP(1);

        //output_obs.set_moisture(reference_point.get_moisture());
        if(point.get_MAP() < 1)
            output_obs.set_moisture(reference_point.get_moisture());
        else 
            output_obs.set_moisture(reference_point.get_moisture()/point.get_MAP());
    } else {
        output_obs.set_MAP(0); 
        Vec_f temps = point.get_temp();
        if(temps.max() == temps.get(index-1))
            output_obs.set_MTWA(1);
        if(temps.min() == temps.get(index-1))
            output_obs.set_MTCO(1);

        output_obs.set_MAT((float)days_in_month[index-1]/ (float)number_of_days_in_year);
        if(temps.get(index-1) > 5/temp_scaler) 
            output_obs.set_GDD5(output_obs.get_MAT());

        float x = temps.get(index-1);
        float h = sqrt(eps) * x;
        volatile float xph = x + h;
        float dx = xph - x;

        temps.set(index-1, xph);

        
        output_obs.set_moisture((mc.true_mi(point.get_MAP(), temps, sta_model_climate.get_sunf(), orb, location)- reference_point.get_moisture())/dx);

        //New derivative/////
        //output_obs.set_moisture(mc.uncorrect_mi_wrt_T(reference_point.get_moisture(), point.get_MAP(), temps.get(index), index-1));
        //printf("This should be much closer to 0: %f\n", output_obs.get_moisture() - old_mi_wrt_T);
    }
    output_obs.set_alpha(budyko_relationship_derivative(reference_point.get_moisture()) * output_obs.get_moisture());

    return output_obs;
}

//////////////////////////////////////////////////////////////////////////////////
//Wang P Budyko Grid//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
Wang_P_Budyko_grid::Wang_P_Budyko_grid(vector<Wang_P_Budyko*> &mini_sto, vector<loc> &location_grid, int state_d, int obs_d, float dist_factor, int number_of_threads) : sto(number_of_threads), grid_sto(state_d, obs_d, dist_factor, mini_sto.size(), location_grid.size()), mini_stos(mini_sto) {
    
     //This weights each observation to its closest state, basd on the great circle distance 
    int lowest_index = 0;
    float lowest_distance = 1e10;
    float cur_distance = 0;

    //Weight the nearest state to the observation only
    for(unsigned int i=0; i<mini_sto.size(); i++) {
        lowest_index = 0;
        lowest_distance = 1e10;

        for(unsigned int j=0; j<location_grid.size(); j++) {
            cur_distance = location_grid[j].angle_distance_between(mini_sto[i]->get_position());

            if(cur_distance < lowest_distance) {
                lowest_distance = cur_distance;
                lowest_index = j;
            }
        }
    
        location_matrix.set(i, lowest_index, 1);
    }
}

//Note here that there are different results for different obs even if they use the same grid cell. One can't just figure the mapping out from one grid cell and then copy it to every observation here (due to the MI calculation being slightly position depedant and the obs being at slightly different places). 
Vec_f Wang_P_Budyko_grid::state_to_observation(const Vec_f &state) {
    int obs_grid_length = location_matrix.get_dim(Dimension::row);
    int state_grid_length = state.get_size()/state_dim;
    Vec_f obs_on_obs_grid = Vec_f(obs_grid_length*obs_dim);
    
    for(int i=0; i<obs_grid_length; i++) {
        //Find the state readings that represent the site we're interested in 
        //Note that this is sort of matrix mutiplication so this function should really be passed to lin maths somehow
        int j=0;
        for(; j<state_grid_length; j++) {
            if(location_matrix.get(i,j) > 0)
                break;
        }
        
        //Convert state readings to obs readings on the obs grid
        obs_on_obs_grid.set(i*obs_dim, mini_stos[i]->state_to_observation(state.get_sub_vec(j*state_dim, (j+1)*state_dim)));
    }

    //Output result grid cells in obs space
    return obs_on_obs_grid;
}

Vec_f Wang_P_Budyko_grid::derivative(const Vec_f &point, const Vec_f &reference_run, int index) {
    int obs_grid_length = location_matrix.get_dim(Dimension::row);
    int full_obs_dim = obs_dim*obs_grid_length;
    Vec_f output = Vec_f(full_obs_dim);
    int state_number = floor(index/state_dim);

    for(int i=0; i<obs_grid_length; i++) {
        if(location_matrix.get(i, state_number) > 0) {
            output.set(i*obs_dim, mini_stos[i]->derivative(
                        point.get_sub_vec(state_number*state_dim, (state_number + 1)*state_dim),
                        reference_run.get_sub_vec(i*obs_dim, (i+1)*obs_dim),
                        index%state_dim
                        ));
        } 
    }

    return output;
}

// Construct jacobian in block matrix form
//TODO: Make sure that all the locations aren't just being filled in, only the blocks where a obserbation is present should be filled in 
Block_f Wang_P_Budyko_grid::block_jacobian(const Vec_f &point) {
    int obs_grid_length = location_matrix.get_dim(Dimension::row);
    int state_grid_length = location_matrix.get_dim(Dimension::column);

    Block_f output = Block_f(obs_dim, state_dim, obs_grid_length, state_grid_length);

    for(int i=0; i<obs_grid_length; i++) {
        for(int j=0; j<state_grid_length; j++) {
            if(location_matrix.get(i, j) > 0) {
                output.set_block(i, j,  mini_stos[i]->jacobian(point.get_sub_vec(j*state_dim, (j+1)*state_dim)));
            }
        } 
    }

    return output;
}

// Construct jacobian in block diagonal matrix form
Block_Diag_f Wang_P_Budyko_grid::block_diag_jacobian(const Vec_f &point) {
    int obs_grid_length = location_matrix.get_dim(Dimension::row);

    Block_Diag_f output = Block_Diag_f(obs_dim, state_dim, obs_grid_length);
    for(int i=0; i<obs_grid_length; i++) {
        if(location_matrix.get(i, i) > 0) {
            output.set(i*obs_dim, i*state_dim, mini_stos[i]->jacobian(point.get_sub_vec(i*state_dim, (i+1)*state_dim)));
        } 
    }

    return output;
}
