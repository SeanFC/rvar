#include "main.h"

// Load information into the SQL database 
// TODO: These must be done in a certain order as several tables are dependant on each other, it maybe take some fiddling to start the tables from scratch
void load_up_tables(Assim_input_parser input_values) {
    sql_database db = sql_database(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address, input_values.tunneling, input_values.tunneling ? input_values.connection_port : input_values.sql_port, input_values.sql_verbose);

    //INFO_accessor inf(db);
    //inf.load_up_tables();
    
    CRU_accessor cru(db);
    //cru.load_up_tables();
    
    PMIP_accessor pmip_full(db);
    //pmip_full.load_in_ice();
    //pmip_full.load_up_tables(cru);
               
    OBS_accessor obs(db);
    obs.load_up_tables(cru);
    
    //ANA_accessor ana(db);
    //ana.load_up_tables();
    
    //Average_PMIP_ensemble pmip = Average_PMIP_ensemble(db, input_values.time_period_id); 
    //pmip.add_to_database();

    db.close_connection();
} 

// Find the analysis for a given run 
void get_analysis_for_given_area(Assim_input_parser input_values) {
    //if(input_values.open_blas_threading)
    //    openblas_set_num_threads(input_values.number_of_threads);

    printf("Starting assimilator\n");

    sql_database db = sql_database(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address, input_values.tunneling, input_values.tunneling ? input_values.connection_port : input_values.sql_port, input_values.sql_verbose);
    CRU_accessor cru(db);
    OBS_accessor obs(db);
    ANA_accessor ana(db);
    
    ANA_experiment ana_exp(db, input_values.experiment_name, input_values.time_period_id, input_values.rh_scale, input_values.temp_scale, input_values.grid_scale, input_values.grid_degrade);
    
    Average_PMIP_ensemble pmip = Average_PMIP_ensemble(db, input_values.time_period_id); 

    //Grab all the data and parse it into the form that's needed for the minimiser
    pmip_seasonal_assimilator p = pmip_seasonal_assimilator(pmip, obs, cru, 
            input_values.grid_degrade, 
            input_values.lat_max,
            input_values.lat_min, 
            input_values.lon_max,
            input_values.lon_min,
            input_values.open_blas_threading ? 1 : input_values.number_of_threads,
            input_values.time_period_id,
            input_values.use_CO2,
            input_values.exclude_PSM_obs
            );
    db.close_connection();

    //Minimise the cost function and find the analysis and as well as other useful analyitics
    Seasonal_correlated_minimiser::os_map analysis = p.find_analysis(input_values.rh_scale, input_values.temp_scale, input_values.grid_scale);

    //Save the analysis as well as the analytics
    db.open_connection(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address);
    p.save_analysis(ana_exp, analysis);

    printf("Ending assimilator\n");
}


// Find and save the analysis for several different values of L_t
// Several flags need to be changed for this to work:
//  The resolution matrix in climate_data_reader.h needs to be state_space_size x state_space_size
//  The result file save is need
//  Compute resolution needs to be on
void run_multiple_runs_L_t(Assim_input_parser input_values) {
    printf("Starting assimilator\n");

    sql_database db = sql_database(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address, input_values.tunneling, input_values.tunneling ? input_values.connection_port : input_values.sql_port, input_values.sql_verbose);
    CRU_accessor cru(db);
    OBS_accessor obs(db);
    ANA_accessor ana(db);
    Average_PMIP_ensemble pmip = Average_PMIP_ensemble(db, input_values.time_period_id); 

    //Grab all the data and parse it into the form that's needed for the minimiser
    pmip_seasonal_assimilator p = pmip_seasonal_assimilator(pmip, obs, cru,
            input_values.grid_degrade, 
            input_values.lat_max,
            input_values.lat_min, 
            input_values.lon_max,
            input_values.lon_min,
            input_values.open_blas_threading ? 1 : input_values.number_of_threads,
            input_values.time_period_id,
            input_values.use_CO2,
            input_values.exclude_PSM_obs
            );

    db.open_connection(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address);

    vector<double> L_ts = {0.1, 1, 3};
    for(auto L_t : L_ts) {
        printf("Running %f\n", L_t);
        ANA_experiment ana_exp(db, input_values.experiment_name, input_values.time_period_id, input_values.rh_scale, L_t, input_values.grid_scale, input_values.grid_degrade);
        
        //Minimise the cost function and find the analysis and as well as other useful analyitics
        Seasonal_correlated_minimiser::os_map analysis = p.find_analysis(input_values.rh_scale, L_t, input_values.grid_scale);

        //Save the analysis as well as the analytics
        p.save_analysis(ana_exp, analysis);
    }

    db.close_connection();

    printf("Ending assimilator\n");
}

// Find and save the BLUE solution (just one step of the minimiser) for several different values of L_s
void run_multiple_runs_L_s(Assim_input_parser input_values) {
    printf("Starting assimilator\n");

    sql_database db = sql_database(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address, input_values.tunneling, input_values.tunneling ? input_values.connection_port : input_values.sql_port, input_values.sql_verbose);
    CRU_accessor cru(db);
    OBS_accessor obs(db);
    ANA_accessor ana(db);
    Average_PMIP_ensemble pmip = Average_PMIP_ensemble(db, input_values.time_period_id); 

    //Grab all the data and parse it into the form that's needed for the minimiser
    pmip_seasonal_assimilator p = pmip_seasonal_assimilator(pmip, obs, cru,
            input_values.grid_degrade, 
            input_values.lat_max,
            input_values.lat_min, 
            input_values.lon_max,
            input_values.lon_min,
            input_values.open_blas_threading ? 1 : input_values.number_of_threads,
            input_values.time_period_id,
            input_values.use_CO2,
            input_values.exclude_PSM_obs
            );

    db.open_connection(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address);

    for(int L_s = 100; L_s < 1000; L_s+=50) {
        printf("Running %i\n", L_s);
        ANA_experiment ana_exp(db, input_values.experiment_name, input_values.time_period_id, input_values.rh_scale, input_values.temp_scale, L_s, input_values.grid_degrade);
        
        //Minimise the cost function and find the analysis and as well as other useful analyitics
        //Seasonal_correlated_minimiser::os_map analysis = p.find_analysis(input_values.rh_scale, input_values.temp_scale, L_s);
        Seasonal_correlated_minimiser::os_map analysis = p.static_analysis(input_values.rh_scale, input_values.temp_scale, L_s);

        //Save the analysis as well as the analytics
        p.save_analysis(ana_exp, analysis);
    }

    db.close_connection();

    printf("Ending assimilator\n");
}

// Finds the resolution matrix for the BLUE 
// TODO: This isn't currently implemented, the resolution matrix calculation is only avaliable for the full minised problem, not the BLUE 
void get_resolution_matrix(Assim_input_parser input_values) {
    sql_database db = sql_database(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address, input_values.tunneling, input_values.tunneling ? input_values.connection_port : input_values.sql_port, input_values.sql_verbose);

    CRU_accessor cru(db);
    OBS_accessor obs(db);
    ANA_accessor ana(db);
    
    ANA_experiment ana_exp(db, input_values.experiment_name, input_values.time_period_id, input_values.rh_scale, input_values.temp_scale, input_values.grid_scale, input_values.grid_degrade);

    Average_PMIP_ensemble pmip = Average_PMIP_ensemble(db, input_values.time_period_id); 

    //Grab all the data and parse it into the form that's needed for the minimiser
    pmip_seasonal_assimilator p = pmip_seasonal_assimilator(pmip, obs, cru,
            input_values.grid_degrade,
            input_values.lat_max,
            input_values.lat_min,
            input_values.lon_max,
            input_values.lon_min,
            input_values.number_of_threads,
            input_values.time_period_id,
            input_values.use_CO2,
            input_values.exclude_PSM_obs
            );
    
    //Create a file to save all the informaiton
    ofstream res_mat_file("res/resolution.mat");
    Mat_f resol = p.get_resolution_matrix(input_values.rh_scale, input_values.temp_scale, input_values.grid_scale);

    resol.print(res_mat_file, false);
    res_mat_file.close(); 
}

void make_CO2_compare_runs(Assim_input_parser input_values) {
    printf("Starting CO2 compare\n");

    sql_database db = sql_database(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address, input_values.tunneling, input_values.tunneling ? input_values.connection_port : input_values.sql_port, input_values.sql_verbose);
    CRU_accessor cru(db);
    OBS_accessor obs(db);
    ANA_accessor ana(db);
    
    Average_PMIP_ensemble pmip = Average_PMIP_ensemble(db, input_values.time_period_id); 

    db.open_connection(input_values.sql_username, input_values.sql_password, input_values.sql_port, input_values.tunnel_address);
    
    //Perform tests/////////////////
    string exp_name = "WO CO2";

    ANA_experiment ana_exp(db, exp_name, input_values.time_period_id, input_values.rh_scale, input_values.temp_scale, input_values.grid_scale, input_values.grid_degrade);

    //Grab all the data and parse it into the form that's needed for the minimiser
    pmip_seasonal_assimilator p = pmip_seasonal_assimilator(pmip, obs, cru, 
            input_values.grid_degrade, 
            input_values.lat_max,
            input_values.lat_min, 
            input_values.lon_max,
            input_values.lon_min,
            input_values.open_blas_threading ? 1 : input_values.number_of_threads,
            input_values.time_period_id,
            0,
            input_values.exclude_PSM_obs
            );

    //Minimise the cost function and find the analysis and as well as other useful analytics
    Seasonal_correlated_minimiser::os_map analysis = p.find_analysis(input_values.rh_scale, input_values.temp_scale, input_values.grid_scale);

    //Save the analysis as well as the analytics
    p.save_analysis(ana_exp, analysis);

    ///////////////////////////////////
    
    exp_name = "With CO2";

    ANA_experiment ana_exp2(db, exp_name, input_values.time_period_id, input_values.rh_scale, input_values.temp_scale, input_values.grid_scale, input_values.grid_degrade);

    //Grab all the data and parse it into the form that's needed for the minimiser
    pmip_seasonal_assimilator p2(pmip, obs, cru, 
            input_values.grid_degrade, 
            input_values.lat_max,
            input_values.lat_min, 
            input_values.lon_max,
            input_values.lon_min,
            input_values.open_blas_threading ? 1 : input_values.number_of_threads,
            input_values.time_period_id,
            1,
            input_values.exclude_PSM_obs
            );

    //Minimise the cost function and find the analysis and as well as other useful analytics
    analysis = p2.find_analysis(input_values.rh_scale, input_values.temp_scale, input_values.grid_scale);

    //Save the analysis as well as the analytics
    p2.save_analysis(ana_exp2, analysis);

    //End tests/////////////////////
    db.close_connection();

    printf("Ending compare CO2\n");
}

int main(int argc, char * argv[]) {
    //Parse all the configuration options from the config file and command line    
    Assim_input_parser input_values;

    if(argc>1)
        for(int i=1; i<argc; i++)
            input_values.parse_line(argv[i]);

    // Use these functions to select with operation you would like 
    
    //load_up_tables(input_values);
    get_analysis_for_given_area(input_values);
    //run_multiple_runs_L_t(input_values);
    //run_multiple_runs_L_s(input_values);
    //get_resolution_matrix(input_values);
    //make_CO2_compare_runs(input_values)
    
    return 0;
}
