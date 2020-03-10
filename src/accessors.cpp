#include "accessors.h"
/////////////////////////////////////////////////////////////////////////////////////
//INFO///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
void INFO_accessor::load_up_tables() {
    //Clear tables
    time_table.delete_table();
    time_table.initialise();

    //Fill up the time period table of the different time periods
    time_table.add_record("LGM");
    time_table.add_record("MH");
   
    position_table.delete_table();
    position_table.initialise();
}

/////////////////////////////////////////////////////////////////////////////////////
//PMIP///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

vector<vector<netcdf_PMIP*>> PMIP_accessor::register_avaliable_pmip_files() {
    //C++17 has a std::filesystem api but most of this code was designed with C++11 in mind so we'll stick to that. This means this will only really work on Linux
    string base_data_path = "data_input/PMIP/2_deg_gridded";

    //Look in the directory, find all the nc files and hold thier names
    DIR *dirp = opendir(base_data_path.c_str());
    class dirent *dp;

    vector<PMIP_nc_file_meta> nc_file_metas = vector<PMIP_nc_file_meta>();

    while ((dp = readdir(dirp)) != NULL) {
        string full_name = dp->d_name;
        if(full_name.size() < 3)
            continue;

        nc_file_metas.push_back(PMIP_nc_file_meta(full_name, base_data_path));
    }

    (void)closedir(dirp);
   
    //Go through all the files picking out the ones which are just different variables of the same model run
    //At the end we have all the groupings of the of each different run in 'model_runs' 
    vector<PMIP_nc_file_meta>::iterator master_meta = nc_file_metas.begin();
    vector<PMIP_nc_file_meta> meta_grouping = vector<PMIP_nc_file_meta>();
    vector<netcdf_PMIP*> model_runs = vector<netcdf_PMIP*>();
    
    //Get the first meta in the list, see the experiment matches any others in the list and group them together 
    while(master_meta != nc_file_metas.end()) {
        meta_grouping.push_back(*master_meta);

        //Go through th rest of the files, if they match then store it and remove it from the list, if not then move on
        for(vector<PMIP_nc_file_meta>::iterator meta = nc_file_metas.begin(); meta != nc_file_metas.end();) {
            if(master_meta != meta && meta->experiment_model_match(*master_meta)) {
                meta_grouping.push_back(*meta);
                meta = nc_file_metas.erase(meta);
            }
            else {
                meta++;
            }
        }
       
        model_runs.push_back(new netcdf_PMIP(meta_grouping));
        meta_grouping.clear();

        master_meta = nc_file_metas.erase(master_meta);
    }

    vector<vector<netcdf_PMIP*>> ordered_model_runs = vector<vector<netcdf_PMIP*>>();   
    
    vector<netcdf_PMIP*>::iterator master_run = model_runs.begin();
    netcdf_PMIP* lgm_experiment = NULL;
    netcdf_PMIP* mh_experiment = NULL;
    netcdf_PMIP* pi_experiment = NULL;

    while(master_run != model_runs.end()) {
        if((*master_run)->experiment_name == "lgm")
            lgm_experiment = *master_run;
        else if((*master_run)->experiment_name == "midHolocene")
            mh_experiment = *master_run;
        else if((*master_run)->experiment_name == "piControl")
            pi_experiment = *master_run;

        for(vector<netcdf_PMIP*>::iterator run = model_runs.begin(); run != model_runs.end();) {
            if(master_run != run && (*run)->model_name == (*master_run)->model_name) {
                if((*run)->experiment_name == "lgm")
                    lgm_experiment = *run;
                else if((*run)->experiment_name == "midHolocene")
                    mh_experiment = *run;
                else if((*run)->experiment_name == "piControl")
                    pi_experiment = *run;

                run = model_runs.erase(run);
            }
            else {
                run++;
            }
        }
       
        if(pi_experiment) {
            //We discount CCSM4 R2 since it isn't a relevent experiment, Also I thint it's given in Kelvin which just makes thins awkward
            if((*master_run)->model_name != "CCSM4_R2I1P1" )
                //&& (*master_run)->model_name != "MPI_ESM_P_R1I1P1"
                {
                    ordered_model_runs.push_back(vector<netcdf_PMIP*>({
                                lgm_experiment, mh_experiment, pi_experiment
                                }));
                }
        } else {
            string exp_name = "";
            if(lgm_experiment) exp_name = lgm_experiment->model_name;
            if(mh_experiment) exp_name = mh_experiment->model_name;
            printf("PI experiemnt not found for %s\n", exp_name.c_str());
        }

        master_run = model_runs.erase(master_run);
        lgm_experiment = NULL; mh_experiment = NULL; pi_experiment = NULL;
    }

    return ordered_model_runs;
}

void PMIP_accessor::add_to_tables(const vector<vector<netcdf_PMIP*>> &files, CRU_accessor mod) {

    //Find out which ids correspond to the time periods we're considering
    int lgm_time_period_id = time_table.get_id_from_time_period("LGM"); 
    int mh_time_period_id = time_table.get_id_from_time_period("MH"); 

    vector<int> lgm_ids = vector<int>();
    vector<int> mh_ids = vector<int>();

    int cur_id = 1; //We're going to assume the sequential id thing here

    //Load up all the files and add the experiments into the database
    for(auto &model_files : files) {
        if(model_files[0]) {
            model_files[0]->load_all_readings();
            model_table.add_record(model_files[0]->model_name, lgm_time_period_id);
            lgm_ids.push_back(cur_id); cur_id++;
        } 

        if(model_files[1]) {
            model_files[1]->load_all_readings();

            model_table.add_record(model_files[1]->model_name, mh_time_period_id);
            mh_ids.push_back(cur_id); cur_id++;
        } 

        if(model_files[2]) model_files[2]->load_all_readings();
    }

    model_table.add_record("AVG", lgm_time_period_id);
    int lgm_average_experiment_id = model_table.get_last_id();    

    model_table.add_record("AVG", mh_time_period_id);
    int mh_average_experiment_id = model_table.get_last_id();    

    //Go through all the experiments and put them in the table///////////////////////////////////////////////////
    Seasonal_Modern_reading means_lgm = Seasonal_Modern_reading();
    Seasonal_Modern_reading sds_lgm = Seasonal_Modern_reading();
    Seasonal_Modern_reading means_mh = Seasonal_Modern_reading();
    Seasonal_Modern_reading sds_mh = Seasonal_Modern_reading();

    vector<Netcdf_PMIP_Info*> current_readings_lgm = vector<Netcdf_PMIP_Info*>();
    vector<Netcdf_PMIP_Info*> current_readings_mh = vector<Netcdf_PMIP_Info*>();
    
    unsigned int position_id = 1;
    unsigned int scaled_lat_size = 180/default_grid_size;
    unsigned int scaled_lon_size = 360/default_grid_size;

    Vec_f amount_of_experiments_with_data_mh = Vec_f(4); //Should be ints but cba, the mark of true quality
    Vec_f amount_of_experiments_with_data_lgm = Vec_f(4); //Should be ints but cba, the mark of true quality

    Netcdf_PMIP_Info *lgm_reading = NULL;
    Netcdf_PMIP_Info *mh_reading = NULL;
    Netcdf_PMIP_Info *pi_reading = NULL;
    Seasonal_Modern_reading* mod_reading = NULL;

    vector<CRU_accessor::CRU_seasonal_reading*> modern_readings = mod.get_area(90, -90, 180, -180, 1);
    for(unsigned int lat_index = 0; lat_index<scaled_lat_size; lat_index++) { 
        if(LOAD_ALL_TABLES_VERBOSE) printf("Lat %i | ID %i\n", lat_index, position_id);
        
        for(unsigned int lon_index = 0; lon_index<scaled_lon_size; lon_index++) {
            amount_of_experiments_with_data_lgm = Vec_f(4);
            amount_of_experiments_with_data_mh = Vec_f(4);
            means_lgm = Seasonal_Modern_reading();
            sds_lgm = Seasonal_Modern_reading();
            means_mh = Seasonal_Modern_reading();
            sds_mh = Seasonal_Modern_reading();

            pi_reading = files[0][2]->get_lat_lon_record(lat_index, lon_index);
            loc cur_location = loc(); 
 
            for(vector<CRU_accessor::CRU_seasonal_reading*>::iterator mod = modern_readings.begin(); mod != modern_readings.end(); mod++) {

                if((*mod)->get_loc().lat == pi_reading->location.lat && (*mod)->get_loc().lon == pi_reading->location.lon) {

                    cur_location= (*mod)->get_loc();
                    mod_reading = *mod; 
                    //modern_readings.erase(mod);
                    break;
                }
            }
            
            if(mod_reading) {
                //printf("Mod %f\n", mod_reading->get_MAP());
                //Open every reading and add it to the total of means
                //TODO:Restucture the modern part, makes this whole loop very confusing
                for(auto &model : files) {
                    if(!pi_reading) pi_reading = model[2]->get_lat_lon_record(lat_index, lon_index); //If no model here then we need to panic
                    

                    //*pi_reading = ignore_add((*pi_reading)*-1, *mod_reading);
                    //*pi_reading = (*pi_reading)*-1 + *mod_reading;
                
                    //delete mod_reading;

                    if(model[0]) lgm_reading = model[0]->get_lat_lon_record(lat_index, lon_index);
                    if(model[1]) mh_reading = model[1]->get_lat_lon_record(lat_index, lon_index);

                    if(lgm_reading) {
                        *lgm_reading = *lgm_reading + (*pi_reading)*-1 + *mod_reading;//*pi_reading;
                        
                        //Can't have negative MAP so we set it very close to 0 here
                        if(lgm_reading->get_MAP()  < 0)
                            lgm_reading->set_MAP(0.1);
                        
                        //Can't have negative sunf so we set it to 0 here 
                        Vec_f tmp_sunf = lgm_reading->get_sunf();
                        for(int i=0; i<12; i++) {
                            if(tmp_sunf.get(i) < 0)
                                tmp_sunf.set(i, 0);
                            lgm_reading->set_sunf(tmp_sunf);
                        }

                        current_readings_lgm.push_back(lgm_reading);
                        means_lgm = ignore_add(*lgm_reading, means_lgm);

                        amount_of_experiments_with_data_lgm += lgm_reading->values_held;
                    }

                    if(mh_reading) {
                        *mh_reading = *mh_reading + (*pi_reading)*-1 + *mod_reading;

                        //Can't have negative MAP so we set to 0 here
                        if(mh_reading->get_MAP()  < 0)
                            mh_reading->set_MAP(0.1);

                        //Can't have negative sunf so we set it to 0 here 
                        Vec_f tmp_sunf = mh_reading->get_sunf();
                        for(int i=0; i<12; i++) {
                            if(tmp_sunf.get(i) < 0)
                                tmp_sunf.set(i, 0);
                            mh_reading->set_sunf(tmp_sunf);
                        }

                        current_readings_mh.push_back(mh_reading);
                        means_mh = ignore_add(*mh_reading, means_mh);
                        amount_of_experiments_with_data_mh += mh_reading->values_held;
                    }

                    delete pi_reading; pi_reading = NULL;
                    lgm_reading = NULL;
                    mh_reading = NULL;
                }

                mod_reading = NULL;

                means_lgm.section_div(amount_of_experiments_with_data_lgm);
                means_mh.section_div(amount_of_experiments_with_data_mh);

                //Calculate the standard deviation of the model
                for(auto &reading : current_readings_mh) {
                    Netcdf_PMIP_Info temp = (*reading - means_mh);
                    //temp = temp.piecewise_square();
                    for(int i=0; i<37; i++)
                        temp.set(i, temp.get(i)*temp.get(i));
                    sds_mh = ignore_add(temp, sds_mh);
                }

                for(auto &reading : current_readings_lgm) {
                    Netcdf_PMIP_Info temp = (*reading - means_lgm);
                    //temp = temp.piecewise_square();
                    for(int i=0; i<37; i++)
                        temp.set(i, temp.get(i)*temp.get(i));
                    sds_lgm = ignore_add(temp, sds_lgm);
                }
               
                sds_lgm.section_div(amount_of_experiments_with_data_lgm);
                sds_mh.section_div(amount_of_experiments_with_data_mh);

                //sds_lgm.piecewise_root();
                //sds_mh.piecewise_root();
                for(int i=0; i<37; i++)
                    sds_lgm.set(i, sqrt(sds_lgm.get(i)));
                for(int i=0; i<37; i++)
                    sds_mh.set(i, sqrt(sds_mh.get(i)));


                //Add the values and sds to the relevent tables; state table, position table, sd table
                position_table.add_record(cur_location); //All readings should have the same locaiton so we just add the first one here

                sd_table.add_record(position_id, mh_time_period_id, sds_mh);
                sd_table.add_record(position_id, lgm_time_period_id, sds_lgm);
                
                for(unsigned int i=0; i<lgm_ids.size(); i++)
                    state_table.add_record(position_id, lgm_ids[i], *current_readings_lgm[i]); //We assume the model id here and sql starts at 1 so we must add 1

                for(unsigned int i=0; i<mh_ids.size(); i++)
                    state_table.add_record(position_id, mh_ids[i], *current_readings_mh[i]); //We assume the model id here and sql starts at 1 so we must add 1

                //TODO:Mod 1 too high
                state_table.add_record(position_id, lgm_average_experiment_id, means_lgm);  //Add the average model in
                state_table.add_record(position_id, mh_average_experiment_id, means_mh);  //Add the average model in
                position_id++;

                //Clean up memory used
                //for(auto it=current_readings_lgm.begin(); it!=current_readings_lgm.end(); it++)
                for(auto &it : current_readings_lgm)
                    delete it;
                current_readings_lgm.clear();
                for(auto &it : current_readings_mh)
                    delete it;
                current_readings_mh.clear();
            } else {
                if(pi_reading) {
                    delete pi_reading; 
                    pi_reading = NULL;
                }
            }
            
        }
    }

    for(auto &mod: modern_readings) 
        delete mod;
}

//We want to load up all five tables that this accessor controls
//  model table
//  state table
//  sd table
void PMIP_accessor::load_up_tables(CRU_accessor mod) {
    //We delete all the tables and reinitalise them
    //Note that these must be done in order of dependance 
    sd_table.delete_table();
    state_table.delete_table(); 
    model_table.delete_table();
    ice_table.delete_table();
    position_table.delete_table();

    position_table.initialise();
    ice_table.initialise();
    model_table.initialise();
    state_table.initialise();
    sd_table.initialise();

    //Fill model table of the different models that we have
    
    vector<vector<netcdf_PMIP*>> all_netcdf_PMIP = register_avaliable_pmip_files();
    
    add_to_tables(all_netcdf_PMIP, mod);

    //Delete netcdf files
    for(auto &grouping : all_netcdf_PMIP) {
        for(auto &model_exp: grouping) {
            if(model_exp)
                delete model_exp;
        }
    }
    
    //Load in the ice data
    load_in_ice();
}

void PMIP_accessor::load_in_ice() {
    Netcdf_ice lgm_ice("data_input/PMIP/ice_sheet/21k_ice2x2.nc");
    lgm_ice.load_in_info();

    Netcdf_ice mh_ice("data_input/PMIP/ice_sheet/6k_ice2x2.nc");
    mh_ice.load_in_info();

    ice_table.delete_table();
    ice_table.initialise();

    int lgm_time_period_id = time_table.get_id_from_time_period("LGM"); 
    int mh_time_period_id = time_table.get_id_from_time_period("MH"); 

    double cur_lat = 0;
    double cur_lon = 0;
    double cur_ice = 0;

    //Add the LGM ice sheet into the database 
    for(int lat_index=0; lat_index<lgm_ice.get_lat_size(); lat_index++) {
        for(int lon_index=0; lon_index<lgm_ice.get_lon_size(); lon_index++) {
            cur_lat = lgm_ice.get_lat(lat_index);
            cur_lon = lgm_ice.get_lon(lon_index);

            cur_ice = lgm_ice.get_has_ice(lat_index, lon_index);
            int position_id = position_table.get_id_from_lat_lon(cur_lat, cur_lon);
            
            if(position_id>=0) {
                if(LOAD_ALL_TABLES_VERBOSE) printf("Adding LGM ice record: %i %f %f %f\n", position_id, cur_lat, cur_lon, cur_ice);
                ice_table.add_record(lgm_time_period_id, position_id, cur_ice);
            }

            //We really want to do something like this however there may be places where there is no such position id, hence we have to ask for the id, then check that it exists and only insert if it does
            //"INSERT INTO PMIP_ice(position_id, time_period_id, ice_percentage) VALUES ((SELECT INFO_position.id FROM INFO_position WHERE INFO_position.lat=+"cur_lat+" AND "+"INFO_position.lon="+cur_lon"),"+lgm_time_period_id+","+cur_ice+");"


        }
    }
    
    //Add the MH ice sheet into the database
    for(int lat_index=0; lat_index<mh_ice.get_lat_size(); lat_index++) {
        for(int lon_index=0; lon_index<mh_ice.get_lon_size(); lon_index++) {
            cur_lat = mh_ice.get_lat(lat_index);
            cur_lon = mh_ice.get_lon(lon_index);

            cur_ice = mh_ice.get_has_ice(lat_index, lon_index);
            int position_id = position_table.get_id_from_lat_lon(cur_lat, cur_lon);
            
            if(position_id>=0) {
                if(LOAD_ALL_TABLES_VERBOSE) printf("Adding MH ice record: %i %f %f %f\n", position_id, cur_lat, cur_lon, cur_ice);
                ice_table.add_record(mh_time_period_id, position_id, cur_ice);
            }
        }
    }
}

PMIP_accessor::PMIP_seasonal_reading* PMIP_accessor::get_state_on_id(unsigned int id) {
    sql::ResultSet* found_state = state_table.get_from_id(id);
    PMIP_accessor::PMIP_seasonal_reading* output_reading = NULL; 

    if(found_state) {
        found_state->next();
        output_reading = get_state_from_state_sql(found_state);
        delete found_state;
    }

    return output_reading;
}

PMIP_accessor::PMIP_seasonal_reading* PMIP_accessor::get_state_from_model_id(unsigned int position_id, unsigned int model_id) {
    sql::ResultSet* found_state = state_table.get_on_condition(string(PMIP_state_table::position_id_name) + "=" + to_string(position_id) + " AND " + PMIP_state_table::model_id_name + "=" + to_string(model_id));

    PMIP_accessor::PMIP_seasonal_reading* output_reading = NULL; 

    if(found_state) {
        found_state->next();
        output_reading = get_state_from_state_sql(found_state);
        delete found_state;
    }

    output_reading->set_position_id(position_id);
    return output_reading;
}

vector<PMIP_accessor::PMIP_seasonal_reading*> PMIP_accessor::get_state_from_model_id(vector<unsigned int> &position_id, unsigned int model_id) {
    sql::ResultSet* found_model = model_table.get_on_condition(string(PMIP_model_table::basic_id_name )+ "=" + to_string(model_id));
    found_model->next(); found_model->next(); found_model->next();

    unsigned int time_period_id = found_model->getInt(1);
    delete found_model;

    //Create a string with all the positions and the model id in it and then get the relevent states
    string position_condition_string = string(PMIP_state_table::model_id_name) + "=" + to_string(model_id) + " AND (";

    for(vector<unsigned int>::iterator it = position_id.begin(); it!=position_id.end(); it++) {
        position_condition_string += string(PMIP_state_table::position_id_name) + "=" + to_string(*it); 
        if(it != position_id.end()-1) 
            position_condition_string += " OR ";
    }

    position_condition_string += ")"; 

    sql::ResultSet* found_state = state_table.get_on_condition(position_condition_string);

    //Parse the states into the output
    vector<PMIP_accessor::PMIP_seasonal_reading*> output_readings = vector<PMIP_accessor::PMIP_seasonal_reading*>(); 

    //Here we get all the state information together into the proper class (PMIP_seasonal_reading)
    //TODO:This gets the sds from the tables 1 by 1 rather than all at once that this method (with the multiple pos ids) is meant to prevent
    if(found_state) {
        int i=0;
        while(found_state->next()) {
            output_readings.push_back(get_state_from_state_sql(found_state, time_period_id));
            output_readings[i]->set_position_id(position_id[i]);
            i++;
        }
    }

    delete found_state;

    return output_readings;
}

PMIP_accessor::PMIP_seasonal_reading* PMIP_accessor::get_state_from_position_id_name(unsigned int position_id, string model_name, int time_period_id) {
    sql::ResultSet* found_model = model_table.get_on_condition(string(PMIP_model_table::model_name_name )+ "=\'" + model_name + "\' AND " + PMIP_model_table::time_period_id_name+ "=\'" + to_string(time_period_id) + "\'");
    found_model->next();
    unsigned int model_id = found_model->getInt(1);
    delete found_model;

    sql::ResultSet* found_state = state_table.get_on_condition(string(PMIP_state_table::position_id_name) + "=" + to_string(position_id) + " AND " + PMIP_state_table::model_id_name + "=" + to_string(model_id));

    PMIP_accessor::PMIP_seasonal_reading* output_reading = NULL; 

    if(found_state) {
        found_state->next();
        output_reading = get_state_from_state_sql(found_state);
        delete found_state;
    }

    output_reading->set_position_id(position_id);
    return output_reading;
}

//Given a vector of positions ids will find all the states associated with said positions for a given model and time period
vector<PMIP_accessor::PMIP_seasonal_reading*> PMIP_accessor::get_state_from_position_id_name(vector<unsigned int> &position_id, string model_name, int time_period_id) {
    //Figure out what the model id number is that we're looking at
    sql::ResultSet* found_model = model_table.get_on_condition(string(PMIP_model_table::model_name_name )+ "=\'" + model_name + "\' AND " + PMIP_model_table::time_period_id_name+ "=\'" + to_string(time_period_id) + "\'");
    found_model->next();
    unsigned int model_id = found_model->getInt(1);
    delete found_model;

    //Create a string with all the positions and the model id in it and then get the relevent states
    string position_condition_string = string(PMIP_state_table::model_id_name) + "=" + to_string(model_id) + " AND (";

    for(vector<unsigned int>::iterator it = position_id.begin(); it!=position_id.end(); it++) {
       position_condition_string += string(PMIP_state_table::position_id_name) + "=" + to_string(*it); 
        if(it != position_id.end()-1) 
            position_condition_string += " OR ";
    }

    position_condition_string += ")"; 

    sql::ResultSet* found_state = state_table.get_on_condition(position_condition_string);

    //Parse the states into the output
    vector<PMIP_accessor::PMIP_seasonal_reading*> output_readings = vector<PMIP_accessor::PMIP_seasonal_reading*>(); 

    //Here we get all the state information together into the proper class (PMIP_seasonal_reading)
    //TODO:This gets the sds from the tables 1 by 1 rather than all at once that this method (with the multiple pos ids) is meant to prevent
    if(found_state) {
        int i=0;
        while(found_state->next()) {
            output_readings.push_back(get_state_from_state_sql(found_state, time_period_id));
            output_readings[i]->set_position_id(position_id[i]);
            i++;
        }
    }

    delete found_state;

    return output_readings;
}

vector<PMIP_accessor::PMIP_seasonal_reading*> PMIP_accessor::get_state_from_position_statement(string position_string, unsigned int model_id, string time_period_name) {

    //int time_period_id = time_table.get_id_from_time_period(time_period_name); 

    //Create a string with all the positions and the model id in it and then get the relevent states
    string position_condition_string = string(PMIP_state_table::model_id_name) + "=" + to_string(model_id) + " AND position_id in (" + position_string + ")";

    sql::ResultSet* found_state = state_table.run_command("SELECT PMIP_state.* FROM PMIP_state WHERE " + position_condition_string + " INNER JOIN INFO_position.lat, INFO_position.lon, INFO_position.elev on PMIP_state.position_id == INFO_position.id");

    //Parse the states into the output
    vector<PMIP_accessor::PMIP_seasonal_reading*> output_readings = vector<PMIP_accessor::PMIP_seasonal_reading*>(); 

    //Here we get all the state information together into the proper class (PMIP_seasonal_reading)
    //TODO:This gets the sds from the tables 1 by 1 rather than all at once that this method (with the multiple pos ids) is meant to prevent
    if(found_state) {
        while(found_state->next()) {
            /*loc location;
            unsigned int id = found_state->getInt(1);
            unsigned int position_id = found_state->getInt(2);
            unsigned int model_id = found_state->getInt(3);

            location = loc(
                found_position->getDouble(2),
                found_position->getDouble(3),
                found_position->getDouble(4)
                );

            output_reading = new PMIP_accessor::PMIP_seasonal_reading(id, location, position_id);
            for(int i=4; i<state_table.modern_state_variable_amount+4; i++) { 
                output_reading->set(i-4, found_state->getDouble(i));
            }
            output_readings.push_back(get_state_from_state_sql(found_state, time_period_id));*/
        }
    }

    delete found_state;

    return output_readings;

    /*PMIP_accessor::PMIP_seasonal_reading* output_reading = NULL; 
    loc location;
    unsigned int id = found_state->getInt(1);
    unsigned int position_id = found_state->getInt(2);
    unsigned int model_id = found_state->getInt(3);

    if(time_period_id <0) {
        sql::ResultSet* rs = model_table.get_from_id(model_id);
        time_period_id = rs->getInt(3);
    }

    sql::ResultSet* found_position = position_table.get_from_id(position_id);

    if(found_position->next()) {
        location = loc(
                found_position->getDouble(2),
                found_position->getDouble(3),
                found_position->getDouble(4)
                );

        if(found_sds->next()) {
            output_reading = new PMIP_accessor::PMIP_seasonal_reading(id, location, position_id);

            for(int i=4; i<state_table.modern_state_variable_amount+4; i++) { 
                output_reading->set(i-4, found_state->getDouble(i));
            }
        }

        delete found_position;
    } 

    return output_reading;*/
}

PMIP_accessor::PMIP_seasonal_reading* PMIP_accessor::get_state_from_state_sql(sql::ResultSet* found_state, int time_period_id) {
    PMIP_accessor::PMIP_seasonal_reading* output_reading = NULL; 
    loc location;
    unsigned int id = found_state->getInt(1);
    unsigned int position_id = found_state->getInt(2);
    unsigned int model_id = found_state->getInt(3);

    if(time_period_id <0) {
        sql::ResultSet* rs = model_table.get_from_id(model_id);
        time_period_id = rs->getInt(3);
    }

    sql::ResultSet* found_position = position_table.get_from_id(position_id);

    string sd_condition = string(PMIP_sd_table::position_id_name) + "=\'" + to_string(position_id) + "\' AND " + PMIP_sd_table::time_period_id_name+ "=\'" + to_string(time_period_id) + "\'"; 
    sql::ResultSet* found_sds = sd_table.get_on_condition(sd_condition);

    if(found_position->next()) {
        location = loc(
                found_position->getDouble(2),
                found_position->getDouble(3),
                found_position->getDouble(4)
                );

        if(found_sds->next()) {
            output_reading = new PMIP_accessor::PMIP_seasonal_reading(id, location, position_id);

            for(int i=4; i<state_table.modern_state_variable_amount+4; i++) { 
                output_reading->set(i-4, found_state->getDouble(i));
                output_reading->get_sd().set(i-4, (double)found_sds->getDouble(i));
            }
            delete found_sds;
        }

        delete found_position;
    } 

    return output_reading;
}

/*vector<PMIP_accessor::Ice_reading*> get_ice_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade, int time_period_id) {

    vector<PMIP_accessor::Ice_reading*> full_state_readings = vector<PMIP_accessor::PMIP_seasonal_reading*>(); 
    sql::ResultSet* found_state = state_table.run_command(
            "SELECT PMIP_ice.*"
            "FROM PMIP_state "
            "INNER JOIN INFO_position on PMIP_ice.position_id = INFO_position.id "
            "WHERE ( INFO_position.lat BETWEEN " + to_string(lat_min) + " AND " + to_string(lat_max) + ") "
            "AND (INFO_position.lon BETWEEN " + to_string(lon_min) + " AND " + to_string(lon_max) + ") "
            "AND ((INFO_position.lat - "+to_string(default_grid_size/2)+")%"+to_string(default_grid_size*grid_degrade) + "=0) "
            "AND ((INFO_position.lon - "+to_string(default_grid_size/2)+")%"+to_string(default_grid_size*grid_degrade) + "=0) "
            "AND (PMIP_state." + string(PMIP_state_table::model_id_name) + "=" + to_string(model_id) + ") "
            "AND (PMIP_state_sd.time_period_id="+ to_string(time_period_id) + ")"
            );

    if(found_state) {
        while(found_state->next()) {
            unsigned int id = found_state->getInt(1);
            unsigned int position_id = found_state->getInt(2);

            int state_and_sd_amount = 2*(3+state_table.modern_state_variable_amount);

            loc location;
            location = loc(
                    found_state->getDouble(state_and_sd_amount + 1),
                    found_state->getDouble(state_and_sd_amount + 2),
                    found_state->getDouble(state_and_sd_amount + 3)
                    );

            //TODO:No one is incharge of deleting these
            PMIP_accessor::PMIP_seasonal_reading *output_reading = new PMIP_accessor::PMIP_seasonal_reading(id, location, position_id);
            for(int i=4; i<state_table.modern_state_variable_amount+4; i++) { 
                output_reading->set(i-4, found_state->getDouble(i));
                output_reading->get_sd().set(i-4, found_state->getDouble(state_table.modern_state_variable_amount+3 + i));
            }

            full_state_readings.push_back(output_reading); 
        }

        delete found_state;
    }

    return full_state_readings;

}*/

vector<int> PMIP_accessor::get_all_model_ids_from_period(int time_period_id) {
    sql::ResultSet *rs = state_table.run_command("SELECT id from PMIP_model where time_period_id=" + std::to_string(time_period_id));
    
    vector<int> output = vector<int>();

    if(rs) {
        while(rs->next())
            output.push_back(rs->getInt(1));

        delete rs;
    }

    return output;
}


int PMIP_accessor::find_time_average(int time_period_id) {
    sql::ResultSet *rs = state_table.run_command("SELECT id from PMIP_model where time_period_id=" + std::to_string(time_period_id) + " and model_name=\'AVG\'");
    
    int output = -1;
    if(rs) {
        while(rs->next())
            output = rs->getInt(1);

        delete rs;
    }

    return output;
}

sql::ResultSet* PMIP_accessor::get_all_positions() {
    return position_table.get_all();
}


//////////////////////////////////////////////////////////////////////////////////////////////////
//Single_Model_PMIP_accessor//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
PMIP_accessor::PMIP_seasonal_reading* Single_Model_PMIP_accessor::get_state_from_position(sql::ResultSet* rs) {
    //We assume that the record that is wanted is already loaded in 
    unsigned int pos_id = -1;
    if(rs) {
        pos_id = rs->getInt(1);
    } else {
        return NULL;
    }

    return get_state_from_model_id(pos_id, model_id);
}

vector<PMIP_accessor::PMIP_seasonal_reading*> Single_Model_PMIP_accessor::get_state_from_positions(sql::ResultSet* rs) {
    vector<unsigned int> position_ids = vector<unsigned int>();
    if(rs)
        while(rs->next())
            position_ids.push_back(rs->getInt(1));

    return get_state_from_model_id(position_ids, model_id);
}

PMIP_accessor::PMIP_seasonal_reading* Single_Model_PMIP_accessor::get_closest_state(loc location) {
    double c_lat = location.lat;
    double c_lon = location.lon;

    //Find a box a small box around the target which contain some records
    float run_size = 0;
    sql::ResultSet *res = NULL; 
    
    //Note here the bottle neck is normally the sql query, not the finding the closest record from the sub set (this however could also be improved)
    while(run_size <2) {
        run_size += 0.25;//TODO: This should be somehow tuned to the grid size
        //We create a box of lat and lon to search for, note that lat and lon wray around the earth so we account for this
        double low_lon = fmod(c_lon - run_size + 180, 360) - 180;
        double high_lon = fmod(c_lon + run_size + 180, 360) - 180; 
        double low_lat = fmod(c_lat - run_size + 90, 180) - 90;
        double high_lat = fmod(c_lat + run_size + 90, 180) - 90;

        try {
            res = position_table.get_on_condition(position_table.get_box_condition(low_lat, high_lat, low_lon, high_lon));
        } catch(sql::SQLException* e) {
            printf("Problem with SQL query: %s\n", e->getSQLState().c_str());
        }

        if(res && res->rowsCount() != 0) {
            break;
        }
    }

    //Find the closest neighbour in the box 
    PMIP_accessor::PMIP_seasonal_reading* closest = NULL;
    float closest_distance = 259201; //2*360**2 + 1

    if(res && res->rowsCount() != 0) {
        vector<PMIP_accessor::PMIP_seasonal_reading*> close_readings = get_state_from_positions(res);
        for(auto &reading_we_think_is_close : close_readings) {
            if(reading_we_think_is_close) {
                float cur_dis = reading_we_think_is_close->get_loc().distance_between(location);

                if(closest_distance > cur_dis) {
                    closest_distance = cur_dis;
                    delete closest;
                    closest = reading_we_think_is_close;
                } else {
                    delete reading_we_think_is_close;
                }
                
                //Since these data sets of often on the same grid this check is worth it
                //TODO:This leaves the rest of the readings 
                if(closest_distance == 0) 
                    break;
            }
        }
    }

    delete res;

    return closest;
}

//Set all the state values for the PMIP model output in a given area
vector<PMIP_accessor::PMIP_seasonal_reading*> Single_Model_PMIP_accessor::get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade) {

    vector<PMIP_accessor::PMIP_seasonal_reading*> full_state_readings = vector<PMIP_accessor::PMIP_seasonal_reading*>(); 
    sql::ResultSet* found_state = state_table.run_command(
            "SELECT PMIP_state.*, PMIP_state_sd.*, INFO_position.lat, INFO_position.lon, INFO_position.elev, PMIP_ice.ice_percentage "
            "FROM PMIP_state "
            "INNER JOIN INFO_position on PMIP_state.position_id = INFO_position.id "
            "INNER JOIN PMIP_state_sd on PMIP_state.position_id = PMIP_state_sd.position_id "
            "INNER JOIN PMIP_ice on PMIP_state.position_id = PMIP_ice.position_id "
            "WHERE ( INFO_position.lat BETWEEN " + to_string(lat_min) + " AND " + to_string(lat_max) + ") "
            "AND (INFO_position.lon BETWEEN " + to_string(lon_min) + " AND " + to_string(lon_max) + ") "
            "AND ((INFO_position.lat - "+to_string(default_grid_size/2)+")%"+to_string(default_grid_size*grid_degrade) + "=0) "
            "AND ((INFO_position.lon - "+to_string(default_grid_size/2)+")%"+to_string(default_grid_size*grid_degrade) + "=0) "
            "AND (PMIP_state." + string(PMIP_state_table::model_id_name) + "=" + to_string(model_id) + ") "
            "AND (PMIP_state_sd.time_period_id="+ to_string(time_period_id) + ") "
            "AND (PMIP_ice.time_period_id="+ to_string(time_period_id) + ") "
            );

    if(found_state) {
        while(found_state->next()) {
            unsigned int id = found_state->getInt(1);
            unsigned int position_id = found_state->getInt(2);

            int state_and_sd_amount = 2*(3+state_table.modern_state_variable_amount);

            loc location;
            location = loc(
                    found_state->getDouble(state_and_sd_amount + 1),
                    found_state->getDouble(state_and_sd_amount + 2),
                    found_state->getDouble(state_and_sd_amount + 3)
                    );

            //TODO:No one is incharge of deleting these
            PMIP_accessor::PMIP_seasonal_reading *output_reading = new PMIP_accessor::PMIP_seasonal_reading(id, location, position_id);
            for(int i=4; i<state_table.modern_state_variable_amount+4; i++) { 
                output_reading->set(i-4, found_state->getDouble(i));
                output_reading->get_sd().set(i-4, found_state->getDouble(state_table.modern_state_variable_amount+3 + i));
            }

            output_reading->set_ice_percentage(found_state->getDouble(state_and_sd_amount + 4));

            full_state_readings.push_back(output_reading); 
        }

        delete found_state;
    }

    return full_state_readings;
}

////////////////////////////////////////////////////////////////////////////
//Average_PMIP_ensemble/////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
Average_PMIP_ensemble::Average_PMIP_ensemble(sql_database &d, int time_per_id) : PMIP_ensemble(d), time_period_id(time_per_id) {
    int time_average_id = find_time_average(time_period_id);

    if(time_average_id>0)
        members.push_back(Single_Model_PMIP_accessor(d, time_average_id, time_period_id));
    else {
        vector<int> useable_ids = get_all_model_ids_from_period(time_period_id);

        for(auto &it : useable_ids) 
            members.push_back(Single_Model_PMIP_accessor(d, it, time_period_id));
    }
}

PMIP_accessor::PMIP_seasonal_reading* Average_PMIP_ensemble::get_closest_state(loc location) {
    if(members.size() > 1) {
        vector<PMIP_accessor::PMIP_seasonal_reading*> every_member_close_readings;
        //You probably only need the position of one of these if everything is on the same grid 
        for(auto& it : members) 
            every_member_close_readings.push_back(it.get_closest_state(location));

        return average_readings(every_member_close_readings); 
    } else if (members.size() == 1) {
        return members[0].get_closest_state(location);
    }

    return NULL;
}

vector<PMIP_accessor::PMIP_seasonal_reading*> Average_PMIP_ensemble::get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade) { 
    
    vector<PMIP_accessor::PMIP_seasonal_reading*> output;
    if(members.size() > 1) {
        //Load together every reading in the area from all the members 
        vector<vector<PMIP_accessor::PMIP_seasonal_reading*>> every_member_area_readings;
        for(auto& it : members) 
            every_member_area_readings.push_back(it.get_area(lat_max, lat_min, lon_max, lon_min, grid_degrade));

        //Go through every grid cell and average all the members
        vector<PMIP_accessor::PMIP_seasonal_reading*> tmp_point_matrix;

        for(unsigned int i=0; i<every_member_area_readings[0].size(); i++) {
            for(auto& it : every_member_area_readings)
               tmp_point_matrix.push_back(it[i]);

            output.push_back(average_readings(tmp_point_matrix));
            tmp_point_matrix.clear();
        }
    } else if(members.size() == 1){
        output = members[0].get_area(lat_max, lat_min, lon_max, lon_min, grid_degrade);
    }

    return output;
}

//Averages all the readings into the first reading
PMIP_accessor::PMIP_seasonal_reading* Average_PMIP_ensemble::average_readings(vector<PMIP_accessor::PMIP_seasonal_reading*> all_readings) { 
    PMIP_accessor::PMIP_seasonal_reading* output = new PMIP_accessor::PMIP_seasonal_reading(all_readings[0]->get_main_id(), all_readings[0]->get_location(), all_readings[0]->get_position_id());
    output->set_sd(all_readings[0]->get_sd());
    double amount = all_readings.size();

    //Add together all the elements to a different reading deleting the elements as we go
    for(auto& it : all_readings) {
        *output += *it;
        delete it;
        it++;
    }

    *output /= amount;
    return output;
}

void Average_PMIP_ensemble::add_to_database() {
    model_table.add_record("AVG", time_period_id);
    int experiment_id = model_table.get_last_id();    

    vector<PMIP_accessor::PMIP_seasonal_reading*> all_avg = get_area(90, -90, 180, -180, 1);
    for(auto &it : all_avg)
        state_table.add_record(it->get_position_id(), experiment_id, *it);
}

////////////////////////////////////////////////////////////////////////////////////////
//PMIP_ice//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
//OBS_accessor//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
//Load a particular bartlein file into the database
void OBS_accessor::load_in_bart(bartlein_csv_file &bartlein_data, int time_period_id, CRU_accessor &modern_values) {
    
    //Start going through the bartlein file
    bartlein_data.start(); //bartlein_data.move_next();
    bartlein_csv_file::bart_format_return *current_record_bart = NULL; 

    //Set up the values to save to the table
    Vec_f full_values = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
    Vec_f full_sds = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
    Seasonal_Modern_reading *modern_state = NULL; 

    do {
        //Go through every reading
        bartlein_data.move_next();
        current_record_bart = bartlein_data.get_next_reading();

        //Reset the variables that will be saved to the database
        full_values = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
        full_sds = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
   
        //Sometimes elevation values aren't given by bart so we set them as 0 
        if(current_record_bart->location.elev != current_record_bart->location.elev)
            current_record_bart->location.elev = 0;

        //Find the modern state and non dimensionalise the P and T values so that they can be used with sto_no_CO2
        modern_state = modern_values.get_closest(current_record_bart->location, 5);
        if(!modern_state) {
            printf("Modern state not found for lat, lon: %f, %f\n", current_record_bart->location.lat, current_record_bart->location.lon);  
            continue;
        }

        Seasonal_space_reading space(modern_state->get_sub_vec(0,13));
        MI_alpha_calc::nondimensionalise_state(space);
        modern_state->set(0, space);

        MI_alpha_calc mi_calc(FIF_orbital, current_record_bart->location);
        am_space_reading modern_obs = mi_calc.sto_no_CO2(*modern_state);
        modern_obs = MI_alpha_calc::dimensionalise_obs(modern_obs);
        delete modern_state;

        double alpha = modern_obs.get(0) + current_record_bart->values.get(0);
        if(alpha<0) {
            printf("Error: Alpha less than 0 observed at %f, %f\n", current_record_bart->location.lat, current_record_bart->location.lon);
            full_values.set(0, NAN); 
            full_sds.set(0, NAN);

        } else {
            full_values.set(0, alpha);
            full_sds.set(0, current_record_bart->sds.get(0));
        }

        //The Bartlein file has alpha values so we have to insert a null reading on the moisture value in the table 
        full_values.set(1, NAN);
        full_sds.set(1, NAN);
        
        double value= 0;
        for(int i=2; i<OBS_obs_table::alpha_moisture_variable_amount; i++) {
            value = modern_obs.get(i) + current_record_bart->values.get(i-1);

            if(value<0 and (i==2 or i==6)) {
                printf("Error: Variable %i less than 0 observed at %f, %f with value %f\n", i, current_record_bart->location.lat, current_record_bart->location.lon, value);
                //TODO:Shouldn't these be i not 0
                full_values.set(i, NAN); 
                full_sds.set(i, NAN);
            } else {
                full_values.set(i, value);
                full_sds.set(i, current_record_bart->sds.get(i-1));
            }
        }


        for(int i=0; i<OBS_obs_table::alpha_moisture_variable_amount; i++) {
            if(full_sds.get(i) <= 0) {
                full_values.set(i, NAN);
                full_sds.set(i, NAN);
            }
        }

        obs.add_record(time_period_id, current_record_bart->location, full_values, full_sds, current_record_bart->co2_corrected);
        delete current_record_bart;
    } while(!bartlein_data.on_last_line());
}

//Load in all the known observations. Both Bartlein time periods and also the Australian info
void OBS_accessor::load_up_tables(CRU_accessor &modern_values) {
    INFO_time_table time_table = INFO_time_table(db); 
    int lgm_time_period_id = time_table.get_id_from_time_period("LGM"); 
    int mh_time_period_id = time_table.get_id_from_time_period("MH"); 
    
    //Flush the table of any previous observations
    obs.reset();

    //Open files for bartlein and load them in
    bartlein_csv_file lgm_bartlein_data("data_input/bartlein_raw/site_anomolies/bartlein_raw_21k.csv");
    bartlein_csv_file mh_bartlein_data("data_input/bartlein_raw/site_anomolies/bartlein_raw_06k.csv");

    load_in_bart(lgm_bartlein_data, lgm_time_period_id, modern_values);
    load_in_bart(mh_bartlein_data, mh_time_period_id, modern_values);

    //Set up vectors to hold Australia data 
    Vec_f full_values = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
    Vec_f full_sds = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
    
    //Load up Australia file and put it into the database
    aus_csv_file australia_data("data_input/aus/aus_21k_data_updated.csv");
    australia_data.start();
    aus_csv_file::aus_format_return *current_record_aus = NULL;
    do{
        australia_data.move_next();
        current_record_aus = australia_data.get_next_reading();
        full_values = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);
        full_sds = Vec_f(OBS_obs_table::alpha_moisture_variable_amount);

        //The australia file has moisture values so we have to insert a null reading on the alpha value in the table 
        full_values.set(0, NAN);
        full_sds.set(0, NAN);
        
        for(int i=0; i<OBS_obs_table::alpha_moisture_variable_amount-1; i++) {
            full_values.set(i+1, current_record_aus->values.get(i));
            full_sds.set(i+1, current_record_aus->sds.get(i));
        }

        obs.add_record(lgm_time_period_id, current_record_aus->location, full_values, full_sds);
        
        delete current_record_aus;

    } while(!australia_data.on_last_line());
}

OBS_accessor::OBS_obs_reading* OBS_accessor::get_from_sql_result(sql::ResultSet* found_obs) {
    OBS_accessor::OBS_obs_reading* output_reading = NULL; 

    if(found_obs) {
        unsigned int id = found_obs->getInt(1);

        loc location = loc(
                found_obs->getDouble(3),
                found_obs->getDouble(4),
                found_obs->getDouble(5)
                );
        output_reading = new OBS_accessor::OBS_obs_reading(id, location, found_obs->getBoolean(6));

        for(int i=0; i<7; i++) {
            output_reading->set(i, found_obs->getDouble(i+7));
            output_reading->get_sd().set(i, found_obs->getDouble(i+14));
        }
    } 

    return output_reading;
}

OBS_accessor::OBS_obs_reading* OBS_accessor::get_obs_on_id(unsigned int id) {
    sql::ResultSet* found_obs = obs.get_from_id(id);
    found_obs->next();
    OBS_accessor::OBS_obs_reading* output_reading = get_from_sql_result(found_obs);
    delete found_obs;
    return output_reading;
}

OBS_accessor::OBS_obs_reading* OBS_accessor::get_closest(loc location) {
    OBS_accessor::OBS_obs_reading* output_full_reading = NULL;

    double c_lat = location.lat;
    double c_lon = location.lon;

    //Find a box a small box around the target which contain some records
    float run_size = 0;
    sql::ResultSet *res = NULL; 
    
    //Note here the bottle neck is normally the sql query, not the finding the closest record from the sub set (this however could also be improved)
    while(run_size <2) {
        run_size += 0.25;//TODO: This should be somehow tuned to the grid size
        //We create a box of lat and lon to search for, note that lat and lon wray around the earth so we account for this
        double low_lon = fmod(c_lon - run_size + 180, 360) - 180;
        double high_lon = fmod(c_lon + run_size + 180, 360) - 180; 
        double low_lat = fmod(c_lat - run_size + 90, 180) - 90;
        double high_lat = fmod(c_lat + run_size + 90, 180) - 90;
        try {
            res = obs.get_on_condition(obs.get_box_condition(low_lat, high_lat, low_lon, high_lon));
        } catch(sql::SQLException* e) {
            printf("Problem with SQL query: %s\n", e->getSQLState().c_str());
        }

        if(res && res->rowsCount() != 0) {
            //Find the closest neighbour in the box 
            float closest_distance = 259201; //2*360**2 + 1

            while(res->next()) {
                OBS_accessor::OBS_obs_reading* reading_we_think_is_close = get_from_sql_result(res);

                if(reading_we_think_is_close) {
                    float cur_dis = reading_we_think_is_close->get_loc().distance_between(location);

                    if(closest_distance > cur_dis) {
                        closest_distance = cur_dis;
                        delete output_full_reading;
                        output_full_reading = reading_we_think_is_close;
                    } else {
                        delete reading_we_think_is_close;
                    }

                    //Since these data sets of often on the same grid this check is worth it
                    if(closest_distance == 0) 
                        break;
                }
            }
            delete res;

            break;
        } else {
            delete res; //TODO:Copy this to the other closest methods
        }
    }

    return output_full_reading;
}

vector<OBS_accessor::OBS_obs_reading*> OBS_accessor::get_area(double lat_max, double lat_min, double lon_max, double lon_min, int time_period_id) {

    vector<OBS_accessor::OBS_obs_reading*> output = vector<OBS_accessor::OBS_obs_reading*>();

    sql::ResultSet *res = NULL; 
    
    try {
        res = obs.get_on_condition(obs.get_box_condition(lat_min, lat_max, lon_min, lon_max) + " AND time_period_id=" + to_string(time_period_id));
    } catch(sql::SQLException* e) {
        printf("Problem with SQL query: %s\n", e->getSQLState().c_str());
    }

    if(res && res->rowsCount() != 0) {
        while(res->next()) {
            output.push_back(get_from_sql_result(res));
        }
    }

    delete res;

    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////
//CRU_accessor///////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
void CRU_accessor::load_up_tables() {
    //Make the cru file handler with all the files we're interested in
    cru_file raw_cru(
            "data_input/CRU_CL_v2/tmp_10min_resort.dat",
            "data_input/CRU_CL_v2/pre_10min_resort.dat",
            "data_input/CRU_CL_v2/sunp_10min_resort.dat",
            "data_input/CRU_CL_v2/reh_10min_resort.dat",
            "data_input/CRU_CL_v2/elv_10min_resort.dat"
            );

    //Restart state table
    state.delete_table();
    state.initialise(); 

    //Start the CRU file accessor
    raw_cru.start(); raw_cru.move_next();

    //Go through every line in the CRU file, convert the line into a modern reading and add it to the table
    int current_id = 1;
    do {
        if(LOAD_ALL_TABLES_VERBOSE) printf("Current ID %i\n", current_id);

        cru_file::cru_format_return* current = raw_cru.get_next_reading();

        //Add the reading to the state table
        state.add_record(current->location, current->values);

        current_id++;
    } while(raw_cru.move_next());
}

CRU_accessor::CRU_seasonal_reading* CRU_accessor::get_state_from_id(unsigned int id) {
    sql::ResultSet* found_state = state.get_from_id(id);
    found_state->next();
    CRU_accessor::CRU_seasonal_reading* output_reading = get_state_from_sql_result(found_state);
    
    delete found_state;
    return output_reading;
}

CRU_accessor::CRU_seasonal_reading* CRU_accessor::get_state_from_sql_result(sql::ResultSet* found_state) {
    CRU_accessor::CRU_seasonal_reading* output_reading = NULL;

    if(found_state){
        loc location = loc(
                found_state->getDouble(2),
                found_state->getDouble(3),
                found_state->getDouble(4)
                );

        output_reading = new CRU_accessor::CRU_seasonal_reading(found_state->getInt(1), location); 

        for(int i=0; i<state.modern_state_variable_amount; i++) {
            output_reading->set(i, found_state->getDouble(5+i));
        }
    }

    return output_reading;
}

CRU_accessor::CRU_seasonal_reading* CRU_accessor::get(loc location) {
    CRU_accessor::CRU_seasonal_reading* output_full_reading = NULL; 

    double c_lat = location.lat;
    double c_lon = location.lon;

    sql::ResultSet *res = NULL; 
    
    try {
        res = state.get_on_condition(string("lat=")+ to_string(c_lat) + " AND " + "lon="+to_string(c_lon));
    } catch(sql::SQLException* e) {
        printf("Problem with SQL query: %s\n", e->getSQLState().c_str());
    }

    if(res && res->rowsCount() != 0) {
        res->next();
        output_full_reading = get_state_from_sql_result(res);
    } 

    delete res;

    return output_full_reading;
}

CRU_accessor::CRU_seasonal_reading* CRU_accessor::get_closest(loc location, double max_grid_radius) {
    CRU_accessor::CRU_seasonal_reading* output_full_reading = NULL; 

    double c_lat = location.lat;
    double c_lon = location.lon;

    //Find a box a small box around the target which contain some records
    float run_size = 0;
    sql::ResultSet *res = NULL; 
    
    //Note here the bottle neck is normally the sql query, not the finding the closest record from the sub set (this however could also be improved)
    while(run_size < max_grid_radius) {
        run_size += 0.2;//TODO: This should be somehow tuned to the grid size
        //We create a box of lat and lon to search for, note that lat and lon wray around the earth so we account for this
        double low_lon = fmod(c_lon - run_size + 180, 360) - 180;
        double high_lon = fmod(c_lon + run_size + 180, 360) - 180; 
        double low_lat = fmod(c_lat - run_size + 90, 180) - 90;
        double high_lat = fmod(c_lat + run_size + 90, 180) - 90;
        try {
            res = state.get_on_condition(state.get_box_condition(low_lat, high_lat, low_lon, high_lon));
        } catch(sql::SQLException* e) {
            printf("Problem with SQL query: %s\n", e->getSQLState().c_str());
        }

        if(res && res->rowsCount() != 0) {
            //Find the closest neighbour in the box 
            float closest_distance = 259201; //2*360**2 + 1
            while(res->next()) {
                CRU_accessor::CRU_seasonal_reading* reading_we_think_is_close = get_state_from_sql_result(res);

                if(reading_we_think_is_close) {
                    float cur_dis = reading_we_think_is_close->get_loc().distance_between(location);

                    if(closest_distance > cur_dis) {
                        closest_distance = cur_dis;
                        delete output_full_reading;
                        output_full_reading = reading_we_think_is_close;
                    } else {
                        delete reading_we_think_is_close;
                    }

                    //Since these data sets of often on the same grid this check is worth it
                    if(closest_distance == 0) 
                        break;
                }
            }

            delete res;
            break;

        } 
        delete res;
    }

    return output_full_reading;
}


//Instead of finding the closest here we're going to cheat a bit and snap to the CRU grid
vector<CRU_accessor::CRU_seasonal_reading*> CRU_accessor::get_closest_set(vector<loc> close_locations) {
    vector<CRU_accessor::CRU_seasonal_reading*> outputs = vector<CRU_accessor::CRU_seasonal_reading*>();
    
    double g = default_grid_size/2;//1.0/12.0;
    for(unsigned int i=0; i<close_locations.size(); i++) {
        double lat_snap, lon_snap;
        lat_snap = g*floor(close_locations[i].lat/g + default_grid_size);
        lon_snap = g*floor(close_locations[i].lon/g + default_grid_size);
        sql::ResultSet* res = state.run_command("SELECT * FROM CRU_state WHERE (lat between " + to_string(lat_snap - g/2) + " AND " + to_string(lat_snap + g/2) +") AND (lon between " + to_string(lon_snap - g/2) + " AND " + to_string(lon_snap + g/2) +")");
                
        if(res && res->rowsCount() != 0) {
            while(res->next()) {
                outputs.push_back(get_state_from_sql_result(res));
            }
        }
    }
    return outputs;
}

vector<CRU_accessor::CRU_seasonal_reading*> CRU_accessor::get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade) {
    vector<CRU_accessor::CRU_seasonal_reading*> full_readings = vector<CRU_accessor::CRU_seasonal_reading*>(); 

    sql::ResultSet *res = NULL; 
    
    try {
        res = state.get_on_condition(state.get_box_condition(lat_min, lat_max, lon_min, lon_max) +" and (lat-"+to_string(default_grid_size/2)+")%"+ to_string(default_grid_size*grid_degrade) + "=0 and (lon-"+ to_string(default_grid_size/2)+")%" + to_string(default_grid_size*grid_degrade)+"=0");
    } catch(sql::SQLException* e) {
        printf("Problem with SQL query: %s\n", e->getSQLState().c_str());
    }

    if(res && res->rowsCount() != 0) {
        while(res->next()) {
            full_readings.push_back(get_state_from_sql_result(res));
        }
    }

    delete res;

    return full_readings;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Analysis/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void ANA_accessor::load_up_tables() {
    seasonal_state_obs.delete_table();
    experiments.delete_table();
    experiments.initialise();
    seasonal_state_obs.initialise();
}

void ANA_accessor::add_seasonal_record(int position_id, int experiment_id, Vec_f state, Vec_f obs, Vec_f obs_sd, Vec_f attrib) {
    seasonal_state_obs.add_record(position_id, experiment_id, state, obs, obs_sd, attrib);
}

//TODO:Should experiment just extend accessor?
ANA_experiment::ANA_experiment(sql_database &d, string experi_name, int time_period_id, double rh_scale, double temp_scale, double grid_scale, int grid_degrade) : ANA_accessor(d) {
    //Create experiment, add it to table and find the id
    experiments.add_record(experi_name, time_period_id, rh_scale, temp_scale, grid_scale, grid_degrade);
    experiment_id = experiments.get_last_id();    
}

void ANA_experiment::add_record(int position_id, Vec_f state, Vec_f obs, Vec_f obs_sd,  Vec_f attrib) {
    ANA_accessor::add_seasonal_record(position_id, experiment_id, state, obs, obs_sd, attrib);
}

void ANA_experiment::add_seasonal_record(int position_id, Vec_f seasonal_state, Vec_f obs, Vec_f obs_sd, Vec_f attrib) {
    ANA_accessor::add_seasonal_record(position_id, experiment_id, seasonal_state, obs, obs_sd, attrib);
}
