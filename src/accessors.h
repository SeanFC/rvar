//How to access all the clim data 

#ifndef accessors_h 
#define accessors_h 

#define LOAD_ALL_TABLES_VERBOSE 1

#include <stdlib.h> //for abs
#include <dirent.h>
#include <string>

#include "db_tables.h"
#include "txt_files.h"
#include "readings.h"
#include "geo_observations.h"

#define default_grid_size 2

class CRU_accessor;

///////////////////////////////////////////////////////////////////////////////////////////////////////
//INFO/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
class INFO_accessor {
    private:
        INFO_time_table time_table;
        INFO_position_table position_table;
        sql_database &db;


    public:
        INFO_accessor(sql_database &d) : time_table(d), position_table(d), db(d) {}

        void load_up_tables();
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//PMIP/////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
class PMIP_accessor {
    private:
        void load_in_timed_readings(
                int time_period_id, 
                vector<netcdf_PMIP*> models, 
                vector<int> model_ids, 
                int average_experiment_id, 
                bool add_position=true);

       vector<vector<netcdf_PMIP*>> register_avaliable_pmip_files();
       void add_to_tables(const vector<vector<netcdf_PMIP*>> &files, CRU_accessor mod);

    public:
        enum Experiments_done {
            LGM_EXPERIMENT,
            MH_EXPERIMENT,
            BOTH_EXPERIMENTS
        };

        sql_database &db; 

        INFO_time_table time_table;
        INFO_position_table position_table;

        PMIP_model_table model_table;
        PMIP_state_table state_table;
        PMIP_sd_table sd_table;
        PMIP_ice_table ice_table;

    public:
        PMIP_accessor(sql_database& d) : db(d), time_table(db), position_table(db), model_table(db), state_table(db), sd_table(db), ice_table(db) {}
        void load_up_tables(CRU_accessor);
        void load_in_ice();

        class PMIP_seasonal_reading : public Modern_sd_loc_reading {
            protected:
                unsigned int position_id;
                double ice_percentage;

            public:
                PMIP_seasonal_reading(unsigned int id, loc location, unsigned int pos_id, double ice_per) : Modern_sd_loc_reading(id, location), position_id(pos_id), ice_percentage(ice_per) {};
                PMIP_seasonal_reading(unsigned int id, loc location, unsigned int pos_id) : Modern_sd_loc_reading(id, location), position_id(pos_id), ice_percentage(0) {}; //We're assuming no ice

                void set_position_id(double pos_id) { position_id = pos_id; };
                unsigned int get_position_id() { return position_id; };
                void set_ice_percentage(double ice) { ice_percentage = ice; };
                double get_ice_percentage() { return ice_percentage; };
        };

        /*class Ice_reading : public id_loc {
            protected:
                unsigned int position_id;
                double ice_percentage;

            public:
                Ice_reading(unsigned int id, loc location, unsigned int pos_id, double ice_per) : id_loc(id, location), position_id(pos_id), ice_percentage(ice_per) {};

                void set_position_id(double pos_id) { position_id = pos_id; };
                unsigned int get_position_id() { return position_id; };
        };*/

        void load_in_model(unsigned int number_number);

        PMIP_seasonal_reading* get_state_on_id(unsigned int id);
        PMIP_seasonal_reading* get_state_from_state_sql(sql::ResultSet* found_state, int time_period_id=-1);

        PMIP_seasonal_reading* get_state_from_model_id(unsigned int position_id, unsigned int model_id);
        PMIP_seasonal_reading* get_state_from_position_id_name(unsigned int position_id, string model_name, int time_period_id);
        vector<PMIP_accessor::PMIP_seasonal_reading*> get_state_from_model_id(vector<unsigned int> &position_id, unsigned int model_id);
        vector<PMIP_accessor::PMIP_seasonal_reading*> get_state_from_position_id_name(vector<unsigned int> &position_id, string model_name, int time_period_id);
        vector<PMIP_accessor::PMIP_seasonal_reading*> get_state_from_position_statement(string position_statement, unsigned int model_name, string time_period_id);

        //vector<PMIP_accessor::Ice_reading*> get_ice_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade, int time_period_id);

        //bool does_model_exist(unsigned int model_num, string time_per);
        vector<int> get_all_model_ids_from_period(int time_per_id);
        int find_time_average(int time_period_id);
        sql::ResultSet* get_all_positions();
};

class Single_PMIP_accessor {
    public:
        virtual PMIP_accessor::PMIP_seasonal_reading* get_closest_state(loc location)=0;

        virtual vector<PMIP_accessor::PMIP_seasonal_reading*> get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade)=0;
        //virtual vector<PMIP_accessor::Ice_reading*> get_ice_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade)=0;
};

class Single_Model_PMIP_accessor : public Single_PMIP_accessor, public PMIP_accessor {
    protected:
        int model_id; 
        int time_period_id;

    public:
        Single_Model_PMIP_accessor(sql_database &d, int mod_id, int time_per_id) : PMIP_accessor(d), model_id(mod_id), time_period_id(time_per_id) {}

        PMIP_accessor::PMIP_seasonal_reading* get_state_from_position(sql::ResultSet* rs);
        vector<PMIP_accessor::PMIP_seasonal_reading*> get_state_from_positions(sql::ResultSet*);

        PMIP_accessor::PMIP_seasonal_reading* get_closest_state(loc location);
        vector<PMIP_accessor::PMIP_seasonal_reading*> get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade);
        //vector<PMIP_accessor::Ice_reading*> get_ice_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade) { return PMIP_accessor::get_ice_area(lat_max, lat_min, lon_max, lon_min, grid_degrade, time_period_id); };
};

class PMIP_ensemble : public PMIP_accessor {
    protected:
        vector<Single_Model_PMIP_accessor> members;

    public:
        PMIP_ensemble(sql_database &d) : PMIP_accessor(d) {}
};

class Average_PMIP_ensemble : public PMIP_ensemble, public Single_PMIP_accessor {
    protected:
        int time_period_id =-1;

    public:
        Average_PMIP_ensemble(sql_database &d, int time_per);

        PMIP_accessor::PMIP_seasonal_reading* get_closest_state(loc location);

        vector<PMIP_accessor::PMIP_seasonal_reading*> get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade);
        //vector<PMIP_accessor::Ice_reading*> get_ice_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade) { return PMIP_accessor::get_ice_area(lat_max, lat_min, lon_max, lon_min, grid_degrade, time_period_id); };

        void add_to_database();

    protected:
        PMIP_accessor::PMIP_seasonal_reading* average_readings(vector<PMIP_accessor::PMIP_seasonal_reading*>);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//CRU//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
class CRU_accessor {
    private:
        sql_database &db; 
        CRU_state_table state;

    public:
        class CRU_seasonal_reading : public Modern_loc_reading {
            public:
                CRU_seasonal_reading(unsigned int id, loc location) : Modern_loc_reading(id, location) {};
        };

        CRU_accessor(sql_database& d) : db(d), state(db) {}
        void load_up_tables();

        CRU_seasonal_reading* get_state_from_sql_result(sql::ResultSet* found_state);
        CRU_seasonal_reading* get_state_from_id(unsigned int id);
        CRU_seasonal_reading* get(loc location);
        CRU_seasonal_reading* get_closest(loc location, double max_grid_radius=1.0/3);
        vector<CRU_seasonal_reading*> get_closest_set(vector<loc>);

        vector<CRU_seasonal_reading*> get_area(double lat_max, double lat_min, double lon_max, double lon_min, double grid_degrade);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//OBS//////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
class OBS_accessor {
    private:
        sql_database &db; 
        OBS_obs_table obs; 

        void load_in_bart(bartlein_csv_file& bartlein_data, int time_period_id, CRU_accessor& modern_values);

    public:
        OBS_accessor(sql_database& d) : db(d), obs(d) { }//obs_result_set = NULL; }
        //~OBS_accessor() { if(obs_result_set) delete obs_result_set; }
        void load_up_tables(CRU_accessor&);

        class OBS_obs_reading : public am_sd_loc_reading {
            public:
                OBS_obs_reading(unsigned int id, loc location, bool co2) : am_sd_loc_reading(id, location, co2) {};
        }; 

        OBS_obs_reading* get_from_sql_result(sql::ResultSet* found_obs);
        OBS_obs_reading* get_obs_on_id(unsigned int id);
        OBS_obs_reading* get_closest(loc location);
        vector<OBS_obs_reading*> get_area(double lat_max, double lat_min, double lon_max, double lon_min, int time_period_id);
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Analysis/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
class ANA_accessor {
    protected:
        sql_database &db; 
        ANA_seasonal_table seasonal_state_obs; 
        ANA_experiment_table experiments;

    public:
        ANA_accessor(sql_database &d) : db(d), seasonal_state_obs(d), experiments(d) {};

        void load_up_tables();
        void add_record(int position_id, int model_id, Vec_f state, Vec_f obs, Vec_f obs_sd, Vec_f attrib);
        void add_seasonal_record(int position_id, int model_id, Vec_f state, Vec_f obs, Vec_f obs_sd, Vec_f attrib);
};

class ANA_experiment : public ANA_accessor {
    private:
        int experiment_id;

    public:
        ANA_experiment(sql_database &d, string experi_name, int time_period_id, double temp_scale, double state_grid_scale, double obs_grid_scale, int grid_degrade);

        void add_record(int position_id, Vec_f state, Vec_f obs, Vec_f obs_sd, Vec_f attrib);
        void add_seasonal_record(int position_id, Vec_f state, Vec_f obs, Vec_f obs_sd, Vec_f attrib);

        int get_experiment_id() { return experiment_id; }
};
#endif
