//All C++ table interfaces in the database, one for each table
#ifndef clim_db_interface_h
#define clim_db_interface_h

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>

#include <cppconn/driver.h>
#include <cppconn/exception.h>
#include <cppconn/resultset.h>
#include <cppconn/statement.h>
#include <cppconn/prepared_statement.h>
#include "mysql_connection.h"

//#include <libssh/libssh.h>

#include "lin_maths.h"
#include "location.h"

using namespace std;

class sql_database {
    private:
        shared_ptr<sql::Connection> con; 
        sql::Driver* driver = NULL;
        bool tunnel; //Are we connecting to a remote server (and so we must set up an ssh tunnel first) 
        int connection_port;
        bool sql_verbose;
    
    public:
        sql_database(string username, string password, int port, string tunnel_address = "", bool tun = false, int con_port = 3307, bool sql_verbose=false);

        virtual ~sql_database();

        void open_connection(string username, string password, int port, string tunnel_address);
        void close_connection();

        shared_ptr<sql::Connection> get_connection()  { return con; }
        bool verbose() { return sql_verbose; }
};

class sql_table {
        sql_database &base;

    public:
        sql_table(sql_database& b) : base(b) {};
        void reset();
        void delete_table();
        void initialise();

        virtual string get_table_name()=0;
        virtual string get_creation_string()=0;

        sql::PreparedStatement* generate_insert_prepared_statement(string variable_list, int amount_of_variables);
        sql::ResultSet* run_command(string command);
        sql::ResultSet* get_on_condition(string condition);
        string get_on_condition_string(string condition);
        sql::ResultSet* get_matching(string to_match, string match_to);
        sql::ResultSet* get_all();
        sql::ResultSet* get_last(string variable_name);

        static inline void sql_set_double_nan_guard(sql::PreparedStatement *ps, int position, double value) {
            if(not isnan(value)) ps->setDouble(position, value); else ps->setNull(position, sql::DataType::DOUBLE);
        }


};

/*Database structure
* ANA_experiment, ANA_seasonal, CRU_state, INFO_position, INFO_time, OBS, PMIP_model, PMIP_state, PMIP_state_sd
 */

class Seasonal_state_table {
    public:
        static const char seasonal_state_name_type[];
        static const char seasonal_state_name[];
        static const int seasonal_state_variable_amount = 13;
};

class Modern_state_table {
    public:
        static const char modern_state_name_type[];
        static const char modern_state_name[];
        static const int modern_state_variable_amount = 37;
};

class alpha_moisture_table {
    public:
        static const char alpha_moisture_name_type[];
        static const char alpha_moisture_name[];
        static const int alpha_moisture_variable_amount = 7;
};

class alpha_moisture_sd_table {
    public:
        static const char alpha_moisture_sd_name_type[];
        static const char alpha_moisture_sd_name[];
        static const int alpha_moisture_sd_variable_amount = 7;
};

class alpha_moisture_attrib_table {
    public:
        static const char alpha_moisture_attrib_name_type[];
        static const char alpha_moisture_attrib_name[];
        static const int alpha_moisture_attrib_variable_amount = 7;
};


class location_table {
    public:
        static const char loc_name_type[];
        static const char loc_name[];
        static const int loc_variable_amount = 3;

        string get_box_condition(double low_lat, double high_lat, double low_lon, double high_lon);
};

//////////////////////////////////////////////////////////////////////////////////////
//INFO tables/////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//INFO_position : id, lat, lon, elev
class INFO_position_table : public sql_table, public location_table {
    public:
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];

        INFO_position_table(sql_database &b) : sql_table(b) {}

        bool add_record(loc location);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                loc_name_type + ", "+ 
                "PRIMARY KEY(" + basic_id_name + ")"; }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
        int get_id_from_lat_lon(double lat, double lon);
};

//INFO_time : id, time_period
class INFO_time_table : public sql_table, public location_table {
    public:
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char time_period_name[];
        static const char time_period_name_type[];

        INFO_time_table(sql_database &b) : sql_table(b) {}

        bool add_record(string time_period);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                time_period_name_type + ", "+ 
                "PRIMARY KEY (" + basic_id_name + ")"; }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
        int get_id_from_time_period(string time_period);
};

///////////////////////////////////////////////////////////////////////////////////
//PMIP/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//PMIP_model : id, model_name
class PMIP_model_table : public sql_table {
    public:
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char model_name_name[];
        static const char model_param_name[];
        static const char time_period_id_name[];
        static const char model_param_name_type[];
        static const int params_length = 2;

        PMIP_model_table(sql_database &b) : sql_table(b) {}

        bool add_record(string given_model_name, int time_period_id);
        int get_last_id();

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                model_param_name_type + ", " +
                "PRIMARY KEY (" + basic_id_name + "), " +
                "FOREIGN KEY (" + time_period_id_name + ") REFERENCES " + INFO_time_table::table_name + "(" + INFO_time_table::basic_id_name + ") ON DELETE CASCADE";
                 }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

//PMIP_state : id, ensamble_id, model_id, state
class PMIP_state_table : public sql_table, public Modern_state_table {
    public:
        static const int no_of_ids = 3;
        static const char basic_id_name_type[];
        static const char extra_id_name_list[];
        static const char position_id_name[];
        static const char position_id_name_type[];
        static const char model_id_name[];
        static const char model_id_name_type[];
        static const char extra_id_name_type[];

        PMIP_state_table(sql_database &b) : sql_table(b) {}

        static const char table_name[];
        static const char basic_id_name[];

        bool add_record(int position_id, int model_id, Vec_f state);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                position_id_name_type + " NOT NULL, " +
                model_id_name_type + " NOT NULL, " +
                modern_state_name_type + ", " + 
                "PRIMARY KEY (" + basic_id_name + "), " +
                "FOREIGN KEY (" + position_id_name + ") REFERENCES " + INFO_position_table::table_name + "(" + INFO_position_table::basic_id_name + ") ON DELETE CASCADE, " +
                "FOREIGN KEY (" + model_id_name + ") REFERENCES " + PMIP_model_table::table_name + "(" + PMIP_model_table::basic_id_name + ") ON DELETE CASCADE"; }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

//PMIP_sd : id, PMIP_position.id, PMIP_time.id, sds
class PMIP_sd_table : public sql_table, public Modern_state_table {
    private:
        static const int no_of_ids = 3;

    public:
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char extra_id_name_list[];
        static const char extra_id_name_type[];
        static const char position_id_name[];
        static const char position_id_name_type[];
        static const char time_period_id_name[];
        static const char time_period_id_name_type[];

        PMIP_sd_table(sql_database &d) : sql_table(d) {}
        bool add_record(int position_id, int time_period_id, Vec_f sds);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                position_id_name_type + " NOT NULL, " +
                time_period_id_name_type + " NOT NULL, " +
                modern_state_name_type + ", " + 
                "PRIMARY KEY(" + basic_id_name + "), " +
                "FOREIGN KEY(" + position_id_name + ") REFERENCES " + INFO_position_table::table_name + "(" + INFO_position_table::basic_id_name + ") ON DELETE CASCADE, " +
                "FOREIGN KEY(" + time_period_id_name + ") REFERENCES " + INFO_time_table::table_name + "(" + INFO_time_table::basic_id_name + ") ON DELETE CASCADE"; }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

//PMIP_ice : id, time_period_id, position_id, ice_percentage 
class PMIP_ice_table : public sql_table {
    public:
        //static const int no_of_ids = 1;

        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char time_period_id_name[];
        static const char time_period_id_name_type[];
        static const char position_id_name[];
        static const char position_id_name_type[];
        static const char ice_variables_name[];
        static const char ice_variables_name_type[];
        
        PMIP_ice_table(sql_database &b) : sql_table(b) {}

        static const char table_name[];

        bool add_record(int time_period_id, int position_id, double ice_percentage);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                position_id_name_type + " NOT NULL, " +
                time_period_id_name_type + " NOT NULL, " + 
                ice_variables_name_type + ", " + 
                "PRIMARY KEY (" + basic_id_name + "), " +
                "FOREIGN KEY(" + time_period_id_name + ") REFERENCES " + INFO_time_table::table_name + "(" + INFO_time_table::basic_id_name + ") ON DELETE CASCADE, " + 
                "FOREIGN KEY (" + position_id_name + ") REFERENCES " + INFO_position_table::table_name + "(" + INFO_position_table::basic_id_name + ") ON DELETE CASCADE";
        }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

////////////////////////////////////////////////////////////////////////////////////
//CRU///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
//CRU_state : id, lat, lon, elev, state
class CRU_state_table : public sql_table, public location_table, public Modern_state_table {
    public:
        static const int no_of_ids = 1;
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];

        CRU_state_table(sql_database &b) : sql_table(b) {}
        bool add_record(loc location, Vec_f sta);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                loc_name_type + ", "+ 
                modern_state_name_type + ", " + 
                "PRIMARY KEY(" + basic_id_name + ")";
                }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

///////////////////////////////////////////////////////////////////////////////////
//Observations/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//OBS: id, PMIP_time(id), loc, obs, sd
class OBS_obs_table : public sql_table, public location_table, public alpha_moisture_table, public alpha_moisture_sd_table {
    private:
        static const int no_of_ids = 2;
        static const int no_of_co2_vars = 1;
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char time_period_id_name[];
        static const char time_period_id_name_type[];
        static const char co2_corrected_name[];
        static const char co2_corrected_name_type[];

    public:
        OBS_obs_table(sql_database &b) : sql_table(b) {}
        bool add_record(int time_per_id, loc location, Vec_f obs, Vec_f sds, bool CO2_corrected=false);  //Note we expect a 7 length obs here

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                time_period_id_name_type + ", " +
                loc_name_type + ", "+ 
                co2_corrected_name_type + ", "+ 
                alpha_moisture_name_type + ", " + 
                alpha_moisture_sd_name_type + ", " +
                "PRIMARY KEY(" + basic_id_name + "), " + 
                "FOREIGN KEY (" + time_period_id_name + ") REFERENCES " + INFO_time_table::table_name + "(" + INFO_time_table::basic_id_name + ") ON DELETE CASCADE";
                }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Analysis/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
//ANA_experiment : id, name, PMIP_time.id
class ANA_experiment_table : public sql_table {
    public:
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char experiment_name_name[];
        static const char experiment_name_name_type[];
        static const char time_period_id_name[];
        static const char time_period_id_name_type[];
        static const char length_scales_name[];
        static const char length_scales_name_type[];
        static const char extra_constants_name[];
        static const char extra_constants_name_type[];

        static const int params_length = 3;

        ANA_experiment_table(sql_database &b) : sql_table(b) {}

        bool add_record(string given_model_name, int time_period_id, double temp_scale, double state_scale, double obs_scale, int grid_degrade);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                experiment_name_name_type + ", " +
                time_period_id_name_type + ", " +
                length_scales_name_type + ", " + 
                extra_constants_name_type + ", " +
                "PRIMARY KEY (" + basic_id_name + "), " +
                "FOREIGN KEY (" + time_period_id_name + ") REFERENCES " + INFO_time_table::table_name + "(" + INFO_time_table::basic_id_name + ") ON DELETE CASCADE";
                 }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
        int get_last_id();
};

//ANA_seasonal : id, experiment_id, PMIP_position.id, Seasonal_state, obs, obs sd, obs attrib
class ANA_seasonal_table : public sql_table, public Seasonal_state_table, public alpha_moisture_table, public alpha_moisture_sd_table, public alpha_moisture_attrib_table {
    public:
        static const int no_of_ids = 3;
        static const char table_name[];
        static const char basic_id_name[];
        static const char basic_id_name_type[];
        static const char position_id_name[];
        static const char position_id_name_type[];
        static const char experiment_id_name[];
        static const char experiment_id_name_type[];

        ANA_seasonal_table(sql_database &b) : sql_table(b) {}

        bool add_record(int position_id, int experiment_id, Vec_f seasonal_state, Vec_f obs, Vec_f obs_sd, Vec_f attrib);

        string get_table_name() { return table_name; }
        string get_creation_string() { 
            return string(basic_id_name_type) + " NOT NULL AUTO_INCREMENT, " +
                position_id_name_type + " NOT NULL, " +
                experiment_id_name_type + " NOT NULL, " +
                seasonal_state_name_type + ", " + 
                alpha_moisture_name_type + ", " + 
                alpha_moisture_sd_name_type + ", " + 
                alpha_moisture_attrib_name_type + ", " + 
                "PRIMARY KEY (" + basic_id_name + "), " +
                "FOREIGN KEY (" + position_id_name + ") REFERENCES " + INFO_position_table::table_name + "(" + INFO_position_table::basic_id_name + ") ON DELETE CASCADE, " +
                "FOREIGN KEY (" + experiment_id_name + ") REFERENCES " + ANA_experiment_table::table_name + "(" + ANA_experiment_table::basic_id_name + ") ON DELETE CASCADE"; }

        sql::ResultSet* get_from_id(int id) { return get_matching(basic_id_name, to_string(id)); } 
};

#endif
