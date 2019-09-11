#ifndef settings_interpreter_h 
#define settings_interpreter_h

#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <unistd.h>

using namespace std;

//Hold all the inputs to a particular experiment 
class Assim_input_parser {
    public:
        //The settings determine an experiment
        double rh_scale = 1;
        double temp_scale = 1;
        double grid_scale = 0.15;
        string sql_username = "";
        string sql_password = "";
        int sql_port = 3306;
        int time_period_id = 1;
        string experiment_name = "";
        int grid_degrade = 1;
        bool tunneling = false;
        int connection_port = 3307;
        bool sql_verbose = false;
        bool use_CO2 = true;
        bool exclude_PSM_obs = false;
        string tunnel_address = "";

        double lat_max = 21.25; 
        double lat_min = 20.8; 
        double lon_max = 18.6; 
        double lon_min = 18.25; 

        int number_of_threads = 0;
        bool open_blas_threading = false;

        //Constructor
        Assim_input_parser();
        Assim_input_parser(string file_path) { parse_config_file(file_path); };

        //Parse a specific line of input
        void parse_line(string line);

    private:
        //Grab all the settings from a particular file
        void parse_config_file(string file_path);

        //Parse specific parts of input
        void parse_segment(string segement);
    
        //Set a particular value
        void set_value(string key, string value);
};


// Trim all the spaces at the start of a string (in place)
static inline std::string& ltrim(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
        return !std::isspace(ch);
    }));
    return s;
}

// Trim all the spaces at the end of a string (in place)
static inline std::string& rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
        return !std::isspace(ch);
    }).base(), s.end());
    return s;
}

#endif
