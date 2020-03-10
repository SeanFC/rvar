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
        //The settings to use in an experiment
        
        string experiment_name = "";        // A name to give to the experiment when saving
        double rh_scale = 1;                // The relative humidity temporal length scale. Deprecated
        double temp_scale = 1;              // The temperature temporal length scale
        double grid_scale = 0.15;           // The spatial length scale

        int time_period_id = 1;             // The time period code for the time period in the database (generally 1 is LGM and 2 is MH)
        bool use_CO2 = true;                // Correct for changes in CO2
        bool exclude_PSM_obs = false;       // Exclude observations made using model inversion
        
        string sql_username = "";           // Username for SQL database
        string sql_password = "";           // Password for SQL database
        int sql_port = 3306;                // Port of SQL database
        int grid_degrade = 1;               // Only use every grid_degrade grid cell to make the background, useful for testing
        bool tunneling = false;             // Tunnel to a remote database. Deprecated
        int connection_port = 3307;         // The port used to connect to the remote database
        bool sql_verbose = false;           // Print out SQL commands
        string tunnel_address = "";         // The address of the SQL database to tunnel to, currently deprecated

        double lat_max = 21.25;             // The maximum latitude to consider
        double lat_min = 20.8;              // The minimum latitude to consider
        double lon_max = 18.6;              // The maximum logitude to consider
        double lon_min = 18.25;             // The minimum logitude to consider

        int number_of_threads = 0;          // How many threads to use if not using open_blas_threading. Deprecated
        bool open_blas_threading = false;   // Leave the multithreading to openBLAS, this is incompatible with using our own multithreading. Deprecated

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
