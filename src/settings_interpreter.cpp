#include "settings_interpreter.h"

//Use the default file name of config/"hostname of computer".cfg and set all the settings
Assim_input_parser::Assim_input_parser() {
    char hostname[1024];
    gethostname(hostname, 1024);
    
    parse_config_file(string("config/") + hostname+ ".cfg");
}

//Process each line and remove anything in the line in a given file after the comment character ('//')
void Assim_input_parser::parse_config_file(string file_path) {
    fstream config_file(file_path);

    if(config_file.is_open()) {
        string line;
        while( getline(config_file, line) )
            parse_line(line.substr(0, line.find("//")));
        
        config_file.close();
    }
}

//Split a line up about the end line character (';')
void Assim_input_parser::parse_line(string line) {
    istringstream is_line(line);

    string segment;
    while( getline(is_line, segment, ';') )
        parse_segment(segment);
}

//Process a segemnt, it should be key=value
void Assim_input_parser::parse_segment(string segment) {
    istringstream is_line(segment);
    string key;

    if( getline(is_line, key, '=') ) {
        string value;
        if( getline(is_line, value) ) 
            set_value(key, value);
    }
}

//Set a particular value, making sure to use the correct type
//TODO:This could be better implemented with a functional and some arrays
void Assim_input_parser::set_value(string key, string value) {
    //Get rid of white space at the beginning and end of the values
    rtrim(ltrim(key));
    rtrim(ltrim(value));

    //TODO:No error checking for sql injection and type
    if(key == "rh_scale") 
        rh_scale = stod(value);
    else if(key == "temp_scale") 
        temp_scale = stod(value);
    else if (key == "grid_scale")
        grid_scale = stod(value);
    else if(key == "sql_username")
        sql_username = value;
    else if(key == "sql_password")
        sql_password = value;
    else if(key == "sql_port")
        sql_port = stoi(value); 
    else if(key == "time_period_id")
        time_period_id = stoi(value);
    else if(key == "experiment_name")
        experiment_name = value;
    else if(key == "grid_degrade")
        grid_degrade = stoi(value);
    else if(key == "lat_min")
        lat_min = stod(value);
    else if(key == "lat_max")
        lat_max = stod(value);
    else if(key == "lon_min")
        lon_min = stod(value);
    else if(key == "lon_max")
        lon_max = stod(value);
    else if(key == "tunneling")
        tunneling = stoi(value);
    else if(key == "connection_port")
        connection_port = stoi(value);
    else if(key == "sql_verbose")
        sql_verbose = stoi(value);
    else if(key == "number_of_threads")
        number_of_threads = stoi(value);
    else if(key == "open_blas_threading")
        open_blas_threading = stoi(value);
    else if(key == "use_CO2")
        use_CO2 = stoi(value);
    else if(key == "exclude_PSM_obs")
        exclude_PSM_obs = stoi(value);
    else if(key == "tunnel_address")
        tunnel_address = value;
    else
        printf("Error: key %s not found for value %s\n", key.c_str(), value.c_str());
}
