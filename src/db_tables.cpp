#include "db_tables.h"

//////////////////////////////////////////////////////////////////////////////////////
// Constants for all the SQL table headers ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

const char Seasonal_state_table::seasonal_state_name_type[] = 
"MAP float, "
"temp_1 float, temp_2 float, temp_3 float, temp_4 float, temp_5 float, temp_6 float, temp_7 float, temp_8 float, temp_9 float, temp_10 float, temp_11 float, temp_12 float";
const char Seasonal_state_table::seasonal_state_name[] =
"MAP, "
"temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8, temp_9, temp_10, temp_11, temp_12";

const char Modern_state_table::modern_state_name_type[] = 
"MAP float, "
"temp_1 float, temp_2 float, temp_3 float, temp_4 float, temp_5 float, temp_6 float, temp_7 float, temp_8 float, temp_9 float, temp_10 float, temp_11 float, temp_12 float, "
"reh_1 float, reh_2 float, reh_3 float, reh_4 float, reh_5 float, reh_6 float, reh_7 float, reh_8 float, reh_9 float, reh_10 float, reh_11 float, reh_12 float, " 
"sunp_1 float, sunp_2 float, sunp_3 float, sunp_4 float, sunp_5 float, sunp_6 float, sunp_7 float, sunp_8 float, sunp_9 float, sunp_10 float, sunp_11 float, sunp_12 float";
const char Modern_state_table::modern_state_name[] =
"MAP, "
"temp_1, temp_2, temp_3, temp_4, temp_5, temp_6, temp_7, temp_8, temp_9, temp_10, temp_11, temp_12, "
"reh_1, reh_2, reh_3, reh_4, reh_5, reh_6, reh_7, reh_8, reh_9, reh_10, reh_11, reh_12, "
"sunp_1, sunp_2, sunp_3, sunp_4, sunp_5, sunp_6, sunp_7, sunp_8, sunp_9, sunp_10, sunp_11, sunp_12";

const char alpha_moisture_table::alpha_moisture_name_type[] = "alpha float, moisture_index float, precipitation float, mat float, mtco float, mtwa float, gdd5 float";
const char alpha_moisture_table::alpha_moisture_name[] = "alpha, moisture_index, precipitation, mat, mtco, mtwa, gdd5";

const char alpha_moisture_sd_table::alpha_moisture_sd_name_type[] = "alpha_sd float, moisture_index_sd float, precipitation_sd float, mat_sd float, mtco_sd float, mtwa_sd float, gdd5_sd float";
const char alpha_moisture_sd_table::alpha_moisture_sd_name[] = "alpha_sd, moisture_index_sd, precipitation_sd, mat_sd, mtco_sd, mtwa_sd, gdd5_sd";

const char alpha_moisture_attrib_table::alpha_moisture_attrib_name_type[] = "alpha_attrib float, moisture_index_attrib float, precipitation_attrib float, mat_attrib float, mtco_attrib float, mtwa_attrib float, gdd5_attrib float";
const char alpha_moisture_attrib_table::alpha_moisture_attrib_name[] = "alpha_attrib, moisture_index_attrib, precipitation_attrib, mat_attrib, mtco_attrib, mtwa_attrib, gdd5_attrib";

const char location_table::loc_name_type[] = "lat float, lon float, elev float";
const char location_table::loc_name[] = "lat, lon, elev";

const char INFO_position_table::table_name[] = "INFO_position";
const char INFO_position_table::basic_id_name[] = "id";
const char INFO_position_table::basic_id_name_type[] = "id int";

const char INFO_time_table::table_name[] = "INFO_time";
const char INFO_time_table::basic_id_name[] = "id";
const char INFO_time_table::basic_id_name_type[] = "id int";
const char INFO_time_table::time_period_name[] = "time_period";
const char INFO_time_table::time_period_name_type[] = "time_period CHAR(10)";

const char PMIP_model_table::table_name[] = "PMIP_model";
const char PMIP_model_table::basic_id_name[] = "id";
const char PMIP_model_table::basic_id_name_type[] = "id int";
const char PMIP_model_table::model_name_name[] = "model_name";
const char PMIP_model_table::model_param_name[] = "model_name, time_period_id";
const char PMIP_model_table::time_period_id_name[] = "time_period_id";
const char PMIP_model_table::model_param_name_type[] = "model_name CHAR(30), time_period_id int";

const char PMIP_state_table::table_name[] = "PMIP_state";
const char PMIP_state_table::basic_id_name[] = "id";
const char PMIP_state_table::basic_id_name_type[] = "id int";
const char PMIP_state_table::extra_id_name_list[] = "position_id, model_id";
const char PMIP_state_table::position_id_name[] = "position_id";
const char PMIP_state_table::position_id_name_type[] = "position_id int";
const char PMIP_state_table::model_id_name[] = "model_id";
const char PMIP_state_table::model_id_name_type[] = "model_id int";
const char PMIP_state_table::extra_id_name_type[] = "position_id int, model_id int";

const char PMIP_sd_table::table_name[] = "PMIP_state_sd";
const char PMIP_sd_table::basic_id_name[] = "id";
const char PMIP_sd_table::basic_id_name_type[] = "id int";
const char PMIP_sd_table::extra_id_name_list[] = "position_id, time_period_id";
const char PMIP_sd_table::position_id_name[] = "position_id";
const char PMIP_sd_table::position_id_name_type[] = "position_id int";
const char PMIP_sd_table::time_period_id_name[] = "time_period_id";
const char PMIP_sd_table::time_period_id_name_type[] = "time_period_id int";
const char PMIP_sd_table::extra_id_name_type[] = "position_id int, time_period_id int";

const char PMIP_ice_table::table_name[] = "PMIP_ice";
const char PMIP_ice_table::basic_id_name[] = "id";
const char PMIP_ice_table::basic_id_name_type[] = "id int";
const char PMIP_ice_table::time_period_id_name[] = "time_period_id";
const char PMIP_ice_table::time_period_id_name_type[] = "time_period_id int";
const char PMIP_ice_table::position_id_name[] = "position_id";
const char PMIP_ice_table::position_id_name_type[] = "position_id int";
const char PMIP_ice_table::ice_variables_name[] = "ice_percentage";
const char PMIP_ice_table::ice_variables_name_type[] = "ice_percentage int";

const char CRU_state_table::table_name[] = "CRU_state";
const char CRU_state_table::basic_id_name[] = "id";
const char CRU_state_table::basic_id_name_type[] = "id int";

const char OBS_obs_table::table_name[] = "OBS";
const char OBS_obs_table::basic_id_name[] = "id";
const char OBS_obs_table::basic_id_name_type[] = "id int";
const char OBS_obs_table::time_period_id_name[] = "time_period_id";
const char OBS_obs_table::time_period_id_name_type[] = "time_period_id int";
const char OBS_obs_table::co2_corrected_name[] = "CO2_corrected";
const char OBS_obs_table::co2_corrected_name_type[] = "CO2_corrected boolean";

const char ANA_experiment_table::table_name[] = "ANA_experiment";
const char ANA_experiment_table::basic_id_name[] = "id";
const char ANA_experiment_table::basic_id_name_type[] = "id int";
const char ANA_experiment_table::experiment_name_name[] = "experiment_name";
const char ANA_experiment_table::experiment_name_name_type[] = "experiment_name CHAR(30)";
const char ANA_experiment_table::time_period_id_name[] = "time_period_id";
const char ANA_experiment_table::time_period_id_name_type[] = "time_period_id int";
const char ANA_experiment_table::length_scales_name[] = "rh_scale, temp_scale, grid_scale";
const char ANA_experiment_table::length_scales_name_type[] = "rh_scale float, temp_scale float, grid_scale float";
const char ANA_experiment_table::extra_constants_name[] = "grid_degrade";
const char ANA_experiment_table::extra_constants_name_type[] = "grid_degrade int";

const char ANA_seasonal_table::table_name[] = "ANA_seasonal";
const char ANA_seasonal_table::basic_id_name[] = "id";
const char ANA_seasonal_table::basic_id_name_type[] = "id int";
const char ANA_seasonal_table::position_id_name[] = "position_id";
const char ANA_seasonal_table::position_id_name_type[] = "position_id int";
const char ANA_seasonal_table::experiment_id_name[] = "experiment_id";
const char ANA_seasonal_table::experiment_id_name_type[] = "experiment_id int";

/* TODO:
 * One can make better use of these commands in the sql db interface
 *  sql::Connection::isValid() checks whether the connection is alive
 *  sql::Connection::reconnect() reconnects if the connection has gone down
 *  OPT_RECONNECT option might be useful
 */

//Connect to a sql clim database
sql_database::sql_database(string username, string password, int port, string tunnel_address, bool tun, int con_port, bool sql_vb) : tunnel(tun), connection_port(con_port), sql_verbose(sql_vb) {
    driver = get_driver_instance();
    
    open_connection(username, password, port, tunnel_address);
}

sql_database::~sql_database() {
    close_connection();
}

void sql_database::open_connection(string username, string password, int port, string tunnel_address) {
    //If the server is remote and the sql interface isn't open we set up an ssh tunnel to the machine at the connection port and interface that way
    if(tunnel) {
        //TODO: Currently unsupported, create a tunnel using bac_rmin instead
    }

    con = shared_ptr<sql::Connection>(driver->connect(("tcp://127.0.0.1:"+to_string(connection_port)).c_str(), username.c_str(), password.c_str())); 
    //con = driver->connect("tcp://localhost", username.c_str(), password.c_str()); 

    //Incase the database doesn't exist we make it and start to use it
    sql::Statement *stmt = con->createStatement();
    if(sql_verbose) printf("Create database: %s\n", "CREATE DATABASE IF NOT EXISTS clim");
    stmt->execute("CREATE DATABASE IF NOT EXISTS clim");
    con->setSchema("clim");
    delete stmt;
}

void sql_database::close_connection() {
    if(tunnel) {
        //TODO: This is currently unsupported
        //if(system(("ssh -S ./sql-socket_"+to_string(connection_port)+" -O exit INSERT HOSTNAME STRING").c_str()))  {
        //    printf("Error: System call to close ssh socket had a problem\n");
        //}
    }
}

//Delete and then initalise the() sql table
void sql_table::reset() {
    delete_table();
    initialise();
}

//Delete the table
void sql_table::delete_table() {
    sql::Statement *stmt = base.get_connection()->createStatement();

    if(base.verbose()) printf("Drop table: %s\n", ("DROP TABLE IF EXISTS "+get_table_name()).c_str());
    stmt->execute("DROP TABLE IF EXISTS "+get_table_name());

    delete stmt;
}


//Initialise the sql table
void sql_table::initialise() {
    sql::Statement *stmt = base.get_connection()->createStatement();
    
    if(base.verbose()) printf("Create table: %s\n", ("CREATE TABLE IF NOT EXISTS "+get_table_name()+"("+get_creation_string()+")").c_str());
    stmt->execute("CREATE TABLE IF NOT EXISTS "+get_table_name()+"("+get_creation_string()+")");

    delete stmt;
}

sql::PreparedStatement* sql_table::generate_insert_prepared_statement(string variable_list, int amount_of_variables) {
    string template_string = "INSERT INTO "+get_table_name()+"(" + variable_list +") VALUES (?";

    //This -1 is from the ? above
    for(int i=0; i<amount_of_variables-1; i++) { 
        template_string += ", ?";
    }
    template_string += ")";

    if(base.verbose()) printf("Preparing: %s\n", template_string.c_str());

    return base.get_connection()->prepareStatement(template_string.c_str());
}

sql::ResultSet* sql_table::run_command(string command) {
    if(base.verbose()) printf("Command: %s\n", (command).c_str());

    sql::Statement *s = base.get_connection()->createStatement();
    sql::ResultSet* rs = s->executeQuery(command);

    delete s;
    return rs;
}

sql::ResultSet* sql_table::get_on_condition(string condition) {
    string sql_ask_string = get_on_condition_string(condition);

    if(base.verbose()) printf("Condition: %s\n", (sql_ask_string).c_str());

    sql::Statement *s = base.get_connection()->createStatement();
    sql::ResultSet* rs = s->executeQuery(sql_ask_string);

    delete s;
    return rs;
}

string sql_table::get_on_condition_string(string condition) {
    return "SELECT * FROM " + get_table_name() + " WHERE "+ condition;
}


sql::ResultSet* sql_table::get_matching(string to_match, string match_to) {
    sql::Statement *s = base.get_connection()->createStatement();
    if(base.verbose()) printf("Matching: %s\n", ("SELECT * FROM " + get_table_name() + " WHERE "+ to_match + "=" + match_to).c_str());
    sql::ResultSet* rs = s->executeQuery("SELECT * FROM " + get_table_name() + " WHERE "+ to_match + "=" + match_to);
    delete s;
    return rs;
}

sql::ResultSet* sql_table::get_all() {
    sql::Statement *s = base.get_connection()->createStatement();
    if(base.verbose()) printf("Get all: %s\n", ("SELECT * FROM " + get_table_name()).c_str());
    sql::ResultSet* rs = s->executeQuery("SELECT * FROM " + get_table_name());
    delete s;
    return rs;
}

sql::ResultSet* sql_table::get_last(string variable_name) {
    sql::Statement *s = base.get_connection()->createStatement();
    if(base.verbose()) printf("Get last: %s\n", ("SELECT " + variable_name + " FROM " + get_table_name() + " ORDER BY " + variable_name + " DESC LIMIT 1").c_str());
    sql::ResultSet* rs = s->executeQuery("SELECT " + variable_name + " FROM " + get_table_name() + " ORDER BY " + variable_name + " DESC LIMIT 1");
    delete s;
    return rs;
}

string location_table::get_box_condition(double low_lat, double high_lat, double low_lon, double high_lon) {
    ostringstream strs; 
    strs << "( lat BETWEEN " << low_lat << " AND " << high_lat << ")";
    strs << " AND (lon BETWEEN " << low_lon << " AND " << high_lon << ")";
    return strs.str();
}

//////////////////////////////////////////////////////////////////////////////////////////
//INFO////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
bool INFO_position_table::add_record(loc location) {
    int variable_amount = loc_variable_amount;

    sql::PreparedStatement *ps = generate_insert_prepared_statement(loc_name, variable_amount);

    ps->setDouble(1, location.lat);
    ps->setDouble(2, location.lon);
    ps->setDouble(3, location.elev);

    ps->execute();

    delete ps;

    return true;
}


int INFO_position_table::get_id_from_lat_lon(double lat, double lon) {
    sql::ResultSet* rs = run_command("SELECT id FROM " + get_table_name() + " WHERE lat = " +to_string(lat)+ " AND lon = " + to_string(lon));
    int id = -1;

    if(rs->next())
        id = rs->getInt(1);

    delete rs;
    return id;
}

bool INFO_time_table::add_record(string time_period) {

    sql::PreparedStatement *ps = generate_insert_prepared_statement(time_period_name, 1);

    ps->setString(1, time_period);
    ps->execute();

    delete ps;
    return true;
}

//TODO:Check for existance of solution
int INFO_time_table::get_id_from_time_period(string time_period) { 
    sql::ResultSet* rs = get_matching(time_period_name, "\"" + time_period + "\""); //Since time period is a string we need to escape it
    int output = -1;

    if(rs->next())
        output = rs->getInt(1);

    delete rs;
    return output;
}

//////////////////////////////////////////////////////////////////////////////////////////
//PMIP////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
bool PMIP_model_table::add_record(string given_model_name, int time_period_id) {
    sql::PreparedStatement *ps = generate_insert_prepared_statement(model_param_name, params_length);

    ps->setString(1, given_model_name);
    ps->setInt(2, time_period_id);

    ps->execute();

    delete ps;
    return true;
}

int PMIP_model_table::get_last_id() {
    sql::ResultSet *rs = get_last(basic_id_name);
    int id=-1;
    if(rs->next()) 
        id = rs->getInt(1);

    delete rs;
    return id;
}

bool PMIP_state_table::add_record(int ensamble_id, int model_id, Vec_f state) {
    int variable_amount = no_of_ids - 1 + modern_state_variable_amount;

    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(extra_id_name_list) + ", " + modern_state_name, variable_amount);

    ps->setInt(1, ensamble_id);
    ps->setInt(2, model_id);

    for(int j=0; j<state.get_size(); j++) {
        if(!isnan(state.get(j)))
            ps->setDouble(j+no_of_ids, state.get(j)); //We take away 1 since we don't include the main id and we start at 1 since the sql connector starts itteratoring at 1 hence we end up not doing anything
        else
            ps->setNull(j+no_of_ids, sql::DataType::DOUBLE);
    }

    ps->execute();

    delete ps;
    return true;
}

bool PMIP_sd_table::add_record(int position_id, int time_period_id, Vec_f sds) {
    int variable_amount = no_of_ids - 1 + modern_state_variable_amount;

    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(extra_id_name_list) + ", " + modern_state_name, variable_amount);

    ps->setInt(1, position_id);
    ps->setInt(2, time_period_id);

    for(int j=0; j<sds.get_size(); j++) {
        ps->setDouble(j+3, sds.get(j)); //2 are from the ints above, the extra is since the sql connector counts from 1
    }

    ps->execute();

    delete ps;
    return true;
}

bool PMIP_ice_table::add_record(int time_period_id, int position_id, double ice_percentage) {
    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(time_period_id_name) + ", " + string(position_id_name) + ", " + string(ice_variables_name), 3);

    ps->setInt(1, time_period_id);
    ps->setInt(2, position_id);
    ps->setDouble(3, ice_percentage); 

    ps->execute();

    delete ps;
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////
//CRU///////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
bool CRU_state_table::add_record(loc location, Vec_f sta) {
    int variable_amount = no_of_ids - 1 + loc_variable_amount + modern_state_variable_amount;
    
    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(loc_name) + ", " + modern_state_name, variable_amount);

    ps->setDouble(1, location.lat);
    ps->setDouble(2, location.lon);
    ps->setDouble(3, location.elev);
    
    int starting_location = no_of_ids - 1 + loc_variable_amount + 1; //Since we miss the normal id, sql starts counting from 1
    for(int i=0; i<sta.get_size(); i++) {
        ps->setDouble(starting_location +i, sta.get(i));
    }

    ps->execute();

    delete ps;

    return true;
}

///////////////////////////////////////////////////////////////////
//Observations/////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
bool OBS_obs_table::add_record(int time_per_id, loc location, Vec_f obs, Vec_f sds, bool CO2_corrected) {
    int variable_amount = no_of_ids -1 + loc_variable_amount + no_of_co2_vars + alpha_moisture_variable_amount + alpha_moisture_sd_variable_amount;

    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(time_period_id_name) + ", " + loc_name + ", " + co2_corrected_name + ", " + alpha_moisture_name + ", " + alpha_moisture_sd_name , variable_amount);

    ps->setInt(1, time_per_id);
    ps->setDouble(2, location.lat);
    ps->setDouble(3, location.lon);
    ps->setDouble(4, location.elev);
    ps->setBoolean(5, CO2_corrected);

    //Go through the record's values and set all the variables in the right places of the statement
    int starting_location = no_of_ids - 1 + loc_variable_amount + no_of_co2_vars + 1; //Since we miss the normal id, sql starts counting from 1
    for(int i=0; i<obs.get_size(); i++) {
        
        //If the value of the obs or the sds isn't NAN then put it in the statement, else set that part of the statement to NULL
        if(not isnan(obs.get(i))) 
            ps->setDouble(starting_location +i, obs.get(i)); 
        else 
            ps->setNull(starting_location +i, sql::DataType::DOUBLE);

        if(not isnan(sds.get(i))) 
            ps->setDouble(starting_location + alpha_moisture_variable_amount + i, sds.get(i)); 
        else 
            ps->setNull(starting_location + alpha_moisture_variable_amount + i, sql::DataType::DOUBLE);
    }

    //Launch the statement so the record is added to the table
    ps->execute();

    //Clean up and finish
    delete ps;
    return true;
}

///////////////////////////////////////////////////////////////////
//Analysis/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
bool ANA_experiment_table::add_record(string given_experiment_name, int time_period_id, double rh_scale, double temp_scale, double grid_scale, int grid_degrade) {
    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(experiment_name_name) + ", " + time_period_id_name + ", " + length_scales_name + ", " + extra_constants_name, 6);
    
    ps->setString(1, given_experiment_name);
    ps->setInt(2, time_period_id);
    ps->setDouble(3, rh_scale);
    ps->setDouble(4, temp_scale);
    ps->setDouble(5, grid_scale);
    ps->setInt(6, grid_degrade);

    ps->execute();

    delete ps;

    return true;
}

int ANA_experiment_table::get_last_id() {
    sql::ResultSet *rs = get_last(basic_id_name);
    int id=-1;
    if(rs->next()) 
        id = rs->getInt(1);

    delete rs;
    return id;
}

bool ANA_seasonal_table::add_record(int position_id, int experiment_id, Vec_f seasonal_state, Vec_f obs, Vec_f obs_sds, Vec_f attrib) {
    int variable_amount = no_of_ids - 1 + seasonal_state_variable_amount + alpha_moisture_variable_amount + alpha_moisture_sd_variable_amount + alpha_moisture_attrib_variable_amount;

    sql::PreparedStatement *ps = generate_insert_prepared_statement(string(position_id_name) + ", " + experiment_id_name + ", " + seasonal_state_name + ", " + alpha_moisture_name + ", " + alpha_moisture_sd_name + ", " + alpha_moisture_attrib_name, variable_amount);

    ps->setInt(1, position_id);
    ps->setInt(2, experiment_id);

    int starting_location = no_of_ids;// - 1 + 1; //Since we miss the normal id, sql starts counting from 1
    int seasonal_state_size = seasonal_state.get_size();

    for(int i=0; i<seasonal_state_size; i++) {
        ps->setDouble(starting_location + i, seasonal_state.get(i));
    }
    
    for(int i=0; i<obs.get_size(); i++) {
        sql_set_double_nan_guard(ps, starting_location + seasonal_state_size + i, obs.get(i));
        sql_set_double_nan_guard(ps, starting_location + seasonal_state_size + alpha_moisture_variable_amount + i, obs_sds.get(i));
        sql_set_double_nan_guard(ps, starting_location + seasonal_state_size + alpha_moisture_variable_amount + alpha_moisture_attrib_variable_amount + i, attrib.get(i));
    }

    ps->execute();

    delete ps;

    return true;
}
