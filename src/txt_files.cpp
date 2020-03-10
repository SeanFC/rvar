#include "txt_files.h"

netcdf_PMIP::netcdf_PMIP(vector<PMIP_nc_file_meta> file_metas) {
    //TODO:The path location + file name uses unix '/' so wont play well with windows etc
    string ordered_metas[4] = {"", "", "", ""}; //Not a very elegant solution, coudl just swap file_metas around

    for(auto it = file_metas.begin(); it != file_metas.end(); it++) {
        if(it->var_name == "tas")
            ordered_metas[0] = it->path_location + "/" + it->file_name;
        else if(it->var_name == "clt")
            ordered_metas[1] = it->path_location + "/" + it->file_name;
        else if(it->var_name == "pr")
            ordered_metas[2] = it->path_location + "/" + it->file_name;
        else if(it->var_name == "hurs")
            ordered_metas[3] = it->path_location + "/" + it->file_name;
    }

    experiment_name = file_metas[0].experiment_name;
    grid_resolution = file_metas[0].grid_resolution;
    project_name = file_metas[0].project_name;
    model_name = file_metas[0].model_name;

    open_up_files(ordered_metas[0], ordered_metas[1], ordered_metas[2], ordered_metas[3]);
}

void netcdf_PMIP::open_up_files(string tas_fp, string cl_fp, string pr_fp, string hur_fp) {
    if(does_file_exist(tas_fp))
        tas_file = new NcFile(tas_fp.c_str(), NcFile::read);

    if(does_file_exist(cl_fp))
        cl_file = new NcFile(cl_fp.c_str(), NcFile::read); 

    if(does_file_exist(pr_fp))
        pr_file = new NcFile(pr_fp.c_str(), NcFile::read);

    if(does_file_exist(hur_fp))
        hur_file = new NcFile(hur_fp.c_str(), NcFile::read);
}

void netcdf_PMIP::load_all_readings() {
    //TODO:This just assumes that the tas file exists, we can use any file we have here
    std::multimap<std::string,NcVar> vars = tas_file->getVars();

    //for(auto it=vars.begin(); it != vars.end(); it++)
    //    printf("Var name %s\n", it->first.c_str());

    //TODO:Check everything here such as the files have the correct dimensions and stuff 
   
    //Get ncvar handles for all the main variables
    NcVar tas_var = tas_file->getVar("tas_ltm");

    //Find the sizes of all the dimensions 
    time_size = tas_var.getDim(0).getSize();
    lat_size = tas_var.getDim(1).getSize();
    lon_size = tas_var.getDim(2).getSize();

    loaded_time = new float[time_size];
    loaded_lat = new float[lat_size];
    loaded_lon = new float[lon_size];

    tas_file->getVar("time").getVar(loaded_time);
    tas_file->getVar("lat").getVar(loaded_lat);
    tas_file->getVar("lon").getVar(loaded_lon);

    vector<NcFile*> files = {tas_file, cl_file, pr_file, hur_file};
    vector<string> var_nc_names= {"tas_ltm", "clt_ltm", "pr_ltm", "hurs_ltm" };
    vector<float*> value_arrays = vector<float*>(); 

    //If we're able to load a file here then load its variables into the given array
    for(unsigned int i=0; i<var_nc_names.size(); i++) {
        if(files[i]) {
            NcVar var = files[i]->getVar(var_nc_names[i].c_str());
            //Runtime sized multi dimensional arrays in C++ and C are absolutely ridiculous, i'm just going to do it this way
            float *hold_array = new float[time_size*lat_size*lon_size];
            var.getVar(hold_array);
            value_arrays.push_back(hold_array);
        } else {
            value_arrays.push_back(NULL);
        }
    }

    loaded_tas = value_arrays[0];
    loaded_cl  = value_arrays[1];
    loaded_pr  = value_arrays[2];
    loaded_hur = value_arrays[3];

    //Figure out how many days are in each month so that we can caluclate the yearly precip correctly
    if(loaded_time[0] != 1) {
        for(int i=0; i<time_size; i++){
            if(i!=time_size-1)
                days_in_given_months[i] = loaded_time[i+1] - loaded_time[i];
            else
                days_in_given_months[i] = 365.25 - (loaded_time[time_size-1] - loaded_time[0]);
        }
    }

}

Netcdf_PMIP_Info* netcdf_PMIP::get_lat_lon_record(int lat_index, int lon_index) {
    if(lat_index-1 > lat_size || lat_index < 0 || lon_index-1> lon_size || lon_index<0) return NULL;
    Netcdf_PMIP_Info* output_info = new Netcdf_PMIP_Info();

    //Make the climate we want to add and put it into the output
    if(loaded_lon[lon_index] < 180)
        output_info->location = loc(loaded_lat[lat_index], loaded_lon[lon_index], 0); //TODO: Not supplied an elevation here
    else 
        output_info->location = loc(loaded_lat[lat_index], loaded_lon[lon_index]-360, 0); //TODO: Not supplied an elevation here

    //output_info->values = Seasonal_Modern_reading();
    output_info->set(NAN);
    double precip_sum = 0;
    for(int time_index=0; time_index<time_size; time_index++) {
        //With this netcdf handler the last variable varies the fastest, the last variable here is longitude and hence we vary if the fastest here, see /usr/include/ncVar.h:437 for how netcdf files are read in all at once
        int data_index = time_index*(lat_size*lon_size) + lon_size*lat_index + lon_index;

        //if(loaded_pr) precip_sum += days_in_given_months[time_index]*loaded_pr[data_index];
        //printf("%f\n", loaded_pr[data_index]);
        if(loaded_pr) precip_sum += loaded_pr[data_index]; 
        if(loaded_tas) {
            //If we're using Kelvin then chagne to celcius
            if(loaded_tas[data_index] > 100)
                output_info->set(time_index +  1, loaded_tas[data_index] - 273.15);  
            else
                output_info->set(time_index +  1, loaded_tas[data_index]); 
        }
        if(loaded_hur) output_info->set(time_index + 13, loaded_hur[data_index]); 
        if(loaded_cl) output_info->set(time_index + 25, 100 - loaded_cl[data_index]); 
    }

    if(loaded_pr) output_info->set_MAP(precip_sum); 

    if(loaded_pr)  output_info->values_held.set(0, 1);
    if(loaded_tas) output_info->values_held.set(1, 1);
    if(loaded_hur) output_info->values_held.set(2, 1);
    if(loaded_cl)  output_info->values_held.set(3, 1);

    return output_info;
}
//////////////////////////////////////////////////////////////////////////////////
//netcdf ice file/////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void Netcdf_ice::load_in_info() {
    //TODO:Load up file into array
    NcVar ice_var = ice_file->getVar("sftgif");

    lat_size = ice_var.getDim(0).getSize();
    lon_size = ice_var.getDim(1).getSize();
    std::multimap<std::string,NcVar> vars = ice_file->getVars();

    //for(auto it=vars.begin(); it != vars.end(); it++)
    //    printf("Var name %s\n", it->first.c_str());

    has_ice = new float[lat_size*lon_size];
    ice_var.getVar(has_ice);
    loaded_lat = new float[lat_size];
    loaded_lon = new float[lon_size];

    //printf("%i\n", ice_file->getVar("lat").getDim(0).getSize());

    ice_file->getVar("lat").getVar(loaded_lat);
    ice_file->getVar("lon").getVar(loaded_lon);

    //printf("Done\n");
}

//////////////////////////////////////////////////////////////////////////////////
//csv file////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
csv_file::csv_file(string f) {
    seperator = ',';
    fp = f;

    data = ifstream(fp.c_str());
}

csv_file::csv_file(string f, char s) {
    seperator = s;
    fp = f;

    data = ifstream(fp.c_str()); //TODO: This was commented out for some reason
}

bool csv_file::move_next() {
    getline(data, current_line);
    return (data.peek() != EOF); //TODO: Does this capture the last line?
}

bool csv_file::on_last_line() {
    return (data.peek() == EOF);
}

vector<double> csv_file::get_next() { 
    vector<string> cdata; 

    stringstream line_stream(current_line);
    string cell;
    while(getline(line_stream, cell, seperator)) {
        cdata.push_back(cell); 
    }

    return strings_to_dubs(cdata);
}

void csv_file::start() {
    if(data.is_open())
        data.close();

    data = ifstream(fp.c_str());
}

vector<double> csv_file::strings_to_dubs(vector<string> input_string) {
    vector<double> values;
    for(unsigned int i=0; i< input_string.size(); i++) {
        if(!input_string[i].compare("\r")) {
            break;
        } else if(input_string[i].empty()) {
            values.push_back(NAN);
        } else {
            try {
                values.push_back(stod(input_string[i])); //TODO: What if these isn't a double in here?
            } catch (exception &e) {
                values.push_back(NAN);
            }
        }
    }
    return values;
}

//CRU///////////////////////////////////////////////////////////////////////////
bool cru_file::move_next() {
    printf("Moving to next file\n");
    return (temp_data.move_next() and pre_data.move_next() and sunp_data.move_next() and reh_data.move_next() and elev_data.move_next());
}

//Open up files, if they're already open then close then first
void cru_file::start() { 
    temp_data.start(); 
    pre_data.start();
    sunp_data.start();
    reh_data.start();
    elev_data.start();
}

//Use the current position in all the files and grab all the items from them, append them to a vector and given them to the formatter, return the observation that the formatter gives
vector<double> cru_file::get_next() {
    vector<double> output = vector<double>(); 

    vector<double> elev_hold = vector<double>();
    vector<double> temp_hold = vector<double>();
    vector<double> prec_hold = vector<double>();
    vector<double> sunp_hold = vector<double>();
    vector<double> reh_hold = vector<double>();

    csv_file* files[5] = { &elev_data, &pre_data, &temp_data, &reh_data, &sunp_data };
    vector<double> tmp_holders[5] = { elev_hold, prec_hold, temp_hold, reh_hold, sunp_hold};

    //Get the current record for everything
    for(int i=0; i<5; i++) {
        tmp_holders[i] = remove_nans(files[i]->get_next());
    }

    //Parse everything together and add it to the input 
    for(int i=0; i<5; i++) {
        output.insert(output.end(), tmp_holders[i].begin(), tmp_holders[i].end()); 
    }

    return remove_nans(output);//Don't think this remove nans is needed anymore
}

cru_file::cru_format_return* cru_file::format_values(vector<double> input_dub) {
    //Set up the location as the first 3 entried
    cru_format_return *output = new cru_format_return(loc(input_dub[0], input_dub[1], input_dub[2]*1000)); //Get the elevation in metres

    //Put all the correct variables in the output holder 
    for(int i=0; i<12; i++) {
        output->values.set(0, output->values.get(0) + input_dub[i+5]);
        output->values.set(i+1, input_dub[i+31]);
        output->values.set(i+13, input_dub[i+45]);
        output->values.set(i+25, input_dub[i+59]);
    }

    return output;
}

//Given an input of values from the file, translate this to a climate reading
bartlein_csv_file::bart_format_return* bartlein_csv_file::format_values(vector<double> input_dub) {

    //Set the location and CO2 correction parameters 
    bart_format_return *output = new bart_format_return(loc(input_dub[2], input_dub[3], input_dub[4]), input_dub[6] == 1 ? true : false );

    //Set the appropriate values for the climate variabiles
    output->values.set(0, input_dub[17]); //This is alpha not moisture
    output->values.set(1, input_dub[13]);
    output->values.set(2, input_dub[11]);
    output->values.set(3, input_dub[7]);
    output->values.set(4, input_dub[9]);
    output->values.set(5, input_dub[15]);

    output->sds.set(0, input_dub[18]); //This is alpha not moisture
    output->sds.set(1, input_dub[14]);
    output->sds.set(2, input_dub[12]);
    output->sds.set(3, input_dub[8]);
    output->sds.set(4, input_dub[10]);
    output->sds.set(5, input_dub[16]);

    return output;
}

vector<double> bartlein_csv_file::strings_to_dubs(vector<string> input_string) {
    vector<double> values;

    for(unsigned int i=0; i< input_string.size(); i++) {
        if(!input_string[i].compare("\r")) {
            break;
        } else if(input_string[i].empty()) {
            values.push_back(NAN);
        } else {
            try {
                if(input_string[i] == "WU" || input_string[i] == "HAIBIN" ) 
                    values.push_back(1.0); 
                else
                    values.push_back(stod(input_string[i])); //TODO: What if there isn't a double in here?

            } catch (exception &e) {
                values.push_back(NAN);
            }
        }
    }
    return values;

}

//Given an input of values from the file, translate this to a climate reading
aus_csv_file::aus_format_return* aus_csv_file::format_values(vector<double> input_dub) {

    aus_format_return *output = new aus_format_return(loc(input_dub[4], input_dub[5], input_dub[3]));

    //Set the appropriate values for the climate varabiles
    output->values.set(0, input_dub[18]); 
    output->values.set(1, NAN);
    output->values.set(2, input_dub[16]);
    output->values.set(3, NAN);
    output->values.set(4, NAN);
    output->values.set(5, NAN);

    //These magic numbers are based on the average values from the bartlein pollen dataset
    output->sds.set(0, 0.0145); //TODO: Update this to 0.218
    output->sds.set(1, NAN);
    output->sds.set(2, 1.82);
    output->sds.set(3, NAN);
    output->sds.set(4, NAN);
    output->sds.set(5, NAN);

    return output;
}

////////////////////////////////////////////////////////////////////
//Utility functions/////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


vector<double> remove_nans(vector<double> input_dubs) {
    vector<double>::iterator it = input_dubs.begin();

    while(it!=input_dubs.end()) {
        if(isnan(*it)) {
            it = input_dubs.erase(it);
        } else {
            ++it;
        }
    }
    return input_dubs;
}
