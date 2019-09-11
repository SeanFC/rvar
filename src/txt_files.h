#ifndef txt_files_h
#define txt_files_h 

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

#include <netcdf>

#include "readings.h"

using namespace netCDF;
using namespace netCDF::exceptions;

class Netcdf_PMIP_Info : public Seasonal_Modern_reading {
    public:
        loc location;
        //Seasonal_Modern_reading values;
        Vec_f values_held = Vec_f(4); //This is double because I can't be bothered to write the int case properly

        virtual Netcdf_PMIP_Info& ignore_add(const Seasonal_Modern_reading& rhs) {
            if(values_held.get(0)) set_MAP  (get_MAP()  + rhs.get_MAP()  ); else set_MAP  ( rhs.get_MAP()  );
            if(values_held.get(1)) set_temp (get_temp() + rhs.get_temp() ); else set_temp ( rhs.get_temp() );
            if(values_held.get(2)) set_RH   (get_RH()   + rhs.get_RH()   ); else set_RH   ( rhs.get_RH()   );
            if(values_held.get(3)) set_sunf (get_sunf() + rhs.get_sunf() ); else set_sunf ( rhs.get_sunf() );

            return *this;
        }

        virtual Netcdf_PMIP_Info& operator+=(const Seasonal_Modern_reading& sq) { this->Seasonal_Modern_reading::operator+=(sq); return *this; }
        virtual Netcdf_PMIP_Info& operator-=(const Seasonal_Modern_reading& sq) { this->Seasonal_Modern_reading::operator-=(sq); return *this; }
        virtual Netcdf_PMIP_Info& operator*=(const double dub) { this->Seasonal_Modern_reading::operator*=(dub); return *this; }
        virtual Netcdf_PMIP_Info& operator/=(const double dub) { this->Seasonal_Modern_reading::operator/=(dub); return *this; }

        virtual Netcdf_PMIP_Info& piecewise_square() { this->Seasonal_Modern_reading::piecewise_square(); return *this; }
};

inline Netcdf_PMIP_Info ignore_add(Netcdf_PMIP_Info lhs, const Seasonal_Modern_reading& rhs) {
    lhs.ignore_add(rhs);
    return lhs;
}

inline Netcdf_PMIP_Info operator+(Netcdf_PMIP_Info lhs, const Seasonal_Modern_reading& rhs) {
    lhs += rhs;
    return lhs;
}

inline Netcdf_PMIP_Info operator-(Netcdf_PMIP_Info lhs, const Seasonal_Modern_reading& rhs) {
    lhs -= rhs;
    return lhs;
}

inline Netcdf_PMIP_Info operator*(Netcdf_PMIP_Info lhs, const double rhs) {
    lhs *= rhs;
    return lhs;
}

inline Netcdf_PMIP_Info operator/(Netcdf_PMIP_Info lhs, const double rhs) {
    lhs /= rhs;
    return lhs;
}

inline bool does_file_exist(string fp) { return ifstream(fp).good(); }

//All the meta data of a particular netCDF file
struct PMIP_nc_file_meta {
    string var_name;
    string experiment_name;
    string grid_resolution;
    string project_name;
    string model_name;
    string file_name;
    string path_location;

    PMIP_nc_file_meta(string fn, string pl) : file_name(fn), path_location(pl) {
        std::stringstream ss(file_name);
        std::string item;
        for(int i=0; i<4; i++) {
            std::getline(ss, item, '_');
            switch(i) {
                case 0:
                    var_name = item; break;
                case 1:
                    experiment_name = item; break;
                case 2:
                    grid_resolution = item; break;
                case 3:
                    project_name = item; break;
            }
        }

        std::getline(ss, item);
        model_name = item;
        model_name.resize(model_name.size() - 3);
    };

    void print() {
        printf(" %s %s %s %s %s\n", var_name.c_str(), experiment_name.c_str(), grid_resolution.c_str(), project_name.c_str(), model_name.c_str());
    };

    bool experiment_model_match(PMIP_nc_file_meta other) {
        return (experiment_name == other.experiment_name) && (project_name == other.project_name) && (model_name == other.model_name);
    }
};


//A PMIP netCDF file
class netcdf_PMIP {
    private:
        NcFile *tas_file = NULL;
        NcFile *cl_file = NULL;
        NcFile *pr_file = NULL;
        NcFile *hur_file = NULL;

        //Find the sizes of all the dimensions 
        int time_size; //TODO: Taken out the const restriction here
        int lat_size; 
        int lon_size;

        //Variably sized Multi dimensional arrays in C++ and C are absolutely ridiculous, i'm just going to do it this way
        float *loaded_tas;
        float *loaded_cl;  
        float *loaded_pr;
        float *loaded_hur;
        float *loaded_lat;
        float *loaded_lon;
        float *loaded_time;
        
        float days_in_given_months[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}; 

        void open_up_files(string tas_fp, string cl_fp, string pr_fp, string hur_fp);

    public:
        string experiment_name;
        string grid_resolution;
        string project_name;
        string model_name;

        netcdf_PMIP(string tas_fp, string cl_fp, string pr_fp, string hur_fp) {
            open_up_files(tas_fp, cl_fp, pr_fp, hur_fp);
        }; 
        netcdf_PMIP(vector<PMIP_nc_file_meta> file_metas);

        ~netcdf_PMIP() { 
            delete []loaded_tas;
            delete []loaded_pr;
            delete []loaded_cl;
            delete []loaded_hur;
            delete []loaded_lat;
            delete []loaded_lon;

            if(tas_file) delete tas_file;
            if(cl_file) delete cl_file;
            if(pr_file) delete pr_file;
            if(hur_file) delete hur_file;
        } 

        void load_all_readings();

        Netcdf_PMIP_Info* get_lat_lon_record(int lat_index, int lon_index);

        void print() { printf("%s, %s, %s, %s\n", experiment_name.c_str(), grid_resolution.c_str(), project_name.c_str(), model_name.c_str()); }
};

class Netcdf_ice {
    private:
        string file_path;

        NcFile *ice_file = NULL;

        float *has_ice;
        float *loaded_lat;
        float *loaded_lon;
        int lat_size;
        int lon_size;

    public:
        Netcdf_ice(string fp) : file_path(fp) {
            if(does_file_exist(file_path))
                ice_file = new NcFile(file_path.c_str(), NcFile::read);
        }

        ~Netcdf_ice() {
            delete []has_ice;

            if(ice_file) delete ice_file;
        }

        void load_in_info();

        double get_lat(int lat_index) { return loaded_lat[lat_index]; }
        double get_lon(int lon_index) { return loaded_lon[lon_index]; }
        double get_has_ice(int lat_index, int lon_index) { return has_ice[lat_index*lon_size + lon_index]; }

        double get_lat_size() { return lat_size; }
        double get_lon_size() { return lon_size; }
};

class csv_file {
        string fp;
        ifstream data;

        char seperator;
        string current_line;

    public:
        csv_file(string);
        csv_file(string, char);
        ~csv_file() { data.close(); } //TODO: Only do this if data is open

        virtual void start();
        virtual bool move_next();
        virtual bool on_last_line();
        virtual void start(int begin_at) { for(int i=0; i<begin_at; i++) { move_next(); } };
        virtual vector<double> get_next();
        virtual vector<double> strings_to_dubs(vector<string>);

        string get_fp()  { return fp; }
        char get_sep()  { return seperator; }

        int current_pos() { return data.tellg(); }
};

//CRU space seperated files
class cru_file {
    private:
        //Actual Files
        csv_file temp_data;
        csv_file pre_data;
        csv_file sunp_data;
        csv_file reh_data;
        csv_file elev_data;

    public:
        bool move_next(); //Move to the next item, don't get it yet
        void start(); //Move to the start and get the first item

        vector<double> get_next(); //Get the next item

        cru_file(string t_fp, string p_fp, string s_fp, string r_fp, string e_fp) : temp_data(t_fp, ' '), pre_data(p_fp, ' '), sunp_data(s_fp, ' '), reh_data(r_fp, ' '), elev_data(e_fp, ' ') {}; //Constructor

        struct cru_format_return {
            loc location;
            Vec_f values;

            cru_format_return(loc lo) : location(lo), values(37) {};
        };

        cru_format_return* format_values(vector<double> input_dub);
        cru_format_return* get_next_reading() { return format_values(get_next()); }
};

class bartlein_csv_file : public csv_file {
    public:
        bartlein_csv_file(string fp) : csv_file(fp) {} 

        void start() { csv_file::start(); move_next(); } //We skip the first line

        struct bart_format_return {
            loc location;
            Vec_f values;
            Vec_f sds; 
            bool co2_corrected;

            bart_format_return(loc lo, bool corrected) : location(lo), values(6), sds(6), co2_corrected(corrected) {};
        };

        bart_format_return* format_values(vector<double> input_dub);
        bart_format_return* get_next_reading() { return format_values(get_next()); }


        virtual vector<double> strings_to_dubs(vector<string>);
};

class aus_csv_file : public csv_file {
    public:
        aus_csv_file(string fp) : csv_file(fp) {}

        void start() { csv_file::start(); move_next(); } //We skip the first line

        struct aus_format_return {
            loc location;
            Vec_f values;
            Vec_f sds; 

            aus_format_return(loc lo) : location(lo), values(6), sds(6) {};
        };

        aus_format_return* format_values(vector<double> input_dub);
        aus_format_return* get_next_reading() { return format_values(get_next()); }
};

/////////////////////////////////////////////////////////////////////////////////////
//Utility Fuctions///////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
vector<double> remove_nans(vector<double> input_dubs);

#endif
