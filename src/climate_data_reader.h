#ifndef climate_data_reader_h 
#define climate_data_reader_h 

#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h> // For system calls
#include <math.h> // For log2

#include "settings_interpreter.h"
#include "readings.h"
#include "geo_observations.h"
#include "accessors.h"
#include "rvar.h"
#include "cblas.h"
#include "observation_function.h"

using namespace std;

// Given input of the observations and background, this class can compuate the analysis along with various properties of the analyis
class Seasonal_correlated_minimiser : public ccm {
    private:
        vector<Wang_P_Budyko*> full_grid_stos;
        int inner_state_dim, inner_obs_dim, amount_of_backgrounds, amount_of_observations, view_size;
        vector<loc> state_location_grid;
        vector<Wang_P_Budyko*> s2o;
        Diagonal_f bac_sd_eye;
        Kroneker_f C; //Background correlation matrix
        Diagonal_f R; //Observation covariance matrix
        Wang_P_Budyko_grid *no_CO2_observer = NULL;

    public:
        Seasonal_correlated_minimiser(
                vector<OBS_accessor::OBS_obs_reading*>&, 
                vector<Seasonal_sd_reading*>&, 
                vector<Wang_P_Budyko*>&, 
                vector<Wang_P_Budyko*>&
                );

        virtual ~Seasonal_correlated_minimiser();

        void set_scales(double, double, double, int);

        struct os_map {
            vector<am_sd_reading> obs;
            vector<Seasonal_sd_reading> states;
            Vec_f x;
            Vec_f w;
            double cost;
            Vec_f gradient;
            Mat_f H_a;
            vector<am_space_reading> attrib;
            Mat_f resolution;
            double cond_number;

            os_map(int state_size, int obs_size) : x(state_size), w(state_size), gradient(state_size), H_a(state_size, obs_size), resolution(1, 1) {} //TODO:Resolution should be dealt with, don't want to be computing it every time really. Change this back to 1,1 for the real plots
        };

        void calculate_uncertainty_shifted(const Vec_f &w, const Block_Diag_f &G_a, const Block_Diag_f &G_b, Vec_f& ana_err, Vec_f& bac_err);

        os_map parse_current_result();
        Mat_f get_hessian(Vec_f);
        double get_condition_number();
        Mat_f get_resolution_matrix();
};

////////////////////////////////////////////////////////////////////////////////////
//Full Data Assimilators////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
// Given raw PMIP and observation readings, compute the analysis and save it, along with several useful properies of the system. This class also takes intout account how the CO_2 should be handled and other user specificed options. 
class pmip_seasonal_assimilator {
    private:
        Single_PMIP_accessor &background;
        OBS_accessor &obs;
        CRU_accessor &modern;
        Seasonal_correlated_minimiser *scm = NULL;
        int grid_degrade;

        vector<OBS_accessor::OBS_obs_reading*> obs_readings;
        vector<PMIP_accessor::PMIP_seasonal_reading*> all_back_readings;
        vector<int> back_reading_position_ids;
        vector<Seasonal_sd_reading*> back_readings;
        vector<CRU_accessor::CRU_seasonal_reading*> modern_readings;
        vector<Wang_P_Budyko*> stos;
        vector<Wang_P_Budyko*> state_grid_stos;

        double lat_max = 21.25; 
        double lat_min = 20.8; 
        double lon_max = 18.6; 
        double lon_min = 18.25; 

        int number_of_threads = 0;

        bool use_CO2;

    public:
        pmip_seasonal_assimilator(
                Single_PMIP_accessor &back,
                OBS_accessor &o,
                CRU_accessor &mod,
                int grid_degrade,
                double lat_max,
                double lat_min,
                double lon_max,
                double lon_min,
                int number_of_threads, 
                int time_period_id,
                bool use_CO2,
                bool exclude_PSM_obs
                ); 

        ~pmip_seasonal_assimilator();

        Seasonal_correlated_minimiser::os_map static_analysis(double rh_scale, double temp_scale, double grid_scale);
        Seasonal_correlated_minimiser::os_map find_analysis(double, double, double);

        void save_analysis(ANA_experiment &analysis, Seasonal_correlated_minimiser::os_map computed_ana);

        void preload_minimisation_matricies();
        double get_condition_number(double rh_scale, double temp_scale, double grid_scale);
        Mat_f get_resolution_matrix(double rh_scale, double temp_scale, double grid_scale);
};

#endif
