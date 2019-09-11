#ifndef geo_observations_h 
#define geo_observations_h 

#include <thread>
#include <mutex>

#include "wang_p_model_correction.h"
#include "readings.h"
#include "observation_function.h"
#include "lin_maths.h"
#include "location.h"


class MI_alpha_calc { 
    protected:
        orbital orb;
        loc location;

        static am_space_reading& process_averages(am_space_reading &output_obs, const Vec_f &temps, float MAP);
        static am_space_reading& process_nondim_averages(am_space_reading &output_obs, const Vec_f &temps, float MAP);

    public:
        MI_alpha_calc(const orbital &o, loc l) : orb(o), location(l) {};
        am_space_reading sto_no_CO2(const Seasonal_Modern_reading&); 

        static am_space_reading& dimensionalise_obs(am_space_reading& input_obs);
        static am_sd_reading& dimensionalise_obs_with_sd(am_sd_reading& input_obs);
        static am_space_reading& nondimensionalise_obs(am_space_reading& input_obs);
        static am_sd_reading& nondimensionalise_obs_with_sd(am_sd_reading& input_obs);

        static Seasonal_space_reading& dimensionalise_state(Seasonal_space_reading& input_state);
        static Seasonal_sd_reading& dimensionalise_state_with_sd(Seasonal_sd_reading& input_state);
        static Seasonal_space_reading& nondimensionalise_state(Seasonal_space_reading& input_state);
        static Seasonal_sd_reading& nondimensionalise_state_with_sd(Seasonal_sd_reading& input_state);

        loc get_position() { return location; }
};

class Wang_P_Budyko : public MI_alpha_calc, public sto {
    protected:
        Seasonal_Modern_reading &sta_modern_climate;
        Seasonal_Modern_reading &sta_model_climate;
        Moisture_Corrector mc;
        bool use_CO2;

    public:
        Wang_P_Budyko(Seasonal_Modern_reading &sta_modern_day_climate, Seasonal_Modern_reading &model_climate, const orbital &orbital_params, loc locatio, bool co2) : MI_alpha_calc(orbital_params, locatio), sta_modern_climate(sta_modern_day_climate), sta_model_climate(model_climate),
        mc(
                sta_modern_climate.get_MAP(),//precip_modern
                sta_modern_climate.get_RH(),//RH_modern 
                sta_model_climate.get_RH(),//RH_past
                sta_modern_climate.get_temp(), //temp_modern
                sta_modern_climate.get_sunf(), //sunf_modern
                location, orbital_params, location, FIF_orbital//loc_past_in, orbital orb_past_in, loc loc_mod_in, orbital orb_mod_in
          ),
        use_CO2(co2)
        {};

        virtual ~Wang_P_Budyko() {}

        //Vec_f full_uncorrection(const Seasonal_Modern_reading &input_state);

        //TODO:These 4 functions should be private (or at least protected), only sto should be public 
        am_space_reading PB_obs_func(const Seasonal_space_reading&);
        am_space_reading PB_obs_func_no_CO2(const Seasonal_space_reading&);

        virtual am_space_reading PB_obs_func_derivative(const Seasonal_space_reading &point, const am_space_reading &reference_point, int index);
        virtual am_space_reading PB_obs_func_no_CO2_derivative(const Seasonal_space_reading &point, const am_space_reading &reference_point, int index);

        virtual am_space_reading state_to_observation(const Seasonal_space_reading &state) { 
            if(use_CO2)
                return PB_obs_func(state); 
            else
                return PB_obs_func_no_CO2(state); 
        }

        virtual Vec_f state_to_observation(const Vec_f &state) { 
            Seasonal_space_reading temp_reading = Seasonal_space_reading();
            temp_reading.set(state);
            return state_to_observation(temp_reading); 
        }
        virtual Vec_f derivative(const Vec_f &point, const Vec_f &reference_point, int index) {
            Seasonal_space_reading temp_reading = Seasonal_space_reading();
            temp_reading.set(point);
            am_space_reading temp_ref = am_space_reading();
            temp_ref.set(reference_point);

            if(use_CO2)
                return PB_obs_func_derivative(temp_reading, temp_ref, index); 
            else
                return PB_obs_func_no_CO2_derivative(temp_reading, temp_ref, index); 
        }

        Seasonal_Modern_reading& get_modern_climate() { return sta_modern_climate; }
        Seasonal_Modern_reading& get_model_climate() { return sta_model_climate; }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
//Grid STOs//////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
class grid_sto {
    protected:
        int state_dim;
        int obs_dim;
        float distance_scale_factor;
        Mat_f location_matrix;

    public:
        grid_sto(int state_d, int obs_d, float dist_factor, int grid_amount, int obs_amount) : state_dim(state_d), obs_dim(obs_d), distance_scale_factor(dist_factor), location_matrix(grid_amount, obs_amount) { }

        int get_state_dim() { return state_dim; }
        int get_obs_dim() { return obs_dim; }
};

class Wang_P_Budyko_grid : public sto, grid_sto {
    protected:
        vector<Wang_P_Budyko*> &mini_stos;

    public:
        Wang_P_Budyko_grid(vector<Wang_P_Budyko*> &mini_sto, vector<loc> &location_grid, int state_d, int obs_d, float dist_factor, int number_of_threads);

        virtual ~Wang_P_Budyko_grid() {}

        //Vec_f state_to_observation_on_grid(vector<Wang_P_Budyko*> &full_stos, Vec_f state);
        virtual Vec_f state_to_observation(const Vec_f &state);
        virtual Vec_f derivative(const Vec_f &point, const Vec_f &reference_run, int index);
        virtual Block_f block_jacobian(const Vec_f &point);
        virtual Block_Diag_f block_diag_jacobian(const Vec_f &point);
};

#endif
