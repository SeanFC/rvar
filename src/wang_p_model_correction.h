#ifndef wang_p_model_correction_h
#define wang_p_model_correction_h

#include <math.h> //For exp
#include <complex> //For complex numbers 
#include <limits> //For infinity (or nan)

#include "orbital.h" //For orbital parameters
#include "location.h" //For location
#include "lin_maths.h" //For taking derivatives

const float I_sc = 1360.8; //W m−2
const float latent_heat_of_vaporisation_of_water = 2.45;//MJ kg^{-1} //lambda
const float psychrometer_constant = 0.067; //kPa K^{-1} at sea level
const float celcius_to_kelvin_addition = 273.15; // Kelvin - this = Celcius
const float delta_H_Gamma = 37830; //J mol^{-1}
const float universal_gas_constant = 8.314; //J mol^{-1} K^{-1}
const float activation_energy_carbon = 79430; //J mol^{-1} 
const float activation_energy_oxgyen = 36380; //J mol^{-1} 
const float atmospheric_concentration_of_oxygen = 210; // permil
const float ratio_of_carboxylation_and_transpiration_25 = 240; //The ratio between carboxylation and transpiration costs at 25C from wang2016towards
const float michaelis_menton_coeff_scale_K = 0.03997; //kPa
const float standard_atmospheric_pressure = 101.325; // kPa, Allen 1973
const int number_of_days_in_year = 365; //Number of days in a year TODO:Assuming no leap years here
const int number_of_months_in_year = 12; 
const float entrainment_factor = 0.26;
const float supply_rate_constant = 1.05; // mm / h 
const float soil_capacity = 150; // mm 
const float pir = M_PI/180.0;
const float viscosity_of_water_ref_num = 138.0; //TODO:This is 183 in prentice2017constraining but could be 138 in other sources
const float viscosity_of_water_at_25C_magic_number = 3.6216047455510463; //580/(25 + celcius_to_kelvin_addition - viscosity_of_water_ref_num)
const float temp_scaler = 5;

const int days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

namespace Wang_Dimension {
    float dim_precip(float non_P);
    float dim_precip_derivative(float non_P);
    float nondim_precip(float P);
    float nondim_precip_derivative(float P);

    float nondim_temp(float temp);
    float dim_temp(float non_temp);
}

class Ground_Radiation {
    private: 
        //Given values
        loc location;
        const orbital &orb;

        //Values calculated for the year
        float xlam;
        float xee;
        float xec;
        float xse;

        //float daytime_ground_radiation[number_of_months_in_year];
        //float nighttime_ground_radiation[number_of_months_in_year];
        float distance_factor[number_of_months_in_year];
        float monthly_hour_sunset_angle[number_of_months_in_year];
        float monthly_ru[number_of_months_in_year];
        float monthly_rw[number_of_months_in_year];
        float monthly_rv[number_of_months_in_year];
        float monthly_hn[number_of_months_in_year];
    
    public:

        Ground_Radiation(loc l, const orbital &o);

        float* get_monthly_ground_radiation_arrays(float *T, float *s_f, float *output);
        float get_single_month_radiation(float temp, float s_f, int month_of_year);
        //float get_average_daily_radiation_for_month(int month);
        //float get_radiation_for_month(int month);
};

float energy_to_water_conversion_with_CO2(
        float temp_past, 
        float relative_humidity_past,
        float c_a_past,
        float temp_modern,
        float relative_humidity_modern,
        float c_a_modern);

float e_s_wrt_temp(
        float temp_past, 
        float relative_humidity_past,
        float c_a_past,
        float temp_modern,
        float relative_humidity_modern,
        float c_a_modern);

float VPD_wang(float temp_past, float relative_humidity_past, float c_a_past, float temp_modern, float c_a_modern, float water_lost);

float water_lost_by_CO2_uptake_nd(float temp, float relative_humidity, float c_a);

//kPa
float saturation_vapour_pressure(float temp);
float saturation_vapour_pressure_nd(float temp);

//kPa //Differentiation of saturaiton_vapour_pressure
float saturation_vapour_pressure_wrt_temp(float temp);
float saturation_vapour_pressure_wrt_temp_nd(float temp_nd);

float wang_compensation_point(float temp);
float wang_compensation_point_nd(float temp);
float wang_compensation_point_wrt_T(float temp);
float wang_compensation_point_wrt_T_nd(float temp_nd);

//Change temperature from degrees Celsius to Kelvin units
float T_C_to_K(float temp);

//The VPD as calculated from temperature and rh with no [CO_2] change calculation needed
float real_VPD(float temp, float relative_humidity);
float real_VPD_nd(float temp, float relative_humidity);

float wang_K(float temp);
float wang_K_nd(float temp_nd);

float wang_K_wrt_T(float temp);
float wang_K_wrt_T_nd(float temp_nd);

//The viscosity of water relative to its value at 25C
float wang_eta(float temp);
float wang_eta_nd(float temp);

float wang_eta_wrt_T(float temp_nd);
float wang_eta_wrt_T_nd(float temp_nd);

float colin_P_model_factor(float eta, float K, float Gamma);
float colin_P_model_factor_wrt_T(float xi, float eta, float K, float Gamma, float eta_wrt_T, float K_wrt_T, float Gamma_wrt_T);
 
//Finds the moisture index as if it were constructed without accounting for changes in [CO_2]
class Moisture_Corrector {
    private:
        float precip_modern;
        const Vec_f relative_humidity_modern, relative_humidity_past, temp_modern, s_f_modern; 
        const loc &loc_past;
        orbital orb_past; //This should be const & but it keeps getting deleted somewhere

        Ground_Radiation radiation_palaeo;
        Vec_f modern_e;
        Vec_f modern_e_wrt_T;

    public:
        Moisture_Corrector(
                float precip_modern_in, 
                const Vec_f relative_humidity_modern_in, 
                const Vec_f relative_humidity_past_in, 
                const Vec_f temp_modern_in, 
                const Vec_f s_f_modern_in,
                const loc &loc_past_in, 
                const orbital &orb_past_in, 
                const loc &loc_mod_in, 
                const orbital &orb_mod_in);    

        float uncorrect_mi(
                float annual_precip, 
                const Vec_f &temp_past
                );

        float uncorrect_mi_wrt_T(
                float mi,
                float annual_precip, 
                float temp_past,
                int month_index
                );

        float true_mi(float annual_precip, const Vec_f &temperature, const Vec_f &sunf, const loc &location) { 
            return true_mi(annual_precip, temperature, sunf, orb_past, location);
        }
        static float true_mi(float annual_precip, const Vec_f &temperature, const Vec_f &sunf, const orbital &orb, const loc &location); 

        orbital get_past_orbital() { return orb_past; }
        Vec_f get_sunf() { return s_f_modern; }
};

float e_s_wrt_temp_alg(
        float temp_past, 
        float relative_humidity_past,
        float c_a_past,
        float mod_e,
        float mod_e_wrt_T
        );

inline float budyko_relationship(float m) {
    if(m >=0)  {
        return 1.26*(1 + m - cbrt(1 + pow(m,3)));
    }else {
        m = -m;
        return -1.26*(1 + m - cbrt(1 + pow(m,3)));
    }
}

inline float budyko_relationship_derivative(float m) {
    if(m >=0)  {
        return 1.26*(1 - pow(m, 2) / pow(cbrt(1 + pow(m, 3)), 2 ));
    }else {
        m = -m;
        return 1.26*(1 - pow(m, 2) / pow(cbrt(1 + pow(m, 3)), 2 ));
    }
}

float nondim_temp(float temp);
float dim_temp(float temp);

float true_mi_simple(float precip, float* temps, float* sunf, float* location, int orb_choice);

#endif
