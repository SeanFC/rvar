#include "wang_p_model_correction.h"

// TODO: This is no longer supported 
float energy_to_water_conversion_with_CO2(
        float temp_past, 
        float relative_humidity_past,
        float c_a_past,
        float temp_modern,
        float relative_humidity_modern,
        float c_a_modern) {

    //convert c_a's into kPa
    c_a_past = standard_atmospheric_pressure*c_a_past*1e-6; 
    c_a_modern = standard_atmospheric_pressure*c_a_modern*1e-6;

    float e_s_diff_temp = e_s_wrt_temp(
            temp_past, 
            relative_humidity_past,
            c_a_past,
            temp_modern,
            relative_humidity_modern,
            c_a_modern);

    float psychrometer_ratio = isinf(e_s_diff_temp) ? 1 : e_s_diff_temp / (e_s_diff_temp + psychrometer_constant);

    //printf("%f %f %f %f %f\n", e_s_diff_temp, temp_past, temp_modern, relative_humidity_past, relative_humidity_modern);
    return psychrometer_ratio / latent_heat_of_vaporisation_of_water;
}

float e_s_wrt_temp(
        float temp_past, 
        float relative_humidity_past,
        float c_a_past,
        float temp_modern,
        float relative_humidity_modern,
        float c_a_modern) {

    float water_lost = water_lost_by_CO2_uptake_nd(
            temp_past, 
            relative_humidity_past, 
            c_a_past
            );

    float delta = 1;//TODO: Configure this 1e-6;
    float VPD_change = 0;
    if(isinf(water_lost)) {
        VPD_change = (sqrt(ratio_of_carboxylation_and_transpiration_25* ((wang_K(temp_modern+delta) +wang_compensation_point(temp_modern+delta)))/(1.6 * wang_eta(temp_modern+delta)))
                -
                sqrt(ratio_of_carboxylation_and_transpiration_25* ((wang_K(temp_modern+delta) +wang_compensation_point(temp_modern+delta)))/(1.6 * wang_eta(temp_modern+delta))));
    } else {
        float pre_VPD = VPD_wang(
                temp_past, 
                relative_humidity_past,
                c_a_past, 
                temp_modern+delta, 
                c_a_modern,
                water_lost
                );
        float post_VPD = VPD_wang(
                temp_past, 
                relative_humidity_past,
                c_a_past, 
                temp_modern-delta, 
                c_a_modern,
                water_lost
                );
        VPD_change = pre_VPD - post_VPD;
        if(isinf(pre_VPD) || isinf(post_VPD))  
            return 0;
    }
    return VPD_change/(2*delta*(1- relative_humidity_modern/100));
}

float VPD_wang(float temp_past, float relative_humidity_past, float c_a_past, float temp_modern, float c_a_modern, float water_lost) {
    float comp_point_mod = wang_compensation_point(temp_modern);

    float d_hat_inv = sqrt(ratio_of_carboxylation_and_transpiration_25* ((wang_K(temp_modern) + comp_point_mod))/(1.6 * wang_eta(temp_modern)));
       
    float delta = pow(d_hat_inv,2) + 4*water_lost/1.6*(c_a_modern - comp_point_mod);
    float D = 0.25*(pow(d_hat_inv,2) + d_hat_inv*sqrt(delta) + delta);

    return D; 
}

//nd here is non-dim
float water_lost_by_CO2_uptake_nd(float temp_nd, float relative_humidity_nd, float c_a){
    float comp_point = wang_compensation_point_nd(temp_nd);

    float D = real_VPD_nd(temp_nd, relative_humidity_nd);
    float t_inv = sqrt(ratio_of_carboxylation_and_transpiration_25*(wang_K_nd(temp_nd) + comp_point)/(1.6 * wang_eta_nd(temp_nd) * D));
    float e = 1.6*D * (t_inv/sqrt(D)+1)/(c_a - comp_point);

    if(c_a < comp_point) {
        //float temp = celcius_to_kelvin_addition * (exp(temp_nd) - 1);
        //printf("Comp_point broken: temp %f temp_nd %f rh_nd %f c_a %f comp_point %f D %f t_inv %f e %f\n", temp, temp_nd, relative_humidity_nd, c_a, comp_point, D, t_inv, e);
        return std::numeric_limits<float>::infinity();
    }

    return e; 
}

//kPa
float saturation_vapour_pressure(float temp) {
    return 0.6108 * exp(17.27 * temp/(temp+237.3)); 
}

//kPa
float saturation_vapour_pressure_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    return 0.6108 * exp(17.27 * temp/(temp+237.3)); 
}

//kPa K^{-1} //Differentiation of saturation_vapour_pressure
float saturation_vapour_pressure_wrt_temp(float temp) {
    return 4098.171*saturation_vapour_pressure(temp)/pow(temp + 237.3, 2);  
}

float saturation_vapour_pressure_wrt_temp_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    return 4098.171*saturation_vapour_pressure_nd(temp_nd)/pow(temp + 237.3, 2); 
}

//kPa
float wang_compensation_point(float temp) {
    return 0.00422 * exp(delta_H_Gamma / universal_gas_constant *( 1.0/298.15 - 1.0/T_C_to_K(temp))); //4.220  42.75 for \mu mol mol^{-1}, 0.00422 for compensation_point
}

//kPa
float wang_compensation_point_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    return 0.00422 * exp(delta_H_Gamma / universal_gas_constant *( 1.0/298.15 - 1.0/T_C_to_K(temp))); 
}

float wang_compensation_point_wrt_T(float temp) {
    return wang_compensation_point(temp) * delta_H_Gamma/universal_gas_constant / pow(T_C_to_K(temp), 2); 
}

float wang_compensation_point_wrt_T_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    return wang_compensation_point(temp) * delta_H_Gamma/universal_gas_constant / pow(T_C_to_K(temp), 2); 
}

//Change temperature from degrees Celsius to Kelvin units
float T_C_to_K(float temp) {
    return temp + celcius_to_kelvin_addition;
}

float real_VPD(float temp, float relative_humidity) {
    return (1 - relative_humidity/100) * saturation_vapour_pressure(temp);
}

float real_VPD_nd(float temp, float relative_humidity) {
    float rh = atan(relative_humidity)*M_PI + 50; 
    return (1 - rh/100) * saturation_vapour_pressure_nd(temp);
}

//kPa
float wang_K(float temp) {
    float farquhar_pre_calc = 1/universal_gas_constant*(1.0/298.15 - 1.0/(float)(T_C_to_K(temp))); 
    /*printf("K: %f %f %f %f\n", 
            temp, 
            404.9 * exp(activation_energy_carbon * farquhar_pre_calc), 
            278.4 * exp(activation_energy_oxgyen * farquhar_pre_calc), 
            404.9 * exp(activation_energy_carbon * farquhar_pre_calc) * ( 1.0 + atmospheric_concentration_of_oxygen/(278.4 * exp(activation_energy_oxgyen * farquhar_pre_calc))));*/
    return michaelis_menton_coeff_scale_K * exp(activation_energy_carbon * farquhar_pre_calc) * ( 1.0 + atmospheric_concentration_of_oxygen/(278.4 * exp(activation_energy_oxgyen * farquhar_pre_calc))); //39.97 404.9 
}

float wang_K_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    float farquhar_pre_calc = 1/universal_gas_constant*(1.0/298.15 - 1.0/(float)(T_C_to_K(temp))); 
    return michaelis_menton_coeff_scale_K * exp(activation_energy_carbon * farquhar_pre_calc) * ( 1.0 + atmospheric_concentration_of_oxygen/(278.4 * exp(activation_energy_oxgyen * farquhar_pre_calc))); //39.97 404.9 
}

float wang_K_wrt_T(float temp) {
    float farquhar_pre_calc = 1/universal_gas_constant*(1.0/298.15 - 1.0/(float)(T_C_to_K(temp))); 

    return michaelis_menton_coeff_scale_K * exp(activation_energy_carbon * farquhar_pre_calc) / (universal_gas_constant * pow(temp + celcius_to_kelvin_addition, 2)) * (activation_energy_carbon + atmospheric_concentration_of_oxygen/(278.4 * exp(activation_energy_oxgyen * farquhar_pre_calc))*(activation_energy_carbon - activation_energy_oxgyen));
}

float wang_K_wrt_T_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    float farquhar_pre_calc = 1/universal_gas_constant*(1.0/298.15 - 1.0/(float)(T_C_to_K(temp))); 

    return michaelis_menton_coeff_scale_K * exp(activation_energy_carbon * farquhar_pre_calc) / (universal_gas_constant * pow(temp + celcius_to_kelvin_addition, 2)) * (activation_energy_carbon + atmospheric_concentration_of_oxygen/(278.4 * exp(activation_energy_oxgyen * farquhar_pre_calc))*(activation_energy_carbon - activation_energy_oxgyen));
}

float colin_P_model_factor(float eta, float K, float Gamma) {
    return sqrt(1.6 * eta / (ratio_of_carboxylation_and_transpiration_25 * (K + Gamma)));
}

float colin_P_model_factor_wrt_T(float xi, float eta, float K, float Gamma, float eta_wrt_T, float K_wrt_T, float Gamma_wrt_T) {
    return  0.8/ ratio_of_carboxylation_and_transpiration_25 / xi * (eta_wrt_T * (K + Gamma) - eta * (K_wrt_T + Gamma_wrt_T)) / pow(K + Gamma, 2); 
}

//The viscosity of water relative to its value at 25C
float wang_eta(float temp) {
    return exp(580.0/(float)(T_C_to_K(temp) - viscosity_of_water_ref_num) - viscosity_of_water_at_25C_magic_number); 
}

//The viscosity of water relative to its value at 25C differentiated wrt T
float wang_eta_wrt_T(float temp) {
    return -580.0/pow(T_C_to_K(temp) - viscosity_of_water_ref_num,2) * wang_eta(temp);
}

//The viscosity of water relative to its value at 25C
float wang_eta_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    
    //if(temp < -100.0)
    //    printf("Eta info dump: %f %f %f %f\n", temp, temp_nd, exp(580.0/(float)(T_C_to_K(temp) - viscosity_of_water_ref_num) - viscosity_of_water_at_25C_magic_number), 580.0/(float)(T_C_to_K(temp) - viscosity_of_water_ref_num) - viscosity_of_water_at_25C_magic_number);
    
    return exp(580.0/(float)(T_C_to_K(temp) - viscosity_of_water_ref_num) - viscosity_of_water_at_25C_magic_number); 
}

float wang_eta_wrt_T_nd(float temp_nd) {
    float temp = Wang_Dimension::dim_temp(temp_nd);
    
    float to_find_eta = 0.024258 * exp(580.0/(float)(T_C_to_K(temp) - viscosity_of_water_ref_num));
    float eta_25 = 0.024258 * exp(580.0/(float)(T_C_to_K(25) - viscosity_of_water_ref_num)); //TODO: Export as const expr
    return to_find_eta/eta_25 * ( -580.0 / pow( T_C_to_K(temp) - viscosity_of_water_ref_num, 2));
}

//Implementation of SPLASH radiation model
//Davis, T. W., Prentice, I. C., Stocker, B. D., Thomas, R. T., Whitley, R. J., Wang, H., Cramer, W. (2017). Simple process-led algorithms for simulating habitats (SPLASH v.1.0): robust indices of radiation, evapotranspiration and plant-available moisture. Geoscientific Model Development, 10 (2), 689–708
Ground_Radiation::Ground_Radiation(loc l, const orbital &o) : location(l), orb(o) {
    //TODO: Check the orb parameters line up with these equations
        
    // Variable substitutes:
    xee = pow(orb.eccentricity, 2.0);
    xec = pow(orb.eccentricity, 3.0);
    xse = sqrt(1.0 - xee);

    // Mean longitude for vernal equinox:
    xlam = (orb.eccentricity/2.0 + xec/8.0)*(1.0 + xse)*sin(pir*orb.perihelion);
    xlam -= xee/4.0*(0.5 + xse)*sin(pir*2.0*orb.perihelion);
    xlam += xec/8.0*(1.0/3.0 + xse)*sin(pir*3.0*orb.perihelion);
    xlam *= 2.0;
    xlam /= pir;
    
    int current_assumed_day = int(days_in_month[0]/2.0);
    for(int month_of_year=0; month_of_year<number_of_months_in_year; month_of_year++) {

        //We want the day right in the middle of the month so we say off the last 'half month' used and then add in the month and then the current half
        if(month_of_year) 
            current_assumed_day = days_in_month[month_of_year - 1]/2.0 + int(1.5*days_in_month[month_of_year]);

        float dlamm = xlam + (current_assumed_day - 80.0)*(360.0/number_of_days_in_year);

        // Mean anomaly:
        float anm = (dlamm - orb.perihelion);
        float ranm = anm*pir;

        // True anomaly:
        float ranv = ranm;
        ranv += (2.0*orb.eccentricity - xec/4.0)*sin(ranm);
        ranv += 5.0/4.0*xee*sin(2.0*ranm);
        ranv += 13.0/12.0*xec*sin(3.0*ranm);
        float anv = ranv/pir;

        // True longitude: 
        float true_longitude = (anv + orb.perihelion); //rad
        if (true_longitude < 0) {
            true_longitude += 360.0;
        } else if (true_longitude > 360) {
            true_longitude -= 360.0;
        }

        float true_anomoly = true_longitude - orb.perihelion; //rad
        if (true_anomoly <0) 
            true_anomoly += 360.0;
        
        //TODO:Be careful about if we're using degrees or rad here
        //Taken from davis2016simple
        distance_factor[month_of_year] = pow((1.0 + cos(pir*true_anomoly)*orb.eccentricity)/(1.0 - pow(orb.eccentricity, 2.0)), 2.0);
        //pow((1+orb.eccentricity*cos(pir*true_anomoly))/(1-orb.eccentricity*orb.eccentricity), 2);
        float declination_angle = asin(sin(pir*true_longitude)*sin(pir*orb.obliquity))/pir;

        //printf("%i %f %f %f\n", month_of_year, declination_angle, true_longitude,orb.obliquity);
        monthly_ru[month_of_year] = sin(pir*declination_angle)*sin(pir*location.lat);
        monthly_rv[month_of_year] = cos(pir*declination_angle)*cos(pir*location.lat);

        float hour_sunset_angle = 0; //rad
        float ru_by_rv = (monthly_ru[month_of_year]/monthly_rv[month_of_year]); 

        if (ru_by_rv >= 1.0) {
            // Polar day (no sunset)
            hour_sunset_angle = 180.0;
        } else if (ru_by_rv <= -1.0) {
            // Polar night (no sunrise)
            hour_sunset_angle = 0.0;
        } else {
            hour_sunset_angle = -1.0*ru_by_rv;
            hour_sunset_angle = acos(hour_sunset_angle);
            hour_sunset_angle /= pir;
        }
        monthly_hour_sunset_angle[month_of_year] = hour_sunset_angle;

        //This is needed for some photosynthetic stuff that we don't use here
        //monthly extraterrestrial solar radiation (ra_d), J/m^2
        //float ra_d = (86400.0/M_PI)*distance_factor[month_of_year]*I_sc;
        //ra_d *= (monthly_ru[month_of_year]*hour_sunset_angle*pir + monthly_rv[month_of_year]*sin(pir*hour_sunset_angle));
    }
}

float* Ground_Radiation::get_monthly_ground_radiation_arrays(float *T, float *s_f, float *output) {
    float rnl = 0;
    float rw = 0;
    float tau =0;
    float hn = 0;

    for(int month_of_year=0; month_of_year<number_of_months_in_year; month_of_year++) {
        //transmittivity (tau), unitless
        tau = (0.25 + 0.5*s_f[month_of_year])*(1.0 + (2.67e-5)*location.elev);

        //net longwave radiation (rnl), W/m^2
        rnl = (0.2 + (1.0 - 0.2)*s_f[month_of_year])*(107 - T[month_of_year]);

        //variable substitute (rw), W/m^2
        rw = (1.0 - 0.17)*tau*distance_factor[month_of_year]*I_sc;

        //Net radiation cross-over hour angle (hn), degrees
        if ((rnl - rw*monthly_ru[month_of_year])/(rw*monthly_rv[month_of_year]) >= 1.0) {
            // Net radiation negative all day
            hn = 0;
        } else if ((rnl - rw*monthly_ru[month_of_year])/(rw*monthly_rv[month_of_year]) <= -1.0) {
            // Net radiation positive all day
            hn = 180.0;
        } else {
            hn = acos((rnl - rw*monthly_ru[month_of_year])/(rw*monthly_rv[month_of_year]));
            hn /= pir;
        }

        //daytime net radiation (rn_d), J/m^2
        float rn_d = pir*hn*(rw*monthly_ru[month_of_year] - rnl) + rw*monthly_rv[month_of_year]*sin(pir*hn);
        rn_d *= (86400.0/M_PI);

        //nighttime net radiation (rnn_d), J/m^2
        float rnn_d = rw*monthly_ru[month_of_year]*(monthly_hour_sunset_angle[month_of_year] - hn)*pir;
        rnn_d += rw*monthly_rv[month_of_year]*(sin(pir*monthly_hour_sunset_angle[month_of_year]) - sin(pir*hn));
        rnn_d += rnl*(M_PI - 2.0*monthly_hour_sunset_angle[month_of_year]*pir + hn*pir);
        rnn_d *= (86400.0/M_PI);

        //printf("T & S_f %i %f %f %f %f %f %f \n", month_of_year, T[month_of_year], s_f[month_of_year], monthtime_ground_radiation[month_of_year], nighttime_ground_radiation[month_of_year], monthly_rw[month_of_year], monthly_hn[month_of_year]);
        //
        
        //output[month_of_year] = 1e-6 * (rn_d - rnn_d) * days_in_month[month_of_year];
        output[month_of_year] = 1e-6 * (rn_d + rnn_d) * days_in_month[month_of_year];
    }
    return output;
}

float Ground_Radiation::get_single_month_radiation(float temp, float s_f, int month_of_year) {
    //transmittivity (tau), unitless
    float tau = (0.25 + 0.5*s_f)*(1.0 + (2.67e-5)*location.elev);

    //net longwave radiation (rnl), W/m^2
    float rnl = (0.2 + (1.0 - 0.2)*s_f)*(107 - temp);

    //variable substitute (rw), W/m^2
    float rw = (1.0 - 0.17)*tau*distance_factor[month_of_year]*I_sc;
    
    float hn = 0;

    //Net radiation cross-over hour angle (hn), degrees
    if ((rnl - rw*monthly_ru[month_of_year])/(rw*monthly_rv[month_of_year]) >= 1.0) {
        // Net radiation negative all day
        hn = 0;
    } else if ((rnl - rw*monthly_ru[month_of_year])/(rw*monthly_rv[month_of_year]) <= -1.0) {
        // Net radiation positive all day
        hn = 180.0;
    } else {
        hn = acos((rnl - rw*monthly_ru[month_of_year])/(rw*monthly_rv[month_of_year]));
        hn /= pir;
    }

    //daytime net radiation (rn_d), J/m^2
    float rn_d = pir*hn*(rw*monthly_ru[month_of_year] - rnl) + rw*monthly_rv[month_of_year]*sin(pir*hn);
    rn_d *= (86400.0/M_PI);

    //nighttime net radiation (rnn_d), J/m^2
    float rnn_d = rw*monthly_ru[month_of_year]*(monthly_hour_sunset_angle[month_of_year] - hn)*pir;
    rnn_d += rw*monthly_rv[month_of_year]*(sin(pir*monthly_hour_sunset_angle[month_of_year]) - sin(pir*hn));
    rnn_d += rnl*(M_PI - 2.0*monthly_hour_sunset_angle[month_of_year]*pir + hn*pir);
    rnn_d *= (86400.0/M_PI);

    //return 1e-6 * (rn_d - rnn_d) * days_in_month[month_of_year];
    return 1e-6 * (rn_d + rnn_d) * days_in_month[month_of_year];
}

////Gets the radiation for a specific day of the year in MJ
//float Ground_Radiation::get_average_daily_radiation_for_month(int month) {
//    return 1e-6 * (daytime_ground_radiation[month] + nighttime_ground_radiation[month]);
//}
//
////Gets the radiation for a specific day of the year in MJ
//float Ground_Radiation::get_radiation_for_month(int month) {
//    return 1e-6 * (daytime_ground_radiation[month] + nighttime_ground_radiation[month]) * days_in_month[month];
//}

///////////////////////////////////////////////////////////////////////////////////////////////////
//Moisture Corrector///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
Moisture_Corrector::Moisture_Corrector(
        float precip_modern_in, 
        const Vec_f relative_humidity_modern_in, 
        const Vec_f relative_humidity_past_in, 
        const Vec_f temp_modern_in, 
        const Vec_f s_f_modern_in,
        const loc &loc_past_in, 
        const orbital &orb_past_in, 
        const loc &loc_mod_in, 
        const orbital &orb_mod_in) : 
    precip_modern(precip_modern_in), 
    relative_humidity_modern(relative_humidity_modern_in), 
    relative_humidity_past(relative_humidity_past_in), 
    temp_modern(temp_modern_in), 
    s_f_modern(s_f_modern_in), 
    loc_past(loc_past_in), 
    orb_past(orb_past_in), 
    radiation_palaeo(loc_past_in, orb_past_in), modern_e(number_of_months_in_year), modern_e_wrt_T(number_of_months_in_year) {

    float c_a      = standard_atmospheric_pressure*orb_mod_in.co2*1e-6;

    for(int month_of_year =0; month_of_year<number_of_months_in_year; month_of_year++) {
        if(relative_humidity_modern.get(month_of_year) >= 100) {
            modern_e.set(month_of_year, 0);
            modern_e_wrt_T.set(month_of_year, 0);
            continue; 
        }

        //The modern e    
        float cur_mod_temp  = temp_modern.get(month_of_year);
        float eta           = wang_eta(cur_mod_temp);
        float eta_wrt_T     = wang_eta_wrt_T(cur_mod_temp);
        float K             = wang_K(cur_mod_temp);
        float K_wrt_T       = wang_K_wrt_T(cur_mod_temp);
        float Gamma         = wang_compensation_point(cur_mod_temp);
        float Gamma_wrt_T   = wang_compensation_point_wrt_T(cur_mod_temp);
        float D             = saturation_vapour_pressure(cur_mod_temp) * (1-relative_humidity_modern.get(month_of_year)/100.0);
        float D_wrt_T       = saturation_vapour_pressure_wrt_temp(cur_mod_temp) * (1-relative_humidity_modern.get(month_of_year)/100.0);
        float xi            = colin_P_model_factor(eta, K, Gamma);
        float xi_wrt_T      = colin_P_model_factor_wrt_T(xi, eta, K, Gamma, eta_wrt_T, K_wrt_T, Gamma_wrt_T);

        //TODO:This seems really high!
        modern_e.set(month_of_year, 1.6*D/(c_a - Gamma) * (1 + 1/xi/sqrt(D))); 
        modern_e_wrt_T.set(month_of_year,   
                //1.6 / (c_a - Gamma) * ((D_wrt_T * ( c_a - Gamma) + D * Gamma_wrt_T) * ( 1 + 1/xi ) / (c_a - Gamma) - D * xi_wrt_T / pow(xi,2))
                1.6 / (c_a - Gamma) * (D_wrt_T * (1 + 1/2.0/xi/sqrt(D)) - xi_wrt_T * sqrt(D)/(xi*xi) + modern_e.get(month_of_year) * Gamma / 1.6)
                );

        /*printf("%f %f %f %f %f %f %f %f %f %f %f\n",
                //1.6*D/(c_a - Gamma),
                //sqrt(ratio_of_carboxylation_and_transpiration_25 * (K + Gamma)/(1.6 * D)),
                cur_mod_temp    , 
                relative_humidity_modern.get(month_of_year),
                K            , 
                K_wrt_T      , 
                c_a          ,
                Gamma        , 
                Gamma_wrt_T  , 
                D            , 
                D_wrt_T      ,
                modern_e.get(month_of_year), 
                modern_e_wrt_T.get(month_of_year)
              );*/
    }     
}

float Moisture_Corrector::true_mi(float annual_precip_nd, const Vec_f &temperature_nd, const Vec_f &sunf, const orbital &orb, const loc &location) {
    Ground_Radiation rad = Ground_Radiation(location, orb);
    
    //Dimensionalise the temperature
    Vec_f temperature_d = Vec_f(temperature_nd);
    for(int i=0; i<temperature_d.get_size(); i++)
        temperature_d.set(i, Wang_Dimension::dim_temp(temperature_nd.get(i)));

    // Calculate radiation
    float *radiation = new float[12];
    rad.get_monthly_ground_radiation_arrays(temperature_d.get_all(), (sunf*((float)1.0/(float)100.0)).get_all(), radiation); 

    // Calculate the equlibrium evapotranspiration
    float e_q_sum = 0;
    float e_s_wrt_T = 0;

    for(int month_of_year=0; month_of_year<number_of_months_in_year; month_of_year++) {
        if(radiation[month_of_year] >0) { 
            e_s_wrt_T = saturation_vapour_pressure_wrt_temp_nd(temperature_nd.get(month_of_year));
            if(isinf(e_s_wrt_T) or e_s_wrt_T != e_s_wrt_T)
                e_q_sum += radiation[month_of_year];
            else
                e_q_sum += radiation[month_of_year]*e_s_wrt_T/(e_s_wrt_T + psychrometer_constant); 
        }
    }
    
    // Calculate the moisture index
    float output_mi = Wang_Dimension::dim_precip(annual_precip_nd)*latent_heat_of_vaporisation_of_water/e_q_sum;
    
    if(Wang_Dimension::dim_precip(annual_precip_nd)>3000000) {

    //if(location.lat == -35.000000 and location.lon == 139.000000) {
        printf("MI %f\n", output_mi);
        printf("P %f\n", Wang_Dimension::dim_precip(annual_precip_nd));
        printf("E_q %f\n", e_q_sum/latent_heat_of_vaporisation_of_water);
        printf("LL %f %f\n", location.lat, location.lon);

        for(int month_of_year=0; month_of_year<number_of_months_in_year; month_of_year++) {
            e_s_wrt_T = saturation_vapour_pressure_wrt_temp_nd(temperature_nd.get(month_of_year));
            printf("%i %f %f %f %f\n", month_of_year, e_s_wrt_T, radiation[month_of_year], e_s_wrt_T/(e_s_wrt_T + psychrometer_constant), temperature_d.get(month_of_year));
        }

        printf("\n");
    }

    //TODO:LABEL
    //if(output_mi != output_mi)
    //    printf("Nan MI, Corrected \n");
    //else if (!isfinite(output_mi)) {
    //    //printf("Infinate MI\n");
    //    output_mi = 1e10;
    //}

    delete[] radiation;
    return output_mi;
}

//Given a climate find the uncorrected MI if the WUE was swap for the one given at the class initialisation 
float Moisture_Corrector::uncorrect_mi(
        float annual_precip_nd, 
        const Vec_f &temp_past) {
    //Initialise variables
    float e_s_wrt_T = 0;
    float e_q_sum = 0;

    Vec_f temperature_d = Vec_f(temp_past); //Dimensionalise the temperature

    for(int i=0; i<temperature_d.get_size(); i++)
        temperature_d.set(i, Wang_Dimension::dim_temp(temp_past.get(i)));

    float *radiation = new float[12];

    radiation_palaeo.get_monthly_ground_radiation_arrays(temperature_d.get_all(), (s_f_modern*((float)1.0/(float)100.0)).get_all(), radiation);  

    for(int month_of_year =0; month_of_year<number_of_months_in_year; month_of_year++) {
        if(radiation[month_of_year] > 0) {
            e_s_wrt_T = e_s_wrt_temp_alg(
                    temp_past.get(month_of_year), 
                    relative_humidity_past.get(month_of_year), 
                    orb_past.co2, 
                    modern_e.get(month_of_year),
                    modern_e_wrt_T.get(month_of_year)
                    );

            if(isinf(e_s_wrt_T) or e_s_wrt_T != e_s_wrt_T)
                e_q_sum += radiation[month_of_year];
            else
                e_q_sum += radiation[month_of_year]*e_s_wrt_T/(e_s_wrt_T + psychrometer_constant); 
        } 
    }

    //printf("e_q_sum %f %f %f\n", annual_precip_nd, I_sc, e_q_sum);
    //printf("precip %f %f %f \n", annual_precip_nd, e_q_sum, exp(annual_precip_nd)*I_sc*latent_heat_of_vaporisation_of_water/e_q_sum);
    
    float output_mi = Wang_Dimension::dim_precip(annual_precip_nd)*latent_heat_of_vaporisation_of_water/e_q_sum;

    //if(annual_precip_nd < 1)
    //    output_mi = exp(annual_precip_nd-1)*I_sc/e_q_sum;
    
    //if(output_mi != output_mi) { TODO: LABEL
    if(false) {
        printf("Nan MI %f %f \n", annual_precip_nd, e_q_sum);

        for(int month_of_year =0; month_of_year<number_of_months_in_year; month_of_year++) {
            printf("MI info dump: %f %f %f %f %f\n",
                    temp_past.get(month_of_year), 
                    relative_humidity_past.get(month_of_year), 
                    orb_past.co2, 
                    modern_e.get(month_of_year),
                    modern_e_wrt_T.get(month_of_year));

            e_s_wrt_T = e_s_wrt_temp_alg(
                    temp_past.get(month_of_year), 
                    relative_humidity_past.get(month_of_year), 
                    orb_past.co2, 
                    modern_e.get(month_of_year),
                    modern_e_wrt_T.get(month_of_year)
                    );

            printf("e_s stuff: %f %f %i %i %i\n", radiation[month_of_year], e_s_wrt_T, isinf(e_s_wrt_T), isnan(e_s_wrt_T), e_s_wrt_T == e_s_wrt_T );

        }
    }
    else if (!isfinite(output_mi)) {
        //printf("Infinate MI\n");
        output_mi = 1e10;
    }

    delete[] radiation;
    return output_mi;
}

//TODO:Currently unsupported
//float Moisture_Corrector::uncorrect_mi_wrt_T(
//        float mi,
//        float annual_precip, 
//        float temp_past,
//        int month_index 
//        ) {
//    //Dimensionalise the temperature
//    float temperature_d = temp_scaler*temp_past;
//    float radiation = radiation_palaeo.get_single_month_radiation(temperature_d, s_f_modern.get(month_index) * 1.0/100.0, month_index);
//
//    float e_s_wrt_T = e_s_wrt_temp_alg(
//                    temp_past,
//                    relative_humidity_past.get(month_index), 
//                    orb_past.co2, 
//                    modern_e.get(month_index),
//                    modern_e_wrt_T.get(month_index)
//                    );
//    float x = temperature_d;
//    float h = sqrt(std::numeric_limits<float>::epsilon()) * x;
//    volatile float xph = x + h;
//    float dx = xph - x;
//
//    float radiation_wrt_T_nd = celcius_to_kelvin_addition*exp(temp_past) * (radiation_palaeo.get_single_month_radiation(xph,s_f_modern.get(month_index) * 1.0/100.0, month_index)- radiation)/dx;
//
//    float output_mi_wrt_T = -pow(mi,2)/(exp(annual_precip)*I_sc) * (radiation_wrt_T_nd * e_s_wrt_T + radiation * psychrometer_constant/(e_s_wrt_T + psychrometer_constant)) /(e_s_wrt_T + psychrometer_constant);
//    
//    //printf("%f %f %f %f\n", radiation_wrt_T_nd, output_mi_wrt_T, e_s_wrt_T, radiation);
//    if(output_mi_wrt_T != output_mi_wrt_T)
//        printf("Nan MI derivative\n");
//    
//
//    return output_mi_wrt_T;
//}

//Algebraic diff
float e_s_wrt_temp_alg(
        float temp_past_nd, 
        float relative_humidity_past,
        float c_a_past,
        float mod_e,
        float mod_e_wrt_T
        ) {
    /*printf("Inputs %f %f %f %f %f\n", 
            temp_past_nd,
            relative_humidity_past,
            c_a_past,
            mod_e,
            mod_e_wrt_T     );*/

    //If relative humidity is really high then we know d e_s/dt is inf
    if(relative_humidity_past >= 100)
        return numeric_limits<float>::infinity();
    
    //convert c_a's into kPa
    float c_a =    standard_atmospheric_pressure*c_a_past*1e-6;  

    float eta          = wang_eta_nd(temp_past_nd);
    float eta_wrt_T    = wang_eta_wrt_T_nd(temp_past_nd);
    float Gamma        = wang_compensation_point_nd(temp_past_nd); 
    float Gamma_wrt_T  = wang_compensation_point_wrt_T_nd(temp_past_nd);
    float K            = wang_K_nd(temp_past_nd);
    float K_wrt_T      = wang_K_wrt_T_nd(temp_past_nd);
    float xi           = colin_P_model_factor(eta, K, Gamma);
    float xi_wrt_T     = colin_P_model_factor_wrt_T(xi, eta, K, Gamma, eta_wrt_T, K_wrt_T, Gamma_wrt_T);
    float delta        = 1/pow(xi, 2) + 4 * mod_e / 1.6 * (c_a - Gamma);
    float delta_wrt_T  = -2 * xi_wrt_T/pow(xi,3) + 4/1.6 * (mod_e_wrt_T * (c_a - Gamma) + mod_e * Gamma_wrt_T);
    float D_wrt_T      = 0.25 * ( delta_wrt_T - xi_wrt_T * 2/pow(xi,3) + delta_wrt_T / (xi * sqrt(delta)) - 2 * sqrt(delta) / pow(xi, 2) * xi_wrt_T );
    //TODO: Cur change: float D_wrt_T      = 0.25 * (xi_wrt_T * (-2)/pow(xi,3) + (0.5*xi*delta_wrt_T*pow(delta, -0.5) - sqrt(delta)*xi_wrt_T)/pow(xi,2) + delta_wrt_T);

    /*printf("%f %f %f %f %f %f %f %f %f %f %f %f\n", 
            eta          ,
            eta_wrt_T    ,
            Gamma        ,
            Gamma_wrt_T  ,
            K            ,
            K_wrt_T      ,
            d            ,
            d_wrt_T      ,
            delta        ,
            delta_wrt_T  ,
            D_wrt_T      ,
            D_wrt_T / (1-relative_humidity_past/100.0));*/

    return D_wrt_T / (1-relative_humidity_past/100.0);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//Dimensionality functions/////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////
float Wang_Dimension::dim_precip(float non_P) {
    if (non_P < p_dim_cutoff_factor/(I_sc/ latent_heat_of_vaporisation_of_water))
        return exp(non_P -1)* I_sc/ latent_heat_of_vaporisation_of_water;
    else
        return non_P * I_sc/ latent_heat_of_vaporisation_of_water;
}

float Wang_Dimension::dim_precip_derivative(float non_P) {
    if (non_P < p_dim_cutoff_factor/(I_sc/ latent_heat_of_vaporisation_of_water))
        return exp(non_P - 1)* I_sc/ latent_heat_of_vaporisation_of_water;
    else 
        return I_sc/ latent_heat_of_vaporisation_of_water;
}

float Wang_Dimension::nondim_precip(float P) {
    if(P < p_dim_cutoff_factor) {
        if(P != 0) 
            return log(P*latent_heat_of_vaporisation_of_water/I_sc) + 1;
        else
            return -1e10;
    } else 
        return P*latent_heat_of_vaporisation_of_water/I_sc;
}

float Wang_Dimension::nondim_precip_derivative(float P) {
    if(P < p_dim_cutoff_factor) 
        if(P != 0) 
            return 1/P; 
        else
            return 1e10;
    else
        return latent_heat_of_vaporisation_of_water/I_sc; 
}

float Wang_Dimension::nondim_temp(float temp) {
    return temp/temp_scaler;
}

float Wang_Dimension::dim_temp(float non_temp) {
    return non_temp * temp_scaler;
}

//Access functions
float true_mi_simple(float precip, float* temps, float* sunf, float* location, int orb_choice) {
    orbital orb_to_use = LGM_orbital;
    switch(orb_choice) {
        case 0:
            orb_to_use = FIF_orbital;
        case 1:
            orb_to_use = MH_orbital;
        case 2:
        default:
            orb_to_use = LGM_orbital;
    }
    
    Vec_f temps_to_use = Vec_f(12);
    for(int i=0; i<12; i++) 
        temps_to_use.set(i, Wang_Dimension::nondim_temp(temps[i]));

    return Moisture_Corrector::true_mi(Wang_Dimension::nondim_precip(precip), temps_to_use, Vec_f(12, sunf), orb_to_use, loc(location[0], location[1], location[2]));
}
