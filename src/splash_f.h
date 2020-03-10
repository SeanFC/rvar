#ifndef splash_f_h 
#define splash_f_h 

#include <stdint.h>
//Note that these types depend highly on how we compile the fortran library, especially in relation to default lengths of REAL etc
//Useful links:
//http://www.math.utah.edu/software/c-with-fortran.html#data-type-diffs
//http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
//http://www.cprogramming.com/tutorial/shared-libraries-linux-gcc.html
//http://www.neurophys.wisc.edu/comp/docs/not017/

//We only need the extern C if we're using cpp, otherwise we can define the functions as normal
#ifdef __cplusplus
extern "C" { 
#endif

typedef struct {
    double sm;     // soil moisture (=water content, mm)
    double ro;     // daily runoff (mm)
} waterbaltype;

// function return value of evap() as derived type
typedef struct {
    double ra;         // daily top-of-atmosphere solar irradiation (J/m2)
    double rn;         // daily net radiation (J/m2)
    double ppfd;       // daily photosynthetic photon flux density (mol/m^2)
    double eet;        // daily equilibrium evapotranspiration (mm)
    double pet;        // daily potential evapotranspiration (mm)
    double aet;        // daily actual evapotranspiration (mm)
    double cn;         // daily condensation (mm)
} outtype_evap;

// function return value of berger_tls() as derived type
typedef struct {
    double nu;
    double lambda;
} outtype_berger;

//double __splash_MOD_outmeet[12];

//double __splash_MOD_outdaet[366];
//double __splash_MOD_outdeet[366];
//

double __splash_MOD_dgsin(double*); 
double __splash_MOD_dgcos(double*); 
//double __splash_MOD_berger_tls
//double __splash_MOD_degrees
//double __splash_MOD_dgcos
//double __splash_MOD_dgsin
//double __splash_MOD_elv2pres
//double __splash_MOD_evap
//double __splash_MOD_get_density_h2o
//double __splash_MOD_get_enthalpy_vap
//double __splash_MOD_get_julian_day
//double __splash_MOD_getout_daily
//double __splash_MOD_getout_monthly
//double __splash_MOD_get_psychro
//double __splash_MOD_get_sat_slope
//double __splash_MOD_initdaily
//double __splash_MOD_initmonthly
//double __splash_MOD_outdaet
//double __splash_MOD_outdcn
//double __splash_MOD_outdeet
//double __splash_MOD_outdpet
//double __splash_MOD_outdppfd
//double __splash_MOD_outdra
//double __splash_MOD_outdrn
//double __splash_MOD_outdro
//double __splash_MOD_outdsm
//double __splash_MOD_out_evap
extern double __splash_MOD_outmaet[12];
//extern double __splash_MOD_outmcpa;
//double __splash_MOD_outmcwd
extern double __splash_MOD_outmeet[12];
//double __splash_MOD_outmpet
//double __splash_MOD_outmppfd
//double __splash_MOD_radians
//double __splash_MOD_run_one_day
//double __splash_MOD_run_one_year
void __splash_MOD_spin_up_sm(int64_t*, double*, double*, double[], double[], double[], int64_t*, double*, double*, double*);
//double __splash_MOD_use_daily_input
extern bool __splash_MOD_verbose;
//double __splash_MOD_waterbal
//double __splash_MOD_write_to_file

struct non_dim_splash_class{
    int no_of_time_segments;
    double elevation, latitude, eccentricity, obliquity, long_of_perihelion;
    double precip_scale=0, temp_shift=0, sunf_shift=0, temp_scale=0, latent_slope_exp_coefficient=0;
    double latent_slope_coefficient=0, latent_slope_temp_shift=0, water_boiling_const=0;
    double sunf_scale[365], incident_lw_shift_const[365], denc_angle[365], heli_long[365];
}; //__non_dim_splash_MOD_splash_class;

struct mi_jac_hess {
    double mi;
    double jacobian[25]; //Should be a 25x1
    double hessian[25*25]; //Should be a 25x25
};

void __non_dim_splash_MOD_splash_class(non_dim_splash_class*, 
        double*, 
        double*, 
        int64_t*, 
        double*, 
        double*, 
        double*);
double __non_dim_splash_MOD_moisture_index_dim_procedural(non_dim_splash_class*, 
        double*, 
        double[], 
        double[]);
mi_jac_hess __non_dim_splash_MOD_moisture_index_dim_jh_procedural(non_dim_splash_class*, 
        double*, 
        double[], 
        double[]);

#ifdef __cplusplus
}
#endif

#endif
