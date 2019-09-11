#include "location.h"

//Give the greater circle distance between two locations in km 
double loc::distance_between(const loc& other) const {
    return get_earth_radius_km()*angle_distance_between(other);
}

//Give the angle between two locations 
//The wikipedia article is a good source of information on this
double loc::angle_distance_between(const loc& other) const {
    //Precalculate latitude and longitude in radians
    double lat_rad_this = lat/360.0 * 2 * M_PI;
    double lon_rad_this = lon/360.0 * 2 * M_PI;
    double lat_rad_othe = other.lat/360.0 * 2 * M_PI;
    double lon_rad_othe = other.lon/360.0 * 2 * M_PI;
    
    //Find the difference between the two 
    double diff_dist = sin(lat_rad_this) * sin(lat_rad_othe)  + cos(lat_rad_this) * cos(lat_rad_othe) * cos(lon_rad_this - lon_rad_othe);
    return acos(diff_dist);
}
