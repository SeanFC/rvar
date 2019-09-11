#ifndef location_h
#define location_h

#include <cmath> //For abs and fmod
#include <ostream> //For ostream

#define EARTH_RADIUS_KM 6371.0

//A specific location on the surface of a sphere using latitude and longitude
struct loc {
    public: 
        //The variables in 3D space
        double lat, lon, elev;

        //Constructors
        loc() { lat = 0; lon = 0; elev = 0;} 
        loc(double la, double lo, double el) { lat = la; lon = lo; elev = el; }

        void print() const { printf("%f, %f, %f\n", lat, lon, elev); }
       
        //Equating operators 
        bool operator==(const loc &o) const { return (lat == o.lat && lon == o.lon && elev == o.elev); }
        bool equal_flat(const loc &o) const { return (lat == o.lat && lon == o.lon); }

        //Linear operators
        loc operator-(const loc &rhs) const { return loc(lat - rhs.lat, lon - rhs.lon, elev - rhs.elev); } 
        loc operator+(const loc &rhs) const { return loc(lat + rhs.lat, lon + rhs.lon, elev + rhs.elev); } 

        //Comparison operators, note that these need to work on a curved surface
        double distance_between(const loc& other) const;
        double angle_distance_between(const loc& other) const;

        //What we're using for the Earth's radius 
        static double get_earth_radius_km() { return EARTH_RADIUS_KM; };
};

//Used to print the variables
inline std::ostream& operator<<(std::ostream& os, const loc& l)
{
    os << l.lat << ' ' << l.lon << ' ' << l.elev;
    return os;
}

#endif
