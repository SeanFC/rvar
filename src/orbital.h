#ifndef orbital_h 
#define orbital_h 

#include <stdio.h>

//Set known characteristics of a climate 
struct orbital {
    //Orbital paramters
    double eccentricity;
    double obliquity; //Degrees
    double perihelion; //Degrees

    //Atmospheric CO_2 concentration in ppm
    double co2;     

    //Contstructor
    orbital(double e, double o, double p, double c) : eccentricity(e), obliquity(o), perihelion(p), co2(c) {}
    
    void print() const { printf("%f %f %f\n", eccentricity, perihelion, obliquity); };
};

//Set values for orbital charactoers
const orbital FIF_orbital(0.016724, 23.446, 102.04+180, 360); //1950CE, fifties
const orbital MH_orbital(0.018682, 24.105, 0.87+180, 264.4); //6ka
const orbital LGM_orbital(0.018994, 22.949, 114.42+180, 180); //21ka

//To be able to easily cycle through the different orbitals
enum Orbitals { FIF, LGM, MH };

//Get a specific orbital based on a code
inline orbital get_orbital(int orb_ident) {
    switch(orb_ident) {
        case FIF:
            return FIF_orbital;
        case LGM:
            return LGM_orbital;
        case MH:
            return MH_orbital;
    }

    return FIF_orbital;
}

#endif
