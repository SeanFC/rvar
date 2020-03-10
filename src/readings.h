/////////////////////////////////////////
// Readings /////////////////////////////
/////////////////////////////////////////
// Classes for readings of the observation, modern or state space, together with additional useful properties such as location and associated ID

//am is for alpha moisture

#ifndef readings_h 
#define readings_h 

#include <vector>
#include <stdlib.h>
#include <cmath>
#include <iostream>

#include "lin_maths.h"
#include "location.h"

using namespace std;


// Standard readings ////////////////////

// An observation space vector at a single site 
class am_space_reading : public Vec_f {
    public:
        static const int default_size = 7;

        am_space_reading() : Vec_f(default_size) {};
        am_space_reading(Vec_f other) : Vec_f(other) {};

        double get_alpha()      const { return get(0); }
        double get_moisture()   const { return get(1); }
        double get_MAP()        const { return get(2); }
        double get_MAT()        const { return get(3); }
        double get_MTCO()       const { return get(4); }
        double get_MTWA()       const { return get(5); }
        double get_GDD5()       const { return get(6); }

        void set_alpha(double given) { set(0, given); }
        void set_moisture(double given) { set(1, given); }
        void set_MAP(double given) { set(2, given); }
        void set_MAT(double given) { set(3, given); }
        void set_MTCO(double given) { set(4, given); }
        void set_MTWA(double given) { set(5, given); }
        void set_GDD5(double given) { set(6, given); }

        am_space_reading& operator+=(const am_space_reading& reading) { 
            this->Vec_f::operator+=(reading);
            return *this;
        }
};

inline am_space_reading operator+(am_space_reading lhs, const am_space_reading& rhs) {
    lhs += rhs;
    return lhs;
}

// A regular state space vector, in a single grid cell
class Seasonal_space_reading : public Vec_f {
    public:
        static const int default_size = 13;

        Seasonal_space_reading() : Vec_f(default_size) {};
        Seasonal_space_reading(Vec_f v) : Vec_f(v) {};
        Seasonal_space_reading(double P, const Vec_f &T) : Seasonal_space_reading() {
            set_MAP(P); set_temp(T);
        };

        //The rule of three
        Seasonal_space_reading(const Seasonal_space_reading &to_clone) : Vec_f(to_clone) {};
        virtual Seasonal_space_reading& operator=(const Seasonal_space_reading& rhs) { 
            Vec_f::operator=(rhs); 
            return *this; 
        };
        virtual ~Seasonal_space_reading() {};

        double get_MAP()        const { return get(0); }
        Vec_f get_temp()  const { return get_sub_vec(1,1+12); }

        void set_MAP(double given)          { set(0, given); }
        void set_temp(const Vec_f &given)    { set(1, given); }

        Seasonal_space_reading& operator+=(const Vec_f& reading) { 
            this->Vec_f::operator+=(reading);
            return *this;
        }
};

inline Seasonal_space_reading operator+(Seasonal_space_reading lhs, const Seasonal_space_reading& rhs) {
    lhs += rhs;
    return lhs;
}

// A modern state space vector in a single grid cell
class Seasonal_Modern_reading : public Vec_f {
    public:
        static const int default_size = 37;

        Seasonal_Modern_reading() : Vec_f(default_size) {};
        Seasonal_Modern_reading(Vec_f v) : Vec_f(v) {};

        double get_MAP()        const { return get(0); }
        Vec_f get_temp()  const { return get_sub_vec(1,1+12); }
        Vec_f get_RH()    const { return get_sub_vec(13,13+12); }
        Vec_f get_sunf()  const { return get_sub_vec(25,25+12); }

        void set_MAP(double given)          { set(0, given); }
        void set_temp(Vec_f given)    { set(1, given); }
        void set_RH(Vec_f given)      { set(13, given); }
        void set_sunf(Vec_f given)    { set(25, given); }

        Seasonal_Modern_reading& operator+=(const Seasonal_Modern_reading& reading) { 
            this->Vec_f::operator+=(reading);
            return *this;
        }

        Seasonal_Modern_reading& section_div(const Vec_f& divider) { 
            this->set_MAP  (this->get_MAP()  / divider.get(0));
            this->set_temp (this->get_temp() / divider.get(1));
            this->set_RH   (this->get_RH()   / divider.get(2));
            this->set_sunf (this->get_sunf() / divider.get(3));

            return *this;
        }
};

inline Seasonal_Modern_reading operator+(Seasonal_Modern_reading lhs, const Seasonal_Modern_reading& rhs) {
    lhs += rhs;
    return lhs;
}

// Readings with standard deviation ////////////////////

// An observation space vector with standard deviations attached
class am_sd_reading : public am_space_reading {
    protected:
        am_space_reading sd;
    
    public:
        am_sd_reading() : am_space_reading() { sd = am_space_reading(); }
        am_sd_reading(am_space_reading vals, am_space_reading sds) : am_space_reading(vals) { sd = sds; }
        am_space_reading& get_sd() { return sd; }
        void set_sd(am_space_reading to_set) { sd = am_space_reading(to_set); }

        am_sd_reading& operator+=(const am_space_reading& reading) { 
            this->am_space_reading::operator+=(reading);
            return *this;
        }
};

inline am_sd_reading operator+(am_sd_reading lhs, const am_space_reading& rhs) {
    lhs += rhs;
    return lhs;
}

// A state space vector with standard deviations attached
class Seasonal_sd_reading : public Seasonal_space_reading {
    protected:
        Seasonal_space_reading sd;
    
    public:
        Seasonal_sd_reading() : Seasonal_space_reading() { sd = Seasonal_space_reading(); } 
        Seasonal_sd_reading(const Seasonal_space_reading &s, const Seasonal_space_reading &sds) : Seasonal_space_reading(s) { sd = sds; }

        Seasonal_space_reading& get_sd() { return sd; }
        void set_sd(Seasonal_space_reading &to_set) { sd = Seasonal_space_reading(to_set); }

        Seasonal_sd_reading& operator+=(const Seasonal_space_reading& reading) { 
            this->Seasonal_space_reading::operator+=(reading);
            return *this;
        }
};

inline Seasonal_sd_reading operator+(Seasonal_sd_reading lhs, const Seasonal_sd_reading& rhs) {
    lhs += rhs;
    return lhs;
}

// A modern state space vector, with standard deviations attached
class Modern_sd_reading : public Seasonal_Modern_reading {
    protected:
        Seasonal_Modern_reading sd;
    
    public:
        Modern_sd_reading() : Seasonal_Modern_reading() { sd = Seasonal_Modern_reading(); }
        Seasonal_Modern_reading& get_sd() { return sd; }
        void set_sd(Seasonal_Modern_reading &to_set) { sd = Seasonal_Modern_reading(to_set); }

        //TODO:Change this back!
        Seasonal_sd_reading* get_seasonal_reading_part() { 
            return new Seasonal_sd_reading(
                    Seasonal_space_reading(get_MAP(), get_temp()),//, get_RH()), 
                    Seasonal_space_reading(sd.get_MAP(), sd.get_temp())//, sd.get_RH())
                    );
        }

        Modern_sd_reading& operator+=(const Seasonal_Modern_reading& reading) { 
            this->Seasonal_Modern_reading::operator+=(reading);
            return *this;
        }
};

inline Modern_sd_reading operator+(Modern_sd_reading lhs, const Modern_sd_reading& rhs) {
    lhs += rhs;
    return lhs;
}

// Readings with ID and location ////////////////////
//
// Hold the ID and location of a certain point
class id_loc {
    protected:
        unsigned int main_id;
        loc location;
    
        id_loc(unsigned int i, loc l) : main_id(i), location(l) {} 

    public:
        unsigned int get_main_id() { return main_id; }
        loc get_location() { return location; }
        loc get_loc() { return get_location(); }
};

// An observation space vector, with an ID and location attached
class am_loc_reading : public am_space_reading, public id_loc {
    public:
        am_loc_reading(unsigned int id, loc l) : am_space_reading(), id_loc(id, l) {};
};

// A state space vector, with an ID and location attached
class Seasonal_loc_reading : public Seasonal_space_reading, public id_loc {
    public:
        Seasonal_loc_reading(unsigned int id, loc l) : Seasonal_space_reading(), id_loc(id, l) {};
};

// A modern state space vector with ID and location attached
class Modern_loc_reading : public Seasonal_Modern_reading, public id_loc {
    public:
        Modern_loc_reading(unsigned int id, loc l) : Seasonal_Modern_reading(), id_loc(id, l) {};
};

// Readings with standard devation, ID and location ////////////////////////////

// An observation space vector with sd, ID and location attached
class am_sd_loc_reading : public am_sd_reading, public id_loc {
    bool CO2_corrected;

    public:
        am_sd_loc_reading(unsigned int id, loc l, bool co2) : am_sd_reading(), id_loc(id, l), CO2_corrected(co2) {};

        bool get_CO2_corrected() { return CO2_corrected; }
        void set_CO2_corrected(bool co2) { CO2_corrected = co2; }
};

// A state space vector with sd, ID and location attached
class Seasonal_sd_loc_reading : public Seasonal_sd_reading, public id_loc {
    public:
        Seasonal_sd_loc_reading(unsigned int id, loc l) : Seasonal_sd_reading(), id_loc(id, l) {};
};

// A modern state space vector with sd, ID and location attached
class Modern_sd_loc_reading : public Modern_sd_reading, public id_loc {
    public:
        Modern_sd_loc_reading(unsigned int id, loc l) : Modern_sd_reading(), id_loc(id, l) {};
};

#endif
