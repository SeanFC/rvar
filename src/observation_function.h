#ifndef observation_function_h
#define observation_function_h

#include <thread>
#include <future>

#include "lin_maths.h"

using namespace std;

const int thread_cap = 14;

//A generic class that takes a state and gives an observation
//State to Observation
class sto {
    protected:
        float eps = std::numeric_limits<float>::epsilon();  // The numeric limits of the system this is being run on 
        int number_of_threads;                              // How many threads we're going to give this
        bool TIMING_VERBOSE = false;                        // uUse to test timings

    public:
        sto(int no_o_threads = 0) : number_of_threads(no_o_threads) {
            //Pick a number of threads that will use the maximum system resources that I can get away with and leave some processing headspace for management
            number_of_threads = number_of_threads == 0 ? 1 : number_of_threads;
            number_of_threads = number_of_threads > 0 ? number_of_threads : std::thread::hardware_concurrency() - 1;
            number_of_threads = number_of_threads > thread_cap ? thread_cap : number_of_threads;
        }

        virtual ~sto() {};

        virtual Vec_f state_to_observation(const Vec_f&)=0;
        Mat_f jacobian(Vec_f);
        Mat_f partial_jacobian(Vec_f point, Mat_f *output, int start_index, int finish_index, Vec_f run_at_point);
        virtual Vec_f derivative(const Vec_f &point, const Vec_f &reference_run, int index);

        virtual Mat_f state_to_observation_variance(Vec_f state, Mat_f covariance);

        //Find the jacobian in block matrix form
        //TODO:Note that this is poor implementation, really the jacobian function should just output a general matrix, then child classes could overload jacobian to place in any matrix needed however this doesn't work well with the current lin_maths setup
        virtual Block_f block_jacobian(const Vec_f&);
};

#endif
