#include "observation_function.h"

Mat_f sto::jacobian(Vec_f point) {
    Vec_f run_at_point = state_to_observation(point);
    Mat_f output(run_at_point.get_size(), point.get_size());

    //Pre save some timings
    struct timespec start_pos, finish_pos;

    if(TIMING_VERBOSE) clock_gettime(CLOCK_REALTIME, &start_pos);
    if(number_of_threads>1) {
        //openblas_set_num_threads(1);
        vector<thread> threads; 
        int amount_of_calcs = point.get_size();
        int calculation_amount = (int)(amount_of_calcs/(1.0 * number_of_threads));
        int end_size = amount_of_calcs - (calculation_amount*number_of_threads) + calculation_amount;

        //Copy a load of points so that they can be sent to the individual threads and messed around with independantly
        //TODO:If these points are const then there isn't really a need to make loads of them
        vector<Vec_f> points = vector<Vec_f>();
        for(int i=0; i<number_of_threads; i++) points.push_back(Vec_f(point));

        //Copy a load of jacobian matrix that can be used as the outputs of the threads 
        vector<Mat_f> jacs = vector<Mat_f>();
        for(int i=0; i<number_of_threads-1; i++) jacs.push_back(Mat_f(run_at_point.get_size(), calculation_amount));
        jacs.push_back(Mat_f(run_at_point.get_size(), end_size));

        for(int i=0; i<number_of_threads; i++) { 
            int start_point = i*calculation_amount;
            int end_point=0;
            if(i == number_of_threads - 1) 
                end_point = amount_of_calcs;
            else
                end_point = (i+1)*calculation_amount;

            threads.push_back(thread(
                        &sto::partial_jacobian,
                        this, 
                        points[i], 
                        &jacs[i],
                        start_point,
                        end_point,
                        run_at_point
                        ));
        }

        for(auto& t: threads) t.join(); //Bring all the threads together to let them finish

        for(int i=0; i<number_of_threads; i++)
            output.set(0, i*calculation_amount, jacs[i]);

        threads.clear();

        //openblas_set_num_threads(number_of_threads);
    } else {
        output = partial_jacobian(point, &output, 0, point.get_size(), run_at_point);
    }
    
    if(TIMING_VERBOSE) {
        clock_gettime(CLOCK_REALTIME, &finish_pos);

        double time_elapsed = (finish_pos.tv_sec - start_pos.tv_sec)/ 1000000000.0;

        printf("POS Timing: %f\n", time_elapsed);
    }

    return output; 
}

Mat_f sto::partial_jacobian(Vec_f point, Mat_f *output, int start_index, int finish_index, Vec_f run_at_point) {
    for(int i=0; i<finish_index-start_index; i++)
        output->set(0, i, derivative(point, run_at_point, i+start_index)); 

    return *output;
}

Vec_f sto::derivative(const Vec_f &point, const Vec_f &reference_run, int index) {
    //TODO:Double here?
    float x = point.get(index);
    float h = sqrt(eps) * x;
    volatile float xph = x + h;
    float dx = xph - x;

    Vec_f point_copy = point;

    point_copy.set(index, xph);
    Vec_f output = (state_to_observation(point_copy) - reference_run)/dx;
    point_copy.set(index, x); 

    return output; 
}

Mat_f sto::state_to_observation_variance(Vec_f state, Mat_f covariance) {
    Mat_f func_jacobian = jacobian(state);

    return func_jacobian*covariance*Transpose(func_jacobian);
}

Block_f sto::block_jacobian(const Vec_f &point) {
    printf("sto::block_jacobian - Implement me\n");

    return Block_f(1,1,1,1);
}
