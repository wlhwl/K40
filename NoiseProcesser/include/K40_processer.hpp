#include <iostream>
#include <random>

class K40_processer {
    public:
        K40_processer();
        ~K40_processer();
        void process_K40(double ,double );

    private:
        double K40poststepdoit(double t0, double t);
        void K40alongstepdoit(int cl, double t);
        int domid;
        int total_noise;
        double T;
        double event_t_min;
        double CL1, CL2, CL3, CL4, CL5;
        std::exponential_distribution<double> exp_CL1;
        std::exponential_distribution<double> exp_CL2;
        std::exponential_distribution<double> exp_CL3;
        std::exponential_distribution<double> exp_CL4;
        std::exponential_distribution<double> exp_CL5;
        std::random_device rd;
};