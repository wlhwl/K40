#include <iostream>


class K40_generation {
    public:
        K40_generation(std::string config_file);
        ~K40_generation();
        void generate_K40(int domid);
        int DomId;

    private:
        double K40poststepdoit(double t0, double t);
        void K40alongstepdoit(int cl, double t);
        double T_e;
        double CL1, CL2, CL3, CL4, CL5;
};