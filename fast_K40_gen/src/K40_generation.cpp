/*
Treat each CL=n event as Poisson process. Select CL=n corresponding to the minimal time every step.
*/
#include "yaml-cpp/yaml.h"
#include "math.h"
#include "iostream"
#include "random"
#include "../include/K40_generation.hpp"

K40_generation::K40_generation(std::string config_file) {
    YAML::Node config = YAML::LoadFile(config_file);
    T_e = config["event_duration"].as<double>();
    CL1 = config["CL_rate"]["CL1"].as<double>();
    CL2 = config["CL_rate"]["CL2"].as<double>();
    CL3 = config["CL_rate"]["CL3"].as<double>();
    CL4 = config["CL_rate"]["CL4"].as<double>();
    CL5 = config["CL_rate"]["CL5"].as<double>();
}

K40_generation::~K40_generation() {}

void K40_generation::generate_K40(int domid) {
    //scale time(ns or s) by 20ns
    double t = T_e/20;
    double t0 = 0;
    DomId = domid;
    while (t0<t){
        t0 = K40poststepdoit(t0,t);
    }
}

double K40_generation::K40poststepdoit(double t0, double t) {
    std::random_device rd;
    std::mt19937 gen(rd());

    std::exponential_distribution<double> exp_CL1(CL1/(5*pow(10,7)));
    std::exponential_distribution<double> exp_CL2(CL2/(5*pow(10,7)));
    std::exponential_distribution<double> exp_CL3(CL3/(5*pow(10,7)));
    std::exponential_distribution<double> exp_CL4(CL4/(5*pow(10,7)));
    std::exponential_distribution<double> exp_CL5(CL5/(5*pow(10,7)));

    double t1 = exp_CL1(gen);
    double t2 = exp_CL2(gen);
    double t3 = exp_CL3(gen);
    double t4 = exp_CL4(gen);
    double t5 = exp_CL5(gen);

    double t_min = std::min({t1, t2, t3, t4, t5});
    if(t_min+t0>t){return t;}
    if (t_min == t1){
        K40alongstepdoit(1, t_min+t0);
    } else if (t_min == t2){
        K40alongstepdoit(2, t_min+t0);
    } else if (t_min == t3){
        K40alongstepdoit(3, t_min+t0);
    } else if (t_min == t4){
        K40alongstepdoit(4, t_min+t0);
    } else if (t_min == t5){
        K40alongstepdoit(5, t_min+t0);
    }

    return t_min+t0;
}

void K40_generation::K40alongstepdoit(int cl, double t) {

    double start_time = floor(t)*20;
    std::cout<<DomId<<","<<cl<<","<<start_time<<std::endl;
    //generate n hits in [start_time, start_time+20ns]

}