/*
Treat each CL=n event as Poisson process. Select CL=n corresponding to the minimal time every step.
*/
#include "yaml-cpp/yaml.h"
#include "math.h"
#include "iostream"
#include "random"
#include "algorithm"
#include "../include/K40_processer.hpp"
#include "../include/Event_processer.hpp"

K40_processer::K40_processer() {
    CL1 = config["K40"]["CL_rate"]["CL1"].as<double>();
    CL2 = config["K40"]["CL_rate"]["CL2"].as<double>();
    CL3 = config["K40"]["CL_rate"]["CL3"].as<double>();
    CL4 = config["K40"]["CL_rate"]["CL4"].as<double>();
    CL5 = config["K40"]["CL_rate"]["CL5"].as<double>();
    exp_CL1 = std::exponential_distribution<double>(CL1/(5*pow(10,7)));
    exp_CL2 = std::exponential_distribution<double>(CL2/(5*pow(10,7)));
    exp_CL3 = std::exponential_distribution<double>(CL3/(5*pow(10,7)));
    exp_CL4 = std::exponential_distribution<double>(CL4/(5*pow(10,7)));
    exp_CL5 = std::exponential_distribution<double>(CL5/(5*pow(10,7)));

}

K40_processer::~K40_processer() {}

void K40_processer::process_K40(double t_min,double time_duration) {
    total_noise=0;
    //scale time(ns or s) by 20ns
    T = time_duration;
    event_t_min = t_min;
    double t = T/20;
    for(const auto& id: DomId){
        domid = id;
        double t0 = 0;
        while (t0<t){
            t0 = K40poststepdoit(t0,t);
        }
    }
    std::cout<<"Event duration: " << T <<" ns. Number of K40 noise hits: "<<total_noise<<std::endl;
}

double K40_processer::K40poststepdoit(double t0, double t) {
    
    std::mt19937 gen(rd());

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

void K40_processer::K40alongstepdoit(int cl, double t) {

    double start_time = floor(t)*20, t_=0;
    total_noise += cl;
    std::vector<double> ghosttime;
    //generate n hits in [start_time, start_time+20ns]
    std::mt19937 gen(rd());
    std::exponential_distribution<double> exp(cl*10 / 20);
    while(t_<20){
        t_ += exp(gen);
        ghosttime.emplace_back(t);
    }
    std::vector<double> hittime(cl);
    if(ghosttime.size()<cl){
        std::cout<<"ghosttime.size(): "<<ghosttime.size()<<"< cl: "<<cl<<std::endl;
        return;
    }
    std::sample(ghosttime.begin(), ghosttime.end(), hittime.begin(), cl, gen);
    for (int i=0; i<cl; i++){
        new_t.emplace_back(event_t_min+start_time+hittime.at(i));
        new_pmt_id.emplace_back(-1);
        new_DOMid.emplace_back(domid);
        new_x.emplace_back(array.at(domid).x());
        new_y.emplace_back(array.at(domid).y());
        new_z.emplace_back(array.at(domid).z());
        new_hit_type.emplace_back(-1);
    }
}