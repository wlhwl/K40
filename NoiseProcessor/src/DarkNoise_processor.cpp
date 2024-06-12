//
// Created by cmo on 2023/12/16.
//
#include <string>
#include <iostream>
#include <vector>

#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "../include/DarkNoise_processor.hpp"
#include "../include/Event_processor.hpp"


void process_pmt_darknoise(double event_t_min, double T){
		NUM_PMTs = config["DarkNoise"]["num_pmt"].as<int>();
        NoiseRate = config["DarkNoise"]["pmtrate"].as<double>() / 1e9; // ns
		int total_noise = 0;
        for(const auto& domId: DomId){
            for(int pmtId=0; pmtId<NUM_PMTs; pmtId++){
                // Sample number of noise hits
                int num_noise = rnd_gen->Poisson(T * NoiseRate);
				total_noise += num_noise;
                for(int iNoise=0; iNoise<num_noise; iNoise++){
                    // Sample noise time
                    auto t_noise = rnd_gen->Uniform(T) + event_t_min;

                    // Fill in noise data
                    new_t.emplace_back(t_noise);
                    new_pmt_id.emplace_back(pmtId);
                    new_DOMid.emplace_back(domId);
                    new_x.emplace_back(array.at(domId).x());
                    new_y.emplace_back(array.at(domId).y());
                    new_z.emplace_back(array.at(domId).z());
                    new_hit_type.emplace_back(0);
                }
            }
        }
		std::cout<<"Event duration: " << T <<" ns. Number of dark noise hits: "<<total_noise<<std::endl;

}

