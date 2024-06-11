#include "../include/Event_processer.hpp"
#include "../include/K40_processer.hpp"
#include "../include/DarkNoise_processer.hpp"
#include "yaml-cpp/yaml.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <numeric>
#include <algorithm>

Event_processer::Event_processer(std::string config_file) {
    
    config = YAML::LoadFile(config_file);
    max_time = config["max_time"].as<double>();
    std::ifstream csvGeometry(config["fileGeometry"].as<std::string>());
    float x, y, z;
    std::string line;
    while (std::getline(csvGeometry, line)) {
        std::sscanf(line.c_str(), "%f,%f,%f", &x, &y, &z);
        array.push_back({x,y,z});
    }
    
    DomId.resize(22420);
    std::iota(DomId.begin(), DomId.end(), 0);
    
    // Environment variable ${TRIDENT_QE_DIR} defines quantum efficiency;
    const char* qe_baseDir = std::getenv("TRIDENT_QE_DIR");
    if (qe_baseDir) {
        std::cout << "Quantum efficiency data directory: " << qe_baseDir << std::endl;
    } else {
        std::cerr << "Environment variable ${TRIDENT_QE_DIR} not set" << std::endl;
    }
    
    // Retrieve quantum efficiency
    auto gPMT = new TGraph(TString() + qe_baseDir + "/PMT_QE.csv", "%lg,%lg");
    // auto gSiPM= new TGraph(TString() + qe_baseDir + "/SiPM_PDE.csv", "%lg,%lg");
    pmtQE = new TSpline3("pmtQE", gPMT);
    // sipmPDE = new TSpline3("sipmPDE", gSiPM);
}

Event_processer::~Event_processer() {
    delete rnd_gen;
    delete pmtQE;
    // delete sipmPDE;
}

void Event_processer::set_up_rootfile(std::string rootfile) {
    fROOT = new TFile(rootfile.c_str(), "UPDATE");
    pmtHits = (TTree*) fROOT->Get("PmtHit");
    // sipmHits = (TTree*) fROOT->Get("SipmHit");
    primaryTree = (TTree*) fROOT->Get("Primary");
    
    
    pmtHits->SetBranchAddress("t0", &pmt_t);
    pmtHits->SetBranchAddress("e0", &pmt_e);
    pmtHits->SetBranchAddress("PmtId", &pmt_id);
    pmtHits->SetBranchAddress("DomId", &pmt_DOMid);
    // sipmHits->SetBranchAddress("t0", &sipm_t);
    // sipmHits->SetBranchAddress("e0", &sipm_e);
    // sipmHits->SetBranchAddress("SipmId", &sipm_id);
    // sipmHits->SetBranchAddress("DomId", &sipm_DOMid);

    primaryTree->SetBranchAddress("x0", &pri_x);
    primaryTree->SetBranchAddress("y0", &pri_y);
    primaryTree->SetBranchAddress("z0", &pri_z);
    primaryTree->SetBranchAddress("px", &pri_nx);
    primaryTree->SetBranchAddress("py", &pri_ny);
    primaryTree->SetBranchAddress("pz", &pri_nz);
    primaryTree->SetBranchAddress("PdgId", &pri_pdgid);

    hitsTree = new TTree("Hits", "Hits");
    hitsTree->Branch("t0", &new_t);
    hitsTree->Branch("PmtId", &new_pmt_id);
    hitsTree->Branch("DomId", &new_DOMid);
    hitsTree->Branch("x0", &new_x);
    hitsTree->Branch("y0", &new_y);
    hitsTree->Branch("z0", &new_z);
    hitsTree->Branch("Type", &new_hit_type);

    const int MAX_CAPACITY = 100000;
    new_t.reserve(MAX_CAPACITY);
    new_pmt_id.reserve(MAX_CAPACITY);
    new_DOMid.reserve(MAX_CAPACITY);
    new_x.reserve(MAX_CAPACITY); new_y.reserve(MAX_CAPACITY); new_z.reserve(MAX_CAPACITY);
    new_hit_type.reserve(MAX_CAPACITY);

}

void Event_processer::Loop_events() {
    const std::string k40 = config["K40"]["process"].as<std::string>(); 
    const std::string dark_noise = config["DarkNoise"]["process"].as<std::string>();
    const std::string event_type = config["event_type"].as<std::string>();
    const double range = config["range"].as<double>();
    const std::string mode = config["noise_mode"].as<std::string>();
    

    for(int iEvent=0; iEvent < pmtHits->GetEntries(); iEvent++){
        pmtHits->GetEntry(iEvent);
        // sipmHits->GetEntry(iEvent);

        primaryTree->GetEntry(iEvent);

        load_pmt();
        // load_sipm();
        if(new_t.size()==0) {
            hitsTree->Fill();
            new_t.clear();
            new_pmt_id.clear();
            new_DOMid.clear();
            new_x.clear();
            new_y.clear();
            new_z.clear();
            new_hit_type.clear();
            continue;
        }
        
        if(mode != "rage"){
            if(event_type=="track"){
                DomId.clear();
                
                auto it = std::find_if(pri_pdgid->begin(), pri_pdgid->end(), [](float i){return i==13 || i == -13;});
                if(it == pri_pdgid->end()) continue;
                int mu_index = std::distance(pri_pdgid->begin(), it);
                TVector3 n(pri_nx->at(mu_index), pri_ny->at(mu_index), pri_nz->at(mu_index));
                
                for(int i=0; i<22420;i++){
                    double lx = array.at(i).x() - pri_x->at(0);
                    double ly = array.at(i).y() - pri_y->at(0);
                    double lz = array.at(i).z() - pri_z->at(0);
                    TVector3 l(lx, ly, lz);
                    TVector3 cross = n.Cross(l);
                    if(cross.Mag() < range){
                        DomId.emplace_back(i);
                    }
                }
            }else if(event_type=="cascade"){
                DomId.clear();
                for(int i=0; i<22420;i++){
                    double lx = array.at(i).x() - pri_x->at(0);
                    double ly = array.at(i).y() - pri_y->at(0);
                    double lz = array.at(i).z() - pri_z->at(0);
                    if(TMath::Sqrt(lx*lx + ly*ly + lz*lz) < range){
                        DomId.emplace_back(i);
                    }
                }
            }
        }
        
        std::cout << "Event "<< iEvent << "/" << pmtHits->GetEntries() << "..." <<std::endl;
        float event_t_min = *std::min_element(new_t.begin(), new_t.end()) - 20;
        event_t_min = event_t_min > 0 ? event_t_min : 0;
        float event_t_max = *std::max_element(new_t.begin(), new_t.end()) + 40;
        float event_duration = event_t_max - event_t_min;
        if(k40=="on"){
            K40_processer k40processer;
            k40processer.process_K40(event_t_min,event_duration);
        }
        if((dark_noise=="both") || (dark_noise=="PMT")){
            process_pmt_darknoise(event_t_min, event_duration);
        }
        // if(dark_noise=="both" or dark_noise=="SiPM"){
            
        // }

        hitsTree->Fill();
        new_t.clear();
        new_pmt_id.clear();
        new_DOMid.clear();
        new_x.clear();
        new_y.clear();
        new_z.clear();
        new_hit_type.clear();
    }
    hitsTree->Write();
    fROOT->Close();
}

void Event_processer::load_pmt() {
    for(size_t iHit=0; iHit<pmt_t->size(); iHit++){
		if(pmt_t->at(iHit) > max_time) continue;
        double wave_length = (1240. / pmt_e->at(iHit));
        // check if pass quantum efficiency test
        if(rnd_gen->Uniform()*100 > pmtQE->Eval(wave_length)){
            new_t.emplace_back(pmt_t->at(iHit));
            new_pmt_id.emplace_back(pmt_id->at(iHit));
            auto domId = pmt_DOMid->at(iHit);
            new_DOMid.emplace_back(domId);
            new_x.emplace_back(array.at(domId).x());
            new_y.emplace_back(array.at(domId).y());
            new_z.emplace_back(array.at(domId).z());
            new_hit_type.emplace_back(1);
        }
    }
}

