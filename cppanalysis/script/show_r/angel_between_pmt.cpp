//
// Created by ineffablord on 23-11-28.
//
#include "TH1.h"
#include "fstream"
#include "iostream"
#include "sstream"
#include "vector"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "random"
#include "analysis_trigger_root_tts.h"
#include "cmath"
#include "set"

Double_t GetTotalTime(TTree* aintree){
    Long64_t hitnum;
    Double_t TotalTime;
    hitnum = aintree->GetEntries();
    TBranch* hittimebranch = aintree->GetBranch("hittime");
    hittimebranch->SetAddress(&TotalTime);
    aintree->GetEntry(hitnum - 1);
    return TotalTime;
}

void ApplyQE(TTree* aintree, const std::vector<QEdata>& qe, TTree* anewtree, std::mt19937 agenerator){
    //calculate wavelength
    Double_t energy;
    TBranch* energybranch = aintree->GetBranch("energy");
    energybranch->SetAddress(&energy);
    Double_t wavelength;
    TBranch* wavelengthbranch = aintree->Branch("wavelength", &wavelength, "wavelength/D");
    for (Long64_t i = 0; i < aintree->GetEntries(); ++i){
        aintree->GetEntry(i);
        wavelength = 1240 / energy;
        wavelengthbranch->Fill();
    }
    aintree->Write("",TObject::kOverwrite);

    anewtree->AddFriend(aintree);
    TBranch* eventid = aintree->GetBranch("eventid");
    Double_t eid;
    eventid->SetAddress(&eid);
    anewtree->Branch("eventid",&eid,"eventid/I");
    TBranch* hittime = aintree->GetBranch("hittime");
    Double_t hitt;
    hittime->SetAddress(&hitt);
    anewtree->Branch("hittime",&hitt,"hittime/D");
    TBranch* PMTid = aintree->GetBranch("PMTid");
    Double_t pmtid;
    PMTid->SetAddress(&pmtid);
    anewtree->Branch("PMTid",&pmtid,"PMTid/D");
    Double_t x,y,z,r;
    TBranch* xb = aintree->GetBranch("x");
    xb->SetAddress(&x);
    TBranch* yb = aintree->GetBranch("y");
    yb->SetAddress(&y);
    TBranch* zb = aintree->GetBranch("z");
    zb->SetAddress(&z);
    anewtree->Branch("r",&r,"r/D");
    std::uniform_real_distribution<double> uniform(0,1);
    for (Long64_t i = 0; i < aintree->GetEntries(); ++i){
        aintree->GetEntry(i);
        r = sqrt(x*x + y*y + z*z);
        Double_t QE = interpolate(qe,wavelength) / 100;
        Double_t probability = uniform(agenerator);
        if (QE >= probability){
            anewtree->Fill();
        }
    }
    anewtree->Write("",TObject::kOverwrite);
}

std::vector<QEdata> Readcsv(const std::string& csvfile) {
    std::vector<QEdata> qedata;
    std::ifstream file(csvfile);
    std::string line;
    //discard the first line
    if (std::getline(file, line)){}
    while (std::getline(file, line)){
        std::istringstream ss(line);
        std::string token;
        QEdata point{};
        std::getline(ss, token, ',');
        std::getline(ss, token, ',');
        point.wl = std::stod(token);
        std::getline(ss, token);
        point.qe = std::stod(token);
        qedata.push_back(point);
    }
    return qedata;
}

Double_t interpolate(const std::vector<QEdata>& qe, Double_t x){
    for (size_t i = 0; i < qe.size() - 1; i++){
        if (x >= qe[i].wl && x <= qe[i +1].wl){
            Double_t x0 = qe[i].wl;
            Double_t x1 = qe[i+1].wl;
            Double_t y0 = qe[i].qe;
            Double_t y1 = qe[i+1].qe;
            Double_t slope = (y1 - y0)/(x1 - x0);
            Double_t intery = y0 + slope*(x - x0);
            return intery;
        }
    }
    if (x <= qe.front().wl){
        return qe.front().qe;
    }
    else if (x >= qe.back().wl){
        return qe.back().qe;
    }
    else return 0;
}

void Addtts(const Double_t tts, TTree* anewtree, std::mt19937 agenerator, Double_t t_win) {
    std::normal_distribution<double> Normal(0,tts/2.3548);
    Double_t hitt;
    TBranch* hittime = anewtree->GetBranch("hittime");
    hittime->SetAddress(&hitt);
    Double_t atts;
    TBranch* addtts = anewtree->Branch("addtts", &atts, "addtts/D");
    Double_t win;
    TBranch* windowbranch = anewtree->Branch("window", &win, "window/D");
    for (Long64_t i = 0; i < hittime->GetEntries(); i++){
        hittime->GetEntry(i);
        atts = hitt + Normal(agenerator);
        win = floor(atts / t_win);
        addtts->Fill();
        windowbranch->Fill();
    }
    anewtree->Write("",TObject::kOverwrite);
}

void CountCoinPmt(TTree* thenewtree, const Int_t coin, std::ofstream* csvfile, std::ofstream* rfile){
    Double_t theta[6]={55,73,107,125,150,180};
    Double_t phi[6][6]={{60,120,180,240,300,360},{30,90,150,210,270,330},{60,120,180,240,300,360},{30,90,150,210,270,330},{60,120,180,240,300,360},{0,0,0,0,0,0}};
    Double_t win,winf;
    TBranch* winbranch = thenewtree->GetBranch("window");
    winbranch->SetAddress(&win);
    Double_t pmtid;
    TBranch* pmtidbranch = thenewtree->GetBranch("PMTid");
    pmtidbranch->SetAddress(&pmtid);
    Double_t r;
    TBranch* rbranch = thenewtree->GetBranch("r");
    rbranch->SetAddress(&r);
    const Long64_t entries = thenewtree->GetEntries();
    Double_t count = 0;
    Long64_t x = 0;
    for ( ; x < entries - coin + 1; x++){
        winbranch->GetEntry(x);
        pmtidbranch->GetEntry(x);
        rbranch->GetEntry(x);
        winf = win;
        winbranch->GetEntry(x+coin-1);
        if (win == winf){
            std::set<Double_t> pmts;
            Int_t number_of_pmts = 0;
            Long64_t I = 0;
            while (win == winf and x+I < entries){
                if (pmts.count(pmtid) == 0){
                    number_of_pmts++;
                }
                pmts.insert(pmtid);
                I++;
                winbranch->GetEntry(x+I);
                pmtidbranch->GetEntry(x+I);
            }
            if (number_of_pmts==coin){
                count++;
                rbranch->GetEntry(x);
                *rfile<<r<<",";
                for (auto i = pmts.begin();i!=pmts.end();++i){
                    *rfile<<*i<<",";
                }
                x+=(I-1);
                Double_t cospsimin = 1;
                Double_t cospsimax = -1;
                for (auto i = pmts.begin(); i!=pmts.end(); ++i){
                    Double_t theta1 = theta[int(floor(*i/6))]*3.1415926/180;
                    Double_t phi1 = phi[int(floor(*i/6))][int(*i)%6]*3.1415926/180;
                    for (auto j = i; j!=pmts.end(); ++j){
                        Double_t theta2 = theta[int(floor(*j/6))]*3.1415926/180;
                        Double_t phi2 = phi[int(floor(*j/6))][int(*j)%6]*3.1415926/180;
                        Double_t _cospsi;
                        Double_t a,b,c;
                        a = pow((sin(theta1)* cos(phi1)- sin(theta2)* cos(phi2)),2);
                        b = pow((sin(theta1)* sin(phi1)- sin(theta2)* sin(phi2)),2);
                        c = pow((cos(theta1)- cos(theta2)),2);
                        _cospsi = (2-a-b-c)/2;
                        if(_cospsi<=cospsimin){cospsimin=_cospsi;}
                        if(_cospsi>=cospsimax and j!=i){cospsimax=_cospsi;}
                    }
                }
                Double_t angmax,angmin;
                angmax = acos(cospsimin)*180/3.1415926;
                angmin = acos(cospsimax)*180/3.1415926;
                std::cout<<cospsimax<<std::endl;
                *rfile<<angmin<<","<<angmax<<std::endl;
            }else{}
        }
        else{
        }
    }
    *csvfile << int(count) << ",";
}

int main(int argc, char** argv){
    if (argc != 2) {
        std::cerr << "Wrong input arguments! The right way is: ./main [path/to/root file]";
        return -1;
    }

    //const Int_t totalevent = 20000000;

    const Char_t* ROOTPath = argv[1];
    const Double_t TTS = 1.3;
    std::random_device rd;
    std::mt19937 generator(rd());

    const double time_win[] = {20};
    const int coin_num[] = {1,2,3,4,5,6,7};
    //open ROOT file. Tree:Hit0 Branch:eventid hittime PMTid energy x y z decaytime
    auto *infile= new TFile(ROOTPath, "UPDATE");
    TTree *intree;
    intree = dynamic_cast<TTree *>(infile->Get("Hit0"));
    TTree* newtree;
    newtree = new TTree("analysis", "fortrigger");
    newtree->SetDirectory(infile);
    //---Get total simulation time
    Double_t sim_time;
    sim_time = GetTotalTime(intree);
    std::cout << sim_time << std::endl;
    //---Get qe
    std::string qepath = "/lustre/neutrino/huangweilun/k40sim/cppanalysis/store/qe.csv";
    std::vector<QEdata> QE = Readcsv(qepath);
    //Apply qe in a new tree
    ApplyQE(intree,QE,newtree,generator);
    //Now the ROOT file should be: Tree:Hit0 Branch:eventid hittime PMTid energy x y z decaytime wavelength
    //                             Tree:analysis Branch:eventid hittime PMTid
    //---Add tts
    Addtts(TTS,newtree,generator,time_win[0]);
    //Now the ROOT file should be: Tree:Hit0 Branch:eventid hittime PMTid energy x y z decaytime wavelength
    //                             Tree:analysis Branch:eventid hittime PMTid addtts
    //---Trigger
    std::ofstream outcsv;
    outcsv.open("/lustre/neutrino/huangweilun/k40sim/k40data/multiruns/new_opt_sim_sak_500/ana_bin3_r/num.csv", std::ios_base::app);

    for (double j : time_win){
        for (int i : coin_num){
            std::ofstream rcsv;
            std::string rpath = "/lustre/neutrino/huangweilun/k40sim/k40data/multiruns/new_opt_sim_sak_500/ana_bin3_r/" + std::to_string(i) + "_" + "r.csv";
            rcsv.open(rpath, std::ios_base::app);
            CountCoinPmt(newtree,i,&outcsv,&rcsv);
            rcsv.close();
        }
    }

    outcsv<<std::endl;
    outcsv.close();
    infile->Write();
    infile->Close();
    delete infile;
    return 0;
}