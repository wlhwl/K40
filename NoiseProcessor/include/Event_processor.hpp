#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"
#include "TSpline.h"
#include "TRandom3.h"
#include "yaml-cpp/yaml.h"

inline YAML::Node config;

inline double max_time;
inline std::vector<float> *pmt_t= nullptr, *pmt_e= nullptr/*, *sipm_t= nullptr, *sipm_e= nullptr*/;
inline std::vector<int> *pmt_id= nullptr, *pmt_DOMid= nullptr/*, *sipm_id= nullptr, *sipm_DOMid= nullptr*/;

inline std::vector<double> new_t;
inline std::vector<int> new_pmt_id, new_DOMid; // 0~NUM_PMTs: PMT; NUM_PMTs~: SiPM
inline std::vector<float> new_x, new_y, new_z;
inline std::vector<int> new_hit_type; //0 for noise, 1 for signal, -1 for K40

inline std::vector<TVector3> array;
inline std::vector<int> DomId;
inline TRandom3* rnd_gen = new TRandom3(0);

class Event_processor {
    public:
    Event_processor(std::string config_file);
    ~Event_processor();

    void set_up_rootfile(std::string rootfile);
    void Loop_events();
    void load_pmt();
    // void load_sipm();

    
    TSpline3 *pmtQE/*,*sipmPDE*/;

    private:
    TFile* fROOT;
    TTree* pmtHits;
    // TTree* sipmHits;
    TTree* hitsTree;
    TTree* primaryTree;
    std::vector<float> *pri_x=nullptr, *pri_y=nullptr, *pri_z=nullptr;
    std::vector<float> *pri_nx=nullptr, *pri_ny=nullptr, *pri_nz=nullptr, *pri_pdgid=nullptr;

};
