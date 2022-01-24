using namespace std;


//////////////////////////////////////////////////////////
// variables of the raw root file tree (siwecaldecoded)
// A. Irles 2021

//pedestal and mapping variables

Float_t map_pointX[15][16][64];
Float_t map_pointY[15][16][64];
Int_t masked[15][16][64];
std::vector<std::vector<std::vector<Double_t> > > ped_mean;
std::vector<std::vector<std::vector<Double_t> > > ped_error;
std::vector<std::vector<std::vector<Double_t> > > ped_width;

std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_mean_slboard;
std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_error_slboard;
std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_width_slboard;

std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_mean_cov_slboard;
std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_w_i_slboard;
std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_w_c1_slboard;
std::vector<std::vector<std::vector<std::vector<Double_t> > > > ped_w_c2_slboard;

//mapping between slboard add and position
int slboard_array_mapping[15]={-1};


// Declaration of leaf types
Int_t           event;
Int_t           acqNumber;
Int_t           n_slboards;
Int_t           slot[15];
Int_t           slboard_id[15];
Int_t           chipid[15][16];
Int_t           nColumns[15][16];
Int_t           startACQ[15];
Int_t           rawTSD[15];
Float_t         TSD[15];
Int_t           rawAVDD0[15];
Int_t           rawAVDD1[15];
Float_t         AVDD0[15];
Float_t         AVDD1[15];
Int_t           bcid[15][16][15];
Int_t           corrected_bcid[15][16][15];
Int_t           badbcid[15][16][15];
Int_t           nhits[15][16][15];
Int_t           charge_lowGain[15][16][15][64];
Int_t           charge_hiGain[15][16][15][64];
Int_t           gain_hit_low[15][16][15][64];
Int_t           gain_hit_high[15][16][15][64];

// List of branches
TBranch        *b_event;   //!
TBranch        *b_acqNumber;   //!
TBranch        *b_n_slboards;   //!
TBranch        *b_slot;   //!
TBranch        *b_slboard_id;   //!
TBranch        *b_chipid;   //!
TBranch        *b_nColumns;   //!
TBranch        *b_startACQ;   //!
TBranch        *b_rawTSD;   //!
TBranch        *b_TSD;   //!
TBranch        *b_rawAVDD0;   //!
TBranch        *b_rawAVDD1;   //!
TBranch        *b_AVDD0;   //!
TBranch        *b_AVDD1;   //!
TBranch        *b_bcid;   //!
TBranch        *b_corrected_bcid;   //!
TBranch        *b_badbcid;   //!
TBranch        *b_nhits;   //!
TBranch        *b_lowGain;   //!
TBranch        *b_highGain;   //!
TBranch        *b_gain_hit_low;   //!
TBranch        *b_gain_hit_high;   //!


void SetBranchAddressFunction(TTree *fChain) {
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("acqNumber", &acqNumber, &b_acqNumber);
  fChain->SetBranchAddress("n_slboards", &n_slboards, &b_n_slboards);
  fChain->SetBranchAddress("slot", slot, &b_slot);
  fChain->SetBranchAddress("slboard_id", slboard_id, &b_slboard_id);
  fChain->SetBranchAddress("chipid", chipid, &b_chipid);
  fChain->SetBranchAddress("nColumns", nColumns, &b_nColumns);
  fChain->SetBranchAddress("startACQ", startACQ, &b_startACQ);
  fChain->SetBranchAddress("rawTSD", rawTSD, &b_rawTSD);
  fChain->SetBranchAddress("TSD", TSD, &b_TSD);
  fChain->SetBranchAddress("rawAVDD0", rawAVDD0, &b_rawAVDD0);
  fChain->SetBranchAddress("rawAVDD1", rawAVDD1, &b_rawAVDD1);
  fChain->SetBranchAddress("AVDD0", AVDD0, &b_AVDD0);
  fChain->SetBranchAddress("AVDD1", AVDD1, &b_AVDD1);
  fChain->SetBranchAddress("bcid", bcid, &b_bcid);
  fChain->SetBranchAddress("corrected_bcid", corrected_bcid, &b_corrected_bcid);
  fChain->SetBranchAddress("badbcid", badbcid, &b_badbcid);
  fChain->SetBranchAddress("nhits", nhits, &b_nhits);
  fChain->SetBranchAddress("charge_lowGain", charge_lowGain, &b_lowGain);
  fChain->SetBranchAddress("charge_hiGain", charge_hiGain, &b_highGain);
  fChain->SetBranchAddress("gain_hit_low", gain_hit_low, &b_gain_hit_low);
  fChain->SetBranchAddress("gain_hit_high", gain_hit_high, &b_gain_hit_high);
  // Notify();
}
