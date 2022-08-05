#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include "TLorentzVector.h"

#include "Functions.h"

using namespace std;

class AK4Jet {

 public:
  
  int  jetid;
  int  jetid_tightlepveto;
  float  pt;
  float  eta;
  float  mass;
  float  phi;
  float  y;
  int    hadronFlavour;                                                                        
  int    partonFlavour;
  float  btag_DeepFlav;
  float  btag_DeepCSV;
  float  puid;
  float  qgl;
  bool   closebymu;
  bool   closebyel;
  float pt_resoup;
  float mass_resoup;
  float pt_resodn;
  float mass_resodn;
  float jesup_total;
  float jesdn_total;
  float btag_DeepFlav_tight_SF;
  float btag_DeepFlav_tight_SF_up;
  float btag_DeepFlav_tight_SF_dn;
  float btag_DeepFlav_med_SF;
  float btag_DeepFlav_med_SF_up;
  float btag_DeepFlav_med_SF_dn;
  float btag_DeepFlav_loose_SF;
  float btag_DeepFlav_loose_SF_up;
  float btag_DeepFlav_loose_SF_dn;
/*    float  reso;
    float  reso_up;
    float  reso_dn;
    float jesup_pu;
    float jesup_rel;
    float jesup_scale;
    float jesdn_pu;
    float jesdn_rel;
    float jesdn_scale;
*/
  //  int   genmatch;
  int AK8_neighbor_index;
  TLorentzVector p4;
};

class AK8Jet {

 public:
  
  int  jetid;
  int  jetid_tightlepveto;
  float  pt;
  float  eta;
  float  mass;
  float  phi;
  float  y;
  float pt_resoup;
  float mass_resoup;
  float pt_resodn;
  float mass_resodn;
  float jesup_total;
  float jesdn_total;
  float chrad;
  float tau21;
  float tau32;
  float DeepTag_TvsQCD;
  float DeepTag_WvsQCD;
  float DeepTag_ZvsQCD;
  float PNet_TvsQCD;
  float PNet_WvsQCD;
  float btag_DeepCSV;
  
  float CHF;
  float NHF;
  float CEMF;
  float NEMF;
  float MUF;
  float PHF;
  float NHadF;
  float EMF;

  int EEM;
  int neucons;
  int chcons;
  
  float neuemfrac;
  float neunhadfrac;

  float sdmass;
  float sub1pt;
  float sub1eta;
  float sub1phi;
  float sub1mass;
  float sub1btag;
  float sub1hadfrac;
  float sub1emfrac;
  float sub2pt;
  float sub2eta;
  float sub2phi;
  float sub2mass;
  float sub2btag;
  float sub2hadfrac;
  float sub2emfrac;
  float subbtag;
  float subhaddiff;
  float subemdiff;
  float subptdiff;
  
  float elinsubpt;
  float elinsubeta;
  float elinsubphi;
  float elinsubjpt;
  float elinsubjeta;
  float elinsubjphi;
  float elinsubjmass;
  float muinsubpt;
  float muinsubeta;
  float muinsubphi;
  float muinsubjpt;
  float muinsubjeta;
  float muinsubjphi;
  float muinsubjmass;
  float muinsubI0;
  float muinsubInear;
  float muinsubIfar;
  int match_muon_index;
  int match_electron_index;
  int match_lepton_index;
  int match_AK4_index;
  float matchAK4deepb;
  float re_tvsb;
  float rmu_tvsb;

  bool haselectron;
  bool hasmuon;
  bool hastau;
  bool hasqg; 
  bool hasb;
  bool hasleptop; 
  bool hashadtop;
  bool hastop;
  bool hasleptop_alldecay;
  bool hashadtop_alldecay;
  bool haspfelectron;
  bool haspfmuon;
  bool hasmatchmu;
  bool hasmatche;

  TLorentzVector p4;
};


class Muon {

 public:

  bool vetoid_tight;
  float pt;
  float eta;
  float phi;
  float y;
  float  mass; //need to ask Suman why he re-added it
  float charge;
  float trkvtx;
  float dz;
  bool ip;
  bool isTRK;
  bool isGL;
  bool isPF;
  bool isLoose;
  bool isGoodGL;
  bool isMed;
  bool isMedPr;
  bool isTight;
  bool isHighPt;
  bool isHighPttrk;
  float minisoall;
  float chi;
  float posmatch;
  float trkink;
  float segcom;
  float hit;
  float mst;
  float pixhit;
  float trklay;
  float valfrac;
  float pfiso;
  float p;
  float mudxy_sv;
  
  TLorentzVector p4;
};

class Electron {

 public:

  bool vetoid_tight;
  float pt;
  float eta;
  float phi;
  float  mass; //need to ask Suman why he re-added it
  //float charge;
  bool mvaWP90id;
  bool Fallv2WP80;
  bool mvaWP90id_noIso;
  bool Fallv2WP80_noIso;
  bool mvaWPloose;
  bool mvaWPloose_noIso;
  bool loosecutid;
  bool vetocutid;
  bool medcutid;
  bool tightcutid;
  float p;
  float y;
  float dxy;
  float dz;
  bool ip;
  float pfiso;
  float eldxy_sv;
  float supcl_eta;
  float supcl_phi;
  float supcl_rawE;
  float sigmaieta;
  float sigmaiphi;
  float r9full;
  float supcl_etaw;
  float supcl_phiw;
  float hcaloverecal;
  float cloctftrkn;
  float cloctftrkchi2;
  float e1x5bye5x5;
  float normchi2;
  float hitsmiss;
  float trkmeasure;
  float ecloverpout;
  float ecaletrkmomentum;
  float deltaetacltrkcalo;
  float supcl_preshvsrawe;
  float pfisolsumphet;
  float pfisolsumchhadpt;
  float pfsiolsumneuhadet;
  float etain;
  float phiin;
  float fbrem;
  float eoverp;
  float hovere;

  TLorentzVector p4;
};

class Lepton {
 
 public:

  float pt;
  float y;
  float eta;
  float phi;
  float mass; //need to ask Suman why he re-added it
  float charge;
  int lepton_id;
  int indexemu;
  int pdgId;
  int AK8_neighbor_index;
  
  TLorentzVector p4;
  
};

class AK4GenJet {

 public:

  float  eta;
  float  y;
  float  mass;
  float  phi;
  float  pt;
  int hadronFlavor;
  int partonFlavor;

  TLorentzVector p4;

};

class AK8GenJet {

 public:

  float  eta;
  float  y;
  float  mass;
  float  phi;
  float  pt;
  int hadronFlavor;
  int partonFlavor;

  TLorentzVector p4;

} ;

class GenParton{

 public:

  float  eta;
  float  y;
  float  mass;
  float  phi;
  float  pt;

  int status;
  int pdgId;
  int mompdgId;
  int grmompdgId;

  bool fromhard;
  bool fromhardbFSR;
  bool isPromptFinalState;
  bool isLastCopyBeforeFSR;

  TLorentzVector p4;
  
} ;

struct LHEparticle{

  float  eta;
  float  y;
  float  mass;
  float  phi;
  float  pt;
  int pdgId;

  TLorentzVector p4;

} ;

class TopQuark{
  // gives 4-momentum of top quark and a vector of its daughters                                    
  //(length of the vector of daughters should be 3)                                                 
  // each daughter is of GenParton-type                                                             
                       
 public:

  TLorentzVector p4;
  vector<GenParton> daughter;

} ;

struct Single_Trigger{
  bool single_hlts;
  float single_pt_cuts;
  int single_pids;
  float single_other_pt_cuts;
  int single_other_pids;
};

struct Double_Trigger{
    bool double_hlts;
    vector<float> double_pt_cuts;
    vector<int> double_pids;
};

struct Cuts{
  bool cut_pass;
  string name;
};

bool AK4Jet_sort_by_pt(AK4Jet i1, AK4Jet i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<AK4Jet> & objs) {
  sort(objs.begin(), objs.end(), AK4Jet_sort_by_pt);
}
bool AK8Jet_sort_by_pt(AK8Jet i1, AK8Jet i2)                                                        
{                                                                                                   
  return (i1.pt > i2.pt);                                                                           
}                                                                                                   
void sorted_by_pt(vector<AK8Jet> & objs) {                                                          
  sort(objs.begin(), objs.end(), AK8Jet_sort_by_pt);
}
bool AK4GenJet_sort_by_pt(AK4GenJet i1, AK4GenJet i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<AK4GenJet> & objs) {
  sort(objs.begin(), objs.end(), AK4GenJet_sort_by_pt);
}
bool Muon_sort_by_pt(Muon i1, Muon i2)                                                                           
{                                                                                                                
  return (i1.pt > i2.pt);                                                                                        
}                                                                                                                
void sorted_by_pt(vector<Muon> & objs) {                                                                         
  sort(objs.begin(), objs.end(), Muon_sort_by_pt);                                                               
}
bool Electron_sort_by_pt(Electron i1, Electron i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<Electron> & objs) {
  sort(objs.begin(), objs.end(), Electron_sort_by_pt);
}
bool Lepton_sort_by_pt(Lepton i1, Lepton i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<Lepton> & objs) {
  sort(objs.begin(), objs.end(), Lepton_sort_by_pt);
}
bool Parton_sort_by_pt(GenParton i1, GenParton i2)
{
  return (i1.pt > i2.pt);
}
void sorted_by_pt(vector<GenParton> & objs) {
  sort(objs.begin(), objs.end(), Parton_sort_by_pt);
}

bool Top_sort_by_pt(TopQuark i1, TopQuark i2)
{
  return (i1.p4.Pt() > i2.p4.Pt());
}
void sorted_by_pt(vector<TopQuark> & objs) {
  sort(objs.begin(), objs.end(), Top_sort_by_pt);
}

float compute_HT(vector<AK4Jet>  & objs, float ptcut, float etacut){
  float HT = 0;
  for(unsigned iobs=0; iobs<objs.size(); iobs++){
    if(objs[iobs].pt > ptcut && abs(objs[iobs].eta)<=etacut){
      HT += objs[iobs].pt;
    }
  }
  return HT;
}


int get_nearest_AK4(vector<AK4Jet>  & objs, TLorentzVector tmp_vec, float minR = 0.8/*0.6*/) {
  // gives the index of AK4 jet nearest to the tmp_vec vector                                                       
  int nearest = -1;
  float btag_score_max=-1;
  for(unsigned iobs=0; iobs<objs.size(); iobs++){
    if(objs[iobs].AK8_neighbor_index>=0 ) continue;
    if(delta2R(objs[iobs].p4,tmp_vec) < minR && objs[iobs].btag_DeepFlav > btag_score_max){
      btag_score_max = objs[iobs].btag_DeepFlav;
      nearest = iobs;
    }
  }
    if(nearest >=0) objs[nearest].AK8_neighbor_index = +1;
  return  nearest;
}

int get_nearest_lepton(vector<Lepton>  & objs, TLorentzVector tmp_vec, int lepton_pdgid=-1, float minR = 0.7) {
  // gives the index of lepton nearest to the tmp_vec vector                                  

  int nearest = -1;
  float ptmax = -99;
  for(unsigned iobs=0; iobs<objs.size(); iobs++){
    if(lepton_pdgid>=0 && objs[iobs].pdgId!=lepton_pdgid) continue;
    if(objs[iobs].AK8_neighbor_index>=0 ) continue;
    float tmpR = delta2R(objs[iobs].p4,tmp_vec);
    if( tmpR < minR && objs[iobs].pt > ptmax){
      ptmax = objs[iobs].pt;
      nearest = iobs;
    }
  }
 /* TString str;
 str = TString::Format("check %i near %i  %f",int(objs.size()),nearest,minR);
	if(gProofServ) gProofServ->SendAsynMessage(str);*/
  if(nearest >=0) objs[nearest].AK8_neighbor_index = +1;
  return  nearest;
}
int get_nearest_AK8Jet(vector<AK8Jet>  & objs, TLorentzVector tmp_vec, float minR = 0.7/*0.8*/) {
  // gives the index of AK8 jet nearest to the tmp_vec vector                                  
                         
  int nearest = -1;
  for(unsigned iobs=0; iobs<objs.size(); iobs++){
    if(delta2R(objs[iobs].p4,tmp_vec) < minR){
      minR = delta2R(objs[iobs].p4,tmp_vec) ;
      nearest = iobs;
    }
  }
  return  nearest;
}

/*
  bool AK4Jet_sort_by_DeepFlav(AK4Jet i1, AK4Jet i2)
  {
  return (i1.btagDeepFlavB > i2.btagDeepFlavB);
  }
  void sorted_by_DeepFlav(vector<AK4Jet> & objs) {
  sort(objs.begin(), objs.end(), AK4Jet_sort_by_DeepFlav);
  }
  bool AK8Jet_sort_by_DeepAK8_Htag(AK8Jet i1, AK8Jet i2)
  {
  return (i1.deepTagMD_bbvsLight > i2.deepTagMD_bbvsLight);
  }
  void sorted_by_DeepAK8_Htag(vector<AK8Jet> & objs) {
  sort(objs.begin(), objs.end(), AK8Jet_sort_by_DeepAK8_Htag);
  }
  bool GenAK4Jet_sort_by_pt(AK4GenJet i1, AK4GenJet i2)
  {
  return (i1.pt > i2.pt);
  }
  void sorted_by_pt(vector<AK4GenJet> & objs) {
  sort(objs.begin(), objs.end(), GenAK4Jet_sort_by_pt);
  }
  bool GenAK8Jet_sort_by_pt(AK8GenJet i1, AK8GenJet i2)
  {
  return (i1.pt > i2.pt);
  }
  void sorted_by_pt(vector<AK8GenJet> & objs) {
  sort(objs.begin(), objs.end(), GenAK8Jet_sort_by_pt);
  }
*/


