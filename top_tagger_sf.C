/*
  weight = event_weight; 
  should be in the same ntuple, event by event we are going to read it.
    
  Order lepton according to the matching of leading Ak8 jets. 
  Lepton tagging
  
*/

#define top_tagger_sf_cxx
//#include "top_tagger_sf.h"
#include "getobjects.h"

#include <TH2.h>
#include <TStyle.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <fstream>
#include <TProofOutputFile.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TProofServ.h>


void top_tagger_sf::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();

}
 
void top_tagger_sf::SlaveBegin(TTree * /*tree*/)
{
  //The SlaveBegin() function is called after the Begin() function.
  //When running with PROOF SlaveBegin() is called on each slave server.
  //The tree argument is deprecated (on PROOF 0 is passed).
  
  TString option = GetOption();
  
  OutFile = new TProofOutputFile("test_output.root");
  
  fileOut = OutFile->OpenFile("RECREATE");
  if ( !(fileOut = OutFile->OpenFile("RECREATE")) )
    {
      Warning("SlaveBegin", "problems opening file: %s/%s",OutFile->GetDir(), OutFile->GetFileName());
    }
   
  isMC = true;
  isTT = true;
  isST = false;
  isDIB = false;
  isWJ = false;
  isDY = false;
  isQCD = false;
  
  isTTX = false;
  isTTsignal = true;

  wp_leptop_tagger[0] = 0.8;
  wp_leptop_tagger[1] = 0.3;

  //GMA
  //Tout = new TTree("leptop","leptop");
  //Tout->Branch("weight",&weight,"weight/F");
  
  Tnewvar = new TTree("newvars","newvars");
  Tnewvar->Branch("irun", &irun, "irun/I");
  Tnewvar->Branch("ilumi", &ilumi, "ilumi/I");
  Tnewvar->Branch("ievt", &ievt, "ievt/i");
  Tnewvar->Branch("weight",&weight,"weight/F");
  Tnewvar->Branch("nmu_veto",&nmu_veto,"nmu_veto/F");
  Tnewvar->Branch("nel_veto",&nel_veto,"nel_veto/F");
  Tnewvar->Branch("weight_mutop",&weight_mutop,"weight_mutop/F");
  Tnewvar->Branch("weight_eltop",&weight_eltop,"weight_eltop/F");
  Tnewvar->Branch("hadtaggersf",&hadtaggersf,"hadtaggersf/F");
  Tnewvar->Branch("chs_met",&chs_met,"chs_met/F");
  Tnewvar->Branch("deltaphi",&deltaphi,"deltaphi/F");
  Tnewvar->Branch("nbjets",&nbjets,"nbjets/F");
  Tnewvar->Branch("tagjetpt",&tagjetpt,"tagjetpt/F");
  Tnewvar->Branch("tagjet_taggerscore",&tagjet_taggerscore,"tagjet_taggerscore/F");
  Tnewvar->Branch("genmatch_tag",&genmatch_tag,"genmatch_tag/O");
  Tnewvar->Branch("genmatch_probe",&genmatch_probe,"genmatch_probe/O");
  Tnewvar->Branch("topjetmatch",&topjetmatch,"topjetmatch/O");
  Tnewvar->Branch("elmatch",&elmatch,"elmatch/O");
  Tnewvar->Branch("mumatch",&mumatch,"mumatch/O");
  Tnewvar->Branch("probejet_score",&probejet_score,"probejet_score/F");
  Tnewvar->Branch("nelectrons",&nelectrons_eltop,"nelectrons/F");
  Tnewvar->Branch("nmuons",&nmuons_mutop,"nmuons/F"); 
  Tnewvar->Branch("bjetmatch",&bjetmatch,"bjetmatch/O");
  Tnewvar->Branch("dirgltrthr",&dirgltrthr,"dirgltrthr/F");
  Tnewvar->Branch("dirglthrmin",&dirglthrmin,"dirglthrmin/F");

  //Tnewvar->Branch("",&,"/F");
  
  hist2d_toptag_sf[0][0] = new TH2D("hist2d_tagpass_eltop","pt of el top jet vs eta of el top after tag jet pass cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[1][0] = new TH2D("hist2d_genmatch_tagpass_eltop","pt of el top jet vs eta of el top after tag jet pass and tag-gen match cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[2][0] = new TH2D("hist2d_probepass_eltop","pt of el top jet vs eta of el top after tag & probe jet pass cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[3][0] = new TH2D("hist2d_genmatch_probepass_eltop","pt of el top jet vs eta of el top after tag & probe jet pass and tag-gen match cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[4][0] = new TH2D("hist2d_genmatch_tagprobepass_eltop","pt of el top jet vs eta of el top after tag & probe jet pass and tag-gen & probe-gen match cut",nptbin,ptbinarray,netabin,etabinarray);
  
  hist2d_toptag_sf[0][1] = new TH2D("hist2d_tagpass_mutop","pt of mu top jet vs eta of mu top after tag jet pass cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[1][1] = new TH2D("hist2d_genmatch_tagpass_mutop","pt of mu top jet vs eta of mu top after tag jet pass and tag-gen match cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[2][1] = new TH2D("hist2d_probepass_mutop","pt of mu top jet vs eta of mu top after tag & probe jet pass cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[3][1] = new TH2D("hist2d_genmatch_probepass_mutop","pt of mu top jet vs eta of mu top after tag & probe jet pass and tag-gen match cut",nptbin,ptbinarray,netabin,etabinarray);
  hist2d_toptag_sf[4][1] = new TH2D("hist2d_genmatch_tagprobepass_mutop","pt of mu top jet vs eta of mu top after tag & probe jet pass and tag-gen & probe-gen match cut",nptbin,ptbinarray,netabin,etabinarray);
  
  hist_ptak8_sf[0][0] = new TH1D("hist_ptak8_tagpass_eltop","pt of el top jet after tag jet pass cut",nptbin,ptbinarray);
  hist_ptak8_sf[1][0] = new TH1D("hist_ptak8_genmatch_tagpass_eltop","pt of el top jet after tag jet pass and tag-gen match cut",nptbin,ptbinarray);
  hist_ptak8_sf[2][0] = new TH1D("hist_ptak8_probepass_eltop","pt of el top jet after tag & probe jet pass cut",nptbin,ptbinarray);
  hist_ptak8_sf[3][0] = new TH1D("hist_ptak8_genmatch_probepass_eltop","pt of el top jet after tag & probe jet pass and tag-gen match cut",nptbin,ptbinarray);
  hist_ptak8_sf[4][0] = new TH1D("hist_ptak8_genmatch_tagprobepass_eltop","pt of el top jet after tag & probe jet pass and tag-gen & probe-gen match cut",nptbin,ptbinarray);
  
  hist_ptak8_sf[0][1] = new TH1D("hist_ptak8_tagpass_mutop","pt of mu top jet after tag jet pass cut",nptbin,ptbinarray);
  hist_ptak8_sf[1][1] = new TH1D("hist_ptak8_genmatch_tagpass_mutop","pt of mu top jet after tag jet pass and tag-gen match cut",nptbin,ptbinarray);
  hist_ptak8_sf[2][1] = new TH1D("hist_ptak8_probepass_mutop","pt of mu top jet after tag & probe jet pass cut",nptbin,ptbinarray);
  hist_ptak8_sf[3][1] = new TH1D("hist_ptak8_genmatch_probepass_mutop","pt of mu top jet after tag & probe jet pass and tag-gen match cut",nptbin,ptbinarray);
  hist_ptak8_sf[4][1] = new TH1D("hist_ptak8_genmatch_tagprobepass_mutop","pt of mu top jet after tag & probe jet pass and tag-gen & probe-gen match cut",nptbin,ptbinarray);

  hist_etaak8_sf[0][0] = new TH1D("hist_etaak8_tagpass_eltop","eta of el top after tag jet pass cut",12,-2.5,2.5);
  hist_etaak8_sf[1][0] = new TH1D("hist_etaak8_genmatch_tagpass_eltop","eta of el top after tag jet pass and tag-gen match cut",12,-2.5,2.5);
  hist_etaak8_sf[2][0] = new TH1D("hist_etaak8_probepass_eltop","eta of el top after tag & probe jet pass cut",12,-2.5,2.5);
  hist_etaak8_sf[3][0] = new TH1D("hist_etaak8_genmatch_probepass_eltop","eta of el top after tag & probe jet pass and tag-gen match cut",12,-2.5,2.5);
  hist_etaak8_sf[4][0] = new TH1D("hist_etaak8_genmatch_tagprobepass_eltop","eta of el top after tag & probe jet pass and tag-gen & probe-gen match cut",12,-2.5,2.5);
  
  hist_etaak8_sf[0][1] = new TH1D("hist_etaak8_tagpass_mutop","eta of mu top after tag jet pass cut",12,-2.5,2.5);
  hist_etaak8_sf[1][1] = new TH1D("hist_etaak8_genmatch_tagpass_mutop","eta of mu top after tag jet pass and tag-gen match cut",12,-2.5,2.5);
  hist_etaak8_sf[2][1] = new TH1D("hist_etaak8_probepass_mutop","eta of mu top after tag & probe jet pass cut",12,-2.5,2.5);
  hist_etaak8_sf[3][1] = new TH1D("hist_etaak8_genmatch_probepass_mutop","eta of mu top after tag & probe jet pass and tag-gen & probe-gen match cut",12,-2.5,2.5);
  hist_etaak8_sf[4][1] = new TH1D("hist_etaak8_genmatch_probepass_mutop","eta of mu top after tag & probe jet pass and tag-gen & probe-gen match cut",12,-2.5,2.5);

  hist_sdmassak8_sf[0][0] = new TH1D("hist_sdmassak8_tagpass_eltop","sdmass of el top after tag jet pass cut",25,0,250);
  hist_sdmassak8_sf[1][0] = new TH1D("hist_sdmassak8_genmatch_tagpass_eltop","sdmass of el top after tag jet pass and tag-gen match cut",25,0,250);
  hist_sdmassak8_sf[2][0] = new TH1D("hist_sdmassak8_probepass_eltop","sdmass of el top after tag & probe jet pass cut",25,0,250);
  hist_sdmassak8_sf[3][0] = new TH1D("hist_sdmassak8_genmatch_probepass_eltop","sdmass of el top after tag & probe jet pass and tag-gen match cut",25,0,250);
  hist_sdmassak8_sf[4][0] = new TH1D("hist_sdmassak8_genmatch_tagprobepass_eltop","sdmass of el top after tag & probe jet pass and tag-gen & probe-gen match cut",25,0,250);
  
  hist_sdmassak8_sf[0][1] = new TH1D("hist_sdmassak8_tagpass_mutop","sdmass of mu top after tag jet pass cut",25,0,250);
  hist_sdmassak8_sf[1][1] = new TH1D("hist_sdmassak8_genmatch_tagpass_mutop","sdmass of mu top after tag jet pass and tag-gen match cut",25,0,250);
  hist_sdmassak8_sf[2][1] = new TH1D("hist_sdmassak8_probepass_mutop","sdmass of mu top after tag & probe jet pass cut",25,0,250);
  hist_sdmassak8_sf[3][1] = new TH1D("hist_sdmassak8_genmatch_probepass_mutop","sdmass of mu top after tag & probe jet pass and tag-gen match cut",25,0,250);
  hist_sdmassak8_sf[4][1] = new TH1D("hist_sdmassak8_genmatch_tagprobepass_mutop","sdmass of mu top after tag & probe jet pass and tag-gen & probe-gen match cut",25,0,250);

  hist_sdmasstagak8_sf[0][0] = new TH1D("hist_tagsdmassak8_tagpass_eltop","sdmass of el top in e-channel after tag jet pass cut",25,0,250);
  hist_sdmasstagak8_sf[1][0] = new TH1D("hist_tagsdmassak8_genmatch_tagpass_eltop","sdmass of el top in e-channel after tag jet pass and tag-gen match cut",25,0,250);
  hist_sdmasstagak8_sf[2][0] = new TH1D("hist_tagsdmassak8_probepass_eltop","sdmass of el top in e-channel after tag & probe jet pass cut",25,0,250);
  hist_sdmasstagak8_sf[3][0] = new TH1D("hist_tagsdmassak8_genmatch_probepass_eltop","sdmass of el top in e-channel after tag & probe jet pass and tag-gen match cut",25,0,250);
  hist_sdmasstagak8_sf[4][0] = new TH1D("hist_tagsdmassak8_genmatch_tagprobepass_eltop","sdmass of el top in e-channel after tag & probe jet pass and tag-gen & probe-gen match cut",25,0,250);
  
  hist_sdmasstagak8_sf[0][1] = new TH1D("hist_tagsdmassak8_tagpass_mutop","sdmass of mu top in #mu-channel after tag jet pass cut",25,0,250);
  hist_sdmasstagak8_sf[1][1] = new TH1D("hist_tagsdmassak8_genmatch_tagpass_mutop","sdmass of mu top in #mu-channel after tag jet pass and tag-gen match cut",25,0,250);
  hist_sdmasstagak8_sf[2][1] = new TH1D("hist_tagsdmassak8_probepass_mutop","sdmass of mu top in #mu-channel after tag & probe jet pass cut",25,0,250);
  hist_sdmasstagak8_sf[3][1] = new TH1D("hist_tagsdmassak8_genmatch_probepass_mutop","sdmass of mu top in #mu-channel after tag & probe jet pass and tag-gen match cut",25,0,250);
  hist_sdmasstagak8_sf[4][1] = new TH1D("hist_tagsdmassak8_genmatch_tagprobepass_mutop","sdmass of mu top in #mu-channel after tag & probe jet pass and tag-gen & probe-gen match cut",25,0,250);
 
  hist_pnettagak8_sf[0][0] = new TH1D("hist_pnettagak8_tagpass_eltop","PNet tvsQCD score of tag jet in e-channel after tag jet pass cut",20,0.96,1);
  hist_pnettagak8_sf[1][0] = new TH1D("hist_pnettagak8_genmatch_tagpass_eltop","PNet tvsQCD score of tag jet in e-channel after tag jet pass and tag-gen match cut",20,0.96,1);
  hist_pnettagak8_sf[2][0] = new TH1D("hist_pnettagak8_probepass_eltop","PNet tvsQCD score of tag jet in e-channel after tag & probe jet pass cut",20,0.96,1);
  hist_pnettagak8_sf[3][0] = new TH1D("hist_pnettagak8_genmatch_probepass_eltop","PNet tvsQCD score of tag jet in e-channel after tag & probe jet pass and tag-gen match cut",20,0.96,1);
  hist_pnettagak8_sf[4][0] = new TH1D("hist_pnettagak8_genmatch_tagprobepass_eltop","PNet tvsQCD score of tag jet in e-channel after tag & probe jet pass and tag-gen & probe-gen match cut",20,0.96,1);
  
  hist_pnettagak8_sf[0][1] = new TH1D("hist_pnettagak8_tagpass_mutop","PNet tvsQCD score of tag jet in #mu-channel after tag jet pass cut",20,0.96,1);
  hist_pnettagak8_sf[1][1] = new TH1D("hist_pnettagak8_genmatch_tagpass_mutop","PNet tvsQCD score of tag jet in #mu-channel after tag jet pass and tag-gen match cut",20,0.96,1);
  hist_pnettagak8_sf[2][1] = new TH1D("hist_pnettagak8_probepass_mutop","PNet tvsQCD score of tag jet in #mu-channel after tag & probe jet pass cut",20,0.96,1);
  hist_pnettagak8_sf[3][1] = new TH1D("hist_pnettagak8_genmatch_probepass_mutop","PNet tvsQCD score of tag jet in #mu-channel after tag & probe jet pass and tag-gen match cut",20,0.96,1);
  hist_pnettagak8_sf[4][1] = new TH1D("hist_pnettagak8_genmatch_tagprobepass_mutop","PNet tvsQCD score of tag jet in #mu-channel after tag & probe jet pass and tag-gen & probe-gen match cut",20,0.96,1);

  hist_probescoreak8_sf[0][0] = new TH1D("hist_probescoreak8_tagpass_eltop","Leptop tagger score of el top after tag jet pass cut",30,-1,1);
  hist_probescoreak8_sf[1][0] = new TH1D("hist_probescoreak8_genmatch_tagpass_eltop","Leptop tagger score of el top after tag jet pass and tag-gen match cut",30,-1,1);
  hist_probescoreak8_sf[2][0] = new TH1D("hist_probescoreak8_probepass_eltop","Leptop tagger score of el top after tag & probe jet pass cut",30,-1,1);
  hist_probescoreak8_sf[3][0] = new TH1D("hist_probescoreak8_genmatch_probepass_eltop","Leptop tagger score of el top after tag & probe jet pass and tag-gen match cut",30,-1,1);
  hist_probescoreak8_sf[4][0] = new TH1D("hist_probescoreak8_genmatch_tagprobepass_eltop","Leptop tagger score of el top after tag & probe jet pass and tag-gen & probe-gen match cut",30,-1,1);
  
  hist_probescoreak8_sf[0][1] = new TH1D("hist_probescoreak8_tagpass_mutop","Leptop tagger score of mu top after tag jet pass cut",30,-1,1);
  hist_probescoreak8_sf[1][1] = new TH1D("hist_probescoreak8_genmatch_tagpass_mutop","Leptop tagger score of mu top after tag jet pass and tag-gen match cut",30,-1,1);
  hist_probescoreak8_sf[2][1] = new TH1D("hist_probescoreak8_probepass_mutop","Leptop tagger score of mu top after tag & probe jet pass cut",30,-1,1);
  hist_probescoreak8_sf[3][1] = new TH1D("hist_probescoreak8_genmatch_probepass_mutop","Leptop tagger score of mu top after tag & probe jet pass and tag-gen match cut",30,-1,1);
  hist_probescoreak8_sf[4][1] = new TH1D("hist_probescoreak8_genmatch_tagprobepass_mutop","Leptop tagger score of mu top after tag & probe jet pass and tag-gen & probe-gen match cut",30,-1,1);
  

  for(int ip =0; ip<2; ip++){
    string ipass;
    if(ip == 0) ipass = "_pass";
    else ipass = "_fail";
    hist_ak8pt_mcdata_allmatch[0][ip] = new TH1D(("hist_ak8pt_eltop_allmatch" + ipass).c_str(),"pt of el top jet",25,300,1300);
    hist_ak8eta_mcdata_allmatch[0][ip] = new TH1D(("hist_ak8eta_eltop_allmatch" + ipass).c_str(),"#eta of el top",12,-2.5,2.5);
    hist_ak8sdmass_mcdata_allmatch[0][ip] = new TH1D(("hist_ak8sdmass_eltop_allmatch" + ipass).c_str(),"sdmass of el top",25,0,250);
    hist_ak8score_mcdata_allmatch[0][ip] = new TH1D(("hist_response_eltop_allmatch" + ipass).c_str(),"Electronic top tagger score",20,-1,1);
    
    hist_ak8pt_mcdata_allmatch[1][ip] = new TH1D(("hist_ak8pt_mutop_allmatch" + ipass).c_str(),"pt of mu top jet",25,300,1300);
    hist_ak8eta_mcdata_allmatch[1][ip] = new TH1D(("hist_ak8eta_mutop_allmatch" + ipass).c_str(),"#eta of mu top",12,-2.5,2.5);
    hist_ak8sdmass_mcdata_allmatch[1][ip] = new TH1D(("hist_ak8sdmass_mutop_allmatch" + ipass).c_str(),"sdmass of mu top",25,0,250);
    hist_ak8score_mcdata_allmatch[1][ip] = new TH1D(("hist_response_mutop_allmatch" + ipass).c_str(),"Muonic top tagger score",20,-1,1);
    
    hist_ak8pt_mcdata_leporbmatch[0][ip] = new TH1D(("hist_ak8pt_eltop_leporbmatch" + ipass).c_str(),"pt of el top jet",25,300,1300);
    hist_ak8eta_mcdata_leporbmatch[0][ip] = new TH1D(("hist_ak8eta_eltop_leporbmatch" + ipass).c_str(),"#eta of el top",12,-2.5,2.5);
    hist_ak8sdmass_mcdata_leporbmatch[0][ip] = new TH1D(("hist_ak8sdmass_eltop_leporbmatch" + ipass).c_str(),"sdmass of el top",25,0,250);
    hist_ak8score_mcdata_leporbmatch[0][ip] = new TH1D(("hist_response_eltop_leporbmatch" + ipass).c_str(),"Electronic top tagger score",20,-1,1);
    
    hist_ak8pt_mcdata_leporbmatch[1][ip] = new TH1D(("hist_ak8pt_mutop_leporbmatch" + ipass).c_str(),"pt of mu top jet",25,300,1300);
    hist_ak8eta_mcdata_leporbmatch[1][ip] = new TH1D(("hist_ak8eta_mutop_leporbmatch" + ipass).c_str(),"#eta of mu top",12,-2.5,2.5);
    hist_ak8sdmass_mcdata_leporbmatch[1][ip] = new TH1D(("hist_ak8sdmass_mutop_leporbmatch" + ipass).c_str(),"sdmass of mu top",25,0,250);
    hist_ak8score_mcdata_leporbmatch[1][ip] = new TH1D(("hist_response_mutop_leporbmatch" + ipass).c_str(),"Muonic top tagger score",20,-1,1);
    
    hist_ak8pt_mcdata_nomatch[0][ip] = new TH1D(("hist_ak8pt_eltop_nomatch" + ipass).c_str(),"pt of el top jet",25,300,1300);
    hist_ak8eta_mcdata_nomatch[0][ip] = new TH1D(("hist_ak8eta_eltop_nomatch" + ipass).c_str(),"#eta of el top",12,-2.5,2.5);
    hist_ak8sdmass_mcdata_nomatch[0][ip] = new TH1D(("hist_ak8sdmass_eltop_nomatch" + ipass).c_str(),"sdmass of el top",25,0,250);
    hist_ak8score_mcdata_nomatch[0][ip] = new TH1D(("hist_response_eltop_nomatch" + ipass).c_str(),"Electronic top tagger score",20,-1,1);
    
    hist_ak8pt_mcdata_nomatch[1][ip] = new TH1D(("hist_ak8pt_mutop_nomatch" + ipass).c_str(),"pt of mu top jet",25,300,1300);
    hist_ak8eta_mcdata_nomatch[1][ip] = new TH1D(("hist_ak8eta_mutop_nomatch" + ipass).c_str(),"#eta of mu top",12,-2.5,2.5);
    hist_ak8sdmass_mcdata_nomatch[1][ip] = new TH1D(("hist_ak8sdmass_mutop_nomatch" + ipass).c_str(),"sdmass of mu top",25,0,250);
    hist_ak8score_mcdata_nomatch[1][ip] = new TH1D(("hist_response_mutop_nomatch" + ipass).c_str(),"Muonic top tagger score",20,-1,1);
  }
  
  for(int ip = 0; ip < nptbin; ip++)
    {
      for(int ie = 0; ie < netabin; ie++)
	{
	  hist_ak8sdmass_sf_allmatch_pass[ip][ie][0] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_eltop_allmatch_pass").c_str(),"sdmass of el top",25,0,250);
	  hist_ak8sdmass_sf_allmatch_pass[ip][ie][1] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_mutop_allmatch_pass").c_str(),"sdmass of mu top",25,0,250);
	  
	  hist_ak8sdmass_sf_leporbmatch_pass[ip][ie][0] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_eltop_leporbmatch_pass").c_str(),"sdmass of el top",25,0,250);
	  hist_ak8sdmass_sf_leporbmatch_pass[ip][ie][1] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_mutop_leporbmatch_pass").c_str(),"sdmass of mu top",25,0,250);

	  hist_ak8sdmass_sf_nomatch_pass[ip][ie][0] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_eltop_nomatch_pass").c_str(),"sdmass of el top",25,0,250);
	  hist_ak8sdmass_sf_nomatch_pass[ip][ie][1] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_mutop_nomatch_pass").c_str(),"sdmass of mu top",25,0,250);

	  hist_ak8sdmass_sf_allmatch_fail[ip][ie][0] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_eltop_allmatch_fail").c_str(),"sdmass of el top",25,0,250);
	  hist_ak8sdmass_sf_allmatch_fail[ip][ie][1] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_mutop_allmatch_fail").c_str(),"sdmass of mu top",25,0,250);
	  
	  hist_ak8sdmass_sf_leporbmatch_fail[ip][ie][0] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_eltop_leporbmatch_fail").c_str(),"sdmass of el top",25,0,250);
	  hist_ak8sdmass_sf_leporbmatch_fail[ip][ie][1] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_mutop_leporbmatch_fail").c_str(),"sdmass of mu top",25,0,250);

	  hist_ak8sdmass_sf_nomatch_fail[ip][ie][0] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_eltop_nomatch_fail").c_str(),"sdmass of el top",25,0,250);
	  hist_ak8sdmass_sf_nomatch_fail[ip][ie][1] = new TH1D(("hist_ak8sdmass_bin_" + to_string(ip) + "_" + to_string(ie) + "_mutop_nomatch_fail").c_str(),"sdmass of mu top",25,0,250);
	}
    }

  for(int lvar=0;lvar<7;lvar++)
    {
      char lnamein[1000];
      sprintf(lnamein,"Var_%s",names_elidvar[lvar]);
      hist_elidvar[lvar] = new TH1D(lnamein,lnamein,elidvar_nbins[lvar],elidvar_low[lvar],elidvar_up[lvar]);
      hist_elidvar[lvar]->Sumw2();
      
      sprintf(lnamein,"Var_%s",names_muidvar[lvar]);
      hist_muidvar[lvar] = new TH1D(lnamein,lnamein,muidvar_nbins[lvar],muidvar_low[lvar],muidvar_up[lvar]);
      hist_muidvar[lvar]->Sumw2();
    }
  
  calib_deepflav = BTagCalibration("DeepJet", (dir + "DeepJet_106XUL18SF_WPonly_V1p1.csv").Data());
  reader_deepflav_tight = BTagCalibrationReader(BTagEntry::OP_TIGHT, "central", {"up", "down"}); 
  reader_deepflav_tight.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
  reader_deepflav_tight.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
  reader_deepflav_tight.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");

  reader_deepflav_med = BTagCalibrationReader(BTagEntry::OP_MEDIUM, "central", {"up", "down"}); 
  reader_deepflav_med.load(calib_deepflav, BTagEntry::FLAV_B, "comb");
  reader_deepflav_med.load(calib_deepflav, BTagEntry::FLAV_C, "comb");
  reader_deepflav_med.load(calib_deepflav, BTagEntry::FLAV_UDSG, "incl");
  
  reader1 = new TMVA::Reader( "BDTG_Re" );
  reader1->AddVariable("selpfjetAK8NHadF", &in_pfjetAK8NHadF);
  reader1->AddVariable("selpfjetAK8neunhadfrac", &in_pfjetAK8neunhadfrac);
  reader1->AddVariable("selpfjetAK8subhaddiff", &in_pfjetAK8subhaddiff);
  reader1->AddVariable("selpfjetAK8tau21", &in_pfjetAK8tau21);
  reader1->AddVariable("selpfjetAK8chrad", &in_pfjetAK8chrad);
  reader1->AddVariable("selpfjetAK8sdmass", &in_pfjetAK8sdmass);
  reader1->AddVariable("selpfjetAK8matchedeldxy_sv", &in_pfjetAK8matchedeldxy_sv);
  reader1->AddVariable("selpfjetAK8matchedelcleta", &in_pfjetAK8matchedelcleta);
  reader1->AddVariable("selpfjetAK8matchedelpt", &in_pfjetAK8matchedelpt);
  reader1->AddVariable("selpfjetAK8matchedelsigmaieta", &in_pfjetAK8matchedelsigmaieta);
  reader1->AddVariable("selpfjetAK8matchedelsigmaiphi", &in_pfjetAK8matchedelsigmaiphi);
  reader1->AddVariable("selpfjetAK8matchedelr9full", &in_pfjetAK8matchedelr9full);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_etaw", &in_pfjetAK8matchedelsupcl_etaw);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_phiw", &in_pfjetAK8matchedelsupcl_phiw);
  reader1->AddVariable("selpfjetAK8matchedelhcaloverecal", &in_pfjetAK8matchedelhcaloverecal);
  reader1->AddVariable("selpfjetAK8matchedelcloctftrkn", &in_pfjetAK8matchedelcloctftrkn);
  reader1->AddVariable("selpfjetAK8matchedelcloctftrkchi2", &in_pfjetAK8matchedelcloctftrkchi2);
  reader1->AddVariable("selpfjetAK8matchedele1x5bye5x5", &in_pfjetAK8matchedele1x5bye5x5);
  reader1->AddVariable("selpfjetAK8matchedelnormchi2", &in_pfjetAK8matchedelnormchi2);
  reader1->AddVariable("selpfjetAK8matchedelhitsmiss", &in_pfjetAK8matchedelhitsmiss);
  reader1->AddVariable("selpfjetAK8matchedeltrkmeasure", &in_pfjetAK8matchedeltrkmeasure);
  reader1->AddVariable("selpfjetAK8matchedelecloverpout", &in_pfjetAK8matchedelecloverpout);
  reader1->AddVariable("selpfjetAK8matchedelecaletrkmomentum", &in_pfjetAK8matchedelecaletrkmomentum);
  reader1->AddVariable("selpfjetAK8matchedeldeltaetacltrkcalo", &in_pfjetAK8matchedeldeltaetacltrkcalo);
  reader1->AddVariable("selpfjetAK8matchedelsupcl_preshvsrawe", &in_pfjetAK8matchedelsupcl_preshvsrawe);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumphet", &in_pfjetAK8matchedelpfisolsumphet);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumchhadpt", &in_pfjetAK8matchedelpfisolsumchhadpt);
  reader1->AddVariable("selpfjetAK8matchedelpfisolsumneuhadet", &in_pfjetAK8matchedelpfisolsumneuhadet);
  reader1->AddVariable("selpfjetAK8matchedeletain", &in_pfjetAK8matchedeletain);
  reader1->AddVariable("selpfjetAK8matchedelphiin", &in_pfjetAK8matchedelphiin);
  reader1->AddVariable("selpfjetAK8matchedelfbrem", &in_pfjetAK8matchedelfbrem);
  reader1->AddVariable("selpfjetAK8matchedeleoverp", &in_pfjetAK8matchedeleoverp);
  reader1->AddVariable("selpfjetAK8matchedelhovere", &in_pfjetAK8matchedelhovere);
  reader1->AddVariable("selpfjetAK8matchedelRho", &in_pfjetAK8matchedelRho);
  reader1->BookMVA("BDTG method", weightfile1);
    
  reader4 = new TMVA::Reader( "BDTG_Rmu" );
  reader4->AddVariable( "selpfjetAK8NHadF", &in_mupfjetAK8NHadF);
  reader4->AddVariable( "selpfjetAK8neunhadfrac", &in_mupfjetAK8neunhadfrac);
  reader4->AddVariable( "selpfjetAK8subhaddiff", &in_mupfjetAK8subhaddiff);
  reader4->AddVariable( "selpfjetAK8tau21", &in_mupfjetAK8tau21);
  reader4->AddVariable( "selpfjetAK8chrad", &in_mupfjetAK8chrad);
  reader4->AddVariable( "selpfjetAK8sdmass", &in_mupfjetAK8sdmass);

  
  reader4->AddVariable("selpfjetAK8matchedmuonchi", &in_pfjetAK8matchedmuonchi);
  reader4->AddVariable("selpfjetAK8matchedmuonposmatch", &in_pfjetAK8matchedmuonposmatch);
  reader4->AddVariable("selpfjetAK8matchedmuontrkink", &in_pfjetAK8matchedmuontrkink);
  reader4->AddVariable("selpfjetAK8matchedmuonsegcom", &in_pfjetAK8matchedmuonsegcom);
  reader4->AddVariable("selpfjetAK8matchedmuonhit", &in_pfjetAK8matchedmuonhit);
  reader4->AddVariable("selpfjetAK8matchedmuonmst", &in_pfjetAK8matchedmuonmst);
  reader4->AddVariable("selpfjetAK8matchedmuontrkvtx", &in_pfjetAK8matchedmuontrkvtx);
  reader4->AddVariable("selpfjetAK8matchedmuondz", &in_pfjetAK8matchedmuondz);
  reader4->AddVariable("selpfjetAK8matchedmuonpixhit", &in_pfjetAK8matchedmuonpixhit);
  reader4->AddVariable("selpfjetAK8matchedmuontrklay", &in_pfjetAK8matchedmuontrklay);
  reader4->AddVariable("selpfjetAK8matchedmuonvalfrac", &in_pfjetAK8matchedmuonvalfrac);
  reader4->AddVariable("selpfjetAK8muinsubptrat", &in_pfjetAK8muinsubptrat);
  reader4->AddVariable("selpfjetAK8muinsubmassrat", &in_pfjetAK8muinsubmassrat);
  reader4->AddVariable("selpfjetAK8muinsubinvmass", &in_pfjetAK8muinsubinvmass);
  reader4->AddVariable("selpfjetAK8muinsubIfarbyI0", &in_pfjetAK8muinsubIfarbyI0);
  reader4->AddVariable("selpfjetAK8muinsubInearbyI0", &in_pfjetAK8muinsubInearbyI0);
  reader4->BookMVA("BDTG method", weightfile4);
}

Bool_t top_tagger_sf::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // It can be passed to either top_tagger_sf::GetEntry() or TBranch::GetEntry()
  // to read either all or the required parts of the data. When processing
  // keyed objects with PROOF, the object is already loaded and is available
  // via the fObject pointer.
  //
  // This function should contain the "body" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  int icheck = GetEntry(entry);
  if(isMC){
    weight = event_weight;
    gen_weight = event_weight;
  }else{
    weight = 1;
    gen_weight = 1;
  }
  //Tout->Fill();
  TString str;
  //str = TString::Format("check start for evt %d ",ievt);          
  //if(gProofServ) gProofServ->SendAsynMessage(str);

  float weight_ch[nchannel] = {1.0,1.0};
  weight_ch[0] = weight;
  weight_ch[1] = weight;
    
  vector<GenParton> genpartons;
  vector<GenParton> LHEtops;
  vector<TopQuark> gentops;
  vector<LHEparticle> lheparticles;
  int cntt = 0; int nel=0; int nmu=0;
  tau_decay_mode = 0;
  toppt_wt = 1;
    
  if( isMC ){
    getPartons(genpartons);
    getLHEParticles(lheparticles);
    // Get GEN-level top quarks
    tt_decay_mode = -1;
    if(isTT || isST || isTTX) {
      getGENTops(gentops,genpartons); // after shower (get top quarks from its daughters) --> will tell details about the signature of ttbar events at GEN level
      getLHETops(LHEtops,genpartons); // before shower (get original top quarks which have decayed) --> will be usedto derive top pt reweighting
      tt_decay_mode = 0;
      for(int ip=0; ip <(int)lheparticles.size();ip++)
	{
	  if(abs(lheparticles[ip].pdgId) == 11)
	    { tt_decay_mode +=10; cntt++; nel++;}
	  else if(abs(lheparticles[ip].pdgId) == 13)
	    { tt_decay_mode +=100; cntt++; nmu++;}
	  else if(abs(lheparticles[ip].pdgId) == 15)
	    { tt_decay_mode +=1000; cntt++;}
	}
      tt_decay_mode += 2 - cntt;
      
      // top pt reweighting //
      if(LHEtops.size()==2 && isTT){
	toppt_wt = SF_TOP(0.0615,0.0005,TMath::Min(float(500),float(LHEtops[0].pt)),TMath::Min(float(500),float(LHEtops[1].pt)));
	weight_ch[0] *= toppt_wt;
	weight_ch[1] *= toppt_wt;
      }
    }
  }    

  // Pile up rewighing //
  if(isMC && npu_vert_true>=0 && npu_vert_true<100){
    float *puweights = Get_PU_Weights(npu_vert_true);
    puWeight = puweights[0];
    puWeightup = puweights[1];
    puWeightdown = puweights[2];
    weight_ch[0] *= puWeight;
    weight_ch[0] *= prefiringweight;
    weight_ch[1] *= puWeight;
    weight_ch[1] *= prefiringweight;
  }
  else if(isMC)
    {
      weight = 0;
      str = TString::Format("npu_vert_true = %d weight = %f nprimi = %d is out of range for evt %d  ;",npu_vert_true,weight,nprimi,ievt);          
      if(gProofServ) gProofServ->SendAsynMessage(str);
      return kFALSE;
    }
  
  
  bool itrig_onlysinglee = (!hlt_Mu37Ele27 && !hlt_Mu27Ele37 && !hlt_Mu37TkMu27 && !hlt_DoubleEle25 && !hlt_Mu50 && hlt_Ele50_PFJet165);
  
  //if(!itrig_onlysinglee) return kFALSE; //event should at least fire a single lepton trigger
  
  vector <Electron> velectrons;
  getelectrons(velectrons,30,absetacut);    // electron_pt_cut & absetacut defined in Proof.h
  
  
  //Here you get muons with your criteria
  vector <Muon> vmuons;
  getmuons(vmuons,30,absetacut);        

  
  //Make lepton collection from electrons & muons (using only common variables)
  vector <Lepton> vleptons;
  getLeptons(vleptons,vmuons,velectrons,30);
    

  //Here you get AK4 jets with your criteria                                     
  vector <AK4Jet> Jets;
  getAK4jets(Jets,AK4jet_pt_cut,absetacut,isMC);     
  
  //Here you get AK8 jets with your criteria 
  vector <AK8Jet> LJets;
  getAK8jets(LJets,300,absetacut,isMC);
  
  if(isMC){
    // assign information from GEN-matching
    AssignGen(LJets,genpartons);
    if(isTT || isST || isTTX) TopAssignment_toJet(LJets,LHEtops,gentops);
  }

  leptonsf_weight = leptonsf_weight_stat = leptonsf_weight_syst = 1;
  if(isMC){
    for(int ie=0; ie<(int)velectrons.size(); ie++) {
      float *sfvalues;

      sfvalues = Electron_SF(velectrons[ie].pt,velectrons[ie].eta,"reco");

      leptonsf_weight *= sfvalues[0];
      leptonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[1]*sfvalues[1] + sfvalues[2]*sfvalues[2]));  // like this for time being
      leptonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4] + sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6]));  // like this for time being
    }
    
    /*    for(int im=0; im<(int)vmuons.size(); im++) {
      float *sfvalues;

      sfvalues = Muon_SF(vmuons[im].pt,vmuons[im].eta,"loose");

      leptonsf_weight *= sfvalues[0];
      leptonsf_weight_stat *= sfvalues[1];
      leptonsf_weight_syst *= sfvalues[2];
      }*/
    weight_ch[0] *= leptonsf_weight;
  }

  leptonsf_weight = leptonsf_weight_stat = leptonsf_weight_syst = 1;
  if(isMC){
    /*for(int ie=0; ie<(int)velectrons.size(); ie++) {
      float *sfvalues;

      sfvalues = Electron_SF(velectrons[ie].pt,velectrons[ie].eta,"mva90noiso");

      leptonsf_weight *= sfvalues[0];
      leptonsf_weight_stat *= (sfvalues[0] + sqrt(sfvalues[1]*sfvalues[1] + sfvalues[2]*sfvalues[2]));  // like this for time being
      leptonsf_weight_syst *= (sfvalues[0] + sqrt(sfvalues[3]*sfvalues[3] + sfvalues[4]*sfvalues[4] + sfvalues[5]*sfvalues[5] + sfvalues[6]*sfvalues[6]));  // like this for time being
      }*/
    
    for(int im=0; im<(int)vmuons.size(); im++) {
      float *sfvalues;

      sfvalues = Muon_SF(vmuons[im].pt,vmuons[im].eta,"reco");

      leptonsf_weight *= sfvalues[0];
      leptonsf_weight_stat *= sfvalues[1];
      leptonsf_weight_syst *= sfvalues[2];
      }
    weight_ch[1] *= leptonsf_weight;
  }
  
  btag_SF = btag_SF_up = btag_SF_dn = 1;
  float p_mc, p_data, p_data_up, p_data_dn;
  p_mc = p_data = p_data_up = p_data_dn = 1;

  float btag_loosecut,btag_medcut,btag_tightcut;
  btag_loosecut = 0.0490;   btag_medcut = 0.2783;   btag_tightcut = 0.7100; //for UL18

  if(isMC){
    for(auto & jet: Jets){
    
      BTagEntry::JetFlavor btv_flav;
      if(abs(jet.hadronFlavour)==5){ btv_flav = BTagEntry::FLAV_B; }
      else if (abs(jet.hadronFlavour)==4){ btv_flav = BTagEntry::FLAV_C; }
      else { btv_flav = BTagEntry::FLAV_UDSG; }
      
      jet.btag_DeepFlav_med_SF = reader_deepflav_med.eval_auto_bounds("central",btv_flav,fabs(jet.eta),jet.pt); 
      jet.btag_DeepFlav_med_SF_up = reader_deepflav_med.eval_auto_bounds("up",btv_flav,fabs(jet.eta),jet.pt);
      jet.btag_DeepFlav_med_SF_dn = reader_deepflav_med.eval_auto_bounds("down",btv_flav,fabs(jet.eta),jet.pt);
      
      //jet.btag_DeepFlav_loose_SF = reader_deepflav_loose.eval_auto_bounds("central",btv_flav,fabs(jet.eta),jet.pt); 
      //jet.btag_DeepFlav_loose_SF_up = reader_deepflav_loose.eval_auto_bounds("up",btv_flav,fabs(jet.eta),jet.pt);
      //jet.btag_DeepFlav_loose_SF_dn = reader_deepflav_loose.eval_auto_bounds("down",btv_flav,fabs(jet.eta),jet.pt);

      double eff_l, eff_m;
      string sample;
      if(isTT) sample = "tt";
      else if(isQCD) sample = "qcd";
      else if(isWJ) sample = "wj";
      else if(isST) sample = "st";
      else if(isTTX) sample = "ttx";
      else if(isDY) sample = "dy";
      else if(isDIB) sample = "dib";
      else sample = "tt";
	
      eff_m = BTag_MCEfficiency_M(sample,abs(jet.hadronFlavour),jet.pt,fabs(jet.eta));
      if(jet.btag_DeepFlav > btag_medcut)
	{
	  p_mc *= max((double)0.00001,eff_m);
	  p_data *= jet.btag_DeepFlav_med_SF*eff_m;
	  p_data_up *= jet.btag_DeepFlav_med_SF_up*eff_m;
	  p_data_dn *= jet.btag_DeepFlav_med_SF_dn*eff_m;
	}
      else
	{
	  p_mc *= max((double)0.00001,1 - eff_m);
	  p_data *= 1 - jet.btag_DeepFlav_med_SF*eff_m;
	  p_data_up *= 1 - jet.btag_DeepFlav_med_SF_up*eff_m;
	  p_data_dn *= 1 - jet.btag_DeepFlav_med_SF_dn*eff_m;
	}
    }
  }

  btag_SF = p_data/p_mc;
  btag_SF_up = p_data_up/p_mc;
  btag_SF_dn = p_data_dn/p_mc;
  weight_ch[0] *= btag_SF;
  weight_ch[1] *= btag_SF;
  
  dirgltrthr = 0;
  dirglthrmin = 0;
  
  if ((int)Jets.size()>=1) {
    std::vector<TLorentzVector> allsjets_4v;
    double Pt_sum(0.);
    for(int ijet=0; ijet< (int)Jets.size(); ijet++){
      TLorentzVector sjv = Jets[ijet].p4;
      allsjets_4v.push_back(sjv);
	Pt_sum += Jets[ijet].pt;
    }
      
    for(int ilep=0; ilep<(int)velectrons.size(); ilep++){
      TLorentzVector sjv = velectrons[ilep].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += velectrons[ilep].p4.Pt();
    }
      
    for(int ilep=0; ilep<(int)vmuons.size(); ilep++){
      TLorentzVector sjv = vmuons[ilep].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += vmuons[ilep].p4.Pt();
    }
    std::vector<double> ThrustAxis;
    std::vector<double> Thrust;
    
    for(unsigned int j =0;j<4;j++){                                                            
      Thrust.push_back(0.01);                                                                    
    }
    Thrust_calculate(allsjets_4v,Thrust);
    dirgltrthr = Thrust[3];
    
    double alpha=atan2(Thrust[1],Thrust[0]);
    for( int i=0; i<(int)allsjets_4v.size(); i++){
      dirglthrmin += fabs(-sin(alpha)*allsjets_4v[i].Px()+cos(alpha)*allsjets_4v[i].Py());
    }
    dirglthrmin = dirglthrmin/max(1.0, Pt_sum);
  }
  
  EventHT = -1;	
  rat_metpt_eventHT = -1;
  
  if ((int)Jets.size()>=1) {
    std::vector<TLorentzVector> allsjets_4v;
    double Pt_sum(0.);
    for(int ijet=0; ijet< (int)Jets.size(); ijet++){
      TLorentzVector sjv = Jets[ijet].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += Jets[ijet].pt;
    }
    
    for(int ilep=0; ilep<(int)velectrons.size(); ilep++){
      TLorentzVector sjv = velectrons[ilep].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += velectrons[ilep].p4.Pt();
    }
    
    for(int ilep=0; ilep<(int)vmuons.size(); ilep++){
      TLorentzVector sjv = vmuons[ilep].p4;
      allsjets_4v.push_back(sjv);
      Pt_sum += vmuons[ilep].p4.Pt();
    }
    
    EventHT = Pt_sum;	
    rat_metpt_eventHT = chs_met/max(float(1.0), EventHT);
  }
  
  bool event_pass[nchannel];
  int ihadjet = -1, ilepjet = -1;

  if((int)LJets.size()>1) {
    if(LJets[0].jetid_tightlepveto && !LJets[1].jetid_tightlepveto)
      {
	ihadjet = 0;
	ilepjet = 1;
      }
    else if(!LJets[0].jetid_tightlepveto && LJets[1].jetid_tightlepveto)
      {
	ihadjet = 1;
	ilepjet = 0;
      }
    else if(LJets[0].DeepTag_TvsQCD > LJets[1].DeepTag_TvsQCD)
      {
	ihadjet = 0;
	ilepjet = 1;
      }
    else if(LJets[0].DeepTag_TvsQCD < LJets[1].DeepTag_TvsQCD)
      {
	ihadjet = 1;
	ilepjet = 0;
      }
    else if(LJets[0].DeepTag_TvsQCD == LJets[1].DeepTag_TvsQCD)
      {
        ihadjet = 0;
        ilepjet = 1;
      }
  }

  vector <AK4Jet> BJets;
  bool bjet_inhadjet = false;
  for(auto & jet: Jets){
    if(isBJet(jet,btag_medcut)){
      BJets.push_back(jet);
      if((int)LJets.size()>1 && delta2R(jet.p4,LJets[ihadjet].p4)<0.8)
	bjet_inhadjet = true;
    }
  }
  
  nbjets = (int)BJets.size();

  nel_veto=0;
  nmu_veto=0;

  for(int ie=0; ie<(int)velectrons.size(); ie++) {
    if(velectrons[ie].mvaWPloose_noIso)
      nel_veto++;
  }
  for(int im=0; im<(int)vmuons.size(); im++) {
    if(vmuons[im].isLoose)
      nmu_veto++;
  }

  //Match lepton with AK8 jets
  int ileppdg = -1;
  Lepton lepcand;
  if((int)LJets.size()>1 && (int)vleptons.size()>0)
    {
      LJets[ilepjet].match_lepton_index = get_nearest_lepton(vleptons,LJets[ilepjet].p4);
      if(LJets[ilepjet].match_lepton_index >= 0)
        {
          ReadTagger(LJets[ilepjet],vleptons,vmuons,velectrons,reader1,reader4);
          if(LJets[ilepjet].match_lepton_index >=0 && LJets[ilepjet].match_lepton_index < (int)vleptons.size())
            {
              lepcand = vleptons[LJets[ilepjet].match_lepton_index];
              ileppdg = abs(lepcand.pdgId);
            }
        }
    }

  bool ttbar_mujet = true;
  bool ttbar_eljet = true;
  if(isMC && isTTsignal)
    {
      if(tt_decay_mode != 11)
	ttbar_eljet = false;
      if(tt_decay_mode != 101)
	ttbar_mujet = false;
    }

  deltaphi = -1.0;
  if((int)LJets.size()>1)  deltaphi = fabs(PhiInRange(LJets[0].phi - LJets[1].phi));
  //if(deltaphi>TMath::Pi()) deltaphi = 2*TMath::Pi() - deltaphi;
  
  event_pass[0] = (int)velectrons.size()>0 && nmu_veto==0 && chs_met > 50 && chs_met <400 && nbjets>0 && (int)LJets.size()>1 && deltaphi > 2.5 && LJets[ihadjet].pt>300 && LJets[ihadjet].jetid_tightlepveto && bjet_inhadjet && ileppdg == 11 && LJets[ihadjet].sdmass > 100 && LJets[ihadjet].sdmass < 240 && EventHT>500 && EventHT <2000 &&  dirgltrthr < 0.35 && dirglthrmin < 0.15;

  event_pass[1] = (int)vmuons.size()>0 && nel_veto==0 && chs_met > 50 && chs_met < 400 && nbjets>0 && (int)LJets.size()>1 && deltaphi > 2.5 && LJets[ihadjet].pt>300 && LJets[ihadjet].jetid_tightlepveto && bjet_inhadjet && ileppdg == 13 && LJets[ihadjet].sdmass > 100 && LJets[ihadjet].sdmass < 240 && EventHT>500 && EventHT <2000 &&  dirgltrthr < 0.35 && dirglthrmin < 0.15;

  float wp_had_tagger = 0.97;
  bool tag_pass = (int)LJets.size()>1 && LJets[ihadjet].PNet_TvsQCD > wp_had_tagger;
  
  genmatch_tag = false;
  hadtaggersf = 1; float hadtaggersf_up = 1, hadtaggersf_down = 1;
  if(isMC && (isTT || isST || isTTX) && (int)LJets.size()>1)
    {
      genmatch_tag = LJets[ihadjet].hashadtop_alldecay;
      if(genmatch_tag)
	{
	  hadtaggersf = PNetHadTopTag_SF(LJets[ihadjet].pt);
	  hadtaggersf_up = PNetHadTopTag_SF_UP(LJets[ihadjet].pt);
	  hadtaggersf_down = PNetHadTopTag_SF_DOWN(LJets[ihadjet].pt);
	  weight_ch[0] *= hadtaggersf;
	  weight_ch[1] *= hadtaggersf;
	}
    }

  bool probe_pass[nchannel] = {false};
  probe_pass[0] = (int)LJets.size()>1 && LJets[ilepjet].re_tvsb > wp_leptop_tagger[0];
  probe_pass[1] = (int)LJets.size()>1 && LJets[ilepjet].re_tvsb > wp_leptop_tagger[1];

  genmatch_probe = false;
  if(isMC && (isTT || isST || isTTX) && (int)LJets.size()>1)
    {
      genmatch_probe = LJets[ilepjet].hasleptop_alldecay;
    }
  
  for(int ich = 0; ich < nchannel; ich++)
    {
      if(event_pass[ich] && tag_pass)
	{
	  hist2d_toptag_sf[0][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),fabs(LJets[ilepjet].eta),weight_ch[ich]);
	  hist_ptak8_sf[0][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	  hist_etaak8_sf[0][ich]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	  hist_sdmassak8_sf[0][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	  hist_sdmasstagak8_sf[0][ich]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[ich]);
	  hist_pnettagak8_sf[0][ich]->Fill(LJets[ihadjet].PNet_TvsQCD,weight_ch[ich]);
	  hist_probescoreak8_sf[0][ich]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	}
    }
  
  for(int ich = 0; ich < nchannel; ich++)
    {
      if(event_pass[ich] && tag_pass && genmatch_tag)
	{
	  hist2d_toptag_sf[1][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),fabs(LJets[ilepjet].eta),weight_ch[ich]);
	  hist_ptak8_sf[1][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	  hist_etaak8_sf[1][ich]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	  hist_sdmassak8_sf[1][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	  hist_sdmasstagak8_sf[1][ich]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[ich]);
	  hist_pnettagak8_sf[1][ich]->Fill(LJets[ihadjet].PNet_TvsQCD,weight_ch[ich]);
	  hist_probescoreak8_sf[1][ich]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	}
    }
  
  for(int ich = 0; ich < nchannel; ich++)
    {
      if(event_pass[ich] && tag_pass && probe_pass[ich])
	{
	  hist2d_toptag_sf[2][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),fabs(LJets[ilepjet].eta),weight_ch[ich]);
	  hist_ptak8_sf[2][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	  hist_etaak8_sf[2][ich]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	  hist_sdmassak8_sf[2][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	  hist_sdmasstagak8_sf[2][ich]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[ich]);
	  hist_pnettagak8_sf[2][ich]->Fill(LJets[ihadjet].PNet_TvsQCD,weight_ch[ich]);
	  hist_probescoreak8_sf[2][ich]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	}
    }
  
  for(int ich = 0; ich < nchannel; ich++)
    {
      if(event_pass[ich] && tag_pass && probe_pass[ich] && genmatch_tag)
	{
	  hist2d_toptag_sf[3][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),fabs(LJets[ilepjet].eta),weight_ch[ich]);
	  hist_ptak8_sf[3][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	  hist_etaak8_sf[3][ich]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	  hist_sdmassak8_sf[3][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	  hist_sdmasstagak8_sf[3][ich]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[ich]);
	  hist_pnettagak8_sf[3][ich]->Fill(LJets[ihadjet].PNet_TvsQCD,weight_ch[ich]);
	  hist_probescoreak8_sf[3][ich]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	}
    }
  
  for(int ich = 0; ich < nchannel; ich++)
    {
      if(event_pass[ich] && tag_pass && probe_pass[ich] && genmatch_tag && genmatch_probe)
	{
	  hist2d_toptag_sf[4][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),fabs(LJets[ilepjet].eta),weight_ch[ich]);
	  hist_ptak8_sf[4][ich]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	  hist_etaak8_sf[4][ich]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	  hist_sdmassak8_sf[4][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	  hist_sdmasstagak8_sf[4][ich]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[ich]);
	  hist_pnettagak8_sf[4][ich]->Fill(LJets[ihadjet].PNet_TvsQCD,weight_ch[ich]);
	  hist_probescoreak8_sf[4][ich]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	}
    }
  
  topjetmatch = true; bjetmatch = true;
  bool lepmatch[2] = {true,true};
  
  if(isMC && (int)LJets.size()>1)
    {
      topjetmatch = genmatch_probe;
      bjetmatch = LJets[ilepjet].hasb;
      lepmatch[0] = LJets[ilepjet].haselectron;
      lepmatch[1] = LJets[ilepjet].hasmuon;
    }
    for(int ich = 0; ich < nchannel; ich++)
    {
      if(event_pass[ich] && tag_pass && probe_pass[ich])
	{
	  int eta_bin_id = (hist2d_toptag_sf[2][ich]->GetYaxis()->FindBin(fabs(LJets[ilepjet].eta))>(hist2d_toptag_sf[2][ich]->GetNbinsY())) ? hist2d_toptag_sf[2][ich]->GetNbinsY() : hist2d_toptag_sf[2][ich]->GetYaxis()->FindBin(fabs(LJets[ilepjet].eta));
	  int pt_bin_id = (hist2d_toptag_sf[2][ich]->GetXaxis()->FindBin(LJets[ilepjet].pt)>(hist2d_toptag_sf[2][ich]->GetNbinsX())) ? hist2d_toptag_sf[2][ich]->GetNbinsX() : hist2d_toptag_sf[2][ich]->GetXaxis()->FindBin(LJets[ilepjet].pt);
	  
	  if(topjetmatch)
	    {
	      hist_ak8pt_mcdata_allmatch[ich][0]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	      hist_ak8eta_mcdata_allmatch[ich][0]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	      hist_ak8sdmass_mcdata_allmatch[ich][0]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	      hist_ak8score_mcdata_allmatch[ich][0]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	      hist_ak8sdmass_sf_allmatch_pass[max(int(0),int(pt_bin_id -1))][max(int(0),int(eta_bin_id -1))][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	    }
	  else if(bjetmatch || lepmatch[ich])
	    {
	      hist_ak8pt_mcdata_leporbmatch[ich][0]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	      hist_ak8eta_mcdata_leporbmatch[ich][0]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	      hist_ak8sdmass_mcdata_leporbmatch[ich][0]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	      hist_ak8score_mcdata_leporbmatch[ich][0]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	      hist_ak8sdmass_sf_leporbmatch_pass[max(int(0),int(pt_bin_id -1))][max(int(0),int(eta_bin_id -1))][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	    }
	  else
	    {
	      hist_ak8pt_mcdata_nomatch[ich][0]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	      hist_ak8eta_mcdata_nomatch[ich][0]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	      hist_ak8sdmass_mcdata_nomatch[ich][0]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	      hist_ak8score_mcdata_nomatch[ich][0]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	      hist_ak8sdmass_sf_nomatch_pass[max(int(0),int(pt_bin_id -1))][max(int(0),int(eta_bin_id -1))][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	    }	  
	}
      
      if(event_pass[ich] && tag_pass && !probe_pass[ich])
	{
	  int eta_bin_id = (hist2d_toptag_sf[2][ich]->GetYaxis()->FindBin(fabs(LJets[ilepjet].eta))>(hist2d_toptag_sf[2][ich]->GetNbinsY())) ? hist2d_toptag_sf[2][ich]->GetNbinsY() : hist2d_toptag_sf[2][ich]->GetYaxis()->FindBin(fabs(LJets[ilepjet].eta));
	  int pt_bin_id = (hist2d_toptag_sf[2][ich]->GetXaxis()->FindBin(LJets[ilepjet].pt)>(hist2d_toptag_sf[2][ich]->GetNbinsX())) ? hist2d_toptag_sf[2][ich]->GetNbinsX() : hist2d_toptag_sf[2][ich]->GetXaxis()->FindBin(LJets[ilepjet].pt);
 
	  if(topjetmatch)
	    {
	      hist_ak8pt_mcdata_allmatch[ich][1]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	      hist_ak8eta_mcdata_allmatch[ich][1]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	      hist_ak8sdmass_mcdata_allmatch[ich][1]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	      hist_ak8score_mcdata_allmatch[ich][1]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	      hist_ak8sdmass_sf_allmatch_fail[max(int(0),int(pt_bin_id -1))][max(int(0),int(eta_bin_id -1))][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	    }
	  else if(bjetmatch || lepmatch[ich])
	    {
	      hist_ak8pt_mcdata_leporbmatch[ich][1]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	      hist_ak8eta_mcdata_leporbmatch[ich][1]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	      hist_ak8sdmass_mcdata_leporbmatch[ich][1]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	      hist_ak8score_mcdata_leporbmatch[ich][1]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	      hist_ak8sdmass_sf_leporbmatch_fail[max(int(0),int(pt_bin_id -1))][max(int(0),int(eta_bin_id -1))][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);	      
	    }
	  else
	    {
	      hist_ak8pt_mcdata_nomatch[ich][1]->Fill(min(float(1299.0),LJets[ilepjet].pt),weight_ch[ich]);
	      hist_ak8eta_mcdata_nomatch[ich][1]->Fill(LJets[ilepjet].eta,weight_ch[ich]);
	      hist_ak8sdmass_mcdata_nomatch[ich][1]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	      hist_ak8score_mcdata_nomatch[ich][1]->Fill(LJets[ilepjet].re_tvsb,weight_ch[ich]);
	      hist_ak8sdmass_sf_nomatch_fail[max(int(0),int(pt_bin_id -1))][max(int(0),int(eta_bin_id -1))][ich]->Fill(min((float)249.0,LJets[ilepjet].sdmass),weight_ch[ich]);
	    }	  
	}
    }

    if((int)nelecs>0 && nmu_veto==0 && nbjets>0 && (int)LJets.size()>1 && deltaphi > 2.5 && LJets[ihadjet].pt>300 && LJets[ihadjet].jetid_tightlepveto && bjet_inhadjet && ileppdg == 11 && LJets[ihadjet].sdmass > 30 && tag_pass){
      hist_elidvar[0]->Fill(min((float)899.0,chs_met),weight_ch[0]);
      hist_elidvar[1]->Fill(min((float)0.3,dirgltrthr),weight_ch[0]);
      hist_elidvar[2]->Fill(min((float)0.6,dirglthrmin),weight_ch[0]);
      hist_elidvar[3]->Fill(deltaphi,weight_ch[0]);
      hist_elidvar[4]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[0]);
      hist_elidvar[5]->Fill(min((float)2499,EventHT),weight_ch[0]);
      hist_elidvar[6]->Fill((int)Jets.size(),weight_ch[0]); 
    }

    if((int)nmuons>0 && nel_veto==0 && chs_met > 50 && nbjets>0 && (int)LJets.size()>1 && deltaphi > 2 && LJets[ihadjet].pt>300 && LJets[ihadjet].jetid_tightlepveto && bjet_inhadjet && ileppdg == 13 && LJets[ihadjet].sdmass > 30 && tag_pass){
      hist_muidvar[0]->Fill(min((float)1499.0,chs_met),weight_ch[1]);
      hist_muidvar[1]->Fill(min((float)0.3,dirgltrthr),weight_ch[1]);
      hist_muidvar[2]->Fill(min((float)0.6,dirglthrmin),weight_ch[1]);
      hist_muidvar[3]->Fill(deltaphi,weight_ch[1]);
      hist_muidvar[4]->Fill(min((float)249.0,LJets[ihadjet].sdmass),weight_ch[1]);
      hist_muidvar[5]->Fill(min((float)2499,EventHT),weight_ch[1]);
      hist_muidvar[6]->Fill((int)Jets.size(),weight_ch[1]); 
    }
       
    if(((int)vmuons.size()>0 || (int)velectrons.size()>0) && (int)LJets.size()>1 && LJets[ihadjet].jetid_tightlepveto && nbjets>0 && bjet_inhadjet && LJets[ihadjet].sdmass > 30 && LJets[ihadjet].PNet_TvsQCD > 0.5)
      {
	nmuons_mutop = (int)vmuons.size();
	nelectrons_eltop = (int)velectrons.size();
	weight_eltop = weight_ch[0];
	weight_mutop = weight_ch[1];
	tagjetpt = LJets[ihadjet].pt;
	tagjet_taggerscore = LJets[ihadjet].DeepTag_TvsQCD;
	probejet_score = LJets[ilepjet].re_tvsb;
	elmatch = lepmatch[0];
	mumatch = lepmatch[1];
	Tnewvar->Fill();
      }
    											
  return kTRUE;     
}
void top_tagger_sf::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  fileOut->cd();
  fileOut->Write();
  fOutput->Add(OutFile);
  fileOut->Close();
}

void top_tagger_sf::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  
}
