/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 28 22:40:34 2019 by ROOT version 5.34/38
// from TTree T1/WPrimeNtuple
// found on file: root://se01.indiacms.res.in//store/user/chatterj/TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8/crab_crab_L5JERC_TTBar_Mtt_700to1000_Autumn18_JECV19/191113_143123/0000/rootuple_jerc_l5_106.root
//////////////////////////////////////////////////////////

#ifndef top_tagger_sf_h
#define top_tagger_sf_h

#include "Objects.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TMath.h>

#include "TLorentzVector.h"
#include <TProofOutputFile.h>
#include <TProofServ.h>

#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

#include "BTagCalibrationStandalone.h"
#include "BTagCalibrationStandalone.cpp"

using namespace std;
using namespace TMVA;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

//Need to ask Suman why some functions are defined in Functions.h and some here in the new structure

double diff_func(double f1, double f2){ 
  double ff = pow(f1-f2,2)*1./pow(f1+f2,2);
  return ff;
}

bool Muon_TightID(bool muonisGL,bool muonisPF, float muonchi, float muonhit, float muonmst,
		  float muontrkvtx, float muondz, float muonpixhit, float muontrklay){
  bool tightid = false;
  if(muonisGL && muonisPF){
    if(muonchi<10 && muonhit>0 && muonmst>1){
      if(fabs(muontrkvtx)<0.2 && fabs(muondz)<0.5){
	if(muonpixhit>0 && muontrklay>5){
	  tightid = true;
	}
      } 
    }
  }	
  return tightid;
}
//https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon

bool Muon_Iso_ID(float muonpfiso)
{
  bool isoid = false;	
  if(muonpfiso<0.15) { isoid = true; } //SR
  //if(muonpfiso>0.15) { isoid = true; }  // CR
  return isoid;
}

double BTag_SF(int flavor, string syst, double pt){
  
  // scale factors taken from the csv file in : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X (for medium WP)
  
  double x = min(pt,999.9);  /// For btagging efficiency need to be updated (not using at the moment)
  
  if(syst ==  "noSyst") {
    if(abs(flavor)==5||abs(flavor)==4){
      if (x >= 20  && x < 1000) return (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))));
    }else{
      if (x >= 20  && x < 1000) return (1.59373+-0.00113028*x+8.66631e-07*x*x+-1.10505/x) ;
    }
  }
  
  if(syst ==  "down") {
    
    if(abs(flavor)==5){
      if(pt >= 20  && pt < 30) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.19459584355354309) );
      if(pt >= 30  && pt < 50) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.04693598672747612) );
      if(pt >= 50  && pt < 70) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.041476961225271225) );
      if(pt >= 70  && pt < 100) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.037213429808616638) );
      if(pt >= 100  && pt < 140) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.033781636506319046) );
      if(pt >= 140  && pt < 200) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.035268638283014297) );
      if(pt >= 200  && pt < 300) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.043516248464584351) );
      if(pt >= 300  && pt < 600) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.10369165241718292) );
      if(pt >= 600  && pt < 1000) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.29925653338432312) );
    }
    if(abs(flavor)==4){
      if(pt >= 20  && pt < 30) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.064865283668041229) );
      if(pt >= 30  && pt < 50) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.015645328909158707) );
      if(pt >= 50  && pt < 70) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.013825654052197933) );
      if(pt >= 70  && pt < 100) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.012404476292431355) );
      if(pt >= 100  && pt < 140) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.011260545812547207) );
      if(pt >= 140  && pt < 200) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.011756212450563908) );
      if(pt >= 200  && pt < 300) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.01450541615486145) );
      if(pt >= 300  && pt < 600) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.034563884139060974) );
      if(pt >= 600  && pt < 1000) return ( 1.0097+((-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18))))))))-0.099752180278301239) );
    }
    if(abs(flavor)!=5 && abs(flavor)!=4){
      if (pt >= 20  && pt < 1000) return ( (1.59373+-0.00113028*x+8.66631e-07*x*x+-1.10505/x)*(1-(0.142253+0.000227323*x+-2.71704e-07*x*x)) );
    }
  }
  
  if(syst ==  "up") {
    
    if(abs(flavor)==5){
      if(pt >= 20  && pt < 30) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.19459584355354309 );
      if(pt >= 30  && pt < 50) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.04693598672747612 );
      if(pt >= 50  && pt < 70) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.041476961225271225 );
      if(pt >= 70  && pt < 100) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.037213429808616638 );
      if(pt >= 100  && pt < 140) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.033781636506319046 );
      if(pt >= 140  && pt < 200) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.035268638283014297 );
      if(pt >= 200  && pt < 300) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.043516248464584351 );
      if(pt >= 300  && pt < 600) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.10369165241718292 );
      if(pt >= 600  && pt < 1000) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.29925653338432312 );
    }
    if(abs(flavor)==4){
      if(pt >= 20  && pt < 30) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.064865283668041229 );
      if(pt >= 30  && pt < 50) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.015645328909158707 );
      if(pt >= 50  && pt < 70) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.013825654052197933 );
      if(pt >= 70  && pt < 100) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.012404476292431355 );
      if(pt >= 100  && pt < 140) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.011260545812547207 );
      if(pt >= 140  && pt < 200) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.011756212450563908 );
      if(pt >= 200  && pt < 300) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.01450541615486145 );
      if(pt >= 300  && pt < 600) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.034563884139060974 );
      if(pt >= 600  && pt < 1000) return ( (1.0097+(-(2.89663e-06*(log(x+19)*(log(x+18)*(3-(-(110.381*log(x+18)))))))))+0.099752180278301239 );
    }
    if(abs(flavor)!=5 && abs(flavor)!=4){
      if (pt >= 20  && pt < 1000) return ( (1.59373+-0.00113028*x+8.66631e-07*x*x+-1.10505/x)*(1+(0.115686+0.000270835*x+-3.2078e-07*x*x)) );
    }
    
  }
  return 1.0;
}

double BTag_MCEfficiency_L(string sample, int flavor, double pt, double eta){
  if(sample == "dy")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.905567 ;
	  if(pt>=50 && pt<70) return 0.92555 ;
	  if(pt>=70 && pt<100) return 0.926837 ;
	  if(pt>=100 && pt<140) return 0.922133 ;
	  if(pt>=140 && pt<200) return 0.928474 ;
	  if(pt>=200 && pt<300) return 0.928943 ;
	  if(pt>=300 && pt<600) return 0.935012 ;
	  if(pt>=600) return 0.94298 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.899925 ;
	  if(pt>=50 && pt<70) return 0.909001 ;
	  if(pt>=70 && pt<100) return 0.921566 ;
	  if(pt>=100 && pt<140) return 0.915961 ;
	  if(pt>=140 && pt<200) return 0.928816 ;
	  if(pt>=200 && pt<300) return 0.924774 ;
	  if(pt>=300 && pt<600) return 0.922638 ;
	  if(pt>=600) return 0.93735 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.874449 ;
	  if(pt>=50 && pt<70) return 0.881698 ;
	  if(pt>=70 && pt<100) return 0.878283 ;
	  if(pt>=100 && pt<140) return 0.907069 ;
	  if(pt>=140 && pt<200) return 0.884462 ;
	  if(pt>=200 && pt<300) return 0.912729 ;
	  if(pt>=300 && pt<600) return 0.897714 ;
	  if(pt>=600) return 0.928375 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.885357 ;
	  if(pt>=50 && pt<70) return 0.865522 ;
	  if(pt>=70 && pt<100) return 0.889644 ;
	  if(pt>=100 && pt<140) return 0.877388 ;
	  if(pt>=140 && pt<200) return 0.870629 ;
	  if(pt>=200 && pt<300) return 0.890291 ;
	  if(pt>=300 && pt<600) return 0.893782 ;
	  if(pt>=600) return 0.91874 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.489011 ;
	  if(pt>=50 && pt<70) return 0.474361 ;
	  if(pt>=70 && pt<100) return 0.478299 ;
	  if(pt>=100 && pt<140) return 0.460146 ;
	  if(pt>=140 && pt<200) return 0.448821 ;
	  if(pt>=200 && pt<300) return 0.486624 ;
	  if(pt>=300 && pt<600) return 0.55235 ;
	  if(pt>=600) return 0.650161 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.482194 ;
	  if(pt>=50 && pt<70) return 0.451984 ;
	  if(pt>=70 && pt<100) return 0.451207 ;
	  if(pt>=100 && pt<140) return 0.458789 ;
	  if(pt>=140 && pt<200) return 0.456792 ;
	  if(pt>=200 && pt<300) return 0.489193 ;
	  if(pt>=300 && pt<600) return 0.54235 ;
	  if(pt>=600) return 0.647703 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.488855 ;
	  if(pt>=50 && pt<70) return 0.467434 ;
	  if(pt>=70 && pt<100) return 0.45497 ;
	  if(pt>=100 && pt<140) return 0.454278 ;
	  if(pt>=140 && pt<200) return 0.466124 ;
	  if(pt>=200 && pt<300) return 0.481173 ;
	  if(pt>=300 && pt<600) return 0.539877 ;
	  if(pt>=600) return 0.693175 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.503897 ;
	  if(pt>=50 && pt<70) return 0.467462 ;
	  if(pt>=70 && pt<100) return 0.441967 ;
	  if(pt>=100 && pt<140) return 0.424543 ;
	  if(pt>=140 && pt<200) return 0.446207 ;
	  if(pt>=200 && pt<300) return 0.474203 ;
	  if(pt>=300 && pt<600) return 0.539517 ;
	  if(pt>=600) return 0.670134 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.200731 ;
	  if(pt>=50 && pt<70) return 0.13241 ;
	  if(pt>=70 && pt<100) return 0.117055 ;
	  if(pt>=100 && pt<140) return 0.102891 ;
	  if(pt>=140 && pt<200) return 0.0888807 ;
	  if(pt>=200 && pt<300) return 0.0907529 ;
	  if(pt>=300 && pt<600) return 0.144651 ;
	  if(pt>=600) return 0.317814 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.220422 ;
	  if(pt>=50 && pt<70) return 0.145122 ;
	  if(pt>=70 && pt<100) return 0.128694 ;
	  if(pt>=100 && pt<140) return 0.115693 ;
	  if(pt>=140 && pt<200) return 0.0984015 ;
	  if(pt>=200 && pt<300) return 0.0990223 ;
	  if(pt>=300 && pt<600) return 0.159046 ;
	  if(pt>=600) return 0.336824 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.273345 ;
	  if(pt>=50 && pt<70) return 0.177947 ;
	  if(pt>=70 && pt<100) return 0.159941 ;
	  if(pt>=100 && pt<140) return 0.147762 ;
	  if(pt>=140 && pt<200) return 0.133582 ;
	  if(pt>=200 && pt<300) return 0.136509 ;
	  if(pt>=300 && pt<600) return 0.225981 ;
	  if(pt>=600) return 0.465439 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.35171 ;
	  if(pt>=50 && pt<70) return 0.215766 ;
	  if(pt>=70 && pt<100) return 0.188233 ;
	  if(pt>=100 && pt<140) return 0.174276 ;
	  if(pt>=140 && pt<200) return 0.167387 ;
	  if(pt>=200 && pt<300) return 0.182848 ;
	  if(pt>=300 && pt<600) return 0.279638 ;
	  if(pt>=600) return 0.454573 ;
	}
      }
    }
  else if(sample == "qcd")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.902352 ;
	  if(pt>=50 && pt<70) return 0.915307 ;
	  if(pt>=70 && pt<100) return 0.920376 ;
	  if(pt>=100 && pt<140) return 0.921623 ;
	  if(pt>=140 && pt<200) return 0.922782 ;
	  if(pt>=200 && pt<300) return 0.930834 ;
	  if(pt>=300 && pt<600) return 0.935952 ;
	  if(pt>=600) return 0.939543 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.901488 ;
	  if(pt>=50 && pt<70) return 0.907403 ;
	  if(pt>=70 && pt<100) return 0.912996 ;
	  if(pt>=100 && pt<140) return 0.914167 ;
	  if(pt>=140 && pt<200) return 0.916191 ;
	  if(pt>=200 && pt<300) return 0.922939 ;
	  if(pt>=300 && pt<600) return 0.923549 ;
	  if(pt>=600) return 0.931809 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.887915 ;
	  if(pt>=50 && pt<70) return 0.898413 ;
	  if(pt>=70 && pt<100) return 0.897696 ;
	  if(pt>=100 && pt<140) return 0.899377 ;
	  if(pt>=140 && pt<200) return 0.898798 ;
	  if(pt>=200 && pt<300) return 0.903549 ;
	  if(pt>=300 && pt<600) return 0.910122 ;
	  if(pt>=600) return 0.926351 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.8881 ;
	  if(pt>=50 && pt<70) return 0.889402 ;
	  if(pt>=70 && pt<100) return 0.881519 ;
	  if(pt>=100 && pt<140) return 0.880708 ;
	  if(pt>=140 && pt<200) return 0.877609 ;
	  if(pt>=200 && pt<300) return 0.88237 ;
	  if(pt>=300 && pt<600) return 0.875434 ;
	  if(pt>=600) return 0.915556 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.481931 ;
	  if(pt>=50 && pt<70) return 0.466843 ;
	  if(pt>=70 && pt<100) return 0.45846 ;
	  if(pt>=100 && pt<140) return 0.461699 ;
	  if(pt>=140 && pt<200) return 0.476007 ;
	  if(pt>=200 && pt<300) return 0.502076 ;
	  if(pt>=300 && pt<600) return 0.556342 ;
	  if(pt>=600) return 0.62951 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.480895 ;
	  if(pt>=50 && pt<70) return 0.46494 ;
	  if(pt>=70 && pt<100) return 0.458564 ;
	  if(pt>=100 && pt<140) return 0.45751 ;
	  if(pt>=140 && pt<200) return 0.468372 ;
	  if(pt>=200 && pt<300) return 0.491343 ;
	  if(pt>=300 && pt<600) return 0.540857 ;
	  if(pt>=600) return 0.622639 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.479873 ;
	  if(pt>=50 && pt<70) return 0.456818 ;
	  if(pt>=70 && pt<100) return 0.448145 ;
	  if(pt>=100 && pt<140) return 0.448193 ;
	  if(pt>=140 && pt<200) return 0.455587 ;
	  if(pt>=200 && pt<300) return 0.475654 ;
	  if(pt>=300 && pt<600) return 0.533171 ;
	  if(pt>=600) return 0.661931 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.501086 ;
	  if(pt>=50 && pt<70) return 0.461148 ;
	  if(pt>=70 && pt<100) return 0.443428 ;
	  if(pt>=100 && pt<140) return 0.434994 ;
	  if(pt>=140 && pt<200) return 0.441925 ;
	  if(pt>=200 && pt<300) return 0.464569 ;
	  if(pt>=300 && pt<600) return 0.543192 ;
	  if(pt>=600) return 0.685836 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.121162 ;
	  if(pt>=50 && pt<70) return 0.0804156 ;
	  if(pt>=70 && pt<100) return 0.0693233 ;
	  if(pt>=100 && pt<140) return 0.0644113 ;
	  if(pt>=140 && pt<200) return 0.0657508 ;
	  if(pt>=200 && pt<300) return 0.0748266 ;
	  if(pt>=300 && pt<600) return 0.123511 ;
	  if(pt>=600) return 0.273714 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.134711 ;
	  if(pt>=50 && pt<70) return 0.0894969 ;
	  if(pt>=70 && pt<100) return 0.0769969 ;
	  if(pt>=100 && pt<140) return 0.0727211 ;
	  if(pt>=140 && pt<200) return 0.072739 ;
	  if(pt>=200 && pt<300) return 0.0797056 ;
	  if(pt>=300 && pt<600) return 0.128969 ;
	  if(pt>=600) return 0.297387 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.17475 ;
	  if(pt>=50 && pt<70) return 0.115469 ;
	  if(pt>=70 && pt<100) return 0.100261 ;
	  if(pt>=100 && pt<140) return 0.0954146 ;
	  if(pt>=140 && pt<200) return 0.0978817 ;
	  if(pt>=200 && pt<300) return 0.109037 ;
	  if(pt>=300 && pt<600) return 0.189974 ;
	  if(pt>=600) return 0.430442 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.25675 ;
	  if(pt>=50 && pt<70) return 0.162866 ;
	  if(pt>=70 && pt<100) return 0.13761 ;
	  if(pt>=100 && pt<140) return 0.132394 ;
	  if(pt>=140 && pt<200) return 0.139593 ;
	  if(pt>=200 && pt<300) return 0.161235 ;
	  if(pt>=300 && pt<600) return 0.262349 ;
	  if(pt>=600) return 0.493307 ;
	}
      }
    }
  else if(sample == "wj")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.899994 ;
	  if(pt>=50 && pt<70) return 0.905468 ;
	  if(pt>=70 && pt<100) return 0.904702 ;
	  if(pt>=100 && pt<140) return 0.909312 ;
	  if(pt>=140 && pt<200) return 0.905616 ;
	  if(pt>=200 && pt<300) return 0.911863 ;
	  if(pt>=300 && pt<600) return 0.922862 ;
	  if(pt>=600) return 0.93567 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.896152 ;
	  if(pt>=50 && pt<70) return 0.904593 ;
	  if(pt>=70 && pt<100) return 0.903827 ;
	  if(pt>=100 && pt<140) return 0.900541 ;
	  if(pt>=140 && pt<200) return 0.897073 ;
	  if(pt>=200 && pt<300) return 0.903928 ;
	  if(pt>=300 && pt<600) return 0.913669 ;
	  if(pt>=600) return 0.93217 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.881679 ;
	  if(pt>=50 && pt<70) return 0.891829 ;
	  if(pt>=70 && pt<100) return 0.884637 ;
	  if(pt>=100 && pt<140) return 0.882495 ;
	  if(pt>=140 && pt<200) return 0.87679 ;
	  if(pt>=200 && pt<300) return 0.882175 ;
	  if(pt>=300 && pt<600) return 0.894753 ;
	  if(pt>=600) return 0.923219 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.881551 ;
	  if(pt>=50 && pt<70) return 0.877085 ;
	  if(pt>=70 && pt<100) return 0.869977 ;
	  if(pt>=100 && pt<140) return 0.860322 ;
	  if(pt>=140 && pt<200) return 0.853857 ;
	  if(pt>=200 && pt<300) return 0.863095 ;
	  if(pt>=300 && pt<600) return 0.879293 ;
	  if(pt>=600) return 0.900473 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.489226 ;
	  if(pt>=50 && pt<70) return 0.474178 ;
	  if(pt>=70 && pt<100) return 0.46603 ;
	  if(pt>=100 && pt<140) return 0.461732 ;
	  if(pt>=140 && pt<200) return 0.471758 ;
	  if(pt>=200 && pt<300) return 0.504238 ;
	  if(pt>=300 && pt<600) return 0.564148 ;
	  if(pt>=600) return 0.643589 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.485791 ;
	  if(pt>=50 && pt<70) return 0.4733 ;
	  if(pt>=70 && pt<100) return 0.465328 ;
	  if(pt>=100 && pt<140) return 0.460108 ;
	  if(pt>=140 && pt<200) return 0.470487 ;
	  if(pt>=200 && pt<300) return 0.494763 ;
	  if(pt>=300 && pt<600) return 0.551211 ;
	  if(pt>=600) return 0.636294 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.48555 ;
	  if(pt>=50 && pt<70) return 0.462082 ;
	  if(pt>=70 && pt<100) return 0.450726 ;
	  if(pt>=100 && pt<140) return 0.455618 ;
	  if(pt>=140 && pt<200) return 0.462374 ;
	  if(pt>=200 && pt<300) return 0.489544 ;
	  if(pt>=300 && pt<600) return 0.554971 ;
	  if(pt>=600) return 0.678982 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.505972 ;
	  if(pt>=50 && pt<70) return 0.45943 ;
	  if(pt>=70 && pt<100) return 0.445447 ;
	  if(pt>=100 && pt<140) return 0.438127 ;
	  if(pt>=140 && pt<200) return 0.448744 ;
	  if(pt>=200 && pt<300) return 0.48258 ;
	  if(pt>=300 && pt<600) return 0.558924 ;
	  if(pt>=600) return 0.704741 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.16463 ;
	  if(pt>=50 && pt<70) return 0.103549 ;
	  if(pt>=70 && pt<100) return 0.0887529 ;
	  if(pt>=100 && pt<140) return 0.0789445 ;
	  if(pt>=140 && pt<200) return 0.0741228 ;
	  if(pt>=200 && pt<300) return 0.085994 ;
	  if(pt>=300 && pt<600) return 0.14279 ;
	  if(pt>=600) return 0.304453 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.185249 ;
	  if(pt>=50 && pt<70) return 0.117692 ;
	  if(pt>=70 && pt<100) return 0.100407 ;
	  if(pt>=100 && pt<140) return 0.0895698 ;
	  if(pt>=140 && pt<200) return 0.0835016 ;
	  if(pt>=200 && pt<300) return 0.0928656 ;
	  if(pt>=300 && pt<600) return 0.15062 ;
	  if(pt>=600) return 0.323542 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.236811 ;
	  if(pt>=50 && pt<70) return 0.150964 ;
	  if(pt>=70 && pt<100) return 0.13019 ;
	  if(pt>=100 && pt<140) return 0.118657 ;
	  if(pt>=140 && pt<200) return 0.115764 ;
	  if(pt>=200 && pt<300) return 0.134158 ;
	  if(pt>=300 && pt<600) return 0.220298 ;
	  if(pt>=600) return 0.449496 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.321382 ;
	  if(pt>=50 && pt<70) return 0.19381 ;
	  if(pt>=70 && pt<100) return 0.158247 ;
	  if(pt>=100 && pt<140) return 0.148209 ;
	  if(pt>=140 && pt<200) return 0.151157 ;
	  if(pt>=200 && pt<300) return 0.179848 ;
	  if(pt>=300 && pt<600) return 0.279109 ;
	  if(pt>=600) return 0.505407 ;
	}
      }
    }
  else if(sample == "ttx")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.925606 ;
	  if(pt>=50 && pt<70) return 0.932672 ;
	  if(pt>=70 && pt<100) return 0.937036 ;
	  if(pt>=100 && pt<140) return 0.939968 ;
	  if(pt>=140 && pt<200) return 0.943322 ;
	  if(pt>=200 && pt<300) return 0.94703 ;
	  if(pt>=300 && pt<600) return 0.950644 ;
	  if(pt>=600) return 0.949843 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.921456 ;
	  if(pt>=50 && pt<70) return 0.928752 ;
	  if(pt>=70 && pt<100) return 0.932756 ;
	  if(pt>=100 && pt<140) return 0.934242 ;
	  if(pt>=140 && pt<200) return 0.937644 ;
	  if(pt>=200 && pt<300) return 0.941689 ;
	  if(pt>=300 && pt<600) return 0.947111 ;
	  if(pt>=600) return 0.946068 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.906871 ;
	  if(pt>=50 && pt<70) return 0.913164 ;
	  if(pt>=70 && pt<100) return 0.915732 ;
	  if(pt>=100 && pt<140) return 0.917757 ;
	  if(pt>=140 && pt<200) return 0.92213 ;
	  if(pt>=200 && pt<300) return 0.927601 ;
	  if(pt>=300 && pt<600) return 0.934946 ;
	  if(pt>=600) return 0.935787 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.899791 ;
	  if(pt>=50 && pt<70) return 0.903967 ;
	  if(pt>=70 && pt<100) return 0.902179 ;
	  if(pt>=100 && pt<140) return 0.905586 ;
	  if(pt>=140 && pt<200) return 0.910687 ;
	  if(pt>=200 && pt<300) return 0.916539 ;
	  if(pt>=300 && pt<600) return 0.92859 ;
	  if(pt>=600) return 0.911436 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.526629 ;
	  if(pt>=50 && pt<70) return 0.49686 ;
	  if(pt>=70 && pt<100) return 0.47841 ;
	  if(pt>=100 && pt<140) return 0.475382 ;
	  if(pt>=140 && pt<200) return 0.490095 ;
	  if(pt>=200 && pt<300) return 0.52829 ;
	  if(pt>=300 && pt<600) return 0.620782 ;
	  if(pt>=600) return 0.694797 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.524497 ;
	  if(pt>=50 && pt<70) return 0.497375 ;
	  if(pt>=70 && pt<100) return 0.48356 ;
	  if(pt>=100 && pt<140) return 0.479597 ;
	  if(pt>=140 && pt<200) return 0.493133 ;
	  if(pt>=200 && pt<300) return 0.526744 ;
	  if(pt>=300 && pt<600) return 0.614025 ;
	  if(pt>=600) return 0.686717 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.514635 ;
	  if(pt>=50 && pt<70) return 0.488406 ;
	  if(pt>=70 && pt<100) return 0.474932 ;
	  if(pt>=100 && pt<140) return 0.475357 ;
	  if(pt>=140 && pt<200) return 0.486594 ;
	  if(pt>=200 && pt<300) return 0.530321 ;
	  if(pt>=300 && pt<600) return 0.623699 ;
	  if(pt>=600) return 0.739486 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.528191 ;
	  if(pt>=50 && pt<70) return 0.485987 ;
	  if(pt>=70 && pt<100) return 0.466245 ;
	  if(pt>=100 && pt<140) return 0.458355 ;
	  if(pt>=140 && pt<200) return 0.476342 ;
	  if(pt>=200 && pt<300) return 0.526803 ;
	  if(pt>=300 && pt<600) return 0.635699 ;
	  if(pt>=600) return 0.781521 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.141655 ;
	  if(pt>=50 && pt<70) return 0.0905282 ;
	  if(pt>=70 && pt<100) return 0.0785527 ;
	  if(pt>=100 && pt<140) return 0.0742288 ;
	  if(pt>=140 && pt<200) return 0.0811064 ;
	  if(pt>=200 && pt<300) return 0.103187 ;
	  if(pt>=300 && pt<600) return 0.183007 ;
	  if(pt>=600) return 0.323246 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.157233 ;
	  if(pt>=50 && pt<70) return 0.101695 ;
	  if(pt>=70 && pt<100) return 0.0882175 ;
	  if(pt>=100 && pt<140) return 0.0834459 ;
	  if(pt>=140 && pt<200) return 0.089645 ;
	  if(pt>=200 && pt<300) return 0.112165 ;
	  if(pt>=300 && pt<600) return 0.194324 ;
	  if(pt>=600) return 0.348737 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.194544 ;
	  if(pt>=50 && pt<70) return 0.126493 ;
	  if(pt>=70 && pt<100) return 0.111574 ;
	  if(pt>=100 && pt<140) return 0.107752 ;
	  if(pt>=140 && pt<200) return 0.121062 ;
	  if(pt>=200 && pt<300) return 0.156634 ;
	  if(pt>=300 && pt<600) return 0.268459 ;
	  if(pt>=600) return 0.488825 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.260974 ;
	  if(pt>=50 && pt<70) return 0.164021 ;
	  if(pt>=70 && pt<100) return 0.137427 ;
	  if(pt>=100 && pt<140) return 0.133262 ;
	  if(pt>=140 && pt<200) return 0.152857 ;
	  if(pt>=200 && pt<300) return 0.202391 ;
	  if(pt>=300 && pt<600) return 0.332196 ;
	  if(pt>=600) return 0.540172 ;
	}
      }
    }
  else if(sample == "st")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.928132 ;
	  if(pt>=50 && pt<70) return 0.936758 ;
	  if(pt>=70 && pt<100) return 0.941327 ;
	  if(pt>=100 && pt<140) return 0.944032 ;
	  if(pt>=140 && pt<200) return 0.947087 ;
	  if(pt>=200 && pt<300) return 0.950895 ;
	  if(pt>=300 && pt<600) return 0.95508 ;
	  if(pt>=600) return 0.953545 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.923448 ;
	  if(pt>=50 && pt<70) return 0.932552 ;
	  if(pt>=70 && pt<100) return 0.936198 ;
	  if(pt>=100 && pt<140) return 0.93933 ;
	  if(pt>=140 && pt<200) return 0.943055 ;
	  if(pt>=200 && pt<300) return 0.945766 ;
	  if(pt>=300 && pt<600) return 0.949781 ;
	  if(pt>=600) return 0.949118 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.909703 ;
	  if(pt>=50 && pt<70) return 0.917972 ;
	  if(pt>=70 && pt<100) return 0.921439 ;
	  if(pt>=100 && pt<140) return 0.923305 ;
	  if(pt>=140 && pt<200) return 0.927091 ;
	  if(pt>=200 && pt<300) return 0.930952 ;
	  if(pt>=300 && pt<600) return 0.936525 ;
	  if(pt>=600) return 0.944971 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.903022 ;
	  if(pt>=50 && pt<70) return 0.905424 ;
	  if(pt>=70 && pt<100) return 0.9077 ;
	  if(pt>=100 && pt<140) return 0.910166 ;
	  if(pt>=140 && pt<200) return 0.913411 ;
	  if(pt>=200 && pt<300) return 0.916477 ;
	  if(pt>=300 && pt<600) return 0.929294 ;
	  if(pt>=600) return 0.933388 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.519589 ;
	  if(pt>=50 && pt<70) return 0.487563 ;
	  if(pt>=70 && pt<100) return 0.472168 ;
	  if(pt>=100 && pt<140) return 0.467855 ;
	  if(pt>=140 && pt<200) return 0.478362 ;
	  if(pt>=200 && pt<300) return 0.508546 ;
	  if(pt>=300 && pt<600) return 0.594934 ;
	  if(pt>=600) return 0.659036 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.515003 ;
	  if(pt>=50 && pt<70) return 0.492236 ;
	  if(pt>=70 && pt<100) return 0.47787 ;
	  if(pt>=100 && pt<140) return 0.46883 ;
	  if(pt>=140 && pt<200) return 0.479881 ;
	  if(pt>=200 && pt<300) return 0.500357 ;
	  if(pt>=300 && pt<600) return 0.581794 ;
	  if(pt>=600) return 0.654855 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.506228 ;
	  if(pt>=50 && pt<70) return 0.485101 ;
	  if(pt>=70 && pt<100) return 0.473494 ;
	  if(pt>=100 && pt<140) return 0.466281 ;
	  if(pt>=140 && pt<200) return 0.479171 ;
	  if(pt>=200 && pt<300) return 0.495515 ;
	  if(pt>=300 && pt<600) return 0.591365 ;
	  if(pt>=600) return 0.686829 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.516708 ;
	  if(pt>=50 && pt<70) return 0.474475 ;
	  if(pt>=70 && pt<100) return 0.455315 ;
	  if(pt>=100 && pt<140) return 0.448262 ;
	  if(pt>=140 && pt<200) return 0.459588 ;
	  if(pt>=200 && pt<300) return 0.478655 ;
	  if(pt>=300 && pt<600) return 0.563391 ;
	  if(pt>=600) return 0.782598 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.163323 ;
	  if(pt>=50 && pt<70) return 0.101563 ;
	  if(pt>=70 && pt<100) return 0.0853963 ;
	  if(pt>=100 && pt<140) return 0.0789168 ;
	  if(pt>=140 && pt<200) return 0.0792538 ;
	  if(pt>=200 && pt<300) return 0.092507 ;
	  if(pt>=300 && pt<600) return 0.16012 ;
	  if(pt>=600) return 0.301679 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.175684 ;
	  if(pt>=50 && pt<70) return 0.110451 ;
	  if(pt>=70 && pt<100) return 0.0932369 ;
	  if(pt>=100 && pt<140) return 0.0859455 ;
	  if(pt>=140 && pt<200) return 0.0862756 ;
	  if(pt>=200 && pt<300) return 0.0960885 ;
	  if(pt>=300 && pt<600) return 0.162116 ;
	  if(pt>=600) return 0.32497 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.210401 ;
	  if(pt>=50 && pt<70) return 0.134437 ;
	  if(pt>=70 && pt<100) return 0.115123 ;
	  if(pt>=100 && pt<140) return 0.108368 ;
	  if(pt>=140 && pt<200) return 0.116546 ;
	  if(pt>=200 && pt<300) return 0.132224 ;
	  if(pt>=300 && pt<600) return 0.225403 ;
	  if(pt>=600) return 0.454104 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.273932 ;
	  if(pt>=50 && pt<70) return 0.168223 ;
	  if(pt>=70 && pt<100) return 0.140219 ;
	  if(pt>=100 && pt<140) return 0.135691 ;
	  if(pt>=140 && pt<200) return 0.152693 ;
	  if(pt>=200 && pt<300) return 0.178629 ;
	  if(pt>=300 && pt<600) return 0.284006 ;
	  if(pt>=600) return 0.51334 ;
	}
      }
    }
  else if(sample == "tt")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.928489 ;
	  if(pt>=50 && pt<70) return 0.936358 ;
	  if(pt>=70 && pt<100) return 0.939329 ;
	  if(pt>=100 && pt<140) return 0.942312 ;
	  if(pt>=140 && pt<200) return 0.946494 ;
	  if(pt>=200 && pt<300) return 0.949071 ;
	  if(pt>=300 && pt<600) return 0.94913 ;
	  if(pt>=600) return 0.94744 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.923875 ;
	  if(pt>=50 && pt<70) return 0.931643 ;
	  if(pt>=70 && pt<100) return 0.934944 ;
	  if(pt>=100 && pt<140) return 0.937767 ;
	  if(pt>=140 && pt<200) return 0.942025 ;
	  if(pt>=200 && pt<300) return 0.944155 ;
	  if(pt>=300 && pt<600) return 0.943651 ;
	  if(pt>=600) return 0.943732 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.908884 ;
	  if(pt>=50 && pt<70) return 0.916452 ;
	  if(pt>=70 && pt<100) return 0.919043 ;
	  if(pt>=100 && pt<140) return 0.921359 ;
	  if(pt>=140 && pt<200) return 0.92691 ;
	  if(pt>=200 && pt<300) return 0.929786 ;
	  if(pt>=300 && pt<600) return 0.929737 ;
	  if(pt>=600) return 0.938173 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.902435 ;
	  if(pt>=50 && pt<70) return 0.905362 ;
	  if(pt>=70 && pt<100) return 0.90628 ;
	  if(pt>=100 && pt<140) return 0.908942 ;
	  if(pt>=140 && pt<200) return 0.914433 ;
	  if(pt>=200 && pt<300) return 0.917965 ;
	  if(pt>=300 && pt<600) return 0.918255 ;
	  if(pt>=600) return 0.922557 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.524975 ;
	  if(pt>=50 && pt<70) return 0.494497 ;
	  if(pt>=70 && pt<100) return 0.477375 ;
	  if(pt>=100 && pt<140) return 0.472183 ;
	  if(pt>=140 && pt<200) return 0.485762 ;
	  if(pt>=200 && pt<300) return 0.517724 ;
	  if(pt>=300 && pt<600) return 0.592044 ;
	  if(pt>=600) return 0.649878 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.521824 ;
	  if(pt>=50 && pt<70) return 0.495354 ;
	  if(pt>=70 && pt<100) return 0.479801 ;
	  if(pt>=100 && pt<140) return 0.475446 ;
	  if(pt>=140 && pt<200) return 0.485749 ;
	  if(pt>=200 && pt<300) return 0.511998 ;
	  if(pt>=300 && pt<600) return 0.583933 ;
	  if(pt>=600) return 0.647665 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.512644 ;
	  if(pt>=50 && pt<70) return 0.485937 ;
	  if(pt>=70 && pt<100) return 0.474412 ;
	  if(pt>=100 && pt<140) return 0.470773 ;
	  if(pt>=140 && pt<200) return 0.484986 ;
	  if(pt>=200 && pt<300) return 0.510647 ;
	  if(pt>=300 && pt<600) return 0.591648 ;
	  if(pt>=600) return 0.685138 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.525138 ;
	  if(pt>=50 && pt<70) return 0.484251 ;
	  if(pt>=70 && pt<100) return 0.461788 ;
	  if(pt>=100 && pt<140) return 0.456025 ;
	  if(pt>=140 && pt<200) return 0.470075 ;
	  if(pt>=200 && pt<300) return 0.503388 ;
	  if(pt>=300 && pt<600) return 0.590552 ;
	  if(pt>=600) return 0.710843 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.14031 ;
	  if(pt>=50 && pt<70) return 0.0894484 ;
	  if(pt>=70 && pt<100) return 0.0761991 ;
	  if(pt>=100 && pt<140) return 0.0709635 ;
	  if(pt>=140 && pt<200) return 0.0751688 ;
	  if(pt>=200 && pt<300) return 0.0908985 ;
	  if(pt>=300 && pt<600) return 0.151593 ;
	  if(pt>=600) return 0.301989 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.154074 ;
	  if(pt>=50 && pt<70) return 0.099789 ;
	  if(pt>=70 && pt<100) return 0.0852958 ;
	  if(pt>=100 && pt<140) return 0.0792641 ;
	  if(pt>=140 && pt<200) return 0.0824022 ;
	  if(pt>=200 && pt<300) return 0.0969287 ;
	  if(pt>=300 && pt<600) return 0.1592 ;
	  if(pt>=600) return 0.325947 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.190491 ;
	  if(pt>=50 && pt<70) return 0.124391 ;
	  if(pt>=70 && pt<100) return 0.107716 ;
	  if(pt>=100 && pt<140) return 0.102422 ;
	  if(pt>=140 && pt<200) return 0.110677 ;
	  if(pt>=200 && pt<300) return 0.133626 ;
	  if(pt>=300 && pt<600) return 0.227186 ;
	  if(pt>=600) return 0.458928 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.257815 ;
	  if(pt>=50 && pt<70) return 0.161815 ;
	  if(pt>=70 && pt<100) return 0.134995 ;
	  if(pt>=100 && pt<140) return 0.129469 ;
	  if(pt>=140 && pt<200) return 0.145057 ;
	  if(pt>=200 && pt<300) return 0.181222 ;
	  if(pt>=300 && pt<600) return 0.291998 ;
	  if(pt>=600) return 0.519698 ;
	}
      }
    }
  else if(sample == "dib")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.919095 ;
	  if(pt>=50 && pt<70) return 0.929909 ;
	  if(pt>=70 && pt<100) return 0.930533 ;
	  if(pt>=100 && pt<140) return 0.926632 ;
	  if(pt>=140 && pt<200) return 0.93491 ;
	  if(pt>=200 && pt<300) return 0.931081 ;
	  if(pt>=300 && pt<600) return 0.947711 ;
	  if(pt>=600) return 0.968265 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.91066 ;
	  if(pt>=50 && pt<70) return 0.929246 ;
	  if(pt>=70 && pt<100) return 0.929913 ;
	  if(pt>=100 && pt<140) return 0.932934 ;
	  if(pt>=140 && pt<200) return 0.926046 ;
	  if(pt>=200 && pt<300) return 0.925153 ;
	  if(pt>=300 && pt<600) return 0.946744 ;
	  if(pt>=600) return 0.975925 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.902229 ;
	  if(pt>=50 && pt<70) return 0.914489 ;
	  if(pt>=70 && pt<100) return 0.909306 ;
	  if(pt>=100 && pt<140) return 0.918275 ;
	  if(pt>=140 && pt<200) return 0.914618 ;
	  if(pt>=200 && pt<300) return 0.918679 ;
	  if(pt>=300 && pt<600) return 0.934546 ;
	  if(pt>=600) return 0.951759 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.891679 ;
	  if(pt>=50 && pt<70) return 0.89782 ;
	  if(pt>=70 && pt<100) return 0.890734 ;
	  if(pt>=100 && pt<140) return 0.899401 ;
	  if(pt>=140 && pt<200) return 0.909105 ;
	  if(pt>=200 && pt<300) return 0.911866 ;
	  if(pt>=300 && pt<600) return 0.93065 ;
	  if(pt>=600) return 1 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.520367 ;
	  if(pt>=50 && pt<70) return 0.493208 ;
	  if(pt>=70 && pt<100) return 0.472977 ;
	  if(pt>=100 && pt<140) return 0.460993 ;
	  if(pt>=140 && pt<200) return 0.478527 ;
	  if(pt>=200 && pt<300) return 0.512678 ;
	  if(pt>=300 && pt<600) return 0.601098 ;
	  if(pt>=600) return 0.658976 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.521075 ;
	  if(pt>=50 && pt<70) return 0.492432 ;
	  if(pt>=70 && pt<100) return 0.483858 ;
	  if(pt>=100 && pt<140) return 0.471868 ;
	  if(pt>=140 && pt<200) return 0.482187 ;
	  if(pt>=200 && pt<300) return 0.501855 ;
	  if(pt>=300 && pt<600) return 0.595792 ;
	  if(pt>=600) return 0.680188 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.514471 ;
	  if(pt>=50 && pt<70) return 0.487848 ;
	  if(pt>=70 && pt<100) return 0.476232 ;
	  if(pt>=100 && pt<140) return 0.476894 ;
	  if(pt>=140 && pt<200) return 0.481134 ;
	  if(pt>=200 && pt<300) return 0.519426 ;
	  if(pt>=300 && pt<600) return 0.614791 ;
	  if(pt>=600) return 0.71972 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.526358 ;
	  if(pt>=50 && pt<70) return 0.497467 ;
	  if(pt>=70 && pt<100) return 0.47235 ;
	  if(pt>=100 && pt<140) return 0.458112 ;
	  if(pt>=140 && pt<200) return 0.483114 ;
	  if(pt>=200 && pt<300) return 0.515958 ;
	  if(pt>=300 && pt<600) return 0.629878 ;
	  if(pt>=600) return 0.713483 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.14036 ;
	  if(pt>=50 && pt<70) return 0.0857675 ;
	  if(pt>=70 && pt<100) return 0.0749535 ;
	  if(pt>=100 && pt<140) return 0.0701135 ;
	  if(pt>=140 && pt<200) return 0.0711584 ;
	  if(pt>=200 && pt<300) return 0.0848574 ;
	  if(pt>=300 && pt<600) return 0.152285 ;
	  if(pt>=600) return 0.310082 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.155529 ;
	  if(pt>=50 && pt<70) return 0.0995451 ;
	  if(pt>=70 && pt<100) return 0.0844271 ;
	  if(pt>=100 && pt<140) return 0.0791756 ;
	  if(pt>=140 && pt<200) return 0.0805337 ;
	  if(pt>=200 && pt<300) return 0.0919421 ;
	  if(pt>=300 && pt<600) return 0.160375 ;
	  if(pt>=600) return 0.319832 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.19542 ;
	  if(pt>=50 && pt<70) return 0.123868 ;
	  if(pt>=70 && pt<100) return 0.111322 ;
	  if(pt>=100 && pt<140) return 0.105929 ;
	  if(pt>=140 && pt<200) return 0.115126 ;
	  if(pt>=200 && pt<300) return 0.132777 ;
	  if(pt>=300 && pt<600) return 0.230198 ;
	  if(pt>=600) return 0.433242 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.26896 ;
	  if(pt>=50 && pt<70) return 0.165431 ;
	  if(pt>=70 && pt<100) return 0.14009 ;
	  if(pt>=100 && pt<140) return 0.135418 ;
	  if(pt>=140 && pt<200) return 0.146241 ;
	  if(pt>=200 && pt<300) return 0.17545 ;
	  if(pt>=300 && pt<600) return 0.287827 ;
	  if(pt>=600) return 0.459712 ;
	}
      }
    }
  return 1.0;
}

double BTag_MCEfficiency_M(string sample, int flavor, double pt, double eta){
  if(sample == "dy")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.788439 ;
	  if(pt>=50 && pt<70) return 0.819981 ;
	  if(pt>=70 && pt<100) return 0.8305 ;
	  if(pt>=100 && pt<140) return 0.829388 ;
	  if(pt>=140 && pt<200) return 0.849152 ;
	  if(pt>=200 && pt<300) return 0.846761 ;
	  if(pt>=300 && pt<600) return 0.836956 ;
	  if(pt>=600) return 0.829924 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.766004 ;
	  if(pt>=50 && pt<70) return 0.801194 ;
	  if(pt>=70 && pt<100) return 0.824395 ;
	  if(pt>=100 && pt<140) return 0.822347 ;
	  if(pt>=140 && pt<200) return 0.83768 ;
	  if(pt>=200 && pt<300) return 0.830082 ;
	  if(pt>=300 && pt<600) return 0.825308 ;
	  if(pt>=600) return 0.805405 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.715652 ;
	  if(pt>=50 && pt<70) return 0.748716 ;
	  if(pt>=70 && pt<100) return 0.750594 ;
	  if(pt>=100 && pt<140) return 0.784126 ;
	  if(pt>=140 && pt<200) return 0.776026 ;
	  if(pt>=200 && pt<300) return 0.793947 ;
	  if(pt>=300 && pt<600) return 0.758987 ;
	  if(pt>=600) return 0.731354 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.699504 ;
	  if(pt>=50 && pt<70) return 0.712544 ;
	  if(pt>=70 && pt<100) return 0.738971 ;
	  if(pt>=100 && pt<140) return 0.704965 ;
	  if(pt>=140 && pt<200) return 0.722872 ;
	  if(pt>=200 && pt<300) return 0.747836 ;
	  if(pt>=300 && pt<600) return 0.735647 ;
	  if(pt>=600) return 0.694327 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.192278 ;
	  if(pt>=50 && pt<70) return 0.175704 ;
	  if(pt>=70 && pt<100) return 0.162206 ;
	  if(pt>=100 && pt<140) return 0.174202 ;
	  if(pt>=140 && pt<200) return 0.162764 ;
	  if(pt>=200 && pt<300) return 0.182075 ;
	  if(pt>=300 && pt<600) return 0.232246 ;
	  if(pt>=600) return 0.252569 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.187912 ;
	  if(pt>=50 && pt<70) return 0.177296 ;
	  if(pt>=70 && pt<100) return 0.158728 ;
	  if(pt>=100 && pt<140) return 0.15943 ;
	  if(pt>=140 && pt<200) return 0.164303 ;
	  if(pt>=200 && pt<300) return 0.182423 ;
	  if(pt>=300 && pt<600) return 0.216468 ;
	  if(pt>=600) return 0.282186 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.174285 ;
	  if(pt>=50 && pt<70) return 0.161508 ;
	  if(pt>=70 && pt<100) return 0.160048 ;
	  if(pt>=100 && pt<140) return 0.156912 ;
	  if(pt>=140 && pt<200) return 0.17182 ;
	  if(pt>=200 && pt<300) return 0.173718 ;
	  if(pt>=300 && pt<600) return 0.195849 ;
	  if(pt>=600) return 0.235279 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.159618 ;
	  if(pt>=50 && pt<70) return 0.14814 ;
	  if(pt>=70 && pt<100) return 0.138957 ;
	  if(pt>=100 && pt<140) return 0.13481 ;
	  if(pt>=140 && pt<200) return 0.138906 ;
	  if(pt>=200 && pt<300) return 0.170971 ;
	  if(pt>=300 && pt<600) return 0.21268 ;
	  if(pt>=600) return 0.220715 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0338363 ;
	  if(pt>=50 && pt<70) return 0.0198584 ;
	  if(pt>=70 && pt<100) return 0.0171191 ;
	  if(pt>=100 && pt<140) return 0.0153923 ;
	  if(pt>=140 && pt<200) return 0.0136957 ;
	  if(pt>=200 && pt<300) return 0.0139145 ;
	  if(pt>=300 && pt<600) return 0.0204948 ;
	  if(pt>=600) return 0.0584709 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0360152 ;
	  if(pt>=50 && pt<70) return 0.0231707 ;
	  if(pt>=70 && pt<100) return 0.019403 ;
	  if(pt>=100 && pt<140) return 0.0175199 ;
	  if(pt>=140 && pt<200) return 0.0153178 ;
	  if(pt>=200 && pt<300) return 0.0142875 ;
	  if(pt>=300 && pt<600) return 0.0219566 ;
	  if(pt>=600) return 0.0603255 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.044138 ;
	  if(pt>=50 && pt<70) return 0.0295545 ;
	  if(pt>=70 && pt<100) return 0.0272254 ;
	  if(pt>=100 && pt<140) return 0.0224177 ;
	  if(pt>=140 && pt<200) return 0.0195233 ;
	  if(pt>=200 && pt<300) return 0.0198841 ;
	  if(pt>=300 && pt<600) return 0.0333882 ;
	  if(pt>=600) return 0.0890321 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0564224 ;
	  if(pt>=50 && pt<70) return 0.0321335 ;
	  if(pt>=70 && pt<100) return 0.0295463 ;
	  if(pt>=100 && pt<140) return 0.0265831 ;
	  if(pt>=140 && pt<200) return 0.0244293 ;
	  if(pt>=200 && pt<300) return 0.0270997 ;
	  if(pt>=300 && pt<600) return 0.043592 ;
	  if(pt>=600) return 0.0747513 ;
	}
      }
    }
  else if(sample == "qcd")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.786142 ;
	  if(pt>=50 && pt<70) return 0.814186 ;
	  if(pt>=70 && pt<100) return 0.826588 ;
	  if(pt>=100 && pt<140) return 0.83243 ;
	  if(pt>=140 && pt<200) return 0.839812 ;
	  if(pt>=200 && pt<300) return 0.852361 ;
	  if(pt>=300 && pt<600) return 0.848049 ;
	  if(pt>=600) return 0.824806 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.774172 ;
	  if(pt>=50 && pt<70) return 0.801184 ;
	  if(pt>=70 && pt<100) return 0.813602 ;
	  if(pt>=100 && pt<140) return 0.817657 ;
	  if(pt>=140 && pt<200) return 0.827111 ;
	  if(pt>=200 && pt<300) return 0.837876 ;
	  if(pt>=300 && pt<600) return 0.828362 ;
	  if(pt>=600) return 0.798402 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.731658 ;
	  if(pt>=50 && pt<70) return 0.764152 ;
	  if(pt>=70 && pt<100) return 0.772209 ;
	  if(pt>=100 && pt<140) return 0.778511 ;
	  if(pt>=140 && pt<200) return 0.784003 ;
	  if(pt>=200 && pt<300) return 0.793285 ;
	  if(pt>=300 && pt<600) return 0.775857 ;
	  if(pt>=600) return 0.735846 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.702141 ;
	  if(pt>=50 && pt<70) return 0.728536 ;
	  if(pt>=70 && pt<100) return 0.732823 ;
	  if(pt>=100 && pt<140) return 0.735174 ;
	  if(pt>=140 && pt<200) return 0.737238 ;
	  if(pt>=200 && pt<300) return 0.742886 ;
	  if(pt>=300 && pt<600) return 0.709444 ;
	  if(pt>=600) return 0.687201 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.188666 ;
	  if(pt>=50 && pt<70) return 0.173627 ;
	  if(pt>=70 && pt<100) return 0.172337 ;
	  if(pt>=100 && pt<140) return 0.176422 ;
	  if(pt>=140 && pt<200) return 0.188635 ;
	  if(pt>=200 && pt<300) return 0.205055 ;
	  if(pt>=300 && pt<600) return 0.240291 ;
	  if(pt>=600) return 0.275866 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.184609 ;
	  if(pt>=50 && pt<70) return 0.1711 ;
	  if(pt>=70 && pt<100) return 0.171769 ;
	  if(pt>=100 && pt<140) return 0.174923 ;
	  if(pt>=140 && pt<200) return 0.184919 ;
	  if(pt>=200 && pt<300) return 0.200295 ;
	  if(pt>=300 && pt<600) return 0.228509 ;
	  if(pt>=600) return 0.257915 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.175878 ;
	  if(pt>=50 && pt<70) return 0.163095 ;
	  if(pt>=70 && pt<100) return 0.160661 ;
	  if(pt>=100 && pt<140) return 0.163074 ;
	  if(pt>=140 && pt<200) return 0.171869 ;
	  if(pt>=200 && pt<300) return 0.18383 ;
	  if(pt>=300 && pt<600) return 0.206844 ;
	  if(pt>=600) return 0.243153 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.166572 ;
	  if(pt>=50 && pt<70) return 0.153152 ;
	  if(pt>=70 && pt<100) return 0.144691 ;
	  if(pt>=100 && pt<140) return 0.142847 ;
	  if(pt>=140 && pt<200) return 0.149761 ;
	  if(pt>=200 && pt<300) return 0.162117 ;
	  if(pt>=300 && pt<600) return 0.185562 ;
	  if(pt>=600) return 0.228545 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0142195 ;
	  if(pt>=50 && pt<70) return 0.0090539 ;
	  if(pt>=70 && pt<100) return 0.00774963 ;
	  if(pt>=100 && pt<140) return 0.00665534 ;
	  if(pt>=140 && pt<200) return 0.00695749 ;
	  if(pt>=200 && pt<300) return 0.008249 ;
	  if(pt>=300 && pt<600) return 0.0142748 ;
	  if(pt>=600) return 0.0479807 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0153347 ;
	  if(pt>=50 && pt<70) return 0.00991585 ;
	  if(pt>=70 && pt<100) return 0.0084633 ;
	  if(pt>=100 && pt<140) return 0.00791183 ;
	  if(pt>=140 && pt<200) return 0.00814272 ;
	  if(pt>=200 && pt<300) return 0.00880325 ;
	  if(pt>=300 && pt<600) return 0.015032 ;
	  if(pt>=600) return 0.0499368 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0192045 ;
	  if(pt>=50 && pt<70) return 0.0132933 ;
	  if(pt>=70 && pt<100) return 0.0120319 ;
	  if(pt>=100 && pt<140) return 0.011021 ;
	  if(pt>=140 && pt<200) return 0.0111536 ;
	  if(pt>=200 && pt<300) return 0.0128595 ;
	  if(pt>=300 && pt<600) return 0.0255242 ;
	  if(pt>=600) return 0.0812452 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0294708 ;
	  if(pt>=50 && pt<70) return 0.0192343 ;
	  if(pt>=70 && pt<100) return 0.0165029 ;
	  if(pt>=100 && pt<140) return 0.0162242 ;
	  if(pt>=140 && pt<200) return 0.017306 ;
	  if(pt>=200 && pt<300) return 0.0210505 ;
	  if(pt>=300 && pt<600) return 0.0380773 ;
	  if(pt>=600) return 0.0823214 ;
	}
      }
    }
  else if(sample == "wj")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.781536 ;
	  if(pt>=50 && pt<70) return 0.805753 ;
	  if(pt>=70 && pt<100) return 0.811077 ;
	  if(pt>=100 && pt<140) return 0.815027 ;
	  if(pt>=140 && pt<200) return 0.811923 ;
	  if(pt>=200 && pt<300) return 0.821861 ;
	  if(pt>=300 && pt<600) return 0.824519 ;
	  if(pt>=600) return 0.815119 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.76887 ;
	  if(pt>=50 && pt<70) return 0.794995 ;
	  if(pt>=70 && pt<100) return 0.804123 ;
	  if(pt>=100 && pt<140) return 0.798826 ;
	  if(pt>=140 && pt<200) return 0.802232 ;
	  if(pt>=200 && pt<300) return 0.805901 ;
	  if(pt>=300 && pt<600) return 0.805633 ;
	  if(pt>=600) return 0.791605 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.726117 ;
	  if(pt>=50 && pt<70) return 0.763066 ;
	  if(pt>=70 && pt<100) return 0.760379 ;
	  if(pt>=100 && pt<140) return 0.757634 ;
	  if(pt>=140 && pt<200) return 0.751598 ;
	  if(pt>=200 && pt<300) return 0.750978 ;
	  if(pt>=300 && pt<600) return 0.747376 ;
	  if(pt>=600) return 0.72565 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.695328 ;
	  if(pt>=50 && pt<70) return 0.716743 ;
	  if(pt>=70 && pt<100) return 0.712366 ;
	  if(pt>=100 && pt<140) return 0.711298 ;
	  if(pt>=140 && pt<200) return 0.702145 ;
	  if(pt>=200 && pt<300) return 0.704393 ;
	  if(pt>=300 && pt<600) return 0.696708 ;
	  if(pt>=600) return 0.668938 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.188022 ;
	  if(pt>=50 && pt<70) return 0.171087 ;
	  if(pt>=70 && pt<100) return 0.166339 ;
	  if(pt>=100 && pt<140) return 0.160366 ;
	  if(pt>=140 && pt<200) return 0.167858 ;
	  if(pt>=200 && pt<300) return 0.188833 ;
	  if(pt>=300 && pt<600) return 0.231108 ;
	  if(pt>=600) return 0.275936 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.184053 ;
	  if(pt>=50 && pt<70) return 0.171958 ;
	  if(pt>=70 && pt<100) return 0.166664 ;
	  if(pt>=100 && pt<140) return 0.16319 ;
	  if(pt>=140 && pt<200) return 0.167898 ;
	  if(pt>=200 && pt<300) return 0.186276 ;
	  if(pt>=300 && pt<600) return 0.221515 ;
	  if(pt>=600) return 0.258639 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.177663 ;
	  if(pt>=50 && pt<70) return 0.161756 ;
	  if(pt>=70 && pt<100) return 0.156033 ;
	  if(pt>=100 && pt<140) return 0.158522 ;
	  if(pt>=140 && pt<200) return 0.162398 ;
	  if(pt>=200 && pt<300) return 0.179004 ;
	  if(pt>=300 && pt<600) return 0.208811 ;
	  if(pt>=600) return 0.246884 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.167256 ;
	  if(pt>=50 && pt<70) return 0.15005 ;
	  if(pt>=70 && pt<100) return 0.142838 ;
	  if(pt>=100 && pt<140) return 0.141105 ;
	  if(pt>=140 && pt<200) return 0.146175 ;
	  if(pt>=200 && pt<300) return 0.161676 ;
	  if(pt>=300 && pt<600) return 0.187246 ;
	  if(pt>=600) return 0.229983 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0238783 ;
	  if(pt>=50 && pt<70) return 0.0130375 ;
	  if(pt>=70 && pt<100) return 0.0104439 ;
	  if(pt>=100 && pt<140) return 0.00912547 ;
	  if(pt>=140 && pt<200) return 0.00830632 ;
	  if(pt>=200 && pt<300) return 0.00993866 ;
	  if(pt>=300 && pt<600) return 0.0193107 ;
	  if(pt>=600) return 0.0567631 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0275034 ;
	  if(pt>=50 && pt<70) return 0.0156884 ;
	  if(pt>=70 && pt<100) return 0.0126489 ;
	  if(pt>=100 && pt<140) return 0.0107182 ;
	  if(pt>=140 && pt<200) return 0.00990751 ;
	  if(pt>=200 && pt<300) return 0.0109258 ;
	  if(pt>=300 && pt<600) return 0.0199511 ;
	  if(pt>=600) return 0.0565818 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0342674 ;
	  if(pt>=50 && pt<70) return 0.0202625 ;
	  if(pt>=70 && pt<100) return 0.0165357 ;
	  if(pt>=100 && pt<140) return 0.0149187 ;
	  if(pt>=140 && pt<200) return 0.0143112 ;
	  if(pt>=200 && pt<300) return 0.0170554 ;
	  if(pt>=300 && pt<600) return 0.0319311 ;
	  if(pt>=600) return 0.0862961 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0464005 ;
	  if(pt>=50 && pt<70) return 0.0261784 ;
	  if(pt>=70 && pt<100) return 0.0200199 ;
	  if(pt>=100 && pt<140) return 0.0193049 ;
	  if(pt>=140 && pt<200) return 0.0201594 ;
	  if(pt>=200 && pt<300) return 0.024464 ;
	  if(pt>=300 && pt<600) return 0.0420911 ;
	  if(pt>=600) return 0.0892878 ;
	}
      }
    }
  else if(sample == "ttx")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.810699 ;
	  if(pt>=50 && pt<70) return 0.83219 ;
	  if(pt>=70 && pt<100) return 0.840355 ;
	  if(pt>=100 && pt<140) return 0.847661 ;
	  if(pt>=140 && pt<200) return 0.854143 ;
	  if(pt>=200 && pt<300) return 0.858623 ;
	  if(pt>=300 && pt<600) return 0.850054 ;
	  if(pt>=600) return 0.804608 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.797587 ;
	  if(pt>=50 && pt<70) return 0.818771 ;
	  if(pt>=70 && pt<100) return 0.829672 ;
	  if(pt>=100 && pt<140) return 0.835647 ;
	  if(pt>=140 && pt<200) return 0.840937 ;
	  if(pt>=200 && pt<300) return 0.844154 ;
	  if(pt>=300 && pt<600) return 0.833636 ;
	  if(pt>=600) return 0.775966 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.754138 ;
	  if(pt>=50 && pt<70) return 0.77777 ;
	  if(pt>=70 && pt<100) return 0.786439 ;
	  if(pt>=100 && pt<140) return 0.790947 ;
	  if(pt>=140 && pt<200) return 0.798017 ;
	  if(pt>=200 && pt<300) return 0.79891 ;
	  if(pt>=300 && pt<600) return 0.780106 ;
	  if(pt>=600) return 0.703957 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.719427 ;
	  if(pt>=50 && pt<70) return 0.746062 ;
	  if(pt>=70 && pt<100) return 0.751058 ;
	  if(pt>=100 && pt<140) return 0.759745 ;
	  if(pt>=140 && pt<200) return 0.765623 ;
	  if(pt>=200 && pt<300) return 0.767898 ;
	  if(pt>=300 && pt<600) return 0.737757 ;
	  if(pt>=600) return 0.648871 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.190312 ;
	  if(pt>=50 && pt<70) return 0.159433 ;
	  if(pt>=70 && pt<100) return 0.150834 ;
	  if(pt>=100 && pt<140) return 0.14886 ;
	  if(pt>=140 && pt<200) return 0.160974 ;
	  if(pt>=200 && pt<300) return 0.185941 ;
	  if(pt>=300 && pt<600) return 0.252934 ;
	  if(pt>=600) return 0.299262 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.193775 ;
	  if(pt>=50 && pt<70) return 0.161194 ;
	  if(pt>=70 && pt<100) return 0.153 ;
	  if(pt>=100 && pt<140) return 0.155073 ;
	  if(pt>=140 && pt<200) return 0.164227 ;
	  if(pt>=200 && pt<300) return 0.187187 ;
	  if(pt>=300 && pt<600) return 0.248341 ;
	  if(pt>=600) return 0.297072 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.186254 ;
	  if(pt>=50 && pt<70) return 0.160366 ;
	  if(pt>=70 && pt<100) return 0.151476 ;
	  if(pt>=100 && pt<140) return 0.152715 ;
	  if(pt>=140 && pt<200) return 0.163174 ;
	  if(pt>=200 && pt<300) return 0.185326 ;
	  if(pt>=300 && pt<600) return 0.239528 ;
	  if(pt>=600) return 0.303139 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.177328 ;
	  if(pt>=50 && pt<70) return 0.14795 ;
	  if(pt>=70 && pt<100) return 0.135348 ;
	  if(pt>=100 && pt<140) return 0.136486 ;
	  if(pt>=140 && pt<200) return 0.152645 ;
	  if(pt>=200 && pt<300) return 0.176127 ;
	  if(pt>=300 && pt<600) return 0.224789 ;
	  if(pt>=600) return 0.242087 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.01917 ;
	  if(pt>=50 && pt<70) return 0.0110702 ;
	  if(pt>=70 && pt<100) return 0.00902598 ;
	  if(pt>=100 && pt<140) return 0.00865529 ;
	  if(pt>=140 && pt<200) return 0.00965998 ;
	  if(pt>=200 && pt<300) return 0.0121528 ;
	  if(pt>=300 && pt<600) return 0.0243258 ;
	  if(pt>=600) return 0.0562305 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0212073 ;
	  if(pt>=50 && pt<70) return 0.0125054 ;
	  if(pt>=70 && pt<100) return 0.0104047 ;
	  if(pt>=100 && pt<140) return 0.0101168 ;
	  if(pt>=140 && pt<200) return 0.0105605 ;
	  if(pt>=200 && pt<300) return 0.0137805 ;
	  if(pt>=300 && pt<600) return 0.0250807 ;
	  if(pt>=600) return 0.0592792 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0241524 ;
	  if(pt>=50 && pt<70) return 0.0146144 ;
	  if(pt>=70 && pt<100) return 0.0125163 ;
	  if(pt>=100 && pt<140) return 0.0123158 ;
	  if(pt>=140 && pt<200) return 0.0144975 ;
	  if(pt>=200 && pt<300) return 0.0195762 ;
	  if(pt>=300 && pt<600) return 0.0367029 ;
	  if(pt>=600) return 0.0878885 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0312554 ;
	  if(pt>=50 && pt<70) return 0.0189097 ;
	  if(pt>=70 && pt<100) return 0.0160487 ;
	  if(pt>=100 && pt<140) return 0.016407 ;
	  if(pt>=140 && pt<200) return 0.0195108 ;
	  if(pt>=200 && pt<300) return 0.0267428 ;
	  if(pt>=300 && pt<600) return 0.0527937 ;
	  if(pt>=600) return 0.105282 ;
	}
      }
    }
  else if(sample == "st")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.815363 ;
	  if(pt>=50 && pt<70) return 0.837827 ;
	  if(pt>=70 && pt<100) return 0.847572 ;
	  if(pt>=100 && pt<140) return 0.855455 ;
	  if(pt>=140 && pt<200) return 0.863129 ;
	  if(pt>=200 && pt<300) return 0.868603 ;
	  if(pt>=300 && pt<600) return 0.864177 ;
	  if(pt>=600) return 0.816161 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.8027 ;
	  if(pt>=50 && pt<70) return 0.825979 ;
	  if(pt>=70 && pt<100) return 0.836009 ;
	  if(pt>=100 && pt<140) return 0.844524 ;
	  if(pt>=140 && pt<200) return 0.853085 ;
	  if(pt>=200 && pt<300) return 0.856609 ;
	  if(pt>=300 && pt<600) return 0.846366 ;
	  if(pt>=600) return 0.787529 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.761572 ;
	  if(pt>=50 && pt<70) return 0.784737 ;
	  if(pt>=70 && pt<100) return 0.794169 ;
	  if(pt>=100 && pt<140) return 0.801162 ;
	  if(pt>=140 && pt<200) return 0.810278 ;
	  if(pt>=200 && pt<300) return 0.811894 ;
	  if(pt>=300 && pt<600) return 0.792393 ;
	  if(pt>=600) return 0.715515 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.725922 ;
	  if(pt>=50 && pt<70) return 0.748433 ;
	  if(pt>=70 && pt<100) return 0.757215 ;
	  if(pt>=100 && pt<140) return 0.765658 ;
	  if(pt>=140 && pt<200) return 0.774759 ;
	  if(pt>=200 && pt<300) return 0.776012 ;
	  if(pt>=300 && pt<600) return 0.763504 ;
	  if(pt>=600) return 0.710885 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.189159 ;
	  if(pt>=50 && pt<70) return 0.155906 ;
	  if(pt>=70 && pt<100) return 0.147914 ;
	  if(pt>=100 && pt<140) return 0.146295 ;
	  if(pt>=140 && pt<200) return 0.15453 ;
	  if(pt>=200 && pt<300) return 0.178051 ;
	  if(pt>=300 && pt<600) return 0.22845 ;
	  if(pt>=600) return 0.266209 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.187647 ;
	  if(pt>=50 && pt<70) return 0.162878 ;
	  if(pt>=70 && pt<100) return 0.152739 ;
	  if(pt>=100 && pt<140) return 0.15007 ;
	  if(pt>=140 && pt<200) return 0.159197 ;
	  if(pt>=200 && pt<300) return 0.178751 ;
	  if(pt>=300 && pt<600) return 0.225035 ;
	  if(pt>=600) return 0.254068 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.184185 ;
	  if(pt>=50 && pt<70) return 0.158081 ;
	  if(pt>=70 && pt<100) return 0.153773 ;
	  if(pt>=100 && pt<140) return 0.151556 ;
	  if(pt>=140 && pt<200) return 0.161687 ;
	  if(pt>=200 && pt<300) return 0.175639 ;
	  if(pt>=300 && pt<600) return 0.216642 ;
	  if(pt>=600) return 0.249795 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.171865 ;
	  if(pt>=50 && pt<70) return 0.145115 ;
	  if(pt>=70 && pt<100) return 0.13597 ;
	  if(pt>=100 && pt<140) return 0.133302 ;
	  if(pt>=140 && pt<200) return 0.144276 ;
	  if(pt>=200 && pt<300) return 0.153882 ;
	  if(pt>=300 && pt<600) return 0.186202 ;
	  if(pt>=600) return 0.275906 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0248063 ;
	  if(pt>=50 && pt<70) return 0.0132907 ;
	  if(pt>=70 && pt<100) return 0.0105736 ;
	  if(pt>=100 && pt<140) return 0.00978675 ;
	  if(pt>=140 && pt<200) return 0.00966648 ;
	  if(pt>=200 && pt<300) return 0.0116475 ;
	  if(pt>=300 && pt<600) return 0.0224327 ;
	  if(pt>=600) return 0.053819 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0260555 ;
	  if(pt>=50 && pt<70) return 0.0148537 ;
	  if(pt>=70 && pt<100) return 0.0116983 ;
	  if(pt>=100 && pt<140) return 0.0108652 ;
	  if(pt>=140 && pt<200) return 0.0108113 ;
	  if(pt>=200 && pt<300) return 0.0119786 ;
	  if(pt>=300 && pt<600) return 0.0216516 ;
	  if(pt>=600) return 0.0545822 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0278645 ;
	  if(pt>=50 && pt<70) return 0.0170498 ;
	  if(pt>=70 && pt<100) return 0.0141185 ;
	  if(pt>=100 && pt<140) return 0.0138236 ;
	  if(pt>=140 && pt<200) return 0.0147357 ;
	  if(pt>=200 && pt<300) return 0.0171958 ;
	  if(pt>=300 && pt<600) return 0.0326425 ;
	  if(pt>=600) return 0.0842763 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0334985 ;
	  if(pt>=50 && pt<70) return 0.0204545 ;
	  if(pt>=70 && pt<100) return 0.0172483 ;
	  if(pt>=100 && pt<140) return 0.0172067 ;
	  if(pt>=140 && pt<200) return 0.0196076 ;
	  if(pt>=200 && pt<300) return 0.0243433 ;
	  if(pt>=300 && pt<600) return 0.0420726 ;
	  if(pt>=600) return 0.0926417 ;
	}
      }
    }
  else if(sample == "tt")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.816965 ;
	  if(pt>=50 && pt<70) return 0.83737 ;
	  if(pt>=70 && pt<100) return 0.846435 ;
	  if(pt>=100 && pt<140) return 0.85388 ;
	  if(pt>=140 && pt<200) return 0.863017 ;
	  if(pt>=200 && pt<300) return 0.864682 ;
	  if(pt>=300 && pt<600) return 0.84824 ;
	  if(pt>=600) return 0.796696 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.803635 ;
	  if(pt>=50 && pt<70) return 0.825482 ;
	  if(pt>=70 && pt<100) return 0.835112 ;
	  if(pt>=100 && pt<140) return 0.842853 ;
	  if(pt>=140 && pt<200) return 0.852502 ;
	  if(pt>=200 && pt<300) return 0.853063 ;
	  if(pt>=300 && pt<600) return 0.832109 ;
	  if(pt>=600) return 0.770352 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.75992 ;
	  if(pt>=50 && pt<70) return 0.78299 ;
	  if(pt>=70 && pt<100) return 0.792012 ;
	  if(pt>=100 && pt<140) return 0.799117 ;
	  if(pt>=140 && pt<200) return 0.810987 ;
	  if(pt>=200 && pt<300) return 0.810666 ;
	  if(pt>=300 && pt<600) return 0.778932 ;
	  if(pt>=600) return 0.701484 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.724582 ;
	  if(pt>=50 && pt<70) return 0.747582 ;
	  if(pt>=70 && pt<100) return 0.756247 ;
	  if(pt>=100 && pt<140) return 0.766018 ;
	  if(pt>=140 && pt<200) return 0.777683 ;
	  if(pt>=200 && pt<300) return 0.776282 ;
	  if(pt>=300 && pt<600) return 0.742312 ;
	  if(pt>=600) return 0.68142 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.189897 ;
	  if(pt>=50 && pt<70) return 0.158656 ;
	  if(pt>=70 && pt<100) return 0.149459 ;
	  if(pt>=100 && pt<140) return 0.149078 ;
	  if(pt>=140 && pt<200) return 0.158698 ;
	  if(pt>=200 && pt<300) return 0.183899 ;
	  if(pt>=300 && pt<600) return 0.233864 ;
	  if(pt>=600) return 0.269221 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.190655 ;
	  if(pt>=50 && pt<70) return 0.162364 ;
	  if(pt>=70 && pt<100) return 0.154252 ;
	  if(pt>=100 && pt<140) return 0.15406 ;
	  if(pt>=140 && pt<200) return 0.162561 ;
	  if(pt>=200 && pt<300) return 0.183967 ;
	  if(pt>=300 && pt<600) return 0.227488 ;
	  if(pt>=600) return 0.255044 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.185327 ;
	  if(pt>=50 && pt<70) return 0.160725 ;
	  if(pt>=70 && pt<100) return 0.153015 ;
	  if(pt>=100 && pt<140) return 0.153086 ;
	  if(pt>=140 && pt<200) return 0.163389 ;
	  if(pt>=200 && pt<300) return 0.181199 ;
	  if(pt>=300 && pt<600) return 0.219281 ;
	  if(pt>=600) return 0.247 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.177768 ;
	  if(pt>=50 && pt<70) return 0.151661 ;
	  if(pt>=70 && pt<100) return 0.138752 ;
	  if(pt>=100 && pt<140) return 0.139292 ;
	  if(pt>=140 && pt<200) return 0.149254 ;
	  if(pt>=200 && pt<300) return 0.167331 ;
	  if(pt>=300 && pt<600) return 0.201055 ;
	  if(pt>=600) return 0.231121 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0189149 ;
	  if(pt>=50 && pt<70) return 0.0105821 ;
	  if(pt>=70 && pt<100) return 0.00874925 ;
	  if(pt>=100 && pt<140) return 0.00809742 ;
	  if(pt>=140 && pt<200) return 0.00855031 ;
	  if(pt>=200 && pt<300) return 0.0104137 ;
	  if(pt>=300 && pt<600) return 0.019268 ;
	  if(pt>=600) return 0.0552555 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0204881 ;
	  if(pt>=50 && pt<70) return 0.0119297 ;
	  if(pt>=70 && pt<100) return 0.00993381 ;
	  if(pt>=100 && pt<140) return 0.00925907 ;
	  if(pt>=140 && pt<200) return 0.00952704 ;
	  if(pt>=200 && pt<300) return 0.0111743 ;
	  if(pt>=300 && pt<600) return 0.0201289 ;
	  if(pt>=600) return 0.0562302 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0232926 ;
	  if(pt>=50 && pt<70) return 0.0143408 ;
	  if(pt>=70 && pt<100) return 0.0122606 ;
	  if(pt>=100 && pt<140) return 0.011937 ;
	  if(pt>=140 && pt<200) return 0.0131115 ;
	  if(pt>=200 && pt<300) return 0.0163477 ;
	  if(pt>=300 && pt<600) return 0.0324909 ;
	  if(pt>=600) return 0.0892736 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0306497 ;
	  if(pt>=50 && pt<70) return 0.0187395 ;
	  if(pt>=70 && pt<100) return 0.0157941 ;
	  if(pt>=100 && pt<140) return 0.0159116 ;
	  if(pt>=140 && pt<200) return 0.0184096 ;
	  if(pt>=200 && pt<300) return 0.0244386 ;
	  if(pt>=300 && pt<600) return 0.0451509 ;
	  if(pt>=600) return 0.0950604 ;
	}
      }
    }
  else if(sample == "dib")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.805297 ;
	  if(pt>=50 && pt<70) return 0.822044 ;
	  if(pt>=70 && pt<100) return 0.838976 ;
	  if(pt>=100 && pt<140) return 0.838412 ;
	  if(pt>=140 && pt<200) return 0.849905 ;
	  if(pt>=200 && pt<300) return 0.845748 ;
	  if(pt>=300 && pt<600) return 0.860119 ;
	  if(pt>=600) return 0.893207 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.784144 ;
	  if(pt>=50 && pt<70) return 0.820133 ;
	  if(pt>=70 && pt<100) return 0.833056 ;
	  if(pt>=100 && pt<140) return 0.84032 ;
	  if(pt>=140 && pt<200) return 0.831514 ;
	  if(pt>=200 && pt<300) return 0.837511 ;
	  if(pt>=300 && pt<600) return 0.853525 ;
	  if(pt>=600) return 0.884336 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.758221 ;
	  if(pt>=50 && pt<70) return 0.783229 ;
	  if(pt>=70 && pt<100) return 0.785464 ;
	  if(pt>=100 && pt<140) return 0.802451 ;
	  if(pt>=140 && pt<200) return 0.800113 ;
	  if(pt>=200 && pt<300) return 0.807111 ;
	  if(pt>=300 && pt<600) return 0.816639 ;
	  if(pt>=600) return 0.760659 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.706545 ;
	  if(pt>=50 && pt<70) return 0.74044 ;
	  if(pt>=70 && pt<100) return 0.749949 ;
	  if(pt>=100 && pt<140) return 0.756609 ;
	  if(pt>=140 && pt<200) return 0.757977 ;
	  if(pt>=200 && pt<300) return 0.774698 ;
	  if(pt>=300 && pt<600) return 0.76279 ;
	  if(pt>=600) return 1 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.186279 ;
	  if(pt>=50 && pt<70) return 0.15766 ;
	  if(pt>=70 && pt<100) return 0.147145 ;
	  if(pt>=100 && pt<140) return 0.14554 ;
	  if(pt>=140 && pt<200) return 0.157733 ;
	  if(pt>=200 && pt<300) return 0.178447 ;
	  if(pt>=300 && pt<600) return 0.241236 ;
	  if(pt>=600) return 0.26208 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.195977 ;
	  if(pt>=50 && pt<70) return 0.159022 ;
	  if(pt>=70 && pt<100) return 0.151617 ;
	  if(pt>=100 && pt<140) return 0.14866 ;
	  if(pt>=140 && pt<200) return 0.156054 ;
	  if(pt>=200 && pt<300) return 0.178598 ;
	  if(pt>=300 && pt<600) return 0.236309 ;
	  if(pt>=600) return 0.284858 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.184481 ;
	  if(pt>=50 && pt<70) return 0.158029 ;
	  if(pt>=70 && pt<100) return 0.151623 ;
	  if(pt>=100 && pt<140) return 0.152472 ;
	  if(pt>=140 && pt<200) return 0.158007 ;
	  if(pt>=200 && pt<300) return 0.184024 ;
	  if(pt>=300 && pt<600) return 0.228013 ;
	  if(pt>=600) return 0.287827 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.175861 ;
	  if(pt>=50 && pt<70) return 0.155588 ;
	  if(pt>=70 && pt<100) return 0.138591 ;
	  if(pt>=100 && pt<140) return 0.132658 ;
	  if(pt>=140 && pt<200) return 0.153439 ;
	  if(pt>=200 && pt<300) return 0.169866 ;
	  if(pt>=300 && pt<600) return 0.215919 ;
	  if(pt>=600) return 0.309754 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0197615 ;
	  if(pt>=50 && pt<70) return 0.010124 ;
	  if(pt>=70 && pt<100) return 0.00880611 ;
	  if(pt>=100 && pt<140) return 0.00835179 ;
	  if(pt>=140 && pt<200) return 0.00849048 ;
	  if(pt>=200 && pt<300) return 0.0106437 ;
	  if(pt>=300 && pt<600) return 0.020581 ;
	  if(pt>=600) return 0.0556525 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0222682 ;
	  if(pt>=50 && pt<70) return 0.0122596 ;
	  if(pt>=70 && pt<100) return 0.00984396 ;
	  if(pt>=100 && pt<140) return 0.00948045 ;
	  if(pt>=140 && pt<200) return 0.00950385 ;
	  if(pt>=200 && pt<300) return 0.0115501 ;
	  if(pt>=300 && pt<600) return 0.0227227 ;
	  if(pt>=600) return 0.0535392 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0253517 ;
	  if(pt>=50 && pt<70) return 0.0151463 ;
	  if(pt>=70 && pt<100) return 0.0133791 ;
	  if(pt>=100 && pt<140) return 0.0133582 ;
	  if(pt>=140 && pt<200) return 0.0140878 ;
	  if(pt>=200 && pt<300) return 0.0162538 ;
	  if(pt>=300 && pt<600) return 0.0328523 ;
	  if(pt>=600) return 0.0742136 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0343988 ;
	  if(pt>=50 && pt<70) return 0.0198165 ;
	  if(pt>=70 && pt<100) return 0.0180723 ;
	  if(pt>=100 && pt<140) return 0.0163887 ;
	  if(pt>=140 && pt<200) return 0.0188143 ;
	  if(pt>=200 && pt<300) return 0.0224715 ;
	  if(pt>=300 && pt<600) return 0.040918 ;
	  if(pt>=600) return 0.0851774 ;
	}
      }
    }
  return 1.0;
}

double BTag_MCEfficiency_T(string sample, int flavor, double pt, double eta){
  if(sample == "dy")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.609562 ;
	  if(pt>=50 && pt<70) return 0.66886 ;
	  if(pt>=70 && pt<100) return 0.668088 ;
	  if(pt>=100 && pt<140) return 0.684057 ;
	  if(pt>=140 && pt<200) return 0.708561 ;
	  if(pt>=200 && pt<300) return 0.720736 ;
	  if(pt>=300 && pt<600) return 0.688987 ;
	  if(pt>=600) return 0.613929 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.579457 ;
	  if(pt>=50 && pt<70) return 0.63145 ;
	  if(pt>=70 && pt<100) return 0.674152 ;
	  if(pt>=100 && pt<140) return 0.676845 ;
	  if(pt>=140 && pt<200) return 0.702885 ;
	  if(pt>=200 && pt<300) return 0.68646 ;
	  if(pt>=300 && pt<600) return 0.675687 ;
	  if(pt>=600) return 0.616382 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.517482 ;
	  if(pt>=50 && pt<70) return 0.553115 ;
	  if(pt>=70 && pt<100) return 0.57192 ;
	  if(pt>=100 && pt<140) return 0.613745 ;
	  if(pt>=140 && pt<200) return 0.607951 ;
	  if(pt>=200 && pt<300) return 0.626562 ;
	  if(pt>=300 && pt<600) return 0.571751 ;
	  if(pt>=600) return 0.491298 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.484882 ;
	  if(pt>=50 && pt<70) return 0.512079 ;
	  if(pt>=70 && pt<100) return 0.547112 ;
	  if(pt>=100 && pt<140) return 0.509671 ;
	  if(pt>=140 && pt<200) return 0.533947 ;
	  if(pt>=200 && pt<300) return 0.583682 ;
	  if(pt>=300 && pt<600) return 0.557144 ;
	  if(pt>=600) return 0.430265 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0385313 ;
	  if(pt>=50 && pt<70) return 0.0406435 ;
	  if(pt>=70 && pt<100) return 0.0352714 ;
	  if(pt>=100 && pt<140) return 0.0344298 ;
	  if(pt>=140 && pt<200) return 0.0353996 ;
	  if(pt>=200 && pt<300) return 0.0379238 ;
	  if(pt>=300 && pt<600) return 0.0545672 ;
	  if(pt>=600) return 0.0599737 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0337207 ;
	  if(pt>=50 && pt<70) return 0.0298835 ;
	  if(pt>=70 && pt<100) return 0.0322891 ;
	  if(pt>=100 && pt<140) return 0.0333392 ;
	  if(pt>=140 && pt<200) return 0.0350568 ;
	  if(pt>=200 && pt<300) return 0.0409501 ;
	  if(pt>=300 && pt<600) return 0.0531835 ;
	  if(pt>=600) return 0.0776586 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0290213 ;
	  if(pt>=50 && pt<70) return 0.0298834 ;
	  if(pt>=70 && pt<100) return 0.0310321 ;
	  if(pt>=100 && pt<140) return 0.0346433 ;
	  if(pt>=140 && pt<200) return 0.0357677 ;
	  if(pt>=200 && pt<300) return 0.0394672 ;
	  if(pt>=300 && pt<600) return 0.0481731 ;
	  if(pt>=600) return 0.0491902 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0277545 ;
	  if(pt>=50 && pt<70) return 0.0337877 ;
	  if(pt>=70 && pt<100) return 0.0288879 ;
	  if(pt>=100 && pt<140) return 0.0295434 ;
	  if(pt>=140 && pt<200) return 0.0268315 ;
	  if(pt>=200 && pt<300) return 0.0378055 ;
	  if(pt>=300 && pt<600) return 0.0455196 ;
	  if(pt>=600) return 0.0406405 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00435117 ;
	  if(pt>=50 && pt<70) return 0.00245451 ;
	  if(pt>=70 && pt<100) return 0.00282819 ;
	  if(pt>=100 && pt<140) return 0.00308405 ;
	  if(pt>=140 && pt<200) return 0.00319033 ;
	  if(pt>=200 && pt<300) return 0.00369907 ;
	  if(pt>=300 && pt<600) return 0.00255475 ;
	  if(pt>=600) return 0.00868971 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00345712 ;
	  if(pt>=50 && pt<70) return 0.00359923 ;
	  if(pt>=70 && pt<100) return 0.00327322 ;
	  if(pt>=100 && pt<140) return 0.00262515 ;
	  if(pt>=140 && pt<200) return 0.00370966 ;
	  if(pt>=200 && pt<300) return 0.00263894 ;
	  if(pt>=300 && pt<600) return 0.00303871 ;
	  if(pt>=600) return 0.00642996 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00459408 ;
	  if(pt>=50 && pt<70) return 0.00351647 ;
	  if(pt>=70 && pt<100) return 0.0047584 ;
	  if(pt>=100 && pt<140) return 0.00424259 ;
	  if(pt>=140 && pt<200) return 0.00344657 ;
	  if(pt>=200 && pt<300) return 0.00323033 ;
	  if(pt>=300 && pt<600) return 0.00447821 ;
	  if(pt>=600) return 0.0142239 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0051219 ;
	  if(pt>=50 && pt<70) return 0.00486096 ;
	  if(pt>=70 && pt<100) return 0.00335606 ;
	  if(pt>=100 && pt<140) return 0.0040784 ;
	  if(pt>=140 && pt<200) return 0.00376852 ;
	  if(pt>=200 && pt<300) return 0.00500413 ;
	  if(pt>=300 && pt<600) return 0.0051441 ;
	  if(pt>=600) return 0.00914939 ;
	}
      }
    }
  else if(sample == "qcd")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.606993 ;
	  if(pt>=50 && pt<70) return 0.657891 ;
	  if(pt>=70 && pt<100) return 0.680853 ;
	  if(pt>=100 && pt<140) return 0.693134 ;
	  if(pt>=140 && pt<200) return 0.710728 ;
	  if(pt>=200 && pt<300) return 0.728514 ;
	  if(pt>=300 && pt<600) return 0.713773 ;
	  if(pt>=600) return 0.648692 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.588006 ;
	  if(pt>=50 && pt<70) return 0.639931 ;
	  if(pt>=70 && pt<100) return 0.661624 ;
	  if(pt>=100 && pt<140) return 0.673706 ;
	  if(pt>=140 && pt<200) return 0.691566 ;
	  if(pt>=200 && pt<300) return 0.70943 ;
	  if(pt>=300 && pt<600) return 0.686956 ;
	  if(pt>=600) return 0.608022 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.519542 ;
	  if(pt>=50 && pt<70) return 0.574792 ;
	  if(pt>=70 && pt<100) return 0.597511 ;
	  if(pt>=100 && pt<140) return 0.608711 ;
	  if(pt>=140 && pt<200) return 0.626355 ;
	  if(pt>=200 && pt<300) return 0.641326 ;
	  if(pt>=300 && pt<600) return 0.604793 ;
	  if(pt>=600) return 0.50543 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.482935 ;
	  if(pt>=50 && pt<70) return 0.530467 ;
	  if(pt>=70 && pt<100) return 0.544468 ;
	  if(pt>=100 && pt<140) return 0.560489 ;
	  if(pt>=140 && pt<200) return 0.567674 ;
	  if(pt>=200 && pt<300) return 0.576709 ;
	  if(pt>=300 && pt<600) return 0.526687 ;
	  if(pt>=600) return 0.462763 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0354337 ;
	  if(pt>=50 && pt<70) return 0.0352486 ;
	  if(pt>=70 && pt<100) return 0.0366913 ;
	  if(pt>=100 && pt<140) return 0.0391796 ;
	  if(pt>=140 && pt<200) return 0.0434168 ;
	  if(pt>=200 && pt<300) return 0.0491295 ;
	  if(pt>=300 && pt<600) return 0.0601161 ;
	  if(pt>=600) return 0.0686708 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0333746 ;
	  if(pt>=50 && pt<70) return 0.0340465 ;
	  if(pt>=70 && pt<100) return 0.0368204 ;
	  if(pt>=100 && pt<140) return 0.0388721 ;
	  if(pt>=140 && pt<200) return 0.0428865 ;
	  if(pt>=200 && pt<300) return 0.0483083 ;
	  if(pt>=300 && pt<600) return 0.0569545 ;
	  if(pt>=600) return 0.0616901 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0322877 ;
	  if(pt>=50 && pt<70) return 0.0320069 ;
	  if(pt>=70 && pt<100) return 0.0329509 ;
	  if(pt>=100 && pt<140) return 0.0348534 ;
	  if(pt>=140 && pt<200) return 0.0396512 ;
	  if(pt>=200 && pt<300) return 0.0434156 ;
	  if(pt>=300 && pt<600) return 0.050365 ;
	  if(pt>=600) return 0.05261 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0302536 ;
	  if(pt>=50 && pt<70) return 0.0283691 ;
	  if(pt>=70 && pt<100) return 0.0283851 ;
	  if(pt>=100 && pt<140) return 0.0318579 ;
	  if(pt>=140 && pt<200) return 0.0353451 ;
	  if(pt>=200 && pt<300) return 0.038864 ;
	  if(pt>=300 && pt<600) return 0.0435706 ;
	  if(pt>=600) return 0.047808 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00225005 ;
	  if(pt>=50 && pt<70) return 0.00142215 ;
	  if(pt>=70 && pt<100) return 0.00132762 ;
	  if(pt>=100 && pt<140) return 0.00112203 ;
	  if(pt>=140 && pt<200) return 0.00104455 ;
	  if(pt>=200 && pt<300) return 0.00148094 ;
	  if(pt>=300 && pt<600) return 0.0017277 ;
	  if(pt>=600) return 0.00631901 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00195784 ;
	  if(pt>=50 && pt<70) return 0.00141935 ;
	  if(pt>=70 && pt<100) return 0.00115079 ;
	  if(pt>=100 && pt<140) return 0.00119895 ;
	  if(pt>=140 && pt<200) return 0.00128774 ;
	  if(pt>=200 && pt<300) return 0.00130869 ;
	  if(pt>=300 && pt<600) return 0.00184315 ;
	  if(pt>=600) return 0.00614823 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00226444 ;
	  if(pt>=50 && pt<70) return 0.00188788 ;
	  if(pt>=70 && pt<100) return 0.00221784 ;
	  if(pt>=100 && pt<140) return 0.00172192 ;
	  if(pt>=140 && pt<200) return 0.00157725 ;
	  if(pt>=200 && pt<300) return 0.00182911 ;
	  if(pt>=300 && pt<600) return 0.00344651 ;
	  if(pt>=600) return 0.0103198 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0030538 ;
	  if(pt>=50 && pt<70) return 0.00241809 ;
	  if(pt>=70 && pt<100) return 0.00239224 ;
	  if(pt>=100 && pt<140) return 0.0024095 ;
	  if(pt>=140 && pt<200) return 0.00248518 ;
	  if(pt>=200 && pt<300) return 0.00326542 ;
	  if(pt>=300 && pt<600) return 0.00520796 ;
	  if(pt>=600) return 0.00956375 ;
	}
      }
    }
  else if(sample == "wj")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.603359 ;
	  if(pt>=50 && pt<70) return 0.652562 ;
	  if(pt>=70 && pt<100) return 0.669954 ;
	  if(pt>=100 && pt<140) return 0.677417 ;
	  if(pt>=140 && pt<200) return 0.677102 ;
	  if(pt>=200 && pt<300) return 0.687615 ;
	  if(pt>=300 && pt<600) return 0.681069 ;
	  if(pt>=600) return 0.633039 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.582496 ;
	  if(pt>=50 && pt<70) return 0.632129 ;
	  if(pt>=70 && pt<100) return 0.654626 ;
	  if(pt>=100 && pt<140) return 0.655224 ;
	  if(pt>=140 && pt<200) return 0.659195 ;
	  if(pt>=200 && pt<300) return 0.665925 ;
	  if(pt>=300 && pt<600) return 0.653157 ;
	  if(pt>=600) return 0.59832 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.518018 ;
	  if(pt>=50 && pt<70) return 0.575565 ;
	  if(pt>=70 && pt<100) return 0.590957 ;
	  if(pt>=100 && pt<140) return 0.591573 ;
	  if(pt>=140 && pt<200) return 0.591465 ;
	  if(pt>=200 && pt<300) return 0.586487 ;
	  if(pt>=300 && pt<600) return 0.567062 ;
	  if(pt>=600) return 0.499226 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.482662 ;
	  if(pt>=50 && pt<70) return 0.524842 ;
	  if(pt>=70 && pt<100) return 0.537956 ;
	  if(pt>=100 && pt<140) return 0.544811 ;
	  if(pt>=140 && pt<200) return 0.532495 ;
	  if(pt>=200 && pt<300) return 0.527679 ;
	  if(pt>=300 && pt<600) return 0.507942 ;
	  if(pt>=600) return 0.445194 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0342263 ;
	  if(pt>=50 && pt<70) return 0.0341014 ;
	  if(pt>=70 && pt<100) return 0.034643 ;
	  if(pt>=100 && pt<140) return 0.0329682 ;
	  if(pt>=140 && pt<200) return 0.0361817 ;
	  if(pt>=200 && pt<300) return 0.0425708 ;
	  if(pt>=300 && pt<600) return 0.0552784 ;
	  if(pt>=600) return 0.0657497 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0329687 ;
	  if(pt>=50 && pt<70) return 0.0340773 ;
	  if(pt>=70 && pt<100) return 0.0342733 ;
	  if(pt>=100 && pt<140) return 0.0339737 ;
	  if(pt>=140 && pt<200) return 0.035726 ;
	  if(pt>=200 && pt<300) return 0.0419052 ;
	  if(pt>=300 && pt<600) return 0.0534383 ;
	  if(pt>=600) return 0.0590889 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0313398 ;
	  if(pt>=50 && pt<70) return 0.0298728 ;
	  if(pt>=70 && pt<100) return 0.030491 ;
	  if(pt>=100 && pt<140) return 0.0323159 ;
	  if(pt>=140 && pt<200) return 0.033676 ;
	  if(pt>=200 && pt<300) return 0.0399848 ;
	  if(pt>=300 && pt<600) return 0.0479233 ;
	  if(pt>=600) return 0.0512813 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0289313 ;
	  if(pt>=50 && pt<70) return 0.0267585 ;
	  if(pt>=70 && pt<100) return 0.0262376 ;
	  if(pt>=100 && pt<140) return 0.028486 ;
	  if(pt>=140 && pt<200) return 0.030746 ;
	  if(pt>=200 && pt<300) return 0.0357853 ;
	  if(pt>=300 && pt<600) return 0.0405556 ;
	  if(pt>=600) return 0.0442072 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0025517 ;
	  if(pt>=50 && pt<70) return 0.00133605 ;
	  if(pt>=70 && pt<100) return 0.00111657 ;
	  if(pt>=100 && pt<140) return 0.0010366 ;
	  if(pt>=140 && pt<200) return 0.000993665 ;
	  if(pt>=200 && pt<300) return 0.00127072 ;
	  if(pt>=300 && pt<600) return 0.0024696 ;
	  if(pt>=600) return 0.00746634 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00269828 ;
	  if(pt>=50 && pt<70) return 0.00157217 ;
	  if(pt>=70 && pt<100) return 0.00128536 ;
	  if(pt>=100 && pt<140) return 0.00125278 ;
	  if(pt>=140 && pt<200) return 0.00128973 ;
	  if(pt>=200 && pt<300) return 0.00143229 ;
	  if(pt>=300 && pt<600) return 0.00252651 ;
	  if(pt>=600) return 0.00685351 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00289002 ;
	  if(pt>=50 && pt<70) return 0.00184071 ;
	  if(pt>=70 && pt<100) return 0.00165596 ;
	  if(pt>=100 && pt<140) return 0.00175765 ;
	  if(pt>=140 && pt<200) return 0.00189095 ;
	  if(pt>=200 && pt<300) return 0.00226578 ;
	  if(pt>=300 && pt<600) return 0.00416797 ;
	  if(pt>=600) return 0.010486 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00381654 ;
	  if(pt>=50 && pt<70) return 0.00246661 ;
	  if(pt>=70 && pt<100) return 0.00220066 ;
	  if(pt>=100 && pt<140) return 0.00254077 ;
	  if(pt>=140 && pt<200) return 0.00285033 ;
	  if(pt>=200 && pt<300) return 0.00350716 ;
	  if(pt>=300 && pt<600) return 0.00563313 ;
	  if(pt>=600) return 0.0105518 ;
	}
      }
    }
  else if(sample == "ttx")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.627654 ;
	  if(pt>=50 && pt<70) return 0.671128 ;
	  if(pt>=70 && pt<100) return 0.689437 ;
	  if(pt>=100 && pt<140) return 0.701529 ;
	  if(pt>=140 && pt<200) return 0.712395 ;
	  if(pt>=200 && pt<300) return 0.712611 ;
	  if(pt>=300 && pt<600) return 0.687257 ;
	  if(pt>=600) return 0.575727 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.605575 ;
	  if(pt>=50 && pt<70) return 0.649829 ;
	  if(pt>=70 && pt<100) return 0.67044 ;
	  if(pt>=100 && pt<140) return 0.683569 ;
	  if(pt>=140 && pt<200) return 0.692684 ;
	  if(pt>=200 && pt<300) return 0.693249 ;
	  if(pt>=300 && pt<600) return 0.658793 ;
	  if(pt>=600) return 0.546297 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.540473 ;
	  if(pt>=50 && pt<70) return 0.583105 ;
	  if(pt>=70 && pt<100) return 0.602155 ;
	  if(pt>=100 && pt<140) return 0.614333 ;
	  if(pt>=140 && pt<200) return 0.625108 ;
	  if(pt>=200 && pt<300) return 0.623971 ;
	  if(pt>=300 && pt<600) return 0.569281 ;
	  if(pt>=600) return 0.45673 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.497495 ;
	  if(pt>=50 && pt<70) return 0.540028 ;
	  if(pt>=70 && pt<100) return 0.554646 ;
	  if(pt>=100 && pt<140) return 0.57635 ;
	  if(pt>=140 && pt<200) return 0.584833 ;
	  if(pt>=200 && pt<300) return 0.580576 ;
	  if(pt>=300 && pt<600) return 0.536363 ;
	  if(pt>=600) return 0.425351 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0306986 ;
	  if(pt>=50 && pt<70) return 0.0253632 ;
	  if(pt>=70 && pt<100) return 0.0259264 ;
	  if(pt>=100 && pt<140) return 0.0269098 ;
	  if(pt>=140 && pt<200) return 0.0321254 ;
	  if(pt>=200 && pt<300) return 0.038077 ;
	  if(pt>=300 && pt<600) return 0.0576638 ;
	  if(pt>=600) return 0.0736834 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0304927 ;
	  if(pt>=50 && pt<70) return 0.0258105 ;
	  if(pt>=70 && pt<100) return 0.026651 ;
	  if(pt>=100 && pt<140) return 0.0291099 ;
	  if(pt>=140 && pt<200) return 0.0319972 ;
	  if(pt>=200 && pt<300) return 0.0391066 ;
	  if(pt>=300 && pt<600) return 0.0533452 ;
	  if(pt>=600) return 0.0674263 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0292608 ;
	  if(pt>=50 && pt<70) return 0.0259198 ;
	  if(pt>=70 && pt<100) return 0.0256836 ;
	  if(pt>=100 && pt<140) return 0.0271594 ;
	  if(pt>=140 && pt<200) return 0.0314313 ;
	  if(pt>=200 && pt<300) return 0.037767 ;
	  if(pt>=300 && pt<600) return 0.0526074 ;
	  if(pt>=600) return 0.0723687 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0294207 ;
	  if(pt>=50 && pt<70) return 0.024925 ;
	  if(pt>=70 && pt<100) return 0.0235528 ;
	  if(pt>=100 && pt<140) return 0.0251182 ;
	  if(pt>=140 && pt<200) return 0.0299068 ;
	  if(pt>=200 && pt<300) return 0.035304 ;
	  if(pt>=300 && pt<600) return 0.0406704 ;
	  if(pt>=600) return 0.0352834 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00208799 ;
	  if(pt>=50 && pt<70) return 0.00109447 ;
	  if(pt>=70 && pt<100) return 0.00107797 ;
	  if(pt>=100 && pt<140) return 0.00102595 ;
	  if(pt>=140 && pt<200) return 0.00135178 ;
	  if(pt>=200 && pt<300) return 0.00172229 ;
	  if(pt>=300 && pt<600) return 0.00305026 ;
	  if(pt>=600) return 0.00760301 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00227443 ;
	  if(pt>=50 && pt<70) return 0.00127177 ;
	  if(pt>=70 && pt<100) return 0.00116439 ;
	  if(pt>=100 && pt<140) return 0.00140324 ;
	  if(pt>=140 && pt<200) return 0.00138572 ;
	  if(pt>=200 && pt<300) return 0.00190477 ;
	  if(pt>=300 && pt<600) return 0.0034383 ;
	  if(pt>=600) return 0.00593384 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00229102 ;
	  if(pt>=50 && pt<70) return 0.00139779 ;
	  if(pt>=70 && pt<100) return 0.00132386 ;
	  if(pt>=100 && pt<140) return 0.00160431 ;
	  if(pt>=140 && pt<200) return 0.00185331 ;
	  if(pt>=200 && pt<300) return 0.00255098 ;
	  if(pt>=300 && pt<600) return 0.00451145 ;
	  if(pt>=600) return 0.0111411 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00264198 ;
	  if(pt>=50 && pt<70) return 0.00184122 ;
	  if(pt>=70 && pt<100) return 0.00181692 ;
	  if(pt>=100 && pt<140) return 0.00219547 ;
	  if(pt>=140 && pt<200) return 0.002465 ;
	  if(pt>=200 && pt<300) return 0.00344826 ;
	  if(pt>=300 && pt<600) return 0.00641116 ;
	  if(pt>=600) return 0.0147716 ;
	}
      }
    }
  else if(sample == "st")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.633743 ;
	  if(pt>=50 && pt<70) return 0.678543 ;
	  if(pt>=70 && pt<100) return 0.698648 ;
	  if(pt>=100 && pt<140) return 0.713832 ;
	  if(pt>=140 && pt<200) return 0.729019 ;
	  if(pt>=200 && pt<300) return 0.734422 ;
	  if(pt>=300 && pt<600) return 0.71115 ;
	  if(pt>=600) return 0.588091 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.612575 ;
	  if(pt>=50 && pt<70) return 0.658401 ;
	  if(pt>=70 && pt<100) return 0.679443 ;
	  if(pt>=100 && pt<140) return 0.69581 ;
	  if(pt>=140 && pt<200) return 0.712492 ;
	  if(pt>=200 && pt<300) return 0.714435 ;
	  if(pt>=300 && pt<600) return 0.683615 ;
	  if(pt>=600) return 0.540342 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.548997 ;
	  if(pt>=50 && pt<70) return 0.589218 ;
	  if(pt>=70 && pt<100) return 0.608676 ;
	  if(pt>=100 && pt<140) return 0.623129 ;
	  if(pt>=140 && pt<200) return 0.641096 ;
	  if(pt>=200 && pt<300) return 0.642261 ;
	  if(pt>=300 && pt<600) return 0.592382 ;
	  if(pt>=600) return 0.417977 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.503975 ;
	  if(pt>=50 && pt<70) return 0.544555 ;
	  if(pt>=70 && pt<100) return 0.562891 ;
	  if(pt>=100 && pt<140) return 0.579535 ;
	  if(pt>=140 && pt<200) return 0.597457 ;
	  if(pt>=200 && pt<300) return 0.59377 ;
	  if(pt>=300 && pt<600) return 0.556541 ;
	  if(pt>=600) return 0.43632 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0308233 ;
	  if(pt>=50 && pt<70) return 0.0261362 ;
	  if(pt>=70 && pt<100) return 0.0262731 ;
	  if(pt>=100 && pt<140) return 0.0271796 ;
	  if(pt>=140 && pt<200) return 0.030467 ;
	  if(pt>=200 && pt<300) return 0.0369805 ;
	  if(pt>=300 && pt<600) return 0.0502294 ;
	  if(pt>=600) return 0.0539919 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0316377 ;
	  if(pt>=50 && pt<70) return 0.0272661 ;
	  if(pt>=70 && pt<100) return 0.0270438 ;
	  if(pt>=100 && pt<140) return 0.0277009 ;
	  if(pt>=140 && pt<200) return 0.0323295 ;
	  if(pt>=200 && pt<300) return 0.0385708 ;
	  if(pt>=300 && pt<600) return 0.0487467 ;
	  if(pt>=600) return 0.0461867 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0302402 ;
	  if(pt>=50 && pt<70) return 0.0255824 ;
	  if(pt>=70 && pt<100) return 0.0266356 ;
	  if(pt>=100 && pt<140) return 0.0276508 ;
	  if(pt>=140 && pt<200) return 0.032506 ;
	  if(pt>=200 && pt<300) return 0.0375771 ;
	  if(pt>=300 && pt<600) return 0.047264 ;
	  if(pt>=600) return 0.0609234 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0286836 ;
	  if(pt>=50 && pt<70) return 0.0239411 ;
	  if(pt>=70 && pt<100) return 0.0246508 ;
	  if(pt>=100 && pt<140) return 0.0249666 ;
	  if(pt>=140 && pt<200) return 0.0290569 ;
	  if(pt>=200 && pt<300) return 0.0331478 ;
	  if(pt>=300 && pt<600) return 0.0408733 ;
	  if(pt>=600) return 0.0529581 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00263941 ;
	  if(pt>=50 && pt<70) return 0.00139648 ;
	  if(pt>=70 && pt<100) return 0.0012187 ;
	  if(pt>=100 && pt<140) return 0.00119532 ;
	  if(pt>=140 && pt<200) return 0.00128022 ;
	  if(pt>=200 && pt<300) return 0.00152341 ;
	  if(pt>=300 && pt<600) return 0.00264386 ;
	  if(pt>=600) return 0.00653632 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00273203 ;
	  if(pt>=50 && pt<70) return 0.00161521 ;
	  if(pt>=70 && pt<100) return 0.00129373 ;
	  if(pt>=100 && pt<140) return 0.00139311 ;
	  if(pt>=140 && pt<200) return 0.00146086 ;
	  if(pt>=200 && pt<300) return 0.00158335 ;
	  if(pt>=300 && pt<600) return 0.00275362 ;
	  if(pt>=600) return 0.00686364 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00260354 ;
	  if(pt>=50 && pt<70) return 0.00177744 ;
	  if(pt>=70 && pt<100) return 0.00153607 ;
	  if(pt>=100 && pt<140) return 0.0016602 ;
	  if(pt>=140 && pt<200) return 0.00197096 ;
	  if(pt>=200 && pt<300) return 0.00235663 ;
	  if(pt>=300 && pt<600) return 0.00442991 ;
	  if(pt>=600) return 0.0120702 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00286177 ;
	  if(pt>=50 && pt<70) return 0.00204757 ;
	  if(pt>=70 && pt<100) return 0.00199336 ;
	  if(pt>=100 && pt<140) return 0.00233169 ;
	  if(pt>=140 && pt<200) return 0.00284861 ;
	  if(pt>=200 && pt<300) return 0.00351394 ;
	  if(pt>=300 && pt<600) return 0.00592714 ;
	  if(pt>=600) return 0.0124033 ;
	}
      }
    }
  else if(sample == "tt")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.635851 ;
	  if(pt>=50 && pt<70) return 0.678729 ;
	  if(pt>=70 && pt<100) return 0.69796 ;
	  if(pt>=100 && pt<140) return 0.712605 ;
	  if(pt>=140 && pt<200) return 0.728733 ;
	  if(pt>=200 && pt<300) return 0.728222 ;
	  if(pt>=300 && pt<600) return 0.687089 ;
	  if(pt>=600) return 0.569871 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.614468 ;
	  if(pt>=50 && pt<70) return 0.658643 ;
	  if(pt>=70 && pt<100) return 0.679534 ;
	  if(pt>=100 && pt<140) return 0.694633 ;
	  if(pt>=140 && pt<200) return 0.71227 ;
	  if(pt>=200 && pt<300) return 0.710596 ;
	  if(pt>=300 && pt<600) return 0.662768 ;
	  if(pt>=600) return 0.532541 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.54648 ;
	  if(pt>=50 && pt<70) return 0.589259 ;
	  if(pt>=70 && pt<100) return 0.608867 ;
	  if(pt>=100 && pt<140) return 0.622579 ;
	  if(pt>=140 && pt<200) return 0.643845 ;
	  if(pt>=200 && pt<300) return 0.641116 ;
	  if(pt>=300 && pt<600) return 0.580178 ;
	  if(pt>=600) return 0.434753 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.502321 ;
	  if(pt>=50 && pt<70) return 0.54417 ;
	  if(pt>=70 && pt<100) return 0.563793 ;
	  if(pt>=100 && pt<140) return 0.581948 ;
	  if(pt>=140 && pt<200) return 0.600118 ;
	  if(pt>=200 && pt<300) return 0.595121 ;
	  if(pt>=300 && pt<600) return 0.534394 ;
	  if(pt>=600) return 0.406315 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0300655 ;
	  if(pt>=50 && pt<70) return 0.0259742 ;
	  if(pt>=70 && pt<100) return 0.0258642 ;
	  if(pt>=100 && pt<140) return 0.0274179 ;
	  if(pt>=140 && pt<200) return 0.0310595 ;
	  if(pt>=200 && pt<300) return 0.03887 ;
	  if(pt>=300 && pt<600) return 0.052695 ;
	  if(pt>=600) return 0.0613515 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0306 ;
	  if(pt>=50 && pt<70) return 0.0267456 ;
	  if(pt>=70 && pt<100) return 0.0269237 ;
	  if(pt>=100 && pt<140) return 0.0288964 ;
	  if(pt>=140 && pt<200) return 0.0323185 ;
	  if(pt>=200 && pt<300) return 0.0391772 ;
	  if(pt>=300 && pt<600) return 0.0512272 ;
	  if(pt>=600) return 0.0561889 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0297573 ;
	  if(pt>=50 && pt<70) return 0.0256131 ;
	  if(pt>=70 && pt<100) return 0.0261067 ;
	  if(pt>=100 && pt<140) return 0.028029 ;
	  if(pt>=140 && pt<200) return 0.0321702 ;
	  if(pt>=200 && pt<300) return 0.0383747 ;
	  if(pt>=300 && pt<600) return 0.0491744 ;
	  if(pt>=600) return 0.0483128 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0293692 ;
	  if(pt>=50 && pt<70) return 0.02469 ;
	  if(pt>=70 && pt<100) return 0.0240683 ;
	  if(pt>=100 && pt<140) return 0.0264314 ;
	  if(pt>=140 && pt<200) return 0.0301735 ;
	  if(pt>=200 && pt<300) return 0.0367625 ;
	  if(pt>=300 && pt<600) return 0.0451669 ;
	  if(pt>=600) return 0.0429804 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00199518 ;
	  if(pt>=50 && pt<70) return 0.00109405 ;
	  if(pt>=70 && pt<100) return 0.000980188 ;
	  if(pt>=100 && pt<140) return 0.000999707 ;
	  if(pt>=140 && pt<200) return 0.00109957 ;
	  if(pt>=200 && pt<300) return 0.00134789 ;
	  if(pt>=300 && pt<600) return 0.00243716 ;
	  if(pt>=600) return 0.00751583 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00207121 ;
	  if(pt>=50 && pt<70) return 0.00123541 ;
	  if(pt>=70 && pt<100) return 0.00112058 ;
	  if(pt>=100 && pt<140) return 0.00114681 ;
	  if(pt>=140 && pt<200) return 0.00127083 ;
	  if(pt>=200 && pt<300) return 0.00146729 ;
	  if(pt>=300 && pt<600) return 0.00257511 ;
	  if(pt>=600) return 0.00710372 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00208823 ;
	  if(pt>=50 && pt<70) return 0.00139718 ;
	  if(pt>=70 && pt<100) return 0.00133492 ;
	  if(pt>=100 && pt<140) return 0.00146313 ;
	  if(pt>=140 && pt<200) return 0.001713 ;
	  if(pt>=200 && pt<300) return 0.00215588 ;
	  if(pt>=300 && pt<600) return 0.00427701 ;
	  if(pt>=600) return 0.0114294 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00268738 ;
	  if(pt>=50 && pt<70) return 0.00183769 ;
	  if(pt>=70 && pt<100) return 0.00177367 ;
	  if(pt>=100 && pt<140) return 0.00205706 ;
	  if(pt>=140 && pt<200) return 0.00257445 ;
	  if(pt>=200 && pt<300) return 0.00346758 ;
	  if(pt>=300 && pt<600) return 0.00636822 ;
	  if(pt>=600) return 0.0115531 ;
	}
      }
    }
  else if(sample == "dib")
    {
      if(flavor==5){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.621745 ;
	  if(pt>=50 && pt<70) return 0.659441 ;
	  if(pt>=70 && pt<100) return 0.693269 ;
	  if(pt>=100 && pt<140) return 0.700953 ;
	  if(pt>=140 && pt<200) return 0.711653 ;
	  if(pt>=200 && pt<300) return 0.715527 ;
	  if(pt>=300 && pt<600) return 0.730902 ;
	  if(pt>=600) return 0.767574 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.598801 ;
	  if(pt>=50 && pt<70) return 0.663723 ;
	  if(pt>=70 && pt<100) return 0.682136 ;
	  if(pt>=100 && pt<140) return 0.70294 ;
	  if(pt>=140 && pt<200) return 0.693904 ;
	  if(pt>=200 && pt<300) return 0.700736 ;
	  if(pt>=300 && pt<600) return 0.714284 ;
	  if(pt>=600) return 0.697737 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.549371 ;
	  if(pt>=50 && pt<70) return 0.592631 ;
	  if(pt>=70 && pt<100) return 0.608454 ;
	  if(pt>=100 && pt<140) return 0.630803 ;
	  if(pt>=140 && pt<200) return 0.63474 ;
	  if(pt>=200 && pt<300) return 0.630934 ;
	  if(pt>=300 && pt<600) return 0.640124 ;
	  if(pt>=600) return 0.591139 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.484292 ;
	  if(pt>=50 && pt<70) return 0.553981 ;
	  if(pt>=70 && pt<100) return 0.550685 ;
	  if(pt>=100 && pt<140) return 0.571858 ;
	  if(pt>=140 && pt<200) return 0.575818 ;
	  if(pt>=200 && pt<300) return 0.598651 ;
	  if(pt>=300 && pt<600) return 0.58943 ;
	  if(pt>=600) return 0.569875 ;
	}
      }
      else if(flavor==4){
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0301811 ;
	  if(pt>=50 && pt<70) return 0.0253151 ;
	  if(pt>=70 && pt<100) return 0.0250534 ;
	  if(pt>=100 && pt<140) return 0.0290466 ;
	  if(pt>=140 && pt<200) return 0.0323913 ;
	  if(pt>=200 && pt<300) return 0.0396841 ;
	  if(pt>=300 && pt<600) return 0.0554936 ;
	  if(pt>=600) return 0.0561431 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0290787 ;
	  if(pt>=50 && pt<70) return 0.0257023 ;
	  if(pt>=70 && pt<100) return 0.0277957 ;
	  if(pt>=100 && pt<140) return 0.0266877 ;
	  if(pt>=140 && pt<200) return 0.0304117 ;
	  if(pt>=200 && pt<300) return 0.0389465 ;
	  if(pt>=300 && pt<600) return 0.0506215 ;
	  if(pt>=600) return 0.0623864 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0287913 ;
	  if(pt>=50 && pt<70) return 0.0269391 ;
	  if(pt>=70 && pt<100) return 0.0252978 ;
	  if(pt>=100 && pt<140) return 0.0265234 ;
	  if(pt>=140 && pt<200) return 0.0315922 ;
	  if(pt>=200 && pt<300) return 0.0384182 ;
	  if(pt>=300 && pt<600) return 0.0524296 ;
	  if(pt>=600) return 0.0453584 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.0254435 ;
	  if(pt>=50 && pt<70) return 0.0254436 ;
	  if(pt>=70 && pt<100) return 0.0251259 ;
	  if(pt>=100 && pt<140) return 0.0237963 ;
	  if(pt>=140 && pt<200) return 0.0312452 ;
	  if(pt>=200 && pt<300) return 0.0323776 ;
	  if(pt>=300 && pt<600) return 0.0485668 ;
	  if(pt>=600) return 0.0523895 ;
	}
      }
      else {
	if(fabs(eta) >= 0&& fabs(eta)< 0.6){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00242421 ;
	  if(pt>=50 && pt<70) return 0.00117227 ;
	  if(pt>=70 && pt<100) return 0.00115708 ;
	  if(pt>=100 && pt<140) return 0.000916735 ;
	  if(pt>=140 && pt<200) return 0.0011464 ;
	  if(pt>=200 && pt<300) return 0.00152216 ;
	  if(pt>=300 && pt<600) return 0.002748 ;
	  if(pt>=600) return 0.00669757 ;
	}
	if(fabs(eta) >= 0.6&& fabs(eta)< 1.2){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00225878 ;
	  if(pt>=50 && pt<70) return 0.00121443 ;
	  if(pt>=70 && pt<100) return 0.00101182 ;
	  if(pt>=100 && pt<140) return 0.00118308 ;
	  if(pt>=140 && pt<200) return 0.00131998 ;
	  if(pt>=200 && pt<300) return 0.00139654 ;
	  if(pt>=300 && pt<600) return 0.00282746 ;
	  if(pt>=600) return 0.00648637 ;
	}
	if(fabs(eta) >= 1.2&& fabs(eta)< 1.9){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00214507 ;
	  if(pt>=50 && pt<70) return 0.00160623 ;
	  if(pt>=70 && pt<100) return 0.00158697 ;
	  if(pt>=100 && pt<140) return 0.00172384 ;
	  if(pt>=140 && pt<200) return 0.00190063 ;
	  if(pt>=200 && pt<300) return 0.00214862 ;
	  if(pt>=300 && pt<600) return 0.00456314 ;
	  if(pt>=600) return 0.0084323 ;
	}
	if(fabs(eta) >= 1.9&& fabs(eta)< 2.5){
	  if(pt>=20 && pt<30) return 0 ;
	  if(pt>=30 && pt<50) return 0.00307277 ;
	  if(pt>=50 && pt<70) return 0.00207509 ;
	  if(pt>=70 && pt<100) return 0.00211169 ;
	  if(pt>=100 && pt<140) return 0.0023582 ;
	  if(pt>=140 && pt<200) return 0.00293513 ;
	  if(pt>=200 && pt<300) return 0.00293857 ;
	  if(pt>=300 && pt<600) return 0.0064019 ;
	  if(pt>=600) return 0.0208246 ;
	}
      }
    }
  return 1.0;
}
double DeepAK8HadTopTag_SF(float pt){
  
// scale factors taken from DAKX twiki https://indico.cern.ch/event/877167/contributions/3744193/attachments/1989744/3379280/DeepAK8_Top_W_SFs_V2.pdf contains sf for all years

/// for WP 0.1%
if(pt>=300 & pt<400) return 0.81;
else if(pt>=400 & pt<480) return 0.93;
else if(pt>=480 & pt<600) return 0.96;
else if(pt>=600) return 0.93;

return 1;
}

double DeepAK8HadTopTag_SF_UP(float pt){

// scale factors taken from DAKX twiki https://indico.cern.ch/event/877167/contributions/3744193/attachments/1989744/3379280/DeepAK8_Top_W_SFs_V2.pdf contains sf for all years

/// for WP 0.1%
if(pt>=300 & pt<400) return 0.88;
else if(pt>=400 & pt<480) return 0.98;
else if(pt>=480 & pt<600) return 1.02;
else if(pt>=600) return 0.98;

return 1;
}

double DeepAK8HadTopTag_SF_DOWN(float pt){

// scale factors taken from DAKX twiki https://indico.cern.ch/event/877167/contributions/3744193/attachments/1989744/3379280/DeepAK8_Top_W_SFs_V2.pdf contains sf for all years

/// for WP 0.1%
if(pt>=300 & pt<400) return 0.74;
else if(pt>=400 & pt<480) return 0.88;
else if(pt>=480 & pt<600) return 0.92;
else if(pt>=600) return 0.88;

return 1;
}

double PNetHadTopTag_SF(float pt){

// scale factors taken from https://indico.cern.ch/event/1103765/contributions/4647556/attachments/2364610/4037250/ParticleNet_2018_ULNanoV9_JMAR_14Dec2021_PK.pdf contains sf for all years

/// for WP 0.1%
if(pt>=300 & pt<400) return 1.04;
else if(pt>=400 & pt<480) return 0.95;
else if(pt>=480 & pt<600) return 0.94;
else if(pt>=600) return 0.96;

return 1;
}

double PNetHadTopTag_SF_UP(float pt){

// scale factors taken from https://indico.cern.ch/event/1103765/contributions/4647556/attachments/2364610/4037250/ParticleNet_2018_ULNanoV9_JMAR_14Dec2021_PK.pdf contains sf for all years

/// for WP 0.1%
if(pt>=300 & pt<400) return 1.13;
else if(pt>=400 & pt<480) return 1.0;
else if(pt>=480 & pt<600) return 1.02;
else if(pt>=600) return 1.03;

return 1;
}

double PNetHadTopTag_SF_DOWN(float pt){

// scale factors taken https://indico.cern.ch/event/1103765/contributions/4647556/attachments/2364610/4037250/ParticleNet_2018_ULNanoV9_JMAR_14Dec2021_PK.pdf contains sf for all years

/// for WP 0.1%
if(pt>=300 & pt<400) return 0.96;
else if(pt>=400 & pt<480) return 0.91;
else if(pt>=480 & pt<600) return 0.84;
else if(pt>=600) return 0.91;

return 1;
}

double SF_TOP(double alpha, double beta, double pt0, double pt1)
{
  double sfwt = sqrt(exp(alpha-beta*pt0) * exp(alpha-beta*pt1));
  return sfwt;
}

void Thrust_calculate(std::vector<TLorentzVector>& Input_4v, std::vector<double>& thrust_meas) {
  
  std::vector<double> Input_Px;
  std::vector<double> Input_Py;
  
  double thrustmax_calc =0;
  double temp_calc =0;
  unsigned int length_thrust_calc =0;
  
  std::vector<double> Thrust_Axis_calc;
  std::vector<double> p_thrust_max_calc;
  std::vector<double> p_dec_1_calc;
  std::vector<double> p_dec_2_calc;
  std::vector<double> p_pt_beam_calc;
  
  if (!thrust_meas.empty()){
    thrust_meas.clear();
    Thrust_Axis_calc.clear();
    p_thrust_max_calc.clear();
    p_dec_1_calc.clear();
    p_dec_2_calc.clear();
    p_pt_beam_calc.clear();
  }
  for(unsigned int j = 0; j < 3; j++){
    p_pt_beam_calc.push_back(0.);
    p_dec_1_calc.push_back(0.);
    p_dec_2_calc.push_back(0.);
    p_thrust_max_calc.push_back(0.);
    Thrust_Axis_calc.push_back(0.);
  }
  
  length_thrust_calc = Input_4v.size();
  for (unsigned k=0; k<length_thrust_calc; k++) {
    Input_Px.push_back(Input_4v[k].Px());
    Input_Py.push_back(Input_4v[k].Py());
  }
  
  double Pt_sum_calc =0;
  
  for(unsigned int k=0;k<length_thrust_calc;k++){
    Pt_sum_calc+=sqrt(pow(Input_Px[k],2)+pow(Input_Py[k],2));
    for(unsigned int j = 0; j < 3; j++){
      p_thrust_max_calc[j]=0;
    }
    //get a vector perpendicular to the beam axis and                            
	//perpendicular to the momentum of particle k                                
	    //per default beam axis b = (0,0,1)                                          
    
		p_pt_beam_calc[0] = Input_Py[k]*1;
    p_pt_beam_calc[1] = - Input_Px[k]*1;
    p_pt_beam_calc[2] = 0.; // GMA p_pt_beam_calc[3] = 0.;                       
    
    for(unsigned int i=0;i<length_thrust_calc;i++){
      if(i!=k){
	if((Input_Px[i]*p_pt_beam_calc[0]+Input_Py[i]*p_pt_beam_calc[1])>=0){
	  p_thrust_max_calc[0]= p_thrust_max_calc[0]+Input_Px[i];
	  p_thrust_max_calc[1]= p_thrust_max_calc[1]+Input_Py[i];
	}
	else{
	  p_thrust_max_calc[0]= p_thrust_max_calc[0]-Input_Px[i];
	  p_thrust_max_calc[1]= p_thrust_max_calc[1]-Input_Py[i];
	}
      }
    }
    
    p_dec_1_calc[0]=p_thrust_max_calc[0]+Input_Px[k];
    p_dec_1_calc[1]=p_thrust_max_calc[1]+Input_Py[k];
    p_dec_1_calc[2]=0;
    
    p_dec_2_calc[0]=p_thrust_max_calc[0]-Input_Px[k];
    p_dec_2_calc[1]=p_thrust_max_calc[1]-Input_Py[k];
    p_dec_2_calc[2]=0;
    
    temp_calc = pow(p_dec_1_calc[0],2)+pow(p_dec_1_calc[1],2);
    
    if(temp_calc>thrustmax_calc){
      thrustmax_calc =temp_calc;
      
      for(unsigned int i=0;i<3;i++){
	Thrust_Axis_calc[i]=p_dec_1_calc[i]/sqrt(thrustmax_calc);
      }
    }
    temp_calc = pow(p_dec_2_calc[0],2)+pow(p_dec_2_calc[1],2);
    
    if(temp_calc>thrustmax_calc){
      thrustmax_calc =temp_calc;
      for(unsigned int i=0;i<3;i++){
	Thrust_Axis_calc[i]=p_dec_2_calc[i]/sqrt(thrustmax_calc);
      }
    }
  }
  for(unsigned int j=0;j<3;j++){
    thrust_meas[j]=Thrust_Axis_calc[j];
  }
  double thrust_calc=0;
  thrust_calc = sqrt(thrustmax_calc)/Pt_sum_calc;
  thrust_meas[3]=1.-thrust_calc;
  if (thrust_meas[3] < 1.e-20) thrust_meas[3] = 1.e-20;
}

class top_tagger_sf : public TSelector {
public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  
  //New more variables stored as ntuple//
  float M_l1l2, rat_l2pt_l1pt, deltaPhi_l1l2, l1pt_nearjet, l2pt_nearjet, met_pt, met_phi, met_pt_puppi, met_phi_puppi, delta_phil1_met, delta_phil2_met, delta_phibl1_met, delta_phibl2_met, rat_metpt_ak4pt, rat_metpt_ak8pt, rat_metpt_eventHT, mt_of_l1met, mt_of_l2met, no_ak4jets, no_ak4bjets, no_ak8jets, EventHT, extra_ak4j, ptsum_extra_ak4, extra_ak4jqgl, extra_ak4jdeepb, rat_extra_ak4jpt_lpt, ak81pt, ak81y, ak81mass, ak81sdmass, ak81deep_tvsqcd, ak81deep_wvsqcd, ak82pt, ak82y, ak82mass, ak82sdmass, ak82deep_tvsqcd, ak82deep_wvsqcd, M_bl1, M_bl2, M_jl1, M_jl2, delta_phibl1bl2, delta_phijl1jl2, deltaR_l1l2, /*deltaR_l1b1, deltaR_l2b1, deltaR_l1b2, deltaR_l2b2,*/ deltaR_l1j1, deltaR_l2j1, deltaR_l1j2, deltaR_l2j2,deltaR_j1j2 , j1_btag_sc, j2_btag_sc, j1_btag_sc_ptsort, j2_btag_sc_ptsort, ak81NHad , ak81chrad ,ak81neuhad , ak81tau21 , ak81subhaddiff , ak81tau32 , ak82NHad , ak82chrad , ak82neuhad , ak82tau21 , ak82subhaddiff , ak82tau32,ak41pt,ak42pt,ak41mass,ak42mass,ak41delrlep,ak42delrlep,ak41inleppt, ak42inleppt,lep1pt,lep2pt, response_ak82,response_ak81,extra_ak4mass,extra_ak4pt,tt_decay_mode, tau_decay_mode,response_eventclass, pt_gen_l1, pt_gen_b1, pt_gen_lb1, deltaR_gen_lb1, pt_gen_l2, pt_gen_b2, pt_gen_lb2, deltaR_gen_lb2, delptbypt_gen_l1, delptbypt_gen_b1, delptbypt_gen_lb1, delptbypt_gen_l2, delptbypt_gen_b2, delptbypt_gen_lb2,puWeight,puWeightup,puWeightdown,btag_SF_dn,btag_SF,btag_SF_up;

  float in_M_l1l2, in_rat_l2pt_l1pt, in_deltaPhi_l1l2, in_l1pt_nearjet, in_l2pt_nearjet, in_met_pt, in_met_phi, in_delta_phil1_met, in_delta_phil2_met, in_delta_phibl1_met, in_delta_phibl2_met, in_rat_metpt_ak4pt, in_rat_metpt_ak8pt, in_rat_metpt_eventHT, in_mt_of_l1met, in_mt_of_l2met, in_no_ak4jets, in_no_ak4bjets, in_no_ak8jets, in_EventHT, in_extra_ak4j, in_ptsum_extra_ak4, in_extra_ak4jqgl, in_extra_ak4jdeepb, in_rat_extra_ak4jpt_lpt, in_ak81pt, in_ak81y, in_ak81mass, in_ak81sdmass, in_ak81deep_tvsqcd, in_ak81deep_wvsqcd, in_ak82pt, in_ak82y, in_ak82mass, in_ak82sdmass, in_ak82deep_tvsqcd, in_ak82deep_wvsqcd, in_M_bl1, in_M_bl2, in_M_jl1, in_M_jl2, in_delta_phibl1bl2, in_delta_phijl1jl2, in_deltaR_l1l2, /*deltaR_l1b1, deltaR_l2b1, deltaR_l1b2, deltaR_l2b2,*/ in_deltaR_l1j1, in_deltaR_l2j1, in_deltaR_l1j2, in_deltaR_l2j2,in_deltaR_j1j2 , in_j1_btag_sc, in_j2_btag_sc, in_j1_btag_sc_ptsort, in_j2_btag_sc_ptsort, in_ak81NHad , in_ak81chrad ,in_ak81neuhad , in_ak81tau21 , in_ak81subhaddiff , in_ak81tau32 , in_ak82NHad , in_ak82chrad , in_ak82neuhad , in_ak82tau21 , in_ak82subhaddiff , in_ak82tau32,in_ak41pt,in_ak42pt,in_ak41mass,in_ak42mass,in_ak41delrlep,in_ak42delrlep,in_ak41inleppt, in_ak42inleppt,in_lep1pt,in_lep2pt, in_response_ak82,in_response_ak81,in_extra_ak4mass,in_extra_ak4pt, in_dirgltrthr, in_dirglthrmin,nmu_veto,nel_veto,weight_mutop,weight_eltop,hadtaggersf,deltaphi,nbjets,tagjetpt,tagjet_taggerscore,probejet_score,nelectrons_eltop,nmuons_mutop;

  bool genmatch_tag,genmatch_probe,topjetmatch,elmatch,mumatch,bjetmatch;

  int ak81hase,ak81hasmu,ak82hase,ak82hasmu,ak81hasgene,ak81hasgenmu,ak81hasgenb,ak81hasgenhasalldecay,ak82hasgene,ak82hasgenmu,ak82hasgenb,ak82hasgenhasalldecay;
  int genmatch_sc2,genmatch_sc1,topmatchvar;

  float dirgltrthr, dirglthrmin;
  
  static const int njetmx = 20;
  static const int njetmxAK8 = 10;
  static const int npartmx = 25;
  static const int ntrigobjmx = 25;
  
  float weight;
  float leptonsf_weight;
  float leptonsf_weight_stat = 1.0;
  float leptonsf_weight_syst = 1.0;
  float weight_puwup;
  float weight_puwdown;
  float weight_btagwup;
  float weight_btagwdown;
  float event_pt_weight;

  bool isMC;
  bool SemiLeptt;
  bool DiLeptt;
  bool Hadtt;
  bool EE, EMU, MUMU, EJets, MUJets, TAUJets, TauTau, ETau, MuTau;
  
  bool isTT;
  bool isTTsignal;
  bool isST;
  bool isDIB;
  bool isWJ;
  bool isDY;
  bool isQCD;
  bool isTTX;
  
  Bool_t hlt_PFMET250,hlt_PFMET300,hlt_PFMET200,hlt_PFMET200_TypeOne ;

  Double_t prefiringweight, prefiringweightup, prefiringweightdown, LHE_weight;

  Int_t nLHEScaleWeights, nLHEPDFWeights, nLHEAlpsWeights,nLHEPSWeights;

  Float_t chs_met,chs_misphi, chs_misetsig,chs_sumEt,miset_UnclusEup,misphi_UnclusEdn,miset_PUPPI,misphi_PUPPI,misetsig_PUPPI,sumEt_PUPPI,miset_PUPPI_JESup,miset_PUPPI_JESdn, miset_PUPPI_JERup, miset_PUPPI_JERdn, miset_PUPPI_UnclusEup, miset_PUPPI_UnclusEdn, misphi_PUPPI_JESup, misphi_PUPPI_JESdn, misphi_PUPPI_JERup, misphi_PUPPI_JERdn, misphi_PUPPI_UnclusEup, misphi_PUPPI_UnclusEdn, misphi_UnclusEup, miset_UnclusEdn;

  Float_t Generator_x1,Generator_x2,Generator_xpdf1,Generator_xpdf2,Generator_id2,Generator_id1,Generator_scalePDF,gen_weight,toppt_wt;

  Float_t Muon_corrected_pt[njetmx], Muon_correctedUp_pt[njetmx], Muon_correctedDown_pt[njetmx], Electron_eccalTrkEnergyPostCorr[njetmx], Electron_energyScaleValue[njetmx], Electron_energyScaleUp[njetmx], Electron_energyScaleDown[njetmx], Electron_energySigmaValue[njetmx], Electron_energySigmaUp[njetmx], Electron_energySigmaDown[njetmx];

  Float_t LHEScaleWeights[npartmx],LHEPDFWeights[npartmx],LHEAlpsWeights[npartmx],LHEPSWeights[npartmx];

  // Declaration of leaf types
  Int_t           irun;
  Int_t           ilumi;
  UInt_t          ievt;
  Int_t           nprim;
  Int_t           nprimi;
  Double_t        Rho;
  Double_t        event_weight;
  Float_t         qscale;
  Int_t           npu_vert;
  Int_t           npu_vert_true;
  Int_t           trig_value;
  
  Bool_t          hlt_IsoMu24;
  Bool_t          hlt_Mu50;
  Bool_t          hlt_Ele50_PFJet165;
  Bool_t          hlt_Ele115;
  Bool_t          hlt_AK8PFJet500;
  Bool_t          hlt_Photon200;
  Bool_t          hlt_Mu37Ele27;
  Bool_t          hlt_Mu27Ele37;
  Bool_t          hlt_Mu37TkMu27;
  Bool_t          hlt_OldMu100;
  Bool_t          hlt_TkMu100;
  Bool_t          hlt_DoubleEle25;

  Int_t           ntrigobjs;
  Float_t         trigobjpt[ntrigobjmx];   //[ntrigobjs]
  Float_t         trigobjeta[ntrigobjmx];   //[ntrigobjs]
  Float_t         trigobjphi[ntrigobjmx];   //[ntrigobjs]
  Float_t         trigobjmass[ntrigobjmx];   //[ntrigobjs]
                 
  Int_t           trigobjpdgId[ntrigobjmx];   //[ntrigobjs]
  Bool_t          trigobjHLT[ntrigobjmx];   //[ntrigobjs]
  Bool_t          trigobjL1[ntrigobjmx];   //[ntrigobjs]
  Bool_t          trigobjBoth[ntrigobjmx];   //[ntrigobjs]
  Int_t           trigobjIhlt[ntrigobjmx];   //[ntrigobjs]
  Float_t         PFMET;
  Float_t         PFMETPhi;
  Float_t         MisEtSig;
  Float_t         sumEt;
  Int_t           npfjetAK8;
  Bool_t          pfjetAK8jetID[njetmxAK8];
  Bool_t          pfjetAK8jetID_tightlepveto[njetmxAK8];
  Float_t         pfjetAK8pt[njetmxAK8];
  Float_t         pfjetAK8y[njetmxAK8];
  Float_t         pfjetAK8eta[njetmxAK8];
  Float_t         pfjetAK8phi[njetmxAK8];
  Float_t         pfjetAK8mass[njetmxAK8];
  Float_t         pfjetAK8JEC[njetmxAK8];
  Float_t         pfjetAK8btag_DeepCSV[njetmxAK8];
  Float_t 	  pfjetAK8matchAK4deepb[njetmxAK8];
  Float_t         pfjetAK8DeepTag_TvsQCD[njetmxAK8];
  Float_t         pfjetAK8DeepTag_WvsQCD[njetmxAK8];
  Float_t         pfjetAK8DeepTag_ZvsQCD[njetmxAK8];
  Float_t         pfjetAK8CHF[njetmxAK8];
  Float_t         pfjetAK8NHF[njetmxAK8];
  Float_t         pfjetAK8CEMF[njetmxAK8];
  Float_t         pfjetAK8NEMF[njetmxAK8];
  Float_t         pfjetAK8MUF[njetmxAK8];
  Float_t         pfjetAK8PHF[njetmxAK8];
  Float_t         pfjetAK8EEF[njetmxAK8];
  Float_t         pfjetAK8HFHF[njetmxAK8];

  Int_t           pfjetAK8CHM[njetmxAK8];
  Int_t           pfjetAK8NHM[njetmxAK8];
  Int_t           pfjetAK8MUM[njetmxAK8];
  Int_t           pfjetAK8PHM[njetmxAK8];
  Int_t           pfjetAK8EEM[njetmxAK8];
  Int_t           pfjetAK8HFHM[njetmxAK8];
  Int_t           pfjetAK8Neucons[njetmxAK8];
  Int_t           pfjetAK8Chcons[njetmxAK8];

  Float_t         pfjetAK8reso[njetmxAK8];
  Float_t         pfjetAK8resoup[njetmxAK8];
  Float_t         pfjetAK8resodn[njetmxAK8];

  Float_t         pfjetAK8jesup_pu[njetmxAK8];
  Float_t         pfjetAK8jesup_rel[njetmxAK8];
  Float_t         pfjetAK8jesup_scale[njetmxAK8];
  Float_t         pfjetAK8jesup_total[njetmxAK8];
  Float_t         pfjetAK8jesdn_pu[njetmxAK8];
  Float_t         pfjetAK8jesdn_rel[njetmxAK8];
  Float_t         pfjetAK8jesdn_scale[njetmxAK8];
  Float_t         pfjetAK8jesdn_total[njetmxAK8];

  Float_t         pfjetAK8chrad[njetmxAK8];
  Float_t         pfjetAK8pTD[njetmxAK8];
  Float_t         pfjetAK8sdmass[njetmxAK8];
  bool		  pfjetAK8haspfelectron[njetmxAK8];
  bool		  pfjetAK8haspfmuon[njetmxAK8];

  Float_t         pfjetAK8tau1[njetmxAK8];
  Float_t         pfjetAK8tau2[njetmxAK8];
  Float_t         pfjetAK8tau3[njetmxAK8];
  Float_t         pfjetAK8tau21[njetmxAK8];
  Float_t         pfjetAK8tau32[njetmxAK8];
  Float_t         pfjetAK8sub1pt[njetmxAK8];
  Float_t         pfjetAK8sub1eta[njetmxAK8];
  Float_t         pfjetAK8sub1phi[njetmxAK8];
  Float_t         pfjetAK8sub1mass[njetmxAK8];
  Float_t         pfjetAK8sub1btag[njetmxAK8];
  
  Float_t         pfjetAK8sub1chhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub1neuhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub1emfrac[njetmxAK8];
  Float_t         pfjetAK8sub1phofrac[njetmxAK8];
  Float_t         pfjetAK8sub1mufrac[njetmxAK8];
 
  Float_t         pfjetAK8sub2pt[njetmxAK8];
  Float_t         pfjetAK8sub2eta[njetmxAK8];
  Float_t         pfjetAK8sub2phi[njetmxAK8];
  Float_t         pfjetAK8sub2mass[njetmxAK8];
  Float_t         pfjetAK8sub2btag[njetmxAK8];
  
  Float_t         pfjetAK8sub2chhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub2neuhadfrac[njetmxAK8];
  Float_t         pfjetAK8sub2emfrac[njetmxAK8];
  Float_t         pfjetAK8sub2phofrac[njetmxAK8];
  Float_t         pfjetAK8sub2mufrac[njetmxAK8];

  Float_t	  pfjetAK8subhaddiff[njetmxAK8];
  Float_t	  pfjetAK8subemdiff[njetmxAK8];
  Float_t         pfjetAK8subptdiff[njetmxAK8];

  Bool_t 	  pfjetAK8_hasmatche[njetmxAK8];
  Bool_t          pfjetAK8_hasmatchmu[njetmxAK8];
  Float_t	  pfjetAK8stau21[njetmxAK8];
  Float_t	  pfjetAK8stau32[njetmxAK8];


  Float_t	  pfjetAK8subbtag[njetmxAK8];
  Bool_t 	  pfjetAK8haselectron[njetmxAK8];
  Bool_t 	  pfjetAK8hasmuon[njetmxAK8];
  Bool_t          pfjetAK8hastau[njetmxAK8];
  Bool_t	  pfjetAK8hasb[njetmxAK8];
  Bool_t	  pfjetAK8hasqg[njetmxAK8];
  Bool_t 	  pfjetAK8hashadtop[njetmxAK8];
  Bool_t 	  pfjetAK8hasleptop[njetmxAK8];
  Bool_t	  pfjetAK8hastop[njetmxAK8];
  Float_t	  pfjetAK8re_tvsb[njetmxAK8];
  Float_t	  pfjetAK8rnu_tvsb[njetmxAK8];
  Float_t	  pfjetAK8rt[njetmxAK8];
  Float_t         pfjetAK8rmu_tvsb[njetmxAK8];

  Float_t         pfjetAK8elinsubpt[njetmxAK8];
  Float_t         pfjetAK8elinsubeta[njetmxAK8];
  Float_t         pfjetAK8elinsubphi[njetmxAK8];
  
  Float_t         pfjetAK8elinsubmass[njetmxAK8];
  Float_t         pfjetAK8elinsubjpt[njetmxAK8];
  Float_t         pfjetAK8elinsubjeta[njetmxAK8];
  Float_t         pfjetAK8elinsubjphi[njetmxAK8];
  
  Float_t         pfjetAK8elinsubjmass[njetmxAK8];

  Float_t         pfjetAK8muinsubpt[njetmxAK8];
  Float_t         pfjetAK8muinsubeta[njetmxAK8];
  Float_t         pfjetAK8muinsubphi[njetmxAK8];
  Float_t         pfjetAK8muinsubmass[njetmxAK8];
  Float_t         pfjetAK8muinsubjpt[njetmxAK8];
  Float_t         pfjetAK8muinsubjeta[njetmxAK8];
  Float_t         pfjetAK8muinsubjphi[njetmxAK8];
  Float_t         pfjetAK8muinsubje[njetmxAK8];
  Float_t         pfjetAK8muinsubjmass[njetmxAK8];
  Float_t         pfjetAK8muinsubI0[njetmxAK8];
  Float_t         pfjetAK8muinsubInear[njetmxAK8];
  Float_t         pfjetAK8muinsubIfar[njetmxAK8];

  Int_t           npfjetAK4;
  Bool_t          pfjetAK4jetID[njetmx];
  Bool_t          pfjetAK4jetID_tightlepveto[njetmx];
  Float_t         pfjetAK4pt[njetmx];
  Float_t         pfjetAK4eta[njetmx];
  Float_t         pfjetAK4y[njetmx];
  Float_t         pfjetAK4phi[njetmx];
  Float_t         pfjetAK4mass[njetmx];
  Float_t         pfjetAK4JEC[njetmx];
  
  Float_t         pfjetAK4inleppt[njetmx];
  Float_t         pfjetAK4delrlep[njetmx];

  Float_t         pfjetAK4btag_DeepCSV[njetmx];
  Float_t         pfjetAK4btag_DeepFlav[njetmx];
  
  Float_t         pfjetAK4reso[njetmx];
  Float_t         pfjetAK4resoup[njetmx];
  Float_t         pfjetAK4resodn[njetmx];
  Float_t         pfjetAK4jesup_pu[njetmx];
  Float_t         pfjetAK4jesup_rel[njetmx];
  Float_t         pfjetAK4jesup_scale[njetmx];
  Float_t         pfjetAK4jesup_total[njetmx];
  Float_t         pfjetAK4jesdn_pu[njetmx];
  Float_t         pfjetAK4jesdn_rel[njetmx];
  Float_t         pfjetAK4jesdn_scale[njetmx];
  Float_t         pfjetAK4jesdn_total[njetmx];
  Int_t           pfjetAK4hadronflav[njetmx];
  Int_t           pfjetAK4partonflav[njetmx];
  Float_t         pfjetAK4qgl[njetmx];
  Float_t         pfjetAK4PUID[njetmx];
  Int_t           pfjetAK4GenMatch;

  Float_t         GENMET;
  Float_t         GENMETPhi;
   
  Int_t           ngenjetAK8;
  Float_t         genjetAK8pt[njetmx];
  Float_t         genjetAK8eta[njetmx];
  Float_t         genjetAK8phi[njetmx];
  Float_t         genjetAK8mass[njetmx];
  Float_t         genjetAK8sdmass[njetmx];
   
   
   
  Int_t           ngenjetAK4;
  Float_t         genjetAK4pt[njetmx];
  Float_t         genjetAK4eta[njetmx];
  Float_t         genjetAK4phi[njetmx];
  Float_t         genjetAK4mass[njetmx];
  Float_t         genjetAK4hadronflav[njetmx];
  Float_t         genjetAK4partonflav[njetmx];

  Int_t           ngenparticles;
  Int_t           genpartstatus[npartmx];
  Int_t           genpartpdg[npartmx];
  Int_t           genpartmompdg[npartmx];
  Int_t           genpartgrmompdg[npartmx];
  Int_t           genpartdaugno[npartmx];
  Bool_t          genpartfromhard[npartmx];
  Bool_t          genpartfromhardbFSR[npartmx];
  Bool_t          genpartisPromptFinalState[npartmx];
  Bool_t          genpartisLastCopyBeforeFSR[npartmx];
  Float_t         genpartpt[npartmx];
  Float_t         genparteta[npartmx];
  Float_t         genpartphi[npartmx];
  Float_t         genpartm[npartmx];
  Int_t           nLHEparticles;
  Int_t           LHEpartpdg[npartmx];
  Float_t         LHEpartpt[npartmx];
  Float_t         LHEparteta[npartmx];
  Float_t         LHEpartphi[npartmx];
  Float_t         LHEpartm[npartmx];

  Int_t           genpartpair[npartmx];  //Not using many of them now
  //Int_t           ngenelectrons;
  Float_t         genelectronpt[npartmx];
  Float_t         genelectroneta[npartmx];
  Float_t         genelectronphi[npartmx];
  Float_t         genelectronm[npartmx];
  Int_t	   nleptop;
  Int_t	   nhadtop;
  //Int_t           ngenmuons;
  Float_t         genmuonpt[npartmx];
  Float_t         genmuoneta[npartmx];
  Float_t         genmuonphi[npartmx];
  Float_t         genmuonm[npartmx];
  //Int_t           ngentaus;
  Float_t         gentaupt[npartmx];
  Float_t         gentaueta[npartmx];
  Float_t         gentauphi[npartmx];
  Float_t         gentaum[npartmx];
  //Int_t           ngenqgs;
  Float_t         genqgpt[npartmx];
  Float_t         genqgeta[npartmx];
  Float_t         genqgphi[npartmx];
  Float_t         genqgm[npartmx];
  //Int_t           ngenbs;
  Float_t         genbpt[npartmx];
  Float_t         genbeta[npartmx];
  Float_t         genbphi[npartmx];
  Float_t         genbm[npartmx];
  //Int_t	   ngentops;
  Float_t	   gentoppt[npartmx];
  Float_t	   gentopeta[npartmx];
  Float_t	   gentopphi[npartmx];
  Float_t	   gentopm[npartmx];
  Int_t	   gentopid[npartmx];
  Int_t           nmuons;
  Bool_t          muonisPF[njetmx];
  Bool_t          muonisGL[njetmx];
  Bool_t          muonisTRK[njetmx];
  Bool_t          muonisLoose[njetmx];
  Bool_t          muonisGoodGL[njetmx];
  Bool_t          muonisMed[njetmx];
  Bool_t          muonisTight[njetmx];
  Bool_t          muonisHighPt[njetmx];
  Bool_t          muonisHighPttrk[njetmx];
  Bool_t          muonisMedPr[njetmx];
  Float_t         muonpt[njetmx];
  Float_t         muonp[njetmx];
   
  Float_t         muoneta[njetmx];
  Float_t         muonphi[njetmx];
  Float_t         muoncharge[njetmx];
   
  Float_t         muontrkvtx[njetmx];
  Float_t         muondz[njetmx];
  Float_t         muonpter[njetmx];
  Float_t         muonchi[njetmx];
  Int_t           muonndf[njetmx];
  Float_t         muonecal[njetmx];
  Float_t         muonhcal[njetmx];
   
  Float_t         muonpfiso[njetmx];
   
  Float_t         muonposmatch[njetmx];
  Float_t         muontrkink[njetmx];
  Float_t         muonsegcom[njetmx];

  Float_t         muonpixhit[njetmx];
  Float_t         muonmst[njetmx];
  Float_t         muontrklay[njetmx];
  Float_t         muonvalfrac[njetmx];
  Float_t         muonchiso[njetmx];
  Float_t         muonnhiso[njetmx];
  Float_t         muonphiso[njetmx];
  Float_t         muonminisoall[njetmx];
  Float_t         mudxy_sv[njetmx];
  Float_t         muonhit[njetmx];

  Int_t           nelecs;
  Float_t         elpt[njetmx];
  Float_t         eldxytrk[njetmx];
  Float_t         eldztrk[njetmx];  
  Float_t         eldxy_sv[njetmx];
  Float_t         eleta[njetmx];
  Float_t         elphi[njetmx];
  Float_t         elp[njetmx];
  Float_t         ele[njetmx];
  Float_t         elcharge[njetmx];
  Bool_t          elmvaid[njetmx];
  Bool_t          elmvaid_Fallv2WP80[njetmx]; 
  Bool_t          elmvaid_noIso[njetmx];
  Bool_t          elmvaid_Fallv2WP80_noIso[njetmx];

  Float_t         elhovere[njetmx];
  Float_t         elchi[njetmx];
  Int_t           elndf[njetmx];
   
  Float_t         eletain[njetmx];
  Float_t         elphiin[njetmx];

  Float_t         elsupcl_eta[njetmx];
  Float_t         elsupcl_phi[njetmx];
  Float_t         elsupcl_rawE[njetmx];
  Float_t         elfbrem[njetmx];
   
  Float_t         eleoverp[njetmx];
  Float_t         elietaieta[njetmx];
  Float_t         elmisshits[njetmx];
   
  Float_t         elpfiso[njetmx];
  Float_t         pfjetAK8PNet_TvsQCD[njetmxAK8];
  Float_t         pfjetAK8PNet_WvsQCD[njetmxAK8];
  Bool_t          elcutid_veto[njetmx];
  Bool_t          elcutid_loose[njetmx];
  Bool_t          elcutid_med[njetmx];
  Bool_t          elcutid_tight[njetmx];
  Bool_t          elmvaid_wpLoose[njetmx];
  Bool_t          elmvaid_noIso_wpLoose[njetmx];
  Bool_t          elsupcl_e[njetmx];
  Int_t           genpartmomstatus[npartmx];
   
  //***al el ID variables****//
  Float_t         elsigmaieta[njetmx];
  Float_t         elsigmaiphi[njetmx];
  Float_t         elr9full[njetmx];
  Float_t         elsupcl_etaw[njetmx];
  Float_t         elsupcl_phiw[njetmx];
  Float_t         elhcaloverecal[njetmx];
  Float_t         elcloctftrkn[njetmx];
  Float_t         elcloctftrkchi2[njetmx];
  Float_t         ele1x5bye5x5[njetmx];
  Float_t         elnormchi2[njetmx];
  Float_t         elhitsmiss[njetmx];
  Float_t         eltrkmeasure[njetmx];
  Float_t         elconvtxprob[njetmx];
  Float_t         elecloverpout[njetmx];
  Float_t         elecaletrkmomentum[njetmx];
  Float_t         eldeltaetacltrkcalo[njetmx];
  Float_t         elsupcl_preshvsrawe[njetmx];
  Float_t         elpfisolsumphet[njetmx];
  Float_t         elpfisolsumchhadpt[njetmx];
  Float_t         elpfsiolsumneuhadet[njetmx];

  float pfjetAK8jesup_AbsoluteStat[njetmxAK8];     TBranch    *b_pfjetAK8jesup_AbsoluteStat;
  float pfjetAK8jesup_AbsoluteScale[njetmxAK8];   TBranch    *b_pfjetAK8jesup_AbsoluteScale;
  float pfjetAK8jesup_AbsoluteMPFBias[njetmxAK8];   TBranch    *b_pfjetAK8jesup_AbsoluteMPFBias;
  float pfjetAK8jesup_FlavorQCD[njetmxAK8];   TBranch    *b_pfjetAK8jesup_FlavorQCD;
  float pfjetAK8jesup_Fragmentation[njetmxAK8];   TBranch    *b_pfjetAK8jesup_Fragmentation;
  float pfjetAK8jesup_PileUpDataMC[njetmxAK8];   TBranch    *b_pfjetAK8jesup_PileUpDataMC;
  float pfjetAK8jesup_PileUpPtBB[njetmxAK8];   TBranch    *b_pfjetAK8jesup_PileUpPtBB;
  float pfjetAK8jesup_PileUpPtEC1[njetmxAK8];   TBranch    *b_pfjetAK8jesup_PileUpPtEC1;
  float pfjetAK8jesup_PileUpPtEC2[njetmxAK8];   TBranch    *b_pfjetAK8jesup_PileUpPtEC2;
  float pfjetAK8jesup_PileUpPtHF[njetmxAK8];   TBranch    *b_pfjetAK8jesup_PileUpPtHF;
  float pfjetAK8jesup_PileUpPtRef[njetmxAK8];   TBranch    *b_pfjetAK8jesup_PileUpPtRef;
  float pfjetAK8jesup_RelativeFSR[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeFSR;
  float pfjetAK8jesup_RelativeJEREC1[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeJEREC1;
  float pfjetAK8jesup_RelativeJEREC2[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeJEREC2;
  float pfjetAK8jesup_RelativeJERHF[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeJERHF;
  float pfjetAK8jesup_RelativePtBB[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativePtBB;
  float pfjetAK8jesup_RelativePtEC1[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativePtEC1;
  float pfjetAK8jesup_RelativePtEC2[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativePtEC2;
  float pfjetAK8jesup_RelativePtHF[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativePtHF;
  float pfjetAK8jesup_RelativeBal[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeBal;
  float pfjetAK8jesup_RelativeSample[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeSample;
  float pfjetAK8jesup_RelativeStatEC[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeStatEC;
  float pfjetAK8jesup_RelativeStatHF[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeStatHF;
  float pfjetAK8jesup_RelativeStatFSR[njetmxAK8];   TBranch    *b_pfjetAK8jesup_RelativeStatFSR;
  float pfjetAK8jesup_SinglePionECAL[njetmxAK8];   TBranch    *b_pfjetAK8jesup_SinglePionECAL;
  float pfjetAK8jesup_SinglePionHCAL[njetmxAK8];   TBranch    *b_pfjetAK8jesup_SinglePionHCAL;
  float pfjetAK8jesup_TimePtEta[njetmxAK8];   TBranch    *b_pfjetAK8jesup_TimePtEta;
  float pfjetAK8jesup_Total[njetmxAK8];   TBranch    *b_pfjetAK8jesup_Total;

  float pfjetAK8jesdn_AbsoluteStat[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_AbsoluteStat;
  float pfjetAK8jesdn_AbsoluteScale[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_AbsoluteScale;
  float pfjetAK8jesdn_AbsoluteMPFBias[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_AbsoluteMPFBias;
  float pfjetAK8jesdn_FlavorQCD[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_FlavorQCD;
  float pfjetAK8jesdn_Fragmentation[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_Fragmentation;
  float pfjetAK8jesdn_PileUpDataMC[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_PileUpDataMC;
  float pfjetAK8jesdn_PileUpPtBB[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_PileUpPtBB;
  float pfjetAK8jesdn_PileUpPtEC1[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_PileUpPtEC1;
  float pfjetAK8jesdn_PileUpPtEC2[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_PileUpPtEC2;
  float pfjetAK8jesdn_PileUpPtHF[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_PileUpPtHF;
  float pfjetAK8jesdn_PileUpPtRef[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_PileUpPtRef;
  float pfjetAK8jesdn_RelativeFSR[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeFSR;
  float pfjetAK8jesdn_RelativeJEREC1[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeJEREC1;
  float pfjetAK8jesdn_RelativeJEREC2[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeJEREC2;
  float pfjetAK8jesdn_RelativeJERHF[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeJERHF;
  float pfjetAK8jesdn_RelativePtBB[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativePtBB;
  float pfjetAK8jesdn_RelativePtEC1[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativePtEC1;
  float pfjetAK8jesdn_RelativePtEC2[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativePtEC2;
  float pfjetAK8jesdn_RelativePtHF[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativePtHF;
  float pfjetAK8jesdn_RelativeBal[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeBal;
  float pfjetAK8jesdn_RelativeSample[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeSample;
  float pfjetAK8jesdn_RelativeStatEC[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeStatEC;
  float pfjetAK8jesdn_RelativeStatHF[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeStatHF;
  float pfjetAK8jesdn_RelativeStatFSR[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_RelativeStatFSR;
  float pfjetAK8jesdn_SinglePionECAL[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_SinglePionECAL;
  float pfjetAK8jesdn_SinglePionHCAL[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_SinglePionHCAL;
  float pfjetAK8jesdn_TimePtEta[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_TimePtEta;
  float pfjetAK8jesdn_Total[njetmxAK8];   TBranch    *b_pfjetAK8jesdn_Total;

  float pfjetAK4jesup_AbsoluteStat[njetmx];     TBranch    *b_pfjetAK4jesup_AbsoluteStat;
  float pfjetAK4jesup_AbsoluteScale[njetmx];   TBranch    *b_pfjetAK4jesup_AbsoluteScale;
  float pfjetAK4jesup_AbsoluteMPFBias[njetmx];   TBranch    *b_pfjetAK4jesup_AbsoluteMPFBias;
  float pfjetAK4jesup_FlavorQCD[njetmx];   TBranch    *b_pfjetAK4jesup_FlavorQCD;
  float pfjetAK4jesup_Fragmentation[njetmx];   TBranch    *b_pfjetAK4jesup_Fragmentation;
  float pfjetAK4jesup_PileUpDataMC[njetmx];   TBranch    *b_pfjetAK4jesup_PileUpDataMC;
  float pfjetAK4jesup_PileUpPtBB[njetmx];   TBranch    *b_pfjetAK4jesup_PileUpPtBB;
  float pfjetAK4jesup_PileUpPtEC1[njetmx];   TBranch    *b_pfjetAK4jesup_PileUpPtEC1;
  float pfjetAK4jesup_PileUpPtEC2[njetmx];   TBranch    *b_pfjetAK4jesup_PileUpPtEC2;
  float pfjetAK4jesup_PileUpPtHF[njetmx];   TBranch    *b_pfjetAK4jesup_PileUpPtHF;
  float pfjetAK4jesup_PileUpPtRef[njetmx];   TBranch    *b_pfjetAK4jesup_PileUpPtRef;
  float pfjetAK4jesup_RelativeFSR[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeFSR;
  float pfjetAK4jesup_RelativeJEREC1[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeJEREC1;
  float pfjetAK4jesup_RelativeJEREC2[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeJEREC2;
  float pfjetAK4jesup_RelativeJERHF[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeJERHF;
  float pfjetAK4jesup_RelativePtBB[njetmx];   TBranch    *b_pfjetAK4jesup_RelativePtBB;
  float pfjetAK4jesup_RelativePtEC1[njetmx];   TBranch    *b_pfjetAK4jesup_RelativePtEC1;
  float pfjetAK4jesup_RelativePtEC2[njetmx];   TBranch    *b_pfjetAK4jesup_RelativePtEC2;
  float pfjetAK4jesup_RelativePtHF[njetmx];   TBranch    *b_pfjetAK4jesup_RelativePtHF;
  float pfjetAK4jesup_RelativeBal[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeBal;
  float pfjetAK4jesup_RelativeSample[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeSample;
  float pfjetAK4jesup_RelativeStatEC[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeStatEC;
  float pfjetAK4jesup_RelativeStatHF[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeStatHF;
  float pfjetAK4jesup_RelativeStatFSR[njetmx];   TBranch    *b_pfjetAK4jesup_RelativeStatFSR;
  float pfjetAK4jesup_SinglePionECAL[njetmx];   TBranch    *b_pfjetAK4jesup_SinglePionECAL;
  float pfjetAK4jesup_SinglePionHCAL[njetmx];   TBranch    *b_pfjetAK4jesup_SinglePionHCAL;
  float pfjetAK4jesup_TimePtEta[njetmx];   TBranch    *b_pfjetAK4jesup_TimePtEta;
  float pfjetAK4jesup_Total[njetmx];   TBranch    *b_pfjetAK4jesup_Total;

  float pfjetAK4jesdn_AbsoluteStat[njetmx];   TBranch    *b_pfjetAK4jesdn_AbsoluteStat;
  float pfjetAK4jesdn_AbsoluteScale[njetmx];   TBranch    *b_pfjetAK4jesdn_AbsoluteScale;
  float pfjetAK4jesdn_AbsoluteMPFBias[njetmx];   TBranch    *b_pfjetAK4jesdn_AbsoluteMPFBias;
  float pfjetAK4jesdn_FlavorQCD[njetmx];   TBranch    *b_pfjetAK4jesdn_FlavorQCD;
  float pfjetAK4jesdn_Fragmentation[njetmx];   TBranch    *b_pfjetAK4jesdn_Fragmentation;
  float pfjetAK4jesdn_PileUpDataMC[njetmx];   TBranch    *b_pfjetAK4jesdn_PileUpDataMC;
  float pfjetAK4jesdn_PileUpPtBB[njetmx];   TBranch    *b_pfjetAK4jesdn_PileUpPtBB;
  float pfjetAK4jesdn_PileUpPtEC1[njetmx];   TBranch    *b_pfjetAK4jesdn_PileUpPtEC1;
  float pfjetAK4jesdn_PileUpPtEC2[njetmx];   TBranch    *b_pfjetAK4jesdn_PileUpPtEC2;
  float pfjetAK4jesdn_PileUpPtHF[njetmx];   TBranch    *b_pfjetAK4jesdn_PileUpPtHF;
  float pfjetAK4jesdn_PileUpPtRef[njetmx];   TBranch    *b_pfjetAK4jesdn_PileUpPtRef;
  float pfjetAK4jesdn_RelativeFSR[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeFSR;
  float pfjetAK4jesdn_RelativeJEREC1[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeJEREC1;
  float pfjetAK4jesdn_RelativeJEREC2[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeJEREC2;
  float pfjetAK4jesdn_RelativeJERHF[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeJERHF;
  float pfjetAK4jesdn_RelativePtBB[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativePtBB;
  float pfjetAK4jesdn_RelativePtEC1[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativePtEC1;
  float pfjetAK4jesdn_RelativePtEC2[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativePtEC2;
  float pfjetAK4jesdn_RelativePtHF[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativePtHF;
  float pfjetAK4jesdn_RelativeBal[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeBal;
  float pfjetAK4jesdn_RelativeSample[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeSample;
  float pfjetAK4jesdn_RelativeStatEC[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeStatEC;
  float pfjetAK4jesdn_RelativeStatHF[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeStatHF;
  float pfjetAK4jesdn_RelativeStatFSR[njetmx];   TBranch    *b_pfjetAK4jesdn_RelativeStatFSR;
  float pfjetAK4jesdn_SinglePionECAL[njetmx];   TBranch    *b_pfjetAK4jesdn_SinglePionECAL;
  float pfjetAK4jesdn_SinglePionHCAL[njetmx];   TBranch    *b_pfjetAK4jesdn_SinglePionHCAL;
  float pfjetAK4jesdn_TimePtEta[njetmx];   TBranch    *b_pfjetAK4jesdn_TimePtEta;
  float pfjetAK4jesdn_Total[njetmx];   TBranch    *b_pfjetAK4jesdn_Total;

  // List of branches
  TBranch        *b_irun;   //!
  TBranch        *b_ilumi;   //!
  TBranch        *b_ievt;   //!
  TBranch        *b_nprim;   //!
  TBranch        *b_nprimi;
  TBranch        *b_Rho;   //!
  TBranch        *b_event_weight;   //!
  TBranch        *b_qscale;   //!
  TBranch        *b_npu_vert;   //!
  TBranch        *b_npu_vert_true;   //!
  TBranch        *b_trig_value;   //!

  TBranch        *b_hlt_IsoMu24;
  TBranch        *b_hlt_Mu50;
  TBranch        *b_hlt_Ele50_PFJet165;
  TBranch        *b_hlt_Ele115;
  TBranch        *b_hlt_AK8PFJet500;
  TBranch        *b_hlt_Photon200;
  TBranch        *b_hlt_Mu37Ele27;
  TBranch        *b_hlt_Mu27Ele37;
  TBranch        *b_hlt_Mu37TkMu27;
  TBranch        *b_hlt_OldMu100;
  TBranch        *b_hlt_TkMu100;
  TBranch        *b_hlt_DoubleEle25;

  TBranch        *b_ntrigobjs;   //!
  TBranch        *b_trigobjpt;   //!
  TBranch        *b_trigobjeta;   //!
  TBranch        *b_trigobjphi;   //!
  TBranch        *b_trigobjmass;   //!                                                                       
  TBranch        *b_trigobjpdgId;   //!
  TBranch        *b_trigobjHLT;   //!
  TBranch        *b_trigobjL1;   //!
  TBranch        *b_trigobjBoth;   //!
  TBranch        *b_trigobjIhlt;   //!
  TBranch        *b_miset;   //!
  TBranch        *b_misphi;   //!
  TBranch        *b_misetsig;   //!
  TBranch        *b_sumEt;   //!
  TBranch        *b_npfjetAK8;   //!
  TBranch        *b_pfjetAK8jetID;   //!                                                    
  TBranch        *b_pfjetAK8jetID_tightlepveto;   //!
  TBranch        *b_pfjetAK8pt;   //!
  TBranch        *b_pfjetAK8y;   //!
  TBranch        *b_pfjetAK8eta; //!
  TBranch        *b_pfjetAK8phi;   //!
  TBranch        *b_pfjetAK8mass;   //!
  TBranch        *b_pfjetAK8JEC;   //!
  TBranch        *b_pfjetAK8btag_DeepCSV;   //!
  TBranch        *b_pfjetAK8DeepTag_TvsQCD;   //!
  TBranch        *b_pfjetAK8DeepTag_WvsQCD;   //!
  TBranch        *b_pfjetAK8DeepTag_ZvsQCD;   //!

  TBranch        *b_pfjetAK8CHF;   //!
  TBranch        *b_pfjetAK8NHF;   //!
  TBranch        *b_pfjetAK8CEMF;   //!
  TBranch        *b_pfjetAK8NEMF;   //!
  TBranch        *b_pfjetAK8MUF;   //!
  TBranch        *b_pfjetAK8PHF;   //!
  TBranch        *b_pfjetAK8EEF;   //!
  TBranch        *b_pfjetAK8HFHF;   //!
  TBranch        *b_pfjetAK8CHM;   //!
  TBranch        *b_pfjetAK8NHM;   //!
  TBranch        *b_pfjetAK8MUM;   //!
  TBranch        *b_pfjetAK8PHM;   //!
  TBranch        *b_pfjetAK8EEM;   //!
  TBranch        *b_pfjetAK8HFHM;   //!
  TBranch        *b_pfjetAK8Neucons;   //!
  TBranch        *b_pfjetAK8Chcons;   //!
  TBranch        *b_pfjetAK8reso;   //!
  TBranch        *b_pfjetAK8resoup;   //!
  TBranch        *b_pfjetAK8resodn;   //!
  TBranch        *b_pfjetAK8jesup_pu;   //!
  TBranch        *b_pfjetAK8jesup_rel;   //!
  TBranch        *b_pfjetAK8jesup_scale;   //!
  TBranch        *b_pfjetAK8jesup_total;   //!
  TBranch        *b_pfjetAK8jesdn_pu;   //!
  TBranch        *b_pfjetAK8jesdn_rel;   //!
  TBranch        *b_pfjetAK8jesdn_scale;   //!
  TBranch        *b_pfjetAK8jesdn_total;   //!
  TBranch        *b_pfjetAK8chrad;   //!
  TBranch        *b_pfjetAK8pTD;   //!
  TBranch        *b_pfjetAK8sdmass;   //!
  TBranch        *b_pfjetAK8tau1;   //!
  TBranch        *b_pfjetAK8tau2;   //!
  TBranch        *b_pfjetAK8tau3;   //!
  TBranch        *b_pfjetAK8sub1pt;   //!
  TBranch        *b_pfjetAK8sub1eta;   //!
  TBranch        *b_pfjetAK8sub1phi;   //!
  TBranch        *b_pfjetAK8sub1mass;   //!
  TBranch        *b_pfjetAK8sub1btag;   //!
  TBranch        *b_pfjetAK8sub1chhadfrac;   //!
  TBranch        *b_pfjetAK8sub1neuhadfrac;   //!
  TBranch        *b_pfjetAK8sub1emfrac;   //!
  TBranch        *b_pfjetAK8sub1phofrac;   //!
  TBranch        *b_pfjetAK8sub1mufrac;   //!
  TBranch        *b_pfjetAK8sub2pt;   //!
  TBranch        *b_pfjetAK8sub2eta;   //!
  TBranch        *b_pfjetAK8sub2phi;   //!
  TBranch        *b_pfjetAK8sub2mass;   //!
  TBranch        *b_pfjetAK8sub2btag;   //!
  TBranch        *b_pfjetAK8sub2chhadfrac;   //!
  TBranch        *b_pfjetAK8sub2neuhadfrac;   //!
  TBranch        *b_pfjetAK8sub2emfrac;   //!
  TBranch        *b_pfjetAK8sub2phofrac;   //!
  TBranch        *b_pfjetAK8sub2mufrac;   //!
   
  TBranch        *b_pfjetAK8subhaddiff;
  TBranch        *b_pfjetAK8subemdiff;
  TBranch        *b_pfjetAK8subptdiff;


  TBranch        *b_pfjetAK8elinsubpt;
  TBranch        *b_pfjetAK8elinsubeta;
  TBranch        *b_pfjetAK8elinsubphi;
   
  TBranch        *b_pfjetAK8elinsubmass;
  TBranch        *b_pfjetAK8elinsubjpt;
  TBranch        *b_pfjetAK8elinsubjeta;
  TBranch        *b_pfjetAK8elinsubjphi;
   
  TBranch        *b_pfjetAK8elinsubjmass;

  TBranch        *b_pfjetAK8muinsubpt;
  TBranch        *b_pfjetAK8muinsubeta;
  TBranch        *b_pfjetAK8muinsubphi;
  TBranch        *b_pfjetAK8muinsubmass;
  TBranch        *b_pfjetAK8muinsubjpt;
  TBranch        *b_pfjetAK8muinsubjeta;
  TBranch        *b_pfjetAK8muinsubjphi;
  TBranch        *b_pfjetAK8muinsubje;
  TBranch        *b_pfjetAK8muinsubjmass;
  TBranch        *b_pfjetAK8muinsubInear;
  TBranch        *b_pfjetAK8muinsubIfar;
  TBranch        *b_pfjetAK8muinsubI0;

  TBranch        *b_npfjetAK4;   //!
  TBranch        *b_pfjetAK4jetID;   //!                                                           
  TBranch        *b_pfjetAK4jetID_tightlepveto;   //!

  TBranch        *b_pfjetAK4pt;   //!
  TBranch        *b_pfjetAK4eta;   //!
  TBranch        *b_pfjetAK4y;   //!
  TBranch        *b_pfjetAK4phi;   //!
  TBranch        *b_pfjetAK4mass;   //!
  TBranch        *b_pfjetAK4JEC;   //!
   
  TBranch        *b_pfjetAK4btag_DeepCSV;   //!
  TBranch        *b_pfjetAK4btag_DeepFlav;   //!
   
  TBranch        *b_pfjetAK4reso;   //!
  TBranch        *b_pfjetAK4resoup;   //!
  TBranch        *b_pfjetAK4resodn;   //!
  TBranch        *b_pfjetAK4jesup_pu;   //!
  TBranch        *b_pfjetAK4jesup_rel;   //!
  TBranch        *b_pfjetAK4jesup_scale;   //!
  TBranch        *b_pfjetAK4jesup_total;   //!
  TBranch        *b_pfjetAK4jesdn_pu;   //!
  TBranch        *b_pfjetAK4jesdn_rel;   //!
  TBranch        *b_pfjetAK4jesdn_scale;   //!
  TBranch        *b_pfjetAK4jesdn_total;   //!
  TBranch        *b_pfjetAK4hadronflav;   //!
  TBranch        *b_pfjetAK4partonflav;   //!
  TBranch        *b_pfjetAK4qgl;   //!
  TBranch        *b_pfjetAK4PUID;   //!
  TBranch        *b_pfjetAK4GenMatch;   //!
  TBranch        *b_genmiset;   //!
  TBranch        *b_genmisphi;   //!
  TBranch        *b_genmisetsig;   //!
  TBranch        *b_ngenjetAK8;   //!
  TBranch        *b_genjetAK8pt;   //!
  TBranch        *b_genjetAK8eta;   //!
  TBranch        *b_genjetAK8phi;   //!
  TBranch        *b_genjetAK8mass;   //!
  TBranch        *b_genjetAK8sdmass;   //!
   
  TBranch        *b_ngenjetAK4;   //!
  TBranch        *b_genjetAK4pt;   //!
  TBranch        *b_genjetAK4phi;   //!
  TBranch        *b_genjetAK4mass;   //!
  TBranch        *b_genjetAK4eta;   //!
  TBranch        *b_genjetAK4hadronflav;   //!
  TBranch        *b_genjetAK4partonflav;   //!
  TBranch        *b_ngenparticles;   //!
  TBranch        *b_genpartstatus;   //!
  TBranch        *b_genpartpdg;   //!
  TBranch        *b_genpartmompdg;   //!
  TBranch        *b_genpartgrmompdg;   //!
  TBranch        *b_genpartdaugno;   //!
  TBranch        *b_genpartfromhard;   //!
  TBranch        *b_genpartfromhardbFSR;   //!
  TBranch        *b_genpartisPromptFinalState;
  TBranch        *b_genpartisLastCopyBeforeFSR;
  TBranch        *b_genpartpt;   //!
  TBranch        *b_genparteta;   //!
  TBranch        *b_genpartphi;   //!
  TBranch        *b_genpartm;   //!

  TBranch        *b_nLHEparticles;   //!
  TBranch        *b_LHEpartpdg;   //!
  TBranch        *b_LHEpartpt;   //!
  TBranch        *b_LHEparteta;   //!
  TBranch        *b_LHEpartphi;   //!
  TBranch        *b_LHEpartm;   //!

  TBranch        *b_nmuons;   //!
  TBranch        *b_muonisPF;   //!
  TBranch        *b_muonisGL;   //!
  TBranch        *b_muonisTRK;   //!
  TBranch        *b_muonisLoose;   //!
  TBranch        *b_muonisGoodGL;   //!
  TBranch        *b_muonisMed;   //!
  TBranch        *b_muonisTight;
  TBranch        *b_muonisHighPt;
  TBranch        *b_muonisHighPttrk;
  TBranch        *b_muonisMedPr;
  TBranch        *b_muonpt;   //!
  TBranch        *b_muonp;   //!
   
  TBranch        *b_muoneta;   //!
  TBranch        *b_muonphi;   //!
  TBranch        *b_muoncharge;
   
  TBranch        *b_muontrkvtx;   //!
  TBranch        *b_muondz;   //!
  TBranch        *b_muonpter;   //!
  TBranch        *b_muonchi;   //!
  TBranch        *b_muonndf;   //!
  TBranch        *b_muonecal;   //!
  TBranch        *b_muonhcal;   //!
   
  TBranch        *b_muonpfiso;   //!
   
  TBranch        *b_muonposmatch;   //!
  TBranch        *b_muontrkink;   //!
  TBranch        *b_muonsegcom;   //!

  TBranch        *b_muonpixhit;   //!
  TBranch        *b_muonmst;   //!
  TBranch        *b_muontrklay;   //!
  TBranch        *b_muonvalfrac;   //!
  TBranch        *b_muonchiso;
  TBranch        *b_muonnhiso;
  TBranch        *b_muonphiso;
  TBranch        *b_muonminisoall;
  TBranch        *b_mudxy_sv;
  TBranch        *b_muonhit;
  TBranch        *b_prefiringweight;
  TBranch        *b_prefiringweightup;
  TBranch        *b_prefiringweightdown;
  TBranch        *b_hlt_PFMET250;
  TBranch        *b_hlt_PFMET300;
  TBranch        *b_hlt_PFMET200;
  TBranch        *b_hlt_PFMET200_TypeOne;
  TBranch        *b_LHE_weight;
  TBranch        *b_nLHEScaleWeights;
  TBranch        *b_nLHEPDFWeights;
  TBranch        *b_nLHEAlpsWeights;
  TBranch        *b_Generator_x1;
  TBranch        *b_Generator_x2;
  TBranch        *b_Generator_xpdf1;
  TBranch        *b_Generator_xpdf2;
  TBranch        *b_Generator_id1;
  TBranch        *b_Generator_id2;
  TBranch        *b_Generator_scalePDF;
  TBranch        *b_Muon_corrected_pt;
  TBranch        *b_Muon_correctedUp_pt;
  TBranch        *b_Muon_correctedDown_pt;
  TBranch        *b_Electron_eccalTrkEnergyPostCorr;
  TBranch        *b_Electron_energyScaleValue;
  TBranch        *b_Electron_energyScaleUp;
  TBranch        *b_Electron_energyScaleDown;
  TBranch        *b_Electron_energySigmaValue;
  TBranch        *b_Electron_energySigmaUp;
  TBranch        *b_Electron_energySigmaDown;
  TBranch        *b_LHEScaleWeights;
  TBranch        *b_LHEPDFWeights;
  TBranch        *b_LHEAlpsWeights;
  TBranch        *b_LHEPSWeights;
  TBranch        *b_miset_CHS_UnclusEup;
  TBranch        *b_miset_CHS_UnclusEdn;
  TBranch        *b_CHSMET_phi_UnclusEup;
  TBranch        *b_CHSMET_phi_UnclusEdn;
  TBranch        *b_miset_PUPPI;
  TBranch        *b_misphi_PUPPI;
  TBranch        *b_misetsig_PUPPI;
  TBranch        *b_sumEt_PUPPI;
  TBranch        *b_miset_PUPPI_JESup;
  TBranch        *b_miset_PUPPI_JESdn;
  TBranch        *b_miset_PUPPI_JERup;
  TBranch        *b_miset_PUPPI_JERdn;
  TBranch        *b_miset_PUPPI_UnclusEup;
  TBranch        *b_miset_PUPPI_UnclusEdn;
  TBranch        *b_misphi_PUPPI_JESup;
  TBranch        *b_misphi_PUPPI_JESdn;
  TBranch        *b_misphi_PUPPI_JERup;
  TBranch        *b_misphi_PUPPI_JERdn;
  TBranch        *b_misphi_PUPPI_UnclusEup;
  TBranch        *b_misphi_PUPPI_UnclusEdn;
  TBranch        *b_nLHEPSWeights;
  TBranch        *b_nelecs;   //!
  TBranch        *b_elpt;   //!
  TBranch        *b_eldxytrk; //!
  TBranch        *b_eldztrk; //!
  TBranch        *b_eldxy_sv;
  TBranch        *b_eleta;   //!
  TBranch        *b_elphi;   //!
  TBranch        *b_elp;   //!
  TBranch        *b_ele;   //!
  TBranch        *b_elcharge;
  TBranch        *b_elmvaid;   //!
  TBranch        *b_elmvaid_Fallv2WP80;  //! 
  TBranch        *b_elmvaid_noIso;   //!
  TBranch        *b_elmvaid_Fallv2WP80_noIso; //!
   
  TBranch        *b_elhovere;   //!
  TBranch        *b_elchi;   //!
  TBranch        *b_elndf;   //!

  TBranch        *b_eletain;   //!
  TBranch        *b_elphiin;   //!

  TBranch        *b_elsupcl_eta;
  TBranch        *b_elsupcl_phi;
  TBranch        *b_elsupcl_rawE;
  TBranch        *b_elfbrem;   //!
   
  TBranch        *b_eleoverp;   //!
  TBranch        *b_elietaieta;   //!
  TBranch        *b_elmisshits;   //!
  TBranch        *b_elpfiso;   //!
   
  TBranch        *b_elsigmaieta;
  TBranch        *b_elsigmaiphi;
  TBranch        *b_elr9full;
  TBranch        *b_elsupcl_etaw;
  TBranch        *b_elsupcl_phiw;
  TBranch        *b_elhcaloverecal;
  TBranch        *b_elcloctftrkn;
  TBranch        *b_elcloctftrkchi2;
  TBranch        *b_ele1x5bye5x5;
  TBranch        *b_elnormchi2;
  TBranch        *b_elhitsmiss;
  TBranch        *b_eltrkmeasure;
  TBranch        *b_elconvtxprob;
  TBranch        *b_elecloverpout;
  TBranch        *b_elecaletrkmomentum;
  TBranch        *b_eldeltaetacltrkcalo;
  TBranch        *b_elsupcl_preshvsrawe;
  TBranch        *b_elpfisolsumphet;
  TBranch        *b_elpfisolsumchhadpt;
  TBranch        *b_elpfsiolsumneuhadet;
  TBranch        *b_pfjetAK8PNet_TvsQCD;   //!
  TBranch        *b_pfjetAK8PNet_WvsQCD;   //!
  TBranch        *b_elcutid_veto;   //!
  TBranch        *b_elcutid_loose;   //!
  TBranch        *b_elcutid_med;   //!
  TBranch        *b_elcutid_tight;   //!
  TBranch        *b_elmvaid_wpLoose;   //!
  TBranch        *b_elmvaid_noIso_wpLoose;   //!
  TBranch        *b_elsupcl_e;   //!
  TBranch        *b_genpartmomstatus;   //!

  float pfjetAK8pt_resoup[njetmxAK8], pfjetAK8mass_resoup[njetmxAK8], pfjetAK8pt_resodown[njetmxAK8], pfjetAK8mass_resodown[njetmxAK8]; 
  float pfjetAK4pt_resoup[njetmx], pfjetAK4mass_resoup[njetmx], pfjetAK4pt_resodown[njetmx], pfjetAK4mass_resodown[njetmx];
  float pfjetAK8sub1hadfrac[njetmxAK8], pfjetAK8sub2hadfrac[njetmxAK8], pfjetAK8HadF[njetmxAK8], pfjetAK8NHadF[njetmxAK8];
  //    //float puWeight, puWeightUp, puWeightDown;
 
  static const int nmaxjet = 10;
 
  double pu_rat18[100] =    {15.0025,9.25619,7.25497,5.06682,3.81133,3.00109,2.47446,2.15125,1.91921,1.76245,1.65799,1.5848,1.53433,1.49373,1.46908,1.44313,1.42046,1.40154,1.37988,1.35679,1.3361,1.31139,1.28604,1.26065,1.23868,1.21481,1.19261,1.17143,1.15268,1.13366,1.11664,1.10102,1.08379,1.06829,1.05251,1.03584,1.01745,1.00004,0.980328,0.958154,0.935215,0.910953,0.885202,0.856835,0.827565,0.797445,0.765605,0.733814,0.701484,0.668361,0.634682,0.60224,0.570207,0.537462,0.505992,0.475353,0.445231,0.417069,0.390368,0.363391,0.339587,0.316478,0.293616,0.272703,0.253076,0.23398,0.216635,0.200789,0.185823,0.171907,0.159236,0.148382,0.13732,0.127015,0.11736,0.108589,0.0994979,0.0923745,0.0861297,0.0803695,0.0741731,0.0689201,0.0636846,0.0574831,0.053418,0.0491854,0.0468157,0.0431244,0.0407697,0.0352331,0.0327207,0.0311201,0.0284833,0.0257933,0.0234082,0.0200422,0.0199332,0.0189272,0.020733,0.0166171};
  double pu_rat18_up[100] = {0,11.3701,49.1593,16.3978,10.4484,7.79227,5.70396,4.15872,3.02768,2.28549,1.82582,1.52983,1.3595,1.2554,1.19605,1.1684,1.16115,1.17185,1.18964,1.20936,1.22873,1.23491,1.23159,1.21107,1.18259,1.14644,1.11133,1.08136,1.05384,1.03331,1.01987,1.01367,1.01107,1.01298,1.01865,1.02593,1.03512,1.0447,1.05099,1.0554,1.05447,1.04466,1.02824,1.00332,0.965566,0.923431,0.871249,0.814665,0.752156,0.689408,0.624858,0.564,0.505617,0.452167,0.402,0.359344,0.321227,0.285921,0.258403,0.233682,0.210464,0.192413,0.174424,0.159861,0.146181,0.131623,0.119227,0.10899,0.0963316,0.086803,0.0773651,0.0712667,0.0629173,0.0552031,0.0481823,0.0455058,0.0376989,0.0339163,0.0298286,0.0264131,0.0255965,0.0179475,0.0169746,0.0136435,0.0117583,0.00988318,0.00674005,0.00661599,0.00316237,0.00149674,0.0010104,0.00106782,0.000384941,0.000591271,0.000423128,0.000165822,7.60044e-05,4.96232e-05,7.51979e-05,1.05862e-05};
  double pu_rat18_dn[100] = {0,15.0557,67.8751,22.3278,14.1211,10.4821,7.88069,5.86513,4.31762,3.35551,2.78627,2.40097,2.16428,2.00485,1.9056,1.85092,1.82051,1.80608,1.78719,1.75544,1.71117,1.64481,1.57234,1.49261,1.42092,1.35612,1.3043,1.26517,1.23118,1.20443,1.18302,1.16596,1.14834,1.13047,1.11055,1.08517,1.05388,1.01479,0.96502,0.907499,0.841466,0.767187,0.68971,0.610695,0.530471,0.45611,0.385995,0.32355,0.268127,0.221267,0.181416,0.149012,0.122387,0.100955,0.0832931,0.0694147,0.0579993,0.0482614,0.0406839,0.0341693,0.0284128,0.0238208,0.0196651,0.0163071,0.0134164,0.0108213,0.00875349,0.00713274,0.00561523,0.00450669,0.00357902,0.00293888,0.00231295,0.00180802,0.00140385,0.00117654,0.000861839,0.000682485,0.000525487,0.000404909,0.00033922,0.000204219,0.000164688,0.000112084,8.12391e-05,5.70485e-05,3.2298e-05,2.61592e-05,1.02574e-05,3.96059e-06,2.16985e-06,1.85204e-06,5.36884e-07,6.60936e-07,3.78607e-07,1.19189e-07,4.4536e-08,2.4673e-08,3.47283e-08,5.35281e-09};
  
 
  double pu_ratUL18[100] = {0.259884,0.223561,0.381037,0.489628,0.522993,0.507155,0.482735,0.461454,0.447782,0.439033,0.449886,0.47743,0.517608,0.564426,0.618908,0.678753,0.742106,0.808915,0.875028,0.935956,0.985031,1.02027,1.03986,1.04695,1.04456,1.04062,1.03793,1.03921,1.04531,1.05651,1.07181,1.09119,1.11281,1.13634,1.16021,1.18376,1.20516,1.22146,1.23095,1.22984,1.21792,1.1951,1.15675,1.10654,1.04226,0.972244,0.892157,0.811123,0.72623,0.645032,0.56595,0.493995,0.427171,0.370178,0.318241,0.274725,0.237453,0.203845,0.176215,0.151949,0.130797,0.112535,0.0969905,0.0826169,0.0692946,0.0586899,0.0492603,0.0407844,0.0336518,0.0273283,0.0224939,0.0180027,0.0140926,0.0113632,0.0090323,0.00707861,0.00547182,0.0042758,0.00333859,0.00185336,0.00142246,0.000955504,0.000678387,0.000578844,0.000357373,0.000265886,0.000180433,0.000167422,0.000091243,0.0000554,0.00003686,0.000029,0.00001185,0.00000481,0.00000443,0.00000657,0.0000017967,0.0000035,0.00000108};

  //TH1D *hist_obs_pu_sys[nobshist][2];
  
  
  TH2D *hist_2D_msd_deepak8=0;
  TH2D *hist_2D_bpass_flavb=0;
  TH2D *hist_2D_bpass_flavc=0;
  TH2D *hist_2D_bpass_flavq=0;
  TH2D *hist_2D_ball_flavb=0;
  TH2D *hist_2D_ball_flavc=0;
  TH2D *hist_2D_ball_flavq=0;
	
  TH1D *hist_count=0;
  //GMA define all these histogrammes. Before you put the crit on these variables, look on those first.
  static const int nprvar=45;	//GMA added these new variables (last three are for M_ll for three cases
  TH1D* hist_prvar[nprvar]={0};

  const char* prvar_name[nprvar] = {"pr_dummy","pr_nprime1","pr_anytrigger2","pr_nlepton3","pr_lepflavmatch4","pr_channeltriggerpass5","pr_leppt1fortrigthreshold6","pr_lepchargeproduct4","pr_nelec8","pr_nmuon8","pr_nak8beforecut9","pr_nak4_10","pr_dr_ak81_ak41_11","pr_dr_ak81_lep1_12","pr_dr_ak82_ak42_13","pr_dr_ak82_lep2_14","pr_nbjets","pr_ptak4", "pr_melel15", "pr_melmu15", "pr_mmumu15","pr_PFMET","pr_trigger_after","pr_extra_muon_pt8","pr_extra_electron_pt8","pr_pt_ak82_9","pr_checknan","pr_trigobjmatch7","pr_pt_ak81_9","pr_eta_ak81_9","pr_phi_ak81_9","pr_eta_ak82_9","pr_phi_ak82_9","pr_pt_lep1_3","pr_eta_lep1_3","pr_phi_lep1_3","pr_pt_lep2_3","pr_eta_lep2_3","pr_phi_lep2_3","pr_pt_elec_3","pr_eta_elec_3","pr_phi_elec_3","pr_pt_muon_3","pr_eta_muon_3","pr_phi_muon_3"};

  int prvar_bins[nprvar] = {200,    101,   65,   12,   6,    2, 60,   35,    10,   10,   10,   20, 30, 30, 30, 30, 10,  60, 60,   60,   60, 60, 65, 40, 40,60, 3, 3, 60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60,60};
  float prvar_low[nprvar] =  {-100.0,   -0.5,  -0.5, -0.5,  -0.5, -0.5, 0.0, -200, -0.5, -0.5, -0.5, -0.5,  0,  0, 0,  0, -0.5,    0,   0,   0,   0, 0, -0.5, 10, 10,0,-0.5, -0.5,0,-2.5,-3.142,-2.5,-3.142,0,-2.5,-3.142,0,-2.5,-3.142,0,-2.5,-3.142,0,-2.5,-3.142};
  float prvar_high[nprvar] ={100.0, 100.5, 64.5, 11.5, 5.5,  1.5, 900,  200,  9.5,  9.5,  9.5, 19.5, 6.5, 6.5, 6.5, 6.5, 9.5, 360, 900, 900,  900, 800, 64.5, 600, 600,1000,1.5, 1.5,1000,2.5,3.142,2.5,3.142,800,2.5,3.142,800,2.5,3.142,800,2.5,3.142,800,2.5,3.142};

  static const int ntypes=5;
  static const int ntcount=3;
  TH2D* hist_prptvar[ntypes][ntcount]={{0}}; //GMA also added on 24th JUly 2021 : Pt-eta plot of all five types for 1st, 2nd and others

  static const int nybin=120;
  double ybins[nybin+1]={10, 10.4232, 10.8643, 11.324, 11.8032, 12.3027, 12.8233, 13.366, 13.9316, 14.5211, 15.1356, 15.7761, 16.4437, 17.1396, 17.8649, 18.6209, 19.4089, 20.2302, 21.0863, 21.9786, 22.9087, 23.8781, 24.8886, 25.9418, 27.0396, 28.1838, 29.3765, 30.6196, 31.9154, 33.266, 34.6737, 36.141, 37.6704, 39.2645, 40.9261, 42.658, 44.4631, 46.3447, 48.3059, 50.3501, 52.4807, 54.7016, 57.0164, 59.4292, 61.9441, 64.5654, 67.2977, 70.1455, 73.1139, 76.2079, 79.4328, 82.7942, 86.2979, 89.9498, 93.7562, 97.7237, 101.859, 106.17, 110.662, 115.345, 120.226, 125.314, 130.617, 136.144, 141.906, 147.911, 154.17, 160.694, 167.494, 174.582, 181.97, 189.671, 197.697, 206.063, 214.783, 223.872, 233.346, 243.22, 253.513, 264.241, 275.423, 287.078, 299.226, 311.889, 325.087, 338.844, 353.183, 368.129, 383.707, 399.945, 416.869, 434.51, 452.898, 472.063, 492.04, 512.861, 534.564, 557.186, 580.764, 605.341, 630.957, 657.658, 685.488, 714.496, 744.732, 776.247, 809.096, 843.335, 879.023, 916.22, 954.993, 995.405, 1037.53, 1081.43, 1127.2, 1174.9, 1224.62, 1276.44, 1330.45, 1386.76, 1445.44};

  const char* ptprvar_name[ntypes] = {"ptpreta_ak8", "ptpreta_ak4", "ptpreta_bjet", "ptpreta_el", "ptpreta_muon"};

  static const int npr_angle=10;
  TH2D* hist_prptangle[npr_angle]={0};
  const char* pr_angle[npr_angle] = {"elec_vs_AK4_dr", "muon_vs_ak4_dr", "ak4_vs_ak8_1pr", "ak4_vs_ak8_2pr", "elec_vs_ak8_1pr", "muon_vs_ak8_1pr", "elec_vs_ak8_2pr", "muon_vs_ak8_2pr", "elec_vs_trig", "muon_vs_trig"};


  static const int noptbins = 32;
  double ptbins[noptbins+1] = {300, 390, 468,
			       507, 548, 592, 638, 686, 737, 790, 846, 905, 967, 1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000, 2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103};

  static const int nobptbins = 9;
  double bptbins[nobptbins+1] = {20, 30, 50, 70, 100, 140, 200, 300, 600, 1000};

  static const int nobetabins = 2;
  double betabins[nobetabins+1] = {0, 1.48, 2.5};

  static const int notptbins = 3;
  double tptbins[notptbins+1] = {300, 480, 600, 1200};

  float DAK8_topcut = 0.470; //0.470 ;//1% mistag rate (0.920 // 0.1% mistag)
  float deep_btag_cut = 0.2783; //0.0494//loose //(0.2770//medium) (0.7264 //tight) => These are for Autumn18
  //for UL18 => 0.0490: loose, 0.2783: medium, 0.7100: tight
  float re_cut = 0.3;
  float rt_cut = 0.7;//0.725;

  float muon_pt_cut = 20;
  float electron_pt_cut = 20;
  float lepton_pt_cut = 30;
  float AK4jet_pt_cut = 30;
  float AK8jet_pt_cut = 300;

  float absetacut = 2.5;

  static const int nre = 100;
  
  TProofOutputFile *OutFile;
  TFile *fileOut;
 
  // TTree *Tout ;
  TTree *Tnewvar;
    
  //TH1D *hist_new_var[15];
  //GMA define nhist_in and use that
  static const int nhist_in=19;
  TH1D *hist_init[nhist_in];

  static const int nobshist = 47;
  TH1D *hist_obs[nobshist];

  const char *obsnames[nobshist] = {
    "pt1","y1","mass1",  //3
    "NHad1","neuhad1","sdmass1","chrad1","subhaddiff1","tau21_1",  //9
    "DAK8_topvsQCD1","DAK8_WvsQCD1","DAK8_ZvsQCD1",   //12
    "re1","rmu1",   //14
    "haspfelectron1", //15
    "haspfmuon1","hasmatchel1","hasmatchmu1", //18
    "deltaR_l1_1","deltaR_l2_1","deltaR_bj1_1","deltaR_bj2_1",  //22
    "deltaR_1_2","pt2","y2","mass2",  //26
    "NHad2","neuhad2","sdmass2","chrad2","subhaddiff2","tau21_2",  //32
    "DAK8_topvsQCD2","DAK8_WvsQCD2","DAK8_ZvsQCD2", //35
    "re2","rmu2",  //37
    "haspfelectron2", //38
    "haspfmuon2","hasmatchel2","hasmatchmu2",  //41
    "deltaR_l1_2","deltaR_l2_2","deltaR_bj1_2","deltaR_bj2_2","tau32_1","tau32_2"  //47
  };

  double obs_low[nobshist] = {240,-2.5,0,0,0,0,-0.2,0,0,0,0,0,-1.0,-1.0,-0.5,-0.5,-0.5,-0.5,0,0,0,0,0,240,-2.5,0,0,0,0,-0.2,0,0,0,0,0,-1.0,-1.0,-0.5,-0.5,-0.5,-0.5,0,0,0,0,0,0};
  double obs_up[nobshist] = {2000,2.5,300,1,1,500,0.2,1,1,1,1,1,1,1,1.5,1.5,1.5,1.5,5,5,5,5,5,2000,2.5,300,1,1,300,0.2,1,1,1,1,1,1,1,1.5,1.5,1.5,1.5,5,5,5,5,1,1};


  int obs_nbins[nobshist] = {30,30,30,30,30,30,30,30,60,30,30,30,60,60,2,2,2,2,30,30,30,30
			     ,30,30,30,30,30,30,30,30,30,60,30,30,30,60,60,2,2,2,2,30,30,30,30,30,30};

  static const int nhistbtag=14;
  TH1D *hist_obs_btag[nhistbtag];
  const char *obsnames_btag[nhistbtag] = {"NHad1_med_btagcut","neuhad1_med_btagcut","sdmass1_med_btagcut","chrad1_med_btagcut","subhaddiff1_med_btagcut","tau21_1_med_btagcut","tau32_1_med_btagcut","NHad2_med_btagcut","neuhad2_med_btagcut","sdmass2_med_btagcut","chrad2_med_btagcut","subhaddiff2_med_btagcut","tau21_2_med_btagcut","tau32_2_med_btagcut"};

  double obs_low_btag[nhistbtag] = {0,0,0,-0.15,0,0,0,0,0,0,-0.15,0,0,0};
  double obs_up_btag[nhistbtag] = {1,1,300,0.15,1,1,1,1,1,300,0.15,1,1,1};
  int obs_nbins_btag[nhistbtag] = {30,30,30,60,30,30,30,30,30,30,60,30,30,30};

  static const int nhistgenmatch=14;
  TH1D *hist_obs_genmatch[nhistgenmatch];
  const char *obsnames_genmatch[nhistgenmatch] = {"NHad1_genmatch","neuhad1_genmatch","sdmass1_genmatch","chrad1_genmatch","subhaddiff1_genmatch","tau21_1_genmatch","tau32_1_genmatch","NHad2_genmatch","neuhad2_genmatch","sdmass2_genmatch","chrad2_genmatch","subhaddiff2_genmatch","tau21_2_genmatch","tau32_2_genmatch"};

  double obs_low_genmatch[nhistgenmatch] = {0,0,0,-0.3,0,0,0,0,0,0,-0.3,0,0,0};
  double obs_up_genmatch[nhistgenmatch] = {1,1,300,0.3,1,1,1,1,1,300,0.3,1,1,1};

  int obs_nbins_genmatch[nhistgenmatch] = {30,30,30,60,30,30,30,30,30,30,60,30,30,30};
  TH1D *hist_init_pu_sys[18][2]={{0}};

  static const int ntreevarhist = 47;
  TH1D *hist_treevar[ntreevarhist];

  const char *treevarnames[ntreevarhist] = {"M_l1l2","rat_l2pt_l1pt","deltaPhi_l1l2",  //3
					    "l1pt_nearjet","l2pt_nearjet","met_pt","met_phi",                                //7
					    "M_bl1","M_bl2",                                                                 //9
					    "delta_phil1_met","delta_phil2_met","delta_phibl1_met","delta_phibl2_met",       //13
					    "rat_metpt_ak4pt","rat_metpt_ak8pt","rat_metpt_eventHT",                         //16
					    "mt_of_l1met","mt_of_l2met","no_ak4jets","no_ak4bjets","no_ak8jets","EventHT",   //22
					    "extra_ak4j","ptsum_extra_ak4","extra_ak4pt","extra_ak4mass","extra_ak4jqgl",    //27
					    "extra_ak4jdeepb","rat_extra_ak4jpt_lpt","ak81deep_tvsqcd","ak81deep_wvsqcd",    //31
					    "ak82deep_tvsqcd","ak82deep_wvsqcd",                                             //33
					    "delta_phibl1bl2","deltaR_l1l2","deltaR_l1j1","deltaR_l2j1",                     //37
					    "deltaR_l1j2","deltaR_l2j2","j1_btag_sc","j2_btag_sc",                           //41
					    "dirgltrthr","dirglthrmin","response_ak81","response_ak82","deltaR_l1l2","response_eventclass"};       //46
  double treevar_low[ntreevarhist] = {100,0,-3.15,0,0,50,-3.15,0,0,-3.15,-3.15,-3.15,-3.15,0,0,0,0,0,-0.5,-0.5,1.5,400,-0.5,0,0,0,0,0,0,0,0,0,0,-3.15,1,0,1,1,0,0,0,0,0,-1,-1,1,-1};
  double treevar_up[ntreevarhist] = {1000,8,3.15,250,250,400,3.15,200,200,3.15,3.15,3.15,3.15,2,0.8,0.3,1000,1000,10.5,5.5,5.5,4000,7.5,500,300,200,1,1,7,1,1,1,1,3.15,5,1,5,5,1,1,1,0.25,0.5,1,1,5,1};
  int treevar_nbins[ntreevarhist] = {30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,11,6,4,30,8,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30};

  const char *initnames[nhist_in] = {"matchN_l1","matchN_l2","N_PV_sel","NJets_AK4","NBJets_AK4","NJets_AK8","mll","l1pt","l1eta","l1phi","l2pt","l2eta","l2phi","ak41jetpt","ak41jeteta","ak41jetphi","ak42jetpt","ak42jeteta","ak42jetphi"};

  const char *titlenames[nhist_in] = {"trigger matching for mu","trigger matching for el","# of Primary Vertices","# of AK4 jets","# of b tagged AK4 jets","# of AK8 jets","mll","l1pt","l1eta","l1phi","l2pt","l2eta","l2phi","pT of leading AK4","eta of leading AK4","phi of leading AK4","pT of subleading AK4","eta of subleading AK4","phi of subleading AK4"};

  //double new_var_low[9] = {0.0,0.0,-5.0,0.0,0.0,0.0,-2.5,-5.0,-5.0};//,0.0,0.0,0.0,0.0,0.0,0.0,0.0,300.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  //double new_var_up[9] = {500.0,2.0,5.0,1000.0,1000.0,500.0,2.5,5.0,5.0};//1000,3.15,3.15,3.15,500,500,3000,500,3.15,3.15,3.15,3.15,500,500,5.};
  //int new_var_nbins[9] = {25,25,50,50,50,25,25,50,50};//30,25,25,25,25,25,30,25,25,25,25,25,25,25,25};

  double ini_low[nhist_in] = {0.0,0.0,-0.1,-0.5,-0.5,-0.5,84.5,15.0,-2.5,-3.15,15.0,-2.5,-3.15,6.0,-2.5,-3.15,6.0,-2.5,-3.15};
  double ini_up[nhist_in] = {1.0,1.0,99.9,10.5,10.5,10.5,1000.0,900.0,2.5,3.15,900.0,2.5,3.15,1400.0,2.5,3.15,1400.0,2.5,3.15};
  int ini_nbins[nhist_in] = {60,60,100,11,11,11,60,60,30,30,60,30,30,60,30,30,60,30,30};

  TH1D *hist_npv;
  TH1D *hist_npv_nopuwt;

  const static int nhistdeltaR=8;
  TH1D *hist_genmatch_deltaR[nhistdeltaR];
  const char *names_genmatch_deltaR[nhistdeltaR] = {"deltaR_genlep_ak8_1","deltaR_genb_ak8_1","deltaR_genblep_ak8_1","deltaR_gentop_ak8_1","deltaR_genlep_ak8_2","deltaR_genb_ak8_2","deltaR_genblep_ak8_2","deltaR_gentop_ak8_2"};  /// Used in the code

  const char *titles_genmatch_deltaR[nhistdeltaR] = {"deltaR(gen lep, leading lep jet)","deltaR(gen b, leading AK4 jet)","deltaR(gen (b +lep), leading AK8 jet)","deltaR(gen top, leading AK8 jet)","deltaR(gen lep, sub leading lep jet)","deltaR(gen b, sub leading AK4 jet)","deltaR(gen (b +lep), sub leading AK8 jet)","deltaR(gen top, sub leading AK8 jet)"};

  const static int nhistdrsc=2;
  TH1D *hist_genmatch_deltaR_score[nhistdrsc];
  const char *names_genmatch_deltaR_score[nhistdrsc] = {"genmatch_score_ak8_1","genmatch_score_ak8_2"};

  const char *titles_genmatch_deltaR_score[nhistdrsc] = {"genmatching score of leading AK8 jet","genmatching score of sub leading AK8 jet"};

  const static int ncutflow=16;
  TH1D *hist_met_cutflow[ncutflow];

  TH1D *hist_elidvar[7];
  const char *names_elidvar[7] = {"chs_met_eltop","dirgltrthr_eltop","dirglthrmin_eltop","deltaphi_eltop","sdmasstagjet_eltop","eventHT_eltop","rat_metpt_eventHT_eltop"};
  double elidvar_low[7] = {0,0,0,0,0,0,0};
  double elidvar_up[7] = {900,0.3,0.6,3.142,250,2500,0.5};
  int elidvar_nbins[7] = {30,30,30,30,30,30,30};

  TH1D *hist_muidvar[7];
  const char *names_muidvar[7] = {"chs_met_mutop","dirgltrthr_mutop","dirglthrmin_mutop","deltaphi_mutop","sdmasstagjet_mutop","eventHT_mutop","rat_metpt_eventHT_mutop"};
  double muidvar_low[7] = {0,0,0};
  double muidvar_up[7] = {900,0.3,0.6,3.142,250,2500,0.5};
  int muidvar_nbins[7] = {30,30,30,30,30,30,30};

  TH2D *hist_2d_deltaR_vsbtagsc[2];
  TH2D *hist_2d_pt_vsbtagsc[4];
  TH2D *hist_2d_deltaR_vspt[2];

  TH1D *hist_npu_vert_true,*hist_npu_vert;
  TH1D *hist_count_lep[2];
  TH1D *hist_count_top,*hist_delptbypt[3];
  TH2D *hist_2d_pt_genlepvsb1, *hist_2d_pt_genlepvsb2, *hist_2d_pt_genlepvsb3, *hist_2d_pt_genlepvsb4, *hist_2d_pt_genlepvsb5, *hist_2d_pt_gentopvsgentop;

    const static int nptbin = 4,netabin = 2;
    double ptbinarray[nptbin+1] = {300.0, 400.0, 500.0, 650.0, 1300.0};
    double etabinarray[netabin+1] = {0.0, 1.48, 2.5};

	const static int nchannel = 2;
	TH2D *hist2d_toptag_sf[5][nchannel];
	TH1D *hist_ptak8_sf[5][nchannel], *hist_etaak8_sf[5][nchannel], *hist_sdmassak8_sf[5][nchannel],*hist_sdmasstagak8_sf[5][nchannel],*hist_pnettagak8_sf[5][nchannel], *hist_probescoreak8_sf[5][nchannel], *hist_ak8score_sf[nchannel];

	const static int nhistmcdata = 8;

	TH1D *hist_ak8pt_mcdata_allmatch[nchannel][2], *hist_ak8eta_mcdata_allmatch[nchannel][2], *hist_ak8sdmass_mcdata_allmatch[nchannel][2], *hist_ak8score_mcdata_allmatch[nchannel][2];

	TH1D *hist_ak8pt_mcdata_leporbmatch[nchannel][2], *hist_ak8eta_mcdata_leporbmatch[nchannel][2], *hist_ak8sdmass_mcdata_leporbmatch[nchannel][2], *hist_ak8score_mcdata_leporbmatch[nchannel][2];

	TH1D *hist_ak8pt_mcdata_nomatch[nchannel][2], *hist_ak8eta_mcdata_nomatch[nchannel][2], *hist_ak8sdmass_mcdata_nomatch[nchannel][2], *hist_ak8score_mcdata_nomatch[nchannel][2];

	TH1D *hist_ak8pt_mcdata_up[nhistmcdata][nchannel], *hist_ak8eta_mcdata_up[nhistmcdata][nchannel], *hist_ak8sdmass_mcdata_up[nhistmcdata][nchannel], *hist_ak8score_mcdata_up[nhistmcdata][nchannel];

	TH1D *hist_ak8pt_mcdata_dn[nhistmcdata][nchannel], *hist_ak8eta_mcdata_dn[nhistmcdata][nchannel], *hist_ak8sdmass_mcdata_dn[nhistmcdata][nchannel], *hist_ak8score_mcdata_dn[nhistmcdata][nchannel];

	TH1D *hist_ak8sdmass_sf_allmatch_pass[nptbin][netabin][nchannel],*hist_ak8sdmass_sf_allmatch_fail[nptbin][netabin][nchannel];
	TH1D *hist_ak8sdmass_sf_leporbmatch_pass[nptbin][netabin][nchannel],*hist_ak8sdmass_sf_leporbmatch_fail[nptbin][netabin][nchannel];
	TH1D *hist_ak8sdmass_sf_nomatch_pass[nptbin][netabin][nchannel],*hist_ak8sdmass_sf_nomatch_fail[nptbin][netabin][nchannel];

	const char *mcdatanames[nhistmcdata] = {"jes","jer","pileup","hadtoptag","pdfscale","btag","met"};

	float wp_leptop_tagger[nchannel] = {-1};
	
  float in_pfjetAK8NHadF;
  float in_pfjetAK8neunhadfrac;
  float in_pfjetAK8subhaddiff;
  float in_pfjetAK8tau21;
  float in_pfjetAK8chrad;
  float in_pfjetAK8sdmass;
  float in_pfjetAK8eldxy_sv;
  float in_pfjetAK8matchedeldxy_sv;
  float in_pfjetAK8matchedelcleta;
  float in_pfjetAK8matchedelpt;
  float in_pfjetAK8matchedelsigmaieta;
  float in_pfjetAK8matchedelsigmaiphi;
  float in_pfjetAK8matchedelr9full;
  float in_pfjetAK8matchedelsupcl_etaw;
  float in_pfjetAK8matchedelsupcl_phiw;
  float in_pfjetAK8matchedelhcaloverecal;
  float in_pfjetAK8matchedelcloctftrkn;
  float in_pfjetAK8matchedelcloctftrkchi2;
  float in_pfjetAK8matchedele1x5bye5x5;
  float in_pfjetAK8matchedelnormchi2;
  float in_pfjetAK8matchedelhitsmiss;
  float in_pfjetAK8matchedeltrkmeasure;
  float in_pfjetAK8matchedelecloverpout;
  float in_pfjetAK8matchedelecaletrkmomentum;
  float in_pfjetAK8matchedeldeltaetacltrkcalo;
  float in_pfjetAK8matchedelsupcl_preshvsrawe;
  float in_pfjetAK8matchedelpfisolsumphet;
  float in_pfjetAK8matchedelpfisolsumchhadpt;
  float in_pfjetAK8matchedelpfisolsumneuhadet;
  float in_pfjetAK8matchedeletain;
  float in_pfjetAK8matchedelphiin;
  float in_pfjetAK8matchedelfbrem;
  float in_pfjetAK8matchedeleoverp;
  float in_pfjetAK8matchedelhovere;
  float in_pfjetAK8matchedelRho;
	
  float in_mupfjetAK8NHadF;
  float in_mupfjetAK8neunhadfrac;
  float in_mupfjetAK8subhaddiff;
  float in_mupfjetAK8tau21;
  float in_mupfjetAK8chrad;
  float in_mupfjetAK8sdmass;
  float in_pfjetAK8matchedmuonchi;
  float in_pfjetAK8matchedmuonposmatch;
  float in_pfjetAK8matchedmuontrkink;
  float in_pfjetAK8matchedmuonsegcom;
  float in_pfjetAK8matchedmuonhit;
  float in_pfjetAK8matchedmuonmst;
  float in_pfjetAK8matchedmuontrkvtx;
  float in_pfjetAK8matchedmuondz;
  float in_pfjetAK8matchedmuonpixhit;
  float in_pfjetAK8matchedmuontrklay;
  float in_pfjetAK8matchedmuonvalfrac;
  float in_pfjetAK8muinsubptrat;
  float in_pfjetAK8muinsubmassrat;
  float in_pfjetAK8muinsubinvmass;
  float in_pfjetAK8muinsubIfarbyI0;
  float in_pfjetAK8muinsubInearbyI0;
	
  TMVA::Reader *reader1;
  TMVA::Reader *reader2;
  TMVA::Reader *reader3;
  TMVA::Reader *reader4;
    
  TString testdir = "/home/deroy/t3store3/CMSSW_10_5_0/src/BDTResponse_validator/Analysis/newvar_sv/Signal/";
  //TString dir = "/home/deroy/t3store3/Muon_MuEl/";
  //TString dir = "/home/rsaxena/t3store3/Muon_MuEl/";
  TString dir = "/home/ritik/Desktop/code2/";
  //TString weightfile1 = testdir + TString("TMVAClassification_BDTG_elIDvarv3.weights.xml");
	
  //TString weightfile1 = dir + TString("TMVAClassification_BDTG_elIDvar_Jan2021Corr_TTbarUL18.weights.xml");
  TString weightfile1 = dir + TString("TMVAClassification_BDTG_for_eltop_pt300_weights.xml"/*"TMVAClassification_BDTG_for_eltop_pt300_withqcdbkg_weights.xml"*/);
  TString weightfile2;
  TString weightfile3 = dir + TString("TMVAClassification_only7varsnomatchel_BDTG.weights.xml");
  TString weightfile4 = dir + TString("TMVAClassification_BDTG_for_mutop_pt300_weights.xml"/*"TMVAClassification_BDTG_for_mutop_pt300_withqcdbkg_weights.xml"*/);

  TFile *electron_reco_sf_file = new TFile(dir + TString("egammaEffi_ptAbove20.txt_EGM2D_UL2018.root"));
  TH2F *h_SF_el_reco = (TH2F*)electron_reco_sf_file->Get("EGamma_SF2D");
  TH2F *h_SF_statData_el_reco = (TH2F*)electron_reco_sf_file->Get("statData");
  TH2F *h_SF_statMC_el_reco = (TH2F*)electron_reco_sf_file->Get("statMC");
  TH2F *h_SF_altBkgModel_el_reco = (TH2F*)electron_reco_sf_file->Get("altBkgModel");
  TH2F *h_SF_altSignalModel_el_reco = (TH2F*)electron_reco_sf_file->Get("altSignalModel");
  TH2F *h_SF_altMCEff_el_reco = (TH2F*)electron_reco_sf_file->Get("altMCEff");
  TH2F *h_SF_altTagSelection_el_reco = (TH2F*)electron_reco_sf_file->Get("altTagSelection");

  TFile *electron_loose_sf_file = new TFile(dir + TString("2018_ElectronLoose.root"));
  TH2F *h_SF_el_loose = (TH2F*)electron_reco_sf_file->Get("EGamma_SF2D");
 /* TH2F *h_SF_statData_el_loose = (TH2F*)electron_reco_sf_file->Get("statData");
  TH2F *h_SF_statMC_el_loose = (TH2F*)electron_reco_sf_file->Get("statMC");
  TH2F *h_SF_altBkgModel_el_loose = (TH2F*)electron_reco_sf_file->Get("altBkgModel");
  TH2F *h_SF_altSignalModel_el_loose = (TH2F*)electron_reco_sf_file->Get("altSignalModel");
  TH2F *h_SF_altMCEff_el_loose = (TH2F*)electron_reco_sf_file->Get("altMCEff");
  TH2F *h_SF_altTagSelection_el_loose = (TH2F*)electron_reco_sf_file->Get("altTagSelection");*/

  TFile *electron_tightid_sf_file = new TFile(dir + TString("egammaEffi.txt_Ele_wp90noiso_EGM2D.root"));
  TH2F *h_SF_el_tight = (TH2F*)electron_tightid_sf_file->Get("EGamma_SF2D");
  TH2F *h_SF_statData_el_tight = (TH2F*)electron_tightid_sf_file->Get("statData");
  TH2F *h_SF_statMC_el_tight = (TH2F*)electron_tightid_sf_file->Get("statMC");
  TH2F *h_SF_altBkgModel_el_tight = (TH2F*)electron_tightid_sf_file->Get("altBkgModel");
  TH2F *h_SF_altSignalModel_el_tight = (TH2F*)electron_tightid_sf_file->Get("altSignalModel");
  TH2F *h_SF_altMCEff_el_tight = (TH2F*)electron_tightid_sf_file->Get("altMCEff");
  TH2F *h_SF_altTagSelection_el_tight = (TH2F*)electron_tightid_sf_file->Get("altTagSelection");

  TFile *muon_reco_sf_file = new TFile(dir + TString("Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root"));
  TH2F *h_SF_mu_reco = (TH2F*)muon_reco_sf_file->Get("NUM_TrackerMuons_DEN_genTracks");

  TFile *muon_tightid_sf_file = new TFile(dir + TString("Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root"));
  TH2F *h_SF_mu_tight = (TH2F*)muon_tightid_sf_file->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt");
  TH2F *h_SF_stat_mu_tight = (TH2F*)muon_tightid_sf_file->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_stat");
  TH2F *h_SF_sys_mu_tight = (TH2F*)muon_tightid_sf_file->Get("NUM_TightID_DEN_TrackerMuons_abseta_pt_syst");

  TH2F *h_SF_mu_loose = (TH2F*)muon_tightid_sf_file->Get("NUM_LooseID_DEN_TrackerMuons_abseta_pt");
  TH2F *h_SF_stat_mu_loose = (TH2F*)muon_tightid_sf_file->Get("NUM_LooseID_DEN_TrackerMuons_abseta_pt_stat");
  TH2F *h_SF_sys_mu_loose = (TH2F*)muon_tightid_sf_file->Get("NUM_LooseID_DEN_TrackerMuons_abseta_pt_syst");

  TFile *file_pu_ratio = new TFile(dir + TString("RatioPileup-UL2018-100bins.root"));
  TH1F *h_pileup_data = (TH1F*)file_pu_ratio->Get("pileup_weight");
  TH1F *h_pileup_data_plus = (TH1F*)file_pu_ratio->Get("pileup_plus_weight");
  TH1F *h_pileup_data_minus = (TH1F*)file_pu_ratio->Get("pileup_minus_weight");

  BTagCalibration calib_deepcsv, calib_deepflav;
  BTagCalibrationReader reader_deepcsv, reader_deepflav_loose,reader_deepflav_med,reader_deepflav_tight;
  //float ptcut = 200; //Need to ask Suman how ptcut of ljet is then set if it is not declared here
  
  top_tagger_sf(TTree * /*tree*/ =0) : fChain(0) { }
  virtual ~top_tagger_sf() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    Begin(TTree *tree);
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify();

  virtual float*  Muon_SF(float pt, float eta, string id);
  virtual float*  Electron_SF(float pt, float eta, string id);
  virtual float*  Get_PU_Weights(int npu);

  // declare object selection functions // 
  virtual void    getmuons(std::vector<Muon> &vmuons, float ptcut, float etacut, int maxsize, float dxy_cut, float dz_cut);
  virtual void    getelectrons(std::vector<Electron> &velectrons, float ptcut, float etacut, int maxsize, float dxy_cut, float dz_cut);
  virtual void    getLeptons(std::vector<Lepton> &vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, float pt_cut);
  virtual void    convertintolep(Muon vmuon, Lepton &vlepton);
  virtual void    convertintolep(Electron velectron, Lepton &vlepton);
  virtual void    getAK4jets(std::vector<AK4Jet> &Jets, float ptcut, float etacut, bool isMC, int maxsize);
  virtual void    getAK8jets(std::vector<AK8Jet> &LJets, float ptcut, float etacut, bool isMC, int maxsize);
  virtual void    getAK4genjets(std::vector<AK4GenJet> &Jets, float ptcut, float etacut, int maxsize);
  //virtual void    LeptonJet_cleaning(std::vector <Electron> velectrons,std::vector <Muon> vmuons, float dR_cut);
  virtual void    getPartons(std::vector<GenParton> &GenPartons, int maxsize);
  virtual void    getLHEParticles(std::vector<LHEparticle> &lheparticles, int maxsize);
  virtual void    getLHETops(std::vector<GenParton> &LHETops, std::vector<GenParton> GenPartons);
  virtual void    getGENTops(vector<TopQuark> &gentops, vector<GenParton> genpartons);
  virtual void    TopAssignment_toJet(std::vector<AK8Jet> &LJets, std::vector<GenParton> lhetops, std::vector<TopQuark> gentops);
  virtual void    AssignGen(std::vector<AK8Jet> &LJets, std::vector<GenParton> GenPartons);
  virtual bool    isBJet(AK4Jet jet, float btag_cut);
  virtual void    ReadTagger(AK8Jet &LJet, std::vector<Lepton> vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, TMVA::Reader *reader_electron, TMVA::Reader *reader_muon);
  virtual void    Match_trigger(vector<Single_Trigger> vsinglelep_trig,vector<Double_Trigger> vdoublelep_trig,
				vector<std::pair<int,TLorentzVector> > TrigRefObj,
				Lepton lepcand_1, Lepton lepcand_2, vector<AK4Jet> Jets,
				bool &trig_threshold_pass,
				bool &trig_matching_pass,
				vector<TH1D*> &hist_init,
				vector<TH2D*> &hist_2d
				);

  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    SlaveTerminate();
  virtual void    Terminate();
  
  ClassDef(top_tagger_sf,0);
};

#endif

#ifdef top_tagger_sf_cxx
void top_tagger_sf::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

   
  fChain->SetBranchAddress("irun", &irun, &b_irun);
  fChain->SetBranchAddress("ilumi", &ilumi, &b_ilumi);
  fChain->SetBranchAddress("ievt", &ievt, &b_ievt);
  fChain->SetBranchAddress("nprim", &nprim, &b_nprim);
  fChain->SetBranchAddress("nprimi", &nprimi, &b_nprimi);
  fChain->SetBranchAddress("Rho", &Rho, &b_Rho);
  fChain->SetBranchAddress("Generator_weight", &event_weight, &b_event_weight);
  fChain->SetBranchAddress("qscale", &qscale, &b_qscale);
  fChain->SetBranchAddress("npu_vert", &npu_vert, &b_npu_vert);
  fChain->SetBranchAddress("npu_vert_true", &npu_vert_true, &b_npu_vert_true);
  fChain->SetBranchAddress("trig_value", &trig_value, &b_trig_value);
  fChain->SetBranchAddress("hlt_IsoMu24",&hlt_IsoMu24,&b_hlt_IsoMu24);
  fChain->SetBranchAddress("hlt_Mu50",&hlt_Mu50,&b_hlt_Mu50);
  fChain->SetBranchAddress("hlt_Ele50_PFJet165",&hlt_Ele50_PFJet165,&b_hlt_Ele50_PFJet165);
  fChain->SetBranchAddress("hlt_Ele115",&hlt_Ele115,&b_hlt_Ele115);
  fChain->SetBranchAddress("hlt_AK8PFJet500", &hlt_AK8PFJet500, &b_hlt_AK8PFJet500);
  fChain->SetBranchAddress("hlt_Photon200", &hlt_Photon200, &b_hlt_Photon200);
  fChain->SetBranchAddress("hlt_Mu37Ele27", &hlt_Mu37Ele27, &b_hlt_Mu37Ele27);
  fChain->SetBranchAddress("hlt_Mu27Ele37", &hlt_Mu27Ele37, &b_hlt_Mu27Ele37);
  fChain->SetBranchAddress("hlt_Mu37TkMu27", &hlt_Mu37TkMu27, &b_hlt_Mu37TkMu27);
  fChain->SetBranchAddress("hlt_OldMu100", &hlt_OldMu100, &b_hlt_OldMu100);
  fChain->SetBranchAddress("hlt_TkMu100", &hlt_TkMu100, &b_hlt_TkMu100);
  fChain->SetBranchAddress("hlt_DoubleEle25", &hlt_DoubleEle25, &b_hlt_DoubleEle25);
  fChain->SetBranchAddress("hlt_PFMET250",&hlt_PFMET250, &b_hlt_PFMET250);
  fChain->SetBranchAddress("hlt_PFMET300",&hlt_PFMET300, &b_hlt_PFMET300);
  fChain->SetBranchAddress("hlt_PFMET200",&hlt_PFMET200, &b_hlt_PFMET200);
  fChain->SetBranchAddress("hlt_PFMET200_TypeOne",&hlt_PFMET200_TypeOne, &b_hlt_PFMET200_TypeOne);
  fChain->SetBranchAddress("ntrigobjs", &ntrigobjs, &b_ntrigobjs);
  fChain->SetBranchAddress("trigobjpt", trigobjpt, &b_trigobjpt);
  fChain->SetBranchAddress("trigobjeta", trigobjeta, &b_trigobjeta);
  fChain->SetBranchAddress("trigobjphi", trigobjphi, &b_trigobjphi);
  fChain->SetBranchAddress("trigobjmass", trigobjmass, &b_trigobjmass);
  fChain->SetBranchAddress("trigobjpdgId", trigobjpdgId, &b_trigobjpdgId);
  fChain->SetBranchAddress("trigobjHLT", trigobjHLT, &b_trigobjHLT);
  fChain->SetBranchAddress("trigobjL1", trigobjL1, &b_trigobjL1);
  fChain->SetBranchAddress("trigobjBoth", trigobjBoth, &b_trigobjBoth);
  fChain->SetBranchAddress("trigobjIhlt", trigobjIhlt, &b_trigobjIhlt);
  fChain->SetBranchAddress("prefiringweight",&prefiringweight, &b_prefiringweight);
  fChain->SetBranchAddress("prefiringweightup",&prefiringweightup, &b_prefiringweightup);
  fChain->SetBranchAddress("prefiringweightdown",&prefiringweightdown, &b_prefiringweightdown);

  fChain->SetBranchAddress("CHSMET_pt",&chs_met, &b_miset);
  fChain->SetBranchAddress("CHSMET_phi",&chs_misphi, &b_misphi);
  fChain->SetBranchAddress("CHSMET_sig",&chs_misetsig, &b_misetsig);
  fChain->SetBranchAddress("CHSMET_sumEt",&chs_sumEt, &b_sumEt);

  fChain->SetBranchAddress("CHSMET_pt_UnclusEup",&miset_UnclusEup, &b_miset_CHS_UnclusEup);
  fChain->SetBranchAddress("CHSMET_pt_UnclusEdn",&miset_UnclusEdn, &b_miset_CHS_UnclusEdn);
  fChain->SetBranchAddress("CHSMET_phi_UnclusEup",&misphi_UnclusEup, &b_CHSMET_phi_UnclusEup);
  fChain->SetBranchAddress("CHSMET_phi_UnclusEdn",&misphi_UnclusEdn, &b_CHSMET_phi_UnclusEdn);

  fChain->SetBranchAddress("PuppiMET_pt",&miset_PUPPI, &b_miset_PUPPI);
  fChain->SetBranchAddress("PuppiMET_phi",&misphi_PUPPI, &b_misphi_PUPPI);
  fChain->SetBranchAddress("PuppiMET_sig",&misetsig_PUPPI, &b_misetsig_PUPPI);
  fChain->SetBranchAddress("PuppiMET_sumEt",&sumEt_PUPPI, &b_sumEt_PUPPI);

  fChain->SetBranchAddress("PuppiMET_pt_JESup",&miset_PUPPI_JESup, &b_miset_PUPPI_JESup);
  fChain->SetBranchAddress("PuppiMET_pt_JESdn",&miset_PUPPI_JESdn, &b_miset_PUPPI_JESdn);
  fChain->SetBranchAddress("PuppiMET_pt_JERup",&miset_PUPPI_JERup, &b_miset_PUPPI_JERup);
  fChain->SetBranchAddress("PuppiMET_pt_JERdn",&miset_PUPPI_JERdn, &b_miset_PUPPI_JERdn);
  fChain->SetBranchAddress("PuppiMET_pt_UnclusEup",&miset_PUPPI_UnclusEup, &b_miset_PUPPI_UnclusEup);
  fChain->SetBranchAddress("PuppiMET_pt_UnclusEdn",&miset_PUPPI_UnclusEdn, &b_miset_PUPPI_UnclusEdn);
  fChain->SetBranchAddress("PuppiMET_phi_JESup",&misphi_PUPPI_JESup, &b_misphi_PUPPI_JESup);
  fChain->SetBranchAddress("PuppiMET_phi_JESdn",&misphi_PUPPI_JESdn, &b_misphi_PUPPI_JESdn);
  fChain->SetBranchAddress("PuppiMET_phi_JERup",&misphi_PUPPI_JERup, &b_misphi_PUPPI_JERup);
  fChain->SetBranchAddress("PuppiMET_phi_JERdn",&misphi_PUPPI_JERdn, &b_misphi_PUPPI_JERdn);
  fChain->SetBranchAddress("PuppiMET_phi_UnclusEup",&misphi_PUPPI_UnclusEup, &b_misphi_PUPPI_UnclusEup);
  fChain->SetBranchAddress("PuppiMET_phi_UnclusEdn",&misphi_PUPPI_UnclusEdn, &b_misphi_PUPPI_UnclusEdn);

  fChain->SetBranchAddress("npfjetAK8", &npfjetAK8, &b_npfjetAK8);
  fChain->SetBranchAddress("pfjetAK8jetID_tightlepveto", pfjetAK8jetID_tightlepveto, &b_pfjetAK8jetID_tightlepveto);
  fChain->SetBranchAddress("pfjetAK8jetID", pfjetAK8jetID, &b_pfjetAK8jetID);
  fChain->SetBranchAddress("pfjetAK8pt", pfjetAK8pt, &b_pfjetAK8pt);
  fChain->SetBranchAddress("pfjetAK8y", pfjetAK8y, &b_pfjetAK8y);
  fChain->SetBranchAddress("pfjetAK8eta", pfjetAK8eta, &b_pfjetAK8eta);
  fChain->SetBranchAddress("pfjetAK8phi", pfjetAK8phi, &b_pfjetAK8phi);
  fChain->SetBranchAddress("pfjetAK8mass", pfjetAK8mass, &b_pfjetAK8mass);
  fChain->SetBranchAddress("pfjetAK8JEC", pfjetAK8JEC, &b_pfjetAK8JEC);
  fChain->SetBranchAddress("pfjetAK8btag_DeepCSV", pfjetAK8btag_DeepCSV, &b_pfjetAK8btag_DeepCSV);

  fChain->SetBranchAddress("pfjetAK8DeepTag_DAK8_TvsQCD", pfjetAK8DeepTag_TvsQCD, &b_pfjetAK8DeepTag_TvsQCD);
  fChain->SetBranchAddress("pfjetAK8DeepTag_DAK8_WvsQCD", pfjetAK8DeepTag_WvsQCD, &b_pfjetAK8DeepTag_WvsQCD);
  fChain->SetBranchAddress("pfjetAK8DeepTag_DAK8_ZvsQCD", pfjetAK8DeepTag_ZvsQCD, &b_pfjetAK8DeepTag_ZvsQCD);
  fChain->SetBranchAddress("pfjetAK8DeepTag_PNet_TvsQCD", pfjetAK8PNet_TvsQCD, &b_pfjetAK8PNet_TvsQCD);
  fChain->SetBranchAddress("pfjetAK8DeepTag_PNet_WvsQCD", pfjetAK8PNet_WvsQCD, &b_pfjetAK8PNet_WvsQCD);
  fChain->SetBranchAddress("pfjetAK8CHF", pfjetAK8CHF, &b_pfjetAK8CHF);
  fChain->SetBranchAddress("pfjetAK8NHF", pfjetAK8NHF, &b_pfjetAK8NHF);
  fChain->SetBranchAddress("pfjetAK8CEMF", pfjetAK8CEMF, &b_pfjetAK8CEMF);
  fChain->SetBranchAddress("pfjetAK8NEMF", pfjetAK8NEMF, &b_pfjetAK8NEMF);
  fChain->SetBranchAddress("pfjetAK8MUF", pfjetAK8MUF, &b_pfjetAK8MUF);
  fChain->SetBranchAddress("pfjetAK8PHF", pfjetAK8PHF, &b_pfjetAK8PHF);
  fChain->SetBranchAddress("pfjetAK8EEF", pfjetAK8EEF, &b_pfjetAK8EEF);
  fChain->SetBranchAddress("pfjetAK8HFHF", pfjetAK8HFHF, &b_pfjetAK8HFHF);
  fChain->SetBranchAddress("pfjetAK8CHM", pfjetAK8CHM, &b_pfjetAK8CHM);
  fChain->SetBranchAddress("pfjetAK8NHM", pfjetAK8NHM, &b_pfjetAK8NHM);
  fChain->SetBranchAddress("pfjetAK8MUM", pfjetAK8MUM, &b_pfjetAK8MUM);
  fChain->SetBranchAddress("pfjetAK8PHM", pfjetAK8PHM, &b_pfjetAK8PHM);
  fChain->SetBranchAddress("pfjetAK8EEM", pfjetAK8EEM, &b_pfjetAK8EEM);
  fChain->SetBranchAddress("pfjetAK8HFHM", pfjetAK8HFHM, &b_pfjetAK8HFHM);
  fChain->SetBranchAddress("pfjetAK8Neucons", pfjetAK8Neucons, &b_pfjetAK8Neucons);
  fChain->SetBranchAddress("pfjetAK8Chcons", pfjetAK8Chcons, &b_pfjetAK8Chcons);

  fChain->SetBranchAddress("pfjetAK8JER", pfjetAK8reso, &b_pfjetAK8reso);
  fChain->SetBranchAddress("pfjetAK8JERup", pfjetAK8resoup, &b_pfjetAK8resoup);
  fChain->SetBranchAddress("pfjetAK8JERdn", pfjetAK8resodn, &b_pfjetAK8resodn);

  fChain->SetBranchAddress("pfjetAK8jesup_AbsoluteStat",pfjetAK8jesup_AbsoluteStat, &b_pfjetAK8jesup_AbsoluteStat);
  fChain->SetBranchAddress("pfjetAK8jesup_AbsoluteScale",pfjetAK8jesup_AbsoluteScale, &b_pfjetAK8jesup_AbsoluteScale);
  fChain->SetBranchAddress("pfjetAK8jesup_AbsoluteMPFBias",pfjetAK8jesup_AbsoluteMPFBias, &b_pfjetAK8jesup_AbsoluteMPFBias);
  fChain->SetBranchAddress("pfjetAK8jesup_FlavorQCD",pfjetAK8jesup_FlavorQCD, &b_pfjetAK8jesup_FlavorQCD);
  fChain->SetBranchAddress("pfjetAK8jesup_Fragmentation",pfjetAK8jesup_Fragmentation, &b_pfjetAK8jesup_Fragmentation);
  fChain->SetBranchAddress("pfjetAK8jesup_PileUpDataMC",pfjetAK8jesup_PileUpDataMC, &b_pfjetAK8jesup_PileUpDataMC);
  fChain->SetBranchAddress("pfjetAK8jesup_PileUpPtBB",pfjetAK8jesup_PileUpPtBB, &b_pfjetAK8jesup_PileUpPtBB);
  fChain->SetBranchAddress("pfjetAK8jesup_PileUpPtEC1",pfjetAK8jesup_PileUpPtEC1, &b_pfjetAK8jesup_PileUpPtEC1);
  fChain->SetBranchAddress("pfjetAK8jesup_PileUpPtEC2",pfjetAK8jesup_PileUpPtEC2, &b_pfjetAK8jesup_PileUpPtEC2);
  fChain->SetBranchAddress("pfjetAK8jesup_PileUpPtHF",pfjetAK8jesup_PileUpPtHF, &b_pfjetAK8jesup_PileUpPtHF);
  fChain->SetBranchAddress("pfjetAK8jesup_PileUpPtRef",pfjetAK8jesup_PileUpPtRef, &b_pfjetAK8jesup_PileUpPtRef);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeFSR",pfjetAK8jesup_RelativeFSR, &b_pfjetAK8jesup_RelativeFSR);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeJEREC1",pfjetAK8jesup_RelativeJEREC1, &b_pfjetAK8jesup_RelativeJEREC1);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeJEREC2",pfjetAK8jesup_RelativeJEREC2, &b_pfjetAK8jesup_RelativeJEREC2);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeJERHF",pfjetAK8jesup_RelativeJERHF, &b_pfjetAK8jesup_RelativeJERHF);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativePtBB",pfjetAK8jesup_RelativePtBB, &b_pfjetAK8jesup_RelativePtBB);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativePtEC1",pfjetAK8jesup_RelativePtEC1, &b_pfjetAK8jesup_RelativePtEC1);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativePtEC2",pfjetAK8jesup_RelativePtEC2, &b_pfjetAK8jesup_RelativePtEC2);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativePtHF",pfjetAK8jesup_RelativePtHF, &b_pfjetAK8jesup_RelativePtHF);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeBal",pfjetAK8jesup_RelativeBal, &b_pfjetAK8jesup_RelativeBal);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeSample",pfjetAK8jesup_RelativeSample, &b_pfjetAK8jesup_RelativeSample);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeStatEC",pfjetAK8jesup_RelativeStatEC, &b_pfjetAK8jesup_RelativeStatEC);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeStatHF",pfjetAK8jesup_RelativeStatHF, &b_pfjetAK8jesup_RelativeStatHF);
  fChain->SetBranchAddress("pfjetAK8jesup_RelativeStatFSR",pfjetAK8jesup_RelativeStatFSR, &b_pfjetAK8jesup_RelativeStatFSR);
  fChain->SetBranchAddress("pfjetAK8jesup_SinglePionECAL",pfjetAK8jesup_SinglePionECAL, &b_pfjetAK8jesup_SinglePionECAL);
  fChain->SetBranchAddress("pfjetAK8jesup_SinglePionHCAL",pfjetAK8jesup_SinglePionHCAL, &b_pfjetAK8jesup_SinglePionHCAL);
  fChain->SetBranchAddress("pfjetAK8jesup_TimePtEta",pfjetAK8jesup_TimePtEta, &b_pfjetAK8jesup_TimePtEta);
  fChain->SetBranchAddress("pfjetAK8jesup_Total",pfjetAK8jesup_Total, &b_pfjetAK8jesup_Total);

  fChain->SetBranchAddress("pfjetAK8jesdn_AbsoluteStat",pfjetAK8jesdn_AbsoluteStat, &b_pfjetAK8jesdn_AbsoluteStat);
  fChain->SetBranchAddress("pfjetAK8jesdn_AbsoluteScale",pfjetAK8jesdn_AbsoluteScale, &b_pfjetAK8jesdn_AbsoluteScale);
  fChain->SetBranchAddress("pfjetAK8jesdn_AbsoluteMPFBias",pfjetAK8jesdn_AbsoluteMPFBias, &b_pfjetAK8jesdn_AbsoluteMPFBias);
  fChain->SetBranchAddress("pfjetAK8jesdn_FlavorQCD",pfjetAK8jesdn_FlavorQCD, &b_pfjetAK8jesdn_FlavorQCD);
  fChain->SetBranchAddress("pfjetAK8jesdn_Fragmentation",pfjetAK8jesdn_Fragmentation, &b_pfjetAK8jesdn_Fragmentation);
  fChain->SetBranchAddress("pfjetAK8jesdn_PileUpDataMC",pfjetAK8jesdn_PileUpDataMC, &b_pfjetAK8jesdn_PileUpDataMC);
  fChain->SetBranchAddress("pfjetAK8jesdn_PileUpPtBB",pfjetAK8jesdn_PileUpPtBB, &b_pfjetAK8jesdn_PileUpPtBB);
  fChain->SetBranchAddress("pfjetAK8jesdn_PileUpPtEC1",pfjetAK8jesdn_PileUpPtEC1, &b_pfjetAK8jesdn_PileUpPtEC1);
  fChain->SetBranchAddress("pfjetAK8jesdn_PileUpPtEC2",pfjetAK8jesdn_PileUpPtEC2, &b_pfjetAK8jesdn_PileUpPtEC2);
  fChain->SetBranchAddress("pfjetAK8jesdn_PileUpPtHF",pfjetAK8jesdn_PileUpPtHF, &b_pfjetAK8jesdn_PileUpPtHF);
  fChain->SetBranchAddress("pfjetAK8jesdn_PileUpPtRef",pfjetAK8jesdn_PileUpPtRef, &b_pfjetAK8jesdn_PileUpPtRef);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeFSR",pfjetAK8jesdn_RelativeFSR, &b_pfjetAK8jesdn_RelativeFSR);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeJEREC1",pfjetAK8jesdn_RelativeJEREC1, &b_pfjetAK8jesdn_RelativeJEREC1);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeJEREC2",pfjetAK8jesdn_RelativeJEREC2, &b_pfjetAK8jesdn_RelativeJEREC2);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeJERHF",pfjetAK8jesdn_RelativeJERHF, &b_pfjetAK8jesdn_RelativeJERHF);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativePtBB",pfjetAK8jesdn_RelativePtBB, &b_pfjetAK8jesdn_RelativePtBB);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativePtEC1",pfjetAK8jesdn_RelativePtEC1, &b_pfjetAK8jesdn_RelativePtEC1);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativePtEC2",pfjetAK8jesdn_RelativePtEC2, &b_pfjetAK8jesdn_RelativePtEC2);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativePtHF",pfjetAK8jesdn_RelativePtHF, &b_pfjetAK8jesdn_RelativePtHF);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeBal",pfjetAK8jesdn_RelativeBal, &b_pfjetAK8jesdn_RelativeBal);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeSample",pfjetAK8jesdn_RelativeSample, &b_pfjetAK8jesdn_RelativeSample);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeStatEC",pfjetAK8jesdn_RelativeStatEC, &b_pfjetAK8jesdn_RelativeStatEC);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeStatHF",pfjetAK8jesdn_RelativeStatHF, &b_pfjetAK8jesdn_RelativeStatHF);
  fChain->SetBranchAddress("pfjetAK8jesdn_RelativeStatFSR",pfjetAK8jesdn_RelativeStatFSR, &b_pfjetAK8jesdn_RelativeStatFSR);
  fChain->SetBranchAddress("pfjetAK8jesdn_SinglePionECAL",pfjetAK8jesdn_SinglePionECAL, &b_pfjetAK8jesdn_SinglePionECAL);
  fChain->SetBranchAddress("pfjetAK8jesdn_SinglePionHCAL",pfjetAK8jesdn_SinglePionHCAL, &b_pfjetAK8jesdn_SinglePionHCAL);
  fChain->SetBranchAddress("pfjetAK8jesdn_TimePtEta",pfjetAK8jesdn_TimePtEta, &b_pfjetAK8jesdn_TimePtEta);
  fChain->SetBranchAddress("pfjetAK8jesdn_Total",pfjetAK8jesdn_Total, &b_pfjetAK8jesdn_Total);

  fChain->SetBranchAddress("pfjetAK8chrad", pfjetAK8chrad, &b_pfjetAK8chrad);
  fChain->SetBranchAddress("pfjetAK8pTD", pfjetAK8pTD, &b_pfjetAK8pTD);
  fChain->SetBranchAddress("pfjetAK8sdmass", pfjetAK8sdmass, &b_pfjetAK8sdmass);
  fChain->SetBranchAddress("pfjetAK8tau1", pfjetAK8tau1, &b_pfjetAK8tau1);
  fChain->SetBranchAddress("pfjetAK8tau2", pfjetAK8tau2, &b_pfjetAK8tau2);
  fChain->SetBranchAddress("pfjetAK8tau3", pfjetAK8tau3, &b_pfjetAK8tau3);
  fChain->SetBranchAddress("pfjetAK8sub1pt", pfjetAK8sub1pt, &b_pfjetAK8sub1pt);
  fChain->SetBranchAddress("pfjetAK8sub1eta", pfjetAK8sub1eta, &b_pfjetAK8sub1eta);
  fChain->SetBranchAddress("pfjetAK8sub1phi", pfjetAK8sub1phi, &b_pfjetAK8sub1phi);
  fChain->SetBranchAddress("pfjetAK8sub1mass", pfjetAK8sub1mass, &b_pfjetAK8sub1mass);
  fChain->SetBranchAddress("pfjetAK8sub1btag", pfjetAK8sub1btag, &b_pfjetAK8sub1btag);
  fChain->SetBranchAddress("pfjetAK8sub1chhadfrac", pfjetAK8sub1chhadfrac, &b_pfjetAK8sub1chhadfrac);
  fChain->SetBranchAddress("pfjetAK8sub1neuhadfrac", pfjetAK8sub1neuhadfrac, &b_pfjetAK8sub1neuhadfrac);
  fChain->SetBranchAddress("pfjetAK8sub1emfrac", pfjetAK8sub1emfrac, &b_pfjetAK8sub1emfrac);
  fChain->SetBranchAddress("pfjetAK8sub1phofrac", pfjetAK8sub1phofrac, &b_pfjetAK8sub1phofrac);
  fChain->SetBranchAddress("pfjetAK8sub1mufrac", pfjetAK8sub1mufrac, &b_pfjetAK8sub1mufrac);
  fChain->SetBranchAddress("pfjetAK8sub2pt", pfjetAK8sub2pt, &b_pfjetAK8sub2pt);
  fChain->SetBranchAddress("pfjetAK8sub2eta", pfjetAK8sub2eta, &b_pfjetAK8sub2eta);
  fChain->SetBranchAddress("pfjetAK8sub2phi", pfjetAK8sub2phi, &b_pfjetAK8sub2phi);
  fChain->SetBranchAddress("pfjetAK8sub2mass", pfjetAK8sub2mass, &b_pfjetAK8sub2mass);
  fChain->SetBranchAddress("pfjetAK8sub2btag", pfjetAK8sub2btag, &b_pfjetAK8sub2btag);
  fChain->SetBranchAddress("pfjetAK8sub2chhadfrac", pfjetAK8sub2chhadfrac, &b_pfjetAK8sub2chhadfrac);
  fChain->SetBranchAddress("pfjetAK8sub2neuhadfrac", pfjetAK8sub2neuhadfrac, &b_pfjetAK8sub2neuhadfrac);
  fChain->SetBranchAddress("pfjetAK8sub2emfrac", pfjetAK8sub2emfrac, &b_pfjetAK8sub2emfrac);
  fChain->SetBranchAddress("pfjetAK8sub2phofrac", pfjetAK8sub2phofrac, &b_pfjetAK8sub2phofrac);
  fChain->SetBranchAddress("pfjetAK8sub2mufrac", pfjetAK8sub2mufrac, &b_pfjetAK8sub2mufrac);
   
  fChain->SetBranchAddress("pfjetAK8subhaddiff", pfjetAK8subhaddiff, &b_pfjetAK8subhaddiff);
  fChain->SetBranchAddress("pfjetAK8subemdiff", pfjetAK8subemdiff, &b_pfjetAK8subemdiff);
  fChain->SetBranchAddress("pfjetAK8subptdiff", pfjetAK8subptdiff, &b_pfjetAK8subptdiff);

  fChain->SetBranchAddress("pfjetAK8elinsubpt", pfjetAK8elinsubpt, &b_pfjetAK8elinsubpt);
  fChain->SetBranchAddress("pfjetAK8elinsubeta", pfjetAK8elinsubeta, &b_pfjetAK8elinsubeta);
  fChain->SetBranchAddress("pfjetAK8elinsubphi", pfjetAK8elinsubphi, &b_pfjetAK8elinsubphi);
  fChain->SetBranchAddress("pfjetAK8elinsubmass", pfjetAK8elinsubmass, &b_pfjetAK8elinsubmass);

  fChain->SetBranchAddress("pfjetAK8elinsubjpt", pfjetAK8elinsubjpt, &b_pfjetAK8elinsubjpt);
  fChain->SetBranchAddress("pfjetAK8elinsubjeta", pfjetAK8elinsubjeta, &b_pfjetAK8elinsubjeta);
  fChain->SetBranchAddress("pfjetAK8elinsubjphi", pfjetAK8elinsubjphi, &b_pfjetAK8elinsubjphi);
  fChain->SetBranchAddress("pfjetAK8elinsubjmass", pfjetAK8elinsubjmass, &b_pfjetAK8elinsubjmass);

  fChain->SetBranchAddress("pfjetAK8muinsubjpt", pfjetAK8muinsubjpt, &b_pfjetAK8muinsubjpt);
  fChain->SetBranchAddress("pfjetAK8muinsubjeta", pfjetAK8muinsubjeta, &b_pfjetAK8muinsubjeta);
  fChain->SetBranchAddress("pfjetAK8muinsubjphi", pfjetAK8muinsubjphi, &b_pfjetAK8muinsubjphi);
  fChain->SetBranchAddress("pfjetAK8muinsubjmass", pfjetAK8muinsubjmass, &b_pfjetAK8muinsubjmass);


  fChain->SetBranchAddress("pfjetAK8muinsubpt", pfjetAK8muinsubpt, &b_pfjetAK8muinsubpt);
  fChain->SetBranchAddress("pfjetAK8muinsubeta", pfjetAK8muinsubeta, &b_pfjetAK8muinsubeta);
  fChain->SetBranchAddress("pfjetAK8muinsubphi", pfjetAK8muinsubphi, &b_pfjetAK8muinsubphi);
  fChain->SetBranchAddress("pfjetAK8muinsubmass", pfjetAK8muinsubmass, &b_pfjetAK8muinsubmass);
  fChain->SetBranchAddress("pfjetAK8muinsubIfar", pfjetAK8muinsubIfar, &b_pfjetAK8muinsubIfar);
  fChain->SetBranchAddress("pfjetAK8muinsubInear", pfjetAK8muinsubInear, &b_pfjetAK8muinsubInear);
  fChain->SetBranchAddress("pfjetAK8muinsubI0", pfjetAK8muinsubI0, &b_pfjetAK8muinsubI0);

  fChain->SetBranchAddress("npfjetAK4", &npfjetAK4, &b_npfjetAK4);
  fChain->SetBranchAddress("pfjetAK4jetID", pfjetAK4jetID, &b_pfjetAK4jetID);
  fChain->SetBranchAddress("pfjetAK4jetID_tightlepveto", pfjetAK4jetID_tightlepveto, &b_pfjetAK4jetID_tightlepveto);
  fChain->SetBranchAddress("pfjetAK4pt", pfjetAK4pt, &b_pfjetAK4pt);
  fChain->SetBranchAddress("pfjetAK4eta", pfjetAK4eta, &b_pfjetAK4eta);
  fChain->SetBranchAddress("pfjetAK4y", pfjetAK4y, &b_pfjetAK4y);
  fChain->SetBranchAddress("pfjetAK4phi", pfjetAK4phi, &b_pfjetAK4phi);
  fChain->SetBranchAddress("pfjetAK4mass", pfjetAK4mass, &b_pfjetAK4mass);
  fChain->SetBranchAddress("pfjetAK4JEC", pfjetAK4JEC, &b_pfjetAK4JEC);
   
  fChain->SetBranchAddress("pfjetAK4btag_DeepCSV", pfjetAK4btag_DeepCSV, &b_pfjetAK4btag_DeepCSV);
  fChain->SetBranchAddress("pfjetAK4btag_DeepFlav", pfjetAK4btag_DeepFlav, &b_pfjetAK4btag_DeepFlav);
   
  fChain->SetBranchAddress("pfjetAK4JER", pfjetAK4reso, &b_pfjetAK4reso);
  fChain->SetBranchAddress("pfjetAK4JERup", pfjetAK4resoup, &b_pfjetAK4resoup);
  fChain->SetBranchAddress("pfjetAK4JERdn", pfjetAK4resodn, &b_pfjetAK4resodn);

  fChain->SetBranchAddress("pfjetAK4jesup_AbsoluteStat",pfjetAK4jesup_AbsoluteStat, &b_pfjetAK4jesup_AbsoluteStat);
  fChain->SetBranchAddress("pfjetAK4jesup_AbsoluteScale",pfjetAK4jesup_AbsoluteScale, &b_pfjetAK4jesup_AbsoluteScale);
  fChain->SetBranchAddress("pfjetAK4jesup_AbsoluteMPFBias",pfjetAK4jesup_AbsoluteMPFBias, &b_pfjetAK4jesup_AbsoluteMPFBias);
  fChain->SetBranchAddress("pfjetAK4jesup_FlavorQCD",pfjetAK4jesup_FlavorQCD, &b_pfjetAK4jesup_FlavorQCD);
  fChain->SetBranchAddress("pfjetAK4jesup_Fragmentation",pfjetAK4jesup_Fragmentation, &b_pfjetAK4jesup_Fragmentation);
  fChain->SetBranchAddress("pfjetAK4jesup_PileUpDataMC",pfjetAK4jesup_PileUpDataMC, &b_pfjetAK4jesup_PileUpDataMC);
  fChain->SetBranchAddress("pfjetAK4jesup_PileUpPtBB",pfjetAK4jesup_PileUpPtBB, &b_pfjetAK4jesup_PileUpPtBB);
  fChain->SetBranchAddress("pfjetAK4jesup_PileUpPtEC1",pfjetAK4jesup_PileUpPtEC1, &b_pfjetAK4jesup_PileUpPtEC1);
  fChain->SetBranchAddress("pfjetAK4jesup_PileUpPtEC2",pfjetAK4jesup_PileUpPtEC2, &b_pfjetAK4jesup_PileUpPtEC2);
  fChain->SetBranchAddress("pfjetAK4jesup_PileUpPtHF",pfjetAK4jesup_PileUpPtHF, &b_pfjetAK4jesup_PileUpPtHF);
  fChain->SetBranchAddress("pfjetAK4jesup_PileUpPtRef",pfjetAK4jesup_PileUpPtRef, &b_pfjetAK4jesup_PileUpPtRef);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeFSR",pfjetAK4jesup_RelativeFSR, &b_pfjetAK4jesup_RelativeFSR);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeJEREC1",pfjetAK4jesup_RelativeJEREC1, &b_pfjetAK4jesup_RelativeJEREC1);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeJEREC2",pfjetAK4jesup_RelativeJEREC2, &b_pfjetAK4jesup_RelativeJEREC2);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeJERHF",pfjetAK4jesup_RelativeJERHF, &b_pfjetAK4jesup_RelativeJERHF);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativePtBB",pfjetAK4jesup_RelativePtBB, &b_pfjetAK4jesup_RelativePtBB);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativePtEC1",pfjetAK4jesup_RelativePtEC1, &b_pfjetAK4jesup_RelativePtEC1);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativePtEC2",pfjetAK4jesup_RelativePtEC2, &b_pfjetAK4jesup_RelativePtEC2);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativePtHF",pfjetAK4jesup_RelativePtHF, &b_pfjetAK4jesup_RelativePtHF);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeBal",pfjetAK4jesup_RelativeBal, &b_pfjetAK4jesup_RelativeBal);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeSample",pfjetAK4jesup_RelativeSample, &b_pfjetAK4jesup_RelativeSample);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeStatEC",pfjetAK4jesup_RelativeStatEC, &b_pfjetAK4jesup_RelativeStatEC);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeStatHF",pfjetAK4jesup_RelativeStatHF, &b_pfjetAK4jesup_RelativeStatHF);
  fChain->SetBranchAddress("pfjetAK4jesup_RelativeStatFSR",pfjetAK4jesup_RelativeStatFSR, &b_pfjetAK4jesup_RelativeStatFSR);
  fChain->SetBranchAddress("pfjetAK4jesup_SinglePionECAL",pfjetAK4jesup_SinglePionECAL, &b_pfjetAK4jesup_SinglePionECAL);
  fChain->SetBranchAddress("pfjetAK4jesup_SinglePionHCAL",pfjetAK4jesup_SinglePionHCAL, &b_pfjetAK4jesup_SinglePionHCAL);
  fChain->SetBranchAddress("pfjetAK4jesup_TimePtEta",pfjetAK4jesup_TimePtEta, &b_pfjetAK4jesup_TimePtEta);
  fChain->SetBranchAddress("pfjetAK4jesup_Total",pfjetAK4jesup_Total, &b_pfjetAK4jesup_Total);

  fChain->SetBranchAddress("pfjetAK4jesdn_AbsoluteStat",pfjetAK4jesdn_AbsoluteStat, &b_pfjetAK4jesdn_AbsoluteStat);
  fChain->SetBranchAddress("pfjetAK4jesdn_AbsoluteScale",pfjetAK4jesdn_AbsoluteScale, &b_pfjetAK4jesdn_AbsoluteScale);
  fChain->SetBranchAddress("pfjetAK4jesdn_AbsoluteMPFBias",pfjetAK4jesdn_AbsoluteMPFBias, &b_pfjetAK4jesdn_AbsoluteMPFBias);
  fChain->SetBranchAddress("pfjetAK4jesdn_FlavorQCD",pfjetAK4jesdn_FlavorQCD, &b_pfjetAK4jesdn_FlavorQCD);
  fChain->SetBranchAddress("pfjetAK4jesdn_Fragmentation",pfjetAK4jesdn_Fragmentation, &b_pfjetAK4jesdn_Fragmentation);
  fChain->SetBranchAddress("pfjetAK4jesdn_PileUpDataMC",pfjetAK4jesdn_PileUpDataMC, &b_pfjetAK4jesdn_PileUpDataMC);
  fChain->SetBranchAddress("pfjetAK4jesdn_PileUpPtBB",pfjetAK4jesdn_PileUpPtBB, &b_pfjetAK4jesdn_PileUpPtBB);
  fChain->SetBranchAddress("pfjetAK4jesdn_PileUpPtEC1",pfjetAK4jesdn_PileUpPtEC1, &b_pfjetAK4jesdn_PileUpPtEC1);
  fChain->SetBranchAddress("pfjetAK4jesdn_PileUpPtEC2",pfjetAK4jesdn_PileUpPtEC2, &b_pfjetAK4jesdn_PileUpPtEC2);
  fChain->SetBranchAddress("pfjetAK4jesdn_PileUpPtHF",pfjetAK4jesdn_PileUpPtHF, &b_pfjetAK4jesdn_PileUpPtHF);
  fChain->SetBranchAddress("pfjetAK4jesdn_PileUpPtRef",pfjetAK4jesdn_PileUpPtRef, &b_pfjetAK4jesdn_PileUpPtRef);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeFSR",pfjetAK4jesdn_RelativeFSR, &b_pfjetAK4jesdn_RelativeFSR);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeJEREC1",pfjetAK4jesdn_RelativeJEREC1, &b_pfjetAK4jesdn_RelativeJEREC1);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeJEREC2",pfjetAK4jesdn_RelativeJEREC2, &b_pfjetAK4jesdn_RelativeJEREC2);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeJERHF",pfjetAK4jesdn_RelativeJERHF, &b_pfjetAK4jesdn_RelativeJERHF);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativePtBB",pfjetAK4jesdn_RelativePtBB, &b_pfjetAK4jesdn_RelativePtBB);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativePtEC1",pfjetAK4jesdn_RelativePtEC1, &b_pfjetAK4jesdn_RelativePtEC1);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativePtEC2",pfjetAK4jesdn_RelativePtEC2, &b_pfjetAK4jesdn_RelativePtEC2);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativePtHF",pfjetAK4jesdn_RelativePtHF, &b_pfjetAK4jesdn_RelativePtHF);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeBal",pfjetAK4jesdn_RelativeBal, &b_pfjetAK4jesdn_RelativeBal);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeSample",pfjetAK4jesdn_RelativeSample, &b_pfjetAK4jesdn_RelativeSample);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeStatEC",pfjetAK4jesdn_RelativeStatEC, &b_pfjetAK4jesdn_RelativeStatEC);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeStatHF",pfjetAK4jesdn_RelativeStatHF, &b_pfjetAK4jesdn_RelativeStatHF);
  fChain->SetBranchAddress("pfjetAK4jesdn_RelativeStatFSR",pfjetAK4jesdn_RelativeStatFSR, &b_pfjetAK4jesdn_RelativeStatFSR);
  fChain->SetBranchAddress("pfjetAK4jesdn_SinglePionECAL",pfjetAK4jesdn_SinglePionECAL, &b_pfjetAK4jesdn_SinglePionECAL);
  fChain->SetBranchAddress("pfjetAK4jesdn_SinglePionHCAL",pfjetAK4jesdn_SinglePionHCAL, &b_pfjetAK4jesdn_SinglePionHCAL);
  fChain->SetBranchAddress("pfjetAK4jesdn_TimePtEta",pfjetAK4jesdn_TimePtEta, &b_pfjetAK4jesdn_TimePtEta);
  fChain->SetBranchAddress("pfjetAK4jesdn_Total",pfjetAK4jesdn_Total, &b_pfjetAK4jesdn_Total);

  fChain->SetBranchAddress("pfjetAK4hadronflav", pfjetAK4hadronflav, &b_pfjetAK4hadronflav);
  fChain->SetBranchAddress("pfjetAK4partonflav", pfjetAK4partonflav, &b_pfjetAK4partonflav);
  fChain->SetBranchAddress("pfjetAK4qgl", pfjetAK4qgl, &b_pfjetAK4qgl);
  fChain->SetBranchAddress("pfjetAK4PUID", pfjetAK4PUID, &b_pfjetAK4PUID);
  fChain->SetBranchAddress("pfjetAK4GenMatch", &pfjetAK4GenMatch, &b_pfjetAK4GenMatch);

  fChain->SetBranchAddress("GENMET", &GENMET, &b_genmiset);
  fChain->SetBranchAddress("GENMETPhi", &GENMETPhi, &b_genmisphi);
   
  fChain->SetBranchAddress("ngenjetAK8", &ngenjetAK8, &b_ngenjetAK8);
  fChain->SetBranchAddress("genjetAK8pt", genjetAK8pt, &b_genjetAK8pt);
  fChain->SetBranchAddress("genjetAK8eta", genjetAK8eta, &b_genjetAK8eta);
  fChain->SetBranchAddress("genjetAK8phi", genjetAK8phi, &b_genjetAK8phi);
  fChain->SetBranchAddress("genjetAK8mass", genjetAK8mass, &b_genjetAK8mass);
  fChain->SetBranchAddress("genjetAK8sdmass", genjetAK8sdmass, &b_genjetAK8sdmass);
   
  fChain->SetBranchAddress("ngenjetAK4", &ngenjetAK4, &b_ngenjetAK4);
  fChain->SetBranchAddress("genjetAK4pt", genjetAK4pt, &b_genjetAK4pt);
  fChain->SetBranchAddress("genjetAK4phi", genjetAK4phi, &b_genjetAK4phi);
  fChain->SetBranchAddress("genjetAK4mass", genjetAK4mass, &b_genjetAK4mass);
  fChain->SetBranchAddress("genjetAK4eta", genjetAK4eta, &b_genjetAK4eta);
  fChain->SetBranchAddress("genjetAK4hadronflav",genjetAK4hadronflav,&b_genjetAK4hadronflav);
  fChain->SetBranchAddress("genjetAK4partonflav",genjetAK4partonflav,&b_genjetAK4partonflav);
  fChain->SetBranchAddress("ngenparticles", &ngenparticles, &b_ngenparticles);
  fChain->SetBranchAddress("genpartstatus", genpartstatus, &b_genpartstatus);
  fChain->SetBranchAddress("genpartpdg", genpartpdg, &b_genpartpdg);
  fChain->SetBranchAddress("genpartmompdg", genpartmompdg, &b_genpartmompdg);
  fChain->SetBranchAddress("genpartmomstatus", genpartmomstatus, &b_genpartmomstatus);
  fChain->SetBranchAddress("genpartgrmompdg", genpartgrmompdg, &b_genpartgrmompdg);
  fChain->SetBranchAddress("genpartdaugno", genpartdaugno, &b_genpartdaugno);
  fChain->SetBranchAddress("genpartfromhard", genpartfromhard, &b_genpartfromhard);
  fChain->SetBranchAddress("genpartfromhardbFSR", genpartfromhardbFSR, &b_genpartfromhardbFSR);
  fChain->SetBranchAddress("genpartisPromptFinalState", genpartisPromptFinalState, &b_genpartisPromptFinalState); //Need to ask Suman as added in new setup
  fChain->SetBranchAddress("genpartisLastCopyBeforeFSR", genpartisLastCopyBeforeFSR, &b_genpartisLastCopyBeforeFSR);
  fChain->SetBranchAddress("genpartpt", genpartpt, &b_genpartpt);
  fChain->SetBranchAddress("genparteta", genparteta, &b_genparteta);
  fChain->SetBranchAddress("genpartphi", genpartphi, &b_genpartphi);
  fChain->SetBranchAddress("genpartm", genpartm, &b_genpartm);

  fChain->SetBranchAddress("Generator_x1",&Generator_x1, &b_Generator_x1);
  fChain->SetBranchAddress("Generator_x2",&Generator_x2, &b_Generator_x2);
  fChain->SetBranchAddress("Generator_xpdf1",&Generator_xpdf1, &b_Generator_xpdf1);
  fChain->SetBranchAddress("Generator_xpdf2",&Generator_xpdf2, &b_Generator_xpdf2);
  fChain->SetBranchAddress("Generator_id1",&Generator_id1, &b_Generator_id1);
  fChain->SetBranchAddress("Generator_id2",&Generator_id2, &b_Generator_id2);
  fChain->SetBranchAddress("Generator_scalePDF",&Generator_scalePDF, &b_Generator_scalePDF);

  fChain->SetBranchAddress("LHE_weight",&LHE_weight, &b_LHE_weight);
  fChain->SetBranchAddress("nLHEScaleWeights",&nLHEScaleWeights, &b_nLHEScaleWeights);
  fChain->SetBranchAddress("LHEScaleWeights",LHEScaleWeights, &b_LHEScaleWeights);
  fChain->SetBranchAddress("nLHEPDFWeights",&nLHEPDFWeights, &b_nLHEPDFWeights);
  //fChain->SetBranchAddress("LHEPDFWeights",LHEPDFWeights, &b_LHEPDFWeights);
  fChain->SetBranchAddress("nLHEAlpsWeights",&nLHEAlpsWeights, &b_nLHEAlpsWeights);
  fChain->SetBranchAddress("LHEAlpsWeights",LHEAlpsWeights, &b_LHEAlpsWeights);
  fChain->SetBranchAddress("nLHEPSWeights",&nLHEPSWeights, &b_nLHEPSWeights);
  fChain->SetBranchAddress("LHEPSWeights",LHEPSWeights, &b_LHEPSWeights);
  fChain->SetBranchAddress("nLHEparticles", &nLHEparticles, &b_nLHEparticles);
  fChain->SetBranchAddress("LHEpartpdg", LHEpartpdg, &b_LHEpartpdg);
  fChain->SetBranchAddress("LHEpartpt", LHEpartpt, &b_LHEpartpt);
  fChain->SetBranchAddress("LHEparteta", LHEparteta, &b_LHEparteta);
  fChain->SetBranchAddress("LHEpartphi", LHEpartphi, &b_LHEpartphi);
  fChain->SetBranchAddress("LHEpartm", LHEpartm, &b_LHEpartm);

  fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
  fChain->SetBranchAddress("muonisPF", muonisPF, &b_muonisPF);
  fChain->SetBranchAddress("muonisGL", muonisGL, &b_muonisGL);
  fChain->SetBranchAddress("muonisTRK", muonisTRK, &b_muonisTRK);
  fChain->SetBranchAddress("muonisLoose", muonisLoose, &b_muonisLoose);
  fChain->SetBranchAddress("muonisGoodGL", muonisGoodGL, &b_muonisGoodGL);
  fChain->SetBranchAddress("muonisMed", muonisMed, &b_muonisMed);
  fChain->SetBranchAddress("muonisTight", muonisTight, &b_muonisTight);
  fChain->SetBranchAddress("muonisHighPt", muonisHighPt, &b_muonisHighPt);
  fChain->SetBranchAddress("muonisHighPttrk", muonisHighPttrk, &b_muonisHighPttrk);
  fChain->SetBranchAddress("muonisMedPr", muonisMedPr, &b_muonisMedPr);
  fChain->SetBranchAddress("muonpt", muonpt, &b_muonpt);
  fChain->SetBranchAddress("muonp", muonp, &b_muonp);
  fChain->SetBranchAddress("muoneta", muoneta, &b_muoneta);
  fChain->SetBranchAddress("muonphi", muonphi, &b_muonphi);
  fChain->SetBranchAddress("muoncharge", muoncharge, &b_muoncharge);
  fChain->SetBranchAddress("muontrkvtx", muontrkvtx, &b_muontrkvtx);
  fChain->SetBranchAddress("muondz", muondz, &b_muondz);
  fChain->SetBranchAddress("muonpter", muonpter, &b_muonpter);
  fChain->SetBranchAddress("muonchi", muonchi, &b_muonchi);
  fChain->SetBranchAddress("muonndf", muonndf, &b_muonndf);
  fChain->SetBranchAddress("muonecal", muonecal, &b_muonecal);
  fChain->SetBranchAddress("muonhcal", muonhcal, &b_muonhcal);
  fChain->SetBranchAddress("muonpfiso", muonpfiso, &b_muonpfiso);
  fChain->SetBranchAddress("muonposmatch", muonposmatch, &b_muonposmatch);
  fChain->SetBranchAddress("muontrkink", muontrkink, &b_muontrkink);
  fChain->SetBranchAddress("muonsegcom", muonsegcom, &b_muonsegcom);
  fChain->SetBranchAddress("muonpixhit", muonpixhit, &b_muonpixhit);
  fChain->SetBranchAddress("muonmst", muonmst, &b_muonmst);
  fChain->SetBranchAddress("muontrklay", muontrklay, &b_muontrklay);
  fChain->SetBranchAddress("muonvalfrac", muonvalfrac, &b_muonvalfrac);
  fChain->SetBranchAddress("muonchiso", muonchiso, &b_muonchiso);
  fChain->SetBranchAddress("muonnhiso", muonnhiso, &b_muonnhiso);
  fChain->SetBranchAddress("muonphiso", muonphiso, &b_muonphiso);
  fChain->SetBranchAddress("muonminisoall", muonminisoall, &b_muonminisoall);
  fChain->SetBranchAddress("mudxy_sv", mudxy_sv, &b_mudxy_sv);
  fChain->SetBranchAddress("muonhit", muonhit, &b_muonhit);
  fChain->SetBranchAddress("Muon_corrected_pt",Muon_corrected_pt, &b_Muon_corrected_pt);
  fChain->SetBranchAddress("Muon_correctedUp_pt",Muon_correctedUp_pt, &b_Muon_correctedUp_pt);
  fChain->SetBranchAddress("Muon_correctedDown_pt",Muon_correctedDown_pt, &b_Muon_correctedDown_pt);

  fChain->SetBranchAddress("nelecs", &nelecs, &b_nelecs);
  fChain->SetBranchAddress("elpt", elpt, &b_elpt);
  fChain->SetBranchAddress("eldxytrk", eldxytrk, &b_eldxytrk);
  fChain->SetBranchAddress("eldztrk", eldztrk, &b_eldztrk);
  fChain->SetBranchAddress("eldxy_sv", eldxy_sv, &b_eldxy_sv);
  fChain->SetBranchAddress("eleta", eleta, &b_eleta);
  fChain->SetBranchAddress("elphi", elphi, &b_elphi);
  fChain->SetBranchAddress("elp", elp, &b_elp);
  fChain->SetBranchAddress("ele", ele, &b_ele);
  fChain->SetBranchAddress("elcharge", elcharge, &b_elcharge);
  fChain->SetBranchAddress("elmvaid", elmvaid, &b_elmvaid);
  fChain->SetBranchAddress("elmvaid_Fallv2WP80", elmvaid_Fallv2WP80, &b_elmvaid_Fallv2WP80);
  fChain->SetBranchAddress("elmvaid_noIso", elmvaid_noIso, &b_elmvaid_noIso);
  fChain->SetBranchAddress("elmvaid_Fallv2WP80_noIso", elmvaid_Fallv2WP80_noIso, &b_elmvaid_Fallv2WP80_noIso);
  fChain->SetBranchAddress("elcutid_veto", elcutid_veto, &b_elcutid_veto);
  fChain->SetBranchAddress("elcutid_loose", elcutid_loose, &b_elcutid_loose);
  fChain->SetBranchAddress("elcutid_med", elcutid_med, &b_elcutid_med);
  fChain->SetBranchAddress("elcutid_tight", elcutid_tight, &b_elcutid_tight);
  fChain->SetBranchAddress("elmvaid_wpLoose", elmvaid_wpLoose, &b_elmvaid_wpLoose);
  fChain->SetBranchAddress("elmvaid_noIso_wpLoose", elmvaid_noIso_wpLoose, &b_elmvaid_noIso_wpLoose);
  fChain->SetBranchAddress("elsupcl_e", elsupcl_e, &b_elsupcl_e);
  fChain->SetBranchAddress("elhovere", elhovere, &b_elhovere);
  fChain->SetBranchAddress("elchi", elchi, &b_elchi);
  fChain->SetBranchAddress("elndf", elndf, &b_elndf);
   
  fChain->SetBranchAddress("eletain", eletain, &b_eletain);
  fChain->SetBranchAddress("elphiin", elphiin, &b_elphiin);
   
  fChain->SetBranchAddress("elsupcl_eta", elsupcl_eta, &b_elsupcl_eta);
  fChain->SetBranchAddress("elsupcl_phi", elsupcl_phi, &b_elsupcl_phi);
  fChain->SetBranchAddress("elsupcl_rawE", elsupcl_rawE, &b_elsupcl_rawE);
  fChain->SetBranchAddress("elfbrem", elfbrem, &b_elfbrem);
   
  fChain->SetBranchAddress("eleoverp", eleoverp, &b_eleoverp);
  fChain->SetBranchAddress("elietaieta", elietaieta, &b_elietaieta);
  fChain->SetBranchAddress("elmisshits", elmisshits, &b_elmisshits);
   
  fChain->SetBranchAddress("elpfiso", elpfiso, &b_elpfiso);
   
  fChain->SetBranchAddress("elsigmaieta", elsigmaieta, &b_elsigmaieta);
  fChain->SetBranchAddress("elsigmaiphi", elsigmaiphi, &b_elsigmaiphi);
  fChain->SetBranchAddress("elr9full", elr9full, &b_elr9full);
  fChain->SetBranchAddress("elsupcl_etaw", elsupcl_etaw, &b_elsupcl_etaw); 
  fChain->SetBranchAddress("elsupcl_phiw", elsupcl_phiw, &b_elsupcl_phiw);
  fChain->SetBranchAddress("elhcaloverecal", elhcaloverecal, &b_elhcaloverecal);
  fChain->SetBranchAddress("elcloctftrkn",  elcloctftrkn, &b_elcloctftrkn);
  fChain->SetBranchAddress("elcloctftrkchi2", elcloctftrkchi2, &b_elcloctftrkchi2);
  fChain->SetBranchAddress("ele1x5bye5x5", ele1x5bye5x5, &b_ele1x5bye5x5);
  fChain->SetBranchAddress("elnormchi2", elnormchi2, &b_elnormchi2);
  fChain->SetBranchAddress("elhitsmiss", elhitsmiss, &b_elhitsmiss);
  fChain->SetBranchAddress("eltrkmeasure", eltrkmeasure, &b_eltrkmeasure);
  fChain->SetBranchAddress("elconvtxprob", elconvtxprob, &b_elconvtxprob);
  fChain->SetBranchAddress("elecloverpout", elecloverpout, &b_elecloverpout);
  fChain->SetBranchAddress("elecaletrkmomentum", elecaletrkmomentum, &b_elecaletrkmomentum);
  fChain->SetBranchAddress("eldeltaetacltrkcalo", eldeltaetacltrkcalo, &b_eldeltaetacltrkcalo);
  fChain->SetBranchAddress("elsupcl_preshvsrawe", elsupcl_preshvsrawe, &b_elsupcl_preshvsrawe);
  fChain->SetBranchAddress("elpfisolsumphet", elpfisolsumphet, &b_elpfisolsumphet);
  fChain->SetBranchAddress("elpfisolsumchhadpt", elpfisolsumchhadpt, &b_elpfisolsumchhadpt);
  fChain->SetBranchAddress("elpfsiolsumneuhadet", elpfsiolsumneuhadet, &b_elpfsiolsumneuhadet);

  fChain->SetBranchAddress("Electron_eccalTrkEnergyPostCorr",Electron_eccalTrkEnergyPostCorr, &b_Electron_eccalTrkEnergyPostCorr);
  fChain->SetBranchAddress("Electron_energyScaleValue",Electron_energyScaleValue, &b_Electron_energyScaleValue);
  fChain->SetBranchAddress("Electron_energyScaleUp",Electron_energyScaleUp, &b_Electron_energyScaleUp);
  fChain->SetBranchAddress("Electron_energyScaleDown",Electron_energyScaleDown, &b_Electron_energyScaleDown);
  fChain->SetBranchAddress("Electron_energySigmaValue",Electron_energySigmaValue, &b_Electron_energySigmaValue);
  fChain->SetBranchAddress("Electron_energySigmaUp",Electron_energySigmaUp, &b_Electron_energySigmaUp);
  fChain->SetBranchAddress("Electron_energySigmaDown",Electron_energySigmaDown, &b_Electron_energySigmaDown);

  /*
    fChain->SetBranchAddress("nphotons", &nphotons, &b_nphotons);
    fChain->SetBranchAddress("phoe", phoe, &b_phoe);
    fChain->SetBranchAddress("phoeta", phoeta, &b_phoeta);
    fChain->SetBranchAddress("phophi", phophi, &b_phophi);
    fChain->SetBranchAddress("phomvaid", phomvaid, &b_phomvaid);
    fChain->SetBranchAddress("phoe1by9", phoe1by9, &b_phoe1by9);
    fChain->SetBranchAddress("phoe9by25", phoe9by25, &b_phoe9by25);
    fChain->SetBranchAddress("photrkiso", photrkiso, &b_photrkiso);
    fChain->SetBranchAddress("phoemiso", phoemiso, &b_phoemiso);
    fChain->SetBranchAddress("phohadiso", phohadiso, &b_phohadiso);
    fChain->SetBranchAddress("phochhadiso", phochhadiso, &b_phochhadiso);
    fChain->SetBranchAddress("phoneuhadiso", phoneuhadiso, &b_phoneuhadiso);
    fChain->SetBranchAddress("phophoiso", phophoiso, &b_phophoiso);
    fChain->SetBranchAddress("phoPUiso", phoPUiso, &b_phoPUiso);
    fChain->SetBranchAddress("phohadbyem", phohadbyem, &b_phohadbyem);
    fChain->SetBranchAddress("phoietaieta", phoietaieta, &b_phoietaieta);
  */
}

Bool_t top_tagger_sf::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif // #ifdef top_tagger_sf_cxx
