#define top_tagger_sf_cxx
#include "top_tagger_sf.h"

float* top_tagger_sf::Muon_SF(float pt, float eta, string id){

 char name[100];
    float sf =0, sf_stat=0, sf_sys=0;

    if(id == "Tight")
    {
	int eta_bin_id = (h_SF_mu_tight->GetXaxis()->FindBin(eta)>(h_SF_mu_tight->GetNbinsX())) ? h_SF_mu_tight->GetNbinsX() : h_SF_mu_tight->GetXaxis()->FindBin(eta) ;
	int pt_bin_id = (h_SF_mu_tight->GetYaxis()->FindBin(pt)>(h_SF_mu_tight->GetNbinsY())) ? h_SF_mu_tight->GetNbinsY() : h_SF_mu_tight->GetYaxis()->FindBin(pt) ;

      if(eta_bin_id>0 && pt_bin_id>0 ){
        sf = h_SF_mu_tight->GetBinContent(eta_bin_id,pt_bin_id);
		sf_stat = h_SF_stat_mu_tight->GetBinContent(eta_bin_id,pt_bin_id);
		sf_sys = h_SF_sys_mu_tight->GetBinContent(eta_bin_id,pt_bin_id);
      }else{
        sf = 1;
		sf_stat = sf_sys = 1;
        }
      }
      else if(id == "reco")
      {
        int eta_bin_id = (h_SF_mu_reco->GetXaxis()->FindBin(eta)>(h_SF_mu_reco->GetNbinsX())) ? h_SF_mu_reco->GetNbinsX() : h_SF_mu_reco->GetXaxis()->FindBin(eta) ;
        int pt_bin_id = (h_SF_mu_reco->GetYaxis()->FindBin(pt)>(h_SF_mu_reco->GetNbinsY())) ? h_SF_mu_reco->GetNbinsY() : h_SF_mu_reco->GetYaxis()->FindBin(pt) ;
        if(eta_bin_id>0 && pt_bin_id>0 ){
          sf = h_SF_mu_reco->GetBinContent(eta_bin_id,pt_bin_id);
          sf_stat = sf + h_SF_mu_reco->GetBinError(eta_bin_id,pt_bin_id);
          sf_sys = 1;
        }else{
          sf = 1;
          sf_stat = sf_sys = 1;
        }
      }
    else if(id == "loose")
      {
        int eta_bin_id = (h_SF_mu_loose->GetXaxis()->FindBin(eta)>(h_SF_mu_loose->GetNbinsX())) ? h_SF_mu_loose->GetNbinsX() : h_SF_mu_loose->GetXaxis()->FindBin(eta) ;
        int pt_bin_id = (h_SF_mu_loose->GetYaxis()->FindBin(pt)>(h_SF_mu_loose->GetNbinsY())) ? h_SF_mu_loose->GetNbinsY() : h_SF_mu_loose->GetYaxis()->FindBin(pt) ;
        if(eta_bin_id>0 && pt_bin_id>0 ){
          sf = h_SF_mu_loose->GetBinContent(eta_bin_id,pt_bin_id);
          sf_stat = sf + h_SF_mu_loose->GetBinError(eta_bin_id,pt_bin_id);
          sf_sys = 1;
        }else{
          sf = 1;
          sf_stat = sf_sys = 1;
        }
      }
	static float sfvalues[3];
	sfvalues[0] = sf;
	sfvalues[1] = sf_stat;
	sfvalues[2] = sf_sys;

	return sfvalues;
}

float* top_tagger_sf::Electron_SF(float pt, float eta, string id){

	char name[100];
	static float sfvalues[7] = {-100,-100,-100,-100,-100,-100,-100};

    if(id == "reco"){
	int eta_bin_id = (h_SF_el_reco->GetXaxis()->FindBin(eta)>(h_SF_el_reco->GetNbinsX())) ? h_SF_el_reco->GetNbinsX() : h_SF_el_reco->GetXaxis()->FindBin(eta) ;
	int pt_bin_id = (h_SF_el_reco->GetYaxis()->FindBin(pt)>(h_SF_el_reco->GetNbinsY())) ? h_SF_el_reco->GetNbinsY() : h_SF_el_reco->GetYaxis()->FindBin(pt) ;

	float sf, sf_stat, sf_sys;

	if(eta_bin_id>0 && pt_bin_id>0 ){
		sfvalues[0] = h_SF_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[1] = h_SF_statData_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[2] = h_SF_statMC_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[3] = h_SF_altBkgModel_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[4] = h_SF_altSignalModel_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[5] = h_SF_altMCEff_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
		sfvalues[6] = h_SF_altTagSelection_el_reco->GetBinContent(eta_bin_id,pt_bin_id);
	}
	else{
		sfvalues[0] = 1;
		sfvalues[1] = sfvalues[2] = sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = 0;
		}
    }
    else if(id == "Tight") {
        int eta_bin_id = (h_SF_el_tight->GetXaxis()->FindBin(eta)>(h_SF_el_tight->GetNbinsX())) ? h_SF_el_tight->GetNbinsX() : h_SF_el_tight->GetXaxis()->FindBin(eta) ;
        int pt_bin_id = (h_SF_el_tight->GetYaxis()->FindBin(pt)>(h_SF_el_tight->GetNbinsY())) ? h_SF_el_tight->GetNbinsY() : h_SF_el_tight->GetYaxis()->FindBin(pt) ;

        float sf, sf_stat, sf_sys;

        if(eta_bin_id>0 && pt_bin_id>0 ){
          sfvalues[0] = h_SF_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[1] = h_SF_statData_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[2] = h_SF_statMC_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[3] = h_SF_altBkgModel_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[4] = h_SF_altSignalModel_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[5] = h_SF_altMCEff_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[6] = h_SF_altTagSelection_el_tight->GetBinContent(eta_bin_id,pt_bin_id);
      }
      else{
          sfvalues[0] = 1;
          sfvalues[1] = sfvalues[2] = sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = 0;
          }}
          else if(id == "loose") {
        int eta_bin_id = (h_SF_el_loose->GetXaxis()->FindBin(eta)>(h_SF_el_loose->GetNbinsX())) ? h_SF_el_loose->GetNbinsX() : h_SF_el_loose->GetXaxis()->FindBin(eta) ;
        int pt_bin_id = (h_SF_el_loose->GetYaxis()->FindBin(pt)>(h_SF_el_loose->GetNbinsY())) ? h_SF_el_loose->GetNbinsY() : h_SF_el_loose->GetYaxis()->FindBin(pt) ;

        float sf, sf_stat, sf_sys;

        if(eta_bin_id>0 && pt_bin_id>0 ){
          sfvalues[0] = h_SF_el_loose->GetBinContent(eta_bin_id,pt_bin_id);
          sfvalues[1] = 0;
          sfvalues[2] = 0;
          sfvalues[3] = 0;
          sfvalues[4] = 0;
          sfvalues[5] = 0;
          sfvalues[6] = 0;
      }
      else{
          sfvalues[0] = 1;
          sfvalues[1] = sfvalues[2] = sfvalues[3] = sfvalues[4] = sfvalues[5] = sfvalues[6] = 0;
          }
    }

	return sfvalues;
}

float* top_tagger_sf::Get_PU_Weights(int npu){

	int bin_id = h_pileup_data->FindBin(npu);

	static float puweight[3] = {0,0,0};
	if(bin_id>=0 && bin_id<100){
		puweight[0] = h_pileup_data->GetBinContent(bin_id);
		puweight[1] = h_pileup_data_plus->GetBinContent(bin_id);
		puweight[2] = h_pileup_data_minus->GetBinContent(bin_id);
	}
	return puweight;
}

void top_tagger_sf::getmuons(std::vector<Muon> &vmuons, float ptcut=25, float etacut=2.5, int maxsize=njetmx, float dxy_cut=0.2, float dz_cut=0.5)
{

  for(int mu=0; mu<(nmuons); mu++){

    if(Muon_corrected_pt[mu]<ptcut) continue;
    if(fabs(muoneta[mu])>etacut/* || (fabs(muoneta[mu])>1.4442 && fabs(muoneta[mu])<1.566)*/)  continue;

    bool mu_id = Muon_TightID(muonisGL[mu],muonisPF[mu],
			      muonchi[mu],muonhit[mu],muonmst[mu],
			      muontrkvtx[mu],muondz[mu],muonpixhit[mu],muontrklay[mu]);
    //    if(!mu_id) continue;
    bool mu_iso = Muon_Iso_ID(muonpfiso[mu]);

    if(!muonisLoose[mu]) continue;

    //    if(fabs(muontrkvtx[mu])>dxy_cut || fabs(muondz[mu])>dz_cut) continue;

    Muon vmuon;

    vmuon.vetoid_tight = mu_id  && (muonminisoall[mu] < 0.1) /*&& fabs(muontrkvtx[mu])<dxy_cut && fabs(muondz[mu])<dz_cut*/;
    vmuon.pt = Muon_corrected_pt[mu];
    vmuon.eta = muoneta[mu];
    vmuon.phi = muonphi[mu];
    vmuon.mass = 0.105658;
    vmuon.p = muonp[mu];
    vmuon.trkvtx = muontrkvtx[mu];
    vmuon.dz = muondz[mu];
    vmuon.isTRK = muonisTRK[mu];
    vmuon.isGL = muonisGL[mu];
    vmuon.isPF = muonisPF[mu];
    vmuon.isLoose = muonisLoose[mu];
    vmuon.isGoodGL = muonisGoodGL[mu];
    vmuon.isMed = muonisMed[mu];
    vmuon.isMedPr = muonisMedPr[mu];
    vmuon.isTight = muonisTight[mu];
    vmuon.isHighPt = muonisHighPt[mu];
    vmuon.isHighPttrk = muonisHighPttrk[mu];
    vmuon.minisoall = muonminisoall[mu];
    vmuon.chi = muonchi[mu];
    vmuon.posmatch = muonposmatch[mu];
    vmuon.trkink = muontrkink[mu];
    vmuon.segcom = muonsegcom[mu];
    vmuon.hit = muonhit[mu];
    vmuon.mst = muonmst[mu];
    vmuon.pixhit = muonpixhit[mu];
    vmuon.trklay = muontrklay[mu];
    vmuon.valfrac = muonvalfrac[mu];
    vmuon.pfiso = muonpfiso[mu];
    vmuon.mudxy_sv = mudxy_sv[mu];

    vmuon.p4.SetPtEtaPhiM(vmuon.pt,vmuon.eta,vmuon.phi,0.105658);

    vmuon.y = vmuon.p4.Rapidity();
    vmuons.push_back(vmuon);

    if(int(vmuons.size())>=maxsize) break;
  }

  sorted_by_pt(vmuons);

}

void top_tagger_sf::getelectrons(std::vector<Electron> &velectrons, float ptcut=25, float etacut=2.5, int maxsize=njetmx, float dxy_cut=0.05, float dz_cut=0.1)
{

  for(int ie=0; ie<(nelecs); ie++) {

    // if(!elmvaid[ie]) continue;
    //  if(!elmvaid_noIso[ie]) continue;

    if(!elmvaid_noIso_wpLoose[ie]) continue;
    Electron velectron;

    velectron.vetoid_tight = elmvaid_noIso[ie];
    velectron.pt = elpt[ie];
    velectron.eta = eleta[ie];
    velectron.phi = elphi[ie];
    velectron.mass = 0.000511;
    velectron.mvaWP90id = elmvaid[ie];
    velectron.Fallv2WP80 = elmvaid_Fallv2WP80[ie];
    velectron.mvaWP90id_noIso = elmvaid_noIso[ie];
    velectron.Fallv2WP80_noIso = elmvaid_Fallv2WP80_noIso[ie];
    velectron.mvaWPloose = elmvaid_wpLoose[ie];
    velectron.mvaWPloose_noIso = elmvaid_noIso_wpLoose[ie];
    velectron.vetocutid = elcutid_veto[ie];
    velectron.loosecutid = elcutid_loose[ie];
    velectron.medcutid = elcutid_med[ie];
    velectron.tightcutid = elcutid_tight[ie];
    velectron.p = elp[ie];
    velectron.dxy = eldxytrk[ie];
    velectron.dz = eldztrk[ie];
    velectron.pfiso = elpfiso[ie];
    velectron.eldxy_sv = eldxy_sv[ie];
    velectron.supcl_eta = elsupcl_eta[ie];
    velectron.supcl_phi = elsupcl_phi[ie];
    velectron.supcl_rawE = elsupcl_rawE[ie];
    velectron.sigmaieta = elsigmaieta[ie];
    velectron.sigmaiphi = elsigmaiphi[ie];
    velectron.r9full = elr9full[ie];
    velectron.supcl_etaw = elsupcl_etaw[ie];
    velectron.supcl_phiw = elsupcl_phiw[ie];
    velectron.hcaloverecal = elhcaloverecal[ie];
    velectron.cloctftrkn = elcloctftrkn[ie];
    velectron.cloctftrkchi2 = elcloctftrkchi2[ie];
    velectron.e1x5bye5x5 = ele1x5bye5x5[ie];
    velectron.normchi2 = elnormchi2[ie];
    velectron.hitsmiss = elhitsmiss[ie];
    velectron.trkmeasure = eltrkmeasure[ie];
    velectron.ecloverpout = elecloverpout[ie];
    velectron.ecaletrkmomentum = elecaletrkmomentum[ie];
    velectron.deltaetacltrkcalo = eldeltaetacltrkcalo[ie];
    velectron.supcl_preshvsrawe = elsupcl_preshvsrawe[ie];
    velectron.pfisolsumphet = elpfisolsumphet[ie];
    velectron.pfisolsumchhadpt = elpfisolsumchhadpt[ie];
    velectron.pfsiolsumneuhadet = elpfsiolsumneuhadet[ie];
    velectron.etain = eletain[ie];
    velectron.phiin = elphiin[ie];
    velectron.fbrem = elfbrem[ie];
    velectron.eoverp = eleoverp[ie];
    velectron.hovere = elhovere[ie];

    velectron.p4.SetPtEtaPhiE(velectron.pt,velectron.eta,velectron.phi, ele[ie]);

/*    TString str = TString::Format(" el pt correction %f = %f charge %f & %f",Electron_eccalTrkEnergyPostCorr[ie] / ele[ie], Electron_eccalTrkEnergyPostCorr[ie] / velectron.p4.E(),elp[ie],velectron.p4.P()*elcharge[ie]);
    if(gProofServ) gProofServ->SendAsynMessage(str);*/


    velectron.p4 = velectron.p4 * (Electron_eccalTrkEnergyPostCorr[ie] / ele[ie]);
    velectron.pt = velectron.p4.Pt();
    velectron.eta = velectron.p4.Eta();
    velectron.phi = velectron.p4.Phi();
//    velectron.p = velectron.p4.P()*elcharge[ie];

    velectron.y = velectron.p4.Rapidity();

    if (elpt[ie]<ptcut) continue;
    if(fabs(eleta[ie])>etacut /*|| (fabs(eleta[ie])>1.4442 && fabs(eleta[ie])<1.566)*/)  continue;

    velectrons.push_back(velectron);

    if(int(velectrons.size()) >= maxsize) break;

  }

  sorted_by_pt(velectrons);

}

void top_tagger_sf::getLeptons(std::vector<Lepton> &vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, float pt_cut=30)
{ 
  for(unsigned imu=0; imu<vmuons.size(); imu++){
    if(vmuons[imu].pt < pt_cut) continue;
    Lepton vlepton;
    vlepton.pt = vmuons[imu].pt;
    vlepton.eta = vmuons[imu].eta;
    vlepton.y = vmuons[imu].y;
    vlepton.phi = vmuons[imu].phi;
    vlepton.mass = vmuons[imu].mass;
    vlepton.pdgId = (vmuons[imu].p>0) ? -13 : 13;
    vlepton.p4 = vmuons[imu].p4;
    vlepton.indexemu = imu;
    vlepton.AK8_neighbor_index=-1;
    vleptons.push_back(vlepton);
  }
  for(unsigned ie=0; ie<velectrons.size(); ie++){
    if(velectrons[ie].pt < pt_cut) continue;
    Lepton vlepton;
    vlepton.pt = velectrons[ie].pt;
    vlepton.eta = velectrons[ie].eta;
    vlepton.y = velectrons[ie].y;
    vlepton.phi = velectrons[ie].phi;
    vlepton.mass = velectrons[ie].mass;
    vlepton.pdgId = (velectrons[ie].p>0) ? -11 : 11;
    vlepton.p4 = velectrons[ie].p4;
    vlepton.indexemu=ie;
    vlepton.AK8_neighbor_index=-1;
    vleptons.push_back(vlepton);
  }
  sorted_by_pt(vleptons);
}

void top_tagger_sf::convertintolep(Muon vmuon, Lepton &vlepton)
{
    vlepton.pt = vmuon.pt;
    vlepton.eta = vmuon.eta;
    vlepton.y = vmuon.y;
    vlepton.phi = vmuon.phi;
    vlepton.mass = vmuon.mass;
    vlepton.pdgId = (vmuon.p>0) ? -13 : 13;
    vlepton.p4 = vmuon.p4;
}

void top_tagger_sf::convertintolep(Electron velectron, Lepton &vlepton)
{
    vlepton.pt = velectron.pt;
    vlepton.eta = velectron.eta;
    vlepton.y = velectron.y;
    vlepton.phi = velectron.phi;
    vlepton.mass = velectron.mass;
    vlepton.pdgId = (velectron.p>0) ? -11 : 11;
    vlepton.p4 = velectron.p4;
}

void top_tagger_sf::getAK4jets(std::vector<AK4Jet> &Jets, float ptcut=30, float etacut=2.5, bool isMC=false, int maxsize=njetmx)
{

  for(int ijet=0; ijet<(npfjetAK4); ijet++){

    AK4Jet sJet;

    if(!pfjetAK4jetID[ijet]) continue;

    pfjetAK4pt[ijet] *= pfjetAK4JEC[ijet] ;
    pfjetAK4mass[ijet] *= pfjetAK4JEC[ijet];

    if(isMC){
      pfjetAK4pt[ijet] *=max(float(0.0),(1+pfjetAK4reso[ijet])) ;
      pfjetAK4mass[ijet] *= max(float(0.0),(1+pfjetAK4reso[ijet])) ;
      pfjetAK4pt_resoup[ijet] = pfjetAK4pt[ijet]*max(float(0.0),(1+pfjetAK4resoup[ijet]));
      pfjetAK4mass_resoup[ijet] = pfjetAK4mass[ijet]*max(float(0.0),(1+pfjetAK4resoup[ijet]));
      pfjetAK4pt_resodown[ijet] = pfjetAK4pt[ijet]*max(float(0.0),(1+pfjetAK4resodn[ijet]));
      pfjetAK4mass_resodown[ijet] = pfjetAK4mass[ijet]*max(float(0.0),(1+pfjetAK4resodn[ijet]));
    }

    if(fabs(pfjetAK4eta[ijet])>etacut /*|| (fabs(pfjetAK4eta[ijet])>1.4442 && fabs(pfjetAK4eta[ijet])<1.566)*/) continue;
    if(pfjetAK4pt[ijet]<ptcut) continue;

    sJet.pt = pfjetAK4pt[ijet];
    sJet.mass = pfjetAK4mass[ijet];
    sJet.eta = pfjetAK4eta[ijet];
    sJet.y = pfjetAK4y[ijet];
    sJet.phi = pfjetAK4phi[ijet];
    sJet.p4.SetPtEtaPhiM(sJet.pt,sJet.eta,sJet.phi,sJet.mass);
    sJet.pt_resoup = pfjetAK4pt_resoup[ijet];
    sJet.mass_resoup = pfjetAK4mass_resoup[ijet];
    sJet.pt_resodn = pfjetAK4pt_resodown[ijet];
    sJet.mass_resodn = pfjetAK4mass_resodown[ijet];
    sJet.jesup_total = pfjetAK4jesup_total[ijet];
    sJet.jesdn_total = pfjetAK4jesdn_total[ijet];
    sJet.jetid = pfjetAK4jetID[ijet];
    sJet.jetid_tightlepveto = pfjetAK4jetID_tightlepveto[ijet];
    sJet.hadronFlavour = pfjetAK4hadronflav[ijet];
    sJet.partonFlavour = pfjetAK4partonflav[ijet];
    sJet.btag_DeepFlav = pfjetAK4btag_DeepFlav[ijet];
    sJet.btag_DeepCSV = pfjetAK4btag_DeepCSV[ijet];
    sJet.puid = pfjetAK4PUID[ijet];
    sJet.qgl = pfjetAK4qgl[ijet];
    sJet.AK8_neighbor_index = -1;
    Jets.push_back(sJet);

    if(int(Jets.size())>=maxsize) break;

  }

  sorted_by_pt(Jets);

}

void top_tagger_sf::getAK8jets(std::vector<AK8Jet> &LJets, float ptcut=200, float etacut=2.5, bool isMC=false, int maxsize=njetmxAK8)
{

  for(int ijet=0; ijet<(npfjetAK8); ijet++){

    if(!pfjetAK8jetID[ijet]) continue;

    pfjetAK8pt[ijet] *= pfjetAK8JEC[ijet] ;
    pfjetAK8mass[ijet] *= pfjetAK8JEC[ijet];

    if(isMC){

      pfjetAK8pt[ijet] *= max(float(0.0),(1+pfjetAK8reso[ijet])) ;
      pfjetAK8mass[ijet] *= max(float(0.0),(1+pfjetAK8reso[ijet])) ;

      pfjetAK8pt_resoup[ijet] = pfjetAK8pt[ijet]*max(float(0.0),(1+pfjetAK8resoup[ijet]));
      pfjetAK8mass_resoup[ijet] = pfjetAK8mass[ijet]*max(float(0.0),(1+pfjetAK8resoup[ijet]));
      pfjetAK8pt_resodown[ijet] = pfjetAK8pt[ijet]*max(float(0.0),(1+pfjetAK8resodn[ijet]));
      pfjetAK8mass_resodown[ijet] = pfjetAK8mass[ijet]*max(float(0.0),(1+pfjetAK8resodn[ijet]));
    }

    if(fabs(pfjetAK8eta[ijet])>etacut /*|| (fabs(pfjetAK8eta[ijet])>1.4442 && fabs(pfjetAK8eta[ijet])<1.566)*/) continue;
    if(pfjetAK8pt[ijet] < ptcut) continue;

    AK8Jet LJet;

    LJet.jetid = pfjetAK8jetID[ijet];
    LJet.jetid_tightlepveto = pfjetAK8jetID_tightlepveto[ijet];
    LJet.pt = pfjetAK8pt[ijet];
    LJet.eta = pfjetAK8eta[ijet];
    LJet.mass = pfjetAK8mass[ijet];
    LJet.phi = pfjetAK8phi[ijet];
    LJet.y = pfjetAK8y[ijet];
    LJet.p4.SetPtEtaPhiM(LJet.pt,LJet.eta,LJet.phi,LJet.mass);
    LJet.pt_resoup = pfjetAK8pt_resoup[ijet];
    LJet.mass_resoup = pfjetAK8mass_resoup[ijet];
    LJet.pt_resodn = pfjetAK8pt_resodown[ijet];
    LJet.mass_resodn = pfjetAK8mass_resodown[ijet];
    LJet.jesup_total = pfjetAK8jesup_total[ijet];
    LJet.jesdn_total = pfjetAK8jesdn_total[ijet];

    LJet.chrad = pfjetAK8chrad[ijet];
    LJet.tau21 = pfjetAK8tau2[ijet]*1./max(float(1.e-6),pfjetAK8tau1[ijet]);
    LJet.tau32 = pfjetAK8tau3[ijet]*1./max(float(1.e-6),pfjetAK8tau2[ijet]);
    LJet.DeepTag_TvsQCD = pfjetAK8DeepTag_TvsQCD[ijet];
    LJet.DeepTag_WvsQCD = pfjetAK8DeepTag_WvsQCD[ijet];
    LJet.DeepTag_ZvsQCD = pfjetAK8DeepTag_ZvsQCD[ijet];
    LJet.PNet_TvsQCD = pfjetAK8PNet_TvsQCD[ijet];
    LJet.PNet_WvsQCD = pfjetAK8PNet_WvsQCD[ijet];
    LJet.btag_DeepCSV = pfjetAK8btag_DeepCSV[ijet];

    LJet.CHF = pfjetAK8CHF[ijet];
    LJet.NHF = pfjetAK8NHF[ijet];
    LJet.CEMF = pfjetAK8CEMF[ijet];
    LJet.NEMF = pfjetAK8NEMF[ijet];
    LJet.MUF = pfjetAK8MUF[ijet];
    LJet.PHF = pfjetAK8PHF[ijet];
    pfjetAK8HadF[ijet] = (pfjetAK8NHF[ijet]+pfjetAK8CHF[ijet]);
    LJet.NHadF = (1.- pfjetAK8HadF[ijet]);

    LJet.EMF = (pfjetAK8NEMF[ijet]+pfjetAK8CEMF[ijet]);
    LJet.EEM = pfjetAK8EEM[ijet];

    LJet.neucons = pfjetAK8Neucons[ijet];
    LJet.chcons = pfjetAK8Chcons[ijet];

    LJet.neuemfrac = (pfjetAK8NEMF[ijet]*1.)/max(float(1.e-6),(pfjetAK8NEMF[ijet]+pfjetAK8CEMF[ijet]));
    LJet.neunhadfrac = (pfjetAK8NEMF[ijet]*1.)/max(1.e-6,(1.- pfjetAK8HadF[ijet]));

    LJet.sdmass = pfjetAK8sdmass[ijet];
    LJet.sub1pt = pfjetAK8sub1pt[ijet];
    LJet.sub1eta = pfjetAK8sub1eta[ijet];
    LJet.sub1phi = pfjetAK8sub1phi[ijet];
    LJet.sub1mass = pfjetAK8sub1mass[ijet];
    LJet.sub1btag = pfjetAK8sub1btag[ijet];
    LJet.sub1hadfrac = pfjetAK8sub1chhadfrac[ijet]+pfjetAK8sub1neuhadfrac[ijet];
    LJet.sub1emfrac = pfjetAK8sub1emfrac[ijet];

    LJet.sub2pt = pfjetAK8sub2pt[ijet];
    LJet.sub2eta = pfjetAK8sub2eta[ijet];
    LJet.sub2phi = pfjetAK8sub2phi[ijet];
    LJet.sub2mass = pfjetAK8sub2mass[ijet];
    LJet.sub2btag = pfjetAK8sub2btag[ijet];
    LJet.sub2hadfrac = pfjetAK8sub2chhadfrac[ijet]+pfjetAK8sub2neuhadfrac[ijet];
    LJet.sub2emfrac = pfjetAK8sub2emfrac[ijet];

    LJet.subbtag = max(pfjetAK8sub1btag[ijet],pfjetAK8sub2btag[ijet]);

    LJet.subhaddiff = pfjetAK8subhaddiff[ijet];
    LJet.subemdiff = pfjetAK8subemdiff[ijet];
    LJet.subptdiff = pfjetAK8subptdiff[ijet];

    LJet.elinsubpt = pfjetAK8elinsubpt[ijet];
    LJet.elinsubeta = pfjetAK8elinsubeta[ijet];
    LJet.elinsubphi = pfjetAK8elinsubphi[ijet];
        
    LJet.elinsubjpt = pfjetAK8elinsubjpt[ijet];
    LJet.elinsubjeta = pfjetAK8elinsubjeta[ijet];
    LJet.elinsubjphi = pfjetAK8elinsubjphi[ijet];
    LJet.elinsubjmass = pfjetAK8elinsubjmass[ijet];

    LJet.muinsubpt = pfjetAK8muinsubpt[ijet];
    LJet.muinsubeta = pfjetAK8muinsubeta[ijet];
    LJet.muinsubphi = pfjetAK8muinsubphi[ijet];
    LJet.muinsubjpt = pfjetAK8muinsubjpt[ijet];
    LJet.muinsubjeta = pfjetAK8muinsubjeta[ijet];
    LJet.muinsubjphi = pfjetAK8muinsubjphi[ijet];
    LJet.muinsubjmass = pfjetAK8muinsubjmass[ijet];

    LJet.muinsubI0 = pfjetAK8muinsubI0[ijet];
    LJet.muinsubInear = pfjetAK8muinsubInear[ijet];
    LJet.muinsubIfar = pfjetAK8muinsubIfar[ijet];

    LJet.haselectron = LJet.hasmuon = LJet.hastau = LJet.hasqg = LJet.hasb = LJet.hasleptop = LJet.hashadtop = LJet.hastop = LJet.hasmatchmu = LJet.hasmatche = false;
    LJet.hasleptop_alldecay = LJet.hashadtop_alldecay = false;

    LJet.haspfelectron = LJet.haspfmuon = false;
    LJet.matchAK4deepb = LJet.re_tvsb = LJet.rmu_tvsb = -99;

    /*TString str;
     *
     *    str = TString::Format(" %f %f %f %f %f %f ",LJet.muinsubpt , LJet.muinsubjpt, LJet.muinsubeta,LJet.muinsubphi,LJet.y,LJet.phi);
     *    if(gProofServ) gProofServ->SendAsynMessage(str);*/

    if(LJet.elinsubpt > 0. && LJet.elinsubjpt > 0.){
      if (delta2R(LJet.elinsubeta,LJet.elinsubphi,LJet.y,LJet.phi) < 0.8) {
	LJet.haspfelectron = true;
      }
    }

    if (LJet.muinsubpt > 0. && LJet.muinsubjpt > 0.) {
      if (delta2R(LJet.muinsubeta,LJet.muinsubphi,LJet.y,LJet.phi) < 0.8) {
	LJet.haspfmuon = true;
      }
    }

    LJets.push_back(LJet);

    if(int(LJets.size())>=maxsize) break;
  }

  sorted_by_pt(LJets);

}

void top_tagger_sf::getAK4genjets(std::vector<AK4GenJet> &Jets, float ptcut=30, float etacut=2.5, int maxsize=njetmx)
{

  for(int ijet=0; ijet<(ngenjetAK4); ijet++){

    AK4GenJet sJet;

    if(fabs(genjetAK4eta[ijet])>etacut /*|| (fabs(pfjetAK4eta[ijet])>1.4442 && fabs(pfjetAK4eta[ijet])<1.566)*/) continue;
    if(genjetAK4pt[ijet]<ptcut) continue;

    sJet.pt = genjetAK4pt[ijet];
    sJet.mass = genjetAK4mass[ijet];
    sJet.eta = genjetAK4eta[ijet];
    sJet.phi = genjetAK4phi[ijet];
    sJet.p4.SetPtEtaPhiM(sJet.pt,sJet.eta,sJet.phi,sJet.mass);
    sJet.y = sJet.p4.Rapidity();
    sJet.hadronFlavor = genjetAK4hadronflav[ijet];
    sJet.partonFlavor = genjetAK4partonflav[ijet];

    Jets.push_back(sJet);

    if(int(Jets.size())>=maxsize) break;

  }

  sorted_by_pt(Jets);

}

/*void top_tagger_sf::LeptonJet_cleaning(std::vector <Electron> velectrons,std::vector <Muon> vmuons, float dR_cut=0.4)
 * {
 * for(int imu=0; imu<(int)vmuons.size(); imu++)
 *    {
 *    TLorentzVector j_mom,clj_mom, mu_mom;
 *    int matchjet=-1;
 *    float delrmin=0.4;
 *    for(int ijet=0; ijet<npfjetAK4; ijet++)
 *    {
 *    float delrtmp = delta2R(pfjetAK4eta[ijet], pfjetAK4phi[ijet], vmuons[imu].eta, vmuons[imu].phi);
 *    hist_prptangle[1]->Fill(pfjetAK4pt[ijet], delrtmp,weight);
 *    if(delrtmp<delrmin)
 *        {
 *        matchjet=ijet;
 *        delrmin = delrtmp;
 *        }
 *    }
 *    if(matchjet>=0)
 *    {
 *    j_mom.SetPtEtaPhiM(pfjetAK4pt[matchjet],pfjetAK4eta[matchjet],pfjetAK4phi[matchjet],pfjetAK4mass[matchjet]);
 *    mu_mom.SetPtEtaPhiM(vmuons[imu].pt, vmuons[imu].eta, vmuons[imu].phi, 0.105658);
 *    clj_mom = j_mom - mu_mom;
 *    pfjetAK4pt[matchjet] = clj_mom.Pt();
 *    pfjetAK4mass[matchjet] = clj_mom.M();
 *    pfjetAK4eta[matchjet] = clj_mom.Eta();
 *    pfjetAK4y[matchjet] = clj_mom.Rapidity();
 *    pfjetAK4phi[matchjet] = clj_mom.Phi();
 *    pfjetAK4inleppt[matchjet] += mu_mom.Pt();
 *    pfjetAK4delrlep[matchjet] = delrmin;
 *    }
 *    }
 *
 * for(int ie=0; ie<(int)velectrons.size(); ie++) {
 *    {
 *    TLorentzVector j_mom,clj_mom, el_mom;
 *    int matchjet=-1;
 *    float delrmin=0.4;
 *    for(int ijet=0; ijet<npfjetAK4; ijet++)
 *    {
 *    float delrtmp = delta2R(pfjetAK4eta[ijet],pfjetAK4phi[ijet], velectrons[ie].y, velectrons[ie].phi);
 *    hist_prptangle[2]->Fill(pfjetAK4pt[ijet], delrtmp,weight);
 *    if(delrtmp<delrmin)
 *        {
 *        delrmin = delrtmp;
 *        matchjet=ijet;
 *        }
 *    }
 *    if(matchjet>=0)
 *    {
 *    j_mom.SetPtEtaPhiM(pfjetAK4pt[matchjet],pfjetAK4eta[matchjet],pfjetAK4phi[matchjet],pfjetAK4mass[matchjet]);
 *    el_mom.SetPtEtaPhiM(velectrons[ie].pt, velectrons[ie].eta, velectrons[ie].phi, 0.000511);
 *    clj_mom = j_mom - el_mom;
 *    pfjetAK4pt[matchjet] = clj_mom.Pt();
 *    pfjetAK4mass[matchjet] = clj_mom.M();
 *    pfjetAK4eta[matchjet] = clj_mom.Eta();
 *    pfjetAK4y[matchjet] = clj_mom.Rapidity();
 *    pfjetAK4phi[matchjet] = clj_mom.Phi();
 *    pfjetAK4inleppt[matchjet] += el_mom.Pt();
 *    pfjetAK4delrlep[matchjet] = delrmin;
 *    }
 *    }
 * }
 * }*/

void top_tagger_sf::getPartons(std::vector<GenParton> &GenPartons, int maxsize=npartmx)
{

  for(int igen=0; igen<(ngenparticles); igen++){

    GenParton parton;

    parton.pt = genpartpt[igen];
    parton.eta = genparteta[igen];
    parton.phi = genpartphi[igen];
    parton.mass = genpartm[igen];
    parton.p4.SetPtEtaPhiM(parton.pt,parton.eta,parton.phi,parton.mass);
    parton.y = parton.p4.Rapidity();
    parton.status = genpartstatus[igen];
    parton.pdgId = genpartpdg[igen];
    parton.mompdgId = genpartmompdg[igen];
    parton.grmompdgId = genpartgrmompdg[igen];

    parton.fromhard = genpartfromhard[igen];
    parton.fromhardbFSR = genpartfromhardbFSR[igen];
    parton.isPromptFinalState = genpartisPromptFinalState[igen];
    parton.isLastCopyBeforeFSR = genpartisLastCopyBeforeFSR[igen];

    GenPartons.push_back(parton);

    if(int(GenPartons.size())>=maxsize) break;

  }

  //	GenPartons.sorted_by_pt();

}

void top_tagger_sf::getLHEParticles(std::vector<LHEparticle> &lheparticles, int maxsize=npartmx)
{

  for(int igen=0; igen<(nLHEparticles); igen++){

    LHEparticle lheparticle;

    lheparticle.pt = LHEpartpt[igen];
    lheparticle.eta = LHEparteta[igen];
    lheparticle.phi = LHEpartphi[igen];
    lheparticle.mass = LHEpartm[igen];
    lheparticle.p4.SetPtEtaPhiM(lheparticle.pt,lheparticle.eta,lheparticle.phi,lheparticle.mass);
    lheparticle.y = lheparticle.p4.Rapidity();
    lheparticle.pdgId = LHEpartpdg[igen];
    lheparticles.push_back(lheparticle);

    if(int(lheparticles.size())>=maxsize) break;

  }

}

void top_tagger_sf::getLHETops(std::vector<GenParton> &LHETops, std::vector<GenParton> GenPartons)
{

  for(auto & part: GenPartons){

    if(abs(part.status)!=22) continue;
    if(!(part.fromhard)) continue;
    if(abs(part.pdgId)!=6) continue;
    LHETops.push_back(part);

  }
}

void top_tagger_sf::getGENTops(vector<TopQuark> &gentops, vector<GenParton> genpartons)  // with daughters after shower
{     
  vector<GenParton> W_dau;
  vector<GenParton> t_bp;

  for(unsigned igen=0; igen<genpartons.size(); igen++){

    if(!(genpartons[igen].status==23 || genpartons[igen].status==1)) continue;
    if(!(genpartons[igen].fromhard)) continue;

    if(!((abs(genpartons[igen].pdgId)>=1 && abs(genpartons[igen].pdgId)<=5)||(abs(genpartons[igen].pdgId)>=11 && abs(genpartons[igen].pdgId)<=16))) continue;
    if(!(abs(genpartons[igen].mompdgId)==6 || abs(genpartons[igen].mompdgId)==24)) continue;
    if(abs(genpartons[igen].mompdgId)==24 && (abs(genpartons[igen].grmompdgId)!=6 && abs(genpartons[igen].grmompdgId)!=24 )) continue;

    if(abs(genpartons[igen].pdgId)>=1 && abs(genpartons[igen].pdgId)<=5 && abs(genpartons[igen].mompdgId)==24 )   {  W_dau.push_back(genpartons[igen]); }
    if(abs(genpartons[igen].pdgId)>=11 && abs(genpartons[igen].pdgId)<=16 && abs(genpartons[igen].mompdgId)==24 ) {  W_dau.push_back(genpartons[igen]); }
    if(abs(genpartons[igen].pdgId)==5 && abs(genpartons[igen].mompdgId)==6) {  t_bp.push_back(genpartons[igen]); }
  }

  for(unsigned ipart=0; ipart<W_dau.size(); ipart++){

    if(!((abs(W_dau[ipart].pdgId)==11 || abs(W_dau[ipart].pdgId)==13 || abs(W_dau[ipart].pdgId)==15 || abs(W_dau[ipart].pdgId)==1 || abs(W_dau[ipart].pdgId)==3 || abs(W_dau[ipart].pdgId)==5))) continue;

    unsigned partner = -1;
    unsigned match_b = -1;

    for(unsigned jpart=0; jpart<W_dau.size(); jpart++){
      if( (W_dau[ipart].mompdgId==W_dau[jpart].mompdgId) && (W_dau[ipart].grmompdgId==W_dau[jpart].grmompdgId) && (W_dau[ipart].pdgId*W_dau[jpart].pdgId<0)){
	partner = jpart;
	break;
      }
    }

    for(unsigned ib=0; ib<t_bp.size(); ib++){
      if(t_bp[ib].pdgId*W_dau[ipart].mompdgId>0){
	match_b = ib;
	break;
      }
    }

    GenParton q1, q2, b;

    q1 = W_dau[ipart];

    if(int(partner)>=0 && partner<W_dau.size()){
      q2 = W_dau[partner];
    }
    if(int(match_b)>=0 && match_b<t_bp.size()){
      b = t_bp[match_b];
    }

    if(int(partner)>=0 && partner<W_dau.size() && int(match_b)>=0 && match_b<t_bp.size()){

      TopQuark topQ;

      topQ.p4 = (b.p4 + q1.p4 + q2.p4);
      topQ.daughter.push_back(q1);
      topQ.daughter.push_back(q2);
      topQ.daughter.push_back(b);

      gentops.push_back(topQ);
    }
  }
}



void top_tagger_sf::TopAssignment_toJet(std::vector<AK8Jet> &LJets, std::vector<GenParton> lhetops, std::vector<TopQuark> gentops)
{
  for(unsigned ijet=0; ijet<LJets.size(); ijet++){
    for(unsigned itop=0; itop<lhetops.size(); itop++)
      {
	if(delta2R(LJets[ijet].p4,lhetops[itop].p4)<0.6){
	  LJets[ijet].hastop = true;
	  break;
	}
      }

    for(unsigned itop=0; itop<gentops.size(); itop++){

      if(abs(gentops[itop].daughter[0].pdgId)==11 || abs(gentops[itop].daughter[0].pdgId)==13 || abs(gentops[itop].daughter[0].pdgId)==15){

	if(delta2R(LJets[ijet].p4,gentops[itop].p4)<0.8){
	  LJets[ijet].hasleptop = true;
	}

	if((delta2R(LJets[ijet].p4,gentops[itop].daughter[0].p4)<0.8)
	   && (delta2R(LJets[ijet].p4,gentops[itop].daughter[2].p4)<0.8))
	  {
	    LJets[ijet].hasleptop_alldecay = true;
	  }

      }

      else if(abs(gentops[itop].daughter[0].pdgId)==1 || abs(gentops[itop].daughter[0].pdgId)==3 || abs(gentops[itop].daughter[0].pdgId)==5){

	if(delta2R(LJets[ijet].p4,gentops[itop].p4)<0.8){
	  LJets[ijet].hashadtop = true;
	}

	if((delta2R(LJets[ijet].p4,gentops[itop].daughter[0].p4)<0.8) && (delta2R(LJets[ijet].p4,gentops[itop].daughter[1].p4)<0.8)
	   && (delta2R(LJets[ijet].p4,gentops[itop].daughter[2].p4)<0.8))
	  {
	    LJets[ijet].hashadtop_alldecay = true;
	  }

      }
    }

  }//ijet

}

void top_tagger_sf::AssignGen(std::vector<AK8Jet> &LJets, std::vector<GenParton> GenPartons){

  for(auto & ljet: LJets){

    for(auto & part: GenPartons){

      if(!(part.status==23 || part.status==1)) continue;
      if(!(part.fromhard)) continue;

      float delr=delta2R(ljet.p4,part.p4);
      if(delr<0.8){
        if(abs(part.pdgId)==11){
	  ljet.haselectron  = true;

        }

        if(abs(part.pdgId)==13){
	  ljet.hasmuon  = true;

        }

        if(abs(part.pdgId)==15){
	  ljet.hastau = true;

        }
      }
      if(delr<0.8){
	if(((abs(part.pdgId)>=1 && abs(part.pdgId)<5) || abs(part.pdgId)==21) ){
	  ljet.hasqg  = true;

        }

        if(abs(part.pdgId)==5 ){
	  ljet.hasb  = true;

        }
      }
    }
  }
}

bool top_tagger_sf::isBJet(AK4Jet jet, float btag_cut)
{
  if(jet.btag_DeepFlav >= btag_cut){
    return true;
  }
  else{
    return false;
  }
}


void top_tagger_sf::ReadTagger(AK8Jet &LJet, std::vector<Lepton> vleptons, std::vector<Muon> vmuons, std::vector<Electron> velectrons, TMVA::Reader *reader_electron, TMVA::Reader *reader_muon){

    //Muon
    in_mupfjetAK8NHadF = -999;
    in_mupfjetAK8neunhadfrac = -999;
    in_mupfjetAK8subhaddiff = -999;
    in_mupfjetAK8tau21 = -999;
    in_mupfjetAK8chrad = -999;
    in_mupfjetAK8sdmass = -999;
    in_pfjetAK8matchedmuonchi = -999;
    in_pfjetAK8matchedmuonposmatch = -999;
    in_pfjetAK8matchedmuontrkink = -999;
    in_pfjetAK8matchedmuonsegcom = -999;
    in_pfjetAK8matchedmuonhit = -999;
    in_pfjetAK8matchedmuonmst = -999;
    in_pfjetAK8matchedmuontrkvtx = -999;
    in_pfjetAK8matchedmuondz = -999;
    in_pfjetAK8matchedmuonpixhit = -999;
    in_pfjetAK8matchedmuontrklay = -999;
    in_pfjetAK8matchedmuonvalfrac = -999;
    in_pfjetAK8muinsubptrat = -999;
    in_pfjetAK8muinsubmassrat = -999;
    in_pfjetAK8muinsubinvmass = -999;
    in_pfjetAK8muinsubIfarbyI0 = -999;
    in_pfjetAK8muinsubInearbyI0 = -999;
    //electron
    in_pfjetAK8NHadF = -999;
    in_pfjetAK8neunhadfrac = -999;
    in_pfjetAK8subhaddiff = -999;
    in_pfjetAK8tau21 = -999;
    in_pfjetAK8chrad = -999;
    in_pfjetAK8sdmass = -999;
    in_pfjetAK8matchedelcleta = -999;
    in_pfjetAK8matchedelpt = -999;
    in_pfjetAK8matchedelsigmaieta = -999;
    in_pfjetAK8matchedelsigmaiphi = -999;
    in_pfjetAK8matchedelr9full = -999;
    in_pfjetAK8matchedelsupcl_etaw = -999;
    in_pfjetAK8matchedelsupcl_phiw = -999;
    in_pfjetAK8matchedelhcaloverecal = -999;
    in_pfjetAK8matchedelcloctftrkn = -999;
    in_pfjetAK8matchedelcloctftrkchi2 = -999;
    in_pfjetAK8matchedele1x5bye5x5 = -999;
    in_pfjetAK8matchedelnormchi2 = -999;
    in_pfjetAK8matchedelhitsmiss = -999;
    in_pfjetAK8matchedeltrkmeasure = -999;
    in_pfjetAK8matchedelecloverpout = -999;
    in_pfjetAK8matchedelecaletrkmomentum = -999;
    in_pfjetAK8matchedeldeltaetacltrkcalo = -999;
    in_pfjetAK8matchedelsupcl_preshvsrawe = -999;
    in_pfjetAK8matchedelpfisolsumphet = -999;
    in_pfjetAK8matchedelpfisolsumchhadpt = -999;
    in_pfjetAK8matchedelpfisolsumneuhadet = -999;
    in_pfjetAK8matchedeletain = -999;
    in_pfjetAK8matchedelphiin = -999;
    in_pfjetAK8matchedelfbrem = -999;
    in_pfjetAK8matchedeleoverp = -999;
    in_pfjetAK8matchedelhovere = -999;
    in_pfjetAK8matchedeldxy_sv = -999;
    in_pfjetAK8matchedelRho = (Rho);

    LJet.hasmatche = false;
    LJet.hasmatchmu = false;

    if (int(vleptons.size()) >0) {

      int lepindex = LJet.match_lepton_index;
      if (abs(vleptons[lepindex].pdgId)==11) { //Electron

	in_pfjetAK8NHadF = LJet.NHadF;
	in_pfjetAK8neunhadfrac = LJet.neunhadfrac;
	in_pfjetAK8subhaddiff = LJet.subhaddiff;
	in_pfjetAK8tau21 = LJet.tau21;
	in_pfjetAK8chrad = LJet.chrad;
	in_pfjetAK8sdmass = LJet.sdmass;

	LJet.hasmatche = true;
	int iel = vleptons[lepindex].indexemu;

	in_pfjetAK8matchedelcleta = velectrons[iel].supcl_eta; //elsupcl_eta[nearest];
	in_pfjetAK8matchedelpt = fabs(velectrons[iel].pt); //elpt[nearest]);
	in_pfjetAK8matchedelsigmaieta = velectrons[iel].sigmaieta; // elsigmaieta[nearest];
	in_pfjetAK8matchedelsigmaiphi = velectrons[iel].sigmaiphi; // elsigmaiphi[nearest];
	in_pfjetAK8matchedelr9full = velectrons[iel].r9full; // elr9full[nearest];
	in_pfjetAK8matchedelsupcl_etaw = velectrons[iel].supcl_etaw; // elsupcl_etaw[nearest];
	in_pfjetAK8matchedelsupcl_phiw = velectrons[iel].supcl_phiw; // elsupcl_phiw[nearest];
	in_pfjetAK8matchedelhcaloverecal = velectrons[iel].hcaloverecal; // elhcaloverecal[nearest];
	in_pfjetAK8matchedelcloctftrkn = velectrons[iel].cloctftrkn; // elcloctftrkn[nearest];
	in_pfjetAK8matchedelcloctftrkchi2 = velectrons[iel].cloctftrkchi2; // elcloctftrkchi2[nearest];
	in_pfjetAK8matchedele1x5bye5x5 = velectrons[iel].e1x5bye5x5; // ele1x5bye5x5[nearest];
	in_pfjetAK8matchedelnormchi2 = velectrons[iel].normchi2; // elnormchi2[nearest];
	in_pfjetAK8matchedelhitsmiss = velectrons[iel].hitsmiss; // elhitsmiss[iel];
	in_pfjetAK8matchedeltrkmeasure = velectrons[iel].trkmeasure; // eltrkmeasure[nearest];
	in_pfjetAK8matchedelecloverpout = velectrons[iel].ecloverpout; // elecloverpout[nearest];
	in_pfjetAK8matchedelecaletrkmomentum = velectrons[iel].ecaletrkmomentum; //elecaletrkmomentum[nearest];
	in_pfjetAK8matchedeldeltaetacltrkcalo = velectrons[iel].deltaetacltrkcalo; // eldeltaetacltrkcalo[nearest];
	in_pfjetAK8matchedelsupcl_preshvsrawe = velectrons[iel].supcl_preshvsrawe; // elsupcl_preshvsrawe[nearest];
	in_pfjetAK8matchedelpfisolsumphet = velectrons[iel].pfisolsumphet; // elpfisolsumphet[nearest];
	in_pfjetAK8matchedelpfisolsumchhadpt = velectrons[iel].pfisolsumchhadpt; // elpfisolsumchhadpt[nearest];
	in_pfjetAK8matchedelpfisolsumneuhadet = velectrons[iel].pfsiolsumneuhadet; //elpfsiolsumneuhadet[nearest];
	in_pfjetAK8matchedeletain = velectrons[iel].etain; // eletain[nearest];
	in_pfjetAK8matchedelphiin = velectrons[iel].phiin; // elphiin[nearest];
	in_pfjetAK8matchedelfbrem = velectrons[iel].fbrem; //elfbrem[nearest];
	in_pfjetAK8matchedeleoverp = velectrons[iel].eoverp; // eleoverp[nearest];
	in_pfjetAK8matchedelhovere = velectrons[iel].hovere; //elhovere[nearest];
	in_pfjetAK8matchedeldxy_sv = velectrons[iel].eldxy_sv;

	LJet.re_tvsb = reader_electron->EvaluateMVA("BDTG method");

      } else { // Muon

	in_mupfjetAK8NHadF = LJet.NHadF;
	in_mupfjetAK8neunhadfrac = LJet.neunhadfrac;
	in_mupfjetAK8subhaddiff = LJet.subhaddiff;
	in_mupfjetAK8tau21 = LJet.tau21;
	in_mupfjetAK8chrad = LJet.chrad;
	in_mupfjetAK8sdmass = LJet.sdmass;

	LJet.hasmatchmu = true;
	int imu = vleptons[lepindex].indexemu;
	in_pfjetAK8matchedmuonchi = vmuons[imu].chi; //GMA //muonchi[nearestmu];
	in_pfjetAK8matchedmuonposmatch = vmuons[imu].posmatch; //muonposmatch[nearestmu];
	in_pfjetAK8matchedmuontrkink = vmuons[imu].trkink; // muontrkink[nearestmu];
	in_pfjetAK8matchedmuonsegcom = vmuons[imu].segcom; // muonsegcom[nearestmu];
	in_pfjetAK8matchedmuonhit = vmuons[imu].hit; // muonhit[nearestmu];
	in_pfjetAK8matchedmuonmst = vmuons[imu].mst; //muonmst[nearestmu];
	in_pfjetAK8matchedmuontrkvtx = vmuons[imu].trkvtx; //muontrkvtx[nearestmu];
	in_pfjetAK8matchedmuondz = vmuons[imu].dz; //muondz[nearestmu];
	in_pfjetAK8matchedmuonpixhit = vmuons[imu].pixhit; // muonpixhit[nearestmu];
	in_pfjetAK8matchedmuontrklay = vmuons[imu].trklay; // muontrklay[nearestmu];
	in_pfjetAK8matchedmuonvalfrac = vmuons[imu].valfrac; // muonvalfrac[nearestmu];
	TLorentzVector lep, lepbj, bj;
	lep.SetPtEtaPhiM(LJet.muinsubpt, LJet.muinsubeta, LJet.muinsubphi, 0.105658);
	bj.SetPtEtaPhiM(LJet.muinsubjpt, LJet.muinsubjeta, LJet.muinsubjphi, LJet.muinsubjmass);
	lepbj = lep + bj;

	in_pfjetAK8muinsubptrat = (lep.Pt())/(max(1.e-6,lepbj.Pt()));
	in_pfjetAK8muinsubmassrat = (bj.M())/(max(1.e-6,lepbj.M()));
	in_pfjetAK8muinsubinvmass = lepbj.M();
	in_pfjetAK8muinsubIfarbyI0 = LJet.muinsubIfar/max(float(1.e-6),LJet.muinsubI0);
	in_pfjetAK8muinsubInearbyI0 = LJet.muinsubInear/max(float(1.e-6),LJet.muinsubI0);

	LJet.re_tvsb = reader_muon->EvaluateMVA("BDTG method");
      } //else of if (vleptons[ij].lepton_id==1)
    } // if (vleptons.size() >=ij)
}

void top_tagger_sf::Match_trigger(vector<Single_Trigger> vsinglelep_trig,vector<Double_Trigger> vdoublelep_trig, vector<std::pair<int,TLorentzVector>> TrigRefObj,Lepton lepcand_1, Lepton lepcand_2, vector<AK4Jet> Jets, bool &trig_threshold_pass, bool &trig_matching_pass, vector<TH1D*> &hist_init,vector<TH2D*> &hist2d)
{

  if(vdoublelep_trig.size()<1 && vsinglelep_trig.size()<1) {
    return;
  }

  else{

    // checking if offline objects passed trigger thresholds //

    bool double_trig_pass(false), single_trig_pass(false);
    bool any_double_hlt_pass = false;

    if(vdoublelep_trig.size()>0){
      for(unsigned ihlt=0; ihlt<vdoublelep_trig.size(); ihlt++){
	if (vdoublelep_trig[ihlt].double_hlts){ any_double_hlt_pass = true; }
	if (vdoublelep_trig[ihlt].double_hlts && ((lepcand_1.pt>vdoublelep_trig[ihlt].double_pt_cuts[0] && vdoublelep_trig[ihlt].double_pids[0]==abs(lepcand_1.pdgId) && lepcand_2.pt>vdoublelep_trig[ihlt].double_pt_cuts[1] && vdoublelep_trig[ihlt].double_pids[1]==abs(lepcand_2.pdgId)) || (lepcand_1.pt>vdoublelep_trig[ihlt].double_pt_cuts[1] && vdoublelep_trig[ihlt].double_pids[1]==abs(lepcand_1.pdgId) && lepcand_2.pt>vdoublelep_trig[ihlt].double_pt_cuts[0] && vdoublelep_trig[ihlt].double_pids[0]==abs(lepcand_1.pdgId)))) {
	  trig_threshold_pass = true;
	  double_trig_pass = true;
	  break;
	}
      }
    }

    int fired_single_trig = -1;

    if(!any_double_hlt_pass){
      if(vsinglelep_trig.size()>0){
	for(unsigned ihlt=0; ihlt<vsinglelep_trig.size(); ihlt++){
	  if (vsinglelep_trig[ihlt].single_hlts && ((lepcand_1.pt>vsinglelep_trig[ihlt].single_pt_cuts && abs(lepcand_1.pdgId)==vsinglelep_trig[ihlt].single_pids) || (lepcand_2.pt>vsinglelep_trig[ihlt].single_pt_cuts && abs(lepcand_2.pdgId)==vsinglelep_trig[ihlt].single_pids)))
	    {
	      if(vsinglelep_trig[ihlt].single_other_pt_cuts>0 && Jets.size()>0 && Jets[0].pt>vsinglelep_trig[ihlt].single_other_pt_cuts)
		{
		  trig_threshold_pass = true;
		  single_trig_pass = true;
		  fired_single_trig = ihlt;
		  break;
		}
	      else if ( vsinglelep_trig[ihlt].single_other_pt_cuts<0)
		{
		  trig_threshold_pass = true;
		  single_trig_pass = true;
		  fired_single_trig = ihlt;
		  break;
		}
	    }
	}
      }
    }

    // check if offline objects match to trigger objects //

    bool lep1_match = false;
    bool lep2_match = false;
    bool jet_match = false;

    float ptratmin_lep1(0.3), ptratmin_lep2(0.3), ptratmin_jet(0.3);

    for (uint trv=0; trv<TrigRefObj.size(); trv++) {
      bool mutrobj(false), eltrobj(false), jettrobj(false);
      if (abs(TrigRefObj[trv].first)==13) mutrobj=true;
      else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() < 1.e-3) eltrobj=true;
      else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.e-3) jettrobj=true;
      //else if (abs(TrigRefObj[trv].first)==0 && TrigRefObj[trv].second.M() > 1.) jettrobj=true;
      //else if (abs(TrigRefObj[trv].first)==0) eltrobj=true;

      TVector3 Trv = TrigRefObj[trv].second.Vect();

      if (mutrobj||eltrobj) {

	TVector3 flep1v = lepcand_1.p4.Vect();
	TVector3 flep2v = lepcand_2.p4.Vect();

	double tmprat1 = fabs((flep1v-Trv).Mag()/max(1.e-6,flep1v.Mag()));
	double tmprat2 = fabs((flep2v-Trv).Mag()/max(1.e-6,flep2v.Mag()));

	hist_init[0]->Fill(tmprat1, weight);
	hist_init[1]->Fill(tmprat2, weight);

	if(abs(lepcand_1.pdgId)==11||abs(lepcand_2.pdgId)==11) {
	  TVector3 felv = (abs(lepcand_1.pdgId)==11) ? flep1v : flep2v;
	  double tmpratel = (abs(lepcand_1.pdgId)==11) ? tmprat1 : tmprat2;
	  hist2d[0]->Fill(felv.Mag(), tmpratel, weight);
	}
	if(abs(lepcand_1.pdgId)==13||abs(lepcand_2.pdgId)==13) {
	  TVector3 fmuv = (abs(lepcand_1.pdgId)==13) ? flep1v : flep2v;
	  double tmpratmu = (abs(lepcand_1.pdgId)==13) ? tmprat1 : tmprat2;
	  hist2d[0]->Fill(fmuv.Mag(), tmpratmu, weight);
	}
	if (tmprat1 < ptratmin_lep1) {
	  if((mutrobj && abs(lepcand_1.pdgId)==13)||(eltrobj && abs(lepcand_1.pdgId)==11)){
	    lep1_match = true;
	    ptratmin_lep1 = tmprat1;
	  }
	}

	if (tmprat2 < ptratmin_lep2) {
	  if((mutrobj && abs(lepcand_2.pdgId)==13)||(eltrobj && abs(lepcand_2.pdgId)==11)){
	    lep2_match = true;
	    ptratmin_lep2 = tmprat2;
	  }
	}
      }

      if(jettrobj && Jets.size()>0){

	double tmprat = fabs((Jets[0].p4.Vect()-Trv).Mag()/max(1.e-6,(Jets[0].p4.Vect().Mag())));
	if (tmprat < ptratmin_jet) {
	  jet_match = true;
	  ptratmin_jet = tmprat;
	}

      }
    }

    //if (((double_hlts.size()>0 && double_trig_pass && lep1_match && lep2_match) || (single_hlts.size()>0 && single_trig_pass && lep1_match))) { trig_matching_pass = true; }
    if (vdoublelep_trig.size()>0 && double_trig_pass && lep1_match && lep2_match)  { trig_matching_pass = true; }
    else if (vsinglelep_trig.size()>0 && single_trig_pass && (lep1_match||lep2_match) && vsinglelep_trig[fired_single_trig].single_other_pt_cuts>0 && jet_match) { trig_matching_pass = true; }
    else if (vsinglelep_trig.size()>0 && single_trig_pass && (lep1_match||lep2_match) && vsinglelep_trig[fired_single_trig].single_other_pt_cuts<0) { trig_matching_pass = true; }

  }

}
