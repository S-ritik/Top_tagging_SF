#include "plotscripts.h"
#include "filesystem.h"
#include "leptop_sf_functions.C"

vector<fileinfo> eltopfiles{jethtdata,ttbar_signal,ttbar,qcd,stop,dyjets,wjets,dib,ttx};
vector<fileinfo> mutopfiles{jethtdata,ttbar_signal,ttbar,qcd,stop,dyjets,wjets,dib,ttx};

void makehistroot()
{

  string eltop_tagger_score = "0.8";
  string mutop_tagger_score = "0.5";

  string dataname_eltop[9] = {"Data","t #bar{t} (el+jets)","t #bar{t} (other)","QCD","Single t","Z+j","W+j","DiBoson","ttX"};
  string dataname_mutop[9] = {"Data","t #bar{t} (#mu+jets)","t #bar{t} (other)","QCD","Single t","Z+j","W+j","DiBoson","ttX"};
  int col_all[8]  = {kRed,kMagenta-7,kGreen+4,kBlue,kOrange+7,kCyan+1,kRed-5,kViolet+6};

  string allsys_names[] = {"JES_AbsoluteStat",
			       "JES_AbsoluteScale",
			       "JES_AbsoluteMPFBias",
			       "JES_FlavorQCD",
			       "JES_Fragmentation",
			       "JES_PileUpDataMC",
			       "JES_PileUpPtBB",
			       "JES_PileUpPtEC1",
			       "JES_PileUpPtEC2",
			       "JES_PileUpPtRef",
			       "JES_RelativeFSR",
			       "JES_RelativeJEREC1",
			       "JES_RelativeJEREC2",
			       "JES_RelativePtBB",
			       "JES_RelativePtEC1",
			       "JES_RelativePtEC2",
			       "JES_RelativeBal",
			       "JES_RelativeSample",
			       "JES_RelativeStatEC",
			       "JES_RelativeStatFSR",
			       "JES_SinglePionECAL",
			       "JES_SinglePionHCAL",
			       "JES_TimePtEta",
			       "JER",
			       "pileup",
			       "btag",
			       "hadtagger",
			       "scale",
			       "partonshower",
      };
      
  int nallsys = sizeof(allsys_names)/sizeof(allsys_names[0]); ;
      
  const int n_sf_ptbins = 5;
  int sf_ptbins[n_sf_ptbins] = {300,400,480,600,3000};
   
  const int nvar = 4;
  const int nvarplot = 2;

  TH1D *hist_basic[nvar];
  int nbins[nvar] = {30,30,30,30};
  float binlow[nvar] = {0,0,0,0};
  float binup[nvar] = {250,250,1500,1500};

  const int nbins_eltop = 8;
  float bins_eltop_fail[nbins_eltop+1] = {30,50,70,85,105,125,150,180,250};
  float bins_eltop_pass[nbins_eltop+1] = {30,55,70,85,100,115,130,150,250};
  string var[nvar] = {"probejet_sdmass","tagjet_sdmass","probejet_pt","tagjet_pt"};
  string vartitle[nvar] = {"Soft drop mass of probe jet (in GeV)","Soft drop mass of tag jet (in GeV)","p_{t} of probe jet (in GeV)","p_{t} of tag jet (in GeV)"};

  string histname_eltop[nvar] = {"hist_eltop_probejet_sdmass","hist_eltop_tagjet_sdmass","hist_eltop_probejet_pt","hist_eltop_tagjet_pt"};
  
  string weight_eltop = "total_weight_eltop/weight_toppt";
  string event_selection_cuts_eltop = "nelectrons_eltop>0 && nvetomuons_eltop==0 && chs_met > 80 && nbjets>0 && nak8jets>1 && tagjet_pt>300 && tagjet_jetidlepveto==1 && bjet_inhadjet < 0.8 && abs(probejet_leppdg) == 11 && tagjet_sdmass > 30 && tagjet_hadtagscore > 0.97 && probejet_sdmass > 30 && jet_trigger_unprescale == 1 && ";

  for(int ir = 0; ir < 3; ir++)
    {
      for(int iv=0; iv< nvarplot; iv++)
	{
	  string evt_sel,hist_name_withregion;
	  if(ir == 0)
	    evt_sel = event_selection_cuts_eltop + "probejet_toptagscore < " + eltop_tagger_score;
	  else if(ir ==1)
	    evt_sel = event_selection_cuts_eltop + "probejet_toptagscore > " + eltop_tagger_score;
	  else
	    evt_sel = event_selection_cuts_eltop + " 10>2";
	  
	  if(ir == 0)
	    hist_name_withregion = histname_eltop[iv] + "_allptbin_fail";
	  else if(ir == 1)
	    hist_name_withregion = histname_eltop[iv] + "_allptbin_pass";
	  else
	    hist_name_withregion = histname_eltop[iv] + "_allptbin";
	  
	  if(ir == 0 && iv == 0 )
	    hist_basic[iv] = new TH1D(hist_name_withregion.c_str(),hist_name_withregion.c_str(),nbins_eltop,bins_eltop_fail);
	  else if(ir == 1 && iv == 0)
	    hist_basic[iv] = new TH1D(hist_name_withregion.c_str(),hist_name_withregion.c_str(),nbins_eltop,bins_eltop_pass);
	  else
	    hist_basic[iv] = new TH1D(hist_name_withregion.c_str(),hist_name_withregion.c_str(),nbins[iv],binlow[iv],binup[iv]);

	  //makerootfile(var[iv],weight_eltop,evt_sel,hist_basic[iv],hist_name_withregion,vartitle[iv],{jethtdata,ttbar_signal,ttbar,qcd,stop,dyjets,wjets,dib,ttx});

	  evt_sel = evt_sel + " && probejet_topmatch == 0 && probejet_top_borelmatch == 0";
	  //makerootfile(var[iv],weight_eltop,evt_sel,hist_basic[iv],(hist_name_withregion + "_nomatch"),vartitle[iv],eltopfiles,"combine_eltop_allptbin.root");

	  plot_datamcratio(9,dataname_eltop,col_all,hist_name_withregion,vartitle[iv],{jethtdata,ttbar_signal,ttbar,qcd,stop,dyjets,wjets,dib,ttx});
	}
    }


  const int nbins_mutop = 8;
  float bins_mutop_fail[nbins_mutop+1] = {30,50,70,85,105,125,150,180,250};
  float bins_mutop_pass[nbins_mutop+1] = {30,55,70,85,100,115,130,150,250};

  string histname_mutop[nvar] = {"hist_mutop_probejet_sdmass","hist_mutop_tagjet_sdmass","hist_mutop_probejet_pt","hist_mutop_tagjet_pt"};
  
  string weight_mutop = "total_weight_mutop/weight_toppt";
  string event_selection_cuts_mutop = "nmuons_mutop>0 && nvetoelectrons_mutop==0 && chs_met > 20 && nbjets>0 && nak8jets>1 && tagjet_pt>300 && tagjet_jetidlepveto==1 && bjet_inhadjet < 0.8 && abs(probejet_leppdg) == 13 && tagjet_sdmass > 30 && tagjet_hadtagscore > 0.97 && probejet_sdmass > 30 && jet_trigger_unprescale == 1 && ";

  for(int ir = 0; ir < 3; ir++)
    {
      for(int iv=0; iv< nvarplot; iv++)
	{
	  string evt_sel,hist_name_withregion;
	  if(ir == 0)
	    evt_sel = event_selection_cuts_mutop + "probejet_toptagscore < " + mutop_tagger_score;
	  else if(ir == 1)
	    evt_sel = event_selection_cuts_mutop + "probejet_toptagscore > " + mutop_tagger_score;
	  else
	    evt_sel = event_selection_cuts_mutop + " 10>2";
	  
	  if(ir == 0)
	    hist_name_withregion = histname_mutop[iv] + "_allptbin_fail";
	  else if(ir==1)
	    hist_name_withregion = histname_mutop[iv] + "_allptbin_pass";
	  else
	    hist_name_withregion = histname_mutop[iv] + "_allptbin";

	  if(ir == 0 && iv == 0 )
	    hist_basic[iv] = new TH1D(hist_name_withregion.c_str(),hist_name_withregion.c_str(),nbins_mutop,bins_mutop_fail);
	  else if(ir == 1 && iv == 0)
	    hist_basic[iv] = new TH1D(hist_name_withregion.c_str(),hist_name_withregion.c_str(),nbins_mutop,bins_mutop_pass);
	  else
	    hist_basic[iv] = new TH1D(hist_name_withregion.c_str(),hist_name_withregion.c_str(),nbins[iv],binlow[iv],binup[iv]);

	  makerootfile(var[iv],weight_mutop,evt_sel,hist_basic[iv],hist_name_withregion,vartitle[iv],{jethtdata,ttbar_signal,ttbar,qcd,stop,dyjets,wjets,dib,ttx});

	  evt_sel = evt_sel + " && probejet_topmatch == 0 && probejet_top_bormumatch == 0";
	  makerootfile(var[iv],weight_mutop,evt_sel,hist_basic[iv],(hist_name_withregion + "_nomatch"),vartitle[iv],mutopfiles,"combine_mutop_allptbin.root");

	  plot_datamcratio(9,dataname_mutop,col_all,hist_name_withregion,vartitle[iv],{jethtdata,ttbar_signal,ttbar,qcd,stop,dyjets,wjets,dib,ttx});
	}
    }
  
  
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////// For eltop hisotgram ///////////////////////////////////////////////////////////////////////////

  /*
  for(int iptbins = 10; iptbins < n_sf_ptbins; iptbins++)
   {
      int ptbin = sf_ptbins[iptbins];
      int minptbin = sf_ptbins[iptbins-1];

      string el_combine_filename = "combine_eltop_maxprobept" + to_string(ptbin) + ".root";
  
      const int nvar_eltop = 8;
  
      string weight_eltop = "total_weight_eltop";
      string event_selection_cuts_eltop = "nelectrons_eltop>0 && nvetomuons_eltop==0 && chs_met > 80 && nbjets>0 && nak8jets>1 && tagjet_pt>300 && tagjet_jetidlepveto==1 && bjet_inhadjet < 0.8 && abs(probejet_leppdg) == 11 && tagjet_sdmass > 30 && tagjet_hadtagscore > 0.97 && probejet_sdmass > 30 && jet_trigger_unprescale == 1 && probejet_pt > " + to_string(minptbin) + " && probejet_pt < " + to_string(ptbin) + " && ";

      string probecut_eltop[nvar_eltop] = {"probejet_toptagscore > " + eltop_tagger_score,"probejet_toptagscore > " + eltop_tagger_score,"probejet_toptagscore > " + eltop_tagger_score,"probejet_toptagscore > " + eltop_tagger_score,"probejet_toptagscore < " + eltop_tagger_score,"probejet_toptagscore < " + eltop_tagger_score,"probejet_toptagscore < " + eltop_tagger_score,"probejet_toptagscore < " + eltop_tagger_score};
      
      string topgencut_eltop[nvar_eltop] = {" && probejet_topmatch == 1"," && probejet_topmatch == 0 && probejet_top_borelmatch == 1"," && probejet_topmatch == 0 && probejet_top_borelmatch == 0"," "," && probejet_topmatch == 1"," && probejet_topmatch == 0 && probejet_top_borelmatch == 1"," && probejet_topmatch == 0 && probejet_top_borelmatch == 0"," "};
    
      TH1D *hist_basic_eltop[nvar_eltop];
      string histname_eltop[nvar_eltop] = {"sdmass_probejet_pass_topmatch_eltop","sdmass_probejet_pass_blepmatch_eltop","sdmass_probejet_pass_nomatch_eltop","sdmass_probejet_pass_all_eltop","sdmass_probejet_fail_topmatch_eltop","sdmass_probejet_fail_blepmatch_eltop","sdmass_probejet_fail_nomatch_eltop","sdmass_probejet_fail_all_eltop"};
      const int nbins_eltop = 8;
      float bins_eltop_fail[nbins_eltop+1] = {30,50,70,85,105,125,150,180,250};
      float bins_eltop_pass[nbins_eltop+1] = {30,55,70,85,100,115,130,150,250};
      //int nbins_eltop[nvar_eltop] = {25,25,25,25,25,25};
      //  float binlow_eltop[nvar_eltop] = {0,0,0,0,0,0};
      //float binup_eltop[nvar_eltop] = {250,250,250,250,250,250};

      string var_eltop[nvar_eltop] = {"probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass"};
      string vartitle_eltop[nvar_eltop] = {"Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Softx drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)"};
  
      for(int iv=0; iv<nvar_eltop; iv++)
	{
	  string evt_sel = event_selection_cuts_eltop + probecut_eltop[iv] + topgencut_eltop[iv];
	  if(iv<nvar_eltop/2)
	    hist_basic_eltop[iv] = new TH1D(histname_eltop[iv].c_str(),histname_eltop[iv].c_str(),nbins_eltop,bins_eltop_pass);
	  else
	    hist_basic_eltop[iv] = new TH1D(histname_eltop[iv].c_str(),histname_eltop[iv].c_str(),nbins_eltop,bins_eltop_fail);

	  string give_el_combine_filename = el_combine_filename;
	  if(iv == 3 || iv == 7) give_el_combine_filename = "";
	  make_theory_syshist(var_eltop[iv],weight_eltop,evt_sel,hist_basic_eltop[iv],histname_eltop[iv],vartitle_eltop[iv],eltopfiles,give_el_combine_filename);
	  make_jes_jer_syshist(var_eltop[iv],weight_eltop,evt_sel,hist_basic_eltop[iv],histname_eltop[iv],vartitle_eltop[iv],eltopfiles,give_el_combine_filename);
	  makerootfile(var_eltop[iv],weight_eltop,evt_sel,hist_basic_eltop[iv],histname_eltop[iv],vartitle_eltop[iv],eltopfiles,give_el_combine_filename);
	}
      
      

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////// For eltop sysmatics hisotgram ///////////////////////////////////////////////////////////////////////////

      const int nweightsys_eltop = 3;

      string nomweight_eltop[nweightsys_eltop] = {"weight_pileup","weight_btag","weight_hadtagger"};
      
      string histname_sys_eltop[nweightsys_eltop] = {"pileup","btag","hadtagger"};

      for(int iv=0; iv<nvar_eltop; iv++)
      {
      for(int isys = 0; isys < nweightsys_eltop; isys++)
      {
	string evt_sel = event_selection_cuts_eltop + probecut_eltop[iv] + topgencut_eltop[iv];
	string hname = histname_eltop[iv] + "_" + histname_sys_eltop[isys];
	if(iv<nvar_eltop/2)
	hist_basic_eltop[iv] = new TH1D(hname.c_str(),hname.c_str(),nbins_eltop,bins_eltop_pass);
	else
	hist_basic_eltop[iv] = new TH1D(hname.c_str(),hname.c_str(),nbins_eltop,bins_eltop_fail);

		string give_el_combine_filename = el_combine_filename;
		if(iv == 3 || iv == 7) give_el_combine_filename = "";
		make_weightsyshist(var_eltop[iv],weight_eltop,nomweight_eltop[isys],evt_sel,hist_basic_eltop[iv],hname,vartitle_eltop[iv],eltopfiles,give_el_combine_filename);
	    }
	} 
      addhisttofile(el_combine_filename,"sdmass_probejet_pass_data_obs_eltop","combined_ToptaggerSF_Summer20UL18_JetHT_2018_output.root","sdmass_probejet_pass_nomatch_eltop");      
      addhisttofile(el_combine_filename,"sdmass_probejet_fail_data_obs_eltop","combined_ToptaggerSF_Summer20UL18_JetHT_2018_output.root","sdmass_probejet_fail_nomatch_eltop");
      
      addhisttofile(el_combine_filename,"sdmass_probejet_fail_qcd_eltop","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_eltop");
      addhisttofile(el_combine_filename,"sdmass_probejet_pass_qcd_eltop","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_eltop");

      for(int ih = 0; ih < nallsys; ih++)
	{
	  string histname = allsys_names[ih];
	  addhisttofile(el_combine_filename,"sdmass_probejet_fail_qcd_eltop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_eltop_" + histname + "Up");
	  addhisttofile(el_combine_filename,"sdmass_probejet_pass_qcd_eltop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_eltop_" + histname + "Up");
	  addhisttofile(el_combine_filename,"sdmass_probejet_fail_qcd_eltop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_eltop_" + histname + "Down");
	  addhisttofile(el_combine_filename,"sdmass_probejet_pass_qcd_eltop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_eltop_" + histname + "Down");
	}

      for(int ih = 1; ih < 103; ih++)
	{
	  string histname = "pdf" + to_string(ih);
	  addhisttofile(el_combine_filename,"sdmass_probejet_fail_qcd_eltop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_eltop_" + histname + "Up");
	  addhisttofile(el_combine_filename,"sdmass_probejet_pass_qcd_eltop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_eltop_" + histname + "Up");
	  addhisttofile(el_combine_filename,"sdmass_probejet_fail_qcd_eltop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_eltop_" + histname + "Down");
	  addhisttofile(el_combine_filename,"sdmass_probejet_pass_qcd_eltop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_eltop_" + histname + "Down");
	}

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////// For mutop hisotgram ///////////////////////////////////////////////////////////////////////////

      string mu_combine_filename = "combine_mutop_maxprobept" + to_string(ptbin) + ".root";
    
      const int nvar_mutop = 8;

      string weight_mutop = "total_weight_mutop";
      string event_selection_cuts_mutop = "nmuons_mutop>0 && nvetoelectrons_mutop==0 && chs_met > 20 && nbjets>0 && nak8jets>1 && tagjet_pt>300 && tagjet_jetidlepveto==1 && bjet_inhadjet < 0.8 && abs(probejet_leppdg) == 13 && tagjet_sdmass > 30 && tagjet_hadtagscore > 0.97 && probejet_sdmass > 30 && jet_trigger_unprescale == 1 && probejet_pt > " + to_string(minptbin) + " && probejet_pt < " + to_string(ptbin) + " && ";

      string probecut_mutop[nvar_mutop] = {"probejet_toptagscore > " + mutop_tagger_score,"probejet_toptagscore > " + mutop_tagger_score,"probejet_toptagscore > " + mutop_tagger_score,"probejet_toptagscore > " + mutop_tagger_score,"probejet_toptagscore < " + mutop_tagger_score,"probejet_toptagscore < " + mutop_tagger_score,"probejet_toptagscore < " + mutop_tagger_score,"probejet_toptagscore < " + mutop_tagger_score};
      string topgencut_mutop[nvar_mutop] = {" && probejet_topmatch == 1"," && probejet_topmatch == 0 && probejet_top_bormumatch == 1"," && probejet_topmatch == 0 && probejet_top_bormumatch == 0"," "," && probejet_topmatch == 1"," && probejet_topmatch == 0 && probejet_top_bormumatch == 1"," && probejet_topmatch == 0 && probejet_top_bormumatch == 0"," "};
    
      TH1D *hist_basic_mutop[nvar_mutop];
      string histname_mutop[nvar_mutop] = {"sdmass_probejet_pass_topmatch_mutop","sdmass_probejet_pass_blepmatch_mutop","sdmass_probejet_pass_nomatch_mutop","sdmass_probejet_pass_all_mutop","sdmass_probejet_fail_topmatch_mutop","sdmass_probejet_fail_blepmatch_mutop","sdmass_probejet_fail_nomatch_mutop","sdmass_probejet_fail_all_mutop"};
      const int nbins_mutop = 8;
      float bins_mutop_fail[nbins_mutop+1] = {30,50,70,85,105,125,150,180,250};
      float bins_mutop_pass[nbins_mutop+1] = {30,55,70,85,100,115,130,150,250};

      //  int nbins_mutop[nvar_mutop] = {25,25,25,25,25,25};
      //float binlow_mutop[nvar_mutop] = {0,0,0,0,0,0};
      //float binup_mutop[nvar_mutop] = {250,250,250,250,250,250};

      string var_mutop[nvar_mutop] = {"probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass","probejet_sdmass"};
      string vartitle_mutop[nvar_mutop] = {"Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)","Soft drop mass of probe jet (GeV)"};
  
      for(int iv=0; iv<nvar_mutop; iv++)
	{
	  string evt_sel = event_selection_cuts_mutop + probecut_mutop[iv] + topgencut_mutop[iv];
	  if(iv<nvar_mutop/2)
	    hist_basic_mutop[iv] = new TH1D(histname_mutop[iv].c_str(),histname_mutop[iv].c_str(),nbins_mutop,bins_mutop_pass);
	  else
	    hist_basic_mutop[iv] = new TH1D(histname_mutop[iv].c_str(),histname_mutop[iv].c_str(),nbins_mutop,bins_mutop_fail);

	  string give_mu_combine_filename = mu_combine_filename;
	  if(iv == 3 || iv == 7) give_mu_combine_filename = "";

	  make_theory_syshist(var_mutop[iv],weight_mutop,evt_sel,hist_basic_mutop[iv],histname_mutop[iv],vartitle_mutop[iv],mutopfiles,give_mu_combine_filename);
	  make_jes_jer_syshist(var_mutop[iv],weight_mutop,evt_sel,hist_basic_mutop[iv],histname_mutop[iv],vartitle_mutop[iv],mutopfiles,give_mu_combine_filename);
	  makerootfile(var_mutop[iv],weight_mutop,evt_sel,hist_basic_mutop[iv],histname_mutop[iv],vartitle_mutop[iv],mutopfiles,give_mu_combine_filename);
	}

      

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////// For mutop sysmatics hisotgram ///////////////////////////////////////////////////////////////////////////

      const int nweightsys_mutop = 3;

      string nomweight_mutop[nweightsys_mutop] = {"weight_pileup","weight_btag","weight_hadtagger"};
      
      string histname_sys_mutop[nweightsys_mutop] = {"pileup","btag","hadtagger"};

      for(int iv=0; iv<nvar_mutop; iv++)
	{
	  for(int isys = 0; isys < nweightsys_mutop; isys++)
	    {
	      string evt_sel = event_selection_cuts_mutop + probecut_mutop[iv] + topgencut_mutop[iv];
	      string hname = histname_mutop[iv] + "_" + histname_sys_mutop[isys];
	      if(iv<nvar_mutop/2)
		hist_basic_mutop[iv] = new TH1D(hname.c_str(),hname.c_str(),nbins_mutop,bins_mutop_pass);
	      else
		hist_basic_mutop[iv] = new TH1D(hname.c_str(),hname.c_str(),nbins_mutop,bins_mutop_fail);

		string give_mu_combine_filename = mu_combine_filename;
		if(iv == 3 || iv == 7) give_mu_combine_filename = "";

	      make_weightsyshist(var_mutop[iv],weight_mutop,nomweight_mutop[isys],evt_sel,hist_basic_mutop[iv],hname,vartitle_mutop[iv],mutopfiles,give_mu_combine_filename);
	    }
	}
      
      addhisttofile(mu_combine_filename,"sdmass_probejet_pass_data_obs_mutop","combined_ToptaggerSF_Summer20UL18_JetHT_2018_output.root","sdmass_probejet_pass_nomatch_mutop");
      addhisttofile(mu_combine_filename,"sdmass_probejet_fail_data_obs_mutop","combined_ToptaggerSF_Summer20UL18_JetHT_2018_output.root","sdmass_probejet_fail_nomatch_mutop");

      addhisttofile(mu_combine_filename,"sdmass_probejet_fail_qcd_mutop","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_mutop");
      addhisttofile(mu_combine_filename,"sdmass_probejet_pass_qcd_mutop","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_mutop");

      for(int ih = 0; ih < nallsys; ih++)
	{
	  string histname = allsys_names[ih];
	  addhisttofile(mu_combine_filename,"sdmass_probejet_fail_qcd_mutop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_mutop_" + histname + "Up");
	  addhisttofile(mu_combine_filename,"sdmass_probejet_pass_qcd_mutop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_mutop_" + histname + "Up");
	  addhisttofile(mu_combine_filename,"sdmass_probejet_fail_qcd_mutop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_mutop_" + histname + "Down");
	  addhisttofile(mu_combine_filename,"sdmass_probejet_pass_qcd_mutop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_mutop_" + histname + "Down");
	}

      for(int ih = 1; ih < 103; ih++)
	{
	  string histname = "pdf" + to_string(ih);
	  addhisttofile(mu_combine_filename,"sdmass_probejet_fail_qcd_mutop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_mutop_" + histname + "Up");
	  addhisttofile(mu_combine_filename,"sdmass_probejet_pass_qcd_mutop_" + histname + "Up","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_mutop_" + histname + "Up");
	  addhisttofile(mu_combine_filename,"sdmass_probejet_fail_qcd_mutop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_fail_all_mutop_" + histname + "Down");
	  addhisttofile(mu_combine_filename,"sdmass_probejet_pass_qcd_mutop_" + histname + "Down","combined_ToptaggerSF_Summer20UL18_QCD_output.root","sdmass_probejet_pass_all_mutop_" + histname + "Down");
	}

      

   }
    
//////////////////////////////////      For data-driven QCD plots           ///////////////////////////////
  
   // string datadriven_event_selection_cuts_eltop = "nelectrons_eltop>0 && nvetomuons_eltop==0 && nak8jets>1 && tagjet_pt>300 && tagjet_jetidlepveto==1 && abs(probejet_leppdg) == 11 && tagjet_sdmass > 30 && tagjet_hadtagscore > 0.97 && probejet_sdmass > 30 && probejet_toptagscore < " + eltop_tagger_score;
  /*  string datadriven_weight_eltop1 = "total_weight_eltop";
  float metbins[3] = {0,80,800};
  float nbjetbins[3] = {0,1,10};
  float hadtopbins[3] = {0,0.97,1};

  string datadriven_event_selection_cuts_eltop2 = "nelectrons_eltop>0 && nvetomuons_eltop==0 && nbjets>0 && nak8jets>1 && tagjet_pt>300 && tagjet_jetidlepveto==1 && jet_trigger_unprescale == 1 && bjet_inhadjet < 0.8 && abs(probejet_leppdg) == 11 && tagjet_sdmass > 30 && probejet_sdmass > 30 && probejet_toptagscore > " + eltop_tagger_score;
  string datadriven_weight_eltop2 = "(tagjet_hadtagscore>0.97) ? total_weight_eltop : total_weight_eltop/weight_hadtagger";
  
  TH2D *hist2d2 = new TH2D("hist2d2","hadtagscore score vs met",2,metbins,2,hadtopbins);
  cout<<"check 1"<<endl;
  makerootfile("tagjet_hadtagscore:chs_met",datadriven_weight_eltop2,datadriven_event_selection_cuts_eltop2,hist2d2,"eltophist2d_hadtopscorevsmet","MET","Hadronic top tagger score of tag jet",{jethtdata,ttbar_signal,qcd});
  
  cout<<"check 2"<<endl;
//plotrootfile("eltophist2d_hadtopscorevsmet",{jethtdata,ttbar_signal,qcd});
  
  string sr_cr_cut[2]; string transfer_factor_cut[2];
  transfer_factor_cut[0] = "tagjet_hadtagscore<0.8"; transfer_factor_cut[1] = "tagjet_hadtagscore>0.97"; sr_cr_cut[0] = "chs_met<50"; sr_cr_cut[1] = "chs_met>80";
  
  compute_QCD_ABCD("eltophist2d_hadtopscorevsmet","sdmass_probejet_fail_topmatch_eltop","probejet_sdmass","Soft drop mass of probe jet (GeV)",datadriven_event_selection_cuts_eltop2, sr_cr_cut, transfer_factor_cut,datadriven_weight_eltop2,jethtdata,qcd,{ttbar_signal,ttbar,stop,dyjets,wjets,dib,ttx});
  */
  // TH2D *hist2d = new TH2D("hist2d","nbjets vs met",6,metbins,2,nbjetbins);
  // makerootfile("nbjets:chs_met",datadriven_weight_eltop1,datadriven_event_selection_cuts_eltop,hist2d,"eltophist2d_nbjetvsmet","MET","No. of bjets",{egammadata,ttbar_signal,qcd});
  // plotrootfile("eltophist2d_nbjetvsmet",{egammadata,ttbar_signal,qcd});
  

  ///////////////////////       For cutflow histograms   //////////////////////////////
/*
  string dataname_eltop[9] = {"Data","t #bar{t} (el+jets)","t #bar{t} (other)","QCD","Single t","Z+j","W+j","DiBoson","ttX"};
  string dataname_mutop[9] = {"Data","t #bar{t} (#mu+jets)","t #bar{t} (other)","QCD","Single t","Z+j","W+j","DiBoson","ttX"};

  const int nhists = 12;
  const int nchannels = 2;

  string prvar_name[nhists][nchannels] = {{"nAK8jets","nAK8jets"},{"leadingAK8pt","leadingAK8pt"},{"nelectrons","nmuons"},{"nmuons_veto","nelectrons_veto"},{"lep_pdg_AK8matched","lep_pdg_AK8matched"},{"tagjet_jetidlepveto","tagjet_jetidlepveto"},{"tagjet_pnettagscore","tagjet_pnettagscore"},{"tagjet_sdmass","tagjet_sdmass"},{"nbjets","nbjets"},{"bjet_inhadjet","bjet_inhadjet"},{"chs_met","chs_met"},{"jet_trigger_unprescale","jet_trigger_unprescale"}};

  string prvar_title[nhists][nchannels] = {{"No. of AK8 jets","No. of AK8 jets"},{"p_{T} of leading AK8 (in GeV)","p_{T} of leading AK8 (in GeV)"},{"No. of electrons","No. of muons"},{"No. of muons for vetoing","No. of electrons for vetoing"},{"pdg id of leption matched to AK8","pdg id of leption matched to AK8"},{"lepveto of tag jet","lepveto of tag jet"},{"PNet hadronic top tagger score of tagjet","PNet hadronic top tagger score of tagjet"},{"SD mass of tagjet (in GeV)","SD mass of tagjet (in GeV)"},{"No. of b-jets","No. of b-jets"},{"bjet_inhadjet","bjet_inhadjet"},{"chs_met","chs_met"},{"jet trigger fired (unprescale)","jet triggers fired (unprescale)"}};

  string channelname[nchannels] = {"eltop","mutop"};
  int col[8]  = {kRed,kMagenta-7,kGreen+4,kBlue,kOrange+7,kCyan+1,kRed-5,kViolet+6};

  for(int ih = 0; ih < nhists+1; ih++)
    {
      string histname;
      if(ih<nhists)histname = "hist_cutvar_" + prvar_name[ih][0] + "_" + channelname[0];
      if(ih<nhists)    compare_datamc_fromallfiles(9,dataname_eltop,col,histname,prvar_title[ih][0],eltopfiles);
      histname = "hist_cutflow_ak8pt_after_" + to_string(ih) + "_cut_" + channelname[0];
      compare_datamc_fromallfiles(9,dataname_eltop,col,histname,"p_{T} of probe AK8 jet (in GeV)",eltopfiles);
      histname = "hist_cutflow_leppt_after_" + to_string(ih) + "_cut_" + channelname[0];
      compare_datamc_fromallfiles(9,dataname_eltop,col,histname,"p_{T} of electron inside probe AK8 jet (in GeV)",eltopfiles);

      if(ih<nhists)histname = "hist_cutvar_" + prvar_name[ih][1] + "_" + channelname[1];
      if(ih<nhists)    compare_datamc_fromallfiles(9,dataname_mutop,col,histname,prvar_title[ih][1],mutopfiles);
      histname = "hist_cutflow_ak8pt_after_" + to_string(ih) + "_cut_" + channelname[1];
      compare_datamc_fromallfiles(9,dataname_mutop,col,histname,"p_{T} of probe AK8 jet (in GeV)",mutopfiles);
      histname = "hist_cutflow_leppt_after_" + to_string(ih) + "_cut_" + channelname[1];
      compare_datamc_fromallfiles(9,dataname_mutop,col,histname,"p_{T} of muon inside probe AK8 jet (in GeV)",mutopfiles);
    }
  */

 
/*  const int nprocs = 3;
  string processes[nprocs] = {"blepmatch","nomatch","topmatch"};
  string eltopchannelname[nprocs] = {"top semi-merge jets","no matched jets","top merge jets"};
  string mutopchannelname[nprocs] = {"top semi-merge jets","no matched jets","top merge jets"};
  
  TGraphAsymmErrors *gr_eltop[nprocs];
  TGraphAsymmErrors *gr_mutop[nprocs];
  
  for(int ip = 0; ip < nprocs; ip++)
    {
      gr_eltop[ip] = new TGraphAsymmErrors(n_sf_ptbins-1);
      gr_mutop[ip] = new TGraphAsymmErrors(n_sf_ptbins-1);
      
     for(int iptbins = 1; iptbins < n_sf_ptbins; iptbins++)
       {
	 string processes_inorder[nprocs];
	 string eltopchannelname_inorder[nprocs];
	 string mutopchannelname_inorder[nprocs];

	 for(int ij =0; ij < nprocs; ij++)
	   {
	     processes_inorder[ij] = processes[int((ip+ij)%nprocs)];
	     eltopchannelname_inorder[ij] = eltopchannelname[int((ip+ij)%nprocs)];
	     mutopchannelname_inorder[ij] = mutopchannelname[int((ip+ij)%nprocs)];
	   }
	 
	 int ptbin = sf_ptbins[iptbins];
	 float ptbin_mid = (sf_ptbins[iptbins] + sf_ptbins[iptbins-1])/2.0;

	 plotfromhiggscombine("eltop","fail",nprocs,processes_inorder,eltopchannelname_inorder,"_maxprobept" + to_string(ptbin));
	 plotfromhiggscombine("eltop","pass",nprocs,processes_inorder,eltopchannelname_inorder,"_maxprobept" + to_string(ptbin));
	 float *sf_witherror_eltop = calculate_sf_fromhiggscombine("eltop",nprocs,processes_inorder,eltopchannelname_inorder,"_maxprobept" + to_string(ptbin),iptbins,gr_eltop[ip]);
      
	 cout<<"eltop SF for "<<processes[ip]<<"=="<<processes_inorder[0]<<" in pt bin "<<ptbin<<" = "<<sf_witherror_eltop[0]<<" "<<sf_witherror_eltop[1]<<"(up) "<<sf_witherror_eltop[2]<<"(down) "<<endl;
	 //gr_eltop[ip]->SetPoint(iptbins-1,ptbin_mid,sf_witherror_eltop[0]);
	 //gr_eltop[ip]->SetPointError(iptbins-1,ptbin_mid - sf_ptbins[iptbins-1],sf_ptbins[iptbins] - ptbin_mid,sf_witherror_eltop[2],sf_witherror_eltop[1]);
	 //gr_eltop[ip]->SetPoint(iptbins-1,iptbins,sf_witherror_eltop[0]);
	 //gr_eltop[ip]->SetPointError(iptbins-1,0,0,sf_witherror_eltop[2],sf_witherror_eltop[1]);
	 
	 
	 plotfromhiggscombine("mutop","fail",nprocs,processes_inorder,mutopchannelname_inorder,"_maxprobept" + to_string(ptbin));
	 plotfromhiggscombine("mutop","pass",nprocs,processes_inorder,mutopchannelname_inorder,"_maxprobept" + to_string(ptbin));
	 float *sf_witherror_mutop = calculate_sf_fromhiggscombine("mutop",nprocs,processes_inorder,mutopchannelname_inorder,"_maxprobept" + to_string(ptbin),iptbins,gr_mutop[ip]);
	  
	  cout<<"mutop SF for "<<processes[ip]<<"=="<<processes_inorder[0]<<" in pt bin "<<ptbin<<" = "<<sf_witherror_mutop[0]<<" "<<sf_witherror_mutop[1]<<"(up) "<<sf_witherror_mutop[2]<<"(down) "<<endl;
	  //gr_mutop[ip]->SetPoint(iptbins-1,ptbin_mid,sf_witherror_mutop[0]);
	  //gr_mutop[ip]->SetPointError(iptbins-1,ptbin_mid - sf_ptbins[iptbins-1],sf_ptbins[iptbins] - ptbin_mid,sf_witherror_mutop[2],sf_witherror_mutop[1]);
	  //	  gr_mutop[ip]->SetPoint(iptbins-1,iptbins,sf_witherror_mutop[0]);
	  //gr_mutop[ip]->SetPointError(iptbins-1,0,0,sf_witherror_mutop[2],sf_witherror_mutop[1]);
	}
	}
 eltop_tagger_score[1] = 'p'; mutop_tagger_score[1]= 'p';  
    string file_name = "toptagging_sf_eltop" + eltop_tagger_score + "_mutop" + mutop_tagger_score + ".root";
    TFile *fout = new TFile(file_name.c_str(),"update");
    for(int ip =0; ip < nprocs; ip++)
    {
     TPaveText *text1 = new TPaveText(0.2,0.6,0.95,0.99,"NDCNB");
     text1->AddText("CMS Preliminary 59.81                                                     fb^{-1} (13 TeV)");
     //text1->SetBorder(0);
     text1->SetFillStyle(0);
     text1->SetTextSize(0.1);
     
     TPaveText *text2 = new TPaveText(0.4,0.6,0.8,0.9,"NDCNB");
     text2->AddText("Electron top tagger (WP > 0.888)");
     text2->AddText("At mistag rate = 0.24% and signal efficiency = 85.4%");
     text2->SetBorderSize(0);
     text2->SetFillStyle(0);
     text2->SetTextSize(0.035);

     TPaveText *text3 = new TPaveText(0.4,0.6,0.8,0.9,"NDCNB");
     text3->AddText("Muon top tagger (WP > 0.9)");
     text3->AddText("At mistag rate = 0.17% and signal efficiency = 91.1%");
     text3->SetBorderSize(0);
     text3->SetFillStyle(0);
     text3->SetTextSize(0.035);
     
     TAxis *ax1 = gr_eltop[ip]->GetHistogram()->GetXaxis();
     Double_t x11 = ax1->GetBinLowEdge(1); 
     Double_t x21 = ax1->GetBinUpEdge(ax1->GetNbins());
     gr_eltop[ip]->GetHistogram()->GetXaxis()->Set(n_sf_ptbins-1,x11,x21);
     
     for(Int_t k=0;k<n_sf_ptbins-1;k++){
       gr_eltop[ip]->GetHistogram()->GetXaxis()->SetBinLabel(k+1,(to_string(sf_ptbins[k])+ "-" + sf_ptbins[k+1]).Data());   
       }
     
     gr_eltop[ip]->SetMarkerColor(2);
     gr_eltop[ip]->SetMarkerStyle(20);
     gr_eltop[ip]->SetMaximum(1.8);
     gr_eltop[ip]->SetMinimum(0.0);
     gr_eltop[ip]->SetTitle(";probe jet p_{T} (in GeV);SF");

     TAxis *ax = gr_mutop[ip]->GetHistogram()->GetXaxis();
     Double_t x1 = ax->GetBinLowEdge(1); 
     Double_t x2 = ax->GetBinUpEdge(ax->GetNbins());
     gr_mutop[ip]->GetHistogram()->GetXaxis()->Set(n_sf_ptbins-1,x1,x2);
     
     for(Int_t k=0;k<n_sf_ptbins-1;k++){
       gr_mutop[ip]->GetHistogram()->GetXaxis()->SetBinLabel(k+1,(to_string(sf_ptbins[k])+ "-" + sf_ptbins[k+1]).Data());   
     }
     
     gr_mutop[ip]->SetMarkerColor(2);
     gr_mutop[ip]->SetMarkerStyle(20);
     gr_mutop[ip]->SetMaximum(1.8);
     gr_mutop[ip]->SetMinimum(0.0);
     gr_mutop[ip]->SetTitle(";probe jet p_{T} (in GeV);SF");
     
     TCanvas *c1 = tdrCanvas("canv_d", (TH1D*)gr_eltop[ip]->GetHistogram(),8,0);
     c1->cd();
     gr_eltop[ip]->Draw("AP");
     //text1->Draw("same");
     text2->Draw("same");
     CMS_lumi( c1, 8, 0 );

     c1->SaveAs(("plots/sf_eltoptagger_" + processes[ip] + "_" + eltop_tagger_score + ".png").c_str());
     
     TCanvas *c2 = tdrCanvas("canv_d", (TH1D*)gr_mutop[ip]->GetHistogram(),8,0);
     c2->cd();
     gr_mutop[ip]->Draw("AP");
     //text1->Draw("same");
     CMS_lumi( c2, 8, 0 );
     text3->Draw("same");
     c2->SaveAs(("plots/sf_mutoptagger_" + processes[ip] + "_" + mutop_tagger_score + ".png").c_str());
    
     gr_eltop[ip]->Write(("toptagger_sf_" + processes[ip] + "_eltop").c_str(),TObject::kOverwrite);
     gr_mutop[ip]->Write(("toptagger_sf_" + processes[ip] + "_mutop").c_str(),TObject::kOverwrite);     
   }
 fout->Close();
*/

  
//////////////////////////// For comparing shapes ////////////////////////////////
/*
 TH1D *hist_compare_eltop[3];
 TH1D *hist_compare_mutop[3];
 
 TFile *fin = new TFile("combined_ToptaggerSF_Summer20UL18_QCD_output.root","read");
 if(!fin->IsOpen())
   {cout<<"Check: file combined_ToptaggerSF_Summer20UL18_QCD_output.root not present"<<endl;exit(0);}
 fin->cd();
 hist_compare_eltop[0] = (TH1D*)fin->Get("hist_eltop_probejet_sdmass_allptbin_fail");
 hist_compare_eltop[0]->SetDirectory(0);
 hist_compare_mutop[0] = (TH1D*)fin->Get("hist_mutop_probejet_sdmass_allptbin_fail");
 hist_compare_mutop[0]->SetDirectory(0);
 fin->Close();
 
 fin = new TFile("combined_Summer20UL18_TTBar_output.root","read");
 if(!fin->IsOpen())
   {cout<<"Check: file combined_Summer20UL18_TTBar_output.root not present"<<endl;exit(0);}
 fin->cd();
 hist_compare_eltop[1] = (TH1D*)fin->Get("hist_eltop_probejet_sdmass_allptbin_fail");
 hist_compare_eltop[1]->SetDirectory(0);
 hist_compare_mutop[1] = (TH1D*)fin->Get("hist_mutop_probejet_sdmass_allptbin_fail");
 hist_compare_mutop[1]->SetDirectory(0);
 fin->Close();

 fin = new TFile("combine_eltop_allptbin.root","read");
 if(!fin->IsOpen())
   {cout<<"Check: file combine_eltop_allptbin.root not present"<<endl;exit(0);}
 fin->cd();
 hist_compare_eltop[2] = (TH1D*)fin->Get("hist_eltop_probejet_sdmass_nomatch");
 hist_compare_eltop[2]->SetDirectory(0);
 fin->Close();

 fin = new TFile("combine_mutop_allptbin.root","read");
 if(!fin->IsOpen())
   {cout<<"Check: file combine_mutop_allptbin.root not present"<<endl;exit(0);}
 fin->cd();
 hist_compare_mutop[2] = (TH1D*)fin->Get("hist_mutop_probejet_sdmass_nomatch");
 hist_compare_mutop[2]->SetDirectory(0);
 fin->Close();

 string legend_names_eltop[3] = {"jets from QCD process","jets from ttbar (non el+jets processes)","nomatch jets"};
 string legend_names_mutop[3] = {"jets from QCD process","jets from ttbar (non #mu+jets processes)","nomatch jets"};

 compareshape(3,hist_compare_eltop,legend_names_eltop,"plot_compare_bkgshapes_eltop.png");
 compareshape(3,hist_compare_mutop,legend_names_mutop,"plot_compare_bkgshapes_mutop.png");
*/
 /////////////// for sysmetaic check plots             ///////////////////////////

 // check_sysematic();
 
}
