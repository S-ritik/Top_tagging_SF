#include</home/ritik/Desktop/new/new_structure/btop_measurement/event_selection_best_yeield/My_Style.C>
/*void savehist(TEfficiency *hist, string plotname)
{
  TH1D *histtemp = new TH1D("tmp","tmp",60,5,10);
  //cout<<"as"<<endl;
  TCanvas *c1 = tdrCanvas("canv_d", histtemp,8,0);
  c1->cd();
  gStyle->SetOptTitle(1);
  
  hist->Draw("histtext, e");

  c1->SaveAs(("./plots/" + plotname + ".png").c_str());  
  }*/
void savehist(TH1D *hist, string plotname)
{
  TCanvas *c1 = tdrCanvas("canv_d", hist,8,0);
  c1->cd();
  gStyle->SetOptTitle(1);

  hist->Draw("histtext, e");

  c1->SaveAs(("./plots/" + plotname + ".png").c_str());  
}
void savehist(TH2D *hist, string plotname)
{
  TCanvas *c1 = tdrCanvas("canv_d", hist->ProjectionX(),8,0);
  c1->cd();
  gStyle->SetOptTitle(1);

  if( plotname.find("lepptvsleppt") != string::npos){
    gPad->SetLogx();
    gPad->SetLogy();
    //hist->Rebin2D(2,2);
  }
  gStyle->SetPaintTextFormat( "1.3f" );
  hist->Draw("colztext, e");

  c1->SaveAs(("./plots/" + plotname + ".png").c_str());
  
}

void savehist(TEfficiency *hist, string plotname)
{
  TH2D *h = new TH2D();
  h = (TH2D*)hist->CreateHistogram();

  int nxbin = h->GetNbinsX();
  int nybin = h->GetNbinsY();

  for(int ij=1; ij<=nxbin; ij++)
    {
      for(int ik=1; ik<=nybin; ik++)
	{
	  h->SetBinError(ij,ik,(hist->GetEfficiencyErrorLow(hist->GetGlobalBin(ij,ik)) + hist->GetEfficiencyErrorUp(hist->GetGlobalBin(ij,ik)))/2);
	  //h->SetBinErrorUp(ij,ik,hist->GetEfficiencyErrorUp(GetEfficiencyErrorUp(ij,ik)));
	}
    }
  
      
  TCanvas *c1 = tdrCanvas("canv_d", h->ProjectionX(),8,0);
  c1->cd();
  gStyle->SetOptTitle(1);

  if( plotname.find("lepptvsleppt") != string::npos){
    gPad->SetLogx();
    gPad->SetLogy();
    //hist->Rebin2D(2,2);
  }

  h->GetXaxis()->SetNdivisions(15);
  h->GetYaxis()->SetNdivisions(15);
  h->GetXaxis()->SetMoreLogLabels();
  h->GetYaxis()->SetMoreLogLabels();

  gStyle->SetPaintTextFormat( "1.3f" );
  h->Draw("colztext, e");
  
  c1->SaveAs(("./plots/" + plotname + ".png").c_str());
  
}

void scale_factor(TH1D* hist_orthomc, TH1D* hist_trigmc, TH1D* hist_orthodata, TH1D* hist_trigdata)
{
  TH1D *hist_sf;
  hist_sf = (TH1D*)hist_orthomc->Clone();

  int nxbin = hist_orthomc->GetNbinsX();

  for(int ij=1; ij<=nxbin; ij++)
    {
      if(hist_trigmc->GetBinContent(ij)!=0 && hist_orthodata->GetBinContent(ij)!=0 && hist_orthomc->GetBinContent(ij)!=0)
	{
	  hist_sf->SetBinContent(ij, hist_trigdata->GetBinContent(ij) * hist_orthomc->GetBinContent(ij) / max(0.0001,( hist_trigmc->GetBinContent(ij) * hist_orthodata->GetBinContent(ij) )) );
	}
      else
	{
	  hist_sf->SetBinContent(ij,0.0);
	}
    }
  string name = "sf_" + string(hist_orthomc->GetName());
  //  savehist(hist_sf,name.c_str());  
}

void scale_factor(TH2D* hist_orthomc, TH2D* hist_trigmc, TH2D* hist_orthodata, TH2D* hist_trigdata)
{
  TH2D *hist_sf;
  hist_sf = (TH2D*)hist_orthomc->Clone();

  int nxbin = hist_orthomc->GetNbinsX();
  int nybin = hist_orthomc->GetNbinsY();

  for(int ij=1; ij<=nxbin; ij++)
    {
      for(int ik=1; ik<=nybin; ik++)
	{
	  if(hist_trigmc->GetBinContent(ij,ik)!=0 && hist_orthodata->GetBinContent(ij,ik)!=0 && hist_orthomc->GetBinContent(ij,ik)!=0)
	    {
	      hist_sf->SetBinContent(ij,ik, hist_trigdata->GetBinContent(ij,ik) * hist_orthomc->GetBinContent(ij,ik) / max(0.0001,( hist_trigmc->GetBinContent(ij,ik) * hist_orthodata->GetBinContent(ij,ik) )) );
	      //hist_sf->SetBinContent(ij,ik,);
	    }
	  else
	    {
	      hist_sf->SetBinContent(ij,ik,0.0);
	    }
	}
    }
  string name = "sf_" + string(hist_trigmc->GetName());
  savehist(hist_sf,name.c_str());
}

void trigg_eff_calculator()
{
  string mcfile = "Trig_Eff_Summer20UL18_TTBar_DiLeptonic_signal_output.root";
  
  string datafile = "mergerd_Trig_Eff_Summer20UL18_MET_output.root"; //"mergerd_EE_Summer20UL18_EGamma_output.root"; 
  
  float data_lumi = 59.83;
  float mc_weight = 0.0000083923462258680;
  
  static const int nhist = 12; static const int nhist2d = 9;
  
  string historthonames[nhist] = {"ee_hist_orthotrig_dileppt","ee_hist_orthotrig_dilepmass","ee_hist_orthotrig_dilepphi","ee_hist_orthotrig_dilepeta","emu_hist_orthotrig_dileppt","emu_hist_orthotrig_dilepmass","emu_hist_orthotrig_dilepphi","emu_hist_orthotrig_dilepeta","mumu_hist_orthotrig_dileppt","mumu_hist_orthotrig_dilepmass","mumu_hist_orthotrig_dilepphi","mumu_hist_orthotrig_dilepeta"};
  string historthotitles[nhist] = {"pt of dilepton","mass of dilepton","phi of dilepton","eta of dilepton","pt of dilepton","mass of dilepton","phi of dilepton","eta of dilepton","pt of dilepton","mass of dilepton","phi of dilepton","eta of dilepton"};
  string histtrignames[nhist] = {"ee_hist_trig_eff_dileppt","ee_hist_trig_eff_dilepmass","ee_hist_trig_eff_dilepphi","ee_hist_trig_eff_dilepeta","emu_hist_trig_eff_dileppt","emu_hist_trig_eff_dilepmass","emu_hist_trig_eff_dilepphi","emu_hist_trig_eff_dilepeta","mumu_hist_trig_eff_dileppt","mumu_hist_trig_eff_dilepmass","mumu_hist_trig_eff_dilepphi","mumu_hist_trig_eff_dilepeta"};
  string histtrigtitles[nhist] = {"pt of dilepton","mass of dilepton","phi of dilepton","eta of dilepton","pt of dilepton","mass of dilepton","phi of dilepton","eta of dilepton","pt of dilepton","mass of dilepton","phi of dilepton","eta of dilepton"};

  TH1D *hist1d_orthomc[nhist];
  TH1D *hist1d_trigmc[nhist];
  TH1D *hist1d_orthodata[nhist];
  TH1D *hist1d_trigdata[nhist];
    
  string historthonames2d[nhist2d] = {"ee_hist_orthotrig_lepptvsleppt","ee_hist_orthotrig_lepetavslepeta","ee_hist_orthotrig_lepphivslepphi","emu_hist_orthotrig_lepptvsleppt","emu_hist_orthotrig_lepetavslepeta","emu_hist_orthotrig_lepphivslepphi","mumu_hist_orthotrig_lepptvsleppt","mumu_hist_orthotrig_lepetavslepeta","mumu_hist_orthotrig_lepphivslepphi"};
  string historthotitles2d[nhist2d] = {"pt of subleading vs pt of leading lepton","eta of subleading vs eta of leading lepton","phi of subleading vs phi of leading lepton","pt of subleading vs pt of leading lepton","eta of subleading vs eta of leading lepton","phi of subleading vs phi of leading lepton","pt of subleading vs pt of leading lepton","eta of subleading vs eta of leading lepton","phi of subleading vs phi of leading lepton"};
  string histtrignames2d[nhist2d] = {"ee_hist_trig_eff_lepptvsleppt","ee_hist_trig_eff_lepetavslepeta","ee_hist_trig_eff_lepphivslepphi","emu_hist_trig_eff_lepptvsleppt","emu_hist_trig_eff_lepetavslepeta","emu_hist_trig_eff_lepphivslepphi","mumu_hist_trig_eff_lepptvsleppt","mumu_hist_trig_eff_lepetavslepeta","mumu_hist_trig_eff_lepphivslepphi"};
  string histtrigtitles2d[nhist2d] = {"pt of subleading vs pt of leading lepton","eta of subleading vs eta of leading lepton","phi of subleading vs phi of leading lepton","pt of subleading vs pt of leading lepton","eta of subleading vs eta of leading lepton","phi of subleading vs phi of leading lepton","pt of subleading vs pt of leading lepton","eta of subleading vs eta of leading lepton","phi of subleading vs phi of leading lepton"};

  string histtrignames2d_unmatch[nhist2d] = {"ee_hist_trig_eff_lepptvsleppt_unmatched","ee_hist_trig_eff_lepetavslepeta_unmatched","ee_hist_trig_eff_lepphivslepphi_unmatched","emu_hist_trig_eff_lepptvsleppt_unmatched","emu_hist_trig_eff_lepetavslepeta_unmatched","emu_hist_trig_eff_lepphivslepphi_unmatched","mumu_hist_trig_eff_lepptvsleppt_unmatched","mumu_hist_trig_eff_lepetavslepeta_unmatched","mumu_hist_trig_eff_lepphivslepphi_unmatched"};
    
  TH2D *hist2d_orthomc[nhist];
  TH2D *hist2d_trigmc[nhist];
  TH2D *hist2d_trigmc_unmatch[nhist];
  TH2D *hist2d_orthodata[nhist];
  TH2D *hist2d_trigdata[nhist];
  TH2D *hist2d_trigdata_unmatch[nhist];
  
  TFile *filemc = new TFile(mcfile.c_str(),"read");
  TFile *filedata = new TFile(datafile.c_str(),"read");
  /*  for(int ihist=0; ihist<nhist;ihist++)
    {
      //if(historthonames[ihist][0] == 'e' && historthonames[ihist][1] == 'm' && historthonames[ihist][2] == 'u') continue;
      hist1d_orthomc[ihist] = (TH1D*)filemc->Get(historthonames[ihist].c_str());
      hist1d_orthomc[ihist]->Scale(data_lumi*mc_weight);
      hist1d_orthomc[ihist]->SetTitle((historthotitles[ihist] + " in mc").c_str());

      hist1d_trigmc[ihist] = (TH1D*)filemc->Get(histtrignames[ihist].c_str());
      hist1d_trigmc[ihist]->Scale(data_lumi*mc_weight);
      hist1d_trigmc[ihist]->SetTitle((histtrigtitles[ihist] + " in mc").c_str());
      
      hist1d_orthodata[ihist] = (TH1D*)filedata->Get(historthonames[ihist].c_str());
      hist1d_orthodata[ihist]->SetTitle((historthotitles[ihist] + " in data").c_str());
      
      hist1d_trigdata[ihist] = (TH1D*)filedata->Get(histtrignames[ihist].c_str());
      hist1d_trigdata[ihist]->SetTitle((histtrigtitles[ihist] + " in data").c_str());

      savehist(hist1d_orthomc[ihist], ("mc_" + historthonames[ihist]).c_str());
      savehist(hist1d_trigmc[ihist], ("mc_" + histtrignames[ihist]).c_str());
      savehist(hist1d_orthodata[ihist], ("data_" + historthonames[ihist]).c_str());
      savehist(hist1d_trigdata[ihist], ("data_" + histtrignames[ihist]).c_str());
      
      TEfficiency *Teff_mc = new TEfficiency(*hist1d_trigmc[ihist],*hist1d_orthomc[ihist]);
      Teff_mc->SetStatisticOption(TEfficiency::kFFC);
      Teff_mc->SetConfidenceLevel(0.683);
     
      //hist_eff_mc = (TH1D*)hist1d_orthomc[ihist]->Clone();
      //hist_eff_mc->Divide(hist1d_trigmc[ihist],hist1d_orthomc[ihist]);
      // cout<<ihist<<endl;
      savehist(Teff_mc, ("ratio_mc_" + historthonames[ihist]).c_str());

      //TH1D *hist_eff_data = (TH1D*)hist1d_orthodata[ihist]->Clone();
      //hist_eff_data->Divide(hist1d_trigdata[ihist],hist1d_orthodata[ihist]);
      //savehist(hist_eff_data, ("ratio_data_" + historthonames[ihist]).c_str());

      TEfficiency *Teff_data = new TEfficiency(*hist1d_trigdata[ihist],*hist1d_orthodata[ihist]);
      Teff_data->SetStatisticOption(TEfficiency::kFFC);
      Teff_data->SetConfidenceLevel(0.683);
      savehist(Teff_data, ("ratio_data_" + historthonames[ihist]).c_str());

      //      TEfficiency *Teff_sf = new TEfficiency(,);
      //Teff_sf->SetStatisticOption(TEfficiency::kFFC);
      //Teff_sf->SetConfidenceLevel(0.683);
      //savehist(Teff_sf, ("sf_" + historthonames[ihist]).c_str());
      //scale_factor(hist1d_orthomc[ihist],hist1d_trigmc[ihist],hist1d_orthodata[ihist],hist1d_trigdata[ihist]);
    }*/
  
  for(int ihist=0; ihist<nhist2d;ihist++)
    {
      if(ihist%3!=0) continue;

      //if(historthonames2d[ihist][0] == 'e' && historthonames2d[ihist][1] == 'm' && historthonames2d[ihist][2] == 'u') continue;
      hist2d_orthomc[ihist] = (TH2D*)filemc->Get(historthonames2d[ihist].c_str());
      hist2d_orthomc[ihist]->Scale(data_lumi*mc_weight);
      hist2d_orthomc[ihist]->SetTitle((historthotitles2d[ihist] + " in mc").c_str());

      hist2d_trigmc[ihist] = (TH2D*)filemc->Get(histtrignames2d[ihist].c_str());
      hist2d_trigmc[ihist]->Scale(data_lumi*mc_weight);
      hist2d_trigmc[ihist]->SetTitle((histtrigtitles2d[ihist] + " in mc").c_str());
      
      hist2d_orthodata[ihist] = (TH2D*)filedata->Get(historthonames2d[ihist].c_str());
      hist2d_orthodata[ihist]->SetTitle((historthotitles2d[ihist] + " in data").c_str());
      
      hist2d_trigdata[ihist] = (TH2D*)filedata->Get(histtrignames2d[ihist].c_str());
      hist2d_trigdata[ihist]->SetTitle((histtrigtitles2d[ihist] + " in data").c_str());

      if( historthonames2d[ihist].find("lepptvsleppt") != string::npos){
	hist2d_orthomc[ihist]->Rebin2D(2,2);
	hist2d_trigmc[ihist]->Rebin2D(2,2);
	hist2d_orthodata[ihist]->Rebin2D(2,2);
	hist2d_trigdata[ihist]->Rebin2D(2,2);
      }
      savehist(hist2d_orthomc[ihist], ("mc_" + historthonames2d[ihist]).c_str());
      savehist(hist2d_trigmc[ihist], ("mc_" + histtrignames2d[ihist]).c_str());
      savehist(hist2d_orthodata[ihist], ("data_" + historthonames2d[ihist]).c_str());
      savehist(hist2d_trigdata[ihist], ("data_" + histtrignames2d[ihist]).c_str());

      TEfficiency *Teff_mc = new TEfficiency(*hist2d_trigmc[ihist],*hist2d_orthomc[ihist]);
      Teff_mc->SetStatisticOption(TEfficiency::kFFC);
      Teff_mc->SetConfidenceLevel(0.683);

      savehist(Teff_mc, ("ratio_mc_" + historthonames2d[ihist]).c_str());

      TEfficiency *Teff_data = new TEfficiency(*hist2d_trigdata[ihist],*hist2d_orthodata[ihist]);
      Teff_data->SetStatisticOption(TEfficiency::kFFC);
      Teff_data->SetConfidenceLevel(0.683);

      savehist(Teff_data, ("ratio_data_" + historthonames2d[ihist]).c_str());

      /* TH2D *hist_sf;
      TEfficiency *Teff_sf = new TEfficiency(*hist_eff_data,*hist_eff_mc);
      Teff_sf->SetStatisticOption(TEfficiency::kFFC);
      Teff_sf->SetConfidenceLevel(0.683);
      hist_sf = (TH2D*)Teff_sf->CreateHistogram();
      //if(TEfficiency::CheckConsistency(*hist_eff_data,*hist_eff_mc))
	savehist(hist_sf, ("sf_" + historthonames2d[ihist]).c_str()); 
	//scale_factor(hist2d_orthomc[ihist],hist2d_trigmc[ihist],hist2d_orthodata[ihist],hist2d_trigdata[ihist]);*/
    }

  for(int ihist=0; ihist<nhist2d;ihist++)
    {
      if(ihist%3!=0) continue;
      
      //hist2d_orthomc[ihist] = (TH2D*)filemc->Get(historthonames2d[ihist].c_str());
      //hist2d_orthomc[ihist]->Scale(data_lumi*mc_weight);
      //hist2d_orthomc[ihist]->SetTitle((historthotitles2d[ihist] + " in mc_unmatch").c_str());
      
      hist2d_trigmc_unmatch[ihist] = (TH2D*)filemc->Get(histtrignames2d_unmatch[ihist].c_str());
      hist2d_trigmc_unmatch[ihist]->Scale(data_lumi*mc_weight);
      hist2d_trigmc_unmatch[ihist]->SetTitle((histtrigtitles2d[ihist] + " in mc_unmatch").c_str());

      //hist2d_orthodata[ihist] = (TH2D*)filedata->Get(historthonames2d[ihist].c_str());
      //hist2d_orthodata[ihist]->SetTitle((historthotitles2d[ihist] + " in data").c_str());
      
      hist2d_trigdata_unmatch[ihist] = (TH2D*)filedata->Get(histtrignames2d_unmatch[ihist].c_str());
      hist2d_trigdata_unmatch[ihist]->SetTitle((histtrigtitles2d[ihist] + " in data").c_str());
      
      if( historthonames2d[ihist].find("lepptvsleppt") != string::npos){
	//hist2d_orthomc[ihist]->Rebin2D(2,2);
	hist2d_trigmc_unmatch[ihist]->Rebin2D(2,2);
	//hist2d_orthodata[ihist]->Rebin2D(2,2);
	hist2d_trigdata_unmatch[ihist]->Rebin2D(2,2);
      }
      savehist(hist2d_orthomc[ihist], ("mc_unmatch_" + historthonames2d[ihist]).c_str());
      savehist(hist2d_trigmc_unmatch[ihist], ("mc_unmatch_" + histtrignames2d_unmatch[ihist]).c_str());

      int nxbin,nybin; 
      /*      TH2D *hist_eff_mc_unmatch = (TH2D*)hist2d_trigmc_unmatch[ihist]->Clone();
      
      nxbin = hist2d_orthomc[ihist]->GetNbinsX();
      nybin = hist2d_orthomc[ihist]->GetNbinsY();

      for(int ij=1; ij<=nxbin; ij++)
	{
	  for(int ik=1; ik<=nybin; ik++)
	    {
	      if(hist2d_orthomc[ihist]->GetBinContent(ij,ik)!=0)
		hist_eff_mc_unmatch->SetBinContent(ij,ik, hist2d_trigmc_unmatch[ihist]->GetBinContent(ij,ik) / hist2d_orthomc[ihist]->GetBinContent(ij,ik) );
	      else
		hist_eff_mc_unmatch->SetBinContent(ij,ik,0);
	    }
	    }
      savehist(hist_eff_mc_unmatch, ("ratio_mc_unmatch_" + historthonames2d[ihist]).c_str());

      nxbin = hist2d_orthodata[ihist]->GetNbinsX();
      nybin = hist2d_orthodata[ihist]->GetNbinsY();

      TH2D *hist_eff_data_unmatch;
      hist_eff_data_unmatch = (TH2D*)hist2d_orthodata[ihist]->Clone();
      for(int ij=1; ij<=nxbin; ij++)
	{
	  for(int ik=1; ik<=nybin; ik++)
	    {
	      if(hist2d_orthodata[ihist]->GetBinContent(ij,ik)!=0)
		hist_eff_data_unmatch->SetBinContent(ij,ik, hist2d_trigdata_unmatch[ihist]->GetBinContent(ij,ik) / hist2d_orthodata[ihist]->GetBinContent(ij,ik) );
	      else
		hist_eff_data_unmatch->SetBinContent(ij,ik,0);
	    }
	}
      savehist(hist_eff_data_unmatch, ("ratio_data_unmatch_" + historthonames2d[ihist]).c_str());

      scale_factor(hist2d_orthomc[ihist],hist2d_trigmc_unmatch[ihist],hist2d_orthodata[ihist],hist2d_trigdata_unmatch[ihist]);*/

      TH2D *hist_eff_mc_unmatch;
      TEfficiency *Teff_mc_unmatch = new TEfficiency(*hist2d_trigmc_unmatch[ihist],*hist2d_orthomc[ihist]);
      Teff_mc_unmatch->SetStatisticOption(TEfficiency::kFFC);
      Teff_mc_unmatch->SetConfidenceLevel(0.683);

      hist_eff_mc_unmatch = (TH2D*)Teff_mc_unmatch->CreateHistogram();
      savehist(hist_eff_mc_unmatch, ("ratio_mc_unmatch_" + historthonames2d[ihist]).c_str());

      nxbin = hist2d_orthodata[ihist]->GetNbinsX();
      nybin = hist2d_orthodata[ihist]->GetNbinsY();

      TH2D *hist_eff_data_unmatch;
      TEfficiency *Teff_data_unmatch = new TEfficiency(*hist2d_trigdata_unmatch[ihist],*hist2d_orthodata[ihist]);
      Teff_data_unmatch->SetStatisticOption(TEfficiency::kFFC);
      Teff_data_unmatch->SetConfidenceLevel(0.683);
      hist_eff_data_unmatch = (TH2D*)Teff_data_unmatch->CreateHistogram();

      savehist(hist_eff_data_unmatch, ("ratio_data_unmatch_" + historthonames2d[ihist]).c_str());

      /*TH2D *hist_sf_unmatch;
      TEfficiency *Teff_sf_unmatch = new TEfficiency(*hist_eff_data_unmatch,*hist_eff_mc_unmatch);
      Teff_sf_unmatch->SetStatisticOption(TEfficiency::kFFC);
      Teff_sf_unmatch->SetConfidenceLevel(0.683);
      hist_sf_unmatch = (TH2D*)Teff_sf_unmatch->CreateHistogram();
      if(TEfficiency::CheckConsistency(*hist_eff_data_unmatch,*hist_eff_mc_unmatch))
      savehist(hist_sf_unmatch, ("sf_unmatch_" + historthonames2d[ihist]).c_str()); */
    }
}



/*hist_eff_mc = (TH2D*)hist2d_orthomc[ihist]->Clone();

      nxbin = hist2d_orthomc[ihist]->GetNbinsX();
      nybin = hist2d_orthomc[ihist]->GetNbinsY();

      for(int ij=1; ij<=nxbin; ij++)
	{
	  for(int ik=1; ik<=nybin; ik++)
	    {
	      if(hist2d_orthomc[ihist]->GetBinContent(ij,ik)!=0)
		hist_eff_mc->SetBinContent(ij,ik, hist2d_trigmc[ihist]->GetBinContent(ij,ik) / hist2d_orthomc[ihist]->GetBinContent(ij,ik) );
	      else
		hist_eff_mc->SetBinContent(ij,ik,0);
	    }
	    }*/
      
