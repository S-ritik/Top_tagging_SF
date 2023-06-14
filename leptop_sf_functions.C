#ifndef LEPTOP_SF_FUNCTIONS_H
#define LEPTOP_SF_FUNCTIONS_H

#include "plotscripts.h"


void plot_1dhists(int nhist, TH1D *hist_obs[nhist], string dataname[nhist],string plotname)
{
  
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",hist_obs[0],8,0);
  cv->cd();
  legv = tdrLeg(0.2,0.65,0.875,0.915);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);

  double max_bincontent = hist_obs[0]->GetMaximum();
  
  for(int ih = 0; ih < nhist; ih++)
    {
      //hist_obs[ih]->Scale(1.0/max(0.00001,hist_obs[ih]->Integral()));
      hist_obs[ih]->SetFillStyle(0);
      hist_obs[ih]->SetFillColor(0);
      hist_obs[ih]->SetLineColor(ih+2);
      //     hist_obs[ih]->SetLineStyle(ih+2);
      hist_obs[ih]->SetLineWidth(2);
      hist_obs[ih]->GetYaxis()->SetTitleOffset(1.5);

      legv->AddEntry(hist_obs[ih],dataname[ih].c_str(),"lep");

      max_bincontent = max(max_bincontent,double(hist_obs[ih]->GetMaximum()));
      hist_obs[ih]->SetMinimum(0.01);
    }

  hist_obs[0]->SetMaximum(1.5*max_bincontent);
  
  //gPad->SetLogy(1);

  hist_obs[0]->Draw("hist");
  for(int ih = 1; ih < nhist; ih++)
    hist_obs[ih]->Draw("histsame");

  legv->Draw("same");
  CMS_lumi( cv, 8, 0 );

  cv->SaveAs(("plots/" + plotname).c_str());
}

void substract_MC(fileinfo mc, string histname_regionA,string histname_regionB, string histname_regionC, string histname_regionD, TH1D *hist_regionA, TH1D *hist_regionB, TH1D *hist_regionC, TH1D *hist_regionD)
{
  TFile *file_mc = new TFile(mc.out_name.c_str(),"read");
  
  if(!file_mc->IsOpen())
    {cout<<"Check: file "<<mc.out_name<<" not present"<<endl;exit(0);}
  
  TH1D *tmp_regionA = (TH1D*)file_mc->Get(histname_regionA.c_str());
  TH1D *tmp_regionB = (TH1D*)file_mc->Get(histname_regionB.c_str());
  TH1D *tmp_regionC = (TH1D*)file_mc->Get(histname_regionC.c_str());
  TH1D *tmp_regionD = (TH1D*)file_mc->Get(histname_regionD.c_str());

  if(tmp_regionA == NULL)
    {cout<<"Check: file "<<mc.out_name<<" does not have 1D Hist "<<histname_regionA<<endl;exit(0);}
  if(tmp_regionB == NULL)
    {cout<<"Check: file "<<mc.out_name<<" does not have 1D Hist "<<histname_regionB<<endl;exit(0);}
  if(tmp_regionC == NULL)
    {cout<<"Check: file "<<mc.out_name<<" does not have 1D Hist "<<histname_regionC<<endl;exit(0);}
  if(tmp_regionD == NULL)
    {cout<<"Check: file "<<mc.out_name<<" does not have 1D Hist "<<histname_regionD<<endl;exit(0);}

  hist_regionA->Add(tmp_regionA,-1);
  hist_regionB->Add(tmp_regionB,-1);
  hist_regionC->Add(tmp_regionC,-1);
  hist_regionD->Add(tmp_regionD,-1);

  file_mc->Close();
}

void compute_QCD_ABCD(string hist2dname, string hist_struct, string fitvar, string fitvartitle, string baseline_sel, string sr_cr_cut[2], string transfer_factor_cut[2], string weight, fileinfo data, fileinfo QCD, vector<fileinfo> remainingMCfiles)
{

  /// Checking corelation of variables used in ABCD method; Should be weakly corelated
  
  TFile *file_data = new TFile(data.out_name.c_str(),"read");

  if(!file_data->IsOpen())
    {cout<<"Check: file "<<data.out_name<<" not present"<<endl;exit(0);}

  TH2D* hist2d_data = (TH2D*)file_data->Get(hist2dname.c_str());

  if(hist2d_data == NULL)
    {cout<<"Check: file "<<data.out_name<<" does not have 2D Hist "<<hist2dname<<endl;exit(0);}

  string plotname = "fordatadivenqcd_" +data.out_name.substr(30,9)+ "_" + hist2dname;
  cout<<"Correlation for hist"<<plotname<<" = "<<hist2d_data->GetCorrelationFactor()<<endl;
  savehist2d(hist2d_data,plotname,false);

  file_data->Close();
  
  TFile *file_qcd = new TFile(qcd.out_name.c_str(),"read");
    
  if(!file_qcd->IsOpen())
    {cout<<"Check: file "<<qcd.out_name<<" not present"<<endl;exit(0);}
  
  TH2D* hist2d_qcd = (TH2D*)file_qcd->Get(hist2dname.c_str());

  if(hist2d_qcd == NULL)
    {cout<<"Check: file "<<qcd.out_name<<" does not have 2D Hist "<<hist2dname<<endl;exit(0);}

  plotname = "fordatadivenqcd_" +qcd.out_name.substr(30,9)+ "_" + hist2dname;
  cout<<"Correlation for hist"<<plotname<<" = "<<hist2d_qcd->GetCorrelationFactor()<<endl;
  savehist2d(hist2d_qcd,plotname,false);


  ///////////////////  Get bining of histogram from another stored histogram
  
  hist_struct = "sdmass_probejet_fail_topmatch_eltop";

  TH1D* hist_copy_struct = (TH1D*)file_qcd->Get(hist_struct.c_str());
  
  if(hist_copy_struct == NULL)
    {cout<<"Check: file "<<qcd.out_name<<" does not have 1D Hist "<<hist_struct<<" which is needed to define the bining of histograms"<<endl;exit(0);}
  hist_copy_struct->SetDirectory(0);  

  file_qcd->Close();


  /*
    ABCD method
     
    For estimate in signal region
    
    |      |
    |  A   |  B
    |______|_____
    |      |
    |  C   |  D
    |______|______ 

    A, B = region for transfer factors
    C = Control region
    D = Signal region = C*B/A
  */

  
  ////////////// For defining histogram in every region in all MC samples and data  /////////////////////////
  
  string evt_sel_regionA = baseline_sel + " && " + transfer_factor_cut[0] + " && " + sr_cr_cut[0];
  string evt_sel_regionB = baseline_sel + " && " + transfer_factor_cut[0] + " && " + sr_cr_cut[1];
  string evt_sel_regionC = baseline_sel + " && " + transfer_factor_cut[1] + " && " + sr_cr_cut[0];
  string evt_sel_regionD = baseline_sel + " && " + transfer_factor_cut[1] + " && " + sr_cr_cut[1];

  cout<<evt_sel_regionA<<endl<<evt_sel_regionB<<endl<<evt_sel_regionC<<endl<<evt_sel_regionD<<endl;
    
  string histname_regionA = "datadrivenABCD_regionA_" + fitvar;
  string histname_regionB = "datadrivenABCD_regionB_" + fitvar;
  string histname_regionC = "datadrivenABCD_regionC_" + fitvar;
  string histname_regionD = "datadrivenABCD_regionD_" + fitvar;
	
  TH1D *hist_regionA = (TH1D*)hist_copy_struct->Clone(histname_regionA.c_str());
  TH1D *hist_regionB = (TH1D*)hist_copy_struct->Clone(histname_regionB.c_str());
  TH1D *hist_regionC = (TH1D*)hist_copy_struct->Clone(histname_regionC.c_str());
  TH1D *hist_regionD = (TH1D*)hist_copy_struct->Clone(histname_regionD.c_str());

  vector<fileinfo> allfiles = {data};
  allfiles.insert(allfiles.end(),qcd);
  allfiles.insert(allfiles.end(),remainingMCfiles.begin(),remainingMCfiles.end());
  
  makerootfile(fitvar,weight,evt_sel_regionA,hist_regionA,histname_regionA,fitvartitle.c_str(),allfiles);

  makerootfile(fitvar,weight,evt_sel_regionB,hist_regionB,histname_regionB,fitvartitle.c_str(),allfiles);
  
  makerootfile(fitvar,weight,evt_sel_regionC,hist_regionC,histname_regionC,fitvartitle.c_str(),allfiles);
  
  makerootfile(fitvar,weight,evt_sel_regionD,hist_regionD,histname_regionD,fitvartitle.c_str(),allfiles);
  

  ////////////////     Get histogram in each region to compare data and QCD yield    /////////////////////////
  
  TH1D *hist_qcd_regionA = (TH1D*)hist_copy_struct->Clone(histname_regionA.c_str());
  TH1D *hist_qcd_regionB = (TH1D*)hist_copy_struct->Clone(histname_regionB.c_str());
  TH1D *hist_qcd_regionC = (TH1D*)hist_copy_struct->Clone(histname_regionC.c_str());
  TH1D *hist_qcd_regionD = (TH1D*)hist_copy_struct->Clone(histname_regionD.c_str());

  TFile* file_qcd2 = new TFile(qcd.out_name.c_str(),"read");
    
  if(!file_qcd2->IsOpen())
    {cout<<"Check: file "<<qcd.out_name<<" not present"<<endl;exit(0);}
  
  hist_qcd_regionA = (TH1D*)file_qcd2->Get(histname_regionA.c_str());
  hist_qcd_regionB = (TH1D*)file_qcd2->Get(histname_regionB.c_str());
  hist_qcd_regionC = (TH1D*)file_qcd2->Get(histname_regionC.c_str());
  hist_qcd_regionD = (TH1D*)file_qcd2->Get(histname_regionD.c_str());
  
  if(hist_qcd_regionA == NULL)
    {cout<<"Check: file "<<qcd.out_name<<" does not have 1D Hist "<<histname_regionA<<endl;exit(0);}
  if(hist_qcd_regionB == NULL)
    {cout<<"Check: file "<<qcd.out_name<<" does not have 1D Hist "<<histname_regionB<<endl;exit(0);}
  if(hist_qcd_regionC == NULL)
    {cout<<"Check: file "<<qcd.out_name<<" does not have 1D Hist "<<histname_regionC<<endl;exit(0);}
  if(hist_qcd_regionD == NULL)
    {cout<<"Check: file "<<qcd.out_name<<" does not have 1D Hist "<<histname_regionD<<endl;exit(0);}
    
  hist_regionA->Reset("ICES");
  hist_regionB->Reset("ICES");
  hist_regionC->Reset("ICES");
  hist_regionD->Reset("ICES");

  hist_qcd_regionA->SetDirectory(0);
  hist_qcd_regionB->SetDirectory(0);
  hist_qcd_regionC->SetDirectory(0);
  hist_qcd_regionD->SetDirectory(0);

  file_qcd2->Close();
  
  TFile* file_data2 = new TFile(data.out_name.c_str(),"read");
  
  if(!file_data2->IsOpen())
    {cout<<"Check: file "<<data.out_name<<" not present"<<endl;exit(0);}

  hist_regionA = (TH1D*)file_data2->Get(histname_regionA.c_str());
  hist_regionB = (TH1D*)file_data2->Get(histname_regionB.c_str());
  hist_regionC = (TH1D*)file_data2->Get(histname_regionC.c_str());
  hist_regionD = (TH1D*)file_data2->Get(histname_regionD.c_str());

  hist_regionA->SetDirectory(0);
  hist_regionB->SetDirectory(0);
  hist_regionC->SetDirectory(0);
  hist_regionD->SetDirectory(0);

  if(hist_regionA == NULL)
    {cout<<"Check: file "<<data.out_name<<" does not have 1D Hist "<<histname_regionA<<endl;exit(0);}
  if(hist_regionB == NULL)
    {cout<<"Check: file "<<data.out_name<<" does not have 1D Hist "<<histname_regionB<<endl;exit(0);}
  if(hist_regionC == NULL)
    {cout<<"Check: file "<<data.out_name<<" does not have 1D Hist "<<histname_regionC<<endl;exit(0);}
  if(hist_regionD == NULL)
    {cout<<"Check: file "<<data.out_name<<" does not have 1D Hist "<<histname_regionD<<endl;exit(0);}

  file_data2->Close();

  cout<<"Check 1";
  for(int ifile=0; ifile < (int)remainingMCfiles.size(); ifile++)
    substract_MC(remainingMCfiles[ifile], histname_regionA, histname_regionB, histname_regionC, histname_regionD, hist_regionA, hist_regionB, hist_regionC, hist_regionD);

  ////////////////     Plot histogram in each region to compare data-allMC and QCD distribution    /////////////////////////

  TH1D *hist_arr[2]; string dataname[2]; int col[1];
  
  for(int i =1; i<hist_regionA->GetNbinsX()+1; i++)
    cout<<"bin "<<i<<" "<<hist_regionA->GetBinLowEdge(i)<<" to "<<hist_regionA->GetBinLowEdge(i)+hist_regionA->GetBinWidth(i)<<" "<<hist_regionA->GetBinContent(i)<<endl;
  cout<<"Check 2";
   
  hist_arr[0] = hist_regionA;   hist_arr[1] = hist_qcd_regionA;
  dataname[0] = "Data - Remaining MC"; dataname[1] = "MC QCD"; col[0] = kRed;
  plot_datamcratio(2, hist_arr, dataname, col, histname_regionA + "_compare", fitvartitle);

  hist_arr[0] = hist_regionB;   hist_arr[1] = hist_qcd_regionB;
  plot_datamcratio(2, hist_arr, dataname, col, histname_regionB + "_compare", fitvartitle);

  hist_arr[0] = hist_regionC;   hist_arr[1] = hist_qcd_regionC;
  plot_datamcratio(2, hist_arr, dataname, col, histname_regionC + "_compare", fitvartitle);

  hist_arr[0] = hist_regionD;   hist_arr[1] = hist_qcd_regionD;
  plot_datamcratio(2, hist_arr, dataname, col, histname_regionD + "_compare", fitvartitle);
  cout<<"Check 3"<<endl;

  string dataname_eltop[9] = {"Data","QCD","t #bar{t} (el+jets)","t #bar{t} (other)","Single t","Z+j","W+j","DiBoson","ttX"};
  int col_all[8]  = {kRed,kMagenta-7,kGreen+4,kBlue,kOrange+7,kCyan+1,kRed-5,kViolet+6};

  plot_datamcratio(9,dataname_eltop,col_all,histname_regionA,fitvartitle,allfiles);
  plot_datamcratio(9,dataname_eltop,col_all,histname_regionB,fitvartitle,allfiles);
  plot_datamcratio(9,dataname_eltop,col_all,histname_regionC,fitvartitle,allfiles);
  plot_datamcratio(9,dataname_eltop,col_all,histname_regionD,fitvartitle,allfiles);


  cout<<"Check 4"<<endl;

  TH1D *hist_datadriven_transfer_factor_plot = (TH1D*)hist_regionB->Clone();
  hist_datadriven_transfer_factor_plot->Divide(hist_regionA);

  TH1D * hist_datadriven_estimated_sr_plot = (TH1D*)hist_regionC->Clone();
  hist_datadriven_estimated_sr_plot->Divide(hist_regionA);
  hist_datadriven_estimated_sr_plot->Multiply(hist_regionB);

  cout<<"Check 5"<<endl;
  TFile *fout = new TFile(qcd.out_name.c_str(),"update"); cout<<"Check 5.5"<<endl;
  fout->cd(); 
  hist_datadriven_estimated_sr_plot->Write(("hist_datadriven_estimated_" + hist_struct).c_str(),TObject::kOverwrite);
  cout<<"Check 6"<<endl;
  fout->Close();
  
  /*  for(int ig=0;ig<qcd.in_names.size();ig++)
      {
      TFile *filein = new TFile(qcd.in_names[ig].c_str(),"update");
      filein->cd();
      float weight_datadriven = 0;
      float fitvar_value = 1;
      
      TTree *tree = (TTree*)filein->Get("newvars");
      tree->SetBranchAddress(fitvar.c_str(),&fitvar_value);

      tree->Branch("weight_datadriven",&weight_datadriven,"weight_datadriven/F");

      Int_t nentries = (Int_t)tree->GetEntries();
      for (Int_t ij=0;ij < nentries; ij++) {
      tree->GetEntry(ij);
      int ibin = 0;
      if(fitvar_value > hist_datadriven_transfer_factor_plot->GetXaxis()->GetBinUpEdge(hist_datadriven_transfer_factor_plot->GetXaxis()->GetNbins() - 1))
      {  ibin = hist_datadriven_transfer_factor_plot->GetXaxis()->GetNbins();}
      else  if(fitvar_value < hist_datadriven_transfer_factor_plot->GetXaxis()->GetBinLowEdge(1))
      {ibin = 1;}
      else
      { ibin = hist_datadriven_transfer_factor_plot->GetXaxis()->FindFixBin(fitvar_value); }
      weight_datadriven = hist_datadriven_transfer_factor_plot->GetBinContent(ibin);
      tree->Fill();
      }
      filein->Close();
      }*/
}

float* calculate_sf_fromhiggscombine(string channel, const int nprocs, string processes[nprocs], string channelname[nprocs], string ptbins, int iptbins, TGraphAsymmErrors *gr_sf)
{
  float sf_witherror[3];  /// 0 = sf;  1 = sf error up;  2 = sf error down
  
  string regions[2] = {"fail","pass"};
  float prefit_failsig_evts[nprocs], postfit_failsig_evts[nprocs], prefit_passsig_evts[nprocs], postfit_passsig_evts[nprocs];
  float eff_mc[nprocs];
  float eff_data[nprocs];
  
  for(int ir = 0 ; ir < 2; ir++)
    {
      
      string region = regions[ir];

      using namespace RooFit;
      
      ///  Open higgs combine output file to get fitted plots
      TFile *fitfile;
      string filename = "fitDiagnostics_" + channel + ptbins + ".root";
      fitfile = new TFile(filename.c_str(),"read");
      if(!fitfile->IsOpen())
        {
	  cout<<"Check: Higgs combine outout file "<<filename<<" not present. Or check its name should be in format fitDiagnostics_*channel*.root"<<endl;exit(0);
        }
      
      ///  Extract histograms with improper bining
      
      string fitdirname = "shapes_fit_s/" + region + channel + "/";
      
      TH1D* hist_result_fit_improperbin = (TH1D*)fitfile->Get((fitdirname + "total").c_str());
      TH1D* hist_sig_fit_improperbin = (TH1D*)fitfile->Get((fitdirname + "total_signal").c_str());
      TH1D* hist_bkg_fit_improperbin = (TH1D*)fitfile->Get((fitdirname + "total_background").c_str());
      
      TH1D *hist_procs_postfit[nprocs];   //// for properly binned hisotgram of each process
      TH1D *hist_procs_postfit_improperbin[nprocs];
      for(int ip = 0; ip < nprocs; ip++)
	{
	  string prochistname = "shapes_fit_s/" + region + channel + "/" + processes[ip];
	  hist_procs_postfit_improperbin[ip] = (TH1D*)fitfile->Get(prochistname.c_str());
	  if(hist_procs_postfit_improperbin[ip] == NULL)
	    {cout<<"Check: file "<<filename<<" does not have histogram for process  "<<processes[ip]<<" with hist name:"<<prochistname<<endl;exit(0);}
	}
      
      if(hist_result_fit_improperbin == NULL || hist_bkg_fit_improperbin == NULL || hist_sig_fit_improperbin == NULL)
	{       cout<<"Check: Either fit dir: ("<<fitdirname<<") not present or the histogram (total) is not present in file "<<filename<<". Make sure that ----saveShapes option was used in higgs combine to get this."<<endl;exit(0);}
      
      //// Use orgional root given to combine to get proper binned histogram
      TFile *filein;
      string datafilename = "combine_" + channel + ptbins + ".root";
      filein = new TFile(datafilename.c_str(),"read");
      if(!filein->IsOpen())
	{cout<<"Check: Input file  to higgs combine  with name "<<datafilename<<" is not present which is needed to get the data histogram (cant taken from combine output as the x asis range will be not correct). Make sure its name is in format combine_*channel*.root"<<endl;exit(0);}
      string datahistname = "sdmass_probejet_" + region + "_data_obs_" + channel;
      filein->cd();
      TH1D *hist_data = (TH1D*)filein->Get(datahistname.c_str());
      cout<<" check data name = "<<datahistname<<" "<<hist_data->GetName()<<" "<<hist_data->GetTitle()<<endl;
      if(hist_data == NULL)
	{cout<<"Check: file "<<datafilename<<" does not have histogram "<<datahistname<<endl;exit(0);}

      string sys_names[] = {"pileup",
			  "btag",
			  "hadtagger",
			  "JER",
			  "alpha",
			  //"pdf",
			  //			  "pdf1","pdf2","pdf3","pdf4","pdf5","pdf6","pdf7","pdf8","pdf9","pdf10","pdf11","pdf12","pdf13","pdf14","pdf15","pdf16","pdf17","pdf18","pdf19","pdf20","pdf21","pdf22","pdf23","pdf24","pdf25","pdf26","pdf27","pdf28","pdf29","pdf30","pdf31","pdf32","pdf33","pdf34","pdf35","pdf36","pdf37","pdf38","pdf39","pdf40","pdf41","pdf42","pdf43","pdf44","pdf45","pdf46","pdf47","pdf48","pdf49","pdf50","pdf51","pdf52","pdf53","pdf54","pdf55","pdf56","pdf57","pdf58","pdf59","pdf60","pdf61","pdf62","pdf63","pdf64","pdf65","pdf66","pdf67","pdf68","pdf69","pdf70","pdf71","pdf72","pdf73","pdf74","pdf75","pdf76","pdf77","pdf78","pdf79","pdf80","pdf81","pdf82","pdf83","pdf84","pdf85","pdf86","pdf87","pdf88","pdf89","pdf90","pdf91","pdf92","pdf93","pdf94","pdf95","pdf96","pdf97","pdf98","pdf99","pdf100","pdf101","pdf102",
			  "partonshower",
			  "JES_AbsoluteStat",
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
			  "JES_TimePtEta"};
      int nsys_err = sizeof(sys_names)/sizeof(sys_names[0]); ;
      
      TH1D *hist_procs_prefit[nprocs];
      TH1D *hist_procs_prefit_sys[nprocs][nsys_err][2];
      TH1D *hist_procs_prefit_relsyserror[nprocs][nsys_err];
      TH1D *hist_procs_prefit_total;
      for(int ip = 0; ip < nprocs; ip++)
	{
	  string prochistname = "sdmass_probejet_" + region + "_" + processes[ip] + "_" + channel;
	  hist_procs_prefit[ip] = (TH1D*)filein->Get(prochistname.c_str());
	  cout<<" check mc name = "<<hist_procs_prefit[ip]->GetName()<<" "<<hist_procs_prefit[ip]->GetTitle()<<endl;
	  if(hist_procs_prefit[ip] == NULL)
	    {cout<<"Check: file "<<datafilename<<" does not have histogram for process  "<<processes[ip]<<" with hist name:"<<prochistname<<endl;exit(0);}
	  for(int isys = 0; isys < nsys_err; isys++)
	    {
	      string prochistname_sys = "sdmass_probejet_" + region + "_" + processes[ip] + "_" + channel + "_" + sys_names[isys];
	      hist_procs_prefit_sys[ip][isys][0] = (TH1D*)filein->Get((prochistname_sys + "Up").c_str());
	      hist_procs_prefit_sys[ip][isys][1] = (TH1D*)filein->Get((prochistname_sys + "Down").c_str());
	      hist_procs_prefit_relsyserror[ip][isys] = (TH1D*)hist_procs_prefit[ip]->Clone();	  
	    }
	  if(ip==0)
	    hist_procs_prefit_total = (TH1D*)hist_procs_prefit[ip]->Clone();
	  //  else
	  //  hist_procs_prefit_total->Add(hist_procs_prefit[ip]);
	}
      
      
      for(int ib = 1; ib < hist_procs_prefit_total->GetNbinsX() +1; ib++)
	{
	  double err = 0;
	  for(int ip = 0; ip < nprocs; ip++)
	    {
	      double nom_yield = hist_procs_prefit[ip]->GetBinContent(ib); 
	      for(int isys = 0; isys < nsys_err; isys++)
		{
		  cout<<region<<" sys  "<<sys_names[isys]<<endl;
		  double sys_error = (hist_procs_prefit_sys[ip][isys][0]->GetBinContent(ib) - hist_procs_prefit_sys[ip][isys][1]->GetBinContent(ib))/nom_yield;
		  err +=  sys_error*sys_error;
		  hist_procs_prefit_relsyserror[ip][isys]->SetBinContent(ib,sys_error*100);
		  hist_procs_prefit_relsyserror[ip][isys]->SetBinError(ib,0);
		  
		}
	    }
	  cout<<"bin = "<<ib<<" content = "<<hist_procs_prefit_total->GetBinContent(ib)<<" error = "<<err;
	  hist_procs_prefit_total->SetBinError(ib,sqrt(err));
	}

      string dataname[3] = {"topmatch","blepmatch","nomatch"};
      for(int isys = 0; isys < nsys_err; isys++)
	{
	  TH1D *giv_hist[3];
	  giv_hist[0] = hist_procs_prefit_relsyserror[0][isys];
	  giv_hist[1] = hist_procs_prefit_relsyserror[1][isys];
	  giv_hist[2] = hist_procs_prefit_relsyserror[2][isys];
	  cout<<endl<<"sys plot "<<isys<<endl<<endl;
	  plot_1dhists(3,giv_hist,dataname,("/sysematic/" + string(hist_procs_prefit[0]->GetTitle()) + "_relerr_" + sys_names[isys] + ".png"));
	}
      string total_MC_name = "total MC with all error"; 
      plot_1dhists(1,&hist_procs_prefit_total,&total_MC_name,("/sysematic/" + string(hist_procs_prefit[0]->GetTitle()) + "_allerr.png"));

      /// correct the visualizing of the histogram (useful only for variable bin width histogram).
      int nbins_data = hist_data->GetNbinsX();
      cout<<nbins_data<<endl;
      if(nbins_data != hist_result_fit_improperbin->GetNbinsX()) {cout<<"Check input file names! The no. of bins in prefit data is not equal to no. of points in roohist data of higgscombine output files"<<endl;exit(0);}
      
      
      cout<<"\n"<<region<<" "<<channel<<" \nData prefit = "<<hist_data->Integral()<<" Prefit signal yield = "<<hist_procs_prefit[0]->Integral()<<" \nPostfit yield (sig+bkg) = "<<hist_result_fit_improperbin->Integral()<<" signal yield = "<<hist_sig_fit_improperbin->Integral()<<endl<<endl;

      /// Get proper bining from orginial histograms
      for(int ip = 0; ip < nprocs; ip++)
	{
	  hist_procs_postfit[ip] = (TH1D*)hist_data->Clone();
	}
      TH1D* hist_result_fit = (TH1D*)hist_data->Clone();
      TH1D* hist_sig_fit = (TH1D*)hist_data->Clone();
      TH1D* hist_bkg_fit = (TH1D*)hist_data->Clone();

      int event_perbin = 10;
      for(int ib = 0; ib < nbins_data; ib++)
	{
	  float ith_databin_width = hist_data->GetBinWidth(ib+1);
	  //float ith_databin_width = 10;
	  
	  cout<<hist_data->GetBinContent(ib+1)<<"*"<<event_perbin<<" / "<<ith_databin_width;

	  hist_data->SetBinContent(ib+1,hist_data->GetBinContent(ib+1)*event_perbin/ith_databin_width);
	  hist_data->SetBinError(ib+1,hist_data->GetBinError(ib+1)*event_perbin/ith_databin_width);
	  cout<<" = "<<hist_data->GetBinContent(ib+1)<<endl;

	  cout<<hist_result_fit_improperbin->GetBinContent(ib+1)<<"*"<<event_perbin<<" / "<<ith_databin_width;

	  hist_result_fit->SetBinContent(ib+1,hist_result_fit_improperbin->GetBinContent(ib+1)*event_perbin/ith_databin_width);
	  hist_result_fit->SetBinError(ib+1,hist_result_fit_improperbin->GetBinError(ib+1)*event_perbin/ith_databin_width);
	  hist_result_fit->SetLineWidth(3);

	  cout<<" = "<<hist_result_fit->GetBinContent(ib+1)<<endl;

	  hist_sig_fit->SetBinContent(ib+1,hist_sig_fit_improperbin->GetBinContent(ib+1)*event_perbin/ith_databin_width);
	  hist_sig_fit->SetBinError(ib+1,hist_sig_fit_improperbin->GetBinError(ib+1)*event_perbin/ith_databin_width);

	  hist_bkg_fit->SetBinContent(ib+1,hist_bkg_fit_improperbin->GetBinContent(ib+1)*event_perbin/ith_databin_width);
	  hist_bkg_fit->SetBinError(ib+1,hist_bkg_fit_improperbin->GetBinError(ib+1)*event_perbin/ith_databin_width);

	  for(int ip = 0; ip < nprocs; ip++)
	    {
	      hist_procs_postfit[ip]->SetBinContent(ib+1,hist_procs_postfit_improperbin[ip]->GetBinContent(ib+1)*event_perbin/ith_databin_width);
	      hist_procs_postfit[ip]->SetBinError(ib+1,hist_procs_postfit_improperbin[ip]->GetBinError(ib+1)*event_perbin/ith_databin_width);
	      hist_procs_postfit[ip]->SetLineStyle(1);
	      hist_procs_postfit[ip]->SetLineWidth(3);
	      hist_procs_postfit[ip]->SetLineColor(ip+6);

	      if(ip == 1)
		cout<<hist_data->GetBinContent(ib+1)<<" = "<<hist_procs_prefit[ip]->GetBinContent(ib+1)<<endl;

	      hist_procs_prefit[ip]->SetBinContent(ib+1,hist_procs_prefit[ip]->GetBinContent(ib+1)*event_perbin/ith_databin_width);

	      hist_procs_prefit[ip]->SetLineStyle(2);
	      hist_procs_prefit[ip]->SetLineWidth(3);
	      hist_procs_prefit[ip]->SetLineColor(ip+6);
	    }
	}

      ofstream f1;
      f1.open("plots/combine_sf_yield.txt",ios::app);
      f1<<"Channel "<<channel<<":"<<processes[0]<<"     Condition: Probe "<<region<<" region for "<<ptbins<<endl;
      f1<<"Data : "<<hist_data->Integral()<<endl;
      
      float nbkg = 0;
      for(int ip = 0; ip < nprocs; ip++)
	nbkg+=hist_procs_prefit[ip]->Integral();
      
      f1<<"Prefit  yeilds        Signal : "<<hist_procs_prefit[0]->Integral()<<" Total background : "<<nbkg<<endl;
      f1<<"Postfit yeilds        Signal : "<<hist_procs_postfit_improperbin[0]->Integral()<<" Total background : "<<hist_bkg_fit_improperbin->Integral()<<endl;
      f1.close();

      if(ir == 0)
	{
	  for(int ip = 0; ip < nprocs; ip++)
	    {
	      prefit_failsig_evts[ip] = hist_procs_prefit[ip]->Integral();
	      postfit_failsig_evts[ip] = hist_procs_postfit_improperbin[ip]->Integral();
	    }
	}
      else if(ir == 1)
	{
	  for(int ip = 0; ip < nprocs; ip++)
	    {
	      prefit_passsig_evts[ip] = hist_procs_prefit[ip]->Integral();
	      postfit_passsig_evts[ip] = hist_procs_postfit_improperbin[ip]->Integral();
	    }
	}
      for(int ip = 0; ip < nprocs; ip++)
	{
	  eff_mc[ip] = prefit_passsig_evts[ip] / max(float(0.00001),(prefit_passsig_evts[ip] + prefit_failsig_evts[ip]));
	  eff_data[ip] = postfit_passsig_evts[ip] / max(float(0.00001),(postfit_passsig_evts[ip] + postfit_failsig_evts[ip]));
	}
      sf_witherror[0] = eff_data[0]/max(float(0.00001),eff_mc[0]);
      sf_witherror[1] = 0;
      sf_witherror[2] = 0;

      if(ir == 1 && processes[0] == "topmatch")
	{
	  TTree * tree_fit = (TTree*)fitfile->Get("tree_fit_sb");
	  double sf,sf_up,sf_down;
	  tree_fit->SetBranchAddress("SF",&sf);
	  tree_fit->SetBranchAddress("SFHiErr",&sf_up);
	  tree_fit->SetBranchAddress("SFLoErr",&sf_down);

	  tree_fit->GetEntry(0);

	  sf_witherror[0] = (float)sf;
	  sf_witherror[1] = (float)sf_up;
	  sf_witherror[2] = (float)sf_down;
	}
      
      hist_data->SetFillStyle(0);
      hist_data->SetFillColor(0);
      hist_data->SetMarkerStyle(kFullCircle);
      hist_data->SetMarkerColor(kBlack);
      hist_data->SetMarkerSize(0.7);
      hist_data->SetLineColor(kBlack);
      hist_data->SetLineWidth(3);
      hist_data->SetLineStyle(1);

      hist_result_fit->SetLineColor(kBlue);

      hist_sig_fit->SetLineColor(kGreen);

      hist_bkg_fit->SetLineColor(kRed);

      TGraphAsymmErrors *h_ratio = new TGraphAsymmErrors(hist_data,hist_result_fit,"pois");

      TLegend *legend = tdrLeg(0.54,0.6,0.84,0.85);
      legend->SetBorderSize(0);
      legend->SetTextFont(42);
      legend->SetTextSize(0.035);
      legend->AddEntry(hist_data, "Data", "6ple1");
      legend->AddEntry(hist_result_fit, "Total Fit", "l");
      for(int ip = 0; ip < nprocs; ip++)
        legend->AddEntry(hist_procs_postfit[ip], channelname[ip].c_str(), "l");
      legend->AddEntry(hist_procs_prefit[0],"Prefit","l");
      legend->AddEntry(hist_procs_postfit[0],"Postfit","l");

      float max_mc_binstat = 0;
      for(int ip = 0; ip < nprocs; ip++)
        max_mc_binstat = max((float)hist_procs_postfit[ip]->GetMaximum(),max((float)hist_procs_prefit[ip]->GetMaximum(),max_mc_binstat));


      TCanvas * can = tdrDiCanvas("can",hist_data,hist_data,8,0);
      can->cd(1);
      gPad->SetLogy(0);
      hist_data->GetYaxis()->SetTitle(("Events / " + to_string(event_perbin) + " GeV").c_str());
      hist_data->GetXaxis()->SetTitle(" ");
      hist_data->SetMinimum(0.001);
      hist_data->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)*max((float)hist_data->GetMaximum(),max_mc_binstat));
      hist_data->GetYaxis()->CenterTitle(1);
      hist_data->GetXaxis()->CenterTitle(1);
      hist_data->GetYaxis()->SetTitleOffset(0.7);
      hist_data->GetYaxis()->SetTitleSize(0.09);
      hist_data->GetYaxis()->SetLabelSize(0.04);

      hist_data->Draw("P");
      hist_result_fit->Draw("histsame");
      for(int ip = 0; ip < nprocs; ip++)
	{
	  hist_procs_postfit[ip]->Draw("histsame");
	  hist_procs_prefit[ip]->Draw("histsame");
	}
      legend->Draw("same");
      hist_data->Draw("PE1SAME");

      can->cd(2);
      h_ratio->SetMarkerColor(1);
      h_ratio->SetMarkerStyle(20);
      h_ratio->SetMarkerSize(0.89);
      h_ratio->SetLineColor(kBlack);
      int nbin = hist_data->GetXaxis()->GetNbins();
      double* bins =const_cast <double *>(hist_data->GetXaxis()->GetXbins()->GetArray());
      TH1D *band = new TH1D("Band", "", nbin, bins);

      for(int i=0;i<nbin; i++)
	{
	  band->SetBinContent(i+1, 1.0);
	  double err = 0;
	  if(hist_data->GetBinContent(i+1)!=0 && hist_result_fit->GetBinContent(i+1)!=0)
	    err=(hist_result_fit->GetBinError(i+1))/(hist_result_fit->GetBinContent(i+1));
	  else if(hist_data->GetBinContent(i+1)!=0 && hist_result_fit->GetBinContent(i+1)==0)
            err = 1;
	  else
            err = 0;
	  band->SetBinError(i+1, err);
	}

      gPad->SetGridx();
      gPad->SetGridy();
      band->SetFillColor(kGray+3);
      band->SetFillStyle(3001);
      band->GetYaxis()->SetTitle("Data/Fit");
      band->GetXaxis()->SetTitle("Soft drop mass of probe jet (in GeV)");
      band->GetYaxis()->CenterTitle(1);
      band->GetYaxis()->SetTitleOffset(0.35);
      band->GetYaxis()->SetTitleSize(0.11);
      band->GetXaxis()->SetTitleSize(0.12);
      band->GetYaxis()->SetLabelSize(0.07);
      band->GetXaxis()->SetLabelSize(0.1);
      band->SetMaximum(1.55);
      band->SetMinimum(0.45);
      band->GetYaxis()->SetNdivisions(10);
      band->GetYaxis()->SetTickSize(0.01);
      band->GetXaxis()->SetNdivisions(10);
      band->GetXaxis()->SetTickSize(0.03);
      band->SetMarkerStyle(1);

      TPaveText *pt = new TPaveText(0.7,0.7,0.85,0.9);
      pt->AddText("Systematic uncertainty band");

      band->Draw("E2C");
      h_ratio->Draw("PE1SAME");
      pt->Draw("same");

      TLine *line = new TLine(band->GetXaxis()->GetXmin(),1,band->GetXaxis()->GetXmax(),1);
      line->SetLineColor(1);
      //line->Draw("sames");

      can->Update();
      can->SaveAs(("plots/Final_combine_fit_"+region+ptbins+channel+".png").c_str());
      
      ///////////////////////////////////////////////////   Plot only prefit histograms   ///////////////////////////////////////////////////////////////////////////////////

      THStack *hist_prefit_stack = new THStack("hist_prefit_stack","");
      THStack *hist_postfit_stack = new THStack("hist_postfit_stack","");

      TH1D *hist_allprocs_prefit = (TH1D*)hist_procs_prefit[0]->Clone();
      for(int ip = nprocs - 1; ip > -1 ; ip--)
	{
	  hist_procs_postfit[ip]->SetFillColor(ip+6);
	  hist_procs_postfit[1]->SetFillColor(7);
	  hist_postfit_stack->Add(hist_procs_postfit[ip]);

	  hist_procs_prefit[ip]->SetFillColor(ip+6);
	  hist_procs_prefit[1]->SetFillColor(7);
	  hist_procs_prefit[ip]->SetLineStyle(1);
	  hist_prefit_stack->Add(hist_procs_prefit[ip]);
	  if(ip >0)   hist_allprocs_prefit->Add(hist_procs_prefit[ip]);
	}

      float max_stack_stack = max(hist_postfit_stack->GetMaximum(),hist_prefit_stack->GetMaximum());

      hist_postfit_stack->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)*max((float)hist_data->GetMaximum(),max_stack_stack));
      hist_postfit_stack->SetMinimum(0.001);

      hist_prefit_stack->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)*max((float)hist_data->GetMaximum(),max_stack_stack));
      hist_prefit_stack->SetMinimum(0.001);

      TGraphAsymmErrors *h_ratio_prefit = new TGraphAsymmErrors(hist_data,hist_allprocs_prefit,"pois");

      TLegend *legend_prefit = tdrLeg(0.54,0.6,0.84,0.85);
      legend_prefit->SetBorderSize(0);
      legend_prefit->SetTextFont(42);
      legend_prefit->SetTextSize(0.035);
      legend_prefit->AddEntry(hist_data, "Data", "6ple1");

      for(int ip = 0; ip < nprocs; ip++)
        legend_prefit->AddEntry(hist_procs_prefit[ip], channelname[ip].c_str(), "l");

      TCanvas * can_prefit = tdrDiCanvas("can_prefit",hist_data,hist_data,8,0);
      can_prefit->cd(1);
      gPad->SetLogy(0);

      TPaveText *mceffstats = new TPaveText(0.2,0.6,0.5,0.9,"NDCNB");
      mceffstats->AddText((string("#epsilon_{") + channelname[0] + "}^{MC} = " + to_string(eff_mc[0])).c_str());
      // mceffstats->AddText((string("#epsilon_{") + channelname[1] + "}^{MC} = " + to_string(eff_mc[1])).c_str());
      //mceffstats->AddText((string("#epsilon_{") + channelname[2] + "}^{MC} = " + to_string(eff_mc[2])).c_str());
      //      mceffstats->AddText((string("#epsilon_{") + channelname[3] + "}^{MC} = " + to_string(eff_mc[3])).c_str());
      
      //mceffstats->AddText((string("#splitline{#splitline{#splitline{#epsilon_{") + processes[0] + "}^{MC} = " + to_string(eff_mc[0]) + string("}{#epsilon_{") + processes[1] + "}^{MC} = " + to_string(eff_mc[1]) + "}}{#epsilon_{" + processes[2] + "}^{MC} = " + to_string(eff_mc[2]) + "}}{#epsilon_{") + processes[3] + "}^{MC} = " + to_string(eff_mc[3]) + string("}"))).c_str());
      mceffstats->SetTextSize(0.03);
      mceffstats->SetFillColor(0);
      
      hist_prefit_stack->Draw("hist");
      hist_prefit_stack->GetYaxis()->SetTitle(("Events / " + to_string(event_perbin) + " GeV").c_str());
      
      legend_prefit->Draw("same");
      hist_data->Draw("PE1SAME");
      if(ir > 0 ) mceffstats->Draw();
	
      can_prefit->cd(2);
      h_ratio_prefit->SetMarkerColor(1);
      h_ratio_prefit->SetMarkerStyle(20);
      h_ratio_prefit->SetMarkerSize(0.89);
      h_ratio_prefit->SetLineColor(kBlack);
      TH1D *band_prefit = new TH1D("Band_prefit", "", nbin, bins);

      for(int i=0;i<nbin; i++)
	{
	  band_prefit->SetBinContent(i+1, 1.0);
	  double err = 0;
	  if(hist_data->GetBinContent(i+1)!=0 && hist_procs_prefit_total->GetBinContent(i+1)!=0)
	    err=hist_procs_prefit_total->GetBinError(i+1); //(hist_procs_prefit_total->GetBinError(i+1))/(hist_procs_prefit_total->GetBinContent(i+1));
	  else if(hist_data->GetBinContent(i+1)!=0 && hist_procs_prefit_total->GetBinContent(i+1)==0)
            err = 1;
	  else
            err = 0;
	    cout<<err<<" ";
	  band_prefit->SetBinError(i+1, err);
		  
	}
    
      gPad->SetGridx();
      gPad->SetGridy();
      band_prefit->SetFillColor(kGray+3);
      band_prefit->SetFillStyle(3001);
      band_prefit->GetYaxis()->SetTitle("Data/MC");
      band_prefit->GetXaxis()->SetTitle("Soft drop mass of probe jet (in GeV)");
      band_prefit->GetYaxis()->CenterTitle(1);
      band_prefit->GetYaxis()->SetTitleOffset(0.35);
      band_prefit->GetYaxis()->SetTitleSize(0.11);
      band_prefit->GetXaxis()->SetTitleSize(0.12);
      band_prefit->GetYaxis()->SetLabelSize(0.07);
      band_prefit->GetXaxis()->SetLabelSize(0.1);
      band_prefit->SetMaximum(1.55);
      band_prefit->SetMinimum(0.45);
      band_prefit->GetYaxis()->SetNdivisions(10);
      band_prefit->GetYaxis()->SetTickSize(0.01);
      band_prefit->GetXaxis()->SetNdivisions(10);
      band_prefit->GetXaxis()->SetTickSize(0.03);
      band_prefit->SetMarkerStyle(1);
      
      band_prefit->Draw("E2C");
      h_ratio_prefit->Draw("PE1SAME");
      //pt->Draw("same");
      line->Draw("same");

      can_prefit->Update();
      can_prefit->SaveAs(("plots/Final_combine_fit_"+region+ptbins+channel+"_prefitonly.png").c_str());

      ///////////////////////////////////////////////////   Plot only postfit histograms   ///////////////////////////////////////////////////////////////////////////////////
      
      hist_result_fit->SetLineStyle(2);

      TLegend *legend_postfit = tdrLeg(0.54,0.6,0.84,0.85);
      legend_postfit->SetBorderSize(0);
      legend_postfit->SetTextFont(42);
      legend_postfit->SetTextSize(0.035);
      legend_postfit->AddEntry(hist_data, "Data", "6ple1");
      legend_postfit->AddEntry(hist_result_fit, "Total Fit", "l");

      for(int ip = 0; ip < nprocs; ip++)
        legend_postfit->AddEntry(hist_procs_postfit[ip], channelname[ip].c_str(), "l");

      TCanvas * can_postfit = tdrDiCanvas("can_postfit",hist_data,hist_data,8,0);
      can_postfit->cd(1);
      gPad->SetLogy(0);

      string text = string("#splitline{#splitline{#epsilon_{") + processes[0] + "}^{Data} = " + to_string(eff_data[0]) + string("}{#epsilon_{") + processes[1] + "}^{Data} = " + to_string(eff_data[1]) + "}}{#epsilon_{" + processes[2] + "}^{Data} = " + to_string(eff_data[2]) + "}";
      cout<<endl<<endl<<endl<<text;
      TPaveText *dataeffstats = new TPaveText(0.2,0.6,0.5,0.9,"NDCNB");
      dataeffstats->AddText((string("#epsilon_{") + channelname[0] + "}^{Data} = " + to_string(eff_data[0])).c_str());
      //dataeffstats->AddText((string("#epsilon_{") + channelname[1] + "}^{Data} = " + to_string(eff_data[1])).c_str());
      //dataeffstats->AddText((string("#epsilon_{") + channelname[2] + "}^{Data} = " + to_string(eff_data[2])).c_str());
      //  dataeffstats->AddText((string("#epsilon_{") + channelname[3] + "}^{Data} = " + to_string(eff_data[3])).c_str());

      dataeffstats->AddText((string("SF_{") + channelname[0] + "} = " + to_string(sf_witherror[0]) + "_{-" + to_string(sf_witherror[2]) + "}^{+" + to_string(sf_witherror[1]) + "}").c_str());

      
      //      dataeffstats->AddText((string("#splitline{#splitline{#splitline{#epsilon_{") + processes[0] + "}^{Data} = " + to_string(eff_data[0]) + string("}{#epsilon_{") + processes[1] + "}^{Data} = " + to_string(eff_data[1]) + "}}{#epsilon_{" + processes[2] + "}^{Data} = " + to_string(eff_data[2]) + "}}{#epsilon_{") + processes[3] + "}^{Data} = " + to_string(eff_data[3]) + string("}")).c_str());
      dataeffstats->SetTextSize(0.03);
      dataeffstats->SetFillColor(0);

      hist_postfit_stack->Draw("hist");
      hist_postfit_stack->GetHistogram()->GetYaxis()->SetTitle(("Events / " + to_string(event_perbin) + " GeV").c_str());

      hist_result_fit->Draw("histsame");
      legend_postfit->Draw("same");
      hist_data->Draw("PE1SAME");
      if(ir > 0) dataeffstats->Draw();
  
      can_postfit->cd(2);

      gPad->SetGridx();
      gPad->SetGridy();

      band->Draw("E2C");
      h_ratio->Draw("PE1SAME");
      pt->Draw("same");

      can_postfit->Update();
      can_postfit->SaveAs(("plots/Final_combine_fit_"+region+ptbins+channel+"_postfitonly.png").c_str());

      delete hist_result_fit_improperbin;
      delete hist_sig_fit_improperbin;
      delete hist_bkg_fit_improperbin;
      delete hist_data;

      for(int ip = nprocs - 1; ip > -1 ; ip--)
	{
	  delete hist_procs_postfit[ip];
	  delete hist_procs_postfit_improperbin[ip];
	  delete hist_procs_prefit[ip];
	}
      
      filein->Close();
      fitfile->Close();  
    }
  
  gr_sf->SetPoint(iptbins-1,iptbins,sf_witherror[0]);
  gr_sf->SetPointError(iptbins-1,0,0,sf_witherror[2],sf_witherror[1]);

  ofstream f1;
  f1.open("plots/combine_sf_yield.txt",ios::app);
  f1<<"Efficiency in MC = "<<eff_mc[0]<<" And Efficiency in Data = "<<eff_data[0]<<endl;
  f1<<"Calculated SF = "<<sf_witherror[0]<<" +- "<<sf_witherror[1]<<"(up) "<<sf_witherror[2]<<"(down)"<<endl<<endl<<endl;
  f1.close();
      
  return sf_witherror;
}


void get_theorysys_hist(string fileinname, string treename, string fillhist, string evt_weight, string event_selection_cuts, float sample_weight, int ntheorysys, TH1D *hist_out[120])
{
    TFile *filein;
    filein = new TFile(fileinname.c_str(),"read");
    if(!filein->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}

    TTree *tree = (TTree*)filein->Get(treename.c_str());
    if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}

    int N_sel_events = tree->Draw("ievt",event_selection_cuts.c_str(),"goff");
    if(N_sel_events < 1)
      return;

    TTree * tree_sel_evts = tree->CopyTree(event_selection_cuts.c_str());
    TH1D *tmp = (TH1D*)hist_out[0]->Clone("tmp");
    cout<<"check start "<<tree_sel_evts<<" entries = "<<tree_sel_evts->GetEntries()<<endl;
    if(ntheorysys > 0)    //////////////   sys for alpha_s
      {
	tree_sel_evts->Draw("nLHEAlpsWeights","nLHEAlpsWeights!=3");
	if(gPad != nullptr)
	  {
	    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
	    if(htemp != nullptr)
	      {
		cout<<"\nCheck file:"<<fileinname<<". It has some events with nLHEAlpsWeights!=3"<<endl;
	      }
	  }
	
	string weight_alphaup = evt_weight + "*LHEAlpsWeights[2]/LHEAlpsWeights[0]";
	string weight_alphadown = evt_weight + "*LHEAlpsWeights[1]/LHEAlpsWeights[0]";
	
	tmp->Reset("ICES");
	tree_sel_evts->Draw(fillhist.c_str(),weight_alphaup.c_str());
	if(tmp != NULL)
	  {
	    tmp->Scale(sample_weight);
	    hist_out[0]->Add(tmp);
	  }

	tmp->Reset("ICES");
	tree_sel_evts->Draw(fillhist.c_str(),weight_alphadown.c_str());
	if(tmp != NULL)
	  {
	    tmp->Scale(sample_weight);
	    hist_out[1]->Add(tmp);
	  }
      }
  
    if(ntheorysys > 1)    //////////////   sys for #mu_r and #mu_f scale unc.
      {
	tree_sel_evts->Draw("nLHEScaleWeights","nLHEScaleWeights!=9");
	if(gPad != nullptr)
	  {
	    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
	    if(htemp != nullptr)
	      {
		cout<<"\nCheck file:"<<fileinname<<". It has some events with nLHEScaleWeights!=9"<<endl;
	      }
	  }
	
	/*  scale ordering of muf and mur
	    facscfact 1    ! factorization scale factor: mufact=muref*facscfact
	    renscfact 1    ! renormalization scale factor: muren=muref*renscfact
	      
	    <weight id="1001"> lhapdf=306000 renscfact=1d0 facscfact=1d0 </weight>
	    <weight id="1002"> lhapdf=306000 renscfact=1d0 facscfact=2d0 </weight>
	    <weight id="1003"> lhapdf=306000 renscfact=1d0 facscfact=0.5d0 </weight>
	    <weight id="1004"> lhapdf=306000 renscfact=2d0 facscfact=1d0 </weight>
	    <weight id="1005"> lhapdf=306000 renscfact=2d0 facscfact=2d0 </weight>
	    <weight id="1006"> lhapdf=306000 renscfact=2d0 facscfact=0.5d0 </weight>
	    <weight id="1007"> lhapdf=306000 renscfact=0.5d0 facscfact=1d0 </weight>
	    <weight id="1008"> lhapdf=306000 renscfact=0.5d0 facscfact=2d0 </weight>
	    <weight id="1009"> lhapdf=306000 renscfact=0.5d0 facscfact=0.5d0 </weight>
	*/
	    
	string weights_scalesys[7];
	int ik=0;
	for(int i = 0; i<9; i++)
	  {
	    if(i != 5 && i != 7)
	      {
		weights_scalesys[ik] = evt_weight + "*LHEScaleWeights[" + to_string(i) + "]/max(0.001,LHEScaleWeights[0])";
		ik++;
	      }
	  }
	for(int ih = 0; ih<7; ih++)
	  {
	    tmp->Reset("ICES");
	    tree_sel_evts->Draw(fillhist.c_str(),weights_scalesys[ih].c_str());
	    if(tmp != NULL)
	      {
		tmp->Scale(sample_weight);
		hist_out[2+ih]->Add(tmp);
	      }
	  }
      }
  
    if(ntheorysys > 2)    //////////////   sys for Parton Shower unc.
      {
	tree_sel_evts->Draw("nLHEPSWeights","nLHEPSWeights!=8");
	if(gPad != nullptr)
	  {
	    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
	    if(htemp != nullptr)
	      {
		cout<<"\nCheck file:"<<fileinname<<". It has some events with nLHEScaleWeights!=9"<<endl;
	      }
	  }
	
	string weights_pssys[8];
	for(int ih = 0; ih<8; ih++)
	  {
	    weights_pssys[ih] = evt_weight + "*LHEPSWeights[" + to_string(ih) + "]";
	    tmp->Reset("ICES");
	    tree_sel_evts->Draw(fillhist.c_str(),weights_pssys[ih].c_str());
	    if(tmp != NULL)
	      {
		tmp->Scale(sample_weight);
		hist_out[9+ih]->Add(tmp);
		
	      }
	  }
      }     
    if(ntheorysys > 3)    //////////////   sys for Parton Shower unc.
      {
	tree_sel_evts->Draw("nLHEPDFWeights","nLHEPDFWeights!=103");
	if(gPad != nullptr)
	  {
	    TH1F *htemp = (TH1F*)gPad->GetPrimitive("htemp");
	    if(htemp != nullptr)
	      {
		cout<<"\nCheck file:"<<fileinname<<". It has some events with nLHEPDFWeights!=103"<<endl;
	      }
	  }

	string weights_pdfsys[103];
	for(int ih = 0; ih<103; ih++)
	  {
	    weights_pdfsys[ih] = evt_weight + "*LHEPDFWeights[" + to_string(ih) + "]";
	    tmp->Reset("ICES");
	    tree_sel_evts->Draw(fillhist.c_str(),weights_pdfsys[ih].c_str());
	    if(tmp != NULL)
	      {
		tmp->Scale(sample_weight);
		hist_out[17+ih]->Add(tmp);
		
	      }
	  }
	tmp->Reset("ICES");
	tree_sel_evts->Draw(fillhist.c_str(),evt_weight.c_str());
	if(tmp != NULL)
	  {
	    tmp->Scale(sample_weight);
	    hist_out[120]->Add(tmp);  
	  }
      }     
    delete tree_sel_evts;
    delete tree;
    filein->Close();
}


void make_theory_syshist(string var, string weight, string event_selection_cuts, TH1D *hist_basic, string histname, string vartitle, vector<fileinfo> allfiles, string higgscombinefilename = "")
{
  string fillhist = var + ">>+" + "tmp";

  TFile *filecombine;

  static int ntheorysys = 106;

  /*
    0 = alpha_s sys ie varying alpha_s of QCD by 0.1180 ± 0.0015 GeV^2
    1 = scale unc.
    2 = PS unc.
    3 = PDF unc.
  */

  TH1D *hist_combine_theorysys_up[ntheorysys];
  TH1D *hist_combine_theorysys_down[ntheorysys];
  
  if(higgscombinefilename != "")
    {
      filecombine = new TFile(higgscombinefilename.c_str(),"update");

      for(int ih = 0; ih < ntheorysys; ih++)
	{
	  hist_combine_theorysys_up[ih] = (TH1D*)hist_basic->Clone();
	  hist_combine_theorysys_up[ih]->Sumw2();
	  hist_combine_theorysys_down[ih] = (TH1D*)hist_basic->Clone();
	  hist_combine_theorysys_down[ih]->Sumw2();
	}
    }
  
  for(int ifile=1;ifile<allfiles.size();ifile++)
    {
      if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}
      
      TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");
      
      TH1D *hist_out[121];
      TH1D *hist_out_theorysys_up[ntheorysys], *hist_out_theorysys_down[ntheorysys];
      for(int ih = 0; ih < ntheorysys; ih++)
	{
	  hist_out_theorysys_up[ih] = (TH1D*)hist_basic->Clone();
	  hist_out_theorysys_up[ih]->Sumw2();
	  hist_out_theorysys_down[ih] = (TH1D*)hist_basic->Clone();
	  hist_out_theorysys_down[ih]->Sumw2();
	}
      for(int ih = 0; ih < 121; ih++)
	{
	  hist_out[ih] = (TH1D*)hist_basic->Clone();
	  hist_out[ih]->Sumw2();
	}
      cout<<" opening file:"<<allfiles[ifile].out_name<<endl;
      for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {
	  get_theorysys_hist(allfiles[ifile].in_names[ig],allfiles[ifile].treename,fillhist,weight,event_selection_cuts,allfiles[ifile].weights[ig],ntheorysys,hist_out);
	}

            
      if(ntheorysys > 0)
	{
	  hist_out_theorysys_up[0] = (TH1D*)hist_out[0]->Clone();;
	  hist_out_theorysys_down[0] = (TH1D*)hist_out[1]->Clone();;
	}

      int nbin = hist_basic->GetNbinsX();
      
      if(ntheorysys > 1)
	{
	  for(int ib = 1; ib<=nbin; ib++)
	    {
	      double min_content = hist_out[2]->GetBinContent(ib), max_content = hist_out[2]->GetBinContent(ib) ;
	      for(int ih = 1; ih<7; ih++)
		{
		  min_content = min(min_content,hist_out[2+ih]->GetBinContent(ib));
		  max_content = max(max_content,hist_out[2+ih]->GetBinContent(ib));
		}
	      cout<<"scale unc "<<max_content<<" min content = "<<min_content<<endl;
	      hist_out_theorysys_up[1]->SetBinContent(ib, max_content);
	      hist_out_theorysys_down[1]->SetBinContent(ib, min_content);
	    }
	}
      
      if(ntheorysys > 2)
	{
	  for(int ib = 1; ib<=nbin; ib++)
	    {
	      double min_content = hist_out[9]->GetBinContent(ib), max_content = hist_out[9]->GetBinContent(ib) ;
	      for(int ih = 2; ih<3; ih++)  ////////////////////// run over The Def(ault) setting (mu_R factor 2) is recommended for now check later for all histogram
		{
		  min_content = min(min_content,hist_out[9+ih]->GetBinContent(ib));
		  max_content = max(max_content,hist_out[9+ih]->GetBinContent(ib));
		}
	      
	      hist_out_theorysys_up[2]->SetBinContent(ib, max_content);
	      hist_out_theorysys_down[2]->SetBinContent(ib, min_content);
	    } 
	}
      
      /*      if(ntheorysys > 3)
	{
	  for(int ib = 1; ib<=nbin; ib++)
	    {
	      double pdf_unc_corr = 0, nomimal_content = hist_out[120]->GetBinContent(ib) ; ///////// CHeck which one is normial
	      for(int ih = 1; ih<103; ih++)
		{
		  pdf_unc_corr += (nomimal_content - hist_out[17+ih]->GetBinContent(ib))*(nomimal_content - hist_out[17+ih]->GetBinContent(ib));
		}
	      cout<<" PDF unc. bin "<<ib<<" nomial = "<<nomimal_content<<" pdf error = "<<pdf_unc_corr<<endl;
	      hist_out_theorysys_up[3]->SetBinContent(ib, nomimal_content + sqrt(pdf_unc_corr));
	      hist_out_theorysys_down[3]->SetBinContent(ib, nomimal_content - sqrt(pdf_unc_corr));
	    } 
	  
	}
      */
           if(ntheorysys > 4)
	{
	  for(int ib = 1; ib<=nbin; ib++)
	    {
	      double pdf_unc_corr = 0, nomimal_content = hist_out[120]->GetBinContent(ib) ; ///////// CHeck which one is normial
	      double pdf0_content = hist_out[17]->GetBinContent(ib);
	      cout<<" PDF unc. bin "<<ib<<" nomial = "<<nomimal_content<<" nomial - 1st egin weight pdf = "<<nomimal_content - pdf0_content<<endl;
	      for(int ih = 0; ih<103; ih++)
		{
		  double diff_pdf = abs(pdf0_content - hist_out[17+ih]->GetBinContent(ib));
		  hist_out_theorysys_up[3+ih]->SetBinContent(ib, nomimal_content + diff_pdf);
		  hist_out_theorysys_down[3+ih]->SetBinContent(ib, nomimal_content - diff_pdf);
		}
	    } 
	  
	}
      if(ifile != 0) {
	
	for(int ih = 0; ih < ntheorysys; ih++)
	  {
	    hist_out_theorysys_up[ih]->Scale(59.83);
	    if(higgscombinefilename != "") hist_combine_theorysys_up[ih]->Add(hist_out_theorysys_up[ih]);
	    hist_out_theorysys_down[ih]->Scale(59.83);
	    if(higgscombinefilename != "") hist_combine_theorysys_down[ih]->Add(hist_out_theorysys_down[ih]);
	  }
      }
      
      fileout->cd();

      for(int ih = 0; ih < ntheorysys; ih++)
	{
	  hist_out_theorysys_up[ih]->GetXaxis()->SetTitle(vartitle.c_str());
	  hist_out_theorysys_down[ih]->GetXaxis()->SetTitle(vartitle.c_str());
	}

      if(ntheorysys > 0)   hist_out_theorysys_up[0]->SetName((histname + "_alphaUp").c_str());
      if(ntheorysys > 1)   hist_out_theorysys_up[1]->SetName((histname + "_scaleUp").c_str());
      if(ntheorysys > 2)   hist_out_theorysys_up[2]->SetName((histname + "_partonshowerUp").c_str());
      if(ntheorysys == 4)  hist_out_theorysys_up[3]->SetName((histname + "_pdfUp").c_str());
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_out_theorysys_up[ih]->SetName((histname + "_pdf" + to_string(ih-3) + "Up").c_str());
	}
    
      if(ntheorysys > 0)   hist_out_theorysys_up[0]->Write((histname + "_alphaUp").c_str(),TObject::kOverwrite);
      if(ntheorysys > 1)   hist_out_theorysys_up[1]->Write((histname + "_scaleUp").c_str(),TObject::kOverwrite);
      if(ntheorysys > 2)   hist_out_theorysys_up[2]->Write((histname + "_partonshowerUp").c_str(),TObject::kOverwrite);
      if(ntheorysys == 4)  hist_out_theorysys_up[3]->Write((histname + "_pdfUp").c_str(),TObject::kOverwrite);
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_out_theorysys_up[ih]->Write((histname + "_pdf" + to_string(ih-3) + "Up").c_str(),TObject::kOverwrite);
	}
    
      if(ntheorysys > 0)   hist_out_theorysys_down[0]->SetName((histname + "_alphaDown").c_str());
      if(ntheorysys > 1)   hist_out_theorysys_down[1]->SetName((histname + "_scaleDown").c_str());
      if(ntheorysys > 2)   hist_out_theorysys_down[2]->SetName((histname + "_partonshowerDown").c_str());
      if(ntheorysys == 4)  hist_out_theorysys_down[3]->SetName((histname + "_pdfDown").c_str());
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_out_theorysys_down[ih]->SetName((histname + "_pdf" + to_string(ih-3) + "Down").c_str());
	}

      if(ntheorysys > 0)   hist_out_theorysys_down[0]->Write((histname + "_alphaDown").c_str(),TObject::kOverwrite);
      if(ntheorysys > 1)   hist_out_theorysys_down[1]->Write((histname + "_scaleDown").c_str(),TObject::kOverwrite);
      if(ntheorysys > 2)   hist_out_theorysys_down[2]->Write((histname + "_partonshowerDown").c_str(),TObject::kOverwrite);
      if(ntheorysys == 4)  hist_out_theorysys_down[3]->Write((histname + "_pdfDown").c_str(),TObject::kOverwrite);
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_out_theorysys_down[ih]->Write((histname + "_pdf" + to_string(ih-3) + "Down").c_str(),TObject::kOverwrite);
	}

      fileout->Close();
    }
  
  if(higgscombinefilename != "")
    {
      filecombine->cd();

      if(ntheorysys > 0)   hist_combine_theorysys_up[0]->SetName((histname + "_alphaUp").c_str());
      if(ntheorysys > 1)   hist_combine_theorysys_up[1]->SetName((histname + "_scaleUp").c_str());
      if(ntheorysys > 2)   hist_combine_theorysys_up[2]->SetName((histname + "_partonshowerUp").c_str());
      if(ntheorysys == 4)  hist_combine_theorysys_up[3]->SetName((histname + "_pdfUp").c_str());
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_combine_theorysys_up[ih]->SetName((histname + "_pdf" + to_string(ih-3) + "Up").c_str());
	}

      if(ntheorysys > 0)   hist_combine_theorysys_up[0]->Write((histname + "_alphaUp").c_str(),TObject::kOverwrite);
      if(ntheorysys > 1)   hist_combine_theorysys_up[1]->Write((histname + "_scaleUp").c_str(),TObject::kOverwrite);
      if(ntheorysys > 2)   hist_combine_theorysys_up[2]->Write((histname + "_partonshowerUp").c_str(),TObject::kOverwrite);
      if(ntheorysys == 4)  hist_combine_theorysys_up[3]->Write((histname + "_pdfUp").c_str(),TObject::kOverwrite);
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_combine_theorysys_up[ih]->Write((histname + "_pdf" + to_string(ih-3) + "Up").c_str(),TObject::kOverwrite);
	}

      if(ntheorysys > 0)   hist_combine_theorysys_down[0]->SetName((histname + "_alphaDown").c_str());
      if(ntheorysys > 1)   hist_combine_theorysys_down[1]->SetName((histname + "_scaleDown").c_str());
      if(ntheorysys > 2)   hist_combine_theorysys_down[2]->SetName((histname + "_partonshowerDown").c_str());
      if(ntheorysys == 4)  hist_combine_theorysys_down[3]->SetName((histname + "_pdfDown").c_str());
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_combine_theorysys_down[ih]->SetName((histname + "_pdf" + to_string(ih-3) + "Down").c_str());
	}

      if(ntheorysys > 0)   hist_combine_theorysys_down[0]->Write((histname + "_alphaDown").c_str(),TObject::kOverwrite);
      if(ntheorysys > 1)   hist_combine_theorysys_down[1]->Write((histname + "_scaleDown").c_str(),TObject::kOverwrite);
      if(ntheorysys > 2)   hist_combine_theorysys_down[2]->Write((histname + "_partonshowerDown").c_str(),TObject::kOverwrite);
      if(ntheorysys == 4)  hist_combine_theorysys_down[3]->Write((histname + "_pdfDown").c_str(),TObject::kOverwrite);
      if(ntheorysys > 4 )
	{
	  for(int ih = 3; ih < ntheorysys; ih++)
	    hist_combine_theorysys_down[ih]->Write((histname + "_pdf" + to_string(ih-3) + "Down").c_str(),TObject::kOverwrite);
	}

      filecombine->Close();
    }
}

void get_jerjes_sys_hist(string fileinname, string treename, string fillhist, string weight_sel, float weight, int njessys, string jes_names[njessys], TH1D *hist_out_up[njessys+1], TH1D *hist_out_down[njessys+1])
{
  TFile *filein;
  filein = new TFile(fileinname.c_str(),"read");
  if(!filein->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}

  TTree *tree = (TTree*)filein->Get(treename.c_str());
  if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}

  TH1D *tmp = (TH1D*)hist_out_up[0]->Clone("tmp");
  tmp->Reset("ICES");
     
  for(int ih = 0; ih < njessys; ih++)
    {
      string weight_sel_up = weight_sel;
      string weight_sel_down = weight_sel;
      string to_replace = "probejet_pt";
      string nameup = "probejet_pt*probejet_pt_jesup_" + jes_names[ih];
      string namedn = "probejet_pt*probejet_pt_jesdn_" + jes_names[ih];
      int index = 0;
      while ((index = weight_sel_up.find(to_replace, index)) != string::npos) {
	weight_sel_up.replace(index, to_replace.length(), nameup);
	weight_sel_down.replace(index, to_replace.length(), namedn);
	index += namedn.length();
      }
      to_replace = "tagjet_pt";
      nameup = "tagjet_pt*tagjet_pt_jesup_" + jes_names[ih];
      namedn = "tagjet_pt*tagjet_pt_jesdn_" + jes_names[ih];
      index = 0;
      while ((index = weight_sel_up.find(to_replace, index)) != string::npos) {
	weight_sel_up.replace(index, to_replace.length(), nameup);
	weight_sel_down.replace(index, to_replace.length(), namedn);
	index += namedn.length();
      }
	
      tmp->Reset("ICES");
      tree->Draw(fillhist.c_str(),weight_sel_up.c_str());
      if(tmp != NULL)
	{
	  tmp->Scale(weight);
	  hist_out_up[ih]->Add(tmp);
	}

      tmp->Reset("ICES");
      tree->Draw(fillhist.c_str(),weight_sel_down.c_str());
      if(tmp != NULL)
	{
	  tmp->Scale(weight);
	  hist_out_down[ih]->Add(tmp);
	}
   
    }

  string weight_sel_up = weight_sel;
  string weight_sel_down = weight_sel;
  string to_replace = "_pt";
  string nameup = "_pt_resoup";
  string namedn = "_pt_resodn";
  int index = 0;
  while ((index = weight_sel_up.find(to_replace, index)) != string::npos) {
    weight_sel_up.replace(index, to_replace.length(), nameup);
    weight_sel_down.replace(index, to_replace.length(), namedn);
    index += namedn.length();
  }

  tmp->Reset("ICES");
  tree->Draw(fillhist.c_str(),weight_sel_up.c_str());
  if(tmp != NULL)
    {
      tmp->Scale(weight);
      hist_out_up[njessys]->Add(tmp);
    }
    
  tmp->Reset("ICES");
  tree->Draw(fillhist.c_str(),weight_sel_down.c_str());
  if(tmp != NULL)
    {
      tmp->Scale(weight);
      hist_out_down[njessys]->Add(tmp);
    }
    
  filein->Close();
}

void make_jes_jer_syshist(string var, string weight, string event_selection_cuts, TH1D *hist_basic, string histname, string vartitle, vector<fileinfo> allfiles, string higgscombinefilename = "")
{
  string weight_sel = weight + "*(" + event_selection_cuts + ")";
  string fillhist = var + ">>+" + "tmp";

  TFile *filecombine;

  static const int njessys = 23;

  string jes_names[njessys] = {"AbsoluteStat",
			       "AbsoluteScale",
			       "AbsoluteMPFBias",
			       "FlavorQCD",
			       "Fragmentation",
			       "PileUpDataMC",
			       "PileUpPtBB",
			       "PileUpPtEC1",
			       "PileUpPtEC2",
			       "PileUpPtRef",
			       "RelativeFSR",
			       "RelativeJEREC1",
			       "RelativeJEREC2",
			       "RelativePtBB",
			       "RelativePtEC1",
			       "RelativePtEC2",
			       "RelativeBal",
			       "RelativeSample",
			       "RelativeStatEC",
			       "RelativeStatFSR",
			       "SinglePionECAL",
			       "SinglePionHCAL",
			       "TimePtEta"};
  string jer_name[2] = {"_pt_resoup","_pt_resodn"};
    
  TH1D *hist_combine_jessys_up[njessys];
  TH1D *hist_combine_jessys_down[njessys];
  TH1D *hist_combine_jersys_up;
  TH1D *hist_combine_jersys_down;
  
  if(higgscombinefilename != "")
    {
      filecombine = new TFile(higgscombinefilename.c_str(),"update");

      for(int ih = 0; ih < njessys; ih++)
	{
	  hist_combine_jessys_up[ih] = (TH1D*)hist_basic->Clone();
	  hist_combine_jessys_up[ih]->Sumw2();
	  hist_combine_jessys_down[ih] = (TH1D*)hist_basic->Clone();
	  hist_combine_jessys_down[ih]->Sumw2();
	}
      hist_combine_jersys_up = (TH1D*)hist_basic->Clone();
      hist_combine_jersys_up->Sumw2();
      hist_combine_jersys_down = (TH1D*)hist_basic->Clone();
      hist_combine_jersys_down->Sumw2();

    }
  
  for(int ifile=1;ifile<allfiles.size();ifile++)
    {
      if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}
      
      TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");
      
      TH1D *hist_out_up[njessys+1], *hist_out_down[njessys+1];

      for(int ih = 0; ih < njessys+1; ih++)
	{
	  hist_out_up[ih] = (TH1D*)hist_basic->Clone();
	  hist_out_up[ih]->Sumw2();
	  hist_out_down[ih] = (TH1D*)hist_basic->Clone();
	  hist_out_down[ih]->Sumw2();
	}
      cout<<" opening file:"<<allfiles[ifile].out_name<<endl;
      for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {
	  get_jerjes_sys_hist(allfiles[ifile].in_names[ig],allfiles[ifile].treename,fillhist,weight_sel,allfiles[ifile].weights[ig],njessys,jes_names,hist_out_up,hist_out_down);
	}
         
      if(ifile != 0) {
	
	for(int ih = 0; ih < njessys; ih++)
	  {
	    hist_out_up[ih]->Scale(59.83);
	    if(higgscombinefilename != "") hist_combine_jessys_up[ih]->Add(hist_out_up[ih]);
	    hist_out_down[ih]->Scale(59.83);
	    if(higgscombinefilename != "") hist_combine_jessys_down[ih]->Add(hist_out_down[ih]);
	  }
	hist_out_up[njessys]->Scale(59.83);
	if(higgscombinefilename != "") hist_combine_jersys_up->Add(hist_out_up[njessys]);
	hist_out_down[njessys]->Scale(59.83);
	if(higgscombinefilename != "") hist_combine_jersys_down->Add(hist_out_down[njessys]);
      }
      
      fileout->cd();

      for(int ih = 0; ih < njessys; ih++)
	{
	  hist_out_up[ih]->GetXaxis()->SetTitle(vartitle.c_str());
	  hist_out_down[ih]->GetXaxis()->SetTitle(vartitle.c_str());
	  hist_out_up[ih]->SetName((histname + "_JES_" + jes_names[ih] + "Up").c_str());
	  hist_out_up[ih]->Write((histname + "_JES_" + jes_names[ih] + "Up").c_str(),TObject::kOverwrite);
	  hist_out_down[ih]->SetName((histname + "_JES_" + jes_names[ih] + "Down").c_str());
	  hist_out_down[ih]->Write((histname + "_JES_" + jes_names[ih] + "Down").c_str(),TObject::kOverwrite);
	}
      
      hist_out_up[njessys]->GetXaxis()->SetTitle(vartitle.c_str());
      hist_out_down[njessys]->GetXaxis()->SetTitle(vartitle.c_str());
      hist_out_up[njessys]->SetName((histname + "_JER" + "Up").c_str());
      hist_out_up[njessys]->Write((histname + "_JER" + "Up").c_str(),TObject::kOverwrite);
      hist_out_down[njessys]->SetName((histname + "_JER" + "Down").c_str());
      hist_out_down[njessys]->Write((histname + "_JER" + "Down").c_str(),TObject::kOverwrite);
      
      fileout->Close();
    }
  
  if(higgscombinefilename != "")
    {
      filecombine->cd();

      for(int ih = 0; ih < njessys; ih++)
	{
	  hist_combine_jessys_up[ih]->GetXaxis()->SetTitle(vartitle.c_str());
	  hist_combine_jessys_down[ih]->GetXaxis()->SetTitle(vartitle.c_str());
	  hist_combine_jessys_up[ih]->SetName((histname + "_JES_" + jes_names[ih] + "Up").c_str());
	  hist_combine_jessys_up[ih]->Write((histname + "_JES_" + jes_names[ih] + "Up").c_str(),TObject::kOverwrite);
	  hist_combine_jessys_down[ih]->SetName((histname + "_JES_" + jes_names[ih] + "Down").c_str());
	  hist_combine_jessys_down[ih]->Write((histname + "_JES_" + jes_names[ih] + "Down").c_str(),TObject::kOverwrite);
	}
      
      hist_combine_jersys_up->GetXaxis()->SetTitle(vartitle.c_str());
      hist_combine_jersys_down->GetXaxis()->SetTitle(vartitle.c_str());
      hist_combine_jersys_up->SetName((histname + "_JER" + "Up").c_str());
      hist_combine_jersys_up->Write((histname + "_JER" + "Up").c_str(),TObject::kOverwrite);
      hist_combine_jersys_down->SetName((histname + "_JER" + "Down").c_str());
      hist_combine_jersys_down->Write((histname + "_JER" + "Down").c_str(),TObject::kOverwrite);

      filecombine->Close();
    }
}


/*
void check_sysematic(int nsys_err, string sys_names, string histname, string filename)
{
  TFile *filein;
  filein = new TFile(filename.c_str(),"read");
  if(!filein->IsOpen())
    {cout<<"Check: Input file  to higgs combine  with name "<<filename<<" is not present which is needed to get the data histogram (cant taken from combine output as the x asis range will be not correct). Make sure its name is in format combine_*channel*.root"<<endl;exit(0);}
  
  filein->cd();
  
  string sys_names[] = {"pileup",
			"btag",
			"hadtagger",
			"JER",
			  "alpha",
			  "pdf",
			  "partonshower",
			  "JES_AbsoluteStat",
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
			  "JES_TimePtEta"};
      int nsys_err = sizeof(sys_names)/sizeof(sys_names[0]); ;
      
      TH1D *hist_procs_prefit[nprocs];
      TH1D *hist_procs_prefit_sys[nprocs][nsys_err][2];
      TH1D *hist_procs_prefit_relsyserror[nprocs][nsys_err];
      TH1D *hist_procs_prefit_total;
      for(int ip = 0; ip < nprocs; ip++)
	{
	  string prochistname = "sdmass_probejet_" + region + "_" + processes[ip] + "_" + channel;
	  hist_procs_prefit[ip] = (TH1D*)filein->Get(prochistname.c_str());
	  cout<<" check mc name = "<<hist_procs_prefit[ip]->GetName()<<" "<<hist_procs_prefit[ip]->GetTitle()<<endl;
	  if(hist_procs_prefit[ip] == NULL)
	    {cout<<"Check: file "<<datafilename<<" does not have histogram for process  "<<processes[ip]<<" with hist name:"<<prochistname<<endl;exit(0);}
	  for(int isys = 0; isys < nsys_err; isys++)
	    {
	      string prochistname_sys = "sdmass_probejet_" + region + "_" + processes[ip] + "_" + channel + "_" + sys_names[isys];
	      hist_procs_prefit_sys[ip][isys][0] = (TH1D*)filein->Get((prochistname_sys + "Up").c_str());
	      hist_procs_prefit_sys[ip][isys][1] = (TH1D*)filein->Get((prochistname_sys + "Down").c_str());
	      hist_procs_prefit_relsyserror[ip][isys] = (TH1D*)hist_procs_prefit[ip]->Clone();
	    }
	  if(ip==0)
	    hist_procs_prefit_total = (TH1D*)hist_procs_prefit[ip]->Clone();
	  //  else
	  //  hist_procs_prefit_total->Add(hist_procs_prefit[ip]);
	}
      
      for(int ib = 1; ib < hist_procs_prefit_total->GetNbinsX() +1; ib++)
	{
	  double err = 0;
	  for(int ip = 0; ip < nprocs; ip++)
	    {
	      double nom_yield = hist_procs_prefit[ip]->GetBinContent(ib); 
	      for(int isys = 0; isys < nsys_err; isys++)
		{
		  double sys_error = (hist_procs_prefit_sys[ip][isys][0]->GetBinContent(ib) - hist_procs_prefit_sys[ip][isys][1]->GetBinContent(ib))/nom_yield;
		  err +=  sys_error*sys_error;
		  hist_procs_prefit_relsyserror[ip][isys]->SetBinContent(ib,sys_error);
		  hist_procs_prefit_relsyserror[ip][isys]->SetBinError(ib,0);
		  
		}
	    }
	  cout<<"bin = "<<ib<<" content = "<<hist_procs_prefit_total->GetBinContent(ib)<<" error = "<<err;
	  hist_procs_prefit_total->SetBinError(ib,err);
	}

      string dataname[3] = {"topmatch","blepmatch","nomatch"};
      for(int isys = 0; isys < nsys_err; isys++)
	{
	  TH1D *giv_hist[3];
	  giv_hist[0] = hist_procs_prefit_relsyserror[0][isys];
	  giv_hist[1] = hist_procs_prefit_relsyserror[1][isys];
	  giv_hist[2] = hist_procs_prefit_relsyserror[2][isys];
	  plot_1dhists(3,giv_hist,dataname,("../sysematic/" + histname + "_relerr_" + sys_names[isys] + ".png"));
	}

      plot_1dhists(1,hist_procs_prefit_total,"total MC",("../sysematic/" + histname + "_allerr.png"));
	
}*/



void compareshape(int nhist, TH1D *hist_obs[nhist], string dataname[nhist],string plotname)
{
  
  TCanvas *cv;
  TLegend *legv;
  
  cv = tdrCanvas("canv_d",hist_obs[0],8,0);
  cv->cd();
  legv = tdrLeg(0.2,0.65,0.875,0.915);

  legv->SetTextFont(42);
  legv->SetTextSize(0.045);
  legv->SetBorderSize(0);
  
  for(int ih = 0; ih < nhist; ih++)
    {
      hist_obs[ih]->Scale(1.0/max(0.00001,hist_obs[ih]->Integral()));
      hist_obs[ih]->SetFillStyle(0);
      hist_obs[ih]->SetFillColor(0);
      hist_obs[ih]->SetLineColor(ih+2);
      //     hist_obs[ih]->SetLineStyle(ih+2);
      hist_obs[ih]->SetLineWidth(2);
      hist_obs[ih]->GetYaxis()->SetTitleOffset(1.5);

      legv->AddEntry(hist_obs[ih],dataname[ih].c_str(),"lep");
   
      hist_obs[ih]->SetMaximum(1);
      hist_obs[ih]->SetMinimum(0.01);
    }
    //gPad->SetLogy(1);

  hist_obs[0]->Draw("hist");
  for(int ih = 1; ih < nhist; ih++)
    hist_obs[ih]->Draw("histsame");

  legv->Draw("same");
  CMS_lumi( cv, 8, 0 );

  cv->SaveAs(("plots/" + plotname).c_str());
}


#endif















 	  /*
	  TH2D *corelation_matrix = (TH2D*)fitfile->Get("covariance_fit_s");

	  TTree * tree_fit = (TTree*)fitfile->Get("tree_fit_sb");
	  double failbin1,failbin2,failbin3,failbin4,failbin5,failbin6,failbin7,failbin8,failbin9,r,passbin1,passbin2,passbin3,passbin4,passbin5,passbin6,passbin7,passbin8,passbin9,btag,hadtagger,pileup;

	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin1").c_str(),&failbin1);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin2").c_str(),&failbin2);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin3").c_str(),&failbin3);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin4").c_str(),&failbin4);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin5").c_str(),&failbin5);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin6").c_str(),&failbin6);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin7").c_str(),&failbin7);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin8").c_str(),&failbin8);
	  tree_fit->SetBranchAddress(("prop_binfail" + channel + "_bin9").c_str(),&failbin9);
	  tree_fit->SetBranchAddress("r",&r);
	  tree_fit->SetBranchAddress("btag",&btag);
	  tree_fit->SetBranchAddress("hadtagger",&hadtagger);
	  tree_fit->SetBranchAddress("pileup",&pileup);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin1").c_str(),&passbin1);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin2").c_str(),&passbin2);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin3").c_str(),&passbin3);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin4").c_str(),&passbin4);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin5").c_str(),&passbin5);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin6").c_str(),&passbin6);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin7").c_str(),&passbin7);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin8").c_str(),&passbin8);
	  tree_fit->SetBranchAddress(("prop_binpass" + channel + "_bin9").c_str(),&passbin9);
      
	  double fail_values_array[13],pass_values_array[13];
	  tree_fit->GetEntry(0);

	  fail_values_array[0] = btag;
	  fail_values_array[1] = hadtagger;
	  fail_values_array[2] = pileup;
	  fail_values_array[3] = failbin1;
	  fail_values_array[4] = failbin2;
	  fail_values_array[5] = failbin3;
	  fail_values_array[6] = failbin4;
	  fail_values_array[7] = failbin5;
	  fail_values_array[8] = failbin6;
	  fail_values_array[9] = failbin7;
	  fail_values_array[10] = failbin8;
	  fail_values_array[11] = failbin9;
	  fail_values_array[12] = r;
      
	  pass_values_array[0] = btag;
	  pass_values_array[1] = hadtagger;
	  pass_values_array[2] = pileup;
	  pass_values_array[3] = passbin1;
	  pass_values_array[4] = passbin2;
	  pass_values_array[5] = passbin3;
	  pass_values_array[6] = passbin4;
	  pass_values_array[7] = passbin5;
	  pass_values_array[8] = passbin6;
	  pass_values_array[9] = passbin7;
	  pass_values_array[10] = passbin8;
	  pass_values_array[11] = passbin9;
	  pass_values_array[12] = r;

	  int n_np = 3, n_bins = 10;

	  double fail_comatrix[13][13],pass_comatrix[13][13];
	  for(int i = 0; i < n_np + n_bins -1; i++)
	    {
	      for(int j = 0; j < n_np + n_bins - 1; j++)
		{
		  fail_comatrix[i][j] = corelation_matrix->GetBinContent(i+1,j+1);
		}
	    }
	  for(int j = 0; j < n_np + n_bins; j++)
	    {
	      fail_comatrix[n_np + n_bins - 1][j] = corelation_matrix->GetBinContent(n_np + 2*n_bins-1,j+1);
	      fail_comatrix[j][n_np + n_bins - 1] = corelation_matrix->GetBinContent(j+1,n_np + 2*n_bins-1);
	    }

	  for(int i = 0; i < n_np + n_bins -1; i++)
	    {
	      for(int j = 0; j < n_np + n_bins - 1; j++)
		{
		  if(j < n_np && i >= n_np)
		    pass_comatrix[i][j] = corelation_matrix->GetBinContent(i+n_bins,j+1);
		  else if(j >= n_np && i < n_np)
		    pass_comatrix[i][j] = corelation_matrix->GetBinContent(i+1,j+n_bins);
		  else if(j >= n_np && i >= n_np)
		    pass_comatrix[i][j] = corelation_matrix->GetBinContent(i+n_bins,j+n_bins);
		  else
		    pass_comatrix[i][j] = corelation_matrix->GetBinContent(i+1,j+1);
		}
	    }
	  for(int j = 0; j < n_np + n_bins; j++)
	    {
	      if(j < n_np)
		{
		  pass_comatrix[n_np + n_bins - 1][j] = corelation_matrix->GetBinContent(n_np + 2*n_bins-1,j+1);
		  pass_comatrix[j][n_np + n_bins - 1] = corelation_matrix->GetBinContent(j+1,n_np + 2*n_bins-1);
		}
	      else
		{
		  pass_comatrix[n_np + n_bins - 1][j] = corelation_matrix->GetBinContent(n_np + 2*n_bins-1,n_bins + j + 1);
		  pass_comatrix[j][n_np + n_bins - 1] = corelation_matrix->GetBinContent(n_bins + j+1,n_np + 2*n_bins-1);
		}
	    }

	  double prod_matrix[13];
	  for(int i = 0; i <n_np + n_bins; i++)
	    {
	      prod_matrix[i] =0;
	      for(int j = 0; j <n_np + n_bins; j++)
		{
		  prod_matrix[i] += pass_comatrix[j][i] *  pass_values_array[j];
		}	      
	    }

	  float pass_unc = 0;
	  for(int j = 0; j <n_np + n_bins; j++)
	    {
	      pass_unc += prod_matrix[j] *  pass_values_array[j];
	    }	      

	  for(int i = 0; i <n_np + n_bins; i++)
	    {
	      prod_matrix[i] =0;
	      for(int j = 0; j <n_np + n_bins; j++)
		{
		  prod_matrix[i] += fail_comatrix[j][i] *  fail_values_array[j];
		}	      
	    }

	  float fail_unc = 0;
	  for(int j = 0; j <n_np + n_bins; j++)
	    {
	      fail_unc += prod_matrix[j] *  fail_values_array[j];
	    }	      

	  float eff_data_up = (postfit_passsig_evts[0] * (1+pass_unc) ) / max(float(0.00001),(postfit_passsig_evts[0]*(1+pass_unc) + postfit_failsig_evts[0]*(1-fail_unc)));
	  float eff_data_down = (postfit_passsig_evts[0] * (1-pass_unc) ) / max(float(0.00001),(postfit_passsig_evts[0]*(1-pass_unc) + postfit_failsig_evts[0]*(1+fail_unc)));
	  
	  sf_witherror[1] = sf_witherror[0]*sqrt(1/max(float(0.00001),prefit_passsig_evts[0]) + 1/(max(float(0.00001),prefit_passsig_evts[0] + prefit_failsig_evts[0])) -2*sqrt(eff_mc[0]) + (eff_data[0] - eff_data_up)*(eff_data[0] - eff_data_up)/(eff_data[0]*eff_data[0]));
	  sf_witherror[2] = sf_witherror[0]*sqrt(1/max(float(0.00001),prefit_passsig_evts[0]) + 1/(max(float(0.00001),prefit_passsig_evts[0] + prefit_failsig_evts[0])) -2*sqrt(eff_mc[0]) + (eff_data[0] - eff_data_down)*(eff_data[0] - eff_data_down)/(eff_data[0]*eff_data[0]));
	  cout<<" pass event unc. = "<<pass_unc<<" fail event unc. = "<<fail_unc;
	  
	}
      else
      {*/
