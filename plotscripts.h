#ifndef PLOTSCRIPTS_H
#define PLOTSCRIPTS_H

#include "filesystem.h"
#include "My_Style.C"
////////////////////////////////////////////////////////////////////////////////////////////////////

void gethist(string fileinname, string treename, string fillhist, string weight_sel, float weight, TH1D *hist_out)
{
    TFile *filein;
    filein = new TFile(fileinname.c_str(),"read");
    if(!filein->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}
    TTree *tree = (TTree*)filein->Get(treename.c_str());
    if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}
    TH1D *tmp = (TH1D*)hist_out->Clone("tmp");
    tmp->Reset("ICES");
    tree->Draw(fillhist.c_str(),weight_sel.c_str());
        //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_basic->GetNbinsX()<<"  "<<hist_basic->GetEntries();
    //TH1D *tmp = (TH1D*)gDirectory->Get(hist_out->GetName());
    //for(int i =1; i<tmp->GetNbinsX()+1; i++)
    //cout<<"bin "<<i<<" "<<tmp->GetBinLowEdge(i)<<" to "<<tmp->GetBinLowEdge(i)+tmp->GetBinWidth(i)<<" "<<tmp->GetBinContent(i)<<endl;
    if(tmp != NULL)
    {
        //tmp->Sumw2();
        tmp->Scale(weight);
        hist_out->Add(tmp);
        //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_out->GetName()<<" "<<tmp->GetNbinsX()<<"  "<<tmp->GetEntries();
    }
    filein->Close();
}


 void makerootfile(string var, string weight, string event_selection_cuts, TH1D *hist_basic, string histname, string vartitle, vector<fileinfo> allfiles, string higgscombinefilename = "")
{
    string weight_sel = weight + "*(" + event_selection_cuts + ")";
    //cout<<hist_basic->GetName()<<" , "<<hist_basic->GetNbinsX()<<endl;
    //string fillhist = var + ">>+" + "tmp" + "(" + to_string(hist_basic->GetNbinsX()) + "," + to_string(hist_basic->GetBinLowEdge(1)) + "," + to_string(hist_basic->GetBinLowEdge(hist_basic->GetNbinsX()) + hist_basic->GetBinWidth(hist_basic->GetNbinsX())) + ")";

    string fillhist = var + ">>+" + "tmp";
  //    for(int i =1; i<tmp->GetNbinsX()+1; i++)
    //cout<<"bin "<<i<<" "<<tmp->GetBinLowEdge(i)<<" to "<<tmp->GetBinLowEdge(i)+tmp->GetBinWidth(i)<<endl;
    TFile *filecombine;
    TH1D *hist_combine;
    if(higgscombinefilename != "")
    {
        filecombine = new TFile(higgscombinefilename.c_str(),"update");
        hist_combine = (TH1D*)hist_basic->Clone();
        hist_combine->Sumw2();
    }

    for(int ifile=0;ifile<allfiles.size();ifile++)
    {
        if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}

        TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");

        TH1D *hist_out = (TH1D*)hist_basic->Clone();
        hist_out->Sumw2();
        for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {
            gethist(allfiles[ifile].in_names[ig],allfiles[ifile].treename,fillhist,weight_sel,allfiles[ifile].weights[ig],hist_out);
            //hist_basic->Scale(allfiles[ifile].weights[ig]);
            //cout<<" = "<<hist_out->GetEntries()<<endl;
            //hist_out->Add(hist_basic);
        }
        if(ifile != 0) {
            hist_out->Scale(59.83);
            if(higgscombinefilename != "") hist_combine->Add(hist_out);
        }
        fileout->cd();
        hist_out->GetXaxis()->SetTitle(vartitle.c_str());
        hist_out->Write(histname.c_str(),TObject::kOverwrite);
        fileout->Close();
    }
    if(higgscombinefilename != "")
    {
        filecombine->cd();
        hist_combine->GetXaxis()->SetTitle(vartitle.c_str());
        hist_combine->Write(histname.c_str(),TObject::kOverwrite);
        filecombine->Close();
    }
}


 void make_weightsyshist(string var, string weight, string sysweight, string event_selection_cuts, TH1D *hist_basic, string histname, string vartitle, vector<fileinfo> allfiles, string higgscombinefilename = "")
{
    string weight_sel_up = weight + "*" + sysweight + "_up" + "*(" + event_selection_cuts + ")/" + sysweight;
    string weight_sel_down = weight + "*" + sysweight + "_down" + "*(" + event_selection_cuts + ")/" + sysweight;

    //string fillhist = var + ">>+" + hist_basic->GetName() + "(" + to_string(hist_basic->GetNbinsX()) + "," + to_string(hist_basic->GetBinLowEdge(1)) + "," + to_string(hist_basic->GetBinLowEdge(hist_basic->GetNbinsX()) + hist_basic->GetBinWidth(hist_basic->GetNbinsX())) + ")";

    string fillhist = var + ">>+" + "tmp";

    TFile *filecombine;
    TH1D *hist_combineUp;
    TH1D *hist_combineDown;

    if(higgscombinefilename != ""){
        filecombine = new TFile(higgscombinefilename.c_str(),"update");
        hist_combineUp = (TH1D*)hist_basic->Clone();
        hist_combineDown = (TH1D*)hist_basic->Clone();
        hist_combineUp->SetTitle((string(hist_combineUp->GetTitle())+"Up").c_str());
        hist_combineDown->SetTitle((string(hist_combineDown->GetTitle())+"Down").c_str());
        hist_combineUp->Sumw2();
        hist_combineDown->Sumw2();
    }

    for(int ifile=0;ifile<allfiles.size();ifile++)
    {
        if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}

        TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");

        TH1D *hist_out_up = (TH1D*)hist_basic->Clone();
        TH1D *hist_out_down = (TH1D*)hist_basic->Clone();
        hist_out_up->Sumw2();
        hist_out_down->Sumw2();

        for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {
            gethist(allfiles[ifile].in_names[ig],allfiles[ifile].treename,fillhist,weight_sel_up,allfiles[ifile].weights[ig],hist_out_up);
            gethist(allfiles[ifile].in_names[ig],allfiles[ifile].treename,fillhist,weight_sel_down,allfiles[ifile].weights[ig],hist_out_down);
            //hist_basic->Scale(allfiles[ifile].weights[ig]);
            //cout<<" = "<<hist_out->GetEntries()<<endl;
            //hist_out->Add(hist_basic);
        }
        if(ifile != 0)
        {
            hist_out_up->Scale(59.83);
            hist_out_down->Scale(59.83);
            if(higgscombinefilename != "")
            {
                hist_combineUp->Add(hist_out_up);
                hist_combineDown->Add(hist_out_down);
            }
        }
        fileout->cd();
        hist_out_up->GetXaxis()->SetTitle(vartitle.c_str());
        hist_out_up->Write((histname + "Up").c_str(),TObject::kOverwrite);
        hist_out_down->GetXaxis()->SetTitle(vartitle.c_str());
        hist_out_down->Write((histname + "Down").c_str(),TObject::kOverwrite);
        fileout->Close();
    }
    if(higgscombinefilename != ""){
        filecombine->cd();
        hist_combineUp->GetXaxis()->SetTitle(vartitle.c_str());
        hist_combineUp->Write((histname + "Up").c_str(),TObject::kOverwrite);
        hist_combineDown->GetXaxis()->SetTitle(vartitle.c_str());
        hist_combineDown->Write((histname + "Down").c_str(),TObject::kOverwrite);
        filecombine->Close();
    }
}


void addhisttofile(string outfile, string outhist, string infile, string inhist)
{
    TFile *fin = new TFile(infile.c_str(),"read");
    if(!fin->IsOpen())
    {cout<<"Check: file "<<infile<<" not present"<<endl;exit(0);}
    TH1D *hist = (TH1D*)fin->Get(inhist.c_str());
    if(hist == NULL)
    {cout<<"Check: file "<<infile<<" does not have Tobject "<<inhist<<endl;exit(0);}
    TFile *fout = new TFile(outfile.c_str(),"update");
    fout->cd();
    hist->SetTitle(outhist.c_str());
    hist->SetName(outhist.c_str());
    hist->Write(outhist.c_str(),TObject::kOverwrite);
    fout->Close();
    delete hist;
    fin->Close();

}



















void get2dhist(string fileinname, string treename, string fillhist, string weight_sel, float weight, TH2D *hist_out)
{
    TFile *filein;
    filein = new TFile(fileinname.c_str(),"read");
    if(!filein->IsOpen())
    {cout<<"Check: file "<<fileinname<<" not present"<<endl;exit(0);}
    TTree *tree = (TTree*)filein->Get(treename.c_str());
    if(tree == NULL)
    {cout<<"Check: file "<<fileinname<<" does not have tree "<<treename<<endl;exit(0);}
    //tree->Draw(fillhist.c_str(),weight_sel.c_str());
        //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_basic->GetNbinsX()<<"  "<<hist_basic->GetEntries();
    //TH2D *tmp = (TH2D*)gDirectory->Get(hist_out->GetName());
    TH1D *tmp = (TH1D*)hist_out->Clone("tmp");
    tmp->Reset("ICES");
    tree->Draw(fillhist.c_str(),weight_sel.c_str());
    if(tmp != NULL)
    {
        tmp->Scale(weight);
        hist_out->Add(tmp);
        //cout<<fileinname<<" "<<treename<<" "<<fillhist<<" , "<<hist_out->GetName()<<" "<<tmp->GetNbinsX()<<"  "<<tmp->GetEntries();
    }
    filein->Close();
}


 void makerootfile(string var, string weight, string event_selection_cuts, TH2D *hist_basic, string histname, string vartitle1, string vartitle2, vector<fileinfo> allfiles, string higgscombinefilename = "")
{
    string weight_sel = weight + "*(" + event_selection_cuts + ")";
    //cout<<hist_basic->GetName()<<" , "<<hist_basic->GetNbinsX()<<endl;
    //string fillhist = var + ">>+" + hist_basic->GetName() + "(" + to_string(hist_basic->GetNbinsX()) + "," + to_string(hist_basic->GetBinLowEdge(1)) + "," + to_string(hist_basic->GetBinLowEdge(hist_basic->GetNbinsX()) + hist_basic->GetBinWidth(hist_basic->GetNbinsX())) + ")";

    string fillhist = var + ">>+" + "tmp";

    TFile *filecombine;
    TH2D *hist_combine;
    if(higgscombinefilename != "")
    {
        filecombine = new TFile(higgscombinefilename.c_str(),"update");
        hist_combine = (TH2D*)hist_basic->Clone();
        hist_combine->Sumw2();
    }
    for(int ifile=0;ifile<allfiles.size();ifile++)
    {
        if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}
  cout<<"check 1: "<<allfiles[ifile].out_name<<endl;
        TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");

        TH2D *hist_out = (TH2D*)hist_basic->Clone();
        for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {
            get2dhist(allfiles[ifile].in_names[ig],allfiles[ifile].treename,fillhist,weight_sel,allfiles[ifile].weights[ig],hist_out);
            //hist_basic->Scale(allfiles[ifile].weights[ig]);
            //cout<<" = "<<hist_out->GetEntries()<<endl;
            //hist_out->Add(hist_basic);
        }
        if(ifile != 0) {
            hist_out->Scale(59.83);
            if(higgscombinefilename != "") hist_combine->Add(hist_out);
        }
        fileout->cd();
        hist_out->GetXaxis()->SetTitle(vartitle1.c_str());
        hist_out->GetYaxis()->SetTitle(vartitle2.c_str());
        hist_out->Write(histname.c_str(),TObject::kOverwrite);
        fileout->Write();
        fileout->Close();
    }
    if(higgscombinefilename != ""){
        filecombine->cd();
        hist_combine->GetXaxis()->SetTitle(vartitle1.c_str());
        hist_combine->GetYaxis()->SetTitle(vartitle2.c_str());
        hist_combine->Write(histname.c_str(),TObject::kOverwrite);
        filecombine->Close();
    }
}

void savehist2d(string filename, string histname, string plotname)
{
  TFile *file = new TFile(filename.c_str(),"read");
  if(!file->IsOpen())
  {cout<<"Check: file "<<filename<<" not present"<<endl;exit(0);}

  TH2D *hist = (TH2D*)file->Get(histname.c_str());
  if(hist == NULL)
  {cout<<"Check: file "<<filename<<" does not have TH2D hist "<<histname<<endl;exit(0);}

  TCanvas *c1 = tdrCanvas("canv_d", hist->ProjectionX(),8,0);
  c1->cd();
  gStyle->SetOptTitle(0);

  gStyle->SetPaintTextFormat( "1.3f" );
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  //  c1->Update();
  hist->GetXaxis()->SetNdivisions(215);
  hist->GetYaxis()->SetNdivisions(215);
  hist->GetXaxis()->SetTickSize(0.03);
  hist->GetYaxis()->SetTickSize(0.03);
  hist->GetXaxis()->SetMoreLogLabels();
  hist->GetYaxis()->SetMoreLogLabels();

  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.06);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetLabelOffset(0.001);

  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.13);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->SetMarkerSize(1.4);
  //hist->GetYaxis()->SetLabelOffset(1.06);

  hist->Draw("coltext, e");
  CMS_lumi( c1, 8, 0 );

  c1->SaveAs(("./plots/" + plotname + ".png").c_str());
  file->Close();
}

void savehist2d(TH2D *hist, string plotname, bool plotwitherror = true)
{
  TCanvas *c1 = tdrCanvas("canv_d", hist->ProjectionX(),8,0);

  c1->cd();
  gStyle->SetOptTitle(0);

  gStyle->SetPaintTextFormat( "1.2f" );
  gPad->SetLogx(0);
  gPad->SetLogy(0);

  //  c1->Update();
  hist->GetXaxis()->SetNdivisions(215);
  hist->GetYaxis()->SetNdivisions(215);
  hist->GetXaxis()->SetTickSize(0.03);
  hist->GetYaxis()->SetTickSize(0.03);
  hist->GetXaxis()->SetMoreLogLabels();
  hist->GetYaxis()->SetMoreLogLabels();

  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(1.06);
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetXaxis()->SetLabelOffset(0.001);

  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(1.13);
  hist->GetYaxis()->SetLabelSize(0.04);
  hist->SetMarkerSize(1.4);
  //hist->GetYaxis()->SetLabelOffset(1.06);
  if(plotwitherror) hist->Draw("coltext, e");
  else hist->Draw("coltext");
  CMS_lumi( c1, 8, 0 );

  c1->SaveAs(("./plots/" + plotname + ".png").c_str());
}

void convertteff2hist2d(TEfficiency *Teff, TH2D* hist_eff)
{
    TH2D *h = new TH2D();
    h = (TH2D*)Teff->CreateHistogram();

    int nxbin = h->GetNbinsX();
    int nybin = h->GetNbinsY();

    for(int ij=1; ij<=nxbin; ij++)
    {
        for(int ik=1; ik<=nybin; ik++)
        {
            h->SetBinError(ij,ik,(Teff->GetEfficiencyErrorLow(Teff->GetGlobalBin(ij,ik)) + Teff->GetEfficiencyErrorUp(Teff->GetGlobalBin(ij,ik)))/2);
        }
    }
//    hist_eff = (TH2D*)
    h->Copy(*hist_eff);
}

void difftriggsf(string datafilename,string mcfilename,string historthonames2d,string  histtrignames2d,string historthonames2d_sys,string histtrignames2d_sys)
{
  TFile *file = new TFile("triggersf.root","update");
  TH2D *hist_sf_nom;
  hist_sf_nom = (TH2D*)file->Get(("sf_" + historthonames2d).c_str());
  if(hist_sf_nom == NULL)
    {cout<<"Check: file triggersf.root does not have TH2D hist "<<("sf_" + historthonames2d).c_str()<<endl;exit(0);}

  TH2D *hist_sf_sys;
  hist_sf_sys = (TH2D*)file->Get(("sf_" + historthonames2d_sys).c_str());
  if(hist_sf_sys == NULL)
    {cout<<"Check: file triggersf.root does not have TH2D hist "<<("sf_" + historthonames2d_sys).c_str()<<endl;exit(0);}

  TH2D *hist_diffsf = (TH2D*)hist_sf_nom->Clone();
  TH2D *hist_diffsf_percent = (TH2D*)hist_sf_nom->Clone();

  int nxbin = hist_diffsf->GetNbinsX();
  int nybin = hist_diffsf->GetNbinsY();

  for(int ij=1; ij<=nxbin; ij++)
	{
	  for(int ik=1; ik<=nybin; ik++)
	    {
	      double norm = hist_sf_nom->GetBinContent(ij,ik);
	      double sys_sf = hist_sf_sys->GetBinContent(ij,ik);
	      double error = sys_sf - norm;
	      hist_diffsf->SetBinContent(ij,ik,error);
          if(norm > 0.0)hist_diffsf_percent->SetBinContent(ij,ik,error/norm * 100);
	    }
	}

  string plotname = "diffsf_" + historthonames2d_sys;
  savehist2d(hist_diffsf,plotname,false);
  savehist2d(hist_diffsf_percent,"percent" + plotname,false);
  hist_diffsf->Write(plotname.c_str(),TObject::kOverwrite);
  hist_diffsf_percent->Write(("percent" + plotname).c_str(),TObject::kOverwrite);
  file->Close();

}

void calculatetriggsf(string datafilename, string mcfilename, string historthonames2d, string histtrignames2d,string syshist[] = {}, int nsyshist = 0)
{
    TFile *fdata = new TFile(datafilename.c_str(),"read");
    if(!fdata->IsOpen())
    {cout<<"Check: file "<<datafilename<<" not present"<<endl;exit(0);}

    TFile *fmc = new TFile(mcfilename.c_str(),"read");
    if(!fmc->IsOpen())
    {cout<<"Check: file "<<mcfilename<<" not present"<<endl;exit(0);}

     TH2D *hist2d_orthodata = (TH2D*)fdata->Get(historthonames2d.c_str());
    if(hist2d_orthodata == NULL)
    {cout<<"Check: file "<<datafilename<<" does not have TH2D hist "<<historthonames2d<<endl;exit(0);}

     TH2D *hist2d_trigdata = (TH2D*)fdata->Get(histtrignames2d.c_str());
    if(hist2d_trigdata == NULL)
    {cout<<"Check: file "<<datafilename<<" does not have TH2D hist "<<histtrignames2d<<endl;exit(0);}


    TH2D *hist2d_orthomc = (TH2D*)fmc->Get(historthonames2d.c_str());
    if(hist2d_orthomc == NULL)
    {cout<<"Check: file "<<mcfilename<<" does not have TH2D hist "<<historthonames2d<<endl;exit(0);}

    TH2D *hist2d_trigmc = (TH2D*)fmc->Get(histtrignames2d.c_str());
    if(hist2d_trigmc == NULL)
    {cout<<"Check: file "<<mcfilename<<" does not have TH2D hist "<<histtrignames2d<<endl;exit(0);}

    savehist2d(hist2d_orthomc, ("mc_" + historthonames2d).c_str());
    savehist2d(hist2d_trigmc, ("mc_" + histtrignames2d).c_str());
    savehist2d(hist2d_orthodata, ("data_" + historthonames2d).c_str());
    savehist2d(hist2d_trigdata, ("data_" + histtrignames2d).c_str());

    TEfficiency *Teff_mc = new TEfficiency(*hist2d_trigmc,*hist2d_orthomc);
    Teff_mc->SetStatisticOption(TEfficiency::kFFC);
    Teff_mc->SetConfidenceLevel(0.683);
    TH2D *hist_effmc = new TH2D();
    convertteff2hist2d(Teff_mc,hist_effmc);
    savehist2d(hist_effmc, ("ratio_mc_" + historthonames2d).c_str());

    TEfficiency *Teff_data = new TEfficiency(*hist2d_trigdata,*hist2d_orthodata);
    Teff_data->SetStatisticOption(TEfficiency::kFFC);
    Teff_data->SetConfidenceLevel(0.683);
    TH2D *hist_effdata = new TH2D();
    convertteff2hist2d(Teff_data,hist_effdata);
    savehist2d(hist_effdata, ("ratio_data_" + historthonames2d).c_str());

    TH2D*hist_sf = (TH2D*)Teff_mc->CreateHistogram();
    int nxbin = hist_sf->GetNbinsX();
    int nybin = hist_sf->GetNbinsY();

    for(int ij=1; ij<=nxbin; ij++)
        {
        for(int ik=1; ik<=nybin; ik++)
            {
                int globalbin = Teff_mc->GetGlobalBin(ij,ik);
                if(Teff_mc->GetEfficiency(globalbin)!=0)
                {
                    float sf = Teff_data->GetEfficiency(globalbin) / Teff_mc->GetEfficiency(globalbin);
                    hist_sf->SetBinContent(ij,ik,sf);
                    double error_mc,error_data;
                    error_mc = (Teff_mc->GetEfficiencyErrorUp(globalbin) + Teff_mc->GetEfficiencyErrorLow(globalbin) ) / 2;
                    error_data = (Teff_data->GetEfficiencyErrorUp(globalbin) + Teff_data->GetEfficiencyErrorLow(globalbin) ) / 2;
                    double error = sf * sqrt(error_mc*error_mc/(Teff_mc->GetEfficiency(globalbin)*Teff_mc->GetEfficiency(globalbin)) + error_data*error_data/(Teff_data->GetEfficiency(globalbin)*Teff_data->GetEfficiency(globalbin)) );
                    hist_sf->SetBinError(ij,ik,error);
                }
                else
                {
                    hist_sf->SetBinContent(ij,ik,0.0);
                    hist_sf->SetBinError(ij,ik,0.0);
                }
            }
    }
  TFile *file = new TFile("triggersf.root","update");
  file->cd();
  string plotname = "sf_" + historthonames2d;
  hist_sf->Write(plotname.c_str(),TObject::kOverwrite);
  savehist2d(hist_sf,plotname);
  file->Close();

  if(nsyshist > 0)
    {
        for(int isys = 0; isys< nsyshist; isys++)
        {
            string historthonames2d_sys = historthonames2d + "_" + syshist[isys];
            string histtrignames2d_sys = histtrignames2d + "_" + syshist[isys];
            calculatetriggsf(datafilename, mcfilename, historthonames2d_sys, histtrignames2d_sys);
            difftriggsf(datafilename, mcfilename, historthonames2d, histtrignames2d, historthonames2d_sys, histtrignames2d_sys);
        }
        //// overall sf calculation
        TFile *file = new TFile("triggersf.root","update");
        TH2D *hist_sf_nom;
        hist_sf_nom = (TH2D*)file->Get(("sf_" + historthonames2d).c_str());
        if(hist_sf_nom == NULL)
        {cout<<"Check: file triggersf.root does not have TH2D hist "<<("sf_" + historthonames2d).c_str()<<endl;exit(0);}

        TH2D *hist_sf_sys[nsyshist];
        for(int isys =0; isys < nsyshist; isys++)
        {
            hist_sf_sys[isys] = (TH2D*)file->Get(("sf_" + historthonames2d + "_" + syshist[isys]).c_str());
            if(hist_sf_sys[isys] == NULL)
            {cout<<"Check: file triggersf.root does not have TH2D hist "<<("sf_" + historthonames2d + "_" + syshist[isys]).c_str()<<endl;exit(0);}
        }

        TH2D *hist_sf_lumiweighted = (TH2D*)file->Get(("sf_" + historthonames2d + "_lumiweighted").c_str());

        TH2D *hist_sfoverall = (TH2D*)hist_sf_nom->Clone();
        TH2D *hist_sferrorpercent = (TH2D*)hist_sf_nom->Clone();

        int nxbin = hist_sfoverall->GetNbinsX();
        int nybin = hist_sfoverall->GetNbinsY();

        for(int ij=1; ij<=nxbin; ij++)
        {
            for(int ik=1; ik<=nybin; ik++)
            {
                double norm = hist_sf_nom->GetBinContent(ij,ik);
                double staterror = hist_sf_nom->GetBinError(ij,ik);
                double topsys_sferror = abs(norm - hist_sf_sys[0]->GetBinContent(ij,ik));
                double primivertex_sferror = max(abs(norm - hist_sf_sys[1]->GetBinContent(ij,ik)),abs(norm - hist_sf_sys[2]->GetBinContent(ij,ik)));
                double sfsys_lumiweighted = abs(norm - hist_sf_lumiweighted->GetBinContent(ij,ik));
                double error = sqrt(staterror*staterror + topsys_sferror*topsys_sferror + primivertex_sferror*primivertex_sferror + sfsys_lumiweighted*sfsys_lumiweighted);
                hist_sfoverall->SetBinContent(ij,ik,norm);
                hist_sfoverall->SetBinError(ij,ik,error);
                if(norm > 0) hist_sferrorpercent->SetBinContent(ij,ik,error/norm*100);
            }
        }

        string plotname = "overallsf_" + historthonames2d;
        savehist2d(hist_sfoverall,plotname);

        hist_sfoverall->Write(plotname.c_str(),TObject::kOverwrite);

        string plotname2 = "allerrorpercent_sf_" + historthonames2d;
        savehist2d(hist_sferrorpercent,plotname2,false);

        hist_sferrorpercent->Write(plotname2.c_str(),TObject::kOverwrite);

        file->Close();
    }
}

void makelumiweightedsf(fileinfo lumiweightedjethtdata, string mcfilename, string historthonames2d, string histtrignames2d)
{
    const int ndatafiles = (int)lumiweightedjethtdata.in_names.size();
    for(int ip=0;ip<ndatafiles;ip++)
    calculatetriggsf(lumiweightedjethtdata.in_names[ip], mcfilename, historthonames2d + "_lumiweighted" + to_string(ip), histtrignames2d + "_lumiweighted" + to_string(ip));

    TFile *file = new TFile("triggersf.root","update");

    TH2D *hist_sf_dataperiods[ndatafiles];
    float totlumi = 0;
    for(int ip=0;ip<ndatafiles;ip++)
    {
        hist_sf_dataperiods[ip] = (TH2D*)file->Get(("sf_" + historthonames2d + "_lumiweighted" + to_string(ip)).c_str());
        totlumi += lumiweightedjethtdata.weights[ip];
        if(hist_sf_dataperiods[ip] == NULL)
        {cout<<"CHeck file triggersf.root does not contain hist"<<("sf_" + historthonames2d + "_lumiweighted" + to_string(ip))<<"  Probably some naming issue."<<endl;exit(0);}
    }
    TH2D *hist_sf_lumiweighted = (TH2D*)hist_sf_dataperiods[0]->Clone();

    int nxbin = hist_sf_lumiweighted->GetNbinsX();
    int nybin = hist_sf_lumiweighted->GetNbinsY();

        for(int ij=1; ij<=nxbin; ij++)
        {
            for(int ik=1; ik<=nybin; ik++)
            {
                float lumiweightedsf = 0;
                for(int ip=0;ip<ndatafiles;ip++)
                    lumiweightedsf += hist_sf_dataperiods[ip]->GetBinContent(ij,ik)*lumiweightedjethtdata.weights[ip];
                lumiweightedsf = lumiweightedsf/totlumi;
                hist_sf_lumiweighted->SetBinContent(ij,ik,lumiweightedsf);
            }
        }
    file->cd();
    string plotname = "sf_" + historthonames2d + "_lumiweighted";
    hist_sf_lumiweighted->Write(plotname.c_str(),TObject::kOverwrite);
    savehist2d(hist_sf_lumiweighted,plotname);
    file->Close();

}

void plotrootfile(string histname, vector<fileinfo> allfiles)
{

    for(int ifile=0;ifile<allfiles.size();ifile++)
    {
        if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}

        TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");

        if(!fileout->IsOpen())
        {cout<<"Check: file "<<allfiles[ifile].out_name<<" not present"<<endl;exit(0);}

        TH2D* hist2d = (TH2D*)fileout->Get(histname.c_str());

        if(hist2d == NULL)
        {cout<<"Check: file "<<allfiles[ifile].out_name<<" does not have 2D Hist "<<histname<<endl;exit(0);}

        string plotname = "fordatadivenqcd_" +allfiles[ifile].out_name.substr(30,9)+ "_" + histname;
        cout<<"Correlation for hist"<<plotname<<" = "<<hist2d->GetCorrelationFactor()<<endl;
        savehist2d(hist2d,plotname,false);
    }
}



void plotfromhiggscombine(string channel, string region, int nprocs, string processes[nprocs], string channelname[nprocs], string ptbins = "")
{
    using namespace RooFit;

    ///  Open higgs combine output file to get fitted plots
    TFile *fitfile;
    string filename = "fitDiagnostics_" + channel + ptbins + ".root";
    fitfile = new TFile(filename.c_str(),"read");
    if(!fitfile->IsOpen())
        {
            cout<<"Check: Higgs combine outout file "<<filename<<" not present. Or check its name should be in format fitDiagnostics_*channel*.root"<<endl;exit(0);
        }

   /*
    Get roocurve for pdfs

    string fitname = region + channel + "_CMS_th1x_fit_s";
    RooPlot *fit = (RooPlot*)fitfile->Get(fitname.c_str());
    if(fit == NULL)
    {       cout<<"Check: roofit plot with name "<<fitname<<" in file "<<filename<<" is not present"<<endl;exit(0);}

    string dataname = "h_" + region + channel;
    string resultfitname = "pdf_bin" + region + channel + "_Norm[CMS_th1x]";
    string sigfitname = "pdf_bin" + region + channel + "_Norm[CMS_th1x]_Comp[shapeSig*]";
    string bkgfitname = "pdf_bin" + region + channel + "_Norm[CMS_th1x]_Comp[shapeBkg*]";

    RooHist * data = (RooHist*)fit->findObject(dataname.c_str());
    RooCurve* result_fit_pdf= (RooCurve*)fit->findObject(resultfitname.c_str());
    RooCurve* sig_fit_pdf= (RooCurve*)fit->findObject(sigfitname.c_str());
    RooCurve* bkg_fit_pdf= (RooCurve*)fit->findObject(bkgfitname.c_str());

    if(data == NULL)
    {       cout<<"Check: data roohist "<<dataname<<"in rooplot ("<<fitname<<") in file "<<filename<<" is not present. Make sure that --plots option was used in higgs combine to get this."<<endl;exit(0);}

    if(result_fit_pdf == NULL)
    {       cout<<"Check: fit pdf roocurve "<<resultfitname<<"in rooplot ("<<fitname<<") in file "<<filename<<" is not present."<<endl;exit(0);}

    if(sig_fit_pdf == NULL)
    {       cout<<"Check: signal pdf roocurve "<<sigfitname<<"in rooplot ("<<fitname<<") in file "<<filename<<" is not present."<<endl;exit(0);}

    if(bkg_fit_pdf == NULL)
    {       cout<<"Check: background pdf roocurve "<<bkgfitname<<"in rooplot ("<<fitname<<") in file "<<filename<<" is not present."<<endl;exit(0);}

    string histfitname = "shapes_fit_s/"+region+channel+"/total";
    TH1D *hist_fit_messedupbins = (TH1D*)fitfile->Get(histfitname.c_str());
    if(hist_fit_messedupbins == NULL)
    {       cout<<"Check: fitted hist "<<histfitname<<"in rooplot ("<<fitname<<") in file "<<filename<<" is not present. Make sure that ----saveShapes option was used in higgs combine to get this."<<endl;exit(0);}

       int nbins_data = hist_data->GetNbinsX();
    cout<<nbins_data<<endl;
    if(nbins_data != hist_result_fit->GetNbinsX()) {cout<<"Check input file names! The no. of bins in prefit data is not equal to no. of points in roohist data of higgscombine output files"<<endl;exit(0);}
    for(int ib = 0; ib < nbins_data; ib++)
    {
        cout<<hist_data->GetBinWidth(ib+1)<<endl;
        data->SetPointY(ib,data->GetPointY(ib)/hist_data->GetBinWidth(ib+1));
        for(int ip = int(ib*result_fit->GetN()/nbins_data); ip <int((ib+1)*result_fit->GetN()/nbins_data); ip++)
            result_fit->SetPointY(ip,result_fit->GetPointY(ip)/hist_data->GetBinWidth(ib+1));
        for(int ip = int(ib*sig_fit->GetN()/nbins_data); ip <int((ib+1)*sig_fit->GetN()/nbins_data); ip++)
            sig_fit->SetPointY(ip,sig_fit->GetPointY(ip)/hist_data->GetBinWidth(ib+1));
        for(int ip = int(ib*bkg_fit->GetN()/nbins_data); ip <int((ib+1)*bkg_fit->GetN()/nbins_data); ip++)
            bkg_fit->SetPointY(ip,bkg_fit->GetPointY(ip)/hist_data->GetBinWidth(ib+1));
    }

      TGraphAsymmErrors *h_ratio = new TGraphAsymmErrors(hist_data,hist_fit,"pois");

    TLegend *legend = tdrLeg(0.72,0.65,0.9,0.915);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    legend->AddEntry(data, "Data", "ple1");
    legend->AddEntry(result_fit, "Fit", "l");
    legend->AddEntry(sig_fit, "sig Pdf", "l");
    legend->AddEntry(bkg_fit, "Bkg Pdf", "l");

    TCanvas * can = tdrDiCanvas("can",(TH1D*)data->GetHistogram(),hist_data,8,0);
    can->cd(1);
    //gPad->SetLogy(1);
    data->GetYaxis()->SetTitle("Events / GeV");
    data->GetXaxis()->SetTitle(" ");
    data->setYAxisLimits(data->getYAxisMin(),1.2*data->getYAxisMax());
    data->GetYaxis()->CenterTitle(1);
    data->GetXaxis()->CenterTitle(1);

    data->Draw("P");
    result_fit->Draw("same");
    sig_fit->Draw("same");
    bkg_fit->Draw("same");
    legend->Draw("same");

    */


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

    TH1D *hist_procs_prefit[nprocs];
    for(int ip = 0; ip < nprocs; ip++)
    {
        string prochistname = "sdmass_probejet_" + region + "_" + processes[ip] + "_" + channel;
        hist_procs_prefit[ip] = (TH1D*)filein->Get(prochistname.c_str());
        cout<<" check mc name = "<<hist_procs_prefit[ip]->GetName()<<" "<<hist_procs_prefit[ip]->GetTitle()<<endl;
        if(hist_procs_prefit[ip] == NULL)
        {cout<<"Check: file "<<datafilename<<" does not have histogram for process  "<<processes[ip]<<" with hist name:"<<prochistname<<endl;exit(0);}
    }

    /// correct the visualizing of the histogram (useful only for variable bin width histogram).
    int nbins_data = hist_data->GetNbinsX();
    cout<<nbins_data<<endl;
    if(nbins_data != hist_result_fit_improperbin->GetNbinsX()) {cout<<"Check input file names! The no. of bins in prefit data is not equal to no. of points in roohist data of higgscombine output files"<<endl;exit(0);}


    cout<<"\n"<<region<<" "<<channel<<" \nData prefit = "<<hist_data->Integral()<<" Prefit signal yield = "<<hist_procs_prefit[0]->Integral()<<" \nPostfit yield (sig+bkg) = "<<hist_result_fit_improperbin->Integral()<<" signal yield = "<<hist_sig_fit_improperbin->Integral()<<endl<<endl;

    /// Get proper bining from orginial histograms
    for(int ip = 0; ip < nprocs; ip++)
    {
        hist_procs_postfit[ip] = (TH1D*)hist_data->Clone();
        hist_procs_prefit[ip] = (TH1D*)hist_data->Clone();
    }
    TH1D* hist_result_fit = (TH1D*)hist_data->Clone();
    TH1D* hist_sig_fit = (TH1D*)hist_data->Clone();
    TH1D* hist_bkg_fit = (TH1D*)hist_data->Clone();

    for(int ib = 0; ib < nbins_data; ib++)
    {
        //float ith_databin_width = hist_data->GetBinWidth(ib+1);
        float ith_databin_width = 100;
        cout<<hist_data->GetBinContent(ib+1)<<"*100 / "<<ith_databin_width;

        hist_data->SetBinContent(ib+1,hist_data->GetBinContent(ib+1)*100/ith_databin_width);
        hist_data->SetBinError(ib+1,hist_data->GetBinError(ib+1)*100/ith_databin_width);
        cout<<" = "<<hist_data->GetBinContent(ib+1)<<endl;

        cout<<hist_result_fit_improperbin->GetBinContent(ib+1)<<"*100 / "<<ith_databin_width;

        hist_result_fit->SetBinContent(ib+1,hist_result_fit_improperbin->GetBinContent(ib+1)*100/ith_databin_width);
        hist_result_fit->SetBinError(ib+1,hist_result_fit_improperbin->GetBinError(ib+1)*100/ith_databin_width);
        hist_result_fit->SetLineWidth(3);

        cout<<" = "<<hist_result_fit->GetBinContent(ib+1)<<endl;

        hist_sig_fit->SetBinContent(ib+1,hist_sig_fit_improperbin->GetBinContent(ib+1)*100/ith_databin_width);
        hist_sig_fit->SetBinError(ib+1,hist_sig_fit_improperbin->GetBinError(ib+1)*100/ith_databin_width);

        hist_bkg_fit->SetBinContent(ib+1,hist_bkg_fit_improperbin->GetBinContent(ib+1)*100/ith_databin_width);
        hist_bkg_fit->SetBinError(ib+1,hist_bkg_fit_improperbin->GetBinError(ib+1)*100/ith_databin_width);

    for(int ip = 0; ip < nprocs; ip++)
    {
        hist_procs_postfit[ip]->SetBinContent(ib+1,hist_procs_postfit_improperbin[ip]->GetBinContent(ib+1)*100/ith_databin_width);
        hist_procs_postfit[ip]->SetBinError(ib+1,hist_procs_postfit_improperbin[ip]->GetBinError(ib+1)*100/ith_databin_width);
        hist_procs_postfit[ip]->SetLineStyle(1);
        hist_procs_postfit[ip]->SetLineWidth(3);
        hist_procs_postfit[ip]->SetLineColor(ip+6);

        if(ip == 1)
                    cout<<hist_data->GetBinContent(ib+1)<<" = "<<hist_procs_prefit[ip]->GetBinContent(ib+1)<<endl;

        hist_procs_prefit[ip]->SetBinContent(ib+1,hist_procs_prefit[ip]->GetBinContent(ib+1)*100/ith_databin_width);
        hist_procs_prefit[ip]->SetBinError(ib+1,hist_procs_prefit[ip]->GetBinError(ib+1)*100/ith_databin_width);
        hist_procs_prefit[ip]->SetLineStyle(2);
        hist_procs_prefit[ip]->SetLineWidth(3);
        hist_procs_prefit[ip]->SetLineColor(ip+6);
    }
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

/*    RooArgSet *norms = (RooArgSet*)fitfile->Get("norm_fit_s");
    if(norms == NULL)
    {       cout<<"Check: normalizations RooArgSet norm_fit_s in file "<<filename<<" is not present. Make sure that --saveNormalizations option was used in higgs combine to get this."<<endl;exit(0);}

    //cout<<" norm = "<<norms->Print();*/

    //TH1D *hist_fit = (TH1D*)hist_data->Clone();

/*    for(int i =1; i<hist_data->GetNbinsX()+1; i++)
    {
        hist_fit->SetBinContent(i,hist_fit_messedupbins->GetBinContent(i));
        hist_fit->SetBinError(i,hist_fit_messedupbins->GetBinError(i));
        cout<<"data bin "<<i<<" "<<hist_data->GetBinLowEdge(i)<<" to "<<hist_data->GetBinLowEdge(i)+hist_data->GetBinWidth(i)<<" "<<hist_data->GetBinContent(i)<<" fit bin "<<i<<" "<<hist_fit->GetBinLowEdge(i)<<" to "<<hist_fit->GetBinLowEdge(i)+hist_fit->GetBinWidth(i)<<" "<<hist_fit_messedupbins->GetBinContent(i)<<" "<<hist_fit->GetBinContent(i)<<endl;
    }*/
    TGraphAsymmErrors *h_ratio = new TGraphAsymmErrors(hist_data,hist_result_fit,"pois");

//    string channelname[3] = {"top matched","leptonic or b matched","no matched"};

    TLegend *legend = tdrLeg(0.54,0.6,0.84,0.85);
    legend->SetBorderSize(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.035);
    legend->AddEntry(hist_data, "Data", "6ple1");
    legend->AddEntry(hist_result_fit, "Total Fit", "l");
//    legend->AddEntry(hist_sig_fit, "Signal", "l");
//    legend->AddEntry(hist_bkg_fit, "Background", "l");
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
    hist_data->GetYaxis()->SetTitle("Events / 100 GeV");
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
//    hist_sig_fit->Draw("histsame");
//    hist_bkg_fit->Draw("histsame");
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
       // cout<<"band bin "<<i+1<<" "<<band->GetBinLowEdge(i+1)<<" to "<<band->GetBinLowEdge(i+1)+hist_data->GetBinWidth(i+1)<<" "<<band->GetBinError(i+1)<<endl;

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

    //hist_data->Draw("P");
    hist_prefit_stack->Draw("hist");
//    hist_sig_fit->Draw("histsame");
//    hist_bkg_fit->Draw("histsame");
    legend_prefit->Draw("same");
    hist_data->Draw("PE1SAME");

    can_prefit->cd(2);
    h_ratio_prefit->SetMarkerColor(1);
    h_ratio_prefit->SetMarkerStyle(20);
    h_ratio_prefit->SetMarkerSize(0.89);
    h_ratio_prefit->SetLineColor(kBlack);
    TH1D *band_prefit = new TH1D("Band_prefit", "", nbin, bins);

    for(int i=0;i<nbin; i++)
    {
        band_prefit->SetBinContent(i+1, 1.0);
        band_prefit->SetBinError(i+1, 0);
    }

    gPad->SetGridx();
    gPad->SetGridy();
//    band_prefit->SetFillColor(kGray+3);
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
    pt->Draw("same");
    line->Draw("same");

    can_prefit->Update();
    can_prefit->SaveAs(("plots/Final_combine_fit_"+region+ptbins+channel+"_prefitonly.png").c_str());






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

    //hist_data->Draw("P");
    hist_postfit_stack->Draw("hist");
    hist_result_fit->Draw("histsame");
//    hist_sig_fit->Draw("histsame");
//    hist_bkg_fit->Draw("histsame");
    legend_postfit->Draw("same");
    hist_data->Draw("PE1SAME");

    can_postfit->cd(2);

    gPad->SetGridx();
    gPad->SetGridy();

    band->Draw("E2C");
    h_ratio->Draw("PE1SAME");
    pt->Draw("same");

    can_postfit->Update();
    can_postfit->SaveAs(("plots/Final_combine_fit_"+region+ptbins+channel+"_postfitonly.png").c_str());


    filein->Close();
    fitfile->Close();

}

void plot_datamcratio( const int nfiles, string dataname[nfiles],  int col[nfiles-1], string obsnames, string vartitle,  vector<fileinfo> allfiles)
{
    const int nfiles_input = (int)allfiles.size();
    if(nfiles != nfiles_input)
    {cout<<"Check file info: No. of datanames given is not equal to no. of input files ie "<<nfiles<<" != "<<nfiles_input<<endl;exit(0);}

  TH1D *hist_obs[nfiles];

  for(int ifile=0;ifile<nfiles;ifile++)
    {
        TFile *filein = new TFile(allfiles[ifile].out_name.c_str(),"read");

        if(!filein->IsOpen())
            {cout<<"Check: file "<<allfiles[ifile].out_name<<" not present"<<endl;exit(0);}

        filein->cd();
        hist_obs[ifile] = (TH1D*)filein->Get(obsnames.c_str());

        if(hist_obs[ifile] == NULL)
        {
            cout<<allfiles[ifile].out_name<<" file does not have or have empty histogram  "<<obsnames<<" hist";
            continue;
        }
        hist_obs[ifile]->SetDirectory(0);
        filein->Close();
    }

    THStack *hist_obs_stack;

  TH1D *hist_obs_bkg;

  TH1D *rat_obs;
  TGraphAsymmErrors *ratioGraph;
  TCanvas *cv;
  TLegend *legv;

  hist_obs[0]->GetXaxis()->SetTitle(" ");
  char name[200];
  sprintf(name,"Stack_%s",obsnames.c_str());
  hist_obs_stack = new THStack(name,"");
  for(int bkg=1; bkg<(nfiles); bkg++){

    hist_obs[nfiles-bkg]->SetFillColor(col[nfiles-bkg-1]);
    hist_obs[nfiles-bkg]->SetLineColor(col[nfiles-bkg-1]);
    hist_obs[nfiles-bkg]->SetLineWidth(3);
    hist_obs[nfiles-bkg]->GetYaxis()->SetTitleOffset(1.2);

    hist_obs_stack->Add(hist_obs[nfiles-bkg]);
    }

    for(int bkg=0; bkg<(nfiles-1); bkg++){
      if(bkg==0) { hist_obs_bkg = (TH1D*)hist_obs[bkg+1]->Clone(); }
      else { hist_obs_bkg->Add(hist_obs[bkg+1]); }
    }
    // canv obs //
    sprintf(name,"Canv_%s",obsnames.c_str());

    cv = tdrDiCanvas(name,hist_obs[0],hist_obs_bkg,8,0);

    legv = tdrLeg(0.65,0.65,0.875,0.915);

    legv->SetTextFont(42);
    legv->SetTextSize(0.045);
    legv->SetBorderSize(0);
    for(int bkg=0; bkg<(nfiles); bkg++){
      if(bkg==0){ legv->AddEntry(hist_obs[bkg],dataname[bkg].c_str(),"lep"); }
      else{
          legv->AddEntry(hist_obs[bkg],dataname[bkg].c_str(),"f");
      }
    }
    cv->cd(1);

    hist_obs[0]->SetFillStyle(0);
    hist_obs[0]->SetFillColor(0);
    hist_obs[0]->SetMarkerStyle(kFullCircle);
    hist_obs[0]->SetMarkerColor(kBlack);
    hist_obs[0]->SetMarkerSize(0.7);
    hist_obs[0]->SetLineColor(kBlack);
    hist_obs[0]->GetYaxis()->SetTitleOffset(1.5);

    //gPad->SetLogy(1);

    hist_obs[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)* max(hist_obs[0]->GetMaximum(),hist_obs_bkg->GetMaximum()));
    hist_obs_bkg->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)*max(hist_obs[0]->GetMaximum(),hist_obs_bkg->GetMaximum()));
    hist_obs[0]->SetMinimum(0.5);
    hist_obs_bkg->SetMinimum(0.5);

    hist_obs_stack->Draw("histSAME");
    hist_obs[0]->Draw("PSAME");

    gPad->RedrawAxis();

    cv->cd(2);
    rat_obs = (TH1D*)hist_obs[0]->Clone();
    rat_obs->Divide(hist_obs_bkg);

/*    if(iv==1 || iv==2){
      gPad->SetLogx(0);
      rat_obs->GetXaxis()->SetMoreLogLabels(kTRUE);
      rat_obs->GetXaxis()->SetNoExponent(kTRUE);
    }
*/
    sprintf(name,"%s",vartitle.c_str());
    rat_obs->GetXaxis()->SetTitle(name);
    rat_obs->GetXaxis()->SetTitleSize(0.1);
    rat_obs->GetXaxis()->SetLabelSize(0.1);
    //    rat_obs->GetXaxis()->CenterTitle();
    rat_obs->GetXaxis()->SetNdivisions(406);

    rat_obs->GetYaxis()->SetTitle("Data / MC");
    rat_obs->GetYaxis()->SetTitleSize(0.1);
    rat_obs->GetYaxis()->SetLabelSize(0.1);
    rat_obs->GetYaxis()->SetTitleOffset(0.75);
    rat_obs->GetYaxis()->SetNdivisions(406);
    rat_obs->GetYaxis()->CenterTitle();

    rat_obs->SetMinimum(0);
    rat_obs->SetMaximum(2.);

    rat_obs->SetFillStyle(0);
    rat_obs->SetFillColor(0);
    rat_obs->SetMarkerStyle(kFullCircle);
    rat_obs->SetMarkerColor(kBlack);
    rat_obs->SetMarkerSize(0.7);
    rat_obs->SetLineColor(kBlack);
    rat_obs->Draw("e,p");

    TLine *line = new TLine(rat_obs->GetXaxis()->GetXmin(),1,rat_obs->GetXaxis()->GetXmax(),1);
    line->SetLineColor(1);
    line->Draw("sames");
    sprintf(name,"plots/%s.png",cv->GetName());
    cv->SaveAs(name);

}

void plot_datamcratio( const int nfiles, TH1D *hist_obs[nfiles], string dataname[nfiles],  int col[nfiles-1], string obsnames, string vartitle)
{
  THStack *hist_obs_stack;

  TH1D *hist_obs_bkg;

  TH1D *rat_obs;
  TGraphAsymmErrors *ratioGraph;
  TCanvas *cv;
  TLegend *legv;

  hist_obs[0]->GetXaxis()->SetTitle(" ");
  char name[200];
  sprintf(name,"Stack_%s",obsnames.c_str());
  hist_obs_stack = new THStack(name,"");
  for(int bkg=1; bkg<(nfiles); bkg++){

    hist_obs[nfiles-bkg]->SetFillColor(col[nfiles-bkg-1]);
    hist_obs[nfiles-bkg]->SetLineColor(col[nfiles-bkg-1]);
    hist_obs[nfiles-bkg]->SetLineWidth(3);
    hist_obs[nfiles-bkg]->GetYaxis()->SetTitleOffset(1.2);

    hist_obs_stack->Add(hist_obs[nfiles-bkg]);
    }

    for(int bkg=0; bkg<(nfiles-1); bkg++){
      if(bkg==0) { hist_obs_bkg = (TH1D*)hist_obs[bkg+1]->Clone(); }
      else { hist_obs_bkg->Add(hist_obs[bkg+1]); }
    }
    // canv obs //
    sprintf(name,"Canv_%s",obsnames.c_str());

    cv = tdrDiCanvas(name,hist_obs[0],hist_obs_bkg,8,0);

    legv = tdrLeg(0.65,0.65,0.875,0.915);

    legv->SetTextFont(42);
    legv->SetTextSize(0.045);
    legv->SetBorderSize(0);
    for(int bkg=0; bkg<(nfiles); bkg++){
      if(bkg==0){ legv->AddEntry(hist_obs[bkg],dataname[bkg].c_str(),"lep"); }
      else{
          legv->AddEntry(hist_obs[bkg],dataname[bkg].c_str(),"f");
      }
    }
    cv->cd(1);

    hist_obs[0]->SetFillStyle(0);
    hist_obs[0]->SetFillColor(0);
    hist_obs[0]->SetMarkerStyle(kFullCircle);
    hist_obs[0]->SetMarkerColor(kBlack);
    hist_obs[0]->SetMarkerSize(0.7);
    hist_obs[0]->SetLineColor(kBlack);
    hist_obs[0]->GetYaxis()->SetTitleOffset(1.5);

    gPad->SetLogy(1);

    hist_obs[0]->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)*max(hist_obs[0]->GetMaximum(),hist_obs_bkg->GetMaximum()));
    hist_obs_bkg->SetMaximum((gPad->GetLogy() == 0 ? 1.4 : 30)*max(hist_obs[0]->GetMaximum(),hist_obs_bkg->GetMaximum()));
    hist_obs[0]->SetMinimum(0.5);
    hist_obs_bkg->SetMinimum(0.5);

    hist_obs_stack->Draw("histSAME");
    hist_obs[0]->Draw("PSAME");

    gPad->RedrawAxis();

    cv->cd(2);
    rat_obs = (TH1D*)hist_obs[0]->Clone();
    rat_obs->Divide(hist_obs_bkg);

/*    if(iv==1 || iv==2){
      gPad->SetLogx(0);
      rat_obs->GetXaxis()->SetMoreLogLabels(kTRUE);
      rat_obs->GetXaxis()->SetNoExponent(kTRUE);
    }
*/
    sprintf(name,"%s",vartitle.c_str());
    rat_obs->GetXaxis()->SetTitle(name);
    rat_obs->GetXaxis()->SetTitleSize(0.1);
    rat_obs->GetXaxis()->SetLabelSize(0.1);
    //    rat_obs->GetXaxis()->CenterTitle();
    rat_obs->GetXaxis()->SetNdivisions(406);

    rat_obs->GetYaxis()->SetTitle("Data / MC");
    rat_obs->GetYaxis()->SetTitleSize(0.1);
    rat_obs->GetYaxis()->SetLabelSize(0.1);
    rat_obs->GetYaxis()->SetTitleOffset(0.75);
    rat_obs->GetYaxis()->SetNdivisions(406);
    rat_obs->GetYaxis()->CenterTitle();

    rat_obs->SetMinimum(0);
    rat_obs->SetMaximum(2.);

    rat_obs->SetFillStyle(0);
    rat_obs->SetFillColor(0);
    rat_obs->SetMarkerStyle(kFullCircle);
    rat_obs->SetMarkerColor(kBlack);
    rat_obs->SetMarkerSize(0.7);
    rat_obs->SetLineColor(kBlack);
    rat_obs->Draw("e,p");

    TLine *line = new TLine(rat_obs->GetXaxis()->GetXmin(),1,rat_obs->GetXaxis()->GetXmax(),1);
    line->SetLineColor(1);
    line->Draw("sames");
    sprintf(name,"plots/%s.png",cv->GetName());
    cv->SaveAs(name);

}

void compare_datamc_fromallfiles(const int nfiles, string dataname[nfiles],  int col[nfiles-1], string obsnames, string vartitle, vector<fileinfo> allfiles)
{
    const int nfiles_input = (int)allfiles.size();
    if(nfiles != nfiles_input)
    {cout<<"Check file info: No. of datanames given is not equal to no. of input files ie "<<nfiles<<" != "<<nfiles_input<<endl;exit(0);}

    TH1D* hist_datamc[nfiles];

    for(int ifile=0;ifile<nfiles;ifile++)
    {
        if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}

        cout<<allfiles[ifile].out_name<<endl;
        TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");
            cout<<"Check 1"<<endl;
            TH1D *hist_out;
        //hist_out->Sumw2();
        for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {   cout<<"Check 2"<<endl;

            TFile *filein = new TFile(allfiles[ifile].in_names[ig].c_str(),"read");
  cout<<"Check 3"<<endl;

            if(!filein->IsOpen())
            {cout<<"Check: file "<<allfiles[ifile].in_names[ig]<<" not present"<<endl;exit(0);}
  cout<<"Check 4"<<endl;

            filein->cd();
            TH1D *hist_in = (TH1D*)filein->Get(obsnames.c_str());
  cout<<"Check 5"<<endl;

            //hist_in->SetDirectory(0);
            if(hist_in == NULL)
            {
                cout<<allfiles[ifile].in_names[ig]<<" file does not have or have empty histogram  "<<obsnames<<" hist";
                continue;
            }
           // hist_in->Scale(allfiles[ifile].weights[ig]);

            if(ig == 0)
            {
                hist_out = (TH1D*)hist_in->Clone();
                hist_out->SetDirectory(0);
            }
        else
            hist_out->Add(hist_in);
  cout<<"Check 6"<<endl;

            filein->Close();   cout<<"Check 7"<<endl;

        }
       /* if(ifile != 0) {
            hist_out->Scale(59.83);
        }*/
        fileout->cd();
        hist_out->Write(obsnames.c_str(),TObject::kOverwrite);

        if(hist_out == NULL)
            {cout<<"Histogram is empty cant clone it. Please check"<<endl; exit(0);}
        else
            hist_datamc[ifile] = (TH1D*)hist_out->Clone();

        hist_datamc[ifile]->SetDirectory(0);
        fileout->Close();
    }
    plot_datamcratio(nfiles, hist_datamc, dataname, col, obsnames, vartitle);
}

void compare_datamc(const int nfiles, string dataname[nfiles],  int col[nfiles-1], string obsnames, string vartitle, vector<fileinfo> allfiles)
{
    const int nfiles_input = (int)allfiles.size();
    if(nfiles != nfiles_input)
    {cout<<"Check file info: No. of datanames given is not equal to no. of input files ie "<<nfiles<<" != "<<nfiles_input<<endl;exit(0);}

    TH1D* hist_datamc[nfiles];

    for(int ifile=0;ifile<nfiles;ifile++)
    {
        if(allfiles[ifile].in_names.size() != allfiles[ifile].weights.size())
        {cout<<"Check file info:";cout<<" "<<allfiles[ifile].in_names[0]<<" no. of files is not equal to no. of given weights"<<endl;exit(0);}

        cout<<allfiles[ifile].out_name<<endl;
        TFile *fileout = new TFile(allfiles[ifile].out_name.c_str(),"update");
        TH1D *hist_out;
        //hist_out->Sumw2();
        for(int ig=0;ig<allfiles[ifile].in_names.size();ig++)
        {
            TFile *filein = new TFile(allfiles[ifile].in_names[ig].c_str(),"read");

            if(!filein->IsOpen())
            {cout<<"Check: file "<<allfiles[ifile].in_names[ig]<<" not present"<<endl;exit(0);}

            filein->cd();
            TH1D *hist_in = (TH1D*)filein->Get(obsnames.c_str());

            //hist_in->SetDirectory(0);

            hist_in->Scale(allfiles[ifile].weights[ig]);

            if(ig == 0)
        {
            hist_out = (TH1D*)hist_in->Clone();
            hist_out->SetDirectory(0);
        }
        else
            hist_out->Add(hist_in);

            filein->Close();
        }
        if(ifile != 0) {
            hist_out->Scale(59.83);
        }
        fileout->cd();
        hist_out->Write(obsnames.c_str(),TObject::kOverwrite);

        if(hist_out == NULL)
            {cout<<"Histogram is empty cant clone it. Please check"<<endl; exit(0);}
        else
            hist_datamc[ifile] = (TH1D*)hist_out->Clone();

        hist_datamc[ifile]->SetDirectory(0);
        fileout->Close();
    }
    plot_datamcratio(nfiles, hist_datamc, dataname, col, obsnames, vartitle);
}

#endif
