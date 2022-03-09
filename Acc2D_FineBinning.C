///////////////////////////////////////////////////////////////////////////
//The purpose of this program is to calculate the 2-Dimensional
//Upsilon y vs Upsilon pT Acceptance spectrum for uneven events and use
//the even events spectrum calculated at Acc2D_even to do the fine binning
//process.
//
//To run, type on console:
//.L Acc2D_FineBinning.C
//Acceptance t
//t.Loop()
///////////////////////////////////////////////////////////////////////////

//defining headers
#define Acceptance_cxx
#include "Acceptance.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//Defining program
void Acceptance::Loop()
{
  //reading inputfile
  TFile* InputFile = new TFile("OutAcc2D_even.root","READ");

  //Loading histograms from TFiles objects
  TH2F* hGenUpsilonAccEven = (TH2F*)InputFile->Get("hGenUpsilonAccEven");

  if (fChain == 0) return;

  //Defining pt binning
  Float_t PtBins[29] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.75,7.5,8.25,9.0,9.75,10.5,11.25,12.0,16.5,21.0,25.5,30.0,35.0,40.0,45.0,50.0};
  Int_t PtnBins= 28;

  //Defining Acceptance Denominator Histogram
  TH1F* hGenUpsilonAccDen = new TH1F ("hGenUpsilonAccDen", "1D Acceptance denominator for uneven events",PtnBins, PtBins); hGenUpsilonAccDen->Sumw2();
  //Defining Acceptance numerator histogram
  TH1F* hGenUpsilonAccNum = new TH1F ("hGenUpsilonAccNum", "1D Acceptance numerator for uneven events",PtnBins, PtBins); hGenUpsilonAccNum->Sumw2();
  //Defining Acceptance corrected numerator histogram
  TH1F* hGenUpsilonAccNumCorr = new TH1F ("hGenUpsilonAccNumCorr", "1D Acceptance Corrected numerator for uneven events",PtnBins, PtBins); hGenUpsilonAccNumCorr->Sumw2();

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)//Loop on ALL events
  //for (Long64_t jentry=0; jentry<100000;jentry++)//Limited to read only 10000 events
    {
      //check if entry is uneven; if not, skip event
      if( jentry % 2 == 0) continue;
      
      cout<<"Reading event "<<jentry<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    
      for(int iQQ=0; iQQ < Gen_QQ_size; iQQ++)//Loop on generated upsilon
	{
	  TLorentzVector *GenQQ = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	  float UpsilonPt = GenQQ->Pt();
	  float UpsilonEta = GenQQ->Eta();
	  float UpsilonRapidity = GenQQ->Rapidity();

	  TLorentzVector *GenMuPl = (TLorentzVector*) Gen_QQ_mupl_4mom->At(iQQ);
	  float MuPlPt = GenMuPl->Pt();
	  float MuPlEta = GenMuPl->Eta();

	  TLorentzVector *GenMuMi = (TLorentzVector*) Gen_QQ_mumi_4mom->At(iQQ);
	  float MuMiPt = GenMuMi->Pt();
	  float MuMiEta = GenMuMi->Eta();

	  //Defining Acceptance conditions as booleans
	  bool UpsilonRapMax = abs(UpsilonRapidity) < 2.4;
       
          bool MuPlPtMin = MuPlPt > 3.5;
	  bool MuPlEtaMax = abs(MuPlEta) < 2.4;
	   
          bool MuMiPtMin = MuMiPt > 3.5;
	  bool MuMiEtaMax = abs(MuMiEta) < 2.4;

	  bool MuonPtMin = (MuPlPtMin && MuMiPtMin);
	  bool MuonEtaMax = (MuPlEtaMax && MuMiEtaMax);

	  if( UpsilonRapMax == false ) continue;
	  //Check for conditions and fill denominator histograms
	  hGenUpsilonAccDen->Fill(UpsilonPt);
	  
	  if( ((MuPlPtMin==false) || (MuMiPtMin==false)) || ((MuPlEtaMax==false) || (MuMiEtaMax==false)) ) continue;
	  //Check for conditions and fill numerator histograms

	  //Get weight value using the corresponding bin value from the
	  //Acceptance even events spectrum
	  Int_t BinNum = hGenUpsilonAccEven->FindBin(UpsilonPt,UpsilonEta);
	  Double_t WeightVal = hGenUpsilonAccEven->GetBinContent(BinNum);

	  hGenUpsilonAccNum->Fill(UpsilonPt);
	  if(WeightVal==0){
	    hGenUpsilonAccNumCorr->Fill(UpsilonPt,0);
	  }
	  else {
	    hGenUpsilonAccNumCorr->Fill(UpsilonPt,(1/WeightVal));
	  }
	  
	  
	}
    }

  //Calculate Acceptance by dividing histograms
  TH1F *hGenUpsilonAccUneven = new TH1F ("hGenUpsilonAccUneven", "Generated Upsilon Pt Even Acceptance", PtnBins, PtBins); hGenUpsilonAccUneven->Divide(hGenUpsilonAccNum,hGenUpsilonAccDen);

  //Calculate ratio by dividing histograms
  TH1F *hGenUpsilonAccRatio = new TH1F ("hGenUpsilonAccRatio", "Acceptance Fine Binning", PtnBins, PtBins); hGenUpsilonAccRatio->Divide(hGenUpsilonAccNumCorr,hGenUpsilonAccDen);


  //Set histograms parametrs and write them on root file
  TFile *fOut = new TFile("OutAcc2D_FineBinning.root","RECREATE");


  hGenUpsilonAccUneven->GetYaxis()->SetRangeUser(0.0,1.0);
  hGenUpsilonAccUneven->Write();

  
  hGenUpsilonAccRatio->SetStats(0);
  hGenUpsilonAccRatio->SetTitle("Acceptance 2D Fine Binning");
  hGenUpsilonAccRatio->SetTitleSize(1);
  hGenUpsilonAccRatio->SetLineColor(kBlue);
  hGenUpsilonAccRatio->SetLineWidth(3);
  hGenUpsilonAccRatio->GetYaxis()->SetRangeUser(0.96,1.04);
  hGenUpsilonAccRatio->GetYaxis()->SetTitle("Ratio");
  hGenUpsilonAccRatio->GetYaxis()->SetTitleSize(0.05);
  hGenUpsilonAccRatio->GetYaxis()->SetTitleOffset(1);
  hGenUpsilonAccRatio->GetYaxis()->SetLabelSize(0.035);
  hGenUpsilonAccRatio->GetXaxis()->SetTitle("Upsilon pT (GeV)");
  hGenUpsilonAccRatio->GetXaxis()->SetTitleSize(0.05);
  hGenUpsilonAccRatio->GetXaxis()->SetTitleOffset(1);
  hGenUpsilonAccRatio->GetXaxis()->SetLabelSize(0.035);
  hGenUpsilonAccRatio->Draw();
  hGenUpsilonAccRatio->Write();


  //Draw results on canvas
  /*
  //Set canvas dimensions
  Double_t w = 800; //width
  Double_t h = 800; //height
  
  //Define new Canvas
  TCanvas *cIntS1 = new TCanvas("cIntS1","cIntS1",w,h);
  gStyle->SetOptStat(0);
  cIntS1->Range(0,0,1,1);
  cIntS1->SetFillColor(0);
  cIntS1->SetBorderSize(1);
  cIntS1->SetBorderMode(1);
  cIntS1->SetFrameBorderMode(1);
  cIntS1->SetFrameLineWidth(0);
  cIntS1->SetTopMargin(0.1);
  cIntS1->SetBottomMargin(0.125);
  cIntS1->SetLeftMargin(0.125);

  cIntS1->SetGridx(0);
  cIntS1->SetGridy(1);

  hGenUpsilonPtRatio->SetStats(0);
  hGenUpsilonPtRatio->SetTitle("Corrected odd Acceptance Numerator over odd Acceptance Denominator");
  hGenUpsilonPtRatio->SetTitle("Acceptance 2D Fine Binning 2");
  hGenUpsilonPtRatio->SetTitleSize(1);
  hGenUpsilonPtRatio->SetLineColor(kBlue);
  hGenUpsilonPtRatio->SetLineWidth(3);
  hGenUpsilonPtRatio->GetYaxis()->SetRangeUser(0.96,1.04);
  hGenUpsilonPtRatio->GetYaxis()->SetTitle("Ratio");
  hGenUpsilonPtRatio->GetYaxis()->SetTitleSize(0.05);
  hGenUpsilonPtRatio->GetYaxis()->SetTitleOffset(1);
  hGenUpsilonPtRatio->GetYaxis()->SetLabelSize(0.035);
  hGenUpsilonPtRatio->GetXaxis()->SetTitle("Upsilon pT (GeV)");
  hGenUpsilonPtRatio->GetXaxis()->SetTitleSize(0.05);
  hGenUpsilonPtRatio->GetXaxis()->SetTitleOffset(1);
  hGenUpsilonPtRatio->GetXaxis()->SetLabelSize(0.035);
  hGenUpsilonPtRatio->Draw();

  //Line in x axis along y = 1
  TLine *line2 = new TLine(0.0,1.0,50.0,1.0);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(1);
  line2->SetLineWidth(1);
  //line->Draw();
  */
}
