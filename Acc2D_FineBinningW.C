///////////////////////////////////////////////////////////////////////////
//The purpose of this program is to do the weighted fine binning procedure
//on the 2D Acceptance using the Acceptance spectrums histograms calculated
//on Acc2D.C
//
//To run, type on console:
//.L Acc2D_FineBinningW.C
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

  //Reading input files
  TFile* InputFileAcc = new TFile("OutAcc2D.root","READ");
  TFile* InputFileEff = new TFile("OutEff2D.root","READ");


  //Loading histograms from TFiles objects
  TH2F* hGenUpsilonAcc = (TH2F*)InputFileAcc->Get("hGenUpsilonAcc");
  TH2F* hGenUpsilonAccW = (TH2F*)InputFileAcc->Get("hGenUpsilonAccW");
  TH2F* hGenUpsilonEff = (TH2F*)InputFileEff->Get("hGenUpsilonEff");
  TH2F* hGenUpsilonEffW = (TH2F*)InputFileEff->Get("hGenUpsilonEffW");

  //Defining pt binning
  Float_t PtBins[8] = {0.0,2.0,4.0,6.0,9.0,12.0,30.0,50.0};
  Int_t PtnBins = 7;

  //Defining Acceptance Numerator Histogram
  TH1F* hGenUpsilonAccNum = new TH1F ("hGenUpsilonAccNum", "Acceptance Numerator",PtnBins, PtBins); hGenUpsilonAccNum->Sumw2();
  //Defining Corrected Acceptance Numerator Histogram
  TH1F* hGenUpsilonAccNumCorr = new TH1F ("hGenUpsilonAccNumCorr", "Acceptance Corrected Numerator",PtnBins, PtBins); hGenUpsilonAccNumCorr->Sumw2();


   if (fChain == 0) return;

   

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)//Loop on ALL events
    //for (Long64_t jentry=0; jentry<100000;jentry++)//Lopp on first 10000 events
    {
      
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

	  if( abs(UpsilonRapidity) > 2.4) continue;

	  if( ((MuPlPtMin==false) || (MuMiPtMin==false)) || ((MuPlEtaMax==false) || (MuMiEtaMax==false)) ) continue;
	  //Check for conditions and fill numerator histograms

	  //Get weight values using unweighted and weighted Acceptance and
	  //Efficiency histograms
	  Int_t BinNumAcc = hGenUpsilonAcc->FindBin(UpsilonPt,UpsilonEta);
	  Int_t BinNumAccW = hGenUpsilonAccW->FindBin(UpsilonPt,UpsilonEta);
	  Int_t BinNumEff = hGenUpsilonEff->FindBin(UpsilonPt,UpsilonEta);
	  Int_t BinNumEffW = hGenUpsilonEffW->FindBin(UpsilonPt,UpsilonEta);
	  Double_t AccVal = hGenUpsilonAcc->GetBinContent(BinNumAcc);
	  Double_t AccValW = hGenUpsilonAccW->GetBinContent(BinNumAccW);
	  Double_t EffVal = hGenUpsilonEff->GetBinContent(BinNumEff);
	  Double_t EffValW = hGenUpsilonEffW->GetBinContent(BinNumEffW);
	  
	  hGenUpsilonAccNum->Fill(UpsilonPt,AccVal*EffVal);
	  hGenUpsilonAccNumCorr->Fill(UpsilonPt,AccValW*EffValW);

	}
    }

  //Calculate ratio between histograms by dividing histograms
  TH1F *hGenUpsilonAccNumRatio  = new TH1F(*hGenUpsilonAccNumCorr); hGenUpsilonAccNumRatio->Divide(hGenUpsilonAccNum);

  //Set histograms parameters and write them out on root file
  TFile *fOut = new TFile("OutAcc2D_FineBinningW.root","RECREATE");

  
  hGenUpsilonAccNumRatio->SetStats(0);
  hGenUpsilonAccNumRatio->SetTitle("Acceptance 2D Weighted Fine Binning");
  hGenUpsilonAccNumRatio->SetName("hGenUpsilonAccNumRatio ");
  hGenUpsilonAccNumRatio->SetTitleSize(1);
  hGenUpsilonAccNumRatio->SetLineColor(kBlue);
  hGenUpsilonAccNumRatio->SetLineWidth(3);
  hGenUpsilonAccNumRatio->GetYaxis()->SetRangeUser(0.96,1.04);
  hGenUpsilonAccNumRatio->GetYaxis()->SetTitleSize(0.05);
  hGenUpsilonAccNumRatio->GetYaxis()->SetTitleOffset(1);
  hGenUpsilonAccNumRatio->GetYaxis()->SetLabelSize(0.035);
  hGenUpsilonAccNumRatio->GetYaxis()->SetTitle("Ratio");
  hGenUpsilonAccNumRatio->GetXaxis()->SetTitle("Upsilon pT (GeV)");
  hGenUpsilonAccNumRatio->GetXaxis()->SetTitleSize(0.05);
  hGenUpsilonAccNumRatio->GetXaxis()->SetTitleOffset(1);
  hGenUpsilonAccNumRatio->GetXaxis()->SetLabelSize(0.035);
   
  hGenUpsilonAccNumRatio ->Write();

  fOut->Close();

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

  hGenUpsilonAccNumRatio->SetStats(0);
  hGenUpsilonAccNumRatio->SetTitle("Acceptance 2D Weighted Fine Binning");
  hGenUpsilonAccNumRatio->SetTitleSize(1);
  hGenUpsilonAccNumRatio->SetLineColor(kBlue);
  hGenUpsilonAccNumRatio->SetLineWidth(3);
  hGenUpsilonAccNumRatio->GetYaxis()->SetRangeUser(0.96,1.04);
  hGenUpsilonAccNumRatio->GetYaxis()->SetTitle("Ratio");
  hGenUpsilonAccNumRatio->GetYaxis()->SetTitleSize(0.05);
  hGenUpsilonAccNumRatio->GetYaxis()->SetTitleOffset(1);
  hGenUpsilonAccNumRatio->GetYaxis()->SetLabelSize(0.035);
  hGenUpsilonAccNumRatio->GetXaxis()->SetTitle("Upsilon pT (GeV)");
  hGenUpsilonAccNumRatio->GetXaxis()->SetTitleSize(0.05);
  hGenUpsilonAccNumRatio->GetXaxis()->SetTitleOffset(1);
  hGenUpsilonAccNumRatio->GetXaxis()->SetLabelSize(0.035);
  hGenUpsilonAccNumRatio->Draw();

  //Line in x axis along y = 1
  TLine *line2 = new TLine(0.0,1.0,50.0,1.0);
  line2->SetLineColor(kBlack);
  line2->SetLineStyle(1);
  line2->SetLineWidth(1);
  //line->Draw();
  */
}
