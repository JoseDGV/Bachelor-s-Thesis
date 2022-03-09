///////////////////////////////////////////////////////////////////////////
//The purpose of this program is to calculate the 2-Dimensional
//Upsilon y vs Upsilon pT Acceptance spectrum for even events to use in
//the fine binning process.
//
//To run, type on console:
//.L Acc2D_even.C
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
  
  if (fChain == 0) return;

  ///////////////////////////////////////////////////////////////////////
  //Old pT binning used for fine binning
  //Float_t PtBins1[8] = {0.0,2.0,4.0,6.0,9.0,12.0,30.0,50.0};
  //Int_t PtnBins1 = 7;

  //Float_t PtBins2[15] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.5,9.0,10.5,12.0,21.0,30.0,40.0,50.0};//,60.0,100.0};
  //Int_t PtnBins2 = 14;
  ////////////////////////////////////////////////////////////////////////
  
  //Defining pt binning
  Float_t PtBins[29] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.75,7.5,8.25,9.0,9.75,10.5,11.25,12.0,16.5,21.0,25.5,30.0,35.0,40.0,45.0,50.0};
  Int_t PtnBins= 28;

  //Defining rapidity binning
  Float_t RapBins[25] = {-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
  Int_t RapnBins = 24;
  

  //Defining Acceptance Denominator Histogram
  TH2F *hGenUpsilonAccDen = new TH2F ("hGenUpsilonAccDen", "2D Acceptance Denominator for even events", PtnBins, PtBins, RapnBins, RapBins);
  //Defining Acceptance Numerator Histogram
  TH2F *hGenUpsilonAccNum = new TH2F ("hGenUpsilonAccNum", "2D Acceptance Numerator for even events", PtnBins, PtBins, RapnBins, RapBins);
   
  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++)//Loop on ALL events
    //for (Long64_t jentry=0; jentry<100000;jentry++)//Loop on first 10000 events
    {
      //Check if entry is even; if not, skip event
      if( jentry % 2 != 0) continue;
      
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

	  if( abs(UpsilonRapidity) < 2.4) {
	    //Check for conditions and fill denominator histograms
	    hGenUpsilonAccDen->Fill(UpsilonPt,UpsilonRapidity);

	    if( ((MuPlPtMin==true) & (MuMiPtMin==true)) & ((MuPlEtaMax==true) & (MuMiEtaMax==true)) ){
	      //Check for conditions and fill numerator histograms
	      hGenUpsilonAccNum->Fill(UpsilonPt,UpsilonRapidity);
	    }
	  }
	}
    }


  //Calculate Acceptance by dividing histograms
  TH2F *hGenUpsilonAccEven  = new TH2F(*hGenUpsilonAccNum); hGenUpsilonAccEven->Divide(hGenUpsilonAccDen);

  
  TFile *fOut = new TFile("OutAcc2D_even.root","RECREATE");

   hGenUpsilonAccEven->SetStats(0);
   hGenUpsilonAccEven->SetContour(1000);
   hGenUpsilonAccEven->SetTitle("Acceptance Upsilon pT vs y distribution for even events");
   hGenUpsilonAccEven->SetTitleSize(1);
   hGenUpsilonAccEven->SetName("hGenUpsilonAccEven");
   hGenUpsilonAccEven->GetXaxis()->SetTitleSize(0.06);
   hGenUpsilonAccEven->GetXaxis()->SetTitleOffset(0.65);
   hGenUpsilonAccEven->GetXaxis()->SetLabelSize(0.04);
   hGenUpsilonAccEven->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonAccEven->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonAccEven->GetYaxis()->SetTitleSize(0.06);
   hGenUpsilonAccEven->GetYaxis()->SetTitleOffset(0.65);
   hGenUpsilonAccEven->GetYaxis()->SetLabelSize(0.04);
   hGenUpsilonAccEven->SetOption("colz");
   
   hGenUpsilonAccEven->Write();

}
