///////////////////////////////////////////////////////////////////////////
//The purpose of this program is to calculate the 2-Dimensional
//Upsilon y vs Upsilon pT Acceptance spectrum for uneven events and use
//the even events spectrum calculated at Eff2D_even to do the fine binning
//process.
//
//To run, type on console:
//.L Eff2D_FineBinning.C
//Efficiency t
//t.Loop()
///////////////////////////////////////////////////////////////////////////

//defining headers
#define Efficiency_cxx
#include "Efficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//Defining program
void Efficiency::Loop()
{

  //reading inputfile
  TFile* InputFile = new TFile("OutEff2Deven.root","READ");

  //Loading histograms from TFiles objects
  TH2F* hGenUpsilonEffEven = (TH2F*)InputFile->Get("hGenUpsilonEffEven");
  
  //Defining pt binning
  Float_t PtBins[29] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.75,7.5,8.25,9.0,9.75,10.5,11.25,12.0,16.5,21.0,25.5,30.0,35.0,40.0,45.0,50.0};
  Int_t PtnBins= 28;


  //Defining histograms

  //Defining Efficiency Denominator Histogram
  TH1F* hGenUpsilonEffDen = new TH1F ("hGenUpsilonEffDen", "1D Efficiency denominator for uneven events",PtnBins, PtBins); hGenUpsilonEffDen->Sumw2();
  //Defining Efficiency numerator histogram
  TH1F* hGenUpsilonEffNum = new TH1F ("hGenUpsilonEffNum", "1D Efficiency numerator for uneven events",PtnBins, PtBins); hGenUpsilonEffNum->Sumw2();
  //Defining Efficiency corrected numerator histogram
  TH1F* hGenUpsilonEffNumCorr = new TH1F ("hGenUpsilonEffNumCorr", "1D Efficiency Corrected numerator for uneven events",PtnBins, PtBins); hGenUpsilonEffNumCorr->Sumw2();
  


  if (fChain == 0) return;


  //Declaring Tree branches
  TFile *File = TFile::Open("/home/danielg/Documents/Tesis/Efficiency.root");
  TTree *MyTree = (TTree*)File->Get("hionia/myTree");
  int nevents= (int)MyTree->GetEntries();
  cout<<"nevents = "<<nevents<<"\n";
   
  TBranch *b_Gen_QQ_size = MyTree->GetBranch("Gen_QQ_size");
  b_Gen_QQ_size->SetAddress(&Gen_QQ_size);

  TBranch *b_Gen_QQ_mupl_idx = MyTree->GetBranch("Gen_QQ_mupl_idx");
  b_Gen_QQ_mupl_idx->SetAddress(&Gen_QQ_mupl_idx);

  TBranch *b_Gen_QQ_mumi_idx = MyTree->GetBranch("Gen_QQ_mumi_idx");
  b_Gen_QQ_mumi_idx->SetAddress(&Gen_QQ_mumi_idx);

  TBranch *b_Gen_QQ_whichRec = MyTree->GetBranch("Gen_QQ_whichRec");
  b_Gen_QQ_whichRec->SetAddress(&Gen_QQ_whichRec);

  TBranch *b_Reco_QQ_size = MyTree->GetBranch("Reco_QQ_size");
  b_Reco_QQ_size->SetAddress(&Reco_QQ_size);

  TBranch *b_Reco_QQ_mupl_idx = MyTree->GetBranch("Reco_QQ_mupl_idx");
  b_Reco_QQ_mupl_idx->SetAddress(&Reco_QQ_mupl_idx);

  TBranch *b_Reco_QQ_mumi_idx = MyTree->GetBranch("Reco_QQ_mumi_idx");
  b_Reco_QQ_mumi_idx->SetAddress(&Reco_QQ_mumi_idx);

  TBranch *b_Reco_mu_SelectionType = MyTree->GetBranch("Reco_mu_SelectionType");
  b_Reco_mu_SelectionType->SetAddress(&Reco_mu_SelectionType);

  TBranch *b_Reco_mu_nTrkWMea = MyTree->GetBranch("Reco_mu_nTrkWMea");
  b_Reco_mu_nTrkWMea->SetAddress(&Reco_mu_nTrkWMea);

  TBranch *b_Reco_mu_nPixWMea = MyTree->GetBranch("Reco_mu_nPixWMea");
  b_Reco_mu_nPixWMea->SetAddress(&Reco_mu_nPixWMea);

  TBranch *b_Reco_mu_dxy = MyTree->GetBranch("Reco_mu_dxy");
  b_Reco_mu_dxy->SetAddress(&Reco_mu_dxy);

  TBranch *b_Reco_mu_dz = MyTree->GetBranch("Reco_mu_dz");
  b_Reco_mu_dz->SetAddress(&Reco_mu_dz);

  TBranch *b_Reco_QQ_sign = MyTree->GetBranch("Reco_QQ_sign");
  b_Reco_QQ_sign->SetAddress(&Reco_QQ_sign);

  TBranch *b_Reco_QQ_trig = MyTree->GetBranch("Reco_QQ_trig");
  b_Reco_QQ_trig->SetAddress(&Reco_QQ_trig);

  TBranch *b_Reco_mu_trig = MyTree->GetBranch("Reco_mu_trig");
  b_Reco_mu_trig->SetAddress(&Reco_mu_trig);

  TBranch *b_HLTriggers = MyTree->GetBranch("HLTriggers");
  b_HLTriggers->SetAddress(&HLTriggers);

  TBranch *b_Reco_QQ_VtxProb = MyTree->GetBranch("Reco_QQ_VtxProb");
  b_Reco_QQ_VtxProb->SetAddress(&Reco_QQ_VtxProb);

  //4-moms
  TBranch *b_Gen_QQ_4mom;
  TClonesArray *Gen_QQ_4mom = new TClonesArray();
  MyTree->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   
  TBranch *b_Gen_mu_4mom;
  TClonesArray *Gen_mu_4mom = new TClonesArray();
  MyTree->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);

  TBranch *b_Reco_QQ_4mom;
  TClonesArray *Reco_QQ_4mom = new TClonesArray();
  MyTree->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);

  TBranch *b_Reco_mu_4mom;
  TClonesArray *Reco_mu_4mom = new TClonesArray();
  MyTree->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++)//Loop on ALL events
    //for (Long64_t jentry=0; jentry<100000;jentry++)//Loop on first 10000 events
    {
      //check if event is even; if not, read event
      if( jentry % 2 == 0) continue;
      
      cout<<"Reading event "<<jentry<<endl;
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      //MyTree->GetEntry(jEntry);
      if(jentry%nentries==0){ cout<<"Part of tree scanned : "<< 100*(double)jentry / (double)nentries<<" % \n";}

      Gen_QQ_4mom->Clear();
      Gen_mu_4mom->Clear();
      Reco_QQ_4mom->Clear();
      Reco_mu_4mom->Clear();

      b_Gen_QQ_size->GetEntry(jentry);
      if(Gen_QQ_size > 0)
	{
	  b_Gen_QQ_mupl_idx->GetEntry(jentry);
	  b_Gen_QQ_mumi_idx->GetEntry(jentry);
	  b_Gen_QQ_4mom->GetEntry(jentry);
	  b_Gen_mu_4mom->GetEntry(jentry);
	  b_Gen_QQ_whichRec->GetEntry(jentry);
	}

      b_Reco_QQ_size->GetEntry(jentry);
      if(Reco_QQ_size > 0)
	{
	  b_Reco_QQ_mupl_idx->GetEntry(jentry);
	  b_Reco_QQ_mumi_idx->GetEntry(jentry);
	  b_Reco_QQ_4mom->GetEntry(jentry);
	  b_Reco_mu_4mom->GetEntry(jentry);
	  b_Reco_mu_SelectionType->GetEntry(jentry);
	  b_Reco_mu_whichGen->GetEntry(jentry);
	  b_Reco_mu_nTrkWMea->GetEntry(jentry);
	  b_Reco_mu_nPixWMea->GetEntry(jentry);
	  b_Reco_mu_dxy->GetEntry(jentry);
	  b_Reco_mu_dz->GetEntry(jentry);
	  b_Reco_QQ_sign->GetEntry(jentry);
	  b_Reco_QQ_trig->GetEntry(jentry);
	  b_Reco_mu_trig->GetEntry(jentry);
	  b_HLTriggers->GetEntry(jentry);
	  b_Reco_QQ_VtxProb->GetEntry(jentry);
	}

      for(int iQQ=0; iQQ < Gen_QQ_size; iQQ++)//Loop on generated upsilon
	{
	 TLorentzVector *GenQQ = (TLorentzVector*) Gen_QQ_4mom->At(iQQ);
	 float GenUpsilonPt = GenQQ->Pt();
	 float GenUpsilonRapidity = GenQQ->Rapidity();
	 float GenUpsilonEta = GenQQ->Eta();

	 if( GenUpsilonPt > 50.0 ) continue;
	 
	 TLorentzVector *GenQQmuPos = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mupl_idx[iQQ]);
	 TLorentzVector *GenQQmuNeg = (TLorentzVector*) Gen_mu_4mom->At(Gen_QQ_mumi_idx[iQQ]);
	 float MuNegGenPt = GenQQmuNeg->Pt();
	 float MuPosGenPt = GenQQmuPos->Pt();
	 
	 float MuNegGenEta = GenQQmuNeg->Eta();
	 float MuPosGenEta = GenQQmuPos->Eta();

	 //Defining EFficiency denominator conditions as booleans
	 bool GenUpsilonRapidityMax = abs(GenUpsilonRapidity) < 2.4;

	 bool MuPosGenPtMin = MuPosGenPt > 3.5;
	 bool MuPosGenEtaMax = abs(MuPosGenEta) < 2.4;

	 bool MuNegGenPtMin = MuNegGenPt > 3.5;
	 bool MuNegGenEtaMax = abs(MuNegGenEta) < 2.4;

	 bool MuonPtMin = ( MuPosGenPtMin && MuNegGenPtMin );
	 bool MuonEtaMax = ( MuPosGenEtaMax && MuNegGenEtaMax );

	 //Check for conditions and fill denominator histograms
	 if( (GenUpsilonRapidityMax == false) || (MuonPtMin == false) || (MuonEtaMax == false) ) continue;

	 hGenUpsilonEffDen->Fill(GenUpsilonPt);

	 //Check reco event
	 int isReco = Gen_QQ_whichRec[iQQ];

	 if( isReco < 0 ) continue;

	 int iMuPos = Reco_QQ_mupl_idx[isReco];
	 int iMuNeg = Reco_QQ_mumi_idx[isReco];

	 int isMuGenPos = Reco_mu_whichGen[iMuPos];
	 int isMuGenNeg = Reco_mu_whichGen[iMuNeg];

	 if( (isMuGenPos < 0) || (isMuGenNeg < 0) ) continue;

	 TLorentzVector *RecoQQ_MuPos = (TLorentzVector*) Reco_mu_4mom->At(iMuPos);
	   TLorentzVector *RecoQQ_MuNeg = (TLorentzVector*) Reco_mu_4mom->At(iMuNeg);

	   float RecoMuPosPt = RecoQQ_MuPos->Pt();
	   float RecoMuPosEta = RecoQQ_MuPos->Eta();
	   float RecoMuNegPt = RecoQQ_MuNeg->Pt();
	   float RecoMuNegEta = RecoQQ_MuNeg->Eta();
	   
	   //Defining muon numerator conditions as booleans
	   
	   //Defining pT conditions as booleans
	   bool RecoMuPosPtMin = RecoMuPosPt > 3.5;
	   bool RecoMuNegPtMin = RecoMuNegPt > 3.5;
	   if( (RecoMuPosPtMin == false) || (RecoMuNegPtMin == false) ) continue;

	   //Defining Eta conditions as booleans
	   bool RecoMuPosEtaMin = abs(RecoMuPosEta) < 2.4;
	   bool RecoMuNegEtaMin = abs(RecoMuNegEta) < 2.4;
	   if( (RecoMuPosEtaMin == false) || (RecoMuNegEtaMin == false) ) continue;
	   
	   //Check if both muons are Global Muons
	   bool IsGlobalMuPos = (Reco_mu_SelectionType[iMuPos]&((ULong64_t)pow(2, 1))) > 0;
	   bool IsGlobalMuNeg = (Reco_mu_SelectionType[iMuNeg]&((ULong64_t)pow(2, 1))) > 0;
	   if ( (IsGlobalMuPos == false) || (IsGlobalMuNeg == false) ) continue;
	   
	   //Check if both muons are Tracker Muons
	   bool IsTrackerMuPos = (Reco_mu_SelectionType[iMuPos]&((ULong64_t)pow(2, 3))) > 0;
	   bool IsTrackerMuNeg = (Reco_mu_SelectionType[iMuNeg]&((ULong64_t)pow(2, 3))) > 0;
	   if ( (IsTrackerMuPos == false) || (IsTrackerMuNeg == false) ) continue;

	   //Check if both muons have at least 6 silicon tracker hits
	   bool nTrkWMeaMinPos = ( Reco_mu_nTrkWMea[iMuPos]  > 5 );
	   bool nTrkWMeaMinNeg = ( Reco_mu_nTrkWMea[iMuNeg] > 5 );
	   if( (nTrkWMeaMinPos == false) || (nTrkWMeaMinNeg == false) ) continue;

	   //Check if both muons have at least 1 silicon pixel hit
	   bool nPixWMeaMinPos = ( Reco_mu_nPixWMea[iMuPos] > 0 );
	   bool nPixWMeaMinNeg = ( Reco_mu_nPixWMea[iMuNeg] > 0 );
	   if( (nPixWMeaMinPos == false) || (nPixWMeaMinNeg == false) ) continue;

	   //Check if both muons' distance from the muon track to the
	   //closest primary vertex is less than 3mm in the transverse direction
	   bool dxyMuPos = ( abs(Reco_mu_dxy[iMuPos]) < 0.3 );
	   bool dxyMuNeg = ( abs(Reco_mu_dxy[iMuNeg]) < 0.3 );
	   if( (dxyMuPos == false) || (dxyMuNeg == false) ) continue;

	   //Check if both muons' distance from the muon track to the
	   //closest primary vertex is less than 20cm in the longitudinal dir.
	   bool dzMuPos = ( abs(Reco_mu_dz[iMuPos]) < 20.0 );
	   bool dzMuNeg = ( abs(Reco_mu_dz[iMuNeg]) < 20.0 );
	   if( ( dzMuPos == false) || (dzMuNeg == false) ) continue;

	   //Check Reconstructed Upsilon
	   TLorentzVector *RecoQQ = (TLorentzVector*) Reco_QQ_4mom->At(isReco);

	   float RecoQQPt = RecoQQ->Pt();
	   float RecoQQRap = RecoQQ->Rapidity();
	   float RecoQQMass = RecoQQ->M();
	   float RecoQQVtxProb = Reco_QQ_VtxProb[isReco];
	   int RecoQQCharge = Reco_QQ_sign[isReco];

	   //Defining Reconstructed Upsilon conditions as a boolean
	   bool RecoQQAcc = ( (RecoQQPt < 100.0) && (RecoQQVtxProb > 0.01) && (RecoQQMass> 8.5) && (RecoQQMass < 11.0) && (RecoQQCharge == 0) );
	   if( RecoQQAcc == false ) continue;

	   int TriggerBit = 3;
	   //Define trigger conditions as booleans
	   bool JpsiHLTtrigers = (HLTriggers&((ULong64_t)pow(2, TriggerBit))) > 0;
	   bool JpsiTrigger = (Reco_QQ_trig[isReco]&((ULong64_t)pow(2, TriggerBit))) > 0;
	   if(JpsiHLTtrigers == false || JpsiTrigger == false) continue;

	   //Get weight value using the corresponding bin value from the
	   //Acceptance even events spectrum
	   Int_t BinNum = hGenUpsilonEffEven->FindBin(GenUpsilonPt,GenUpsilonEta);
	   Double_t WeightVal = hGenUpsilonEffEven->GetBinContent(BinNum);

	   //Check for conditions and fill numerator histograms
	   hGenUpsilonEffNum->Fill(GenUpsilonPt);
	   if(WeightVal==0){
	     hGenUpsilonEffNumCorr->Fill(GenUpsilonPt,0);
	   }
	   else {
	     hGenUpsilonEffNumCorr->Fill(GenUpsilonPt,(1/WeightVal));
	   }
	 
	 
	}
    }

  //Calculate Efficiency  by dividing histograms
  TH1F *hGenUpsilonEffUneven = new TH1F ("hGenUpsilonEffUneven", "Generated Upsilon Pt Even Acceptance", PtnBins, PtBins); hGenUpsilonEffUneven->Divide(hGenUpsilonEffNum,hGenUpsilonEffDen);

  //Calculate ratio by dividing histograms
  TH1F *hGenUpsilonEffRatio = new TH1F ("hGenUpsilonEffRatio", "Efficiency Fine Binning", PtnBins, PtBins); hGenUpsilonEffRatio->Divide(hGenUpsilonEffNumCorr,hGenUpsilonEffDen);

  
  
  TFile *fOut = new TFile("OutEff2D_FineBinning.root","RECREATE");

  
  hGenUpsilonEffUneven->GetYaxis()->SetRangeUser(0.0,1.0);
  hGenUpsilonEffUneven->Write();

  
  hGenUpsilonEffRatio->SetStats(0);
  hGenUpsilonEffRatio->SetTitle("Efficiency 2D Fine Binning");
  hGenUpsilonEffRatio->SetTitleSize(1);
  hGenUpsilonEffRatio->SetLineColor(kBlue);
  hGenUpsilonEffRatio->SetLineWidth(3);
  hGenUpsilonEffRatio->GetYaxis()->SetRangeUser(0.96,1.04);
  hGenUpsilonEffRatio->GetYaxis()->SetTitle("Ratio");
  hGenUpsilonEffRatio->GetYaxis()->SetTitleSize(0.05);
  hGenUpsilonEffRatio->GetYaxis()->SetTitleOffset(1);
  hGenUpsilonEffRatio->GetYaxis()->SetLabelSize(0.035);
  hGenUpsilonEffRatio->GetXaxis()->SetTitle("Upsilon pT (GeV)");
  hGenUpsilonEffRatio->GetXaxis()->SetTitleSize(0.05);
  hGenUpsilonEffRatio->GetXaxis()->SetTitleOffset(1);
  hGenUpsilonEffRatio->GetXaxis()->SetLabelSize(0.035);
  hGenUpsilonEffRatio->Draw();
  hGenUpsilonEffRatio->Write();


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
  hGenUpsilonPtRatio->SetTitle("Efficiency 2D Fine Binning 3");
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
