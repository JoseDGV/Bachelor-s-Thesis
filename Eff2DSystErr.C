///////////////////////////////////////////////////////////////////////////
//The purpose of this program is to calculate the 2-Dimensional
//Upsilon y vs Upsilon pT Efficiency systematic errors corrections.
//
//To run, type on console:
//.L Eff2DSystErr.C
//Efficiency t
//t.Loop()
///////////////////////////////////////////////////////////////////////////

//defining headers
#define Efficiency_cxx
#include "Efficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//systematic errors
#include "tnp_weight.h"

/*
Systematic errors Index
i=1 : central_value
i=2 : statistical error
i=3 : systematic error
i=4 : total error
 */


//Defining function to get Global Muon Tight Acceptance error values
Double_t tnp_weight_GlobalMuon_TightAcceptance_pp_func(Double_t pt,Double_t eta,Int_t i){
 auto SFwithError = tnp_weight_GlobalMuon_TightAcceptance_pp(pt,eta );
 double central_value = std::get<0>(SFwithError);
 double stat_error = std::get<1>(SFwithError);
 double syst_error = std::get<2>(SFwithError);
 double tot_error = std::get<3>(SFwithError);
 if (i==1) return central_value;
 if (i==2) return stat_error;
 if (i==3) return syst_error;
 if (i==4) return tot_error;
}

//Defining function to get Hybrid Hard Soft Muons error values
Double_t tnp_weight_HybridSoftID_LooseAcceptance_pp_func(Double_t pt,Double_t eta,Int_t i){
 auto SFwithError = tnp_weight_HybridSoftID_LooseAcceptance_pp(pt,eta );
 double central_value = std::get<0>(SFwithError);
 double stat_error = std::get<1>(SFwithError);
 double syst_error = std::get<2>(SFwithError);
 double tot_error = std::get<3>(SFwithError);
 if (i==1) return central_value;
 if (i==2) return stat_error;
 if (i==3) return syst_error;
 if (i==4) return tot_error;
}

//Defining weight function
Double_t weightF(Double_t x) {
     //Defining function parameters
     Double_t A = 251.102;
     Double_t B = 18.5373;
     Double_t C = 23.5803;
     Double_t D = -4.97668;

     //Defining function
     Double_t arg = ( A + B*x + C*pow(x,2) )/pow(x-D,3);
     return arg;
   }

//Defining program
void Efficiency::Loop()
{
  //Defining pt binning
  Float_t PtBins[29] = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.75,7.5,8.25,9.0,9.75,10.5,11.25,12.0,16.5,21.0,25.5,30.0,35.0,40.0,45.0,50.0};
  Int_t PtnBins= 28;

  //Defining rapidity binning
  Float_t RapBins[25] = {-2.4,-2.2,-2.0,-1.8,-1.6,-1.4,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4};
  Int_t RapnBins = 24;

  //Defining Efficiency Denominator Histogram
   TH2F *hGenUpsilonEffDen = new TH2F ("hGenUpsilonEffDen", "2D Efficiency Nominal Error Correction Denominator", PtnBins, PtBins, RapnBins, RapBins);
   
  //Defining Efficiency Nominal Correction Numerator Histogram
   TH2F *hGenUpsilonEffNumNominal = new TH2F ("hGenUpsilonEffNumNominal", "2D Efficiency Nominal Error Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);

   //Defining Efficiency Systematic Errors Corrections Histograms
   TH2F *hGenUpsilonEffNumHSPlus1SystErr = new TH2F ("hGenUpsilonEffNumHSPlus1SystErr", "2D Efficiency HSPlus1SystErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   TH2F *hGenUpsilonEffNumHSMinus1SystErr = new TH2F ("hGenUpsilonEffNumHSMinus1SystErr", "2D Efficiency  HSMinus1SystErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   TH2F *hGenUpsilonEffNumHSPlus1StatErr = new TH2F ("hGenUpsilonEffNumHSPlus1StatErr", "2D Efficiency  HSPlus1StatErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   TH2F *hGenUpsilonEffNumHSMinus1StatErr = new TH2F ("hGenUpsilonEffNumHSMinus1StatErr", "2D Efficiency  HSMinus1StatErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);

   TH2F *hGenUpsilonEffNumGMPlus1SystErr = new TH2F ("hGenUpsilonEffNumGMPlus1SystErr", "2D Efficiency GMPlus1SystErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   TH2F *hGenUpsilonEffNumGMMinus1SystErr = new TH2F ("hGenUpsilonEffNumGMMinus1SystErr", "2D Efficiency GMMinus1SystErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   TH2F *hGenUpsilonEffNumGMPlus1StatErr = new TH2F ("hGenUpsilonEffNumGMPlus1StatErr", "2D Efficiency GMPlus1StatErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   TH2F *hGenUpsilonEffNumGMMinus1StatErr = new TH2F ("hGenUpsilonEffNumGMMinus1StatErr", "2D Efficiency GMPlus1StatErr Correction Numerator", PtnBins, PtBins, RapnBins, RapBins);
   
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

	 if(MuPosGenPt>30 || MuNegGenPt>30)continue;
	 
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

	 //Check for conditions and fill denominator histogram
	 if( (GenUpsilonRapidityMax == false) || (MuonPtMin == false) || (MuonEtaMax == false) ) continue;

	 hGenUpsilonEffDen->Fill(GenUpsilonPt,GenUpsilonRapidity);
	 

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

	   //Use functions to get systematic error values
	   //////////SYSTEMATIC ERRORS/////////////////
	   //i=1 : central_value
	   //i=2 : statistical error
	   //i=3 : systematic error
	   //i=4 : total error
	   
	   ///POSITIVE MUON///
	   Double_t MuPosGM_Ctrl_Val = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),1);
	   Double_t MuPosGM_Stat_Err = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),2);
	   Double_t MuPosGM_Syst_Err = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),3);
	   Double_t MuPosGM_Tot_Err = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),4);

	   Double_t MuPosHS_Ctrl_Val = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),1);
	   Double_t MuPosHS_Stat_Err = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),2);
	   Double_t MuPosHS_Syst_Err = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),3);
	   Double_t MuPosHS_Tot_Err = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuPosGenPt,abs(MuPosGenEta),4);
	   cout<<"("<<MuPosGenPt<<","<<MuPosGenEta<<")"<<","<<MuPosGM_Ctrl_Val<<","<<MuPosGM_Stat_Err<<","<<MuPosGM_Syst_Err<<","<<MuPosGM_Tot_Err<<endl;
	   cout<<"("<<MuPosGenPt<<","<<MuPosGenEta<<")"<<","<<MuPosHS_Ctrl_Val<<","<<MuPosHS_Stat_Err<<","<<MuPosHS_Syst_Err<<","<<MuPosHS_Tot_Err<<endl;
	   ///POSITVE MUON///

	   ///NEGATIVE MUON///
	   Double_t MuNegGM_Ctrl_Val = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),1);
	   Double_t MuNegGM_Stat_Err = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),2);
	   Double_t MuNegGM_Syst_Err = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),3);
	   Double_t MuNegGM_Tot_Err = tnp_weight_GlobalMuon_TightAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),4);

	   Double_t MuNegHS_Ctrl_Val = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),1);
	   Double_t MuNegHS_Stat_Err = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),2);
	   Double_t MuNegHS_Syst_Err = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),3);
	   Double_t MuNegHS_Tot_Err = tnp_weight_HybridSoftID_LooseAcceptance_pp_func(MuNegGenPt,abs(MuNegGenEta),4);
	   cout<<"("<<MuNegGenPt<<","<<MuNegGenEta<<")"<<","<<MuNegGM_Ctrl_Val<<","<<MuNegGM_Stat_Err<<","<<MuNegGM_Syst_Err<<","<<MuNegGM_Tot_Err<<endl;
	   cout<<"("<<MuNegGenPt<<","<<MuNegGenEta<<")"<<","<<MuNegHS_Ctrl_Val<<","<<MuNegHS_Stat_Err<<","<<MuNegHS_Syst_Err<<","<<MuNegHS_Tot_Err<<endl;
	   ///NEGATIVE MUON///

	   //Calculate systematic errors corrections

	   //Nominal
	   Double_t tnp_weight_nominal= (MuPosGM_Ctrl_Val * MuPosHS_Ctrl_Val)*(MuNegGM_Ctrl_Val*MuNegHS_Ctrl_Val);
	   hGenUpsilonEffNumNominal->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_nominal);

	   //HYBRID SOFT ID//

	   //Plus 1 systematic error
	   Double_t tnp_weight_HSplus1syst = (MuPosGM_Ctrl_Val * (MuPosHS_Ctrl_Val+MuPosHS_Syst_Err)) * (MuNegGM_Ctrl_Val * (MuNegHS_Ctrl_Val+MuNegHS_Syst_Err));
	   hGenUpsilonEffNumHSPlus1SystErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_HSplus1syst);

	   //Minus 1 systematic error
	   Double_t tnp_weight_HSminus1syst = (MuPosGM_Ctrl_Val * (MuPosHS_Ctrl_Val-MuPosHS_Syst_Err)) * (MuNegGM_Ctrl_Val * (MuNegHS_Ctrl_Val-MuNegHS_Syst_Err));
	   hGenUpsilonEffNumHSMinus1SystErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_HSminus1syst);

	   //Plus 1 statistical error
	   Double_t tnp_weight_HSplus1stat = (MuPosGM_Ctrl_Val * (MuPosHS_Ctrl_Val+MuPosHS_Stat_Err)) * (MuNegGM_Ctrl_Val * (MuNegHS_Ctrl_Val+MuNegHS_Stat_Err));
	   hGenUpsilonEffNumHSPlus1StatErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_HSplus1stat);

	   //Minus 1 statistical error
	   Double_t tnp_weight_HSminus1stat = (MuPosGM_Ctrl_Val * (MuPosHS_Ctrl_Val-MuPosHS_Stat_Err)) * (MuNegGM_Ctrl_Val * (MuNegHS_Ctrl_Val-MuNegHS_Stat_Err));
	   hGenUpsilonEffNumHSMinus1StatErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_HSminus1stat);
	    //HYBRID SOFT ID//

	   //GLOBAL MUON//

	   //Plus 1 systematic error
	   Double_t tnp_weight_GMplus1syst = ((MuPosGM_Ctrl_Val+MuPosGM_Syst_Err) * MuPosHS_Ctrl_Val) * ((MuNegGM_Ctrl_Val+MuNegGM_Syst_Err) * MuNegHS_Ctrl_Val);
	   hGenUpsilonEffNumGMPlus1SystErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_GMplus1syst);

	   //Minus 1 systematic error
	   Double_t tnp_weight_GMminus1syst = ((MuPosGM_Ctrl_Val-MuPosGM_Syst_Err) * MuPosHS_Ctrl_Val) * ((MuNegGM_Ctrl_Val-MuNegGM_Syst_Err) * MuNegHS_Ctrl_Val);
	   hGenUpsilonEffNumGMMinus1SystErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_GMminus1syst);

	   //Plus 1 statistical error
	   Double_t tnp_weight_GMplus1stat = ((MuPosGM_Ctrl_Val+MuPosGM_Stat_Err) * MuPosHS_Ctrl_Val) * ((MuNegGM_Ctrl_Val+MuNegGM_Stat_Err) * MuNegHS_Ctrl_Val);
	   hGenUpsilonEffNumGMPlus1StatErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_GMplus1stat);

	   //Minus 1 statitical error
	   Double_t tnp_weight_GMminus1stat = ((MuPosGM_Ctrl_Val-MuPosGM_Stat_Err) * MuPosHS_Ctrl_Val) * ((MuNegGM_Ctrl_Val-MuNegGM_Stat_Err) * MuNegHS_Ctrl_Val);
	   hGenUpsilonEffNumGMMinus1StatErr->Fill(GenUpsilonPt,GenUpsilonRapidity,tnp_weight_GMminus1stat);
	    //GLOBAL MUON//

	   //Print error values
	   cout<<tnp_weight_nominal<<","<<tnp_weight_HSplus1syst<<","<<tnp_weight_HSminus1syst<<","<<tnp_weight_HSplus1stat<<","<<tnp_weight_HSminus1stat<<endl;
	   cout<<tnp_weight_nominal<<","<<tnp_weight_GMplus1syst<<","<<tnp_weight_GMminus1syst<<","<<tnp_weight_GMplus1stat<<","<<tnp_weight_GMminus1stat<<endl;
	   
	   /////////SYSTEMATIC ERRORS////////////////
	   
           
	 
	 
	}
    }
  
  //Calculate Efficiency by dividing histograms
  TH2F *hGenUpsilonEffNominal  = new TH2F(*hGenUpsilonEffNumNominal); hGenUpsilonEffNominal->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffHSPlus1SystErr  = new TH2F(*hGenUpsilonEffNumHSPlus1SystErr); hGenUpsilonEffHSPlus1SystErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffHSMinus1SystErr  = new TH2F(*hGenUpsilonEffNumHSMinus1SystErr); hGenUpsilonEffHSMinus1SystErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffHSPlus1StatErr  = new TH2F(*hGenUpsilonEffNumHSPlus1StatErr); hGenUpsilonEffHSPlus1StatErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffHSMinus1StatErr  = new TH2F(*hGenUpsilonEffNumHSMinus1StatErr); hGenUpsilonEffHSMinus1StatErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffGMPlus1SystErr  = new TH2F(*hGenUpsilonEffNumGMPlus1SystErr); hGenUpsilonEffGMPlus1SystErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffGMMinus1SystErr  = new TH2F(*hGenUpsilonEffNumGMMinus1SystErr); hGenUpsilonEffGMMinus1SystErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffGMPlus1StatErr  = new TH2F(*hGenUpsilonEffNumGMPlus1StatErr); hGenUpsilonEffGMPlus1StatErr->Divide(hGenUpsilonEffDen);
  TH2F *hGenUpsilonEffGMMinus1StatErr  = new TH2F(*hGenUpsilonEffNumGMMinus1StatErr); hGenUpsilonEffGMMinus1StatErr->Divide(hGenUpsilonEffDen);
  
  TFile *fOut = new TFile("OutEff2DSyst.root","RECREATE");


   hGenUpsilonEffDen->SetStats(0);
   hGenUpsilonEffDen->SetContour(1000);
   hGenUpsilonEffDen->SetTitle("Efficiency denominator");
   hGenUpsilonEffDen->SetName("hGenUpsilonEffDen");
   hGenUpsilonEffDen->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffDen->GetYaxis()->SetTitle("y");
   hGenUpsilonEffDen->SetOption("colz");
   
   hGenUpsilonEffDen->Write();
   
   hGenUpsilonEffNumNominal->SetStats(0);
   hGenUpsilonEffNumNominal->SetContour(1000);
   hGenUpsilonEffNumNominal->SetTitle("Efficiency numerator Nominal");
   hGenUpsilonEffNumNominal->SetName("hGenUpsilonEffNumNominal");
   hGenUpsilonEffNumNominal->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumNominal->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumNominal->SetOption("colz");
   
   hGenUpsilonEffNumNominal->Write();

   hGenUpsilonEffNominal->SetStats(0);
   hGenUpsilonEffNominal->SetContour(1000);
   hGenUpsilonEffNominal->SetTitle("Efficiency with Nominal Error Correction");
   hGenUpsilonEffNominal->SetTitleSize(1);
   hGenUpsilonEffNominal->SetName("hGenUpsilonEffNominal");
   hGenUpsilonEffNominal->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffNominal->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffNominal->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffNominal->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffNominal->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffNominal->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffNominal->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffNominal->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffNominal->SetOption("colz");
   
   hGenUpsilonEffNominal->Write();
  
   hGenUpsilonEffNumHSPlus1SystErr->SetStats(0);
   hGenUpsilonEffNumHSPlus1SystErr->SetContour(1000);
   hGenUpsilonEffNumHSPlus1SystErr->SetTitle("Efficiency numerator HSPlus1SystErr");
   hGenUpsilonEffNumHSPlus1SystErr->SetName("hGenUpsilonEffNumHSPlus1SystErr");
   hGenUpsilonEffNumHSPlus1SystErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumHSPlus1SystErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumHSPlus1SystErr->SetOption("colz");
   
   hGenUpsilonEffNumHSPlus1SystErr->Write();

   hGenUpsilonEffHSPlus1SystErr->SetStats(0);
   hGenUpsilonEffHSPlus1SystErr->SetContour(1000);
   hGenUpsilonEffHSPlus1SystErr->SetName("hGenUpsilonEffHSPlus1SystErr");
   hGenUpsilonEffHSPlus1SystErr->SetTitle("Efficiency with Plus 1 Syst. Err. Hard Soft Muon Correction");
   hGenUpsilonEffHSPlus1SystErr->SetTitleSize(1);
   hGenUpsilonEffHSPlus1SystErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSPlus1SystErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSPlus1SystErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSPlus1SystErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffHSPlus1SystErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffHSPlus1SystErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSPlus1SystErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSPlus1SystErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSPlus1SystErr->SetOption("colz");
   
   hGenUpsilonEffHSPlus1SystErr->Write();
   
   hGenUpsilonEffNumHSMinus1SystErr->SetStats(0);
   hGenUpsilonEffNumHSMinus1SystErr->SetContour(1000);
   hGenUpsilonEffNumHSMinus1SystErr->SetTitle("Efficiency numerator HSMinus1SystErr");
   hGenUpsilonEffNumHSMinus1SystErr->SetName("hGenUpsilonEffNumHSMinus1SystErr");
   hGenUpsilonEffNumHSMinus1SystErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumHSMinus1SystErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumHSMinus1SystErr->SetOption("colz");
   
   hGenUpsilonEffNumHSMinus1SystErr->Write();

   hGenUpsilonEffHSMinus1SystErr->SetStats(0);
   hGenUpsilonEffHSMinus1SystErr->SetContour(1000);
   hGenUpsilonEffHSMinus1SystErr->SetName("hGenUpsilonEffHSMinus1SystErr");
   hGenUpsilonEffHSMinus1SystErr->SetTitle("Efficiency with Minus 1 Syst. Err. Hard Soft Muon Correction");
   hGenUpsilonEffHSMinus1SystErr->SetTitleSize(1);
   hGenUpsilonEffHSMinus1SystErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSMinus1SystErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSMinus1SystErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSMinus1SystErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffHSMinus1SystErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffHSMinus1SystErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSMinus1SystErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSMinus1SystErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSMinus1SystErr->SetOption("colz");
   
   hGenUpsilonEffHSMinus1SystErr->Write();

   hGenUpsilonEffNumHSPlus1StatErr->SetStats(0);
   hGenUpsilonEffNumHSPlus1StatErr->SetContour(1000);
   hGenUpsilonEffNumHSPlus1StatErr->SetTitle("Efficiency numerator HSPlus1StatErr");
   hGenUpsilonEffNumHSPlus1StatErr->SetName("hGenUpsilonEffNumHSPlus1StatErr");
   hGenUpsilonEffNumHSPlus1StatErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumHSPlus1StatErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumHSPlus1StatErr->SetOption("colz");
   
   hGenUpsilonEffNumHSPlus1StatErr->Write();

   hGenUpsilonEffHSPlus1StatErr->SetStats(0);
   hGenUpsilonEffHSPlus1StatErr->SetContour(1000);
   hGenUpsilonEffHSPlus1StatErr->SetName("hGenUpsilonEffHSPlus1StatErr");
   hGenUpsilonEffHSPlus1StatErr->SetTitle("Efficiency with Plus 1 Stat. Err. Hard Soft Muon Correction");
   hGenUpsilonEffHSPlus1StatErr->SetTitleSize(1);
   hGenUpsilonEffHSPlus1StatErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffHSPlus1StatErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffHSPlus1StatErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSPlus1StatErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSPlus1StatErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSPlus1StatErr->SetOption("colz");
   
   hGenUpsilonEffHSPlus1StatErr->Write();
   

   hGenUpsilonEffNumHSMinus1StatErr->SetStats(0);
   hGenUpsilonEffNumHSMinus1StatErr->SetContour(1000);
   hGenUpsilonEffNumHSMinus1StatErr->SetTitle("Efficiency numerator HSMinus1StatErr");
   hGenUpsilonEffNumHSMinus1StatErr->SetName("hGenUpsilonEffNumHSMinus1StatErr");
   hGenUpsilonEffNumHSMinus1StatErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumHSMinus1StatErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumHSMinus1StatErr->SetOption("colz");
   
   hGenUpsilonEffNumHSMinus1StatErr->Write();

   hGenUpsilonEffHSMinus1StatErr->SetStats(0);
   hGenUpsilonEffHSMinus1StatErr->SetContour(1000);
   hGenUpsilonEffHSMinus1StatErr->SetName("hGenUpsilonEffHSMinus1StatErr");
   hGenUpsilonEffHSMinus1StatErr->SetTitle("Efficiency with Minus 1 Stat. Err. Hard Soft Muon Correction");
   hGenUpsilonEffHSMinus1StatErr->SetTitleSize(1);
   hGenUpsilonEffHSMinus1StatErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSMinus1StatErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSMinus1StatErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSMinus1StatErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffHSMinus1StatErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffHSMinus1StatErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffHSMinus1StatErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffHSMinus1StatErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffHSMinus1StatErr->SetOption("colz");
   
   hGenUpsilonEffHSMinus1StatErr->Write();

   
   hGenUpsilonEffNumGMPlus1SystErr->SetStats(0);
   hGenUpsilonEffNumGMPlus1SystErr->SetContour(1000);
   hGenUpsilonEffNumGMPlus1SystErr->SetTitle("Efficiency numerator GMPlus1SystErr");
   hGenUpsilonEffNumGMPlus1SystErr->SetName("hGenUpsilonEffNumGMPlus1SystErr");
   hGenUpsilonEffNumGMPlus1SystErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumGMPlus1SystErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumGMPlus1SystErr->SetOption("colz");
   
   hGenUpsilonEffNumGMPlus1SystErr->Write();

   hGenUpsilonEffGMPlus1SystErr->SetStats(0);
   hGenUpsilonEffGMPlus1SystErr->SetContour(1000);
   hGenUpsilonEffGMPlus1SystErr->SetName("hGenUpsilonEffGMPlus1SystErr");
   hGenUpsilonEffGMPlus1SystErr->SetTitle("Efficiency with Plus 1 Syst. Err. Global Muon Correction");
   hGenUpsilonEffGMPlus1SystErr->SetTitleSize(1);
   hGenUpsilonEffGMPlus1SystErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMPlus1SystErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMPlus1SystErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMPlus1SystErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffGMPlus1SystErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffGMPlus1SystErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMPlus1SystErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMPlus1SystErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMPlus1SystErr->SetOption("colz");
   
   hGenUpsilonEffGMPlus1SystErr->Write();
   
   hGenUpsilonEffNumGMMinus1SystErr->SetStats(0);
   hGenUpsilonEffNumGMMinus1SystErr->SetContour(1000);
   hGenUpsilonEffNumGMMinus1SystErr->SetTitle("Efficiency numerator GMMinus1SystErr");
   hGenUpsilonEffNumGMMinus1SystErr->SetName("hGenUpsilonEffNumGMMinus1SystErr");
   hGenUpsilonEffNumGMMinus1SystErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumGMMinus1SystErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumGMMinus1SystErr->SetOption("colz");
   
   hGenUpsilonEffNumGMMinus1SystErr->Write();

   hGenUpsilonEffGMMinus1SystErr->SetStats(0);
   hGenUpsilonEffGMMinus1SystErr->SetContour(1000);
   hGenUpsilonEffGMMinus1SystErr->SetName("hGenUpsilonEffGMMinus1SystErr");
   hGenUpsilonEffGMMinus1SystErr->SetTitle("Efficiency with Minus 1 Syst. Err. Global Muon Correction");
   hGenUpsilonEffGMMinus1SystErr->SetTitleSize(1);
   hGenUpsilonEffGMMinus1SystErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMMinus1SystErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMMinus1SystErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMMinus1SystErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffGMMinus1SystErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffGMMinus1SystErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMMinus1SystErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMMinus1SystErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMMinus1SystErr->SetOption("colz");
   
   hGenUpsilonEffGMMinus1SystErr->Write();

   hGenUpsilonEffNumGMPlus1StatErr->SetStats(0);
   hGenUpsilonEffNumGMPlus1StatErr->SetContour(1000);
   hGenUpsilonEffNumGMPlus1StatErr->SetTitle("Efficiency numerator GMPlus1StatErr");
   hGenUpsilonEffNumGMPlus1StatErr->SetName("hGenUpsilonEffNumGMPlus1StatErr");
   hGenUpsilonEffNumGMPlus1StatErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumGMPlus1StatErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumGMPlus1StatErr->SetOption("colz");
   
   hGenUpsilonEffNumGMPlus1StatErr->Write();

   hGenUpsilonEffGMPlus1StatErr->SetStats(0);
   hGenUpsilonEffGMPlus1StatErr->SetContour(1000);
   hGenUpsilonEffGMPlus1StatErr->SetName("hGenUpsilonEffGMPlus1StatErr");
   hGenUpsilonEffGMPlus1StatErr->SetTitle("Efficiency with Plus 1 Stat. Err. Global Muon Correction");
   hGenUpsilonEffGMPlus1StatErr->SetTitleSize(1);
   hGenUpsilonEffGMPlus1StatErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMPlus1StatErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMPlus1StatErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMPlus1StatErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffGMPlus1StatErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffGMPlus1StatErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMPlus1StatErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMPlus1StatErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMPlus1StatErr->SetOption("colz");
   
   hGenUpsilonEffGMPlus1StatErr->Write();

   hGenUpsilonEffNumGMMinus1StatErr->SetStats(0);
   hGenUpsilonEffNumGMMinus1StatErr->SetContour(1000);
   hGenUpsilonEffNumGMMinus1StatErr->SetTitle("Efficiency numerator GMMinus1StatErr");
   hGenUpsilonEffNumGMMinus1StatErr->SetName("hGenUpsilonEffNumGMMinus1StatErr");
   hGenUpsilonEffNumGMMinus1StatErr->GetXaxis()->SetTitle("pT (GeV)");
   hGenUpsilonEffNumGMMinus1StatErr->GetYaxis()->SetTitle("y");
   hGenUpsilonEffNumGMMinus1StatErr->SetOption("colz");
   
   hGenUpsilonEffNumGMMinus1StatErr->Write();

   hGenUpsilonEffGMMinus1StatErr->SetStats(0);
   hGenUpsilonEffGMMinus1StatErr->SetContour(1000);
   hGenUpsilonEffGMMinus1StatErr->SetName("hGenUpsilonEffGMMinus1StatErr");
   hGenUpsilonEffGMMinus1StatErr->SetTitle("Efficiency with Minus 1 Stat. Err. Global Muon Correction");
   hGenUpsilonEffGMMinus1StatErr->SetTitleSize(1);
   hGenUpsilonEffGMMinus1StatErr->GetXaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMMinus1StatErr->GetXaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMMinus1StatErr->GetXaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMMinus1StatErr->GetXaxis()->SetTitle("Upsilon pT (GeV)");
   hGenUpsilonEffGMMinus1StatErr->GetYaxis()->SetTitle("Upsilon y");
   hGenUpsilonEffGMMinus1StatErr->GetYaxis()->SetTitleSize(0.05);
   hGenUpsilonEffGMMinus1StatErr->GetYaxis()->SetTitleOffset(1);
   hGenUpsilonEffGMMinus1StatErr->GetYaxis()->SetLabelSize(0.035);
   hGenUpsilonEffGMMinus1StatErr->SetOption("colz");
   
   hGenUpsilonEffGMMinus1StatErr->Write();
   
  fOut->Close();

}
