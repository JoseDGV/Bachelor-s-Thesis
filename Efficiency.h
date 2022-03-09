/////////////////////////////////////////////////////////////////
//The purpose of this header is to read the Efficiency.root file
//////////////////////////////////////////////////////////////////

#ifndef Efficiency_h
#define Efficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class Efficiency {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   TClonesArray    *Gen_QQ_4mom;
   TClonesArray    *Gen_mu_4mom;
   TClonesArray    *Reco_QQ_4mom;
   TClonesArray    *Reco_mu_4mom;
   Int_t            Gen_QQ_size;
   Int_t            Gen_QQ_mupl_idx[1000];
   Int_t            Gen_QQ_mumi_idx[1000];
   Int_t            Gen_QQ_whichRec[1000];
   Int_t            Gen_mu_whichRec[1000];
   Int_t            Reco_QQ_size;
   Int_t            Reco_QQ_sign[1000];
   ULong64_t        Reco_QQ_trig[1000];
   Int_t            Reco_QQ_mupl_idx[1000];
   Int_t            Reco_QQ_mumi_idx[1000];
   Float_t          Reco_QQ_VtxProb[1000];
   Int_t            Reco_mu_size;
   Int_t            Reco_mu_whichGen[1000];
   Int_t            Reco_mu_nTrkWMea[1000];
   Int_t            Reco_mu_nPixWMea[1000];
   Float_t          Reco_mu_dxy[1000];
   Float_t          Reco_mu_dxyErr[1000];  
   Float_t          Reco_mu_dz[1000];  
   Float_t          Reco_mu_dzErr[1000];
   ULong64_t        Reco_mu_trig[1000];
   Int_t            Reco_mu_SelectionType[1000];
   ULong64_t        HLTriggers;
   Int_t            nTrig;
   Int_t            trigPrescale[15];          

   // List of branches
   TBranch        *b_Gen_QQ_size;       //!
   TBranch        *b_Gen_QQ_4mom;       //!
   TBranch        *b_Gen_QQ_whichRec;   //!
   TBranch        *b_Gen_QQ_mupl_idx;   //!
   TBranch        *b_Gen_QQ_mumi_idx;   //!
   TBranch        *b_Gen_mu_4mom;       //!
   TBranch        *b_Gen_mu_whichRec;   //!
   TBranch        *b_Reco_QQ_4mom;      //!
   TBranch        *b_Reco_QQ_size;      //!
   TBranch        *b_Reco_QQ_sign;      //!
   TBranch        *b_Reco_QQ_mupl_idx;  //!
   TBranch        *b_Reco_QQ_mumi_idx;  //!
   TBranch        *b_Reco_QQ_trig;      //!
   TBranch        *b_Reco_QQ_VtxProb;   //!
   TBranch        *b_Reco_mu_size;      //!
   TBranch        *b_Reco_mu_whichGen;  //!
   TBranch        *b_Reco_mu_4mom;      //!
   TBranch        *b_Reco_mu_nTrkWMea;  //!
   TBranch        *b_Reco_mu_nPixWMea;  //!
   TBranch        *b_Reco_mu_dxy;       //!
   TBranch        *b_Reco_mu_dxyErr;    //!
   TBranch        *b_Reco_mu_dz;        //!
   TBranch        *b_Reco_mu_dzErr;     //!
   TBranch        *b_Reco_mu_trig;      //!
   TBranch        *b_Reco_mu_SelectionType; 
   TBranch        *b_HLTriggers;        //!
   TBranch        *b_nTrig;             //!
   TBranch        *b_trigPrescale;      //!
   
   
   
   Efficiency(TTree *tree=0);
   virtual ~Efficiency();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Efficiency_cxx
Efficiency::Efficiency(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Efficiency.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Efficiency.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Efficiency.root:/hionia");
      dir->GetObject("myTree",tree);

   }
   Init(tree);
}

Efficiency::~Efficiency()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Efficiency::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Efficiency::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Efficiency::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Gen_QQ_4mom = 0;
   Gen_mu_4mom = 0;
   Reco_QQ_4mom = 0;
   Reco_mu_4mom = 0;

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Gen_QQ_size", &Gen_QQ_size, &b_Gen_QQ_size);
   fChain->SetBranchAddress("Gen_QQ_4mom", &Gen_QQ_4mom, &b_Gen_QQ_4mom);
   fChain->SetBranchAddress("Gen_mu_4mom", &Gen_mu_4mom, &b_Gen_mu_4mom);
   fChain->SetBranchAddress("Gen_mu_whichRec", &Gen_mu_whichRec, &b_Gen_mu_whichRec);
   fChain->SetBranchAddress("Gen_QQ_whichRec", &Gen_QQ_whichRec, &b_Gen_QQ_whichRec);
   fChain->SetBranchAddress("Gen_QQ_mupl_idx", &Gen_QQ_mupl_idx, &b_Gen_QQ_mupl_idx);
   fChain->SetBranchAddress("Gen_QQ_mumi_idx", &Gen_QQ_mumi_idx, &b_Gen_QQ_mumi_idx);
   fChain->SetBranchAddress("Reco_mu_whichGen", &Reco_mu_whichGen, &b_Reco_mu_whichGen);
   fChain->SetBranchAddress("Reco_QQ_size", &Reco_QQ_size, &b_Reco_QQ_size);
   fChain->SetBranchAddress("Reco_QQ_sign", &Reco_QQ_sign, &b_Reco_QQ_sign);
   fChain->SetBranchAddress("Reco_mu_size", &Reco_mu_size, &b_Reco_mu_size);
   fChain->SetBranchAddress("Reco_QQ_4mom", &Reco_QQ_4mom, &b_Reco_QQ_4mom);
   fChain->SetBranchAddress("Reco_mu_4mom", &Reco_mu_4mom, &b_Reco_mu_4mom);
   fChain->SetBranchAddress("Reco_QQ_VtxProb", &Reco_QQ_VtxProb, &b_Reco_QQ_VtxProb);
   fChain->SetBranchAddress("Reco_mu_nTrkWMea", &Reco_mu_nTrkWMea, &b_Reco_mu_nTrkWMea);
   fChain->SetBranchAddress("Reco_mu_nPixWMea", &Reco_mu_nPixWMea, &b_Reco_mu_nPixWMea);
   fChain->SetBranchAddress("Reco_mu_dxy", &Reco_mu_dxy, &b_Reco_mu_dxy);
   fChain->SetBranchAddress("Reco_mu_dxyErr", &Reco_mu_dxyErr, &b_Reco_mu_dxyErr);
   fChain->SetBranchAddress("Reco_mu_dz", &Reco_mu_dz, &b_Reco_mu_dz);
   fChain->SetBranchAddress("Reco_mu_dzErr", &Reco_mu_dzErr, &b_Reco_mu_dzErr);
   fChain->SetBranchAddress("Reco_mu_trig", &Reco_mu_trig, &b_Reco_mu_trig);
   fChain->SetBranchAddress("Reco_mu_SelectionType", &Reco_mu_SelectionType, &b_Reco_mu_SelectionType);
   fChain->SetBranchAddress("HLTriggers", &HLTriggers, &b_HLTriggers);
   fChain->SetBranchAddress("nTrig", &nTrig, &b_nTrig);
   fChain->SetBranchAddress("trigPrescale", &trigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("Reco_QQ_mupl_idx", &Reco_QQ_mupl_idx, &b_Reco_QQ_mupl_idx);
   fChain->SetBranchAddress("Reco_QQ_mumi_idx", &Reco_QQ_mumi_idx, &b_Reco_QQ_mumi_idx);
   fChain->SetBranchAddress("Reco_QQ_trig", &Reco_QQ_trig, &b_Reco_QQ_trig);
   Notify();
}

Bool_t Efficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Efficiency::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Efficiency::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Efficiency_cxx
