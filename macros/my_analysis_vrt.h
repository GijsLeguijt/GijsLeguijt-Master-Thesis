//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan 27 13:34:34 2020 by ROOT version 6.14/04
// from TTree evt/Event Data
// found on file: events.root
//////////////////////////////////////////////////////////

#ifndef my_analysis_vrt_h
#define my_analysis_vrt_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class my_analysis_vrt {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           eventid;
   Float_t         xp_pri;
   Float_t         yp_pri;
   Float_t         zp_pri;
   Float_t         e_pri;
   Float_t         w_pri;
   vector<string>  *edproc;
   vector<float>   *xp;
   vector<float>   *yp;
   vector<float>   *zp;
   vector<float>   *ed;

   // List of branches
   TBranch        *b_eventid;   //!
   TBranch        *b_xp_pri;   //!
   TBranch        *b_yp_pri;   //!
   TBranch        *b_zp_pri;   //!
   TBranch        *b_e_pri;   //!
   TBranch        *b_w_pri;   //!
   TBranch        *b_edproc;   //!
   TBranch        *b_xp;   //!
   TBranch        *b_yp;   //!
   TBranch        *b_zp;   //!
   TBranch        *b_ed;   //!

   my_analysis_vrt(TTree *tree=0);
   virtual ~my_analysis_vrt();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef my_analysis_vrt_cxx
my_analysis_vrt::my_analysis_vrt(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("events.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("events.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("events.root:/events");
      dir->GetObject("evt",tree);

   }
   Init(tree);
}

my_analysis_vrt::~my_analysis_vrt()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t my_analysis_vrt::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t my_analysis_vrt::LoadTree(Long64_t entry)
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

void my_analysis_vrt::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   edproc = 0;
   xp = 0;
   yp = 0;
   zp = 0;
   ed = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("xp_pri", &xp_pri, &b_xp_pri);
   fChain->SetBranchAddress("yp_pri", &yp_pri, &b_yp_pri);
   fChain->SetBranchAddress("zp_pri", &zp_pri, &b_zp_pri);
   fChain->SetBranchAddress("e_pri", &e_pri, &b_e_pri);
   fChain->SetBranchAddress("w_pri", &w_pri, &b_w_pri);
   fChain->SetBranchAddress("edproc", &edproc, &b_edproc);
   fChain->SetBranchAddress("xp", &xp, &b_xp);
   fChain->SetBranchAddress("yp", &yp, &b_yp);
   fChain->SetBranchAddress("zp", &zp, &b_zp);
   fChain->SetBranchAddress("ed", &ed, &b_ed);
   Notify();
}

Bool_t my_analysis_vrt::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void my_analysis_vrt::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t my_analysis_vrt::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef my_analysis_vrt_cxx