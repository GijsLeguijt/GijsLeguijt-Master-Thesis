#define my_analysis_vrt_cxx
#include "my_analysis_vrt.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void my_analysis_vrt::Loop()
{
//   In a ROOT session, you can do:
//      root> .L my_analysis_vrt.C
//      root> my_analysis_vrt t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TCanvas *myCanvas = new TCanvas();
   myCanvas->Divide(2,2);
   TH1F *myHisto   = new TH1F("myHisto", "Energy spectrum all;energy (keV);counts",  100,0,1050);
   TH1F *myHisto1  = new TH1F("myHisto1","Energy spectrum 1;energy (keV);counts",    100,0,1050);
   TH1F *myHisto2  = new TH1F("myHisto2","Energy spectrum 2;energy (keV);counts",    100,0,1050);
   TH1F *myHisto3  = new TH1F("myHisto3","Energy spectrum more;energy (keV);counts", 100,0,1050);
   gPad->SetLogy();

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   int nhits1 = 0, nhits2 = 0, nhits3 = 0; //number of hits in an event, nhits3 is for all >2
   int maxhits = 0;

   for (Long64_t jentry=0; jentry < nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      nb = fChain->GetEntry(jentry);   nbytes += nb; //gets to next event

      //
      // Getting data from the simulation
      //
      vector<float> x = * xp;
      vector<float> y = * yp;
      vector<float> z = * zp;
      vector<float> e = * ed;
      vector<string> dep_proc = * edproc;

      struct Hits //data-structure of individual hits
      {
         float x, y, z, e;
         string process;
      };

      std::vector<Hits> hits, clusters; //hits are individual hits, clusters are clustered hits

      for (int i = 0; i < x.size(); i++){
         hits.push_back({x[i],y[i],z[i],e[i],dep_proc[i]});
      }

      //
      // Data cuts, cutting on clustered data
      //
      int nhits = clusters.size();
      
      if (nhits > 0)
      {  if (nhits > maxhits) {maxhits = nhits;}
         
         if (abs(clusters[0].z) > 670 ) {continue;}

         double r = sqrt( pow(clusters[0].x,2) + pow(clusters[0].y,2) );
         
         if (r > 570) {continue;}

         if (nhits == 1) {nhits1++; myHisto1->Fill(clusters[0].e);}
         if (nhits == 2) {nhits2++; myHisto2->Fill(clusters[0].e);}
         if (nhits >= 3) {nhits3++; myHisto3->Fill(clusters[0].e);}

         myHisto->Fill(clusters[0].e);
      }

   }
   cout << maxhits << " and: " << nentries << endl;

   cout << "Times 1 hit: "     << nhits1 << endl;
   cout << "Times 2 hits: "    << nhits2 << endl;
   cout << "Times more hits: " << nhits3 << endl;
   
   myCanvas->cd(1);
   gPad->SetLogy();
   myHisto->Draw("c1");
   myCanvas->cd(2);
   gPad->SetLogy();
   myHisto1->Draw("c2");
   myCanvas->cd(3);
   gPad->SetLogy();
   myHisto2->Draw("c3");
   myCanvas->cd(4);
   gPad->SetLogy();
   myHisto3->Draw("c4");
}