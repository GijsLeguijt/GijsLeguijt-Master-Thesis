#define my_analysis_cxx
#include "my_analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void my_analysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L my_analysis.C
//      root> my_analysis t
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
      vector<float> t = * time;
      vector<string> par_type = *parenttype;
      vector<string> dep_proc = *edproc;

      struct Hits //data-structure of individual hits
      {
         float x, y, z, e, t;
         string parent, process;
      };

      std::vector<Hits> hits, clusters; //hits are individual hits, clusters are clustered hits

      for (int i = 0; i < x.size(); i++){
         hits.push_back({x[i],y[i],z[i],e[i],t[i],par_type[i],dep_proc[i]});
      }
   
      //sorting the hits in ascending time order
      sort(hits.begin(), hits.end(), [](const Hits& h1, const Hits& h2) { return h1.t < h2.t; });


      //
      // Clustering of hits
      //
      for (Long64_t i = 0; i < hits.size(); i++) {

         if (hits[i].parent != "none"){//not primary particle
            continue;
         }         

         if (hits[i].e == 0){//G4 also adds hits without energy deposits, e.g. transport
            continue;
         }

         float xx   = 0;
         float yy   = 0;
         float zz   = 0;
         float ee   = 0;
         float tt   = 0;
         int n_prim = 0;

         for (Long64_t j=i; j<hits.size();j++) {

            if (hits[j].parent == "none"){//if primary particle
               n_prim += 1;
            }
            
            if (n_prim > 1){//everything happening after next primary reaction belongs to that reaction
               break;
            }

            //
            // Averaging using energy as weight
            //
            double x2 = hits[j].x;
            double y2 = hits[j].y;
            double z2 = hits[j].z;
            double e2 = hits[j].e;
            double t2 = hits[j].t;

            xx += x2 * e2;
            yy += y2 * e2;
            zz += z2 * e2;
            tt += t2 * e2;
            ee += e2;
         }

         xx /= ee;
         yy /= ee;
         zz /= ee;
         tt /= ee;

         clusters.push_back({xx,yy,zz,ee,tt,par_type[i],dep_proc[i]});

      }

      //
      // Data cuts, cutting on clustered data
      //
      int nhits = clusters.size();
      
      if (nhits > 0)      
      {  if (abs(clusters[0].z) > 670 ) {continue;}

         double r = sqrt( pow(clusters[0].x,2) + pow(clusters[0].y,2) );
         
         if (r > 570) {continue;}

         if (nhits == 1) {nhits1++; myHisto1->Fill(clusters[0].e);}
         if (nhits == 2) {nhits2++; myHisto2->Fill(clusters[0].e);}
         if (nhits >= 3) {nhits3++; myHisto3->Fill(clusters[0].e);}

         myHisto->Fill(clusters[0].e);
      }

   }

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