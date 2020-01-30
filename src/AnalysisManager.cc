#include <G4SDManager.hh>
#include <G4Run.hh>
#include <G4Event.hh>
#include <G4HCofThisEvent.hh>
#include <G4EmCalculator.hh>
#include <G4Material.hh>
#include <G4HadronicProcessStore.hh>
#include <G4ParticleTable.hh>
#include <G4NistManager.hh>
#include <G4ElementTable.hh>
#include <G4Version.hh>
#include <numeric>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>
#include <TDirectory.h>
#include <TH1.h>

#include "stdHit.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventData.hh"

#include "AnalysisMessenger.hh"
#include "AnalysisManager.hh"

#include "G4UImanager.hh"
#include "DetectorConstruction.hh"
#include <algorithm>
#include "Particle.hh"
#include <cmath>

using namespace CLHEP;

//__________________________________________________________________________________________________________

AnalysisManager::AnalysisManager(PrimaryGeneratorAction *pPrimaryGeneratorAction)
{
    m_pAnalysisMessenger = new AnalysisMessenger(this);
    m_hTreeType          = "raw";
    
    runTime = new G4Timer();
    m_CollectionIDs.clear();
    
    m_LXeCollectionID = -1;
    m_hDataFilename   = "events.root";
    
    m_pPrimaryGeneratorAction = pPrimaryGeneratorAction;
    
    m_pEventData = new EventData();
    //plotPhysics      = kTRUE;
    // writeEmptyEvents = kTRUE;
}

//__________________________________________________________________________________________________________

AnalysisManager::~AnalysisManager()
{
    delete m_pAnalysisMessenger;
}

//__________________________________________________________________________________________________________

void
AnalysisManager::BeginOfRun(const G4Run *)
{
    // start a timer for this run....
    runTime->Start();
    // do we write empty events or not?
    // writeEmptyEvents = kTRUE;
    
    m_pTreeFile = new TFile(m_hDataFilename.c_str(), "RECREATE");
    // make tree structure
    TNamed *G4version = new TNamed("G4VERSION_TAG",G4VERSION_TAG);
    G4version->Write();
    
    _events = m_pTreeFile->mkdir("events");
    _events->cd();
    
    G4cout << "AnalysisManager:: Init data tree ..." << G4endl;
    m_pTree = new TTree("evt", "Event Data");
    m_pTree_vrt = new TTree("evt_vrt", "Event Data vrt");
    
    gROOT->ProcessLine("#include <vector>");
    
    m_pTree->Branch("eventid", &m_pEventData->m_iEventId, "eventid/I");
    
    m_pTree->Branch("etot",   &m_pEventData->m_fTotalEnergyDeposited, "etot/F");
    m_pTree->Branch("nsteps", &m_pEventData->m_iNbSteps, "nsteps/I");
    m_pTree->Branch("type_pri", "vector<string>", &m_pEventData->m_pPrimaryParticleType);
    m_pTree->Branch("xp_pri", &m_pEventData->m_fPrimaryX, "xp_pri/F");
    m_pTree->Branch("yp_pri", &m_pEventData->m_fPrimaryY, "yp_pri/F");
    m_pTree->Branch("zp_pri", &m_pEventData->m_fPrimaryZ, "zp_pri/F");
    m_pTree->Branch("e_pri",  &m_pEventData->m_fPrimaryE, "e_pri/F");
    m_pTree->Branch("w_pri",  &m_pEventData->m_fPrimaryW, "w_pri/F");
    
    if(m_hTreeType == "raw"){ // all hits
        m_pTree->Branch("trackid",    "vector<int>",    &m_pEventData->m_pTrackId);
        m_pTree->Branch("type",       "vector<string>", &m_pEventData->m_pParticleType);
        m_pTree->Branch("parentid",   "vector<int>",    &m_pEventData->m_pParentId);
        m_pTree->Branch("collid",     "vector<int>",    &m_pEventData->m_pCollectionId);
        m_pTree->Branch("parenttype", "vector<string>", &m_pEventData->m_pParentType);
        m_pTree->Branch("creaproc",   "vector<string>", &m_pEventData->m_pCreatorProcess);
        m_pTree->Branch("edproc",     "vector<string>", &m_pEventData->m_pDepositingProcess);
        m_pTree->Branch("xp",         "vector<float>",  &m_pEventData->m_pX);
        m_pTree->Branch("yp",         "vector<float>",  &m_pEventData->m_pY);
        m_pTree->Branch("zp",         "vector<float>",  &m_pEventData->m_pZ);
        m_pTree->Branch("ed",         "vector<float>",  &m_pEventData->m_pEnergyDeposited);
        m_pTree->Branch("time",       "vector<float>",  &m_pEventData->m_pTime);
    } else if (m_hTreeType == "compact"){ // only average per detector
        m_pTree->Branch("collid",     "vector<int>",    &m_pEventData->m_pCollectionId);
        m_pTree->Branch("xp",         "vector<float>",  &m_pEventData->m_pX);
        m_pTree->Branch("yp",         "vector<float>",  &m_pEventData->m_pY);
        m_pTree->Branch("zp",         "vector<float>",  &m_pEventData->m_pZ);
        m_pTree->Branch("ed",         "vector<float>",  &m_pEventData->m_pEnergyDeposited);
        m_pTree->Branch("time",       "vector<float>",  &m_pEventData->m_pTime);
    }
    
    //Data part of vrt MC
    m_pTree_vrt->Branch("eventid", &m_pEventData->m_iEventId_vrt,  "eventid/I");
    m_pTree_vrt->Branch("xp_pri",  &m_pEventData->m_fPrimaryX_vrt, "xp_pri/F");
    m_pTree_vrt->Branch("yp_pri",  &m_pEventData->m_fPrimaryY_vrt, "yp_pri/F");
    m_pTree_vrt->Branch("zp_pri",  &m_pEventData->m_fPrimaryZ_vrt, "zp_pri/F");
    m_pTree_vrt->Branch("e_pri",   &m_pEventData->m_fPrimaryE_vrt, "e_pri/F");
    m_pTree_vrt->Branch("w_pri",   &m_pEventData->m_fPrimaryW_vrt, "w_pri/F");
    m_pTree_vrt->Branch("xp",     "vector<float>",  &m_pEventData->m_pX_vrt);
    m_pTree_vrt->Branch("yp",     "vector<float>",  &m_pEventData->m_pY_vrt);
    m_pTree_vrt->Branch("zp",     "vector<float>",  &m_pEventData->m_pZ_vrt);
    m_pTree_vrt->Branch("ed",     "vector<float>",  &m_pEventData->m_pEnergyDeposited_vrt);
    m_pTree_vrt->Branch("edproc", "vector<string>", &m_pEventData->m_pDepositingProcess_vrt);


    m_pNbEventsToSimulateParameter = new TParameter<int>("nbevents", m_iNbEventsToSimulate);
    m_pNbEventsToSimulateParameter->Write();

    m_pTreeFile->cd();
}

//__________________________________________________________________________________________________________

void
AnalysisManager::EndOfRun(const G4Run *)
{
    runTime->Stop();
    G4double dt = runTime->GetRealElapsed();
    // make tree structure
    TParameter<G4double> *dtPar = new TParameter<G4double>("G4RUNTIME",dt);
    dtPar->Write();
    
    m_pTreeFile->cd();
    
    m_pTreeFile->Write();
    m_pTreeFile->Close();
}

//__________________________________________________________________________________________________________

void
AnalysisManager::BeginOfEvent(const G4Event *pEvent)
{   
    // only do this if the collection has not yet been defined yet
    if(m_LXeCollectionID == -1)
    {   G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
        m_LXeCollectionID = pSDManager->GetCollectionID("LXe/HitsCollection");
        m_CollectionIDs.push_back(m_LXeCollectionID);
        G4cout << "AnalysisManager::BeginOfEvent  Found LXe collection at ID = " << m_LXeCollectionID << G4endl;
    }
    
    /*
     * Following code checks whether the initial path of the primary would intersect with the FV if it
     * would not scatter. Only particles that would hit the FV are kept
     */
    /*
    G4PrimaryVertex   * primaryVertex   = pEvent->GetPrimaryVertex();
    G4PrimaryParticle * primaryParticle = primaryVertex->GetPrimary();

    // Getting position and momentum of the initial particle
    G4ThreeVector pos = primaryVertex->GetPosition();
    G4ThreeVector mom = primaryParticle->GetMomentumDirection();
    G4double      ene = primaryParticle->GetKineticEnergy();
  
    Particle myParticle;

    myParticle.setX0(pos);
    myParticle.setDirection(mom);
    myParticle.setEnergy(ene);
    myParticle.setX0start(pos);
    myParticle.setVrt("fiducial_scatter");
    
    myParticle.Propagate();
    
    m_pEventData->Clear();
    
    _events->cd();

    m_pEventData->m_iEventId_vrt  = pEvent->GetEventID();
    m_pEventData->m_fPrimaryX_vrt = myParticle.getX0start()[0] / mm;
    m_pEventData->m_fPrimaryY_vrt = myParticle.getX0start()[1] / mm;
    m_pEventData->m_fPrimaryZ_vrt = myParticle.getX0start()[2] / mm;
    m_pEventData->m_fPrimaryE_vrt = (myParticle.getEnergy() + myParticle.getEdep()) / keV;
    m_pEventData->m_fPrimaryW_vrt = myParticle.getWeight();
    
    for (G4int i; i < myParticle.getEdep_int().size(); i++)
    {   
        m_pEventData->m_pDepositingProcess_vrt->push_back(myParticle.getPro_int()[i]);
    
        m_pEventData->m_pX_vrt->push_back(myParticle.getPos_int()[i][0] / mm);
        m_pEventData->m_pY_vrt->push_back(myParticle.getPos_int()[i][1] / mm);
        m_pEventData->m_pZ_vrt->push_back(myParticle.getPos_int()[i][2] / mm);
                        
        m_pEventData->m_pEnergyDeposited_vrt->push_back(myParticle.getEdep_int()[i] / keV);
    }

    m_pTree_vrt->Fill(); // write all events to the tree
    
    m_pEventData->Clear();
    m_pTreeFile->cd();
    */
}

//__________________________________________________________________________________________________________

void
AnalysisManager::EndOfEvent(const G4Event *pEvent)
{
    m_pEventData->Clear();
    
    _events->cd();
    
    G4HCofThisEvent* pHCofThisEvent = pEvent->GetHCofThisEvent();
    
    // get the event ID and primary particle information
    m_pEventData->m_iEventId = pEvent->GetEventID();
    m_pEventData->m_pPrimaryParticleType->push_back(m_pPrimaryGeneratorAction->GetParticleTypeOfPrimary());
    
    m_pEventData->m_fPrimaryX = m_pPrimaryGeneratorAction->GetPositionOfPrimary().x();
    m_pEventData->m_fPrimaryY = m_pPrimaryGeneratorAction->GetPositionOfPrimary().y();
    m_pEventData->m_fPrimaryZ = m_pPrimaryGeneratorAction->GetPositionOfPrimary().z();
    m_pEventData->m_fPrimaryE = m_pPrimaryGeneratorAction->GetEnergyOfPrimary() / keV;
    m_pEventData->m_fPrimaryW = pEvent->GetPrimaryVertex()->GetWeight();
    
    //G4cout  << "Primary particle ID: " << m_pEventData->m_iEventId << G4endl;
    
    // unpack the hit collections
    G4int    iNbHits               = 0;
    G4int    iNbSteps              = 0;
    G4double fTotalEnergyDeposited = 0;
    
    if(pHCofThisEvent) {
        //
        // loop over all our hit collections
        //
        
        for(G4int icol = 0; icol < (G4int)m_CollectionIDs.size(); icol++){
            // check if the ID of the collection is OK
            if(m_CollectionIDs[icol] != -1){
                stdHitsCollection* pHitsCollection = 0;
                // get the hits
                pHitsCollection = static_cast<stdHitsCollection*>(pHCofThisEvent->GetHC(m_CollectionIDs[icol]));
                               
                pHitsCollection->DrawAllHits();
                
                // the number of hits
                iNbHits = (pHitsCollection)?(pHitsCollection->entries()):(0);
                
                //G4cout << icol << " Nb hits = " << iNbHits<< " ID = " << m_CollectionIDs[icol] << " (AnalysisManager)" << G4endl;
                
                if(iNbHits) {
                    // hits
                    if (m_hTreeType == "raw"){
                        // write all the GEANT4 steps and hits to the output tree
                        for(G4int i = 0; i < iNbHits; i++) {
                            stdHit *pHit = (*pHitsCollection)[i];
                            if(pHit->GetParticleType() != "opticalphoton"){
                                m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
                                m_pEventData->m_pParentId->push_back(pHit->GetParentId());
                                m_pEventData->m_pCollectionId->push_back(icol);
                                
                                m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
                                m_pEventData->m_pParentType->push_back(pHit->GetParentType());
                                m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
                                m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());
                                
                                m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
                                m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
                                m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);
                                
                                fTotalEnergyDeposited += pHit->GetEnergyDeposited()/keV;
                                m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);
                                
                                m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
                                m_pEventData->m_pTime->push_back(pHit->GetTime()/second);
                                
                                iNbSteps++;
                            } // if !opticalphoton
                        }
                    } else if (m_hTreeType == "compact"){
                        //
                        // calculate the average position of a hit in each detector
                        //
                                                
                        for(G4int i = 0; i < iNbHits; i++) {

                            G4double xx = 0;
                            G4double yy = 0;
                            G4double zz = 0;
                            G4double ee = 0;
                            G4double tt = 0;

                            stdHit *pHit = (*pHitsCollection)[i];
                            if (pHit->GetTrackBanner()){
                                continue;
                            }
                            
                            G4double x1 = pHit->GetPosition().x();
                            G4double y1 = pHit->GetPosition().y();
                            G4double z1 = pHit->GetPosition().z();
                            G4double t1 = pHit->GetTime();

                            for(G4int j = 0; j < iNbHits; j++) { 

                                stdHit *pHit2 = (*pHitsCollection)[j];
                                if (pHit2->GetTrackBanner()){
                                    continue;
                                }
                                
                                G4double ed2 = pHit2->GetEnergyDeposited()/keV;

                                G4double x2 = pHit2->GetPosition().x();
                                G4double y2 = pHit2->GetPosition().y();
                                G4double z2 = pHit2->GetPosition().z();
                                G4double t2 = pHit2->GetTime();

                                G4double distance = sqrt( pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2) );
                                G4double time_dif = abs ( t2 - t1 );

                                if (distance <= 0 * mm && time_dif <= 0 * ns){//doesnt cluster now, cause doing so in root file
                                    pHit2->SetTrackBanner(true);
                                    xx += x2/mm * ed2;
                                    yy += y2/mm * ed2;
                                    zz += z2/mm * ed2;
                                    tt += t2/second     * ed2;
                                    ee += ed2;
                                }

                            }

                            if(ee!=0){
                                xx /= ee;
                                yy /= ee;
                                zz /= ee;
                                tt /= ee;
                            }
                        
                            // fill the tree variables
                            m_pEventData->m_pCollectionId->push_back(icol);
                            m_pEventData->m_pEnergyDeposited->push_back(ee);
                            m_pEventData->m_pX->push_back(xx);
                            m_pEventData->m_pY->push_back(yy);
                            m_pEventData->m_pZ->push_back(zz);
                            m_pEventData->m_pTime->push_back(tt);

                            fTotalEnergyDeposited += ee;
                        }                       
                                                
                    } else {
                        G4cout << "AnalysisManager::EndOfEvent ERROR: wrong Tree type selected" << G4endl;
                        return;
                    } // if m_hTreeType
                } // if Hits
            }
        }
    }
    
    //G4cout << " size = " << m_pEventData->m_pX->size() <<G4endl;
    // also write the header information + primary vertex of the empty events....
    m_pEventData->m_iNbSteps              = iNbSteps;
    m_pEventData->m_fTotalEnergyDeposited = fTotalEnergyDeposited;
    
    // save only energy depositing events
    m_pTree->Fill(); // write all events to the tree
    
    m_pEventData->Clear();
    m_pTreeFile->cd();
    
    //  delete pHitsCollection;
}

//__________________________________________________________________________________________________________

void
AnalysisManager::Step(const G4Step *)
{   
}

//__________________________________________________________________________________________________________