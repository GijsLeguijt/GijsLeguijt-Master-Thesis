#include <G4HCofThisEvent.hh>
#include <G4Step.hh>
#include <G4VProcess.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>

#include <map>

using std::map;

#include "stdHit.hh"
#include "SensitiveDetector.hh"

//__________________________________________________________________________________________________________

SensitiveDetector::SensitiveDetector(G4String hName): G4VSensitiveDetector(hName), m_pHitsCollection(0), iHitsCollectionID(-1)
{
	collectionName.insert("HitsCollection");
}

//__________________________________________________________________________________________________________

SensitiveDetector::~SensitiveDetector()
{
}

//__________________________________________________________________________________________________________

void SensitiveDetector::Initialize(G4HCofThisEvent* pHitsCollectionOfThisEvent)
{
	m_pHitsCollection = new stdHitsCollection(SensitiveDetectorName, collectionName[0]);
	//m_pHitsCollection->PrintAllHits();
	if(iHitsCollectionID < 0)
		iHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(m_pHitsCollection);

    	//G4cout << "SensitiveDetector::Initialize ID = " << iHitsCollectionID << " (SensitiveDetector)" << G4endl;
	
	pHitsCollectionOfThisEvent->AddHitsCollection(iHitsCollectionID, m_pHitsCollection);

	//m_hParticleTypes.clear();
}

//__________________________________________________________________________________________________________

G4bool SensitiveDetector::ProcessHits(G4Step* pStep, G4TouchableHistory *)
{
	G4double dEnergyDeposited = pStep->GetTotalEnergyDeposit();
	G4Track *pTrack 		  = pStep->GetTrack();
	
	stdHit* pHit = new stdHit();

	pHit->SetTrackId(pTrack->GetTrackID());
	//	G4cout <<" Track weight = "<<pTrack->GetWeight()<<G4endl;

	if(!m_hParticleTypes.count(pTrack->GetTrackID()))
		m_hParticleTypes[pTrack->GetTrackID()] = pTrack->GetDefinition()->GetParticleName();

	pHit->SetParentId(pTrack->GetParentID());
	pHit->SetParticleType(pTrack->GetDefinition()->GetParticleName());

	if(pTrack->GetParentID())
		pHit->SetParentType(m_hParticleTypes[pTrack->GetParentID()]);
	else
		pHit->SetParentType(G4String("none"));

	if(pTrack->GetCreatorProcess())
		pHit->SetCreatorProcess(pTrack->GetCreatorProcess()->GetProcessName());
	else
		pHit->SetCreatorProcess(G4String("Null"));

	pHit->SetDepositingProcess(pStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName());
	pHit->SetPosition(pStep->GetPostStepPoint()->GetPosition());
	pHit->SetEnergyDeposited(dEnergyDeposited);
	pHit->SetKineticEnergy(pTrack->GetKineticEnergy());
	pHit->SetTime(pTrack->GetGlobalTime());

	m_pHitsCollection->insert(pHit);

	return true;
}

//__________________________________________________________________________________________________________

void SensitiveDetector::EndOfEvent(G4HCofThisEvent *)
{
//  if (verboseLevel>0) { 
//     G4int NbHits = trackerCollection->entries();
//     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
//            << " hits in the tracker chambers: " << G4endl;
//     for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
//    } 
}

//__________________________________________________________________________________________________________