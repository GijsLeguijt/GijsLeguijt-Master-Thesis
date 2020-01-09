#include "SteppingAction.hh"
#include "AnalysisManager.hh"

#include "G4SteppingManager.hh"

#include <string.h>
#include <cmath>

//__________________________________________________________________________________________________________

SteppingAction::SteppingAction(AnalysisManager *myAM):myAnalysisManager(myAM)
{
}

//__________________________________________________________________________________________________________

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{   //
    // Part below here identifies new materials
    //

    /*G4StepPoint* preStepPoint = aStep->GetPreStepPoint();

    if(preStepPoint->GetStepStatus() == fGeomBoundary)
    {
        G4Material* newMaterial = aStep->GetPreStepPoint()->GetMaterial();
        G4cout << "New material:" << newMaterial->GetName() << G4endl;
    }*/
}

//__________________________________________________________________________________________________________