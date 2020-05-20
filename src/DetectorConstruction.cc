#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Cons.hh>
#include <G4Sphere.hh>
#include <G4Torus.hh>
#include <G4Orb.hh>
#include <G4Polyhedra.hh>
#include <G4Polycone.hh>
#include <G4Ellipsoid.hh>
#include <G4ExtrudedSolid.hh>  // SERENA
//#include <G4Trd.hh>
#include <G4UnionSolid.hh>
//#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4SDManager.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <globals.hh>

#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using std::vector;
using std::stringstream;
using std::max;

#include "DetectorMessenger.hh"

#include "SensitiveDetector.hh"
#include "DetectorConstruction.hh"

#include "G4PhysicalVolumeStore.hh"

#include "TFile.h"
#include "TParameter.h"

using namespace CLHEP;

map<G4String, G4double> DetectorConstruction::m_hGeometryParameters;

//__________________________________________________________________________________________________________

DetectorConstruction::DetectorConstruction(G4String fName)
{
    m_pDetectorMessenger = new DetectorMessenger(this);
    detRootFile          = fName;
}

//__________________________________________________________________________________________________________

DetectorConstruction::~DetectorConstruction()
{
    delete m_pDetectorMessenger;
}

//__________________________________________________________________________________________________________

G4VPhysicalVolume*
DetectorConstruction::Construct()
{
    DefineMaterials();
    
    DefineGeometryParameters();
    
    // construction of the geometry
    ConstructLaboratory();
    // TPC
    ConstructCollimatorSystem();
    // make the appropriate components active
    ConstructSD();
    
    PrintGeometryInformation();
    
    MakeDetectorPlots();
    
    return m_pLabPhysicalVolume;
}

//__________________________________________________________________________________________________________

void
DetectorConstruction::DefineMaterials()
{
    G4NistManager* pNistManager = G4NistManager::Instance();
    
    //========================================= elements ==========================================
    
    G4Element *H  = new G4Element("Hydrogen",  "H",   1.,  1.0079 *g/mole);
    G4Element *C  = new G4Element("Carbon",    "C",   6., 12.011  *g/mole);
    G4Element *N  = new G4Element("Nitrogen",  "N",   7., 14.007  *g/mole);
    G4Element *O  = new G4Element("Oxygen",    "O",   8., 15.999  *g/mole);
    G4Element *F  = new G4Element("Fluorine",  "F",   9., 18.998  *g/mole);
    G4Element *Si = new G4Element("Silicon",   "Si", 14., 28.086  *g/mole);
    G4Element *Cr = new G4Element("Chromium",  "Cr", 24., 51.996  *g/mole);
    G4Element *Mn = new G4Element("Manganese", "Mn", 25., 54.938  *g/mole);
    G4Element *Fe = new G4Element("Iron",      "Fe", 26., 55.85   *g/mole);
    G4Element *Ni = new G4Element("Nickel",    "Ni", 28., 58.693  *g/mole);
    
    //========================================= materials =========================================
    
    //------------------------------------------ vacuum -------------------------------------------
    G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3, 2, kStateGas);
    Vacuum->AddElement(N, 0.755);
    Vacuum->AddElement(O, 0.245);
    
    //-------------------------------------------- air --------------------------------------------
    G4Material *Air = pNistManager->FindOrBuildMaterial("G4_AIR");

    //------------------------------------------- xenon -------------------------------------------
    //G4Isotope *iso_Xe = new G4Isotope("Xe1313", 54, 128, 128*g/mole);//131.3*g/mole);
    //G4Element  *Xe = new G4Element("Xe", "Xe", 1);
    //Xe->AddIsotope(iso_Xe, 100.*perCent);
    
    G4Element *Xe = pNistManager->FindOrBuildElement("Xe");
    G4Material *LXe = new G4Material("LXe", 3.0*g/cm3, 1, kStateLiquid);
    
    LXe->AddElement(Xe, 1.000);

    //-------------------------------------- Xenon-Isotopes ---------------------------------------
    G4Isotope *iso_Xe124 = new G4Isotope("Xe124", 54, 124, 124*g/mole);
    G4Isotope *iso_Xe126 = new G4Isotope("Xe126", 54, 126, 126*g/mole);
    G4Isotope *iso_Xe128 = new G4Isotope("Xe128", 54, 128, 128*g/mole);
    G4Isotope *iso_Xe129 = new G4Isotope("Xe129", 54, 129, 129*g/mole);
    G4Isotope *iso_Xe130 = new G4Isotope("Xe130", 54, 130, 130*g/mole);
    G4Isotope *iso_Xe131 = new G4Isotope("Xe131", 54, 131, 131*g/mole);
    G4Isotope *iso_Xe132 = new G4Isotope("Xe132", 54, 132, 132*g/mole);
    G4Isotope *iso_Xe134 = new G4Isotope("Xe134", 54, 134, 134*g/mole);
    G4Isotope *iso_Xe136 = new G4Isotope("Xe136", 54, 136, 136*g/mole);

    // Adding the isotopes as elements such that cross-sections are plotted
    G4Element  *ele_Xe124 = new G4Element("element Xe124", "Xe", 1);
    ele_Xe124->AddIsotope(iso_Xe124, 100.*perCent);
    G4Element  *ele_Xe126 = new G4Element("element Xe126", "Xe", 1);
    ele_Xe126->AddIsotope(iso_Xe126, 100.*perCent);
    G4Element  *ele_Xe128 = new G4Element("element Xe128", "Xe", 1);
    ele_Xe128->AddIsotope(iso_Xe128, 100.*perCent);
    G4Element  *ele_Xe129 = new G4Element("element Xe129", "Xe", 1);
    ele_Xe129->AddIsotope(iso_Xe129, 100.*perCent);
    G4Element  *ele_Xe130 = new G4Element("element Xe130", "Xe", 1);
    ele_Xe130->AddIsotope(iso_Xe130, 100.*perCent);
    G4Element  *ele_Xe131 = new G4Element("element Xe131", "Xe", 1);
    ele_Xe131->AddIsotope(iso_Xe131, 100.*perCent);
    G4Element  *ele_Xe132 = new G4Element("element Xe132", "Xe", 1);
    ele_Xe132->AddIsotope(iso_Xe132, 100.*perCent);
    G4Element  *ele_Xe134 = new G4Element("element Xe134", "Xe", 1);
    ele_Xe134->AddIsotope(iso_Xe134, 100.*perCent);
    G4Element  *ele_Xe136 = new G4Element("element Xe136", "Xe", 1);
    ele_Xe136->AddIsotope(iso_Xe136, 100.*perCent);

    //-------------------------------------- stainless steel --------------------------------------
    /*G4Material *SS304LSteel = new G4Material("SS304LSteel", 8.00*g/cm3, 5, kStateSolid);
    SS304LSteel->AddElement(Fe, 0.65);
    SS304LSteel->AddElement(Cr, 0.20);
    SS304LSteel->AddElement(Ni, 0.12);
    SS304LSteel->AddElement(Mn, 0.02);
    SS304LSteel->AddElement(Si, 0.01);*/

    //------------------------------------------ teflon -------------------------------------------
    // G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
    // Teflon->AddElement(C, 0.240183);
    // Teflon->AddElement(F, 0.759817);
    
}

//__________________________________________________________________________________________________________

void
DetectorConstruction::DefineGeometryParameters()
{
    //======================================== Laboratory =========================================
    m_hGeometryParameters["LabSize"]    = 250   * cm;
    //============================================ TPC ============================================
    m_hGeometryParameters["TPC_outerR"] =  65.5 * cm;
    m_hGeometryParameters["TPC_Z"]      = 151   * cm;
    //============================================ LXe ============================================
    m_hGeometryParameters["LXe_outerR"] =  65   * cm;
    m_hGeometryParameters["LXe_Z"]      = 150   * cm;
    //============================================= FV ============================================
    m_hGeometryParameters["FV_outerR"] =  57   * cm;
    m_hGeometryParameters["FV_Z"]      = 134   * cm;
}

//__________________________________________________________________________________________________________

G4double
DetectorConstruction::GetGeometryParameter(const char *szParameter)
{
    return m_hGeometryParameters[szParameter];
}

//__________________________________________________________________________________________________________

void
DetectorConstruction::ConstructLaboratory()
{
    //======================================== Laboratory =========================================
    const G4double dLabHalfSize = 0.5*GetGeometryParameter("LabSize");
    
    G4Material *Air = G4Material::GetMaterial("G4_AIR");

    G4Box *pLabBox         = new G4Box("LabBox", dLabHalfSize, dLabHalfSize, dLabHalfSize);
    m_pLabLogicalVolume    = new G4LogicalVolume(pLabBox, Air, "LabVolume", 0, 0, 0);
    m_pLabPhysicalVolume   = new G4PVPlacement(0, G4ThreeVector(), m_pLabLogicalVolume, "Lab", 0, false, 0);
    m_pMotherLogicalVolume = m_pLabLogicalVolume;
    
    m_pLabLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

}

//__________________________________________________________________________________________________________

void
DetectorConstruction::ConstructSD()
{
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    // make LXe sensitive
    G4VSensitiveDetector* LXe_SD = new SensitiveDetector(SDname="/LXe");
    SDman->AddNewDetector(LXe_SD);
    m_pLXe_LogicalVolume->SetSensitiveDetector(LXe_SD);
    m_pTPC_LogicalVolume->SetSensitiveDetector(LXe_SD);

}

//__________________________________________________________________________________________________________

void
DetectorConstruction::ConstructCollimatorSystem()
{
    
    //========================================= Materials =========================================
    G4Material *LXe    = G4Material::GetMaterial("LXe");
    //G4Material *SSteel = G4Material::GetMaterial("SS304LSteel");

    // rotation matrix for the TPC
    G4RotationMatrix *pRot      = new G4RotationMatrix();
    G4RotationMatrix *pRot_idle = new G4RotationMatrix();
    //pRot->rotateY(90.*deg);
       
    //================================ Construction Xenon and TPC =================================
    const G4double dTPCouterRadius =     GetGeometryParameter("TPC_outerR");
    const G4double dTPCHalfZ       = 0.5*GetGeometryParameter("TPC_Z");
    const G4double dLXeouterRadius =     GetGeometryParameter("LXe_outerR");
    const G4double dLXeHalfZ       = 0.5*GetGeometryParameter("LXe_Z");
    
    G4Tubs *pTPC   			= new G4Tubs("TPC",	0.*cm,
												dTPCouterRadius,
												dTPCHalfZ,
												0.*deg, 360.*deg);
    //m_pTPC_LogicalVolume	= new G4LogicalVolume(pTPC, SSteel, "TPCTUBS", 0, 0, 0);
    m_pTPC_LogicalVolume	= new G4LogicalVolume(pTPC, LXe, "TPCTUBS", 0, 0, 0);
    m_pTPC_PhysicalVolume 	= new G4PVPlacement(pRot, G4ThreeVector(0,0,0),
                                                        m_pTPC_LogicalVolume,
                                                        "TPC",
                                                        m_pMotherLogicalVolume,
                                                        false, 0);

    G4Tubs *pLXe_volume     = new G4Tubs("LXe",	0.*cm,
    											dLXeouterRadius,
    											dLXeHalfZ,
    											0.*deg, 360.*deg);
    m_pLXe_LogicalVolume 	= new G4LogicalVolume(pLXe_volume, LXe, "LXe", 0, 0, 0);
    m_pLXe_PhysicalVolume	= new G4PVPlacement(pRot_idle, G4ThreeVector(0,0,0),
                                                        m_pLXe_LogicalVolume,
                                                        "LXe",
                                                        m_pTPC_LogicalVolume,
                                                        false, 0);
                                                

    G4cout << "end sysdef" << G4endl;
}

//__________________________________________________________________________________________________________

void
DetectorConstruction::PrintGeometryInformation()
{
}

//__________________________________________________________________________________________________________

void DetectorConstruction::MakeDetectorPlots()
{
    _fGeom    = new TFile(detRootFile,"RECREATE");
    _detector = _fGeom->mkdir("detector");
    
    // store properties of the materials that are used
    StoreMaterialParameters();
    // store geometry parameters
    StoreGeometryParameters();
    
    _fGeom->Write();
    _fGeom->Close();
    
}

//__________________________________________________________________________________________________________

void DetectorConstruction::StoreGeometryParameters()
{
    // TDirectory for storage of teh geometry parameters
    TDirectory *_geometry = _detector->mkdir("geometry");
    _geometry->cd();
        
    _fGeom->cd();
}

//__________________________________________________________________________________________________________

void DetectorConstruction::StoreMaterialParameters()
{
    
    // make a list of materials for graphs
    G4int nmaterial = G4Material::GetNumberOfMaterials();
    G4cout << "MakeDetectorPlots:: Number of materials = " << nmaterial << G4endl;
    
    TDirectory *_materials = _detector->mkdir("materials");
    _materials->cd();
    

    vector<TDirectory*> matdirs;
    
    for(G4int imat=0; imat<nmaterial; imat++){
        G4Material *mat  = G4NistManager::Instance()->GetMaterial(imat);
        G4String matname = mat->GetName();
        G4cout << "MakeDetectorPlots:: " << matname << G4endl;
        G4double T       = mat->GetTemperature();
        G4double rho     = mat->GetDensity();
        G4double P       = mat->GetPressure();
        
        matdirs.push_back(_materials->mkdir(matname));
        matdirs[imat]->cd();
        TParameter<double> *TemperaturePar = new TParameter<double>("Temperature",T);
        TemperaturePar->Write();
        TParameter<double> *DensityPar     = new TParameter<double>("Density",rho / (g/cm3));
        DensityPar->Write();
        TParameter<double> *PressurePar    = new TParameter<double>("Pressure",P / bar);
        PressurePar->Write();
        // disect the material
        size_t nele = mat->GetNumberOfElements();
        G4ElementVector *elems       = (G4ElementVector*)mat->GetElementVector();
        G4double        *fractionVec = (G4double*)mat->GetFractionVector();
        
        for(size_t iele = 0; iele < nele; iele++){
            G4String elname = (*elems)[iele]->GetName();
            G4double frac   = fractionVec[iele];
            //      G4cout <<iele<<" elem = "<<(*elems)[iele]->GetName()<<" f = "<<fractionVec[iele]<<G4endl;
            char  elFrac[100];
            sprintf(elFrac,"f_%s",elname.c_str());
            TParameter<double> *_fracPar = new TParameter<double>((G4String)elFrac,frac);
            _fracPar->Write();
        }
        
        
        _materials->cd();
    }
    
    _fGeom->cd();
}

//__________________________________________________________________________________________________________