#include "EventData.hh"

//__________________________________________________________________________________________________________

EventData::EventData()
{
    m_iEventId = 0;                               m_iEventId_vrt = 0;
    
    m_fTotalEnergyDeposited = 0.;
    m_iNbSteps = 0;
    
    m_pEtot                = new vector<float>;   
    m_pTrackId             = new vector<int>;     
    m_pParentId            = new vector<int>;     
    m_pCollectionId        = new vector<int>;     
    m_pParticleType        = new vector<string>;  
    m_pParentType          = new vector<string>;  
    m_pCreatorProcess      = new vector<string>;  
    m_pDepositingProcess   = new vector<string>;  m_pDepositingProcess_vrt   = new vector<string>;
    m_pX                   = new vector<float>;   m_pX_vrt                   = new vector<float>;
    m_pY                   = new vector<float>;   m_pY_vrt                   = new vector<float>;
    m_pZ                   = new vector<float>;   m_pZ_vrt                   = new vector<float>;
    m_pEnergyDeposited     = new vector<float>;   m_pEnergyDeposited_vrt     = new vector<float>;
    m_pKineticEnergy       = new vector<float>;   
    m_pTime                = new vector<float>; 
    
    m_pPrimaryParticleType = new vector<string>;  
    m_fPrimaryX = 0.;                             m_fPrimaryX_vrt            = 0.;
    m_fPrimaryY = 0.;                             m_fPrimaryY_vrt            = 0.;
    m_fPrimaryZ = 0.;                             m_fPrimaryZ_vrt            = 0.;
    m_fPrimaryE = 0.;                             m_fPrimaryE_vrt            = 0.;
                                                  m_fPrimaryW_vrt            = 0.; 
    
}

//__________________________________________________________________________________________________________

EventData::~EventData()
{
    delete m_pEtot;                 
    delete m_pTrackId;              
    delete m_pParentId;             
    delete m_pCollectionId;         
    delete m_pParticleType;         
    delete m_pParentType;           
    delete m_pCreatorProcess;       
    delete m_pDepositingProcess;    delete m_pDepositingProcess_vrt;
    delete m_pX;                    delete m_pX_vrt;
    delete m_pY;                    delete m_pY_vrt;
    delete m_pZ;                    delete m_pZ_vrt;
    delete m_pEnergyDeposited;      delete m_pEnergyDeposited_vrt;
    delete m_pKineticEnergy;        
    delete m_pTime;                 
                                    
    delete m_pPrimaryParticleType;  
    
}

//__________________________________________________________________________________________________________

void
EventData::Clear()
{
    m_iEventId = 0;                     m_iEventId_vrt = 0; 
    
    m_fTotalEnergyDeposited = 0.0;
    m_iNbSteps = 0;
    
    m_pEtot->clear();                   
    m_pTrackId->clear();
    m_pParentId->clear();
    m_pCollectionId->clear();
    m_pParticleType->clear();           
    m_pParentType->clear();             
    m_pCreatorProcess->clear();         
    m_pDepositingProcess->clear();      m_pDepositingProcess_vrt->clear();
    m_pX->clear();                      m_pX_vrt->clear();
    m_pY->clear();                      m_pY_vrt->clear();
    m_pZ->clear();                      m_pZ_vrt->clear();
    m_pEnergyDeposited->clear();        m_pEnergyDeposited_vrt->clear();
    m_pKineticEnergy->clear();          
    m_pTime->clear();
    
    m_pPrimaryParticleType->clear();    
    m_fPrimaryX = 0.;                   m_fPrimaryX_vrt = 0.;
    m_fPrimaryY = 0.;                   m_fPrimaryY_vrt = 0.;
    m_fPrimaryZ = 0.;                   m_fPrimaryZ_vrt = 0.;	
    m_fPrimaryE = 0.;                   m_fPrimaryE_vrt = 0.;
    m_fPrimaryW = 0.;                   m_fPrimaryW_vrt = 0.;
    
}

//__________________________________________________________________________________________________________