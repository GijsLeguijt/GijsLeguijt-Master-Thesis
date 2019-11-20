#ifndef _Particle_
#define _Particle_

#include <globals.hh>
#include <G4ThreeVector.hh>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4EmCalculator.hh>

using namespace CLHEP;

class Particle
{
public:
    Particle();
    void Print();
    std::vector<G4double> intersect(G4double cyl_outerRadius, G4double cyl_halfZ);
    G4double Get_att_probability(G4double distance);
    void Propagate();
    G4double Generate_interaction_point(G4double smax = -1);
    std::string Scatter(G4double s);
    std::string Select_scatter_process();
private:
    std::vector<G4double> intersect_side(G4ThreeVector x0, G4ThreeVector p, G4double r_cyl, G4double z_cyl);
    std::vector<G4double> intersect_plane(G4ThreeVector x0, G4ThreeVector p, G4double r_cyl, G4double z_cyl);
    std::vector<G4double> sort_vector(std::vector<G4double> s1, std::vector<G4double> s2);
    G4double Random_uniform(G4double min, G4double max);
    G4EmCalculator emCalc;
    
private:
    //Parameters
    std::string           m_type         = "gamma";
    std::string           m_vrt          = "";
    G4bool                m_debug        = false;
    G4double              m_weight       = 1;
    G4double              m_energy       = 0;
    G4double              m_edep         = 0;
    G4double              m_edep_max     = 100000 ;//* MeV;
    G4int                 m_nscatter     = 0;
    G4int                 m_nscatter_max = 1;
    G4ThreeVector         m_x0           = G4ThreeVector(0, 0, 0);
    G4ThreeVector         m_x0start      = G4ThreeVector(0, 0, 0);
    G4ThreeVector         m_direction    = G4ThreeVector(0, 0, 0);
    std::vector<G4double> m_xint         = {0,0};
    std::vector<G4double> m_sint         = {0,0};
    G4Material *          m_material     = G4Material::GetMaterial("LXe");
  
public:
    //Sets
    void setType(std::string type)             {m_type         = type;}
    void setVrt(std::string vrt)               {m_vrt          = vrt;}
    void setDebug(G4bool debug)                {m_debug        = debug;}
    void setWeight(G4double weight)            {m_weight       = weight;}
    void setEnergy(G4double energy)            {m_energy       = energy;}
    void setEdep(G4double edep)                {m_edep         = edep;}
    void setEdep_max(G4double edep_max)        {m_edep_max     = edep_max;}
    void setNscatter(G4int nscatter)           {m_nscatter     = nscatter;}
    void setNscatter_max(G4int nscatter_max)   {m_nscatter_max = nscatter_max;}
    void setX0(G4ThreeVector x0)               {m_x0           = x0;}
    void setX0start(G4ThreeVector x0start)     {m_x0start      = x0start;}
    void setdirection(G4ThreeVector direction) {m_direction    = direction;}
    void setXint(std::vector<G4double> xint)   {m_xint         = xint;}
    void setSint(std::vector<G4double> sint)   {m_sint         = sint;}
    void setMaterial(G4Material * material)    {m_material     = material;}
    
    
    //Gets
    std::string           getType()         { return m_type; }
    std::string           getVrt()          { return m_vrt; }
    G4bool                getDebug()        { return m_debug; }
    G4double              getWeight()       { return m_weight; }
    G4double              getEnergy()       { return m_energy; }
    G4double              getEdep()         { return m_edep; }
    G4double              getEdep_max()     { return m_edep_max; }    
    G4int                 getNscatter()     { return m_nscatter; }
    G4int                 getNscatter_max() { return m_nscatter_max; }
    G4ThreeVector         getX0()           { return m_x0; }
    G4ThreeVector         getX0start()      { return m_x0start; }
    G4ThreeVector         getdirection()    { return m_direction; }
    std::vector<G4double> getXint()         { return m_xint; }
    std::vector<G4double> getSint()         { return m_sint; }
    G4Material *          getMaterial()     { return m_material; }
    
};

#endif