#ifndef _FMC_
#define _FMC_

#include <globals.hh>
#include <G4ThreeVector.hh>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>

class FMC
{
public:
    FMC();

public:
    std::vector<G4double> intersect_side(G4ThreeVector x0, G4ThreeVector p, G4bool &keepevent, G4double r_cyl, G4double z_cyl);
    std::vector<G4double> intersect_plane(G4ThreeVector x0, G4ThreeVector p, G4bool &keepevent, G4double r_cyl, G4double z_cyl);
    std::vector<G4double> sort_vector(std::vector<G4double> s1, std::vector<G4double> s2);

};

#endif