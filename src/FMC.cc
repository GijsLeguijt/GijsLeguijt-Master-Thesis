/*
 * Fast Monte Carlo-code. Reduces the number of events that needs to be simulated,
 * by forcing particles to show the wanted behaviour. Statistics are taken care of
 * by adding weights to the particles.
 */

#include "FMC.hh"
#include "G4UImanager.hh"
#include "DetectorConstruction.hh"


using namespace CLHEP;

//__________________________________________________________________________________________________________

FMC::FMC()
{

}

//__________________________________________________________________________________________________________

vector<G4double> FMC::intersect_side(G4ThreeVector x0, G4ThreeVector p, G4bool &keepevent, G4double r_cyl, G4double z_cyl)
{   /*
     * Checks whether a trajectory intersects with the side (NOT TOP/BOTTOM) of a cylinder.
     * Returns the path length to the closest intersection point
     * 
     * Trajectory is parametrised with x0 + s * p
     */

    vector<G4double> s_int =  {0,0}; // Path length to the intersection points

    
    // Checks intersection with circle
    G4double A = pow(p.x(),2) + pow(p.y(),2);
    G4double B = 2 * (x0.x() * p.x() + x0.y() * p.y());
    G4double C = pow(x0.x(),2) + pow(x0.y(),2) - pow(r_cyl,2);

    G4double discriminant = pow(B,2) - 4 * A * C;

    vector<G4int> sign = {1,-1};

    // Checks whether the height of the intersection is within the cylinder
    if (discriminant >= 0)
    {   for (G4int i = 0; i < 2; i++)
        {
            G4double s  = (-B + sign[i] * sqrt(discriminant)) / (2 * A);
            G4double z  = x0.z() + s * p.z();
            if (abs(z) <= z_cyl && s > 0)
            {   
                // Keep the particle to track
                keepevent = true;
                
                if (!s_int[0])
                {   s_int[0] = s;
                }
                else
                {   s_int[1] = s;
                }
            }
        }
    }

    return s_int;
}

//__________________________________________________________________________________________________________

vector<G4double> FMC::intersect_plane(G4ThreeVector x0, G4ThreeVector p, G4bool &keepevent, G4double r_cyl, G4double z_cyl)
{   /*
     * Checks whether a trajectory intersects with the top or bottom (NOT SIDE) of a cylinder.
     * Returns the path length to the closest intersection point
     * 
     * Trajectory is parametrised with x0 + s * p
     */

    vector<G4int> sign = {1,-1};

    vector<G4double> s_int =  {0,0}; // Path length to the intersection points

    for (G4int i = 0; i < 2; i++)
    {   // Path length to the top(/bottom) of the cylinder
        G4double s = (sign[i] * z_cyl - x0.z()) / p.z();

        G4double x = x0.x() + s * p.x();
        G4double y = x0.y() + s * p.y();
        G4double r = sqrt( pow(x,2) + pow(y,2) );

        // Checks whether track, at correct z, lies within the circle
        if (r <= r_cyl && s > 0)
        {
            // Keep the particle to track
            keepevent = true;
         
            if (!s_int[0])
            {   s_int[0] = s;
            }
            else
            {   s_int[1] = s;
            }
            
        }
    }

    return s_int;
}

//__________________________________________________________________________________________________________

vector<G4double> FMC::sort_vector(vector<G4double> s1, vector<G4double> s2)
{   /*
     * Combines two arrays of path lengths, sorts them (descending) and removes the zeros
     */

    s1.insert(end(s1), begin(s2), end(s2));
	sort(s1.begin(), s1.end(), std::greater<G4double>());

	s1.erase(remove(s1.begin(), s1.end(), 0), s1.end());
    
    return s1;
}

//__________________________________________________________________________________________________________