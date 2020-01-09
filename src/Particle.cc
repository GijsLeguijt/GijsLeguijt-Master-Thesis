#include "Particle.hh"
#include <G4NistManager.hh>
#include "DetectorConstruction.hh"
#include <G4EmCalculator.hh>
#include <cmath>
#include <G4Material.hh>
#include <RunAction.hh>
#include <Randomize.hh>
#include <G4KleinNishinaModel.hh>
#include <G4MaterialCutsCouple.hh>
#include <G4DynamicParticle.hh>
#include <G4Track.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTable.hh>
#include <G4DataVector.hh>

#include "G4ParticleChangeForGamma.hh"
#include "G4VEmModel.hh"
#include "G4ProductionCutsTable.hh"

//__________________________________________________________________________________________________________

Particle::Particle()
{   
}

//__________________________________________________________________________________________________________

void Particle::Print(G4bool extended) //default = false
{   
    G4cout << "particle::print PARTICLE STATUS"                              << G4endl
           << "particle::print type               = " << m_type              << G4endl
           << "particle::print energy             = " << m_energy  << " MeV" << G4endl
           << "particle::print deposited energy   = " << m_edep    << " MeV" << G4endl
           << "particle::print origin             = " << m_x0start << " mm"  << G4endl
           << "particle::print final position     = " << m_x0      << " mm"  << G4endl
           << "particle::print direction          = " << m_direction         << G4endl
           << "particle::print weight             = " << m_weight            << G4endl;
    if (extended)
    {
    G4cout << "particle::print variance reduction = " << m_vrt               << G4endl
           << "particle::print n scatter max      = " << m_nscatter_max      << G4endl
           << "particle::print e deposit max      = " << m_edep_max          << G4endl;
    }
}

//__________________________________________________________________________________________________________

vector<G4double> Particle::intersect(G4double cyl_outerRadius, G4double cyl_halfZ)
{   /*
     * Checks whether a trajectory intersects with a cylinder. Split into intersections
     * with the side and with the top/bottom. Returns the path lengths to the inter-
     * sections in ascending order
     */

    vector<G4double> s_side  = intersect_side( m_x0, m_direction, cyl_outerRadius, cyl_halfZ);
    vector<G4double> s_plane = intersect_plane(m_x0, m_direction, cyl_outerRadius, cyl_halfZ);       
    
    vector<G4double> s       = sort_vector(s_side, s_plane);

    return s;
}

//__________________________________________________________________________________________________________

vector<G4double> Particle::intersect_side(G4ThreeVector x0, G4ThreeVector p, G4double r_cyl, G4double z_cyl)
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
            {   // We have an interaction   
                if (!s_int[0]) // 1st spot empty? save there
                {   s_int[0] = s;
                }
                else // Save on 2nd spot
                {   s_int[1] = s;
                }
            }
        }
    }

    return s_int;
}

//__________________________________________________________________________________________________________

vector<G4double> Particle::intersect_plane(G4ThreeVector x0, G4ThreeVector p, G4double r_cyl, G4double z_cyl)
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
        {   // We have an interaction
            if (!s_int[0]) // 1st spot empty? save there
            {   s_int[0] = s;
            }
            else // Save on 2nd spot
            {   s_int[1] = s;
            }
            
        }
    }

    return s_int;
}

//__________________________________________________________________________________________________________

vector<G4double> Particle::sort_vector(vector<G4double> s1, vector<G4double> s2)
{   /*
     * Combines two arrays of path lengths, sorts them (ascending) and removes the zeros
     */

    s1.insert(end(s1), begin(s2), end(s2));               //Two arrays glued together
	sort(s1.begin(), s1.end());                           //Sorted (ascending)

	s1.erase(remove(s1.begin(), s1.end(), 0), s1.end());  //Removes zeros
    
    return s1;
}

//__________________________________________________________________________________________________________

G4double Particle::Get_att_probability(G4double distance)
{   /*
     * Get the probability for a particle of energy "m_energy", to travel distance
     */

    G4double mu = emCalc.ComputeGammaAttenuationLength(m_energy, m_material);

    G4double prob = exp(- distance / mu);

    return prob;
}

//__________________________________________________________________________________________________________


void Particle::Propagate()
{   /*
     * Propagates the particle
     */

    // Get the dimensions of the cylinders, the LXe-volume and the fiducial volume (FV)
    const G4double LXeouterRadius =       DetectorConstruction::GetGeometryParameter("LXe_outerR");
    const G4double LXeHalfZ       = 0.5 * DetectorConstruction::GetGeometryParameter("LXe_Z");
    const G4double FVouterRadius  =       DetectorConstruction::GetGeometryParameter("FV_outerR");
    const G4double FVHalfZ        = 0.5 * DetectorConstruction::GetGeometryParameter("FV_Z");
    
    G4bool terminate = false;
    G4double s_max;
    m_nscatter = 0; // redundant?

    if (m_debug)
    {
        G4cout << "particle::propagate Next event" << G4endl;
    }

    while (!terminate) // While terminate is false
    {   
        // intersections with the LXe
        vector<G4double> cryostat_intersections = intersect(LXeouterRadius, LXeHalfZ);

        if (cryostat_intersections.size() > 0) // if we intersect with the LXe
        {   
            if (cryostat_intersections.size() == 2) // particle is still outside the LXe
            {   s_max = cryostat_intersections[1];  // maximum length until trajectory leaves the LXe
            }
            else                                    // particle is already inside the LXe
            {   s_max = cryostat_intersections[0];  // maximum length until trajectory leaves the LXe                
            }
            
            G4double s_max_fiducial = -1.0;

            if (m_vrt == "fiducial_scatter") // if we require scatters in the FV
            {
                if (m_nscatter == 0) // if particle didn't scatter yet
                {   
                    if (m_debug)
                    {
                        G4cout << "particle::propagate VRT - transport to fiducial volume is ON" << G4endl;
                    }

                    // intersections with the FV
                    vector<G4double> fiducial_intersections = intersect(FVouterRadius, FVHalfZ);    
                    
                    if (fiducial_intersections.size() > 0) // if we intersect with the FV
                    {
                        G4double s_fiducial_entry = fiducial_intersections[0];
                        G4double s_fiducial_exit  = fiducial_intersections[1];

                        // update the weight for forcing to reach the FV
                        m_weight *= Get_att_probability(s_fiducial_entry);
                        Update_particle(s_fiducial_entry,"transport");

                        // maximal s will be lower, as we reached the FV, subtract the travelled distance
                        s_max -= s_fiducial_entry;
                        s_max_fiducial = s_fiducial_exit-s_fiducial_entry; // distance of path through FV
                    }
                    else // we do not intersect with the FV
                    {
                        terminate = true;
                        continue;
                    }
                }
                else if (m_nscatter > 0 && m_nscatter < m_nscatter_max) // particle scattered, but allowed to scatter more
                {
                    // intersections with the FV
                    vector<G4double> fiducial_intersections = intersect(FVouterRadius, FVHalfZ);
                    
                    if (fiducial_intersections.size() > 0) // if we intersect with the FV
                    {   // particle is already in the FV, so will find 1 (and only 1) intersection with the FV
                        s_max_fiducial = fiducial_intersections[0];
                    }
                    else // shouldn't occur, as particle is in the FV, so should intersect
                    {
                        G4cout << "particle::propagate ERROR.... bad intersection. discard event." << G4endl;
                        terminate = true;
                        continue;
                    }
                }
                else if (m_nscatter == m_nscatter_max) // particle not allowed to scatter more, calculate chance to exit
                {        
                    //s_max = cryostat_intersections[0]; // was in python code, seems redundant

                    m_weight *= Get_att_probability(s_max); // chance to leave the detector without more interactions
                        
                    terminate = true;
                    continue;
                }
            }
            
            // generate a distance to the next interaction, s_max_fiducial is maximal distance to interaction
            // if we do not require interaction within the FV, s_max_fiducial = -1, and ignored
            G4double s_gen = Generate_interaction_point(s_max_fiducial);

            if(s_gen < s_max) // interaction within the detector
            {   // scatter the particle, either PE or compton
                std::string process = Update_particle(s_gen);
                 
                if (process == "pho") // energy deposit is too high, not interested
                {
                    terminate = true;
                }
            }
            else // no interaction within the detector
            {
                terminate = true;
            }
        }
        else // no intersection of trajectory with the LXe
        {
            terminate = true;
        }
        

    }

    if (m_debug)
    {
        G4cout << "particle::propagate exit propagator" << G4endl;
    }
}

//__________________________________________________________________________________________________________

G4double Particle::Generate_interaction_point(G4double smax) //default for smax = -1
{   /*
     * Picks the length of the track until the next interaction
     * argument smax = maximal path length to generate, default is "-1" (= infinite)
     */

    // Get the attenuation length for the photon
    G4double   mu = emCalc.ComputeGammaAttenuationLength(m_energy, m_material);
    G4double rmax = 1; // sets the range for the chosen random number, updated if there is a limit on s
    
    if (smax > 0.0) // calculate the range to fulfill smax, also update the weight
    {   rmax = 1.0 - exp(- smax / mu);
        m_weight *= rmax;
    }

    // generate the path length
    G4double r = Random_uniform(0.0, rmax);
    G4double L = - log(1-r) * mu;

    return L;
}

//__________________________________________________________________________________________________________

std::string Particle::Update_particle(G4double s_scatter, std::string process)//Default of process = ""
{   /*
     * Combination of "scatter" and "update_particle" from python code
     * 
     * Will transport a "time" s_scatter and scatter via inc or PE, unless process equals
     * transport, in which case only the transportation takes place
     */

    if (process == "") //if no process pre-selected, get either compton or PE
    {   process = Select_scatter_process();
    }

    m_x0 += s_scatter * m_direction; //get to the new position
    G4double edep;

    if (process == "inc")
    {   //Inverse Compton scattering
        edep = Do_compton();

        Save_interaction(m_x0, edep, process); //save interaction data
    }
    else if (process == "pho")
    {   //Photo-electric effect, all energy is deposited
        edep        = m_energy;
        m_nscatter += 1;
        m_edep     += edep;
        m_energy    = 0;

        Save_interaction(m_x0, edep, process); //save interaction data
    }
    else if (process == "transport")
    {   //Just transportation, nothing additional happens
    }
    else
    {
        G4cout << "Invalid process name, leave empty or use transport, pho or inc" << G4endl;
    }
    
    return process;
}

//__________________________________________________________________________________________________________

std::string Particle::Select_scatter_process()
{   /*
     * Selects the process based on the relative cross-sections
     * "inc" = inverse compton
     * "pho" = photo-electric effect (PE)
     */

    std::string process;

    // Cross-sections for compton and PE, units don't matter as we need a fraction    
    G4double sigma_compt = emCalc.ComputeCrossSectionPerVolume(m_energy,m_type,"compt",m_material->GetName(),0);
    G4double sigma_phot  = emCalc.ComputeCrossSectionPerVolume(m_energy,m_type,"phot", m_material->GetName(),0);
    G4double sigma_total = sigma_compt + sigma_phot;
    
    G4double frac        = sigma_compt / sigma_total;

    if (m_edep_max < m_energy)
    // if the maximum energy deposit is smaller than the total kinetic energy, PE is off and the weight updated
    {    
        process   = "inc";
        m_weight *= frac;
    }
    else
    {
        G4double r = Random_uniform(0.0, 1.0); // select process based on relative cross-section
        
        if (r < frac)
        {   process  = "inc";
        }
        else
        {    process = "pho";
        }

        
        if (m_vrt == "fiducial_scatter" && m_nscatter < m_nscatter_max - 1)
        // if we require multiple scatters, PE is not allowed until required number of scatters is reached
        {   process   = "inc";
            m_weight *= frac;
        }
    }

    return process;
}

//__________________________________________________________________________________________________________

G4double Particle::Random_uniform(G4double min, G4double max)
{   /*
     * Returns a number randomly picked from a uniform distribution on interval [min, max]
     */
    
    G4double rand = G4UniformRand();
    
    return min + (max - min) * rand;
}

//__________________________________________________________________________________________________________

G4double Particle::Do_compton()
{   /*
     * look into this!
     * 
     * Does compton scattering using G4KleinNishinaModel
     * 
     * Initial part is declaration of (irrelevant) objects needed for G4KleinNishinaModel to work
     */

    const G4ParticleDefinition* part_def = G4ParticleTable::GetParticleTable()->FindParticle(m_type);

    //
    //Make the photon that will scatter + a storage for the scattered electron
    //

    G4DynamicParticle * ppart_dyn = new G4DynamicParticle(part_def,m_direction,m_energy);
    std::vector< G4DynamicParticle * >    vec_elec = {};
    std::vector< G4DynamicParticle * > * pvec_elec = &vec_elec;
    
    //
    //G4KleinNishinaModel's SampleSecondaries requires following 2 variables but doesn't use them
    //
    G4double tmin      = 0;    
    G4double maxEnergy = 0; 
    

    //G4MaterialCutsCouple    mat_cuts = G4MaterialCutsCouple(m_material);
    //G4MaterialCutsCouple * pmat_cuts = &mat_cuts; //pointer to material cuts
    //pmat_cuts->SetIndex(0);//is NOT how to do this, but otherwise it tries to get -1th element of list -> look into this!
    const G4MaterialCutsCouple * pmat_cuts = G4ProductionCutsTable::GetProductionCutsTable()->GetMaterialCutsCouple(0); //wrong material cut -> look into this!
    
    size_t x = 1; //think this is: number of cuts, can't be zero as compiler will crash, trying to find something that's not there -> look into this!
    G4DataVector cuts = {x, 0}; //thought it was: number of cuts, energy cuts, but 2nd number doesnt change anything -> look into this!


    G4ParticleChangeForGamma    Gammainfo;                          //storage for new direction and energy
    G4ParticleChangeForGamma * pGammainfo = &Gammainfo;
    
    //
    // Actual KN model, initialisation and afterwards compton scattering
    //
    G4KleinNishinaModel KNmodel = G4KleinNishinaModel();
    KNmodel.SetParticleChange(pGammainfo);                          //make sure it uses correct storage
    KNmodel.Initialise(part_def,cuts);
    
    //do the scattering
    KNmodel.SampleSecondaries(pvec_elec, pmat_cuts, ppart_dyn, tmin, maxEnergy);

    //Hook to the scattered electron, in case it's needed
    //G4DynamicParticle * electron = pvec_elec->at(0);    

    //saving data to the particle    
    m_direction         = Gammainfo.GetProposedMomentumDirection();
    G4double new_energy = Gammainfo.GetProposedKineticEnergy();

    G4double edep = (m_energy - new_energy);
    m_edep     += edep;
    m_energy    = new_energy;          
    m_nscatter += 1;

    return edep;
}

//__________________________________________________________________________________________________________

void Particle::Save_interaction(G4ThreeVector position, G4double deposit, std::string process)
{   /*
     * Saves position, deposited energy and the process of an interaction
     */

    m_pos_int.push_back(position);
    m_edep_int.push_back(deposit);
    m_pro_int.push_back(process);
}


/* To do:
Zoeken op "look into this!"



Zoeken op: "LXe" en zorgen dat dit correcte materiaal is (ook in .hh)
Zorgen bij aanroepen functie dat elk deeltje nieuw seed krijgt + 1e zelfde als run messenger
Zoeken op "std::string Particle::Select_scatter_process()", compton vs inverse?
Zoeken op "redundant"
*/