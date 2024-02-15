#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <vector>
#include <cassert>
#include <math.h>
#include <chrono>

const double deg2rad = acos(-1)/180.0; // pi/180 for changing degs to radians
double accumulated_forces_bond  = 0.;     // Checksum: accumulated size of forces
double accumulated_forces_angle = 0.;     // Checksum: accumulated size of forces
double accumulated_forces_non_bond = 0.;  // Checksum: accumulated size of forces


/* a class for the bond between two atoms U = 0.5k(r12-L0)^2 */
class Bond {
public:
    double K;    // force constant
    double L0;   // relaxed length
    int a1, a2;  // the indexes of the atoms at either end
};

/* a class for the angle between three atoms  U=0.5K(phi123-phi0)^2 */
class Angle {
public:
    double K;
    double Phi0;
    int a1, a2, a3; // the indexes of the three atoms, with a2 being the centre atom
};

// ===============================================================================
// Two new classes arranging Atoms in a Structure-of-Array data structure
// ===============================================================================

/* atom class, represent N instances of identical atoms */
class Atoms {
public:
    // The mass of the atom in (U)
    double mass;
    double ep;            // epsilon for LJ potential
    double sigma;         // Sigma, somehow the size of the atom
    double charge;        // charge of the atom (partial charge)
    std::string name;     // Name of the atom
    // the position in (nm), velocity (nm/ps) and forces (k_BT/nm) of the atom
    std::vector<double> px,py,pz,vx,vy,vz,fx,fy,fz;
    // constructor, takes parameters and allocates p, v and f properly to have N_identical elements
    Atoms(double mass, double ep, double sigma, double charge, std::string name, size_t N_identical) 
    : mass{mass}, ep{ep}, sigma{sigma}, charge{charge}, name{name}, 
    px(N_identical),py(N_identical),pz(N_identical), vx(N_identical), vy(N_identical), vz(N_identical), fx(N_identical), fy(N_identical), fz(N_identical)
    {}
};

/* molecule class for no_mol identical molecules */
class Molecules {
public:
    std::vector<Atoms> atoms;         // list of atoms in the N identical molecule
    std::vector<Bond> bonds;          // the bond potentials, eg for water the left and right bonds
    std::vector<Angle> angles;        // the angle potentials, for water just the single one, but keep it a list for generality
    int no_mol;
};

// ===============================================================================


/* system class */
class System {
public:
    Molecules molecules;          // all the molecules in the system
    double time = 0;                          // current simulation time
};

class Sim_Configuration {
public:
    int steps = 10000;     // number of steps
    int no_mol = 4;        // number of molecules
    double dt = 0.0005;    // integrator time step
    int data_period = 100; // how often to save coordinate to trajectory
    std::string filename = "trajectory.txt";   // name of the output file with trajectory
    // system box size. for this code these values are only used for vmd, but in general md codes, period boundary conditions exist

    // simulation configurations: number of step, number of the molecules in the system, 
    // IO frequency, time step and file name
    Sim_Configuration(std::vector <std::string> argument){
        for (long unsigned int i = 1; i<argument.size() ; i += 2){
            std::string arg = argument.at(i);
            if(arg=="-h"){ // Write help
                std::cout << "MD -steps <number of steps> -no_mol <number of molecules>"
                          << " -fwrite <io frequency> -dt <size of timestep> -ofile <filename> \n";
                exit(0);
                break;
            } else if(arg=="-steps"){
                steps = std::stoi(argument[i+1]);
            } else if(arg=="-no_mol"){
                no_mol = std::stoi(argument[i+1]);
            } else if(arg=="-fwrite"){
                data_period = std::stoi(argument[i+1]);
            } else if(arg=="-dt"){
                dt = std::stof(argument[i+1]);
            } else if(arg=="-ofile"){
                filename = argument[i+1];
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }

        dt /= 1.57350; /// convert to ps based on having energy in k_BT, and length in nm
    }
};
double inline mag(double a, double b, double c){
    return sqrt(a*a+b*b+c*c);
}
// Given a bond, updates the force on all atoms correspondingly
void UpdateBondForces(System& sys){
    Molecules& molecule = sys.molecules;
    // Loops over the (2 for water) bond constraints
    for (Bond& bond : molecule.bonds){
        auto& atom1s = molecule.atoms[bond.a1];
        auto& atom2s = molecule.atoms[bond.a2];
        for (int i = 0; i<molecule.no_mol; i++ ){
            double dpx = atom1s.px[i]-atom2s.px[i];
            double dpy = atom1s.py[i]-atom2s.py[i];
            double dpz = atom1s.pz[i]-atom2s.pz[i];
            double dpmag = mag(dpx,dpy,dpz);
            double fac = -bond.K*(1-bond.L0/dpmag);
            double fx   = fac*dpx;
            double fy   = fac*dpy;
            double fz   = fac*dpz;
            atom1s.fx[i] += fx;
            atom1s.fy[i] += fy;
            atom1s.fz[i] += fz;
            atom2s.fx[i] -= fx;
            atom2s.fy[i] -= fy;
            atom2s.fz[i] -= fz;
            accumulated_forces_bond += mag(fx,fy,fz);
        }
    }
}

// Iterates over all bonds in molecules (for water only 2: the left and right)
// And updates forces on atoms correpondingly
void UpdateAngleForces(System& sys){
    Molecules& molecule = sys.molecules;
    for (Angle& angle : molecule.angles){
        auto& atom1s = molecule.atoms[angle.a1];
        auto& atom2s = molecule.atoms[angle.a2];
        auto& atom3s = molecule.atoms[angle.a3];
        #pragma omp simd reduction(+:accumulated_forces_angle)
        for (int i = 0; i<molecule.no_mol; i++ ){
            //====  angle forces  (H--O---H bonds) U_angle = 0.5*k_a(phi-phi_0)^2
            //f_H1 =  K(phi-ph0)/|H1O|*Ta
            // f_H2 =  K(phi-ph0)/|H2O|*Tc
            // f_O = -f1 - f2
            // Ta = norm(H1O x (H1O x H2O))
            // Tc = norm(H2O x (H2O x H1O))
            //=============================================================
            double d21x = atom2s.px[i]-atom1s.px[i];     
            double d21y = atom2s.py[i]-atom1s.py[i];     
            double d21z = atom2s.pz[i]-atom1s.pz[i];     
            double d23x = atom2s.px[i]-atom3s.px[i];     
            double d23y = atom2s.py[i]-atom3s.py[i];     
            double d23z = atom2s.pz[i]-atom3s.pz[i];     

            // phi = d21 dot d23 / |d21| |d23|
            double norm_d21 = mag(d21x, d21y, d21z);
            double norm_d23 = mag(d23x, d23y, d23z);
            double dot = d21x*d23x + d21y*d23y + d21z*d23z;
            double phi = acos(dot / (norm_d21*norm_d23));

            double c21_23x = d21y*d23z-d21z*d23y;
            double c21_23y = d21z*d23x-d21x*d23z;
            double c21_23z = d21x*d23y-d21y*d23x;

            double Tax = d21y*c21_23z-d21z*c21_23y;
            double Tay = d21z*c21_23x-d21x*c21_23z;
            double Taz = d21x*c21_23y-d21y*c21_23x;
            double TaMag = mag(Tax, Tay, Taz);
            Tax /= TaMag;
            Tay /= TaMag;
            Taz /= TaMag;

            // d23 cross (d23 cross d21) = - d23 cross (d21 cross d23) = c21_23 cross d23
            double Tcx = c21_23y*d23z-c21_23z*d23y;
            double Tcy = c21_23z*d23x-c21_23x*d23z;
            double Tcz = c21_23x*d23y-c21_23y*d23x;
            double TcMag = mag(Tcx, Tcy, Tcz);
            Tcx /= TcMag;
            Tcy /= TcMag;
            Tcz /= TcMag;

            double f1fac = angle.K*(phi-angle.Phi0)/norm_d21;
            double f1x = Tax*f1fac;
            double f1y = Tay*f1fac;
            double f1z = Taz*f1fac;
            double f3fac = (angle.K*(phi-angle.Phi0)/norm_d23); 
            double f3x = Tcx*f3fac;
            double f3y = Tcy*f3fac;
            double f3z = Tcz*f3fac;

            atom1s.fx[i] += f1x;
            atom1s.fy[i] += f1y;
            atom1s.fz[i] += f1z;
            atom2s.fx[i] -= f1x+f3x;
            atom2s.fy[i] -= f1y+f3y;
            atom2s.fz[i] -= f1z+f3z;
            atom3s.fx[i] += f3x;
            atom3s.fy[i] += f3y;
            atom3s.fz[i] += f3z;
            accumulated_forces_angle += mag(f1x,f1y,f1z) + mag(f3x,f3y,f3z);
        }
    }
}

// Iterates over all atoms in both molecules
// And updates forces on atoms correpondingly
void UpdateNonBondedForces(System& sys){
    /* nonbonded forces: only a force between atoms in different molecules
       The total non-bonded forces come from Lennard Jones (LJ) and coulomb interactions
       U = ep[(sigma/r)^12-(sigma/r)^6] + C*q1*q2/r */
    std::vector<double> fxs(sys.molecules.no_mol);
    std::vector<double> fys(sys.molecules.no_mol);
    std::vector<double> fzs(sys.molecules.no_mol);
    for (size_t typei = 0; typei<sys.molecules.atoms.size(); typei++){
        auto& atom1s = sys.molecules.atoms[typei];
        for (size_t typej = 0; typej<sys.molecules.atoms.size(); typej++){
            auto& atom2s = sys.molecules.atoms[typej];
            double const ep = sqrt(atom1s.ep*atom2s.ep); // ep = sqrt(ep1*ep2)
            double const sigma = 0.5*(atom1s.sigma+atom2s.sigma);  // sigma = (sigma1+sigma2)/2
            double const sigma2 = sigma*sigma;  // sigma = (sigma1+sigma2)/2
            double const q1 = atom1s.charge;
            double const q2 = atom2s.charge;
            for (int i = 0;   i < sys.molecules.no_mol; i++){
                #pragma omp simd
                for (int j = i+1;   j < sys.molecules.no_mol; j++){
                    double dpx = atom1s.px[i]-atom2s.px[j];
                    double dpy = atom1s.py[i]-atom2s.py[j];
                    double dpz = atom1s.pz[i]-atom2s.pz[j];

                    double r  = mag(dpx,dpy,dpz);       
                    double r2 = r*r;
                    double KC = 80*0.7;          // Coulomb prefactor
                    double sir = sigma2/r2; // crossection**2 times inverse squared distance
                    double pow2 = sir*sir;
                    double pow3 = pow2*sir;
                    double pow6 = pow3*pow3;
                    double fac1 = ep*6*(2*pow6-pow3)*sir;
                    double fac2 = KC*q1*q2/(r2*r);
                    double fx = fac1*dpx + fac2*dpx; // LJ + Coulomb forces
                    double fy = fac1*dpy + fac2*dpy; 
                    double fz = fac1*dpz + fac2*dpz; 
                    atom2s.fx[j] -= fx;
                    atom2s.fy[j] -= fy;
                    atom2s.fz[j] -= fz;
                    fxs[j] = fx;
                    fys[j] = fy;
                    fzs[j] = fz;
                }
                double x = 0;
                double y = 0;
                double z = 0;
                #pragma omp simd reduction(+:x,y,z,accumulated_forces_non_bond)
                for (int j = i+1;   j < sys.molecules.no_mol; j++){
                    x += fxs[j];
                    y += fys[j];
                    z += fzs[j];
                    accumulated_forces_non_bond += mag(fxs[j],fys[j],fzs[j]);
                }
                atom1s.fx[i] += x;
                atom1s.fy[i] += y;
                atom1s.fz[i] += z;
            }
        }
    }
}

// integrating the system for one time step using Leapfrog symplectic integration
void Evolve(System &sys, Sim_Configuration &sc){

    // Kick velocities and zero forces for next update
    // Drift positions: Loop over molecules and atoms inside the molecules
    for (auto& atomis : sys.molecules.atoms)
    for (int j = 0; j < sys.molecules.no_mol; j++){
        double fac = sc.dt/atomis.mass;
        atomis.vx[j] += fac*atomis.fx[j];    // Update the velocities
        atomis.vy[j] += fac*atomis.fy[j];
        atomis.vz[j] += fac*atomis.fz[j];
        atomis.fx[j]  = 0;// set the forces zero to prepare for next potential calculation
        atomis.fy[j]  = 0;
        atomis.fz[j]  = 0;
        atomis.px[j] += sc.dt* atomis.vx[j];               // update position
        atomis.py[j] += sc.dt* atomis.vy[j];
        atomis.pz[j] += sc.dt* atomis.vz[j];
    }

    // Update the forces on each particle based on the particles positions
    // Calculate the intermolecular forces in all molecules
    UpdateBondForces(sys);
    UpdateAngleForces(sys);
    // Calculate the intramolecular LJ and Coulomb potential forces between all molecules
    UpdateNonBondedForces(sys);

    sys.time += sc.dt; // update time
}

// Setup one water molecule
System MakeWater(int N_molecules){
    //===========================================================
    // creating water molecules at position X0,Y0,Z0. 3 atoms
    //                        H---O---H
    // The angle is 104.45 degrees and bond length is 0.09584 nm
    //===========================================================
    // mass units of dalton
    // initial velocity and force is set to zero for all the atoms by the constructor
    const double L0 = 0.09584;
    const double angle = 104.45*deg2rad;    

    System sys;
    // bonds beetween first H-O and second H-O respectively
    sys.molecules.bonds  = {{ .K = 20000, .L0 = L0, .a1 = 0, .a2 = 1},
                            { .K = 20000, .L0 = L0, .a1 = 0, .a2 = 2}};
    sys.molecules.angles = {{ .K = 1000, .Phi0 = angle, .a1 = 1, .a2 = 0, .a3 = 2 }};
    sys.molecules.atoms  = {Atoms(16, 0.65,    0.31, -0.82, "O", N_molecules), 
                            Atoms( 1, 0.18828, 0.238, 0.41, "H", N_molecules), 
                            Atoms( 1, 0.18828, 0.238, 0.41, "H", N_molecules)};
    for (int i = 0; i < N_molecules; i++){
        double P0x = i * 0.2;
        double P0y = i * 0.2;
        double P0z = 0;
        sys.molecules.atoms[0].px[i]=P0x;
        sys.molecules.atoms[0].py[i]=P0y;
        sys.molecules.atoms[0].pz[i]=P0z;
        sys.molecules.atoms[1].px[i]=P0x+L0*sin(angle/2);
        sys.molecules.atoms[1].py[i]=P0y+L0*cos(angle/2);
        sys.molecules.atoms[1].pz[i]=P0z;
        sys.molecules.atoms[2].px[i]=P0x-L0*sin(angle/2);
        sys.molecules.atoms[2].py[i]=P0y+L0*cos(angle/2);
        sys.molecules.atoms[2].pz[i]=P0z;
    }
    // Store atoms, bonds and angles in Water class and return
    sys.molecules.no_mol = N_molecules;
    return sys;
}

// Write the system configurations in the trajectory file.
void WriteOutput(System& sys, std::ofstream& file){  
    // Loop over all atoms in model one molecule at a time and write out position
    Molecules& molecule = sys.molecules;
    for (int i = 0; i<molecule.no_mol; i++)
    for (size_t j = 0; j < molecule.atoms.size(); j++){
        file << sys.time << " " << molecule.atoms[j].name << " " 
            << molecule.atoms[j].px[i] << " " 
            << molecule.atoms[j].py[i] << " " 
            << molecule.atoms[j].pz[i] << '\n';
    }
}

//======================================================================================================
//======================== Main function ===============================================================
//======================================================================================================
int main(int argc, char* argv[]){    
    Sim_Configuration sc({argv, argv+argc}); // Load the system configuration from command line data
    
    System sys  = MakeWater(sc.no_mol);   // this will create a system containing sc.no_mol water molecules
    std::ofstream file(sc.filename); // open file

    WriteOutput(sys, file);    // writing the initial configuration in the trajectory file
    
    auto tstart = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)
    
    // Molecular dynamics simulation
    for (int step = 0;step<sc.steps ; step++){

        Evolve(sys, sc); // evolving the system by one step
        if (step % sc.data_period == 0){
            //writing the configuration in the trajectory file
            WriteOutput(sys, file);
        }
    }

    auto tend = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)

    std::cout <<  "Elapsed time:" << std::setw(9) << std::setprecision(4)
              << (tend - tstart).count()*1e-9 << "\n";
    std::cout <<  "Accumulated forces Bonds   : "  << std::setw(9) << std::setprecision(5) 
              << accumulated_forces_bond << "\n";
    std::cout <<  "Accumulated forces Angles  : "  << std::setw(9) << std::setprecision(5)
              << accumulated_forces_angle << "\n";
    std::cout <<  "Accumulated forces Non-bond: "  << std::setw(9) << std::setprecision(5)
              << accumulated_forces_non_bond << "\n";
}
