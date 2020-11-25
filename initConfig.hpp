/*
    File: initConfig.hpp
    Function: Initialize chromosomes and nuclear envelope
    Model: chromoCell
    Created: 7 December, 2017 (XF)
*/

#ifndef INITCONFIG_HPP_INCLUDED
#define INITCONFIG_HPP_INCLUDED

#include <stdio.h>
//#include <math.h>
#include <iostream>
#include <fstream>

#include <cstring>
#include <string>
#include <sstream>
#include <iomanip>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cmath>

#include<random>
#include<chrono>

#include<omp.h>

#include <unistd.h>     // for getpid

using namespace std;

typedef struct {double x,y,z;} dVec;
typedef struct {int i,j,k;} iVec;
typedef struct {bool x,y,z;} bVec;    // index for Octree leaf
typedef struct {int i,j;} iPair;

// system settings
const double DT                 = 1E-4;     // s; 1.69E-5 s in Cheng2015 paper; 1E-3 used for RECONSTRUCTION
const double T                  = 1200+DT;  // 1200+DT //s; 1E5 used for RECONSTRUCTION
//const double T                  = 10+DT;       // s; used for Relaxation simulation only
const double T_INIT             = 0;        // s; temporarily freeze fluorophores
const double T_PREP             = 10;      // s; 120 relaxation stage
const int N_ITER_TOTAL          = T/DT;
const int CNT_START_INTERACT    = 0.0*N_ITER_TOTAL;
//const int FREQ_WRITE            = N_ITER_TOTAL/500;
//const int FREQ_PRINT            = N_ITER_TOTAL/1000;
const int FREQ_WRITE            = 1E4;   //1E4 default, 1E4*2 previously
const int FREQ_PRINT            = 5E3;
const int MSD_FREQ_WRITE        = 1E5;      // 3e3; 3E5 ~ every 30 seconds
const int MSD_WRITE_INTERVAL    = 200;      // 200   x DT = 20 ms time window
const int MSD_WRITE_DURATION    = 20000;    // 20000 x DT = 2  s

/////////////////////////////////////////////////////////////////////////               !!!!!!!!!!!!!!!!!!!!!!!!!
const int DC_FREQ_WRITE        = 3E5;      // every 30 seconds
const int DC_WRITE_INTERVAL    = 200;      // 200   x DT = 20 ms time window
const int DC_WRITE_DURATION    = 20000;    // 20000 x DT = 2  s

// ---------------------------------------- 01-08-2018 ---------------------------
const int DB_FREQ_WRITE        = FREQ_WRITE*1E-5;      // every 30 seconds; frequency for writing of dynamics Binders
// -------------------------------------------------------------------------------


const int DIM_VOXMAP            = 30;   // n x n x n
const int NUM_THREADS           = 4;
const bool OMP_ON               = true;
const pid_t PROC_ID             = getpid();           // process id of the current program

// constants
const double PI                 = 3.1415926535897;
const double AVOG               = 6.02214E+23;        // 1/mol, Avogadro number
const double NM_PER_BP          = 0.34;               // nm/bp; old value used in RECONSTRUCTION: 4.6875E-4 (where does it come from?)
const double RAD_NUCLEUS        = 1500;               // nm; Haering paper S. pombe
const double DEPTH_NUCLEOLUS    = 500;                // nm;
const double RAD_NUCLEOLUS      = 2*RAD_NUCLEUS - DEPTH_NUCLEOLUS;
const double K_BOLTZMANN        = 1.38E-2;                   // (pN*nm/K)
const double TEMPERATURE        = 310;                       // (K)  body temperature
const double COARSE_FACTOR      = 10;

/* -------- Parameters Brownian Simulation --------- */
const string SIM_TYPE           = "SIMULATION";

// structure
const int 	 INTERPHASE_TAD_SIZE      = 133; // in beads; ~ 267kb Yasu's paper Nature genetics 2017 
const int    MITOSIS_TAD_SIZE         = 240; // in beads; ~ 480kb Yasu's paper Nature genetics 2017
const int    BINDER_SPACING           = 10;             // every 10th bead is a binder

const double BP_SEP_NUCLEOSOME        = 154;                        // base pairs; 154  (constant) separation between nucleosomes in S. pombe (may include bps wrapping the nucleosome; read the paper Lantermann2010, Moyle-heyrmann2013 to double check!)
const double BP_SEP_CONDENSIN_BINDING = 15000;                      // base pairs; 15000 (default)  separation between condensin binding sites in S. pombe
//const int    NUM_BEAD = (int) 1876*200/BP_SEP_NUCLEOSOME + 1;      // ~2436.5 with ~ 2kb resolution = 1bead = 10 nucleosomes; NGS-8595.2000.I_ArmL.500.noself.normMatrix has 1876 rows (columns), representing 1876x2000 base pairs
//const int    NUM_BEAD = (int) 1876*2000/BP_SEP_NUCLEOSOME + 1;      // ~24365; NGS-8595.2000.I_ArmL.500.noself.normMatrix has 1876 rows (columns), representing 1876x2000 base pairs
const int    NUM_BEAD = 1880; //1880*2kb = 3.76 Mb  						// 2436 for testing dynamics
const double WHOLE_GENOME       = 12.0E+6;
const double VOL_FRAC_OF_GENOME = NUM_BEAD *COARSE_FACTOR * BP_SEP_NUCLEOSOME / WHOLE_GENOME;
const double MOLW_PER_BP_DNA    = 650;      // 650 g/mol or dalton
const double MOLW_PER_NUCLEO    = 2.1E+5;   // 2.1E+5 g/mol or dalton
const int    N_BP_PER_NUCLEO    = 146;      // 146 (constant) base pairs wrap around nucleosome

const int    NUM_TAD_INTERPHASE = NUM_BEAD/INTERPHASE_TAD_SIZE;
const int    NUM_TAD_MITOSIS    = NUM_BEAD/MITOSIS_TAD_SIZE;
const int SWITCH_INTERPHASE_MITOSIS   = 0;

//const double MASS_BEAD          = (MOLW_PER_BP_DNA+MOLW_PER_NUCLEO/N_BP_PER_NUCLEO)/AVOG*BP_SEP_NUCLEOSOME *1E-3;    // kg, ~
const double MASS_BEAD          = (MOLW_PER_NUCLEO + MOLW_PER_BP_DNA*BP_SEP_NUCLEOSOME)*COARSE_FACTOR/AVOG*1E-3;    // kg, ~3.6668E-22
//const double RAD_BEAD           = BP_SEP_NUCLEOSOME*NM_PER_BP;  // is this correct?
const double RAD_BEAD           = 5*5;        // nm; ~ 5nm 1 nucleosome, 10 nucleosomes = 50 nm or 15nm when better packing

// avoid floating point number error
const double L_DIFF_RESOLUTION  = 1E-3;        // nm

// dynamics
// ...viscous
const double GAMMA_VISCOUS      = 3E-8;       // kg/s;
// ...entropic
//const double F_ENTROPIC         = 1E-18;//1E-18;    // pN; 24.5 in Cheng2015 model
const double D_BROWNIAN_SQRT	= sqrt(2*K_BOLTZMANN*TEMPERATURE/GAMMA_VISCOUS*1E-3);	// nm^2/s sqrt(2*D)
const double F_ENTROPIC         = sqrt(2*GAMMA_VISCOUS*K_BOLTZMANN*TEMPERATURE*DT);  // ~1E-15
const double F_ENTROPIC_UNIT_MASS = F_ENTROPIC*1E-3 / MASS_BEAD;     //nN/kg
const double F_ENTROPIC_VARIANCE = 2*GAMMA_VISCOUS*K_BOLTZMANN*TEMPERATURE*1E3;  // pN*pN   *DT?
// ...tension [affects mean inter-bead gap size which influences chain passing]
const double K_TENSION          = 1E-3*COARSE_FACTOR*10;        // pN/nm; 50 in Cheng2015 model
const double K_TENSION_UNIT_MASS = K_TENSION*1E-3 / MASS_BEAD;       // nN/nm/kg
const double L0_TENSION         = 2*RAD_BEAD + 2;        // 12.5 nm; 15 in Cheng2015 model of budding yeast; use NM_PER_BP*BP_SEP_NUCLEOSOME ?
// ...repulsion
const double F_REPULSION        = 5E-1; //5E-1 previously but 5 appears to be better       //5E-4 default pN; 10 in Cheng2015 model
const double F_REPULSION_UNIT_MASS = F_REPULSION*1E-3 / MASS_BEAD;     // nN/kg
const double LTHR1_REPULSION    = 2*RAD_BEAD + 2;        // nm; 10 in Cheng2015 model
const double LTHR2_REPULSION    = 2*RAD_BEAD + 2;          // nm; 15 in Cheng2015 model
// ...attraction [affects mean angle between 3 beads]
const double K_ATTRACTION       = 1E-5*10;        // pN/nm (not given in Cheng2015 paper)
const double K_ATTRACTION_UNIT_MASS = K_ATTRACTION*1E-3 / MASS_BEAD;    // nN/nm/kg
// ...condensin tension (diffusion capture mechanism)
const double DIFFUSION_CAPTURE_ON   = true;
const double K_TENSION_CONDENSIN    = K_TENSION;    // pN/nm; 50 in Cheng2015 model
const double K_TENSION_CONDENSIN_UNIT_MASS = K_TENSION_CONDENSIN*1E-3 / MASS_BEAD;  // nN/nm/kg
const double LTHR_TENSION_CONDENSIN = (3*RAD_BEAD + 2);      // 40 default value; nm;
const double L0_TENSION_CONDENSIN   = (2*RAD_BEAD + 2);      // 30 default value; nm; 30 in Cheng2015 model
const double P_DISSOCIATE_CONDENSIN = 1E-2;    // a probability per simulation step that an existing interaction is lost (eventually force dependent?). 1E-3 for interphase; 1E-4 for condensin

const double P_DISSOCIATE_MITOSIS   = 1E-2; 
// ... binder tension (loop extrusion mechanism)
const double LOOP_EXTRUSION_ON      = false;
const double K_TENSION_BINDER       = K_TENSION;    // pN/nm
const double L0_TENSION_BINDER      = (2*RAD_BEAD + 2);    // 12.5 default value; nm   (check literature)
const double K_TENSION_BINDER_BEAD  = K_TENSION;
const double L0_TENSION_BINDER_BEAD = (2*RAD_BEAD + 2);    // 12.5 default value
const double P_BINDER_SLIDING       = 1E-4;    // 2E-5 default value; 1E-6 (too slow!) ~ 1E-4 (too fast?) ; in LE only, translocation rate: 30 bp/s (1E-5) 65~70 bp/s (2.5E-5) 75 bp/s (5E-5)
const double P_BINDER_STOP_UPON_CAPTURE = 0.9;
const double STATIC_BINDER_BLOCK_ON = false;
const string STATIC_BINDER_BLOCK_TYPE = "Interphase-TAD-Mitosis";  // "chipseq", "random", "TAD-Interphase", "TAD-Mitosis", "domains", "Interphase-domains-Mitosis", "Interphase-TAD-Mitosis"

// ... dynamic loading and unloading  1st
const bool SEPARATE_DYNAMICS_BINDER = false; // true 
const double P_ACTIVATION_BINDER      = 1E-6;;//1E-3;  // !!! need estimation !!  ~ a probability for loading of a binder
const double P_DISSOCIATE_BINDER      = 1E-3;//1E-3;  // !!! need estimation !!  ~ a probability for unloading of a binder

// ... dynamic loading and unloading  2nd
const bool COUPLED_DYNAMICS_BINDER  = false;  // ~ Mirny's way of doing it
const double P_COUPLED_BINDER       = 1E-2; //

// ... regulate concentration of Binders
const string DELETION_BINDER_TYPE     = "random"; // "uniform", "random" 
const double DELETE_BINDER            = 1.0;  // *100 = per cent of Binders to be deleted
//const int    Nth_REMOVE               = NUM_BEAD    ;
/* ----------------------08-08-2018-------------------------- */
const string BINDER_POSITION          =   "CHIP-corrected"; //"uniform"; // "random1"; "uniform"; "CHIP-seq";
const int NUM_CONDENSIN_PEAKS         = 154; // number of condensins to be loaded based on number of peaks in ChIP-seq data
// --------------------------------------------------------------

const int T_RESIDENCE                 = 500; // every 1st iteration do load/unload of Binders
// ---------------------12-08-2018-------------------------------

const string DIFFUSION_CAPTURE   = "multivalent"; //"probabilistic", "bivalent", "multivalent" - NO probab for bivalent and multivalent
const string DC_MECHANISM        = "BF-BF"; //"CSB-CSB", "CSB-BF"
const string DIFFUSION_CAPTURE_M = "multivalent";
const double P_DC_PROBAB         = 0.5; // needs an estimation!
const int    DC_VALENCE          = 8;

const string LOOP_EXTRUSION      = "symmetric"; //"asymmetricF"; "asymmetricR" "symmetric"; "probabilistic"; "asymprobabilistic" - NO probab for symmetric and asymmetric
const double P_LE_PROBAB         = 0.5; // needs an estimation!
//---------------------------------------------------------------


//const int N_BINDERS_ON_CHAIN     = int( double ((BP_SEP_NUCLEOSOME*10*NUM_BEAD/BP_SEP_CONDENSIN_BINDING)*(1-DELETE_BINDER)));


//----------------------------------------------------------------


/* ------------ Parameters Reconstruct ------------- */
// const string SIM_TYPE           = "RECONSTRUCT";

// structure
const double GENOMIC_RESOLUTION = 2000;     // bp (base pair ) "RECONSTRUCT" 3000 for budding data; 2000 for fission data;
/*
const double RAD_BEAD           = GENOMIC_RESOLUTION*NM_PER_BP;
const double MOLW_PER_BP_DNA    = 650;      // g/mol or dalton
const double MOLW_PER_NUCLEO    = 2.1E+5;   // g/mol or dalton
const int    N_NUCLEO_PER_BP    = 146;      // 146 base pairs wrap around nucleosome
const double MASS_BEAD          = (MOLW_PER_BP_DNA+MOLW_PER_NUCLEO/N_NUCLEO_PER_BP)/AVOG*GENOMIC_RESOLUTION *1E-3;    // kg, ~1E-23 (this is pure DNA. how to account for nucleosome?)
*/
// dynamics
const double CUTOFF_NORMMAT     = 0.001;      // "inter_NGS-9316_3kb_500" dataset has a max of 0.144361; 173 over 0.06; 675 over 0.05; 59823 over 0.005
const int    GENOMIC_LEAST_SEP_INTRA  = 0;  // make it 0 when nearest interactions are accounted for according to HiC map
const double K_NEAR             = 5E-21;    // nN/nm
const double K_INTRA            = 0.1*K_NEAR;
const double K_INTER            = 0.1*K_NEAR;

// data input
const string SUFFIX_normMatrix = ".normMatrix";
const string SUFFIX_numberBead = ".numberBead";
const string SUFFIX_chipseq    = ".csv";
const string SUFFIX_domains = ".csv";
const string SUFFIX_domains1 = ".Icsv";
const string SUFFIX_domains2 = ".Mcsv";
const string SUFFIX_PDB        = ".txt"; // ".pdb"
const string SUFFIX_random     = ".distrib" ; 

/* -------------------------------------------------- */

// [NOT USED] octree
const int OCTREE_DEPTH_MAX      = 4;        // at most 4096 smallest boxes
const int OCTREE_BRANCH_TRIGGER = 10;       // critical number of beads to trigger subdivision (branching)
const vector<bVec> OCTREE_POSTIIONS_CHILD   =   {{false,false,false}, {false,false,true}, {false,true,false}, {false,true,true},
                                                 {true,false,false},  {true,false,true},  {true,true,false},  {true,true,true}  };

// ===================== Classes =========================

/* -------------- Node ---------------- */
class Node
{
    int id, hostChromoId;
    dVec coord; // nm
    dVec veloc; // nm/s
    dVec accel, accelPrev; // nm/s/s

public:
    Node ();
    Node (int id, int hostChromoId);
    Node (int id, int hostChromoId, dVec coord, dVec veloc, dVec accel, dVec accelPrev);
    ~Node ();

    // setters
    void set_id(int id) {this->id=id;}
    void set_hostChromoId(int hostChromoId) {this->hostChromoId=hostChromoId;}
    void set_coord(dVec coord) {this->coord=coord;}
    void set_veloc(dVec veloc) {this->veloc=veloc;}
    void set_accel(dVec accel) {this->accel=accel;}
    void set_accelPrev(dVec accelPrev) {this->accelPrev=accelPrev;}

    // getters
    int get_id() const {return this->id;}
    int get_hostChromoId() const {return this->hostChromoId;}
    dVec get_coord() const {return this->coord;}
    dVec get_veloc() const {return this->veloc;}
    dVec get_accel() const {return this->accel;}
    dVec get_accelPrev() const {return this->accelPrev;}
};

/* -------------- Chromosomes ---------------- */
class Chromosome
{
    int id;
    vector<int> vecBeadIds; // may change to vector of pointer to Node objects in vecNode
    vector<int> vecGeneIds; // may change to vector of iPairs {start, end}
    vector<int> vecBinderIds;
    //iPair lociTelomere;
    //iPair lociCentromere;
public:
    Chromosome ();
    Chromosome (int id, vector<int> vecBeadIds);
    Chromosome (int id, vector<int> vecBeadIds, vector<int> vecGeneIds, vector<int> vecBinderIds);
    ~Chromosome();

    // setters
    void set_id(int id) {this->id = id;}
    void set_vecBeadIds(vector<int> vecBeadIds) {this->vecBeadIds = vecBeadIds;}
    void set_vecGeneIds(vector<int> vecGeneIds) {this->vecGeneIds = vecGeneIds;}
    void set_vecBinderIds(vector<int> vecBinderIds) {this->vecBinderIds = vecBinderIds;}

    // getters
    int get_id() const {return this->id;}
    vector<int> get_vecBeadIds() const {return this->vecBeadIds;}
    vector<int> get_vecGeneIds() const {return this->vecGeneIds;}
    vector<int> get_vecBinderIds() const {return this->vecBinderIds;}
};

/* -------------- Binders ---------------- */
class Binder
{
    int id, attachChromoId;
    vector<bool> canSlide;
    vector<int> vecBeadIds;         // should include just 2 beads (front and rear)
    vector<int> vecAttachBeadIds;   // should include just 2 beads (front attachment and rear attachment)

public:
    Binder ();
    Binder (int id, int attachChromoId, vector<bool> canSlide, vector<int> vecBeadIds, vector<int> vecAttachBeadIds);
    ~Binder ();

    // setters
    void set_id(int id) {this->id = id;}
    void set_attachChromoId(int attachChromoId) {this->attachChromoId = attachChromoId;}
    void set_canSlide(vector<bool> canSlide) {this->canSlide = canSlide;}
    void set_vecBeadIds(vector<int> vecBeadIds) {this->vecBeadIds = vecBeadIds;}
    void set_vecAttachBeadIds(vector<int> vecAttachBeadIds) {this->vecAttachBeadIds = vecAttachBeadIds;}

    // getters
    int get_id() const {return this->id;}
    int get_attachChromoId() const {return this->attachChromoId;}
    vector<bool> get_canSlide() const {return this->canSlide;}
    vector<int> get_vecBeadIds() const {return this->vecBeadIds;}
    vector<int> get_vecAttachBeadIds() const {return this->vecAttachBeadIds;}
};

/* -------------- Domains ---------------- */
class Domain
{
    string name;        // e.g., "telomere_left", "telomere_right", "centromere", "gene_XXXXXX"
    int hostChromoId;
    iPair location;     // {start, end}; if start > end, it's on the reverse strand

public:
    Domain ();
    Domain (string name, int hostChromoId, iPair location);
    ~Domain ();

    // setters
    void set_name(string name) {this->name = name;}
    void set_hostChromoId(int hostChromoId) {this->hostChromoId = hostChromoId;}
    void set_location(iPair location) {this->location = location;}

    // getters
    string get_name() const {return this->name;}
    int get_hostChromoId() const {return this->hostChromoId;}
    iPair get_location() const {return this->location;}
};

/* -------------- VoxMap ----------------- */
class VoxMap
{
    iVec dim;
    unordered_map<int, vector<int>> mapBeadIds;         // contains Bead ids in VoxBox; key is collapsed from 3d indexing
    unordered_map<int, vector<int>> mapBeadIdsHalo;     // contains Bead ids in halo of VoxBox;

public:
    VoxMap ();
    VoxMap (iVec dim, unordered_map<int, vector<int>> mapBeadIds, unordered_map<int, vector<int>> mapBeadIdsHalo);
    ~VoxMap ();

    // setters
    void set_dim(iVec dim) {this->dim = dim;}
    void set_mapBeadIds(unordered_map<int, vector<int>> mapBeadIds) {this->mapBeadIds = mapBeadIds;}
    void set_mapBeadIdsHalo(unordered_map<int, vector<int>> mapBeadIdsHalo) {this->mapBeadIdsHalo = mapBeadIdsHalo;}

    // getters
    iVec get_dim() const {return this->dim;}
    unordered_map<int, vector<int>> get_mapBeadIds() const {return this->mapBeadIds;}
    unordered_map<int, vector<int>> get_mapBeadIdsHalo() const {return this->mapBeadIdsHalo;}
};

/* -------------- Octree (NOT FINISHED) ---------------- */
class Leaf
{
    vector<bVec> position;   // size indicates depth, {0, 1} indicate location at that depth
    vector<int> vecBeadIds;  // vector of bead ids located at this leaf
    vector<Leaf> vecLeaf;    // vector of Leafs located at this leaf
    bool flagBranched;       // indicates if this is an "end" leaf

public:
    Leaf ();
    Leaf (vector<bVec> position, vector<int> vecBeadIds, vector<Leaf> vecLeaf, bool flagBranched);
    ~Leaf ();

    // setters
    void set_position(vector<bVec> position) {this->position = position;}
    void set_vecBeadIds(vector<int> vecBeadIds) {this->vecBeadIds = vecBeadIds;}
    void set_vecLeaf(vector<Leaf> vecLeaf) {this->vecLeaf = vecLeaf;}
    void set_flagBranched(bool flagBranched) {this->flagBranched = flagBranched;}

    // getters
    vector<bVec> get_position() const {return this->position;}
    vector<int> get_vecBeadIds() const {return this->vecBeadIds;}
    vector<Leaf> get_vecLeaf() const {return this->vecLeaf;}
    bool get_flagBranched() const {return this->flagBranched;}
};
class Octree
{
    int depth;
    vector<Leaf> vecLeaf;

public:
    Octree ();
    Octree (int depth, vector<Leaf>);
    ~Octree ();

    // setters
    void set_vecLeaf(vector<Leaf> vecLeaf) {this->vecLeaf = vecLeaf;}
    void set_depth(int depth) {this->depth;}

    // getters
    int get_depth() const {return this->depth;}
    vector<Leaf> get_vecLeaf() const {return this->vecLeaf;}
};

// ===================== Functions =========================
// --- following functions are for RECONSTRUCTION objective ---
void initChromosome(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<iPair> & vecBeadPairInteract, \
                    vector<double> & vecBeadPairInteractFreq, int argc, char ** argv);
void initDynamics(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, int argc, char ** argv, string initType, vector<double> & randNum01);
// ------------------------------------------------------------
// --- following functions are for Brownian SIMULATION objective ---
void initChromosomeBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode);
void initDynamicsBrownianSimulation(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecFluoroSite, string initType, vector<double> & randNum01, int argc, char ** argv);
void initBindersBrownianSimulation(vector<Binder> & vecBinder, vector<Node> & vecNode, vector<int> & vecCondensinSite, vector<Chromosome> & vecChromosome);
void initBindersLoadingBrownianSimulation(vector<Binder> & vecBinder, vector<Node> & vecNode, vector<int> & vecCondensinSite, vector<Chromosome> & vecChromosome,  \
                                          vector<int> & vecLoadedBinders, vector<int> & vecCondensinSiteEmpty, vector<int> & vecBinderIds2);   // ------------- 21-06-2018 --------------- k_on

// ------------------------ 14-08-2018 -----------------------------------------------------------
void initBindersUnloadingBrownianSimulation(vector<Binder> & vecBinder, vector<Node> & vecNode, vector<int> & vecCondensinSite, vector<Chromosome> & vecChromosome, vector<int> & vecLoadedBinders);



// -----------------------------------------------------------------
void printChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<Binder> & vecBinder);
void writeChromoInfo(vector<Chromosome> & vecChromosome, vector<Node> & vecNode, vector<int> & vecFluoroSite, \
                     vector<int> & vecCondensinSite, vector<int> & vecBinderBlockSite, vector<Binder> & vecBinder, \
                     VoxMap & voxMap, double t_now);

void writeMSD(vector<Node> & vecNode, double t_now);

/////////////////////////////////////////////////////////////////////////

//void writeDC(vector<iPair> & capturer, double t_now);
void writeDC2(vector<iPair> & capturer, double t_now);

// -------------------------------- 03 - 08 - 2018 -----------------------------------------------
void writeDC2_stats(vector<iPair> & capturer, vector<int> & hit_valence, double t_now);


/////////////////////////////////////////////////////////////////////////
void writeDynamicsBinders(vector<Binder> & vecBinder, vector<int> & vecLoadedBinders, vector<int> & vecUnloadedBinders, vector<int> & vecCondensinSiteEmpty, double t_now, vector<int> & vecUnLoadedBinders);
/////////////////////////////////////////////////////////////////////////

#endif // INITCONFIG_HPP_INCLUDED
