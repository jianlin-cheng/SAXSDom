////  Define header for Foxs, be sure to put this part first, otherwise, get error
#include <IMP.h>
#include <IMP/foxs/internal/Gnuplot.h>
#include <IMP/foxs/internal/JmolWriter.h>

#include <IMP/saxs/Profile.h>
#include <IMP/saxs/ProfileFitter.h>
#include <IMP/saxs/ChiScoreLog.h>
#include <IMP/saxs/ChiFreeScore.h>
#include <IMP/saxs/RatioVolatilityScore.h>
#include <IMP/saxs/FormFactorTable.h>
#include <IMP/saxs/utility.h>
#include <IMP/saxs/utility_SAXSDom.h>

#include <IMP/benchmark/Profiler.h>
#include <IMP/saxs/Distribution.h>
#include <IMP/saxs/RadiusOfGyrationRestraint.h>
#include <boost/program_options.hpp>
#include "boost/filesystem.hpp" 
namespace po = boost::program_options;

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

struct stat st = {0};

using namespace IMP::saxs;
using namespace IMP::foxs::internal;

//// Define header for Unicon3D
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <cassert>
#include <functional>
#include <string>
#include <vector>
#include <cfloat>
#include <cstring>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <unistd.h>
#include "mocapy.h"
#include "vonmises2d/vonmises2dess.h"
#include "vonmises2d/vonmises2ddensities.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "MatVec.h"
#include <ctime> 
#define X 0
#define Y 1
#define Z 2

// import C++ Standard Library, you can use 'vector' instead of using 'std::vector'

//// Define headers for pulchra
#include <sys/timeb.h>
#define COMPILE_BB
#define COMPILE_ROT
#define uchar unsigned char
#define uint unsigned int
#define real double

#include "pulchra_306/pulchra_common.h"
#include "pulchra_306/nco_data_new.h"
#include "pulchra_306/rot_data_coords.h"
#include "pulchra_306/rot_data_idx.h"


#define PULCHRA_VERSION 3.06
#define MAX_BUF_SIZE 1000

#define FILE_SUCCESS     0
#define FILE_NOT_FOUND  -1
#define FILE_WARNING    -2

#define FATAL_MAE -1

#define FLAG_BACKBONE  1
#define FLAG_CALPHA    2
#define FLAG_SIDECHAIN 4
#define FLAG_SCM       8
#define FLAG_INITIAL  16

#define FLAG_PROTEIN  1
#define FLAG_DNA      2
#define FLAG_RNA      4
#define FLAG_CHYDRO   8

#define RADDEG 180./M_PI
#define DEGRAD M_PI/180.

int _VERBOSE = 0;
int _BB_REARRANGE = 1;
int _BB_OPTIMIZE = 0;
int _CA_OPTIMIZE = 1;
int _CA_RANDOM = 0;
int _CA_ITER = 100;
int _CA_TRAJECTORY = 0;
int _CISPRO = 0;
int _CHIRAL = 1;
int _CENTER_CHAIN = 0;
int _REBUILD_BB = 1;
int _REBUILD_SC = 1;
int _REBUILD_H = 0;
int _PDB_SG = 0;
int _TIME_SEED = 0;
int _XVOLUME = 1;
int _XVOL_ITER = 3;
int _PRESERVE = 1;
real _CA_START_DIST = 3.0;
real _CA_XVOL_DIST = 3.5;
real _SG_XVOL_DIST = 1.6;


#define CALC_C_ALPHA
#define CALC_C_ALPHA_ANGLES
#define CALC_C_ALPHA_START
#define CALC_C_ALPHA_XVOL

real CA_K=10.0;
real CA_ANGLE_K=20.0;
real CA_START_K=0.01;
real CA_XVOL_K=10.00;

#define CA_DIST 3.8
#define CA_DIST_TOL 0.1
#define CA_DIST_CISPRO 2.9
#define CA_DIST_CISPRO_TOL 0.1
#define E_EPS 1e-10

#ifndef bool
#define bool int
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

int **RBINS = NULL;
real **X_COORDS = NULL;
real **C_ALPHA = NULL;

// grid resolution (used for fast clash detection)
#define GRID_RES 6.0

int chain_length = 0;

char AA_NAMES[21][4] =
  { "GLY", "ALA", "SER", "CYS", "VAL",
    "THR", "ILE", "PRO", "MET", "ASP",
    "ASN", "LEU", "LYS", "GLU", "GLN",
    "ARG", "HIS", "PHE", "TYR", "TRP",
    "UNK" };

char SHORT_AA_NAMES[22] = { "GASCVTIPMDNLKEQRHFYWX" };

int AA_NUMS[256];

int nheavy[20] = { 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 4, 5, 5, 5, 7, 6, 7, 8, 10};

char *backbone_atoms[4] = { "N  ", "CA ", "C  ", "O  " };

char *heavy_atoms[1000]= {
/* GLY */  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ALA */ "CB ", NULL,   NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* SER */ "CB ", "OG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* CYS */ "CB ", "SG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* VAL */ "CB ", "CG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* THR */ "CB ", "OG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ILE */ "CB ", "CG1", "CG2", "CD1",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* PRO */ "CB ", "CG ", "CD ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* MET */ "CB ", "CG ", "SD ", "CE ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASP */ "CB ", "CG ", "OD1", "OD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASN */ "CB ", "CG ", "OD1", "ND2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LEU */ "CB ", "CG ", "CD1", "CD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LYS */ "CB ", "CG ", "CD ", "CE ", "NZ ",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLU */ "CB ", "CG ", "CD ", "OE1", "OE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLN */ "CB ", "CG ", "CD ", "OE1", "NE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* ARG */ "CB ", "CG ", "CD ", "NE ", "CZ ", "NH1", "NH2",  NULL,  NULL,  NULL,
/* HIS */ "CB ", "CG ", "ND1", "CD2", "CE1", "NE2",  NULL,  NULL,  NULL,  NULL,
/* PHE */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ",  NULL,  NULL,  NULL,
/* TYR */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ", "OH ",  NULL,  NULL,
/* TRP */ "CB ", "CG ", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};

/* reads full-atom pdb file */

struct _res_type;

typedef struct _atom_type {
  struct _atom_type *next;
  real x, y, z;
  char *name;
  int num, locnum;
  int flag;
  char cispro;
  int gx, gy, gz;
  struct _res_type *res;
  struct _atom_type *prev;
} atom_type;

typedef struct _res_type {
  struct _res_type *next;
  atom_type *atoms;
  int num, locnum, natoms;
  int type;
  char pdbsg;
  char protein;
  char *name;
  char chain;
  real sgx, sgy, sgz;
  real cmx, cmy, cmz;
  struct _res_type *prev;
} res_type;

typedef struct _mol_type {
  struct _mol_type *next;
  res_type *residua;
  int nres;
  unsigned char *r14;
  char *name;
  uchar *seq;
  char **contacts;
  real **cutoffs;
  struct _mol_type *prev;
} mol_type;

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

mol_type *chain = NULL;


real rnd(void)
{
  return 0.001*(real)(rand()%1000);
}

// atom/res/mol manipulation functions

atom_type *new_atom(void)
{
  atom_type *tmpatom;

    tmpatom = (atom_type*) calloc(sizeof(atom_type),1);
    if (tmpatom) {
      tmpatom->x=tmpatom->y=tmpatom->z=0.;
      tmpatom->name=NULL;
      tmpatom->num=tmpatom->locnum=tmpatom->flag=0;
      tmpatom->next=tmpatom->prev=NULL;
    }

  return tmpatom;
}

res_type* new_res(void)
{
  res_type *tmpres;

    tmpres = (res_type*) calloc(sizeof(res_type),1);
    if (tmpres) {
      tmpres->num=0;
      tmpres->name=NULL;
      tmpres->atoms=NULL;
      tmpres->chain=' ';
      tmpres->next=tmpres->prev=NULL;
    }

  return tmpres;
}

mol_type *new_mol(void)
{
  mol_type *tmpmol;

    tmpmol = (mol_type*) calloc(sizeof(mol_type),1);
    if (tmpmol) {
      tmpmol->name=NULL;
      tmpmol->residua=NULL;
      tmpmol->next=tmpmol->prev=NULL;
    }

  return tmpmol;
}


typedef struct _atom_list {
  atom_type *atom;
  struct _atom_list *next;
} atom_list;

int read_pdb_file(char* filename, mol_type* molecules, char *realname);
int read_pdb_file_saxs(char* filename, mol_type* molecules, char *realname, std::string pdbString);
real calc_ca_energy(atom_type **c_alpha, real **new_c_alpha, real **init_c_alpha, real **gradient, real alpha, real *ene, bool calc_gradient);
void ca_optimize(char *tname, char *iname);
void center_chain(mol_type *mol);
void write_pdb(char *name, mol_type *mol);
void write_pdb_insideSAXS(char *name, std::string &pdbString_pulchar, mol_type *mol);
void write_pdb_sg(char *name, mol_type *mol);
// distance
real calc_distance(real x1, real y1, real z1,real x2, real y2, real z2);
real calc_r14(real x1, real y1, real z1, real x2, real y2, real z2, real x3, real y3, real z3, real x4, real y4, real z4);
real superimpose2(real **coords1, real **coords2, int npoints, real **tpoints, int ntpoints);
atom_type *find_atom(res_type *res, char *aname);
void add_replace(res_type *res, char *aname, real x, real y, real z, int flags);
void prepare_rbins(void);
void rebuild_backbone(void);
int read_rotamers(void);
void cross(real *v1, real *v2, real *v3);
void norm(real *v);
int check_xvol(res_type *res);
void rebuild_sidechains(void);
int get_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid);
int display_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid);
void allocate_grid(atom_list *****grid_, int *xgrid_, int *ygrid_, int *zgrid_);
void optimize_exvol(void);
void optimize_backbone(mol_type *chain);
int chirality_check(void);


//########################################### end define pulchar
using namespace mocapy;
using namespace std;
using std::vector;

const double PI = std::atan(1.0)*4;
const double RAD2DEG = 180/PI;
const double DEG2RAD = PI/180;
const double CA2CA = 3.8;
const double START_TEMP = 1000.0;
const double FINAL_TEMP = 298.0;
const double BOLTZMANN_CONSTANT = 0.0019872041;
//const double RR_CONTACT_THRESHOLD = 8.0;
const double RR_CONTACT_THRESHOLD = 25.0; // increate threshold for domain
const int RR_SEQ_SEP = 1;
const int MIN_FOLDON_LEN = 20;
const int NUM_DAT = 7;
const int NUM_MIS = 6;

// Define a base random number generator and initialize it with a seed.
boost::minstd_rand baseGen(std::time(0));

// Define distribution U[0,1) [double values]
boost::uniform_real<> uniDblUnit(0,1);

// Define a random variate generator using our base generator and distribution
boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uniDblGen(baseGen, uniDblUnit);

// amino acid residue types
const int ALA = 0; // Alanine	Ala	   A
const int CYS = 1; // Cysteine	Cys	   C
const int ASP = 2; // Aspartic Acid	Asp	   D
const int GLU = 3; // Glutamic Acid	Glu	   E
const int PHE = 4; // Phenylalanine	Phe	   F
const int GLY = 5; // Glycine	Gly	   G
const int HIS = 6; // Histidine	His	   H
const int ILE = 7; // Isoleucine	Ile	   I
const int LYS = 8; // Lysine	Lys	   K
const int LEU = 9; // Leucine	Leu	   L
const int MET = 10; // Methionine	Met	   M
const int ASN = 11; // Asparagine	Asn	   N
const int PRO = 12; // Proline	Pro	   P
const int GLN = 13; // Glutamine	Gln	   Q
const int ARG = 14; // Arginine	Arg	   R
const int SER = 15; // Serine	Ser	   S
const int THR = 16; // Threonine	Thr	   T
const int VAL = 17; // Valine	Val	   V
const int TRP = 18; // Tryptophan	Trp	   W
const int TYR = 19; // Tyrosine	Tyr	   Y


// amino acide residue code to three letter format
string seq3[20] = {"ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR"};

// amino acide residue code to one letter format
string seq[20] = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

// secondary structure types
const int ALPHA_HELIX = 0;
const int THREE_HELIX = 1;
const int FIVE_HELIX = 2;
const int ISOLATED_STRAND = 3;
const int EXTENDED_STRAND = 4;
const int HBONDED_TURN = 5;
const int NONHBONDED_BEND = 6;
const int RANDOM_COIL = 7;

// secondary structure code to one letter format
string sec[8] = {"H", "G", "I", "B", "E", "T", "S", "C"};

// point3d object
struct point3d {
    double x;
    double y;
    double z;
};

// pdbInfo object
struct pdbInfo {
    int id;
    int aa;
    point3d ca;
    point3d sc;
};

// poseInfo object
struct poseInfo {
    int id;
    int aa;
    int ss;
    double bca;
    double tao;
    double theta;
    double bsc;
    double phi;
    double delta;
};

// contactInfo object
struct contactInfo {
    int ri;
    int rj;
    double wt;
};

// vector to store the conatcts
vector<contactInfo> rrContact;

//max and min fragment size
int MIN_FRG_LEN = 1;
int MAX_FRG_LEN = 5; 
int adjustLinker = 1; 

// weights of energy terms
double w_sc_sc = 1.0;
double w_sc_bb = 2.73684;
double w_bb_bb = 0.06833;
double w_ri_rj = 0.0;


//double w_ri_rj = 0.0;

double saxs_chi_global_minima = 20.0; 
int total_sampling_cycle = 0; 

double w_saxs_chi_penalty = 10.0;   
double w_saxs_chi_penalty_initial = 0.0001; 
double w_saxs_chi_penalty_final =10;
double w_saxs_chi_final =  1; 
double w_saxs_chi_initial = 30.0; 
double w_saxs_chi = 30; 


double w_saxs_KL_penalty = 10.0;   
double w_saxs_KL_penalty_initial = 0.0001; 
double w_saxs_KL_penalty_final =10;
double w_saxs_KL_final =  1; 
double w_saxs_KL_initial = 700.0; 
double w_saxs_KL = 700; 


double w_saxs_score2_penalty = 10.0;   
double w_saxs_score2_penalty_initial = 0.0001; 
double w_saxs_score2_penalty_final =10;
double w_saxs_score2_final =  1; 
double w_saxs_score2_initial = 700.0; 
double w_saxs_score2 = 700; 



double w_saxs_RG_normalize_penalty = 10.0;   
double w_saxs_RG_normalize_penalty_initial = 0.0001; 
double w_saxs_RG_normalize_penalty_final =10;
double w_saxs_RG_normalize_final =  1; 
double w_saxs_RG_normalize_initial = 800.0; 
double w_saxs_RG_normalize = 800; 


int scoreFunction = 1;
bool scoreCombine = false;

// energyInfo object
struct energyInfo {
    double sc_sc;
    double sc_bb;
    double bb_bb;
    double ri_rj;
    double saxs_chi_global;
    double saxs_KL;
    double saxs_score2;
    double saxs_RG_normalize;
    double structure_energy;
    double saxs_energy;
    double saxs_penalty;
    double total;
};

// sseInfo object
struct sseInfo {
    int start;
    int end;
    int type;
};

// amino acid string defined by Liwo et. al. 1997 J. Comp. Chem. I
string liwoAminoStr = "CMFILVWYAGTSQNEDHRKP";

// epsilon0 parameter (lower triangle and digonal) for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoEpsilon0[20][20] = {
    {1.050, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.030, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020},
    {1.260, 1.450, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.030, 0.030, 0.030, 0.020, 0.020, 0.020, 0.030, 0.020},
    {1.190, 1.340, 1.270, 0.010, 0.010, 0.010, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.010, 0.010, 0.020, 0.010, 0.010, 0.010},
    {1.300, 1.470, 1.410, 1.580, 0.010, 0.010, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {1.250, 1.510, 1.400, 1.590, 1.550, 0.010, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {1.170, 1.380, 1.310, 1.520, 1.500, 1.400, 0.020, 0.010, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {0.990, 1.170, 1.150, 1.210, 1.180, 1.100, 0.970, 0.020, 0.020, 0.020, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020},
    {0.920, 1.150, 1.050, 1.220, 1.180, 1.040, 0.870, 0.810, 0.010, 0.010, 0.010, 0.010, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {0.980, 1.200, 0.990, 1.240, 1.260, 1.190, 0.770, 0.810, 1.020, 0.010, 0.010, 0.010, 0.030, 0.020, 0.020, 0.020, 0.020, 0.020, 0.020, 0.010},
    {0.980, 1.030, 0.840, 1.060, 1.130, 1.010, 0.710, 0.720, 0.820, 0.560, 0.010, 0.010, 0.020, 0.020, 0.010, 0.020, 0.020, 0.020, 0.000, 0.010},
    {0.800, 0.910, 0.760, 0.980, 0.900, 0.870, 0.560, 0.580, 0.640, 0.550, 0.430, 0.010, 0.020, 0.020, 0.020, 0.010, 0.020, 0.020, 0.000, 0.010},
    {0.770, 0.860, 0.680, 0.900, 0.930, 0.830, 0.510, 0.520, 0.590, 0.470, 0.450, 0.280, 0.020, 0.020, 0.010, 0.020, 0.020, 0.020, 0.000, 0.010},
    {0.810, 1.020, 0.720, 0.950, 0.980, 0.850, 0.590, 0.630, 0.750, 0.330, 0.350, 0.260, -0.280, 0.060, 0.030, 0.020, 0.030, 0.040, 0.010, 0.020},
    {0.730, 0.950, 0.700, 0.870, 1.000, 0.890, 0.600, 0.600, 0.770, 0.490, 0.380, 0.380, 0.530, 0.660, 0.010, 0.040, 0.030, 0.050, 0.000, 0.020},
    {0.640, 0.810, 0.530, 0.880, 0.790, 0.720, 0.520, 0.510, 0.470, -0.060, 0.200, 0.040, -0.230, -0.020, -1.580, 0.100, 0.030, 0.040, 0.040, 0.020},
    {0.680, 0.640, 0.520, 0.790, 0.680, 0.620, 0.530, 0.570, 0.510, 0.230, 0.290, 0.120, -0.120, 0.270, -0.930, -0.660, 0.030, 0.030, 0.050, 0.020},
    {0.910, 1.050, 0.920, 0.940, 0.980, 0.830, 0.820, 0.760, 0.650, 0.560, 0.570, 0.490, 0.380, 0.590, 0.420, 0.620, 0.800, 0.030, 0.000, 0.020},
    {0.580, 0.870, 0.690, 0.960, 0.950, 0.750, 0.660, 0.670, 0.530, 0.380, 0.430, 0.400, 0.360, 0.330, 1.010, 1.000, 0.490, -0.020, 0.090, 0.020},
    {0.590, 0.810, 0.550, 0.960, 0.970, 0.850, 0.500, 0.620, 0.680, -0.010, 0.000, -0.010, -0.020, 0.000, 1.300, 1.090, -0.010, -0.480, -11.960, 0.020},
    {0.820, 0.950, 0.810, 0.980, 1.000, 0.920, 0.770, 0.790, 0.740, 0.690, 0.570, 0.580, 0.620, 0.620, 0.420, 0.420, 0.610, 0.530, 0.560, 0.820}
};

// sigma_0 parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoSigma0[20] = {2.4366, 2.3204, 2.5098, 2.5933, 2.2823, 2.3359, 2.3409, 2.5919, 2.7249, 2.8905, 2.4984, 2.6954, 2.723, 2.6269, 2.3694, 2.4471, 2.6047, 2.7251, 1.6947, 2.1346};

// r_0 parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoR0[20] = {2.1574, 5.7866, 1.4061, 1.2819, 6.3367, 2.5197, 3.357, 4.4859, 0.2712, 3.3121, 3.5449, 0.7634, 3.332, 1.1813, 1.8119, 2.2432, 3.0723, 3.777, 9.2904, 4.8607};

// sigma_double_bar_over_sigma_tee_square parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoSigmaDoubleBarOverSigmaTeeSquare[20] = {1.809, 2.6006, 1.9262, 2.5707, 3.964, 1.0429, 3.6263, 3.2406, 8.0078, 2.3636, 4.4303, 2.0433, 1.7905, 2.6172, 6.6061, 1.6795, 2.2451, 2.0347, 7.5089, 5.9976};

// chi_prime parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoChiPrime[20] = {0.0333, -0.0025, 0.309, 0.4904, 0.0491, 0.2238, 0.1351, 0.0897, 0.579, 0.0749, 0.0968, 0.2732, -0.1105, 0.296, 0.2624, -0.0029, 0.0236, 0.077, 0.0731, 0.1177};

// alpha parameter for GBV Potential by Liwo et. al. 1997 J. Comp. Chem. I
double liwoAlpha[20] = {0.1052, 0.0299, 0.025, -0.0266, 0.0801, -0.1277, 0.0589, 0.0664, 0.0115, 0.1108, 0.0878, -0.0064, -0.019, 0.0505, 0.0062, -0.0348, -0.0264, 0.0679, 0.0549, 0.0438};

// bsc parameter for side chain length by Levitt M. 1976 JMB (Table 1)
double levittBsc[20] = {0.77, 1.38, 1.99, 2.63, 2.97, 0.0, 2.76, 1.83, 2.94, 2.08, 2.34, 1.98, 1.42, 2.58, 3.72, 1.28, 1.43, 1.49, 3.58, 3.36};

//global variables
char jobId[1000] = "";
char faFile[1000] = "";
char ssFile[1000] = "";
char alnFile[1000] = ""; 
char cmFile[1000] = "";
char dmFile[1000] = "";
char moFile[1000] = "";
int numCycles = 300;
int genInitial_numDecoys = 1;
int numDecoys = 20;
char natFile[1000] = "";
char simulationStatFile[1000] = "";
char statFile[1000] = "";
char statScoreFile[1000] = "";
char saxsFile[1000] = ""; 
char zeroFile[1000] = ""; 
char outputdir[1000] = ""; 
int saxsFit = 0; 
int saxsCAFit = 0; 
int saxsfoxsonly = 0; 
int saxspulcharonly = 0; 
vector<std::string> FullFastaArray;
vector<std::string> FullSSArray;

Profiles exp_profiles_custom; 
bool exp_profiles_loaded = false;



char pdbTempFile[1000];
char pdbTempFile_prefix[1000];
char pdbTempFile_pulchra[1000]; 
char pdbTempFile_pose[1000];
char bufposeStat[1000];
char pdbfile_initial[1000];

char pdbfile_initial_pulchra[1000];

char pdbSimulatedFile[1000];
char simpdbfile_initial_pulchra[1000];	


char pdbdomainFile[1000];
char dompdbfile_initial_pulchra[1000];	

char pdbFoldonFile[1000];
// buffer statistics
char bufStat[1000];
char pdbFoldonFile_pulchra[1000];
char bufsimStat[1000];
string pdbString;
string pdbString_pulchar;


#include <math.h>
#define LN_2      0.69314718055994530942      /* ln(2) */
#define log2(x)   (log(x)/LN_2)

// pr_dist object
struct pr_dist {
    double radius;
    double dist;
};
// domLinkerInfo object
struct domLinkerInfo {
    double start;
    double end;
};
// saxsInfo object
struct saxsInfo {
  string expfile;
  string pdfile;
  double chi_score;
  double expRG;
  double modelRG;
  double RGdiff;
  double RGdiff_normalize;
  double KL;
  double KL_symetric;
  double fun_score1;
  double fun_score2;
  double fun_score3;
  double fun_score4;
  double fun_score5;
  double fun_score6;
  double fun_score7;
  double fun_score8;
  double fun_score9;
  double fun_score10;
  double fun_score11;
			  
			  
};
// fitting_profile object
struct fitting_profile {
    double q;
    double expInt;
    double modelInt;
    double error;
};
// radius of gyration object
struct RG_saxs {
    double expRG;
    double modelRG;
    double RGdiff;
    double RGdiff_normalize;
};
// radius of gyration object
struct sample_region {
    int start;
    int end;
};


//indicators
bool jId = false;
bool fFile = false;
bool sFile = false;
bool aFile = false;
bool cFile = false;
bool dFile = false;
bool mFile = false;
bool nFile = false;
bool eFile = false; 
bool odir = false; 
bool zFile = false; 
bool tofitting = false; 
bool CAfitting = false; 
bool genInitial = false; 
bool regularize = false; 
bool coolingDown = false;
bool foxsonly = false;
bool pulcharonly = false;
bool evaluateonly = false;
std::string testonly = "";
int allowJump = 0;

void parseNextItem(int argc, char ** argv, int & i);
void parseNextItem_special(int argc, char ** argv, int & i);
void parseCommandLine(int argc, char ** argv);

int getAA(const char * aa);
int getSS(const char * ss);

double getDistance(point3d & p1, point3d &p2);
double getDotProduct(point3d & p1, point3d &p2);
point3d getCrossProduct(point3d & p1, point3d &p2);
double getNorm(point3d & p);
point3d getDifference(point3d & p1, point3d &p2);
point3d getMidpoint(point3d & p1, point3d &p2);
point3d getUnit(point3d & p);
double getAngle(point3d & p1, point3d &p2, point3d &p3);
double getDihedral(point3d & p1, point3d &p2, point3d &p3, point3d &p4);
point3d getCentroid(vector<point3d> &pointCloud);
void setCoordinate(point3d &c0, point3d &c1, point3d &c2, point3d &c3, double alpha, double tau, double normw);

void loadPdb(char *filename, vector<pdbInfo> &pdb);
void loadFasta(char *filename, vector<int> &aa);
void loadSec(char *filename, vector<int> &ss);
void loadAlign(char *filename, vector<int> &alnCode); 
void loadCmap(char *filename, MDArray<double> &cm);

void pdb2pose(vector<pdbInfo> &pdb, vector<int> &ss, vector<poseInfo> &pose);
void pose2pdb(vector<poseInfo> &pose, vector<pdbInfo> &pdb);
void sample2pose(MDArray<double> &sample, vector<poseInfo> &pose);
void pose2sample(vector<poseInfo> &pose, MDArray<double> &sample);

bool isHydrophobic(int aa);
void getRrContact(MDArray<double> &cm, int end);
double getUpperBoundHarmonic(double d, double bound);
double getFade(double d, double lb, double ub, double z, double w);
double getLiwoEpsilon0(int ai, int aj);

double getScScEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
double getScBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
double getBbBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
double getRiRjEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight);
energyInfo getWeightedEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn);
void showEnergy(energyInfo &e);

double getRmsd(vector<pdbInfo> &pdb1, vector<pdbInfo> &pdb2, string mode = "ca");
void getAlignment(vector<vector<double> > &  A, vector<vector<double> > & B, double rot[3][3], double trans[3]);
void quat2mat(double q[4], double mat[3][3]);
void vApply (double m[3][3], const double a[3], double b[3]);

void assembleFoldon(vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow);
void assembleDomain(vector<MDArray<double> > &foldonSample, int i,  vector<int> &alnCode, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow);
void doSimulatedAnneling(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp);
bool acceptMove(energyInfo &prev, energyInfo &curr, double temperature);

void writePdb(vector<pdbInfo> &pdb, char * filename);


void runPulcha(char * pdbTempFile, char * pdbTempFile_pulchra);
void runPulcha2(char * pdbTempFile, char * pdbTempFile_pulchra);
void runPulcha3(char * pdbTempFile, char * pdbTempFile_pulchra, string pdbString, string &pdbString_pulchar, int savefile);
double getChiEnergy(char *filename, char *SAXSfilename, double weight);  
void showEnergy_Print(energyInfo &e, vector<pdbInfo> &pdb, int decoyId, int startPos, int endPos);
void showEnergy_Print_regularize(energyInfo &e, vector<pdbInfo> &pdb, int decoyId, int startPos, int endPos);
void assembleDomainSAXS_multiDomain(vector<MDArray<double> > &foldonSample, int i,  vector<int> &alnCode, vector<domLinkerInfo> linkers_info_array , DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow, int decoyId, int scoreType, bool combineScore); 
void doSimulatedAnnelingSAXS_multiDomain(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, vector<domLinkerInfo> linkers_info_array , int minLength, int maxLength, double initTemp,  double finalTemp, int decoyId, int scoreType, bool combineScore);
void doSimulatedAnnelingSAXS_multiDomain_regularized(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, vector<domLinkerInfo> linkers_info_array , int minLength, int maxLength, double initTemp,  double finalTemp, int decoyId, int scoreType, bool combineScore);
void doSimulatedAnnelingSAXS_multiDomain_regularized_pro(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, vector<domLinkerInfo> linkers_info_array , int minLength, int maxLength, double initTemp,  double finalTemp, int decoyId, int scoreType, bool combineScore);
void removeTempFile (char *Prefix); 
energyInfo getWeightedEnergy_saxs_unknown_inLoop(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn, vector<double> &SXS_chi_scores,double proceeding_ratio);
energyInfo getWeightedEnergy_saxs_unknown_inLoop_pro(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn, vector<double> &SXS_chi_scores,double proceeding_ratio, int scoreType, bool combineScore);
energyInfo getWeightedEnergy_saxs_unknown(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn); 
energyInfo getWeightedEnergy_saxs_pro(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn, int scoreType, bool combineScore); 
energyInfo getWeightedEnergy_saxs_regularize_pro(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn, int scoreType, bool combineScore);
void sample2file(MDArray<double> &sample,  char * filename) ;
double run_foxs(int foxsargc, char** foxsargv);
double getChiEnergy_inside(char *filename, char *SAXSfilename, double weight);
double getChiEnergy_inside2(char *filename, char *SAXSfilename, string pdbString, double weight) ;
bool acceptMove_jump(energyInfo &prev, energyInfo &curr, double temperature, double weight);
double run_foxs_saxs(int foxsargc, char** foxsargv, string inputString);
int run_Pulcha_int(int pulchaargc, char** pulchaargv); 
int run_Pulcha_int2(int pulchaargc, char** pulchaargv,string pdbString, string &pdbString_pulchar,int savefile);  
void Pdb2String(vector<pdbInfo> &pdb, string &inputString);
double run_foxs_only(int foxsargc, char** foxsargv);
void printHelpInfo(int argc, char ** argv);
void parseOption(int argc, char ** argv);
void parseOption_special(int argc, char ** argv);
void gen_initial_model(vector<int> &aa,vector<int> &ss, MDArray<double> &cm,vector<pdbInfo> &native,DBN &dbn, vector<pdbInfo> &pdbFoldon);
void backboneSampling(vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow);
void backboneModelling(vector<int> &aa,vector<int> &ss,MDArray<double> &sample,DBN &dbn, vector<pdbInfo> &pdbFoldon);
saxsInfo getChiEnergy_inside_pro(char *filename, char *SAXSfilename, string pdbString,double weight);
saxsInfo run_foxs_saxs_pro(int foxsargc, char** foxsargv, string inputString);
/**
*	compares two vector lengths
*/
template<class T>
void comp(vector<T> v1, vector<T> v2){
	if(v1.size() != v2.size()){
		cout << "vectors not the same size\n";
		exit(1);
	}
}

/**
*	Euclidean distance between two vectors of type T such that T has binary +,-,* 
*/
template<class T>
double euclidean(vector<T> v1, vector<T> v2){
	comp(v1, v2);
	T diff, sum;

	diff = v1[0] - v2[0];
	sum = diff * diff;

	for (unsigned int i=1; i < v1.size(); i++){
		diff = v1[i] - v2[i];
		sum += diff * diff;
	}
	return sqrt(double(sum));
}

/**
*	Jaccard Coefficient.	Use for asymetric binary values
*/
template<class T>
double jaccard(vector<T> v1, vector<T> v2){
	comp(v1, v2);
	int f11 = 0, f00 = 0;

	for (unsigned int i=0; i < v1.size(); i++){
		if(v1[i] == v2[i]){
			if(v1[i])
				f11++;
			else
				f00++;
		}
	}
	return double(f11) / double(v1.size() - (f11+f00));
}

/**
*	Cosine
*/
template<class T>
double cosine(vector<T> v1, vector<T> v2){
	comp(v1, v2);

	T lv1 = v1[0] * v1[0];
	T lv2 = v2[0] * v2[0];
	T dot12 = v1[0] * v2[0];

	for (unsigned int i=0; i < v1.size(); i++){
		lv1 += v1[i] * v1[i];
		lv2 += v2[i] * v2[i];
		dot12 += v1[i] * v2[i];
	}
	return double(dot12) / ( sqrt(double(lv1)) * sqrt(double(lv2)) );
}

/**
*	The mean of a vector
*/
template<class T>
double mean(vector<T> v1){
	T sum = v1[0];
	for (unsigned int i=1; i < v1.size(); i++)
		sum += v1[i];
	return double(sum) / double(v1.size());
}

/**
*	The Covariance
*/
template<class T>
double covariance(vector<T> v1, vector<T> v2){
	comp(v1, v2);
	double mean1 = mean(v1), mean2 = mean(v2);
	double sum = (double(v1[0]) - mean1) * (double(v2[0]) - mean2);
	for (unsigned int i=1; i < v1.size(); i++){
		sum += (double(v1[i]) - mean1) * (double(v2[i]) - mean2);
	}
	return double(sum) / double(v1.size()-1);
}

/**
*	standard deviation the covariance where both vectors are the same.
*/
template<class T>
double std_dev(vector<T> v1){
	return sqrt(covariance(v1, v1));
}

/**
*	Pearson Correlation
*/
template<class T>
double pearson(vector<T> v1, vector<T> v2){
	if (std_dev(v1) * std_dev(v2) == 0){
		cout << "( a standard deviaton was 0 )";
		return -2; // I dont know what to do here???
	}
	return covariance(v1,v2) / ( std_dev(v1) * std_dev(v2));
}

/*************************************************************************
 * Name        : main
 * Purpose     : United residue protein folding via stepwise, conditional sampling
 * Arguments   : int argc, char ** argv
 * Return Type : int
 *************************************************************************/
int main(int argc, char ** argv) {
    cout << endl;
    cout << "###########################################################################" << endl;
    cout << "#                                UniCon3D                                 #" << endl;
    cout << "#                              Version: 1.0                               #" << endl;
    cout << "#         De novo protein structure prediction using united-residue       #" << endl;
    cout << "#        conformational search via stepwise, probabilistic sampling       #" << endl;
    cout << "#         Copyright (C) Debswapna Bhattacharya and Jianlin Cheng          #" << endl;
    cout << "#     Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory      #" << endl;
    cout << "#              University of Missouri, Columbia MO 65211                  #" << endl;
    cout << "###########################################################################" << endl;
    
    parseCommandLine(argc, argv);
    
	// if run fox only 
	if(foxsonly) {
		
		return EXIT_SUCCESS;
	}
	// if run pulchar only 
	if(pulcharonly) {
		
		return EXIT_SUCCESS;
	}
	
	
	if(testonly.compare("checkparameter") == 0 )
	{
		cout << "Testing checkparameter finish! Stop here" << endl;
		exit(0);
		
	}
    // load model
    DBN dbn;
    dbn.load(moFile);
    dbn.randomGen->get_rand();
    
    // load amino acid
    vector<int> aa;
    loadFasta(faFile, aa);
	
	cout << "Fasta:  ";
	string FullFastaString="";
	for (int i = 0; i < aa.size(); i++) {
		//cout << seq[aa[i]];
		FullFastaString.append(seq[aa[i]]);
		FullFastaArray.push_back(seq[aa[i]]);
	}
	cout << FullFastaString << endl;

    // load secondary structure
    vector<int> ss;
    loadSec(ssFile, ss);
	
	cout << "ss:     ";
	for (int i = 0; i < ss.size(); i++) {
		cout << sec[ss[i]];
		FullSSArray.push_back(seq[aa[i]]);
	}
	cout << endl;
	
	
    // load alignment code
    vector<int> DomainAlnCode;
	for (int i = 0; i < aa.size(); i++) {
		//cout << seq[aa[i]];
		DomainAlnCode.push_back(0);
	}
	
 
    //load contact matrix
    MDArray<double> cm;
    cm.set_shape(aa.size(), aa.size());
    loadCmap(cmFile, cm);
    
	
	
    // load native pdb if provided
    vector<pdbInfo> native;
    vector<poseInfo> nativePose;
    //MDArray<double> pdbSampleFile;
    if (nFile) {
        loadPdb(natFile, native);
        if (aa.size() != native.size()) {
            cout << "Error! size mismatch in native structure" << endl;
            exit(0);
        }
        pdb2pose(native, ss, nativePose);
        //pose2sample(nativePose,pdbSampleFile);
		
		if(odir)
		{
			sprintf(pdbfile_initial, "%s/%s_initial_native_model.pdb", outputdir, jobId);
			sprintf(pdbfile_initial_pulchra, "%s/%s_initial_native_model.rebuilt.pdb", outputdir, jobId);
		}else{
			sprintf(pdbfile_initial, "%s_initial_native_model.pdb", jobId);
			sprintf(pdbfile_initial_pulchra, "%s_initial_native_model.rebuilt.pdb", jobId);
		}
		
		writePdb(native, pdbfile_initial);	

		
		// run pulcha on global structure
		//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);		
		string dompdbString;
		string dompdbString_pulchar;    
		//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
		Pdb2String(native, dompdbString);	
		runPulcha3(pdbfile_initial, pdbfile_initial_pulchra,dompdbString,dompdbString_pulchar,1);	
		//runPulcha2(pdbdomainFile, dompdbfile_initial_pulchra);	
		//vector<int> alnCode = DomainAlnCode		
		
		vector<pdbInfo> pdb_native_structure_tmp;
		pose2pdb(nativePose, pdb_native_structure_tmp);
		double rmsd;
        rmsd = getRmsd(pdb_native_structure_tmp, native);
        cout << "UniCon3D.io: Ca_rmsd to native structure after loading structure = " << rmsd << endl;
        if(rmsd> 2)
		{
			cout << "Warning: Large RMSD variation after loading native, need refine first!" << endl;
			//exit(0);
		}
   
    }
	
	
	
	// if run evaluate only 
	if(evaluateonly) {
		
		if(eFile)
		{		
			energyInfo e;
			e.sc_sc = getScScEnergy(native, nativePose, 1.0);
			e.sc_bb = getScBbEnergy(native, nativePose, 1.0);
			e.bb_bb = getBbBbEnergy(native, nativePose, 1.0);
			e.ri_rj = getRiRjEnergy(native, nativePose, 1.0);
			
			int idrand = rand();
			
			
			if(odir)
			{
				sprintf(pdbTempFile, "%s/%s%d.pdb", outputdir, "GlobalFoldon_",idrand);
				sprintf(pdbTempFile_prefix, "%s/%s%d", outputdir, "GlobalFoldon_",idrand);
				sprintf(pdbTempFile_pulchra, "%s/%s.rebuilt.pdb", outputdir,pdbTempFile_prefix);

			}else{
				sprintf(pdbTempFile, "%s%d.pdb", "GlobalFoldon_",idrand);
				sprintf(pdbTempFile_prefix, "%s%d", "GlobalFoldon_",idrand);
				sprintf(pdbTempFile_pulchra, "%s.rebuilt.pdb",pdbTempFile_prefix);

			}
			//writePdb(pdb, pdbTempFile);
			pdbString.clear();
			pdbString_pulchar.clear();
			Pdb2String(native, pdbString); // convert pdb file to string
			runPulcha3(pdbTempFile, pdbTempFile_pulchra,pdbString,pdbString_pulchar,1);

			
			saxsInfo saxs_score=getChiEnergy_inside_pro(pdbTempFile,saxsFile,pdbString_pulchar,1.0);  
			e.saxs_chi_global = saxs_score.chi_score;  	
			e.saxs_KL =saxs_score.KL;  	
			e.saxs_score2 = saxs_score.fun_score2;  	
			e.saxs_RG_normalize = saxs_score.RGdiff_normalize;  
			
			
		
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_chi * e.saxs_chi_global+ w_saxs_KL * e.saxs_KL + w_saxs_RG_normalize * e.saxs_RG_normalize; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy =  w_saxs_chi * e.saxs_chi_global+ w_saxs_KL * e.saxs_KL+ w_saxs_RG_normalize * e.saxs_RG_normalize; 
		
			e.saxs_penalty=0;
				
			if(odir)
			{
				sprintf(statScoreFile, "%s/energy_score.txt", outputdir);

			}else{
				sprintf(statScoreFile, "energy_score.txt");

			}	
			
			//cout << "Writing score to file: " << statScoreFile << endl;	  
			ofstream statsaxs(statScoreFile);
			if (!statsaxs.is_open()) {
				cout << "Error! folding statistics file can not open " << statScoreFile << endl;
				exit(0);
			}
			
			// buffer statistics
			char bufStat[1000];
			sprintf(bufStat,"Total\tEnergy_structure\tEnergy_SAXS\tChi-score\tSAXS_KL\tSAXS_RG");
			statsaxs << bufStat << endl;

			sprintf(bufStat,"%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f",e.total, e.structure_energy, e.saxs_energy, e.saxs_chi_global, e.saxs_KL, e.saxs_RG_normalize);
			statsaxs << bufStat << endl;
				
		}else{
			energyInfo energyNext = getWeightedEnergy(native, nativePose, dbn);
				
			if(odir)
			{
				sprintf(statScoreFile, "%s/energy_score.txt", outputdir);

			}else{
				sprintf(statScoreFile, "energy_score.txt");

			}	

			
			//cout << "Writing score to file: " << statScoreFile << endl;	  
			ofstream statsaxs(statScoreFile);
			if (!statsaxs.is_open()) {
				cout << "Error! folding statistics file can not open " << statScoreFile << endl;
				exit(0);
			}

			// buffer statistics
			char bufStat[1000];
			sprintf(bufStat,"Total");
			statsaxs << bufStat << endl;

			sprintf(bufStat,"%8.3f",energyNext.total);
			statsaxs << bufStat << endl;
		  
		
		}
		
		return EXIT_SUCCESS;
	
	}
	
	
		   
    // load initial pdb if provided
    vector<pdbInfo> initialstr;
    vector<poseInfo> initialstrPose;
    MDArray<double> pdbSampleFile;
    if (zFile) {
        loadPdb(zeroFile, initialstr);
        if (aa.size() != initialstr.size()) {
            cout << "Error! size mismatch in initial structure" << endl;
            exit(0);
        }
        pdb2pose(initialstr, ss, initialstrPose);
        pose2sample(initialstrPose,pdbSampleFile);
   
    }else{  
        cout << "The initial modeling file couldn't be found, need generate by simulation" << endl;
        //exit(0);
    }	
	
   // extract secondary structural elements
    vector<sseInfo> sse;
    bool sseStarted = false;
    sseInfo sseData;
    for (int i = 0; i < ss.size() - 1; i++) {
        if (!sseStarted) {
            sseData.start = i;
            sseData.type = ss[i];
            sseStarted = true;
        }
        if (sseStarted && ss[i] != ss[i+1]) {
            sseData.end = i;
            sse.push_back(sseData);
            sseStarted = false;
        }
    }
    // for the last sse
    sseData.type = ss[ss.size() - 1];
    sseData.end = ss.size() - 1;
    sse.push_back(sseData);
    		
		
/*************************************************************************
 * Added by Jie for initial model generation
 *************************************************************************/

    // Initial model generation
    vector<sseInfo> foldon;
    sseInfo foldonData_initial;
	foldonData_initial.start = sse[0].start;
	foldonData_initial.end = sse[sse.size() - 1].end;
	foldonData_initial.type = sse[sse.size() - 1].type;
	foldon.push_back(foldonData_initial);
    foldon[foldon.size() - 1].end = sse[sse.size()-1].end;
 /*************************************************************************
 * Added by Jie for initial input model
 *************************************************************************/
    vector<pdbInfo> inital_full_model;
	MDArray<double> pdbSampleFile_initial;
	
    if (zFile) {
		pdbSampleFile_initial = pdbSampleFile;
		vector<poseInfo> pose_initial;
		sample2pose(pdbSampleFile_initial, pose_initial);
		vector<pdbInfo> pdb_initial;

		if(odir)
		{
			sprintf(pdbfile_initial, "%s/%s_%.6d_initial_check.pdb",outputdir, jobId, 1);
			sprintf(pdbfile_initial_pulchra, "%s/%s_%.6d_initial_check.rebuilt.pdb",outputdir, jobId, 1);
		}else{
			sprintf(pdbfile_initial, "%s_%.6d_initial_check.pdb", jobId, 1);
			sprintf(pdbfile_initial_pulchra, "%s_%.6d_initial_check.rebuilt.pdb", jobId, 1);
		}

		
		writePdb(pdb_initial, pdbfile_initial);
	   	

		
		// run pulcha on global structure
		//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);
		pdbString.clear();
		pdbString_pulchar.clear();
		//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
		Pdb2String(pdb_initial, pdbString);	
		runPulcha3(pdbfile_initial, pdbfile_initial_pulchra,pdbString,pdbString_pulchar,1);	

		inital_full_model = pdb_initial;
		
	
	
    }else{
	
	
			
		 /*************************************************************************
		 * Added by Jie for initial input model, calculating RMSD between native and initial model if native provided
		 *************************************************************************/
		 
		
		 /*************************************************************************
		 * Added by Jie for generating initial model 
		 *************************************************************************/
		int numCycles_default = numCycles;
		if(genInitial)
		{
			numCycles = genInitial_numDecoys;  // set to temporary epoch to fast generate initial model
		}else{
			
			numCycles = 1;  // set to temporary epoch to fast generate initial model
		}
		
		vector<pdbInfo> simulated_inital_model;
		
		gen_initial_model(aa,ss,cm, native, dbn, simulated_inital_model);
		
		
	

		if(odir)
		{
			
			sprintf(pdbSimulatedFile, "%s/%s_simulated_initial_model.pdb",outputdir, jobId);
			sprintf(simpdbfile_initial_pulchra, "%s/%s_simulated_initial_model.rebuilt.pdb",outputdir, jobId);
		}else{
			sprintf(pdbSimulatedFile, "%s_simulated_initial_model.pdb", jobId);
			sprintf(simpdbfile_initial_pulchra, "%s_simulated_initial_model.rebuilt.pdb", jobId);
		}
					
		writePdb(simulated_inital_model, pdbSimulatedFile);


		
		// run pulcha on global structure
		//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);		
		string simpdbString;
		string simpdbString_pulchar;    
		//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
		Pdb2String(simulated_inital_model, simpdbString);	
		runPulcha3(pdbSimulatedFile, simpdbfile_initial_pulchra,simpdbString,simpdbString_pulchar,1);
		//runPulcha2(pdbSimulatedFile, simpdbfile_initial_pulchra);
		//runPulcha(pdbSimulatedFile, simpdbfile_initial_pulchra);
		numCycles  = numCycles_default;	
		
		inital_full_model = simulated_inital_model;
		
		if(genInitial)
		{
			cout << "Initial model generation finished! Saved to " << pdbSimulatedFile << endl;
			return EXIT_SUCCESS;
		}
	}
	
	 /*************************************************************************
	 * Added by Jie for loading domain models
	 *************************************************************************/
	 
	
	vector<poseInfo> full_domainPose;
	pdb2pose(inital_full_model, ss, full_domainPose);
	
			
	
	string line, str;
    string header (">");
    ifstream fin (dmFile);
	int domain_id=0;
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            if(line.length() != 0 && line.compare(0, header.length(), header) != 0) {
				ifstream fin2( line.c_str());
				if( fin2.fail() ) {
					cout << endl;
					cout << "Error! Domain file " <<  line <<" not present" << endl << endl;
					exit(0);
				}else{
					
					cout << "Decting domain file " <<  line << endl;
					domain_id++;
					// start loading the domain pdb
					// first, extract sequence from domain pdb, get domain start compared to fasta, then assign angles to initial model
					vector<int> domain_aa;
					
					vector<pdbInfo> domainstr;
					char * dmpdbpath = const_cast<char*> (line.c_str());
					loadPdb(dmpdbpath, domainstr);
					for (int i = 0; i < domainstr.size(); i++) {
						int aa_code = domainstr[i].aa;
						domain_aa.push_back(aa_code);
					}
					
					
					
					
					char pdbdomainFile_tmp0[1000];	
					char dompdbfile_initial_pulchra_tmp0[1000];	
					if(odir)
					{
						sprintf(pdbdomainFile_tmp0, "%s/%s_initial_domain_D%i_model_before.pdb",outputdir, jobId,domain_id);
						sprintf(dompdbfile_initial_pulchra_tmp0, "%s/%s_initial_domain_D%i_model_before.rebuilt.pdb",outputdir, jobId,domain_id);
					}else{
						sprintf(pdbdomainFile_tmp0, "%s_initial_domain_D%i_model_before.pdb", jobId,domain_id);
						sprintf(dompdbfile_initial_pulchra_tmp0, "%s_initial_domain_D%i_model_before.rebuilt.pdb", jobId,domain_id);
					}
					writePdb(domainstr, pdbdomainFile_tmp0);
						
					// run pulcha on global structure
					//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);		
					string dompdbString_tmp0;
					string dompdbString_pulchar_tmp0;    
					//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
					Pdb2String(domainstr, dompdbString_tmp0);	
					runPulcha3(pdbdomainFile_tmp0, dompdbfile_initial_pulchra_tmp0,dompdbString_tmp0,dompdbString_pulchar_tmp0,1);	
					
					cout << "Domain fasta:  ";
					string domainString="";
					for (int i = 0; i < domain_aa.size(); i++) {
						//cout << seq[domain_aa[i]];
						domainString.append(seq[domain_aa[i]]);
					}
					cout << domainString << endl;
					
					std::size_t domain_found = FullFastaString.find(domainString);
					if (domain_found!=std::string::npos)
					{
					 	std::cout << "Found domain start at: " << domain_found << " end at : "<< domain_found + domainString.length() -1 <<endl;
					}else{
						
						std::cout << "Error! The fasta sequence in domain "<< line <<" is not subsequence of full fasta"  << endl;
					}
					
					cout << endl;
					
					int domain_start = domain_found;
					int domain_end= domain_found + domainString.length() -1;
					
					
					cout << "Domain ss:  ";
					string domainSS="";
					vector<int> domainss_arr;
					for (int i = domain_start; i < domain_end+1; i++) {
						//cout << seq[domain_aa[i]];
						domainSS.append(sec[ss[i]]);
						domainss_arr.push_back(ss[i]);
						 
						DomainAlnCode[i]=1;
					}
					if(domain_end!=FullFastaString.size()-1)
					{
						DomainAlnCode[domain_end]=0; // in case two domain is completely connected  1-56, 57,90
						
					}
					cout << domainSS << endl;
					
					
					// start load domain model into initial model using pose function 
					vector<poseInfo> domainPose;
					pdb2pose(domainstr, domainss_arr, domainPose);
					
					
					vector<pdbInfo> pdb_assembly_structure_tmp;
					pose2pdb(domainPose, pdb_assembly_structure_tmp);
					
					char pdbdomainFile_tmp[1000];	
					char dompdbfile_initial_pulchra_tmp[1000];	
					
					
					if(odir)
					{
						sprintf(pdbdomainFile_tmp, "%s/%s_initial_domain_D%i_model.pdb",outputdir, jobId,domain_id);
						sprintf(dompdbfile_initial_pulchra_tmp, "%s/%s_initial_domain_D%i_model.rebuilt.pdb",outputdir, jobId,domain_id);
					}else{
						sprintf(pdbdomainFile_tmp, "%s_initial_domain_D%i_model.pdb", jobId,domain_id);
						sprintf(dompdbfile_initial_pulchra_tmp, "%s_initial_domain_D%i_model.rebuilt.pdb", jobId,domain_id);
					}
					writePdb(pdb_assembly_structure_tmp, pdbdomainFile_tmp);
						
					// run pulcha on global structure
					//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);		
					string dompdbString_tmp;
					string dompdbString_pulchar_tmp;    
					//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
					Pdb2String(pdb_assembly_structure_tmp, dompdbString_tmp);	
					runPulcha3(pdbdomainFile_tmp, dompdbfile_initial_pulchra_tmp,dompdbString_tmp,dompdbString_pulchar_tmp,1);		
					//runPulcha2(pdbdomainFile_tmp, dompdbfile_initial_pulchra_tmp);		
					//runPulcha(pdbdomainFile_tmp, dompdbfile_initial_pulchra_tmp);		
					
					if(false) // not needed for backbone samplning, the same
					{
						/*************************************************************************
						 * Added by Jie for modeling the backbone distance
						 *************************************************************************/
						int numCycles_default_tmp = numCycles;
						numCycles = 1;  // set to temporary epoch to fast generate initial model
						vector<pdbInfo> simulated_inital_model_tmp;
						MDArray<double> pdbSampleFile_tmp;
						pose2sample(domainPose,pdbSampleFile_tmp);
						backboneModelling(domain_aa,domainss_arr, pdbSampleFile_tmp, dbn, simulated_inital_model_tmp);
						
						
						char pdbSimulatedFile_tmp[1000];
						char simpdbfile_initial_pulchra_tmp[1000];	
						if(odir)
						{
							sprintf(pdbSimulatedFile_tmp, "%s/%s_initial_domain_D%i_model_resampled.pdb",outputdir, jobId,domain_id);
							sprintf(simpdbfile_initial_pulchra_tmp, "%s/%s_initial_domain_D%i_model_resampled.rebuilt.pdb",outputdir, jobId,domain_id);
						}else{
							sprintf(pdbSimulatedFile_tmp, "%s_initial_domain_D%i_model_resampled.pdb", jobId,domain_id);
							sprintf(simpdbfile_initial_pulchra_tmp, "%s_initial_domain_D%i_model_resampled.rebuilt.pdb", jobId,domain_id);
						}
						
						writePdb(simulated_inital_model_tmp, pdbSimulatedFile_tmp);


						
						// run pulcha on global structure
						//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);		
						string simpdbString_tmp;
						string simpdbString_pulchar_tmp;    
						//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
						Pdb2String(simulated_inital_model_tmp, simpdbString_tmp);	
						runPulcha3(pdbSimulatedFile_tmp, simpdbfile_initial_pulchra_tmp,simpdbString_tmp,simpdbString_pulchar_tmp,1);		
						numCycles  = numCycles_default_tmp;
					}
					
					// update domains 
					for (int i = 0; i < domainPose.size(); i++) {
						
						if(full_domainPose[i+domain_start].aa == domainPose[i].aa and full_domainPose[i+domain_start].ss == domainPose[i].ss)
						{
							//full_domainPose[i+domain_start].id = domainPose[i].id; // no id, need use default full pdb id
							full_domainPose[i+domain_start].bca = domainPose[i].bca;
							full_domainPose[i+domain_start].tao = domainPose[i].tao;
							full_domainPose[i+domain_start].theta = domainPose[i].theta;
							full_domainPose[i+domain_start].bsc = domainPose[i].bsc;
							full_domainPose[i+domain_start].phi = domainPose[i].phi;
							full_domainPose[i+domain_start].delta = domainPose[i].delta;
						}else{
							
							cout << "Error: The aa in pdb and domain not match" << endl;
							exit(0);
						}
					}
				}
            }
        }
        fin.close();
    }
    else {
        cout << "Error! domain file can not open " << dmFile << endl;
        exit(0);
    }	 
	
	
	
	vector<pdbInfo> pdb_assembly_structure;
	pose2pdb(full_domainPose, pdb_assembly_structure);
	if(odir)
	{
		sprintf(pdbdomainFile, "%s/%s_initial_domain_assembly_model.pdb", outputdir, jobId);
		sprintf(dompdbfile_initial_pulchra, "%s/%s_initial_domain_assembly_model.rebuilt.pdb", outputdir, jobId);
	}else{
		sprintf(pdbdomainFile, "%s_initial_domain_assembly_model.pdb", jobId);
		sprintf(dompdbfile_initial_pulchra, "%s_initial_domain_assembly_model.rebuilt.pdb", jobId);
	}
    
    writePdb(pdb_assembly_structure, pdbdomainFile);	

	
	// run pulcha on global structure
	//runPulcha2(pdbfile_initial, pdbfile_initial_pulchra);		
	string dompdbString;
	string dompdbString_pulchar;    
	//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
	Pdb2String(pdb_assembly_structure, dompdbString);	
	runPulcha3(pdbdomainFile, dompdbfile_initial_pulchra,dompdbString,dompdbString_pulchar,1);	
	//runPulcha2(pdbdomainFile, dompdbfile_initial_pulchra);	
	//vector<int> alnCode = DomainAlnCode		
		
	/*
	char** Files;
	Files = new char*[1000];// initialize the double pointer
	Files[0]=new char[1000];// initialize 1st char*, with capacity of 4 chars
	Files[1]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[2]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
		
	strcpy(Files[0],"pulchra");//copy some data to 1st string
	strcpy(Files[1],"-g");//copy some data to 2nd string
	strcpy(Files[2],pdbdomainFile);//copy some data to 2nd string
	run_Pulcha_int(3,Files);
	*/
	
	
	
	
	cout << "DomainAlnCode: " << endl;
	int linker_num=0;
	for (int i = 0; i < DomainAlnCode.size(); i++) {
		cout << DomainAlnCode[i];
		if(DomainAlnCode[i] ==0)
		{
			linker_num++;
		}
	}
	cout << endl;
	
	if(linker_num >0)
	{
		cout << linker_num <<" residues are not aligned!" << endl;
	}else{
		cout << "Warning: None unaligned residues are found! Not sure where is the domain linker!" << endl;
		exit(0);
	}


    if (DomainAlnCode.size() != aa.size()) {
        cout << "Error! size mismatch in fasta, alignment" << endl;
        exit(0);
    }
	
	// define domain linker start and end, by Jie 
    // extract domain linker region, 1111111111111111111111000000000000011111111111111111
	
	domLinkerInfo linker_info;
	vector<domLinkerInfo> linkers_info_array;
	
	int link_start = 0;
    int link_end = 0;
	bool link_start_done = false;
    for (int i = 0; i < DomainAlnCode.size(); i++) {
		if(DomainAlnCode[i] == 0)
		{
			if (!link_start_done) {
				link_start = i;
				link_start_done = true;
			}
			if (link_start_done && DomainAlnCode[i] != DomainAlnCode[i+1]) {
				link_end = i;
				link_start_done = false;
				linker_info.start = link_start;
				linker_info.end = link_end;
				linkers_info_array.push_back(linker_info);
			}
		}	
    }
	
	
	
	cout << linkers_info_array.size() +1 << " domains are found for assembling" << endl;
	
	for (int i = 0; i < linkers_info_array.size(); i++) {
		link_start = linkers_info_array[i].start;
		link_end = linkers_info_array[i].end;
		cout << "Linker " << i << " start: " << link_start << endl;;
		cout << "Linker " << i << " end:   " << link_end << endl;
		
		if(link_end - link_start + 1 == 1 )
		{
			// set 5 redidues  for sampling
			//for (int j = link_start-adjustLinker; j < link_start+adjustLinker; j++) {
			//	int cur_ind=j;
			//	if(cur_ind<0)
			//	{
			//		cur_ind=0;
			//	}
			//	DomainAlnCode[cur_ind]=0;
			//}
			//for (int j = link_end-adjustLinker; j < link_end+adjustLinker+1; j++) {
			//	int cur_ind=j;
			//	if(cur_ind<0)
			//	{
			//		cur_ind=0;
			//	}
			//	DomainAlnCode[cur_ind]=0;
			//}
			DomainAlnCode[link_start-2]=0;
			DomainAlnCode[link_start-1]=0;
			DomainAlnCode[link_start]=0;
			DomainAlnCode[link_end+1]=0;
			DomainAlnCode[link_end+2]=0;
		}else if(link_end - link_start + 1 == 2 )
		{
			DomainAlnCode[link_start-2]=0;
			DomainAlnCode[link_start-1]=0;
			DomainAlnCode[link_start]=0;
			DomainAlnCode[link_end+1]=0;
		}else if(link_end - link_start + 1 == 3 )
		{
			DomainAlnCode[link_start-1]=0;
			DomainAlnCode[link_start]=0;
			DomainAlnCode[link_end+1]=0;
		}else if(link_end - link_start + 1 == 4 )
		{
			DomainAlnCode[link_start-1]=0;
			DomainAlnCode[link_start]=0;
		}
	}
	
	
	cout << "Extended DomainAlnCode for sampling: " << endl;
	for (int i = 0; i < DomainAlnCode.size(); i++) {
		cout << DomainAlnCode[i];
	}
	cout << endl;
	
	
	domLinkerInfo linker_info_new;
	vector<domLinkerInfo> linkers_info_array_new;
	
	link_start = 0;
    link_end = 0;
	link_start_done = false;
    for (int i = 0; i < DomainAlnCode.size(); i++) {
		if(DomainAlnCode[i] == 0)
		{
			if (!link_start_done) {
				link_start = i;
				link_start_done = true;
			}
			if (link_start_done && DomainAlnCode[i] != DomainAlnCode[i+1]) {
				link_end = i;
				link_start_done = false;
				linker_info_new.start = link_start;
				linker_info_new.end = link_end;
				linkers_info_array_new.push_back(linker_info_new);
			}
		}	
    }
	
	cout << "Before ss:     "<< endl;
	for (int i = 0; i < ss.size(); i++) {
		cout << sec[ss[i]];
	}
	cout << endl;
	for (int i = 0; i < linkers_info_array_new.size(); i++) {
		link_start = linkers_info_array_new[i].start;
		link_end = linkers_info_array_new[i].end;
		cout << "Linker " << i << " start: " << link_start << endl;;
		cout << "Linker " << i << " end:   " << link_end << endl;
		
		// set the linker secondary structure as loop, set ss code to 7
		 for (int t = link_start; t <= link_end; t++) {
			ss[t] = 7;
			cout << "converting initial ss from " <<full_domainPose[t].ss << " to 7"<<endl;
			full_domainPose[t].ss = 7;
			cout << "converting native ss from " <<nativePose[t].ss << " to 7"<<endl;
			nativePose[t].ss=7;
		 }
	}
	
	cout << "updated ss:     "<< endl;
	for (int i = 0; i < ss.size(); i++) {
		cout << sec[ss[i]];
	}
	cout << endl;
	
	
	
	MDArray<double> FullpdbSampleFile;
    pose2sample(full_domainPose,FullpdbSampleFile);
		
	
 /*************************************************************************
 * Added by Jie for assigning initial sample two first structure for folding
 *************************************************************************/      
    // populate sample corresponding to foldon units
    vector<MDArray<double> > foldonSample;
    for (int i = 0; i < foldon.size(); i++) {
        MDArray<double> foldonSampleData;
        int foldonSize = foldon[i].end + 1;
        foldonSampleData.set_shape(2 * foldonSize, NUM_DAT);
        for (int j = 0; j <= foldon[i].end; j++) {
            // for backbone
            foldonSampleData.set(2 * j, 0, FullpdbSampleFile.get(2 * j, 0));
            foldonSampleData.set(2 * j, 1, FullpdbSampleFile.get(2 * j, 1));
			if(FullpdbSampleFile.get(2 * j, 2) != aa[j])
			{
				cout << " assign initial model: " << FullpdbSampleFile.get(2 * j, 2) << " not equal to "<< aa[j] << " at residue aa[j] " << j << endl;
				exit(0);	
			}
            foldonSampleData.set(2 * j, 2, FullpdbSampleFile.get(2 * j, 2));
			if(FullpdbSampleFile.get(2 * j, 3) != ss[j])
			{
				cout << " assign initial model: " << FullpdbSampleFile.get(2 * j, 3) << " not equal to "<< ss[j] << " at residue ss[j] " << j << endl;
				exit(0);	
			}
            foldonSampleData.set(2 * j, 3, FullpdbSampleFile.get(2 * j, 3));
            foldonSampleData.set(2 * j, 4, FullpdbSampleFile.get(2 * j, 4));
            foldonSampleData.set(2 * j, 5, FullpdbSampleFile.get(2 * j, 5));
            foldonSampleData.set(2 * j, 6, FullpdbSampleFile.get(2 * j, 6));
            // for side chain
            foldonSampleData.set(2 * j + 1, 0, FullpdbSampleFile.get(2 * j + 1, 0));
            foldonSampleData.set(2 * j + 1, 1, FullpdbSampleFile.get(2 * j + 1, 1));
            foldonSampleData.set(2 * j + 1, 2, FullpdbSampleFile.get(2 * j + 1, 2));
            foldonSampleData.set(2 * j + 1, 3, FullpdbSampleFile.get(2 * j + 1, 3));
            foldonSampleData.set(2 * j + 1, 4, FullpdbSampleFile.get(2 * j + 1, 4));
            foldonSampleData.set(2 * j + 1, 5, FullpdbSampleFile.get(2 * j + 1, 5));
            foldonSampleData.set(2 * j + 1, 6, FullpdbSampleFile.get(2 * j + 1, 6));
        }
        foldonSample.push_back(foldonSampleData);
    }		
    // output the native information 
	//cout << "check-0" << endl;
	cout << "Comparing the initial structure to native: " << endl;
	vector<poseInfo> pose_initial_structure;
	//cout << "check0" << endl;
	sample2pose(FullpdbSampleFile, pose_initial_structure);
	//cout << "check00" << endl;
	//cout << "check1" << endl;
	vector<pdbInfo> pdb_initial_structure;
	pose2pdb(pose_initial_structure, pdb_initial_structure);
	//cout << "check2" << endl;
	//energyInfo energyNext = getWeightedEnergy(pdbNext, poseNext, dbn);
	if(eFile)
	{
		//cout <<  "saxs1" << endl;
		energyInfo energy_initial_structure= getWeightedEnergy_saxs_pro(pdb_initial_structure, pose_initial_structure, dbn,1,false);
		//cout <<  "saxs2" << endl;
		showEnergy_Print(energy_initial_structure,pdb_initial_structure, 0,0,0);
		//cout <<  "saxs3" << endl;
		// run again to see if get same results
		energy_initial_structure= getWeightedEnergy_saxs_pro(pdb_initial_structure, pose_initial_structure, dbn,1,false);
		//cout <<  "saxs4" << endl;
		showEnergy_Print(energy_initial_structure,pdb_initial_structure, 0,0,0);
		//cout <<  "saxs5" << endl;
	}else{
		energyInfo energyNext = getWeightedEnergy(pdb_initial_structure, pose_initial_structure, dbn);
		showEnergy(energyNext);
	}
	
    
    // assemble foldon units via stepwise addition
	
    vector<pdbInfo> pdbFoldon;
    vector<poseInfo> poseFoldon;
    
	
	if(testonly.compare("test_input") == 0 )
	{
		cout << "Testing finish! Stop here" << endl;
        exit(0);
		
	}
	
	

	if(odir)
	{
		sprintf(statFile, "%s/%s_stats.txt", outputdir,jobId);
	}else{
		sprintf(statFile, "%s_stats.txt", jobId);
    }
    cout << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    cout << "Job Id                   : " << jobId << endl;
    cout << "Protein Length           : " << aa.size() << endl;
    cout << "Amino Acid Sequence      : ";
    for (int i = 0; i < aa.size(); i++) {
        cout << seq[aa[i]];
    }
    cout << endl;
    cout << "Secondary Structure      : ";
    for (int i = 0; i < ss.size(); i++) {
        cout << sec[ss[i]];
    }
    cout << endl;
    cout << "Alignment file           : ";
    for (int i = 0; i < DomainAlnCode.size(); i++) {
        cout << DomainAlnCode[i];
    }
    cout << endl;
    //cout << "Domain linker sites      : " << domain_start << " ---- " << domain_end << endl;
    cout << "Number of Decoys         : " << numDecoys << endl;
    cout << "Monte Carlo Cycle Factor : " << numCycles/300 << endl;
    cout << "Native Structure         : ";
    if (nFile) {
        cout << "Available" << endl;
    }
    else {
        cout << "Not Available" << endl;
    }
    cout << "Number of Foldon Units   : " << foldon.size() << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
    
    ofstream stat(statFile);
    if (!stat.is_open()) {
        cout << "Error! folding statistics file can not open " << statFile << endl;
        exit(0);
    }

    cout << "Folden size is: "<< foldon.size() <<endl;
    for (int n = 1; n <= numDecoys; n++) {
        cout << "UniCon3D.io: Generating decoy " << n << " of " << numDecoys << endl;
		clock_t Decoy_start = clock();
        for (int i = 0; i < foldon.size(); i++) {
            
            // get residue-residue contacts from native pdb
            //getRrContact(cm, foldon[i].end);
            
            if (i == 0) {
                cout << endl <<"UniCon3D.sampling: Sampling foldon unit " << i + 1 << " [residue : 1-" << foldon[i].end + 1 << "]"<< endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            else {
                cout << endl << "UniCon3D.sampling: Conditional sampling foldon unit " << i + 1 << " [residue : " << foldon[i-1].end + 2 << "-" << foldon[i].end + 1 << "]" << endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            
            //assembleFoldon(foldonSample, i, dbn, pdbFoldon, poseFoldon);
            //assembleDomain(foldonSample, i, alnCode, dbn, pdbFoldon, poseFoldon);
			
			if(eFile)
			{
				cout << "UniCon3D.sampling: Performing simulated annealing energy minimization using SAXS" << endl;
				// to do: need consider inputted al or  from domain inside
				//assembleDomainSAXS(foldonSample, i, DomainAlnCode, domain_start, domain_end, dbn, pdbFoldon, poseFoldon,n);
				assembleDomainSAXS_multiDomain(foldonSample, i, DomainAlnCode, linkers_info_array_new, dbn, pdbFoldon, poseFoldon,n,scoreFunction,scoreCombine);
			}else{
				cout << "UniCon3D.sampling: Performing simulated annealing energy minimization" << endl;
				assembleDomain(foldonSample, i, DomainAlnCode, dbn, pdbFoldon, poseFoldon);
			}

            cout << endl;
            //sprintf(pdbFoldonFile, "%s_D%.6d_F%.6d.pdb", jobId, n, i+1);
            //if (saveInd == 0) {
            //    writePdb(pdbFoldon, pdbFoldonFile);
            //}
        }
		
		   
			
		if(odir)
		{
			sprintf(pdbFoldonFile, "%s/%s_%.6d.pdb", outputdir,jobId, n);
			sprintf(pdbFoldonFile_pulchra, "%s/%s_%.6d.rebuilt.pdb", outputdir,jobId, 1);
		}else{
			sprintf(pdbFoldonFile, "%s_%.6d.pdb", jobId, n);
			sprintf(pdbFoldonFile_pulchra, "%s_%.6d.rebuilt.pdb", jobId, 1);
		}
        
        writePdb(pdbFoldon, pdbFoldonFile);

		
		// run pulcha on global structure	    
		//runPulcha(pdbFoldonFile, pdbFoldonFile_pulchra);	
        
		string FinalpdbString;
		string FinalpdbString_pulchar;    
		//runPulcha(pdbfile_initial, pdbfile_initial_pulchra);	
		Pdb2String(pdbFoldon, FinalpdbString);	
		runPulcha3(pdbFoldonFile, pdbFoldonFile_pulchra,FinalpdbString,FinalpdbString_pulchar,1);	
		
			
        cout << endl << "UniCon3D.io: Finished generating decoy " << pdbFoldonFile << endl;
        
        
        double rmsd;
        if (nFile) {
            rmsd = getRmsd(pdbFoldon, native);
            cout << "UniCon3D.io: Ca_rmsd to native structure = " << rmsd << endl;
        }
        // prepare to write to folding statistics file
        energyInfo e = getWeightedEnergy(pdbFoldon, poseFoldon, dbn);
        if (nFile) {
            sprintf(bufStat,"Decoy_name %s E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f Ca_rmsd %8.3f", pdbFoldonFile, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total, rmsd);
            stat << bufStat << endl;
        }
        else {
            sprintf(bufStat,"Decoy_name %s E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f", pdbFoldonFile, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total);
            stat << bufStat << endl;
        }
        cout << "UniCon3D.io: Folding statistics written to " << statFile << endl;
        cout << "---------------------------------------------------------------------------" << endl;
		
		clock_t Decoy_end = clock();
		double elapsedecoy_secs1 = double(Decoy_end - Decoy_start) / CLOCKS_PER_SEC;
		cout << "UniCon3D.io: Generating decoy " << n << " of " << numDecoys << " finished within "<< elapsedecoy_secs1 << " sec!" << endl << endl;
		
    }
	
	
	if(testonly.compare("test_assembly") == 0 )
	{
		cout << "Testing assembly finish! Stop here" << endl;
        exit(0);
		
	}
    return EXIT_SUCCESS;
}

void printHelpInfo(int argc, char ** argv)
{
	cout << "Purpose: United residue protein folding via stepwise, probabilistic sampling" << endl;
	cout << "Usage: " << argv[0] << " -i id -f fasta -s ss -c cmap -m model" << endl;
	cout << "   -i id     : job id                                     : mandatory" << endl;
	cout << "   -f fasta  : fasta file                                 : mandatory" << endl;
	cout << "   -s ss     : secondary structure file                   : mandatory" << endl;
	cout << "   -c cmap   : contact map file                           : mandatory" << endl;
	cout << "   -m model  : model file                                 : mandatory" << endl;
	cout << "   -d decoy  : number of decoys to generate               : optional (default is 100)" << endl;
	cout << "   -x factor : increase monte carlo cycles by this factor : optional (default is 1)" << endl;
	cout << "   -n native : native structure for comparison            : optional" << endl;
	cout << "   -e profile : SAXS experimental profile for fitting     : optional" << endl;
	cout << "   -z initalmodel : initial structure to model            : optional" << endl;
	cout << "   -t tofitting : Fitting the SAXS profile for chi        : optional" << endl;
	cout << "   -b CAfitting : Fitting the SAXS profile for chi only CA atom        : optional" << endl;
	cout << "   -j pulchar : run pulchar program only        : optional" << endl;
	cout << "   -k foxs : run foxs program only         : optional" << endl;
	cout << "   -h help   : this message" << endl;
}
/*************************************************************************
 * Name        : parseNextItem
 * Purpose     : parse next item in command line argument
 * Arguments   : int argc, char ** argv, int & i
 * Return Type : void
 *************************************************************************/
void parseNextItem(int argc, char ** argv, int & i) {
    if (strncmp(argv[i], "-i", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No job id provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        jId = true;
        strcpy(jobId, argv[++i]);
    }
    else if (strncmp(argv[i], "-f", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No fasta file provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        fFile = true;
        strcpy(faFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-s", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No secondary structure file provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        sFile = true;
        strcpy(ssFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-a", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No alignment file provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        aFile = true;
        strcpy(alnFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-c", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No contact map file provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        cFile = true;
        strcpy(cmFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-m", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No model file provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        mFile = true;
        strcpy(moFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-d", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No number of decoys provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        numDecoys = atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-x", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No monte carlo cycle factor provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        numCycles = numCycles * atoi(argv[++i]);
    }
    else if (strncmp(argv[i], "-n", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No native structure provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        nFile = true;
        strcpy(natFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-z", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No initial structure provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        zFile = true;
        strcpy(zeroFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-h", 2) == 0) {
        cout << endl;
        printHelpInfo(argc,argv);
        exit(0);
    }
    else if (strncmp(argv[i], "-e", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! No SAXS experimental profile provided" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        eFile = true;
        strcpy(saxsFile, argv[++i]);
    }
    else if (strncmp(argv[i], "-t", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! Fitting parameter is required! (1/0)" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        
        saxsFit = atoi(argv[++i]);
		if(saxsFit == 1)
		{
			tofitting = true;
		}else if(saxsFit == 0){
			tofitting = false;
		}else{
			cout << "Error! Fitting parameter is required! (1/0)" << endl << endl;
			exit(-1);
		}
    }
    else if (strncmp(argv[i], "-b", 2) == 0) {
        if (argc < i + 2) {
            cout << endl;
            cout << "Error! Fitting using CA atom parameter is required! (1/0)" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        
        saxsCAFit = atoi(argv[++i]);
		if(saxsCAFit == 1)
		{
			CAfitting = true;
		}else if(saxsCAFit == 0){
			CAfitting = false;
		}else{
			cout << "Error! Fitting CA atom parameter is required! (1/0)" << endl << endl;
			exit(-1);
		}
    }
    else {
        cout << endl;
        cout << "Error! Invalid option" << endl << endl;
        printHelpInfo(argc,argv);
        exit(0);
    }
    i++;
}


void parseOption(int argc, char ** argv) {
	
	
	std::string jobId_s;
	std::string fasta_s;
	std::string ss_s;
	std::string aln_s;
	std::string cm_s;
	std::string mod_s;
	std::string nat_s;
	std::string init_s;
	std::string outdir_s;
	std::string saxs_s;
	std::string score_s;
	std::string score_w;
	std::string score_w_init;
	std::string score_w_final;
	std::string score_w_p_init;
	std::string score_w_p_final;
	double wsaxs_s;
	std::string test_s;
	std::string domain_s;
	int geninit_s;
	int de_s=0;
	int cy_s=0;
	int lk_s=0;
	
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("verbose,v", "print words with verbosity") 
      ("id,i", po::value<std::string>(&jobId_s), "job id (required)")
      ("fasta,f", po::value<std::string>(&fasta_s), "fasta file (required)")
      ("ss,s", po::value<std::string>(&ss_s), "secondary structure file (required)")
      ("cm,c", po::value<std::string>(&cm_s), "contact map file (required)")
      ("domainlist,l", po::value<std::string>(&domain_s), "domain path list (required)")
      ("model,m", po::value<std::string>(&mod_s), "IOHMM model file (required)")
      ("aln,a", po::value<std::string>(&aln_s), "alignment file")
	  ("numDecoys,d", po::value<int>(&de_s),"number of decoys")
	  ("numCycles,x", po::value<int>(&cy_s),"monte carlo cycle factor")
	  ("native,n", po::value<std::string>(&nat_s), "native structure file")
	  ("init,z", po::value<std::string>(&init_s), "initial structure file")
	  ("outputdir,o", po::value<std::string>(&outdir_s), "Output directory")
	  ("saxs,e", po::value<std::string>(&saxs_s), "SAXS experimental profile file")
	  ("scoreFun,F", po::value<std::string>(&score_s), "SAXS score function (chi, KL, score2, RG)")
	  ("scoreWeight", po::value<std::string>(&score_w), "weight for SAXS score function (chi, KL, score2, RG):  default: 100-100-100-100")
	  ("scoreWeightInitial", po::value<std::string>(&score_w_init), "intial weight for SAXS score function (chi, KL, score2, RG):  default: 100-100-100-100")
	  ("scoreWeightFinal", po::value<std::string>(&score_w_final), "Final weight for SAXS score function (chi, KL, score2, RG):  default: 1-1-1-1")
	  ("scoreWeightPenalyInitial", po::value<std::string>(&score_w_p_init), "intial penalty weight for SAXS score function (chi, KL, score2, RG):  default: 0.001-0.001-0.001-0.001")
	  ("scoreWeightPenalyFinal", po::value<std::string>(&score_w_p_final), "Final penalty weight for SAXS score function (chi, KL, score2, RG):  default: 10-10-10-10")
	  ("scoreCombine,B", "SAXS score function combination (chi, KL, score2, RG)")
	  ("wsaxs,w", po::value<double>(&wsaxs_s), "SAXS weight (default=50)")
	  ("testprog,g", po::value<std::string>(&test_s), "test the program(Options: checkparameter, test_input)")
      ("reg,r", "perform saxs regularization fitting (default = false)")
      ("saxsFit,t", "perform saxs fitting (default = false)")
      ("saxsCaAtomFit,b", "Fitting using CA atom (default = false)")
      ("saxsFullAtomFit,u", "Fitting using full atom (default = false)")
      ("saxsHeavyAtomFit,h", "Fitting using heavy atom (default = true)") 
      ("genInitial,p",po::value<int>(&geninit_s), "generating simulated initial model (default = false)")
      ("foxs,k", "running foxs program")
      ("pulchar,j", "running pulchar program")		
      ("evaluate,E", "running evaluation mode")
      ("adjust,A", po::value<int>(&lk_s), "setting linker adjust"); 
	
   
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files")
	("min_c1", "running foxs program with min_c1") 
	("max_c1", "running foxs program with max_c1") 
	("min_c2", "running foxs program with min_c2") 
	("max_c2", "running foxs program with max_c2") 
	("pr_dmax", "running foxs program with pr_dmax") ;
	
	po::options_description cmdline_options;
	cmdline_options.add(desc).add(hidden);

	po::options_description visible(
	  "Purpose: Hidden markov model based domain assembly via stepwise, probabilistic sampling\nUsage: prog  -i id -f fasta -s ss -c cmap -m model");
	visible.add(desc);
  
	
	  po::positional_options_description p;
	  p.add("input-files", -1);
	  po::variables_map vm;
	  po::store(po::command_line_parser(argc, argv)
					.options(cmdline_options)
					.positional(p)
					.run(),
				vm);
	  po::notify(vm);
	  
	  std::vector<std::string> files;
	  if (vm.count("input-files")) {
		files = vm["input-files"].as<std::vector<std::string> >();
	  }
	  
	  /** --help option 
       */ 
	
	  //if (vm.count("help") || files.size() == 0) {
	  if (vm.count("help")) {
		std::cout << visible << "\n";
		exit(0);
	  }
	
	if ( vm.count("id") ) 
    { 
		jId = true;
        std::cout << "Setting job name: " 
                << vm["id"].as<std::string>() << std::endl; 
		strcpy(jobId, vm["id"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No job id provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		std::cout << visible << "\n";
		exit(0);
	}
	
	if ( vm.count("fasta") ) 
    { 
		fFile = true;
        std::cout << "Setting fasta file: " 
                << vm["fasta"].as<std::string>() << std::endl; 
		strcpy(faFile, vm["fasta"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No fasta file provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		std::cout << visible << "\n";
		exit(0);
	}
	
	
	if ( vm.count("ss") ) 
    { 
		sFile = true;
        std::cout << "Setting ss file: " 
                << vm["ss"].as<std::string>() << std::endl; 
		strcpy(ssFile, vm["ss"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No secondary structure file provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		std::cout << visible << "\n";
		exit(0);
	}
	
	if ( vm.count("aln") ) 
    { 
		aFile = true;
        std::cout << "Setting alignment file: " 
                << vm["aln"].as<std::string>() << std::endl; 
		strcpy(alnFile, vm["aln"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No alignment file provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		//std::cout << visible << "\n";
		//exit(0);
	}
	
	if ( vm.count("cm") ) 
    { 
		cFile = true;
        std::cout << "Setting alignment file: " 
                << vm["cm"].as<std::string>() << std::endl; 
		strcpy(cmFile, vm["cm"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No contact map file provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		std::cout << visible << "\n";
		exit(0);
	}
	if ( vm.count("domainlist") ) 
    { 
		dFile = true;
        std::cout << "Setting domain list file: " 
                << vm["domainlist"].as<std::string>() << std::endl; 
		strcpy(dmFile, vm["domainlist"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No domain list file provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		std::cout << visible << "\n";
		exit(0);
	}
	
	if ( vm.count("model") ) 
    { 
		mFile = true;
        std::cout << "Setting alignment file: " 
                << vm["model"].as<std::string>() << std::endl; 
		strcpy(moFile, vm["model"].as<std::string>().c_str());
    }else{
		cout << endl;
		cout << "Error! No model file provided" << endl << endl;
		
        //printHelpInfo(argc,argv);
		std::cout << visible << "\n";
		exit(0);
	}
	
	if ( vm.count("numDecoys") ) 
    { 
        std::cout << "Setting number of decoys: " 
                << vm["numDecoys"].as<int>()  << std::endl;  
		numDecoys = vm["numDecoys"].as<int>() ;
    }
	if ( vm.count("numCycles") ) 
    { 
        std::cout << "Setting monte carlo cycle factor: " 
                << vm["numCycles"].as<int>() << std::endl; 
		numCycles = numCycles * vm["numCycles"].as<int>();
    }
	if ( vm.count("adjust") ) 
    { 
        std::cout << "Setting linker adjustment: " 
                << vm["adjust"].as<int>() << std::endl; 
		adjustLinker = vm["adjust"].as<int>();
    }
	if ( vm.count("native") ) 
    { 
        nFile = true;
        std::cout << "Setting native structure: " 
                << vm["native"].as<std::string>() << std::endl; 
		strcpy(natFile, vm["native"].as<std::string>().c_str());
    }
	
	if ( vm.count("init") ) 
    { 
        zFile = true;
        std::cout << "Setting initial structure: " 
                << vm["init"].as<std::string>() << std::endl; 
		strcpy(zeroFile, vm["init"].as<std::string>().c_str());
    }
	if ( vm.count("saxs") ) 
    { 
        eFile = true;
        std::cout << "Setting SAXS experimental profile: " 
                << vm["saxs"].as<std::string>() << std::endl; 
				
		std::string pdb_subfix (".pdb");
		if (vm["saxs"].as<std::string>().find(pdb_subfix) != std::string::npos) 
		{
			std::cout << "SAXS experimental profile file has <.pdb> substring, better to remove!"  << std::endl; 
			exit(0);
		}
		
		strcpy(saxsFile, vm["saxs"].as<std::string>().c_str());
    }
	if ( vm.count("outputdir") ) 
    { 
        odir = true;
        std::cout << "Setting output directory: " 
                << vm["outputdir"].as<std::string>() << std::endl; 
				
		
		strcpy(outputdir, vm["outputdir"].as<std::string>().c_str());

		if (stat(vm["outputdir"].as<std::string>().c_str(), &st) == -1) {
			mkdir(vm["outputdir"].as<std::string>().c_str(), 0700);
		}
    }
	if ( vm.count("wsaxs") ) 
    { 
        std::cout << "Setting SAXS weight: " << vm["wsaxs"].as<double>() << std::endl; 
				
		w_saxs_chi_initial = vm["wsaxs"].as<double>();
    }
	if ( vm.count("scoreFun") ) 
    { 
        
		std::string score1 ("chi");
		std::string score2 ("KL");		
		std::string score3 ("score2");		
		std::string score4 ("RG");		
		
		if(vm["scoreFun"].as<std::string>().compare(score1) ==0 )
		{
			std::cout << "Setting SAXS scorefunction: chi-score" << std::endl; 
			scoreFunction = 1;
		}else if(vm["scoreFun"].as<std::string>().compare(score2) ==0 )
		{
			std::cout << "Setting SAXS scorefunction: KL" << std::endl; 
			scoreFunction = 2;
		}else if(vm["scoreFun"].as<std::string>().compare(score3) ==0 )
		{
			std::cout << "Setting SAXS scorefunction: score2" << std::endl;
			scoreFunction = 3;
		}else if(vm["scoreFun"].as<std::string>().compare(score4) ==0 )
		{
			std::cout << "Setting SAXS scorefunction: RG" << std::endl;
			scoreFunction = 4;
		}else{
			std::cout << "Only 4 functions are supported (chi, KL, RG, score2) " << std::endl; 
			exit(0);
		}
    }
	
	if ( vm.count("scoreWeight") ) 
    { 
        string line=vm["scoreWeight"].as<std::string>();
		double arr[4];
		int i = 0;
		string item;
		stringstream split(line);
		while (std::getline(split, item, '_') && i < 4){
			arr[i]=atof(item.c_str());
			++i;
		}
		for(i = 0; i < 4; i++){
			cout << arr[i] << endl;
		}
	
		
		 w_saxs_chi = arr[0];
		 w_saxs_KL = arr[1];
		 w_saxs_score2 = arr[2];
		 w_saxs_RG_normalize = arr[3];
		 
		 cout << "Setting saxs chi weight to "<< w_saxs_chi<< endl;
		 cout << "Setting saxs KL weight to "<< w_saxs_KL<< endl;
		 cout << "Setting saxs score2 weight to "<< w_saxs_score2<< endl;
		 cout << "Setting saxs RG weight to "<< w_saxs_RG_normalize<< endl;

    }
	
	if ( vm.count("scoreWeightInitial") ) 
    { 
        string line=vm["scoreWeightInitial"].as<std::string>();
		double arr[4];
		int i = 0;
		string item;
		stringstream split(line);
		while (std::getline(split, item, '_') && i < 4){
			arr[i]=atof(item.c_str());
			++i;
		}
		for(i = 0; i < 4; i++){
			cout << arr[i] << endl;
		}
		
		
		 w_saxs_chi_initial = arr[0];
		 w_saxs_KL_initial = arr[1];
		 w_saxs_score2_initial = arr[2];
		 w_saxs_RG_normalize_initial = arr[3];
		 
		 cout << "Setting saxs chi initial weight to "<< w_saxs_chi_initial<< endl;
		 cout << "Setting saxs KL initial weight to "<< w_saxs_KL_initial<< endl;
		 cout << "Setting saxs score2 initial weight to "<< w_saxs_score2_initial<< endl;
		 cout << "Setting saxs RG initial weight to "<< w_saxs_RG_normalize_initial<< endl;

    }
	
	if ( vm.count("scoreWeightFinal") ) 
    { 
        string line=vm["scoreWeightFinal"].as<std::string>();
		double arr[4];
		int i = 0;
		string item;
		stringstream split(line);
		while (std::getline(split, item, '_') && i < 4){
			arr[i]=atof(item.c_str());
			++i;
		}
		for(i = 0; i < 4; i++){
			cout << arr[i] << endl;
		}
	
		
		 w_saxs_chi_final = arr[0];
		 w_saxs_KL_final = arr[1];
		 w_saxs_score2_final = arr[2];
		 w_saxs_RG_normalize_final = arr[3];
		 
		 cout << "Setting saxs chi final weight to "<< w_saxs_chi_final<< endl;
		 cout << "Setting saxs KL final weight to "<< w_saxs_KL_final<< endl;
		 cout << "Setting saxs score2 final weight to "<< w_saxs_score2_final<< endl;
		 cout << "Setting saxs RG final weight to "<< w_saxs_RG_normalize_final<< endl;

    }
	if ( vm.count("scoreWeightPenaltyInitial") ) 
    { 
        string line=vm["scoreWeightPenaltyInitial"].as<std::string>();
		double arr[4];
		int i = 0;
		string item;
		stringstream split(line);
		while (std::getline(split, item, '_') && i < 4){
			arr[i]=atof(item.c_str());
			++i;
		}
		for(i = 0; i < 4; i++){
			cout << arr[i] << endl;
		}
	
		
		 w_saxs_chi_penalty_initial = arr[0];  
		 w_saxs_KL_penalty_initial = arr[1];  
		 w_saxs_score2_penalty_initial = arr[2];  
		 w_saxs_RG_normalize_penalty_initial = arr[3]; 
		 
		 cout << "Setting saxs chi penalty_initial weight to "<< w_saxs_chi_penalty_initial<< endl;
		 cout << "Setting saxs KL penalty_initial weight to "<< w_saxs_KL_penalty_initial<< endl;
		 cout << "Setting saxs score2 penalty_initial weight to "<< w_saxs_score2_penalty_initial<< endl;
		 cout << "Setting saxs RG penalty_initial weight to "<< w_saxs_RG_normalize_penalty_initial<< endl;

    }
	
	if ( vm.count("scoreWeightPenaltyFinal") ) 
    { 
        string line=vm["scoreWeightPenaltyFinal"].as<std::string>();
		double arr[4];
		int i = 0;
		string item;
		stringstream split(line);
		while (std::getline(split, item, '_') && i < 4){
			arr[i]=atof(item.c_str());
			++i;
		}
		for(i = 0; i < 4; i++){
			cout << arr[i] << endl;
		}
	
		
		 w_saxs_chi_penalty_final = arr[0]; 
		 w_saxs_KL_penalty_final = arr[1];
		 w_saxs_score2_penalty_final = arr[2];  
		 w_saxs_RG_normalize_penalty_final = arr[3];  
		 
		 cout << "Setting saxs chi penalty_final weight to "<< w_saxs_chi_penalty_final<< endl;
		 cout << "Setting saxs KL penalty_final weight to "<< w_saxs_KL_penalty_final<< endl;
		 cout << "Setting saxs score2 penalty_final weight to "<< w_saxs_score2_penalty_final<< endl;
		 cout << "Setting saxs RG penalty_final weight to "<< w_saxs_RG_normalize_penalty_final<< endl;

    }
	if ( vm.count("scoreCombine") ) 
    { 
        std::cout << "Setting SAXS score combination function: True" << std::endl; 
		scoreCombine = true;
    }
	if ( vm.count("saxsFit") ) 
    { 
        tofitting = true;
        std::cout << "Setting Fitting parameter: True"  << std::endl; 
    }
	if ( vm.count("saxsCaAtomFit") ) 
    { 
        CAfitting = true;
        std::cout << "Setting CA atom Fitting parameter: True"  << std::endl; 
    }
	if ( vm.count("genInitial") ) 
    { 
        genInitial = true;
        std::cout << "Setting Initial model generation parameter (number of decoys): "  << vm["genInitial"].as<int>() << std::endl;  
		genInitial_numDecoys = vm["genInitial"].as<int>();
    }
	if ( vm.count("reg") ) 
    { 
        regularize = true;
        std::cout << "Setting SAXS regularization: True" << std::endl;  
    }
	 
	if ( vm.count("foxs") ) 
    { 
        foxsonly = true;
        std::cout << "Setting to foxs program: True" << std::endl; 
    }
	if ( vm.count("pulchar") ) 
    { 
        pulcharonly = true;
        std::cout << "Setting to pulchar program: True" << std::endl; 
    }
	
	if ( vm.count("testprog") ) 
    { 
        testonly = vm["testprog"].as<std::string>();
        std::cout << "Setting to testing program: " << testonly << std::endl; 
    }
	if ( vm.count("evaluate") ) 
    { 
        evaluateonly = true;
        std::cout << "Setting to evaluation mode: " << std::endl; 
    }
}




/*************************************************************************
 * Name        : parseNextItem_special
 * Purpose     : parse next item in command line argument, only for fox and pulchar
 * Arguments   : int argc, char ** argv, int & i
 * Return Type : void
 *************************************************************************/
void parseNextItem_special(int argc, char ** argv, int & i) {
    if (strncmp(argv[i], "-k", 2) == 0) {
           if (argc < i + 2) {
            cout << endl;
            cout << "Error! Running fox program only (1/0)" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        saxsfoxsonly = atoi(argv[++i]);
		if(saxsfoxsonly == 1)
		{
			foxsonly = true;
		}else if(saxsfoxsonly == 0){
			foxsonly = false;
		}else{
			cout << "Error! Fitting CA atom parameter is required! (1/0)" << endl << endl;
			exit(-1);
		}
		
    }
    else if (strncmp(argv[i], "-j", 2) == 0) {
           if (argc < i + 2) {
            cout << endl;
            cout << "Error! Running pulchar program only (1/0)" << endl << endl;
            printHelpInfo(argc,argv);
            exit(0);
        }
        saxspulcharonly = atoi(argv[++i]);
		if(saxspulcharonly == 1)
		{
			pulcharonly = true;
		}else if(saxspulcharonly == 0){
			pulcharonly = false;
		}else{
			cout << "Error! Fitting CA atom parameter is required! (1/0)" << endl << endl;
			exit(-1);
		}
		
    }
    i++;
}



void parseOption_special(int argc, char ** argv) {
	
	
	std::string jobId_s;
	std::string fasta_s;
	std::string ss_s;
	std::string aln_s;
	std::string cm_s;
	std::string mod_s;
	std::string nat_s;
	std::string init_s;
	std::string saxs_s;
	int de_s=0;
	int cy_s=0;
	
    po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("verbose,v", "print words with verbosity") 
      ("id,i", "job id (required)")
      ("fasta,f", "fasta file (required)")
      ("ss,s",  "secondary structure file (required)")
      ("aln,a",  "alignment file (required)")
      ("cm,c", "contact map file (required)")
      ("domainlist,l",  "domain path list (required)")
      ("model,m",  "IOHMM model file (required)")
	  ("numDecoys,d", "number of decoys")
	  ("numCycles,x", "monte carlo cycle factor")
	  ("native,n",  "native structure file")
	  ("init,z", "initial structure file")
	  ("saxs,e", "SAXS experimental profile file")
	  ("scoreFun,F", "SAXS score function (chi, KL, score2, RG)")
	  ("scoreWeight",  "weight for SAXS score function (chi, KL, score2, RG):  default: 100-100-100-100")
	  ("scoreWeightInitial",  "intial weight for SAXS score function (chi, KL, score2, RG):  default: 100-100-100-100")
	  ("scoreWeightFinal", "Final weight for SAXS score function (chi, KL, score2, RG):  default: 1-1-1-1")
	  ("scoreWeightPenalyInitial", "intial penalty weight for SAXS score function (chi, KL, score2, RG):  default: 0.001-0.001-0.001-0.001")
	  ("scoreCombine,B", "SAXS score function combination (chi, RG, score2, KL)")
	  ("wsaxs,w", "SAXS weight (default=50)")
	  ("testprog,g",  "test the program")
      ("reg,r", "perform saxs regularization fitting (default = false)")
      ("saxsFit,t", "perform saxs fitting (default = false)")
      ("saxsCaAtomFit,b", "Fitting using CA atom (default = false)")
      ("saxsFullAtomFit,u", "Fitting using full atom (default = false)")
      ("saxsHeavyAtomFit,h", "Fitting using heavy atom (default = true)")
      ("genInitial,p", "generating simulated initial model (default = false)")
      ("foxs,k", "running foxs program")
      ("pulchar,j", "running pulchar program")	
      ("evaluate,E", "running evaluation mode");
	
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files")
	("min_c1", "running foxs program with min_c1") 
	("max_c1", "running foxs program with max_c1") 
	("min_c2", "running foxs program with min_c2") 
	("max_c2", "running foxs program with max_c2") 
	("pr_dmax", "running foxs program with pr_dmax") ;
	
	po::options_description cmdline_options;
	cmdline_options.add(desc).add(hidden);

	po::options_description visible(
	  "Purpose: Hidden markov model based domain assembly via stepwise, probabilistic sampling\nUsage: prog  -i id -f fasta -s ss -c cmap -m model");
	visible.add(desc);
  
	
	  po::positional_options_description p;
	  p.add("input-files", -1);
	  po::variables_map vm;
	  po::store(po::command_line_parser(argc, argv)
					.options(cmdline_options)
					.positional(p)
					.run(),
				vm);
	  po::notify(vm);
	  
	  std::vector<std::string> files;
	  if (vm.count("input-files")) {
		files = vm["input-files"].as<std::vector<std::string> >();
	  }
	  
	  /** --help option 
       */ 
	
	  //if (vm.count("help") || files.size() == 0) {
	  if (vm.count("help")) {
		std::cout << visible << "\n";
		exit(0);
	  }
	
	if ( vm.count("foxs") ) 
    { 
        foxsonly = true;
        std::cout << "Setting to foxs program: " << std::endl; 
    }
	if ( vm.count("pulchar") ) 
    { 
        pulcharonly = true;
        std::cout << "Setting to pulchar program: " << std::endl; 
    }
	if ( vm.count("evaluate") ) 
    { 
        evaluateonly = true;
        std::cout << "Setting to evaluation mode: " << std::endl; 
    }
	
}
/*************************************************************************
 * Name        : parseCommandLine
 * Purpose     : parse command line arguments
 * Arguments   : int argc, char ** argv
 * Return Type : void
 *************************************************************************/
void parseCommandLine(int argc, char ** argv) {
    int i = 1;
    while (i < argc)
        parseNextItem_special(argc, argv, i);
	//parseOption_special(argc,argv);
	
	
	po::options_description desc("Options"); 
    desc.add_options() 
      ("help,h", "Print help messages") 
      ("verbose,v", "print words with verbosity") 
      ("id,i", "job id (required)")
      ("fasta,f", "fasta file (required)")
      ("ss,s",  "secondary structure file (required)")
      ("aln,a",  "alignment file (required)")
      ("cm,c", "contact map file (required)")
      ("domainlist,l",  "domain path list (required)")
      ("model,m",  "IOHMM model file (required)")
	  ("numDecoys,d", "number of decoys")
	  ("numCycles,x", "monte carlo cycle factor")
	  ("native,n",  "native structure file")
	  ("init,z", "initial structure file")
	  ("saxs,e", "SAXS experimental profile file")
	  ("scoreFun,F", "SAXS score function (chi, RG, score2, KL)")
	  ("scoreCombine,B", "SAXS score function combination (chi, RG, score2, KL)")
	  ("wsaxs,w", "SAXS weight (default=50)")
	  ("testprog,g",  "test the program")
      ("reg,r", "perform saxs regularization fitting (default = false)")
      ("saxsFit,t", "perform saxs fitting (default = false)")
      ("saxsCaAtomFit,b", "Fitting using CA atom (default = false)")
      ("saxsFullAtomFit,u", "Fitting using full atom (default = false)")
      ("saxsHeavyAtomFit,h", "Fitting using heavy atom (default = true)")
      ("genInitial,p", "generating simulated initial model (default = false)")
      ("foxs,k", "running foxs program")
      ("pulchar,j", "running pulchar program")
      ("evaluate,E", "running evaluation mode")
      ("adjust,A", "setting linker adjust"); 	
	
    po::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files");
	
	po::options_description cmdline_options;
	cmdline_options.add(desc).add(hidden);

	po::options_description visible(
	  "Purpose: Hidden markov model based domain assembly via stepwise, probabilistic sampling\nUsage: prog  -i id -f fasta -s ss -c cmap -m model");
	visible.add(desc);
	
		// if run fox only 
	if (pulcharonly && foxsonly)
	{
            cout << endl;
            cout << "Error! <-j> and <-k> can not use the same time!" << endl << endl;
            std::cout << visible << "\n";
            exit(0);
		
	}else if(foxsonly) {
		cout << "Running fox program only" << endl;
		char ** foxsargv;
		int foxsargc=0;
	    foxsargv = new char*[5000];// initialize the double pointer

		for (int i = 0; i < argc; i++) 
		{
			if(strncmp(argv[i], "-k", 2) != 0)
			{
				foxsargv[foxsargc]=new char[1000];
				strcpy(foxsargv[foxsargc], argv[i]);
				foxsargc++;
			}else{
				i++; // ignore the next value  -j 1
			}
		}
		run_foxs_only(foxsargc,foxsargv);
	}else if(pulcharonly) {
		cout << "Running pulchar program only" << endl;
		char ** pulcharargv;
		int pulcharargc=0;
	    pulcharargv = new char*[1000];// initialize the double pointer

		for (int i = 0; i < argc; i++) 
		{
			if(strncmp(argv[i], "-j", 2) != 0)
			{
				pulcharargv[pulcharargc]=new char[1000];
				strcpy(pulcharargv[pulcharargc], argv[i]);
				pulcharargc++;
			}else{
				i++; // ignore the next value  -k 1
			}
		}
		
		run_Pulcha_int(pulcharargc,pulcharargv);
	}else{
	
		//int i = 1;
		//while (i < argc)
		//	parseNextItem(argc, argv, i);
		parseOption(argc,argv);
		
		if(genInitial and zFile)
		{
            cout << endl;
            cout << "Error! <-z> and <-p> can not use the same time for initial model generation!" << endl << endl;
            std::cout << visible << "\n";
            exit(0);
		}
	
		if (!jId) {
			cout << endl;
			cout << "Error! Job id must be provided" << endl << endl;
            std::cout << visible << "\n";
			exit(0);
		}
		if (!fFile) {
			cout << endl;
			cout << "Error! Fasta file must be provided" << endl << endl;
            std::cout << visible << "\n";
			exit(0);
		}
		else {
			ifstream fin( faFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Fasta file not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (!sFile) {
			cout << endl;
			cout << "Error! Secondary structure file must be provided" << endl << endl;
            std::cout << visible << "\n";
			exit(0);
		}
		else {
			ifstream fin( ssFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Secondary structure file not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (!aFile) {  // alignment is not required because if multiple domain is provide, can detect by algorithm
			cout << endl;
			//cout << "Error! Alignment file must be provided" << endl << endl;
            //std::cout << visible << "\n";
			//exit(0);
		}
		else {
			ifstream fin( alnFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Alignment file not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (!cFile) {
			cout << endl;
			cout << "Error! Contact map file must be provided" << endl << endl;
            std::cout << visible << "\n";
			exit(0);
		}
		else {
			ifstream fin( cmFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Contact map file not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (!dFile) {
			cout << endl;
			cout << "Error! Domain list file must be provided" << endl << endl;
            std::cout << visible << "\n";
			exit(0);
		}
		else {
			ifstream fin( dmFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Domain list file not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (!mFile) {
			cout << endl;
			cout << "Error! Model file must be provided" << endl << endl;
            std::cout << visible << "\n";
			exit(0);
		}
		else {
			ifstream fin( moFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Model file not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (nFile) {
			ifstream fin( natFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! Native structure not present" << endl << endl;
				std::cout << visible << "\n";
				exit(0);
			}
		}
		if (zFile) {
			ifstream fin( zeroFile);
			if( fin.fail() ) {
				cout << endl;
				cout << "Error! initial structure not present" << endl << endl;
				//printHelpInfo(argc,argv);
				std::cout << visible << "\n";
				exit(0);
			}
		}
	}
}

/*************************************************************************
 * Name        : getAA
 * Purpose     : convert AA name to a numerical code
 * Arguments   : const char * aa
 * Return Type : int
 *************************************************************************/
int getAA(const char * aa) {
    if (strlen(aa) == 3) {
        if (strcmp(aa, "ALA") == 0)
            return (ALA);
        else if (strcmp(aa, "ARG") == 0)
            return (ARG);
        else if (strcmp(aa, "ASN") == 0)
            return (ASN);
        else if (strcmp(aa, "ASP") == 0)
            return (ASP);
        else if (strcmp(aa, "CYS") == 0)
            return (CYS);
        else if (strcmp(aa, "GLN") == 0)
            return (GLN);
        else if (strcmp(aa, "GLU") == 0)
            return (GLU);
        else if (strcmp(aa, "GLY") == 0)
            return (GLY);
        else if (strcmp(aa, "HIS") == 0)
            return (HIS);
        else if (strcmp(aa, "ILE") == 0)
            return (ILE);
        else if (strcmp(aa, "LEU") == 0)
            return (LEU);
        else if (strcmp(aa, "LYS") == 0)
            return (LYS);
        else if (strcmp(aa, "MET") == 0)
            return (MET);
        else if (strcmp(aa, "PHE") == 0)
            return (PHE);
        else if (strcmp(aa, "PRO") == 0)
            return (PRO);
        else if (strcmp(aa, "SER") == 0)
            return (SER);
        else if (strcmp(aa, "THR") == 0)
            return (THR);
        else if (strcmp(aa, "TRP") == 0)
            return (TRP);
        else if (strcmp(aa, "TYR") == 0)
            return (TYR);
        else if (strcmp(aa, "VAL") == 0)
            return (VAL);
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else if (strlen(aa) == 1) {
        if (aa[0] == 'A')
            return ALA;
        else if (aa[0] == 'C')
            return CYS;
        else if (aa[0] == 'D')
            return ASP;
        else if (aa[0] == 'E')
            return GLU;
        else if (aa[0] == 'F')
            return PHE;
        else if (aa[0] == 'G')
            return GLY;
        else if (aa[0] == 'H')
            return HIS;
        else if (aa[0] == 'I')
            return ILE;
        else if (aa[0] == 'K')
            return LYS;
        else if (aa[0] == 'L')
            return LEU;
        else if (aa[0] == 'M')
            return MET;
        else if (aa[0] == 'N')
            return ASN;
        else if (aa[0] == 'P')
            return PRO;
        else if (aa[0] == 'Q')
            return GLN;
        else if (aa[0] == 'R')
            return ARG;
        else if (aa[0] == 'S')
            return SER;
        else if (aa[0] == 'T')
            return THR;
        else if (aa[0] == 'V')
            return VAL;
        else if (aa[0] == 'W')
            return TRP;
        else if (aa[0] == 'Y')
            return TYR;
        else {
            cout << "Error! Invalid amino acid " << aa << endl;
            exit(0);
        }
    }
    else {
        cout << "Error! Invalid amino acid " << aa << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : getSS
 * Purpose     : convert SS name to a numerical code
 * Arguments   : const char * ss
 * Return Type : int
 *************************************************************************/
int getSS(const char * ss) {
    if (strlen(ss) == 1) {
        if (strcmp(ss, "H") == 0)
            return (ALPHA_HELIX);
        else if (strcmp(ss, "G") == 0)
            return (THREE_HELIX);
        else if (strcmp(ss, "I") == 0)
            return (FIVE_HELIX);
        else if (strcmp(ss, "B") == 0)
            return (ISOLATED_STRAND);
        else if (strcmp(ss, "E") == 0)
            return (EXTENDED_STRAND);
        else if (strcmp(ss, "T") == 0)
            return (HBONDED_TURN);
        else if (strcmp(ss, "S") == 0)
            return (NONHBONDED_BEND);
        else if (strcmp(ss, "C") == 0)
            return (RANDOM_COIL);
        else {
            cout << "Error! Invalid secondary structure " << ss << endl;
            exit(0);
        }
    }
    else {
        cout << "Error! Invalid secondary structure " << ss << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadPdb
 * Purpose     : loads a pdb file into a vector of pdbInfo object
 * Arguments   : char *filename, vector<pdbInfo> &pdb
 * Return Type : void
 *************************************************************************/
void loadPdb(char *filename, vector<pdbInfo> &pdb) {
    string line, str;
    string atom ("ATOM ");
    int prevRes = -999999;
    point3d caAtom;
    point3d scAtom;
    vector<point3d> scAtomCloud;
    pdbInfo pdbData;
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            if(line.compare(0, atom.length(), atom)==0) {
                int res = atoi(line.substr(22, 4).c_str());
                int anmae = getAA(line.substr(17, 3).c_str());
                // seek for the next residue
                if (res != prevRes) {
                    if (prevRes != -999999) {
                        // inser the data collected so far
                        pdbData.ca = caAtom;
                        pdbData.sc = getCentroid(scAtomCloud);
                        pdb.push_back(pdbData);
                    }
                    prevRes = res;
                    scAtomCloud.clear();
                    pdbData.id = res;
                    pdbData.aa = anmae;
                }
                // consider the first alternate location id (i.e. A) if present
                if( line.compare(16, 1, " ") == 0 || line.compare(16, 1, "A") == 0 ) {
                    // get the CA atom coordinate
                    if( line.compare(12, 4, "CA  ") == 0 || line.compare(12, 4, " CA ") == 0 || line.compare(12, 4, "  CA") == 0 ) {
                        
                        caAtom.x = atof(line.substr(30, 8).c_str());
                        caAtom.y = atof(line.substr(38, 8).c_str());
                        caAtom.z = atof(line.substr(46, 8).c_str());
                        scAtomCloud.push_back(caAtom);
                    }
                    //cout << "line "<< line <<" length: "<<line.length() << endl;
                    if(line.length()<77)// sometimes native structure don't have label at position 77
                    {
                      // get the sidechain heavy atoms
                      if(line.compare(12, 4, "N   ") != 0 && line.compare(12, 4, " N  ") != 0 && line.compare(12, 4, "  N ") != 0 &&
                         line.compare(12, 4, "C   ") != 0 && line.compare(12, 4, " C  ") != 0 && line.compare(12, 4, "  C ") != 0 &&
                         line.compare(12, 4, "O   ") != 0 && line.compare(12, 4, " O  ") != 0 && line.compare(12, 4, "  O ") != 0 &&
                         line.compare(12, 4, "CA  ") != 0 && line.compare(12, 4, " CA ") != 0 && line.compare(12, 4, "  CA") != 0) { 
                          
                          scAtom.x = atof(line.substr(30, 8).c_str());
                          scAtom.y = atof(line.substr(38, 8).c_str());
                          scAtom.z = atof(line.substr(46, 8).c_str());
                          scAtomCloud.push_back(scAtom);
                      }
                    
                    }else{
                      // get the sidechain heavy atoms
                      if(line.compare(12, 4, "N   ") != 0 && line.compare(12, 4, " N  ") != 0 && line.compare(12, 4, "  N ") != 0 &&
                         line.compare(12, 4, "C   ") != 0 && line.compare(12, 4, " C  ") != 0 && line.compare(12, 4, "  C ") != 0 &&
                         line.compare(12, 4, "O   ") != 0 && line.compare(12, 4, " O  ") != 0 && line.compare(12, 4, "  O ") != 0 &&
                         line.compare(12, 4, "CA  ") != 0 && line.compare(12, 4, " CA ") != 0 && line.compare(12, 4, "  CA") != 0 &&
                         line.compare(77, 1, "H") != 0) { 
                          
                          scAtom.x = atof(line.substr(30, 8).c_str());
                          scAtom.y = atof(line.substr(38, 8).c_str());
                          scAtom.z = atof(line.substr(46, 8).c_str());
                          scAtomCloud.push_back(scAtom);
                      }
                    }
                }
            }
        }
        fin.close();
        // for the last residue
        pdbData.ca = caAtom;
        pdbData.sc = getCentroid(scAtomCloud);
        pdb.push_back(pdbData);
    }
    else {
        cout << "Error! pdb file can not open " << filename << endl;
        exit(0);
    }
    
}

/*************************************************************************
 * Name        : loadFasta
 * Purpose     : loads a fasta file into a vector of int object
 * Arguments   : char *filename, vector<int> &aa
 * Return Type : void
 *************************************************************************/
void loadFasta(char *filename, vector<int> &aa) {
    string line, str;
    string header (">");
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            if(line.length() != 0 && line.compare(0, header.length(), header) != 0) {
                for (int i = 0; i < line.length(); i++) {
                    int aa_code = getAA(line.substr(i,1).c_str());
                    aa.push_back(aa_code);
                }
            }
        }
        fin.close();
    }
    else {
        cout << "Error! fasta file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadSec
 * Purpose     : loads a secondary structure file into a vector of int object
 * Arguments   : char *filename, vector<int> &ss
 * Return Type : void
 *************************************************************************/
void loadSec(char *filename, vector<int> &ss) {
    string line, str;
    string header (">");
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            if(line.length() != 0 && line.compare(0, header.length(), header) != 0) {
                for (int i = 0; i < line.length(); i++) {
                    int ss_code = getSS(line.substr(i,1).c_str());
                    ss.push_back(ss_code);
                }
            }
        }
        fin.close();
    }
    else {
        cout << "Error! secondary structure file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadAlign
 * Purpose     : loads a Alignment Code into a vector of int object
 * Arguments   : char *filename, vector<int> &alnCode
 * Return Type : void
 *************************************************************************/
void loadAlign(char *filename, vector<int> &alnCode) {
    string line, str;
    string header (">");
    ifstream fin (filename);
    if (fin.is_open()) {
        while ( fin.good() ) {
            getline(fin, line);
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            if(line.length() != 0 && line.compare(0, header.length(), header) != 0) {
                for (int i = 0; i < line.length(); i++) {
                    int aln_code =  atoi(line.substr(i,1).c_str()); //  1: align or 0:unalign
                    alnCode.push_back(aln_code);
                }
            }
        }
        fin.close();
    }
    else {
        cout << "Error! alignment file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : loadCmap
 * Purpose     : loads a contact map into a MDArray object
 * Arguments   : char *filename, vector<int> &sa
 * Return Type : void
 *************************************************************************/
void loadCmap(char *filename, MDArray<double> &cm) {
    ifstream fin (filename);
    double val;
    if (fin.is_open()) {
        for (int i = 0; i < cm.get_shape()[0]; i++) {
            for (int j = 0; j < cm.get_shape()[0]; j++) {
                fin >> val;
                cm.set(i, j, val);
            }
        }
        fin.close();
    }
    else {
        cout << "Error! contact map file can not open " << filename << endl;
        exit(0);
    }
}

/*************************************************************************
 * Name        : getDistance
 * Purpose     : gets distance between two points
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : double
 *************************************************************************/
double getDistance(point3d & p1, point3d &p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

/*************************************************************************
 * Name        : getDotProduct
 * Purpose     : gets the dot product (i.e. scalar product) for two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : double
 *************************************************************************/
double getDotProduct(point3d & p1, point3d &p2) {
    return (p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
}

/*************************************************************************
 * Name        : getCrossProduct
 * Purpose     : gets the cross product (i.e. vector product) for two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getCrossProduct(point3d & p1, point3d &p2) {
    point3d p;
    p.x = p1.y * p2.z - p1.z * p2.y;
    p.y = p1.z * p2.x - p1.x * p2.z;
    p.z = p1.x * p2.y - p1.y * p2.x;
    return p;
}

/*************************************************************************
 * Name        : getNorm
 * Purpose     : gets the norm (i.e. length) of a vector
 * Arguments   : point3d & p
 * Return Type : double
 *************************************************************************/
double getNorm(point3d & p) {
    return sqrt( p.x * p.x + p.y * p.y + p.z * p.z);
}

/*************************************************************************
 * Name        : getDifference
 * Purpose     : gets the difference vector between two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getDifference(point3d & p1, point3d &p2) {
    point3d p;
    p.x = p2.x - p1.x;
    p.y = p2.y - p1.y;
    p.z = p2.z - p1.z;
    return p;
}

/*************************************************************************
 * Name        : getMidpoint
 * Purpose     : gets the midpoint vector between two vectors
 * Arguments   : point3d & p1, point3d &p2
 * Return Type : point3d
 *************************************************************************/
point3d getMidpoint(point3d & p1, point3d &p2) {
    point3d p;
    p.x = (p2.x + p1.x) / 2.0;
    p.y = (p2.y + p1.y) / 2.0;
    p.z = (p2.z + p1.z) / 2.0;
    return p;
}

/*************************************************************************
 * Name        : getUnit
 * Purpose     : gets the unit vector for a vector
 * Arguments   : point3d & p
 * Return Type : point3d
 *************************************************************************/
point3d getUnit(point3d & p) {
    point3d u;
    double norm = getNorm(p);
    if (norm == 0.0) {
        u.x = 0.0;
        u.y = 0.0;
        u.z = 0.0;
    }
    else {
        u.x = p.x / norm;
        u.y = p.y / norm;
        u.z = p.z / norm;
    }
    return u;
}

/*************************************************************************
 * Name        : getAngle
 * Purpose     : gets the angle between three points
 * Arguments   : point3d & p1, point3d &p2, point3d &p3
 * Return Type : double
 *************************************************************************/
double getAngle(point3d & p1, point3d &p2, point3d &p3) {
    double acc = 0.0;
    acc = (p2.x - p1.x) * (p2.x - p3.x) + (p2.y - p1.y) * (p2.y - p3.y) + (p2.z - p1.z) * (p2.z - p3.z);
    double d1 = getDistance(p1, p2);
    double d2 = getDistance(p3, p2);
    if (d1 == 0 || d2 == 0) {
        return 0;
    }
    acc = acc / (d1 * d2);
    if (acc > 1.0)
        acc = 1.0;
    else if (acc < -1.0)
        acc = -1.0;
    acc = acos(acc);
    return acc;
}

/*************************************************************************
 * Name        : getDihedral
 * Purpose     : gets the dihedral between four points
 * Arguments   : point3d & p1, point3d &p2, point3d &p3, point3d &p4
 * Return Type : double
 *************************************************************************/
double getDihedral(point3d & p1, point3d &p2, point3d &p3, point3d &p4) {
    point3d q, r, s, t, u, v, z;
    double acc, w;
    z.x = z.y = z.z = 0.0;
    q = getDifference(p1, p2);
    r = getDifference(p3, p2);
    s = getDifference(p4, p3);
    t = getCrossProduct(q, r);
    u = getCrossProduct(s, r);
    v = getCrossProduct(u, t);
    w = getDotProduct(v, r);
    acc = getAngle(t, z, u);
    if (w < 0)
        acc = -acc;
    return (acc);
}

/*************************************************************************
 * Name        : getCentroid
 * Purpose     : gets centroid of a cloud of points
 * Arguments   : vector<point3d> &pointCloud
 * Return Type : point3d
 *************************************************************************/
point3d getCentroid(vector<point3d> &pointCloud) {
    point3d centroid;
    centroid.x = 0.0;
    centroid.y = 0.0;
    centroid.z = 0.0;
    for (int i = 0; i < pointCloud.size(); i++) {
        centroid.x += pointCloud[i].x / pointCloud.size();
        centroid.y += pointCloud[i].y / pointCloud.size();
        centroid.z += pointCloud[i].z / pointCloud.size();
    }
    return centroid;
}

/*************************************************************************
 * Name        : setCoordinate
 * Purpose     : set coordinate of a point based on dihedral alpha, bond angle tau, bond length normw
 * Arguments   : point3d &c0, point3d &c1, point3d &c2, point3d &c3, double alpha, double tau, double normw
 * Return Type : void
 *************************************************************************/
void setCoordinate(point3d &c0, point3d &c1, point3d &c2, point3d &c3, double alpha, double tau, double normw) {
    double	u1, u2, u3, v1, v2, v3, norm, pvuv1, pvuv2, pvuv3, pvvuv1, pvvuv2, pvvuv3, nsa, nca, nct;
    u1 = (c1.x - c0.x);
    u2 = (c1.y - c0.y);
    u3 = (c1.z - c0.z);
    v1 = (c2.x - c1.x);
    v2 = (c2.y - c1.y);
    v3 = (c2.z - c1.z);
    norm = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
    v1 /=  norm;
    v2 /=  norm;
    v3 /=  norm;
    pvuv1 = u2 * v3 - u3 * v2;
    pvuv2 = u3 * v1 - u1 * v3;
    pvuv3 = u1 * v2 - u2 * v1;
    norm = sqrt(pvuv1 * pvuv1 + pvuv2 * pvuv2 + pvuv3 * pvuv3);
    pvuv1 /= norm;
    pvuv2 /= norm;
    pvuv3 /= norm;
    pvvuv1 = v2 * pvuv3 - v3 * pvuv2;
    pvvuv2 = v3 * pvuv1 - v1 * pvuv3;
    pvvuv3 = v1 * pvuv2 - v2 * pvuv1;
    norm = sqrt(pvvuv1 * pvvuv1 + pvvuv2 * pvvuv2 + pvvuv3 * pvvuv3);
    pvvuv1 /= norm;
    pvvuv2 /= norm;
    pvvuv3 /= norm;
    nca = cos(alpha);
    nsa = sin(alpha);
    nct = tan(tau - PI/2);
    u1 = nca * (-pvvuv1) + nsa * pvuv1 + v1 * nct;
    u2 = nca * (-pvvuv2) + nsa * pvuv2 + v2 * nct;
    u3 = nca * (-pvvuv3) + nsa * pvuv3 + v3 * nct;
    norm = sqrt(u1 * u1 + u2 * u2 + u3 * u3);
    u1 = u1 * normw/norm;
    u2 = u2 * normw/norm;
    u3 = u3 * normw/norm;
    c3.x = u1 + c2.x;
    c3.y = u2 + c2.y;
    c3.z = u3 + c2.z;
}

/*************************************************************************
 * Name        : pdb2pose
 * Purpose     : converts pdb to pose
 * Arguments   : vector<pdbInfo> &pdb, vector<int> &ss, vector<poseInfo> &pose
 * Return Type : void
 *************************************************************************/
void pdb2pose(vector<pdbInfo> &pdb, vector<int> &ss, vector<poseInfo> &pose) {
    pose.clear();
    poseInfo poseData;
    //define a virtual N terminus
    point3d N;
    N.x = pdb[0].ca.x - cos(PI) * CA2CA;
    N.y = pdb[0].ca.y - sin(PI) * CA2CA;
    N.z = pdb[0].ca.z;
    //define a virtual C terminus
    point3d C;
    C.x = pdb[pdb.size()-1].ca.x + cos(PI) * CA2CA;
    C.y = pdb[pdb.size()-1].ca.y + sin(PI) * CA2CA;
    C.z = pdb[pdb.size()-1].ca.z;
    for (int i = 0; i < pdb.size(); i++) {
        poseData.id = pdb[i].id;
        poseData.aa = pdb[i].aa;
        poseData.ss = ss[i];
        // calculate bca
        if (i == 0) {
            poseData.bca = CA2CA;
        }
        else {
            poseData.bca = getDistance(pdb[i-1].ca, pdb[i].ca);
        }
        // calculate tao
        if (i == 0 || i == pdb.size() - 1 || i == pdb.size() - 2) {
            poseData.tao = 2 * PI;
        }
        else {
            poseData.tao = getDihedral(pdb[i-1].ca, pdb[i].ca, pdb[i+1].ca, pdb[i+2].ca);
        }
        // calculate theta
        if (i == 0 || i == pdb.size() - 1) {
            poseData.theta = 2 * PI;
        }
        else {
            poseData.theta = getAngle(pdb[i-1].ca, pdb[i].ca, pdb[i+1].ca);
        }
        //calculate bsc
        poseData.bsc = getDistance(pdb[i].ca, pdb[i].sc);
        //calculate phi
        if (i == 0) {
            poseData.phi = getDihedral(pdb[i+1].ca, N, pdb[i].ca, pdb[i].sc);
        }
        else if (i == pdb.size() - 1) {
            poseData.phi = getDihedral(C, pdb[i-1].ca, pdb[i].ca, pdb[i].sc);
        }
        else {
            poseData.phi = getDihedral(pdb[i+1].ca, pdb[i-1].ca, pdb[i].ca, pdb[i].sc);
        }
        //calculate delta
        if (i == 0) {
            poseData.delta = getAngle(N, pdb[i].ca, pdb[i].sc);
        }
        else {
            poseData.delta = getAngle(pdb[i-1].ca,  pdb[i].ca, pdb[i].sc);
        }
        // populate the pose
        pose.push_back(poseData);
    }
    
}

/*************************************************************************
 * Name        : pose2pdb
 * Purpose     : converts pose to pdb
 * Arguments   : vector<poseInfo> &pose, vector<pdbInfo> &pdb
 * Return Type : void
 *************************************************************************/
void pose2pdb(vector<poseInfo> &pose, vector<pdbInfo> &pdb) {
    pdb.clear();
    pdbInfo pdbData;
    //define a virtual N terminus
    point3d N;
    N.x = 0.0 - cos(PI) * CA2CA;
    N.y = 0.0 - sin(PI) * CA2CA;
    N.z = 0.0;
    // set the first residue's CA atom at the origin
    pdbData.id = pose[0].id;
    pdbData.aa = pose[0].aa;
    pdbData.ca.x = 0.0;
    pdbData.ca.y = 0.0;
    pdbData.ca.z = 0.0;
    pdbData.sc = pdbData.ca;
    pdb.push_back(pdbData);
    // set the second residue's CA atom by bond angle and bond length
    pdbData.id = pose[1].id;
    pdbData.aa = pose[1].aa;
    pdbData.ca.x = pose[1].bca;
    pdbData.ca.y = 0.0;
    pdbData.ca.z = 0.0;
    pdbData.sc = pdbData.ca;
    pdb.push_back(pdbData);
    // set the third residue's CA atom by dihedral, bond angle and bond length
    pdbData.id = pose[2].id;
    pdbData.aa = pose[2].aa;
    setCoordinate(N, pdb[0].ca, pdb[1].ca, pdbData.ca, pose[1].tao, pose[1].theta, pose[1].bca);
    pdbData.sc = pdbData.ca;
    pdb.push_back(pdbData);
    // set the rest of the residues' CA atoms by dihedral, bond angle and bond length
    for (int i = 3; i < pose.size(); i++) {
        pdbData.id = pose[i].id;
        pdbData.aa = pose[i].aa;
        setCoordinate(pdb[i-3].ca, pdb[i-2].ca, pdb[i-1].ca, pdbData.ca, pose[i-2].tao, pose[i-1].theta, pose[i-1].bca);
        pdbData.sc = pdbData.ca;
        pdb.push_back(pdbData);
    }
    // define a virtual C terminus
    point3d C;
    C.x = pdb[pdb.size()-1].ca.x + cos(PI) * CA2CA;
    C.y = pdb[pdb.size()-1].ca.y + sin(PI) * CA2CA;
    C.z = pdb[pdb.size()-1].ca.z;
    // set all residues' SC atoms based on CA atoms by dihedral, bond angle and bond length
    for (int i = 0; i < pose.size(); i++) {
        if (pose[i].aa != GLY) {
            if (i == 0) {
                setCoordinate(pdb[i+1].ca, N, pdb[i].ca, pdb[i].sc, pose[i].phi, pose[i].delta, pose[i].bsc);
            }
            else if (i == pose.size() - 1) {
                setCoordinate(C, pdb[i-1].ca, pdb[i].ca, pdb[i].sc, pose[i].phi, pose[i].delta, pose[i].bsc);
            }
            else {
                setCoordinate(pdb[i+1].ca, pdb[i-1].ca, pdb[i].ca, pdb[i].sc, pose[i].phi, pose[i].delta, pose[i].bsc);
            }
        }
    }
}

/*************************************************************************
 * Name        : sample2pose
 * Purpose     : converts sample to pose
 * Arguments   : MDArray<double> &sample, vector<poseInfo> &pose
 * Return Type : void
 *************************************************************************/
void sample2pose(MDArray<double> &sample, vector<poseInfo> &pose) {
    pose.clear();
    for (int i = 0; i < sample.get_shape()[0] / 2; i++) {
        poseInfo poseData;
        poseData.id = i + 1;
        poseData.aa = (int)sample.get(2 * i, 2);
        poseData.ss = (int)sample.get(2 * i, 3);
        poseData.bca = sample.get(2 * i, 4);
        poseData.tao = sample.get(2 * i, 5);
        poseData.theta = sample.get(2 * i, 6);
        poseData.bsc = sample.get(2 * i + 1, 4);
        poseData.phi = sample.get(2 * i + 1, 5);
        poseData.delta = sample.get(2 * i + 1, 6);
        pose.push_back(poseData);
    }
}

/*************************************************************************
 * Name        : pose2sample
 * Purpose     : converts pose to sample
 * Arguments   : vector<poseInfo> &pose, MDArray<double> &sample
 * Return Type : void
 *************************************************************************/
void pose2sample(vector<poseInfo> &pose, MDArray<double> &sample) {
    sample.set_shape(pose.size() * 2, NUM_DAT);
    for (int i = 0; i < pose.size(); i++) {
        sample.set(2 * i, 0, 0);
        sample.set(2 * i, 1, 0);
        sample.set(2 * i, 2, pose[i].aa);
        sample.set(2 * i, 3, pose[i].ss);
        sample.set(2 * i, 4, pose[i].bca);
        sample.set(2 * i, 5, pose[i].tao);
        sample.set(2 * i, 6, pose[i].theta); // deb's original version has problem

        sample.set(2 * i + 1, 0, 1);
        sample.set(2 * i + 1, 1, 0);
        sample.set(2 * i + 1, 2, pose[i].aa);
        sample.set(2 * i + 1, 3, pose[i].ss);
        sample.set(2 * i + 1, 4, pose[i].bsc);
        sample.set(2 * i + 1, 5, pose[i].phi);
        sample.set(2 * i + 1, 6, pose[i].delta);
   }
}



/*************************************************************************
 * Name        : gen_initial_model
 * Purpose     : converts pose to pdb
 * Arguments   : vector<int> &aa,vector<int> &ss, DBN &dbn, vector<pdbInfo> &pdbFoldon
 * Return Type : void
 *************************************************************************/
void gen_initial_model(vector<int> &aa,vector<int> &ss, MDArray<double> &cm,vector<pdbInfo> &native,DBN &dbn, vector<pdbInfo> &pdbFoldon) {
		
    // extract secondary structural elements
    vector<sseInfo> sse;
    bool sseStarted = false;
    sseInfo sseData;
    for (int i = 0; i < ss.size() - 1; i++) {
        if (!sseStarted) {
            sseData.start = i;
            sseData.type = ss[i];
            sseStarted = true;
        }
        if (sseStarted && ss[i] != ss[i+1]) {
            sseData.end = i;
            sse.push_back(sseData);
            sseStarted = false;
        }
    }
    // for the last sse
    sseData.type = ss[ss.size() - 1];
    sseData.end = ss.size() - 1;
    sse.push_back(sseData);
    
    
    // extract foldon units
    vector<sseInfo> foldon;
    bool foldonStarted = false;
    sseInfo foldonData;
	
	foldonData.start = sse[0].start;
	foldonData.end = sse[sse.size() - 1].end;
	foldonData.type = sse[sse.size() - 1].type;
	foldon.push_back(foldonData);
    foldon[foldon.size() - 1].end = sse[sse.size()-1].end;
    
    // populate sample corresponding to foldon units
    vector<MDArray<double> > foldonSample;
    for (int i = 0; i < foldon.size(); i++) {
        MDArray<double> foldonSampleData;
        int foldonSize = foldon[i].end + 1;
        foldonSampleData.set_shape(2 * foldonSize, NUM_DAT);
        for (int j = 0; j <= foldon[i].end; j++) {
            // for backbone
            foldonSampleData.set(2 * j, 0, 0);
            foldonSampleData.set(2 * j, 1, 0);
            foldonSampleData.set(2 * j, 2, aa[j]);
            foldonSampleData.set(2 * j, 3, ss[j]);
            foldonSampleData.set(2 * j, 4, CA2CA);
            foldonSampleData.set(2 * j, 5, PI);
            foldonSampleData.set(2 * j, 6, PI);
            // for side chain
            foldonSampleData.set(2 * j + 1, 0, 1);
            foldonSampleData.set(2 * j + 1, 1, 0);
            foldonSampleData.set(2 * j + 1, 2, aa[j]);
            foldonSampleData.set(2 * j + 1, 3, ss[j]);
            foldonSampleData.set(2 * j + 1, 4, levittBsc[aa[j]]);
            foldonSampleData.set(2 * j + 1, 5, 0);
            foldonSampleData.set(2 * j + 1, 6, 0);
        }
        foldonSample.push_back(foldonSampleData);
    }
    
    // assemble foldon units via stepwise addition
    char pdbFoldonFile[1000];
    //vector<pdbInfo> pdbFoldon;
    vector<poseInfo> poseFoldon;
    
    
    cout << endl;
    cout << "-------------------- initial model generation -------------------------------------------------------" << endl;
    cout << "Job Id                   : " << jobId << endl;
    cout << "Protein Length           : " << aa.size() << endl;
    cout << "Amino Acid Sequence      : ";
    for (int i = 0; i < aa.size(); i++) {
        cout << seq[aa[i]];
    }
    cout << endl;
    cout << "Secondary Structure      : ";
    for (int i = 0; i < ss.size(); i++) {
        cout << sec[ss[i]];
    }
    cout << endl;
	int epoch=1;
    cout << "Number of Decoys         : " << epoch << endl;
    cout << "Monte Carlo Cycle Factor : " << numCycles/100 << endl;
    cout << "Native Structure         : ";
    if (nFile) {
        cout << "Available" << endl;
    }
    else {
        cout << "Not Available" << endl;
    }
    cout << "Number of Foldon Units   : " << foldon.size() << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
    
    // buffer statistics
    char bufStat[1000];

    for (int n = 1; n <= epoch; n++) {
        cout << "UniCon3D.io: Generating decoy " << n << " of " << epoch << endl;
        for (int i = 0; i < foldon.size(); i++) {
            
            // get residue-residue contacts from native pdb
            getRrContact(cm, foldon[i].end);
            
            if (i == 0) {
                cout << endl <<"UniCon3D.sampling: Sampling foldon unit " << i + 1 << " [residue : 1-" << foldon[i].end + 1 << "]"<< endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            else {
                cout << endl << "UniCon3D.sampling: Conditional sampling foldon unit " << i + 1 << " [residue : " << foldon[i-1].end + 2 << "-" << foldon[i].end + 1 << "]" << endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            cout << "UniCon3D.sampling: Performing simulated annealing energy minimization" << endl;
            assembleFoldon(foldonSample, i, dbn, pdbFoldon, poseFoldon);

            cout << endl;
            //sprintf(pdbFoldonFile, "%s_D%.6d_F%.6d.pdb", jobId, n, i+1);
            //if (saveInd == 0) {
            //    writePdb(pdbFoldon, pdbFoldonFile);
            //}
        }
        //sprintf(pdbFoldonFile, "%s_%.6d.pdb", jobId, n);
        //writePdb(pdbFoldon, pdbFoldonFile);
        
        cout << endl << "UniCon3D.io: Finished generating initial model " << n << endl;
        
        double rmsd;
        if (nFile) {
            rmsd = getRmsd(pdbFoldon, native);
            cout << "UniCon3D.io: Ca_rmsd to native structure = " << rmsd << endl;
        }
        // prepare to write to folding statistics file
        energyInfo e = getWeightedEnergy(pdbFoldon, poseFoldon, dbn);
        if (nFile) {
            sprintf(bufStat,"Decoy_name initial_model E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f Ca_rmsd %8.3f", e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total, rmsd);
            cout << bufStat << endl << endl;
        }
        else {
            sprintf(bufStat,"Decoy_name initial_model E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f", e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total);
            cout << bufStat << endl << endl;
        }
        cout << "---------------------------------------------------------------------------" << endl;
    }
}
/*************************************************************************
 * Name        : backboneModelling
 * Purpose     : converts pose to pdb
 * Arguments   : vector<int> &aa,vector<int> &ss, DBN &dbn, vector<pdbInfo> &pdbFoldon
 * Return Type : void
 *************************************************************************/
void backboneModelling(vector<int> &aa,vector<int> &ss,MDArray<double> &sample,DBN &dbn, vector<pdbInfo> &pdbFoldon) {
		
    // extract secondary structural elements
    vector<sseInfo> sse;
    bool sseStarted = false;
    sseInfo sseData;
    for (int i = 0; i < ss.size() - 1; i++) {
        if (!sseStarted) {
            sseData.start = i;
            sseData.type = ss[i];
            sseStarted = true;
        }
        if (sseStarted && ss[i] != ss[i+1]) {
            sseData.end = i;
            sse.push_back(sseData);
            sseStarted = false;
        }
    }
    // for the last sse
    sseData.type = ss[ss.size() - 1];
    sseData.end = ss.size() - 1;
    sse.push_back(sseData);
    
    
    // extract foldon units
    vector<sseInfo> foldon;
    bool foldonStarted = false;
    sseInfo foldonData;
	
	foldonData.start = sse[0].start;
	foldonData.end = sse[sse.size() - 1].end;
	foldonData.type = sse[sse.size() - 1].type;
	foldon.push_back(foldonData);
    foldon[foldon.size() - 1].end = sse[sse.size()-1].end;
    
    // populate sample corresponding to foldon units
    vector<MDArray<double> > foldonSample;
	
	cout << " check here: " << endl;
	for (int i = 0; i < foldon.size(); i++) {
        MDArray<double> foldonSampleData;
        int foldonSize = foldon[i].end + 1;
        foldonSampleData.set_shape(2 * foldonSize, NUM_DAT);
        for (int j = 0; j <= foldon[i].end; j++) {
            // for backbone
            foldonSampleData.set(2 * j, 0, sample.get(2 * j, 0));
            foldonSampleData.set(2 * j, 1, sample.get(2 * j, 1));
			if(sample.get(2 * j, 2) != aa[j])
			{
				cout << " assign initial model: " << sample.get(2 * j, 2) << " not equal to "<< aa[j] << " at residue aa[j] " << j << endl;
				exit(0);	
			}
            foldonSampleData.set(2 * j, 2, sample.get(2 * j, 2));
			if(sample.get(2 * j, 3) != ss[j])
			{
				cout << " assign initial model: " << sample.get(2 * j, 3) << " not equal to "<< ss[j] << " at residue ss[j] " << j << endl;
				exit(0);	
			}
            foldonSampleData.set(2 * j, 3, sample.get(2 * j, 3));
            foldonSampleData.set(2 * j, 4, sample.get(2 * j, 4));
            foldonSampleData.set(2 * j, 5, sample.get(2 * j, 5));
            foldonSampleData.set(2 * j, 6, sample.get(2 * j, 6));
            // for side chain
            foldonSampleData.set(2 * j + 1, 0, sample.get(2 * j + 1, 0));
            foldonSampleData.set(2 * j + 1, 1, sample.get(2 * j + 1, 1));
            foldonSampleData.set(2 * j + 1, 2, sample.get(2 * j + 1, 2));
            foldonSampleData.set(2 * j + 1, 3, sample.get(2 * j + 1, 3));
            foldonSampleData.set(2 * j + 1, 4, sample.get(2 * j + 1, 4));
            foldonSampleData.set(2 * j + 1, 5, sample.get(2 * j + 1, 5));
            foldonSampleData.set(2 * j + 1, 6, sample.get(2 * j + 1, 6));
        }
        foldonSample.push_back(foldonSampleData);
    }
	
	
    
    // assemble foldon units via stepwise addition
    char pdbFoldonFile[1000];
    //vector<pdbInfo> pdbFoldon;
    vector<poseInfo> poseFoldon;
    
    
    cout << endl;
    cout << "-------------------- Backbone modeling -------------------------------------------------------" << endl;
    cout << "Job Id                   : " << jobId << endl;
    cout << "Protein Length           : " << aa.size() << endl;
    cout << "Amino Acid Sequence      : ";
    for (int i = 0; i < aa.size(); i++) {
        cout << seq[aa[i]];
    }
    cout << endl;
    cout << "Secondary Structure      : ";
    for (int i = 0; i < ss.size(); i++) {
        cout << sec[ss[i]];
    }
    cout << endl;
	int epoch=1;
    cout << "Number of Decoys         : " << epoch << endl;
    cout << "Monte Carlo Cycle Factor : " << numCycles/100 << endl;
    cout << "Native Structure         : ";
    if (nFile) {
        cout << "Available" << endl;
    }
    else {
        cout << "Not Available" << endl;
    }
    cout << "Number of Foldon Units   : " << foldon.size() << endl;
    cout << "---------------------------------------------------------------------------" << endl;
    
    
    // buffer statistics
    char bufStat[1000];

    for (int n = 1; n <= epoch; n++) {
        cout << "UniCon3D.io: Generating decoy " << n << " of " << epoch << endl;
        for (int i = 0; i < foldon.size(); i++) {
            
            // get residue-residue contacts from native pdb
            //getRrContact(cm, foldon[i].end);
            
            if (i == 0) {
                cout << endl <<"UniCon3D.sampling: Sampling foldon unit " << i + 1 << " [residue : 1-" << foldon[i].end + 1 << "]"<< endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            else {
                cout << endl << "UniCon3D.sampling: Conditional sampling foldon unit " << i + 1 << " [residue : " << foldon[i-1].end + 2 << "-" << foldon[i].end + 1 << "]" << endl << endl;
                //cout << "---------------------------------------------------------------------------" << endl;
            }
            cout << "UniCon3D.sampling: Performing simulated annealing energy minimization" << endl;
            backboneSampling(foldonSample, i, dbn, pdbFoldon, poseFoldon);

            cout << endl;
            //sprintf(pdbFoldonFile, "%s_D%.6d_F%.6d.pdb", jobId, n, i+1);
            //if (saveInd == 0) {
            //    writePdb(pdbFoldon, pdbFoldonFile);
            //}
        }
        //sprintf(pdbFoldonFile, "%s_%.6d.pdb", jobId, n);
        //writePdb(pdbFoldon, pdbFoldonFile);
        
        cout << endl << "UniCon3D.io: Finished generating decoy " << pdbFoldonFile << endl;
        
        
        // prepare to write to folding statistics file
        energyInfo e = getWeightedEnergy(pdbFoldon, poseFoldon, dbn);
        
        sprintf(bufStat,"Decoy_name %s E_sc_sc %8.3f E_sc_bb %8.3f E_bb_bb %8.3f E_ri_rj %8.3f E_total %8.3f", pdbFoldonFile, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total);
            
        
        cout << "---------------------------------------------------------------------------" << endl;
    }
}
/*************************************************************************
 * Name        : isHydrophobic
 * Purpose     : returns true if residue is hydrophobic, false otherwise
 * Arguments   : int aa
 * Return Type : bool
 *************************************************************************/
bool isHydrophobic(int aa) {
    if ( aa == CYS || aa == MET || aa == PHE || aa == ILE || aa == LEU || aa == VAL || aa == TRP || aa == TYR || aa == ALA) {
        return true;
    }
    else {
        return false;
    }
}

/*************************************************************************
 * Name        : getRrContact
 * Purpose     : gets the residue-residue contacts upto a specific residue
 * Arguments   : MDArray<double> &cm, int end
 * Return Type : void
 *************************************************************************/
void getRrContact(MDArray<double> &cm, int end) {
    rrContact.clear();
    for (int i = 0; i < end; i++) {
        for (int j = i + RR_SEQ_SEP; j < end; j++) {
            contactInfo contactData;
            contactData.ri = i;
            contactData.rj = j;
            contactData.wt = cm.get(i, j);
            rrContact.push_back(contactData);
        }
    }
}

/*************************************************************************
 * Name        : getUpperBoundHarmonic
 * Purpose     : gets upper bound harmonic potential
 * Arguments   : double d, double bound
 * Return Type : double
 *************************************************************************/
double getUpperBoundHarmonic(double d, double bound) {
    double potential = 0.0;
    if (d > bound) {
        potential = ((d - bound) * (d - bound));
    }
    return potential;
}

/*************************************************************************
 * Name        : getFade
 * Purpose     : gets FADE potential
 * Arguments   : double d, double lb, double ub, double z, double w
 * Return Type : double
 *************************************************************************/
double getFade(double d, double lb, double ub, double z, double w) {
    double potential = 0.0;
    double lf = lb + z;
    double uf = ub - z;
    if (d < lb || d > ub) {
        potential = 0.0;
    }
    else if (d < lf) {
        potential = -2.0 * pow(((d - lf) / z), 3.0) - 3.0 * pow(((d - lf) / z), 2.0);
    }
    else if (d < uf) {
        potential = -2.0 * pow(((d - uf) / z), 3.0) - 3.0 * pow(((d - uf) / z), 2.0);
    }
    else {
        potential = w;
    }
    return potential;
}

/*************************************************************************
 * Name        : getLiwoEpsilon0
 * Purpose     : gets the epsilon 0 values for amino acid pairs
 * Arguments   : int ai, int aj
 * Return Type : double
 *************************************************************************/
double getLiwoEpsilon0(int ai, int aj) {
    int id1 = liwoAminoStr.find(seq[ai]);
    int id2 = liwoAminoStr.find(seq[aj]);
    if (id1 < 0 || id1 > 19) {
        return 0;
    }
    if (id2 < 0 || id2 > 19) {
        return 0;
    }
    // swap in case id1 is less than id2
    if (id1 < id2)
    {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
    }
    return liwoEpsilon0[id1][id2];
}

/*************************************************************************
 * Name        : getScScEnergy
 * Purpose     : gets the side chain - side chain interaction energy (GBV)
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getScScEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < pdb.size() - 1; i++) {
        for (int j = i + 1; j < pdb.size(); j++) {
            // consider only secondary structural elements
            if (pose[i].ss != RANDOM_COIL && pose[j].ss != RANDOM_COIL) {
            //if (pose[i].ss != RANDOM_COIL && pose[j].ss != RANDOM_COIL) {
                double rij = getDistance(pdb[i].sc, pdb[j].sc);
                double rij0 = liwoR0[pose[i].aa] + liwoR0[pose[j].aa];
                double sigmaij0 = sqrt(liwoSigma0[pdb[i].aa] * liwoSigma0[pdb[i].aa] + liwoSigma0[pdb[j].aa] * liwoSigma0[pdb[j].aa]);
                point3d dij1 = getDifference(pdb[i].ca, pdb[i].sc);
                point3d dij2 = getDifference(pdb[j].ca, pdb[j].sc);
                point3d drij = getDifference(pdb[i].sc, pdb[j].sc);
                point3d uij1 = getUnit(dij1);
                point3d uij2 = getUnit(dij2);
                point3d urij2 = getUnit(drij);
                double omegaij1 = getDotProduct(uij1, urij2);
                double omegaij2 = getDotProduct(uij2, urij2);
                double omegaij12 = getDotProduct(uij1, uij2);
                double chiij1Nr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[i].aa] - 1.0;
                double chiij1Dr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[i].aa] + pow((liwoSigma0[pdb[j].aa]/liwoSigma0[pdb[i].aa]), 2.0);
                double chiij1 = chiij1Nr / chiij1Dr;
                double chiij2Nr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[j].aa] - 1.0;
                double chiij2Dr = liwoSigmaDoubleBarOverSigmaTeeSquare[pdb[j].aa] + pow((liwoSigma0[pdb[i].aa]/liwoSigma0[pdb[j].aa]), 2.0);
                double chiij2 = chiij2Nr / chiij2Dr;
                double epsilonij0 = getLiwoEpsilon0(pdb[i].aa, pdb[j].aa);
                double epsilonij1Dr = (1.0 - chiij1 * chiij2 * omegaij12);
                double epsilonij1 = sqrt(epsilonij1Dr);
                double epsilonij2Nr = (liwoChiPrime[pdb[i].aa] * omegaij1 * omegaij1 + liwoChiPrime[pdb[j].aa] * omegaij2 * omegaij2 - 2.0 * liwoChiPrime[pdb[i].aa] * liwoChiPrime[pdb[j].aa] * omegaij1 * omegaij2 * omegaij12);
                double epsilonij2Dr = (1.0 - liwoChiPrime[pose[i].aa] * liwoChiPrime[pose[j].aa] * omegaij12 * omegaij12);
                double epsilonij2 = (1.0 - epsilonij2Nr / epsilonij2Dr) * (1.0 - epsilonij2Nr / epsilonij2Dr);
                double epsilonij3Nr = (1.0 - liwoAlpha[pdb[i].aa] * omegaij1 + liwoAlpha[pdb[j].aa] * omegaij2 - 0.5 * (liwoAlpha[pdb[i].aa] + liwoAlpha[pdb[j].aa]) * omegaij12);
                double epsilonij3 = epsilonij3Nr * epsilonij3Nr;
                double sigmaijNr = (chiij1 * omegaij1 * omegaij1 + chiij2 * omegaij2 * omegaij2 - 2.0 * chiij1 * chiij2 * omegaij1 * omegaij2 * omegaij12);
                double sigmaijDr = (1.0 - chiij1 * chiij2 * omegaij12 * omegaij12);
                double sigmaij = sigmaij0 * sqrt(1.0 - sigmaijNr / sigmaijDr);;
                double epsilonij = epsilonij0 * epsilonij1 * epsilonij2 * epsilonij3;
                double xijNr = rij0;
                double xijDr = (rij - sigmaij + rij0);
                double xij = 0.0;
                if (xijDr != 0.0) {
                    xij = xijNr / xijDr;
                }
                energy += 4 * (abs(epsilonij) * pow(xij, 12.0) - epsilonij * pow(xij, 6.0));
            }
        }
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getScBbEnergy
 * Purpose     : gets the side chain - backbone interaction energy (repulsive only)
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getScBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < pdb.size(); i++) {
        for (int j = 0; j < pdb.size() - 1; j++) {
            if (j != i - 1 && j != i) {
                point3d pj = getMidpoint(pdb[j].ca, pdb[j+1].ca);
                double rij = getDistance(pdb[i].sc, pj);
                energy += 0.3 * pow((4.0 / rij), 6.0);
            }
        }
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getBbBbEnergy
 * Purpose     : gets the backnone - backbone interaction energy (electrostatic)
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getBbBbEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < pdb.size() - 3; i++) {
        for (int j = i + 2; j < pdb.size() - 1; j++) {
            double A = 0.0;
            double B = 0.0;
            double epsilon = 0.0;
            double r0 = 0.0;
            // proline-proline interaction
            if (pdb[i+1].aa == PRO && pdb[j+1].aa == PRO) {
                epsilon = 0.574;
                r0 = 4.48;
                A = 5.13;
                B = 335.0;
            }
            // ordinary-proline interaction
            if ((pdb[i+1].aa == PRO && pdb[j+1].aa != PRO) || (pdb[i+1].aa != PRO && pdb[j+1].aa == PRO)) {
                epsilon = 0.365;
                r0 = 4.54;
                A = 0.0;
                B = 1129.0;
            }
            // ordinary-ordinary interaction
            if (pdb[i+1].aa != PRO && pdb[j+1].aa != PRO) {
                epsilon = 0.305;
                r0 = 4.51;
                A = 3.73;
                B = 1306.0;
            }
            point3d vi = getDifference(pdb[i].ca, pdb[i+1].ca);
            point3d vj = getDifference(pdb[j].ca, pdb[j+1].ca);
            point3d uvi = getUnit(vi);
            point3d uvj = getUnit(vj);
            point3d pi = getMidpoint(pdb[i].ca, pdb[i+1].ca);
            point3d pj = getMidpoint(pdb[j].ca, pdb[j+1].ca);
            point3d pij = getDifference(pi, pj);
            point3d upij = getUnit(pij);
            double rij = getDistance(pi, pj);
            double cosAlphaij = getDotProduct(uvi, uvj);
            double cosBetaij = getDotProduct(uvi, upij);
            double cosGammaij = getDotProduct(uvj, upij);
            double term1 = (A / pow(rij, 3.0)) * (cosAlphaij - 3.0 * cosBetaij * cosGammaij);
            double term2 = (B / pow(rij, 6.0)) * (4.0 + ((cosAlphaij - 3.0 * cosBetaij * cosGammaij) * (cosAlphaij - 3.0 * cosBetaij * cosGammaij)) - 3.0 * (cosBetaij * cosBetaij + cosGammaij * cosGammaij));
            double term3 = epsilon * (pow((r0 / rij), 12.0) - 2.0 * pow((r0 / rij), 6.0));
            double total = term1 - term2 + term3;
            energy += total;
        }
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getRiRjEnergy
 * Purpose     : gets the residue-residue contact energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight
 * Return Type : double
 *************************************************************************/
double getRiRjEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, double weight) {
    double energy = 0.0;
    for (int i = 0; i < rrContact.size(); i++) {
        contactInfo contactData;
        int ci = rrContact[i].ri;
        int cj = rrContact[i].rj;
        double dij = getDistance(pdb[ci].sc, pdb[cj].sc);
        double pij = rrContact[i].wt;
        if (dij <= RR_CONTACT_THRESHOLD ) {
            energy += -pij;
        }
        else {
            energy += (-pij * exp(-pow((dij - RR_CONTACT_THRESHOLD), 2.0))) + (pij * ((dij - RR_CONTACT_THRESHOLD)/dij));
        }
        //energy += (log(abs(ci - cj)) * pij * (dij - RR_CONTACT_THRESHOLD));
    }
    return energy * weight;
}

/*************************************************************************
 * Name        : getWeightedEnergy
 * Purpose     : gets total weighted energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn
 * Return Type : energyInfo
 *************************************************************************/
energyInfo getWeightedEnergy(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn) {
    energyInfo e;
    e.sc_sc = getScScEnergy(pdb, pose, 1.0);
    e.sc_bb = getScBbEnergy(pdb, pose, 1.0);
    e.bb_bb = getBbBbEnergy(pdb, pose, 1.0);
    e.ri_rj = getRiRjEnergy(pdb, pose, 1.0);
    e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj;
    return e;
}

/*************************************************************************
 * Name        : showEnergy
 * Purpose     : shows (prints) the energy breakdown and total weighted energy
 * Arguments   : energyInfo
 * Return Type : void
 *************************************************************************/
void showEnergy(energyInfo &e) {
    char buf[1000];
    sprintf(buf,"UniCon3D.scoring: E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_total = %8.3f", e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.total);
    printf("%s\r",buf);
    //fflush(stdout);
    //cout << buf << endl;
}

/*************************************************************************
 * Name        : getRmsd
 * Purpose     : gets the root mean square deviation between two pdb objects
 * Arguments   : vector<pdbInfo> &pdb1, vector<pdbInfo> &pdb2, string mode
 * Return Type : double
 *************************************************************************/
double getRmsd(vector<pdbInfo> &pdb1, vector<pdbInfo> &pdb2, string mode) {
    if (pdb1.size() != pdb2.size()) {
        cout << "Error! The two compared structures are of different length" << endl;
        exit(0);
    }
    double rot[3][3];
    double trans[3];
    vector<vector<double> > A;
    for (int i = 0; i < pdb1.size(); i++) {
        vector<double> a(3);
        //load CA only
        if (mode.compare("ca") == 0) {
            a[0] = pdb1[i].ca.x;
            a[1] = pdb1[i].ca.y;
            a[2] = pdb1[i].ca.z;
            A.push_back(a);
        }
        //load SC only
        else if (mode.compare("sc") == 0) {
            a[0] = pdb1[i].sc.x;
            a[1] = pdb1[i].sc.y;
            a[2] = pdb1[i].sc.z;
            A.push_back(a);
        }
        //load all atoms (CA + SC)
        else if (mode.compare("aa") == 0) {
            a[0] = pdb1[i].ca.x;
            a[1] = pdb1[i].ca.y;
            a[2] = pdb1[i].ca.z;
            A.push_back(a);
            a[0] = pdb1[i].sc.x;
            a[1] = pdb1[i].sc.y;
            a[2] = pdb1[i].sc.z;
            A.push_back(a);
        }
        else {
            cout << "Error! Invalid mode in rmsd calculation" << endl;
            exit(0);
        }
    }
    vector<vector<double> > B;
    for (int i = 0; i < pdb2.size(); i++) {
        vector<double> b(3);
        //load CA only
        if (mode.compare("ca") == 0) {
            b[0] = pdb2[i].ca.x;
            b[1] = pdb2[i].ca.y;
            b[2] = pdb2[i].ca.z;
            B.push_back(b);
        }
        //load SC only
        else if (mode.compare("sc") == 0) {
            b[0] = pdb2[i].sc.x;
            b[1] = pdb2[i].sc.y;
            b[2] = pdb2[i].sc.z;
            B.push_back(b);
        }
        //load all atoms (CA + SC)
        else if (mode.compare("aa") == 0) {
            b[0] = pdb2[i].ca.x;
            b[1] = pdb2[i].ca.y;
            b[2] = pdb2[i].ca.z;
            B.push_back(b);
            b[0] = pdb2[i].sc.x;
            b[1] = pdb2[i].sc.y;
            b[2] = pdb2[i].sc.z;
            B.push_back(b);
        }
        else {
            cout << "Error! Invalid mode in rmsd calculation" << endl;
            exit(0);
        }
    }
    double Bt[A.size()][3];
    getAlignment(A, B, rot, trans);
    // Transform the second chain to optimally align with the first.
    for (int k = 0; k < A.size(); k++) {
        Bt[k][X] = B[k][X] * rot[0][0] + B[k][Y] * rot[1][0] +
        B[k][Z] * rot[2][0] + trans[0];
        Bt[k][Y] = B[k][X] * rot[0][1] + B[k][Y] * rot[1][1] +
        B[k][Z] * rot[2][1] + trans[1];
        Bt[k][Z] = B[k][X] * rot[0][2] + B[k][Y] * rot[1][2] +
        B[k][Z] * rot[2][2] + trans[2];
    }
    double rmsd = 0;
    for (int i = 0; i < A.size(); i++) {
        double a0 = Bt[i][X] - A[i][X];
        double a1 = Bt[i][Y] - A[i][Y];
        double a2 = Bt[i][Z] - A[i][Z];
        
        rmsd += (a0*a0 + a1*a1 + a2*a2);
    }
    return sqrt(rmsd / A.size());
}

/*************************************************************************
 * Name        : getAlignment
 * Purpose     : gets the optimal alignemnt between two sets of point clouds
 * Arguments   : vector<vector<double> > &  A, vector<vector<double> > & B, double rot[3][3], double trans[3]
 * Return Type : void
 *************************************************************************/
void getAlignment(vector<vector<double> > &  A, vector<vector<double> > & B, double rot[3][3], double trans[3]) {
    int i,j;
    double c1[3],c2[3];   /* center of mass for two point collections */
    double v1[3],v2[3];
    double recip;
    double tr;
    double m[4][4], q[4][4];
    double v[4];
    double cov[3][3];
    double aij[3][3];
    double quat[4];
    // find the center of mass for the two collections
    c1[X] = c1[Y] = c1[Z] = 0;
    c2[X] = c2[Y] = c2[Z] = 0;
    
    for (int i = 0; i < A.size(); i++) {
        c1[X] += A[i][X];
        c1[Y] += A[i][Y];
        c1[Z] += A[i][Z];
        
        c2[X] += B[i][X];
        c2[Y] += B[i][Y];
        c2[Z] += B[i][Z];
    }
    recip = 1.0 / A.size();
    c1[X] *= recip;
    c1[Y] *= recip;
    c1[Z] *= recip;
    c2[X] *= recip;
    c2[Y] *= recip;
    c2[Z] *= recip;
    // create the cross-covariance matrix
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            cov[i][j] = 0;
    for (int i = 0; i < A.size(); i++) {
        v1[X] = A[i][X] - c1[X];
        v1[Y] = A[i][Y] - c1[Y];
        v1[Z] = A[i][Z] - c1[Z];
        v2[X] = B[i][X] - c2[X];
        v2[Y] = B[i][Y] - c2[Y];
        v2[Z] = B[i][Z] - c2[Z];
        cov[X][X] += v1[X] * v2[X];
        cov[X][Y] += v1[X] * v2[Y];
        cov[X][Z] += v1[X] * v2[Z];
        cov[Y][X] += v1[Y] * v2[X];
        cov[Y][Y] += v1[Y] * v2[Y];
        cov[Y][Z] += v1[Y] * v2[Z];
        cov[Z][X] += v1[Z] * v2[X];
        cov[Z][Y] += v1[Z] * v2[Y];
        cov[Z][Z] += v1[Z] * v2[Z];
    }
    // aij = cov - transpose(cov)
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            aij[i][j] = cov[i][j] - cov[j][i];
    // find the trace of the covariance matrix
    tr = cov[X][X] + cov[Y][Y] + cov[Z][Z];
    m[0][0] = tr;
    m[1][0] = m[0][1] = aij[1][2];
    m[2][0] = m[0][2] = aij[2][0];
    m[3][0] = m[0][3] = aij[0][1];
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            m[i+1][j+1] = cov[i][j] + cov[j][i] - (i == j) * tr;
    // find the eigenvector corresponding to the largest eigenvalue of this matrix
    Meigen4(q, v, m);
    if( v[0] > v[1] ) {
        if( v[0] > v[2] ) {
            if( v[0] > v[3] )
                i = 0;
            else
                i = 3;
        }
        else{
            if( v[2] > v[3] )
                i = 2;
            else
                i =3;
        }
    }
    else {
        if( v[1] > v[2] ){
            if( v[1] > v[3] )
                i = 1;
            else
                i = 3;
        }
        else {
            if( v[2] > v[3] )
                i = 2;
            else
                i =3;
        }
    }
    quat[0] = q[0][i];
    quat[1] = q[1][i];
    quat[2] = q[2][i];
    quat[3] = q[3][i];
    quat2mat(quat, rot);
    double c3[3];
    // determine best translation
    vApply (rot, c2, c3);
    trans[0] = c1[0] - c3[0];
    trans[1] = c1[1] - c3[1];
    trans[2] = c1[2] - c3[2];
}

/*************************************************************************
 * Name        : quat2mat
 * Purpose     : convert quat to mat
 * Arguments   : double q[4], double mat[3][3]
 * Return Type : void
 *************************************************************************/
void quat2mat(double q[4], double mat[3][3]) {
    double q00,q01,q02,q03;
    double q11,q12,q13;
    double q22,q23;
    double q33;
    q00 = q[0] * q[0];
    q01 = q[0] * q[1];
    q02 = q[0] * q[2];
    q03 = q[0] * q[3];
    q11 = q[1] * q[1];
    q12 = q[1] * q[2];
    q13 = q[1] * q[3];
    q22 = q[2] * q[2];
    q23 = q[2] * q[3];
    q33 = q[3] * q[3];
    mat[X][X] = q00 + q11 - q22 - q33;
    mat[X][Y] = 2 * (q12 - q03);
    mat[X][Z] = 2 * (q13 + q02);
    mat[Y][X] =  2 * (q12 + q03);
    mat[Y][Y] =  q00 + q22 - q11 - q33;
    mat[Y][Z] = 2 * (q23 - q01);
    mat[Z][X] = 2 * (q13 - q02);
    mat[Z][Y] = 2 * (q23 + q01);
    mat[Z][Z] = q00 + q33 - q11 - q22;
}

/*************************************************************************
 * Name        : vApply
 * Purpose     : apply translation
 * Arguments   : double m[3][3], const double a[3], double b[3]
 * Return Type : void
 *************************************************************************/
void vApply (double m[3][3], const double a[3], double b[3]) {
    int j;
    for (j = 0; j <= 2; j++)
        b[j] = a[0] * m[0][j] + a[1] * m [1][j] +
        a[2] * m[2][j];
}

/*************************************************************************
 * Name        : assembleFoldon
 * Purpose     : assemble foldon unit via simulated anneling
 * Arguments   : vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow
 * Return Type : void
 *************************************************************************/
void assembleFoldon(vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow) {
    int start = 0;
    if (i != 0) {
        start = foldonSample[i-1].get_shape()[0];
    }
    int end = foldonSample[i].get_shape()[0];
    // populate the current foldonSample based on the previos other than the first foldonSample
    if (i != 0) {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0]; j++) {
            foldonSample[i].set(j, 0, foldonSample[i-1].get(j, 0));
            foldonSample[i].set(j, 1, foldonSample[i-1].get(j, 1));
            foldonSample[i].set(j, 2, foldonSample[i-1].get(j, 2));
            foldonSample[i].set(j, 3, foldonSample[i-1].get(j, 3));
            foldonSample[i].set(j, 4, foldonSample[i-1].get(j, 4));
            foldonSample[i].set(j, 5, foldonSample[i-1].get(j, 5));
            foldonSample[i].set(j, 6, foldonSample[i-1].get(j, 6));
        }
    }
    // populate data and mismask from foldonSample
    Sequence data;
    data.set_shape(foldonSample[i].get_shape()[0], NUM_DAT);
    MDArray<eMISMASK> mism;
    mism.set_shape(foldonSample[i].get_shape()[0], NUM_MIS);
    if (i == 0) {
        for (int j = 0; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            mism.set(j * 2, 4, MOCAPY_HIDDEN);
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            mism.set(j * 2, 5, MOCAPY_HIDDEN);
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
            }
        }
    }
    else {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
                mism.set(j * 2, 4, MOCAPY_HIDDEN);
            }
            else {
                mism.set(j * 2, 4, MOCAPY_OBSERVED);
            }
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
                mism.set(j * 2, 5, MOCAPY_HIDDEN);
            }
            else {
                mism.set(j * 2, 5, MOCAPY_OBSERVED);
            }
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
            }
        }
        for (int j = foldonSample[i-1].get_shape()[0] / 2; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            mism.set(j * 2, 4, MOCAPY_HIDDEN);
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            mism.set(j * 2, 5, MOCAPY_HIDDEN);
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
            }
        }
    }
    // setup the sampler
    SampleInfEngineHMM sampler(&dbn, data, mism, 1);
    //MDArray<double> sample = sampler.sample_next();
    //sampler.undo();
    sampler.set_start_end(start, end);
    MDArray<double> sample = sampler.sample_next();

    double Ti = START_TEMP;
    double Tf = FINAL_TEMP;
    
    // perform simulated anneling
    doSimulatedAnneling(sampler, dbn, data, mism, sample, poseLow, MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf);
    // populate the lowest scoring sample into current foldonSample
    for (int j = 0; j < sample.get_shape()[0]; j++) {
        foldonSample[i].set(j, 0, sample.get(j, 0));
        foldonSample[i].set(j, 1, sample.get(j, 1));
        foldonSample[i].set(j, 2, sample.get(j, 2));
        foldonSample[i].set(j, 3, sample.get(j, 3));
        foldonSample[i].set(j, 4, sample.get(j, 4));
        foldonSample[i].set(j, 5, sample.get(j, 5));
        foldonSample[i].set(j, 6, sample.get(j, 6));
    }
    pose2pdb(poseLow, pdbLow);
}


/*************************************************************************
 * Name        : assembleFoldon
 * Purpose     : assemble foldon unit via simulated anneling
 * Arguments   : vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow
 * Return Type : void
 *************************************************************************/
void backboneSampling(vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow) {
    int start = 0;
    if (i != 0) {
        start = foldonSample[i-1].get_shape()[0];
    }
    int end = foldonSample[i].get_shape()[0];
    // populate the current foldonSample based on the previos other than the first foldonSample
    if (i != 0) {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0]; j++) {
            foldonSample[i].set(j, 0, foldonSample[i-1].get(j, 0));
            foldonSample[i].set(j, 1, foldonSample[i-1].get(j, 1));
            foldonSample[i].set(j, 2, foldonSample[i-1].get(j, 2));
            foldonSample[i].set(j, 3, foldonSample[i-1].get(j, 3));
            foldonSample[i].set(j, 4, foldonSample[i-1].get(j, 4));
            foldonSample[i].set(j, 5, foldonSample[i-1].get(j, 5));
            foldonSample[i].set(j, 6, foldonSample[i-1].get(j, 6));
        }
    }
    // populate data and mismask from foldonSample
    Sequence data;
    data.set_shape(foldonSample[i].get_shape()[0], NUM_DAT);
    MDArray<eMISMASK> mism;
    mism.set_shape(foldonSample[i].get_shape()[0], NUM_MIS);
    if (i == 0) {
        for (int j = 0; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            mism.set(j * 2, 4, MOCAPY_HIDDEN); 
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            mism.set(j * 2, 5, MOCAPY_OBSERVED);// not sample angle
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED); // not sample angle
            }
        }
    }
    else {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
                mism.set(j * 2, 4, MOCAPY_HIDDEN);
            }
            else {
                mism.set(j * 2, 4, MOCAPY_OBSERVED);
            }
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
                mism.set(j * 2, 5, MOCAPY_OBSERVED); // not sample angle
            }
            else {
                mism.set(j * 2, 5, MOCAPY_OBSERVED);
            }
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED); // not sample angle
            }
        }
        for (int j = foldonSample[i-1].get_shape()[0] / 2; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
            // for backbone
            data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
            mism.set(j * 2, 0, MOCAPY_OBSERVED);
            data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
            mism.set(j * 2, 1, MOCAPY_HIDDEN);
            data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
            mism.set(j * 2, 2, MOCAPY_OBSERVED);
            data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
            mism.set(j * 2, 3, MOCAPY_OBSERVED);
            data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
            mism.set(j * 2, 4, MOCAPY_HIDDEN);
            data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
            data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
            mism.set(j * 2, 5, MOCAPY_OBSERVED); // not sample angle
            // for side chain
            data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
            mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
            mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
            data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
            mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
            mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
            data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
            }
            data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
            data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
            if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
            }
            else {
                mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED); // not sample angle
            }
        }
    }
    // setup the sampler
    SampleInfEngineHMM sampler(&dbn, data, mism, 1);
    //MDArray<double> sample = sampler.sample_next();
    //sampler.undo();
    sampler.set_start_end(start, end);
    MDArray<double> sample = sampler.sample_next();

    double Ti = START_TEMP;
    double Tf = FINAL_TEMP;
    
    // perform simulated anneling
    doSimulatedAnneling(sampler, dbn, data, mism, sample, poseLow, MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf);
    // populate the lowest scoring sample into current foldonSample
    for (int j = 0; j < sample.get_shape()[0]; j++) {
        foldonSample[i].set(j, 0, sample.get(j, 0));
        foldonSample[i].set(j, 1, sample.get(j, 1));
        foldonSample[i].set(j, 2, sample.get(j, 2));
        foldonSample[i].set(j, 3, sample.get(j, 3));
        foldonSample[i].set(j, 4, sample.get(j, 4));
        foldonSample[i].set(j, 5, sample.get(j, 5));
        foldonSample[i].set(j, 6, sample.get(j, 6));
    }
    pose2pdb(poseLow, pdbLow);
}



/*************************************************************************
 * Name        : sample2pose
 * Purpose     : converts sample to pose
 * Arguments   : MDArray<double> &sample, vector<poseInfo> &pose
 * Return Type : void
 *************************************************************************/
void sample2file(MDArray<double> &sample,  char * filename) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cout << "Error! output pdb file can not open " << filename << endl;
        return;
    }
    for (int i = 0; i < sample.get_shape()[0] / 2; i++) {
        char bufCa[1000];
        char bufSC[1000];
        sprintf(bufCa,"%6d CA %6d %6d %6.3f %6.3f %6.3f",i + 1,(int)sample.get(2 * i, 2),(int)sample.get(2 * i, 3),sample.get(2 * i, 4),sample.get(2 * i, 5),sample.get(2 * i, 6));
        sprintf(bufSC,"%6d SC %6d %6d %6.3f %6.3f %6.3f",i + 1,(int)sample.get(2 * i, 2),(int)sample.get(2 * i, 3),sample.get(2 * i+1, 4),sample.get(2 * i+1, 5),sample.get(2 * i+1, 6));
        fout << bufCa << endl;
		fout << bufSC << endl;
    }
}


/*************************************************************************
 * Name        : assembleDomain
 * Purpose     : assemble foldon unit via simulated anneling
 * Arguments   : vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow
 * Return Type : void
 *************************************************************************/
void assembleDomain(vector<MDArray<double> > &foldonSample, int i,  vector<int> &alnCode, DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow) {
    int start = 0;
    if (i != 0) {
        start = foldonSample[i-1].get_shape()[0];
    }
    int end = foldonSample[i].get_shape()[0];
    // populate the current foldonSample based on the previos other than the first foldonSample
    if (i != 0) {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0]; j++) {
            foldonSample[i].set(j, 0, foldonSample[i-1].get(j, 0));
            foldonSample[i].set(j, 1, foldonSample[i-1].get(j, 1));
            foldonSample[i].set(j, 2, foldonSample[i-1].get(j, 2));
            foldonSample[i].set(j, 3, foldonSample[i-1].get(j, 3));
            foldonSample[i].set(j, 4, foldonSample[i-1].get(j, 4));
            foldonSample[i].set(j, 5, foldonSample[i-1].get(j, 5));
            foldonSample[i].set(j, 6, foldonSample[i-1].get(j, 6));
        }
    }
    // populate data and mismask from foldonSample
    Sequence data;
    data.set_shape(foldonSample[i].get_shape()[0], NUM_DAT);
    MDArray<eMISMASK> mism;
    mism.set_shape(foldonSample[i].get_shape()[0], NUM_MIS);
	for (int j = 0; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
		// for backbone
		data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
		mism.set(j * 2, 0, MOCAPY_OBSERVED);
		data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
		mism.set(j * 2, 1, MOCAPY_HIDDEN);
		data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
		mism.set(j * 2, 2, MOCAPY_OBSERVED);
		data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
		mism.set(j * 2, 3, MOCAPY_OBSERVED);
		data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
		//if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
		//	mism.set(j * 2, 4, MOCAPY_HIDDEN);
		//}
		//else {
		//	mism.set(j * 2, 4, MOCAPY_OBSERVED);
		//}
		data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
		data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
		//if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
		if ((int)alnCode[j] == 0) {
			mism.set(j * 2, 4, MOCAPY_HIDDEN);
			mism.set(j * 2, 5, MOCAPY_HIDDEN);
		}
		else {
			mism.set(j * 2, 4, MOCAPY_OBSERVED);
			mism.set(j * 2, 5, MOCAPY_OBSERVED);
		}
		// for side chain
		data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
		mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
		data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
		mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
		data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
		mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
		data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
		mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
		data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
		if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
			mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
		}
		else {
			mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
		}
		data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
		data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
		if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
			mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
		}
		else {
			mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
		}
	}
        
    
    // setup the sampler
    SampleInfEngineHMM sampler(&dbn, data, mism, 1);
    //MDArray<double> sample = sampler.sample_next();
    //sampler.undo();
    sampler.set_start_end(start, end);
    MDArray<double> sample = sampler.sample_next();

    double Ti = START_TEMP;
    double Tf = FINAL_TEMP;
    //cout << "Start simulated annealing" <<endl;
    // perform simulated anneling
    doSimulatedAnneling(sampler, dbn, data, mism, sample, poseLow, MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf);
    //cout << "finish simulated annealing" <<endl;
    // populate the lowest scoring sample into current foldonSample
    for (int j = 0; j < sample.get_shape()[0]; j++) {
        foldonSample[i].set(j, 0, sample.get(j, 0));
        foldonSample[i].set(j, 1, sample.get(j, 1));
        foldonSample[i].set(j, 2, sample.get(j, 2));
        foldonSample[i].set(j, 3, sample.get(j, 3));
        foldonSample[i].set(j, 4, sample.get(j, 4));
        foldonSample[i].set(j, 5, sample.get(j, 5));
        foldonSample[i].set(j, 6, sample.get(j, 6));
    }
    pose2pdb(poseLow, pdbLow);
}

/*************************************************************************
 * Name        : assembleDomainSAXS
 * Purpose     : assemble foldon unit via simulated anneling
 * Arguments   : vector<MDArray<double> > &foldonSample, int i, DBN &dbn, vector<pdbInfo> &pdbLow
 * Return Type : void
 *************************************************************************/
void assembleDomainSAXS_multiDomain(vector<MDArray<double> > &foldonSample, int i,  vector<int> &alnCode, vector<domLinkerInfo> linkers_info_array , DBN &dbn, vector<pdbInfo> &pdbLow, vector<poseInfo> &poseLow, int decoyId, int scoreType, bool combineScore) {
    int start = 0;
    if (i != 0) {
        start = foldonSample[i-1].get_shape()[0];
    }
    int end = foldonSample[i].get_shape()[0];
    // populate the current foldonSample based on the previos other than the first foldonSample
    if (i != 0) {
        for (int j = 0; j < foldonSample[i-1].get_shape()[0]; j++) {
            foldonSample[i].set(j, 0, foldonSample[i-1].get(j, 0));
            foldonSample[i].set(j, 1, foldonSample[i-1].get(j, 1));
            foldonSample[i].set(j, 2, foldonSample[i-1].get(j, 2));
            foldonSample[i].set(j, 3, foldonSample[i-1].get(j, 3));
            foldonSample[i].set(j, 4, foldonSample[i-1].get(j, 4));
            foldonSample[i].set(j, 5, foldonSample[i-1].get(j, 5));
            foldonSample[i].set(j, 6, foldonSample[i-1].get(j, 6));
        }
    }
    // populate data and mismask from foldonSample
    Sequence data;
    data.set_shape(foldonSample[i].get_shape()[0], NUM_DAT);
    MDArray<eMISMASK> mism;
    mism.set_shape(foldonSample[i].get_shape()[0], NUM_MIS);
	

	for (int j = 0; j < foldonSample[i].get_shape()[0] / 2 ; j++) {
		// for backbone
		//if ((int)foldonSample[i].get(j * 2 , 3) == HBONDED_TURN || (int)foldonSample[i].get(j * 2 , 3) == NONHBONDED_BEND || (int)foldonSample[i].get(j * 2 , 3) == RANDOM_COIL) {
		if ((int)alnCode[j] == 0) {
			data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
			mism.set(j * 2, 0, MOCAPY_OBSERVED);
			data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
			mism.set(j * 2, 1, MOCAPY_HIDDEN);
			data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
			mism.set(j * 2, 2, MOCAPY_OBSERVED);
			data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
			mism.set(j * 2, 3, MOCAPY_OBSERVED);
			data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
			data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
			data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
			mism.set(j * 2, 4, MOCAPY_HIDDEN);
			mism.set(j * 2, 5, MOCAPY_HIDDEN);
		}
		else {
			data.set(j * 2, 0, foldonSample[i].get(j * 2, 0));
			mism.set(j * 2, 0, MOCAPY_OBSERVED);
			data.set(j * 2, 1, foldonSample[i].get(j * 2, 1));
			mism.set(j * 2, 1, MOCAPY_HIDDEN);
			data.set(j * 2, 2, foldonSample[i].get(j * 2, 2));
			mism.set(j * 2, 2, MOCAPY_OBSERVED);
			data.set(j * 2, 3, foldonSample[i].get(j * 2, 3));
			mism.set(j * 2, 3, MOCAPY_OBSERVED);
			data.set(j * 2, 4, foldonSample[i].get(j * 2, 4));
			data.set(j * 2, 5, foldonSample[i].get(j * 2, 5));
			data.set(j * 2, 6, foldonSample[i].get(j * 2, 6));
			mism.set(j * 2, 4, MOCAPY_OBSERVED);
			mism.set(j * 2, 5, MOCAPY_OBSERVED);
		}
		// for side chain

		if ((int)alnCode[j] == 0) {
			data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
			mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
			data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
			mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
			data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
			mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
			data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
			mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
			data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
			data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
			data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
			if ((int)foldonSample[i].get(j * 2 + 1, 2) == GLY) {
				mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
				mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
			}else{
				mism.set(j * 2 + 1, 4, MOCAPY_HIDDEN);
				mism.set(j * 2 + 1, 5, MOCAPY_HIDDEN);
			}
		}
		else {
			data.set(j * 2 + 1, 0, foldonSample[i].get(j * 2 + 1, 0));
			mism.set(j * 2 + 1, 0, MOCAPY_OBSERVED);
			data.set(j * 2 + 1, 1, foldonSample[i].get(j * 2 + 1, 1));
			mism.set(j * 2 + 1, 1, MOCAPY_HIDDEN);
			data.set(j * 2 + 1, 2, foldonSample[i].get(j * 2 + 1, 2));
			mism.set(j * 2 + 1, 2, MOCAPY_OBSERVED);
			data.set(j * 2 + 1, 3, foldonSample[i].get(j * 2 + 1, 3));
			mism.set(j * 2 + 1, 3, MOCAPY_OBSERVED);
			data.set(j * 2 + 1, 4, foldonSample[i].get(j * 2 + 1, 4));
			data.set(j * 2 + 1, 5, foldonSample[i].get(j * 2 + 1, 5));
			data.set(j * 2 + 1, 6, foldonSample[i].get(j * 2 + 1, 6));
			mism.set(j * 2 + 1, 4, MOCAPY_OBSERVED);
			mism.set(j * 2 + 1, 5, MOCAPY_OBSERVED);
		}
	}
        
    
    // setup the sampler
    SampleInfEngineHMM sampler(&dbn, data, mism, 1);
    //MDArray<double> sample = sampler.sample_next();
    //sampler.undo();
    sampler.set_start_end(start, end);
    MDArray<double> sample = sampler.sample_next();

    double Ti = START_TEMP;
    double Tf = FINAL_TEMP;
    //cout << "Start simulated annealing" <<endl;
    // perform simulated anneling
    //doSimulatedAnneling(sampler, dbn, data, mism, sample, poseLow, MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf);
	if(regularize)
	{ 
		cout << "---------------- Starting use SAXS regularization method for sampling! --------------------"<< endl;
		//doSimulatedAnnelingSAXS_multiDomain_regularized(sampler, dbn, data, mism, sample, poseLow,linkers_info_array,MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf, decoyId,scoreType, combineScore);
		doSimulatedAnnelingSAXS_multiDomain_regularized_pro(sampler, dbn, data, mism, sample, poseLow,linkers_info_array,MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf, decoyId,scoreType, combineScore);
	}else{
		doSimulatedAnnelingSAXS_multiDomain(sampler, dbn, data, mism, sample, poseLow,linkers_info_array,MIN_FRG_LEN, MAX_FRG_LEN, Ti, Tf, decoyId,scoreType, combineScore);
    }
    //cout << "finish simulated annealing" <<endl;
    // populate the lowest scoring sample into current foldonSample
    for (int j = 0; j < sample.get_shape()[0]; j++) {
        foldonSample[i].set(j, 0, sample.get(j, 0));
        foldonSample[i].set(j, 1, sample.get(j, 1));
        foldonSample[i].set(j, 2, sample.get(j, 2));
        foldonSample[i].set(j, 3, sample.get(j, 3));
        foldonSample[i].set(j, 4, sample.get(j, 4));
        foldonSample[i].set(j, 5, sample.get(j, 5));
        foldonSample[i].set(j, 6, sample.get(j, 6));
    }
    pose2pdb(poseLow, pdbLow);
	
	// release the memory, add by jie
	mism.clear(); 
	
}

/*************************************************************************
 * Name        : doSimulatedAnneling
 * Purpose     : perform simulated anneling
 * Arguments   : SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp, CImgDisplay &dispLow, CImgDisplay &dispLive
 * Return Type : void
 *************************************************************************/
void doSimulatedAnneling(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp) {
    vector<poseInfo> pose;
    sample2pose(sample, pose);
    vector<pdbInfo> pdb;
    pose2pdb(pose, pdb);
    energyInfo energy = getWeightedEnergy(pdb, pose, dbn);
    
    // assign lowest energy pose to current pose
    poseLow = pose;
    vector<pdbInfo> pdbLow;
    pdbLow = pdb;
    energyInfo energyLow = getWeightedEnergy(pdbLow, poseLow, dbn);
    
    // define integer, discrete uniform distribution for position
    boost::uniform_int<> uniIntPos(0, sample.get_shape()[0]/2 - maxLength);
    // define a random variate generator using our base generator and distribution
    boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenPos(baseGen, uniIntPos);
    
    // define integer, discrete uniform distribution for length
    boost::uniform_int<> uniIntLen(minLength, maxLength);
    // define a random variate generator using our base generator and distribution
    boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenLen(baseGen, uniIntLen);
    
    // cycles and temperature schedule
    int outerCycles = sample.get_shape()[0]/2;
    int innerCycles = numCycles;
    //double initTemp = 1000.0;
    //double finalTemp = 298.0;
    double gamma = pow((finalTemp / initTemp),(double)(1.0 / (outerCycles * innerCycles)));
    double T = initTemp;
    
    // simulated anneling
    for (int i = 0; i < outerCycles; i++) {
        //sampler.set_seq_mismask(sample, mism);
        for (int j = 0; j < innerCycles; j++) {
            
            // randomly sample start and end
            int startPos = uniIntGenPos();
            int endPos;
            int len;
            while(true) {
                len = uniIntGenLen();
                endPos = startPos + len;
                if (endPos < sample.get_shape()[0]/2) {
                    break;
                }
            }
            sampler.set_start_end(startPos * 2, endPos * 2);
            MDArray<double> sampleNext = sampler.sample_next();
            
            for (int k = startPos; k < endPos; k++) {
                if (mism.get(k * 2, 5) == MOCAPY_OBSERVED) {
                	sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
                  sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
                }
            }
            vector<poseInfo> poseNext;
            sample2pose(sampleNext, poseNext);
            vector<pdbInfo> pdbNext;
            pose2pdb(poseNext, pdbNext);
            energyInfo energyNext = getWeightedEnergy(pdbNext, poseNext, dbn);
            // accept the move
            if (acceptMove(energy, energyNext, T)) {
                sampler.set_seq_mismask(sampleNext, mism);
                energy = energyNext;

                //if (((i + 1) * ( j + 1)) % innerCycles == 0) {
                    showEnergy(energyNext);
                //}
                
            }
            // reject the move
            else {
				//cout << "The current energy is: " << energyNext.total << " The previous energy is: " << energy.total << endl;
                sampler.undo();
            }
            // in case energy is lower than lowest-energy conformation so far
            if (energyNext.total < energyLow.total) {
                poseLow = poseNext;
                pdbLow = pdbNext;
                energyLow = energyNext;
                sample = sampleNext;
            }
            T = T * gamma;
        }
        
    }
}

/*************************************************************************
 * Name        : doSimulatedAnnelingSAXS
 * Purpose     : perform simulated anneling
 * Arguments   : SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp, CImgDisplay &dispLow, CImgDisplay &dispLive
 * Return Type : void
 *************************************************************************/
void doSimulatedAnnelingSAXS_multiDomain(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, vector<domLinkerInfo> linkers_info_array , int minLength, int maxLength, double initTemp,  double finalTemp, int decoyId, int scoreType, bool combineScore) {
    vector<poseInfo> pose;
    sample2pose(sample, pose);
    vector<pdbInfo> pdb;
    pose2pdb(pose, pdb);
    //energyInfo energy = getWeightedEnergy(pdb, pose, dbn);
    energyInfo energy = getWeightedEnergy_saxs_pro(pdb, pose, dbn,scoreType,combineScore);
    
    // assign lowest energy pose to current pose
    poseLow = pose;
    vector<pdbInfo> pdbLow;
    pdbLow = pdb;
    //energyInfo energyLow = getWeightedEnergy(pdbLow, poseLow, dbn);
    energyInfo energyLow = getWeightedEnergy_saxs_pro(pdbLow, poseLow, dbn,scoreType,combineScore);
    
	int linker_start = 0;
	int linker_end = 0;
	
	
	vector<double> Energy_total;
	vector<double> Energy_structure;
	vector<double> Energy_saxs;
	vector<int> Energy_selected;
	
	vector<double> SXS_chi_scores;
	double saxs_chi_global_initial = 100.0;
	SXS_chi_scores.push_back(saxs_chi_global_initial);
	int outerCycles = 1;
	// simulated anneling
	for (int i = 0; i < outerCycles; i++) { 
		//sampler.set_seq_mismask(sample, mism);
		clock_t OutCycle_start = clock();
		
		vector<sample_region> regions_array;
		int innerCycles=0;
		for(int linker=0;linker<linkers_info_array.size();linker++)
		{
			linker_start = linkers_info_array[linker].start;
			linker_end = linkers_info_array[linker].end;
			cout << "Start modeling linker " << linker << endl;
			cout << "linker start: " << linker_start << endl;
			cout << "linker end: " << linker_end << endl;
			// define integer, discrete uniform distribution for position, only for linkers, added by jie
			//boost::uniform_int<> uniIntPos(0, sample.get_shape()[0]/2);
			//boost::uniform_int<> uniIntPos(linker_start+1, linker_end-maxLength);  // should be set to linker_end-maxLength, otherwise, if the start is sampled at endpoint, the while loop won't out
			boost::uniform_int<> uniIntPos(linker_start+1, linker_end-1);  // should be set to linker_end-1 and linker_start+1, otherwise, if the start is sampled at endpoint, the while loop won't out

			//boost::uniform_int<> uniIntPos(linker_start, linker_end-minLength); // 120 ~ 140, more site to sample
			// define a random variate generator using our base generator and distribution
			boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenPos(baseGen, uniIntPos);
			
			// define integer, discrete uniform distribution for length
			boost::uniform_int<> uniIntLen(minLength, maxLength); //~ 3,6 
			// define a random variate generator using our base generator and distribution
			boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenLen(baseGen, uniIntLen);
			
			// cycles and temperature schedule
			//int outerCycles = sample.get_shape()[0]/2;
			
			int length_factor = ceil((linker_end-linker_start)/5);// set outer as 30
			if(length_factor < 1 )
			{
				length_factor = 1;
			}
			cout << "Setting length factor to " << length_factor << endl;
			innerCycles += numCycles*length_factor;//1000
			//double initTemp = 1000.0;
			//double finalTemp = 298.0;

			
			// set recording saxs score 
			int sample_num = outerCycles * innerCycles;
			total_sampling_cycle = sample_num;
			coolingDown = false;
			w_saxs_chi = w_saxs_chi_initial;
		

			for (int j = 0; j < innerCycles; j++) {
				// randomly sample start and end
				int startPos_tmp = uniIntGenPos();
		
				unsigned int startPos;
				unsigned int endPos;
				int len;
				int tmpPos;
				//right side
				while(true) {
					len = uniIntGenLen();
					endPos = startPos_tmp + len;
					if (endPos <= linker_end ) {
						break;
					}else{
						endPos = linker_end;
						break;
					}
					cout << "in sample " << "start: " << startPos_tmp << "end: "<< endPos<< endl;
				}
				//if(endPos > sample.get_shape()[0]/2)
				if(endPos > linker_end)
				{
					//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
					endPos = linker_end;
				}
				//if(endPos > sample.get_shape()[0]/2)
				if(startPos_tmp < linker_start)
				{
					//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
					startPos = linker_start;
				}else{
					startPos = startPos_tmp;
				}
				
				sample_region region;
				region.start = startPos;
				region.end = endPos;
				regions_array.push_back(region);
				//left side
				while(true) {
					len = uniIntGenLen();
					tmpPos = startPos_tmp-len;
					if (tmpPos >= linker_start) {
						endPos=startPos_tmp;
						startPos=tmpPos;
						break;
					}else{
						endPos=startPos_tmp;
						startPos=linker_start;
						break;
					}
					cout << "in sample " << "start: " << startPos << "end: "<< endPos<< endl;
				}
				//if(endPos > sample.get_shape()[0]/2)
				if(endPos > linker_end)
				{
					//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
					endPos = linker_end;
				}
				//if(endPos > sample.get_shape()[0]/2)
				if(startPos < linker_start)
				{
					//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
					startPos = linker_start;
				}
				
				region.start = startPos;
				region.end = endPos;
				regions_array.push_back(region);
			}
				
			cout << "Intermediate sample region size: " << regions_array.size() << " in region (" << linker_start << " ---- "<<linker_end <<")" << endl;
				
		}
		double gamma = pow((finalTemp / initTemp),(double)(1.0 / (outerCycles * innerCycles)));
		double T = initTemp;
		cout << "Total sample region size: " << regions_array.size() << endl;
		std::random_shuffle ( regions_array.begin(), regions_array.end() );
		for (int j = 0; j < regions_array.size(); j++) {
			int startPos = regions_array[j].start;
			int endPos = regions_array[j].end;
			
			if(j % 300 == 0)
			{
				cout << "Cycle " << j << " finished!" << endl;
			}
			sampler.set_start_end(startPos * 2, endPos * 2);
			MDArray<double> sampleNext = sampler.sample_next();
			
			for (unsigned int k = startPos; k < endPos; k++) { // even though won't have observed part, but still keep here 
				if (mism.get(k * 2, 5) == MOCAPY_OBSERVED) {\
					cout << "Intersting, won't happen here, why? The position is: " << k << endl;
					sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
				  sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
				}
			}
			//char sampleTempFile[1000];
			//sprintf(sampleTempFile, "%s%d_%d.txt", "samplefile_",i,j);
			//sample2file(sampleNext, sampleTempFile);
			
			
			vector<poseInfo> poseNext;
			sample2pose(sampleNext, poseNext);
			vector<pdbInfo> pdbNext;
			pose2pdb(poseNext, pdbNext);
			
			//char pdbTempFile[1000];
			//char pdbTempFile_pulchra[1000];
			//sprintf(pdbTempFile, "%s%d_%d.pdb", "samplefile_",i,j);
			//sprintf(pdbTempFile_pulchra, "%s%d_%d.rebuilt.pdb", "samplefile_",i,j);
			//writePdb(pdbNext, pdbTempFile);
			//runPulcha(pdbTempFile, pdbTempFile_pulchra);
			
			//energyInfo energyNext = getWeightedEnergy(pdbNext, poseNext, dbn);
			//energyInfo energyNext = getWeightedEnergy_saxs_unknown_inLoop(pdbNext, poseNext, dbn,SXS_chi_scores, proceeding_ratio);
			energyInfo energyNext = getWeightedEnergy_saxs_pro(pdbNext, poseNext, dbn,scoreType,combineScore);
			
			// record the energy 
			Energy_total.push_back(energyNext.total);
			Energy_structure.push_back(energyNext.structure_energy);
			Energy_saxs.push_back(energyNext.saxs_energy);
			SXS_chi_scores.push_back(energyNext.saxs_chi_global); //save the chi score
			
			// accept the move
			if (acceptMove(energy, energyNext, T)) {
				sampler.set_seq_mismask(sampleNext, mism);
				energy = energyNext;
				
				Energy_selected.push_back(1);

				//if (((i + 1) * ( j + 1)) % innerCycles == 0) {
					//showEnergy(energyNext);
				//}
				
				showEnergy_Print(energyNext,pdbNext, decoyId,startPos,endPos);
				cout<<endl; // make each accept energy showing on screen
				if(odir)
				{
					sprintf(pdbTempFile, "%s/%s%d_%d.pdb", outputdir, "samplefile_",i,j);
					sprintf(pdbTempFile_pulchra, "%s/%s%d_%d.rebuilt.pdb", outputdir, "samplefile_",i,j);
					sprintf(pdbTempFile_pose, "%s/%s%d_%d_pose.txt", outputdir, "samplefile_",i,j);
				}else{
					sprintf(pdbTempFile, "%s%d_%d.pdb", "samplefile_",i,j);
					sprintf(pdbTempFile_pulchra, "%s%d_%d.rebuilt.pdb", "samplefile_",i,j);
					sprintf(pdbTempFile_pose, "%s%d_%d_pose.txt", "samplefile_",i,j);
				}
				//writePdb(pdbNext, pdbTempFile);
				//runPulcha(pdbTempFile, pdbTempFile_pulchra);
				pdbString.clear();
				pdbString_pulchar.clear();
				Pdb2String(pdbNext, pdbString);				
				runPulcha3(pdbTempFile, pdbTempFile_pulchra,pdbString,pdbString_pulchar,1);	
				
				//output the pose information
				
				ofstream posestat(pdbTempFile_pose);
				if (!posestat.is_open()) {
					cout << "Error! folding statistics file can not open " << pdbTempFile_pose << endl;
					exit(0);
				}
				char bufposeStat[1000];
				sprintf(bufposeStat,"id\taa\tss\tbca\ttao\ttheta\tbsc\tphi\tdelta");
				posestat << bufposeStat << endl;
				for(int h=0;h<poseNext.size();h++)
				{
					// buffer statistics
					char bufposeStat[1000]; 
					sprintf(bufposeStat,"%i\t%i\t%i\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f", poseNext[h].id, poseNext[h].aa,poseNext[h].ss,poseNext[h].bca, poseNext[h].tao, poseNext[h].theta,poseNext[h].bsc,poseNext[h].phi,poseNext[h].delta);
				
					posestat << bufposeStat << endl;
				}
				
				
			}
			// reject the move
			else {
				Energy_selected.push_back(0);
				sampler.undo();
			}
			// in case energy is lower than lowest-energy conformation so far
			if (energyNext.total < energyLow.total) {
				poseLow = poseNext;
				pdbLow = pdbNext;
				energyLow = energyNext;
				sample = sampleNext;
			}
			T = T * gamma;
				
		}
		clock_t OutCycle_end = clock();
		double elapsecycle_secs1 = double(OutCycle_end - OutCycle_start) / CLOCKS_PER_SEC;
		cout << "UniCon3D.io: Outcycle " << i <<" finished within "<< elapsecycle_secs1 << " sec!" << endl << endl;
		
		
	
	}
	

	
	if(Energy_total.size() != Energy_structure.size() and Energy_total.size() != Energy_saxs.size()  and Energy_total.size() != Energy_selected.size() and Energy_total.size() != SXS_chi_scores.size())
	{
		cout << "Error happens!" <<endl;
	}
	if(odir)
	{
		
		sprintf(simulationStatFile, "%s/%s_w%d_simulationStats.txt", outputdir,jobId,w_saxs_chi_initial);
	}else{
		sprintf(simulationStatFile, "%s_w%d_simulationStats.txt", jobId,w_saxs_chi_initial);
	}
	
	ofstream simustat(simulationStatFile);
    if (!simustat.is_open()) {
        cout << "Error! folding statistics file can not open " << simulationStatFile << endl;
        exit(0);
    }
	char bufsimStat[1000];
	sprintf(bufsimStat,"Epoch\tAccept\tEnergy_structure\tEnergy_saxs\tEnergy_total\tChi-score");
	simustat << bufsimStat << endl;
	for(int h=0;h<Energy_total.size();h++)
	{
		// buffer statistics
		char bufsimStat[1000];
		sprintf(bufsimStat,"%i\t%i\t%8.3f\t%8.3f\t%8.3f\t%8.3f", h, Energy_selected[h], Energy_structure[h],Energy_saxs[h],Energy_total[h],SXS_chi_scores[h]);
		simustat << bufsimStat << endl;
	}
	
}

/*************************************************************************
 * Name        : doSimulatedAnnelingSAXS
 * Purpose     : perform simulated anneling
 * Arguments   : SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, int minLength, int maxLength, double initTemp,  double finalTemp, CImgDisplay &dispLow, CImgDisplay &dispLive
 * Return Type : void
 *************************************************************************/
void doSimulatedAnnelingSAXS_multiDomain_regularized_pro(SampleInfEngineHMM &sampler, DBN &dbn, Sequence &data, MDArray<eMISMASK> &mism, MDArray<double> &sample, vector<poseInfo> &poseLow, vector<domLinkerInfo> linkers_info_array , int minLength, int maxLength, double initTemp,  double finalTemp, int decoyId, int scoreType, bool combineScore) {
    vector<poseInfo> pose;
    sample2pose(sample, pose);
    vector<pdbInfo> pdb;
    pose2pdb(pose, pdb);
    //energyInfo energy = getWeightedEnergy(pdb, pose, dbn);
    //energyInfo energy = getWeightedEnergy_saxs_pro(pdb, pose, dbn,scoreType,combineScore);
    w_saxs_chi = w_saxs_chi_initial;
	w_saxs_chi_penalty = w_saxs_chi_penalty_initial;
	
	w_saxs_KL = w_saxs_KL_initial;	
	w_saxs_KL_penalty = w_saxs_KL_penalty_initial;
	
	
	w_saxs_score2 = w_saxs_score2_initial;	
	w_saxs_score2_penalty = w_saxs_score2_penalty_initial;
	
	
	w_saxs_RG_normalize = w_saxs_RG_normalize_initial;	
	w_saxs_RG_normalize_penalty = w_saxs_RG_normalize_penalty_initial;
    energyInfo energy = getWeightedEnergy_saxs_regularize_pro(pdb, pose, dbn,scoreType,combineScore);
	cout << "The initial energy is: " << endl;
    showEnergy_Print_regularize(energy,pdb, 0,0,0);
    // assign lowest energy pose to current pose
    poseLow = pose;
    vector<pdbInfo> pdbLow;
    pdbLow = pdb;
    //energyInfo energyLow = getWeightedEnergy(pdbLow, poseLow, dbn);
    //energyInfo energyLow = getWeightedEnergy_saxs_pro(pdbLow, poseLow, dbn,scoreType,combineScore);
    energyInfo energyLow = getWeightedEnergy_saxs_regularize_pro(pdbLow, poseLow, dbn,scoreType,combineScore);
    
	int linker_start = 0;
	int linker_end = 0;
	
	
	vector<double> Energy_total;
	vector<double> Energy_structure;
	vector<double> Energy_saxs;
	vector<int> Energy_selected;
	
	vector<double> SXS_chi_scores;
	double saxs_chi_global_initial = 100.0;
	SXS_chi_scores.push_back(saxs_chi_global_initial);
	
		

		// simulated anneling
		// cycles and temperature schedule
		//int outerCycles = sample.get_shape()[0]/2;
		int outerCycles = 1;// set outer as 30
		int innerCycles = numCycles;//1000
		//double initTemp = 1000.0;
		//double finalTemp = 298.0;
		double gamma = pow((finalTemp / initTemp),(double)(1.0 / (outerCycles * innerCycles*2)));// *2 because in my sampling, each time generate 2 sample regions
		double T = initTemp;
		
		double gamma_w_chi = pow((w_saxs_chi_final / w_saxs_chi_initial),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_chi = w_saxs_chi_initial;
		double gamma_w_chi_penalty = pow((w_saxs_chi_penalty_final/w_saxs_chi_penalty_initial ),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_chi_penalty = w_saxs_chi_penalty_initial;
		
		double gamma_w_KL = pow((w_saxs_KL_final / w_saxs_KL_initial),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_KL = w_saxs_KL_initial;	
		double gamma_w_KL_penalty = pow((w_saxs_KL_penalty_final/w_saxs_KL_penalty_initial ),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_KL_penalty = w_saxs_KL_penalty_initial;
		
		double gamma_w_score2 = pow((w_saxs_score2_final / w_saxs_score2_initial),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_score2 = w_saxs_score2_initial;	
		double gamma_w_score2_penalty = pow((w_saxs_score2_penalty_final/w_saxs_score2_penalty_initial ),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_score2_penalty = w_saxs_score2_penalty_initial;
		
		double gamma_w_RG_normalize = pow((w_saxs_RG_normalize_final / w_saxs_RG_normalize_initial),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_RG_normalize = w_saxs_RG_normalize_initial;	
		double gamma_w_RG_normalize_penalty = pow((w_saxs_RG_normalize_penalty_final/w_saxs_RG_normalize_penalty_initial ),(double)(1.0 / (outerCycles * innerCycles*2)));
		w_saxs_RG_normalize_penalty = w_saxs_RG_normalize_penalty_initial;
		
		int epoch=0;
		MDArray<double> sampleNext;
		for (int i = 0; i < outerCycles+1; i++) { // plus additional cycle for final approximation 
			
			clock_t OutCycle_start = clock();
			vector<sample_region> regions_array;
			
			for(int linker=0;linker<linkers_info_array.size();linker++)
			{
				linker_start = linkers_info_array[linker].start;
				linker_end = linkers_info_array[linker].end;
				cout << "Start modeling linker " << linker << endl;
				cout << "linker start: " << linker_start << endl;
				cout << "linker end: " << linker_end << endl;
				// define integer, discrete uniform distribution for position, only for linkers, added by jie
				//boost::uniform_int<> uniIntPos(0, sample.get_shape()[0]/2);
				//boost::uniform_int<> uniIntPos(linker_start+1, linker_end-maxLength);  // should be set to linker_end-maxLength, otherwise, if the start is sampled at endpoint, the while loop won't out
				boost::uniform_int<> uniIntPos(linker_start+1, linker_end-1);  // should be set to linker_end-1 and linker_start+1, otherwise, if the start is sampled at endpoint, the while loop won't out
				//boost::uniform_int<> uniIntPos(linker_start, linker_end-minLength); // 120 ~ 140, more site to sample
				// define a random variate generator using our base generator and distribution
				boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenPos(baseGen, uniIntPos);
				
				// define integer, discrete uniform distribution for length
				boost::uniform_int<> uniIntLen(minLength, maxLength); //~ 3,6 
				// define a random variate generator using our base generator and distribution
				boost::variate_generator<boost::minstd_rand&, boost::uniform_int<> > uniIntGenLen(baseGen, uniIntLen);
				


				// set recording saxs score 
				int sample_num = outerCycles * innerCycles;
				total_sampling_cycle = sample_num;
			
				for (int j = 0; j < innerCycles; j++) {
					int startPos_tmp = uniIntGenPos();
			
					unsigned int startPos;
					unsigned int endPos;
					int len;
					int tmpPos;
					//right side
					while(true) {
						len = uniIntGenLen();
						endPos = startPos_tmp + len;
						if (endPos <= linker_end ) {
							break;
						}else{
							endPos = linker_end;
							break;
						}
						cout << "in sample " << "start: " << startPos_tmp << "end: "<< endPos<< endl;
					}
					//if(endPos > sample.get_shape()[0]/2)
					if(endPos > linker_end)
					{
						//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
						endPos = linker_end;
					}
					//if(endPos > sample.get_shape()[0]/2)
					if(startPos_tmp < linker_start)
					{
						//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
						startPos = linker_start;
					}else{
						startPos = startPos_tmp;
					}
					
					sample_region region;
					region.start = startPos;
					region.end = endPos;
					regions_array.push_back(region);
					//left side
					while(true) {
						len = uniIntGenLen();
						tmpPos = startPos_tmp-len;
						if (tmpPos >= linker_start) {
							endPos=startPos_tmp;
							startPos=tmpPos;
							break;
						}else{
							endPos=startPos_tmp;
							startPos=linker_start;
							break;
						}
						cout << "in sample " << "start: " << startPos << "end: "<< endPos<< endl;
					}
					//if(endPos > sample.get_shape()[0]/2)
					if(endPos > linker_end)
					{
						//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
						endPos = linker_end;
					}
					//if(endPos > sample.get_shape()[0]/2)
					if(startPos < linker_start)
					{
						//cout << " The endpos is exceeding the margin (" << endPos << ")" << endl;
						startPos = linker_start;
					}
					
					region.start = startPos;
					region.end = endPos;
					regions_array.push_back(region);
				}
				
				cout << "Intermediate sample region size: " << regions_array.size() << " in region (" << linker_start << " ---- "<<linker_end <<")" << endl;
				
			}
			cout << "Total sample region size: " << regions_array.size() << endl;
			std::random_shuffle ( regions_array.begin(), regions_array.end() );
			for (int j = 0; j < regions_array.size(); j++) {
				
				// randomly sample start and end
				epoch++;
				//if((epoch % 100) ==0)
				//{
				//	cout << endl;
				//	cout << "epoch: "<< epoch  << endl;
				//	cout << "T: "<< T  << endl;
				//	cout << "w_saxs_chi: "<< w_saxs_chi  << endl;
				//	cout<< "w_saxs_chi_penalty: " <<  w_saxs_chi_penalty << endl;
				//	cout << "gamma: "<< gamma  << " gamma_w: " <<  gamma_w << " gamma_w_penalty: " <<  gamma_w_penalty << endl;
				//}
				
				int startPos = regions_array[j].start;
				int endPos = regions_array[j].end;
				//cout << "In program: " << "start: " << startPos << "end: "<< endPos<< endl;
				sampler.set_start_end(startPos * 2, endPos * 2);
				sampleNext.clear();
				sampleNext = sampler.sample_next();
				
				for (unsigned int k = startPos; k < endPos; k++) { // even though won't have observed part, but still keep here 
				    if (mism.get(k * 2, 5) == MOCAPY_OBSERVED) {\
						cout << "Intersting, won't happen here, why? The position is: " << k << endl;
				    	sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
				      sampleNext.set(k * 2, 5, DEG2RAD * (RAD2DEG * sampleNext.get(k * 2, 5) - 1.0 + uniDblGen() * 2.0));
				    }
				}
				//char sampleTempFile[1000];
				//sprintf(sampleTempFile, "%s%d_%d.txt", "samplefile_",i,j);
				//sample2file(sampleNext, sampleTempFile);
				
				
				vector<poseInfo> poseNext;
				sample2pose(sampleNext, poseNext);
				vector<pdbInfo> pdbNext;
				pose2pdb(poseNext, pdbNext);
				
				//char pdbTempFile[1000];
				//char pdbTempFile_pulchra[1000];
				//sprintf(pdbTempFile, "%s%d_%d.pdb", "samplefile_",i,j);
				//sprintf(pdbTempFile_pulchra, "%s%d_%d.rebuilt.pdb", "samplefile_",i,j);
				//writePdb(pdbNext, pdbTempFile);
				//runPulcha(pdbTempFile, pdbTempFile_pulchra);
				
				//energyInfo energyNext = getWeightedEnergy(pdbNext, poseNext, dbn);
				//energyInfo energyNext = getWeightedEnergy_saxs_unknown_inLoop(pdbNext, poseNext, dbn,SXS_chi_scores, proceeding_ratio);
				//energyInfo energyNext = getWeightedEnergy_saxs_pro(pdbNext, poseNext, dbn,scoreType,combineScore);
				energyInfo energyNext = getWeightedEnergy_saxs_regularize_pro(pdbNext, poseNext, dbn,scoreType,combineScore);
				
				// record the energy 
				Energy_total.push_back(energyNext.total);
				Energy_structure.push_back(energyNext.structure_energy);
				Energy_saxs.push_back(energyNext.saxs_energy);
				SXS_chi_scores.push_back(energyNext.saxs_chi_global); //save the chi score
				//showEnergy_Print_regularize(energyNext,pdbNext, decoyId,startPos,endPos);
				// accept the move
				if (acceptMove(energy, energyNext, T)) {
					sampler.set_seq_mismask(sampleNext, mism);
					energy = energyNext;
					
					Energy_selected.push_back(1);

					//if (((i + 1) * ( j + 1)) % innerCycles == 0) {
						//showEnergy(energyNext);
					//}
					
					// assign the global chi minima
					if(energyNext.saxs_chi_global <= saxs_chi_global_minima)
					{
						saxs_chi_global_minima = energyNext.saxs_chi_global;
					}
					
					
					showEnergy_Print_regularize(energyNext,pdbNext, decoyId,startPos,endPos);
					cout<<endl; // make each accept energy showing on screen
					if(odir)
					{
						sprintf(pdbTempFile, "%s/%s%d_%d.pdb", outputdir, "samplefile_",i,j);
						sprintf(pdbTempFile_pulchra, "%s/%s%d_%d.rebuilt.pdb", outputdir, "samplefile_",i,j);
						sprintf(pdbTempFile_pose, "%s/%s%d_%d_pose.txt", outputdir, "samplefile_",i,j);
					}else{
						sprintf(pdbTempFile, "%s%d_%d.pdb", "samplefile_",i,j);
						sprintf(pdbTempFile_pulchra, "%s%d_%d.rebuilt.pdb", "samplefile_",i,j);
						sprintf(pdbTempFile_pose, "%s%d_%d_pose.txt", "samplefile_",i,j);
					}
					//writePdb(pdbNext, pdbTempFile);
					//runPulcha(pdbTempFile, pdbTempFile_pulchra);	
					pdbString.clear();
					pdbString_pulchar.clear();
					Pdb2String(pdbNext, pdbString);				
					runPulcha3(pdbTempFile, pdbTempFile_pulchra,pdbString,pdbString_pulchar,1);	
					
					//output the pose information
					
					ofstream posestat(pdbTempFile_pose);
					if (!posestat.is_open()) {
						cout << "Error! folding statistics file can not open " << pdbTempFile_pose << endl;
						exit(0);
					}
					sprintf(bufposeStat,"id\taa\tss\tbca\ttao\ttheta\tbsc\tphi\tdelta");
					posestat << bufposeStat << endl;
					for(int h=0;h<poseNext.size();h++)
					{
						// buffer statistics
						sprintf(bufposeStat,"%i\t%i\t%i\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f", poseNext[h].id, poseNext[h].aa,poseNext[h].ss,poseNext[h].bca, poseNext[h].tao, poseNext[h].theta,poseNext[h].bsc,poseNext[h].phi,poseNext[h].delta);
					
						posestat << bufposeStat << endl;
					}
					
					
				}
				// reject the move
				else {
					Energy_selected.push_back(0);
					sampler.undo();
				}
				// in case energy is lower than lowest-energy conformation so far
				if (energyNext.total < energyLow.total) {
					poseLow = poseNext;
					pdbLow = pdbNext;
					energyLow = energyNext;
					sample = sampleNext;
				}
				T = T * gamma;
				
				// update chi
				if(w_saxs_chi>w_saxs_chi_final)
				{
					w_saxs_chi = w_saxs_chi*gamma_w_chi;
				}else{
					w_saxs_chi = w_saxs_chi_final;
				}
				if(w_saxs_chi_penalty<w_saxs_chi_penalty_final)
				{
					w_saxs_chi_penalty = w_saxs_chi_penalty*gamma_w_chi_penalty;
				}else{
					w_saxs_chi_penalty = w_saxs_chi_penalty_final;
				}
				// update KL
				if(w_saxs_KL>w_saxs_KL_final)
				{
					w_saxs_KL = w_saxs_KL*gamma_w_KL;
				}else{
					w_saxs_KL = w_saxs_KL_final;
				}
				if(w_saxs_KL_penalty<w_saxs_KL_penalty_final)
				{
					w_saxs_KL_penalty = w_saxs_KL_penalty*gamma_w_KL_penalty;
				}else{
					w_saxs_KL_penalty = w_saxs_KL_penalty_final;
				}
				
				// update score2
				if(w_saxs_score2>w_saxs_score2_final)
				{
					w_saxs_score2 = w_saxs_score2*gamma_w_score2;
				}else{
					w_saxs_score2 = w_saxs_score2_final;
				}
				if(w_saxs_score2_penalty<w_saxs_score2_penalty_final)
				{
					w_saxs_score2_penalty = w_saxs_score2_penalty*gamma_w_score2_penalty;
				}else{
					w_saxs_score2_penalty = w_saxs_score2_penalty_final;
				}
				
				// update RG_normalize
				if(w_saxs_RG_normalize>w_saxs_RG_normalize_final)
				{
					w_saxs_RG_normalize = w_saxs_RG_normalize*gamma_w_RG_normalize;
				}else{
					w_saxs_RG_normalize = w_saxs_RG_normalize_final;
				}
				if(w_saxs_RG_normalize_penalty<w_saxs_RG_normalize_penalty_final)
				{
					w_saxs_RG_normalize_penalty = w_saxs_RG_normalize_penalty*gamma_w_RG_normalize_penalty;
				}else{
					w_saxs_RG_normalize_penalty = w_saxs_RG_normalize_penalty_final;
				}
				
				
				// clear the memory
				poseNext.clear();
				pdbNext.clear();
				vector<poseInfo>().swap(poseNext);
				vector<pdbInfo>().swap(pdbNext);
					
			}
			clock_t OutCycle_end = clock();
			double elapsecycle_secs1 = double(OutCycle_end - OutCycle_start) / CLOCKS_PER_SEC;
			cout << "UniCon3D.io: Outcycle " << i <<" finished within "<< elapsecycle_secs1 << " sec!" << endl << endl;
			
			// clear the memory
			vector<sample_region>().swap(regions_array);
			
		}
	
	

	
	if(Energy_total.size() != Energy_structure.size() and Energy_total.size() != Energy_saxs.size()  and Energy_total.size() != Energy_selected.size()  and Energy_total.size() != SXS_chi_scores.size())
	{
		cout << "Error happens!" <<endl;
	}
	if(odir)
	{
		
		sprintf(simulationStatFile, "%s/%s_w%d_simulationStats.txt", outputdir,jobId,w_saxs_chi_initial);
	}else{
		sprintf(simulationStatFile, "%s_w%d_simulationStats.txt", jobId,w_saxs_chi_initial);
	}
	
	ofstream simustat(simulationStatFile);
    if (!simustat.is_open()) {
        cout << "Error! folding statistics file can not open " << simulationStatFile << endl;
        exit(0);
    }
	
	sprintf(bufsimStat,"Epoch\tAccept\tEnergy_structure\tEnergy_saxs\tEnergy_total\tChi-score");
	simustat << bufsimStat << endl;
	for(int h=0;h<Energy_total.size();h++)
	{
		// buffer statistics
		sprintf(bufsimStat,"%i\t%i\t%8.3f\t%8.3f\t%8.3f\t%8.3f", h, Energy_selected[h], Energy_structure[h],Energy_saxs[h],Energy_total[h],SXS_chi_scores[h]);
		simustat << bufsimStat << endl;
	}
	
	//release memory
	pose.clear();
	pdb.clear();
	pdbLow.clear();
	vector<poseInfo>().swap(pose);
	vector<pdbInfo>().swap(pdb);
	vector<pdbInfo>().swap(pdbLow);
	
	
	vector<double>().swap(Energy_total);
	vector<double>().swap(Energy_structure);
	vector<double>().swap(Energy_saxs);
	vector<int>().swap(Energy_selected);
	
	vector<double>().swap(SXS_chi_scores);
	
	
}

/*************************************************************************
 * Name        : getWeightedEnergy_saxs_unknown
 * Purpose     : gets total weighted energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn
 * Return Type : energyInfo
 *************************************************************************/
energyInfo getWeightedEnergy_saxs_unknown(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn) {
    energyInfo e;
    e.sc_sc = getScScEnergy(pdb, pose, 1.0);
    e.sc_bb = getScBbEnergy(pdb, pose, 1.0);
    e.bb_bb = getBbBbEnergy(pdb, pose, 1.0);
    e.ri_rj = getRiRjEnergy(pdb, pose, 1.0);
	
	int idrand = rand();
	

	if(odir)
	{
		sprintf(pdbTempFile, "%s/%s%d.pdb", outputdir, "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_prefix, "%s/%s%d", outputdir, "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_pulchra, "%s/%s.rebuilt.pdb", outputdir,pdbTempFile_prefix);

	}else{
		sprintf(pdbTempFile, "%s%d.pdb", "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_prefix, "%s%d", "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_pulchra, "%s.rebuilt.pdb",pdbTempFile_prefix);

	}	
	
	//writePdb(pdb, pdbTempFile);
	pdbString.clear();
	pdbString_pulchar.clear();
	Pdb2String(pdb, pdbString); // convert pdb file to string
	runPulcha3(pdbTempFile, pdbTempFile_pulchra,pdbString,pdbString_pulchar,1);

  //cout << "Output pdb "<<endl;
   //std::string line;
  //std::istringstream split(pdbString);
  //while(std::getline(split, line, ';')) {
  // std::cout << line << '\n';
  //}  //cout << "Output pdb "<<endl;
  // std::string line2;
  //std::istringstream split2(pdbString_pulchar);
 // while(std::getline(split2, line2, ';')) {
  // std::cout << line2 << '\n';
  //}
	
	if(CAfitting)
	{
		//e.saxs_chi_global = getChiEnergy(pdbTempFile_pulchra,saxsFile,1.0);  
		//e.saxs_chi_global = getChiEnergy_inside(pdbTempFile,saxsFile,1.0);  	
		e.saxs_chi_global = getChiEnergy_inside2(pdbTempFile,saxsFile,pdbString_pulchar,1.0);  	
	}else{
		// run pulcha on global structure
		//runPulcha(pdbTempFile, pdbTempFile_pulchra);
		//e.saxs_chi_global = getChiEnergy(pdbTempFile_pulchra,saxsFile,1.0);  
		//e.saxs_chi_global = getChiEnergy_inside(pdbTempFile_pulchra,saxsFile,1.0);  		
        //cout << "here" << endl;			
		e.saxs_chi_global = getChiEnergy_inside2(pdbTempFile_pulchra,saxsFile,pdbString_pulchar,1.0);  				
                //cout << "finish" << endl;	
	}
	
	// Remove temporary file 
	//removeTempFile(pdbTempFile_prefix);
 
	
	if((e.saxs_chi_global - 0.0) < 0.001 )
	{ 
		cout << "Error: saxs score is incorrect: " << e.saxs_chi_global << endl; 
		exit(0);
	}
	//cout << "e.saxs_chi: " << e.saxs_chi << endl; 
    //e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj;
    e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_chi * e.saxs_chi_global; 
    return e;
}


/*************************************************************************
 * Name        : getWeightedEnergy_saxs_unknown
 * Purpose     : gets total weighted energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn
 * Return Type : energyInfo
 *************************************************************************/
energyInfo getWeightedEnergy_saxs_pro(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn, int scoreType, bool combineScore) {
    energyInfo e;
    e.sc_sc = getScScEnergy(pdb, pose, 1.0);
    e.sc_bb = getScBbEnergy(pdb, pose, 1.0);
    e.bb_bb = getBbBbEnergy(pdb, pose, 1.0);
    e.ri_rj = getRiRjEnergy(pdb, pose, 1.0);
	
	int idrand = rand();
	
	
	if(odir)
	{
		sprintf(pdbTempFile, "%s/%s%d.pdb", outputdir, "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_prefix, "%s/%s%d", outputdir, "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_pulchra, "%s/%s.rebuilt.pdb", outputdir,pdbTempFile_prefix);

	}else{
		sprintf(pdbTempFile, "%s%d.pdb", "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_prefix, "%s%d", "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_pulchra, "%s.rebuilt.pdb",pdbTempFile_prefix);

	}
	//writePdb(pdb, pdbTempFile);
	pdbString.clear();
	pdbString_pulchar.clear();
	Pdb2String(pdb, pdbString); // convert pdb file to string
	runPulcha3(pdbTempFile, pdbTempFile_pulchra,pdbString,pdbString_pulchar,1);

	
	
	// sometimes sample failed, pdb contains nan, ignore 
	if (pdbString_pulchar.find("nan") != std::string::npos) {
		std::cout << "found nan, ignore!" << '\n';
		e.saxs_chi_global = 100.0;  	
		e.saxs_KL =100.0;  	
		e.saxs_score2 = 100.0;  	
		e.saxs_RG_normalize = 100.0;  
	}else{
		saxsInfo saxs_score=getChiEnergy_inside_pro(pdbTempFile,saxsFile,pdbString_pulchar,1.0);  
		e.saxs_chi_global = saxs_score.chi_score;  	
		e.saxs_KL =saxs_score.KL;  	
		e.saxs_score2 = saxs_score.fun_score2;  	
		e.saxs_RG_normalize = saxs_score.RGdiff_normalize;  
			
	}
	
	
	
	if((e.saxs_chi_global - 0.0) < 0.001 )
	{ 
		cout << "Error: saxs score is incorrect: " << e.saxs_chi_global << endl; 
		exit(0);
	}
	
	
	if(combineScore)
	{
		e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_chi * e.saxs_chi_global+ w_saxs_KL * e.saxs_KL+ w_saxs_score2 * e.saxs_score2+ w_saxs_RG_normalize * e.saxs_RG_normalize; 
		e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
		e.saxs_energy =  w_saxs_chi * e.saxs_chi_global+ w_saxs_KL * e.saxs_KL+ w_saxs_score2 * e.saxs_score2+ w_saxs_RG_normalize * e.saxs_RG_normalize; 
	}else{
		if(scoreType == 1)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_chi * e.saxs_chi_global; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_chi * e.saxs_chi_global; 
		}else if(scoreType == 2)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_KL * e.saxs_KL; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_KL * e.saxs_KL; 
		}else if(scoreType == 3)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_score2 * e.saxs_score2; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_score2 * e.saxs_score2; 
		}else if(scoreType == 4)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_RG_normalize * e.saxs_RG_normalize; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_RG_normalize * e.saxs_RG_normalize; 
		}
	}
	e.saxs_penalty=0;
	
	// Remove temporary file 
	//removeTempFile(pdbTempFile_prefix);
 
	
	//cout << "e.saxs_chi: " << e.saxs_chi << endl; 
    //e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj;
    
    return e;
}


/*************************************************************************
 * Name        : getWeightedEnergy_saxs_unknown
 * Purpose     : gets total weighted energy
 * Arguments   : vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn
 * Return Type : energyInfo
 *************************************************************************/
energyInfo getWeightedEnergy_saxs_regularize_pro(vector<pdbInfo> &pdb, vector<poseInfo> &pose, DBN &dbn, int scoreType, bool combineScore) {
    energyInfo e;
    e.sc_sc = getScScEnergy(pdb, pose, 1.0);
    e.sc_bb = getScBbEnergy(pdb, pose, 1.0);
    e.bb_bb = getBbBbEnergy(pdb, pose, 1.0);
    e.ri_rj = getRiRjEnergy(pdb, pose, 1.0);
	
	int idrand = rand();
	
	
	if(odir)
	{
		sprintf(pdbTempFile, "%s/%s%d.pdb", outputdir, "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_prefix, "%s/%s%d", outputdir, "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_pulchra, "%s/%s.rebuilt.pdb", outputdir,pdbTempFile_prefix);

	}else{
		sprintf(pdbTempFile, "%s%d.pdb", "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_prefix, "%s%d", "GlobalFoldon_",idrand);
		sprintf(pdbTempFile_pulchra, "%s.rebuilt.pdb",pdbTempFile_prefix);

	}
	//writePdb(pdb, pdbTempFile);
	pdbString.clear();
	pdbString_pulchar.clear();
	Pdb2String(pdb, pdbString); // convert pdb file to string
	runPulcha3(pdbTempFile, pdbTempFile_pulchra,pdbString,pdbString_pulchar,1);

	
	saxsInfo saxs_score=getChiEnergy_inside_pro(pdbTempFile,saxsFile,pdbString_pulchar,1.0);
	e.saxs_chi_global = saxs_score.chi_score;  	
	e.saxs_KL =saxs_score.KL;  	
	e.saxs_score2 = saxs_score.fun_score2;  	
	e.saxs_RG_normalize = saxs_score.RGdiff_normalize;  
	
	
	
	double deviation = abs(e.saxs_chi_global-saxs_chi_global_minima);
	e.saxs_penalty = deviation*deviation; 

 
 
	if((e.saxs_chi_global - 0.0) < 0.001 )
	{ 
		cout << "Error: saxs score is incorrect: " << e.saxs_chi_global << endl; 
		exit(0);
	}
	
	
	if(combineScore)
	{
		e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_chi * e.saxs_chi_global+ w_saxs_KL * e.saxs_KL+ w_saxs_score2 * e.saxs_score2+ w_saxs_RG_normalize * e.saxs_RG_normalize + w_saxs_chi_penalty*e.saxs_penalty; 
		e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
		e.saxs_energy =  w_saxs_chi * e.saxs_chi_global+ w_saxs_KL * e.saxs_KL+ w_saxs_score2 * e.saxs_score2+ w_saxs_RG_normalize * e.saxs_RG_normalize + w_saxs_chi_penalty*e.saxs_penalty; 
	}else{
		if(scoreType == 1)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_chi * e.saxs_chi_global + w_saxs_chi_penalty*e.saxs_penalty; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_chi * e.saxs_chi_global + w_saxs_chi_penalty*e.saxs_penalty; 
		}else if(scoreType == 2)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_KL * e.saxs_KL + w_saxs_chi_penalty*e.saxs_penalty; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_KL * e.saxs_KL + w_saxs_chi_penalty*e.saxs_penalty; 
		}else if(scoreType == 3)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_score2 * e.saxs_score2 + w_saxs_chi_penalty*e.saxs_penalty; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_score2 * e.saxs_score2 + w_saxs_chi_penalty*e.saxs_penalty; 
		}else if(scoreType == 4)
		{
			e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj + w_saxs_RG_normalize * e.saxs_RG_normalize; 
			e.structure_energy = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj; 
			e.saxs_energy = w_saxs_RG_normalize * e.saxs_RG_normalize + w_saxs_chi_penalty*e.saxs_penalty; 
			
		}
		
	}
	
	// Remove temporary file 
	//removeTempFile(pdbTempFile_prefix);
 
	
	//cout << "e.saxs_chi: " << e.saxs_chi << endl; 
    //e.total = w_sc_sc * e.sc_sc + w_sc_bb * e.sc_bb + w_bb_bb * e.bb_bb + w_ri_rj * e.ri_rj;
    
    return e;
}




/*************************************************************************
 * Name        : getChiEnergy
 * Purpose     : gets the chi-score energy
 * Arguments   : char *filename, char *SAXSfilename, double weight
 * Return Type : double
 *************************************************************************/
double getChiEnergy(char *filename, char *SAXSfilename, double weight) {
    double energy = 0.0;
	string line;
	FILE * stream;
    const int max_buffer = 256;
    char buffer[max_buffer];
	string Chi ("Chi = ");
	string C1 ("c1 =");
	//char foxlog[1000];
	//int idrand = rand();
	//sprintf(foxlog, "%s%d.log", "fox",idrand);
	string cmd="";
	if(CAfitting)
	{
		 cmd = "/home/casp11/CASP12/SAXS-CASP12/tools/FoXS/foxs -q 0.3 -s 500 -r   ";
	}else{
		 cmd = "/home/casp11/CASP12/SAXS-CASP12/tools/FoXS/foxs -q 0.3 -s 500   ";
	}
	
	cmd.append(filename);
	cmd.append("  ");
	cmd.append(SAXSfilename);
	cmd.append(" 2>&1 ");
	//cmd.append(foxlog);
    //cout << "command: " << cmd << endl;
	
	
    stream = popen(cmd.c_str(), "r");
    if (stream) {
		//while (!feof(stream))
		//{
			while (fgets(buffer, max_buffer, stream) != NULL){
				line="";
				// assign the values in the array to a string
				line += buffer;
				//cout << "Chi-score: " <<  line << endl;
				if(line.find(Chi,0) != string::npos)
				{
					if( line.compare(line.find(Chi), 6, "Chi = ") == 0) {
						//cout << "line: " << line << endl;
						string chi_score = line.substr(line.find(Chi)+6, line.find(C1)-line.find(Chi)-6);
						chi_score.erase(remove_if(chi_score.begin(), chi_score.end(), (int(*)(int))isspace), chi_score.end());
						energy = atof(chi_score.c_str());
						//cout << "Chi-score: " <<  energy << endl;
					} 
					//fflush(stdout);
				} 
			} 
					
		//}
	}
	pclose(stream);	

/*
	stream = popen(cmd.c_str(), "r");
	if (stream) {
		ifstream fin ("fox.log");
		if (fin.is_open()) {
			while ( fin.good() ) {
				getline(fin, line);
				if(line.length() != 0)
				{
					if(line.find(Chi,0) != string::npos)
					{
						if( line.compare(line.find(Chi), 6, "Chi = ") == 0) {
							cout << "line: " << line << endl;
							string chi_score = line.substr(line.find(Chi)+6, line.find(C1)-line.find(Chi)-6);
							chi_score.erase(remove_if(chi_score.begin(), chi_score.end(), (int(*)(int))isspace), chi_score.end());
							energy = atof(chi_score.c_str());
							cout << "Chi-score: " <<  energy << endl;
						}
					}
				}
			}
		}
		fin.close();
	}else{
		cout << "Failed to run: " << cmd << endl;
		exit(0);
		
	}
*/	
    return energy * weight;

}



/*************************************************************************
 * Name        : getChiEnergy inside program
 * Purpose     : gets the chi-score energy
 * Arguments   : char *filename, char *SAXSfilename, double weight
 * Return Type : double
 *************************************************************************/
double getChiEnergy_inside(char *filename, char *SAXSfilename, double weight) {
    double energy = 0.0;
	char** Files;
	Files = new char*[1000];// initialize the double pointer
	Files[0]=new char[1000];// initialize 1st char*, with capacity of 4 chars
	Files[1]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[2]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[3]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[4]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[5]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[6]=new char[1000];// initialize 2nd char*, with capacity of 5 chars

	strcpy(Files[0],"foxs");//copy some data to 1st string
	strcpy(Files[1],"-q");//copy some data to 2nd string
	strcpy(Files[2],"0.3");//copy some data to 2nd string
	strcpy(Files[3],"-s");//copy some data to 2nd string
	strcpy(Files[4],"500");//copy some data to 2nd string
	strcpy(Files[5],filename);//copy some data to 2nd string
	strcpy(Files[6],SAXSfilename);//copy some data to 2nd string
	

	double chi_score=0;
	
	if(tofitting)
	{
		if(CAfitting)
		{
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[8]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"-n");
			strcpy(Files[8],"-r");//copy some data to 2nd string
			//for(int i = 0; i < 9;i++)
			//{
			//	cout << i << " : " << Files[i] << endl;
			//}
			 chi_score = run_foxs(9,Files);
			 //cout <<  " Fitting SAXS profile using CA atom only" << endl;

		}else{
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"-n");//copy some data to 2nd string
			//for(int i = 0; i < 8;i++)
			//{
			//	cout << i << " : " << Files[i] << endl;
			//}
			 chi_score = run_foxs(8,Files);
			 //cout <<  " Fitting SAXS profile" << endl;
		}
	}else{
		//for(int i = 0; i < 7;i++)
		//{
		//	cout << i << " : " << Files[i] << endl;
		//}
		if(CAfitting)
		{
			Files[7]=new char[1000];
			strcpy(Files[7],"-r");
			chi_score = run_foxs(8,Files);
			//cout <<  " Calculating default SAXS using CA atom only" << endl;
		}else{
			chi_score = run_foxs(7,Files);
			//cout <<  " Calculating default SAXS using CA atom only" << endl;
		}
		
	}
    energy =chi_score;
	delete[] Files;
	double energy_outside = getChiEnergy(filename,SAXSfilename,weight);
	//cout << "The final external chi score of "<< filename <<" is: " <<energy_outside<< endl;
	//cout << "The final inside chi score of "<< filename <<" is: " <<energy<< endl;
	//cout << "The final internal score: " <<chi_score<< endl;
    return energy * weight;

}




/*************************************************************************
 * Name        : getChiEnergy inside program
 * Purpose     : gets the chi-score energy
 * Arguments   : char *filename, char *SAXSfilename, double weight
 * Return Type : double
 *************************************************************************/
double getChiEnergy_inside2(char *filename, char *SAXSfilename, string pdbString, double weight) {
    double energy = 0.0;
	char** Files;
	Files = new char*[1000];// initialize the double pointer
	Files[0]=new char[1000];// initialize 1st char*, with capacity of 4 chars
	Files[1]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[2]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[3]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[4]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[5]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[6]=new char[1000];// initialize 2nd char*, with capacity of 5 chars

	strcpy(Files[0],"foxs");//copy some data to 1st string
	strcpy(Files[1],"-q");//copy some data to 2nd string
	strcpy(Files[2],"0.3");//copy some data to 2nd string
	strcpy(Files[3],"-s");//copy some data to 2nd string
	strcpy(Files[4],"500");//copy some data to 2nd string
	strcpy(Files[5],filename);//copy some data to 2nd string
	strcpy(Files[6],SAXSfilename);//copy some data to 2nd string
	

	double chi_score=0;
	
	if(tofitting)
	{
		if(CAfitting)
		{
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[8]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"-n");
			strcpy(Files[8],"-r");//copy some data to 2nd string
			//for(int i = 0; i < 9;i++)
			//{
			//	cout << i << " : " << Files[i] << endl;
			//}
			 //chi_score = run_foxs(9,Files);
			 chi_score = run_foxs_saxs(9,Files,pdbString);
			 //cout <<  " Fitting SAXS profile using CA atom only" << endl;

		}else{
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"-n");//copy some data to 2nd string
			//for(int i = 0; i < 8;i++)
			//{
			//	cout << i << " : " << Files[i] << endl;
			//}
			 //chi_score = run_foxs(8,Files);
			 chi_score = run_foxs_saxs(8,Files,pdbString);
			 //cout <<  " Fitting SAXS profile" << endl;
		}
	}else{
		//for(int i = 0; i < 7;i++)
		//{
		//	cout << i << " : " << Files[i] << endl;
		//}
		if(CAfitting)
		{
			Files[7]=new char[1000];
			strcpy(Files[7],"-r");
			//chi_score = run_foxs(8,Files);
			chi_score = run_foxs_saxs(8,Files,pdbString);
			//cout <<  " Calculating default SAXS using CA atom only" << endl;
		}else{
			//chi_score = run_foxs(7,Files);
			chi_score = run_foxs_saxs(7,Files,pdbString);
			//cout <<  " Calculating default SAXS using CA atom only" << endl;
		}
		
	}
    energy =chi_score;
	delete[] Files;
	double energy_outside = getChiEnergy(filename,SAXSfilename,weight);
	//cout << "The final external chi score of "<< filename <<" is: " <<energy_outside<< endl;
	//cout << "The final inside chi score of "<< filename <<" is: " <<energy<< endl;
	//cout << "The final internal score: " <<chi_score<< endl;
    return energy * weight;

}




//int run_foxs(int foxsargc, char** foxsargv) {
double run_foxs_only(int foxsargc, char** foxsargv) {
  // output arguments
  // comment by Jie
  clock_t foxs_start = clock();

  for (int i = 0; i < foxsargc; i++) std::cerr << foxsargv[i] << " ";
  std::cerr << std::endl;

  int profile_size = 500;
  float max_q = 0.0; // change after read
  float min_c1 = 0.99;
  float max_c1 = 1.05;
  float min_c2 = -0.5;
  float max_c2 = 2.0;
  bool heavy_atoms_only = true;
  bool residue_level = false;
  float background_adjustment_q = 0.0;
  bool use_offset = false;
  bool write_partial_profile = false;
  int multi_model_pdb = 1;
  bool vr_score = false;
  bool score_log = false;
  bool gnuplot_script = false;
  bool tofitting = false; 
  float pr_dmax = 0.0;  //  jie put this here from hidden option
  std::string jobname = "SAXS_fitting";  //  jie put this here from hidden option

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Any number of input PDBs and profiles is supported. \
Each PDB will be fitted against each profile.")
    ("version", "FoXS (IMP applications)\nCopyright 2007-2016 IMP Inventors.\n\
All rights reserved. \nLicense: GNU LGPL version 2.1 or later\n\
<http://gnu.org/licenses/lgpl.html>.\n\
Written by Dina Schneidman.")
    ("profile_size,s", po::value<int>(&profile_size)->default_value(500, "500"),
     "number of points in the profile")
    ("jobname,w", po::value<std::string>(&jobname)->default_value("SAXS_fitting", "SAXS_fitting"), "job id")
    ("max_q,q", po::value<float>(&max_q)->default_value(0.5, "0.50"), "max q value")
    ("min_c1", po::value<float>(&min_c1)->default_value(0.99, "0.99"), "min c1 value")
    ("max_c1", po::value<float>(&max_c1)->default_value(1.05, "1.05"), "max c1 value")
    ("min_c2", po::value<float>(&min_c2)->default_value(-2.0, "-2.00"), "min c2 value")
    ("max_c2", po::value<float>(&max_c2)->default_value(4.0, "4.00"), "max c2 value")
    ("tofitting,n", "only default chi is calculated (default = false)") 
    ("hydrogens,h", "explicitly consider hydrogens in PDB files (default = false)")
    ("residues,r", "fast coarse grained calculation using CA atoms only (default = false)")
    ("background_q,b", po::value<float>(&background_adjustment_q)->default_value(0.0),
     "background adjustment, not used by default. if enabled, recommended q value is 0.2")
    ("offset,o", "use offset in fitting (default = false)")
    ("write-partial-profile,p", "write partial profile file (default = false)")
    ("multi-model-pdb,m", po::value<int>(&multi_model_pdb)->default_value(1),
     "1 - read the first MODEL only (default), \
2 - read each MODEL into a separate structure, \
3 - read all models into a single structure")
    ("volatility_ratio,v","calculate volatility ratio score (default = false)")
    ("score_log,l", "use log(intensity) in fitting and scoring (default = false)")
    ("gnuplot_script,g", "print gnuplot script for gnuplot viewing (default = false)")
	("pr_dmax", po::value<float>(&pr_dmax)->default_value(0.0, "0.0"),
     "Dmax value for P(r) calculation. P(r) is calculated only is pr_dmax > 0");

  std::string form_factor_table_file;
  std::string beam_profile_file;
  bool ab_initio = false;
  bool vacuum = false;
  bool javascript = false;
  int chi_free = 0;
  /*float pr_dmax = 0.0;*/
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files")
    ("form_factor_table,f", po::value<std::string>(&form_factor_table_file),
     "ff table name")
    ("beam_profile", po::value<std::string>(&beam_profile_file),
     "beam profile file name for desmearing")
    ("ab_initio,a", "compute profile for a bead model with \
constant form factor (default = false)")
    ("vacuum", "compute profile in vacuum (default = false)")
    ("javascript,j",
     "output javascript for browser viewing of the results (default = false)")
    ("chi_free,x", po::value<int>(&chi_free)->default_value(0),
     "compute chi-free instead of chi, specify iteration number (default = 0)");
    /*("pr_dmax", po::value<float>(&pr_dmax)->default_value(0.0, "0.0"),
     "Dmax value for P(r) calculation. P(r) is calculated only is pr_dmax > 0");*/

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::options_description visible(
      "Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ... ");
  visible.add(desc);

  po::positional_options_description p;
  p.add("input-files", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(foxsargc, foxsargv)
                .options(cmdline_options)
                .positional(p)
                .run(),
            vm);
  po::notify(vm);

  bool fit = true;
  std::vector<std::string> files, pdb_files, dat_files;
  if (vm.count("input-files")) {
    files = vm["input-files"].as<std::vector<std::string> >();
  }
  if (vm.count("help") || files.size() == 0) {
    std::cout << visible << "\n";
    return 0;
  }
  if (vm.count("tofitting")) tofitting = true; 
  if (vm.count("hydrogens")) heavy_atoms_only = false;
  if (vm.count("residues")) residue_level = true;
  if (vm.count("offset")) use_offset = true;
  if (vm.count("write-partial-profile")) write_partial_profile = true;
  if (vm.count("score_log")) score_log = true;
  if (vm.count("gnuplot_script")) gnuplot_script = true;

  // no water layer or fitting in ab initio mode for now
  if (vm.count("ab_initio")) {
    ab_initio = true;
    fit = false;
  }
  if (vm.count("vacuum")) {
    vacuum = true;
  }
  if (vm.count("javascript")) {
    javascript = true;
  }
  if (vm.count("volatility_ratio")) {
    vr_score = true;
  }

  if (multi_model_pdb != 1 && multi_model_pdb != 2 && multi_model_pdb != 3) {
    std::cerr << "Incorrect option for multi_model_pdb " << multi_model_pdb
              << std::endl;
    std::cerr << "Use 1 to read first MODEL only\n"
              << "    2 to read each MODEL into a separate structure,\n"
              << "    3 to read all models into a single structure\n";
    std::cerr << "Default value of 1 is used\n";
    multi_model_pdb = 1;
  }

  //IMP::benchmark::Profiler pp("prof_out");

  // determine form factor type
  FormFactorType ff_type = HEAVY_ATOMS;
  if (!heavy_atoms_only) ff_type = ALL_ATOMS;
  if (residue_level) ff_type = CA_ATOMS;

  // 1. read pdbs and profiles, prepare particles
  std::vector<IMP::Particles> particles_vec;
  Profiles exp_profiles;
  //cout << "Check here" << endl;
  read_files(files, pdb_files, dat_files, particles_vec, exp_profiles,
             residue_level, heavy_atoms_only, multi_model_pdb, max_q);

  if (background_adjustment_q > 0.0) {
    for (unsigned int i = 0; i < exp_profiles.size(); i++)
      exp_profiles[i]->background_adjust(background_adjustment_q);
  }

  if (exp_profiles.size() == 0 && !write_partial_profile) fit = false;

  if (max_q == 0.0) { // determine max_q
    if (exp_profiles.size() > 0) {
      for (unsigned int i = 0; i < exp_profiles.size(); i++)
        if (exp_profiles[i]->get_max_q() > max_q)
          max_q = exp_profiles[i]->get_max_q();
    } else {
      max_q = 0.5;
    }
  }
  float delta_q = max_q / profile_size;

  // read in or use default form factor table
  bool reciprocal = false;
  FormFactorTable* ft = NULL;
  if (form_factor_table_file.length() > 0) {
    // reciprocal space calculation, requires form factor file
    ft = new FormFactorTable(form_factor_table_file, 0.0, max_q, delta_q);
    reciprocal = true;
  } else {
    ft = get_default_form_factor_table();
  }

   // 2. get pr from experimental profiles , added by Jie
    std::vector<std::vector<pr_dist> > exp_prInfo;
    for (unsigned int j = 0; j < dat_files.size(); j++) {
      Profile* exp_saxs_profile = exp_profiles[j];
      std::string pr_file_name = trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) +
          ".dat.pr";
		// calculate P(r)
		if(pr_dmax > 0.0) {
		  RadialDistributionFunction pr(0.5);
		  exp_saxs_profile->profile_2_distribution(pr, pr_dmax);
		  pr.normalize();
		  std::ofstream pr_file(pr_file_name.c_str());
		  pr.show(pr_file);
		  
		  std::vector<pr_dist> prInfo;
		  prInfo.clear();
		  pr_dist prData;
		  for (unsigned int i = 0; i < pr.size(); i++) {
			//std::cout << pr.get_distance_from_index(i) << " " << (pr)[i] << std::endl;
			
			prData.radius = pr.get_distance_from_index(i);
			prData.dist = (pr)[i];
			
			prInfo.push_back(prData);

		  }
		  exp_prInfo.push_back(prInfo);
		}
	}
  // 2. compute profiles for input pdbs
  Profiles profiles;
  std::vector<FitParameters> fps;
  std::vector<std::vector<pr_dist> > pdb_prInfo;  // add by jie
  std::vector<std::string > pdb_fittingInfo;  // add by jie
  std::vector<RG_saxs> RG_fittingInfo;  // add by jie
  std::vector<double> chi_fittingInfo;  // add by jie
  
  double chi_score = 0.0;
  for (unsigned int i = 0; i < particles_vec.size(); i++) {
    //std::cerr << "Computing profile for " << pdb_files[i] << " "
     //         << particles_vec[i].size() << " atoms " << std::endl; // comment by Jie
    IMP::Pointer<Profile> profile =
        compute_profile(particles_vec[i], 0.0, max_q, delta_q, ft, ff_type,
                        fit, fit, reciprocal, ab_initio, vacuum,
                        beam_profile_file);

    // save the profile
    profiles.push_back(profile);
    // write profile file
/*  comment by Jie*/
    std::string profile_file_name = std::string(pdb_files[i]) + ".dat";
    if (write_partial_profile)
      profile->write_partial_profiles(profile_file_name);
    else {  // write normal profile
      profile->add_errors();
      profile->write_SAXS_file(profile_file_name);
      if (gnuplot_script) Gnuplot::print_profile_script(pdb_files[i]);
    }

    // calculate P(r)
    if(pr_dmax > 0.0) {
      RadialDistributionFunction pr(0.5);
      profile->profile_2_distribution(pr, pr_dmax);
      pr.normalize();
      std::string pr_file_name = std::string(pdb_files[i]) + ".pr";
      std::ofstream pr_file(pr_file_name.c_str());
      pr.show(pr_file);  
	  
	  std::vector<pr_dist> prInfo;
	  prInfo.clear();
	  pr_dist prData;
	  for (unsigned int m = 0; m < pr.size(); m++) {
		//std::cout << pr.get_distance_from_index(i) << " " << (pr)[i] << std::endl;
		
		prData.radius = pr.get_distance_from_index(m);
		prData.dist = (pr)[m];
		
		prInfo.push_back(prData);

	  }
	  pdb_prInfo.push_back(prInfo);
    }
	
	
 
	
    // 3. fit experimental profiles
    for (unsigned int j = 0; j < dat_files.size(); j++) {
	  //cout << "datfile: " << j << " : " << dat_files[j] << endl;
      Profile* exp_saxs_profile = exp_profiles[j];
      std::string fit_file_name2 =
        trim_extension(pdb_files[i]) + "_" +
        trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) +
          ".dat";
         
		//string fit_profile_array;
		std::vector<std::string> fit_profile_array;
		//Compute default score directly
		if(!tofitting)
		{		
			IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
			//chi_score = pf->compute_score(profile, use_offset, fit_file_name2);
			chi_score = pf->compute_score_custom(profile,fit_profile_array, use_offset, fit_file_name2);  // jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
			cout << "The no fitting score is: " << chi_score << endl; 
			
			chi_fittingInfo.push_back(chi_score);
		}else{
		
		  FitParameters fp;
		  if (score_log) {
			IMP_NEW(ProfileFitter<ChiScoreLog>, pf, (exp_saxs_profile));
			//fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
			//					 use_offset, fit_file_name2);
			fp = pf->fit_profile_custom(profile,fit_profile_array, min_c1, max_c1, min_c2, max_c2,// jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
								 use_offset, fit_file_name2);
		  } else {
			if (vr_score) {
			  IMP_NEW(ProfileFitter<RatioVolatilityScore>, pf, (exp_saxs_profile));
			  //fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
				//				   use_offset, fit_file_name2);
			  fp = pf->fit_profile_custom(profile,fit_profile_array, min_c1, max_c1, min_c2, max_c2,// jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
								   use_offset, fit_file_name2);
			} else {
			  IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
			  //fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
				//				   use_offset, fit_file_name2);
			  fp = pf->fit_profile_custom(profile,fit_profile_array, min_c1, max_c1, min_c2, max_c2,// jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
								   use_offset, fit_file_name2);
			  if (chi_free > 0) {
				float dmax = compute_max_distance(particles_vec[i]);
				unsigned int ns = IMP::algebra::get_rounded(
							   exp_saxs_profile->get_max_q() * dmax / IMP::PI);
				int K = chi_free;
				IMP_NEW(ChiFreeScore, cfs, (ns, K));
				cfs->set_was_used(true);
				// IMP_NEW(RatioVolatilityScore, rvs, ());
				// rvs->set_was_used(true);
				// resample the profile
				IMP_NEW(Profile, resampled_profile,
						(exp_saxs_profile->get_min_q(), exp_saxs_profile->get_max_q(),
						 exp_saxs_profile->get_delta_q()));
				pf->resample(profile, resampled_profile);
				float chi_free =
				  cfs->compute_score(exp_saxs_profile, resampled_profile);
				fp.set_chi(chi_free);
			  }
			}
		  }
		  chi_score = fp.get_chi();
		  double default_chi_score = fp.get_default_chi();
		  //cout << "The score is: " << chi_score << endl;
		  //cout << "The default score is: " << default_chi_score << endl;
			std::cout << pdb_files[i] << " " << dat_files[j]
			   << " Chi = " << fp.get_chi() << " c1 = " << fp.get_c1()
				<< " c2 = " << fp.get_c2()
				<< " default chi = " << fp.get_default_chi() << std::endl;
			
			chi_fittingInfo.push_back(chi_score);
		}
		
		std::string fit_profile_array_str="";
		for(int h = 0; h< fit_profile_array.size();h++)
		{
			//cout << fit_profile_array[h].c_str() << endl;
			fit_profile_array_str += fit_profile_array[h];
			//char charInfo[1000];
			//std::strcpy(charInfo, fit_profile_array[h].c_str());
			//fit_profile_array_str.append(charInfo);  
		}
		
		
	
		pdb_fittingInfo.push_back(fit_profile_array_str);
		
		
		
		
		// calculate RG from dat, add by Jie
		double end_q_rg = 1; // set to default 1
		double expRG = exp_saxs_profile->radius_of_gyration(end_q_rg);
		
		
		// calcualte RG from pdb , added by jie
		IMP::Particles pdb_particles = particles_vec[i];
		//RadiusOfGyrationRestraint pdbRG(pdb_particles,exp_saxs_profile,1.3);   // set q_end_rg as default 1.3
		IMP::algebra::Vector3D centroid(0.0, 0.0, 0.0);
		IMP::Vector<IMP::algebra::Vector3D> coordinates(pdb_particles.size());
		get_coordinates(pdb_particles, coordinates);
		for (unsigned int m = 0; m < pdb_particles.size(); m++) {
			centroid += coordinates[m];
		  }
		  centroid /= pdb_particles.size();
		  double radg = 0;
		  for (unsigned int m = 0; m < pdb_particles.size(); m++) {
			radg += get_squared_distance(coordinates[m], centroid);
		  }
		  radg /= pdb_particles.size();
		  radg = sqrt(radg);

		  double RGdiff = abs(radg - expRG);  // TODO: improve  
		  double RGdiff_normalize = abs(radg - expRG) / expRG;  // TODO: improve  

		  
		  
		//cout << "RG for " <<  trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) << ".dat: " <<  expRG << endl;
		//cout << "RG for " <<  trim_extension(pdb_files[i]) << ".pdb: " <<  radg << endl;
		//cout << "RG diff for " <<  fit_file_name2  << ": " <<  RGdiff << endl;
		//cout << "RG diff normalize for " <<  fit_file_name2 << ": " <<  RGdiff_normalize << endl;
		
		
		  std::vector<RG_saxs> rgInfo;
		  rgInfo.clear();
		  RG_saxs rgData;
		  rgData.expRG = expRG;
		  rgData.modelRG = radg;
		  rgData.RGdiff = RGdiff;
		  rgData.RGdiff_normalize = RGdiff_normalize;
		  RG_fittingInfo.push_back(rgData);
	  
	  
		  
/*
      std::cout << pdb_files[i] << " " << dat_files[j]
               << " Chi = " << fp.get_chi() << " c1 = " << fp.get_c1()
                << " c2 = " << fp.get_c2()
                << " default chi = " << fp.get_default_chi() << std::endl;
      fp.set_pdb_file_name(pdb_files[i]);
      fp.set_profile_file_name(dat_files[j]);
      fp.set_mol_index(i);
      if (gnuplot_script) Gnuplot::print_fit_script(fp);
      fps.push_back(fp);
*/
    }
  }

	//cout << "Writing score to file! " << endl;	  
	// calculate profile functions, KL divergence, added by Jie 
	//cout << "pdb_prInfo.size(): " << pdb_prInfo.size() << endl;
	//cout << "exp_prInfo.size(): " << exp_prInfo.size() << endl;
	
	if(pdb_prInfo.size()!=0)
	{
		
		
		if(odir)
		{
			sprintf(statScoreFile, "%s/%s_saxsfitting_score.txt", outputdir,const_cast<char*>(jobname.c_str()));

		}else{
			sprintf(statScoreFile, "%s_saxsfitting_score.txt",const_cast<char*>(jobname.c_str()));

		}	

		
		//cout << "Writing score to file: " << statScoreFile << endl;	  
		ofstream statsaxs(statScoreFile);
		if (!statsaxs.is_open()) {
			cout << "Error! folding statistics file can not open " << statScoreFile << endl;
			exit(0);
		}
		
		// buffer statistics
		char bufStat[1000];
		sprintf(bufStat,"Exp_dat\tDecoy_name\tChi_score\tRG_exp\tRG_model\tRG_diff\tRG_diff_normalize\tKL_div\tKL_div_sym\tScore1\tScore2\tScore3\tScore4\tScore5\tScore6\tScore7\tScore8\tScore9\tScore10\tScore11");
		statsaxs << bufStat << endl;
	
		for (int j = 0; j < pdb_prInfo.size(); j++) {
			for (int m = 0; m < exp_prInfo.size(); m++) {
				string fit_profile_array = pdb_fittingInfo[j+m];
				RG_saxs rgData = RG_fittingInfo[j+m];
				double chi_score = chi_fittingInfo[j+m];
				
				// calculate profile function used in SAXSTER http://www.cell.com/cms/attachment/2073830349/2068835456/mmc1.pdf
				std::string line;
				std::istringstream split(fit_profile_array);
				std::vector<fitting_profile> fittingdatInfo;
				fittingdatInfo.clear();
				while(std::getline(split, line, ';')) {
					//std::cout  << "line:" << line << endl;
					
					std::string line2;
					std::istringstream split2(line);
					fitting_profile datData;
					std::vector<double> vect;
					while(std::getline(split2, line2, ' ')) {
						//vect.push_back(std::strtod(line2.c_str() ));
						//cout << "--> "<< line2<<endl;
						vect.push_back(atof(line2.c_str()));
						
					}
					if(vect.size() !=4)
					{
						std::cout<<"The fitting data has incorrect column number!" <<std::endl;
						exit(-1); 
					}
					datData.q = vect[0];
					datData.expInt = vect[1];
					datData.modelInt = vect[2];
					datData.error = vect[3];
					fittingdatInfo.push_back(datData);
				}
				
				double fun_score1=0;
				double fun_score2=0;
				double fun_score2_nom=0;
				double fun_score2_denom=0;
				double fun_score3=0;
				double fun_score3_nom=0;
				double fun_score3_denom=0;
				double fun_score4=0;
				double fun_score5=0;
				double fun_score6=0;
				double fun_score6_nom=0;
				double fun_score6_denom=0;
				
				vector<double> fun_score8_expInt;
				vector<double> fun_score8_modelInt;
				for (unsigned int i = 0; i < fittingdatInfo.size(); i++) { 
					  //std::cout.precision(10);
					//std::cout << fittingdatInfo[i].q << " " << fittingdatInfo[i].expInt << " " << fittingdatInfo[i].modelInt<< " " << fittingdatInfo[i].error << std::endl;
					
					if(fittingdatInfo[i].expInt < 10e-10)
					{
						fittingdatInfo[i].expInt = 10e-10;
					}
					if(fittingdatInfo[i].modelInt < 10e-10)
					{
						fittingdatInfo[i].modelInt = 10e-10;
					}
						
					
					//function 1
					fun_score1 +=((fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt)*(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt))/(fittingdatInfo[i].error*fittingdatInfo[i].error);
					//function 2
					fun_score2_nom +=std::abs(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt);
					fun_score2_denom +=std::abs(fittingdatInfo[i].expInt);
					//function 3
					fun_score3_nom +=std::abs(log2(fittingdatInfo[i].expInt) - log2(fittingdatInfo[i].modelInt));
					fun_score3_denom +=std::abs(log2(fittingdatInfo[i].expInt));
					//function4
					fun_score4 += fittingdatInfo[i].q * ((fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt)*(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt));
					//function5
					fun_score5 += fittingdatInfo[i].q * ((log2(fittingdatInfo[i].expInt) - log2(fittingdatInfo[i].modelInt))*(log2(fittingdatInfo[i].expInt) - log2(fittingdatInfo[i].modelInt)));
					
					//function 6
					fun_score6_nom +=fittingdatInfo[i].q *fittingdatInfo[i].q *std::abs(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt);
					fun_score6_denom +=fittingdatInfo[i].q *fittingdatInfo[i].q*std::abs(fittingdatInfo[i].expInt);
					
					//function 8
					fun_score8_expInt.push_back(fittingdatInfo[i].expInt);
					fun_score8_modelInt.push_back(fittingdatInfo[i].modelInt);
				}
				
				fun_score1  /= fittingdatInfo.size();
				fun_score2 = fun_score2_nom /fun_score2_denom;
				fun_score3 = fun_score3_nom  / fun_score3_denom;
				fun_score6 = fun_score6_nom  / fun_score6_denom;
				
				double fun_score8 = log2(1-pearson(fun_score8_expInt,fun_score8_modelInt));
				double fun_score10 = log2(1-cosine(fun_score8_expInt,fun_score8_modelInt));
				double IntPCC = pearson(fun_score8_expInt,fun_score8_modelInt);
				double Intcosine = cosine(fun_score8_expInt,fun_score8_modelInt);
				double Inteuclidean = euclidean(fun_score8_expInt,fun_score8_modelInt);
				double Intjaccard = jaccard(fun_score8_expInt,fun_score8_modelInt);
				
			
				
				std::string exp_file_name = trim_extension(basename(const_cast<char*>(dat_files[m].c_str()))) +
			  ".dat";
			  std::string pdb_file_name = std::string(pdb_files[j]);
			  std::vector<pr_dist> e_prInfo=exp_prInfo[m];
			  std::vector<pr_dist> p_prInfo=pdb_prInfo[j];
			  int min_size=0;
			  if(e_prInfo.size() < p_prInfo.size())
			  {
				  min_size = e_prInfo.size();
			  }else{
				  min_size = p_prInfo.size();
			  }
			  double KL=0;
			  double KL_symetric=0;
			  double fun_score7=0;
			  vector<double> fun_score9_expPr;
			  vector<double> fun_score9_modelPr;
			  //std::cout << "Particle size for KL: "<< min_size << std::endl << std::endl;
			  for (int k = 0; k < min_size; k++) {
				if(e_prInfo[k].radius != p_prInfo[k].radius)
				{ 
					std::cout << "Error! the pr format is not correct!" << std::endl << std::endl;
					exit(-1);
				}
				if(e_prInfo[k].dist !=0 and p_prInfo[k].dist !=0)
				{
					if(e_prInfo[k].dist < 10e-10)
					{
						e_prInfo[k].dist = 10e-10;
					}
					if(p_prInfo[k].dist < 10e-10)
					{
						p_prInfo[k].dist = 10e-10;
					}
					KL += p_prInfo[k].dist *log2(p_prInfo[k].dist /e_prInfo[k].dist );
					KL_symetric += p_prInfo[k].dist *log2(p_prInfo[k].dist /e_prInfo[k].dist ) + e_prInfo[k].dist *log2(e_prInfo[k].dist /p_prInfo[k].dist );
					fun_score7 += (p_prInfo[k].dist - e_prInfo[k].dist)*(p_prInfo[k].dist - e_prInfo[k].dist);
					
					fun_score9_expPr.push_back(p_prInfo[k].dist);
					fun_score9_modelPr.push_back(e_prInfo[k].dist);
					
				}
			  }
			  double fun_score9 = log2(1-pearson(fun_score9_expPr,fun_score9_modelPr));
			  double fun_score11 = log2(1-cosine(fun_score9_expPr,fun_score9_modelPr));
			  double PrPCC = pearson(fun_score9_expPr,fun_score9_modelPr);
			  double Prcosine = cosine(fun_score9_expPr,fun_score9_modelPr);
			  double Preuclidean = euclidean(fun_score9_expPr,fun_score9_modelPr);
			  double Prjaccard = jaccard(fun_score9_expPr,fun_score9_modelPr);
				
			  if (KL  <0)
			  {
				  KL = 0.0;
			  }
			  if (KL_symetric  <0)
			  {
				  KL_symetric = 0.0;
			  }
			  
			  //std::cout << exp_file_name << "&" << pdb_file_name<< " KL: " << KL << " score1: " << fun_score1 << " score2: " << fun_score2 << " score3: " << fun_score3;
			  //std::cout << " score3: " << fun_score3 << " score4: " << fun_score4 << " score5: " << fun_score5 << " score6: " << fun_score6;
			  //std::cout << " score7: " << fun_score7 << " score8: " << fun_score8 << " score9: " << fun_score9 << std::endl;
			  //std::cout << " PrPCC: " << PrPCC << " Prcosine: " << Prcosine << " Preuclidean: " << Preuclidean << " Prjaccard: " << Prjaccard << std::endl;
			  //std::cout << " IntPCC: " << IntPCC << " Intcosine: " << Intcosine << " Inteuclidean: " << Inteuclidean << " Intjaccard: " << Intjaccard << std::endl;
			  
			  
			  char bufStat[1000];
			  sprintf(bufStat,"%s\t%s\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f",basename(const_cast<char*>(exp_file_name.c_str())),basename(const_cast<char*>(pdb_file_name.c_str())),chi_score,rgData.expRG,rgData.modelRG,rgData.RGdiff,rgData.RGdiff_normalize,KL,KL_symetric,fun_score1,fun_score2,fun_score3,fun_score4,fun_score5,fun_score6,fun_score7,fun_score8,fun_score9,fun_score10,fun_score11 );
			  statsaxs << bufStat << endl;
		
		  }
		}
	}
  //std::sort(fps.begin(), fps.end(), FitParameters::compare_fit_parameters());
 /*
  if (pdb_files.size() > 1 && gnuplot_script) {
    Gnuplot::print_profile_script(pdb_files);
    if (dat_files.size() > 0) Gnuplot::print_fit_script(fps);
  }

  if (javascript) {
    if (dat_files.size() > 0) {
      Gnuplot::print_canvas_script(fps, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(fps, particles_vec, "jmoltable");
    } else {
      Gnuplot::print_canvas_script(pdb_files, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(pdb_files, particles_vec, "jmoltable");
    }
  }
 */
 
	clock_t foxs_end = clock();
		double elapse_secs1 = double(foxs_end - foxs_start) / CLOCKS_PER_SEC;
		cout << "Foxs finished within "<< elapse_secs1 << " sec!" << endl << endl;
  return 0;
  //return chi_score;
}




/*************************************************************************
 * Name        : getChiEnergy inside program
 * Purpose     : gets the chi-score energy
 * Arguments   : char *filename, char *SAXSfilename, double weight
 * Return Type : double
 *************************************************************************/
saxsInfo getChiEnergy_inside_pro(char *filename, char *SAXSfilename, string pdbString,double weight) {
    double energy = 0.0;
	char** Files;
	Files = new char*[1000];// initialize the double pointer
	Files[0]=new char[1000];// initialize 1st char*, with capacity of 4 chars
	Files[1]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[2]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[3]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[4]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[5]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[6]=new char[1000];// initialize 2nd char*, with capacity of 5 chars

	strcpy(Files[0],"foxs");//copy some data to 1st string
	strcpy(Files[1],"-q");//copy some data to 2nd string
	strcpy(Files[2],"0.3");//copy some data to 2nd string
	strcpy(Files[3],"-s");//copy some data to 2nd string
	strcpy(Files[4],"500");//copy some data to 2nd string
	strcpy(Files[5],filename);//copy some data to 2nd string
	strcpy(Files[6],SAXSfilename);//copy some data to 2nd string
	

	//double saxs_score=0;
	saxsInfo saxsScore;
	if(tofitting)
	{
		if(CAfitting)
		{
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[8]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[9]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[10]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"-n");
			strcpy(Files[8],"-r");//copy some data to 2nd string
			strcpy(Files[9],"--pr_dmax");//copy some data to 2nd string
			strcpy(Files[10],"100");//copy some data to 2nd string
			//for(int i = 0; i < 9;i++)
			//{
			//	cout << i << " : " << Files[i] << endl;
			//}
			 //chi_score = run_foxs(9,Files);
			 //chi_score = run_foxs_saxs(9,Files,pdbString);
			 saxsScore = run_foxs_saxs_pro(11,Files,pdbString);
			 //cout <<  " Fitting SAXS profile using CA atom only" << endl;

		}else{
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"-n");//copy some data to 2nd string
			Files[8]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[9]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[8],"--pr_dmax");//copy some data to 2nd string
			strcpy(Files[9],"100");//copy some data to 2nd string
			//for(int i = 0; i < 8;i++)
			//{
			//	cout << i << " : " << Files[i] << endl;
			//}
			 //chi_score = run_foxs(8,Files);
			 //chi_score = run_foxs_saxs(8,Files,pdbString);
			 saxsScore = run_foxs_saxs_pro(10,Files,pdbString);
			 //cout <<  " Fitting SAXS profile" << endl;
		}
	}else{
		//for(int i = 0; i < 7;i++)
		//{
		//	cout << i << " : " << Files[i] << endl;
		//}
		if(CAfitting)
		{
			Files[7]=new char[1000];
			Files[8]=new char[1000];
			Files[9]=new char[1000];
			strcpy(Files[7],"-r");
			strcpy(Files[8],"--pr_dmax");//copy some data to 2nd string
			strcpy(Files[9],"100");//copy some data to 2nd string
			//chi_score = run_foxs(8,Files);
			//chi_score = run_foxs_saxs(8,Files,pdbString);
			saxsScore = run_foxs_saxs_pro(10,Files,pdbString);
			//cout <<  " Calculating default SAXS using CA atom only" << endl;
		}else{
			//chi_score = run_foxs(7,Files);
			Files[7]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			Files[8]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
			strcpy(Files[7],"--pr_dmax");//copy some data to 2nd string
			strcpy(Files[8],"100");//copy some data to 2nd string
			//chi_score = run_foxs_saxs(7,Files,pdbString);
			saxsScore = run_foxs_saxs_pro(9,Files,pdbString);
			//cout <<  " Calculating default SAXS using CA atom only" << endl;
		}
		
	}
	
	delete[] Files;
    //energy =saxs_score;
	
	//double energy_outside = getChiEnergy(filename,SAXSfilename,weight);
	//cout << "The final external chi score of "<< filename <<" is: " <<energy_outside<< endl;
	//cout << "The final inside chi score of "<< filename <<" is: " <<energy<< endl;
	//cout << "The final internal score: " <<chi_score<< endl;
    //return energy * weight;
    return saxsScore;

}



//int run_foxs(int foxsargc, char** foxsargv) {
saxsInfo run_foxs_saxs_pro(int foxsargc, char** foxsargv, string inputString) {
  // output arguments
  // comment by Jie
  //clock_t foxs_start = clock();
  //for (int i = 0; i < foxsargc; i++) std::cerr << foxsargv[i] << " ";
 // std::cerr << std::endl;
  
  int profile_size = 500;
  float max_q = 0.0; // change after read
  float min_c1 = 0.99;
  float max_c1 = 1.05;
  float min_c2 = -0.5;
  float max_c2 = 2.0;
  bool heavy_atoms_only = true;
  bool residue_level = false;
  float background_adjustment_q = 0.0;
  bool use_offset = false;
  bool write_partial_profile = false;
  int multi_model_pdb = 1;
  bool vr_score = false;
  bool score_log = false;
  bool gnuplot_script = false;
  bool tofitting = false; 
  float pr_dmax = 0.0;  //  jie put this here from hidden option
  std::string jobname = "SAXS_fitting";  //  jie put this here from hidden option

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Any number of input PDBs and profiles is supported. \
Each PDB will be fitted against each profile.")
    ("version", "FoXS (IMP applications)\nCopyright 2007-2016 IMP Inventors.\n\
All rights reserved. \nLicense: GNU LGPL version 2.1 or later\n\
<http://gnu.org/licenses/lgpl.html>.\n\
Written by Dina Schneidman.")
    ("profile_size,s", po::value<int>(&profile_size)->default_value(500, "500"),
     "number of points in the profile")
    ("jobname,w", po::value<std::string>(&jobname)->default_value("SAXS_fitting", "SAXS_fitting"), "job id")
    ("max_q,q", po::value<float>(&max_q)->default_value(0.5, "0.50"), "max q value")
    ("min_c1", po::value<float>(&min_c1)->default_value(0.99, "0.99"), "min c1 value")
    ("max_c1", po::value<float>(&max_c1)->default_value(1.05, "1.05"), "max c1 value")
    ("min_c2", po::value<float>(&min_c2)->default_value(-2.0, "-2.00"), "min c2 value")
    ("max_c2", po::value<float>(&max_c2)->default_value(4.0, "4.00"), "max c2 value")
    ("tofitting,n", "only default chi is calculated (default = false)") 
    ("hydrogens,h", "explicitly consider hydrogens in PDB files (default = false)")
    ("residues,r", "fast coarse grained calculation using CA atoms only (default = false)")
    ("background_q,b", po::value<float>(&background_adjustment_q)->default_value(0.0),
     "background adjustment, not used by default. if enabled, recommended q value is 0.2")
    ("offset,o", "use offset in fitting (default = false)")
    ("write-partial-profile,p", "write partial profile file (default = false)")
    ("multi-model-pdb,m", po::value<int>(&multi_model_pdb)->default_value(1),
     "1 - read the first MODEL only (default), \
2 - read each MODEL into a separate structure, \
3 - read all models into a single structure")
    ("volatility_ratio,v","calculate volatility ratio score (default = false)")
    ("score_log,l", "use log(intensity) in fitting and scoring (default = false)")
    ("gnuplot_script,g", "print gnuplot script for gnuplot viewing (default = false)")
	("pr_dmax", po::value<float>(&pr_dmax)->default_value(0.0, "0.0"),
     "Dmax value for P(r) calculation. P(r) is calculated only is pr_dmax > 0");

  std::string form_factor_table_file;
  std::string beam_profile_file;
  bool ab_initio = false;
  bool vacuum = false;
  bool javascript = false;
  int chi_free = 0;
  /*float pr_dmax = 0.0;*/
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files")
    ("form_factor_table,f", po::value<std::string>(&form_factor_table_file),
     "ff table name")
    ("beam_profile", po::value<std::string>(&beam_profile_file),
     "beam profile file name for desmearing")
    ("ab_initio,a", "compute profile for a bead model with \
constant form factor (default = false)")
    ("vacuum", "compute profile in vacuum (default = false)")
    ("javascript,j",
     "output javascript for browser viewing of the results (default = false)")
    ("chi_free,x", po::value<int>(&chi_free)->default_value(0),
     "compute chi-free instead of chi, specify iteration number (default = 0)");
    /*("pr_dmax", po::value<float>(&pr_dmax)->default_value(0.0, "0.0"),
     "Dmax value for P(r) calculation. P(r) is calculated only is pr_dmax > 0");*/

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::options_description visible(
      "Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ... ");
  visible.add(desc);

  po::positional_options_description p;
  p.add("input-files", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(foxsargc, foxsargv)
                .options(cmdline_options)
                .positional(p)
                .run(),
            vm);
  po::notify(vm);

  bool fit = true;
  std::vector<std::string> files, pdb_files, dat_files;
  if (vm.count("input-files")) {
    files = vm["input-files"].as<std::vector<std::string> >();
  }
  if (vm.count("help") || files.size() == 0) {
    std::cout << visible << "\n";
	exit(0);
    //return 0;
  }
  if (vm.count("tofitting")) tofitting = true; 
  if (vm.count("hydrogens")) heavy_atoms_only = false;
  if (vm.count("residues")) residue_level = true;
  if (vm.count("offset")) use_offset = true;
  if (vm.count("write-partial-profile")) write_partial_profile = true;
  if (vm.count("score_log")) score_log = true;
  if (vm.count("gnuplot_script")) gnuplot_script = true;

  // no water layer or fitting in ab initio mode for now
  if (vm.count("ab_initio")) {
    ab_initio = true;
    fit = false;
  }
  if (vm.count("vacuum")) {
    vacuum = true;
  }
  if (vm.count("javascript")) {
    javascript = true;
  }
  if (vm.count("volatility_ratio")) {
    vr_score = true;
  }

  if (multi_model_pdb != 1 && multi_model_pdb != 2 && multi_model_pdb != 3) {
    std::cerr << "Incorrect option for multi_model_pdb " << multi_model_pdb
              << std::endl;
    std::cerr << "Use 1 to read first MODEL only\n"
              << "    2 to read each MODEL into a separate structure,\n"
              << "    3 to read all models into a single structure\n";
    std::cerr << "Default value of 1 is used\n";
    multi_model_pdb = 1;
  }

  //IMP::benchmark::Profiler pp("prof_out");

  // determine form factor type
  FormFactorType ff_type = HEAVY_ATOMS;
  if (!heavy_atoms_only) ff_type = ALL_ATOMS;
  if (residue_level) ff_type = CA_ATOMS;

  // 1. read pdbs and profiles, prepare particles
  std::vector<IMP::Particles> particles_vec;

  Profiles exp_profiles_custom;  // set this as global variable, so only need be loaded in the first time
  //Profiles exp_profiles;  // set this as global variable, so only need be loaded in the first time
  
  
  
  //cout << "Check here" << endl;
  //read_files(files, pdb_files, dat_files, particles_vec, exp_profiles,
   //         residue_level, heavy_atoms_only, multi_model_pdb, max_q);
   
   
   read_files_saxs(files, pdb_files, dat_files,inputString, particles_vec, exp_profiles_custom,
             residue_level, heavy_atoms_only, multi_model_pdb, max_q);
   /*
  if(exp_profiles_loaded)
  {
	 read_files_saxs_pdbOnly(files, pdb_files, dat_files,inputString, particles_vec,
             residue_level, heavy_atoms_only, multi_model_pdb, max_q);
	  
			cout << "already loaded profile" << exp_profiles_custom.size()<< endl;
  }else{
	  read_files_saxs(files, pdb_files, dat_files,inputString, particles_vec, exp_profiles_custom,
             residue_level, heavy_atoms_only, multi_model_pdb, max_q);
			 
			cout << "start  loaded profile" << exp_profiles_custom.size() << endl;
		exp_profiles_loaded=true;
	  
  }
  */
  
  
  if (background_adjustment_q > 0.0) {
    for (unsigned int i = 0; i < exp_profiles_custom.size(); i++)
      exp_profiles_custom[i]->background_adjust(background_adjustment_q);
  }

  if (exp_profiles_custom.size() == 0 && !write_partial_profile) fit = false;

  if (max_q == 0.0) { // determine max_q
    if (exp_profiles_custom.size() > 0) {
      for (unsigned int i = 0; i < exp_profiles_custom.size(); i++)
        if (exp_profiles_custom[i]->get_max_q() > max_q)
          max_q = exp_profiles_custom[i]->get_max_q();
    } else {
      max_q = 0.5;
    }
  }
  float delta_q = max_q / profile_size;

  // read in or use default form factor table
  bool reciprocal = false;
  FormFactorTable* ft = NULL;
  if (form_factor_table_file.length() > 0) {
    // reciprocal space calculation, requires form factor file
    ft = new FormFactorTable(form_factor_table_file, 0.0, max_q, delta_q);
    reciprocal = true;
  } else {
    ft = get_default_form_factor_table();
  }

   // 2. get pr from experimental profiles , added by Jie
    std::vector<std::vector<pr_dist> > exp_prInfo;
    for (unsigned int j = 0; j < dat_files.size(); j++) {
      Profile* exp_saxs_profile = exp_profiles_custom[j];
      std::string pr_file_name = trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) +
          ".dat.pr";
		// calculate P(r)
		if(pr_dmax > 0.0) {
		  RadialDistributionFunction pr(0.5);
		  exp_saxs_profile->profile_2_distribution(pr, pr_dmax);
		  pr.normalize();
		  std::ofstream pr_file(pr_file_name.c_str());
		  pr.show(pr_file);
		  
		  std::vector<pr_dist> prInfo;
		  prInfo.clear();
		  pr_dist prData;
		  for (unsigned int i = 0; i < pr.size(); i++) {
			//std::cout << pr.get_distance_from_index(i) << " " << (pr)[i] << std::endl;
			
			prData.radius = pr.get_distance_from_index(i);
			prData.dist = (pr)[i];
			
			prInfo.push_back(prData);

		  }
		  exp_prInfo.push_back(prInfo);
		}
	}
	
  // 2. compute profiles for input pdbs
  Profiles profiles;
  std::vector<FitParameters> fps;
  std::vector<std::vector<pr_dist> > pdb_prInfo;  // add by jie
  std::vector<std::string > pdb_fittingInfo;  // add by jie
  std::vector<RG_saxs> RG_fittingInfo;  // add by jie
  std::vector<double> chi_fittingInfo;  // add by jie
  
  double chi_score = 0.0;
  for (unsigned int i = 0; i < particles_vec.size(); i++) {
    //std::cerr << "Computing profile for " << pdb_files[i] << " "
     //         << particles_vec[i].size() << " atoms " << std::endl; // comment by Jie
    IMP::Pointer<Profile> profile =
        compute_profile(particles_vec[i], 0.0, max_q, delta_q, ft, ff_type,
                        fit, fit, reciprocal, ab_initio, vacuum,
                        beam_profile_file);

    // save the profile
    profiles.push_back(profile);
    // write profile file
/*  comment by Jie
    std::string profile_file_name = std::string(pdb_files[i]) + ".dat";
    if (write_partial_profile)
      profile->write_partial_profiles(profile_file_name);
    else {  // write normal profile
      profile->add_errors();
      profile->write_SAXS_file(profile_file_name);
      if (gnuplot_script) Gnuplot::print_profile_script(pdb_files[i]);
    }
*/
    // calculate P(r)
    if(pr_dmax > 0.0) {
      RadialDistributionFunction pr(0.5);
      profile->profile_2_distribution(pr, pr_dmax);
      pr.normalize();
      //std::string pr_file_name = std::string(pdb_files[i]) + ".pr";
      //std::ofstream pr_file(pr_file_name.c_str());
      //pr.show(pr_file);  
	  
	  std::vector<pr_dist> prInfo;
	  prInfo.clear();
	  pr_dist prData;
	  for (unsigned int m = 0; m < pr.size(); m++) {
		//std::cout << pr.get_distance_from_index(i) << " " << (pr)[i] << std::endl;
		
		prData.radius = pr.get_distance_from_index(m);
		prData.dist = (pr)[m];
		
		prInfo.push_back(prData);

	  }
	  pdb_prInfo.push_back(prInfo);
    }
	

	
    // 3. fit experimental profiles
    for (unsigned int j = 0; j < dat_files.size(); j++) {
	  //cout << "datfile: " << j << " : " << dat_files[j] << endl;
      Profile* exp_saxs_profile = exp_profiles_custom[j];
      std::string fit_file_name2 =
        trim_extension(pdb_files[i]) + "_" +
        trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) +
          ".dat";
         
		//string fit_profile_array;
		std::vector<std::string> fit_profile_array;
		//Compute default score directly
		if(!tofitting)
		{		
			IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
			//chi_score = pf->compute_score(profile, use_offset, fit_file_name2);
			chi_score = pf->compute_score_custom(profile,fit_profile_array, use_offset, fit_file_name2);  // jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
			cout << "The no fitting score is: " << chi_score << endl; 
			
			chi_fittingInfo.push_back(chi_score);
		}else{
		
		  FitParameters fp;
		  if (score_log) {
			IMP_NEW(ProfileFitter<ChiScoreLog>, pf, (exp_saxs_profile));
			//fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
			//					 use_offset, fit_file_name2);
			fp = pf->fit_profile_custom(profile,fit_profile_array, min_c1, max_c1, min_c2, max_c2,// jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
								 use_offset, fit_file_name2);
		  } else {
			if (vr_score) {
			  IMP_NEW(ProfileFitter<RatioVolatilityScore>, pf, (exp_saxs_profile));
			  //fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
				//				   use_offset, fit_file_name2);
			  fp = pf->fit_profile_custom(profile,fit_profile_array, min_c1, max_c1, min_c2, max_c2,// jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
								   use_offset, fit_file_name2);
			} else {
			  IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
			  //fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
				//				   use_offset, fit_file_name2);
			  fp = pf->fit_profile_custom(profile,fit_profile_array, min_c1, max_c1, min_c2, max_c2,// jie add new function on 06/16/2017, into saxs/src/include/ProfileFitter.h
								   use_offset, fit_file_name2);
			  if (chi_free > 0) {
				float dmax = compute_max_distance(particles_vec[i]);
				unsigned int ns = IMP::algebra::get_rounded(
							   exp_saxs_profile->get_max_q() * dmax / IMP::PI);
				int K = chi_free;
				IMP_NEW(ChiFreeScore, cfs, (ns, K));
				cfs->set_was_used(true);
				// IMP_NEW(RatioVolatilityScore, rvs, ());
				// rvs->set_was_used(true);
				// resample the profile
				IMP_NEW(Profile, resampled_profile,
						(exp_saxs_profile->get_min_q(), exp_saxs_profile->get_max_q(),
						 exp_saxs_profile->get_delta_q()));
				pf->resample(profile, resampled_profile);
				float chi_free =
				  cfs->compute_score(exp_saxs_profile, resampled_profile);
				fp.set_chi(chi_free);
			  }
			}
		  }
		  chi_score = fp.get_chi();
		  double default_chi_score = fp.get_default_chi();
		  //cout << "The score is: " << chi_score << endl;
		  //cout << "The default score is: " << default_chi_score << endl;
		  /*
			std::cout << pdb_files[i] << " " << dat_files[j]
			   << " Chi = " << fp.get_chi() << " c1 = " << fp.get_c1()
				<< " c2 = " << fp.get_c2()
				<< " default chi = " << fp.get_default_chi() << std::endl;
			*/
			chi_fittingInfo.push_back(chi_score);
		
		}
		
		
		std::string fit_profile_array_str="";
		for(int h = 0; h< fit_profile_array.size();h++)
		{
			//cout << fit_profile_array[h].c_str() << endl;
			fit_profile_array_str += fit_profile_array[h];
			//char charInfo[1000];
			//std::strcpy(charInfo, fit_profile_array[h].c_str());
			//fit_profile_array_str.append(charInfo);  
		}
		pdb_fittingInfo.push_back(fit_profile_array_str);
		
		
		// calculate RG from dat, add by Jie
		double end_q_rg = 1; // set to default 1
		double expRG = exp_saxs_profile->radius_of_gyration(end_q_rg);
		
		
		// calcualte RG from pdb , added by jie
		IMP::Particles pdb_particles = particles_vec[i];
		//RadiusOfGyrationRestraint pdbRG(pdb_particles,exp_saxs_profile,1.3);   // set q_end_rg as default 1.3
		IMP::algebra::Vector3D centroid(0.0, 0.0, 0.0);
		IMP::Vector<IMP::algebra::Vector3D> coordinates(pdb_particles.size());
		get_coordinates(pdb_particles, coordinates);
		for (unsigned int m = 0; m < pdb_particles.size(); m++) {
			centroid += coordinates[m];
		  }
		  centroid /= pdb_particles.size();
		  double radg = 0;
		  for (unsigned int m = 0; m < pdb_particles.size(); m++) {
			radg += get_squared_distance(coordinates[m], centroid);
		  }
		  radg /= pdb_particles.size();
		  radg = sqrt(radg);

		  double RGdiff = abs(radg - expRG);  // TODO: improve  
		  double RGdiff_normalize = abs(radg - expRG) / expRG;  // TODO: improve  

		  
		  
		//cout << "RG for " <<  trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) << ".dat: " <<  expRG << endl;
		//cout << "RG for " <<  trim_extension(pdb_files[i]) << ".pdb: " <<  radg << endl;
		//cout << "RG diff for " <<  fit_file_name2  << ": " <<  RGdiff << endl;
		//cout << "RG diff normalize for " <<  fit_file_name2 << ": " <<  RGdiff_normalize << endl;
		
		
		  std::vector<RG_saxs> rgInfo;
		  rgInfo.clear();
		  RG_saxs rgData;
		  rgData.expRG = expRG;
		  rgData.modelRG = radg;
		  rgData.RGdiff = RGdiff;
		  rgData.RGdiff_normalize = RGdiff_normalize;
		  RG_fittingInfo.push_back(rgData);
	  
	  
		  
/*
      std::cout << pdb_files[i] << " " << dat_files[j]
               << " Chi = " << fp.get_chi() << " c1 = " << fp.get_c1()
                << " c2 = " << fp.get_c2()
                << " default chi = " << fp.get_default_chi() << std::endl;
      fp.set_pdb_file_name(pdb_files[i]);
      fp.set_profile_file_name(dat_files[j]);
      fp.set_mol_index(i);
      if (gnuplot_script) Gnuplot::print_fit_script(fp);
      fps.push_back(fp);
*/
    }
  }
	

	//cout << "Writing score to file! " << endl;	  
	// calculate profile functions, KL divergence, added by Jie 
	//cout << "pdb_prInfo.size(): " << pdb_prInfo.size() << endl;
	//cout << "exp_prInfo.size(): " << exp_prInfo.size() << endl;
	saxsInfo saxsScore;
	if(pdb_prInfo.size()!=0)
	{
		//sprintf(statScoreFile, "%s_score.txt",const_cast<char*>(jobname.c_str()));
		//cout << "Writing score to file: " << statScoreFile << endl;	  
		//ofstream statsaxs(statScoreFile);
		//if (!statsaxs.is_open()) {
		//	cout << "Error! folding statistics file can not open " << statScoreFile << endl;
		//	exit(0);
		//}
		
		// buffer statistics
		//char bufStat[1000];
		//sprintf(bufStat,"Exp_dat\tDecoy_name\tChi_score\tRG_exp\tRG_model\tRG_diff\tRG_diff_normalize\tKL_div\tKL_div_sym\tScore1\tScore2\tScore3\tScore4\tScore5\tScore6\tScore7\tScore8\tScore9\tScore10\tScore11");
		//statsaxs << bufStat << endl;
	
		for (int j = 0; j < pdb_prInfo.size(); j++) {
			for (int m = 0; m < exp_prInfo.size(); m++) {
				string fit_profile_array = pdb_fittingInfo[j+m];
				RG_saxs rgData = RG_fittingInfo[j+m];
				double chi_score = chi_fittingInfo[j+m];
				// calculate profile function used in SAXSTER http://www.cell.com/cms/attachment/2073830349/2068835456/mmc1.pdf
				std::string line;
				std::istringstream split(fit_profile_array);
				std::vector<fitting_profile> fittingdatInfo;
				fittingdatInfo.clear();
				while(std::getline(split, line, ';')) {
					//std::cout  << "line:" << line << endl;
					
					std::string line2;
					std::istringstream split2(line);
					fitting_profile datData;
					std::vector<double> vect;
					while(std::getline(split2, line2, ' ')) {
						//vect.push_back(std::strtod(line2.c_str() ));
						//cout << "--> "<< line2<<endl;
						vect.push_back(atof(line2.c_str()));
						
					}
					if(vect.size() !=4)
					{
						std::cout<<"The fitting data has incorrect column number!" <<std::endl;
						exit(-1); 
					}
					datData.q = vect[0];
					datData.expInt = vect[1];
					datData.modelInt = vect[2];
					datData.error = vect[3];
					fittingdatInfo.push_back(datData);
				}
				
				double fun_score1=0;
				double fun_score2=0;
				double fun_score2_nom=0;
				double fun_score2_denom=0;
				double fun_score3=0;
				double fun_score3_nom=0;
				double fun_score3_denom=0;
				double fun_score4=0;
				double fun_score5=0;
				double fun_score6=0;
				double fun_score6_nom=0;
				double fun_score6_denom=0;
				
				vector<double> fun_score8_expInt;
				vector<double> fun_score8_modelInt;
				for (unsigned int i = 0; i < fittingdatInfo.size(); i++) { 
					  //std::cout.precision(10);
					//std::cout << fittingdatInfo[i].q << " " << fittingdatInfo[i].expInt << " " << fittingdatInfo[i].modelInt<< " " << fittingdatInfo[i].error << std::endl;
					
					if(fittingdatInfo[i].expInt < 10e-10)
					{
						fittingdatInfo[i].expInt = 10e-10;
					}
					if(fittingdatInfo[i].modelInt < 10e-10)
					{
						fittingdatInfo[i].modelInt = 10e-10;
					}
						
					
					//function 1
					fun_score1 +=((fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt)*(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt))/(fittingdatInfo[i].error*fittingdatInfo[i].error);
					//function 2
					fun_score2_nom +=std::abs(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt);
					fun_score2_denom +=std::abs(fittingdatInfo[i].expInt);
					//function 3
					fun_score3_nom +=std::abs(log2(fittingdatInfo[i].expInt) - log2(fittingdatInfo[i].modelInt));
					fun_score3_denom +=std::abs(log2(fittingdatInfo[i].expInt));
					//function4
					fun_score4 += fittingdatInfo[i].q * ((fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt)*(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt));
					//function5
					fun_score5 += fittingdatInfo[i].q * ((log2(fittingdatInfo[i].expInt) - log2(fittingdatInfo[i].modelInt))*(log2(fittingdatInfo[i].expInt) - log2(fittingdatInfo[i].modelInt)));
					
					//function 6
					fun_score6_nom +=fittingdatInfo[i].q *fittingdatInfo[i].q *std::abs(fittingdatInfo[i].expInt - fittingdatInfo[i].modelInt);
					fun_score6_denom +=fittingdatInfo[i].q *fittingdatInfo[i].q*std::abs(fittingdatInfo[i].expInt);
					
					//function 8
					fun_score8_expInt.push_back(fittingdatInfo[i].expInt);
					fun_score8_modelInt.push_back(fittingdatInfo[i].modelInt);
				}
				fun_score1  /= fittingdatInfo.size();
				fun_score2 = fun_score2_nom /fun_score2_denom;
				fun_score3 = fun_score3_nom  / fun_score3_denom;
				fun_score6 = fun_score6_nom  / fun_score6_denom;
				
				double fun_score8 = log2(1-pearson(fun_score8_expInt,fun_score8_modelInt));
				double fun_score10 = log2(1-cosine(fun_score8_expInt,fun_score8_modelInt));
				double IntPCC = pearson(fun_score8_expInt,fun_score8_modelInt);
				double Intcosine = cosine(fun_score8_expInt,fun_score8_modelInt);
				double Inteuclidean = euclidean(fun_score8_expInt,fun_score8_modelInt);
				double Intjaccard = jaccard(fun_score8_expInt,fun_score8_modelInt);
				
			
				
				std::string exp_file_name = trim_extension(basename(const_cast<char*>(dat_files[m].c_str()))) +
			  ".dat";
			  std::string pdb_file_name = std::string(pdb_files[j]);
			  std::vector<pr_dist> e_prInfo=exp_prInfo[m];
			  std::vector<pr_dist> p_prInfo=pdb_prInfo[j];
			  int min_size=0;
			  if(e_prInfo.size() < p_prInfo.size())
			  {
				  min_size = e_prInfo.size();
			  }else{
				  min_size = p_prInfo.size();
			  }
			  double KL=0;
			  double KL_symetric=0;
			  double fun_score7=0;
			  vector<double> fun_score9_expPr;
			  vector<double> fun_score9_modelPr;
			  //std::cout << "Particle size for KL: "<< min_size << std::endl << std::endl;
			  for (int k = 0; k < min_size; k++) {
				if(e_prInfo[k].radius != p_prInfo[k].radius)
				{ 
					std::cout << "Error! the pr format is not correct!" << std::endl << std::endl;
					exit(-1);
				}
				if(e_prInfo[k].dist !=0 and p_prInfo[k].dist !=0)
				{
					if(e_prInfo[k].dist < 10e-10)
					{
						e_prInfo[k].dist = 10e-10;
					}
					if(p_prInfo[k].dist < 10e-10)
					{
						p_prInfo[k].dist = 10e-10;
					}
					KL += p_prInfo[k].dist *log2(p_prInfo[k].dist /e_prInfo[k].dist );
					KL_symetric += p_prInfo[k].dist *log2(p_prInfo[k].dist /e_prInfo[k].dist ) + e_prInfo[k].dist *log2(e_prInfo[k].dist /p_prInfo[k].dist );
					fun_score7 += (p_prInfo[k].dist - e_prInfo[k].dist)*(p_prInfo[k].dist - e_prInfo[k].dist);
					
					fun_score9_expPr.push_back(p_prInfo[k].dist);
					fun_score9_modelPr.push_back(e_prInfo[k].dist);
					
				}
			  }
			  double fun_score9 = log2(1-pearson(fun_score9_expPr,fun_score9_modelPr));
			  double fun_score11 = log2(1-cosine(fun_score9_expPr,fun_score9_modelPr));
			  double PrPCC = pearson(fun_score9_expPr,fun_score9_modelPr);
			  double Prcosine = cosine(fun_score9_expPr,fun_score9_modelPr);
			  double Preuclidean = euclidean(fun_score9_expPr,fun_score9_modelPr);
			  double Prjaccard = jaccard(fun_score9_expPr,fun_score9_modelPr);
				
			  if (KL  <0)
			  {
				  KL = 0.0;
			  }
			  if (KL_symetric  <0)
			  {
				  KL_symetric = 0.0;
			  }

			  
			  //char bufStat[1000];
			  //sprintf(bufStat,"%s\t%s\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f\t%8.3f",basename(const_cast<char*>(exp_file_name.c_str())),basename(const_cast<char*>(pdb_file_name.c_str())),chi_score,rgData.expRG,rgData.modelRG,rgData.RGdiff,rgData.RGdiff_normalize,KL,KL_symetric,fun_score1,fun_score2,fun_score3,fun_score4,fun_score5,fun_score6,fun_score7,fun_score8,fun_score9,fun_score10,fun_score11 );
			  //statsaxs << bufStat << endl;
			  
			 
			  saxsScore.expfile=exp_file_name;
			  saxsScore.pdfile=pdb_file_name;
			  saxsScore.chi_score=chi_score;
			  saxsScore.expRG=rgData.expRG;
			  saxsScore.modelRG=rgData.modelRG;
			  saxsScore.RGdiff=rgData.RGdiff;
			  saxsScore.RGdiff_normalize=rgData.RGdiff_normalize;
			  saxsScore.KL=KL;
			  saxsScore.KL_symetric=KL_symetric;
			  saxsScore.fun_score1=fun_score1;
			  saxsScore.fun_score2=fun_score2;
			  saxsScore.fun_score3=fun_score3;
			  saxsScore.fun_score4=fun_score4;
			  saxsScore.fun_score5=fun_score5;
			  saxsScore.fun_score6=fun_score6;
			  saxsScore.fun_score7=fun_score7;
			  saxsScore.fun_score8=fun_score8;
			  saxsScore.fun_score9=fun_score9;
			  saxsScore.fun_score10=fun_score10;
			  saxsScore.fun_score11=fun_score11;
		  }
		}
	}
  //std::sort(fps.begin(), fps.end(), FitParameters::compare_fit_parameters());
 /*
  if (pdb_files.size() > 1 && gnuplot_script) {
    Gnuplot::print_profile_script(pdb_files);
    if (dat_files.size() > 0) Gnuplot::print_fit_script(fps);
  }

  if (javascript) {
    if (dat_files.size() > 0) {
      Gnuplot::print_canvas_script(fps, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(fps, particles_vec, "jmoltable");
    } else {
      Gnuplot::print_canvas_script(pdb_files, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(pdb_files, particles_vec, "jmoltable");
    }
  }
 */
 
	//clock_t foxs_end = clock();
	//	double elapse_secs1 = double(foxs_end - foxs_start) / CLOCKS_PER_SEC;
		//cout << "Foxs finished within "<< elapse_secs1 << " sec!" << endl << endl;
  //return 0;
  // free memory, add by jie
  std::vector<FitParameters>().swap(fps);
  std::vector<std::vector<pr_dist> >().swap(pdb_prInfo);  // add by jie
  std::vector<std::string >().swap(pdb_fittingInfo);  // add by jie
  std::vector<RG_saxs>().swap(RG_fittingInfo);  // add by jie
  std::vector<double>().swap(chi_fittingInfo);  // add by jie
  std::vector<IMP::Particles>().swap(particles_vec);
  
  return saxsScore;
}


//int run_foxs(int foxsargc, char** foxsargv) {
double run_foxs(int foxsargc, char** foxsargv) {
  // output arguments
  // comment by Jie
  //for (int i = 0; i < foxsargc; i++) std::cerr << foxsargv[i] << " ";
  //std::cerr << std::endl;

  int profile_size = 500;
  float max_q = 0.0; // change after read
  float min_c1 = 0.99;
  float max_c1 = 1.05;
  float min_c2 = -0.5;
  float max_c2 = 2.0;
  bool heavy_atoms_only = true;
  bool residue_level = false;
  float background_adjustment_q = 0.0;
  bool use_offset = false;
  bool write_partial_profile = false;
  int multi_model_pdb = 1;
  bool vr_score = false;
  bool score_log = false;
  bool gnuplot_script = false;
  bool tofitting = false; 

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Any number of input PDBs and profiles is supported. \
Each PDB will be fitted against each profile.")
    ("version", "FoXS (IMP applications)\nCopyright 2007-2016 IMP Inventors.\n\
All rights reserved. \nLicense: GNU LGPL version 2.1 or later\n\
<http://gnu.org/licenses/lgpl.html>.\n\
Written by Dina Schneidman.")
    ("profile_size,s", po::value<int>(&profile_size)->default_value(500, "500"),
     "number of points in the profile")
    ("max_q,q", po::value<float>(&max_q)->default_value(0.5, "0.50"), "max q value")
    ("min_c1", po::value<float>(&min_c1)->default_value(0.99, "0.99"), "min c1 value")
    ("max_c1", po::value<float>(&max_c1)->default_value(1.05, "1.05"), "max c1 value")
    ("min_c2", po::value<float>(&min_c2)->default_value(-2.0, "-2.00"), "min c2 value")
    ("max_c2", po::value<float>(&max_c2)->default_value(4.0, "4.00"), "max c2 value")
    ("tofitting,n", "only default chi is calculated (default = false)") 
    ("hydrogens,h", "explicitly consider hydrogens in PDB files (default = false)")
    ("residues,r", "fast coarse grained calculation using CA atoms only (default = false)")
    ("background_q,b", po::value<float>(&background_adjustment_q)->default_value(0.0),
     "background adjustment, not used by default. if enabled, recommended q value is 0.2")
    ("offset,o", "use offset in fitting (default = false)")
    ("write-partial-profile,p", "write partial profile file (default = false)")
    ("multi-model-pdb,m", po::value<int>(&multi_model_pdb)->default_value(1),
     "1 - read the first MODEL only (default), \
2 - read each MODEL into a separate structure, \
3 - read all models into a single structure")
    ("volatility_ratio,v","calculate volatility ratio score (default = false)")
    ("score_log,l", "use log(intensity) in fitting and scoring (default = false)")
    ("gnuplot_script,g", "print gnuplot script for gnuplot viewing (default = false)");

  std::string form_factor_table_file;
  std::string beam_profile_file;
  bool ab_initio = false;
  bool vacuum = false;
  bool javascript = false;
  int chi_free = 0;
  float pr_dmax = 0.0;
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files")
    ("form_factor_table,f", po::value<std::string>(&form_factor_table_file),
     "ff table name")
    ("beam_profile", po::value<std::string>(&beam_profile_file),
     "beam profile file name for desmearing")
    ("ab_initio,a", "compute profile for a bead model with \
constant form factor (default = false)")
    ("vacuum", "compute profile in vacuum (default = false)")
    ("javascript,j",
     "output javascript for browser viewing of the results (default = false)")
    ("chi_free,x", po::value<int>(&chi_free)->default_value(0),
     "compute chi-free instead of chi, specify iteration number (default = 0)")
    ("pr_dmax", po::value<float>(&pr_dmax)->default_value(0.0, "0.0"),
     "Dmax value for P(r) calculation. P(r) is calculated only is pr_dmax > 0");

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::options_description visible(
      "Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ... ");
  visible.add(desc);

  po::positional_options_description p;
  p.add("input-files", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(foxsargc, foxsargv)
                .options(cmdline_options)
                .positional(p)
                .run(),
            vm);
  po::notify(vm);

  bool fit = true;
  std::vector<std::string> files, pdb_files, dat_files;
  if (vm.count("input-files")) {
    files = vm["input-files"].as<std::vector<std::string> >();
  }
  if (vm.count("help") || files.size() == 0) {
    std::cout << visible << "\n";
    return 0;
  }
  if (vm.count("tofitting")) tofitting = true; 
  if (vm.count("hydrogens")) heavy_atoms_only = false;
  if (vm.count("residues")) residue_level = true;
  if (vm.count("offset")) use_offset = true;
  if (vm.count("write-partial-profile")) write_partial_profile = true;
  if (vm.count("score_log")) score_log = true;
  if (vm.count("gnuplot_script")) gnuplot_script = true;

  // no water layer or fitting in ab initio mode for now
  if (vm.count("ab_initio")) {
    ab_initio = true;
    fit = false;
  }
  if (vm.count("vacuum")) {
    vacuum = true;
  }
  if (vm.count("javascript")) {
    javascript = true;
  }
  if (vm.count("volatility_ratio")) {
    vr_score = true;
  }

  if (multi_model_pdb != 1 && multi_model_pdb != 2 && multi_model_pdb != 3) {
    std::cerr << "Incorrect option for multi_model_pdb " << multi_model_pdb
              << std::endl;
    std::cerr << "Use 1 to read first MODEL only\n"
              << "    2 to read each MODEL into a separate structure,\n"
              << "    3 to read all models into a single structure\n";
    std::cerr << "Default value of 1 is used\n";
    multi_model_pdb = 1;
  }

  //IMP::benchmark::Profiler pp("prof_out");

  // determine form factor type
  FormFactorType ff_type = HEAVY_ATOMS;
  if (!heavy_atoms_only) ff_type = ALL_ATOMS;
  if (residue_level) ff_type = CA_ATOMS;

  // 1. read pdbs and profiles, prepare particles
  std::vector<IMP::Particles> particles_vec;
  Profiles exp_profiles;
  //cout << "Check here" << endl;
  read_files(files, pdb_files, dat_files, particles_vec, exp_profiles,
             residue_level, heavy_atoms_only, multi_model_pdb, max_q);

  if (background_adjustment_q > 0.0) {
    for (unsigned int i = 0; i < exp_profiles.size(); i++)
      exp_profiles[i]->background_adjust(background_adjustment_q);
  }

  if (exp_profiles.size() == 0 && !write_partial_profile) fit = false;

  if (max_q == 0.0) { // determine max_q
    if (exp_profiles.size() > 0) {
      for (unsigned int i = 0; i < exp_profiles.size(); i++)
        if (exp_profiles[i]->get_max_q() > max_q)
          max_q = exp_profiles[i]->get_max_q();
    } else {
      max_q = 0.5;
    }
  }
  float delta_q = max_q / profile_size;

  // read in or use default form factor table
  bool reciprocal = false;
  FormFactorTable* ft = NULL;
  if (form_factor_table_file.length() > 0) {
    // reciprocal space calculation, requires form factor file
    ft = new FormFactorTable(form_factor_table_file, 0.0, max_q, delta_q);
    reciprocal = true;
  } else {
    ft = get_default_form_factor_table();
  }

  // 2. compute profiles for input pdbs
  Profiles profiles;
  std::vector<FitParameters> fps;
  
  double chi_score = 0.0;
  for (unsigned int i = 0; i < particles_vec.size(); i++) {
    //std::cerr << "Computing profile for " << pdb_files[i] << " "
     //         << particles_vec[i].size() << " atoms " << std::endl; // comment by Jie
    IMP::Pointer<Profile> profile =
        compute_profile(particles_vec[i], 0.0, max_q, delta_q, ft, ff_type,
                        fit, fit, reciprocal, ab_initio, vacuum,
                        beam_profile_file);

    // save the profile
    profiles.push_back(profile);
    // write profile file
/*  comment by Jie
    std::string profile_file_name = std::string(pdb_files[i]) + ".dat";
    if (write_partial_profile)
      profile->write_partial_profiles(profile_file_name);
    else {  // write normal profile
      profile->add_errors();
      profile->write_SAXS_file(profile_file_name);
      if (gnuplot_script) Gnuplot::print_profile_script(pdb_files[i]);
    }

    // calculate P(r)
    if(pr_dmax > 0.0) {
      RadialDistributionFunction pr(0.5);
      profile->profile_2_distribution(pr, pr_dmax);
      pr.normalize();
      std::string pr_file_name = std::string(pdb_files[i]) + ".pr";
      std::ofstream pr_file(pr_file_name.c_str());
      pr.show(pr_file);
    }
*/


    // 3. fit experimental profiles
    for (unsigned int j = 0; j < dat_files.size(); j++) {
	  //cout << "datfile: " << j << " : " << dat_files[j] << endl;
      Profile* exp_saxs_profile = exp_profiles[j];
      std::string fit_file_name2 =
        trim_extension(pdb_files[i]) + "_" +
        trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) +
          ".dat";
     
    //Compute default score directly
	if(!tofitting)
	{		
		IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
		chi_score = pf->compute_score(profile, use_offset, fit_file_name2);
		//cout << "The no fitting score is: " << default_chi << endl; 
	}else{
	
      FitParameters fp;
      if (score_log) {
        IMP_NEW(ProfileFitter<ChiScoreLog>, pf, (exp_saxs_profile));
        fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                             use_offset, fit_file_name2);
      } else {
        if (vr_score) {
          IMP_NEW(ProfileFitter<RatioVolatilityScore>, pf, (exp_saxs_profile));
          fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                               use_offset, fit_file_name2);
        } else {
          IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
          fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                               use_offset, fit_file_name2);
          if (chi_free > 0) {
            float dmax = compute_max_distance(particles_vec[i]);
            unsigned int ns = IMP::algebra::get_rounded(
                           exp_saxs_profile->get_max_q() * dmax / IMP::PI);
            int K = chi_free;
            IMP_NEW(ChiFreeScore, cfs, (ns, K));
            cfs->set_was_used(true);
            // IMP_NEW(RatioVolatilityScore, rvs, ());
            // rvs->set_was_used(true);
            // resample the profile
            IMP_NEW(Profile, resampled_profile,
                    (exp_saxs_profile->get_min_q(), exp_saxs_profile->get_max_q(),
                     exp_saxs_profile->get_delta_q()));
            pf->resample(profile, resampled_profile);
            float chi_free =
              cfs->compute_score(exp_saxs_profile, resampled_profile);
            fp.set_chi(chi_free);
          }
        }
      }
	  chi_score = fp.get_chi();
	  double default_chi_score = fp.get_default_chi();
	  //cout << "The score is: " << chi_score << endl;
	  //cout << "The default score is: " << default_chi_score << endl;
	}

/*
      std::cout << pdb_files[i] << " " << dat_files[j]
               << " Chi = " << fp.get_chi() << " c1 = " << fp.get_c1()
                << " c2 = " << fp.get_c2()
                << " default chi = " << fp.get_default_chi() << std::endl;
      fp.set_pdb_file_name(pdb_files[i]);
      fp.set_profile_file_name(dat_files[j]);
      fp.set_mol_index(i);
      if (gnuplot_script) Gnuplot::print_fit_script(fp);
      fps.push_back(fp);
*/
    }
  }

  //std::sort(fps.begin(), fps.end(), FitParameters::compare_fit_parameters());
 /*
  if (pdb_files.size() > 1 && gnuplot_script) {
    Gnuplot::print_profile_script(pdb_files);
    if (dat_files.size() > 0) Gnuplot::print_fit_script(fps);
  }

  if (javascript) {
    if (dat_files.size() > 0) {
      Gnuplot::print_canvas_script(fps, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(fps, particles_vec, "jmoltable");
    } else {
      Gnuplot::print_canvas_script(pdb_files, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(pdb_files, particles_vec, "jmoltable");
    }
  }
 */
  //return 0;
  return chi_score;
}




//int run_foxs(int foxsargc, char** foxsargv) {
double run_foxs_saxs(int foxsargc, char** foxsargv, string inputString) {
  // output arguments
  // comment by Jie
  //for (int i = 0; i < foxsargc; i++) std::cerr << foxsargv[i] << " ";
  //std::cerr << std::endl;

  int profile_size = 500;
  float max_q = 0.0; // change after read
  float min_c1 = 0.99;
  float max_c1 = 1.05;
  float min_c2 = -0.5;
  float max_c2 = 2.0;
  bool heavy_atoms_only = true;
  bool residue_level = false;
  float background_adjustment_q = 0.0;
  bool use_offset = false;
  bool write_partial_profile = false;
  int multi_model_pdb = 1;
  bool vr_score = false;
  bool score_log = false;
  bool gnuplot_script = false;
  bool tofitting = false; 

  po::options_description desc("Options");
  desc.add_options()
    ("help", "Any number of input PDBs and profiles is supported. \
Each PDB will be fitted against each profile.")
    ("version", "FoXS (IMP applications)\nCopyright 2007-2016 IMP Inventors.\n\
All rights reserved. \nLicense: GNU LGPL version 2.1 or later\n\
<http://gnu.org/licenses/lgpl.html>.\n\
Written by Dina Schneidman.")
    ("profile_size,s", po::value<int>(&profile_size)->default_value(500, "500"),
     "number of points in the profile")
    ("max_q,q", po::value<float>(&max_q)->default_value(0.5, "0.50"), "max q value")
    ("min_c1", po::value<float>(&min_c1)->default_value(0.99, "0.99"), "min c1 value")
    ("max_c1", po::value<float>(&max_c1)->default_value(1.05, "1.05"), "max c1 value")
    ("min_c2", po::value<float>(&min_c2)->default_value(-2.0, "-2.00"), "min c2 value")
    ("max_c2", po::value<float>(&max_c2)->default_value(4.0, "4.00"), "max c2 value")
    ("tofitting,n", "only default chi is calculated (default = false)") 
    ("hydrogens,h", "explicitly consider hydrogens in PDB files (default = false)")
    ("residues,r", "fast coarse grained calculation using CA atoms only (default = false)")
    ("background_q,b", po::value<float>(&background_adjustment_q)->default_value(0.0),
     "background adjustment, not used by default. if enabled, recommended q value is 0.2")
    ("offset,o", "use offset in fitting (default = false)")
    ("write-partial-profile,p", "write partial profile file (default = false)")
    ("multi-model-pdb,m", po::value<int>(&multi_model_pdb)->default_value(1),
     "1 - read the first MODEL only (default), \
2 - read each MODEL into a separate structure, \
3 - read all models into a single structure")
    ("volatility_ratio,v","calculate volatility ratio score (default = false)")
    ("score_log,l", "use log(intensity) in fitting and scoring (default = false)")
    ("gnuplot_script,g", "print gnuplot script for gnuplot viewing (default = false)");

  std::string form_factor_table_file;
  std::string beam_profile_file;
  bool ab_initio = false;
  bool vacuum = false;
  bool javascript = false;
  int chi_free = 0;
  float pr_dmax = 0.0;
  po::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-files", po::value<std::vector<std::string> >(),
     "input PDB and profile files")
    ("form_factor_table,f", po::value<std::string>(&form_factor_table_file),
     "ff table name")
    ("beam_profile", po::value<std::string>(&beam_profile_file),
     "beam profile file name for desmearing")
    ("ab_initio,a", "compute profile for a bead model with \
constant form factor (default = false)")
    ("vacuum", "compute profile in vacuum (default = false)")
    ("javascript,j",
     "output javascript for browser viewing of the results (default = false)")
    ("chi_free,x", po::value<int>(&chi_free)->default_value(0),
     "compute chi-free instead of chi, specify iteration number (default = 0)")
    ("pr_dmax", po::value<float>(&pr_dmax)->default_value(0.0, "0.0"),
     "Dmax value for P(r) calculation. P(r) is calculated only is pr_dmax > 0");

  po::options_description cmdline_options;
  cmdline_options.add(desc).add(hidden);

  po::options_description visible(
      "Usage: <pdb_file1> <pdb_file2> ... <profile_file1> <profile_file2> ... ");
  visible.add(desc);

  po::positional_options_description p;
  p.add("input-files", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(foxsargc, foxsargv)
                .options(cmdline_options)
                .positional(p)
                .run(),
            vm);
  po::notify(vm);

  bool fit = true;
  std::vector<std::string> files, pdb_files, dat_files;
  if (vm.count("input-files")) {
    files = vm["input-files"].as<std::vector<std::string> >();
  }
  if (vm.count("help") || files.size() == 0) {
    std::cout << visible << "\n";
    return 0;
  }
  if (vm.count("tofitting")) tofitting = true; 
  if (vm.count("hydrogens")) heavy_atoms_only = false;
  if (vm.count("residues")) residue_level = true;
  if (vm.count("offset")) use_offset = true;
  if (vm.count("write-partial-profile")) write_partial_profile = true;
  if (vm.count("score_log")) score_log = true;
  if (vm.count("gnuplot_script")) gnuplot_script = true;

  // no water layer or fitting in ab initio mode for now
  if (vm.count("ab_initio")) {
    ab_initio = true;
    fit = false;
  }
  if (vm.count("vacuum")) {
    vacuum = true;
  }
  if (vm.count("javascript")) {
    javascript = true;
  }
  if (vm.count("volatility_ratio")) {
    vr_score = true;
  }

  if (multi_model_pdb != 1 && multi_model_pdb != 2 && multi_model_pdb != 3) {
    std::cerr << "Incorrect option for multi_model_pdb " << multi_model_pdb
              << std::endl;
    std::cerr << "Use 1 to read first MODEL only\n"
              << "    2 to read each MODEL into a separate structure,\n"
              << "    3 to read all models into a single structure\n";
    std::cerr << "Default value of 1 is used\n";
    multi_model_pdb = 1;
  }

  //IMP::benchmark::Profiler pp("prof_out");

  // determine form factor type
  FormFactorType ff_type = HEAVY_ATOMS;
  if (!heavy_atoms_only) ff_type = ALL_ATOMS;
  if (residue_level) ff_type = CA_ATOMS;

  // 1. read pdbs and profiles, prepare particles
  std::vector<IMP::Particles> particles_vec;
  Profiles exp_profiles;
  //cout << "Check here" << endl;
  //read_files(files, pdb_files, dat_files, particles_vec, exp_profiles,
  //           residue_level, heavy_atoms_only, multi_model_pdb, max_q);
  read_files_saxs(files, pdb_files, dat_files,inputString, particles_vec, exp_profiles,
             residue_level, heavy_atoms_only, multi_model_pdb, max_q);

  if (background_adjustment_q > 0.0) {
    for (unsigned int i = 0; i < exp_profiles.size(); i++)
      exp_profiles[i]->background_adjust(background_adjustment_q);
  }

  if (exp_profiles.size() == 0 && !write_partial_profile) fit = false;

  if (max_q == 0.0) { // determine max_q
    if (exp_profiles.size() > 0) {
      for (unsigned int i = 0; i < exp_profiles.size(); i++)
        if (exp_profiles[i]->get_max_q() > max_q)
          max_q = exp_profiles[i]->get_max_q();
    } else {
      max_q = 0.5;
    }
  }
  float delta_q = max_q / profile_size;

  // read in or use default form factor table
  bool reciprocal = false;
  FormFactorTable* ft = NULL;
  if (form_factor_table_file.length() > 0) {
    // reciprocal space calculation, requires form factor file
    ft = new FormFactorTable(form_factor_table_file, 0.0, max_q, delta_q);
    reciprocal = true;
  } else {
    ft = get_default_form_factor_table();
  }

  // 2. compute profiles for input pdbs
  Profiles profiles;
  std::vector<FitParameters> fps;
  
  double chi_score = 0.0;
  for (unsigned int i = 0; i < particles_vec.size(); i++) {
    //std::cerr << "Computing profile for " << pdb_files[i] << " "
     //         << particles_vec[i].size() << " atoms " << std::endl; // comment by Jie
    IMP::Pointer<Profile> profile =
        compute_profile(particles_vec[i], 0.0, max_q, delta_q, ft, ff_type,
                        fit, fit, reciprocal, ab_initio, vacuum,
                        beam_profile_file);

    // save the profile
    profiles.push_back(profile);
    // write profile file
/*  comment by Jie
    std::string profile_file_name = std::string(pdb_files[i]) + ".dat";
    if (write_partial_profile)
      profile->write_partial_profiles(profile_file_name);
    else {  // write normal profile
      profile->add_errors();
      profile->write_SAXS_file(profile_file_name);
      if (gnuplot_script) Gnuplot::print_profile_script(pdb_files[i]);
    }

    // calculate P(r)
    if(pr_dmax > 0.0) {
      RadialDistributionFunction pr(0.5);
      profile->profile_2_distribution(pr, pr_dmax);
      pr.normalize();
      std::string pr_file_name = std::string(pdb_files[i]) + ".pr";
      std::ofstream pr_file(pr_file_name.c_str());
      pr.show(pr_file);
    }
*/


    // 3. fit experimental profiles
    for (unsigned int j = 0; j < dat_files.size(); j++) {
	  //cout << "datfile: " << j << " : " << dat_files[j] << endl;
      Profile* exp_saxs_profile = exp_profiles[j];
      std::string fit_file_name2 =
        trim_extension(pdb_files[i]) + "_" +
        trim_extension(basename(const_cast<char*>(dat_files[j].c_str()))) +
          ".dat";
     
    //Compute default score directly
	if(!tofitting)
	{		
		IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
		chi_score = pf->compute_score(profile, use_offset, fit_file_name2);
		//cout << "The no fitting score is: " << default_chi << endl; 
	}else{
	
      FitParameters fp;
      if (score_log) {
        IMP_NEW(ProfileFitter<ChiScoreLog>, pf, (exp_saxs_profile));
        fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                             use_offset, fit_file_name2);
      } else {
        if (vr_score) {
          IMP_NEW(ProfileFitter<RatioVolatilityScore>, pf, (exp_saxs_profile));
          fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                               use_offset, fit_file_name2);
        } else {
          IMP_NEW(ProfileFitter<ChiScore>, pf, (exp_saxs_profile));
          fp = pf->fit_profile(profile, min_c1, max_c1, min_c2, max_c2,
                               use_offset, fit_file_name2);
          if (chi_free > 0) {
            float dmax = compute_max_distance(particles_vec[i]);
            unsigned int ns = IMP::algebra::get_rounded(
                           exp_saxs_profile->get_max_q() * dmax / IMP::PI);
            int K = chi_free;
            IMP_NEW(ChiFreeScore, cfs, (ns, K));
            cfs->set_was_used(true);
            // IMP_NEW(RatioVolatilityScore, rvs, ());
            // rvs->set_was_used(true);
            // resample the profile
            IMP_NEW(Profile, resampled_profile,
                    (exp_saxs_profile->get_min_q(), exp_saxs_profile->get_max_q(),
                     exp_saxs_profile->get_delta_q()));
            pf->resample(profile, resampled_profile);
            float chi_free =
              cfs->compute_score(exp_saxs_profile, resampled_profile);
            fp.set_chi(chi_free);
          }
        }
      }
	  chi_score = fp.get_chi();
	  double default_chi_score = fp.get_default_chi();
	  //cout << "The score is: " << chi_score << endl;
	  //cout << "The default score is: " << default_chi_score << endl;
	}

/*
      std::cout << pdb_files[i] << " " << dat_files[j]
               << " Chi = " << fp.get_chi() << " c1 = " << fp.get_c1()
                << " c2 = " << fp.get_c2()
                << " default chi = " << fp.get_default_chi() << std::endl;
      fp.set_pdb_file_name(pdb_files[i]);
      fp.set_profile_file_name(dat_files[j]);
      fp.set_mol_index(i);
      if (gnuplot_script) Gnuplot::print_fit_script(fp);
      fps.push_back(fp);
*/
    }
  }

  //std::sort(fps.begin(), fps.end(), FitParameters::compare_fit_parameters());
 /*
  if (pdb_files.size() > 1 && gnuplot_script) {
    Gnuplot::print_profile_script(pdb_files);
    if (dat_files.size() > 0) Gnuplot::print_fit_script(fps);
  }

  if (javascript) {
    if (dat_files.size() > 0) {
      Gnuplot::print_canvas_script(fps, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(fps, particles_vec, "jmoltable");
    } else {
      Gnuplot::print_canvas_script(pdb_files, JmolWriter::MAX_DISPLAY_NUM_);
      JmolWriter::prepare_jmol_script(pdb_files, particles_vec, "jmoltable");
    }
  }
 */
  //return 0;
  return chi_score;
}





/*************************************************************************
 * Name        : showEnergy_Print
 * Purpose     : shows (prints) the energy breakdown and total weighted energy
 * Arguments   : energyInfo
 * Return Type : void
 *************************************************************************/
void showEnergy_Print(energyInfo &e, vector<pdbInfo> &pdb, int decoyId, int startPos, int endPos) {
	

	char buf[1000];
	if (nFile) {
		// load native pdb if provided
		vector<pdbInfo> native;      
		double rmsd;
        loadPdb(natFile, native);    // load amino acid
		vector<int> aa;
		loadFasta(faFile, aa);
        if (aa.size() != native.size()) {
            cout << "Error! size mismatch in native structure" << endl;
            exit(0);
        }
		rmsd = getRmsd(pdb, native);
		cout << "get rmsd: " << rmsd << endl;
		//sprintf(buf,"UniCon3D.scoring: E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_chi_score_global = %8.3f   RMSD = %8.3f   E_total = %8.3f", e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,rmsd, e.total);  
		if(coolingDown)
		{
			//sprintf(buf,"UniCon3D.scoring(CoolingDown): Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f ( %.3f )    E_sc_bb = %8.3f ( %.3f )   E_bb_bb = %8.3f ( %.3f )   E_ri_rj = %8.3f ( %.3f )   E_chi_score_global = %8.3f ( %.3f )   RMSD = %8.3f  Min_Chi = %8.3f  E_total = %8.3f  ( Sampling cycles: %d looping_averaged_chi_dev: %.3f )", startPos, endPos, e.sc_sc,w_sc_sc, e.sc_bb,w_sc_bb, e.bb_bb,w_bb_bb, e.ri_rj,w_ri_rj, e.saxs_chi_global,w_saxs_chi,rmsd,saxs_chi_global_minima, e.total, total_sampling_cycle, looping_averaged_chi_dev);  
			sprintf(buf,"UniCon3D.scoring(CoolingDown): Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_chi_score_global = %8.3f ( %.3f )   E_KL = %8.3f ( %.3f )   E_score2 = %8.3f ( %.3f )   E_RG_normalize = %8.3f ( %.3f )   RMSD = %8.3f  Min_Chi = %8.3f  E_structure = %8.3f  E_saxs_energy = %8.3f  E_penalty = %8.3f  E_total = %8.3f", startPos, endPos, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,w_saxs_chi, e.saxs_KL,w_saxs_KL, e.saxs_score2,w_saxs_score2, e.saxs_RG_normalize,w_saxs_RG_normalize,rmsd,saxs_chi_global_minima, e.structure_energy, e.saxs_energy, e.saxs_penalty, e.total);  
		}else{
			//sprintf(buf,"UniCon3D.scoring: Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f ( %.3f )    E_sc_bb = %8.3f ( %.3f )   E_bb_bb = %8.3f ( %.3f )   E_ri_rj = %8.3f ( %.3f )   E_chi_score_global = %8.3f ( %.3f )   RMSD = %8.3f  Min_Chi = %8.3f  E_total = %8.3f  ( Sampling cycles: %d looping_averaged_chi_dev: %.3f )", startPos, endPos,  e.sc_sc,w_sc_sc, e.sc_bb,w_sc_bb, e.bb_bb,w_bb_bb, e.ri_rj,w_ri_rj, e.saxs_chi_global,w_saxs_chi,rmsd,saxs_chi_global_minima, e.total, total_sampling_cycle, looping_averaged_chi_dev);  
			sprintf(buf,"UniCon3D.scoring: Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_chi_score_global = %8.3f ( %.3f )   E_KL = %8.3f ( %.3f )   E_score2 = %8.3f ( %.3f )   E_RG_normalize = %8.3f ( %.3f )   RMSD = %8.3f  Min_Chi = %8.3f  E_structure = %8.3f  E_saxs_energy = %8.3f  E_penalty = %8.3f  E_total = %8.3f", startPos, endPos,  e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,w_saxs_chi, e.saxs_KL,w_saxs_KL, e.saxs_score2,w_saxs_score2, e.saxs_RG_normalize,w_saxs_RG_normalize,rmsd,saxs_chi_global_minima, e.structure_energy, e.saxs_energy, e.saxs_penalty, e.total);  
		}   
	}else{
		//sprintf(buf,"UniCon3D.scoring: E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_chi_score_global = %8.3f   E_total = %8.3f", e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global, e.total);  
 		if(coolingDown)
		{
			//sprintf(buf,"UniCon3D.scoring(CoolingDown):  Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f ( %.3f )    E_sc_bb = %8.3f ( %.3f )   E_bb_bb = %8.3f ( %.3f )   E_ri_rj = %8.3f ( %.3f )   E_chi_score_global = %8.3f ( %.3f )  Min_Chi = %8.3f  E_total = %8.3f  ( Sampling cycles: %d looping_averaged_chi_dev: %.3f )", startPos, endPos,  e.sc_sc,w_sc_sc, e.sc_bb,w_sc_bb, e.bb_bb,w_bb_bb, e.ri_rj,w_ri_rj, e.saxs_chi_global,w_saxs_chi,saxs_chi_global_minima, e.total, total_sampling_cycle, looping_averaged_chi_dev);  
			sprintf(buf,"UniCon3D.scoring(CoolingDown):  Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f    E_sc_bb = %8.3f  E_bb_bb = %8.3f   E_ri_rj = %8.3f    E_chi_score_global = %8.3f ( %.3f )   E_KL = %8.3f ( %.3f )   E_score2 = %8.3f ( %.3f )   E_RG_normalize = %8.3f ( %.3f )  Min_Chi = %8.3f  E_structure = %8.3f  E_saxs_energy = %8.3f  E_penalty = %8.3f  E_total = %8.3f", startPos, endPos,  e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,w_saxs_chi, e.saxs_KL,w_saxs_KL, e.saxs_score2,w_saxs_score2, e.saxs_RG_normalize,w_saxs_RG_normalize,saxs_chi_global_minima, e.structure_energy, e.saxs_energy, e.saxs_penalty, e.total);  
		}else{
			//sprintf(buf,"UniCon3D.scoring:  Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f ( %.3f )    E_sc_bb = %8.3f ( %.3f )   E_bb_bb = %8.3f ( %.3f )   E_ri_rj = %8.3f ( %.3f )   E_chi_score_global = %8.3f ( %.3f )  Min_Chi = %8.3f  E_total = %8.3f  ( Sampling cycles: %d looping_averaged_chi_dev: %.3f )",  startPos, endPos, e.sc_sc,w_sc_sc, e.sc_bb,w_sc_bb, e.bb_bb,w_bb_bb, e.ri_rj,w_ri_rj, e.saxs_chi_global,w_saxs_chi,saxs_chi_global_minima, e.total, total_sampling_cycle, looping_averaged_chi_dev);  
			sprintf(buf,"UniCon3D.scoring:  Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_chi_score_global = %8.3f ( %.3f )   E_KL = %8.3f ( %.3f )   E_score2 = %8.3f ( %.3f )   E_RG_normalize = %8.3f ( %.3f )  Min_Chi = %8.3f  E_structure = %8.3f  E_saxs_energy = %8.3f  E_penalty = %8.3f  E_total = %8.3f",  startPos, endPos, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,w_saxs_chi, e.saxs_KL,w_saxs_KL, e.saxs_score2,w_saxs_score2, e.saxs_RG_normalize,w_saxs_RG_normalize,saxs_chi_global_minima, e.structure_energy, e.saxs_energy, e.saxs_penalty, e.total);  
		}
	}
	cout << "energy done" << endl;
    printf("%s\r",buf);
    fflush(stdout);
	
	char simulatefile[1000] = "";
	int foldonsize = static_cast<int>(pdb.size());
	
		
	if(odir)
	{
		sprintf(simulatefile, "%s/simulate_decoy%d_%d_stats.txt", outputdir,decoyId, foldonsize );

	}else{
		sprintf(simulatefile, "simulate_decoy%d_%d_stats.txt",decoyId, foldonsize );

	}	
	
    
	ofstream ofs;
	ofs.open (simulatefile, ofstream::out | ofstream::app);
	ofs << buf << endl;  
	ofs.close();
	
	 
    //cout << buf << endl;
} 

/*************************************************************************
 * Name        : showEnergy_Print
 * Purpose     : shows (prints) the energy breakdown and total weighted energy
 * Arguments   : energyInfo
 * Return Type : void
 *************************************************************************/
void showEnergy_Print_regularize(energyInfo &e, vector<pdbInfo> &pdb, int decoyId, int startPos, int endPos) {
	

	char buf[1000];
	if (nFile) {
		// load native pdb if provided
		vector<pdbInfo> native;      
		double rmsd;
        loadPdb(natFile, native);    // load amino acid
		vector<int> aa;
		loadFasta(faFile, aa);
        if (aa.size() != native.size()) {
            cout << "Error! size mismatch in native structure" << endl;
            exit(0);
        }
		rmsd = getRmsd(pdb, native);
		sprintf(buf,"UniCon3D.scoring(Regularization): Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f    E_sc_bb = %8.3f   E_bb_bb = %8.3f   E_ri_rj = %8.3f   E_chi_score_global = %8.3f ( %.3f )   E_KL = %8.3f ( %.3f )   E_score2 = %8.3f ( %.3f )   E_RG_normalize = %8.3f ( %.3f )   E_saxs_penalty = %8.3f ( %.3f )    E_structure = %8.3f    E_saxs_energy = %8.3f    RMSD = %8.3f  Min_Chi = %8.3f  E_total = %8.3f", startPos, endPos, e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,w_saxs_chi, e.saxs_KL,w_saxs_KL, e.saxs_score2,w_saxs_score2, e.saxs_RG_normalize,w_saxs_RG_normalize, e.saxs_penalty,w_saxs_chi_penalty, e.structure_energy, e.saxs_energy,rmsd,saxs_chi_global_minima, e.total);  
		 
	}else{
		sprintf(buf,"UniCon3D.scoring(Regularization): Sample_start = %d Sample_end = %d  E_sc_sc = %8.3f    E_sc_bb = %8.3f  E_bb_bb = %8.3f   E_ri_rj = %8.3f    E_chi_score_global = %8.3f ( %.3f )   E_KL = %8.3f ( %.3f )   E_score2 = %8.3f ( %.3f )   E_RG_normalize = %8.3f ( %.3f )   E_saxs_penalty = %8.3f ( %.3f )    E_structure = %8.3f    E_saxs_energy = %8.3f    Min_Chi = %8.3f  E_total = %8.3f", startPos, endPos,  e.sc_sc, e.sc_bb, e.bb_bb, e.ri_rj, e.saxs_chi_global,w_saxs_chi, e.saxs_KL,w_saxs_KL, e.saxs_score2,w_saxs_score2, e.saxs_RG_normalize,w_saxs_RG_normalize, e.saxs_penalty,w_saxs_chi_penalty, e.structure_energy, e.saxs_energy,saxs_chi_global_minima, e.total);  
		
	}
    printf("%s\r",buf);
    fflush(stdout);
	
	char simulatefile[1000] = "";
	int foldonsize = static_cast<int>(pdb.size());
	
		
	if(odir)
	{
		sprintf(simulatefile, "%s/simulate_decoy%d_%d_stats.txt", outputdir,decoyId, foldonsize );

	}else{
		sprintf(simulatefile, "simulate_decoy%d_%d_stats.txt",decoyId, foldonsize );

	}	
	
    
	ofstream ofs;
	ofs.open (simulatefile, ofstream::out | ofstream::app);
	ofs << buf << endl;  
	ofs.close();
	
	 
    //cout << buf << endl;
} 

/*************************************************************************
 * Name        : acceptMove
 * Purpose     : accept or reject a move
 * Arguments   : energyInfo &prev, energyInfo &curr, double temperature
 * Return Type : bool
 *************************************************************************/
bool acceptMove(energyInfo &prev, energyInfo &curr, double temperature) {
    double diff = prev.total - curr.total;
    if (diff > 0)
        return true;
    else {
        double e = exp(diff / (temperature * BOLTZMANN_CONSTANT));
        double r = uniDblGen();
        if (e > r)
            return true;
        else
            return false;
    }
}

bool acceptMove_jump(energyInfo &prev, energyInfo &curr, double temperature, double weight) {
    double diff = weight*prev.total - curr.total;
    if (diff > 0)
        return true;
    else {
        double e = exp(diff / (temperature * BOLTZMANN_CONSTANT));
        double r = uniDblGen();
        if (e > r)
            return true;
        else
            return false;
    }
}

void removeTempFile (char *Prefix){
	FILE * stream;
	char buffer[1028];	
	string cmd = "rm ";
	cmd.append(Prefix);
	cmd.append("*");
    stream = popen(cmd.c_str(), "r");
	if (stream == NULL) 
	{ 
        cout << "Failed to remove " << cmd << endl; 
		exit(0); 
	 
    }else{
		while (fgets(buffer, 1028, stream) != NULL) 
		{ 
			//printf(buffer); 
			fflush(stdout);
		} 
	}
	pclose(stream);

}
/*************************************************************************
 * Name        : writePdb
 * Purpose     : writes a pdb file in  united residue representation into a file
 * Arguments   : vector<pdbInfo> &pdb, char * filename
 * Return Type : void
 *************************************************************************/
void writePdb(vector<pdbInfo> &pdb, char * filename) {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cout << "Error! output pdb file can not open " << filename << endl;
        return;
    }
    
    for (int i=0; i < pdb.size(); i++) {
        // buffer CA atoms
        char bufCa[1000];
        sprintf(bufCa,"ATOM %6d  CA  %-3s %5d    %8.3f%8.3f%8.3f  1.00  0.00", 2 * i + 1, seq3[pdb[i].aa].c_str(), pdb[i].id, pdb[i].ca.x, pdb[i].ca.y, pdb[i].ca.z);
        // buffer SC atoms
        char bufSc[1000];
        sprintf(bufSc,"ATOM %6d  SC  %-3s %5d    %8.3f%8.3f%8.3f  1.00  0.00", 2 * i + 2, seq3[pdb[i].aa].c_str(), pdb[i].id, pdb[i].sc.x, pdb[i].sc.y, pdb[i].sc.z);
        // write CA and SC atoms to file
        fout << bufCa << endl;
        fout << bufSc << endl;
    }
}



/*************************************************************************
 * Name        : Pdb2String
 * Purpose     : writes a pdb file in  united residue representation into a file
 * Arguments   : vector<pdbInfo> &pdb, char * filename
 * Return Type : void
 *************************************************************************/
void Pdb2String(vector<pdbInfo> &pdb, string &inputString) {
    for (int i=0; i < pdb.size(); i++) {
        // buffer CA atoms
        char bufCa[1000];
        sprintf(bufCa,"ATOM %6d  CA  %-3s %5d    %8.3f%8.3f%8.3f  1.00  0.00;", 2 * i + 1, seq3[pdb[i].aa].c_str(), pdb[i].id, pdb[i].ca.x, pdb[i].ca.y, pdb[i].ca.z);
        // buffer SC atoms
        char bufSc[1000];
        sprintf(bufSc,"ATOM %6d  SC  %-3s %5d    %8.3f%8.3f%8.3f  1.00  0.00;", 2 * i + 2, seq3[pdb[i].aa].c_str(), pdb[i].id, pdb[i].sc.x, pdb[i].sc.y, pdb[i].sc.z);
        // write CA and SC atoms to file
        inputString.append(bufCa);
        inputString.append(bufSc); 
    }
}





/*************************************************************************
 * Name        : runPulcha
 * Purpose     : run pulcha on pose structure to get heavy atoms
 * Arguments   : char * filename, char * filename
 * Return Type : void
 *************************************************************************/
void runPulcha(char * pdbTempFile, char * pdbTempFile_pulchra) {
	
	FILE * stream;
	char buffer[1028]; 
	string cmd = "/home/jh7x3/tools/pulchra_306/pulchra  -g -v   ";
	cmd.append(pdbTempFile);
	cmd.append(" 2>&1");
	//cmd.append("  &>  PULCHRA_");
	//cmd.append(pdbTempFile_prefix);
	//cmd.append(".log");
	//cmd.append(foxlog);
    //cout << "command: " << cmd << endl;
	
    stream = popen(cmd.c_str(), "r");
	if (stream == NULL) 
	{ 
        cout << "Failed to run pulcha with " << cmd << endl; 
		exit(0); 
	 
    }else{
		while (fgets(buffer, 1028, stream) != NULL) 
		{ 
			//printf(buffer); 
			fflush(stdout);
		} 
		ifstream fin( pdbTempFile_pulchra);
		if( fin.fail() ) {
			cout << endl;
			cout << "Error! Pulcha file " << pdbTempFile_pulchra << " not present" << endl << endl;
			exit(0);
		}
		fin.close();
	}
	pclose(stream);	
}



/*************************************************************************
 * Name        : runPulcha
 * Purpose     : run pulcha on pose structure to get heavy atoms
 * Arguments   : char * filename, char * filename
 * Return Type : void
 *************************************************************************/
void runPulcha2(char * pdbTempFile, char * pdbTempFile_pulchra) {
	
	char** Files;
	Files = new char*[1000];// initialize the double pointer
	Files[0]=new char[1000];// initialize 1st char*, with capacity of 4 chars
	Files[1]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[2]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
		
	strcpy(Files[0],"pulchra");//copy some data to 1st string
	strcpy(Files[1],"-g");//copy some data to 2nd string
	strcpy(Files[2],pdbTempFile);//copy some data to 2nd string
	ifstream fin1( Files[2]);
	if( fin1.fail() ) {
		cout << endl;
		cout << "Error! Pulcha file " << Files[2] << " not present" << endl << endl;
		exit(0);
	}	
	fin1.close();
	int pulcha_back = run_Pulcha_int(3,Files);
	ifstream fin( pdbTempFile_pulchra);
	if( fin.fail() ) {
		cout << endl;
		cout << "Error! Pulcha file " << pdbTempFile_pulchra << " not present" << endl << endl;
		exit(0);
	}
	fin.close();
}


/*************************************************************************
 * Name        : runPulcha
 * Purpose     : run pulcha on pose structure to get heavy atoms
 * Arguments   : char * filename, char * filename
 * Return Type : void
 *************************************************************************/
void runPulcha3(char * pdbTempFile, char * pdbTempFile_pulchra, string pdbString, string &pdbString_pulchar, int savefile) {
	
	char** Files;
	Files = new char*[1000];// initialize the double pointer
	Files[0]=new char[1000];// initialize 1st char*, with capacity of 4 chars
	Files[1]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
	Files[2]=new char[1000];// initialize 2nd char*, with capacity of 5 chars
		
	strcpy(Files[0],"pulchra");//copy some data to 1st string
	strcpy(Files[1],"-g");//copy some data to 2nd string
	strcpy(Files[2],pdbTempFile);//copy some data to 2nd string
	//ifstream fin1( Files[2]);
	//if( fin1.fail() ) {
	//	cout << endl;
	//	cout << "Error! Pulcha file " << Files[2] << " not present" << endl << endl;
	//	exit(0);
	//}	
	//fin1.close();
	int pulcha_back = run_Pulcha_int2(3,Files,pdbString, pdbString_pulchar,savefile);
	//int pulcha_back = run_Pulcha_int(3,Files);
	delete[] Files;
			
	/*		
	cout << "in pulchar3 finish" << endl;
	ifstream fin( pdbTempFile_pulchra);
	if( fin.fail() ) {
		cout << endl;
		cout << "Error! Pulcha file " << pdbTempFile_pulchra << " not present" << endl << endl;
		exit(0);
	}
	fin.close();
	
	string line, str;
    ifstream fin2 (pdbTempFile_pulchra);
    if (fin2.is_open()) {
        while ( fin2.good() ) {
            getline(fin2, line);
			pdbString_pulchar.append(line);
		}
	}
	fin2.close();
	*/
}

/************************************************************************
 * Name        : run_Pulcha_int
 * Purpose     : run pulcha on pose structure to get heavy atoms
 * Arguments   : char * filename, char * filename
 * Return Type : void
 *************************************************************************/
int run_Pulcha_int(int pulchaargc, char** pulchaargv) {


	// initialize pulchar default parameter 
	_VERBOSE = 0;
	_BB_REARRANGE = 1;
	_BB_OPTIMIZE = 0;
	_CA_OPTIMIZE = 1;
	_CA_RANDOM = 0;
	_CA_ITER = 100;
	_CA_TRAJECTORY = 0;
	_CISPRO = 0;
	_CHIRAL = 1;
	_CENTER_CHAIN = 0;
	_REBUILD_BB = 1;
	_REBUILD_SC = 1;
	_REBUILD_H = 0;
	_PDB_SG = 0;
	_TIME_SEED = 0;
	_XVOLUME = 1;
	_XVOL_ITER = 3;
	_PRESERVE = 1;
	_CA_START_DIST = 3.0;
	_CA_XVOL_DIST = 3.5;
	_SG_XVOL_DIST = 1.6;


	CA_K=10.0;
	CA_ANGLE_K=20.0;
	CA_START_K=0.01;
	CA_XVOL_K=10.00;



	RBINS = NULL;
	X_COORDS = NULL;
	C_ALPHA = NULL;
	// initialize finished


  int i, j, next;
  char buf[100];
  char *name=NULL, *ini_name=NULL;
  char *ptr, out_name[1000];
  real f;
  mol_type *mol;
  struct timeb time0, time1;

    for (i=1; i<pulchaargc; i++) {
		
      if (pulchaargv[i][0]=='-') {
        next=0;
        for (j=1; j<(int)strlen(pulchaargv[i]); j++) {
          switch(pulchaargv[i][j]) {
            case 'v': _VERBOSE=1; break;
            case 'c': _CA_OPTIMIZE=0; break;
            case 'e': _BB_REARRANGE=1; break;
            case 'r': _CA_RANDOM=1; break;
            case 'z': _CHIRAL=0; break;
            case 't': _CA_TRAJECTORY=1; break;
            case 'n': _CENTER_CHAIN=1; break;
            case 'b': _REBUILD_BB=0; break;
            case 's': _REBUILD_SC=0; break;
            case 'i': ini_name = pulchaargv[++i]; next=1; break;
            case 'g': _PDB_SG=1; break;
            case 'x': _TIME_SEED=1; break;
            case 'o': _XVOLUME=0; break;
            case 'h': _REBUILD_H=0; break;
            case 'q': _BB_OPTIMIZE=1; break;
            case 'p': _CISPRO=1; break;
            case 'f': _PRESERVE=1; break;
            case 'u':
              if (sscanf(pulchaargv[++i],"%lf",&f)==1) {
                _CA_START_DIST = f;
              }
              next=1;
            break;
            default: {
              printf("Unknown option: %c\n", pulchaargv[i][j]);
              return -1;
            }
          }
          if (next) break;
        }
      } else {
        if (!name) name=pulchaargv[i];
      }
    }

    if (!name) {
      printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);
      printf("Usage: %s [options] <pdb_file>\n", pulchaargv[0]);
      printf("The program default input is a PDB file.\n");
      printf("Output file <pdb_file.rebuild.pdb> will be created as a result.\n");
      printf("Valid options are:\n\n");
      printf("  -v : verbose output (default: off)\n");
      printf("  -n : center chain (default: off)\n");
      printf("  -x : time-seed random number generator (default: off)\n");
      printf("  -g : use PDBSG as an input format (CA=C-alpha, SC or CM=side chain c.m.)\n\n");
      printf("  -c : skip C-alpha positions optimization (default: on)\n");
      printf("  -p : detect cis-prolins (default: off)\n");
      printf("  -r : start from a random chain (default: off)\n");
      printf("  -i pdbfile : read the initial C-alpha coordinates from a PDB file\n");
      printf("  -t : save chain optimization trajectory to file <pdb_file.pdb.trajectory>\n");
      printf("  -u value : maximum shift from the restraint coordinates (default: 0.5A)\n\n");
      printf("  -e : rearrange backbone atoms (C, O are output after side chain) (default: off)\n");
      printf("  -f : preserve initial coordinates (default: off, implies '-c' on and '-n' off)\n");

#ifdef COMPILE_BB
      printf("  -b : skip backbone reconstruction (default: on)\n");
      printf("  -q : optimize backbone hydrogen bonds pattern (default: off)\n");
      printf("  -h : outputs hydrogen atoms (default: off)\n");
#endif

#ifdef COMPILE_ROT
      printf("  -s : skip side chains reconstruction (default: on)\n");
      printf("  -o : don't attempt to fix excluded volume conflicts (default: on)\n");
      printf("  -z : don't check amino acid chirality (default: on)\n");
#endif
	  printf("\n");
      return -1;
    }

    for (i=0; i<255; i++) /* prepare hash table*/
      AA_NUMS[i] = 20; /* dummy aa code */
    for (i=0; i<20; i++)
      AA_NUMS[(int)SHORT_AA_NAMES[i]] = i;

    setbuf(stdout,0);

    if (_TIME_SEED) srand(time(NULL)); else srand(1234);

    if (_VERBOSE) printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);

    ftime(&time0);

	
    chain = new_mol();

    if (read_pdb_file(name,chain,"chain")==FILE_NOT_FOUND) {
      if (_VERBOSE) printf("Can't read the input file!\n");
      return -1;
    }

    if (_VERBOSE) printf("%d residua read.\n", chain->nres);

    //if (_PRESERVE) printf("Initial coordinates will be preserved.\n");

    chain_length = chain->nres;

    if (_CA_OPTIMIZE && !_PRESERVE) {
      snprintf(out_name,1000,"%s.tra",name);
      ca_optimize(out_name, ini_name);
    }

#ifdef COMPILE_BB
    if (_REBUILD_BB) {
      rebuild_backbone();
      if (_BB_OPTIMIZE) {
        optimize_backbone(chain);
      }
    }
#endif

#ifdef COMPILE_ROT
    if (_REBUILD_SC) {
      rebuild_sidechains();
      if (_XVOLUME)
        optimize_exvol();
      if (_CHIRAL)
        chirality_check();
    }
#endif

    if (_CENTER_CHAIN) {
      center_chain(chain);
    }


    if (_BB_REARRANGE) {
      if (_VERBOSE) printf("Rearranging backbone atoms...\n", out_name);
    }

    ptr = strstr(name,".pdb");
    if (ptr) ptr[0]=0;
      
    snprintf(out_name,1000,"%s.rebuilt.pdb",name);
    
    if (_VERBOSE) printf("Writing output file %s...\n", out_name);
    write_pdb(out_name, chain);

    ftime(&time1);

    if (_VERBOSE) printf("Done. Reconstruction finished in %.3f s.\n", (real)0.001*(1000.*(time1.time-time0.time)+(time1.millitm-time0.millitm)));
  return 0;
}



/*************************************************************************
 * Name        : run_Pulcha_int
 * Purpose     : run pulcha on pose structure to get heavy atoms
 * Arguments   : char * filename, char * filename
 * Return Type : void
 *************************************************************************/
int run_Pulcha_int2(int pulchaargc, char** pulchaargv,string pdbString, string &pdbString_pulchar, int savefile) {


	// initialize pulchar default parameter 
	_VERBOSE = 0;
	_BB_REARRANGE = 1;
	_BB_OPTIMIZE = 0;
	_CA_OPTIMIZE = 1;
	_CA_RANDOM = 0;
	_CA_ITER = 100;
	_CA_TRAJECTORY = 0;
	_CISPRO = 0;
	_CHIRAL = 1;
	_CENTER_CHAIN = 0;
	_REBUILD_BB = 1;
	_REBUILD_SC = 1;
	_REBUILD_H = 0;
	_PDB_SG = 0;
	_TIME_SEED = 0;
	_XVOLUME = 1;
	_XVOL_ITER = 3;
	_PRESERVE = 1;
	_CA_START_DIST = 3.0;
	_CA_XVOL_DIST = 3.5;
	_SG_XVOL_DIST = 1.6;


	CA_K=10.0;
	CA_ANGLE_K=20.0;
	CA_START_K=0.01;
	CA_XVOL_K=10.00;



	RBINS = NULL;
	X_COORDS = NULL;
	C_ALPHA = NULL;
	// initialize finished

  int i, j, next;
  char buf[100];
  char *name=NULL, *ini_name=NULL;
  char *ptr, out_name[1000];
  real f;
  mol_type *mol;
  struct timeb time0, time1;

    for (i=1; i<pulchaargc; i++) {
      if (pulchaargv[i][0]=='-') {
        next=0;
        for (j=1; j<(int)strlen(pulchaargv[i]); j++) {
          switch(pulchaargv[i][j]) {
            case 'v': _VERBOSE=1; break;
            case 'c': _CA_OPTIMIZE=0; break;
            case 'e': _BB_REARRANGE=1; break;
            case 'r': _CA_RANDOM=1; break;
            case 'z': _CHIRAL=0; break;
            case 't': _CA_TRAJECTORY=1; break;
            case 'n': _CENTER_CHAIN=1; break;
            case 'b': _REBUILD_BB=0; break;
            case 's': _REBUILD_SC=0; break;
            case 'i': ini_name = pulchaargv[++i]; next=1; break;
            case 'g': _PDB_SG=1; break;
            case 'x': _TIME_SEED=1; break;
            case 'o': _XVOLUME=0; break;
            case 'h': _REBUILD_H=0; break;
            case 'q': _BB_OPTIMIZE=1; break;
            case 'p': _CISPRO=1; break;
            case 'f': _PRESERVE=1; break;
            case 'u':
              if (sscanf(pulchaargv[++i],"%lf",&f)==1) {
                _CA_START_DIST = f;
              }
              next=1;
            break;
            default: {
              printf("Unknown option: %c\n", pulchaargv[i][j]);
              return -1;
            }
          }
          if (next) break;
        }
      } else {
        if (!name) name=pulchaargv[i];
      }
    }

    if (!name) {
      printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);
      printf("Usage: %s [options] <pdb_file>\n", pulchaargv[0]);
      printf("The program default input is a PDB file.\n");
      printf("Output file <pdb_file.rebuild.pdb> will be created as a result.\n");
      printf("Valid options are:\n\n");
      printf("  -v : verbose output (default: off)\n");
      printf("  -n : center chain (default: off)\n");
      printf("  -x : time-seed random number generator (default: off)\n");
      printf("  -g : use PDBSG as an input format (CA=C-alpha, SC or CM=side chain c.m.)\n\n");
      printf("  -c : skip C-alpha positions optimization (default: on)\n");
      printf("  -p : detect cis-prolins (default: off)\n");
      printf("  -r : start from a random chain (default: off)\n");
      printf("  -i pdbfile : read the initial C-alpha coordinates from a PDB file\n");
      printf("  -t : save chain optimization trajectory to file <pdb_file.pdb.trajectory>\n");
      printf("  -u value : maximum shift from the restraint coordinates (default: 0.5A)\n\n");
      printf("  -e : rearrange backbone atoms (C, O are output after side chain) (default: off)\n");
      printf("  -f : preserve initial coordinates (default: off, implies '-c' on and '-n' off)\n");

#ifdef COMPILE_BB
      printf("  -b : skip backbone reconstruction (default: on)\n");
      printf("  -q : optimize backbone hydrogen bonds pattern (default: off)\n");
      printf("  -h : outputs hydrogen atoms (default: off)\n");
#endif

#ifdef COMPILE_ROT
      printf("  -s : skip side chains reconstruction (default: on)\n");
      printf("  -o : don't attempt to fix excluded volume conflicts (default: on)\n");
      printf("  -z : don't check amino acid chirality (default: on)\n");
#endif
	  printf("\n");
      return -1;
    }

    for (i=0; i<255; i++) /* prepare hash table*/
      AA_NUMS[i] = 20; /* dummy aa code */
    for (i=0; i<20; i++)
      AA_NUMS[(int)SHORT_AA_NAMES[i]] = i;

    setbuf(stdout,0);

    if (_TIME_SEED) srand(time(NULL)); else srand(1234);

    if (_VERBOSE) printf("PULCHRA Protein Chain Restoration Algorithm version %4.2f\n", PULCHRA_VERSION);

    ftime(&time0);

	
    chain = new_mol();

    if (read_pdb_file_saxs(name,chain,"chain",pdbString)==FILE_NOT_FOUND) {
      if (_VERBOSE) printf("Can't read the input file!\n");
      return -1;
    }

    if (_VERBOSE) printf("%d residua read.\n", chain->nres);

    //if (_PRESERVE) printf("Initial coordinates will be preserved.\n");

    chain_length = chain->nres;

    if (_CA_OPTIMIZE && !_PRESERVE) {
      snprintf(out_name,1000,"%s.tra",name);
      ca_optimize(out_name, ini_name);
    }

#ifdef COMPILE_BB
    if (_REBUILD_BB) {
      rebuild_backbone();
      if (_BB_OPTIMIZE) {
        optimize_backbone(chain);
      }
    }
#endif

#ifdef COMPILE_ROT
    if (_REBUILD_SC) {
      rebuild_sidechains();
      if (_XVOLUME)
        optimize_exvol();
      if (_CHIRAL)
        chirality_check();
    }
#endif

    if (_CENTER_CHAIN) {
      center_chain(chain);
    }


    if (_BB_REARRANGE) {
      if (_VERBOSE) printf("Rearranging backbone atoms...\n", out_name);
    }

    ptr = strstr(name,".pdb");
    if (ptr) ptr[0]=0;
      
    snprintf(out_name,1000,"%s.rebuilt.pdb",name);
    
    if (_VERBOSE) printf("Writing output file %s...\n", out_name);
    if(savefile)
	{
		write_pdb(out_name, chain);
	}
	write_pdb_insideSAXS(out_name,pdbString_pulchar, chain);
	//std::string line;
	//std::istringstream split(pdbString_pulchar);
	//while(std::getline(split, line, ';')) {
		//std::cout << "chi inside string: " << line << endl;
	//}
	//std::string line2;
	//std::istringstream split2(pdbString);
	//while(std::getline(split2, line2, ';')) {
	//	std::cout  << line2 << endl;
	//}
    ftime(&time1);

    if (_VERBOSE) printf("Done. Reconstruction finished in %.3f s.\n", (real)0.001*(1000.*(time1.time-time0.time)+(time1.millitm-time0.millitm)));
  return 0;
}

//// Define pulcha functions 


void add_atom(atom_type* atomlist, atom_type* newatom)
{
  atom_type *tmpatom;

    if (!atomlist)
      atomlist=newatom;
    else {
      tmpatom=atomlist->next;
      atomlist->next=newatom;
      newatom->prev=atomlist;
      newatom->next=tmpatom;
      if (tmpatom) tmpatom->prev=newatom;
    }
}

void add_res(res_type* reslist, res_type* newres)
{
  res_type *tmpres;

    if (!reslist)
      reslist=newres;
    else {
      tmpres=reslist->next;
      reslist->next=newres;
      newres->prev=reslist;
      newres->next=tmpres;
      if (tmpres) tmpres->prev=newres;
    }
}

void add_mol(mol_type* mollist, mol_type* newmol)
{
  mol_type *tmpmol;

    if (!mollist)
      mollist=newmol;
    else {
      tmpmol=mollist->next;
      mollist->next=newmol;
      newmol->prev=mollist;
      newmol->next=tmpmol;
      if (tmpmol) tmpmol->prev=newmol;
    }
}

void delete_atom(atom_type* atom)
{
  atom_type *tmpatom;

    if (atom->prev) atom->prev->next=atom->next;
    if (atom->next) atom->next->prev=atom->prev;
    if (atom->name) free(atom->name);
    free(atom);
    atom=NULL;
}

void delete_res(res_type* res)
{
  res_type *tmpres;
  atom_type *tmpatom;

    if (res->prev) res->prev->next=res->next;
    if (res->next) res->next->prev=res->prev;
    if (res->name) free(res->name);
    if (res->atoms) {
      while (res->atoms) {
        tmpatom = res->atoms->next;
        delete_atom(res->atoms);
        res->atoms=tmpatom;
      }
    }
    free(res);
    res=NULL;
}

void delete_mol(mol_type* mol)
{
  mol_type *tmpmol;
  res_type *tmpres;
  int i;

    if (mol->prev) mol->prev->next=mol->next;
    if (mol->next) mol->next->prev=mol->prev;
    if (mol->name) free(mol->name);
    if (mol->residua) {
      while (mol->residua) {
        tmpres = mol->residua->next;
        delete_res(mol->residua);
        mol->residua=tmpres;
      }
    }
    if (mol->contacts) {
      for (i=0; i<mol->nres; i++) free(mol->contacts[i]);
      free(mol->contacts);
    }
    if (mol->cutoffs) {
      for (i=0; i<mol->nres; i++) free(mol->cutoffs[i]);
      free(mol->cutoffs);
    }
    free(mol);
    mol=NULL;
}


atom_type* get_last_atom(atom_type* atom)
{
    while (atom->next) atom=atom->next;

  return atom;
}

res_type* get_last_res(res_type* res)
{
    while (res->next) res=res->next;

  return res;
}

mol_type *get_last_mol(mol_type* mol)
{
    while (mol->next) mol=mol->next;

  return mol;
}

// single-aa from 3-letter code
char setseq(char* aaname)
{
  int i;

    for (i=0; i<21; i++) {
      if ((aaname[0]==AA_NAMES[i][0]) &&
          (aaname[1]==AA_NAMES[i][1]) &&
          (aaname[2]==AA_NAMES[i][2]))
         break;
    }
    
    if (i==21) {
      if (!strcmp(aaname, "GLX"))
        return 'E';
      if (!strcmp(aaname, "ASX"))
        return 'D';
      if (!strcmp(aaname, "HID"))
        return 'H';
      if (!strcmp(aaname, "MSE"))
        return 'M';
      if (!strcmp(aaname, "SEP"))
        return 'S';
      if (!strcmp(aaname, "TPO"))
        return 'T';
      if (!strcmp(aaname, "PTR"))
        return 'Y';
      i--;
    }

  return SHORT_AA_NAMES[i];
}

// side chain - side chain orientation
int orient(res_type *res1, res_type *res2)
{
  real x1, y1, z1;
  real x2, y2, z2;
  real cax, cay, caz;
  real len, vect, angle;
  atom_type *atom;

    if (!res1 || !res2) return 0;

    atom=res1->atoms;
    cax=cay=caz=0.;
    while (atom) {
      if (!strncmp(atom->name,"CA",2)) {
        cax=atom->x; cay=atom->y; caz=atom->z;
      }
      atom=atom->next;
    }
    x1=res1->sgx-cax; y1=res1->sgy-cay; z1=res1->sgz-caz;
    if (x1==0. && y1==0. && z1==0.) x1+=1.0;

    atom=res2->atoms;
    cax=cay=caz=0.;
    while (atom) {
      if (!strncmp(atom->name,"CA",2)) {
        cax=atom->x; cay=atom->y; caz=atom->z;
      }
      atom=atom->next;
    }
    x2=res2->sgx-cax; y2=res2->sgy-cay; z2=res2->sgz-caz;
    if (x2==0. && y2==0. && z2==0.) x2+=1.0;

    vect = x1*x2+y1*y2+z1*z2;
    len = sqrt(x1*x1+y1*y1+z1*z1)*sqrt(x2*x2+y2*y2+z2*z2);
    if (len) vect /= len;

    angle=RADDEG*acos(vect);

    if (angle>120.) return 1; /*anti*/
    if (angle>60.) return 2;  /*mid*/

  return 3; /*par*/
}

int res_contact(res_type *res1, res_type *res2) {
  atom_type *atoms1, *atoms2;
  real dx, dy, dz;

    atoms1 = res1->atoms;
    while (atoms1) {
      atoms2 = res2->atoms;
      while (atoms2) {
        dx=atoms1->x-atoms2->x;
        dy=atoms1->y-atoms2->y;
        dz=atoms1->z-atoms2->z;
        if ((atoms1->flag & FLAG_SIDECHAIN) && (atoms2->flag & FLAG_SIDECHAIN) && (dx*dx+dy*dy+dz*dz<20.25)) {
          return 1;
        }
        atoms2=atoms2->next;
      }
      atoms1=atoms1->next;
    }

  return 0;
}


int read_pdb_file(char* filename, mol_type* molecules, char *realname)
{
  FILE *inp;
  char buffer[1000];
  char atmname[10];
  char resname[10];
  char version;
  int prevresnum, resnum, atmnum, locatmnum, num, locnum=0, i, j;
  atom_type *prevatom1, *prevatom2, *prevatom3, *prevatom4;
  int sgnum, cc, nres, ok, natom;
  real sgx, sgy, sgz;
  res_type *res, *test1, *test2;
  atom_type *atom;
  real x, y, z;
  real dist;
  unsigned char bin;
  int warn=0;
  real cutoff;

    if (_VERBOSE) printf("Reading input file %s...\n", filename);

    inp = fopen(filename, "r");
    if (!inp) {
      if (_VERBOSE) printf("ERROR: can't open %s !!!\n", filename);
      return FILE_NOT_FOUND;
    }

    molecules->nres=0;
    molecules->name=(char*)calloc(strlen(realname)+1,1);
    strcpy(molecules->name, realname);

    atmname[3]=0;
    resname[3]=0;
    prevresnum=-666;
    locatmnum=0;
    sgnum=0;
    sgx=sgy=sgz=0.;
    res=NULL;
    while (!feof(inp)) {
      if (fgets(buffer, 1000, inp)!=buffer) break;
      if (!strncmp(buffer, "END", 3) || !strncmp(buffer, "TER", 3)) break; // end of file; only single molecule read
      if (!strncmp(buffer, "ATOM", 4) || !strncmp(buffer, "HETATM", 6)) {
        if (buffer[16]!=' ' && buffer[16]!='A') continue;
        sscanf(&buffer[22], "%d", &resnum);
        strncpy(resname, &buffer[17], 3);
        strncpy(atmname, &buffer[13], 3);
        if (resnum==prevresnum && !strncmp(atmname, "N ", 2)) {
          if (_VERBOSE) printf("WARNING: fault in numeration at residuum %s[%d]\n", resname, resnum);
          warn=1;
        }
        if (atmname[0]=='H') continue;
        if (resnum!=prevresnum || !strncmp(atmname, "N ", 2)) {
          prevresnum=resnum;
          if (res)
            if (sgnum) {
              res->sgx=sgx/(real)sgnum;
              res->sgy=sgy/(real)sgnum;
              res->sgz=sgz/(real)sgnum;
            } else {
              res->sgx=res->sgy=res->sgz=0.;
            }
          locatmnum=0;
          version=' ';
          res = new_res();
          sgnum=0;
          sgx=sgy=sgz=0.;
          molecules->nres++;
          res->name = (char *)calloc(strlen(resname)+1, 1);  (char *)
          res->type = AA_NUMS[setseq(resname)];
          res->locnum=locnum++;
          res->num = resnum;
          res->natoms=0;
          res->chain = buffer[21];
          strcpy(res->name, resname);
          if (molecules->residua) {
            add_res(get_last_res(molecules->residua), res);
          } else {
            molecules->residua = res;
          }
        }
        atom = new_atom();
        atom->res = res;
        atom->flag |= FLAG_INITIAL;
        res->natoms++;
        locatmnum++;
        sscanf(&buffer[7], "%d", &atmnum);
        sscanf(&buffer[30], "%lf", &x);
        sscanf(&buffer[38], "%lf", &y);
        sscanf(&buffer[46], "%lf", &z);
        version = buffer[16];
        atom->name = (char *) calloc(strlen(atmname)+1,1);  (char *)
        strcpy(atom->name, atmname);
        atom->x=x; atom->y=y; atom->z=z;
        atom->num = atmnum;
        atom->locnum = locatmnum;
        if ((atmname[0]=='S' && atmname[1]=='C')||(atmname[0]=='C' && atmname[1]=='M')) {
          res->cmx = x;
          res->cmy = y;
          res->cmz = z;
          res->pdbsg=1;
          if (res->type<20) {
            res->protein=1;
          }
        } else
        if (!( ((atmname[0]=='C' || atmname[0]=='N' || atmname[0]=='O') && atmname[1]==' ') ||
               (atmname[0]=='H') ||
               (atmname[0]=='C' && atmname[1]=='A') ||
               (atmname[0]=='O' && atmname[1]=='X' && atmname[2]=='T') ) ) {
          sgx+=x;
          sgy+=y;
          sgz+=z;
          sgnum++;
          atom->flag |= FLAG_SIDECHAIN;
        } else
          atom->flag |= FLAG_BACKBONE;
        if (atmname[0]=='C' && atmname[1]=='A') {
          atom->flag |= FLAG_BACKBONE;
          if (res->type<20) {
            res->protein=1;
          }
          if (!res->pdbsg) {
            res->cmx = x;
            res->cmy = y;
            res->cmz = z;
          }
        }
        if (atmname[0]=='C' && atmname[1]=='M') {
          atom->flag |= FLAG_SCM;
        }
        if (atmname[0]=='S' && atmname[1]=='C') {
          atom->flag |= FLAG_SCM;
        }
        if (res->atoms) {
          add_atom(get_last_atom(res->atoms), atom);
        } else {
          res->atoms = atom;
        }
      }
    }

    if (res)
      if (sgnum) {
        res->sgx=sgx/(real)sgnum;
        res->sgy=sgy/(real)sgnum;
        res->sgz=sgz/(real)sgnum;
      } else {
        res->sgx=res->sgy=res->sgz=0.;
      }

    fclose(inp);

    molecules->seq = (uchar*)calloc(sizeof(uchar)*molecules->nres+1,1);
    res=molecules->residua; i=0;
    while (res) {
      molecules->seq[i++]=(uchar)AA_NUMS[(int)setseq(res->name)];
      res=res->next;
    }

  if (!warn) return FILE_SUCCESS; else return FILE_WARNING;
}


int read_pdb_file_saxs(char* filename, mol_type* molecules, char *realname, string pdbString)
{
  //FILE *inp;
  char buffer[1000];
  char atmname[10];
  char resname[10];
  char version;
  int prevresnum, resnum, atmnum, locatmnum, num, locnum=0, i, j;
  atom_type *prevatom1, *prevatom2, *prevatom3, *prevatom4;
  int sgnum, cc, nres, ok, natom;
  real sgx, sgy, sgz;
  res_type *res, *test1, *test2;
  atom_type *atom;
  real x, y, z;
  real dist;
  unsigned char bin;
  int warn=0;
  real cutoff;

    if (_VERBOSE) printf("Reading input file %s...\n", filename);

    //inp = fopen(filename, "r");
    //if (!inp) {
    //  if (_VERBOSE) printf("ERROR: can't open %s !!!\n", filename);
    //  return FILE_NOT_FOUND;
    //}

  
    molecules->nres=0;
    molecules->name=(char*)calloc(strlen(realname)+1,1);
    strcpy(molecules->name, realname);

    atmname[3]=0;
    resname[3]=0;
    prevresnum=-666;
    locatmnum=0;
    sgnum=0;
    sgx=sgy=sgz=0.;
    res=NULL;
    //while (!feof(inp)) {
    //  if (fgets(buffer, 1000, inp)!=buffer) break;
	std::string line;
	std::istringstream split(pdbString);
	while(std::getline(split, line, ';')) {
		//std::cout << "pulchar inside: " << line << endl;;
		strcpy(buffer, line.c_str());
	//}
	  
      if (!strncmp(buffer, "END", 3) || !strncmp(buffer, "TER", 3)) break; // end of file; only single molecule read
      if (!strncmp(buffer, "ATOM", 4) || !strncmp(buffer, "HETATM", 6)) {
		  //cout << "start read: "<< buffer << endl;
        if (buffer[16]!=' ' && buffer[16]!='A') continue;
        sscanf(&buffer[22], "%d", &resnum);
        strncpy(resname, &buffer[17], 3);
        strncpy(atmname, &buffer[13], 3);
        if (resnum==prevresnum && !strncmp(atmname, "N ", 2)) {
          if (_VERBOSE) printf("WARNING: fault in numeration at residuum %s[%d]\n", resname, resnum);
          warn=1;
        }
        if (atmname[0]=='H') continue;
        if (resnum!=prevresnum || !strncmp(atmname, "N ", 2)) {
          prevresnum=resnum;
          if (res)
            if (sgnum) {
              res->sgx=sgx/(real)sgnum;
              res->sgy=sgy/(real)sgnum;
              res->sgz=sgz/(real)sgnum;
            } else {
              res->sgx=res->sgy=res->sgz=0.;
            }
          locatmnum=0;
          version=' ';
          res = new_res();
          sgnum=0;
          sgx=sgy=sgz=0.;
          molecules->nres++;
          res->name = (char *)calloc(strlen(resname)+1, 1);  (char *)
          res->type = AA_NUMS[setseq(resname)];
          res->locnum=locnum++;
          res->num = resnum;
          res->natoms=0;
          res->chain = buffer[21];
          strcpy(res->name, resname);
          if (molecules->residua) {
            add_res(get_last_res(molecules->residua), res);
          } else {
            molecules->residua = res;
          }
        }
        atom = new_atom();
        atom->res = res;
        atom->flag |= FLAG_INITIAL;
        res->natoms++;
        locatmnum++;
        sscanf(&buffer[7], "%d", &atmnum);
        sscanf(&buffer[30], "%lf", &x);
        sscanf(&buffer[38], "%lf", &y);
        sscanf(&buffer[46], "%lf", &z);
        version = buffer[16];
        atom->name = (char *) calloc(strlen(atmname)+1,1);  (char *)
        strcpy(atom->name, atmname);
        atom->x=x; atom->y=y; atom->z=z;
        atom->num = atmnum;
        atom->locnum = locatmnum;
        if ((atmname[0]=='S' && atmname[1]=='C')||(atmname[0]=='C' && atmname[1]=='M')) {
          res->cmx = x;
          res->cmy = y;
          res->cmz = z;
          res->pdbsg=1;
          if (res->type<20) {
            res->protein=1;
          }
        } else
        if (!( ((atmname[0]=='C' || atmname[0]=='N' || atmname[0]=='O') && atmname[1]==' ') ||
               (atmname[0]=='H') ||
               (atmname[0]=='C' && atmname[1]=='A') ||
               (atmname[0]=='O' && atmname[1]=='X' && atmname[2]=='T') ) ) {
          sgx+=x;
          sgy+=y;
          sgz+=z;
          sgnum++;
          atom->flag |= FLAG_SIDECHAIN;
        } else
          atom->flag |= FLAG_BACKBONE;
        if (atmname[0]=='C' && atmname[1]=='A') {
          atom->flag |= FLAG_BACKBONE;
          if (res->type<20) {
            res->protein=1;
          }
          if (!res->pdbsg) {
            res->cmx = x;
            res->cmy = y;
            res->cmz = z;
          }
        }
        if (atmname[0]=='C' && atmname[1]=='M') {
          atom->flag |= FLAG_SCM;
        }
        if (atmname[0]=='S' && atmname[1]=='C') {
          atom->flag |= FLAG_SCM;
        }
        if (res->atoms) {
          add_atom(get_last_atom(res->atoms), atom);
        } else {
          res->atoms = atom;
        }
      }
    }
	
	//cout << "pulchar inside1" << endl;

    if (res)
      if (sgnum) {
        res->sgx=sgx/(real)sgnum;
        res->sgy=sgy/(real)sgnum;
        res->sgz=sgz/(real)sgnum;
      } else {
        res->sgx=res->sgy=res->sgz=0.;
      }

    //fclose(inp);
	//cout << "pulchar inside2" << endl;

    molecules->seq = (uchar*)calloc(sizeof(uchar)*molecules->nres+1,1);
    res=molecules->residua; i=0;
    while (res) {
      molecules->seq[i++]=(uchar)AA_NUMS[(int)setseq(res->name)];
      res=res->next;
    }
	//cout << "pulchar inside3" << endl;

  if (!warn) return FILE_SUCCESS; else return FILE_WARNING;
}

// energy calculation for C-alpha optimizer
real calc_ca_energy(atom_type **c_alpha, real **new_c_alpha, real **init_c_alpha, real **gradient, real alpha, real *ene, bool calc_gradient)
{
  int i, j;
  real dx, dy, dz;
  real dist, ddist, ddist2;
  real new_e_pot;
  real theta0, tdif, th, aa, bb, ab;
  real ff0, ff2, dth, m0, m2, grad, f0[3], f2[3];
  real adiff[3], bdiff[3];
  real deriv, theta, dtheta, len1, len2, cos_theta, sin_theta;
  real dx1, dy1, dz1;
  real dx2, dy2, dz2;
  real dx3, dy3, dz3;
  real vx1, vy1, vz1;
  real vx2, vy2, vz2;
  real vx3, vy3, vz3;

  real r12x, r12y, r12z;
  real r32x, r32y, r32z;
  real d12, d32, d12inv, d32inv, c1, c2, diff;
  real f1x, f1y, f1z;
  real f2x, f2y, f2z;
  real f3x, f3y, f3z;

        for (i=0; i<chain_length; i++) {
          new_c_alpha[i][0]=c_alpha[i]->x+alpha*gradient[i][0];
          new_c_alpha[i][1]=c_alpha[i]->y+alpha*gradient[i][1];
          new_c_alpha[i][2]=c_alpha[i]->z+alpha*gradient[i][2];
        }

        new_e_pot = 0.0;

        ene[0]=ene[1]=ene[2]=ene[3]=0.0;

        for (i=0; i<chain_length; i++) {
#ifdef CALC_C_ALPHA_START
          dx=new_c_alpha[i][0]-init_c_alpha[i][0];
          dy=new_c_alpha[i][1]-init_c_alpha[i][1];
          dz=new_c_alpha[i][2]-init_c_alpha[i][2];
          dist=sqrt(dx*dx+dy*dy+dz*dz);
          ddist = -dist;
          if (dist>_CA_START_DIST) {
            ddist2=dist*dist;
            new_e_pot+=CA_START_K*ddist2;
            ene[1] += CA_START_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_START_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
            }
          }

#endif


#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=new_c_alpha[i][0]-new_c_alpha[i-1][0];
            dy=new_c_alpha[i][1]-new_c_alpha[i-1][1];
            dz=new_c_alpha[i][2]-new_c_alpha[i-1][2];
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;
//              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;
//              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }
            ddist2=ddist*ddist;
            new_e_pot+=CA_K*ddist2;
            ene[0] += CA_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
              gradient[i-1][0]+=grad*dx;
              gradient[i-1][1]+=grad*dy;
              gradient[i-1][2]+=grad*dz;
            }
          }
#endif

#ifdef CALC_C_ALPHA_XVOL
          for (j=0;j<i;j++) {
            if (abs(i-j)>2) {
              dx=new_c_alpha[i][0]-new_c_alpha[j][0];
              dy=new_c_alpha[i][1]-new_c_alpha[j][1];
              dz=new_c_alpha[i][2]-new_c_alpha[j][2];
              dist=sqrt(dx*dx+dy*dy+dz*dz);
              ddist = dist-_CA_XVOL_DIST;
              if (dist<_CA_XVOL_DIST) {
                ddist2 = dist*dist;
                new_e_pot+=CA_XVOL_K*ddist2;
                ene[3] += CA_XVOL_K*ddist2;
                if (calc_gradient) {
                  grad = ddist*(8.0*CA_XVOL_K)/dist;
                  gradient[i][0]-=grad*dx;
                  gradient[i][1]-=grad*dy;
                  gradient[i][2]-=grad*dz;
                  gradient[j][0]+=grad*dx;
                  gradient[j][1]+=grad*dy;
                  gradient[j][2]+=grad*dz;
                }
              }
            }
          }
#endif

#ifdef CALC_C_ALPHA_ANGLES

        if (i>0 && i<chain_length-1) {
          r12x=new_c_alpha[i-1][0]-new_c_alpha[i][0];
          r12y=new_c_alpha[i-1][1]-new_c_alpha[i][1];
          r12z=new_c_alpha[i-1][2]-new_c_alpha[i][2];
          r32x=new_c_alpha[i+1][0]-new_c_alpha[i][0];
          r32y=new_c_alpha[i+1][1]-new_c_alpha[i][1];
          r32z=new_c_alpha[i+1][2]-new_c_alpha[i][2];
          d12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);
          d32 = sqrt(r32x*r32x+r32y*r32y+r32z*r32z);
          cos_theta = (r12x*r32x+r12y*r32y+r12z*r32z)/(d12*d32);
          if (cos_theta>1.0)
            cos_theta = 1.0;
          else
          if (cos_theta<-1.0)
            cos_theta = -1.0;
          sin_theta = sqrt(1.0-cos_theta*cos_theta);
          theta = acos(cos_theta);

          if (RADDEG*theta<80.)
            diff = theta-80.*DEGRAD;
          else
          if (RADDEG*theta>150.)
            diff = theta-150.*DEGRAD;
          else
            diff = 0.0;

          new_e_pot += CA_ANGLE_K*diff*diff;
          ene[2] += CA_ANGLE_K*diff*diff;
          if (calc_gradient) {
            d12inv = 1.0/d12;
            d32inv = 1.0/d32;
            diff *= (-2.0*CA_ANGLE_K)/sin_theta;
            c1 = diff*d12inv;
            c2 = diff*d32inv;
            f1x = c1*(r12x*(d12inv*cos_theta)-r32x*d32inv);
            f1y = c1*(r12y*(d12inv*cos_theta)-r32y*d32inv);
            f1z = c1*(r12z*(d12inv*cos_theta)-r32z*d32inv);
            f3x = c2*(r32x*(d32inv*cos_theta)-r12x*d12inv);
            f3y = c2*(r32y*(d32inv*cos_theta)-r12y*d12inv);
            f3z = c2*(r32z*(d32inv*cos_theta)-r12z*d12inv);
            f2x = -f1x-f3x;
            f2y = -f1y-f3y;
            f2z = -f1z-f3z;
            gradient[i-1][0]+=f1x;
            gradient[i-1][1]+=f1y;
            gradient[i-1][2]+=f1z;
            gradient[i][0]+=f2x;
            gradient[i][1]+=f2y;
            gradient[i][2]+=f2z;
            gradient[i+1][0]+=f3x;
            gradient[i+1][1]+=f3y;
            gradient[i+1][2]+=f3z;
          }
        }

#endif

        }

//printf("ene[3] = %f\n", ene[3]);

  return new_e_pot;
}

/*
 *  Steepest gradient optimization using v=k*(r0-r)^2
 *  k = CA_K, r0 = CA_DIST
 */
 

void ca_optimize(char *tname, char *iname)
{
  char buf[1000];
  int i, j, hx, my_iter;
  real dx, dy, dz, dd, dist, dist2, dist3, ddist, ddist2;
  real e_pot, new_e_pot, grad, alpha, e_pot1, e_pot2, e_pot3;
  real adiff[3], bdiff[3];
  real ff0, ff2, aa, ab, bb, th, tdif, dth, m0, m2;
  real theta0, deg_th, maxgrad, sum;
  real f0[3], f2[3];
  real x, y, z;
  int numsteps, numsteps2, msteps;
  int *sec;
  real **new_c_alpha, **gradient, **init_c_alpha, last_alpha, tmp, last_good_alpha, d_alpha, last_e_pot;
  atom_type *atom, **c_alpha;
  res_type *res;
  FILE *inp, *out;
  int mnum, init, ok;
  real alpha1, alpha2, alpha3, a0;
  real ene1, ene2, ene3, e0;
  real energies[4];
  real w1, w2, w3, eps;
  real gnorm, last_gnorm;
  int mode, fcnt;


    if (_CA_TRAJECTORY) {
      out = fopen(tname,"w");
      if (out) fclose(out);
    }

    if (_VERBOSE) printf("Alpha carbons optimization...\n");
      
    new_c_alpha = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    init_c_alpha = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    for (i=0;i<=chain_length;i++) {
      new_c_alpha[i] = (real*)calloc(sizeof(real)*3,1);
      init_c_alpha[i] = (real*)calloc(sizeof(real)*3,1);
    }
    gradient = (real**)calloc(sizeof(real*)*(chain_length+1),1);
    for (i=0;i<=chain_length;i++) {
      gradient[i] = (real*)calloc(sizeof(real)*3,1);
    }

    c_alpha = (atom_type**)calloc(sizeof(atom_type*)*(chain_length+1),1);

    i = 0;
    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A') {
          if (i<chain_length) {
            c_alpha[i] = atom;
            i++;
            break;
          } else {
            if (_VERBOSE) printf("WARNING: number of C-alpha atoms exceeds the chain length!\n");
            break;
          }
        }
        atom = atom->next;
      }
      res = res->next;
    }

    if (i<chain_length) chain_length = i;
      
    for (i=0; i<chain_length; i++) {
      init_c_alpha[i][0] = c_alpha[i]->x;
      init_c_alpha[i][1] = c_alpha[i]->y;
      init_c_alpha[i][2] = c_alpha[i]->z;
    }

    if (_CISPRO) {
      for (i=1; i<chain_length; i++) {
        dx = c_alpha[i]->x-c_alpha[i-1]->x;
        dy = c_alpha[i]->y-c_alpha[i-1]->y;
        dz = c_alpha[i]->z-c_alpha[i-1]->z;
        dd = sqrt(dx*dx+dy*dy+dz*dz);
        if ((setseq(c_alpha[i]->res->name)=='P') && (dd>CA_DIST_CISPRO-5*CA_DIST_CISPRO_TOL) && (dd<CA_DIST_CISPRO+5*CA_DIST_CISPRO_TOL)) {
          if (_VERBOSE) printf("Probable cis-proline found at postion %d\n", c_alpha[i]->res->num);
          c_alpha[i]->cispro = 1;
        }
      }
    }

    if (_CA_RANDOM) {
      if (_VERBOSE) printf("Generating random C-alpha coordinates...\n");
      c_alpha[0]->x = 0.0;
      c_alpha[0]->y = 0.0;
      c_alpha[0]->z = 0.0;
      for (i=1;i<chain_length;i++) {
        dx = 0.01*(100-rand()%200);
        dy = 0.01*(100-rand()%200);
        dz = 0.01*(100-rand()%200);
        dd = 3.8/sqrt(dx*dx+dy*dy+dz*dz);
        dx *= dd;
        dy *= dd;
        dz *= dd;
        c_alpha[i]->x = c_alpha[i-1]->x+dx;
        c_alpha[i]->y = c_alpha[i-1]->y+dy;
        c_alpha[i]->z = c_alpha[i-1]->z+dz;
      }
    }

    if (iname) {
      inp = fopen(iname,"r");
      if (inp) {
        if (_VERBOSE) printf("Reading initial structure %s...\n", iname);
        i = 0;
        while (!feof(inp)) {
          if (fgets(buf,1000,inp)==buf && buf[13]=='C' && buf[14]=='A') {
            if (i<chain_length) {
              if (sscanf(&buf[30],"%lf%lf%lf",&x,&y,&z)==3) {
                c_alpha[i]->x = x;
                c_alpha[i]->y = y;
                c_alpha[i]->z = z;
                i++;
              }
            } else {
              if (_VERBOSE) printf("WARNING: number of ini-file C-alpha atoms exceeds the chain length!\n");
              break;
            }
          }
        }
        fclose(inp);
      } else
        if (_VERBOSE) printf("WARNING: can't read initial corrdinates %s\n", iname);
    }

    mnum = 1;
    mode = 0;
    init = 0;
    numsteps=numsteps2=0;
    last_alpha = 0.0;


    if (_VERBOSE) printf("Optimizing alpha carbons...\n");

    eps = 0.5;

    fcnt=0;

    last_gnorm = 1000.;

    do {
      last_e_pot = e_pot;

      if (_CA_TRAJECTORY) {
        out = fopen(tname,"a");
        if (out) {
          fprintf(out,"MODEL  %d\n",mnum++);
          for (i=0; i<chain_length; i++) {
            fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                    i+1, "CA ", c_alpha[i]->res->name, ' ', c_alpha[i]->res->num,
  				          c_alpha[i]->x, c_alpha[i]->y, c_alpha[i]->z);

          }
          fprintf(out,"ENDMDL\n");
          fclose(out);
        }
      }

// calculate gradients

      e_pot=e_pot1=e_pot2=e_pot3=0.;

      for (i=0; i<chain_length; i++)
        gradient[i][0]=gradient[i][1]=gradient[i][2]=0.;

      e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, true);

      if (_VERBOSE && !init) {
        printf("Initial energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", energies[0], energies[2], energies[1], energies[3], e_pot);
      }

      if (!init) init=1;

// LINE SEARCH

      alpha1 = -1.0;
      alpha2 = 0.0;
      alpha3 = 1.0;

      ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
      ene2 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha2, energies, false);
      ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);

      msteps = 0;
      while (ene2>MIN(ene1,ene3) && msteps<_CA_ITER) {
        msteps++;
        alpha1 *= 2.0;
        alpha3 *= 2.0;
        ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
        ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);
      }

      msteps = 0;
      do {
        if (alpha3-alpha2>alpha2-alpha1) {
          a0 = 0.5*(alpha2+alpha3);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha1 = alpha2;
            alpha2 = a0;
            ene1 = ene2;
            ene2 = e0;
          } else {
            alpha3 = a0;
            ene3 = e0;
          }
        } else {
          a0 = 0.5*(alpha1+alpha2);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          e0 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha3 = alpha2;
            alpha2 = a0;
            ene3 = ene2;
            ene2 = e0;
          } else {
            alpha1 = a0;
            ene1 = e0;
          }
        }
        msteps++;
      } while (alpha3-alpha1>1e-6 && msteps<20);

      last_alpha = alpha2;
      e_pot = ene2;

      for (i=0; i<chain_length; i++) {
        c_alpha[i]->x=c_alpha[i]->x+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][0];
        c_alpha[i]->y=c_alpha[i]->y+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][1];
        c_alpha[i]->z=c_alpha[i]->z+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][2];
      }

      e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, false);

      eps *= 0.75;
      if (eps<1e-3) eps=0.0;

      numsteps++;

      gnorm = 0.0;

	  for (i=0; i<chain_length; i++) {
        gnorm += gradient[i][0]*gradient[i][0] + gradient[i][1]*gradient[i][1] + gradient[i][2]*gradient[i][2];
      }

      gnorm = sqrt(gnorm/(double)chain_length);

      if (last_gnorm-gnorm<1e-3) fcnt++;

      last_gnorm = gnorm;

    } while ( (fcnt<3) &&  (gnorm>0.01) && (numsteps<_CA_ITER));


     if (_VERBOSE) {
        for (i=0; i<chain_length; i++) {

#ifdef CALC_C_ALPHA
          if (i>0) {
            dx=c_alpha[i]->x-c_alpha[i-1]->x;
            dy=c_alpha[i]->y-c_alpha[i-1]->y;
            dz=c_alpha[i]->z-c_alpha[i-1]->z;
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            if (c_alpha[i]->cispro) {
              ddist=CA_DIST_CISPRO-dist;
              if (fabs(ddist)<CA_DIST_CISPRO_TOL) ddist=0.0;
            } else {
              ddist=CA_DIST-dist;
              if (fabs(ddist)<CA_DIST_TOL) ddist=0.0;
            }
            ddist2=ddist*ddist;
       	    if (fabs(ddist)>=CA_DIST_TOL) printf("WARNING: distance %d = %.3lf A\n", i, dist);
          }
#endif
        }

        for (i=0; i<chain_length; i++) {
#ifdef CALC_C_ALPHA_ANGLES
          if (i>0 && i<chain_length-1) {
            aa=ab=bb=0.0;
            adiff[0]=c_alpha[i-1]->x-c_alpha[i]->x;
            bdiff[0]=c_alpha[i+1]->x-c_alpha[i]->x;
            aa+=adiff[0]*adiff[0];
            ab+=adiff[0]*bdiff[0];
            bb+=bdiff[0]*bdiff[0];
            adiff[1]=c_alpha[i-1]->y-c_alpha[i]->y;
            bdiff[1]=c_alpha[i+1]->y-c_alpha[i]->y;
            aa+=adiff[1]*adiff[1];
            ab+=adiff[1]*bdiff[1];
            bb+=bdiff[1]*bdiff[1];
            adiff[2]=c_alpha[i-1]->z-c_alpha[i]->z;
            bdiff[2]=c_alpha[i+1]->z-c_alpha[i]->z;
            aa+=adiff[2]*adiff[2];
            ab+=adiff[2]*bdiff[2];
            bb+=bdiff[2]*bdiff[2];

            th=ab/sqrt(aa*bb);
            if (th<-1.0) th=-1.0;
            if (th>1.0) th=1.0;
            th=acos(th);
            deg_th=RADDEG*th;
            if (deg_th>150.) theta0=DEGRAD*150.; else
            if (deg_th<75.) theta0=DEGRAD*75.; else
            theta0=th;
       	    if (fabs(deg_th-RADDEG*theta0)>1.0) printf("WARNING: angle %d = %.3lf degrees\n", i, deg_th);
          }
#endif
        }
      }

    if (_VERBOSE) printf("Optimization done after %d step(s).\nFinal energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", numsteps, energies[0], energies[2], energies[1], energies[3], e_pot);

    if (_CA_TRAJECTORY) {
      out = fopen(tname,"a");
      if (out) {
        fprintf(out,"END\n");
      }
    }

    for (i=0;i<chain_length+1;i++) {
      free(init_c_alpha[i]);
      free(new_c_alpha[i]);
      free(gradient[i]);
    }
    free(new_c_alpha);
    free(gradient);
    free(c_alpha);
    free(init_c_alpha);
}



void center_chain(mol_type *mol)
{
  real cx, cy, cz;
  int natom;
  res_type *res;
  atom_type *atom;

    cx = cy = cz = 0.0;
    natom = 0;

    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        cx += atom->x;
        cy += atom->y;
        cz += atom->z;
        natom++;
  			atom=atom->next;
  		}
      res = res->next;
    }

    cx /= (real)natom;
    cy /= (real)natom;
    cz /= (real)natom;

    if (_VERBOSE) printf("Molecule center: %8.3f %8.3f %8.3f -> 0.000 0.000 0.000\n", cx, cy, cz);

    res = mol->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (!(_PRESERVE && (atom->flag & FLAG_INITIAL))) {        
          atom->x -= cx;
          atom->y -= cy;
          atom->z -= cz;
        }
        natom++;
  		atom=atom->next;
  	  }
      res = res->next;
    }

}



void write_pdb(char *name, mol_type *mol)
{
  FILE *out;
  res_type *res;
  atom_type *atom, *oxt;
  int anum;

    oxt = NULL;
    out = fopen(name,"w");
    if (!out) {
      if (_VERBOSE) printf("Can't open output file!\n");
      return;
    }
    fprintf(out,"REMARK 999 REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    anum=1;
    res = mol->residua;
    while (res) {
      if (res->protein) {
        if (!_BB_REARRANGE) {
          atom = res->atoms;
          while (atom) {
            if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
                !(atom->name[0]=='S' && atom->name[1]=='C') &&
                !(atom->name[0]=='C' && atom->name[1]=='M') &&
                !(atom->name[0]=='O' && atom->name[1]=='X') &&
                !(atom->name[0]=='H' && !_REBUILD_H))
              fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                            anum++, atom->name, res->name, ' ', res->num,
    	    				          atom->x, atom->y, atom->z);
            if (atom->name[0]=='O' && atom->name[1]=='X') oxt = atom;      	  
      	    atom=atom->next;
	  }
      	} else {
          atom = res->atoms;
          while (atom) {
            if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
                !(atom->name[0]=='S' && atom->name[1]=='C') &&
                !(atom->name[0]=='C' && atom->name[1]==' ') &&
                !(atom->name[0]=='O' && atom->name[1]==' ') &&
                !(atom->name[0]=='C' && atom->name[1]=='M') &&
                !(atom->name[0]=='O' && atom->name[1]=='X') &&
                !(atom->name[0]=='H' && !_REBUILD_H))
              fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                            anum++, atom->name, res->name, ' ', res->num,
    	    				          atom->x, atom->y, atom->z);
            if (atom->name[0]=='O' && atom->name[1]=='X') oxt = atom;      	  
      	    atom=atom->next;
      	  }
          atom = res->atoms;
          while (atom) {
            if (((atom->name[0]=='C' && atom->name[1]==' ') ||
                (atom->name[0]=='O' && atom->name[1]==' ')) &&
               !(atom->name[0]=='H' && !_REBUILD_H))
              fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                            anum++, atom->name, res->name, ' ', res->num,
    	    				          atom->x, atom->y, atom->z);
      			atom=atom->next;
      		}
      	}
      }
      if (!res->next && oxt) {
        atom = oxt;
        fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                     anum++, atom->name, res->name, ' ', res->num, atom->x, atom->y, atom->z);
      }
      res = res->next;
    }
    fprintf(out,"TER\nEND\n");
    fclose(out);
}


void write_pdb_insideSAXS(char *name, string &pdbString_pulchar, mol_type *mol)
{
  //FILE *out;
  res_type *res;
  atom_type *atom, *oxt;
  int anum;
  int anum2;
    oxt = NULL;
    //out = fopen(name,"w");
    //if (!out) {
    //  if (_VERBOSE) printf("Can't open output file!\n");
    //  return;
    //}
	char content[1000];
    //fprintf(out,"REMARK 999 REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    sprintf(content,"REMARK 999 REBUILT BY PULCHRA V.%.2f;", PULCHRA_VERSION);
	//cout << "In wirtepdb: " << content << endl;
	pdbString_pulchar.append(content);
	
    anum=1;
    anum2=1;
    res = mol->residua;
    while (res) {
      if (res->protein) {
        if (!_BB_REARRANGE) {
          atom = res->atoms;
          while (atom) {
            if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
                !(atom->name[0]=='S' && atom->name[1]=='C') &&
                !(atom->name[0]=='C' && atom->name[1]=='M') &&
                !(atom->name[0]=='O' && atom->name[1]=='X') &&
                !(atom->name[0]=='H' && !_REBUILD_H))
				{
				  //fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
					//			anum++, atom->name, res->name, ' ', res->num,
					//					  atom->x, atom->y, atom->z);
				  char content[1000];
				  sprintf(content, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f;",
							   anum2++, atom->name, res->name, ' ', res->num,
										  atom->x, atom->y, atom->z);
				  pdbString_pulchar.append(content);
				  
				  //cout << "In wirtepdb res: " << content << endl;
				}
			  
            if (atom->name[0]=='O' && atom->name[1]=='X') oxt = atom;      	  
      	    atom=atom->next;
	  }
      	} else {
          atom = res->atoms;
          while (atom) {
            if (!(atom->name[0]=='D' && atom->name[1]=='U') &&
                !(atom->name[0]=='S' && atom->name[1]=='C') &&
                !(atom->name[0]=='C' && atom->name[1]==' ') &&
                !(atom->name[0]=='O' && atom->name[1]==' ') &&
                !(atom->name[0]=='C' && atom->name[1]=='M') &&
                !(atom->name[0]=='O' && atom->name[1]=='X') &&
                !(atom->name[0]=='H' && !_REBUILD_H))
				{
					//fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
						//		anum++, atom->name, res->name, ' ', res->num,
						//				  atom->x, atom->y, atom->z);
					char content[1000];
					sprintf(content, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f;",
								anum2++, atom->name, res->name, ' ', res->num,
										  atom->x, atom->y, atom->z);
					pdbString_pulchar.append(content);
					//cout << "In wirtepdb unres: " << content << endl;
			   }
			  
			  
            if (atom->name[0]=='O' && atom->name[1]=='X') oxt = atom;      	  
      	    atom=atom->next;
      	  }
          atom = res->atoms;
          while (atom) {
            if (((atom->name[0]=='C' && atom->name[1]==' ') ||
                (atom->name[0]=='O' && atom->name[1]==' ')) &&
               !(atom->name[0]=='H' && !_REBUILD_H))
			   {
				  //fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
						//		anum++, atom->name, res->name, ' ', res->num,
						//				  atom->x, atom->y, atom->z);
				  char content[1000];
				  sprintf(content, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f;",
								anum2++, atom->name, res->name, ' ', res->num,
										  atom->x, atom->y, atom->z);
				  pdbString_pulchar.append(content);
				  //cout << "In wirtepdb unbb: " << content << endl;
			   }
			  
      			atom=atom->next;
      		}
      	}
      }
      if (!res->next && oxt) {
        atom = oxt;
        //fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                  //   anum++, atom->name, res->name, ' ', res->num, atom->x, atom->y, atom->z);
 
        char content[1000];
        sprintf(content, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f;",
                     anum2++, atom->name, res->name, ' ', res->num, atom->x, atom->y, atom->z);
        pdbString_pulchar.append(content);
		//cout << "In wirtepdb ok: " << content << endl;
		
      }
      res = res->next;
    }
    //fprintf(out,"TER\nEND\n");
	char content2[1000];
    sprintf(content2,"TER;END;");
	pdbString_pulchar.append(content2);
    //fclose(out);
}

void write_pdb_sg(char *name, mol_type *mol)
{
  FILE *out;
  res_type *res;
  atom_type *atom;
  int anum;

    out = fopen(name,"w");
    if (!out) {
      if (_VERBOSE) printf("Can't open output file!\n");
      return;
    }
    fprintf(out,"REMARK 999 REBUILT BY PULCHRA V.%.2f\n", PULCHRA_VERSION);
    anum=1;
    res = mol->residua;
    while (res) {
      if (res->protein) {
        atom = res->atoms;
        while (atom) {
          if ((atom->name[0]=='C' && atom->name[1]=='A'))
            fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                          anum++, atom->name, res->name, ' ', res->num,
    	  				          atom->x, atom->y, atom->z);
    			atom=atom->next;
    		}
        fprintf(out, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
                      anum++, "CM ", res->name, ' ', res->num,
    				          res->cmx, res->cmy, res->cmz);
    	}
      res = res->next;
    }
    fprintf(out,"TER\nEND\n");
    fclose(out);
}


real calc_distance(real x1, real y1, real z1,real x2, real y2, real z2)
{
  real dx,dy,dz;
  real dist2;

    dx = (x1) - (x2);
    dy = (y1) - (y2);
    dz = (z1) - (z2);
    if (dx || dy || dz ) {
      dist2 = dx*dx+dy*dy+dz*dz;
      return (sqrt(dist2));
    } else
      return 0.0;
}

// r14 chiral distance

real calc_r14(real x1, real y1, real z1,
							 real x2, real y2, real z2,
							 real x3, real y3, real z3,
							 real x4, real y4, real z4)
{
  real r, dx, dy, dz;
  real vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3;
  real hand;

    dx = x4-x1;
    dy = y4-y1;
    dz = z4-z1;

    r = sqrt(dx*dx+dy*dy+dz*dz);

    vx1=x2-x1;
    vy1=y2-y1;
    vz1=z2-z1;
    vx2=x3-x2;
    vy2=y3-y2;
    vz2=z3-z2;
    vx3=x4-x3;
    vy3=y4-y3;
    vz3=z4-z3;

    hand = (vy1*vz2-vy2*vz1)*vx3+
           (vz1*vx2-vz2*vx1)*vy3+
           (vx1*vy2-vx2*vy1)*vz3;

    if (hand<0) r=-r;

  return r;
}

// superimposition of two sets for coordinates + optional transformation of tpoints

real superimpose2(real **coords1, real **coords2, int npoints, real **tpoints, int ntpoints)
{
  real mat_s[3][3], mat_a[3][3], mat_b[3][3], mat_g[3][3];
  real mat_u[3][3], tmp_mat[3][3];
  real val, d, alpha, beta, gamma, x, y, z;
  real cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz;
  int i, j, k, n;

    cx1=cy1=cz1=cx2=cy2=cz2=0.;

    for (i=0; i<npoints; i++) {
      cx1+=coords1[i][0];
      cy1+=coords1[i][1];
      cz1+=coords1[i][2];
      cx2+=coords2[i][0];
      cy2+=coords2[i][1];
      cz2+=coords2[i][2];
    }

    cx1/=(real)npoints;
    cy1/=(real)npoints;
    cz1/=(real)npoints;

    cx2/=(real)npoints;
    cy2/=(real)npoints;
    cz2/=(real)npoints;

    for (i=0; i<npoints; i++) {
      coords1[i][0]-=cx1;
      coords1[i][1]-=cy1;
      coords1[i][2]-=cz1;
      coords2[i][0]-=cx2;
      coords2[i][1]-=cy2;
      coords2[i][2]-=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]-=cx2;
      tpoints[i][1]-=cy2;
      tpoints[i][2]-=cz2;
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
        if (i==j)
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=1.0;
        else
          mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=0.0;
        mat_u[i][j]=0.;
      }

    for (n=0; n<npoints; n++) {
      mat_u[0][0]+=coords1[n][0]*coords2[n][0];
      mat_u[0][1]+=coords1[n][0]*coords2[n][1];
      mat_u[0][2]+=coords1[n][0]*coords2[n][2];
      mat_u[1][0]+=coords1[n][1]*coords2[n][0];
      mat_u[1][1]+=coords1[n][1]*coords2[n][1];
      mat_u[1][2]+=coords1[n][1]*coords2[n][2];
      mat_u[2][0]+=coords1[n][2]*coords2[n][0];
      mat_u[2][1]+=coords1[n][2]*coords2[n][1];
      mat_u[2][2]+=coords1[n][2]*coords2[n][2];
    }

    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        tmp_mat[i][j]=0.;

    do {
      d=mat_u[2][1]-mat_u[1][2];
      if (d==0) alpha=0; else alpha=atan(d/(mat_u[1][1]+mat_u[2][2]));
      if (cos(alpha)*(mat_u[1][1]+mat_u[2][2])+sin(alpha)*(mat_u[2][1]-mat_u[1][2])<0.0)       alpha+=M_PI;
      mat_a[1][1]=mat_a[2][2]=cos(alpha);
      mat_a[2][1]=sin(alpha);
      mat_a[1][2]=-mat_a[2][1];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_a[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_a[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[0][2]-mat_u[2][0];
      if (d==0) beta=0; else beta=atan(d/(mat_u[0][0]+mat_u[2][2]));
      if (cos(beta)*(mat_u[0][0]+mat_u[2][2])+sin(beta)*(mat_u[0][2]-mat_u[2][0])<0.0) beta+=M_PI;
      mat_b[0][0]=mat_b[2][2]=cos(beta);
      mat_b[0][2]=sin(beta);
      mat_b[2][0]=-mat_b[0][2];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_b[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_b[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[1][0]-mat_u[0][1];
      if (d==0) gamma=0; else gamma=atan(d/(mat_u[0][0]+mat_u[1][1]));
      if (cos(gamma)*(mat_u[0][0]+mat_u[1][1])+sin(gamma)*(mat_u[1][0]-mat_u[0][1])<0.0)
        gamma+=M_PI;
      mat_g[0][0]=mat_g[1][1]=cos(gamma);
      mat_g[1][0]=sin(gamma);
      mat_g[0][1]=-mat_g[1][0];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_g[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_g[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      val=fabs(alpha)+fabs(beta)+fabs(gamma);
    } while (val>0.001);

    val=0.;
    for (i=0; i<npoints; i++) {
      x=coords2[i][0];
      y=coords2[i][1];
      z=coords2[i][2];
      tmpx=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tmpy=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tmpz=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
      x=coords1[i][0]-tmpx;
      y=coords1[i][1]-tmpy;
      z=coords1[i][2]-tmpz;
      val+=x*x+y*y+z*z;
    }

    for (i=0; i<ntpoints; i++) {
      x=tpoints[i][0];
      y=tpoints[i][1];
      z=tpoints[i][2];
      tpoints[i][0]=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tpoints[i][1]=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tpoints[i][2]=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
    }

    for (i=0; i<npoints; i++) {
      coords1[i][0]+=cx1;
      coords1[i][1]+=cy1;
      coords1[i][2]+=cz1;
      coords2[i][0]+=cx2;
      coords2[i][1]+=cy2;
      coords2[i][2]+=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]+=cx1;
      tpoints[i][1]+=cy1;
      tpoints[i][2]+=cz1;
    }

  return sqrt(val/(real)npoints);
}


atom_type *find_atom(res_type *res, char *aname)
{
  atom_type *atom;

    atom = res->atoms;
    while (atom) {
      if (atom->name[0]==aname[0] && atom->name[1]==aname[1] && atom->name[2]==aname[2]) {
        return atom;
        break;
      }
      atom = atom->next;
    }

  return NULL;
}

void add_replace(res_type *res, char *aname, real x, real y, real z, int flags)
{
  atom_type *atom, *newatom;

    atom = res->atoms;
    while (atom) {
      if (atom->name[0]==aname[0] && atom->name[1]==aname[1] && atom->name[2]==aname[2]) {
        if (!(_PRESERVE && (atom->flag & FLAG_INITIAL))) { atom->x = x; atom->y = y; atom->z = z; }
        atom->flag |= flags;
        break;
      }
      atom = atom->next;
    }

    if (!atom) {
      newatom = (atom_type*)calloc(sizeof(atom_type),1);
      newatom->x = x;
      newatom->y = y;
      newatom->z = z;
      newatom->flag |= flags;
      newatom->res = res;
      newatom->name = (char*)calloc(4,1);
      strcpy(newatom->name,aname);

      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A')
          break;
        atom = atom->next;
      }
      if (aname[0]=='N' && aname[1]==' ') {
        newatom->next = res->atoms;
        res->atoms = newatom;
      } else {
        while (atom->next) atom=atom->next;
        atom->next = newatom;
      }
    }
}




void prepare_rbins(void)
{
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real r13_1, r13_2, r14;
  real **cacoords, **tmpcoords, **tmpstat;
  res_type *res, *prevres;
  atom_type *atom;

  if (!RBINS) {
    RBINS = (int**)calloc(sizeof(int*)*(chain_length+1),1);
    for (i=0;i<chain_length+1;i++)
      RBINS[i] = (int*)calloc(sizeof(int)*3,1);

    X_COORDS = (real**)calloc(sizeof(real*)*(chain_length+10),1);
    for (i=0;i<chain_length+10;i++)
      X_COORDS[i] = (real*)calloc(sizeof(real)*3,1);

    i = 5;

    res = chain->residua;
    while (res) {
      atom = res->atoms;
      while (atom) {
        if (atom->name[0]=='C' && atom->name[1]=='A') {
          X_COORDS[i][0] = atom->x;
          X_COORDS[i][1] = atom->y;
          X_COORDS[i][2] = atom->z;
          i++;
        }
        atom = atom->next;
      }
      res = res->next;
    }

    C_ALPHA = &X_COORDS[5];

    cacoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpcoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpstat = (real**)calloc(sizeof(real*)*(8),1);
    for (i=0;i<8;i++) {
      cacoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpcoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpstat[i] = (real*)calloc(sizeof(real)*3,1);;
    }

    // rebuild ends...

    for (i=0,j=0;i<5;i++,j++)
      for (k=0;k<3;k++)
        tmpcoords[j][k] = C_ALPHA[i][k];
    for (i=2,j=0;i<5;i++,j++)
      for (k=0;k<3;k++)
        cacoords[j][k] = C_ALPHA[i][k];
    for (i=0,j=0;i<3;i++,j++)
      for (k=0;k<3;k++)
        tmpstat[j][k] = C_ALPHA[i][k];

    superimpose2(tmpstat,cacoords,3,tmpcoords,5);

    for (i=-2,j=0;i<0;i++,j++)
      for (k=0;k<3;k++)
        C_ALPHA[i][k] = tmpcoords[j][k];

    for (i=chain_length-5,j=0;i<chain_length;i++,j++)
      for (k=0;k<3;k++)
        tmpcoords[j][k] = C_ALPHA[i][k];
    for (i=chain_length-5,j=0;i<chain_length-2;i++,j++)
      for (k=0;k<3;k++)
        cacoords[j][k] = C_ALPHA[i][k];
    for (i=chain_length-3,j=0;i<chain_length;i++,j++)
      for (k=0;k<3;k++)
        tmpstat[j][k] = C_ALPHA[i][k];

    superimpose2(tmpstat,cacoords,3,tmpcoords,5);

    for (i=chain_length-3,j=0;i<chain_length;i++,j++)
      for (k=0;k<3;k++)
        C_ALPHA[i+3][k] = tmpcoords[j+3][k];

    for (i=0;i<chain_length+1;i++) {
    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];

    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];

    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];

    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

    	r13_1 = calc_distance(x1, y1, z1, x3, y3, z3);
    	r13_2 = calc_distance(x2, y2, z2, x4, y4, z4);
    	r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

      bin13_1 = (int)((r13_1-4.6)/0.3);
      bin13_2 = (int)((r13_2-4.6)/0.3);
      bin14 = (int)((r14+11.)/0.3);

      if (bin13_1<0) bin13_1=0;
      if (bin13_2<0) bin13_2=0;
      if (bin14<0) bin14=0;
      if (bin13_1>9) bin13_1=9;
      if (bin13_2>9) bin13_2=9;
      if (bin14>73) bin14=73;

      RBINS[i][0] = bin13_1;
      RBINS[i][1] = bin13_2;
      RBINS[i][2] = bin14;
    }
  }
}

#ifdef COMPILE_BB

void rebuild_backbone(void)
{

  res_type *res, *prevres;
  atom_type *atom;
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real besthit, hit;
  int bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real rmsd, total, maxrms;
  FILE *debug, *out;

    if (_VERBOSE) printf("Rebuilding backbone...\n");

    prepare_rbins();

    cacoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpcoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpstat = (real**)calloc(sizeof(real*)*(8),1);
    for (i=0;i<8;i++) {
      cacoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpcoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpstat[i] = (real*)calloc(sizeof(real)*3,1);;
    }


    prevres = NULL;
    res = chain->residua;


    total = maxrms = 0.0;

    for (i=0;i<chain_length+1;i++) {
    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];

    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];

    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];

    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

    	cacoords[0][0] = x1;
    	cacoords[0][1] = y1;
    	cacoords[0][2] = z1;

    	cacoords[1][0] = x2;
     	cacoords[1][1] = y2;
     	cacoords[1][2] = z2;

     	cacoords[2][0] = x3;
     	cacoords[2][1] = y3;
     	cacoords[2][2] = z3;

     	cacoords[3][0] = x4;
     	cacoords[3][1] = y4;
     	cacoords[3][2] = z4;

      bin13_1 = RBINS[i][0];
      bin13_2 = RBINS[i][1];
      bin14 = RBINS[i][2];

      pro = 0;

      if (prevres && !strncmp(prevres->name,"PRO",3)) {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  }
     	  }
      } else {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);
        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
        for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat[bestpos].data[j][k];
       	  }
     	  }
      }

     	rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 8);

     	total += rmsd;
     	if (rmsd>maxrms) maxrms=rmsd;

// add-or-replace

      if (prevres) {
        add_replace(prevres, "C  ", tmpcoords[4][0], tmpcoords[4][1], tmpcoords[4][2], FLAG_BACKBONE);
        add_replace(prevres, "O  ", tmpcoords[5][0], tmpcoords[5][1], tmpcoords[5][2], FLAG_BACKBONE);
      }

      if (res) {
        add_replace(res, "N  ", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      } else { // terminal oxygen instead of nitrogen
        add_replace(prevres, "OXT", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      }

      prevres = res;
      if (res)
        res = res->next;
    }

    if (_VERBOSE) printf("Backbone rebuilding deviation: average = %.3f, max = %.3f\n", total/(real)chain_length, maxrms);

}

#endif


#ifdef COMPILE_ROT

typedef struct _rot_struct {
  int r13_1, r13_2, r14;
  int nc;
  real ***coords;
  struct _rot_struct *next;
} rot_struct;

rot_struct *rotamers[20];

/* this is obsolete in a standalone version of PULCHRA */

int read_rotamers(void)
{
  FILE *inp;
  char buf[1000];
  char dum[100];
  int aa, i, j, k, l, n;
  rot_struct *new_rot, *last_rot;
  real x, y, z;

    if (_VERBOSE) printf("Reading rotamer library...\n");

    inp = fopen("NEWROT","r");
    last_rot=NULL;
    while (!feof(inp)) {
      if (fgets(buf,1000,inp)==buf) {
        if (buf[0]=='A') {
          sscanf(buf,"%s %d", dum, &aa);
          if (last_rot) last_rot->next = NULL;
          last_rot = NULL;
          if (fgets(buf,1000,inp)!=buf) break;
        }
//        printf("aa: %d\n", aa);
        if (aa==20) break;
        sscanf(buf,"%d %d %d %s %d", &i, &j, &k, dum, &l);
        new_rot = (rot_struct*)calloc(sizeof(rot_struct),1);
//        printf("%d %d %d nc: %d\n", i, j, k, l);
        new_rot->r13_1 = i;
        new_rot->r13_2 = j;
        new_rot->r14 = k;
        new_rot->nc = l;
        new_rot->next = NULL;
        new_rot->coords = (real***)calloc(sizeof(real**)*l,1);
        for (i=0;i<l;i++) {
          new_rot->coords[i]=(real**)calloc(sizeof(real*)*(nheavy[aa]+1),1);
          for (j=0;j<(nheavy[aa]+1);j++) {
            new_rot->coords[i][j]=(real*)calloc(sizeof(real)*3,1);
          }
        }
        for (i=0;i<l;i++) {
          fgets(buf,1000,inp);
          for (j=0;j<(nheavy[aa]+1);j++) {
            fgets(buf,1000,inp);
            sscanf(buf,"%lf%lf%lf",&x, &y, &z);
            new_rot->coords[i][j][0]=x;
            new_rot->coords[i][j][1]=y;
            new_rot->coords[i][j][2]=z;
          }
          if (last_rot) {
            last_rot->next = new_rot;
          } else {
            rotamers[aa] = new_rot;
          }
          last_rot = new_rot;
        }
      }
    }
    fclose(inp);
}


void cross(real *v1, real *v2, real *v3)
{
  v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}


void norm(real *v)
{
  real d;

    d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}


int check_xvol(res_type *res)
{
  res_type *res2;
  atom_type *atom1, *atom2;
  real dx, dy, dz, dd;

    res2 = chain->residua;

    while (res2) {
      atom2 = res2->atoms;
      if (res!=res2) {
        while (atom2) {
          atom1 = res->atoms;
          while (atom1) {
            if (atom1->flag & FLAG_SIDECHAIN) {
              dx = atom1->x-atom2->x;
              dy = atom1->y-atom2->y;
              dz = atom1->z-atom2->z;
              dd = dx*dx+dy*dy+dz*dz;
              if (dd<(1.7*1.7)) {
                return 1;
              }
            }
            atom1=atom1->next;
          }
          atom2=atom2->next;
        }
      }
      res2=res2->next;
    }

  return 0;
}


real ***SORTED_ROTAMERS;


void rebuild_sidechains(void)
{
  FILE *out;
  res_type *res, *prevres, *testres;
  atom_type *atom, *atom1, *atom2;
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real x5, y5, z5;
  real r14, r13_1, r13_2;
  real dx, dy, dz, dd;
  real hit, besthit;
  int exvol, bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14;
  real rmsd, total;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  real **lsys, **vv, **sc;
  char scn[12][4];
  rot_struct *rot;
  int ok, last_a, last_b, last_c, last_d, jpos;
  int jx, jy, jz, jxi, jyi, jzi, b13_1, b13_2, b14, jm;
  int crot, bestrot, minexvol, totexvol, rtried, pos, cpos;
  real cmx, cmy, cmz, ddx, ddy, ddz, ddd, bestdd;
  real sort_rot[100][2];

    if (_VERBOSE) printf("Rebuilding side chains...\n");

    prepare_rbins();

    lsys = (real**)calloc(sizeof(real*)*3,1);
    vv = (real**)calloc(sizeof(real*)*3,1);
    sc = (real**)calloc(sizeof(real*)*12,1);
    for (i=0;i<12;i++)
      sc[i] = (real*)calloc(sizeof(real)*3,1);
    for (i=0;i<3;i++) {
      lsys[i] = (real*)calloc(sizeof(real)*3,1);
      vv[i] = (real*)calloc(sizeof(real)*3,1);
    }

    SORTED_ROTAMERS = (real***)calloc(sizeof(real**)*(chain_length+1),1);
    for (i=0;i<chain_length+1;i++) {
      SORTED_ROTAMERS[i] = (real**)calloc(sizeof(real*)*10,1);
      for (j=0;j<10;j++) {
        SORTED_ROTAMERS[i][j] = (real*)calloc(sizeof(real)*2,1);
      }
    }

    prevres = NULL;
    res = chain->residua;
    totexvol = 0;

    for (i=0;i<chain_length;i++) {
      if (!strncmp(res->name,"GLY",3) || !res->protein) {
        if (res->next) res = res->next;
        continue;
      }

  		x1 = C_ALPHA[i-2][0];
  		y1 = C_ALPHA[i-2][1];
  		z1 = C_ALPHA[i-2][2];
  		x2 = C_ALPHA[i-1][0];
  		y2 = C_ALPHA[i-1][1];
  		z2 = C_ALPHA[i-1][2];
  		x3 = C_ALPHA[i][0];
  		y3 = C_ALPHA[i][1];
  		z3 = C_ALPHA[i][2];
  		x4 = C_ALPHA[i+1][0];
  		y4 = C_ALPHA[i+1][1];
  		z4 = C_ALPHA[i+1][2];

      bin13_1 = RBINS[i][0];
      bin13_2 = RBINS[i][1];
      bin14 = RBINS[i][2];

      v1[0] = x4-x2;
      v1[1] = y4-y2;
      v1[2] = z4-z2;

      v2a[0] = x4-x3;
      v2a[1] = y4-y3;
      v2a[2] = z4-z3;

      v2b[0] = x3-x2;
      v2b[1] = y3-y2;
      v2b[2] = z3-z2;

      cross(v2a, v2b, v2);
      cross(v1, v2, v3);

      norm(v1);
      norm(v2);
      norm(v3);

// gather 10 closest rotamer conformations...

      for (j=0;j<10;j++)
        SORTED_ROTAMERS[i][j][0] = 500.;

      j = 0;
      besthit = 1000.;
      bestpos = 0;
      do {
        if (rot_stat_idx[j][0]==res->type) {
          hit = fabs(rot_stat_idx[j][1]-bin13_1)+fabs(rot_stat_idx[j][2]-bin13_2)+0.2*fabs(rot_stat_idx[j][3]-bin14);
          if (hit<SORTED_ROTAMERS[i][9][0]) {
            k = 9;
            while (k>=0 && hit<SORTED_ROTAMERS[i][k][0]) {
              k--;
            }
            k++;
            // k = hit
            for (l=9;l>k;l--) {
              SORTED_ROTAMERS[i][l][0]=SORTED_ROTAMERS[i][l-1][0];
              SORTED_ROTAMERS[i][l][1]=SORTED_ROTAMERS[i][l-1][1];
            }
            SORTED_ROTAMERS[i][k][0]=hit;
            SORTED_ROTAMERS[i][k][1]=j;
          }
        }
        j++;
      } while (rot_stat_idx[j][0]>=0);


      besthit = SORTED_ROTAMERS[i][0][0];
      bestpos = SORTED_ROTAMERS[i][0][1];


// new rebuild...

      pos = rot_stat_idx[bestpos][5];
      nsc = nheavy[res->type]+1;

      if (_PDB_SG) { // more than one rotamer - check SC
        bestdd = 100.; crot = 0;
        for (l=0;l<2;l++) { // check two closest conformations
          cpos = SORTED_ROTAMERS[i][l][1];
          for (m=0;m<rot_stat_idx[cpos][4];m++) {
            for (j=0;j<3;j++) {
              vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
              for (k=0;k<3;k++) {
                if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
              }
            }
            pos = rot_stat_idx[cpos][5]+nsc*m;
            for (j=0;j<nsc;j++) {
              for (k=0;k<3;k++) {
                sc[j][k] = rot_stat_coords[pos+j][k];
              }
            }
            superimpose2(vv,lsys,3,sc,nsc);
            for (j=0;j<nsc;j++) {
              sc[j][0] += x3;
              sc[j][1] += y3;
              sc[j][2] += z3;
            }
            cmx = 0.; cmy = 0.; cmz = 0.;
            for (j=0;j<nsc;j++) {
              cmx += sc[j][0];
              cmy += sc[j][1];
              cmz += sc[j][2];
            }
            cmx /= (real) nsc;
            cmy /= (real) nsc;
            cmz /= (real) nsc;
            ddx = res->cmx-cmx;
            ddy = res->cmy-cmy;
            ddz = res->cmz-cmz;
            ddx *= ddx;
            ddy *= ddy;
            ddz *= ddz;
            ddd = ddx+ddy+ddz;
            if (ddd<bestdd) {
              bestdd = ddd;
              crot = pos; // closest rotamer position
            }
          }
        }
        pos = crot;
      } // PDB_SG

      for (j=0;j<3;j++) {
        vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
        for (k=0;k<3;k++) {
          if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
        }
      }

      for (j=0;j<nsc;j++) {
        for (k=0;k<3;k++) {
          sc[j][k] = rot_stat_coords[pos+j][k];
        }
      }

      superimpose2(vv,lsys,3,sc,nsc);

      for (j=0;j<nsc;j++) {
        sc[j][0] += x3;
        sc[j][1] += y3;
        sc[j][2] += z3;
      }

      for (j=1;j<nsc;j++) {
        add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
      }

      if (res->next) res = res->next;

    } // i++, next res

    for (i=0;i<12;i++)
      free(sc[i]);
    for (i=0;i<3;i++) {
      free(lsys[i]);
      free(vv[i]);
    }
    free(sc);
	free(lsys);
	free(vv);
}



// finds steric conflicts

int get_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
  atom_list *llist;
  atom_type *atom, *atom2;
  int i, j, k, x, y, z;
  int ii, jj, kk, con, iter, maxcon, merged;
  real dx, dy, dz, dd;

    con = 0;
    atom = res->atoms;
    while (atom) {
      i = atom->gx;
      j = atom->gy;
      k = atom->gz;
      for (ii=i-2;ii<=i+2;ii++)
        for (jj=j-2;jj<=j+2;jj++)
          for (kk=k-2;kk<=k+2;kk++) {
            if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<zgrid) {
              llist = grid[ii][jj][kk];
              while (llist) {
                atom2 = llist->atom;
                if (atom && atom2 && res && atom2->res) {
                  merged=0;
                  if (res==atom2->res) { // self-xvol
                    if (atom->flag & FLAG_SIDECHAIN && atom2->flag & FLAG_SIDECHAIN) merged=1;
                    if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='A' && atom2->name[0]=='C' && atom2->name[1]=='B') merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='B' && atom2->name[0]=='C' && atom2->name[1]=='A') merged=1;
                    if (res->name[0]=='P') {
                      if (atom->name[0]=='C' && atom->name[1]=='D' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                      if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]=='D') merged=1;
                    }

                    if (!merged) {
//                      printf("merged: %s[%d] %s-%s %d %d\n", res->name,res->num,atom->name,atom2->name,atom->flag,atom2->flag);
                    }
                  } else
                  if (res->next==atom2->res || res==atom2->res->next) {
                    if (atom->name[0]=='C' && atom->name[1]==' ' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                    if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]==' ') merged=1;
                  }
                  if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1; // for now
                  if (atom->flag & FLAG_SCM || atom2->flag & FLAG_SCM) merged=1; // for now
                  if (!merged) {
                    dx = atom->x-atom2->x;
                    dx*=dx;
                    dy = atom->y-atom2->y;
                    dy*=dy;
                    dz = atom->z-atom2->z;
                    dz*=dz;
                    dd = dx+dy+dz;
                    if (dd<_SG_XVOL_DIST*_SG_XVOL_DIST) {
                      con++;
                    }
                  }
                }
                llist = llist->next;
              }
            }
          }
      atom = atom->next;
    }

  return con;
}

int display_conflicts(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
  atom_list *llist;
  atom_type *atom, *atom2;
  int i, j, k, x, y, z;
  int ii, jj, kk, con, iter, maxcon, merged;
  real dx, dy, dz, dd;

    con = 0;
    atom = res->atoms;
    while (atom) {
      i = atom->gx;
      j = atom->gy;
      k = atom->gz;
      for (ii=i-2;ii<=i+2;ii++)
        for (jj=j-2;jj<=j+2;jj++)
          for (kk=k-2;kk<=k+2;kk++) {
            if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<zgrid) {
              llist = grid[ii][jj][kk];
              while (llist) {
                atom2 = llist->atom;
                if (atom && atom2 && res && atom2->res) {
                  merged=0;
                  if (res==atom2->res) { // self-xvol
                    if (atom->flag & FLAG_SIDECHAIN && atom2->flag & FLAG_SIDECHAIN) merged=1;
                    if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='A' && atom2->name[0]=='C' && atom2->name[1]=='B') merged=1;
                    if (atom->name[0]=='C' && atom->name[1]=='B' && atom2->name[0]=='C' && atom2->name[1]=='A') merged=1;
                    if (res->name[0]=='P') {
                      if (atom->name[0]=='C' && atom->name[1]=='D' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                      if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]=='D') merged=1;
                    }
                    if (!merged) {
//                      printf("merged: %s[%d] %s-%s %d %d\n", res->name,res->num,atom->name,atom2->name,atom->flag,atom2->flag);
                    }
                  } else
                  if (res->next==atom2->res || res==atom2->res->next) {
                    if (atom->name[0]=='C' && atom->name[1]==' ' && atom2->name[0]=='N' && atom2->name[1]==' ') merged=1;
                    if (atom->name[0]=='N' && atom->name[1]==' ' && atom2->name[0]=='C' && atom2->name[1]==' ') merged=1;
                  }
                  if (atom->flag & FLAG_BACKBONE && atom2->flag & FLAG_BACKBONE) merged=1; // for now
                  if (atom->flag & FLAG_SCM || atom2->flag & FLAG_SCM) merged=1; // for now
                  if (!merged) {
                    dx = atom->x-atom2->x;
                    dx*=dx;
                    dy = atom->y-atom2->y;
                    dy*=dy;
                    dz = atom->z-atom2->z;
                    dz*=dz;
                    dd = dx+dy+dz;
                    if (dd<1.6*1.6) {
                      printf("STERIC CONFLICT: %s[%d]%s-%s[%d]%s\n", atom->res->name,atom->res->num,atom->name,atom2->res->name,atom2->res->num,atom2->name);
                      con++;
                    }
                  }
                }
                llist = llist->next;
              }
            }
          }
      atom = atom->next;
    }

  return con;
}


void allocate_grid(atom_list *****grid_, int *xgrid_, int *ygrid_, int *zgrid_)
{
  static int xgrid, ygrid, zgrid;
  static atom_list ****grid = NULL;
  atom_list *llist, *alist;
  real min[3], max[3];
  res_type *res, *worst;
  atom_type *atom, *atom2;
  int i, j, x, y, z;

    if (!grid && chain->residua && chain->residua->atoms) {
 	    res = chain->residua;
      min[0]=max[0]=res->atoms->x;
      min[1]=max[1]=res->atoms->y;
      min[2]=max[2]=res->atoms->z;
	    while (res) {
	      atom = res->atoms;
	      while (atom) {
	        if (atom->x<min[0]) min[0]=atom->x;
	        if (atom->y<min[1]) min[1]=atom->y;
	        if (atom->z<min[2]) min[2]=atom->z;
	        if (atom->x>max[0]) max[0]=atom->x;
	        if (atom->y>max[1]) max[1]=atom->y;
	        if (atom->z>max[2]) max[2]=atom->z;
	        atom = atom->next;
	      }
	      res = res->next;
	    }

	    xgrid = (max[0]-min[0])/GRID_RES;
	    ygrid = (max[1]-min[1])/GRID_RES;
	    zgrid = (max[2]-min[2])/GRID_RES;

	    if (_VERBOSE) printf("Allocating grid (%d %d %d)...\n", xgrid, ygrid, zgrid);

	   grid = (atom_list****)calloc(sizeof(atom_list***)*(xgrid+1),1);
	   for (i=0;i<xgrid+1;i++) {
	     grid[i] = (atom_list***)calloc(sizeof(atom_list**)*(ygrid+1),1);
	     for (j=0;j<ygrid+1;j++) {
	       grid[i][j] = (atom_list**)calloc(sizeof(atom_list*)*(zgrid+1),1);
	     }
	   }

	   res = chain->residua;
	   while (res) {
	     atom = res->atoms;
	     while (atom) {
	       x = xgrid*(atom->x-min[0])/(max[0]-min[0]);
	       y = ygrid*(atom->y-min[1])/(max[1]-min[1]);
	       z = zgrid*(atom->z-min[2])/(max[2]-min[2]);
	       alist = (atom_list*)calloc(sizeof(atom_list),1);
	       alist->atom = atom;
	       atom->gx = x;
	       atom->gy = y;
	       atom->gz = z;
	       if (grid[x][y][z]!=NULL) {
	         llist = grid[x][y][z];
	         while (llist->next) llist=llist->next;
	         llist->next = alist;
	       } else {
	         grid[x][y][z]=alist;
	       }
	       atom = atom->next;
	     }
	     res = res->next;
	   }
	} else {
	   if (_VERBOSE) printf("Grid already allocated (%d %d %d)\n", xgrid, ygrid, zgrid);
	}

	*grid_ = grid;
	*xgrid_ = xgrid;
	*ygrid_ = ygrid;
	*zgrid_ = zgrid;
}

void optimize_exvol(void)
{
  real min[3], max[3];
  res_type *res, *worst;
  atom_type *atom, *atom2;
  int xgrid, ygrid, zgrid;
  atom_list ****grid, *llist, *alist;
  int i, j, k, l, m, x, y, z;
  int ii, jj, kk, con, iter, maxcon, totcon;
  int cpos, bestpos, pos, con0;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  real **lsys, **vv, **sc;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;

    min[0]=1e5;
    min[1]=1e5;
    min[2]=1e5;
    max[0]=-1e5;
    max[1]=-1e5;
    max[2]=-1e5;

    lsys = (real**)calloc(sizeof(real*)*3,1);
    vv = (real**)calloc(sizeof(real*)*3,1);
    sc = (real**)calloc(sizeof(real*)*12,1);
    for (i=0;i<12;i++)
      sc[i] = (real*)calloc(sizeof(real)*3,1);
    for (i=0;i<3;i++) {
      lsys[i] = (real*)calloc(sizeof(real)*3,1);
      vv[i] = (real*)calloc(sizeof(real)*3,1);
    }

  	allocate_grid(&grid, &xgrid, &ygrid, &zgrid);

    if (_VERBOSE) printf("Finding excluded volume conflicts...\n", xgrid, ygrid, zgrid);

	iter = 0;

	do {
//printf("ITER: %d\n", iter);

      maxcon = 0;
      totcon=0;

      res = chain->residua;
      while (res) {
        if (res->protein) {
          con = get_conflicts(res, grid, xgrid, ygrid, zgrid);
          if (con>0) {
            totcon+=con;
            if (con>maxcon) {
              maxcon = con;
              worst = res;
            }
          }
        }
        res = res->next;
      }

      if (_VERBOSE && iter==0) {
        printf("Total number of conflicts: %d\n", totcon);
      }

      if (totcon==0) break;

      if (_VERBOSE && iter==0) {
        printf("Maximum number of conflicts: %s[%d] : %d\n", worst->name, worst->num, maxcon);
      }

      totcon=0;

      if (maxcon>0) {

// try to fix...

    res = chain->residua;
    for (i=0;i<chain_length;i++) {
      if (!strncmp(res->name,"GLY",3) || !res->protein) {
        if (res->next) res = res->next;
        continue;
      }

      nsc = nheavy[res->type]+1;

    	x1 = C_ALPHA[i-2][0];
    	y1 = C_ALPHA[i-2][1];
    	z1 = C_ALPHA[i-2][2];
    	x2 = C_ALPHA[i-1][0];
    	y2 = C_ALPHA[i-1][1];
    	z2 = C_ALPHA[i-1][2];
    	x3 = C_ALPHA[i][0];
    	y3 = C_ALPHA[i][1];
    	z3 = C_ALPHA[i][2];
    	x4 = C_ALPHA[i+1][0];
    	y4 = C_ALPHA[i+1][1];
    	z4 = C_ALPHA[i+1][2];

      v1[0] = x4-x2;
      v1[1] = y4-y2;
      v1[2] = z4-z2;

      v2a[0] = x4-x3;
      v2a[1] = y4-y3;
      v2a[2] = z4-z3;

      v2b[0] = x3-x2;
      v2b[1] = y3-y2;
      v2b[2] = z3-z2;

      cross(v2a, v2b, v2);
      cross(v1, v2, v3);

      norm(v1);
      norm(v2);
      norm(v3);

      con = get_conflicts(res, grid, xgrid, ygrid, zgrid);

      if (con>0) {

        bestpos=0;
        con0 = 100;
        for (l=0;l<10;l++) { // check two closest conformations
          cpos = SORTED_ROTAMERS[i][l][1];
          for (m=0;m<rot_stat_idx[cpos][4];m++) {
            for (j=0;j<3;j++) {
              vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
              for (k=0;k<3;k++) {
                if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
              }
            }
            pos = rot_stat_idx[cpos][5]+nsc*m;
            for (j=0;j<nsc;j++) {
              for (k=0;k<3;k++) {
                sc[j][k] = rot_stat_coords[pos+j][k];
              }
            }
            superimpose2(vv,lsys,3,sc,nsc);
            for (j=0;j<nsc;j++) {
              sc[j][0] += x3;
              sc[j][1] += y3;
              sc[j][2] += z3;
            }
            for (j=1;j<nsc;j++) {
              add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
            }
            con = get_conflicts(res, grid, xgrid, ygrid, zgrid);
//printf("test: %d\n", con);

            if (con<con0) {
              con0 = con;
              bestpos = pos;
            }
            if (con==0) break;
          }
          if (con==0) break;
        }

		totcon += con0;

        for (j=0;j<3;j++) {
          vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
          for (k=0;k<3;k++) {
            if (j==k) lsys[j][k]=1.; else lsys[j][k]=0.;
          }
        }
        pos = bestpos;
        for (j=0;j<nsc;j++) {
          for (k=0;k<3;k++) {
            sc[j][k] = rot_stat_coords[pos+j][k];
          }
        }
        superimpose2(vv,lsys,3,sc,nsc);
        for (j=0;j<nsc;j++) {
          sc[j][0] += x3;
          sc[j][1] += y3;
          sc[j][2] += z3;
        }
        for (j=1;j<nsc;j++) {
          add_replace(res, heavy_atoms[10*res->type+j-1], sc[j][0], sc[j][1], sc[j][2], FLAG_SIDECHAIN);
        }
      }



      	res=res->next;

    	} // i
	}

	iter++;

	} while (iter<_XVOL_ITER);


    if (_VERBOSE) {
      if (totcon>0)
        printf("WARNING: %d steric conflict(s) are still there.\n", totcon);
      else
        printf("All steric conflicts removed.\n");
    }

    for (i=0;i<12;i++)
      free(sc[i]);
    for (i=0;i<3;i++) {
      free(lsys[i]);
      free(vv[i]);
    }
    free(sc);
	free(lsys);
	free(vv);


}


void vcross(real ax,real ay,real az,real bx,real by,real bz,real *cx,real *cy,real *cz)
{
    *cx = ay * bz - by * az;
    *cy = az * bx - bz * ax;
    *cz = ax * by - bx * ay;
}

real vdot(real ax,real ay,real az,real bx,real by,real bz)
{
    return ax*bx+ay*by+az*bz;
}

real calc_torsion(atom_type *a1, atom_type *a2, atom_type *a3, atom_type *a4)
{
  real v12x, v12y, v12z;
  real v43x, v43y, v43z;
  real zx, zy, zz;
  real px, py, pz;
  real xx, xy, xz;
  real yx, yy, yz;
  real u, v, angle;

    v12x = a1->x-a2->x;
    v12y = a1->y-a2->y;
    v12z = a1->z-a2->z;

    v43x = a4->x-a3->x;
    v43y = a4->y-a3->y;
    v43z = a4->z-a3->z;

    zx = a2->x-a3->x;
    zy = a2->y-a3->y;
    zz = a2->z-a3->z;

    vcross(zx,zy,zz,v12x,v12y,v12z,&px,&py,&pz);
    vcross(zx,zy,zz,v43x,v43y,v43z,&xx,&xy,&xz);
    vcross(zx,zy,zz,xx,xy,xz,&yx,&yy,&yz);

    u = vdot(xx,xy,xz,xx,xy,xz);
    v = vdot(yx,yy,yz,yx,yy,yz);

    angle = 360.;

    if (u<0. || v<0.) return angle;

    u = vdot(px,py,pz,xx,xy,xz) / sqrt(u);
    v = vdot(px,py,pz,yx,yy,yz) / sqrt(v);

    if (u != 0.0 || v != 0.0) angle = atan2(v, u) * RADDEG;


  return angle;

}


// Ca-N-C-Cb angle should be close to 34 deg
// check and fix 

int chirality_check(void)
{
  int i;
  atom_type *a_ca, *a_n, *a_c, *a_cb;
  atom_type *atom;
  res_type *res;
  real angle;
  real nx, ny, nz;
  real px, py, pz;
  real qx, qy, qz;
  real rx, ry, rz;
  real xx, xy, xz;
  real yx, yy, yz;
  real dd, costheta, sintheta;

    if (_VERBOSE) printf("Checking chirality...\n");
    res = chain->residua;
    while (res) {
      a_ca = a_n = a_c = a_cb = NULL;
      a_ca = find_atom(res,"CA ");
      a_n = find_atom(res,"N  ");
      a_c = find_atom(res,"C  ");
      a_cb = find_atom(res,"CB ");
      if (a_ca && a_n && a_c && a_cb) {
        angle = calc_torsion(a_ca, a_n, a_c, a_cb);
        if (angle<0.) {
          if (_VERBOSE) printf("WARNING: D-aa detected at %s %3d : %5.2f", res->name, res->num, angle);
          xx = a_ca->x-a_n->x;
          xy = a_ca->y-a_n->y;
          xz = a_ca->z-a_n->z;
          yx = a_c->x-a_ca->x;
          yy = a_c->y-a_ca->y;
          yz = a_c->z-a_ca->z;
          vcross(xx,xy,xz,yx,yy,yz,&nx,&ny,&nz);
          dd = sqrt(nx*nx+ny*ny+nz*nz);
          nx /= dd;
          ny /= dd;
          nz /= dd;
          // nx, ny, nz = reflection plane normal
          rx = xx-yx;
          ry = xy-yy;
          rz = xz-yz;
          dd = sqrt(rx*rx+ry*ry+rz*rz);
          rx /= dd;
          ry /= dd;
          rz /= dd;
          costheta = -1.;
          sintheta = 0.;
          atom = res->atoms;
          while (atom) {
            if (atom->flag & FLAG_SIDECHAIN) {
              px = atom->x-a_ca->x;
              py = atom->y-a_ca->y;
              pz = atom->z-a_ca->z;
              qx = qy = qz = 0.;
              qx += (costheta + (1 - costheta) * rx * rx) * px;
              qx += ((1 - costheta) * rx * ry - rz * sintheta) * py;
              qx += ((1 - costheta) * rx * rz + ry * sintheta) * pz;
              qy += ((1 - costheta) * rx * ry + rz * sintheta) * px;
              qy += (costheta + (1 - costheta) * ry * ry) * py;
              qy += ((1 - costheta) * ry * rz - rx * sintheta) * pz;
              qz += ((1 - costheta) * rx * rz - ry * sintheta) * px;
              qz += ((1 - costheta) * ry * rz + rx * sintheta) * py;
              qz += (costheta + (1 - costheta) * rz * rz) * pz;
              qx += a_ca->x;
              qy += a_ca->y;
              qz += a_ca->z;
              atom->x = qx;
              atom->y = qy;
              atom->z = qz;
            }
            atom = atom->next;
          }
          angle = calc_torsion(a_ca, a_n, a_c, a_cb);
          if (_VERBOSE) printf(", fixed : %5.2f\n", angle);
        }
      }
      res = res->next;
    }
}


#endif

// DSSP energy of petide-peptide HB

real hb_energy(res_type *res, atom_list ****grid, int xgrid, int ygrid, int zgrid)
{
	atom_type *atom, *c_atom1, *o_atom1, *n_atom1, *c_atom2, *o_atom2, *n_atom2, *tmp_atom;
	atom_type h_atom;
	int i, j, k, ii, jj, kk;
  atom_list *llist, *alist;
  real dx, dy, dz, dist, min_dist1, min_dist2;
  real hx1, hy1, hz1, dd;
  real dno, dnc, dho, dhc;
  real ene, Q;

    ene = 1e3;

    if (!res || !res->prev) return ene;
            
    Q = -27888.0; // DSSP h-bond energy constant

		c_atom1 = o_atom1 = n_atom1 = NULL;

  	atom = res->prev->atoms;
  	while (atom) {
			if (atom->name[0]=='C' && atom->name[1]==' ') c_atom1 = atom;
			if (atom->name[0]=='O' && atom->name[1]==' ') o_atom1 = atom;
			atom = atom->next;
		}

  	atom = res->atoms;
  	while (atom) {
			if (atom->name[0]=='N' && atom->name[1]==' ') { n_atom1 = atom; break; }
			atom = atom->next;
		}

// first bond

    min_dist2 = 1e10;
    o_atom2 = c_atom2 = NULL;
		if (n_atom1) {
			i = n_atom1->gx;
			j = n_atom1->gy;
			k = n_atom1->gz;
			for (ii=i-1;ii<=i+1;ii++) {
				for (jj=j-1;jj<=j+1;jj++) {
					for (kk=k-1;kk<=k+1;kk++) {
						if (ii>=0 && ii<xgrid && jj>=0 && jj<ygrid && kk>=0 && kk<=zgrid) {
							llist = grid[ii][jj][kk];
							while (llist) {
							  if (llist->atom->name[0]=='O' && llist->atom->name[1]==' ' && abs(llist->atom->res->locnum-n_atom1->res->locnum)>2)  {
							    tmp_atom = llist->atom;
							    dx = n_atom1->x-tmp_atom->x;
							    dy = n_atom1->y-tmp_atom->y;
							    dz = n_atom1->z-tmp_atom->z;
							    dist = dx*dx+dy*dy+dz*dz;
							    if (dist<min_dist2 && dist<25.0) {
							      o_atom2=tmp_atom;
							      min_dist2 = dist;
   						    } 
							  } 
							  llist = llist->next;
		          }					
						}
					}
				}
			}
		}

    if (o_atom2) {
      atom = o_atom2->res->atoms;
    	while (atom) {
  			if (atom->name[0]=='C' && atom->name[1]==' ') { c_atom2 = atom; break; }
  			atom = atom->next;
  		}                      
      if (c_atom2) {    
    		hx1 = o_atom1->x-c_atom1->x;
    		hy1 = o_atom1->y-c_atom1->y;
    		hz1 = o_atom1->z-c_atom1->z;
    		dd = -1.081f/sqrt(hx1*hx1+hy1*hy1+hz1*hz1);
    		hx1 *= dd;
    		hy1 *= dd;
    		hz1 *= dd;
    		
    		hx1 += n_atom1->x;
    		hy1 += n_atom1->y;
    		hz1 += n_atom1->z;
        
        add_replace(n_atom1->res, "H  ", hx1, hy1, hz1, FLAG_BACKBONE);

  // dno
        dx = n_atom1->x-o_atom2->x;
        dy = n_atom1->y-o_atom2->y;
        dz = n_atom1->z-o_atom2->z;
        dno = sqrt(dx*dx+dy*dy+dz*dz);
  
  // dnc
        dx = n_atom1->x-c_atom2->x;
        dy = n_atom1->y-c_atom2->y;
        dz = n_atom1->z-c_atom2->z;
        dnc = sqrt(dx*dx+dy*dy+dz*dz);
  
  // dho
        dx = hx1-o_atom2->x;
        dy = hy1-o_atom2->y;
        dz = hz1-o_atom2->z;
        dho = sqrt(dx*dx+dy*dy+dz*dz);
  
  // dhc
        dx = hx1-c_atom2->x;
        dy = hy1-c_atom2->y;
        dz = hz1-c_atom2->z;
        dhc = sqrt(dx*dx+dy*dy+dz*dz);
        if (dho<0.01F || dhc<0.01F || dnc<0.01F || dno<0.01F) {
          ene = -10.0;
        } else {
          ene = 0.001*(Q/dho - Q/dhc + Q/dnc - Q/dno);
        }       
      }
    }

  return ene;
}

// rotates a point around a vector
void rot_point_vector(real *x, real *y, real *z, real u, real v, real w, real angle)
{
  real ux, uy, uz, vx, vy, vz, wx, wy, wz, sa, ca;
  
    sa = sinf(10.0*M_PI*angle/180.0);
    ca = cosf(10.0*M_PI*angle/180.0);
    
    ux = u**x;
    uy = u**y;
    uz = u**z;
    vx = v**x;
    vy = v**y;
    vz = v**z;
    wx = w**x;
    wy = w**y;
    wz = w**z;

    *x = u*(ux+vy+wz)+(*x*(v*v+w*w)-u*(vy+wz))*ca+(-wy+vz)*sa;
    *y = v*(ux+vy+wz)+(*y*(u*u+w*w)-v*(ux+wz))*ca+( wx-uz)*sa;
    *z = w*(ux+vy+wz)+(*z*(u*u+v*v)-w*(ux+vy))*ca+(-vx+uy)*sa;
}


// rotates a peptide plate

void rot_peptide(res_type *res, real angle)
{
	atom_type *atom, *c_atom, *o_atom, *n_atom, *ca_atom1, *ca_atom2;
  real u, v, w, x, y, z, dd;

    if (!res || !res->prev) return;
      
    c_atom = o_atom = n_atom = ca_atom1 = ca_atom2 = NULL;
    
  	atom = res->prev->atoms;
  	while (atom) {
			if (atom->name[0]=='C' && atom->name[1]=='A') ca_atom1 = atom;
			if (atom->name[0]=='C' && atom->name[1]==' ') c_atom = atom;
			if (atom->name[0]=='O' && atom->name[1]==' ') o_atom = atom;
			atom = atom->next;
		}
		
    atom = res->atoms;
  	while (atom) {
			if (atom->name[0]=='C' && atom->name[1]=='A') ca_atom2 = atom;
			if (atom->name[0]=='N' && atom->name[1]==' ') n_atom = atom;
			atom = atom->next;
		}
		
    if (c_atom && o_atom && n_atom && ca_atom1 && ca_atom2) {
      u = ca_atom2->x-ca_atom1->x;
      v = ca_atom2->y-ca_atom1->y;
      w = ca_atom2->z-ca_atom1->z;
      dd = 1.0f/sqrt(u*u+v*v+w*w); 
      u*=dd; v*=dd; w*=dd; // normalize ca-ca vector
      x = n_atom->x-ca_atom1->x;
      y = n_atom->y-ca_atom1->y;
      z = n_atom->z-ca_atom1->z;
      rot_point_vector(&x, &y, &z, u, v, w, angle);
      n_atom->x = x+ca_atom1->x;
      n_atom->y = y+ca_atom1->y;
      n_atom->z = z+ca_atom1->z;
      x = c_atom->x-ca_atom1->x;
      y = c_atom->y-ca_atom1->y;
      z = c_atom->z-ca_atom1->z;
      rot_point_vector(&x, &y, &z, u, v, w, angle);
      c_atom->x = x+ca_atom1->x;
      c_atom->y = y+ca_atom1->y;
      c_atom->z = z+ca_atom1->z;
      x = o_atom->x-ca_atom1->x;
      y = o_atom->y-ca_atom1->y;
      z = o_atom->z-ca_atom1->z;
      rot_point_vector(&x, &y, &z, u, v, w, angle);
      o_atom->x = x+ca_atom1->x;
      o_atom->y = y+ca_atom1->y;
      o_atom->z = z+ca_atom1->z;
    }  		
  
}

void optimize_backbone(mol_type *chain)
{
  int xgrid, ygrid, zgrid;
  atom_list ****grid;
  atom_type *atom;
	res_type *res;
  real ene, min_ene, tot1, tot2;
  int i, k, best;
FILE *out;
   
  	if (_VERBOSE) printf("Optimizing backbone...\n", xgrid, ygrid, zgrid);

		allocate_grid(&grid, &xgrid, &ygrid, &zgrid);

    tot1 = tot2 = 0.0;
    
		res = chain->residua;
		while (res) {
			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
			if (ene<-0.5) tot1 += ene;
			res = res->next;
		}

		res = chain->residua;
		while (res) {
		  if (res->type!=7) {
  			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			if (ene<1.0) { // try to optimize
  			  min_ene = ene;
  			  rot_peptide(res, -1.1);
  			  best = 0;
  			  for (i=-10;i<10;i++) {
  			    rot_peptide(res, 0.1);
      			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			    if (ene<min_ene) {
  			      best = i;
  			      min_ene = ene;
  			    }			  
  			  }	
  			  rot_peptide(res,-0.9);		
    			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			  if (min_ene<ene) {
  			    rot_peptide(res,0.1*best);
      			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
  			  } 
  			}			
  		}
			res = res->next;
		}

		res = chain->residua;
		while (res) {
			ene = hb_energy(res, grid, xgrid, ygrid, zgrid);
			if (ene<-0.5) tot2 += ene;
			res = res->next;
		}

    if (_VERBOSE) printf("Backbone HB energy: before %g, after: %g, difference: %g\n", tot1, tot2, tot2-tot1);

}