// analysis.h

#ifndef ANALYSIS_H
#define ANALYSIS_H

#define MAXCELL 9000000
#define rcutoff 1.122462048

#include "../../LD/LD-cpp/src/library.h"
#include "../../LD/LD-cpp/src/interactions.h"
#include "library.h"

using namespace program;

namespace analysis {

	class atomsXYZ {

	public:

	char element;
	int atom_id, type, si;
	float rxt1, ryt1;
	float rxt2, ryt2;		// *) Used in MSD calculations
	float vx, vy, vz;
	float vxth, vyth, vzth;
	float fx, fy;
	float fx_int, fy_int;
	int jumpx, jumpy;

	atomsXYZ();
	~atomsXYZ();
	
	};

	typedef class molecularXYZ : public atomsXYZ {

	public:

	int mol_id;
	float rzt1;
	int jumpz;

	molecularXYZ();
	~molecularXYZ();

	}atom_style;

	class System {

	public:

	int dimension;
	float Lx, Ly, Lz, rcellx, rcelly, rcellz;
	int Ncellx, Ncelly, Ncellz, ncells;
	int nAtoms, nAtomTypes, nMols;

	int MAPS[MAXCELL], HEAD[MAXCELL], LIST[MAXCELL];

	System(float Lx, float Ly, float Lz, int nAtoms, int nAtomTypes, bool buildMaps = false, float rcellx = rcutoff, float rcelly = rcutoff, float rcellz = rcutoff);
	~System();

	int cellindex(int ix, int iy);
	void buildCellMaps();
	void buildCellList(atom_style *ATOMS);
	void checkMinImage(float *dx = NULL, float *dy = NULL, float *dz = NULL);

	};

	class velocityDist {

	public:

	int nDims, nBins;
	float **Bins;
	float *min_val;
	float *max_val;
	float *binW;
	long *ctr;

	velocityDist(int nDims, int nBins, System *BOX);
	~velocityDist(){};

	void createBins(System *BOX);
	void scanVelocities(atom_style *ATOMS, System *BOX);
	void binVelocities(atom_style *ATOMS, System *BOX, bool subtract_vcom = false);
	void normalize();

	};

	class Bin1D {

	public:

	float Lx, Ly, binWidth;
	int nBins;
	float *bin;
	float *binavg;

	Bin1D(float Lx, float Ly, float binWidth);
	Bin1D();
	~Bin1D();

	void addToBin(float rx, float val = 1.0);
	void normalize(float *binTonorm, float *normv = NULL, float normc = 1.0);
	void zero(float *binTozero);
	void addBins(float *bin1, float *bin2);

	};

	class Hist2D {

		public:

		float rcut_x, rcut_y;
		float binW_x, binW_y;
		int nBin_x, nBin_y;
		float **bin;

		Hist2D(float rcut_x, float rcut_y, float binW_x, float binW_y);
		~Hist2D();

		void addToBin(char iid, char jid, float dxij, float dyij);
		void normalize(float fac = 1.0);
	};

	class Trajectory {

	public:

	long step; 
	int frame_nr, totalFrames, frameWidth;
	float timeStep, time;
	int nAtoms, nAtomTypes, nMols;
	int dimension;
	float Lx, Ly, Lz;
	float xCom, yCom, zCom; 
	char *format;
	char line_fmt[50];
	int line_fmt_mode;

	char *fpathI, *fpathO, *pipeString, *pipeChar;
	FILE *fileI, *fileO;
	fpos_t curr_pos;

	Trajectory(float timeStep = 1.0, int frameWidth = 1, char fileFormat[10] = "xyz");
	~Trajectory();

	void openTrajectory(bool count = false);
	void createOutputFile(char line[] = "", bool newfile = true);
	void closeTrajectory(bool closeI = true, bool closeO = true);
	void loadTrajectory(atom_style **ATOMS, System *BOX, int frameStart, int frameEnd, bool unwrap_pbc = false);
	void countFrames();
	void readThisFrame(atom_style *ATOMS);
	void readNextFrame(atom_style *ATOMS);
	void writeThisFrame(atom_style *ATOMS, System *BOX, long add_step = 0);
	void copyThisFrame(atom_style *From, atom_style *To, bool unwrap_pbc = false);
	void computeCom(atom_style *ATOMS);
	long **sortAtomsByType(atom_style *ATOMS, System *BOX);

	void write2file();
	void write2file(Bin1D *bin1, int ctr = 0);
	void write2file(Bin1D *binA, Bin1D *binB, int ctr = 0, float timeToAvg = 0.0);
	void write2file(Bin1D **Pin, Bin1D **Pkin, Bin1D *Pswim, int ctr = 0);
	void write2file(double **MSD, long *count, int *delFrames, int nSample);
	void write2file(float time, float order, int step = -1);
	void write2file(float *rn, float **RDF, int nRDF, int nBins);
	void write2file(float ***RDF_x_y, int nRDF, int nBins[], float binW[], char option[]);
	void write2file(float ***RDFxy, int nRDFtypes, int nBins[]);
	void write2file(float **Sk_x_y, float kx[], float ky[], int Nk[], int tag = 0);
	void write2file(float **Hist, int nBins[]);
	void write2file(float pe);
	void write2file(velocityDist *Dist);
	void write2file(float *power);
	};

	float computeNonBondedInteractions(atom_style *ATOMS, System *BOX, WCA_2P *INTERACTIONS);
	int bounds(float a);
}

#endif /*ANALYSIS_H*/