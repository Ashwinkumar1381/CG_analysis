/*
	Code to compute the Mean-Squared-Displacement of a group of atoms
	by averaging over each particle's trajectory and over all particles
	in the group.

	Author			: Ashwin Kumar
	Date created 	: 15.03.25
	Last modified	: 20.02.26
*/

#include "analysis.h"

using namespace analysis;

void computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, char group[], int delFrames[], int nSample);

int main(int argc, char *argv[])
{
	// ----------- System Params -----------
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	// ----------- Trajectory Params -----------
	int frameStart = int(15e3), frameEnd = int(5e4);
	float dt = 1e-4;
	int frameW = int(2e4);

	char group[5] = "all";
	int nSample = int(43);

	int delFrames[nSample] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600,
							  700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 35000, 40000};

	// int delFrames[nSample];
	// for(int i = 0; i < 199; i++) delFrames[i] = int(1 + 1*i);
	// for(int i = 0; i < 180; i++) delFrames[i + 199] = int(200 + 10*i);
	// for(int i = 0; i <= 80; i++) delFrames[i + 379] = int(2000 + 100*i);

	// for(int i = 0; i < nSample; i++) printf("\ndelFrame = %d", delFrames[i]);

	// exit(0);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");

	if(argc == 1)
	{
		sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_5e-1/Pe30/traj.cfg");
		sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_5e-1/Pe30/msd.dat");	
	}
	else
	{
		int temp;

		sprintf(TRAJ -> fpathI, argv[1]);
		sprintf(TRAJ -> fpathO, argv[2]);
		
		sscanf(argv[3], "%d", &temp);
		frameStart = temp;
		
		sscanf(argv[4], "%d", &temp);
		frameEnd = temp;

		sscanf(argv[5], "%d", &temp);
		TRAJ->frameWidth = temp;
	}

	TRAJ -> openTrajectory();

	if(frameEnd < frameStart)
	{
		printf("Error: Ending frame number cannot be smaller than the starting frame. Exiting...\n");
		exit(-1);
	}
	else
	{
		TRAJ->totalFrames = frameEnd - frameStart + 1;

		atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
		for(int i = 0; i < TRAJ->totalFrames; i++) 
			ATOMS[i] = new atom_style [TRAJ->nAtoms];

		TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, true); 
		computeMeanSquaredDisplacement(TRAJ, BOX, ATOMS, group, delFrames, nSample);

		delete[] ATOMS;
	}

	TRAJ -> closeTrajectory(true, true);

	return(0);
}

void computeMeanSquaredDisplacement(Trajectory *TRAJ, System *BOX, atom_style **ATOMS, char group[], int delFrames[], int nSample)
{
	double **MSD = new double*[4];
	long *count = new long[nSample];

	for(int i = 0; i < 4; i++)
	{
		MSD[i] = new double[nSample];
		for(int j = 0; j < nSample; j++)
			MSD[i][j] = 0.0;
	}

	for(int i = 0; i < nSample; i++)
		count[i] = 0;

	long **group_lists = TRAJ->sortAtomsByType(ATOMS[0], BOX);

	int loopMax, group_id = 0;

	if(strcmp(group, "all") == 0)
		loopMax = TRAJ->nAtoms;
	else
	{
		sscanf(group, "%d\n", &group_id);
		loopMax = group_lists[0][group_id - 1];
	}

	printf("\nComputing MSD...");

	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		// printf("\nComputing MSD for frame %d...", currFrame);

		int ctr = 0;

		while(ctr < nSample)
		{
			int nextFrame = currFrame + delFrames[ctr];
			
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int j = 0; j < loopMax; j++)
				{
					int k = j;
					if(group_id != 0)
						k = group_lists[group_id][j] - 1;

					double dx = ATOMS[nextFrame][k].rxt1 - ATOMS[currFrame][k].rxt1 + BOX->Lx*(ATOMS[nextFrame][k].jumpx - ATOMS[currFrame][k].jumpx);
					double dy = ATOMS[nextFrame][k].ryt1 - ATOMS[currFrame][k].ryt1 + BOX->Ly*(ATOMS[nextFrame][k].jumpy - ATOMS[currFrame][k].jumpy);
					double dz = 0.0;
					if(TRAJ->dimension == 3)
						dz = ATOMS[nextFrame][k].rzt1 - ATOMS[currFrame][k].rzt1 + BOX->Lz*(ATOMS[nextFrame][k].jumpz - ATOMS[currFrame][k].jumpz);
					double dx2 = dx*dx;
					double dy2 = dy*dy;
					double dxdy = abs(dx)*abs(dy);

					if(TRAJ->dimension == 2)
					{
						MSD[0][ctr] += (dx2 + dy2);
					}

					else if(TRAJ->dimension == 3)
					{
						double dz2 = dz*dz;
						MSD[0][ctr] += (dx2 + dy2 + dz2);
					}

					MSD[1][ctr] += dx2;
					MSD[2][ctr] += dy2;
					MSD[3][ctr] += dxdy;
					
					count[ctr] += 1;
				}
			}

			ctr++;
		}
	}

	for(int i = 0; i < nSample; i++)
	{
		for(int j = 0; j < 4; j++)
			MSD[j][i] /= double(count[i]);
	}

	TRAJ -> createOutputFile("time count MSD_all MSD_xx MSD_yy MSD_xy");
	TRAJ -> write2file(MSD, count, delFrames, nSample);
}

void analysis::Trajectory::write2file(double **MSD, long *count, int *delFrames, int nSample)
{
	for(int i = 0; i < nSample; i++)
	{
		fprintf(fileO, "\n%g %ld %f %f %f %f", delFrames[i]*frameWidth*timeStep, count[i], MSD[0][i], MSD[1][i], MSD[2][i], MSD[3][i]);
		// printf("\n%g %ld %f %f %f %f", delFrames[i]*frameWidth*timeStep, count[i], MSD[0][i], MSD[1][i], MSD[2][i], MSD[3][i]);
	}

	fprintf(fileO, "\n");
}