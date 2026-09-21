/*
	Computes the Lane order parameter as defined by the following articles:
	*) Wächtler, C.W. et al. (2016) Physical Review E, 94(5) - computeLaneOrder1()
	*) Dzubiella, J. et al. (2002) Physical Review E, 65(2) - computeLaneOrder2(), computeLaneOrder2B()

	Author 		  : Ashwin Kumar
	Last Modified : 13.03.26 
*/

#define NUM_THREADS 1
#include "analysis.h"

using namespace analysis;

float computeLaneOrder1(atom_style *ATOMS, System *BOX, float slabW_y);
float computeLaneOrder2(atom_style *ATOMS, System *BOX, float acut);
float computeLaneOrder2B(atom_style *ATOMS, System *BOX, float acut);

int main(int argc, char*argv[])
{
	char *option = new char[20];
	sprintf(option, "time_evolve_traj");

	// Trajectory params
	float dt = 1e-4;
	int frameW = int(1e3);
	int frameStart = int(0), frameEnd = int(5e4);

	// System params
	int nAtomTypes = 2;
	float Lx = 150.0, Ly = 30.0; 
	float slabW = 0.5;
	float rho = 0.45;
	float acut = 1.0/sqrt(rho);

	// Keywords and values for option = time_evolve_traj 
	int Nfreq = 1; 						// Obtain averages for every Nfreq frames
	int Nsample = 1; 					// Use Nsample frames to compute the average  
	int Nevery = 1;						// Obtain Nsample frames at intervals of Nevery before the current frame whose average we desire to compute

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_8e-1/Pe800/traj.cfg");

	if(strcmp(option, "time_evolve_traj") == 0)
	{
		sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_8e-1/Pe800/lane_order.dat");

		TRAJ -> openTrajectory();
		TRAJ -> createOutputFile("step order");

		atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

		int ctr = 0;
		float avg = 0.0;

		int nextFrame = frameStart + ctr*Nfreq;
		int frameStartAvg = nextFrame;

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				if(TRAJ->frame_nr == frameStart)
				{
					ctr++;
					nextFrame = frameStart + ctr*Nfreq;
					frameStartAvg = nextFrame - Nsample*Nevery;

					float order = computeLaneOrder1(ATOMS, BOX, slabW);
					TRAJ -> write2file(TRAJ -> time, order);
				}

				if((TRAJ->frame_nr > frameStartAvg) and (TRAJ->frame_nr <= nextFrame))
				{
					if((TRAJ->frame_nr - frameStartAvg)%Nevery == 0)
						avg += computeLaneOrder1(ATOMS, BOX, slabW);
				}

				if(TRAJ->frame_nr == nextFrame)
				{
					ctr++;
					nextFrame = frameStart + ctr*Nfreq;
					frameStartAvg = nextFrame - Nsample*Nevery;

					float order = float(avg/Nsample);
					TRAJ -> write2file(TRAJ -> time, order);

					avg = 0.0;
				}
			}

			if(TRAJ->frame_nr == frameEnd)
				break;
		}

		printf("\nComputed order parameter for %d frames.\n", ctr);

		TRAJ -> closeTrajectory();
	}

	else if(strcmp(option, "time_evolve_frames") == 0)
	{
		sprintf(TRAJ->fpathO, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data28/laneOrder.dat");
		TRAJ -> write2file(0.0, 0.0, 1);

		int total_frames1 = 944;
		int total_frames2 = 57;

		for(int i = 0; i < total_frames1; i++)
		{
			long curr_step = long(1e7)*(i + 1);

			printf("\nProcessing file %d / %d...", i+1, total_frames1 + total_frames2);

			sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data28/frames/frame.%ld.cfg", curr_step);
			TRAJ -> openTrajectory();

			atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
			System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

			TRAJ -> readThisFrame(ATOMS);
			float order = computeLaneOrder1(ATOMS, BOX, slabW);

			TRAJ -> write2file(0.0, order);

			delete[] ATOMS;
			delete BOX;
		}

		long add_step = TRAJ->step;

		for(int i = 0; i < total_frames2; i++)
		{
			long curr_step = long(1e7)*(i + 1);

			printf("\nProcessing file %d / %d...", i+1+total_frames1, total_frames1 + total_frames2);

			sprintf(TRAJ->fpathI, "//media/ashwin/One Touch/ashwin_md/lane/May2025/lmp/Data28/frames_res1/frame.%ld.cfg", curr_step);
			TRAJ -> openTrajectory();

			atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
			System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

			TRAJ -> readThisFrame(ATOMS);
			float order = computeLaneOrder1(ATOMS, BOX, slabW);

			TRAJ->step = curr_step + add_step;
			TRAJ -> write2file(0.0, order);
		}

		TRAJ -> write2file(0.0, 0.0, 0);
	}

	else if(strcmp(option, "time_avg_traj") == 0)
	{
		TRAJ -> openTrajectory();
		atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
		System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

		int ctr = 0;
		float avg = 0.0, max = 0.0, min = 1.0;

		while( !feof(TRAJ -> fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				float order = computeLaneOrder1(ATOMS, BOX, slabW); 
				avg += order;

				if(order > max) max = order;
				if(order < min) min = order;

				ctr++;
			}

			if(TRAJ->frame_nr == frameEnd)
			{
				avg /= ctr;
				break;
			}
		}

		printf("\nLane order parameter averaged over %d frames.\n\nMax Min Avg\n%f %f %f\n", ctr, max, min, avg);
		TRAJ -> closeTrajectory();
	}
}

/* Based on the lane order parameter as defined by Wächtler, C.W. et al. (2016) Physical Review E, 94(5) */
float computeLaneOrder1(atom_style *ATOMS, System *BOX, float slabW_y)
{
	int nBins = int(BOX->Ly/slabW_y);

	float phi[nBins];
	int bin_ctr[nBins];
	for(int i = 0; i < nBins; i++) 
	{
		phi[i] = 0.0;
		bin_ctr[i] = 0;
	}

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		int bin_id = int(ATOMS[i].ryt1 / slabW_y);

		if(ATOMS[i].element == 'N')
			phi[bin_id] += -1.0;

		else if(ATOMS[i].element == 'O')
			phi[bin_id] += 1.0;

		bin_ctr[bin_id] += 1;
	}

	float order = 0.0;
	for(int i = 0; i < nBins; i++)
	{
		if(bin_ctr[i] != 0)
			order += abs(phi[i] / bin_ctr[i]);
		else
			order += 0.0;
	}

	return(order/nBins);
}

/* Based on the lane order parameter as defined by Dzubiella, J. et al. (2002) Physical Review E, 65(2) */
float computeLaneOrder2(atom_style *ATOMS, System *BOX, float acut)
{
	float avg = 0.0;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		int add = 1;

		for(int j = 0; j < BOX->nAtoms; j++)
		{
			if(ATOMS[i].element != ATOMS[j].element)
			{
				float dxij = ATOMS[i].rxt1 - ATOMS[j].rxt1;
				float dyij = ATOMS[i].ryt1 - ATOMS[j].ryt1;
				BOX -> checkMinImage(&dxij, &dyij);

				if(dxij*dxij + dyij*dyij <= acut*acut)
				{
					// printf("Atoms %d(%c) and %d(%c) are separated by %f\n", i, ATOMS[i].element, j, ATOMS[j].element, dyij);
					add = 0;
					break;
				}

			}
		}

		avg += float(add);
	}

	return(avg/BOX->nAtoms);
}

float computeLaneOrder2B(atom_style *ATOMS, System *BOX, float acut)
{
	int **atomList = new int*[BOX->nAtomTypes];
	int *ctr = new int[BOX->nAtomTypes];

	int nAtomsPerType = int(BOX->nAtoms/BOX->nAtomTypes);

	for(int i = 0; i < BOX->nAtomTypes; i++)
	{
		atomList[i] = new int[nAtomsPerType];
		for(int j = 0; j < nAtomsPerType; j++)
			atomList[i][j] = 0;

		ctr[i] = 0;
	}

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		if(ATOMS[i].element == 'N')
		{
			atomList[0][ctr[0]] = i;
			ctr[0]++;
		}
		else if(ATOMS[i].element == 'O')
		{
			atomList[1][ctr[1]] = i;
			ctr[1]++;
		}
	}

	float avg = 0.0;
	#pragma omp parallel num_threads(NUM_THREADS)
	{
		int sum = 0;

		#pragma omp for
		for(int i = 0; i < nAtomsPerType; i++)
		{
			int iid = atomList[0][i];
			int add = 1;

			for(int j = 0; j < nAtomsPerType; j++)
			{
				int jid = atomList[1][j];

				float dyij = ATOMS[iid].ryt1 - ATOMS[jid].ryt1;
				BOX -> checkMinImage(NULL, &dyij);
				// printf("%f\n", dyij);

				if(abs(dyij) < acut) 
				{
					// printf("Hi %d\n", i);
					add = 0;
					break;
				}
			} 

			sum += add;
		}

		#pragma omp critical
		avg += sum;	
	}

	// exit(-1);

	#pragma omp parallel num_threads(NUM_THREADS)
	{
		int sum = 0;

		#pragma omp for
		for(int i = 0; i < nAtomsPerType; i++)
		{
			int iid = atomList[1][i];
			int add = 1;

			for(int j = 0; j < nAtomsPerType; j++)
			{
				int jid = atomList[0][j];

				if(abs(ATOMS[iid].ryt1 - ATOMS[jid].ryt1) < acut)
				{
					add = 0;
					break;
				}
			}

			sum += add;
		}

		#pragma omp critical
		avg += sum;
	}

	printf("%f\n", avg);

	return(avg/BOX->nAtoms);
}

void analysis::Trajectory::write2file(float time, float order, int option)
{
	fprintf(fileO, "\n%ld %f", step, order);
}