/*
	Extract specific frames from multiple trajectories and combine into one trajectory

	Last modified : 14.09.25
*/

#include "analysis.h"

using namespace analysis;

int main(int argc, char* argv[])
{
	/* ---------- System Params ---------- */
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	/* ---------- Trajectory Params ---------- */
	float dt = 1e-4;
	int frameW = int(2e4);
	
	int nFiles = 2;
	long startStep[nFiles] = {long(1e7+1.0e9), long(0)};

	long endStep[nFiles]   = {long(1e7+2.0e9), long(2e9)};

	string filenames[nFiles] = {"long_traj", "long_traj_res"};

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathO, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/lmp/jamming/tau_4e-1/Pe32/long_traj_merged.%s", TRAJ->format);
	TRAJ -> createOutputFile();

	atom_style *ATOMS;
	System *BOX;

	for(int i = 0; i < nFiles; i++)
	{
		sprintf(TRAJ->fpathI, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/lmp/jamming/tau_4e-1/Pe32/%s.%s", (filenames[i]).c_str(), TRAJ->format);
		TRAJ -> openTrajectory();

		long addStep = 0;

		if(i == 0)
		{
			ATOMS = new atom_style[TRAJ->nAtoms];
			BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

			TRAJ->totalFrames = 0;
		}
		else
			addStep = endStep[i - 1];

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->step >= startStep[i] and TRAJ->step <= endStep[i])
				TRAJ -> writeThisFrame(ATOMS, BOX, addStep);

			if(TRAJ->step == endStep[i]) 
			{
				endStep[i] += addStep;
				TRAJ->totalFrames += TRAJ->frame_nr;

				printf("\nSource: %s\nDestination: %s\nFrames extracted this step: %d\nTotal frames extracted: %d\n", TRAJ->fpathI, TRAJ->fpathO, TRAJ->frame_nr, TRAJ->totalFrames);
				break;
			}
		}

		if(i < nFiles - 1)
			TRAJ -> closeTrajectory(true, false);
		else
			TRAJ -> closeTrajectory(true, true);
	}
}