// temperature.cpp

#include "analysis.h"
#define NUM_THREADS 1

using namespace analysis;

void getVelocities(Trajectory *TRAJ, atom_style *ATOMS, System *BOX);
float* computeKineticTemperature(atom_style *ATOMS, System *BOX);
float computeThermalTemperature(atom_style *ATOMS, System *BOX);

int main(int argc, char *argv[])
{
	char *option = new char[20];
	sprintf(option, "get_velocity_frame");

	/* -------- System Params -------- */
	float Lx = 150.0, Ly = 30.0;
	int nAtomTypes = 2;

	/* -------- Trajectory Params -------- */
	float dt = 1e-5;
	int frameW = int(5e4);
	int frameStart = int(1468), frameEnd = int(5e4);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_1e0/Pe5/dt_1e-5/traj.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/tau_1e0/Pe5/dt_1e-5/");

	TRAJ -> openTrajectory();

	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);

	if(strcmp(option, "time_evolve_traj") == 0)
	{
		sprintf(TRAJ->fpathO, "%stemp_all.dat", TRAJ->fpathO);
		TRAJ -> createOutputFile("step temp temp_x temp_y");

		int ctr = 0;
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				float *temp = computeKineticTemperature(ATOMS, BOX);
				TRAJ -> write2file(temp);
				
				ctr++;
			}

			if(TRAJ->frame_nr == frameEnd)
				break;
		}

		printf("Computed kinetic temperature for %d frames.\n", ctr);
	}

	if(strcmp(option, "time_avg_traj") == 0)
	{
		int choice = 0;
		float min = 1e4, max = 0.0, avg = 0.0;

		int ctr = 0;
		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				float *temp = computeKineticTemperature(ATOMS, BOX); 
				avg += temp[choice];

				if(temp[choice] < min) min = temp[choice];
				if(temp[choice] > max) max = temp[choice];

				ctr++;
			}

			if(TRAJ->frame_nr == frameEnd)
				break;
		}

		avg /= ctr;

		printf("\nTemperature averaged over %d frames.\nMin Max Avg\n%f %f %f\n", ctr, min, max, avg);
	}

	else if(strcmp(option, "get_velocity_frame") == 0)
	{
		sprintf(TRAJ->fpathO, "%svelocities_frame%d.dat", TRAJ->fpathO, frameStart);

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr == frameStart)
			{
				getVelocities(TRAJ, ATOMS, BOX);

				printf("Frame %d exported!\n", frameStart);	
				exit(-1);	
			}
		}
	}
	
	TRAJ -> closeTrajectory(true, true);

	delete[] ATOMS;
	delete TRAJ;
	return(0);
}

void getVelocities(Trajectory *TRAJ, atom_style *ATOMS, System *BOX)
{
	TRAJ->fileO = fopen(TRAJ->fpathO, "w");
	if(TRAJ->fileO == NULL)
	{
		printf("Cannot open file. exiting...");
		exit(-1);
	}
	else
		fprintf(TRAJ->fileO, "id type rx ry vx vy\n");

	for(int i = 0; i < TRAJ->nAtoms; i++)
		fprintf(TRAJ->fileO, "%d %d %g %g %g %g\n", ATOMS[i].atom_id, ATOMS[i].type, ATOMS[i].rxt1, ATOMS[i].ryt1, ATOMS[i].vx, ATOMS[i].vy);

	fclose(TRAJ->fileO);
}

float* computeKineticTemperature(atom_style *ATOMS, System *BOX)
{
	float *vCom = new float[2];

	for(int i = 0; i < 2; i++)
		vCom[i] = 0.0;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		vCom[0] += ATOMS[i].vx;
		vCom[1] += ATOMS[i].vy; 
	}

	for(int i = 0; i < 2; i++)
		vCom[i] /= BOX->nAtoms;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		ATOMS[i].vxth = ATOMS[i].vx - vCom[0];
		ATOMS[i].vyth = ATOMS[i].vy - vCom[1];
	}

	float *temp = new float[3];
	for(int i = 0; i < 3; i++)
		temp[i] = 0.0;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		float vx2 = ATOMS[i].vxth*ATOMS[i].vxth;
		float vy2 = ATOMS[i].vyth*ATOMS[i].vyth;
		temp[0] += vx2 + vy2;
		temp[1] += vx2;
		temp[2] += vy2;
	}

	temp[0] /= (2.0*BOX->nAtoms);
	temp[1] /= (BOX->nAtoms);
	temp[2] /= (BOX->nAtoms);
	
	return(temp);
}

float computeThermalTemperature(atom_style *ATOMS, System *BOX)
{
	float **vCom = new float*[BOX->nAtomTypes];
	for(int i = 0; i < BOX->nAtomTypes; i++)
	{
		vCom[i] = new float[2];
		for(int j = 0; j < 2; j++)
			vCom[i][j] = 0.0;
	}

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		int ptype = ATOMS[i].type - 1;
		vCom[ptype][0] += ATOMS[i].vx;
		vCom[ptype][1] += ATOMS[i].vy;
	}

	for(int i = 0; i < BOX->nAtomTypes; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			vCom[i][j] /= (BOX->nAtoms/BOX->nAtomTypes);
		}
	}

	float ke = 0.0;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		int ptype = ATOMS[i].type - 1;

		ATOMS[i].vxth = ATOMS[i].vx - vCom[ptype][0];
		ATOMS[i].vyth = ATOMS[i].vy - vCom[ptype][1];

		ke += (ATOMS[i].vxth*ATOMS[i].vxth + ATOMS[i].vyth*ATOMS[i].vyth);
	}

	ke /= (2.0*BOX->nAtoms);
	
	return(ke);
}

void analysis::Trajectory::write2file(float *temp)
{
	fprintf(fileO, "\n%ld %g %g %g", step, temp[0], temp[1], temp[2]);
}