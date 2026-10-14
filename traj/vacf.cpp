#include "analysis.h"

using namespace analysis;

void computeVACF(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples);
void computeVACF2(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples);
void compute_Directional_Corr(Trajectory *TRAJ, atom_style **ATOMS, System *BOX, int *frame_sep, int nSamples);

int main(int argc, char* argv[])
{
	float dt = 5e-4;
	int frameW = 10;
	int frameStart = int(0);
	int frameEnd = int(5e4);
	int nAtomTypes = 2;
	char *option = new char[50];
	sprintf(option, "vacf");

	int nSamples = 10000;
	int *frame_sep = new int[nSamples];
	for(int i = 0; i < nSamples; i++)
		frame_sep[i] = int(i*1);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0/R3/short_traj_res.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0/R3/vacf.dat");
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
			ATOMS[i] = new atom_style[TRAJ->nAtoms];

		TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, true);

		if(strcmp(option, "directional_corr") == 0)
			compute_Directional_Corr(TRAJ, ATOMS, BOX, frame_sep, nSamples);

		if(strcmp(option, "vacf") == 0)
			computeVACF2(TRAJ, ATOMS, frame_sep, nSamples);

		delete[] ATOMS;
	}

	TRAJ -> closeTrajectory(true, true);

}

void computeVACF(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples)
{
	double *correlations = new double[nSamples];
	long *counts = new long[nSamples];
	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		printf("\nProcessing frame %d", currFrame);

	// int currFrame = 0;

		for(int ctr = 0; ctr < nSamples; ctr++)
		{
			int nextFrame = currFrame + frame_sep[ctr];
			if(nextFrame < TRAJ->totalFrames)
			{
				for(int i = 0; i < TRAJ->nAtoms; i++)
				{
					double curr_velocity = sqrt(ATOMS[currFrame][i].vx*ATOMS[currFrame][i].vx + ATOMS[currFrame][i].vy*ATOMS[currFrame][i].vy);
					double next_velocity = sqrt(ATOMS[nextFrame][i].vx*ATOMS[nextFrame][i].vx + ATOMS[nextFrame][i].vy*ATOMS[nextFrame][i].vy);

					double prod = (ATOMS[currFrame][i].vx*ATOMS[nextFrame][i].vx + ATOMS[currFrame][i].vy*ATOMS[nextFrame][i].vy)/(curr_velocity*next_velocity);

					correlations[ctr] += prod;
					counts[ctr] += 1;
				}
			}
		}
	}

	for(int ctr = 0; ctr < nSamples; ctr++)
		correlations[ctr] /= double(counts[ctr]);

	TRAJ -> createOutputFile("frame_sep time count vacf");

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		// correlations[ctr] /= correlations[0];
		// printf("Frame_sep = %d, counts = %ld, vacf = %f\n", frame_sep[ctr], counts[ctr], correlations[ctr]);
		fprintf(TRAJ->fileO, "\n%d %f %ld %f", frame_sep[ctr], float(frame_sep[ctr]*TRAJ->frameWidth*TRAJ->timeStep), counts[ctr], correlations[ctr]);
	}

}

void computeVACF2(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples)
{
	double *correlations = new double[nSamples];
	long *counts = new long[nSamples];
	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		printf("\nProcessing frame %d", currFrame);

	// int currFrame = 0;

		#pragma omp parallel num_threads(10)
		{
			#pragma omp for
			for(int ctr = 0; ctr < nSamples; ctr++)
			{
				int nextFrame = currFrame + frame_sep[ctr];
				if(nextFrame < TRAJ->totalFrames)
				{
					for(int i = 0; i < TRAJ->nAtoms; i++)
					{
						double curr_velocity = sqrt(ATOMS[currFrame][i].vx*ATOMS[currFrame][i].vx + ATOMS[currFrame][i].vy*ATOMS[currFrame][i].vy);
						double next_velocity = sqrt(ATOMS[nextFrame][i].vx*ATOMS[nextFrame][i].vx + ATOMS[nextFrame][i].vy*ATOMS[nextFrame][i].vy);

						double prod = (ATOMS[currFrame][i].vx*ATOMS[nextFrame][i].vx + ATOMS[currFrame][i].vy*ATOMS[nextFrame][i].vy)/(curr_velocity*next_velocity);

						correlations[ctr] += prod;
						counts[ctr] += 1;
					}
				}
			}
		}
	}

	for(int ctr = 0; ctr < nSamples; ctr++)
		correlations[ctr] /= double(counts[ctr]);

	TRAJ -> createOutputFile("frame_sep time count vacf");

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		// correlations[ctr] /= correlations[0];
		// printf("Frame_sep = %d, counts = %ld, vacf = %f\n", frame_sep[ctr], counts[ctr], correlations[ctr]);
		fprintf(TRAJ->fileO, "\n%d %f %ld %f", frame_sep[ctr], float(frame_sep[ctr]*TRAJ->frameWidth*TRAJ->timeStep), counts[ctr], correlations[ctr]);
	}

}

void compute_Directional_Corr(Trajectory *TRAJ, atom_style **ATOMS, System *BOX, int *frame_sep, int nSamples)
{
	float *correlations = new float[nSamples];
	long *counts = new long[nSamples];

	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		int nextFrame = 0 + frame_sep[ctr];

		for(int i = 0; i < TRAJ->nAtoms; i++)
		{
			float dx_t1 = ATOMS[1][i].rxt1 - ATOMS[0][i].rxt1 + BOX->Lx*(ATOMS[1][i].jumpx - ATOMS[0][i].jumpx);
			float dy_t1 = ATOMS[1][i].ryt1 - ATOMS[0][i].ryt1 + BOX->Ly*(ATOMS[1][i].jumpy - ATOMS[0][i].jumpy);

			float dx_t2 = ATOMS[nextFrame + 1][i].rxt1 - ATOMS[nextFrame][i].rxt1 + BOX->Lx*(ATOMS[nextFrame + 1][i].jumpx - ATOMS[nextFrame][i].jumpx);
			float dy_t2 = ATOMS[nextFrame + 1][i].ryt1 - ATOMS[nextFrame][i].ryt1 + BOX->Ly*(ATOMS[nextFrame + 1][i].jumpy - ATOMS[nextFrame][i].jumpy);

			float disp_t1 = sqrt(dx_t1*dx_t1 + dy_t1*dy_t1);
			float disp_t2 = sqrt(dx_t2*dx_t2 + dy_t2*dy_t2);
			float prod = (dx_t1*dx_t2 + dy_t1*dy_t2)/(disp_t1*disp_t2);

			correlations[ctr] += prod;
			counts[ctr] += 1;
		}
	}

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		correlations[ctr] /= counts[ctr];
		printf("Frame_sep = %d, counts = %ld, vacf = %f\n", frame_sep[ctr], counts[ctr], correlations[ctr]);
		// fprintf(TRAJ->fileO, "\n%d %f %ld %f", frame_sep[ctr], float(frame_sep[ctr]*TRAJ->frameWidth*TRAJ->timeStep), counts[ctr], correlations[ctr]);
	}
}