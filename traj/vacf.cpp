#include "analysis.h"

using namespace analysis;

void computeVACF(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples, int frameMax, int total_threads = 1);
void computeVACF2(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples, int frameMax, int total_threads = 1);
void compute_Directional_Corr(Trajectory *TRAJ, atom_style **ATOMS, System *BOX, int *frame_sep, int nSamples, int frameMax, int total_threads = 1);
void compute_Residual_VACF(Trajectory *TRAJ, int frameStart, int frameEnd, float vd);

int main(int argc, char* argv[])
{
	int total_threads = 1;

	if(argc > 1)
	{
		sscanf(argv[1], "%d", &total_threads);
	}

	float dt = 5e-4;
	float tau_in = 1e-2, Ped = 100.0;
	int frameW = int(1);
	int frameStart = int(0);
	int frameEnd = int(4e4);
	int frameMax = int(4e4);
	int nAtomTypes = 2;
	char *option = new char[50];
	sprintf(option, "residual");

	int nSamples = int(5e3);
	int *frame_sep = new int[nSamples];
	for(int i = 0; i < nSamples; i++)
		frame_sep[i] = int(i*2);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/vacf/tau_1e-2/R2/short_traj.cfg");
	// sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/vacf/tau_1e-2/R2/vacf.dat");
	// sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/vacf/tau_1e-2/R2/pcorr.dat");
	TRAJ -> openTrajectory();

	if(frameEnd < frameStart)
	{
		printf("Error: Ending frame number cannot be smaller than the starting frame. Exiting...\n");
		exit(-1);
	}
	else
	{
		TRAJ->totalFrames = frameEnd - frameStart + 1;

		if(strcmp(option, "residual") == 0)
			compute_Residual_VACF(TRAJ, frameStart, frameEnd, float(Ped*tau_in));
		else
		{
			atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
			System *BOX = new System(TRAJ->Lx, TRAJ->Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
			for(int i = 0; i < TRAJ->totalFrames; i++)
				ATOMS[i] = new atom_style[TRAJ->nAtoms];

			TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd, true);

			if(strcmp(option, "pcorr") == 0)
				compute_Directional_Corr(TRAJ, ATOMS, BOX, frame_sep, nSamples, frameMax, total_threads);

			if(strcmp(option, "vacf") == 0)
				computeVACF(TRAJ, ATOMS, frame_sep, nSamples, frameMax, total_threads);

			delete[] ATOMS;
		}
	}

	TRAJ -> closeTrajectory(true, true);

}

void computeVACF(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples, int frameMax, int total_threads)
{
	double *correlations = new double[nSamples];
	long *counts = new long[nSamples];
	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int currFrame = 0; currFrame < frameMax; currFrame++)
	{
		printf("\nProcessing frame %d", currFrame);

	// int currFrame = 0;

		#pragma omp parallel num_threads(total_threads)
		{
			#pragma omp parallel for
			for(int ctr = 0; ctr < nSamples; ctr++)
			{
				int nextFrame = currFrame + frame_sep[ctr];
				if(nextFrame < TRAJ->totalFrames)
				{
					for(int i = 0; i < TRAJ->nAtoms; i++)
					{
						double prod = (ATOMS[currFrame][i].vx*ATOMS[nextFrame][i].vx + ATOMS[currFrame][i].vy*ATOMS[nextFrame][i].vy);

						correlations[ctr] += prod;
						counts[ctr] += 1;
					}
				}
			}
		}
	}

	for(int ctr = 0; ctr < nSamples; ctr++)
		correlations[ctr] /= double(counts[ctr]);

	double norm = correlations[0];

	for(int ctr = 0; ctr < nSamples; ctr++)
		correlations[ctr] /= norm;

	TRAJ -> createOutputFile("frame_sep time count vacf");

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		// correlations[ctr] /= correlations[0];
		// printf("Frame_sep = %d, counts = %ld, vacf = %f\n", frame_sep[ctr], counts[ctr], correlations[ctr]);
		fprintf(TRAJ->fileO, "\n%d %f %ld %f", frame_sep[ctr], float(frame_sep[ctr]*TRAJ->frameWidth*TRAJ->timeStep), counts[ctr], correlations[ctr]);
	}

}

void computeVACF2(Trajectory *TRAJ, atom_style **ATOMS, int *frame_sep, int nSamples, int frameMax, int total_threads)
{
	double *correlations = new double[nSamples];
	long *counts = new long[nSamples];
	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int currFrame = 0; currFrame < frameMax; currFrame++)
	{
		printf("\nProcessing frame %d", currFrame);

	// int currFrame = 0;

		#pragma omp parallel num_threads(total_threads)
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

						if(curr_velocity!= 0 and next_velocity != 0)
						{
							double prod = (ATOMS[currFrame][i].vx*ATOMS[nextFrame][i].vx + ATOMS[currFrame][i].vy*ATOMS[nextFrame][i].vy)/(curr_velocity*next_velocity);

							correlations[ctr] += prod;
							counts[ctr] += 1;
						}
					}
				}
			}
		}
	}

	TRAJ -> createOutputFile("frame_sep time count vacf");

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		if(counts[ctr] > 0)
			correlations[ctr] /= double(counts[ctr]);

		fprintf(TRAJ->fileO, "\n%d %f %ld %f", frame_sep[ctr], float(frame_sep[ctr]*TRAJ->frameWidth*TRAJ->timeStep), counts[ctr], correlations[ctr]);
	}

}

void compute_Directional_Corr(Trajectory *TRAJ, atom_style **ATOMS, System *BOX, int *frame_sep, int nSamples, int frameMax, int total_threads)
{
	double *correlations = new double[nSamples];
	long *counts = new long[nSamples];
	for(int i = 0; i < nSamples; i++)
	{
		correlations[i] = 0.0;
		counts[i] = 0;
	}

	for(int currFrame = 0; currFrame < frameMax - 1; currFrame++)
	{
		printf("\nProcessing frame %d", currFrame);
	
	// int currFrame = 0;

		#pragma omp parallel num_threads(total_threads)
		{
			#pragma omp for
			for(int ctr = 0; ctr < nSamples; ctr++)
			{
				int nextFrame = currFrame + frame_sep[ctr];
				if(nextFrame < TRAJ->totalFrames - 1)
				{
					for(int i = 0; i < TRAJ->nAtoms; i++)
					{
						double dx_t1 = ATOMS[currFrame + 1][i].rxt1 - ATOMS[currFrame][i].rxt1 + BOX->Lx*(ATOMS[currFrame + 1][i].jumpx - ATOMS[currFrame][i].jumpx);
						double dy_t1 = ATOMS[currFrame + 1][i].ryt1 - ATOMS[currFrame][i].ryt1 + BOX->Ly*(ATOMS[currFrame + 1][i].jumpy - ATOMS[currFrame][i].jumpy);

						double dx_t2 = ATOMS[nextFrame + 1][i].rxt1 - ATOMS[nextFrame][i].rxt1 + BOX->Lx*(ATOMS[nextFrame + 1][i].jumpx - ATOMS[nextFrame][i].jumpx);
						double dy_t2 = ATOMS[nextFrame + 1][i].ryt1 - ATOMS[nextFrame][i].ryt1 + BOX->Ly*(ATOMS[nextFrame + 1][i].jumpy - ATOMS[nextFrame][i].jumpy);

						double disp_t1 = sqrt(dx_t1*dx_t1 + dy_t1*dy_t1);
						double disp_t2 = sqrt(dx_t2*dx_t2 + dy_t2*dy_t2);

						if(disp_t1 != 0 and disp_t2 != 0)
						{
							double prod = (dx_t1*dx_t2 + dy_t1*dy_t2)/(disp_t1*disp_t2);

							correlations[ctr] += prod;
							counts[ctr] += 1;
						}
					}
				}
			}
		}
	}

	TRAJ -> createOutputFile("frame_sep time count pcorr");

	for(int ctr = 0; ctr < nSamples; ctr++)
	{
		if(counts[ctr] > 0)
			correlations[ctr] /= counts[ctr];
		
		fprintf(TRAJ->fileO, "\n%d %f %ld %f", frame_sep[ctr], float(frame_sep[ctr]*TRAJ->frameWidth*TRAJ->timeStep), counts[ctr], correlations[ctr]);
	}
}

void compute_Residual_VACF(Trajectory *TRAJ, int frameStart, int frameEnd, float vd)
{
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];
	float prod1 = 0.0;
	float prod2 = 0.0;
	long counts = 0;

	while( !feof(TRAJ->fileI) )
	{
		TRAJ -> readThisFrame(ATOMS);

		if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd and TRAJ->frame_nr%10 == 0)
		{
			printf("\nProcessing frame %d", TRAJ->frame_nr);

			for(int j = 0; j < TRAJ->nAtoms; j++)
			{
				prod1 += ATOMS[j].si*ATOMS[j].vx;
				prod2 += ATOMS[j].vx*ATOMS[j].vx + ATOMS[j].vy*ATOMS[j].vy;
				counts += 1;
			}
		}

		if(TRAJ->frame_nr == frameEnd)
			break;
	}

	prod1 /= counts;
	prod2 /= counts;

	printf("\nResidual correlations : %f, %f, %f, total_counts = %ld\n", prod1, prod2, (vd*prod1)/prod2, counts);

	delete[] ATOMS;
}