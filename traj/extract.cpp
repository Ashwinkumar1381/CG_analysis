#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	char *option = new char[10];
	sprintf(option, "single");

	int num_particles = 1;
	int frameStart = int(3820), frameEnd = int(4940);
	int frameW = int(50);
	float dt = 5e-4;

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");

	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0_R2/short_traj.cfg");
	// sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/coords_tau_1e-3.dat");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0_R2/coords");
	
	if(argc != 1)
		num_particles = argc - 1;

	int *particle_nums = new int[num_particles];
	if(argc != 1)
	{
		int temp;
		for(int i = 1; i <= num_particles; i++)
		{
			sscanf(argv[i], "%d", &temp);
			particle_nums[i - 1] = temp;
		}
	}
	
	TRAJ -> openTrajectory();

	TRAJ -> totalFrames = int(frameEnd - frameStart + 1);
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];

	if(strcmp(option, "single") == 0)
	{
		printf("\nExtracting coordinates for %d particles from frame %d to frame %d\n", num_particles, frameStart, frameEnd);

		char *fmt = new char[50];
		sprintf(fmt, "%%%s %%%s %%%s %%%s %%%s\n", "ld", "d", "c", "g", "g");

		char *fpathO = new char[500];
		sprintf(fpathO, "%s", TRAJ->fpathO);

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr == frameStart)
			{
				for(int i = 0; i < num_particles; i++)
				{
					sprintf(TRAJ->fpathO, "%s/coords_p%d.dat", fpathO, particle_nums[i]);
					TRAJ -> createOutputFile("step id element x y\n");
					fclose(TRAJ->fileO);
				}
			}

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				printf("\nProcessing frame %d", TRAJ->frame_nr - frameStart + 1);

				for(int i = 0; i < num_particles; i++)
				{
					sprintf(TRAJ->fpathO, "%s/coords_p%d.dat", fpathO, particle_nums[i]);
					TRAJ->fileO = fopen(TRAJ->fpathO, "a+");

					int pid = particle_nums[i] - 1;

					fprintf(TRAJ->fileO, fmt, TRAJ->step, ATOMS[pid].atom_id, ATOMS[pid].element, ATOMS[pid].rxt1, ATOMS[pid].ryt1, ATOMS[pid].jumpx, ATOMS[pid].jumpy);

					fclose(TRAJ->fileO);
				}
			}

			if(TRAJ->frame_nr == frameEnd)
				break;
		}

		// printf("\n\nSucessfully imported coordinates for %d frames\n", TRAJ->totalFrames);
	}

	if(strcmp(option, "binary") == 0)
	{
		TRAJ -> createOutputFile("step x0 y0 jumpx0 jumpy0 vx vy x1 y1 jumpx1 jumpy1");

		for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
		{
			TRAJ -> readThisFrame(ATOMS);
			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
				fprintf(TRAJ->fileO, "\n%ld %f %f %d %d %f %f %f %f %d %d", TRAJ->step, ATOMS[0].rxt1, ATOMS[0].ryt1, ATOMS[0].jumpx, ATOMS[0].jumpy, ATOMS[0].vx, ATOMS[0].vy, ATOMS[1].rxt1, ATOMS[1].ryt1, ATOMS[1].jumpx, ATOMS[1].jumpy);
		}	
	}

	if(TRAJ -> fileO == NULL)
		TRAJ -> closeTrajectory(true, false);
	else
		TRAJ -> closeTrajectory(true, true);

	return(0);
}