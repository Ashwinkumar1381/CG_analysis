#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	char *option = new char[10];
	sprintf(option, "single");

	int particle_num = 825;
	int frameStart = int(3820), frameEnd = int(4940);
	int frameW = int(50);
	float dt = 5e-4;

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0_R2/short_traj.cfg");
	// sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/coords_tau_1e-3.dat");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/jamming/Fd100/tau_1e0_R2/coords_p%d.dat", particle_num);

	TRAJ -> openTrajectory();

	TRAJ -> totalFrames = int(frameEnd - frameStart + 1);
	atom_style *ATOMS = new atom_style[TRAJ->nAtoms];

	if(strcmp(option, "single") == 0)
	{
		TRAJ -> createOutputFile("step id element x y vx vy fx fy ix iy\n");

		char *fmt = new char[50];
		sprintf(fmt, "%%%s %s", "ld", TRAJ->line_fmt);

		while( !feof(TRAJ->fileI) )
		{
			TRAJ -> readThisFrame(ATOMS);

			if(TRAJ->frame_nr >= frameStart and TRAJ->frame_nr <= frameEnd)
			{
				int pid = particle_num - 1;

				fprintf(TRAJ->fileO, fmt, TRAJ->step, ATOMS[pid].atom_id, ATOMS[pid].element, ATOMS[pid].rxt1, ATOMS[pid].ryt1, ATOMS[pid].vx, ATOMS[pid].vy, ATOMS[pid].fx, ATOMS[pid].fy, ATOMS[pid].jumpx, ATOMS[pid].jumpy);
			}
		}

		printf("\n\nSucessfully imported coordinates for %d frames\n", TRAJ->totalFrames);
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

	TRAJ -> closeTrajectory(true, true);

	return(0);
}