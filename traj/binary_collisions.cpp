#include "analysis.h"

using namespace std;
using namespace analysis;

int main(int argc, char *argv[])
{
	float Lx = 100.0, Ly = 100.0;
	int nAtomTypes = 2;

	int frameStart = int(0), frameEnd = int(1e4);
	float dt = 5e-4;
	int frameW = int(1);

	int run_nr = 0;
	sscanf(argv[1], "%d", &run_nr);

	Trajectory *TRAJ = new Trajectory(dt, frameW, "cfg");
	sprintf(TRAJ->fpathI, "//media/ashwin/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/traj.cfg");
	sprintf(TRAJ->fpathO, "//media/ashwin/ASH_DRIVE_3/ashwin_md/binary_collisions/tests/collision_stats.dat");

	TRAJ -> openTrajectory();

	if(run_nr == 1)
		TRAJ -> createOutputFile("run_nr tau_col dely_init delx_dist dely_dist delx_disp dely_disp vel_corr", true);
	else
		TRAJ -> createOutputFile("", false);

	TRAJ->totalFrames = frameEnd - frameStart + 1;
	atom_style **ATOMS = new atom_style*[TRAJ->totalFrames];
	System *BOX = new System(Lx, Ly, TRAJ->Lz, TRAJ->nAtoms, nAtomTypes);
	for(int i = 0; i < TRAJ->totalFrames; i++)
		ATOMS[i] = new atom_style[TRAJ->nAtoms];

	TRAJ -> loadTrajectory(ATOMS, BOX, frameStart, frameEnd);

	float init_pos_x = 0.0, init_pos_y = 0.0;
	float dely_init = 0.0, dist_x = 0.0, dist_y = 0.0, disp_x = 0.0, disp_y = 0.0;

	bool pos_track = false, vel_track = true;
	int nFrames = 0, tau_col_frames = 0;

	float *u_init = new float[2];
	float *u_final = new float[2];
	for(int i = 0; i < 2; i++)
	{
		u_init[i] = 0.0;
		u_final[i] = 0.0;
	}
	
	for(int currFrame = 0; currFrame < TRAJ->totalFrames; currFrame++)
	{
		float dxij = ATOMS[currFrame][1].rxt1 - ATOMS[currFrame][0].rxt1 + BOX->Lx*(ATOMS[currFrame][1].jumpx - ATOMS[currFrame][0].jumpx);
		float dyij = ATOMS[currFrame][1].ryt1 - ATOMS[currFrame][0].ryt1 + BOX->Ly*(ATOMS[currFrame][1].jumpy - ATOMS[currFrame][0].jumpy);

		if(pos_track == false and vel_track == true)
		{
			if(currFrame == 0)
			{
				init_pos_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx;
				init_pos_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy;
			}

			if(currFrame != 0)
			{
				float curr_pos_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx;
				float curr_pos_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy;
				nFrames += 1;

				u_init[0] += float((curr_pos_x - init_pos_x)/(nFrames*frameW*dt));
				u_init[1] += float((curr_pos_y - init_pos_y)/(nFrames*frameW*dt));
			}

			if(dxij < 1.0 and abs(dyij) < 1.0)
			{
				pos_track = true;
				vel_track = false;

				u_init[0] /= nFrames;
				u_init[1] /= nFrames;
				nFrames = 0;

				disp_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx;
				disp_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy;
				dely_init = dyij;
			}
		}

		if(pos_track == true and vel_track == false)
		{
			float dxTracer = ATOMS[currFrame][0].rxt1 - ATOMS[currFrame - 1][0].rxt1 + BOX->Lx*(ATOMS[currFrame][0].jumpx - ATOMS[currFrame - 1][0].jumpx);
			float dyTracer = ATOMS[currFrame][0].ryt1 - ATOMS[currFrame - 1][0].ryt1 + BOX->Ly*(ATOMS[currFrame][0].jumpy - ATOMS[currFrame - 1][0].jumpy);

			dist_x += abs(dxTracer);
			dist_y += abs(dyTracer);
			nFrames += 1;

			if(abs(dyij) >= 1.0)
			{
				vel_track = true;

				disp_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx - disp_x;
				disp_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy - disp_y;
				tau_col_frames = nFrames;
				nFrames = 0;

				init_pos_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx;
				init_pos_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy;
			}
		}

		if(pos_track == true and vel_track == true)
		{
			float curr_pos_x = ATOMS[currFrame][0].rxt1 + BOX->Lx*ATOMS[currFrame][0].jumpx;
			float curr_pos_y = ATOMS[currFrame][0].ryt1 + BOX->Ly*ATOMS[currFrame][0].jumpy;
			nFrames += 1;

			u_final[0] += float((curr_pos_x - init_pos_x)/(nFrames*frameW*dt));
			u_final[1] += float((curr_pos_y - init_pos_y)/(nFrames*frameW*dt));

			if(nFrames == int(2*tau_col_frames))
			{
				u_final[0] /= nFrames;
				u_final[1] /= nFrames;

				float vel_init = sqrt(u_init[0]*u_init[0] + u_init[1]*u_init[1]);
				float vel_final = sqrt(u_final[0]*u_final[0] + u_final[1]*u_final[1]);
				float vel_corr = (u_init[0]*u_final[0] + u_init[1]*u_final[1])/(vel_init*vel_final);

				float tau_col = float(tau_col_frames*frameW*dt);

				fprintf(TRAJ->fileO, "\n%d %g %g %g %g %g %g %g", run_nr, tau_col, dely_init, dist_x, dist_y, disp_x, disp_y, vel_corr);

				delete[] ATOMS;

				TRAJ->closeTrajectory(true, true);
				break;
			}
		}
	}

	return(0);
}