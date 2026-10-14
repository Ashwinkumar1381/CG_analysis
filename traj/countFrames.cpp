#include "analysis.h"

using namespace analysis;

int main(int argc, char *argv[])
{
	float dt = 1e-4;
	int frameW = int(1e3);

	Trajectory *TRAJ = new Trajectory(dt ,frameW, "cfg");
	sprintf(TRAJ->fpathI, "/home/ethaya_lab_c1/Desktop/ashwin_md/LD/lmp/jamming/tau_9e-1/Pe600/long_traj_res.cfg");
	
	TRAJ -> openTrajectory(true);

	TRAJ -> closeTrajectory();

	return(0);
}