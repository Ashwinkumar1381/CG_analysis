/*
	analysis.cpp

	Some utility classes to perform analysis on Trajectories

	Date created  : 28.11.24
	Last modified : 20.02.26
*/

#include "analysis.h"

using namespace analysis;

/* ----------------- atomsXYZ members ----------------- */

analysis::atomsXYZ::atomsXYZ()
{
	rxt1 = ryt1 = rxt2 = ryt2 = 0.0;
	vx = vy = vz = vxth = vyth = vzth = 0.0;
	fx = fy = fx_int = fy_int = 0.0;
	si = 1;
	jumpx = jumpy = 0;
}

analysis::atomsXYZ::~atomsXYZ(){}

/* ----------------- molecularXYZ members ----------------- */

analysis::molecularXYZ::molecularXYZ()
{
	mol_id = 1;
	rzt1 = 0.0;
	jumpz = 0;
}

analysis::molecularXYZ::~molecularXYZ(){}

/* ----------------- System members ----------------- */

analysis::System::System(float Lx, float Ly, float Lz, int nAtoms, int nAtomTypes, bool buildMaps, float rcellx, float rcelly, float rcellz)
{
	this -> Lx = Lx;
	this -> Ly = Ly;
	this -> Lz = Lz;
	this -> nAtoms = nAtoms;
	this -> nAtomTypes = nAtomTypes;
	this -> rcellx = rcellx;
	this -> rcelly = rcelly;
	this -> rcellz = rcellz;

	Ncellx = int(this->Lx/this->rcellx);
	Ncelly = int(this->Ly/this->rcelly);

	this->rcellx = this->Lx/Ncellx;
	this->rcelly = this->Ly/Ncelly;
	
	ncells = int(Ncellx*Ncelly);

	for(int i = 0; i < MAXCELL; i++)
	{
		MAPS[i] = 0;
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	if(buildMaps == true)
		buildCellMaps();
}

analysis::System::~System(){}

int analysis::System::cellindex(int ix, int iy)
{   
    if (ix >= Ncellx) ix = ix - Ncellx; 
	else if (ix <= -1) ix = ix + Ncellx; 
	if (iy >= Ncelly) iy = iy - Ncelly;
	else if (iy <= -1) iy = iy + Ncelly;

	if(Ncellx >= Ncelly)
		return (1 + ix + iy*Ncellx);
	else
		return (1 + ix + iy*Ncelly);
}

void analysis::System::buildCellMaps()
{
	for (int ix = 0; ix < Ncellx; ix++)
	{
		for (int iy = 0; iy < Ncelly; iy++)
		{
			int imap = 4*(cellindex(ix,iy)-1);
			MAPS[imap+1] = cellindex(ix-1,iy);
			MAPS[imap+2] = cellindex(ix-1,iy+1);
			MAPS[imap+3] = cellindex(ix,iy+1);
			MAPS[imap+4] = cellindex(ix+1,iy+1);
		}
	}
    
    printf("\nSuccessfully constructed MAPS array with %d cells.\n", int(Ncellx * Ncelly));
}

void analysis::System::buildCellList(atom_style *ATOMS)
{
	for(int i = 0; i < MAXCELL; i++)
	{
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	for (int i = 1; i <= nAtoms; i++)
	{
		int ii = i - 1;
		int ix = int(ATOMS[ii].rxt1/rcellx);
		int iy = int(ATOMS[ii].ryt1/rcelly);
		
		int icell = cellindex(ix, iy);	        
		LIST[i] = HEAD[icell];
		HEAD[icell] = i;		    
	}	
}

void analysis::System::checkMinImage(float *dx, float *dy, float *dz)
{
	if(dx != NULL)
	{
		if(*dx >= 0.5*Lx) *dx -= Lx;
		else if(*dx <= -0.5*Lx) *dx += Lx;
	}
	if(dy != NULL)
	{
		if(*dy >= 0.5*Ly) *dy -= Ly;
		else if(*dy <= -0.5*Ly) *dy += Ly; 
	}
	if(dz != NULL)
	{
		if(*dz >= 0.5*Lz) *dz -= Lz;
		else if(*dz <= -0.5*Lz) *dz += Lz;
	}
}

/* ----------------- Bin1D members -----------------*/

analysis::Bin1D::Bin1D(float Lx, float Ly, float binWidth)
{
	this -> Lx = Lx;
	this -> Ly = Ly;
	this -> binWidth = binWidth;

	nBins = int(Lx/binWidth);

	bin = new float[nBins];
	binavg = new float[nBins];

	zero(bin);
	zero(binavg);
}

analysis::Bin1D::~Bin1D(){}

void analysis::Bin1D::addToBin(float rx, float val)
{
	if(val == 1.0)
		bin[int(rx/binWidth)] += 1;
	else
		bin[int(rx/binWidth)] += val;	
}

void analysis::Bin1D::normalize(float *binTonorm, float *normv, float normc)
{
	if(normv != NULL)
	{
		for(int i = 0; i < nBins; i++)
			binTonorm[i] /= normv[i];
	}
	else
	{
		if(normc == 1.0) normc = Ly*binWidth;

		for(int i = 0; i < nBins; i++)
			binTonorm[i] /= normc;	
	}
}

void analysis::Bin1D::zero(float *binTozero)
{
	for(int i = 0; i < nBins; i++)
		binTozero[i] = 0.0;
}

void analysis::Bin1D::addBins(float *bin1, float *bin2)
{
	for(int i = 0; i < nBins; i++)
		bin1[i] += bin2[i];
}

/* ----------------- Hist2D members -----------------*/

analysis::Hist2D::Hist2D(float rcut_x, float rcut_y, float binW_x, float binW_y)
{
	this -> rcut_x = rcut_x;
	this -> rcut_y = rcut_y;
	this -> binW_x = binW_x;
	this -> binW_y = binW_y;

	nBin_x = 2*int(rcut_x / binW_x);
	nBin_y = 2*int(rcut_y / binW_y);

	bin = new float*[nBin_x];
	for(int i = 0; i < nBin_x; i++)
	{
		bin[i] = new float[nBin_y];

		for(int j = 0; j < nBin_y; j++)
			bin[i][j] = 0.0;
	}
} 

analysis::Hist2D::~Hist2D(){}

void analysis::Hist2D::addToBin(char iid, char jid, float dxij, float dyij)
{
	int bin_xi = int((rcut_x + dxij)/binW_x);
	int bin_xj = int((rcut_x - dxij)/binW_x);
	int bin_yi = int((rcut_y + dyij)/binW_y);
	int bin_yj = int((rcut_y - dyij)/binW_y);

	if(bin_xi == nBin_x) bin_xi--;
	if(bin_xj == nBin_x) bin_xj--;
	if(bin_yi == nBin_y) bin_yi--;
	if(bin_yj == nBin_y) bin_yj--;

	if(iid == jid)
	{
		bin[bin_xi][bin_yi] += 1.0;
		bin[bin_xj][bin_yj] += 1.0;
	}
	else
		bin[bin_xi][bin_yi] += 1.0;	
}

void analysis::Hist2D::normalize(float fac)
{
	for(int i = 0; i < nBin_x; i++)
	{
		for(int j = 0; j < nBin_y; j++)
			bin[i][j] /= (fac * binW_x * binW_y);
	}
}

/* ----------------- Trajectory members -----------------*/

analysis::Trajectory::Trajectory(float timeStep, int frameWidth, char fileFormat[10])
{
	fpathI = new char [500];
	fpathO = new char [500];
	pipeString = new char [500];
	pipeChar = new char [500];
	fileI = nullptr;
	fileO = NULL;

	frame_nr = -1;
	this->timeStep = timeStep;
	this->frameWidth = frameWidth;
	format = fileFormat;
	nAtoms = nAtomTypes = nMols = 0;
	Lx = Ly = Lz = 0.0;
	xCom = yCom = zCom = 0.0;
}

analysis::Trajectory::~Trajectory(){}

void analysis::Trajectory::openTrajectory(bool count)
{
	fileI = fopen(fpathI, "r");
	if(fileI == NULL)
	{
		printf("Could not open %s. Exiting ...\n", fpathI);
		exit(-1);
	}
	else
	{
		if(strcmp(format, "xyz") == 0)
		{
			fgets(pipeString, 500, fileI);
			sscanf(pipeString, "%d", &nAtoms);			
		}
		else if(strcmp(format, "cfg") == 0)
		{
			for(int i = 0; i < 4; i++) 
				fgets(pipeString, 500, fileI);

			sscanf(pipeString, "%d", &nAtoms);

			fgets(pipeString, 500, fileI);

			float x1, x2 = 0.0;

			fgets(pipeString, 500, fileI);
			sscanf(pipeString, "%f %f", &x1, &x2);
			Lx = abs(float(x2) - float(x1));

			fgets(pipeString, 500, fileI);
			sscanf(pipeString, "%f %f", &x1, &x2);
			Ly = abs(float(x2) - float(x1));

			fgets(pipeString, 500, fileI);
			sscanf(pipeString, "%f %f", &x1, &x2);
			Lz = abs(float(x2) - float(x1));

			if(Lz < 0.1)
			{
				Lz = 0.0;
				dimension = 2;
			}
			else
				dimension = 3;

			fgets(pipeString, 500, fileI);

			if(strcmp(pipeString, "ITEM: ATOMS id element x y ix iy\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s %%%s %%%s\n", "d", "c", "f", "f", "d", "d");
				line_fmt_mode = 1;
			}
			else if(strcmp(pipeString, "ITEM: ATOMS id element x y vx vy ix iy\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s %%%s %%%s %%%s %%%s\n", "d", "c", "f", "f", "f", "f", "d", "d");
				line_fmt_mode = 2;
			}
			else if(strcmp(pipeString, "ITEM: ATOMS id element x y vx vy fx fy ix iy\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s %%%s %%%s %%%s %%%s %%%s %%%s\n", "d", "c", "f", "f", "f", "f", "f", "f", "d", "d");
				line_fmt_mode = 3;
			}
			else if(strcmp(pipeString, "ITEM: ATOMS id mol element x y z ix iy iz\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s %%%s %%%s %%%s %%%s %%%s\n", "d", "d", "c", "f", "f", "f", "d", "d", "d");
				line_fmt_mode = 4;

				int mol_nr;
				for(int i = 0; i < nAtoms; i++)
				{
					fgets(pipeString, 500, fileI);
					sscanf(pipeString, "%*d %d", &mol_nr);

					if(mol_nr > nMols) nMols = mol_nr;
				}
			}
			else if(strcmp(pipeString, "ITEM: ATOMS id element vx vy vz\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s %%%s\n", "d", "c", "f", "f", "f");
				line_fmt_mode = 5;
			}
			else if(strcmp(pipeString, "ITEM: ATOMS id element x y\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s\n", "d", "c", "f", "f");
				line_fmt_mode = 6;
			}
			else if(strcmp(pipeString, "ITEM: ATOMS id element vx vy\n") == 0)
			{
				sprintf(line_fmt, "%%%s %%%s %%%s %%%s\n", "d", "c", "f", "f");
				line_fmt_mode = 7;
			}
			else
			{
				printf("Invalid line format in trajectory %s. Exiting...\n", fpathI);
				exit(-1);
			}
		}

		rewind(fileI);
		frame_nr = -1;

		printf("\nInput trajectory %s file opened and ready to be read...\n", fpathI);
		printf("\nTrajectory Stats:\n");
		printf("\nTotal number of atoms: %d", nAtoms);

		if(nMols > 0)
			printf("\nTotal number of molecules: %d", nMols);

		if(strcmp(format, "cfg") == 0)
		{
			if(Lz == 0.0)
				printf("\n2D simulation box | Lx = %.2f | Ly = %.2f\n", Lx, Ly);
			else
				printf("\n3D simulation box | Lx = %.2f | Ly = %.2f | Lz = %.2f\n", Lx, Ly, Lz);
		}
	}

	if(count == true)
		countFrames();
}

void analysis::Trajectory::createOutputFile(char line[], bool newfile)
{
	if(newfile == true)
		remove(fpathO);

	fileO = fopen(fpathO, "a+");
	if(fileO == NULL)
	{
		printf("Cannot create new file or open existing file %s for writing. Exiting...\n", fpathO);
		exit(-1);
	}
	if(strcmp(line, "") != 0)
		fprintf(fileO, "%s", line);		
}

void analysis::Trajectory::closeTrajectory(bool closeI, bool closeO)
{
	if(closeI == true)
	{
		printf("\nClosing input trajectory file.\n");
		fclose(fileI);
	}

	if(closeO == true)
	{
		if(fileO != NULL)
		{
			printf("\nClosing output file.\n");
			fclose(fileO);	
		}
	}
}

void analysis::Trajectory::loadTrajectory(atom_style **ATOMS, System *BOX, int frameStart, int frameEnd, bool unwrap_pbc)
{
	atom_style *tempATOMS = new atom_style[nAtoms];

	frame_nr = -1;
	while( !feof(fileI) )
	{
		readThisFrame(tempATOMS);

		if(frame_nr >= frameStart and frame_nr <= frameEnd)
		{
			int new_frame = frame_nr - frameStart;

			copyThisFrame(tempATOMS, ATOMS[new_frame], unwrap_pbc);

			if(new_frame == 0) 
				printf("\nFirst frame: Step %ld\n", step);

			if(new_frame > 0 and unwrap_pbc == true)
			{
				for(int i = 0; i < nAtoms; i++)
				{
					ATOMS[new_frame][i].jumpx = ATOMS[new_frame - 1][i].jumpx;
					ATOMS[new_frame][i].jumpy = ATOMS[new_frame - 1][i].jumpy;

					float dx = ATOMS[new_frame][i].rxt1 - ATOMS[new_frame - 1][i].rxt1;
					float dy = ATOMS[new_frame][i].ryt1 - ATOMS[new_frame - 1][i].ryt1;

					if(dx <= -0.5*BOX->Lx) ATOMS[new_frame][i].jumpx++;
					else if(dx >= 0.5*BOX->Lx) ATOMS[new_frame][i].jumpx--;

					if(dy <= -0.5*BOX->Ly) ATOMS[new_frame][i].jumpy++;
					else if(dy >= 0.5*BOX->Ly) ATOMS[new_frame][i].jumpy--;
				}	
			}

			if(frame_nr == frameEnd)
			{
				printf("Last frame : Step %ld\n", step);
				break;
			} 
		}
	}

	printf("\nCoordinates loaded successfully for %d frames!\n", frame_nr - frameStart + 1);
	rewind(fileI);

	delete[] tempATOMS;
}

void analysis::Trajectory::copyThisFrame(atom_style *From, atom_style *To, bool unwrap_pbc)
{
	for(int i = 0; i < nAtoms; i++)
	{
		To[i].atom_id = From[i].atom_id;
		To[i].type = From[i].type;
		To[i].element = From[i].element;
		To[i].rxt1 = From[i].rxt1;
		To[i].ryt1 = From[i].ryt1;
		To[i].vx = From[i].vx;
		To[i].vy = From[i].vy;
		
		if(dimension == 3)
		{
			To[i].rzt1 = From[i].rzt1;
			To[i].vz = From[i].vz;
		}

		if(unwrap_pbc == false)
		{
			To[i].jumpx = From[i].jumpx;
			To[i].jumpy = From[i].jumpy;
			if(dimension == 3) 
				To[i].jumpz = From[i].jumpz;
		}
	}
}

void analysis::Trajectory::countFrames()
{
	atom_style *ATOMS = new atom_style [nAtoms];

	frame_nr = 0;
	while( !(feof(fileI)) )
		readThisFrame(ATOMS);

	totalFrames = frame_nr - 1;

	printf("Counting %d frames, %d atoms\n", totalFrames, nAtoms);

	frame_nr = -1;

	rewind(fileI);
	delete[] ATOMS;
}

void analysis::Trajectory::readThisFrame(atom_style *ATOMS)
{
	if(strcmp(format, "xyz") == 0)
	{
		for(int i = 0; i < 2; i++)
			fgets(pipeString, 500, fileI);

		sscanf(pipeString, "%*s %*s %ld", &step);

		for(int i = 0; i < nAtoms; i++)
		{
			fgets(pipeString, 500, fileI);
			sscanf(pipeString, "%c %f %f %*f %*f %*f", &ATOMS[i].element, &ATOMS[i].rxt1, &ATOMS[i].ryt1);
		}		
	}

	else if(strcmp(format, "cfg") == 0)
	{
		int pid, temp_id;

		for(int i = 0; i < 2; i++)
			fgets(pipeString, 500, fileI);

		sscanf(pipeString, "%ld", &step);

		for(int i = 0; i < 7; i++)
			fgets(pipeString, 500, fileI);

		for(int i = 0; i < nAtoms; i++)
		{
			fgets(pipeString, 500, fileI);

			sscanf(pipeString, "%d", &pid);
			pid -= 1;

			if(line_fmt_mode == 1)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].element, &ATOMS[pid].rxt1, &ATOMS[pid].ryt1, &ATOMS[pid].jumpx, &ATOMS[pid].jumpy);

			else if(line_fmt_mode == 2)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].element, &ATOMS[pid].rxt1, &ATOMS[pid].ryt1, &ATOMS[pid].vx, &ATOMS[pid].vy, &ATOMS[pid].jumpx, &ATOMS[pid].jumpy);

			else if(line_fmt_mode == 3)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].element, &ATOMS[pid].rxt1, &ATOMS[pid].ryt1, &ATOMS[pid].vx, &ATOMS[pid].vy, &ATOMS[pid].fx, &ATOMS[pid].fy, &ATOMS[pid].jumpx, &ATOMS[pid].jumpy);

			else if(line_fmt_mode == 4)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].mol_id, &ATOMS[pid].element, &ATOMS[pid].rxt1, &ATOMS[pid].ryt1, &ATOMS[pid].rzt1, &ATOMS[pid].jumpx, &ATOMS[pid].jumpy, &ATOMS[pid].jumpz);

			else if(line_fmt_mode == 5)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].element, &ATOMS[pid].vx, &ATOMS[pid].vy, &ATOMS[pid].vz);

			else if(line_fmt_mode == 6)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].element, &ATOMS[pid].rxt1, &ATOMS[pid].ryt1);

			else if(line_fmt_mode == 7)
				sscanf(pipeString, line_fmt, &ATOMS[pid].atom_id, &ATOMS[pid].element, &ATOMS[pid].vx, &ATOMS[pid].vy);

			if(ATOMS[pid].element == 'O')
			{
				ATOMS[pid].si = +1;
				ATOMS[pid].type = 1;
			}
			else if(ATOMS[pid].element == 'N')
			{
				ATOMS[pid].si = -1;
				ATOMS[pid].type = 2;
			}
		}
	}

	frame_nr++;
}

void analysis::Trajectory::writeThisFrame(atom_style *ATOMS, System *BOX, long add_step)
{
	if(strcmp(format, "xyz") == 0)
	{
		fprintf(fileO, "%d\n", nAtoms);
		fprintf(fileO, " Atoms. Timestep: %ld\n", step + add_step);

		for(int i = 0; i < nAtoms; i++)
			fprintf(fileO, "%c %g %g 0\n", ATOMS[i].element, ATOMS[i].rxt1, ATOMS[i].ryt1);
	}

	else if(strcmp(format, "cfg") == 0)
	{
		fprintf(fileO, "ITEM: TIMESTEP\n");
		fprintf(fileO, "%ld\n", step + add_step);
		fprintf(fileO, "ITEM: NUMBER OF ATOMS\n");
		fprintf(fileO, "%d\n", nAtoms);
		fprintf(fileO, "ITEM: BOX BOUNDS pp pp pp\n");
		fprintf(fileO, "%g %g\n", 0.0, BOX->Lx);
		fprintf(fileO, "%g %g\n", 0.0, BOX->Ly);
		fprintf(fileO, "0 0\n");

		if(line_fmt_mode == 1)
		{
			fprintf(fileO, "ITEM: ATOMS id element x y ix iy\n");

			for(int i = 0; i < nAtoms; i++)
				fprintf(fileO, line_fmt, ATOMS[i].atom_id, ATOMS[i].element, ATOMS[i].rxt1, ATOMS[i].ryt1, ATOMS[i].jumpx, ATOMS[i].jumpy);	
		}

		else if(line_fmt_mode == 2)
		{
			fprintf(fileO, "ITEM: ATOMS id element x y vx vy ix iy\n");

			for(int i = 0; i < nAtoms; i++)
				fprintf(fileO, line_fmt, ATOMS[i].atom_id, ATOMS[i].element, ATOMS[i].rxt1, ATOMS[i].ryt1, ATOMS[i].vx, ATOMS[i].vy, ATOMS[i].jumpx, ATOMS[i].jumpy);	
		}

		else if(line_fmt_mode == 3)
		{
			fprintf(fileO, "ITEM: ATOMS id element x y vx vy fx fy ix iy\n");

			for(int i = 0; i < nAtoms; i++)
				fprintf(fileO, line_fmt, ATOMS[i].atom_id, ATOMS[i].element, ATOMS[i].rxt1, ATOMS[i].ryt1, ATOMS[i].vx, ATOMS[i].vy, ATOMS[i].fx, ATOMS[i].fy, ATOMS[i].jumpx, ATOMS[i].jumpy);	
		}
	}
}

void analysis::Trajectory::computeCom(atom_style *ATOMS)
{
	xCom = yCom = 0.0;

	for(int i = 0; i < nAtoms; i++)
	{
		xCom += ATOMS[i].rxt1;
		yCom += ATOMS[i].ryt1;
	}	

	xCom /= nAtoms;
	yCom /= nAtoms;
}

long** analysis::Trajectory::sortAtomsByType(atom_style *ATOMS, System *BOX)
{
	if(nAtoms == 0)
	{
		printf("Error: Zero atoms present in trajectory! Exiting...\n");
		exit(-1);
	}
	if(BOX->nAtomTypes == 0)
	{
		printf("Error: Zero atom types present in trajectory! Exiting...\n");
		exit(-1);
	}

	long **group_lists = new long*[1 + BOX->nAtomTypes];
	int ctr[BOX->nAtomTypes];

	group_lists[0] = new long[BOX->nAtomTypes];
	for(int i = 0; i < BOX->nAtomTypes; i++)
	{
		group_lists[0][i] = 0;
		ctr[i] = 0;
	}

	// Counts the number of atoms belonging to a specific group
	for(int i = 0; i < nAtoms; i++)
	{
		if(ATOMS[i].type == 0 or ATOMS[i].type > BOX->nAtomTypes)
		{
			printf("Error: Invalid group type for Atom %d! Exiting ...\n", i + 1);
			exit(-1);
		}

		group_lists[0][ATOMS[i].type - 1]++;
	}

	// Creates lists to categorize atoms by type 
	for(int i = 0; i < BOX->nAtomTypes; i++)
	{
		group_lists[i + 1] = new long[group_lists[0][i]];
		for(int j = 0; j < group_lists[0][i]; j++)
			group_lists[i + 1][j] = 0;
	}

	// Sort atoms into their respective lists
	for(int i = 0; i < nAtoms; i++)
	{
		group_lists[ATOMS[i].type][ctr[ATOMS[i].type - 1]] = i + 1;
		ctr[ATOMS[i].type - 1]++;
	}

	return(group_lists);
}

int analysis::bounds(float a)
{
	float num = abs(a);
	int d = 0;
	while(int(num) > 0)
	{
		num = num/10;
		d += 1;
	}

	d = int(pow(10, d-1));

	int sign;
	if(a < 0) sign = -1;
	else sign = 1;

	int b = int(a) - sign*(int(sign*a)%d) + sign*d;
	return(b);
}