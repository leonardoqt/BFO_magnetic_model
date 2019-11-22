#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "afd_cell.h"

using namespace std;

void afd_cell :: init(int Nx, int Ny, int Nz, double Tild, double V0, double J0)
{
	nx=Nx; ny=Ny; nz=Nz;
	tild=Tild; v0=V0; j0=J0;
	num_tot = nx*ny*nz;
	theta.resize(num_tot);
	phi.resize(num_tot);
	spin.resize(num_tot);
}

void afd_cell :: rand_angle()
{
	int index;
	for(size_t t1=0; t1<nx; t1++)
	for(size_t t2=0; t2<ny; t2++)
	for(size_t t3=0; t3<nz; t3++)
	{
		index = (t1*ny+t2)*nz+t3;
		theta[index] = rand()/(double)RAND_MAX * M_PI;
		phi[index] = rand()/(double)RAND_MAX * 2 * M_PI;
	}
}

void afd_cell :: calculate_spin()
{
	int index;
	double tmp[3];
	for(size_t t1=0; t1<nx; t1++)
	for(size_t t2=0; t2<ny; t2++)
	for(size_t t3=0; t3<nz; t3++)
	{
		index = (t1*ny+t2)*nz+t3;
		tmp[0] = sin(theta[index])*cos(phi[index]);
		tmp[1] = sin(theta[index])*sin(phi[index]);
		tmp[2] = cos(theta[index]);
		spin[index] = tmp;
	}
}

void afd_cell :: send_param(vector <vector <double> >& Mc_param)
{
	Mc_param.resize(2);
	Mc_param[0].resize(num_tot);
	Mc_param[1].resize(num_tot);
	Mc_param[0] = theta;
	Mc_param[1] = phi;
}

void afd_cell :: receive_param(vector <vector <double> >& Mc_param)
{
	theta = Mc_param[0];
	phi = Mc_param[1];
}

double afd_cell :: ene()
{
	double energy;
	int index1,index2,i_sum;
	vec r1,r2,s1,s2;
	double tmp[3];

	energy = 0;
	for(size_t t1=0; t1<nx; t1++)
	for(size_t t2=0; t2<ny; t2++)
	for(size_t t3=0; t3<nz; t3++)
	{
		i_sum = (t1+t2+t3)%2;
		index1 = (t1*ny+t2)*nz+t3;
		s1 = spin[index1];
		// x direction
		index2 = (((t1+1)%nx)*ny+t2)*nz+t3;
		tmp[0] = 1; tmp[1] = tild*pow(-1,i_sum); tmp[2] = tild*pow(-1,i_sum);
		r1 = tmp;
		tmp[0] = -1;
		r2 = tmp;
		s2 = spin[index2];
		energy += v0*((r1^r2)*(s1^s2)) - j0*(s1*s2);
		// y direction
		index2 = (t1*ny+(t2+1)%ny)*nz+t3;
		tmp[1] = 1; tmp[0] = tild*pow(-1,i_sum); tmp[2] = tild*pow(-1,i_sum);
		r1 = tmp;
		tmp[1] = -1;
		r2 = tmp;
		s2 = spin[index2];
		energy += v0*((r1^r2)*(s1^s2)) - j0*(s1*s2);
		// z direction
		index2 = (t1*ny+t2)*nz+(t3+1)%nz;
		tmp[0] = 1; tmp[1] = tild*pow(-1,i_sum); tmp[2] = tild*pow(-1,i_sum);
		r1 = tmp;
		tmp[0] = -1;
		r2 = tmp;
		s2 = spin[index2];
		energy += v0*((r1^r2)*(s1^s2)) - j0*(s1*s2);
	}
	return energy;
}

void afd_cell :: print_spin()
{
	double index;
	for(size_t t1=0; t1<nx; t1++)
	for(size_t t2=0; t2<ny; t2++)
	for(size_t t3=0; t3<nz; t3++)
	{
		index = (t1*ny+t2)*nz+t3;
		cout<<t1<<'\t'<<t2<<'\t'<<t3<<'\t'<<spin[index]<<endl;
	}
}
