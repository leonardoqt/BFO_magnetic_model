#ifndef __AFD_CELL__
#define __AFD_CELL__

#include <vector>
#include <vec.h>

class afd_cell;

class afd_cell
{
private:
	int nx,ny,nz;
	int num_tot;
	std::vector < double > theta;
	std::vector < double > phi;
	std::vector < vec > spin;
public:
	double tild, v0, j0;
	void init(int Nx, int Ny, int Nz, double Tild, double V0, double J0);
	void rand_angle();
	void calculate_spin();
	void send_param(std::vector <std::vector <double> >& Mc_param);
	void receive_param(std::vector <std::vector <double> >& Mc_param);
	double ene();
	void print_spin();
};

#endif
