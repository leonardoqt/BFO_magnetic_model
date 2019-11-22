#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <chrono>
#include <vec.h>
#include <mc.h>
#include "afd_cell.h"

using namespace std;

int main()
{
	chrono::high_resolution_clock::time_point now = chrono::high_resolution_clock::now();
	srand(now.time_since_epoch().count());
	ifstream input;

	// mc temperature
	double T_max, T_min;
	int period,repeat;
	thermo_profile TT;
	// mc control
	int num_kind, check_point;
	vector<int> num_param;
	vector<double> lambda;
	vector<double> max_param, min_param;
	vector<vector <double> > param,param_opt;
	mc mc_run;
	// cell 
	int nx,ny,nz;
	double tild,v0,j0;
	afd_cell cell;
	// other
	double ene_tmp, ene_opt, ene_new;

	// pass in parameters
	input.open("param.in");
	input>>T_min>>T_max;
	input>>period>>repeat>>check_point;
	input>>nx>>ny>>nz;
	input>>tild>>v0>>j0;
	input.close();
	// initialize objects
	TT.init(period,T_max,T_min);
	//
	num_kind = 2;
	num_param.resize(num_kind);
	num_param[0] = nx*ny*nz;
	num_param[1] = nx*ny*nz;
	lambda.resize(num_kind);
	lambda[0] = 1.5;
	lambda[1] = 3;
	max_param.resize(num_kind);
	max_param[0] = M_PI;
	max_param[1] = 2*M_PI;
	min_param.resize(num_kind);
	min_param[0] = 0;
	min_param[1] = 0;
	//
	cell.init(nx,ny,nz,tild,v0,j0);
	cell.rand_angle();
	cell.calculate_spin();
	ene_tmp = cell.ene();
	cell.send_param(param);
	mc_run.init(num_kind,check_point,num_param,lambda,max_param,min_param,param,ene_tmp);

	//run mc
	for(size_t t1=0; t1<repeat*period; t1++)
	{
		mc_run.gen_param_kind(param);
		cell.receive_param(param);
		cell.calculate_spin();
		ene_tmp = cell.ene();
		TT.gen_T(t1);
		mc_run.evaluate(TT,ene_tmp);
		mc_run.get_param_ene(param,param_opt,ene_tmp,ene_opt,ene_new);
		if ((t1+1)%period==0)
			cout<<ene_tmp<<'\t'<<ene_opt<<'\t'<<ene_new<<endl;
	}
	cell.receive_param(param_opt);
	cell.calculate_spin();
	cell.print_spin();

	return 0;
}
