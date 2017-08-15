// Desc :
// Xinxin's article: https://zhuanlan.zhihu.com/p/26882619
#include <stdio.h>
#include <random>
#include <math.h>

using namespace std;

//define the mollify radius
#define EPS 0.01

//return a random number \in [a,b]
double frand(double a, double b)
{
	double  r = (double)rand()/(double)RAND_MAX;
	return (b-a)*r+a;
}

//a simple 2D vortexc struct
struct vortex2D
{
	double x; double y; double vort;
	vortex2D(){}
	vortex2D(const double& _x, const double& _y, const double& _vort)
	{
		x=_x; y=_y; vort=_vort;
	}
	vortex2D(const vortex2D &vortex)
	{
		x = vortex.x; y=vortex.y; vort=vortex.vort;
	}
};

//u component of velocity from a single vortex
double compute_u_from_single_vortex(double x, double y, vortex2D &vortex)
{
	double r_ij2 = (x-vortex.x)*(x-vortex.x) + (y-vortex.y)*(y-vortex.y);
	return vortex.vort*(vortex.y-y)/(r_ij2*3.1415926) * 0.5 *(1.0 - exp(-r_ij2/(EPS*EPS)));
}

//v component of velocity from a single vortex
double compute_v_from_single_vortex(double x, double y, vortex2D &vortex)
{
	double r_ij2 = (x-vortex.x)*(x-vortex.x) + (y-vortex.y)*(y-vortex.y);
	return vortex.vort*(x-vortex.x)/(r_ij2*3.1415926) * 0.5 *(1.0 - exp(-r_ij2/(EPS*EPS)));
}

//u component of velocity from all vortex influence
double compute_u_full_influence(double x, double y, std::vector<vortex2D> &vortex)
{
	double u=0;
	for (int i=0;i<vortex.size();i++)
	{
		u+=compute_u_from_single_vortex(x,y,vortex[i]);
	}
	return u;
}

//v component of velocity from all vortex influence
double compute_v_full_influence(double x, double y, std::vector<vortex2D> &vortex)
{
	double v=0;
	for (int i=0;i<vortex.size();i++)
	{
		v+=compute_v_from_single_vortex(x,y,vortex[i]);
	}
	return v;
}

//a simple rk3 integrator
void rk3_integrate_pos(double &x, double &y, std::vector<vortex2D> &vortex, double dt)
{
	double x0 = x;double y0 = y;
	double u0 = compute_u_full_influence(x0,y0,vortex);
	double v0 = compute_v_full_influence(x0,y0,vortex);
	double x1 = x0 + 0.5*dt*u0; double y1 = y0 + 0.5*dt*v0;

	double u1 = compute_u_full_influence(x1,y1,vortex);
	double v1 = compute_v_full_influence(x1,y1,vortex);
	double x2 = x0 + 0.75*dt*u1; double y2 = y0 + 0.75*dt*v1;

	double u2 = compute_u_full_influence(x2,y2,vortex);
	double v2 = compute_v_full_influence(x2,y2,vortex);
	x += 0.222222222222*dt*u0 + 0.333333333333*dt*u1 + 0.444444444444*dt*u2;
	y += 0.222222222222*dt*v0 + 0.333333333333*dt*v1 + 0.444444444444*dt*v2;
}

//write tracer particle data to a file so that can be visulized by matlab
void write_file(double *pos_x, double *pos_y, int num, int frame)
{
	char filename[256];
	sprintf(filename,"../../outputs/sims/Particle_data%04d.bin",frame);
	float *data;
	data = new float[num*4];
	for (int i=0;i<num;i++)
	{
		data[i*4+0] = pos_x[i];
		data[i*4+1] = pos_y[i];
		data[i*4+2] = 0;
		data[i*4+3] = 1;
	}
	FILE *f;
	f = fopen(filename,"wb");
	fwrite(&num, sizeof(int), 1, f);
	fwrite(data,sizeof(float),4*num,f);

	delete[]data;
	fclose(f);
}

// our main function computes two vortex ring leapforgging
int main(int argc, char * argv[])
{
	//initialize computation : 2 vortex ring pair and a lot tracer particles used to visualize

	//2 vortex ring pairs;
	std::vector<vortex2D> vortex_particles;
	vortex_particles.resize(0);
	vortex_particles.push_back(vortex2D(0.0,  1.0,  1.0));
	vortex_particles.push_back(vortex2D(0.0, -1.0, -1.0));
	vortex_particles.push_back(vortex2D(0.0,  0.3,  1.0));
	vortex_particles.push_back(vortex2D(0.0, -0.3, -1.0));


	//200000 tracer particles 
	int num_tracer = 200000;
	double *pos_x = new double[num_tracer];
	double *pos_y = new double[num_tracer];
	
	int num=0;
	while(num<num_tracer)
	{
		double x=frand(-0.5,0.5);
		double y=frand(-1.5,1.5);
		pos_x[num] = x;
		pos_y[num] = y;
		num++;
	}

	//our simulation
	double dt = 0.1;
	for (int T=0;T<300;T++)
	{
		write_file(pos_x,pos_y,num_tracer,T);

		//few substeps, not necessary
		for(int substep=0;substep<4;substep++)
		{
		//integrate tracers, we are going to use rk3 integrator
			#pragma omp parallel for
			for (int i=0;i<num_tracer;i++)
			{
				rk3_integrate_pos(pos_x[i],pos_y[i],vortex_particles, dt);
			}
		
			//integrate vortex particles
			std::vector<vortex2D> vortex_particles_temp;
			vortex_particles_temp.resize(vortex_particles.size());
			//before we integrate, we should copy vortex particles to
			//freeze them in space.
			for (int i=0;i<vortex_particles.size();i++)
			{
				vortex_particles_temp[i] = vortex_particles[i];
			}
			for(int i=0;i<vortex_particles_temp.size();i++)
			{
			
				double u=0;double v = 0;
				for (int j=0;j<vortex_particles.size();j++)
				{
					if(j!=i){

						u += compute_u_from_single_vortex(vortex_particles_temp[i].x,vortex_particles_temp[i].y,vortex_particles[j]);
						v += compute_v_from_single_vortex(vortex_particles_temp[i].x,vortex_particles_temp[i].y,vortex_particles[j]);
					}
				}
				vortex_particles_temp[i].x += dt*u;
				vortex_particles_temp[i].y += dt*v;
			}
			//copy back data
			for (int i=0;i<vortex_particles.size();i++)
			{
				vortex_particles[i] = vortex_particles_temp[i];
			}
		
		}
		
		printf("step %d done\n",T);
	}
	
	delete[]pos_x; delete[]pos_y;
	return 0;
}