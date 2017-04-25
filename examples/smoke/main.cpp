/*******************************************************************************
    VFXEPOCH - Physically based simulation VFX

    Copyright (c) 2016 Snow Tsui <trevor.miscellaneous@gmail.com>

    All rights reserved. Use of this source code is governed by
    the MIT license as written in the LICENSE file.
*******************************************************************************/
#include "Helpers.h"

VFXEpoch::Grid2DVector2DfField v0;
VFXEpoch::Grid2DVector2DfField v;
VFXEpoch::Grid2DVector2DfField grav;
VFXEpoch::Grid2DVector2DfField buoy;
VFXEpoch::Grid2DfScalarField d;
VFXEpoch::Grid2DfScalarField d0;
VFXEpoch::Grid2DfScalarField t;
VFXEpoch::Grid2DfScalarField t0;
VFXEpoch::Grid2DfScalarField pressure;
VFXEpoch::Grid2DfScalarField divergence;

// TODO: When apply IVOCK, use the following variables.
VFXEpoch::Grid2DfScalarField wn;
VFXEpoch::Grid2DfScalarField wBar;
VFXEpoch::Grid2DfScalarField wStar;
VFXEpoch::Grid2DfScalarField dw;
VFXEpoch::Grid2DfScalarField psi;
VFXEpoch::Grid2DVector2DfField dvel;
std::vector<VFXEpoch::Particle2D> particles;

VFXEpoch::Solvers::SL2D *sl2D_solver = NULL;
Helper::SimulationParameters simParams;
int ID = 0;
int width = 640;
int height = 720;
int mouse_status[3];
int mx0, my0, mx, my;
bool bVel = false;
bool bSmoke = false;
bool bParticles = true;
bool bPause = false;
int stopFrame = -1;
int frame_counter = 0;

void Init(int argc, char** argv);
void WindowShowup(int width, int height);
void PreDisplay();
void PostDisplay();
void Display();
void DisplayParticles();
void DisplayVelocityField();
void DispolayDensityField();
void GetUserOperations(VFXEpoch::Grid2DfScalarField& density, VFXEpoch::Grid2DVector2DfField& vel);
void Advance();
void IVOCKAdvance();
void ParticlesAdvector_RKII();
void Reset();
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void Reshape(int width, int height);
void Idle();
void Keys(unsigned char key, int x, int y);
void Loop();
void KeepSource();
void Close();

void
Init(int argc, char **argv)
{
	// Simulation parameters
	simParams.nx = 166;	simParams.ny = 128;
	simParams.dt = 0.01f;
	simParams.diff = 0.0f; simParams.visc = 0.0f;
	simParams.src = 100.0f;	simParams.src_rate = 1.0f;
	simParams.user_force = 0.0f;
	simParams.heat_source = 2000.0f;
	simParams.streamer_len = 5.0f;
	simParams.num_particles = 10000;
	simParams.vort_conf_eps = 0.95f;
	simParams.particle_life_span_rev = 0.9f;
	simParams.linear_solver_iterations = 30;
	// Simulation parameters End

	// Related field
	particles.resize(simParams.num_particles);
	v.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	v0.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	grav.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	buoy.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	d.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	d0.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	t.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	t0.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	wn.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	psi.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	dw.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	wBar.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	wStar.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	dvel.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	pressure.ResetDimension(simParams.nx + 2, simParams.ny + 2);
	divergence.ResetDimension(simParams.nx + 2, simParams.ny + 2);

	VFXEpoch::Zeros(v);
	VFXEpoch::Zeros(v0);
	VFXEpoch::Zeros(grav);
	VFXEpoch::Zeros(buoy);
	VFXEpoch::Zeros(d);
	VFXEpoch::Zeros(d0);
	VFXEpoch::Zeros(t);
	VFXEpoch::Zeros(t0);
	VFXEpoch::Zeros(wn);
	VFXEpoch::Zeros(wBar);
	VFXEpoch::Zeros(wStar);
	VFXEpoch::Zeros(dw);
	VFXEpoch::Zeros(dvel);
	VFXEpoch::Zeros(psi);
	VFXEpoch::Zeros(pressure);
	VFXEpoch::Zeros(divergence);

	sl2D_solver = new VFXEpoch::Solvers::SL2D();
	if(!sl2D_solver) exit(-1);

	// Setup particles
	float r(0.0f), g(0.0f), b(0.0f);
	float x(0.0f), y(0.0f);
	for (std::vector<VFXEpoch::Particle2D>::iterator ite = particles.begin(); ite != particles.end(); ite++){
		x = (simParams.nx / 2 + VFXEpoch::RandomI(-40, 40)) * 1.0f / simParams.nx;
		y = (VFXEpoch::RandomI(0, 30)) * 1.0f / simParams.nx;
		ite->pos = VFXEpoch::Vector2Df(x, y);
		ite->vel = VFXEpoch::Vector2Df(0.0f, 0.0f);
		ite->color = VFXEpoch::Vector3Df(0.0f, 0.0f, 0.0f);
	}

	// Initialize gas solver
	if (!sl2D_solver->Initialize(VFXEpoch::Vector2Di(simParams.nx + 2, simParams.nx + 2), VFXEpoch::Vector2Df(1.0f / simParams.nx, 1.0f / simParams.nx),
		simParams.linear_solver_iterations, simParams.dt, simParams.diff, simParams.visc, simParams.src_rate)) {
		cout << "Solver initialization failed" << endl;
		exit(-1);
	}
	else {
		cout << "Solver was initialized successfully" << endl;
	}

	// If anx field get new value, it requires call the following
	// interfaces to transport data to the solver.
	sl2D_solver->SetField(v, VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL);
	sl2D_solver->SetField(v0, VFXEpoch::COMPUTATIONAL_VECTOR_FIELD_2D::VEL_PREV);
	sl2D_solver->SetField(d, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::DENSITY);
	sl2D_solver->SetField(d0, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::DENSITY_PREV);
	sl2D_solver->SetField(t, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::TEMPERATURE);
	sl2D_solver->SetField(t0, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::TEMPERATURE_PREV);
	sl2D_solver->SetField(pressure, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::PRESSURE);
	sl2D_solver->SetField(divergence, VFXEpoch::COMPUTATIONAL_SCALAR_FIELD_2D::DIVERGENCE);
	sl2D_solver->SetFieldBoundary(VFXEpoch::BOUNDARY::NEUMANN_OPEN, VFXEpoch::EDGES_2DSIM::TOP);
	sl2D_solver->SetFieldBoundary(VFXEpoch::BOUNDARY::NEUMANN_OPEN, VFXEpoch::EDGES_2DSIM::BOTTOM);
	sl2D_solver->SetFieldBoundary(VFXEpoch::BOUNDARY::NEUMANN_OPEN, VFXEpoch::EDGES_2DSIM::RIGHT);
	sl2D_solver->SetFieldBoundary(VFXEpoch::BOUNDARY::NEUMANN_OPEN, VFXEpoch::EDGES_2DSIM::LEFT);

	// Set gravity
	for (int i = 1; i != grav.getDimY() - 1; i++) {
		for (int j = 1; j != grav.getDimX() - 1; j++) {
			grav.setData(VFXEpoch::Vector2Df(0.0f, -9.8f), i, j);
		}
	}

	glutInit(&argc, argv);
}

void
WindowShowup(int width, int height)
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition((glutGet(GLUT_SCREEN_WIDTH) - width) / 2,
		(glutGet(GLUT_SCREEN_HEIGHT) - height) / 2);
	glutInitWindowSize(width, height);
	ID = glutCreateWindow("Smoke Sim VFXEpoch");
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	PreDisplay();
	glutKeyboardFunc(Keys);
	glutMouseFunc(mouse_func);
	glutMotionFunc(motion_func);
	glutReshapeFunc(Reshape);
	glutIdleFunc(Idle);
	glutDisplayFunc(Display);
}

void
PreDisplay()
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0f, 1.0, 0.0f, 1.0f);
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glEnable(GL_POINT_SMOOTH);
}

void
PostDisplay()
{
	glutSwapBuffers();
}

void
Display()
{
	PreDisplay();

	if (bVel)
		DisplayVelocityField();
	else if (bSmoke)
		DispolayDensityField();
	else if (bParticles)
		DisplayParticles();
	else
		DisplayParticles();

	PostDisplay();
}

void
DisplayVelocityField()
{
	glColor3f(0.0f, 0.0f, 0.0f);
	glLineWidth(1.0f);
	float x, y;
	float dx = 1.0f / (simParams.nx), dy = 1.0f / (simParams.nx);
	float streamerxRate = simParams.streamer_len / simParams.nx, streameryRate = simParams.streamer_len / simParams.nx;

	glBegin(GL_LINES);
	for (int i = 1; i <= simParams.nx + 1; i++){
		for (int j = 1; j <= simParams.nx + 1; j++){
			x = (i - 0.5f) * dx;
			y = (j - 0.5f) * dy;
			glVertex2f(x, y);
			glVertex2f(x + v(i, j).m_x * streamerxRate, y + v(i, j).m_y * streameryRate);
		}
	}
	glEnd();

	//glColor3f(0.2f, 0.6f, 1.0f);
	//glPointSize(1.0f);
	//glBegin(GL_POINTS);
	//for (int i = 1; i <= simParams.nx + 1; i++){
	//	for (int j = 1; j <= simParams.nx + 1; j++){
	//		x = (j - 0.5f) * dx;
	//		y = (i - 0.5f) * dy;
	//		glVertex2f(x, y);
	//	}
	//}
	//glEnd();
}

void
DisplayParticles()
{
	glPointSize(1.0f);

	glBegin(GL_POINTS);
	for (std::vector<VFXEpoch::Particle2D>::iterator ite = particles.begin(); ite != particles.end(); ite++){
		if (ite->pos.m_x < 0 || ite->pos.m_x > 1 ||
			ite->pos.m_y < 0 || ite->pos.m_y > 1) {
			ite->pos.m_x = (simParams.nx / 2 + VFXEpoch::RandomI(-40, 40)) * 1.0f / simParams.nx;
			ite->pos.m_y = (VFXEpoch::RandomI(0, 30)) * 1.0f / simParams.nx;
		}
		glColor3f(ite->color.m_x, ite->color.m_y, ite->color.m_z);
		glVertex2f(ite->pos.m_x, ite->pos.m_y);
	}
	glEnd();
}

void
DispolayDensityField()
{
	float x, y, hx, hy, d00, d01, d10, d11;
	hx = 1.0f / simParams.nx;
	hy = 1.0f / simParams.nx;
	glBegin(GL_QUADS);
	for (int i = 0; i <= simParams.nx; i++)	{
		x = (i - 0.5f) * hx;
		for (int j = 0; j <= simParams.nx; j++)	{
			y = (j - 0.5f) * hy;
			d00 = d(i, j);
			d01 = d(i, j + 1);
			d10 = d(i + 1, j);
			d11 = d(i + 1, j + 1);

			glColor3f(1.f - d00, 1.f - d00, 1.f - d00); glVertex2f(x, y);
			glColor3f(1.f - d10, 1.f - d10, 1.f - d10); glVertex2f(x + hx, y);
			glColor3f(1.f - d11, 1.f - d11, 1.f - d11); glVertex2f(x + hx, y + hy);
			glColor3f(1.f - d01, 1.f - d01, 1.f - d01); glVertex2f(x, y + hy);
		}
	}
	glEnd();
}

void
GetUserOperations(VFXEpoch::Grid2DfScalarField& density, VFXEpoch::Grid2DVector2DfField& vel)
{
	VFXEpoch::Vector2Df vec(0.0f, 0.0f);
	VFXEpoch::Zeros(vel);
	VFXEpoch::Zeros(density);
	if (!mouse_status[0] && !mouse_status[2]) return;

	int i = (int)((mx / (float)width) * simParams.nx + 1);
	int j = (int)(((height - my) / (float)height) * simParams.nx + 1);
	if (i < 1 || i > simParams.nx || j < 1 || j > simParams.nx) return;
	if (mouse_status[0]) {
		//vec.m_x = simParams.user_force * (mx - mx0);
		//vec.m_y = simParams.user_force * (my0 - my);
		//vel.setData(vec, i, j);
	}
	if (mouse_status[2]) {
		density(i, j) = simParams.src;
	}

	mx0 = mx;
	my0 = my;
}

void
ParticlesAdvector_RKII()
{
	VFXEpoch::Grid2DfScalarField du(v.getDimY(), v.getDimX(), 1.0f / v.getDimX(), 1.0f / v.getDimY());
	VFXEpoch::Grid2DfScalarField dv(v.getDimY(), v.getDimX(), 1.0f / v.getDimX(), 1.0f / v.getDimY());
	VFXEpoch::ExtractComponents(du, v, VFXEpoch::VECTOR_COMPONENTS::X);
	VFXEpoch::ExtractComponents(dv, v, VFXEpoch::VECTOR_COMPONENTS::Y);
	sl2D_solver->_advect_particles_rk2(du, dv, particles);
	du.clear();
	dv.clear();
}

void
Advance()
{
	/*-------------------------------- Advect Velocity Field --------------------------------*/
	VFXEpoch::Zeros(wn);
	VFXEpoch::Zeros(wBar);
	VFXEpoch::Zeros(wStar);
	VFXEpoch::Zeros(psi);
	VFXEpoch::Zeros(dvel);
	VFXEpoch::Zeros(dw);
	VFXEpoch::Zeros(buoy);

	sl2D_solver->_set_source(v0, grav);
	sl2D_solver->_set_source(v, v0);

	ParticlesAdvector_RKII();

	VFXEpoch::Swap(v, v0);
	sl2D_solver->_diffuse(v, v0);
	sl2D_solver->_project(v, pressure, divergence);
	// VFXEpoch::Swap(v, v0);
	sl2D_solver->_advect(v, v0, v0);
	sl2D_solver->_get_buoyancy(d, t, buoy, 0.1f, 0.4f);
	v += buoy;
	//sl2D_solver->AddVortConf(v, simParams.vort_conf_eps, VFXEpoch::VORT_METHODS::LEAST_SQUARE);
	sl2D_solver->_project(v, pressure, divergence);

	///*------------------------------ Advect temperature Field ------------------------------*/
	sl2D_solver->_set_source(t, t0);
	VFXEpoch::Swap(t, t0);
	sl2D_solver->_diffuse(t, t0);
	VFXEpoch::Swap(t, t0);
	sl2D_solver->_advect(t, t0, v);


	/*-------------------------------- Advect density Field --------------------------------*/
	sl2D_solver->_set_source(d, d0);
	VFXEpoch::Swap(d, d0);
	sl2D_solver->_diffuse(d, d0);
	VFXEpoch::Swap(d, d0);
	sl2D_solver->_advect(d, d0, v);
}

void
IVOCKAdvance()
{
	/*-------------------------------- Advect Velocity Field --------------------------------*/
	VFXEpoch::Zeros(wn);
	VFXEpoch::Zeros(wBar);
	VFXEpoch::Zeros(wStar);
	VFXEpoch::Zeros(psi);
	VFXEpoch::Zeros(dvel);
	VFXEpoch::Zeros(dw);
	VFXEpoch::Zeros(buoy);

	sl2D_solver->_set_source(v0, grav);
	sl2D_solver->_set_source(v, v0);

	ParticlesAdvector_RKII();

	VFXEpoch::Swap(v, v0);
	sl2D_solver->_diffuse(v, v0);
	sl2D_solver->_project(v, pressure, divergence);
	VFXEpoch::Swap(v, v0);

	// Compute vorticity from the original velocity field
	// VFXEpoch::Analysis::computeCurl_uniform_Richardson(wn, v0);
	VFXEpoch::Analysis::computeCurl_uniform_LS(wn, v0);
	// VFXEpoch::Analysis::computeCurl_uniform_Stokes(wn, v0);

	// Advect vorticity by following velocity field
	sl2D_solver->_advect(wBar, wn, v0);

	sl2D_solver->_advect(v, v0, v0);

	// Compute vorticity from the velocity field which has been advected
	// VFXEpoch::Analysis::computeCurl_uniform_Richardson(wStar, v);
	VFXEpoch::Analysis::computeCurl_uniform_LS(wStar, v);
	// VFXEpoch::Analysis::computeCurl_uniform_Stokes(wStar, v);

	// Get the difference from 2 vorticity field and scaled by -1.0f
	dw = wBar - wStar;
	dw.scale(-1.0f);

	// Linearly solve the psi and deduce velocity from vorticity
	//VFXEpoch::LinearSolver::GSSolve(psi, dw, sl2D_solver->getFieldBoundaries(), 1, 4, simParams.linear_solver_iterations);
	VFXEpoch::LinearSolver::JacobiSolve(psi, dw, sl2D_solver->getFieldBoundaries(), 1, 4, simParams.linear_solver_iterations);
	//VFXEpoch::LinearSolver::MultigridSolve_V_Cycle(1.f / (simParams.nx - 1), psi, dw, sl2D_solver->getFieldBoundaries(), 1, 4, 30);
	VFXEpoch::Analysis::find_vector_from_vector_potential_2D(dvel, psi);
	// Combine the differences
	v += dvel;

	// Get buoyancy
	sl2D_solver->_get_buoyancy(d, t, buoy, 0.1f, 0.4f);
	v += buoy;
	//sl2D_solver->AddVortConf(v, simParams.vort_conf_eps, VFXEpoch::VORT_METHODS::LEAST_SQUARE);
	sl2D_solver->_project(v, pressure, divergence);

	/*------------------------------ Advect temperature Field ------------------------------*/
	sl2D_solver->_set_source(t, t0);
	VFXEpoch::Swap(t, t0);
	sl2D_solver->_diffuse(t, t0);
	VFXEpoch::Swap(t, t0);
	sl2D_solver->_advect(t, t0, v);


	/*-------------------------------- Advect density Field --------------------------------*/
	sl2D_solver->_set_source(d, d0);
	VFXEpoch::Swap(d, d0);
	sl2D_solver->_diffuse(d, d0);
	VFXEpoch::Swap(d, d0);
	sl2D_solver->_advect(d, d0, v);
}

void
Reset()
{
	float r(0.0f), g(0.0f), b(0.0f);
	float x(0.0f), y(0.0f);
	for (std::vector<VFXEpoch::Particle2D>::iterator ite = particles.begin(); ite != particles.end(); ite++){
		x = (simParams.nx / 2 + VFXEpoch::RandomI(-40, 40)) * 1.0f / simParams.nx;
		y = (VFXEpoch::RandomI(0, 30)) * 1.0f / simParams.nx;
		ite->pos = VFXEpoch::Vector2Df(x, y);
		ite->vel = VFXEpoch::Vector2Df(0.0f, 0.0f);
		ite->color = VFXEpoch::Vector3Df(0.0f, 0.0f, 0.0f);
	}

	v.zeroVectors(); v0.zeroVectors();
	d.zeroScalars(); d0.zeroScalars();
	t.zeroScalars(); t0.zeroScalars();
	pressure.zeroScalars();	divergence.zeroScalars();
	sl2D_solver->Reset();

	frame_counter = 0;
}

void mouse_func(int button, int state, int x, int y)
{
	mx0 = mx = x;
	my0 = my = y;

	mouse_status[button] = state == GLUT_DOWN;
}

void
motion_func(int x, int y)
{
	mx = x;
	my = y;
}

void
Reshape(int width, int height)
{
	glutSetWindow(ID);
	glutReshapeWindow(width, height);
	::width = width;
	::height = height;
}

void
Keys(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 's':
		bSmoke = !bSmoke;
		break;
	case 'v':
		bVel = !bVel;
		break;
	case 'p':
		bParticles = !bParticles;
		break;
	case 'x':
		exit(0);
		break;
	case 'c':
		Reset();
		break;
	case ' ':
		bPause = !bPause;
		break;
	default:
		break;
	}
}

void
Idle()
{
	GetUserOperations(d0, v0);
	if (!bPause)
	{
		if (frame_counter != stopFrame)
		{
			KeepSource();
			//Advance();
			IVOCKAdvance();
			frame_counter++;
		}
		else
		{
			Reset();
		}
	}
	glutSetWindow(ID);
	glutPostRedisplay();
}

void
Loop()
{
	glutMainLoop();
}

void
KeepSource()
{
	int idxi = simParams.nx / 2;
	int idxj = 10;
	v0(idxi, idxj).m_y	= simParams.user_force;
	d0(idxi, idxj) = simParams.src;
	t0(idxi, idxj) = simParams.heat_source;
	t0(idxi + 1, idxj) = simParams.heat_source;
	t0(idxi - 1, idxj) = simParams.heat_source;
	t0(idxi + 2, idxj) = simParams.heat_source;
	t0(idxi - 2, idxj) = simParams.heat_source;
}

void
Close()
{
	particles.clear();
	v.clear();	v0.clear();
	d.clear();	d0.clear();
	t.clear(); t0.clear();
	grav.clear(); buoy.clear();
	wn.clear(); wBar.clear();
	wStar.clear(); dvel.clear();
	psi.clear();
	pressure.clear(); divergence.clear();

	if(sl2D_solver){
		sl2D_solver->Shutdown();
		delete [] sl2D_solver;
	}
}

int
main(int argc, char** argv)
{
	cout << "-------------- Smoke Simulation Info --------------" << endl;
	cout << "Press 'v' switching to display velocity field" << endl;
	cout << "Press 's' switching to display smoke" << endl;
	cout << "Press 'c' to reset simulation" << endl;
	cout << "Press 'x' to exit" << endl;

	cout << endl << "-------------- Simulation Setup --------------" << endl;
	Init(argc, argv);
	WindowShowup(width, height);
	cout << endl << "Simulation parameters:";
	cout << endl << simParams << endl;

	cout << endl << "-------------- Simulation Start --------------" << endl;
	Loop();
	Close();

	return 0;
}
