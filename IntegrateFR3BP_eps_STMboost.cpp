// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

//#include <stdafx.h>
#define _CRT_SECURE_NO_WARNINGS
//#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <mex.h>
#include <matrix.h>
#include <math.h>


#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif


using namespace boost::numeric::odeint;
typedef std::vector< double > state_type;
using std::vector;

// CR3BP Equations of Motion
class eomCR3BP {
	double nu;
	double m;
	double e;
	double phase;

public:
	eomCR3BP(double nu_i,double m_i,double e_i,double phase_i){
		nu = nu_i;
		m = m_i;
		e = e_i;
		phase = phase_i;
	}

	void operator() (const state_type &x, state_type &dxdt, const double t)
	{
		double a1,a_1,a2,a_2,a3,a_3,xi,eeta,m2_a03,invR1_v3,invR_v3,invR1_v5,invR_v5;
		a1 = 0.1875*pow(m,2) + 0.5*pow(m,3) + pow(m,4)*7/12 + pow(m,5)*11/36 - pow(m,6)*30749/110592;
		a_1 = -pow(m,2)*19/16 -pow(m,3)*5/3 -pow(m,4)*43/36 -pow(m,5)*14/27 -pow(m,6)*7381/82944;
		a2 = pow(m,4)*25/256 + pow(m,5)*803/1920 + pow(m,6)*6109/7200;
		a_2 = pow(m,5)*23/640 + pow(m,6)*299/2400;
		a3 = pow(m,6)*833/12288;
		a_3 = pow(m,6)/192;
		// xi,eeta relative position of moon wrt earth in rotating frame 2(1+1/m), gives same result as variation_orbit_rot_2(m,tau)
		xi = ((a1+a_1)*cos(2*(t+phase))+(a2+a_2)*cos(4*(t+phase))+(a3+a_3)*cos(6*(t+phase)));
		eeta = ((a1-a_1)*sin(2*(t+phase))+(a2-a_2)*sin(4*(t+phase))+(a3-a_3)*sin(6*(t+phase)));

		m2_a03 = 1/pow(1 - 2*m/3 + 7*pow(m,2)/18 - 4*pow(m,3)/81,3);
		invR1_v3 = 1/pow(sqrt(pow(x[0]+nu*(1+xi),2) + pow(x[1]+nu*eeta,2) + pow(x[2],2)),3);
		invR_v3 = 1/pow(sqrt(pow(x[0]-(1-nu)*(1+xi),2) + pow(x[1]-(1-nu)*eeta,2) + pow(x[2],2)),3);
		invR1_v5 = 1/pow(sqrt(pow(x[0]+nu*(1+xi),2) + pow(x[1]+nu*eeta,2) + pow(x[2],2)),5);
		invR_v5 = 1/pow(sqrt(pow(x[0]-(1-nu)*(1+xi),2) + pow(x[1]-(1-nu)*eeta,2) + pow(x[2],2)),5);

		double Vx,Vy,Vz,Vxx,Vxy,Vxz,Vyy,Vyz,Vzz;

		Vx = (1+2*m+m*m)*x[0] + e*(0.75*m*m*(2*x[0]*cos(2*t)-2*x[1]*sin(2*t)) + 0.5*m*m*x[0]) + m2_a03*(-(1-nu)*(x[0]+nu*(1+xi))*invR1_v3 -nu*(x[0]-(1-nu)*(1+xi))*invR_v3);
		Vy = (1+2*m+m*m)*x[1] + e*(0.75*m*m*(-2*x[1]*cos(2*t)-2*x[0]*sin(2*t))+ 0.5*m*m*x[1]) + m2_a03*(-(1-nu)*(x[1]+nu*eeta)*invR1_v3 -nu*(x[1]-(1-nu)*eeta)*invR_v3);
		Vz = -e*m*m*x[2] + m2_a03*(-(1-nu)*x[2]*invR1_v3 -nu*x[2]*invR_v3);
		
		Vxx = 1+2*m+m*m + e*(1.5*m*m*cos(2*t) + 0.5*m*m) + m2_a03*(3*(1-nu)*pow(x[0]+nu*(1+xi),2)*invR1_v5 + 3*nu*pow(x[0]-(1-nu)*(1+xi),2)*invR_v5 - (1-nu)*invR1_v3 - nu*invR_v3);
		Vxy = e*-1.5*m*m*sin(2*t) + m2_a03*(3*(1-nu)*(x[0]+nu*(1+xi))*(x[1]+nu*eeta)*invR1_v5 + 3*nu*(x[0]-(1-nu)*(1+xi))*(x[1]-(1-nu)*eeta)*invR_v5);
		Vxz = m2_a03*(3*(1-nu)*(x[0]+nu*(1+xi))*x[2]*invR1_v5 + 3*nu*(x[0]-(1-nu)*(1+xi))*x[2]*invR_v5);
		Vyy = 1+2*m+m*m + e*(-1.5*m*m*cos(2*t) + 0.5*m*m) + m2_a03*(3*(1-nu)*pow(x[1]+nu*eeta,2)*invR1_v5 + 3*nu*pow(x[1]-(1-nu)*eeta,2)*invR_v5 - (1-nu)*invR1_v3 - nu*invR_v3);
		Vyz = m2_a03*(3*(1-nu)*(x[1]+nu*eeta)*x[2]*invR1_v5 + 3*nu*(x[1]-(1-nu)*eeta)*x[2]*invR_v5);
		Vzz = -e*m*m + m2_a03*(3*(1-nu)*pow(x[2],2)*invR1_v5 + 3*nu*pow(x[2],2)*invR_v5 -(1-nu)*invR1_v3 - nu*invR_v3);

		double phi[6][6];
		double A[6][6];
		double phidot[6][6];

		
		dxdt[0] = x[3];
		dxdt[1] = x[4];
		dxdt[2] = x[5];
		dxdt[3] = 2*(1+m)*x[4] + Vx;
		dxdt[4] = -2*(1+m)*x[3] + Vy;
		dxdt[5] = Vz;
	//define elements of STM, phi
	for (int k=0; k<6; k++)
	{
		for (int j=0; j<6; j++)
		{
			phi[j][k] = x[6+(6*k)+j]; 
		}
	}
	//compute terms for A, phidot = A*phi
	A[0][0] = 0; A[0][1] = 0; A[0][2] = 0;
	A[1][0] = 0; A[1][1] = 0; A[1][2] = 0;
	A[2][0] = 0; A[2][1] = 0; A[2][2] = 0;
    
	A[0][3] = 1; A[0][4] = 0; A[0][5] = 0;
	A[1][3] = 0; A[1][4] = 1; A[1][5] = 0;
	A[2][3] = 0; A[2][4] = 0; A[2][5] = 1;
    
	A[3][0] = Vxx; A[3][1] = Vxy; A[3][2] = Vxz;
	A[4][0] = Vxy; A[4][1] = Vyy; A[4][2] = Vyz;
	A[5][0] = Vxz; A[5][1] = Vyz; A[5][2] = Vzz;
    
	A[3][3] = 0; A[3][4] = 2*(1+m); A[3][5] = 0;
	A[4][3] = -2*(1+m); A[4][4] = 0; A[4][5] = 0;
	A[5][3] = 0; A[5][4] = 0; A[5][5] = 0;

	for (int k=0; k<6; k++)
	{
		for (int j=0; j<6; j++)
		{
			phidot[k][j] = A[k][0]*phi[0][j]+A[k][1]*phi[1][j]+A[k][2]*phi[2][j]+A[k][3]*phi[3][j]+A[k][4]*phi[4][j]+A[k][5]*phi[5][j];
		}
	}
    
	for (int k=0; k<6; k++)
	{
		for (int j=0; j<6; j++)
		{
			dxdt[6+(6*k)+j] = phidot[j][k];
		}
	}

	}
};

// Used to pull states out of the integration
struct getStateAndTime
{
	std::vector< state_type >& m_states;
	std::vector< double >& m_times;
	
	getStateAndTime(std::vector< state_type > &states, std::vector< double > &times) : m_states(states), m_times(times) {}

	void operator() (const state_type &x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}

};


// Propagation function similar to ATD that mexFunction uses to convert data to MATLAB readable
std::vector < vector < double >> prop(double ic[], double in_times[], double nu[], double m[],double e[],double phase[],int state_dim, int t_dim) {
	using namespace std;
	using namespace boost::numeric::odeint;

	typedef vector<double> state_type;
	vector < state_type > statesAndTimes;
	// Set vectors intermediate steps during integration
	vector<double> tOut;
	vector < state_type > statesOut;

	// Set initial integration condition for CR3BP
	state_type IC(42, 0);
	state_type t(t_dim, 0);
	for (int k=0; k<42; k++)
	{
		IC[k] = ic[k];
	}
	//IC[0] = ic[0];
	//IC[1] = ic[1];
	//IC[2] = ic[2];
	//IC[3] = ic[3];
	//IC[4] = ic[4];
	//IC[5] = ic[5];


	

	// Set integrator type -> Currently set at RK78
	typedef runge_kutta_fehlberg78<state_type> rk78;
	double relTol = 1e-12;
	double absTol = 1e-12;

	// Create eom to integrate
	eomCR3BP eom(nu[0],m[0],e[0],phase[0]);
	size_t steps;
	
	if (t_dim > 2) {
		statesOut.resize(t_dim);
		tOut.resize(t_dim);

		statesOut[0].resize(42);
		for (int k=0; k<42; k++)
		{
			statesOut[0][k] = IC[k];
		}

		for (int i = 1; i < t_dim; i++) {
			tOut[i-1] = in_times[i-1];
			tOut[i] = in_times[i];
			steps = integrate_adaptive(make_controlled <rk78>(absTol, relTol), eom, IC, tOut[i - 1], tOut[i], tOut[i] - tOut[i - 1]);
			statesOut[i].resize(42);
			for (int k=0; k<42; k++)
			{
				statesOut[i][k] = IC[k];
			}
		}
		steps = tOut.size()-1;
	}
	else {
		t[0] = in_times[0];
		t[1] = in_times[1];

		double h = t[1] > t[0] ? 1e-5 : -1e-5;

		// Integrate EOMs
		steps = integrate_adaptive(make_controlled <rk78>(absTol, relTol), eom, IC, t[0], t[1], h, getStateAndTime(statesOut, tOut));

	}

	// Determine step size (forward or backward) and set initial step size
	

	
	// Insert IC into list of state vectors
	statesAndTimes.resize(statesOut.size());
	// Put integrated states into state vectors
	for (int i = 0; i <= steps; i++) {
		statesAndTimes[i].resize(IC.size() + 1);
		for (int j = 0; j < statesAndTimes[i].size(); j++) {
			if (j == 0) {
				statesAndTimes[i][j] = tOut[i];
			}
			else {
				statesAndTimes[i][j] = statesOut[i][j-1];
			}
		}
	}

	return statesAndTimes;

}


// FOR TESTING PURPOSES
//int main() {
//
//	double *mu;
//	double mmu = 0.0122;
//
//	mu = &mmu;
//
//	double IC[6];
//	double t[2];
//	double stateDim = 6;
//	double tDim = 2;
//	vector<vector<double>> outTest;
//
//	IC[0] = -0.2700;
//	IC[1] = -0.4200;
//	IC[2] = 0.0;
//	IC[3] = 0.300;
//	IC[4] = -1.000;
//	IC[5] = 0.0;
//	t[0] = 0;
//	//t[1] = 0.01;
//	//t[2] = 0.02;
//	//t[3] = 0.03;
//	//t[4] = 0.04;
//	//t[5] = 0.05;
//	//t[6] = 0.06;
//	//t[7] = 0.07;
//	//t[8] = 0.08;
//	//t[9] = 0.09;
//	t[1] = 200;
//
//	outTest = prop(IC, t, mu, 6, 2);
//
//	double test = 0;
//
//}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	//declare variables
	mxArray *c_out_m, *t_out_m;
	const mwSize *dims, *time_dims;
	double *a, *c, *t, *times, *nu,*m,*e,*phase;
	int dim, time_dim;
	std::vector<vector<double>> state;

	//figure out dimensions
	dims = mxGetDimensions(prhs[0]); //get dimensions of state input vector
	time_dims = mxGetDimensions(prhs[1]); //get dimensions of time input vector

	if ((int)dims[0] != 1) {
		printf("input vector should be 1xn\n");
		return;
	}
	if ((int)time_dims[0] != 1) {
		printf("input vector should be 1xm\n");
		return;
	}

	dim = (int)dims[1]; //get dimension of state vector
	time_dim = (int)time_dims[1]; //get number of time inputs -> 2: [t0,tf], >2: [t0,t1,...,tf]

								  //associate pointers
	a = mxGetPr(prhs[0]);
	times = mxGetPr(prhs[1]);
	nu = mxGetPr(prhs[2]);
	m = mxGetPr(prhs[3]);
	e = mxGetPr(prhs[4]);
	phase = mxGetPr(prhs[5]);


	if (dim == 42) {
		state = prop(a, times, nu,m, e,phase,dim, time_dim);
	}
	else {
		printf("State vector should have 42 columns\n");
	}

	int nsteps = state.size();

	//associate outputs
	c_out_m = plhs[0] = mxCreateDoubleMatrix(nsteps, dim, mxREAL);
	t_out_m = plhs[1] = mxCreateDoubleMatrix(nsteps, 1, mxREAL);
	
	//associate pointers
	c = mxGetPr(c_out_m);
	t = mxGetPr(t_out_m);


	for (int k = 0; k<nsteps; k++) {
		for (int j = 0; j<=dim; j++) {
			if (j == 0) {
				t[k] = state[k][j];
			}
			else if (j>0) {
				c[nsteps*(j-1) + k] = state[k][j];
			}
		}
	}

	return;
}
