#include <math.h>
#include <iostream>
#include <cstdlib>

// Numerical Integrator Function

using namespace std;

/*
N -> number of steps
F, G, H, I-> unkowns
z-> value to integrate over
params-> other numerical constants
-> lower boundary in 
*/

double add(double x, double y)
{
	return x+y;
}

double *Runge_Kutta_4ODE(int N, double z_0, double z_max, double F_0, double G_0, double H_0, double I_0, double (*Ffunc)(double,double,double,double,double,void*), double (*Gfunc)(double,double,double,double,double,void*), double (*Hfunc)(double,double,double,double,double,void*), double (*Ifunc)(double,double,double,double,double,void*), void *params)
{
	// STEP
	double h = (z_max-z_0)/N;
	double k_0, k_1, k_2, k_3;
	double l_0, l_1, l_2, l_3;
	double m_0, m_1, m_2, m_3;
	double n_0, n_1, n_2, n_3;

	double F_i = F_0;
	double G_i = G_0;
	double H_i = H_0;
	double I_i = I_0;

	double funcGrid[4];

	if (h > 0)
	{
		for (int i = 0; i<=N; ++i)
		{
			k_0 = h*(*Ffunc(z_0+h*i, F_i, G_i, H_i, I_i, params));
			l_0 = h*(*Gfunc(z_0+h*i, F_i, G_i, H_i, I_i, params));
			m_0 = h*(*Hfunc(z_0+h*i, F_i, G_i, H_i, I_i, params));
			n_0 = h*(*Ifunc(z_0+h*i, F_i, G_i, H_i, I_i, params));

			k_1 = h*(*Ffunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));
			l_1 = h*(*Gfunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));
			m_1 = h*(*Hfunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));
			n_1 = h*(*Ifunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));

			k_2 = h*(*Ffunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));
			l_2 = h*(*Gfunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));
			m_2 = h*(*Hfunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));
			n_2 = h*(*Ifunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));

			k_3 = h*(*Ffunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));
			l_3 = h*(*Gfunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));
			m_3 = h*(*Hfunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));
			n_3 = h*(*Ifunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));

			F_i = F_i + 1.0/6.*(k_0+2*k_1+2*k_2+k_3);
			G_i = G_i + 1.0/6.*(l_0+2*l_1+2*l_2+l_3);
			H_i = H_i + 1.0/6.*(m_0+2*m_1+2*m_2+m_3);
			I_i = I_i + 1.0/6.*(n_0+2*n_1+2*n_2+n_3);
		}
	}
	else
	{
		for (int i = 0; i<=N; --i)
		{
			k_0 = h*(*Ffunc(z_0+h*i, F_i, G_i, H_i, I_i, params));
			l_0 = h*(*Gfunc(z_0+h*i, F_i, G_i, H_i, I_i, params));
			m_0 = h*(*Hfunc(z_0+h*i, F_i, G_i, H_i, I_i, params));
			n_0 = h*(*Ifunc(z_0+h*i, F_i, G_i, H_i, I_i, params));

			k_1 = h*(*Ffunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));
			l_1 = h*(*Gfunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));
			m_1 = h*(*Hfunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));
			n_1 = h*(*Ifunc(z_0+h*i+0.5*h, F_i+0.5*k_0, G_i+0.5*l_0, H_i+0.5*m_0, I_i+0.5*n_0, params));

			k_2 = h*(*Ffunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));
			l_2 = h*(*Gfunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));
			m_2 = h*(*Hfunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));
			n_2 = h*(*Ifunc(z_0+h*i+0.5*h, F_i+0.5*k_1, G_i+0.5*l_1, H_i+0.5*m_1, I_i+0.5*n_1, params));

			k_3 = h*(*Ffunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));
			l_3 = h*(*Gfunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));
			m_3 = h*(*Hfunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));
			n_3 = h*(*Ifunc(z_0+h*i+h, F_i+k_2, G_i+l_2, H_i+m_2, I_i+n_2, params));

			F_i = F_i + 1.0/6.*(k_0+2*k_1+2*k_2+k_3);
			G_i = G_i + 1.0/6.*(l_0+2*l_1+2*l_2+l_3);
			H_i = H_i + 1.0/6.*(m_0+2*m_1+2*m_2+m_3);
			I_i = I_i + 1.0/6.*(n_0+2*n_1+2*n_2+n_3);
		}
	}

	funcGrid[0] = F_i;
	funcGrid[1] = G_i;
	funcGrid[2] = H_i;
	funcGrid[3] = I_i;

	return funcGrid;
} 
