#include "stdafx.h"
#include <iostream>
#include <cmath>

#define ALBEDO .136		//Albedo of Lunar Surface
#define DX .0325		//Depth-direction Grid Element Size, m
#define DT 2360.448		//Time-direction Grid Element Size, m
#define EPSILON 0.92	//Emissivity of Lunar Regolith
#define ERRMAX 1.0E-6	//Upper Limit of Error
#define LA .8978		//Latitude Coefficient = cos(Lunar Local Latitude)
#define NX 1.3			//Constant-temperature Depth, m
#define NT 2360448.0	//Lunar Period 27.32 days = 2360448s
#define PI 3.14159
#define QE .12668		//Earth Radiation, W/m2
#define S 1353.			//Solar Constant, W/m2
#define SIGMA 5.67E-8	//Stefan-Boltzmann Constant
#define U0 200.			//Initial Temperature, K

double u[101][2001];
int n, m;
double fomax;			//Fourier Number

//Compute Specific Heat Capacity of Lunar Regolith
double CompC(double temp){
	const double d0 = -.05277;
	const double d1 = .15899E-2;
	const double d2 = -.03366E-4;
	const double d3 = .03142E-7;
	double ct = d0 + d1 * temp + d2 * pow(temp, 2) + d3 * pow(temp, 3);
	
	return(ct * 3600.);
}

//Compute Density of Lunar Regolith
double CompRho(double depth){
	return(1000. * 1.92 * ((depth * 100. + 12.2) / (depth * 100. + 18.)));
}

//Compute Thermal Conduvtivity of Lunar Regolith
double CompK(double temp, double dens){
	const double c0 = .40627783E-2;
	const double c1 = -.10295491E-4;
	const double c2 = .91660767E-8;
	const double c3 = -.22511580E-11;
	const double b0 = .30821872E-10;
	const double b1 = -.18565704E-13;
	const double b2 = -.12893852E-16;
	const double b3 = .17879447E-19;
	double k0 = c0 + c1 * dens + c2 * pow(dens, 2) + c3 * pow(dens, 3);
	double k3 = b0 + b1 * dens + b2 * pow(dens, 2) + b3 * pow(dens, 3);

	return(k0 + k3 * pow(temp, 3));
}

int main(){
	using namespace std;
	double a, mu, fo, uold;
	double err, errt, min, max;
	double c;					//Specific Heat Capacity, J/(kg.K)
	double rho;					//Density, kg/m3
	double k;					//Thermal Conduvtivity, W/(m.K)
	int times;					//Iteration Times

	n = (int)(NX / DX + .5) + 1;
	m = (int)(NT / DT + .5) + 1;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			u[i][j] = U0;
	times = 0;
	fomax = .0;
	err = 1.;
	while (err > ERRMAX){
		err = 0.;
		times++;
		for (int j = 0; j < m; j++){
			mu = S * sin(2 * PI * j * DT / NT) * LA;
			if (mu < 0) mu = 0;
			mu = (mu + QE) * (1. - ALBEDO) - EPSILON * SIGMA * pow(u[0][j], 4);
			int jj = (j == m - 1) ? 0 : j + 1;
			for (int i = 1; i < n; i++){
				uold = u[i][jj];
				c = CompC(uold);
				rho = CompRho(i * DX);
				k = CompK(uold, rho);
				a = k / rho / c;
				fo = a * DT / (DX * DX);
				if (fo > fomax) fomax = fo;
				if (i == n - 1)
					u[i][jj] = (1. - 2 * fo) * u[i][j] + 2 * fo * u[i - 1][j];
				else
					u[i][jj] = fo * (u[i + 1][j] + u[i - 1][j]) + (1. - 2 * fo) * u[i][j];
				errt = fabs(u[i][jj] - uold) / uold;
				if (errt > err) err = errt;
			}

			//Lunar surface where there is radiation input, is defined as heat flux boundary.
			uold = u[0][jj];
			c = CompC(uold);
			rho = CompRho(0.);
			k = CompK(uold, rho);
			a = k / rho / c;
			fo = a * DT / (DX * DX);
			if (fo > fomax) fomax = fo;
			u[0][jj] = (1. - 2 * fo) * u[0][j] + 2 * fo * u[1][j] + 2 * DT / (rho * c * DX) * mu;
			errt = fabs(u[0][jj] - uold) / uold;
			if (errt > err) err = errt;
		}
	}
	cout << "Iterate " << times << " Times, MaxFo = " << fomax << endl;
	min = 1000.;
	max = 0.;
	for (int j = 0; j < m; j++){
		if (u[0][j] < min) min = u[0][j];
		if (u[0][j] > max) max = u[0][j];
	}
	cout << "MinTemp = " << min << " K, MaxTemp = " << max << " K" << endl;

	return 0;
}
