#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

#pragma warning(disable : 4996)



double norm(double* x, int n)
{
	double res = 0;
	for (int i = 0; i < n; i++)
		res += fabs(x[i]);
	return res;
}

const int n = 12;

double F = 96483;//Faraday constant (C/mol)
double A = 1;//cross-sectional area (cm2)
double a = 490000;//specific area (cm2/cm3)

double I = 1e-3;//discharge current (A)
double tau = 1;//time step (s)
double epsilon = 1e-3;//error for Newton method
double Vt = 0.026;//thermal voltage (V)
double L = 0.01;//cathode thickness (cm)

double k1 = 1;
double K1 = 1.9e-5;
double V1 = 123.6;
double k2 = 5e8;
double K2 = 5e-15;
double V2 = 27.58;
double i01 = 1e-4 * a / F;
double i02 = 1e-5 * a / F;
double i03 = 1e-6 * a / F;
double i04 = 1e-7 * a / F;
double i05 = 1e-10 * a / F;
double i06 = 1e-3;
double U0 = 2.45;

double fSulfur = 0.1;//initial fraction of S8
double fLi2S = 1e-7;//initial fraction of Li2S
double fCarbon = 0.1;//initial fraction of carbon
double c_S80 = 1.9e-5;//reference concentration for S80
double c_S8 = 5e-7;//reference concentration for S8
double c_S6 = 1e-6;//reference concentration for S6
double c_S4 = 1e-8;//reference concentration for S4
double c_S2 = 1e-14;//reference concentration for S2
double c_S = 1e-17;//reference concentration for S
double c_Li = 1e-3;//reference concentration for Li

double* x_old, * x, * y, * rhs;
double** J;


void ComputeFunct(double* rhs)
{
	rhs[0] = (x[0] - x_old[0]) / tau + k1 * V1 * x[0] * (x[2] * c_S80 - K1);
	rhs[1] = (x[1] - x_old[1]) / tau - k2 * V2 * x[1] * (c_S * c_Li * c_Li * x[8] * x[8] - K2);
	rhs[2] = (c_S80 * (x[9] * x[2] - x_old[9] * x_old[2]) / tau) + k1 * x[0] * (c_S80 * x[2] - K1) - 0, 5 * i01 * a * (sqrt(x[3]) * exp(x[10]) - sqrt(x[2]) * exp(-x[10]));
	rhs[3] = (c_S8 * (x[9] * x[3] - x_old[9] * x_old[3]) / tau) + 0.5 * i01 * a / F * (sqrt(x[3]) * exp(x[10]) + sqrt(x[2]) * exp(-x[10]));
	rhs[4] = (c_S6 * (x[9] * x[4] - x_old[9] * x_old[4]) / tau) + 2 * i02 * a / F * (x[4] * x[4] * exp(x[10]) - pow(x[3], 1.5) * exp(-x[10])) - i03 * a / F * (pow(x[5], 1.5) * exp(x[10]) - x[4] * exp(-x[10]));
	rhs[5] = (c_S4 * (x[9] * x[5] - x_old[9] * x_old[5]) / tau) + 1.5 * i03 * a / F * (pow(x[5], 1.5) * exp(x[10]) - x[4] * exp(-x[10])) - 0.5 * i04 * a * (x[6] * exp(x[10]) - pow(x[5], 0.5) * exp(-x[10]));
	rhs[6] = (c_S2 * (x[9] * x[6] - x_old[9] * x_old[6]) / tau) + i04 * a / F * (x[6] * exp(x[10]) - pow(x[4], 0.5) * exp(-x[10])) - 0.5 * i05 * a / F * (x[7] * exp(x[10]) - pow(x[6], 0.5) * exp(-x[10]));
	rhs[7] = (c_S * (x[9] * x[7] - x_old[9] * x_old[7]) / tau) + i05 * a / F * (x[7] * exp(x[10]) + sqrt(x[6]) * exp(-x[10])) + k2 * x[1] * (c_S * c_Li * c_Li * x[7] * x[8] * x[8] - K2);
	rhs[8] = (c_Li * (x[9] * x[8] - x_old[9] * x_old[8]) / tau) + 2 * k2 * x[1] * (c_S * pow(c_Li, 2) * x[7] * x[8] * x[8] - K2) - i06 / (L * F) * (exp(x[11]) - x[8] * exp(-x[11]));
	rhs[9] = x[9] + x[0] + x[1] + fCarbon - 1;
	rhs[10] = I - i06 * A * (exp(x[11]) - x[8] * exp(-x[11]));
	rhs[11] = c_Li * (x[8] - 1) + 2 * c_S8 * (1 - x[3]) + 2 * c_S6 * (1 - x[4]) + 2 * c_S4 * (1 - x[5]) + 2 * c_S2 * (1 - x[6]);
}

void ComputeJacob(double** J)
{
	for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) J[i][j] = 0;

	J[0][0] = 1 / tau - k1 * V1 * c_S80 * x[3] + k1 * K1 * V1;
	J[0][2] = -K1 * V1 * c_S80 * x[0];
	J[1][1] = 1 / tau - k2 * V2 * c_S * pow(c_Li, 2) * x[7] * x[8] + k2 * V2 * K2;
	J[1][7] = -k2 * V2 * x[1] * c_S * pow(c_Li, 2) * x[8];
	J[1][8] = -k2 * V2 * x[1] * c_S * pow(c_Li, 2) * x[7] * 2 * x[8];

	J[2][0] = k1 * (x[2] * c_S80 - K1);
	J[2][2] = c_S80 * x[9] / tau + k1 * x[0] * c_S80 + 0, 25 * i01 * a / F * pow(x[2], 0.5) * exp(-x[10]);
	J[2][3] = -0.25 * i01 * a / F * pow(x[3], -0.5) * exp(-x[10]);
	J[2][9] = c_S80 * x[2] / tau;
	J[2][10] = -0.5 * i01 * a / F * (sqrt(x[3]) * exp(x[10]) + sqrt(x[2]) * exp(-x[10]));

	J[3][2] = -0.25 * i01 * a / F * pow(x[2], -0.5) * exp(-x[10]);
	J[3][3] = c_S8 * x[9] / tau + 0.25 * i01 * a / F * pow(x[3], -0.5) * exp(x[10]) + 2.25 * i02 * a / F * sqrt(x[3]) * exp(-x[10]);
	J[3][4] = -3 * i02 * a / F * x[4] * exp(x[10]);
	J[3][9] = c_S8 * x[3] / tau;
	J[3][10] = 0.5 * i01 * a / F * (sqrt(x[3]) * exp(x[10]) + sqrt(x[2]) * exp(-x[10]));

	J[4][3] = -3 * i02 * a / F * sqrt(x[3]) * exp(-x[10]);
	J[4][4] = (c_S8 * x[9] / tau) + (4 * i02 * a / F) * x[4] * exp(x[10]) + i03 * a / F * exp(-x[10]);
	J[4][5] = (-1.5 * i03 * a / F) * (sqrt(x[5]) * exp(x[10]));
	J[4][9] = c_S6 * x[4] / tau;
	J[4][10] = 2 * i02 * a / F * (x[4] * x[4] * exp(x[10]) + pow(x[3], 1.5) * exp(-x[10]));

	J[5][4] = -1.5 * i03 * a * exp(-x[10]) / F;
	J[5][5] = c_S * x[9] / tau + 2.25 * i03 * a / F * sqrt(x[5]) * exp(x[10]) + 0.25 * i04 * a / F * pow(x[5], -0.5) * exp(-x[10]);
	J[5][6] = -0.5 * i04 * a * exp(x[10]) / F;
	J[5][9] = c_S4 * x[5] / tau;
	J[5][10] = 1.5 * i03 * a / F * (pow(x[5], 1.5) * exp(x[10]) + x[4] * exp(-x[10])) - 0.5 * i04 * a * (x[6] * exp(x[10]) + pow(x[5], 0.5) * exp(-x[10])) / F;

	J[6][5] = -0.5 * i04 * a / F * pow(x[5], -0.5) * exp(-x[10]);
	J[6][6] = c_S2 * x[9] / tau + i04 * a / F * exp(x[10]) + 0.25 * i05 * a / F * (pow(x[6], -0.5) * exp(-x[10]));
	J[6][7] = -0.5 * i05 * a / F * exp(x[10]);
	J[6][9] = c_S2 * x[6] / tau;
	J[6][10] = i04 * a / F * (x[6] * exp(x[10]) + pow(x[4], 0.5) * exp(-x[10])) - 0.5 * i05 * a / F * (x[7] * exp(x[10]) + pow(x[6], 0.5) * exp(-x[10]));

	J[7][1] = k2 * c_S * pow(c_Li, 2) * x[7] * x[8] * x[8] - k2 * K2;
	J[7][6] = -0.5 * i05 * a / F * pow(x[6], -0.5) * exp(-x[10]);
	J[7][7] = c_S * x[9] / tau + i05 * a / F * exp(x[10]) + k2 * x[1] * c_S * c_Li * x[8] * x[8];
	J[7][8] = 2 * k2 * x[1] * c_S * pow(c_Li, 2) * x[7] * x[8];
	J[7][9] = c_S * x[7] / tau;
	J[7][10] = i05 * a / F * (x[7] * exp(x[10]) + sqrt(x[6]) * exp(-x[10]));

	J[8][1] = 2 * K2 * c_S * pow(c_Li, 2) * x[7] * x[8] - 2 * k2 * K2;
	J[8][7] = 2 * k2 * x[1] * c_S * pow(c_Li, 2) * x[8] * x[8];
	J[8][8] = c_Li * x[9] / tau + 4 * K2 * x[1] * c_S * pow(c_Li, 2) * x[7] * x[8] + (i06 / L * F) * exp(-x[11]);
	J[8][9] = c_Li * x[8] / tau;
	J[8][11] = (i06 / L * F) * exp(-x[11]) - (i06 / L * F) * x[8] * exp(-x[11]);

	J[9][0] = 1;
	J[9][1] = 1;
	J[9][9] = 1;

	J[10][8] = i06 * A * exp(-x[11]);
	J[10][11] = -i06 * A * exp(x[11]) + x[8] * exp(-x[11]);

	J[11][3] = -2 * c_S8;
	J[11][4] = -2 * c_S6;
	J[11][5] = -2 * c_S4;
	J[11][6] = -2 * c_S2;
	J[11][7] = -2 * c_S;
	J[11][8] = c_Li;
}

void Solve(double** A, double* rhs, double* X, int n)
{
	double** w = new double* [n];
	for (int i = 0; i < n; i++) w[i] = new double[n];
	int** index = new int* [n];
	for (int i = 0; i < n; i++) index[i] = new int[3];

	int nv = 1;
	int irow = 0;
	int icol = 0;
	double pivot;

	for (int i = 0; i < n; i++) {
		w[i][0] = rhs[i];
		index[i][2] = 0;
	}

	double determ = 1;
	for (int i = 0; i < n; i++) {
		double big = 0;
		for (int j = 0; j < n; j++) {
			if (index[j][2] != 1)
				for (int k = 0; k < n; k++)
				{
					if (index[k][2] > 1)
					{
						printf("Error at Gauss-Jordan: Matrix is singular\n");
						exit(1);
					}
					if (index[k][2] < 1)
						if (fabs(A[j][k]) > big)
						{
							irow = j;
							icol = k;
							big = fabs(A[j][k]);
						}
				}
		}

		index[icol][2] += 1;
		index[icol][0] = irow;
		index[icol][1] = icol;

		if (irow != icol)
		{
			determ = -determ;
			for (int l = 0; l < n; l++)
			{
				double temp = A[irow][l];
				A[irow][l] = A[icol][l];
				A[icol][l] = temp;
			}
			if (nv > 0)
				for (int l = 0; l < nv; l++)
				{
					double temp = w[irow][l];
					w[irow][l] = w[icol][l];
					w[icol][l] = temp;
				}
		}

		pivot = A[icol][icol];
		determ *= pivot;
		A[icol][icol] = 1;

		for (int l = 0; l < n; l++)
			A[icol][l] /= pivot;
		if (nv != 0)
			for (int l = 0; l < nv; l++)
				w[icol][l] /= pivot;
		for (int l1 = 0; l1 < n; l1++) {
			if (l1 != icol)
			{
				double t = A[l1][icol];
				A[l1][icol] = 0;
				for (int l = 0; l < n; l++)
					A[l1][l] -= A[icol][l] * t;
				if (nv != 0)
					for (int l = 0; l < nv; l++)
						w[l1][l] -= w[icol][l] * t;
			}
		}
	}


	for (int i = 0; i < n; i++) {
		int l = n - i - 1;
		if (index[l][0] != index[l][1]) {
			irow = index[l][0];
			icol = index[l][1];
			for (int k = 0; k < n; k++)
			{
				double temp = A[k][irow];
				A[k][irow] = A[k][icol];
				A[k][icol] = temp;
			}
		}
	}

	for (int k = 0; k < n; k++)
		if (index[k][2] != 1) {
			printf("Error at Gauss-Jordan: Matrix is singular\n");
			exit(1);
		}

	for (int i = 0; i < n; i++)
		X[i] = w[i][0];

	for (int i = 0; i < n; i++) delete[](w[i]);
	for (int i = 0; i < n; i++) delete[](index[i]);
	delete[]w;
	delete[]index;

}




void main()
{
	char temp;
	FILE* fp = fopen("voltage.csv", "w");//save the data in this file
	printf("SIMULATION STARTED...\n");

	//allocate memory
	J = new double* [n];
	for (int i = 0; i < n; i++) J[i] = new double[n];
	x = new double[n];
	y = new double[n];
	rhs = new double[n];
	x_old = new double[n];

	//initial conditions
	x[0] = fSulfur;// fraction of solid S8
	x[1] = fLi2S;// fraction of solid Li2S
	x[2] = 1;// relative concentration of S80
	x[3] = 1;// relative concentration of S8
	x[4] = 1;// relative concentration of S6
	x[5] = 1;// relative concentration of S4
	x[6] = 1;// relative concentration of S2
	x[7] = 1;// relative concentration of S
	x[8] = 1;// relative concentration of Li
	x[9] = 1 - fCarbon - fSulfur - fLi2S;// porosity
	x[10] = 0;// VC
	x[11] = 0;// VA


	double tmax = 3600 * 100;//10 hours
	for (int k = 1; k <= tmax / tau; k++)//time iterations
	{
		for (int j = 0; j < n; j++)
			x_old[j] = x[j];

		double time = k * tau;

		for (int i = 1; i <= 100; i++)//Newton iterations
		{
			ComputeFunct(rhs);
			ComputeJacob(J);
			Solve(J, rhs, y, n);//use the same solve function that you used in project 2
			double error = norm(y, n);

			for (int j = 0; j < n; j++)
				x[j] -= y[j];

			//printf("   error = %e\n", error);
			if (error < epsilon)
				break;

			if (i == 100)
			{
				printf("Newton did not converge after 100 interations. Make sure your Jacobian matrix is correctly set");
				scanf("%c",&temp); exit(1);
			}
		}

		double V = U0 - 2 * Vt * (x[11] - x[10]);
		printf("time = %.4f h\t Capacity = %.2f mAh/g\t Voltage = %.4f\n", time / 3600, I * time / 3.6 / (2.07 * fSulfur * L * A), V);
		fprintf(fp, "%e, %e\n", k * tau, V);

		if (V < 1.5) break;
	}

	fclose(fp);

	//free the memory allocated at the beginning
	for (int i = 0; i < n; i++) delete[] J[i];
	delete[]J;
	delete[]x;
	delete[]x_old;
	delete[]rhs;
	delete[]y;

	printf("...SIMULATION FINSIHED!\n");
	scanf("%c", &temp);

	
}

