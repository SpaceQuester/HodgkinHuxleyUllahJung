#define _CRT_SECURE_NO_WARNINGS

#define _USE_MATH_DEFINES
#define USE_OMP

#include "math.h"
#include <stdlib.h>
#include <stdio.h>
#include <locale.h>
#include <time.h>
#include <stdbool.h>
#include <list>
#include <vector>

#ifdef USE_OMP
#include <omp.h>
#endif

using namespace std;

#define Node_count 200 // 200 Default

int MaxDeep = 50; // 50 Default

#define Equations_per_node 12 // !!! 12 - Don't change !!!
#define Equations_count Node_count * Equations_per_node

double* f;
double* f_diff;

double** k;
double* phi_k1;
double* phi_k2;
double* phi_k3;

bool enable_I_syn_out = false;

bool write_file = false;
bool omp_enable = false;

double c_0 = 2; // uM
double c_1 = 0.185;
double v_1 = 6; // s^-1
double v_2 = 0.11; // s^-1
double v_3 = 2.2; // uM/s
double* v_4; // uM/s - Controling parameter // 0.5 //double v_4[Node_count]; // uM/s - Controling parameter //0.495
double v_5 = 0.025; // uM/s
double v_6 = 0.2;  // uM/s
double k_1 = 0.5; // s^-1
double k_2 = 1; // uM
double k_3 = 0.1;
double k_4 = 1.1; // uM/s
double a_2 = 0.14; // uM/s
double d_1 = 0.13; // uM
double d_2 = 1.049; // uM
double d_3 = 0.9434; // uM
double d_5 = 0.082; // uM
double alpha = 0.8;
double tau_IP3 = 7.143; // s
double IP3_star = 0.16; // uM
double d_Ca = 0.001; // 0.001
double d_IP3 = 0.2; // 0.12
double alpha_Glu = 2; // 2
double g_astro; // 3

// https://neuronaldynamics.epfl.ch/online/Ch2.S2.html
double C_m = 1; // muF/cm^2
double g_K = 35; // mS/cm^2
double g_Na = 40; // mS/cm^2
double g_L = 0.3; // mS/cm^2
double E_K = -77; // mV
double E_Na = 55; // mV
double E_L = -65; // mV

double C_m_P = 1; // muF/cm^2
double g_K_P = 35; // mS/cm^2
double g_Na_P = 40; // mS/cm^2
double g_L_P = 0.3; // mS/cm^2
double E_K_P = -77; // mV
double E_Na_P = 55; // mV
double E_L_P = -65; // mV

double p_rewir;
double p_inhib;

double I_app_min;
double I_app_max;

double g_syn;
double k_syn = 0.2;
double* E_syn;

double g_syn_P;
double k_syn_P = 0.2;
double* E_syn_P;

double alpha_G_P = 25; //s^-1
double beta_G_P = 500; //s^-1

double* I_app;
double* I_app_P;

double** A_A;
double** B_A;
double* C_A;

double** A_N;
double** B_N;
double* C_N;

double** A_N_P;
double** B_N_P;
double* C_N_P;

list<double>* V_spikes;
list<double>* V_spikes_Freq;

FILE* fp_I_syn;

double** tau;

#define tau_min  2 // ms
#define tau_max 12 // ms

#define ms_to_step 40 // (0.001 / dt) !!! Don't forget !!!

#define Max_delay tau_max * ms_to_step
double** V_old_array;

double Poisson_Freq; // Hz
//const double Min_magintude = -0.20; // muA/cm^2
double Max_magnitude = 2.5; // muA/cm^2 // 0.20
const double Duration = 0.002; // sec

double* Meander_start_from_zero;
double* Meander_width;
double* Meander_height;
double* Meander_interval;
double* last_meander_end;

struct Interval
{
	double T_start;
	double T_end;
};

enum Ca_state
{
	Wait_for_up,
	Wait_for_down
};

Ca_state ca_all_state[Node_count];

vector<Interval> Ca_all_intervals[Node_count];
vector<Interval> Ca_intervals;

//bool thread_count_printed = false;

double I_stim(int i, double t)
{
	if (t < Meander_start_from_zero[i])
		return 0;

	t -= Meander_start_from_zero[i];
	t = fmod(t, Meander_width[i] + Meander_interval[i]);

	return t < Meander_width[i] ? Meander_height[i] : 0;
}

double Ca(int i)
{
	return f[i * Equations_per_node];
}

void SetCa(int i, double value)
{
	f[i * Equations_per_node] = value;
}

double IP3(int i)
{
	return f[i * Equations_per_node + 1];
}

void SetIP3(int i, double value)
{
	f[i * Equations_per_node + 1] = value;
}

double z(int i)
{
	return f[i * Equations_per_node + 2];
}

void Setz(int i, double value)
{
	f[i * Equations_per_node + 2] = value;
}

double G_P(int i)
{
	return f[i * Equations_per_node + 3];
}

void SetG_P(int i, double value)
{
	f[i * Equations_per_node + 3] = value;
}

double V(int i)
{
	return f[i * Equations_per_node + 4];
}

void SetV(int i, double value)
{
	f[i * Equations_per_node + 4] = value;
}

double m(int i)
{
	return f[i * Equations_per_node + 5];
}

void Setm(int i, double value)
{
	f[i * Equations_per_node + 5] = value;
}

double n(int i)
{
	return f[i * Equations_per_node + 6];
}

void Setn(int i, double value)
{
	f[i * Equations_per_node + 6] = value;
}

double h(int i)
{
	return f[i * Equations_per_node + 7];
}

void Seth(int i, double value)
{
	f[i * Equations_per_node + 7] = value;
}

///

double V_P(int i)
{
	return f[i * Equations_per_node + 8];
}

void SetV_P(int i, double value)
{
	f[i * Equations_per_node + 8] = value;
}

double m_P(int i)
{
	return f[i * Equations_per_node + 9];
}

void Setm_P(int i, double value)
{
	f[i * Equations_per_node + 9] = value;
}

double n_P(int i)
{
	return f[i * Equations_per_node + 10];
}

void Setn_P(int i, double value)
{
	f[i * Equations_per_node + 10] = value;
}

double h_P(int i)
{
	return f[i * Equations_per_node + 11];
}

void Seth_P(int i, double value)
{
	f[i * Equations_per_node + 11] = value;
}

double V_old(int i, int delay)
{
	return V_old_array[i][Max_delay - 1 - delay];
}

int RandomI(int min, int max)
{
	return ((double)rand() / (RAND_MAX - 1)) * (max - min) + min;
}

double RandomD(double min, double max)
{
	return ((double)rand() / RAND_MAX) * (max - min) + min;
}

double J_channel(double* f, int i)
{
	return c_1 * v_1 * pow(IP3(i), 3) * pow(Ca(i), 3) * pow(z(i), 3) * (c_0 / c_1 - (1 + 1 / c_1) * Ca(i)) / pow((IP3(i) + d_1) * (Ca(i) + d_5), 3);
}

double J_PLC(double* f, int i)
{
	return v_4[i] * (Ca(i) + (1 - alpha) * k_4) / (Ca(i) + k_4);
}

double J_leak(double* f, int i)
{
	return c_1 * v_2 * (c_0 / c_1 - (1 + 1 / c_1) * Ca(i));
}

double J_pump(double* f, int i)
{
	return v_3 * pow(Ca(i), 2) / (pow(k_3, 2) + pow(Ca(i), 2));
}

double J_in(double* f, int i)
{
	return v_5 + v_6 * pow(IP3(i), 2) / (pow(k_2, 2) + pow(IP3(i), 2));
}

double J_out(double* f, int i)
{
	return k_1 * Ca(i);
}

double J_Glu(double* f, int i)
{
	/*double J = 0;
	if (E_syn_P[i] == 0)
	{
		J += alpha_Glu / (1 + exp(-(G_P(i) - 0.25) / 0.01));
	}
	return J;*/

	return alpha_Glu / (1 + exp(-(G_P(i) - 0.25) / 0.01));
}

double alpha_m(double* f, int i)
{
	return 0.182 * (V(i) + 35) / (1 - exp(-(V(i) + 35) / 9));
}

double beta_m(double* f, int i)
{
	return -0.124 * (V(i) + 35) / (1 - exp((V(i) + 35) / 9));
}

double alpha_n(double* f, int i)
{
	return 0.02 * (V(i) - 25) / (1 - exp(-(V(i) - 25) / 9));
}

double beta_n(double* f, int i)
{
	return -0.002 * (V(i) - 25) / (1 - exp((V(i) - 25) / 9));
}

double alpha_h(double* f, int i)
{
	return 0.25 * exp(-(V(i) + 90) / 12);
}

double beta_h(double* f, int i)
{
	return 0.25 * exp((V(i) + 62) / 6) / exp((V(i) + 90) / 12);
}

//

double alpha_m_P(double* f, int i)
{
	return 0.182 * (V_P(i) + 35) / (1 - exp(-(V_P(i) + 35) / 9));
}

double beta_m_P(double* f, int i)
{
	return -0.124 * (V_P(i) + 35) / (1 - exp((V_P(i) + 35) / 9));
}

double alpha_n_P(double* f, int i)
{
	return 0.02 * (V_P(i) - 25) / (1 - exp(-(V_P(i) - 25) / 9));
}

double beta_n_P(double* f, int i)
{
	return -0.002 * (V_P(i) - 25) / (1 - exp((V_P(i) - 25) / 9));
}

double alpha_h_P(double* f, int i)
{
	return 0.25 * exp(-(V_P(i) + 90) / 12);
}

double beta_h_P(double* f, int i)
{
	return 0.25 * exp((V_P(i) + 62) / 6) / exp((V_P(i) + 90) / 12);
}

double UllahJung_HodgkinHuxley(int i, double* f, double t)
{
	int in = i / Equations_per_node;
	int il = i % Equations_per_node;

	switch (il)
	{
	case 0: // Ca
	{
		double sum_1 = 0;

		/*for (int j = 0; j < Node_count; j++)
		{
		  sum_1 += d_Ca * (Ca(j) - Ca(in));
		}*/

		for (int j = 0; j < C_A[in]; j++)
		{
			sum_1 += d_Ca * (Ca((int)B_A[in][j]) - Ca(in));
		}

		return J_channel(f, in) - J_pump(f, in) + J_leak(f, in) + J_in(f, in) - J_out(f, in) + sum_1;
	}

	case 1: // IP3
	{
		double sum_2 = 0;

		/*for (int j = 0; j < Node_count; j++)
		{
		sum_2 += d_IP3 * (IP3(j) - IP3(in));
		}*/

		for (int j = 0; j < C_A[in]; j++)
		{
			sum_2 += d_IP3 * (IP3((int)B_A[in][j]) - IP3(in));
		}

		return (IP3_star - IP3(in)) / tau_IP3 + J_PLC(f, in) + sum_2 + J_Glu(f, in);
	}

	case 2: // z
	{
		return a_2 * (d_2 * (IP3(in) + d_1) / (IP3(in) + d_3) * (1 - z(in)) - Ca(in) * z(in));
	}

	case 3: // G_P
	{
		return -alpha_G_P * G_P(in) + beta_G_P * (1 / (1 + exp(-V_P(in) / 0.5)));
	}

	case 4: // V
	{
		double I_syn = 0;
		double I_syn_P = 0;

		/*for (int j = 0; j < Node_count; j++)
		{
		//sum += A[in][j] * g_syn * (V(in) - V_old(j, tau[in][j]));
		//sum += A[in][j] * g_syn * (V(j) - V(in));
		//sum += A[in][j] * g_syn * (V(in) - E_syn[in]) / (1 + exp(-V_old(j, tau[in][j]) / k_syn));
		//sum += 1 / (0.2 * Node_count_half) * A[in][j] * g_syn * (V(in) - E_syn[in]) / (1 + exp(-V(j) / k_syn)); // i up, j down
		//sum += 1 / (0.2 * Node_count_half) * A[in][j] * g_syn * (V(j) - E_syn[j]) / (1 + exp(-V(in) / k_syn)); // j up, i down
		  I_syn += A_N[in][j] * g_syn * (E_syn[in] - V(in)) / (1 + exp(-(V(j) / k_syn)));
		//printf("in = %d\t Node_count = %d\t A[in][j] = %f\t V(in) = %f\t E_syn[in] = %f\t sum = %f\n", in, j, A[in][j], V(in), E_syn[in], sum);
	  }*/

	  /*for (int j = 0; j < C[in]; j++)
	  {
	  sum += sigma[in][(int)B[in][j]] * (V((int)B[in][j]) - V(in));
	  }*/

		for (int j = 0; j < C_N[in]; j++)
		{
			/*if ((1 + g_astro * Ca(in)) > 0)
			{
				if (Ca(in) >= 0.3)
				{
					I_syn += g_syn * (1 + g_astro * Ca(in)) * (V(in) - E_syn[(int)B_N[in][j]]) / (1 + exp(-(V((int)B_N[in][j]) / k_syn))); // версия с V (без V_old)
				}
				else
				{
					I_syn += g_syn * (V(in) - E_syn[(int)B_N[in][j]]) / (1 + exp(-(V((int)B_N[in][j]) / k_syn))); // версия с V (без V_old)
				}
			}*/
			I_syn += g_syn * (V(in) - E_syn[(int)B_N[in][j]]) / (1 + exp(-(V((int)B_N[in][j]) / k_syn))); // версия с V (без V_old)
			//I_syn += g_syn * (1 + g_astro * Ca(in)) * (E_syn[in] - V(in)) / (1 + exp(-(V_old((int)B_N[in][j], tau[in][(int)B_N[in][j]]) / k_syn)); // версия с V_old
			//I_syn += g_syn * (1 + g_astro * Ca(in)) * (E_syn[in] - V(in)) / (1 + exp(-(V(j) / k_syn))); // версия без с V (без V_old)
			//I_syn += g_syn * (E_syn[in] - V(in)) / (1 + exp(-(V(j) / k_syn))); // версия без с V (без V_old), упрощенная версия
			// sum_3 += g_syn * (1 + g_astro * Ca(in)) * (E_syn[i] - V(i)) / (1 + exp(-(V(j) / k_syn))); // образец из старой версии
			//I_syn += g_syn * (E_syn[in] - V(in)) / (1 + exp(-(V_old((int)B_N[in][j], tau[in][(int)B_N[in][j]])) / k_syn));*/
		  //I_syn += g_syn * (E_syn[(int)B_N[in][j]] - V(in)) / (1 + exp(-(V(j) / k_syn))); // версия с V (без V_old)
		  //I_syn += g_syn * (E_syn[(int)B_N[in][j]] - V(in)) / (1 + exp(-(V((int)B_N[in][j]) / k_syn))); // версия с V (без V_old), упрощенная версия !!!
		  //I_syn += g_syn * (E_syn[(int)B_N[in][j]] - V(in)) / (1 + exp(-(V_old((int)B_N[in][j], tau[in][(int)B_N[in][j]]) / k_syn))); // версия с V (без V_old), упрощенная версия !!!
		//sum += g_syn * (V((int)B[in][j]) - V(in)); // устаревшая часть, нужна для проверки разностной схемы
		//sum += A[in][j] * g_syn * (V(in) - E_syn[in]) / (1 + exp(-V_old((int)B[in][j]) / k_syn)); // устаревшая часть
		//sum += A[in][(int)B[in][j]] * g_syn * (V(in) - E_syn[in]) / (1 + exp(-V_old((int)B[in][j], tau[in][(int)B[in][j]]) / k_syn));
		//sum += A[in][(int)B[in][j]] * g_syn * (V((int)B[in][j]) - E_syn[(int)B[in][j]]) / (1 + exp(-V_old(in, tau[in][(int)B[in][j]]) / k_syn));
		//sum += 1 / (0.2 * Node_count_half) * A[in][(int)B[in][j]] * g_syn * (V((int)B[in][j]) - E_syn[(int)B[in][j]]) / (1 + exp(-V(in) / k_syn)); // j up, i down
		//sum += 1 / (0.2 * Node_count) * /*(int)A_N[in][(int)B_N[in][j]] * */ g_syn * (1 + g_astro * Ca(in)) * (V(in) - E_syn[in]) / (1 + exp(-V((int)B_N[in][j]) / k_syn)); // i up, j down
		  // i up, j down
		  /*printf("i = %d\t j = %d\t A[i, j] = %d\n", in, (int)B_N[in][j], (int)A_N[in][(int)B_N[in][j]]);*/
		/*printf("i = %d\t V_old = %f\t exp = %f\n", in, Vold, ee);*/
		}
		/*printf("i = %d\t sum = %f\n", in, sum);*/

		if ((1 + g_astro * Ca(in)) > 0)
		{
			if (Ca(in) >= 0.3)
			{
				I_syn_P += g_syn_P * (1 + g_astro * Ca(in)) * (V_P(in) - E_syn_P[in]) / (1 + exp(-(V_P(in) / k_syn_P))); // версия с V (без V_old)
			}
			else
			{
				I_syn_P += g_syn_P * (V_P(in) - E_syn_P[in]) / (1 + exp(-(V_P(in) / k_syn_P))); // версия с V (без V_old)
			}
		}
		//I_syn_P += g_syn_P * (V_P(in) - E_syn_P[in]) / (1 + exp(-(V_P(in) / k_syn_P))); // версия с V (без V_old)

		if (enable_I_syn_out)
			fprintf(fp_I_syn, i == Equations_count - 1 ? "%f" : "%f\t", I_syn);

		return 1000 * ((g_Na * pow(m(in), 3) * h(in) * (E_Na - V(in)) + g_K * n(in) * (E_K - V(in)) + g_L * (E_L - V(in)) + I_app[in] + I_syn + I_syn_P) / C_m); // V
	}

	case 5: // m
	{
		return 1000 * (alpha_m(f, in) * (1 - m(in)) - beta_m(f, in) * m(in)); // m
	}

	case 6: // n
	{
		return 1000 * (alpha_n(f, in) * (1 - n(in)) - beta_n(f, in) * n(in)); // n
	}

	case 7: // h
	{
		return 1000 * (alpha_h(f, in) * (1 - h(in)) - beta_h(f, in) * h(in)); // h
	}

	case 8: // V_P
	{
		return 1000 * ((g_Na_P * pow(m_P(in), 3) * h_P(in) * (E_Na_P - V_P(in)) + g_K_P * n_P(in) * (E_K_P - V_P(in)) + g_L_P * (E_L_P - V_P(in)) + I_app_P[in] + I_stim(in, t)) / C_m_P); // V_P
	}

	case 9: // m_P
	{
		return 1000 * (alpha_m_P(f, in) * (1 - m_P(in)) - beta_m_P(f, in) * m_P(in)); // m_P
	}

	case 10: // n_P
	{
		return 1000 * (alpha_n_P(f, in) * (1 - n_P(in)) - beta_n_P(f, in) * n_P(in)); // n_P
	}

	case 11: // h_P
	{
		return 1000 * (alpha_h_P(f, in) * (1 - h_P(in)) - beta_h_P(f, in) * h_P(in)); // h_P
	}
	}

	return 0;
}

void RungeKutta(double t, double dt, double* f, double* f_next)
{
	// k1
#ifdef USE_OMP
#pragma omp parallel for
#endif
	for (int i = 0; i < Equations_count; i++)
	{
		//if (!thread_count_printed)
		//{
		//	thread_count_printed = true;
		//	printf("Threads = %d\n", omp_get_num_threads());
		//}

		k[i][0] = UllahJung_HodgkinHuxley(i, f, t) * dt;
		phi_k1[i] = f[i] + k[i][0] / 2;
		k[i][1] = UllahJung_HodgkinHuxley(i, phi_k1, t) * dt;
		phi_k2[i] = f[i] + k[i][1] / 2;
		k[i][2] = UllahJung_HodgkinHuxley(i, phi_k2, t) * dt;
		phi_k3[i] = f[i] + k[i][2] / 2;
		k[i][3] = UllahJung_HodgkinHuxley(i, phi_k3, t) * dt;
		f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
	}

	//for (int i = 0; i < Equations_count; i++)
	//	phi_k1[i] = f[i] + k[i][0] / 2;

	// k2
	//for (int i = 0; i < Equations_count; i++)
	//	k[i][1] = UllahJung_HodgkinHuxley(i, phi_k1, t) * dt;


	//for (int i = 0; i < Equations_count; i++)
	//	phi_k2[i] = f[i] + k[i][1] / 2;

	// k3
	//for (int i = 0; i < Equations_count; i++)
	//	k[i][2] = UllahJung_HodgkinHuxley(i, phi_k2, t) * dt;


	//for (int i = 0; i < Equations_count; i++)
	//	phi_k3[i] = f[i] + k[i][2] / 2;

	//enable_I_syn_out = true;

	// k4
	//for (int i = 0; i < Equations_count; i++)
	//	k[i][3] = UllahJung_HodgkinHuxley(i, phi_k3, t) * dt;

	//enable_I_syn_out = false;

	//for (int i = 0; i < Equations_count; i++)
	//	f_next[i] = f[i] + (k[i][0] + 2 * k[i][1] + 2 * k[i][2] + k[i][3]) / 6;
}

void CopyArray(double* source, double* target, int N)
{
	for (int i = 0; i < N; i++)
		target[i] = source[i];
}

bool Approximately(double a, double b)
{
	if (a < 0)
		a = -a;

	if (b < 0)
		b = -b;

	return a - b <= 0.000001;
}

//bool CheckSameLine(int i, int j)
//{
//	return i / Node_wire_width == j / Node_wire_width;
//}
//
//bool IsWireNeighbors(int i, int j)
//{
//	if (CheckSameLine(i, j) && (i == j - 1 || i == j + 1))
//		return true;
//
//	if (i == j - Node_wire_width || i == j + Node_wire_width)
//		return true;
//
//	return false;
//}

// http://preshing.com/20111007/how-to-generate-random-timings-for-a-poisson-process/
double nextTime(double rateParameter)
{
	return -log(1.0 - (double)rand() / (RAND_MAX)) / rateParameter;
}

void GenerateRandomMeander(int i, double min_start_time)
{
	double offset = nextTime(Poisson_Freq);

	if (offset < 0)
	{
		int a = 0;
	}

	Meander_start_from_zero[i] = min_start_time + offset;
	Meander_width[i] = Duration;
	//Meander_height[i] = RandomD(-Max_magnitude, Max_magnitude);
	Meander_height[i] = RandomD(0, Max_magnitude);
	//Meander_height[i] = Max_magnitude;
}

void FillAMatrixZero()
{
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			A_A[i][j] = 0;
			A_N[i][j] = 0;
			A_N_P[i][j] = 0;
		}
	}
}

void FillAstrociteMatrix()
{
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			if (i == j)
			{
				A_A[i][j] = 0;
				continue;
			}

			if (i > j)
			{
				A_A[i][j] = A_A[j][i];
				continue;
			}

			if (i == 0 && j == Node_count - 1)
			{
				A_A[i][j] = 1;
				continue;
			}

			if (i == Node_count - 1 && j == 0)
			{
				A_A[i][j] = 1;
				continue;
			}

			if (i == j - 1 || i == j + 1)
			{
				A_A[i][j] = 1;
				continue;
			}
		}
	}
	//A_A[0][1] = 0; // only for debug. diffusion Ca test
	//A_A[1][0] = 0; // only for debug. diffusion Ca test
	//A_A[1][3] = 0;
	//A_A[3][1] = 0;
	//A_A[0][2] = 0;
	//A_A[2][0] = 0;
}

void FillBCMatrix_A()
{
	for (int i = 0; i < Node_count; i++)
	{
		int bIndex = 0;
		C_A[i] = 0;
		for (int j = 0; j < Node_count; j++)
		{
			if (A_A[i][j] == 1)
			{
				B_A[i][bIndex] = j;
				bIndex++;
				C_A[i]++;
			}
		}
	}
}

void FillBCMatrix_N()
{
	for (int i = 0; i < Node_count; i++)
	{
		int bIndex = 0;
		C_N[i] = 0;
		for (int j = 0; j < Node_count; j++)
		{
			if (A_N[i][j] == 1)
			{
				B_N[i][bIndex] = j;
				bIndex++;
				C_N[i]++;
			}
		}
	}
}

void FillBCMatrix_N_P()
{
	for (int i = 0; i < Node_count; i++)
	{
		int bIndex = 0;
		C_N_P[i] = 0;
		for (int j = 0; j < Node_count; j++)
		{
			if (A_N_P[i][j] == 1)
			{
				B_N_P[i][bIndex] = j;
				bIndex++;
				C_N_P[i]++;
			}
		}
	}
}

bool IsWireNeighbors(int i, int j, int deep)
{
	int j_border_left = j - deep < 0 ? j + Node_count : j;
	int j_border_right = j + deep >= Node_count ? j - Node_count : j;

	if (i == j_border_left - deep || i == j_border_right + deep)
	{
		return true;
	}

	return false;
}

void FillNeuronMatrix()
{
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			if (i == j)
			{
				A_N[i][j] = 0;
				continue;
			}

			if (i > j)
			{
				A_N[i][j] = A_N[j][i];
				continue;
			}

			for (int deep = 1; deep <= MaxDeep; deep++)
			{
				if (IsWireNeighbors(i, j, deep))
					A_N[i][j] = 1;
			}
		}
	}
}

void RandomizeNeuronMatrix()
{
	srand(time(NULL));

	for (int i = 0; i < Node_count; i++)
	{
		//if (i == Node_count / 2)
		//	srand(time(NULL));

		for (int link = 0; link < MaxDeep * 2; link++)
		{
			double x = RandomD(0, 1);

			if (x > p_rewir)
				continue;

			int rndJ;

			do
			{
				rndJ = RandomI(0, Node_count);
			} while (i == rndJ || A_N[i][rndJ] == 1);

			int rndJ_last;

			do
			{
				rndJ_last = RandomI(i - MaxDeep - 1, i + MaxDeep + 1);

				if (rndJ_last < 0)
					rndJ_last += Node_count;
				else if (rndJ_last >= Node_count)
					rndJ_last -= Node_count;

			} while (i == rndJ_last || A_N[i][rndJ_last] == 0);

			A_N[i][rndJ_last] = 0;

			A_N[i][rndJ] = 1;
		}
	}
}

void FillNeuronPoissonMatrix()
{
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			A_N_P[i][j] = 0;
		}
	}
}

void FillVOldFromCurrent()
{
	for (int i = 0; i < Node_count; i++)
		for (int j = 0; j < Max_delay; j++)
			V_old_array[i][j] = V(i);
}

void UpdateVOld()
{
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 1; j < Max_delay; j++)
			V_old_array[i][j - 1] = V_old_array[i][j];

		V_old_array[i][Max_delay - 1] = V(i);
	}
}

//void FillFullTauMatrix()
//{
//	for (int i = 0; i < Node_count; i++)
//	{
//		for (int j = 0; j < Node_count; j++)
//		{
//			if (i < Node_count || j < Node_count)
//			{
//				tau[i][j] = 0;
//				continue;
//			}
//
//			int i_neuron = i - Node_count;
//			int j_neuron = j - Node_count;
//
//			int i_wire_x = i_neuron / Node_wire_width;
//			int i_wire_y = i_neuron % Node_wire_width;
//
//			int j_wire_x = j_neuron / Node_wire_width;
//			int j_wire_y = j_neuron % Node_wire_width;
//
//			double distance_max = sqrt(2.) * (Node_wire_width - 1);
//			double distance = sqrt((i_wire_x - j_wire_x) * (i_wire_x - j_wire_x) + (i_wire_y - j_wire_y) * (i_wire_y - j_wire_y));
//
//			tau[i][j] = (tau_min + distance / (distance_max) * (tau_max - tau_min)) * ms_to_step;
//		}
//	}
//}
//
//void FillTauMatrix()
//{
//	for (int i = 0; i < Node_count; i++)
//	{
//		for (int j = 0; j < Node_count; j++)
//		{
//			if (i == j || A_N[i][j] == 0)
//			{
//				tau[i][j] = 0;
//				continue;
//			}
//
//			int i_neuron = i;
//			int j_neuron = j;
//
//			int i_wire_x = i_neuron / Node_wire_width;
//			int i_wire_y = i_neuron % Node_wire_width;
//
//			int j_wire_x = j_neuron / Node_wire_width;
//			int j_wire_y = j_neuron % Node_wire_width;
//
//			double distance_max = sqrt(2.) * (Node_wire_width - 1);
//			double distance = sqrt((i_wire_x - j_wire_x) * (i_wire_x - j_wire_x) + (i_wire_y - j_wire_y) * (i_wire_y - j_wire_y));
//
//			double t = (distance - 1) / (distance_max - 1);
//			tau[i][j] = (tau_min + t * (tau_max - tau_min)) * ms_to_step;
//		}
//	}
//}

int main(int argc, char* argv[])
{
	// run like: UJ_HH_Ring_acc.out 250 0.2 0.4 0.05 3.0 1.05 1.50 // (1) p_rewir (2) g_syn (3) g_astro (4) Max_magnitude (5) g_syn_P (6) Poisson_Freq
	/*sscanf(argv[1], "%d", &Node_count);
	FILE* fp_Node_count;
	fp_Node_count = fopen("Node_count.txt", "w");
	fprintf(fp_Node_count, "%d\t", Node_count);
	fclose(fp_Node_count);
	printf("Node_count = %d\n", Node_count);*/

	sscanf(argv[1], "%lf", &p_rewir);
	FILE* fp_p_rewir;
	fp_p_rewir = fopen("p_rewir.txt", "w");
	fprintf(fp_p_rewir, "%f\t", p_rewir);
	fclose(fp_p_rewir);
	printf("p_rewir = %f\n", p_rewir);

	/*sscanf(argv[3], "%lf", &p_inhib);
	FILE* fp_p_inhib;
	fp_p_inhib = fopen("p_inhib.txt", "w");
	fprintf(fp_p_inhib, "%f\t", p_inhib);
	fclose(fp_p_inhib);
	printf("p_inhib = %f\n", p_inhib);*/

	sscanf(argv[2], "%lf", &g_syn);
	FILE* fp_g_syn;
	fp_g_syn = fopen("g_syn.txt", "w");
	fprintf(fp_g_syn, "%f\t", g_syn);
	fclose(fp_g_syn);
	printf("g_syn = %f\n", g_syn);

	sscanf(argv[3], "%lf", &g_astro);
	FILE* fp_g_astro;
	fp_g_astro = fopen("g_astro.txt", "w");
	fprintf(fp_g_astro, "%f\t", g_astro);
	fclose(fp_g_astro);
	printf("g_astro = %f\n", g_astro);

	sscanf(argv[4], "%lf", &g_syn_P);
	FILE* fp_g_syn_P;
	fp_g_syn_P = fopen("g_syn_P.txt", "w");
	fprintf(fp_g_syn_P, "%f\t", g_syn_P);
	fclose(fp_g_syn_P);
	printf("g_syn_P = %f\n", g_syn_P);

	sscanf(argv[5], "%lf", &Poisson_Freq);
	FILE* fp_Poisson_Freq;
	fp_Poisson_Freq = fopen("Poisson_Freq.txt", "w");
	fprintf(fp_Poisson_Freq, "%f\t", Poisson_Freq);
	fclose(fp_Poisson_Freq);
	printf("Poisson_Freq = %f\n", Poisson_Freq);

	int write = 0;
	sscanf(argv[6], "%i", &write);
	write_file = write;
	printf("write_file = %i\n", write_file);

	int omp = 0;
	sscanf(argv[7], "%i", &omp);
	omp_enable = omp;
	printf("omp_enable = %i\n", omp_enable);


#ifdef USE_OMP
	printf("Macro USE_OMP seted\n");
#endif

#ifndef USE_OMP
	printf("Macro USE_OMP not seted\n");
#endif

#ifndef USE_OMP
	if (omp_enable)
	{
		printf("WARNING! omp_enable == true but macro #USE_OMP not seted, OMP will not be used \n");
		omp_enable = false;
	}
#endif

	//sscanf(argv[6], "%lf", &I_app_min);
	//sscanf(argv[7], "%lf", &I_app_max);

	for (int i = 0; i < Node_count; i++)
		ca_all_state[i] = Wait_for_up;

	f = new double[Equations_count];
	f_diff = new double[Equations_count];
	v_4 = new double[Node_count];
	E_syn = new double[Node_count];
	I_app = new double[Node_count];
	E_syn_P = new double[Node_count];
	I_app_P = new double[Node_count];
	V_spikes = new list<double>[Node_count];
	V_spikes_Freq = new list<double>[Node_count];

	Meander_start_from_zero = new double[Node_count];
	Meander_width = new double[Node_count];
	Meander_height = new double[Node_count];
	Meander_interval = new double[Node_count];
	last_meander_end = new double[Node_count];

	tau = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		tau[i] = new double[Node_count];

	V_old_array = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		V_old_array[i] = new double[Max_delay];

	FILE* fp0;
	FILE* fp_I_stim = 0;
	FILE* fp_Ca = 0;
	FILE* fp_IP3 = 0;
	//FILE *fp_z;
	//FILE* fp_G_P;
	FILE* fp_V = 0;
	FILE* fp_V_P = 0;
	//FILE *fp_m;
	//FILE *fp_n;
	//FILE *fp_h;
	FILE* fp_V_spikes = 0;
	FILE* fp_Esyn;
	FILE* fp_Esyn_P;

	//FILE* fp_res;
	srand(time(NULL));

	//for (int i = 0; i < Node_count; i++)
	//	V_old_length[i] = 0;

	A_A = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		A_A[i] = new double[Node_count];

	B_A = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		B_A[i] = new double[Node_count];

	C_A = new double[Node_count];

	A_N = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		A_N[i] = new double[Node_count];

	B_N = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		B_N[i] = new double[Node_count];

	C_N = new double[Node_count];

	A_N_P = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		A_N_P[i] = new double[Node_count];

	B_N_P = new double* [Node_count];
	for (int i = 0; i < Node_count; i++)
		B_N_P[i] = new double[Node_count];

	C_N_P = new double[Node_count];

	FillAMatrixZero();
	FillAstrociteMatrix();
	FillNeuronMatrix();
	RandomizeNeuronMatrix();
	FillNeuronPoissonMatrix();
	FillBCMatrix_A();
	FillBCMatrix_N();
	FillBCMatrix_N_P();
	//FillTauMatrix();

	fp0 = fopen("A_A.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			fprintf(fp0, "%d\t", (int)A_A[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("A_N.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			fprintf(fp0, "%d\t", (int)A_N[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("tau.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < Node_count; j++)
		{
			fprintf(fp0, "%f\t", tau[i][j] / ms_to_step);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	// Write to file number of links for each neuron
	fp0 = fopen("links.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		int links_count = 0;
		for (int j = 0; j < Node_count; j++)
		{
			if (A_N[i][j] == 1)
			{
				links_count++;
			}
		}
		fprintf(fp0, "%d\n", (int)links_count);
	}
	fclose(fp0);

	//setlocale(LC_NUMERIC, "French_Canada.1252");
	fp0 = fopen("test_Poisson.txt", "w+");
	for (int i = 0; i < 1000; i++)
		fprintf(fp0, "%f\n", nextTime(Poisson_Freq));
	fclose(fp0);

	fp0 = fopen("B_A.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < C_A[i]; j++)
		{
			fprintf(fp0, "%d\t", (int)B_A[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("B_N.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		for (int j = 0; j < C_N[i]; j++)
		{
			fprintf(fp0, "%d\t", (int)B_N[i][j]);
		}
		fprintf(fp0, "\n");
	}
	fclose(fp0);

	fp0 = fopen("C_A.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		fprintf(fp0, "%d\n", (int)C_A[i]);
	}
	fclose(fp0);

	fp0 = fopen("C_N.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		fprintf(fp0, "%d\n", (int)C_N[i]);
	}
	fclose(fp0);

	/*for (int i = 0; i < 6; i++)
	{
		v_4[i] = 0.6;
	}*/
	for (int i = 0; i < Node_count; i++)
	{
		v_4[i] = 0.4; // 0.4
	}

	// Initial values
	/*for (int i = 0; i < Equations_count; i++)
	{
	  f[i] = 0;
  }*/

  //I_app_min = 1.1;
  //I_app_max = 1.5;

	FILE* fp_I_app;
	fp_I_app = fopen("I_app.txt", "w");
	for (int i = 0; i < Node_count; i++)
	{
		I_app[i] = 0.7; //RandomD(I_app_min, I_app_max);
		fprintf(fp_I_app, "%f\n", I_app[i]);
	}
	fclose(fp_I_app);

	FILE* fp_I_app_P;
	fp_I_app_P = fopen("I_app_P.txt", "w");
	for (int i = 0; i < Node_count; i++)
	{
		I_app_P[i] = 0.7;
		fprintf(fp_I_app_P, "%f\n", I_app_P[i]);
	}
	fclose(fp_I_app_P);

	for (int i = 0; i < Node_count; i++) // init array for all nodes
	{
		SetG_P(i, 0); // G_P
	}

	double percent_stable_state = 0.50; // 0.40
	double eps_persent = 0.05; //0.05

	double Ca0 = 0.07;
	double IP30 = 0.16;
	double z0 = 0.67;

	for (int i = 0; i < Node_count; i++)
	{
		/*SetCa(i, Ca0 + RandomD(-Ca0 * eps_persent, Ca0 * eps_persent)); // Ca
		SetIP3(i, IP30 + RandomD(-IP30 * eps_persent, IP30 * eps_persent)); // IP3
		Setz(i, z0 + RandomD(-z0 * eps_persent, z0 * eps_persent)); // z */
		SetCa(i, Ca0); // Ca
		SetIP3(i, IP30); // IP3
		Setz(i, z0); // z
	}

	double V0 = -58.7085;
	double m0 = 0.0953;
	double n0 = 0.000913;
	double h0 = 0.3662;

	double V1 = 14.8409;
	double m1 = 0.9174;
	double n1 = 0.0140;
	double h1 = 0.0539;

	/*for (int i = 0; i < Node_count; i++) // init only for neurons
	{
	  double random = RandomD(0, 1);

	  SetV(i, random < percent_stable_state ? V0 + RandomD(-V0 * eps_persent, V0 * eps_persent) : V1 + RandomD(-V1 * eps_persent, V1 * eps_persent)); // V
	  Setm(i, random < percent_stable_state ? m0 + RandomD(-m0 * eps_persent, m0 * eps_persent) : m1 + RandomD(-m1 * eps_persent, m1 * eps_persent)); // m
	  Setn(i, random < percent_stable_state ? n0 + RandomD(-n0 * eps_persent, n0 * eps_persent) : n1 + RandomD(-n1 * eps_persent, n1 * eps_persent)); // n
	  Seth(i, random < percent_stable_state ? h0 + RandomD(-h0 * eps_persent, h0 * eps_persent) : h1 + RandomD(-h1 * eps_persent, h1 * eps_persent)); // h
	}*/

	for (int i = 0; i < Node_count; i++) // init only for neurons
	{
		double random = RandomD(0, 1);

		/*SetV(i, random < percent_stable_state ? V0 : V1); // V
		Setm(i, random < percent_stable_state ? m0 : m1); // m
		Setn(i, random < percent_stable_state ? n0 : n1); // n
		Seth(i, random < percent_stable_state ? h0 : h1); // h*/

		SetV(i, V1); // V
		Setm(i, m1); // m
		Setn(i, n1); // n
		Seth(i, h1); // h

		/*SetV_P(i, random < percent_stable_state ? V0 : V1); // V
		Setm_P(i, random < percent_stable_state ? m0 : m1); // m
		Setn_P(i, random < percent_stable_state ? n0 : n1); // n
		Seth_P(i, random < percent_stable_state ? h0 : h1); // h*/

		SetV_P(i, V0); // V
		Setm_P(i, m0); // m
		Setn_P(i, n0); // n
		Seth_P(i, h0); // h
	}

	/*for (int i = 0; i < Node_count; i++) // init only for neurons
	{*/
	/*SetV(i, RandomD(-80, 20)); // V
	Setm(i, RandomD(0, 1)); // m
	Setn(i, RandomD(0, 1)); // n
	Seth(i, RandomD(0, 1)); // h*/
	/*SetV(i, V1); // V
	Setm(i, m1); // m
	Setn(i, n1); // n
	Seth(i, h1); // h
}*/

	double E_syn0 = 0; // Excitatory neuron
	double E_syn1 = -90; // Inhibitory neuron

	fp_Esyn = fopen("results_E_syn.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		/*E_syn[i] = E_syn0;

		double x = RandomD(0, 1);

		if (x > p_inhib)
		{
			fprintf(fp_Esyn, "%f\n", E_syn[i]);
			continue;
		}*/

		E_syn[i] = E_syn1;

		fprintf(fp_Esyn, "%f\n", E_syn[i]);

		E_syn_P[i] = E_syn0;
	}
	fclose(fp_Esyn);

	for (int i = 0; i < Node_count; i++)
	{
		GenerateRandomMeander(i, 0);
		last_meander_end[i] = Meander_start_from_zero[i] + Duration;
	}

	const double t_start = 0;
	const double t_max = 34; // 100 msec = 0.1 sec // 240 // 270
	//const double dt = 0.000002; // 0.01 msec = 0.00001 sec; 0.1 msec = 0.0001 sec; 1 msec = 0.001 sec // 0.000025
	const double dt = 0.00002; // 0.01 msec = 0.00001 sec; 0.1 msec = 0.0001 sec; 1 msec = 0.001 sec // 0.000025

	double t = t_start;

	//fp_res = fopen("matlab_res.txt", "a+");
	//fprintf(fp_res, "%f\t%f\t", Max_magintude, g_syn);

	//fp0 = fopen("results.txt", "w+");
	//setlocale(LC_NUMERIC, "French_Canada.1252");

	double start_rk4, end_rk4;
	clock_t omp_start_rk4, omp_end_rk4;

	if (omp_enable)
	{
#ifdef USE_OMP
		omp_start_rk4 = omp_get_wtime();
#endif
	}
	else
		start_rk4 = clock();

	int lastPercent = -1;

	//FillVOldFromCurrent();

	k = new double* [Equations_count];
	for (int i = 0; i < Equations_count; i++)
		k[i] = new double[4];

	phi_k1 = new double[Equations_count];
	phi_k2 = new double[Equations_count];
	phi_k3 = new double[Equations_count];

	if (write_file)
	{
		fp_I_stim = fopen("results_I_stim.txt", "w+");
		//fp_I_syn = fopen("results_I_syn.txt", "w+");
		fp_Ca = fopen("results_Ca.txt", "w+");
		//fp_IP3 = fopen("results_IP3.txt", "w+");
		//fp_z    = fopen("results_z.txt", "w+");
		//fp_G_P = fopen("results_G_P.txt", "w+");
		fp_V = fopen("results_V.txt", "w+");
		fp_V_P = fopen("results_V_P.txt", "w+");
		//fp_m = fopen("results_m.txt", "w+");
		//fp_n = fopen("results_n.txt", "w+");
		//fp_h = fopen("results_h.txt", "w+");
		fp_V_spikes = fopen("results_V_spikes.txt", "w+");
	}

	//

	double* f_next = new double[Equations_count];

	while (t < t_max || Approximately(t, t_max))
	{
		if (write_file)
		{
			fprintf(fp_I_stim, "%f\t", t);
			fprintf(fp_Ca, "%f\t", t);
			//fprintf(fp_IP3, "%f\t", t);
			//fprintf(fp_z, "%f\t", t);
			//fprintf(fp_G, "%f\t", t);
			fprintf(fp_V, "%f\t", t);
			fprintf(fp_V_P, "%f\t", t);
			//fprintf(fp_m, "%f\t", t);
			//fprintf(fp_n, "%f\t", t);
			//fprintf(fp_h, "%f\t", t);
			//fprintf(fp_V_spikes, "%f\t", t);
		}

		for (int i = 0; i < Node_count; i++)
		{
			if (t > last_meander_end[i])
			{
				GenerateRandomMeander(i, t);
				last_meander_end[i] = Meander_start_from_zero[i] + Duration;
			}

			if (write_file)
				fprintf(fp_I_stim, "%f\t", I_stim(i, t));
		}

		if (write_file)
		{
			fprintf(fp_I_stim, "\n");

			for (int i = 0; i < Equations_count; i += Equations_per_node)
				fprintf(fp_Ca, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // Ca

			//for (int i = 1; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_IP3, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // IP3

			//for (int i = 2; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_z, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // z

			//for (int i = 3; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_G_P, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // G

			for (int i = 4; i < Equations_count; i += Equations_per_node)
			{
				if (isnan(f[i]))
					return 1;

				fprintf(fp_V, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // V
			}

			//for (int i = 5; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_m, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // m

			//for (int i = 6; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_n, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // n

			//for (int i = 7; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_h, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // h

			for (int i = 8; i < Equations_count; i += Equations_per_node)
				fprintf(fp_V_P, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // V_P

			//for (int i = 9; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_m_P, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // m_P

			//for (int i = 10; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_n_P, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // n_P

			//for (int i = 11; i < Equations_count; i += Equations_per_node)
			//	fprintf(fp_h_P, i == Equations_count - 1 ? "%f" : "%f\t", f[i]); // h_P

			fprintf(fp_Ca, "\n");
			//fprintf(fp_IP3, "\n");
			//fprintf(fp_z, "\n");
			//fprintf(fp_G, "\n");
			fprintf(fp_V, "\n");
			//fprintf(fp_m, "\n");
			//fprintf(fp_n, "\n");
			//fprintf(fp_h, "\n");

			//fprintf(fp_G_P, "\n");
			fprintf(fp_V_P, "\n");
			//fprintf(fp_m_P, "\n");
			//fprintf(fp_n_P, "\n");
			//fprintf(fp_h_P, "\n");
		}

		RungeKutta(t, dt, f, f_next);

#ifdef USE_OMP
#pragma omp parallel for
#endif
		for (int i = 0; i < Node_count; i++)
		{
			int index = Equations_per_node * i + 4;
			double diff = f_next[index] - f[index];

			//fprintf(fp_V_spikes, i == Equations_count - 1 ? "%d" : "%d\t", diff < 0 && f_diff[index] > 0 && f[index] > -10 && (V_spikes[i].size() == 0 || t - V_spikes[i].back() > 0.001) ? 1 : 0);

			if (diff < 0 && f_diff[index] > 0 && f[index] > -10 && (V_spikes[i].size() == 0 || t - V_spikes[i].back() > 0.001))
			{
				V_spikes[i].push_back(t);
			}

			f_diff[index] = diff;
		}

		//fprintf(fp_V_spikes, "\n");

		CopyArray(f_next, f, Equations_count);

		//printf("V(24) = %f\t V_old(24) = %f\n", f[24*4], V_old(24));
		//UpdateVOld();

		//fprintf(fp_I_syn, "\n");

		for (int i = 0; i < Node_count; i++)
		{
			if (ca_all_state[i] == Wait_for_up)
			{
				if (Ca(i) > 0.3)
				{
					ca_all_state[i] = Wait_for_down;
					Interval interval;
					interval.T_start = t;
					interval.T_end = -1;
					Ca_all_intervals[i].push_back(interval);
				}
			}
			else if (ca_all_state[i] == Wait_for_down)
			{
				if (Ca(i) < 0.3)
				{
					Ca_all_intervals[i].back().T_end = t;
					ca_all_state[i] = Wait_for_up;
				}
			}
		}

		t += dt;

		/*int percent = (int)(100 * (t - t_start) / (t_max - t_start));
		if (percent != lastPercent)
		{
			printf("Progress: %d%%\n", percent);
			lastPercent = percent;
		}*/
	}

	delete[] f_next;

	double* V_mean_freqs = new double[Node_count];
	double* V_STD_freqs = new double[Node_count];
	//list<double> V_mean_freqs;

	int max_interval_count = 0;

	for (int i = 0; i < Node_count; i++)
	{
		if (Ca_all_intervals[i].size() > max_interval_count)
		{
			max_interval_count = Ca_all_intervals[i].size();
		}
	}

	for (int i = 0; i < max_interval_count; i++)
	{
		Interval mean_interval;
		mean_interval.T_start = 0;
		mean_interval.T_end = 0;

		int mean_count = 0;

		for (int j = 0; j < Node_count; j++)
		{
			if (Ca_all_intervals[j].size() > i && Ca_all_intervals[j][i].T_end > 0)
			{
				mean_interval.T_start += Ca_all_intervals[j][i].T_start;
				mean_interval.T_end += Ca_all_intervals[j][i].T_end;

				mean_count++;
			}
		}

		mean_interval.T_start /= mean_count;
		mean_interval.T_end /= mean_count;

		Ca_intervals.push_back(mean_interval);
	}

#ifdef USE_OMP
#pragma omp parallel for
#endif
	for (int i = 0; i < Node_count; i++)
	{
		list<double>::iterator it_V_spikes = V_spikes[i].begin();

		while (it_V_spikes != V_spikes[i].end() && *it_V_spikes < 0.25 * t_max)
		{
			V_spikes[i].pop_front();
			it_V_spikes = V_spikes[i].begin();
		}

		list<double> V_freqs;

		it_V_spikes = V_spikes[i].begin();

		V_mean_freqs[i] = 0;
		V_STD_freqs[i] = 0;

		if (V_spikes[i].size() <= 1)
		{
			continue;
		}
		else
		{
			for (int j = 1; j < V_spikes[i].size(); j++)
			{
				double first = *it_V_spikes;
				advance(it_V_spikes, 1);
				double next = *it_V_spikes;

				double T = next - first;
				V_freqs.push_back(1 / T);
			}
		}

		list<double>::iterator it_Freq = V_freqs.begin();

		for (int j = 0; j < V_freqs.size(); j++)
		{
			V_mean_freqs[i] += *it_Freq;
			advance(it_Freq, 1);
		}

		V_mean_freqs[i] /= V_freqs.size();

		it_Freq = V_freqs.begin();

		for (int j = 0; j < V_freqs.size(); j++)
		{
			V_STD_freqs[i] += pow(*it_Freq - V_mean_freqs[i], 2);
			advance(it_Freq, 1);
		}

		V_STD_freqs[i] /= V_freqs.size();
		V_STD_freqs[i] = sqrt(V_STD_freqs[i]);
	}

	fp0 = fopen("V_STD_freqs.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		fprintf(fp0, "%f\t", V_STD_freqs[i]);
	}
	fclose(fp0);

	double V_STD_mean_Freq = 0;
	double V_STD_mean_Freq_not_zero = 0;
	double V_STD_mean_Freq_not_zero_count = 0;

	for (int j = 0; j < Node_count; j++)
	{
		if (V_STD_freqs[j] != 0)
		{
			V_STD_mean_Freq_not_zero_count++;
			V_STD_mean_Freq += V_STD_freqs[j];
		}
	}

	V_STD_mean_Freq_not_zero = V_STD_mean_Freq / V_STD_mean_Freq_not_zero_count;
	V_STD_mean_Freq /= Node_count;

	fp0 = fopen("V_STD_mean_Freq.txt", "w+");
	fprintf(fp0, "%f\n", V_STD_mean_Freq);
	//fprintf(fp_res, "%f\t", V_STD_mean_Freq);
	fclose(fp0);

	fp0 = fopen("V_STD_mean_Freq_not_zero.txt", "w+");
	fprintf(fp0, "%f\n", V_STD_mean_Freq_not_zero);
	//fprintf(fp_res, "%f\t", V_STD_mean_Freq_not_zero);
	fclose(fp0);

	double V_mean_mean_Freq = 0;
	double V_mean_mean_Freq_not_zero = 0;
	double V_mean_mean_Freq_not_zero_count = 0;

	for (int j = 0; j < Node_count; j++)
	{
		if (V_mean_freqs[j] != 0)
		{
			V_mean_mean_Freq_not_zero_count++;
			V_mean_mean_Freq += V_mean_freqs[j];
		}
	}

	delete[] V_mean_freqs;
	delete[] V_STD_freqs;

	V_mean_mean_Freq_not_zero = V_mean_mean_Freq / V_mean_mean_Freq_not_zero_count;
	V_mean_mean_Freq /= Node_count;

	fp0 = fopen("V_mean_mean_Freq.txt", "w+");
	fprintf(fp0, "%f\n", V_mean_mean_Freq);
	fclose(fp0);

	fp0 = fopen("V_mean_mean_Freq_not_zero.txt", "w+");
	fprintf(fp0, "%f\n", V_mean_mean_Freq_not_zero);
	fclose(fp0);

	// CHUNKS
	double chunk_t_start = 0.25 * t_max;
	double chunk_t_step = 0.5;
	int chunk_step_count = (0.75 * t_max / chunk_t_step);

	double dt_chunk = 0.1 / V_mean_mean_Freq_not_zero; // 0.25
	int chunks_count = chunk_t_step / dt_chunk;

	double corr_aver_mean = 0, corr_aver_mean_not_zero = 0;

	FILE* fp_corr_not_zero = fopen("corr_not_zero.txt", "w+");
	FILE* fp_V_spikes_chunks = fopen("V_spikes_chunks.txt", "w+");

	vector<double> k_syn;
	vector<double> k_syn_time;

	for (int s = 0; s < chunk_step_count; s++)
	{
		if (chunks_count > 0)
		{
			int** V_spikes_chunks = new int* [Node_count];

#ifdef USE_OMP
#pragma omp parallel for
#endif
			for (int i = 0; i < Node_count; i++)
			{
				V_spikes_chunks[i] = new int[chunks_count];

				if (V_spikes[i].size() == 0)
					for (int ch = 0; ch < chunks_count; ch++)
						V_spikes_chunks[i][ch] = 0;

				double ch_start = chunk_t_start + chunk_t_step * s;
				double ch_end = ch_start + dt_chunk;
				list<double>::iterator currentSpike = V_spikes[i].begin();

				for (int ch = 0; ch < chunks_count; ch++)
				{
					while (currentSpike != V_spikes[i].end() && *currentSpike < ch_start)
						currentSpike++;

					if (currentSpike == V_spikes[i].end())
						V_spikes_chunks[i][ch] = 0;
					else
						V_spikes_chunks[i][ch] = *currentSpike >= ch_start && *currentSpike <= ch_end;

					ch_start += dt_chunk;
					ch_end += dt_chunk;
				}
			}

			char buffer[50];
			//sprintf(buffer, "V_spikes_chunks_%d.txt", s);
			//fp0 = fopen(buffer, "w+");

			for (int ch = 0; ch < chunks_count; ch++)
			{
				for (int i = 0; i < Node_count; i++)
					fprintf(fp_V_spikes_chunks, "%d\t", V_spikes_chunks[i][ch]);

				fprintf(fp_V_spikes_chunks, "\n");
			}

			//fclose(fp0);

			double corr = 0;
			double corr_not_zero = 0;
			double counter = 0;
			double counter_not_zero = 0;

			for (int i = 0; i < Node_count; i++)
			{
				for (int j = 0; j < Node_count; j++)
				{
					if (i == j)
						continue;

					int k1 = 0;
					int k2 = 0;
					int k3 = 0;

					for (int l = 0; l < chunks_count; l++)
					{
						if (V_spikes_chunks[i][l] == 1 && V_spikes_chunks[j][l] == 1)
							k1++;

						k2 += V_spikes_chunks[i][l];
						k3 += V_spikes_chunks[j][l];
					}

					if (k2 != 0 && k3 != 0)
					{
						corr += (double)k1 / sqrt((double)k2 * (double)k3);
						counter_not_zero++;
					}

					counter++;
				}
			}

			double corr_aver = corr / counter;

			double corr_aver_not_zero = corr;

			if (counter_not_zero != 0)
				corr_aver_not_zero /= counter_not_zero;

			corr_aver_mean += corr_aver;
			corr_aver_mean_not_zero += corr_aver_not_zero;

			//sprintf(buffer, "corr_aver_%d.txt", s);

			//fp0 = fopen(buffer, "w+");
			//fprintf(fp0, "%f\n", corr_aver);
			//fclose(fp0);

			double ch_start = chunk_t_start + chunk_t_step * s;
			fprintf(fp_corr_not_zero, "%f\t%f\n", ch_start, corr_aver_not_zero);

			for (int i = 0; i < Node_count; i++)
				delete[] V_spikes_chunks[i];

			delete[] V_spikes_chunks;

			k_syn.push_back(corr_aver_not_zero);
			k_syn_time.push_back(ch_start);
		}
		else
		{
			double ch_start = chunk_t_start + chunk_t_step * s;
			fprintf(fp_corr_not_zero, "%f\t%f\n", ch_start, 0.0);
		}
	}

	fclose(fp_corr_not_zero);
	fclose(fp_V_spikes_chunks);

	corr_aver_mean /= chunk_step_count;
	corr_aver_mean_not_zero /= chunk_step_count;

	double k_syn_local_max_mean = 0;
	int k_syn_local_max_mean_count = 0;

	double k_syn_local_min_mean = 0;
	int k_syn_local_min_mean_count = 0;

	FILE* fp_max = fopen("k_syn_local_max.txt", "w+");
	FILE* fp_min = fopen("k_syn_local_min.txt", "w+");

	for (int i = 0; i < Ca_intervals.size(); i++)
	{
		if (Ca_intervals[i].T_end < 0)
			break;

		{
			double k_syn_local_max = 0;
			bool firstValue = true;

			for (int j = 0; j < k_syn.size(); j++)
			{
				if (k_syn_time[j] >= Ca_intervals[i].T_start && k_syn_time[j] <= Ca_intervals[i].T_end)
				{
					if (firstValue)
					{
						k_syn_local_max = k_syn[j];
						firstValue = false;
					}
					else if (k_syn[j] > k_syn_local_max)
					{
						k_syn_local_max = k_syn[j];
					}
				}
			}

			if (!firstValue)
			{
				k_syn_local_max_mean += k_syn_local_max;
				k_syn_local_max_mean_count++;

				fprintf(fp_max, "%f\n", k_syn_local_max);
			}
		}

		{
			double k_syn_local_min = 1;
			bool firstValue = true;

			for (int j = 0; j < k_syn.size(); j++)
			{
				if (k_syn_time[j] >= Ca_intervals[i].T_start && k_syn_time[j] <= Ca_intervals[i].T_end)
				{
					if (firstValue)
					{
						k_syn_local_min = k_syn[j];
						firstValue = false;
					}
					else if (k_syn[j] < k_syn_local_min)
					{
						k_syn_local_min = k_syn[j];
					}
				}
			}

			if (!firstValue)
			{
				k_syn_local_min_mean += k_syn_local_min;
				k_syn_local_min_mean_count++;

				fprintf(fp_min, "%f\n", k_syn_local_min);
			}
		}
	}

	fclose(fp_max);
	fclose(fp_min);

	if (k_syn_local_max_mean_count != 0)
		k_syn_local_max_mean /= k_syn_local_max_mean_count;

	if (k_syn_local_min_mean_count != 0)
		k_syn_local_min_mean /= k_syn_local_min_mean_count;

	//fp0 = fopen("corr_aver_mean.txt", "w+");
	//fprintf(fp0, "%f\n", corr_aver_mean);
	//fclose(fp0);

	//fp0 = fopen("corr_aver_mean_not_zero.txt", "w+");
	//fprintf(fp0, "%f\n", corr_aver_mean_not_zero);
	//fclose(fp0);
	////// CHUNKS END

	//2
	double* V_freq_sync_time = new double[Node_count];
	double* V_freq_sync_time_relative = new double[Node_count];

#ifdef USE_OMP
#pragma omp parallel for
#endif
	for (int i = 0; i < Node_count; i++)
	{
		list<double>::iterator it_V_spikes = V_spikes[i].begin();

		while (it_V_spikes != V_spikes[i].end() && *it_V_spikes < 0.25 * t_max)
		{
			V_spikes[i].pop_front();
			it_V_spikes = V_spikes[i].begin();
		}

		list<double> V_freqs_normalized;
		list<double> V_freqs_time;

		it_V_spikes = V_spikes[i].begin();

		if (V_spikes[i].size() <= 1)
		{
			continue;
		}
		else
		{
			for (int j = 1; j < V_spikes[i].size(); j++)
			{
				double first = *it_V_spikes;
				advance(it_V_spikes, 1);
				double next = *it_V_spikes;

				double T = next - first;
				V_freqs_normalized.push_back(1 / T - V_mean_mean_Freq);
				V_freqs_time.push_back(next);
			}
		}

		list<double>::iterator it_Freq_normalized = V_freqs_normalized.begin();
		list<double>::iterator it_Freq_time = V_freqs_time.begin();

		V_freq_sync_time[i] = 0;

		for (int j = 1; j < V_freqs_normalized.size(); j++)
		{
			double Freq_normalized_last = *it_Freq_normalized;
			advance(it_Freq_normalized, 1);
			double Freq_normalized_next = *it_Freq_normalized;

			double Freq_time_last = *it_Freq_time;
			advance(it_Freq_time, 1);
			double Freq_time_next = *it_Freq_time;

			if (abs(Freq_normalized_last) <= 0.5 && abs(Freq_normalized_next) <= 0.5 && (Freq_time_next - Freq_time_last) <= 0.035)
				V_freq_sync_time[i] += Freq_time_next - Freq_time_last;
		}

		V_freq_sync_time_relative[i] = V_freq_sync_time[i] / (t_max - (0.25 * t_max));
	}

	double V_mean_freq_sync_time_relative = 0;

	for (int j = 0; j < Node_count; j++)
	{
		V_mean_freq_sync_time_relative += V_freq_sync_time_relative[j];
	}

	V_mean_freq_sync_time_relative /= Node_count;

	fp0 = fopen("V_freq_sync_time_relative.txt", "w+");
	for (int i = 0; i < Node_count; i++)
	{
		fprintf(fp0, "%f\t", V_freq_sync_time_relative[i]);
	}
	fclose(fp0);

	fp0 = fopen("V_mean_freq_sync_time_relative.txt", "w+");
	fprintf(fp0, "%f\n", V_mean_freq_sync_time_relative);
	//fprintf(fp_res, "%f\n", V_mean_freq_sync_time_relative);
	fclose(fp0);

	if (write_file)
	{
		//fclose(fp_res);
		fclose(fp_I_stim);
		///fclose(fp_I_syn);
		fclose(fp_Ca);
		//fclose(fp_IP3);
		//fclose(fp_z);
		fclose(fp_V);
		//fclose(fp_m);
		//fclose(fp_n);
		//fclose(fp_h);
		//fclose(fp_G_P);
		fclose(fp_V_P);
		//fclose(fp_m_P);
		//fclose(fp_n_P);
		//fclose(fp_h_P);

		fclose(fp_V_spikes);
	}

	double extime_rk4 = 0;

	if (omp_enable)
	{
#ifdef USE_OMP
		omp_end_rk4 = omp_get_wtime();
		extime_rk4 = (double)(omp_end_rk4 - omp_start_rk4);// / CLOCKS_PER_SEC;
#endif
	}
	else
	{
		end_rk4 = clock();
		extime_rk4 = (double)(end_rk4 - start_rk4);// / CLOCKS_PER_SEC;
	}

	int minutes = (int)extime_rk4 / 60;
	int seconds = (int)extime_rk4 % 60;
	printf("\nExecution time is: %d minutes %d seconds\n ", minutes, seconds);

	/*int nth;
	#pragma omp parallel
	{
		#pragma omp master
		nth = omp_get_num_threads();
	}*/

	fp0 = fopen("time_exec.txt", "a");
	//fprintf(fp0, "%d %lf\n", nth, extime_rk4);
	fprintf(fp0, "%lf\n", extime_rk4);
	fclose(fp0);

	fp0 = fopen("results_g_syn_P_g_syn_k_syn.txt", "a");
	fprintf(fp0, "%f\t%f\t%f\t%f\t%f\n", g_syn_P, g_syn, k_syn_local_min_mean, corr_aver_mean_not_zero, k_syn_local_max_mean);
	fclose(fp0);

	return 0;

	for (int i = 0; i < Node_count; i++)
		delete[] A_A[i];

	delete[] A_A;

	for (int i = 0; i < Node_count; i++)
		delete[] B_A[i];

	delete[] B_A;

	delete[] C_A;

	for (int i = 0; i < Node_count; i++)
		delete[] A_N[i];

	delete[] A_N;

	for (int i = 0; i < Node_count; i++)
		delete[] B_N[i];

	delete[] B_N;

	delete[] C_N;

	delete[] A_N_P;

	for (int i = 0; i < Node_count; i++)
		delete[] B_N_P[i];

	delete[] B_N_P;

	delete[] C_N_P;

	delete[] f;
	delete[] f_diff;
	delete[] v_4;
	delete[] E_syn;
	delete[] I_app;
	delete[] V_spikes;
	delete[] V_spikes_Freq;

	delete[] Meander_start_from_zero;
	delete[] Meander_width;
	delete[] Meander_height;
	delete[] Meander_interval;
	delete[] last_meander_end;
	delete[] V_freq_sync_time;

	for (int i = 0; i < Node_count; i++)
		delete[] tau[i];

	delete[] tau;

	for (int i = 0; i < Node_count; i++)
		delete[] V_old_array[i];

	delete[] V_old_array;

	for (int i = 0; i < Equations_count; i++)
		delete[] k[i];

	delete[] k;

	delete[] phi_k1;
	delete[] phi_k2;
	delete[] phi_k3;
}
