// DU 8 Gragg–Bulirsch–Stoer.cpp : Defines the entry point for the console application.
//
//prevzate z http://web.mit.edu/ehliu/Public/Spring2006/18.304/implementation_bulirsch_stoer.pdf
//a z numerical recepies

#include "stdafx.h"
#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>

using namespace std;

const double pi = 3.14159265358979323;

double function_f1(double t, double y1, double y2)
{
	return y2;
}

double function_f2(double t, double y1, double y2)
{
	double omega = 1.0 / sqrt(2);
	return -omega*omega*sin(y1);
}

vector<double> modified_midpoint(int i, double time, int n_steps, double h_small, double H_big, double y1, double y2)
{
	int m;
	vector<double> z1, z2;
	vector<double> y;

	//priradit hodnoty na kraji intervalu y(i) -> z(0)
	z1.push_back(y1);
	z2.push_back(y2);

	z1.push_back(z1[0] + h_small*function_f1(time, z1[0], z2[0]));
	z2.push_back(z2[0] + h_small*function_f2(time, z1[0], z2[0]));

	//itegrovat so zjemnenym krokom na intervale od z(1) -> z(n)
	for (m = 1; m < n_steps; m++)
	{
		z1.push_back(z1[m -1] + 2*h_small*function_f1(time * m, z1[m], z2[m]));
		z2.push_back(z2[m -1] + 2*h_small*function_f2(time * m, z1[m], z2[m]));
	}

	//spocitat dalsi krok -> y(i + H) pomocou z(n) a z(n - 1)
	y.push_back((z1[n_steps] + z1[n_steps - 1] + h_small*function_f1(time + H_big, z1[n_steps], z2[n_steps])) / 2);
	y.push_back((z2[n_steps] + z2[n_steps - 1] + h_small*function_f2(time + H_big, z1[n_steps], z2[n_steps])) / 2);

	return y;
}

vector<double> TDMA(const int N, vector<vector<double>> matrix, vector<double> vec) //pre tridiagonalne matice (Thomas algorithm)
{
	int i;
	vector<double> x;
	double help;
	x.resize(N);

	for (i = 1; i < N; i++)
	{
		help = matrix[i][i - 1] / matrix[i - 1][i - 1];
		matrix[i][i - 1] -= matrix[i - 1][i - 1] * help;
		matrix[i][i] -= matrix[i - 1][i] * help;;
		vec[i] -= vec[i - 1] * help;
	}

	x[N - 1] = vec[N - 1] / matrix[N - 1][N - 1];

	for (i = N - 2; i > -1; i--)
		x[i] = (vec[i] - x[i + 1] * matrix[i][i + 1]) / matrix[i][i];

	return x;
}

double interpolation(vector< double > y1, vector< double > time, int iter, const int interpolation_points)
{
	int i, j, k;
	double x;
	double spline;
	double min_value = 100;
	double quarter_period;
	vector< double > h_i, y_i, a_i, sigma, u_i;
	vector<vector< double >> matrix;

	cout.precision(13);
	cout.setf(std::ios::fixed, std::ios::floatfield);

	//nastavit velkost poli
	a_i.resize(interpolation_points);
	y_i.resize(interpolation_points);
	h_i.resize(interpolation_points - 1);
	u_i.resize(interpolation_points - 2);
	sigma.resize(interpolation_points - 2);

	matrix.resize(interpolation_points - 2);
	for (i = 0; i < interpolation_points - 2; i++)
		matrix[i].resize(interpolation_points - 2);

	//ekvidistantne delenie
	for (i = 0; i < interpolation_points; i++)
	{
		a_i[i] = time[iter - interpolation_points + i + 1];
		y_i[i] = y1[iter - interpolation_points + i + 1];
	}

	//pocitanie dlzky intervalov
	for (i = 0; i < interpolation_points - 1; i++)
		h_i[i] = a_i[i + 1] - a_i[i];

	//kontrolne vypisy
	cout << "a" << endl;
	for (i = 0; i < interpolation_points; i++)
	cout << a_i[i] << " ";
	cout << endl;

	cout << "y" << endl;
	for (i = 0; i < interpolation_points; i++)
	cout << y_i[i] << " ";
	cout << endl;

	cout << "h" << endl;
	for (i = 0; i < interpolation_points - 1; i++)
	cout << h_i[i] << " ";
	cout << endl;

	cout << "sigma" << endl;
	for (i = 0; i < interpolation_points - 2; i++)
	cout << sigma[i] << " ";
	cout << endl;

	//naplnanie matice
	for (i = 0; i < interpolation_points - 2; i++)
		matrix[i][i] = 2 * (h_i[i] + h_i[i + 1]);

	for (i = 0; i < interpolation_points - 3; i++)
		matrix[i][i + 1] = h_i[i];

	for (i = 1; i < interpolation_points - 2; i++)
		matrix[i][i - 1] = h_i[i];

	//naplnanie pravej strany
	for (i = 1; i < interpolation_points - 1; i++)
		u_i[i - 1] = (y_i[i + 1] - y_i[i]) / h_i[i] - (y_i[i] - y_i[i - 1]) / h_i[i - 1];

	//kontrolny vypis
	/*cout << "matica" << endl;
	for (i = 0; i < interpolation_points - 2; i++)
	{
	for (j = 0; j < interpolation_points - 2; j++)
	cout << matrix[i][j] << " ";
	cout << "\t" << u_i[i] << endl;
	}*/

	//riesenie sustavy trojdiagonalnym algoritmom
	sigma.swap(TDMA(interpolation_points - 2, matrix, u_i));

	//doplnit nuly na zaciatok a koniec vektora sigma
	sigma.insert(sigma.begin(), 0);
	sigma.push_back(0);

	//kontrolny vypis
	/*cout << "sigma" << endl;
	for (i = 0; i < interpolation_points; i++)
	cout << sigma[i] << " ";
	cout << endl;*/

	//vykreslenie aproximacie a najdenie bodu najblizsie k nule
	ofstream outfile2;
	outfile2.open("interpolation.dat");
	outfile2.precision(13);
	outfile2.setf(std::ios::fixed, std::ios::floatfield);
	for (i = 0; i < interpolation_points - 1; i++)
	{
		for (j = 0; j < 10000; j++)
		{
			x = a_i[i] + j*(h_i[i]) / 9999;
			spline = (pow(x - a_i[i], 3)*sigma[i + 1]) / h_i[i]
				+ (pow(a_i[i + 1] - x, 3)*sigma[i]) / h_i[i]
				+ ((y_i[i] - pow(h_i[i], 2)*sigma[i])*(a_i[i + 1] - x)) / h_i[i]
				+ ((y_i[i + 1] - pow(h_i[i], 2)*sigma[i + 1])* (x - a_i[i])) / h_i[i];

			outfile2 << x << " " << spline << endl;

			//hladanie bodu kde funkcia pretina 0
			if (abs(spline) < min_value)
			{
				min_value = spline;
				quarter_period = x;
			}
		}
	}
	outfile2.close();

	return quarter_period;
}

int main()
{
	bool repeat = true;
	unsigned int i, j, k, m;
	unsigned int n_steps;
	int interpolation_points;
	double H_big, h_small;
	double local_error;
	vector<double> scale;
	vector<double> H_evolution;
	vector<double> time, temp;
	vector<double> delta;
	vector<double> y1, y2;
	vector<double> z1, z2;
	vector<vector<double>> Richardson_1, Richardson_2;

	//const double start = 0;
	//const double end = 13;
	const double epsilon = 1E-14;
	const double max_n_steps = 8;
	const double H_max = 0.00005;
	const double absolute_tolerance = epsilon;
	const double relative_tolerance = epsilon;

	Richardson_1.resize(max_n_steps);
	for (i = 0; i <	max_n_steps; i++)
		Richardson_1[i].resize(i + 1);

	Richardson_2.resize(max_n_steps);
	for (i = 0; i < max_n_steps; i++)
		Richardson_2[i].resize(i + 1);

	delta.resize(2);
	scale.resize(2);

	cout.precision(13);
	cout.setf(std::ios::fixed, std::ios::floatfield);

#pragma region kontrolny vypis

	//cout << "initial richardson coefficients" << endl;

	////kontrolny vypis
	//cout << "Richardson 1" << endl;
	//for (j = 0; j < max_n_steps; j++)
	//{
	//	for (m = 0; m < j; m++)
	//	{
	//		cout << Richardson_1[j][m] << "\t";
	//	}
	//	cout << endl;
	//}
	//cout << endl;
	//cout << "Richardson 2" << endl;
	//for (j = 0; j < max_n_steps; j++)
	//{
	//	for (m = 0; m < j; m++)
	//	{
	//		cout << Richardson_2[j][m] << "\t";
	//	}
	//	cout << endl;
	//}
	//cout << endl;

	//cout << "repeat " << repeat << endl;

#pragma endregion

	H_big = 0.00005;//velkost kroku

	i = 0;
	time.push_back(0);//zaciatok casu

	y1.push_back(2 * pi / 3);//pocitocna podmienka x(0)
	y2.push_back(0);//pocitocna podmienka x´(0)

	//while (time[i] < end)//integrovanie cez urcity cas
	while(y1[i] > -2*H_big)// hladanie periody
	{
		//cout << "----------------------- step " << i + 1 << "-----------------------------------" << endl;
		//cout << "stepsize H: " << H_big << endl;

		//max_n_steps urcuje max pocet podintervalov na intervale H
		for (k = 0; k < max_n_steps; k++)
		{

			n_steps = 2 * (k + 1); //pocet podintevalov v jednom kroku H
			h_small = H_big / n_steps;

			//cout << "# of substeps: " << n_steps << endl;

			//spocitat dalsi krok pomocou modified midpoint rule -> y(i + H) pomocou z(n) a z(n - 1)
			temp = modified_midpoint(i, time[i], n_steps, h_small, H_big, y1[i], y2[i]);//ak by som priradoval postupne tak by sa funkcia musela vykonat 2x
			Richardson_1[k][0] = temp[0];
			Richardson_2[k][0] = temp[1];

			//cout << "temp " << temp[0] << "\t" << temp[1] << endl;

			//ak niesom v prvom kroku tak pomocou Richardsona dopocitat aproximacie vyssich radov
			if (k > 0)
			{
				for (j = 0; j < k; j++)
				{
					//cout << "nsteps " << n_steps << " / " << (2 * (k - j + 1)) << endl;
					//cout << "menovatel " << (pow(double(n_steps) / (2 * (k - j)), 2) - 1) << endl;
					Richardson_1[k][j + 1] = Richardson_1[k][j] + (Richardson_1[k][j] - Richardson_1[k - 1][j]) / (pow(double(n_steps) / (2 * (k - j)), 2) - 1);
					Richardson_2[k][j + 1] = Richardson_2[k][j] + (Richardson_2[k][j] - Richardson_2[k - 1][j]) / (pow(double(n_steps) / (2 * (k - j)), 2) - 1);
				}

				//spocitat rozdiely aproximacii
				delta[0] = abs(Richardson_1[k][k] - Richardson_1[k][k - 1]);
				delta[1] = abs(Richardson_2[k][k] - Richardson_2[k][k - 1]);

				scale[0] = absolute_tolerance + abs(Richardson_1[k][k])*relative_tolerance;
				scale[1] = absolute_tolerance + abs(Richardson_2[k][k])*relative_tolerance;

				local_error = 0;
				for (m = 0; m < 2; m++)
					local_error += pow(delta[m] / scale[m], 2);

				local_error = sqrt(local_error / 2);

				//cout << "local_error" << local_error << endl;

				//ked sa dosiahne pozadovana lokalna presnost -> break
				if (local_error <= 1)
				{
					//cout << "**********************************  break  *************************************" << endl;
					repeat = false;
					break;
				}
			}

#pragma region kontrolne vypisy

			//kontrolny vypis
			/*cout << "Richardson 1" << endl;
			for (j = 0; j < k; j++)
			{
				for (m = 0; m < j + 1; m++)
				{
					if (Richardson_1[j][m] != 0)
						cout << Richardson_1[j][m] << "\t";
				}
				cout << endl;
			}
			cout << endl;
			cout << "Richardson 2" << endl;
			for (j = 0; j < k; j++)
			{
				for (m = 0; m < j + 1; m++)
				{
					if (Richardson_2[j][m] != 0)
						cout << Richardson_2[j][m] << "\t";
				}
				cout << endl;
			}
			cout << endl;*/

#pragma endregion

		}
		k--; //k treba o 1 zmensit lebo po skonceni cyklu sa este pricita 1 naviac

		//cout << "repeat " << repeat << endl;

		// repeat = false ked sa v danom kroku nepodarilo pomocou richardsonovej extrapolacie
		//skonvergovat na pozadovanu presnost a je potrebne opakovat krok s inym H
		if (repeat == false)
		{
			//priradit novy krok do trajektorie
			y1.push_back(Richardson_1[k][k]);
			y2.push_back(Richardson_2[k][k]);

			//posunut sa do dalsieho casoveho kroku
			time.push_back(time[i] + H_big);

			i++;
		}
		repeat = true;

		//spocitat velkost noveho kroku H - adaptive stepsize control
		H_big = H_big*0.98*pow(0.98 / local_error, 1.0 / (2 * k + 1));

		if (H_big > H_max) H_big = H_max;

		H_evolution.push_back(H_big);

		//cout << "----------------------- end of step " << i + 1 << "------------------------------" << endl;
	}

	//zistit kolko cisel za 0 sa naslo pravdepodobne cca 3
	interpolation_points = 0;
	while (y1[y1.size() - interpolation_points - 1] < 0)
	{
		interpolation_points++;
	}
	
	if (interpolation_points < 3) cout << "warning!! to few interpolation points" << endl;

	cout << "perioda kyvadla je: " << 4 * interpolation(y1, time, i, interpolation_points * 2) << endl;

#pragma region zápis na disk

	ofstream outfile1; //zapis do suboru
	outfile1.precision(13);
	outfile1.setf(std::ios::fixed, std::ios::floatfield);

	outfile1.open("x(t).dat");
	for (i = 0; i < y1.size(); i++)
	{
		outfile1 << time[i] << "\t";
		outfile1 << y1[i] << endl;
	}
	outfile1.close();

	ofstream outfile2; //zapis do suboru
	outfile2.precision(13);
	outfile2.setf(std::ios::fixed, std::ios::floatfield);

	outfile2.open("v(t).dat");
	for (i = 0; i < y2.size(); i++)
	{
		outfile2 << time[i] << "\t";
		outfile2 << y2[i] << endl;
	}
	outfile2.close();

	ofstream outfile3; //zapis do suboru
	outfile3.precision(13);
	outfile3.setf(std::ios::fixed, std::ios::floatfield);

	outfile3.open("x(v).dat");
	for (i = 1; i < y1.size(); i++)
	{

		outfile3 << y1[i] << "\t";
		outfile3 << y2[i] << endl;
	}
	outfile3.close();

	ofstream outfile4; //zapis do suboru
	outfile4.precision(13);
	outfile4.setf(std::ios::fixed, std::ios::floatfield);

	outfile4.open("H(t).dat");
	for (i = 1; i < y1.size(); i++)
	{

		outfile4 << time[i] << "\t";
		outfile4 << H_evolution[i] << endl;
	}
	outfile4.close();

#pragma endregion

	system("PAUSE");

	return 0;
}

