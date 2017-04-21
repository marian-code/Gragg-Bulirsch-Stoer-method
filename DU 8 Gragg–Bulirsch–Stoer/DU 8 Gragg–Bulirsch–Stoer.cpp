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

using namespace std;

double function_f1(double t, double y1, double y2)
{
	return y2;
}

double function_f2(double t, double y1, double y2)
{
	double omega = 2 * 3.14159;
	return -omega*omega*sin(y1);
}

vector<double> modified_midpoint(int i, double time, int n_steps, double h_small, double H_big, double y1, double y2 )
{
	int m;
	vector<double> z1, z2;
	vector<double> y;

	//kôli tomu aby sa indexy nemuseli posuvat v cykle a bolo to prehladnejsie
	z1.push_back(0);
	z2.push_back(0);

	//priradit hodnoty na kraji intervalu y(i) -> z(0)
	z1.push_back(y1);
	z2.push_back(y2);

	//itegrovat so zjemnenym krokom na intervale od z(1) -> z(n)
	for (m = 1; m < n_steps; m++)
	{
		z1.push_back(z1[m] + h_small*function_f1(time * m, z1[m], z2[m]));
		z2.push_back(z2[m] + h_small*function_f2(time * m, z1[m], z2[m]));
	}

	//spocitat dalsi krok -> y(i + H) pomocou z(n) a z(n - 1)
	y.push_back((z1[n_steps] + z1[n_steps - 1] + h_small*function_f1(time + H_big, z1[n_steps], z2[n_steps])) / 2);
	y.push_back((z2[n_steps] + z2[n_steps - 1] + h_small*function_f2(time + H_big, z1[n_steps], z2[n_steps])) / 2);

	return y;
}

int main()
{
	unsigned int i, j, k, m;
	double H_big, h_small;
	int n_steps;
	vector<double> time, temp;
	vector<double> local_error;
	vector<double> y1, y2;
	vector<double> z1, z2;
	vector<vector<double>> Richardson_1, Richardson_2;

	const double start = 0;
	const double end = 1.5;
	const double epsilon = 0.0001;
	const double max_n_steps = 8;

	Richardson_1.resize(8);
	for (i = 0; i < 8; i++)
		Richardson_1[i].resize(i + 1);

	Richardson_2.resize(8);
	for (i = 0; i < 8; i++)
		Richardson_2[i].resize(i + 1);

	cout.precision(5);
	cout.setf(std::ios::fixed, std::ios::floatfield);

#pragma region kontrolny vypis

	cout << "initial richardson coefficients" << endl;

	//kontrolny vypis
	cout << "Richardson 1" << endl;
	for (j = 0; j < max_n_steps; j++)
	{
		for (m = 0; m < j; m++)
		{
			cout << Richardson_1[j][m] << "\t";
		}
		cout << endl;
	}
	cout << endl;

	//kontrolny vypis
	cout << "Richardson 2" << endl;
	for (j = 0; j < max_n_steps; j++)
	{
		for (m = 0; m < j; m++)
		{
			cout << Richardson_2[j][m] << "\t";
		}
		cout << endl;
	}
	cout << endl;

#pragma endregion
	
	H_big = 0.005;//velkost kroku

	i = 0;
	time.push_back(0);//zaciatok casu

	y1.push_back(0);//pocitocna podmienka x(0)
	y2.push_back(2);//pocitocna podmienka x´(0)

	while (time[i] < end)
	{
		cout << "----------------------- step " << i + 1 << "-----------------------------------" << endl;

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

			cout << "temp " << temp[0] << "\t" << temp[1] << endl;

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
				local_error.push_back(abs(Richardson_1[k][k] - Richardson_1[k - 1][k - 1]));
				local_error.push_back(abs(Richardson_2[k][k] - Richardson_2[k - 1][k - 1]));

				//najst najväèšiu chybu
				sort(local_error.begin(), local_error.end());

				//ked sa dosiahne pozadovana lokalna presnost -> break
				if (local_error[local_error.size() - 1] < epsilon)
				{
					cout << "break-------------------------" << endl;
					break;
				}
			}

#pragma region kontrolne vypissy

			//kontrolny vypis
			cout << "Richardson 1" << endl;
			for (j = 0; j < max_n_steps; j++)
			{
				for (m = 0; m < j + 1; m++)
				{
					if(Richardson_1[j][m] != 0)
						cout << Richardson_1[j][m] << "\t";
				}
				cout << endl;
			}
			cout << endl;
			cout << "Richardson 2" << endl;
			for (j = 0; j < max_n_steps; j++)
			{
				for (m = 0; m < j + 1; m++)
				{
					if(Richardson_2[j][m] != 0)
						cout << Richardson_2[j][m] << "\t";
				}
				cout << endl;
			}
			cout << endl;

#pragma endregion

		}
		k--; //k treba o 1 zmensit lebo po skonceni cyklu sa este pricita 1 naviac
		
		//priradit novy krok do trajektorie
		y1.push_back(Richardson_1[k][k]);
		y2.push_back(Richardson_2[k][k]);

		//posunut sa do dalsieho casoveho kroku
		time.push_back(time[i] + H_big);;

		//if (i > 100) break;

		cout << "----------------------- end of step " << i + 1 << "------------------------------" << endl;

		i++;
	}

	ofstream outfile1; //zapis do suboru

	//nastavenie presnosti vypisu
	outfile1.precision(13);
	outfile1.setf(std::ios::fixed, std::ios::floatfield);

	outfile1.open("graf.dat");
	for (i = 0; i < y1.size(); i++)
	{
		outfile1 << time[i] << "\t";
		outfile1 << y1[i] << endl;
	}
	outfile1.close();

	ofstream outfile2; //zapis do suboru

	//nastavenie presnosti vypisu
	outfile2.precision(13);
	outfile2.setf(std::ios::fixed, std::ios::floatfield);

	outfile2.open("graf1.dat");
	for (i = 0; i < y2.size(); i++)
	{
		outfile2 << time[i] << "\t";
		outfile2 << y2[i] << endl;
	}
	outfile2.close();

	ofstream outfile3; //zapis do suboru

	//nastavenie presnosti vypisu
	outfile3.precision(13);
	outfile3.setf(std::ios::fixed, std::ios::floatfield);

	outfile3.open("graf2.dat");
	for (i = 1; i < y1.size(); i++)
	{

		outfile3 << y1[i] << "\t";
		outfile3 << y2[i] << endl;
	}
	outfile3.close();

	system("PAUSE");

    return 0;
}

