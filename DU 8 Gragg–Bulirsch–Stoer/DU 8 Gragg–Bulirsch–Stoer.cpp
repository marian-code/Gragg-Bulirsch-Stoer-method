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

vector<double> modified_midpoint(int i, double time, double n_steps, double h_small, double H_big, double y1, double y2 )
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
	int i, m;
	double H_big, h_small;
	double n_steps;
	vector<double> time, temp;
	vector<double> y1, y2;
	vector<double> z1, z2;
	vector<vector<double>> Richardson;

	const double start = 0;
	const double end = 1.5;

	//modified midpoint method

	n_steps = 12;//pocet podintevalov v jednom kroku H
	H_big = 0.0005;//velkost kroku
	h_small = H_big / n_steps;

	i = 0;
	time.push_back(0);//zaciatok casu

	y1.push_back(0);//pocitocna podmienka x(0)
	y2.push_back(2);//pocitocna podmienka x´(0)

	while (time[i] < end)
	{
		
		//spocitat dalsi krok pomocou modified midpoint rule -> y(i + H) pomocou z(n) a z(n - 1)
		temp = modified_midpoint(i, time[i], n_steps, h_small, H_big, y1[i], y2[i]);//ak by som priradoval postupne tak by sa funkcia musela vykonat 2x
		y1.push_back(temp[0]);
		y2.push_back(temp[1]);

		//posunut sa do dalsieho casoveho kroku
		time.push_back(time[i] + H_big);;

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

