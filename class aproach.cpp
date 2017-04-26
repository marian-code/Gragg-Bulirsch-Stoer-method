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

class Modified_midpoint
{
public:

	double n_steps;
	double H_big, h_small;
	double time;
	double y1, y2;
	vector<double> z1, z2;

	double function_f1(double t, double y1, double y2)
	{
		return y2;
	}

	double function_f2(double t, double y1, double y2)
	{
		double omega = 2 * 3.14159;
		return -omega*omega*sin(y1);
	}

	void cycle()
	{
		int m;

		//itegrovat so zjemnenym krokom na intervale od z(1) -> z(n)
		for (m = 1; m < n_steps; m++)
		{
			z1.push_back(z1[m] + h_small*function_f1(time * m, z1[m], z2[m]));
			z2.push_back(z2[m] + h_small*function_f2(time * m, z1[m], z2[m]));
			cout << "ok" << endl;
		}
	}

	void new_step()
	{
		//spocitat dalsi krok -> y(i + H) pomocou z(n) a z(n - 1)
		y1 = (z1[n_steps] + z1[n_steps - 1] + h_small*function_f1(time + H_big, z1[n_steps], z2[n_steps])) / 2;
		y2 = (z2[n_steps] + z2[n_steps - 1] + h_small*function_f2(time + H_big, z1[n_steps], z2[n_steps])) / 2;
	}

};



int main()
{
	Modified_midpoint MP;

	int i, m;
	double H_big, h_small;
	double n_steps;
	vector<double> time;
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
		////kôli tomu aby sa indexy nemuseli posuvat v cykle a bolo to prehladnejsie
		//z1.push_back(0);
		//z2.push_back(0);

		////priradit hodnoty na kraji intervalu y(i) -> z(0)
		//z1.push_back(y1[i]);
		//z2.push_back(y2[i]);

		MP.time = time[i];
		MP.H_big = H_big;
		MP.h_small = h_small;

		MP.z1.resize(2);
		MP.z1[0] = 0;
		MP.z1[1] = y1[i];

		MP.z2.resize(2);
		MP.z2[0] = 0;
		MP.z2[1] = y2[i];

		MP.cycle();
		MP.new_step();

		////itegrovat so zjemnenym krokom na intervale od z(1) -> z(n)
		//for (m = 1; m < n_steps; m++)
		//{
		//	z1.push_back(z1[m] + h_small*function_f1(time[i] * m, z1[m], z2[m]));
		//	z2.push_back(z2[m] + h_small*function_f2(time[i] * m, z1[m], z2[m]));
		//}

		////spocitat dalsi krok -> y(i + H) pomocou z(n) a z(n - 1)
		//y1.push_back((z1[n_steps] + z1[n_steps - 1] + h_small*function_f1(time[i] + H_big, z1[n_steps], z2[n_steps])) / 2);
		//y2.push_back((z2[n_steps] + z2[n_steps - 1] + h_small*function_f2(time[i] + H_big, z1[n_steps], z2[n_steps])) / 2);

		y1.push_back(MP.y1);
		y2.push_back(MP.y2);

		//posunut sa do dalsieho casoveho kroku
		time.push_back(time[i] + H_big);

		////zmazat hodnoty na podintervaloch z(i)
		//z1.clear();
		//z2.clear();

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

