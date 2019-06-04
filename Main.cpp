#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>
#include <cmath>
#include <vector>
#include <iomanip> 
#include "CaseSpecificInputs.h"
#include "SlugFlowFilmProfile.h"
#include "SlugFlowHeatTransfer.h"
#include "dPCalculation.h"
#include "OveralHeatTransfer.h"
#include <omp.h>
//namespace MultiphaseTemperature {


int main() {
	double a;
	user::PipeGeometry Pipe;
	int Nx=50, NTime=1000;
	double T = 29+273.15;
	double P = (14.7+350)*6894.76;
	double VSG =1*0.305515;
	double VSL =4*0.305515;
	double ID = Pipe.IDPipe; //ID should be an input parameter for slugfilm but not for heat transfer
	double Roughness = 0.00004572;
	double IT = 0.025;
	double TInitial = 302;
	double InitialGuess =0.1, dx = 0.01, Epsilon = 0.01, hMax =0.1, EpsilonZh = 0.001;
	user::FluidProperties Fluid(T, P);
	std::vector<double> FilmInputs;
	std::vector<double> HydroInputs;
	std::vector<double> HeatInputs;

	std::vector<std::vector<double>> WaxThickness;
	std::vector<std::vector<double>> WaxFraction;

	std::vector<double> q;
	std::vector<double> TS;
	std::vector<double> TL;
	std::vector<double> TG;

	std::vector<bool> SlugFlow;
	std::vector<double> HLTbLocation;

	//Initial condition assignment
	std::vector<std::vector<std::vector<double>>>Tw;
	Tw.resize(Nx + 1);
	for (int xx = 0; xx < Nx; xx++) //To include the initial condtion as contant temperature
	{
		Tw[xx].resize(4 + 1);
		for (int WL = 1; WL < 4 + 1; WL++) //starts from 1 to 3
		{
			Tw[xx][WL].resize(4 + 1);
			for (int RS = 1; RS < 4 + 1; RS++) //starts from 1 to 4
			{
				//***************************Possible modifications are requried
				Tw[xx][WL][RS] = TInitial;  /*It is assumed that inlet temperature is constant and it is equal to T at t=0 everywhere.
											Later this assumption should be substitute with pipe-in-pipe counter courent flow */
			}
		}

		TL.push_back(TInitial);
		TG.push_back(TInitial);
		TS.push_back(TInitial);
		q.push_back(8.0);
	} 


	WaxThickness.resize(Nx + 1);
	WaxFraction.resize(Nx + 1);
	for (int xx = 0; xx < Nx; xx++)
	{
		WaxFraction[xx].resize(5);
		WaxFraction[xx][1] = 0.001/* + double(xx) / double(Nx) * 0.5*0.04*/;
		WaxFraction[xx][2] = 0.001/*WaxFraction[xx][1] - WaxFraction[xx][1] * 0.1*/;
		WaxFraction[xx][3] = 0.001/* WaxFraction[xx][1] - WaxFraction[xx][1] * 0.2*/;
		WaxFraction[xx][4] = WaxFraction[xx][2];

		WaxThickness[xx].resize(5);
		WaxThickness[xx][1] = 0.0005 /*+ double(xx) / double (Nx) * 0.5*0.001*/;
		WaxThickness[xx][2] = WaxThickness[xx][1] /*- WaxThickness[xx][1] * 0.1*/;
		WaxThickness[xx][3]= WaxThickness[xx][1] /*- WaxThickness[xx][1] * 0.2*/;
		WaxThickness[xx][4] = WaxThickness[xx][2];

	}
	
	//FilmInputs, WaxWidthVector,  WaxSolidFraction, HydroInputs, TwInitial,
	//	Ts,  Tf,  Tg,  FilmAssign,
	//	 SlugAssign,  q, NxInput, Tm, TInlet


		FilmInputs.clear();
		HydroInputs.clear();


		FilmInputs.push_back(VSL); FilmInputs.push_back(VSG); FilmInputs.push_back(ID); FilmInputs.push_back(T-273.15); //T has to be in C
		FilmInputs.push_back(P); FilmInputs.push_back(Roughness); FilmInputs.push_back(IT);

		HydroInputs.push_back(InitialGuess); HydroInputs.push_back(dx); HydroInputs.push_back(Epsilon); HydroInputs.push_back(hMax); HydroInputs.push_back(EpsilonZh);
		double TCorr;
	//MultiphaseTemperature::PressureCal Press(FilmInputs, HydroInputs); 

		MultiphaseTemperature::SlugFLowHeatTransferCal OBJ;

		OBJ.HydroUpdate(FilmInputs, HydroInputs);
		SlugFlow.resize(Nx+1);
		HLTbLocation.resize(Nx + 1);
		for (int xx = 0; xx < Nx; xx++)
		{
			SlugFlow[xx] = false;
			HLTbLocation[xx] = false;
			if (xx <= (OBJ.NxF - 1))
			{
				SlugFlow[xx] = false;
				HLTbLocation[xx] = OBJ.HLTBOrgNew[xx];

			}
			if (xx >(OBJ.NxF - 1))
			{
				SlugFlow[xx] = true;
			}
			//s1 = HLTbLocation[xx]; s2 = SlugFlow[xx];
		}
		
		double ch;
		MultiphaseTemperature::SlugFlowOveralHeatTranfer obj2;
		ch=obj2.AvgTempForHydroUpdate(FilmInputs, WaxThickness, WaxFraction, HydroInputs, Tw, TS, TL, TG, HLTbLocation, SlugFlow, q, 302.15,0.001);
		
		//std::vector<std::vector<double>> Save;
		//Save.resize(1000);
		//for (int j = 0; j < 1000; j++) {
		//	Save[j].resize(Nx + 1);
		//}

		//for (int i = 0; i < 1000; i++)
		//{
		//	OBJ.SlugFlowHeatModel(FilmInputs, WaxThickness, WaxFraction, HydroInputs, Tw, TS, TL, TG, HLTbLocation, SlugFlow, q, 302.15);
		//	Tw = OBJ.Tw;
		//	TS = OBJ.TS;
		//	TL = OBJ.TL;
		//	TG = OBJ.TG;
		//	OBJ.SlugFlowBooleanLater(SlugFlow, HLTbLocation);
		//	SlugFlow = OBJ.SlugFlow;
		//	HLTbLocation= OBJ.HLTbLocation;
		//	Save[i] = OBJ.TGActual;
		//	std::cout << OBJ.SumhCheck << " " << OBJ.NuNumCal << std::endl;

		//}
		

	//std::cout << VSG << " " << VSL << " " << Film.LfFinal<<" "<<Film.HLLS<<" "<<Film.HLTBAvg<<" "<<Film.err<< std::endl;
		//for (int i = 0; i < 90; i++)
		//{
		//	for (int j = 0; j < Nx; j++)
		//	{
		//		std::cout<<std::setprecision(16) << Save[i][j] << " ";
		//	}
		//	std::cout << std::endl;
		//}
	

	//MultiphaseTemperature::SlugFLowHeatTransferCal SFHTC(FilmInputs);
	//VS,VL,ID,Temp,Pressure (FilmProfile)
	auto start_time = std::chrono::high_resolution_clock::now();
	//MultiphaseTemperature::FilmProfile Film(FilmInputs);
	//Film.HyrodynamicProperties(InitialGuess, dx, Epsilon, hMax);
	
	//std::cout << Press.dP;
	auto end_time = std::chrono::high_resolution_clock::now();
	std::cout << "The total program run-time in milliseconds is ";
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(
		end_time - start_time).count() << '\n';

	

	
	//for (int i = 0; i < Film.ZSize+1; i++)
	//{
	//	std::cout << Film.XSave[i] << " " << Film.HLTBSave[i] << std::endl;
	//}
	std::system("pause");

	//a=Film.hFDetermination(10);
	//a = Film.LfDetermination();
	//cout << Film.hFDetermination(a);

	//std::cin >> a;





}
