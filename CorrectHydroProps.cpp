#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <iomanip> 
#include "CaseSpecificInputs.h"
#include "StratifiedFlowGeometry.h"
#include "SlugFlowFilmProfile.h"
#include "SlugFlowHeatTransfer.h"
#include "dPCalculation.h"
#include "CorrectHydroProps.h"

namespace MultiphaseTemperature
{
	double SPModel::TForHydroUpdate(std::vector<double>& SlugFLowHeatTransferCalInputs, std::vector<std::vector<double>>& WaxVector,
		std::vector<std::vector<double>>& WaxSolidFraction,
		std::vector<double>& HydroInputs, std::vector<std::vector<std::vector<double>>>& TwInitial,
		std::vector<double>& Ts, std::vector<double>& Tf,  std::vector<double>& Tg, std::vector<double>& q, int Nx, double Tm)
	{

		double LfPre, LfCurr = 1000;
		double Temp=302;
		double AvgTBK = Tm;
		double EpsLf;
		
		for (int kL = 0; kL < 100; kL++)
		{
			//Hydrodynamic properties calculation
			LfPre = LfCurr;
			SlugFLowHeatTransferCal SFHT;
			SFHT.SlugFlowHeatModel(SlugFLowHeatTransferCalInputs, WaxVector, WaxSolidFraction,
				HydroInputs,TwInitial,Ts,Tf,Tg, q,Nx, Nx, Temp);
			Temp = SFHT.AvgTBK;
			SFHT.HydroUpdate(HydroInputs, Temp);
			LfCurr = SFHT.Lf;
			EpsLf = abs(LfCurr - LfPre);

			if (EpsLf < 0.001) {
				break;
			}
		}
		return(Temp);
	}

	void SPModel::SlugHeatTransferWithHydroUpdate(std::vector<double>& SlugFLowHeatTransferCalInputs, std::vector<std::vector<double>>& WaxVector,
		std::vector<std::vector<double>>& WaxSolidFraction,
		std::vector<double>& HydroInputs, std::vector<std::vector<std::vector<double>>>& TwInitial,
		std::vector<double>& Ts, std::vector<double>& Tf, std::vector<double>& Tg, std::vector<double>& q,
		int Nx, int NTime, double Tm)
	{ 
		TwWithHydroUpdate.clear();
		TSWithHydroUpdate.clear();
		TGWithHydroUpdate.clear();
		TLWithHydroUpdate.clear();
		KDepo.clear();
		hHot.clear();

			double CorrectTemp;
			CorrectTemp = TForHydroUpdate(SlugFLowHeatTransferCalInputs, WaxVector, WaxSolidFraction, HydroInputs,
				TwInitial, Ts, Tf, Tg, q, Nx, Tm);

			SlugFLowHeatTransferCal SFHT;
			SFHT.SlugFlowHeatModel(SlugFLowHeatTransferCalInputs, WaxVector, WaxSolidFraction,
				HydroInputs, TwInitial, Ts, Tf, Tg, q, Nx, Nx, CorrectTemp);
			TwWithHydroUpdate = SFHT.TwSave;
			TSWithHydroUpdate = SFHT.TSSave;
			TGWithHydroUpdate = SFHT.TGSave;
			TLWithHydroUpdate = SFHT.TLSave;
			KDepo=SFHT.KDepoSave;
			hHot=SFHT.hHotSave;
			TwInitial = SFHT.Tw;
			Ts = SFHT.TS;
			Tf = SFHT.TL;
			Tg = SFHT.TG;


		//TForHydroUpdate(SlugFLowHeatTransferCalInputs, WaxVector,WaxSolidFraction, HydroInputs, Nx, Tm)
	}

	//Gnielinski(2009), friction factor in the annulus
	double SPModel::FrictionFactorAnnulus(double Rho, double V, double dpOut, double danIn, double Mu) {
		double Dhyd = danIn - dpOut;
		double Re,ReStar, a;
		double ans; 
		a = dpOut / danIn;
		Re = Rho * V*Dhyd / Mu;
		ReStar = Re * ((1 + a * a)*log(a) + 1 - a * a) / ((1 - a * a)*log(a));
		ans = pow(1.8*log10(ReStar) - 1.5, -2.0);
		return(ans);
	} 

	//k is assumed to be 1, Not accurate
	double SPModel::NuxGn(double Rho, double v, double dpOut, double danIn,
		double Mu, double Cp, double K, double dx)
	{
		double ans;
		double Dhyd = danIn - dpOut;
		double Re = Rho * v*Dhyd / Mu;
		double Pr = Mu * Cp / K;
		double a = dpOut / danIn;
		double Fan = 0.75*pow(a, -0.17);
		double k1; 
		double friction;
		k1 = 1.07 + 900.0 / Re - 0.63 / (1 + 10 * Pr);
		friction = FrictionFactorAnnulus(Rho, v,dpOut,danIn,Mu);

		ans = (friction / 8.0*Pr*Re)*(1 + pow(Dhyd / dx, 2.0 / 3.0))*(Fan) / (k1 + 12.7*pow(friction / 8.0, 0.5)*(pow(Pr, 2.0 / 3.0) - 1));
		return(ans);
	}

	void SPModel::RrefUref(std::vector<double>& SlugFLowHeatTransferCalInputs, std::vector<std::vector<double>>& WaxVector,
		std::vector<std::vector<double>>& WaxSolidFraction,
		std::vector<double>& HydroInputs, std::vector<std::vector<std::vector<double>>>& TwInitial,
		std::vector<double>& Ts, std::vector<double>& Tf, std::vector<double>& Tg, std::vector<double>& q,
		int Nx, int NTime, double Tm)
		{
		//By this method, the slug flow heat tranfser calcualtion is performed and 
		//TS, TG, TL, Tw and nusslet number is calcuated with correct 
		SlugHeatTransferWithHydroUpdate(SlugFLowHeatTransferCalInputs, WaxVector,
			WaxSolidFraction, HydroInputs, TwInitial, Ts, Tf, Tg, q, Nx, NTime, Tm);


		}

} //end of namespace MultiphaseTemperature
















	//{
	//std::ofstream TLOut, TGOut, TW1Out, TW2Out, TW3Out, NuCal, NuCheck,
	//	TW4Out, TW5Out, TW6Out, TW7Out, TW8Out, TWax1, TWax2, TWax3,
	//	TW9Out, Time, Axial, Tbulk;
	//TLOut.open("TL.txt");
	//TGOut.open("TG.txt");
	//TWax1.open("TWax1.txt");
	//TWax2.open("TWax2.txt");
	//TWax3.open("TWax3.txt");
	//TW1Out.open("TW1.txt");
	//TW2Out.open("TW2.txt");
	//TW3Out.open("TW3.txt");
	//TW4Out.open("TW4.txt");
	//TW5Out.open("TW5.txt");
	//TW6Out.open("TW6.txt");
	//TW7Out.open("TW7.txt");
	//TW8Out.open("TW8.txt");
	//TW9Out.open("TW9.txt");
	//Tbulk.open("TBulk.txt");
	//NuCal.open("NuCal.txt");
	//NuCheck.open("NuCheck.txt");
	//Time.open("Time.txt");
	//Axial.open("Axial.txt"); //Axial discritization
	//bool flag = 0;

	//for (int Count = 0; Count < 20; Count++)
	//{
	//	double CorrectTemp;
	//	CorrectTemp = TForHydroUpdate(SlugFLowHeatTransferCalInputs, WaxVector, WaxSolidFraction, HydroInputs, TwInitial, Ts, Tf, Tg, Nx, Tm);

	//	SlugFLowHeatTransferCal SFHT;
	//	SFHT.SlugFlowHeatModel(SlugFLowHeatTransferCalInputs, WaxVector, WaxSolidFraction,
	//		HydroInputs, TwInitial, Ts, Tf, Tg, Nx, Nx, CorrectTemp);
	//	TwInitial = SFHT.Tw;
	//	Ts = SFHT.TS;
	//	Tf = SFHT.TL;
	//	Tg = SFHT.TG;

	//	if (flag == 0) {
	//		for (int xx = 0; xx < Nx; xx++) {
	//			Axial << SFHT.dX*double(xx) << " ";
	//		}
	//		flag = 1;
	//	}
	//	for (int tt = 0; tt < NTime; tt++) {
	//		Time << double(Count)*NTime*SFHT.dt + SFHT.dt*tt << " ";
	//		for (int xx = 0; xx < Nx; xx++) {
	//			TLOut << std::setprecision(16) << SFHT.TLSave[tt][xx] << " ";
	//			TGOut << std::setprecision(16) << SFHT.TGSave[tt][xx] << " ";
	//			Tbulk << std::setprecision(16) << SFHT.TSSave[tt][xx] << " ";
	//			TWax1 << std::setprecision(16) << SFHT.TwSave[tt][xx][1][1] << " ";
	//			TWax2 << std::setprecision(16) << SFHT.TwSave[tt][xx][1][2] << " ";
	//			TWax3 << std::setprecision(16) << SFHT.TwSave[tt][xx][1][3] << " ";
	//			TW1Out << std::setprecision(16) << SFHT.TwSave[tt][xx][2][1] << " ";
	//			TW2Out << std::setprecision(16) << SFHT.TwSave[tt][xx][2][2] << " ";
	//			TW3Out << std::setprecision(16) << SFHT.TwSave[tt][xx][2][3] << " ";
	//			TW4Out << std::setprecision(16) << SFHT.TwSave[tt][xx][3][1] << " ";
	//			TW5Out << std::setprecision(16) << SFHT.TwSave[tt][xx][3][2] << " ";
	//			TW6Out << std::setprecision(16) << SFHT.TwSave[tt][xx][3][3] << " ";
	//			TW7Out << std::setprecision(16) << SFHT.TwSave[tt][xx][4][1] << " ";
	//			TW8Out << std::setprecision(16) << SFHT.TwSave[tt][xx][4][2] << " ";
	//			TW9Out << std::setprecision(16) << SFHT.TwSave[tt][xx][4][3] << " ";

	//		}
	//		TWax1 << std::endl;
	//		TWax2 << std::endl;
	//		TWax3 << std::endl;
	//		NuCheck << std::endl;
	//		NuCal << std::endl;
	//		TLOut << std::endl;
	//		TGOut << std::endl;
	//		Tbulk << std::endl;
	//		TW1Out << std::endl;
	//		TW2Out << std::endl;
	//		TW3Out << std::endl;
	//		TW4Out << std::endl;
	//		TW5Out << std::endl;
	//		TW6Out << std::endl;
	//		TW7Out << std::endl;
	//		TW8Out << std::endl;
	//		TW9Out << std::endl;
	//		if (Count == 19) {
	//			std::cout << SFHT.NuNumCheck << " " << SFHT.NuNumCal << std::endl;
	//		}
	//	}
	//}
	////TForHydroUpdate(SlugFLowHeatTransferCalInputs, WaxVector,WaxSolidFraction, HydroInputs, Nx, Tm)
	//}