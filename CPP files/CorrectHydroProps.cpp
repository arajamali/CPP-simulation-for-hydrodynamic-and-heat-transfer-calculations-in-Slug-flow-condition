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
	double SPModel::TForHydroUpdate(std::vector<double>& SlugFLowHeatTransferCalInputs, 
                                  std::vector<std::vector<double>>& WaxVector,
		                              std::vector<std::vector<double>>& WaxSolidFraction,
		                              std::vector<double>& HydroInputs, 
                                  std::vector<std::vector<std::vector<double>>>& TwInitial,
	                              	std::vector<double>& Ts, std::vector<double>& Tf,  
                                  std::vector<double>& Tg, std::vector<double>& q, 
                                  int Nx, double Tm)
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

	void SPModel::SlugHeatTransferWithHydroUpdate(std::vector<double>& SlugFLowHeatTransferCalInputs, 
                                                std::vector<std::vector<double>>& WaxVector,
		                                            std::vector<std::vector<double>>& WaxSolidFraction,
		                                            std::vector<double>& HydroInputs, 
                                                std::vector<std::vector<std::vector<double>>>& TwInitial,
		                                            std::vector<double>& Ts, std::vector<double>& Tf, 
                                                std::vector<double>& Tg, std::vector<double>& q,
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

	}

	//Gnielinski(2009), friction factor in the annulus
	double SPModel::FrictionFactorAnnulus(double Rho, double V, double dpOut, 
                                        double danIn, double Mu) {
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

		ans = (friction / 8.0*Pr*Re)*(1 + pow(Dhyd / dx, 2.0 / 3.0))*(Fan) /
          (k1 + 12.7*pow(friction / 8.0, 0.5)*(pow(Pr, 2.0 / 3.0) - 1));
		return(ans);
	}

	void SPModel::RrefUref(std::vector<double>& SlugFLowHeatTransferCalInputs, 
                         std::vector<std::vector<double>>& WaxVector,
		                     std::vector<std::vector<double>>& WaxSolidFraction,
		                     std::vector<double>& HydroInputs, 
                         std::vector<std::vector<std::vector<double>>>& TwInitial,
		                     std::vector<double>& Ts, std::vector<double>& Tf, 
                         std::vector<double>& Tg, std::vector<double>& q,
		                     int Nx, int NTime, double Tm)
		{
		//By this method, the slug flow heat tranfser calcualtion is performed and 
		//TS, TG, TL, Tw and nusslet number is calcuated with correct 
		SlugHeatTransferWithHydroUpdate(SlugFLowHeatTransferCalInputs, WaxVector,
			WaxSolidFraction, HydroInputs, TwInitial, Ts, Tf, Tg, q, Nx, NTime, Tm);

		}

} //end of namespace MultiphaseTemperature












