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
#include "OveralHeatTransfer.h"

namespace MultiphaseTemperature {

	/*								IMPORTANT NOTE
	In this method, the objective is to find correct hydrodynamic properties
	for t=Tu. The hydrodynamic properties are calcualted based on the 
	thermophysical properties of liquid and gas and they are heavility dependent
	on temperature. However, the temperature anywhere in the slug unit changes
	with time and with location. Therefore, thermophysical properties change too.
	In addition to thermophysical properties, wax characteristics also change
	including the thickness and ultimately the effective pipe's radius which is
	very important for hydrodynamic property estimation. This inconsostency needed
	to get fixed. In the following method, we find the most appropeiate temperature
	for hydrodynamic property calculation. This iterative process is to check Lf
	with temperature. Since wax deposition calculation is costly (timely), we 
	did not included the wax thickness update. In other word, wax characteristics
	at tIniial is considered. This assumption should note jepordize the accuracy much*/ 

	/*The following method was checked for many flow rates. While following method's
	importance might not be obvious at the first seen but it is!!
	For VSG=10 and VSL=1, the average temperautre calculated to be ALOT different than the
	intiial Temperature */
	double SlugFlowOveralHeatTranfer::AvgTempForHydroUpdate(std::vector<double>& FilmInputs,
		std::vector<std::vector<double>>& WaxWidthVector, std::vector<std::vector<double>>& WaxSolidFraction,
		std::vector<double>& HydroInputs, std::vector<std::vector<std::vector<double>>>& TwInitial,
		std::vector<double>& Ts, std::vector<double>& Tf, std::vector<double>& Tg, std::vector<double>& FilmAssign,
		std::vector<bool>& SlugAssign, std::vector<double>& q, double TInlet, double EpsLimit)
	{
		
		double Lf1=100, Lf2;
		double Eps =100;
		std::vector<std::vector<double>> Save;
		Save.resize(1000);
		for (int j = 0; j < 1000; j++) {
			Save[j].resize(50 + 1);
		}
		double AvgTb= 0;
		int count = 0;
		
		//t=Tu 
		while (Eps > EpsLimit)
		{
			//Counter of the number of iteration
			count = count + 1;
			AvgTb = 0;
			MultiphaseTemperature::SlugFLowHeatTransferCal OBJ;
			Lf2 = Lf1;
			OBJ.HydroUpdate(FilmInputs, HydroInputs);
			Lf1 = OBJ.Lf;
			for (int t = 0; t < OBJ.Nx; t++)
			{
				OBJ.SlugFlowHeatModel(FilmInputs, WaxWidthVector, WaxSolidFraction, HydroInputs, TwInitial, Ts, Tf, Tg, FilmAssign, SlugAssign, q, TInlet);
				AvgTb = AvgTb + OBJ.AvgTb;
				TwInitial = OBJ.Tw;
				Ts = OBJ.TS;
				Tf = OBJ.TL;
				Tg = OBJ.TG;
				OBJ.SlugFlowBooleanLater(SlugAssign, FilmAssign);
				SlugAssign = OBJ.SlugFlow;
				FilmAssign = OBJ.HLTbLocation;
			}
			std::cout << OBJ.NuNumCal << " " << OBJ.SumhCheck<<std::endl;
			AvgTb = AvgTb / double(OBJ.Nx);
			FilmInputs[3] = AvgTb - 273.15;
			Eps = abs(Lf1 - Lf2);
		}
		return(AvgTb);
	}

	/* From the previous method, the AvgTb temperature for the correct
	calcualtion of hydrodynamic properties is done. This temperature can
	be used now for further calcualtion and larger time
	*/


} //end of MultiphaseTemperature namespace