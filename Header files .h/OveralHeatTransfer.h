#ifndef OverlHeatTransfer_H
#define OverlHeatTransfer_H
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "CaseSpecificInputs.h"
#include "StratifiedFlowGeometry.h"
#include "SlugFlowFilmProfile.h"
#include "dPCalculation.h"
#include "SlugFlowHeatTransfer.h"
#include "OveralHeatTransfer.h"
namespace MultiphaseTemperature {

	//In this class, we will try to calucalte the temperature distribution in counter-current flow
	class SlugFlowOveralHeatTranfer {
	public:
	
		//In this method, the correct average temperature for hydrodynamic property update will be calculated
		double AvgTempForHydroUpdate(std::vector<double>& FilmInputs,
			std::vector<std::vector<double>>& WaxWidthVector, std::vector<std::vector<double>>& WaxSolidFraction,
			std::vector<double>& HydroInputs, std::vector<std::vector<std::vector<double>>>& TwInitial,
			std::vector<double>& Ts, std::vector<double>& Tf, std::vector<double>& Tg, std::vector<double>& FilmAssign,
			std::vector<bool>& SlugAssign, std::vector<double>& q, double TInlet, double EpsLimit);
	};

} //end of MultiphaseTemperature namespace
#endif