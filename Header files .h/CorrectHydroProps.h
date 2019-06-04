#ifndef CorrectHydroProps_H
#define CorrectHydroProps_H
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "CaseSpecificInputs.h"
#include "StratifiedFlowGeometry.h"
#include "SlugFlowFilmProfile.h"
#include "dPCalculation.h"
#include "SlugFlowHeatTransfer.h"

namespace MultiphaseTemperature
{
	class SPModel 
	{
	public:
		std::vector<std::vector<std::vector<std::vector<double>>>>TwWithHydroUpdate;
		std::vector<std::vector<double>> TSWithHydroUpdate;
		std::vector<std::vector<double>> TLWithHydroUpdate;
		std::vector<std::vector<double>> TGWithHydroUpdate;
		std::vector<std::vector<double>> hHot;
		std::vector<std::vector<std::vector<double>>> KDepo;
		//This vector contains the r_ref*U_ref parameter for the over heat transfer calculation
		std::vector<double> ru;

		/*In this method, the temperature distribution is estimated using
		the correct hydrodynamic update after one T_u*/
		double TForHydroUpdate(std::vector<double>& SlugFLowHeatTransferCalInputs, 
                           std::vector<std::vector<double>>& WaxVector,
			                     std::vector<std::vector<double>>& WaxSolidFraction,
			                     std::vector<double>& HydroInputs, 
                           std::vector<std::vector<std::vector<double>>>& TwInitial,
			                     std::vector<double>& Ts, std::vector<double>& Tf, 
                           std::vector<double>& Tg,
                           std::vector<std::vector<double>> SlugFilmAssign, 
                           std::vector<double>& q, int Nx, double Tm);

		void SlugHeatTransferWithHydroUpdate(std::vector<double>& SlugFLowHeatTransferCalInputs, 
                                         std::vector<std::vector<double>>& WaxVector,
			                                   std::vector<std::vector<double>>& WaxSolidFraction,
			                                   std::vector<double>& HydroInputs, 
                                         std::vector<std::vector<std::vector<double>>>& TwInitial,
			                                   std::vector<double>& Ts, std::vector<double>& Tf, 
                                         std::vector<double>& Tg, 
                                         std::vector<std::vector<double>> SlugFilmAssign,
			                                   std::vector<double>& q, int Nx, int NTime,  double Tm);

		void RrefUref(std::vector<double>& SlugFLowHeatTransferCalInputs, 
                  std::vector<std::vector<double>>& WaxVector,
			            std::vector<std::vector<double>>& WaxSolidFraction,
		            	std::vector<double>& HydroInputs, 
                  std::vector<std::vector<std::vector<double>>>& TwInitial,
			            std::vector<double>& Ts, std::vector<double>& Tf,
                  std::vector<double>& Tg, std::vector<std::vector<double>> SlugFilmAssign,
		            	std::vector<double>& q, int Nx, int NTime, double Tm);

		double FrictionFactorAnnulus(double Rho, double V, double dpOut, double danIn, double Mu);
		double NuxGn(double Rho, double v, double dpOut, double danIn,
			           double Mu, double Cp, double K, double dx);
	};

} //End of namespace MultiphaseTemperature
#endif
