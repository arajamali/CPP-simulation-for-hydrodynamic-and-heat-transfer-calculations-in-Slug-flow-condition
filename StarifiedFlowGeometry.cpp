#include "StratifiedFlowGeometry.h"
#include "CaseSpecificInputs.h"
#include <iostream>
#include <cmath>

namespace MultiphaseTemperature {
	//In this constructor please note that we want to calculate liquid height (hL).
	//h is the liquid holdup. One should make sure to know the differences between the two
	//h=AL/AP, by knowing this relationship, we can calculate AL and then calculate hL
	StratifiedFlowGeometry::StratifiedFlowGeometry(double H,double IDEffective) {
		user::PipeGeometry Pipe;
		ID = IDEffective;
		double AP = 0.25*ID * ID*Pi;
		double hh = H * ID;
		double AL = H * AP;
		double hLTi;
		double Eps = 10;
		double dx = 0.000000001;
		double F, FPlusdx, FPrime;
		while (Eps > 0.00000001)
		{
			F = ID * ID*(0.25*(Pi - acos(2 * (hh / ID) - 1) + (2 * (hh / ID) - 1)*pow(1 - pow(2 * (hh / ID) - 1, 2), 0.5))) - H * AP;
			FPlusdx = ID * ID*(0.25*(Pi - acos(2 * ((hh + dx) / ID) - 1) + (2 * ((hh + dx) / ID) - 1)*pow(1 - pow(2 * ((hh + dx) / ID) - 1, 2), 0.5))) - H * AP;
			FPrime = ((FPlusdx)-(F)) / dx;
			hh = hh - F / FPrime;
			Eps = F;
		}
		H = hh;
		hL = H;

		A = Pi * ID*ID*0.25;
		hLTilt = hL / ID;
		SGTilt = acos(2 * hLTilt - 1);
		SLTilt = Pi - acos(2 * hLTilt - 1);
		SITilt = pow(1 - pow(2 * hLTilt - 1, 2), 0.5);
		AGTilt = 0.25*(acos(2 * hLTilt - 1) - SITilt * (2 * hLTilt - 1));
		ALTilt = 0.25*(Pi - acos(2 * hLTilt - 1) + (2 * hLTilt - 1)*pow(1 - pow(2 * hLTilt - 1, 2), 0.5));
		dL = 4.0 * ALTilt*ID*ID / (SLTilt*ID);
		dG = 4.0 * AGTilt*ID*ID / (SITilt*ID + SGTilt * ID);


		if (hL < 0.5*ID)
		{
			Theta = 2 * acos((0.5*ID - hL) / (0.5*ID));
		}
		else if (hL == 0.5*ID) {
			Theta = Pi;
		}
		else {
			Theta = Pi + 2 * asin((hL - 0.5*ID) / (0.5*ID));
		}

	}
}//namespace MultiphaseTemperature {