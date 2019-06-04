#ifndef dPCalculation_H
#define dPCalculation_H
#include <iostream>
#include <cmath>
#include <vector>
#include "CaseSpecificInputs.h"
#include "SlugFlowFilmProfile.h"
#include "StratifiedFlowGeometry.h"

namespace MultiphaseTemperature {
	class PressureCal {
	public:
		PressureCal(std::vector<double>& Inputs, std::vector<double>& HydroInputs);
		double dP;
		double ShearS, PressSlug /*Pressure contribtion from slug*/, fS /*friction factor*/;
		double ShearF, ShearSF, PressSFilm, PressFilm, fF,fSF;
		double ShearG,ShearSG, PressSGas, PressGas, fG, fSG;
	private: 
		double VSL, VSG, ID, T, P, PipeRoughness, IT;
		double Vm;
		double RhoG, RhoL, MuL, MuG, RhoS, MuS;
		double HLLS, Lf, Ls,HLTBAvg;
		double AF, AG, SG, SF, SI;
		double VF, VG, VTb, VLLS, VGLS, VGTB,VLTB;
		double LfInitial, dx, Epsilon, hMax, EpsilonZh;

		//Tells the number of points in RK45 discritization
		int ZSize;
		std::vector<double> zArray, HLTBArray;
		//There is one shear associated with the slug body section

	};

}//End of MultiphaseTemperature namespace

#endif
