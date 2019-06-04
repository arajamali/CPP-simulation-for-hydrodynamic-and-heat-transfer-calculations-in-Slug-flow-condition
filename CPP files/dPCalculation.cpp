#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <vector>
#include "CaseSpecificInputs.h"
#include "SlugFlowFilmProfile.h"
#include "StratifiedFlowGeometry.h"
#include "dPCalculation.h"

namespace MultiphaseTemperature {
	PressureCal::PressureCal(std::vector<double>& Inputs, std::vector<double>& HydroInputs)
	{

		VSL = Inputs[0];
		VSG = Inputs[1];
		ID = Inputs[2];
		T = Inputs[3];
		P = Inputs[4];
		PipeRoughness = Inputs[5];
		IT = Inputs[6];//InterfactialTension

		LfInitial = HydroInputs[0];
		dx = HydroInputs[1];
		Epsilon = HydroInputs[2];
		hMax = HydroInputs[3];
		EpsilonZh = HydroInputs[4];
		zArray.clear();
		HLTBArray.clear();
		user::PipeGeometry PipeProps;
		user::FluidProperties Fluid(T, P);
		FilmProfile Film(Inputs);

		//Fluid Properties
		RhoG = Fluid.Gas.Rho(Fluid); 
		RhoL = Fluid.Liquid.Rho(Fluid);
		MuL = Fluid.Liquid.Mu(Fluid);
		MuG = Fluid.Gas.Mu(Fluid);
		//Now lets call the hydrodynamic properties
		
		Film.HyrodynamicProperties(LfInitial, dx, Epsilon, hMax, EpsilonZh);
		//Slug section properties
		HLLS=Film.HLLS;
		ZSize = Film.ZSize;
		RhoS = RhoL * HLLS + RhoG * (1 - HLLS);
		MuS = MuL * HLLS + MuG * (1 - HLLS);
		Lf = Film.LfFinal;
		Ls= Film.Ls;
		Vm = Film.Vm;
		HLTBAvg= Film.HLTBAvg;
		VTb = Film.VTb;
		VLLS = Film.VLLS;
		VGLS = Film.VGLS;

		//Assigning HLTB in z direction calculated by RKF45
		for (int i = 0; i < ZSize + 1; i++)
		{
			zArray.push_back(Film.XSave[i]);
			HLTBArray.push_back(Film.HLTBSave[i]);
		}
		fS = Film.FannningFrictionFactor(RhoS, Vm, ID, MuS, PipeRoughness);
		ShearS = fS * RhoS*Vm*Vm / 2.0;
		PressSlug = (ShearS * M_PI*ID / (M_PI*0.25*ID*ID))*Ls;
		double z;
		double HLTB;
		PressFilm = 0;
		PressGas = 0;
		PressSFilm = 0;
		PressSGas = 0;
		for (int i = 1; i < ZSize + 1; i++) //Check this for 2
		{
			z= (zArray[i] - zArray[i-1]);
			HLTB=(HLTBArray[i] + HLTBArray[i-1])/2.0;
			StratifiedFlowGeometry Geometry(HLTB,ID);
			AF = Geometry.ALTilt*ID*ID;
			AG = Geometry.AGTilt*ID*ID;
			SG = Geometry.SGTilt*ID;
			SF = Geometry.SLTilt*ID;
			SI = Geometry.SITilt*ID;
			VF = (VTb - VLLS)*HLLS / HLTB;
			VG = (VTb - VGLS)*(1 - HLLS) / (1 - HLTB);
			VLTB= VTb - VF;
			VGTB= VTb - VG;
			fF = Film.FannningFrictionFactor(RhoL, abs(VLTB), 
                                       4 * AF / SF, MuL, PipeRoughness);
			fG = Film.FannningFrictionFactor(RhoG, abs(VGTB), 
                                       4 * AG / (SI + SG), MuG, PipeRoughness);

			ShearF = fF * RhoL*abs(VLTB)*VLTB / 2.0;
			ShearG = fG * RhoG*abs(VGTB)*VGTB / 2.0;
			PressFilm = PressFilm+ShearF * SF*z / (M_PI*0.25*ID*ID);
			PressGas = PressGas+ShearG * SG*z / (M_PI*0.25*ID*ID);

			//For hibiki X-Parameter, these factors are needed
			fSF = Film.FannningFrictionFactor(RhoL, VSL, M_PI*ID*ID*0.25, MuL, PipeRoughness);
			fSG = Film.FannningFrictionFactor(RhoG, VSG, M_PI*ID*ID*0.25, MuG, PipeRoughness);
			ShearSF = fSF * RhoL*VSL*VSL / 2.0;
			ShearSG = fSG * RhoG*VSG*VSG / 2.0;
			PressSFilm = PressSFilm + ShearSF * (SF + SF)*z / (M_PI*0.25*ID*ID);
			PressSGas = PressSGas + ShearSG * (SF + SF)*z / (M_PI*0.25*ID*ID);
		}

		dP = (PressGas + PressFilm + PressSlug)*0.000145038;
	}

}//End of MultiphaseTemperature namespace
