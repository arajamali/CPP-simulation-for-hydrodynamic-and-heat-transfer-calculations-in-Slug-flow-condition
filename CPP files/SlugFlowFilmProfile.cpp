//In this source file, slug flow hydrodynamic properties are calculated
//with full film profile calculation. In this source file, mainly Taitel and Barnea
//model is developed. 
//In this source file, all units are in SI units
#include <iostream>
#include <chrono>
#include <cmath>
#include "CaseSpecificInputs.h"
#include "SlugFlowFilmProfile.h"
#include <fstream>


namespace MultiphaseTemperature {

	//Final results for Film profile is calculated in "HyrodynamicProperties" method	
	FilmProfile::FilmProfile(std::vector<double>& FilmProfileInputs) {

		//???????????????
		Flag2 = 0;
		g = 9.8;
		VSL = FilmProfileInputs[0];
		VSG = FilmProfileInputs[1];
		ID = FilmProfileInputs[2];
		T = FilmProfileInputs[3];
		P = FilmProfileInputs[4];
		PipeRoughness = FilmProfileInputs[5];
		IT = FilmProfileInputs[6];//InterfactialTension
		//Mixture velocity
		Vm = VSL + VSG;
		//Distribution variable
		C0 = 1.2;

		user::PipeGeometry PipeProps;
		OD = PipeProps.ODPipe;

		user::FluidProperties Fluid(T, P);
		//This value is slightly overestimated in compare with what Ake 
		//reported in EXCEL file. The answer will be noticably different
		//if the values from Ake EXCEL file is used
		RhoG = Fluid.Gas.Rho(Fluid); // or RhoG = Fluid.Gas.Rho2(Fluid) (for P=350Psi)

		//From fitting correlations
		RhoL = Fluid.Liquid.Rho(Fluid);
		MuL = Fluid.Liquid.Mu(Fluid);
		MuG = Fluid.Gas.Mu(Fluid);

		// Gregory - suggested by Zheng 2003 for intial gueess
		HLLS = 1 / (1 + pow(Vm / 8.66, 1.39));

		/*********Some other HLLS relations***********/

		//HLLS = exp(-(2.48*pow(10, -6)*ReLS));
		/*HLLS = 1 - 0.058*(2 * pow(0.4*IT / 
     ((RhoL - RhoG)*9.8), 0.5)*pow(RhoL / 
     IT, 0.6)*pow(2.0/ID* 0.007*pow(VSG+VSL,3), 0.4) - 0.725, 2); */



		//Mixture properties
		RhoS = RhoL * HLLS + RhoG * (1 - HLLS);
		MuS = MuL * HLLS + MuG * (1 - HLLS);

		//From (Pan, 2010), the average slug length
		//for oil/Air experiments was reported to be 24D for Mu = 4 Cp 
		Ls = 30 * ID;

		//Translational velocity
		VTb = C0 * (VSL + VSG) + 0.54*pow(g*ID, 0.5);

		//Velocity of gas and liquid phase in slug section
		VGLS = C0 * (VSL + VSG);
		VLLS = ((VSL + VSG) - VGLS * (1 - HLLS)) / HLLS;
		Flag2 = 0;
	}

	//Fanning friction factor (main relation-implicit)
	double FilmProfile::FannningFrictionFactor(double Rho, double V, 
                                             double ID, double Mu, 
                                             double Roughness) {
		double Eps = 10, FPlusdx, FPrime;
		double F, Re;
		double ans = 0.0001;
		double dx = 1E-10;
		while (Eps > 1E-15)
		{
			Re = Rho * V*ID / Mu;
			F = 1 / pow(ans, 0.5) + 4 * log10(((Roughness / ID) / (3.7)) + (1.256 / (Re*pow(ans, 0.5))));
			FPlusdx = 1 / pow((ans + dx), 0.5) + 4 * log10(((Roughness / ID) / (3.7)) + (1.256 / (Re*pow((ans + dx), 0.5))));
			FPrime = ((FPlusdx)-(F)) / dx;
			ans = ans - F / FPrime;
			Eps = F;
		}
		return(ans);
	}

	//ODE of dhF/dz is calculated in this method
	double FilmProfile::hFdZ(double hF) {

		//Stratified Geometric parameters 
		double SG, SF, SI;
		double AF, AG;

		//Shear terms from liqiud and gas to wall and from gas to liquid interface
		double fF, fG, fI;
		double ShearF, ShearG, ShearI;


		double checking;
		//Relative velocities to VTb
		double VF, VG;
		//other variables
		double hLTilt = hF / ID;//This is not liquid phase holdup
		double HL; //This is liquid holdup
		double A;
		double dHLTBdhF;// dHLTB/dhF
		double ans;

		// H should be liquid holdup!
		//At this point, we need to calculate the liquid holdup
		//Becase we need liquid holdup (not liquid height) for stratified geometry class
		hLTilt = hF / ID;

		//Now we can calculate the liquid phase holdup
		HL = (1 / Pi)*(Pi - acos(2 * hF / ID - 1) + (2 * hF / ID - 1)*
			pow(1 - pow(2 * hF / ID - 1, 2), 0.5));

		//Parameters of stratified flow gemoetry are calculated for the given lqidui holdup
		StratifiedFlowGeometry Geometry(HL,ID);
		AF = Geometry.ALTilt*ID*ID;
		AG = Geometry.AGTilt*ID*ID;
		SG = Geometry.SGTilt*ID;
		SF = Geometry.SLTilt*ID;
		SI = Geometry.SITilt*ID;

		//Relative velocities (referneced to VTb)
		VF = (VTb - VLLS)*HLLS / HL;
		VG = (VTb - VGLS)*(1 - HLLS) / (1 - HL);
		
		//Actual velocities of film and gas core in stratified region
		VLTB = VTb - VF;
		VGTB = VTb - VG;

		dHLTBdhF = (4.0 / (Pi*ID))*pow(1 - pow(2 * hF / ID - 1, 2), 0.5);

		//Friction factors 
		fF = FannningFrictionFactor(RhoL, abs(VLTB), 4 * AF / SF, MuL, PipeRoughness);
		fG = FannningFrictionFactor(RhoG, abs(VGTB), 4 * AG / (SI + SG),
			MuG, PipeRoughness);

		//shear terms
		ShearF = fF * RhoL*abs(VLTB)*VLTB / 2.0;
		ShearG = fG * RhoG*abs(VGTB)*VGTB / 2.0;
		ShearI = fG * RhoG*abs(VGTB - VLTB)*(VGTB - VLTB) / 2.0;

		//If interface shear is neglected, it affects the the results massively
		ans = ((ShearF*SF / AF) - (ShearG*SG / AG) - (ShearI*SI*(1 / AF + 1 / AG)))
			/ ((RhoL - RhoG)*g - (RhoL*VF*(VTb - VLLS)*HLLS*dHLTBdhF / HL / HL)
				- (RhoG*VG*(VTb - VGLS)*(1 - HLLS)*dHLTBdhF / (1 - HL) / (1 - HL)));
		
		return(ans);
	}

	double FilmProfile::hFCritical(double hF)
	{
		//Stratified Geometric parameters 
		double SG, SF, SI;
		double AF, AG;

		//Relative velocities to VTb
		double VF, VG;
		//other variables
		double hLTilt = hF / ID;//This is not liquid phase holdup
		double HL; //This is liquid holdup
		double dHLTBdhF;// dHLTB/dhF
		hLTilt = hF / ID;

		//Now we can calculate the liquid phase holdup
		HL = (1 / Pi)*(Pi - acos(2 * hF / ID - 1) + (2 * hF / ID - 1)*
			pow(1 - pow(2 * hF / ID - 1, 2), 0.5));
		StratifiedFlowGeometry Geometry(HL,ID);
		AF = Geometry.ALTilt*ID*ID;
		AG = Geometry.AGTilt*ID*ID;
		SG = Geometry.SGTilt*ID;
		SF = Geometry.SLTilt*ID;
		SI = Geometry.SITilt*ID;

		//Relative velocities (referneced to VTb)
		VF = (VTb - VLLS)*HLLS / HL;
		VG = (VTb - VGLS)*(1 - HLLS) / (1 - HL);

		//Actual velocities of film and gas core in stratified region
		VLTB = VTb - VF;
		VGTB = VTb - VG;

		dHLTBdhF = (4.0 / (Pi*ID))*pow(1 - pow(2 * hF / ID - 1, 2), 0.5);

		return(((RhoL - RhoG)*g - (RhoL*VF*(VTb - VLLS)*HLLS*dHLTBdhF / HL / HL)
			- (RhoG*VG*(VTb - VGLS)*(1 - HLLS)*dHLTBdhF / (1 - HL) / (1 - HL))));
	}

	double FilmProfile::hFCriticalDetermination()
	{
		double HL;
		double dx = 0.00001;
		double x, FPrime = 0;
		double Error = 10;
		HLLS = 1 / (1 + pow(Vm / 8.66, 1.39));
		StratifiedFlowGeometry Geometry(HLLS,ID);
		x= Geometry.hLTilt *ID*0.81;
		while (Error > 0.00001)
		{
			FPrime = (hFCritical(x +dx) - hFCritical(x)) / dx;
			x = x - hFCritical(x) / FPrime;
			Error = abs(hFCritical(x));
		}
		//err = hFCritical(x, hMax, EpsilonZh);
		HL=(1 / Pi)*(Pi - acos(2 * x / ID - 1) + (2 * x / ID - 1)*
			pow(1 - pow(2 * x / ID - 1, 2), 0.5));
		return(x);
	}
	//In this method, the the film profile is calculated and 
	//mass balanced is checked from a given film length
	double FilmProfile::hFDetermination(double Lf, double hMax,double EpsilonZh) {
		hFArray.clear();
		zArray.clear();
		std::ofstream output;
		output.open("Output.txt");
		int zLoop;
		//Average HLTB based on predicted film profile

		//Unit slug length
		double Lu;
		Lu = Lf + Ls;

		//RK45 variables
		double z = 0;
		double h = 0.000001;
		double RKFEps = 0.00001;
		double TotalZ = Lf;
		double StepRatio = 1;
		double hMin = 0.000001;
		double TotalCount, F;
		double XxAtZ, hminCondition;
		double xx;
		double KK[6];
		double ERK4, ERKFVal;
		double HLTBPrev;
		//Error
		double Eps = 10, Eps2 = 10;
		//location dependent film liquid holdup and velociy 
		


		//Other hydrodynamic parameters
		double x;
		double hF;
		double hFLLS; //hF for HLLS 
		double fS;
		double HLTBComparison = 10;

		//liquid mass flow rate
		double WLCal; //Calculated liquid mass flow rate to be checked
		double WL = VSL * Pi*0.25*ID*ID*RhoL; //input liquid mass flow rate

		//other variables
		double gEnglish = 32.2;
		double Sum = 0, Sum2 = 0, Sum3 = 0, SumHfAvg = 0;
		double Flag = 0;
		double Tsm;
		double PreVal = 0;


		k = 0;
		HLLS = 1 / (1 + pow(Vm / 8.66, 1.39));



		//In this while loop, the HLLS from Zhang 2003 is calculated
		//His approach is in implicit form and that is why it required 
		//iterative aprroach to be olved
		while (HLTBComparison > EpsilonZh)
		{
			hFArray.clear();
			zArray.clear();
			z = 0;
			h = 0.0000001;
			RKFEps = 0.00001;
			TotalZ = Lf;
			double StepRatio = 1;
			hMin = 0.000001;
			k = k + 1;
			//liquid height for HLLS is calculated
			StratifiedFlowGeometry Geometry(HLLS,ID);
			hFLLS = Geometry.hLTilt*ID*0.99;
			if (hFdZ(hFLLS) > 0)
			{
				hF = hFCriticalDetermination()*0.99;
			}
			else {
				hF = hFLLS;
			}

			hFArray.push_back(hF);
			zArray.push_back(0);
			//This is boundary condition for ODE of dhF/dz 

			//These hydrodynamic parameters are recalculated using the new HLLS
			x = (VTb - VLLS)*RhoL*0.25*Pi*ID*ID*HLLS;
			VLLS = ((VSL + VSG) - VGLS * (1 - HLLS)) / HLLS;
			RhoS = RhoL * HLLS + RhoG * (1 - HLLS);
			MuS = MuL * HLLS + MuG * (1 - HLLS);

			//In this "for loop", axial locations are iterated
			for (zLoop = 0; zLoop < 4000; zLoop++) {
				if (z >= TotalZ) {

					break;
				}
				XxAtZ = hF;
				F = hFdZ(hF);

				hminCondition = 0;

				//This is RK45 calculation
				for (int RKFLoop = 1; RKFLoop < 1000; RKFLoop++) {
					h = h * StepRatio;
					if (h > hMax) {
						h = hMax;
					}
					if (h < hMin) {
						h = hMin;
						hminCondition = 1;
					}

					if (z + h >= TotalZ) {
						h = TotalZ - z;
					}

					//Cal KK1
					xx = hF;
					F = hFdZ(hF);
					KK[0] = F;

					//Cal KK2
					xx = XxAtZ + h / 4.0*KK[0];
					F = hFdZ(xx);
					KK[1] = F;

					//Cal KK3
					xx = XxAtZ + h * (3.0 / 32.0*KK[0] + 9.0 / 32.0*KK[1]);
					F = hFdZ(xx);
					KK[2] = F;

					//Cal KK4
					xx = XxAtZ + h * (1932.0 / 2197.0*KK[0] - 7200.0 /
                            2197.0*KK[1] + 7296.0 / 2197.0*KK[2]);
					F = hFdZ(xx);
					KK[3] = F;

					//Cal KK5
					xx = XxAtZ + h * (439.0 / 216.0*KK[0] - 8.0*KK[1] + 3680.0 / 
                            513.0*KK[2] - 845.0 / 4104.0*KK[3]);
					F = hFdZ(xx);
					KK[4] = F;

					//Cal KK6
					xx = XxAtZ + h * (-8.0 / 27.0*KK[0] + 2.0*KK[1] - 3544.0 / 2565.0*KK[2]
                            + 1859.0 / 4104.0*KK[3] - 11.0 / 40.0*KK[4]);
					F = hFdZ(xx);
					KK[5] = F;

					ERK4 = abs(h*(KK[0] / 360.0 - 128.0 / 4275.0*KK[2] - 2197.0 / 75240.0*
                 KK[3] + KK[4] / 50.0 + 2.0 / 55.0*KK[5]));

					ERKFVal = ERK4;
					if (ERKFVal != 0) {
						StepRatio = pow((abs(RKFEps*h / 2.0 / ERKFVal)), 0.25);
					}
					else {
						StepRatio = 1.25;
					}

					if (StepRatio < 0.05) {
						StepRatio = 0.05;
					}

					if (StepRatio >= 1.0) {
						break;
					}

					if (RKFLoop > 999) {
						std::cout << "WRONG";
					}
				}

				//At one specific axial location, hF is calculated
				XxAtZ = XxAtZ + h * (16.0 / 135.0*KK[0] + 6656.0 / 12825.0*KK[2] +
					28561.0 / 56430.0*KK[3] - 9.0 / 50.0*KK[4] + 2.0 / 55.0*KK[5]);
				hF = XxAtZ;

	
				z = z + h;
				hFArray.push_back(hF);
				zArray.push_back(z);

				Sum = 0;
				Flag = 0;
				Sum2 = 0;
				Sum3 = 0;
				SumHfAvg = 0;
				for (int jj = 1; jj < zLoop + 2; jj++)
				{
					if (zLoop == 10)
					{
						double a;
						a = 1;
					}
					//HLTB is calculated for each predicted hF
					if (jj > 0) {
						HLTBPrev = (1 / Pi)*(Pi - acos(2 * (hFArray[jj - 1]) / ID - 1) + (2 * (hFArray[jj - 1]) / ID - 1)*pow(1 - pow(2 * (hFArray[jj - 1]) / ID - 1, 2), 0.5));;
					}
					else {
						HLTBPrev = 0;
					}
					HLTB = (1 / Pi)*(Pi - acos(2 * (hFArray[jj]) / ID - 1) + (2 * (hFArray[jj]) / ID - 1)*pow(1 - pow(2 * (hFArray[jj]) / ID - 1, 2), 0.5));
					HLTB = (HLTB + HLTBPrev) / 2.0;
					hFAvg = (hFArray[jj] + hFArray[jj - 1]) / 2.0;


					//Integral in mass balance relation
					Sum = Sum + RhoL * Pi * ID*ID*0.25*HLTB*(zArray[jj] - zArray[jj - 1]);

					//To calculate the average HLTB
					Sum2 = Sum2 + HLTB * (zArray[jj] - zArray[jj - 1]);

					SumHfAvg = SumHfAvg + hFAvg * (zArray[jj] - zArray[jj - 1]);

					//Sum3 = Sum3 + (VTb / Lu)*(1 - HLTB)*(zArray[jj] - zArray[jj - 1]);
				}

				//Trapezoidal method for calculation of average HLTB
				double SumNew = 0;
				for (int jj = 1; jj < zLoop +2; jj++)
				{
					double HLTBNew = (1 / Pi)*(Pi - acos(2 * (hFArray[jj]) / ID - 1) +
                           (2 * (hFArray[jj]) / ID - 1)*pow(1 - 
                           pow(2 * (hFArray[jj]) / ID - 1, 2), 0.5));
					double HLTBPast = (1 / Pi)*(Pi - acos(2 * (hFArray[jj-1]) / 
                             ID - 1) + (2 * (hFArray[jj-1]) / ID - 1)*
                            pow(1 - pow(2 * (hFArray[jj-1]) / ID - 1, 2), 0.5));
					HLTB = (HLTBPast + HLTBNew) / 2.0;
					double m = (HLTBNew - HLTBPast) / (zArray[jj] - zArray[jj-1]);
					double x1, x2, x3;
					double y1, y2, y3;
					x1 = zArray[jj-1] + (zArray[jj] - zArray[jj - 1]) / 3.0;
					y1 = m * x1 - m * zArray[jj-1] + HLTBPast;

					x2 = x1 + (zArray[jj] - zArray[jj - 1]) / 3.0;
					y2 = m * x2 - m * zArray[jj-1] + HLTBPast;

					Sum3 = Sum3 + (VTb / Lu)*(1 - HLTB)*(zArray[jj] - zArray[jj - 1]);
					SumNew = SumNew + ((zArray[jj] - zArray[jj - 1]) / 6.0)*
                            (HLTBPast + 2 * y1 + 2 * y2 + HLTBNew);
				}
						//Trapezoidal method HLTBAvg
				HLTBAvg=SumNew / Lf;
				StratifiedFlowGeometry Geom(HLTBAvg,ID);
				hFAvg = Geom.hLTilt*ID;
				//Calculated mass flow rate
				WLCal = (RhoL * (Ls)*Pi * ID*ID*0.25*HLLS + RhoL * Pi * 
                 ID*ID*0.25*HLTBAvg*Lf)*(VTb / Lu) - x;

				//Error calculation of calculated mass flow rate and input mass flow rate
				Eps = WLCal - WL;
				Eps2 = VSL - VLLS * HLLS - VTb * (1 - HLLS)*Lf / Lu + Sum3;

				//Average liquid film holdup 

				//HLTBAvg = Sum2 / Lf;
				
				//hFAvg = SumHfAvg / Lf;
				if (Con == 1) {
					std::cout << std::endl;
				}
			}
			VLTB = VTb - (VTb - VLLS)*HLLS / HLTBAvg;
			VGTB = VTb - (VTb - VGLS)*(1 - HLLS) / (1 - HLTBAvg);


			//HLLS calculation from Zhang 2003
			fS = FannningFrictionFactor(RhoS, Vm, ID, MuL, PipeRoughness);
			Tsm = (2.0 / 2.5)*(0.5*fS*RhoS*Vm*Vm + 0.25*ID*RhoL*HLTBAvg*
            (VTb - VLTB)*(Vm - VLTB) / Ls);
			HLLS = 1 / (1 + (Tsm / (3.16*pow((RhoL - RhoG)*g*IT, 0.5))));

			//Comparing newlly calculated HLTBAvg from the new HLLS with HLTBAvg
			//from previous iteration and HLLS 
			HLTBComparison = abs(PreVal - HLTBAvg);
			PreVal = HLTBAvg;
			z = 0;
		}


		ZSize = zLoop;

		return(Eps);
	}
	//In this method, the right Lf will be chosen so the mass is conserved
   //In other words, by this method, correct film length is calculated
   //I used bisection method for soliving for the root
   //the two values of negative and positive values shoudl be checked for each case
	double FilmProfile::LfDetermination(double InitialGuess,double dx, double Epsilon,
                                      double hMax, double EpsilonZh) {
		double Lf = 16, Eps = 10;;
		double F = 10;
		double err1 = 10;
		double check;
		double x,FPrime = 0;

		

		x = InitialGuess;
		err=hFDetermination(x, hMax, EpsilonZh);
		err = 10;
		while (err > Epsilon)
		{
			FPrime = (hFDetermination(x + dx,hMax, EpsilonZh) - 
                hFDetermination(x, hMax, EpsilonZh)) / dx;
			x = x - hFDetermination(x, hMax, EpsilonZh) / FPrime;
			err = abs(hFDetermination(x, hMax, EpsilonZh));
		}
		err = hFDetermination(x, hMax, EpsilonZh);
		return(x);
	}

	void FilmProfile::HyrodynamicProperties(double InitialGuess, double dx,
                                          double Epsilon, double hMax, 
                                          double EpsilonZh) {
		XSave.clear();
		HLTBSave.clear();
		std::ofstream output;
		output.open("Output.txt");
		LfFinal = LfDetermination(InitialGuess,dx,Epsilon, hMax, EpsilonZh);
		hFDetermination(LfFinal, hMax,EpsilonZh);

			StratifiedFlowGeometry Geo(HLTBAvg,ID);
			output << "Calculated Error:	" << err << std::endl;
			output << "Pressure:	" << P / 6894.76 << "	Psia" << std::endl;
			output << "Liquid superficial velocity (VSL):	" << VSL << "	m/s" << std::endl;
			output << "Gas superficial velocity (VSG):	" << VSG << "	m/s" << std::endl;
			output << "Liquid density:	" << RhoL << "	kg/m3" << std::endl;
			output << "Gas density:	" << RhoG << "	kg/m3" << std::endl;
			output << "Liquid viscosity:	" << MuL << "	Pa.s" << std::endl;
			output << "Gas viscosity:	" << MuG << "	Pa.s" << std::endl;
			output << "Liquid film holdup:	" << HLTBAvg << std::endl;
			output << "Liquid slug holdup:	" << HLLS << std::endl;
			output << "SL:	" << (Geo.SLTilt)*ID << "	m" << std::endl; //0.16485 0.06954
			output << "SI:	" << Geo.SITilt*ID << "	m" << std::endl;
			output << "Film velocity (Lf):	" << VLTB << "	m/s" << std::endl;
			output << "Gas core velocity (Vg):	" << VTb - (VTb - VGLS)*(1 - HLLS) /
                 (1 - HLTBAvg) << "	m/s" << std::endl;
			output << "Unit slug length (Lu):	" << LfFinal+Ls << "	m" << std::endl;
			output << "Film length (Lf):	" << LfFinal << "	m" << std::endl;
			output << "Slug length (Ls):	" << Ls << "	m" << std::endl;
			output << std::endl;
			output << "Z	" << "HLTB:	" << std::endl;

			for (int j = 0; j < ZSize + 1; j++) {
				HLTB = (1 / Pi)*(Pi - acos(2 * (hFArray[j]) / ID - 1) + (2 * (hFArray[j])
               / ID - 1)*pow(1 - pow(2 * (hFArray[j]) / ID - 1, 2), 0.5));
				output << zArray[j] << "	" << HLTB << std::endl;
					XSave.push_back(zArray[j]);
					HLTBSave.push_back(HLTB);
			}

	}


}//namespace MultiphaseTemperature




//Some formulas for Frequency::
//Freq = 0.0226*(pow((VSL / gEnglish / ID), (1.2)))*(pow((212.6 / (VSL + VSG) / 3.28084 + (VSL + VSG)*3.28084), (1.2)))*(0.836);
//double Freq = 0.19;
//Freq = exp(0.8 + 1.53*log(VSL / 0.36) + 0.27*((VSG / (1 - 0.36)) - VSL / 0.36) / Vm - 34.1*ID);
//Freq = 0.0226*pow((VSL / g / ID)*(19.75 / Vm + Vm), 1.2);
//Freq = 0.8428*pow((VSL / g / (ID * 1000))*(19.75 / Vm + Vm), 0.25);
//Freq = 0.0364*(VSL / Vm)*pow(2.02 / ID + Vm * Vm / g / ID, 1.06); //Not good
//Freq = 0.47*pow(VSL, 0.75) / (pow(ID, 1.2)*pow(0.3048 * 53, 0.55)); //Not good
//Freq = exp(0.8 + 1.53*log(VSL / HLTB) + 0.27*(VTb / Vm) - 34.1*ID); //Not too bad
//Freq = 2.623*(VSL / ID)*(1 / pow(pow(ID, 3.0 / 2.0)*pow(RhoL*(RhoL - RhoG)*g, 0.5) / MuL, 0.612)); //hugely underestimated
//Freq = 0.61 - (RhoG - VSG / (1 - 0.36)) / (RhoL*(ID - 0.36*ID)); //Not good
//Freq = 0.088*(VSL + 1.5)*(VSL + 1.5) / g / ID; //Overestimated
//Freq = 1.2*VSL / (32 * ID);
//Lu = VTb / Freq;

