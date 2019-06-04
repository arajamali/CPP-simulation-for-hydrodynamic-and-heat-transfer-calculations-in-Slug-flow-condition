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

namespace MultiphaseTemperature {

	double s1 = 0; bool s2;

	void SlugFLowHeatTransferCal::SlugFlowHeatModel(std::vector<double>& FilmInputs,
		std::vector<std::vector<double>>& WaxWidthVector, std::vector<std::vector<double>>& WaxSolidFraction ,
		std::vector<double>& HydroInputs, std::vector<std::vector<std::vector<double>>>& TwInitial,
		std::vector<double>& Ts, std::vector<double>& Tf, std::vector<double>& Tg, std::vector<double>& FilmAssign,
		std::vector<bool>& SlugAssign, std::vector<double>& q,double TInlet) {


		SetVariablesFromIputFiles(FilmInputs,WaxWidthVector, WaxSolidFraction,
			TwInitial,Ts,Tf,Tg, FilmAssign, SlugAssign);


		//The requried vector sizes are assigned  
		VectorResize();
		//In this simplified case, Inlet temperature is assumed to be constant and the same as 
	
		AvgTb = 0;
		//Pipe specs
		const user::PipeGeometry PipeProps;
		PipeSpecs(PipeProps);

		//This method uses the pipe's ID from the above method (PipeSpecs(PipeProps);)
		IDEf = IDEff(WaxWidthVector);

			Update(AvgTempForHydroUpdate, HLLS);
	


				for (int xx = 0; xx < Nx; xx++)
				{
					/*Effective surface area of pipe is different axially so, and
					   it is calculated by averaging the 4 peaces surface areas of the pipe*/
					ApEff[xx] = 0.25*M_PI*(pow((ri - WaxThickness[xx][1]), 2)
						+ pow((ri - WaxThickness[xx][2]), 2) + pow((ri - WaxThickness[xx][3]), 2)
						+ pow((ri - WaxThickness[xx][4]), 2));


					//ApEff[xx] = M_PI * ri*ri;


					if (SlugFlow[xx] == false)
					{
						//flow geometries are calculated and reported by the following method
						HLTBAvg = HLTbLocation[xx];
						const StratifiedFlowGeometry FlowGeom(HLTBAvg, IDEf);
						UpdateFlowGeomParams(FlowGeom, HLTBAvg);
					}

					for (int WL = 1; WL < WallLay + 1; WL++) //starts from 1 to 4
					{

						for (int RS = 1; RS < RingSec; RS++)  //starts from 1 to 3. Due to symmetry, Ring 2 is assigned to Ring 4
						{

							//Thermal conductivity coefficient is assingned for wax and pipe layers
							if (WL == 1) {
								KThermal = WaxDepoK(Tw[xx][WL][RS], WaxSolidFrac[xx][RS]);
								KDepoSave[xx][RS] = KThermal;
							}
							else {
								KThermal = KPipe;
							}

							if (WL == 1) //Inner layer contacting the fluid
							{
								//This method assignes the AG and AL arrays
								SurfaceAreaAGAL(thetahF, WaxThickness[xx][3], WaxThickness[xx][2], WaxThickness[xx][1]);
							}

							//Es calcualations
							if (WL == WallLay)
							{
								Es[xx][WL][RS] = q[xx] * Aout(WL)*dt;
								check = Aout(WL);

							}
							else
							{
								Es[xx][WL][RS] = Aout(WL) * KThermal*(Tw[xx][WL + 1][RS] - Tw[xx][WL][RS])*dt*3.0 / w;//double T[NTime + 1][Nx + 1][WLay + 1][RingSec + 1];

							}

							//Er calcualations
							//For the most inner layer, the caculations are calculated based on theta
							if (WL == 1) //Heat is convectively transfered from liquid to the wall
							{

								if (SlugFlow[xx] == true)
								{
									//The following if, else is to take care of the x=0 and x=N-1

									if (xx == Nx - 1)
									{
										//The boundary condition of dT/dx  is assumed for the end part
										Update(TS[xx], HLLS);
										Er[xx][WL][RS] = -Ain(WL, WaxThickness[xx][RS])*HeatCoe(KS, IDEf, RhoS, Vm, MuS, CpS) *dt*(Tw[xx][WL][RS] - TS[xx]);
									}
									else //Other xx other than the last one
									{
										Update(TS[xx], HLLS);

										Er[xx][WL][RS] = -Ain(WL, WaxThickness[xx][RS])*HeatCoe(KS, IDEf, RhoS, Vm, MuS, CpS) *dt*((Tw[xx + 1][WL][RS] + Tw[xx][WL][RS]) / 2.0 - TS[xx]);

									}
									EWsLocal[RS] = -Er[xx][WL][RS];
									if (RS == 2)
										EWsLocal[4] = EWsLocal[2];
								}
								else //Film section
								{
									//the following if and else statements are to care of x=N-1
									if (xx == Nx - 1)
									{
										AvgTemp = Tw[xx][WL][RS];//This variable is used in Er calculations
									}
									else
									{
										AvgTemp = (Tw[xx + 1][WL][RS] + Tw[xx][WL][RS]) / 2.0;//This variable is used in Er calculations
									}

									//EwfLocal and EwgLocal 
									Update(TL[xx], HLLS); EwfLocal[RS] = AL[RS] * HeatCoe(KL, dL, RhoL, VLTB,
										MuL, CpL)*dt*(AvgTemp - TL[xx]);
									Update(TG[xx], HLLS); EwgLocal[RS] = AG[RS] * HeatCoe(KG, dG, RhoG, VGTB,
										MuG, CpG)*dt*(AvgTemp - TG[xx]);
									//Er calculation
									Er[xx][WL][RS] = -(EwfLocal[RS] + EwgLocal[RS]);

									//Assigning RS=2 to RS=4 (sidewall section are the same)
									if (RS == 2) { EwfLocal[4] = EwfLocal[2]; EwgLocal[4] = EwgLocal[2]; }

								}
							}
							else //Wall layers other than the first one
							{
								Er[xx][WL][RS] = -Ain(WL, WaxThickness[xx][RS]) * KThermal*(Tw[xx][WL][RS] - Tw[xx][WL - 1][RS])*dt*3.0 / w;
							}
							//the following if and else statements are to care of x=N-1
							if (xx == Nx - 1)
							{
								Ex[xx][WL][RS] = Ax(WL, WaxThickness[xx][RS]) * KThermal*dt*(1 / dX)*(Tw[xx][WL][RS] -
									2 * Tw[xx][WL][RS] + Tw[xx - 1][WL][RS]);
							}
							else if (xx == 0)
							{
								Ex[xx][WL][RS] = Ax(WL, WaxThickness[xx][RS]) * KThermal*dt*(1 / dX)*(Tw[xx + 1][WL][RS] -
									2 * Tw[xx][WL][RS] + TInlet);
							}
							else // xx other than xx==0 and xx=Nx-1
							{
								Ex[xx][WL][RS] = Ax(WL, WaxThickness[xx][RS]) * KThermal*dt*(1 / dX)*(Tw[xx + 1][WL][RS] -
									2 * Tw[xx][WL][RS] + Tw[xx - 1][WL][RS]);
							}

							if (RS == 1)
							{
								ETheta[xx][WL][RS] = ATheta(WL, WaxThickness[xx][RS]) * KThermal*dt*(1 / dxTheta(WL, WaxThickness[xx][RS]))*
									(Tw[xx][WL][RS + 1] - 2 * Tw[xx][WL][RS]
										+ Tw[xx][WL][4]);
							}
							else { //Ring section other than the first one 

								ETheta[xx][WL][RS] = ATheta(WL, WaxThickness[xx][RS]) * KThermal*dt*(1 / dxTheta(WL, WaxThickness[xx][RS]))*
									(Tw[xx][WL][RS + 1] - 2 * Tw[xx][WL][RS]
										+ Tw[xx][WL][RS - 1]);
							}

							//balance equations for ring elements
							if (WL == 1)
							{

								//CPL was used for the heat capacity of the wax deposition, PLEASE CHECK!!
								Tw[xx][WL][RS] = Tw[xx][WL][RS] + (Es[xx][WL][RS] +
									Er[xx][WL][RS] + Ex[xx][WL][RS] + ETheta[xx][WL][RS])
									/ (Ax(WL, WaxThickness[xx][RS])*dX*RhoDepo*CpL);
							}
							else
							{

								Tw[xx][WL][RS] = Tw[xx][WL][RS] + (Es[xx][WL][RS] +
									Er[xx][WL][RS] + Ex[xx][WL][RS] + ETheta[xx][WL][RS])
									/ (Ax(WL, WaxThickness[xx][RS])*dX*RhoPipe*CPPipe);
							}

							//std::cout << "Es: " << Es[xx][WL][RS] << " Er: " << Er[xx][WL][RS]
								//<< " Ex: " << Ex[xx][WL][RS] << " ETheta: " << ETheta[xx][WL][RS] << std::endl;
							//std::cout << std::setprecision(10)<< Tw[xx][WL][RS] << " ";


							//std::cout <<std::setprecision(10)<< Tw[tt + 1][xx][WL][RS] << " ";

							if (RS == 2)
							{
								//In this section, the energy of section 4 will be assigned by the energies of section 2
								Es[xx][WL][4] = Es[xx][WL][RS];
								Er[xx][WL][4] = Er[xx][WL][RS];
								Ex[xx][WL][4] = Ex[xx][WL][RS];
								ETheta[xx][WL][4] = ETheta[xx][WL][RS];
								Tw[xx][WL][4] = Tw[xx][WL][RS];
								//std::cout << std::setprecision(10) << Tw[tt + 1][xx][WL][4] << " ";
							}
						}//End of Ring section of one ring

					}//End of Wall layer

					if (SlugFlow[xx] == false)
					{
						sumWF = 0;
						sumWG = 0;
						for (int itr = 1; itr < RingSec + 1; itr++) //In here the forth ring section should also be included
						{
							sumWF = EwfLocal[itr] + sumWF;
							sumWG = EwgLocal[itr] + sumWG;
						}
						EWf[xx] = sumWF;
						EWG[xx] = sumWG;
					}
					else
					{
						sumWS = 0;
						for (int itr = 1; itr < RingSec + 1; itr++)//In here the forth ring section should also be included
						{
							sumWS = EWsLocal[itr] + sumWS;
						}
						EWS[xx] = sumWS;
					}


					//the following if and else statements are to take care of x=N-1 Based on dT/dx
					if (xx == Nx - 1)
					{
						if (SlugFlow[xx] == false)
						{
							ECF[xx] = 0;
							ECG[xx] = 0;
						}
						else // slug
						{
							ECS[xx] = 0;
						}
					}
					else
					{
						if (SlugFlow[xx] == false)
						{
							Update(TL[xx], HLLS);

							if (HLTbLocation[xx + 1] > 0.000001) {
								HLTB2 = HLTbLocation[xx + 1];
							}
							else {
								HLTB2 = HLTBAvg;
							}


							Update(TL[xx + 1], HLLS);
							//ECF[xx] = ECF2- ECF1;
							Update((TL[xx + 1] + TL[xx]) / 2.0, HLLS);
							ECF[xx] = (VTb - VLTB)*ApEff[xx] * HLTBAvg*RhoL*CpL*dt*(TL[xx + 1] - TL[xx]);
							//ECF[xx] = ((VTb - VLTB)*Ap*HLTBAvg*RhoL*CpL*dt)*(TL[xx + 1] - TL[xx]);//check the Ap. I chose the Ap as the inner surface of the pipe. It could be the sirface of the pipe itself

							Update(TG[xx], HLLS);
							ECG1 = Cs * ApEff[xx] * Vm*(1 - HLLS)*dt*CpG*RhoG*TG[xx];

							Update(TG[xx + 1], HLLS);
							ECG2 = Cs * ApEff[xx + 1] * Vm*(1 - HLLS)*dt*CpG*RhoG*TG[xx + 1];
							//ECG[xx] = ECG2 - ECG1;
							Update((TG[xx + 1] + TG[xx]) / 2.0, HLLS);
							ECG[xx] = Cs * ApEff[xx] * Vm*(1 - HLLS)*dt*CpG*RhoG*(TG[xx + 1] - TG[xx]);
							//ECG[xx] = Cs * Ap*Vm*(1 - HLLS)*dt*CpG*RhoG*(TG[xx + 1] - TG[xx]);
						}
						else // Slug
						{
							Update(TS[xx], HLLS);

							//ECS[xx] = ECS2 - ECS1;
							Update((TS[xx + 1] + TS[xx]) / 2.0, HLLS);
							ECS[xx] = Cs * ApEff[xx] * Vm*RhoS*CpS*dt*(TS[xx + 1] - TS[xx]);
						}
					}

					TwAv[xx] = (Tw[xx][1][1] + Tw[xx][1][2] + Tw[xx][1][3] + Tw[xx][1][4]) / 4.0;

					//Now the balance equations can be written for TL,TG and TS
					if (SlugFlow[xx] == false)
					{
						//needs fix
						UInt = 1 / (1 / (HeatCoe(KG, dG, RhoG, abs(VGTB - VLTB), MuG, CpG)) + 1 / (HeatCoe(KL, dL, RhoL, abs(VGTB - VLTB), MuL, CpL)));
						//Future TG
						TG[xx] = TG[xx] + (EWG[xx] + ECG[xx] + SI * dX* UInt*dt*(TL[xx] - TG[xx])) / ((1 - HLTBAvg)*ApEff[xx] * dX*CpG*RhoG);
						//Future TL
						TL[xx] = TL[xx] + (EWf[xx] + ECF[xx] + SI * dX * UInt*dt*(TG[xx] - TL[xx])) / (HLTBAvg*ApEff[xx] * dX*CpL*RhoL);
						TLActual[xx] = TL[xx];
						TGActual[xx] = TG[xx];
						//Future Tb
						Tb[xx] = (TL[xx] * VLTB*ALP + TG[xx] * VGTB*AGP) / (VLTB*ALP + VGTB * AGP);
						AvgTb = AvgTb + Tb[xx];
						Kb[xx] = (KLSave * VLTB*ALP + KG * KGSave*AGP) / (VLTB*ALP + VGTB * AGP);
						//Kb[xx] = KLSave *HLTBAvg+(1- HLTBAvg)*KGSave;
						TS[xx] = Tb[xx];

						//For the next timestep
						QFSum = QFSum + (EWf[xx] + EWG[xx])*dX / (dX*dt*M_PI*IDEf);
						TfSum = TfSum + (-Tb[xx] + TwAv[xx])*dX;
						hHotSave[xx] = QFSum / TfSum;
						
					}

					else //Slug
					{

						Update(TS[xx], HLLS);
						TG[xx] = TG[xx] + (1 - HLLS)*(ECS[xx] + EWS[xx]) / (ApEff[xx] * dX*RhoS*CpS);
						TL[xx] = TL[xx] + HLLS * (ECS[xx] + EWS[xx]) / (ApEff[xx] * dX*RhoS*CpS);
						TS[xx] = TS[xx] + (ECS[xx] + EWS[xx]) / (ApEff[xx] * dX*RhoS*CpS);
						TLActual[xx] = TS[xx];
						TGActual[xx] = TS[xx];
						Tb[xx] = TS[xx];
						AvgTb = AvgTb + Tb[xx];
						Kb[xx] = KS;

						//For convective heat tranfser coefficient calcualtion
						QSSum = QSSum + (EWS[xx])*dX / (dX*dt*M_PI*IDEf);
						TSSum = TSSum + (-Tb[xx] + TwAv[xx])*dX;
						hHotSave[xx] = QSSum / TSSum;


					}
					Update(Tw[xx][1][1], HLLS);
					double help = MuL;
					Update(Tb[xx], HLLS);
					if (xx != 0)
					{
						SumhCheck = SumhCheck + (1 / IDEf)*KL*SinglePhaseNu(RhoL, VSL, IDEf, MuL, CpL, KL, xx*dX, PipeRoughness, VSL / (VSL + VSG));
					}

				}//End of axial loop

				AvgTb = AvgTb / double(Nx);

				NuNumCal = (QSSum + QFSum) / (TSSum + TfSum);
				SumhCheck = SumhCheck / double(Nx);
		
	}//End of constructor


	void SlugFLowHeatTransferCal::SetVariablesFromIputFiles(std::vector<double>& FilmInputs,
		std::vector<std::vector<double>>& WaxWidthVector, std::vector<std::vector<double>>& WaxSolidFraction,
		std::vector<std::vector<std::vector<double>>>& TwInitial,
		std::vector<double>& Ts, std::vector<double>& Tf, std::vector<double>& Tg, std::vector<double>& FilmAssign,
		std::vector<bool>& SlugAssign)
	{


		//Reading input vectors and assiging them into double variables
		VSL = FilmInputs[0];
		VSG = FilmInputs[1];
		//T = SlugFLowHeatTransferCalInputs[2];
		AvgTempForHydroUpdate= FilmInputs[3];
		/*	Please note that this temperature is used as a defualt
		temperature value for themorpysical property calcualtion*/
		P = FilmInputs[4];
		PipeRoughness = FilmInputs[5];
		IT = FilmInputs[6];//InterfactialTension

		WaxThickness = WaxWidthVector;
		WaxSolidFrac = WaxSolidFraction;
		Tw = TwInitial;
		TL = Tf;
		TG = Tg;
		TS = Ts;
		HLTbLocation = FilmAssign;
		SlugFlow = SlugAssign;

	}

	//In this method, heat transfer coefficient is calculated 
	double SlugFLowHeatTransferCal::HeatCoe(double K, double d, double Rho,
		double V, double Mu, double Cp)
	{
		double ans;
		ans = (0.023*K / d) * pow((Rho*V * d / Mu), 0.8) * pow((Cp*Mu / K), 1.0 / 3.0);
		return(ans);
	}

	double SlugFLowHeatTransferCal::DistCoe(double Rho, double V, double d, double Mu)
	{
		double Re, co, c;
		Re = Rho * V*d / Mu;
		if (Re <= 2100) {
			co = 2;
			c = co - 1;
		}
		if (Re > 2100) {
			co = 1.2;
			c = co - 1;
		}
		return(c);
	}

	void SlugFLowHeatTransferCal::Update(double T,double HLLS)
	{
		user::FluidProperties Fluid(T-273.15, P);

		RhoL = Fluid.Liquid.Rho(Fluid);
		RhoG = Fluid.Gas.Rho(Fluid);
		MuL = Fluid.Liquid.Mu(Fluid);
		MuG = Fluid.Gas.Mu(Fluid);
		KL = Fluid.Liquid.K(Fluid);
		KG = Fluid.Gas.K(Fluid);
		CpL = Fluid.Liquid.Cp(Fluid);
		CpG = Fluid.Gas.Cp(Fluid);

		//Slug properties
		KS = HLLS * KL + (1 - HLLS)*KG;
		RhoS = HLLS * RhoL + (1 - HLLS)*RhoG;
		MuS = HLLS * MuL + (1 - HLLS)*MuG;
		CpS = HLLS * CpL + (1 - HLLS)*CpG;
	}


	//Please note that the temperature in FilmInputs needs to be in in Centigrade 
	void SlugFLowHeatTransferCal::HydroUpdate(std::vector<double>&FimInputs, std::vector<double>& HydroInputs)
	{	
		double  InitialGuess, dxx, Epsilon, hMax, EpsilonZh;

		std::vector<double> xAxisOrg, HLTBOrg;
		double Temp;
		//HtdroInputs for final hydrodynamic properties of slug flow
		InitialGuess = HydroInputs[0];
		dxx = HydroInputs[1];
		Epsilon = HydroInputs[2];
		hMax = HydroInputs[3];
		EpsilonZh = HydroInputs[4];

		//Pipe specs
		const user::PipeGeometry PipeProps;
		PipeSpecs(PipeProps);

		FilmProfile Film(FimInputs);
		PressureCal Press(FimInputs, HydroInputs);
		//These two variables are for VSL and VSG for X parameter
		PressFilm = Press.PressSFilm;
		PressGas = Press.PressSGas;
		//Calculating the film profile
		Film.HyrodynamicProperties(InitialGuess, dxx, Epsilon, hMax, EpsilonZh);
		VTb = Film.VTb;
		HLLS = Film.HLLS;//NP
		VLLS = Film.VLLS;//NP
		VGLS = Film.VGLS;//NP
		VGTB = Film.VGTB;
		VLTB = Film.VLTB;
		Ls = Film.Ls;//NP
		Vm = Film.Vm;//NP
		Lf = Film.LfFinal;//NP
		Lu = Lf + Ls;
		//HLTBAvg = Film.HLTBAvg;
		ro = OD / 2.0;
		ri = IDPipe / 2.0;
		w = ro - ri;
		dX = Lu / double(Nx);
		dt = dX / VTb;
		NxS = round(Ls / dX);
		NxF = round(Lf / dX);

		//Film =sdf
		for (int i = 0; i < Film.ZSize + 1; i++)
		{
			xAxisOrg.push_back(Film.XSave[i]);
			HLTBOrg.push_back(Film.HLTBSave[i]);
		}
		
		for (int i = 0; i < NxF; i++) //check for Nx+1
		{
			xAxisNew.push_back(i*dX);
		}
		//In this section, the adjusted film profile is calculated
		double x1, x2,y1,y2;
		double m;
		for (int i = 0; i < NxF; i++) 
		{
			for (int j = 0; j < Film.ZSize + 1; j++)
			{
				if (xAxisNew[i] < Film.XSave[j])
				{
					x1 = Film.XSave[j-1];
					x2 = Film.XSave[j];
					y1 = Film.HLTBSave[j-1];
					y2 = Film.HLTBSave[j];
					m = (y2 - y1) / (x2 - x1);
					HLTBOrgNew.push_back(y1 + m * (xAxisNew[i] - x1));
					//std::cout << xAxisNew[i] << " " << HLTBOrgNew[i] << std::endl;
					break;
				}

			}
		}
		double n,X_middle, Y_middle;
		for (int i = 0; i < NxF-1; i++)
		{
			HLTBOrgNew[i] = (HLTBOrgNew[i + 1] + HLTBOrgNew[i])*0.5;
			xAxisNew[i] = (xAxisNew[i] + xAxisNew[i + 1])*0.5;
		}

	}

	double SlugFLowHeatTransferCal::SinglePhaseNu(double Rho, double v, double D, 
		double Mu, double Cp, double K,double dx,double E,double H)
	{
		double Nu1, f, Ref, Prf;
		double RHS, Alpha, Phi2,m;
		double X = pow((PressFilm / PressGas), 0.5);
		double RHS2;
		RHS2 = pow(H, -0.194)*(1 + 0.687*pow(X, -0.7));
		m = 2;
		Ref = Rho * v*D / Mu;
		Prf = Mu * Cp / K;
		f = 1.325 / pow(log(E / (3.7*IDEf) + 5.74 / (pow(Ref, 0.9))), 2);
		Nu1 = (((f / 8.0)*(Ref - 1000)*Prf) / (1 + 12.7*pow(f / 8.0, 0.5)*(pow(Prf, 2.0 / 3.0) - 1)))*(1 + pow(D / dx, 2.0 / 3.0));
		Phi2 = 1 / (pow(H, m));
		RHS = pow(H, 1.28)*Phi2;
		return(Nu1*RHS2);
	}
		
	void SlugFLowHeatTransferCal::SurfaceAreaAGAL(double Theta, double Top,double Side, double Bot)
	{

		if (thetahF <= M_PI * 0.5)
		{
			AL[1] = (ri- Bot) * thetahF*dX;
			AG[1] = (ri- Bot) * (0.5*M_PI - thetahF)*dX;

			AL[2] = 0;
			AG[2] = (ri- Side) * 0.5*M_PI*dX;

			AL[3] = 0;
			AG[3] = (ri- Top) * 0.5*M_PI*dX;

			AL[4] = 0;
			AG[4] = (ri- Side) * 0.5*M_PI*dX;
		}//thetahF <= M_PI * 0.5

		if (thetahF > M_PI * 0.5 && thetahF <= 3.0*M_PI * 0.5)
		{
			AL[1] = 0.5*M_PI*(ri- Bot)*dX;
			AG[1] = 0;

			AL[2] = 0.5*(ri- Side)*(thetahF - 0.5*M_PI)*dX;
			AG[2] = 0.5*(ri- Side)*(3.0*M_PI / 2.0 - thetahF)*dX;

			AL[3] = 0;
			AG[3] = (ri- Top) * 0.5*M_PI*dX;

			AL[4] = 0.5*(ri- Side)*(thetahF - 0.5*M_PI)*dX;
			AG[4] = 0.5*(ri- Side)*(3.0*M_PI / 2.0 - thetahF)*dX;
		}//(thetahF > M_PI * 0.5 && thetahF <= 3.0*M_PI * 0.5)

		if (thetahF > 3.0*M_PI * 0.5)
		{
			AL[1] = 0.5*M_PI*(ri- Bot)*dX;
			AG[1] = 0;

			AL[2] = 0.5*M_PI*(ri- Side)*dX;
			AG[2] = 0;

			AL[3] = (ri- Top) * (thetahF - 3.0*M_PI / 2.0)*dX;
			AG[3] = (ri- Top) * (2 * M_PI - thetahF)*dX;

			AL[4] = 0.5*M_PI*(ri- Side)*dX;
			AG[4] = 0;
		}//(thetahF > M_PI * 0.5 && thetahF <= 3.0*M_PI * 0.5)
	}



	void SlugFLowHeatTransferCal::UpdateFlowGeomParams(const StratifiedFlowGeometry& FlowGeom,double HL)
	{
		thetahF = FlowGeom.Theta;
		AGP = FlowGeom.AGTilt*IDEf*IDEf;
		ALP = FlowGeom.ALTilt*IDEf*IDEf;
		SI = FlowGeom.SITilt*IDEf;
		dL = FlowGeom.dL;
		dG = FlowGeom.dG;
		VLTB = VTb - (VTb - VLLS)*HLLS / HL;
		VGTB = VTb - (VTb - VGLS)*(1 - HLLS) / (1 - HL);
		VG = VGTB;
		VL = VLTB;

		Cg = DistCoe(RhoG, VGTB, dG, MuG);
		Cs = DistCoe(RhoS, Vm, IDEf, MuS);
	}

	void SlugFLowHeatTransferCal::SlugFlowBooleanInitial()
	{
		for (int xx = 0; xx < Nx; xx++)
		{
			SlugFlow[xx] = false;
			HLTbLocation[xx] = false;
			if (xx <= (NxF - 1))
			{
				SlugFlow[xx] = false;
				HLTbLocation[xx] = HLTBOrgNew[xx];

			}
			if (xx >(NxF - 1))
			{
				SlugFlow[xx] = true;
			}
			//s1 = HLTbLocation[xx]; s2 = SlugFlow[xx];
		}
	}

	void SlugFLowHeatTransferCal::SlugFlowBooleanLater(std::vector<bool> SlugFlowPre, std::vector<double> HLTbLocationPre)
	{

		/*
		An example for SlugFlow[xx] boolean pre-adjustment
		-----------------------------------SlugFlow[xx]--------------------------------------
		0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1	1
		1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	1
		1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
		1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		0	0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		0	0	0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		0	0	0	0	0	0	1	1	1	0	0	0	0	0	0	0	0	0	0	0	0	0
		*/
		
		std::vector<double> HLTbLocation2;
		std::vector<bool> SlugFlow2;

		SlugFlow2.resize(Nx+1);
		HLTbLocation2.resize(Nx + 1);


		for (int xx = 0; xx < Nx; xx++)
		{
				if (xx != 0)
				{
					SlugFlow2[xx] = SlugFlowPre[xx - 1];
					HLTbLocation2[xx] = HLTbLocationPre[xx - 1];

					//std::cout << HLTbLocation2[xx] << " ";

				}
				else
				{
					SlugFlow2[xx] = SlugFlowPre[Nx - 1];
					HLTbLocation2[xx] = HLTbLocationPre[Nx - 1];
					//std::cout << HLTbLocation2[xx] << " ";
				}
		}
			for (int xx = 0; xx < Nx; xx++)
			{
				SlugFlow[xx] = SlugFlow2[xx];
				HLTbLocation[xx] = HLTbLocation2[xx];
				//std::cout << SlugFlow[xx] << " ";
			}
		
	}

	//In this method, the requried multi-dimensional vectors are resized as required 
	void SlugFLowHeatTransferCal::VectorResize()
	{

		hHotSave.clear();
		Es.clear();
		Er.clear();
		Ex.clear();
		ETheta.clear();
		KDepoSave.clear();
		riEff.clear();
		EWf.clear();
		EWG.clear();
		EWS.clear();
		ECF.clear();
		ECG.clear();
		ECS.clear();
		Tb.clear();
		Kb.clear();
		TwAv.clear();
		Qtot.clear();
		hBulk.clear();
		ApEff.clear();
		EwfLocal.clear();
		EwgLocal.clear();
		EWsLocal.clear();
		AL.clear();
		AG.clear();		
		TLActual.clear();
		TGActual.clear();

		TLActual.resize(Nx + 1);
		TGActual.resize(Nx + 1);
		EWf.resize(Nx + 1);
		EWS.resize(Nx + 1);
		EWG.resize(Nx + 1);
		ECF.resize(Nx + 1);
		ECG.resize(Nx + 1);
		ECS.resize(Nx + 1);
		Tb.resize(Nx + 1);
		Kb.resize(Nx + 1);
		TwAv.resize(Nx + 1);
		Qtot.resize(Nx + 1);
		hBulk.resize(Nx + 1);
		ApEff.resize(Nx + 1);

		EwfLocal.resize(RingSec + 1);
		EwgLocal.resize(RingSec + 1);
		EWsLocal.resize(RingSec + 1);
		AL.resize(RingSec + 1);
		AG.resize(RingSec + 1);
		

		hHotSave.resize(Nx + 1);
		KDepoSave.resize(Nx + 1);
		riEff.resize(Nx + 1);
		Es.resize(Nx + 1);
		Er.resize(Nx + 1);
		Ex.resize(Nx + 1);
		ETheta.resize(Nx + 1);			
			for (int xx = 0; xx < Nx; xx++)
			{
				KDepoSave[xx].resize(RingSec + 1);
				riEff[xx].resize(RingSec + 1);
				Es[xx].resize(WallLay + 1);
				Er[xx].resize(WallLay + 1);
				Ex[xx].resize(WallLay + 1);
				ETheta[xx].resize(WallLay + 1);
			
				for (int WL = 1; WL < WallLay + 1; WL++) //starts from 1 to 3
				{
					Es[xx][WL].resize(RingSec + 1);
					Er[xx][WL].resize(RingSec + 1);
					Ex[xx][WL].resize(RingSec + 1);
					ETheta[xx][WL].resize(RingSec + 1);
				}
			}
	}

	void SlugFLowHeatTransferCal::PipeSpecs(const  user::PipeGeometry& PipeProps) {
		IDPipe = PipeProps.IDPipe;
		OD = PipeProps.ODPipe;
		KPipe = PipeProps.KPipe;
		CPPipe = PipeProps.CPPipe;
		RhoPipe = PipeProps.RhoPipe;
	}
		
	double SlugFLowHeatTransferCal::Aout(int WallLayer)
	{
		double ans;
		ans = M_PI * 0.5*(ro - (WallLay - double(WallLayer))*w / 3.0)*dX;
		return(ans);
	}

	double SlugFLowHeatTransferCal::Ain(int WallLayer, double WaxThickness) {
		double ans;
		if (WallLayer = !1)
		{
			ans = M_PI * 0.5*(ro - (WallLay - double(WallLayer) + 1)*w / 3.0)*dX;
		}
		else {
			ans = M_PI * 0.5*(ro - w - WaxThickness)*dX;
		}
		return(ans);
	}

	double SlugFLowHeatTransferCal::Ax(int WallLayer, double WaxThickness) {
		
		double ans;
		if (WallLayer = !1)
		{
			ans = M_PI * 0.5*(ro - (WallLay - double(WallLayer))* (w / 3.0) - w / (2 * 3.0))*w / 3.0;//
		}
		else {
			ans = M_PI * 0.5*(ro - w - WaxThickness / 2.0)*WaxThickness;
		}
		return(ans);
	}

	double SlugFLowHeatTransferCal::ATheta(int WallLayer, double WaxThickness) {
		double ans;
		if (WallLayer = !1)
		{
			ans = w / 3.0*dX;
		}
		else {
			ans = WaxThickness * dX;
		}
		return(ans);
	}

	double SlugFLowHeatTransferCal::dxTheta(int WallLayer, double WaxThickness)
	{
		double ans;
		if (WallLayer = !1)
		{
			ans = M_PI * 0.5*(ro - (WallLay - double(WallLayer))* (w / 3.0) - w / (2 * 3.0));
		}
		else {
			ans = M_PI * 0.5*(ro - w - WaxThickness / 2.0);
		}
		return(ans);
	}

	double SlugFLowHeatTransferCal::WaxDepoK(double T,double Fw)
	{
		double Kdepo;
		//This updates the thermal conductivity of the liquid phase
		Update(T, HLLS);
		Kdepo = KL * (1 + (3 * Fw) / ((KWax + 2 * KL) / (KWax - KL) - Fw));
		return(Kdepo);
	}

	double SlugFLowHeatTransferCal::IDEff(std::vector<std::vector<double>>& WaxVector)
	{
		double ans;
		double sum=0;
		double sum2=0;
		double AvgThickness;
		for (int j = 0; j < Nx; j++) 
		{
			sum = 0;
			for (int i = 1; i < RingSec +1; i++)
			{
				sum = sum + WaxVector[j][i];
			}
			sum2 = sum2 + sum / 4.0;
		}
		AvgThickness = sum2 / double(Nx);
		ans = IDPipe - 2 * AvgThickness;
		return(ans);
	}

} //MultiphaseTemperature namespace






