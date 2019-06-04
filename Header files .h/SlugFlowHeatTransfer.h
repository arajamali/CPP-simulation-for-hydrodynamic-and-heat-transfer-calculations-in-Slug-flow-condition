#ifndef SlugFlowHeatTransfer_H
#define SlugFlowHeatTransfer_H
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "CaseSpecificInputs.h"
#include "StratifiedFlowGeometry.h"
#include "SlugFlowFilmProfile.h"
#include "dPCalculation.h"
#include "OveralHeatTransfer.h"
namespace MultiphaseTemperature {
	//For those flow conditions where Lu is large, bigger Nx can help the accuracy
	class SlugFLowHeatTransferCal {
	public:
		//static const int NXSize = 50;
		double NuNumCal;
		//double NuNumCheck;
		void SlugFlowHeatModel(std::vector<double>& FilmInputs,
			                     std::vector<std::vector<double>>& WaxWidthVector, 
                           std::vector<std::vector<double>>& WaxSolidFraction,
			                     std::vector<double>& HydroInputs, 
                           std::vector<std::vector<std::vector<double>>>& TwInitial,
			                     std::vector<double>& Ts, std::vector<double>& Tf, 
                           std::vector<double>& Tg, std::vector<double>& FilmAssign,
			                     std::vector<bool>& SlugAssign, 
                           std::vector<double>& q, double TInlet);

	
		//static const int NTime = 1000;
		static const int Nx = 50;
		static const int WallLay = 4;
		static const int RingSec = 4;

		//Pipe's inlet temperature
		//In the case, the Inlet temperature is assumed to he same as the initial condition
		//Number of sections in slug and in film sections
		double UInt;
		double AvgTb=0;
		double AvgTempForHydroUpdate;
		double NxS, NxF;
		double g = 9.8;
		double RhoL, RhoG, MuL, MuG, KG, 
           KL, CpL, CpG, IT, KS, RhoS, MuS, CpS;
		double VSG, VSL, VTb, VLLS, VLTB, 
           VGTB, VLTB2, VGTB2, VGLS
				, VL, VG, VL2, VG2, Vm;
		double T, AvgTemp,P;
		double Cg, Cs;
		double Ls, Lf, Lu;
		double ECF1, ECF2, ECG1, ECG2, ECS1, ECS2;
		double dL, dG, SI, AGP, ALP;
		double HLTBAvg, hFAvg,HLTB2,HLLS;
		double thetahF;
		double w, KPipe,KThermal, CPPipe, RhoPipe, 
           PipeRoughness, CpPipe, IDPipe,IDEf, OD, ro, ri;
		double dX, dt;
		double PressGas, PressFilm;
		double KGSave, KLSave;
		double sumWF = 0, sumWG = 0, sumWS = 0,
           QSSum = 0, QFSum = 0, QGSum = 0,
			     SumhCheck = 0, TSSum = 0, TfSum = 0;
		double WLayD = double(WallLay);
		double hFilmGas, hSlug, hFilmLiquid;
		bool Flag = 1;
		double RhoDepo = 855; //Assumed
		double KWax=0.25;//Assumed

		std::vector<double> EWf, EWG, EWS, ECF, ECG,ECS, Tb, Kb, TwAv, Qtot, hBulk,
                        EwfLocal, EwgLocal, EWsLocal,AL, AG,ApEff;
		double check;
		double RhoJacket = 1080; //should be checked, it is just a prototype
		double QCoolant=0.004600455; //This should be input, it is just a protoype
		double WcJacket; 

		std::vector<double> q;
		std::vector<double> TS;
		std::vector<double> TL,TLActual;
		std::vector<double> TG, TGActual;
		std::vector<std::vector<double>> WaxThickness, WaxSolidFrac;
		std::vector<double> xAxisNew, HLTBOrgNew;
		std::vector<std::vector<std::vector<double>>>Es, Er, Ex, ETheta;
		std::vector<std::vector<std::vector<double>>>Tw;
		std::vector<std::vector<double>>KDepoSave;
		std::vector<double> hHotSave;
		std::vector<std::vector<double>>riEff;
		std::vector<double> HLTbLocation;
		std::vector<bool> SlugFlow;
			
		//std::ofstream TLOut, TGOut, TW1Out, TW2Out, TW3Out, NuCal, NuCheck,
		//	TW4Out, TW5Out, TW6Out, TW7Out, TW8Out,TWax1, TWax2, TWax3,
		//	TW9Out, Time, Axial, Tbulk;

		/*Input variables are read from the input vectors and are stored to 
		appropriate class parameters*/
		void SetVariablesFromIputFiles(std::vector<double>& FilmInputs,
			                             std::vector<std::vector<double>>& WaxWidthVector, 
                                   std::vector<std::vector<double>>& WaxSolidFraction,
			                             std::vector<std::vector<std::vector<double>>>& TwInitial,
			                             std::vector<double>& Ts, std::vector<double>& Tf, 
                                   std::vector<double>& Tg, std::vector<double>& FilmAssign,
			                             std::vector<bool>& SlugAssign);


		//Hydrodynamic properties are calculated and assigned to class variables
		void HydroUpdate(std::vector<double>&FimInputs,std::vector<double>& HydroInputs);

		/*Convective heat transfer coefficient is calculated. hL,
		hG and hS are calculated in the program*/
		double HeatCoe(double K, double d, double Rho, double V, double Mu, double Cp);

		/*Distribution coefficient is calcualted. This parameter is different for
		laminar and for turbulent flow */
		double DistCoe(double Rho, double V, double d, double Mu);

		/*In this method, fluid properties are updated for a given Temperature 
		and slug holdup*/
		void Update(double T, double HLLS);

		/* Nusselt number for liquid single phase flow is calculated for Two-Phase 
		  Nusselt number coefficient correlation */
		double SinglePhaseNu(double Rho, double v, double D, double Mu,
			double Cp, double K, double dx, double E, double H);



		//AG and AL for Er calculation and EwLocal calcualtions
		void SurfaceAreaAGAL(double Theta, double Top, double Side, double Bot);

		//Stratified flow geormetries plus some other related values are calculated
		void UpdateFlowGeomParams(const StratifiedFlowGeometry& FlowGeom, double HL);

		/*The slugFlow boolean is set appropriately for each time step from the
		following two methods*/
		void SlugFlowBooleanInitial();
		void SlugFlowBooleanLater(std::vector<bool> SlugFlowPre,            
                              std::vector<double> HLTbLocationPre);

		//Corrct sizes are assigned to the vectors
		void VectorResize();



		//The pipe's specifications are obtained through this method
		void PipeSpecs(const  user::PipeGeometry& PipeProps);
		
		//Aout in the calcualtions

		//Ax, ArInn, ArOut, ATheta are calculated in this method
		double Aout(int WallLayer);
		double Ain(int WallLayer, double WaxThickness);
		double Ax(int WallLayer, double WaxThickness);
		double ATheta(int WallLayer, double WaxThickness);
		double dxTheta(int WallLayer, double WaxThickness);

		double WaxDepoK(double T, double Fw);
		double IDEff(std::vector<std::vector<double>>& WaxVector);
	};




}//namespace MultiphaseTemperature 
#endif