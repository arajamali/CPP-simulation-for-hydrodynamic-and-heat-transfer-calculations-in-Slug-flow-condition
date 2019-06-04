#ifndef SlugFlowFilmProfile_H
#define SlugFlowFilmProfile_H
#include <cmath>
#include <iostream>
#include <vector>
#include "CaseSpecificInputs.h"
#include "StratifiedFlowGeometry.h"
namespace MultiphaseTemperature {

	//Please note the following class is for horizontal case
	//however, it can easily be modified for inclined case

	class FilmProfile {
	private:

		double Pi = 3.14159265358979;//Pi number
		double ID, OD, VSG, VSL, T, P, PipeRoughness, RhoL, RhoG, MuL, MuG,C0, g, RhoS, MuS, IT;
		
		int k;
		int Con;
		double Tol = 0.001;
		int Flag2;
		double hFdZ(double hF);
		double hFDetermination(double Lf, double hMax, double EpsilonZh);
		double LfDetermination(double InitialGuess, double dx, double Err,
                           double hMax, double EpsilonZh);
	public:
		//These arrays are used to store hf and z values
		std::vector<double>zArray, hFArray;
		//double hFArray[1000], zArray[1000];
		//Final error
		double err = 10;
		double VLTB, VGTB, VLLS, VGLS, Vm, VTb;
		double HLLS, LfFinal, Ls, HLTB,HLTBAvg;
		double  hFAvg;
		int ZSize;

		FilmProfile(std::vector < double > & FilmProfileInputs);
		double FannningFrictionFactor(double Rho, double V, double d, double Mu, 
                                  double Roughness);
		void HyrodynamicProperties(double InitialGuess, double dx, double Epsilon,
                               double hMax, double EpsilonZh);
		double hFCritical(double hF);
		double hFCriticalDetermination();
		//void ReportParams();
		friend class SlugFLowHeatTransferCal;
		std::vector<double> XSave;
		std::vector<double> HLTBSave;
		double CalculatedError;
		double AvgFilmHoldUp_HLTB;
		double AvgSlugHoldUp_HLLS;
		double FilmVelocity_VLTB;
		double GasCoreVelocity_Vg;
		double UnitSlugLength_Lu;
		double FilmLength_Lf;
		double SlugLength_Ls;




	};

}//namespace MultiphaseTemperature 
#endif
