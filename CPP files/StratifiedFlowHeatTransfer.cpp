/*In this cpp file, different fluid properties are
calculated based on temperature and pressure*/

//Glycol properties are not included
//Properties of gas and oil are listed here


#include "StratifiedFlowGeometry.h" //It provides hydraulic parameters of both phases
#include "CaseSpecificInputs.h" //Pipe and fluid properties
#include <iostream>
#include <cmath>

using namespace std;


SFConstWallTemp::SFConstWallTemp(double Length, double SuperVelocityGas,
                                 double SuperVelocityLiquid, double LiqHoldUp,
                                 double TempWall, double TempInlet, double Pressure)
{
	L = Length, TW = TempWall; TLI = TempInlet; VSG = SuperVelocityGas;
  VSL = SuperVelocityLiquid; P = Pressure; h = LiqHoldUp;
}
void SFConstWallTemp::Liquid::SetContainer(SFConstWallTemp &container)
{
	this->container = &container;
}
void SFConstWallTemp::Liquid::Render(void)
{
	VSL = container->VSL;
	L = container->L;
	TW = container->TW;
	TLI = container->TLI;
	P = container->P;
	h = container->h;
}


double SFConstWallTemp::Liquid::Nu(double T, SFConstWallTemp &container) {
	this->container = &container;
	Render();
	PipeGeometry PipeGeometry;
  // Temp and pressure are needed for properties class
	FluidProperties Properties(T, P); 
	StratifiedFlowGeometry FlowGeo(h);
	ID = PipeGeometry.IDPipe;
	double ans;
	double AL = (FlowGeo.ALTilt)*ID*ID;
	double dL = pow(AL * 4 / Pi, 0.5);
	double Rho = Properties.Liquid.Rho(Properties);
	double Mu = Properties.Liquid.Mu(Properties);
	double Cp = Properties.Liquid.Cp(Properties);;
	double K = 0.15;//For now ******
	double VS = VSL / h;
	double Re = dL * VS*Rho / Mu;
	double Pr = Cp * Mu / K;
	ans = 0.023*pow(Re, 0.8)*pow(Pr, 1.0 / 3.0);
	return(ans);
}

double SFConstWallTemp::Liquid::TempProfile(SFConstWallTemp &container) {
	this->container = &container;
	Render();
	StratifiedFlowGeometry FlowGeo(h);
	PipeGeometry PipeGeometry;
	ID = PipeGeometry.IDPipe;
	double X = L / 200.0;
	double hAverage;
	double NuL;
	double T = TLI;
	double CP;
	double RHO;
	double SL = FlowGeo.SLTilt*ID;
	double AL = (FlowGeo.ALTilt)*ID*ID;
	double dL = pow(AL * 4 / Pi, 0.5);
	double W;
	double TLO = T - 5;//Initial guess
	double F, Fdx, Fprime;
	double dx = 0.0000001;
	double Eps = 10;
	for (int i = 0; i < 200; i++) {

		//This while basically does NR
		while (Eps > 10E-12) {
			FluidProperties Properties(T, P);
			CP = Properties.Liquid.Cp(Properties);
			RHO = Properties.Liquid.Rho(Properties);
			W = (VSL / h)*AL*RHO;
			NuL = Nu(T, container);
			hAverage = NuL * 0.15 / dL;
			F = hAverage * SL*X*((TW - T) - (TW - TLO)) / log((TW - T) /
         (TW - TLO)) - W * CP*(TLO - T);
			Eps = F;
			Fdx = hAverage * SL*X*((TW - T) - (TW - (TLO + dx))) /
         log((TW - T) / (TW - (TLO + dx))) - W * CP*((TLO + dx) - T);
			Fprime = (Fdx - F) / dx;
			TLO = TLO - F / Fprime;
		}
		T = TLO;
		cout << TLO << endl;
		TLO = T - 1;
		Eps = 10;
	}
	return(TLO);
}

void SFConstWallTemp::Gas::SetContainer(SFConstWallTemp &container)
{
	this->container = &container;
}
void SFConstWallTemp::Gas::Render(void)//
{

	VSG = container->VSG;
	L = container->L;
	TW = container->TW;
	TLI = container->TLI;
	P = container->P;
	h = container->h;

}

double SFConstWallTemp::Gas::Nu(double T, SFConstWallTemp &container) {
	this->container = &container;
	Render();
	PipeGeometry PipeGeometry;
  // Temp and pressure are needed for properties class
	FluidProperties Properties(T, P); 
	StratifiedFlowGeometry FlowGeo(h);
	ID = PipeGeometry.IDPipe;
	double ans;
	double AG = (FlowGeo.AGTilt)*ID*ID;
	double dG = pow(AG * 4 / Pi, 0.5);
	double Rho = Properties.Gas.Rho(Properties);
	double Mu = Properties.Gas.Mu(Properties);
	double Cp = Properties.Gas.Cp(Properties);
	double K = 0.035854;//For now ******
	double VG = VSG / h;
	double Re = dG * VG*Rho / Mu;
	double Pr = Cp * Mu / K;
	ans = 0.023*pow(Re, 0.8)*pow(Pr, 1.0 / 3.0);
	return(ans);
}

double SFConstWallTemp::Gas::TempProfile(SFConstWallTemp &container) {
	this->container = &container;
	Render();
	StratifiedFlowGeometry FlowGeo(h);
	PipeGeometry PipeGeometry;
	ID = PipeGeometry.IDPipe;
	double X = L / 200.0;
	double hAverage;
	double NuG;
	double T = TLI;
	double CP;
	double RHO;
	double SG = FlowGeo.SGTilt*ID;
	double AG = (FlowGeo.AGTilt)*ID*ID;
	double dG = pow(AG * 4 / Pi, 0.5);
	double W;
	double TLO = T - 5;//Initial guess
	double F, Fdx, Fprime;
	double dx = 0.0000001;
	double Eps = 10;
	for (int i = 0; i < 200; i++) {

		//This while basically does NR
		while (Eps > 10E-12) {
			FluidProperties Properties(T, P);
			CP = Properties.Gas.Cp(Properties);
			RHO = Properties.Gas.Rho(Properties);
			W = (VSG / h)*AG*RHO;
			NuG = Nu(T, container);
			hAverage = NuG * 0.15 / dG;
			F = hAverage * SG*X*((TW - T) - (TW - TLO)) / log((TW - T) / 
            (TW - TLO)) - W * CP*(TLO - T);
			Eps = F;
			Fdx = hAverage * SG*X*((TW - T) - (TW - (TLO + dx))) / 
            log((TW - T) / (TW - (TLO + dx))) - W * CP*((TLO + dx) - T);
			Fprime = (Fdx - F) / dx;
			TLO = TLO - F / Fprime;
		}
		T = TLO;
		cout << TLO << endl;
		TLO = T - 1;
		Eps = 10;
	}
	return(TLO);
}


SFConstHeatFlux::
    SFConstHeatFlux(double Length, double SuperVelocityGas, 
                    double SuperVelocityLiquid, double LiqHoldUp,
                    double InitialTempWall, double TempInlet, double Pressure)
{
	L = Length, TW = InitialTempWall; TLI = TempInlet; VSG = SuperVelocityGas;
  VSL = SuperVelocityLiquid; P = Pressure; h = LiqHoldUp;
}
void SFConstHeatFlux::Liquid::SetContainer(SFConstHeatFlux &container)
{
	this->container = &container;
}
void SFConstHeatFlux::Liquid::Render(void)
{
	VSL = container->VSL;
	L = container->L;
	TWI = container->TW;
	TLI = container->TLI;
	P = container->P;
	h = container->h;
}


double SFConstHeatFlux::Liquid::Nu(double T, SFConstHeatFlux &container) {
	this->container = &container;
	Render();
	PipeGeometry PipeGeometry;
  // Temp and pressure are needed for properties class
	FluidProperties Properties(T, P); 
	StratifiedFlowGeometry FlowGeo(h);
	ID = PipeGeometry.IDPipe;
	double ans;
	double AL = (FlowGeo.ALTilt)*ID*ID;
	double dL = pow(AL * 4 / Pi, 0.5);
	double Rho = Properties.Liquid.Rho(Properties);
	double Mu = Properties.Liquid.Mu(Properties);
	double Cp = Properties.Liquid.Cp(Properties);;
	double K = 0.15;//For now ******
	double VS = VSL / h;
	double Re = dL * VS*Rho / Mu;
	double Pr = Cp * Mu / K;
	ans = 0.023*pow(Re, 0.8)*pow(Pr, 1.0 / 3.0);
	return(ans);
}

double SFConstHeatFlux::Liquid::TempProfile(SFConstHeatFlux &container) {
	this->container = &container;
	Render();
	StratifiedFlowGeometry FlowGeo(h);
	PipeGeometry PipeGeometry;
	ID = PipeGeometry.IDPipe;
  //This is one of the assumptions of constant heat flux approach
	double ConstTwTo = TWI - TLI;
	double X = L / 200.0;
	double hAverage;
	double NuL;
	double T = TLI;
	double AL = (FlowGeo.ALTilt)*ID*ID;
	double dL = pow(AL * 4 / Pi, 0.5);
	double TW = TWI;
	double CP;
	double RHO;
	double SL = FlowGeo.SLTilt*ID;
	double W;
	double ans;

	for (int i = 0; i < 200; i++) {

		FluidProperties Properties(T, P);
		CP = Properties.Liquid.Cp(Properties);
		RHO = Properties.Liquid.Rho(Properties);
		W = (VSL / h)*AL*RHO;
		NuL = Nu(T, container);
		hAverage = NuL * 0.15 / dL;
		ans = hAverage * SL*X*ConstTwTo / (W*CP);
    // Plus sign here, I am not so sure!!
		T = T + hAverage * SL*X*ConstTwTo / (W*CP); 
		TW = T + ConstTwTo;
		cout << T << "	" << TW << endl;

	}
	return(T);

}

void SFConstHeatFlux::Gas::SetContainer(SFConstHeatFlux &container)
{
	this->container = &container;
}
void SFConstHeatFlux::Gas::Render(void)
{

	VSG = container->VSG;
	L = container->L;
	TWI = container->TW;
	TLI = container->TLI;
	P = container->P;
	h = container->h;

}

double SFConstHeatFlux::Gas::Nu(double T, SFConstHeatFlux &container) {
	this->container = &container;
	Render();
	PipeGeometry PipeGeometry;
  // Temp and pressure are needed for properties class
	FluidProperties Properties(T, P); 
	StratifiedFlowGeometry FlowGeo(h);
	ID = PipeGeometry.IDPipe;
	double ans;
	double AG = (FlowGeo.AGTilt)*ID*ID;
	double dG = pow(AG * 4 / Pi, 0.5);
	double Rho = Properties.Gas.Rho(Properties);
	double Mu = Properties.Gas.Mu(Properties);
	double Cp = Properties.Gas.Cp(Properties);;
	double K = 0.035854;//For now ******
	double VG = VSG / h;
	double Re = dG * VG*Rho / Mu;
	double Pr = Cp * Mu / K;
	ans = 0.023*pow(Re, 0.8)*pow(Pr, 1.0 / 3.0);
	return(ans);
}

double SFConstHeatFlux::Gas::TempProfile(SFConstHeatFlux &container) {
	this->container = &container;
	Render();
	StratifiedFlowGeometry FlowGeo(h);
	PipeGeometry PipeGeometry;
	ID = PipeGeometry.IDPipe;
  //This is one of the assumptions of constant heat flux approach
	double ConstTwTo = TWI - TLI;
	double X = L / 200.0;
	double hAverage;
	double NuG;
	double T = TLI;
	double AG = (FlowGeo.AGTilt)*ID*ID;
	double dG = pow(AG * 4 / Pi, 0.5);
	double TW = TWI;
	double CP;
	double RHO;
	double SG = FlowGeo.SGTilt*ID;
	double W;
	double ans;

	for (int i = 0; i < 200; i++) {

		FluidProperties Properties(T, P);
		CP = Properties.Gas.Cp(Properties);
		RHO = Properties.Gas.Rho(Properties);
		W = (VSG / h)*AG*RHO;
		NuG = Nu(T, container);
		hAverage = NuG * 0.15 / dG;
		ans = hAverage * SG*X*ConstTwTo / (W*CP);
    // Plus sign here, I am not so sure!!
		T = T + hAverage * SG*X*ConstTwTo / (W*CP); 
		TW = T + ConstTwTo;
		cout << T << "	" << TW << endl;
	}
	return(T);
}
