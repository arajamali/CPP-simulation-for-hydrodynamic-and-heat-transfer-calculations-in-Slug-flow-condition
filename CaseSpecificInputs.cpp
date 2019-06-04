#include <iostream>
#include <cmath>
#include "CaseSpecificInputs.h"

namespace user {
	FluidProperties::FluidProperties(double Temp, double Press) {
		T = Temp + 273.15;
		P = Press;
	}

	//Liquid properties
	void FluidProperties::Liquid::SetContainer(FluidProperties &container)
	{
		this->container = &container;
	}
	void FluidProperties::Liquid::Render(void)
	{
		T = container->T;
		P = container->P;
	}

	//Checked
	double FluidProperties::Liquid::Cp(FluidProperties &container)
	{
		this->container = &container;
		Render();
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(2.0819*TF + 1.8858 * 1000); //very good accuracy
	} //T should be Kelvin, P is 350 Psig

	  //Checked
	double FluidProperties::Liquid::Rho(FluidProperties &container)//At 350 Psig ONLY
	{//[kg/m3]
		this->container = &container;
		Render();
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(-3.963*0.1*TF + 8.586 * 100);
	}

	//Checked
	double FluidProperties::Liquid::Mu(FluidProperties &container)
	{
		//kg/(m s)
		this->container = &container;
		Render();
		double TR = T * 9.0 / 5.0;//Rankin
		double a = -3.9932*log10(TR) + 10.7408;
		double ans = (Rho(container) / 1000.0)*(pow(10, pow(10, a)) - 0.7);
		return(ans / 1000.0);
	}

	//Surface tension between gas and liquid (RR)
	//Not reliable
	double FluidProperties::Liquid::STG(FluidProperties &container)
	{
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(pow(10, -3)*(-1.5667*pow(10, -2)*TF + 1.7666 * 10));
	}

	//Thermal conductivit of liquid phase
	double FluidProperties::Liquid::K(FluidProperties &container)
	{
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(-2.998*pow(10, -4.0)*TF + 2.0954*0.1);
	}

	void FluidProperties::Gas::SetContainer(FluidProperties &container)
	{
		this->container = &container;
	}
	void FluidProperties::Gas::Render(void)
	{
		T = container->T;
		P = container->P;
	}

	//Checked!
	double FluidProperties::Gas::Z(FluidProperties &container)
	{
		this->container = &container;
		Render();
		double gamma = MW / 28.9625;
		double Ppsi = (P) / 6894.7572932;//PsiA
		Ppsi = Ppsi - 14.7;//psig
		return(1 / (1 + (Ppsi * 344400 * pow(10, 1.785*gamma)) / (pow(T * 9.0 / 5.0, 3.825))));
	}

	//Checked 
	double FluidProperties::Gas::Cp(FluidProperties &container)
	{
		this->container = &container;
		Render();
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(8.2976*0.1*TF + 2.2822 * 1000);
	}

	//Checked
	double FluidProperties::Gas::Rho(FluidProperties &container) // [kg / m3]
	{//[J/kg/K]
		this->container = &container;
		double ans;
		Render();
		double Ppsi = P / 6894.7572932;
		double R = 10.7316;//psi.ft3/mol/R
		double TR = T * 9.0 / 5.0;//Rankin
		ans = Z(container);
		return(Ppsi*MW / (Z(container)*R*TR) / 62.4279606 * 1000);
	}
	//Not prefered, use RhoG instead
	double FluidProperties::Gas::Rho2(FluidProperties &container)
	{
		this->container = &container;
		Render();
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(-4.3*0.01*TF + 2.23 * 10);
	}

	//Checked
	double FluidProperties::Gas::Mu(FluidProperties &container)
	{//[kg/(m s)]
		this->container = &container;
		Render();
		double TR = T * 9.0 / 5.0;//Rankin
		double X = 3.448 + 986.4 / (TR)+0.01009*MW;
		double Y = 2.447 - 0.2224*X;
		double K = (9.379 + 0.01607*MW)*pow(TR, 1.5) / (209.2 + 19.26*MW + TR);
		double ans = pow(10, -7.0)*K*exp(X*pow(Rho(container)*0.001, Y));
		return(ans);
	}

	double FluidProperties::Gas::K(FluidProperties &container)
	{
		double TF = (9.0 / 5.0)*(T - 273.15) + 32;
		return(5.963*pow(10, -5)*TF + 3.2221*0.01);
	}

}//usr

