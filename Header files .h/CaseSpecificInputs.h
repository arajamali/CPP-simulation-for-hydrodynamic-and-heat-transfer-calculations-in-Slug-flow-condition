#ifndef CaseSpecificInputs_H
#define CaseSpecificInputs_H
#include<cmath>
#include <iostream>

namespace user {
	class FluidProperties {
		double T;
		double P;
		double MW = 17.68358;

	public:

		FluidProperties(double T, double P);

		class Liquid {

		public:
			double T;
			double P;
			double MW;
			FluidProperties *container;

			void SetContainer(FluidProperties &container);

			void Render(void);
      //T should be Kelvin, P is 350 Psig
			double Cp(FluidProperties &container); 
      //At 350 Psig ONLY
			double Rho(FluidProperties &container);
			double Mu(FluidProperties &container);
			double STG(FluidProperties &container);
			double K(FluidProperties &container);


		};

		class Gas {

		public:
			double T;
			double P;
      //In this parameter is manually put here. (Be carefull)
			double MW = 17.68358; 
			FluidProperties *container;

			void SetContainer(FluidProperties &container);
			void Render(void);
			double Z(FluidProperties &container);
			double Cp(FluidProperties &container);
			double Rho(FluidProperties &container);
      //At 350 psig pressure ONLY
			double Rho2(FluidProperties &container);
			double Mu(FluidProperties &container);
			double K(FluidProperties &container);

		};

		Liquid Liquid;
		Gas Gas;

	};

	class PipeGeometry {

	public:

		double Roughness = 0.0018*0.0254; //Meter
		double IDPipe = 2.067*0.0254; //Meter
		double ODPipe = 2.375*0.0254; //Meter
		double KPipe = 16;//w/m/k
    //https://www.engineersedge.com/materials/specific_heat_capacity_of_metals_13259.htm
		double CPPipe = 502.416; //J/Kg/K 
    //https://www.engineeringtoolbox.com/steel-pipes-dimensions-d_43.html
		double RhoPipe = 7835.34367; //Kg/m3 (3.65 lb/ft) / (0.007462024 ft2) 
		double IDJacket = 3.826*0.0254;
		double ODJacket = 4.5*0.0254;
		//Seems test length
		double Length = 0.3048 * 53;
	};
}//namespace user {

#endif
