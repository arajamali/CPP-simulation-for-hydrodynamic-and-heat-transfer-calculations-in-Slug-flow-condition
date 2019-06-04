#ifndef StratifiedFlowGeometry_H
#define StratifiedFlowGeometry_H
#include<cmath>
#include <iostream>
#include "CaseSpecificInputs.h"

namespace MultiphaseTemperature {
	class StratifiedFlowGeometry {
		//hL should be known. hL is the height of liquid in the pipe
		//Always ID is equal or larger than hL

		double hL;
		double Pi = 3.14159265358979;//Pi number


	public:

		StratifiedFlowGeometry(double h, double IDEffective);

		double ID;
		double A;
		double hLTilt;
		double SGTilt;
		double SLTilt;
		double SITilt;
		double AGTilt;
		double ALTilt;
		double dL;
		double dG;
		double Theta;
		friend class PipeGeometry;
	};
}

#endif
