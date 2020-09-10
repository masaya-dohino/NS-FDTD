#ifndef _LIGHT_SOURCE_H
#define _LIGHT_SOURCE_H

#include "Object.h"
#include "PhysicalConstant.h"
#include "Field.h"

class LightSource{
	const int angle;
	const double lambda_s;	//”g’·
	const double k_s;			//”g”
	const double w_s;			//Špü”g”
	const double T_s;			//üŠú
public:
	LightSource(int ang, double lambda);
	int getAngle();
	double getAngleRad();
	Complex lightSource(int time);

	//void pointSource(complex<double> *p);
	//void linearSource(complex<double> *p);
};

inline int LightSource::getAngle(){
	return angle;
}

inline double LightSource::getAngleRad(){
	return M_PI*angle/180.0;
}

/*
//“_ŒõŒ¹‚Ì“üË
void LightSource::pointSource(complex<double> *p){
	double ray_coef = exp(-0.0001*(time-500)*(time-500));
	p[index(mField->getNx()/2, mField->getNy()/2)] = ray_coef*polar(1.0, w_s*time);
}

//üŒõŒ¹‚Ì“üË
void LightSource::linearLightSource(complex<double> *p){

	double ray_coef = exp(-0.0001*(time-500)*(time-500));
	for(int i=1; i<mField->getNy()-1; i++)
		p[index(5,i, +1)]  +=ray_coef*polar(1.0, w_s*time);	//todo ‹«ŠEã‚É‚à“ü‚ê‚Ä‚¢‚¢‚Ì‚©?
}

*/

#endif //_LIGHT_SOURCE_H