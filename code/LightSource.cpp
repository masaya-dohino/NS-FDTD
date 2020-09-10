#include "LightSource.h"

LightSource::LightSource(int ang, double lambda)
:angle(ang), lambda_s(lambda), k_s(2*M_PI/lambda_s),w_s(LIGHT_SPEED_S*k_s), T_s(lambda_s/LIGHT_SPEED_S)
{
}
