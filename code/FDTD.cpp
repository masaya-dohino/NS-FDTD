#include "FDTD.h"
static double min_lambda = 200*_pow(10,-9);		//調べる波長の範囲
static double max_lambda = 800*_pow(10,-9);

FDTD::FDTD()
:Solver()
{
	phi  = new complex<double>[3*mField->getNcel()];		//領域確保
	np   = new double[mField->getNcel()];		//計算用定数

	for(int i=0; i<3*mField->getNcel(); i++)	phi[i] = 0;				//初期化

	for(int i=0; i<mField->getNcel(); i++)	np[i] = 0;

	cout << "FDTD Constructor" << endl;
}

FDTD::~FDTD(){
	delete[] phi;
	delete[] np;
	cout << "FDTD Destructor" << endl;
}

void FDTD::Initialize(){
	for(int i=0; i<3*mField->getNcel(); i++)	phi[i] = 0;				//初期化
	time = 0;
}

void FDTD::Initialize(double _lambda)
{
	for(int i=0; i<3*mField->getNcel(); i++)
		phi[i] = 0;				//初期化
	time = 0;
}

bool FDTD::calc(){
	return true;
}

//描画
void FDTD::draw(){
	Solver::draw(phi);
}