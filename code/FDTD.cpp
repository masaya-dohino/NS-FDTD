#include "FDTD.h"
static double min_lambda = 200*_pow(10,-9);		//’²‚×‚é”g’·‚Ì”ÍˆÍ
static double max_lambda = 800*_pow(10,-9);

FDTD::FDTD()
:Solver()
{
	phi  = new complex<double>[3*mField->getNcel()];		//—ÌˆæŠm•Û
	np   = new double[mField->getNcel()];		//ŒvZ—p’è”

	for(int i=0; i<3*mField->getNcel(); i++)	phi[i] = 0;				//‰Šú‰»

	for(int i=0; i<mField->getNcel(); i++)	np[i] = 0;

	cout << "FDTD Constructor" << endl;
}

FDTD::~FDTD(){
	delete[] phi;
	delete[] np;
	cout << "FDTD Destructor" << endl;
}

void FDTD::Initialize(){
	for(int i=0; i<3*mField->getNcel(); i++)	phi[i] = 0;				//‰Šú‰»
	time = 0;
}

void FDTD::Initialize(double _lambda)
{
	for(int i=0; i<3*mField->getNcel(); i++)
		phi[i] = 0;				//‰Šú‰»
	time = 0;
}

bool FDTD::calc(){
	return true;
}

//•`‰æ
void FDTD::draw(){
	Solver::draw(phi);
}