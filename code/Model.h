#ifndef _MODEL_H
#define _MODEL_H
#include"PhysicalConstant.h"
#include<string>
#include<direct.h>
#include"Object.h"
#include"Field.h"
#include<iostream>

using namespace std;
enum INTEG{
	D_XY = 0,
	D_X,
	D_Y
};

class FazzyModel {
protected:
	Field * mField;

public:
	FazzyModel(Field *field) {
		mField = field;
	};
	virtual string mkdir(string root) = 0;
	virtual double calcEPS(const double&, const double&, enum INTEG = D_XY) = 0;
	virtual double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) = 0;
	virtual bool update(int) = 0;
	virtual void Initialize() = 0;
	virtual void InitializeLambda(double a) {
	}
};

class FazzySlabModel :public FazzyModel{
	const double ep1, ep2;
	const int width1, width2;
public:
	FazzySlabModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	bool update(int){
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyMieModel :public FazzyModel{
	double ep;
	double r;
public:
	FazzyMieModel(Field* f);
	string mkdir(string root);	//ディレクトリ作成
	virtual void InitializeLambda(double a) {
		this->r = mField->nanoToCell(a);
		cout << "r=" + to_s((int)mField->cellToNano(r)) << endl;
	}

	//誘電率計算
	double calcEPS(const double& x, const double& y, enum INTEG f);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY) {
		return 0;
	}
	bool update(int a){
		 // return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_incidenceModel :public FazzyModel {
	const double 
		ep1, ep2;
	const int alpha;
	double alphaR, length, ln, lx, ly, rn, cwidth, r;
public:
	FazzyHair_incidenceModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()


	{
	}
};

class FazzyHair_incidenceLayerModel :public FazzyModel {
	const double ep1, ep2;
	const int alpha;
	double alphaR, r, rn, cwidth, cn, cmc, mn, length, ln, lx, ly;
	/*
	alpha:傾き		alphaR:傾き(ラジアン)
	r:毛皮質範囲の半径(μm)							rn:毛皮質範囲の半径(nmシミュレーション値)
	cwidth:キューティクル厚さ(μm)					cn:キューティクル厚さ(nmシミュレーション値)
	cmc:CMC幅(μm)									mn:CMC幅(nmシミュレーション値)
	length:キューティクル長さ(μm)					ln:キューティクル長さ(nmシミュレーション値)
	ly:キューティクル範囲(nmシミュレーション値)		lx:x方向長さ(nmシミュレーション値)
	*/
public:
	FazzyHair_incidenceLayerModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_incidenceLayerModel_try :public FazzyModel {
	const double ep1, ep2;
	const int alpha;
	double alphaR, r, rn, cwidth, cn, cmc, mn, length, ln, lx, ly;
	/*
	alpha:傾き		alphaR:傾き(ラジアン)
	r:毛皮質範囲の半径(μm)							rn:毛皮質範囲の半径(nmシミュレーション値)
	cwidth:キューティクル厚さ(μm)					cn:キューティクル厚さ(nmシミュレーション値)
	cmc:CMC幅(μm)									mn:CMC幅(nmシミュレーション値)
	length:キューティクル長さ(μm)					ln:キューティクル長さ(nmシミュレーション値)
	ly:キューティクル範囲(nmシミュレーション値)		lx:x方向長さ(nmシミュレーション値)
	*/
public:
	FazzyHair_incidenceLayerModel_try(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_normalModel :public FazzyModel {
	const double ep1, ep2, e;
	const int r;
	double rn, ax, by;
public:
	FazzyHair_normalModel(Field*);
	string mkdir(string root);
	double calcEPS(const double&, const double&, enum INTEG f);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};

class FazzyHair_NONcuticleModel :public FazzyModel {
	const double ep1, ep2;
	const int r;
	double rn;
public:
	FazzyHair_NONcuticleModel(Field*);
	string mkdir(string root);
    double calcEPS(const double&, const double&, enum INTEG = D_XY);
	double calcSIG(const double&, const double&, const double lam, enum INTEG = D_XY);
	bool update(int) {
		//return true;
		return false;
	}
	void Initialize()
	{
	}
};


//the simulation parameter of BuprestidaeModel is (2038, 2038, 2, 20),the parameter of the field
class BuprestidaeModel :public FazzyModel {
	double epOfWhite, epOfBlack;
	double sigOfWhite, sigOfBlack;
	double realPartOfWhite[7];
	double imaginaryPartOfWhite[7];
	double realPartOfBlack[7];
	double imaginaryPartOfBlack[7];
	double lambda, w;
public:
	BuprestidaeModel(Field*, double);
	string mkdir(string root);
	double calcEPS(const double& x, const double& y, enum INTEG f = D_XY);
	double calcSIG(const double& x, const double& y, const double lam, enum INTEG f = D_XY);
	double calcEPSFromN(bool white);
	double calcSIGFromN(bool white);
	double initializeRoW(int i);
	double initializeRoB(int i);
	double initializeIoW(int i);
	double initializeIoB(int i);
	virtual void InitializeLambda(double a) {
        #define WHITE true
        #define BLACK false
		setLambda(a);
		epOfWhite = calcEPSFromN(WHITE);
		epOfBlack = calcEPSFromN(BLACK);
		sigOfWhite = calcSIGFromN(WHITE);
		sigOfBlack = calcSIGFromN(BLACK);
		cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
		cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
	}
	double getNByLambda(double a[]) {
		if (lambda >= 650)
			return a[5] + (lambda - 650)*(a[6] - a[5]) / 50;
		if (lambda >= 600)
			return a[4] + (lambda - 600)*(a[5] - a[4]) / 50;
		if (lambda >= 550)
			return a[3] + (lambda - 550)*(a[4] - a[3]) / 50;
		if (lambda >= 500)
			return a[2] + (lambda - 500)*(a[3] - a[2]) / 50;
		if (lambda >= 450)
			return a[1] + (lambda - 450)*(a[2] - a[1]) / 50;
		if (lambda >= 400)
			return a[0] + (lambda - 400)*(a[1] - a[0]) / 50;

		return 0.0;
	}
	void setLambda(double i) {
		this->lambda = i;
		w = LIGHT_SPEED_S * (2 * PI / lambda);
		
	}
	
	bool update(int) {
		return false;
	}

	void Initialize() {
	}
};

/*
simulation field (1700,1700,5,20)
*/

class BuprestidaeModelWithNoise :public FazzyModel {

	double epOfWhite, epOfBlack;
	double sigOfWhite, sigOfBlack;
	double realPartOfWhite[7];
	double imaginaryPartOfWhite[7];
	double realPartOfBlack[7];
	double imaginaryPartOfBlack[7];
	double lambda, w;
public:
	BuprestidaeModelWithNoise(Field*, double);
	string mkdir(string root);
	double calcEPS(const double& x, const double& y, enum INTEG f = D_XY);
	double calcSIG(const double& x, const double& y, const double lam, enum INTEG f = D_XY);
	double calcEPSFromN(bool white);
	double calcSIGFromN(bool white);
	double initializeRoW(int i);
	double initializeRoB(int i);
	double initializeIoW(int i);
	double initializeIoB(int i);
	double calcMultilayer(double x, double y,double start);
	double selectLayer(double x, double start);
	virtual void InitializeLambda(double a) {
#define WHITE true
#define BLACK false
		setLambda(a);
		epOfWhite = calcEPSFromN(WHITE);
		epOfBlack = calcEPSFromN(BLACK);
		sigOfWhite = calcSIGFromN(WHITE);
		sigOfBlack = calcSIGFromN(BLACK);
		cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
		cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
	}
	double getNByLambda(double a[]) {
		if (lambda >= 650)
			return a[5] + (lambda - 650)*(a[6] - a[5]) / 50;
		if (lambda >= 600)
			return a[4] + (lambda - 600)*(a[5] - a[4]) / 50;
		if (lambda >= 550)
			return a[3] + (lambda - 550)*(a[4] - a[3]) / 50;
		if (lambda >= 500)
			return a[2] + (lambda - 500)*(a[3] - a[2]) / 50;
		if (lambda >= 450)
			return a[1] + (lambda - 450)*(a[2] - a[1]) / 50;
		if (lambda >= 400)
			return a[0] + (lambda - 400)*(a[1] - a[0]) / 50;

		return 0.0;
	}
	void setLambda(double i) {
		this->lambda = i;
		w = LIGHT_SPEED_S * (2 * PI / lambda);

	}

	bool update(int) {
		return false;
	}

	void Initialize() {
	}

};

//1865,1865,5,20
class BuprestidaeModelsmooth2nd :public FazzyModel {
	double epOfWhite, epOfBlack;
	double sigOfWhite, sigOfBlack;
	double realPartOfWhite[7];
	double imaginaryPartOfWhite[7];
	double realPartOfBlack[7];
	double imaginaryPartOfBlack[7];
	double lambda, w;
public:
	BuprestidaeModelsmooth2nd(Field*, double);
	string mkdir(string root);
	double calcEPS(const double& x, const double& y, enum INTEG f = D_XY);
	double calcSIG(const double& x, const double& y, const double lam, enum INTEG f = D_XY);
	double calcEPSFromN(bool white);
	double calcSIGFromN(bool white);
	double initializeRoW(int i);
	double initializeRoB(int i);
	double initializeIoW(int i);
	double initializeIoB(int i);
	double calcMultilayer(double x, double y, double start);
	double selectLayer(double x, double start);
	virtual void InitializeLambda(double a) {
#define WHITE true
#define BLACK false
		setLambda(a);
		epOfWhite = calcEPSFromN(WHITE);
		epOfBlack = calcEPSFromN(BLACK);
		sigOfWhite = calcSIGFromN(WHITE);
		sigOfBlack = calcSIGFromN(BLACK);
		cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
		cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
	}
	double getNByLambda(double a[]) {
		if (lambda >= 650)
			return a[5] + (lambda - 650)*(a[6] - a[5]) / 50;
		if (lambda >= 600)
			return a[4] + (lambda - 600)*(a[5] - a[4]) / 50;
		if (lambda >= 550)
			return a[3] + (lambda - 550)*(a[4] - a[3]) / 50;
		if (lambda >= 500)
			return a[2] + (lambda - 500)*(a[3] - a[2]) / 50;
		if (lambda >= 450)
			return a[1] + (lambda - 450)*(a[2] - a[1]) / 50;
		if (lambda >= 400)
			return a[0] + (lambda - 400)*(a[1] - a[0]) / 50;

		return 0.0;
	}
	void setLambda(double i) {
		this->lambda = i;
		w = LIGHT_SPEED_S * (2 * PI / lambda);

	}

	bool update(int) {
		return false;
	}

	void Initialize() {
	}
};

//field(1845,1845,5,20)

class BuprestidaeModelWithNoise2nd :public FazzyModel {

	double epOfWhite, epOfBlack;
	double sigOfWhite, sigOfBlack;
	double realPartOfWhite[7];
	double imaginaryPartOfWhite[7];
	double realPartOfBlack[7];
	double imaginaryPartOfBlack[7];
	double lambda, w;
public:
	BuprestidaeModelWithNoise2nd(Field*, double);
	string mkdir(string root);
	double calcEPS(const double& x, const double& y, enum INTEG f = D_XY);
	double calcSIG(const double& x, const double& y, const double lam, enum INTEG f = D_XY);
	double calcEPSFromN(bool white);
	double calcSIGFromN(bool white);
	double initializeRoW(int i);
	double initializeRoB(int i);
	double initializeIoW(int i);
	double initializeIoB(int i);
	double calcMultilayer(double x, double y, double start);
	double selectLayer(double x, double start);
	virtual void InitializeLambda(double a) {
#define WHITE true
#define BLACK false
		setLambda(a);
		epOfWhite = calcEPSFromN(WHITE);
		epOfBlack = calcEPSFromN(BLACK);
		sigOfWhite = calcSIGFromN(WHITE);
		sigOfBlack = calcSIGFromN(BLACK);
		cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
		cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
	}
	double getNByLambda(double a[]) {
		if (lambda >= 650)
			return a[5] + (lambda - 650)*(a[6] - a[5]) / 50;
		if (lambda >= 600)
			return a[4] + (lambda - 600)*(a[5] - a[4]) / 50;
		if (lambda >= 550)
			return a[3] + (lambda - 550)*(a[4] - a[3]) / 50;
		if (lambda >= 500)
			return a[2] + (lambda - 500)*(a[3] - a[2]) / 50;
		if (lambda >= 450)
			return a[1] + (lambda - 450)*(a[2] - a[1]) / 50;
		if (lambda >= 400)
			return a[0] + (lambda - 400)*(a[1] - a[0]) / 50;

		return 0.0;
	}
	void setLambda(double i) {
		this->lambda = i;
		w = LIGHT_SPEED_S * (2 * PI / lambda);

	}

	bool update(int) {
		return false;
	}

	void Initialize() {
	}

};


//field(2345,2345,5)
class BuprestidaeModelSmooth24 :public FazzyModel {

	double epOfWhite, epOfBlack;
	double sigOfWhite, sigOfBlack;
	double realPartOfWhite[7];
	double imaginaryPartOfWhite[7];
	double realPartOfBlack[7];
	double imaginaryPartOfBlack[7];
	double lambda, w;
public:
	BuprestidaeModelSmooth24(Field*, double);
	string mkdir(string root);
	double calcEPS(const double& x, const double& y, enum INTEG f = D_XY);
	double calcSIG(const double& x, const double& y, const double lam, enum INTEG f = D_XY);
	double calcEPSFromN(bool white);
	double calcSIGFromN(bool white);
	double initializeRoW(int i);
	double initializeRoB(int i);
	double initializeIoW(int i);
	double initializeIoB(int i);
	double calcMultilayer(double x, double y, double start);
	double selectLayer(double x, double start);
	virtual void InitializeLambda(double a) {
#define WHITE true
#define BLACK false
		setLambda(a);
		epOfWhite = calcEPSFromN(WHITE);
		epOfBlack = calcEPSFromN(BLACK);
		sigOfWhite = calcSIGFromN(WHITE);
		sigOfBlack = calcSIGFromN(BLACK);
		cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
		cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
	}
	double getNByLambda(double a[]) {
		if (lambda >= 650)
			return a[5] + (lambda - 650)*(a[6] - a[5]) / 50;
		if (lambda >= 600)
			return a[4] + (lambda - 600)*(a[5] - a[4]) / 50;
		if (lambda >= 550)
			return a[3] + (lambda - 550)*(a[4] - a[3]) / 50;
		if (lambda >= 500)
			return a[2] + (lambda - 500)*(a[3] - a[2]) / 50;
		if (lambda >= 450)
			return a[1] + (lambda - 450)*(a[2] - a[1]) / 50;
		if (lambda >= 400)
			return a[0] + (lambda - 400)*(a[1] - a[0]) / 50;

		return 0.0;
	}
	void setLambda(double i) {
		this->lambda = i;
		w = LIGHT_SPEED_S * (2 * PI / lambda);

	}

	bool update(int) {
		return false;
	}

	void Initialize() {
	}

};



#endif //_MODEL_H