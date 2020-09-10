#ifndef _FDTD_H
#define _FDTD_H

#include "Object.h"
#include "Solver.h"
#include<complex>


class FDTD:public Solver{
	typedef complex<double> Complex;
protected:
	complex<double> *phi;	//ƒÓ
	double *np;
public:
	FDTD();
	virtual ~FDTD();
	virtual bool calc();
	virtual void draw();
	virtual void field() = 0;

protected:
	void Mie_Cylinder_Incidence();	//MiewU—ƒ‚ƒfƒ‹‚Ì“üË”g
	void Mie_Slub_Incidence();
	
	void Initialize(double);
	void Initialize();
};

class StFDTD:public FDTD{
protected:

public:
	StFDTD();
	~StFDTD();
	bool calc();
	void field();

protected:
	//‹zû‹«ŠE
	void absorbing();

	//üŠú‹«ŠE
	void cycle();

};

class NsFDTD: public FDTD{
public:
	NsFDTD();
	~NsFDTD();
	bool calc();
	void field();
private:
	double r_s, ky_s, kx_s;
	double DxDy2(double *p, int i, int j, int t);
	double D0_2(double *p, int i, int j, int t);
	complex<double> DxDy2(complex<double> *p, int i, int j, int t);
	complex<double>  D0_2(complex<double> *p, int i, int j, int t);

	//‹zû‹«ŠE
	virtual void absorbing_left();
	virtual void absorbing_right();
	virtual void absorbing_up();
	virtual void absorbing_down();
	virtual void absorbing();

	void Mie_Cylinder_Incidence();	//MiewU—ƒ‚ƒfƒ‹‚Ì“üË”g
};

#endif //_FDTD_H