#ifndef _ST_FDTD_TM_H
#define _ST_FDTD_TM_H
#include "FDTD_TM.h"

class StFDTD_TM: public FDTD_TM{
	typedef FDTD_TM super;
public:
	StFDTD_TM();
	~StFDTD_TM();
	bool calc();
	void field();
private:
	void absorbing();
	void cycle();
	bool EndTask();		//1��̃V�~�����[�V�������I������Ƃ��̏���
	void ReStart() {
		super::Initialize();
		field();
	}

	//todo ���E�ߖT��S-FDTD���g���ĂȂ�
	void CalcE() {
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i < mField->getNpx() - 1; i++)
				for (int j = 1; j < mField->getNpy() - 1; j++)
					EZX(i, j, +1) = CEZX(i, j)*EZX(i, j, 0) + CEZXLX(i, j)*(HY(i, j, 0) - HY(i - 1, j, 0));

		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i < mField->getNpx() - 1; i++)
				for (int j = 1; j < mField->getNpy() - 1; j++)
					EZY(i, j, +1) = CEZY(i, j)*EZY(i, j, 0) - CEZYLY(i, j)*(HX(i, j, 0) - HX(i, j - 1, 0));

		#ifdef _OPENMP
		#pragma omp for
		#endif
			for (int i = 1; i < mField->getNpx() - 1; i++)
				for (int j = 1; j < mField->getNpy() - 1; j++)
					EZ(i, j, +1) = EZX(i, j, +1) + EZY(i, j, +1);
		}
	}

	//�z�����E��Ez�ɂ����K�p���Ȃ�����,H�͗̈�̒[�����ʂɌv�Z����(�ł��镪��)
	void CalcH(){	//todo �v�Z�̈� i=0?, 1? j=0?, 1?
	#ifdef _OPENMP
	#pragma omp parallel
	#endif
		{
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=0; i<mField->getNpx(); i++)
			for(int j=0; j<mField->getNpy()-1; j++)
				HX(i,j,+1) = CHX(i,j)*HX(i,j,0) 
						- CHXLY(i,j)*( EZX(i,j+1,+1)-EZX(i,j,+1) + EZY(i,j+1,+1)-EZY(i,j,+1));	//Hx�̌v�Z Hx(i, j+1/2) -> Hx[i,j]

				
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i=0; i<mField->getNpx()-1; i++)
			for(int j=0; j<mField->getNpy(); j++)
				HY(i,j,+1) = CHY(i,j)*HY(i,j, 0) 
						+ CHYLX(i,j)*( EZX(i+1,j,+1)-EZX(i,j,+1) + EZY(i+1,j,+1)-EZY(i,j,+1) );	//Hy�̌v�Z Hy(i+1/2, j) -> Hy[i,j]

		}
	};

};
#endif //_ST_FDTD_TM_H