#include"Model.h"
#include"Field.h"

/*---------------------------------------------*/
/*--------------�~��Mie�U��--------------------*/
/*---------------------------------------------*/
FazzyMieModel::FazzyMieModel(Field *f):
FazzyModel(f)
{
	ep = 1.6*1.6*EPSILON_0_S;			//�U�d�� = (���ܗ�)^2
		
}

string FazzyMieModel::mkdir(string root){
	_mkdir((root + "Mie").c_str());

	string name = "Mie/" + to_s((int)(mField->cellToNano(r))) +"nm,"+ mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

double FazzyMieModel::calcEPS(const double& x, const double& y, enum INTEG f){

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)�����_�ɃV�t�g
	double _y = my - 0.5*mField->getNy();

	//���S�Ƃ̋��������a+��2/2�Z���ȏ�Ȃ�, ���S�ɔ}���̊O(�O�̂���, ���a+1 �ȏォ���ׂĂ���)
	if(_x*_x + _y*_y >= pow(r+1, 2.0))

		return EPSILON_0_S;

	//���S�Ƃ̋�����, ���a-��2/2�Z���ȓ��Ȃ�, ���S�ɔ}���̒�
	if(_x*_x + _y*_y <= pow(r-1, 2.0))
		return ep;

	double s=0;

	double a = 1.0,b=1.0;
	if(f == D_X) b = 0;
	if(f == D_Y) a = 0;
	for(double i=-16+0.5; i<16; i+=1)
		for(double j=-16+0.5; j<16; j+=1)
			if(pow(_x+a*i/32.0, 2.0) + pow(_y+b*j/32.0, 2.0) <= r*r)
				s+=1; 
	s /= 32.0*32.0;
	return ep*s + EPSILON_0_S*(1-s);
}

/*---------------------------------------------*/
/*--------------���w��-------------------------*/
/*---------------------------------------------*/
FazzySlabModel::FazzySlabModel(Field* f):
FazzyModel(f), ep1(2.0*2.0*EPSILON_0_S), ep2(EPSILON_0_S), width1(250), width2(50)
{
}

double FazzySlabModel::calcEPS(const double& x, const double& y, enum INTEG f){
//��100nm����,250nm�Ԋu��50nm�̃X���u�����Ă���  **��250nm����(L70.71)10nm�X���u�ɕύX(L73)
//���w��
	
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();

	if(mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy() ) return EPSILON_0_S;

	int k    = (int)(mField->cellToNano(mx) - 250)%250;
	double l =      (mField->cellToNano(mx) - 250)/250;

	if( k > 0 && k <=10 && l < 5)
		return ep1;
	else
		return ep2;

}

string FazzySlabModel::mkdir(string root){
	_mkdir((root + "SlabModel").c_str());

	string name = "SlabModel/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

/*---------------------------------------------*/
/*---------------�є�--------------------------*/
/*---------------------------------------------*/

/*----------------------------*/
/*-----------�c�f��-----------*/
/*----------------------------*/
FazzyHair_incidenceModel::FazzyHair_incidenceModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(5), cwidth(0.5), r(32+(4.36-cwidth))
	//alpha:�L���[�e�B�N���̊p�x(deg)  cwidth:�L���[�e�B�N���̕�(��m)  r:�т̔��a(��m)(���a+�L���[�e�B�N�����d�Ȃ�̈�)
{
	alphaR = alpha * PI / 180;
	length = cwidth / sin(alphaR);
	cout << "�L���[�e�B�N���̊p�x : " + to_s(alpha) + "deg" << endl;
	cout << "�L���[�e�B�N���� : " + to_s(cwidth) + "micro" << endl;
	cout << "�L���[�e�B�N��1���̘I�o�� : " + to_s(length) + "micro" << endl;
}
/*
	alpha:�X��		alphaR:�X��(���W�A��)
	r:�є玿�͈͂̔��a(��m)							rn:�є玿�͈͂̔��a(nm�V�~�����[�V�����l)
	cwidth:�L���[�e�B�N������(��m)					cn:�L���[�e�B�N������(nm�V�~�����[�V�����l)
	cmc:CMC��(��m)									mn:CMC��(nm�V�~�����[�V�����l)
	length:�L���[�e�B�N������(��m)					ln:�L���[�e�B�N������(nm�V�~�����[�V�����l)
	ly:�L���[�e�B�N���͈�(nm�V�~�����[�V�����l)		lx:x��������(nm�V�~�����[�V�����l)
	*/

double FazzyHair_incidenceModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	rn = mField->nanoToCell(r * 1000);
	
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;
	
	double h = mField->nanoToCell(0*1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�


	/***************************************************/
	//��ʂ����L���[�e�B�N���Ȃ�
//	if (y - mField->getNpml() >= mField->nanoToCell(32 * 1000) + cy)	return ep2;
	/***************************************************/


	int c = mField->getNx() / lx + 1;		//�v�Z�͈͓��̃L���[�e�B�N���̐�
	for (int i = 0; i < c; i++) {
		if (mx > i * lx + h && mx < (i + 1) * lx + h && mx < mField->getNx() - h) {
			//			if (my > tan(alphaR) * (mx - lx*i) + cy + rn)	return ep2;
			//			else return ep1;		//Fuzzy�Ȃ�(Staircase���f��)

			double dy1 = my - (tan(alphaR) * (mx - lx*i - h) + cy + rn);
			double dy2 = my - (tan(alphaR) * ((mx - lx*i - h) + 1) + cy + rn);
			double s;
			if (dy1 > 0 && dy2 > 0) return ep2;		//�L���[�e�B�N�������̊O�� (1)
			if (fabs(dy1) > 1 && fabs(dy2) > 1) return ep1;		//�L���[�e�B�N�������̓��� (2)

			if (dy1 <= 0 && dy2 <= 0) {
				if (fabs(dy1) <= 1 && fabs(dy2) <= 1) {
					s = (fabs(dy1) + fabs(dy2)) * 1.0 / 2.0;
					return ep1 * s + ep2 * (1 - s);		// (3)
				}
				if (fabs(dy1) < 1 && fabs(dy2) > 1) {
					s = (1 - fabs(dy1)) * ((my - cy - rn) / tan(alphaR) - (mx - lx*i - h)) / 2;
					return ep2 * s + ep1 * (1 - s);		// (4)
				}
			}
			if (dy1 > 0 && dy2 < 0) {
				s = fabs(dy2) * (((mx - lx*i - h) + 1) - (my - cy - rn) / tan(alphaR)) / 2;
				return ep1 * s + ep2 * (1 - s);		// (5)
			}
		}
		else
			continue;
			//break;
	}


	return ep2;
}

double FazzyHair_incidenceModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//���F�̐F�f
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceplane").c_str());				//�z���W���Ȃ��̏ꍇ
		name = "HairModel/incidenceplane/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceplane_withSig").c_str());		//�z���W������̏ꍇ
		name = "HairModel/incidenceplane_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	
	return name + "/";
}

/*--------------------------------------------------*/
/*-----------�c�f��(���w���L���[�e�B�N��)-----------*/
/*--------------------------------------------------*/
FazzyHair_incidenceLayerModel::FazzyHair_incidenceLayerModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(5), cwidth(0.5), length(50), r(32), cmc(0.06)
	//alpha:�L���[�e�B�N���̊p�x(deg)  cwidth:�L���[�e�B�N���̌���(��m)  length:�L���[�e�B�N���̒���(��m)	r:�є玿�͈͂̔��a(��m)  cmc:CMC�͈�(��m)
{
	alphaR = alpha * PI / 180;
	int n = length * sin(alphaR) / (cmc + cwidth);

	cout << "�L���[�e�B�N���̊p�x : " + to_s(alpha) + "deg" << endl;
	cout << "�L���[�e�B�N������ : " + to_s(cwidth) + "micro" << endl;
	cout << "�L���[�e�B�N������ : " + to_s(length) + "micro" << endl;
	cout << "�L���[�e�B�N��1���̘I�o�� : " + to_s((cmc+cwidth)/sin(alphaR)) + "micro" << endl;
	cout << "�L���[�e�B�N���͈͕� : " + to_s(length*sin(alphaR)) + "micro" << endl;
	cout << "�L���[�e�B�N���d�Ȃ薇�� : " + to_s(n) + "��" << endl;
}

double FazzyHair_incidenceLayerModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 1000);
	cn = mField->nanoToCell(cwidth * 1000);

	mn = mField->nanoToCell(cmc * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w

	/**** �є��S�� ****/
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�
	if (my <= rn + cy)		return ep1;		//�є玿����

	my = my - rn - cy;		//���ړ�

	/****************** �L���[�e�B�N�������̂� ******************/
	/* Field�����T�C�Y                                          */
	/* Field(8000, 8000, 5, 10) Field(16000, 8000, 10, 10) �Ȃ� */
	/************************************************************/
//	if (my <= 100)		return ep1;		//�є玿����
//	my = my - 100;		//���ړ�


	if (my > ly)	return ep2;
	for (int i = -(ly / (cn + mn)); (i*mn + i*cn) / tan(alphaR) < mField->getNx(); i++) {
		if (my <= ly - cn) {
			if (my >= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + i*cn) / tan(alphaR)) - 1)
				return ep1;

			else if (my >= tan(alphaR) * (mx - ((i + 1)*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1)
				return ep2;

			else if ((my > tan(alphaR) * (mx - (i*mn + i*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + i*cn) / tan(alphaR)) + 1)
				|| (my > tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1)) {
				
				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (my + b / 32.0 >= tan(alphaR) * (mx + a / 32.0 - (i*mn + (i + 1)*cn) / tan(alphaR)) && my + b / 32.0 <= tan(alphaR) * (mx + a / 32.0 - (i*mn + i*cn) / tan(alphaR)) && mx <= (ly + i*mn + i*cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1*s + ep2*(1 - s);
				
			}
			
		}
		else if (my > ly - cn) {
			if (mx >= (my + i*mn + i*cn + 1) / tan(alphaR) && mx <= (ly + i*mn + i*cn) / tan(alphaR))
				return ep1;

			else if (mx < (my + i*mn + i*cn + 1) / tan(alphaR) && mx >(my + i*mn + i*cn - 1) / tan(alphaR)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (mx + a / 32.0 >= (my + b / 32.0 + i*mn + i*cn) / tan(alphaR) && mx <= (ly + i*mn + i*cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1*s + ep2*(1 - s);
			}	
		}
	}

	return ep2;
}

double FazzyHair_incidenceLayerModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//���F�̐F�f
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceLayerModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceLayer").c_str());				//�z���W���Ȃ��̏ꍇ
		name = "HairModel/incidenceLayer/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceLayer_withSig").c_str());		//�z���W������̏ꍇ
		name = "HairModel/incidenceLayer_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬

	return name + "/";
}


FazzyHair_incidenceLayerModel_try::FazzyHair_incidenceLayerModel_try(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), alpha(5), cwidth(0.5), length(30), r(4), cmc(0.06)
	//alpha:�L���[�e�B�N���̊p�x(deg)  cwidth:�L���[�e�B�N���̌���(��m)  length:�L���[�e�B�N���̒���(��m)	r:�є玿�͈͂̔��a(��m)  cmc:CMC�͈�(��m)
{
	alphaR = alpha * PI / 180;
	int n = length * sin(alphaR) / (cmc + cwidth);

	cout << "�L���[�e�B�N���̊p�x : " + to_s(alpha) + "deg" << endl;
	cout << "�L���[�e�B�N������ : " + to_s(cwidth) + "micro" << endl;
	cout << "�L���[�e�B�N������ : " + to_s(length) + "micro" << endl;
	cout << "�L���[�e�B�N��1���̘I�o�� : " + to_s((cmc + cwidth) / sin(alphaR)) + "micro" << endl;
	cout << "�L���[�e�B�N���͈͕� : " + to_s(length*sin(alphaR)) + "micro" << endl;
	cout << "�L���[�e�B�N���d�Ȃ薇�� : " + to_s(n) + "��" << endl;
}

double FazzyHair_incidenceLayerModel_try::calcEPS(const double& x, const double& y, enum INTEG f) {
	alphaR = alpha * PI / 180;
	ln = mField->nanoToCell(length * 1000);
	lx = ln * cos(alphaR);
	ly = ln * sin(alphaR);
	rn = mField->nanoToCell(r * 1000);
	cn = mField->nanoToCell(cwidth * 1000);

	mn = mField->nanoToCell(cmc * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w

	/**** �є��S�� ****/
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�
	if (my <= rn + cy)		return ep1;		//�є玿����

	my = my - rn - cy;		//���ړ�

	/****************** �L���[�e�B�N�������̂� ******************/
	/* Field�����T�C�Y                                          */
	/* Field(8000, 8000, 5, 10) Field(16000, 8000, 10, 10) �Ȃ� */
	/************************************************************/
//	if (my <= 100)		return ep1;		//�є玿����
//	my = my - 100;		//���ړ�


	if (my > ly)	return ep2;
	for (int i = -(ly / (cn + mn)); (i*mn + i * cn) / tan(alphaR) < mField->getNx(); i++) {
		if (my <= ly - cn) {
			if (my >= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) - 1)
				return ep1;

			else if (my >= tan(alphaR) * (mx - ((i + 1)*mn + (i + 1)*cn) / tan(alphaR)) + 1 && my <= tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1)
				return ep2;

			else if ((my > tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + i * cn) / tan(alphaR)) + 1)
				|| (my > tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) - 1 && my < tan(alphaR) * (mx - (i*mn + (i + 1)*cn) / tan(alphaR)) + 1)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (my + b / 32.0 >= tan(alphaR) * (mx + a / 32.0 - (i*mn + (i + 1)*cn) / tan(alphaR)) && my + b / 32.0 <= tan(alphaR) * (mx + a / 32.0 - (i*mn + i * cn) / tan(alphaR)) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1 * s + ep2 * (1 - s);

			}

		}
		else if (my > ly - cn) {
			if (mx >= (my + i * mn + i * cn + 1) / tan(alphaR) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
				return ep1;

			else if (mx < (my + i * mn + i * cn + 1) / tan(alphaR) && mx >(my + i * mn + i * cn - 1) / tan(alphaR)) {

				double s = 0;
				for (double a = -16; a < 16; a += 1)
					for (double b = -16; b < 16; b += 1)
						if (mx + a / 32.0 >= (my + b / 32.0 + i * mn + i * cn) / tan(alphaR) && mx <= (ly + i * mn + i * cn) / tan(alphaR))
							s += 1;
				s /= 32.0*32.0;
				return ep1 * s + ep2 * (1 - s);
			}
		}
	}

	return ep2;
}

double FazzyHair_incidenceLayerModel_try::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(32 * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy) {
		
		//int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		//double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;
		
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) ) % 8000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) ) / 8000;
		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//���F�̐F�f
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_incidenceLayerModel_try::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/incidenceLayer").c_str());				//�z���W���Ȃ��̏ꍇ
		name = "HairModel/incidenceLayer/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/incidenceLayer_withSig").c_str());		//�z���W������̏ꍇ
		name = "HairModel/incidenceLayer_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬

	return name + "/";
}

/*------------------------------------------------*/
/*-----------�c�f��(�L���[�e�B�N���Ȃ�)-----------*/
/*------------------------------------------------*/
FazzyHair_NONcuticleModel::FazzyHair_NONcuticleModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), r(32)
	//r:�т̔��a(��m)
{
	cout << "�L���[�e�B�N�� : �Ȃ�"<< endl;
}

double FazzyHair_NONcuticleModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy)	return ep1;
	else  return ep2;
}

double FazzyHair_NONcuticleModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	double h = mField->nanoToCell(0 * 1000);

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;	//PML�w
	if (my < cy)	my = 2 * cy - my;		//x���ɑ΂��Đ��Ώ�

	if (my <= rn + cy) {
		int k = (int)(mField->cellToNano(my) - mField->cellToNano(cy) - 1500) % 4000;
		double l = (mField->cellToNano(my) - mField->cellToNano(cy) - 1500) / 4000;

		if (k > 0 && k <= 1000 && l < 8)	return 1.0;		//���F�̐F�f
		else	return 0;
	}
	else  return 0;
}

string FazzyHair_NONcuticleModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	string name;

	if (mField->sig == false) {
		_mkdir((root + "HairModel/NONcuticle").c_str());				//�z���W���Ȃ��̏ꍇ
		name = "HairModel/NONcuticle/" + mField->getStringCellInfo();
	}
	else if (mField->sig == true) {
		_mkdir((root + "HairModel/NONcuticle_withSig").c_str());		//�z���W������̏ꍇ
		name = "HairModel/NONcuticle_withSig/" + mField->getStringCellInfo();
	}

	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

/*----------------------------*/
/*-----------���f��-----------*/
/*----------------------------*/
FazzyHair_normalModel::FazzyHair_normalModel(Field* f) :
	FazzyModel(f), ep1(1.55*1.55*EPSILON_0_S), ep2(EPSILON_0_S), e(0.6), r(32)
	//a:���S��  r:�т̔��a(��m)
{
	cout << "�ȉ~�̗��S�� = " + to_s((double)e) << endl;
}

double FazzyHair_normalModel::calcEPS(const double& x, const double& y, enum INTEG f) {
	rn = mField->nanoToCell(r * 1000);
	ax = rn;
	by = ax * sqrt(1 - e*e);

	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	double cx = mField->getNx() / 2;
	double cy = mField->getNy() / 2;

	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return ep2;	//PML�w

	double _x = mx - 0.5*mField->getNx();//(N_X/2, N_Y/2)�����_�ɃV�t�g
	double _y = my - 0.5*mField->getNy();

	//���S�Ƃ̋��������a+��2/2�Z���ȏ�Ȃ�, ���S�ɔ}���̊O(�O�̂���, ���a+1 �ȏォ���ׂĂ���)
	double _ax = ax+1, _by = by+1;
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) >= 1)
		return ep2;

	_ax = ax - 1;
	_by = by - 1;
	//���S�Ƃ̋�����, ���a-��2/2�Z���ȓ��Ȃ�, ���S�ɔ}���̊O
	if ((_x*_x) / (_ax*_ax) + (_y*_y) / (_by*_by) <= 1)
		return ep1;

	double s = 0;

	double a = 1.0, b = 1.0;
	if (f == D_X) b = 0;
	if (f == D_Y) a = 0;
	for (double i = -16 + 0.5; i<16; i += 1)
		for (double j = -16 + 0.5; j<16; j += 1)
			if (pow(_x + a*i / 32.0, 2.0) / (ax*ax) + pow(_y + b*j / 32.0, 2.0) / (by*by) <= 1)
				s += 1;
	s /= 32.0*32.0;
	return ep1*s + ep2*(1 - s);
}

string FazzyHair_normalModel::mkdir(string root) {
	_mkdir((root + "HairModel").c_str());
	_mkdir((root + "HairModel/normalplane").c_str());
	
	string name = "HairModel/normalplane/e=" + to_s((double)e);
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	
	name = "HairModel/normalplane/e=" + to_s((double)e) + "/" + mField->getStringCellInfo();
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "/";
}

BuprestidaeModel::BuprestidaeModel(Field* f,double lambda):FazzyModel(f){
	#define WHITE true
    #define BLACK false
	//this->lambda = lambda;
	//this->w = w;
	setLambda(lambda);
	
	for (int i = 0; i < 7; i++) {
		realPartOfWhite[i] = initializeRoW(i);
		realPartOfBlack[i] = initializeRoB(i);
		imaginaryPartOfWhite[i] = initializeIoW(i);
		imaginaryPartOfBlack[i] = initializeIoB(i);
	}

	epOfWhite = calcEPSFromN(WHITE);
	epOfBlack = calcEPSFromN(BLACK);
	sigOfWhite = calcSIGFromN(WHITE);
	sigOfBlack = calcSIGFromN(BLACK);
	cout <<"epw="<< epOfWhite << " epb=" << epOfBlack << endl;
	cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
}
double BuprestidaeModel::initializeRoW(int i) {
	switch (i) {
	case 0:
		return 1.606;
	case 1:
		return 1.588;
	case 2:
		return 1.573;
	case 3:
		return 1.562;
	case 4:
		return 1.550;
	case 5:
		return 1.542;
	case 6:
		return 1.531;

	}
	return 0.0;
}

double BuprestidaeModel::initializeRoB(int i) {
	switch (i) {
	case 0:
		return 1.784;
	case 1:
		return 1.745;
	case 2:
		return 1.718;
	case 3:
		return 1.683;
	case 4:
		return 1.669;
	case 5:
		return 1.658;
	case 6:
		return 1.640;

	}
	return 0.0;
}

double BuprestidaeModel::initializeIoW(int i) {
	switch (i) {
	case 0:
		return 0.017;
	case 1:
		return 0.003;
	default:
		break;
	}
	return 0.0;
}

double BuprestidaeModel::initializeIoB(int i) {
	switch (i) {
	case 0:
		return 0.096;
	case 1:
		return 0.069;
	case 2:
		return 0.050;
	case 3:
		return 0.034;
	case 4:
		return 0.025;
	case 5:
		return 0.021;
	case 6:
		return 0.013;

	}
	return 0.0;
}

string BuprestidaeModel::mkdir(string root) {
	_mkdir((root + "Buprestidae").c_str());

	string name = "Buprestidae\\"+to_s(int(lambda));
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name+"\\" ;
}

double BuprestidaeModel::calcEPS(const double& x,const double& y,enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0&&mx - start_x <length*10) {
		int calc_x = (int)(mx - start_x) % length;
		if (calc_x < w_length)
			return epOfWhite;
		return epOfBlack;
		
	}
	
	return EPSILON_0_S;
}

double BuprestidaeModel::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;

	int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0&&mx - start_x<length*10) {
		int calc_x = (int)(mx - start_x) % length;
		if (calc_x < w_length)
			return sigOfWhite;
		return sigOfBlack;
		
	}
	

	return 0;
}

double BuprestidaeModel::calcEPSFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return realP * realP - imaginaryP * imaginaryP;
	
}

double BuprestidaeModel::calcSIGFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return 2 * realP*imaginaryP*w;
	
}

BuprestidaeModelWithNoise::BuprestidaeModelWithNoise(Field* f, double lambda) :FazzyModel(f) {
#define WHITE true
#define BLACK false
	//this->lambda = lambda;
	//this->w = w;
	setLambda(lambda);

	for (int i = 0; i < 7; i++) {
		realPartOfWhite[i] = initializeRoW(i);
		realPartOfBlack[i] = initializeRoB(i);
		imaginaryPartOfWhite[i] = initializeIoW(i);
		imaginaryPartOfBlack[i] = initializeIoB(i);
	}

	epOfWhite = calcEPSFromN(WHITE);
	epOfBlack = calcEPSFromN(BLACK);
	sigOfWhite = calcSIGFromN(WHITE);
	sigOfBlack = calcSIGFromN(BLACK);
	cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
	cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
}
double BuprestidaeModelWithNoise::initializeRoW(int i) {
	switch (i) {
	case 0:
		return 1.606;
	case 1:
		return 1.588;
	case 2:
		return 1.573;
	case 3:
		return 1.562;
	case 4:
		return 1.550;
	case 5:
		return 1.542;
	case 6:
		return 1.531;

	}
	return 0.0;
}

double BuprestidaeModelWithNoise::initializeRoB(int i) {
	switch (i) {
	case 0:
		return 1.784;
	case 1:
		return 1.745;
	case 2:
		return 1.718;
	case 3:
		return 1.683;
	case 4:
		return 1.669;
	case 5:
		return 1.658;
	case 6:
		return 1.640;

	}
	return 0.0;
}

double BuprestidaeModelWithNoise::initializeIoW(int i) {
	switch (i) {
	case 0:
		return 0.017;
	case 1:
		return 0.003;
	default:
		break;
	}
	return 0.0;
}

double BuprestidaeModelWithNoise::initializeIoB(int i) {
	switch (i) {
	case 0:
		return 0.096;
	case 1:
		return 0.069;
	case 2:
		return 0.050;
	case 3:
		return 0.034;
	case 4:
		return 0.025;
	case 5:
		return 0.021;
	case 6:
		return 0.013;

	}
	return 0.0;
}

string BuprestidaeModelWithNoise::mkdir(string root) {
	_mkdir((root + "BuprestidaeWithNoise").c_str());

	string name = "BuprestidaeWithNoise\\" + to_s(int(lambda));
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "\\";
}

double BuprestidaeModelWithNoise::calcEPS(const double& x, const double& y, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x <length * 10) {
		int calc_x = (int)(mx - start_x) % length;
		if (calc_x < w_length)
			return epOfWhite;
		return epOfBlack;

	}*/

	double startLine = 60;
	double weight_w = calcMultilayer(mx, my, startLine);
	if(weight_w<0)
       return EPSILON_0_S;
	return weight_w * epOfWhite + (1 - weight_w)*epOfBlack;
}

double BuprestidaeModelWithNoise::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x<length * 10) {
		int calc_x = (int)(mx - start_x) % length;
		if (calc_x < w_length)
			return sigOfWhite;
		return sigOfBlack;

	}*/
	double startLine = 60;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w<0)
		return 0;
	return weight_w * sigOfWhite + (1 - weight_w)*sigOfBlack;

	
}

double BuprestidaeModelWithNoise::calcEPSFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return realP * realP - imaginaryP * imaginaryP;

}

double BuprestidaeModelWithNoise::calcSIGFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return 2 * realP*imaginaryP*w;

}


/*���G�ȑ��w��

���v18�w

in filed()
���f��
1�F���@150nm     30
1.5  5nm �@�@ (3��,2.05 ��)
2�F���@30nm      6
2.5  5nm�@    (2.55��,2.05��)�@
3�F���@115nm     23
3.5  5nm      (2.55��,0.55��)
4�F���@55nm      11
4.5  5nm      (1.75��,0.55��)
5�F���@95nm      19
5.5  5nm      (1.75��,2.413��)
6�F���@50nm      10
6.5  5nm      (1��,2.413��)
7�F���@85nm      17
7.5  5nm      (1��,1.85��)
8�F���@40nm      8
8.5  5nm      (0��,1.85��)
9�F���@90nm      18
9.5  5nm      (0��,1��)
10�F���@45nm     9
10.5 5nm      (1.225��,1��)

11�F���@80nm
11.5 5nm      (1.225��,1.86��)
12�F��  55nm
12.5 5nm      (1.225��,1.86��)



*/
double BuprestidaeModelWithNoise::calcMultilayer(double x, double y,double start) {
	double layer = selectLayer(x, start);
	if (layer < 0)
		return -1;
	int layer_s = layer * 2 - 1;
	int a = layer_s % 4;
	if (a == 3)
		return 0;
	if (a == 1)
		return 1;

	switch (layer_s) {
	case 2:
		return 3 / (3 + 2.05);
	case 4:
		return 2.55 / (2.55 + 2.05);
	case 6:
		return 2.55 / (2.55 + 0.55);
	case 8:
		return 1.75 / (1.75 + 0.55);
	case 10:
		return 1.75 / (1.75 + 2.413);
	case 12:
		return 1.0 / (1.0 + 2.413);
	case 14:
		return 1.0 / (1.0 + 1.85);
	case 16:
		return 0 / (0 + 1.85);
	case 18:
		return 0 / (0 + 1.0);
	case 20:
		return 1.225 / (1.225 + 1.0);
	case 22:
		return 1.225 / (1.225 + 1.86);
	case 24:
		return 1.225 / (1.225 + 1.86);
	case 26:
		return 1.225 / (1.225 + 1.86);
	case 28:
		return 1.225 / (1.225 + 1.86);
	case 30:
		return 1.225 / (1.225 + 1.86);
	case 32:
		return 1.225 / (1.225 + 1.86);
	case 34:
		return 1.225 / (1.225 + 1.86);
	default:
		return 1.225 / (1.225 + 1.86);
	

	}
}



double BuprestidaeModelWithNoise::selectLayer(double x, double start) {
	double length = x - start;
	if (length < 0)
		return -1;
	if (length > 275)//-11 18
		return -1;
	if (length > 264)//-1
		return 18;
	if (length > 263)//-16 17
		return 17.5;
	if (length > 247)//-1
		return 17;
	if (length > 246)//-11 16
		return 16.5;
	if (length > 235)//-1
		return 16;
	if (length > 234)//-16 15
		return 15.5;
	if (length > 218)//-1
		return 15;
	if (length > 217)//-11 14
		return 14.5;
	if (length > 206)//-1
		return 14;
	if (length > 205)//-16 13
		return 13.5;
	if (length > 189)//-1
		return 13;
	if (length > 188)//-11 12
		return 12.5;
	if (length > 177)//-1
		return 12;
	if (length > 176)//-16 11
		return 11.5;
	if (length > 160)//-1
		return 11;
	if (length > 159)//-9  10
		return 10.5;
	if (length > 150)//-1
		return 10;
	if (length > 149)//-18 9
		return 9.5;
	if (length > 131)//-1
		return 9;
	if (length > 130)//-8  8
		return 8.5;
	if (length > 122)//-1
		return 8;
	if (length > 121)//-17  7
		return 7.5;
	if (length > 104)//-1
		return 7;
	if (length > 103)//-10  6
		return 6.5;
	if (length > 93)//-1
		return 6;
	if (length > 92)//-19  5
		return 5.5;
	if (length > 73)//-1
		return 5;
	if (length > 72)//-11  4
		return 4.5;
	if (length > 61)//-1
		return 4;
	if (length > 60)//-23  3
		return 3.5;
	if (length > 37)//-1
		return 3;
	if (length > 36)//-6 2
		return 2.5;
	if (length > 30)//-1
		return 2;
	if (length > 29)//-30
		return 1.5;
	return 1;
}

BuprestidaeModelsmooth2nd::BuprestidaeModelsmooth2nd(Field* f, double lambda) :FazzyModel(f) {
#define WHITE true
#define BLACK false
	//this->lambda = lambda;
	//this->w = w;
	setLambda(lambda);

	for (int i = 0; i < 7; i++) {
		realPartOfWhite[i] = initializeRoW(i);
		realPartOfBlack[i] = initializeRoB(i);
		imaginaryPartOfWhite[i] = initializeIoW(i);
		imaginaryPartOfBlack[i] = initializeIoB(i);
	}

	epOfWhite = calcEPSFromN(WHITE);
	epOfBlack = calcEPSFromN(BLACK);
	sigOfWhite = calcSIGFromN(WHITE);
	sigOfBlack = calcSIGFromN(BLACK);
	cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
	cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
}
double BuprestidaeModelsmooth2nd::initializeRoW(int i) {
	switch (i) {
	case 0:
		return 1.606;
	case 1:
		return 1.588;
	case 2:
		return 1.573;
	case 3:
		return 1.562;
	case 4:
		return 1.550;
	case 5:
		return 1.542;
	case 6:
		return 1.531;

	}
	return 0.0;
}

double BuprestidaeModelsmooth2nd::initializeRoB(int i) {
	switch (i) {
	case 0:
		return 1.784;
	case 1:
		return 1.745;
	case 2:
		return 1.718;
	case 3:
		return 1.683;
	case 4:
		return 1.669;
	case 5:
		return 1.658;
	case 6:
		return 1.640;

	}
	return 0.0;
}

double BuprestidaeModelsmooth2nd::initializeIoW(int i) {
	switch (i) {
	case 0:
		return 0.017;
	case 1:
		return 0.003;
	default:
		break;
	}
	return 0.0;
}

double BuprestidaeModelsmooth2nd::initializeIoB(int i) {
	switch (i) {
	case 0:
		return 0.096;
	case 1:
		return 0.069;
	case 2:
		return 0.050;
	case 3:
		return 0.034;
	case 4:
		return 0.025;
	case 5:
		return 0.021;
	case 6:
		return 0.013;

	}
	return 0.0;
}

string BuprestidaeModelsmooth2nd::mkdir(string root) {
	_mkdir((root + "BuprestidaeModelsmooth2nd").c_str());

	string name = "BuprestidaeModelsmooth2nd\\" + to_s(int(lambda));
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "\\";
}

double BuprestidaeModelsmooth2nd::calcEPS(const double& x, const double& y, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x <length * 10) {
		int calc_x = (int)(mx - start_x) % length;
		if (calc_x < w_length)
			return epOfWhite;
		return epOfBlack;

	}*/

	double startLine = 62;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w < 0)
		return EPSILON_0_S;
	return weight_w * epOfWhite + (1 - weight_w)*epOfBlack;
}

double BuprestidaeModelsmooth2nd::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x<length * 10) {
		int calc_x = (int)(mx - start_x) % length;
		if (calc_x < w_length)
			return sigOfWhite;
		return sigOfBlack;

	}*/
	double startLine = 62;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w < 0)
		return 0;
	return weight_w * sigOfWhite + (1 - weight_w)*sigOfBlack;

}

double BuprestidaeModelsmooth2nd::calcEPSFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return realP * realP - imaginaryP * imaginaryP;

}

double BuprestidaeModelsmooth2nd::calcSIGFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return 2 * realP*imaginaryP*w;

}

double BuprestidaeModelsmooth2nd::calcMultilayer(double x, double y, double start) {
	double layer = selectLayer(x, start);
	if (layer < 0)
		return -1;
	int layer_s = int(layer) % 2;
	
	if (layer_s == 1)
		return 1;
	if (layer_s == 0)
		return 0;

	
}

double BuprestidaeModelsmooth2nd::selectLayer(double x, double start) {
	double length = x - start;
	if (length < 0)
		return -1;
	if (length > 311)//-13
		return -1;
	if (length > 298)//-17
		return 20;
	if (length > 281)//-13
		return 19;
	if (length > 268)//-17
		return 18;
	if (length > 251)//-13
		return 17;
	if (length > 238)//-17
		return 16;
	if (length > 221)//-13
		return 15;
	if (length > 208)//-17
		return 14;
	if (length > 191)//-13
		return 13;
	if (length > 178)//-17
		return 12;
	if (length > 161)//-11
		return 11;
	if (length > 150)//-18
		return 10;
	if (length > 132)//-9
		return 9;
	if (length > 123)//-18
		return 8;
	if (length > 105)//-12
		return 7;
	if (length > 93)//-19
		return 6;
	if (length > 74)//-12
		return 5;
	if (length > 62)//-22
		return 4;
	if (length > 40)//-9
		return 3;
	if (length > 31)//-31
		return 2;
	 
	return 1;
}

BuprestidaeModelWithNoise2nd::BuprestidaeModelWithNoise2nd(Field* f, double lambda) :FazzyModel(f) {
#define WHITE true
#define BLACK false
	//this->lambda = lambda;
	//this->w = w;
	setLambda(lambda);

	for (int i = 0; i < 7; i++) {
		realPartOfWhite[i] = initializeRoW(i);
		realPartOfBlack[i] = initializeRoB(i);
		imaginaryPartOfWhite[i] = initializeIoW(i);
		imaginaryPartOfBlack[i] = initializeIoB(i);
	}

	epOfWhite = calcEPSFromN(WHITE);
	epOfBlack = calcEPSFromN(BLACK);
	sigOfWhite = calcSIGFromN(WHITE);
	sigOfBlack = calcSIGFromN(BLACK);
	cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
	cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
}
double BuprestidaeModelWithNoise2nd::initializeRoW(int i) {
	switch (i) {
	case 0:
		return 1.606;
	case 1:
		return 1.588;
	case 2:
		return 1.573;
	case 3:
		return 1.562;
	case 4:
		return 1.550;
	case 5:
		return 1.542;
	case 6:
		return 1.531;

	}
	return 0.0;
}

double BuprestidaeModelWithNoise2nd::initializeRoB(int i) {
	switch (i) {
	case 0:
		return 1.784;
	case 1:
		return 1.745;
	case 2:
		return 1.718;
	case 3:
		return 1.683;
	case 4:
		return 1.669;
	case 5:
		return 1.658;
	case 6:
		return 1.640;

	}
	return 0.0;
}

double BuprestidaeModelWithNoise2nd::initializeIoW(int i) {
	switch (i) {
	case 0:
		return 0.017;
	case 1:
		return 0.003;
	default:
		break;
	}
	return 0.0;
}

double BuprestidaeModelWithNoise2nd::initializeIoB(int i) {
	switch (i) {
	case 0:
		return 0.096;
	case 1:
		return 0.069;
	case 2:
		return 0.050;
	case 3:
		return 0.034;
	case 4:
		return 0.025;
	case 5:
		return 0.021;
	case 6:
		return 0.013;

	}
	return 0.0;
}

string BuprestidaeModelWithNoise2nd::mkdir(string root) {
	_mkdir((root + "BuprestidaeWithNoise2nd").c_str());

	string name = "BuprestidaeWithNoise2nd\\" + to_s(int(lambda));
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "\\";
}

double BuprestidaeModelWithNoise2nd::calcEPS(const double& x, const double& y, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x <length * 10) {
	int calc_x = (int)(mx - start_x) % length;
	if (calc_x < w_length)
	return epOfWhite;
	return epOfBlack;

	}*/

	double startLine = 60;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w<0)
		return EPSILON_0_S;
	return weight_w * epOfWhite + (1 - weight_w)*epOfBlack;
}

double BuprestidaeModelWithNoise2nd::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x<length * 10) {
	int calc_x = (int)(mx - start_x) % length;
	if (calc_x < w_length)
	return sigOfWhite;
	return sigOfBlack;

	}*/
	double startLine = 60;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w<0)
		return 0;
	return weight_w * sigOfWhite + (1 - weight_w)*sigOfBlack;


}

double BuprestidaeModelWithNoise2nd::calcEPSFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return realP * realP - imaginaryP * imaginaryP;

}

double BuprestidaeModelWithNoise2nd::calcSIGFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return 2 * realP*imaginaryP*w;

}


/*���G�ȑ��w��

���v20�w

in filed()
���f��
1�F���@150nm     30
1.5  5nm �@�@ (3��,2.05 ��)
2�F���@30nm      6
2.5  5nm�@    (2.55��,2.05��)�@
3�F���@115nm     23
3.5  5nm      (2.55��,0.55��)
4�F���@55nm      11
4.5  5nm      (1.75��,0.55��)
5�F���@95nm      19
5.5  5nm      (1.75��,2.413��)
6�F���@50nm      10
6.5  5nm      (1��,2.413��)
7�F���@85nm      17
7.5  5nm      (1��,1.85��)
8�F���@40nm      8
8.5  5nm      (0��,1.85��)
9�F���@90nm      18
9.5  5nm      (0��,1��)
10�F���@45nm     9
10.5 5nm      (1.225��,1��)

11�F���@80nm
11.5 5nm      (1.225��,1.86��)
12�F��  55nm
12.5 5nm      (1.225��,1.86��)



*/
double BuprestidaeModelWithNoise2nd::calcMultilayer(double x, double y, double start) {
	double layer = selectLayer(x, start);
	if (layer < 0)
		return -1;
	int layer_s = layer * 2 - 1;
	int a = layer_s % 4;
	if (a == 3)
		return 0;
	if (a == 1)
		return 1;

	switch (layer_s) {
	case 2:
		return 3 / (3 + 2.05);
	case 4:
		return 2.55 / (2.55 + 2.05);
	case 6:
		return 2.55 / (2.55 + 0.55);
	case 8:
		return 1.75 / (1.75 + 0.55);
	case 10:
		return 1.75 / (1.75 + 2.413);
	case 12:
		return 1.0 / (1.0 + 2.413);
	case 14:
		return 1.0 / (1.0 + 1.85);
	case 16:
		return 0 / (0 + 1.85);
	case 18:
		return 0 / (0 + 1.0);
	case 20:
		return 1.225 / (1.225 + 1.0);
	case 22:
		return 1.225 / (1.225 + 1.86);
	case 24:
		return 1.225 / (1.225 + 1.86);
	case 26:
		return 1.225 / (1.225 + 1.86);
	case 28:
		return 1.225 / (1.225 + 1.86);
	case 30:
		return 1.225 / (1.225 + 1.86);
	case 32:
		return 1.225 / (1.225 + 1.86);
	case 34:
		return 1.225 / (1.225 + 1.86);
	default:
		return 1.225 / (1.225 + 1.86);


	}
}



double BuprestidaeModelWithNoise2nd::selectLayer(double x, double start) {
	double length = x - start;
	if (length < 0)
		return -1;
	if (length > 304)//-11 20
		return -1;
	if (length > 293)//-1 
		return 20;
	if (length > 292)//-16 19
		return 19.5;
	if (length > 276)//-1
		return 19;
	if (length > 275)//-11 18
		return 19.5;
	if (length > 264)//-1
		return 18;
	if (length > 263)//-16 17
		return 17.5;
	if (length > 247)//-1
		return 17;
	if (length > 246)//-11 16
		return 16.5;
	if (length > 235)//-1
		return 16;
	if (length > 234)//-16 15
		return 15.5;
	if (length > 218)//-1
		return 15;
	if (length > 217)//-11 14
		return 14.5;
	if (length > 206)//-1
		return 14;
	if (length > 205)//-16 13
		return 13.5;
	if (length > 189)//-1
		return 13;
	if (length > 188)//-11 12
		return 12.5;
	if (length > 177)//-1
		return 12;
	if (length > 176)//-16 11
		return 11.5;
	if (length > 160)//-1
		return 11;
	if (length > 159)//-9  10
		return 10.5;
	if (length > 150)//-1
		return 10;
	if (length > 149)//-18 9
		return 9.5;
	if (length > 131)//-1
		return 9;
	if (length > 130)//-8  8
		return 8.5;
	if (length > 122)//-1
		return 8;
	if (length > 121)//-17  7
		return 7.5;
	if (length > 104)//-1
		return 7;
	if (length > 103)//-10  6
		return 6.5;
	if (length > 93)//-1
		return 6;
	if (length > 92)//-19  5
		return 5.5;
	if (length > 73)//-1
		return 5;
	if (length > 72)//-11  4
		return 4.5;
	if (length > 61)//-1
		return 4;
	if (length > 60)//-23  3
		return 3.5;
	if (length > 37)//-1
		return 3;
	if (length > 36)//-6 2
		return 2.5;
	if (length > 30)//-1
		return 2;
	if (length > 29)//-30
		return 1.5;
	return 1;
}

BuprestidaeModelSmooth24::BuprestidaeModelSmooth24(Field* f, double lambda) :FazzyModel(f) {
#define WHITE true
#define BLACK false
	//this->lambda = lambda;
	//this->w = w;
	setLambda(lambda);

	for (int i = 0; i < 7; i++) {
		realPartOfWhite[i] = initializeRoW(i);
		realPartOfBlack[i] = initializeRoB(i);
		imaginaryPartOfWhite[i] = initializeIoW(i);
		imaginaryPartOfBlack[i] = initializeIoB(i);
	}

	epOfWhite = calcEPSFromN(WHITE);
	epOfBlack = calcEPSFromN(BLACK);
	sigOfWhite = calcSIGFromN(WHITE);
	sigOfBlack = calcSIGFromN(BLACK);
	cout << "epw=" << epOfWhite << " epb=" << epOfBlack << endl;
	cout << "sigw=" << sigOfWhite << " sigb=" << sigOfBlack << endl;
}
double BuprestidaeModelSmooth24::initializeRoW(int i) {
	switch (i) {
	case 0:
		return 1.606;
	case 1:
		return 1.588;
	case 2:
		return 1.573;
	case 3:
		return 1.562;
	case 4:
		return 1.550;
	case 5:
		return 1.542;
	case 6:
		return 1.531;

	}
	return 0.0;
}

double BuprestidaeModelSmooth24::initializeRoB(int i) {
	switch (i) {
	case 0:
		return 1.784;
	case 1:
		return 1.745;
	case 2:
		return 1.718;
	case 3:
		return 1.683;
	case 4:
		return 1.669;
	case 5:
		return 1.658;
	case 6:
		return 1.640;

	}
	return 0.0;
}

double BuprestidaeModelSmooth24::initializeIoW(int i) {
	switch (i) {
	case 0:
		return 0.017;
	case 1:
		return 0.003;
	default:
		break;
	}
	return 0.0;
}

double BuprestidaeModelSmooth24::initializeIoB(int i) {
	switch (i) {
	case 0:
		return 0.096;
	case 1:
		return 0.069;
	case 2:
		return 0.050;
	case 3:
		return 0.034;
	case 4:
		return 0.025;
	case 5:
		return 0.021;
	case 6:
		return 0.013;

	}
	return 0.0;
}

string BuprestidaeModelSmooth24::mkdir(string root) {
	_mkdir((root + "BuprestidaeSmooth24").c_str());

	string name = "BuprestidaeSmooth24\\" + to_s(int(lambda));
	_mkdir((root + name).c_str());	//�f�B���N�g���̍쐬
	return name + "\\";
}

double BuprestidaeModelSmooth24::calcEPS(const double& x, const double& y, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return EPSILON_0_S;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x <length * 10) {
	int calc_x = (int)(mx - start_x) % length;
	if (calc_x < w_length)
	return epOfWhite;
	return epOfBlack;

	}*/

	double startLine = 60;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w < 0)
		return EPSILON_0_S;
	return weight_w * epOfWhite + (1 - weight_w)*epOfBlack;
}

double BuprestidaeModelSmooth24::calcSIG(const double& x, const double& y, const double lam, enum INTEG f) {
	double mx = x - mField->getNpml(); //�v�Z�̈���֎ʑ�
	double my = y - mField->getNpml();
	if (mx < 0 || my < 0 || mx >= mField->getNx() || my >= mField->getNy()) return 0;

	/*int w_length, b_length;
	w_length = mField->nanoToCell(97);
	b_length = mField->nanoToCell(66);
	int length = w_length + b_length;

	double length_x = mField->getNx();
	double start_x = length_x / 5;
	if (mx - start_x >= 0 && mx - start_x<length * 10) {
	int calc_x = (int)(mx - start_x) % length;
	if (calc_x < w_length)
	return sigOfWhite;
	return sigOfBlack;

	}*/
	double startLine = 60;
	double weight_w = calcMultilayer(mx, my, startLine);
	if (weight_w < 0)
		return 0;
	return weight_w * sigOfWhite + (1 - weight_w)*sigOfBlack;


}

double BuprestidaeModelSmooth24::calcEPSFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return realP * realP - imaginaryP * imaginaryP;

}

double BuprestidaeModelSmooth24::calcSIGFromN(bool white) {
	double realP, imaginaryP;
	if (white) {
		realP = getNByLambda(realPartOfWhite);
		imaginaryP = getNByLambda(imaginaryPartOfWhite);
	}
	else {
		realP = getNByLambda(realPartOfBlack);
		imaginaryP = getNByLambda(imaginaryPartOfBlack);
	}

	return 2 * realP*imaginaryP*w;

}


double BuprestidaeModelSmooth24::calcMultilayer(double x, double y, double start) {
	double layer = selectLayer(x, start);
	if (layer < 0)
		return -1;
	return (double)((int)layer % 2);

	

	
}



double BuprestidaeModelSmooth24::selectLayer(double x, double start) {
	double length = x - start;

	if (length < 0)
		return -1;
	if (length > 409)
		return -1;
	int lengthEx12 = int(length - 46);
	if (lengthEx12 > 0) {

		if (lengthEx12 % 33 < 19)
			return 3;
		return 4;
	}

	if ((int)length % 46 < 32)
		return 1;
	return 2;
	
}