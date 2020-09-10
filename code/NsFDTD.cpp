#include "FDTD.h"

NsFDTD::NsFDTD()
:FDTD(){
		cout << "NsFDTD Constructor" << endl;
}

NsFDTD::~NsFDTD(){
	cout << "NsFDTD Destructor" << endl;
}

//NsFDTD‚Ì‚Æ‚«‚Ì·•ª–@
complex<double> NsFDTD::DxDy2(complex<double> *p, int i, int j, int t){
	return     p[index(i+1,j-1, t)] + p[index(i+1,j+1, t)] + p[index(i-1,j-1, t)] + p[index(i-1,j+1, t)]
	      -2.0*( p[index(i+1,j  , t)] + p[index(i-1,j  , t)] + p[index(i  ,j+1, t)] + p[index(i  ,j-1, t)] )
		  +4.0*p[index(i,j, t)];
}

complex<double> NsFDTD::D0_2(complex<double> *p, int i, int j, int t){
	return Dx2(p, i,j,t) + Dy2(p, i,j,t) + r_s*DxDy2(p, i,j,t);
}

void NsFDTD::field(){
	//‹üÜ—¦‚Ìİ’è
	//ŒvZ—p’è”‚Ìİ’è
	kx_s = 1/sqrt(sqrt(2.0)) * k_s;
	ky_s = sqrt(1 - 1/sqrt(2.0) ) * k_s;

	double sin2_kx = pow(sin(kx_s*H_S/2), 2);
	double sin2_ky = pow(sin(ky_s*H_S/2), 2);
	double sin2_k  = pow(sin(k_s *H_S/2), 2);
	r_s = (sin2_kx + sin2_ky - sin2_k)/(4*sin2_kx*sin2_ky);  //(1-ƒÁ0)/2

	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	//double kx_b = k_b * cos(PI/4);
	//double ky_b = k_b * sin(PI/4)

	for(int i=0; i<mField->getNx(); i++)
		for(int j=0; j<mField->getNy();j++)
			np[index(i,j)] = _pow(sin(w_b/n_s[index(i,j)]) / sin(k_b),2);
}

bool NsFDTD::calc(){
	time += DT_S;		//Ÿ‚ÌƒXƒeƒbƒv‚ÌŒvZ
	//printf("%lf \n", time/T_s);

	//‹«ŠEˆÈŠO
	for(int i=1; i<mField->getNx()-1; i++){
		for(int j=1; j<mField->getNy()-1; j++){
			if( (i==1 || i==mField->getNx()-2) && (j==1 || j==mField->getNy()-2) ) //l‹÷‚ğQÆ‚µ‚Ä‚µ‚Ü‚¤êŠ‚ÍS-FDTD‚Å‰ğ‚­
				phi[index(i,j, +1)] = np[index(i,j)]*( Dx2(phi, i,j, 0) + Dy2(phi, i,j, 0) ) - phi[index(i,j, -1)] + 2.0*phi[index(i,j,0)];
			else 
				phi[index(i,j, +1)] = np[index(i,j)]*D0_2(phi, i,j, 0) - phi[index(i,j, -1)] + 2.0*phi[index(i,j,0)];
		}
	}

//‹zû‹«ŠE
	absorbing();

	//phi[index(mField->getNx()/2, mField->getNy()/2, +1)] += 5*(1-exp(-0.01*time*time))*sin(- w_s*time);
	pointLightSource(phi);
	ButtonFactory::setButton("time", time);
	if(time > 2500){
		MiePrint(phi, "Mie_Nsphi2");
		return false;
	}

	return true;
}


//-------------------------------------------------------//
//-----------«---‹zû‹«ŠE‚±‚Á‚©‚ç «--------------------//
//-------------------------------------------------------//
//¶•Ç‚Ì‹zû‹«ŠE
void NsFDTD::absorbing_left(){
	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int i=1; i<mField->getNx()-1; i++){

		u1 = tan(w_b/n_s[index(0,i)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/n_s[index(0,i)]), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(i==1 || i==mField->getNx()-2)	//l‹÷‚Ì‰¡‚ÍˆêŸŒ³‹zû‹«ŠE
			phi[index(0,i, +1)]    = phi[index(1,i, 0)]    + (1- u1)/(1+u1)*(phi[index(0,i, 0)]    - phi[index(1,i, +1)]);

		else				//“ñŸŒ³‹zû‹«ŠE
			phi[index(0,i, +1)] = -phi[index(1,i, -1)]     
								 -  (1-u1)/(1+u1)*(phi[index(0,i, -1)] + phi[index(1,i, +1)])    
							     +       2/(1+u1)*(phi[index(0,i,  0)] + phi[index(1,i,  0)]) 
								 + u2*u2/(1+u1)/2*(   Dy2(phi,0,i, 0)  + Dy2(phi,1,i, 0)	);
												//  dy^2 ƒÓn     +   dy^2 ƒÓb
	}			
}

//‰E•Ç‚Ì‹zû‹«ŠE
void NsFDTD::absorbing_right(){

	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int i=1; i<mField->getNx()-1; i++){
		u1 = tan(w_b/n_s[index(mField->getNx()-1,i)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/n_s[index(mField->getNx()-1,i)]), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(i == 1 || i == mField->getNx()-2)		// l‹÷‚Ì‰¡‚ÍˆêŸŒ³‹zû‹«ŠE
			phi[index(mField->getNx()-1,i, +1)] = phi[index(mField->getNx()-2,i, 0)] + (1- u1)/(1+u1)*(phi[index(mField->getNx()-1,i, 0)] - phi[index(mField->getNx()-2,i, +1)]);

		else						//‚»‚êˆÈŠO‚Í“ñŸŒ³‹zû‹«ŠE
			phi[index(mField->getNx()-1,i, +1)] = - phi[index(mField->getNx()-2,i, -1)] 
								     -   (1-u1)/(1+u1)*(phi[index(mField->getNx()-1,i, -1)] + phi[index(mField->getNx()-2,i, +1)]) 
								     +        2/(1+u1)*(phi[index(mField->getNx()-1,i,  0)] + phi[index(mField->getNx()-2,i,  0)]) 
								     + u2*u2/(1+u1)/2*(   Dy2(phi, mField->getNx()-1,i, 0) + Dy2(phi, mField->getNx()-2,i, 0)	);
												        //  dy^2 ƒÓn     +   dy^2 ƒÓb
	}
}

//ã•Ç‚Ì‹zû‹«ŠE
void NsFDTD::absorbing_up(){

	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int i=1; i<mField->getNx()-1; i++){
		u1 = tan(w_b/n_s[index(i,0)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/n_s[index(i,0)]), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(i==1 || i==mField->getNx()-2)	//l‹÷‚Ì‰¡‚ÍˆêŸŒ³‹zû‹«ŠE
			phi[index(i,0, +1)]    = phi[index(i,1, 0)]    + (1- u1)/(1+u1)*(phi[index(i,0, 0)]    - phi[index(i,1, +1)]);

		else				//“ñŸŒ³‹zû‹«ŠE
			phi[index(i,0, +1)]  = -phi[index(i,1, -1)]    
								  -  (1-u1)/(1+u1)*(phi[index(i,0, -1)]+phi[index(i,1, +1)]) 
							      +       2/(1+u1)*(phi[index(i,0, 0)]+phi[index(i,1, 0)])     
								  + u2*u2/(1+u1)/2*( Dx2(phi, i,0, 0) + Dx2(phi, i,1, 0) 	);
												  //  dx^2 ƒÓn     +   dx^2 ƒÓb
	}		
}

//‰º•Ç‚Ì‹zû‹«ŠE
void NsFDTD::absorbing_down(){

	double u1, u2;
	double w_b  = w_s*DT_S/2;
	double k_b  = k_s*H_S/2;
	double kx_b = kx_s*H_S/2;
	double ky_b = ky_s*H_S/2;

	for(int i=1; i<mField->getNx()-1; i++){

		u1 = tan(w_b/n_s[index(i,mField->getNx()-1)]) / tan(k_b);
		u2 = 2 * _pow(sin(w_b/n_s[index(i,mField->getNx()-1)]), 2) / _pow(sin(ky_b),2) * (1 - tan(kx_b)/tan(k_b));

		if(i==1 || i==mField->getNx()-2)	//l‹÷‚Ì‰¡‚ÍˆêŸŒ³‹zû‹«ŠE
			phi[index(i,mField->getNx()-1, +1)] = phi[index(i,mField->getNx()-2, 0)] + (1- u1)/(1+u1)*(phi[index(i,mField->getNx()-1, 0)] - phi[index(i,mField->getNx()-2, +1)]);



		else				//“ñŸŒ³‹zû‹«ŠE
			phi[index(i,mField->getNx()-1, +1)] = -phi[index(i,mField->getNx()-2, -1)] 
			                      - (1-u1)/(1+u1)*(phi[index(i,mField->getNx()-1, -1)]+phi[index(i,mField->getNx()-2, +1)]) 
							      +     2/(1+u1)*(phi[index(i,mField->getNx()-1,  0)]+phi[index(i,mField->getNx()-2,  0)]) 
								  + u2*u2/(1+u1)/2*( Dx2(phi, i,mField->getNx()-1, 0) + Dx2(phi, i,mField->getNx()-2, 0)	 );
													 //  dx^2 ƒÓn     +   dx^2 ƒÓb
	}		
}

//‹zû‹«ŠE
void NsFDTD::absorbing(){
	/*
	absorbing_left();
	absorbing_right();
	absorbing_up();
	absorbing_down();	
	*/
	absorbing_nsRL(phi, 0, LEFT);
	absorbing_nsRL(phi, mField->getNx()-1, RIGHT);
	absorbing_nsTB(phi, 0,    BOTTOM);
	absorbing_nsTB(phi, mField->getNy()-1, TOP);
}

//-------------------------------------------------------//
//-------------ª---‹zû‹«ŠE‚±‚±‚Ü‚Å---ª----------------//
//-------------------------------------------------------//


void NsFDTD::Mie_Cylinder_Incidence(){
	//’†S‚Ì‰~‚É,U—”g”­¶
	for(int i=0; i<mField->getNx(); i++)
		for(int j=0; j<mField->getNx(); j++)	
	//		if((i-mField->getNx()/2)*(i-mField->getNx()/2) + (j-mField->getNx()/2)*(j-mField->getNx()/2) <= lambda_s*lambda_s)
				phi[index(i,j, +1)] += (1-exp(-0.01*time*time))*(1/_pow(n_s[index(i,j)], 2) - 1) * ( sin(k_s*i - w_s*(time + DT_S) ) + sin(k_s*i - w_s*(time - DT_S)) - 2*sin(k_s*i - w_s*time)); 
	
}
