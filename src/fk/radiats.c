/****************************************************************
   compute horizontal radiation pattens for
	double-couple	specified by az-strike,dip,rake
	single force	specified by az-strike and dip
	moment tnesor	specified by the moment tensor and az

*****************************************************************/
#include "radiats.h"

/*******************************************************************
horizontal radiation coefficients of a double-couple, Haskell'64
with tangential corrected

INPUT:	az of the station measured from the strike of the fault clockwise,
	dip, and rake of the fault-plane solution
OUTPUT:	rad[i][j] -> summation coef. for i-th azimuthal order and j-th component
			(0-> vertical, 1-> radial, 2->transverse)

Algorithm:
 V/R = f3n3*Z0
     +((f1n3+f3n1)*cos(theta)+(f2n3+f3n2)*sin(theta))*Z1
     +((f2n2-f1n1)*cos(2theta)+(-f1n2-f2n1)*sin(2theta))*Z2
  T =-((f1n3+f3n1)*sin(theta)-(f2n3+f3n2)*cos(theta))*T1
     -((f2n2-f1n1)*sin(2theta)-(-f1n2-f2n1)*cos(2theta))*T2

where theta=pi/2-az.
  n = (sin(delta),0,cos(delta))
  F = (-sin(lamda)*cos(delta), cos(lamda), sin(lamda)*sin(delta))
where delta is dip from the horizontal, lambda is the rake from the
strike CCW.

********************************************************************/
void	dc_radiat(float stk,float dip,float rak,float rad[4][3]) {
   float sstk,sdip,srak,sstk2,sdip2;
   float cstk,cdip,crak,cstk2,cdip2;
   stk*=DEG2RAD; dip*=DEG2RAD; rak*=DEG2RAD;
   sstk=sin(stk);cstk=cos(stk);
   sdip=sin(dip);cdip=cos(dip);
   srak=sin(rak);crak=cos(rak);
   sstk2=2*sstk*cstk; cstk2=cstk*cstk-sstk*sstk;
   sdip2=2*sdip*cdip; cdip2=cdip*cdip-sdip*sdip;
   rad[0][0]=0.5*srak*sdip2;
   rad[0][1]=rad[0][0];
   rad[0][2]=0.;
   rad[1][0]=-sstk*srak*cdip2+cstk*crak*cdip;
   rad[1][1]=rad[1][0];
   rad[1][2]= cstk*srak*cdip2+sstk*crak*cdip;
   rad[2][0]=-sstk2*crak*sdip-0.5*cstk2*srak*sdip2;
   rad[2][1]=rad[2][0];
   rad[2][2]=cstk2*crak*sdip-0.5*sstk2*srak*sdip2;
}

/******************************************************
horizontal radiation coefficients of a single-force
   In:
	stk: az_of_obs w.r.t to strike of the force
		 measured clockwise
	dip: dip of the force, from horizontal down

   algorithm:
	vertical (UP) = f3*Z0 + (f1*cos(theta)+f2*sin(theta))*Z1
	radial  (OUT) = f3*R0 + (f1*cos(theta)+f2*sin(theta))*R1
	tangen   (CW) =       - (f1*sin(theta)-f2*cos(theta))*T1
    where F = (0,cos(dip),-sin(dip))
******************************************************/
void	sf_radiat(float stk,float dip,float rad[4][3]) {
   float sstk,sdip,cstk,cdip;
   stk*=DEG2RAD; dip*=DEG2RAD;
   sstk=sin(stk);cstk=cos(stk);
   sdip=sin(dip);cdip=cos(dip);
   rad[0][0]=-sdip;
   rad[0][1]=rad[0][0];
   rad[0][2]=0.;
   rad[1][0]=cdip*cstk;
   rad[1][1]=rad[1][0];
   rad[1][2]=cdip*sstk;
}

/*****************************************************************
 horizontal radiation coefficients from a moment-tensor m
   see Jost and Herrmann, 1989 (note an error in Eq A5.4-A5.6)
*****************************************************************/
void	mt_radiat(float az, float m[3][3], float rad[4][3]) {
   float saz, caz, saz2, caz2;
   az *= DEG2RAD;
   saz = sin(az); caz = cos(az);
   saz2 = 2*saz*caz; caz2 = caz*caz-saz*saz;
   rad[2][0] = rad[2][1] = -0.5*(m[0][0]-m[1][1])*caz2 - m[0][1]*saz2; 
   rad[1][0] = rad[1][1] = -m[0][2]*caz - m[1][2]*saz;
   rad[0][0] = rad[0][1] =  (2*m[2][2]-m[0][0]-m[1][1])/6.;
   rad[2][2] = -0.5*(m[0][0]-m[1][1])*saz2 + m[0][1]*caz2;
   rad[1][2] = -m[0][2]*saz + m[1][2]*caz;
   rad[0][2] = 0.;
   /* contribution from explosion: */
   rad[3][0] = rad[3][1] = (m[0][0]+m[1][1]+m[2][2])/3.;
   rad[3][2] = 0;
}

/******************************************************************
 construct a normalized moment tensor (i.e. without the scalar M_0)
   M = sqrt(2/3)*iso*I + sqrt(1-iso*iso)*[sqrt(1-clvd*clvd)*DC+clvd*CLVD]
 where
   DC_ij = n_i v_j + v_i n_j,
   CLVD_ij = (2 N_i N_j - v_i v_j - n_i n_j )/sqrt(3),
   1>=iso>=-1,
   0.5>=clvd>=-0.5,
 from strike, dip, and rake (all in degrees). See A&R P117
******************************************************************/
void nmtensor(float iso, float clvd, float str, float dip, float rake, float tensor[3][3]) {
  float cstr,cdip,crak,sstr,sdip,srak,sstr2,cstr2,sdip2,cdip2,dum,dev;
  float n[3], v[3], N[3];
  dum = 0.8164966*iso;
  tensor[0][0] = tensor[1][1] = tensor[2][2] = dum;
  tensor[0][1] = tensor[0][2] = tensor[1][2] = 0.;
  dev = 1.-iso*iso;
  if (dev>0.) {
     dev = sqrt(dev);
     dum = dev*sqrt(1.-clvd*clvd);
     str  *= DEG2RAD; dip  *= DEG2RAD; rake *= DEG2RAD;
     sstr=sin(str);cstr=cos(str);sstr2=2*sstr*cstr;cstr2=1-2*sstr*sstr;
     sdip=sin(dip);cdip=cos(dip);sdip2=2*sdip*cdip;cdip2=1-2*sdip*sdip;
     crak=cos(rake);srak=sin(rake);
     tensor[0][0] += -dum*(sdip*crak*sstr2+sdip2*srak*sstr*sstr);
     tensor[0][1] +=  dum*(sdip*crak*cstr2+0.5*sdip2*srak*sstr2);
     tensor[0][2] += -dum*(cdip*crak*cstr+cdip2*srak*sstr);
     tensor[1][1] +=  dum*(sdip*crak*sstr2-sdip2*srak*cstr*cstr);
     tensor[1][2] +=  dum*(cdip2*srak*cstr-cdip*crak*sstr);
     tensor[2][2] +=  dum*sdip2*srak;
     if (clvd>0.0001 || clvd<-0.0001) {
        n[0] = -sdip*sstr; n[1] = sdip*cstr; n[2] = -cdip;
        v[0] =  crak*cstr+srak*cdip*sstr; v[1] = crak*sstr-srak*cdip*cstr; v[2] = -srak*sdip;
        N[0] = n[1]*v[2]-n[2]*v[1]; N[1] = n[2]*v[0]-n[0]*v[2]; N[2] = n[0]*v[1]-n[1]*v[0];
        //fprintf(stderr," n= %f %f %f  v = %f %f %f\n",n[0],n[1],n[2],v[0],v[1],v[2]);
        dum = dev*clvd/sqrt(3.);
        tensor[0][0] += dum*(2*N[0]*N[0]-n[0]*n[0]-v[0]*v[0]);
        tensor[0][1] += dum*(2*N[0]*N[1]-n[0]*n[1]-v[0]*v[1]);
        tensor[0][2] += dum*(2*N[0]*N[2]-n[0]*n[2]-v[0]*v[2]);
        tensor[1][1] += dum*(2*N[1]*N[1]-n[1]*n[1]-v[1]*v[1]);
        tensor[1][2] += dum*(2*N[1]*N[2]-n[1]*n[2]-v[1]*v[2]);
        tensor[2][2] += dum*(2*N[2]*N[2]-n[2]*n[2]-v[2]*v[2]);
     }
  }
  tensor[1][0] = tensor[0][1];
  tensor[2][0] = tensor[0][2];
  tensor[2][1] = tensor[1][2];
}

/*****************************************************************
  radiation pattern from a double couple,  A & R, P118
******************************************************************/
float radpmt(float mom[3][3], float alpha, float az, int type)
{
  float	dir[3], wave[3], sth, cth, cphi, sphi;
  int	m,n;
  
  sth=sin(alpha*DEG2RAD);
  cth=cos(alpha*DEG2RAD);
  cphi=cos(az*DEG2RAD);
  sphi=sin(az*DEG2RAD);
  dir[0]=sth*cphi;
  dir[1]=sth*sphi;
  dir[2]=cth;
  switch(type) {
  case 1:	/* P wave */
    wave[0]=dir[0];
    wave[1]=dir[1];
    wave[2]=dir[2];
    break;
  case 2:	/* SV wave */
    wave[0]=cth*cphi;
    wave[1]=cth*sphi;
    wave[2]=-sth;
    break;
  case 3:	/* SH wave */
    wave[0]=-sphi;
    wave[1]=cphi;
    wave[2]=0;
  break;
  default:
    fprintf(stderr,"wrong phase type %d\n",type);
    return 0.;
  }
  
  sth=0;
  for(m=0;m<3;m++){
    for(n=0;n<3;n++){
      sth+=mom[m][n]*dir[m]*wave[n];
    }
  }
  return sth;
}
