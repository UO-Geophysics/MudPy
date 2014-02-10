/******************************************************************
   F-K: @(#) syn.c             1.0 4/29/2000
 
   Copyright (c) 200 by L. Zhu
   See README file for copying and redistribution conditions.

 Synthesize the Green's functions by adding source radiation pattern

 Written by Lupei Zhu, 1996, seismo lab, Caltech

 The Green's functions are computed by fk.f

 Revision History:
   05/05/1996	Lupei Zhu	Initial coding.
   04/29/2000	Lupei Zhu	documatation.
   07/12/2000	Lupei Zhu	add component orientations in SAC.
   07/18/2000	Lupei Zhu	add -I option for integration.
   07/22/2000	Lupei Zhu	add option for convolving a trapezoid.
   03/29/2002	Lupei Zhu	add -P option for computing static disp.
   04/16/2004	Lupei Zhu	add -Mm for double-couple moment-tensor source.
   03/13/2006   Lupei Zhu	add -Me for explosion source.
   11/07/2006	Lupei Zhu	modify to input a general moment tensor source.
   11/01/2007	Lupei Zhu	add band-pass filtering (-F) and attenuation (-T) options.
   05/15/2008   Lupei Zhu	correct a bug introduced when using a general MT.
   02/13/2012	Lupei Zhu	correct a bug of writing sac head info.
******************************************************************/
#include <stdio.h> 
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "sac.h"
#include "Complex.h"
#include "radiats.h"

int main(int argc, char **argv) {
  int	i,j,k,nn,ns,npt,error,intg,diff,src_type, filter;
  char	nam[128],outnm[128],*ccc,com[3]={'z','r','t'};
  int	dynamic=1;
  float	coef,rad[4][3],m0,az,*grn,*syn[3],*src,*pt,disp[3],mt[3][3];
  float cmpinc[3]={0,90.,90.}, cmpaz[3];
  float dt, dura, rise, tp, ts, tmp, dist, shift;
  float	*trapezoid(float, float, float, int *);
  SACHEAD	hd;
#ifdef SAC_LIB
  char type[2] = {'B','P'}, proto[2] = {'B','U'};
  float	sn[30], sd[30];
  double f1, f2;
  long int order, nsects;
#endif
  void  fttq_(float *, float *, int *, int *, float *);
  int	mftm=2048, nftm;
  float	tstar=0., ftm[2048];

  /* input parameters */
  ns = 0;
  dura = 0.;
  src_type=0;
  intg=0;
  diff=0;
  filter=0;
  shift=0.;
  error = 0;
  for (i=1; !error && i < argc; i++) {
     if (argv[i][0] == '-') {
	switch(argv[i][1]) {

 	   case 'A':
	      sscanf(&argv[i][2], "%f",&az);
	      cmpaz[0] = 0.;
	      cmpaz[1] = az;
	      cmpaz[2] = az+90.;
	      if (cmpaz[2] > 360.) cmpaz[2] -= 360.;
	      break;

 	   case 'D':
	      j = sscanf(&argv[i][2], "%f/%f",&dura,&rise);
	      if (j<2) rise = 0.5;
	      break;

#ifdef SAC_LIB
 	   case 'F':
	      filter = 1;
	      j = sscanf(&argv[i][2], "%lf/%lf/%ld",&f1,&f2,&order);
	      if (j<3) order = 4;
	      break;
#endif

 	   case 'G':
	      strcpy(nam,&argv[i][2]);
	      break;

	   case 'I':
	      intg = 1;
	      break;

	   case 'J':
	      diff = 1;
	      break;

 	   case 'M':
	      src_type = sscanf(&argv[i][2], "%f/%f/%f/%f/%f/%f/%f",&m0,&mt[0][0],&mt[0][1],&mt[0][2],&mt[1][1],&mt[1][2],&mt[2][2]);
	      break;

 	   case 'O':
	      strcpy(outnm, &argv[i][2]);
	      break;

	   case 'P':
	      dynamic = 0;
	      break;

	   case 'Q':
	      sscanf(&argv[i][2], "%f",&tstar);
	      break;

 	   case 'S':
	      if ( (src=read_sac(&argv[i][2],&hd)) != NULL ) {
                 ns = hd.npts;
		 shift = -hd.b;
	      }
	      break;

	   default:
	      error = 1;
	      break;
	}
     }
        else error = 1;
  }

  switch (src_type) {
  case 1:
     nn = 1;
     m0 = m0*1.0e-20;
     break;
  case 3:
     nn = 2;
     sf_radiat(az-mt[0][0],mt[0][1],rad);
     m0 = m0*1.0e-15;
     break;
  case 4:
     nn = 3;
     dc_radiat(az-mt[0][0],mt[0][1],mt[0][2],rad);
     m0 = pow(10.,1.5*m0+16.1-20); 
     break;
  case 7:
     nn = 4;
     mt_radiat(az,mt,rad);
     m0 = m0*1.0e-20;
     break;
  default:
     error = 1;
  }

  if (dynamic) ccc = strrchr(nam, (int) '.') + 1;

  if(argc == 1 || error || (dynamic && src_type == 1 && (*ccc) != 'a') ) {
     fprintf(stderr,"Usage: %s -Mmag([[/Strike/Dip]/Rake]|/Mxx/Mxy/Mxz/Myy/Myz/Mzz) -Aazimuth ([-SsrcFunctionName | -Ddura[/rise]] [-Ff1/f2[/n]] [-I | -J] -OoutName.z -GFirstCompOfGreen | -P)\n\
   Compute displacements in cm in the up, radial, and transverse (clockwise) directions produced by difference seismic sources\n\
   -M Specify source magnitude and orientation or moment-tensor\n\
      For double-couple, mag is Mw, strike/dip/rake are in A&R convention\n\
      For explosion; mag in in dyne-cm, no strike, dip, and rake needed\n\
      For single-force source; mag is in dyne, only strike and dip are needed\n\
      For moment-tensor; mag in dyne-cm, x=N,y=E,z=Down\n\
   -A Set station azimuth in degree measured from the North\n\
   -D Specify the source time function as a trapezoid,\n\
      give the total duration and rise-time (0-0.5, default 0.5=triangle)\n\
   -F apply n-th order Butterworth band-pass filter, SAC lib required (off, n=4, must be < 10)\n\
   -G Give the name of the first component of the FK Green function\n\
   -I Integration once\n\
   -J Differentiate the synthetics\n\
   -O Output SAC file name\n\
   -P Compute static displacement, input Green functions from stdin in the form\n\
	distance Z45 R45 T45 ZDD RDD TDD ZSS RSS TSS [distance ZEX REX TEX]\n\
      The displacements will be output to stdout in the form of\n\
	distance azimuth z r t\n\
   -Q Convolve a Futterman Q operator of tstar (no)\n\
   -S Specify the SAC file name of the source time function (its sum. must be 1)\n\
   Examples:\n\
   * To compute three-component velocity at N33.5E azimuth from a Mw 4.5\n\
earthquake (strike 355, dip 80, rake -70), use:\n\
	syn -M4.5/355/80/-70 -D1 -A33.5 -OPAS.z -Ghk_15/50.grn.0\n\
   * To compute the static displacements from the same earthquake, use:\n\
	nawk \'$1==50\' st.out | syn -M4.5/355/80/-70 -A33.5 -P\n\
   * To compute displacement from an explosion, use:\n\
   	syn -M3.3e20 -D1 -A33.5 -OPAS.z -Ghk_15/50.grn.a\n\
      or\n\
        syn -M3.3e20/1/0/0/1/0/1 -D1 -A33.5 -OPAS.z -Ghk_15/50.grn.0\n\
	\n",argv[0]);
    return -1;
  }

  for(j=0; j<3; j++) disp[j] = 0.;
  for(i=0; i<nn; i++) {	/* sum up contribution from a DD, DS, SS, and EX */
     for(j=0; j<3; j++) {
	coef = m0; if (nn>1) coef *= rad[i][j];
	if (!dynamic) {
	   if ((i==0||i==3) && j==0) scanf("%f",&dist);
	   scanf("%f",&tmp);
	   disp[j] += coef*tmp;
	   continue;
	}
        if ( (grn=read_sac(nam,&hd)) == NULL ) continue;
	if ( i==0 && j== 0 ) {
           npt = hd.npts;
	   dt = hd.delta;
	   tp = hd.t1; ts = hd.t2;
  	   for(k=0; k<3; k++) syn[k]=(float *) calloc(npt, sizeof(float));
	   if (dura>0.) src = trapezoid(dura,rise,dt,&ns);
	}
	else if (hd.npts != npt) {
	   fprintf(stderr,"number points in %s not agree with %d\n",nam,npt);
	}
        for(pt=syn[j],k=0;k<npt;k++,pt++) (*pt) += coef*grn[k];
	free(grn);
        (*ccc)++;
	if (*ccc == '9') (*ccc)+=40;	/* explosion components start at 'a' instead of '0' */
     }
  }

  if (!dynamic) {
     printf("%8.2f %8.2f %10.3e %10.3e %10.3e\n",dist,az,disp[0],disp[1],disp[2]);
     return 0;
  }

  /* convolve a source time function. integrate or filtering if needed */
#ifdef SAC_LIB
  if (filter) design(order, type, proto, 1., 1., f1, f2, (double) dt, sn, sd, &nsects);
#endif
  if (tstar>0.) fttq_(&dt,&tstar,&mftm,&nftm,ftm);
  for(j=0;j<3;j++) {
     if (ns > 0) conv(src,ns,syn[j],npt);
     if (intg) cumsum(syn[j],npt,dt);
     if (diff) diffrt(syn[j],npt,dt);
#ifdef SAC_LIB
     if (filter) apply(syn[j],(long int)npt,0,sn,sd,nsects);
#endif
     if (tstar>0.) conv(ftm,nftm,syn[j],npt);
  }

  /* output */
  ccc = outnm + strlen(outnm) - 1;
  if (ns == 0) strcat(outnm,"i");
  hd.b -= shift;
  hd.e -= shift;
  hd.t1 = tp; hd.t2 = ts;
  hd.az = az;
  for(j=0;j<3;j++) {
     *ccc = com[j];
     hd.cmpinc = cmpinc[j];
     hd.cmpaz  = cmpaz[j];
     write_sac(outnm,hd,syn[j]);
  }

  return 0;

}

/* construct a trapezoid of duration dura with rise portion rise */
/* and an area of 1*dt */
float *trapezoid(float dura, float rise, float dt, int *ns) {
   int i, nr, trap_len;
   float amp, *src;
   FILE *filePtr;   /*DM Added */
   *ns = rint(dura/dt);
   trap_len=0; /*DM*/ 
   if (*ns<2) *ns = 2;src = malloc((1+*ns)*sizeof(float));
   if (src == NULL) {*ns=0; return src;}
   nr = rint(rise*(*ns)); if (nr<1) nr = 1; if (2*nr>*ns) nr=*ns/2;
   amp = 1./(nr*(*ns-nr));
   for(i=0;i<nr;i++) {
       src[i] = i*amp;
       trap_len++;
   }
   for(;i<*ns-nr;i++) {
       src[i] = nr*amp;
       trap_len++;
   }
   for(;i<=*ns;i++)   {
       src[i] = (*ns-i)*amp;
       trap_len++;
   }
   (*ns)++;
   /* Diego says: Write the thing to a file so you can look at it yo*/
   filePtr = fopen("floatArray.txt","w");
   for (i = 0; i < trap_len; i++) {
      fprintf(filePtr, "%.5g\n", src[i]);
   }
   fclose(filePtr);
   /**/
   return src;
}
