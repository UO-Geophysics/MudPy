/*******************************************************************
*			sacio.c
* SAC I/O functions:
*	read_sachead	read SAC header
*	read_sac	read SAC binary data
*	read_sac2	read SAC data with cut option
*	write_sac	write SAC binary data
*	swab4		reverse byte order for integer/float
*	wrtsac2		write 2 1D arrays as XY SAC data
*	sachdr		creat new sac header
*	rdsac0_		fortran wraper for read_sac
*	my_brsac_	fortran binary sac data+head reader.
*	wrtsac0_	fortran write 1D array as SAC binary data
*	wrtsac2_	fortran wraper for wrtsac2
*	wrtsac3_	wrtsac0 with component orientation cmpaz/cmpinc
*	ResetSacTime	reset reference time in the head
*	sac_head_index	find the index of a sac head field.

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
	08/31/04	Lupei Zhu	read_sac2 uses t2 instead of number of samples;
			fill zeros if it's not completely in the window.
	05/22/12	Lupei Zhu	read_sac auto detect byte swap; always write in local byte order.
	02/08/13	Lupei Zhu	correct a bug in wrtsac2 (SAC puts y first followed by x).
*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac.h"

/***********************************************************

  read_sachead

  Description:	read binary SAC header from file.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled

  Return:	0 if success, -1 if failed

************************************************************/

int	read_sachead(const char *name, SACHEAD *hd) {
  FILE		*strm;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return -1;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     fclose(strm);
     return -1;
  }

  if ( testByte((char *)hd) ) swab4((char *) hd, HD_SIZE);

  fclose(strm);
  return 0;

}


/***********************************************************

  read_sac

  Description:	read binary data from file. If succeed, it will return
		a float pointer to the read data array. The SAC header
		is also filled. A NULL pointer is returned if failed.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled

  Return:	float pointer to the data array, NULL if failed

************************************************************/

float*	read_sac(const char *name, SACHEAD *hd) {
  FILE		*strm;
  float		*ar;
  unsigned	sz;
  int		byteSwap;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

  if ( (byteSwap=testByte((char *)hd)) ) swab4((char *) hd, HD_SIZE);

  sz = hd->npts*sizeof(float);
  if (hd->iftype == IRLIM || hd->iftype == IAMPH || hd->iftype == IXY) sz *= 2;
  if ((ar = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     return NULL;
  }

  if (fread((char *) ar, sz, 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }

  fclose(strm);

  if ( byteSwap ) swab4((char *) ar, sz);

  return ar;

}


/***********************************************************

  write_sac

  Description:	write SAC binary data.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD hd		SAC header
		const float *ar		pointer to the data

  Return:	0 if succeed; -1 if failed

************************************************************/

int	write_sac(const char *name, SACHEAD hd, const float *ar) {
  FILE		*strm;
  int		error = 0;
  unsigned	sz;
  sz = hd.npts*sizeof(float);
  if (hd.iftype == IRLIM || hd.iftype == IAMPH || hd.iftype == IXY) sz *= 2;

  if ( !error && (strm = fopen(name, "w")) == NULL ) {
     fprintf(stderr,"Error in opening file for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(&hd, sizeof(SACHEAD), 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC header for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(ar, sz, 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC data for writing %s\n",name);
     error = 1;
  }
  fclose(strm);

  return (error==0) ? 0 : -1;

}


/*****************************************************

  swab4

  Description:	reverse byte order for float/integer

  Author:	Lupei Zhu

  Arguments:	char *pt	pointer to byte array
		int    n	number of bytes

  Return:	none

*******************************************************/

void	swab4(char *pt, int n) {
  int i;
  char temp;
  for(i=0;i<n;i+=4) {
    temp = pt[i+3];
    pt[i+3] = pt[i];
    pt[i] = temp;
    temp = pt[i+2];
    pt[i+2] = pt[i+1];
    pt[i+1] = temp;
  }
}


/***********************************************************

  sachdr

  Description:	creat a new SAC header

  Author:	Lupei Zhu

  Arguments:	float	dt		sampling interval
		int	ns		number of points
		float	b0		starting time

  Return:	SACHEAD

************************************************************/

SACHEAD    sachdr(float dt, int ns, float b0) {
  SACHEAD	hd = sac_null;
  hd.npts = ns;
  hd.delta = dt;
  hd.b = b0;
  hd.o = 0.;
  hd.e = b0+(ns-1)*hd.delta;
  hd.iztype = IO;
  hd.iftype = ITIME;
  hd.leven = TRUE;
  return hd;
}


/***********************************************************

  wrtsac2

  Description:	write 2 arrays into XY SAC data.

  Author:	Lupei Zhu

  Arguments:	const char *name	file name
		int	ns		number of points
		const float *x		x data array
		const float *y		y data array

  Return:	0 succeed, -1 fail

************************************************************/

int	wrtsac2(const char	*name,
		int		n,
		const float	*x,
		const float	*y
	) {
  SACHEAD	hd = sac_null;
  float		*ar;
  unsigned	sz;
  int		exit_code;

  hd.npts = n;
  hd.iftype = IXY;
  hd.leven = FALSE;

  sz = n*sizeof(float);

  if ( (ar = (float *) malloc(2*sz)) == NULL ) {
     fprintf(stderr, "error in allocating memory%s\n",name);
     return -1;
  }

  if (memcpy(ar, y, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }
  if (memcpy(ar+n, x, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }

  exit_code = write_sac(name, hd, ar);
  
  free(ar);
  
  return exit_code;

}


/*for fortran--read evenly-spaced data */
void	rdsac0_(const char *name, float *dt, int *ns, float *b0, float *ar) {
   int i;
   SACHEAD hd;
   float *temp;
   temp = read_sac(name,&hd);
   *dt = hd.delta;
   *ns = hd.npts;
   *b0 = hd.b;
   for(i=0;i<*ns;i++) ar[i]=temp[i];
   free(temp);
}


/*for fortran--read evely-spaced data with header, make sure hdr size >= 158 */
void	my_brsac_(char *name, float *hdr, int *hdi, char *hdc, float *ar, int *err) {
   int i, *ipt;
   float *temp;
   char *cpt;
   cpt = strchr(name, (int) ' '); *cpt = 0;
   *err = 0;
   temp = read_sac(name, (SACHEAD *)hdr);
   if (temp == NULL) {*err = -1; return;}
   ipt = (int *) (hdr + 70);
   for(i=0;i<30;i++) hdi[i] = ipt[i];
   cpt = (char *)  (hdr + 110);
   for(i=0;i<192;i++) hdc[i] = cpt[i];
   for(i=0;i<hdi[9];i++) ar[i]=temp[i];
   free(temp);
}


/* for fortran--write evenly-spaced data */
void    wrtsac0_(const char *name, float *dt, int *ns, float *b0, float *dist, const float *ar) {
  SACHEAD hd;
  hd = sachdr(*dt,*ns,*b0);
  hd.dist = *dist;
  write_sac(name, hd, ar);
}


/* for fortran--write x-y data */
void    wrtsac2_(const char *name, int n, const float *x, const float *y) {
  wrtsac2(name, n, x, y);
}


/* for fortran--write evenly-spaced data with comp orientation */
void    wrtsac3_(const char *name, float dt, int ns, float b0, float dist, float cmpaz, float cmpinc, const float *ar) {
  SACHEAD hd;
  hd = sachdr(dt,ns,b0);
  hd.dist = dist;
  hd.cmpaz=cmpaz;
  hd.cmpinc=cmpinc;
  write_sac(name, hd, ar);
}


/***********************************************************

  read_sac2

  Description:	read portion of data from file. If succeed, it will return
		a float pointer to the read data array which may be
		zero-filled if the SAC trace is partially or not in the 
                specified window. The SAC head is also filled.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled
		int	tmark,		time mark in sac header
					0-9 -> tn; -5 -> b; -3 -> o; -2 -> a
					others refer to t=0.
		float	t1		begin time is tmark + t1
		float	t2		end time is tmark + t2

  Return:	float pointer to the data array, NULL if failed

************************************************************/

float*	read_sac2(const char	*name,
		SACHEAD		*hd,
		int		tmark,
		float		t1,
		float		t2
	) {
  FILE		*strm;
  int		nn, nt1, nt2, npts, byteSwap;
  float		tref, *ar, *fpt;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

  if ( (byteSwap=testByte((char *)hd)) ) swab4((char *) hd, HD_SIZE);

  nn = (int) rint((t2-t1)/hd->delta);
  if (nn<=0 || (ar = (float *) calloc(nn,sizeof(float)))==NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s n=%d\n",name,nn);
     return NULL;
  }
  tref = 0.;
  if (tmark==-5 || tmark==-3 || tmark==-2 || (tmark>=0&&tmark<10) ) {
     tref = *( (float *) hd + 10 + tmark);
     if (tref==-12345.) {
        fprintf(stderr,"Time mark undefined in %s\n",name);
        return NULL;
     }
  }
  t1 += tref;
  nt1 = (int) rint((t1-hd->b)/hd->delta);
  nt2 = nt1+nn;
  npts = hd->npts;
  hd->npts = nn;
  hd->b = t1;
  hd->e = t1+nn*hd->delta;

  if ( nt1>=npts || nt2<0 ) return ar;

  if (nt1<0) {
     fpt = ar-nt1;
     nt1 = 0;
  } else {
     if (fseek(strm,nt1*sizeof(float),SEEK_CUR) < 0) {
	fprintf(stderr, "error in seek %s\n",name);
	return NULL;
     }
     fpt = ar;
  }
  if (nt2>npts) nt2=npts;
  nn = nt2-nt1;
  if (fread((char *) fpt, sizeof(float), nn, strm) != nn) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }
  fclose(strm);
  if (byteSwap) swab4((char *) fpt, nn*sizeof(float));
  return ar;

}


/*reset reference time in sac header*/
void	ResetSacTime(SACHEAD *hd) {
     hd->o = 0.;
     hd->a = 0.;
     hd->nzyear = -12345;
     hd->nzjday = -12345;
     hd->nzhour = -12345;
     hd->nzmin = -12345;
     hd->nzsec = -12345;
     hd->nzmsec = -12345;
}


/*find the index of sac head field name, return -1 if not found*/
int	sac_head_index(const char *name) {
  if(strcmp(name,"delta")==0) return(0);
  else if (strcmp(name,"depmin")==0) return(1);
  else if (strcmp(name,"depmax")==0) return(2);
  else if (strcmp(name,"b")==0) return(5);
  else if (strcmp(name,"e")==0) return(6);
  else if (strcmp(name,"o")==0) return(7);
  else if (strcmp(name,"a")==0) return(8);
  else if (strcmp(name,"t0")==0) return(10);
  else if (strcmp(name,"t1")==0) return(11);
  else if (strcmp(name,"t2")==0) return(12);
  else if (strcmp(name,"t3")==0) return(13);
  else if (strcmp(name,"t4")==0) return(14);
  else if (strcmp(name,"t5")==0) return(15);
  else if (strcmp(name,"t6")==0) return(16);
  else if (strcmp(name,"t7")==0) return(17);
  else if (strcmp(name,"t8")==0) return(18);
  else if (strcmp(name,"t9")==0) return(19);
  else if (strcmp(name,"stla")==0) return(31);
  else if (strcmp(name,"stlo")==0) return(32);
  else if (strcmp(name,"stel")==0) return(33);
  else if (strcmp(name,"stdp")==0) return(34);
  else if (strcmp(name,"evla")==0) return(35);
  else if (strcmp(name,"evlo")==0) return(36);
  else if (strcmp(name,"evel")==0) return(37);
  else if (strcmp(name,"evdp")==0) return(38);
  else if (strcmp(name,"user0")==0) return(40);
  else if (strcmp(name,"user1")==0) return(41);
  else if (strcmp(name,"user2")==0) return(42);
  else if (strcmp(name,"user3")==0) return(43);
  else if (strcmp(name,"user4")==0) return(44);
  else if (strcmp(name,"user5")==0) return(45);
  else if (strcmp(name,"user6")==0) return(46);
  else if (strcmp(name,"user7")==0) return(47);
  else if (strcmp(name,"user8")==0) return(48);
  else if (strcmp(name,"user9")==0) return(49);
  else if (strcmp(name,"dist")==0) return(50);
  else if (strcmp(name,"az")==0) return(51);
  else if (strcmp(name,"baz")==0) return(52);
  else if (strcmp(name,"gcarc")==0) return(53);
  else if (strcmp(name,"depmen")==0) return(56);
  else if (strcmp(name,"cmpaz")==0) return(57);
  else if (strcmp(name,"cmpinc")==0) return(58);
  else if (strcmp(name,"kztime")==0) return(70);
  else if (strcmp(name,"npts")==0) return(79);
  else if (strcmp(name,"kstnm")==0) return(110);
  else return -1;
}


int testByte(char *hd) {
   int ver = *((int *) (hd+304));
   if (ver<1||ver>6) {
	   swab4((char *)&ver,4);
	   if (ver<1||ver>6) {
		   fprintf(stderr,"****** not in sac format, version = %d\n",ver);
	   }
	   return 1;
   }
   return 0;
}
