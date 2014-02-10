/****************************************************************
   Head file for radiats.c

   rad[4][3] is the horizontal radiation rad[DD,DS,SS,EX][Z,R,T]
*****************************************************************/

#ifndef __RADIATS_HEAD__
  #define __RADIATS_HEAD__

#include<stdio.h>
#include<math.h>

#define DEG2RAD  1.745329252e-2  /*degree to rad*/

void    dc_radiat(float stk,float dip,float rak,float rad[4][3]);
void	sf_radiat(float, float, float rad[4][3]);
void	mt_radiat(float, float mt[3][3], float rad[4][3]);
void	nmtensor(float iso, float clvd, float stk, float dip, float rak, float mt[3][3]);
float	radpmt(float mt[3][3], float alpha, float az, int type);

#endif
