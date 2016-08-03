c***************************************************************
c  F-K: @(#) layer.h			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c  include header for transfering parameters of any given layer
c  from one subroutines to another
c***************************************************************
      include 'constants.h'
      real*8 kd,exa,exb,mu2,y,y1
      complex*16 ra,rb,Ca,Cb,Ya,Xa,Yb,Xb,r,r1
      common /layer/kd,mu2,exa,exb,ra,rb,Ca,Cb,Ya,Xa,Yb,Xb,r,r1,y,y1
