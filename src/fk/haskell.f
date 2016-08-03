c***************************************************************
c  F-K: @(#) haskell.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c  This file contains the following subroutines for computing Haskell propagator
c  matrix and other related matrices/vectors for the dynamic case.
c	layerParameter()	compute some parameters for each layer.
c       sh_ch()			compute Cosh(x), Sinh(x)/x, and x*Sinh(x).
c	haskell()		compute the P-SV Haskell matrix.
c       compoundMatrix()	compute the compound matrix of Haskell matrix.
c	eVector()		compute the compound matrix E|_12^ij.
c       initialG()		initialize g vector in the bottom half space.
c***************************************************************


      subroutine layerParameter(k, lay)
c***************************************************************
c compute some parameters for this layer
c	IN:
c		k   --- wave-number
c		lay --- layer number
c		volocity model passed in common/model/
c	Out:
c		common/layer/
c called by: kernel()		in kernel.f
c***************************************************************
      IMPLICIT NONE
      integer	lay
      real	k
      real*8	k2
      complex*16	kka, kkb
      include 'layer.h'
      include 'model.h'
      k2 = k*k
      kka = ka(lay)/k2
      kkb = kb(lay)/k2
      r = two/kkb
      kd = k*d(lay)
      mu2 = 2.*mu(lay)
      ra = zsqrt(one - kka)
      rb = zsqrt(one - kkb)
      r1 = one - one/r
      return
      end


      subroutine sh_ch(c,y,x,ex,a,kd)
c***************************************************************
c compute c=cosh(a*kd); y=sinh(a*kd)/a; x=sinh(a*kd)*a
c and multiply them by ex=exp(-Real(a*kd)) to supress overflow
c
c called by: compoundMatrix()		in compound.f
c***************************************************************
      IMPLICIT NONE
      complex*16 c,y,x,a
      real*8 kd,r,i,ex
      y = kd*a
      r = dreal(y)
      i = dimag(y)
      ex = dexp(-r)
      y = 0.5d0*dcmplx(dcos(i),dsin(i))
      x = ex*ex*dconjg(y)
      c = y + x
      x = y - x
      y = x/a
      x = x*a
      return
      end


      subroutine haskellMatrix(a)
c***************************************************************
c compute 4x4 P-SV Haskell a for the layer, the layer parameter
c is passed in by common /layer/.
c***************************************************************
      IMPLICIT NONE
      complex*16 a(5,5)
      include 'layer.h'

      Ca = Ca*exb
      Xa = Xa*exb
      Ya = Ya*exb
      Cb = Cb*exa
      Yb = Yb*exa
      Xb = Xb*exa

c p-sv, scaled by exa*exb, see p381/Haskell1964 or EQ 17 of ZR
      a(1,1) = r*(Ca-r1*Cb)
      a(1,2) = r*(r1*Ya-Xb)
      a(1,3) = (Cb-Ca)*r/mu2
      a(1,4) = (Xb-Ya)*r/mu2

      a(2,1) = r*(r1*Yb-Xa)
      a(2,2) = r*(Cb-r1*Ca)
      a(2,3) = (Xa-Yb)*r/mu2
      a(2,4) =-a(1,3)

      a(3,1) = mu2*r*r1*(Ca-Cb)
      a(3,2) = mu2*r*(r1*r1*Ya-Xb)
      a(3,3) = a(2,2)
      a(3,4) =-a(1,2)

      a(4,1) = mu2*r*(r1*r1*Yb-Xa)
      a(4,2) =-a(3,1)
      a(4,3) =-a(2,1)
      a(4,4) = a(1,1)

c sh, the Haskell matrix is not needed. it is replaced by exb
      a(5,5) = exb

      return
      end


      subroutine compoundMatrix(a)
c***************************************************************
c The upper-left 5x5 is the 6x6 compound matrix of the P-SV Haskell matrix,
c	a(ij,kl) = A|_kl^ij, ij=12,13,14,23,24,34,
c after dropping the 3rd row and colummn and replacing the 4th row
c by (2A41, 2A42, 2A44-1,2A45,2A46) (see W&H, P1035).
c The lower-right c 2x2 is the SH part of the Haskell matrix.
c Input: layer parameters passed in by /layer/.
c Output: compound matrix a, scaled by exa*exb for the P-SV and exb for the SH.
c***************************************************************
      IMPLICIT NONE
      real*8 ex
      complex*16 a(7,7)
      complex*16 CaCb, CaXb, CaYb, XaCb, XaXb, YaCb, YaYb, r2, r3
      include 'layer.h'

      call sh_ch(Ca,Ya,Xa,exa,ra,kd)
      call sh_ch(Cb,Yb,Xb,exb,rb,kd)
      CaCb=Ca*Cb
      CaYb=Ca*Yb
      CaXb=Ca*Xb
      XaCb=Xa*Cb
      XaXb=Xa*Xb
      YaCb=Ya*Cb
      YaYb=Ya*Yb
      ex = exa*exb
      r2 = r*r
      r3 = r1*r1

c p-sv, scaled by exa*exb to supress overflow
      a(1,1) = ((one+r3)*CaCb-XaXb-r3*YaYb-two*r1*ex)*r2
      a(1,2) = (XaCb-CaYb)*r/mu2
      a(1,3) = ((one+r1)*(CaCb-ex)-XaXb-r1*YaYb)*r2/mu2
      a(1,4) = (YaCb-CaXb)*r/mu2
      a(1,5) = (two*(CaCb-ex)-XaXb-YaYb)*r2/(mu2*mu2)

      a(2,1) = (r3*YaCb-CaXb)*r*mu2
      a(2,2) = CaCb
      a(2,3) = (r1*YaCb-CaXb)*r
      a(2,4) =-Ya*Xb
      a(2,5) = a(1,4)

      a(3,1) = two*mu2*r2*(r1*r3*YaYb-(CaCb-ex)*(r3+r1)+XaXb)
      a(3,2) = two*r*(r1*CaYb-XaCb)
      a(3,3) = two*(CaCb - a(1,1)) + ex
      a(3,4) =-two*a(2,3)
      a(3,5) =-two*a(1,3)

      a(4,1) = mu2*r*(XaCb-r3*CaYb)
      a(4,2) =-Xa*Yb
      a(4,3) =-a(3,2)/two
      a(4,4) = a(2,2)
      a(4,5) = a(1,2)

      a(5,1) = mu2*mu2*r2*(two*(CaCb-ex)*r3-XaXb-r3*r3*YaYb)
      a(5,2) = a(4,1)
      a(5,3) =-a(3,1)/two
      a(5,4) = a(2,1)
      a(5,5) = a(1,1)

c sh, scaled by exb
      a(6,6) = Cb
      a(6,7) =-two*Yb/mu2
      a(7,6) =-mu2*Xb/two
      a(7,7) = Cb

      return
      end


      subroutine eVector(e)
c***************************************************************
c The first 5 members are E|_12^ij, ij=12,13,23,24,34.
c The last two are the first column of SH E matrix.
c***************************************************************
      IMPLICIT NONE
      include 'layer.h'
      complex*16 e(7)
c For p-sv, compute E|_(12)^(ij), ij=12, 13, 23, 24, 34.
      e(1) = ra*rb-one
      e(2) = mu2*rb*(one-r1)
      e(3) = mu2*(r1-ra*rb)
      e(4) = mu2*ra*(r1-one)
      e(5) = mu2*mu2*(ra*rb-r1*r1)
c sh part
      e(6) =-one
      e(7) = mu2*rb/two
      return
      end


      subroutine initialG(g)
c***************************************************************
c Initialize the g row-vector. The first 5 elements are the
c inverse(E)|_{ij}^{12}, ij=12,13,23,24,34.
c g14 is omitted because g14 = -g23.
c The last two are the 5th row of E^-1.
c***************************************************************
      IMPLICIT NONE
      complex*16 g(7), delta
      include 'layer.h'
c p-sv, see EQ 33 on ZR/p623, constant omitted.
      delta  = r*(one-ra*rb)-one
      g(1) = mu2*(delta-r1)
      g(2) = ra
      g(3) = delta
      g(4) =-rb
      g(5) = (one+delta)/mu2
c sh, use the 5th row of E^-1, see EQ A4 on ZR/p625, 1/2 omitted
      g(6) =-one
      g(7) = two/(rb*mu2)
      return
      end
