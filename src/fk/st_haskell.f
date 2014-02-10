c***************************************************************
c  F-K: @(#) st_haskell.f                  1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c  This file contains the following subroutines for computing static Haskell propagator
c  matrix and other related matrices/vectors for the static case.
c       layerParameter()        compute some parameters for each layer.
c       haskell()               compute the 6x6 Haskell matrix.
c	compoundMatrix()	compute the compound matrix.
c	eVector()		compute E_{12}^{ij}.
c       initialG()              initialize g matrix in the bottom half space.
c***************************************************************

      subroutine layerParameter(k, lay)
      IMPLICIT NONE
      include   'model.h'
      include   'layer.h'
      integer	lay
      real	k
      mu2 = two*mu(lay)
      kd  = k*d(lay)
      y   = xi(lay)
      y1  = one-y
      return
      end


      subroutine haskellMatrix(a)
c***************************************************************
c p-sv haskell matrix, scaled by exp(-2kd) to suppress overflow.
c***************************************************************
      IMPLICIT NONE
      real*8 C, S, x, y2
      complex*16 a(5,5)
      include 'layer.h'

      C = exb*(one+exa)/two
      S = exb*(one-exa)/two
      y2 = one+y
      x = kd*y1

c for p-sv
      a(1,1) = C+S*x
      a(1,2) = C*x+S*y
      a(1,3) =-S*x/mu2
      a(1,4) =-(S*y2+C*x)/mu2

      a(2,1) = S*y-C*x
      a(2,2) = C-S*x
      a(2,3) = (C*x-S*y2)/mu2
      a(2,4) =-a(1,3)

      a(3,1) = mu2*S*x
      a(3,2) = mu2*(C*x-S*y1)
      a(3,3) = a(2,2)
      a(3,4) =-a(1,2)

      a(4,1) =-mu2*(C*x+S*y1)
      a(4,2) =-a(3,1)
      a(4,3) =-a(2,1)
      a(4,4) = a(1,1)

c for sh, not needed other than the exb factor
      a(5,5) = exb

      return
      end


      subroutine compoundMatrix(a)
c***************************************************************
c p-sv compound matrix scaled by exb*exb; sh haskell matrix scaled by exb
c***************************************************************
      IMPLICIT NONE
      real*8 C, S, C2, S2, CS, x, x2, y2
      complex*16 a(7,7)
      include 'layer.h'

      exb = dexp(-kd)
      exa = exb*exb
      C = (one+exa)/two
      S = (one-exa)/two
      y2 = one+y
      x = kd*y1
      C2 = C*C
      S2 = S*S
      CS = C*S
      x2 = x*x

c for p-sv
      a(1,1) = (x2+one)*exa+y1*y2*S2
      a(1,2) = (x*exa-y2*CS)/mu2
      a(1,3) = (x2*exa-y*y2*S2)/mu2
      a(1,4) = (x*exa+y2*CS)/mu2
      a(1,5) = (x2*exa-y2*y2*S2)/(mu2*mu2)

      a(2,1) = (x*exa-y1*CS)*mu2
      a(2,2) = C2
      a(2,3) = x*exa+y*CS
      a(2,4) =-S2
      a(2,5) = a(1,4)
      
      a(3,1) =-two*mu2*(x2*exa+y*y1*S2)
      a(3,2) = two*(y*CS-x*exa)
      a(3,3) = two*((one-x2)*exa+y*y*S2)-exa
      a(3,4) =-two*a(2,3)
      a(3,5) =-two*a(1,3)
      
      a(4,1) = mu2*(x*exa+y1*CS)
      a(4,2) = a(2,4)
      a(4,3) =-a(3,2)/two
      a(4,4) = a(2,2)
      a(4,5) = a(1,2)

      a(5,1) = mu2*mu2*(x2*exa-y1*y1*S2)
      a(5,2) = a(4,1)
      a(5,3) =-a(3,1)/two
      a(5,4) = a(2,1)
      a(5,5) = a(1,1)

c for sh
      a(6,6) = C
      a(6,7) =-two*S/mu2
      a(7,6) =-mu2*S/two
      a(7,7) = C

      return
      end


      subroutine eVector(e)
c***************************************************************
c e(1:5) = E_12^ij; e(6:7) is the 1st column of E_sh.
c***************************************************************
      IMPLICIT NONE
      include 'layer.h'
      complex*16 e(7)
c for p-sv, they are E|_(12)^(ij), ij=12, 13, 23, 24, 34.
      e(1) =-(one+y)
      e(2) = mu2
      e(3) = mu2*y
      e(4) =-mu2
      e(5) = mu2*mu2*y1
c for sh
      e(6) =-one
      e(7) = mu2/two
      return
      end


      subroutine initialG(g)
c**************************************************************
c Initialize the g row-vector. The first 5 elements are the
c inverse(E)|_{ij}^{12}, ij=12,13,23,24,34.
c The last two are the 5th row of E^-1, see Eq. 34 and Eq. A9 in ZR.
c**************************************************************
      IMPLICIT NONE
      include 'layer.h'
      complex*16 g(7)
c p-sv
      g(1) =-mu2*y1
      g(2) = one
      g(3) = y
      g(4) =-one
      g(5) = (one+y)/mu2
c sh
      g(6) = g(4)
      g(7) = two/mu2
      return
      end
