c***************************************************************
c  F-K: @(#) kernel.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c Compute dynamic displacement kernels from a point source in a multi-layer medium
c using the Haskell propagator matrix.
c
c	reference:
c		Haskell (1964), BSSA
c		Wang and Herrmann (1980), BSSA
c
c	with some modifications:
c	(1) modify the definition of B and K and put back w^2 so all matrix
c	    elements are either dimensionless or in in shear modular unit.
c	(2) suppress overflow by factoring out e^(Re(ra+rb)*d) in matrix
c	    elements.
c	(3) the cordinate is Z up, R outward and T counterclockwise
c Input:
c	k	wavenumber
c	/model/	velocity model (see model.h)
c Output:
c	u(i,3)	kernels for azimuthal mode i = 0, 1, 2 in the order of Z, R, T.
c		Note that 1/2*pi*mu is omitted.
c Called by
c	main() in fk.f.
c
c Subroutines called are in haskell.f (st_haskell.f in static case) and prop.f.
c***************************************************************
      subroutine kernel(k, u)
      IMPLICIT NONE
      include	'layer.h'
      include	'model.h'
      integer	i, j
      real	k
      complex u(3,3)
      complex*16 rayl, love, dum
      complex*16 a(5,5), b(7,7), c(7,7), e(7), g(7), z(3,5), ss(3,6)
c Explanations:
c a --- 4x4 p-sv Haskell matrix and a(5,5)=exb.
c b --- product of compound matrices from the receiver to the surface.
c c --- compound matrix.
c e --- vector, the first 5 members are E0|_{12}^{ij}, ij=12,13,23,24,34;
c	the last two are the 1st column of the 2x2 SH E0 (or unit vector if the top is elastic).
c g --- vector containing the Rayleigh and Love denominators. It is initialized in the
c	bottom half-space with  (E^-1)|_{ij}^{12}, ij=12,13,23,24,34, and
c	the 1st row of the 2x2 SH E^-1 (or a unit vector if the bottom is vacume).
c z --- z(n,j)=s(i)*X|_{ij}^{12} for p-sv and s(i)*X_5i for sh.

      call initialB(b, e, g)
c propagation - start from the bottom
      do j = mb, 1, -1
	 call layerParameter(k, j)
	 if ( j.EQ.mb .AND. d(j).LT.epsilon ) then
            call initialG(g)
         else if ( j.EQ.1 .AND. d(1).LT.epsilon ) then
	    call eVector(e)
	    exit
         else
	    call compoundMatrix(c)
	    call propagateG(c, g)
         endif
	 if ( j.EQ.src ) then
	    call separatS(ss)
	    call initialZ(ss, g, z)
	 endif
	 if ( j.LT.src ) then
	    if ( j.GE.rcv ) then
	       call haskellMatrix(a)
	       call propagateZ(a, z)
	    else
	       call propagateB(c, b)
	    endif
	 endif
      enddo

c add the top halfspace boundary condition
      e(3) = two*e(3)
      rayl = g(1)*e(1)+g(2)*e(2)+g(3)*e(3)+g(4)*e(4)+g(5)*e(5)
      love = g(6)*e(6)+g(7)*e(7)
      do i=1, 4
	 g(i) = zero
	 do j=1, 5
	    g(i) = g(i) + b(i,j)*e(j)
	 enddo
      enddo
      g(3) = g(3)/two
      g(6) = b(6,6)*e(6) + b(6,7)*e(7)
      do i=1, 3
         dum    = z(i,2)*g(1)+z(i,3)*g(2)-z(i,4)*g(3)
         z(i,2) =-z(i,1)*g(1)+z(i,3)*g(3)+z(i,4)*g(4)
         z(i,1) = dum
         z(i,5) = z(i,5)*g(6)
      enddo

c displacement kernels at the receiver
      dum = k
      if ( stype.EQ.1 ) dum = one
      do i = 1, 3
         u(i,1) = dum*z(i,2)/rayl
	 u(i,2) = dum*z(i,1)/rayl
         u(i,3) = dum*z(i,5)/love
      enddo

      return
      end


      subroutine separatS(ss)
      IMPLICIT NONE
      include	'layer.h'
      include	'model.h'
      integer i, j, ii, jj
      complex*16 temp(4,4), temp_sh, ss(3,6), ra1, rb1, dum

      ii = stype+1
      if (updn.EQ.0) then
         do i=1, ii
            do j=1,6
               ss(i,j) = si(i,j)
            enddo
         enddo
         return
      endif

c down-going (updn=1) or up-going(updn=-1) matrix: E*diag(...)*inv(E), without the 1/2 factor
      ra1 = one/ra
      rb1 = one/rb
      dum = updn*r
      temp(1,1) = one
      temp(1,2) = dum*(rb-r1*ra1)
      temp(1,3) = zero
      temp(1,4) = dum*(ra1-rb)/mu2
      temp(2,1) = dum*(ra-r1*rb1)
      temp(2,2) = one
      temp(2,3) = dum*(rb1-ra)/mu2
      temp(2,4) = zero
      temp(3,1) = zero
      temp(3,2) = dum*(rb-r1*r1*ra1)*mu2
      temp(3,3) = one
      temp(3,4) = dum*(r1*ra1-rb)
      temp(4,1) = dum*(ra-r1*r1*rb1)*mu2
      temp(4,2) = zero
      temp(4,3) = dum*(r1*rb1-ra)
      temp(4,4) = one
      temp_sh   = (updn*two/mu2)*rb1
c
      do i=1, ii
         do j=1,4
	    dum = zero
            do jj = 1,4
               dum = dum + temp(j,jj)*si(i,jj)
            enddo
	    ss(i,j) = dum/two
         enddo
	 ss(i,5) = (si(i,5) + temp_sh*si(i,6))/two
	 ss(i,6) = (si(i,6) + si(i,5)/temp_sh)/two
      enddo
      return
      end
