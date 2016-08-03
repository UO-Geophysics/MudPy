c***************************************************************
c  F-K: @(#) prop.f			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c  This file contains the following subroutines for propagating Haskell matrices.
c       propagateG()		propoagate the g vector by the compound matrix.
c       initialZ()		initialize z vector at the source.
c	propagateZ()		propagate the z vector by mulitiplying the Haskell matrix.
c	initialB()		initialize the B matrix and  e/g vectors.
c	propagateB()		compound matrix product.
c***************************************************************


      subroutine propagateG(a, g)
c***************************************************************
c propagate g vector upward using the compound matrix
c	g = g*a
c***************************************************************
      IMPLICIT NONE
      include 'constants.h'
      integer i, j
      complex*16 g(7), a(7,7), temp(7)
c p-sv
      do i = 1, 5
         temp(i) = zero
         do j = 1, 5
            temp(i) = temp(i) + g(j)*a(j,i)
         enddo
      enddo
      do i = 1, 5
         g(i) = temp(i)
      enddo
c sh
      do i = 6, 7
         temp(i) = g(6)*a(6,i) + g(7)*a(7,i)
      enddo
      do i = 6, 7
         g(i) = temp(i)
      enddo
      return
      end


      subroutine initialZ(s, g, z)
c***************************************************************
c initialize the row-vector z at the source z(j)=s(i)*X|_ij^12
c for P-SV and z(j)=s(i)*X(5,i) for SH.
c  input:
c	s(3,6)	---- source coef. for n=0,1,2
c	g(7)	---- g vector used to construct matrix X|_ij^12
c		     |	0   g1  g2 -g3 |
c	 X|_ij^12 =  | -g1  0   g3  g4 | for P-SV.
c		     | -g2 -g3  0   g5 |
c		     |  g3 -g4 -g5  0  |
c	 X(5,i) = ( g6 g7 )	for SH.
c  output:
c	z(3,5)  ---- z vector for n=0,1,2
c***************************************************************
      IMPLICIT NONE
      integer	i
c     real	s(3,6)
      complex*16 g(7), z(3,5), s(3,6)
      do i=1, 3
c for p-sv, see WH p1018
	 z(i,1) =-s(i,2)*g(1)-s(i,3)*g(2)+s(i,4)*g(3)
	 z(i,2) = s(i,1)*g(1)-s(i,3)*g(3)-s(i,4)*g(4)
	 z(i,3) = s(i,1)*g(2)+s(i,2)*g(3)-s(i,4)*g(5)
	 z(i,4) =-s(i,1)*g(3)+s(i,2)*g(4)+s(i,3)*g(5)
c for sh
	 z(i,5) = s(i,5)*g(6)+s(i,6)*g(7)
      enddo
      return
      end


      subroutine propagateZ(a, z)
c***************************************************************
c  apply the Haskell matrix a to the z vector
c	z = z*a
c***************************************************************
      IMPLICIT NONE
      include 'constants.h'
      integer i, j, l
      complex*16 z(3,5), a(5,5), temp(4)
      do i = 1, 3
c p-sv
         do j = 1, 4
            temp(j) = zero
            do l = 1, 4
               temp(j) = temp(j) + z(i,l)*a(l,j)
            enddo
         enddo
         do j = 1, 4
            z(i,j) = temp(j)
         enddo
c sh, due to the exb scaling of the sh haskell matrix.
         z(i,5) = z(i,5)*a(5,5)
      enddo
      return
      end


      subroutine initialB(b, e, g)
c***************************************************************
c Initialize b as an unit matrix; e as an unit vector;
c e = (1 0 0 0 0 1 0) for top halfspace boundary condition;
c g = (0 0 0 0 1 0 1) for free bottom boundary condition.
c***************************************************************
      IMPLICIT NONE
      include 'constants.h'
      integer	i, j
      complex*16 b(7,7), e(7), g(7)
      do i=1, 7
	 do j=1, 7
	    b(i,j) = zero
	 enddo
	 b(i,i) = one
	 e(i) = zero
	 g(i) = zero
      enddo
      e(1) = one
      g(5) = one
      e(6) = one
      g(7) = one
      return
      end


      subroutine propagateB(c, b)
c***************************************************************
c	b = b*c
c***************************************************************
      IMPLICIT NONE
      include 'constants.h'
      integer i, j, l
      complex*16 b(7,7), c(7,7), temp(7)
c p-sv
      do i = 1, 5
         do j = 1, 5
            temp(j) = zero
            do l = 1, 5
               temp(j) = temp(j) + b(i,l)*c(l,j)
            enddo
         enddo
         do j = 1, 5
            b(i,j) = temp(j)
         enddo
      enddo
c sh
      do i = 6, 7
         do j = 6, 7
            temp(j) = b(i,6)*c(6,j) + b(i,7)*c(7,j)
         enddo
         do j = 6, 7
            b(i,j) = temp(j)
         enddo
      enddo
      return
      end
