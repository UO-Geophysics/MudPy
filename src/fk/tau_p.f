      SUBROUTINE findp0(x,p0)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     solve d(tau)/dp = 0 for p0, tau=p*x+eta*z
c     input:
c	 x --- distance
c	 p0 -- the largest possible p
c     output: p0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL ZERO, dtdp0, x
      COMPLEX dtdp, p0, p1, p2
      ZERO = 1.E-7
      p1 = CMPLX(0.,AIMAG(p0))
      DO WHILE ( p1.NE.p0 )
	 p2 = p0
	 p0 = 0.5*(p1+p2)
	 dtdp0 = dtdp(x,p0)
	 IF ( ABS(dtdp0).LT.ZERO .OR. p0.EQ.p1 .OR. p0.EQ.p2 ) RETURN
	 IF( dtdp0 .GT. 0. ) THEN
	    p1 = p0
	    p0 = p2
         END IF
      ENDDO
      RETURN
      END
 
      COMPLEX FUNCTION taup(p,x)
c define function tau(p) = p x + eta h
      IMPLICIT NONE
      INCLUDE'aseries.h'
      INTEGER i
      REAL x
      COMPLEX p, pp
      taup = p*x
      pp = p*p
      DO i = topp, bttm
	 taup=taup+SQRT(vps(1,i)-pp)*ray_len(1,i)
     &		  +SQRT(vps(2,i)-pp)*ray_len(2,i)
      ENDDO
      RETURN
      END
 
      COMPLEX FUNCTION dtdp(x,p)
c define d(tau)/dp
      IMPLICIT NONE
      INCLUDE'aseries.h'
      INTEGER j
      REAL x
      COMPLEX p, pp
      pp = p*p
      dtdp = 0.0
      DO j = topp, bttm
	 dtdp=dtdp-ray_len(1,j)/SQRT(vps(1,j)-pp)
     &		  -ray_len(2,j)/SQRT(vps(2,j)-pp)
      ENDDO
      dtdp = x + p*dtdp
      RETURN
      END
 
      SUBROUTINE dtdp23(p, dt2dp2, dt3dp3)
c calculate the second and 3rd derivatives of tau(p)
      IMPLICIT NONE
      INCLUDE'aseries.h'
      COMPLEX e1, e2, p, pp, dt2dp2, dt3dp3
      INTEGER j
      pp = p*p
      dt2dp2 = 0.
      dt3dp3 = 0.
      DO j = topp, bttm
	 e1 = vps(1,j)-pp
	 e2 = vps(2,j)-pp
	 dt2dp2=dt2dp2-ray_len(1,j)*vps(1,j)/e1/SQRT(e1)
     &	              -ray_len(2,j)*vps(2,j)/e2/SQRT(e2)
	 dt3dp3=dt3dp3-3*p*ray_len(1,j)*vps(1,j)/e1/e1/SQRT(e1)
     &	              -3*p*ray_len(2,j)*vps(2,j)/e2/e2/SQRT(e2)
      ENDDO
      RETURN
      END
      
      COMPLEX FUNCTION time2(x,p,dpdt)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c find the p(tau) contour which makes Im(tau) = 0. Return the tau.
c INPUT:
c	x: distance
c	p: complex ray parameter on the contour
c	dpdt: dp/dt at this point p
c OUT:
c	p: new point on the contour with time advanced by ~dt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      IMPLICIT NONE
      REAL ti, ti1, ti2, x, ZERO, delpr
      COMPLEX taup, p, p1, p2, dpdt
c
      ZERO = 1.E-7		! precision
c
      p1 = p+CMPLX(0.,AIMAG(dpdt))	! move the point up above the contour
      time2 = taup(p1,x)
      ti1 = AIMAG(time2)		! this should be > 0.
      IF (ti1 .LT. 0) THEN
         WRITE(0,*)'Contour search failed, p1 is not on the left, p=',p1
	 CALL EXIT(1)
      END IF
      delpr = 2*REAL(dpdt)
      p = p1
      ti = ti1
c first make sure p and p1 are on the opposite sides of the contour
      DO WHILE ( ti*ti1 .GT. 0. )
c        WRITE(0,*)REAL(p1),ti1,REAL(p),ti
	 p1 = p
	 ti1 = ti
	 p = p+delpr
	 time2 = taup(p,x)
	 ti = AIMAG(time2)
      ENDDO
c begin to find the point on the contour bracketed by [p1 p]
      DO WHILE ( ABS(ti).GT.ZERO )
         p2 = p
         ti2 = ti
	 p = p1 + REAL(p-p1)*ti1/(ti1-ti)
	 time2 = taup(p,x)
	 IF ( p.EQ.p1 .OR. p.EQ.p2 ) RETURN
	 ti = AIMAG(time2)
	 IF( ti1*ti.GT.0. )THEN
	    p1 = p2
	    ti1 = ti2
	 ENDIF
      ENDDO
      RETURN
      END
 
      REAL FUNCTION tstar(p)
c tstar = traveltime/Q, not working for head-wave
      IMPLICIT NONE
      INCLUDE'aseries.h'
      INTEGER i
      REAL pp
      COMPLEX p
      pp = REAL(p)*REAL(p)
      tstar = 0.
      DO i = topp, bttm
         IF (ray_len(1,i) .GT. 0.) THEN
            tstar=tstar+ray_len(1,i)*vps(1,i)/SQRT(vps(1,i)-pp)/q(1,i)
	 ENDIF
         IF (ray_len(2,i) .GT. 0.) THEN
            tstar=tstar+ray_len(2,i)*vps(2,i)/SQRT(vps(2,i)-pp)/q(2,i)
         ENDIF
      ENDDO
      RETURN
      END
