      SUBROUTINE FTTQ(dt,ts,nfft,n,f)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Compute FUTTERMAN'S Q OPERATOR F(N) in time domain.
C  See A&R P165-170
C  IN:
C  dt	TIME SEPARATION IN THE TIME DOMAIN
C  ts	RATIO T/Q (t-star)
C  nfft	Max. number of points, must be 2**integer
C  OUT:
C  n	number of points
C  f()	the Futterman filter
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
      INTEGER 	i, n, nfft, nfft2
      REAL	dt, ts, w, dw, Pi, a, ep, f(nfft)
      Pi = 3.1415926
      nfft2 = nfft/2
      dw=Pi/(nfft2*dt)
      w = dw
      f(1) = 1.
      f(2) = 0.
      DO i=2,nfft2
         a=w*(ts*ALOG(w)/Pi-nfft2*dt)
	 ep = EXP(-0.5*w*ts)
         f(2*i-1) = ep*cos(a)
	 f(2*i  ) = ep*sin(a)
	 w = w + dw
      ENDDO
      CALL fftr(f,nfft2,-dt)
C cut f(t) so that only values larger enough are returned
      a = 0.
      DO i=1,nfft
	 f(i) = f(i)*dt
	 a = a + f(i)
	 IF (a .GT. 0.02) EXIT
      ENDDO
      f(1) = f(i)
      DO n=2,nfft-i
	 f(n) = f(n+i-1)*dt
	 a = a + f(n)
	 IF (a .GT. 0.98) EXIT
      ENDDO
      RETURN
      END
