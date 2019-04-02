# 1 "bessel.FF"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 360 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "bessel.FF" 2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c bessel.FF:	Compute Bessel function Jn(z) for n=0,1,2
c Reivsion History
c	03/05/1996  Lupei Zhu   Initial coding.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine besselFn(z, aj0, aj1, aj2)
	IMPLICIT NONE
        real z, aj0, aj1, aj2








	aj0 = BesJ0(z)
	aj1 = BesJ1(z)
	aj2 = BesJN(2,z)

# 30 "bessel.FF"
        return
      end

