/********************************************************************************
 * Complex.c	Complex number operations.
 * Revision History:
 * 	Lupei Zhu	06/20/94	Initial coding
*********************************************************************************/
#include <math.h>
#include "Complex.h"

complex cplus(complex a, complex b) {
  a.x += b.x;
  a.y += b.y;
  return(a);
}

complex cmltp(complex a, complex b) {
  complex c;
  c.x = a.x*b.x - a.y*b.y;
  c.y = a.x*b.y + a.y*b.x;
  return(c);
}

complex cngtv(complex a) {
  a.x = -a.x;
  a.y = -a.y;
  return(a);
}

complex cinvs(complex a) {
  complex dmltp(float, complex);
  complex conjg(complex a);
  return(dmltp(1./(a.x*a.x+a.y*a.y), conjg(a)));
}

complex conjg(complex a) {
  a.y = -a.y;
  return(a);
}

complex dmltp(float a, complex b) {
  b.x *= a;
  b.y *= a;
  return(b);
}

complex Csqrt(complex a) {
  double mo, ar;
  double ccabs(complex);
  mo = sqrt(ccabs(a));
  ar = 0.5*atan2(a.y, a.x);
  a.x = mo*cos(ar);
  a.y = mo*sin(ar);
  return(a);
}

complex cmplx(float x, float y) {
  complex a;
  a.x = x;
  a.y = y;
  return(a);
}

complex cphase(complex w) {
  double mo;
  mo = exp(w.x);
  return cmplx(mo*cos(w.y), mo*sin(w.y));
}

double ccabs(complex a) {
  return(sqrt(a.x*a.x+a.y*a.y));
}
