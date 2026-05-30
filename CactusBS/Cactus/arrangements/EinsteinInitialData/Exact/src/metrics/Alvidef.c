#ifdef Append
#define Alvireg1 Alvireg1_
#define Alvireg2 Alvireg2_
#define Alvireg3 Alvireg3_
#define Alvireg4 Alvireg4_
#endif

#ifdef Append2
#define Alvireg1 Alvireg1__
#define Alvireg2 Alvireg2__
#define Alvireg3 Alvireg3__
#define Alvireg4 Alvireg4__
#endif

#ifdef Uppercase
#define Alvireg1 ALVIREG1
#define Alvireg2 ALVIREG2
#define Alvireg3 ALVIREG3
#define Alvireg4 ALVIREG4
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Power(x,y) pow((double) (x), (double) (y))
#define Sqrt(x)    pow((double) (x), 0.5)
#define Log(x)     log((double) (x)))
#define pow2(x)    ((x)*(x))
#define pow2inv(x) (1.0/((x)*(x)))

void Alvireg1(double *x, double *y, double *z, double *zm1, double *zm2, double *zr12, double *results) {

double g00,g0x,g0y,g0z,gxx,gxy,gxz,gyy,gyz,gzz;

double m1;
double m2;
double b;
double xp;
double yp;
double zp;

m1 = *zm1;
m2 = *zm2;
b=*zr12;
xp = *x;
yp = *y;
zp = *z;

/* finding the 4-metric */

g00 =  Power(1 - (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* ((m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2* Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)) ;
g0x =        (-2*Power(m2,2)*Sqrt((m1 + m2)/b)*(1 - (m2*(1 + m2/(2.*(m1 +
       m2))))/b)* (1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)* (1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2))*yp*Power(zp,2)*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,5)
       + (1 - (m2*(1 + m2/(2.*(m1 + m2))))/b)* (1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2))* ((-2*m2*Sqrt((m1 +
       m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       ((-2*Power(m1,2))/ (Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)) +
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3) )) +
       (m2*Sqrt((m1 + m2)/b)*(1 - (m2*(1 + m2/(2.*(m1 + m2))))/b)* (1 +
       m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b)*yp*
       ((m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,2) ;
g0y =        -((m2*(1 - (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((-2*m2*Sqrt((m1
       + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2)) + (1 - (m2*(1 + m2/(2.*(m1 +
       m2))))/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)* ((2*m2*Sqrt((m1
       + m2)/b)*(-(Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2)) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/
       (2.*Power(b,2)))*Power(1 + m1/(2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)
       )) + (1 - (m2*(1 + m2/(2.*(m1 + m2))))/b)* (-(Sqrt((m1 +
       m2)/b)*(m2/(m1 + m2) + (m2*(3 + Power(m2,2)/(2.*Power(m1 +
       m2,2)) + m1/(m1 + m2)))/b)) + (m2*Sqrt((m1 + m2)/b)*((1 + m1/(m1
       + m2))*(-((b*m2)/(m1 + m2)) + xp) - (3*Power(-((b*m2)/(m1 + m2))
       + xp,2) - 3*Power(yp,2) - Power(zp,2))/(2.*b)))/Power(b,2))*
       ((m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)) ;
g0z =        (2*m2*Sqrt((m1 + m2)/b)*(1 - (m2*(1 + m2/(2.*(m1 + m2))))/b)* (1
       + (m2*(1 + m2/(2.*(m1 + m2))))/b)* Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*yp*zp* Power(1 -
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + (m2*(1 - (m2*(1 + m2/(2.*(m1 + m2))))/b)*zp* ((-2*m2*Sqrt((m1
       + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + (m2*Sqrt((m1 + m2)/b)*(1 -
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp*zp* ((m2*(-(Power(1 + (m2*(1
       + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,3) ;
gxx =         (4*m1*Power(m2,2)*Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2)*((1 + m2/b)*(-((b*m2)/(m1 + m2))
       + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (4*Power(m2,3)*(m1 + m2)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*
       (1 + m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b)* (1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2))*Power(yp,2)*
       Power(zp,2)*Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,8)
       + (Power(m2,2)*Power(zp,2)*Power(1 + m1/(2.*Sqrt(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2)* (1 +
       Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)))))
       /Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)))))/ Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) +
       (2*m2*Sqrt((m1 + m2)/b)*(1 + m1/(m1 + m2) - (3*(-((b*m2)/(m1 +
       m2)) + xp))/b)* (1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))*yp* ((-2*m2*Sqrt((m1 + m2)/b)*(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) +
       xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/ (2.*Sqrt(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + (Power(m2,2)*(m1 + m2)*Power(1 +
       m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b,2)*Power(yp,2)*
       ((m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,5) ;
gxy =         (2*m1*Power(m2,2)*Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       (1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))*yp*Power(zp,2)* (-(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* (1 + Power(m1,2)/
       (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*m1*m2*Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* (1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2))*yp* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/
       (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,3)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*m1*Power(m2,3)*(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))* yp*Power(zp,2)*((1 + m2/b)*(-((b*m2)/(m1 +
       m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* (1 + Power(m1,2)/
       (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,7)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*Power(m2,2)*Sqrt((m1 + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b)* (1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))*yp*Power(zp,2)* (-(Sqrt((m1 + m2)/b)*(m2/(m1 +
       m2) + (m2*(3 + Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 +
       m2)))/b)) + (m2*Sqrt((m1 + m2)/b)* ((1 + m1/(m1 +
       m2))*(-((b*m2)/(m1 + m2)) + xp) - (3*Power(-((b*m2)/(m1 + m2)) +
       xp,2) - 3*Power(yp,2) - Power(zp,2))/(2.*b)))/Power(b,2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,5)
       - (m2*(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))*yp* Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,2)
       - (Power(m2,2)*Sqrt((m1 + m2)/b)* (1 + m1/(m1 + m2) -
       (3*(-((b*m2)/(m1 + m2)) + xp))/b)*Power(yp,2)* ((-2*m2*Sqrt((m1
       + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,4) + (1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2))* (-(Sqrt((m1 + m2)/b)*(m2/(m1 + m2) +
       (m2*(3 + Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 + m2)))/b))
       + (m2*Sqrt((m1 + m2)/b)*((1 + m1/(m1 + m2))*(-((b*m2)/(m1 + m2))
       + xp) - (3*Power(-((b*m2)/(m1 + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((-2*m2*Sqrt((m1 + m2)/b)*(1
       + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       ((-2*Power(m1,2))/ (Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)) +
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3) )) +
       (m2*Sqrt((m1 + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)* (1 +
       m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b)*yp*
       ((2*m2*Sqrt((m1 + m2)/b)*(-(Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2)) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + (m2*Sqrt((m1 + m2)/b)*(1 +
       m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b)*yp* (-(Sqrt((m1
       + m2)/b)*(m2/(m1 + m2) + (m2*(3 + Power(m2,2)/(2.*Power(m1 +
       m2,2)) + m1/(m1 + m2)))/b)) + (m2*Sqrt((m1 + m2)/b)* ((1 +
       m1/(m1 + m2))*(-((b*m2)/(m1 + m2)) + xp) - (3*Power(-((b*m2)/(m1
       + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((m2*(-(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,2);
gxz =         (-2*m1*m2*Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),3)*zp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,3)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*m1*Power(m2,3)*(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))* Power(zp,3)*((1 + m2/b)*(-((b*m2)/(m1 + m2)) +
       xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,7)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*Power(m2,2)*(m1 + m2)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*
       (1 + m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b)* Power(1 +
       m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(yp,2)*zp*Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,6)
       - (2*Power(m2,3)*(m1 + m2)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*
       (1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))*Power(yp,2)* Power(zp,3)*Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,9) - (m2*(1 +
       m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2))*zp* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2)* (1 +
       Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)))))
       /Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,2)
       + (m2*(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2))*zp* Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,2)
       + (Power(m2,2)*Sqrt((m1 + m2)/b)* (1 + m1/(m1 + m2) -
       (3*(-((b*m2)/(m1 + m2)) + xp))/b)*yp*zp* ((-2*m2*Sqrt((m1 +
       m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,4) + (m2*Sqrt((m1 + m2)/b)*(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2))*yp*zp*
       ((-2*m2*Sqrt((m1 + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/ (2.*Sqrt(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,3) + (Power(m2,2)*(m1 + m2)*(1 +
       m1/(m1 + m2) - (3*(-((b*m2)/(m1 + m2)) + xp))/b)*Power(yp,2)*zp*
       ((m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,6) ;
gyy =        (4*m1*Power(m2,2)*Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2)* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp)
       - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) + Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       ((-2*m1)/ (Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)) -
       (2*m1*Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)))))/ Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) +
       (Power(m2,2)*Power(yp,2)*Power(1 + m1/(2.*Sqrt(Power(1 + (m2*(1
       + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       - (2*m2*yp*(-(Sqrt((m1 + m2)/b)* (m2/(m1 + m2) + (m2*(3 +
       Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 + m2)))/b)) +
       (m2*Sqrt((m1 + m2)/b)*((1 + m1/(m1 + m2))* (-((b*m2)/(m1 + m2))
       + xp) - (3*Power(-((b*m2)/(m1 + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((-2*m2*Sqrt((m1 + m2)/b)*(1
       + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + 2*(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b)* (-(Sqrt((m1 + m2)/b)*(m2/(m1 + m2) + (m2*(3 +
       Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 + m2)))/b)) +
       (m2*Sqrt((m1 + m2)/b)*((1 + m1/(m1 + m2))*(-((b*m2)/(m1 + m2)) +
       xp) - (3*Power(-((b*m2)/(m1 + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((2*m2*Sqrt((m1 +
       m2)/b)*(-(Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2)) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/
       (2.*Power(b,2)))*Power(1 + m1/(2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)
       )) + Power(-(Sqrt((m1 + m2)/b)* (m2/(m1 + m2) + (m2*(3 +
       Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 + m2)))/b)) +
       (m2*Sqrt((m1 + m2)/b)*((1 + m1/(m1 + m2))*(-((b*m2)/(m1 + m2)) +
       xp) - (3*Power(-((b*m2)/(m1 + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2),2)* ((m2*(-(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)) ;
gyz =        (-2*m1*m2*Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(1 +
       m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*yp*zp*
       (-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) -
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* (1 + Power(m1,2)/
       (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,3)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*m1*Power(m2,2)*Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*yp*zp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*m1*Power(m2,2)*Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*yp*zp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp)
       - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*m2*Sqrt((m1 + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*yp*zp* (-(Sqrt((m1 + m2)/b)*(m2/(m1 + m2) +
       (m2*(3 + Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 + m2)))/b))
       + (m2*Sqrt((m1 + m2)/b)* ((1 + m1/(m1 + m2))*(-((b*m2)/(m1 +
       m2)) + xp) - (3*Power(-((b*m2)/(m1 + m2)) + xp,2) -
       3*Power(yp,2) - Power(zp,2))/(2.*b)))/Power(b,2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) -
       (Power(m2,2)*yp*zp*Power(1 + m1/(2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       - (Power(m2,2)*Sqrt((m1 + m2)/b)*Power(yp,2)*zp*
       ((-2*m2*Sqrt((m1 + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/ (2.*Sqrt(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,5) + (m2*zp*(-(Sqrt((m1 +
       m2)/b)*(m2/(m1 + m2) + (m2*(3 + Power(m2,2)/(2.*Power(m1 +
       m2,2)) + m1/(m1 + m2)))/b)) + (m2*Sqrt((m1 + m2)/b)* ((1 +
       m1/(m1 + m2))*(-((b*m2)/(m1 + m2)) + xp) - (3*Power(-((b*m2)/(m1
       + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((-2*m2*Sqrt((m1 + m2)/b)*(1
       + (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + (m2*Sqrt((m1 + m2)/b)*(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp*zp* ((2*m2*Sqrt((m1 +
       m2)/b)*(-(Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2)) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,3) + (m2*Sqrt((m1 +
       m2)/b)*yp*zp*(-(Sqrt((m1 + m2)/b)* (m2/(m1 + m2) + (m2*(3 +
       Power(m2,2)/(2.*Power(m1 + m2,2)) + m1/(m1 + m2)))/b)) +
       (m2*Sqrt((m1 + m2)/b)*((1 + m1/(m1 + m2))* (-((b*m2)/(m1 + m2))
       + xp) - (3*Power(-((b*m2)/(m1 + m2)) + xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((m2*(-(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,3) ;
gzz =        (-4*m1*Power(m2,2)*Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2)* ((1 + m2/b)*(-((b*m2)/(m1 + m2))
       + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (4*Power(m2,2)*(m1 + m2)*(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b)*
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(yp,2)*Power(zp,2)*Power(1 - m1/(2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,7)
       + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2)*(1 +
       Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)))))/
       Power(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) +
       (Power(m2,2)*Power(zp,2)*Power(1 + m1/(2.*Sqrt(Power(1 + (m2*(1
       + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m1)/ (Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m1*Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m1,2)/ (4.*(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       + (2*Power(m2,2)*Sqrt((m1 + m2)/b)*yp*Power(zp,2)*
       ((-2*m2*Sqrt((m1 + m2)/b)*(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b)*yp* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m1/ (2.*Sqrt(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b)*yp* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m2*(-(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,5) + (Power(m2,2)*(m1 +
       m2)*Power(yp,2)*Power(zp,2)* ((m2*(-(Power(1 + (m2*(1 +
       m2/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m2*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m2*(1 + m2/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m2/b - (m2*(-((b*m2)/(m1
       + m2)) + xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m1/ (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 +
       m2)) + xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m2/b)*(-((b*m2)/(m1 + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 +
       m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m2*(-(Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m2/b
       - (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m1,2))/ (Power(1 +
       (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m2/b -
       (m2*(-((b*m2)/(m1 + m2)) + xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m2/b)*(-((b*m2)/(m1 + m2)) + xp) -
       (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m1/
       (2.*Sqrt(Power(1 + (m2*(1 + m2/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m2/b - (m2*(-((b*m2)/(m1 + m2)) +
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m2/b)*(-((b*m2)/(m1
       + m2)) + xp) - (m2*(Power(-((b*m2)/(m1 + m2)) + xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,7) ;

*results=g00;
*(results+1)=g0x;
*(results+2)=g0y;
*(results+3)=g0z;
*(results+4)=gxx;
*(results+5)=gxy;
*(results+6)=gxz;
*(results+7)=gyy;
*(results+8)=gyz;
*(results+9)=gzz; 

}

void Alvireg2(double *x, double *y, double *z, double *zm1, double *zm2, double *zr12, double *results ) {

double g00,g0x,g0y,g0z,gxx,gxy,gxz,gyy,gyz,gzz;

double m1;
double m2;
double b;
double xp;
double yp;
double zp;

m1 = *zm1;
m2 = *zm2;
b=*zr12;
xp = *x;
yp = *y;
zp = *z;

/* finding the 4-metric */

g00 =        Power(1 - (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* ((m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2* Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)) ;
g0x =        (-2*Power(m1,2)*Sqrt((m1 + m2)/b)*(1 - (m1*(1 + m1/(2.*(m1 +
       m2))))/b)* (1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)* (1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2))*yp*Power(zp,2)*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,5)
       - (1 - (m1*(1 + m1/(2.*(m1 + m2))))/b)* (1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2))* ((2*m1*Sqrt((m1 +
       m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       ((-2*Power(m2,2))/ (Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)) +
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3) )) +
       (m1*Sqrt((m1 + m2)/b)*(1 - (m1*(1 + m1/(2.*(m1 + m2))))/b)* (1 +
       m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b)*yp*
       ((m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,2) ;
g0y =        -((m1*(1 - (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((2*m1*Sqrt((m1 +
       m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2)) - (1 - (m1*(1 + m1/(2.*(m1 +
       m2))))/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)* ((2*m1*Sqrt((m1
       + m2)/b)*(-(Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2)) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/
       (2.*Power(b,2)))*Power(1 + m2/(2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)
       )) - (1 - (m1*(1 + m1/(2.*(m1 + m2))))/b)* (-(Sqrt((m1 +
       m2)/b)*(m1/(m1 + m2) + (m1*(3 + Power(m1,2)/(2.*Power(m1 +
       m2,2)) + m2/(m1 + m2)))/b)) + (m1*Sqrt((m1 + m2)/b)*((1 + m2/(m1
       + m2))*(-((b*m1)/(m1 + m2)) - xp) - (3*Power(-((b*m1)/(m1 + m2))
       - xp,2) - 3*Power(yp,2) - Power(zp,2))/(2.*b)))/Power(b,2))*
       ((m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)) ;
g0z =        (-2*m1*Sqrt((m1 + m2)/b)*(1 - (m1*(1 + m1/(2.*(m1 + m2))))/b)*
       (1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)* Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*yp*zp* Power(1 -
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + (m1*(1 - (m1*(1 + m1/(2.*(m1 + m2))))/b)*zp* ((2*m1*Sqrt((m1 +
       m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) - (m1*Sqrt((m1 + m2)/b)*(1 -
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp*zp* ((m1*(-(Power(1 + (m1*(1
       + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,3) ;
gxx =        (4*Power(m1,2)*m2*Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2)*((1 + m1/b)*(-((b*m1)/(m1 + m2))
       - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (4*Power(m1,3)*(m1 + m2)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*
       (1 + m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b)* (1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2))*Power(yp,2)*
       Power(zp,2)*Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,8)
       + (Power(m1,2)*Power(zp,2)*Power(1 + m2/(2.*Sqrt(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2)* (1 +
       Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)))))
       /Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)))))/ Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) -
       (2*m1*Sqrt((m1 + m2)/b)*(1 + m2/(m1 + m2) - (3*(-((b*m1)/(m1 +
       m2)) - xp))/b)* (1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))*yp* ((2*m1*Sqrt((m1 + m2)/b)*(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b)*yp* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) -
       xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m2/ (2.*Sqrt(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + (Power(m1,2)*(m1 + m2)*Power(1 +
       m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b,2)*Power(yp,2)*
       ((m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,5) ;

gxy =        (-2*Power(m1,2)*m2*Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       (1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))*yp*Power(zp,2)* (-(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* (1 + Power(m2,2)/
       (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*m1*m2*Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* (1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2))*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/
       (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,3)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*Power(m1,3)*m2*(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))* yp*Power(zp,2)*((1 + m1/b)*(-((b*m1)/(m1 +
       m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* (1 + Power(m2,2)/
       (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,7)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*Power(m1,2)*Sqrt((m1 + m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b)* (1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))*yp*Power(zp,2)* (-(Sqrt((m1 + m2)/b)*(m1/(m1 +
       m2) + (m1*(3 + Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 +
       m2)))/b)) + (m1*Sqrt((m1 + m2)/b)* ((1 + m2/(m1 +
       m2))*(-((b*m1)/(m1 + m2)) - xp) - (3*Power(-((b*m1)/(m1 + m2)) -
       xp,2) - 3*Power(yp,2) - Power(zp,2))/(2.*b)))/Power(b,2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,5)
       + (m1*(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))*yp* Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,2)
       - (Power(m1,2)*Sqrt((m1 + m2)/b)* (1 + m2/(m1 + m2) -
       (3*(-((b*m1)/(m1 + m2)) - xp))/b)*Power(yp,2)* ((2*m1*Sqrt((m1 +
       m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,4) + (1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2))* (-(Sqrt((m1 + m2)/b)*(m1/(m1 + m2) +
       (m1*(3 + Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 + m2)))/b))
       + (m1*Sqrt((m1 + m2)/b)*((1 + m2/(m1 + m2))*(-((b*m1)/(m1 + m2))
       - xp) - (3*Power(-((b*m1)/(m1 + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((2*m1*Sqrt((m1 + m2)/b)*(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 + m1/b)*(-((b*m1)/(m1 +
       m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       ((-2*Power(m2,2))/ (Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)) +
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3) )) -
       (m1*Sqrt((m1 + m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)* (1 +
       m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b)*yp*
       ((2*m1*Sqrt((m1 + m2)/b)*(-(Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2)) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) - (m1*Sqrt((m1 + m2)/b)*(1 +
       m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b)*yp* (-(Sqrt((m1
       + m2)/b)*(m1/(m1 + m2) + (m1*(3 + Power(m1,2)/(2.*Power(m1 +
       m2,2)) + m2/(m1 + m2)))/b)) + (m1*Sqrt((m1 + m2)/b)* ((1 +
       m2/(m1 + m2))*(-((b*m1)/(m1 + m2)) - xp) - (3*Power(-((b*m1)/(m1
       + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((m1*(-(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,2) ;

gxz =        (2*m1*m2*Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),3)*zp* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,3)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*Power(m1,3)*m2*(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))* Power(zp,3)*((1 + m1/b)*(-((b*m1)/(m1 + m2)) -
       xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,7)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*Power(m1,2)*(m1 + m2)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*
       (1 + m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b)* Power(1 +
       m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(yp,2)*zp*Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,6)
       + (2*Power(m1,3)*(m1 + m2)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*
       (1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))*Power(yp,2)* Power(zp,3)*Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,9) + (m1*(1 +
       m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2))*zp* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2)* (1 +
       Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)))))
       /Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,2)
       - (m1*(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2))*zp* Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,2)
       + (Power(m1,2)*Sqrt((m1 + m2)/b)* (1 + m2/(m1 + m2) -
       (3*(-((b*m1)/(m1 + m2)) - xp))/b)*yp*zp* ((2*m1*Sqrt((m1 +
       m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,4) + (m1*Sqrt((m1 + m2)/b)*(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2))*yp*zp*
       ((2*m1*Sqrt((m1 + m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b)*yp* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m2/ (2.*Sqrt(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,3) - (Power(m1,2)*(m1 + m2)*(1 +
       m2/(m1 + m2) - (3*(-((b*m1)/(m1 + m2)) - xp))/b)*Power(yp,2)*zp*
       ((m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) - Power(1 -
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,6) ;
gyy =        (4*Power(m1,2)*m2*Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2)* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp)
       - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) + Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       ((-2*m2)/ (Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)) -
       (2*m2*Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)))))/ Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) +
       (Power(m1,2)*Power(yp,2)*Power(1 + m2/(2.*Sqrt(Power(1 + (m1*(1
       + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       + (2*m1*yp*(-(Sqrt((m1 + m2)/b)* (m1/(m1 + m2) + (m1*(3 +
       Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 + m2)))/b)) +
       (m1*Sqrt((m1 + m2)/b)*((1 + m2/(m1 + m2))* (-((b*m1)/(m1 + m2))
       - xp) - (3*Power(-((b*m1)/(m1 + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((2*m1*Sqrt((m1 + m2)/b)*(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 + m1/b)*(-((b*m1)/(m1 +
       m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + 2*(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b)* (-(Sqrt((m1 + m2)/b)*(m1/(m1 + m2) + (m1*(3 +
       Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 + m2)))/b)) +
       (m1*Sqrt((m1 + m2)/b)*((1 + m2/(m1 + m2))*(-((b*m1)/(m1 + m2)) -
       xp) - (3*Power(-((b*m1)/(m1 + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((2*m1*Sqrt((m1 +
       m2)/b)*(-(Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2)) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/
       (2.*Power(b,2)))*Power(1 + m2/(2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)
       )) + Power(-(Sqrt((m1 + m2)/b)* (m1/(m1 + m2) + (m1*(3 +
       Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 + m2)))/b)) +
       (m1*Sqrt((m1 + m2)/b)*((1 + m2/(m1 + m2))*(-((b*m1)/(m1 + m2)) -
       xp) - (3*Power(-((b*m1)/(m1 + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2),2)* ((m1*(-(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)) ;
gyz =         (-2*m1*m2*Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(1 +
       m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*yp*zp*
       (-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) -
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* (1 + Power(m2,2)/
       (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,3)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) -
       (2*Power(m1,2)*m2*Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*yp*zp* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*Power(m1,2)*m2*Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*yp*zp* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp)
       - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (2*m1*Sqrt((m1 + m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*yp*zp* (-(Sqrt((m1 + m2)/b)*(m1/(m1 + m2) +
       (m1*(3 + Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 + m2)))/b))
       + (m1*Sqrt((m1 + m2)/b)* ((1 + m2/(m1 + m2))*(-((b*m1)/(m1 +
       m2)) - xp) - (3*Power(-((b*m1)/(m1 + m2)) - xp,2) -
       3*Power(yp,2) - Power(zp,2))/(2.*b)))/Power(b,2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3) -
       (Power(m1,2)*yp*zp*Power(1 + m2/(2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       + (Power(m1,2)*Sqrt((m1 + m2)/b)*Power(yp,2)*zp* ((2*m1*Sqrt((m1
       + m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,5) - (m1*zp*(-(Sqrt((m1 +
       m2)/b)*(m1/(m1 + m2) + (m1*(3 + Power(m1,2)/(2.*Power(m1 +
       m2,2)) + m2/(m1 + m2)))/b)) + (m1*Sqrt((m1 + m2)/b)* ((1 +
       m2/(m1 + m2))*(-((b*m1)/(m1 + m2)) - xp) - (3*Power(-((b*m1)/(m1
       + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((2*m1*Sqrt((m1 + m2)/b)*(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* ((1 + m1/b)*(-((b*m1)/(m1 +
       m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,2) + (m1*Sqrt((m1 + m2)/b)*(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp*zp* ((2*m1*Sqrt((m1 +
       m2)/b)*(-(Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2)) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 +
       m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* ((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,3) + (m1*Sqrt((m1 +
       m2)/b)*yp*zp*(-(Sqrt((m1 + m2)/b)* (m1/(m1 + m2) + (m1*(3 +
       Power(m1,2)/(2.*Power(m1 + m2,2)) + m2/(m1 + m2)))/b)) +
       (m1*Sqrt((m1 + m2)/b)*((1 + m2/(m1 + m2))* (-((b*m1)/(m1 + m2))
       - xp) - (3*Power(-((b*m1)/(m1 + m2)) - xp,2) - 3*Power(yp,2) -
       Power(zp,2))/(2.*b)))/Power(b,2))* ((m1*(-(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,3) ;

gzz =        (-4*Power(m1,2)*m2*Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2)* ((1 + m1/b)*(-((b*m1)/(m1 + m2))
       - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/ (2.*Power(b,2)))*(-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)*Power(zp,2) + 2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)*Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/
       (Power(b,5)*Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2),1.5)) +
       (4*Power(m1,2)*(m1 + m2)*(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b)*
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(yp,2)*Power(zp,2)*Power(1 - m2/(2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,7)
       + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2)*(1 +
       Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)))))/
       Power(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) +
       Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*
       Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) +
       (Power(m1,2)*Power(zp,2)*Power(1 + m2/(2.*Sqrt(Power(1 + (m1*(1
       + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*m2)/ (Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) - (2*m2*Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2)*
       (1 + Power(m2,2)/ (4.*(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)* Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))))) /Power(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2),1.5) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)))/Power(b,3)) )/Power(b,4)
       - (2*Power(m1,2)*Sqrt((m1 + m2)/b)*yp*Power(zp,2)*
       ((2*m1*Sqrt((m1 + m2)/b)*(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b)*yp* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* Power(1 - m2/ (2.*Sqrt(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2) + Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + Sqrt((m1 + m2)/Power(b,3))*(1 - (m1*m2)/(b*(m1 + m2)))* (1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b)*yp* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 +
       (m1*(-(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2))
       - Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),
       2)*Power(zp,2) + 2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3))))/Power(b,5) + (Power(m1,2)*(m1 +
       m2)*Power(yp,2)*Power(zp,2)* ((m1*(-(Power(1 + (m1*(1 +
       m1/(2.*(m1 + m2))))/b,2)*Power(yp,2)) - Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* Power(1 - m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       - Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2),2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)/ Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),2) -
       (4*m1*Sqrt((m1 + m2)/Power(b,3))*Sqrt((m1 + m2)/b)* (1 -
       (m1*m2)/(b*(m1 + m2)))* ((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)))* (-(Power(1 + (m1*(1 + m1/(2.*(m1
       + m2))))/b,2)*Power(yp,2)) + Power(1 + m1/b - (m1*(-((b*m1)/(m1
       + m2)) - xp))/Power(b,2),2)* Power(zp,2) - Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 - m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),2)* Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2), 2)*Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),4))/Power(b,3)
       + ((m1 + m2)*Power(1 - (m1*m2)/(b*(m1 + m2)),2)* (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))*
       Power(1 + m2/ (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 +
       m2))))/b,2)* Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 +
       m2)) - xp))/Power(b,2), 2)*Power(zp,2) + Power((1 +
       m1/b)*(-((b*m1)/(m1 + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 +
       m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))),4)* (1 + (m1*(-(Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)* Power(yp,2)) - Power(1 + m1/b
       - (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2), 2)*Power(zp,2) +
       2*Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2))* ((-2*Power(m2,2))/ (Power(1 +
       (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*Power(yp,2) + Power(1 + m1/b -
       (m1*(-((b*m1)/(m1 + m2)) - xp))/Power(b,2),2)* Power(zp,2) +
       Power((1 + m1/b)*(-((b*m1)/(m1 + m2)) - xp) -
       (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) + Power(yp,2) -
       Power(zp,2)))/(2.*Power(b,2)),2)) + Power(1 + m2/
       (2.*Sqrt(Power(1 + (m1*(1 + m1/(2.*(m1 + m2))))/b,2)*
       Power(yp,2) + Power(1 + m1/b - (m1*(-((b*m1)/(m1 + m2)) -
       xp))/Power(b,2),2)* Power(zp,2) + Power((1 + m1/b)*(-((b*m1)/(m1
       + m2)) - xp) - (m1*(Power(-((b*m1)/(m1 + m2)) - xp,2) +
       Power(yp,2) - Power(zp,2)))/(2.*Power(b,2)),2))),
       4)))/Power(b,3)))/Power(b,3)))/Power(b,7) ;

*results=g00;
*(results+1)=g0x;
*(results+2)=g0y;
*(results+3)=g0z;
*(results+4)=gxx;
*(results+5)=gxy;
*(results+6)=gxz;
*(results+7)=gyy;
*(results+8)=gyz;
*(results+9)=gzz; 


}

void Alvireg3(double *x, double *y, double *z, double *zm1, double *zm2, double *zr12, double *results) {

double g00,g0x,g0y,g0z,gxx,gxy,gxz,gyy,gyz,gzz;

double m1;
double m2;
double b;
double xp;
double yp;
double zp;

m1 = *zm1;
m2 = *zm2;
b=*zr12;
xp = *x;
yp = *y;
zp = *z;

g00 = 
       -1 + (2*m1)/Sqrt(Power(-b/2. + xp,2) + Power(yp,2) +
       Power(zp,2)) + (2*m2)/Sqrt(Power(b/2. + xp,2) + Power(yp,2) +
       Power(zp,2)) - (m1*m2*Power(yp,2)*(m2/ Power(Power(-b/2. + xp,2)
       + Power(yp,2) + Power(zp,2),1.5) + m1/Power(Power(b/2. + xp,2) +
       Power(yp,2) + Power(zp,2),1.5)))/ (b*(m1 + m2)) - (7*m1*m2*(1/
       Sqrt(Power(-b/2. + xp,2) + Power(yp,2) + Power(zp,2)) -
       1/Sqrt(Power(b/2. + xp,2) + Power(yp,2) + Power(zp,2))))/b -
       (2*m1*m2*(1/Sqrt(Power(-b/2. + xp,2) + Power(yp,2) +
       Power(zp,2)) + 1/Sqrt(Power(b/2. + xp,2) + Power(yp,2) +
       Power(zp,2))))/b + (3*m1*m2*(m2/Sqrt(Power(-b/2. + xp,2) +
       Power(yp,2) + Power(zp,2)) + m1/Sqrt(Power(b/2. + xp,2) +
       Power(yp,2) + Power(zp,2))))/ (b*(m1 + m2)) - 2*Power(m1/
       Sqrt(Power(-b/2. + xp,2) + Power(yp,2) + Power(zp,2)) +
       m2/Sqrt(Power(b/2. + xp,2) + Power(yp,2) + Power(zp,2)),2) +
       ((m1 + m2)*(Power(xp,2) + Power(yp,2))* (1 +
       (2*m1)/Sqrt(Power(-b/2. + xp,2) + Power(yp,2) + Power(zp,2)) +
       (2*m2)/Sqrt(Power(b/2. + xp,2) + Power(yp,2) + Power(zp,2))))/
       Power(b,3) ;
g0x = 
       Sqrt((m1 + m2)/Power(b,3))*yp*(1 + (2*m1)/Sqrt(Power(-b/2. +
       xp,2) + Power(yp,2) + Power(zp,2)) + (2*m2)/Sqrt(Power(b/2. +
       xp,2) + Power(yp,2) + Power(zp,2))) ;
g0y = 
       (-4*m1*m2*Sqrt((m1 + m2)/b)*(1/ Sqrt(Power(-b/2. + xp,2) +
       Power(yp,2) + Power(zp,2)) - 1/Sqrt(Power(b/2. + xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2) + Sqrt((m1 +
       m2)/Power(b,3))*xp* (1 + (2*m1)/Sqrt(Power(-b/2. + xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*m2)/Sqrt(Power(b/2. + xp,2) +
       Power(yp,2) + Power(zp,2))) ;
g0z = 0;
gxx = 
       1 + (2*m1)/Sqrt(Power(-b/2. + xp,2) + Power(yp,2) + Power(zp,2))
       + (2*m2)/Sqrt(Power(b/2. + xp,2) + Power(yp,2) + Power(zp,2)) ;
gxy = 0;
gxz = 0;
gyy = 
       1 + (2*m1)/Sqrt(Power(-b/2. + xp,2) + Power(yp,2) + Power(zp,2))
       + (2*m2)/Sqrt(Power(b/2. + xp,2) + Power(yp,2) + Power(zp,2)) ;
gyz = 0;
gzz = 
       1 + (2*m1)/Sqrt(Power(-b/2. + xp,2) + Power(yp,2) + Power(zp,2))
       + (2*m2)/Sqrt(Power(b/2. + xp,2) + Power(yp,2) + Power(zp,2)) ;

*results=g00;
*(results+1)=g0x;
*(results+2)=g0y;
*(results+3)=g0z;
*(results+4)=gxx;
*(results+5)=gxy;
*(results+6)=gxz;
*(results+7)=gyy;
*(results+8)=gyz;
*(results+9)=gzz; 

}

void Alvireg4(double *x, double *y, double *z, double *zm1, double *zm2, double *zr12, double *results) {

double g00,g0x,g0y,g0z,gxx,gxy,gxz,gyy,gyz,gzz;

double m1;
double m2;
double b;
double xp;
double yp;
double zp;

m1 = *zm1;
m2 = *zm2;
b=*zr12;
xp = *x;
yp = *y;
zp = *z;

g00 = 
       -1 - (2*Power(m1 + m2,2))/(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (2*(m1 + m2)*(1 - (m1*m2)/(2.*b*(m1 + m2))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) - (12*m1*(m1 -
       m2)*m2*Power((m1 + m2)/b,2.5)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ (Power(b,2)*Power(m1 +
       m2,2)*(Power(xp,2) + Power(yp,2) + Power(zp,2))) + (4*m1*m2*(m1
       + m2)*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2)
       + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)* (1 - ((m1 - m2)*((Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/ (Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2)))/
       (Power(b,4)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) +
       (2*m1*m2*(m1 + m2)*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)* (-2/b + ((m1 - m2)*((Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/b +
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       (Power(b,3)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) +
       (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (2*(m1 + m2)*(Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) - Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/b + (Power(b,2)*(-Power(xp,2) -
       Power(yp,2) - Power(zp,2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ ((m1 + m2)*Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),1.5)) + (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 +
       m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)*(yp* cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2)) + (2*(m1 + m2)*(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* ((-4*m1*m2*(Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*(Power(xp,2) + Power(yp,2) + Power(zp,2))) + (2*m1*(m1 -
       m2)*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,2)*(-Power(xp,2) - Power(yp,2) -
       Power(zp,2) + 3*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))))/ (Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       ((m1 + m2)*(2*Power(yp* cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) - 3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/b) )/(Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2), 1.5))))/Power(b,2) -
       (4*m1*m2*Sqrt((m1 + m2)/b)*(yp* cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* ((2*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/b - ((m1 - m2)*((4*Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))*
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/b +
       (-Power(xp,2) - Power(yp,2) - Power(zp,2) +
       3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       (b*(Power(xp,2) + Power(yp,2) + Power(zp,2))) + ((m1 +
       m2)*(Power(xp,2) + Power(yp,2))* (1 + Power(m1 +
       m2,2)/(Power(xp,2) + Power(yp,2) + Power(zp,2)) + (2*(m1 +
       m2)*(1 - (m1*m2)/(2.*b*(m1 + m2))))/ Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*m1*(m1 - m2)*m2*(Sqrt((m1 +
       m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ (b*(m1 +
       m2)*(Power(xp,2) + Power(yp,2) + Power(zp,2))) +
       (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (2*(m1 + m2)*(Power(yp* cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) - Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/b + (Power(b,2)*(-Power(xp,2) -
       Power(yp,2) - Power(zp,2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ ((m1 + m2)*Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),1.5)) + (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 +
       m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2)))) /Power(b,3) ;
g0x =
       (4*m1*m2*Power((m1 + m2)/b,1.5)* sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))* (xp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (1 - ((m1 -
       m2)*((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/ (Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2)))/ (b*(m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) -
       (2*m1*m2*Sqrt((m1 + m2)/b)*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))* (yp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-2/b + ((m1 -
       m2)*((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/b + (xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       ((m1 + m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       (b*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + Sqrt((m1 +
       m2)/b)*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))* ((-4*m1*m2*(Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*(Power(xp,2) + Power(yp,2) + Power(zp,2))) + (2*m1*(m1 -
       m2)*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,2)*(-Power(xp,2) - Power(yp,2) -
       Power(zp,2) + 3*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))))/ (Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       ((m1 + m2)*(2*Power(yp* cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) - 3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/b))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),1.5))) + (2*m1*m2*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))*
       ((2*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/b - ((m1 -
       m2)*((4*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/b + (-Power(xp,2) - Power(yp,2) - Power(zp,2) +
       3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       (Power(xp,2) + Power(yp,2) + Power(zp,2)) - Sqrt((m1 +
       m2)/Power(b,3))*yp* (1 + Power(m1 + m2,2)/(Power(xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*(m1 + m2)*(1 - (m1*m2)/(2.*b*(m1
       + m2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (2*m1*(m1 - m2)*m2*(Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))) + (b*(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ (b*(m1 + m2)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))*
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) + (2*(m1 +
       m2)*(Power(yp* cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2)
       + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) - Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/b + (Power(b,2)*(-Power(xp,2) - Power(yp,2) -
       Power(zp,2) + 3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       ((m1 + m2)*Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5)) +
       (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2))) + (6*m1*(m1 -
       m2)*m2*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))*
       (xp*cos(2*Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(2*Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ (Power(b,3)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) ;
g0y =
       (4*m1*m2*Power((m1 + m2)/b,1.5)* cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))* (xp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (1 - ((m1 -
       m2)*((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/ (Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2)))/ (b*(m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) +
       (2*m1*m2*Sqrt((m1 + m2)/b)*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))* (yp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-2/b + ((m1 -
       m2)*((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/b + (xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       ((m1 + m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       (b*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + Sqrt((m1 +
       m2)/b)*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))* ((-4*m1*m2*(Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*(Power(xp,2) + Power(yp,2) + Power(zp,2))) + (2*m1*(m1 -
       m2)*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,2)*(-Power(xp,2) - Power(yp,2) -
       Power(zp,2) + 3*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))))/ (Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       ((m1 + m2)*(2*Power(yp* cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) - 3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/b))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),1.5))) + (2*m1*m2*sin(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))*
       ((2*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/b - ((m1 -
       m2)*((4*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/b + (-Power(xp,2) - Power(yp,2) - Power(zp,2) +
       3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       (Power(xp,2) + Power(yp,2) + Power(zp,2)) + Sqrt((m1 +
       m2)/Power(b,3))*xp* (1 + Power(m1 + m2,2)/(Power(xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*(m1 + m2)*(1 - (m1*m2)/(2.*b*(m1
       + m2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (2*m1*(m1 - m2)*m2*(Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))) + (b*(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ (b*(m1 + m2)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))*
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) + (2*(m1 +
       m2)*(Power(yp* cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2)
       + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) - Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/b + (Power(b,2)*(-Power(xp,2) - Power(yp,2) -
       Power(zp,2) + 3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       ((m1 + m2)*Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5)) +
       (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2))) - (6*m1*(m1 -
       m2)*m2*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))*
       (yp*cos(2*Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(2*Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ (Power(b,3)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) ;
g0z = 0;
gxx = 
       1 + (Power(m1 + m2,2)*Power(xp,2))/ Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2) + Power(m1 + m2,2)/(Power(xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*(m1 + m2)*(1 - (m1*m2)/(2.*b*(m1
       + m2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (2*m1*(m1 - m2)*m2*(Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))) + (b*(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ (b*(m1 + m2)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + (4*m1*m2*Power(sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)* (1 - ((m1 -
       m2)*((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/ (Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2)))/ (b*Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + (2*m1*m2*Power(cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)* (-2/b + ((m1 - m2)*((Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/b +
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (2*(m1 + m2)*(Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) - Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/b + (Power(b,2)*(-Power(xp,2) -
       Power(yp,2) - Power(zp,2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ ((m1 + m2)*Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),1.5)) + (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 +
       m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)*(yp* cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2)) + (6*m1*(m1 - m2)*m2*Power((m1 +
       m2)/b,1.5)* sin(2*Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/ (Power(m1 + m2,2)*(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) ;
gxy =
       (Power(m1 + m2,2)*xp*yp)/Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),2) + (6*m1*(m1 - m2)*m2*Power((m1 + m2)/b,1.5)*
       cos(2*Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)))* (xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ (Power(m1 + m2,2)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + (((4*m1*m2*(1 - ((m1 - m2)*((Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/ (Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2)))/ (b*Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - (2*b*m1*m2*(-2/b + ((m1 - m2)*
       ((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/b + (xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       ((m1 + m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))* sin(2*Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/2. ;
gxz =
       (Power(m1 + m2,2)*xp*zp)/Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),2) ;
gyy =
       1 + (Power(m1 + m2,2)*Power(yp,2))/ Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2) + Power(m1 + m2,2)/(Power(xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*(m1 + m2)*(1 - (m1*m2)/(2.*b*(m1
       + m2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (2*m1*(m1 - m2)*m2*(Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))) + (b*(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ (b*(m1 + m2)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + (4*m1*m2*Power(cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)* (1 - ((m1 -
       m2)*((Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (b*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/ (Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/(m1 + m2)))/ (b*Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + (2*m1*m2*Power(sin(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2)* (-2/b + ((m1 - m2)*((Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/b +
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (2*(m1 + m2)*(Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) - Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/b + (Power(b,2)*(-Power(xp,2) -
       Power(yp,2) - Power(zp,2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ ((m1 + m2)*Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),1.5)) + (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 +
       m2)/b)* (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)*(yp* cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2)) - (6*m1*(m1 - m2)*m2*Power((m1 +
       m2)/b,1.5)* sin(2*Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))/ (Power(m1 + m2,2)*(Power(xp,2) +
       Power(yp,2) + Power(zp,2))); 
gyz =
       (Power(m1 + m2,2)*yp*zp)/Power(Power(xp,2) + Power(yp,2) +
       Power(zp,2),2); 
gzz =  
       1 + (Power(m1 + m2,2)*Power(zp,2))/ Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2) + Power(m1 + m2,2)/(Power(xp,2) +
       Power(yp,2) + Power(zp,2)) + (2*(m1 + m2)*(1 - (m1*m2)/(2.*b*(m1
       + m2))))/ Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) +
       (2*m1*(m1 - m2)*m2*(Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))) + (b*(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2)))))/ Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))/ (b*(m1 + m2)*(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + (m1*m2*((6*b*Sqrt((m1 + m2)/b)* (yp*cos(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))*
       (xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))))/
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)) + (2*(m1 +
       m2)*(Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2) -
       Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/b +
       (Power(b,2)*(-Power(xp,2) - Power(yp,2) - Power(zp,2) +
       3*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       (Power(xp,2) + Power(yp,2) + Power(zp,2))))/ ((m1 +
       m2)*Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5)) +
       (m1*(m1 - m2)*m2*((3*Power(b,2)*Sqrt((m1 + m2)/b)*
       (yp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (Power(xp,2) +
       Power(yp,2) + Power(zp,2) - 5*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ (Power(xp,2) + Power(yp,2) +
       Power(zp,2)) + (Power(b,3)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))))* (3*(Power(xp,2) + Power(yp,2) + Power(zp,2)) -
       5*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))),2)))/
       Power(Power(xp,2) + Power(yp,2) + Power(zp,2),1.5) + (2*(m1 +
       m2)*(xp*cos(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))))* (-Power(xp,2) -
       Power(yp,2) - Power(zp,2) - 6*Power(yp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2) + 3*Power(xp*cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       - yp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))),2)))/ Sqrt(Power(xp,2) + Power(yp,2)
       + Power(zp,2)) + Power((m1 + m2)/b,1.5)*(yp* cos(Sqrt((m1 +
       m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2)))
       + xp*sin(Sqrt((m1 + m2)/Power(b,3))* Sqrt(Power(xp,2) +
       Power(yp,2) + Power(zp,2))))* (-Power(xp,2) - Power(yp,2) -
       Power(zp,2) - 2*Power(yp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) + xp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2) + 7*Power(xp*cos(Sqrt((m1 + m2)/Power(b,3))*
       Sqrt(Power(xp,2) + Power(yp,2) + Power(zp,2))) - yp*sin(Sqrt((m1
       + m2)/Power(b,3))* Sqrt(Power(xp,2) + Power(yp,2) +
       Power(zp,2))),2))))/ (Power(m1 + m2,2)*Power(Power(xp,2) +
       Power(yp,2) + Power(zp,2),2)) ;

*results=g00;
*(results+1)=g0x;
*(results+2)=g0y;
*(results+3)=g0z;
*(results+4)=gxx;
*(results+5)=gxy;
*(results+6)=gxz;
*(results+7)=gyy;
*(results+8)=gyz;
*(results+9)=gzz; 


}
