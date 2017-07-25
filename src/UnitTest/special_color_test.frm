#include ../color.h
#include ../special_color.h



Off statistics;

Indices <i1=NR>,...,<i8=NR>;
Indices <j1=NA>,...,<j8=NA>;
S   cR,cA; 
Dimension 3;

Local colErr=f(j1,j4,j5)*f(j2,j5,j7)*f(j3,j7,j4)*T(i3,i4,j1)*T(i4,i5,j2)*T(i5,i3,j3);
Local colErr2=T(i1,i2,j1)*T(i2,i3,j1)*T(i3,i4,j2)*T(i4,i5,j2)*T(i5,i6,j3)*T(i6,i1,j3);

* #call SU2fundexplicit
#call SU3fundexplicit

#call checkSUNfund
* #call checkSONfund
* #call checkSPNfund
* #call color

id I2R=1/2;

* id NR=2;
* id NA=3;
* id cR=3/4;
* id cA=2;

id NR=3;
id NA=8;
id cR=4/3;
id cA=3;
.sort


print;
.end



