********************************************************************************
* gamma5 handling before Dirac trace in d dimensions using Larin scheme (and (10) in 1506.04517)
* input: single gammastring as gamma(...)
* output: expression to be contracted and then actually traced
* #define FastGamma5Tracedneq4 "1" for faster implementation for single gamma5
* required definitions:
* CFunction FTxTr, FTxTr2, FTxTr3, FTxeInt
********************************************************************************
#procedure PreDiracTrace5

id once gamma (?mom1, spi1?, spi1?) = FTxTr(?mom1);


* replace squares of gamma5s
Repeat;
id FTxTr(?lorx1,g5,g5,?lorx2)=FTxTr(?lorx1,?lorx2);
EndRepeat;

*remove gamma4 from first position equivalent to reading the trace from a certain position
id once FTxTr(g5,?lorx1)=FTxTr(?lorx1,g5); 

*remove another set of gamma^5s
id FTxTr(?lorx1,g5,g5,?lorx2)=FTxTr(?lorx1,?lorx2);

* replace all but one gamma5 using Larin scheme (removed I for Eucl conventions)
Repeat;
if(match(FTxTr(?lorx1,g5,?lorx2,g5,?lorx3)));
id once FTxTr(g5,?lorx1)=FTxTr(?lorx1,g5);
id once FTxTr(?lorx1,lorx3?!setg5,g5,?lorx2)=FTxeInt(lorx3,lorx4,lorx5,lorx6)*FTxTr(?lorx1,lorx4,lorx5,lorx6,?lorx2)/6;
Sum lorx4,lorx5,lorx6;
Sum lorx4,lorx5,lorx6;
endif;
EndRepeat;

* eliminating products of an odd number of gammas
id once FTxTr(?lorx1,g5,?lorx2)=FTxTr2(?lorx2,?lorx1,g5)*FTxTr(?lorx2,?lorx1,g5);
#define MAXGAM "10"
#do i = `MAXGAM'-1,3,-2
id FTxTr2(lor1?,...,lor{`i'}?,g5) = 0;
#enddo
id FTxTr2(lor1?,lor2?,g5) = 0;
id FTxTr2(lor1?,g5) = 0;
id FTxTr2(g5) = 0;
id FTxTr2(?lorx1)=1;

#ifndef `FastGamma5Tracedneq4'
* default: apply Larin again
id once FTxTr(?lorx1,lorx3?,g5,?lorx2)=FTxeInt(lorx3,lorx4,lorx5,lorx6)*FTxTr(?lorx1,lorx4,lorx5,lorx6,?lorx2)/6;
Sum lorx4,lorx5,lorx6;
#else
* handle the remaining gamma5 (implementing eq 10)
repeat;
id once FTxTr(?a,lorx2?,g5) = distrib_(-2,3,FTxTr,FTxTr2,?a)*FTxTr3(lorx2,g5);
id FTxTr2(lorx1?,lorx2?,lorx3?)*FTxTr3(lorx4?,g5) = FTxeInt(lorx1,...,lorx4);
endrepeat;
#endif

#endprocedure


