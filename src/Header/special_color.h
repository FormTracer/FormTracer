* protected variables: T,f,I2R,NR,NA(,cR,cA)
* for xxxfund
Tensor T, f(antisymmetric);
Tensor scOlTp;
CFunction scOlTr(cyclic);
Symbols I2R,NR,NA;
Dimension NA;
AutoDeclare Index scOlj;
Dimension NR;
AutoDeclare Index scOli;

* for checkxxx
CFunction scOlcolor;
Symbols  cR,cA;
Symbols scOlx,scOlx1,scOlx2; 

* for SU3fundexplicit
CFunction scOld(symmetric);
 

#procedure SUNfund
*	USES COLOR.H CONVENTIONS FOR GENERATORS
*	Procedure to compute color traces for the SU(NR) groups
*	We follow the article by Cvitanovic (Phys.Rev.D14(1976)1536)
*
*	We use [T(i),T(j)] = i_*f(i,j,k)*T(k)    (f is the C of Cvitanovic)
*
*	We use the indices i in the space of the fundamental representation
*	The indices j are in the space of the adjoint.
*	The dimension should be the dimension of the fundamental representation.
*
*	We need: (C)Function Tr(cyclic). Indicates the traces.
*	CFunction T, scOlTp, f(antisymmetric);
*	Symbols I2R,NR,NA;
*	Indices scOli1=NR,scOli2=NR,scOli3=NR,scOli4=NR;
*	Indices scOlj1=NA,scOlj2=NA,scOlj3=NA;
*	Dimension NR;
*
*	Usually the value of I2R is taken to be 1/2;
*	NR is the dimension of the fundamental representation.
*	NA is the dimension of the adjoint representation.
*
*	based on routine by J.Vermaseren
*

* step 1:express f via generators
* step 2: identify contracted generators and apply completeness relation
* these are essentially two while loops
#do scOlpri = 1,1
if ( count(f,1) || match(T(scOli1?,scOli2?,scOlj1?!fixed_)*T(scOli3?,scOli4?,scOlj1?!fixed_)) )
		redefine scOlpri "0";
		
id,once,f(scOlj1?,scOlj2?,scOlj3?) = 1/I2R/i_*T(scOli1,scOli2,scOlj1)*T(scOli2,scOli3,scOlj2)*T(scOli3,scOli1,scOlj3)
			-1/I2R/i_*T(scOli1,scOli2,scOlj3)*T(scOli2,scOli3,scOlj2)*T(scOli3,scOli1,scOlj1);
sum scOli1,scOli2,scOli3;
id	T(scOli1?,scOli2?,scOlj1?!fixed_)*T(scOli3?,scOli4?,scOlj1?!fixed_) = scOlTp(scOli1,scOli2,scOli3,scOli4);

#do scOlprj = 1,1
if ( count(scOlTp,1) ) redefine scOlprj "0";
.sort
id,once,scOlTp(scOli1?,scOli2?,scOli3?,scOli4?) =
			I2R*(d_(scOli1,scOli4)*d_(scOli2,scOli3)-d_(scOli1,scOli2)*d_(scOli3,scOli4)/NR);
#enddo
#enddo

*transform everything into traces (consistent with other routines below)
repeat;
	id	T(scOli1?,scOli2?!fixed_,scOlj1?,?a)*T(scOli2?!fixed_,scOli3?,scOlj2?,?b) = T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
endrepeat;
id	T(scOli1?!fixed_,scOli1?!fixed_,?a) = scOlTr(?a);
id	scOlTr(scOlj1?) = 0;
id	scOlTr(scOlj1?,scOlj2?) = I2R*d_(scOlj1,scOlj2);
.sort

#endprocedure

#procedure SONfund
*	USES COLOR.H CONVENTIONS FOR GENERATORS
*	Procedure to compute color traces for the SO(NF) groups
*	We follow the article by Cvitanovic (Phys.Rev.D14(1976)1536
*
*	We use [T(i),T(j)] = i_*f(i,j,k)*T(k)    (f is the C of Cvitanovic)
*
*	We use the indices i in the space of the fundamental representation
*	The indices j are in the space of the adjoint.
*	The dimension should be the dimension of the fundamental representation.
*
*	We need: (C)Function Tr(cyclic). Indicates the traces.
*	CFunction T, scOlTp, f(antisymmetric);
*	Symbols I2R,NR,NA;
*	Indices scOli1=NR,scOli2=NR,scOli3=NR,scOli4=NR;
*	Indices scOlj1=NA,scOlj2=NA,scOlj3=NA;
*	Dimension NR;
*
*	Usually the value of I2R is taken to be 1/2;
*	NR is the dimension of the fundamental representation.
*	NA is the dimension of the adjoint representation.
*
*	based on routine by J.Vermaseren
*

#do scOlpri = 1,1
if ( count(f,1) || match(T(scOli1?,scOli2?,scOlj1?!fixed_)*T(scOli3?,scOli4?,scOlj1?!fixed_)) )
		redefine scOlpri "0";
		
id,once,f(scOlj1?,scOlj2?,scOlj3?) = 2/I2R/i_*T(scOli1,scOli2,scOlj1)*T(scOli2,scOli3,scOlj2)*T(scOli3,scOli1,scOlj3);
sum scOli1,scOli2,scOli3;

id	T(scOli1?,scOli2?,scOlj1?!fixed_)*T(scOli3?,scOli4?,scOlj1?!fixed_) = scOlTp(scOli1,scOli2,scOli3,scOli4);
#do scOlpri = 1,1
if ( count(scOlTp,1) ) redefine scOlpri "0";
.sort
id,once,scOlTp(scOli1?,scOli2?,scOli3?,scOli4?) =
			I2R/2*(d_(scOli1,scOli4)*d_(scOli2,scOli3)-d_(scOli1,scOli3)*d_(scOli2,scOli4));
#enddo
#enddo

repeat;
	id	T(scOli1?,scOli2?!fixed_,scOlj1?,?a)*T(scOli2?!fixed_,scOli3?,scOlj2?,?b) =  T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
	id	T(scOli1?,scOli2?!fixed_,scOlj1?,?a)*T(scOli3?,scOli2?!fixed_,scOlj2?,?b) = -T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
	id	T(scOli2?!fixed_,scOli1?,scOlj1?,?a)*T(scOli2?!fixed_,scOli3?,scOlj2?,?b) = -T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
	id	T(scOli2?!fixed_,scOli1?,scOlj1?,?a)*T(scOli3?,scOli2?!fixed_,scOlj2?,?b) =  T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
endrepeat;
id	T(scOli1?!fixed_,scOli1?!fixed_,?a) = scOlTr(?a);
id	scOlTr(scOlj1?,scOlj2?) = I2R*d_(scOlj1,scOlj2);
id	scOlTr(scOlj1?,scOlj2?,scOlj3?) = i_*I2R/2*f(scOlj1,scOlj2,scOlj3);
.sort
#endprocedure

#procedure SPNfund
*	USES COLOR.H CONVENTIONS FOR GENERATORS
*	Procedure to compute color traces for the SP(NF) groups
*	We follow the article by Cvitanovic (Phys.Rev.D14(1976)1536
*
*	We use [T(i),T(j)] = i_*f(i,j,k)*T(k)    (f is the C of Cvitanovic)
*
*	We use the indices i in the space of the fundamental representation
*	The indices j are in the space of the adjoint.
*	The dimension should be the dimension of the fundamental representation.
*
*	We need: (C)Function Tr(cyclic). Indicates the traces.
*	CFunction T, scOlTp, f(antisymmetric);
*	Symbols I2R,NR,NA;
*	Indices scOli1=NR,scOli2=NR,scOli3=NR,scOli4=NR;
*	Indices scOlj1=NA,scOlj2=NA,scOlj3=NA;
*	Dimension NR;
*
*	Usually the value of I2R is taken to be 1/2;
*	NR is the dimension of the fundamental representation.
*	NA is the dimension of the adjoint representation.
*
*	based on routine by J.Vermaseren, 7-jan-1997
*

#do scOlpri = 1,1
if ( count(f,1) || match(T(scOli1?,scOli2?,scOlj1?!fixed_)*T(scOli3?,scOli4?,scOlj1?!fixed_)) )
		redefine scOlpri "0";
		
id,once,f(scOlj1?,scOlj2?,scOlj3?) = 2/I2R/i_*T(scOli1,scOli2,scOlj1)*T(scOli2,scOli3,scOlj2)*T(scOli3,scOli1,scOlj3);
sum scOli1,scOli2,scOli3;

id	T(scOli1?,scOli2?,scOlj1?!fixed_)*T(scOli3?,scOli4?,scOlj1?!fixed_) = scOlTp(scOli1,scOli2,scOli3,scOli4);
#do scOlpri = 1,1
if ( count(scOlTp,1) ) redefine scOlpri "0";
.sort
id,once,scOlTp(scOli1?,scOli2?,scOli3?,scOli4?) =
			I2R/2*(d_(scOli1,scOli4)*d_(scOli2,scOli3)+f(scOli1,scOli3)*f(scOli2,scOli4));
repeat id	f(scOli1?,scOli2?)*f(scOli1?,scOli3?) = -d_(scOli2,scOli3);
#enddo
#enddo

repeat;
	id	T(scOli1?,scOli2?!fixed_,scOlj1?,?a)*T(scOli2?!fixed_,scOli3?,scOlj2?,?b) =  T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
	id	T(scOli1?,scOli2?!fixed_,scOlj1?,?a)*T(scOli3?,scOli2?!fixed_,scOlj2?,?b) = -T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
	id	T(scOli2?!fixed_,scOli1?,scOlj1?,?a)*T(scOli2?!fixed_,scOli3?,scOlj2?,?b) = -T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
	id	T(scOli2?!fixed_,scOli1?,scOlj1?,?a)*T(scOli3?,scOli2?!fixed_,scOlj2?,?b) =  T(scOli1,scOli3,scOlj1,?a,scOlj2,?b);
endrepeat;
id	T(scOli1?!fixed_,scOli1?!fixed_,?a) = scOlTr(?a);
id	scOlTr(scOlj1?) = 0;
id	scOlTr(scOlj1?,scOlj2?) = I2R*d_(scOlj1,scOlj2);
id	scOlTr(scOlj1?,scOlj2?,scOlj3?) = i_*I2R/2*f(scOlj1,scOlj2,scOlj3);
.sort
#endprocedure

#procedure checkSUNfund
*performs an SUn color trace (fundamental rep) and tries to express some constants via Casimirs
*	cA is the quadratic Casimir of the adjoint rep (2*I2R*NR)
*	cR is the quadratic Casimir of the fundamental rep I2R*(NR^2-1)/NR
* 	cR*NR = I2R*NA e.g. from hep-ph/9701390 (valid for all groups)

#call SUNfund
B   I2R;
.sort
Collect scOlcolor;
FactArg scOlcolor;
ChainOut,scOlcolor;
id  scOlcolor(scOlx?number_) = scOlx;
id  scOlcolor(scOlx?symbol_) = scOlx;
id  scOlcolor(NR-1)*scOlcolor(NR+1) = cR*NR/I2R;
id  scOlcolor(scOlx?) = scOlx;
id  cR*NR = I2R*NA;
id  I2R*NR = cA/2;
#endprocedure

#procedure checkSONfund
*performs an SOn color trace (fundamental rep) and tries to express some constants via Casimirs
*	cA is the quadratic Casimir of the adjoint rep I2R*(NR-2)
*	cR is the quadratic Casimir of the fundamental rep I2R/2*(NR-1)

#call SONfund
B   I2R;
.sort
Collect scOlcolor;
FactArg scOlcolor;
ChainOut,scOlcolor;
id  scOlcolor(scOlx?number_) = scOlx;
id  scOlcolor(scOlx?symbol_) = scOlx;
id  scOlcolor(NR-2) = cA/I2R;
id  scOlcolor(NR-1) = 2*cR/I2R;
id  scOlcolor(scOlx?) = scOlx;
id  cR*NR = I2R*NA;
#endprocedure

#procedure checkSPNfund
* performs an Spn color trace (fundamental rep) tries to express some constants via Casimirs
*	cA is the quadratic Casimir of the adjoint rep I2R*(NR+2)
*	cR is the quadratic Casimir of the fundamental rep I2R/2*(NR+1)
#call SPNfund
B   I2R;
.sort
Collect scOlcolor;
FactArg scOlcolor;
ChainOut,scOlcolor;
id  scOlcolor(scOlx?number_) = scOlx;
id  scOlcolor(scOlx?symbol_) = scOlx;
id  scOlcolor(NR+2) = cA/I2R;
id  scOlcolor(NR+1) = 2*cR/I2R;
id  scOlcolor(scOlx?) = scOlx;
id  cR*NR = I2R*NA;
#endprocedure

* carries out a fundamental SU(3) trace (allowing explicit indices, partial traces, transposed generators)
* dimension has to be set to NR (for SUNfund)
#procedure SU3fundexplicit

* handle transposed generators
Repeat;
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,1)=T(scOli1,scOli3,scOlj1,?scOlj2,1);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,2)=-T(scOli1,scOli3,scOlj1,?scOlj2,2);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,3)=T(scOli1,scOli3,scOlj1,?scOlj2,3);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,4)=T(scOli1,scOli3,scOlj1,?scOlj2,4);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,5)=-T(scOli1,scOli3,scOlj1,?scOlj2,5);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,6)=T(scOli1,scOli3,scOlj1,?scOlj2,6);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,7)=-T(scOli1,scOli3,scOlj1,?scOlj2,7);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,8)=T(scOli1,scOli3,scOlj1,?scOlj2,8);

id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,1)=T(scOli3,scOli1,1,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,2)=-T(scOli3,scOli1,2,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,3)=T(scOli3,scOli1,3,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,4)=T(scOli3,scOli1,4,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,5)=-T(scOli3,scOli1,5,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,6)=T(scOli3,scOli1,6,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,7)=-T(scOli3,scOli1,7,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,8)=T(scOli3,scOli1,8,scOlj1,?scOlj2);

id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,scOlj3?!fixed_)=-4*T(scOli1,scOli3,scOlj1,?scOlj2,2,scOlj3,2)-4*T(scOli1,scOli3,scOlj1,?scOlj2,5,scOlj3,5)-4*T(scOli1,scOli3,scOlj1,?scOlj2,7,scOlj3,7);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,scOlj3?!fixed_)=-4*T(scOli3,scOli1,2,scOlj3,2,scOlj1,?scOlj2)-4*T(scOli3,scOli1,5,scOlj3,5,scOlj1,?scOlj2)-4*T(scOli3,scOli1,7,scOlj3,7,scOlj1,?scOlj2);
EndRepeat;
.sort

* undo combined T for SUNfund
Repeat;
id,once T(scOli1?,scOli2?,scOlj1?,scOlj2?,?scOlj3)=T(scOli1,scOli3,scOlj1)*T(scOli3,scOli2,scOlj2,?scOlj3);
sum scOli3;
EndRepeat;
.sort

* call standard SUNfund trace
#call SUNfund
.sort
* set dimension to NR=3 for sum below
Dimension 3;
Multiply replace_(NR,3,NA,8);

* undo scOlTr
Repeat;
id,once scOlTr(?scOljx)=T(scOlix,scOlix,?scOljx);
sum scOlix;
Endrepeat;

*explicitly reduce T
Repeat;
id,once T(scOli1?,scOli3?,scOlj1?,scOlj2?,?scOlj3)=1/2/3*d_(scOlj1,scOlj2)*T(scOli1,scOli3,?scOlj3)+1/2*sum_(scOlj4, 1, 8, (I_*f(scOlj1,scOlj2,scOlj4)+scOld(scOlj1,scOlj2,scOlj4))*T(scOli1,scOli3,scOlj4,?scOlj3));
id T(scOli1?,scOli1?,scOlj?)=0;
id T(scOli1?,scOli2?)=d_(scOli1,scOli2);

id f(1,2,3)=1;
id f(1,4,7)=1/2;
id f(2,4,6)=1/2;
id f(2,5,7)=1/2;
id f(3,4,5)=1/2;
id f(1,5,6)=-1/2;
id f(3,6,7)=-1/2;
id f(4,5,8)=sqrt_(3)/2;
id f(6,7,8)=sqrt_(3)/2;
id scOld(1,1,8)=sqrt_(3)/3;
id scOld(2,2,8)=sqrt_(3)/3;
id scOld(3,3,8)=sqrt_(3)/3;
id scOld(8,8,8)=-sqrt_(3)/3;
id scOld(4,4,8)=-sqrt_(3)/6;
id scOld(5,5,8)=-sqrt_(3)/6;
id scOld(6,6,8)=-sqrt_(3)/6;
id scOld(7,7,8)=-sqrt_(3)/6;
id scOld(1,4,6)=1/2;
id scOld(1,5,7)=1/2;
id scOld(2,4,7)=-1/2;
id scOld(2,5,6)=1/2;
id scOld(3,4,4)=1/2;
id scOld(3,5,5)=1/2;
id scOld(3,6,6)=-1/2;
id scOld(3,7,7)=-1/2;

id f(scOlj1?fixed_,scOlj2?fixed_,scOlj3?fixed_)=0;
id scOld(scOlj1?fixed_,scOlj2?fixed_,scOlj3?fixed_)=0;
id sqrt_(3)^2=3;
Endrepeat;
.sort;
#endprocedure

* explicitly carries out a SU(2) trace (allowing explicit indices, partial traces, transposed generators)
* Dimension has to be set to NA=3
#procedure SU2fundexplicit

Multiply replace_(NR,2,NA,3);

id f(scOlj1?,scOlj2?,scOlj3?)=e_(scOlj1,scOlj2,scOlj3);
Contract 0;

* handle transposed generators
Repeat;
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,1)=T(scOli1,scOli3,scOlj1,?scOlj2,1);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,2)=-T(scOli1,scOli3,scOlj1,?scOlj2,2);
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,3)=T(scOli1,scOli3,scOlj1,?scOlj2,3);

id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,1)=T(scOli3,scOli1,1,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,2)=-T(scOli3,scOli1,2,scOlj1,?scOlj2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,3)=-T(scOli3,scOli1,3,scOlj1,?scOlj2);

id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli3?,scOli2?!fixed_,scOlj3?!fixed_)=-4*T(scOli1,scOli3,scOlj1,?scOlj2,2,scOlj3,2);
id T(scOli2?!fixed_,scOli1?,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,scOlj3?!fixed_)=-4*T(scOli3,scOli1,2,scOlj3,2,scOlj1,?scOlj2);
EndRepeat;

Repeat;
id T(scOli1?,scOli2?!fixed_,scOlj1?,?scOlj2)*T(scOli2?!fixed_,scOli3?,scOlj3?,?scOlj4)=T(scOli1,scOli3,scOlj1,?scOlj2,scOlj3,?scOlj4);
EndRepeat;

Repeat;
id T(scOli1?!fixed_,scOli1?!fixed_,scOlj1?)=0;
id,once T(scOli1?,scOli2?,scOlj1?,scOlj2?,?scOlj3)=I_/2*e_(scOlj1,scOlj2,scOlj4)*T(scOli1,scOli2,scOlj4,?scOlj3)+1/4*d_(scOlj1,scOlj2)*T(scOli1,scOli2,?scOlj3);
Sum scOlj4;
id T(scOli1?,scOli2?)=d_(scOli1,scOli2);
Contract;
EndRepeat;
.sort
#endprocedure
