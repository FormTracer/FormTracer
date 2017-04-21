*workaround to replace epsilon symbols with an arbitrary number of arguments
#procedure replaceepsilon
AntiBracket e_;
.sort 
Tensor FTxeps;
Vector epsilonvec;
CFunction epsilonguard;
Symbol epsilontmp;
Collect epsilonguard;
FactArg epsilonguard;
Chainout epsilonguard;
Argument;
ToVector e_,epsilonvec;
ToTensor epsilonvec,FTxeps;
EndArgument;
id epsilonguard(FTxeps(?x))=FTxeps(?x); 
id epsilonguard(epsilontmp?number_)=epsilontmp;
#endprocedure
