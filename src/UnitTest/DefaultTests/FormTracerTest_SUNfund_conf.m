Hold[If[versionStringToNumber[$FormTracerVersionNumber] > 1.7, 
  {DefineGroupTensors[{SUNfund, {color, Nc}, delta[adj, a, b], SUNF[a, b, c], 
      delta[FundCol, a, b], SUNT[a, i, j]}]; }, 
  If[versionStringToNumber[$FormTracerVersionNumber] >= 1.5, 
   {DefineGroupTensors[{SUN, {color, Ncolor}, delta[adj, a, b], 
       SUNF[a, b, c], delta[FundCol, a, b], SUNT[a, i, j]}]; }, 
   {DefineGroupTensors[{{color, delta[adj, a, b], SUNF[a, b, c], 
        delta[FundCol, a, b], SUNT[a, i, j]}}]; }]]]
