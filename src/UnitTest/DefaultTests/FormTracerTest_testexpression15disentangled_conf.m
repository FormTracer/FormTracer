{Hold[Join[{If[versionStringToNumber[$FormTracerVersionNumber] > 1.7, 
       DefineLorentzTensors[delta[lor, a, b], {vec[p, mu], p, mu}, 
         {sp[p, q], p, q}, eps[], {delta[Spinor, i, j], i, j}, 
         gamma[mu, i, j], gamma5[i, j]]; , 
       DefineLorentzTensors[delta[lor, a, b], {vec[p, mu], p, mu}, 
         {sp[p, q], p, q}, {delta[Spinor, i, j], i, j}, gamma[mu, i, j], 
         gamma5[i, j]]; ]*DisentangleLorentzStructures[False]; 
     DefineLorentzTensorIdentities[{{PiT[p, mu, rho]*PiT[p, rho, nu], 
        PiT[p, mu, nu]}, {PiT[p, mu, rho]*PiT[-p, rho, nu], 
        PiT[p, mu, nu]}}]; DefineCombinedLorentzTensors[
      {{PiT[p, mu, nu], delta[lor, mu, nu] - vec[p, mu]*
          (vec[p, nu]/sp[p, p])}, {PiL[p, mu, nu], 
        vec[p, mu]*(vec[p, nu]/sp[p, p])}}]; DefineFormAutoDeclareFunctions[
      {dR, R, Z, dZ, eta, h, sqrM, lambda, M, alphaS, V}]; }, 
   If[versionStringToNumber[$FormTracerVersionNumber] > 1.7, 
    DefineGroupTensors[{{SUNfund, {color, Nc}, delta[adj, a, b], 
        SUNF[a, b, c], delta[FundCol, a, b], SUNT[a, i, j]}, 
       {SUNfund, {flavor, Nf}, delta[MesonFlav, a, b], SUNFF[a, b, c], 
        delta[FundFlav, a, b], SUNTF[a, i, j]}}]; , 
    If[versionStringToNumber[$FormTracerVersionNumber] > 1.5, 
     {DefineGroupTensors[{{SUN, {color, Nc}, delta[adj, a, b], SUNF[a, b, c], 
          delta[FundCol, a, b], SUNT[a, i, j]}, {SUN, {flavor, Nf}, 
          delta[MesonFlav, a, b], SUNFF[a, b, c], delta[FundFlav, a, b], 
          SUNTF[a, i, j]}}]; }, 
     {DefineGroupTensors[{{color, delta[adj, a, b], SUNF[a, b, c], 
          delta[FundCol, a, b], SUNT[a, i, j]}, {flavor, delta[MesonFlav, a, 
           b], SUNFF[a, b, c], delta[FundFlav, a, b], SUNTF[a, i, 
           j]}}]; }]]]], Hold[DisentangleLorentzStructures[True]]}
