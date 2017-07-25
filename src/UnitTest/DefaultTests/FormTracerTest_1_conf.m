Hold[{DefineLorentzTensors[deltaLorentz[a, b], vec[p, mu], sp[p, q]], 
  DefineCombinedLorentzTensors[
   {{PiT[p, mu, nu], deltaLorentz[mu, nu] - (vec[p, mu]*vec[p, nu])/
       sp[p, p]}, {PiL[p, mu, nu], (vec[p, mu]*vec[p, nu])/sp[p, p]}}]}]
