#!/usr/bin/env wolframscript
(* ---------------------------------------------------------------------------
   Regression test: contraction of Levi-Civita tensors from separate Dirac
   traces.

   Two independent closed Dirac traces each evaluate to an epsilon tensor
   eps(l1, p1, m, n) (via Tr[ (g.l1)(g.p1) g_m g_n g5 ]). The free Lorentz
   indices m, n are shared, so the product must contract to
       2 (sp[l1,l1] sp[p1,p1] - sp[l1,p1]^2)   (times the trace normalisation).

   In disentangled mode the two traces live in separate Lorentz structures and
   are only multiplied together after the contraction rounds. FORM's per-trace
   (Trace4/Tracen) epsilon contraction does NOT see the cross-structure product,
   and historically FormTracer issued no global `contract` afterwards, so the
   result was left as the un-contractible  16*(e_)[l1,p1,lor1,lor2]^2 , which
   then leaked into generated code as `e_` / Pattern[..] / Blank[..].

   The fix adds a single `contract 0;` (+ .sort) right before the
   "Rewriting scalar products and momentum labels" stage in generateFormFile,
   so the cross-trace e_*e_ is contracted and the resulting dot products flow
   through the existing FTxsp conversion.

   UNPATCHED FormTracer -> prints FAIL (leftover e_).
   PATCHED   FormTracer -> prints PASS (clean scalar products).
--------------------------------------------------------------------------- *)

Get["FormTracer`"];

(* adjust if FORM/TFORM is elsewhere; "-w16" is just a worker-count flag *)
DefineFormExecutable["/usr/bin/tform -w16"];

DefineLorentzDimensions[4];
DefineLorentzTensors[deltaLorentz[mu, nu], vec[p, mu], sp[p, q], eps[],
  deltaDirac[i, j], gamma[mu, i, j], gamma5[i, j]];
DisentangleLorentzStructures[True];

(* keep each trace on ONE line: a newline would terminate the implicit product *)
t1 = gamma[u1, i1, i2] vec[l1, u1] gamma[u2, i2, i3] vec[p1, u2] gamma[m, i3, i4] gamma[n, i4, i5] gamma5[i5, i1];
t2 = gamma[w1, j1, j2] vec[l1, w1] gamma[w2, j2, j3] vec[p1, w2] gamma[m, j3, j4] gamma[n, j4, j5] gamma5[j5, j1];

res = FormTrace[t1 t2];
Print["result = ", InputForm[res]];

leftoverEps = (! FreeQ[res, eps]) || StringContainsQ[ToString[res, InputForm], "e_"];

If[leftoverEps,
  Print["FAIL: leftover Levi-Civita (e_) -- FormTracer is UNPATCHED."],
  Print["PASS: epsilons contracted into scalar products (no e_)."]
];
