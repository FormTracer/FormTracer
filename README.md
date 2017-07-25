# FormTracer
FormTracer is a high-performance, general purpose, easy-to-use Mathematica tracing package which uses [FORM](https://github.com/vermaseren/form). It supports arbitrary space
and spinor dimensions as well as an arbitrary number of simple compact Lie groups. While keeping the usability of the Mathematica interface, it
relies on the efficiency of FORM. An additional performance gain is achieved by a decomposition algorithm that avoids redundant traces in the product
tensors spaces. FormTracer supports a wide range of syntaxes which endows it with a high flexibility. Mathematica notebooks that automatically install
the package and guide the user through performing standard traces in spacetime, spinor and gauge-group spaces are provided.

FormTracer has been developed by Anton K. Cyrol, Mario Mitter, Jan M. Pawlowski, and Nils Strodthoff.
If used in scientific publications, please acknowledge our work by citing

> **A. K. Cyrol, M. Mitter, and N. Strodthoff, [Comput. Phys. Commun. 219C (2017) 346-352](https://doi.org/10.1016/j.cpc.2017.05.024), [arXiv:1610.09331 [hep-ph]](https://arxiv.org/abs/1610.09331)**

FormTracer is maintained by [selected members](https://github.com/orgs/FormTracer/people) of the fQCD collaboration.

## Features
* evaluation of (Euclidean) Lorentz/Dirac traces in arbitrary dimensions and traces over an arbitrary number of group product spaces
* intuitive, easy-to-use and highly customizable Mathematica frontend
* high performance due to FORM backend combined with an efficient decomposition algorithm in Mathematica
* supports
  * a special time-like direction for (Euclidean) finite temperature and density applications
  * partial traces involving open indices
  * creation of optimized output (including bracketing) using FORM’s optimization algorithm for further numerical processing in C/C++/Fortran
  * user-defined combined Lorentz tensors and corresponding identities, e.g. (transverse and longitudinal) projectors and their orthogonality relations, for speedup
  * convenient installation and update procedure within Mathematica

## Build instructions
FormTracer is equipped with a fully automated installation script in Mathematica, which can be downloaded and started by evaluating<br>
`Import["https://raw.githubusercontent.com/FormTracer/FormTracer/master/src/FormTracerInstaller.m"]`<br>
in a Mathematica input cell. Example notebooks are available for download:<br>
* [FormTracerShowcase](https://raw.githubusercontent.com/FormTracer/FormTracer/master/src/Examples/FormTracerShowcase.nb)
* [FormTracerMinimalExample](https://raw.githubusercontent.com/FormTracer/FormTracer/master/src/Examples/FormTracerMinimalExample.nb)
* [FourQuarkInteraction](https://raw.githubusercontent.com/FormTracer/FormTracer/master/src/Examples/FourQuarkInteraction.nb)

Both can also be used to install FormTracer. For further information we refer to  [arXiv:1610.09331 [hep-ph]](https://arxiv.org/abs/1610.09331).

## Bug reports
Please reports bugs via the [issue tracker](https://github.com/FormTracer/FormTracer/issues) on github.

## Links to related software
* [FORM](http://www.nikhef.nl/~form/) - project for symbolic manipulation of very big expressions
* [DoFun](http://physik.uni-graz.at/~mqh/DoFun/) - Derivation of Functional Equations
* [FeynCalc](https://feyncalc.github.io/) - Tools and Tables for QFT Calculations 

