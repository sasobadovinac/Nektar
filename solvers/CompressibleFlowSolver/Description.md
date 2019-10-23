Nektar++/feature/scalar-transport
========================
## Short description
Compressible flow solver with scalar transport equation.
Supported Riemann solvers are
- AUSM0 (not tested)
- AUSM1 (not tested)
- AUSM2 (tested)
- AUSM3 (fails even without scalar)
- Average (not tested)
- HLL (not tested)
- HLLC (tested)
- LaxFriedrichs (tested)
(Roe and Exact Toro are not supported)

To enable scalar transport equation simply add a variable(s) to xml file and
place them after energy (E). For more details see folder Examples.

