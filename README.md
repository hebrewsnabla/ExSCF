# ExSCF
## Pre-requisites
* sympy
* PySCF
* [MOKIT](https://gitlab.com/jxzou/mokit) (optional)

## Features
* SUHF
* SUHF+DFT, SUHF+*f*DFT
* guess: 
  + mix (`guess.mix`)
  + fragment (`guess.from_frag`)
* converging strategy
  + DIIS (`diis_on = True`)
  + `level_shift = ` (in a.u.)
* SUHF relaxed density, natural orbitals
* interface (require MOKIT)
  + read guess MO from fch
  + dump MO, NO to fch (`tofch = True`)

## TODO
* SU-PDFT
* SUPT2
