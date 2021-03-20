# ExSCF
## Pre-requisites
* sympy
* PySCF
* libxc
* [MOKIT](https://gitlab.com/jxzou/mokit) (optional)

## Build
```
cd ExSCF/pyphf
cmake . -DCMAKE_INSTALL_PREFIX:PATH=$HOME/pyscf_deps # the directory above libxc's lib
make xc_itrf
```

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
