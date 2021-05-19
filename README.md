# ExSCF
## Pre-requisites
* sympy
* PySCF (1.7.6 recommended)
* libxc
* [MOKIT](https://gitlab.com/jxzou/mokit) (optional, for read/write fch and auto CAS calculation)

## Installation
* git clone
* add `/path/to/ExSCF` to your `PYTHONPATH`
* build(deprecated, only if you need modified DFT functional)
```
cd ExSCF/pyphf
cmake . -DCMAKE_INSTALL_PREFIX:PATH=$HOME/pyscf_deps # the directory above libxc's lib
make xc_itrf_mod
```

## Features
* SUHF
* SUHF+DFT, SUHF+*f*DFT, SUHF+*fc*DFT
* guess: 
  + mix (`guess.mix`)
  + fragment (`guess.from_frag`)
  + stablize (`guess.check_stab`)
* converging strategy
  + DIIS (`diis_on = True`)
  + `level_shift = ` (in a.u.)
* SUHF relaxed density, natural orbitals
* interface (require MOKIT)
  + read guess MO from fch (`guess.from_fch_simp`)
  + dump MO, NO to fch (`tofch = True`)

## TODO
* SU-PDFT
* SUPT2
* TD-SUHF
