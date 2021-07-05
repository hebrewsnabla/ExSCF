# ExSCF
## Pre-requisites
* sympy
* PySCF (1.7.6 recommended)
* libxc
* [MOKIT](https://gitlab.com/jxzou/mokit) (optional, for read/write fch and CASDFT)
* [pyAutoMR](https://github.com/hebrewsnabla/pyAutoMR) (optional, for CASDFT)

Of course we need numpy, scipy, etc., but these are also required by PySCF, so they are not listed here.
libxc is optional in PySCF, but it's required here. Usually you've installed that when installing PySCF.

## Installation
* git clone
* add `/path/to/ExSCF` to your `PYTHONPATH`

## Features
**Theoretial features**:
* SUHF
  + energy
  + density, natural orbital

not fully tested:
* SUHF+DFT, SUHF+*f*DFT, SUHF+*fc*DFT
* CAS+*f*DFT, CAS-DFT2
* DeltaSCF with MOM

**Technical features**:
* guess
  + mix (`guess.mix`)
  + fragment (`guess.from_frag`)
  + stablize UHF (`guess.check_stab`)
* converging strategy
  + DIIS (`diis_on = True`)
  + `level_shift = ` (in a.u.)
* interface (require MOKIT)
  + read guess MO from fch (`guess.from_fch_simp`)
  + dump MO, NO to fch (`tofch = True`)

## Quick Start
```
from pyphf import suscf, guess

xyz = 'H 0.0 0.0 0.0; H 0.0 0.0 2.0'''
bas = '3-21g'
mf = guess.mix(xyz, bas, conv='tight')

mf2 = suscf.SUHF(mf)
mf2.kernel()
```

## TODO
* SU-PDFT
* SUPT2
* TD-SUHF
* SUHF gradient
