# ExSCF

## Features
* SUHF
* SUHF+DFT
* guess: 
  + mix (read from Gaussian)
  + fragment (`frag.guess_frag`)
* converging strategy
  + DIIS (`diis_on = True`)
  + `level_shift = ` (in a.u.)
* SUHF relaxed density, natural orbitals
* interface
  + read guess MO from fch
  + dump MO, NO to fch (`tofch = True`)

