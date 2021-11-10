!> @file
!!
!!
!!   Copyright (C) 2019-2020
!!   @authors: Ioannis Deretzis, Giuseppe Fisicaro and Antonino La Magna
!!   This file is part of the MulSKIPS code.
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/licenses/gpl-3.0.txt .
!!   For the list of contributors, see ~/AUTHORS

!!   MulSKIPS is a free software: you can redistribute it and/or modify
!!   it under the terms of the GNU General Public License as published by
!!   the Free Software Foundation, either version 3 of the License, or
!!   (at your option) any later version.

!!   MulSKIPS is distributed in the hope that it will be useful,
!!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!   GNU General Public License for more details
      MODULE Definitions

       IMPLICIT NONE

! Random number generation: Initialization variable for Ran3
       INTEGER :: Idum

       REAL(8) :: sf ! [nm] superlattice parameter, to convert KMC units

! Energy
       REAL(8) :: Temp,Beta
       REAL(8) :: Tm, P0, sigma, Tflex, LenVac, LenNuc
       INTEGER :: tmax

       INTEGER     IPF,OPF0,OPF1,OPF2,OPF3,OPF4,OPF5,OPF6,OPF7,OPF8,OPF9
       INTEGER     OPF10,OPF11,OPF12,OPF13,OPF14,OPF15,counter
       CHARACTER(Len=2)  :: InitSt,RunType,SaveFinalState,SaveCoo
       CHARACTER(Len=2)  :: Homogeneous ! T or F (for LA)
       CHARACTER(Len=4)  :: ExitStrategy

       CHARACTER(Len=200) :: cadfilename
       CHARACTER(Len=200) :: tempfilename
       CHARACTER(Len=200) :: coofilename
       CHARACTER(Len=200) :: restartfilename

      END MODULE Definitions

