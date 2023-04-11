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
******************************************************************************
      SUBROUTINE InterpProbFromField(temperature)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER temperature, ilower, iupper, wlower, wupper
      INTEGER :: diff(ntemp)

      diff = reftempvalues(:ntemp) - temperature
      IF(ANY(diff.EQ.0))THEN
        iupper = FINDLOC(diff, 0, DIM=1)
        PtransE = refPtransE(iupper,:,:,:,:,:,:,:)
        
        PtransD = refPtransD(iupper,:,:)
        
        PtransAbs = refPtransAbs(iupper,:,:)
        PtransDes = refPtransDes(iupper,:,:)
      ELSE
        iupper = MINLOC(diff, DIM=1, MASK=(diff.GT.0))
        ilower = iupper-1 ! assuming that T is monotonously increasing
        wlower = ABS(diff(iupper)) ! lower bound weights more if upper bound is far from temperature
        wupper = ABS(diff(ilower))
        PtransE = (wlower*refPtransE(ilower,:,:,:,:,:,:,:) + 
     >             wupper*refPtransE(iupper,:,:,:,:,:,:,:) ) 
     >              / (wlower+wupper)
        PtransD = (wlower*refPtransD(ilower,:,:) + 
     >             wupper*refPtransD(iupper,:,:) ) 
     >              / (wlower+wupper)
        PtransAbs = (wlower*refPtransAbs(ilower,:,:) + 
     >               wupper*refPtransAbs(iupper,:,:) ) 
     >              / (wlower+wupper)
        PtransDes = (wlower*refPtransDes(ilower,:,:) + 
     >               wupper*refPtransDes(iupper,:,:) ) 
     >              / (wlower+wupper)
      END IF
      
      END SUBROUTINE InterpProbFromField
