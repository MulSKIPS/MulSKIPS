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
***********************************************************************
       SUBROUTINE UpdateCoor_1(PosFrom,LenFrom,PosTo,LenTo,CoorNN)

        USE Definitions

        IMPLICIT NONE
        INTEGER :: PosFrom,LenFrom,PosTo,LenTo,CoorNN

        INTEGER    :: CoorFrom,CoorTo

        CoorTo   = IBITS(CoorNN,PosTo,LenTo) + 1
        CoorFrom = IBITS(CoorNN,PosFrom,LenFrom) - 1
        CALL MVBITS(CoorFrom,0,LenFrom,CoorNN,PosFrom)
        CALL MVBITS(CoorTo,0,LenTo,CoorNN,PosTo)
       END SUBROUTINE UpdateCoor_1
******|****************************************************************

***********************************************************************
       SUBROUTINE UpdateCoor_2(PosFrom,LenFrom,CoorNN)

        USE Definitions

        IMPLICIT NONE
        INTEGER :: PosFrom,LenFrom,CoorNN

        INTEGER    :: CoorFrom

        CoorFrom = IBITS(CoorNN,PosFrom,LenFrom) - 1
        CALL MVBITS(CoorFrom,0,LenFrom,CoorNN,PosFrom)
       END SUBROUTINE UpdateCoor_2
******|****************************************************************

***********************************************************************
       SUBROUTINE UpdateCoor_3(PosTo,LenTo,CoorNN)

        USE Definitions

        IMPLICIT NONE
        INTEGER :: PosTo,LenTo,CoorNN

        INTEGER   :: CoorTo

        CoorTo   = IBITS(CoorNN,PosTo,LenTo) + 1
        CALL MVBITS(CoorTo,0,LenTo,CoorNN,PosTo)
       END SUBROUTINE UpdateCoor_3
******|***************************************************************



***********************************************************************
       SUBROUTINE ExtrCoor(Site,Coor)

        USE DefSystem

        IMPLICIT NONE
        INTEGER :: Site,Coor(3)
        INTEGER :: res

        Coor(3) =Site/LenX/LenY
        res = Site - Coor(3)*LenX*LenY
        Coor(2) =res/LenX
	    Coor(1) = Res - Coor(2)*LenX

       END SUBROUTINE ExtrCoor
******|****************************************************************
