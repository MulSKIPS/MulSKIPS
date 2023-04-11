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

      MODULE DefDerType



      TYPE BulkParticle
      INTEGER, DIMENSION(3)  :: AtomXYZ
      INTEGER, DIMENSION(3,4):: NextNXYZ
      END TYPE BulkParticle


      TYPE MCParticle
      INTEGER, DIMENSION(3)  :: AtomXYZ
      INTEGER, DIMENSION(3,4):: NextNXYZ
      INTEGER :: Ind_Event
      REAL(8) :: ProbTrans
      END TYPE MCParticle

      TYPE Vacancy

      INTEGER, DIMENSION(3)  :: AtomXYZ
      INTEGER, DIMENSION(3,4):: NextNXYZ
!	la dimensione dell'array è importante perche specifica
!	i possibili eventi inclusi nella cinetica di diffusione
!	bisogna aggiungere sempre +1 per il desorbimento
!	NUMTRANS=8 deve essere settato qui a mano, non in rif a Definitions

      END TYPE Vacancy

      TYPE AtBuffer
      INTEGER, DIMENSION(3)  :: AtomXYZ
      INTEGER :: Ind_Atype
      END TYPE AtBuffer


      TYPE Position
      INTEGER :: posX
      INTEGER :: posY
      INTEGER :: posZ
      END TYPE Position

      END MODULE DefDerType
