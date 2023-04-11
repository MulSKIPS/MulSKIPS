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
      FUNCTION IF_Arm_Choice(SiteOld,SiteCOld,NZig,NArm)
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE

      INTEGER i,j,IndArm,IndZig,NumArm,NumZig,CoorArm,CoorZig
      INTEGER NArm(3,3),NZig(3,3)
      INTEGER NNArm(3,3),NNZig(3,3),SNSite(3,3)
      INTEGER SiteOld(3),SiteCOld(3),Site(3),SiteC(3)
      LOGICAL IF_Arm_Choice
      REAL(8) random
      EXTERNAL random

      CoorArm=0
      CoorZig=0
      SiteC=SiteOld
      DO i=1,3
       CoorArm=CoorArm+IBITS(LattCoo(NArm(1,i)
     >          ,NArm(2,i),NArm(3,i)),PosCoor,LenCoor) ! ???
       CoorZig=CoorZig+IBITS(LattCoo(NZig(1,i)
     >          ,NZig(2,i),NZig(3,i)),PosCoor,LenCoor) ! ???
      END DO
      IF(CoorZig.GT.CoorArm)THEN
         IF_Arm_Choice=.FALSE.
      ELSE IF(CoorZig.EQ.CoorArm)THEN
         IF(random(idum).le.PTransZig)THEN
           IF_Arm_Choice=.FALSE.
         ELSE
           IF_Arm_Choice=.TRUE.
         END IF
      ELSE
         IF_Arm_Choice=.TRUE.
      END IF
      !write(*,*)'IF_Arm_choice',CoorZig,CoorArm,IF_Arm_Choice
      END FUNCTION IF_Arm_Choice
