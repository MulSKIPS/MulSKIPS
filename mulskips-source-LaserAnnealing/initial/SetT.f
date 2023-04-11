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
**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetT(b)  ! Called only if InitSt = 'LA'
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER x,y,z
      INTEGER :: IPF51,temperature
      
      REAL(8) :: LenIn(LenX*LenY)

      INTEGER, DIMENSION(LenZ,LenX*LenY) :: b
!      INTEGER, DIMENSION(30,LenX*LenY) :: b

      IF(tempfilename.EQ.'None')THEN ! array received from sockets

        write(*,*) 'b(1,1)', b(1,1)
        write(*,*)'Setting temperatures in LattCoo from array received via sockets'
        DO z=0,LenZ-1
!        DO z=0,30-1
          DO y=0,LenY-1
            DO x=0,LenX-1
              temperature = b(1+z,1+x+LenX*y)
              IF(temperature.NE.10)THEN  ! set T only if not wall 
                CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
              END IF
            END DO
          END DO
        END DO        
        write(*,*)'Done setting temperatures in LattCoo'

      ELSE
      
        write(*,*)'Reading input temperature file: ', tempfilename
        IPF51=51
        OPEN(IPF51,FILE=tempfilename,STATUS='OLD')  ! KELVIN
        DO z=0,LenZ-1
          READ(IPF51,*)LenIn
          DO y=0,LenY-1
            DO x=0,LenX-1
              temperature = INT(LenIn(1+x+LenX*y))
              IF(temperature.NE.10)THEN  ! set T only if not wall 
                CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
              END IF
            END DO
          END DO
        END DO
        CLOSE(IPF51)
        write(*,*)'Done reading input temperatures file'
      
      END IF
      

      END SUBROUTINE SetT

**************************************************************************************
