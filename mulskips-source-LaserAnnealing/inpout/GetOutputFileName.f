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
      SUBROUTINE GetOutputFileName(Time, Iter, FileNameBase)

      IMPLICIT    NONE
      INTEGER     Iter,T
      REAL(8)     Time
      CHARACTER(Len=9)   FileNameBase
      CHARACTER(LEN=1)   C1,C1A,C1B,C1C,C1D,C1E,C1F,C1G,C1H

      Iter=Iter*10

      IF (Iter.GE.1000000000)
     &  STOP 'GetOutputFileName:: Iter can not be represented'
C    Character representation offset: 48 !!!
      T = ITER
      IF (T.GE.100000000) THEN
        C1 = CHAR(INT(T/100000000) + 48)
        T = MOD(T,100000000)
      ELSE
        C1 = '0'
      ENDIF
      IF (T.GE.10000000) THEN
        C1A = CHAR(INT(T/10000000) + 48)
        T = MOD(T,10000000)
      ELSE
        C1A = '0'
      ENDIF
      IF (T.GE.1000000) THEN
        C1B = CHAR(INT(T/1000000) + 48)
        T = MOD(T,1000000)
      ELSE
        C1B = '0'
      ENDIF
      IF (T.GE.100000) THEN
        C1C = CHAR(INT(T/100000) + 48)
        T = MOD(T,100000)
      ELSE
        C1C = '0'
      ENDIF
      IF (T.GE.10000) THEN
        C1D = CHAR(INT(T/10000) + 48)
        T = MOD(T,10000)
      ELSE
        C1D = '0'
      ENDIF
      IF (T.GE.1000) THEN
        C1E = CHAR(INT(T/1000) + 48)
        T = MOD(T,1000)
      ELSE
        C1E = '0'
      ENDIF
      IF (T.GE.100) THEN
        C1F = CHAR(INT(T/100) + 48)
        T = MOD(T,100)
      ELSE
        C1F = '0'
      ENDIF
      C1G = CHAR(INT(T/10) + 48)
      C1H = CHAR(MOD(T,10) + 48)

      FileNameBase = 'I'//C1//C1A//C1B//C1C//C1D//C1E//C1F//C1G//C1H

      END SUBROUTINE GetOutputFileName
