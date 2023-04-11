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
      SUBROUTINE FileOpen()

      USE Definitions
      IMPLICIT    NONE

C   File descriptor 18 used in 'SUBROUTINE ReadDefectData'
      IPF = 19
      OPF0 = 20
      OPF1 = 21
      OPF2 = 22
      OPF3 = 23
      OPF4 = 24
      OPF5 = 25
      OPF6 = 26
      OPF7 = 27
      OPF8 = 28
      OPF9 = 29
      OPF10 = 30
      OPF11 = 31
      OPF12 = 32
      OPF13 = 33
      OPF14 = 34
      OPF15 = 35

      OPEN(IPF,FILE='start.dat',STATUS='OLD')
      OPEN(OPF0,FILE='Run_time.dat',STATUS='REPLACE')
!      OPEN(OPF1,FILE='check01.dat',STATUS='REPLACE')
!      OPEN(OPF2,FILE='check02.dat',STATUS='REPLACE')
!      OPEN(OPF3,FILE='check03.dat',STATUS='REPLACE')
!      OPEN(OPF4,FILE='check04.dat',STATUS='REPLACE')
!      OPEN(OPF5,FILE='check05.dat',STATUS='REPLACE')
!      OPEN(OPF6,FILE='check06.dat',STATUS='REPLACE')
!      OPEN(OPF7,FILE='check07.dat',STATUS='REPLACE')
!      OPEN(OPF8,FILE='check08.dat',STATUS='REPLACE')
!      OPEN(OPF9,FILE='check09.dat',STATUS='REPLACE')
!      OPEN(OPF10,FILE='check10.dat',STATUS='REPLACE')
!      OPEN(OPF11,FILE='check11.dat',STATUS='REPLACE')
!      OPEN(OPF12,FILE='check12.dat',STATUS='REPLACE')
!      OPEN(OPF13,FILE='check13.dat',STATUS='REPLACE')
!      OPEN(OPF14,FILE='check14.dat',STATUS='REPLACE')
!      OPEN(OPF15,FILE='check15.dat',STATUS='REPLACE')

      END SUBROUTINE FileOpen
******|****************************************************************
