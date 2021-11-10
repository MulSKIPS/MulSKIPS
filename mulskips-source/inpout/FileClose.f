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
      SUBROUTINE FileClose()

      USE Definitions

      IMPLICIT    NONE

      CLOSE(IPF)
      CLOSE(OPF0)
!      CLOSE(OPF1)
!      CLOSE(OPF2)
!      CLOSE(OPF3)
!      CLOSE(OPF4)
!      CLOSE(OPF5)
!      CLOSE(OPF6)
!      CLOSE(OPF7)
!      CLOSE(OPF8)
!      CLOSE(OPF9)
!      CLOSE(OPF10)
!      CLOSE(OPF11)
!      CLOSE(OPF12)
!      CLOSE(OPF13)
!      CLOSE(OPF14)
!      CLOSE(OPF15)

      END SUBROUTINE FileClose
******|****************************************************************
