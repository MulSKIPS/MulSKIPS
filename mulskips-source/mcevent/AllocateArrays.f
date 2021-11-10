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
*******************************************************
*  Allocate Binary Tree for the Transition probability
*******************************************************
       SUBROUTINE AllocateArrays()
       USE DefSystem
       IMPLICIT NONE
       INTEGER ind

!!!!!!it calculates necessary number of levels !!!!!!!
       DO ind=0,31
         IF(BTEST(NumMcMax,ind))Levels=ind ! level necessary in the tree
       END DO
       IF(NumMcMax.GT.ISHFT(1,Levels))THEN
         Levels = Levels + 2
         SizeTree = ISHFT(1,Levels) - 1
       ELSE
         Levels = Levels + 1
         SizeTree = ISHFT(1,Levels) - 1
       END IF
       WRITE(*,*) 'numParticelle=',NumMcMax,
     >     'Levels =',Levels,', SizeTree =',SizeTree

       ALLOCATE(Tree(SizeTree))

       END SUBROUTINE AllocateArrays
******************************************************
