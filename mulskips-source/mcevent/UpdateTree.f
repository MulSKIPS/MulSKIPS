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
**********************************************************************
*     It modifies to Patom the probability related to the particle
*     and recursively update the tree
***********************************************************************
      SUBROUTINE UpdateTree(Index,Patom)

!       modifica la particella Index, assegnadogli la probabilita
!       Patom, ATTENTO che Patom rappresenta prob totale
!       aggiorna poi l'albero di conseguenza

	USE DefSystem
	USE Definitions
	IMPLICIT NONE
	INTEGER :: Index
	Real(8)  Patom
	INTEGER Pos,Parent,NewPos,NeighPos,ControlBit

	IF(Patom.EQ.0.and.Index.lt.NumAdAtom) THEN
		write(*,*) Index,ListAdatom(Index) % AtomXYZ,Patom
		!CALL ERR("stabilizzato  da dentro Update!")
		write(*,*)"error stabilizzato  da dentro Update!"
		STOP
	ELSE if (Index.gt.NumAdAtom)THEN
		write(*,*) Index,NumAdAtom
		!CALL ERR("stabilizzato  da dentro Update!")
		write(*,*)"Index > NumAdAtom"
                STOP
	ENDIF
	Pos = Index + ISHFT(1,Levels-1)-1
	Tree(Pos) = Patom
	Parent=Levels-1
	DO
	   IF (Parent.LT.1) EXIT
	   ControlBit = IAND(Pos,1)
	   SELECT CASE (ControlBit)
	   CASE (0)
	      NeighPos = Pos + 1
	   CASE (1)
	      NeighPos = Pos - 1
	   END SELECT
	   NewPos = ISHFT(Pos,-1)
	   Tree(NewPos) = Tree(Pos) + Tree(NeighPos)
	   Pos = NewPos
	   Parent = Parent - 1
	END DO
       END SUBROUTINE UpdateTree
