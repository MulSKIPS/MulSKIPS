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
********************************************************************
*  It adds a MC particle to the lists and in the lattice
********************************************************************
       SUBROUTINE AddMC(Site,NextN,Index_Event,Prob)
!    Questo metodo aggiorna l'array delle cinetiche: AtprobList
!    inserisce la particella nell'albero e ne aggiorna le probabilita
!    MODIFICA: AtprobList, Tree (attraverso UpdateTree), LattInd, Latt

	USE DefDerType
	USE DefSystem
	IMPLICIT NONE
	Real*8 :: Prob
	INTEGER Index_Event
	INTEGER, DIMENSION(3) :: Site
	INTEGER, DIMENSION(3,4) :: NextN

	IF(Prob.EQ.0) THEN
           WRITE(*,*)"provo ad aggiungere prob nulla"
           STOP
        END IF
	NumAdAtom = NumAdAtom + 1
	IF(NumAdAtom.GE.NumMCMax) then
	  write(*,*) NumAdAtom,NumMCMax
	  STOP "troppe particelle cinetiche"
	endif
	LattInd(Site(1),Site(2),Site(3)) = NumAdAtom
	ListAdAtom(NumAdAtom) % AtomXYZ = Site
	ListAdAtom(NumAdAtom) % NextNXYZ = NextN
	ListAdAtom(NumAdAtom) % Ind_Event = Index_Event
	ListAdAtom(NumAdAtom) % ProbTrans = Prob
	CALL UpdateTree(NumAdAtom,Prob)
	END SUBROUTINE AddMC

********************************************************************
*  It adds a bulk particle to the lists and in the lattice
********************************************************************
       SUBROUTINE AddBulk(Site,NextN)
!    Questo metodo aggiorna l'array delle cinetiche: AtprobList
!    inserisce la particella nell'albero e ne aggiorna le probabilita
!    MODIFICA: AtprobList, Tree (attraverso UpdateTree), LattInd, Latt

	USE DefDerType
	USE DefSystem
	IMPLICIT NONE
	INTEGER, DIMENSION(3) :: Site
	INTEGER, DIMENSION(3,4) :: NextN

	NumAtoms = NumAtoms + 1
	IF(NumAtoms.GE.NumAtMax) then
	  write(*,*) NumAtoms,NumAtMax
	  STOP "troppe particelle bulk"
	endif
	LattInd(Site(1),Site(2),Site(3)) = NumAtoms
	ListAtom(NumAtoms) % AtomXYZ = Site
	ListAtom(NumAtoms) % NextNXYZ = NextN
	END SUBROUTINE AddBulk


********************************************************************
*  It adds a Vacancy to the lists and in the lattice
********************************************************************
       SUBROUTINE AddVoid(Site,NextN)
!    Questo metodo aggiorna l'array delle cinetiche: AtprobList
!    inserisce la particella nell'albero e ne aggiorna le probabilita
!    MODIFICA: AtprobList, Tree (attraverso UpdateTree), LattInd, Latt

	USE DefDerType
	USE DefSystem
	IMPLICIT NONE
	INTEGER, DIMENSION(3) :: Site
	INTEGER, DIMENSION(3,4) :: NextN

	NumVoids = NumVoids + 1
	IF(NumVoids.GE.NumVdMax) then
	  write(*,*) NumVoids,NumVdMax
	  STOP "troppe vacanze"
	endif
	LattInd(Site(1),Site(2),Site(3)) = NumVoids
	ListVoid(NumVoids) % AtomXYZ = Site
	ListVoid(NumVoids) % NextNXYZ = NextN
	END SUBROUTINE AddVoid
