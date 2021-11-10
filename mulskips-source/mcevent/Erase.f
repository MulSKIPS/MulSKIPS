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
*    Erase a particle from the MC list
************************************************************************
      SUBROUTINE EraseMC(Site)
!	It modifies:
!	AtProbList(LattInd(x,y,z)), LattInd(x,y,z)


      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Site(3),Itype
      INTEGER, DIMENSION(3)   :: SiteTemp
      INTEGER, DIMENSION(3,4) :: NextNTemp
      INTEGER x,y,z,IndNN
      INTEGER n,ind
      x=Site(1)
      y=Site(2)
      z=Site(3)
      IF(z.eq.-1)THEN
       !CALL ERR("accesso a latt fuori bond")
       write(*,*)"accesso a latt fuori bond"
       STOP
      END IF
      indNN=LattInd(x,y,z)
      IF(IndNN.GT.NumMcMax.or.IndNN.EQ.0) THEN
        write(*,*) "err: EraseMC",IndNN,x,y,z
        STOP
      ENDIF

      IF(IndNN.GE.NumAdAtom)THEN
            ListAdAtom(IndNN) % AtomXYZ = 0
            ListAdAtom(IndNN) % NextNXYZ = 0
            ListAdAtom(IndNN) % ProbTrans = 0.d0
            ListAdAtom(IndNN) % Ind_Event = 0
            LattInd(Site(1),Site(2),Site(3))=0
            CALL UpdateTree(IndNN,0.d0)
            NumAdAtom = NumAdAtom - 1
      ELSE
            ListAdAtom(IndNN) = ListAdAtom(NumAdAtom)
       	    ind = NumAdAtom + ISHFT(1,Levels-1)-1
       	    CALL UpdateTree(IndNN,Tree(ind))
       	    CALL UpdateTree(NumAdAtom,0.d0)
       	    SiteTemp =  ListAdAtom(NumAdAtom) % AtomXYZ
!            NextNTemp = AtprobList(NumAdAtom) % NextNXYZ
       	    LattInd(SiteTemp(1),SiteTemp(2),SiteTemp(3)) = IndNN
            ListAdAtom(NumAdAtom) % ProbTrans = 0.d0
            ListAdAtom(NumAdAtom) % AtomXYZ = 0
            ListAdAtom(NumAdAtom) % NextNXYZ = 0
            ListAdAtom(NumAdAtom) % Ind_Event = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumAdAtom = NumAdAtom - 1
      END IF
c	      IF(NumMc.EQ.2) STOP "sono finite le particelle cinetiche...."
      END SUBROUTINE EraseMC
***********************************************************************

***********************************************************************
*    Erase a particle from the bulk particle list
************************************************************************
      SUBROUTINE EraseBulk(Site)
!	It modifies:
!	AtProbList(LattInd(x,y,z)), LattInd(x,y,z)


      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Site(3),Itype
      INTEGER, DIMENSION(3)   :: SiteTemp
      INTEGER, DIMENSION(3,4) :: NextNTemp
      INTEGER x,y,z,IndNN
      INTEGER n,ind
      x=Site(1)
      y=Site(2)
      z=Site(3)
      IF(z.eq.-1)THEN
       !CALL ERR("accesso a latt fuori bond")
       write(*,*)"accesso a latt fuori bond"
       STOP
      END IF
      indNN=LattInd(x,y,z)
      IF(IndNN.GT.NumAtMax.or.IndNN.EQ.0) THEN
        write(*,*) "err: EraseBulk",IndNN,x,y,z
        STOP
      ENDIF
      !write(*,*)'erase bulk indnn site',indnn,site
      IF(IndNN.EQ.NumAtoms)THEN
            ListAtom(NumAtoms) % AtomXYZ = 0
            ListAtom(NumAtoms) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumAtoms = NumAtoms - 1
      ELSE
            ListAtom(IndNN) = ListAtom(NumAtoms)
       	    SiteTemp =  ListAtom(NumAtoms) % AtomXYZ
       	    !write(*,*)SiteTemp
!            NextNTemp = AtprobList(NumAdAtom) % NextNXYZ
       	    LattInd(SiteTemp(1),SiteTemp(2),SiteTemp(3)) = IndNN
            ListAtom(NumAtoms) % AtomXYZ = 0
            ListAtom(NumAtoms) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumAtoms = NumAtoms - 1
      END IF
c	      IF(NumMc.EQ.2) STOP "sono finite le particelle cinetiche...."
      END SUBROUTINE EraseBulk
***********************************************************************


***********************************************************************
*    Erase a particle from the bulk particle list
************************************************************************
      SUBROUTINE EraseVoid(Site)
!	It modifies:
!	AtProbList(LattInd(x,y,z)), LattInd(x,y,z)


      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Site(3),Itype
      INTEGER, DIMENSION(3)   :: SiteTemp
      INTEGER, DIMENSION(3,4) :: NextNTemp
      INTEGER x,y,z,IndNN
      INTEGER n,ind
      x=Site(1)
      y=Site(2)
      z=Site(3)
      IF(z.eq.-1)THEN
       !CALL ERR("accesso a latt fuori bond")
       write(*,*)"accesso a latt fuori bond"
       STOP
      END IF
      indNN=LattInd(x,y,z)
      IF(IndNN.GT.NumVdMax.or.IndNN.EQ.0) THEN
        write(*,*) "err: EraseVoid",IndNN,x,y,z
        STOP
      ENDIF
      !write(*,*)'erase void indnn site',indnn,site
      IF(IndNN.EQ.NumVoids)THEN
            ListVoid(NumVoids) % AtomXYZ = 0
            ListVoid(NumVoids) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumVoids = NumVoids - 1
      ELSE
            ListVoid(IndNN) = ListVoid(NumVoids)
       	    SiteTemp =  ListVoid(NumVoids) % AtomXYZ
       	    !write(*,*)SiteTemp
!            NextNTemp = AtprobList(NumAdAtom) % NextNXYZ
       	    LattInd(SiteTemp(1),SiteTemp(2),SiteTemp(3)) = IndNN
            ListVoid(NumVoids) % AtomXYZ = 0
            ListVoid(NumVoids) % NextNXYZ = 0
            LattInd(Site(1),Site(2),Site(3))=0
            NumVoids = NumVoids - 1
      END IF
c	      IF(NumMc.EQ.2) STOP "sono finite le particelle cinetiche...."
      END SUBROUTINE EraseVoid
***********************************************************************
