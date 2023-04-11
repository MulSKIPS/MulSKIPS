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
***************************************************************************
!     Absorption SUBROUTINE
!     it updates the systems's state after the absorption of an AdActiveAtom of the type CovInd (.ne.113)
!     the IndVoid points to a empty site which has one of its NN already filled
      SUBROUTINE Absorption(IndVoid,CovInd) ! H: CovInd=1
      USE DefDerType
      USE DefSystem
      USE Definitions

      IMPLICIT NONE
      INTEGER IndVoid,CovInd,IndNN,Index_Event,CovIndN
      INTEGER i,j,Coor,OccN,CoorN,SiOrC
      INTEGER Site(3),NextN(3,4),NextNN(3,4)
      REAL(8) Prob
      LOGICAL Recipr

      Site=ListAdAtom(IndVoid) % AtomXYZ
      Coor=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor)
      CALL MVBITS(CovInd,0,LenIndex,
     >            LattCoo(Site(1),Site(2),Site(3)),PosIndex)  ! Set coverage

      NextN = ListAdAtom(IndVoid) % NextNXYZ
      DO i=1,4
        OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
        IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
        CoorN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)     
        NextNN=0
        IF (OccN.NE.0)THEN
          CovIndN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosIndex,LenIndex)
          IF(CovIndN.NE.113)THEN
            IF(CoorN.LE.3)THEN 
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the site belongs to the NN of its NN
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3))) Recipr=.TRUE.
              END DO
              IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                SiOrC=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CovIndN,CoorN,NextN(1:3,i),      ! set evaporation probability (varies with coor and coverage)
     >                    NextNN,Index_Event,Prob)
                ListAdAtom(IndNN) % Ind_Event = Index_Event
                ListAdAtom(IndNN) % ProbTrans = Prob
                CALL Updatetree(IndNN,Prob)
              END IF
            END IF
          END IF
        END IF
      END DO

      CALL Get_Prob(0,0,CovInd,Coor,Site,NextN,Index_Event,Prob)        ! set desorption probability 
      ListAdAtom(IndVoid) % Ind_Event = Index_Event
      ListAdAtom(IndVoid) % ProbTrans = Prob
      CALL Updatetree(IndVoid,Prob)
      END SUBROUTINE Absorption
