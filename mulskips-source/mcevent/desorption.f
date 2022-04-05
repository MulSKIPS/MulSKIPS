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
!     Desorption SUBROUTINE
!     it updates the systems's state after the desorption of a coverage atom
!     the IndAtom points to a full coverage site before the update procedure
      SUBROUTINE desorption(IndAtom) ! H: CovInd=1 
      USE DefDerType
      USE DefSystem
      USE Definitions

      IMPLICIT NONE
      INTEGER IndAtom,IndNN,Index_Event,CovIndN,CovInd
      INTEGER i,j,Coor,OccN,CoorN,SiOrC
      INTEGER Site(3),NextN(3,4),NextNN(3,4)
      REAL(8) Prob
      LOGICAL Recipr

      Site=ListAdAtom(IndAtom) % AtomXYZ

      CovInd=IBITS(LattCoo(Site(1),Site(2),Site(3)),
     >                  PosIndex,LenIndex)
      ! Reduce number of corresponding coverage atom
      IF(CovInd.EQ.ListCov(1))THEN
        CountCov(1) = CountCov(1) - 1
      ELSE IF(CovInd.EQ.ListCov(2))THEN
        CountCov(2) = CountCov(2) - 1
      ELSE IF(CovInd.EQ.ListCov(3))THEN
        CountCov(3) = CountCov(3) - 1
      END IF

      Coor=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor)
      CALL MVBITS(0,0,LenIndex,
     >            LattCoo(Site(1),Site(2),Site(3)),PosIndex)  ! Unset coverage

      IF(Coor.EQ.1)THEN ! there is only one NN to update. The others were not fixed! 
        NextN = ListAdAtom(IndAtom) % NextNXYZ 
        OccN = IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosOcc,LenOcc)
        IndNN = LattInd(NextN(1,1),NextN(2,1),NextN(3,1))
        CoorN=IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosCoor,LenCoor)
        NextNN=0
        IF (OccN.NE.0)THEN
          CovIndN=IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosIndex,LenIndex)
          IF(CovIndN.NE.113)THEN
            NextNN = ListAdAtom(IndNN) % NextNXYZ
            Recipr=.FALSE.
            DO j=1,4 ! check is the Site belongs to the NN of Next
              IF(ALL(NextNN(1:3,j).EQ.Site(1:3))) Recipr=.TRUE.
            END DO
            IF(Recipr)THEN ! action only if the site belongs to the NN of its NN
              SiOrC = IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosSiC,LenSiC)
              CALL Get_Prob(1,SiOrC,CovIndN,CoorN,NextN(1:3,1),NextNN,  ! set evaporation probability (varies with coor and coverage)
     >                     Index_Event,Prob)
              ListAdAtom(IndNN) % Ind_Event = Index_Event
              ListAdAtom(IndNN) % ProbTrans = Prob
              CALL Updatetree(IndNN,Prob)
            END IF
          END IF
        END IF
      ELSE
        NextN = ListAdAtom(IndAtom) % NextNXYZ 
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
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3))) Recipr=.TRUE.
              END DO
              IF(Recipr)THEN ! action only if the site belongs to the NN of its NN
                SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CovIndN,CoorN,NextN(1:3,i),NextNN,  ! set evaporation probability (varies with coor and coverage)
     >                     Index_Event,Prob)
                ListAdAtom(IndNN) % Ind_Event = Index_Event
                ListAdAtom(IndNN) % ProbTrans = Prob
                CALL Updatetree(IndNN,Prob)
              END IF
            END IF
          END IF
        END DO
      END IF
      CALL Get_Prob(0,0,0,Coor,Site,NextN,Index_Event,Prob)             ! set probability for deposition OR absorption
      ListAdAtom(IndAtom) % Ind_Event = Index_Event
      ListAdAtom(IndAtom) % ProbTrans = Prob
      CALL Updatetree(IndAtom,Prob)
      END SUBROUTINE desorption
