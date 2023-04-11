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
!     Evaporation SUBROUTINE
!     it updates the systems's state after the evaporation of an atom
!     the IndAtom points to a full site before the update procedure
      SUBROUTINE Evaporation(IndAtom) ! IAtom 1=01 Si 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions

      IMPLICIT NONE
      INTEGER IndAtom,IAtom,IndNN,IndNNN,Index_Event,CovInd,CovIndNN      ! WALL ! added CovInd and CovIndNN
      INTEGER i,j,jj,Coor,Occ,OccN,CoorNN,CoorNNN,SiOrC
      INTEGER Site(3),SiteC(3),NextN(3,4),NextNN(3,4),SNSite(3,3)
      INTEGER NNArm(3,3),NNZig(3,3),NBuff(3),NextNNN(3,4)
      REAL(8)Prob
      LOGICAL Recipr
      INTEGER temperature

      Site=ListAdAtom(IndAtom) % AtomXYZ
      IAtom=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosSiC,LenSiC)
      CountCrystal(IAtom) = CountCrystal(IAtom) - 1

      Coor=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor) ! not used
      !write(*,*)'Evaporation Site',Site,Coor
      CALL MVBITS(0,0,LenOcc,
     >            LattCoo(Site(1),Site(2),Site(3)),PosOcc)  ! No occupancy
      NumAdAtomOcc = NumAdAtomOcc - 1   ! Update Count of occupied AdAtoms

      CALL MVBITS(0,0,LenSiC,
     >            LattCoo(Site(1),Site(2),Site(3)),PosSiC)  ! No atom type

      NextN = ListAdAtom(IndAtom) % NextNXYZ ! All the NextN are fixed for
      DO i=1,4  ! The coordination of the next neighbour sites decreases of 1
        Occ = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
        CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
        IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
        CovInd=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),         ! WALL
     >                  PosIndex,LenIndex)                              ! WALL
        !write(*,*)"Ev. NextN,Occ,CoorNN",NextN(1:3,i),Occ,CoorNN
        NextNN=0
        IF (Occ.EQ.0)THEN  ! empty site
            IF (CoorNN.EQ.1) THEN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
              END DO
              IF(Recipr)THEN ! this empty site has not now any NN is must be deleted from the MC list
                  CALL EraseMC(NextN(1:3,i)) ! out of the MC particles
                  IF((InitSt.EQ.'LA').OR.(InitSt.EQ.'ST'))THEN
                    temperature = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >             NextN(3,i)),PosT,LenT) ! K
                    LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))=0
                    CALL MVBITS(temperature,0,LenT,LattCoo(NextN(1,i),
     >               NextN(2,i),NextN(3,i)),PosT)
                  ELSE  
                    LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))=0
                  END IF
              END IF
            ELSE IF (CoorNN.EQ.2) THEN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                 CoorNN=1
                 CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                 DO j=1,4
                   OccN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j),
     >                  NextNN(3,j)),PosOcc,LenOcc)
                   IF(OccN.EQ.0)THEN
                       ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                   ELSE
                     CoorNNN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j)
     >                  ,NextNN(3,j)),PosCoor,LenCoor)
                     CovIndNN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j)       ! WALL
     >                          ,NextNN(3,j)),PosIndex,LenIndex)          ! WALL 
                     IF((CoorNNN.LE.3).AND.(CovIndNN.NE.113))THEN         ! WALL NextNN could be a wall here (having OccN.NE.0 and CoorNNN.EQ.0). Go on looking for NextNNN only if NextNN is not a wall (otherwise it would give checkbounds error at line 92)
                        IndNNN=LattInd(NextNN(1,j),NextNN(2,j),
     >                                 NextNN(3,j))
                        NextNNN = ListAdAtom(IndNNN) % NextNXYZ
                        Recipr=.FALSE.
                        DO jj=1,4
                          IF(ALL(NextNNN(1:3,jj).EQ.NextN(1:3,i)))
     >                        Recipr=.TRUE.
                        END DO
                        IF(Recipr)THEN
                          NBuff(1:3)=
     >                        ListAdAtom(IndNN) % NextNXYZ(1:3,j)
                          ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                        ELSE
                          ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                        END IF
                     ELSE ! If CoorNNN=4 NextNN Cannot be NN of Next   ! Do this also if CoorNNN=0 and CovInd=113
                        ListAdAtom(IndNN) % NextNXYZ(1:3,j)=0
                     END IF
                   END IF
                 END DO
                 ListAdAtom(IndNN) % NextNXYZ(1:3,1)=NBuff(1:3)
                 CALL Get_Prob(0,0,CovInd,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ListAdAtom(IndNN) % Ind_Event = Index_Event
                 ListAdAtom(IndNN) % ProbTrans = Prob
                 CALL Updatetree(IndNN,Prob) !
              END IF
            ELSE IF (CoorNN.EQ.3) THEN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=2
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                CALL Get_Prob(0,0,CovInd,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                ListAdAtom(IndNN) % Ind_Event = Index_Event
                ListAdAtom(IndNN) % ProbTrans = Prob
                CALL Updatetree(IndNN,Prob) !
              END IF
            ELSE IF (CoorNN.EQ.4) THEN ! CoorNN.EQ.4 A previously Vacancy site now is a MC particle
              NextNN = ListVoid(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
C                 IF(Coor.EQ.0)THEN ! vacancy was trapped in vacuum/liquid
C                   ! eliminate vacancy and its neighbours
C                 ELSE                
                CoorNN=3
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                CALL Get_Prob(0,0,CovInd,CoorNN,NextN(1:3,i),NextNN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL EraseVoid(NextN(1:3,i))
                CALL AddMC(NextN(1:3,i),NextNN,Index_Event,Prob)
              END IF
            END IF
        ELSE
          IF(CovInd.NE.113)THEN                                         ! WALL:  Go on only if NextN is not a coverage site (note that it cannot be wall because Occ=0 here)
            IF(CovInd.NE.0)THEN
              write(*,*)'ERROR in deposition: Site can only
     >          can only NN with CovInd=0 or 113 when its Occ=1'
              STOP
            END IF
            IF (CoorNN.EQ.4)THEN! full site
              NextNN = ListAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=3
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CovInd,CoorNN,NextN(1:3,i),NextNN,     ! Set evaporation probability (can vary depending on coverage and coor) (Here CovInd can only be zero)
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL EraseBulk(NextN(1:3,i))
                CALL AddMC(NextN(1:3,i),NextNN,Index_Event,Prob)
                NumAdAtomOcc = NumAdAtomOcc + 1   ! Update Count of occupied AdAtoms
              END IF
            ELSE
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))
     >             Recipr=.TRUE.
              END DO
              IF(Recipr)THEN
                CoorNN=CoorNN-1
                CALL MVBITS(CoorNN,0,LenCoor,
     >              LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosCoor)
                SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CovInd,CoorNN,NextN(1:3,i),NextNN,     ! Set evaporation probability (can vary depending on coverage and coor) (Here CovInd can only be zero)
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                ListAdAtom(IndNN) % Ind_Event = Index_Event               ! NOTE: if CoorNN = 0, then Index_Event will be set to 3 with huge probability so it will evaporate immediately (Case N_Si+N_C+N_wall+N_cov=0 in Get_Prob). So no need to Erase it now. (Sicuri?)
                ListAdAtom(IndNN) % ProbTrans = Prob
                CALL Updatetree(IndNN,Prob) !
              END IF
            END IF
          ELSE                                                            ! WALL
            Coor = Coor - 1                                               ! WALL
          END IF
        END IF
      END DO
      IndAtom=LattInd(Site(1),Site(2),Site(3)) ! it can change due to the Next Neigh. update
      CALL MVBITS(Coor,0,LenCoor,LattCoo(Site(1),                         ! WALL ! Coor might have changed due to a CovInd.NE.0 
     >                        Site(2),Site(3)),PosCoor)
      IF(Coor.EQ.0)THEN  ! Isolated Atom Evaporation
         CALL EraseMC(Site)
      ELSE
        CALL Get_Prob(0,0,0,Coor,Site,NextN,Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C ! Here CovInd=0 by definition 
        ListAdAtom(IndAtom) % Ind_Event = Index_Event
        ListAdAtom(IndAtom) % ProbTrans = Prob
        IF(Coor.EQ.1)THEN ! the new empty site will not have a fixed coordination
          NextN = ListAdAtom(IndAtom) % NextNXYZ
          DO i=1,4
           OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
           IF(OccN.EQ.0)THEN
               ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
           ELSE ! Note that Coor
             CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
             CovIndNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),    ! WALL
     >                  PosIndex,LenIndex)                                ! WALL
             IF((CoorNN.LE.3).AND.(CovIndNN.NE.113))THEN                  ! WALL NextN could be a wall here (having OccN.NE.0 and CoorNN.EQ.0). Go on looking for NextNN only if NextN is not a wall (otherwise it would give checkbounds error at line 221)
               IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
               NextNN = ListAdAtom(IndNN) % NextNXYZ
               Recipr=.FALSE.
               DO j=1,4 ! check is the Site belongs to the NN of its NN
                 IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
               END DO
             !write(*,*)'Evap 2 Recipr',Recipr
               IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                 NBuff(1:3)=
     >             ListAdAtom(IndAtom) % NextNXYZ(1:3,i)
                 ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
               ELSE
                 ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
               END IF
             ELSE
               !write(*,*)'Evap Coordination = 4'
               ListAdAtom(IndAtom) % NextNXYZ(1:3,i)=0
             END IF
           END IF
          END DO
          ListAdAtom(IndAtom) % NextNXYZ(1:3,1)=NBuff(1:3)
        END IF
        CALL Updatetree(IndAtom,Prob)
      END IF
      END SUBROUTINE Evaporation
