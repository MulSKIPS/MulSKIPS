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
!     Deposition SUBROUTINE
!     it updates the systems's state after the deposition of an atom of the type IAtom
!     the IndVoid points to a empty site which have one of its NN already filled
      SUBROUTINE Deposition(IndVoid,IAtom) ! IAtom 1=01 Si 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions

      IMPLICIT NONE
      INTEGER IndVoid,IAtom,IndNN,IndNNN,Index_Event,CovInd,CovIndNNN     ! WALL added CovInd and CovIndNNN 
      INTEGER DisNN2,Dx,Dy
      INTEGER i,j,jj,Coor,CoorOld,Occ,OccN,CoorNN,CoorNNN,SiOrC
      INTEGER Site(3),SiteC(3),NextN(3,4),NextNN(3,4),NextNNN(3,4)
      INTEGER Newsite(3,2),NNArm(3,3),NNZig(3,3),SNSite(3,3)
      REAL(8):: Prob
      LOGICAL IF_Arm,IF_Arm_Choice,Recipr
      EXTERNAL IF_Arm_Choice

      Site=ListAdAtom(IndVoid) % AtomXYZ
      !write(*,*)'Deposition Site',Site
      Coor=IBITS(LattCoo(Site(1),Site(2),Site(3)),PosCoor,LenCoor) ! Get the previous coordination
      CoorOld=Coor
      CALL MVBITS(1,0,LenOcc,
     >            LattCoo(Site(1),Site(2),Site(3)),PosOcc) ! Set occupancy
      CALL MVBITS(IAtom,0,LenSiC,
     >            LattCoo(Site(1),Site(2),Site(3)),PosSiC)  ! Set atom kind
      !write(*,*)" "
      !write(*,*)"Deposition Coor IndVoid",Coor
      IF(Coor.EQ.1)THEN ! the new full site has not get yet a fixed coordination
        SiteC = ListAdAtom(IndVoid) % NextNXYZ(1:3,1)
        IndNN = LattInd(SiteC(1),SiteC(2),SiteC(3))
        CoorNN=IBITS(LattCoo(SiteC(1),SiteC(2),SiteC(3)),
     >                  PosCoor,LenCoor)  ! not used
        SNSite=0
        j=0
        DO i=1,4
         IF(.NOT.ALL(ListAdAtom(IndNN) % NextNXYZ(1:3,i).EQ.
     >      Site(1:3)))THEN
            j=j+1  ! should reaches the value = 3
            SNSite(1:3,j)= ListAdAtom(IndNN) % NextNXYZ(1:3,i)
         END IF
        END DO
        !write(*,*)'Dep. Call SN Routine Site,SiteC',Site,SiteC
        !write(*,*)SNSite,IndNN
        CALL FIND_SN(Site,SiteC,SNSite,NNZig,NNArm,LenX,LenY)           ! This could also find coverage neighbours...
        !write(*,*)NNZiG
        !write(*,*)NNArm
        IF_Arm = IF_Arm_Choice(Site,SiteC,NNZig,NNArm)  ! Fix the coordination type and store the NN positions
        IF(IF_Arm)THEN
          ListAdAtom(IndVoid) % NextNXYZ(1:3,2:4)=NNArm
        ELSE
          ListAdAtom(IndVoid) % NextNXYZ(1:3,2:4)=NNZig
        END IF
      END IF

      NextN = ListAdAtom(IndVoid) % NextNXYZ
      !write(*,*)"Deposition NextN",NextN
      DO i=1,4
        Occ = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosOcc,LenOcc)
        IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
        CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosCoor,LenCoor)
        CovInd=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),         ! WALL
     >                  PosIndex,LenIndex)                              ! WALL
        !write(*,*)"Deposition IndNN,CoorNN,Occ",IndNN,CoorNN,Occ
        NextNN=0
        IF (Occ.EQ.0)THEN  ! empty site
            IF (CoorNN.EQ.0) THEN  ! This empty site is a new MC particle its Coordination is now 1 the type is not fixed
              CoorNN=1
              CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
              NextNN(1:3,1)=Site ! Only one NN is fixed
              CALL Get_Prob(0,0,CovInd,CoorNN,NextN(1:3,i),NextNN,
     >                      Index_Event,Prob) ! Occ=1/0, NSiC 1=01 Si, 2=10 C, If Occ=0 NSiC=0
              IF(Index_Event.NE.0)THEN
                 !write(*,*)"Dep. Coor1 Ind_Ev,Prob",Index_Event,Prob
                 CALL AddMC(NextN(1:3,i),NextNN,Index_Event,Prob)
              END IF
            ELSE IF (CoorNN.EQ.1) THEN ! This empty site now fixes its next neighbours if a tetraedral configuration can be find
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Dx=ABS(NextNN(1,1)-Site(1))
              Dx=MIN(Dx,LenX-Dx)
              Dy=ABS(NextNN(2,1)-Site(2))
              Dy=MIN(Dy,LenY-Dy)
              DisNN2=Dx*Dx+Dy*Dy+(NextNN(3,1)-Site(3))**2
              IF(DisNN2.EQ.72)THEN  ! 72 is the second neighbour distance in tetrahedral configuration
                CoorNN = 2
                NextNN(1:3,2) = Site
                CALL FIND_NN(NextNN(1:3,1),NextNN(1:3,2),
     >                     NextN(1:3,i),NewSite,LenX,LenY) ! find remanent nn
                DO j=3,4
                   NextNN(1:3,j)=NewSite(1:3,j-2)
                   OccN = IBITS(LattCoo(NextNN(1,j),NextNN(2,j),
     >                          NextNN(3,j)),PosOcc,LenOcc)
                   IF(OccN.NE.0)THEN
                     !write(*,*)'NN occupied',NextN(1:3,i)
                     !write(*,*)'NN occupied',NextNN(1:3,j)
                     CoorNNN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j)
     >                  ,NextNN(3,j)),PosCoor,LenCoor)
                     CovIndNNN=IBITS(LattCoo(NextNN(1,j),NextNN(2,j)      ! WALL 
     >                     ,NextNN(3,j)),PosIndex,LenIndex)               ! WALL 
                     IF((CoorNNN.LE.3).AND.(CovIndNNN.NE.113))THEN        ! WALL NextNN could be a wall here (having OccN.NE.0 and CoorNNN.EQ.0). Go on looking for NextNNN only if NextNN is not a wall (otherwise it would give checkbounds error at line 124)
                        IndNNN=LattInd(NextNN(1,j),NextNN(2,j),
     >                                   NextNN(3,j))
                        NextNNN = ListAdAtom(IndNNN) % NextNXYZ
                        Recipr=.FALSE.
                        DO jj=1,4 ! check is the site belongs to the NN of its NN
                           IF(ALL(NextNNN(1:3,jj).EQ.NextN(1:3,i)))
     >                          Recipr=.TRUE.
                        END DO
                        IF(Recipr)THEN
                          CoorNN = CoorNN+1
                        END IF
                     ELSE
                        Recipr=.FALSE.
                     END IF
                     !write(*,*)Recipr
                   END IF
                END DO
                CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
                ListAdAtom(IndNN) % NextNXYZ = NextNN
                CALL Get_Prob(0,0,CovInd,CoorNN,NextN(1:3,i),NextNN,
     >                      Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C

                IF(Index_Event.EQ.0)THEN
                  CALL EraseMC(NextN(1:3,i))
                  CALL AddVoid(NextN(1:3,i),NextNN)
                ELSE ! 72 is the NN in tetrahedral configuration Site can belong to NextNN
                  ListAdAtom(IndNN) % Ind_Event = Index_Event
                  ListAdAtom(IndNN) % ProbTrans = Prob
                  CALL Updatetree(IndNN,Prob) !
                END IF
              !ELSE ! distance between NextNN and Site is not seconf NN one Site cannot be added to NextNN and CoorNN=1
                  !write(*,*)'DisNN2',DisNN2
                  !write(*,*)NextNN(1:3,1)
                  !write(*,*)Site
              !    STOP
              END IF
            ELSE IF(CoorNN.EQ.2.OR.CoorNN.EQ.3)THEN ! This empty site now increases the coordination if Site is in its NN
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the site belongs to the NN of its NN
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
              END DO
              !write(*,*)'Deposition Recipr',Recipr
              IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                CoorNN = CoorNN + 1
                CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
C                 IF((CovInd.GT.0).AND.(CoorNN.EQ.4)) at this point, get_prob should put it in ListVoid! It bahaves as a vacancy 
                CALL Get_Prob(0,0,CovInd,CoorNN,NextN(1:3,i),NextNN,
     >                      Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C

                IF(Index_Event.EQ.0)THEN
                  !write(*,*)"Dep. Vacancy Generation",NextN(1:3,i)
                  CALL EraseMC(NextN(1:3,i))
                  CALL AddVoid(NextN(1:3,i),NextNN)
                ELSE ! 72 is the NN in tetrahedral configuration Site can belong to NextNN
                  ListAdAtom(IndNN) % Ind_Event = Index_Event
                  ListAdAtom(IndNN) % ProbTrans = Prob
                  !write(*,*)"Dep.C3 Index_Event,Prob",Index_Event,Prob
                  CALL Updatetree(IndNN,Prob) !
                END IF
              ELSE ! Site in not in the NN list of Next Coord cannot increase
                !write(*,*)NextNN
                !write(*,*)Site
               ! STOP
              END IF
            END IF
        ELSE ! full site
          IF(CovInd.NE.113)THEN                                         ! Go on only if NextN is not a wall site, otherwise increase Coor (note that it cannot be coverage, because Occ=1 here)
            IF(CovInd.NE.0)THEN
              write(*,*)'ERROR in deposition: Site can only
     >          can only have NN with CovInd=0 or 113 when its Occ=1'
              STOP
            END IF
            IF(CoorNN.LE.3)THEN ! Only if belong to MC particles
              NextNN = ListAdAtom(IndNN) % NextNXYZ
              Recipr=.FALSE.
              DO j=1,4 ! check is the site belongs to the NN of its NN
                IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
              END DO
            !write(*,*)'Deposition Recipr',Recipr
              IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                CoorNN=CoorNN + 1
                CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i)
     >              ,NextN(2,i),NextN(3,i)),PosCoor)
                SiOrC=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
                CALL Get_Prob(1,SiOrC,CovInd,CoorNN,NextN(1:3,i),NextNN,     ! Set evaporation probability (can vary depending on coverage and coor) (Here CovInd can only be zero!!!)
     >                    Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                
                IF(Index_Event.EQ.0)THEN                                ! Set Bulk also If CoorNN reaches 4 
                !write(*,*)"Dep. Bulk Site Generation",Index_Event,Prob
                  CALL EraseMC(NextN(1:3,i))
                  CALL AddBulk(NextN(1:3,i),NextNN)
                ELSE
                  ListAdAtom(IndNN) % Ind_Event = Index_Event
                  ListAdAtom(IndNN) % ProbTrans = Prob
                !write(*,*)"Dep.full NN",Index_Event,Prob,SiOrC,CoorNN
                  CALL Updatetree(IndNN,Prob) !
                END IF
                IF(CoorOld.EQ.1.AND.i.GE.2)THEN ! if finds other full sites apart the 1st after the coordination fixing at lines 49-74 
                !write(*,*)'Increase Coor',Site,Coor
                  Coor=Coor+1
                END IF
              END IF
            END IF
          ELSE                                                            ! WALL
            Coor = Coor + 1                                               ! WALL
          END IF
        END IF
      END DO
      IndVoid=LattInd(Site(1),Site(2),Site(3)) ! it can change due to the Next Neigh. update
      IF(Coor.NE.CoorOld)CALL MVBITS(Coor,0,LenCoor,LattCoo(Site(1),
     >                               Site(2),Site(3)),PosCoor)
      CALL Get_Prob(1,IAtom,0,Coor,Site,NextN,Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C ! Here CovInd=0 by definition
      
      IF(Index_Event.EQ.0)THEN
                !write(*,*)"Dep. Bulk Site Generation",Index_Event,Prob
        CALL EraseMC(Site(1:3))
        CALL AddBulk(Site(1:3),NextN)
      ELSE
        ListAdAtom(IndVoid) % Ind_Event = Index_Event
        ListAdAtom(IndVoid) % ProbTrans = Prob
      !write(*,*)"Dep. IndVoid Index_Event,Prob",Index_Event,Prob
        CALL Updatetree(IndVoid,Prob)
      END IF
      END SUBROUTINE Deposition
