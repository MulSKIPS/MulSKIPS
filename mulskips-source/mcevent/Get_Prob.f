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
******************************************************************************
      FUNCTION IndexEvent_Abs(Coordination,Ptmp,PDepo)
        USE DefDerType
        USE DefSystem
        USE Definitions
        IMPLICIT NONE
        INTEGER, INTENT(In) :: Coordination
        REAL(8), INTENT(In) :: Ptmp,PDepo
        REAL(8) :: random
        EXTERNAL random
        REAL(8) :: pmin,pmax
        INTEGER :: i,IndexEvent_Abs        
        
        IF(NCov.EQ.1)THEN
          IndexEvent_Abs=5
        ELSE IF(Ncov.GE.2)THEN
          DO i=1,NCov  ! loop over various segments
            IF(i.EQ.1)THEN
              pmin=PDepo
              pmax=PDepo+PtransAbs(1,Coordination)
            ELSE
              pmin=PDepo+SUM(PtransAbs(:i-1,Coordination))
              pmax=PDepo+SUM(PtransAbs(:i,Coordination))
            END IF
            IF((Ptmp.GT.pmin).AND.(Ptmp.LE.pmax))THEN
              IndexEvent_Abs = 4+i
C               write(*,*)'Absorbing',IndexEvent_Abs
            END IF
          END DO
        END IF
        RETURN
      END FUNCTION

      SUBROUTINE Get_Prob(Occ,NSiC,CovInd,Coor,Site,NextN,
     >                    Index_Event,Prob_Event) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Occ,OccN,NdiC,Coor,CoorNN,IndNN
      INTEGER Site(3),NextN(3,4),NextNN(3,4)
      INTEGER i,j,Index_Event,SiOrC,NSiC,N_C,N_Si,N_wall,CovIndN         ! WALL added CovInd
      INTEGER CovInd
      INTEGER IndexEvent_Abs
      REAL(8) Prob_Event
      LOGICAL Recipr
      REAL(8):: Ptot,Pcurr
      REAL(8) :: random
      EXTERNAL random
      REAL(8) :: kb=0.00008617330350d0 ! [eV/K]
      INTEGER temperature
      INTEGER, DIMENSION(Ncov) :: covNN 
      INTEGER iCov,covNNsum
      REAL(8) :: PDep,PAbs

      IF((CovInd.GT.0).AND.(CovInd.LT.113))THEN ! desorption 
        IF(Occ.NE.0)THEN
          write(*,*)'ERROR Covind!=0 with Occ!=0 in get_prob'
          STOP
        END IF
        IF(Coor.LE.3)THEN
          Index_Event = 4
          iCov=0
          DO j=1,Ncov
!            write(*,*)j, ListCov(j)
            IF (ListCov(j).EQ.CovInd)THEN
              iCov=j
              EXIT
            END IF
          END DO
          IF(iCov.EQ.0)THEN
            write(*,*)'ERROR: could not find CovInd=', CovInd, 
     >           ' inside ListCov'
!            write(*,*)Occ,NSiC,CovInd,Coor,Site,NextN
            STOP
          END IF
          Prob_Event = PTransDes(iCov,Coor)
        ELSE IF(Coor.GE.4)THEN ! Behave like Vacancy Site
          Index_Event = 0
          Prob_Event = 0.d0
          IF(Coor.GT.4)THEN
            write(*,*)'anomalous coordination in desorption'
          END IF
        ELSE
           write(*,*)'anomalous coordination in desorption',Coor
        END IF
      ELSE ! if wall or crystal

        ! LA case
        IF(InitSt.EQ.'LA')THEN
          ! Ref. La Magna etal PRB 75, 235201 (2007)
          ! PtransE here is n_s * (Psi_L - Psi_S) (it has to be divided by kb*T)
          ! PtransD here is 2 * (Psi_L - Psi_S) / Tm (it has to be divided by kb)
          temperature = IBITS(LattCoo(Site(1),Site(2),Site(3)),
     >        PosT,LenT) ! K
          ! Melting
!          write(*,*)"temperature ", temperature
!          write(*,*)"PtransD ", PtransD
          PE = P0 * EXP(-PtransE/(kb*temperature)) 
     >            * merge(1.d0,0.d0,PtransE.NE.0.d0)        ! Set PE to zero for all array elements not specified in start.dat (e.g. (1,0,0) or (1,3,1))
          ! Solidification
          PD = 0.5d0 * P0 * EXP(-PtransD/(kb*Tm)) * 
     >      (1.0d0 + ERF((temperature - Tflex) / sigma)) 
     >            * merge(1.d0,0.d0,PtransD.NE.0.d0)        ! Set PD to zero for all array elements not specified in start.dat (so this determines whether 2nd or third crystal species will be deposited)
!          write(*,*)"PD ", PD
!         PD = 0.5d0 * P0 * EXP(-PtransD/kb)    
        ELSE
          PE = PtransE
          PD = PtransD
        END IF      

        IF (Occ.EQ.0)THEN ! Depositions or absorptions
          IF (Coor.EQ.1)THEN
             CoorNN=IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosCoor,LenCoor)
             IF(CoorNN.GE.2)THEN ! in the case Coor=1 the single first neighbor should have CoorNN>=2
               PDep=PD(1,Coor)+PD(2,Coor)
               PAbs=SUM(PtransAbs(:,Coor)) ! for LA this will always be zero. Depo will always be selected
               Ptot=PDep+PAbs
!               write(*,*)'Probs',Ptot,PD(1,Coor),PD(2,Coor),PAbs
               Pcurr=random(idum)*Ptot
               IF(Pcurr.LE.PDep)THEN  ! deposition
                 IF(Pcurr.LE.PD(1,Coor))THEN
!                   write(*,*)'depositing Si'
                   Index_Event = 1
                 ELSE
!                   write(*,*)'depositing C'
                   Index_Event = 2
                 END IF
               ELSE  ! absorption
                 Index_Event = IndexEvent_Abs(Coor,Pcurr,PDep)
               END IF
               Prob_Event = Ptot
             ELSE ! Event to be excluded but MC particle is present in the list
               Index_Event = 1
               Prob_Event = 1.e-8                                          ! WALL TODO: might be decreased to e.g. 1e-20 to avoid Tree errors due to very small Ptrans... 
             END IF
          ELSE IF ((Coor.EQ.2).OR.(Coor.EQ.3))THEN
             PDep=PD(1,Coor)+PD(2,Coor)
             PAbs=SUM(PtransAbs(:,Coor))  ! for LA this will always be zero. Depo will always be selected
             Ptot=PDep+PAbs
             Pcurr=random(idum)*Ptot
             IF(Pcurr.LE.PDep)THEN  ! deposition
               IF(Pcurr.LE.PD(1,Coor))THEN
!                 write(*,*)'depositing Si'
                 Index_Event = 1
               ELSE
!                 write(*,*)'depositing C'
                 Index_Event = 2
               END IF
             ELSE  ! absorption
               Index_Event = IndexEvent_Abs(Coor,Pcurr,PDep)
             END IF
             Prob_Event = Ptot
          ELSE IF(Coor.GE.4)THEN! Vacancy Site?  Error
             Index_Event = 0
             Prob_Event = 0.d0
             !write(*,*)Coor,'Get Prob Vacancy Generation'
             IF(Coor.GT.4)THEN
                write(*,*)'anomalous coordination'
             END IF
          ELSE
             write(*,*)'anomalous coordination'
             write(*,*)Coor
          END IF
        ELSE ! Occ=1 Evaporation
          N_Si=0
          N_C =0
          N_wall = 0                                                        ! WALL ! TODO: eventually there should be an array including N_Si, N_C, N_H, N_walls etc in case of more species. 
          covNN = 0  ! array of Len = Ncov
          !write(*,*)'Get Prob Site',Site,LattInd(Site(1),Site(2),Site(3))
          DO i=1,4
            OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                PosOcc,LenOcc)
            CovIndN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))     ! WALL
     >            ,PosIndex,LenIndex)                                     ! WALL
            CoorNN=IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                    PosCoor,LenCoor)
            IF(OccN.NE.0)THEN
              IF(CovIndN.NE.113)THEN                                          ! WALL Go on only if NextN is a crystal site, otherwise increase N_wall at line 128
                IF(CoorNN.GE.1.AND.CoorNN.LE.3)THEN
                  IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
                  NextNN = ListAdAtom(IndNN) % NextNXYZ
                  Recipr=.FALSE.
                  DO j=1,4 ! check is the site belongs to the NN of its NN
                    IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
                  END DO
!             write(*,*)'Get Prob Recipr',Recipr
                  IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                    SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i)
     >                   ,NextN(3,i)),PosSiC,LenSiC)
                    IF(SiOrC.EQ.1)N_Si=N_Si+1
                    IF(SiOrC.EQ.2)N_C=N_C+1   ! When Ncrystal=1 this should always be zero, because Ind_Event will never be =2 (see above)
!              write(*,*)'Get Prob SiOrC',NextN(1:3,i),SiOrC
                  END IF
                ELSE ! ripete la stessa cosa anche per CoorNN.EQ.4...per ora lasciamo così
                  IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
                  NextNN = ListAtom(IndNN) % NextNXYZ
                  Recipr=.FALSE.
                  DO j=1,4 ! check is the site belongs to the NN of its NN
                    IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
                  END DO
!             write(*,*)'Get Prob Recipr',Recipr
                  IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                    SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                   NextN(3,i)),PosSiC,LenSiC)
                    IF(SiOrC.EQ.1)N_Si=N_Si+1
                    IF(SiOrC.EQ.2)N_C=N_C+1
                  END IF
                END IF
              ELSE                                                         ! WALL 
                N_wall=N_wall+1                                              ! WALL 
              END IF
            ELSE ! I should also look for coverage neighbours, which have Occ=0
              IF(CovIndN.NE.0)THEN  ! only if coverage (this will never occur if Ncov=0 in defsystem.f)
                IndNN = LattInd(NextN(1,i),NextN(2,i),NextN(3,i))
                NextNN = ListAdAtom(IndNN) % NextNXYZ
                Recipr=.FALSE.
                DO j=1,4 ! check is the site belongs to the NN of its NN
                  IF(ALL(NextNN(1:3,j).EQ.Site(1:3)))Recipr=.TRUE.
                END DO
                IF (Recipr) THEN ! action only if the site belongs to the NN of its NN
                  iCov=0
                  DO j=1,Ncov
                    IF (ListCov(j).EQ.CovIndN)THEN
                      iCov=j
                      EXIT
                    END IF
                  END DO
                  IF(iCov.EQ.0)THEN
                    write(*,*)'ERROR: could not find CovIndN=', CovIndN, 
     >                ' inside ListCov'
                    STOP
                  END IF
                  covNN(iCov)=covNN(iCov)+1
                END IF               
              END IF
            END IF
          END DO
          covNNsum = SUM(covNN) ! in case there is no coverage atoms, covNNsum=0 and usual behaviour is retrieved. In case there is only H, covNNsum = covNN(1)
          IF(N_Si+N_C+N_wall.NE.Coor)THEN                              ! note: Coor is not changed by coverage!
             write(*,*)Coor,N_Si,N_C,N_wall,covNNsum               
     >          ,'Get Prob Error N_Si+N_C+N_wall.NE.Coor'
             write(*,*)Site,Occ,NSiC,CovInd,Coor,covNN
             write(*,*)covNN
             STOP
          END IF
          IF (N_Si+N_C+N_wall.GT.3)THEN                                ! note: Coor is not changed by coverage!
             Index_Event = 0
             Prob_Event = 0.d0
!             write(*,*)Coor,'Get Prob Bulk Site'
          ELSE IF(N_Si+N_C+N_wall.EQ.0)THEN                            ! note: Coor is not changed by coverage! An isolated Si should evaporate even if it has a coverage NN
             Index_Event = 3
             Prob_Event = 100000.*Tree(1)  ! it makes its evaporation a fast event
!             write(*,*)Coor,'Get Prob Isolated Atom'
          ELSE
! We'll need to fix the 4th arg in PE below if we want to include a third crystal species! We'll count NC1, NC2, NC3 instead of N_Si and N_C
             Index_Event = 3                                           ! for now it works only with H or Cl. What if there are other species? Will every species be Hard coded? Should we define an effective probability scaling factor as an input array in start.dat?
             IF(NCov.EQ.1)THEN
               Prob_Event = PE(NSiC,N_Si+N_wall,N_C,0,covNN(1),0,0)    ! WALL: reduce probability if a NN is a wall. If more than 1 NN is wall, the effect sums up
             ELSE IF(NCov.EQ.2)THEN
               Prob_Event = PE(NSiC,N_Si+N_wall,N_C,0,covNN(1),
     >                                                     covNN(2),0)    ! WALL: reduce probability if a NN is a wall. If more than 1 NN is wall, the effect sums up
             ELSE IF(NCov.EQ.3)THEN
               Prob_Event = PE(NSiC,N_Si+N_wall,N_C,0,covNN(1),
     >                                                     covNN(2),
     >                                                     covNN(3))    ! WALL: reduce probability if a NN is a wall. If more than 1 NN is wall, the effect sums up
             ELSE
               Prob_Event = PE(NSiC,N_Si+N_wall,N_C,0,0,0,0)    ! WALL: reduce probability if a NN is a wall. If more than 1 NN is wall, the effect sums up
!               write(*,*)NSiC,N_Si,N_C,NCov,N_wall,'Get Prob PE()'
             END IF             
          END IF
        END IF
      END IF

      IF(Index_Event.NE.0.AND.Prob_Event.LT.1e-12)THEN  ! DO NOT INCREASE THIS THRESHOLD VALUE. Otherwise the tree might mess up 
        write(*,*)Index_Event,Prob_Event,Coor
      END IF

      END SUBROUTINE Get_Prob
