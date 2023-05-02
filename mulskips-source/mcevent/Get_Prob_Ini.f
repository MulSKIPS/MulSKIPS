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
**************************************************************************8

      SUBROUTINE Get_Prob_Ini(Occ,NSiC,Coor,Site,NextN,  ! Only for pure 3C Configurations
     >                    Index_Event,Prob_Event) ! Occ=1/0 NSiC 1=001 Si, 2=010 C 3=011 B
      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Occ,NdiC,Coor,CoorNN,CovInd
      INTEGER Site(3),NextN(3,4)
      INTEGER i,Index_Event,SiOrC,NSiC,N_C,N_Si,N_B,N_cov
      REAL(8) Prob_Event
      REAL(8):: Ptot,Pcurr
      REAL(8) :: random
      EXTERNAL random
      REAL(8) :: kb=0.00008617330350d0 ! [eV/K]
      INTEGER temperature

      

      ! Only if Silicon trench, interpolate probabilities from external table based on local temperature value 
      IF((InitSt.EQ.'ST').AND.(.NOT.fixedT))THEN
!      IF(InitSt.EQ.'ST')THEN
        temperature = IBITS(LattCoo(Site(1),Site(2),Site(3)),PosT,LenT) ! K
        IF(temperature.EQ.0)THEN
           write(*,*)'Temperature is zero in get_prob_ini().' 
           write(*,*)'Are you sure that it should be like this?'
           write(*,*)'IMPORTANT: check that the temperature is reset'
           write(*,*)'every time LattCoo is put =0 in the other '
           write(*,*)'routines (e.g. evaporation.f)!!!'
           write(*,*)'Stopping now...'
           STOP
        END IF
        CALL InterpProbFromField(temperature)
      END IF




      IF(InitSt.EQ.'LA')THEN
        ! Ref. La Magna etal PRB 75, 235201 (2007)
        ! PtransE here is n_s * (Psi_L - Psi_S) (it has to be divided by kb*T)
        ! PtransD here is 2 * (Psi_L - Psi_S) / Tm (it has to be divided by kb)
        temperature = IBITS(LattCoo(Site(1),Site(2),Site(3)),PosT,LenT) ! K
        IF(temperature.EQ.0)THEN
           write(*,*)'Temperature is zero in get_prob_ini().' 
           write(*,*)'Are you sure that it should be like this?'
           write(*,*)'IMPORTANT: check that the temperature is reset'
           write(*,*)'every time LattCoo is put =0 in the other '
           write(*,*)'routines (e.g. evaporation.f)!!!'
           write(*,*)'Stopping now...'
           STOP
        END IF
        ! Melting
        PE = P0 * EXP(-PtransE/(kb*temperature)) 
     >          * merge(1.d0,0.d0,PtransE.NE.0.d0)        ! Set PE to zero for all array elements not specified in start.dat (e.g. (1,0,0) or (1,3,1))
        ! Solidification
        PD = 0.5d0 * P0 * EXP(-PtransD/(kb*Tm)) * 
     >    (1.0d0 + ERF((temperature - Tflex) / sigma))
     >            * merge(1.d0,0.d0,PtransD.NE.0.d0)        ! Set PD to zero for all array elements not specified in start.dat (so this determines whether 2nd or third crystal species will be deposited)
C         PD = 0.5d0 * P0 * EXP(-PtransD/kb)    
      ELSE
        PE = PtransE
        PD = PtransD
      END IF

      IF (Occ.EQ.0)THEN ! Depositions
         IF (Coor.EQ.1)THEN
            CoorNN=IBITS(LattCoo(NextN(1,1),NextN(2,1),NextN(3,1)),
     >                  PosCoor,LenCoor)
            IF(CoorNN.GE.2)THEN ! in the case Coor=1 the single first neighbor should have CoorNN=3
              Ptot=PD(1,Coor)+PD(2,Coor)+PD(3,Coor)
              Pcurr=random(idum)*Ptot
              IF(Pcurr.LE.PD(1,Coor))THEN
                Index_Event = 1
              ELSE IF(Pcurr.GT.(PD(1,Coor)+PD(2,Coor)))THEN
                Index_Event = 3        
              ELSE
                Index_Event = 2
              END IF
              Prob_Event = Ptot
            ELSE ! Event to be excluded but MC particle is present in the list
              Index_Event = 1
              Prob_Event = 1.e-6
            END IF
          ELSE IF ((Coor.EQ.2).OR.(Coor.EQ.3))THEN
            Ptot=PD(1,Coor)+PD(2,Coor)+PD(3,Coor)
            Pcurr=random(idum)*Ptot
            IF(Pcurr.LE.PD(1,Coor))THEN
               Index_Event = 1
            ELSE IF(Pcurr.GT.(PD(1,Coor)+PD(2,Coor)))THEN
               Index_Event = 3     
            ELSE
               Index_Event = 2
            END IF
            Prob_Event = Ptot
         ELSE IF(Coor.EQ.4)THEN! Vacancy Site?  Error
            Index_Event = 0
            Prob_Event = 0.d0
            write(*,*)Coor,'Get Prob Vacancy Generation'
         END IF
      ELSE ! Occ=1 Evaporation
         N_Si=0
         N_C =0
         N_B =0
         N_cov = 0  ! TODO: eventually there should be an array including N_Si, N_C, N_B, N_walls etc in case of more species. 
         !write(*,*)'Get Prob Site',Site,Occ,Coor
         DO i=1,4
           CovInd = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >            ,PosIndex,LenIndex)
           IF(CovInd.EQ.0)THEN ! only if the NN is not a coverage or a wall site
             SiOrC = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),
     >                  PosSiC,LenSiC)
             IF(SiOrC.EQ.1)N_Si=N_Si+1
             IF(SiOrC.EQ.2)N_C=N_C+1
             IF(SiOrC.EQ.3)N_B=N_B+1
C             write(*,*)'Get Prob SiOrC',NextN(1:3,i),SiOrC
           ELSE
             N_cov=N_cov+1
           END IF
         END DO
         IF(N_Si+N_C+N_B+N_cov.NE.Coor)THEN
            write(*,*)Coor,N_Si,N_C,N_B,N_cov,
     >           'Get Prob Error N_Si+N_C++N_B+N_cov.NE.Coor'
            write(*,*)Site,NextN            
            STOP
         END IF
         IF (N_Si+N_C+N_B+N_cov.GT.3)THEN
            Index_Event = 0
            Prob_Event = 0.d0
            !write(*,*)Coor,'Get Prob Bulk Site'
         ELSE IF(N_Si+N_C+N_B+N_cov.EQ.0)THEN
            Index_Event = 4
            Prob_Event = 100000.*Tree(1)  ! it makes its evaporation a fast event
            write(*,*)Coor,'Get Prob Isolated Atom'
         ELSE
            Index_Event = 4
            Prob_Event = PE(NSic,N_Si+N_cov,N_C,N_B,0,0,0)
         END IF
      END IF
      END SUBROUTINE Get_Prob_Ini
