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
      SUBROUTINE SetProbabilityFromField()

      USE DefDerType
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER i,j,errcheck,maxElines,indsp,tableidx
      REAL(8) probvalue
      INTEGER elind
      INTEGER indices(6)
  
      INTEGER :: IPF51

      ! Read table filename from start.dat
      READ(IPF,'(A)')tempfilename
      write(*,*)'Reading input probability table: ', tempfilename
      
      ! Open file
      IPF51=51
      OPEN(IPF51,FILE=tempfilename,STATUS='OLD')  
      READ(IPF51,*) ntemp ! number of temperature values at which probabilities were evaluated
      READ(IPF51,*) nprob ! number of probabilities to read in every block 

      ! Get number of blocks in tempfilename
      DO tableidx=1,ntemp

        ! Store tabulated temperature value 
        READ(IPF51,*) reftempvalues(tableidx) ! temperature value at which the next nprob probabilities were evaluated
!        WRITE(*,*)tableidx, reftempvalues(tableidx)

        ! Store tabulated PtransD, PtransE, PtransDes and PtransAbs as arrays of max dimension ntempMax, but effectively ntemp

        !!!!!!!! EVAPORATION !!!!!!!!

        !refPtransE=0.d0 ! not allowed case (e.g. 1,0,0) or (1,3,1) the default is 0

        ! Set Max number of evaporation prob entries 
        IF((NCrystal.EQ.1).AND.(NCov.EQ.0))THEN        
          maxElines=3
        ELSE IF((NCrystal.EQ.2).AND.(NCov.EQ.0))THEN
          maxElines=18
        ELSE IF((NCrystal.EQ.3).AND.(NCov.EQ.0))THEN
          maxElines=57
        ELSE IF((NCrystal.EQ.1).AND.(NCov.EQ.1))THEN
          maxElines=9
        ELSE IF((NCrystal.EQ.2).AND.(NCov.EQ.1))THEN
          maxElines=50
        ELSE IF((NCrystal.EQ.3).AND.(NCov.EQ.1))THEN
          maxElines=150
        ELSE IF((NCrystal.EQ.1).AND.(NCov.EQ.2))THEN
          maxElines=19
        ELSE IF((NCrystal.EQ.2).AND.(NCov.EQ.2))THEN
          maxElines=100
        ELSE IF((NCrystal.EQ.3).AND.(NCov.EQ.2))THEN
          maxElines=288
        ELSE IF((NCrystal.EQ.1).AND.(NCov.EQ.3))THEN
          maxElines=34
        ELSE IF((NCrystal.EQ.2).AND.(NCov.EQ.3))THEN
          maxElines=172
        ELSE IF((NCrystal.EQ.3).AND.(NCov.EQ.3))THEN
          maxElines=480
        ELSE
          write(*,*)'ERROR: NCrystal and NCov must be <=',
     >              NCrystalMax, NCovMax
          STOP
        END IF
        
        IF(tableidx.EQ.1)THEN
          write(*,*)'maxElines',maxElines
        END IF

        ! Read evaporation probabilities
        DO i=1,maxElines
          READ(IPF51,*)indsp,elind,probvalue
          IF(indsp.GT.NCrystal)THEN
            write(*,*)'ERROR: 1st index in PtransE must be <= NCrystal'
            STOP
          END IF
          indices=0
          DO j=1,NCrystal
            ! below is a simple way to select the jth digit from a number of length NCrystal+NCov 
            indices(j) = MOD(elind/(10**(NCrystal+NCov-j)),10)
          END DO
          IF(NCov.GT.0)THEN
            DO j=NCrystal+1,NCrystal+NCov
              indices(NCrystalMax+j-NCrystal)=MOD(elind/(10**(NCrystal
     >         +NCov-j)),10)
            END DO
          END IF
          refPtransE(tableidx,indsp,indices(1),indices(2),indices(3),
     >     indices(4),indices(5),indices(6)) = probvalue
C  C        write(*,*)indsp,indices(:NCrystal),indices(NCrystalMax+1:
C  C     >     NCrystalMax+NCov),PtransE(indsp,indices(1),indices(2),
C  C     >     indices(3),indices(4),indices(5),indices(6))
        END DO

        !!!!!!!! DEPOSITION !!!!!!!!

        DO j=1,NCrystal*3
          READ(IPF51,*)indsp,elind,refPtransD(tableidx,indsp,elind)
          IF(indsp.GT.NCrystal)THEN
            write(*,*)'ERROR: 1st index in depo prob block must be 
     > <= NCrystal'
            STOP
          END IF
C  C         write(*,*)indsp,elind,PtransD(indsp,elind)
        END DO


        !!!!!!!! ABSORPTION !!!!!!!!

        DO j=1,NCov*3
          READ(IPF51,*)indsp,elind,refPtransAbs(tableidx,indsp-NCrystal,
     >                                         elind)
          IF((indsp.LE.NCrystal).OR.(indsp.GT.NCrystal+NCov))THEN
            write(*,*)'ERROR: 1st index in absorption prob block must be 
     > NCrystal < index <= NCrystal+NCov'
            STOP
          END IF
C  C         write(*,*)indsp,elind,PtransAbs(indsp-NCrystal,elind)
        END DO


        !!!!!!!! DESORPTION !!!!!!!!

        DO j=1,NCov*3
          READ(IPF51,*)indsp,elind,refPtransDes(tableidx,indsp-NCrystal,
     >                                         elind)
          IF((indsp.LE.NCrystal).OR.(indsp.GT.NCrystal+NCov))THEN
            write(*,*)'ERROR: 1st index in desorption prob block must be 
     > NCrystal < index <= NCrystal+NCov'
            STOP
          END IF
C  C         write(*,*)indsp,elind,PtransDes(indsp-NCrystal,elind)
        END DO

      END DO

      CLOSE(IPF51)
      write(*,*)'Done reading input probability table'


      ! The actual PtransD, PtransE, PtransDes and PtransAbs will be set within a InterpProb.f routine every time get_prob or get_prob_ini are called. 
      ! Call like this:
      ! write(*,*) refPtransE(tableidx,:,:,:,:,:,:,:)
      ! write(*,*) refPtransD(tableidx,:,:)
      ! write(*,*) refPtransAbs(tableidx,:,:)
      ! write(*,*) refPtransDes(tableidx,:,:)


      END SUBROUTINE SetProbabilityFromField
