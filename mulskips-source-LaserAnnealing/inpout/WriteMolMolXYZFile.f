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
      SUBROUTINE WriteMolMolXYZFile(FN,Time,Iter)

      USE  DefSystem
      USE  Definitions
      IMPLICIT NONE
      INTEGER IAt,IfSi,IfC,IcellSiC,JcellSiC,IGr,JGr,indAtom
      INTEGER x,y,z,x1,y1,z1,IfOcc,Coo,NCoorCorr,CovInd
      INTEGER NumAtTotal,NoEpiBulk,NoEpiSurf,NvoidMerged,NumAtWrong
      INTEGER NumIsWrong,ind_atotype,NumCovTotal
      INTEGER FN,NumSurfAt,Site(3),i,j,ic,sitec(3),NextN(3,4)
      REAL(8)     Time
      INTEGER, PARAMETER :: k10 = selected_int_kind(10)
      INTEGER(kind=k10) :: Iter
      REAL(8) rad,a1(3),a2(3),a3(3),cr(3)
      LOGICAL noepi
      REAL(8), dimension(3) :: alat
      INTEGER ll

!      DATA a1 / 1,  1, -1/  !  BCC
!      DATA a2 /-1,  1,  1/
!      DATA a3 / 1, -1,  1/

!      DATA a1 / 1,  1,  0/   ! FCC
!      DATA a2 / 0,  1,  1/
!      DATA a3 / 1,  0,  1/

       DATA a1 / 1.,  0.,   0./   ! HEX
       DATA a2 / 0.5,  0.866025037844386,  0./
!       DATA a3 / 0.,  0.,   1./
       DATA a3 / 0.,  0.,   1.414213562373/

C       sf = 4.63/12.0  ! valid only for SiC

      alat(1)=sf*real(LenX,kind=8)
      alat(2)=sf*real(LenY,kind=8)
      alat(3)=sf*real(LenZ,kind=8)

100   FORMAT('Si ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
101   FORMAT('B ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
102   FORMAT('C ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
103   FORMAT('Ge ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
104   FORMAT('O ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
105   FORMAT('H ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
106   FORMAT('Cl ',F10.5,' ',F10.5,' ',F10.5,'  ',F6.2)
121   FORMAT('Mg ',F10.5,' ',F10.5,' ',F10.5,' # SV (Si vacancy)')
122   FORMAT('He ',F10.5,' ',F10.5,' ',F10.5,' # CV (C vacancy)')
123   FORMAT('Ne ',F10.5,' ',F10.5,' ',F10.5,' # SAV (Si antisite vac)')
124   FORMAT('Be ',F10.5,' ',F10.5,' ',F10.5,' # CAV (C antisite vac)')
125   FORMAT('Ar ',F10.5,' ',F10.5,' ',F10.5,' # surrounded by 2C-2Si')

      ! Get total number of various particle types

      NumAtTotal=0
      NumCovTotal=0
      DO z=1,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            CovInd=IBITS(LattCoo(x,y,z),PosIndex,LenIndex)
            IF(CovInd.NE.113)THEN ! only if the NN is not a wall site
              IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
              Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
              IF(IfOcc.EQ.1.AND.Coo.LE.3)THEN 
                  NumAtTotal=NumAtTotal+1
              ELSE IF ((IfOcc.EQ.0).AND.(CovInd.NE.0))THEN
                  NumCovTotal=NumCovTotal+1
              END IF
            END IF 
         END DO
       END DO
      END DO

      NoEpiBulk=0
      DO z=0,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfOcc.EQ.1.AND.Coo.EQ.4)THEN
              noepi=.TRUE.
              IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.
     >           AND.MOD(z,3).EQ.0)THEN
                x1=x/3
                y1=y/3
                z1=z/3
                IF((MOD(x1,2).EQ.MOD(y1,2)).AND.
     >             (MOD(x1,2).EQ.MOD(z1,2)))THEN
                  IF(MOD(x1+y1+z1,4).EQ.0.OR.
     >               MOD(x1+y1+z1,4).EQ.3)THEN ! Si sites
                     noepi=.FALSE.
                  END IF
                END IF
              END IF
              IF(noepi)NoEpiBulk=NoEpiBulk+1
            END IF
         END DO
        END DO
      END DO



      ! +++ Write defective (non-epitaxial) +++

      write(*,*)NoEpiBulk
      WRITE(FN+3,'(1x,i8,1x,a,1x,a,1x,e15.8,1x,a,1x,i12)') NoEpiBulk,
     > 'angstroem','KMC-time:',Time,'Iter:',Iter
      WRITE(FN+3,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)

      DO z=0,LenZ-1
        DO y=0,LenY-1
         DO x=0,LenX-1
            IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            IF(IfOcc.EQ.1.AND.Coo.EQ.4)THEN
              noepi=.TRUE.
              IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.
     >           AND.MOD(z,3).EQ.0)THEN
                x1=x/3
                y1=y/3
                z1=z/3
                IF((MOD(x1,2).EQ.MOD(y1,2)).AND.
     >             (MOD(x1,2).EQ.MOD(z1,2)))THEN
                  IF(MOD(x1+y1+z1,4).EQ.0.OR.
     >               MOD(x1+y1+z1,4).EQ.3)THEN ! Si sites
                     noepi=.FALSE.
                  END IF
                END IF
              END IF
              IF(noepi)THEN
                IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                IF(IfSi.EQ.1)WRITE(FN+3,103)sf*x,sf*y,sf*z
                IF(IfC.EQ.1)WRITE(FN+3,101)sf*x,sf*y,sf*z
              END IF
            END IF
         END DO
        END DO
      END DO


      ! +++ Write undercoordinated + vacancies +++

      NvoidMerged=NumVoids
!      WRITE(FN,'(1x,i8,1x,a)') NumAtTotal+NvoidMerged,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
      WRITE(FN,'(1x,i8,1x,a,1x,a,1x,e15.8,1x,a,1x,i12)') 
     > NumAtTotal+NvoidMerged+NumCovTotal,'angstroem','KMC-time:',
     > Time,'Iter:',Iter
      WRITE(FN,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)

      DO z=0,LenZ-1
        DO y=0,LenY-1
          DO x=0,LenX-1
            IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
            IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
            Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
            CovInd = IBITS(LattCoo(x,y,z),PosIndex,LenIndex)
            IF(CovInd.EQ.1)WRITE(FN,105)sf*x,sf*y,sf*z !, rad
            IF(CovInd.EQ.17)WRITE(FN,106)sf*x,sf*y,sf*z !, rad
            IF(IfSi.EQ.1.AND.Coo.LE.3)WRITE(FN,100)sf*x,sf*y,sf*z !, rad
            IF(IfC.EQ.1.AND.Coo.LE.3)THEN
              IF(NCrystal.EQ.1)THEN
                WRITE(*,*)'ERR: IfC=1 should not happen, NCrystal=1!!!'
                STOP
              END IF
              IF(ANY(ListCrystal.EQ.6))THEN
                WRITE(FN,102)sf*x,sf*y,sf*z !, rad
              ELSE IF(ANY(ListCrystal.EQ.32))THEN
                WRITE(FN,103)sf*x,sf*y,sf*z !, rad
              ELSE
                WRITE(*,*)'ERROR: 2nd crystal specie can only be C/Ge'
                STOP
              END IF
            END IF
          END DO
        END DO
      END DO
      DO i=1,NvoidMerged
        x = ListVoid (i) % AtomXYZ(1)
        y = ListVoid (i) % AtomXYZ(2)
        z = ListVoid (i) % AtomXYZ(3)
        WRITE(FN,104)sf*x,sf*y,sf*z !, rad
      END DO




      ! +++ Write Vacancies +++

      WRITE(FN+4,'(1x,i8,1x,a)') NumVoids,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
      WRITE(FN+4,'(1x,a,3(1x,e15.8))')'periodic',alat(1),alat(2),alat(3)

      DO i=1,NumVoids
        x = ListVoid (i) % AtomXYZ(1)
        y = ListVoid (i) % AtomXYZ(2)
        z = ListVoid (i) % AtomXYZ(3)
        ! only for SiC PVD. In CVD or doped cases we do not distinguish point defects for now...
        IF(NCov.EQ.0.AND.NCrystal.EQ.2.AND.ListCrystal(2).EQ.6)THEN 
          NextN = ListVoid(i) % NextNXYZ
          IAt=0
          DO j=1,4
            IAt=IAt+IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),
     >        PosSiC,LenSiC)
          END DO
          IF(IAt.EQ.4)THEN
            WRITE(FN+4,121)sf*x,sf*y,sf*z
          ELSE IF(IAt.EQ.8)THEN
            WRITE(FN+4,122)sf*x,sf*y,sf*z
          ELSE IF(IAt.EQ.7)THEN
            WRITE(FN+4,123)sf*x,sf*y,sf*z
          ELSE IF(IAt.EQ.5)THEN
            WRITE(FN+4,124)sf*x,sf*y,sf*z
          ELSE IF(IAt.EQ.6)THEN
            WRITE(FN+4,125)sf*x,sf*y,sf*z
          ELSE
            WRITE(*,*)'ERROR: Weird Defect type was formed...ABORTING'
            STOP
          END IF
        ELSE
          WRITE(FN+4,104)sf*x,sf*y,sf*z !, rad
        END IF
      END DO



      ! +++ Write wrong occupations (SiC only) +++

      ! Analyze and write "wrong" sites only in SiC epi growth case 
      IF((NCrystal.EQ.2).AND.(ANY(ListCrystal.EQ.6)))THEN

        NumIsWrong=0
        DO z=1,LenZ-1
          DO y=0,LenY-1
           DO x=0,LenX-1
              CovInd=IBITS(LattCoo(x,y,z),PosIndex,LenIndex)
              IF(CovInd.NE.113)THEN ! only if the NN is not a wall site
                IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
                Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
                IF (IfOcc.EQ.1.AND.Coo.EQ.4)THEN
                    indAtom=LattInd(x,y,z)
                    IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                    IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                    NextN=ListAtom(indAtom) % NextNXYZ
                    NCoorCorr=0
                    IF(IfSi.EQ.1)THEN
                        DO i=1,4
                            NCoorCorr=NCoorCorr+
     >    IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosC,LenC)
                        END DO
                    ELSE IF(IfC.EQ.1)THEN
                        DO i=1,4
                            NCoorCorr=NCoorCorr+
     >    IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosSi,LenSi)
                        END DO
                    END IF
                    IF(NCoorCorr.EQ.0)THEN ! only isolated 
                        NumIsWrong=NumIsWrong+1
                        IF(IfSi.EQ.1)THEN
                            LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosSi)
                            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                            ListAtBuffer(NumIsWrong) % AtomXYZ(1)=x
                            ListAtBuffer(NumIsWrong) % AtomXYZ(2)=y
                            ListAtBuffer(NumIsWrong) % AtomXYZ(3)=z
                            ListAtBuffer(NumIsWrong) % Ind_Atype = 1
                        ELSE 
                            LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosC)
                            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                            ListAtBuffer(NumIsWrong) % AtomXYZ(1)=x
                            ListAtBuffer(NumIsWrong) % AtomXYZ(2)=y
                            ListAtBuffer(NumIsWrong) % AtomXYZ(3)=z
                            ListAtBuffer(NumIsWrong) % Ind_Atype = 2
                        END IF
                    END IF
                END IF
              END IF 
           END DO
         END DO
        END DO

        NumAtWrong=0
        DO z=1,LenZ-1
          DO y=0,LenY-1
           DO x=0,LenX-1
              IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
              Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
              IF (IfOcc.EQ.1.AND.Coo.EQ.4)THEN
                  indAtom=LattInd(x,y,z)
                  IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                  IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                  NextN=ListAtom(indAtom) % NextNXYZ
                  NCoorCorr=0
                  IF(IfSi.EQ.1)THEN
                      DO i=1,4
                          NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosC,LenC)
                      END DO
                  ELSE IF(IfC.EQ.1)THEN
                      DO i=1,4
                          NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosSi,LenSi)
                      END DO
                  END IF
                  IF(NCoorCorr.EQ.3)NumAtWrong=NumAtWrong+1
              END IF
           END DO
          END DO
        END DO

        write(*,*)'NumAtWrong: ', NumAtWrong 

        WRITE(FN+5,'(1x,i8,1x,a)') NumAtWrong,'angstroem'  !+Lenx*Leny*Lenz+(Lenx*Leny*Lenz*8)/9
        WRITE(FN+5,'(1x,a,3(1x,e15.8))')'periodic',
     >                  alat(1),alat(2),alat(3)
        DO z=1,LenZ-1
          DO y=0,LenY-1
           DO x=0,LenX-1
              IfOcc = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
              Coo = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
              IF(IfOcc.EQ.1.AND.Coo.LE.3)THEN
                  NumAtTotal=NumAtTotal+1
              ELSE IF (IfOcc.EQ.1.AND.Coo.EQ.4)THEN
                  indAtom=LattInd(x,y,z)
                  IfSi=IBITS(LattCoo(x,y,z),PosSi,LenSi)
                  IfC=IBITS(LattCoo(x,y,z),PosC,LenC)
                  NextN=ListAtom(indAtom) % NextNXYZ
                  NCoorCorr=0
                  IF(IfSi.EQ.1)THEN
                      DO i=1,4
                          NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosC,LenC)
                      END DO
                  ELSE IF(IfC.EQ.1)THEN
                      DO i=1,4
                          NCoorCorr=NCoorCorr+
     > IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)),PosSi,LenSi)
                      END DO
                  END IF
                  IF(NCoorCorr.EQ.3)THEN 
                      IF(IfSi.EQ.1)WRITE(FN+5,100)sf*x,sf*y,sf*z
                      IF(IfC.EQ.1)WRITE(FN+5,102)sf*x,sf*y,sf*z
                  END IF 
              END IF
           END DO
          END DO
        END DO

        DO i=1,NumIsWrong
          x=ListAtBuffer(i) % AtomXYZ(1)
          y=ListAtBuffer(i) % AtomXYZ(2)
          z=ListAtBuffer(i) % AtomXYZ(3)
          ind_atotype=ListAtBuffer(NumIsWrong) % Ind_Atype 
          IF(ind_atotype.EQ.1)THEN
              LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
          ELSE
              LattCoo(x,y,z)=IBCLR(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
          END IF
        END DO 
        
        CLOSE(FN+5)
      
      END IF



      CLOSE(FN)
      CLOSE(FN+3)
      CLOSE(FN+4)

      END SUBROUTINE WriteMolMolXYZFile
******|****************************************************************
