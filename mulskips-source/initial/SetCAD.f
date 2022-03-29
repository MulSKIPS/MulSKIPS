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
**************************************************************************************
**
**
**************************************************************************************
      SUBROUTINE SetCAD()
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER i,j,x,y,z,x1,y1,z1,x0,y0,z0
      INTEGER Id,CovInd,OccN,Occ,IAt,CovIndN,Coor,CoorNN
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob

      INTEGER :: temperature

      INTEGER :: IPF50, IPF60
      REAL(8) :: LenIn(LenX*LenY)



      ! CASE 1
      IF(RunType.EQ.'C')THEN

        write(*,*)'Reading checkpoint file', cadfilename
        IPF60=60
        OPEN(IPF60,FILE=cadfilename,STATUS='OLD')     
        DO z=0,LenZ-1
          READ(IPF60,*)LenIn
          DO y=0,LenY-1
            DO x=0,LenX-1
              IF(InitSt.EQ.'LA')THEN
                ! Needed, otherwise T will be overwritten!
                temperature = IBITS(LattCoo(x,y,z),PosT,LenT) ! K
                LattCoo(x,y,z) = LenIn(1+x+LenX*y)
                CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
                ! Check max T
                IF(temperature.GT.tmax)THEN
                  tmax = temperature
                END IF
              ELSE  
                LattCoo(x,y,z) = LenIn(1+x+LenX*y)
              END IF
            END DO
          END DO
        END DO
        CLOSE(IPF60)
        write(*,*)'Done reading checkpoint file'

      ELSE

        !  Filled wall Sites at z=0 (solves boundary vacancies bug and overwrites geo setting)
        DO z=0,3
         DO y=0,LenY-1
          DO x=0,LenX-1
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            CovInd=113 ! coverage (for now it indicates walls)
            CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
          END DO
         END DO
        END DO            

        !  Filled wall Sites at z=LenZ (solves boundary vacancies bug and overwrites geo setting)
        DO z=LenZ-4,LenZ-1
         DO y=0,LenY-1
          DO x=0,LenX-1
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            CovInd=113 ! coverage (for now it indicates walls)
            CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
          END DO
         END DO
        END DO                   

        write(*,*)'Reading input geometry file: ', cadfilename
        IPF50=50
        OPEN(IPF50,FILE=cadfilename,STATUS='OLD')     
        ! read and ignore lines in cadfilename which relate to z of wall sites
        DO z=0,3
          READ(IPF50,*)LenIn 
        END DO  
        ! Now start reading properly cadfilename           
        DO z=4,LenZ-5
          READ(IPF50,*)LenIn
          DO y=0,LenY-1
            DO x=0,LenX-1
              Id = LenIn(1+x+LenX*y)
              ! Set walls
              IF(Id.EQ.10)THEN ! WALL, AIR, SiO2
                LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
                CovInd=113 
                CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
              ! Set species
              ELSE IF(Id.NE.0)THEN! Si  
                IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0
     >                .AND.MOD(z,3).EQ.0)THEN
                  x1=x/3
                  y1=y/3
                  z1=z/3
                  IF((MOD(x1,2).EQ.MOD(y1,2)).AND.
     >                 (MOD(x1,2).EQ.MOD(z1,2)))THEN
                    IF((MOD(x1+y1+z1,4).EQ.0).OR.
     >                   (MOD(x1+y1+z1,4).EQ.3))THEN ! Si sites                    
                      IF(InitSt.EQ.'LA')THEN
                        ! we need to nucleate at least 2 nm of liquid
                        IF(Homogeneous.EQ.'T')THEN
                        ! here we do it homogeneously as a thin layer
                          IF(z.LE.(LenZ-INT((LenVac+LenNuc)/sf)))THEN 
                            ! LenVac + LenNucleus in [ang]
                            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
                            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)                      
                          END IF
                        ELSE IF(Homogeneous.EQ.'F')THEN
                        ! here instead we do it dishomogeneously, starting from a semi-spherical nucleus of radius LenNuc (Ang)
                          x0 = INT(LenX/2)
                          y0 = INT(LenY/2)
                          z0 = LenZ-INT(LenVac/sf)
                          IF(((x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0))
     >                      .GT.(INT((LenNuc/sf)*(LenNuc/sf))))THEN 
                            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
                            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)                      
                          END IF
                        ELSE
                          WRITE(*,*)'ERROR: Homogeneous needs to be eith
     >er T or F'
                        END IF
                      ELSE IF(InitSt.EQ.'IN')THEN
                        LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
                        LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)                      
                      END IF
                    END IF
                  END IF  
                END IF  
              END IF 
            END DO
          END DO
        END DO
        CLOSE(IPF50)
        write(*,*)'Done reading input geometry file'
        
        !  Set the coordination
        DO z=0,LenZ-1
         DO y=0,LenY-1
          DO x=0,LenX-1
           CovInd=IBITS(LattCoo(x,y,z),PosIndex,LenIndex)
           IF(CovInd.NE.113)THEN ! only if it is not a wall site 
            IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
             x1=x/3
             y1=y/3
             z1=z/3
             IF((MOD(x1,2).EQ.MOD(y1,2)).AND.
     >         (MOD(x1,2).EQ.MOD(z1,2)))THEN
              IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
               Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
               Coor=0
               Site(1)=x
               Site(2)=y
               Site(3)=z
               DO i=1,4
                  NextN(1:3,i) = Site + JumpDia(1:3,i)
                  IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
                  IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
                  IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
                  IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
                  OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
                  IF(Occ.GT.0)THEN
                    IF(OccN.GT.0)Coor=Coor+1 ! this should include walls!             
                  ELSE
                    CovIndN = IBITS(LattCoo(NextN(1,i),NextN(2,i)
     >                     ,NextN(3,i)),PosIndex,LenIndex)
                    IF((OccN.GT.0).AND.(CovIndN.NE.113))Coor=Coor+1 ! if a NN of an empty site is wall, then its coor should not increase 
                  END IF
               END DO
               IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
                 IF(InitSt.EQ.'LA')THEN
                   temperature = IBITS(LattCoo(x,y,z),PosT,LenT) ! K
                   LattCoo(x,y,z)=0
                   CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
                 ELSE  
                   LattCoo(x,y,z)=0
                 END IF               
                 DO i=1,4
                   CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                   CoorNN=CoorNN-1
                   IF(CoorNN.LT.0)CoorNN=0
                   CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
                 END DO
               ELSE
                 CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
               END IF
              ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! Si sites
               Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
               Coor=0
               Site(1)=x
               Site(2)=y
               Site(3)=z
               DO i=1,4
                  NextN(1:3,i) = Site + JumpDia(1:3,i+4)
                  IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
                  IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
                  IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
                  IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
                  OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
                  IF(Occ.GT.0)THEN
                    IF(OccN.GT.0)Coor=Coor+1 ! this should include walls!             
                  ELSE
                    CovIndN = IBITS(LattCoo(NextN(1,i),NextN(2,i)
     >                     ,NextN(3,i)),PosIndex,LenIndex)
                    IF((OccN.GT.0).AND.(CovIndN.NE.113))Coor=Coor+1 ! if a NN of an empty site is wall, then its coor should not increase 
                  END IF
               END DO
               IF(Occ.NE.0.AND.Coor.EQ.0)THEN ! Remove isolated atom
                 IF(InitSt.EQ.'LA')THEN
                   temperature = IBITS(LattCoo(x,y,z),PosT,LenT) ! K
                   LattCoo(x,y,z)=0
                   CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
                 ELSE  
                   LattCoo(x,y,z)=0
                 END IF    
                 DO i=1,4
                   CoorNN = IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                          NextN(3,i)),PosCoor,LenCoor)
                   CoorNN=CoorNN-1
                   IF(CoorNN.LT.0)CoorNN=0
                   CALL MVBITS(CoorNN,0,LenCoor,LattCoo(NextN(1,i),
     >                       NextN(2,i),NextN(3,i)),PosCoor)
                 END DO
               ELSE
                 CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
               END IF
              END IF
             END IF
            END IF
           END IF
          END DO
         END DO
        END DO

      END IF


!     Set the probability
      DO z=0,LenZ-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         CovInd=IBITS(LattCoo(x,y,z),PosIndex,LenIndex)
         IF(CovInd.NE.113)THEN ! only if it is not a wall site 
          IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
           x1=x/3
           y1=y/3
           z1=z/3
           IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
            IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
                NextN(1:3,i) = Site + JumpDia(1:3,i)
                IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
                IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
                IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
                IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied or walls
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    CovIndN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosIndex,LenIndex)
                    IF((OccN.EQ.0).OR.(CovIndN.EQ.113))THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 IF(Coor.EQ.4)THEN
                   write(*,*)'Initializing a vacancy'
                   CALL AddVoid(Site,NextN)
                 ELSE
                   CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                   CALL AddMC(Site,NextN,Index_Event,Prob)
                 END IF
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
            ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! Si sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
                NextN(1:3,i) = Site + JumpDia(1:3,i+4)
                IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
                IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
                IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
                IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied or walls
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    CovIndN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosIndex,LenIndex)
                    IF((OccN.EQ.0).OR.(CovIndN.EQ.113))THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 IF(Coor.EQ.4)THEN
                   write(*,*)'Initializing a vacancy'
                   CALL AddVoid(Site,NextN)
                 ELSE
                   CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                   CALL AddMC(Site,NextN,Index_Event,Prob)
                 END IF
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
            END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      END SUBROUTINE SetCAD

**************************************************************************************
