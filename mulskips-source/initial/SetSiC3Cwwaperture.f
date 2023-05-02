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
      SUBROUTINE SetSiC3Cwwaperture(Len,Lenww,Laperture)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len,Lenww,Laperture
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,CovInd
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),Buf(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled wall Sites at z=0 (helps with migrating vacancies issue)
      DO z=0,5
       DO y=0,LenY-1
        DO x=0,LenX-1
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        END DO
       END DO
      END DO      
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      DO z=6,Len-4
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4 ! coordination
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
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
            CALL AddBulk(Site,NextN)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=4
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
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
            CALL AddBulk(Site,NextN)
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
!  Filled  Sites in 100 SiC lattice MC particles..NumAdatoms
      DO z=Len-3,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si Atoms
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >            Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
            Prob = PtransE(1,2,0,0,0,0,0)! Si atom on SiC Surface
            Index_Event = 4
            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C Atoms
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Coo=0
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
     >           Coo=Coo+1
            END DO
!            write(*,*)Site,Coo,LattCoo(Site(1),Site(2),Site(3))
            Prob = PtransE(2,2,0,0,0,0,0)! C atom on SiC Surface
            Index_Event = 4
            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAdAtom",NumAdAtom,Len
      END DO
!  Empty  Sites in 100 SiC lattice (to be also consider MC particles in NumAdatoms)
      DO z=Len,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
            CALL Get_Prob(0,0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
            Buf=0
!            write(*,*)'Site',Site
            DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
!               write(*,*)NextN(1:3,i)
!               IF(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i)).GT.0)
!     >           write(*,*)'nn found'
            END DO
            CALL Get_Prob(0,0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO
!  Filled Coverage/walls Sites
      DO z=Lenww,LenZ-7
       DO y=0,LenY-1
        DO x=0,INT(0.5*(LenX-Laperture))-1
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        END DO
        DO x=INT(0.5*(LenX+Laperture)),LenX-1
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        END DO
       END DO
      END DO
!  Filled wall Sites at z=LenZ (helps with migrating vacancies issue)
      DO z=LenZ-6,LenZ-1
       DO y=0,LenY-1
        DO x=0,LenX-1
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        END DO
       END DO
      END DO     
      END SUBROUTINE SetSiC3Cwwaperture

**************************************************************************************
