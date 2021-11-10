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
      SUBROUTINE SetAFB(Len)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len
      INTEGER i,j,x,y,z,x1,y1,z1,Bon,Jum,Coo,Coor,CoorNN,Occ,OccN,IAt
      INTEGER Bon1
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NBuff(3)

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      DO z=0,Len-4
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN
            IF(x+y.LE.LenY)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
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
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF(x+y.LE.LenY)THEN  ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
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
            IF(x+y.LE.LenY)THEN ! Si sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
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
!            Prob = PtransE(1,2,0)! Si atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            IF(x+y.LE.LenY)THEN  ! C sites
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            END IF
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
!            Prob = PtransE(2,2,0)! C atom on SiC Surface
!            Index_Event = 3
!            CALL AddMC(Site,NextN,Index_Event,Prob)
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
 !           CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
 !           CALL AddMC(Site,NextN,Index_Event,Prob)
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN
            Bon=2
            CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
            Site(1)=x
            Site(2)=y
            Site(3)=z
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
!            CALL Get_Prob(0,0,2,Site,NextN,Index_Event,Prob) ! Occ= 0 NSiC =0 Coor=2
!            CALL AddMC(Site,NextN,Index_Event,Prob)
           END IF
          END IF
         END IF
        END DO
       END DO
       !write(*,*)"NumAtoms,NumAdAtom",NumAtoms,NumAdAtom
      END DO

! Set the probability
      DO z=Len-3,Len
       DO y=0,LenY-1
        DO x=0,LenX-1
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
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminate NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! Si sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             Site(1)=x
             Site(2)=y
             Site(3)=z
             DO i=1,4
               NextN(1:3,i) = Site + JumpDia(1:3,i+4)
               IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
               IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
               IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
               IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             END DO
             Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
             Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
             IF(Occ.NE.0.AND.Coor.EQ.0)write(*,*)'error',x,y,z,Occ,Coor
             IF(Occ.EQ.0)THEN ! Not Occupied
               IF(Coor.GT.0)THEN ! MC Particle
                 IF(Coor.EQ.1)THEN  ! Eliminete NN not occupied
                   DO i=1,4
                    OccN=IBITS(LattCoo(NextN(1,i),NextN(2,i),
     >                         NextN(3,i)),PosOcc,LenOcc)
                    IF(OccN.EQ.0)THEN
                     NextN(1:3,i)=0
                    ELSE
                     NBuff(1:3)= NextN(1:3,i)
                     NextN(1:3,i)=0
                    END IF
                   END DO
                   NextN(1:3,1)=NBuff(1:3)
                 END IF
                 CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               END IF
             ELSE ! Occupied
               IF(Coor.LE.3)THEN
                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                 IF(x+y.LE.LenY)THEN ! C sites
                   CALL Get_Prob_Ini(2,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 ELSE ! Si Sites
                   CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                 END IF
                 CALL AddMC(Site,NextN,Index_Event,Prob)
               ELSE
                 CALL AddBulk(Site,NextN)
               END IF
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO

      END SUBROUTINE SetAFB
