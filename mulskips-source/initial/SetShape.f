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
      SUBROUTINE SetShape(NFilled)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER i,j,x,y,z,x1,y1,z1,x0,y0,z0
      INTEGER Id,CovInd,OccN,Occ,IAt,Coor,Coord,CoorNN,Coorold,ia
      INTEGER SiOrC,DSN,Ichange,NumAdAtomold,IndS,IndNN,Nfilled,nlist
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER NextN(3,4),NextNN(3,4),NextND(3,4),Buf(3)
      INTEGER ListNN(3,32), NewSites(3,2) 
      LOGICAL Check,Recipr
      INTEGER ListXYZ(3,NFilled)
      CHARACTER(Len=2) :: sym

      INTEGER :: Index_Event
      REAL(8) :: Prob

      INTEGER :: temperature
      REAL(8) :: rr, LenSiGe
      INTEGER :: k

      INTEGER :: IPF50, IPF60



! Read and set occupations
      IPF50=50
      write(*,*)'Reading input geometry file (SHAPE): ', cadfilename
      OPEN(IPF50,FILE=cadfilename,STATUS='OLD')     
      ! read number of elements to be filled
      READ(IPF50,*)NFilled
      ! skip second line
      READ(IPF50,*)
      ! Now start reading indices and filling with Si or C           
      DO ia=1,NFilled
        READ(IPF50,*)sym,x,y,z
        LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
        !  Filled wall Sites at z=0 and z=LenZ (solves boundary vacancies bug and overwrites geo setting)
        IF(z.GE.(LenZ-4).OR.z.LE.3)THEN
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        ELSE            
          ! Otherwise set Si or C/Ge
          IF(sym.EQ.'Si')THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
          ELSE IF(sym.EQ.'C'.OR.sym.EQ.'Ge')THEN
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
          END IF
        END IF
        ! Store xyz for later usage
        ListXYZ(1,ia)=x
        ListXYZ(2,ia)=y
        ListXYZ(3,ia)=z 
      END DO
      CLOSE(IPF50)
      write(*,*)'Done reading input geometry file'




      nlist=0
      DO k=-6,6
        DO j=-6,6
          DO i=-6,6
           IF(i*i+j*j+k*k.EQ.27)THEN
             nlist=nlist+1
             ListNN(1,nlist)=i
             ListNN(2,nlist)=j
             ListNN(3,nlist)=k
           END IF
          END DO
        END DO
      END DO




      DO z=4,LenZ-5
       DO y=0,LenY-1
        DO x=0,LenX-1
          Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
          CovInd=IBITS(LattCoo(x,y,z),PosIndex,LenIndex)

          ! First set coord and prob on all occupied sites 
          IF(Occ.NE.0)THEN
             Coor=0
             Coord=0
             Site(1)=x
             Site(2)=y
             Site(3)=z

             ! Find nearest neighbors by checking the entire cubic/hex shell 
             ! NB: I've removed 0- and 1-coordinated ones in SISL, so Coord can only be 2,3,4 ! 
             j=0
             NextND=0
             DO i=1,32 !full NN list check
               NNSite=Site+ListNN(1:3,i)
               IF(NNSite(1).GE.LenX)NNSite(1)=NNSite(1)-LenX ! No Action If Len1<<LenX
               IF(NNSite(1).LT.0)NNSite(1)=NNSite(1)+LenX    ! No Action If Len1<<LenX
               IF(NNSite(2).GE.LenY)NNSite(2)=NNSite(2)-LenY ! No Action If Len2<<LenY
               IF(NNSite(2).LT.0)NNSite(2)=NNSite(2)+LenY    ! No Action If Len2<<LenY
               OccN = IBITS(LattCoo(NNSite(1),NNSite(2),NNSite(3))
     >                     ,PosOcc,LenOcc)
               IF(OccN.GT.0)THEN
                 Coord=Coord+1
                 j=j+1
                 NextND(1:3,j)=NNSite
               END IF
             END DO

             CALL MVBITS(Coord,0,LenCoor,LattCoo(x,y,z),PosCoor)

             ! Set prob 
             IF(Coord.LE.3)THEN
              
                ! Various checks
                Buf=NextND(1:3,3)
                CALL FIND_NN(NextND(1:3,1),NextND(1:3,2),Site,NewSites,
     >                       LenX,LenY)
                NextND(1:3,3)=NewSites(1:3,1)
                NextND(1:3,4)=NewSites(1:3,2)
                IF(Coord.EQ.3)THEN
                 Check=.TRUE.
                 DO i=1,4
                   IF(ALL(NextND(1:3,i).EQ.Buf))THEN
                    Check=.FALSE.
                   END IF
                 END DO
                 IF(Check)write(*,*)'Coord.eq.3 wrong check Buf',Buf
                END IF
C                 DO i=1,3
C                  DO k=i+1,4
C                   DSN=  (NextND(1,i)-NextND(1,k))**2+
C      >               (NextND(2,i)-NextND(2,k))**2+
C      >               (NextND(3,i)-NextND(3,k))**2
C                   IF(DSN.NE.72)THEN
C                    write(*,*)'not correct bonds direction'
C                    write(*,*)NextND(1:3,i),NextND(1:3,k),Site,Coord,
C      >                             CovInd
C                   END IF
C                  END DO
C                 END DO

                IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                CALL Get_Prob_Ini(1,IAt,Coord,Site,NextND,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C

                CALL AddMC(Site,NextND,Index_Event,Prob)

             ELSE
                CALL AddBulk(Site,NextND)
             END IF

          END IF

        END DO
       END DO
      END DO








C ! not efficient loop to eliminate occupied sites with coor less and equal to one
C       DO z=0,LenZ-1
C        DO y=0,LenY-1
C         DO x=0,LenX-1
C         Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
C         IF(Occ.NE.0)THEN
C            Coor=0
C            Coord=0
C            Site(1)=x
C            Site(2)=y
C            Site(3)=z
C            x1=x/3
C            y1=y/3
C            z1=z/3
C            SiOrC=0
C            IF(MOD(x1+y1+z1,4).EQ.3)SiOrC=4
C            DO i=1,4 ! Diamond lattice only check
C              NextN(1:3,i) = Site + JumpDia(1:3,i+SiOrC)
C              IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
C              IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
C              IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
C              IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
C              OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
C      >                     ,PosOcc,LenOcc)
C              IF(OccN.GT.0)Coor=Coor+1
C            END DO
C            IF(Coor.LT.2)THEN ! Possible defective site
C              j=0
C              NextND=0
C              DO i=1,32 !full NN list check
C                NNSite=Site+ListNN(1:3,i)
C                IF(NNSite(1).GE.LenX)NNSite(1)=NNSite(1)-LenX ! No Action If Len1<<LenX
C                IF(NNSite(1).LT.0)NNSite(1)=NNSite(1)+LenX    ! No Action If Len1<<LenX
C                IF(NNSite(2).GE.LenY)NNSite(2)=NNSite(2)-LenY ! No Action If Len2<<LenY
C                IF(NNSite(2).LT.0)NNSite(2)=NNSite(2)+LenY    ! No Action If Len2<<LenY
C                OccN = IBITS(LattCoo(NNSite(1),NNSite(2),NNSite(3))
C      >                     ,PosOcc,LenOcc)
C                IF(OccN.GT.0)THEN
C                  Coord=Coord+1
C                  j=j+1
C                  NextND(1:3,j)=NNSite
C                END IF
C              END DO
C            END IF
C            IF(Coord.GT.Coor)THEN
C              DO i=1,j-1
C                DO k=i+1,j
C                DSN=  (NextND(1,i)-NextND(1,k))**2+
C      >               (NextND(2,i)-NextND(2,k))**2+
C      >               (NextND(3,i)-NextND(3,k))**2
C                IF(DSN.NE.72)THEN
C                  write(*,*)'not correct bonds direction'
C                  write(*,*)NextND(1:3,i),NextND(1:3,k)
C                END IF
C                END DO
C              END DO
C              IF(Coord.LE.1)LattCoo(x,y,z)=0
C            ELSE
C              IF(Coor.LE.1)LattCoo(x,y,z)=0
C            END IF
C         END IF
C         END DO
C        END DO
C       END DO

C       DO z=0,LenZ-1
C        DO y=0,LenY-1
C         DO x=0,LenX-1
C         Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
C         IF(Occ.NE.0)THEN
C            Coor=0
C            Coord=0
C            Site(1)=x
C            Site(2)=y
C            Site(3)=z
C            x1=x/3
C            y1=y/3
C            z1=z/3
C            SiOrC=0
C            IF(MOD(x1+y1+z1,4).EQ.3)SiOrC=4
C            DO i=1,4 ! Diamond lattice only check
C              NextN(1:3,i) = Site + JumpDia(1:3,i+SiOrC)
C              IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
C              IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
C              IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
C              IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
C              OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
C      >                     ,PosOcc,LenOcc)
C              IF(OccN.GT.0)Coor=Coor+1
C            END DO
C            IF(Coor.LT.2)THEN ! Possible defective site
C              j=0
C              NextND=0
C              DO i=1,32 !full NN list check
C                NNSite=Site+ListNN(1:3,i)
C                IF(NNSite(1).GE.LenX)NNSite(1)=NNSite(1)-LenX ! No Action If Len1<<LenX
C                IF(NNSite(1).LT.0)NNSite(1)=NNSite(1)+LenX    ! No Action If Len1<<LenX
C                IF(NNSite(2).GE.LenY)NNSite(2)=NNSite(2)-LenY ! No Action If Len2<<LenY
C                IF(NNSite(2).LT.0)NNSite(2)=NNSite(2)+LenY    ! No Action If Len2<<LenY
C                OccN = IBITS(LattCoo(NNSite(1),NNSite(2),NNSite(3))
C      >                     ,PosOcc,LenOcc)
C                IF(OccN.GT.0)THEN
C                  Coord=Coord+1
C                  j=j+1
C                  NextND(1:3,j)=NNSite
C                END IF
C              END DO
C            END IF
C            IF(Coord.GT.Coor)THEN
C C              write(*,*)Site,Coor,Coord
C              CALL MVBITS(Coord,0,LenCoor,LattCoo(x,y,z),PosCoor)
C              IF(Coord.LE.3)THEN
C                 Buf=NextND(1:3,3)
C                 CALL FIND_NN(NextND(1:3,1),NextND(1:3,2),Site,NewSites,
C      >                       LenX,LenY)
C                 NextND(1:3,3)=NewSites(1:3,1)
C                 NextND(1:3,4)=NewSites(1:3,2)
C                 IF(Coord.EQ.3)THEN
C                  Check=.TRUE.
C                  DO i=1,4
C                    IF(ALL(NextND(1:3,i).EQ.Buf))THEN
C                     Check=.FALSE.
C                    END IF
C                  END DO
C                  IF(Check)write(*,*)'Coord.eq.3 wrong check Buf',Buf
C                 END IF
C                 DO i=1,3
C                  DO k=i+1,4
C                   DSN=  (NextND(1,i)-NextND(1,k))**2+
C      >               (NextND(2,i)-NextND(2,k))**2+
C      >               (NextND(3,i)-NextND(3,k))**2
C                   IF(DSN.NE.72)THEN
C                    write(*,*)'not correct bonds direction'
C                    write(*,*)NextND(1:3,i),NextND(1:3,k)
C                   END IF
C                  END DO
C                 END DO
C                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
C                 CALL Get_Prob_Ini(1,IAt,Coord,Site,NextND,
C      >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C

C                 CALL AddMC(Site,NextND,Index_Event,Prob)
C              ELSE
C                 CALL AddBulk(Site,NextND)
C              END IF
C            ELSE
C              CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
C              IF(Coor.LE.3)THEN
C !                IF(Coor.EQ.2)write(*,*)x,y,z,Coor
C                 IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
C                 CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
C      >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
C                 CALL AddMC(Site,NextN,Index_Event,Prob)
C              ELSE
C                 CALL AddBulk(Site,NextN)
C              END IF
C            END IF
C         END IF
C         END DO
C        END DO
C       END DO


! The 3 loops below are identify the actual advoids and adatoms, and are needed to set prob correctly on advoids 

      Ichange=0
! Restore Coordination for Reciprocal bonds only
      DO i=1,NumAdAtom
        Site = ListAdAtom(i) % AtomXYZ
        x=Site(1)
        y=Site(2)
        z=Site(3)
        NextN = ListAdAtom(i) % NextNXYZ
        Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
        Coorold=Coor
        DO j=1,4
          OccN = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),
     >                  PosOcc,LenOcc)
        !write(*,*)"Ev. NextN,Occ,CoorNN",NextN(1:3,i),Occ,CoorNN
          IF(OccN.NE.0)THEN
            CoorNN = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),
     >          PosCoor,LenCoor)
            IndNN = LattInd(NextN(1,j),NextN(2,j),NextN(3,j))
            IF(CoorNN.EQ.4)THEN
              NextNN = ListAtom(IndNN) % NextNXYZ
            ELSE
              NextNN = ListAdAtom(IndNN) % NextNXYZ
            END IF
            Recipr=.FALSE.
            DO k=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,k).EQ.Site(1:3)))Recipr=.TRUE.
            END DO
            IF(Recipr.EQV..FALSE.)Coor=Coor-1
          END IF
        END DO
        IF(Coor.NE.Coorold)THEN
          ichange=ichange+1
C           write(*,*)'change coor',Site,Coor,Coorold
          IndS=LattInd(x,y,z)
          IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
          CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
          CALL Get_Prob(1,IAt,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
          ListAdAtom(IndS) % Ind_Event = Index_Event  ! redundant
          ListAdAtom(IndS) % ProbTrans = Prob
          CALL Updatetree(IndS,Prob) !
        END IF
      END DO


! Now do the same, but for bulk atoms
C       write(*,*)NumAdAtom,ichange
      ichange=0
      DO i=1,NumAtoms
        Site = ListAtom(i) % AtomXYZ
        x=Site(1)
        y=Site(2)
        z=Site(3)
        IF(z.GE.13)THEN
        NextN = ListAtom(i) % NextNXYZ
        Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
        Coorold=Coor
        DO j=1,4
          OccN = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),
     >                  PosOcc,LenOcc)
        !write(*,*)"Ev. NextN,Occ,CoorNN",NextN(1:3,i),Occ,CoorNN
          IF(OccN.NE.0)THEN
            CoorNN = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),
     >          PosCoor,LenCoor)
            IndNN = LattInd(NextN(1,j),NextN(2,j),NextN(3,j))
            IF(CoorNN.EQ.4)THEN
              NextNN = ListAtom(IndNN) % NextNXYZ
            ELSE
              NextNN = ListAdAtom(IndNN) % NextNXYZ
            END IF
            Recipr=.FALSE.
            DO k=1,4 ! check is the Site belongs to the NN of Next
                IF(ALL(NextNN(1:3,k).EQ.Site(1:3)))Recipr=.TRUE.
            END DO
            IF(Recipr.EQV..FALSE.)Coor=Coor-1
          END IF
        END DO
        IF(Coor.NE.Coorold)THEN
          ichange=ichange+1
C           write(*,*)'change coor',Site,Coor,Coorold
          IndS=LattInd(x,y,z)
          IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
          CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
          CALL Get_Prob(1,IAt,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
          CALL EraseBulk(Site)
          CALL AddMC(Site,NextN,Index_Event,Prob)

!         ListAdAtom(IndS) % Ind_Event = Index_Event  ! redundant
!          ListAdAtom(IndS) % ProbTrans = Prob
!          CALL Updatetree(IndS,Prob) !
        END IF
        END IF
      END DO


! Now set prob on advoids
C       write(*,*)NumAtoms,ichange
! Set the Ad-Voids from the current ListAdAtom (containing proper ad-Atoms)
      NumAdAtomold=NumAdAtom
      DO i=1,NumAdAtomold
        NNSite = ListAdAtom(i) % AtomXYZ
        DO j=1,4
           Site = ListAdAtom(i) % NextNXYZ(1:3,j)
           x=Site(1)
           y=Site(2)
           z=Site(3)
           Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
           IF(Occ.EQ.0)THEN
              NextN = 0
              Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
              Coor = Coor + 1
              IF(Coor.EQ.1)THEN ! New MC Particle
                NextN(1:3,1) = NNSite
                CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
                CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL AddMC(Site,NextN,Index_Event,Prob)
              ELSE IF(Coor.EQ.2)THEN ! Already A MC
                IndS=LattInd(x,y,z)
                NextN = ListAdAtom(IndS)% NextNXYZ
                NextN(1:3,2) = NNSite
                CALL FIND_NN(NextN(1:3,1),NextN(1:3,2),Site,NewSites,
     >                       LenX,LenY)
                NextN(1:3,3)=NewSites(1:3,1)
                NextN(1:3,4)=NewSites(1:3,2)
                CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
                CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                ListAdAtom(IndS) % NextNXYZ = NextN
                ListAdAtom(IndS) % Ind_Event = Index_Event
                ListAdAtom(IndS) % ProbTrans = Prob
                CALL Updatetree(IndS,Prob) !
              ELSE IF(Coor.EQ.3)THEN ! Already A MC with NN fixed
                IndS=LattInd(x,y,z)
                NextN = ListAdAtom(IndS) % NextNXYZ
                CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
                CALL Get_Prob_Ini(0,0,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                ListAdAtom(IndS) % Ind_Event = Index_Event
                ListAdAtom(IndS) % ProbTrans = Prob
                CALL Updatetree(IndS,Prob) !
              ELSE
                 write(*,*)'Vacancy in the SF setting'
                 STOP
              END IF
           END IF
        END DO
      END DO      


      END SUBROUTINE SetShape

**************************************************************************************
