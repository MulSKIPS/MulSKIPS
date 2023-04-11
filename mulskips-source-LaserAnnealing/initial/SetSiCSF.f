**************************************************************************************
**   Set a triple SF in 100 sustrate with triangular shape diplacing a set of atoms in 111 plane
**   Len z=cost position of the 100 surface
**   VSF position of the triangle's vertex in the substrate bulk (one Si sites of the conventional cubic
**   cells)
**************************************************************************************
      SUBROUTINE SetSiCSF(Len,VSF)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len
      INTEGER i,j,k,x,y,z,x1,y1,z1,Bon,Jum,Coor,Coord,CoorNN,IndNN
      INTEGER BoxSF,plan1Si,plan2Si,plan1C,plan2C,SiOrC,DSN,Ichange
      INTEGER Bon1,nlist,Occ,OccN,IAt,NumAdAtomold,IndS,Coorold
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER VSF(3)
      INTEGER NextN(3,4),NextNN(3,4),NextND(3,4),Buf(3)
      INTEGER ListNN(3,32),NewSites(3,2)
      LOGICAL VPos,Check,Recipr

      INTEGER :: Index_Event
      REAL(8) :: Prob
!  Consistency check
      BoxSF = Len-VSF(3)
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
      IF(VSF(1)+BoxSF.GT.Lenx.OR.VSF(2)+BoxSF.GT.Leny.
     >   OR.BoxSF.LE.0)THEN
        write(*,*)'error SF dimension is not correctly set'
        STOP
      END IF
      Vpos=.FALSE.
      IF(MOD(VSF(1),3).EQ.0.AND.MOD(VSF(2),3).EQ.0.
     >   AND.MOD(VSF(3),3).EQ.0)THEN
        x1=VSF(1)/3
        y1=VSF(2)/3
        z1=VSF(3)/3
        IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
         IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
           Vpos=.TRUE.
         END IF
        END IF
      END IF
      IF(Vpos)THEN
         write(*,*)'SF Vertex position in a Si Site',VSF,Len
      ELSE
         write(*,*)'SF Vertex position in wrong Site',VSF,Len
         STOP
      END IF
!  First layer Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      DO z=0,12
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

!  Filled Bulk Sites in 100 SiC lattice and a triple ST
      DO z=13,Len-1
       DO y=0,LenY-1
        DO x=0,LenX-1
         IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
          x1=x/3
          y1=y/3
          z1=z/3
          IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
           IF(MOD(x1+y1+z1,4).EQ.0)THEN ! Si sites
             plan1Si=-1
             plan2Si=-1
             IF(x.GE.VSF(1).AND.y.GE.VSF(2).AND.z.GE.VSF(3))THEN
               plan1Si = (x-VSF(1))+(y-VSF(2))-(z-VSF(3))-12
               plan2Si = (x-VSF(1))+(y-VSF(2))-(z-VSF(3))-24
             END IF
             IF(plan1Si.EQ.0)THEN
              IF(x+2.LT.Lenx.AND.y+2.LT.Leny.AND.z+4.LT.Len)THEN
               LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosSi)
               LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosOcc)
              END IF
             ELSE IF(plan2Si.EQ.0)THEN
              IF(x-2.GT.0.AND.y-2.GT.0.AND.z-4.GT.0)THEN
               LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosSi)
               LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosOcc)
              END IF
             ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           ELSE IF(MOD(x1+y1+z1,4).EQ.3)THEN ! C sites
             plan1C=-1
             plan2C=-1
             IF(x.GE.VSF(1).AND.y.GE.VSF(2).AND.z.GE.VSF(3))THEN
               plan1C  = (x-VSF(1))+(y-VSF(2))-(z-VSF(3))-3
               plan2C  = (x-VSF(1))+(y-VSF(2))-(z-VSF(3))-15
             END IF
             IF(plan1C.EQ.0)THEN
              IF(x+2.LT.Lenx.AND.y+2.LT.Leny.AND.z+4.LT.Len)THEN
               LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosC)
               LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosOcc)
              END IF
             ELSE IF(plan2C.EQ.0)THEN
              IF(x-2.GT.0.AND.y-2.GT.0.AND.z-4.GT.0)THEN
               LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosC)
               LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosOcc)
              END IF
             ELSE
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
              LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
             END IF
           END IF
          END IF
         END IF
        END DO
       END DO
      END DO
      

! not efficient loop to eliminate occupied sites with coor less and equal to one
      DO z=13,Len+6
       DO y=0,LenY-1
        DO x=0,LenX-1
        Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
        IF(Occ.NE.0)THEN
           Coor=0
           Coord=0
           Site(1)=x
           Site(2)=y
           Site(3)=z
           x1=x/3
           y1=y/3
           z1=z/3
           SiOrC=0
           IF(MOD(x1+y1+z1,4).EQ.3)SiOrC=4
           DO i=1,4 ! Diamond lattice only check
             NextN(1:3,i) = Site + JumpDia(1:3,i+SiOrC)
             IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
             IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
             IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
             IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
             IF(OccN.GT.0)Coor=Coor+1
           END DO
           IF(Coor.LT.2)THEN ! Possible defective site
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
           END IF
           IF(Coord.GT.Coor)THEN
             DO i=1,j-1
               DO k=i+1,j
               DSN=  (NextND(1,i)-NextND(1,k))**2+
     >               (NextND(2,i)-NextND(2,k))**2+
     >               (NextND(3,i)-NextND(3,k))**2
               IF(DSN.NE.72)THEN
                 write(*,*)'not correct bonds direction'
                 write(*,*)NextND(1:3,i),NextND(1:3,k)
               END IF
               END DO
             END DO
             IF(Coord.LE.1)LattCoo(x,y,z)=0
           ELSE
             IF(Coor.LE.1)LattCoo(x,y,z)=0
           END IF
        END IF
        END DO
       END DO
      END DO

      DO z=13,Len+6
       DO y=0,LenY-1
        DO x=0,LenX-1
        Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
        IF(Occ.NE.0)THEN
           Coor=0
           Coord=0
           Site(1)=x
           Site(2)=y
           Site(3)=z
           x1=x/3
           y1=y/3
           z1=z/3
           SiOrC=0
           IF(MOD(x1+y1+z1,4).EQ.3)SiOrC=4
           DO i=1,4 ! Diamond lattice only check
             NextN(1:3,i) = Site + JumpDia(1:3,i+SiOrC)
             IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX ! No Action If Len1<<LenX
             IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX    ! No Action If Len1<<LenX
             IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY ! No Action If Len2<<LenY
             IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY    ! No Action If Len2<<LenY
             OccN = IBITS(LattCoo(NextN(1,i),NextN(2,i),NextN(3,i))
     >                     ,PosOcc,LenOcc)
             IF(OccN.GT.0)Coor=Coor+1
           END DO
           IF(Coor.LT.2)THEN ! Possible defective site
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
           END IF
           IF(Coord.GT.Coor)THEN
C              write(*,*)Site,Coor,Coord
             CALL MVBITS(Coord,0,LenCoor,LattCoo(x,y,z),PosCoor)
             IF(Coord.LE.3)THEN
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
                DO i=1,3
                 DO k=i+1,4
                  DSN=  (NextND(1,i)-NextND(1,k))**2+
     >               (NextND(2,i)-NextND(2,k))**2+
     >               (NextND(3,i)-NextND(3,k))**2
                  IF(DSN.NE.72)THEN
                   write(*,*)'not correct bonds direction'
                   write(*,*)NextND(1:3,i),NextND(1:3,k)
                  END IF
                 END DO
                END DO
                IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                CALL Get_Prob_Ini(1,IAt,Coord,Site,NextND,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C

                CALL AddMC(Site,NextND,Index_Event,Prob)
             ELSE
                CALL AddBulk(Site,NextND)
             END IF
           ELSE
             CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             IF(Coor.LE.3)THEN
!                IF(Coor.EQ.2)write(*,*)x,y,z,Coor
                IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL AddMC(Site,NextN,Index_Event,Prob)
             ELSE
                CALL AddBulk(Site,NextN)
             END IF
           END IF
        END IF
        END DO
       END DO
      END DO
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
      END SUBROUTINE SetSiCSF

**************************************************************************************
