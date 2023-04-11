**************************************************************************************
**   Set a triple SF in 100 sustrate with triangular shape diplacing a set of atoms in 111 plane
**   Len z=cost position of the 100 surface
**   VSF position of the triangle's vertex in the substrate bulk (one Si sites of the conventional cubic
**   cells)
**************************************************************************************
      SUBROUTINE SetSiGeSF(Len,VSF,LenOffset,a)
      USE DefSystem
      USE Definitions
      IMPLICIT NONE
      INTEGER Len, LenNew, LenOffset
      INTEGER i,j,k,x,y,z,x1,y1,z1,Bon,Jum,Coor,Coord,CoorNN,IndNN,CovInd
      INTEGER BoxSF,plan1Si,plan2Si,plan1C,plan2C,SiOrC,DSN,Ichange
      INTEGER Bon1,nlist,Occ,OccN,IAt,NumAdAtomold,IndS,Coorold
      INTEGER SiteSiC(3),Site(3),SiteGr(3),NNSite(3)
      INTEGER VSF(3)
      INTEGER NextN(3,4),NextNN(3,4),NextND(3,4),Buf(3)
      INTEGER ListNN(3,32),NewSites(3,2)
      LOGICAL VPos,Check,Recipr

      INTEGER :: Index_Event, indsp, Id
      REAL(8) :: Prob, rr
      INTEGER :: temperature

      INTEGER, DIMENSION(LenZ,LenX*LenY) :: a

      IF(alloyfraction.GT.0)THEN
        write(*,*) 'Setting up a', LenSiGe, 'Angstrom-thick SiGe layer with fraction = ', alloyfraction
      ELSE
        write(*,*) 'Setting up Silicon'
      END IF 
      write(*,*) 'a(1,1)', a(1,1)

!  Filled wall Sites at z=0 (helps with migrating vacancies issue)
      DO z=0,3
       DO y=0,LenY-1
        DO x=0,LenX-1
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        END DO
       END DO
      END DO      

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
      IF(LenOffset.EQ.0)THEN
        IF(VSF(1)+BoxSF.GT.Lenx.OR.VSF(2)+BoxSF.GT.Leny.
     >     OR.BoxSF.LE.0)THEN
          write(*,*)'error SF dimension is not correctly set'
          STOP
        END IF
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

      IF(INT(LenNuc/sf).GT.LenOffset)THEN
        write(*,*)'LenNuc cannot be larger than LenOffSet'
        STOP
      END IF

      write(*,*)'NumAtoms 1', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0,SUM(CountCrystal) ! <--------------------


!  First layer Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      DO z=4,12
       DO y=0,LenY-1
        DO x=0,LenX-1

         ! Check geom stored in socket
         Id = a(1+z, 1+x+LenX*y)  ! this is the correct one!
         
         IF(Id.EQ.10)THEN ! WALL, AIR, SiO2
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
         
         ELSE IF(Id.EQ.1)THEN ! Set species 
         
          IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
           x1=x/3
           y1=y/3
           z1=z/3
           IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
            IF((MOD(x1+y1+z1,4).EQ.0).OR.(MOD(x1+y1+z1,4).EQ.3))THEN ! Si/Ge sites

             ! Increase number of solids existing before nucleation
             NumOcc0=NumOcc0+1  

             IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
               rr = 1.0   ! choose Si
             ELSE
               rr = RAND()
             END IF
             IF(rr.GE.alloyfraction)THEN
                 LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                 CountCrystal(1) = CountCrystal(1) + 1
             ELSE
                 LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                 CountCrystal(2) = CountCrystal(2) + 1
             END IF
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
            END IF
           END IF
          END IF
         
         END IF
        END DO
       END DO
      END DO


      write(*,*)'NumAtoms 2', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0,SUM(CountCrystal) ! <--------------------

!  Filled Bulk Sites in 100 SiC lattice and a triple ST
      LenNew = Len-LenOffset 

      DO z=13,LenNew-1
       DO y=0,LenY-1
        DO x=0,LenX-1

         ! Check geom stored in socket
         Id = a(1+z, 1+x+LenX*y)  ! this is the correct one!
         
         IF(Id.EQ.10)THEN ! WALL, AIR, SiO2
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
         
         ELSE IF(Id.EQ.1)THEN ! Set species 

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
               IF(x+2.LT.Lenx.AND.y+2.LT.Leny.AND.z+4.LT.LenNew)THEN
                
                ! Increase number of solids existing before nucleation
                NumOcc0=NumOcc0+1  
                IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                  rr = 1.0   ! choose Si
                ELSE
                  rr = RAND()
                END IF
                IF(rr.GE.alloyfraction)THEN
                    LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosSi)
                    CountCrystal(1) = CountCrystal(1) + 1
                ELSE
                    LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosC)
                    CountCrystal(2) = CountCrystal(2) + 1
                END IF
                LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosOcc)
               END IF
              
              ELSE IF(plan2Si.EQ.0)THEN
               IF(x-2.GT.0.AND.y-2.GT.0.AND.z-4.GT.0)THEN
                
                ! Increase number of solids existing before nucleation
                NumOcc0=NumOcc0+1  
                IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                  rr = 1.0   ! choose Si
                ELSE
                  rr = RAND()
                END IF
                IF(rr.GE.alloyfraction)THEN
                    LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosSi)
                    CountCrystal(1) = CountCrystal(1) + 1
                ELSE
                    LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosC)
                    CountCrystal(2) = CountCrystal(2) + 1
                END IF
                LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosOcc)
               END IF

              ELSE
               
               ! Increase number of solids existing before nucleation
               NumOcc0=NumOcc0+1  
               IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                 rr = 1.0   ! choose Si
               ELSE
                 rr = RAND()
               END IF
               IF(rr.GE.alloyfraction)THEN
                   LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                   CountCrystal(1) = CountCrystal(1) + 1
               ELSE
                   LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                   CountCrystal(2) = CountCrystal(2) + 1
               END IF
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
               IF(x+2.LT.Lenx.AND.y+2.LT.Leny.AND.z+4.LT.LenNew)THEN
                
                ! Increase number of solids existing before nucleation
                NumOcc0=NumOcc0+1  
                IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                  rr = 1.0   ! choose Si
                ELSE
                  rr = RAND()
                END IF
                IF(rr.GE.alloyfraction)THEN
                    LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosSi)
                    CountCrystal(1) = CountCrystal(1) + 1
                ELSE
                    LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosC)
                    CountCrystal(2) = CountCrystal(2) + 1
                END IF
                LattCoo(x+2,y+2,z+4)=IBSET(LattCoo(x+2,y+2,z+4),PosOcc)
               END IF
              ELSE IF(plan2C.EQ.0)THEN
               IF(x-2.GT.0.AND.y-2.GT.0.AND.z-4.GT.0)THEN
                
                ! Increase number of solids existing before nucleation
                NumOcc0=NumOcc0+1  
                IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                  rr = 1.0   ! choose Si
                ELSE
                  rr = RAND()
                END IF
                IF(rr.GE.alloyfraction)THEN
                    LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosSi)
                    CountCrystal(1) = CountCrystal(1) + 1
                ELSE
                    LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosC)
                    CountCrystal(2) = CountCrystal(2) + 1
                END IF
                LattCoo(x-2,y-2,z-4)=IBSET(LattCoo(x-2,y-2,z-4),PosOcc)
               END IF

              ELSE
               
               ! Increase number of solids existing before nucleation
               NumOcc0=NumOcc0+1  
               IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                 rr = 1.0   ! choose Si
               ELSE
                 rr = RAND()
               END IF
               IF(rr.GE.alloyfraction)THEN
                   LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                   CountCrystal(1) = CountCrystal(1) + 1
               ELSE
                   LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                   CountCrystal(2) = CountCrystal(2) + 1
               END IF
               LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
              END IF
          
            END IF
           
           END IF
          END IF

         END IF
        
        END DO
       END DO
      END DO
      

      write(*,*)'NumAtoms 3', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0,SUM(CountCrystal) ! <--------------------

!  Top layer Filled Bulk Sites in 100 SiC lattice
!      write(*,*)'15 9 45',LattCoo(15,9,45)
      IF(LenOffset.GT.0)THEN
        DO z=LenNew,Len-1
         DO y=0,LenY-1
          DO x=0,LenX-1

           ! Check geom stored in socket
           Id = a(1+z, 1+x+LenX*y)  ! this is the correct one!
           
           IF(Id.EQ.10)THEN ! WALL, AIR, SiO2
            LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
            CovInd=113 
            CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
           
           ELSE IF(Id.EQ.1)THEN ! Set species 
           
            IF(MOD(x,3).EQ.0.AND.MOD(y,3).EQ.0.AND.MOD(z,3).EQ.0)THEN
             x1=x/3
             y1=y/3
             z1=z/3
             IF((MOD(x1,2).EQ.MOD(y1,2)).AND.(MOD(x1,2).EQ.MOD(z1,2)))THEN
              IF((MOD(x1+y1+z1,4).EQ.0).OR.(MOD(x1+y1+z1,4).EQ.3))THEN ! Si/Ge sites

               ! Increase number of solids existing before nucleation
               NumOcc0=NumOcc0+1  

               ! here we nucleate homogeneously a thin liquid layer
               IF(z.LE.(Len-INT(LenNuc/sf)))THEN 

                 IF(z.LE.(Len-INT((LenSiGe)/sf)))THEN
                   rr = 1.0   ! choose Si
                 ELSE
                   rr = RAND()
                 END IF
                 IF(rr.GE.alloyfraction)THEN
                     LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosSi)
                     CountCrystal(1) = CountCrystal(1) + 1
                 ELSE
                     LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosC)
                     CountCrystal(2) = CountCrystal(2) + 1
                 END IF
                 LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
               
                 Bon=4 ! coordination
                 CALL MVBITS(Bon,0,LenCoor,LattCoo(x,y,z),PosCoor)
C                  Site(1)=x
C                  Site(2)=y
C                  Site(3)=z
C                  DO i=1,4
C                     NextN(1:3,i) = Site + JumpDia(1:3,i)
C                     IF(NextN(1,i).GE.LenX)NextN(1,i)=NextN(1,i)-LenX
C                     IF(NextN(1,i).LT.0)NextN(1,i)=NextN(1,i)+LenX
C                     IF(NextN(2,i).GE.LenY)NextN(2,i)=NextN(2,i)-LenY
C                     IF(NextN(2,i).LT.0)NextN(2,i)=NextN(2,i)+LenY
C                  END DO
C                  CALL AddBulk(Site,NextN)
!                ! Probabilities and AddBulk will be set further below, once the correct neighbours are found                  

               END IF

              END IF
             END IF
            END IF
           
           END IF
          
          END DO
         END DO
        END DO
      END IF

      write(*,*)'NumAtoms 4', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0,SUM(CountCrystal) ! <--------------------


! not efficient loop to eliminate occupied sites with coor less and equal to one
      DO z=13,Len+6
       DO y=0,LenY-1
        DO x=0,LenX-1
        Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
        IF(Occ.NE.0)THEN
         CovInd = IBITS(LattCoo(x,y,z),PosIndex,LenIndex) ! WALL 
         IF(CovInd.NE.113)THEN
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

           IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
           IF(Coord.GT.Coor)THEN
             DO i=1,j-1
               DO k=i+1,j
               DSN=  (NextND(1,i)-NextND(1,k))**2+
     >               (NextND(2,i)-NextND(2,k))**2+
     >               (NextND(3,i)-NextND(3,k))**2
               IF(DSN.NE.72)THEN
                 write(*,*)'not correct bonds direction 1'
                 write(*,*)NextND(1:3,i),NextND(1:3,k)
               END IF
               END DO
             END DO
             IF(Coord.LE.1)THEN
               IF(InitSt.EQ.'LA')THEN
                 temperature = IBITS(LattCoo(x,y,z),PosT,LenT) ! K
                 LattCoo(x,y,z)=0
                 CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
               ELSE  
                 LattCoo(x,y,z)=0
               END IF
               IF(IAt.GT.0)THEN ! <--------------------------------------------------------------------------
                CountCrystal(IAt) = CountCrystal(IAt) - 1
                NumOcc0 = NumOcc0 - 1
               END IF

             END IF
           ELSE
             IF(Coor.LE.1)THEN
               IF(InitSt.EQ.'LA')THEN
                 temperature = IBITS(LattCoo(x,y,z),PosT,LenT) ! K
                 LattCoo(x,y,z)=0
                 CALL MVBITS(temperature,0,LenT,LattCoo(x,y,z),PosT)
               ELSE  
                 LattCoo(x,y,z)=0
               END IF
               IF(IAt.GT.0)THEN  ! <--------------------------------------------------------------------------
                CountCrystal(IAt) = CountCrystal(IAt) - 1
                NumOcc0 = NumOcc0 - 1
               END IF
             END IF
           END IF
         END IF
        END IF
        END DO
       END DO
      END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Don't move this block below!!! CountLiquid is used in get_prob so it should be defined before its first call!
      DO k=1,NCrystal
        IF(ListCrystal(k).EQ.14) THEN
          InitCountCrystal(k) = CountCrystal(k) + INT((NumOcc0 - SUM(CountCrystal))*(1-alloyfraction))
        ELSE IF(ListCrystal(k).EQ.32) THEN
          InitCountCrystal(k) = CountCrystal(k) + INT((NumOcc0 - SUM(CountCrystal))*alloyfraction)
        END IF
      END DO
      CountLiquid = InitCountCrystal - CountCrystal 
      write(*,*)' '
      write(*,*)'SetSiGeSF.f > Initial number of occupied sites (before nucleation): ', NumOcc0 ! Counted inside SetCAD.f
      write(*,*)'SetSiGeSF.f > Initial CountCrystal (before nucleation): ', InitCountCrystal ! Counted inside SetCAD.f
      write(*,*)'SetSiGeSF.f > Initial number of occupied sites (after nucleation): ', SUM(CountCrystal)
      write(*,*)'SetSiGeSF.f > Initial CountCrystal (after nucleation): ', CountCrystal ! Counted inside SetCAD.f
      write(*,*)'SetSiGeSF.f > Initial number of liquid sites (after nucleation): ', SUM(CountLiquid)
      write(*,*)'SetSiGeSF.f > Initial CountLiquid (after nucleation): ', CountLiquid ! Counted inside SetCAD.f
      write(*,*)' '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      write(*,*)'NumAtoms 5', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0 ! <--------------------


      DO z=13,Len+6
       DO y=0,LenY-1
        DO x=0,LenX-1
        Occ = IBITS(LattCoo(x,y,z),PosOcc,LenOcc)
        IF(Occ.NE.0)THEN
         CovInd = IBITS(LattCoo(x,y,z),PosIndex,LenIndex) ! WALL 
         IF(CovInd.NE.113)THEN
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
                   write(*,*)'not correct bonds direction 2'
                   write(*,*)NextND(1:3,i),NextND(1:3,k)
                  END IF
                 END DO
                END DO
                IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                CALL Get_Prob_Ini(1,IAt,Coord,Site,NextND,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL AddMC(Site,NextND,Index_Event,Prob)
                NumAdAtomOcc = NumAdAtomOcc + 1   ! Update Count of occupied AdAtoms
             ELSE
                CALL AddBulk(Site,NextND) ! This if is to avoid double counting of the top layer in case of LenOffset != 0
                !IF(z.LE.LenNew-1) CALL AddBulk(Site,NextND) ! This if is to avoid double counting of the top layer in case of LenOffset != 0
             END IF
           ELSE
             CALL MVBITS(Coor,0,LenCoor,LattCoo(x,y,z),PosCoor)
             IF(Coor.LE.3)THEN
!                IF(Coor.EQ.2)write(*,*)x,y,z,Coor
                IAt=IBITS(LattCoo(x,y,z),PosSiC,LenSiC)
                CALL Get_Prob_Ini(1,IAt,Coor,Site,NextN,
     >                     Index_Event,Prob) ! Occ=1/0 NSiC 1=01 Si, 2=10 C
                CALL AddMC(Site,NextN,Index_Event,Prob)
                NumAdAtomOcc = NumAdAtomOcc + 1   ! Update Count of occupied AdAtoms
             ELSE
                CALL AddBulk(Site,NextN) ! This if is to avoid double counting of the top layer in case of LenOffset != 0
                !IF(z.LE.LenNew-1) CALL AddBulk(Site,NextN) ! This if is to avoid double counting of the top layer in case of LenOffset != 0
             END IF
           END IF
         END IF
        END IF
        END DO
       END DO
      END DO
      Ichange=0

      write(*,*)'NumAtoms 6', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0 ! <--------------------

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
           CovInd = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),PosIndex,LenIndex) ! WALL 
           IF(CovInd.NE.113)THEN
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

      write(*,*)'NumAtoms 7', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0 ! <--------------------

      ichange=0
      DO i=1,NumAtoms
        Site = ListAtom(i) % AtomXYZ
        x=Site(1)
        y=Site(2)
        z=Site(3)
!        IF((z.GE.13).AND.(z.LE.Len-1))THEN ! This if is to avoid double counting of the top layer in case of LenOffset != 0
!        IF((z.GE.13).AND.(z.LE.LenNew-1))THEN ! This if is to avoid double counting of the top layer in case of LenOffset != 0
        IF(z.GE.13)THEN ! This if is to avoid double counting of the top layer in case of LenOffset != 0
          NextN = ListAtom(i) % NextNXYZ
          Coor = IBITS(LattCoo(x,y,z),PosCoor,LenCoor)
          Coorold=Coor
          DO j=1,4
            OccN = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),
     >                  PosOcc,LenOcc)
          !write(*,*)"Ev. NextN,Occ,CoorNN",NextN(1:3,i),Occ,CoorNN
            IF(OccN.NE.0)THEN
             CovInd = IBITS(LattCoo(NextN(1,j),NextN(2,j),NextN(3,j)),PosIndex,LenIndex) ! WALL 
             IF(CovInd.NE.113)THEN
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
            NumAdAtomOcc = NumAdAtomOcc + 1   ! Update Count of occupied AdAtoms 

  !         ListAdAtom(IndS) % Ind_Event = Index_Event  ! redundant
  !          ListAdAtom(IndS) % ProbTrans = Prob
  !          CALL Updatetree(IndS,Prob) !
          END IF
        END IF
      END DO
!      write(*,*)NumAtoms,ichange  

      write(*,*)'NumAtoms 8', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0 ! <--------------------

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
      
!  Filled wall Sites at z=LenZ (helps with migrating vacancies issue)
      DO z=LenZ-4,LenZ-1
       DO y=0,LenY-1
        DO x=0,LenX-1
          LattCoo(x,y,z)=IBSET(LattCoo(x,y,z),PosOcc)
          CovInd=113 ! coverage (for now it indicates walls)
          CALL MVBITS(CovInd,0,LenIndex,LattCoo(x,y,z),PosIndex)
        END DO
       END DO
      END DO      


      write(*,*)'NumAtoms 9', NumAtoms,NumAdAtom,NumAdAtomOcc,NumOcc0 ! <--------------------

      END SUBROUTINE SetSiGeSF

**************************************************************************************
