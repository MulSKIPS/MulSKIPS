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
      PROGRAM KLMCSiC3C
      USE Definitions
      USE DefSystem
      IMPLICIT NONE
      INTEGER :: i,j,k,l,in,irot,K1,Itrans
      INTEGER, PARAMETER :: k10 = selected_int_kind(10)
      INTEGER(kind=k10) :: Iter=1, IterOut=0
      INTEGER(kind=k10) :: MaxIter=2100000000
      INTEGER(kind=k10) :: OutMolMol=100000000
      INTEGER(kind=k10) :: framecounter=0
      REAL(8) :: MaxTime=100000000.0
      REAL(8) :: OutTime=10000000.0
      INTEGER :: Ind,IndTest
      INTEGER :: MGigaCycle=400,IGigaCycle=0,Kilo=1000,Mega=1000000
      INTEGER :: MaxZeta=10000
      INTEGER :: Len=LenZ-1,Len1=120,Len2=120,Len3=120,Len4=120
      INTEGER :: OPF16,OPF17,IQstat,Coo,Occ
      INTEGER :: LenOut(LenX*LenY)
      INTEGER :: LenSaveOut(LenX*LenY)
      REAL(8) :: angle ! degrees
      REAL(8) :: Ttop, Tbottom ! Kelvin
      INTEGER :: NFilled

! Random NumBer
      REAL(8) :: Pi=3.1415926535897932384626433832795028841971
      REAL(8)  :: Time=0.,Dt,Patom
      INTEGER :: Site(3),Coord,Nbon,Njum,Jump,ijump,CovInd
      INTEGER :: Site1(3),Site2(3),SiteC(3),Newsite(3,2)
      INTEGER :: SNSite(3,3),SNZigZag(3,3),SNArmChair(3,3)
      LOGICAL   MolmOutpStat
      REAL(8) Random
      EXTERNAL Random
! Executable statments
8     FORMAT(I4,I5,' ',8ES10.3)
9     FORMAT(I4,' ',8ES10.3)
10    FORMAT('Jump: ',3I3,'   from: ',3I4)
11    FORMAT(I7,I5,I3,' ',ES10.3)
12    FORMAT(I8,ES10.3,' ',2F8.3)
13    FORMAT('o_site:',3I4,'  bon:',I2,'  jum: ',I2, ' nn:',8I2)
14    FORMAT('n_site:',3I4,'  bon:',I2,'  jum: ',I2, ' nn:',8I2)
115   FORMAT('JumpSt:',12I9)

      LattCoo=0
      LattInd=0
      counter=0
      NumAtoms=0
      NumVoids=0
      NumAdAtom=0
      CountCrystal=0
      CountCov=0
     
      write(*,*)'KMC box size: ',LenX, LenY, LenZ
 
      ! ------------------------------------------- INPUT
      ! Read input start.dat
      CALL Fileopen()
      
      ! Initialize crystal and possible coverage species
      READ(IPF,*)NCrystal,NCov
      IF(Ncrystal.GT.NcrystalMax)THEN
        write(*,*)'ERR NOT IMPLEMENTED: NCrystal must be <=',NCrystalMax
        STOP
      END IF
      READ(IPF,*)ListCrystal(:NCrystal)
      IF((NCov.GT.0).AND.(NCov.LE.NCovMax))THEN
        write(*,*)'           ***** Simulating CVD *****'
        READ(IPF,*)ListCov(:NCov)
      ELSE IF(NCov.EQ.0)THEN
        write(*,*)'           ***** Simulating PVD or LA *****'
      ELSE
        write(*,*)'ERR NOT IMPLEMENTED: NCov must be <=',NCovMax
        STOP
      END IF
      write(*,*)'Crystal species (NCrystal=',NCrystal,') --> Z:',
     > ListCrystal(:NCrystal)
      IF(NCov.GT.0)THEN
        write(*,*)'Coverage species (NCov=',NCov,') --> Z:',
     >   ListCov(:NCov)
      END IF
      
      ! Initialization state ! ONLY CHAR ! 
      READ(IPF,*)InitSt
      IF(InitSt.EQ.'FF')THEN
        READ(IPF,*)Len1,Len2,Len3,Len4
      ELSE IF(InitSt.EQ.'ST')THEN
        READ(IPF,*)Len1,Len2,Len3,angle,Ttop,Tbottom
        write(*,*)'Temperature on top in Kelvin: ', Ttop
        write(*,*)'Temperature on bottom in Kelvin: ', Tbottom
        IF(Ttop.EQ.Tbottom)THEN
          fixedT = .TRUE.
        END IF
      ELSE IF(InitSt.EQ.'IN'.OR.InitSt.EQ.'SH')THEN
        READ(IPF,'(A)')cadfilename
      ELSE IF(InitSt.EQ.'LA')THEN
        write(*,*)'****** MulSKIPs run for laser annealing ******'
        IF(NCov.GT.0)THEN
          write(*,*)'ERR NOT IMPLEMENTED: NCov=0 for LA'
          STOP
        END IF
        READ(IPF,'(A)')cadfilename
        READ(IPF,'(A)')tempfilename
        READ(IPF,*)Tm ! Kelvin
        READ(IPF,*)P0 
        READ(IPF,*)sigma
        READ(IPF,*)Tflex ! Kelvin
        write(*,*)'Melting temperature in Kelvin: ', Tm
        READ(IPF,*)LenVac
        write(*,*)'Top air thickness (ang): ', LenVac
        READ(IPF,*)LenNuc
        write(*,*)'Nucleus thickness (ang): ', LenNuc ! = radius when nucleation is Inhomogeneous
        READ(IPF,*)Homogeneous
        write(*,*)'Model homogeneous nucleation: ', Homogeneous        
      ELSE
        READ(IPF,*)Len1,Len2,Len3
      END IF
      
      ! If SiGe read alloy fraction
      IF((ALL(ListCrystal(:NCrystal).EQ.(/14, 32/))).OR.
     >   (ALL(ListCrystal(:NCrystal).EQ.(/14, 32, 5/))))THEN ! SiGex or SiGex:B (NEEDS GENERALIZATION IF THE ALLOY HAS DIFFERENT SPECIES)
        READ(IPF,*) alloyfraction
        write(*,*)'Ge fraction in substrate: ', alloyfraction
      ELSE
        alloyfraction = 0.0 ! Ge concentration in SiGe (set to zero for Si)
      END IF

      ! Probability for defect formation
      READ(IPF,*)PtransZig

      ! KMC Super-Lattice parameter (ang)
      READ(IPF,*)sf
      write(*,*)'KMC Super-Lattice parameter (ang): ', sf
      
      ! Exit strategy
      READ(IPF,*)ExitStrategy ! 'Iter' OR 'Time'
      write(*,*)'Exit strategy: ', ExitStrategy
      IF(ExitStrategy.EQ.'Iter')THEN
        READ(IPF,*)MaxIter
        READ(IPF,*)OutMolMol
        write(*,*)'MaxIter: ', MaxIter
        write(*,*)'OutMolMol: ', OutMolMol
      ELSE IF(ExitStrategy.EQ.'Time')THEN
        READ(IPF,*)MaxTime
        READ(IPF,*)OutTime
        write(*,*)'Total simulation time in seconds: ', MaxTime
        write(*,*)'Output frequency in seconds: ', OutTime
      END IF
!      READ(IPF,*)MaxZeta
      
      ! Define run type
      READ(IPF,*)RunType ! random R, test T or continuation C
      IF(RunType.EQ.'T') THEN
       write(*,*)'ATTENTION! You are running mulskips in test modality,
     > that is random numbers are not so random!!'
      ELSE IF(RunType.EQ.'C')THEN
        write(*,*)'ATTENTION!!! You have chosen a continuation run. 
     >LattCoo will be read from the following checkpoint 
     >file: ', cadfilename     
      END IF

      ! Random seed
      READ(IPF,*)IDUM
      WRITE(*,*)'IDUM',IDUM

      ! Do we want to save final state as checkpoint for a future continuation run?
      READ(IPF,*)SaveFinalState
      IF(SaveFinalState.EQ.'T') THEN
        READ(IPF,'(A)')restartfilename
        write(*,*)'ATTENTION! Storing checkpoint file ', restartfilename
        write(*,*)'This will occupy a lot of disk space...'
      ELSE 
        write(*,*)'SaveFinalState flag is not "T": checkpoint file will
     > not be written.'
      END IF 

      ! Do we want to save the coordinations at the final state? When InitSt='LA' it's redundant, do not read it 
      READ(IPF,*)SaveCoo
      IF(SaveCoo.EQ.'T') THEN
        READ(IPF,'(A)')coofilename
        write(*,*)'ATTENTION! Storing final Coor file ', coofilename
        write(*,*)'This might occupy a a lot of disk space...'
      ELSE 
        IF(InitSt.EQ.'LA')THEN
          write(*,*)'ERR: SaveCoo flag must be "T" for LA simulations!'
          STOP
        ELSE
          write(*,*)'SaveCoo flag is not "T": Coor file not stored'
        END IF
      END IF 

      ! Setup Probabilities for KMC events
      IF(InitSt.EQ.'ST')THEN
        CALL SetProbabilityFromField()
      ELSE
        CALL SetProbability()
      END IF
      CALL AllocateArrays()

      ! ------------------------------------------END INPUT

      ! Setup input structure
      IF(InitSt.EQ.'F')THEN
         Len1=Len1-MOD(Len1,24) ! Len1: init Slab height.
         CALL SetSiC3C(Len1)
      ELSE IF(InitSt.EQ.'Fw')THEN
         Len1=Len1-MOD(Len1,24)
         Len2=Len2-MOD(Len2,24) 
         CALL SetSiC3Cww(Len1,Len2) ! Len1: init Slab height. Len2: wall height. 
      ELSE IF(InitSt.EQ.'Fa')THEN
         Len1=Len1-MOD(Len1,24)
         Len2=Len2-MOD(Len2,24)
         Len3=Len3-MOD(Len3,24)
         CALL SetSiC3Cwwaperture(Len1,Len2,Len3) ! Len1: init Slab height. Len2: wall height. Len3: aperture size along x
      ELSE IF(InitSt.EQ.'FS')THEN
         Len1=Len1-MOD(Len1,24)
         CALL SetSi(Len1) ! Len1: init Slab height. 
      ELSE IF(InitSt.EQ.'SG')THEN
         Len1=Len1-MOD(Len1,24)
         CALL SetSiGex(Len1) ! Len1: init Slab height. 
      ELSE IF(InitSt.EQ.'A')THEN
         CALL SetAFBSymZimb1(Len1)
      ELSE IF(InitSt.EQ.'I')THEN
         CALL SetInvPOff(Len1,Len2)
      ELSE IF(InitSt.EQ.'D')THEN
         CALL SetInvPC(Len1,Len2)
      ELSE IF(InitSt.EQ.'Z')THEN
         CALL SetInvPSi(Len1,Len2)
      ELSE IF(InitSt.EQ.'J')THEN
         CALL SetInvPwAPB(Len1,Len2)
      ELSE IF(InitSt.EQ.'K')THEN
         CALL SetInvPwOnlySiorC(Len1,Len2)
      ELSE IF(InitSt.EQ.'L')THEN
         CALL SetInvPwOnlyCorSi(Len1,Len2)         
      ELSE IF(InitSt.EQ.'S')THEN
         CALL SetNCSp(Len1)
      ELSE IF(InitSt.EQ.'SS')THEN
         CALL SetNCSpSi(Len1)
      ELSE IF(InitSt.EQ.'C')THEN
         CALL SetNC(Len1,Len2,Len3)
      ELSE IF(InitSt.EQ.'T')THEN
         CALL SetTrench(Len1,Len2,Len3)
      ELSE IF(InitSt.EQ.'ST')THEN
         CALL SetTrenchSi(Len1,Len2,Len3,angle,Ttop,Tbottom) ! angle float in degrees
      ELSE IF(InitSt.EQ.'E')THEN
         SiteC(1)=Len2
         SiteC(2)=Len3
         SiteC(3)=Len1-Len2
         write(*,*)Len1,SiteC
         CALL SetSiCSF(Len1,SiteC)
      ELSE IF(InitSt.EQ.'FF')THEN ! silicon solid-liquid finfet shape 
         Len1=Len1-MOD(Len1,24)
         Len2=Len2-MOD(Len2,24) 
         Len3=Len3-MOD(Len3,24) 
         Len4=Len4-MOD(Len4,24) 
         CALL SetSiFINFET(Len1,Len2,Len3,Len4) ! Len1: init Slab height. Len2: wall bottom. Len3: wall top. Len4: aperture size along x.
      ELSE IF(InitSt.EQ.'IN')THEN
         CALL SetCAD()
      ELSE IF(InitSt.EQ.'SH')THEN
         ! read number of filled sites in shape and pass it to SetShape
         OPEN(50,FILE=cadfilename,STATUS='OLD')     
         READ(50,*)NFilled
         CLOSE(50)
         CALL SetShape(NFilled)
      ELSE IF(InitSt.EQ.'LA')THEN
         tmax=0
         CALL SetT()  ! do this first to allocate T, which will be used in Get_Prob_Ini in SetCAD
         CALL SetCAD()
         write(*,*)'Max temperature is: ', tmax, 'K'
      ELSE
         write(*,*)'initialization not_implemented'
         STOP
      END IF


      write(*,*)'Atoms,Voids,AdAtoms',NumAtoms,NumVoids,NumAdAtom
      CountCrystalOld = CountCrystal
      CALL WriteMolMolSource(Time, IterOut, framecounter)
      
      ! --------------------------- START KMC loop
      ! Run until MaxIter or MaxTime are achieved
      DO
       IterOut=Iter

       IF(NumAdAtom.LE.20)THEN
         write(*,*)"few MC particles: stop MC",Iter,NumAdAtom
         exit
       END IF

       CALL PickTreeEvent(Ind,Patom,RunType)
       Time=Time+1.0/Tree(1)
       Itrans=ListAdAtom(Ind) % Ind_Event
       
       SELECT CASE (ExitStrategy)
         CASE ('Iter') 
           IF(MOD(Iter,OutMolMol).EQ.0)THEN
              framecounter = IterOut/OutMolMol
              CALL WriteMolMolSource(Time, IterOut, framecounter)
              Site=ListAdAtom(Ind) % AtomXYZ
              write(*,*)'Iter, Time, Site',Iter,Time,Site
              write(*,*)'NatNvoidNad',NumAtoms,NumVoids,NumAdAtom
           END IF
         CASE ('Time') 
           IF(Time.GE.((framecounter+1)*OutTime))THEN
              CALL WriteMolMolSource(Time, IterOut, framecounter+1)
              Site=ListAdAtom(Ind) % AtomXYZ
              write(*,*)'Iter, Time, Site',Iter,Time,Site
              write(*,*)'NatNvoidNad',NumAtoms,NumVoids,NumAdAtom
              framecounter = framecounter+1
           END IF         
       END SELECT       

       SELECT CASE (Itrans)
         CASE (:3) ! deposition crystal species (e.g. 1=Si, 2=C); will need to shift all Itrans for other events by 1 if we want to include a third crystal species
C           IF(MOD(Iter,OutMolMol).EQ.0)THEN
C             write(*,*)'Iter',Iter,' Deposition',Ind,Patom,Itrans
C           END IF
           CALL Deposition(Ind,Itrans)
           CountCrystal(Itrans) = CountCrystal(Itrans) + 1
         CASE (4) ! evaporation of crystal species
C            IF(MOD(Iter,OutMolMol).EQ.0)THEN
C              write(*,*)'Iter',Iter,'Evaporation',Ind,Patom
C            END IF
           CALL Evaporation(Ind)
         CASE (5) ! desorption of coverage species
C            Site=ListAdAtom(Ind) % AtomXYZ
C            write(*,*)'Iter',Iter,Ind,Patom,'   desorption of',
C      >      IBITS(LattCoo(Site(1),Site(2),Site(3)),PosIndex,LenIndex)
           CALL desorption(Ind)
         CASE (6:) ! absorption of coverage species (Itrans >= 6)
           CovInd = ListCov(Itrans-5)
C            write(*,*)'Iter',Iter,Ind,Patom,'absorption of',CovInd
           CALL absorption(Ind,CovInd)
           CountCov(Itrans-5) = CountCov(Itrans-5) + 1
       END SELECT

       Iter = Iter +1

       SELECT CASE (ExitStrategy)
         CASE ('Iter') 
           IF (Iter.GT.MaxIter) exit
         CASE ('Time') 
           IF (Time.GT.MaxTime) exit
C          CASE ('Zeta') 
C            IF (Site(3).GT.MaxZeta) exit
       END SELECT


      END DO
      ! --------------------------- END KMC loop

      SELECT CASE (ExitStrategy)
        CASE ('Iter') 
          CALL WriteMolMolSource(Time, Iterout, Iterout/OutMolMol)
        CASE ('Time') 
          CALL WriteMolMolSource(Time, Iterout, framecounter)
      END SELECT      

      CALL Fileclose()


      ! Writing checkpoint file for a restart run
      IF(SaveFinalState.EQ.'T') THEN
        write(*,*)'Writing checkpoint file for a restart run', 
     >         restartfilename
        OPF17=37
        OPEN(OPF17,FILE=restartfilename,STATUS='REPLACE')
  !16    FORMAT(500000I4) ! ATT ATT ATT format must be > lenX*Leny!!!
        DO k=0,LenZ-1
          LenSaveOut=0
          DO j=0,LenY-1
            DO i=0,LenX-1
              LenSaveOut(1+i+LenX*j) = LattCoo(i,j,k)
            END DO
          END DO
          write(OPF17,*)LenSaveOut
        END DO
        CLOSE(OPF17)
        write(*,*)'Done writing checkpoint file for a restart run', 
     >         restartfilename      
      END IF


      ! Store final Coor. If LA mode is selected, then this is mandatory
      IF(SaveCoo.EQ.'T')THEN
        write(*,*)'Writing Coor (LA phases) file ', coofilename
        OPF16=36
        OPEN(OPF16,FILE=coofilename,STATUS='REPLACE')
  16    FORMAT(500000I4) ! ATT ATT ATT format must be > lenX*Leny!!!
        DO k=0,LenZ-1
          LenOut=0
          DO j=0,LenY-1
            DO i=0,LenX-1
              IQstat = IBITS(LattCoo(i,j,k),PosIndex,LenIndex)
              Coo    = IBITS(LattCoo(i,j,k),PosCoor,LenCoor)
              Occ    = IBITS(LattCoo(i,j,k),PosOcc,LenOcc)
              IF(IQStat.NE.113)THEN
                IF(Occ.EQ.1) LenOut(1+i+LenX*j)=Coo
              ELSE
                LenOut(1+i+LenX*j)=10 ! Write 10 if site = wall
              END IF
            END DO
          END DO
          write(OPF16,16)LenOut
        END DO
        CLOSE(OPF16)
        write(*,*)'Done writing Coor (LA phases) file ', coofilename
      END IF


      END PROGRAM KLMCSiC3C
