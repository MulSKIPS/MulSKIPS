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
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! for usage with sockets
      use f90sockets, only : open_socket, writebuffer, readbuffer
      use, intrinsic :: iso_c_binding
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      INTEGER :: i,j,k,l,in,irot,K1,Itrans
      INTEGER, PARAMETER :: k10 = selected_int_kind(10)
      INTEGER(kind=k10) :: Iter=1, IterOut=0
      INTEGER(kind=k10) :: MaxIter=2100000000
      INTEGER(kind=k10) :: OutMolMol=100000000
      INTEGER :: framecounter=0
      INTEGER :: LAcounter=0
      REAL(8) :: MaxTime=100000000.0
      REAL(8) :: OutTime=10000000.0
      INTEGER :: Ind,IndTest
      INTEGER :: MGigaCycle=400,IGigaCycle=0,Kilo=1000,Mega=1000000
      INTEGER :: MaxZeta=10000
      INTEGER :: Len=LenZ-1,Len1=120,Len2=120,Len3=120,Len4=120,Len5=120
      INTEGER :: OPF16,OPF17,IQstat,Coo,Occ
      INTEGER :: LenOut(LenX*LenY)
      INTEGER :: LenSaveOut(LenX*LenY)
      INTEGER :: OPF46, IAt, zz, deltaz, ctot, cGe, cSi
      REAL(8) :: angle ! degrees
      REAL(8) :: Ttop, Tbottom ! Kelvin

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


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! socket variables
      integer, parameter :: msglen=12   ! length of the headers of the driver/wrapper communication protocol
      integer socket, inet, port        ! socket id & address of the server
      character(len=1024) :: host

      ! command line parsing
      character(len=1024) :: cmdbuffer
      logical :: hflag=.false., pflag=.false.

      ! socket communication buffers
      character(len=12) :: header
      logical :: isinit=.false.   ! The driver has been initialised by the server
      logical :: phasesready=.false.   ! The driver has finished computing and can send data to server
      logical :: XGeready=.false.   ! The driver has finished computing and can send data to server
      logical :: geoisinit=.false.   ! The driver has finished computing and can send data to server
      logical :: fieldisinit=.false.   ! The driver has finished computing and can send data to server

C       real(kind=4) , allocatable :: msgbuffer(:)
      INTEGER, allocatable :: msgbuffer(:) ! if you change type of "a" you need to change this type as well

      ! data to send and receive (REMEMBER TO CHANGE ALSO type(msgbuffer)!!!)
      ! IS THERE A WAY TO AVOID THESE ALLOCATIONS BELOW?
      INTEGER, DIMENSION(LenZ,LenX*LenY) :: a=0 ! geo 
      INTEGER, DIMENSION(LenZ,LenX*LenY) :: b=0 ! field 
      INTEGER, DIMENSION(LenZ,LenX*LenY) :: c=0 ! phases 
      INTEGER, DIMENSION(LenZ,LenX*LenY) :: d=0 ! Ge fraction 

      integer :: ii
!       integer, allocatable :: sh(:)

      ! intialize defaults
      inet = 1
      host = "localhost"//achar(0)
      port = 31475

      ! read command arguments
      if (mod(command_argument_count(), 2) /= 0) then
          call helpmessage
          stop "ended"
      end if
      
      do ii = 1, command_argument_count()
          call get_command_argument(ii, cmdbuffer)
          if (cmdbuffer == "-h") then ! read the hostname
              hflag = .true.
          elseif (cmdbuffer == "-p") then ! reads the port number
              pflag = .true.
          elseif (hflag .and. mod(ii, 2) == 0) then
              host = trim(cmdbuffer)//achar(0)
              hflag = .false.
          elseif (pflag .and. mod(ii, 2) == 0) then
              read(cmdbuffer,*) port
              pflag = .false.
          else
              write(*,*) " unrecognized command line argument", ii
              call helpmessage
              stop "ended"
          endif
       enddo

      ! open port
      write(*,*) " driver - connecting to host ", trim(host)
      write(*,*) " on port ", port, " using an internet socket."
      call open_socket(socket, inet, port, host)


      ! main loop
C       a = reshape([1,2,3,4,5,6,7,8], shape(a)) ! we define a as LattCoo inside defsystem!!!
!       dims = size(shape(a))
!       sh = shape(a)

      ! remember that for testing purposes LenZ is set to 30 in the python!!!!!
      write(*,*) "Initial array is:"   
      write(*,*) a(1,1), shape(a), size(shape(a))
      write(*,*) b(1,1), shape(b), size(shape(b))
      write(*,*) c(1,1), shape(c), size(shape(c))
      write(*,*) d(1,1), shape(d), size(shape(d))
C       do ii=1,sh(1)
C           write(*,*) a(ii,:)
C       end do

      


      do while (.true.) ! loops forever (or until the wrapper ends!)

          write(*,*) 'Current status: geoisinit, fieldisinit, phasesready, XGeready', geoisinit, fieldisinit, phasesready, XGeready

          ! reads from the socket one message header
          ! I will Stop here till python sends me a message 
          call readbuffer(socket, header, msglen)

          write(*,*) "@ message from server: ", trim(header)



          if (trim(header) == "STATUS") then ! the wrapper is inquiring on what we are doing
              if (.not. isinit) then
                  call writebuffer(socket, "NEEDINIT    ", msglen)  ! signals that we need initialization
                  write(*,*) "@ message to server: NEEDINIT"
              elseif (phasesready.OR.XGeready) then
                  call writebuffer(socket, "HAVEDATA    ", msglen)  ! signals that we are done computing and can data
                  write(*,*) "@ message to server: HAVEDATA"
              else
                  call writebuffer(socket, "READY       ", msglen)  ! we are idling and eager to compute something
                  write(*,*) "@ message to server: READY"
              endif


          elseif (trim(header) == "INIT") then     ! the driver is kindly sending a string for initialization
              write(*,*) " Initializing system from server"
              isinit=.true. ! we actually do nothing with this string, thanks anyway. could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver


          elseif (trim(header) == "SENDGEO") then  ! Server wants to send data to the driver
              if (.not. isinit) then
                  write(*,*) "Driver not iniliasied."
              elseif (phasesready.OR.XGeready) then
                  write(*,*) "Driver has data to send back to server"
              else ! Driver is ready to receive data
                  if (.not. allocated(msgbuffer)) allocate(msgbuffer(size(a)))
                  call readbuffer(socket, msgbuffer, size(a))

                  a = reshape(msgbuffer, shape(a))

                  write(*,*) "Received array from server:"
                  write(*,*) a(1,1)
                  write(*,*) b(1,1)
                  write(*,*) c(1,1)
                  write(*,*) d(1,1)
                  
                  geoisinit = .true.

              end if

          elseif (trim(header) == "SENDFIELD") then  ! Server wants to send data to the driver
              if (.not. isinit) then
                  write(*,*) "Driver not iniliasied."
              elseif (phasesready.OR.XGeready) then
                  write(*,*) "Driver has data to send back to server"
              else ! Driver is ready to receive data
                  if (.not. allocated(msgbuffer)) allocate(msgbuffer(size(b)))
                  call readbuffer(socket, msgbuffer, size(b))

                  b = reshape(msgbuffer, shape(b))

                  write(*,*) "Received array from server:"
                  write(*,*) a(1,1)
                  write(*,*) b(1,1)
                  write(*,*) c(1,1)
                  write(*,*) d(1,1)

                  fieldisinit = .true.

              end if

          elseif (trim(header) == "GETGEO") then  ! Server signaling driver to send data
              if (.not. isinit) then
                  write(*,*) "Driver not iniliasied."
              elseif (.not. geoisinit) then
                  write(*,*) "Driver does not have data to send"
              else
                  call writebuffer(socket, "DATAREADY   ", msglen)
                  write(*,*) "@ message to server: DATAREADY"

                  msgbuffer = reshape(a, [size(a)])   ! flatten data
                  call writebuffer(socket, msgbuffer, size(a)) ! writing data

                  write(*,*) "Sent array to server:"
                  write(*,*) a(1,1)
                  write(*,*) b(1,1)
                  write(*,*) c(1,1)
                  write(*,*) d(1,1)

              end if

          elseif (trim(header) == "GETFIELD") then  ! Server signaling driver to send data
              if (.not. isinit) then
                  write(*,*) "Driver not iniliasied."
              elseif (.not. fieldisinit) then
                  write(*,*) "Driver does not have data to send"
              else
                  call writebuffer(socket, "DATAREADY   ", msglen)
                  write(*,*) "@ message to server: DATAREADY"

                  msgbuffer = reshape(b, [size(b)])   ! flatten data
                  call writebuffer(socket, msgbuffer, size(b)) ! writing data

                  write(*,*) "Sent array to server:"
                  write(*,*) a(1,1)
                  write(*,*) b(1,1)
                  write(*,*) c(1,1)
                  write(*,*) d(1,1)

              end if


          elseif (trim(header) == "GETPHASES") then  ! Server signaling driver to send data
              if (.not. isinit) then
                  write(*,*) "Driver not iniliasied."
              elseif (.not. phasesready) then
                  write(*,*) "Driver does not have phases to send"
              else
                  call writebuffer(socket, "DATAREADY   ", msglen)
                  write(*,*) "@ message to server: DATAREADY"

C                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                   OPF16=36
C                   OPEN(OPF16,FILE='LattStat.dat',STATUS='REPLACE')
C 26                FORMAT(500000I4) ! ATT ATT ATT format must be > lenX*Leny!!!
C                   DO k=0,LenZ-1
C                     write(OPF16,26) c(1+k,:)
C                   END DO
C                   CLOSE(OPF16)
C                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  msgbuffer = reshape(c, [size(c)])   ! flatten data
                  call writebuffer(socket, msgbuffer, size(c)) ! writing data

                  write(*,*) "Sent array to server:"
                  write(*,*) a(1,1)
                  write(*,*) b(1,1)
                  write(*,*) c(1,1)
                  write(*,*) d(1,1)

                  phasesready = .false.

              end if
                  

          elseif (trim(header) == "GETXGE") then  ! Server signaling driver to send data
              if (.not. isinit) then
                  write(*,*) "Driver not iniliasied."
              elseif (.not. XGeready) then
                  write(*,*) "Driver does not have Ge fraction to send"
              else
                  call writebuffer(socket, "DATAREADY   ", msglen)
                  write(*,*) "@ message to server: DATAREADY"

C                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                   OPF16=36
C                   OPEN(OPF16,FILE='LattStat.dat',STATUS='REPLACE')
C 26                FORMAT(500000I4) ! ATT ATT ATT format must be > lenX*Leny!!!
C                   DO k=0,LenZ-1
C                     write(OPF16,26) c(1+k,:)
C                   END DO
C                   CLOSE(OPF16)
C                   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                  msgbuffer = reshape(d, [size(d)])   ! flatten data
                  call writebuffer(socket, msgbuffer, size(d)) ! writing data

                  write(*,*) "Sent array to server:"
                  write(*,*) a(1,1)
                  write(*,*) b(1,1)
                  write(*,*) c(1,1)
                  write(*,*) d(1,1)

                  XGeready = .false.

              end if


          else
              write(*,*) " unexpected header ", header
              stop "ended"
          endif

        ! CYCLE
        IF ((.NOT.geoisinit).OR.(.NOT.fieldisinit)) CYCLE   ! cycle until geo and temperature have been set for the first time
        write(*,*) 'KMC ACTUALLY STARTING NOW: geoisinit, fieldisinit, phasesready, XGeready', geoisinit, fieldisinit, phasesready, XGeready



        ! ------------------------------------------- INPUT (only for LAcounter=0)
        IF(LAcounter.EQ.0)THEN

          LattCoo=0
          LattInd=0
          counter=0
          NumAtoms=0
          NumVoids=0
          NumAdAtom=0
          NumAdAtomOcc=0
          NumOcc0=0
          NumOcc1=0
          NumOcc=0
          CountCrystal=0
          CountCov=0
         
          write(*,*)'KMC box size: ',LenX, LenY, LenZ
     
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
     >    ListCrystal(:NCrystal)
          IF(NCov.GT.0)THEN
            write(*,*)'Coverage species (NCov=',NCov,') --> Z:',
     >      ListCov(:NCov)
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
          ELSE IF(InitSt.EQ.'IN')THEN
            READ(IPF,'(A)')cadfilename
          ELSE IF(InitSt.EQ.'LA')THEN
            write(*,*)'****** MulSKIPS run for laser annealing ******'
            IF(NCov.GT.0)THEN
              write(*,*)'ERR NOT IMPLEMENTED: NCov=0 for LA'
              STOP
            END IF
            READ(IPF,'(A)')cadfilename
            READ(IPF,'(A)')tempfilename
            READ(IPF,*)Tm(:NCrystal)  ! Kelvin
            READ(IPF,*)P0(:NCrystal)
            READ(IPF,*)sigma(:NCrystal)
            READ(IPF,*)Tflex(:NCrystal) ! Kelvin
            write(*,*)'Melting temperature in Kelvin: ', Tm(:NCrystal)
            READ(IPF,*)LenVac
            write(*,*)'Top air thickness (ang): ', LenVac
            READ(IPF,*)LenNuc
            write(*,*)'Nucleus thickness (ang): ', LenNuc ! = radius when nucleation is Inhomogeneous
            READ(IPF,*)Homogeneous
            write(*,*)'Model homogeneous nucleation: ', Homogeneous
            IF(Homogeneous.EQ.'D')THEN ! set stacking fault with these coordinates
             READ(IPF,*)Len1,Len2,Len3,Len4,Len5
            END IF
          ELSE
            READ(IPF,*)Len1,Len2,Len3
          END IF
          
          ! If SiGe read alloy fraction
          IF(ALL(ListCrystal(:NCrystal).EQ.(/14, 32/)))THEN ! SiGex
            READ(IPF,*) alloyfraction
            write(*,*)'Ge fraction in substrate: ', alloyfraction        
            READ(IPF,*)LenSiGe
            write(*,*)'SiGe thickness (ang): ', LenSiGe ! = radius when nucleation is Inhomogeneous
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
!          READ(IPF,*)MaxZeta
          
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
            IF(coofilename.EQ.'None')THEN
              write(*,*)'Phases passed through sockets...'
            ELSE
              write(*,*)'ATTENTION! Storing final Coor file ', coofilename
              write(*,*)'This might occupy a a lot of disk space...'
            END IF
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
             CALL SetSiGex(Len1, alloyfraction) ! Len1: init Slab height. 
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
          ELSE IF(InitSt.EQ.'LA')THEN
             tmax=0
             CALL SetT(b)  ! do this first to allocate T, which will be used in Get_Prob_Ini in SetCAD
             IF(Homogeneous.EQ.'D')THEN ! initiate stacking fault
              ! Len1 : surface height
              ! Len2,3,4 : coordinates of SF vertex
              ! Len5 : offset between top line of SF and the surface (may also be 0)
              SiteC(1)=Len2
              SiteC(2)=Len3
              SiteC(3)=Len4
              write(*,*)'Len1,SiteC,Len5',Len1,SiteC,Len5
              CALL SetSiGeSF(Len1,SiteC,Len5,a,alloyfraction)
             ELSE
              CALL SetCAD(a, alloyfraction)
             END IF
             write(*,*)'Max temperature is: ', tmax, 'K'
          ELSE
             write(*,*)'initialization not_implemented'
             STOP
          END IF

          ! Store initial number of occupied sites
          NumOcc1 = NumAtoms + NumAdAtomOcc
          !write(*,*)'CountCrystal',CountCrystal
          DO k=1,NCrystal
            IF(k.EQ.1) THEN
              InitCountCrystal(k) = CountCrystal(k) + INT((NumOcc0 - NumOcc1)*(1-alloyfraction))
            ELSE IF(k.EQ.2) THEN
              InitCountCrystal(k) = CountCrystal(k) + INT((NumOcc0 - NumOcc1)*alloyfraction)
            END IF
          END DO
          write(*,*)' '
          write(*,*)'Initial number of occupied sites (before nucleation): ', NumOcc0 ! Counted inside SetCAD.f
          write(*,*)'Initial CountCrystal (before nucleation): ', InitCountCrystal ! Counted inside SetCAD.f
          write(*,*)'Initial number of occupied sites (after nucleation): ', NumOcc1
          write(*,*)'Initial CountCrystal (after nucleation): ', CountCrystal ! Counted inside SetCAD.f
          write(*,*)'Initial number of liquid sites (after nucleation): ', NumOcc0 - NumOcc1
          write(*,*)'Initial CountLiquid (after nucleation): ', CountLiquid ! Counted inside SetCAD.f
          write(*,*)' '

        ELSE

          CALL SetT(b)  ! Update T from server
          ! maybe it would be safer to reinitialize b=0 here
        END IF
        ! ------------------------------------------END INPUT (only for LAcounter=0)

        Time=0.
        Iter=1
        IterOut=0
        framecounter=0

        write(*,*)'Atoms,Voids,AdAtoms,AdAtomsOcc',NumAtoms,NumVoids,NumAdAtom,NumAdAtomOcc
        CountCrystalOld = CountCrystal
        CALL WriteMolMolSource(Time, IterOut, framecounter)
        
        ! --------------------------- START KMC loop
        ! Run until MaxIter or MaxTime are achieved
        DO
         IterOut=Iter

         CountLiquid = InitCountCrystal - CountCrystal

         IF(NumAdAtom.LE.20)THEN
           write(*,*)"few MC particles: stop MC",Iter,NumAdAtom,NumAdAtomOcc
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
                write(*,*)'NatNvoidNadNadocc',NumAtoms,NumVoids,NumAdAtom,NumAdAtomOcc
                NumOcc = NumAtoms + NumAdAtomOcc
                write(*,*)'Number of occupied sites: ', NumOcc
             END IF
           CASE ('Time') 
             IF(Time.GE.((framecounter+1)*OutTime))THEN
                CALL WriteMolMolSource(Time, IterOut, framecounter+1)
                Site=ListAdAtom(Ind) % AtomXYZ
                write(*,*)'Iter, Time, Site',Iter,Time,Site
                write(*,*)'NatNvoidNadNadocc',NumAtoms,NumVoids,NumAdAtom,NumAdAtomOcc
                NumOcc = NumAtoms + NumAdAtomOcc
                write(*,*)'Number of occupied sites: ', NumOcc
                framecounter = framecounter+1
             END IF         
         END SELECT       

         SELECT CASE (Itrans)
           CASE (:2) ! deposition crystal species (e.g. 1=Si, 2=C); will need to shift all Itrans for other events by 1 if we want to include a third crystal species
C             IF(MOD(Iter,OutMolMol).EQ.0)THEN
C               write(*,*)'Iter',Iter,' Deposition',Ind,Patom,Itrans
C             END IF
             IF((InitSt.EQ.'LA').AND.
     >                (ALL(ListCrystal(:NCrystal).EQ.(/14, 32/)))) THEN
               IF(CountCrystal(Itrans).GE.InitCountCrystal(Itrans))THEN 
                 Itrans = MOD(Itrans,2)+1
               END IF
             END IF
             CALL Deposition(Ind,Itrans)
             CountCrystal(Itrans) = CountCrystal(Itrans) + 1
           CASE (3) ! evaporation of crystal species
C              IF(MOD(Iter,OutMolMol).EQ.0)THEN
C                write(*,*)'Iter',Iter,'Evaporation',Ind,Patom
C              END IF
             CALL Evaporation(Ind)
           CASE (4) ! desorption of coverage species
C              Site=ListAdAtom(Ind) % AtomXYZ
C              write(*,*)'Iter',Iter,Ind,Patom,'   desorption of',
C        >      IBITS(LattCoo(Site(1),Site(2),Site(3)),PosIndex,LenIndex)
             CALL desorption(Ind)
           CASE (5:) ! absorption of coverage species (Itrans >= 5)
             CovInd = ListCov(Itrans-4)
C              write(*,*)'Iter',Iter,Ind,Patom,'absorption of',CovInd
             CALL absorption(Ind,CovInd)
             CountCov(Itrans-4) = CountCov(Itrans-4) + 1
         END SELECT

         Iter = Iter +1

         SELECT CASE (ExitStrategy)
           CASE ('Iter') 
             IF (Iter.GT.MaxIter) exit
           CASE ('Time') 
             IF (Time.GT.MaxTime) exit
C            CASE ('Zeta') 
C              IF (Site(3).GT.MaxZeta) exit
         END SELECT


         ! Check if initial number of occupied atoms was recovered 
         NumOcc = NumAtoms + NumAdAtomOcc
         IF (NumOcc.GT.NumOcc0)THEN
           write(*,*)' '
           write(*,*)'Current number of occupied sites: ', NumOcc
           write(*,*)'   > initial number of occupied sites: ', NumOcc0
           write(*,*)'... STOPPING MulSKIPS ...'
           write(*,*)' '
           exit  ! pass phases and stop mulskips
           ! remember that here we make an error because Delta t is not over!
         END IF


        END DO
        ! --------------------------- END KMC loop

        !!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! If SiGe read alloy fraction
        IF(ALL(ListCrystal(:NCrystal).EQ.(/14, 32/)))THEN ! SiGex
          deltaz = 24-MOD(LenZ,24)
          OPF46=46
          OPEN(OPF46,FILE='zint_xGeS.dat',STATUS='REPLACE')
          DO zz=0,LenZ-1, deltaz
            cGe = 0
            cSi = 0
            DO k=zz,zz+deltaz-1
              DO j=0,LenY-1
                DO i=0,LenX-1
                  Occ    = IBITS(LattCoo(i,j,k),PosOcc,LenOcc)
                  IF(Occ.EQ.1)THEN
                    IAt    = IBITS(LattCoo(i,j,k),PosSiC,LenSiC)
                    IF(IAt.EQ.2)THEN
                      cGe = cGe + 1
                    ELSE IF(IAt.EQ.1)THEN
                      cSi = cSi + 1
                    END IF
                  END IF
                END DO
              END DO
            END DO
            IF((cGe.EQ.0).AND.(cSi.EQ.0))THEN
              write(OPF46,*) k, 0, 0, 0, 0
            ELSE
              write(OPF46,*) k, cGe, cSi, cGe+cSi, FLOAT(cGe)/(cGe+cSi)
            END IF
          END DO
          CLOSE(OPF46)
          write(*,*)'Done writing Ge conc. vs. z '        
        END IF
        !!!!!!!!!!!!!!!!!!!!!!!!!!!

        SELECT CASE (ExitStrategy)
          CASE ('Iter') 
            CALL WriteMolMolSource(Time, Iterout, Iterout/OutMolMol)
          CASE ('Time') 
            CALL WriteMolMolSource(Time, Iterout, framecounter)
        END SELECT      

        ! Writing checkpoint file for a restart run
        IF(SaveFinalState.EQ.'T') THEN
          write(*,*)'Writing checkpoint file for a restart run', 
     >         restartfilename
          OPF17=37
          OPEN(OPF17,FILE=restartfilename,STATUS='REPLACE')
!16        FORMAT(500000I4) ! ATT ATT ATT format must be > lenX*Leny!!!
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
          IF(coofilename.EQ.'None')THEN
            write(*,*)'Sending phases and Ge fraction back to LA code '
            c=0  ! this is important, otherwise the phases passed to python will be wrong! (atoms melt during the current iter will still be considered solids! Because of IF Occ.EQ.1. Solids would be then overestimated, and even remain constant during the LA simulation!)
            d=0  ! this is important, otherwise the phases passed to python will be wrong!
            DO k=0,LenZ-1
              DO j=0,LenY-1
                DO i=0,LenX-1
                  IQstat = IBITS(LattCoo(i,j,k),PosIndex,LenIndex)
                  Coo    = IBITS(LattCoo(i,j,k),PosCoor,LenCoor)
                  IAt    = IBITS(LattCoo(i,j,k),PosSiC,LenSiC)
                  Occ    = IBITS(LattCoo(i,j,k),PosOcc,LenOcc)
                  IF(IQStat.NE.113)THEN
                    IF(Occ.EQ.1)THEN
                      c(1+k,1+i+LenX*j)=Coo
                      d(1+k,1+i+LenX*j)=IAt
                    ENDIF
                  ELSE
                    c(1+k,1+i+LenX*j)=10 ! Write 10 if site = wall
                    d(1+k,1+i+LenX*j)=10 ! Write 10 if site = wall
                  END IF
                END DO
              END DO
            END DO            
            write(*,*)'Done sending phases and Ge fraction back to LA code '

          ELSE

            write(*,*)'Writing Coor (LA phases) file ', coofilename
            OPF16=36
            OPEN(OPF16,FILE=coofilename,STATUS='REPLACE')
16          FORMAT(500000I4) ! ATT ATT ATT format must be > lenX*Leny!!!
            DO k=0,LenZ-1
              LenOut=0  ! this is important, otherwise the phases passed to python will be wrong! (Solids would be overestimated, remaining constant after the first few iter)
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
        END IF
        
        fieldisinit=.false.  ! now wait again until python passes the new temperature map... 
        phasesready = .true.   
        XGeready = .true.   
        LAcounter = LAcounter +1 

      enddo

      CALL Fileclose()



      CONTAINS

C         SUBROUTINE do_something(a)
C             real(kind=4), intent(inout) :: a(:,:)
C             call sleep(5)
C             a = a+2
C         end SUBROUTINE do_something

        SUBROUTINE helpmessage ! Help banner
            write(*,*) " syntax: driver.x -h hostname -p port "
        end SUBROUTINE helpmessage

      END PROGRAM KLMCSiC3C
