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
      MODULE DefSystem
       USE DefDerType
       IMPLICIT  NONE
!     LenX,LenY,LenZ:    Dimensioni del sistema
!     LattCoor      :    Numero Di Coordinazione Relativo al Sito
!     LattInd       :    Indice del Sito nel vettore Listporb
       INTEGER, PARAMETER :: LenX=1920, LenY=120, LenZ=19200
       INTEGER, PARAMETER :: LenXGr=4*LenX/3,LenYGr=2*LenY/3
       INTEGER :: Nedge=-1
       INTEGER, PARAMETER :: NSites=LenX*LenY*LenZ
       INTEGER, PARAMETER :: NeigFCC=8,NeigGRA=6,NextNeig=18
       INTEGER, PARAMETER :: NeigBCC=8
       INTEGER, PARAMETER :: NTransD=6,NTransE1=2,NTransE2=3
       INTEGER, PARAMETER :: NumMcMax=500000,NumAtMax=100000000
       INTEGER, PARAMETER :: NumVdMax=100000,NumBfMax=500000
       INTEGER  NumAtoms,NumVoids,NumAdAtom

       INTEGER, DIMENSION(0:LenX-1,0:LenY-1,0:LenZ) :: LattCoo
       INTEGER, DIMENSION(0:LenX-1,0:LenY-1,0:LenZ) :: LattInd

       INTEGER, DIMENSION(3), PARAMETER :: Box=(/LenX,LenY,LenZ/)
       INTEGER, DIMENSION(3), PARAMETER :: Box2D=(/LenX,LenY,0/)

       TYPE(MCParticle), DIMENSION(NumMcMax) :: ListAdAtom
       TYPE(BulkParticle), DIMENSION(NumAtMax) :: ListAtom
       TYPE(Vacancy), DIMENSION(NumVdMax) :: ListVoid
       TYPE(AtBuffer), DIMENSION(NumBfMax) :: ListAtBuffer
       REAL(8), DIMENSION(1:NTransD) :: PTransD
       REAL(8), DIMENSION(NTransE1,0:NTransE2,0:NtransE2) :: PTransE
       REAL(8) :: PtransZig


       REAL(8) :: EBind=0.8,EBind1=1.,EbindG=2.5,Energ,Enemed
       REAL(8) :: EevapC=2.5,EevapSi=0.5


       INTEGER,PARAMETER:: PosIndex=0,LenIndex=22,PosOcc=23,LenOcc=1
       INTEGER,PARAMETER:: PosCoor=24,LenCoor=3
       INTEGER,PARAMETER:: PosSi=28,PosC=29,LenSiC=2
       INTEGER,PARAMETER:: PosSiC=28,LenSi=1,LenC=1
       INTEGER, DIMENSION(8)  :: NComplBCC
       INTEGER, DIMENSION(8)  :: NComplFCC
       INTEGER, DIMENSION(3,8)  :: NComplDia
       INTEGER, DIMENSION(12) :: NBit
       INTEGER  JumpStat(NeigFCC)
       INTEGER, DIMENSION(3, NeigBCC) :: JumpDia       ! Parameters
       INTEGER, DIMENSION(3, NeigBCC) :: JumpBCC       ! Parameters
       INTEGER, DIMENSION(3, NeigFCC) :: JumpFCC       ! Parameters
       INTEGER, DIMENSION(3, NextNeig) :: Jump2NN       ! Parameters
       INTEGER, DIMENSION(3, NeigGRA) :: JumpGRA       ! Parameters
       INTEGER, DIMENSION(3,3,8) :: IMatRotat

       DATA JumpDia /3, 3, 3,   3,-3,-3, -3,-3, 3, -3, 3,-3,
     >               3, 3,-3,   3,-3, 3, -3,-3,-3, -3, 3, 3 /
       DATA JumpBCC /1, 0, 0,  -1, 0, 0,  0, 1, 0,  0,-1, 0,
     >               0, 0, 1,   0, 0,-1,  1, 1, 1, -1,-1,-1 /
       DATA JumpFCC /1, 0, 0,   0, 1, 0, -1, 1, 0, -1, 0, 0,
     >               0,-1, 0,   1,-1, 0,  0, 0, 1,  0, 0,-1 /
       DATA JumpGRA /1, 0, 0,  -1, 0, 0,  1,-1, 0, ! pari
     >               1, 0, 0,  -1, 0, 0, -1, 1, 0/ ! dispari
       DATA Jump2NN /1, 1, 0,  -1, 2, 0, -2, 1, 0, -1,-1, 0,
     >               1,-2, 0,   2,-1, 0,
     >               1, 0, 1,   0, 1, 1, -1, 1, 1, -1, 0, 1,
     >               0,-1, 1,   1,-1, 1,
     >               1, 0,-1,   0, 1,-1, -1, 1,-1, -1, 0,-1,
     >               0,-1,-1,   1,-1,-1/
       DATA NComplBCC /2, 1, 4, 3, 6, 5, 8, 7 /
       DATA NComplDia /2, 3, 4,
     >                 3, 4, 1,
     >                 4, 1, 2,
     >                 1, 2, 3,
     >                 6, 7, 8,
     >                 7, 8, 5,
     >                 8, 5, 6,
     >                 5, 6, 7  /

!       DATA NComplFCC /2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11/
       DATA NComplFCC /4, 5, 6, 1, 2, 3, 8, 7/
       DATA NBit /1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048/
       DATA IMatRotat / 2, 2,-1,  -1, 2, 2,  2,-1, 2,
     >                  2,-2, 1,   1, 2, 2, -2,-1, 2,
     >                  2, 2, 1,  -1, 2,-2, -2, 1, 2,
     >                  2,-2,-1,   1, 2,-2,  2, 1, 2,
     >                  2,-1,-2,   2, 2, 1,  1,-2, 2,
     >                  2, 1, 2,  -2, 2, 1, -1,-2, 2,
     >                  2,-1, 2,   2, 2,-1, -1, 2, 2,
     >                  2, 1,-2,  -2, 2,-1,  1, 2, 2 /

!!!!!!!!!!!!!!!
!     Tree for the event picking
      INTEGER :: Levels, SizeTree
      REAL(8), DIMENSION(:), ALLOCATABLE :: Tree
!!!!!!!!!!!!!!!


      END MODULE DefSystem


