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
      SUBROUTINE WriteMolMolSource(Time,Iter,Iout)

      USE  Definitions
      USE  DefSystem
      IMPLICIT    NONE
      INTEGER, PARAMETER :: k10 = selected_int_kind(10)
      INTEGER(kind=k10) :: Iter
      INTEGER :: Iout
      REAL(8)     Time
      REAL(8), DIMENSION(NCrystalMax) :: diff

      INTEGER FN
      CHARACTER(Len=13) :: XYZFileName,SourceFileName,GifFileName
      CHARACTER(Len=15) :: XYZ_dFileName,XYZ_vFileName,XYZ_wFileName
      CHARACTER(Len=9)  :: FileNameBase
      CHARACTER(Len=4)  :: XYZExt='.xyz',SourceExt='.src'
      CHARACTER(Len=6)  :: DefectExt='_d.xyz',VacancyExt='_v.xyz'
      CHARACTER(Len=6)  :: WrongExt='_w.xyz'
      CHARACTER(Len=4)  :: GifExt='.gif'

      IF(Iter.NE.0)THEN
        write(*,*)'SUM(CountCrystal)', SUM(CountCrystal)
        write(*,*)'CountCrystal', CountCrystal
        write(*,*)'CountCrystalOld', CountCrystalOld
        diff=FLOAT(CountCrystal-CountCrystalOld)
        write(*,*)'diff', diff
        write(*,*)'xGesolid', diff/
     >                     (SUM(CountCrystal)-SUM(CountCrystalOld))
        ! Update Old counts
        CountCrystalOld = CountCrystal
        ! Write out coverage counts 
        IF(NCov.NE.0)THEN
          write(*,*)'CountCov', CountCov
        END IF
      END IF

      FN = 97
      CALL GetOutputFileName(Time, Iout, FileNameBase)
      write(*,*)FileNameBase
      XYZFileName = FileNameBase//XYZExt
      XYZ_dFileName = FileNameBase//DefectExt
      XYZ_vFileName = FileNameBase//VacancyExt
      XYZ_wFileName = FileNameBase//WrongExt
      SourceFileName = FileNameBase//SourceExt
      GifFileName = FileNameBase//GifExt
      IF((NCrystal.EQ.2).AND.(ANY(ListCrystal.EQ.6)))THEN
        OPEN(FN+6,FILE=XYZ_wFileName(:LEN_TRIM(XYZ_wFileName)))
      END IF
      OPEN(FN+5,FILE=XYZ_vFileName(:LEN_TRIM(XYZ_vFileName)))
      OPEN(FN+4,FILE=XYZ_dFileName(:LEN_TRIM(XYZ_dFileName)))
      OPEN(FN,FILE=SourceFileName(:LEN_TRIM(SourceFileName)))
      OPEN(FN+1,FILE=XYZFileName(:LEN_TRIM(XYZFileName)))
      OPEN(FN+2,FILE='movie.mos')
      OPEN(FN+3,FILE='DefAtEnergy.dat')
      CALL WriteMolMolXYZFile(FN+1,Time,Iter)
      write(OPF10,*)time

      WRITE(FN,*) 'load xyz ', XYZFileName
      WRITE(FN,*) 'color temperature'
      WRITE(FN,*) 'set background white'

      WRITE(FN,*) 'set boundbox true'
      WRITE(FN,*) 'zoom 140'
      WRITE(FN,*) 'rotate x 41'
      WRITE(FN,*) 'rotate y 4'
C      WRITE(FN,*) 'spacefill 45'
      WRITE(FN,*) 'spacefill 2'
      WRITE(FN,*) 'color axes black'
      WRITE(FN,*) 'wireframe 2'
!      WRITE(FN,*) 'select atomno<',NumS
      WRITE(FN,*) 'color atom green'
!      WRITE(FN,*) 'select atomno>=',NumS
      WRITE(FN,*) 'wireframe 12'
      WRITE(FN,*) 'color atom purple'
      WRITE(FN,*) 'set bondmode and'

32    FORMAT('select atomno = ',I5)
33    FORMAT('select atomno = ',I5,', atomno = ',I5)
34    FORMAT(I5,2X,5I5,2X,F8.4)

      CLOSE(FN)
      WRITE(FN+3,*)''
      CLOSE(FN+3)
      WRITE(FN+2,*)'source    ',SourceFileName
      WRITE(FN+2,*)'write gif ',GifFileName
      WRITE(FN+2,*)'zap'
      CLOSE(FN+2)
      OPEN(FN+2,FILE='movie.mos',STATUS='OLD',POSITION='APPEND')
      OPEN(FN+3,FILE='DefAtEnergy.dat',STATUS='OLD',POSITION='APPEND')

      END SUBROUTINE WriteMolMolSource
******|****************************************************************
