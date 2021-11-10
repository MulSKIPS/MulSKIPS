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
**    This subroutine find Arm-Chair and Zig-Zag configurations which next-neighbour
**    of the site of Site and Second Neighbour of the Site SiteC. Site C has a
**    fixed Coordination
**************************************************************************************
      SUBROUTINE FIND_SN(Site,SiteC_inp,SNSite_inp,SNZiGZag,
     >                   SNArmChair,LX,LY)
      IMPLICIT NONE
      REAL(8)  :: snorm,ax(3),Cent(3)
      REAL(8)  :: RSite(3),RSiteC(3)
      REAL(8)  :: RSNZiGZag(3),RSNArmChair(3)
      INTEGER :: i,LX,LY,Site(3),SiteC_inp(3)
      INTEGER :: SNSite_inp(3,3),SNZiGZag(3,3),SNArmChair(3,3)
      INTEGER :: SNSite(3,3),SiteC(3)
      LOGICAL debug_bc

      SNSite=SNSite_inp
      SiteC=SiteC_inp
!      write(*,*)'standard',Site,SiteC,SNSite
      debug_bc=.FALSE.
      IF(IABS(Site(1)-SiteC(1)).GT.LX/2)THEN
         SiteC(1)=SiteC(1)+ISIGN(LX,Site(1)-SiteC(1))
         debug_bc=.TRUE.
      END IF
      IF(IABS(Site(2)-SiteC(2)).GT.LY/2)THEN
         SiteC(2)=SiteC(2)+ISIGN(LY,Site(2)-SiteC(2))
         debug_bc=.TRUE.
      END IF
      DO i=1,3
        IF(IABS(Site(1)-SNSite(1,i)).GT.LX/2)THEN
           SNSite(1,i)=SNSite(1,i)+ISIGN(LX,Site(1)-SNSite(1,i))
           debug_bc=.TRUE.
        END IF
        IF(IABS(Site(2)-SNSite(2,i)).GT.LY/2)THEN
           SNSite(2,i)=SNSite(2,i)+ISIGN(LY,Site(2)-SNSite(2,i))
           debug_bc=.TRUE.
        END IF
      END DO
!      if(debug_bc)write(*,*)'debug_bc',Site,SiteC,SNSite
      Cent=0.5*(DFLOAT(SiteC)+FLOAT(Site))
      ax=DFLOAT(Site)-DFLOAT(SiteC)
      snorm=DSQRT(ax(1)**2+ax(2)**2+ax(3)**2)
      ax=ax/snorm
!      write(*,*)ax
      DO i=1,3
         RSNZigZag=2*Cent-DFLOAT(SNsite(1:3,i))
         RSNArmChair=DFLOAT(SNsite(1:3,i))+8.6602254037*ax
         SNZiGZag(1:3,i)=NINT(RSNZigZag)
         SNArmChair(1:3,i)=NINT(RSNArmChair)
         IF(SNZiGZag(1,i).GE.LX)SNZiGZag(1,i)=SNZiGZag(1,i)-LX
         IF(SNZiGZag(1,i).LT.0)SNZiGZag(1,i)=SNZiGZag(1,i)+LX
         IF(SNZiGZag(2,i).GE.LY)SNZiGZag(2,i)=SNZiGZag(2,i)-LY
         IF(SNZiGZag(2,i).LT.0)SNZiGZag(2,i)=SNZiGZag(2,i)+LY
         IF(SNArmChair(1,i).GE.LX)SNArmChair(1,i)=SNArmChair(1,i)-LX
         IF(SNArmChair(1,i).LT.0)SNArmChair(1,i)=SNArmChair(1,i)+LX
         IF(SNArmChair(2,i).GE.LY)SNArmChair(2,i)=SNArmChair(2,i)-LY
         IF(SNArmChair(2,i).LT.0)SNArmChair(2,i)=SNArmChair(2,i)+LY
      END DO
!      write(*,*)'standard Arm',SNArmChair
!      write(*,*)'standard Zig',SNZiGZag
!      if(debug_bc)write(*,*)'debug_bc Arm',SNArmChair
!      if(debug_bc)write(*,*)'debug_bc Zig',SNZiGZag

      END SUBROUTINE FIND_SN
