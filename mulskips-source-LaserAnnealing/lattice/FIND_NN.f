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
**    This subroutine find remanent two next-neighbour site of SiteC when
**    two next-neighbour are given
**************************************************************************************
      SUBROUTINE FIND_NN(Site1,Site2,SiteC,NewSite,LX,LY)
      IMPLICIT NONE
      REAL(8) :: Pi=3.14159265358979323
      REAL(8)  :: Cth,Sth,snorm,Rot(3,3),ax(3)
      REAL(8)  :: RNewsite(3)
      INTEGER :: i,k,l,irot,LX,LY
      INTEGER :: Site1(3),Site2(3),SiteC(3),Newsite(3,2)
      INTEGER :: Site11(3),Site22(3)
      Site11=Site1
      Site22=Site2
! rotation axis as versor of the direction SiteC-Site1
      IF(IABS(SiteC(1)-Site11(1)).GT.LX/2)
     >   Site11(1)=Site11(1)+ISIGN(LX,SiteC(1)-Site11(1))
      IF(IABS(SiteC(2)-Site11(2)).GT.LY/2)
     >   Site11(2)=Site11(2)+ISIGN(LY,SiteC(2)-Site11(2))

      IF(IABS(SiteC(1)-Site22(1)).GT.LX/2)
     >   Site22(1)=Site22(1)+ISIGN(LX,SiteC(1)-Site22(1))
      IF(IABS(SiteC(2)-Site22(2)).GT.LY/2)
     >   Site22(2)=Site22(2)+ISIGN(LY,SiteC(2)-Site22(2))

      ax=DFLOAT(SiteC)-DFLOAT(Site11)
      snorm=DSQRT(ax(1)**2+ax(2)**2+ax(3)**2)
      ax=ax/snorm
      DO irot=1,2
! rotation matrix for an angle 2/3 Pi or 4/3 Pi
         Cth=DCOS(DFLOAT(2*irot)*pi/3.)
         Sth=DSIN(DFLOAT(2*irot)*pi/3.)
         DO i=1,3
            Rot(i,i)=ax(i)**2*(1-Cth)+Cth
         END DO
         Rot(2,1)=ax(1)*ax(2)*(1-Cth)-ax(3)*Sth
         Rot(1,2)=ax(1)*ax(2)*(1-Cth)+ax(3)*Sth
         Rot(3,1)=ax(1)*ax(3)*(1-Cth)+ax(2)*Sth
         Rot(1,3)=ax(1)*ax(3)*(1-Cth)-ax(2)*Sth
         Rot(3,2)=ax(2)*ax(3)*(1-Cth)-ax(1)*Sth
         Rot(2,3)=ax(2)*ax(3)*(1-Cth)+ax(1)*Sth
! rotate Site2
         DO k=1,3
            RNewSite(k)=0.0
            DO l=1,3
               RNewSite(k)=RNewSite(k)+Rot(l,k)*
     >                   (DFLOAT(Site22(l))-DFLOAT(SiteC(l)))
            END DO
         END DO
         NewSite(1:3,irot)=NINT(RNewSite)+SiteC
         IF(NewSite(1,irot).GE.LX)NewSite(1,irot)=NewSite(1,irot)-LX
         IF(NewSite(1,irot).LT.0)NewSite(1,irot)=NewSite(1,irot)+LX
         IF(NewSite(2,irot).GE.LY)NewSite(2,irot)=NewSite(2,irot)-LY
         IF(NewSite(2,irot).LT.0)NewSite(2,irot)=NewSite(2,irot)+LY
!         write(*,*)'FindNN',NewSite(1:3,irot),RNewSite
      END DO
      END SUBROUTINE FIND_NN
