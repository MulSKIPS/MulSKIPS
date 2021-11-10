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
      real(8) function Random(idum)

       Real(8) :: Fac
       INTEGER :: Idum
       INTEGER :: Mbig,Mseed,Mz,Iff,Mj,Mk,I,II,Ma,K,Inext,Inextp


         Parameter (mbig=1000000000,Mseed=161803398,Mz=0,fac=1./Mbig)

         Dimension MA(55)  !the dimension of MA is  special
         save
!       inizialize MA(55) using idum and mseed
         if (idum.lt.0.or.iff.eq.0) then
          iff=1
          mj=mseed-iabs(idum)
          mj=mod(mj,mbig)
          ma(55)=mj
          mk=1
!        inizialize the rest of MA in slightly random order with
!        Number that are not especially random
          do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if (mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
11        continue
!         randomize MA NumBer by the "warming up generator"
          do 13 k=1,4
            do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if (ma(i).lt.mz)ma(i)=ma(i)+mbig
12          continue
13        continue
          inext=0
          inextp=31     ! 31 is special
          idum=1
        end if
        inext=inext+1
        if (inext.eq.56) inext=1
         inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)     ! generate a new random NumBer
        if (mj.lt.mz)mj=mj+mbig     ! Be sure that it is in the range
        ma(inext)=mj                ! Store it
        random=mj*fac                 ! outpup ran3 belongs (0,1)

      end function random
