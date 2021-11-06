!   
!    This program calculates the matrixes of kernels.
! This is a new version of my code, only four different type of particles are distinguished:
!        water, pristine ice, snow and densely rimed graupel (with constan and size dependent density). 
!
!    Collision between the water drops (both with and without break up)
!    also collision of 'clean ' and 'polluted' drops for tracking the aerosols inside of the drops
!     break up kernels are also calculated here  
!    Collision between the ice crystals and water drops  
!    Aggregation of pristine ice crystals (results in snow) 
!    Collision between the ice crystals and the snow particles
!    Collision between the snow particles  and water drops
!    Self aggregation of snow crystals  
!    Collision between the high density graupel particles and water drops 
!    Collision between the high density graupel particles and pristine ice
!    Collision between the high density graupel particles and snow    
!    Terminal velocities of the diferent type of particles to calculate the fall out
!    terms in time dependent models, the value of the velocities are given in the midle of the bins.
!    Height corrections with (ro0/rp)**0.5 is necessery 
!    Precalculations to give the g(R) functions in the integral for ice multiplication.
!         int1=int(E(Mg,md),md=mk..mk+1), int2=int(md*E(Mg,md),md=mk..mk+1)
!    scavenging of the AgI particles by different processes for water drops, snow flakes, pristine ice and graupel particles  
!    The number of bins is 36, pk=2  
! 
      implicit real*8 (a-h,o-z)            
      real*8 m(40),r(40),m1(55),raer(55),PK,PI,R1,RHOAER,RHOW
      character filen*40
      parameter(k4=36,k8=55,rhow=1000.,r1=1.5625d-6)
      integer jgraupel,jwater
!  k4 = number of categories.
      common/m/ m
      common/m1/m1
      COMMON /RHOAER/ RHOAER
      
!  m(k) = mass of drops in the k'th category (kg)
      pk=2.d0
!  pk = width parameter
!  r(1) - radius of drop in first category in meters
!  m(1) - mass of drop in first category in kgr.
! 
      RHOAER = 1750.0
      pi=asin(1.d0)*2.d0
      r(1)=r1
      m(1)=4.d0*pi*rhow*r(1)**3/3.d0
      do k=1,k4+1
         km1=k-1
         m(k)=m(1)*pk**km1
         r(k)=r(1)*pk**(km1/3.d0)
         write(*,*)k, m(k), r(k)
      enddo
      m1(1) = m(1)/(pk**19)
      raer (1) = (3.0*m1(1)/(4.0*RHOAER*pi))**(1.0/3.0)
      do k = 1, k8
       km1=k-1
       m1(k) = m1(1)*pk**km1
       raer(k)= raer(1)*pk**(km1/3.d0)
      enddo
      open(unit=98,file='debug_kernel.txt')
!  
!  calculation of the collection kernels of water drop - water drop collision
!
!      coefdd.dat contains only two arrays ccc(36,36,12) and cc(36,3)
!      coefdd1.dat contains one more arry c5(36,36,4), this array is used in the case of coal. of two differnt type of water 
      
      jwater =1
      filen='coefdd.dat'
      open(3,file=filen,form='unformatted')        
      rewind 3
      call gencoef(k4,jwater)
      close(3)
      jwater=2
      filen='coefdd1.dat'                             ! this file contains c5(36,36,4) array
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoef(k4,jwater)
      close(3)


!      calculation of the collection (collision*coalescence) kernels of water drop - water drop collision
!      coalescence efficiency is given by Low and List
!      coefdd.dat contains only two arrays ccc(36,36,12) and cc(36,3)
!      coefdd1.dat contains one more arry c5(36,36,4)   , it is for tracking the aerosol particles inside of the water drops 
!      
!      coefddbr1.dat  Low and List eq. is used if dsmall > 200 and  otherwise ecoal =1.0 (same is used by Axel)
!      coefddbr1_new  break up based on the Staub et al paper
      write(*,*) 'open coefddbr1_new'
      jwater =1
      filen='coefddbr1_new.dat'
      open(3,file=filen,form='unformatted')        
      rewind 3
      write(*,*) 'before call gencoefbr1'
      call gencoefbr1(k4,jwater)
      close(3)
      write(*,*) 'open coefddbr2_new'
      jwater =1
      filen='coefddbr2_new.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoefbr2
      close(3)
      jwater =2
      filen='coefddbr11_new.dat'                    ! this file contains c5(36,36,4) array
      open(3,file=filen,form='unformatted')
      rewind 3
      write(*,*) 'before call gencoefbr1'
      call gencoefbr1(k4,jwater)
      close(3)
!
!  calculation of the collection kernels  graupel particle - water drop collision
!
!  
!   high density graupel particle - water drop collision
!
     
       filen='coefgd.dat'      
       open(3,file=filen,form='unformatted')
       rewind 3
       call gencoef2(k4)
       close(3)
!
!  the density of the graupel particles depends on it size, if r < about .1 mm rhogr = 450.
!                                                           if r > about 1.0 mm rhogr = 800
       filen='coefgd1.dat'
       open(3,file=filen,form='unformatted')
       rewind 3
       call gencoef21(k4)
       close(3)
	
!
!     calculation of the collection kernels of pristine ice - water drop collision
!                   
      filen='coeficd.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoef4(k4)
      close(3)
!	  
!  calculation of the collection kernels of snow flakes - water drop collision 
!     
      filen='coefsnowd1.dat'      
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoef5(k4)
      close(3)    
!  
!  calculation of the collection kernels of pristine ice crystal - pristine ice crystal collision
!     
      filen='coefsnow1.dat'      
      open(3,file=filen,form='unformatted')
      call gencoef6(k4)
      close(3) 
!  
!  calculation of the collection kernels of snow flakes - snow flakes collision
!     
      filen='coefsnow2.dat'      
      open(3,file=filen,form='unformatted')  
      call gencoef7(k4)
      close(3)      
!  
!  calculation of the collection kernels of snow flakes - pristine ice collision
!      
      filen='coefsnowic.dat'      
      open(3,file=filen,form='unformatted')
      call gencoef8(k4)
      close(3)
!  
!  calculation of the collection kernels of graupel particle - pristine ice collision
!     
      filen='coefgric.dat'      
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoef9(k4)
      close(3) 
!  
!  calculation of the collection kernels of graupel particle - snow flake collision
!     
      filen='coefgrsnow.dat'      
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoef10(k4)
      close(3)
!
!  calculation of the terminal velocities with sizes equal to the mass at the middle of the bin
!   
      filen='termv1.dat'      
      open(3,file=filen,form='unformatted')
      rewind 3
      call terminalv(k4)
      close(3) 
! ice multiplication due to the collison of high dens. graupel and water drops
      filen='icmulp.dat'      
      open(3,file=filen,form='unformatted')
      rewind 3
      call icemulp(k4)
      close(3)
! ice multiplication due to the collison of variable density graupel and water drops
      filen='icmulp1.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call icemulp1(k4)
      close(3)
! scavenging of the AgI particles  by water drops
      filen='coefda.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call gencoef1(k4)
      close(3)
      continue
      filen='brownian.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call genbrown(k4)
      close(3)
      filen='phoresis.dat'
      open(3,file=filen,form='unformatted')
      call genphoresis(k4)
      rewind 3
      close(3)
      filen='turbcoag.dat'
      open(3,file=filen,form='unformatted')
      call genturbcoag(k4)
      rewind 3
      close(3)
! scavenging of the AgI particles  by pristine ice
      filen='brownian_pice.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call genbrown_pice(k4)
      close(3)
      filen='phoresis_pice.dat'
      open(3,file=filen,form='unformatted')
      call genphoresis_pice(k4)
      rewind 3
      close(3)
      filen='turbcoag_pice.dat'
      open(3,file=filen,form='unformatted')
      call genturbcoag_pice(k4)
      rewind 3
      close(3)
! scavenging of the AgI particles  by snowflakes
      filen='brownian_snow.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call genbrown_snow(k4)
      close(3)
      filen='phoresis_snow.dat'
      open(3,file=filen,form='unformatted')
      call genphoresis_snow(k4)
      rewind 3
      close(3)
      filen='turbcoag_snow.dat'
      open(3,file=filen,form='unformatted')
      call genturbcoag_snow(k4)
      rewind 3
      close(3)
! scavenging of the AgI particles  by graupel
      filen='brownian_gr.dat'
      open(3,file=filen,form='unformatted')
      rewind 3
      call genbrown_gr(k4)
      close(3)
      filen='phoresis_gr.dat'
      open(3,file=filen,form='unformatted')
      call genphoresis_gr(k4)
      rewind 3
      close(3)
      filen='turbcoag_gr.dat'
      open(3,file=filen,form='unformatted')
      call genturbcoag_gr(k4)
      rewind 3
      close(3)
      end          

!         END of the main program      
      subroutine gencoef(k4,jwater)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,3)
!    and if it is necessary one more array c5(k,i,4)
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,3),ccc(36,36,12)
      real*8 c5(36,36,4)
      real*8 m(40)
      integer jwater
      common/m/ m
      common/k/ k
      external f1,f2,f3,f4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,3
        do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo
      do i=1,12
        do j=1,k4
          do k=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do j=1,k4
          do k=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4 
      errest=0.d0
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1) 
      write(*,*) '11',k,i
      call dtwodq(f1,bot,top,g1,h1,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(f2,bot,top,g1,h1,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,2)=result     
      call dtwodq(f3,bot,top,g1,h1,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,3)=result     
      call dtwodq(f4,bot,top,g1,h1,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*) '21',k,i           
      call dtwodq(f1,bot,top,g2,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(f2,bot,top,g2,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(f3,bot,top,g2,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(f4,bot,top,g2,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo      
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)
      write(*,*) '31',k,i
      call dtwodq(f1,bot,top,h1,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,9)=result    
      call dtwodq(f2,bot,top,h1,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,10)=result 
      call dtwodq(f3,bot,top,h1,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,11)=result     
      call dtwodq(f4,bot,top,h1,h2,errabs,errrel,irule,result,          &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*) '21',k,i           
      call dtwodq(f1,bot,top,h1,g2,errabs,errrel,irule,result,          &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(f2,bot,top,h1,g2,errabs,errrel,irule,result,          &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(f3,bot,top,h1,g2,errabs,errrel,irule,result,          &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(f4,bot,top,h1,g2,errabs,errrel,irule,result,          &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo   
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(f1,bot,top,g3,h1,errabs,errrel,irule,result,          &
     &            errest)
      cc(k,1)=result
      call dtwodq(f2,bot,top,g3,h1,errabs,errrel,irule,result,          &
     &            errest)
      cc(k,2)=result
      call dtwodq(f4,bot,top,g3,h1,errabs,errrel,irule,result,          &
     &            errest)
      cc(k,3)=result
      enddo
      if (jwater.eq.1) then
       write(3)cc
       write(3)ccc
      else
       write(3)cc
       write(3)ccc
       write(3)c5
      endif 
      close (3) 
      write(*,*) 'end of calc. of kernel for water -water'
      return
      end
!
      real*8 function f1(x,y)
      implicit none
      real*8 x,y 
      real*8 pi
      real*8 vx,rx,vy,ry,rn,rk,rnt,rkt,p,ec,e10,e20,dedp1,dedp2
      integer jump,k1,k2  
      real*8 e0(0:12,0:20)
      common /collwd/e0
!    The collection efficiencies are calculated on the base of the work of
!    Hall , 1980:A detailed microphysical model within a two-dimensional
!    dynamic framework: model description and preliminary results.
!    J. Atm. Sci. 37, 2486-2507
!
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419
      jump=2
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvwd(y,vy,ry)
      ec=0.d0
      rn=rx*1.d6
      rk=ry*1.d6
      if (rn.lt.rk) then
      rnt=rn
      rkt=rk
      rn=rkt
      rk=rnt
      endif
      if (jump.eq.1) then
       if(rn.ge.50.d0)then
        ec=1.d0
        go to 40
       endif
       if(rk.le.3.d0.and.rn.lt.50.d0)then
        ec=0.d0
        go to 40
       endif
       ec=4.5d-4*rn*rn*(1.d0-3.d0/rk)
       goto 40
      endif
!    if jump=2 then
      p=rk/rn
      if (rn.ge.300.) then
         ec=1.d0
         goto 40
      endif 
      if (rn.le.10.d0.and.p.le.0.05d0) then
         ec=0.0001d0
         goto 40
      endif
      if (rn.lt.5.d0) then
         ec=0.d0
         goto  40
      endif  
      if (p.lt.0.05) then
        k1=13
34      k1=k1-1
        if (rn.gt.e0(k1,0)) goto 34
        dedp1=(e0(k1,2)-e0(k1,1))/0.05d0
        dedp2=(e0(k1+1,2)-e0(k1+1,1))/0.05d0
        e10=e0(k1,1)-dedp1*(0.05d0-p)
        e20=e0(k1+1,1)-dedp2*(0.05d0-p)
        if(e10.lt.0.d0) e10=0.d0
        if(e20.lt.0.d0) e20=0.d0
        ec=(e10-e20)/(e0(k1,0)-e0(k1+1,0))*(rn-e0(k1,0))+e10
        goto 40
      endif
      k1=13
38    k1=k1-1
      if (rn.gt.e0(k1,0)) goto 38
      k2=0
39    k2=k2+1
      if (p.gt.e0(0,k2)) goto 39  
      e10=(p-e0(0,k2-1))*(e0(k1,k2)-e0(k1,k2-1))/0.05d0+e0(k1,k2-1)
      e20=(p-e0(0,k2-1))*(e0(k1+1,k2)-e0(k1+1,k2-1))/0.05d0+            &
     &      e0(k1+1,k2-1)
      ec=(e10-e20)/(e0(k1,0)-e0(k1+1,0))*(rn-e0(k1,0))+e10
40    continue
      if(ec.lt.0.d0)ec=0.d0
      f1=pi*(rx+ry)**2.d0*ec*dabs(vx-vy)
      return
      end
      function f2(x,y)
      implicit real*8 (a-h,o-z)
      f2=y*f1(x,y)
      return
      end
      function f3(x,y)
      implicit real*8 (a-h,o-z)
      f3=x*f1(x,y)
      return
      end
      function f4(x,y)
      implicit real*8 (a-h,o-z)
      f4=y*x*f1(x,y)
      return
      end
!
!    calculation of the kernels  coalesence efficiency is not equal to 1  
!
      subroutine gencoefbr1(k4,jwater)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,3)
!    and if it is necessary one more array c5(k,i,4)
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable.
!
      real*8 cc(36,3),ccc(36,36,12)
      real*8 c5(36,36,4)
      real*8 m(40)
      integer jwater
      common/m/ m
      common/k/ k
      external fbrl1,fbrl2,fbrl3,fbrl4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,3
        do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo
      do i=1,12
        do j=1,k4
          do k=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do j=1,k4
          do k=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      errest=0.d0
      irule=3
      do k=3,k4
      do i=1,k-2
      bot=m(i)
      top=m(i+1)
      write(*,*) '11br',k,i
      call dtwodq(fbrl1,bot,top,g1,h1,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fbrl2,bot,top,g1,h1,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,2)=result     
      call dtwodq(fbrl3,bot,top,g1,h1,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,3)=result     
      call dtwodq(fbrl4,bot,top,g1,h1,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*) '21br',k,i           
      call dtwodq(fbrl1,bot,top,g2,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fbrl2,bot,top,g2,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fbrl3,bot,top,g2,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fbrl4,bot,top,g2,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo      
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)
      write(*,*) '31',k,i
      call dtwodq(fbrl1,bot,top,h1,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,9)=result    
      call dtwodq(fbrl2,bot,top,h1,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,10)=result 
      call dtwodq(fbrl3,bot,top,h1,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,11)=result     
      call dtwodq(fbrl4,bot,top,h1,h2,errabs,errrel,irule,result,       &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*) '21',k,i           
      call dtwodq(fbrl1,bot,top,h1,g2,errabs,errrel,irule,result,       &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fbrl2,bot,top,h1,g2,errabs,errrel,irule,result,       &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fbrl3,bot,top,h1,g2,errabs,errrel,irule,result,       &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fbrl4,bot,top,h1,g2,errabs,errrel,irule,result,       &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo   
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(fbrl1,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest)
      cc(k,1)=result
      call dtwodq(fbrl2,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest)
      cc(k,2)=result
      call dtwodq(fbrl4,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest)
      cc(k,3)=result
      enddo
      if (jwater.eq.1) then
       write(3)cc
       write(3)ccc
      else
       write(3)cc
       write(3)ccc
       write(3)c5
      endif 
      close (3) 
      write(*,*) 'end for calc breauk mod. kernels'
      return
      end
      real*8 function fbrl1(x,y)
      implicit none
      real*8 x,y
      real*8 pi,tk,sigma
      real*8 vx,rx,vy,ry,rn,rk,rnt,rkt,p,ec,e10,e20,dedp1,dedp2
      real*8 cke,sc,web,deltast,etotal,ecoal
      integer jump,k1,k2
      real*8 e0(0:12,0:20)
      common /collwd/e0
!    The collection efficiencies are calculated on the base of the work of
!    Hall , 1980:A detailed microphysical model within a two-dimensional
!    dynamic framework: model description and preliminary results.
!    J. Atm. Sci. 37, 2486-2507
      jump=2
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvwd(y,vy,ry)
      ec=0.d0
      tk=273.15+20.
      sigma=(75.93-0.115*(tk-273.15))*1e-3
      rn=rx*1.d6
      rk=ry*1.d6
      if (rn.lt.rk) then
      rnt=rn
      rkt=rk
      rn=rkt
      rk=rnt
      endif
      if (jump.eq.1) then
       if(rn.ge.50.d0)then
        ec=1.d0
        go to 40
       endif
       if(rk.le.3.d0.and.rn.lt.50.d0)then
        ec=0.d0
        go to 40
       endif
        ec=4.5d-4*rn*rn*(1.d0-3.d0/rk)
       goto 40
      endif
!    if jump=2 then
      p=rk/rn
      if (rn.ge.300.) then
         ec=1.d0
         goto 40
      endif
      if (rn.le.10.d0.and.p.le.0.05d0) then
         ec=0.0001d0
         goto 40
      endif
      if (rn.lt.5.d0) then
         ec=0.d0
         goto  40
      endif
      if (p.lt.0.05) then
        k1=13
34      k1=k1-1
        if (rn.gt.e0(k1,0)) goto 34
        dedp1=(e0(k1,2)-e0(k1,1))/0.05d0
        dedp2=(e0(k1+1,2)-e0(k1+1,1))/0.05d0
        e10=e0(k1,1)-dedp1*(0.05d0-p)
        e20=e0(k1+1,1)-dedp2*(0.05d0-p)
        if(e10.lt.0.d0) e10=0.d0
        if(e20.lt.0.d0) e20=0.d0
        ec=(e10-e20)/(e0(k1,0)-e0(k1+1,0))*(rn-e0(k1,0))+e10
        goto 40
       endif
       k1=13
38     k1=k1-1
       if (rn.gt.e0(k1,0)) goto 38
       k2=0
39     k2=k2+1
       if (p.gt.e0(0,k2)) goto 39  
       e10=(p-e0(0,k2-1))*(e0(k1,k2)-e0(k1,k2-1))/0.05d0+e0(k1,k2-1)    
       e20=(p-e0(0,k2-1))*(e0(k1+1,k2)-e0(k1+1,k2-1))/0.05d0+           &
     &      e0(k1+1,k2-1)
       ec=(e10-e20)/(e0(k1,0)-e0(k1+1,0))*(rn-e0(k1,0))+e10
40    continue
      if(ec.lt.0.d0)ec=0.d0
      cke=0.5*x*y/(x+y)*(vx-vy)**2
      sc = 4.*pi*sigma*(rx**3 + ry**3)**(2./3.)
      web=cke/sc
      deltast = 4*pi*sigma*(rx**2 + ry**2) - sc
      etotal = cke + deltast
!      if (etotal.lt.5e-6) then
!       ecoal = 0.778*(rn/(rk+rn))**2*
!     1                       dexp(-2.61e6*sigma*etotal*etotal/sc)
!      else
!       ecoal = 0.0
!      endif
      ecoal = dexp(-1.15d0*web)      
      if (rk.lt.100) then   ! observation is valid outside of this size interval   (same condition is used by Axel)
!       ecoal =1.0
      endif
      if (ecoal.gt.1) ecoal =1.0
      fbrl1=pi*(rx+ry)**2.d0*ec*ecoal*dabs(vx-vy)
      return
      end
      function fbrl2(x,y)
      implicit real*8 (a-h,o-z)
      fbrl2=y*fbrl1(x,y)
      return
      end
      function fbrl3(x,y)
      implicit real*8 (a-h,o-z)
      fbrl3=x*fbrl1(x,y)
      return
      end
      function fbrl4(x,y)
      implicit real*8 (a-h,o-z)
      fbrl4=y*x*fbrl1(x,y)
      return
      end
!
!   This program calculates the probality function for the break up process the based on paper of Straub_et_al. 2009 and Schlottke et al 2010 .
! 
      SUBROUTINE GENCOEFBR2
      IMPLICIT NONE
      REAL*8 M(40),R(40),D(40),R1,PK,RHOW,PI,RX,RY,RX1,RY1,DSMALL,      &
     &       DLARGE,MSMALL,MLARGE,DDELTA
      REAL*8 QTN(36,36,36),PFN(36,36,36),PTSEGF(20,20,37)
      REAL*8 CCCMEAN(36,36)
      REAL*8 SIGMAW,TK,VX,VY,X1,X2,Y1,CKE,CKE0,ST,SC,WEB,WEBSTAR,CW,    &
     &      DELTAST,ETOTAL,ECOAL,XV,NRE
! These parameters are to calculate the fragmet concentration in different break up processes
      REAL*8 DCOAL,MCOAL,MNEW,DS,DL,SIGMA0,SIGMA1,SIGMA2,SIGMA3,DELTAD1,&
     & DELTAD2,DELTAD3,VAR,E1,N1,N2,N3,MU1,MU2,MU3,M31,M32,M33,M34,GAMMA
      REAL*8 SEG0,SEG1,SEG2,SEGG,DK1,DK2,SUMF,SUMF1,SUMF2,SUMF3,SUMM1,  &
     &  SUMM2,SUMM3,SUMM4,DMASSI,DMASSJ,MIST,MJST,MIV,MJV,DIV,DJV
      REAL*8 ERF
      INTEGER K4,K,KM1,MINBN,MAXBN,IMAX,J,I,IFLAGF1,IFLAGF2,IFLAGF3,    &
     &        IPRINT,J2,I2,IV,JV,ISEG,JSEG,IVMAX,JVMAX
      CHARACTER FILEN*40
! dummies
      real*4 pfk1(36),pfk2(36),pfk3(36),pfk4(36),fnt
      real*8 summ00,segf1,segf2,segf3,segf4,QTN1(36,36,36)
!
      PARAMETER(K4=36,RHOW=1000.)
!
!   formats
!
 600  FORMAT(5X,'FRAGMENTS SIZE DISTRIBUTIONS AT DIFFERENT DS AND DL')
 610  FORMAT(/2X,'DIAMETER OF THE SMALLER DROP: ',F7.4,' (cm)',2X,      &
     &    'DIAMETER OF THE LARGER DROP: ',F7.4,' (cm)'/2X,              &
     &    'CKE: ',F6.3,' (uJ)',2X,'ST: ',F6.3,' (uJ)',2X,'SC: ',F6.3,   &
     &     2X,' (uJ)',2X,'ETOTAL: ',F6.3,' (uJ)',2X,'ECOAL: ',F5.2,     &
     &    ' WEB. NUM',F6.3,' CW',F6.3/)
 620  FORMAT(2X,'PARAMETRES FOR THE  FRAGMENTS'/2X,                     &
     &    'NUMBER OF FRAGMENTS (FIRST, SECOND AND THIRD SIZE RANGE):',  &
     &    3(2X,F5.2)/2X,' MOD1, MOD2, MOD3',3(1X,F6.3),                 &
     &    ' SIGMA1, SIGMA2, SIGMA3',3(1X,F6.3))
 625  FORMAT(2X,'MASS OF COLL. DROPS ',E13.4,' MASS OF THE FRAGMENTS ', &
     &      3(E13.4,2X),' SUM MASS OF FRAGMENTS ',E13.4)
 630  FORMAT (/2X,'SIZE DISTRIBUTION OF THE FRAGMENTS'/                 &
     &      2X,'BIN N',2X,' D (cm)',9X,'FI1',11X,'FI2',11X,'FI3',11X,   &
     &      'FI4')
 640  FORMAT (3X,I3,4X,F7.4,2X,4(2X,E12.5))
 650  FORMAT (2X,'CHECK OF MASS CONSERVATION'/2X,'FIRST SIZE RANGE ',   &
     &     2(E13.4,2x),' SECOND SIZE RANGE ', 2(E13.4,2X),              &
     &     ' THIRD SIZE RANGE ', 2(E13.4,2X),                           &
     &     ' FOURTH SIZE RANGE ', 2(E13.4,2X))
 660  FORMAT (/5X,'SUM OF NUMBER OF FRAGMENTS: ',e13.4 )
!  K4 = Number of categories.
      COMMON/M/ M
      OPEN(unit=13,file='fragmentsd.txt')
      OPEN(unit=77,file='debug.txt')
      open(unit=99,file='test.txt')
!  M(k) = mass of drops in the k'th category (kg)
      PK=2.d0
!  pk = width parameter
!  R(K) - radius of drop in meters
!  M(K) - mass of drop in mass
      PI=ASIN(1.d0)*2.d0
      R(1)= (3.d0*M(1)/(4.d0*PI*RHOW))**(1.d0/3d0)
      M(1)=4.d0*PI*RHOW*R(1)**3/3.d0

      DO K=1,K4+1
       KM1 = K-1
       R(K)=R(1)*PK**(KM1/3.d0)
       D(K)=2.*R(K)
      ENDDO
!       MINBN= 19    ! minimum drop size (r ~ 0.01 cm) could be involved in  break up process
!       MAXBN = K4-MINBN
!       IMAX = MAXBN + 1
      IVMAX = 20
      JVMAX = 20
      TK = 293.15
      SIGMAW=(75.93-0.115*(TK-273.15))*1e-3
      SEG0=(6.d0/(PI*RHOW))**(1./3.)
      SEG1=PI/6.d0        ! rhow = 1 g/cm3
! clean the three dimensional array
      DO I = 1, K4
       DO J = 1, K4
        DO K = 1, K4
         QTN(I,J,K) = 0.d0
        ENDDO
       ENDDO
      ENDDO
      SUMM1=0.d0
      SUMM2=0.d0
      SUMM3=0.d0
      SUMM4=0.d0
      summ00=0.0
      DO 500 I = 1, K4-1     ! loop for the larger drops
!      ISEG = I-1+MINBN
       MIST = M(I)
       DMASSI = (M(I+1) - M(I))/IVMAX
       DO 550 J = 1, I      ! loop for the smaller drops
!         JSEG = J - 1 + MINBN
        IF (J.LT.19) GOTO 550      ! if the radius of the drop is smaller than 100 um nobreakup occures
         MJST = M(J)
         DMASSJ = (M(J+1) - M(J))/JVMAX
! start of the numerical integration of the QTN at i,j
         DO IV = 1,IVMAX
          DO JV = 1,JVMAX
           DO K = 1,K4
            PTSEGF(IV,JV,K) = 0.d0
           ENDDO
          ENDDO
         ENDDO
         DO 200 JV =1,JVMAX
          MJV = MJST + 0.5*(2*JV -1)*DMASSJ    
          CALL TVWD (MJV,VX,RX)
          DJV = 2.d0*RX
          DO 300 IV = 1,IVMAX
           IPRINT=0
           IF (I.EQ.28.AND.J.EQ.24) THEN
            IF (IV.EQ.10.AND.JV.EQ.10) THEN
             IPRINT=1
            ENDIF
           ENDIF
           MIV = MIST + 0.5*(2*IV - 1)*DMASSI
           CALL TVWD (MIV,VY,RY)
           DIV = 2.d0*RY
           IFLAGF1 = 0
           IFLAGF2 = 0
           IFLAGF3 = 0
           CKE=0.5d0*MIV*MJV/(MIV+MJV)*(VX-VY)**2.0d0
           SC = 4.d0*PI*SIGMAW*(RX**3.d0 + RY**3.d0)**(2.d0/3.d0)
           WEB = CKE/SC
           ST = 4.d0*PI*SIGMAW*(RX**2.d0 + RY**2.d0)
           DELTAST = ST - SC
           ETOTAL = CKE + DELTAST
           ECOAL = DEXP(-1.15d0*WEB)
           IF(DIV .GT. DJV) THEN
            DSMALL = DJV
            DLARGE = DIV
            MSMALL = MJV
            MLARGE = MIV
           ELSE
            DSMALL = DIV
            DLARGE = DJV
            MSMALL = MIV
            MLARGE = MJV
           ENDIF
           DS = 100.*DSMALL        ! drop size in cm
           DL = 100.*DLARGE        ! drop size in cm
           DCOAL = (DS**3 + DL**3)**(1./3.)       ! size of the drop formed due to the collision
           GAMMA = DL/DS
           CW = CKE*WEB*1d6                            ! the unit of the cw is uJ
           IFLAGF1 =0
           IFLAGF2 =0
           IFLAGF3 =0
!   fragment formed in the first size range (smallest fragments)
           DELTAD1 = 0.0125d0*DSQRT(CW)
           E1 = 0.04
           VAR = DELTAD1**2.d0/12.d0
           SIGMA1=DSQRT(DLOG(VAR/E1**2d0+1.d0))
           MU1 = DLOG(E1)-0.5d0*SIGMA1**2.d0
           IF (GAMMA*CW.GT.7.0d0) THEN
            N1 = 0.088*(GAMMA*CW - 7.0d0)
            IFLAGF1 =1
           ELSE
            N1 = 0.0
           ENDIF
!   fragment formed in the second size range (slightly larger fragments)
           IF (CW.GT.21.d0) THEN
            N2 = 0.22d0*(CW - 21.d0)
            DELTAD2 = 7.0d-3*(CW - 21.d0)
            IFLAGF2 =1
           ELSE
            N2 = 0.d0
            DELTAD2 = 0.d0
           ENDIF
           MU2 = 9.5d-2
           SIGMA2 = DSQRT(DELTAD2**2.d0/12.d0)
!  fragment formed in the third size range (size of the fragments are near to that of the smaller drop)
           IF (CW.LT.21.d0) THEN
             N3 = 1.0d0
             IFLAGF3 =1
           ELSEIF (CW.GT.21.d0.AND.CW.LT.46.d0) THEN
             N3 = 0.04d0*(46.d0 - CW)
             IFLAGF3=1
           ELSE
             N3 = 0.d0
           ENDIF
           DELTAD3 = 0.01d0*(1.d0 + 0.76d0*DSQRT(CW))
           SIGMA3 = DSQRT(DELTAD3**2.d0/12.d0)
           MU3 =0.9*DS
!  fragment formed in the fourth size range (near to the largeer drop)
           M31 = N1*DEXP(3.d0*MU1 + 9.d0*SIGMA1**2.d0)*SEG1
           M32 = N2*(MU2**3.d0 + 3.d0*MU2*SIGMA2**2.d0)*SEG1
           M33 = N3*(MU3**3.d0 + 3.d0*MU3*SIGMA3**2.d0)*SEG1
           M34 = SEG1*(DL**3.d0+DS**3.d0) - M31 - M32 - M33
           MCOAL = MSMALL + MLARGE
           SUMF = M31 + M32 + M33
           IF (IPRINT.EQ.1) THEN
            WRITE(13,610) DS,DL,CKE*1e6,ST*1e6,SC*1e6,ETOTAL*1e6,       &
     &                    ECOAL,WEB,CW
            WRITE(13,620) N1,N2,N3,MU1,MU2,MU3,SIGMA1,SIGMA2,SIGMA3
            WRITE(13,625) MCOAL,M31,M32,M33,M34,SUMF
            segg=(1+gamma**3)**(11.0/3.0)/(1+gamma)**2/gamma**6
            segg=segg*(1+gamma**2-(1+gamma**3)**(2./3.0))
            segg=segg*24./5.
            webstar=12.*web*(1+gamma**(-3.))**(5./3.)/gamma**(-2.0)
            write(13,*) 'webstar', webstar
            write(13,*)'gamma,E*w**,f(g)',gamma,ecoal*webstar,segg
            webstar=rhow*dsmall*(VX-VY)**2.0d0/sigmaw
            write(13,*) 'webstar2',webstar
           ENDIF
           SUMF1=0.0
           SUMF2=0.0
           SUMF3=0.0
           do k = 1, 36
            pfk1(k)= 0.0
            pfk2(k)= 0.0
            pfk3(k)= 0.0
            pfk4(k)= 0.0
           enddo
           SEG2 = (1. - ECOAL)*DMASSI*DMASSJ
           DO 400 K = 1, K4-1       ! loop for the fragments
            DK1 = SEG0*M(K)**0.33333*100.            ! bin boundary for fragments  in cm
            DK2 = SEG0*M(K+1)**0.33333*100.            ! bin boundary for fragments
!       distribution of the fragments from the first size range into the bins
            IF (IFLAGF1 .EQ. 0) GOTO 402
            X1=(DLOG(DK1)-MU1)/SQRT(2.0)/SIGMA1
            X2=(DLOG(DK2)-MU1)/SQRT(2.0)/SIGMA1
            PFK1(K)=N1*0.5d0*(ERF(X2) - ERF(X1))
            PTSEGF(IV,JV,K) = SEG2*PFK1(K)

!       distribution of the fragments from the second size range into the bins
402         CONTINUE
            IF (IFLAGF2 .EQ. 0) GOTO 403
            X1=(MU2-DK1)/SQRT(2.0)/SIGMA2
            X2=(MU2-DK2)/SQRT(2.0)/SIGMA2
            PFK2(K)=N2*0.5d0*(ERF(X1) - ERF(X2))
            PTSEGF(IV,JV,K) = PTSEGF(IV,JV,K)+ SEG2*PFK2(K)
!        distribution of the fragments from the third size range into the bins
403         CONTINUE
            IF (IFLAGF3.EQ.0) GOTO 404
            X1=(MU3-DK1)/SQRT(2.0)/SIGMA3
            X2=(MU3-DK2)/SQRT(2.0)/SIGMA3
            PFK3(K)=N3*0.5d0*(ERF(X1) - ERF(X2))
            PTSEGF(IV,JV,K) = PTSEGF(IV,JV,K)+ SEG2*PFK3(K)
404         CONTINUE
            SUMF1=SUMF1 + PFK1(K)
            SUMF2=SUMF2 + PFK2(K)
            SUMF3=SUMF3 + PFK3(K)
400        CONTINUE
           K = 1
410        K = K+1
           IF (M34*1e-3.GT.M(K)) GOTO 410
           if (iprint.eq.1) then
             write(*,*) i,j,iv,jv
             write(*,*)'k,d4',k,(6*m34/3.14159)**0.3333
           endif
           
           PFK4(K-1) = 1
!           PTSEGF(IV,JV,K) = PTSEGF(IV,JV,K)+ SEG2*PFK4(K)
           PTSEGF(IV,JV,K-1) = PTSEGF(IV,JV,K-1)+ SEG2*PFK4(K-1)
           DO 421 K =1,K4
             SUMM1 = SUMM1 + PFK1(K)*m31
             SUMM2 = SUMM2 + PFK2(K)*m32
             SUMM3 = SUMM3 + PFK3(K)*m33
             SUMM4 = SUMM4 + PFK4(K)*m34
             summ00=summ00+(1-ecoal)*(miv+mjv)
 421       continue
           
           IF (IPRINT.EQ.-1) THEN
            WRITE(13,630)
            SUMM1=0.d0
            SUMM2=0.d0
            SUMM3=0.d0
            SUMM4=0.d0
            DO 420 K =1,K4
             WRITE(13,640)K,D(K)*1e2,PFK1(K),PFK2(K),PFK3(K),PFK4(K)
             SUMM1 = SUMM1 + PFK1(K)*1.5d0*M(K)*1e3
             SUMM2 = SUMM2 + PFK2(K)*1.5d0*M(K)*1e3
             SUMM3 = SUMM3 + PFK3(K)*1.5d0*M(K)*1e3
             SUMM4 = SUMM4 + PFK4(K)*1.5d0*M(K)*1e3
420         CONTINUE
            WRITE(13,650) M31,SUMM1,M32,SUMM2,M33,SUMM3,M34,SUMM4
            WRITE(13,660) SUMF1+SUMF2+SUMF3
            write(13,*)
            do k=1,k4
             DDELTA = SEG0*(M(K+1)**0.3333d0 - M(K)**0.3333d0)*100.
             segf1=PFK1(K)/ddelta
             segf2=PFK2(K)/ddelta
             segf3=PFK3(K)/ddelta
             segf4=PFK4(K)/ddelta
             write(13,*) k,d(k)*1e2,segf1,segf2,segf3,segf4
            enddo
           ENDIF
300       CONTINUE      ! loop for IV (large drops)
200      CONTINUE      ! lop for JV (small drops)
! end of numreical integration over the domain given by the drop size D(I) and D(J)
         DO K = 1,K4
          DDELTA = SEG0*(M(K+1)**0.3333d0 - M(K)**0.3333d0)*100.        ! diameter difference in cm
          DO IV = 1,IVMAX
           DO JV = 1,JVMAX
            PFN(I,J,K) =PFN(I,J,K) + PTSEGF(IV,JV,K)
           ENDDO
          ENDDO
          QTN(I,J,K)=1.d0/M(J)/M(I)/M(K)*PFN(I,J,K)
          qtn1(i,j,k) =PFN(I,J,K)/M(J)/M(I)
         ENDDO
550    CONTINUE       ! loop for J small drops
500   CONTINUE         ! loop for I large drops
       write(13,*) 'summ1,summ2,summ3,summ4',summ1,summ2,summ3,summ4
       write(13,*) 'summ00',summ00*1e3
      FNT =0.0
      DO I = 1,K4
       DO J = 1,K4
        DO K=1,k4
         FNT = FNT + PFN(I,J,K)
        ENDDO
       ENDDO
      ENDDO
      WRITE(13,660) FNT
!          stop
!
!  calculation of the mean collection collision kernels of water drop - water drop collision
!           (K(x,y)= Kmean)
     
      CALL GENCOEF_CCMEAN(K4,CCCMEAN)
      WRITE(3) QTN,CCCMEAN
      CLOSE (3)
       write(*,*) 'at the end of the code for break up'
!
      END
      FUNCTION ERF(X)
      IMPLICIT NONE
      REAL*8 ERF,X
      REAL*8 GAMMP
      IF (X.LT.0) THEN
       ERF = -GAMMP(0.5d0,X**2)
      ELSE
       ERF = GAMMP(0.5d0,X**2)
      ENDIF
      RETURN
      END
      FUNCTION GAMMP(A,X)
      IMPLICIT NONE
      REAL*8 A,GAMMP,X
! Returns of the incomplete gamma function P(a,x)
      REAL*8 GAMMCF,GAMSER,GLN
      IF (X.LT.0..OR.A.LT.0) THEN
       WRITE(*,*) 'Bad aRgumeNts iN gammp'
      ENDIF
      IF (X.LT.A+1) THEN
       CALL GSER(GAMSER,A,X,GLN)
       GAMMP = GAMSER
      ELSE
       CALL GCF(GAMMCF,A,X,GLN)
       GAMMP = 1. - GAMMCF
      ENDIF
      RETURN
      END
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 A,GAMSER,GLN,X,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
!
! Returns the incomplete gamma function P(a,x) evaulated by its series representation as gamser.
! Also returns ln(Gamma(a)) as gln
      INTEGER N
      REAL*8 AP,DEL,SUM,GAMMLN
      GLN = GAMMLN(A)
      IF (X .LE. 0) THEN
       WRITE(*,*) 'X < 0 iN gseR'
       GAMSER =0.0
       RETURN
      ENDIF
      AP = A
      SUM = 1./A
      DEL = SUM
      DO 10 N =1,ITMAX
       AP = AP + 1
       DEL = DEL *X/AP
       SUM = SUM + DEL
       IF (DABS(DEL).LT. DABS(SUM)*EPS) GOTO 1
10    CONTINUE
1     GAMSER = SUM*DEXP(-X + A *DLOG(X) -GLN)
      RETURN
      END
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      IMPLICIT NONE
      INTEGER ITMAX
      REAL*8 A,GAMMCF, GLN,X,EPS,FPMIN
      PARAMETER(ITMAX = 100,EPS=3.e-7,FPMIN=1.e-30)
      INTEGER I
      REAL*8 AN,B,C,D,DEL,H,GAMMLN
      GLN = GAMMLN(A)
      B = X + 1.0 -A
      C =1./FPMIN
      D = 1./B
      H = D
      DO 10 I =1,ITMAX
       AN =-I*(I-A)
       B = B + 2.0
       D = AN*D + B
       IF (DABS(D).LT.FPMIN) D = FPMIN
       C = B + AN/C
       IF (DABS(C) .LT. FPMIN) C = FPMIN
       D = 1./D
       DEL = D * C
       H = H * DEL
       IF (DABS(DEL-1).LT.EPS) GOTO 1
10    CONTINUE
!  WRITE(*,*) 'A too large, itmax too small in gcf'
1     GAMMCF = DEXP(-X + A*DLOG(X) - GLN) * H
      RETURN
      END
      FUNCTION GAMMLN(XX)
      IMPLICIT NONE
      REAL*8 GAMMLN,XX,GAMM
      INTEGER J
      REAL*8 SER,STP,TMP,X,Y,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146d0,-86.50532032941677d0,            &
     &     24.01409824083091d0,-1.231739572450155d0,                    &
     &     0.1208650973866179d-2,-0.5395239384953d-5,                   &
     &     2.50662872746310005d0/
      X = XX
      Y = X
      TMP = X + 5.5d0
      TMP = (X + 0.5d0)*DLOG(TMP) - TMP
      SER = 1.000000000190015d0
      DO 10 J = 1,6
       Y = Y + 1.d0
       SER = SER + COF(J)/Y
10    CONTINUE
      GAMMLN = TMP + DLOG(STP*SER/X)
      GAMM = DEXP(GAMMLN)
      RETURN
      END
!
!   subroutine for calculation of the mean collison kernel
!
!   subroutine for calculation of the mean collison kernel
      SUBROUTINE GENCOEF_CCMEAN(K4,CCCMEAN)
      IMPLICIT NONE
!
!  calculates the analytical expressions for the following coef. :
!    CCMEAN(K,I)
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable.
!
      REAL*8 CCCMEAN(36,36)
      REAL*8 M(40),FBRUP2,H1,H2,RESULT
      REAL*8 ERRABS,ERRREL,ERREST,BOT,TOP
      INTEGER K4,K,I,IRULE
      COMMON/M/ M
      COMMON/K/ K
      EXTERNAL FBRUP2
      EXTERNAL H1,H2
!   dummies
      real*4 xx
!
!  cleaning  matrix
!
      DO K=1,K4
       DO I=1,K4
        CCCMEAN(K,I)=0.d0
       ENDDO
      ENDDO
      errabs=0.d0
      errrel=5.d-4
      errest=0.d0
      irule=3
      DO K=1,K4
       DO I=1,K4
        BOT=M(I)
        TOP=M(I+1)
        write(*,*) '31',k,i
        if (k.eq.18.and.i.eq.19) then
          xx=1
        endif
       CALL dtwodq(fbrup2,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &                                                        errest)
       CCCMEAN(K,I)=RESULT/M(K)/M(I)
       ENDDO
      ENDDO
      RETURN
      END

      real*8 function fbrup2(x,y)
      implicit none
      real*8 x,y
      real*8 pi
      real*8 vx,rx,vy,ry,rn,rk,rnt,rkt,p,ec,e10,e20,dedp1,dedp2
      integer jump,k1,k2
      real*8 e0(0:12,0:20)
      common /collwd/e0
      jump=2
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvwd(y,vy,ry)
      ec=0.d0
      rn=rx*1.d6
      rk=ry*1.d6
      if (rn.lt.rk) then
      rnt=rn
      rkt=rk
      rn=rkt
      rk=rnt
      endif
      if (jump.eq.1) then
       if(rn.ge.50.d0)then
        ec=1.d0
        go to 40
       endif
       if(rk.le.3.d0.and.rn.lt.50.d0)then
        ec=0.d0
        go to 40
       endif
       ec=4.5d-4*rn*rn*(1.d0-3.d0/rk)
       goto 40
      endif
!    if jump=2 then
      p=rk/rn
      if (rn.ge.300.) then
         ec=1.d0
         goto 40
      endif
      if (rn.le.10.d0.and.p.le.0.05d0) then
         ec=0.0001d0
         goto 40
      endif
      if (rn.lt.5.d0) then
         ec=0.d0
         goto  40
      endif
      if (p.lt.0.05) then
        k1=13
 34     k1=k1-1
        if (rn.gt.e0(k1,0)) goto 34
        dedp1=(e0(k1,2)-e0(k1,1))/0.05d0
        dedp2=(e0(k1+1,2)-e0(k1+1,1))/0.05d0
        e10=e0(k1,1)-dedp1*(0.05d0-p)
        e20=e0(k1+1,1)-dedp2*(0.05d0-p)
        if(e10.lt.0.d0) e10=0.d0
        if(e20.lt.0.d0) e20=0.d0
        ec=(e10-e20)/(e0(k1,0)-e0(k1+1,0))*(rn-e0(k1,0))+e10
        goto 40
      endif
      k1=13
 38   k1=k1-1
      if (rn.gt.e0(k1,0)) goto 38
      k2=0
 39   k2=k2+1
      if (p.gt.e0(0,k2)) goto 39
      e10=(p-e0(0,k2-1))*(e0(k1,k2)-e0(k1,k2-1))/0.05d0+e0(k1,k2-1)
      e20=(p-e0(0,k2-1))*(e0(k1+1,k2)-e0(k1+1,k2-1))/0.05d0+            &
     &      e0(k1+1,k2-1)
      ec=(e10-e20)/(e0(k1,0)-e0(k1+1,0))*(rn-e0(k1,0))+e10
 40   continue
      if(ec.lt.0.d0)ec=0.d0
      fbrup2=pi*(rx+ry)**2.d0*ec*dabs(vx-vy)
      return
      end

! end   for break up kernels


!
!  the h1,h2,g1,g2,g3 functions define the inner integral bounderies 
!  for the imsl integration function dtwodq
!
      function g1(x)
      implicit real*8 (a-h,o-z)
      real*8 m(40)
      common/m/ m
      common/k/ k
      g1=m(k)-x
      return
      end
      function g2(x)
      implicit real*8 (a-h,o-z)
      real*8 m(40)
      common/m/ m
      common/k/ k
      g2=m(k+1)-x
      return
      end
      function g3(x)
      implicit real*8 (a-h,o-z)
      real*8 m(40)
      common/m/ m
      common/k/ k
      g3=m(k-1)
      return
      end
      function h1(x)
      implicit real*8 (a-h,o-z)
      real*8 m(40)
      common/m/ m
      common/k/ k
      h1=m(k)
      return
      end
      function h2(x)
      implicit real*8 (a-h,o-z)
      real*8 m(40)
      common/m/ m
      common/k/ k
      h2=m(k+1)
      return
      end
      SUBROUTINE dtwodq(ff,a,b,y1,y2,errabs,errrel,irule,s,errest)
      implicit real*8 (a-h,o-z)
      external ff
      PARAMETER (EPS=5.d-4, JMAX=20)
      errabs=0.d0
      errrel=0.d0
!      irule=0
!      errest=0.d0
      ost=-1.d30
      os= -1.d30
      a1=a
      b1=b
      it=0 
      do 11 j=1,JMAX 
       if (j.eq.1) then
        y11=y1(a)
        y22=y2(a)
        call qsimpy(ff,y11,y22,a,fa)
        y11=y1(b)
        y22=y2(b)        
        call qsimpy(ff,y11,y22,b,fb)
        st=0.5d0*(b-a)*(fa+fb)        
       else
        it=2**(j-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 15 k=1,it
          y11=y1(x)
          y22=y2(x) 
          call qsimpy (ff,y11,y22,x,fx) 
          seg3=fx
          sum=sum+fx
          x=x+del
15      continue
        st=0.5d0*(st+(b-a)*sum/tnm)
       endif
       if(s.lt.1d-100.and.it.gt.1000) then
         s=0.d0
         return
       endif        
       s=(4.d0*st-ost)/3.d0 
       if (dabs(s-os).le.EPS*dabs(os)) then
         return 
       endif 
       os=s
       ost=st
11    continue
      return
      END
      SUBROUTINE qsimpy(ff,a,b,xx,s)
      implicit real*8 (a-h,o-z)
      external ff
      PARAMETER (EPS=5.d-4, JMAX=20)
      if (dabs(b-a).lt.1.d-40) then
       s=0.d0
       return
      endif
      it=0
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
       if (j.eq.1) then
        fa=ff(xx,a)
        fb=ff(xx,b)
        st=0.5d0*(b-a)*(fa+fb)
       else
        it=2**(j-2)
        tnm=it 
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 15 k=1,it
          fx=ff(xx,x)
          sum=sum+fx
          x=x+del
15      continue
        st=0.5d0*(st+(b-a)*sum/tnm)
       endif
       if(s.lt.1d-100.and.it.gt.1000) then
         s=0.d0
         return
       endif 
       s=(4.d0*st-ost)/3.d0  
       if (dabs(s-os).le.EPS*dabs(os))then
         return 
       endif 
       os=s
       ost=st
11    continue
      return  
      END
      SUBROUTINE qsimp(func,a,b,ra,s)
      INTEGER JMAX
      DOUBLE PRECISION a,b,func,ra,s,EPS
      EXTERNAL func
      PARAMETER (EPS=1.d-6, JMAX=20)
!     USES trapzd
      INTEGER j
      DOUBLE PRECISION os,ost,st
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
        call trapzd(func,a,b,ra,st,j)
        s=(4.d0*st-ost)/3.d0
        if (abs(s-os).le.EPS*abs(os)) return
        os=s
        ost=st
11    continue
      write(*,*) 'too many steps in qsimp',os,s
      END
!  (C) Copr. 1986-92 Numerical Recipes Software 'LME.
      SUBROUTINE trapzd(func,a,b,ra,s,n)
      INTEGER n
      DOUBLE PRECISION a,b,s,func,ra
      EXTERNAL func
      INTEGER it,j
      DOUBLE PRECISION del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a,ra)+func(b,ra))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x,ra)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
!
!   This routine calculates the collection kernel of the water drop - high density particle
!
      subroutine gencoef2(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,12)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external fd1,fd2,fd3,fd4,fdv1,fdv2,fdv3,fdv4
      external g1,g2,g3
      external h1,h2
      
!
!  cleaning all matrixes
!
      do i=1,4
        do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo
      do i=1,12
        do j=1,k4
          do k=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
       do j=1,k4
        do k=1,k4
         cccv(k,j,i)=0.d0 
         c5v(k,j,i)=0.d0
        enddo
       enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'12',k,i 
      call dtwodq(fd1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(fd2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result        
      call dtwodq(fd3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result        
      call dtwodq(fd4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result  
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'22',k,i      
      call dtwodq(fd1,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(fd2,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(fd3,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(fd4,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)     
      write(*,*)'32',k,i  
      call dtwodq(fd1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(fd2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(fd3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(fd4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(fd1,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,1)=result
      call dtwodq(fd2,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,2)=result   
      call dtwodq(fd3,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,3)=result      
      call dtwodq(fd4,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,4)=result
      enddo
!        next,  the water drops are larger than the graupel particles 
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'12v',k,i
      call dtwodq(fdv1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,1)=result
      call dtwodq(fdv2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,2)=result          
      call dtwodq(fdv3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,3)=result    
      call dtwodq(fdv4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,4)=result 
        
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)     
      write(*,*)'52v',k,i           
      call dtwodq(fdv1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,1)=result
      call dtwodq(fdv2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,2)=result
      call dtwodq(fdv3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,3)=result
      call dtwodq(fdv4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)cccv
      write(3)c5v
      close (3)  
      write(*,*) 'end of graupel - water drop coag'
      return
      end
      real*8 function fd1(x,y)
      implicit none 
      real*8 x,y
      real*8 vx,rx,vy,ry,rey,rhogr,ec,pi
      real*8 kk,sege1,renh,ev,ea,tk,visc
      real*8 rx1,ry1,dedy,dedry,dedrx,ddy1,ddy2,ddx1,ddx2,dedx,ey1,     &
     &       ey2,e10
      real*8 e0(0:18,0:43)
      integer jump, k1,k2
      common/collhdgrwd/e0      ! collision effcienci  given by Khain et al.
!	     
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419      
!
!           
      jump=2     ! collision effcienci  given by Khain et al.;jump=1  collision efficiency given by Langmuir
      rhogr=800.
      tk=273.15+20.
      call tvwd (x,vx,rx)
      call tvgr (y,rhogr,vy,ry,rey)
      pi=asin(1.d0)*2.d0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      ec=0.0
      if (jump.eq.1) then
       sege1=2000.d0/9.d0*rx**2.d0
       renh=rey/60.d0              ! Nre/60.
       kk=sege1*vy/(visc*ry)             
       if (kk.gt.1.214d0) then
        ev=(1.d0+0.75*dlog(2.d0*kk)/(kk-1.214d0))**(-2.d0)
       else
        ev=0.d0
       endif
       if (kk.gt.0.0833d0) then           
        ea=(kk/(kk+0.5d0))**2.d0
       else
        ea=0.d0
       endif
       ec=(ev+ea*renh)/(1.d0+renh)
       if (ec.gt.1.) ec=1.0
       fd1=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif
      if (jump.eq.2) then
       rx1=rx*1e6                     ! water drop radius in um
       ry1=ry*1e6                     ! graupel radius in um 
       if (rx1.ge.260.) then
	if (ry1.ge.5.00) then
         ec=1.d0
        else
	 dedy=0.25
	 ec=dedy*(ry1-1.0)
        endif	   
        goto 40
       endif
       if (rx1.lt.1.0) then
	   ec=0.0
	   goto 40
       endif
       if (ry1.ge.320.0) then
	dedry=0.01/20.
	k1=0
55      k1=k1+1
        if(rx1.gt.e0(0,k1)) goto 55
	dedrx=(e0(18,k1)-e0(18,k1-1))/(e0(0,k1)-e0(0,k1-1))
	e10=dedrx*(rx1-e0(0,k1-1))+e0(18,k1-1)
	ec=(ry1-e0(18,0))*dedry+e10
	if (ec.gt.1.0) ec=1.0
	goto 40
       endif
       k1=0
58     k1=k1+1
       if (rx1.gt.e0(0,k1)) goto 58
       k2=0
59     k2=k2+1
       if (ry1.gt.e0(k2,0)) goto 59  
       ddy1=(ry1-e0(k2-1,0))
       ddy2=(e0(k2,0)-ry1)
       ddx1=(rx1-e0(0,k1-1))
       ddx2=(e0(0,k1)-rx1)
       ey1=(e0(k2-1,k1-1)*ddy2+e0(k2,k1-1)*ddy1)/(ddy1+ddy2)
       ey2=(e0(k2-1,k1)*ddy2+e0(k2,k1)*ddy1)/(ddy1+ddy2)
       dedx=(ey2-ey1)/(ddx1+ddx2)
       ec=dedx*ddx1+ey1
40     if (ec.lt.0.0) ec=0.0
       fd1=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif     ! endif for jump=2
      return
      end      
      function fd2(x,y)
      implicit real*8 (a-h,o-z)
      fd2=y*fd1(x,y)
      return
      end
      function fd3(x,y)
      implicit real*8 (a-h,o-z)
      fd3=x*fd1(x,y)
      return
      end
      function fd4(x,y)
      implicit real*8 (a-h,o-z)
      fd4=y*x*fd1(x,y)
      return
      end
      real*8 function fdv1(y,x) 
      implicit none
      real*8 x,y
      real*8 vx,rx,vy,ry,rey,rhogr,ec,pi
      real*8 kk,sege1,renh,ev,ea,tk,visc
      real*8 rx1,ry1,dedy,dedry,dedrx,ddy1,ddy2,ddx1,ddx2,dedx,ey1,     &
     &       ey2,e10
      integer k1,k2,jump
      real*8 e0(0:18,0:43)
      common /collhdgrwd/e0
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419      
      jump=2     ! collision effcienci  given by Khain et al.;jump=1  collision efficiency given by Langmuir
      rhogr=800.d0
      call tvwd (x,vx,rx)
      call tvgr (y,rhogr,vy,ry,rey)
      pi=asin(1.d0)*2.d0
      ec=0.0
      tk=273.15+20.
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      if (jump.eq.1) then
       sege1=2000.d0/9.d0*rx**2.d0
       renh=rey/60.d0              ! Nre/60
       kk=sege1*vy/(visc*ry)             
       if (kk.gt.1.214d0) then
        ev=(1.d0+0.75*dlog(2.d0*kk)/(kk-1.214d0))**(-2.d0)
       else
        ev=0.d0
       endif
       if (kk.gt.0.0833d0) then           
        ea=(kk/(kk+0.5d0))**2.d0
       else
        ea=0.d0
       endif
       ec=(ev+ea*renh)/(1.d0+renh)
       if (ec.gt.1.) ec=1.0
       fdv1=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif
      if (jump.eq.2) then
       rx1=rx*1e6                     ! water drop radius in um
       ry1=ry*1e6                     ! graupel radius in um 
       if (rx1.ge.260.) then
        if (ry1.ge.5.00) then
         ec=1.d0
        else
	 dedy=0.25
	 ec=dedy*(ry1-1.0)
        endif	   
        goto 40
       endif
       if (rx1.lt.1.0) then
	ec=0.0
	goto 40
       endif
       if (ry1.ge.320.0) then
	dedry=0.01/20.
	k1=0
55      k1=k1+1
        if(rx1.gt.e0(0,k1)) goto 55
	dedrx=(e0(12,k1)-e0(12,k1-1))/(e0(0,k1)-e0(0,k1-1))
	e10=dedrx*(rx1-e0(0,k1-1))+e0(12,k1-1)
	ec=(ry1-e0(12,0))*dedry+e10
	if (ec.gt.1.0) ec=1.0
	goto 40
       endif
       k1=0
58     k1=k1+1
       if (rx1.gt.e0(0,k1)) goto 58
       k2=0
59     k2=k2+1
       if (ry1.gt.e0(k2,0)) goto 59  
       ddy1=(ry1-e0(k2-1,0))
       ddy2=(e0(k2,0)-ry1)
       ddx1=(rx1-e0(0,k1-1))
       ddx2=(e0(0,k1)-rx1)
       ey1=(e0(k2-1,k1-1)*ddy2+e0(k2,k1-1)*ddy1)/(ddy1+ddy2)
       ey2=(e0(k2-1,k1)*ddy2+e0(k2,k1)*ddy1)/(ddy1+ddy2)
       dedx=(ey2-ey1)/(ddx1+ddx2)
       ec=dedx*ddx1+ey1	
40     continue
       fdv1=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif    ! endif for jump=2
      return
      end      
      function fdv2(x,y)
      implicit real*8 (a-h,o-z)
      fdv2=y*fdv1(x,y)
      return
      end
      function fdv3(x,y)
      implicit real*8 (a-h,o-z)
      fdv3=x*fdv1(x,y)
      return
      end
      function fdv4(x,y)
      implicit real*8 (a-h,o-z)
      fdv4=y*x*fdv1(x,y)
      return
      end

!
!  this subroutine calculates the kernel of graupel - water drop coalescence if the density of the graupel depends on its mass
!
      subroutine gencoef21(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable.
!
      real*8 cc(36,4),ccc(36,36,12)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external fd11,fd12,fd13,fd14,fdv11,fdv12,fdv13,fdv14
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,4
       do k=1,k4
        cc(k,i)=0.d0
       enddo
      enddo
      do i=1,12
       do j=1,k4
        do k=1,k4
          ccc(k,j,i)=0.d0
        enddo
       enddo
      enddo
      do i=1,4
       do j=1,k4
        do k=1,k4
         cccv(k,j,i)=0.d0
         c5v(k,j,i)=0.d0
        enddo
       enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2
      bot=m(i)
      top=m(i+1)
      write(*,*)'12',k,i
      call dtwodq(fd11,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(fd12,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result
      call dtwodq(fd13,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result
      call dtwodq(fd14,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'22',k,i
      call dtwodq(fd11,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(fd12,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(fd13,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(fd14,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)
      write(*,*)'32',k,i
      call dtwodq(fd11,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(fd12,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(fd13,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(fd14,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)
      call dtwodq(fd11,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,1)=result
      call dtwodq(fd12,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,2)=result
      call dtwodq(fd13,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,3)=result
      call dtwodq(fd14,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,4)=result
      enddo
!                 the water drops are larger than the graupel particles
      do k=3,k4
      do i=1,k-2
      bot=m(i)
      top=m(i+1)
      write(*,*)'12v',k,i
      call dtwodq(fdv11,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,1)=result
      call dtwodq(fdv12,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,2)=result
      call dtwodq(fdv13,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,3)=result
      call dtwodq(fdv14,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,4)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'52v',k,i
      call dtwodq(fdv11,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,1)=result
      call dtwodq(fdv12,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,2)=result
      call dtwodq(fdv13,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,3)=result
      call dtwodq(fdv14,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)cccv
      write(3)c5v
      close (3)
      write(*,*) 'end of graupel (w. changing density) - water drop'
      return
      end
!----+-----------------------------------------------------------------+
      real*8 function fd11(x,y)
      implicit none
      real*8 x,y
      real*8 vx,rx,vy,ry,rey,rhogr,ec,pi
      real*8 kk,sege1,renh,ev,ea,tk,visc
      real*8 rx1,ry1,dedy,dedry,dedrx,ddy1,ddy2,ddx1,ddx2,dedx,ey1,     &
     &       ey2,e10
      integer k1,k2,jump
      real*8 e0(0:18,0:43)
      common /collgrwd/e0
! the collison coefficients for rhogr = 400 kg/m3 is used for graupel size less than 260 um, and
! the collison coefficients for rhogr = 800 kg/m3 is used for graupel size larger and equal than 260 um
      real*8 m450,m800
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
!
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419
!
!
      call tvwd (x,vx,rx)
      call tvgr (y,rhogr,vy,ry,rey)
      pi=asin(1.d0)*2.d0
      ec=0.0
      tk=273.15+20.
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      jump=2    ! collision effcienci  given by Khain1 et al.;jump=1  collision efficiency given by Langmuir
      if (jump.eq.1) then
       sege1=2000.d0/9.d0*rx**2.d0
       renh=rey/60.d0              ! Nre/60.
       kk=sege1*vy/(visc*ry)
       if (kk.gt.1.214d0) then
        ev=(1.d0+0.75*dlog(2.d0*kk)/(kk-1.214d0))**(-2.d0)
       else
        ev=0.d0
       endif
       if (kk.gt.0.0833d0) then
        ea=(kk/(kk+0.5d0))**2.d0
       else
        ea=0.d0
       endif
       ec=(ev+ea*rey)/(1.d0+rey)
       if (ec.gt.1.) ec=1.0
       fd11=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif
      if (jump.eq.2) then
       rx1=rx*1e6                     ! water drop radius in um
       ry1=ry*1e6                     ! graupel radius in um
       if (rx1.ge.260.) then
        if (ry1.ge.5.00) then
         ec=1.d0
        else
         dedy=0.25
         ec=dedy*(ry1-1.0)
        endif
        goto 40
       endif
       if (rx1.lt.1.0) then
        ec=0.0
        goto 40
       endif
       if (ry1.ge.320.0) then
        dedry=0.01/20.
        k1=0
55      k1=k1+1
        if(rx1.gt.e0(0,k1)) goto 55
        dedrx=(e0(18,k1)-e0(18,k1-1))/(e0(0,k1)-e0(0,k1-1))
        e10=dedrx*(rx1-e0(0,k1-1))+e0(18,k1-1)
        ec=(ry1-e0(18,0))*dedry+e10
        if (ec.gt.1.0) ec=1.0
        goto 40
       endif
       k1=0
58     k1=k1+1
       if (rx1.gt.e0(0,k1)) goto 58
       k2=0
59     k2=k2+1
       if (ry1.gt.e0(k2,0)) goto 59
       ddy1=(ry1-e0(k2-1,0))
       ddy2=(e0(k2,0)-ry1)
       ddx1=(rx1-e0(0,k1-1))
       ddx2=(e0(0,k1)-rx1)
       ey1=(e0(k2-1,k1-1)*ddy2+e0(k2,k1-1)*ddy1)/(ddy1+ddy2)
       ey2=(e0(k2-1,k1)*ddy2+e0(k2,k1)*ddy1)/(ddy1+ddy2)
       dedx=(ey2-ey1)/(ddx1+ddx2)
       ec=dedx*ddx1+ey1
40     if (ec.lt.0.0) ec=0.0
       fd11=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif     ! endif for jump=2
      return
      end
      function fd12(x,y)
      implicit real*8 (a-h,o-z)
      fd12=y*fd11(x,y)
      return
      end
      function fd13(x,y)
      implicit real*8 (a-h,o-z)
      fd13=x*fd11(x,y)
      return
      end
      function fd14(x,y)
      implicit real*8 (a-h,o-z)
      fd14=y*x*fd11(x,y)
      return
      end
      real*8 function fdv11(y,x)
      implicit none
      real*8 x,y
      real*8 vx,rx,vy,ry,rey,rhogr,ec,pi
      real*8 kk,sege1,renh,ev,ea,tk,visc
      real*8 rx1,ry1,dedy,dedry,dedrx,ddy1,ddy2,ddx1,ddx2,dedx,ey1,     &
     &       ey2,e10
      integer k1,k2,jump
      real*8 e0(0:18,0:43)
      common /collgrwd/ e0

! the collison coefficients for rhogr = 400 kg/m3 is used for graupel size less than 260 um, and
! the collison coefficients for rhogr = 800 kg/m3 is used for graupel size larger and equal than 260 um
      real*8 m450,m800
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419
      jump=2     ! collision effcienci  given by Khain et al.;jump=1  collision efficiency given by Langmuir
      call tvwd (x,vx,rx)
      call tvgr (y,rhogr,vy,ry,rey)
      pi=asin(1.d0)*2.d0
      ec=0.0
      tk=273.15+20.
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      if (jump.eq.1) then
       sege1=2000.d0/9.d0*rx**2.d0
       renh=rey/60.d0              ! Nre/60.
       kk=sege1*vy/(visc*ry)
       if (kk.gt.1.214d0) then
        ev=(1.d0+0.75*dlog(2.d0*kk)/(kk-1.214d0))**(-2.d0)
       else
        ev=0.d0
       endif
       if (kk.gt.0.0833d0) then
        ea=(kk/(kk+0.5d0))**2.d0
       else
        ea=0.d0
       endif
       ec=(ev+ea*renh)/(1.d0+renh)
       if (ec.gt.1.) ec=1.0
       fdv11=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif
      if (jump.eq.2) then
       rx1=rx*1e6                     ! water drop radius in um
       ry1=ry*1e6                     ! graupel radius in um
       if (rx1.ge.260.) then
        if (ry1.ge.5.00) then
         ec=1.d0
        else
         dedy=0.25
         ec=dedy*(ry1-1.0)
        endif
        goto 40
       endif
       if (rx1.lt.1.0) then
        ec=0.0
        goto 40
       endif
       if (ry1.ge.320.0) then
        dedry=0.01/20.
        k1=0
55      k1=k1+1
        if(rx1.gt.e0(0,k1)) goto 55
        dedrx=(e0(18,k1)-e0(18,k1-1))/(e0(0,k1)-e0(0,k1-1))
        e10=dedrx*(rx1-e0(0,k1-1))+e0(18,k1-1)
        ec=(ry1-e0(18,0))*dedry+e10
        if (ec.gt.1.0) ec=1.0
        goto 40
       endif
       k1=0
58     k1=k1+1
       if (rx1.gt.e0(0,k1)) goto 58
       k2=0
59     k2=k2+1
       if (ry1.gt.e0(k2,0)) goto 59
       ddy1=(ry1-e0(k2-1,0))
       ddy2=(e0(k2,0)-ry1)
       ddx1=(rx1-e0(0,k1-1))
       ddx2=(e0(0,k1)-rx1)
       ey1=(e0(k2-1,k1-1)*ddy2+e0(k2,k1-1)*ddy1)/(ddy1+ddy2)
       ey2=(e0(k2-1,k1)*ddy2+e0(k2,k1)*ddy1)/(ddy1+ddy2)
       dedx=(ey2-ey1)/(ddx1+ddx2)
       ec=dedx*ddx1+ey1
40     continue
       fdv11=ec*pi*(rx+ry)**2.d0*dabs(vx-vy)
      endif    ! endif for jump=2
      return
      end
      function fdv12(x,y)
      implicit real*8 (a-h,o-z)
      fdv12=y*fdv11(x,y)
      return
      end
      function fdv13(x,y)
      implicit real*8 (a-h,o-z)
      fdv13=x*fdv11(x,y)
      return
      end
      function fdv14(x,y)
      implicit real*8 (a-h,o-z)
      fdv14=y*x*fdv11(x,y)
      return
      end



!
!   this subroutine calculates the kernel for pristine ice crystal - water drop collision
!                       
      subroutine gencoef4(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!    c5(k,i,1)..c5(k,i,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner integral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,8),c5(36,36,4)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40) 
!      character filen*40
      common/m/ m
      common/k/ k
      external fcid1,fcid2,fcid3,fcid4,fdci1,fdci2,fdci3,fdci4
      external g1,g2,g3
      external h1,h2 
!
! cleaning all matrixes
!
      do i=1,4 
       do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo
      do i=1,8        
       do k=1,k4
        do j=1,k4
         ccc(k,j,i)=0.d0
        enddo
       enddo 
      enddo
      do i=1,4
       do k=1,k4 
        do j=1,k4
         cccv(k,j,i)=0.d0 
         c5v(k,j,i)=0.d0
        enddo
       enddo   
      enddo
!
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'13',k,i 
      call dtwodq(fcid1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(fcid2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result        
      call dtwodq(fcid3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result        
      call dtwodq(fcid4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result  
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)     
      write(*,*)'23',k,i    
      call dtwodq(fcid1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(fcid2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(fcid3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(fcid4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)      
      call dtwodq(fcid1,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,1)=result
      call dtwodq(fcid2,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,2)=result 
      call dtwodq(fcid3,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,3)=result       
      call dtwodq(fcid4,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,4)=result
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)     
      write(*,*)'53',k,i           
      call dtwodq(fcid1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,1)=result
      call dtwodq(fcid2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,2)=result
      call dtwodq(fcid3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,3)=result
      call dtwodq(fcid4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,4)=result
      enddo
      enddo
!                 the mass of water drops is larger than that of ice crystals   
      do k=3,k4 
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'13v',k,i
      call dtwodq(fdci1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,1)=result
      call dtwodq(fdci2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,2)=result         
      call dtwodq(fdci3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,3)=result       
      call dtwodq(fdci4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,4)=result    
      enddo
      enddo
      do k=2,k4 
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'53v',k                   
      call dtwodq(fdci1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,1)=result
      call dtwodq(fdci2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,2)=result
      call dtwodq(fdci3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,3)=result
      call dtwodq(fdci4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)c5
      write(3)cccv
      write(3)c5v
      close (3)
      write(*,*) 'end of  p. ice  - water drop col.'
      return
      end
      real*8 function fcid1(x,y)
      implicit none 
      real*8 x,y
      real*8 vx,rx,vy,dy
      real*8 e,dym,rxm,e1,e2,pi      
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419      
  
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvpice(y,vy,dy)                 
      if (rx.lt.5.d-6.or.dy.lt.1.d-4) then
       fcid1=0.d0
       return
      endif
      if (rx.gt.51.6666d-6.and.dy.gt.500.d-6) then
       e=0.5
       goto 30
      endif            
      dym=dy*1.d6
      rxm=rx*1.d6
      if (dym.le.160.d0) then   
       e1=0.d0
       e2=rxm*(-0.0175d0*rxm+0.5775d0)-4.41d0
       if(e2.lt.0.d0) e2=0.d0
       e=(e2-e1)/60.d0*(dym-100.d0)
       goto 30
      endif
      if(dym.gt.160.d0.and.dym.le.194.d0) then
       e1=rxm*(-0.0175d0*rxm+0.5775d0)-4.41d0
       e2=rxm*(-0.005d0*rxm+0.205d0)-1.55d0 
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0
       e=(e2-e1)/(194.d0-160.d0)*(dym-160.d0)+e1
       goto 30
      endif
      if(dym.gt.194.d0.and.dym.lt.213.d0) then
       e1=rxm*(-0.005d0*rxm+0.205d0)-1.55d0
       e2=rxm*(-0.0038d0*rxm+0.1692d0)-1.21d0 
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0       
       e=(e2-e1)/(213.d0-194.d0)*(dym-194.d0)+e1
       goto 30    
      endif
      if(dym.gt.213.d0.and.dym.lt.289.d0) then
       e1=rxm*(-0.0038d0*rxm+0.1692d0)-1.21d0
       e2=rxm*(-0.0019d0*rxm+0.1047d0)-0.64d0 
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0      
       e=(e2-e1)/(289.d0-213.d0)*(dym-213.d0)+e1
       goto 30    
      endif
      if(dym.gt.289.d0.and.dym.lt.404.d0) then
       e1=rxm*(-0.0019d0*rxm+0.1047d0)-0.64d0 
       e2=rxm*(-0.0012d0*rxm+0.080d0)-0.43d0 
       if (rxm.gt.51.6667d0) e2=0.5d0
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0      
       e=(e2-e1)/(404.d0-289.d0)*(dym-289.d0)+e1
       goto 30    
      endif
      if (dym.gt.404.d0.and.dym.le.500.d0) then
       e1=rxm*(-0.0012d0*rxm+0.080d0)-0.43d0 
       if (rxm.gt.51.666d0) e1=0.5d0 
       if (rxm.gt.10.d0) then
        e2=0.5d0
       else
        e2=0.5d0/5.d0*(rxm-5.d0) 
       endif  
       if(e1.lt.0.d0) e1=0.d0
       e=(e2-e1)/(500.d0-404.d0)*(dym-404.d0)+e1
      endif 
      if (dym.gt.500.d0) then
       if (rxm.gt.10.d0) then
        e=0.5d0
        goto 30
       endif
       e=0.5d0/5.d0*(rxm-5.d0) 
      endif  
30    fcid1=e*pi*(rx+dy/2.d0)**2.d0*dabs(vx-vy)        
      return
      end      
      function fcid2(x,y)
      implicit real*8 (a-h,o-z)
      fcid2=y*fcid1(x,y)
      return
      end
      function fcid3(x,y)
      implicit real*8 (a-h,o-z)
      fcid3=x*fcid1(x,y)
      return
      end
      function fcid4(x,y)
      implicit real*8 (a-h,o-z)
      fcid4=y*x*fcid1(x,y)
      return
      end
      real*8 function fdci1(y,x)
      implicit none
      real*8 x,y
      real*8 vx,rx,vy,dy
      real*8 e,dym,rxm,e1,e2,pi   
!	   
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419      
!
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvpice(y,vy,dy)                 
      if (rx.lt.5.d-6.or.dy.lt.1d-4) then
       fdci1=0.d0
       return
      endif	
      if (rx.gt.51.6666d-6.and.dy.gt.500.d-6) then
       e=0.5
       goto 30
      endif            
      dym=dy*1.d6
      rxm=rx*1.d6
      if (dym.le.160.d0) then   
       e1=0.d0
       e2=rxm*(-0.0175d0*rxm+0.5775d0)-4.41d0
       if(e2.lt.0.d0) e2=0.d0
       e=(e2-e1)/60.d0*(dym-100.d0)
       goto 30
      endif
      if(dym.gt.160.d0.and.dym.le.194.d0) then
       e1=rxm*(-0.0175d0*rxm+0.5775d0)-4.41d0
       e2=rxm*(-0.005d0*rxm+0.205d0)-1.55d0 
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0
       e=(e2-e1)/(194.d0-160.d0)*(dym-160.d0)+e1
       goto 30
      endif
      if(dym.gt.194.d0.and.dym.lt.213.d0) then
       e1=rxm*(-0.005d0*rxm+0.205d0)-1.55d0
       e2=rxm*(-0.0038d0*rxm+0.1692d0)-1.21d0 
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0       
       e=(e2-e1)/(213.d0-194.d0)*(dym-194.d0)+e1
       goto 30    
      endif
      if(dym.gt.213.d0.and.dym.lt.289.d0) then
       e1=rxm*(-0.0038d0*rxm+0.1692d0)-1.21d0
       e2=rxm*(-0.0019d0*rxm+0.1047d0)-0.64d0 
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0      
       e=(e2-e1)/(289.d0-213.d0)*(dym-213.d0)+e1
       goto 30    
      endif
      if(dym.gt.289.d0.and.dym.lt.404.d0) then
       e1=rxm*(-0.0019d0*rxm+0.1047d0)-0.64d0 
       e2=rxm*(-0.0012d0*rxm+0.080d0)-0.43d0   
       if (rxm.gt.51.6666d0) e2=0.5d0
       if(e1.lt.0.d0) e1=0.d0
       if(e2.lt.0.d0) e2=0.d0      
       e=(e2-e1)/(404.d0-289.d0)*(dym-289.d0)+e1
       goto 30    
      endif
      if (dym.gt.404.d0.and.dym.lt.500.d0) then
       e1=rxm*(-0.0012d0*rxm+0.080d0)-0.43d0   
       if(rxm.gt.51.6666d0)  e1=0.5d0
       if (rxm.gt.10.d0) then
        e2=0.5d0
       else
        e2=0.5d0/5.d0*(rxm-5.d0) 
       endif  
       if(e1.lt.0.d0) e1=0.d0
       e=(e2-e1)/(500.d0-404.d0)*(dym-404.d0)+e1
      endif
      if (dym.gt.500.d0) then
       if (rxm.gt.10.d0) then
        e=0.5d0
        goto 30
       endif
       e=0.5d0/5.d0*(rxm-5.d0) 
      endif      
30    fdci1=e*pi*(rx+dy/2.d0)**2.d0*dabs(vx-vy)        
      return
      end      
      function fdci2(x,y)
      implicit real*8 (a-h,o-z)
      fdci2=y*fdci1(x,y)
      return
      end
      function fdci3(x,y)
      implicit real*8 (a-h,o-z)
      fdci3=x*fdci1(x,y)
      return
      end
      function fdci4(x,y)
      implicit real*8 (a-h,o-z)
      fdci4=y*x*fdci1(x,y)
      return
      end   
! This routine calculates the collection kernel of the rimed ice - crystal water drop collision     
      subroutine gencoef5(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!    c5(k,i,1)..c5(k,i,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,12)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external fri1,fri2,fri3,fri4,friv1,friv2,friv3,friv4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,4
        do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo 
      do i=1,12
        do j=1,k4
          do k=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo     
      do i=1,4
        do j=1,k4
          do k=1,k4   
          cccv(j,k,i)=0.d0
          c5v(j,k,i)=0.d0
          enddo
        enddo
      enddo
!
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'14',k,i 
      call dtwodq(fri1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(fri2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result        
      call dtwodq(fri3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result        
      call dtwodq(fri4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result  
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'24',k,i      
      call dtwodq(fri1,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(fri2,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(fri3,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(fri4,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)  
      write(*,*)'34',k,i  
      call dtwodq(fri1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(fri2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(fri3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(fri4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(fri1,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,1)=result
      call dtwodq(fri2,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,2)=result   
      call dtwodq(fri3,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,3)=result      
      call dtwodq(fri4,bot,top,g3,h1,errabs,errrel,irule,result,errest)
      cc(k,4)=result
      enddo
!                 the water drops are larger than the rimed ice crystals
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'14v',k,i 
      call dtwodq(friv1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,1)=result
      call dtwodq(friv2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,2)=result          
      call dtwodq(friv3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,3)=result    
      call dtwodq(friv4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      cccv(k,i,4)=result   
      enddo
      enddo
      do k=2,k4 
      do i=1,k-1
      bot=m(i)
      top=m(i+1)     
      write(*,*)'54v',k,i  
      call dtwodq(friv1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,1)=result
      call dtwodq(friv2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,2)=result
      call dtwodq(friv3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,3)=result
      call dtwodq(friv4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)cccv
      write(3)c5v
      close (3) 
      write(*,*) 'end of snow flake water coll.'
      return
      end
      real*8 function fri1(x,y)
      implicit none
      real*8 x,y
      real*8 vx,vy,rx,dy,rhosnw 
      real*8 pi,y1,y0,y2,x1,x0,x2,a,b,e,ta,bv,yv
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419            
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvsnow(y,vy,dy,rhosnw) 
      if (rx.lt.5.0d-6) then
       fri1=0.d0
       return
      endif 
      y1=1.d-4
      y0=5.d-4
      y2=1.d-5
      x1=5.d-6
      x0=5.d-5
      x2=5.d-4
      if(rx.ge.x0.and.dy.ge.y0) then
       e=1.d0
       goto 30 
      endif 
      if (rx.lt.x0.and.dy.lt.y0) then
       a=-1.d0/(x1-x0)**2.d0
       b=-1.d0/(y1-y0)**2.d0
       e=a*(rx-x0)**2.d0+b*(dy-y0)**2.d0+1.d0
       goto 30
      endif  
      if (rx.lt.x0) then
       a=-1.d0/(x1-x0)**2.d0
       e=a*(rx-x0)**2.d0+1.d0
       goto 30
      endif
      if (rx.gt.x0.and.rx.le.x2) then
       ta=(y1-y2)/(x0-x2)
       bv=y1-ta*x0
       yv=ta*rx+bv
       b=-1.d0/(yv-y0)**2.d0
       e=b*(dy-y0)**2.d0+1.d0
       goto 30
      endif
      b=-1.d0/(y2-y0)**2.d0
      e=b*(dy-y0)**2.d0+1.d0 
30    if (e.lt.0.d0) e=0.d0 
      if (e.gt.1.d0) e=1.d0
       fri1=e*(rx+dy/2.0d0)**2.d0*pi*dabs(vx-vy)      
      return
      end      
      function fri2(x,y)
      implicit real*8 (a-h,o-z)
      fri2=y*fri1(x,y)
      return
      end
      function fri3(x,y)
      implicit real*8 (a-h,o-z)
      fri3=x*fri1(x,y)
      return
      end
      function fri4(x,y)
      implicit real*8 (a-h,o-z)
      fri4=y*x*fri1(x,y)
      return
      end
      real*8 function friv1(y,x)
      implicit none
      real*8 x,y
      real*8 vx,vy,rx,dy,rhosnw
      real*8 pi,y1,y0,y2,x1,x0,x2,a,b,e,ta,bv,yv
      pi=asin(1.d0)*2.d0
      call tvwd(x,vx,rx)
      call tvsnow(y,vy,dy,rhosnw)  
      if (rx.le.5.d-6) then
       friv1=0.d0
       return
      endif
      y1=1.d-4
      y0=5.d-4
      y2=1.d-5
      x1=5.d-6
      x0=5.d-5
      x2=5.d-4
      if(rx.ge.x0.and.dy.ge.y0) then
       e=1.d0
       goto 30 
      endif 
      if (rx.lt.x0.and.dy.lt.y0) then
       a=-1.d0/(x1-x0)**2.d0
       b=-1.d0/(y1-y0)**2.d0
       e=a*(rx-x0)**2.d0+b*(dy-y0)**2.d0+1.d0
       goto 30
      endif  
      if (rx.lt.x0) then
       a=-1.d0/(x1-x0)**2.d0
       e=a*(rx-x0)**2.d0+1.d0
       goto 30
      endif
      if (rx.gt.x0.and.rx.le.x2) then
       ta=(y1-y2)/(x0-x2)
       bv=y1-ta*x0
       yv=ta*rx+bv
       b=-1.d0/(yv-y0)**2.d0
       e=b*(dy-y0)**2.d0+1.d0       
       goto 30
      endif
      b=-1.d0/(y2-y0)**2.d0
      e=b*(dy-y0)**2.d0+1.d0 
30    if (e.lt.0.d0) e=0.d0 
      if (e.gt.1.d0) e=1.d0
       friv1=e*(rx+dy/2.0d0)**2.d0*pi*dabs(vx-vy)      
      return
      end      
      function friv2(x,y)
      implicit real*8 (a-h,o-z)
      friv2=y*friv1(x,y)
      return
      end
      function friv3(x,y)
      implicit real*8 (a-h,o-z)
      friv3=x*friv1(x,y)
      return
      end
      function friv4(x,y)
      implicit real*8 (a-h,o-z)
      friv4=y*x*friv1(x,y)
      return
      end

!    calculation of kernel of ice crystal - ice crystal  collision, collision efficiency 
!   is  equal to unity if the crysatls are larger than 50 um, is zero if the one of he crystals less
!    20 um, in other cases the collison efficiency changes linearly between 0 and 1, depending on the
!    crystal size. 
!    The coalesence efficiency depends on the temperature
!          
      subroutine gencoef6(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,8)
!    cc(k,1)...cc(k,4)
!    c5(k,i,1)..c5(k,i,4)
!  instead of ccc(k,i,52)... ccc(k,i,55)
!  we take    c5(k,i,1).... c5(k,i,4)  to save memory space
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,8),c5(36,36,4)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external ficic1,ficic2,ficic3,ficic4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,4
       do k=1,k4
        cc(k,i)=0.d0
       enddo
      enddo
      do i=1,8
       do j=1,k4
        do k=1,k4
          ccc(k,j,i)=0.d0
        enddo
       enddo
      enddo
      do i=1,4
        do j=1,k4
          do k=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
!
      errabs=0.d0
      errrel=5.d-4 
      errest=0.d0
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1) 
      write(*,*) '1',k,i
      call dtwodq(ficic1,bot,top,g1,h1,errabs,errrel,irule,result,      &          
     &             errest)
      ccc(k,i,1)=result
      call dtwodq(ficic2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result     
      call dtwodq(ficic3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result     
      call dtwodq(ficic4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1) 
      write(*,*) '3',k,i
      call dtwodq(ficic1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(ficic2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result 
      call dtwodq(ficic3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result     
      call dtwodq(ficic4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k) 
      write(*,*) '41',k      
      call dtwodq(ficic1,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest) 
      cc(k,1)=result
      call dtwodq(ficic2,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest)
      cc(k,2)=result
      call dtwodq(ficic3,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest)
      cc(k,3)=result
      call dtwodq(ficic4,bot,top,g3,h1,errabs,errrel,irule,result,       &
     &            errest)
      cc(k,4)=result
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1) 
      write(*,*) '4',k,i            
      call dtwodq(ficic1,bot,top,h1,g2,errabs,errrel,irule,result,        &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(ficic2,bot,top,h1,g2,errabs,errrel,irule,result,        &
     &            errest)
      c5(k,i,2)=result 
      call dtwodq(ficic3,bot,top,h1,g2,errabs,errrel,irule,result,        &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(ficic4,bot,top,h1,g2,errabs,errrel,irule,result,        &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)c5
      close (3) 
      write(*,*) 'end of p. ice - p. ice aggregation '
      return
      end

      real*8 function ficic1(x,y)
      implicit none 
      real*8 x,y 
      real*8 vx,vy,dx,dy
      real*8 dx1,dx2,dy1,dy2,d,tgalfa,xv,yv,ec,pi
      pi=asin(1.d0)*2.d0
      call tvpice (x,vx,dx)
      call tvpice (y,vy,dy) 
      dx1=2.d-5
      dx2=1.d-4    
      dy1=2.d-5
      dy2=1.d-4        
      ec=0.d0              
      if (dabs(dx-dy).lt.dy*1.d-12) then
       d=dsqrt((dx2-dx1)**2.d0+(dy2-dy1)**2.d0)
       ec=dsqrt((dx1-dx)**2.d0+(dy1-dy)**2.d0)/d 
       vy=vx
       goto 15
      endif
      if (dx.le.dx1.or.dy.le.dy1) then
       ficic1=0.d0
       return
      endif
      if (dx.ge.dx2.and.dy.ge.dy2) then
       ec=1.d0
       goto 15
      endif     
      if (dy.ge.dy2.and.dx.le.dx2) then
       ec=1.d0/(dx2-dx1)*(dx-dx1)  
       goto 15
      endif     
      if (dx.ge.dx2.and.dy.le.dy2) then
       ec=1.d0/(dy2-dy1)*(dy-dy1)  
       goto 15
      endif  
      if (dx.gt.dy) then
       tgalfa=(dx2-dx)/(dy2-dy) 
       xv=dx2-tgalfa*(dy2-dy1)
       yv=dy1
       ec=(yv*dx-xv*dy)/(yv*dx2-xv*dy2)
       goto 15      
      endif 
      if(dx.lt.dy) then 
       tgalfa=(dy2-dy)/(dx2-dx)
       xv=dx1 
       yv=dy2-tgalfa*(dx2-dx1)
       ec=(yv*dx-xv*dy)/(yv*dx2-xv*dy2)
      endif
15    ficic1=pi*(dx+dy)**2.d0/4.d0*ec*dabs(vx-vy)
      return
      end
      function ficic2(x,y)
      implicit real*8 (a-h,o-z)
      ficic2=y*ficic1(x,y)
      return
      end
      function ficic3(x,y)
      implicit real*8 (a-h,o-z)
      ficic3=x*ficic1(x,y)
      return
      end
      function ficic4(x,y)
      implicit real*8 (a-h,o-z)
      ficic4=y*x*ficic1(x,y)
      return
      end
!
!    calculation of kernel of snow flake crystal - snow flake crystal  collision, collision efficiency 
!   is supposed to be equal to unity, the coalesence efficiency depends on the temperature
!          
      subroutine gencoef7(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,8)
!    cc(k,1)...cc(k,4)
!    c5(k,i,1)..c5(k,i,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,3),ccc(36,36,12)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external fricric1,fricric2,fricric3,fricric4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,3
        do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo
      do i=1,12
        do j=1,k4
          do k=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4 
      errest=0.d0
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1) 
      write(*,*) '1',k,i
      call dtwodq(fricric1,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fricric2,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,2)=result     
      call dtwodq(fricric3,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,3)=result     
      call dtwodq(fricric4,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1) 
      write(*,*) '3',k,i
      call dtwodq(fricric1,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fricric2,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,6)=result 
      call dtwodq(fricric3,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,7)=result     
      call dtwodq(fricric4,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1) 
      write(*,*) '3',k,i
      call dtwodq(fricric1,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fricric2,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,10)=result 
      call dtwodq(fricric3,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,11)=result     
      call dtwodq(fricric4,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)
      write(*,*) '41',k      
      call dtwodq(fricric1,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,1)=result
      call dtwodq(fricric2,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,2)=result
      call dtwodq(fricric4,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,3)=result
      enddo
      write(3)cc
      write(3)ccc
      close (3) 
      write(*,*) 'end of snow flake -snow flake '
      return
      end
!
      real*8 function fricric1(x,y)
      implicit none
      real*8 x,y
      real*8 vx,vy,dx,dy,pi,ec,rhosnw 
      pi=asin(1.d0)*2.d0
      call tvsnow(x,vx,dx,rhosnw)
      call tvsnow(y,vy,dy,rhosnw)
      if(dabs(dx-dy).lt.dx*1.d-12) then
       fricric1=0.d0
       return
      endif
      ec=1.0
      fricric1=pi*(dx+dy)**2.d0/4.d0*ec*dabs(vx-vy)
      return
      end
      function fricric2(x,y)
      implicit real*8 (a-h,o-z)
      fricric2=y*fricric1(x,y)
      return
      end
      function fricric3(x,y)
      implicit real*8 (a-h,o-z)
      fricric3=x*fricric1(x,y)
      return
      end
      function fricric4(x,y)
      implicit real*8 (a-h,o-z)
      fricric4=y*x*fricric1(x,y)
      return
      end
!
!   this subroutine calculates the kernel for snow flakes  - pristine ice collision
!                       
      subroutine gencoef8(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!    c5(k,i,1)..c5(k,i,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner integral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,12)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40) 
!      character filen*40
      common/m/ m
      common/k/ k
      external fsnowic1,fsnowic2,fsnowic3,fsnowic4,ficsnow1,            &
     &   ficsnow2,ficsnow3,ficsnow4
      external g1,g2,g3
      external h1,h2 
!
!  cleaning all matrixes
!
      do i=1,4 
       do k=1,k4
        cc(k,i)=0.d0
        enddo
      enddo
      do i=1,12        
        do k=1,k4
          do j=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo 
      enddo
      do i=1,4
        do k=1,k4 
          do j=1,k4
           cccv(k,j,i)=0.d0 
           c5v(k,j,i)=0.d0
          enddo
        enddo   
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      call dtwodq(fsnowic1,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fsnowic2,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,2)=result        
      call dtwodq(fsnowic3,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,3)=result        
      call dtwodq(fsnowic4,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,4)=result  
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'22',k,i      
      call dtwodq(fsnowic1,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fsnowic2,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fsnowic3,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fsnowic4,bot,top,g2,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)     
      write(*,*)'32',k,i  
      call dtwodq(fsnowic1,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fsnowic2,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fsnowic3,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fsnowic4,bot,top,h1,h2,errabs,errrel,irule,result,    &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(fsnowic1,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,1)=result
      call dtwodq(fsnowic2,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,2)=result   
      call dtwodq(fsnowic3,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,3)=result      
      call dtwodq(fsnowic4,bot,top,g3,h1,errabs,errrel,irule,result,    &
     &            errest)
      cc(k,4)=result
      enddo
!                 the pristine ice crystals are larger than the snow flakes      
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'12v',k,i 
      call dtwodq(ficsnow1,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      cccv(k,i,1)=result
      call dtwodq(ficsnow2,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      cccv(k,i,2)=result          
      call dtwodq(ficsnow3,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      cccv(k,i,3)=result    
      call dtwodq(ficsnow4,bot,top,g1,h1,errabs,errrel,irule,result,    &
     &            errest)
      cccv(k,i,4)=result   
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)     
      write(*,*)'52v',k,i           
      call dtwodq(ficsnow1,bot,top,h1,g2,errabs,errrel,irule,result,    &
     &            errest)
      c5v(k,i,1)=result
      call dtwodq(ficsnow2,bot,top,h1,g2,errabs,errrel,irule,result,    &
     &            errest)
      c5v(k,i,2)=result
      call dtwodq(ficsnow3,bot,top,h1,g2,errabs,errrel,irule,result,    &
     &            errest)
      c5v(k,i,3)=result
      call dtwodq(ficsnow4,bot,top,h1,g2,errabs,errrel,irule,result,    &
     &            errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)cccv
      write(3)c5v
      close (3)
      write(*,*) 'end of snow flakes - pristine ice '
      return  
      end
      real*8 function fsnowic1(x,y)
      implicit real*8 (a-h,o-z)
      real*8 msnow100,msnow500,msnowtv,rhosnw 
      pi=asin(1.d0)*2.d0
      call tvpice (x,vx,dx)
      call tvsnow (y,vy,dy,rhosnw)
      e=1.d0
      fsnowic1=e*pi*(dx/2.d0+dy/2.d0)**2.d0*dabs(vx-vy)        
      return
      end      
      function fsnowic2(x,y)
      implicit real*8 (a-h,o-z)
      fsnowic2=y*fsnowic1(x,y)
      return
      end
      function fsnowic3(x,y)
      implicit real*8 (a-h,o-z)
      fsnowic3=x*fsnowic1(x,y)
      return
      end
      function fsnowic4(x,y)
      implicit real*8 (a-h,o-z)
      fsnowic4=y*x*fsnowic1(x,y)
      return
      end
      real*8 function ficsnow1(y,x)
      implicit real*8 (a-h,o-z)
      real*8 msnow100,msnow500,msnowtv,rhosnw
      pi=asin(1.d0)*2.d0
      call tvpice (x,vx,dx)
      call tvsnow (y,vy,dy,rhosnw)
      e=1.d0
      ficsnow1=e*pi*(dx/2.d0+dy/2.d0)**2.d0*dabs(vx-vy)        
      return
      end      
      function ficsnow2(x,y)
      implicit real*8 (a-h,o-z)
      ficsnow2=y*ficsnow1(x,y)
      return
      end
      function ficsnow3(x,y)
      implicit real*8 (a-h,o-z)
      ficsnow3=x*ficsnow1(x,y)
      return
      end
      function ficsnow4(x,y)
      implicit real*8 (a-h,o-z)
      ficsnow4=y*x*ficsnow1(x,y)
      return
      end  
!
!   This routine calculates the collection kernel of the pristine ice -  high density graupel particle
!
      subroutine gencoef9(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,12)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external fgrci1,fgrci2,fgrci3,fgrci4,fcigr1,fcigr2,fcigr3,fcigr4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,4
       do k=1,k4
        cc(k,i)=0.d0
       enddo
      enddo
      do i=1,12
        do j=1,k4
          do k=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
       do j=1,k4
        do k=1,k4
         cccv(k,j,i)=0.d0 
         c5v(k,j,i)=0.d0
        enddo
       enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      call dtwodq(fgrci1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fgrci2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result        
      call dtwodq(fgrci3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result        
      call dtwodq(fgrci4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result  
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'22',k,i      
      call dtwodq(fgrci1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fgrci2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fgrci3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fgrci4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)     
      write(*,*)'32',k,i  
      call dtwodq(fgrci1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fgrci2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fgrci3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fgrci4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(fgrci1,bot,top,g3,h1,errabs,errrel,irule,result,      &
     &            errest)
      cc(k,1)=result
      call dtwodq(fgrci2,bot,top,g3,h1,errabs,errrel,irule,result,      &
     &            errest)
      cc(k,2)=result   
      call dtwodq(fgrci3,bot,top,g3,h1,errabs,errrel,irule,result,      &
     &            errest)
      cc(k,3)=result      
      call dtwodq(fgrci4,bot,top,g3,h1,errabs,errrel,irule,result,      &
     &            errest)
      cc(k,4)=result
      enddo
!                 the pristine ice crystals are larger than the graupel particles      
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'12v',k,i 
      call dtwodq(fcigr1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      cccv(k,i,1)=result
      call dtwodq(fcigr2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      cccv(k,i,2)=result          
      call dtwodq(fcigr3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      cccv(k,i,3)=result    
      call dtwodq(fcigr4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      cccv(k,i,4)=result   
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)     
      write(*,*)'52v',k,i           
      call dtwodq(fcigr1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5v(k,i,1)=result
      call dtwodq(fcigr2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5v(k,i,2)=result
      call dtwodq(fcigr3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5v(k,i,3)=result
      call dtwodq(fcigr4,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)cccv
      write(3)c5v
      close (3)  
      write(*,*) 'end of   graupel - pristine ice'
      return
      end
      real*8 function fgrci1(x,y)
      implicit none
      real*8 x,y 
      real*8 vx,dx,vy,ry,rhogr,rey,pi,e
      real*8 m450,m800
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      pi=asin(1.d0)*2.d0
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      call tvpice(x,vx,dx)
      call tvgr(y,rhogr,vy,ry,rey)
      e=1.0
      fgrci1=e*pi*(dx/2.d0+ry)**2.d0*dabs(vx-vy)
      return
      end      
      function fgrci2(x,y)
      implicit real*8 (a-h,o-z)
      fgrci2=y*fgrci1(x,y)
      return
      end
      function fgrci3(x,y)
      implicit real*8 (a-h,o-z)
      fgrci3=x*fgrci1(x,y)
      return
      end
      function fgrci4(x,y)
      implicit real*8 (a-h,o-z)
      fgrci4=y*x*fgrci1(x,y)
      return
      end
      real*8 function fcigr1(y,x)
      implicit none
      real*8 x,y
      real*8 vx,dx,vy,ry,rhogr,rey,pi,e
      real*8 m450,m800
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      pi=asin(1.d0)*2.d0
      call tvpice(x,vx,dx)
      call tvgr(y,rhogr,vy,ry,rey)
      e=1.0
      fcigr1=e*pi*(dx/2.d0+ry)**2.d0*dabs(vx-vy)
      return
      end      
      function fcigr2(x,y)
      implicit real*8 (a-h,o-z)
      fcigr2=y*fcigr1(x,y)
      return
      end
      function fcigr3(x,y)
      implicit real*8 (a-h,o-z)
      fcigr3=x*fcigr1(x,y)
      return
      end
      function fcigr4(x,y)
      implicit real*8 (a-h,o-z)
      fcigr4=y*x*fcigr1(x,y)
      return
      end	 
!
!   This routine calculates the collection kernel of the snow flakes -  high density graupel particle
!
      subroutine gencoef10(k4)
      implicit real*8 (a-h,o-z)
!
!  calculates the analytical expressions for the following coef. :
!    ccc(k,i,1)...ccc(k,i,12)
!    cc(k,1)...cc(k,4)
!
!  This routine uses my routine dtwodq to perform the double
!  integration by Simpson's rule. The inner inegral bounderies are in many cases function
!  of the outer integration variable. 
!
      real*8 cc(36,4),ccc(36,36,12)
      real*8 cccv(36,36,4),c5v(36,36,4)
      real*8 m(40)
      common/m/ m
      common/k/ k
      external fgrrim1,fgrrim2,fgrrim3,fgrrim4,frimgr1,frimgr2,         &
     &	frimgr3,frimgr4
      external g1,g2,g3
      external h1,h2
!
!  cleaning all matrixes
!
      do i=1,4
       do k=1,k4
        cc(k,i)=0.d0
       enddo
      enddo
      do i=1,12
       do j=1,k4
        do k=1,k4
          ccc(k,j,i)=0.d0
        enddo
       enddo
      enddo
      do i=1,4
       do j=1,k4
        do k=1,k4
         cccv(k,j,i)=0.d0 
         c5v(k,j,i)=0.d0
        enddo
       enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      call dtwodq(fgrrim1,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fgrrim2,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,2)=result                                                 
      call dtwodq(fgrrim3,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,3)=result        
      call dtwodq(fgrrim4,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,4)=result  
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)
      write(*,*)'22',k,i      
      call dtwodq(fgrrim1,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fgrrim2,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fgrrim3,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fgrrim4,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m(i)
      top=m(i+1)     
      write(*,*)'32',k,i  
      call dtwodq(fgrrim1,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fgrrim2,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fgrrim3,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fgrrim4,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=2,k4
      bot=m(k-1)
      top=m(k)   
      call dtwodq(fgrrim1,bot,top,g3,h1,errabs,errrel,irule,result,     &
     &            errest)
      cc(k,1)=result
      call dtwodq(fgrrim2,bot,top,g3,h1,errabs,errrel,irule,result,     &
     &            errest)
      cc(k,2)=result   
      call dtwodq(fgrrim3,bot,top,g3,h1,errabs,errrel,irule,result,     &
     &            errest)
      cc(k,3)=result      
      call dtwodq(fgrrim4,bot,top,g3,h1,errabs,errrel,irule,result,     &
     &            errest)
      cc(k,4)=result
      enddo
!                 the snowflakes  are larger than the graupel particles      
      do k=3,k4
      do i=1,k-2 
      bot=m(i)
      top=m(i+1)
      write(*,*)'12v',k,i 
      call dtwodq(frimgr1,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      cccv(k,i,1)=result
      call dtwodq(frimgr2,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      cccv(k,i,2)=result          
      call dtwodq(frimgr3,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      cccv(k,i,3)=result    
      call dtwodq(frimgr4,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      cccv(k,i,4)=result   
      enddo
      enddo
      do k=2,k4
      do i=1,k-1
      bot=m(i)
      top=m(i+1)     
      write(*,*)'52v',k,i           
      call dtwodq(frimgr1,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5v(k,i,1)=result
      call dtwodq(frimgr2,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5v(k,i,2)=result
      call dtwodq(frimgr3,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5v(k,i,3)=result
      call dtwodq(frimgr4,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5v(k,i,4)=result
      enddo
      enddo
      write(3)cc
      write(3)ccc
      write(3)cccv
      write(3)c5v
      close (3) 
      write(*,*) 'end of  graupel - rimed ice coll.'
      return
      end
      real*8 function fgrrim1(x,y)
      implicit none
      real*8 x,y
      real*8 vx,dx,vy,ry,rhogr,rhosnw,rey,pi,e
      real*8 m450,m800
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      pi=asin(1.d0)*2.d0
      call tvsnow(x,vx,dx,rhosnw)
      call tvgr (y,rhogr,vy,ry,rey)
      e=1.0
      fgrrim1=e*pi*(dx/2.d0+ry)**2.d0*dabs(vx-vy)
      return
      end      
      function fgrrim2(x,y)
      implicit real*8 (a-h,o-z)
      fgrrim2=y*fgrrim1(x,y)
      return
      end
      function fgrrim3(x,y)
      implicit real*8 (a-h,o-z)
      fgrrim3=x*fgrrim1(x,y)
      return
      end
      function fgrrim4(x,y)
      implicit real*8 (a-h,o-z)
      fgrrim4=y*x*fgrrim1(x,y)
      return
      end
      real*8 function frimgr1(y,x)
      implicit none
      real*8 x,y
      real*8 vx,dx,vy,ry,rhogr,rhosnw,rey,pi,e
      real*8 m450,m800
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      pi=asin(1.d0)*2.d0
      call tvsnow(x,vx,dx,rhosnw)
      call tvgr (y,rhogr,vy,ry,rey)      
      e=1.0
      frimgr1=e*pi*(dx/2.d0+ry)**2.d0*dabs(vx-vy)
      return
      end      
      function frimgr2(x,y)
      implicit real*8 (a-h,o-z)
      frimgr2=y*frimgr1(x,y)
      return
      end
      function frimgr3(x,y)
      implicit real*8 (a-h,o-z)
      frimgr3=x*frimgr1(x,y)
      return
      end
      function frimgr4(x,y)
      implicit real*8 (a-h,o-z)
      frimgr4=y*x*frimgr1(x,y)
      return
      end	 
!
!     This subroutine calculates the terminal velocities of the water drops,
!       pristine ice crystals, rimed ice crystals, aggregates and graupel.
!       The terminal velocities depends on the particle mass at the midle of
!       bins.
!
      subroutine terminalv(k4)
      implicit none
      real*8 m(40),mmean,vx,rx,dx,rhogr,m450,m800,rey,rhosnw
      real*4 vc(36),vic(36),vric(36),vgra(36),vgra1(36)
      real*4 rd,dic,dsnow,rgr,rgr1
      integer k,k4	 
      common /m/ m              
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
      open(unit=29,file='termvel1.txt')
 90   FORMAT('K',6x,'mass',12x,'rd',6x,'vtwd',12x,'dic',6x,'vtpice',    &
     & 12x,'dsnow',6x,'vtsnow',12x,'rgr',6x,'vtgra',12x,'rgr1',6x,      &
     & 'vtgra1')
 100  FORMAT(I3,11(2x,e10.3))
      write(29,90)   
      do k=1,k4
       mmean=1.5*m(k)
       call tvwd(mmean,vx,rx)
       vc(k) = vx
       rd=rx
       call tvpice(mmean,vx,dx)
       vic(k) =vx
       dic=dx
       call tvsnow(mmean,vx,dx,rhosnw)
       vric(k) = vx
       dsnow =dx
       rhogr=800.0
       call tvgr(mmean,rhogr,vx,rx,rey)
       vgra(k) = vx
       rgr = rx
       if (mmean.lt.m450) then
        rhogr = 450.0
       endif
       if (mmean.gt.m800) then
        rhogr = 800.0
       endif
       if (mmean.gt.m450.and.mmean.lt.m800) then
        rhogr=(800.0 - 450.0)/(m800 - m450)*(mmean - m450)+ 450.0
       endif
       call tvgr(mmean,rhogr,vx,rx,rey)
       vgra1(k)= vx
       rgr1 = rx
       write(29,100) k,mmean,rd,vc(k),dic,vic(k),dsnow,vric(k),         &
     &             rgr,vgra(k),rgr1,vgra1(k)
      enddo 
      write(3) vc,vic,vric,vgra,vgra1
      return
      end
!
!   This subroutine calculates the kernel for the ice multiplication
!    ice fragments are produced due to the high density graupel - drop collision
!      
      subroutine icemulp(k4)
      implicit real*8 (a-h,o-z)
      real*8 c(36,36,2)
!  This routine uses qsimp to perform the single integration by Simpson's rule. 
!               (Numerical recipies routine)
!
      real*8 m(40),mmean
      common/m/ m
      common/k/ k
      external fm1,fm2
!
!  cleaning  matrix
!
      do i=1,k4
       do k=1,k4
	do l=1,2
         c(k,i,l)=0.d0
	enddo
       enddo
      enddo
      do i=1,k4
       do k=1,k4
       bot=m(k)
       top=m(k+1)
       mmean=1.5d0*m(i)
       if(i.eq.8.and.k.eq.1) then
	 xx=0.0
       endif
       write(*,*) i,k
       call qsimp(fm1,bot,top,mmean,result)
       c(k,i,1)=result
       call qsimp(fm2,bot,top,mmean,result)
       c(k,i,2)=result
       enddo
      enddo
      write(*,*) 'end of icemultipl.'
      write(3) c
      return
      end
      real*8 function fm1(x,y)
      implicit none
      real*8 x,y
      real*8 vx,rx,vy,ry,rhogr,rey
      real*8 sege1,renh,kk,tk,visc
      real*8 ev,ea,ec
      real*8 rx1,ry1,dedx,dedrx,dedy,dedry,ddy1,ddx1,ddy2,ddx2,e10,     &
     &       ey1,ey2
      real*8 e0(0:18,0:43) 
      common /collhdgrwd/ e0
      integer jump,k1,k2
!	     
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419      
!
!           
      
      jump=2     ! collision effcienci  given by Khain et al.;jump=1  collision efficiency given by Langmuir
      tk=273.15+20.
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      rhogr=800.0
      call tvwd(x,vx,rx)
      call tvgr(y,rhogr,vy,ry,rey)
!        collection effciency given by Langmuir (1948)
      if (jump.eq.1) then
       sege1=2000.d0/9.d0*rx**2.d0
       renh=rey/60.d0              ! Nre/60.
       kk=sege1*vy/(visc*ry)             
       if (kk.gt.1.214d0) then
        ev=(1.d0+0.75*dlog(2.d0*kk)/(kk-1.214d0))**(-2.d0)
       else
        ev=0.d0
       endif
       if (kk.gt.0.0833d0) then           
        ea=(kk/(kk+0.5d0))**2.d0
       else
        ea=0.d0
       endif
       ec=(ev+ea*renh)/(1.d0+renh)
       if (ec.gt.1.) ec=1.0
       fm1=ec
      endif
      if (jump.eq.2) then
       rx1=rx*1e6                     ! water drop radius in um
       ry1=ry*1e6                     ! graupel radius in um 
       if (rx1.ge.260.) then
        if (ry1.ge.5.00) then
         ec=1.d0
        else
	 dedy=0.25
         ec=dedy*(ry1-1.0)
        endif	   
        goto 40
       endif
       if (rx1.lt.1.0) then
        ec=0.0
        goto 40
       endif
       if (ry1.ge.320.0) then
        dedry=0.01/20.
        k1=0
55      k1=k1+1
        if(rx1.gt.e0(0,k1)) goto 55
        dedrx=(e0(18,k1)-e0(18,k1-1))/(e0(0,k1)-e0(0,k1-1))
        e10=dedrx*(rx1-e0(0,k1-1))+e0(18,k1-1)
        ec=(ry1-e0(18,0))*dedry+e10
        if (ec.gt.1.0) ec=1.0
        goto 40
       endif
       k1=0
58     k1=k1+1
       if (rx1.gt.e0(0,k1)) goto 58
       k2=0
59     k2=k2+1
       if (ry1.gt.e0(k2,0)) goto 59  
!          write(*,*) p,k2,k1,rn,rk
       ddy1=(ry1-e0(k2-1,0))
       ddy2=(e0(k2,0)-ry1)
       ddx1=(rx1-e0(0,k1-1))
       ddx2=(e0(0,k1)-rx1)
       ey1=(e0(k2-1,k1-1)*ddy2+e0(k2,k1-1)*ddy1)/(ddy1+ddy2)
       ey2=(e0(k2-1,k1)*ddy2+e0(k2,k1)*ddy1)/(ddy1+ddy2)
       dedx=(ey2-ey1)/(ddx1+ddx2)
       ec=dedx*ddx1+ey1
40     if (ec.lt.0.0) ec=0.0
       fm1=ec*x**(2.d0/3.d0)
      endif
      return
      end      
      function fm2(x,y)
      implicit real*8 (a-h,o-z)
      fm2=x*fm1(x,y)
      return
      end            			     		
!
!   This subroutine calculates the kernel for the ice multiplication
!    ice fragments are produced due to the variable density graupel - drop collision
!
      subroutine icemulp1(k4)
      implicit real*8 (a-h,o-z)
      real*8 c(36,36,2)
!  This routine uses qsimp to perform the single integration by Simpson's rule.
!               (Numerical recipies routine)
!
      real*8 m(40),mmean
      common/m/ m
      common/k/ k
      external fm11,fm12
!
!  cleaning  matrix
!
      do i=1,k4
       do k=1,k4
        do l=1,2
         c(k,i,l)=0.d0
        enddo
       enddo
      enddo
      do i=1,k4
       do k=1,k4
       bot=m(k)
       top=m(k+1)
       mmean=1.5d0*m(i)
       if(i.eq.8.and.k.eq.1) then
         xx=0.0
       endif
       write(*,*) i,k
       call qsimp(fm11,bot,top,mmean,result)
       c(k,i,1)=result
       call qsimp(fm12,bot,top,mmean,result)
       c(k,i,2)=result
       enddo
      enddo
      write(*,*) 'end of icemultipl.'
      write(3) c
      return
      end
      real*8 function fm11(x,y)
      implicit none
      real*8 x,y
      real*8 vx,rx,vy,ry,rhogr,rey
      real*8 sege1,renh,kk,tk,visc,m450,m800
      real*8 ev,ea,ec
      real*8 rx1,ry1,dedx,dedrx,dedy,dedry,ddy1,ddx1,ddy2,ddx2,e10,     &
     &       ey1,ey2
      real*8 e0(0:18,0:43)
      common /collgrwd/ e0
      integer jump,k1,k2
      parameter ( m450=4.888790e-9,m800=4.289321169e-6)
!
! The terminal velocity of the drops are given by the book of Pruppacher and Klett p417-419
!
!

      jump=2     ! collision effcienci  given by Khain et al.;jump=1  collision efficiency given by Langmuir
      tk=273.15+20.
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      call tvwd(x,vx,rx)
      call tvgr(y,rhogr,vy,ry,rey)
!        collection effciency given by Langmuir (1948)
      if (jump.eq.1) then
       sege1=2000.d0/9.d0*rx**2.d0
       renh=rey/60.d0              ! Nre/60.
       kk=sege1*vy/(visc*ry)
       if (kk.gt.1.214d0) then
        ev=(1.d0+0.75*dlog(2.d0*kk)/(kk-1.214d0))**(-2.d0)
       else
        ev=0.d0
       endif
       if (kk.gt.0.0833d0) then
        ea=(kk/(kk+0.5d0))**2.d0
       else
        ea=0.d0
       endif
       ec=(ev+ea*renh)/(1.d0+renh)
       if (ec.gt.1.) ec=1.0
       fm11=ec
      endif
      if (jump.eq.2) then
       rx1=rx*1e6                     ! water drop radius in um^M
       ry1=ry*1e6                     ! graupel radius in um ^M
       if (rx1.ge.260.) then
        if (ry1.ge.5.00) then
         ec=1.d0
        else
         dedy=0.25
         ec=dedy*(ry1-1.0)
        endif
        goto 40
       endif
       if (rx1.lt.1.0) then
        ec=0.0
        goto 40
       endif
       if (ry1.ge.320.0) then
        dedry=0.01/20.
        k1=0
55      k1=k1+1
        if(rx1.gt.e0(0,k1)) goto 55
        dedrx=(e0(18,k1)-e0(18,k1-1))/(e0(0,k1)-e0(0,k1-1))
        e10=dedrx*(rx1-e0(0,k1-1))+e0(18,k1-1)
        ec=(ry1-e0(18,0))*dedry+e10
        if (ec.gt.1.0) ec=1.0
        goto 40
       endif
       k1=0
58     k1=k1+1
       if (rx1.gt.e0(0,k1)) goto 58
       k2=0
59     k2=k2+1
       if (ry1.gt.e0(k2,0)) goto 59
!          write(*,*) p,k2,k1,rn,rk
       ddy1=(ry1-e0(k2-1,0))
       ddy2=(e0(k2,0)-ry1)
       ddx1=(rx1-e0(0,k1-1))
       ddx2=(e0(0,k1)-rx1)
       ey1=(e0(k2-1,k1-1)*ddy2+e0(k2,k1-1)*ddy1)/(ddy1+ddy2)
       ey2=(e0(k2-1,k1)*ddy2+e0(k2,k1)*ddy1)/(ddy1+ddy2)
       dedx=(ey2-ey1)/(ddx1+ddx2)
       ec=dedx*ddx1+ey1
40     if (ec.lt.0.0) ec=0.0
       fm11=ec*x**(2.d0/3.d0)
      endif
      return
      end
      function fm12(x,y)
      implicit real*8 (a-h,o-z)
      fm12=x*fm11(x,y)
      return
      end
  
!
!   This routine calculates the collection kernel of the water drop - aerosol particles gravitational collison
!        (it is supposed that the size of the drops is much lager than that of the aerosol particles)
!
      subroutine gencoef1(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,36,12),c5(36,36,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external ff1,ff2,ff3,ff4,ffv1,ffv2,ffv3,ffv4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      imax=min(36,k+17)       !  17=k8-k4-2
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 grav ',k,i
      call dtwodq(ff1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(ff2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result
      call dtwodq(ff3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result
      call dtwodq(ff4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      imax=min(36,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 grav',k,i
      call dtwodq(ff1,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(ff2,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(ff3,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(ff4,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 grav',k,i
      call dtwodq(ff1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(ff2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(ff3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(ff4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      imax=min(k4,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(ff1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,1)=result
      call dtwodq(ff2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,2)=result
      call dtwodq(ff3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,3)=result
      call dtwodq(ff4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of the water drop -aerosol praticle grav. coag.'
      return
      end
      real*8 function ff1(x,y)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 RD,RAP,V
      REAL*8 nre,k1,k0,GS,CSC,Z,H,YC0,ECOL,SNRE,Y1,RY1,PSEG
      REAL*8 B0,B1,B2,B3,A0,A1,A2
      REAL*8 RG,RA,KB,RHOW,RHOAER,RHOA,PI,TK,VISC,SIGMA,PRES,LAMB
      parameter(b0=0.1465d0,b1=1.302d0,b2=-0.607d0,b3=0.293d0,          &
     &  rg=287.0,RA=1.87e-10,KB=1.38e-23)
      parameter(a0=-.1007,a1=-0.358,a2=0.0261)
      COMMON /RHOAER/ RHOAER
      
      CALL tvwd(y,v,rd)
      if (rd.lt.4.d-5) then         ! the Ecol=0.0 in this case
       ff1=0.d0
       return
      endif
      pi=dasin(1.d0)*2.d0
      RAP=(3.d0*x/(4.d0*PI*RHOAER))**(1.d0/3.d0)
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      sigma=(75.93-0.115*(tk-273.15))*1e-3
      RHOA=pres/tk/rg
      lamb=KB*TK/(4.d0*1.414*pi*ra*ra*pres)
      snre=visc/2./RHOA             
      nre =v*rd/snre
      ecol=0.d0
      y1=dlog(nre)
      gs=a0+a1*y1+a2*y1*y1
      k0=dexp(gs)
      csc=1.+1.26*lamb/RAP
      pseg=RAP/rd
      if(pseg.gt.1.) then
        ecol=0.0
        goto 40
      endif
      k1=pseg*pseg*RHOAER*nre*csc/9./RHOA
      z=dlog(k1/k0)
      h=b0+b1*z+b2*z*z+b3*z*z*z
      yc0=2./pi*datan(h)
      if (yc0+pseg.gt.0.0) then
        ecol=((yc0+pseg)/(1.+pseg))**2
      else
        ecol=0.0
      endif
40    ff1=pi*rd*rd*ecol*v
      return
      end
      function ff2(x,y)
      implicit real*8 (a-h,o-z)
      ff2=y*ff1(x,y)
      return
      end
      function ff3(x,y)
      implicit real*8 (a-h,o-z)
      ff3=x*ff1(x,y)
      return
      end
      function ff4(x,y)
      implicit real*8 (a-h,o-z)
      ff4=x*y*ff1(x,y)
      return
      end
      real*8 function ffv1(y,x)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 RD,RAP,V
      REAL*8 lamb,nre,k1,k0,GS,CSC,Z,H,YC0,ECOL
      REAL*8 B0,B1,B2,B3,A0,A1,A2
      REAL*8 G,RG,RA,KB,RHOAER,RHOA,PI,TK,VISC,SIGMA,PRES
      parameter(b0=0.1465d0,b1=1.302d0,b2=-0.607d0,b3=0.293d0,rg=287.0, &
     &           RA=1.87e-10,KB=1.38e-23)
      parameter(a0=-.1007,a1=-0.358,a2=0.0261)
      REAL*8 SEG,SNRE,XV,X1,Y1,RY1,PSEG
      COMMON /RHOAER/ RHOAER
! All the parameters are given at temperature of 20 C and pressure of 1013mb
      CALL tvwd(y,v,rd)
      if (rd.lt.4.d-5) then         ! the Ecol=0.0 in this case
       ffv1=0.d0
       return
      endif
      pi=dasin(1.d0)*2.d0
      RAP=(3.d0*x/(4.d0*PI*RHOAER))**(1.d0/3.d0)
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      sigma=(75.93-0.115*(tk-273.15))*1e-3
      RHOA=pres/tk/rg
      lamb=KB*TK/(4.d0*1.414*pi*ra*ra*pres)
      snre=visc/2.d0/RHOA
      nre =v*rd/snre
      ecol=0.d0
20    y1=dlog(nre)
      gs=a0+a1*y1+a2*y1*y1
      k0=dexp(gs)
      csc=1.+1.26*lamb/RAP
      pseg=ra/rd
      if(pseg.gt.1.) then
        ecol=0.0
        goto 40
      endif
      k1=pseg*pseg*RHOAER*nre*csc/9./RHOA
      z=dlog(k1/k0)
      h=b0+b1*z+b2*z*z+b3*z*z*z
      yc0=2./pi*datan(h)
      if (yc0+pseg.gt.0.0) then
        ecol=((yc0+pseg)/(1.+pseg))**2
      else
        ecol=0.0
      endif
40    ffv1=pi*rd*rd*ecol*v
      return
      end
      function ffv2(x,y)
      implicit real*8 (a-h,o-z)
      ffv2=y*ffv1(x,y)
      return
      end
      function ffv3(x,y)
      implicit real*8 (a-h,o-z)
      ffv3=x*ffv1(x,y)
      return
      end
      function ffv4(x,y)
      implicit real*8 (a-h,o-z)
      ffv4=x*y*ffv1(x,y)
      return
      end
!
!   This subroutine calculates the kernel for tutbulenece generated collison
!
      subroutine genturbcoag(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,36,12),c5(36,36,4)
      real*8 m(40),m1(55)
      common/m/ m
      common/m1/ m1
      common/k/ k
      external ftr1,ftr2,ftr3,ftr4,ftrv1,ftrv2,ftrv3,ftrv4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
      do i=1,12
        do k=1,k4
          do j=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      imax=min(k4,k+17)       !  17=k8-k4-2
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 grav ',k,i
      call dtwodq(ftr1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(ftr2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result
      call dtwodq(ftr3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result
      call dtwodq(ftr4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      imax=min(36,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 grav',k,i
      call dtwodq(ftr1,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(ftr2,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(ftr3,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(ftr4,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 grav',k,i
      call dtwodq(ftr1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(ftr2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(ftr3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(ftr4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      imax=min(k4,k+18)        ! 18=k8-k4-1                       
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(ftr1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,1)=result
      call dtwodq(ftr2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,2)=result
      call dtwodq(ftr3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,3)=result
      call dtwodq(ftr4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of water aerosol scavenging by turbulence'
      return
      end
      real*8 function ftr1(x,y)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 RHOAER,RD,RAP,PI
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      RD=(3.d0*Y/(4.d0*3.14159*1000.d0))**(1.d0/3.d0)
      RAP=(3.d0*X/(4.d0*3.14159*RHOAER))**(1.d0/3.d0)
      FTR1=(2*RD+2*RAP)**3.d0
      return
      end
      function ftr2(x,y)
      implicit real*8 (a-h,o-z)
      ftr2=y*ftr1(x,y)
      return
      end
      function ftr3(x,y)
      implicit real*8 (a-h,o-z)
      ftr3=x*ftr1(x,y)
      return
      end
      function ftr4(x,y)
      implicit real*8 (a-h,o-z)
      ftr4=x*y*ftr1(x,y)
      return
      end
      real*8 function ftrv1(y,x)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 RHOAER,RD,RAP,PI
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      RD=(3.d0*Y/(4.d0*3.14159*1000.d0))**(1.d0/3.d0)
      RAP=(3.d0*X/(4.d0*3.14159*RHOAER))**(1.d0/3.d0)
      FTRV1=(2*RD+2*RAP)**3.d0
      return
      end
      function ftrv2(x,y)
      implicit real*8 (a-h,o-z)
      ftrv2=y*ftrv1(x,y)
      return
      end
      function ftrv3(x,y)
      implicit real*8 (a-h,o-z)
      ftrv3=x*ftrv1(x,y)
      return
      end
      function ftrv4(x,y)
      implicit real*8 (a-h,o-z)
      ftrv4=x*y*ftrv1(x,y)
      return
      end
!
!   This subroutine calculates the kernel for tutbulenece generated collison with p. ice
!
      subroutine genturbcoag_pice(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,36,12),c5(36,36,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external ftrpi1,ftrpi2,ftrpi3,ftrpi4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
      do i=1,12
        do k=1,k4
          do j=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      imax=min(k4,k+17)       !  17=k8-k4-2
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 grav ',k,i
      call dtwodq(ftrpi1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &     errest)
      ccc(k,i,1)=result
      call dtwodq(ftrpi2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(ftrpi3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(ftrpi4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      imax=min(36,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 grav',k,i
      call dtwodq(ftrpi1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(ftrpi2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(ftrpi3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(ftrpi4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 grav',k,i
      call dtwodq(ftrpi1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(ftrpi2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(ftrpi3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(ftrpi4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      imax=min(k4,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(ftrpi1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(ftrpi2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(ftrpi3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(ftrpi4,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*)'end of aerosol -  pristine ice due to the turb'
      return
      end
      real*8 function ftrpi1(x,y)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 SEG,RHOAER,DY,RAP,PI
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      DY =16.28d0*DSQRT(Y)
      RAP=(3.d0*X/(4.d0*PI*RHOAER))**(1.d0/3.d0)
      FTRPI1=(DY/PI+2*RAP)**3.d0                    ! for the p. ice particles the capacitance is used
      return
      end
      function ftrpi2(x,y)
      implicit real*8 (a-h,o-z)
      ftrpi2=y*ftrpi1(x,y)
      return
      end
      function ftrpi3(x,y)
      implicit real*8 (a-h,o-z)
      ftrpi3=x*ftrpi1(x,y)
      return
      end
      function ftrpi4(x,y)
      implicit real*8 (a-h,o-z)
      ftrpi4=x*y*ftrpi1(x,y)
      return
      end
!
!   This subroutine calculates the kernel for tutbulenece generated collison with snowflakes
!
      subroutine genturbcoag_snow(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,36,12),c5(36,36,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external ftrsnw1,ftrsnw2,ftrsnw3,ftrsnw4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
      do i=1,12
        do k=1,k4
          do j=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      imax=min(k4,k+17)       !  17=k8-k4-2
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 trb_snow ',k,i
      call dtwodq(ftrsnw1,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(ftrsnw2,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(ftrsnw3,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(ftrsnw4,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      imax=min(36,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 trb_snw',k,i
      call dtwodq(ftrsnw1,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(ftrsnw2,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(ftrsnw3,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(ftrsnw4,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 trb_snw',k,i
      call dtwodq(ftrsnw1,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(ftrsnw2,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(ftrsnw3,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(ftrsnw4,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      imax=min(k4,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53 trb_snw',k,i
      call dtwodq(ftrsnw1,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(ftrsnw2,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(ftrsnw3,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(ftrsnw4,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of snow - aerosol scav. by turbulence'
      return
      end
      real*8 function ftrsnw1(x,y)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 RHOAER,VY,DY,RAP,PI,axrmin,axrmax,msnow100,                   &
     &       msnow500,capac,a1,a2,b1,b2,axr,rhosnow,epsilon
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      call tvsnow (y,vy,dy,rhosnow)
      axrmin=0.08d0
      axrmax=1.0d0
      msnow100=3.773e-11
      msnow500=2.22529e-8
      if (y.le.msnow100) then
       axr = axrmin
      endif
      if (y.ge.msnow500) then
       axr = axrmax
      endif
      if (y.gt.msnow100.and.y.lt.msnow500) then
       axr=6.0d0*y/(dy**3.d0*PI*RHOSNOW)
      endif
      if (axr.gt.0.9999) then
       capac=0.5
      else
       EPSILON=DSQRT(1.d0-AXR*AXR)
       CAPAC=EPSILON/2.d0/ASIN(EPSILON)
      endif
      RAP=(3.d0*X/(4.d0*PI*RHOAER))**(1.d0/3.d0)
      FTRSNW1=(CAPAC*DY+2*RAP)**3.d0
      return
      end
      function ftrsnw2(x,y)
      implicit real*8 (a-h,o-z)
      ftrsnw2=y*ftrsnw1(x,y)
      return
      end
      function ftrsnw3(x,y)
      implicit real*8 (a-h,o-z)
      ftrsnw3=x*ftrsnw1(x,y)
      return
      end
      function ftrsnw4(x,y)
      implicit real*8 (a-h,o-z)
      ftrsnw4=x*y*ftrsnw1(x,y)
      return
      end
!
!   This subroutine calculates the kernel for tutbulenece generated collison with graupel
!
      subroutine genturbcoag_gr(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,36,12),c5(36,36,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external ftrgr1,ftrgr2,ftrgr3,ftrgr4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
      pk=2.d0
      do i=1,12
        do k=1,k4
          do j=1,k4
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,k4
          c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      imax=min(k4,k+17)       !  17=k8-k4-2
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 trbgr ',k,i
      call dtwodq(ftrgr1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(ftrgr2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(ftrgr3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(ftrgr4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      imax=min(36,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 trbgr',k,i
      call dtwodq(ftrgr1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(ftrgr2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(ftrgr3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(ftrgr4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,k4
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 trbgr',k,i
      call dtwodq(ftrgr1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(ftrgr2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(ftrgr3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(ftrgr4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      imax=min(k4,k+18)        ! 18=k8-k4-1
      do i=1,imax
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(ftrgr1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(ftrgr2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(ftrgr3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(ftrgr4,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of  graupel - aerosol turbulence scavanging'
      return
      end
      real*8 function ftrgr1(x,y)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 RHOAER,RD,RAP,PI,RHOGR,M450,M800
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      m450=4.888790e-9
      m800=4.289321169e-6
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      RD=(3.d0*X/(4.d0*PI*RHOGR))**(1.d0/3.d0)
      RAP=(3.d0*X/(4.d0*PI*RHOAER))**(1.d0/3.d0)
      FTRGR1=(2*RD+2*RAP)**3.d0
      return
      end
      function ftrgr2(x,y)
      implicit real*8 (a-h,o-z)
      ftrgr2=y*ftrgr1(x,y)
      return
      end
      function ftrgr3(x,y)
      implicit real*8 (a-h,o-z)
      ftrgr3=x*ftrgr1(x,y)
      return
      end
      function ftrgr4(x,y)
      implicit real*8 (a-h,o-z)
      ftrgr4=x*y*ftrgr1(x,y)
      return
      end
!   This routine calculates the collection kernel of the water drop - aerosol particles
!             due to the Brownian collison
!
      subroutine genbrown(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fbr1,fbr2,fbr3,fbr4,fbrv1,fbrv2,fbrv3,fbrv4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
!
      errabs=0.d0
      errrel=5.d-4
      irule=3
!      do k=3,k4
      do K=2,K4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 br',k,i
      call dtwodq(fbr1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(fbr2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result
      call dtwodq(fbr3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result
      call dtwodq(fbr4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 br',k,i
      call dtwodq(fbr1,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(fbr2,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(fbr3,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(fbr4,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 br',k,i
      call dtwodq(fbr1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(fbr2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(fbr3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(fbr4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,20
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(fbr1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,1)=result
      call dtwodq(fbr2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,2)=result
      call dtwodq(fbr3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,3)=result
      call dtwodq(fbr4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of water drops  - aerosol scav. by Brownian'
      return
      end
      real*8 function fbr1(x,y)
      IMPLICIT NONE
      real*8 X,Y
      REAL*8 V,RD,RAP
      REAL*8 lamb,nre,kn,nsch,al,bl,cl,alfa
      REAL*8 PI,RG,RHOAER,KB,RA,TK,PRES,VISC,RHOA,DIFF,FV
      parameter(RG=287.0,KB=1.38e-23,RA=1.87e-10)
      parameter(al=1.257,bl=0.4,cl=1.1)  !parameters for alfa  in calc. of  diff
      REAL*8 SEG,SEG1,SNRE,XV,X1,Y1,RY1
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      call tvwd (y,v,rd)
      rap=(3.d0*x/(4.d0*3.14159*RHOAER))**(1.d0/3.d0)
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5     
      RHOA=pres/tk/rg
      lamb=KB*TK/(4.d0*1.414*pi*ra*ra*pres)            ! free path of air moleculas
      KN=LAMB/RAP
      alfa = al + bl*dexp(-cl/kn)
      diff=kb*tk*(1.+alfa*kn)/6./rap/visc/pi
!      nsch=visc/diff/rhoa
      nsch=0.63
      snre=visc/2./RHOA
      NRE =v*rd/snre
      seg1=nsch**0.33333*dsqrt(nre)
      if (seg1.lt.1.4) then
        fv = 1. + 0.108*seg1**2
      else
        fv = 0.78 + 0.308*seg1
      endif
      fbr1=4*pi*rd*diff*fv
!      if (rd.gt.10.0e-6) then
!         write(99,*) 'rp,rd,diff,nsch,fv,v,nre,fbr1'
!         write(99,*) rp,rd,diff,nsch,fv,v,nre,fbr1
!        stop
!      endif
!      fbr1=(y/x)**(1.d0/3.d0)
      return
      end
      function fbr2(x,y)
      implicit real*8 (a-h,o-z)
      fbr2=y*fbr1(x,y)
      return
      end
      function fbr3(x,y)
      implicit real*8 (a-h,o-z)
      fbr3=x*fbr1(x,y)
      return
      end
      function fbr4(x,y)
      implicit real*8 (a-h,o-z)
      fbr4=x*y*fbr1(x,y)
      return
      end
      real*8 function fbrv1(y,x)
      implicit real*8 (a-h,o-z)
      fbrv1=(y/x)**(1.d0/3.d0)
      return
      end
      function fbrv2(x,y)
      implicit real*8 (a-h,o-z)
      fbrv2=y*fbrv1(x,y)
      return
      end
      function fbrv3(x,y)
      implicit real*8 (a-h,o-z)
      fbrv3=x*fbrv1(x,y)
      return
      end
      function fbrv4(x,y)
      implicit real*8 (a-h,o-z)
      fbrv4=x*y*fbrv1(x,y)
      return
      end
!   This routine calculates the collection kernel of the pristine ice - aerosol particles
!             due to the Brownian collison
!
      subroutine genbrown_pice(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fbrpi1,fbrpi2,fbrpi3,fbrpi4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do K=2,K4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 brpi',k,i
      call dtwodq(fbrpi1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fbrpi2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(fbrpi3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(fbrpi4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 brpi',k,i
      call dtwodq(fbrpi1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fbrpi2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fbrpi3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fbrpi4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 brpi',k,i
      call dtwodq(fbrpi1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fbrpi2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fbrpi3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fbrpi4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,20
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53 brpi',k,i
      call dtwodq(fbrpi1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fbrpi2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fbrpi3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fbrpi4,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of  p ice - aerosol scav. by Brownian'
      return
      end
      real*8 function fbrpi1(x,y)
      IMPLICIT NONE
      real*8 X,Y
      REAL*8 VY,DY,DY1,RAP,lamb,kn,nsch,al,bl,cl,alfa
      REAL*8 PI,RG,RHOAER,KB,RA,TK,PRES,VISC,RHOA,DIFF,FV,HCI0,RHOI,MV         
      parameter(RG=287.0,KB=1.38e-23,MV=1.075d-13,                      &
     &         RA=1.87e-10,HCI0=5.33776d-6,RHOI=900.0)
      parameter(al=1.257,bl=0.4,cl=1.1)  !parameters for alfa  in calc. of  diff
! All the parameters are given at temperature of 20 C and pressure of 1013mb
      REAL*8 CORR,HCI,RAXI,EPSILON,CAPAC,REYCI,KHI
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      CALL TVPICE (y,vy,dy1)
      RAP=(3.d0*x/(4.d0*3.14159*RHOAER))**(1.d0/3.d0)
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      RHOA=pres/tk/rg
!      lamb=8.89919d-8
      lamb=KB*TK/(4.d0*1.414*pi*ra*ra*pres)            ! free path of air moleculas
      KN=LAMB/RAP
      alfa = al + bl*dexp(-cl/kn)
      diff=kb*tk*(1.+alfa*kn)/6./RAP/visc/pi
!      nsch=(visc/diff/rho)**0.3333
      nsch=0.63**0.33333
      IF (Y.LT.MV) THEN
       CORR = 2.0/(4*Y)**0.166667/(PI*RHOI)**0.33333
       DY = 2*(4*Y/PI/RHOI)**0.333333
       HCI = DY
      ELSE
       DY = DY1
       HCI =HCI0
      ENDIF
      RAXI = HCI/DY
      IF (RAXI.GT.0.99) THEN
        CAPAC=0.5
      ELSE
        EPSILON = SQRT(1.- RAXI*RAXI)
        CAPAC = EPSILON/2./ASIN(EPSILON)
      ENDIF
      REYCI = (0.5*DY + HCI)*VY*RHOA/VISC
      KHI = NSCH*SQRT(REYCI)
      IF (KHI.LT.1.0) THEN
        FV = 1.0 + 0.14*KHI*KHI
      ELSE
        FV = 0.86 + 0.28*KHI
      ENDIF
      fbrpi1=4*pi*capac*dy*diff*fv
      return
      end
      function fbrpi2(x,y)
      implicit real*8 (a-h,o-z)
      fbrpi2=y*fbrpi1(x,y)
      return
      end
      function fbrpi3(x,y)
      implicit real*8 (a-h,o-z)
      fbrpi3=x*fbrpi1(x,y)
      return
      end
      function fbrpi4(x,y)
      implicit real*8 (a-h,o-z)
      fbrpi4=x*y*fbrpi1(x,y)
      return
      end
!   This routine calculates the collection kernel of the snow - aerosol particles
!             due to the Brownian collison
!
      subroutine genbrown_snow(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fbrsnw1,fbrsnw2,fbrsnw3,fbrsnw4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do K=2,K4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 br snw',k,i
      call dtwodq(fbrsnw1,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fbrsnw2,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(fbrsnw3,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(fbrsnw4,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 br snw',k,i
      call dtwodq(fbrsnw1,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fbrsnw2,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fbrsnw3,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fbrsnw4,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 br snw',k,i
      call dtwodq(fbrsnw1,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fbrsnw2,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fbrsnw3,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fbrsnw4,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,20
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53 br snw',k,i
      call dtwodq(fbrsnw1,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fbrsnw2,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fbrsnw3,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fbrsnw4,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of  snow - aerosol scavenging by Brownian'
      return
      end
      real*8 function fbrsnw1(x,y)
      IMPLICIT NONE
      real*8 X,Y
      REAL*8 DY,V,RHOSNOW,RAP
      REAL*8 lamb,k1,k0,kn,nsch,al,bl,cl,alfa
      REAL*8 PI,RG,RHOAER,KB,RA,TK,PRES,VISC,RHOA,DIFF
      parameter(RG=287.0,KB=1.38e-23,RA=1.87e-10)
      parameter(al=1.257,bl=0.4,cl=1.1)  !parameters for alfa  in calc. of  diff
      REAL*8 axrmin,axrmax,msnow100,msnow500,axr,capac,omega,epsilon,   &
     &       khi,pker,reycir,fv
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      CALL TVSNOW (y,v,dy,rhosnow)
      axrmin=0.08d0
      axrmax=1.0d0
      msnow100=3.773e-11
      msnow500=2.22529e-8
      RAP=(3.d0*x/(4.d0*3.14159*RHOAER))**(1.d0/3.d0)
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      RHOA=pres/tk/rg
!      lamb=8.89919d-8
      lamb=KB*TK/(4.d0*1.414*pi*ra*ra*pres)            ! free path of air moleculas
      KN=LAMB/RAP
      alfa = al + bl*dexp(-cl/kn)
      diff=kb*tk*(1.+alfa*kn)/6./RAP/visc/pi
!      nsch=visc/diff/roa
      NSCH = (0.632)**0.33333
      if (y.le.msnow100) then
       axr=axrmin
      endif
      if (y.ge.msnow500) then
       axr = axrmax
      endif
      if (y.gt.msnow100.and.y.lt.msnow500) then
       axr=6.d0*y/(dy**3.d0*PI*RHOSNOW)
      endif
      IF (axr.GT.0.999) THEN
       CAPAC=0.5d0
       OMEGA=PI*DY*DY
      ELSE
       EPSILON=DSQRT(1.d0-AXR*AXR)
       CAPAC=EPSILON/2.d0/ASIN(EPSILON)
       OMEGA=0.5*PI*DY*DY*(1.d0+                                        &
     &        AXR*AXR*DLOG((1.d0+EPSILON)/(1.d0-EPSILON))/              &
     &      2.d0/EPSILON)
      ENDIF      
      PKER=DY*PI
      REYCIR = OMEGA/PKER*V*RHOA/VISC
      KHI = NSCH*SQRT(REYCIR)
      IF (KHI.LT.1.d0) THEN
        FV = 1.0 + 0.14*KHI*KHI
      ELSE
        FV = 0.86 + 0.28*KHI
      ENDIF
      fbrsnw1=4*pi*capac*dy*diff*fv
      RETURN
      END
      function fbrsnw2(x,y)
      implicit real*8 (a-h,o-z)
      fbrsnw2=y*fbrsnw1(x,y)
      return
      end
      function fbrsnw3(x,y)
      implicit real*8 (a-h,o-z)
      fbrsnw3=x*fbrsnw1(x,y)
      return
      end
      function fbrsnw4(x,y)
      implicit real*8 (a-h,o-z)
      fbrsnw4=x*y*fbrsnw1(x,y)
      return
      end
!   This routine calculates the collection kernel of the graupel- aerosol particles
!             due to the Brownian collison
!
      subroutine genbrown_gr(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fbrgr1,fbrgr2,fbrgr3,fbrgr4
      external g1,g2,g3
      external h1,h2
      external h11,h12
      pk=2.d0
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do K=2,K4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 br_gr',k,i
      call dtwodq(fbrgr1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fbrgr2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(fbrgr3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(fbrgr4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 br_gr',k,i
      call dtwodq(fbrgr1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fbrgr2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fbrgr3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fbrgr4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 br_gr',k,i
      call dtwodq(fbrgr1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fbrgr2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fbrgr3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fbrgr4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,20
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(fbrgr1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fbrgr2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fbrgr3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fbrgr4,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of  graupel - aerosol scavenging by Brownian'
      return
      end
      real*8 function fbrgr1(x,y)
      IMPLICIT NONE
      real*8 X,Y
      REAL*8 RY,RHOGR,VY,REY,RAP
      REAL*8 lamb,kn,nsch,al,bl,cl,alfa
      REAL*8 PI,TK,PRES,VISC,RHOAER,KB,RA,DIFF,FV
      parameter(KB=1.38e-23,RA=1.87e-10)
      parameter(al=1.257,bl=0.4,cl=1.1)  !parameters for alfa  in calc. of  diff
      REAL*8 M450,M800
      COMMON /RHOAER/ RHOAER
      PI=dasin(1.d0)*2.d0
      m450=4.888790e-9
      m800=4.289321169e-6
      TK = 273.15 + 20.0
      PRES = 101300.00
      VISC=(1.718+0.0049*(tk-273.15))*1e-5
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      CALL TVGR (y,rhogr,vy,ry,rey)
      nsch=0.63d0**0.33333d0
      rap=(3.d0*x/(4.d0*PI*RHOAER))**(1.d0/3.d0)
!      lamb=8.89919d-8
      lamb=KB*TK/(4.d0*1.414*pi*ra*ra*pres)            ! free path of air moleculas
      KN=LAMB/RAP
      alfa = al + bl*dexp(-cl/kn)
      diff=kb*tk*(1.+alfa*kn)/6./RAP/visc/pi 
      fv=0.78d0+0.308d0*NSCH*sqrt(rey)
      fbrgr1=4*pi*ry*diff*fv
      RETURN
      END
      function fbrgr2(x,y)
      implicit real*8 (a-h,o-z)
      fbrgr2=y*fbrgr1(x,y)
      return
      end
      function fbrgr3(x,y)
      implicit real*8 (a-h,o-z)
      fbrgr3=x*fbrgr1(x,y)
      return
      end
      function fbrgr4(x,y)
      implicit real*8 (a-h,o-z)
      fbrgr4=x*y*fbrgr1(x,y)
      return
      end
!   This routine calculates the collection kernel of the water drop - aerosol particles
!             due to the phoretic transport
!      the size intervall of aerosol particles is between 0.062 and 0.246 um
      subroutine genphoresis(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fph1,fph2,fph3,fph4,fphv1,fphv2,fphv3,fphv4
      external g1,g2,g3
      external h1,h2
      external h11,h12
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo

      errabs=0.d0
      errrel=5.d-4
      irule=3
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 ph',k,i
      call dtwodq(fph1,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,1)=result
      call dtwodq(fph2,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,2)=result
      call dtwodq(fph3,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,3)=result
      call dtwodq(fph4,bot,top,g1,h1,errabs,errrel,irule,result,errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 ph',k,i
      call dtwodq(fph1,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,5)=result
      call dtwodq(fph2,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,6)=result
      call dtwodq(fph3,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,7)=result
      call dtwodq(fph4,bot,top,g2,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 ph',k,i
      call dtwodq(fph1,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,9)=result
      call dtwodq(fph2,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,10)=result
      call dtwodq(fph3,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,11)=result
      call dtwodq(fph4,bot,top,h1,h2,errabs,errrel,irule,result,errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(fph1,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,1)=result
      call dtwodq(fph2,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,2)=result
      call dtwodq(fph3,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,3)=result
      call dtwodq(fph4,bot,top,h1,g2,errabs,errrel,irule,result,errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of water  - aerosol scav. by phoretic transport'
      return
      end
      real*8 function fph1(x,y)
      implicit real*8 (a-h,o-z)
      fph1=y**(1.d0/3.d0)
      return
      end
      function fph2(x,y)
      implicit real*8 (a-h,o-z)
      fph2=y*fph1(x,y)
      return
      end
      function fph3(x,y)
      implicit real*8 (a-h,o-z)
      fph3=x*fph1(x,y)
      return
      end
      function fph4(x,y)
      implicit real*8 (a-h,o-z)
      fph4=x*y*fph1(x,y)
      return
      end
      real*8 function fphv1(y,x)
      implicit real*8 (a-h,o-z)
      fphv1=y**(1.d0/3.d0)
      return
      end
      function fphv2(x,y)
      implicit real*8 (a-h,o-z)
      fphv2=y*fphv1(x,y)
      return
      end
      function fphv3(x,y)
      implicit real*8 (a-h,o-z)
      fphv3=x*fphv1(x,y)
      return
      end
      function fphv4(x,y)
      implicit real*8 (a-h,o-z)
      fphv4=x*y*fphv1(x,y)
      return
      end
!   This routine calculates the collection kernel of the pristine ice - aerosol particles
!             due to the phoretic transport
!      the size intervall of aerosol particles is between 0.062 and 0.246 um
      subroutine genphoresis_pice(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fphpi1,fphpi2,fphpi3,fphpi4
      external g1,g2,g3
      external h1,h2
      external h11,h12
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 ph',k,i
      call dtwodq(fphpi1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fphpi2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(fphpi3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(fphpi4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 ph',k,i
      call dtwodq(fphpi1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fphpi2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fphpi3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fphpi4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 ph',k,i
      call dtwodq(fphpi1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fphpi2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fphpi3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fphpi4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(fphpi1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fphpi2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fphpi3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fphpi4,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*)'end of p. ice  - aerosol scav. by phoretic transport'
      return
      end
      real*8 function fphpi1(x,y)
      implicit real*8 (a-h,o-z)
      DY =16.28*DSQRT(Y)
      fphpi1=DY/3.14159
      return
      end
      function fphpi2(x,y)
      implicit real*8 (a-h,o-z)
      fphpi2=y*fphpi1(x,y)
      return
      end
      function fphpi3(x,y)
      implicit real*8 (a-h,o-z)
      fphpi3=x*fphpi1(x,y)
      return
      end
      function fphpi4(x,y)
      implicit real*8 (a-h,o-z)
      fphpi4=x*y*fphpi1(x,y)
      return
      end

!   This routine calculates the collection kernel of the snowflakes - aerosol particles
!             due to the phoretic transport
!      the size intervall of aerosol particles is between 0.062 and 0.246 um
      subroutine genphoresis_snow(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fphsnw1,fphsnw2,fphsnw3,fphsnw4
      external g1,g2,g3
      external h1,h2
      external h11,h12
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 ph_snw',k,i
      call dtwodq(fphsnw1,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fphsnw2,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(fphsnw3,bot,top,g1,h1,errabs,errrel,irule,result,     & 
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(fphsnw4,bot,top,g1,h1,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 ph_snw',k,i
      call dtwodq(fphsnw1,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fphsnw2,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fphsnw3,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fphsnw4,bot,top,g2,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 ph_snw',k,i
      call dtwodq(fphsnw1,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fphsnw2,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fphsnw3,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fphsnw4,bot,top,h1,h2,errabs,errrel,irule,result,     &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53_snw',k,i
      call dtwodq(fphsnw1,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fphsnw2,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fphsnw3,bot,top,h1,g2,errabs,errrel,irule,result,     & 
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fphsnw4,bot,top,h1,g2,errabs,errrel,irule,result,     &
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*) 'end of  snow-aerosol scav. by phoretic transport'
      return
      end
      real*8 function fphsnw1(x,y)
      IMPLICIT NONE
      real*8 X,Y
      REAL*8 VY,DY,RHOSNOW
      REAL*8 PI,axrmin,axrmax, msnow100,msnow500,axr,capac,epsilon
      PI=dasin(1.d0)*2.d0
      CALL TVSNOW (Y,VY,DY,rhosnow)
      axrmin=0.08d0
      axrmax=1.0d0
      msnow100=3.773e-11
      msnow500=2.22529e-8
      if (y.le.msnow100) then
       AXR = AXRMIN
      endif
      if (y.ge.msnow500) then
       AXR = AXRMAX
      endif
      if (y.gt.msnow100.and.y.lt.msnow500) then
       AXR=6.d0*Y/(DY**3.d0*PI*RHOSNOW)
      endif
      IF (axr.GT.0.999) THEN
       CAPAC=0.5d0
      ELSE
       EPSILON=DSQRT(1.d0-AXR*AXR)
       CAPAC=EPSILON/2.d0/ASIN(EPSILON)
      ENDIF
      fphsnw1=capac*dy
      return
      end
      function fphsnw2(x,y)
      implicit real*8 (a-h,o-z)
      fphsnw2=y*fphsnw1(x,y)
      return
      end
      function fphsnw3(x,y)
      implicit real*8 (a-h,o-z)
      fphsnw3=x*fphsnw1(x,y)
      return
      end
      function fphsnw4(x,y)
      implicit real*8 (a-h,o-z)
      fphsnw4=x*y*fphsnw1(x,y)
      return
      end
!   This routine calculates the collection kernel of the graupel - aerosol particles
!             due to the phoretic transport
!      the size intervall of aerosol particles is between 0.062 and 0.246 um
      subroutine genphoresis_gr(k4)
      implicit real*8 (a-h,o-z)
      real*8 ccc(36,20,12),c5(36,20,4)
      real*8 m(40),m1(55)
!      character filen*40
      common/m/ m
      common/m1/ m1
      common/k/ k
      external fphgr1,fphgr2,fphgr3,fphgr4
      external g1,g2,g3
      external h1,h2
      external h11,h12
!
!  cleaning all matrixes
!
      do i=1,12
        do k=1,k4
          do j=1,20
          ccc(k,j,i)=0.d0
          enddo
        enddo
      enddo
      do i=1,4
        do k=1,k4
          do j=1,20
           c5(k,j,i)=0.d0
          enddo
        enddo
      enddo
      errabs=0.d0
      errrel=5.d-4
      irule=3
      do k=2,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'13 ph_gr',k,i
      call dtwodq(fphgr1,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,1)=result
      call dtwodq(fphgr2,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,2)=result
      call dtwodq(fphgr3,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,3)=result
      call dtwodq(fphgr4,bot,top,g1,h1,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,4)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'23 ph',k,i
      call dtwodq(fphgr1,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,5)=result
      call dtwodq(fphgr2,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,6)=result
      call dtwodq(fphgr3,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,7)=result
      call dtwodq(fphgr4,bot,top,g2,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,8)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'33 ph',k,i
      call dtwodq(fphgr1,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,9)=result
      call dtwodq(fphgr2,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,10)=result
      call dtwodq(fphgr3,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,11)=result
      call dtwodq(fphgr4,bot,top,h1,h2,errabs,errrel,irule,result,      &
     &            errest)
      ccc(k,i,12)=result
      enddo
      enddo
      do k=1,k4
      do i=1,19
      bot=m1(i)
      top=m1(i+1)
      write(*,*)'53',k,i
      call dtwodq(fphgr1,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,1)=result
      call dtwodq(fphgr2,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,2)=result
      call dtwodq(fphgr3,bot,top,h1,g2,errabs,errrel,irule,result,      &
     &            errest)
      c5(k,i,3)=result
      call dtwodq(fphgr4,bot,top,h1,g2,errabs,errrel,irule,result,      &                          
     &            errest)
      c5(k,i,4)=result
      enddo
      enddo
      write(3)ccc
      write(3)c5
      close (3)
      write(*,*)'end of graupel  - aerosol scav. by phoretic transport'
      return
      end
      real*8 function fphgr1(x,y)
      IMPLICIT NONE
      REAL*8 X,Y
      REAL*8 SEG,RHOGR,M450,M800,PI
      m450=4.888790e-9
      m800=4.289321169e-6
      if (y.lt.m450) then
       rhogr = 450.0
      endif
      if (y.gt.m800) then
       rhogr = 800.0
      endif
      if (y.gt.m450.and.y.lt.m800) then
       rhogr=(800.0 - 450.0)/(m800 - m450)*(y - m450)+ 450.0
      endif
      PI=3.14159
      SEG=1./(4./3.*pi*RHOGR)**(1./3.)
      fphgr1=y**(1.d0/3.d0)*SEG
      return
      end
      function fphgr2(x,y)
      implicit real*8 (a-h,o-z)
      fphgr2=y*fphgr1(x,y)
      return
      end
      function fphgr3(x,y)
      implicit real*8 (a-h,o-z)
      fphgr3=x*fphgr1(x,y)
      return
      end
      function fphgr4(x,y)
      implicit real*8 (a-h,o-z)
      fphgr4=x*y*fphgr1(x,y)
      return
      end


! 
!    block data for the collision efficency
!
      block data
      implicit none
      real*8 ewd0(0:12,0:20)
      real*8 ehdgrwd0(0:18,0:43)
      real*8 egrwd0(0:18,0:43)
      common /collwd/ ewd0                    ! water - water drop coal. eff
      common /collhdgrwd/ ehdgrwd0            ! high density gra. - water drop coal. eff
      common /collgrwd/ egrwd0                ! gra. - water drop coal. eff
!    The collection efficiencies are calculated on the base of the work of
!    Hall , 1980:A detailed microphysical model within a two-dimensional
!    dynamic framework: model description and preliminary results.
!    J. Atm. Sci. 37, 2486-2507

!----+------------------------------------------------------------------+----
      data ewd0/0,300,200,150,100,70,60,50,40,30,20,10,5.0,0.05,0.97,   &
     & 0.87 ,0.77,0.50,0.20,0.05,0.005,0.001,3*0.0001,0.0,0.1,1.,0.96,  &
     & 0.93,0.79,0.58,0.43,0.40,0.07,0.002,2*0.0001,0.0,0.15,1.0,0.98,  &
     & 0.97,0.91,0.75,0.64,0.60,0.28,0.02,0.005,0.0001,0.0,0.2,2*1.,0.97&
     & ,0.95,0.84,0.77,0.70,0.50,0.04,0.016,0.014,0.0,0.25,3*1.0,0.95,  &
     & 0.88,0.84,0.78,0.62,0.085,0.022,0.017,0.0,0.30,4*1.0,0.90,0.87,  &
     & 0.83,0.68,0.17,0.03,0.019,0.0,0.35,4*1.0,0.92,0.89,0.86,0.74,    &
     & 0.27,0.043,0.022,0.0,0.40,4*1.0,0.94,0.90,0.88,0.78,0.40,0.052,  &
     & 0.027,0.0,0.45,4*1.0,0.95,0.91,0.90,0.80,0.50,0.064,0.03,0.0,    &
     & 0.5,4*1.0,0.95,0.91,0.90,0.80,0.55,0.072,0.033,0.0,0.55,4*1.0,   &
     & 0.95,0.91,0.90,0.80,0.58,0.079,0.035,0.0,0.60,4*1.0,0.95,0.91,   &
     & 0.90,0.78,0.59,0.082,0.037,0.0,0.65,4*1.0,0.95,0.91,0.89,0.77,   &
     & 0.58,0.08,0.038,0.0,0.7,4*1.0,0.95,0.92,0.88,0.76,0.54,0.076,    &
     & 0.038,0.0,0.75,4*1.0,0.97,0.93,0.88,0.77,0.51,0.067,0.037,0.0,   &
     & 0.8,5*1.0,0.95,0.89,0.77,0.49,0.057,0.036,0.0,0.85,4*1.0,1.02,   &
     & 1.0,0.92,0.78,0.47,0.048,0.035,0.0,0.90,4*1.0,1.04,1.03,1.01,    &
     & 0.79,0.45,0.04,0.032,0.0,0.95,4*1.0,2.3,1.7,1.3,0.95,0.47,       &
     &  0.033,0.029,0.0,1.0,4*1.0,4.0,3.0,2.3,1.4,                      &
     &   0.52,0.027,0.027,0.0/
! collision efficienci  given by Khain et al.
!----+------------------------------------------------------------------+----
      data ehdgrwd0/0,1.,5.,10.,20.,50.,80.,100,120,140,160,180,200,220,&
     &       240,260,280,300,320,                                       &
     & 1.,18*0.0,                                                       &
     & 2.,11*0.0,0.01,0.03,0.06,0.08,0.10,0.13,0.15,                    &
     & 4.,5*0.0,.12,.19,.26,.32,.37,.41,.45,.48,.51,.54,.57,            &
     &    .59,.62,                                                      &
     & 6.,4*.00,.18,.31,.47,.54,.59,.63,.66,.68,.71,.73,.75,.77,        &
     &    .79,.81,                                                      &
     & 8.,4*.00,.36,.45,.64,.69,.73,.76,.78,.80,.82,.83,.84,.86,        &
     &    .87,.89,                                                      &
     & 10.,3*.00,.24,.48,.54,.75,.78,.81,.83,.85,.86,.88,.88,.90,.91,   &
     &     .92,.93,                                                     &
     & 12.,2*.00,.12,.36,.56,.60,.81,.84,.86,.88,.89,.90,.91,.92,.93,   &
     &     .93,.94,.95,                                                 &
     & 14.,2*.00,.24,.45,.62,.66,.85,.88,.90,.91,.92,.93,.93,.94,.94,   &
     &     .95,.96,.96,                                                 &
     & 16.,2*.00,.33,.52,.68,.70,.88,.90,.92,.93,.94,.94,.95,.95,.96,   &
     &     .96,.96,.97,                                                 &
     & 18.,.00,.17,.40,.57,.72,.73,.90,.92,.93,.94,.95,.95,.96,.96,.96, &
     &     .97,.97,.97,                                                 &
     & 20.,.00,.25,.46,.62,.75,.75,.92,.93,.94,.95,.96,.96,.96,.97,.97, &
     &     .97,.98,.98,                                                 &
     & 22.,.00,.31,.51,.66,.77,.77,.93,.94,.95,.96,.96,.97,.97,.97,.97, &
     &     .98,.98,.98,                                                 &
     & 24.,.00,.37,.55,.69,.80,.80,.94,.95,.96,.97,.97,.97,.97,.98,.98, &
     &     .98,.98,.99,                                                 &
     & 26.,.00,.41,.58,.72,.82,.85,.94,.96,.96,.97,.97,.97,.98,.98,.98, &
     &     .98,.99,.99,                                                 &
     & 28.,.00,.46,.62,.74,.83,.86,.95,.96,.97,.97,.97,.98,.98,.98,.98, &
     &     .99,.99,.99,                                                 &
     & 30.,.00,.49,.65,.76,.85,.90,.95,.96,.97,.97,.98,.98,.98,.99,.99, &
     &     .99,.99,.99,                                                 &
     & 32.,.00,.52,.67,.78,.86,.91,.96,.97,.97,.98,.98,.98,.99,.99,.99, &
     &     .99,.99,.99,                                                 &
     & 34.,.00,.55,.69,.80,.87,.92,.96,.97,.98,.98,.98,.99,.99,.99,.99, &
     &     .99,.99,.99,                                                 &
     & 36.,.00,.58,.72,.82,.88,.93,.96,.97,.98,.98,.99,.99,.99,.99,.99, &
     &     .99,.99,.99,                                                 &
     & 38.,.00,.60,.73,.83,.88,.93,.96,.97,.98,.98,.99,.99,.99,.99,.99, &
     &     .99,.99,.99,                                                 &
     & 40.,.00,.62,.75,.84,.89,.94,.96,.97,.98,.99,.99,.99,.99,.99,.99, &
     &     .99,.99,.99,                                                 &
     & 50.,.00,.71,.82,.88,.92,.93,.96,.98,.99,.99,.99,.99,.99,.99,.99, &
     &     .99,.99,1.00,                                                &
     & 60.,.00,.77,.86,.91,.93,.96,.94,.98,.99,.99,.99,.99,.99,.99,.99, &
     &     .99,1.00,1.00,                                               &
     & 70.,.00,.81,.89,.93,.94,.90,.74,.97,.99,.99,.99,.99,.99,.99,     &
     &    1.00,1.00,1.00,1.00,                                          &
     & 80.,.00,.85,.91,.94,.95,.99,1.16,.95,.98,.99,.99,.99,.99,1.00,   &
     &    1.00,1.00,1.00,1.00,                                          &
     & 90.,.00,.87,.92,.95,.95,.98,.99,.83,.98,.99,.99,.99,.99,1.00,    &
     &    1.00,1.00,1.00,1.00,                                          &
     & 100.,.00,.89,.93,.96,.96,.98,.99,1.36,.97,.99,.99,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 110.,.00,.90,.94,.96,.96,.98,.99,1.00,.93,.99,.99,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 120.,.00,.92,.95,.97,.96,.98,.99,0.99,.00,.98,.99,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 130.,.00,.92,.96,.97,.98,.98,.99,.99,1.02,.97,.99,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 140.,.00,.93,.96,.98,.98,.99,.99,.99,0.99,.91,.99,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 150.,.00,.94,.96,.98,.99,.99,.99,.99,.99,1.15,.98,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 160.,.00,.95,.97,.98,.99,.99,.99,.99,.99,1.00,.99,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 170.,.00,.95,.97,.98,.99,.99,.99,.99,.99,.99,3.59,.99,.99,1.00,  &
     &     1.00,1.00,1.00,1.00,                                         &
     & 180.,.00,.96,.97,.98,.99,.99,.99,.99,.99,.99,1.03,1.01,.99,1.00, &
     &     1.00,1.00,1.00,1.00,                                         &
     & 190.,.00,.96,.98,.99,.99,.99,.99,.99,.99,.99,1.00,2.14,1.00,1.00,&
     &     1.00,1.00,1.00,1.00,                                         &
     & 200.,.00,.96,.98,.99,.99,.99,.99,.99,.99,.99,1.00,1.12,1.02,1.00,&
     &     1.00,1.00,1.00,1.00,                                         &
     & 210.,.00,.97,.98,.99,.99,.99,.99,1.00,1.00,0.99,1.00,1.01,1.42,  &
     &     1.01,4*1.00,                                                 &
     & 220.,.00,.97,.98,.99,.99,.99,.99,1.00,1.00,1.00,1.00,1.00,1.34,  &
     &     1.04,1.01,3*1.00,                                            &
     & 230.,.00,.97,.98,.99,.99,.99,6*1.00,1.04,1.52,1.02,1.01,1.00,    &
     &     1.00,                                                        &
     & 240.,.00,.97,.98,.99,.99,.99,6*1.00,1.01,1.49,1.09,1.01,1.01,    &
     &     1.00,                                                        &
     & 250.,.00,.97,.99,.99,.99,.99,7*1.00,1.06,3.75,1.06,1.02,1.00,    &
     & 260.,.00,17*1.00/      
! ! collision effcienci  given by Khain1 et al.  (the density of the graupel changes with its size
! the collison coefficients for rhogr = 400 kg/m3 is used for graupel size less than 260 um, and
! the collison coefficients for rhogr = 800 kg/m3 is used for graupel size larger and equal than 260 um
       data egrwd0/0,1.,5.,10.,20.,50.,80.,100,120,140,160,180,200,220, &
     &     240,260,280,300,320,                                         &
     &     1.,18*0.0,                                                   &
     &     2.,11*0.0,0.00,0.00,0.04,0.08,0.10,0.13,0.15,                &
     &     4.,5*0.0,.00,.02,.07,.12,.16,.20,.24,.35,.51,.54,.57,        &
     &   .59,.62,                                                       &
     &     6.,4*.00,.05,.12,.22,.29,.35,.40,.44,.48,.52,.65,.75,.77,    &
     &   .79,.81,                                                       &
     &     8.,4*.00,.16,.31,.41,.48,.54,.58,.62,.65,.67,.73,.84,.86,    &
     &   .87,.89,                                                       &
     &     10.,3*.00,.01,.31,.43,.55,.61,.66,.70,.73,.75,.77,.85,.90,   &
     &   .91,.92,.93,                                                   &
     &     12.,2*.00,.002,.12,.41,.51,.64,.70,.74,.77,.80,.82,.83,.90,  &
     &   .93,.93,.94,.95,                                               &
     &     14.,2*.00,.002,.24,.49,.59,.71,.76,.79,.82,.84,.86,.87,.90,  &
     &   .94,.95,.96,.96,                                               &
     &     16.,2*.00,.07,.33,.56,.65,.75,.80,.83,.86,.87,.89,.90,.93,   &
     &   .96,.96,.96,.97,                                               &
     &     18.,.00,.00,.17,.40,.61,.70,.78,.83,.86,.88,.90,.91,.91,.93, &
     &   .96,.97,.97,.97,                                               &
     &     20.,.00,.00,.25,.45,.65,.73,.81,.85,.88,.90,.91,.92,.93,.95, &
     &   .97,.97,.98,.98,                                               &
     &     22.,.00,.05,.31,.50,.68,.75,.83,.87,.90,.91,.93,.93,.94,.96, &
     &   .97,.98,.98,.98,                                               &
     &     24.,.00,.13,.37,.55,.71,.78,.84,.88,.91,.93,.94,.94,.95,.96, &
     &   .98,.98,.98,.99,                                               &
     &     26.,.00,.19,.41,.58,.74,.78,.84,.89,.92,.93,.94,.95,.96,.97, &
     &   .98,.98,.99,.99,                                               &
     &     28.,.00,.24,.45,.62,.76,.79,.85,.90,.93,.94,.95,.96,.96,.97, &
     &   .98,.99,.99,.99,                                               &
     &     30.,.00,.29,.49,.64,.78,.80,.85,.90,.93,.94,.95,.96,.96,.97, &
     &   .99,.99,.99,.99,                                               &
     &     32.,.00,.33,.52,.67,.80,.82,.84,.91,.93,.95,.96,.96,.97,.98, &
     &   .99,.99,.99,.99,                                               &
     &     34.,.00,.36,.55,.69,.81,.83,.84,.91,.94,.95,.96,.97,.97,.98, &
     &   .99,.99,.99,.99,                                               &
     &     36.,.00,.40,.58,.71,.83,.83,.82,.91,.94,.96,.96,.97,.97,.98, &
     &   .99,.99,.99,.99,                                               &
     &     38.,.00,.43,.60,.73,.83,.85,.79,.91,.94,.96,.97,.97,.97,.98, &
     &   .99,.99,.99,.99,                                               &
     &     40.,.00,.46,.62,.74,.85,.86,.74,.90,.94,.96,.97,.97,.97,.98, &
     &   .99,.99,.99,.99,                                               &
     &     50.,.00,.57,.71,.82,.89,.80,.00,.81,.93,.96,.97,.98,.98,.98, &
     &   .99,.99,.99,1.00,                                              &
     &     60.,.00,.65,.77,.85,.91,.96,1.03,.00,.87,.95,.97,.98,.98,    &
     &   .99,.99,.99,1.00,1.00,                                         &
     &     70.,.00,.71,.81,.88,.93,.90,.99,1.09,.00,.92,.96,.98,.98,    &
     &   .99,1.00,1.00,1.00,1.00,                                       &
     &     80.,.00,.75,.84,.90,.94,.99,.99,1.0,1.28,.70,.94,.97,.98,    &
     &     1.00,1.00,1.00,1.00,1.00,                                    &
     &     90.,.00,.79,.87,.92,.95,.98,.99,.99,1.02,3.85,.87,.96,.98,   &
     &   .99,1.00,1.00,1.00,1.00,                                       &
     &     100.,.00,.82,.89,.93,.96,.98,.99,.99,1.00,1.07,.00,.93,.97,  &
     &   .98,1.00,1.00,1.00,1.00,                                       &
     &     110.,.00,.83,.90,.94,.96,.98,.99,.99,.99,1.01,1.27,.73,.95,  &
     &   .98,1.00,1.00,1.00,1.00,                                       &
     &     120.,.00,.85,.91,.95,.96,.98,.99,.99,.99,1.00,1.03,5.54,.88, &
     &   .96,1.00,1.00,1.00,1.00,                                       &
     &     130.,.00,.87,.92,.96,.98,.98,.99,.99,.99,1.00,1.01,1.11,.00, &
     &   .93,1.00,1.00,1.00,1.00,                                       &
     &     140.,.00,.89,.93,.96,.97,.99,.99,.99,.99,.99,1.00,1.02,1.04, &
     &   .85,1.00,1.00,1.00,1.00,                                       &
     &     150.,.00,.90,.94,.96,.97,.99,.99,.99,.99,.99,1.00,1.01,1.05, &
     &     3.50,1.00,1.00,1.00,1.00,                                    &
     &     160.,.00,.91,.95,.97,.97,.99,.99,.99,.99,.99,1.00,1.00,1.01, &
     &     1.00,1.00,1.00,1.00,1.00,                                    &
     &     170.,.00,.92,.95,.97,.97,.99,.99,.99,.99,.99,1.00,1.00,1.01, &
     &     1.00,1.00,1.00,1.00,1.00,                                    &
     &     180.,.00,.92,.95,.97,.97,.99,.99,.99,.99,1.00,1.00,1.00,1.00,&
     &     1.00,1.00,1.00,1.00,1.00,                                    &
     &     190.,.00,.93,.96,.98,.98,.99,.99,.99,.99,1.00,1.00,1.00,1.00,&
     &     1.00,1.00,1.00,1.00,1.00,                                    &
     &     200.,.00,.93,.96,.98,.99,.99,.99,.99,1.0,1.0,1.0,1.0,1.0,    &
     &     1.00,1.00,1.00,1.00,1.00,                                    &
     &     210.,.00,.94,.97,.98,.98,.99,.99,.99,1.0,1.0,1.0,1.0,1.0,    &
     &     1.00,4*1.00,                                                 &
     &     220.,.00,.95,.96,.98,.99,.99,.99,.99,1.0,1.0,1.0,1.0,1.0,    &
     &     1.00,1.00,3*1.00,                                            &
     &     230.,.00,.95,.97,.98,.99,.99,.99,.99,6*1.00,1.02,1.01,1.00,  &
     &     1.00,                                                        &
     &     240.,.00,.95,.97,.98,.99,.99,.99,7*1.00,1.09,1.01,1.01,1.00, &                                                      
     &     250.,.00,.95,.97,.98,.99,.99,.99,7*1.00,2.75,1.06,1.02,1.00, &
     &     260.,.00,17*1.00/


      end      
!   
!   terminal velocities for different species
!
      subroutine tvwd(x,vx,rx)
      implicit none
      real*8 x,vx,rx
      real*8 lamb,nre,nbo,nbo1,np,nst,nbe,snre
      real*8 rhoa,rhow,g,rg,tk,pres,visc,sigma,pi,seg1,xv,x1,y1,rx1
      rhow=1000.0
      g=9.81
      rg=287.0
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      sigma=(75.93-0.115*(tk-273.15))*1e-3
      pi=asin(1.d0)*2.d0
      seg1=1./(4./3.*pi*1000)**(1./3.)
      rx=x**0.33333d0*seg1      
      rhoa=pres/tk/rg
      lamb=8.89919d-8
      nbe=32.d0*(rhow-rhoa)*rhoa*g/3d0/visc/visc        ! 3.163114206d14  Best num./rd**3 d
      nbo=g*(rhow-rhoa)/sigma                          ! Bond number/rd**2 129082.5717d0
      np=(sigma**3.d0*rhoa*rhoa/(visc**4*g*(rhow-rhoa)))**(1.d0/6.d0)   ! 86.15124908d0
      snre=visc/2./rhoa               !9.619260920d-6           !
      nst=2./9./visc*(rhow-rhoa)*g    ! 1.267784202d8             ! 2./9./visc*(row-roa)*g
      if (rx.le.1d-5) then
       vx=nst*rx*rx*(1.d0+1.26*lamb/rx)
       goto 20
      endif
      if (rx.gt.1.d-5.and.rx.lt.535.d-6) then
       xv=rx**3*nbe
       x1=dlog(xv)
      y1=-3.18657d0+0.992696d0*x1-0.153193d-2*x1*x1-0.987059d-3*x1*x1*x1&
     &   -0.578878d-3*x1**4.d0+0.855176d-4*x1**5.d0-0.327815d-5*x1**6.0
       nre=dexp(y1)
       vx=nre*snre/rx
       goto 20
      endif
      if (rx.gt.535.d-6.and.rx.lt.2.5d-3) then
       nbo1=nbo*rx*rx
       x1=dlog(16.d0/3.d0*nbo1*np)
       y1=-5.00015d0+5.23778d0*x1-2.04914d0*x1*x1+0.475294d0*x1**3-     &
     &    0.0542819d0*x1**4.d0+0.00238449d0*x1**5.d0
       nre=np*dexp(y1)
       vx=nre*snre/rx
       goto 20
      endif
      if (rx.gt.2.5d-3) then
       rx1=2.5d-3
       nbo1=nbo*rx1*rx1
       x1=dlog(16.d0/3.d0*nbo1*np)
       y1=-5.00015d0+5.23778d0*x1-2.04914d0*x1*x1+0.475294d0*x1**3-     &
     &    0.0542819d0*x1**4.d0+0.00238449d0*x1**5.d0
       nre=np*dexp(y1)
       vx=nre*snre/rx1
      endif
 20   continue
      return
      end
      subroutine tvgr (y,rhogr,vy,ry,rey)
      implicit none
      real*8 y,rhogr,vy,ry,rey
      real*8 tk,pres,visc,sigma,rhoa,lamb,rg,g
      real*8 nbe1,nbe2,snre,nst1
      real*8 pi,seg2,xv,x1,y1,bnx,w,lgre
      rg=287.0
      pi=asin(1.d0)*2.d0
      lamb=8.89919d-8
      g=9.81
      tk=273.15+20.
      pres=101300.0
      visc=(1.718+0.0049*(tk-273.15))*1e-5
      sigma=(75.93-0.115*(tk-273.15))*1e-3
      seg2=1./(4./3.*pi*rhogr)**(1./3.)
      ry=y**0.33333d0*seg2
      rhoa=pres/tk/rg
      nbe1=32.d0*(rhogr-rhoa)*rhoa*g/3d0/visc/visc       ! Best num. /rd**3
      nbe2=8*g*rhoa/pi/visc/visc                      ! Best num/m        7.558128329d10
      snre=visc/2./rhoa               !9.619260920d-6           !
      nst1=2./9./visc*(rhogr-rhoa)*g    !1.140892468d8
      if(ry.le.1.d-5) then
       vy=nst1*ry*ry*(1.d0+1.26*lamb/ry)
       rey=2.*ry*vy*rhoa/visc
       goto 35
      endif
      if(ry.le.63.5d-6) then
       xv=ry**3*nbe1
       x1=dlog(xv)
      y1=-3.18657d0+0.992696d0*x1-0.153193d-2*x1*x1-0.987059d-3*x1*x1*x1&
     &   -0.578878d-3*x1**4.d0+0.855176d-4*x1**5.d0-0.327815d-5*x1**6.0
       rey=dexp(y1)
       vy=rey*snre/ry
       goto 35
      endif
      bnx=nbe2*y
      if (bnx.le.562.d0) then
       w=dlog10(bnx)
       lgre=-1.7095d0+1.33438d0*w-0.11591d0*w*w
       rey=10.d0**lgre                                                ! Reynolds number
       goto 30
      endif
      if (bnx.le.1830.d0) then
       w=dlog10(bnx)
       lgre=-1.81391d0+1.34671d0*w-0.12427d0*w*w+.0063d0*w*w*w
       rey=10.d0**lgre
       goto 30
      endif
      if (bnx.le.3.46d8) then
       rey=0.4487d0*bnx**0.5536d0
       goto 30
      endif
      rey=sqrt(bnx/.6d0)
30    vy=visc*rey/(2.d0*ry*rhoa)
35    continue
      return
      end
!
!     this subroutine calculates the terminal velocity of the pristine ice
!
      subroutine tvpice(x,vx,dx)
      implicit none
      real*8 x,vx,dx
      dx=16.28d0*dsqrt(x)
      vx =304.d0*dx
      if (vx.gt.1.d0) then
       vx=1.d0
      endif 
      return
      end
!
!   terminal velocity and mass - size relation for snowflakes
!
      subroutine tvsnow(x,vx,dx,rhosnow)
      implicit none
      real*8 x,vx,dx
      real*8 msnow100,msnow500,msnowtv,axrmin,axrmax,rhomax,rhomin
      real*8 axr,rhosnow,a1,a2,b1,b2,pi
      pi=asin(1.d0)*2.d0
      axrmin=0.08d0
      axrmax=1.0d0
      rhomin=900.d0
      rhomax=340.d0
      msnow100=3.773e-11
      msnow500=2.22529e-8
      msnowtv=3.2906e-8
      if (x.le.msnow100) then
       dx=16.28*dsqrt(x)
       rhosnow=rhomin
      endif
      if (x.ge.msnow500) then
       dx=dsqrt(6.0*x/pi/0.17)
       RHOSNOW=0.17/DX
       IF (RHOSNOW.LT.20) RHOSNOW=20.0
       if ( (0.17/dx) .lt. 20.0) then         
        dx = (6.0*x/3.14159/1.0/20.0)**0.33333333
       endif
      endif
      if (x.gt.msnow100.and.x.lt.msnow500) then
       a1=(axrmax-axrmin)/dlog10(msnow500/msnow100)
       b1=axrmax-a1*dlog10(msnow500)
       a2=(rhomax-rhomin)/dlog10(msnow500/msnow100)
       b2=rhomax-a2*dlog10(msnow500)
       axr=a1*dlog10(x)+b1
       rhosnow=a2*dlog10(x)+b2
       dx=(6.*x/pi/axr/rhosnow)**0.333333
      endif
      if (dx .lt. 432.85E-6) then
       vx=1250.d0*dx
      elseif (dx .lt. 4.347E-3) then
       vx=40.0*dx**0.55 * dexp(-100.0*dx)
      else
       vx=1.301
      endif
      return
      end
