!##############################################################################
! module vacphys - mhd

!##############################################################################
! module vacphys.mhd0 - common subroutines for mhd and mhdiso





!=============================================================================
subroutine physini

! Tell VAC which variables are vectors, set default entropy coefficients

include 'vacdef.f'
integer:: il
!-----------------------------------------------------------------------------

iw_vector(1)=m0_; iw_vector(2)=b0_

! The values of the constants are taken from Ryu & Jones ApJ 442, 228
do il=1,nw
   select case(il)
   case(fastRW_,fastLW_,slowRW_,slowLW_)
      entropycoef(il)=0.2
   case(alfvRW_,alfvLW_)
      entropycoef(il)=0.4
   case default
      entropycoef(il)= -one
   end select
end do

return
end

!=============================================================================
subroutine process(count,idimmin,idimmax,w)

! Process w before it is advected in directions idim^LIM, or before save
! count=1 and 2 for first and second (half step) processing during advection
! count=ifile+2 for saving results into the file indexed by ifile

include 'vacdef.f'

integer:: count,idimmin,idimmax
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

logical:: oktime
double precision:: cputime,time1,timeproc
data timeproc /0.D0/

! The processing should eliminate divergence of B.
!-----------------------------------------------------------------------------

oktest=index(teststr,'process')>=1
oktime=index(teststr,'timeproc')>=1

if(oktest)write(*,*)'Process it,idim^LIM,count',it,idimmin,idimmax,count

if(oktime)time1=cputime()


if(count==0)then
   if(divbconstrain)then
      call die('CT module is OFF: setvac -on=ct; make vac')
   endif
else
   ! Use the projection scheme
   call die('Poisson module is OFF: setvac -on=poisson;make vac')
endif


if(oktime)then
   time1=cputime()-time1
   timeproc=timeproc+time1
   write(*,*)'Time.Proc:',time1,timeproc
endif

return
end

!=============================================================================
subroutine getdt(w,ixmin1,ixmin2,ixmax1,ixmax2)

! If resistivity is  not zero, check diffusion time limit for dt

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1
if(oktest)write(*,*)'GetDt'

if(eqpar(eta_)==zero)return


       write(*,*)'Error: Resistive MHD module is OFF'
call die('Recompile with setvac -on=resist or set eqpar(eta_)=0')

return
end

!=============================================================================
subroutine getdivb(w,ixOmin1,ixOmin2,ixOmax1,ixOmax2,divb)

! Calculate div B within ixO

include 'vacdef.f'

integer:: ix1,ix2
integer::          ixOmin1,ixOmin2,ixOmax1,ixOmax2,ixmin1,ixmin2,ixmax1,&
   ixmax2,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),divb(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

oktest=index(teststr,'getdivb')>=1
if(oktest)write(*,*)'getdivb ixO=',ixOmin1,ixOmin2,ixOmax1,ixOmax2

if(fourthorder)then
   ixmin1=ixOmin1-2;ixmin2=ixOmin2-2;ixmax1=ixOmax1+2;ixmax2=ixOmax2+2;
else
   ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;
endif
!divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=zero
!$OMP DO
      do ix1=ixOmin1,ixOmax1
        do ix2=ixOmin2,ixOmax2
           divb(ix1,ix2)=zero
       enddo
      enddo
!$OMP ENDDO





do idim=1,ndim
!   tmp(ixmin1:ixmax1,ixmin2:ixmax2)=w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)&
!      +w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)
!$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
            tmp(ix1,ix2)=w(ix1,ix2,b0_+idim)&
                +w(ix1,ix2,bg0_+idim)
       enddo
      enddo
!$OMP ENDDO
      call gradient4(.false.,tmp,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,tmp2)

!   divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=divb(ixOmin1:ixOmax1,&
!      ixOmin2:ixOmax2)+tmp2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)

!$OMP DO
      do ix1=ixOmin1,ixOmax1
        do ix2=ixOmin2,ixOmax2
             divb(ix1,ix2)=divb(ix1,&
      ix2)+tmp2(ix1,ix2)

       enddo
      enddo
!$OMP ENDDO


enddo

if(oktest)then
   write(*,*)'divb:',divb(ixtest1,ixtest2)
!   write(*,*)'bx=',w(ixtest1-1:ixtest1+1,ixtest2,b1_)
!   write(*,*)'by=',w(ixtest1,ixtest2-1:ixtest2+1,b2_)
!   write(*,*)'x =',x(ixtest1-1:ixtest1+1,ixtest2,1)
!   write(*,*)'y =',x(ixtest1,ixtest2-1:ixtest2+1,2)
!   write(*,*)'dx=',dx(ixtest1,ixtest2,1)
!   write(*,*)'dy=',dx(ixtest1,ixtest2,2)
endif

return
end

!=============================================================================
subroutine getflux(w,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,f,transport)

! Calculate non-transport flux f_idim[iw] within ix^L.
! Set transport=.true. if a transport flux should be added

include 'vacdef.f'

integer:: ix1, ix2, ix3
integer::          ixmin1,ixmin2,ixmax1,ixmax2,iw,idim
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),f(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2), fb(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
logical::          transport
!-----------------------------------------------------------------------------

oktest= index(teststr,'getflux')>=1
if(oktest.and.iw==iwtest)write(*,*)'Getflux idim,w:',idim,w(ixtest1,ixtest2,&
   iwtest)





transport=.true.


select case(iw)
   case(rho_)
!      f(ixmin1:ixmax1,ixmin2:ixmax2)=w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)&
!         *w(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/(w(ixmin1:ixmax1,&
!         ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))
      !$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
            f(ix1,ix2)=w(ix1,ix2,rhob_)&
                *w(ix1,ix2,m0_+idim)/(w(ix1,&
                ix2,rho_)+w(ix1,ix2,rhob_))
       enddo
      enddo
      !$OMP ENDDO
   case(m1_)
      if(idim==1)then

      call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,f)
	  call getptotal_bg(w,ixmin1,ixmin2,ixmax1,ixmax2,fb)
!          fb(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0
    !$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
            fb(ix1,ix2)=0.d0
       enddo
      enddo
      !$OMP ENDDO



!          f(ixmin1:ixmax1,ixmin2:ixmax2)=f(ixmin1:ixmax1,ixmin2:ixmax2)&
!             +fb(ixmin1:ixmax1,ixmin2:ixmax2)-(w(ixmin1:ixmax1,ixmin2:ixmax2,&
!             b1_)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)+w(ixmin1:ixmax1,&
!             ixmin2:ixmax2,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))-&
!                           w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)&
!                              *w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim) !-&
			 !w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_)    !remove for perturbed
		!$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
            f(ix1,ix2)=f(ix1,ix2)&
             +fb(ix1,ix2)-(w(ix1,ix2,&
             b1_)*w(ix1,ix2,bg0_+idim)+w(ix1,&
             ix2,b0_+idim)*w(ix1,ix2,bg1_))-&
                           w(ix1,ix2,b1_)&
                              *w(ix1,ix2,b0_+idim) !-&
			 !w(ix1,ix2,bg0_+idim)*w(ix1,ix2,bg1_)    !remove for perturbed

       enddo
      enddo
      !$OMP ENDDO


      else
!          f(ixmin1:ixmax1,ixmin2:ixmax2)=-(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)&
!             *w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)+w(ixmin1:ixmax1,&
!             ixmin2:ixmax2,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))-&
!                    w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)*w(ixmin1:ixmax1,&
!                       ixmin2:ixmax2,b0_+idim) !-&
		 !-w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_)  !remove for perturbed
!$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
          f(ix1,ix2)=-(w(ix1,ix2,b1_)&
             *w(ix1,ix2,bg0_+idim)+w(ix1,&
             ix2,b0_+idim)*w(ix1,ix2,bg1_))-&
                    w(ix1,ix2,b1_)*w(ix1,&
                       ix2,b0_+idim) !-&
		 !-w(ix1,ix2,bg0_+idim)*w(ix1,ix2,bg1_)  !remove for perturbed

       enddo
      enddo
      !$OMP ENDDO
      endif
   case(m2_)
      if(idim==2)then

          call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,f)
	  call getptotal_bg(w,ixmin1,ixmin2,ixmax1,ixmax2,fb)
!          fb(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0

    !$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
            fb(ix1,ix2)=0.d0
       enddo
      enddo
      !$OMP ENDDO

!          f(ixmin1:ixmax1,ixmin2:ixmax2)=f(ixmin1:ixmax1,ixmin2:ixmax2)&
!             +fb(ixmin1:ixmax1,ixmin2:ixmax2)-(w(ixmin1:ixmax1,ixmin2:ixmax2,&
!             b2_)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)+w(ixmin1:ixmax1,&
!             ixmin2:ixmax2,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_))-&
!                           w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)&
!                              *w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim) !-&
			 !w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)    !remove for perturbed

      !$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
          f(ix1,ix2)=f(ix1,ix2)&
             +fb(ix1,ix2)-(w(ix1,ix2,&
             b2_)*w(ix1,ix2,bg0_+idim)+w(ix1,&
             ix2,b0_+idim)*w(ix1,ix2,bg2_))-&
                           w(ix1,ix2,b2_)&
                              *w(ix1,ix2,b0_+idim) !-&
			 !w(ix1,ix2,bg0_+idim)*w(ix1,ix2,bg2_)    !remove for perturbed
       enddo
      enddo
      !$OMP ENDDO

      else
!          f(ixmin1:ixmax1,ixmin2:ixmax2)=-(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)&
!             *w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)+w(ixmin1:ixmax1,&
!             ixmin2:ixmax2,b0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_))-&
!                    w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)*w(ixmin1:ixmax1,&
!                       ixmin2:ixmax2,b0_+idim) !-&
		 !-w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)  !remove for perturbed

		 !$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
          f(ix1,ix2)=-(w(ix1,ix2,b2_)&
             *w(ix1,ix2,bg0_+idim)+w(ix1,&
             ix2,b0_+idim)*w(ix1,ix2,bg2_))-&
                    w(ix1,ix2,b2_)*w(ix1,&
                       ix2,b0_+idim) !-&
		 !-w(ix1,ix2,bg0_+idim)*w(ix1,ix2,bg2_)  !remove for perturbed

       enddo
      enddo
      !$OMP ENDDO

      endif

   case(e_)

      call getptotal(w,ixmin1,ixmin2,ixmax1,ixmax2,f)
      call getptotal_bg(w,ixmin1,ixmin2,ixmax1,ixmax2,fb)

     !$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
            fb(ix1,ix2)=0.d0
       enddo
      enddo
      !$OMP ENDDO

      !fb(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0

!      f(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
!         +idim)*(f(ixmin1:ixmax1,ixmin2:ixmax2)+fb(ixmin1:ixmax1,&
!         ixmin2:ixmax2))-w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)&
!         *( (w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_))*w(ixmin1:ixmax1,&
!         ixmin2:ixmax2,m1_)+(w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_))&
!         *w(ixmin1:ixmax1,ixmin2:ixmax2,m2_) )-w(ixmin1:ixmax1,ixmin2:ixmax2,&
!         bg0_+idim)*( (w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))*w(ixmin1:ixmax1,&
!         ixmin2:ixmax2,m1_)+(w(ixmin1:ixmax1,ixmin2:ixmax2,b2_))&
!         *w(ixmin1:ixmax1,ixmin2:ixmax2,m2_) ))/(w(ixmin1:ixmax1,&
!         ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))+&
!              !  -w(ix^S,bg0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_))  ! remove for perturbed
!               w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)*w(ixmin1:ixmax1,&
!                  ixmin2:ixmax2,m0_+idim)/(w(ixmin1:ixmax1,ixmin2:ixmax2,&
!                  rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))&
!                  -              w(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
!                  +idim)*( (w(ixmin1:ixmax1,ixmin2:ixmax2,b1_))&
!                  *w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)+(w(ixmin1:ixmax1,&
!                  ixmin2:ixmax2,b2_))*w(ixmin1:ixmax1,ixmin2:ixmax2,m2_) )&
!                  /(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,&
!                  ixmin2:ixmax2,rhob_)) !

!$OMP DO
      do ix1=ixmin1,ixmax1
        do ix2=ixmin2,ixmax2
      f(ix1,ix2)=(w(ix1,ix2,m0_&
         +idim)*(f(ix1,ix2)+fb(ix1,&
         ix2))-w(ix1,ix2,b0_+idim)&
         *( (w(ix1,ix2,bg1_))*w(ix1,&
         ix2,m1_)+(w(ix1,ix2,bg2_))&
         *w(ix1,ix2,m2_) )-w(ix1,ix2,&
         bg0_+idim)*( (w(ix1,ix2,b1_))*w(ix1,&
         ix2,m1_)+(w(ix1,ix2,b2_))&
         *w(ix1,ix2,m2_) ))/(w(ix1,&
         ix2,rho_)+w(ix1,ix2,rhob_))+&
              !  -w(ix^S,bg0_+idim)*( ^C&(w(ix^S,bg^C_))*w(ix^S,m^C_)+ )/(w(ix^S,rho_)+w(ix^S,rhob_))  ! remove for perturbed
               w(ix1,ix2,eb_)*w(ix1,&
                  ix2,m0_+idim)/(w(ix1,ix2,&
                  rho_)+w(ix1,ix2,rhob_))&
                  -              w(ix1,ix2,b0_&
                  +idim)*( (w(ix1,ix2,b1_))&
                  *w(ix1,ix2,m1_)+(w(ix1,&
                  ix2,b2_))*w(ix1,ix2,m2_) )&
                  /(w(ix1,ix2,rho_)+w(ix1,&
                  ix2,rhob_)) !
       enddo
      enddo
      !$OMP ENDDO

   case(b1_)
      if(idim==1) then
!         f(ixmin1:ixmax1,ixmin2:ixmax2)= zero
         transport=.false.
        !$OMP DO
        do ix1=ixmin1,ixmax1
            do ix2=ixmin2,ixmax2
                f(ix1,ix2)=zero
        enddo
        enddo
        !$OMP ENDDO
      else

!         f(ixmin1:ixmax1,ixmin2:ixmax2)= -w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)&
!            /(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,&
!            ixmin2:ixmax2,rhob_))*(w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)&
!            +w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim))+ &
!                  w(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/(w(ixmin1:ixmax1,&
!                     ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
!                     rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_)
          !$OMP DO
          do ix1=ixmin1,ixmax1
            do ix2=ixmin2,ixmax2
               f(ix1,ix2)= -w(ix1,ix2,m1_)&
                        /(w(ix1,ix2,rho_)+w(ix1,&
                        ix2,rhob_))*(w(ix1,ix2,b0_+idim)&
                        +w(ix1,ix2,bg0_+idim))+ &
                              w(ix1,ix2,m0_+idim)/(w(ix1,&
                                 ix2,rho_)+w(ix1,ix2,&
                                 rhob_))*w(ix1,ix2,bg1_)
           enddo
          enddo
          !$OMP ENDDO




      endif
   case(b2_)
      if(idim==2) then
!         f(ixmin1:ixmax1,ixmin2:ixmax2)= zero
        !$OMP DO
        do ix1=ixmin1,ixmax1
            do ix2=ixmin2,ixmax2
                f(ix1,ix2)=zero
            enddo
        enddo
        !$OMP ENDDO
         transport=.false.
      else

 !        f(ixmin1:ixmax1,ixmin2:ixmax2)= -w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)&
 !           /(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,&
 !           ixmin2:ixmax2,rhob_))*(w(ixmin1:ixmax1,ixmin2:ixmax2,b0_+idim)&
 !           +w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim))+ &
 !                 w(ixmin1:ixmax1,ixmin2:ixmax2,m0_+idim)/(w(ixmin1:ixmax1,&
 !                    ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
 !                    rhob_))*w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)

          !$OMP DO
          do ix1=ixmin1,ixmax1
            do ix2=ixmin2,ixmax2
             f(ix1,ix2)= -w(ix1,ix2,m2_)&
                /(w(ix1,ix2,rho_)+w(ix1,&
                ix2,rhob_))*(w(ix1,ix2,b0_+idim)&
                +w(ix1,ix2,bg0_+idim))+ &
                      w(ix1,ix2,m0_+idim)/(w(ix1,&
                         ix2,rho_)+w(ix1,ix2,&
                         rhob_))*w(ix1,ix2,bg2_)

           enddo
          enddo
          !$OMP ENDDO





      endif

   case default
      call die('Error in getflux: unknown flow variable')
end select

if(oktest.and.iw==iwtest)write(*,*)'transport,f:',transport,f(ixtest1,ixtest2)

return
end

!=============================================================================
subroutine addsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

! Add sources from resistivity and Powell solver

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
!-----------------------------------------------------------------------------

oktest=index(teststr,'addsource')>=1
if(oktest)write(*,*)'Addsource, compactres,divbfix:',compactres,divbfix
if(oktest)write(*,*)'Before adding source:',wnew(ixtest1,ixtest2,iwtest)

! Sources for resistivity in eqs. for e, B1, B2 and B3
if(abs(eqpar(eta_))>smalldouble)then

          write(*,*)'Error: Resistive MHD module is OFF'
   call die('Recompile with setvac -on=resist or set eqpar(eta_)=0')
endif


! Sources related to div B in the Powell solver
 if(divbfix) call addsource_divb(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
    ixOmin2,ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

if(oktest)write(*,*)'After adding source:',wnew(ixtest1,ixtest2,iwtest)

return
end

!=============================================================================
subroutine addsource_divb(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

! Add Powell's divB related sources to wnew within ixO if possible,
! otherwise shrink ixO

include 'vacdef.f'

integer:: ix1,ix2
integer::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,iws(niw_),iiw,iw
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision:: divb(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
!-----------------------------------------------------------------------------

! Calculating div B involves first derivatives
call ensurebound(1,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,qtC,w)

! We calculate now div B
call getdivb(w,ixOmin1,ixOmin2,ixOmax1,ixOmax2,divb)
!divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=qdt*divb(ixOmin1:ixOmax1,&
!   ixOmin2:ixOmax2)

!$OMP DO
do ix1=ixOmin1,ixOmax1
        do ix2=ixOmin2,ixOmax2
            divb(ix1,ix2)=qdt*divb(ix1,ix2)
       enddo
enddo
!$OMP ENDDO


do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
      case(m1_)
!         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,iw)-(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)&
!            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bg1_))*divb(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2)
      !$OMP DO
        do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
         wnew(ix1,ix2,iw)=wnew(ix1,&
            ix2,iw)-(w(ix1,ix2,b1_)&
            +w(ix1,ix2,bg1_))*divb(ix1,&
            ix2)
            enddo
        enddo
      !$OMP ENDDO


      case(m2_)
!         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,iw)-(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)&
!            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bg2_))*divb(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2)
      !$OMP DO
        do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
         wnew(ix1,ix2,iw)=wnew(ix1,&
            ix2,iw)-(w(ix1,ix2,b2_)&
            +w(ix1,ix2,bg2_))*divb(ix1,&
            ix2)
            enddo
        enddo
      !$OMP ENDDO



      case(b1_)
!         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,iw)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
!            /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,rhob_))*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
       !$OMP DO
        do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
                wnew(ix1,ix2,iw)=wnew(ix1,&
                ix2,iw)-w(ix1,ix2,m1_)&
                /(w(ix1,ix2,rho_)+w(ix1,&
                ix2,rhob_))*divb(ix1,ix2)
            enddo
        enddo
      !$OMP ENDDO




      case(b2_)
!         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,iw)-w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
!            /(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,rhob_))*divb(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      !$OMP DO
        do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
                wnew(ix1,ix2,iw)=wnew(ix1,&
                    ix2,iw)-w(ix1,ix2,m2_)&
                    /(w(ix1,ix2,rho_)+w(ix1,&
                    ix2,rhob_))*divb(ix1,ix2)
            enddo
        enddo
      !$OMP ENDDO
      case(e_)
!         wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw)=wnew(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,iw)-(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m1_)&
!            *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b1_)+w(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,bg1_))+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m2_)&
!            *(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,b2_)+w(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2,bg2_)) )/(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
!            +w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rhob_))*divb(ixOmin1:ixOmax1,&
!            ixOmin2:ixOmax2)
      !$OMP DO
        do ix1=ixOmin1,ixOmax1
            do ix2=ixOmin2,ixOmax2
                wnew(ix1,ix2,iw)=wnew(ix1,&
                ix2,iw)-(w(ix1,ix2,m1_)&
                *(w(ix1,ix2,b1_)+w(ix1,&
                ix2,bg1_))+w(ix1,ix2,m2_)&
                *(w(ix1,ix2,b2_)+w(ix1,&
                ix2,bg2_)) )/(w(ix1,ix2,rho_)&
                +w(ix1,ix2,rhob_))*divb(ix1,&
                ix2)
            enddo
        enddo
      !$OMP ENDDO


   end select
end do

return
end


!=============================================================================
! end module vacphys.mhd0
!##############################################################################

!=============================================================================

subroutine keeppositive(ixmin1,ixmin2,ixmax1,ixmax2,w)

! Keep pressure and density positive (following D.Ryu)

include 'vacdef.f'

integer::          ixmin1,ixmin2,ixmax1,ixmax2
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
logical:: toosmallp
!-----------------------------------------------------------------------------
   ! Where rho is small use vacuum state: rho=vacuumrho, v=0, p=smallp, same B
   where((w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,&
      rhob_))<smallrho)
      w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)=zero
      w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)=zero;
!!!      ^C&w(ix^S,m^C_)=w(ix^S,m^C_)/w(ix^S,rho_)*vacuumrho;
      w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=vacuumrho-w(ixmin1:ixmax1,&
         ixmin2:ixmax2,rhob_)
      w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=smallp/(eqpar(gamma_)-one)&
         +half*(w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2+w(ixmin1:ixmax1,&
         ixmin2:ixmax2,b2_)**2)-w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)
   endwhere
! Calculate pressure without clipping toosmall values (.false.)
call getpthermal(w,ixmin1,ixmin2,ixmax1,ixmax2,tmp)
toosmallp=any(tmp(ixmin1:ixmax1,ixmin2:ixmax2)<max(zero,smallp))
if(toosmallp)then
   nerror(toosmallp_)=nerror(toosmallp_)+1
   if(nerror(toosmallp_)==1)then
      write(*,'(a,i2,a,i7)')'Too small pressure (code=',toosmallp_,') at it=',&
         it
      write(*,*)'Value < smallp: ',minval(tmp(ixmin1:ixmax1,ixmin2:ixmax2)),&
         smallp
!     write(*,*)'Location: ',minloc(tmp(ix^S)) !F77_
   endif
   if(smallp>zero)w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=max(tmp(ixmin1:ixmax1,&
      ixmin2:ixmax2),smallp)/(eqpar(gamma_)-1)+half*((w(ixmin1:ixmax1,&
      ixmin2:ixmax2,m1_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)**2)&
      /w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+(w(ixmin1:ixmax1,ixmin2:ixmax2,&
      b1_)**2+w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)**2))-w(ixmin1:ixmax1,&
      ixmin2:ixmax2,eb_)
endif








return
end

!=============================================================================
! end module vacphys - mhd
!##############################################################################
