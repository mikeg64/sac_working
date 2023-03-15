!#############################################################################
! module vaccd
! Centered difference scheme
!=============================================================================

!=============================================================================
subroutine centdiff4(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,idimmin,idimmax,qtC,wCT,qt,w)

! Advance the iws flow variables from t to t+qdt within ixO^L by 
! fourth order centered  differencing in space the dw/dt+dF_i(w)/dx_i=S 
! type equation. 
! wCT contains the time centered variables at time qtC for flux and source.
! w is the old value at qt on input and the new value at qt+qdt on output.

include 'vacdef.f'

double precision:: qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   iws(niw_),idimmin,idimmax
logical :: transport

double precision:: v(ixGlo1:ixGhi1,ixGlo2:ixGhi2),f(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2), fb(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
integer:: iiw,iw,ixmin1,ixmin2,ixmax1,ixmax2,idim,idir
!-----------------------------------------------------------------------------


! Two extra layers are needed in each direction for which fluxes are added.
ixmin1=ixOmin1;ixmin2=ixOmin2;ixmax1=ixOmax1;ixmax2=ixOmax2;
do idim= idimmin,idimmax
   ixmin1=ixmin1-2*kr(idim,1);ixmin2=ixmin2-2*kr(idim,2)
   ixmax1=ixmax1+2*kr(idim,1);ixmax2=ixmax2+2*kr(idim,2);
enddo
if(ixImin1>ixmin1.or.ixImin2>ixmin2.or.ixImax1<ixmax1.or.ixImax2&
   <ixmax2) call die( 'Error in CentDiff4: Non-conforming input limits')

! Add fluxes to w
do idim= idimmin,idimmax
   ixmin1=ixOmin1-2*kr(idim,1);ixmin2=ixOmin2-2*kr(idim,2)
   ixmax1=ixOmax1+2*kr(idim,1);ixmax2=ixOmax2+2*kr(idim,2);

   call getv(wCT,ixmin1,ixmin2,ixmax1,ixmax2,idim,v)

   do iiw=1,iws(niw_); iw=iws(iiw)
!   print*,'iiw', iiw,idim,idir
      ! Get non-transported flux
      call getflux(wCT,ixmin1,ixmin2,ixmax1,ixmax2,iw,idim,f,transport)

      ! Add transport flux
      if(transport)f(ixmin1:ixmax1,ixmin2:ixmax2)=f(ixmin1:ixmax1,&
         ixmin2:ixmax2)+v(ixmin1:ixmax1,ixmin2:ixmax2)*wCT(ixmin1:ixmax1,&
         ixmin2:ixmax2,iw)

      ! Add divergence of flux
      call gradient4(.false.,f,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,tmp)
      w(ixmin1:ixmax1,ixmin2:ixmax2,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,iw)&
         -qdt*tmp(ixmin1:ixmax1,ixmin2:ixmax2)

   select case(iw)

    case(e_)

         call gradient4(.false.,v,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idim,tmp)   
         call getptotal_bg(w,ixmin1,ixmin2,ixmax1,ixmax2,fb)
         
         w(ixmin1:ixmax1,ixmin2:ixmax2,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,iw)&
            -qdt*tmp(ixmin1:ixmax1,ixmin2:ixmax2)*fb(ixmin1:ixmax1,&
            ixmin2:ixmax2)
         
        do idir= idimmin,idimmax 
             call gradient4(.false.,v,ixOmin1,ixOmin2,ixOmax1,ixOmax2,idir,&
                tmp)   
             w(ixmin1:ixmax1,ixmin2:ixmax2,iw)=w(ixmin1:ixmax1,ixmin2:ixmax2,&
                iw)+qdt*w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idir)&
                *w(ixmin1:ixmax1,ixmin2:ixmax2,bg0_+idim)*tmp(ixmin1:ixmax1,&
                ixmin2:ixmax2)
        enddo

   endselect        


   end do    !next iw
end do       !next idim


if(sourceunsplit) call addsource2(qdt*(idimmax-idimmin+one)&
   /ndim, ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,iws,&
   qtC,wCT,qt,w)

return
end

!=============================================================================
! end module vaccd
!#############################################################################
