!##############################################################################
! module vacusr - sim1 ! setvac -d=22 -g=204,204 -p=hdadiab -u=sim1


!==============================================================================
!
!    THE FOLLOWING SUBROUTINES ADD GRAVITATIONAL SOURCE TERMS, SET GRAVITY
!
!------------------------------------------------------------------------------
!    See vacusr.t.gravity and vacusrpar.t.gravity for an example of usage
!
!    Gravitational force is added to the momentum equation:
!
!    d m_i/dt += rho*eqpar(grav0_+i)
!
!    Gravitational work is added to the energy equation (if present):
!
!    de/dt += Sum_i m_i*eqpar(grav0_+i)
!
!    The eqpar(grav1_),eqpar(grav2_),... coefficients are the components of
!    the gravitational acceleration in each dimension. Set them to 0 for no
!    gravity in that direction.
!    The !!! comments show how a grav array could be used for a spatially
!    (and maybe temporally) varying gravitational field.
!    The setgrav subroutine has to be completed then.
!
!============================================================================
subroutine addsource_grav(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

! Add gravity source calculated from w to wnew within ixO for all variables
! in iws. w is at time qtC, wnew is advanced from qt to qt+qdt.

include 'vacdef.f'

integer::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
   ixOmax2,iws(niw_)
double precision:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: iiw,iw,idim
integer:: ix1, ix2
integer:: ixmin1, ixmin2, ixmax1, ixmax2
!!! ! For a spatially varying gravity define the common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common /gravity/ grav

!-----------------------------------------------------------------------------

!!! ! If grav needs to be calculated only once do it for the whole grid
!!! if(it==itmin)call setgrav(w,ixG^L,ixG^L,grav)
!!! ! Otherwise call setgrav in every time step
!!! call setgrav(w,ixI^L,ixO^L,grav)

!!$OMP DO
!      do ix1=ixmin1,ixmax1
!        do ix2=ixmin2,ixmax2!
!
!       enddo
!      enddo
!!$OMP ENDDO


ixmin1=ixOmin1
ixmax1=ixOmax1
ixmin2=ixOmin2
ixmax2=ixOmax2


! add sources from gravity
do iiw=1,iws(niw_); iw=iws(iiw)
   select case(iw)
   case(m1_,m2_)
      ! dm_i/dt= +rho*g_i
      idim=iw-m0_
      !if(abs(eqpar(grav0_+idim))>smalldouble) wnew(ixOmin1:ixOmax1,&
      !   ixOmin2:ixOmax2,m0_+idim)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_&
      !   +idim)+ qdt*eqpar(grav0_+idim)*(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      !   rho_))

        if(abs(eqpar(grav0_+idim))>smalldouble)
            !$OMP DO
                  do ix1=ixmin1,ixmax1
                    do ix2=ixmin2,ixmax2!
                        wnew(ix1,ix2,m0_+idim)=wnew(ix1,ix2,m0_&
                            +idim)+ qdt*eqpar(grav0_+idim)*(w(ix1,ix2,rho_))
                   enddo
                  enddo
            !$OMP ENDDO
        endif

!          wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
!              qdt*eqpar(grav0_+idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

      !!! ! For a spatially varying gravity use instead of the above lines
      !!! wnew(ixO^S,m0_+idim)=wnew(ixO^S,m0_+idim)+ &
      !!!    qdt*grav(ixO^S,idim)*(w(ixO^S,rho_)+w(ixO^S,rhob_))

   case(e_)
      ! de/dt= +g_i*m_i
      !do idim=1,ndim
      !   if(abs(eqpar(grav0_+idim))>smalldouble) wnew(ixOmin1:ixOmax1,&
      !      ixOmin2:ixOmax2,ee_)=wnew(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ee_)&
      !      + qdt*eqpar(grav0_+idim)*w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)&
      !      *w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,m0_+idim)/(w(ixOmin1:ixOmax1,&
      !      ixOmin2:ixOmax2,rho_)+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rhob_

        do idim=1,ndim
            if(abs(eqpar(grav0_+idim))>smalldouble)
                !$OMP DO
                      do ix1=ixmin1,ixmax1
                        do ix2=ixmin2,ixmax2!
                            wnew(ix1,ix2,ee_)=wnew(ix1,ix2,ee_)&
                                + qdt*eqpar(grav0_+idim)*w(ix1,ix2,rho_)&
                                *w(ix1,ix2,m0_+idim)/(w(ix1,&
                                ix2,rho_)+w(ix1,ix2,rhob_

                       enddo
                      enddo
                !$OMP ENDDO
            endif
        enddo



!            wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
!               qdt*eqpar(grav0_+idim)*w(ixO^S,m0_+idim)

         !!! ! For a spatially varying gravity use instead of the above lines
         !!! wnew(ixO^S,ee_)=wnew(ixO^S,ee_)+ &
         !!!    qdt*grav(ixO^S,idim)*w(ixO^S,m0_+idim)

      end do
   end select ! iw
end do        ! iiw

return
end
!=============================================================================
!!! subroutine setgrav(w,ixI^L,ixO^L,grav)

! Set the gravitational acceleration within ixO based on x(ixI,ndim)
! and/or w(ixI,nw)

!!! include 'vacdef.f'

!!! double precision:: w(ixG^T,nw),grav(ixG^T,ndim)
!!! integer:: ixI^L,ixO^L
!----------------------------------------------------------------------------
!!! return
!!! end
!=============================================================================

subroutine getdt_grav(w,ixmin1,ixmin2,ixmax1,ixmax2)

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idim
double precision:: dtgrav
save dtgrav

!!! ! For spatially varying gravity you need a common grav array
!!! double precision:: grav(ixG^T,ndim)
!!! common/gravity/grav

!----------------------------------------------------------------------------

oktest=index(teststr,'getdt')>=1

!global mixima here!
if(it==itmin)then
   ! If gravity is descibed by the equation parameters, use this:
   dtgrav=bigdouble
   do idim=1,ndim
      if(abs(eqpar(grav0_+idim))>zero)dtgrav=min(dtgrav,one&
         /sqrt(maxval(abs(eqpar(grav0_+idim))/dx(ixMmin1:ixMmax1,&
         ixMmin2:ixMmax2,1:ndim))))
   enddo
   !!! ! For spatially varying gravity use this instead of the lines above:
   !!! call setgrav(w,ixG^L,ixM^L,grav)
   !!! ! If gravity does not change with time, calculate dtgrav here:
   !!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))
endif

!!! ! If gravity changes with time, calculate dtgrav here:
!!! dtgrav=one/sqrt(maxval(abs(grav(ixM^S,1:ndim))/dx(ixM^S,1:ndim)))



! limit the time step
dt=min(dt,dtgrav)
if(oktest)write(*,*)'Gravity limit for dt:',dtgrav

return
end

!=============================================================================
!INCLUDE:vacnul.specialini.t
!INCLUDE:vacnul.specialbound.t
! INCLUDE:vacnul.specialsource.t
!INCLUDE:vacnul.specialio.t
!==============================================================================
SUBROUTINE addsource_visc(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,w,qt,wnew)

  ! Add viscosity source to wnew within ixO

  INCLUDE 'vacdef.f'

  INTEGER::          ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,iws(niw_)
  DOUBLE PRECISION:: qdt,qtC,qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
     wnew(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)

  INTEGER:: ix,ixmin1,ixmin2,ixmax1,ixmax2,idim,idir,jdir,iiw,iw

  !already declared in vacusr.f
  !double precision:: tmp2(ixG^T)
  DOUBLE PRECISION:: nushk(ixGlo1:ixGhi1,ixGlo2:ixGhi2,ndim)



  DOUBLE PRECISION:: tmprhoL(ixGlo1:ixGhi1,ixGlo2:ixGhi2), &
     tmprhoR(ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmprhoC(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)
  DOUBLE PRECISION:: tmpVL(ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmpVR(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2), tmpVC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
  DOUBLE PRECISION:: tmpBL(ixGlo1:ixGhi1,ixGlo2:ixGhi2), tmpBR(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2), tmpBC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

  DOUBLE PRECISION:: tmpL(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tmpR(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2), tmpC(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

  DOUBLE PRECISION:: nuL(ixGlo1:ixGhi1,ixGlo2:ixGhi2),nuR(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)

  INTEGER:: jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,hxmin2,hxmax1,hxmax2, hxOmin1,&
     hxOmin2,hxOmax1,hxOmax2

  DOUBLE PRECISION:: c_ene,c_shk

  INTEGER:: i,j,k,l,m,ii0,ii1,t00

  DOUBLE PRECISION:: sB

  !-----------------------------------------------------------------------------

  ! Calculating viscosity sources
  ! involves second derivatives, two extra layers
  CALL ensurebound(2,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,qtC,w)
  ixmin1=ixOmin1-1;ixmin2=ixOmin2-1;ixmax1=ixOmax1+1;ixmax2=ixOmax2+1;

  !sehr wichtig
  CALL setnushk(w,ixmin1,ixmin2,ixmax1,ixmax2,nushk)

  DO idim=1,ndim
     tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
        rho_)
     CALL setnu(w,rho_,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)
     CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
     tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
     CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
     tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
     wnew(ixImin1:ixImax1,ixImin2:ixImax2,rho_)=wnew(ixImin1:ixImax1,&
        ixImin2:ixImax2,rho_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
        -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
        ixImin2:ixImax2,idim)*qdt
  ENDDO


  DO idim=1,ndim
     tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
        e_)-half*((w(ixImin1:ixImax1,ixImin2:ixImax2,b1_)**2&
        +w(ixImin1:ixImax1,ixImin2:ixImax2,b2_)**2)+(w(ixImin1:ixImax1,&
        ixImin2:ixImax2,m1_)**2+w(ixImin1:ixImax1,ixImin2:ixImax2,m2_)**2)&
        /(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)+w(ixImin1:ixImax1,&
        ixImin2:ixImax2,rhob_)))
     CALL setnu(w,173,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)
     CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
     tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
     CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
     tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
     wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
        ixImin2:ixImax2,e_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
        -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
        ixImin2:ixImax2,idim)*qdt
  ENDDO




  tmprhoC(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
     rho_)+w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_)



  DO k=1,ndim
     jxmin1=ixmin1+kr(k,1);jxmin2=ixmin2+kr(k,2);jxmax1=ixmax1+kr(k,1)
     jxmax2=ixmax2+kr(k,2);
     hxmin1=ixmin1-kr(k,1);hxmin2=ixmin2-kr(k,2);hxmax1=ixmax1-kr(k,1)
     hxmax2=ixmax2-kr(k,2);
     tmprhoL(ixmin1:ixmax1,ixmin2:ixmax2)=((w(ixmin1:ixmax1,ixmin2:ixmax2,&
        rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))+(w(hxmin1:hxmax1,&
        hxmin2:hxmax2,rho_)+w(hxmin1:hxmax1,hxmin2:hxmax2,rhob_)))/two
     tmprhoR(ixmin1:ixmax1,ixmin2:ixmax2)=((w(jxmin1:jxmax1,jxmin2:jxmax2,&
        rho_)+w(jxmin1:jxmax1,jxmin2:jxmax2,rhob_))+(w(ixmin1:ixmax1,&
        ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)))/two

     DO l=1,ndim
        CALL setnu(w,l+m0_,k,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)
        tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,m0_+l)/(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)&
           +w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_))


        DO ii1=0,1
           IF (ii1 .EQ. 0) THEN
              i=k
              ii0=l
           ELSE
              i=l
              ii0=k
           ENDIF



           IF (i .EQ. k) THEN
              tmpVL(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,&
                 ixmin2:ixmax2,m0_+ii0)+w(hxmin1:hxmax1,hxmin2:hxmax2,m0_&
                 +ii0))/two
              tmpVR(ixmin1:ixmax1,ixmin2:ixmax2)=(w(jxmin1:jxmax1,&
                 jxmin2:jxmax2,m0_+ii0)+w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
                 +ii0))/two

              CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
              tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
                 ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,k))&
                 *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
              CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
              tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
                 ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,k))&
                 *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

              tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=(tmprhoR(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
                 -tmprhoL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpL&
                 (ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
                 ixImin2:ixImax2,k)/two

              wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 =wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 +tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*qdt

              tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=(tmpVR(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
                 -tmpVL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpL(ixImin1:ixImax1,&
                 ixImin2:ixImax2))/dx(ixImin1:ixImax1,ixImin2:ixImax2,k)/two

              wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
                 ixImin2:ixImax2,e_)+tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*qdt
           ENDIF




           IF (i .NE. k) THEN
              CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
              tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*(nuL(ixImin1:ixImax1,ixImin2:ixImax2)&
                 +nuR(ixImin1:ixImax1,ixImin2:ixImax2)+two*nushk&
                 (ixImin1:ixImax1,ixImin2:ixImax2,k))/two/two

              tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmprhoC(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
              CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,i,tmpC)

              wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 =wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 +tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt

              tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
                 ixImin2:ixImax2,m0_+ii0)*tmp2(ixImin1:ixImax1,&
                 ixImin2:ixImax2)
              CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,i,tmpC)

              wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
                 ixImin2:ixImax2,e_)+tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt
           ENDIF

        ENDDO
     ENDDO
  ENDDO





  DO k=1,ndim
     DO l=1,ndim

        IF (k .NE. l) THEN

           CALL setnu(w,b0_+l,k,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)

           DO ii1=0,1

              IF (ii1 .EQ. 0) THEN
                 ii0=k
                 m=l
                 sB=-1.d0
                 j=k
              ENDIF

              IF (ii1 .EQ. 1) THEN
                 ii0=l    !ii0 is index B
                 m=k      !first derivative
                 sB=1.d0  !sign B
                 j=l      !first B in energy
              ENDIF



              !print*,'k,l,m,j,ii0,ii1=',k,l,m,j,ii0,ii1



              IF (m .EQ. k) THEN

                 jxmin1=ixmin1+kr(m,1);jxmin2=ixmin2+kr(m,2)
                 jxmax1=ixmax1+kr(m,1);jxmax2=ixmax2+kr(m,2);
                 hxmin1=ixmin1-kr(m,1);hxmin2=ixmin2-kr(m,2)
                 hxmax1=ixmax1-kr(m,1);hxmax2=ixmax2-kr(m,2);
                 tmpBL(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,&
                    ixmin2:ixmax2,b0_+j)+w(hxmin1:hxmax1,hxmin2:hxmax2,b0_&
                    +j))/two
                 tmpBR(ixmin1:ixmax1,ixmin2:ixmax2)=(w(jxmin1:jxmax1,&
                    jxmin2:jxmax2,b0_+j)+w(ixmin1:ixmax1,ixmin2:ixmax2,b0_&
                    +j))/two

                 tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
                    ixImin2:ixImax2,b0_+l)

                 CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
                 tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
                    ixImin2:ixImax2))*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
                 CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)
                 tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
                    ixImin2:ixImax2))*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)

                 wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)&
                    =wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)&
                    +sB*(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
                    -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))&
                    /dx(ixImin1:ixImax1,ixImin2:ixImax2,k)*qdt

                 wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew&
                    (ixImin1:ixImax1,ixImin2:ixImax2,e_)+sB&
                    *(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
                    *tmpBR(ixImin1:ixImax1,ixImin2:ixImax2)&
                    -tmpL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpBL&
                    (ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
                    ixImin2:ixImax2,k)*qdt


              ENDIF



              IF (m .NE. k) THEN

                 tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
                    ixImin2:ixImax2,b0_+l)

                 CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)

                 tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
                    ixImin2:ixImax2)*(nuL(ixImin1:ixImax1,ixImin2:ixImax2)&
                    +nuR(ixImin1:ixImax1,ixImin2:ixImax2))/two

                 CALL gradient1(tmp2,ixmin1,ixmin2,ixmax1,ixmax2,m,tmpC)

                 wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)&
                    =wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)&
                    +sB*tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt

                 tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
                    ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,b0_+j)

                 CALL gradient1(tmp2,ixmin1,ixmin2,ixmax1,ixmax2,m,tmpC)

                 wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew&
                    (ixImin1:ixImax1,ixImin2:ixImax2,e_)+sB&
                    *tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt

              ENDIF


           ENDDO
        ENDIF
     ENDDO
  ENDDO




  RETURN
END SUBROUTINE addsource_visc

!=============================================================================
SUBROUTINE setnu(w,iw,idim,ixmin1,ixmin2,ixmax1,ixmax2,nuR,nuL)

  ! Set the viscosity coefficient nu within ixO based on w(ixI).

  INCLUDE 'vacdef.f'

  INTEGER:: iximin1,iximin2,iximax1,iximax2
  DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
  DOUBLE PRECISION:: d1R(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1),d1L(ixGlo1:ixGhi1&
     +1,ixGlo2:ixGhi2+1)
  DOUBLE PRECISION:: d3R(ixGlo1:ixGhi1+1,ixGlo2:ixGhi2+1),d3L(ixGlo1:ixGhi1&
     +1,ixGlo2:ixGhi2+1)
  DOUBLE PRECISION:: md3R(ixGlo1:ixGhi1,ixGlo2:ixGhi2),md3L(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)
  DOUBLE PRECISION:: md1R(ixGlo1:ixGhi1,ixGlo2:ixGhi2),md1L(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)
  DOUBLE PRECISION:: nuR(ixGlo1:ixGhi1,ixGlo2:ixGhi2),nuL(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)

  DOUBLE PRECISION:: c_tot, c_hyp,cmax(ixGlo1:ixGhi1,ixGlo2:ixGhi2),&
      tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
  INTEGER:: ixmin1,ixmin2,ixmax1,ixmax2,idim, iw
  INTEGER:: kxmin1,kxmin2,kxmax1,kxmax2,jxmin1,jxmin2,jxmax1,jxmax2,hxmin1,&
     hxmin2,hxmax1,hxmax2,gxmin1,gxmin2,gxmax1,gxmax2,ixFFmin1,ixFFmin2,&
     ixFFmax1,ixFFmax2,jxFFmin1,jxFFmin2,jxFFmax1,jxFFmax2,hxFFmin1,hxFFmin2,&
     hxFFmax1,hxFFmax2
  INTEGER:: ix_1,ix_2,ix_3

  INTEGER:: ixFlo1,ixFlo2,ixFhi1,ixFhi2,ixFmin1,ixFmin2,ixFmax1,ixFmax2,&
     ixYlo1,ixYlo2,ixYhi1,ixYhi2

  LOGICAL:: new_cmax

  DOUBLE PRECISION:: tmp_nuI(ixGlo1:ixGhi1+2,ixGlo2:ixGhi2+2)

  INTEGER:: k,iwc

  INTEGER:: ix,ixe



  !----------------------------------------------------------------------------

  new_cmax=.TRUE.

  CALL getcmax(new_cmax,w,ixmin1,ixmin2,ixmax1,ixmax2,idim,cmax)
  c_tot=MAXVAL(cmax(ixmin1:ixmax1,ixmin2:ixmax2))



  !---------------------------------------------
  ! Set HyperVis coefficients here:
  !---------------------------------------------

  c_hyp=0.4d0 ! 1.4d0 ! 0.6

  IF (iw.EQ.b1_.OR.iw.EQ.b2_) c_hyp=0.02d0 ! 2d0

  IF (iw .EQ. rho_) c_hyp=0.02d0 !5d0
!  IF (iw .EQ. rho_) c_hyp=0.045d0 !5d0  used for orszag-tang test



  IF (iw .EQ. 173) c_hyp=0.02d0 !2d0


  !---------------------------------------------


  IF (iw .NE. 173) THEN
     tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw)
     IF (iw.EQ.m1_.OR.iw.EQ.m2_) tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)&
        =w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,iw)/(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
        rho_)+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,rhob_))
  ENDIF

  IF (iw .EQ. 173) tmp_nu(ixGlo1:ixGhi1,ixGlo2:ixGhi2)=w(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2,e_)-half*((w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,b1_)**2&
     +w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,b2_)**2)+(w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,&
     m1_)**2+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,m2_)**2)/(w(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2,rho_)+w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,rhob_)))


  ixYlo1=ixmin1-2;ixYlo2=ixmin2-2;ixYhi1=ixmax1+2;ixYhi2=ixmax2+2;

  ixFlo1=ixYlo1+1;ixFlo2=ixYlo2+1;ixFhi1=ixYhi1+1;ixFhi2=ixYhi2+1;

  tmp_nuI(ixFlo1:ixFhi1,ixFlo2:ixFhi2)=tmp_nu(ixYlo1:ixYhi1,ixYlo2:ixYhi2)





  IF (iw .EQ. 173) THEN
     iwc=e_
  ELSE
     iwc=iw
  ENDIF

  DO k=0,1  !left-right bc

     IF (typeB(iwc,2*idim-1+k) .NE. 'mpi') THEN
        IF (upperB(2*idim-1+k)) THEN

           SELECT CASE(idim)
                 CASE(1)
              tmp_nuI(ixFhi1+1,ixFlo2:ixFhi2)=tmp_nuI(ixFhi1-5,ixFlo2:ixFhi2)
                 CASE(2)
              tmp_nuI(ixFlo1:ixFhi1,ixFhi2+1)=tmp_nuI(ixFlo1:ixFhi1,ixFhi2-5)

           END SELECT

        ELSE

           SELECT CASE(idim)
                 CASE(1)
              tmp_nuI(ixFlo1-1,ixFlo2:ixFhi2)=tmp_nuI(ixFlo1+5,ixFlo2:ixFhi2)
                 CASE(2)
              tmp_nuI(ixFlo1:ixFhi1,ixFlo2-1)=tmp_nuI(ixFlo1:ixFhi1,ixFlo2+5)

           END SELECT

        ENDIF
     ENDIF

  ENDDO

  ixFmin1=ixFlo1+1;ixFmin2=ixFlo2+1;ixFmax1=ixFhi1-1;ixFmax2=ixFhi2-1;

  kxmin1=ixFmin1+2*kr(idim,1);kxmin2=ixFmin2+2*kr(idim,2)
  kxmax1=ixFmax1+2*kr(idim,1);kxmax2=ixFmax2+2*kr(idim,2); !5:66
  jxmin1=ixFmin1+kr(idim,1);jxmin2=ixFmin2+kr(idim,2)
  jxmax1=ixFmax1+kr(idim,1);jxmax2=ixFmax2+kr(idim,2); !4:65
  hxmin1=ixFmin1-kr(idim,1);hxmin2=ixFmin2-kr(idim,2)
  hxmax1=ixFmax1-kr(idim,1);hxmax2=ixFmax2-kr(idim,2); !2:63
  gxmin1=ixFmin1-2*kr(idim,1);gxmin2=ixFmin2-2*kr(idim,2)
  gxmax1=ixFmax1-2*kr(idim,1);gxmax2=ixFmax2-2*kr(idim,2); !1:62

  ixFFmin1=ixFlo1;ixFFmin2=ixFlo2;ixFFmax1=ixFhi1;ixFFmax2=ixFhi2;   !2:65
  jxFFmin1=ixFlo1+kr(idim,1);jxFFmin2=ixFlo2+kr(idim,2)
  jxFFmax1=ixFhi1+kr(idim,1);jxFFmax2=ixFhi2+kr(idim,2); !3:66
  hxFFmin1=ixFlo1-kr(idim,1);hxFFmin2=ixFlo2-kr(idim,2)
  hxFFmax1=ixFhi1-kr(idim,1);hxFFmax2=ixFhi2-kr(idim,2); !1:64

  d3R(ixFmin1:ixFmax1,ixFmin2:ixFmax2)=ABS(3.d0*(tmp_nuI(jxmin1:jxmax1,&
     jxmin2:jxmax2)-tmp_nuI(ixFmin1:ixFmax1,ixFmin2:ixFmax2))&
     -(tmp_nuI(kxmin1:kxmax1,kxmin2:kxmax2)-tmp_nuI(hxmin1:hxmax1,&
     hxmin2:hxmax2))) !3:64
  d1R(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2)=ABS(tmp_nuI(jxFFmin1:jxFFmax1,&
     jxFFmin2:jxFFmax2)-tmp_nuI(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2)) !2:65

  DO ix_1=ixmin1,ixmax1
  DO ix_2=ixmin2,ixmax2    !3:62  +1=4:63

  md3R(ix_1,ix_2)=MAXVAL(d3R(ix_1+1-kr(idim,1):ix_1+1+kr(idim,1),ix_2&
     +1-kr(idim,2):ix_2+1+kr(idim,2)))
  md1R(ix_1,ix_2)=MAXVAL(d1R(ix_1+1-2*kr(idim,1):ix_1+1+2*kr(idim,1),ix_2&
     +1-2*kr(idim,2):ix_2+1+2*kr(idim,2)))

  ENDDO
  ENDDO

  WHERE (md1R(ixmin1:ixmax1,ixmin2:ixmax2).GT.0.d0)
     nuR(ixmin1:ixmax1,ixmin2:ixmax2)=c_tot*c_hyp*md3R(ixmin1:ixmax1,&
        ixmin2:ixmax2)/md1R(ixmin1:ixmax1,ixmin2:ixmax2)*dx(ixmin1:ixmax1,&
        ixmin2:ixmax2,idim)
  ELSEWHERE
     nuR(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0
  END WHERE

  maxviscoef=MAX(MAXVAL(nuR(ixmin1:ixmax1,ixmin2:ixmax2)), maxviscoef)


  !************

  d3L(ixFmin1:ixFmax1,ixFmin2:ixFmax2)=ABS(3.d0*(tmp_nuI(ixFmin1:ixFmax1,&
     ixFmin2:ixFmax2)-tmp_nuI(hxmin1:hxmax1,hxmin2:hxmax2))&
     -(tmp_nuI(jxmin1:jxmax1,jxmin2:jxmax2)-tmp_nuI(gxmin1:gxmax1,&
     gxmin2:gxmax2)))
  d1L(ixFFmin1:ixFFmax1,ixFFmin2:ixFFmax2)=ABS(tmp_nuI(ixFFmin1:ixFFmax1,&
     ixFFmin2:ixFFmax2)-tmp_nuI(hxFFmin1:hxFFmax1,hxFFmin2:hxFFmax2))

  DO ix_1=ixmin1,ixmax1
  DO ix_2=ixmin2,ixmax2

  md3L(ix_1,ix_2)=MAXVAL(d3L(ix_1+1-kr(idim,1):ix_1+1+kr(idim,1),ix_2&
     +1-kr(idim,2):ix_2+1+kr(idim,2)))
  md1L(ix_1,ix_2)=MAXVAL(d1L(ix_1+1-2*kr(idim,1):ix_1+1+2*kr(idim,1),ix_2&
     +1-2*kr(idim,2):ix_2+1+2*kr(idim,2)))

  ENDDO
  ENDDO

  WHERE (md1L(ixmin1:ixmax1,ixmin2:ixmax2).GT.0.d0)
     nuL(ixmin1:ixmax1,ixmin2:ixmax2)=c_tot*c_hyp*md3L(ixmin1:ixmax1,&
        ixmin2:ixmax2)/md1L(ixmin1:ixmax1,ixmin2:ixmax2)*dx(ixmin1:ixmax1,&
        ixmin2:ixmax2,idim)
  ELSEWHERE
     nuL(ixmin1:ixmax1,ixmin2:ixmax2)=0.d0
  END WHERE

  maxviscoef=MAX(MAXVAL(nuL(ixmin1:ixmax1,ixmin2:ixmax2)), maxviscoef)



  RETURN
END SUBROUTINE setnu


!=============================================================================
!=============================================================================
SUBROUTINE setnushk(w,ixmin1,ixmin2,ixmax1,ixmax2,nushk)

  INCLUDE 'vacdef.f'

  !double precision:: w(ixG^T,nw),tmp2(ixG^T),nushk(ixG^T,ndim)
  DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),nushk(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2,ndim)

  DOUBLE PRECISION:: c_shk

  DOUBLE PRECISION:: tmp3(ixGlo1:ixGhi1,ixGlo2:ixGhi2)

  INTEGER:: ixmin1,ixmin2,ixmax1,ixmax2,idim, iw,i

  INTEGER:: ix_1,ix_2

  DO idim=1,ndim
     nushk(ixmin1:ixmax1,ixmin2:ixmax2,idim)=0.d0
  ENDDO


  !--------------------------------------------------
  ! Comment this out and NOT USE A DAMN GOTO!
  !--------------------------------------------------
!!$c_shk=0.5d0
!!$
!!$tmp3(ix^S)=0.d0
!!$
!!$!**************************BEGIN shock viscosity*******************************
!!$      do idim=1,ndim
!!$         tmp(ix^S)=w(ix^S,m0_+idim)/(w(ix^S,rho_)+w(ix^S,rhob_))
!!$         call gradient1(tmp,ix^L,idim,tmp2)
!!$         tmp3(ix^S)=tmp3(ix^S)+tmp2(ix^S)
!!$       enddo
!!$      do idim=1,ndim
!!$        nushk(ix^S,idim)=tmp3(ix^S)*(dx(ix^S,idim)**2.d0)*c_shk
!!$	WHERE (tmp3(ix^S) .ge. 0.d0)
!!$!	  nushk(ix^S,idim)=0.d0
!!$	END WHERE
!!$	nushk(ix^S,idim)=abs(nushk(ix^S,idim))
!!$      enddo
!!$!****************************END shock viscosity*******************************


  RETURN
END SUBROUTINE setnushk



!=============================================================================
SUBROUTINE getdt_visc(w,ixmin1,ixmin2,ixmax1,ixmax2)

  ! Check diffusion time limit for dt < dtdiffpar * dx**2 / (nu/rho)

  ! Based on Hirsch volume 2, p.631, eq.23.2.17

  INCLUDE 'vacdef.f'

  DOUBLE PRECISION:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),dtdiff_visc
  INTEGER:: ixmin1,ixmin2,ixmax1,ixmax2,idim, ix_1,ix_2

  INTEGER:: aa

  ! For spatially varying nu you need a common nu array
  DOUBLE PRECISION::tmpdt(ixGlo1:ixGhi1,ixGlo2:ixGhi2), nuL(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2),nuR(ixGlo1:ixGhi1,ixGlo2:ixGhi2), nushk(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2,ndim)
  COMMON/visc/nuL
  COMMON/visc/nuR
  !-----------------------------------------------------------------------------

  CALL setnushk(w,ixmin1,ixmin2,ixmax1,ixmax2,nushk)

  dtdiffpar=0.25d0

  DO idim=1,ndim
     tmpdt(ixmin1:ixmax1,ixmin2:ixmax2)=(maxviscoef+nushk(ixmin1:ixmax1,&
        ixmin2:ixmax2,idim)) !/(w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))   ! 1/dt
     dtdiff_visc=dtdiffpar/MAXVAL(tmpdt(ixmin1:ixmax1,ixmin2:ixmax2)&
        /(dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)**2))

     dt=MIN(dt,dtdiff_visc)
  END DO

  maxviscoef=0.d0

  RETURN
END SUBROUTINE getdt_visc


!***** 2-point central finite difference gradient******

SUBROUTINE gradient1(q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)
  INCLUDE 'vacdef.f'
  INTEGER:: ixmin1,ixmin2,ixmax1,ixmax2,idim
  DOUBLE PRECISION:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)
  INTEGER:: hxmin1,hxmin2,hxmax1,hxmax2,kxmin1,kxmin2,kxmax1,kxmax2
  INTEGER:: minx11,minx12,maxx11,maxx12,k
  !-----------------------------------------------------------------------------

  hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmax1=ixmax1-kr(idim,1)
  hxmax2=ixmax2-kr(idim,2);
  kxmin1=ixmin1+kr(idim,1);kxmin2=ixmin2+kr(idim,2);kxmax1=ixmax1+kr(idim,1)
  kxmax2=ixmax2+kr(idim,2);
  gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(q(kxmin1:kxmax1,kxmin2:kxmax2)&
     -q(hxmin1:hxmax1,hxmin2:hxmax2))/dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)/two

  minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);
  maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);

  DO k=0,1  !left-right bc
     IF (typeB(1,2*idim-1+k) .NE. 'periodic') THEN
        IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
           IF (upperB(2*idim-1+k)) THEN
              SELECT CASE(idim)
                    CASE(1)
                 gradq(ixmax1,ixmin2:ixmax2)=0.d0
                 gradq(maxx11,ixmin2:ixmax2)=0.d0
                    CASE(2)
                 gradq(ixmin1:ixmax1,ixmax2)=0.d0
                 gradq(ixmin1:ixmax1,maxx12)=0.d0

              END SELECT
           ELSE
              SELECT CASE(idim)
                    CASE(1)
                 gradq(ixmin1,ixmin2:ixmax2)=0.d0
                 gradq(minx11,ixmin2:ixmax2)=0.d0
                    CASE(2)
                 gradq(ixmin1:ixmax1,ixmin2)=0.d0
                 gradq(ixmin1:ixmax1,minx12)=0.d0

              END SELECT
           ENDIF
        ENDIF
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE gradient1

!=============================================================================


!*****left upwind forward 2-point non-central finite difference gradient******

SUBROUTINE gradient1L(q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)
  INCLUDE 'vacdef.f'
  INTEGER:: ixmin1,ixmin2,ixmax1,ixmax2,idim
  DOUBLE PRECISION:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)
  INTEGER:: hxmin1,hxmin2,hxmax1,hxmax2
  INTEGER:: minx11,minx12,maxx11,maxx12,k
  !-----------------------------------------------------------------------------

  hxmin1=ixmin1-kr(idim,1);hxmin2=ixmin2-kr(idim,2);hxmax1=ixmax1-kr(idim,1)
  hxmax2=ixmax2-kr(idim,2);
  gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(q(ixmin1:ixmax1,ixmin2:ixmax2)&
     -q(hxmin1:hxmax1,hxmin2:hxmax2))/dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)

  minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);
  maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);

  DO k=0,1  !left-right bc
     IF (typeB(1,2*idim-1+k) .NE. 'periodic') THEN
        IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
           IF (upperB(2*idim-1+k)) THEN
              SELECT CASE(idim)
                    CASE(1)
                 gradq(ixmax1,ixmin2:ixmax2)=0.d0
                 gradq(maxx11,ixmin2:ixmax2)=0.d0
                    CASE(2)
                 gradq(ixmin1:ixmax1,ixmax2)=0.d0
                 gradq(ixmin1:ixmax1,maxx12)=0.d0

              END SELECT
           ELSE
              SELECT CASE(idim)
                    CASE(1)
                 gradq(ixmin1,ixmin2:ixmax2)=0.d0
                 gradq(minx11,ixmin2:ixmax2)=0.d0
                    CASE(2)
                 gradq(ixmin1:ixmax1,ixmin2)=0.d0
                 gradq(ixmin1:ixmax1,minx12)=0.d0

              END SELECT
           ENDIF
        ENDIF
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE gradient1L

!=============================================================================

!*****right upwind forward 2-point non-central finite difference gradient*****

SUBROUTINE gradient1R(q,ixmin1,ixmin2,ixmax1,ixmax2,idim,gradq)
  INCLUDE 'vacdef.f'
  INTEGER:: ixmin1,ixmin2,ixmax1,ixmax2,idim
  DOUBLE PRECISION:: q(ixGlo1:ixGhi1,ixGlo2:ixGhi2),gradq(ixGlo1:ixGhi1,&
     ixGlo2:ixGhi2)
  INTEGER:: hxmin1,hxmin2,hxmax1,hxmax2
  INTEGER:: minx11,minx12,maxx11,maxx12,k
  !-----------------------------------------------------------------------------

  hxmin1=ixmin1+kr(idim,1);hxmin2=ixmin2+kr(idim,2);hxmax1=ixmax1+kr(idim,1)
  hxmax2=ixmax2+kr(idim,2);
  gradq(ixmin1:ixmax1,ixmin2:ixmax2)=(q(hxmin1:hxmax1,hxmin2:hxmax2)&
     -q(ixmin1:ixmax1,ixmin2:ixmax2))/dx(ixmin1:ixmax1,ixmin2:ixmax2,idim)

  minx11=ixmin1+kr(idim,1);minx12=ixmin2+kr(idim,2);
  maxx11=ixmax1-kr(idim,1);maxx12=ixmax2-kr(idim,2);

  DO k=0,1  !left-right bc
     IF (typeB(1,2*idim-1+k) .NE. 'periodic') THEN
        IF (typeB(1,2*idim-1+k) .NE. 'mpi') THEN
           IF (upperB(2*idim-1+k)) THEN
              SELECT CASE(idim)
                    CASE(1)
                 gradq(ixmax1,ixmin2:ixmax2)=0.d0
                 gradq(maxx11,ixmin2:ixmax2)=0.d0
                    CASE(2)
                 gradq(ixmin1:ixmax1,ixmax2)=0.d0
                 gradq(ixmin1:ixmax1,maxx12)=0.d0

              END SELECT
           ELSE
              SELECT CASE(idim)
                    CASE(1)
                 gradq(ixmin1,ixmin2:ixmax2)=0.d0
                 gradq(minx11,ixmin2:ixmax2)=0.d0
                    CASE(2)
                 gradq(ixmin1:ixmax1,ixmin2)=0.d0
                 gradq(ixmin1:ixmax1,minx12)=0.d0

              END SELECT
           ENDIF
        ENDIF
     ENDIF
  ENDDO


  RETURN
END SUBROUTINE gradient1R
!INCLUDE:vacusr.diffusion.t
!=============================================================================
subroutine specialini(ixmin1,ixmin2,ixmax1,ixmax2,w)

include 'vacdef.f'

integer:: ixmin1,ixmin2,ixmax1,ixmax2
integer:: ix_1,ix_2,ix_3
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)

double precision:: rhoin,xcent1,xcent2,radius
double precision:: inirho,iniene
double precision:: onemor,inix,ddx
double precision:: p_1,p_2

integer:: iii_,iix_1,info,i,j
!double precision:: pi,comi,eneu,sum,mode,bmax,l
double precision:: comi,eneu,sum,mode,bmax,l

character*79 atmfilename

double precision:: p1,p2,rho1,rho2,v1,v2,T1,T2, b1_1,b1_2,b2_1,b2_2
double precision:: vx,vz,bx,bz

!-----------------------------------------------------------------------------

p1=5.d0/12.d0/Pi
rho1=25.d0/36.d0/Pi

vx=1.d0
vz=1.d0

bx=1.d0/2.0/sqrt(Pi)
bz=1.d0/2.0/sqrt(Pi)

  w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=0.d0
  w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)=rho1

  w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)=-vx*rho1*sin(2.d0*Pi*x(ixmin1:ixmax1,&
     ixmin2:ixmax2,1)/x(ixmax1,10,1))

  w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)=vz*rho1*sin(2.d0*Pi*x(ixmin1:ixmax1,&
     ixmin2:ixmax2,2)/x(10,ixmax2,2))

  w(ixmin1:ixmax1,ixmin2:ixmax2,b2_)=-bx*sin(2.d0*Pi*x(ixmin1:ixmax1,&
     ixmin2:ixmax2,1)/x(ixmax1,10,1))
  w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)=bz*sin(4.d0*Pi*x(ixmin1:ixmax1,&
     ixmin2:ixmax2,2)/x(10,ixmax2,2))

  w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=0.d0
  w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)=p1/(eqpar(gamma_)-1.d0)


!  w(ix^S,e_)=half*(^C&w(ix^S,m^C_)**2.d0+)/(w(ix^S,rho_)+w(ix^S,rhob_))

  w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)=w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)&
     +half*((w(ixmin1:ixmax1,ixmin2:ixmax2,b1_)**2.d0+w(ixmin1:ixmax1,&
     ixmin2:ixmax2,b2_)**2.d0))+half*(w(ixmin1:ixmax1,ixmin2:ixmax2,m1_)&
     **2.d0+w(ixmin1:ixmax1,ixmin2:ixmax2,m2_)**2.d0)/(w(ixmin1:ixmax1,&
     ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))


  w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)=w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)
  w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)=0.d0

  w(ixmin1:ixmax1,ixmin2:ixmax2,e_)=w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)
  w(ixmin1:ixmax1,ixmin2:ixmax2,eb_)=0.d0

  w(ixmin1:ixmax1,ixmin2:ixmax2,bg2_)=0.d0
  w(ixmin1:ixmax1,ixmin2:ixmax2,bg1_)=0.d0

!  w(ix^S,bg2_)=w(ix^S,b2_)
!  w(ix^S,b2_)=0.d0
!  w(ix^S,bg1_)=w(ix^S,b1_)
!  w(ix^S,b1_)=0.d0


return
end


!=============================================================================
subroutine specialsource(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
   ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)


include 'vacdef.f'

integer:: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,&
   iws(niw_)
double precision:: qdt,qtC,qt,wCT(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw),&
   w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
double precision:: fdt,fdthalf2

double precision:: pre(ixGlo1:ixGhi1,ixGlo2:ixGhi2),tem(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),kapr(ixGlo1:ixGhi1,ixGlo2:ixGhi2),so(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2),flux(ixGlo1:ixGhi1,ixGlo2:ixGhi2)
double precision:: tau(ixGlo1:ixGhi1,ixGlo2:ixGhi2),ine(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

double precision:: preg(ixGlo1:ixGhi1,ixGlo2:ixGhi2),pret(ixGlo1:ixGhi1,&
   ixGlo2:ixGhi2)

integer:: rix_1,i,j
double precision:: mol_0, rrr_

double precision:: fsokr,avgflux

integer:: iw,iiw,iix_1



!point rast source**********
double precision:: rad,sigma1,sigma2,Q0,qt0,xc1,xc2,zzs,qmin,qra
double precision:: rfc,tdep,sdep
double precision:: xs(100),zs(100),ts(100),qs(100)
double precision:: tlast,rn
integer:: ns
!logical:: filexist
integer:: singl
!***************

integer:: ix_1,ix_2,idim, ixmin1,ixmin2,ixmax1,ixmax2

!*****************
!double precision:: t01,t02,a1,a2,s1,s2,sf,xc1,xc2,rad,rfc,sdep,tdep,sigma2


!-----------------------------------------------------------------------------



eqpar(nu_)=1.d0
!eqpar(nu_)=0.d0



if(abs(eqpar(nu_))>smalldouble)&
   call addsource_visc(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,iws,qtC,wCT,qt,w)

!call addsource_grav(qdt,ixI^L,ixO^L,iws,qtC,wCT,qt,w)


write(*,*) '***time=',qt



end



!=============================================================================
subroutine specialbound(qt,ixmin1,ixmin2,ixmax1,ixmax2,iw,iB,w)
include 'vacdef.f'



integer:: iwmin,iwmax,idimmin,idimmax
double precision:: qt,w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,1:nw)
integer:: ix,ix1,ix2,ixe,ixf,ixmin1,ixmin2,ixmax1,ixmax2,ixpairmin1,&
   ixpairmin2,ixpairmax1,ixpairmax2,idim,iw,iB
integer:: iwv,jdim

integer:: Ns,i,j
double precision:: ki

integer:: ix_1,ix_2

double precision:: tmpp1,tmpp2



return
end

!=============================================================================
subroutine getdt_special(w,ixmin1,ixmin2,ixmax1,ixmax2)

! If the Coriolis force is made very strong it may require time step limiting,
! but this is not implemented here.

include 'vacdef.f'
double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2
!-----------------------------------------------------------------------------

!call getdt_diff(w,ix^L)


if(abs(eqpar(nu_))>smalldouble)&
   call getdt_visc(w,ixmin1,ixmin2,ixmax1,ixmax2)

!call getdt_res(w,ix^L)

!call getdt_grav(w,ix^L)

return
end


subroutine specialeta(w,ixmin1,ixmin2,ixmax1,ixmax2,idirmin)

include 'vacdef.f'

double precision:: w(ixGlo1:ixGhi1,ixGlo2:ixGhi2,nw)
integer:: ixmin1,ixmin2,ixmax1,ixmax2,idirmin
!-----------------------------------------------------------------------------

stop 'specialeta is not defined'

end

