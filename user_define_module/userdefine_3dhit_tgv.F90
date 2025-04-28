!+---------------------------------------------------------------------+
!| This module contains user defined subroutines to interfere  program |
!+---------------------------------------------------------------------+
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!| 18-08-2023  | Created by J. Fang                                    |
!+---------------------------------------------------------------------+
module userdefine
  !
  implicit none
  !
  real(8) :: hsource
  !
  contains
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to set flow environment, such as, incoming     |
  !| free stream variables.                                            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_setflowenv
    !
!     use commvar,  only: roinf,uinf,vinf,winf,pinf,tinf,spcinf,num_species
!     use fludyna,  only: thermal
!     !
! #ifdef COMB
!     use thermchem,only : tranco,spcindex,mixture,convertxiyi
!     use cantera 
!     !
!     real(8) :: specr(num_species)
!     ! 
!     specr(:)=0.d0
!     specr(spcindex('H2'))=0.0173
!     specr(spcindex('O2'))=0.2289
!     specr(spcindex('N2'))=1.d0-sum(specr)
!     !
!     ! pinf=5.d0*pinf
!     uinf=0.97d0
!     vinf=0.d0
!     winf=0.d0
!     tinf=300.d0
!     spcinf(:)=specr(:)
!     roinf=thermal(pressure=pinf,temperature=tinf,species=spcinf(:))
!     !
! #endif
    !
  end subroutine udf_setflowenv
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_setflowenv.                         |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate fluctuations for inflow            |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 05-Oct-2023: Created by by Jian Fang @ Daresbury                  |
  !+-------------------------------------------------------------------+
  subroutine udf_inflow_fluc(umean,uinst)
    !
    use commvar, only : jm,km
    !
    real(8),intent(in) ::  umean(0:jm,1:3)  ! inflow mean velocity
    real(8),intent(out) :: uinst(0:jm,0:km,1:3)  ! velocity with fluctuations
    !
  end subroutine udf_inflow_fluc
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_inflow_fluc.                        |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to initialise flowfield by a user              |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 25-May-2023: Created by Yifan Xu @ Peking University              |
  !| 18-Aug-2023: Rename and relocated by Jian Fang @ Daresbury        |
  !+-------------------------------------------------------------------+
  subroutine udf_flowinit
    !
!     use commvar,  only: im,jm,km,ndims,roinf,uinf,nondimen,xmax,pinf,  &
!                         ia,num_species
!     use commarray,only: x,vel,rho,prs,spc,tmp,q
!     use parallel, only: lio
!     use fludyna,  only: thermal
!     !
! #ifdef COMB
!     !
!     use thermchem,only : tranco,spcindex,mixture,convertxiyi
!     use cantera 
!     !
!     ! local data
!     integer :: i,j,k
!     real(8) ::  xc,yc,zc,tmpr,tmpp,xloc,xwid,specr(num_species),  &
!       specp(num_species),arg,prgvar,masflx,specx(num_species)
!     real(8) :: pthick
!     !
!     tmpr=300.d0
!     xloc=3.d0*xmax/4.d0
!     xwid=xmax/(12.d0*5.3d0*2.d0)
!     !
!     !reactants
!     specr(:)=0.d0
!     specr(spcindex('H2'))=0.0173
!     specr(spcindex('O2'))=0.2289
!     specr(spcindex('N2'))=1.d0-sum(specr)
!     !
!     !products
!     tmpp=1814.32d0
!     !
!     ! pthick=1.d-4
!     !
!     do k=0,km
!     do j=0,jm
!     do i=0,im
!       !
!       xc=x(i,j,k,1)
!       !
!       !prgvar=0.5d0*(1.d0+tanh(10.d0*(xc-xloc)/xloc))
!       ! if(xc-xloc<xwid*0.5d0*1.2d0) then 
!       !   prgvar=0.d0
!       !   if(xc-xloc>xwid*0.5d0) &
!       !   prgvar=1.d0-(xc-xloc-(xwid*0.5d0))/(xwid*0.5d0*0.2d0)
!       ! else
!       !   prgvar=1.d0
!       ! endif
!       !
!       prgvar=1.d0*exp(-0.5d0*((xc-xloc)/xwid)**2)
!       !
!       spc(i,j,k,:)=specr(:)
!       !
!       vel(i,j,k,1)=uinf
!       !
!       vel(i,j,k,2)=0.d0
!       vel(i,j,k,3)=0.d0
!       !
!       tmp(i,j,k)=tmpr+prgvar*(tmpp-tmpr)
!       !
!       prs(i,j,k)=pinf
!       !
!       rho(i,j,k)=thermal(pressure=prs(i,j,k),temperature=tmp(i,j,k), &
!                           species=spc(i,j,k,:))
!     enddo
!     enddo
!     enddo
!     !
!     !
!     if(lio)  write(*,'(A,I1,A)')'  ** HIT flame initialised.'
!     !
! #endif
    !
  end subroutine udf_flowinit
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_flowinit.                           |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to generate grid.                              | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 23-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_grid
  end subroutine udf_grid
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_grid.                               |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to list something during a computation.        | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !| 12-Oct-2023: adapted by Chensheng Luo @ Beihang University        |
  !+-------------------------------------------------------------------+
  subroutine udf_stalist
    !
    use constdef
    use commvar,  only : reynolds,lrestart,mach,ia,ja,ka,im,jm,km,roinf,const6
    use commarray,only : vel,rho,tmp,dvel,q,vorbis,dvor,prs
    use fludyna,  only : miucal,sos
    use comsolver,only : solvrinit,grad
    use utility,  only : listinit,listwrite
    use parallel, only : dataswap,psum,lio,pmax,pmin
    use constdef, only : pi,num1d3
    !
    integer :: i,j,k,ns
    real(8) :: rsamples,R
    real(8) :: urms,Kall,Kw,Krho,pavg,eavg,ptheta,rhocotheta,uutheta,rhouugrad,dissps,disspd,dissp
    real(8) :: divU
    real(8) :: divUdivU,traceSS,OmegaOmega,du2,m11m11,m22m22,m33m33, &
              m12m12,m13m13,m21m21,m23m23,m31m31,m32m32
    real(8) :: traceSSth,wsw,s3,divUdivUdivU,OmegaOmegadivU,du3, &
              m11m11m11,m22m22m22,m33m33m33,m11m12m12,m11m12m21,m11m11m22
    real(8) :: kolmloc,mfpath,machrms,csavg,niuavg,miuavg,rhoavg,rho2nd,w2drho
    real(8) :: miu,niu,rhoprime,u1,u2,u3,v2,cs,&
              du11,du12,du13,du21,du22,du23,du31,du32,du33,&
              s11,s22,s33,s12,s13,s23,div,omegax,omegay,omegaz,dPhi,dPsi,dissloc
    real(8) :: ens,macht,skewness,ufluc,Kollength,Taylength,ReTay,&
                Intlength,ReInt
    !
    logical,save :: linit=.true.
    integer,save :: hand_a,hand_b,hand_c,hand_d,hand_e
    !
    R = 8.31446261815324d0
    if(ka .eq. 0) then
      !
      stop "Error! Not 3D"
      ! 3D part
    endif
    !
    !!!!! File initialization
    if(linit) then
      !
      if(lio) then
        call listinit(filename='log/stat2d_ener.dat',handle=hand_a, &
            firstline='ns ti urms Kall Kw Krho pavg eavg ptheta rhocotheta uutheta rhouugrad dissps disspd dissp')
        call listinit(filename='log/stat2d_2nd.dat',handle=hand_b, &
                      firstline='ns ti th2 o2 trss 2ndx 2ndy 2ndz 12^2 13^2 21^2 23^2 31^2 32^2')
        call listinit(filename='log/stat2d_3rd.dat',handle=hand_c, &
                      firstline='ns ti th3 o2th trssth wsw s3 A111 A222 A333 A111212 A111221 A111122')
        call listinit(filename='log/stat2d_di.dat',handle=hand_d, &
                      firstline='ns ti th kolmloc mfpath marms csavg nuav muav roav rho2nd w2drho')
        call listinit(filename='log/stat2d_scale.dat',handle=hand_e, &
                      firstline='ns ti ens macht skew ufluc Kol Tay ReTay Int ReInt')
      endif
      !
      linit=.false.
      !
    endif
    !
    rsamples=dble(ia*ja*ka)
    !
    !!!!! Quantities initialization
    !
    ! u-2nd-order
    urms=0.d0
    !
    ! Energy terms
    Kall=0.d0      ! Kinetic energy
    Kw=0.d0        ! Weakly compressible kinetic energy
    Krho=0.d0      ! Density fluctuation kinetic energy
    pavg=0.d0      ! Average pressure
    eavg=0.d0      ! Average internal energy
    !
    ! Transfer terms
    ptheta=0.d0    ! Transfer to p
    rhocotheta=0.d0    ! Transfer to p(considered as 0)
    uutheta=0.d0   ! Kw-Krho transfer
    rhouugrad=0.d0 ! Kw-Krho transfer ignored ny w.c.a.
    dissps=0.d0    ! Solenoidal dissipation (with density)
    disspd=0.d0    ! dilatational dissipation (with density)
    dissp=0.d0     ! total dissipation
    !
    ! Velgrad 1st order
    divU=0.d0
    !
    ! Velgrad 2nd order 
    OmegaOmega=0.d0
    divUdivU= 0.d0
    traceSS = 0.d0
    m11m11  = 0.d0
    m22m22  = 0.d0
    m33m33  = 0.d0
    m12m12  = 0.d0
    m13m13  = 0.d0
    m21m21  = 0.d0
    m23m23  = 0.d0
    m31m31  = 0.d0
    m32m32  = 0.d0
    du2=0.d0     ! Total second order(m11m11+m22m22)
    !
    ! Velgrad 3rd order
    OmegaOmegadivU=0.d0
    divUdivUdivU=0.d0
    traceSSth=0.d0
    wsw = 0.d0
    s3 = 0.d0
    m11m11m11 = 0.d0
    m22m22m22 = 0.d0
    m33m33m33 = 0.d0
    m11m12m12 = 0.d0
    m11m12m21 = 0.d0
    m11m11m22 = 0.d0
    du3=0.d0      ! Total third order(m11m11m11+m22m22m22)
    !
    ! Lengths
    kolmloc=2*pi  ! min local Kolmogorov length
    mfpath=0.d0   ! max mean free path
    !
    ! Sound speed related
    machrms=0.d0  ! root-mean-square mach number (turbulent mach number)
    csavg=0.d0    ! average speed of sound
    !
    ! Density related
    niuavg=0.d0  ! average miu/rho
    miuavg=0.d0  ! average miu
    rhoavg=0.d0  ! average rho
    rho2nd=0.d0  ! average rho**2
    w2drho=0.d0  ! average w**2/rho
    !
    !
    !!!!! Point calculation
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      !! Values within this loop
      miu = miucal(tmp(i,j,k))/Reynolds
      niu = miu/rho(i,j,k)
      rhoprime = rho(i,j,k) - roinf
      u1 = vel(i,j,k,1)
      u2 = vel(i,j,k,2)
      u3 = vel(i,j,k,3)
      v2 = u1**2+u2**2+u3**2
      cs = sos(tmp(i,j,k))
      !
      du11=dvel(i,j,k,1,1)
      du12=dvel(i,j,k,1,2)
      du13=dvel(i,j,k,1,3)
      du21=dvel(i,j,k,2,1)
      du22=dvel(i,j,k,2,2)
      du23=dvel(i,j,k,2,3)
      du31=dvel(i,j,k,3,1)
      du32=dvel(i,j,k,3,2)
      du33=dvel(i,j,k,3,3)
      !
      s11=du11
      s22=du22
      s33=du33
      s12=0.5d0*(du12+du21)
      s13=0.5d0*(du13+du31)
      s23=0.5d0*(du23+du32)
      !
      div    = du11+du22+du33
      omegax = du32-du23
      omegay = du13-du31
      omegaz = du21-du12
      dPsi   = du11-du22
      dPhi   = du11-du33
      !
      dissloc   = 2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)- num1d3*div**2) 
      !! Statistic quantities
      !
      urms = urms+v2
      !
      ! Energy terms
      Kall = Kall + 0.5d0* rho(i,j,k) *v2
      Kw   = Kw   + 0.5d0* roinf      *v2
      Krho = Krho + 0.5d0* rhoprime   *v2
      pavg = pavg + prs(i,j,k)
      eavg = eavg + prs(i,j,k)*const6
      !
      ! Transfer terms
      ptheta     = ptheta     + prs(i,j,k)*div
      rhocotheta = rhocotheta + rho(i,j,k)*(cs**2)*div
      uutheta    = uutheta    + 0.5d0*roinf*v2*div
      rhouugrad  = rhouugrad  + rhoprime * (u1*u1*du11 + u1*u2*du12 + u1*u3*du13 + &
                                            u2*u1*du21 + u2*u2*du22 + u2*u3*du23 + &
                                            u3*u1*du31 + u3*u2*du32 + u3*u3*du33)
      dissps    = dissps + 1.d0*(omegax**2 + omegay**2 + omegaz**2) ! The multiplication by miuavg will be outside the loop
      disspd    = disspd + 4.d0/3.d0*(div**2) ! The multiplication by miuavg will be outside the loop
      dissp     = dissp  + dissloc
      !
      ! Velgrad 1st order
      divU  = divU  + div
      !
      s12=du12+du21
      s13=du13+du31
      s23=du23+du32
      ! Velgrad 2nd order 
      OmegaOmega = OmegaOmega + omegax**2 + omegay**2 + omegaz**2
      divUdivU   = divUdivU   + div*div
      traceSS = traceSS + s12*s12/2.d0 + s13*s13/2.d0 + s23*s23/2.d0 + 2.d0*dPsi*dPsi/3.d0  &
                        - 2.d0*dPsi*dPhi/3.d0 + 2.d0*dPhi*dPhi/3.d0
      m11m11  = m11m11 + du11**2
      m22m22  = m22m22 + du22**2
      m33m33  = m33m33 + du33**2
      m12m12  = m12m12 + du12*du12
      m13m13  = m13m13 + du13*du13
      m21m21  = m21m21 + du21*du21
      m23m23  = m23m23 + du23*du23
      m31m31  = m31m31 + du31*du31
      m32m32  = m32m32 + du32*du32
      du2 = du2 + du11**2 + du22**2 + du33**2
      !
      ! Velgrad 3rd order 
      OmegaOmegadivU = OmegaOmegadivU + (omegax**2  + omegay**2 + omegaz**2)*div
      divUdivUdivU   = divUdivUdivU   + div*div*div
      traceSSth = traceSSth + (s12*s12/2.d0 + s13*s13/2.d0 + s23*s23/2.d0 + 2.d0*dPsi*dPsi/3.d0 &
                              - 2*dPsi*dPhi/3.d0 + 2*dPhi*dPhi/3.d0)*div
      wsw = wsw + s12*omegay*omegax/4.d0 + s13*omegaz*omegax/4.d0 + s23*omegay*omegaz/4.d0   &
                + dPsi*omegax*omegax/12.d0 - dPsi*omegay*omegay/6.d0  + dPsi*omegaz*omegaz/12.d0 &
                + dPhi*omegax*omegax/12.d0 + dPhi*omegay*omegay/12.d0 - dPhi*omegaz*omegaz/6.d0 
      s3 = s3 - s12*s12*dPsi/4.d0 + s13*s13*dPsi/2.d0 - s23*s23*dPsi/4.d0 +3.d0/4.d0*s12*s13*s23 &
              + s12*s12*dPhi/2.d0 - s13*s13*dPhi/4.d0 - s23*s23*dPhi/4.d0 &
              - 2.d0/9.d0 * dPsi**3 + dPsi**2*dPhi /3.d0 + dPsi*dPhi**2/3.d0 - 2.d0/9.d0 * dPhi**3
      m11m11m11 = m11m11m11 + du11**3
      m22m22m22 = m22m22m22 + du22**3
      m33m33m33 = m33m33m33 + du33**3
      m11m12m12 = m11m12m12 + du11*du12*du12
      m11m12m21 = m11m12m21 + du11*du12*du21
      m11m11m22 = m11m11m22 + du11*du11*du22
      du3  = du3 + (du11**3+du22**3+du33**3)
      !
      ! Lengths
      kolmloc = min(kolmloc,sqrt(sqrt((niu**3)/dissloc)))
      mfpath  = max(mfpath,2.d0*niu/0.921/sqrt(3.d0*R*tmp(i,j,k)))
      !
      machrms = machrms+v2/(cs*cs)
      csavg   = csavg+cs
      !
      niuavg = niuavg + niu
      miuavg = miuavg + miu
      rhoavg = rhoavg + rho(i,j,k)
      rho2nd = rho2nd + rho(i,j,k)*rho(i,j,k)
      w2drho = w2drho + (omegax**2 + omegay**2 + omegaz**2)/rho(i,j,k)
      !
      !
    enddo
    enddo
    enddo
    !
    !!!! Summation and average
    !
    urms = sqrt(psum(urms)/rsamples)
    !
    Kall = psum(Kall)/rsamples
    Kw   = psum(Kw)/rsamples
    Krho = psum(Krho)/rsamples
    pavg = psum(pavg)/rsamples
    eavg = psum(eavg)/rsamples
    !
    ptheta     = psum(ptheta)/rsamples
    rhocotheta = psum(rhocotheta)/rsamples
    uutheta    = psum(uutheta)/rsamples
    rhouugrad  = psum(rhouugrad)/rsamples
    dissps     = psum(dissps)/rsamples
    disspd     = psum(disspd)/rsamples
    dissp      = psum(dissp)/rsamples
    !
    divU  = psum(divU)/rsamples
    !
    OmegaOmega  = psum(OmegaOmega)/rsamples
    divUdivU    = psum(divUdivU)/rsamples
    traceSS     = psum(traceSS)/rsamples
    m11m11  = psum(m11m11)/rsamples
    m22m22  = psum(m22m22)/rsamples
    m33m33  = psum(m22m22)/rsamples
    m12m12  = psum(m12m12)/rsamples
    m13m13  = psum(m13m13)/rsamples
    m21m21  = psum(m21m21)/rsamples
    m23m23  = psum(m23m23)/rsamples
    m31m31  = psum(m31m31)/rsamples
    m32m32  = psum(m32m32)/rsamples
    du2     = psum(du2)/rsamples
    !
    OmegaOmegadivU  = psum(OmegaOmegadivU)/rsamples
    divUdivUdivU  = psum(divUdivUdivU)/rsamples
    traceSSth  = psum(traceSSth)/rsamples
    wsw = psum(wsw)/rsamples
    s3 = psum(s3)/rsamples
    m11m11m11  = psum(m11m11m11)/rsamples
    m22m22m22  = psum(m22m22m22)/rsamples
    m33m33m33  = psum(m33m33m33)/rsamples
    m11m11m22  = psum(m11m11m22)/rsamples
    m11m12m12  = psum(m11m12m12)/rsamples
    m11m12m21  = psum(m11m12m21)/rsamples
    du3        = psum(du3)/rsamples
    !
    kolmloc    = pmin(kolmloc)
    mfpath     = pmax(mfpath)
    !
    machrms = sqrt(psum(machrms)/rsamples)
    csavg   = psum(csavg)/rsamples
    !
    niuavg  = psum(niuavg)/rsamples
    miuavg  = psum(miuavg)/rsamples
    rhoavg = psum(rhoavg)/rsamples
    rho2nd = psum(rho2nd)/rsamples
    w2drho = psum(w2drho)/rsamples
    !
    dissps = dissps * miuavg
    disspd = disspd * miuavg
    !
    ens = 0.5d0 * OmegaOmega
    !
    ! Other statistics
    macht         = urms/csavg
    skewness      = du3/3.d0/sqrt((du2/3.d0)**3)
    !
    ! Scales and Reynolds number
    ufluc     = urms/sqrt(3.d0) ! Prepare work
    ! Energy based
    Kollength = sqrt(sqrt(niuavg**3/dissp))
    Taylength = ufluc/sqrt(du2/3.d0)
    ReTay     = ufluc * Taylength / niuavg
    Intlength = sqrt(Kall**3)/dissp
    ReInt     = Kall**2/(niuavg*dissp)
    !
    !
    if(lio) then 
      call listwrite(hand_a,urms,Kall,Kw,Krho,pavg,eavg,ptheta,&
                    rhocotheta,uutheta,rhouugrad,dissps,disspd,dissp)
      call listwrite(hand_b,divUdivU,OmegaOmega,traceSS,&
                    m11m11,m22m22,m33m33,m12m12,m13m13,m21m21,m23m23,m31m31,m32m32)
      call listwrite(hand_c,divUdivUdivU,OmegaOmegadivU,traceSSth,wsw,s3,&
                    m11m11m11,m22m22m22,m33m33m33,m11m12m12,m11m12m21,m11m11m22)
      call listwrite(hand_d,divU,kolmloc,mfpath,&
                    machrms,csavg,niuavg,miuavg,rhoavg,rho2nd,w2drho)
      call listwrite(hand_e,ens,macht,skewness,ufluc,Kollength,Taylength,&
                    ReTay,Intlength,ReInt)
    endif
    !
  end subroutine udf_stalist
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_stalist.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine add a source term to the rsd of the equation to   |
  !| hit flame.                                                        |
  !| a random force acting like fans to input energy at largest scale  |
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 13-06-2023: Created by Yifan Xu @ Peking University               |
  !+-------------------------------------------------------------------+
  subroutine udf_src
    !
    use commvar,  only : im,jm,km,ndims,deltat,ia,ja,ka,rkstep,xmax,ymax,zmax,&
                         lforce,nstep, forcenum,lhyper,roinf,flowtype
    use parallel, only : lio,psum,bcast
    use commarray,only : rho,tmp,vel,qrhs,x,jacob,forcep,forcek,forcespes
    use utility,  only : listinit,listwrite
    use constdef, only : pi
    use statistic,only : diss_rate_cal
    !
    logical,save :: linit=.true.
    integer,save :: hand_force
    ! Random iniforce generation
    integer :: NumTheta, n, i,j,k,t
    real(8) :: theta
    real(8) :: power,rsamples,Tpower
    real(8) :: kappaT
    real(8) :: k0,Am,xx,yy,zz
    real(8), save :: dissp=0.d0
    character(len=4) :: forcename
    !
    !
    !
    if(lforce .and. trim(flowtype)=='hittgf') then
      !
      !
      if(linit) then
        !
        if(lio) then
          !
          call listinit(filename='log/forcestat.dat',handle =hand_force, firstline='nstep time rkstep Am k0 forcekT')
          !
        endif
        !
        linit = .false.
        !
      endif
      !
      k0 = dble(forcek(1))
      Am = forcespes(1)
      do k=0,km
      do j=0,jm
      do i=0,im
        xx = x(i,j,k,1)
        yy = x(i,j,k,2)
        zz = x(i,j,k,3)
        forcep(i,j,k,1) =   Am * ( sin(k0*xx)*cos(k0*yy)*sin(k0*zz) )
        forcep(i,j,k,2) = - Am * ( cos(k0*xx)*sin(k0*yy)*cos(k0*zz) )
        forcep(i,j,k,3) = 0.d0
      enddo
      enddo
      enddo
      !
      ! Calculate power
      power = 0.0d0
      Tpower = 0.0d0
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        power  = power  + roinf*(forcep(i,j,k,1)*vel(i,j,k,1) + &
                                  forcep(i,j,k,2)*vel(i,j,k,2) + &
                                  forcep(i,j,k,3)*vel(i,j,k,3))
        Tpower = Tpower + tmp(i,j,k)**4
      enddo
      enddo
      enddo
      rsamples=dble(ia*ja*ka)
      !
      power = psum(power)/rsamples
      Tpower = psum(Tpower)/rsamples
      kappaT = power/Tpower
      !
      if(lio) call listwrite(hand_force,dble(rkstep),Am,k0,kappaT)
      !
      ! Add in qrhs and calculate power
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        !
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+roinf*forcep(i,j,k,1)*jacob(i,j,k)
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+roinf*forcep(i,j,k,2)*jacob(i,j,k)
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+roinf*forcep(i,j,k,3)*jacob(i,j,k)
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+roinf*(forcep(i,j,k,1)*vel(i,j,k,1) + &
                                            forcep(i,j,k,2)*vel(i,j,k,2) + &
                                            forcep(i,j,k,3)*vel(i,j,k,3) )*jacob(i,j,k)
        !
        !
        ! temperation dissipation
        qrhs(i,j,k,5)=qrhs(i,j,k,5)-kappaT*(tmp(i,j,k)**4)*jacob(i,j,k)
        !
      end do
      end do
      end do
      !
      !
    endif
    !
  end subroutine udf_src
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_src.                                |
  !+-------------------------------------------------------------------+
  !
  subroutine udf_write
    !
    
    !
  end subroutine udf_write
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_write.                              |
  !+-------------------------------------------------------------------+
  !
    !+-------------------------------------------------------------------+
  !| This subroutine is to manipulate data solver as one likes at the  |
  !| end of each loop.                                                 | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 30-Oct-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  subroutine udf_eom_set
    !
  end subroutine udf_eom_set
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_eom_set.                            |
  !+-------------------------------------------------------------------+
  !
end module userdefine
!+---------------------------------------------------------------------+
!| The end of the module userdefine.                                   |
!+---------------------------------------------------------------------+
