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
    use commvar,  only : reynolds,lrestart,mach,ia,ja,ka,im,jm,km,roinf, const6
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
    real(8) :: divU,Omega,Psi,Phi
    real(8) :: divUdivU,dPsidPsi,dPhidPhi,OmegadivU,OmegaOmega, &
              m11m11,m22m22,m11m22,m12m21,m12m12,m21m21,du2
    real(8) :: dPsidPsidivU,dPhidPhidivU,divUdivUdivU,OmegaOmegadivU, &
              m11m11m11,m22m22m22,m11m11m22,m22m22m11,m11m12m12,&
              m22m21m21,m22m21m12,m11m12m21,m11m21m21,m22m12m12,du3
    real(8) :: kolmloc,mfpath,machrms,csavg,niuavg,miuavg,rhoavg,rho2nd,w2drho,ensdissp
    real(8) :: miu,niu,rhoprime,u1,u2,v2,cs,du11,du12,du21,du22,s11,s12,s22,div,omegaz,dPhi,dPsi,dissloc
    real(8) :: ens,macht,skewness,ufluc,Kollength,Taylength,ReTay,&
                Intlength,ReInt,EnsKol,EnsLarge,ReEnsLar,EnsMicro,&
                ReEnsMic ! mfpath = mean free path
    !
    logical,save :: linit=.true.
    integer,save :: hand_a,hand_b,hand_c,hand_d,hand_e
    !
    R = 8.31446261815324d0
    if(ka .ne. 0) then
      !
      stop "Error! Not 2D"
      ! 2D part
    endif
    !
    !!!!! File initialization
    if(linit) then
      !
      if(lio) then
        call listinit(filename='log/stat2d_ener.dat',handle=hand_a, &
            firstline='ns ti urms Kall Kw Krho pavg eavg ptheta rhocotheta uutheta rhouugrad dissps disspd dissp')
        call listinit(filename='log/stat2d_2nd.dat',handle=hand_b, &
                      firstline='ns ti th2 ps2 f2 oth o2 A1111 A2222 A1122 A1221 A1212 A2121')
        call listinit(filename='log/stat2d_3rd.dat',handle=hand_c, &
                      firstline='ns ti ps2th f2th th3 o2th A111 A222 A112 A221 A111212 A222121 A222112 A111221 A112121 A221212')
        call listinit(filename='log/stat2d_di.dat',handle=hand_d, &
                      firstline='ns ti th o ps f kolmloc mfpath marms csavg nuav muav roav rho2nd w2drho ensdis')
        call listinit(filename='log/stat2d_scale.dat',handle=hand_e, &
                      firstline='ns ti ens macht skew ufluc Kol Tay ReTay Int ReInt EnsKol EnsLar ReLar EnsMicro ReMic')
      endif
      !
      linit=.false.
      !
    endif
    !
    rsamples=dble(ia*ja)
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
    Omega=0.d0
    Psi=0.d0
    Phi=0.d0
    divU=0.d0
    !
    ! Velgrad 2nd order 
    OmegaOmega=0.d0
    OmegadivU=0.d0
    divUdivU=0.d0
    dPsidPsi=0.d0
    dPhidPhi=0.d0
    m11m11=0.d0
    m22m22=0.d0
    m11m22=0.d0
    m12m21=0.d0
    m12m12=0.d0
    m21m21=0.d0
    du2=0.d0     ! Total second order(m11m11+m22m22)
    !
    ! Velgrad 3rd order
    OmegaOmegadivU=0.d0
    divUdivUdivU=0.d0
    dPsidPsidivU=0.d0
    dPhidPhidivU=0.d0
    m11m11m11=0.d0
    m22m22m22=0.d0
    m11m11m22=0.d0
    m22m22m11=0.d0
    m11m12m12=0.d0
    m22m12m12=0.d0
    m11m21m21=0.d0
    m22m21m21=0.d0
    m11m12m21=0.d0
    m22m21m12=0.d0
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
    ensdissp=0.d0
    !!!!! Vorticity Gradient Calculation
    k=0
    do j=0,jm
    do i=0,im
      vorbis(i,j,k)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
    enddo
    enddo
    call dataswap(vorbis)
    dvor=grad(vorbis)
    !
    !!!!! Point calculation
    k=0
    do j=1,jm
    do i=1,im
      !
      !! Values within this loop
      miu = miucal(tmp(i,j,k))/Reynolds
      niu = miu/rho(i,j,k)
      rhoprime = rho(i,j,k) - roinf
      u1 = vel(i,j,k,1)
      u2 = vel(i,j,k,2)
      v2 = u1**2+u2**2
      cs = sos(tmp(i,j,k))
      !
      du11=dvel(i,j,k,1,1)
      du12=dvel(i,j,k,1,2)
      du21=dvel(i,j,k,2,1)
      du22=dvel(i,j,k,2,2)
      !
      s11=du11
      s12=0.5d0*(du12+du21)
      s22=du22
      !
      div    = du11+du22
      omegaz = du21-du12
      dPsi   = du12+du21
      dPhi   = du11-du22
      !
      dissloc   = 2.d0*miu*(s11**2+s22**2+2.d0*(s12**2)-num1d3*div**2)
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
      rhouugrad  = rhouugrad  + rhoprime * (u1*u1*du11 + u1*u2*du12 + &
                                            u1*u2*du21 + u2*u2*du22)
      dissps    = dissps + 1.d0*(omegaz**2) ! The multiplication by miuavg will be outside the loop
      disspd    = disspd + 4.d0/3.d0*(div**2) ! The multiplication by miuavg will be outside the loop
      dissp     = dissp  + dissloc
      !
      ! Velgrad 1st order
      divU  = divU  + div
      Omega = Omega + omegaz
      Psi   = Psi   + dPsi
      Phi   = Phi   + dPhi
      !
      ! Velgrad 2nd order 
      OmegaOmega = OmegaOmega + omegaz*omegaz
      OmegadivU  = OmegadivU  + omegaz*div
      divUdivU   = divUdivU   + div*div
      dPsidPsi   = dPsidPsi   + dPsi*dPsi
      dPhidPhi   = dPhidPhi   + dPhi*dPhi
      m11m11     = m11m11     + du11*du11
      m22m22     = m22m22     + du22*du22
      m11m22     = m11m22     + du11*du22
      m12m21     = m12m21     + du12*du21
      m12m12     = m12m12     + du12*du12
      m21m21     = m21m21     + du21*du21
      du2        = du2        + du11**2+du22**2
      !
      ! Velgrad 3rd order 
      OmegaOmegadivU = OmegaOmegadivU + omegaz*omegaz*div
      divUdivUdivU   = divUdivUdivU   + div*div*div
      dPsidPsidivU   = dPsidPsidivU   + dPsi*dPsi*div
      dPhidPhidivU   = dPhidPhidivU   + dPhi*dPhi*div
      m11m11m11   = m11m11m11 + du11*du11*du11
      m22m22m22   = m22m22m22 + du22*du22*du22
      m11m11m22   = m11m11m22 + du11*du11*du22
      m22m22m11   = m22m22m11 + du22*du22*du11
      m11m12m12   = m11m12m12 + du11*du12*du12
      m22m12m12   = m22m12m12 + du22*du12*du12
      m11m21m21   = m11m21m21 + du11*du21*du21
      m22m21m21   = m22m21m21 + du22*du21*du21
      m11m12m21   = m11m12m21 + du11*du12*du21
      m22m21m12   = m22m21m12 + du22*du21*du12
      du3         = du3 + (du11**3+du22**3)
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
      w2drho = w2drho + (omegaz*omegaz)/rho(i,j,k)
      !
      ensdissp = ensdissp + niu * (dvor(i,j,k,1)**2+dvor(i,j,k,2)**2)
      !
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
    Omega  = psum(Omega)/rsamples
    Psi  = psum(Psi)/rsamples
    Phi  = psum(Phi)/rsamples
    !
    OmegaOmega  = psum(OmegaOmega)/rsamples
    OmegadivU   = psum(OmegadivU)/rsamples
    divUdivU    = psum(divUdivU)/rsamples
    dPsidPsi    = psum(dPsidPsi)/rsamples
    dPhidPhi    = psum(dPhidPhi)/rsamples
    m11m11  = psum(m11m11)/rsamples
    m22m22  = psum(m22m22)/rsamples
    m11m22  = psum(m11m22)/rsamples
    m12m21  = psum(m12m21)/rsamples
    m12m12  = psum(m12m12)/rsamples
    m21m21  = psum(m21m21)/rsamples
    du2     = psum(du2)/rsamples
    !
    OmegaOmegadivU  = psum(OmegaOmegadivU)/rsamples
    divUdivUdivU  = psum(divUdivUdivU)/rsamples
    dPsidPsidivU  = psum(dPsidPsidivU)/rsamples
    dPhidPhidivU  = psum(dPhidPhidivU)/rsamples
    m11m11m11  = psum(m11m11m11)/rsamples
    m22m22m22  = psum(m22m22m22)/rsamples
    m11m11m22  = psum(m11m11m22)/rsamples
    m22m22m11  = psum(m22m22m11)/rsamples
    m11m12m12  = psum(m11m12m12)/rsamples
    m22m12m12  = psum(m22m12m12)/rsamples
    m11m21m21  = psum(m11m21m21)/rsamples
    m22m21m21  = psum(m22m21m21)/rsamples
    m11m12m21  = psum(m11m12m21)/rsamples
    m22m21m12  = psum(m22m21m12)/rsamples
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
    ensdissp = psum(ensdissp)/rsamples
    !
    ! Other statistics
    macht         = urms/csavg
    skewness      = du3/2.d0/sqrt((du2/2.d0)**3)
    !
    ! Scales and Reynolds number
    ufluc     = urms/sqrt(2.d0) ! Prepare work
    ! Energy based
    Kollength = sqrt(sqrt(niuavg**3/dissp))
    Taylength = ufluc/sqrt(du2/2.d0)
    ReTay     = ufluc * Taylength / niuavg
    Intlength = sqrt(Kall**3)/dissp
    ReInt     = Kall**2/(niuavg*dissp)
    ! Enstrophy based based on Herring 1974
    EnsKol    = (niu**3/ensdissp)**(1.d0/6.d0)
    EnsLarge  = urms/sqrt(2.d0)/ensdissp**(1.d0/3.d0)
    ReEnsLar  = 0.5d0*urms**2/(niuavg*ensdissp**(1.d0/3.d0))
    EnsMicro  = sqrt(niuavg*ens/ensdissp)
    ReEnsMic  = sqrt(ens**3)/(ensdissp)
    !
    !
    if(lio) then 
      call listwrite(hand_a,urms,Kall,Kw,Krho,pavg,eavg,ptheta,&
                    rhocotheta,uutheta,rhouugrad,dissps,disspd,dissp)
      call listwrite(hand_b,divUdivU,dPsidPsi,dPhidPhi,OmegadivU,&
                      OmegaOmega, m11m11,m22m22,m11m22,m12m21,m12m12,m21m21)
      call listwrite(hand_c,dPsidPsidivU,dPhidPhidivU,divUdivUdivU, &
                     OmegaOmegadivU,m11m11m11,m22m22m22,m11m11m22,m22m22m11,&
                     m11m12m12,m22m21m21,m22m21m12,m11m12m21,m11m21m21,m22m12m12)
      call listwrite(hand_d,divU,Omega,Psi,Phi,kolmloc,mfpath,&
                    machrms,csavg,niuavg,miuavg,rhoavg,rho2nd,w2drho,ensdissp)
      call listwrite(hand_e,ens,macht,skewness,ufluc,Kollength,Taylength,&
                    ReTay,Intlength,ReInt,EnsKol,EnsLarge,ReEnsLar,EnsMicro,&
                    ReEnsMic)
    endif
    !
  end subroutine udf_stalist
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_stalist.                            |
  !+-------------------------------------------------------------------+
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is to add vortical fluctuations to initial field  | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
  ! subroutine addvortex(xc,yc,radius,amp)
    !
    ! use commvar,  only: im,jm,km,ndims,roinf,uinf
    ! use parallel, only: lio
    ! use commarray,only: x,vel,rho,prs,spc,tmp,q
    ! use fludyna,  only: thermal
    ! !
    ! ! local data
    ! real(8),intent(in) :: xc,yc,radius,amp
    ! !
    ! integer :: i,j,k
    ! real(8) :: var1,radi2,cvor
    ! !
    ! cvor=amp*uinf*radius
    ! !
    ! do k=0,km
    ! do j=0,jm
    ! do i=0,im
    !   radi2=((x(i,j,k,1)-xc)**2+(x(i,j,k,2)-yc)**2)/radius/radius
    !   var1=cvor/radius/radius*exp(-0.5d0*radi2)
    !   !
    !   vel(i,j,k,1)=vel(i,j,k,1)-var1*(x(i,j,k,2)-yc)
    !   if(ndims>=2) vel(i,j,k,2)=vel(i,j,k,2)+var1*(x(i,j,k,1)-xc)
    !   if(ndims==3) vel(i,j,k,3)=0.d0
    !   prs(i,j,k)  =prs(i,j,k)-0.5d0*roinf*cvor*cvor/radi2/radi2*exp(-radi2)
    !   !
    !   tmp(i,j,k)=thermal(density=rho(i,j,k),pressure=prs(i,j,k),species=spc(i,j,k,:))
    !   !
    ! enddo
    ! enddo
    ! enddo
    !
  ! end subroutine addvortex
  !+-------------------------------------------------------------------+
  !| The end of the subroutine addvortex.                              |
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
                         lforce,nstep, forcenum,lhyper,roinf
    use parallel, only : lio,psum,bcast
    use commarray,only : rho,tmp,vel,qrhs,x,jacob,forcep,forcek
    use utility,  only : listinit,listwrite
    use constdef, only : pi
    use statistic,only : diss_rate_cal
    !
    logical,save :: linit=.true.
    integer,save :: hand_force
    integer, allocatable, dimension(:), save :: hand_forcea
    ! Random iniforce generation
    integer :: NumTheta, n, i,j,k,t
    real(8) :: theta
    real(8) :: power,rsamples,Tpower
    real(8) :: kappaT,kappaF
    real(8), save :: dissp=0.d0
    real(8), allocatable,dimension(:) :: alphas, alphad
    character(len=4) :: forcename
    !
    !
    !
    if(lforce) then
      !
      if(.not. allocated(hand_forcea))then
        allocate(hand_forcea(1:forcenum))
      else
        if(size(hand_forcea) .ne. forcenum) stop 'Error in hand_forcea @ udf_src'
      endif
      !
      if(linit) then
        !
        if(lio) then
          !
          call listinit(filename='log/forcestat.dat',handle =hand_force, firstline='nstep time rkstep forcekT power disspation')
          !
          do t=1,forcenum
            write(forcename,'(i4.4)')t
            call listinit(filename='log/forcestat'//forcename//'.dat',handle =hand_forcea(t), &
             firstline='nstep time kappaT kappaF alphad')
          enddo
          !
        endif
        !
        linit = .false.
        !
      endif
      !
      dissp = diss_rate_cal()
      !
      if(rkstep==1)then
        !
        allocate(alphas(1:forcenum),alphad(1:forcenum))
        !
        call udf_generate_force(dissp,alphas,alphad)
        !
        if(lio) then 
          !
          do t=1,forcenum
            call listwrite(hand_forcea(t),real(forcek(t),8),alphas(t),alphad(t))
          enddo
        endif
      endif
      !
      !
      ! Calculate power
      power = 0.0d0
      Tpower = 0.0d0
      !
      if(ndims == 2)then
        do j=1,jm
        do i=1,im
          power  = power  + roinf*(forcep(i,j,0,1)*vel(i,j,0,1) + forcep(i,j,0,2)*vel(i,j,0,2))
          Tpower = Tpower + tmp(i,j,0)**4
        enddo
        enddo
        !
        rsamples=dble(ia*ja)
        !
      endif
      !
      power = psum(power)/rsamples
      Tpower = psum(Tpower)/rsamples
      kappaF = dissp/power
      kappaT = dissp/Tpower
      !
      if(lio) call listwrite(hand_force,dble(rkstep),kappaT,kappaF,dissp)
      !
      ! Add in qrhs and calculate power
      do k=0,km
      do j=0,jm
      do i=0,im
        !
        !
        qrhs(i,j,k,2)=qrhs(i,j,k,2)+kappaF*roinf*forcep(i,j,k,1)*jacob(i,j,k)
        qrhs(i,j,k,3)=qrhs(i,j,k,3)+kappaF*roinf*forcep(i,j,k,2)*jacob(i,j,k)
        qrhs(i,j,k,4)=qrhs(i,j,k,4)+kappaF*roinf*forcep(i,j,k,3)*jacob(i,j,k)
        qrhs(i,j,k,5)=qrhs(i,j,k,5)+kappaF*roinf*(forcep(i,j,k,1)*vel(i,j,k,1) + &
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
      deallocate(alphas,alphad)
    endif
    !
  end subroutine udf_src
  !+-------------------------------------------------------------------+
  !| The end of the subroutine udf_src.                                |
  !+-------------------------------------------------------------------+
  !
  subroutine udf_generate_force(dissp,alphas, alphad)
    !
    use commvar, only: ndims
    real(8), intent(in) ::  dissp
    real(8), allocatable, dimension(:), intent(out) :: alphas, alphad
    !   
    if(ndims == 2) then
      call udf_generate_force_2D(dissp,alphas, alphad)
    else
      print *, "ndims = ", ndims
      stop "Not implemented error! udf_generate_force unvalid ndims"
    endif
    !
  end subroutine udf_generate_force
    !
  subroutine udf_generate_force_2D(dissipation,alphas, alphad)
    !
    use, intrinsic :: iso_c_binding
    use commvar,        only : forcenum,hypervisk,hypervismiu,im,jm,ia,ja,deltat,lhyper
    use commarray,      only : vel, forcep,forcek,forcespes,forcesped
    use fftwlink,       only : jmf, alloc_local, iafftw, jmfftw, fftw_grid_fence, fftw_fence_grid, jafftw, j0f
    use parallel,       only : MPI_COMM_WORLD, psum, mpiright, mpiup, mpitag, mpileft, mpidown
    use udf_tool,       only : GenerateWave, kint
    use mpi
    use statistic,      only : diss_rate_cal
    include 'fftw3-mpi.f03'
    !
    real(8), intent(in) ::  dissipation
    real(8), allocatable, dimension(:), intent(out):: alphas, alphad
    real(8), allocatable, dimension(:,:) :: localvel1t, localvel2t, force1t, force2t
    real(8), allocatable, dimension(:,:) :: fftvel1, fftvel2, fftforce1, fftforce2
    type(C_PTR) :: forward_plan, backward_plan, c_u1spe, c_u2spe, c_fftforce1, c_fftforce2
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1spe,u2spe
    real(8), allocatable, dimension(:,:) :: k1,k2
    complex(8), allocatable, dimension(:,:) :: usspe,udspe,u1s,u2s,u1d,u2d
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: fftforce1c,fftforce2c
    real(8) ::  dk, kk
    integer :: i,j,ierr,allkmax,kOrdinal,t
    real(8), allocatable, dimension(:) :: Ed, Es
    real(8), allocatable, dimension(:) :: sendim,recvim, sendjm,recvjm
    integer :: status(mpi_status_size) 
    !
    dk = 1.d0
    !
    allkmax=ceiling(real(sqrt(2.d0)/3*min(ia,ja))/dk)
    !
    if(.not. allocated(alphas)) allocate(alphas(1:forcenum))
    if(.not. allocated(alphad)) allocate(alphad(1:forcenum))
    allocate(localvel1t(1:jm,1:im),localvel2t(1:jm,1:im))
    allocate(force1t(1:jm,1:im),force2t(1:jm,1:im))
    allocate(fftvel1(1:ia,1:jmf),fftvel2(1:ia,1:jmf))
    allocate(fftforce1(1:ia,1:jmf),fftforce2(1:ia,1:jmf))
    !
    do j=1,jm
    do i=1,im
      !
      localvel1t(j,i)=vel(i,j,0,1)
      localvel2t(j,i)=vel(i,j,0,2)
      !
    enddo
    enddo
    !
    !
    call fftw_grid_fence(localvel1t,fftvel1)
    call fftw_grid_fence(localvel2t,fftvel2)
    !
    ! Begin FFTW
      !
    c_u1spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1spe, u1spe, [iafftw,jmfftw])
    c_u2spe = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2spe, u2spe, [iafftw,jmfftw])
    !
    !!!! Do S-C decomposition
    c_fftforce1 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_fftforce1, fftforce1c, [iafftw,jmfftw])
    c_fftforce2 = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_fftforce2, fftforce2c, [iafftw,jmfftw])
    !
    ! planning
    forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
    backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1spe,u1spe, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
    !
    !
    do j=1,jmf
    do i=1,ia
        !
        u1spe(i,j)=CMPLX(fftvel1(i,j),0.d0,C_INTPTR_T);
        u2spe(i,j)=CMPLX(fftvel2(i,j),0.d0,C_INTPTR_T);
        !
    end do
    end do
    !
    !
    !!!! Do 2d FFT
    call fftw_mpi_execute_dft(forward_plan,u1spe,u1spe)
    call fftw_mpi_execute_dft(forward_plan,u2spe,u2spe)
    !
    do j=1,jmf
    do i=1,ia
        !
        u1spe(i,j)=u1spe(i,j)/(1.d0*ia*ja)
        u2spe(i,j)=u2spe(i,j)/(1.d0*ia*ja)
        !
    end do
    end do
    !
    ! Wavenumber calculation
    allocate(k1(1:ia,1:jmf),k2(1:ia,1:jmf))
    call GenerateWave(ia,jmf,ia,ja,j0f,k1,k2)
    !
    allocate(usspe(1:ia,1:jmf),udspe(1:ia,1:jmf))
    allocate(u1s(1:ia,1:jmf),u1d(1:ia,1:jmf),u2s(1:ia,1:jmf),u2d(1:ia,1:jmf))
    allocate(Es(0:allkmax),Ed(0:allkmax))
    Ed = 0.d0
    Es = 0.d0
    !
    do j=1,jmf
    do i=1,ia
        kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
        if(kk > 1.d-5)then
          usspe(i,j) = u1spe(i,j)*k2(i,j)/kk - u2spe(i,j)*k1(i,j)/kk
          udspe(i,j) = u1spe(i,j)*k1(i,j)/kk + u2spe(i,j)*k2(i,j)/kk
          u1d(i,j)=  udspe(i,j)*k1(i,j)/kk
          u2d(i,j)=  udspe(i,j)*k2(i,j)/kk
          u1s(i,j)=  usspe(i,j)*k2(i,j)/kk 
          u2s(i,j)= -usspe(i,j)*k1(i,j)/kk
        else
          usspe(i,j) = 0
          udspe(i,j) = 0
          u1d(i,j) = 0
          u2d(i,j) = 0
          u1s(i,j) = 0
          u2s(i,j) = 0
        endif
        kOrdinal = kint(kk,dk,2,1.d0)
        if(kOrdinal <= allkmax)then
          Es(kOrdinal) = Es(kOrdinal) + usspe(i,j)*conjg(usspe(i,j))
          Ed(kOrdinal) = Ed(kOrdinal) + udspe(i,j)*conjg(udspe(i,j))
        endif
        !
    end do
    end do
    !
    !
    do i=1,allkmax
      Es(i) = psum(Es(i))
      Ed(i) = psum(Ed(i))
    enddo
    !
    do t=1,forcenum
      alphas(t) = min(max(forcespes(t)*dissipation/Es(forcek(t)), 0.d0),10.d0)
      alphad(t) = min(max(forcesped(t)*dissipation/Ed(forcek(t)), 0.d0),10.d0)
    enddo
    !
    !
    fftforce1c = 0.d0
    fftforce2c = 0.d0
    do j=1,jmf
    do i=1,ia
      kk=dsqrt(k1(i,j)**2+k2(i,j)**2)
      do t=1,forcenum
        if(kint(kk,dk,2,1.d0)==forcek(t))then
          fftforce1c(i,j) = fftforce1c(i,j) + alphas(t) * u1s(i,j) + alphad(t) * u1d(i,j)
          fftforce2c(i,j) = fftforce2c(i,j) + alphas(t) * u2s(i,j) + alphad(t) * u2d(i,j)
        endif
      enddo
      !
    enddo
    enddo
    !
    call fftw_mpi_execute_dft(backward_plan,fftforce1c,fftforce1c)
    call fftw_mpi_execute_dft(backward_plan,fftforce2c,fftforce2c)
    !
    do j=1,jmf
    do i=1,ia
      fftforce1(i,j) = real(fftforce1c(i,j))
      fftforce2(i,j) = real(fftforce2c(i,j))
    enddo
    enddo
    !
    call fftw_fence_grid(fftforce1,force1t)
    call fftw_fence_grid(fftforce2,force2t)
    !
    forcep = 0.d0
    !
    !
    do j=1,jm
    do i=1,im
      forcep(i,j,0,1) = force1t(j,i)
      forcep(i,j,0,2) = force2t(j,i)
    enddo
    enddo
    !
    call mpi_sendrecv(forcep(im,jm,0,1),1,mpi_real8,mpiright,mpitag, &
                      forcep(0,jm,0,1),1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(forcep(im,jm,0,2),1,mpi_real8,mpiright,mpitag, &
                      forcep(0,jm,0,2),1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(forcep(im,jm,0,1),1,mpi_real8,mpiup,mpitag,    &
                      forcep(im,0,0,1),1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    call mpi_sendrecv(forcep(im,jm,0,2),1,mpi_real8,mpiup,mpitag,    &
                      forcep(im,0,0,2),1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    allocate(sendim(0:im),recvim(0:im),sendjm(0:jm),recvjm(0:jm))
    !
    sendim(0:im) = forcep(0:im,jm,0,1)
    call mpi_sendrecv(sendim,im+1,mpi_real8,mpiup,mpitag,     &
                      recvim,im+1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    forcep(0:im,0,0,1) = recvim(0:im)
    mpitag=mpitag+1
    !
    sendim(0:im) = forcep(0:im,jm,0,2)
    call mpi_sendrecv(sendim,im+1,mpi_real8,mpiup,mpitag,     &
                      recvim,im+1,mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    forcep(0:im,0,0,2) = recvim(0:im)
    !
    sendjm(0:jm) = forcep(im,0:jm,0,1)
    call mpi_sendrecv(sendjm,jm+1,mpi_real8,mpiright,mpitag,  &
                      recvjm,jm+1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    forcep(0,0:jm,0,1) = recvjm(0:jm)
    !
    sendjm(0:jm) = forcep(im,0:jm,0,2)
    call mpi_sendrecv(sendjm,jm+1,mpi_real8,mpiright,mpitag,  &
                      recvjm,jm+1,mpi_real8,mpileft,mpitag,   &
                      mpi_comm_world,status,ierr)
    forcep(0,0:jm,0,2) = recvjm(0:jm)
    mpitag=mpitag+1
    !
    call fftw_destroy_plan(forward_plan)
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1spe)
    call fftw_free(c_u2spe)
    call fftw_free(c_fftforce1)
    call fftw_free(c_fftforce2)
    !
    deallocate(localvel1t, localvel2t, force1t, force2t)
    deallocate(fftvel1, fftvel2, fftforce1, fftforce2)
    deallocate(k1,k2,usspe,udspe, u1d, u1s, u2d, u2s)
    deallocate(Es,Ed)
    deallocate(sendim,sendjm,recvim,recvjm)
    !
  end subroutine udf_generate_force_2D
  !+-------------------------------------------------------------------+
  !| This subroutine is to defined an output by a user.                | 
  !+-------------------------------------------------------------------+
  !| CHANGE RECORD                                                     |
  !| -------------                                                     |
  !| 18-Aug-2023: created by Jian Fang @ Daresbury                     |
  !+-------------------------------------------------------------------+
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
