!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| initial field generating.                                           |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_hitgen
  !
  !
  use constdef
  use stlaio,  only: get_unit
  use udf_tool, only: GenerateWave
  !
  implicit none
  !
  contains
  !
  subroutine ppHitgenentrance
    !
    use cmdefne
    use parallel,        only : mpirank,bcast,mpisize,lio
    !
    ! local data
    character(len=64) :: casefolder,inputfile,outputfile,viewmode, &
                          flowfieldfile, readmode,method
    integer :: velmethod,thermomethod,filetype
    !
    !
    if(mpirank == 0) then
      call readkeyboad(readmode)
    endif
    call bcast(readmode)
    
    !
    if(trim(readmode)=='3D') then
      !
      if(mpirank == 0)then
        print *, ' ** 3D mode'
        call readkeyboad(method)
        read(method,'(i1)')velmethod
        call readkeyboad(method)
        read(method,'(i1)')thermomethod
      endif
      call bcast(velmethod)
      call bcast(thermomethod)
      !
      call hitvelgen3d(velmethod,thermomethod)
      !
    elseif(trim(readmode)=='2D') then
      !
      if(mpirank == 0)then
        print *, ' ** 2D mode'
        call readkeyboad(method)
        read(method,'(i1)')velmethod
        call readkeyboad(method)
        read(method,'(i1)')thermomethod
      endif
      call bcast(velmethod)
      call bcast(thermomethod)
      !
      call hitvelgen2d(velmethod,thermomethod)
    elseif(trim(readmode)=='p2D') then
      !
      if(mpirank == 0)then
        print *, ' ** 2D mode'
        call readkeyboad(method)
        read(method,'(i1)')thermomethod
        call readkeyboad(method)
        read(method,'(i1)')filetype
      endif
      call bcast(thermomethod)
      call bcast(filetype)
      !
      call hitgenmodifyp2d(thermomethod,filetype)
      !
    elseif(trim(readmode)=='hitstat2D') then
      !
      if(mpirank == 0)   print *, ' ** hitstat 2D'
      call hitstat2d
    elseif(trim(readmode)=='hitstat3D') then
      !
      if(mpirank == 0)   print *, ' ** hitstat 3D'
      call hitstat3d
    elseif(trim(readmode)=='scale3D') then
      !
      if(mpirank == 0) then
        print *, ' ** scale3D'
        !
        call readkeyboad(flowfieldfile)
        print*,' ** flowfieldfile command: ',flowfieldfile
      endif
      call bcast(flowfieldfile)
      !
      call scale3D(flowfieldfile)
    else
      print* ,"Readmode is not defined!", readmode
    endif
    !
  end subroutine ppHitgenentrance
  !
  !
  !
  subroutine hitvelgen3d(velmethod,thermomethod)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readgrid, readic, readinput
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          roinf,pinf,spcinf,nondimen,&
                          ickmax,iomode,icurms,icsolenoidal,icdilatational
    use bc,        only : twall
    use commarray, only : vel,rho,tmp,prs,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpi_ikgroup,mpirank, psum, pmax,&
                          mpitag,mpifront,mpiback,mpi_comm_world,status,mpi_real8
    use fludyna,   only : thermal
    use hdf5io
    use tecio
    include 'fftw3-mpi.f03'
    !
    integer, intent(in) :: velmethod,thermomethod
    integer :: i,j,k,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2,wn3, wn12, wna, var1, var2
    real(8) :: ran1, ran2,ran3,ran4,rn3
    complex(8) :: vac1, vac2, vac3, crn1, crn2, crn4
    real(8) :: amplitude, urms, uav, vav, wav
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: u1c,u2c,u3c
    real(C_DOUBLE), pointer, dimension(:,:,:) :: u1r,u2r,u3r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u3c, c_u1r, c_u2r, c_u3r
    character(len=1) :: modeio
    real(8), allocatable, dimension(:,:) :: sendkm,recvkm
    !
    modeio = 'h'
    !
    !
    call readinput
    call readic
    if(mpirank==0)  print *, " ** ia:",ia,",ja:",ja,",ka:",ka
    !
    call refcal
    if(mpirank==0)  print*, ' ** refcal done!'
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, " ** fftw_mpi initialized"
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, ' ** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, ' ** parallelini done!'
    !
    if(mpirank==0)  print*, ' ** velmethod = ', velmethod, 'thermomethod = ', thermomethod
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:km),tmp(0:(2*im),0:jm,0:km),prs(0:(2*im),0:jm,0:km))
    !
    ! Generate field
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
    call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw,kmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw,kmfftw])
    c_u3c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u3c, u3c, [imfftw,jmfftw,kmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw,kmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw,kmfftw])
    c_u3r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u3r, u3r, [2*imfftw,jmfftw,kmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_3d(kafftw,jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      if(k1(i,j,k)==0 .and. k2(i,j,k)==0) then
        u1c(i,j,k)=0.d0
        u2c(i,j,k)=0.d0
        u3c(i,j,k)=0.d0
      else
        !
        call random_number(ran1)
        call random_number(ran2)
        call random_number(ran3)
        call random_number(ran4)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        rn3 =ran3*2.d0*pi
        crn4=ran4*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the module of the wavenumber in each direction
        wn1=dble(k1(i,j,k))
        wn2=dble(k2(i,j,k))
        wn3=dble(k3(i,j,k))
        wn12=sqrt(wn1**2+wn2**2)
        wna=sqrt(wn1**2+wn2**2+wn3**2)
        !
        var1=IniEnergDis(dble(ickmax),wna,velmethod)
        var2=sqrt(var1/4.d0/pi/wna**2)
        !
        vac1=var2*cdexp(crn1)*dcos(rn3)
        vac2=var2*cdexp(crn2)*dsin(rn3)
        vac3=var2*cdexp(crn4)
        !
        if(k1(i,j,k)==0 .and. k2(i,j,k)==0) then
          u1c(i,j,k)=vac1*icsolenoidal + vac3*wn1/wna*icdilatational
          u2c(i,j,k)=vac2*icsolenoidal + vac3*wn2/wna*icdilatational
          u3c(i,j,k)=vac3*wn3/wna*icdilatational
        else
          u1c(i,j,k)=(vac1*wna*wn2+vac2*wn1*wn3)/(wna*wn12)*icsolenoidal + vac3*wn1/wna*icdilatational
          u2c(i,j,k)=(vac2*wn2*wn3-vac1*wna*wn1)/(wna*wn12)*icsolenoidal + vac3*wn2/wna*icdilatational
          u3c(i,j,k)=-vac2*wn12/wna*icsolenoidal + vac3*wn3/wna*icdilatational
        end if
        !
      endif
    enddo
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, ' ** field generated!'
    !
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u3c,u3r)
    !
    if(mpirank==0)  print*,' ** Inverse fft, project to physical space. '
    !
    im = im*2-2
    !
    vel(1:im,1:jm,1:km,1)=u1r(1:im,1:jm,1:km)
    vel(1:im,1:jm,1:km,2)=u2r(1:im,1:jm,1:km)
    vel(1:im,1:jm,1:km,3)=u3r(1:im,1:jm,1:km)
    !
    uav = 0.d0
    vav = 0.d0
    wav = 0.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      uav  = uav + vel(i,j,k,1)
      vav  = vav + vel(i,j,k,2)
      wav  = wav + vel(i,j,k,3)
      !
    end do
    end do
    end do
    uav = psum(uav) /dble(ia*ja*ka)
    vav = psum(vav) /dble(ia*ja*ka)
    wav = psum(wav) /dble(ia*ja*ka)
    !
    vel(1:im,1:jm,1:km,1) = vel(1:im,1:jm,1:km,1) - uav
    vel(1:im,1:jm,1:km,2) = vel(1:im,1:jm,1:km,2) - vav
    vel(1:im,1:jm,1:km,3) = vel(1:im,1:jm,1:km,3) - wav
    !
    if(mpirank == 0) print *, " ** Average velocity set to 0."
    !
    urms = 0.d0
    !
    do k=1,km
    do j=1,jm
    do i=1,im
      !
      urms = urms + vel(i,j,k,1)**2 + vel(i,j,k,2)**2 + vel(i,j,k,3)**2
      !
    end do
    end do
    end do
    !
    amplitude = icurms / sqrt(psum(urms)/dble(ia*ja*ka))
    !
    vel(1:im,1:jm,1:km,1) = amplitude * vel(1:im,1:jm,1:km,1)
    vel(1:im,1:jm,1:km,2) = amplitude * vel(1:im,1:jm,1:km,2)
    vel(1:im,1:jm,1:km,3) = amplitude * vel(1:im,1:jm,1:km,3)
    !
    if(mpirank == 0) print *, " ** Urms set to ", icurms
    !
    vel(0,1:jm,1:km,1) = vel(im,1:jm,1:km,1)
    vel(0,1:jm,1:km,2) = vel(im,1:jm,1:km,2)
    vel(0,1:jm,1:km,3) = vel(im,1:jm,1:km,3)
    !
    vel(0:im,0,1:km,1) = vel(0:im,jm,1:km,1)
    vel(0:im,0,1:km,2) = vel(0:im,jm,1:km,2)
    vel(0:im,0,1:km,3) = vel(0:im,jm,1:km,3)
    !
    allocate(sendkm(0:im,0:jm),recvkm(0:im,0:jm))
    !
    sendkm(0:im,0:jm) = vel(0:im,0:jm,km,1)
    call mpi_sendrecv(sendkm,(im+1)*(jm+1),mpi_real8,mpifront,mpitag,  &
                      recvkm,(im+1)*(jm+1),mpi_real8,mpiback,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    vel(0:im,0:jm,0,1) = recvkm(0:im,0:jm)
    !
    sendkm(0:im,0:jm) = vel(0:im,0:jm,km,2)
    call mpi_sendrecv(sendkm,(im+1)*(jm+1),mpi_real8,mpifront,mpitag,  &
                      recvkm,(im+1)*(jm+1),mpi_real8,mpiback,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    vel(0:im,0:jm,0,2) = recvkm(0:im,0:jm)
    !
    sendkm(0:im,0:jm) = vel(0:im,0:jm,km,3)
    call mpi_sendrecv(sendkm,(im+1)*(jm+1),mpi_real8,mpifront,mpitag,  &
                      recvkm,(im+1)*(jm+1),mpi_real8,mpiback,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    vel(0:im,0:jm,0,3) = recvkm(0:im,0:jm)
    !
    if(mpirank == 0) print *, " ** Periodic boundary condition and parallel information transfer"
    !
    rho(0:im,0:jm,0:km)  = roinf
    select case (thermomethod)
      case(1)
        if(mpirank == 0) print *, " ** Incompressible Poisson solver pressure"
        call incompressuresolve3d
      case default
        if(mpirank == 0) print *, " ** Uniform pressure"
        prs(0:im,0:jm,0:km)  = pinf
    end select
    !
    do k=0,km
    do j=0,jm
    do i=0,im
      !
      if(nondimen) then
        tmp(i,j,k)  = thermal(density=rho(i,j,k),pressure=prs(i,j,k))
      else
        spc(i,j,k,:)= spcinf(:)
        tmp(i,j,k)  = thermal(density=rho(i,j,k),pressure=prs(i,j,k),species=spc(i,j,k,:))
      endif
    enddo
    enddo
    enddo
    !
    if(mpirank == 0) print *, " ** Density and temperature allocation"
    !
    call h5io_init(trim('datin/flowini3d.h5'),mode='write')
    call h5write(var=rho(0:im,0:jm,0:km),  varname='ro', mode = modeio) 
    call h5write(var=vel(0:im,0:jm,0:km,1),varname='u1', mode = modeio)
    call h5write(var=vel(0:im,0:jm,0:km,2),varname='u2', mode = modeio)
    call h5write(var=vel(0:im,0:jm,0:km,3),varname='u3', mode = modeio)
    call h5write(var=prs(0:im,0:jm,0:km),  varname='p', mode = modeio)
    call h5write(var=tmp(0:im,0:jm,0:km),  varname='t', mode = modeio)
    call h5io_end
    !
    call fftw_destroy_plan(backward_plan)
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u3c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    call fftw_free(c_u3r)
    deallocate(sendkm,recvkm)
    deallocate(k1,k2,k3)
    deallocate(vel)
    deallocate(rho,tmp,prs)
    !
  end subroutine hitvelgen3d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen.                                 |
  !+-------------------------------------------------------------------+
  !
  !
  subroutine hitvelgen2d(velmethod,thermomethod)
    !
    use, intrinsic :: iso_c_binding
    use readwrite, only : readgrid, readic, readinput
    use fftwlink
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm,Mach,Reynolds, &
                          roinf,pinf,spcinf,nondimen,&
                          ickmax,iomode,icurms,icsolenoidal,icdilatational
    use bc,        only : twall
    use commarray, only : vel,rho,tmp,prs,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpi_ikgroup,mpirank, psum, pmax, &
                          mpitag,mpiup,mpidown,mpi_comm_world,status,mpi_real8
    use fludyna,   only : thermal
    use hdf5io
    use tecio
    use solver,    only : refcal
    include 'fftw3-mpi.f03'
    !
    integer, intent(in) :: velmethod,thermomethod
    integer :: i,j,n,clock,irandom,total_m,proc_m,m
    real(8), allocatable, dimension(:,:) :: k1,k2
    integer,allocatable :: seed(:)
    real(8) :: wn1, wn2, wna, var1, var2, ran1, ran2
    complex(8) :: vac1, vac2, crn1, crn2
    real(8) :: amplitude, urms, uav, vav
    complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1c,u2c
    real(C_DOUBLE), pointer, dimension(:,:) :: u1r,u2r
    type(C_PTR) ::  backward_plan, c_u1c, c_u2c, c_u1r, c_u2r
    real(8), allocatable, dimension(:) :: sendjm,recvjm
    !
    call readinput
    call readic
    if(mpirank==0)  print *, " ** ia:",ia,",ja:",ja
    !
    call refcal
    if(mpirank==0)  print*, ' ** refcal done!'
    !
    call fftw_mpi_init()
    if(mpirank==0)  print *, " ** fftw_mpi initialized"
    !
    call mpisizedis_half_fftw
    if(mpirank==0)  print*, ' ** mpisizedis & parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, ' ** parallelini done!'
    !
    if(mpirank==0)  print*, ' ** velmethod = ', velmethod, 'thermomethod = ', thermomethod
    !
    !
    allocate( vel(-hm:2*im+hm,-hm:jm+hm,-hm:hm,1:3) )
    allocate(rho(0:(2*im),0:jm,0:0),tmp(0:(2*im),0:jm,0:0),prs(0:(2*im),0:jm,0:0))
    if(mpirank==0)  print*, ' ** allocation finished!'
    !
    ! Generate field
    !
    !! random seed
    call random_seed(size=n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock  +  37  *  (/ (irandom  -  1, irandom = 1, n) /)
    call random_seed(put=seed)
    deallocate(seed)
    !
    !! wavenumber generation
    allocate(k1(1:im,1:jm),k2(1:im,1:jm))
    call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
    !
    !! complex speed allocation
    c_u1c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u1c, u1c, [imfftw,jmfftw])
    c_u2c = fftw_alloc_complex(alloc_local)
    call c_f_pointer(c_u2c, u2c, [imfftw,jmfftw])
    !
    c_u1r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u1r, u1r, [2*imfftw,jmfftw])
    c_u2r = fftw_alloc_complex(2*alloc_local)
    call c_f_pointer(c_u2r, u2r, [2*imfftw,jmfftw])
    !
    backward_plan = fftw_mpi_plan_dft_c2r_2d(jafftw,iafftw, u1c,u1r, MPI_COMM_WORLD,FFTW_MEASURE)
    !
    !! half spectral generation
    !
    do j=1,jm
    do i=1,im
      if(k1(i,j)==0 .and. k2(i,j)==0) then
        u1c(i,j)=0.d0
        u2c(i,j)=0.d0
      else
        call random_number(ran1)
        call random_number(ran2)
        ! ran1: random number distributied in (0,1)
        !
        crn1=ran1*2.d0*pi*(0.d0,1.d0)
        crn2=ran2*2.d0*pi*(0.d0,1.d0)
        !
        ! Calculate the modul of the wavenumber in each direction
        wn1=dble(k1(i,j))
        wn2=dble(k2(i,j))
        wna=sqrt(wn1**2+wn2**2)
        !
        var1=IniEnergDis(dble(ickmax),wna,velmethod)
        var2=sqrt(var1/2.d0/pi/wna)
        !
        vac1=var2*cdexp(crn1)
        vac2=var2*cdexp(crn2)
        !
        u1c(i,j)=vac1*wn2/wna*icsolenoidal + vac2*wn1/wna*icdilatational
        u2c(i,j)=-vac1*wn1/wna*icsolenoidal + vac2*wn2/wna*icdilatational
      end if
    enddo
    enddo
    !
    call mpi_barrier(mpi_comm_world,ierr)
    !
    !
    if(mpirank==0)  print*, ' ** Field generated!'
    !
    call fftw_mpi_execute_dft_c2r(backward_plan,u1c,u1r)
    call fftw_mpi_execute_dft_c2r(backward_plan,u2c,u2r)
    !
    if(mpirank==0)  print*,' ** Inverse fft, project to physical space. '
    !
    im = im*2-2
    !
    vel(1:im,1:jm,0,1)=u1r(1:im,1:jm)
    vel(1:im,1:jm,0,2)=u2r(1:im,1:jm)
    !
    uav = 0.d0
    vav = 0.d0
    !
    do j=1,jm
    do i=1,im
      !
      uav  = uav + vel(i,j,0,1)
      vav  = vav + vel(i,j,0,2)
      !
    end do
    end do
    !
    uav = psum(uav) /dble(ia*ja)
    vav = psum(vav) /dble(ia*ja)
    !
    vel(1:im,1:jm,0,1)=vel(1:im,1:jm,0,1) - uav
    vel(1:im,1:jm,0,2)=vel(1:im,1:jm,0,2) - vav
    !
    if(mpirank == 0) print *, " ** Average velocity set to 0."
    !
    urms = 0.d0
    !
    do j=1,jm
    do i=1,im
      !
      urms = urms + vel(i,j,0,1)**2 + vel(i,j,0,2)**2
      !
    end do
    end do
    !
    amplitude = icurms / sqrt(psum(urms)/dble(ia*ja))
    !
    vel(1:im,1:jm,0,1) = amplitude * vel(1:im,1:jm,0,1)
    vel(1:im,1:jm,0,2) = amplitude * vel(1:im,1:jm,0,2)
    !
    if(mpirank == 0) print *, " ** Urms set to ", icurms
    !
    vel(0,1:jm,0,1)=vel(im,1:jm,0,1)
    vel(0,1:jm,0,2)=vel(im,1:jm,0,2)
    !
    allocate(sendjm(0:im),recvjm(0:im))
    !
    sendjm(0:im) = vel(0:im,jm,0,1)
    call mpi_sendrecv(sendjm,(im+1),mpi_real8,mpiup,mpitag,  &
                      recvjm,(im+1),mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    vel(0:im,0,0,1) = recvjm(0:im)
    !
    sendjm(0:im) = vel(0:im,jm,0,2)
    call mpi_sendrecv(sendjm,(im+1),mpi_real8,mpiup,mpitag,  &
                      recvjm,(im+1),mpi_real8,mpidown,mpitag,   &
                      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    vel(0:im,0,0,2) = recvjm(0:im)
    !
    !
    if(mpirank == 0) print *, " ** Periodic boundary condition and parallel information transfer"
    !
    call fftw_destroy_plan(backward_plan)
    !
    vel(0:im,0:jm,0,3)= 0.d0
    rho(0:im,0:jm,0)  = roinf
    select case (thermomethod)
    case(1)
      if(mpirank == 0) print *, " ** Incompressible Poisson solver pressure"
      call incompressuresolve2d
    case default
      if(mpirank == 0) print *, " ** Uniform pressure"
      prs(0:im,0:jm,0:km)  = pinf
    end select
    do j=0,jm
    do i=0,im
      !
      if(nondimen) then
        tmp(i,j,0)  = thermal(density=rho(i,j,0),pressure=prs(i,j,0))
      else
        spc(i,j,0,:)= spcinf(:)
        tmp(i,j,0)  = thermal(density=rho(i,j,0),pressure=prs(i,j,0),species=spc(i,j,0,:))
      endif
    enddo
    enddo
    !
    if(mpirank == 0) print *, " ** Temparature and uniform density allocation"
    !
    ! Output
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
    call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
    call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
    call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
    call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    !
    call h5io_end
    !
    !
    call fftw_mpi_cleanup()
    call fftw_free(c_u1c)
    call fftw_free(c_u2c)
    call fftw_free(c_u1r)
    call fftw_free(c_u2r)
    deallocate(sendjm,recvjm)
    deallocate(k1,k2)
    deallocate(vel)
    deallocate(rho,tmp,prs)
    ! 
    !
  end subroutine hitvelgen2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitgen2d_parallel.                      |
  !+-------------------------------------------------------------------+
  !
  !
  subroutine hitgenmodifyp2d(thermomethod,filetype)
    !
    use readwrite, only : readgrid, readic, readinput
    use commvar,   only : gridfile,im,jm,km,ia,ja,ka,hm, &
                          roinf,pinf,spcinf,nondimen
    use commarray, only : vel,rho,tmp,prs,spc
    use solver,    only : refcal
    use parallel,  only : parallelini,mpirank,psum,pmax,mpisizedis,parapp, &
                          mpi_comm_world,status
    use fludyna,   only : thermal
    use hdf5io
    use tecio
    use solver,    only : refcal
    include 'fftw3-mpi.f03'
    !
    integer, intent(in) :: thermomethod,filetype
    integer :: i,j,n
    !
    call readinput
    call readic
    if(mpirank==0)  print *, " ** ia:",ia,",ja:",ja
    !
    call refcal
    if(mpirank==0)  print*, ' ** refcal done!'
    !
    if(mpirank==0)  print*, '==== Examination begins:'
    !
    call mpisizedis
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parapp
    if(mpirank==0) print*, '** parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3), rho(0:im,0:jm,0:km))
    allocate(prs(0:im,0:jm,0:km), tmp(0:im,0:jm,0:km))
    !
    call h5io_init(trim('datin/flowini2d.h5'),mode='read')
    if(filetype == 0)then
      call h5read(varname='u1', var=vel(0:im,0:jm,0,1),dir='k')
      call h5read(varname='u2', var=vel(0:im,0:jm,0,2),dir='k')
    else
      call h5read(varname='u1', var=vel(0:im,0:jm,0:0,1),mode = 'h')
      call h5read(varname='u2', var=vel(0:im,0:jm,0:0,2),mode = 'h')
    endif
    call h5io_end
    !
    !
    rho(0:im,0:jm,0)  = roinf
    select case (thermomethod)
    case(1)
      if(mpirank == 0) print *, " ** Incompressible Poisson solver pressure"
      call incompressuresolve2d
    case default
      if(mpirank == 0) print *, " ** Uniform pressure"
      prs(0:im,0:jm,0:km)  = pinf
    end select
    do j=0,jm
    do i=0,im
      !
      if(nondimen) then
        tmp(i,j,0)  = thermal(density=rho(i,j,0),pressure=prs(i,j,0))
      else
        spc(i,j,0,:)= spcinf(:)
        tmp(i,j,0)  = thermal(density=rho(i,j,0),pressure=prs(i,j,0),species=spc(i,j,0,:))
      endif
    enddo
    enddo
    !
    if(mpirank == 0) print *, " ** Temparature and uniform density allocation"
    !
    ! Output
    call h5io_init(trim('datin/flowini2d.h5'),mode='write')
    !
    call h5wa2d_r8(varname='ro',var=rho(0:im,0:jm,0),  dir='k')
    call h5wa2d_r8(varname='u1',var=vel(0:im,0:jm,0,1),dir='k')
    call h5wa2d_r8(varname='u2',var=vel(0:im,0:jm,0,2),dir='k')
    call h5wa2d_r8(varname='p', var=prs(0:im,0:jm,0),  dir='k')
    call h5wa2d_r8(varname='t', var=tmp(0:im,0:jm,0),  dir='k')
    !
    call h5io_end
    !
    !
    deallocate(vel)
    deallocate(rho,tmp,prs)
    ! 
    !
  end subroutine hitgenmodifyp2d
  !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit.         |
  !+-------------------------------------------------------------------+
  subroutine hitstat3d
    !
    use readwrite, only : readgrid, readic, readinput
    use constdef
    use commvar,   only : reynolds,mach,im,jm,km,ia,ja,ka,hm,tinf
    use comsolver, only : solvrinit, grad
    use commarray, only : vel,dvel,rho,prs,tmp,x
    use fludyna,   only : miucal
    use parallel,  only : mpisizedis, mpirank, parapp,dataswap, psum, pmax, parallelini
    use hdf5io
    use utility,   only : listinit,listwrite
    use gridgeneration, only : gridcube
    use geom,      only : geomcal
    use solver,    only : refcal
    !
    real(8) :: vort1,vort2,vort3,vorts,s11,s22,s33,s12,s13,s23,div,var1,miu
    real(8) :: urms,roav,uav,vav,wav,tav,pav,umax
    real(8) :: kenergy,enstrophy,dissp,miuav,divav,divvari
    real(8) :: kolmlength,ukolm,timekolm,taylorlength,retaylor,machrms,intlength,tau,macht,dudx2
    integer :: i,j,k,hand
    character(len=1) :: modeio
    !
    modeio = 'h'
    !
    if(mpirank==0)  print*, '==== Examination begins:'
    !
    !
    call readinput
    call readic
    if(mpirank==0)  print *, " ** ia:",ia,",ja:",ja,",ka:",ka
    !
    call refcal
    if(mpirank==0)  print*, ' ** refcal done!'
    !
    call mpisizedis
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parapp
    if(mpirank==0) print*, '** parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3), rho(0:im,0:jm,0:km))
    allocate(prs(0:im,0:jm,0:km), tmp(0:im,0:jm,0:km))
    allocate(dvel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3))
    call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
    call geomcal
    !
    call h5io_init(trim('datin/flowini3d.h5'),mode='read')
    call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',  var=prs(0:im,0:jm,0:km),  mode = modeio)
    call h5read(varname='t',  var=tmp(0:im,0:jm,0:km),  mode = modeio)
    call h5io_end
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Physical word calculation
    !
    call dataswap(vel)
    !
    call solvrinit
    !
    dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
    dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
    dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
    !
    urms = 0.d0
    roav = 0.d0
    uav = 0.d0
    vav = 0.d0
    wav = 0.d0
    tav = 0.d0
    pav = 0.d0
    umax = 0.d0
    !
    kenergy=0.d0
    enstrophy=0.d0
    dissp=0.d0
    miuav=0.d0
    divav=0.d0
    divvari=0.d0
    dudx2=0.d0
    !
    do i=1,im
    do j=1,jm
    do k=1,km
      var1    = vel(i,j,k,1)**2 + vel(i,j,k,2)**2 + vel(i,j,k,3)**2
      urms    = urms + var1
      roav    = roav + rho(i,j,k)
      uav     = uav  + vel(i,j,k,1)
      vav     = vav  + vel(i,j,k,2)
      wav     = wav  + vel(i,j,k,3)
      tav     = tav  + tmp(i,j,k)
      pav     = pav  + prs(i,j,k)
      umax    = max(umax,sqrt(var1))
      !
      vort1    = dvel(i,j,k,2,3)-dvel(i,j,k,3,2)
      vort2    = dvel(i,j,k,3,1)-dvel(i,j,k,1,3)
      vort3    = dvel(i,j,k,1,2)-dvel(i,j,k,2,1)
      miu      = miucal(tmp(i,j,k))/reynolds
      s11      = dvel(i,j,k,1,1)
      s22      = dvel(i,j,k,2,2)
      s33      = dvel(i,j,k,3,3)
      div      = s11+s22+s33
      s12      = 0.5d0*(dvel(i,j,k,2,1)+dvel(i,j,k,1,2))
      s13      = 0.5d0*(dvel(i,j,k,3,1)+dvel(i,j,k,1,3))
      s23      = 0.5d0*(dvel(i,j,k,3,2)+dvel(i,j,k,2,3))
      dudx2    = dudx2+s11**2+s22**2+s33**2
      !
      kenergy  = kenergy + rho(i,j,k)*var1 ! Multiplication is later
      enstrophy= enstrophy + rho(i,j,k)*(vort1**2+vort2**2+vort3**2) ! Multiplication is later
      miuav    = miuav + miu
      divav    = divav + div
      divvari  = divvari + div**2
      dissp    = dissp + 2.d0*miu*(s11**2+s22**2+s33**2+2.d0*(s12**2+s13**2+s23**2)-num1d3*div**2)
    enddo
    enddo
    enddo
    !
    var1=dble(ia*ja*ka)
    urms= sqrt(psum(urms)/var1)
    roav= psum(roav)/var1
    uav = psum(uav) /var1
    vav = psum(vav) /var1
    wav = psum(wav) /var1
    tav = psum(tav) /var1
    pav = psum(pav) /var1
    umax= pmax(umax)
    dudx2=psum(dudx2)/var1
    !
    kenergy   = psum(kenergy)  /var1*0.5d0
    enstrophy = psum(enstrophy)/var1*0.5d0
    dissp     = psum(dissp)    /var1
    miuav     = psum(miuav)    /var1
    divav     = psum(divav)    /var1
    divvari   = psum(divvari)  /var1
    !
    !
    kolmlength = sqrt(sqrt((miuav/roav)**3/dissp)) ! Kolmogorov length
    ukolm = sqrt(sqrt(dissp*miuav/roav)) ! Kolmogorov velocity
    timekolm = sqrt(miuav/roav/dissp) ! Kolmogorov time
    !
    taylorlength = sqrt(urms/dudx2) ! Taylor length
    !
    ! if(mpirank == 0) print *, "** Test usage Taylor length 1", taylorlength 
    ! taylorlength = sqrt(10.d0*kenergy*miuav/roav/dissp) ! Taylor length
    ! if(mpirank == 0) print *, "** Test usage Taylor length 2", taylorlength
    !
    retaylor = taylorlength*roav*urms/miuav/1.7320508075688773d0 ! Taylor Reynolds number
    ! intlength = sqrt(kenergy**3)/dissp ! Integral length
    !
    macht = urms/sqrt(tav)*mach ! Turbulent mach number
    ! tau = intlength / urms! Integral scale time
    !
    !
    if(mpirank == 0) then
      !
      open(unit=20,file='datin/initfield_data.dat',status='REPLACE')
      write(20,'(A10,X,E20.13E2)')'urms', urms
      write(20,'(A10,X,E20.13E2)')'roav', roav
      write(20,'(A10,X,E20.13E2)')'uav', uav
      write(20,'(A10,X,E20.13E2)')'vav', vav
      write(20,'(A10,X,E20.13E2)')'wav', wav
      write(20,'(A10,X,E20.13E2)')'tav', tav
      write(20,'(A10,X,E20.13E2)')'pav', pav
      write(20,'(A10,X,E20.13E2)')'umax', umax
      write(20,'(A10,X,E20.13E2)')'kenergy', kenergy
      write(20,'(A10,X,E20.13E2)')'ens', enstrophy
      write(20,'(A10,X,E20.13E2)')'dissp', dissp
      write(20,'(A10,X,E20.13E2)')'miuav', miuav
      write(20,'(A10,X,E20.13E2)')'Reequiv', miucal(tinf)/miuav
      write(20,'(A10,X,E20.13E2)')'divav', divav
      write(20,'(A10,X,E20.13E2)')'divvari', divvari
      write(20,'(A10,X,E20.13E2)')'eta', kolmlength
      write(20,'(A10,X,E20.13E2)')'eta/dx', kolmlength/(2*pi/ia)
      write(20,'(A10,X,E20.13E2)')'ukolm', ukolm
      write(20,'(A10,X,E20.13E2)')'tkolm', timekolm
      write(20,'(A10,X,E20.13E2)')'lambda', taylorlength
      write(20,'(A10,X,E20.13E2)')'retaylor', retaylor
      ! write(20,'(A10,X,E20.13E2)')'intlength', intlength
      ! write(20,'(A10,X,E20.13E2)')'tau', tau
      write(20,'(A10,X,E20.13E2)')'macht', macht
      CLOSE(20)
      !
      print*, " ** Print statistics"
      !
    endif
  end subroutine hitstat3d
  ! !+-------------------------------------------------------------------+
  ! !| The end of the subroutine hitstat3d.                              |
  ! !+-------------------------------------------------------------------+
  ! !
  !+-------------------------------------------------------------------+
  !| This subroutine is used calcualted the statistics of hit2d.       |
  !+-------------------------------------------------------------------+
  subroutine hitstat2d
    !
    use constdef
    use readwrite, only : readgrid, readic, readinput
    use commvar,   only : reynolds,mach,im,jm,km,ia,ja,ka,hm, tinf
    use comsolver, only : solvrinit, grad
    use commarray, only : vel,dvel,rho,prs,tmp,x
    use fludyna,   only : miucal
    use parallel,  only : mpisizedis, mpirank, parapp,dataswap, psum, pmax, parallelini
    use hdf5io
    use utility,   only : listinit,listwrite
    use gridgeneration, only : gridsquare
    use geom,      only : geomcal
    use solver,    only : refcal
    !
    real(8) :: vort3,vorts,s11,s22,s12,div,var1,miu
    real(8) :: urms,roav,uav,vav,tav,pav,umax
    real(8) :: kenergy,enstrophy,dissp,miuav,divav,divvari
    real(8) :: kolmlength,ukolm,timekolm,taylorlength,retaylor,machrms,intlength,tau,macht,dudx2
    integer :: i,j,k,hand
    character(len=1) :: modeio
    !
    modeio = 'h'
    !
    call readinput
    call readic
    if(mpirank==0)  print *, " ** ia:",ia,",ja:",ja
    !
    call refcal
    if(mpirank==0)  print*, ' ** refcal done!'
    !
    if(mpirank==0)  print*, '==== Examination begins:'
    !
    call mpisizedis
    if(mpirank==0)  print*, '** mpisizedis & parapp done!'
    !
    call parapp
    if(mpirank==0) print*, '** parapp done!'
    !
    call parallelini
    if(mpirank==0)  print*, '** parallelini done!'
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3), rho(0:im,0:jm,0:km))
    allocate(prs(0:im,0:jm,0:km), tmp(0:im,0:jm,0:km))
    allocate(dvel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3))
    call gridsquare(2.d0*pi,2.d0*pi)
    call geomcal
    !
    call h5io_init(trim('datin/flowini2d.h5'),mode='read')
    call h5read(varname='ro', var=rho(0:im,0:jm,0),  dir='k')
    call h5read(varname='u1', var=vel(0:im,0:jm,0,1),dir='k')
    call h5read(varname='u2', var=vel(0:im,0:jm,0,2),dir='k')
    call h5read(varname='p',  var=prs(0:im,0:jm,0),  dir='k')
    call h5read(varname='t',  var=tmp(0:im,0:jm,0),  dir='k')
    call h5io_end
    !
    if(mpirank==0)  print *, "Field read finish!"
    !
    ! Physical word calculation
    !
    call dataswap(vel)
    !
    call solvrinit
    !
    dvel(0:im,0:jm,0:km,1,:)=grad(vel(:,:,:,1))
    dvel(0:im,0:jm,0:km,2,:)=grad(vel(:,:,:,2))
    !
    urms = 0.d0
    roav = 0.d0
    uav = 0.d0
    vav = 0.d0
    tav = 0.d0
    pav = 0.d0
    umax = 0.d0
    !
    kenergy=0.d0
    enstrophy=0.d0
    dissp=0.d0
    miuav=0.d0
    divav=0.d0
    divvari=0.d0
    dudx2=0.d0
    !
    do i=1,im
    do j=1,jm
      var1    = vel(i,j,0,1)**2 + vel(i,j,0,2)**2
      urms    = urms + var1
      roav    = roav + rho(i,j,0)
      uav     = uav  + vel(i,j,0,1)
      vav     = vav  + vel(i,j,0,2)
      tav     = tav  + tmp(i,j,0)
      pav     = pav  + prs(i,j,0)
      umax    = max(umax,sqrt(var1))
      !
      vort3    = dvel(i,j,0,1,2)-dvel(i,j,0,2,1)
      miu      = miucal(tmp(i,j,0))/reynolds
      s11      = dvel(i,j,0,1,1)
      s22      = dvel(i,j,0,2,2)
      div      = s11+s22
      s12      = 0.5d0*(dvel(i,j,0,2,1)+dvel(i,j,0,1,2))
      dudx2    = dudx2+s11**2+s22**2
      !
      kenergy  = kenergy + rho(i,j,0)*var1 ! Multiplication is later
      enstrophy= enstrophy + rho(i,j,0)*vort3**2 ! Multiplication is later
      miuav    = miuav + miu
      divav    = divav + div
      divvari  = divvari + div**2
      dissp    = dissp + 2.d0*miu*(s11**2+s22**2+2.d0*s12**2-num1d3*div**2)
    enddo
    enddo
    !
    var1=dble(ia*ja)
    urms= sqrt(psum(urms)/var1)
    roav= psum(roav)/var1
    uav = psum(uav) /var1
    vav = psum(vav) /var1
    tav = psum(tav) /var1
    pav = psum(pav) /var1
    umax= pmax(umax)
    dudx2=psum(dudx2)/var1
    !
    kenergy   = psum(kenergy)  /var1*0.5d0
    enstrophy = psum(enstrophy)/var1*0.5d0
    dissp     = psum(dissp)    /var1
    miuav     = psum(miuav)    /var1
    divav     = psum(divav)    /var1
    divvari   = psum(divvari)  /var1
    !
    !
    kolmlength = sqrt(sqrt((miuav/roav)**3/dissp)) ! Kolmogorov length
    ukolm = sqrt(sqrt(dissp*miuav/roav)) ! Kolmogorov velocity
    timekolm = sqrt(miuav/roav/dissp) ! Kolmogorov time
    !
    taylorlength = sqrt(urms**2/dudx2) ! Taylor length
    !
    ! if(mpirank == 0) print *, "** Test usage Taylor length 1", taylorlength 
    ! taylorlength = sqrt(10.d0*kenergy*miuav/roav/dissp) ! Taylor length
    ! if(mpirank == 0) print *, "** Test usage Taylor length 2", taylorlength
    !
    retaylor = taylorlength*roav*urms/miuav/sqrt(2.d0) ! Taylor Reynolds number
    !
    macht = urms/sqrt(tav)*mach ! Turbulent mach number
    !
    !
    if(mpirank == 0) then
      !
      open(unit=20,file='datin/initfield_data.dat',status='REPLACE')
      write(20,'(A10,X,E20.13E2)')'urms', urms
      write(20,'(A10,X,E20.13E2)')'roav', roav
      write(20,'(A10,X,E20.13E2)')'uav', uav
      write(20,'(A10,X,E20.13E2)')'vav', vav
      write(20,'(A10,X,E20.13E2)')'tav', tav
      write(20,'(A10,X,E20.13E2)')'pav', pav
      write(20,'(A10,X,E20.13E2)')'umax', umax
      write(20,'(A10,X,E20.13E2)')'kenergy', kenergy
      write(20,'(A10,X,E20.13E2)')'ens', enstrophy
      write(20,'(A10,X,E20.13E2)')'dissp', dissp
      write(20,'(A10,X,E20.13E2)')'miuav', miuav
      write(20,'(A10,X,E20.13E2)')'Reequiv', miucal(tinf)/miuav
      write(20,'(A10,X,E20.13E2)')'divav', divav
      write(20,'(A10,X,E20.13E2)')'divvari', divvari
      write(20,'(A10,X,E20.13E2)')'eta', kolmlength
      write(20,'(A10,X,E20.13E2)')'eta/dx', kolmlength/(2*pi/ia)
      write(20,'(A10,X,E20.13E2)')'ukolm', ukolm
      write(20,'(A10,X,E20.13E2)')'tkolm', timekolm
      write(20,'(A10,X,E20.13E2)')'lambda', taylorlength
      write(20,'(A10,X,E20.13E2)')'retaylor', retaylor
      ! write(20,'(A10,X,E20.13E2)')'intlength', intlength
      ! write(20,'(A10,X,E20.13E2)')'tau', tau
      write(20,'(A10,X,E20.13E2)')'macht', macht
      CLOSE(20)
      !
      print*, " ** Print statistics"
      !
    endif
  end subroutine hitstat2d
  !+-------------------------------------------------------------------+
  !| The end of the subroutine hitstat2d.                              |
  !+-------------------------------------------------------------------+
  !
  function IniEnergDis(k0,wnb,method)
    !
    real(8),intent(in) :: k0,wnb
    real(8) :: Ac,var1,IniEnergDis
    integer,intent(in) :: method
    !
    !
    select case(method)
    case(1)
      ! IC4
      Ac = 2.d0/224.7699d0
      var1= -2.d0*(wnb/k0)**2
      IniEnergDis = Ac*wnb**4*exp(var1)
    case(2)
      ! pic
      if(abs(wnb - k0)< 0.5d0)then
        IniEnergDis = 1.d0
      else
        IniEnergDis = 0.d0
      endif
    case(3)
      ! Herring
      Ac = 2.d0/224.7699d0
      if (k0 /= 0.d0) then
        IniEnergDis = Ac*wnb/k0*exp(-wnb/1.d0*k0)
      else
        IniEnergDis = 0.d0
      endif
    case(4)
      ! 
      Ac = 1.d0
      if(wnb<k0)then
        IniEnergDis = Ac*((k0/wnb)**(5.d0/3.d0))
      else
        IniEnergDis = Ac*((k0/wnb)**(3.d0))
      endif
      !
    case(5)
      !
      Ac = 1.d0
      IniEnergDis = Ac*((k0/wnb)**(2.d0))
      !
    case(6)
      !
      if(wnb<k0)then
        IniEnergDis = 1.d0
      else
        IniEnergDis = 0.d0
      endif
      !
    case default
      stop "Undefined IniEnergDis method"
    end select
    !
    return
    !
  end function IniEnergDis
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of the function Ek.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  subroutine incompressuresolve3d
    !
    use commarray, only : x, vel, dvel, prs, rho
    use parallel,  only : dataswap, psum, pmin,mpirank,pmax
    use comsolver, only : grad,solvrinit
    use commvar,   only : pinf,ia,ja,ka,im,jm,km,hm,roinf,const2,reynolds
    use gridgeneration, only : gridcube
    use geom,      only : geomcal
    use fludyna,   only : thermal,miucal
    use solver,    only : refcal
    !
    integer :: i,j,k,l,m
    real(8) :: dx, pav, prsmin, equmax, miu, div, oldequ
    real(8), allocatable, dimension(:,:,:) :: P
    real(8), allocatable, dimension(:,:,:,:,:) :: sigma
    real(8), allocatable, dimension(:,:,:,:) :: dsigma11, dsigma22, dsigma33, dsigma12, dsigma13, dsigma23
    real(8), allocatable, dimension(:,:,:,:) :: ddsigma11dx1, ddsigma12dx1,ddsigma13dx1,ddsigma22dx2,ddsigma23dx2,ddsigma33dx3
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
    deallocate(prs)
    allocate(prs(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
    !
    call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
    call geomcal
    call dataswap(vel)
    allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
    allocate(P(0:im,0:jm,0:km))
    !
    if(ia == ja .and. ja==ka)then
      dx = 2.d0 * pi/ real(ia)
    else
      stop "error!"
    endif
    !
    call solvrinit
    !
    dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
    dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
    dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
    if(mpirank==0) print *, " **** dvel calculation"
    !
    allocate(sigma(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3))
    allocate(dsigma11(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
              dsigma22(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
              dsigma33(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
              dsigma12(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
              dsigma13(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
              dsigma23(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    allocate(ddsigma11dx1(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
            ddsigma12dx1(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
            ddsigma13dx1(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
            ddsigma22dx2(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
            ddsigma23dx2(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),&
            ddsigma33dx3(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
    !
    prs = pinf
    oldequ = 1000000000.d0
    equmax = 100000000.d0
    !
    ! Using iteration
    do while(equmax>10.d0 .and. (oldequ-equmax)/equmax>0.001d0)
      !
      oldequ = equmax
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        P(i,j,k) = roinf * (dvel(i,j,0,1,1) ** 2 + dvel(i,j,0,2,2) ** 2 + dvel(i,j,0,3,3) ** 2 &
                          + 2 * dvel(i,j,0,1,2) * dvel(i,j,0,2,1) + 2 * dvel(i,j,0,1,3) * dvel(i,j,0,3,1) &
                          + 2 * dvel(i,j,0,2,3) * dvel(i,j,0,3,2)) !TODO
        miu = miucal(thermal(density=rho(i,j,0),pressure=prs(i,j,0)))/Reynolds
        div = dvel(i,j,0,1,1) + dvel(i,j,0,2,2) + dvel(i,j,0,3,3)
        sigma(i,j,0,1,1) = 2 * miu * dvel(i,j,0,1,1) - 2.d0/3.d0 * miu * div
        sigma(i,j,0,2,2) = 2 * miu * dvel(i,j,0,2,2) - 2.d0/3.d0 * miu * div
        sigma(i,j,0,3,3) = 2 * miu * dvel(i,j,0,3,3) - 2.d0/3.d0 * miu * div
        sigma(i,j,0,1,2) = miu * (dvel(i,j,0,1,2) + dvel(i,j,0,2,1))
        sigma(i,j,0,1,3) = miu * (dvel(i,j,0,1,3) + dvel(i,j,0,3,1))
        sigma(i,j,0,2,3) = miu * (dvel(i,j,0,2,3) + dvel(i,j,0,3,2))
      enddo
      enddo
      enddo
      !
      call dataswap(sigma)
      !
      dsigma11(0:im,0:jm,0:km,:) = grad(sigma(:,:,:,1,1))
      dsigma22(0:im,0:jm,0:km,:) = grad(sigma(:,:,:,2,2))
      dsigma33(0:im,0:jm,0:km,:) = grad(sigma(:,:,:,3,3))
      dsigma12(0:im,0:jm,0:km,:) = grad(sigma(:,:,:,1,2))
      dsigma13(0:im,0:jm,0:km,:) = grad(sigma(:,:,:,1,3))
      dsigma23(0:im,0:jm,0:km,:) = grad(sigma(:,:,:,2,3))
      !
      call dataswap(dsigma11)
      call dataswap(dsigma22)
      call dataswap(dsigma33)
      call dataswap(dsigma12)
      call dataswap(dsigma13)
      call dataswap(dsigma23)
      !
      ddsigma11dx1(0:im,0:jm,0:km,:) = grad(dsigma11(:,:,0,1))
      ddsigma22dx2(0:im,0:jm,0:km,:) = grad(dsigma22(:,:,0,2))
      ddsigma33dx3(0:im,0:jm,0:km,:) = grad(dsigma33(:,:,0,3))
      ddsigma12dx1(0:im,0:jm,0:km,:) = grad(dsigma12(:,:,0,1))
      ddsigma13dx1(0:im,0:jm,0:km,:) = grad(dsigma13(:,:,0,1))
      ddsigma23dx2(0:im,0:jm,0:km,:) = grad(dsigma23(:,:,0,2))
      !
      do i=0,im
      do j=0,jm
      do k=0,km
        P(i,j,k) = P(i,j,k) - ddsigma11dx1(i,j,0,1) - ddsigma22dx2(i,j,0,2)  - ddsigma33dx3(i,j,0,3) &
                - 2* ddsigma12dx1(i,j,0,2)  - 2* ddsigma13dx1(i,j,0,3)  - 2* ddsigma23dx2(i,j,0,3) 
      enddo
      enddo
      enddo
      P = - P * const2 
      !
      do k=0,km
      do j=0,jm
      do i=0,im
        prs(i,j,k) =  ( - P(i,j,k) * (dx**2) &
                        + 1.d0/90.d0*prs(i-3,j,k)  - 3.d0/20.d0*prs(i-2,j,k) &
                        + 3.d0/2.d0 *prs(i-1,j,k)  + 3.d0/2.d0 *prs(i+1,j,k) &
                        - 3.d0/20.d0*prs(i+2,j,k)  + 1.d0/90.d0*prs(i+3,j,k) &
                        + 1.d0/90.d0*prs(i,j-3,k)  - 3.d0/20.d0*prs(i,j-2,k) &
                        + 3.d0/2.d0 *prs(i,j-1,k)  + 3.d0/2.d0 *prs(i,j+1,k) &
                        - 3.d0/20.d0*prs(i,j+2,k)  + 1.d0/90.d0*prs(i,j+3,k) &
                        + 1.d0/90.d0*prs(i,j,k-3)  - 3.d0/20.d0*prs(i,j,k-2) &
                        + 3.d0/2.d0 *prs(i,j,k-1)  + 3.d0/2.d0 *prs(i,j,k+1) &
                        - 3.d0/20.d0*prs(i,j,k+2)  + 1.d0/90.d0*prs(i,j,k+3))&
                        /(49.d0/6.d0)
      enddo
      enddo
      enddo
      !
      call dataswap(prs)
      !
      pav = 0.d0
      do k=1,km
      do j=1,jm
      do i=1,im
        pav =  pav + prs(i,j,k)
        prsmin = min(prsmin,prs(i,j,k))
      enddo
      enddo
      enddo
      !
      pav = psum(pav) / real(ia*ja*ka)
      prsmin = pmin(prsmin)
      !
      if(prsmin - pav + pinf>0.d0)then
        prs(:,:,:) = prs(:,:,:) - pav + pinf
      endif
      !
      equmax = 0.d0
      do k=0,km
      do j=0,jm
      do i=0,im
        equmax = max(abs(- P(i,j,k)+ &
                    + (1.d0/90.d0*prs(i-3,j,k) - 3.d0/20.d0*prs(i-2,j,k) &
                    + 3.d0/2.d0 *prs(i-1,j,k)  + 3.d0/2.d0 *prs(i+1,j,k) &
                    - 3.d0/20.d0*prs(i+2,j,k)  + 1.d0/90.d0*prs(i+3,j,k) &
                    + 1.d0/90.d0*prs(i,j-3,k)  - 3.d0/20.d0*prs(i,j-2,k) &
                    + 3.d0/2.d0 *prs(i,j-1,k)  + 3.d0/2.d0 *prs(i,j+1,k) &
                    - 3.d0/20.d0*prs(i,j+2,k)  + 1.d0/90.d0*prs(i,j+3,k) &
                    + 1.d0/90.d0*prs(i,j,k-3)  - 3.d0/20.d0*prs(i,j,k-2) &
                    + 3.d0/2.d0 *prs(i,j,k-1)  + 3.d0/2.d0 *prs(i,j,k+1) &
                    - 3.d0/20.d0*prs(i,j,k+2)  + 1.d0/90.d0*prs(i,j,k+3) &
                    - 49.d0/6.d0*prs(i,j,k))/dx/dx),equmax)
        
      enddo
      enddo
      enddo
      !
      !
      equmax = pmax(equmax)
      !
    enddo
    !
    if(mpirank==0) print *, ' **** Iteration end! equmax = ', equmax, 'oldequ=', oldequ
    !
    !
    prsmin = pmin(prsmin)
    if(prsmin< 0.d0)then
      print *, "Error! pmin<0!", prsmin
    endif
    !
    deallocate(x)
    deallocate(dvel)
    deallocate(sigma)
    deallocate(dsigma11,dsigma22,dsigma33,dsigma12,dsigma13,dsigma23)
    deallocate(ddsigma11dx1,ddsigma22dx2,ddsigma33dx3,ddsigma12dx1,ddsigma13dx1,ddsigma23dx2)
  end subroutine incompressuresolve3d
  !
  !
  subroutine incompressuresolve2d
    !
    use commarray, only : x, vel, dvel, prs, rho
    use parallel,  only : dataswap,psum,pmin,mpirank,pmax
    use comsolver, only : grad,solvrinit
    use commvar,   only : pinf,ia,ja,ka,ia,ja,ka,im,jm,km,hm,roinf,const2,reynolds
    use gridgeneration, only : gridsquare
    use geom,      only : geomcal
    use fludyna,   only : thermal,miucal
    use solver,    only : refcal
    !
    integer :: i,j,k,l,m
    real(8) :: dx, pav, prsmin, equmax, miu, div, oldequ
    real(8), allocatable, dimension(:,:) :: P
    real(8), allocatable, dimension(:,:,:,:,:) :: sigma
    real(8), allocatable, dimension(:,:,:,:) :: dsigma11, dsigma12, dsigma22, ddsigma11dx1, ddsigma12dx1, ddsigma22dx2
    !
    allocate(x(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3) )
    deallocate(prs)
    allocate(prs(-hm:im+hm,-hm:jm+hm,-hm:hm))
    !
    call gridsquare(2.d0*pi,2.d0*pi)
    call geomcal
    call dataswap(vel)
    allocate(dvel(0:im,0:jm,0:km,1:2,1:3))
    !
    allocate(P(0:im,0:jm))
    !
    if(ia == ja)then
      dx = 2.d0 * pi/ real(ia)
    else
      stop "error!"
    endif
    !
    call solvrinit
    !
    dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
    dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
    if(mpirank==0) print *, " **** dvel calculation"
    !
    allocate(sigma(-hm:im+hm,-hm:jm+hm,-hm:hm,1:2,1:2))
    allocate(dsigma11(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3),dsigma12(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3),&
            dsigma22(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3),ddsigma11dx1(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3),&
            ddsigma12dx1(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3),ddsigma22dx2(-hm:im+hm,-hm:jm+hm,-hm:hm,1:3))
    !
    prs = pinf
    oldequ = 1000000000.d0
    equmax = 100000000.d0
    !
    ! Using iteration
    do while(equmax>10.d0 .and. (oldequ-equmax)/equmax>0.001d0)
      !
      oldequ = equmax
      !
      do i=0,im
      do j=0,jm
        P(i,j) = roinf * (dvel(i,j,0,1,1) ** 2 + dvel(i,j,0,2,2) ** 2 + 2 * dvel(i,j,0,1,2) * dvel(i,j,0,2,1))
        miu = miucal(thermal(density=rho(i,j,0),pressure=prs(i,j,0)))/Reynolds
        div = dvel(i,j,0,1,1) + dvel(i,j,0,2,2)
        sigma(i,j,0,1,1) = 2 * miu * dvel(i,j,0,1,1) - 2.d0/3.d0 * miu * div
        sigma(i,j,0,1,2) = miu * (dvel(i,j,0,1,2) + dvel(i,j,0,2,1))
        sigma(i,j,0,2,2) = 2 * miu * dvel(i,j,0,2,2) - 2.d0/3.d0 * miu * div
      enddo
      enddo
      !
      call dataswap(sigma)
      !
      dsigma11(0:im,0:jm,0:0,:) = grad(sigma(:,:,:,1,1))
      dsigma12(0:im,0:jm,0:0,:) = grad(sigma(:,:,:,1,2))
      dsigma22(0:im,0:jm,0:0,:) = grad(sigma(:,:,:,2,2))
      !
      call dataswap(dsigma11)
      call dataswap(dsigma12)
      call dataswap(dsigma22)
      !
      ddsigma11dx1(0:im,0:jm,0:0,:) = grad(dsigma11(:,:,0,1))
      ddsigma12dx1(0:im,0:jm,0:0,:) = grad(dsigma12(:,:,0,1))
      ddsigma22dx2(0:im,0:jm,0:0,:) = grad(dsigma22(:,:,0,2))
      !
      !
      do i=0,im
      do j=0,jm
        P(i,j) = P(i,j) - ddsigma11dx1(i,j,0,1)  - 2* ddsigma12dx1(i,j,0,2)  - ddsigma22dx2(i,j,0,2)
      enddo
      enddo
      P = - P * const2 
      !
      !
      do j=0,jm
      do i=0,im
        prs(i,j,0) =  ( - P(i,j) * (dx**2) &
                        + 1.d0/90.d0*prs(i-3,j,0) - 3.d0/20.d0*prs(i-2,j,0) &
                        + 3.d0/2.d0 *prs(i-1,j,0) + 3.d0/2.d0 *prs(i+1,j,0) &
                        - 3.d0/20.d0*prs(i+2,j,0) + 1.d0/90.d0*prs(i+3,j,0) &
                        + 1.d0/90.d0*prs(i,j-3,0) - 3.d0/20.d0*prs(i,j-2,0) &
                        + 3.d0/2.d0 *prs(i,j-1,0) + 3.d0/2.d0 *prs(i,j+1,0) &
                        - 3.d0/20.d0*prs(i,j+2,0) + 1.d0/90.d0*prs(i,j+3,0))&
                        /(49.d0/9.d0)
      enddo
      enddo
      !
      call dataswap(prs)
      !
      pav = 0.d0
      do j=1,jm
      do i=1,im
        pav =  pav + prs(i,j,0)
        prsmin = min(prsmin,prs(i,j,0))
      enddo
      enddo
      !
      pav = psum(pav) / real(ia*ja)
      prsmin = pmin(prsmin)
      !
      if(prsmin - pav + pinf>0.d0)then
        prs(:,:,:) = prs(:,:,:) - pav + pinf
      endif
      !
      equmax = 0.d0
      do j=0,jm
      do i=0,im
        equmax = max(abs(- P(i,j)+ &
                    + (1.d0/90.d0*prs(i-3,j,0) - 3.d0/20.d0*prs(i-2,j,0) &
                    + 3.d0/2.d0 *prs(i-1,j,0)  + 3.d0/2.d0 *prs(i+1,j,0) &
                    - 3.d0/20.d0*prs(i+2,j,0)  + 1.d0/90.d0*prs(i+3,j,0) &
                    + 1.d0/90.d0*prs(i,j-3,0)  - 3.d0/20.d0*prs(i,j-2,0) &
                    + 3.d0/2.d0 *prs(i,j-1,0)  + 3.d0/2.d0 *prs(i,j+1,0) &
                    - 3.d0/20.d0*prs(i,j+2,0)  + 1.d0/90.d0*prs(i,j+3,0) &
                    - 49.d0/9.d0*prs(i,j,0))/dx/dx),equmax)
        
      enddo
      enddo
      !
      !
      equmax = pmax(equmax)
      !
    enddo
    !
    !
    if(mpirank==0) print *, ' **** Iteration end! equmax = ', equmax, 'oldequ=', oldequ
    !
    !
    prsmin = pmin(prsmin)
    if(prsmin< 0.d0)then
      print *, "Error! pmin<0!", prsmin
    endif
    !
    deallocate(x)
    deallocate(dvel)
    deallocate(sigma,dsigma11,dsigma12,dsigma22)
    deallocate(ddsigma11dx1,ddsigma12dx1,ddsigma22dx2)
  end subroutine incompressuresolve2d
  !
  subroutine scale3D(flowfile)
    !!
    use hdf5io
    use readwrite, only : readinput
    use commvar,   only : im,jm,km,ia,ja,ka
    use commarray, only : rho,vel,prs,tmp
    use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpirank, mpistop, &
                          mpi_ikgroup,mpi_kgroup,mpitag,mpiright, mpiup, mpileft, mpidown, &
                          mpifront, mpiback
    use solver,    only : refcal
    use mpi
    !
    character(len=*),intent(in) :: flowfile
    !
    real(8), allocatable, dimension(:,:,:) :: rhon,prsn,tmpn
    real(8), allocatable, dimension(:,:,:,:) :: veln
    real(8), allocatable, dimension(:,:) :: sendimjm,recvimjm,sendimkm,&
                                            recvimkm,sendjmkm,recvjmkm
    !
    integer :: i,j,k,ratio1,ratio2,ratio3,l,m,n,imn,jmn,kmn,ian,jan,kan, ierr
    integer :: status(mpi_status_size) 
    !
    character(len=1) :: modeio
    character(len=128) :: outfilename
    !
    ian = 1024
    jan = 1024
    kan = 1024
    call readinput
    !
    call mpisizedis
    if(mpirank==0) print*, '** mpisizedis done!'
    !
    call parapp
    if(mpirank==0) print*, '** parapp done!'
    !
    call parallelini
    if(mpirank==0) print*, '** parallelini done!'
    !
    call refcal
    if(mpirank==0) print*, '** refcal done!'
    !
    modeio='h'
    !
    if(mpirank==0)then
      if(ka==0)then
        print *,"2D, ia:",ia,",ja:",ja
      else
        print *, "3D grid dimensions: ia=", ia, ", ja=", ja, ", ka=", ka
      endif
    endif
    !
    !
    !allocate(x(0:im,0:jm,0:km,1:3) )
    allocate(vel(0:im,0:jm,0:km,1:3))
    allocate(rho(0:im,0:jm,0:km),prs(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km))
    !
    !call geomcal
    !
    call h5io_init(filename=flowfile,mode='read')
    !
    call h5read(varname='ro',var=rho(0:im,0:jm,0:km), mode = modeio)
    call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
    call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
    call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
    call h5read(varname='p',var=prs(0:im,0:jm,0:km), mode = modeio)
    call h5read(varname='t',var=tmp(0:im,0:jm,0:km), mode = modeio)
    !
    call h5io_end
    !
    !!!! Global size verification assumption
    if(mpirank==0)then
      print *,'ia=',ia,'ian=',ian
      print *,'ja=',ja,'jan=',jan
      print *,'ka=',ka,'kan=',kan
    endif
    !
    ! 
    !!!! Do scale
    !! 3D case
    imn = im*2
    jmn = jm*2
    kmn = km*2
    allocate(veln(0:imn,0:jmn,0:kmn,1:3))
    allocate(rhon(0:imn,0:jmn,0:kmn), &
              prsn(0:imn,0:jmn,0:kmn), &
              tmpn(0:imn,0:jmn,0:kmn))
    !
    do i=0,(im-1)
    do j=0,(jm-1)
    do k=0,(km-1)
      ! 
      veln(2*i,2*j,2*k,1)       =  vel(i,j,k,1)
      veln(2*i+1,2*j,2*k,1)     = (vel(i,j,k,1) + vel(i+1,j,k,1)    )/2
      veln(2*i+1,2*j+1,2*k,1)   = (vel(i,j,k,1) + vel(i+1,j+1,k,1)  )/2
      veln(2*i+1,2*j,2*k+1,1)   = (vel(i,j,k,1) + vel(i+1,j,k+1,1)  )/2
      veln(2*i+1,2*j+1,2*k+1,1) = (vel(i,j,k,1) + vel(i+1,j+1,k+1,1))/2
      veln(2*i,2*j+1,2*k,1)     = (vel(i,j,k,1) + vel(i,j+1,k,1)    )/2
      veln(2*i,2*j,2*k+1,1)     = (vel(i,j,k,1) + vel(i,j,k+1,1)    )/2
      veln(2*i,2*j+1,2*k+1,1)   = (vel(i,j,k,1) + vel(i,j+1,k+1,1)  )/2

      veln(2*i,2*j,2*k,2)       =  vel(i,j,k,2)
      veln(2*i+1,2*j,2*k,2)     = (vel(i,j,k,2) + vel(i+1,j,k,2)    )/2
      veln(2*i+1,2*j+1,2*k,2)   = (vel(i,j,k,2) + vel(i+1,j+1,k,2)  )/2
      veln(2*i+1,2*j,2*k+1,2)   = (vel(i,j,k,2) + vel(i+1,j,k+1,2)  )/2
      veln(2*i+1,2*j+1,2*k+1,2) = (vel(i,j,k,2) + vel(i+1,j+1,k+1,2))/2
      veln(2*i,2*j+1,2*k,2)     = (vel(i,j,k,2) + vel(i,j+1,k,2)    )/2
      veln(2*i,2*j,2*k+1,2)     = (vel(i,j,k,2) + vel(i,j,k+1,2)    )/2
      veln(2*i,2*j+1,2*k+1,2)   = (vel(i,j,k,2) + vel(i,j+1,k+1,2)  )/2

      veln(2*i,2*j,2*k,3)       =  vel(i,j,k,3)
      veln(2*i+1,2*j,2*k,3)     = (vel(i,j,k,3) + vel(i+1,j,k,3)    )/2
      veln(2*i+1,2*j+1,2*k,3)   = (vel(i,j,k,3) + vel(i+1,j+1,k,3)  )/2
      veln(2*i+1,2*j,2*k+1,3)   = (vel(i,j,k,3) + vel(i+1,j,k+1,3)  )/2
      veln(2*i+1,2*j+1,2*k+1,3) = (vel(i,j,k,3) + vel(i+1,j+1,k+1,3))/2
      veln(2*i,2*j+1,2*k,3)     = (vel(i,j,k,3) + vel(i,j+1,k,3)    )/2
      veln(2*i,2*j,2*k+1,3)     = (vel(i,j,k,3) + vel(i,j,k+1,3)    )/2
      veln(2*i,2*j+1,2*k+1,3)   = (vel(i,j,k,3) + vel(i,j+1,k+1,3)  )/2

      prsn(2*i,2*j,2*k)        = prs(i,j,k)
      prsn(2*i+1,2*j,2*k)     = (prs(i,j,k) + prs(i+1,j,k)    )/2
      prsn(2*i+1,2*j+1,2*k)   = (prs(i,j,k) + prs(i+1,j+1,k)  )/2
      prsn(2*i+1,2*j,2*k+1)   = (prs(i,j,k) + prs(i+1,j,k+1)  )/2
      prsn(2*i+1,2*j+1,2*k+1) = (prs(i,j,k) + prs(i+1,j+1,k+1))/2
      prsn(2*i,2*j+1,2*k)     = (prs(i,j,k) + prs(i,j+1,k)    )/2
      prsn(2*i,2*j,2*k+1)     = (prs(i,j,k) + prs(i,j,k+1)    )/2
      prsn(2*i,2*j+1,2*k+1)   = (prs(i,j,k) + prs(i,j+1,k+1)  )/2

      rhon(2*i,2*j,2*k)       = rho(i,j,k)
      rhon(2*i+1,2*j,2*k)     = (rho(i,j,k) + rho(i+1,j,k)     )/2
      rhon(2*i+1,2*j+1,2*k)   = (rho(i,j,k) + rho(i+1,j+1,k)   )/2
      rhon(2*i+1,2*j,2*k+1)   = (rho(i,j,k) + rho(i+1,j,k+1)   )/2
      rhon(2*i+1,2*j+1,2*k+1) = (rho(i,j,k) + rho(i+1,j+1,k+1) )/2
      rhon(2*i,2*j+1,2*k)     = (rho(i,j,k) + rho(i,j+1,k)     )/2
      rhon(2*i,2*j,2*k+1)     = (rho(i,j,k) + rho(i,j,k+1)     )/2
      rhon(2*i,2*j+1,2*k+1)   = (rho(i,j,k) + rho(i,j+1,k+1)   )/2

      tmpn(2*i,2*j,2*k)        = tmp(i,j,k)
      tmpn(2*i+1,2*j,2*k)     = (tmp(i,j,k) + tmp(i+1,j,k)     )/2
      tmpn(2*i+1,2*j+1,2*k)   = (tmp(i,j,k) + tmp(i+1,j+1,k)   )/2
      tmpn(2*i+1,2*j,2*k+1)   = (tmp(i,j,k) + tmp(i+1,j,k+1)   )/2
      tmpn(2*i+1,2*j+1,2*k+1) = (tmp(i,j,k) + tmp(i+1,j+1,k+1) )/2
      tmpn(2*i,2*j+1,2*k)     = (tmp(i,j,k) + tmp(i,j+1,k)     )/2
      tmpn(2*i,2*j,2*k+1)     = (tmp(i,j,k) + tmp(i,j,k+1)     )/2
      tmpn(2*i,2*j+1,2*k+1)   = (tmp(i,j,k) + tmp(i,j+1,k+1)   )/2
    end do
    end do
    end do
    !
    !
    !
    call mpi_sendrecv(veln(0,0,0,1) ,1,mpi_real8,mpileft ,mpitag, &
          veln(imn,0,0,1),1,mpi_real8,mpiright,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,1) ,1,mpi_real8,mpidown,mpitag, &
          veln(0,jmn,0,1),1,mpi_real8,mpiup  ,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,1) ,1,mpi_real8,mpiback ,mpitag, &
          veln(0,0,kmn,1),1,mpi_real8,mpifront,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,2) ,1,mpi_real8,mpileft ,mpitag, &
          veln(imn,0,0,2),1,mpi_real8,mpiright,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,2) ,1,mpi_real8,mpidown,mpitag, &
          veln(0,jmn,0,2),1,mpi_real8,mpiup  ,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,2) ,1,mpi_real8,mpiback ,mpitag, &
          veln(0,0,kmn,2),1,mpi_real8,mpifront,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,3) ,1,mpi_real8,mpileft ,mpitag, &
          veln(imn,0,0,3),1,mpi_real8,mpiright,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,3) ,1,mpi_real8,mpidown,mpitag, &
          veln(0,jmn,0,3),1,mpi_real8,mpiup  ,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(veln(0,0,0,3) ,1,mpi_real8,mpiback ,mpitag, &
          veln(0,0,kmn,3),1,mpi_real8,mpifront,mpitag, &
          mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    ! For prsn
    call mpi_sendrecv(prsn(0,0,0) ,1,mpi_real8,mpileft ,mpitag, &
    prsn(imn,0,0),1,mpi_real8,mpiright,mpitag, &
    mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(prsn(0,0,0) ,1,mpi_real8,mpidown,mpitag, &
      prsn(0,jmn,0),1,mpi_real8,mpiup  ,mpitag, &
      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(prsn(0,0,0) ,1,mpi_real8,mpiback ,mpitag, &
      prsn(0,0,kmn),1,mpi_real8,mpifront,mpitag, &
      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    ! For rhon
    call mpi_sendrecv(rhon(0,0,0) ,1,mpi_real8,mpileft ,mpitag, &
    rhon(imn,0,0),1,mpi_real8,mpiright,mpitag, &
    mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(rhon(0,0,0) ,1,mpi_real8,mpidown,mpitag, &
      rhon(0,jmn,0),1,mpi_real8,mpiup  ,mpitag, &
      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(rhon(0,0,0) ,1,mpi_real8,mpiback ,mpitag, &
      rhon(0,0,kmn),1,mpi_real8,mpifront,mpitag, &
      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    ! For tmpn
    call mpi_sendrecv(tmpn(0,0,0) ,1,mpi_real8,mpileft ,mpitag, &
    tmpn(imn,0,0),1,mpi_real8,mpiright,mpitag, &
    mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(tmpn(0,0,0) ,1,mpi_real8,mpidown,mpitag, &
      tmpn(0,jmn,0),1,mpi_real8,mpiup  ,mpitag, &
      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    call mpi_sendrecv(tmpn(0,0,0) ,1,mpi_real8,mpiback ,mpitag, &
      tmpn(0,0,kmn),1,mpi_real8,mpifront,mpitag, &
      mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    !
    allocate(sendimjm(0:imn,0:jmn),recvimjm(0:imn,0:jmn), &
         sendimkm(0:imn,0:kmn),recvimkm(0:imn,0:kmn), &
         sendjmkm(0:jmn,0:kmn),recvjmkm(0:jmn,0:kmn))
    !v1
    sendimjm(0:imn,0:jmn) = veln(0:imn,0:jmn,0,1)
    call mpi_sendrecv(sendimjm,(imn+1)*(jmn+1),mpi_real8,mpiback ,mpitag, &
              recvimjm,(imn+1)*(jmn+1),mpi_real8,mpifront,mpitag, &
              mpi_comm_world,status,ierr)
    veln(0:imn,0:jmn,kmn,1) = recvimjm(0:imn,0:jmn)
    mpitag=mpitag+1
    !
    sendimkm(0:imn,0:kmn) = veln(0:imn,0,0:kmn,1)
    call mpi_sendrecv(sendimkm,(imn+1)*(kmn+1),mpi_real8,mpidown,mpitag, &
              recvimkm,(imn+1)*(kmn+1),mpi_real8,mpiup  ,mpitag, &
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    veln(0:imn,jmn,0:kmn,1) = recvimkm(0:imn,0:kmn)
    !
    sendjmkm(0:jmn,0:kmn) = veln(0,0:jmn,0:kmn,1)
    call mpi_sendrecv(sendjmkm,(jmn+1)*(kmn+1),mpi_real8,mpileft ,mpitag,&
              recvjmkm,(jmn+1)*(kmn+1),mpi_real8,mpiright,mpitag,&
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    veln(imn,0:jmn,0:kmn,1) = recvjmkm(0:jmn,0:kmn)
    ! v2
    sendimjm(0:imn,0:jmn) = veln(0:imn,0:jmn,0,2)
    call mpi_sendrecv(sendimjm,(imn+1)*(jmn+1),mpi_real8,mpiback ,mpitag, &
              recvimjm,(imn+1)*(jmn+1),mpi_real8,mpifront,mpitag, &
              mpi_comm_world,status,ierr)
    veln(0:imn,0:jmn,kmn,2) = recvimjm(0:imn,0:jmn)
    mpitag=mpitag+1
    !
    sendimkm(0:imn,0:kmn) = veln(0:imn,0,0:kmn,2)
    call mpi_sendrecv(sendimkm,(imn+1)*(kmn+1),mpi_real8,mpidown,mpitag, &
              recvimkm,(imn+1)*(kmn+1),mpi_real8,mpiup  ,mpitag, &
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    veln(0:imn,jmn,0:kmn,2) = recvimkm(0:imn,0:kmn)
    !
    sendjmkm(0:jmn,0:kmn) = veln(0,0:jmn,0:kmn,2)
    call mpi_sendrecv(sendjmkm,(jmn+1)*(kmn+1),mpi_real8,mpileft ,mpitag,&
              recvjmkm,(jmn+1)*(kmn+1),mpi_real8,mpiright,mpitag,&
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    veln(imn,0:jmn,0:kmn,2) = recvjmkm(0:jmn,0:kmn)
    ! v3
    sendimjm(0:imn,0:jmn) = veln(0:imn,0:jmn,0,3)
    call mpi_sendrecv(sendimjm,(imn+1)*(jmn+1),mpi_real8,mpiback ,mpitag, &
              recvimjm,(imn+1)*(jmn+1),mpi_real8,mpifront,mpitag, &
              mpi_comm_world,status,ierr)
    veln(0:imn,0:jmn,kmn,3) = recvimjm(0:imn,0:jmn)
    mpitag=mpitag+1
    !
    sendimkm(0:imn,0:kmn) = veln(0:imn,0,0:kmn,3)
    call mpi_sendrecv(sendimkm,(imn+1)*(kmn+1),mpi_real8,mpidown,mpitag, &
              recvimkm,(imn+1)*(kmn+1),mpi_real8,mpiup  ,mpitag, &
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    veln(0:imn,jmn,0:kmn,3) = recvimkm(0:imn,0:kmn)
    !
    sendjmkm(0:jmn,0:kmn) = veln(0,0:jmn,0:kmn,3)
    call mpi_sendrecv(sendjmkm,(jmn+1)*(kmn+1),mpi_real8,mpileft ,mpitag,&
              recvjmkm,(jmn+1)*(kmn+1),mpi_real8,mpiright,mpitag,&
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    veln(imn,0:jmn,0:kmn,3) = recvjmkm(0:jmn,0:kmn)
    !
    ! prs
    sendimjm(0:imn,0:jmn) = prsn(0:imn,0:jmn,0)
    call mpi_sendrecv(sendimjm,(imn+1)*(jmn+1),mpi_real8,mpiback ,mpitag, &
              recvimjm,(imn+1)*(jmn+1),mpi_real8,mpifront,mpitag, &
              mpi_comm_world,status,ierr)
    prsn(0:imn,0:jmn,kmn) = recvimjm(0:imn,0:jmn)
    mpitag=mpitag+1
    !
    sendimkm(0:imn,0:kmn) = prsn(0:imn,0,0:kmn)
    call mpi_sendrecv(sendimkm,(imn+1)*(kmn+1),mpi_real8,mpidown,mpitag, &
              recvimkm,(imn+1)*(kmn+1),mpi_real8,mpiup  ,mpitag, &
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    prsn(0:imn,jmn,0:kmn) = recvimkm(0:imn,0:kmn)
    !
    sendjmkm(0:jmn,0:kmn) = prsn(0,0:jmn,0:kmn)
    call mpi_sendrecv(sendjmkm,(jmn+1)*(kmn+1),mpi_real8,mpileft ,mpitag,&
              recvjmkm,(jmn+1)*(kmn+1),mpi_real8,mpiright,mpitag,&
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    prsn(imn,0:jmn,0:kmn) = recvjmkm(0:jmn,0:kmn)
    ! rho
    sendimjm(0:imn,0:jmn) = rhon(0:imn,0:jmn,0)
    call mpi_sendrecv(sendimjm,(imn+1)*(jmn+1),mpi_real8,mpiback ,mpitag, &
              recvimjm,(imn+1)*(jmn+1),mpi_real8,mpifront,mpitag, &
              mpi_comm_world,status,ierr)
    rhon(0:imn,0:jmn,kmn) = recvimjm(0:imn,0:jmn)
    mpitag=mpitag+1
    !
    sendimkm(0:imn,0:kmn) = rhon(0:imn,0,0:kmn)
    call mpi_sendrecv(sendimkm,(imn+1)*(kmn+1),mpi_real8,mpidown,mpitag, &
              recvimkm,(imn+1)*(kmn+1),mpi_real8,mpiup  ,mpitag, &
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    rhon(0:imn,jmn,0:kmn) = recvimkm(0:imn,0:kmn)
    !
    sendjmkm(0:jmn,0:kmn) = rhon(0,0:jmn,0:kmn)
    call mpi_sendrecv(sendjmkm,(jmn+1)*(kmn+1),mpi_real8,mpileft ,mpitag,&
              recvjmkm,(jmn+1)*(kmn+1),mpi_real8,mpiright,mpitag,&
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    rhon(imn,0:jmn,0:kmn) = recvjmkm(0:jmn,0:kmn)
    ! tmp
    sendimjm(0:imn,0:jmn) = tmpn(0:imn,0:jmn,0)
    call mpi_sendrecv(sendimjm,(imn+1)*(jmn+1),mpi_real8,mpiback ,mpitag, &
              recvimjm,(imn+1)*(jmn+1),mpi_real8,mpifront,mpitag, &
              mpi_comm_world,status,ierr)
    tmpn(0:imn,0:jmn,kmn) = recvimjm(0:imn,0:jmn)
    mpitag=mpitag+1
    !
    sendimkm(0:imn,0:kmn) = tmpn(0:imn,0,0:kmn)
    call mpi_sendrecv(sendimkm,(imn+1)*(kmn+1),mpi_real8,mpidown,mpitag, &
              recvimkm,(imn+1)*(kmn+1),mpi_real8,mpiup  ,mpitag, &
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    tmpn(0:imn,jmn,0:kmn) = recvimkm(0:imn,0:kmn)
    !
    sendjmkm(0:jmn,0:kmn) = tmpn(0,0:jmn,0:kmn)
    call mpi_sendrecv(sendjmkm,(jmn+1)*(kmn+1),mpi_real8,mpileft ,mpitag,&
              recvjmkm,(jmn+1)*(kmn+1),mpi_real8,mpiright,mpitag,&
              mpi_comm_world,status,ierr)
    mpitag=mpitag+1
    tmpn(imn,0:jmn,0:kmn) = recvjmkm(0:jmn,0:kmn)
    !
    !
    ia = ia * 2
    ja = ja * 2
    ka = ka * 2
    !
    call parapp
    !
    deallocate(mpi_ikgroup)
    !
    call parallelini
    !
    if((imn .ne. im) .or. (jmn .ne. jm) .or. (kmn .ne. km))then
      stop "Problem emerge when redistributing the array"
    end if
    !
    !
    outfilename = flowfile//'_scaled.'//modeio//'5'
    !
    call h5io_init(trim(outfilename),mode='write')
    call h5write(varname='ro',var=rhon(0:im,0:jm,0:km), mode = modeio)
    call h5write(varname='u1',var=veln(0:im,0:jm,0:km,1), mode = modeio)
    call h5write(varname='u2',var=veln(0:im,0:jm,0:km,2), mode = modeio)
    call h5write(varname='u3',var=veln(0:im,0:jm,0:km,3), mode = modeio)
    call h5write(varname='p', var=prsn(0:im,0:jm,0:km), mode = modeio)
    call h5write(varname='t', var=tmpn(0:im,0:jm,0:km), mode = modeio)
    call h5io_end
    !
    if(mpirank==0) print*,' <<< ', outfilename, '... done.'
    !
    !
    call mpi_barrier(mpi_comm_world, ierr)
    call mpistop
    !
    deallocate(rho,vel,prs,tmp)
    deallocate(veln,rhon,prsn,tmpn)
    deallocate(recvimjm,recvimkm,recvjmkm,sendimjm,sendimkm,sendjmkm)
    !
  end subroutine scale3D
  !
  end module udf_pp_hitgen