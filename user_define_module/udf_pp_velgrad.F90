!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| spectrual calculations.                                             |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_velgrad       !
    !
    !
    use constdef
    use stlaio,  only: get_unit
    use udf_tool
    !
    implicit none
    !
    contains
    !
    subroutine ppVelgradentrance
        !
        use cmdefne
        use parallel,        only : mpirank,bcast,mpisize,lio
        !
        ! local data
        character(len=64) :: casefolder,inputfile,outputfile,viewmode, &
                             flowfieldfile, readmode
        integer :: filenumb
        !
        !
        if(mpirank == 0) then
            call readkeyboad(readmode)
            print*,' ** readmode command: ',readmode
        endif
        !
        call bcast(readmode)
        !
        if(trim(readmode)=='instant') then
        ! 
        if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        call instantvelgradient(filenumb)
        !
        elseif(trim(readmode)=='PQR') then
        ! TODO: write PQR
        if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        stop ' !! PQR not defined'
        !call PQR
        elseif(trim(readmode)=='ScaleLen') then
        !
        if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
        endif
        call bcast(filenumb)
        !
        call velgradient_scale_lengths(filenumb) 
        !
        elseif(trim(readmode)=='vortex2D') then
            !
            if(mpirank == 0) then
            call readkeyboad(inputfile) 
            read(inputfile,'(i4)') filenumb
            print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            !
            call instantvortex2D(filenumb) 
            !
        elseif(trim(readmode)=='theta3D') then
            !
            if(mpirank == 0) then
              call readkeyboad(inputfile) 
              read(inputfile,'(i4)') filenumb
              print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            !
            call instanttheta3D(filenumb) 
            !
        else
        print* ,"Readmode is not defined!", readmode
        endif
        !
    end subroutine ppVelgradentrance
    !
    subroutine instantvelgradient(thefilenumb)
        !
        !
        use readwrite, only : readinput
        use commvar,only : time,nstep,im,jm,km,hm,ia,ja,ka
        use commarray, only : x,vel,dvel
        use hdf5io
        use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank
        use comsolver, only : solvrinit,grad
        use solver,    only : refcal
        use geom,      only : geomcal
        use gridgeneration
        !
        ! arguments
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=128) :: outfilename
        character(len=1) :: modeio
        integer :: i,j
        !
        real(8), allocatable, dimension(:,:,:) :: omega,theta,psi,phi
        real(8), allocatable, dimension(:,:,:) :: alpha,beta,s12,s13,s23,omega1,omega2,omega3
        !
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
            print *,"3D, ia:",ia,",ja:",ja, ",ka:", ka
          endif
        endif
        !
        allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
        allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
        allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
        !
        if(ka==0)then
          call gridsquare(2.d0*pi,2.d0*pi)
          allocate(omega(0:im,0:jm,0:0),theta(0:im,0:jm,0:0),psi(0:im,0:jm,0:0),phi(0:im,0:jm,0:0))
        else
          call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
          allocate(theta(0:im,0:jm,0:km),alpha(0:im,0:jm,0:km),beta(0:im,0:jm,0:km), &
                    s12(0:im,0:jm,0:km),s13(0:im,0:jm,0:km),s23(0:im,0:jm,0:km),     &
                    omega1(0:im,0:jm,0:km),omega2(0:im,0:jm,0:km),omega3(0:im,0:jm,0:km))
        endif
        !
        call geomcal
        !
        if (thefilenumb .ne. 0) then
          write(stepname,'(i4.4)')thefilenumb
          infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
        else
          infilename='outdat/flowfield.'//modeio//'5'
        endif
        !
        call h5io_init(filename=infilename,mode='read')
        !
        call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
        call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
        call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        if(mpirank==0) print *, "Swap velocity"
        !
        call dataswap(vel)
        !
        call solvrinit
        !
        if(mpirank==0) print *, "Calculate gradient"
        !
        dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
        dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
        if(ka .ne. 0)then
          dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
        endif
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/velgrad'//stepname//'.'//modeio//'5'
        else
          outfilename = 'pp/velgrad.'//modeio//'5'
        endif
        
        ! !
        if(ka==0)then
          !
          theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2)
          omega(:,:,:) = dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
          psi(:,:,:)   = dvel(:,:,:,2,1) + dvel(:,:,:,1,2)
          phi(:,:,:)   = dvel(:,:,:,1,1) - dvel(:,:,:,2,2)
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          call h5write(varname='m11',var=dvel(1:im,1:jm,0:0,1,1),mode=modeio)
          call h5write(varname='m12',var=dvel(1:im,1:jm,0:0,1,2),mode=modeio)
          call h5write(varname='m21',var=dvel(1:im,1:jm,0:0,2,1),mode=modeio)
          call h5write(varname='m22',var=dvel(1:im,1:jm,0:0,2,2),mode=modeio)
          !
          call h5write(varname='theta',var=theta(1:im,1:jm,0:0),mode=modeio)
          call h5write(varname='omega',var=omega(1:im,1:jm,0:0),mode=modeio)
          call h5write(varname='psi',var=psi(1:im,1:jm,0:0),mode=modeio)
          call h5write(varname='phi',var=phi(1:im,1:jm,0:0),mode=modeio)
          !
          call h5io_end
          !
          if(mpirank==0) then
            call h5srite(varname='time',var=time,filename=trim(outfilename))
            call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
          endif
          !
          !
        else
          !
          theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2) + dvel(:,:,:,3,3)
          alpha(:,:,:) = dvel(:,:,:,1,1) - dvel(:,:,:,2,2)
          beta(:,:,:)  = dvel(:,:,:,1,1) - dvel(:,:,:,3,3)
          s12(:,:,:)   = dvel(:,:,:,1,2) + dvel(:,:,:,2,1)
          s13(:,:,:)   = dvel(:,:,:,1,3) + dvel(:,:,:,3,1)
          s23(:,:,:)   = dvel(:,:,:,2,3) + dvel(:,:,:,3,2)
          omega1(:,:,:)= dvel(:,:,:,3,2) - dvel(:,:,:,2,3)
          omega2(:,:,:)= dvel(:,:,:,1,3) - dvel(:,:,:,3,1)
          omega3(:,:,:)= dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          call h5write(varname='m11',var=dvel(1:im,1:jm,1:km,1,1),mode=modeio)
          call h5write(varname='m12',var=dvel(1:im,1:jm,1:km,1,2),mode=modeio)
          call h5write(varname='m13',var=dvel(1:im,1:jm,1:km,1,3),mode=modeio)
          call h5write(varname='m21',var=dvel(1:im,1:jm,1:km,2,1),mode=modeio)
          call h5write(varname='m22',var=dvel(1:im,1:jm,1:km,2,2),mode=modeio)
          call h5write(varname='m23',var=dvel(1:im,1:jm,1:km,2,3),mode=modeio)
          call h5write(varname='m31',var=dvel(1:im,1:jm,1:km,3,1),mode=modeio)
          call h5write(varname='m32',var=dvel(1:im,1:jm,1:km,3,2),mode=modeio)
          call h5write(varname='m33',var=dvel(1:im,1:jm,1:km,3,3),mode=modeio)
          !
          call h5write(varname='theta',var=theta(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='alpha',var=alpha(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='beta',var=beta(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='s12',var=s12(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='s23',var=s23(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='s13',var=s13(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='omega1',var=omega1(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='omega2',var=omega2(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='omega3',var=omega3(1:im,1:jm,1:km),mode=modeio)
          !
          !
          call h5io_end
          !
          if(mpirank==0) then
            call h5srite(varname='time',var=time,filename=trim(outfilename))
            call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
          endif
          !
        endif
        !
        call mpistop
        !
        deallocate(x,vel,dvel)
        if(ka==0)then
          deallocate(omega,theta,psi,phi)
        else
          deallocate(theta,alpha,beta,s12,s13,s23,omega1,omega2,omega3)
        endif
        !
      end subroutine instantvelgradient
      !
      subroutine instantvortex2D(thefilenumb)
        !
        use, intrinsic :: iso_c_binding
        use readwrite, only : readinput
        use commvar,only : time,nstep,im,jm,km,hm,ia,ja,ka,reynolds
        use commarray, only : x,vel,dvel,vorbis,dvor,rho,tmp
        use hdf5io
        use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank,psum
        use comsolver, only : solvrinit,grad
        use solver,    only : refcal
        use geom,      only : geomcal
        use gridgeneration
        use fludyna,   only : miucal
        use fftwlink
        include 'fftw3-mpi.f03'
        !
        ! arguments
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=128) :: outfilename
        character(len=1) :: modeio
        integer :: i,j,k,allkmax
        real(8) :: miu, dissa, rsamples,beta,ens1,ens2
        type(C_PTR) :: c_u1,c_u2,c_w,forward_plan,backward_plan,c_wx1,c_wx2
        complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1,u2,w,wx1,wx2
        real(8), allocatable, dimension(:,:) :: k1,k2
        complex(8) :: imag
        !
        call readinput
        !
        call refcal
        if(mpirank==0) print*, '** refcal done!'
        !
        call mpisizedis_fftw
        if(mpirank==0) print*, '** mpisizedis_fftw done!'
        !
        call parapp
        if(mpirank==0) print*, '** parapp done!'
        !
        call parallelini
        if(mpirank==0) print*, '** parallelini done!'
        !
        modeio='h'
        !
        if(mpirank==0) print *,"2D, ia:",ia,",ja:",ja
        !
        allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
        !
        allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
        allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
        allocate(dvel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3,1:3))
        allocate(vorbis(-hm:im+hm,-hm:jm+hm,-hm:km+hm))
        allocate(dvor(0:im,0:jm,0:km,1:3))
        allocate(rho(0:im,0:jm,0:km))
        allocate(tmp(0:im,0:jm,0:km))
        !
        call gridsquare(2.d0*pi,2.d0*pi)
        !
        call geomcal
        !
        if (thefilenumb .ne. 0) then
          write(stepname,'(i4.4)')thefilenumb
          infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
        else
          infilename='outdat/flowfield.'//modeio//'5'
        endif
        !
        call h5io_init(filename=infilename,mode='read')
        !
        call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
        call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
        call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
        call h5read(varname='ro', var=rho(0:im,0:jm,0:km),mode = modeio)
        call h5read(varname='t', var=tmp(0:im,0:jm,0:km),mode = modeio)
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        if(mpirank==0) print *, "Swap velocity"
        !
        call dataswap(vel)
        !
        call solvrinit
        !
        if(mpirank==0) print *, "Calculate gradient"
        !
        dvel(0:im,0:jm,0:km,1,:)=grad(vel(:,:,:,1))
        dvel(0:im,0:jm,0:km,2,:)=grad(vel(:,:,:,2))
        !
        k=0
        do j=0,jm
        do i=0,im
          vorbis(i,j,k)=dvel(i,j,k,2,1)-dvel(i,j,k,1,2)
        enddo
        enddo
        !
        call dataswap(vorbis)
        !
        !
        dvor=grad(vorbis)
        !
        !
        !! wavenumber
        allocate(k1(1:im,1:jm),k2(1:im,1:jm))
        call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
        !
        !! Imaginary number prepare
        imag = CMPLX(0.d0,1.d0,8)
        !
        !!!! Prepare initial field in Fourier space
        !! velocity
        c_u1 = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u1, u1, [imfftw,jmfftw])
        c_u2 = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_u2, u2, [imfftw,jmfftw])
        c_w = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_w, w, [imfftw,jmfftw])
        c_wx1 = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_wx1, wx1, [imfftw,jmfftw])
        c_wx2 = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_wx2, wx2, [imfftw,jmfftw])
        !
        !
        forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1,u1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
        backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, u1,u1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
        !
        do j=1,jm
        do i=1,im
          !
          u1(i,j)=CMPLX(vel(i,j,0,1),0.d0,C_INTPTR_T);
          u2(i,j)=CMPLX(vel(i,j,0,2),0.d0,C_INTPTR_T);
          !
        end do
        end do
        !
        !After this bloc, u1,u2,w are in spectral space
        call fftw_mpi_execute_dft(forward_plan,u1,u1)
        call fftw_mpi_execute_dft(forward_plan,u2,u2)
        !
        do j=1,jm
        do i=1,im
          !
          u1(i,j)=u1(i,j)/(1.d0*ia*ja)
          u2(i,j)=u2(i,j)/(1.d0*ia*ja)
          !
          w(i,j)=imag*k1(i,j)*u2(i,j)-imag*k2(i,j)*u1(i,j)
          if(sqrt(k1(i,j)**2+k2(i,j)**2) > allkmax)then
            w(i,j) = 0.d0
          endif
          !
          wx1(i,j)=imag*k1(i,j)*w(i,j)
          wx2(i,j)=imag*k2(i,j)*w(i,j)
        end do
        end do
        !
        !After this bloc,u1,u2,w,wx1,wx2 are in physical space
        call fftw_mpi_execute_dft(backward_plan,w,w)
        call fftw_mpi_execute_dft(backward_plan,wx1,wx1)
        call fftw_mpi_execute_dft(backward_plan,wx2,wx2)
        !
        beta = 0.d0
        dissa = 0.d0
        ens1 = 0.d0
        ens2 = 0.d0
        rsamples=dble(ia*ja)
        !
        k = 0
        do j=1,jm
        do i=1,im
          !
          !
          miu = miucal(tmp(i,j,k))/Reynolds
          beta = beta + miu/rho(i,j,0) * (dreal(wx1(i,j))**2 + dreal(wx2(i,j))**2)
          dissa = dissa + miu/rho(i,j,k)*(dvor(i,j,k,1)**2+dvor(i,j,k,2)**2)
          ens1 = ens1 + dreal(w(i,j))**2
          ens2 = ens2 + vorbis(i,j,k)**2
          !
        end do
        end do
        !
        beta  = psum(beta) / rsamples
        dissa= psum(dissa)/rsamples
        ens1 = psum(ens1)/rsamples
        ens2 = psum(ens2)/rsamples
        !
        if(mpirank==0)then
          print *, 'dissp: spectral, physical', dissa, beta
          print *, 'ens: spectral, physical', ens1, ens2
        endif
        !
        call mpistop
        !
        deallocate(x,vel,dvel,dvor,tmp,rho)
        !
      end subroutine instantvortex2D
      !
      subroutine velgradient_scale_lengths(thefilenumb)
        !
        !
        use singleton
        use readwrite, only : readinput
        use commvar,   only : time,nstep,im,jm,km,hm,ia,ja,ka,Reynolds
        use commarray, only : x,vel,dvel,tmp,rho
        use hdf5io
        use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank,psum
        use fludyna,   only : miucal
        use comsolver, only : solvrinit,grad
        use solver,    only : refcal
        use geom,      only : geomcal
        use gridgeneration
        !
        ! arguments
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=128) :: outfilename
        character(len=1) :: modeio
        !
        real(8), allocatable, dimension(:,:,:) :: omega,theta,psi,phi,miu
        real(8), allocatable, dimension(:,:,:) :: alpha,beta,s12,s13,s23,omega1,omega2,omega3
        real(8), allocatable, dimension(:,:,:) :: ds, eta, domega
        real(8) :: etaavg,nuavg,disspavg, rsamples
        integer ::i,j,k
        !
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
            print *,"3D, ia:",ia,",ja:",ja, ",ka:", ka
          endif
        endif
        !
        allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
        allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3),tmp(0:im,0:jm,0:km),rho(0:im,0:jm,0:km),miu(0:im,0:jm,0:km))
        allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
        !
        if(ka==0)then
          call gridsquare(2.d0*pi,2.d0*pi)
          allocate(omega(0:im,0:jm,0:0),theta(0:im,0:jm,0:0),psi(0:im,0:jm,0:0),phi(0:im,0:jm,0:0))
          allocate(ds(0:im,0:jm,0:0),eta(0:im,0:jm,0:0),domega(0:im,0:jm,0:0))
        else
          call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
          allocate(theta(0:im,0:jm,0:km),alpha(0:im,0:jm,0:km),beta(0:im,0:jm,0:km), &
                    s12(0:im,0:jm,0:km),s13(0:im,0:jm,0:km),s23(0:im,0:jm,0:km),     &
                    omega1(0:im,0:jm,0:km),omega2(0:im,0:jm,0:km),omega3(0:im,0:jm,0:km))
          allocate(ds(0:im,0:jm,0:km),eta(0:im,0:jm,0:km),domega(0:im,0:jm,0:km))
        endif
        !
        call geomcal
        !
        if (thefilenumb .ne. 0) then
          write(stepname,'(i4.4)')thefilenumb
          infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
        else
          infilename='outdat/flowfield.'//modeio//'5'
        endif
        !
        call h5io_init(filename=infilename,mode='read')
        !
        call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
        call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
        call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
        call h5read(varname='t',  var=tmp(0:im,0:jm,0:km),  mode = modeio)
        call h5read(varname='ro',var=rho(0:im,0:jm,0:km), mode = modeio)
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        if(mpirank==0) print *, "Swap velocity"
        !
        call dataswap(vel)
        !
        call solvrinit
        !
        if(mpirank==0)print *, "Calculate gradient"
        !
        dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
        dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
        if(ka .ne. 0)then
          dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
        endif
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/velgrad_scale'//stepname//'.'//modeio//'5'
        else
          outfilename = 'pp/velgrad_scale.'//modeio//'5'
        endif
        
        ! !
        !
        do k=0,km
        do j=0,jm
        do i=0,im
          miu(i,j,k) = miucal(tmp(i,j,k))/Reynolds
        enddo
        enddo
        enddo
        !
        if(ka==0)then
          !
          theta(:,:,0) = dvel(:,:,0,1,1) + dvel(:,:,0,2,2)
          omega(:,:,0) = dvel(:,:,0,2,1) - dvel(:,:,0,1,2)
          ds = sqrt(miu/rho/abs(theta))
          eta = sqrt(miu/rho/sqrt(theta**2+omega**2))
          domega = sqrt(miu/rho/abs(omega))
          !
          nuavg = 0.d0
          disspavg = 0.d0
          !
          k = 0
          do i=1,im
          do j=1,jm
            nuavg = nuavg + miu(i,j,k)/rho(i,j,k)
            disspavg = disspavg + miu(i,j,k)/rho(i,j,k) * (theta(i,j,k)**2 + omega(i,j,k)**2)
          enddo
          enddo
          !
          rsamples = dble(ia*ja)
          nuavg = psum(nuavg)/rsamples
          disspavg = psum(disspavg)/rsamples
          !
          etaavg = sqrt(sqrt(nuavg**3/disspavg))
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          call h5write(varname='ds',var=ds(1:im,1:jm,0:0),mode=modeio)
          call h5write(varname='eta',var=eta(1:im,1:jm,0:0),mode=modeio)
          call h5write(varname='domega',var=domega(1:im,1:jm,0:0),mode=modeio)
          !
          call h5io_end
          !
          if(mpirank==0) then
            call h5srite(varname='time',var=time,filename=trim(outfilename))
            call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
            call h5srite(varname='etaavg',var=etaavg,filename=trim(outfilename))
          endif
          !
          !
        else
          !
          theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2) + dvel(:,:,:,3,3)
          omega1(:,:,:)= dvel(:,:,:,3,2) - dvel(:,:,:,2,3)
          omega2(:,:,:)= dvel(:,:,:,1,3) - dvel(:,:,:,3,1)
          omega3(:,:,:)= dvel(:,:,:,2,1) - dvel(:,:,:,1,2)
          !
          ds = sqrt(miu/rho/abs(theta))
          eta = sqrt(miu/rho/sqrt(theta**2+omega1**2+omega2**2+omega3**2))
          domega = sqrt(miu/rho/sqrt(omega1**2+omega2**2+omega3**2))
          !
          nuavg = 0.d0
          disspavg = 0.d0
          !
          do i=1,im
          do j=1,jm
          do k=1,km
            nuavg = nuavg + miu(i,j,k)/rho(i,j,k)
            disspavg = disspavg + miu(i,j,k)/rho(i,j,k) * (theta(i,j,k)**2 + omega1(i,j,k)**2 + omega2(i,j,k)**2 + omega3(i,j,k)**2)
          enddo
          enddo
          enddo
          !
          rsamples = dble(ia*ja*ka)
          nuavg = psum(nuavg)/rsamples
          disspavg = psum(disspavg)/rsamples
          !
          etaavg = sqrt(sqrt(nuavg**3/disspavg))
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          call h5write(varname='ds',var=ds(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='eta',var=eta(1:im,1:jm,1:km),mode=modeio)
          call h5write(varname='domega',var=domega(1:im,1:jm,1:km),mode=modeio)
          !
          call h5io_end
          !
          if(mpirank==0) then
            call h5srite(varname='time',var=time,filename=trim(outfilename))
            call h5srite(varname='nstep',var=nstep,filename=trim(outfilename))
            call h5srite(varname='etaavg',var=etaavg,filename=trim(outfilename))
          endif
          !
        endif
        !
        call mpistop
        !
        deallocate(x,vel,dvel,tmp,rho,miu)
        if(ka==0)then
          deallocate(omega,theta,psi,phi)
        else
          deallocate(theta,alpha,beta,s12,s13,s23,omega1,omega2,omega3)
        endif
        deallocate(ds,eta,domega)
        !
      end subroutine velgradient_scale_lengths
      !
      subroutine instanttheta3D(thefilenumb)
        !
        use, intrinsic :: iso_c_binding
        use readwrite, only : readinput
        use fftwlink
        use commvar,only : time,nstep,im,jm,km,hm,ia,ja,ka
        use commarray, only : x,vel,dvel
        use hdf5io
        use parallel,  only : dataswap, mpisizedis,parapp,parallelini,mpistop,mpirank,bcast, pmax, pmin, psum, lio
        use comsolver, only : solvrinit,grad
        use solver,    only : refcal
        use geom,      only : geomcal
        use gridgeneration
        use utility,  only : listinit,listwrite
        use udf_tool, only : GenerateWave, kint
        include 'fftw3-mpi.f03'
        !
        ! arguments
        integer,intent(in) :: thefilenumb
        character(len=128) :: infilename
        character(len=4) :: stepname
        character(len=128) :: outfilename
        character(len=1) :: modeio
        integer,save :: hand_a
        integer :: i,j,k
        type(C_PTR) :: forward_plan, backward_plan, c_theta
        real(8), allocatable, dimension(:,:,:) :: theta
        complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: thetac
        real(8), allocatable, dimension(:) :: Etheta, kn
        real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
        integer :: allkmax,kOrdinal
        real(8) :: kk,dk,lambda
        real(8), allocatable, dimension(:,:) :: Sx,Sy
        real(8), allocatable, dimension(:) :: r
        integer :: n,p,in,jn,nend,pend
        real(8) :: dtheta_x, dtheta_y
        integer :: nbins,ibin
        real(8) :: thetamax, thetamin, dtheta
        real(8), allocatable, dimension(:) :: pdf, centers
        !
        call readinput
        !
        call refcal
        if(mpirank==0) print*, '** refcal done!'
        !
        call fftw_mpi_init()
        if(mpirank==0)  print *, "** fftw_mpi initialized"
        !
        call mpisizedis_fftw
        if(mpirank==0)  print*, '** mpisizedis & parapp done!'
        !
        call parapp
        if(mpirank==0) print*, '** parapp done!'
        !
        call parallelini
        if(mpirank==0) print*, '** parallelini done!'
        !
        !
        modeio='h'
        !
        if(mpirank==0) print *,"ia:",ia,",ja:",ja, ",ka:", ka
        !
        allocate(x(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3) )
        allocate(vel(-hm:im+hm,-hm:jm+hm,-hm:km+hm,1:3))
        allocate(dvel(0:im,0:jm,0:km,1:3,1:3))
        !
        call gridcube(2.d0*pi,2.d0*pi,2.d0*pi)
        allocate(theta(0:im,0:jm,0:km))
        !
        call geomcal
        !
        if (thefilenumb .ne. 0) then
          write(stepname,'(i4.4)')thefilenumb
          infilename='outdat/flowfield'//stepname//'.'//modeio//'5'
        else
          infilename='outdat/flowfield.'//modeio//'5'
        endif
        !
        call h5io_init(filename=infilename,mode='read')
        !
        call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
        call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
        call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
        call h5read(varname='time',var=time)
        call h5read(varname='nstep',var=nstep)
        !
        call h5io_end
        !
        if(mpirank==0)  print *, "** Field read finish!"
        !
        call dataswap(vel)
        !
        if(mpirank==0) print *, "** Swap velocity"
        !
        call solvrinit
        !
        if(mpirank==0) print *, "** Calculate gradient"
        !
        dvel(:,:,:,1,:)=grad(vel(:,:,:,1))
        dvel(:,:,:,2,:)=grad(vel(:,:,:,2))
        dvel(:,:,:,3,:)=grad(vel(:,:,:,3))
        !
        
        theta(:,:,:) = dvel(:,:,:,1,1) + dvel(:,:,:,2,2) + dvel(:,:,:,3,3)
        
        !
        !
        ! Calculate theta spectra
        !
        if(mpirank==0) print *, "** Spectra calculation"
        dk = 1.d0
        allkmax=ceiling(real(sqrt(2.d0)/3*min(min(ia,ja),ka))/dk)
        !
        if(mpirank==0)  print *, "**** knumber:",allkmax
        !
        !
        allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
        call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)

        c_theta = fftw_alloc_complex(alloc_local)
        call c_f_pointer(c_theta, thetac, [imfftw,jmfftw,kmfftw])
        forward_plan = fftw_mpi_plan_dft_3d(kafftw, jafftw, iafftw, thetac,thetac, &
                    MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          thetac(i,j,k)=CMPLX(theta(i,j,k),0.d0,C_INTPTR_T)
        enddo
        enddo
        enddo
        !
        call fftw_mpi_execute_dft(forward_plan,thetac,thetac)
        thetac(:,:,:)=thetac(:,:,:)/(1.d0*ia*ja*ka)
        if(mpirank==0)  print *, "**** fft done"
        !
        allocate(Etheta(0:allkmax),kn(0:allkmax))
        Etheta = 0.d0
        do k=1,km
        do j=1,jm
        do i=1,im
          kk=dsqrt(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2+1.d-15)
          kOrdinal = kint(kk,dk,3,lambda)
          if (kOrdinal<=allkmax) then
            Etheta(kOrdinal) = Etheta(kOrdinal) + thetac(i,j,k) * conjg(thetac(i,j,k))
          endif
        enddo
        enddo
        enddo
        !
        do i=0,allkmax
          Etheta(i) = psum(Etheta(i))
          kn(i) = real(i)
        enddo
        if(mpirank==0)  print *, "**** spectra calculation done"
        !
        if(mpirank == 0) then
          if (thefilenumb > 0) then
            outfilename = 'pp/Etheta'//stepname//'.dat'
          elseif(thefilenumb == 0) then
            outfilename = 'pp/Etheta.dat'
          endif
          !
          call listinit(filename=outfilename,handle=hand_a, &
                            firstline='nstep time k Etheta')
          do i=0,allkmax
            call listwrite(hand_a,kn(i),Etheta(i))
          end do
          print*,' <<< '//outfilename//'... done.'
        endif
        !
        call fftw_destroy_plan(forward_plan)
        call fftw_destroy_plan(backward_plan)
        call fftw_free(c_theta)
        deallocate(Etheta,kn,k1,k2,k3)
        ! Calculate theta structure function
        !
        if(mpirank==0) print *, "** Theta structure function"
        !
        nend = ia/4
        pend = 4
        allocate(Sx(1:pend,1:nend),Sy(1:pend,1:nend),r(1:nend))
        Sx = 0.d0
        Sy = 0.d0
        r = 0.d0
        do k=1,km
        do j=1,jm
        do i=1,im
          do n = 1,nend
            in = mod(i+n,ia)
            jn = mod(j+n,ja)
            dtheta_x = theta(in,j,k) - theta(i,j,k)
            dtheta_y = theta(i,jn,k) - theta(i,j,k)
            do p=1,pend
              Sx(p,n) = Sx(p,n) + dtheta_x ** p
              Sy(p,n) = Sy(p,n) + dtheta_y ** p
            enddo
          enddo
        enddo
        enddo
        enddo
        !
        do n = 1,nend
          do p=1,pend
            Sx(p,n) = psum(Sx(p,n))/(1.d0*ia*ja*ka)
            Sy(p,n) = psum(Sy(p,n))/(1.d0*ia*ja*ka)
          enddo
          r(n) = 2*pi/real(ia,8) * n
        enddo
        !
        if(mpirank==0)  print *, "** Structure function calculation done"
        !
        if(mpirank == 0) then
          if (thefilenumb > 0) then
            outfilename = 'pp/SFtheta'//stepname//'.dat'
          elseif(thefilenumb == 0) then
            outfilename = 'pp/SFtheta.dat'
          endif
          !
          call listinit(filename=outfilename,handle=hand_a, &
                            firstline='nstep time r Sx1 Sy1 Sx2 Sy2 Sx3 Sy3 Sx4 Sy4')
          do n=1,nend
            call listwrite(hand_a,r(n),Sx(1,n),Sy(1,n),Sx(2,n),Sy(2,n),Sx(3,n),Sy(3,n),Sx(4,n),Sy(4,n))
          end do
          print*,' <<< '//outfilename//'... done.'
        endif
        !
        deallocate(Sx,Sy,r)
        ! Calculate theta pdf
        !
        thetamax = 0.d0
        thetamin = 0.d0
        nbins = 200
        do k=1,km
        do j=1,jm
        do i=1,im
          if(thetamax<theta(i,j,k))then
            thetamax = theta(i,j,k)
          endif
          if(thetamin>theta(i,j,k))then
            thetamin = theta(i,j,k)
          endif
        enddo
        enddo
        enddo
        !
        thetamax = pmax(thetamax)
        thetamin = pmin(thetamin)
        dtheta = (thetamax - thetamin) / real(nbins-1,8)
        !
        allocate(pdf(1:nbins),centers(1:nbins))
        !
        pdf = 0.d0
        centers = 0.d0
        do k=1,km
        do j=1,jm
        do i=1,im
          ibin = int( (theta(i,j,k) - thetamin) / dtheta ) + 1
          pdf(ibin) = pdf(ibin) + 1
        enddo
        enddo
        enddo
        !
        do ibin = 1, nbins
          centers(ibin) = thetamin + (real(ibin,8) - 0.5d0) * dtheta
          pdf(ibin) = psum(pdf(ibin))/(1.d0*ia*ja*ka)/ dtheta
        end do
        !
        if(mpirank==0)  print *, "** PDF calculation done"
        !
        if(mpirank == 0) then
          if (thefilenumb > 0) then
            outfilename = 'pp/pdftheta'//stepname//'.dat'
          elseif(thefilenumb == 0) then
            outfilename = 'pp/pdftheta.dat'
          endif
          !
          call listinit(filename=outfilename,handle=hand_a, &
                            firstline='nstep time center pdf')
          do ibin=1,nbins
            call listwrite(hand_a,centers(ibin),pdf(ibin))
          end do
          print*,' <<< '//outfilename//'... done.'
        endif
        !
        deallocate(centers,pdf)
        !
        call fftw_mpi_cleanup()
        call mpistop
        !
        deallocate(x,vel,dvel,theta)
      end subroutine
end module udf_pp_velgrad