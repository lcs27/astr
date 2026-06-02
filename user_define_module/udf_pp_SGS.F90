!+---------------------------------------------------------------------+
!| This module contains subroutines for post-process concerning        |
!| physical SGS method for energy flux calculation.                    |
!+---------------------------------------------------------------------+
!| ==============                                                      |
!| CHANGE RECORD                                                       |
!| -------------                                                       |
!|  30-10-2024  | Created by C.S.Luo @ Beihang                         |
!+---------------------------------------------------------------------+
module udf_pp_SGS
    !
    !
    use constdef
    use stlaio,  only: get_unit
    use udf_tool, only: GenerateWave
    !
    implicit none
    !
    interface tensor_scale_3d
      module procedure tensor_scale_3d_real
      module procedure tensor_scale_3d_complex
    end interface
    interface tensor_multi_3d_scalar
      module procedure tensor_multi_3d_scalar_real
      module procedure tensor_multi_3d_scalar_auto_real
      module procedure tensor_multi_3d_scalar_auto_complex
      module procedure tensor_multi_3d_scalar_complex
    end interface
    contains
    !
    subroutine ppSGSentrance
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
            print*,' ** ppSGS readmode command: ',readmode
        endif
        !
        call bcast(readmode)
        !
        if(trim(readmode)=='Pi2Dint') then
        ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi2Dint"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi2Dint(filenumb)
        !
        elseif(trim(readmode)=='Pi2Dtot') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi2Dtot"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi2Dtot(filenumb)
            !
        elseif(trim(readmode)=='PiOmega2D') then
            ! 
            if(mpirank == 0) then
              print* ," ** Use SGSPiOmega2D"
              call readkeyboad(inputfile) 
              read(inputfile,'(i4)') filenumb
              print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPiOmega2D(filenumb)
            !
        elseif(trim(readmode)=='Pi2Dlocal') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi2Dlocal"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi2Dlocal(filenumb)
            !
        elseif(trim(readmode)=='E2D') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSE2D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSE2D(filenumb)
            !
        elseif(trim(readmode)=='E3D') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSE3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSE3D(filenumb)
            !
        elseif(trim(readmode)=='ET3D') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSET3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSET3D(filenumb)
            !
        elseif(trim(readmode)=='ET2D') then
          ! 
            if(mpirank == 0) then
                print* ," ** Use SGSET2D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSET2D(filenumb)
            !
        elseif(trim(readmode)=='Pi3Dint') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi3Dint"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi3Dint(filenumb)
            !
        elseif(trim(readmode)=='PiB3Dint') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPiB3Dint"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPiB3Dint(filenumb)
            !
        elseif(trim(readmode)=='Pi3Dtot') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi3Dtot"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi3Dtot(filenumb)
            !
        elseif(trim(readmode)=='Pi3Dlocal') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSPi3Dlocal"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSPi3Dlocal(filenumb)
            !
        elseif(trim(readmode)=='LES3D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSLES3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSLES3D(filenumb)
            !
        elseif(trim(readmode)=='T3D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGST3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGST3D(filenumb)
            !
        elseif(trim(readmode)=='stress2D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSstress2D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSstress2D(filenumb)
            !
        elseif(trim(readmode)=='stress3D') then
            ! 
            if(mpirank == 0) then
                print* ," ** Use SGSstress3D"
                call readkeyboad(inputfile) 
                read(inputfile,'(i4)') filenumb
                print*,' ** Filenumb: ',filenumb
            endif
            call bcast(filenumb)
            call SGSstress3D(filenumb)
            !
        else
            print* ,"Readmode is not defined!", readmode
        endif
    end subroutine ppSGSentrance
    !
    !
  subroutine SGSPiOmega2D(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka,reynolds
      use commarray, only: vel,tmp,rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use fludyna,   only : miucal
      use solver,    only : refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1,u2,w,wu1,wu2
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      integer :: allkmax
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl, beta,roav,miu,miuav,miudrho
      real(8), allocatable, dimension(:) :: Pi_omega
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: u1_filted,u2_filted,w_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wu1_filted,wu2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wx1_filted,wx2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: wx1,wx2
      !
      !
      type(C_PTR) :: c_u1,c_u2,c_w,forward_plan,backward_plan
      type(C_PTR) :: c_wx1,c_wx2,c_wu1,c_wu2
      type(C_PTR) :: c_u1_filted,c_u2_filted,c_w_filted
      type(C_PTR) :: c_wu1_filted,c_wu2_filted
      type(C_PTR) :: c_wx1_filted,c_wx2_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja, 'Reynolds:',reynolds
      allkmax=ceiling(sqrt(2.d0)/3*min(ia,ja))
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2),rho(0:im,0:jm,0:km),tmp(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),mode = modeio)
      call h5read(varname='t', var=tmp(0:im,0:jm,0:km),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
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
      c_wu1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu1, wu1, [imfftw,jmfftw])
      c_wu2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu2, wu2, [imfftw,jmfftw])
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
        !
        if(sqrt(k1(i,j)**2+k2(i,j)**2) > allkmax)then
          w(i,j) = 0.d0
        endif
        !
        wx1(i,j)=imag*k1(i,j)*w(i,j)
        wx2(i,j)=imag*k2(i,j)*w(i,j)
        !
      end do
      end do
      !
      !After this bloc,u1,u2,w,wx1,wx2 are in physical space
      call fftw_mpi_execute_dft(backward_plan,u1,u1)
      call fftw_mpi_execute_dft(backward_plan,u2,u2)
      call fftw_mpi_execute_dft(backward_plan,w,w)
      call fftw_mpi_execute_dft(backward_plan,wx1,wx1)
      call fftw_mpi_execute_dft(backward_plan,wx2,wx2)
      !
      beta = 0.d0
      roav = 0.d0
      miu = 0.d0
      miuav = 0.d0
      miudrho = 0.d0
      !!! 
      do j=1,jm
      do i=1,im
        !
        !
        wu1(i,j)=w(i,j)*u1(i,j)
        wu2(i,j)=w(i,j)*u2(i,j)
        !
        !
        miu=miucal(tmp(i,j,0))/reynolds
        beta = beta + miu/rho(i,j,0) * (dreal(wx1(i,j))**2 + dreal(wx2(i,j))**2)
        !
        roav=roav+rho(i,j,0)
        miuav=miuav+miu
        miudrho = miudrho+miu/rho(i,j,0)
        !
      end do
      end do
      !
      beta  = psum(beta) / (ia*ja)
      roav  = psum(roav) / (ia*ja)
      miuav = psum(miuav)/ (ia*ja)
      miudrho = psum(miudrho) / (ia*ja)
      !
      !After this bloc, u1,u2,w,wu1,wu2 are in spectral space
      call fftw_mpi_execute_dft(forward_plan,u1,u1)
      call fftw_mpi_execute_dft(forward_plan,u2,u2)
      call fftw_mpi_execute_dft(forward_plan,w,w)
      call fftw_mpi_execute_dft(forward_plan,wu1,wu1)
      call fftw_mpi_execute_dft(forward_plan,wu2,wu2)
      !
      do j=1,jm
      do i=1,im
        !
        u1(i,j)=u1(i,j)/(1.d0*ia*ja)
        u2(i,j)=u2(i,j)/(1.d0*ia*ja)
        w(i,j)=w(i,j)/(1.d0*ia*ja)
        wu1(i,j)=wu1(i,j)/(1.d0*ia*ja)
        wu2(i,j)=wu2(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi_omega(1:num_l))
      !
      c_u1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_u1_filted, u1_filted, [imfftw,jmfftw])
      c_u2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_u2_filted, u2_filted, [imfftw,jmfftw])
      c_w_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w_filted, w_filted, [imfftw,jmfftw])
      c_wu1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu1_filted, wu1_filted, [imfftw,jmfftw])
      c_wu2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wu2_filted, wu2_filted, [imfftw,jmfftw])
      c_wx1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wx1_filted, wx1_filted, [imfftw,jmfftw])
      c_wx2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_wx2_filted, wx2_filted, [imfftw,jmfftw])
      !
      Pi_omega = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, u1_filted is in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) !  Filtre scale :l
          !
          u1_filted(i,j) = u1(i,j) *Gl
          u2_filted(i,j) = u2(i,j) *Gl
          w_filted(i,j)  = w(i,j)  *Gl
          wu1_filted(i,j)= wu1(i,j)*Gl
          wu2_filted(i,j)= wu2(i,j)*Gl
          wx1_filted(i,j)= w(i,j)  *Gl*imag*k1(i,j)
          wx2_filted(i,j)= w(i,j)  *Gl*imag*k2(i,j)
          !
        enddo
        enddo
        !
        ! After this bloc, u1_filted is in physical space
        call fftw_mpi_execute_dft(backward_plan,u1_filted,u1_filted)
        call fftw_mpi_execute_dft(backward_plan,u2_filted,u2_filted)
        call fftw_mpi_execute_dft(backward_plan,w_filted,w_filted)
        call fftw_mpi_execute_dft(backward_plan,wu1_filted,wu1_filted)
        call fftw_mpi_execute_dft(backward_plan,wu2_filted,wu2_filted)
        call fftw_mpi_execute_dft(backward_plan,wx1_filted,wx1_filted)
        call fftw_mpi_execute_dft(backward_plan,wx2_filted,wx2_filted)
        !
        ! 
        do j=1,jm
        do i=1,im
          !
          Pi_omega(m) = Pi_omega(m) + &
              dreal(wu1_filted(i,j) - w_filted(i,j) * u1_filted(i,j))*dreal(wx1_filted(i,j)) + &
              dreal(wu2_filted(i,j) - w_filted(i,j) * u2_filted(i,j))*dreal(wx2_filted(i,j)) 
          !
        enddo
        enddo
        !
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pi_omega(m) =	 psum(Pi_omega(m)) / (ia*ja)
        !
        !
      enddo
      !
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Piomega_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Piomega.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell piomega beta miudrho miudrho2 lens')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi_omega(m), beta,miuav/roav,miudrho,((miudrho)**3/beta)**(1.d0/6.d0))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_u1)
      call fftw_free(c_u2)
      call fftw_free(c_w)
      call fftw_free(c_wx1)
      call fftw_free(c_wx2)
      call fftw_free(c_wu1)
      call fftw_free(c_wu2)
      call fftw_free(c_u1_filted)
      call fftw_free(c_u2_filted)
      call fftw_free(c_w_filted)
      call fftw_free(c_wu1_filted)
      call fftw_free(c_wu2_filted)
      call fftw_free(c_wx1_filted)
      call fftw_free(c_wx2_filted)
      call mpistop
      deallocate(k1,k2)
      deallocate(l_lim)
      deallocate(Pi_omega)
      !
    end subroutine SGSPiOmega2D
    !
    subroutine SGSPi2Dtot(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pi_tot
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1,w1w2,w2w1,w2w2
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1_filted,w1w2_filted,w2w1_filted,w2w2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
      real(8), allocatable, dimension(:,:,:,:) :: tau ! 1:2,1:2,1:im,1:jm
      !
      !
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w2w1,c_w2w2
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w2w1_filted,c_w2w2_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        w1w1(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w1w2(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w1(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w2(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        w1w1(i,j)=w1w1(i,j)/(1.d0*ia*ja)
        w1w2(i,j)=w1w2(i,j)/(1.d0*ia*ja)
        w2w1(i,j)=w2w1(i,j)/(1.d0*ia*ja)
        w2w2(i,j)=w2w2(i,j)/(1.d0*ia*ja)
        !
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      allocate(tau(1:2,1:2,1:im,1:jm))
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi_tot(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted, [imfftw,jmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted, [imfftw,jmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted, [imfftw,jmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted, [imfftw,jmfftw])
      !
      Pi_tot = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) !  Filtre scale :l
          !
          w1_filted(i,j)    = w1(i,j)   *Gl
          w2_filted(i,j)    = w2(i,j)   *Gl
          !
          w1w1_filted(i,j)  = w1w1(i,j) *Gl
          w1w2_filted(i,j)  = w1w2(i,j) *Gl
          w2w1_filted(i,j)  = w2w1(i,j) *Gl
          w2w2_filted(i,j)  = w2w2(i,j) *Gl
          !
          rho_filted(i,j)   = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do j=1,jm
        do i=1,im
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          !
          tau(1,1,i,j) = dreal(w1w1_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w1_filted(i,j))
          tau(1,2,i,j) = dreal(w1w2_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w2_filted(i,j))
          tau(2,1,i,j) = dreal(w2w1_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w1_filted(i,j))
          tau(2,2,i,j) = dreal(w2w2_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w2_filted(i,j))
          !
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        !
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          !
          A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        !
        do j=1,jm
        do i=1,im
          !
          Pi_tot(m) = Pi_tot(m) + tau(1,1,i,j) * A11_filted(i,j) + &
                                  tau(1,2,i,j) * A12_filted(i,j) + &
                                  tau(2,1,i,j) * A21_filted(i,j) + &
                                  tau(2,2,i,j) * A22_filted(i,j)
          !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pi_tot(m) =	 psum(Pi_tot(m)) / (ia*ja)
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pitot_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pitot.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pitot')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi_tot(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call mpistop
      deallocate(k1,k2,tau)
      deallocate(l_lim)
      deallocate(Pi_tot)
      !
    end subroutine SGSPi2Dtot
    !
    subroutine SGSPi2Dlocal(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pis1,Pis2,Pim2,Pim3,Pid
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1f,w2f,rhof
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
      !
      complex(8), allocatable, dimension(:,:) :: All
      complex(8), allocatable, dimension(:,:) :: S11,S12,S21,S22
      complex(8), allocatable, dimension(:,:) :: W12,W21
      !
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1f,c_w2f,c_rho
      type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      allocate(All(1:im,1:jm),&
              S11(1:im,1:jm),S12(1:im,1:jm),&
              S21(1:im,1:jm),S22(1:im,1:jm),&
              W12(1:im,1:jm),W21(1:im,1:jm))
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pis1(1:num_l),Pis2(1:num_l),Pim2(1:num_l),Pim3(1:num_l),Pid(1:num_l))
      !
      c_w1f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw])
      c_w2f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw])
      c_rho = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho, rhof,[imfftw,jmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11,[imfftw,jmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12,[imfftw,jmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21,[imfftw,jmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22,[imfftw,jmfftw])
      !
      Pis1=0.d0
      Pis2=0.d0
      Pim2=0.d0
      Pim3=0.d0
      Pid=0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1 is (rho*u1) in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1f(i,j)    = w1(i,j)    *Gl
          w2f(i,j)    = w2(i,j)    *Gl
          !
          rhof(i,j)   = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1 is (rho*u1) in physical space
        call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
        !
        ! After this bloc, w1 is u1 in physical space
        do j=1,jm
        do i=1,im
          w1f(i,j) = w1f(i,j)/rhof(i,j)
          w2f(i,j) = w2f(i,j)/rhof(i,j)
          !
          !
        enddo
        enddo
        !
        ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
        !
        do j=1,jm
        do i=1,im
          !
          w1f(i,j)  = w1f(i,j)/(1.d0*ia*ja)
          w2f(i,j)  = w2f(i,j)/(1.d0*ia*ja)
          !
          A11(i,j) = imag*w1f(i,j)*k1(i,j)
          A21(i,j) = imag*w2f(i,j)*k1(i,j)
          A12(i,j) = imag*w1f(i,j)*k2(i,j)
          A22(i,j) = imag*w2f(i,j)*k2(i,j)
          !
        end do
        end do
        !
        !
        !
        ! After this bloc, A11 is A11 in physical space
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        !
        do j=1,jm
        do i=1,im
          !
          All(i,j) = A11(i,j)+A22(i,j)
          !
          S11(i,j) = A11(i,j) - 1.d0/2.d0 * All(i,j)
          S22(i,j) = A22(i,j) - 1.d0/2.d0 * All(i,j)
          S12(i,j) = (A12(i,j) + A21(i,j))*0.5d0
          S21(i,j) = S12(i,j)
          !
          W12(i,j) = (A12(i,j)-A21(i,j))*0.5d0
          W21(i,j) = -1.d0*W12(i,j)
          !
        end do
        end do
        !
        do j=1,jm
        do i=1,im
          !
          Pis1(m) = Pis1(m) + rhof(i,j) *(S11(i,j)*S11(i,j)*S11(i,j)+ &
                                          S12(i,j)*S12(i,j)*S11(i,j)+ &
                                          S11(i,j)*S21(i,j)*S12(i,j)+ &
                                          S12(i,j)*S22(i,j)*S12(i,j)+ &
                                          S21(i,j)*S11(i,j)*S21(i,j)+ &
                                          S22(i,j)*S12(i,j)*S21(i,j)+ &
                                          S21(i,j)*S21(i,j)*S22(i,j)+ &
                                          S22(i,j)*S22(i,j)*S22(i,j)) * &
                                          l_lim(m) * l_lim(m)
          Pis2(m) = Pis2(m) + rhof(i,j) *(W12(i,j)*W21(i,j)*S11(i,j)+ &
                                            W21(i,j)*W12(i,j)*S22(i,j)) * &
                                          l_lim(m) * l_lim(m)
          Pim2(m) = Pim2(m) + rhof(i,j) *(S11(i,j)*S11(i,j)*All(i,j)+ &
                                            S12(i,j)*S12(i,j)*All(i,j)+ &
                                            S21(i,j)*S21(i,j)*All(i,j)+ &
                                            S22(i,j)*S22(i,j)*All(i,j)) * &
                                          l_lim(m) * l_lim(m)
          Pim3(m) = Pim3(m) + rhof(i,j) *(W12(i,j)*W21(i,j)*All(i,j)+ &
                                            W21(i,j)*W12(i,j)*All(i,j)) * &
                                        l_lim(m) * l_lim(m)
          Pid(m) = Pid(m) + rhof(i,j) * All(i,j) * All(i,j) * All(i,j) * &
                                      l_lim(m) * l_lim(m)
          !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pis1(m) =	 psum(Pis1(m)) / (ia*ja)
        Pis2(m) =	 - psum(Pis2(m)) / (ia*ja)
        Pim2(m) =	 psum(Pim2(m)) / (ia*ja) /2
        Pim3(m) =	 - psum(Pim3(m)) / (ia*ja) /2
        Pid(m) =	 psum(Pid(m)) / (ia*ja) /4
        !
        !
      enddo
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pilocal_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pilocal.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pis1 pis2 pim2 pim3 pid')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pis1(m), Pis2(m), Pim2(m), Pim3(m), Pid(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1f)
      call fftw_free(c_w2f)
      call fftw_free(c_rho)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call mpistop
      deallocate(k1,k2)
      deallocate(l_lim)
      deallocate(All,S11,S12,S21,S22,W12,W21)
      deallocate(Pis1,Pis2,Pim2,Pim3,Pid)
      !
    end subroutine SGSPi2Dlocal
    !
    subroutine SGSE2D(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: ES,EW,ED,E
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1f,w2f,rhof
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
      !
      complex(8), allocatable, dimension(:,:) :: All
      complex(8), allocatable, dimension(:,:) :: S11,S12,S21,S22
      complex(8), allocatable, dimension(:,:) :: W12,W21
      !
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1f,c_w2f,c_rho
      type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      allocate(All(1:im,1:jm),&
              S11(1:im,1:jm),S12(1:im,1:jm),&
              S21(1:im,1:jm),S22(1:im,1:jm),&
              W12(1:im,1:jm),W21(1:im,1:jm))
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(ES(1:num_l),EW(1:num_l),ED(1:num_l),E(1:num_l))
      !
      c_w1f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw])
      c_w2f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw])
      c_rho = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho, rhof,[imfftw,jmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11,[imfftw,jmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12,[imfftw,jmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21,[imfftw,jmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22,[imfftw,jmfftw])
      !
      ES=0.d0
      EW=0.d0
      ED=0.d0
      E=0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1 is (rho*u1) in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1f(i,j)    = w1(i,j)    *Gl
          w2f(i,j)    = w2(i,j)    *Gl
          !
          rhof(i,j)   = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1 is (rho*u1) in physical space
        call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
        !
        ! After this bloc, w1 is u1 in physical space
        do j=1,jm
        do i=1,im
          !
          w1f(i,j) = w1f(i,j)/rhof(i,j)
          w2f(i,j) = w2f(i,j)/rhof(i,j)
          !
        enddo
        enddo
        !
        ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
        !
        do j=1,jm
        do i=1,im
          !
          w1f(i,j)  = w1f(i,j)/(1.d0*ia*ja)
          w2f(i,j)  = w2f(i,j)/(1.d0*ia*ja)
          !
          A11(i,j) = imag*w1f(i,j)*k1(i,j)
          A21(i,j) = imag*w2f(i,j)*k1(i,j)
          A12(i,j) = imag*w1f(i,j)*k2(i,j)
          A22(i,j) = imag*w2f(i,j)*k2(i,j)
          !
        end do
        end do
        !
        ! After this bloc, A11 is A11 in physical space
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        !
        do j=1,jm
        do i=1,im
          !
          All(i,j) = A11(i,j)+A22(i,j)
          !
          S11(i,j) = A11(i,j) - 1.d0/2.d0 * All(i,j)
          S22(i,j) = A22(i,j) - 1.d0/2.d0 * All(i,j)
          S12(i,j) = (A12(i,j) + A21(i,j))*0.5d0
          S21(i,j) = S12(i,j)
          !
          W12(i,j) = (A12(i,j)-A21(i,j))*0.5d0
          W21(i,j) = -1.d0*W12(i,j)
          !
        end do
        end do
        !
        do j=1,jm
        do i=1,im
          !
          ES(m) = ES(m) + rhof(i,j) * (S11(i,j)*S11(i,j) + S12(i,j)*S12(i,j) + &
                                      S21(i,j)*S21(i,j) + S22(i,j)*S22(i,j)) / 2.d0
          !
          EW(m) = EW(m) + rhof(i,j) * (W12(i,j)*W12(i,j) + W21(i,j)*W21(i,j) ) / 2.d0
          !
          ED(m) = ED(m) + rhof(i,j) * (All(i,j)*All(i,j)) / 4.d0
          !
          E(m)  = E(m) +  rhof(i,j) * (A11(i,j)*A11(i,j) + A12(i,j)*A12(i,j) + &
                                       A21(i,j)*A21(i,j) + A22(i,j)*A22(i,j)) / 2.d0
          !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        ES(m) = psum(ES(m)) / (ia*ja)
        EW(m) = psum(EW(m)) / (ia*ja)
        ED(m) = psum(ED(m)) / (ia*ja)
        E(m)  = psum(E(m))  / (ia*ja)
        !
      enddo
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_E_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_E.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell ES EW ED E')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), ES(m), EW(m), ED(m), E(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1f)
      call fftw_free(c_w2f)
      call fftw_free(c_rho)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call mpistop
      deallocate(k1,k2)
      deallocate(l_lim)
      deallocate(All,S11,S12,S21,S22,W12,W21)
      deallocate(ES,EW,ED,E)
      !
    end subroutine SGSE2D
    !
    subroutine SGSET2D(thefilenumb)
      ! ! TODO : Improve need / Test need
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh,dim
      integer :: i,j,k,m,n,mmm
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_aS,hand_aW,hand_aD,hand_bS,hand_bW,hand_bD
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      real(8), allocatable, dimension(:) :: Pi1S,Pi2S,Pi3S,Pi4S,Pi5S, &
                                            Pi1W,Pi2W,Pi3W,Pi4W,Pi5W, &
                                            Pi1D,Pi2D,Pi3D,Pi4D,Pi5D
      real(8) ::  Pi1Sint,Pi2Sint,Pi3Sint,Pi4Sint,Pi5Sint,&
                  Pi1Wint,Pi2Wint,Pi3Wint,Pi4Wint,Pi5Wint,&
                  Pi1Dint,Pi2Dint,Pi3Dint,Pi4Dint,Pi5Dint
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: S1mm1_filted_l, S1mm2_filted_l, S2mm1_filted_l, S2mm2_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: W1mm1_filted_l, W1mm2_filted_l,W2mm1_filted_l, W2mm2_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: Amm11_filted_l, Amm12_filted_l,Amm21_filted_l, Amm22_filted_l
      complex(8), allocatable, dimension(:,:) :: All_filted,W12_filted,W21_filted
      complex(8), allocatable, dimension(:,:) :: S11_filted,S12_filted,S21_filted,S22_filted
      !
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: termSS_11,termSS_12,termSS_21,termSS_22,&
                                                              termWW_11,termWW_12,termWW_21,termWW_22,&
                                                              termWS_11,termWS_12,termWS_21,termWS_22,&
                                                              termDD,&
                                                              termSD_11,termSD_12,termSD_21,termSD_22
      real(8) :: vxr_D1S,vxr_D2S,vxr_D3S,vxr_D4S,vxr_D5S, &
                vxr_D1W,vxr_D2W,vxr_D3W,vxr_D4W,vxr_D5W, &
                vxr_D1D,vxr_D2D,vxr_D3D,vxr_D4D,vxr_D5D
      !
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
      type(C_PTR) :: c_termSS_11,c_termSS_12,c_termSS_21,c_termSS_22
      type(C_PTR) :: c_termWW_11,c_termWW_12,c_termWW_21,c_termWW_22
      type(C_PTR) :: c_termWS_11,c_termWS_12,c_termWS_21,c_termWS_22
      type(C_PTR) :: c_termDD
      type(C_PTR) :: c_termSD_11,c_termSD_12,c_termSD_21,c_termSD_22
      type(C_PTR) :: c_S1mm1_filted_l, c_S1mm2_filted_l,c_S2mm1_filted_l, c_S2mm2_filted_l
      type(C_PTR) :: c_W1mm1_filted_l, c_W1mm2_filted_l,c_W2mm1_filted_l, c_W2mm2_filted_l
      type(C_PTR) :: c_Amm11_filted_l, c_Amm12_filted_l,c_Amm21_filted_l, c_Amm22_filted_l
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      dim = 2
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi1S(1:num_l), Pi2S(1:num_l), Pi3S(1:num_l), Pi4S(1:num_l), Pi5S(1:num_l))
      allocate(Pi1W(1:num_l), Pi2W(1:num_l), Pi3W(1:num_l), Pi4W(1:num_l), Pi5W(1:num_l))
      allocate(Pi1D(1:num_l), Pi2D(1:num_l), Pi3D(1:num_l), Pi4D(1:num_l), Pi5D(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw])
      !
      c_S1mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S1mm1_filted_l, S1mm1_filted_l, [imfftw,jmfftw])
      c_S1mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S1mm2_filted_l, S1mm2_filted_l, [imfftw,jmfftw])
      c_S2mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S2mm1_filted_l, S2mm1_filted_l, [imfftw,jmfftw])
      c_S2mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S2mm2_filted_l, S2mm2_filted_l, [imfftw,jmfftw])
      !
      c_W1mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W1mm1_filted_l, W1mm1_filted_l, [imfftw,jmfftw])
      c_W1mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W1mm2_filted_l, W1mm2_filted_l, [imfftw,jmfftw])
      c_W2mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W2mm1_filted_l, W2mm1_filted_l, [imfftw,jmfftw])
      c_W2mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W2mm2_filted_l, W2mm2_filted_l, [imfftw,jmfftw])
      !
      c_Amm11_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm11_filted_l, Amm11_filted_l, [imfftw,jmfftw])
      c_Amm12_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm12_filted_l, Amm12_filted_l, [imfftw,jmfftw])
      c_Amm21_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm21_filted_l, Amm21_filted_l, [imfftw,jmfftw])
      c_Amm22_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm22_filted_l, Amm22_filted_l, [imfftw,jmfftw])
      !
      allocate(All_filted(1:im,1:jm),W12_filted(1:im,1:jm),W21_filted(1:im,1:jm),&
              S11_filted(1:im,1:jm),S12_filted(1:im,1:jm),S21_filted(1:im,1:jm),S22_filted(1:im,1:jm))
      !
      c_termSS_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_11, termSS_11, [imfftw,jmfftw])
      c_termSS_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_12, termSS_12, [imfftw,jmfftw])
      c_termSS_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_21, termSS_21, [imfftw,jmfftw])
      c_termSS_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_22, termSS_22, [imfftw,jmfftw])
      c_termWW_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_11, termWW_11, [imfftw,jmfftw])
      c_termWW_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_12, termWW_12, [imfftw,jmfftw])
      c_termWW_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_21, termWW_21, [imfftw,jmfftw])
      c_termWW_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_22, termWW_22, [imfftw,jmfftw])
      c_termWS_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_11, termWS_11, [imfftw,jmfftw])
      c_termWS_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_12, termWS_12, [imfftw,jmfftw])
      c_termWS_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_21, termWS_21, [imfftw,jmfftw])
      c_termWS_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_22, termWS_22, [imfftw,jmfftw])
      c_termSD_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_11, termSD_11, [imfftw,jmfftw])
      c_termSD_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_12, termSD_12, [imfftw,jmfftw])
      c_termSD_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_21, termSD_21, [imfftw,jmfftw])
      c_termSD_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_22, termSD_22, [imfftw,jmfftw])
      c_termDD = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termDD, termDD, [imfftw,jmfftw])
      !
      Pi1S =0.d0
      Pi2S =0.d0
      Pi3S =0.d0
      Pi4S =0.d0
      Pi5S =0.d0
      Pi1W =0.d0
      Pi2W =0.d0
      Pi3W =0.d0
      Pi4W =0.d0
      Pi5W =0.d0
      Pi1D =0.d0
      Pi2D =0.d0
      Pi3D =0.d0
      Pi4D =0.d0
      Pi5D =0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
            if (thefilenumb .ne. 0) then
              outfilename2 = 'pp/SGS_ETS_precise_'//stepname//'_'//mname//'.dat'
            else
              outfilename2 = 'pp/SGS_ETS_precise_'//mname//'.dat'
            endif
            call listinit(filename=outfilename2,handle=hand_bS, &
                  firstline='nstep time sqrtalpha pi1S pi2S pi3S pi4S pi5S')
            if (thefilenumb .ne. 0) then
              outfilename2 = 'pp/SGS_ETW_precise_'//stepname//'_'//mname//'.dat'
            else
              outfilename2 = 'pp/SGS_ETW_precise_'//mname//'.dat'
            endif
            call listinit(filename=outfilename2,handle=hand_bW, &
                  firstline='nstep time sqrtalpha pi1W pi2W pi3W pi4W pi5W')
            if (thefilenumb .ne. 0) then
              outfilename2 = 'pp/SGS_ETD_precise_'//stepname//'_'//mname//'.dat'
            else
              outfilename2 = 'pp/SGS_ETD_precise_'//mname//'.dat'
            endif
            call listinit(filename=outfilename2,handle=hand_bD, &
                  firstline='nstep time sqrtalpha pi1D pi2D pi3D pi4D pi5D')
        endif
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j)    = w1(i,j)    *Gl
          w2_filted(i,j)    = w2(i,j)    *Gl
          !
          rho_filted(i,j)   = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do j=1,jm
        do i=1,im
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        !
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          !
          A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
          All_filted(i,j) = (A11_filted(i,j)+A22_filted(i,j))
          !
          S11_filted(i,j) = (A11_filted(i,j)) - 1.d0/real(dim) * All_filted(i,j)
          S22_filted(i,j) = (A22_filted(i,j)) - 1.d0/real(dim) * All_filted(i,j)
          S12_filted(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
          S21_filted(i,j) = S12_filted(i,j)
          !
          W12_filted(i,j) = (A12_filted(i,j) - A21_filted(i,j))*0.5d0
          W21_filted(i,j) = -1.d0 * W12_filted(i,j)
          !
          ! 
          !SAmmB_filted(i,j) = - k1(i,j)*KB(i,j)*SA1_filted(i,j) - k2(i,j)*KB(i,j)*SA2_filted(i,j)
          !
          S1mm1_filted_l(i,j) = - k1(i,j)*k1(i,j)*S11_filted(i,j) - k2(i,j)*k1(i,j)*S12_filted(i,j)
          S1mm2_filted_l(i,j) = - k1(i,j)*k2(i,j)*S11_filted(i,j) - k2(i,j)*k2(i,j)*S12_filted(i,j)
          S2mm1_filted_l(i,j) = - k1(i,j)*k1(i,j)*S21_filted(i,j) - k2(i,j)*k1(i,j)*S22_filted(i,j)
          S2mm2_filted_l(i,j) = - k1(i,j)*k2(i,j)*S21_filted(i,j) - k2(i,j)*k2(i,j)*S22_filted(i,j)
          !
          W1mm1_filted_l(i,j) = - k2(i,j)*k1(i,j)*W12_filted(i,j)
          W1mm2_filted_l(i,j) = - k2(i,j)*k2(i,j)*W12_filted(i,j)
          W2mm1_filted_l(i,j) = - k1(i,j)*k1(i,j)*W21_filted(i,j)
          W2mm2_filted_l(i,j) = - k1(i,j)*k2(i,j)*W21_filted(i,j)
          !
          Amm11_filted_l(i,j) = - k1(i,j)*k1(i,j)*All_filted(i,j)
          Amm12_filted_l(i,j) = - k1(i,j)*k2(i,j)*All_filted(i,j)
          Amm21_filted_l(i,j) = - k2(i,j)*k1(i,j)*All_filted(i,j)
          Amm22_filted_l(i,j) = - k2(i,j)*k2(i,j)*All_filted(i,j)
          !
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,S1mm1_filted_l,S1mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S1mm2_filted_l,S1mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S2mm1_filted_l,S2mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S2mm2_filted_l,S2mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W1mm1_filted_l,W1mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W1mm2_filted_l,W1mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W2mm1_filted_l,W2mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W2mm2_filted_l,W2mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm11_filted_l,Amm11_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm12_filted_l,Amm12_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm21_filted_l,Amm21_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm22_filted_l,Amm22_filted_l)
        !
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !!!! Velocity Favre average and density average
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j)  = w1(i,j)    *Galpha
            w2_filted(i,j)  = w2(i,j)    *Galpha
            rho_filted(i,j) = rhocom(i,j)*Galpha
          enddo
          enddo
          !
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          do j=1,jm
          do i=1,im
            w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
            w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          enddo
          enddo
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
            w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
            !
            A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
            A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
            A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
            A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
            !
          end do
          end do
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
          call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
          call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
          call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
          !
          !
          do j=1,jm
          do i=1,im
            !
            All_filted(i,j) = dreal(A11_filted(i,j)+A22_filted(i,j))
            !
            S11_filted(i,j) = dreal(A11_filted(i,j)) - 1.d0/real(dim) * All_filted(i,j)
            S22_filted(i,j) = dreal(A22_filted(i,j)) - 1.d0/real(dim) * All_filted(i,j)
            S12_filted(i,j) = dreal(A12_filted(i,j) + A21_filted(i,j))*0.5d0
            S21_filted(i,j) = S12_filted(i,j)
            !
            W12_filted(i,j) = dreal(A12_filted(i,j)-A21_filted(i,j))*0.5d0
            W21_filted(i,j) = -1.d0*W12_filted(i,j)
            !
          end do
          end do
          !
          !!!! Pi terms
          !
          do j=1,jm
          do i=1,im
            rho_filted(i,j) = dreal(rho_filted(i,j))
            !termSS_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted
            termSS_11(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j) + &
                              rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)
            termSS_12(i,j) = rho_filted(i,j)*S11_filted(i,j)*S21_filted(i,j) + &
                              rho_filted(i,j)*S12_filted(i,j)*S22_filted(i,j)
            termSS_21(i,j) = rho_filted(i,j)*S21_filted(i,j)*S11_filted(i,j) + &
                              rho_filted(i,j)*S22_filted(i,j)*S12_filted(i,j)
            termSS_22(i,j) = rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j) + &
                              rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
            ! 
            ! termSD_IJ = rho_filted*All_filted*SIJ_filted
            termSD_11(i,j) = rho_filted(i,j)*All_filted(i,j)*S11_filted(i,j)
            termSD_12(i,j) = rho_filted(i,j)*All_filted(i,j)*S12_filted(i,j)
            termSD_21(i,j) = rho_filted(i,j)*All_filted(i,j)*S21_filted(i,j)
            termSD_22(i,j) = rho_filted(i,j)*All_filted(i,j)*S22_filted(i,j)
            !
            !termWW_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted
            termWW_11(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j)
            termWW_12(i,j) = 0.d0
            termWW_21(i,j) = 0.d0
            termWW_22(i,j) = rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
            !
            !
            termWS_11(i,j)= rho_filted(i,j)*S21_filted(i,j)*W12_filted(i,j) &
                            -rho_filted(i,j)*S12_filted(i,j)*W21_filted(i,j)
            termWS_21(i,j)= rho_filted(i,j)*S11_filted(i,j)*W21_filted(i,j) &
                            -rho_filted(i,j)*S22_filted(i,j)*W21_filted(i,j)
            termWS_12(i,j)=-rho_filted(i,j)*S11_filted(i,j)*W12_filted(i,j) &
                            +rho_filted(i,j)*S22_filted(i,j)*W12_filted(i,j)
            termWS_22(i,j)= rho_filted(i,j)*S12_filted(i,j)*W21_filted(i,j) &
                            -rho_filted(i,j)*S21_filted(i,j)*W12_filted(i,j)
            !
            ! termDD
            termDD(i,j) = rho_filted(i,j)*All_filted(i,j)*All_filted(i,j)
          enddo
          enddo
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,termSS_11,termSS_11)
          call fftw_mpi_execute_dft(forward_plan,termSS_12,termSS_12)
          call fftw_mpi_execute_dft(forward_plan,termSS_21,termSS_21)
          call fftw_mpi_execute_dft(forward_plan,termSS_22,termSS_22)
          !
          call fftw_mpi_execute_dft(forward_plan,termWW_11,termWW_11)
          call fftw_mpi_execute_dft(forward_plan,termWW_12,termWW_12)
          call fftw_mpi_execute_dft(forward_plan,termWW_21,termWW_21)
          call fftw_mpi_execute_dft(forward_plan,termWW_22,termWW_22)
          !
          call fftw_mpi_execute_dft(forward_plan,termWS_11,termWS_11)
          call fftw_mpi_execute_dft(forward_plan,termWS_12,termWS_12)
          call fftw_mpi_execute_dft(forward_plan,termWS_21,termWS_21)
          call fftw_mpi_execute_dft(forward_plan,termWS_22,termWS_22)
          ! !
          ! call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
          ! !
          call fftw_mpi_execute_dft(forward_plan,termSD_11,termSD_11)
          call fftw_mpi_execute_dft(forward_plan,termSD_12,termSD_12)
          call fftw_mpi_execute_dft(forward_plan,termSD_21,termSD_21)
          call fftw_mpi_execute_dft(forward_plan,termSD_22,termSD_22)
          !
          call fftw_mpi_execute_dft(forward_plan,termDD   ,termDD   )
          !
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            termSS_11(i,j) = termSS_11(i,j)*Gphi/(1.d0*ia*ja)
            termSS_12(i,j) = termSS_12(i,j)*Gphi/(1.d0*ia*ja)
            termSS_21(i,j) = termSS_21(i,j)*Gphi/(1.d0*ia*ja)
            termSS_22(i,j) = termSS_22(i,j)*Gphi/(1.d0*ia*ja)
            !
            termWW_11(i,j) = termWW_11(i,j)*Gphi/(1.d0*ia*ja)
            termWW_12(i,j) = termWW_12(i,j)*Gphi/(1.d0*ia*ja)
            termWW_21(i,j) = termWW_21(i,j)*Gphi/(1.d0*ia*ja)
            termWW_22(i,j) = termWW_22(i,j)*Gphi/(1.d0*ia*ja)
            !
            termWS_11(i,j) = termWS_11(i,j)*Gphi/(1.d0*ia*ja)
            termWS_12(i,j) = termWS_12(i,j)*Gphi/(1.d0*ia*ja)
            termWS_21(i,j) = termWS_21(i,j)*Gphi/(1.d0*ia*ja)
            termWS_22(i,j) = termWS_22(i,j)*Gphi/(1.d0*ia*ja)
            !
            termSD_11(i,j) = termSD_11(i,j)*Gphi/(1.d0*ia*ja)
            termSD_12(i,j) = termSD_12(i,j)*Gphi/(1.d0*ia*ja)
            termSD_21(i,j) = termSD_21(i,j)*Gphi/(1.d0*ia*ja)
            termSD_22(i,j) = termSD_22(i,j)*Gphi/(1.d0*ia*ja)
            !
            termDD(i,j)    = termDD(i,j)   *Gphi/(1.d0*ia*ja)
            !
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,termSS_11,termSS_11)
          call fftw_mpi_execute_dft(backward_plan,termSS_12,termSS_12)
          call fftw_mpi_execute_dft(backward_plan,termSS_21,termSS_21)
          call fftw_mpi_execute_dft(backward_plan,termSS_22,termSS_22)
          !
          call fftw_mpi_execute_dft(backward_plan,termWW_11,termWW_11)
          call fftw_mpi_execute_dft(backward_plan,termWW_12,termWW_12)
          call fftw_mpi_execute_dft(backward_plan,termWW_21,termWW_21)
          call fftw_mpi_execute_dft(backward_plan,termWW_22,termWW_22)
          !
          call fftw_mpi_execute_dft(backward_plan,termWS_11,termWS_11)
          call fftw_mpi_execute_dft(backward_plan,termWS_12,termWS_12)
          call fftw_mpi_execute_dft(backward_plan,termWS_21,termWS_21)
          call fftw_mpi_execute_dft(backward_plan,termWS_22,termWS_22)
          !
          call fftw_mpi_execute_dft(backward_plan,termSD_11,termSD_11)
          call fftw_mpi_execute_dft(backward_plan,termSD_12,termSD_12)
          call fftw_mpi_execute_dft(backward_plan,termSD_21,termSD_21)
          call fftw_mpi_execute_dft(backward_plan,termSD_22,termSD_22)
          !
          call fftw_mpi_execute_dft(backward_plan,termDD   ,termDD   )
          !
          !
          Pi1Sint = 0.d0
          Pi2Sint = 0.d0
          Pi3Sint = 0.d0
          Pi4Sint = 0.d0
          Pi5Sint = 0.d0
          Pi1Wint = 0.d0
          Pi2Wint = 0.d0
          Pi3Wint = 0.d0
          Pi4Wint = 0.d0
          Pi5Wint = 0.d0
          Pi1Dint = 0.d0
          Pi2Dint = 0.d0
          Pi3Dint = 0.d0
          Pi4Dint = 0.d0
          Pi5Dint = 0.d0
          !
          do j=1,jm
          do i=1,im
            vxr_D1S = dreal(termSS_11(i,j) * S1mm1_filted_l(i,j) + &
                          termSS_12(i,j) * S1mm2_filted_l(i,j) + &
                          termSS_21(i,j) * S2mm1_filted_l(i,j) + &
                          termSS_22(i,j) * S2mm2_filted_l(i,j))
            Pi1S(m) = Pi1S(m) + vxr_D1S * dl_alpha(m,n)
            Pi1Sint = Pi1Sint + vxr_D1S * dl_alpha(m,n)
            !
            vxr_D1W = dreal(termSS_11(i,j) * W1mm1_filted_l(i,j) + &
                          termSS_12(i,j) * W1mm2_filted_l(i,j) + &
                          termSS_21(i,j) * W2mm1_filted_l(i,j) + &
                          termSS_22(i,j) * W2mm2_filted_l(i,j))
            Pi1W(m) = Pi1W(m) + vxr_D1W * dl_alpha(m,n)
            Pi1Wint = Pi1Wint + vxr_D1W * dl_alpha(m,n)
            !
            vxr_D1D = dreal(termSS_11(i,j) * Amm11_filted_l(i,j) + &
                            termSS_12(i,j) * Amm12_filted_l(i,j) + &
                            termSS_21(i,j) * Amm21_filted_l(i,j) + &
                            termSS_22(i,j) * Amm22_filted_l(i,j))
            Pi1D(m) = Pi1D(m) + 1.d0/real(dim)*vxr_D1D * dl_alpha(m,n)
            Pi1Dint = Pi1Dint + 1.d0/real(dim)*vxr_D1D * dl_alpha(m,n)
            !
            vxr_D2S = dreal(termWW_11(i,j) * S1mm1_filted_l(i,j) + &
                    termWW_12(i,j) * S1mm2_filted_l(i,j) + &
                    termWW_21(i,j) * S2mm1_filted_l(i,j) + &
                    termWW_22(i,j) * S2mm2_filted_l(i,j))
            Pi2S(m) = Pi2S(m) - vxr_D2S * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            Pi2Sint = Pi2Sint - vxr_D2S * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            !
            vxr_D2W = dreal(termWW_11(i,j) * W1mm1_filted_l(i,j) + &
                    termWW_12(i,j) * W1mm2_filted_l(i,j) + &
                    termWW_21(i,j) * W2mm1_filted_l(i,j) + &
                    termWW_22(i,j) * W2mm2_filted_l(i,j))
            Pi2W(m) = Pi2W(m) - vxr_D2W * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            Pi2Wint = Pi2Wint - vxr_D2W * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            !
            vxr_D2D = dreal(termWW_11(i,j) * Amm11_filted_l(i,j) + &
                    termWW_12(i,j) * Amm12_filted_l(i,j) + &
                    termWW_21(i,j) * Amm21_filted_l(i,j) + &
                    termWW_22(i,j) * Amm22_filted_l(i,j))
            Pi2D(m) = Pi2D(m) - 1.d0/real(dim)*vxr_D2D * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            Pi2Dint = Pi2Dint - 1.d0/real(dim)*vxr_D2D * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            !
            vxr_D3S = dreal(termWS_11(i,j) * S1mm1_filted_l(i,j) + &
                          termWS_12(i,j) * S1mm2_filted_l(i,j) + &
                          termWS_21(i,j) * S2mm1_filted_l(i,j) + &
                          termWS_22(i,j) * S2mm2_filted_l(i,j))
            Pi3S(m) = Pi3S(m) + vxr_D3S * dl_alpha(m,n)
            Pi3Sint = Pi3Sint + vxr_D3S * dl_alpha(m,n)
            !
            vxr_D3W = dreal(termWS_11(i,j) * W1mm1_filted_l(i,j) + &
                          termWS_12(i,j) * W1mm2_filted_l(i,j) + &
                          termWS_21(i,j) * W2mm1_filted_l(i,j) + &
                          termWS_22(i,j) * W2mm2_filted_l(i,j))
            Pi3W(m) = Pi3W(m) + vxr_D3W * dl_alpha(m,n)
            Pi3Wint = Pi3Wint + vxr_D3W * dl_alpha(m,n)
            !
            !
            vxr_D3D = dreal(termWS_11(i,j) * Amm11_filted_l(i,j) + &
                    termWS_12(i,j) * Amm12_filted_l(i,j) + &
                    termWS_21(i,j) * Amm21_filted_l(i,j) + &
                    termWS_22(i,j) * Amm22_filted_l(i,j))
            Pi3D(m) = Pi3D(m) + 1.d0/real(dim)*vxr_D3D * dl_alpha(m,n)
            Pi3Dint = Pi3Dint + 1.d0/real(dim)*vxr_D3D * dl_alpha(m,n)
            !
            vxr_D4S = dreal(termSD_11(i,j) * S1mm1_filted_l(i,j) + &
                          termSD_12(i,j) * S1mm2_filted_l(i,j) + &
                          termSD_21(i,j) * S2mm1_filted_l(i,j) + &
                          termSD_22(i,j) * S2mm2_filted_l(i,j))
            Pi4S(m) = Pi4S(m) + 2.d0/real(dim) * vxr_D4S * dl_alpha(m,n)
            Pi4Sint = Pi4Sint + 2.d0/real(dim) * vxr_D4S * dl_alpha(m,n)
            !
            vxr_D4W = dreal(termSD_11(i,j) * W1mm1_filted_l(i,j) + &
                          termSD_12(i,j) * W1mm2_filted_l(i,j) + &
                          termSD_21(i,j) * W2mm1_filted_l(i,j) + &
                          termSD_22(i,j) * W2mm2_filted_l(i,j))
            Pi4W(m) = Pi4W(m) + 2.d0/real(dim) * vxr_D4W * dl_alpha(m,n)
            Pi4Wint = Pi4Wint + 2.d0/real(dim) * vxr_D4W * dl_alpha(m,n)
            !
            !
            vxr_D4D = dreal(termSD_11(i,j) * Amm11_filted_l(i,j) + &
                    termSD_12(i,j) * Amm12_filted_l(i,j) + &
                    termSD_21(i,j) * Amm21_filted_l(i,j) + &
                    termSD_22(i,j) * Amm22_filted_l(i,j))
            Pi4D(m) = Pi4D(m) + 2.d0/real(dim)/real(dim)* vxr_D4D * dl_alpha(m,n)
            Pi4Dint = Pi4Dint + 2.d0/real(dim)/real(dim)* vxr_D4D * dl_alpha(m,n)
            !
            !
            vxr_D5S = dreal(termDD(i,j) * (S1mm1_filted_l(i,j) + S2mm2_filted_l(i,j)))
            Pi5S(m) = Pi5S(m) + 1.d0/real(dim)/real(dim)*vxr_D5S * dl_alpha(m,n)
            Pi5Sint = Pi5Sint + 1.d0/real(dim)/real(dim)*vxr_D5S * dl_alpha(m,n)
            !
            vxr_D5W = dreal(termDD(i,j) * (W1mm1_filted_l(i,j) + W2mm2_filted_l(i,j)))
            Pi5W(m) = Pi5W(m) + 1.d0/real(dim)/real(dim)*vxr_D5W * dl_alpha(m,n)
            Pi5Wint = Pi5Wint + 1.d0/real(dim)/real(dim)*vxr_D5W * dl_alpha(m,n)
            !
            !
            vxr_D5D = dreal(termDD(i,j) * (Amm11_filted_l(i,j) + Amm22_filted_l(i,j)))
            Pi5D(m) = Pi5D(m) + 1.d0/real(dim)/real(dim)/real(dim)*vxr_D5D * dl_alpha(m,n)
            Pi5Dint = Pi5Dint + 1.d0/real(dim)/real(dim)/real(dim)*vxr_D5D * dl_alpha(m,n)
            !
          enddo
          enddo
          !
          Pi1Sint = psum(Pi1Sint) / (ia*ja)
          Pi2Sint = psum(Pi2Sint) / (ia*ja)
          Pi3Sint = psum(Pi3Sint) / (ia*ja)
          Pi4Sint = psum(Pi4Sint) / (ia*ja)
          Pi5Sint = psum(Pi5Sint) / (ia*ja)
          Pi1Wint = psum(Pi1Wint) / (ia*ja)
          Pi2Wint = psum(Pi2Wint) / (ia*ja)
          Pi3Wint = psum(Pi3Wint) / (ia*ja)
          Pi4Wint = psum(Pi4Wint) / (ia*ja)
          Pi5Wint = psum(Pi5Wint) / (ia*ja)
          Pi1Dint = psum(Pi1Dint) / (ia*ja)
          Pi2Dint = psum(Pi2Dint) / (ia*ja)
          Pi3Dint = psum(Pi3Dint) / (ia*ja)
          Pi4Dint = psum(Pi4Dint) / (ia*ja)
          Pi5Dint = psum(Pi5Dint) / (ia*ja)
          !
          if(mpirank==0) then
            call listwrite(hand_bS,l_sqrtalpha(m,n),Pi1Sint, Pi2Sint,Pi3Sint, &
              Pi4Sint, Pi5Sint)
            call listwrite(hand_bD,l_sqrtalpha(m,n),Pi1Dint, Pi2Dint,Pi3Dint, &
              Pi4Dint, Pi5Dint)
            call listwrite(hand_bW,l_sqrtalpha(m,n),Pi1Wint, Pi2Wint,Pi3Wint, &
              Pi4Wint, Pi5Wint)
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        !
        Pi1S(m) =  psum(Pi1S(m)) / (ia*ja)
        Pi2S(m) =  psum(Pi2S(m)) / (ia*ja)
        Pi3S(m) =  psum(Pi3S(m)) / (ia*ja)
        Pi4S(m) =  psum(Pi4S(m)) / (ia*ja)
        Pi5S(m) =  psum(Pi5S(m)) / (ia*ja)
        Pi1W(m) =  psum(Pi1W(m)) / (ia*ja)
        Pi2W(m) =  psum(Pi2W(m)) / (ia*ja)
        Pi3W(m) =  psum(Pi3W(m)) / (ia*ja)
        Pi4W(m) =  psum(Pi4W(m)) / (ia*ja)
        Pi5W(m) =  psum(Pi5W(m)) / (ia*ja)
        Pi1D(m) =  psum(Pi1D(m)) / (ia*ja)
        Pi2D(m) =  psum(Pi2D(m)) / (ia*ja)
        Pi3D(m) =  psum(Pi3D(m)) / (ia*ja)
        Pi4D(m) =  psum(Pi4D(m)) / (ia*ja)
        Pi5D(m) =  psum(Pi5D(m)) / (ia*ja)
        !
        !
        !
        if(mpirank==0) then
            call listwrite(hand_bS, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
            call listwrite(hand_bS,Pi1S(m)+Pi2S(m)+Pi3S(m)+Pi4S(m)+Pi5S(m), & 
            Pi1S(m), Pi2S(m),Pi3S(m),Pi4S(m),Pi5S(m))
            call listwrite(hand_bD, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
            call listwrite(hand_bD,Pi1D(m)+Pi2D(m)+Pi3D(m)+Pi4D(m)+Pi5D(m), & 
            Pi1D(m), Pi2D(m),Pi3D(m),Pi4D(m),Pi5D(m))
            call listwrite(hand_bW, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
            call listwrite(hand_bW,Pi1W(m)+Pi2W(m)+Pi3W(m)+Pi4W(m)+Pi5W(m), & 
            Pi1W(m), Pi2W(m),Pi3W(m),Pi4W(m),Pi5W(m))
          !
          close(unit=hand_bS)
          close(unit=hand_bD)
          close(unit=hand_bW)
          !
          print *, '>>>>', outfilename2
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_ETS_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_ETS.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_aS, &
                      firstline='nstep time ell pi1S pi2S pi3S pi4S pi5S')
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_ETW_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_ETW.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_aW, &
                      firstline='nstep time ell pi1W pi2W pi3W pi4W pi5W')
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_ETD_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_ETD.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_aD, &
                      firstline='nstep time ell pi1D pi2D pi3D pi4D pi5D')
        !
        do m=1,num_l
          call listwrite(hand_aS,l_lim(m),Pi1S(m), Pi2S(m),&
            Pi3S(m), Pi4S(m), Pi5S(m))
          call listwrite(hand_aW,l_lim(m),Pi1W(m), Pi2W(m),&
            Pi3W(m), Pi4W(m), Pi5W(m))
          call listwrite(hand_aD,l_lim(m),Pi1D(m), Pi2D(m),&
            Pi3D(m), Pi4D(m), Pi5D(m))
        enddo
        !
        close(unit=hand_aS)
        close(unit=hand_aW)
        close(unit=hand_aD)
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_termSS_11)
      call fftw_free(c_termSS_12)
      call fftw_free(c_termSS_21)
      call fftw_free(c_termSS_22)
      call fftw_free(c_termSD_11)
      call fftw_free(c_termSD_12)
      call fftw_free(c_termSD_21)
      call fftw_free(c_termSD_22)
      call fftw_free(c_termWW_11)
      call fftw_free(c_termWW_12)
      call fftw_free(c_termWW_21)
      call fftw_free(c_termWW_22)
      call fftw_free(c_termWS_11)
      call fftw_free(c_termWS_12)
      call fftw_free(c_termWS_21)
      call fftw_free(c_termWS_22)
      call fftw_free(c_termDD)
      call fftw_free(c_S1mm1_filted_l)
      call fftw_free(c_S1mm2_filted_l)
      call fftw_free(c_S2mm1_filted_l)
      call fftw_free(c_S2mm2_filted_l)
      call fftw_free(c_W1mm1_filted_l)
      call fftw_free(c_W1mm2_filted_l)
      call fftw_free(c_W2mm1_filted_l)
      call fftw_free(c_W2mm2_filted_l)
      call mpistop
      deallocate(All_filted,S11_filted,S12_filted,S21_filted,S22_filted)
      deallocate(W12_filted,W21_filted)
      deallocate(k1,k2,k3)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
      deallocate(Pi1S,Pi2S,Pi3S,Pi4S,Pi5S,Pi1W,Pi2W,Pi3W,Pi4W,Pi5W,Pi1D,Pi2D,Pi3D,Pi4D,Pi5D)
      !
    end subroutine SGSET2D
    !
    subroutine SGSE3D(thefilenumb)
      ! ! TODO : Improve need / Test need
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: ES,EW,ED,E
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1f,w2f,w3f,rhof
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
      !
      real(8), allocatable, dimension(:,:,:) :: All
      real(8), allocatable, dimension(:,:,:) :: S11,S12,S13
      real(8), allocatable, dimension(:,:,:) :: S21,S22,S23
      real(8), allocatable, dimension(:,:,:) :: S31,S32,S33
      real(8), allocatable, dimension(:,:,:) :: W12,W13,W21,W23,W31,W32
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1f,c_w2f,c_w3f,c_rho
      type(C_PTR) :: c_A11,c_A12,c_A13
      type(C_PTR) :: c_A21,c_A22,c_A23
      type(C_PTR) :: c_A31,c_A32,c_A33
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      allocate(All(1:im,1:jm,1:km),&
          S11(1:im,1:jm,1:km),S12(1:im,1:jm,1:km),S13(1:im,1:jm,1:km),&
          S21(1:im,1:jm,1:km),S22(1:im,1:jm,1:km),S23(1:im,1:jm,1:km),&
          S31(1:im,1:jm,1:km),S32(1:im,1:jm,1:km),S33(1:im,1:jm,1:km),&
          W12(1:im,1:jm,1:km),W13(1:im,1:jm,1:km),W21(1:im,1:jm,1:km),&
          W23(1:im,1:jm,1:km),W31(1:im,1:jm,1:km),W32(1:im,1:jm,1:km))
      !
      do k=1,km
      do j=1,jm
      do i=1,im
      w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T)
      w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T)
      w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T)
      rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T)
      end do
      end do
      end do
      !
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
      w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
      w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
      w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
      rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
      end do
      end do
      end do
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      if(mpirank==0)  print *, "Integrate point allocated"
      call mpi_barrier(mpi_comm_world,ierr)
      !
      allocate(ES(1:num_l),EW(1:num_l),ED(1:num_l),E(1:num_l))
      !
      c_w1f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw,kmfftw])
      c_w2f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw,kmfftw])
      c_w3f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3f, w3f,  [imfftw,jmfftw,kmfftw])
      c_rho = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho, rhof,[imfftw,jmfftw,kmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11,[imfftw,jmfftw,kmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12,[imfftw,jmfftw,kmfftw])
      c_A13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13, A13,[imfftw,jmfftw,kmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21,[imfftw,jmfftw,kmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22,[imfftw,jmfftw,kmfftw])
      c_A23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23, A23,[imfftw,jmfftw,kmfftw])
      c_A31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31, A31,[imfftw,jmfftw,kmfftw])
      c_A32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32, A32,[imfftw,jmfftw,kmfftw])
      c_A33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33, A33,[imfftw,jmfftw,kmfftw])
      !
      ES=0.d0
      EW=0.d0
      ED=0.d0
      E=0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0)
          w1f(i,j,k)    = w1(i,j,k)    *Gl
          w2f(i,j,k)    = w2(i,j,k)    *Gl
          w3f(i,j,k)    = w3(i,j,k)    *Gl
          rhof(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(backward_plan,w3f,w3f)
        call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          w1f(i,j,k) = w1f(i,j,k)/rhof(i,j,k)
          w2f(i,j,k) = w2f(i,j,k)/rhof(i,j,k)
          w3f(i,j,k) = w3f(i,j,k)/rhof(i,j,k)
        enddo
        enddo
        enddo
        !
        call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(forward_plan,w3f,w3f)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          w1f(i,j,k)  = w1f(i,j,k)/(1.d0*ia*ja*ka)
          w2f(i,j,k)  = w2f(i,j,k)/(1.d0*ia*ja*ka)
          w3f(i,j,k)  = w3f(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11(i,j,k) = imag*w1f(i,j,k)*k1(i,j,k)
          A21(i,j,k) = imag*w2f(i,j,k)*k1(i,j,k)
          A31(i,j,k) = imag*w3f(i,j,k)*k1(i,j,k)
          A12(i,j,k) = imag*w1f(i,j,k)*k2(i,j,k)
          A22(i,j,k) = imag*w2f(i,j,k)*k2(i,j,k)
          A32(i,j,k) = imag*w3f(i,j,k)*k2(i,j,k)
          A13(i,j,k) = imag*w1f(i,j,k)*k3(i,j,k)
          A23(i,j,k) = imag*w2f(i,j,k)*k3(i,j,k)
          A33(i,j,k) = imag*w3f(i,j,k)*k3(i,j,k)
        end do
        end do
        end do
        !
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A31,A31)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        call fftw_mpi_execute_dft(backward_plan,A32,A32)
        call fftw_mpi_execute_dft(backward_plan,A13,A13)
        call fftw_mpi_execute_dft(backward_plan,A23,A23)
        call fftw_mpi_execute_dft(backward_plan,A33,A33)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          All(i,j,k) = dreal(A11(i,j,k)+A22(i,j,k)+A33(i,j,k))
          S11(i,j,k) = dreal(A11(i,j,k)) - 1.d0/3.d0 * All(i,j,k)
          S22(i,j,k) = dreal(A22(i,j,k)) - 1.d0/3.d0 * All(i,j,k)
          S33(i,j,k) = dreal(A33(i,j,k)) - 1.d0/3.d0 * All(i,j,k)
          S12(i,j,k) = dreal(A12(i,j,k) + A21(i,j,k))*0.5d0
          S21(i,j,k) = S12(i,j,k)
          S13(i,j,k) = dreal(A13(i,j,k) + A31(i,j,k))*0.5d0
          S31(i,j,k) = S13(i,j,k)
          S23(i,j,k) = dreal(A23(i,j,k) + A32(i,j,k))*0.5d0
          S32(i,j,k) = S23(i,j,k)
          W12(i,j,k) = dreal(A12(i,j,k)-A21(i,j,k))*0.5d0
          W21(i,j,k) = -1.d0*W12(i,j,k)
          W13(i,j,k) = dreal(A13(i,j,k)-A31(i,j,k))*0.5d0
          W31(i,j,k) = -1.d0*W13(i,j,k)
          W23(i,j,k) = dreal(A23(i,j,k)-A32(i,j,k))*0.5d0
          W32(i,j,k) = -1.d0*W23(i,j,k)
        end do
        end do
        end do
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          ES(m) = ES(m) + dreal(rhof(i,j,k)) * (S11(i,j,k)*S11(i,j,k) + S12(i,j,k)*S12(i,j,k) + S13(i,j,k)*S13(i,j,k) + &
                            S21(i,j,k)*S21(i,j,k) + S22(i,j,k)*S22(i,j,k) + S23(i,j,k)*S23(i,j,k) + &
                            S31(i,j,k)*S31(i,j,k) + S32(i,j,k)*S32(i,j,k) + S33(i,j,k)*S33(i,j,k)) / 2.d0
          EW(m) = EW(m) + dreal(rhof(i,j,k)) * (W12(i,j,k)*W12(i,j,k) + W13(i,j,k)*W13(i,j,k) + &
                            W21(i,j,k)*W21(i,j,k) + W23(i,j,k)*W23(i,j,k) + &
                            W31(i,j,k)*W31(i,j,k) + W32(i,j,k)*W32(i,j,k)) / 2.d0
          ED(m) = ED(m) + dreal(rhof(i,j,k)) * (All(i,j,k)*All(i,j,k)) / 6.d0
          E(m)  = E(m)  + dreal(rhof(i,j,k)) * (dreal(A11(i,j,k))*dreal(A11(i,j,k)) + dreal(A12(i,j,k))*dreal(A12(i,j,k)) + &
                            dreal(A13(i,j,k))*dreal(A13(i,j,k)) + dreal(A21(i,j,k))*dreal(A21(i,j,k)) + &
                            dreal(A22(i,j,k))*dreal(A22(i,j,k)) + dreal(A23(i,j,k))*dreal(A23(i,j,k)) + &
                            dreal(A31(i,j,k))*dreal(A31(i,j,k)) + dreal(A32(i,j,k))*dreal(A32(i,j,k)) + &
                            dreal(A33(i,j,k))*dreal(A33(i,j,k))) / 2.d0
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        ES(m) = psum(ES(m)) / (ia*ja*ka)
        EW(m) = psum(EW(m)) / (ia*ja*ka)
        ED(m) = psum(ED(m)) / (ia*ja*ka)
        E(m)  = psum(E(m))  / (ia*ja*ka)
      enddo
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
      if (thefilenumb .ne. 0) then
        outfilename = 'pp/SGS_E_'//stepname//'.dat'
      else
        outfilename = 'pp/SGS_E.dat'
      endif
      call listinit(filename=outfilename,handle=hand_a, &
              firstline='nstep time ell ES EW ED E')
      do m=1,num_l
        call listwrite(hand_a,l_lim(m), ES(m), EW(m), ED(m), E(m))
      enddo
      print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1f)
      call fftw_free(c_w2f)
      call fftw_free(c_w3f)
      call fftw_free(c_rho)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A13)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_A23)
      call fftw_free(c_A31)
      call fftw_free(c_A32)
      call fftw_free(c_A33)
      call mpistop
      deallocate(k1,k2,k3)
      deallocate(l_lim)
      deallocate(All,S11,S12,S13,S21,S22,S23,S31,S32,S33)
      deallocate(W12,W13,W21,W23,W31,W32)
      deallocate(ES,EW,ED,E)
      !
    end subroutine SGSE3D
    !
    subroutine SGSET3D(thefilenumb)
      ! ! TODO : Improve need / Test need
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh,dim
      integer :: i,j,k,m,n,mmm
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_aS,hand_aW,hand_aD,hand_bS,hand_bW,hand_bD
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      real(8), allocatable, dimension(:) :: Pi1S,Pi2S,Pi3S,Pi4S,Pi5S, &
                                            Pi1W,Pi2W,Pi3W,Pi4W,Pi5W, &
                                            Pi1D,Pi2D,Pi3D,Pi4D,Pi5D
      real(8) ::  Pi1Sint,Pi2Sint,Pi3Sint,Pi4Sint,Pi5Sint,&
                  Pi1Wint,Pi2Wint,Pi3Wint,Pi4Wint,Pi5Wint,&
                  Pi1Dint,Pi2Dint,Pi3Dint,Pi4Dint,Pi5Dint
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: S1mm1_filted_l, S1mm2_filted_l, S1mm3_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: S2mm1_filted_l, S2mm2_filted_l, S2mm3_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: S3mm1_filted_l, S3mm2_filted_l, S3mm3_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: W1mm1_filted_l, W1mm2_filted_l, W1mm3_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: W2mm1_filted_l, W2mm2_filted_l, W2mm3_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: W3mm1_filted_l, W3mm2_filted_l, W3mm3_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: Amm11_filted_l, Amm12_filted_l, Amm13_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: Amm21_filted_l, Amm22_filted_l, Amm23_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: Amm31_filted_l, Amm32_filted_l, Amm33_filted_l
      complex(8), allocatable, dimension(:,:,:) :: All_filted
      complex(8), allocatable, dimension(:,:,:) :: S11_filted,S12_filted,S13_filted
      complex(8), allocatable, dimension(:,:,:) :: S21_filted,S22_filted,S23_filted
      complex(8), allocatable, dimension(:,:,:) :: S31_filted,S32_filted,S33_filted
      complex(8), allocatable, dimension(:,:,:) :: W12_filted,W21_filted
      complex(8), allocatable, dimension(:,:,:) :: W13_filted,W31_filted
      complex(8), allocatable, dimension(:,:,:) :: W23_filted,W32_filted
      !
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: termSS_11,termSS_12,termSS_13,&
                                                              termSS_21,termSS_22,termSS_23,&
                                                              termSS_31,termSS_32,termSS_33,&
                                                              termWW_11,termWW_12,termWW_13,&
                                                              termWW_21,termWW_22,termWW_23,&
                                                              termWW_31,termWW_32,termWW_33,&
                                                              termWS_11,termWS_12,termWS_13,&
                                                              termWS_21,termWS_22,termWS_23,&
                                                              termWS_31,termWS_32,termWS_33,&
                                                              termDD,&
                                                              termSD_11,termSD_12,termSD_13,&
                                                              termSD_21,termSD_22,termSD_23,&
                                                              termSD_31,termSD_32,termSD_33
      real(8) :: vxr_D1S,vxr_D2S,vxr_D3S,vxr_D4S,vxr_D5S, &
                vxr_D1W,vxr_D2W,vxr_D3W,vxr_D4W,vxr_D5W, &
                vxr_D1D,vxr_D2D,vxr_D3D,vxr_D4D,vxr_D5D
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      type(C_PTR) :: c_termSS_11,c_termSS_12,c_termSS_13,c_termSS_21,c_termSS_22,c_termSS_23,c_termSS_31,c_termSS_32,c_termSS_33
      type(C_PTR) :: c_termWW_11,c_termWW_12,c_termWW_13,c_termWW_21,c_termWW_22,c_termWW_23,c_termWW_31,c_termWW_32,c_termWW_33
      type(C_PTR) :: c_termWS_11,c_termWS_12,c_termWS_13,c_termWS_21,c_termWS_22,c_termWS_23,c_termWS_31,c_termWS_32,c_termWS_33
      type(C_PTR) :: c_termDD
      type(C_PTR) :: c_termSD_11,c_termSD_12,c_termSD_13,c_termSD_21,c_termSD_22,c_termSD_23,c_termSD_31,c_termSD_32,c_termSD_33
      type(C_PTR) :: c_S1mm1_filted_l, c_S1mm2_filted_l, c_S1mm3_filted_l
      type(C_PTR) :: c_S2mm1_filted_l, c_S2mm2_filted_l, c_S2mm3_filted_l
      type(C_PTR) :: c_S3mm1_filted_l, c_S3mm2_filted_l, c_S3mm3_filted_l
      type(C_PTR) :: c_W1mm1_filted_l, c_W1mm2_filted_l, c_W1mm3_filted_l
      type(C_PTR) :: c_W2mm1_filted_l, c_W2mm2_filted_l, c_W2mm3_filted_l
      type(C_PTR) :: c_W3mm1_filted_l, c_W3mm2_filted_l, c_W3mm3_filted_l
      type(C_PTR) :: c_Amm11_filted_l, c_Amm12_filted_l, c_Amm13_filted_l
      type(C_PTR) :: c_Amm21_filted_l, c_Amm22_filted_l, c_Amm23_filted_l
      type(C_PTR) :: c_Amm31_filted_l, c_Amm32_filted_l, c_Amm33_filted_l
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      dim = 3
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi1S(1:num_l), Pi2S(1:num_l), Pi3S(1:num_l), Pi4S(1:num_l), Pi5S(1:num_l))
      allocate(Pi1W(1:num_l), Pi2W(1:num_l), Pi3W(1:num_l), Pi4W(1:num_l), Pi5W(1:num_l))
      allocate(Pi1D(1:num_l), Pi2D(1:num_l), Pi3D(1:num_l), Pi4D(1:num_l), Pi5D(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      c_S1mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S1mm1_filted_l, S1mm1_filted_l, [imfftw,jmfftw,kmfftw])
      c_S1mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S1mm2_filted_l, S1mm2_filted_l, [imfftw,jmfftw,kmfftw])
      c_S1mm3_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S1mm3_filted_l, S1mm3_filted_l, [imfftw,jmfftw,kmfftw])
      c_S2mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S2mm1_filted_l, S2mm1_filted_l, [imfftw,jmfftw,kmfftw])
      c_S2mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S2mm2_filted_l, S2mm2_filted_l, [imfftw,jmfftw,kmfftw])
      c_S2mm3_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S2mm3_filted_l, S2mm3_filted_l, [imfftw,jmfftw,kmfftw])
      c_S3mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S3mm1_filted_l, S3mm1_filted_l, [imfftw,jmfftw,kmfftw])
      c_S3mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S3mm2_filted_l, S3mm2_filted_l, [imfftw,jmfftw,kmfftw])
      c_S3mm3_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_S3mm3_filted_l, S3mm3_filted_l, [imfftw,jmfftw,kmfftw])
      !
      c_W1mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W1mm1_filted_l, W1mm1_filted_l, [imfftw,jmfftw,kmfftw])
      c_W1mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W1mm2_filted_l, W1mm2_filted_l, [imfftw,jmfftw,kmfftw])
      c_W1mm3_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W1mm3_filted_l, W1mm3_filted_l, [imfftw,jmfftw,kmfftw])
      c_W2mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W2mm1_filted_l, W2mm1_filted_l, [imfftw,jmfftw,kmfftw])
      c_W2mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W2mm2_filted_l, W2mm2_filted_l, [imfftw,jmfftw,kmfftw])
      c_W2mm3_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W2mm3_filted_l, W2mm3_filted_l, [imfftw,jmfftw,kmfftw])
      c_W3mm1_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W3mm1_filted_l, W3mm1_filted_l, [imfftw,jmfftw,kmfftw])
      c_W3mm2_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W3mm2_filted_l, W3mm2_filted_l, [imfftw,jmfftw,kmfftw])
      c_W3mm3_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_W3mm3_filted_l, W3mm3_filted_l, [imfftw,jmfftw,kmfftw])
      !
      c_Amm11_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm11_filted_l, Amm11_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm12_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm12_filted_l, Amm12_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm13_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm13_filted_l, Amm13_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm21_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm21_filted_l, Amm21_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm22_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm22_filted_l, Amm22_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm23_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm23_filted_l, Amm23_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm31_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm31_filted_l, Amm31_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm32_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm32_filted_l, Amm32_filted_l, [imfftw,jmfftw,kmfftw])
      c_Amm33_filted_l = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_Amm33_filted_l, Amm33_filted_l, [imfftw,jmfftw,kmfftw])
      !
      allocate(All_filted(1:im,1:jm,1:km),&
              S11_filted(1:im,1:jm,1:km),S12_filted(1:im,1:jm,1:km),S13_filted(1:im,1:jm,1:km),&
              S21_filted(1:im,1:jm,1:km),S22_filted(1:im,1:jm,1:km),S23_filted(1:im,1:jm,1:km),&
              S31_filted(1:im,1:jm,1:km),S32_filted(1:im,1:jm,1:km),S33_filted(1:im,1:jm,1:km),&
              W12_filted(1:im,1:jm,1:km),W21_filted(1:im,1:jm,1:km),&
              W13_filted(1:im,1:jm,1:km),W31_filted(1:im,1:jm,1:km),&
              W23_filted(1:im,1:jm,1:km),W32_filted(1:im,1:jm,1:km))
      !
      c_termSS_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_11, termSS_11, [imfftw,jmfftw,kmfftw])
      c_termSS_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_12, termSS_12, [imfftw,jmfftw,kmfftw])
      c_termSS_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_13, termSS_13, [imfftw,jmfftw,kmfftw])
      c_termSS_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_21, termSS_21, [imfftw,jmfftw,kmfftw])
      c_termSS_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_22, termSS_22, [imfftw,jmfftw,kmfftw])
      c_termSS_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_23, termSS_23, [imfftw,jmfftw,kmfftw])
      c_termSS_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_31, termSS_31, [imfftw,jmfftw,kmfftw])
      c_termSS_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_32, termSS_32, [imfftw,jmfftw,kmfftw])
      c_termSS_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSS_33, termSS_33, [imfftw,jmfftw,kmfftw])
      c_termWW_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_11, termWW_11, [imfftw,jmfftw,kmfftw])
      c_termWW_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_12, termWW_12, [imfftw,jmfftw,kmfftw])
      c_termWW_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_13, termWW_13, [imfftw,jmfftw,kmfftw])
      c_termWW_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_21, termWW_21, [imfftw,jmfftw,kmfftw])
      c_termWW_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_22, termWW_22, [imfftw,jmfftw,kmfftw])
      c_termWW_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_23, termWW_23, [imfftw,jmfftw,kmfftw])
      c_termWW_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_31, termWW_31, [imfftw,jmfftw,kmfftw])
      c_termWW_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_32, termWW_32, [imfftw,jmfftw,kmfftw])
      c_termWW_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWW_33, termWW_33, [imfftw,jmfftw,kmfftw])
      c_termWS_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_11, termWS_11, [imfftw,jmfftw,kmfftw])
      c_termWS_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_12, termWS_12, [imfftw,jmfftw,kmfftw])
      c_termWS_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_13, termWS_13, [imfftw,jmfftw,kmfftw])
      c_termWS_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_21, termWS_21, [imfftw,jmfftw,kmfftw])
      c_termWS_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_22, termWS_22, [imfftw,jmfftw,kmfftw])
      c_termWS_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_23, termWS_23, [imfftw,jmfftw,kmfftw])
      c_termWS_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_31, termWS_31, [imfftw,jmfftw,kmfftw])
      c_termWS_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_32, termWS_32, [imfftw,jmfftw,kmfftw])
      c_termWS_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termWS_33, termWS_33, [imfftw,jmfftw,kmfftw])
      c_termSD_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_11, termSD_11, [imfftw,jmfftw,kmfftw])
      c_termSD_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_12, termSD_12, [imfftw,jmfftw,kmfftw])
      c_termSD_13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_13, termSD_13, [imfftw,jmfftw,kmfftw])
      c_termSD_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_21, termSD_21, [imfftw,jmfftw,kmfftw])
      c_termSD_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_22, termSD_22, [imfftw,jmfftw,kmfftw])
      c_termSD_23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_23, termSD_23, [imfftw,jmfftw,kmfftw])
      c_termSD_31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_31, termSD_31, [imfftw,jmfftw,kmfftw])
      c_termSD_32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_32, termSD_32, [imfftw,jmfftw,kmfftw])
      c_termSD_33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termSD_33, termSD_33, [imfftw,jmfftw,kmfftw])
      c_termDD = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termDD, termDD, [imfftw,jmfftw,kmfftw])
      !
      !
      Pi1S =0.d0
      Pi2S =0.d0
      Pi3S =0.d0
      Pi4S =0.d0
      Pi5S =0.d0
      Pi1W =0.d0
      Pi2W =0.d0
      Pi3W =0.d0
      Pi4W =0.d0
      Pi5W =0.d0
      Pi1D =0.d0
      Pi2D =0.d0
      Pi3D =0.d0
      Pi4D =0.d0
      Pi5D =0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
            if (thefilenumb .ne. 0) then
              outfilename2 = 'pp/SGS_ETS_precise_'//stepname//'_'//mname//'.dat'
            else
              outfilename2 = 'pp/SGS_ETS_precise_'//mname//'.dat'
            endif
            call listinit(filename=outfilename2,handle=hand_bS, &
                  firstline='nstep time sqrtalpha pi1S pi2S pi3S pi4S pi5S')
            if (thefilenumb .ne. 0) then
              outfilename2 = 'pp/SGS_ETW_precise_'//stepname//'_'//mname//'.dat'
            else
              outfilename2 = 'pp/SGS_ETW_precise_'//mname//'.dat'
            endif
            call listinit(filename=outfilename2,handle=hand_bW, &
                  firstline='nstep time sqrtalpha pi1W pi2W pi3W pi4W pi5W')
            if (thefilenumb .ne. 0) then
              outfilename2 = 'pp/SGS_ETD_precise_'//stepname//'_'//mname//'.dat'
            else
              outfilename2 = 'pp/SGS_ETD_precise_'//mname//'.dat'
            endif
            call listinit(filename=outfilename2,handle=hand_bD, &
                  firstline='nstep time sqrtalpha pi1D pi2D pi3D pi4D pi5D')
        endif
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j,k)    = w1(i,j,k)    *Gl
          w2_filted(i,j,k)    = w2(i,j,k)    *Gl
          w3_filted(i,j,k)    = w3(i,j,k)    *Gl
          !
          rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
          All_filted(i,j,k) = (A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k))
          !
          S11_filted(i,j,k) = (A11_filted(i,j,k)) - 1.d0/real(dim) * All_filted(i,j,k)
          S22_filted(i,j,k) = (A22_filted(i,j,k)) - 1.d0/real(dim) * All_filted(i,j,k)
          S33_filted(i,j,k) = (A33_filted(i,j,k)) - 1.d0/real(dim) * All_filted(i,j,k)
          S12_filted(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
          S21_filted(i,j,k) = S12_filted(i,j,k)
          S13_filted(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
          S31_filted(i,j,k) = S13_filted(i,j,k)
          S23_filted(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
          S32_filted(i,j,k) = S23_filted(i,j,k)
          !
          W12_filted(i,j,k) = (A12_filted(i,j,k) - A21_filted(i,j,k))*0.5d0
          W21_filted(i,j,k) = -1.d0 * W12_filted(i,j,k)
          W13_filted(i,j,k) = (A13_filted(i,j,k) - A31_filted(i,j,k))*0.5d0
          W31_filted(i,j,k) = -1.d0 * W13_filted(i,j,k)
          W23_filted(i,j,k) = (A23_filted(i,j,k) - A32_filted(i,j,k))*0.5d0
          W32_filted(i,j,k) = -1.d0 * W23_filted(i,j,k)
          !
          ! 
          !SAmmB_filted(i,j,k) = - k1(i,j,k)*KB(i,j,k)*SA1_filted(i,j,k) - k2(i,j,k)*KB(i,j,k)*SA2_filted(i,j,k) - k3(i,j,k)*KB(i,j,k)*SA3_filted(i,j,k)
          !
          S1mm1_filted_l(i,j,k) = - k1(i,j,k)*k1(i,j,k)*S11_filted(i,j,k) - k2(i,j,k)*k1(i,j,k)*S12_filted(i,j,k) &
          - k3(i,j,k)*k1(i,j,k)*S13_filted(i,j,k)
          S1mm2_filted_l(i,j,k) = - k1(i,j,k)*k2(i,j,k)*S11_filted(i,j,k) - k2(i,j,k)*k2(i,j,k)*S12_filted(i,j,k) &
          - k3(i,j,k)*k2(i,j,k)*S13_filted(i,j,k)
          S1mm3_filted_l(i,j,k) = - k1(i,j,k)*k3(i,j,k)*S11_filted(i,j,k) - k2(i,j,k)*k3(i,j,k)*S12_filted(i,j,k) &
          - k3(i,j,k)*k3(i,j,k)*S13_filted(i,j,k)
          S2mm1_filted_l(i,j,k) = - k1(i,j,k)*k1(i,j,k)*S21_filted(i,j,k) - k2(i,j,k)*k1(i,j,k)*S22_filted(i,j,k) &
          - k3(i,j,k)*k1(i,j,k)*S23_filted(i,j,k)
          S2mm2_filted_l(i,j,k) = - k1(i,j,k)*k2(i,j,k)*S21_filted(i,j,k) - k2(i,j,k)*k2(i,j,k)*S22_filted(i,j,k) &
          - k3(i,j,k)*k2(i,j,k)*S23_filted(i,j,k)
          S2mm3_filted_l(i,j,k) = - k1(i,j,k)*k3(i,j,k)*S21_filted(i,j,k) - k2(i,j,k)*k3(i,j,k)*S22_filted(i,j,k) &
          - k3(i,j,k)*k3(i,j,k)*S23_filted(i,j,k)
          S3mm1_filted_l(i,j,k) = - k1(i,j,k)*k1(i,j,k)*S31_filted(i,j,k) - k2(i,j,k)*k1(i,j,k)*S32_filted(i,j,k) &
          - k3(i,j,k)*k1(i,j,k)*S33_filted(i,j,k)
          S3mm2_filted_l(i,j,k) = - k1(i,j,k)*k2(i,j,k)*S31_filted(i,j,k) - k2(i,j,k)*k2(i,j,k)*S32_filted(i,j,k) &
          - k3(i,j,k)*k2(i,j,k)*S33_filted(i,j,k)
          S3mm3_filted_l(i,j,k) = - k1(i,j,k)*k3(i,j,k)*S31_filted(i,j,k) - k2(i,j,k)*k3(i,j,k)*S32_filted(i,j,k) &
          - k3(i,j,k)*k3(i,j,k)*S33_filted(i,j,k)
          !
          W1mm1_filted_l(i,j,k) = - k2(i,j,k)*k1(i,j,k)*W12_filted(i,j,k) - k3(i,j,k)*k1(i,j,k)*W13_filted(i,j,k)
          W1mm2_filted_l(i,j,k) = - k2(i,j,k)*k2(i,j,k)*W12_filted(i,j,k) - k3(i,j,k)*k2(i,j,k)*W13_filted(i,j,k)
          W1mm3_filted_l(i,j,k) = - k2(i,j,k)*k3(i,j,k)*W12_filted(i,j,k) - k3(i,j,k)*k3(i,j,k)*W13_filted(i,j,k)
          W2mm1_filted_l(i,j,k) = - k1(i,j,k)*k1(i,j,k)*W21_filted(i,j,k) - k3(i,j,k)*k1(i,j,k)*W23_filted(i,j,k)
          W2mm2_filted_l(i,j,k) = - k1(i,j,k)*k2(i,j,k)*W21_filted(i,j,k) - k3(i,j,k)*k2(i,j,k)*W23_filted(i,j,k)
          W2mm3_filted_l(i,j,k) = - k1(i,j,k)*k3(i,j,k)*W21_filted(i,j,k) - k3(i,j,k)*k3(i,j,k)*W23_filted(i,j,k)
          W3mm1_filted_l(i,j,k) = - k1(i,j,k)*k1(i,j,k)*W31_filted(i,j,k) - k2(i,j,k)*k1(i,j,k)*W32_filted(i,j,k) 
          W3mm2_filted_l(i,j,k) = - k1(i,j,k)*k2(i,j,k)*W31_filted(i,j,k) - k2(i,j,k)*k2(i,j,k)*W32_filted(i,j,k) 
          W3mm3_filted_l(i,j,k) = - k1(i,j,k)*k3(i,j,k)*W31_filted(i,j,k) - k2(i,j,k)*k3(i,j,k)*W32_filted(i,j,k)
          !
          Amm11_filted_l(i,j,k) = - k1(i,j,k)*k1(i,j,k)*All_filted(i,j,k)
          Amm12_filted_l(i,j,k) = - k1(i,j,k)*k2(i,j,k)*All_filted(i,j,k)
          Amm13_filted_l(i,j,k) = - k1(i,j,k)*k3(i,j,k)*All_filted(i,j,k)
          Amm21_filted_l(i,j,k) = - k2(i,j,k)*k1(i,j,k)*All_filted(i,j,k)
          Amm22_filted_l(i,j,k) = - k2(i,j,k)*k2(i,j,k)*All_filted(i,j,k)
          Amm23_filted_l(i,j,k) = - k2(i,j,k)*k3(i,j,k)*All_filted(i,j,k)
          Amm31_filted_l(i,j,k) = - k3(i,j,k)*k1(i,j,k)*All_filted(i,j,k)
          Amm32_filted_l(i,j,k) = - k3(i,j,k)*k2(i,j,k)*All_filted(i,j,k)
          Amm33_filted_l(i,j,k) = - k3(i,j,k)*k3(i,j,k)*All_filted(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,S1mm1_filted_l,S1mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S1mm2_filted_l,S1mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S1mm3_filted_l,S1mm3_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S2mm1_filted_l,S2mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S2mm2_filted_l,S2mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S2mm3_filted_l,S2mm3_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S3mm1_filted_l,S3mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S3mm2_filted_l,S3mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,S3mm3_filted_l,S3mm3_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W1mm1_filted_l,W1mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W1mm2_filted_l,W1mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W1mm3_filted_l,W1mm3_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W2mm1_filted_l,W2mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W2mm2_filted_l,W2mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W2mm3_filted_l,W2mm3_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W3mm1_filted_l,W3mm1_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W3mm2_filted_l,W3mm2_filted_l)
        call fftw_mpi_execute_dft(backward_plan,W3mm3_filted_l,W3mm3_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm11_filted_l,Amm11_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm12_filted_l,Amm12_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm13_filted_l,Amm13_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm21_filted_l,Amm21_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm22_filted_l,Amm22_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm23_filted_l,Amm23_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm31_filted_l,Amm31_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm32_filted_l,Amm32_filted_l)
        call fftw_mpi_execute_dft(backward_plan,Amm33_filted_l,Amm33_filted_l)
        !
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !!!! Velocity Favre average and density average
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          do k=1,km
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j,k)  = w1(i,j,k)    *Galpha
            w2_filted(i,j,k)  = w2(i,j,k)    *Galpha
            w3_filted(i,j,k)  = w3(i,j,k)    *Galpha
            rho_filted(i,j,k) = rhocom(i,j,k)*Galpha
          enddo
          enddo
          enddo
          !
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          do k=1,km
          do j=1,jm
          do i=1,im
            w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
            w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
            w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          enddo
          enddo
          enddo
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
            w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
            w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
            !
            A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
            A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
            A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
            A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
            A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
            A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
            A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
            A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
            A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
            !
          end do
          end do
          end do
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
          call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
          call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
          call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
          call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
          call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
          call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
          call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
          call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
          !
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            All_filted(i,j,k) = dreal(A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k))
            !
            S11_filted(i,j,k) = dreal(A11_filted(i,j,k)) - 1.d0/real(dim) * All_filted(i,j,k)
            S22_filted(i,j,k) = dreal(A22_filted(i,j,k)) - 1.d0/real(dim) * All_filted(i,j,k)
            S33_filted(i,j,k) = dreal(A33_filted(i,j,k)) - 1.d0/real(dim) * All_filted(i,j,k)
            S12_filted(i,j,k) = dreal(A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
            S21_filted(i,j,k) = S12_filted(i,j,k)
            S13_filted(i,j,k) = dreal(A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
            S31_filted(i,j,k) = S13_filted(i,j,k)
            S23_filted(i,j,k) = dreal(A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
            S32_filted(i,j,k) = S23_filted(i,j,k)
            !
            W12_filted(i,j,k) = dreal(A12_filted(i,j,k)-A21_filted(i,j,k))*0.5d0
            W21_filted(i,j,k) = -1.d0*W12_filted(i,j,k)
            W13_filted(i,j,k) = dreal(A13_filted(i,j,k)-A31_filted(i,j,k))*0.5d0
            W31_filted(i,j,k) = -1.d0*W13_filted(i,j,k)
            W23_filted(i,j,k) = dreal(A23_filted(i,j,k)-A32_filted(i,j,k))*0.5d0
            W32_filted(i,j,k) = -1.d0*W23_filted(i,j,k)
            !
          end do
          end do
          end do
          !
          !!!! Pi terms
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            rho_filted(i,j,k) = dreal(rho_filted(i,j,k))
            !termSS_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted + rho_filted*SI3_filted*SJ3_filted
            termSS_11(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S11_filted(i,j,k) + &
                              rho_filted(i,j,k)*S12_filted(i,j,k)*S12_filted(i,j,k) + &
                              rho_filted(i,j,k)*S13_filted(i,j,k)*S13_filted(i,j,k)
            termSS_12(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S21_filted(i,j,k) + &
                              rho_filted(i,j,k)*S12_filted(i,j,k)*S22_filted(i,j,k) + &
                              rho_filted(i,j,k)*S13_filted(i,j,k)*S23_filted(i,j,k)
            termSS_13(i,j,k) = rho_filted(i,j,k)*S11_filted(i,j,k)*S31_filted(i,j,k) + &
                              rho_filted(i,j,k)*S12_filted(i,j,k)*S32_filted(i,j,k) + &
                              rho_filted(i,j,k)*S13_filted(i,j,k)*S33_filted(i,j,k)
            termSS_21(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S11_filted(i,j,k) + &
                              rho_filted(i,j,k)*S22_filted(i,j,k)*S12_filted(i,j,k) + &
                              rho_filted(i,j,k)*S23_filted(i,j,k)*S13_filted(i,j,k)
            termSS_22(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S21_filted(i,j,k) + &
                              rho_filted(i,j,k)*S22_filted(i,j,k)*S22_filted(i,j,k) + &
                              rho_filted(i,j,k)*S23_filted(i,j,k)*S23_filted(i,j,k)
            termSS_23(i,j,k) = rho_filted(i,j,k)*S21_filted(i,j,k)*S31_filted(i,j,k) + &
                              rho_filted(i,j,k)*S22_filted(i,j,k)*S32_filted(i,j,k) + &
                              rho_filted(i,j,k)*S23_filted(i,j,k)*S33_filted(i,j,k)
            termSS_31(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S11_filted(i,j,k) + &
                              rho_filted(i,j,k)*S32_filted(i,j,k)*S12_filted(i,j,k) + &
                              rho_filted(i,j,k)*S33_filted(i,j,k)*S13_filted(i,j,k)
            termSS_32(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S21_filted(i,j,k) + &
                              rho_filted(i,j,k)*S32_filted(i,j,k)*S22_filted(i,j,k) + &
                              rho_filted(i,j,k)*S33_filted(i,j,k)*S23_filted(i,j,k)
            termSS_33(i,j,k) = rho_filted(i,j,k)*S31_filted(i,j,k)*S31_filted(i,j,k) + &
                              rho_filted(i,j,k)*S32_filted(i,j,k)*S32_filted(i,j,k) + &
                              rho_filted(i,j,k)*S33_filted(i,j,k)*S33_filted(i,j,k)
            ! 
            ! termSD_IJ = rho_filted*All_filted*SIJ_filted
            termSD_11(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S11_filted(i,j,k)
            termSD_12(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S12_filted(i,j,k)
            termSD_13(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S13_filted(i,j,k)
            termSD_21(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S21_filted(i,j,k)
            termSD_22(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S22_filted(i,j,k)
            termSD_23(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S23_filted(i,j,k)
            termSD_31(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S31_filted(i,j,k)
            termSD_32(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S32_filted(i,j,k)
            termSD_33(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*S33_filted(i,j,k)
            !
            !termWW_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted + &
            !rho_filted*WI3_filted*W3J_filted 
            termWW_11(i,j,k) = rho_filted(i,j,k)*W12_filted(i,j,k)*W21_filted(i,j,k) + &
                              rho_filted(i,j,k)*W13_filted(i,j,k)*W31_filted(i,j,k) 
            termWW_21(i,j,k) = rho_filted(i,j,k)*W23_filted(i,j,k)*W31_filted(i,j,k) 
            termWW_31(i,j,k) = rho_filted(i,j,k)*W32_filted(i,j,k)*W21_filted(i,j,k)
            !
            termWW_12(i,j,k) = rho_filted(i,j,k)*W13_filted(i,j,k)*W32_filted(i,j,k) 
            termWW_22(i,j,k) = rho_filted(i,j,k)*W21_filted(i,j,k)*W12_filted(i,j,k) + &
                              rho_filted(i,j,k)*W23_filted(i,j,k)*W32_filted(i,j,k) 
            termWW_32(i,j,k) = rho_filted(i,j,k)*W31_filted(i,j,k)*W12_filted(i,j,k)
            !
            termWW_13(i,j,k) = rho_filted(i,j,k)*W12_filted(i,j,k)*W23_filted(i,j,k)
            termWW_23(i,j,k) = rho_filted(i,j,k)*W21_filted(i,j,k)*W13_filted(i,j,k)
            termWW_33(i,j,k) = rho_filted(i,j,k)*W31_filted(i,j,k)*W13_filted(i,j,k) + &
                              rho_filted(i,j,k)*W32_filted(i,j,k)*W23_filted(i,j,k)
            !
            !
            termWS_11(i,j,k)= rho_filted(i,j,k)*S21_filted(i,j,k)*W12_filted(i,j,k) &
                            -rho_filted(i,j,k)*S12_filted(i,j,k)*W21_filted(i,j,k) &
                            +rho_filted(i,j,k)*S31_filted(i,j,k)*W13_filted(i,j,k) &
                            -rho_filted(i,j,k)*S13_filted(i,j,k)*W31_filted(i,j,k)
            termWS_21(i,j,k)= rho_filted(i,j,k)*S11_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S22_filted(i,j,k)*W21_filted(i,j,k) &
                            +rho_filted(i,j,k)*S31_filted(i,j,k)*W23_filted(i,j,k) &
                            -rho_filted(i,j,k)*S23_filted(i,j,k)*W31_filted(i,j,k)
            termWS_31(i,j,k)= rho_filted(i,j,k)*S11_filted(i,j,k)*W31_filted(i,j,k) &
                            +rho_filted(i,j,k)*S21_filted(i,j,k)*W32_filted(i,j,k) &
                            -rho_filted(i,j,k)*S32_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S33_filted(i,j,k)*W31_filted(i,j,k) 
            !
            termWS_12(i,j,k)=-rho_filted(i,j,k)*S11_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S22_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S32_filted(i,j,k)*W13_filted(i,j,k) &
                            -rho_filted(i,j,k)*S13_filted(i,j,k)*W32_filted(i,j,k)
            termWS_22(i,j,k)= rho_filted(i,j,k)*S12_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S21_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S32_filted(i,j,k)*W23_filted(i,j,k) &
                            -rho_filted(i,j,k)*S23_filted(i,j,k)*W32_filted(i,j,k)
            termWS_32(i,j,k)= rho_filted(i,j,k)*S12_filted(i,j,k)*W31_filted(i,j,k) &
                            -rho_filted(i,j,k)*S31_filted(i,j,k)*W12_filted(i,j,k) &
                            +rho_filted(i,j,k)*S22_filted(i,j,k)*W32_filted(i,j,k) &
                            -rho_filted(i,j,k)*S33_filted(i,j,k)*W32_filted(i,j,k)
            !
            termWS_13(i,j,k)=-rho_filted(i,j,k)*S11_filted(i,j,k)*W13_filted(i,j,k) &
                            +rho_filted(i,j,k)*S23_filted(i,j,k)*W12_filted(i,j,k) &
                            -rho_filted(i,j,k)*S12_filted(i,j,k)*W23_filted(i,j,k) &
                            +rho_filted(i,j,k)*S33_filted(i,j,k)*W13_filted(i,j,k)
            termWS_23(i,j,k)= rho_filted(i,j,k)*S13_filted(i,j,k)*W21_filted(i,j,k) &
                            -rho_filted(i,j,k)*S21_filted(i,j,k)*W13_filted(i,j,k) &
                            -rho_filted(i,j,k)*S22_filted(i,j,k)*W23_filted(i,j,k) &
                            +rho_filted(i,j,k)*S33_filted(i,j,k)*W23_filted(i,j,k)
            termWS_33(i,j,k)= rho_filted(i,j,k)*S13_filted(i,j,k)*W31_filted(i,j,k) &
                            -rho_filted(i,j,k)*S31_filted(i,j,k)*W13_filted(i,j,k) &
                            +rho_filted(i,j,k)*S23_filted(i,j,k)*W32_filted(i,j,k) &
                            -rho_filted(i,j,k)*S32_filted(i,j,k)*W23_filted(i,j,k)
            !
            ! termDD
            termDD(i,j,k) = rho_filted(i,j,k)*All_filted(i,j,k)*All_filted(i,j,k)
          enddo
          enddo
          enddo
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,termSS_11,termSS_11)
          call fftw_mpi_execute_dft(forward_plan,termSS_12,termSS_12)
          call fftw_mpi_execute_dft(forward_plan,termSS_13,termSS_13)
          call fftw_mpi_execute_dft(forward_plan,termSS_21,termSS_21)
          call fftw_mpi_execute_dft(forward_plan,termSS_22,termSS_22)
          call fftw_mpi_execute_dft(forward_plan,termSS_23,termSS_23)
          call fftw_mpi_execute_dft(forward_plan,termSS_31,termSS_31)
          call fftw_mpi_execute_dft(forward_plan,termSS_32,termSS_32)
          call fftw_mpi_execute_dft(forward_plan,termSS_33,termSS_33)
          !
          call fftw_mpi_execute_dft(forward_plan,termWW_11,termWW_11)
          call fftw_mpi_execute_dft(forward_plan,termWW_12,termWW_12)
          call fftw_mpi_execute_dft(forward_plan,termWW_13,termWW_13)
          call fftw_mpi_execute_dft(forward_plan,termWW_21,termWW_21)
          call fftw_mpi_execute_dft(forward_plan,termWW_22,termWW_22)
          call fftw_mpi_execute_dft(forward_plan,termWW_23,termWW_23)
          call fftw_mpi_execute_dft(forward_plan,termWW_31,termWW_31)
          call fftw_mpi_execute_dft(forward_plan,termWW_32,termWW_32)
          call fftw_mpi_execute_dft(forward_plan,termWW_33,termWW_33)
          !
          call fftw_mpi_execute_dft(forward_plan,termWS_11,termWS_11)
          call fftw_mpi_execute_dft(forward_plan,termWS_12,termWS_12)
          call fftw_mpi_execute_dft(forward_plan,termWS_13,termWS_13)
          call fftw_mpi_execute_dft(forward_plan,termWS_21,termWS_21)
          call fftw_mpi_execute_dft(forward_plan,termWS_22,termWS_22)
          call fftw_mpi_execute_dft(forward_plan,termWS_23,termWS_23)
          call fftw_mpi_execute_dft(forward_plan,termWS_31,termWS_31)
          call fftw_mpi_execute_dft(forward_plan,termWS_32,termWS_32)
          call fftw_mpi_execute_dft(forward_plan,termWS_33,termWS_33)
          ! !
          ! call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
          ! !
          call fftw_mpi_execute_dft(forward_plan,termSD_11,termSD_11)
          call fftw_mpi_execute_dft(forward_plan,termSD_12,termSD_12)
          call fftw_mpi_execute_dft(forward_plan,termSD_13,termSD_13)
          call fftw_mpi_execute_dft(forward_plan,termSD_21,termSD_21)
          call fftw_mpi_execute_dft(forward_plan,termSD_22,termSD_22)
          call fftw_mpi_execute_dft(forward_plan,termSD_23,termSD_23)
          call fftw_mpi_execute_dft(forward_plan,termSD_31,termSD_31)
          call fftw_mpi_execute_dft(forward_plan,termSD_32,termSD_32)
          call fftw_mpi_execute_dft(forward_plan,termSD_33,termSD_33)
          !
          call fftw_mpi_execute_dft(forward_plan,termDD   ,termDD   )
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            termSS_11(i,j,k) = termSS_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_12(i,j,k) = termSS_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_13(i,j,k) = termSS_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_21(i,j,k) = termSS_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_22(i,j,k) = termSS_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_23(i,j,k) = termSS_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_31(i,j,k) = termSS_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_32(i,j,k) = termSS_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSS_33(i,j,k) = termSS_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            termWW_11(i,j,k) = termWW_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_12(i,j,k) = termWW_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_13(i,j,k) = termWW_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_21(i,j,k) = termWW_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_22(i,j,k) = termWW_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_23(i,j,k) = termWW_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_31(i,j,k) = termWW_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_32(i,j,k) = termWW_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWW_33(i,j,k) = termWW_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            termWS_11(i,j,k) = termWS_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_12(i,j,k) = termWS_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_13(i,j,k) = termWS_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_21(i,j,k) = termWS_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_22(i,j,k) = termWS_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_23(i,j,k) = termWS_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_31(i,j,k) = termWS_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_32(i,j,k) = termWS_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termWS_33(i,j,k) = termWS_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            termSD_11(i,j,k) = termSD_11(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_12(i,j,k) = termSD_12(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_13(i,j,k) = termSD_13(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_21(i,j,k) = termSD_21(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_22(i,j,k) = termSD_22(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_23(i,j,k) = termSD_23(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_31(i,j,k) = termSD_31(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_32(i,j,k) = termSD_32(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            termSD_33(i,j,k) = termSD_33(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
            termDD(i,j,k)    = termDD(i,j,k)   *Gphi/(1.d0*ia*ja*ka)
            !
          enddo
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,termSS_11,termSS_11)
          call fftw_mpi_execute_dft(backward_plan,termSS_12,termSS_12)
          call fftw_mpi_execute_dft(backward_plan,termSS_13,termSS_13)
          call fftw_mpi_execute_dft(backward_plan,termSS_21,termSS_21)
          call fftw_mpi_execute_dft(backward_plan,termSS_22,termSS_22)
          call fftw_mpi_execute_dft(backward_plan,termSS_23,termSS_23)
          call fftw_mpi_execute_dft(backward_plan,termSS_31,termSS_31)
          call fftw_mpi_execute_dft(backward_plan,termSS_32,termSS_32)
          call fftw_mpi_execute_dft(backward_plan,termSS_33,termSS_33)
          !
          call fftw_mpi_execute_dft(backward_plan,termWW_11,termWW_11)
          call fftw_mpi_execute_dft(backward_plan,termWW_12,termWW_12)
          call fftw_mpi_execute_dft(backward_plan,termWW_13,termWW_13)
          call fftw_mpi_execute_dft(backward_plan,termWW_21,termWW_21)
          call fftw_mpi_execute_dft(backward_plan,termWW_22,termWW_22)
          call fftw_mpi_execute_dft(backward_plan,termWW_23,termWW_23)
          call fftw_mpi_execute_dft(backward_plan,termWW_31,termWW_31)
          call fftw_mpi_execute_dft(backward_plan,termWW_32,termWW_32)
          call fftw_mpi_execute_dft(backward_plan,termWW_33,termWW_33)
          !
          call fftw_mpi_execute_dft(backward_plan,termWS_11,termWS_11)
          call fftw_mpi_execute_dft(backward_plan,termWS_12,termWS_12)
          call fftw_mpi_execute_dft(backward_plan,termWS_13,termWS_13)
          call fftw_mpi_execute_dft(backward_plan,termWS_21,termWS_21)
          call fftw_mpi_execute_dft(backward_plan,termWS_22,termWS_22)
          call fftw_mpi_execute_dft(backward_plan,termWS_23,termWS_23)
          call fftw_mpi_execute_dft(backward_plan,termWS_31,termWS_31)
          call fftw_mpi_execute_dft(backward_plan,termWS_32,termWS_32)
          call fftw_mpi_execute_dft(backward_plan,termWS_33,termWS_33)
          !
          call fftw_mpi_execute_dft(backward_plan,termSD_11,termSD_11)
          call fftw_mpi_execute_dft(backward_plan,termSD_12,termSD_12)
          call fftw_mpi_execute_dft(backward_plan,termSD_13,termSD_13)
          call fftw_mpi_execute_dft(backward_plan,termSD_21,termSD_21)
          call fftw_mpi_execute_dft(backward_plan,termSD_22,termSD_22)
          call fftw_mpi_execute_dft(backward_plan,termSD_23,termSD_23)
          call fftw_mpi_execute_dft(backward_plan,termSD_31,termSD_31)
          call fftw_mpi_execute_dft(backward_plan,termSD_32,termSD_32)
          call fftw_mpi_execute_dft(backward_plan,termSD_33,termSD_33)
          !
          call fftw_mpi_execute_dft(backward_plan,termDD   ,termDD   )
          !
          !
          Pi1Sint = 0.d0
          Pi2Sint = 0.d0
          Pi3Sint = 0.d0
          Pi4Sint = 0.d0
          Pi5Sint = 0.d0
          Pi1Wint = 0.d0
          Pi2Wint = 0.d0
          Pi3Wint = 0.d0
          Pi4Wint = 0.d0
          Pi5Wint = 0.d0
          Pi1Dint = 0.d0
          Pi2Dint = 0.d0
          Pi3Dint = 0.d0
          Pi4Dint = 0.d0
          Pi5Dint = 0.d0
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            vxr_D1S = dreal(termSS_11(i,j,k) * S1mm1_filted_l(i,j,k) + &
                          termSS_12(i,j,k) * S1mm2_filted_l(i,j,k) + &
                          termSS_13(i,j,k) * S1mm3_filted_l(i,j,k) + &
                          termSS_21(i,j,k) * S2mm1_filted_l(i,j,k) + &
                          termSS_22(i,j,k) * S2mm2_filted_l(i,j,k) + &
                          termSS_23(i,j,k) * S2mm3_filted_l(i,j,k) + &
                          termSS_31(i,j,k) * S3mm1_filted_l(i,j,k) + &
                          termSS_32(i,j,k) * S3mm2_filted_l(i,j,k) + &
                          termSS_33(i,j,k) * S3mm3_filted_l(i,j,k))
            Pi1S(m) = Pi1S(m) + vxr_D1S * dl_alpha(m,n)
            Pi1Sint = Pi1Sint + vxr_D1S * dl_alpha(m,n)
            !
            vxr_D1W = dreal(termSS_11(i,j,k) * W1mm1_filted_l(i,j,k) + &
                          termSS_12(i,j,k) * W1mm2_filted_l(i,j,k) + &
                          termSS_13(i,j,k) * W1mm3_filted_l(i,j,k) + &
                          termSS_21(i,j,k) * W2mm1_filted_l(i,j,k) + &
                          termSS_22(i,j,k) * W2mm2_filted_l(i,j,k) + &
                          termSS_23(i,j,k) * W2mm3_filted_l(i,j,k) + &
                          termSS_31(i,j,k) * W3mm1_filted_l(i,j,k) + &
                          termSS_32(i,j,k) * W3mm2_filted_l(i,j,k) + &
                          termSS_33(i,j,k) * W3mm3_filted_l(i,j,k))
            Pi1W(m) = Pi1W(m) + vxr_D1W * dl_alpha(m,n)
            Pi1Wint = Pi1Wint + vxr_D1W * dl_alpha(m,n)
            !
            vxr_D1D = dreal(termSS_11(i,j,k) * Amm11_filted_l(i,j,k) + &
                            termSS_12(i,j,k) * Amm12_filted_l(i,j,k) + &
                            termSS_13(i,j,k) * Amm13_filted_l(i,j,k) + &
                            termSS_21(i,j,k) * Amm21_filted_l(i,j,k) + &
                            termSS_22(i,j,k) * Amm22_filted_l(i,j,k) + &
                            termSS_23(i,j,k) * Amm23_filted_l(i,j,k) + &
                            termSS_31(i,j,k) * Amm31_filted_l(i,j,k) + &
                            termSS_32(i,j,k) * Amm32_filted_l(i,j,k) + &
                            termSS_33(i,j,k) * Amm33_filted_l(i,j,k))
            Pi1D(m) = Pi1D(m) + 1.d0/real(dim)*vxr_D1D * dl_alpha(m,n)
            Pi1Dint = Pi1Dint + 1.d0/real(dim)*vxr_D1D * dl_alpha(m,n)
            !
            vxr_D2S = dreal(termWW_11(i,j,k) * S1mm1_filted_l(i,j,k) + &
                    termWW_12(i,j,k) * S1mm2_filted_l(i,j,k) + &
                    termWW_13(i,j,k) * S1mm3_filted_l(i,j,k) + &
                    termWW_21(i,j,k) * S2mm1_filted_l(i,j,k) + &
                    termWW_22(i,j,k) * S2mm2_filted_l(i,j,k) + &
                    termWW_23(i,j,k) * S2mm3_filted_l(i,j,k) + &
                    termWW_31(i,j,k) * S3mm1_filted_l(i,j,k) + &
                    termWW_32(i,j,k) * S3mm2_filted_l(i,j,k) + &
                    termWW_33(i,j,k) * S3mm3_filted_l(i,j,k))
            Pi2S(m) = Pi2S(m) - vxr_D2S * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            Pi2Sint = Pi2Sint - vxr_D2S * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            !
            vxr_D2W = dreal(termWW_11(i,j,k) * W1mm1_filted_l(i,j,k) + &
                    termWW_12(i,j,k) * W1mm2_filted_l(i,j,k) + &
                    termWW_13(i,j,k) * W1mm3_filted_l(i,j,k) + &
                    termWW_21(i,j,k) * W2mm1_filted_l(i,j,k) + &
                    termWW_22(i,j,k) * W2mm2_filted_l(i,j,k) + &
                    termWW_23(i,j,k) * W2mm3_filted_l(i,j,k) + &
                    termWW_31(i,j,k) * W3mm1_filted_l(i,j,k) + &
                    termWW_32(i,j,k) * W3mm2_filted_l(i,j,k) + &
                    termWW_33(i,j,k) * W3mm3_filted_l(i,j,k))
            Pi2W(m) = Pi2W(m) - vxr_D2W * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            Pi2Wint = Pi2Wint - vxr_D2W * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            !
            vxr_D2D = dreal(termWW_11(i,j,k) * Amm11_filted_l(i,j,k) + &
                    termWW_12(i,j,k) * Amm12_filted_l(i,j,k) + &
                    termWW_13(i,j,k) * Amm13_filted_l(i,j,k) + &
                    termWW_21(i,j,k) * Amm21_filted_l(i,j,k) + &
                    termWW_22(i,j,k) * Amm22_filted_l(i,j,k) + &
                    termWW_23(i,j,k) * Amm23_filted_l(i,j,k) + &
                    termWW_31(i,j,k) * Amm31_filted_l(i,j,k) + &
                    termWW_32(i,j,k) * Amm32_filted_l(i,j,k) + &
                    termWW_33(i,j,k) * Amm33_filted_l(i,j,k))
            Pi2D(m) = Pi2D(m) - 1.d0/real(dim)*vxr_D2D * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            Pi2Dint = Pi2Dint - 1.d0/real(dim)*vxr_D2D * dl_alpha(m,n) ! Negative because of W convention, this is not a mistake
            !
            vxr_D3S = dreal(termWS_11(i,j,k) * S1mm1_filted_l(i,j,k) + &
                          termWS_12(i,j,k) * S1mm2_filted_l(i,j,k) + &
                          termWS_13(i,j,k) * S1mm3_filted_l(i,j,k) + &
                          termWS_21(i,j,k) * S2mm1_filted_l(i,j,k) + &
                          termWS_22(i,j,k) * S2mm2_filted_l(i,j,k) + &
                          termWS_23(i,j,k) * S2mm3_filted_l(i,j,k) + &
                          termWS_31(i,j,k) * S3mm1_filted_l(i,j,k) + &
                          termWS_32(i,j,k) * S3mm2_filted_l(i,j,k) + &
                          termWS_33(i,j,k) * S3mm3_filted_l(i,j,k))
            Pi3S(m) = Pi3S(m) + vxr_D3S * dl_alpha(m,n)
            Pi3Sint = Pi3Sint + vxr_D3S * dl_alpha(m,n)
            !
            vxr_D3W = dreal(termWS_11(i,j,k) * W1mm1_filted_l(i,j,k) + &
                          termWS_12(i,j,k) * W1mm2_filted_l(i,j,k) + &
                          termWS_13(i,j,k) * W1mm3_filted_l(i,j,k) + &
                          termWS_21(i,j,k) * W2mm1_filted_l(i,j,k) + &
                          termWS_22(i,j,k) * W2mm2_filted_l(i,j,k) + &
                          termWS_23(i,j,k) * W2mm3_filted_l(i,j,k) + &
                          termWS_31(i,j,k) * W3mm1_filted_l(i,j,k) + &
                          termWS_32(i,j,k) * W3mm2_filted_l(i,j,k) + &
                          termWS_33(i,j,k) * W3mm3_filted_l(i,j,k))
            Pi3W(m) = Pi3W(m) + vxr_D3W * dl_alpha(m,n)
            Pi3Wint = Pi3Wint + vxr_D3W * dl_alpha(m,n)
            !
            !
            vxr_D3D = dreal(termWS_11(i,j,k) * Amm11_filted_l(i,j,k) + &
                    termWS_12(i,j,k) * Amm12_filted_l(i,j,k) + &
                    termWS_13(i,j,k) * Amm13_filted_l(i,j,k) + &
                    termWS_21(i,j,k) * Amm21_filted_l(i,j,k) + &
                    termWS_22(i,j,k) * Amm22_filted_l(i,j,k) + &
                    termWS_23(i,j,k) * Amm23_filted_l(i,j,k) + &
                    termWS_31(i,j,k) * Amm31_filted_l(i,j,k) + &
                    termWS_32(i,j,k) * Amm32_filted_l(i,j,k) + &
                    termWS_33(i,j,k) * Amm33_filted_l(i,j,k))
            Pi3D(m) = Pi3D(m) + 1.d0/real(dim)*vxr_D3D * dl_alpha(m,n)
            Pi3Dint = Pi3Dint + 1.d0/real(dim)*vxr_D3D * dl_alpha(m,n)
            !
            vxr_D4S = dreal(termSD_11(i,j,k) * S1mm1_filted_l(i,j,k) + &
                          termSD_12(i,j,k) * S1mm2_filted_l(i,j,k) + &
                          termSD_13(i,j,k) * S1mm3_filted_l(i,j,k) + &
                          termSD_21(i,j,k) * S2mm1_filted_l(i,j,k) + &
                          termSD_22(i,j,k) * S2mm2_filted_l(i,j,k) + &
                          termSD_23(i,j,k) * S2mm3_filted_l(i,j,k) + &
                          termSD_31(i,j,k) * S3mm1_filted_l(i,j,k) + &
                          termSD_32(i,j,k) * S3mm2_filted_l(i,j,k) + &
                          termSD_33(i,j,k) * S3mm3_filted_l(i,j,k))
            Pi4S(m) = Pi4S(m) + 2.d0/real(dim) * vxr_D4S * dl_alpha(m,n)
            Pi4Sint = Pi4Sint + 2.d0/real(dim) * vxr_D4S * dl_alpha(m,n)
            !
            vxr_D4W = dreal(termSD_11(i,j,k) * W1mm1_filted_l(i,j,k) + &
                          termSD_12(i,j,k) * W1mm2_filted_l(i,j,k) + &
                          termSD_13(i,j,k) * W1mm3_filted_l(i,j,k) + &
                          termSD_21(i,j,k) * W2mm1_filted_l(i,j,k) + &
                          termSD_22(i,j,k) * W2mm2_filted_l(i,j,k) + &
                          termSD_23(i,j,k) * W2mm3_filted_l(i,j,k) + &
                          termSD_31(i,j,k) * W3mm1_filted_l(i,j,k) + &
                          termSD_32(i,j,k) * W3mm2_filted_l(i,j,k) + &
                          termSD_33(i,j,k) * W3mm3_filted_l(i,j,k))
            Pi4W(m) = Pi4W(m) + 2.d0/real(dim) * vxr_D4W * dl_alpha(m,n)
            Pi4Wint = Pi4Wint + 2.d0/real(dim) * vxr_D4W * dl_alpha(m,n)
            !
            !
            vxr_D4D = dreal(termSD_11(i,j,k) * Amm11_filted_l(i,j,k) + &
                    termSD_12(i,j,k) * Amm12_filted_l(i,j,k) + &
                    termSD_13(i,j,k) * Amm13_filted_l(i,j,k) + &
                    termSD_21(i,j,k) * Amm21_filted_l(i,j,k) + &
                    termSD_22(i,j,k) * Amm22_filted_l(i,j,k) + &
                    termSD_23(i,j,k) * Amm23_filted_l(i,j,k) + &
                    termSD_31(i,j,k) * Amm31_filted_l(i,j,k) + &
                    termSD_32(i,j,k) * Amm32_filted_l(i,j,k) + &
                    termSD_33(i,j,k) * Amm33_filted_l(i,j,k))
            Pi4D(m) = Pi4D(m) + 2.d0/real(dim)/real(dim)* vxr_D4D * dl_alpha(m,n)
            Pi4Dint = Pi4Dint + 2.d0/real(dim)/real(dim)* vxr_D4D * dl_alpha(m,n)
            !
            !
            vxr_D5S = dreal(termDD(i,j,k) * (S1mm1_filted_l(i,j,k) + &
                      S2mm2_filted_l(i,j,k) + S3mm3_filted_l(i,j,k)))
            Pi5S(m) = Pi5S(m) + 1.d0/real(dim)/real(dim)*vxr_D5S * dl_alpha(m,n)
            Pi5Sint = Pi5Sint + 1.d0/real(dim)/real(dim)*vxr_D5S * dl_alpha(m,n)
            !
            vxr_D5W = dreal(termDD(i,j,k) * (W1mm1_filted_l(i,j,k) + &
                      W2mm2_filted_l(i,j,k) + W3mm3_filted_l(i,j,k)))
            Pi5W(m) = Pi5W(m) + 1.d0/real(dim)/real(dim)*vxr_D5W * dl_alpha(m,n)
            Pi5Wint = Pi5Wint + 1.d0/real(dim)/real(dim)*vxr_D5W * dl_alpha(m,n)
            !
            !
            vxr_D5D = dreal(termDD(i,j,k) * (Amm11_filted_l(i,j,k) + &
                      Amm22_filted_l(i,j,k) + Amm33_filted_l(i,j,k)))
            Pi5D(m) = Pi5D(m) + 1.d0/real(dim)/real(dim)/real(dim)*vxr_D5D * dl_alpha(m,n)
            Pi5Dint = Pi5Dint + 1.d0/real(dim)/real(dim)/real(dim)*vxr_D5D * dl_alpha(m,n)
            !
          enddo
          enddo
          enddo
          !
          Pi1Sint = psum(Pi1Sint) / (ia*ja*ka)
          Pi2Sint = psum(Pi2Sint) / (ia*ja*ka)
          Pi3Sint = psum(Pi3Sint) / (ia*ja*ka)
          Pi4Sint = psum(Pi4Sint) / (ia*ja*ka)
          Pi5Sint = psum(Pi5Sint) / (ia*ja*ka)
          Pi1Wint = psum(Pi1Wint) / (ia*ja*ka)
          Pi2Wint = psum(Pi2Wint) / (ia*ja*ka)
          Pi3Wint = psum(Pi3Wint) / (ia*ja*ka)
          Pi4Wint = psum(Pi4Wint) / (ia*ja*ka)
          Pi5Wint = psum(Pi5Wint) / (ia*ja*ka)
          Pi1Dint = psum(Pi1Dint) / (ia*ja*ka)
          Pi2Dint = psum(Pi2Dint) / (ia*ja*ka)
          Pi3Dint = psum(Pi3Dint) / (ia*ja*ka)
          Pi4Dint = psum(Pi4Dint) / (ia*ja*ka)
          Pi5Dint = psum(Pi5Dint) / (ia*ja*ka)
          !
          if(mpirank==0) then
            call listwrite(hand_bS,l_sqrtalpha(m,n),Pi1Sint, Pi2Sint,Pi3Sint, &
              Pi4Sint, Pi5Sint)
            call listwrite(hand_bD,l_sqrtalpha(m,n),Pi1Dint, Pi2Dint,Pi3Dint, &
              Pi4Dint, Pi5Dint)
            call listwrite(hand_bW,l_sqrtalpha(m,n),Pi1Wint, Pi2Wint,Pi3Wint, &
              Pi4Wint, Pi5Wint)
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        !
        Pi1S(m) =  psum(Pi1S(m)) / (ia*ja*ka)
        Pi2S(m) =  psum(Pi2S(m)) / (ia*ja*ka)
        Pi3S(m) =  psum(Pi3S(m)) / (ia*ja*ka)
        Pi4S(m) =  psum(Pi4S(m)) / (ia*ja*ka)
        Pi5S(m) =  psum(Pi5S(m)) / (ia*ja*ka)
        Pi1W(m) =  psum(Pi1W(m)) / (ia*ja*ka)
        Pi2W(m) =  psum(Pi2W(m)) / (ia*ja*ka)
        Pi3W(m) =  psum(Pi3W(m)) / (ia*ja*ka)
        Pi4W(m) =  psum(Pi4W(m)) / (ia*ja*ka)
        Pi5W(m) =  psum(Pi5W(m)) / (ia*ja*ka)
        Pi1D(m) =  psum(Pi1D(m)) / (ia*ja*ka)
        Pi2D(m) =  psum(Pi2D(m)) / (ia*ja*ka)
        Pi3D(m) =  psum(Pi3D(m)) / (ia*ja*ka)
        Pi4D(m) =  psum(Pi4D(m)) / (ia*ja*ka)
        Pi5D(m) =  psum(Pi5D(m)) / (ia*ja*ka)
        !
        !
        !
        if(mpirank==0) then
            call listwrite(hand_bS, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
            call listwrite(hand_bS,Pi1S(m)+Pi2S(m)+Pi3S(m)+Pi4S(m)+Pi5S(m), & 
            Pi1S(m), Pi2S(m),Pi3S(m),Pi4S(m),Pi5S(m))
            call listwrite(hand_bD, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
            call listwrite(hand_bD,Pi1D(m)+Pi2D(m)+Pi3D(m)+Pi4D(m)+Pi5D(m), & 
            Pi1D(m), Pi2D(m),Pi3D(m),Pi4D(m),Pi5D(m))
            call listwrite(hand_bW, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
            call listwrite(hand_bW,Pi1W(m)+Pi2W(m)+Pi3W(m)+Pi4W(m)+Pi5W(m), & 
            Pi1W(m), Pi2W(m),Pi3W(m),Pi4W(m),Pi5W(m))
          !
          close(unit=hand_bS)
          close(unit=hand_bD)
          close(unit=hand_bW)
          !
          print *, '>>>>', outfilename2
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_ETS_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_ETS.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_aS, &
                      firstline='nstep time ell pi1S pi2S pi3S pi4S pi5S')
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_ETW_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_ETW.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_aW, &
                      firstline='nstep time ell pi1W pi2W pi3W pi4W pi5W')
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_ETD_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_ETD.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_aD, &
                      firstline='nstep time ell pi1D pi2D pi3D pi4D pi5D')
        !
        do m=1,num_l
          call listwrite(hand_aS,l_lim(m),Pi1S(m), Pi2S(m),&
            Pi3S(m), Pi4S(m), Pi5S(m))
          call listwrite(hand_aW,l_lim(m),Pi1W(m), Pi2W(m),&
            Pi3W(m), Pi4W(m), Pi5W(m))
          call listwrite(hand_aD,l_lim(m),Pi1D(m), Pi2D(m),&
            Pi3D(m), Pi4D(m), Pi5D(m))
        enddo
        !
        close(unit=hand_aS)
        close(unit=hand_aW)
        close(unit=hand_aD)
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call fftw_free(c_termSS_11)
      call fftw_free(c_termSS_12)
      call fftw_free(c_termSS_13)
      call fftw_free(c_termSS_21)
      call fftw_free(c_termSS_22)
      call fftw_free(c_termSS_23)
      call fftw_free(c_termSS_31)
      call fftw_free(c_termSS_32)
      call fftw_free(c_termSS_33)
      call fftw_free(c_termDD)
      call fftw_free(c_termSD_11)
      call fftw_free(c_termSD_12)
      call fftw_free(c_termSD_13)
      call fftw_free(c_termSD_21)
      call fftw_free(c_termSD_22)
      call fftw_free(c_termSD_23)
      call fftw_free(c_termSD_31)
      call fftw_free(c_termSD_32)
      call fftw_free(c_termSD_33)
      call fftw_free(c_termWW_11)
      call fftw_free(c_termWW_12)
      call fftw_free(c_termWW_13)
      call fftw_free(c_termWW_21)
      call fftw_free(c_termWW_22)
      call fftw_free(c_termWW_23)
      call fftw_free(c_termWW_31)
      call fftw_free(c_termWW_32)
      call fftw_free(c_termWW_33)
      call fftw_free(c_termWS_11)
      call fftw_free(c_termWS_12)
      call fftw_free(c_termWS_13)
      call fftw_free(c_termWS_21)
      call fftw_free(c_termWS_22)
      call fftw_free(c_termWS_23)
      call fftw_free(c_termWS_31)
      call fftw_free(c_termWS_32)
      call fftw_free(c_termWS_33)
      call fftw_free(c_S1mm1_filted_l)
      call fftw_free(c_S1mm2_filted_l)
      call fftw_free(c_S1mm3_filted_l)
      call fftw_free(c_S2mm1_filted_l)
      call fftw_free(c_S2mm2_filted_l)
      call fftw_free(c_S2mm3_filted_l)
      call fftw_free(c_S3mm1_filted_l)
      call fftw_free(c_S3mm2_filted_l)
      call fftw_free(c_S3mm3_filted_l)
      call fftw_free(c_W1mm1_filted_l)
      call fftw_free(c_W1mm2_filted_l)
      call fftw_free(c_W1mm3_filted_l)
      call fftw_free(c_W2mm1_filted_l)
      call fftw_free(c_W2mm2_filted_l)
      call fftw_free(c_W2mm3_filted_l)
      call fftw_free(c_W3mm1_filted_l)
      call fftw_free(c_W3mm2_filted_l)
      call fftw_free(c_W3mm3_filted_l)
      call mpistop
      deallocate(All_filted,S11_filted,S12_filted,S13_filted)
      deallocate(S21_filted,S22_filted,S23_filted)
      deallocate(S31_filted,S32_filted,S33_filted)
      deallocate(W12_filted,W21_filted,W13_filted,W31_filted,W23_filted,W32_filted)
      deallocate(k1,k2,k3)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
      deallocate(Pi1S,Pi2S,Pi3S,Pi4S,Pi5S,Pi1W,Pi2W,Pi3W,Pi4W,Pi5W,Pi1D,Pi2D,Pi3D,Pi4D,Pi5D)
      !
    end subroutine SGSET3D
    !
    subroutine SGSPi2Dint(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini, mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhocom
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11_filted,A12_filted,A21_filted,A22_filted
      complex(8), allocatable, dimension(:,:) :: All_filted,S11_filted,S12_filted,S21_filted,S22_filted,W12_filted,W21_filted
      complex(8), allocatable, dimension(:,:) :: All_filted_l,S11_filted_l,S12_filted_l,S21_filted_l,S22_filted_l
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: term1_11,term1_12,term1_21,term1_22,term2, &
                                                            term3_11,term3_12,term3_21,term3_22, &
                                                            term4_11,term4_22,term5,term6_11,term6_12, &
                                                            term6_21,term6_22,term7
      real(8) :: vxr_D1,vxr_D2,vxr_D3,vxr_D4,vxr_D5,vxr_D6,vxr_D7
      type(C_PTR) :: c_w1,c_w2,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A21_filted,c_A22_filted
      type(C_PTR) :: c_term1_11,c_term1_12,c_term1_21,c_term1_22,c_term2,c_term3_11,c_term3_12,c_term3_21,c_term3_22, &
                                                  c_term4_11,c_term4_22,c_term5,c_term6_11,c_term6_12,c_term6_21,c_term6_22,c_term7
      real(8), allocatable, dimension(:) :: Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhocom(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        rhocom(i,j)=rhocom(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!! allocation
      allocate(Pi1(1:num_l), Pi2(1:num_l), Pi3(1:num_l), Pi4(1:num_l), Pi5(1:num_l), Pi6(1:num_l), Pi7(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted, [imfftw,jmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted, [imfftw,jmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted, [imfftw,jmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted, [imfftw,jmfftw])
      !
      allocate(All_filted(1:im,1:jm),S11_filted(1:im,1:jm),S12_filted(1:im,1:jm),S21_filted(1:im,1:jm),S22_filted(1:im,1:jm), &
                W12_filted(1:im,1:jm),W21_filted(1:im,1:jm))
      !
      allocate(All_filted_l(1:im,1:jm),S11_filted_l(1:im,1:jm),S12_filted_l(1:im,1:jm),S21_filted_l(1:im,1:jm),  &
              S22_filted_l(1:im,1:jm))
      !
      c_term1_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_11, term1_11, [imfftw,jmfftw])
      c_term1_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_12, term1_12, [imfftw,jmfftw])
      c_term1_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_21, term1_21, [imfftw,jmfftw])
      c_term1_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term1_22, term1_22, [imfftw,jmfftw])
      c_term2    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term2, term2, [imfftw,jmfftw])
      c_term3_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_11, term3_11, [imfftw,jmfftw])
      c_term3_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_12, term3_12, [imfftw,jmfftw])
      c_term3_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_21, term3_21, [imfftw,jmfftw])
      c_term3_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term3_22, term3_22, [imfftw,jmfftw])
      c_term4_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_11, term4_11, [imfftw,jmfftw])
      c_term4_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term4_22, term4_22, [imfftw,jmfftw])
      c_term5    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term5, term5, [imfftw,jmfftw])
      c_term6_11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_11, term6_11, [imfftw,jmfftw])
      c_term6_12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_12, term6_12, [imfftw,jmfftw])
      c_term6_21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_21, term6_21, [imfftw,jmfftw])
      c_term6_22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term6_22, term6_22, [imfftw,jmfftw])
      c_term7    = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_term7, term7, [imfftw,jmfftw])
      !
      !
      Pi1 = 0.d0
      Pi2 =	0.d0
      Pi3 =	0.d0
      Pi4 =	0.d0
      Pi5 =	0.d0
      Pi6 =	0.d0
      Pi7 =	0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename2 = 'pp/SGS_Pi_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename2 = 'pp/SGS_Pi_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename2,handle=hand_b, &
                      firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        endif
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! Filtre scale :1
          w1_filted(i,j) = w1(i,j)     *Gl
          w2_filted(i,j) = w2(i,j)     *Gl
          rho_filted(i,j) = rhocom(i,j)*Gl
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do j=1,jm
        do i=1,im
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
          w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
          A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
          A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
          A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
          A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
          !
        end do
        end do
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        !
        do j=1,jm
        do i=1,im
          All_filted_l(i,j) = A11_filted(i,j)+A22_filted(i,j)
        !
          S11_filted_l(i,j) = A11_filted(i,j) - 1.d0/2.d0 * All_filted_l(i,j)
          S12_filted_l(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
          S21_filted_l(i,j) = S12_filted_l(i,j)
          S22_filted_l(i,j) = A22_filted(i,j) - 1.d0/2.d0 * All_filted_l(i,j)
        !
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !!!! Velocity Favre average and density average
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          do i=1,im
          do j=1,jm
            Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j)  = w1(i,j)    *Galpha
            w2_filted(i,j)  = w2(i,j)    *Galpha
            rho_filted(i,j) = rhocom(i,j)*Galpha
          enddo
          enddo
          !
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          do i=1,im
          do j=1,jm
            w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
            w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          enddo
          enddo
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
            w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
            A11_filted(i,j) = imag*w1_filted(i,j)*k1(i,j)
            A21_filted(i,j) = imag*w2_filted(i,j)*k1(i,j)
            A12_filted(i,j) = imag*w1_filted(i,j)*k2(i,j)
            A22_filted(i,j) = imag*w2_filted(i,j)*k2(i,j)
            !
          end do
          end do
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
          call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
          call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
          call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
          !
          !
          do j=1,jm
          do i=1,im
            All_filted(i,j) = A11_filted(i,j)+A22_filted(i,j)
            ! 
            S11_filted(i,j) = A11_filted(i,j) - 1.d0/2.d0 * All_filted(i,j)
            S12_filted(i,j) = (A12_filted(i,j) + A21_filted(i,j))*0.5d0
            S21_filted(i,j) = S12_filted(i,j)
            S22_filted(i,j) = A22_filted(i,j) - 1.d0/2.d0 * All_filted(i,j)
            !
            W12_filted(i,j) = (A12_filted(i,j)-A21_filted(i,j))*0.5d0
            W21_filted(i,j) = -1.d0*W12_filted(i,j)
          enddo
          enddo
          !
          !!!! Pi terms
          !
          do j=1,jm
          do i=1,im
            !term1_IJ = rho_filted*SI1_filted*SJ1_filted + rho_filted*SI2_filted*SJ2_filted
            term1_11(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j) + rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)
            term1_21(i,j) = rho_filted(i,j)*S21_filted(i,j)*S11_filted(i,j) + rho_filted(i,j)*S22_filted(i,j)*S12_filted(i,j)
            term1_12(i,j) = rho_filted(i,j)*S11_filted(i,j)*S21_filted(i,j) + rho_filted(i,j)*S12_filted(i,j)*S22_filted(i,j)
            term1_22(i,j) = rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j) + rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
            ! 
            !
            term2(i,j) = rho_filted(i,j)*S11_filted(i,j)*S11_filted(i,j)+rho_filted(i,j)*S12_filted(i,j)*S12_filted(i,j)+ &
                    rho_filted(i,j)*S21_filted(i,j)*S21_filted(i,j)+rho_filted(i,j)*S22_filted(i,j)*S22_filted(i,j)
            !
            ! term3_IJ = rho_filted*All_filted*SIJ_filted
            term3_11(i,j) = rho_filted(i,j)*All_filted(i,j)*S11_filted(i,j)
            term3_12(i,j) = rho_filted(i,j)*All_filted(i,j)*S12_filted(i,j)
            term3_21(i,j) = rho_filted(i,j)*All_filted(i,j)*S21_filted(i,j)
            term3_22(i,j) = rho_filted(i,j)*All_filted(i,j)*S22_filted(i,j)
            !
            !term4_IJ = rho_filted*WI1_filted*W1J_filted+rho_filted*WI2_filted*W2J_filted
            term4_11(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j)
            term4_22(i,j) = rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
            !
            term5(i,j) = rho_filted(i,j)*W12_filted(i,j)*W21_filted(i,j) + rho_filted(i,j)*W21_filted(i,j)*W12_filted(i,j)
            !
            !term6_IJ= rho_filted*(S1J_filted*WI1_filted-SI1_filted*W1J_filted) + rho_filted*(S2J_filted*WI2_filted-SI2_filted*W2J_filted)
            term6_11(i,j)=   rho_filted(i,j)*(S21_filted(i,j)*W12_filted(i,j)-S12_filted(i,j)*W21_filted(i,j))
            term6_12(i,j)= - rho_filted(i,j)*(S11_filted(i,j)*W12_filted(i,j)) + rho_filted(i,j)*(S22_filted(i,j)*W12_filted(i,j))
            term6_21(i,j)=   rho_filted(i,j)*(S11_filted(i,j)*W21_filted(i,j)) - rho_filted(i,j)*(S22_filted(i,j)*W21_filted(i,j))
            term6_22(i,j)=   rho_filted(i,j)*(S12_filted(i,j)*W21_filted(i,j)-S21_filted(i,j)*W12_filted(i,j))
            !
            term7(i,j) = rho_filted(i,j)*All_filted(i,j)*All_filted(i,j)
          enddo
          enddo
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,term1_11,term1_11)
          call fftw_mpi_execute_dft(forward_plan,term1_12,term1_12)
          call fftw_mpi_execute_dft(forward_plan,term1_21,term1_21)
          call fftw_mpi_execute_dft(forward_plan,term1_22,term1_22)
          call fftw_mpi_execute_dft(forward_plan,term2   ,term2   )
          call fftw_mpi_execute_dft(forward_plan,term3_11,term3_11)
          call fftw_mpi_execute_dft(forward_plan,term3_12,term3_12)
          call fftw_mpi_execute_dft(forward_plan,term3_21,term3_21)
          call fftw_mpi_execute_dft(forward_plan,term3_22,term3_22)
          call fftw_mpi_execute_dft(forward_plan,term4_11,term4_11)
          call fftw_mpi_execute_dft(forward_plan,term4_22,term4_22)
          call fftw_mpi_execute_dft(forward_plan,term5   ,term5   )
          call fftw_mpi_execute_dft(forward_plan,term6_11,term6_11)
          call fftw_mpi_execute_dft(forward_plan,term6_12,term6_12)
          call fftw_mpi_execute_dft(forward_plan,term6_21,term6_21)
          call fftw_mpi_execute_dft(forward_plan,term6_22,term6_22)
          call fftw_mpi_execute_dft(forward_plan,term7   ,term7   )
          !
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            term1_11(i,j) = term1_11(i,j)*Gphi/(1.d0*ia*ja)
            term1_12(i,j) = term1_12(i,j)*Gphi/(1.d0*ia*ja)
            term1_21(i,j) = term1_21(i,j)*Gphi/(1.d0*ia*ja)
            term1_22(i,j) = term1_22(i,j)*Gphi/(1.d0*ia*ja)
            term2(i,j)    = term2(i,j)   *Gphi/(1.d0*ia*ja)
            term3_11(i,j) = term3_11(i,j)*Gphi/(1.d0*ia*ja)
            term3_12(i,j) = term3_12(i,j)*Gphi/(1.d0*ia*ja)
            term3_21(i,j) = term3_21(i,j)*Gphi/(1.d0*ia*ja)
            term3_22(i,j) = term3_22(i,j)*Gphi/(1.d0*ia*ja)
            term4_11(i,j) = term4_11(i,j)*Gphi/(1.d0*ia*ja)
            term4_22(i,j) = term4_22(i,j)*Gphi/(1.d0*ia*ja)
            term5(i,j)    = term5(i,j)   *Gphi/(1.d0*ia*ja)
            term6_11(i,j) = term6_11(i,j)*Gphi/(1.d0*ia*ja)
            term6_12(i,j) = term6_12(i,j)*Gphi/(1.d0*ia*ja)
            term6_21(i,j) = term6_21(i,j)*Gphi/(1.d0*ia*ja)
            term6_22(i,j) = term6_22(i,j)*Gphi/(1.d0*ia*ja)
            term7(i,j)    = term7(i,j)   *Gphi/(1.d0*ia*ja)
            !
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,term1_11,term1_11)
          call fftw_mpi_execute_dft(backward_plan,term1_12,term1_12)
          call fftw_mpi_execute_dft(backward_plan,term1_21,term1_21)
          call fftw_mpi_execute_dft(backward_plan,term1_22,term1_22)
          call fftw_mpi_execute_dft(backward_plan,term2   ,term2   )
          call fftw_mpi_execute_dft(backward_plan,term3_11,term3_11)
          call fftw_mpi_execute_dft(backward_plan,term3_12,term3_12)
          call fftw_mpi_execute_dft(backward_plan,term3_21,term3_21)
          call fftw_mpi_execute_dft(backward_plan,term3_22,term3_22)
          call fftw_mpi_execute_dft(backward_plan,term4_11,term4_11)
          call fftw_mpi_execute_dft(backward_plan,term4_22,term4_22)
          call fftw_mpi_execute_dft(backward_plan,term5   ,term5   )
          call fftw_mpi_execute_dft(backward_plan,term6_11,term6_11)
          call fftw_mpi_execute_dft(backward_plan,term6_12,term6_12)
          call fftw_mpi_execute_dft(backward_plan,term6_21,term6_21)
          call fftw_mpi_execute_dft(backward_plan,term6_22,term6_22)
          call fftw_mpi_execute_dft(backward_plan,term7   ,term7   )
          !
          !
          do j=1,jm
          do i=1,im
            vxr_D1 = term1_11(i,j) * S11_filted_l(i,j) + term1_12(i,j) * S12_filted_l(i,j) + &
                    term1_21(i,j) * S21_filted_l(i,j) + term1_22(i,j) * S22_filted_l(i,j)
            Pi1(m) = Pi1(m) + vxr_D1 * dl_alpha(m,n)
            !
            vxr_D2 = term2(i,j) * All_filted_l(i,j)
            Pi2(m) = Pi2(m) + vxr_D2 * dl_alpha(m,n) * 0.5d0 
            !
            vxr_D3 = term3_11(i,j) * S11_filted_l(i,j) + term3_12(i,j) * S12_filted_l(i,j) + &
                    term3_21(i,j) * S21_filted_l(i,j) + term3_22(i,j) * S22_filted_l(i,j)
            Pi3(m) = Pi3(m) + vxr_D3 * dl_alpha(m,n)
            !
            vxr_D4 = term4_11(i,j) * S11_filted_l(i,j) + term4_22(i,j) * S22_filted_l(i,j)
            Pi4(m) = Pi4(m) - vxr_D4 * dl_alpha(m,n)
            !
            vxr_D5 = term5(i,j) * All_filted_l(i,j)
            Pi5(m) = Pi5(m) - vxr_D5 * dl_alpha(m,n) * 0.5d0
            !
            vxr_D6 = term6_11(i,j) * S11_filted_l(i,j) + term6_12(i,j) * S12_filted_l(i,j) + &
                    term6_21(i,j) * S21_filted_l(i,j) + term6_22(i,j) * S22_filted_l(i,j)
            Pi6(m) = Pi6(m) + vxr_D6 * dl_alpha(m,n)
            !
            vxr_D7 = term7(i,j) * All_filted_l(i,j)
            Pi7(m) = Pi7(m) + vxr_D7 * dl_alpha(m,n) * 0.25d0
          enddo
          enddo
          !
          vxr_D1 = psum(vxr_D1)
          vxr_D2 = psum(vxr_D2)
          vxr_D3 = psum(vxr_D3)
          vxr_D4 = psum(vxr_D4)
          vxr_D5 = psum(vxr_D5)
          vxr_D6 = psum(vxr_D6)
          vxr_D7 = psum(vxr_D7)
          !
          if(mpirank==0) then
            call listwrite(hand_b,l_sqrtalpha(m,n),vxr_D1 * dl_alpha(m,n), &
                          vxr_D2 * dl_alpha(m,n) / 2.d0, &
                          vxr_D3 * dl_alpha(m,n) , &
                          - vxr_D4 * dl_alpha(m,n), &
                          - vxr_D5 * dl_alpha(m,n) * 0.5d0, &
                          vxr_D6 * dl_alpha(m,n), &
                          vxr_D7 * dl_alpha(m,n) * 0.25d0)
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        !
        Pi1(m) =	 psum(Pi1(m)) / (ia*ja)
        Pi2(m) =	 psum(Pi2(m)) / (ia*ja)
        Pi3(m) =	 psum(Pi3(m)) / (ia*ja)
        Pi4(m) =	 psum(Pi4(m)) / (ia*ja)
        Pi5(m) =	 psum(Pi5(m)) / (ia*ja)
        Pi6(m) =	 psum(Pi6(m)) / (ia*ja)
        Pi7(m) =	 psum(Pi7(m)) / (ia*ja)
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pi_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pi.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi1(m), Pi2(m), &
                          Pi3(m), Pi4(m), Pi5(m),     &
                          Pi6(m), Pi7(m))
        end do
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_term1_11)
      call fftw_free(c_term1_12)
      call fftw_free(c_term1_21)
      call fftw_free(c_term1_22)
      call fftw_free(c_term2)
      call fftw_free(c_term3_11)
      call fftw_free(c_term3_12)
      call fftw_free(c_term3_21)
      call fftw_free(c_term3_22)
      call fftw_free(c_term4_11)
      call fftw_free(c_term4_22)
      call fftw_free(c_term5)
      call fftw_free(c_term6_11)
      call fftw_free(c_term6_12)
      call fftw_free(c_term6_21)
      call fftw_free(c_term6_22)
      call fftw_free(c_term7)
      call mpistop
      deallocate(All_filted,S11_filted,S12_filted,S21_filted,S22_filted,W12_filted,W21_filted)
      deallocate(All_filted_l,S11_filted_l,S12_filted_l,S21_filted_l,S22_filted_l)
      deallocate(k1,k2)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
      deallocate(Pi1,Pi2,Pi3,Pi4,Pi5,Pi6,Pi7)
      !
    end subroutine SGSPi2Dint
    !
    subroutine SGSPi3Dtot(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pi_tot
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1,w1w2,w1w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1,w2w2,w2w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1,w3w2,w3w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1_filted,w1w2_filted,w1w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1_filted,w2w2_filted,w2w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1_filted,w3w2_filted,w3w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      real(8), allocatable, dimension(:,:,:,:,:) :: tau ! 1:3,1:3,1:im,1:jm,1:km
      !
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w1w3
      type(C_PTR) :: c_w2w1,c_w2w2,c_w2w3
      type(C_PTR) :: c_w3w1,c_w3w2,c_w3w3
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w1w3_filted
      type(C_PTR) :: c_w2w1_filted,c_w2w2_filted,c_w2w3_filted
      type(C_PTR) :: c_w3w1_filted,c_w3w2_filted,c_w3w3_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw,kmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw,kmfftw])
      c_w1w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3, w1w3, [imfftw,jmfftw,kmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw,kmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw,kmfftw])
      c_w2w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3, w2w3, [imfftw,jmfftw,kmfftw])
      c_w3w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1, w3w1, [imfftw,jmfftw,kmfftw])
      c_w3w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2, w3w2, [imfftw,jmfftw,kmfftw])
      c_w3w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3, w3w3, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw,kmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw,kmfftw])
      c_w1w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3_filted, w1w3_filted, [imfftw,jmfftw,kmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw,kmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw,kmfftw])
      c_w2w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3_filted, w2w3_filted, [imfftw,jmfftw,kmfftw])
      c_w3w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1_filted, w3w1_filted, [imfftw,jmfftw,kmfftw])
      c_w3w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2_filted, w3w2_filted, [imfftw,jmfftw,kmfftw])
      c_w3w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3_filted, w3w3_filted, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        w1w1(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w2(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w3(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w1(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w2(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w3(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w1(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w2(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w3(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w1w3,w1w3)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,w2w3,w2w3)
      call fftw_mpi_execute_dft(forward_plan,w3w1,w3w1)
      call fftw_mpi_execute_dft(forward_plan,w3w2,w3w2)
      call fftw_mpi_execute_dft(forward_plan,w3w3,w3w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        w1w1(i,j,k)=w1w1(i,j,k)/(1.d0*ia*ja*ka)
        w1w2(i,j,k)=w1w2(i,j,k)/(1.d0*ia*ja*ka)
        w1w3(i,j,k)=w1w3(i,j,k)/(1.d0*ia*ja*ka)
        w2w1(i,j,k)=w2w1(i,j,k)/(1.d0*ia*ja*ka)
        w2w2(i,j,k)=w2w2(i,j,k)/(1.d0*ia*ja*ka)
        w2w3(i,j,k)=w2w3(i,j,k)/(1.d0*ia*ja*ka)
        w3w1(i,j,k)=w3w1(i,j,k)/(1.d0*ia*ja*ka)
        w3w2(i,j,k)=w3w2(i,j,k)/(1.d0*ia*ja*ka)
        w3w3(i,j,k)=w3w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      allocate(tau(1:3,1:3,1:im,1:jm,1:km))
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pi_tot(1:num_l))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      Pi_tot = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j,k)    = w1(i,j,k)    *Gl
          w2_filted(i,j,k)    = w2(i,j,k)    *Gl
          w3_filted(i,j,k)    = w3(i,j,k)    *Gl
          !
          w1w1_filted(i,j,k)  = w1w1(i,j,k) *Gl
          w1w2_filted(i,j,k)  = w1w2(i,j,k) *Gl
          w1w3_filted(i,j,k)  = w1w3(i,j,k) *Gl
          w2w1_filted(i,j,k)  = w2w1(i,j,k) *Gl
          w2w2_filted(i,j,k)  = w2w2(i,j,k) *Gl
          w2w3_filted(i,j,k)  = w2w3(i,j,k) *Gl
          w3w1_filted(i,j,k)  = w3w1(i,j,k) *Gl
          w3w2_filted(i,j,k)  = w3w2(i,j,k) *Gl
          w3w3_filted(i,j,k)  = w3w3(i,j,k) *Gl
          !
          rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w3_filted,w1w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w3_filted,w2w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w1_filted,w3w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w2_filted,w3w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w3_filted,w3w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
          tau(1,1,i,j,k) = dreal(w1w1_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w1_filted(i,j,k)) * dreal(w1_filted(i,j,k))
          tau(1,2,i,j,k) = dreal(w1w2_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w1_filted(i,j,k)) * dreal(w2_filted(i,j,k))
          tau(1,3,i,j,k) = dreal(w1w3_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w1_filted(i,j,k)) * dreal(w3_filted(i,j,k))
          tau(2,1,i,j,k) = dreal(w2w1_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w2_filted(i,j,k)) * dreal(w1_filted(i,j,k))
          tau(2,2,i,j,k) = dreal(w2w2_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w2_filted(i,j,k)) * dreal(w2_filted(i,j,k))
          tau(2,3,i,j,k) = dreal(w2w3_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w2_filted(i,j,k)) * dreal(w3_filted(i,j,k))
          tau(3,1,i,j,k) = dreal(w3w1_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w3_filted(i,j,k)) * dreal(w1_filted(i,j,k))
          tau(3,2,i,j,k) = dreal(w3w2_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w3_filted(i,j,k)) * dreal(w2_filted(i,j,k))
          tau(3,3,i,j,k) = dreal(w3w3_filted(i,j,k)) - dreal(rho_filted(i,j,k)) * dreal(w3_filted(i,j,k)) * dreal(w3_filted(i,j,k))
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          Pi_tot(m) = Pi_tot(m) + tau(1,1,i,j,k) * dreal(A11_filted(i,j,k)) + &
                                  tau(1,2,i,j,k) * dreal(A12_filted(i,j,k)) + &
                                  tau(1,3,i,j,k) * dreal(A13_filted(i,j,k)) + &
                                  tau(2,1,i,j,k) * dreal(A21_filted(i,j,k)) + &
                                  tau(2,2,i,j,k) * dreal(A22_filted(i,j,k)) + &
                                  tau(2,3,i,j,k) * dreal(A23_filted(i,j,k)) + &
                                  tau(3,1,i,j,k) * dreal(A31_filted(i,j,k)) + &
                                  tau(3,2,i,j,k) * dreal(A32_filted(i,j,k)) + &
                                  tau(3,3,i,j,k) * dreal(A33_filted(i,j,k))
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pi_tot(m) =	 psum(Pi_tot(m)) / (ia*ja*ka)
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pitot_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pitot.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pitot')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pi_tot(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w1w3)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w2w3)
      call fftw_free(c_w3w1)
      call fftw_free(c_w3w2)
      call fftw_free(c_w3w3)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w1w3_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_w2w3_filted)
      call fftw_free(c_w3w1_filted)
      call fftw_free(c_w3w2_filted)
      call fftw_free(c_w3w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call mpistop
      deallocate(k1,k2,k3,tau)
      deallocate(l_lim)
      deallocate(Pi_tot)
      !
    end subroutine SGSPi3Dtot
    !
    subroutine SGSPi3Dlocal(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: Pis1,Pis2,Pim2,Pim3,Pid
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1f,w2f,w3f,rhof
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
      !
      real(8), allocatable, dimension(:,:,:) :: All
      real(8), allocatable, dimension(:,:,:) :: S11,S12,S13
      real(8), allocatable, dimension(:,:,:) :: S21,S22,S23
      real(8), allocatable, dimension(:,:,:) :: S31,S32,S33
      real(8), allocatable, dimension(:,:,:) :: W12,W21
      real(8), allocatable, dimension(:,:,:) :: W13,W31
      real(8), allocatable, dimension(:,:,:) :: W23,W32
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1f,c_w2f,c_w3f,c_rho
      type(C_PTR) :: c_A11,c_A12,c_A13
      type(C_PTR) :: c_A21,c_A22,c_A23
      type(C_PTR) :: c_A31,c_A32,c_A33
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      allocate(All(1:im,1:jm,1:km),&
              S11(1:im,1:jm,1:km),S12(1:im,1:jm,1:km),S13(1:im,1:jm,1:km),&
              S21(1:im,1:jm,1:km),S22(1:im,1:jm,1:km),S23(1:im,1:jm,1:km),&
              S31(1:im,1:jm,1:km),S32(1:im,1:jm,1:km),S33(1:im,1:jm,1:km),&
              W12(1:im,1:jm,1:km),W21(1:im,1:jm,1:km),&
              W13(1:im,1:jm,1:km),W31(1:im,1:jm,1:km),&
              W23(1:im,1:jm,1:km),W32(1:im,1:jm,1:km))
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(Pis1(1:num_l),Pis2(1:num_l),Pim2(1:num_l),Pim3(1:num_l),Pid(1:num_l))
      !
      c_w1f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1f, w1f,  [imfftw,jmfftw,kmfftw])
      c_w2f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2f, w2f,  [imfftw,jmfftw,kmfftw])
      c_w3f = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3f, w3f,  [imfftw,jmfftw,kmfftw])
      c_rho = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho, rhof,[imfftw,jmfftw,kmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11,[imfftw,jmfftw,kmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12,[imfftw,jmfftw,kmfftw])
      c_A13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13, A13,[imfftw,jmfftw,kmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21,[imfftw,jmfftw,kmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22,[imfftw,jmfftw,kmfftw])
      c_A23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23, A23,[imfftw,jmfftw,kmfftw])
      c_A31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31, A31,[imfftw,jmfftw,kmfftw])
      c_A32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32, A32,[imfftw,jmfftw,kmfftw])
      c_A33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33, A33,[imfftw,jmfftw,kmfftw])
      !
      Pis1=0.d0
      Pis2=0.d0
      Pim2=0.d0
      Pim3=0.d0
      Pid=0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1 is (rho*u1) in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1f(i,j,k)    = w1(i,j,k)    *Gl
          w2f(i,j,k)    = w2(i,j,k)    *Gl
          w3f(i,j,k)    = w3(i,j,k)    *Gl
          !
          rhof(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1 is (rho*u1) in physical space
        call fftw_mpi_execute_dft(backward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(backward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(backward_plan,w3f,w3f)
        call fftw_mpi_execute_dft(backward_plan,rhof,rhof)
        !
        ! After this bloc, w1 is u1 in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1f(i,j,k) = w1f(i,j,k)/rhof(i,j,k)
          w2f(i,j,k) = w2f(i,j,k)/rhof(i,j,k)
          w3f(i,j,k) = w3f(i,j,k)/rhof(i,j,k)
          !
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1 is u1 in fourier space, A11 is A11 in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1f,w1f)
        call fftw_mpi_execute_dft(forward_plan,w2f,w2f)
        call fftw_mpi_execute_dft(forward_plan,w3f,w3f)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1f(i,j,k)  = w1f(i,j,k)/(1.d0*ia*ja*ka)
          w2f(i,j,k)  = w2f(i,j,k)/(1.d0*ia*ja*ka)
          w3f(i,j,k)  = w3f(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11(i,j,k) = imag*w1f(i,j,k)*k1(i,j,k)
          A21(i,j,k) = imag*w2f(i,j,k)*k1(i,j,k)
          A31(i,j,k) = imag*w3f(i,j,k)*k1(i,j,k)
          A12(i,j,k) = imag*w1f(i,j,k)*k2(i,j,k)
          A22(i,j,k) = imag*w2f(i,j,k)*k2(i,j,k)
          A32(i,j,k) = imag*w3f(i,j,k)*k2(i,j,k)
          A13(i,j,k) = imag*w1f(i,j,k)*k3(i,j,k)
          A23(i,j,k) = imag*w2f(i,j,k)*k3(i,j,k)
          A33(i,j,k) = imag*w3f(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11 is A11 in physical space
        call fftw_mpi_execute_dft(backward_plan,A11,A11)
        call fftw_mpi_execute_dft(backward_plan,A21,A21)
        call fftw_mpi_execute_dft(backward_plan,A31,A31)
        call fftw_mpi_execute_dft(backward_plan,A12,A12)
        call fftw_mpi_execute_dft(backward_plan,A22,A22)
        call fftw_mpi_execute_dft(backward_plan,A32,A32)
        call fftw_mpi_execute_dft(backward_plan,A13,A13)
        call fftw_mpi_execute_dft(backward_plan,A23,A23)
        call fftw_mpi_execute_dft(backward_plan,A33,A33)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          All(i,j,k) = dreal(A11(i,j,k)+A22(i,j,k)+A33(i,j,k))
          !
          S11(i,j,k) = dreal(A11(i,j,k)) - 1.d0/3.d0 * All(i,j,k)
          S22(i,j,k) = dreal(A22(i,j,k)) - 1.d0/3.d0 * All(i,j,k)
          S33(i,j,k) = dreal(A33(i,j,k)) - 1.d0/3.d0 * All(i,j,k)
          S12(i,j,k) = dreal(A12(i,j,k) + A21(i,j,k))*0.5d0
          S21(i,j,k) = S12(i,j,k)
          S13(i,j,k) = dreal(A13(i,j,k) + A31(i,j,k))*0.5d0
          S31(i,j,k) = S13(i,j,k)
          S23(i,j,k) = dreal(A23(i,j,k) + A32(i,j,k))*0.5d0
          S32(i,j,k) = S23(i,j,k)
          !
          W12(i,j,k) = dreal(A12(i,j,k)-A21(i,j,k))*0.5d0
          W21(i,j,k) = -1.d0*W12(i,j,k)
          W13(i,j,k) = dreal(A13(i,j,k)-A31(i,j,k))*0.5d0
          W31(i,j,k) = -1.d0*W13(i,j,k)
          W23(i,j,k) = dreal(A23(i,j,k)-A32(i,j,k))*0.5d0
          W32(i,j,k) = -1.d0*W23(i,j,k)
          !
        end do
        end do
        end do
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          Pis1(m) = Pis1(m) + dreal(rhof(i,j,k)) *(S11(i,j,k)*S11(i,j,k)*S11(i,j,k)+ &
                                                  S12(i,j,k)*S12(i,j,k)*S11(i,j,k)+ &
                                                  S13(i,j,k)*S13(i,j,k)*S11(i,j,k)+ &
                                                  S11(i,j,k)*S21(i,j,k)*S12(i,j,k)+ &
                                                  S12(i,j,k)*S22(i,j,k)*S12(i,j,k)+ &
                                                  S13(i,j,k)*S23(i,j,k)*S12(i,j,k)+ &
                                                  S11(i,j,k)*S31(i,j,k)*S13(i,j,k)+ &
                                                  S12(i,j,k)*S32(i,j,k)*S13(i,j,k)+ &
                                                  S13(i,j,k)*S33(i,j,k)*S13(i,j,k)+ &
                                                  S21(i,j,k)*S11(i,j,k)*S21(i,j,k)+ &
                                                  S22(i,j,k)*S12(i,j,k)*S21(i,j,k)+ &
                                                  S23(i,j,k)*S13(i,j,k)*S21(i,j,k)+ &
                                                  S21(i,j,k)*S21(i,j,k)*S22(i,j,k)+ &
                                                  S22(i,j,k)*S22(i,j,k)*S22(i,j,k)+ &
                                                  S23(i,j,k)*S23(i,j,k)*S22(i,j,k)+ &
                                                  S21(i,j,k)*S31(i,j,k)*S23(i,j,k)+ &
                                                  S22(i,j,k)*S32(i,j,k)*S23(i,j,k)+ &
                                                  S23(i,j,k)*S33(i,j,k)*S23(i,j,k)+ &
                                                  S31(i,j,k)*S11(i,j,k)*S31(i,j,k)+ &
                                                  S32(i,j,k)*S12(i,j,k)*S31(i,j,k)+ &
                                                  S33(i,j,k)*S13(i,j,k)*S31(i,j,k)+ &
                                                  S31(i,j,k)*S21(i,j,k)*S32(i,j,k)+ &
                                                  S32(i,j,k)*S22(i,j,k)*S32(i,j,k)+ &
                                                  S33(i,j,k)*S23(i,j,k)*S32(i,j,k)+ &
                                                  S31(i,j,k)*S31(i,j,k)*S33(i,j,k)+ &
                                                  S32(i,j,k)*S32(i,j,k)*S33(i,j,k)+ &
                                                  S33(i,j,k)*S33(i,j,k)*S33(i,j,k)) * &
                                                  l_lim(m) * l_lim(m)
          Pis2(m) = Pis2(m) + dreal(rhof(i,j,k)) *(W12(i,j,k)*W21(i,j,k)*S11(i,j,k)+ &
                                                  W13(i,j,k)*W31(i,j,k)*S11(i,j,k)+ &
                                                  W13(i,j,k)*W32(i,j,k)*S12(i,j,k)+ &
                                                  W12(i,j,k)*W23(i,j,k)*S13(i,j,k)+ &
                                                  W23(i,j,k)*W31(i,j,k)*S21(i,j,k)+ &
                                                  W21(i,j,k)*W12(i,j,k)*S22(i,j,k)+ &
                                                  W23(i,j,k)*W32(i,j,k)*S22(i,j,k)+ &
                                                  W21(i,j,k)*W13(i,j,k)*S23(i,j,k)+ &
                                                  W32(i,j,k)*W21(i,j,k)*S31(i,j,k)+ &
                                                  W31(i,j,k)*W12(i,j,k)*S32(i,j,k)+ &
                                                  W31(i,j,k)*W13(i,j,k)*S33(i,j,k)+ &
                                                  W32(i,j,k)*W23(i,j,k)*S33(i,j,k)) * &
                                                l_lim(m) * l_lim(m)
          Pim2(m) = Pim2(m) + dreal(rhof(i,j,k)) *(S11(i,j,k)*S11(i,j,k)*All(i,j,k)+ &
                                                  S12(i,j,k)*S12(i,j,k)*All(i,j,k)+ &
                                                  S13(i,j,k)*S13(i,j,k)*All(i,j,k)+ &
                                                  S21(i,j,k)*S21(i,j,k)*All(i,j,k)+ &
                                                  S22(i,j,k)*S22(i,j,k)*All(i,j,k)+ &
                                                  S23(i,j,k)*S23(i,j,k)*All(i,j,k)+ &
                                                  S31(i,j,k)*S31(i,j,k)*All(i,j,k)+ &
                                                  S32(i,j,k)*S32(i,j,k)*All(i,j,k)+ &
                                                  S33(i,j,k)*S33(i,j,k)*All(i,j,k)) * &
                                                l_lim(m) * l_lim(m)
          Pim3(m) = Pim3(m) + dreal(rhof(i,j,k)) *(W12(i,j,k)*W21(i,j,k)*All(i,j,k)+ &
                                                  W13(i,j,k)*W31(i,j,k)*All(i,j,k)+ &
                                                  W21(i,j,k)*W12(i,j,k)*All(i,j,k)+ &
                                                  W23(i,j,k)*W32(i,j,k)*All(i,j,k)+ &
                                                  W31(i,j,k)*W13(i,j,k)*All(i,j,k)+ &
                                                  W32(i,j,k)*W23(i,j,k)*All(i,j,k)) * &
                                              l_lim(m) * l_lim(m)
          Pid(m) = Pid(m) + dreal(rhof(i,j,k)) * All(i,j,k) * All(i,j,k) * All(i,j,k) * &
                                              l_lim(m) * l_lim(m)
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        Pis1(m) =	 psum(Pis1(m)) / (ia*ja*ka)
        Pis2(m) =	 - psum(Pis2(m)) / (ia*ja*ka)
        Pim2(m) =	 psum(Pim2(m)) / (ia*ja*ka) /3
        Pim3(m) =	 - psum(Pim3(m)) / (ia*ja*ka) /3
        Pid(m) =	 psum(Pid(m)) / (ia*ja*ka) /9
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pilocal_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pilocal.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pis1 pis2 pim2 pim3 pid')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), Pis1(m), Pis2(m), Pim2(m), Pim3(m), Pid(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1f)
      call fftw_free(c_w2f)
      call fftw_free(c_w3f)
      call fftw_free(c_rho)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A13)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_A23)
      call fftw_free(c_A31)
      call fftw_free(c_A32)
      call fftw_free(c_A33)
      call mpistop
      deallocate(k1,k2,k3)
      deallocate(l_lim)
      deallocate(All,S11,S12,S13,S21,S22,S23,S31,S32,S33)
      deallocate(W12,W13,W21,W23,W31,W32)
      deallocate(Pis1,Pis2,Pim2,Pim3,Pid)
      !
    end subroutine SGSPi3Dlocal
    !
    subroutine SGSLES3D(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl
      real(8), allocatable, dimension(:) :: AllAll,SijSij
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      real(8) :: All_filted, SijSij_this, Smodulu
      complex(8), allocatable, dimension(:,:,:) :: S11_filted,S12_filted,S13_filted
      complex(8), allocatable, dimension(:,:,:) :: S21_filted,S22_filted,S23_filted
      complex(8), allocatable, dimension(:,:,:) :: S31_filted,S32_filted,S33_filted
      !
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      !
      allocate(S11_filted(1:im,1:jm,1:km),S12_filted(1:im,1:jm,1:km),S13_filted(1:im,1:jm,1:km),&
                S21_filted(1:im,1:jm,1:km),S22_filted(1:im,1:jm,1:km),S23_filted(1:im,1:jm,1:km),&
                S31_filted(1:im,1:jm,1:km),S32_filted(1:im,1:jm,1:km),S33_filted(1:im,1:jm,1:km))
      !
      allocate(AllAll(1:num_l),SijSij(1:num_l))
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      AllAll = 0.d0
      SijSij = 0.d0
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        !
        !!!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0) ! Filtre scale :l
          !
          w1_filted(i,j,k)    = w1(i,j,k)    *Gl
          w2_filted(i,j,k)    = w2(i,j,k)    *Gl
          w3_filted(i,j,k)    = w3(i,j,k)    *Gl
          rho_filted(i,j,k)   = rhocom(i,j,k)*Gl
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do k=1,km
        do j=1,jm
        do i=1,im
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        !
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          All_filted = dreal(A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k))
          !
          S11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0 * All_filted
          S22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0 * All_filted
          S33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0 * All_filted
          S12_filted(i,j,k) = (A12_filted(i,j,k) + A21_filted(i,j,k))*0.5d0
          S21_filted(i,j,k) = S12_filted(i,j,k)
          S13_filted(i,j,k) = (A13_filted(i,j,k) + A31_filted(i,j,k))*0.5d0
          S31_filted(i,j,k) = S13_filted(i,j,k)
          S23_filted(i,j,k) = (A23_filted(i,j,k) + A32_filted(i,j,k))*0.5d0
          S32_filted(i,j,k) = S23_filted(i,j,k)
          !
          SijSij_this = dreal(S11_filted(i,j,k))**2 + dreal(S12_filted(i,j,k))**2 + dreal(S13_filted(i,j,k))**2 + &
                        dreal(S21_filted(i,j,k))**2 + dreal(S22_filted(i,j,k))**2 + dreal(S23_filted(i,j,k))**2 + &
                        dreal(S31_filted(i,j,k))**2 + dreal(S32_filted(i,j,k))**2 + dreal(S33_filted(i,j,k))**2
          Smodulu = dsqrt(2.d0*(SijSij_this+ 1.d0/3.d0 * All_filted**2))
          !
          AllAll(m) = AllAll(m) + dreal(rho_filted(i,j,k)) * All_filted * Smodulu
          SijSij(m) = SijSij(m) + dreal(rho_filted(i,j,k)) * SijSij_this
          !
        end do
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        AllAll(m) =	 psum(AllAll(m)) / (ia*ja*ka)
        SijSij(m) =	 psum(SijSij(m)) / (ia*ja*ka)
        !
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_LES_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_LES.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell AllAll SijSij')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), AllAll(m),SijSij(m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call mpistop
      deallocate(S11_filted,S12_filted,S13_filted,&
                  S21_filted,S22_filted,S23_filted,&
                  S31_filted,S32_filted,S33_filted)
      deallocate(k1,k2,k3)
      deallocate(l_lim)
      deallocate(AllAll,SijSij)
      !
    end subroutine SGSLES3D
    !
    subroutine SGSPi3Dint(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n,mmm
      character(len=128) :: infilename,outfilename,outfilename2
      character(len=4) :: stepname,mname
      
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3,ksq,Gl,Galpha,Gphi
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_b
      real(8) :: l_min, ratio_max, ratio_min
      real(8), allocatable, dimension(:,:) :: PiI
      real(8), allocatable, dimension(:) :: Pirank, Pisum
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:,:) :: A_filted
      real(8), allocatable, dimension(:,:,:) :: All_filted_l,All_filted
      real(8), allocatable, dimension(:,:,:,:,:) :: S_filted_l,S_filted,W_filted
      !
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:,:) :: term
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: termM
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,forward_plan,backward_plan!,forward_tensorplan,backward_tensorplan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_A_filted,c_term,c_termM
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      ! Fill spectral arrays with density-weighted velocity and density
      w1(1:im,1:jm,1:km) = CMPLX(vel(1:im,1:jm,1:km,1) * rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      w2(1:im,1:jm,1:km) = CMPLX(vel(1:im,1:jm,1:km,2) * rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      w3(1:im,1:jm,1:km) = CMPLX(vel(1:im,1:jm,1:km,3) * rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      rhocom(1:im,1:jm,1:km) = CMPLX(rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      deallocate(vel,rho)
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fft3d(w1,forward_plan)
      call fft3d(w2,forward_plan)
      call fft3d(w3,forward_plan)
      call fft3d(rhocom,forward_plan)
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km),ksq(1:im,1:jm,1:km))
      allocate(Gl(1:im,1:jm,1:km),Galpha(1:im,1:jm,1:km),Gphi(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      ksq = k1*k1 + k2*k2 + k3*k3
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(PiI(1:7,1:num_l), Pirank(1:7), Pisum(1:7))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A_filted = fftw_alloc_complex(9*alloc_local)
      call c_f_pointer(c_A_filted, A_filted,[imfftw,jmfftw,kmfftw,3_C_SIZE_T,3_C_SIZE_T])
      !
      allocate(All_filted_l(1:im,1:jm,1:km),S_filted_l(1:im,1:jm,1:km,1:3,1:3))
      !
      allocate(All_filted(1:im,1:jm,1:km),S_filted(1:im,1:jm,1:km,1:3,1:3),&
              W_filted(1:im,1:jm,1:km,1:3,1:3))
      !
      c_term = fftw_alloc_complex(9*alloc_local)
      call c_f_pointer(c_term, term, [imfftw,jmfftw,kmfftw,3_C_SIZE_T,3_C_SIZE_T])
      c_termM   = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termM,termM,[imfftw,jmfftw,kmfftw])
      !
      !
      Pirank = 0.d0
      PiI =	0.d0
      Pisum =	0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename2 = 'pp/SGS_Pi_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename2 = 'pp/SGS_Pi_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename2,handle=hand_b, &
                      firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        endif
        !
        !!!! l filetering --> outside
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        Gl = exp(-ksq*l_lim(m)**2/2.d0) ! Filtre scale :l
        w1_filted    = w1    *Gl
        w2_filted    = w2    *Gl
        w3_filted    = w3    *Gl
        rho_filted   = rhocom*Gl
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call ifft3d(w1_filted,backward_plan)
        call ifft3d(w2_filted,backward_plan)
        call ifft3d(w3_filted,backward_plan)
        call ifft3d(rho_filted,backward_plan)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        w1_filted = w1_filted/rho_filted
        w2_filted = w2_filted/rho_filted
        w3_filted = w3_filted/rho_filted
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fft3d(w1_filted,forward_plan)
        call fft3d(w2_filted,forward_plan)
        call fft3d(w3_filted,forward_plan)
        !
        call vector_gradient(A_filted, w1_filted, w2_filted, w3_filted, &
                                 k1, k2, k3)
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call ifft3dtensor(A_filted,backward_plan)
        !
        All_filted_l(:,:,:) = dreal(A_filted(:,:,:,1,1)+A_filted(:,:,:,2,2)+A_filted(:,:,:,3,3))
          !
        do j=1,3
        do i=1,3
          if(i==j)then
            S_filted_l(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j)) - 1.d0/3.d0 * All_filted_l(:,:,:)
          else
            S_filted_l(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j) + A_filted(:,:,:,j,i))*0.5d0
          endif
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          ! 
          !!!! Alpha filetering ---> inside
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          Galpha = exp(-ksq*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
          Gphi = exp(-ksq*l_phi(m,n)**2/2.d0)
          w1_filted  = w1    *Galpha
          w2_filted  = w2    *Galpha
          w3_filted  = w3    *Galpha
          rho_filted = rhocom*Galpha
          !
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call ifft3d(w1_filted,backward_plan)
          call ifft3d(w2_filted,backward_plan)
          call ifft3d(w3_filted,backward_plan)
          call ifft3d(rho_filted,backward_plan)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          w1_filted = w1_filted/rho_filted
          w2_filted = w2_filted/rho_filted
          w3_filted = w3_filted/rho_filted
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fft3d(w1_filted,forward_plan)
          call fft3d(w2_filted,forward_plan)
          call fft3d(w3_filted,forward_plan)
          !
          call vector_gradient(A_filted, w1_filted, w2_filted, w3_filted, &
                                 k1, k2, k3)
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call ifft3dtensor(A_filted,backward_plan)
          !
          !
          All_filted(:,:,:) = dreal(A_filted(:,:,:,1,1)+A_filted(:,:,:,2,2)+A_filted(:,:,:,3,3))
          do j=1,3
          do i=1,3
            if(i==j)then
              S_filted(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j)) - 1.d0/3.d0 * All_filted_l(:,:,:)
              W_filted(:,:,:,i,j)=0.d0
            else
              S_filted(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j) + A_filted(:,:,:,j,i))*0.5d0
              W_filted(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j) - A_filted(:,:,:,j,i))*0.5d0
            endif
          end do
          end do
          !
          !!!! Pi terms
          rho_filted = dreal(rho_filted)
          !
          !! Action I: SS
          call tensor_multi_3d_rhoABT(term,rho_filted,S_filted,S_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(1) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n) ! --> Pi1 = (SS)S
          Pirank(2) = sum(real(term(:,:,:,1,1),8)*All_filted_l)*dl_alpha(m,n)/3.d0 + &
                      sum(real(term(:,:,:,2,2),8)*All_filted_l)*dl_alpha(m,n)/3.d0 + &
                      sum(real(term(:,:,:,3,3),8)*All_filted_l)*dl_alpha(m,n)/3.d0
          !
          !! Action II: WW
          call tensor_multi_3d_rhoABT(term,rho_filted,W_filted,W_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(4) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n) ! --> Pi4 = (WW)S
          Pirank(5) = sum(real(term(:,:,:,1,1),8)*All_filted_l)*dl_alpha(m,n)*1.d0/3.d0 + &
                      sum(real(term(:,:,:,2,2),8)*All_filted_l)*dl_alpha(m,n)*1.d0/3.d0 + &
                      sum(real(term(:,:,:,3,3),8)*All_filted_l)*dl_alpha(m,n)*1.d0/3.d0 
                      ! --> Pi5 = (WW)Theta
          !
          !! Action III: SW
          call tensor_multi_3d_rhoABT(term,rho_filted,S_filted,W_filted,1)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(6) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n) !--> Pi6 (SW)S
          !
          !! Action IV: (STheta)
          call tensor_multi_3d_scalar(term,S_filted,All_filted)
          call tensor_multi_3d_scalar(term,rho_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(3) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)*2.d0/3.d0
          !
          !! Action V:(ThetaTheta)
          termM = rho_filted * All_filted*All_filted
          call fft3d(termM, forward_plan)
          termM = termM * Gphi
          call ifft3d(termM, backward_plan)
          Pirank(7) = sum(real(termM,8)*All_filted_l)*dl_alpha(m,n)/9.d0
          !
          !
          do i=1,7
            Pisum(i)=psum(Pirank(i))/(ia*ja*ka)
          enddo
          !
          PiI(:,m) = PiI(:,m) + Pisum(:)
          !
          if(mpirank==0) then
            call listwrite(hand_b,l_sqrtalpha(m,n),Pisum(1), Pisum(2),Pisum(3), &
                          Pisum(4), Pisum(5),Pisum(6), Pisum(7))
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        !
        if(mpirank==0) then
          call listwrite(hand_b,0.d0, 0.d0, 0.d0, &
                      0.d0, 0.d0, 0.d0,&
                      0.d0, 0.d0)
          call listwrite(hand_b,sum(PiI(:,m)), & 
          PiI(1,m), PiI(2,m),PiI(3,m),PiI(4,m), PiI(5,m),PiI(6,m), PiI(7,m))
          !
          close(unit=hand_b)
          !
          print *, '>>>>', outfilename2
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_Pi_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_Pi.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5 pi6 pi7')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), PiI(1,m), PiI(2,m),PiI(3,m),PiI(4,m), PiI(5,m),PiI(6,m), PiI(7,m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A_filted)
      call fftw_free(c_term)
      call mpistop
      deallocate(All_filted_l,S_filted_l)
      deallocate(All_filted,S_filted,W_filted)
      deallocate(k1,k2,k3,ksq,Galpha,Gl,Gphi)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha)
      deallocate(PiI,Pirank,Pisum)
      !
    end subroutine SGSPi3Dint
    !
    subroutine SGSPiB3Dint(thefilenumb)
      ! 
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n,mmm
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname,mname
      
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3,ksq,Gl,Galpha,Gphi
      real(8), allocatable, dimension(:,:,:,:) :: mag
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a,hand_pipI,hand_pipM,hand_pipA,hand_pipD
      real(8) :: l_min, ratio_max, ratio_min
      real(8), allocatable, dimension(:,:) :: PiI,PiM,PiA,PiD
      real(8), allocatable, dimension(:) :: Pirank, Pisum
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,b1,b2,b3,rhocom
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted,&
                                                              b1_filted,b2_filted,b3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:,:) :: A_filted, C_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:) :: H_filted
      real(8), allocatable, dimension(:,:,:) :: All_filted_l,All_filted
      real(8), allocatable, dimension(:,:,:,:,:) :: S_filted_l,S_filted,W_filted
      real(8), allocatable, dimension(:,:,:,:,:) :: C_filted_l,Sigma_filted_l,J_filted_l,Sigma_filted,J_filted
      !
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:,:,:) :: term
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: termM
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,c_b1,c_b2,c_b3
      type(C_PTR) :: forward_plan,backward_plan!,forward_tensorplan,backward_tensorplan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted,&
                     c_b1_filted,c_b2_filted,c_b3_filted
      type(C_PTR) :: c_A_filted,c_C_filted,c_H_filted,c_term,c_termM
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), mag(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='b1', var=mag(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='b2', var=mag(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='b3', var=mag(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      c_b1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_b1, b1, [imfftw,jmfftw,kmfftw])
      c_b2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_b2, b2, [imfftw,jmfftw,kmfftw])
      c_b3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_b3, b3, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      ! Fill spectral arrays with density-weighted velocity and density
      w1(1:im,1:jm,1:km) = CMPLX(vel(1:im,1:jm,1:km,1) * rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      w2(1:im,1:jm,1:km) = CMPLX(vel(1:im,1:jm,1:km,2) * rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      w3(1:im,1:jm,1:km) = CMPLX(vel(1:im,1:jm,1:km,3) * rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      rhocom(1:im,1:jm,1:km) = CMPLX(rho(1:im,1:jm,1:km), 0.d0, C_INTPTR_T)
      b1(1:im,1:jm,1:km) = CMPLX(mag(1:im,1:jm,1:km,1), 0.d0, C_INTPTR_T)
      b2(1:im,1:jm,1:km) = CMPLX(mag(1:im,1:jm,1:km,2), 0.d0, C_INTPTR_T)
      b3(1:im,1:jm,1:km) = CMPLX(mag(1:im,1:jm,1:km,3), 0.d0, C_INTPTR_T)
      deallocate(vel,mag,rho)
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fft3d(w1,forward_plan)
      call fft3d(w2,forward_plan)
      call fft3d(w3,forward_plan)
      call fft3d(rhocom,forward_plan)
      call fft3d(b1,forward_plan)
      call fft3d(b2,forward_plan)
      call fft3d(b3,forward_plan)
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km),ksq(1:im,1:jm,1:km))
      allocate(Gl(1:im,1:jm,1:km),Galpha(1:im,1:jm,1:km),Gphi(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      ksq = k1*k1 + k2*k2 + k3*k3
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      !!!!
      allocate(PiI(1:7,1:num_l),PiM(1:5,1:num_l),PiA(1:13,1:num_l),PiD(1:13,1:num_l), Pirank(1:13), Pisum(1:13))
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      c_b1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_b1_filted, b1_filted,  [imfftw,jmfftw,kmfftw])
      c_b2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_b2_filted, b2_filted,  [imfftw,jmfftw,kmfftw])
      c_b3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_b3_filted, b3_filted,  [imfftw,jmfftw,kmfftw])
      !
      c_A_filted = fftw_alloc_complex(9*alloc_local)
      call c_f_pointer(c_A_filted, A_filted,[imfftw,jmfftw,kmfftw,3_C_SIZE_T,3_C_SIZE_T])
      c_C_filted = fftw_alloc_complex(9*alloc_local)
      call c_f_pointer(c_C_filted, C_filted,[imfftw,jmfftw,kmfftw,3_C_SIZE_T,3_C_SIZE_T])
      c_H_filted = fftw_alloc_complex(3*alloc_local)
      call c_f_pointer(c_H_filted, H_filted,[imfftw,jmfftw,kmfftw,3_C_SIZE_T])
      !
      allocate(All_filted_l(1:im,1:jm,1:km),S_filted_l(1:im,1:jm,1:km,1:3,1:3))
      allocate(C_filted_l(1:im,1:jm,1:km,1:3,1:3),Sigma_filted_l(1:im,1:jm,1:km,1:3,1:3),&
               J_filted_l(1:im,1:jm,1:km,1:3,1:3))
      !
      allocate(All_filted(1:im,1:jm,1:km),S_filted(1:im,1:jm,1:km,1:3,1:3),&
              W_filted(1:im,1:jm,1:km,1:3,1:3))
      allocate(Sigma_filted(1:im,1:jm,1:km,1:3,1:3),J_filted(1:im,1:jm,1:km,1:3,1:3))
      !
      c_term = fftw_alloc_complex(9*alloc_local)
      call c_f_pointer(c_term, term, [imfftw,jmfftw,kmfftw,3_C_SIZE_T,3_C_SIZE_T])
      c_termM   = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_termM,termM,[imfftw,jmfftw,kmfftw])
      !
      !
      PiI =	0.d0
      PiM = 0.d0
      PiA = 0.d0
      PiD = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do m=1,num_l
        !
        !!!!!! Filter to get Sij filted by l
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        if(mpirank == 0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_PiI_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_PiI_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename,handle=hand_pipI, &
                      firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7')
          !
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_PiM_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_PiM_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename,handle=hand_pipM, &
                      firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5')
          !
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_PiA_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_PiA_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename,handle=hand_pipA, &
          firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7 pi8 pi9 pi10 pi11 pi12 pi13')
          !
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_PiD_precise_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_PiD_precise_'//mname//'.dat'
          endif
          call listinit(filename=outfilename,handle=hand_pipD, &
          firstline='nstep time sqrtalpha pi1 pi2 pi3 pi4 pi5 pi6 pi7 pi8 pi9 pi10 pi11 pi12 pi13')
          !
        endif
        !
        !!!! l filetering --> outside
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        Gl = exp(-ksq*l_lim(m)**2/2.d0) ! Filtre scale :l
        w1_filted    = w1    *Gl
        w2_filted    = w2    *Gl
        w3_filted    = w3    *Gl
        rho_filted   = rhocom*Gl
        b1_filted    = b1    *Gl
        b2_filted    = b2    *Gl
        b3_filted    = b3    *Gl
        !
        ! Only velocity do Favre filtering
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call ifft3d(w1_filted,backward_plan)
        call ifft3d(w2_filted,backward_plan)
        call ifft3d(w3_filted,backward_plan)
        call ifft3d(rho_filted,backward_plan)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        w1_filted = w1_filted/rho_filted
        w2_filted = w2_filted/rho_filted
        w3_filted = w3_filted/rho_filted
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fft3d(w1_filted,forward_plan)
        call fft3d(w2_filted,forward_plan)
        call fft3d(w3_filted,forward_plan)
        !
        call vector_gradient(A_filted, w1_filted, w2_filted, w3_filted, &
                                 k1, k2, k3)
        call vector_gradient(C_filted, b1_filted, b2_filted, b3_filted, &
                                 k1, k2, k3)
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call ifft3dtensor(A_filted,backward_plan)
        call ifft3dtensor(C_filted,backward_plan)
        !
        All_filted_l(:,:,:) = dreal(A_filted(:,:,:,1,1)+A_filted(:,:,:,2,2)+A_filted(:,:,:,3,3))
        !
        do j=1,3
        do i=1,3
          if(i==j)then
            S_filted_l(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j)) - 1.d0/3.d0 * All_filted_l(:,:,:)
          else
            S_filted_l(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j) + A_filted(:,:,:,j,i))*0.5d0
          endif
          Sigma_filted_l(:,:,:,i,j)=dreal(C_filted(:,:,:,i,j) + C_filted(:,:,:,j,i))*0.5d0
          J_filted_l(:,:,:,i,j)=dreal(C_filted(:,:,:,i,j) - C_filted(:,:,:,j,i))*0.5d0
          C_filted_l(:,:,:,i,j)=C_filted(:,:,:,i,j)
        end do
        end do
        !
        if(mpirank==0)  print *, '** l filted!'
        !
        !!!!!! Begin integral
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          ! 
          !!!! Alpha filetering ---> inside
          ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
          Galpha = exp(-ksq*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
          Gphi = exp(-ksq*l_phi(m,n)**2/2.d0)
          !
          w1_filted  = w1    *Galpha
          w2_filted  = w2    *Galpha
          w3_filted  = w3    *Galpha
          rho_filted = rhocom*Galpha
          b1_filted  = b1    *Galpha
          b2_filted  = b2    *Galpha
          b3_filted  = b3    *Galpha
          !
          ! Only velocity do Favre filtering
          ! After this bloc, w1_filted is (rho*u1)_filted in physical space
          call ifft3d(w1_filted,backward_plan)
          call ifft3d(w2_filted,backward_plan)
          call ifft3d(w3_filted,backward_plan)
          call ifft3d(rho_filted,backward_plan)
          ! B need to be in physical space
          call ifft3d(b1_filted,backward_plan)
          call ifft3d(b2_filted,backward_plan)
          call ifft3d(b3_filted,backward_plan)
          !
          ! After this bloc, w1_filted is u1_filted in physical space
          w1_filted = w1_filted/rho_filted
          w2_filted = w2_filted/rho_filted
          w3_filted = w3_filted/rho_filted
          !
          ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
          call fft3d(w1_filted,forward_plan)
          call fft3d(w2_filted,forward_plan)
          call fft3d(w3_filted,forward_plan)
          call fft3d(rho_filted,forward_plan)
          !
          call vector_gradient(A_filted, w1_filted, w2_filted, w3_filted, &
                                 k1, k2, k3)
          call vector_gradient(C_filted, b1_filted, b2_filted, b3_filted, &
                                 k1, k2, k3)
          call scalar_gradient(H_filted, rho_filted, k1,k2,k3)
          !
          ! After this bloc, A11_filted is A11_filted in physical space
          call ifft3dtensor(A_filted,backward_plan)
          call ifft3dtensor(C_filted,backward_plan)
          call ifft3dvector(H_filted,backward_plan)
          call ifft3d(rho_filted,backward_plan)
          !
          All_filted(:,:,:) = dreal(A_filted(:,:,:,1,1)+A_filted(:,:,:,2,2)+A_filted(:,:,:,3,3))
          do j=1,3
            do i=1,3
              if(i==j)then
                S_filted(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j)) - 1.d0/3.d0 * All_filted_l(:,:,:)
                W_filted(:,:,:,i,j)=0.d0
                J_filted(:,:,:,i,j)=0.d0
              else
                S_filted(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j) + A_filted(:,:,:,j,i))*0.5d0
                W_filted(:,:,:,i,j)=dreal(A_filted(:,:,:,i,j) - A_filted(:,:,:,j,i))*0.5d0
                J_filted(:,:,:,i,j)=dreal(C_filted(:,:,:,i,j) - C_filted(:,:,:,j,i))*0.5d0
              endif
              Sigma_filted(:,:,:,i,j)=dreal(C_filted(:,:,:,i,j) + C_filted(:,:,:,j,i))*0.5d0
            end do
            H_filted(:,:,:,j)=H_filted(:,:,:,j)/rho_filted
          end do
          !
          !!!! Pi terms
          rho_filted = dreal(rho_filted)
          !
          !!! Mechanism 1:Inertial - I
          Pirank = 0.d0
          Pisum =	0.d0
          !! Action I: SS
          call tensor_multi_3d_rhoABT(term,rho_filted,S_filted,S_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(1) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)
          Pirank(2) = sum(real(term(:,:,:,1,1),8)*All_filted_l)*dl_alpha(m,n)+ &
                      sum(real(term(:,:,:,2,2),8)*All_filted_l)*dl_alpha(m,n)+ &
                      sum(real(term(:,:,:,3,3),8)*All_filted_l)*dl_alpha(m,n)
          Pirank(2) = Pirank(2)/3.d0
          !
          !! Action II: WW
          call tensor_multi_3d_rhoABT(term,rho_filted,W_filted,W_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(3) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)
          Pirank(4) = sum(real(term(:,:,:,1,1),8)*All_filted_l)*dl_alpha(m,n) + &
                      sum(real(term(:,:,:,2,2),8)*All_filted_l)*dl_alpha(m,n) + &
                      sum(real(term(:,:,:,3,3),8)*All_filted_l)*dl_alpha(m,n) 
          Pirank(4) = Pirank(4)/3.d0
          !
          !! Action III: SW
          call tensor_multi_3d_rhoABT(term,rho_filted,S_filted,W_filted,1)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(5) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)
          !
          !! Action IV: (STheta)
          call tensor_multi_3d_scalar(term,S_filted,All_filted)
          call tensor_multi_3d_scalar(term,rho_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(6) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)*2.d0/3.d0
          !
          !! Action V:(ThetaTheta)
          termM = rho_filted * All_filted*All_filted
          call fft3d(termM, forward_plan)
          termM = termM * Gphi
          call ifft3d(termM, backward_plan)
          Pirank(7) = sum(real(termM,8)*All_filted_l)*dl_alpha(m,n)/9.d0
          !
          !! Summation
          do i=1,7
            Pisum(i)=psum(Pirank(i))/(ia*ja*ka)
          enddo
          !
          PiI(:,m) = PiI(:,m) + Pisum(1:7)
          !
          if(mpirank==0) then
            call listwrite(hand_pipI,l_sqrtalpha(m,n),Pisum(1), Pisum(2),Pisum(3), &
                          Pisum(4), Pisum(5),Pisum(6), Pisum(7))
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
          !!! Mechanism B:Velocities
          Pirank = 0.d0
          Pisum =	0.d0
          !
          !! Action I: SigmaSigma
          call tensor_multi_3d_ABT(term,Sigma_filted,Sigma_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(1) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)
          Pirank(2) = sum(real(term(:,:,:,1,1),8)*All_filted_l)*dl_alpha(m,n) + &
                      sum(real(term(:,:,:,2,2),8)*All_filted_l)*dl_alpha(m,n) + &
                      sum(real(term(:,:,:,3,3),8)*All_filted_l)*dl_alpha(m,n)
          Pirank(2) = -Pirank(2)/6.d0
          !
          !! Action II: JJ
          call tensor_multi_3d_ABT(term,J_filted,J_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(3) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)
          Pirank(4) = sum(real(term(:,:,:,1,1),8)*All_filted_l)*dl_alpha(m,n) + &
                      sum(real(term(:,:,:,2,2),8)*All_filted_l)*dl_alpha(m,n) + &
                      sum(real(term(:,:,:,3,3),8)*All_filted_l)*dl_alpha(m,n)
          Pirank(4) = -Pirank(2)/6.d0
          !
          !! Action III: SigmaJ
          call tensor_multi_3d_ABT(term,J_filted,Sigma_filted,sym=.true.)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,  backward_plan)
          Pirank(5) = sum(real(term,8)*S_filted_l)*dl_alpha(m,n)
          !
          !! Summation
          do i=1,5
            Pisum(i)=psum(Pirank(i))/(ia*ja*ka)
          enddo
          !
          PiM(:,m) = PiM(:,m) + Pisum(1:5)
          !
          if(mpirank==0) then
            call listwrite(hand_pipM,l_sqrtalpha(m,n),Pisum(1), Pisum(2),Pisum(3), &
                          Pisum(4), Pisum(5))
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
          !!! Mechanism C:Advections
          Pirank = 0.d0
          Pisum =	0.d0
          !
          !! Action I: Sigma S
          call tensor_multi_3d_ABT(term,Sigma_filted,S_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(1) = - sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(2) = - sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action II: Sigma W
          call tensor_multi_3d_ABT(term,Sigma_filted,W_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(3) = - sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(4) = - sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action III: Sigma Theta
          call tensor_multi_3d_scalar(term,Sigma_filted,All_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(5) = - sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)/3.d0
          Pirank(6) = - sum(real(term,8)*J_filted_l)*dl_alpha(m,n)/3.d0
          !
          !! Action IIII: J S
          call tensor_multi_3d_ABT(term,J_filted,S_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(7) = - sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(8) = - sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action V: J W
          call tensor_multi_3d_ABT(term,J_filted,W_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(9) = - sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(10) = - sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action VI: J Theta
          call tensor_multi_3d_scalar(term,J_filted,All_filted)
          call fft3dtensor(term,forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(11) = - sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)/3.d0
          Pirank(12) = - sum(real(term,8)*J_filted_l)*dl_alpha(m,n)/3.d0
          !
          !! Action VII: density
          !
          term = 0.d0
          do j=1,3
          do k=1,3
            term(:,:,:,1,j)= term(:,:,:,1,j) + dreal(b1_filted) * dreal(H_filted(:,:,:,k)) * A_filted(:,:,:,j,k)
            term(:,:,:,2,j)= term(:,:,:,2,j) + dreal(b2_filted) * dreal(H_filted(:,:,:,k)) * A_filted(:,:,:,j,k)
            term(:,:,:,3,j)= term(:,:,:,3,j) + dreal(b3_filted) * dreal(H_filted(:,:,:,k)) * A_filted(:,:,:,j,k)
          enddo
          enddo
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(13) = sum(real(term,8)*real(C_filted_l,8))*dl_alpha(m,n)
          !
          !! Summation
          do i=1,13
            Pisum(i)=psum(Pirank(i))/(ia*ja*ka)
          enddo
          !
          PiA(:,m) = PiA(:,m) + Pisum(:)
          !
          if(mpirank==0) then
            call listwrite(hand_pipA,l_sqrtalpha(m,n),Pisum(1), Pisum(2),Pisum(3), &
                          Pisum(4), Pisum(5), Pisum(6), Pisum(7), Pisum(8), &
                          Pisum(9), Pisum(10), Pisum(11), Pisum(12), Pisum(13))
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
          !!! Mechanism D:Dynamo
          Pirank = 0.d0
          Pisum =	0.d0
          !
          !! Action I: S Sigma 
          call tensor_multi_3d_ABT(term,S_filted,Sigma_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(1) = sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(2) = sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action II: W Sigma 
          call tensor_multi_3d_ABT(term,W_filted,Sigma_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(3) = sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(4) = sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action III: Theta Sigma 
          call tensor_multi_3d_scalar(term,Sigma_filted,All_filted,rev=.true.)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(5) = sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)/3.d0
          Pirank(6) = sum(real(term,8)*J_filted_l)*dl_alpha(m,n)/3.d0
          !
          !! Action IIII: S J 
          call tensor_multi_3d_ABT(term,S_filted,J_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(7) = sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(8) = sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action V: W J 
          call tensor_multi_3d_ABT(term,W_filted,J_filted)
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(9) = sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)
          Pirank(10) =sum(real(term,8)*J_filted_l)*dl_alpha(m,n)
          !
          !! Action VI: Theta J 
          call tensor_multi_3d_scalar(term,J_filted,All_filted,rev=.true.)
          call fft3dtensor(term,forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(11) = sum(real(term,8)*Sigma_filted_l)*dl_alpha(m,n)/3.d0
          Pirank(12) = sum(real(term,8)*J_filted_l)*dl_alpha(m,n)/3.d0
          !
          !! Action VII: density
          term = 0.d0
          do i=1,3
          do k=1,3
            term(:,:,:,i,1)= term(:,:,:,i,1) + dreal(b1_filted) * dreal(H_filted(:,:,:,k)) * A_filted(:,:,:,i,k)
            term(:,:,:,i,2)= term(:,:,:,i,2) + dreal(b2_filted) * dreal(H_filted(:,:,:,k)) * A_filted(:,:,:,i,k)
            term(:,:,:,i,3)= term(:,:,:,i,3) + dreal(b3_filted) * dreal(H_filted(:,:,:,k)) * A_filted(:,:,:,i,k)
          enddo
          enddo
          call fft3dtensor(term, forward_plan)
          call tensor_scale_3d(term,Gphi)
          call ifft3dtensor(term,backward_plan)
          Pirank(13) = - sum(real(term,8)*real(C_filted_l,8))*dl_alpha(m,n)
          !
          !! Summation
          do i=1,13
            Pisum(i)=psum(Pirank(i))/(ia*ja*ka)
          enddo
          !
          PiD(:,m) = PiD(:,m) + Pisum(:)
          !
          if(mpirank==0) then
            call listwrite(hand_pipD,l_sqrtalpha(m,n),Pisum(1), Pisum(2),Pisum(3), &
                          Pisum(4), Pisum(5), Pisum(6), Pisum(7), Pisum(8), &
                          Pisum(9), Pisum(10), Pisum(11), Pisum(12), Pisum(13))
          endif
          !
          call mpi_barrier(mpi_comm_world,ierr)
          !
        enddo
        ! 
        if(mpirank==0) then
          call listwrite(hand_pipI,0.d0, 0.d0, 0.d0, &
                      0.d0, 0.d0, 0.d0,&
                      0.d0, 0.d0)
          call listwrite(hand_pipI,sum(PiI(:,m)), & 
          PiI(1,m), PiI(2,m),PiI(3,m),PiI(4,m), PiI(5,m),PiI(6,m), PiI(7,m))
          !
          call listwrite(hand_pipM, 0.d0, 0.d0, 0.d0, &
                      0.d0, 0.d0, 0.d0)
          call listwrite(hand_pipM,sum(PiM(:,m)), & 
          PiM(1,m), PiM(2,m),PiM(3,m),PiM(4,m), PiM(5,m))
          !
          call listwrite(hand_pipA, 0.d0, &
          0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
          0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
          call listwrite(hand_pipA,sum(PiA(:,m)), & 
          PiA(1,m), PiA(2,m),PiA(3,m),PiA(4,m), PiA(5,m), PiA(6,m), PiA(7,m),&
          PiA(8,m), PiA(9,m),PiA(10,m),PiA(11,m), PiA(12,m), PiA(13,m))
          !
          call listwrite(hand_pipD, 0.d0, &
          0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,&
          0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0)
          call listwrite(hand_pipD,sum(PiD(:,m)), & 
          PiD(1,m), PiD(2,m),PiD(3,m),PiD(4,m), PiD(5,m), PiD(6,m), PiD(7,m),&
          PiD(8,m), PiD(9,m),PiD(10,m),PiD(11,m), PiD(12,m), PiD(13,m))
          !
          close(unit=hand_pipI)
          close(unit=hand_pipM)
          close(unit=hand_pipA)
          close(unit=hand_pipD)
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
      enddo
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_PiI_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_PiI.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pis1 pim2 pis2 pim3 pis3 pim1 pid')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), PiI(1,m), PiI(2,m),PiI(3,m),PiI(4,m), PiI(5,m),PiI(6,m), PiI(7,m))
        enddo
        !
        print *, '>>>>', outfilename
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_PiM_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_PiM.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), PiM(1,m), PiM(2,m),PiM(3,m),PiM(4,m), PiM(5,m))
        enddo
        !
        print *, '>>>>', outfilename
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_PiA_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_PiA.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5 pi6 pi7 pi8 pi9 pi10 pi11 pi12 pi13')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), PiA(1,m), PiA(2,m),PiA(3,m),PiA(4,m), PiA(5,m),PiA(6,m), PiA(7,m),&
                                          PiA(8,m), PiA(9,m),PiA(10,m),PiA(11,m), PiA(12,m),PiA(13,m))
        enddo
        !
        print *, '>>>>', outfilename
        !
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_PiD_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_PiD.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time ell pi1 pi2 pi3 pi4 pi5 pi6 pi7 pi8 pi9 pi10 pi11 pi12 pi13')
        do m=1,num_l
          call listwrite(hand_a,l_lim(m), PiD(1,m), PiD(2,m),PiD(3,m),PiD(4,m), PiD(5,m),PiD(6,m), PiD(7,m),&
                                          PiD(8,m), PiD(9,m),PiD(10,m),PiD(11,m), PiD(12,m),PiD(13,m))
        enddo
        !
        print *, '>>>>', outfilename
      endif
      !
      ! TODO: clean
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_b1)
      call fftw_free(c_b2)
      call fftw_free(c_b3)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_b1_filted)
      call fftw_free(c_b2_filted)
      call fftw_free(c_b3_filted)
      call fftw_free(c_A_filted)
      call fftw_free(c_C_filted)
      call fftw_free(c_H_filted)
      call fftw_free(c_term)
      call mpistop
      deallocate(All_filted_l,S_filted_l)
      deallocate(All_filted,S_filted,W_filted)
      deallocate(C_filted_l,Sigma_filted_l,J_filted_l,Sigma_filted,J_filted)
      deallocate(k1,k2,k3,ksq,Galpha,Gl,Gphi)
      deallocate(l_lim,l_sqrtalpha,l_phi,dl_alpha,num_alphas)
      deallocate(PiI,PiM,PiD,PiA,Pirank,Pisum)
      !
    end subroutine SGSPiB3Dint
    !
    subroutine SGSstress2D(thefilenumb)
      ! ! TODO : Improve need / Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,m,n,p,q
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname,mname
      character(len=10):: termname
      real(8), allocatable, dimension(:,:) :: k1,k2
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1,w2,rhol
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1_filted,w2_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1,w1w2,w2w1,w2w2
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: w1w1_filted,w1w2_filted,w2w1_filted,w2w2_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: A11,A12,A21,A22
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:) :: tau11_term,tau12_term,tau21_term,tau22_term
      real(8), allocatable, dimension(:,:,:,:) :: tau ! 1:2,1:2,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:,:,:) :: tau_bis ! 1:2,1:2,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:) :: errormax,erroravg, errorgtr10,errorgtr100
      real(8) :: result,norm2,norm2bis
      real(8) :: errornorm2max,errornorm2avg,errornorm2gtr10,errornorm2gtr100
      !
      !
      type(C_PTR) :: forward_plan,backward_plan
      type(C_PTR) :: c_w1,c_w2,c_rhol
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w2w1,c_w2w2
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w2w1_filted,c_w2w2_filted
      type(C_PTR) :: c_A11,c_A12,c_A21,c_A22
      type(C_PTR) :: c_tau11_term,c_tau12_term,c_tau21_term,c_tau22_term
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      allocate(erroravg(1:3,1:3),errormax(1:3,1:3),errorgtr10(1:3,1:3),errorgtr100(1:3,1:3))
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:2), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !! wavenumber
      allocate(k1(1:im,1:jm),k2(1:im,1:jm))
      call GenerateWave(im,jm,ia,ja,j0f,k1,k2)
      !
      allocate(tau(1:2,1:2,1:im,1:jm),tau_bis(1:2,1:2,1:im,1:jm))
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      do i=1,num_l
        dl_alpha(i,1) = l_sqrtalpha(i,1)**2 
        !
        do j=2,num_alphas(i)
          dl_alpha(i,j) = l_sqrtalpha(i,j)**2 -l_sqrtalpha(i,j-1)**2 
        enddo
      enddo
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw])
      c_rhol = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhol, rhol, [imfftw,jmfftw])
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11, [imfftw,jmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12, [imfftw,jmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21, [imfftw,jmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22, [imfftw,jmfftw])
      !
      c_tau11_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau11_term, tau11_term, [imfftw,jmfftw])
      c_tau12_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau12_term, tau12_term, [imfftw,jmfftw])
      c_tau21_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau21_term, tau21_term, [imfftw,jmfftw])
      c_tau22_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau22_term, tau22_term, [imfftw,jmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_2d(jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      ! wi,wiwj,rhol in physical space, not filted
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=CMPLX(vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2(i,j)=CMPLX(vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        w1w1(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w1w2(i,j)=CMPLX(vel(i,j,0,1)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w1(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,1)*rho(i,j,0),0.d0,C_INTPTR_T);
        w2w2(i,j)=CMPLX(vel(i,j,0,2)*vel(i,j,0,2)*rho(i,j,0),0.d0,C_INTPTR_T);
        rhol(i,j)=CMPLX(rho(i,j,0),0.d0,C_INTPTR_T);
        !
      end do
      end do
      !
      ! wi=rho ui,wiwj =rho ui uj,rhol in spectral space, not filted
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,rhol,rhol)
      !
      do j=1,jm
      do i=1,im
        !
        w1(i,j)=w1(i,j)/(1.d0*ia*ja)
        w2(i,j)=w2(i,j)/(1.d0*ia*ja)
        !
        w1w1(i,j)=w1w1(i,j)/(1.d0*ia*ja)
        w1w2(i,j)=w1w2(i,j)/(1.d0*ia*ja)
        w2w1(i,j)=w2w1(i,j)/(1.d0*ia*ja)
        w2w2(i,j)=w2w2(i,j)/(1.d0*ia*ja)
        !
        rhol(i,j)=rhol(i,j)/(1.d0*ia*ja)
        !
      end do
      end do
      !
      do m=1,num_l
        !
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        ! Method 1: tauij = rho uiuj - rho ui uj
        ! Filter scale: l
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,wiwj=rho ui uj,rho_filted in spectral space, filted by l
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j)**2+k2(i,j)**2)*l_lim(m)**2/2.d0) ! 
          !
          w1_filted(i,j)    = w1(i,j)   *Gl
          w2_filted(i,j)    = w2(i,j)   *Gl
          !
          w1w1_filted(i,j)  = w1w1(i,j) *Gl
          w1w2_filted(i,j)  = w1w2(i,j) *Gl
          w2w1_filted(i,j)  = w2w1(i,j) *Gl
          w2w2_filted(i,j)  = w2w2(i,j) *Gl
          !
          rho_filted(i,j)   = rhol(i,j) *Gl
        enddo
        enddo
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(rho ui)_filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(ui)~filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
          w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          !
          tau(1,1,i,j) = dreal(w1w1_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w1_filted(i,j))
          tau(1,2,i,j) = dreal(w1w2_filted(i,j) - rho_filted(i,j) * w1_filted(i,j) * w2_filted(i,j))
          tau(2,1,i,j) = dreal(w2w1_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w1_filted(i,j))
          tau(2,2,i,j) = dreal(w2w2_filted(i,j) - rho_filted(i,j) * w2_filted(i,j) * w2_filted(i,j))
          !
          do p=1,2
          do q=1,2
            tau_bis(p,q,i,j) = 0.d0
          enddo
          enddo
          !
        enddo
        enddo
        !
        ! Method 2: tauij = int_0^l2 （rho_√α  Aik_√α  Ajk_√α_√l2-α 
        ! Filter scale: l
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in spectral space, filted by sqrtalpha
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j)**2+k2(i,j)**2)*l_sqrtalpha(m,n)**2/2.d0) ! Filtre scale :sqrtalpha
            w1_filted(i,j)  = w1(i,j)  *Galpha
            w2_filted(i,j)  = w2(i,j)  *Galpha
            rho_filted(i,j) = rhol(i,j)*Galpha
          enddo
          enddo
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in physical space, filted by sqrtalpha
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! wi_filted=(ui)~filted, rho_filted in physical space, filted by sqrtalpha 
          do j=1,jm
          do i=1,im
            w1_filted(i,j) = w1_filted(i,j)/rho_filted(i,j)
            w2_filted(i,j) = w2_filted(i,j)/rho_filted(i,j)
          enddo
          enddo
          !
          ! wi_filted=(ui)~filted, Aij = Aij~filted in spectral space, filted by sqrtalpha 
          ! rho_filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j)  = w1_filted(i,j)/(1.d0*ia*ja)
            w2_filted(i,j)  = w2_filted(i,j)/(1.d0*ia*ja)
            !
            A11(i,j) = imag*w1_filted(i,j)*k1(i,j)
            A21(i,j) = imag*w2_filted(i,j)*k1(i,j)
            A12(i,j) = imag*w1_filted(i,j)*k2(i,j)
            A22(i,j) = imag*w2_filted(i,j)*k2(i,j)
            !
          end do
          end do
          !
          ! wi_filted=(ui)~filted in spectral space, filted by sqrtalpha 
          ! rho_filted, Aij = Aij~filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(backward_plan,A11,A11)
          call fftw_mpi_execute_dft(backward_plan,A21,A21)
          call fftw_mpi_execute_dft(backward_plan,A12,A12)
          call fftw_mpi_execute_dft(backward_plan,A22,A22)
          !
          do j=1,jm
          do i=1,im
            !
            tau11_term(i,j) = rho_filted(i,j) * (A11(i,j)*A11(i,j)+A12(i,j)*A12(i,j))
            tau12_term(i,j) = rho_filted(i,j) * (A11(i,j)*A21(i,j)+A12(i,j)*A22(i,j))
            tau21_term(i,j) = rho_filted(i,j) * (A21(i,j)*A11(i,j)+A22(i,j)*A12(i,j))
            tau22_term(i,j) = rho_filted(i,j) * (A21(i,j)*A21(i,j)+A22(i,j)*A22(i,j))
            !
          end do
          end do
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(forward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(forward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(forward_plan,tau22_term,tau22_term)
          !
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j)**2+k2(i,j)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            tau11_term(i,j) = tau11_term(i,j)*Gphi/(1.d0*ia*ja)
            tau12_term(i,j) = tau12_term(i,j)*Gphi/(1.d0*ia*ja)
            tau21_term(i,j) = tau21_term(i,j)*Gphi/(1.d0*ia*ja)
            tau22_term(i,j) = tau22_term(i,j)*Gphi/(1.d0*ia*ja)
            !
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(backward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(backward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(backward_plan,tau22_term,tau22_term)
          !
          !
          do j=1,jm
          do i=1,im
            !
            tau_bis(1,1,i,j) = tau_bis(1,1,i,j) + dreal(tau11_term(i,j)) * dl_alpha(m,n)
            tau_bis(1,2,i,j) = tau_bis(1,2,i,j) + dreal(tau12_term(i,j)) * dl_alpha(m,n)
            tau_bis(2,1,i,j) = tau_bis(2,1,i,j) + dreal(tau21_term(i,j)) * dl_alpha(m,n)
            tau_bis(2,2,i,j) = tau_bis(2,2,i,j) + dreal(tau22_term(i,j)) * dl_alpha(m,n)
            !
          enddo
          enddo
          !
        enddo ! loop of integral (alpha)
        !
        do p=1,2
        do q=1,2
          erroravg(p,q)=0.d0
          errormax(p,q)=0.d0
          errorgtr10(p,q)=0.d0
          errorgtr100(p,q)=0.d0
          errornorm2max=0.d0
          errornorm2avg=0.d0
          errornorm2gtr10=0.d0
          errornorm2gtr100=0.d0
        enddo
        enddo
        !
        ! Output comparaison results
        do j=1,jm
        do i=1,im
          !
          norm2 = 0.d0
          norm2bis = 0.d0
          do p=1,2
          do q=1,2
            if(abs(tau(p,q,i,j))>1.d-6)then
              result = abs(tau_bis(p,q,i,j)-tau(p,q,i,j))/(abs(tau_bis(p,q,i,j))+abs(tau(p,q,i,j)))
              norm2 = norm2 + tau(p,q,i,j)**2
              norm2bis = norm2bis + tau_bis(p,q,i,j)**2
              errormax(p,q) = max(errormax(p,q),result)
              erroravg(p,q) = erroravg(p,q) + result
              if(result > 0.1)then
                errorgtr10(p,q) = errorgtr10(p,q) + 1.d0
              endif
              if(result > 1)then
                errorgtr100(p,q) = errorgtr100(p,q) + 1.d0
              endif
            endif
          enddo
          enddo
          result = abs(norm2bis-norm2)/(abs(norm2)+abs(norm2bis))
          errornorm2max = max(errornorm2max,result)
          errornorm2avg = errornorm2avg + result
          if(result > 0.1)then
            errornorm2gtr10 = errornorm2gtr10 + 1.d0
          endif
          if(result > 1)then
            errornorm2gtr100 = errornorm2gtr100 + 1.d0
          endif
          !
        enddo
        enddo
        !
        !
        do p=1,2
        do q=1,2
          errormax(p,q) = pmax(errormax(p,q))
          erroravg(p,q) = psum(erroravg(p,q))/(1.d0*ia*ja)
          errorgtr10(p,q) = psum(errorgtr10(p,q))/(1.d0*ia*ja)
          errorgtr100(p,q) = psum(errorgtr100(p,q))/(1.d0*ia*ja)
          errornorm2max = psum(errornorm2max)
          errornorm2avg = psum(errornorm2avg)/(1.d0*ia*ja)
          errornorm2gtr10 = psum(errornorm2gtr10)/(1.d0*ia*ja)
          errornorm2gtr100 = psum(errornorm2gtr100)/(1.d0*ia*ja)
        enddo
        enddo
        !
        if(mpirank==0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_relative_error_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_stress_relative_error_'//mname//'.dat'
          endif
          !
          if(mpirank == 0)then
            open(fh,file=outfilename,form='formatted')
            write(fh,"(A7,1x,2(A20,1x))")'nstep','time','l'
            write(fh,"(I7,1x,2(E20.13E2,1x))")nstep,time,l_lim(m)
            write(fh,"(A8,1x,5(A20,1x))")'type','tau11','tau12','tau21','tau22','norm2'
            write(fh,"(A8,1x,5(E20.13E2,1x))")'max',((errormax(p,q),q=1,2),p=1,2),errornorm2max
            write(fh,"(A8,1x,5(E20.13E2,1x))")'avg',((erroravg(p,q),q=1,2),p=1,2),errornorm2avg
            write(fh,"(A8,1x,5(E20.13E2,1x))")'gtr0.1',((errorgtr10(p,q),q=1,2),p=1,2),errornorm2gtr10
            write(fh,"(A8,1x,5(E20.13E2,1x))")'gtr1',((errorgtr100(p,q),q=1,2),p=1,2),errornorm2gtr100
            close(fh)
            print *, '>>>>', outfilename
          endif
          !
          !
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(loutput)then
          !
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_'//stepname//'_'//mname//'.h5'
          else
            outfilename = 'pp/SGS_stress_'//mname//'.h5'
          endif
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          do p=1,2
          do q= 1,2
            write (termname, "(A3,I1,I1)") "tau",p,q
            call h5wa2d_r8(varname=termname,var=tau(p,q,1:im,1:jm),    dir='k')
            write (termname, "(A3,I1,I1,A3)") "tau",p,q,"bis"
            call h5wa2d_r8(varname=termname,var=tau_bis(p,q,1:im,1:jm),dir='k')
          enddo
          enddo
          !
          call h5io_end
          !
        endif
        !
      enddo ! loop of filter point l
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_rhol)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_tau11_term)
      call fftw_free(c_tau12_term)
      call fftw_free(c_tau21_term)
      call fftw_free(c_tau22_term)
      call mpistop
      deallocate(tau,tau_bis,erroravg,errormax,errorgtr10,errorgtr100)
      !
    end subroutine SGSstress2D
    !
    subroutine SGSstress3D(thefilenumb)
      ! ! TODO : Improve need/ Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: fh
      integer :: i,j,k,m,n,p,q
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname,mname
      character(len=10):: termname
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: l_lim
      real(8),allocatable,dimension(:,:) :: l_sqrtalpha,l_phi,dl_alpha
      integer,allocatable,dimension(:) :: num_alphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max, ratio_min
      real(8) :: Gl,Galpha,Gphi
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhol
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1,w1w2,w1w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1,w2w2,w2w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1,w3w2,w3w3
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1w1_filted,w1w2_filted,w1w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w2w1_filted,w2w2_filted,w2w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w3w1_filted,w3w2_filted,w3w3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11,A12,A13
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21,A22,A23
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31,A32,A33
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau11_term,tau12_term,tau13_term
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau21_term,tau22_term,tau23_term
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: tau31_term,tau32_term,tau33_term
      real(8), allocatable, dimension(:,:,:,:,:) :: tau ! 1:3,1:3,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:,:,:,:) :: tau_bis ! 1:3,1:3,1:im,1:jm,1:km
      real(8), allocatable, dimension(:,:) :: errormax,erroravg, errorgtr10,errorgtr100
      real(8) :: result,norm2,norm2bis
      real(8) :: errornorm2max,errornorm2avg,errornorm2gtr10,errornorm2gtr100
      !
      !
      type(C_PTR) :: forward_plan,backward_plan
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhol
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted
      type(C_PTR) :: c_w1w1,c_w1w2,c_w1w3
      type(C_PTR) :: c_w2w1,c_w2w2,c_w2w3
      type(C_PTR) :: c_w3w1,c_w3w2,c_w3w3
      type(C_PTR) :: c_w1w1_filted,c_w1w2_filted,c_w1w3_filted
      type(C_PTR) :: c_w2w1_filted,c_w2w2_filted,c_w2w3_filted
      type(C_PTR) :: c_w3w1_filted,c_w3w2_filted,c_w3w3_filted
      type(C_PTR) :: c_A11,c_A12,c_A13
      type(C_PTR) :: c_A21,c_A22,c_A23
      type(C_PTR) :: c_A31,c_A32,c_A33
      type(C_PTR) :: c_tau11_term,c_tau12_term,c_tau13_term,c_tau21_term,c_tau22_term
      type(C_PTR) :: c_tau23_term,c_tau31_term,c_tau32_term,c_tau33_term
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      allocate(erroravg(1:3,1:3),errormax(1:3,1:3),errorgtr10(1:3,1:3),errorgtr100(1:3,1:3))
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      allocate(tau(1:3,1:3,1:im,1:jm,1:km),tau_bis(1:3,1:3,1:im,1:jm,1:km))
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare l,alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(l_lim(1:num_l),num_alphas(1:num_l),l_sqrtalpha(1:num_l,1:num_alpha))
      allocate(l_phi(1:num_l,1:num_alpha),dl_alpha(1:num_l,1:num_alpha))
      !
      call SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      if(mpirank==0) then
        open(fh,file='pp/SGSintegral.info',form='formatted')
        write(fh,"(2(A9,1x))")'NumL','NumAlpha'
        write(fh,"(2(I9,1x))")num_l,num_alpha
        write(fh,"(2(A9,1x),2(A15,1x))")'i','j','l_lim','l_sqrtalpha'
        do i=1,num_l
          do j=1,num_alphas(i)
          ! Output file of rank information.
            write(fh,"(2(I9,1x),2(E15.7E3,1x))")i,j,l_lim(i),l_sqrtalpha(i,j)
          enddo
        enddo
        !
        close(fh)
        print*,' << SGSintegral.info ... done !'
      endif
      !
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhol = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhol, rhol, [imfftw,jmfftw,kmfftw])
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted, [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted, [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted, [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1, w1w1, [imfftw,jmfftw,kmfftw])
      c_w1w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2, w1w2, [imfftw,jmfftw,kmfftw])
      c_w1w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3, w1w3, [imfftw,jmfftw,kmfftw])
      c_w2w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1, w2w1, [imfftw,jmfftw,kmfftw])
      c_w2w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2, w2w2, [imfftw,jmfftw,kmfftw])
      c_w2w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3, w2w3, [imfftw,jmfftw,kmfftw])
      c_w3w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1, w3w1, [imfftw,jmfftw,kmfftw])
      c_w3w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2, w3w2, [imfftw,jmfftw,kmfftw])
      c_w3w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3, w3w3, [imfftw,jmfftw,kmfftw])
      !
      c_w1w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w1_filted, w1w1_filted, [imfftw,jmfftw,kmfftw])
      c_w1w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w2_filted, w1w2_filted, [imfftw,jmfftw,kmfftw])
      c_w1w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1w3_filted, w1w3_filted, [imfftw,jmfftw,kmfftw])
      c_w2w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w1_filted, w2w1_filted, [imfftw,jmfftw,kmfftw])
      c_w2w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w2_filted, w2w2_filted, [imfftw,jmfftw,kmfftw])
      c_w2w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2w3_filted, w2w3_filted, [imfftw,jmfftw,kmfftw])
      c_w3w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w1_filted, w3w1_filted, [imfftw,jmfftw,kmfftw])
      c_w3w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w2_filted, w3w2_filted, [imfftw,jmfftw,kmfftw])
      c_w3w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3w3_filted, w3w3_filted, [imfftw,jmfftw,kmfftw])
      !
      c_A11 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11, A11, [imfftw,jmfftw,kmfftw])
      c_A12 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12, A12, [imfftw,jmfftw,kmfftw])
      c_A13 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13, A13, [imfftw,jmfftw,kmfftw])
      c_A21 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21, A21, [imfftw,jmfftw,kmfftw])
      c_A22 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22, A22, [imfftw,jmfftw,kmfftw])
      c_A23 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23, A23, [imfftw,jmfftw,kmfftw])
      c_A31 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31, A31, [imfftw,jmfftw,kmfftw])
      c_A32 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32, A32, [imfftw,jmfftw,kmfftw])
      c_A33 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33, A33, [imfftw,jmfftw,kmfftw])
      !
      c_tau11_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau11_term, tau11_term, [imfftw,jmfftw,kmfftw])
      c_tau12_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau12_term, tau12_term, [imfftw,jmfftw,kmfftw])
      c_tau13_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau13_term, tau13_term, [imfftw,jmfftw,kmfftw])
      c_tau21_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau21_term, tau21_term, [imfftw,jmfftw,kmfftw])
      c_tau22_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau22_term, tau22_term, [imfftw,jmfftw,kmfftw])
      c_tau23_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau23_term, tau23_term, [imfftw,jmfftw,kmfftw])
      c_tau31_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau31_term, tau31_term, [imfftw,jmfftw,kmfftw])
      c_tau32_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau32_term, tau32_term, [imfftw,jmfftw,kmfftw])
      c_tau33_term = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_tau33_term, tau33_term, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      ! wi,wiwj,rhol in physical space, not filted
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w1(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w2(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w1w3(i,j,k)=CMPLX(vel(i,j,k,1)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w1(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w2(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w2w3(i,j,k)=CMPLX(vel(i,j,k,2)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w1(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w2(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T);
        w3w3(i,j,k)=CMPLX(vel(i,j,k,3)*vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T);
        rhol(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T);
        !
      end do
      end do
      end do
      !
      ! wi=rho ui,wiwj =rho ui uj,rhol in spectral space, not filted
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,w1w1,w1w1)
      call fftw_mpi_execute_dft(forward_plan,w1w2,w1w2)
      call fftw_mpi_execute_dft(forward_plan,w1w3,w1w3)
      call fftw_mpi_execute_dft(forward_plan,w2w1,w2w1)
      call fftw_mpi_execute_dft(forward_plan,w2w2,w2w2)
      call fftw_mpi_execute_dft(forward_plan,w2w3,w2w3)
      call fftw_mpi_execute_dft(forward_plan,w3w1,w3w1)
      call fftw_mpi_execute_dft(forward_plan,w3w2,w3w2)
      call fftw_mpi_execute_dft(forward_plan,w3w3,w3w3)
      call fftw_mpi_execute_dft(forward_plan,rhol,rhol)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        w1w1(i,j,k)=w1w1(i,j,k)/(1.d0*ia*ja*ka)
        w1w2(i,j,k)=w1w2(i,j,k)/(1.d0*ia*ja*ka)
        w1w3(i,j,k)=w1w3(i,j,k)/(1.d0*ia*ja*ka)
        w2w1(i,j,k)=w2w1(i,j,k)/(1.d0*ia*ja*ka)
        w2w2(i,j,k)=w2w2(i,j,k)/(1.d0*ia*ja*ka)
        w2w3(i,j,k)=w2w3(i,j,k)/(1.d0*ia*ja*ka)
        w3w1(i,j,k)=w3w1(i,j,k)/(1.d0*ia*ja*ka)
        w3w2(i,j,k)=w3w2(i,j,k)/(1.d0*ia*ja*ka)
        w3w3(i,j,k)=w3w3(i,j,k)/(1.d0*ia*ja*ka)
        !
        rhol(i,j,k)=rhol(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
      !
      do m=1,num_l
        !
        if(mpirank==0)  print *, '* l = ', l_lim(m) ,' at', m, '/', num_l
        !
        ! Method 1: tauij = rho uiuj - rho ui uj
        ! Filter scale: l
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=rho ui,wiwj=rho ui uj,rho_filted in spectral space, filted by l
        do k=1,km
        do j=1,jm
        do i=1,im
          Gl = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_lim(m)**2/2.d0)
          !
          w1_filted(i,j,k)    = w1(i,j,k)   *Gl
          w2_filted(i,j,k)    = w2(i,j,k)   *Gl
          w3_filted(i,j,k)    = w3(i,j,k)   *Gl
          !
          w1w1_filted(i,j,k)  = w1w1(i,j,k) *Gl
          w1w2_filted(i,j,k)  = w1w2(i,j,k) *Gl
          w1w3_filted(i,j,k)  = w1w3(i,j,k) *Gl
          w2w1_filted(i,j,k)  = w2w1(i,j,k) *Gl
          w2w2_filted(i,j,k)  = w2w2(i,j,k) *Gl
          w2w3_filted(i,j,k)  = w2w3(i,j,k) *Gl
          w3w1_filted(i,j,k)  = w3w1(i,j,k) *Gl
          w3w2_filted(i,j,k)  = w3w2(i,j,k) *Gl
          w3w3_filted(i,j,k)  = w3w3(i,j,k) *Gl
          !
          rho_filted(i,j,k)   = rhol(i,j,k) *Gl
        enddo
        enddo
        enddo
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(rho ui)_filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w1_filted,w1w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w2_filted,w1w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w1w3_filted,w1w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w1_filted,w2w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w2_filted,w2w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w2w3_filted,w2w3_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w1_filted,w3w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w2_filted,w3w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3w3_filted,w3w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        !
        ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
        ! wi_filted=(ui)~filted,wiwj=(rho ui uj)_filted,rho_filted in physical space, filted by l 
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
          tau(1,1,i,j,k) = dreal(w1w1_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w1_filted(i,j,k))
          tau(1,2,i,j,k) = dreal(w1w2_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w2_filted(i,j,k))
          tau(1,3,i,j,k) = dreal(w1w3_filted(i,j,k) - rho_filted(i,j,k) * w1_filted(i,j,k) * w3_filted(i,j,k))
          tau(2,1,i,j,k) = dreal(w2w1_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w1_filted(i,j,k))
          tau(2,2,i,j,k) = dreal(w2w2_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w2_filted(i,j,k))
          tau(2,3,i,j,k) = dreal(w2w3_filted(i,j,k) - rho_filted(i,j,k) * w2_filted(i,j,k) * w3_filted(i,j,k))
          tau(3,1,i,j,k) = dreal(w3w1_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w1_filted(i,j,k))
          tau(3,2,i,j,k) = dreal(w3w2_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w2_filted(i,j,k))
          tau(3,3,i,j,k) = dreal(w3w3_filted(i,j,k) - rho_filted(i,j,k) * w3_filted(i,j,k) * w3_filted(i,j,k))
          !
          do p=1,3
          do q=1,3
            tau_bis(p,q,i,j,k) = 0.d0
          enddo
          enddo
          !
        enddo
        enddo
        enddo
        !
        ! Method 2: tauij = int_0^l2 （rho_√α  Aik_√α  Ajk_√α_√l2-α 
        ! Filter scale: l
        !
        do n=1,num_alphas(m)
          !
          call date_and_time(values=value) 
          !
          if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alphas(m),',now is ',&
                                  value(5), ':', value(6),':',value(7)
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in spectral space, filted by sqrtalpha
          do k=1,km
          do j=1,jm
          do i=1,im
            Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_sqrtalpha(m,n)**2/2.d0)
            w1_filted(i,j,k)  = w1(i,j,k)  *Galpha
            w2_filted(i,j,k)  = w2(i,j,k)  *Galpha
            w3_filted(i,j,k)  = w3(i,j,k)  *Galpha
            rho_filted(i,j,k) = rhol(i,j,k)*Galpha
          enddo
          enddo
          enddo
          !
          ! wi=rho ui,wiwj=rho ui uj in spectral space, not filted
          ! wi_filted=rho ui,rho_filted in physical space, filted by sqrtalpha
          call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
          call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
          !
          ! wi_filted=(ui)~filted, rho_filted in physical space, filted by sqrtalpha 
          do k=1,km
          do j=1,jm
          do i=1,im
            w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
            w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
            w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          enddo
          enddo
          enddo
          !
          ! wi_filted=(ui)~filted, Aij = Aij~filted in spectral space, filted by sqrtalpha 
          ! rho_filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
          call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
          call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            w1_filted(i,j,k)  = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
            w2_filted(i,j,k)  = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
            w3_filted(i,j,k)  = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
            !
            A11(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
            A21(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
            A31(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
            A12(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
            A22(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
            A32(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
            A13(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
            A23(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
            A33(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
            !
          end do
          end do
          end do
          !
          ! wi_filted=(ui)~filted in spectral space, filted by sqrtalpha 
          ! rho_filted, Aij = Aij~filted in physical space, filted by sqrtalpha 
          call fftw_mpi_execute_dft(backward_plan,A11,A11)
          call fftw_mpi_execute_dft(backward_plan,A21,A21)
          call fftw_mpi_execute_dft(backward_plan,A31,A31)
          call fftw_mpi_execute_dft(backward_plan,A12,A12)
          call fftw_mpi_execute_dft(backward_plan,A22,A22)
          call fftw_mpi_execute_dft(backward_plan,A32,A32)
          call fftw_mpi_execute_dft(backward_plan,A13,A13)
          call fftw_mpi_execute_dft(backward_plan,A23,A23)
          call fftw_mpi_execute_dft(backward_plan,A33,A33)
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            tau11_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A11(i,j,k)+A12(i,j,k)*A12(i,j,k)+A13(i,j,k)*A13(i,j,k))
            tau12_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A21(i,j,k)+A12(i,j,k)*A22(i,j,k)+A13(i,j,k)*A23(i,j,k))
            tau13_term(i,j,k) = rho_filted(i,j,k) * (A11(i,j,k)*A31(i,j,k)+A12(i,j,k)*A32(i,j,k)+A13(i,j,k)*A33(i,j,k))
            tau21_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A11(i,j,k)+A22(i,j,k)*A12(i,j,k)+A23(i,j,k)*A13(i,j,k))
            tau22_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A21(i,j,k)+A22(i,j,k)*A22(i,j,k)+A23(i,j,k)*A23(i,j,k))
            tau23_term(i,j,k) = rho_filted(i,j,k) * (A21(i,j,k)*A31(i,j,k)+A22(i,j,k)*A32(i,j,k)+A23(i,j,k)*A33(i,j,k))
            tau31_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A11(i,j,k)+A32(i,j,k)*A12(i,j,k)+A33(i,j,k)*A13(i,j,k))
            tau32_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A21(i,j,k)+A32(i,j,k)*A22(i,j,k)+A33(i,j,k)*A23(i,j,k))
            tau33_term(i,j,k) = rho_filted(i,j,k) * (A31(i,j,k)*A31(i,j,k)+A32(i,j,k)*A32(i,j,k)+A33(i,j,k)*A33(i,j,k))
            !
          end do
          end do
          end do
          !
          ! Do filter phi:
          ! F -> product -> F inverse
          call fftw_mpi_execute_dft(forward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(forward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(forward_plan,tau13_term,tau13_term)
          call fftw_mpi_execute_dft(forward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(forward_plan,tau22_term,tau22_term)
          call fftw_mpi_execute_dft(forward_plan,tau23_term,tau23_term)
          call fftw_mpi_execute_dft(forward_plan,tau31_term,tau31_term)
          call fftw_mpi_execute_dft(forward_plan,tau32_term,tau32_term)
          call fftw_mpi_execute_dft(forward_plan,tau33_term,tau33_term)
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            Gphi = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*l_phi(m,n)**2/2.d0) ! Filtre scale :phi
            tau11_term(i,j,k) = tau11_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau12_term(i,j,k) = tau12_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau13_term(i,j,k) = tau13_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau21_term(i,j,k) = tau21_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau22_term(i,j,k) = tau22_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau23_term(i,j,k) = tau23_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau31_term(i,j,k) = tau31_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau32_term(i,j,k) = tau32_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            tau33_term(i,j,k) = tau33_term(i,j,k)*Gphi/(1.d0*ia*ja*ka)
            !
          enddo
          enddo
          enddo
          !
          !
          call fftw_mpi_execute_dft(backward_plan,tau11_term,tau11_term)
          call fftw_mpi_execute_dft(backward_plan,tau12_term,tau12_term)
          call fftw_mpi_execute_dft(backward_plan,tau13_term,tau13_term)
          call fftw_mpi_execute_dft(backward_plan,tau21_term,tau21_term)
          call fftw_mpi_execute_dft(backward_plan,tau22_term,tau22_term)
          call fftw_mpi_execute_dft(backward_plan,tau23_term,tau23_term)
          call fftw_mpi_execute_dft(backward_plan,tau31_term,tau31_term)
          call fftw_mpi_execute_dft(backward_plan,tau32_term,tau32_term)
          call fftw_mpi_execute_dft(backward_plan,tau33_term,tau33_term)
          !
          !
          do k=1,km
          do j=1,jm
          do i=1,im
            !
            tau_bis(1,1,i,j,k) = tau_bis(1,1,i,j,k) + dreal(tau11_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(1,2,i,j,k) = tau_bis(1,2,i,j,k) + dreal(tau12_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(1,3,i,j,k) = tau_bis(1,3,i,j,k) + dreal(tau13_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(2,1,i,j,k) = tau_bis(2,1,i,j,k) + dreal(tau21_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(2,2,i,j,k) = tau_bis(2,2,i,j,k) + dreal(tau22_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(2,3,i,j,k) = tau_bis(2,3,i,j,k) + dreal(tau23_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(3,1,i,j,k) = tau_bis(3,1,i,j,k) + dreal(tau31_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(3,2,i,j,k) = tau_bis(3,2,i,j,k) + dreal(tau32_term(i,j,k)) * dl_alpha(m,n)
            tau_bis(3,3,i,j,k) = tau_bis(3,3,i,j,k) + dreal(tau33_term(i,j,k)) * dl_alpha(m,n)
            !
          enddo
          enddo
          enddo
          !
        enddo ! loop of integral (alpha)
        !
        do p=1,3
        do q=1,3
          erroravg(p,q)=0.d0
          errormax(p,q)=0.d0
          errorgtr10(p,q)=0.d0
          errorgtr100(p,q)=0.d0
        enddo
        enddo
        errornorm2max=0.d0
        errornorm2avg=0.d0
        errornorm2gtr10=0.d0
        errornorm2gtr100=0.d0
        !
        ! Output comparaison results
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          norm2 = 0.d0
          norm2bis = 0.d0
          do p=1,3
          do q=1,3
            if(abs(tau(p,q,i,j,k))>2.d-4)then
              result = 2*abs(tau_bis(p,q,i,j,k)-tau(p,q,i,j,k))/(abs(tau_bis(p,q,i,j,k))+abs(tau(p,q,i,j,k)))
              norm2 = norm2 + tau(p,q,i,j,k)**2
              norm2bis = norm2bis + tau_bis(p,q,i,j,k)**2
              errormax(p,q) = max(errormax(p,q),result)
              erroravg(p,q) = erroravg(p,q) + result
              if(result > 0.1)then
                errorgtr10(p,q) = errorgtr10(p,q) + 1.d0
              endif
              if(result > 1)then
                errorgtr100(p,q) = errorgtr100(p,q) + 1.d0
                print *, 'p',p,'q',q,tau_bis(p,q,i,j,k),tau(p,q,i,j,k),result, '*'
              endif
            endif
          enddo
          enddo
          result = 2*abs(norm2bis-norm2)/(abs(norm2)+abs(norm2bis))
          errornorm2max = max(errornorm2max,result)
          errornorm2avg = errornorm2avg + result
          if(result > 0.1)then
            errornorm2gtr10 = errornorm2gtr10 + 1.d0
          endif
          if(result > 1)then
            errornorm2gtr100 = errornorm2gtr100 + 1.d0
            print *, norm2, norm2bis, result, '*'
          endif
          !
        enddo
        enddo
        enddo
        !
        !
        do p=1,3
        do q=1,3
          errormax(p,q) = pmax(errormax(p,q))
          erroravg(p,q) = psum(erroravg(p,q))/(1.d0*ia*ja*ka)
          errorgtr10(p,q) = psum(errorgtr10(p,q))/(1.d0*ia*ja*ka)
          errorgtr100(p,q) = psum(errorgtr100(p,q))/(1.d0*ia*ja*ka)
        enddo
        enddo
        errornorm2max = pmax(errornorm2max)
        errornorm2avg = psum(errornorm2avg)/(1.d0*ia*ja*ka)
        errornorm2gtr10 = psum(errornorm2gtr10)/(1.d0*ia*ja*ka)
        errornorm2gtr100 = psum(errornorm2gtr100)/(1.d0*ia*ja*ka)
        !
        if(mpirank==0) then
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_relative_error_'//stepname//'_'//mname//'.dat'
          else
            outfilename = 'pp/SGS_stress_relative_error_'//mname//'.dat'
          endif
          !
          if(mpirank == 0)then
            open(fh,file=outfilename,form='formatted')
            write(fh,"(A7,1x,2(A20,1x))")'nstep','time','l'
            write(fh,"(I7,1x,2(E20.13E2,1x))")nstep,time,l_lim(m)
            write(fh,"(A8,1x,10(A20,1x))")'type','tau11','tau12','tau13','tau21','tau22','tau23','tau31','tau32','tau33','norm2'
            write(fh,"(A8,1x,10(E20.13E2,1x))")'max',((errormax(p,q),q=1,3),p=1,3),errornorm2max
            write(fh,"(A8,1x,10(E20.13E2,1x))")'avg',((erroravg(p,q),q=1,3),p=1,3),errornorm2avg
            write(fh,"(A8,1x,10(E20.13E2,1x))")'gtr0.1',((errorgtr10(p,q),q=1,3),p=1,3),errornorm2gtr10
            write(fh,"(A8,1x,10(E20.13E2,1x))")'gtr1',((errorgtr100(p,q),q=1,3),p=1,3),errornorm2gtr100
            close(fh)
            print *, '>>>>', outfilename
          endif
          !
          !
          !
        endif
        !
        call mpi_barrier(mpi_comm_world,ierr)
        !
        if(loutput)then
          !
          write(mname,'(i4.4)')m
          if (thefilenumb .ne. 0) then
            outfilename = 'pp/SGS_stress_'//stepname//'_'//mname//'.h5'
          else
            outfilename = 'pp/SGS_stress_'//mname//'.h5'
          endif
          !
          call h5io_init(trim(outfilename),mode='write')
          !
          do p=1,3
            do q= 1,3
              write (termname, "(A3,I1,I1)") "tau",p,q
              call h5write(var=tau(p,q,1:im,1:jm,1:km),      varname=termname,    mode = modeio) 
              write (termname, "(A3,I1,I1,A3)") "tau",p,q,"bis"
              call h5write(var=tau_bis(p,q,1:im,1:jm,1:km),  varname=termname,    mode = modeio)
            enddo
          enddo
          !
          call h5io_end
          !
        endif
        !
      enddo ! loop of filter point l
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhol)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_w1w1)
      call fftw_free(c_w1w2)
      call fftw_free(c_w1w3)
      call fftw_free(c_w2w1)
      call fftw_free(c_w2w2)
      call fftw_free(c_w2w3)
      call fftw_free(c_w3w1)
      call fftw_free(c_w3w2)
      call fftw_free(c_w3w3)
      call fftw_free(c_w1w1_filted)
      call fftw_free(c_w1w2_filted)
      call fftw_free(c_w1w3_filted)
      call fftw_free(c_w2w1_filted)
      call fftw_free(c_w2w2_filted)
      call fftw_free(c_w2w3_filted)
      call fftw_free(c_w3w1_filted)
      call fftw_free(c_w3w2_filted)
      call fftw_free(c_w3w3_filted)
      call fftw_free(c_A11)
      call fftw_free(c_A12)
      call fftw_free(c_A13)
      call fftw_free(c_A21)
      call fftw_free(c_A22)
      call fftw_free(c_A23)
      call fftw_free(c_A31)
      call fftw_free(c_A32)
      call fftw_free(c_A33)
      call fftw_free(c_tau11_term)
      call fftw_free(c_tau12_term)
      call fftw_free(c_tau13_term)
      call fftw_free(c_tau21_term)
      call fftw_free(c_tau22_term)
      call fftw_free(c_tau23_term)
      call fftw_free(c_tau31_term)
      call fftw_free(c_tau32_term)
      call fftw_free(c_tau33_term)
      call mpistop
      deallocate(tau,tau_bis,erroravg,errormax,errorgtr10,errorgtr100)
      !
    end subroutine SGSstress3D
    !
    subroutine SGST3D(thefilenumb)
      ! ! TODO : Improve need/ Test need
      !
      use, intrinsic :: iso_c_binding
      use readwrite, only : readinput
      use fftwlink
      use commvar,only : time,nstep,im,jm,km,ia,ja,ka
      use commarray, only: vel, rho, prs
      use hdf5io
      use utility,  only : listinit,listwrite
      use parallel, only : bcast, pmax, pmin, psum, lio, parallelini,mpistop
      use solver, only: refcal
      include 'fftw3-mpi.f03'
      !
      integer,intent(in) :: thefilenumb
      integer :: i,j,k,m,n
      character(len=128) :: infilename,outfilename
      character(len=4) :: stepname
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1,w2,w3,rhocom,pcom
      real(8), allocatable, dimension(:,:,:) :: k1,k2,k3
      complex(8) :: imag
      real(8),allocatable,dimension(:) :: sqrtalphas,dalphas
      integer :: num_l,num_alpha,num_alphamin
      integer :: hand_a
      real(8) :: l_min, ratio_max,ratio_min
      real(8) :: Galpha
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: w1_filted,w2_filted,w3_filted,rho_filted,p_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A11_filted,A12_filted,A13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A21_filted,A22_filted,A23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: A31_filted,A32_filted,A33_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p11_filted,p12_filted,p13_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p21_filted,p22_filted,p23_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p31_filted,p32_filted,p33_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: p1_filted,p2_filted,p3_filted
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:,:,:) :: rho1_filted,rho2_filted,rho3_filted
      complex(8), allocatable, dimension(:,:,:) :: theta_filted,M11_filted,M22_filted,M33_filted
      complex(8), allocatable, dimension(:,:,:) :: M12_filted,M13_filted,M21_filted
      complex(8), allocatable, dimension(:,:,:) :: M23_filted,M31_filted,M32_filted
      real(8) :: Es,Ec,Ts,Tc,Ps1,Pc1,Pc2
      !
      complex(8) :: vxr_D
      !
      type(C_PTR) :: c_w1,c_w2,c_w3,c_rhocom,c_pcom,forward_plan,backward_plan
      type(C_PTR) :: c_w1_filted,c_w2_filted,c_w3_filted,c_rho_filted,c_p_filted
      type(C_PTR) :: c_A11_filted,c_A12_filted,c_A13_filted
      type(C_PTR) :: c_A21_filted,c_A22_filted,c_A23_filted
      type(C_PTR) :: c_A31_filted,c_A32_filted,c_A33_filted
      type(C_PTR) :: c_p11_filted,c_p12_filted,c_p13_filted
      type(C_PTR) :: c_p21_filted,c_p22_filted,c_p23_filted
      type(C_PTR) :: c_p31_filted,c_p32_filted,c_p33_filted
      type(C_PTR) :: c_p1_filted,c_p2_filted,c_p3_filted
      type(C_PTR) :: c_rho1_filted,c_rho2_filted,c_rho3_filted
      !
      integer,dimension(8) :: value
      character(len=1) :: modeio
      logical :: loutput
      !
      call readinput
      call refcal
      if(mpirank==0)  print*, '** refcal done!'
      !
      modeio='h'
      ! Initialization
      call fftw_mpi_init()
      if(mpirank==0)  print *, "fftw_mpi initialized"
      !
      if(mpirank==0)  print *, "ia:",ia,",ja:",ja,",ka:",ka
      !
      call mpisizedis_fftw
      if(mpirank==0)  print*, '** mpisizedis & parapp done!'
      !
      call parallelini
      if(mpirank==0)  print*, '** parallelini done!'
      !
      !!!! Read velocity and density field
      allocate(vel(0:im,0:jm,0:km,1:3), rho(0:im,0:jm,0:km),prs(0:im,0:jm,0:km))
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
      call h5read(varname='ro', var=rho(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='u1', var=vel(0:im,0:jm,0:km,1),mode = modeio)
      call h5read(varname='u2', var=vel(0:im,0:jm,0:km,2),mode = modeio)
      call h5read(varname='u3', var=vel(0:im,0:jm,0:km,3),mode = modeio)
      call h5read(varname='p',  var=prs(0:im,0:jm,0:km),  mode = modeio)
      call h5read(varname='time',var=time)
      call h5read(varname='nstep',var=nstep)
      !
      call h5io_end
      !
      call mpi_barrier(mpi_comm_world,ierr)
      !
      if(mpirank==0)  print *, "Field read finish!"
      !
      !!!! Prepare initial field in Fourier space
      !! velocity
      c_w1 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1, w1, [imfftw,jmfftw,kmfftw])
      c_w2 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2, w2, [imfftw,jmfftw,kmfftw])
      c_w3 = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3, w3, [imfftw,jmfftw,kmfftw])
      c_rhocom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rhocom, rhocom, [imfftw,jmfftw,kmfftw])
      c_pcom = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_pcom, pcom, [imfftw,jmfftw,kmfftw])
      !
      forward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE)
      backward_plan = fftw_mpi_plan_dft_3d(kafftw,jafftw,iafftw, w1,w1, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=CMPLX(vel(i,j,k,1)*rho(i,j,k),0.d0,C_INTPTR_T)
        w2(i,j,k)=CMPLX(vel(i,j,k,2)*rho(i,j,k),0.d0,C_INTPTR_T)
        w3(i,j,k)=CMPLX(vel(i,j,k,3)*rho(i,j,k),0.d0,C_INTPTR_T)
        rhocom(i,j,k)=CMPLX(rho(i,j,k),0.d0,C_INTPTR_T)
        pcom(i,j,k)=CMPLX(prs(i,j,k),0.d0,C_INTPTR_T)
        !
      end do
      end do
      end do
      !
      !After this bloc, w1 is (rho*u1) in spectral space
      call fftw_mpi_execute_dft(forward_plan,w1,w1)
      call fftw_mpi_execute_dft(forward_plan,w2,w2)
      call fftw_mpi_execute_dft(forward_plan,w3,w3)
      call fftw_mpi_execute_dft(forward_plan,rhocom,rhocom)
      call fftw_mpi_execute_dft(forward_plan,pcom,pcom)
      !
      do k=1,km
      do j=1,jm
      do i=1,im
        !
        w1(i,j,k)=w1(i,j,k)/(1.d0*ia*ja*ka)
        w2(i,j,k)=w2(i,j,k)/(1.d0*ia*ja*ka)
        w3(i,j,k)=w3(i,j,k)/(1.d0*ia*ja*ka)
        rhocom(i,j,k)=rhocom(i,j,k)/(1.d0*ia*ja*ka)
        pcom(i,j,k)=pcom(i,j,k)/(1.d0*ia*ja*ka)
        !
      end do
      end do
      end do
  
      !
      !
      !! wavenumber
      allocate(k1(1:im,1:jm,1:km),k2(1:im,1:jm,1:km),k3(1:im,1:jm,1:km))
      call GenerateWave(im,jm,km,ia,ja,ka,k0f,k1,k2,k3)
      !
      !! Imaginary number prepare
      imag = CMPLX(0.d0,1.d0,8)
      !
      if(mpirank==0)  print *, "Velocity field and wavenum prepare finish"
      !!!! Prepare alpha and others
      call readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      l_min = 2*pi/ia
      allocate(sqrtalphas(num_alpha),dalphas(num_alpha))
      !
      do i=1,num_alpha
        sqrtalphas(i) = sqrt( exp(log(ratio_max**2) * (i-1) / (num_alpha-1)) ) * l_min
      enddo
      !
      dalphas(1) = sqrtalphas(1)**2 
      !
      do i=2,num_alpha
        dalphas(i) = sqrtalphas(i)**2 - sqrtalphas(i-1)**2 
      enddo
      !
      if(mpirank==0)  print *, "Integrate point allocated"
      !
      !!!!
      !
      c_w1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w1_filted, w1_filted,  [imfftw,jmfftw,kmfftw])
      c_w2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w2_filted, w2_filted,  [imfftw,jmfftw,kmfftw])
      c_w3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_w3_filted, w3_filted,  [imfftw,jmfftw,kmfftw])
      c_rho_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho_filted, rho_filted,[imfftw,jmfftw,kmfftw])
      c_p_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p_filted, p_filted,[imfftw,jmfftw,kmfftw])
      !
      c_A11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A11_filted, A11_filted,[imfftw,jmfftw,kmfftw])
      c_A12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A12_filted, A12_filted,[imfftw,jmfftw,kmfftw])
      c_A13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A13_filted, A13_filted,[imfftw,jmfftw,kmfftw])
      c_A21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A21_filted, A21_filted,[imfftw,jmfftw,kmfftw])
      c_A22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A22_filted, A22_filted,[imfftw,jmfftw,kmfftw])
      c_A23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A23_filted, A23_filted,[imfftw,jmfftw,kmfftw])
      c_A31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A31_filted, A31_filted,[imfftw,jmfftw,kmfftw])
      c_A32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A32_filted, A32_filted,[imfftw,jmfftw,kmfftw])
      c_A33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_A33_filted, A33_filted,[imfftw,jmfftw,kmfftw])
      !
      c_p11_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p11_filted, p11_filted,[imfftw,jmfftw,kmfftw])
      c_p12_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p12_filted, p12_filted,[imfftw,jmfftw,kmfftw])
      c_p13_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p13_filted, p13_filted,[imfftw,jmfftw,kmfftw])
      c_p21_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p21_filted, p21_filted,[imfftw,jmfftw,kmfftw])
      c_p22_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p22_filted, p22_filted,[imfftw,jmfftw,kmfftw])
      c_p23_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p23_filted, p23_filted,[imfftw,jmfftw,kmfftw])
      c_p31_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p31_filted, p31_filted,[imfftw,jmfftw,kmfftw])
      c_p32_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p32_filted, p32_filted,[imfftw,jmfftw,kmfftw])
      c_p33_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p33_filted, p33_filted,[imfftw,jmfftw,kmfftw])
      c_p1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p1_filted, p1_filted,[imfftw,jmfftw,kmfftw])
      c_p2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p2_filted, p2_filted,[imfftw,jmfftw,kmfftw])
      c_p3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_p3_filted, p3_filted,[imfftw,jmfftw,kmfftw])
      c_rho1_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho1_filted, rho1_filted,[imfftw,jmfftw,kmfftw])
      c_rho2_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho2_filted, rho2_filted,[imfftw,jmfftw,kmfftw])
      c_rho3_filted = fftw_alloc_complex(alloc_local)
      call c_f_pointer(c_rho3_filted, rho3_filted,[imfftw,jmfftw,kmfftw])
      !
      !
      allocate(theta_filted(1:im,1:jm,1:km),&
      M11_filted(1:im,1:jm,1:km),M22_filted(1:im,1:jm,1:km),M33_filted(1:im,1:jm,1:km),&
      M12_filted(1:im,1:jm,1:km),M13_filted(1:im,1:jm,1:km),M21_filted(1:im,1:jm,1:km),&
      M23_filted(1:im,1:jm,1:km),M31_filted(1:im,1:jm,1:km),M32_filted(1:im,1:jm,1:km))
      !
      !
      !
      Es = 0.d0
      Ec = 0.d0
      Ts = 0.d0
      Tc = 0.d0
      Ps1 = 0.d0
      Pc1 = 0.d0
      Pc2 = 0.d0
      !
      if(mpirank==0)  print *, "Array allocated and initialized"
      !
      do n=1,num_alpha
        !
        call date_and_time(values=value) 
        !
        if(mpirank==0)  print *, '** Integrate for ',n,'/',num_alpha,',now is ',&
                                value(5), ':', value(6),':',value(7)
        !
        !!! Velocity Favre average and density average
        ! After this bloc, w1_filted is (rho*u1)_filted in spectral space
        do i=1,im
        do j=1,jm
        do k=1,km
          Galpha = exp(-(k1(i,j,k)**2+k2(i,j,k)**2+k3(i,j,k)**2)*sqrtalphas(n)**2/2.d0) ! Filtre scale :sqrtalpha
          w1_filted(i,j,k)  = w1(i,j,k)    *Galpha
          w2_filted(i,j,k)  = w2(i,j,k)    *Galpha
          w3_filted(i,j,k)  = w3(i,j,k)    *Galpha
          rho_filted(i,j,k) = rhocom(i,j,k)*Galpha
          p_filted(i,j,k)   = pcom(i,j,k)  *Galpha
          !
          rho1_filted(i,j,k) = imag*rho_filted(i,j,k)*k1(i,j,k)
          rho2_filted(i,j,k) = imag*rho_filted(i,j,k)*k2(i,j,k)
          rho3_filted(i,j,k) = imag*rho_filted(i,j,k)*k3(i,j,k)
          !
          p11_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k1(i,j,k)
          p21_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k1(i,j,k)
          p31_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k1(i,j,k)
          p12_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k2(i,j,k)
          p22_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k2(i,j,k)
          p32_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k2(i,j,k)
          p13_filted(i,j,k)  = -p_filted(i,j,k)*k1(i,j,k)*k3(i,j,k)
          p23_filted(i,j,k)  = -p_filted(i,j,k)*k2(i,j,k)*k3(i,j,k)
          p33_filted(i,j,k)  = -p_filted(i,j,k)*k3(i,j,k)*k3(i,j,k)
          p1_filted(i,j,k)   = imag*p_filted(i,j,k)*k1(i,j,k)
          p2_filted(i,j,k)   = imag*p_filted(i,j,k)*k2(i,j,k)
          p3_filted(i,j,k)   = imag*p_filted(i,j,k)*k3(i,j,k)
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is (rho*u1)_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(backward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(backward_plan,w3_filted,w3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho_filted,rho_filted)
        call fftw_mpi_execute_dft(backward_plan,p_filted,p_filted)
        !
        call fftw_mpi_execute_dft(backward_plan,p11_filted,p11_filted)
        call fftw_mpi_execute_dft(backward_plan,p21_filted,p21_filted)
        call fftw_mpi_execute_dft(backward_plan,p31_filted,p31_filted)
        call fftw_mpi_execute_dft(backward_plan,p12_filted,p12_filted)
        call fftw_mpi_execute_dft(backward_plan,p22_filted,p22_filted)
        call fftw_mpi_execute_dft(backward_plan,p32_filted,p32_filted)
        call fftw_mpi_execute_dft(backward_plan,p13_filted,p13_filted)
        call fftw_mpi_execute_dft(backward_plan,p23_filted,p23_filted)
        call fftw_mpi_execute_dft(backward_plan,p33_filted,p33_filted)
        call fftw_mpi_execute_dft(backward_plan,p1_filted,p1_filted)
        call fftw_mpi_execute_dft(backward_plan,p2_filted,p2_filted)
        call fftw_mpi_execute_dft(backward_plan,p3_filted,p3_filted)
        call fftw_mpi_execute_dft(backward_plan,rho1_filted,rho1_filted)
        call fftw_mpi_execute_dft(backward_plan,rho2_filted,rho2_filted)
        call fftw_mpi_execute_dft(backward_plan,rho3_filted,rho3_filted)
        !
        ! After this bloc, w1_filted is u1_filted in physical space
        do i=1,im
        do j=1,jm
        do k=1,km
          !
          w1_filted(i,j,k) = w1_filted(i,j,k)/rho_filted(i,j,k)
          w2_filted(i,j,k) = w2_filted(i,j,k)/rho_filted(i,j,k)
          w3_filted(i,j,k) = w3_filted(i,j,k)/rho_filted(i,j,k)
          !
        enddo
        enddo
        enddo
        !
        ! After this bloc, w1_filted is u1_filted in fourier space, A11_filted is A11_filted in fourier space
        call fftw_mpi_execute_dft(forward_plan,w1_filted,w1_filted)
        call fftw_mpi_execute_dft(forward_plan,w2_filted,w2_filted)
        call fftw_mpi_execute_dft(forward_plan,w3_filted,w3_filted)
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          w1_filted(i,j,k)   = w1_filted(i,j,k)/(1.d0*ia*ja*ka)
          w2_filted(i,j,k)   = w2_filted(i,j,k)/(1.d0*ia*ja*ka)
          w3_filted(i,j,k)   = w3_filted(i,j,k)/(1.d0*ia*ja*ka)
          !
          A11_filted(i,j,k) = imag*w1_filted(i,j,k)*k1(i,j,k)
          A21_filted(i,j,k) = imag*w2_filted(i,j,k)*k1(i,j,k)
          A31_filted(i,j,k) = imag*w3_filted(i,j,k)*k1(i,j,k)
          A12_filted(i,j,k) = imag*w1_filted(i,j,k)*k2(i,j,k)
          A22_filted(i,j,k) = imag*w2_filted(i,j,k)*k2(i,j,k)
          A32_filted(i,j,k) = imag*w3_filted(i,j,k)*k2(i,j,k)
          A13_filted(i,j,k) = imag*w1_filted(i,j,k)*k3(i,j,k)
          A23_filted(i,j,k) = imag*w2_filted(i,j,k)*k3(i,j,k)
          A33_filted(i,j,k) = imag*w3_filted(i,j,k)*k3(i,j,k)
          !
        end do
        end do
        end do
        !
        ! After this bloc, A11_filted is A11_filted in physical space
        call fftw_mpi_execute_dft(backward_plan,A11_filted,A11_filted)
        call fftw_mpi_execute_dft(backward_plan,A21_filted,A21_filted)
        call fftw_mpi_execute_dft(backward_plan,A31_filted,A31_filted)
        call fftw_mpi_execute_dft(backward_plan,A12_filted,A12_filted)
        call fftw_mpi_execute_dft(backward_plan,A22_filted,A22_filted)
        call fftw_mpi_execute_dft(backward_plan,A32_filted,A32_filted)
        call fftw_mpi_execute_dft(backward_plan,A13_filted,A13_filted)
        call fftw_mpi_execute_dft(backward_plan,A23_filted,A23_filted)
        call fftw_mpi_execute_dft(backward_plan,A33_filted,A33_filted)
        !
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          theta_filted(i,j,k) = A11_filted(i,j,k)+A22_filted(i,j,k)+A33_filted(i,j,k)
          !
          M11_filted(i,j,k) = A11_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
          M22_filted(i,j,k) = A22_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
          M33_filted(i,j,k) = A33_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)
          !
          M12_filted(i,j,k) = A12_filted(i,j,k)
          M13_filted(i,j,k) = A13_filted(i,j,k)
          M21_filted(i,j,k) = A21_filted(i,j,k)
          M23_filted(i,j,k) = A23_filted(i,j,k)
          M31_filted(i,j,k) = A31_filted(i,j,k)
          M32_filted(i,j,k) = A32_filted(i,j,k)
          !
        end do
        end do
        end do
        !
        !!!! T terms
        !
        do k=1,km
        do j=1,jm
        do i=1,im
          !
          Es = Es + dreal(0.5d0 * rho_filted(i,j,k)*&
          (M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + &
          M13_filted(i,j,k)*M13_filted(i,j,k) + M21_filted(i,j,k)*M21_filted(i,j,k) + &
          M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
          M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + &
          M33_filted(i,j,k)*M33_filted(i,j,k))*dalphas(n))
          !
          Ec = Ec + dreal(1.d0/6.d0 * rho_filted(i,j,k)* theta_filted(i,j,k)*theta_filted(i,j,k) *dalphas(n))
          !
          Ts = Ts + dreal(rho_filted(i,j,k)*&
          (M11_filted(i,j,k)*M11_filted(i,j,k)*M11_filted(i,j,k)+&
          M12_filted(i,j,k)*M21_filted(i,j,k)*M11_filted(i,j,k)+&
          M13_filted(i,j,k)*M31_filted(i,j,k)*M11_filted(i,j,k)+&
          M11_filted(i,j,k)*M12_filted(i,j,k)*M12_filted(i,j,k)+&
          M12_filted(i,j,k)*M22_filted(i,j,k)*M12_filted(i,j,k)+&
          M13_filted(i,j,k)*M32_filted(i,j,k)*M12_filted(i,j,k)+&
          M11_filted(i,j,k)*M13_filted(i,j,k)*M13_filted(i,j,k)+&
          M12_filted(i,j,k)*M23_filted(i,j,k)*M13_filted(i,j,k)+&
          M13_filted(i,j,k)*M33_filted(i,j,k)*M13_filted(i,j,k)+&
          M21_filted(i,j,k)*M11_filted(i,j,k)*M21_filted(i,j,k)+&
          M22_filted(i,j,k)*M21_filted(i,j,k)*M21_filted(i,j,k)+&
          M23_filted(i,j,k)*M31_filted(i,j,k)*M21_filted(i,j,k)+&
          M21_filted(i,j,k)*M12_filted(i,j,k)*M22_filted(i,j,k)+&
          M22_filted(i,j,k)*M22_filted(i,j,k)*M22_filted(i,j,k)+&
          M23_filted(i,j,k)*M32_filted(i,j,k)*M22_filted(i,j,k)+&
          M21_filted(i,j,k)*M13_filted(i,j,k)*M23_filted(i,j,k)+&
          M22_filted(i,j,k)*M23_filted(i,j,k)*M23_filted(i,j,k)+&
          M23_filted(i,j,k)*M33_filted(i,j,k)*M23_filted(i,j,k)+&
          M31_filted(i,j,k)*M11_filted(i,j,k)*M31_filted(i,j,k)+&
          M32_filted(i,j,k)*M21_filted(i,j,k)*M31_filted(i,j,k)+&
          M33_filted(i,j,k)*M31_filted(i,j,k)*M31_filted(i,j,k)+&
          M31_filted(i,j,k)*M12_filted(i,j,k)*M32_filted(i,j,k)+&
          M32_filted(i,j,k)*M22_filted(i,j,k)*M32_filted(i,j,k)+&
          M33_filted(i,j,k)*M32_filted(i,j,k)*M32_filted(i,j,k)+&
          M31_filted(i,j,k)*M13_filted(i,j,k)*M33_filted(i,j,k)+&
          M32_filted(i,j,k)*M23_filted(i,j,k)*M33_filted(i,j,k)+&
          M33_filted(i,j,k)*M33_filted(i,j,k)*M33_filted(i,j,k)+&
          2.d0/3.d0 * &
          (M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k)+ &
          M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
          M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k)) &
          *theta_filted(i,j,k))*dalphas(n))
          !
          Tc = Tc + dreal(rho_filted(i,j,k)*&
          ((M11_filted(i,j,k)*M11_filted(i,j,k) + M12_filted(i,j,k)*M12_filted(i,j,k) + M13_filted(i,j,k)*M13_filted(i,j,k) + &
          M21_filted(i,j,k)*M21_filted(i,j,k) + M22_filted(i,j,k)*M22_filted(i,j,k) + M23_filted(i,j,k)*M23_filted(i,j,k) + &
          M31_filted(i,j,k)*M31_filted(i,j,k) + M32_filted(i,j,k)*M32_filted(i,j,k) + M33_filted(i,j,k)*M33_filted(i,j,k)) &
          *theta_filted(i,j,k) - 1.d0/3.d0*theta_filted(i,j,k)**3)*dalphas(n))
          !
          Ps1 = Ps1 + dreal((&
          (p11_filted(i,j,k)*A11_filted(i,j,k)+p12_filted(i,j,k)*A12_filted(i,j,k)+p13_filted(i,j,k)*A13_filted(i,j,k) &
          +p21_filted(i,j,k)*A21_filted(i,j,k)+p22_filted(i,j,k)*A22_filted(i,j,k)+p23_filted(i,j,k)*A23_filted(i,j,k) &
          +p31_filted(i,j,k)*A31_filted(i,j,k)+p32_filted(i,j,k)*A32_filted(i,j,k)+p33_filted(i,j,k)*A33_filted(i,j,k))& 
          - 1.d0/3.d0 * (p11_filted(i,j,k)+p22_filted(i,j,k)+p33_filted(i,j,k)) * theta_filted(i,j,k)&
          )*dalphas(n))
          !
          Pc1 = Pc1 + dreal(1.d0/3.d0 * (p11_filted(i,j,k)+p22_filted(i,j,k)+p33_filted(i,j,k)) * &
                theta_filted(i,j,k)*dalphas(n))
          !
          Pc2 = Pc2 + dreal(1.d0/3.d0 / rho_filted(i,j,k) * (rho1_filted(i,j,k)*p1_filted(i,j,k) + &
          rho2_filted(i,j,k)*p2_filted(i,j,k) + rho3_filted(i,j,k)*p3_filted(i,j,k))* theta_filted(i,j,k) *dalphas(n))
          !
        enddo
        enddo
        enddo
        !
      enddo
      !
      Es = psum(Es)/(1.d0*ia*ja*ka)
      Ec = psum(Ec)/(1.d0*ia*ja*ka)
      Ts = psum(Ts)/(1.d0*ia*ja*ka)
      Tc = psum(Tc)/(1.d0*ia*ja*ka)
      Ps1 = psum(Ps1)/(1.d0*ia*ja*ka)
      Pc1 = psum(Pc1)/(1.d0*ia*ja*ka)
      Pc2 = psum(Pc2)/(1.d0*ia*ja*ka)
      !
      if(mpirank==0)  print *, 'Job finish'
      !
      if(mpirank==0) then
        if (thefilenumb .ne. 0) then
          outfilename = 'pp/SGS_T_'//stepname//'.dat'
        else
          outfilename = 'pp/SGS_T.dat'
        endif
        
        call listinit(filename=outfilename,handle=hand_a, &
                      firstline='nstep time Es Ec Ts Tc Ps1 Pc1 Pc2')
        call listwrite(hand_a,Es,Ec,Ts,Tc,Ps1,Pc1,Pc2)
        !
        print *, '>>>>', outfilename
      endif
      !
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)
      call fftw_mpi_cleanup()
      call fftw_free(c_w1)
      call fftw_free(c_w2)
      call fftw_free(c_w3)
      call fftw_free(c_rhocom)
      call fftw_free(c_w1_filted)
      call fftw_free(c_w2_filted)
      call fftw_free(c_w3_filted)
      call fftw_free(c_rho_filted)
      call fftw_free(c_A11_filted)
      call fftw_free(c_A12_filted)
      call fftw_free(c_A13_filted)
      call fftw_free(c_A21_filted)
      call fftw_free(c_A22_filted)
      call fftw_free(c_A23_filted)
      call fftw_free(c_A31_filted)
      call fftw_free(c_A32_filted)
      call fftw_free(c_A33_filted)
      call fftw_free(c_p11_filted)
      call fftw_free(c_p12_filted)
      call fftw_free(c_p13_filted)
      call fftw_free(c_p21_filted)
      call fftw_free(c_p22_filted)
      call fftw_free(c_p23_filted)
      call fftw_free(c_p31_filted)
      call fftw_free(c_p32_filted)
      call fftw_free(c_p33_filted)
      call fftw_free(c_p1_filted)
      call fftw_free(c_p2_filted)
      call fftw_free(c_p3_filted)
      call fftw_free(c_rho1_filted)
      call fftw_free(c_rho2_filted)
      call fftw_free(c_rho3_filted)
      call mpistop
      deallocate(theta_filted,M11_filted,M22_filted,M33_filted)
      deallocate(M12_filted,M13_filted,M21_filted)
      deallocate(M23_filted,M31_filted,M32_filted)
      deallocate(k1,k2,k3)
      deallocate(sqrtalphas,dalphas)
      !
    end subroutine SGST3D
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This function is used to read the SGS postprocess 
    ! calculation file at datin/SGSinput
    ! Outputs: 
    !   num_l        : number of ell
    !   num_alpha    : number of maximum alpha for integration
    !   num_alphamin : number of minimum alpha for integration
    !   ratio_max    : maximum ell/lmin 
    !   ratio_min    : minimum ell/lmin
    !   loutput      : logical, output stress hdf5 file or not
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine readSGSinput(num_l,num_alpha,num_alphamin,ratio_max,ratio_min,loutput)
      !
      use parallel,only: bcast,mpirank
      !
      ! local data
      integer, intent(out) :: num_l,num_alpha,num_alphamin
      real(8), intent(out) :: ratio_max,ratio_min
      logical, intent(out) :: loutput
      character(len=64) :: inputfile
      integer :: fh
      !
      inputfile='datin/SGSinput'
      !
      if(mpirank==0) then
        !
        fh=get_unit()
        !
        open(fh,file=trim(inputfile),action='read')
        read(fh,'(//)')
        read(fh,*)num_l,num_alpha,num_alphamin
        read(fh,'(/)')
        read(fh,*)ratio_max,ratio_min
        read(fh,'(/)')
        read(fh,*)loutput
        close(fh)
        print*,' >> ',trim(inputfile),' ... done'
        print*,' >>> Get: Number of l is',num_l,'Number of alpha is',num_alpha,'Minimum number of alpha is',num_alphamin
        print*,' >>> Ratio max:',ratio_max,'Ratio min',ratio_min
        if(loutput)then
          print *, ' >>> Output stress'
        else
          print *, ' >>> No output stress'
        endif
        !
      endif
      !
      call bcast(num_l)
      call bcast(num_alpha)
      call bcast(num_alphamin)
      call bcast(ratio_max)
      call bcast(ratio_min)
      call bcast(loutput)
      !
    end subroutine readSGSinput
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This function is used to generate ell and alpha used for integration 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine SGSscale_allocate(num_l,l_min,ratio_max,ratio_min,l_lim,num_alpha,num_alphamin,&
                                num_alphas,l_sqrtalpha,l_phi,dl_alpha)
      ! 
      integer, intent(in) :: num_l
      real(8), intent(in) :: l_min
      real(8), intent(in) :: ratio_max,ratio_min
      real(8), intent(out) :: l_lim(:)
      integer, intent(in),optional :: num_alpha,num_alphamin
      real(8), intent(out),optional :: l_sqrtalpha(:,:),l_phi(:,:),dl_alpha(:,:)
      integer, intent(out),optional :: num_alphas(:)
      integer :: i,j
      !
      do i=1,num_l
          l_lim(i) = exp(log(ratio_min)+(i-1) * (log(ratio_max)-log(ratio_min)) / (num_l-1)) * l_min
          if(present(num_alpha)) then
              num_alphas(i) = num_alphamin + int(sqrt(l_lim(i)/ratio_max/l_min)*real(num_alpha-num_alphamin))
              !
              do j=1,num_alphas(i)
                  l_sqrtalpha(i,j) = sqrt( exp(log(0.2**2) + &
                  (log((l_lim(i)/l_min)**2) - log(0.2**2))* (j-1) / (num_alphas(i)-1)) )*l_min
                  l_phi(i,j) = sqrt( abs(l_lim(i)**2 - l_sqrtalpha(i,j)**2) )
              enddo
          endif
      enddo
      !
      if(present(num_alpha)) then
          do i=1,num_l
            dl_alpha(i,1) = 0.5d0 * (l_sqrtalpha(i,2)**2)
            !
            do j=2,num_alphas(i)-1
              dl_alpha(i,j) = 0.5d0 * (l_sqrtalpha(i,j+1)**2 -l_sqrtalpha(i,j-1)**2)
            enddo
            !
            dl_alpha(i,num_alphas(i)) = 0.5d0 * (l_sqrtalpha(i,num_alphas(i))**2 -l_sqrtalpha(i,num_alphas(i)-1)**2 )
          enddo
      endif
    end subroutine SGSscale_allocate
    !
    subroutine fft2d(array,plan)
        !
        use, intrinsic :: iso_c_binding
        use commvar, only: im,jm,ia,ja
        type(C_PTR), intent(in) :: plan
        complex(C_DOUBLE_COMPLEX), intent(inout) :: array(:,:,:)
        integer :: i,j
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(plan,array,array)
        array=array/(1.d0*ia*ja)
    end subroutine fft2d
    !
    subroutine fft3d(array,plan)
        !
        use, intrinsic :: iso_c_binding
        use commvar, only: im,jm,km,ia,ja,ka
        !
        complex(C_DOUBLE_COMPLEX), intent(inout) :: array(:,:,:)
        type(C_PTR), intent(in) :: plan
        integer :: i,j,k
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(plan,array,array)
        array=array/(1.d0*ia*ja*ka)
    end subroutine fft3d
    !
    subroutine fft3dtensor(tensor,plan)
        !
        use, intrinsic :: iso_c_binding
        use commvar, only: im,jm,km,ia,ja,ka
        !
        complex(C_DOUBLE_COMPLEX), intent(inout) :: tensor(:,:,:,:,:)
        type(C_PTR), intent(in) :: plan
        integer :: i,j,k
        !
        include 'fftw3-mpi.f03'
        !
        do j=1,3
        do i=1,3
          call fftw_mpi_execute_dft(plan,tensor(:,:,:,i,j),tensor(:,:,:,i,j))
        enddo
        enddo
        tensor=tensor/(1.d0*ia*ja*ka)
    end subroutine fft3dtensor
    !
     subroutine ifft2d(array,plan)
        !
        use, intrinsic :: iso_c_binding
        complex(C_DOUBLE_COMPLEX), intent(inout) :: array(:,:,:)
        type(C_PTR), intent(in) :: plan
        !
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(plan,array,array)
        !
    end subroutine ifft2d
        !
     subroutine ifft3d(array,plan)
        !
        use, intrinsic :: iso_c_binding
        complex(C_DOUBLE_COMPLEX), intent(inout) :: array(:,:,:)
        type(C_PTR), intent(in) :: plan
        include 'fftw3-mpi.f03'
        !
        call fftw_mpi_execute_dft(plan,array,array)
        !
    end subroutine ifft3d
    !
    subroutine ifft3dvector(vector,plan)
        !
        use, intrinsic :: iso_c_binding
        use commvar, only: im,jm,km,ia,ja,ka
        !
        complex(C_DOUBLE_COMPLEX), intent(inout) :: vector(:,:,:,:)
        type(C_PTR), intent(in) :: plan
        integer :: i,j,k
        !
        include 'fftw3-mpi.f03'
        !
        do i=1,3
          call fftw_mpi_execute_dft(plan,vector(:,:,:,i),vector(:,:,:,i))
        enddo
        !
    end subroutine ifft3dvector
    !
    subroutine ifft3dtensor(tensor,plan)
        !
        use, intrinsic :: iso_c_binding
        use commvar, only: im,jm,km,ia,ja,ka
        !
        complex(C_DOUBLE_COMPLEX), intent(inout) :: tensor(:,:,:,:,:)
        type(C_PTR), intent(in) :: plan
        integer :: i,j,k
        !
        include 'fftw3-mpi.f03'
        !
        do j=1,3
        do i=1,3
          call fftw_mpi_execute_dft(plan,tensor(:,:,:,i,j),tensor(:,:,:,i,j))
        enddo
        enddo
        !
    end subroutine ifft3dtensor
    !
    subroutine tensor_multi_3d_rhoABT(term,rho,A,B,sym)
      ! Attention: ik*jk = ij
      !
      implicit none
      !
      complex(8), intent(out) :: term(:,:,:,:,:)
      real(8), intent(in)     :: A(:,:,:,:,:),B(:,:,:,:,:)
      complex(8), intent(in)  :: rho(:,:,:)
      integer,    intent(in), optional :: sym
      !
      integer :: i,j,k
      !
      !
      term = 0.d0
      !
      if (present(sym))then
        !
        do j=1,3
        do i=1,3
        do k=1,3
          term(:,:,:,i,j) = term(:,:,:,i,j) + &
                            dreal(rho(:,:,:)) * (A(:,:,:,i,k) * B(:,:,:,j,k) + B(:,:,:,i,k) * A(:,:,:,j,k))
        enddo
        enddo
        enddo
        !
      else
        do j=1,3
        do i=1,3
        do k=1,3
          term(:,:,:,i,j) = term(:,:,:,i,j) + &
                            dreal(rho(:,:,:)) * A(:,:,:,i,k) * B(:,:,:,j,k)
        enddo
        enddo
        enddo
        !
      endif

    end subroutine tensor_multi_3d_rhoABT
    !
    subroutine tensor_multi_3d_ABT(term,A,B,sym)
      ! Attention: ik*jk = ij
      !
      implicit none
      !
      complex(8), intent(out) :: term(:,:,:,:,:)
      real(8), intent(in)     :: A(:,:,:,:,:),B(:,:,:,:,:)
      logical, intent(in), optional :: sym
      !
      integer :: i,j,k
      !
      !
      term = 0.d0
      !
      if (merge(sym, .false., present(sym))) then
        !
        do j=1,3
        do i=1,3
        do k=1,3
          term(:,:,:,i,j) = term(:,:,:,i,j) + &
                            (A(:,:,:,i,k) * B(:,:,:,j,k) + B(:,:,:,i,k) * A(:,:,:,j,k))
        enddo
        enddo
        enddo
        !
      else
        do j=1,3
        do i=1,3
        do k=1,3
          term(:,:,:,i,j) = term(:,:,:,i,j) + A(:,:,:,i,k) * B(:,:,:,j,k)
        enddo
        enddo
        enddo
        !
      endif

    end subroutine tensor_multi_3d_ABT
    !
    subroutine scalar_gradient(A, w, k1, k2, k3)

      complex(8), intent(out) :: A(:,:,:,:)
      complex(8), intent(in)  :: w(:,:,:)

      real(8),    intent(in)  :: k1(:,:,:)
      real(8),    intent(in)  :: k2(:,:,:)
      real(8),    intent(in)  :: k3(:,:,:)

      complex(8), parameter :: imag = CMPLX(0.d0,1.d0,8)

      ! d()/dx
      A(:,:,:,1) = imag*w*k1
      A(:,:,:,2) = imag*w*k2
      A(:,:,:,3) = imag*w*k3

    end subroutine scalar_gradient
    !
    subroutine vector_gradient(A, w1, w2, w3, &
                                 k1, k2, k3)

      complex(8), intent(out) :: A(:,:,:,:,:)
      complex(8), intent(in)  :: w1(:,:,:)
      complex(8), intent(in)  :: w2(:,:,:)
      complex(8), intent(in)  :: w3(:,:,:)

      real(8),    intent(in)  :: k1(:,:,:)
      real(8),    intent(in)  :: k2(:,:,:)
      real(8),    intent(in)  :: k3(:,:,:)

      complex(8), parameter :: imag = CMPLX(0.d0,1.d0,8)

      ! d()/dx
      A(:,:,:,1,1) = imag*w1*k1
      A(:,:,:,2,1) = imag*w2*k1
      A(:,:,:,3,1) = imag*w3*k1
      A(:,:,:,1,2) = imag*w1*k2
      A(:,:,:,2,2) = imag*w2*k2
      A(:,:,:,3,2) = imag*w3*k2
      A(:,:,:,1,3) = imag*w1*k3
      A(:,:,:,2,3) = imag*w2*k3
      A(:,:,:,3,3) = imag*w3*k3

    end subroutine vector_gradient
    !
    subroutine tensor_scale_3d_real(T,scale)
      implicit none

      complex(8), intent(inout), dimension(:,:,:,:,:)  :: T
      real(8),    intent(in)   , dimension(:,:,:)      :: scale

      integer :: i,j

      do j=1,3
      do i=1,3
        T(:,:,:,i,j) = T(:,:,:,i,j) * scale
      enddo
      enddo

    end subroutine tensor_scale_3d_real
    !
    subroutine tensor_scale_3d_complex(T,scale)
      implicit none

      complex(8), intent(inout), dimension(:,:,:,:,:)  :: T
      complex(8), intent(in)   , dimension(:,:,:)      :: scale

      integer :: i,j

      do j=1,3
      do i=1,3
        T(:,:,:,i,j) = T(:,:,:,i,j) * scale
      enddo
      enddo

    end subroutine tensor_scale_3d_complex
    !
    subroutine tensor_multi_3d_scalar_real(term,S,f1,rev)
      implicit none

      complex(8), intent(out), dimension (:,:,:,:,:) :: term
      real(8),    intent(in) , dimension (:,:,:,:,:)  :: S
      real(8),    intent(in) , dimension (:,:,:)  :: f1
      logical,    intent(in) , optional :: rev
      integer :: i,j
      if (merge(rev, .false., present(rev)))then
        do j=1,3
        do i=1,3
          term(:,:,:,i,j) = S(:,:,:,j,i)*f1
        end do
        end do
      else
        do j=1,3
        do i=1,3
          term(:,:,:,i,j) = S(:,:,:,i,j)*f1
        end do
        end do
      endif

    end subroutine tensor_multi_3d_scalar_real
    !
    !
    subroutine tensor_multi_3d_scalar_auto_real(term,f1)
      implicit none

      complex(8), intent(inout), dimension (:,:,:,:,:) :: term
      real(8),    intent(in)   , dimension (:,:,:)     :: f1

      integer :: i,j

      do j=1,3
      do i=1,3
        term(:,:,:,i,j) = term(:,:,:,i,j)*f1
      end do
      end do

    end subroutine tensor_multi_3d_scalar_auto_real
    !
    subroutine tensor_multi_3d_scalar_auto_complex(term,f1)
      implicit none

      complex(8), intent(inout), dimension (:,:,:,:,:) :: term
      complex(8), intent(in)   , dimension (:,:,:)     :: f1

      integer :: i,j

      do j=1,3
      do i=1,3
        term(:,:,:,i,j) = term(:,:,:,i,j)*f1
      end do
      end do

    end subroutine tensor_multi_3d_scalar_auto_complex
    !
    subroutine tensor_multi_3d_scalar_complex(term,S,f1)
      implicit none

      complex(8), intent(out), dimension (:,:,:,:,:) :: term
      real(8),    intent(in),  dimension (:,:,:,:,:)  :: S
      complex(8), intent(in),  dimension (:,:,:)  :: f1

      integer :: i,j

      do j=1,3
      do i=1,3
        term(:,:,:,i,j) = S(:,:,:,i,j)*f1
      end do
      end do

    end subroutine tensor_multi_3d_scalar_complex
end module udf_pp_SGS