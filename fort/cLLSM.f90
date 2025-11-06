
program call_cLLSM

  implicit none

  complex*16, parameter         :: Jcomp = (0.d0,1.d0)

  integer                       :: ng
  integer                       :: l_max, nk, npt
  integer                       :: i, j, k, l

  logical                       :: cusp

  character(len = 40)           :: name_file
  character(len = 5 )           :: charI

  double precision              :: expo_re, expo_im, rr
  complex*16                    :: tmp

  double precision, allocatable :: rgrid(:)
  complex*16,       allocatable :: expo(:,:), coeff(:,:,:)
  complex*16,       allocatable :: F(:,:,:), Fit(:,:,:), er(:,:,:)

  read *, ng 
  read *, npt
  read *, l_max 
  read *, nk
  read *, i
  cusp = .false.
  if(i .eq. 1) cusp = .true. ! add r^{l+1} is cGTOs
  print *, " add r^{l+1} ?", cusp

  allocate(expo(ng,0:l_max))
  do l = 0, l_max
    write(charI, '(I5)') l
    write(name_file, '("../cGTOs/expo/set_3/expo_l", A, ".txt")') trim(adjustl(charI))
    open(unit = 11, file=name_file, action="read")
      do i = 1, ng
        read(11,*) expo_re, expo_im
        expo(i,l) = expo_re + Jcomp*expo_im
      enddo
    close(11)
  enddo

  allocate(rgrid(npt))
  open(unit = 11, file="my_grid.dat", action="read")
    do i = 1, npt
      read(11,*) rgrid(i)
    enddo
  close(11)

  allocate(F(npt,nk,0:l_max))
  do l = 0, l_max
    write(charI, '(I5)') l
    write(name_file, '("radfunc_l", A, ".dat")') trim(adjustl(charI))
    open(unit = 11, file=name_file, action="read")
      do k = 1, nk
        do i = 1, npt
          read(11,*) F(i,k,l)
        enddo
      enddo
    close(11)
  enddo


  allocate(coeff(ng,nk,0:l_max))
  ! u_l = r^{nl+1} sum_i c_i exp(-alpha_i r^2)
  do l = 0, l_max
    call cLLSM(nk, ng, npt, l, rgrid, expo(:,l), F(:,:,l), coeff(:,:,l), cusp)
  enddo


  allocate(Fit(npt,nk,0:l_max), er(npt,nk,0:l_max))
  do l = 0, l_max
    do k = 1, nk
      do j = 1, npt
        rr  = rgrid(j)
        tmp = (0.d0,0.d0)
        do i = 1, ng
          tmp = tmp + coeff(i,k,l) * zexp(-expo(i,l)* rr*rr)
        enddo
        Fit(j,k,l) = tmp * rr**dble(l+1)
        er (j,k,l) = dabs( real(F(j,k,l)) - real(Fit(j,k,l)) ) + &
                               Jcomp*dabs( aimag(Fit(j,k,l)) )
      enddo
    enddo
  enddo


  do l = 0, l_max
    write(charI, '(I5)') l

    write(name_file, '("coeff_l", A, ".dat")') trim(adjustl(charI))
    open(unit=100, file=name_file, action="write")
      do k = 1, nk
        do i = 1, ng
          write(100,'(2(d35.28,1x))') coeff(i,k,l)
        enddo
      enddo
    close(100)

    write(name_file, '("F_l", A, ".dat")') trim(adjustl(charI))
    open(unit=100, file=name_file, action="write")
      do i = 1, npt
        write(100,'(21(e15.8,1x))') rgrid(i), (F(i,1:nk,l))
      end do
    close(100)

    write(name_file, '("Fit_l", A, ".dat")') trim(adjustl(charI))
    open(unit=100, file=name_file, action="write")
      do i = 1, npt
        write(100,'(21(e15.8,1x))') rgrid(i), (Fit(i,1:nk,l))
      end do
    close(100)

    write(name_file, '("Err_l", A, ".dat")') trim(adjustl(charI))
    open(unit=100, file=name_file, action="write")
      do i = 1, npt
        write(100,'(21(e15.8,1x))') rgrid(i), (er(i,1:nk,l))
      end do
    close(100)

  enddo

  deallocate(rgrid, expo, coeff, F, Fit, er)

end program call_cLLSM

! ---

subroutine cLLSM(nk, ng, npt, nl, rgrid, expo, F, coeff, cusp)

  implicit none

  logical,          intent(in)  :: cusp
  integer,          intent(in)  :: nk, ng, npt, nl
  double precision, intent(in)  :: rgrid(npt)
  complex*16,       intent(in)  :: F(npt,nk), expo(ng)

  complex*16, intent(out)       :: coeff(ng,nk)

  integer                       :: i, j, k

  ! used by Lapack
  integer                       :: rank, info, nlvl, LworK, Lrwork, Liwork
  double precision              :: rcond = -1.d0, nlvlr
  integer, allocatable          :: Iwork(:)
  double precision, allocatable :: S(:), RWork(:)
  complex*16, allocatable       :: AA(:,:), Work(:), Fout(:,:)

  !-----------------------------------------------------------------------------
  ! set Lapack parameters
  nlvlr  = dlog(dble(ng)/26.d0) / dlog(2.d0)
  nlvl   = int(nlvlr) + 1
  Lwork  = (2+ng)*npt
  Lrwork = 60*npt + 8*npt*nlvl + 75*nk + 676
  Liwork = ng*(3*nlvl+11)
  !-----------------------------------------------------------------------------

  ! functions to Fit
  allocate( Fout(npt,nk) )
  Fout(1:npt,1:nk) = F(1:npt,1:nk)
  ! On exit from ZGELSD, Fout is overwritten by optimized linear coefficients

  ! minimize |AA X - B | where AA is given by
  allocate( AA(npt,ng) )
  if(cusp) then
    do i = 1, ng
      do j = 1, npt
        AA(j,i) = rgrid(j)**dble(nl+1) * zexp(-expo(i)*rgrid(j)*rgrid(j))
      enddo
    enddo
  else
    do i = 1, ng
      do j = 1, npt
        AA(j,i) = zexp(-expo(i)*rgrid(j)*rgrid(j))
      enddo
    enddo
  endif

  ! call lapack 
  allocate( S(ng), Work(LworK), Rwork(Lrwork), Iwork(Liwork) )
  call ZGELSD(npt, ng, nk, AA, npt, Fout, npt, S, rcond, rank, Work, LworK, &
       RWork, Iwork, info)
  deallocate(AA, S, Work, Rwork, Iwork)

  if( info .ne. 0 ) then
    write(*,*) 'ZGELSD in cLLSM subroutine: info != 0'
    return
  end if

  ! optimised linear coefficients
  coeff(1:ng,1:nk) = Fout(1:ng,1:nk)
  deallocate(Fout)

end subroutine cLLSM

! ---

