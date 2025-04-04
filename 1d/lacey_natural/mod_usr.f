module mod_usr
  use mod_hd

      implicit none

      character(32) :: init_cond_file

contains

  subroutine usr_init()

      usr_init_one_grid => rm1d_init_one_grid

      call set_coordinate_system("Cartesian")


      unit_length         = 335312549104.5808 !cm
      unit_numberdensity  = 170818214.24 
      unit_temperature    = 1.143d6

      !usr_source          => special_source

      call params_read(par_files)
      call hd_activate()

  end subroutine usr_init

    !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

      namelist /init_cond_list/ init_cond_file

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, init_cond_list, end=113)
113    close(unitpar)
      end do
      
  end subroutine params_read


  function read_initial_conditions(r_in,index) result(var)
      integer, intent(in) :: index
      double precision, intent(in) :: r_in
      double precision :: var

      double precision :: w(1:4), w_mo(1:4), w_po(1:4)
      integer :: ll, ioerr
      integer :: total_lines
      logical :: first = .true.  

      total_lines = get_number_lines(init_cond_file)
      
      w(:) = 0.d0
      
      if (mype==0 .and. first .eqv. .true.) then
         print*, 'Reading initial condition data from file:   ', trim(init_cond_file)
         first = .false.
      endif
      open(unit=1, file=trim(init_cond_file), status="old")!,&
      ioerr = 1
      if (ioerr/=0) then
         read(1,*) w            !> first line of data
         do ll = 1,total_lines-1
            w_mo = w
            read(1,*) w
            if (w(1) .gt. r_in) then
               w_po = w
               goto 8765
            endif
            w_po = w
         enddo
      endif
      
8765 CLOSE(1)
      var = w_mo(index)
      
  end function read_initial_conditions


      
  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
      use mod_init_datafromfile

      integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
      double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
      double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
      
      integer :: ii
      
      do ii = ixgmin1, ixgmax1
         w(ii,rho_) = read_initial_conditions(x(ii,1),2)
         w(ii,mom(1)) = read_initial_conditions(x(ii,1),3)
         w(ii,p_) = read_initial_conditions(x(ii,1),4)
         !print*, x(ii,1),  w(ii,rho_), w(ii,mom(1)), w(ii,p_)
      enddo
    
    call hd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
  end subroutine rm1d_init_one_grid
!
!      
! subroutine special_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,w,x)
!      integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax
!      double precision, intent(in) :: qdt, qtC, qt
!      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim),wCT(ixImin1:ixImax1,1:nw)
!      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
!      double precision :: lQgrid(ixImin1:ixImax1)
!
!! add local heating lQ
!      
!      call getlQ(lQgrid,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,wCT,x)
!      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)+qdt*lQgrid(ixOmin1:ixOmax1)
!      end subroutine special_source
!
!
!      
!
!   subroutine getlQ(lQgrid,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x)
!  ! calculate localized heating lQ
!  !nueva subrutina MLB
!    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
!    double precision, intent(in) :: qt, x(ixImin1:ixImax1,&
!       1:ndim), w(ixImin1:ixImax1,1:nw)
!
!    double precision :: lQgrid(ixImin1:ixImax1),&
!        A(ixImin1:ixImax1), lQd, lQ0, hp, lQt,Atop,Alow
!    double precision :: auxradius1(ixImin1:ixImax1),xcenter1,auxm1
!       !ycenter1,y11(ixImin1:ixImax1),y12(ixImin1:ixImax1)
!    double precision :: radiuscentral(ixImin1:ixImax1)
!    double precision :: auxradius2(ixImin1:ixImax1),xcenter2,auxm2
!       !ycenter2,,y21(ixImin1:ixImax1,ixImin2:ixImax2),y22(ixImin1:ixImax1,&
!       !ixImin2:ixImax2)
!    double precision :: lambdamin,lambdamax,lambdamax2,lambda1(ixImin1:ixImax1),lambda11(ixImin1:ixImax1)
!    double precision :: lambda2(ixImin1:ixImax1),&
!       lambda22(ixImin1:ixImax1),xprimaL(ixImin1:ixImax1),zprimaL(ixImin1:ixImax1)
!    double precision :: xprimaR(ixImin1:ixImax1), zprimaR(ixImin1:ixImax1)
!    double precision :: thetaangle01,thetaangle1(ixImin1:ixImax1),thetaangle02,thetaangle2(ixImin1:ixImax1)
!    double precision :: lQgridleft(ixImin1:ixImax1),lQgridright(ixImin1:ixImax1)
!    double precision :: recta1(ixImin1:ixImax1),recta2(ixImin1:ixImax1)
!    double precision :: sigma_heating,sigma_heating_impulsive,y1,lam1
!    double precision :: impulsiveheatingfactor,impulsiveheatingduration,&
!       timeimpulse
!       double precision :: impulsiveheatingfactoramplitude, A0
!      
!       double precision :: heating_half, heating_start, heating_end
!
!    lQ0= 5.d-1/heatunit ! 8.d-3/heatunit !OSP!5.d0*1.d-2/heatunit!MLB 1.d-2/heatunit
!    lQt=10.d0!5.d2/unit_time !MLB 5.d2/unit_time
!    lQd=0.2d0 !MLB 0.2d0
!    hp=0.22d0
!    lQgrid=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! lambda
!    lambdamin=0.1d0
!    lambdamax=0.320d0
!    lambdamax2=500.d0 !es esencialmente infinita
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Calentamiento tipo Zhou et al. 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      sigma_heating=0.3d0
!      sigma_heating_impulsive=0.15d0
!      y1=0.22d0
!      lam1=1.0d0!0.2d0!0.32d0
!      xcenter1=0.0!-5.2288454d0 !LEFT
!      !xcenter2= 8.5744060!5.2288454d0  !RIGHT
!      A0=-132.0d0 !vector potential at the center of the field lines bundle
!      !
!      ! !impulsive heating
!      impulsiveheatingfactor=1.d0
!      impulsiveheatingduration=1.d0 ! every unit is 85 secs
!      timeimpulse=3.001d0
!      impulsiveheatingfactoramplitude=0.d0
!      if(qt>timeimpulse .and. qt<timeimpulse + impulsiveheatingduration) then
!        if(qt<timeimpulse+impulsiveheatingduration/4.d0) then
!          impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude* &
!             (4.d0*(qt-timeimpulse)/impulsiveheatingduration)
!        endif
!        if(qt>=timeimpulse+impulsiveheatingduration/4.d0 .and. qt < &
!           timeimpulse+impulsiveheatingduration*3.d0/4.d0) then
!         impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude
!        endif
!        if(qt>=timeimpulse+impulsiveheatingduration*3.d0/4.d0 .and. &
!           qt<timeimpulse+impulsiveheatingduration) then
!         impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude* &
!            (-4.d0/impulsiveheatingduration*(qt-timeimpulse-&
!            3.d0/4.d0*impulsiveheatingduration)+1.d0)
!        endif
!      
!      endif
!
!      !!CUIDADO HAY UNA ASIMETRÍA EN EL HEATING DEL 0.75 POR UNO.
!      lQgrid(ixOmin1:ixOmax1)=lQ0*(impulsiveheatingfactor*dexp(&
!      -(x(ixOmin1:ixOmax1,1)-xcenter1)**2/sigma_heating**2))
!
!! !MLB
!!OSP dimensional time in s for beginning and end of heating
!      heating_start = 100d0 / unit_time
!      heating_end = 200d0 / unit_time
!      heating_half = (heating_end + heating_start)/2
!      
!      if(qt<heating_start) then
!         lQgrid(ixOmin1:ixOmax1)=0.d0
!      endif
!      if(qt>heating_start .and. qt< heating_half) then
!         lQgrid(ixOmin1:ixOmax1)=(qt-heating_start)/(heating_half-heating_start)*lQgrid(ixOmin1:ixOmax1)
!      endif
!
!      if(qt>heating_half .and. qt<heating_end) then
!         lQgrid(ixOmin1:ixOmax1)=(heating_end-qt)/(heating_end-heating_half)*lQgrid(ixOmin1:ixOmax1)
!      endif
!      if(qt>heating_end) then
!         lQgrid(ixOmin1:ixOmax1)=0.d0
!      endif      
!  end subroutine getlQ

  function get_number_lines(filename)
    integer, parameter :: MAX_LINE_LENGTH=1024
    character(len=MAX_LINE_LENGTH) line
    character(len=*), intent(in) :: filename
    integer, parameter :: lunit=338
    integer :: get_number_lines

    get_number_lines = 0
    OPEN (lunit, file = filename)
    DO
        READ (lunit,*, END=10) line
        if (line(1:1)=='#') cycle
        get_number_lines = get_number_lines + 1
    END DO
    10 CLOSE (lunit)
  end function get_number_lines

      
end module mod_usr
