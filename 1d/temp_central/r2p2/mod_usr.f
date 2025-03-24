module mod_usr
  use mod_hd

      implicit none

      double precision :: heatunit, bQ0

contains

  subroutine usr_init()

      usr_init_one_grid => rm1d_init_one_grid

      !usr_create_particles => generate_particles

      call set_coordinate_system("Cartesian")


      unit_length         = 1d9
      unit_numberdensity  = 1d9
      unit_temperature    = 1d6

      usr_source          => special_source

      
      call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
    double precision :: rho0, p0, v0

      heatunit=unit_pressure/unit_time !3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s
      bQ0=6.0d-5/heatunit       !background heating power density !MLB lo he dividido entre 2, 1d-4
      
      rho0 = 2                  ! 0.5 quizas mejor  quiero 3e-16, con 1 tengo 2.3e-15
      p0 = 2                    ! 0.5 quizas mejor quiero 7e-2, con 1 tengo 3e-1
      v0 = 0                      ! 0.2 quizas mejor
      
      w(ixmin1:ixmax1,rho_)   = rho0
      w(ixmin1:ixmax1,mom(1)) = v0
      w(ixmin1:ixmax1,p_)     = p0
    
    call hd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)

  end subroutine rm1d_init_one_grid

      
 subroutine special_source(qdt,ixImin1,ixImax1,ixOmin1,ixOmax1,iwmin,iwmax,qtC,wCT,qt,w,x)
      integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1, iwmin,iwmax
      double precision, intent(in) :: qdt, qtC, qt
      double precision, intent(in) :: x(ixImin1:ixImax1,1:ndim),wCT(ixImin1:ixImax1,1:nw)
      double precision, intent(inout) :: w(ixImin1:ixImax1,1:nw)
      double precision :: lQgrid(ixImin1:ixImax1)

! add local heating lQ
      
      call getlQ(lQgrid,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,wCT,x)
      w(ixOmin1:ixOmax1,e_)=w(ixOmin1:ixOmax1,e_)+qdt*lQgrid(ixOmin1:ixOmax1)
      end subroutine special_source


      

   subroutine getlQ(lQgrid,ixImin1,ixImax1,ixOmin1,ixOmax1,qt,w,x)
  ! calculate localized heating lQ
  !nueva subrutina MLB
    integer, intent(in) :: ixImin1,ixImax1, ixOmin1,ixOmax1
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,&
       1:ndim), w(ixImin1:ixImax1,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1),&
        A(ixImin1:ixImax1), lQd, lQ0, hp, lQt,Atop,Alow
    double precision :: auxradius1(ixImin1:ixImax1),xcenter1,auxm1
       !ycenter1,y11(ixImin1:ixImax1),y12(ixImin1:ixImax1)
    double precision :: radiuscentral(ixImin1:ixImax1)
    double precision :: auxradius2(ixImin1:ixImax1),xcenter2,auxm2
       !ycenter2,,y21(ixImin1:ixImax1,ixImin2:ixImax2),y22(ixImin1:ixImax1,&
       !ixImin2:ixImax2)
    double precision :: lambdamin,lambdamax,lambdamax2,lambda1(ixImin1:ixImax1),lambda11(ixImin1:ixImax1)
    double precision :: lambda2(ixImin1:ixImax1),&
       lambda22(ixImin1:ixImax1),xprimaL(ixImin1:ixImax1),zprimaL(ixImin1:ixImax1)
    double precision :: xprimaR(ixImin1:ixImax1), zprimaR(ixImin1:ixImax1)
    double precision :: thetaangle01,thetaangle1(ixImin1:ixImax1),thetaangle02,thetaangle2(ixImin1:ixImax1)
    double precision :: lQgridleft(ixImin1:ixImax1),lQgridright(ixImin1:ixImax1)
    double precision :: recta1(ixImin1:ixImax1),recta2(ixImin1:ixImax1)
    double precision :: sigma_heating,sigma_heating_impulsive,y1,lam1
    double precision :: impulsiveheatingfactor,impulsiveheatingduration,&
       timeimpulse
       double precision :: impulsiveheatingfactoramplitude, A0
      
       double precision :: heating_half, heating_start, heating_end

    lQ0= 1.d-1/heatunit ! 8.d-3/heatunit !OSP!5.d0*1.d-2/heatunit!MLB 1.d-2/heatunit
    lQt=10.d0!5.d2/unit_time !MLB 5.d2/unit_time
    lQd=0.2d0 !MLB 0.2d0
    hp=0.22d0
    lQgrid=0.d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! lambda
    lambdamin=0.1d0
    lambdamax=0.320d0
    lambdamax2=500.d0 !es esencialmente infinita


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calentamiento tipo Zhou et al. 2023
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      sigma_heating=0.3d0
      sigma_heating_impulsive=0.15d0
      y1=0.22d0
      lam1=1.0d0!0.2d0!0.32d0
      xcenter1=0.0!-5.2288454d0 !LEFT
      !xcenter2= 8.5744060!5.2288454d0  !RIGHT
      A0=-132.0d0 !vector potential at the center of the field lines bundle
      !
      ! !impulsive heating
      impulsiveheatingfactor=1.d0
      impulsiveheatingduration=1.d0 ! every unit is 85 secs
      timeimpulse=3.001d0
      impulsiveheatingfactoramplitude=0.d0
      if(qt>timeimpulse .and. qt<timeimpulse + impulsiveheatingduration) then
        if(qt<timeimpulse+impulsiveheatingduration/4.d0) then
          impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude* &
             (4.d0*(qt-timeimpulse)/impulsiveheatingduration)
        endif
        if(qt>=timeimpulse+impulsiveheatingduration/4.d0 .and. qt < &
           timeimpulse+impulsiveheatingduration*3.d0/4.d0) then
         impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude
        endif
        if(qt>=timeimpulse+impulsiveheatingduration*3.d0/4.d0 .and. &
           qt<timeimpulse+impulsiveheatingduration) then
         impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude* &
            (-4.d0/impulsiveheatingduration*(qt-timeimpulse-&
            3.d0/4.d0*impulsiveheatingduration)+1.d0)
        endif
      
      endif

      !!CUIDADO HAY UNA ASIMETRÍA EN EL HEATING DEL 0.75 POR UNO.
      lQgrid(ixOmin1:ixOmax1)=lQ0*(impulsiveheatingfactor*dexp(&
      -(x(ixOmin1:ixOmax1,1)-xcenter1)**2/sigma_heating**2))

! !MLB
!OSP dimensional time in s for beginning and end of heating
      heating_start = 100d0 / unit_time
      heating_end = 200d0 / unit_time
      heating_half = (heating_end + heating_start)/2
      
      if(qt<heating_start) then
         lQgrid(ixOmin1:ixOmax1)=0.d0
      endif
      if(qt>heating_start .and. qt< heating_half) then
         lQgrid(ixOmin1:ixOmax1)=(qt-heating_start)/(heating_half-heating_start)*lQgrid(ixOmin1:ixOmax1)
      endif

      if(qt>heating_half .and. qt<heating_end) then
         lQgrid(ixOmin1:ixOmax1)=(heating_end-qt)/(heating_end-heating_half)*lQgrid(ixOmin1:ixOmax1)
      endif
      if(qt>heating_end) then
         lQgrid(ixOmin1:ixOmax1)=0.d0
      endif      
  end subroutine getlQ

      
end module mod_usr
