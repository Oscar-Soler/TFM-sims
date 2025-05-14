module mod_usr
  use mod_mhd
  implicit none
  double precision, allocatable :: pbc(:),rbc(:)
  double precision :: usr_grav,vmax,La
  double precision :: heatunit,gzone,B0,theta,SRadius,kx,ly,bQ0,dya
  double precision :: md0, xc1, xc2, xc3, hb, hb2, hb3
  double precision :: md1, md2, md3, adipo, factorcentralarcade
  double precision, allocatable :: pa(:),ra(:)

  integer, parameter :: jmax=8000

  character(2) :: name_bt='bt'
  integer, dimension(3) :: bt_
  integer :: i

contains

  subroutine usr_init()
    call set_coordinate_system("Cartesian_2.5D")

    unit_length        = 1.d9                                         ! cm
    unit_temperature   = 1.d6                                         ! K
    unit_numberdensity = 1.d9                                         ! cm^-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    usr_aux_output      => specialvar_output
    usr_add_aux_names   => specialvarnames_output
    usr_init_vector_potential=>initvecpot_usr
    
    usr_modify_output => set_output_vars

    ! usr_set_wLR         => boundary_wLR
    !usr_refine_threshold=> specialthreshold
    !usr_set_electric_field => boundary_electric_field
    
    call mhd_activate()
    do i=1, ndir
      bt_(i) =  var_set_extravar(name_bt, name_bt, i)
    end do

  end subroutine usr_init

  subroutine initglobaldata_usr()
    heatunit=unit_pressure/unit_time          ! 3.697693390805347E-003 erg*cm^-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=6.0d-5/heatunit ! background heating power density !MLB lo he dividido entre 2, 1d-4
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary !MLB
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) ! cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 ! the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade
    kx=dpi/(xprobmax1-xprobmin1)
    ly=kx*dcos(theta)
    SRadius=69.61d0 ! Solar radius
    vmax=0.1d0!7.d5/unit_velocity ! maximal driven velocity
    La=1.5d9/unit_length
    ! hydrostatic vertical stratification of density, temperature, pressure

    !!!!!!!!Magnetic field
    ! B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    hb = 49.09696d8/unit_length !31.7d8/unit_length
    hb2= hb!15.0d8/unit_length
    hb3= 25.d8/unit_length
    ! md0=-B0*hb**2
    ! md1= 0.6*md0
    ! md2= 0.2d0*md0*1.d0 !0.88
    ! md3=-0.4*md0
    xc1= 60.d8/unit_length ! 0.d0*30.d8/unit_length
    xc2= -60.d8/unit_length !-0.d0*30.d8/unit_length
    xc3= 0.d0 !-0.d0*32.d8/unit_length
    adipo=0.1d0
    factorcentralarcade=0.23d0


    call inithdstatic

  end subroutine initglobaldata_usr

  subroutine inithdstatic
    use mod_solar_atmosphere
    ! initialize the table in a vertical line through the global domain
    integer :: j,na,ibc
    double precision, allocatable :: Ta(:),gg(:),ya(:)
    double precision :: rpho,Ttop,Tpho,wtra,res,rhob,pb,htra,Ttr,Fc,invT,kappa
    double precision :: rhohc,hc
    logical :: simple_temperature_curve

    simple_temperature_curve=.true.

    allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))

    if(simple_temperature_curve) then
      rpho=1.151d15/unit_numberdensity ! number density at the bottom of height table
      Tpho=8.d3/unit_temperature ! temperature of chromosphere
      Ttop=1.5d6/unit_temperature ! estimated temperature in the top
      htra=0.2d0 ! height of initial transition region
      wtra=0.02d0 ! width of initial transition region 
      Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
      Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
      kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      do j=1,jmax
         ya(j)=(dble(j)-0.5d0)*dya-gzone
         if(ya(j)>htra) then
           Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
         else
           Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
         endif
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      ! rpho=1.151d15/unit_numberdensity ! number density at the bottom of height table !MLB paso de 15 a 12
      ! Tpho=8.d3/unit_temperature ! temperature of chromosphere !MLB de 8 a 40
      ! Ttop=3.0d6/unit_temperature ! estimated temperature in the top !MLB de 1 M a 0.8 0.8d6/unit_temperature
      ! htra=0.2d0 ! height of initial transition region !MLB  de 0.2 a 0.1
      ! wtra=0.02d0 ! width of initial transition region
      ! Ttr=0.8d6/unit_temperature ! lowest temperature of upper profile !MLB de 1.6d5 a 1.d6 1.6d5/unit_temperature
      ! Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
      ! kappa=8.d-7*unit_temperature**3.5d0/unit_length/unit_density/unit_velocity**3
      ! do j=1,jmax
      !    ya(j)=(dble(j)-0.5d0)*dya-gzone
      !    if(ya(j)>htra+0.1d0*wtra) then
      !      ! Ta(j)=(3.5d0*Fc/kappa*(ya(j)-htra)+Ttr**3.5d0)**(2.d0/7.d0)
      !      ! Ta(j)=Ttr*((htra/ya(j)))**(2.d0/7.d0) !MLB 1 NO VALIDA
      !      ! Ta(j)=((htra/ya(j))*(ttr**7.d0/2.d0-Ttop**7.d0/2.d0)+Ttop**7.d0/2.d0)**(2.d0/7.d0) !MLB 2 NO VALIDA
      !      Ta(j)=(((Ttop**3.5d0-Ttr**3.5d0)/&
      !      (xprobmax2-htra-2.d0*wtra))*(ya(j)-htra-2.d0*wtra)+Ttr**3.5d0)**(2.d0/7.d0) !MLB 3 ojo con la altura máxima
      !    else
      !      ! Ta(j)=Tpho+0.5d0*(Ttop-Tpho)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
      !      Ta(j)=Tpho+0.5d0*(Ttr-Tpho)*(tanh((ya(j)-htra)/wtra)+1.d0)
      !    endif
      !    gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      ! enddo
      ! solution of hydrostatic equation
      ra(1)=rpho
      pa(1)=rpho*Tpho
      invT=gg(1)/Ta(1)
      invT=0.d0
      do j=2,jmax
         invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
         pa(j)=pa(1)*dexp(invT*dya)
         ra(j)=pa(j)/Ta(j)
      end do
    else
      do j=1,jmax
         ! get height table
         ya(j)=(dble(j)-0.5d0)*dya-gzone
         ! get gravity table
         gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
      enddo
      ! a coronal height at 10 Mm
      hc=1.d9/unit_length
      ! the number density at the coronal height
      rhohc=6.2d8/unit_numberdensity
      ! get density and pressure table in hydrostatic state with a preset temperature table
      ! call get_atm_para(ya,ra,pa,gg,jmax,'AL-C7',hc,rhohc)
      call get_atm_para(ya,ra,pa,gg,jmax,'VAL-C',hc,rhohc) !MLB voy a probar otros perfiles de temperatura
    end if
    deallocate(ya,gg,Ta)
    !! initialized rho and p in the fixed bottom boundary
    na=floor(gzone/dya+0.5d0)
    res=gzone-(dble(na)-0.5d0)*dya
    rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
    pb=pa(na)+res/dya*(pa(na+1)-pa(na))
    allocate(rbc(nghostcells))
    allocate(pbc(nghostcells))
    do ibc=nghostcells,1,-1
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  

  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir) !MLB
    double precision :: res
    integer :: ix^D,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating prominence perturbation'
      endif
      first=.false.
    endif
    {do ix^DB=ixOmin^DB,ixOmax^DB\}
        na=floor((x(ix^D,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix^D,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix^D,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix^D,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    {end do\}
    w(ixO^S,mom(:))=zero
    if(B0field) then
      w(ixO^S,mag(:))=zero
    else if(stagger_grid) then
      call b_from_vector_potential(ixGs^LL,ixI^L,ixO^L,block%ws,x)
      call mhd_face_to_center(ixO^L,block)
      w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
    else
      !MLB
      ! w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
      ! w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
      ! w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      call specialset_B0(ixI^L,ixO^L,x,Bf)!MLB
      ! w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir) !MLB
      w(ixO^S,mag(1))=Bf(ixO^S,1)!MLB
      w(ixO^S,mag(2))=Bf(ixO^S,2)!MLB
      w(ixO^S,mag(3))=Bf(ixO^S,3)!MLB
    endif

    call mhd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixI^L, ixC^L, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixI^L, ixC^L,idir
    double precision, intent(in)       :: xC(ixI^S,1:ndim)
    double precision, intent(out)      :: A(ixI^S)
    double precision :: rho1(ixI^S),rho2(ixI^S),rho3(ixI^S)
        rho1(ixC^S)=sqrt((xC(ixC^S,1)-xc1)**2+(xC(ixC^S,2)+hb)**2)
        rho2(ixC^S)=sqrt((xC(ixC^S,1)-xc2)**2+(xC(ixC^S,2)+hb2)**2)
        rho3(ixC^S)=sqrt((xC(ixC^S,1)-xc3)**2+(xC(ixC^S,2)+hb3)**2)
    if (idir==3) then
      ! A(ixC^S) = B0/ly*dcos(kx*xC(ixC^S,1))*dexp(-ly*xC(ixC^S,2))*dcos(theta)
      A(ixC^S) = -B0*(adipo**2-hb**2)/4.d0/adipo*log((adipo**2+rho1(ixC^S)**2-&
        2.d0*adipo*(xC(ixC^S,2)+hb))/(adipo**2+rho1(ixC^S)**2+&
        2.d0*adipo*(xC(ixC^S,2)+hb)))-B0*(adipo**2-hb2**2)/4.d0/adipo*&
        log((adipo**2+rho2(ixC^S)**2-2.d0*adipo*(xC(ixC^S,2)+hb2))/&
          (adipo**2+rho2(ixC^S)**2+2.d0*adipo*(xC(ixC^S,2)+hb2)))-&
        B0*factorcentralarcade*(adipo**2-hb3**2)/4.d0/adipo*&
        log((adipo**2+rho3(ixC^S)**2-2.d0*adipo*(xC(ixC^S,2)+hb3))&
          /(adipo**2+rho3(ixC^S)**2+2.d0*adipo*(xC(ixC^S,2)+hb3)))
      !  A(ixC^S) = -B0*(adipo**2-hb**2)/4.d0/adipo*log((adipo**2+((xC(ixC^S,1)-xc1)**2+(xC(ixC^S,2)+hb)**2)-&
      !  2.d0*adipo*(xC(ixC^S,2)+hb))/(adipo**2+((xC(ixC^S,1)-xc1)**2+(xC(ixC^S,2)+hb)**2)+&
      !  2.d0*adipo*(xC(ixC^S,2)+hb)))-B0*(adipo**2-hb2**2)/4.d0/adipo*&
      !  log((adipo**2+((xC(ixC^S,1)-xc2)**2+(xC(ixC^S,2)+hb2)**2)-2.d0*adipo*(xC(ixC^S,2)+hb2))/&
      !    (adipo**2+((xC(ixC^S,1)-xc2)**2+(xC(ixC^S,2)+hb2)**2)+2.d0*adipo*(xC(ixC^S,2)+hb2)))-&
      !  B0*factorcentralarcade*(adipo**2-hb3**2)/4.d0/adipo*&
      !  log((adipo**2+((xC(ixC^S,1)-xc3)**2+(xC(ixC^S,2)+hb3)**2)-2.d0*adipo*(xC(ixC^S,2)+hb3))&
      !    /(adipo**2+((xC(ixC^S,1)-xc3)**2+(xC(ixC^S,2)+hb3)**2)+2.d0*adipo*(xC(ixC^S,2)+hb3)))
    else
      A(ixC^S) = 0.d0
    end if

  end subroutine initvecpot_usr


  


  subroutine driven_velocity(ixI^L,ixO^L,qt,x,vdx,vdz)
   !Imposing a velocity on the bottom boundary to simulate convergence and shearing at the PIL.
   ! Here the shearing has been selected to be positive
   ! i.e. to increase the shear along the loops in time.
   use mod_global_parameters
   integer, intent(in)             :: ixI^L, ixO^L
   double precision, intent(in) :: qt, x(ixI^S,1:ndim)
   double precision, intent(out) :: vdx(ixI^S),vdz(ixI^S)
   double precision :: v_0t1,tb1,tr1,t1,t2, xcc1, ghost_region, resol, sigmam, conv_region
   double precision :: v_0t1cancellation,tb1cancellation,tr1cancellation,t1cancellation,t2cancellation
   double precision :: v_0t3, t3, t4, v_00, v_00cancellation



!!!!!!!!!!!!!!!Converging and  shearing velocity control for two FRs formation

  ! write(*,*)'**************TIME=',qt
   tb1=85000000000000.d0     !0.d0/unit_time+relax
   tr1=10.d0+tb1  !50.0d0/unit_time+tb1

   t1=10000000000000000000.d0+tr1 !(4000000.d0/unit_time)+tb1
   t2=50000.d0+t1    !(4000000.d0/unit_time)+tb1

   tb1cancellation=1000000000000.d0/unit_time!+relax !MLB
   tr1cancellation=tr1 !10.0d0/unit_time+relax !MLB

   t1cancellation=t1 !(90.d0/unit_time)+relax !MLB
   t2cancellation=t2 !(100.d0/unit_time)+relax !MLB



 !   print*,'XXXXXXXXX',unit_velocity


   xcc1=-6.4651938d0!3.2318599d0 ! debería ser -6.4651938
   v_00=0.7d0 !10.d5/unit_velocity !MLB
   v_00cancellation=0.0d0 !10.d5/unit_velocity

   resol=(xprobmax1-xprobmin1)/domain_nx1
   sigmam=0.1d0!0.12d0 !MLB
   conv_region=2.d0*sigmam !ahora no lo estoy usando

   if (qt>=tb1) then
     if (qt>=tb1.and.qt<tr1) then
       v_0t1=((qt-tb1)/(tr1-tb1))*v_00

     endif
     if (qt>=tr1.and.qt<t1) then
       v_0t1=v_00
     endif
     if (qt>=t1.and.qt<t2) then
       v_0t1=v_00*((t2-qt)/(t2-t1))
     endif
     if (qt>=t2) then
       v_0t1=0.0d0
     end if


       ! vdz(ixO^S)=v_0t1*dexp((-(x(ixO^S,1)-xcc1)**2)/sigmam**2)*dsin(2.d0*dpi*(x(ixO^S,1)-xcc1)/(conv_region))


!    perfil de cizallado usado anteriormente
       !vdz(ixO^S)=v_0t1*(dexp((-(x(ixO^S,1)-xcc1-sigmam)**2)/sigmam**2)-dexp((-(x(ixO^S,1)-xcc1+sigmam)**2)/sigmam**2))   !*dsin(2.d0*dpi*(x(ixO^S,1)-xcc1)/(conv_region))
!    nuevo perfil de cizallado
!El numero que aparece en la tanh me da lo abrupto que es el cizallado. Grande significa muy abrupto.
       vdz(ixO^S)=v_0t1*dexp((-(x(ixO^S,1)-xcc1)**2)/sigmam**2)*tanh(20.d0*(x(ixO^S,1)-xcc1))
       ! vdx(ixO^S)=v_0t1*(-dexp((-(x(ixO^S,1)-xcc1-sigmam)**2)/sigmam**2)+dexp((-(x(ixO^S,1)-xcc1+sigmam)**2)/sigmam**2))   !*dsin(2.d0*dpi*(x(ixO^S,1)-xcc1)/(conv_region))
    end if
!!!cancellation
    if (qt>=tb1cancellation) then
      if (qt>=tb1cancellation.and.qt<tr1cancellation) then
        v_0t1cancellation=((qt-tb1cancellation)/(tr1cancellation-tb1cancellation))*v_00cancellation
      endif
      if (qt>=tr1cancellation.and.qt<t1cancellation) then
        v_0t1cancellation=v_00cancellation
      endif
      if (qt>=t1cancellation.and.qt<t2cancellation) then
        v_0t1cancellation=v_00cancellation*((t2cancellation-qt)/(t2cancellation-t1cancellation))
      endif
      if (qt>=t2cancellation) then
        v_0t1cancellation=0.0d0
      end if

        vdx(ixO^S)=v_0t1cancellation*(-dexp((-(x(ixO^S,1)-xcc1-sigmam*2.d0)**2)/(sigmam*2.d0)**2)+dexp((-(x(ixO^S,1)-xcc1+sigmam*2.d0)**2)/(sigmam*2.d0)*2))   !*dsin(2.d0*dpi*(x(ixO^S,1)-xcc1)/(conv_region))
     end if



 end subroutine driven_velocity

  subroutine specialbound_usr(qt,ixI^L,ixO^L,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixO^L, iB, ixI^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: Bf(ixI^S,1:ndir) !MLB
    double precision :: pth(ixI^S),tmp(ixI^S),ggrid(ixI^S),invT(ixI^S)
    double precision :: Q(ixI^S),Qp(ixI^S)
    integer :: ix^D,ixOs^L,ixC^L,hxC^L,jxO^L,idir

    select case(iB)
    case(3)
      if(iprob<3) then
        !! fixed zero velocity
        do idir=1,ndir
          w(ixO^S,mom(idir))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        end do
      else
        ! call driven_velocity(ixI^L,ixO^L,qt,x,w(ixI^S,mom(1)))
        ! w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(2))&
        !             /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        ! w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(3))&
        !             /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        call driven_velocity(ixI^L,ixO^L,qt,x,w(ixI^S,mom(1)),w(ixI^S,mom(3)))
        !w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(1))/w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_) !MLB
        w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(2))/w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end if
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixO^S,mag(:))=0.d0
      else if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          ! 2nd order one-sided zero-gradient extrapolation
          !do ix2=ixOsmax2,ixOsmin2,-1
          !   block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
          !         (-block%ws(ix2+2^%2ixOs^S,idir)&
          !     +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          !end do
          ! 4th order one-sided equal-gradient extrapolation
          do ix2=ixOsmax2,ixOsmin2,-1
            block%ws(ix2^%2ixOs^S,idir)= &
              0.12d0*block%ws(ix2+5^%2ixOs^S,idir) &
             -0.76d0*block%ws(ix2+4^%2ixOs^S,idir) &
             +2.08d0*block%ws(ix2+3^%2ixOs^S,idir) &
             -3.36d0*block%ws(ix2+2^%2ixOs^S,idir) &
             +2.92d0*block%ws(ix2+1^%2ixOs^S,idir)
          end do
        end do
        ixOs^L=ixO^L-kr(2,^D);
        jxO^L=ixO^L+nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=Qp(ix2+1^%2ixO^S)*block%dvolume(ix2+1^%2ixO^S)&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        if(iprob<3) then
          !MLB
          ! w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
          call specialset_B0(ixI^L,ixO^L,x,Bf)!MLB
          w(ixO^S,mag(3))=Bf(ixO^S,3)!MLB
        else
          do ix2=ixOmax2,ixOmin2,-1
            w(ix2^%2ixO^S,mag(3))= &
              0.12d0*w(ix2+5^%2ixO^S,mag(3)) &
             -0.76d0*w(ix2+4^%2ixO^S,mag(3)) &
             +2.08d0*w(ix2+3^%2ixO^S,mag(3)) &
             -3.36d0*w(ix2+2^%2ixO^S,mag(3)) &
             +2.92d0*w(ix2+1^%2ixO^S,mag(3))
          end do
        end if
      else
        !MLB
        call specialset_B0(ixI^L,ixO^L,x,Bf)
        w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir) !MLB
        ! w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
        ! w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
        ! w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
        ! w(ixO^S,mag(1))= Bf(ixO^S,1)
        ! w(ixO^S,mag(2))= Bf(ixO^S,2)
        ! w(ixO^S,mag(3))= Bf(ixO^S,3)
      endif
      !! fixed gravity stratification of density and pressure pre-determined in initial condition
      do ix2=ixOmin2,ixOmax2
        w(ixOmin1:ixOmax1,ix2,rho_)=rbc(ix2)
        w(ixOmin1:ixOmax1,ix2,p_)=pbc(ix2)
      enddo
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case(4)
      ixOs^L=ixO^L;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixI^L,ixOs^L,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixI^L,ixOs^L,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin2-1^%2ixO^S)=w(ixOmin2-1^%2ixO^S,rho_)/pth(ixOmin2-1^%2ixO^S)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin2-1^%2ixO^S)=tmp(ixOmin2-1^%2ixO^S)+0.5d0*&
            (ggrid(ix2^%2ixO^S)+ggrid(ix2-1^%2ixO^S))*invT(ixOmin2-1^%2ixO^S)
        w(ix2^%2ixO^S,p_)=pth(ixOmin2-1^%2ixO^S)*dexp(tmp(ixOmin2-1^%2ixO^S)*dxlevel(2))
        w(ix2^%2ixO^S,rho_)=w(ix2^%2ixO^S,p_)*invT(ixOmin2-1^%2ixO^S)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixO^S,mom(idir)) =-w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))&
                     /w(ixOmin1:ixOmax1,ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax^D=ixOmax^D;
          ixOsmin^D=ixOmin^D-kr(^D,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
                   (-block%ws(ix2-2^%2ixOs^S,idir)&
               +4.d0*block%ws(ix2-1^%2ixOs^S,idir))
          end do
        end do
        ixOs^L=ixO^L;
        jxO^L=ixO^L-nghostcells*kr(2,^D);
        block%ws(ixOs^S,2)=zero
        call get_divb(w,ixI^L,jxO^L,Q)
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixI^L,ixO^L,Qp)
          block%ws(ix2^%2ixOs^S,2)=&
            (Q(jxOmax2^%2jxO^S)*block%dvolume(jxOmax2^%2jxO^S)&
           -Qp(ix2^%2ixO^S)*block%dvolume(ix2^%2ixO^S))&
            /block%surfaceC(ix2^%2ixOs^S,2)
        end do
        call mhd_face_to_center(ixO^L,block)
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(3))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2-2,mag(3))&
                 +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(3)))
        enddo
      else
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* &
                      (-w(ixOmin1:ixOmax1,ix2-2,mag(:))&
                 +4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
        enddo
      end if
      call mhd_to_conserved(ixI^L,ixO^L,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine gravity(ixI^L,ixO^L,wCT,x,gravity_field)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision                :: ggrid(ixI^S)

    gravity_field=0.d0
    call getggrav(ggrid,ixI^L,ixO^L,x)
    gravity_field(ixO^S,2)=ggrid(ixO^S)

  end subroutine gravity

  subroutine getggrav(ggrid,ixI^L,ixO^L,x)
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: ggrid(ixI^S)

    ggrid(ixO^S)=usr_grav*(SRadius/(SRadius+x(ixO^S,2)))**2
  end subroutine

  subroutine special_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S),bQgrid(ixI^S)

    ! add global background heating bQ
    call getbQ(bQgrid,ixI^L,ixO^L,qtC,wCT,x)
    w(ixO^S,e_)=w(ixO^S,e_)+qdt*bQgrid(ixO^S)
    if(iprob>1) then !if(iprob==2) then !MLB
      call getlQ(lQgrid,ixI^L,ixO^L,qt,wCT,x)
      w(ixO^S,e_)=w(ixO^S,e_)+qdt*lQgrid(ixO^S)
    end if

  end subroutine special_source

  subroutine getbQ(bQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate background heating bQ
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: bQgrid(ixI^S)
    double precision :: bQgridspecial1(ixI^S),bQgridspecial2(ixI^S)
    double precision :: envolvente(ixI^S),aenvolvente
    double precision :: bQgridleft(ixI^S),bQgridright(ixI^S),bQgridcenter(ixI^S)
    double precision :: rraau(ixI^S)
    double precision :: xcenterleft,xcenterright,ycenterleft,ycenterright
    double precision :: xcentercenter,ycentercenter
    double precision :: rcleft(ixI^S),rcright(ixI^S),rccenter(ixI^S),A(ixI^S)
    double precision :: radiusthreshold,radiusthreshold2,radiusthreshold3
    double precision :: lambdabackground
    ! xcenterleft= -2.5d0!-4.4d0
    ! xcenterright= 2.5d0!4.4d0
    ! xcentercenter=0.d0
    ! ycenterleft= -0.6d0!-3.9d0
    ! ycenterright=-0.6d0!-3.9d0
    ! ycentercenter=-0.5d0
    ! radiusthreshold=1.5d0!5.2d0 ;right
    ! radiusthreshold2=1.5d0!5.2d0 ;left
    ! radiusthreshold3=1.1d0
    ! bQgridleft=0.d0
    ! bQgridright=0.d0
    ! bQgridcenter=0.d0

    lambdabackground=5.d0





      call initvecpot_usr(ixI^L, ixO^L, x, A, 3)
      
      
      aenvolvente=-140.+8.
      envolvente=0.5d0*(tanh(-0.7*(A(ixO^S)-aenvolvente))+1.d0)

      bQgridspecial1(ixO^S)=8.d0*bQ0*(abs(A(ixO^S))/171.73493d0)**(6.d0)*envolvente
      
      bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/lambdabackground)+ &
      bQgridspecial1(ixO^S)

  end subroutine getbQ

  

  subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
  ! calculate localized heating lQ
  !nueva subrutina MLB
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)

    double precision :: lQgrid(ixI^S), A(ixI^S), lQd, lQ0, hp, lQt,Atop,Alow
    double precision :: auxradius1(ixI^S),xcenter1,ycenter1,auxm1,y11(ixI^S),y12(ixI^S)
    double precision :: radiuscentral(ixI^S)
    double precision :: auxradius2(ixI^S),xcenter2,ycenter2,auxm2,y21(ixI^S),y22(ixI^S)
    double precision :: lambdamin,lambdamax,lambdamax2,lambda1(ixI^S),lambda11(ixI^S)
    double precision :: lambda2(ixI^S),lambda22(ixI^S),xprimaL(ixI^S),zprimaL(ixI^S)
    double precision :: xprimaR(ixI^S),zprimaR(ixI^S)
    double precision :: thetaangle01,thetaangle1(ixI^S),thetaangle02,thetaangle2(ixI^S)
    double precision :: lQgridleft(ixI^S),lQgridright(ixI^S)
    double precision :: recta1(ixI^S),recta2(ixI^S)
    double precision :: sigma_heating,sigma_heating_impulsive,y1,lam1
    double precision :: impulsiveheatingfactor,impulsiveheatingduration,timeimpulse
    double precision :: impulsiveheatingfactoramplitude, A0

    lQ0=5.d0*1.d-2/heatunit!MLB 1.d-2/heatunit
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
      xcenter1=-8.5744060!-5.2288454d0 !LEFT
      xcenter2= 8.5744060!5.2288454d0  !RIGHT
      A0=-132.0d0 !vector potential at the center of the field lines bundle
      !
      ! !impulsive heating
      impulsiveheatingfactor=1.d0
      impulsiveheatingduration=20.d0 ! every unit is 85 secs
      timeimpulse=10.001d0
      impulsiveheatingfactoramplitude=5.d0
      if(qt>timeimpulse .and. qt<timeimpulse + impulsiveheatingduration) then
        if(qt<timeimpulse+impulsiveheatingduration/3.d0) then
          impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude* &
          (3.d0*(qt-timeimpulse)/impulsiveheatingduration)
        endif
        if(qt>=timeimpulse+impulsiveheatingduration/3.d0 .and. &
         qt < timeimpulse+impulsiveheatingduration*2.d0/3.d0) then
         impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude
        endif
        if(qt>=timeimpulse+impulsiveheatingduration*2.d0/3.d0 .and. &
        qt<timeimpulse+impulsiveheatingduration) then
         impulsiveheatingfactor=1.d0+impulsiveheatingfactoramplitude* &
         (-3.d0/impulsiveheatingduration*&
         (qt-timeimpulse-2.d0/3.d0*impulsiveheatingduration)+1.d0)
        endif
      
      endif
      
      !El calentamiento en la derecha lo he quitado: *0.d0
      !El impulsivehatingfactor vale 1 fuera de timeimpulse y timeimpulse+impulsiveheatingduration
      !Y dentro de esos tiempos vale 1+impulsiveheatingfactoramplitude
      !Para que no haya calentamiento lo que hacemos es que la función de calentamiento de abajo
      !vaya multiplicada por impulsiveheatingfactor-1.d0
      
      lQgrid(ixO^S)=lQ0*((impulsiveheatingfactor-1.d0)*&
        dexp(-(x(ixO^S,1)-xcenter1)**2/sigma_heating**2)+&
        0.d0*dexp(-(x(ixO^S,1)-xcenter2)**2/sigma_heating**2) )
      lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(x(ixO^S,2)-y1)**2/lam1**2)

      call initvecpot_usr(ixI^L, ixO^L, x, A, 3)
      lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(A(ixO^S)-A0)**2/2.d0/(0.50d0)**2)

                         !
                         ! +&
                         ! lQ0*(impulsiveheatingfactor-1.d0)*&
                         ! dexp(-(x(ixO^S,1)-xcenter1)**2/sigma_heating_impulsive**2)*&
                         ! dexp(-(x(ixO^S,2)-y1)**2/lam1**2)


      ! where(A(ixO^S)<Atop .and. A(ixO^S)>Alow)
      !   lQgrid(ixO^S)=lQ0
      ! end where

    ! !MLB
    if(qt<10.d0) then
      lQgrid(ixO^S)=0.d0
    endif
    if(qt>10.d0 .and. qt<60.d0) then
      lQgrid(ixO^S)=(qt-10.d0)/50.d0*lQgrid(ixO^S)
    endif
    


  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixI^L,ixO^L,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixI^L, ixO^L
    double precision, intent(in) :: qt, w(ixI^S,1:nw), x(ixI^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixO^S,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine special_refine_grid

subroutine specialvar_output(ixI^L,ixO^L,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixI^L,ixO^L
    double precision, intent(in)       :: x(ixI^S,1:ndim)
    double precision                   :: w(ixI^S,nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixI^S),B2(ixI^S),tmp2(ixI^S),dRdT(ixI^S)
    double precision :: ens(ixI^S),divb(ixI^S),wlocal(ixI^S,1:nw)
    double precision :: Btotal(ixI^S,1:ndir),curlvec(ixI^S,1:ndir)
    double precision :: Te(ixI^S),tco_local
    double precision, dimension(ixI^S,1:ndim) :: gradT, bunitvec
    integer :: idirmin,idir,ix^D
    logical :: lrlt(ixI^S)

    wlocal(ixI^S,1:nw)=w(ixI^S,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixI^L,ixI^L,pth)
    Te(ixI^S)=pth(ixI^S)/w(ixI^S,rho_)
    w(ixO^S,nw+1)=Te(ixO^S)

    do idir=1,ndir
      if(B0field) then
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))+block%B0(ixI^S,idir,0)
      else
        Btotal(ixI^S,idir)=w(ixI^S,mag(idir))
      endif
    end do
    ! B^2
    B2(ixO^S)=sum((Btotal(ixO^S,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixO^S,nw+2)=dsqrt(B2(ixO^S)/w(ixO^S,rho_))

    ! output divB1
    call get_divb(wlocal,ixI^L,ixO^L,divb)
    w(ixO^S,nw+3)=divb(ixO^S)
    ! output the plasma beta p*2/B**2
    w(ixO^S,nw+4)=pth(ixO^S)*two/B2(ixO^S)
    ! output heating rate
    call getbQ(ens,ixI^L,ixO^L,global_time,wlocal,x)
    w(ixO^S,nw+5)=ens(ixO^S)
    ! store the cooling rate
    if(mhd_radiative_cooling)call getvar_cooling(ixI^L,ixO^L,wlocal,x,ens,rc_fl)
    w(ixO^S,nw+6)=ens(ixO^S)

    ! store current
    call get_current(wlocal,ixI^L,ixO^L,idirmin,curlvec)
    do idir=1,ndir
      w(ixO^S,nw+6+idir)=curlvec(ixO^S,idir)
    end do

    if(nw_extra>0) w(ixO^S,nw+10)=block%wextra(ixO^S,iw_tcoff)

  end subroutine specialvar_output

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
       character(len=*) :: varnames
       varnames = 'bt1 bt2 bt3'
  end subroutine specialvarnames_output


  subroutine set_output_vars(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x)
      integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
      double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
      double precision, intent(inout)    :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw+nwauxio)
      double precision, intent(in) :: qt
      integer :: idir
      double precision :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
      
      do idir=1,ndir
         if(B0field) then
            Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
            ixImin2:ixImax2,mag(idir))+block%B0(ixImin1:ixImax1,ixImin2:ixImax2,&
            idir,0)
         else
            Btotal(ixImin1:ixImax1,ixImin2:ixImax2,idir)=w(ixImin1:ixImax1,&
            ixImin2:ixImax2,mag(idir))
         endif
      end do
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bt_(1)) = Btotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bt_(2)) = Btotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,bt_(3)) = Btotal(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)

  end subroutine set_output_vars

  ! subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! ! Here add a steady (time-independent) potential or
  ! ! linear force-free background field
  !   integer, intent(in)           :: ixI^L,ixO^L
  !   double precision, intent(in)  :: x(ixI^S,1:ndim)
  !   double precision, intent(inout) :: wB0(ixI^S,1:ndir)
  !
  !   wB0(ixO^S,1)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
  !   wB0(ixO^S,2)=+B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
  !   wB0(ixO^S,3)=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
  !
  ! end subroutine specialset_B0

  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
 ! Here add a time-independent background magnetic field
   use mod_global_parameters
   integer, intent(in)           :: ixI^L,ixO^L
   double precision, intent(in)  :: x(ixI^S,1:ndim)
   double precision, intent(inout) :: wB0(ixI^S,1:ndir)
   double precision :: Bx1(ixI^S),By1(ixI^S),Bx2(ixI^S),By2(ixI^S),Bx3(ixI^S),By3(ixI^S)
   double precision :: Bradial1(ixI^S),Btheta1(ixI^S),Bradial2(ixI^S),Btheta2(ixI^S)
   double precision :: rho1(ixI^S),cosphi1(ixI^S),sinphi1(ixI^S),cos2phi1(ixI^S)
   double precision :: rho2(ixI^S),cosphi2(ixI^S),sinphi2(ixI^S),cos2phi2(ixI^S)
   ! double precision :: Bx1(ixI^S),By1(ixI^S)
   double precision :: xn1(ixI^S),yn1(ixI^S),xn2(ixI^S),yn2(ixI^S),xn3(ixI^S),yn3(ixI^S)

   xn1(ixO^S)=(x(ixO^S,1)-xc1)
   yn1(ixO^S)=(x(ixO^S,2)+hb)

   xn2(ixO^S)=(x(ixO^S,1)-xc2)
   yn2(ixO^S)=(x(ixO^S,2)+hb2)

   xn3(ixO^S)=(x(ixO^S,1)-xc3)
   yn3(ixO^S)=(x(ixO^S,2)+hb3)

   Bx1(ixO^S)=-B0*((adipo**2-hb**2)*(adipo**2+xn1(ixO^S)**2-yn1(ixO^S)**2))/((adipo**2+xn1(ixO^S)**2+yn1(ixO^S)**2-2.d0*adipo*yn1(ixO^S))*(adipo**2+xn1(ixO^S)**2+yn1(ixO^S)**2+2.d0*adipo*yn1(ixO^S)))
   By1(ixO^S)=-B0*(2.d0*(adipo**2-hb**2)*xn1(ixO^S)*yn1(ixO^S))/((adipo**2+xn1(ixO^S)**2+yn1(ixO^S)**2-2.d0*adipo*yn1(ixO^S))*(adipo**2+xn1(ixO^S)**2+yn1(ixO^S)**2+2.d0*adipo*yn1(ixO^S)))

   Bx2(ixO^S)=-B0*((adipo**2-hb2**2)*(adipo**2+xn2(ixO^S)**2-yn2(ixO^S)**2))/((adipo**2+xn2(ixO^S)**2+yn2(ixO^S)**2-2.d0*adipo*yn2(ixO^S))*(adipo**2+xn2(ixO^S)**2+yn2(ixO^S)**2+2.d0*adipo*yn2(ixO^S)))
   By2(ixO^S)=-B0*(2.d0*(adipo**2-hb2**2)*xn2(ixO^S)*yn2(ixO^S))/((adipo**2+xn2(ixO^S)**2+yn2(ixO^S)**2-2.d0*adipo*yn2(ixO^S))*(adipo**2+xn2(ixO^S)**2+yn2(ixO^S)**2+2.d0*adipo*yn2(ixO^S)))

   Bx3(ixO^S)=-B0*factorcentralarcade*((adipo**2-hb3**2)*(adipo**2+xn3(ixO^S)**2-yn3(ixO^S)**2))/((adipo**2+xn3(ixO^S)**2+yn3(ixO^S)**2-2.d0*adipo*yn3(ixO^S))*(adipo**2+xn3(ixO^S)**2+yn3(ixO^S)**2+2.d0*adipo*yn3(ixO^S)))
   By3(ixO^S)=-B0*factorcentralarcade*(2.d0*(adipo**2-hb3**2)*xn3(ixO^S)*yn3(ixO^S))/((adipo**2+xn3(ixO^S)**2+yn3(ixO^S)**2-2.d0*adipo*yn3(ixO^S))*(adipo**2+xn3(ixO^S)**2+yn3(ixO^S)**2+2.d0*adipo*yn3(ixO^S)))

   ! rho1(ixO^S)=sqrt((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)+hb)**2)
   ! cosphi1(ixO^S)=(x(ixO^S,1)-xc1)/rho1(ixO^S)
   ! sinphi1(ixO^S)=(x(ixO^S,2)+hb)/rho1(ixO^S)
   ! cos2phi1(ixO^S)=1.d0-2.d0*(sinphi1(ixO^S))**2
   ! Bradial1(ixO^S)=B0*((adipo**2-hb**2)*rho1(ixO^S)*(adipo**2+(rho1(ixO^S))**2)* &
   !    cosphi1(ixO^S))/(adipo**4+(rho1(ixO^S))**4+2.d0*adipo**2*(rho1(ixO^S))**2*cos2phi1(ixO^S))
   ! Btheta1(ixO^S)=-B0*((adipo**2-hb**2)*(adipo**2-(rho1(ixO^S))**2)*sinphi1(ixO^S))/&
   !  (adipo**4+(rho1(ixO^S))**4+2.d0*adipo**2*(rho1(ixO^S))**2*cos2phi1(ixO^S))
   !
   !  rho2(ixO^S)=sqrt((x(ixO^S,1)-xc2)**2+(x(ixO^S,2)+hb)**2)
   !  cosphi2(ixO^S)=(x(ixO^S,1)-xc2)/rho2(ixO^S)
   !  sinphi2(ixO^S)=(x(ixO^S,2)+hb)/rho2(ixO^S)
   !  cos2phi2(ixO^S)=1.d0-2.d0*(sinphi2(ixO^S))**2
   !  Bradial2(ixO^S)=B0*((adipo**2-hb**2)*rho2(ixO^S)*(adipo**2+(rho2(ixO^S))**2)* &
   !     cosphi2(ixO^S))/(adipo**4+(rho2(ixO^S))**4+2.d0*adipo**2*(rho2(ixO^S))**2*cos2phi2(ixO^S))
   !  Btheta2(ixO^S)=-B0*((adipo**2-hb**2)*(adipo**2-(rho2(ixO^S))**2)*sinphi2(ixO^S))/&
   !   (adipo**4+(rho2(ixO^S))**4+2.d0*adipo**2*(rho2(ixO^S))**2*cos2phi2(ixO^S))


   ! Bx1(ixO^S)=md1*((x(ixO^S,1)-xc1)**2-(x(ixO^S,2)+hb)**2)/(((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)+hb)**2)**2)
   !
   ! By1(ixO^S)=2.d0*md1*(x(ixO^S,1)-xc1)*(x(ixO^S,2)+hb)/(((x(ixO^S,1)-xc1)**2+(x(ixO^S,2)+hb)**2)**2)
   !
   ! Bx2(ixO^S)=md2*((x(ixO^S,1)-xc2)**2-(x(ixO^S,2)+hb)**2)/(((x(ixO^S,1)-xc2)**2+(x(ixO^S,2)+hb)**2)**2)
   !
   ! By2(ixO^S)=2.d0*md2*(x(ixO^S,1)-xc2)*(x(ixO^S,2)+hb)/(((x(ixO^S,1)-xc2)**2+(x(ixO^S,2)+hb)**2)**2)
   !
   ! Bx3(ixO^S)=md3*((x(ixO^S,1)-xc3)**2-(x(ixO^S,2)+hb3)**2)/(((x(ixO^S,1)-xc3)**2+(x(ixO^S,2)+hb3)**2)**2)
   !
   ! By3(ixO^S)=2.d0*md3*(x(ixO^S,1)-xc3)*(x(ixO^S,2)+hb3)/(((x(ixO^S,1)-xc3)**2+(x(ixO^S,2)+hb3)**2)**2)

   wB0(ixO^S,1)=Bx1(ixO^S)+Bx2(ixO^S)+Bx3(ixO^S)
   wB0(ixO^S,2)=By1(ixO^S)+By2(ixO^S)+By3(ixO^S)
   wB0(ixO^S,3)=0.0!*x(ixO^S,2)
 end subroutine specialset_B0







  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    !PURPOSE: use different tolerance in special regions for AMR to
    !reduce/increase resolution there where nothing/something interesting happens.
    use mod_global_parameters

    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    integer, intent(in) :: level

    double precision :: bczone^D,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.6d0
    ! thickness of near-boundary region
    bczone1=0.8d0
    bczone2=2.d0
    ! linear changing of additional tolerance
    if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-min(xlocal(1)-xprobmin1,xprobmax1-xlocal(1))/bczone1)*addtol
    endif
    if(xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-(xprobmax2-xlocal(2))/bczone2)*addtol
    endif
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

end module mod_usr
