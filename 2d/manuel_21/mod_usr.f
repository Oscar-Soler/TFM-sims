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
    unit_numberdensity = 1.d9 !cm-3,cm-3

    usr_set_parameters  => initglobaldata_usr
    usr_init_one_grid   => initonegrid_usr
    usr_special_bc      => specialbound_usr
    usr_source          => special_source
    usr_gravity         => gravity
    usr_refine_grid     => special_refine_grid
    usr_set_B0          => specialset_B0
    !usr_aux_output      => specialvar_output
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
    heatunit=unit_pressure/unit_time !3.697693390805347E-003 erg*cm-3/s,erg*cm-3/s

    usr_grav=-2.74d4*unit_length/unit_velocity**2 ! solar gravity
    bQ0=6.0d-5/heatunit !background heating power density !MLB lo he dividido entre 2, 1d-4
    gzone=0.2d0 ! thickness of a ghostzone below the bottom boundary !MLB
    dya=(2.d0*gzone+xprobmax2-xprobmin2)/dble(jmax) !cells size of high-resolution 1D solar atmosphere
    B0=Busr/unit_magneticfield ! magnetic field strength at the bottom
    theta=60.d0*dpi/180.d0 !the angle to the plane xy, 90-theta is the angle to the polarity inversion line of the arcade
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
      rpho=1.151d15/unit_numberdensity !number density at the bottom of height table
      Tpho=8.d3/unit_temperature ! temperature of chromosphere
      Ttop=1.5d6/unit_temperature ! estimated temperature in the top
      htra=0.2d0 ! height of initial transition region
      wtra=0.02d0 ! width of initial transition region 
      Ttr=1.6d5/unit_temperature ! lowest temperature of upper profile
      Fc=2.d5/heatunit/unit_length ! constant thermal conduction flux
      kappa=8.d-7*unit_temperature&
         **3.5d0/unit_length/unit_density/unit_velocity**3
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
      !! solution of hydrostatic equation
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
      na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0))/dya+0.5d0)
      res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+&
         1)-0.5d0)-(dble(na)-0.5d0)*dya
      rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
      pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
    end do

    if (mype==0) then
     print*,'minra',minval(ra)
     print*,'rhob',rhob
     print*,'pb',pb
    endif

  end subroutine inithdstatic

  ! subroutine inithdstatic !Valeriia's code
  !   use mod_solar_atmosphere
  !   ! initialize the table in a vertical line through the global domain
  !   integer :: j,na,ibc
  !   double precision, allocatable :: Ta(:),gg(:),ya(:)
  !   double precision :: rch,Ttop,Tch,wtra,res,rhob,pb,htra,invT
  !
  !
  !   allocate(ya(jmax),Ta(jmax),gg(jmax),pa(jmax),ra(jmax))
  !
  !     rch=1.15d15/unit_numberdensity ! number density at the bottom of height table
  !     Tch=8.d3/unit_temperature ! temperature of chromosphere
  !     Ttop=1.d6/unit_temperature ! estimated temperature in the top
  !     htra=0.1d0 ! height of initial transition region
  !     wtra=0.02d0 ! width of initial transition region
  !     do j=1,jmax
  !        ya(j)=(dble(j)-0.5d0)*dya-gzone
  !        Ta(j)=Tch+0.5d0*(Ttop-Tch)*(tanh((ya(j)-htra-0.027d0)/wtra)+1.d0)
  !        gg(j)=usr_grav*(SRadius/(SRadius+ya(j)))**2
  !     enddo
  !     !! solution of hydrostatic equation
  !     ra(1)=rch
  !     pa(1)=rch*Tch
  !     invT=gg(1)/Ta(1)
  !     invT=0.d0
  !     do j=2,jmax
  !        invT=invT+(gg(j)/Ta(j)+gg(j-1)/Ta(j-1))*0.5d0
  !        pa(j)=pa(1)*dexp(invT*dya)
  !        ra(j)=pa(j)/Ta(j)
  !     end do
  !   deallocate(ya,gg,Ta)
  !   !! initialized rho and p in the fixed bottom boundary
  !   na=floor(gzone/dya+0.5d0)
  !   res=gzone-(dble(na)-0.5d0)*dya
  !   rhob=ra(na)+res/dya*(ra(na+1)-ra(na))
  !   pb=pa(na)+res/dya*(pa(na+1)-pa(na))
  !   allocate(rbc(nghostcells))
  !   allocate(pbc(nghostcells))
  !   do ibc=nghostcells,1,-1
  !     na=floor((gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0))/dya+0.5d0)
  !     res=gzone-dx(2,refine_max_level)*(dble(nghostcells-ibc+1)-0.5d0)-(dble(na)-0.5d0)*dya
  !     rbc(ibc)=ra(na)+res/dya*(ra(na+1)-ra(na))
  !     pbc(ibc)=pa(na)+res/dya*(pa(na+1)-pa(na))
  !   end do
  !
  !   if (mype==0) then
  !    print*,'minra',minval(ra)
  !    print*,'rhob',rhob
  !    print*,'pb',pb
  !   endif
  !
  ! end subroutine inithdstatic

  subroutine initonegrid_usr(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x)
    ! initialize one grid
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir) !MLB
    double precision :: res
    integer :: ix1,ix2,na
    logical, save :: first=.true.

    if(first)then
      if(mype==0) then
        write(*,*)'Simulating prominence perturbation'
      endif
      first=.false.
    endif
    do ix2=ixOmin2,ixOmax2
    do ix1=ixOmin1,ixOmax1
        na=floor((x(ix1,ix2,2)-xprobmin2+gzone)/dya+0.5d0)
        res=x(ix1,ix2,2)-xprobmin2+gzone-(dble(na)-0.5d0)*dya
        w(ix1,ix2,rho_)=ra(na)+(one-cos(dpi*res/dya))/two*(ra(na+1)-ra(na))
        w(ix1,ix2,p_)  =pa(na)+(one-cos(dpi*res/dya))/two*(pa(na+1)-pa(na))
    end do
    end do
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(:))=zero
    if(B0field) then
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=zero
    else if(stagger_grid) then
      call b_from_vector_potential(ixGslo1,ixGslo2,ixGshi1,ixGshi2,ixImin1,&
         ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,block%ws,x)
      call mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,block)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=-B0*dcos(kx*x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1))*dexp(-ly*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2))*dsin(theta)
    else
      !MLB
      ! w(ixO^S,mag(1))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dcos(theta)
      ! w(ixO^S,mag(2))= B0*dsin(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))
      ! w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
      call specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,x,Bf) !MLB
      ! w(ixO^S,mag(1:ndir))=Bf(ixO^S,1:ndir) !MLB
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1))=Bf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1) !MLB
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(2))=Bf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,2) !MLB
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=Bf(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,3) !MLB
    endif

    call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,w,x)

  end subroutine initonegrid_usr

  subroutine initvecpot_usr(ixImin1,ixImin2,ixImax1,ixImax2, ixCmin1,ixCmin2,&
     ixCmax1,ixCmax2, xC, A, idir)
    ! initialize the vectorpotential on the edges
    ! used by b_from_vectorpotential()
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixCmin1,ixCmin2,ixCmax1,ixCmax2,idir
    double precision, intent(in)       :: xC(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)      :: A(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: rho1(ixImin1:ixImax1,ixImin2:ixImax2),&
       rho2(ixImin1:ixImax1,ixImin2:ixImax2),rho3(ixImin1:ixImax1,&
       ixImin2:ixImax2)
        rho1(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt((xC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,1)-xc1)**2+(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           2)+hb)**2)
        rho2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt((xC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,1)-xc2)**2+(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           2)+hb2)**2)
        rho3(ixCmin1:ixCmax1,ixCmin2:ixCmax2)=sqrt((xC(ixCmin1:ixCmax1,&
           ixCmin2:ixCmax2,1)-xc3)**2+(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
           2)+hb3)**2)
    if (idir==3) then
      ! A(ixC^S) = B0/ly*dcos(kx*xC(ixC^S,1))*dexp(-ly*xC(ixC^S,2))*dcos(theta)
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = &
         -B0*(adipo**2-hb**2)/4.d0/adipo*log((adipo**2+rho1(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2-2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+hb))/(adipo**2+rho1(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2+2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+hb)))-B0*(adipo**2-hb2**2)/4.d0/adipo*log((adipo**2+&
         rho2(ixCmin1:ixCmax1,ixCmin2:ixCmax2)**2-&
         2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+hb2))/(adipo**2+rho2(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2+2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+hb2)))-B0*factorcentralarcade*(adipo**2-&
         hb3**2)/4.d0/adipo*log((adipo**2+rho3(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2-2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+hb3))/(adipo**2+rho3(ixCmin1:ixCmax1,&
         ixCmin2:ixCmax2)**2+2.d0*adipo*(xC(ixCmin1:ixCmax1,ixCmin2:ixCmax2,&
         2)+hb3)))
      !  A(ixC^S) = -B0*(adipo**2-hb**2)/4.d0/adipo*log((adipo**2+((xC(ixC^S,1)-xc1)**2+(xC(ixC^S,2)+hb)**2)-&
      !  2.d0*adipo*(xC(ixC^S,2)+hb))/(adipo**2+((xC(ixC^S,1)-xc1)**2+(xC(ixC^S,2)+hb)**2)+&
      !  2.d0*adipo*(xC(ixC^S,2)+hb)))-B0*(adipo**2-hb2**2)/4.d0/adipo*&
      !  log((adipo**2+((xC(ixC^S,1)-xc2)**2+(xC(ixC^S,2)+hb2)**2)-2.d0*adipo*(xC(ixC^S,2)+hb2))/&
      !    (adipo**2+((xC(ixC^S,1)-xc2)**2+(xC(ixC^S,2)+hb2)**2)+2.d0*adipo*(xC(ixC^S,2)+hb2)))-&
      !  B0*factorcentralarcade*(adipo**2-hb3**2)/4.d0/adipo*&
      !  log((adipo**2+((xC(ixC^S,1)-xc3)**2+(xC(ixC^S,2)+hb3)**2)-2.d0*adipo*(xC(ixC^S,2)+hb3))&
      !    /(adipo**2+((xC(ixC^S,1)-xc3)**2+(xC(ixC^S,2)+hb3)**2)+2.d0*adipo*(xC(ixC^S,2)+hb3)))
    else
      A(ixCmin1:ixCmax1,ixCmin2:ixCmax2) = 0.d0
    end if

  end subroutine initvecpot_usr














  ! subroutine boundary_electric_field(ixI^L,ixO^L,qt,qdt,fE,s)
  !   ! specify tangential electric field at physical boundaries
  !   ! to fix or drive normal magnetic field
  !   integer, intent(in)                :: ixI^L, ixO^L
  !   double precision, intent(in)       :: qt,qdt
  !   type(state)                        :: s
  !   double precision, intent(inout)    :: fE(ixI^S,7-2*ndim:3)
  !
  !   double precision :: xC(ixI^S,1:ndim), tmp(ixI^S)
  !   integer :: idir,ixC^L,ixA^L
  !
  !   if(s%is_physical_boundary(3)) then
  !     if(iprob<3) then
  !       ! to fix normal magnetic field at bottom boundary surface
  !       idir=3
  !       ixCmin^D=ixOmin^D+kr(idir,^D)-1;
  !       ixCmax^D=ixOmax^D;
  !       fE(ixCmin2^%2ixC^S,3)=0.d0
  !     else
  !       ! horizontal flow driven boundary
  !       ixCmin^D=ixOmin^D-1;
  !       ixCmax^D=ixOmax^D;
  !       ixAmin^D=ixCmin^D;
  !       ixAmax^D=ixCmax^D+1;
  !       !! cell center electric field
  !       tmp(ixAmin2^%2ixA^S)=-s%ws(ixAmin2^%2ixA^S,2)*s%w(ixAmin2^%2ixA^S,mom(1))/s%w(ixAmin2^%2ixA^S,rho_)
  !       ! neighor cell center average to get cell edge
  !       ixAmin^D=ixCmin^D+kr(1,^D);
  !       ixAmax^D=ixCmax^D+kr(1,^D);
  !       !print*,'difference',maxval(abs(fE(ixCmin2^%2ixC^S,3)-vx(ixCmin2^%2ixC^S)))
  !       fE(ixCmin2^%2ixC^S,3)=0.5d0*(tmp(ixCmin2^%2ixC^S)+tmp(ixAmin2^%2ixA^S))
  !     end if
  !   end if

  ! end subroutine boundary_electric_field






















  !> allow user to specify variables' left and right state at physical boundaries to control flux through the boundary surface
  ! subroutine boundary_wLR(ixI^L,ixO^L,qt,wLC,wRC,wLp,wRp,s,idir)
  !   integer, intent(in)             :: ixI^L, ixO^L, idir
  !   double precision, intent(in)    :: qt
  !   double precision, intent(inout) :: wLC(ixI^S,1:nw), wRC(ixI^S,1:nw)
  !   double precision, intent(inout) :: wLp(ixI^S,1:nw), wRp(ixI^S,1:nw)
  !   type(state)                     :: s
  !
  !   double precision :: vx(ixI^S)
  !
  !   if(iprob<3) then
  !     if(s%is_physical_boundary(3).and.idir==2) then
  !       wLp(ixOmin2^%2ixO^S,mom(:))=0.d0
  !       wRp(ixOmin2^%2ixO^S,mom(:))=0.d0
  !       wLC(ixOmin2^%2ixO^S,mom(:))=0.d0
  !       wRC(ixOmin2^%2ixO^S,mom(:))=0.d0
  !     end if
  !   else
  !     if(s%is_physical_boundary(3).and.idir==2) then
  !       call driven_velocity(ixI^L,ixO^L,qt,s%x,vx)
  !       wLp(ixOmin2^%2ixO^S,mom(1))=vx(ixOmin2^%2ixO^S)
  !       wRp(ixOmin2^%2ixO^S,mom(1))=wRp(ixOmin2^%2ixO^S,mom(1))
  !       wLC(ixOmin2^%2ixO^S,mom(1))=wLp(ixOmin2^%2ixO^S,mom(1))*wLp(ixOmin2^%2ixO^S,rho_)
  !       wRC(ixOmin2^%2ixO^S,mom(1))=wRp(ixOmin2^%2ixO^S,mom(1))*wLp(ixOmin2^%2ixO^S,rho_)
  !       wLp(ixOmin2^%2ixO^S,mom(2:3))=0.d0
  !       wRp(ixOmin2^%2ixO^S,mom(2:3))=0.d0
  !       wLC(ixOmin2^%2ixO^S,mom(2:3))=0.d0
  !       wRC(ixOmin2^%2ixO^S,mom(2:3))=0.d0
  !     end if
  !   end if
  !
  ! end subroutine boundary_wLR

  ! subroutine driven_velocity(ixI^L,ixO^L,qt,x,vdr)
  !   integer, intent(in)             :: ixI^L, ixO^L
  !   double precision, intent(in) :: qt, x(ixI^S,1:ndim)
  !   double precision, intent(out) :: vdr(ixI^S)
  !
  !   double precision :: tramp1,ft
  !
  !   tramp1=3.d0
  !   ! if(qt<tramp1) then !MLB
  !   !   ft=qt/tramp1 !MLB
  !   ! else
  !     ft=1.d0
  !   !end if !MLB
  !   where(abs(x(ixO^S,1)+6.5d0)<=La)
  !     vdr(ixO^S)=-ft*vmax*sin(dpi*x(ixO^S,1)/La)
  !   else where
  !     vdr(ixO^S)=0.d0
  !   end where
  !
  ! end subroutine driven_velocity


  subroutine driven_velocity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,x,vdx,vdz)
   !Imposing a velocity on the bottom boundary to simulate convergence and shearing at the PIL.
   ! Here the shearing has been selected to be positive
   ! i.e. to increase the shear along the loops in time.
   use mod_global_parameters
   integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2
   double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
      1:ndim)
   double precision, intent(out) :: vdx(ixImin1:ixImax1,ixImin2:ixImax2),&
      vdz(ixImin1:ixImax1,ixImin2:ixImax2)
   double precision :: v_0t1,tb1,tr1,t1,t2, xcc1, ghost_region, resol, sigmam,&
       conv_region
   double precision :: v_0t1cancellation,tb1cancellation,tr1cancellation,&
      t1cancellation,t2cancellation
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


   xcc1=-6.4651938d0!3.2318599d0 ! deberC-a ser -6.4651938
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
       vdz(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=v_0t1*dexp((-(x(ixOmin1:ixOmax1,&
          ixOmin2:ixOmax2,1)-xcc1)**2)/sigmam**&
          2)*tanh(20.d0*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-xcc1))
       ! vdx(ixO^S)=v_0t1*(-dexp((-(x(ixO^S,1)-xcc1-sigmam)**2)/sigmam**2)+dexp((-(x(ixO^S,1)-xcc1+sigmam)**2)/sigmam**2))   !*dsin(2.d0*dpi*(x(ixO^S,1)-xcc1)/(conv_region))
    end if
!!!cancellation
    if (qt>=tb1cancellation) then
      if (qt>=tb1cancellation.and.qt<tr1cancellation) then
        v_0t1cancellation=((qt-tb1cancellation)/(tr1cancellation-&
           tb1cancellation))*v_00cancellation
      endif
      if (qt>=tr1cancellation.and.qt<t1cancellation) then
        v_0t1cancellation=v_00cancellation
      endif
      if (qt>=t1cancellation.and.qt<t2cancellation) then
        v_0t1cancellation=v_00cancellation*((t2cancellation-&
           qt)/(t2cancellation-t1cancellation))
      endif
      if (qt>=t2cancellation) then
        v_0t1cancellation=0.0d0
      end if

        vdx(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=v_0t1cancellation*(-dexp((-&
           (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)-xcc1-sigmam*2.d0)**2)/(sigmam*2.d0)**2)+&
           dexp((-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
           1)-xcc1+sigmam*2.d0)**2)/(sigmam*2.d0)*2)) !*dsin(2.d0*dpi*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)-xcc1)/(conv_region))
     end if



 end subroutine driven_velocity

  subroutine specialbound_usr(qt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iB,w,x)
    ! special boundary types, user defined
    integer, intent(in) :: ixOmin1,ixOmin2,ixOmax1,ixOmax2, iB, ixImin1,&
       ixImin2,ixImax1,ixImax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision :: Bf(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir) !MLB
    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       tmp(ixImin1:ixImax1,ixImin2:ixImax2),ggrid(ixImin1:ixImax1,&
       ixImin2:ixImax2),invT(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: Q(ixImin1:ixImax1,ixImin2:ixImax2),Qp(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    integer :: ix1,ix2,ixOsmin1,ixOsmin2,ixOsmax1,ixOsmax2,ixCmin1,ixCmin2,&
       ixCmax1,ixCmax2,hxCmin1,hxCmin2,hxCmax1,hxCmax2,jxOmin1,jxOmin2,jxOmax1,&
       jxOmax2,idir

    select case(iB)
    case(3)
      if(iprob<3) then
        !! fixed zero velocity
        do idir=1,ndir
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir))=-w(ixOmin1:ixOmax1,&
             ixOmax2+nghostcells:ixOmax2+1:-1,mom(idir))/w(ixOmin1:ixOmax1,&
             ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        end do
      else
        ! call driven_velocity(ixI^L,ixO^L,qt,x,w(ixI^S,mom(1)))
        ! w(ixO^S,mom(2))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(2))&
        !             /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        ! w(ixO^S,mom(3))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(3))&
        !             /w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
        call driven_velocity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,qt,x,w(ixImin1:ixImax1,ixImin2:ixImax2,mom(1)),&
           w(ixImin1:ixImax1,ixImin2:ixImax2,mom(3)))
        !w(ixO^S,mom(1))=-w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,mom(1))/w(ixOmin1:ixOmax1,ixOmax2+nghostcells:ixOmax2+1:-1,rho_) !MLB
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2))=-w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,mom(2))/w(ixOmin1:ixOmax1,&
           ixOmax2+nghostcells:ixOmax2+1:-1,rho_)
      end if
      !! fixed b1 b2 b3
      if(iprob==0 .or. B0field) then
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(:))=0.d0
      else if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
          ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
          ! 2nd order one-sided zero-gradient extrapolation
          !do ix2=ixOsmax2,ixOsmin2,-1
          !   block%ws(ix2^%2ixOs^S,idir)=1.d0/3.d0*&
          !         (-block%ws(ix2+2^%2ixOs^S,idir)&
          !     +4.d0*block%ws(ix2+1^%2ixOs^S,idir))
          !end do
          ! 4th order one-sided equal-gradient extrapolation
          do ix2=ixOsmax2,ixOsmin2,-1
            block%ws(ixOsmin1:ixOsmax1,ix2,&
               idir)= 0.12d0*block%ws(ixOsmin1:ixOsmax1,ix2+5,&
               idir) -0.76d0*block%ws(ixOsmin1:ixOsmax1,ix2+4,&
               idir) +2.08d0*block%ws(ixOsmin1:ixOsmax1,ix2+3,&
               idir) -3.36d0*block%ws(ixOsmin1:ixOsmax1,ix2+2,&
               idir) +2.92d0*block%ws(ixOsmin1:ixOsmax1,ix2+1,idir)
          end do
        end do
        ixOsmin1=ixOmin1-kr(2,1);ixOsmin2=ixOmin2-kr(2,2)
        ixOsmax1=ixOmax1-kr(2,1);ixOsmax2=ixOmax2-kr(2,2);
        jxOmin1=ixOmin1+nghostcells*kr(2,1)
        jxOmin2=ixOmin2+nghostcells*kr(2,2)
        jxOmax1=ixOmax1+nghostcells*kr(2,1)
        jxOmax2=ixOmax2+nghostcells*kr(2,2);
        block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,2)=zero
        do ix2=ixOsmax2,ixOsmin2,-1
          call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,Qp)
          block%ws(ixOsmin1:ixOsmax1,ix2,2)=Qp(ixOmin1:ixOmax1,&
             ix2+1)*block%dvolume(ixOmin1:ixOmax1,&
             ix2+1)/block%surfaceC(ixOsmin1:ixOsmax1,ix2,2)
        end do
        call mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,block)
        if(iprob<3) then
          !MLB
          ! w(ixO^S,mag(3))=-B0*dcos(kx*x(ixO^S,1))*dexp(-ly*x(ixO^S,2))*dsin(theta)
          call specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,x,Bf) !MLB
          w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(3))=Bf(ixOmin1:ixOmax1,&
             ixOmin2:ixOmax2,3) !MLB
        else
          do ix2=ixOmax2,ixOmin2,-1
            w(ixOmin1:ixOmax1,ix2,mag(3))= 0.12d0*w(ixOmin1:ixOmax1,ix2+5,&
               mag(3)) -0.76d0*w(ixOmin1:ixOmax1,ix2+4,&
               mag(3)) +2.08d0*w(ixOmin1:ixOmax1,ix2+3,&
               mag(3)) -3.36d0*w(ixOmin1:ixOmax1,ix2+2,&
               mag(3)) +2.92d0*w(ixOmin1:ixOmax1,ix2+1,mag(3))
          end do
        end if
      else
        !MLB
        call specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
           ixOmax1,ixOmax2,x,Bf)
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mag(1:ndir))=Bf(ixOmin1:ixOmax1,&
           ixOmin2:ixOmax2,1:ndir) !MLB
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
      call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case(4)
      ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmin2-1;
      call mhd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOsmin1,&
         ixOsmin2,ixOsmax1,ixOsmax2,pth)
      ixOsmin2=ixOmin2-1;ixOsmax2=ixOmax2;
      call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOsmin1,ixOsmin2,&
         ixOsmax1,ixOsmax2,x)
      !> fill pth, rho ghost layers according to gravity stratification
      invT(ixOmin1:ixOmax1,ixOmin2-1)=w(ixOmin1:ixOmax1,ixOmin2-1,&
         rho_)/pth(ixOmin1:ixOmax1,ixOmin2-1)
      tmp=0.d0
      do ix2=ixOmin2,ixOmax2
        tmp(ixOmin1:ixOmax1,ixOmin2-1)=tmp(ixOmin1:ixOmax1,&
           ixOmin2-1)+0.5d0*(ggrid(ixOmin1:ixOmax1,ix2)+ggrid(ixOmin1:ixOmax1,&
           ix2-1))*invT(ixOmin1:ixOmax1,ixOmin2-1)
        w(ixOmin1:ixOmax1,ix2,p_)=pth(ixOmin1:ixOmax1,&
           ixOmin2-1)*dexp(tmp(ixOmin1:ixOmax1,ixOmin2-1)*dxlevel(2))
        w(ixOmin1:ixOmax1,ix2,rho_)=w(ixOmin1:ixOmax1,ix2,&
           p_)*invT(ixOmin1:ixOmax1,ixOmin2-1)
      enddo
      !> fixed zero velocity
      do idir=1,ndir
        w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(idir)) =-w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,mom(idir))/w(ixOmin1:ixOmax1,&
           ixOmin2-1:ixOmin2-nghostcells:-1,rho_)
      end do
      !> zero normal gradient extrapolation
      if(stagger_grid) then
        do idir=1,nws
          if(idir==2) cycle
          ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
          ixOsmin1=ixOmin1-kr(1,idir);ixOsmin2=ixOmin2-kr(2,idir);
          do ix2=ixOsmin2,ixOsmax2
             block%ws(ixOsmin1:ixOsmax1,ix2,&
                idir)=1.d0/3.d0*(-block%ws(ixOsmin1:ixOsmax1,ix2-2,&
                idir)+4.d0*block%ws(ixOsmin1:ixOsmax1,ix2-1,idir))
          end do
        end do
        ixOsmin1=ixOmin1;ixOsmin2=ixOmin2;ixOsmax1=ixOmax1;ixOsmax2=ixOmax2;
        jxOmin1=ixOmin1-nghostcells*kr(2,1)
        jxOmin2=ixOmin2-nghostcells*kr(2,2)
        jxOmax1=ixOmax1-nghostcells*kr(2,1)
        jxOmax2=ixOmax2-nghostcells*kr(2,2);
        block%ws(ixOsmin1:ixOsmax1,ixOsmin2:ixOsmax2,2)=zero
        call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,jxOmin1,jxOmin2,&
           jxOmax1,jxOmax2,Q)
        do ix2=ixOsmin2,ixOsmax2
          call get_divb(w,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
             ixOmax1,ixOmax2,Qp)
          block%ws(ixOsmin1:ixOsmax1,ix2,2)=(Q(jxOmin1:jxOmax1,&
             jxOmax2)*block%dvolume(jxOmin1:jxOmax1,&
             jxOmax2)-Qp(ixOmin1:ixOmax1,ix2)*block%dvolume(ixOmin1:ixOmax1,&
             ix2))/block%surfaceC(ixOsmin1:ixOsmax1,ix2,2)
        end do
        call mhd_face_to_center(ixOmin1,ixOmin2,ixOmax1,ixOmax2,block)
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(3))=(1.0d0/3.0d0)* (-w(ixOmin1:ixOmax1,&
             ix2-2,mag(3))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(3)))
        enddo
      else
        do ix2=ixOmin2,ixOmax2
          w(ixOmin1:ixOmax1,ix2,mag(:))=(1.0d0/3.0d0)* (-w(ixOmin1:ixOmax1,&
             ix2-2,mag(:))+4.0d0*w(ixOmin1:ixOmax1,ix2-1,mag(:)))
        enddo
      end if
      call mhd_to_conserved(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,w,x)
    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  subroutine gravity(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
     ixOmax2,wCT,x,gravity_field)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(in)    :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw)
    double precision, intent(out)   :: gravity_field(ixImin1:ixImax1,&
       ixImin2:ixImax2,ndim)

    double precision                :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    gravity_field=0.d0
    call getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,x)
    gravity_field(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=ggrid(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

  end subroutine gravity

  subroutine getggrav(ggrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x)
    integer, intent(in)             :: ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)    :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision, intent(out)   :: ggrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ggrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=usr_grav*(SRadius/(SRadius+&
       x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))**2
  end subroutine

  subroutine special_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
        wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgrid(ixImin1:ixImax1,ixImin2:ixImax2)

    ! add global background heating bQ
    call getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,qtC,wCT,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       e_)+qdt*bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    if(iprob>1) then !if(iprob==2) then !MLB
      call getlQ(lQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2,qt,wCT,x)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,e_)=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         e_)+qdt*lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    end if

  end subroutine special_source

  subroutine getbQ(bQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
  ! calculate background heating bQ
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: bQgrid(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: bQgridspecial1(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgridspecial2(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: envolvente(ixImin1:ixImax1,ixImin2:ixImax2),&
       aenvolvente
    double precision :: bQgridleft(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgridright(ixImin1:ixImax1,ixImin2:ixImax2),&
       bQgridcenter(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: rraau(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: xcenterleft,xcenterright,ycenterleft,ycenterright
    double precision :: xcentercenter,ycentercenter
    double precision :: rcleft(ixImin1:ixImax1,ixImin2:ixImax2),&
       rcright(ixImin1:ixImax1,ixImin2:ixImax2),rccenter(ixImin1:ixImax1,&
       ixImin2:ixImax2),A(ixImin1:ixImax1,ixImin2:ixImax2)
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



    ! rcleft(ixO^S)=sqrt((x(ixO^S,1)-xcenterleft)**2+(x(ixO^S,2)-ycenterleft)**2)
    ! rcright(ixO^S)=sqrt((x(ixO^S,1)-xcenterright)**2+(x(ixO^S,2)-ycenterright)**2)
    ! rccenter(ixO^S)=sqrt((x(ixO^S,1)-xcentercenter)**2+(x(ixO^S,2)-ycentercenter)**2)
    ! where(rcleft(ixO^S)<radiusthreshold2)
    !      ! bQgrid(ixO^S)=10.d0*bQ0*dexp(-x(ixO^S,2)/(lambdabackground/3.d0))+bQ0*5.0d0!+rcleft(ixO^S)-rcleft(ixO^S)
    !      bQgridleft(ixO^S)=5.d0*bQ0*dexp(-rcleft(ixO^S)/(lambdabackground/3.d0))+bQ0*5.0d0!+rcleft(ixO^S)-rcleft(ixO^S)
    ! end where
    ! where(rcright(ixO^S)<radiusthreshold)
    !      ! bQgrid(ixO^S)=10.d0*bQ0*dexp(-x(ixO^S,2)/(lambdabackground/3.d0))+bQ0*5.0d0!+rcright(ixO^S)-rcright(ixO^S)
    !      bQgridright(ixO^S)=5.d0*bQ0*dexp(-rcright(ixO^S)/(lambdabackground/3.d0))+bQ0*5.0d0!+rcright(ixO^S)-rcright(ixO^S)
    ! end where
    ! where(rccenter(ixO^S)<radiusthreshold3)
    !      bQgridcenter(ixO^S)=25.d0*bQ0*dexp(-rccenter(ixO^S)/(lambdabackground/3.d0))+bQ0*10.0d0!+rcright(ixO^S)-rcright(ixO^S)
    ! end where



      call initvecpot_usr(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, x, A, 3)
      ! bQgridspecial1(ixO^S)=0.d0
      ! bQgridspecial2(ixO^S)=0.d0
      ! !*abs(A(ixO^S))!/minval(A(ixO^S)))!**16.d0
      ! rraau(ixO^S)=sqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)

      ! where(A(ixO^S)>-58.d0 .and. rraau(ixO^S)<1.85) !usamos el codigo IDL para determinar estos valores
            !el factor 4 es por que incrementamos en 4 el valor de B0
            ! bQgrid(ixO^S)=10.d0*bQ0*dexp(-x(ixO^S,2)/lambdabackground/1.d0)
      
      !hay un problema de cizalla entre los loops de la derecha y
      !el canal de filamento. Creo que es que los loops están poco
      !poblados con lo que aumento la extensión del calentamiento (oscar)
      !bQgridspecial1(ixO^S)=15.d0*bQ0*(abs(A(ixO^S))/171.73493d0)**(26.d0)!16.d0*bQ0*(abs(A(ixO^S))/maxval(abs(A(ixO^S))))**8
      
      aenvolvente=-140.+8.
      envolvente=0.5d0*(tanh(-0.7*(A(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-aenvolvente))+1.d0)

      bQgridspecial1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=8.d0*bQ0*(abs(A(&
         ixOmin1:ixOmax1,ixOmin2:ixOmax2))/171.73493d0)**(6.d0)*envolvente
      ! bQgridspecial1(ixO^S)=15.d0*bQ0*(abs(A(ixO^S))/171.73493d0)**(6.d0)
      ! end where
      ! bQgridspecial2(ixO^S)=abs(1.d0*bQ0*((abs(A(ixO^S)/(94.395597d0)))**(-1.d0)))
      !! where(A(ixO^S)>-88.d0 .and. rraau(ixO^S)<5.d0)
      !!   bQgridspecial2(ixO^S)=1.d0*bQ0*(abs(A(ixO^S))/94.395597d0)**(-1.d0)
      !! end where
      !!                     where(A(ixO^S)<Atop .and. A(ixO^S)>Alow)
                          !   lQgrid(ixO^S)=lQ0
                          ! end where
      ! bQgrid(ixO^S)=bQ0*dexp(-x(ixO^S,2)/lambdabackground)+ bQgridspecial(ixO^S)!+bQgridcenter(ixO^S) !he quitado el central
      bQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=bQ0*dexp(-x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,2)/lambdabackground)+ bQgridspecial1(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2) !+bQgridleft(ixOmin1:ixOmax1,ixOmin2:ixOmax2)+bQgridright(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
      !+ bQgridspecial2(ixO^S)!+bQgridcenter(ixO^S) !he quitado el central
      !  bQgrid(ixO^S)=bQgridspecial2(ixO^S)
  end subroutine getbQ

  ! subroutine getlQ(lQgrid,ixI^L,ixO^L,qt,w,x)
  ! ! calculate localized heating lQ
  !   integer, intent(in) :: ixI^L, ixO^L
  !   double precision, intent(in) :: qt, x(ixI^S,1:ndim), w(ixI^S,1:nw)
  !
  !   double precision :: lQgrid(ixI^S), A(ixI^S), lQd, lQ0, hp, lQt,Atop,Alow
  !
  !   lQ0=1.d-2/heatunit
  !   lQt=5.d2/unit_time !MLB 5.d2/unit_time
  !   lQd=0.2d0
  !   hp=0.3d0
  !
  !   call initvecpot_usr(ixI^L, ixO^L, x, A, 3)
  !
  !   lQgrid=0.d0
  !   Atop=B0/kx*cos(kx*1.35d0)
  !   Alow=B0/kx*cos(kx*1.99d0)
  !   !MLB
  !   ! Atop=B0/kx*cos(kx*1.35d0)
  !   ! Alow=B0/kx*cos(kx*1.75d0)
  !   where(A(ixO^S)<Atop .and. A(ixO^S)>Alow)
  !     lQgrid(ixO^S)=lQ0
  !   end where
  !   where(x(ixO^S,2)>hp)
  !     lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(x(ixO^S,2)-hp)**2/lQd)
  !   end where
  !   !MLB
  !   ! if(qt<lQt) then
  !   !   lQgrid(ixO^S)=qt/lQt*lQgrid(ixO^S)
  !   ! endif
  !
  ! end subroutine getlQ

  subroutine getlQ(lQgrid,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,qt,w,x)
  ! calculate localized heating lQ
  !nueva subrutina MLB
    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2
    double precision, intent(in) :: qt, x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim), w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    double precision :: lQgrid(ixImin1:ixImax1,ixImin2:ixImax2),&
        A(ixImin1:ixImax1,ixImin2:ixImax2), lQd, lQ0, hp, lQt,Atop,Alow
    double precision :: auxradius1(ixImin1:ixImax1,ixImin2:ixImax2),xcenter1,&
       ycenter1,auxm1,y11(ixImin1:ixImax1,ixImin2:ixImax2),y12(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: radiuscentral(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: auxradius2(ixImin1:ixImax1,ixImin2:ixImax2),xcenter2,&
       ycenter2,auxm2,y21(ixImin1:ixImax1,ixImin2:ixImax2),y22(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: lambdamin,lambdamax,lambdamax2,lambda1(ixImin1:ixImax1,&
       ixImin2:ixImax2),lambda11(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: lambda2(ixImin1:ixImax1,ixImin2:ixImax2),&
       lambda22(ixImin1:ixImax1,ixImin2:ixImax2),xprimaL(ixImin1:ixImax1,&
       ixImin2:ixImax2),zprimaL(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: xprimaR(ixImin1:ixImax1,ixImin2:ixImax2),&
       zprimaR(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: thetaangle01,thetaangle1(ixImin1:ixImax1,&
       ixImin2:ixImax2),thetaangle02,thetaangle2(ixImin1:ixImax1,&
       ixImin2:ixImax2)
    double precision :: lQgridleft(ixImin1:ixImax1,ixImin2:ixImax2),&
       lQgridright(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: recta1(ixImin1:ixImax1,ixImin2:ixImax2),&
       recta2(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: sigma_heating,sigma_heating_impulsive,y1,lam1
    double precision :: impulsiveheatingfactor,impulsiveheatingduration,&
       timeimpulse
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !calentamiento en el pie izquierdo
    !xcenter1=-5.2288454d0!+0.d0!-8.4243435d0+0.3d0!-8.1d0!-3.935d0*1.d0 !lo he puesto para que coincidad con el otro extremo
!    ycenter1=hp*0.d0
!    auxm1=dtan(52.006568d0*dpi/180.d0)
!    thetaangle01=52.006568d0*dpi/180.d0
!    thetaangle1(ixO^S)=atan((x(ixO^S,2)-ycenter1)/(x(ixO^S,1)-xcenter1))
!    lambda1(ixO^S)=lambdamin+(lambdamax-lambdamin)*cos(thetaangle1(ixO^S)-thetaangle01)**2
!    lambda11(ixO^S)=lambdamin+(lambdamax2-lambdamin)*cos(thetaangle1(ixO^S)-thetaangle01)**2
!    auxradius1(ixO^S)=sqrt((x(ixO^S,1)-xcenter1)**2+(x(ixO^S,2)-ycenter1)**2)
    ! recta1(ixO^S)=-1.d0/auxm1*(x(ixO^S,1)-xcenter1)+ycenter1
    ! call initvecpot_usr(ixI^L, ixO^L, x, A, 3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !izquierda
    !rotamos las coordenadas
!        xprimaL(ixO^S)=(x(ixO^S,1)-xcenter1)*cos(thetaangle01)+&
!                      (x(ixO^S,2)-ycenter1)*sin(thetaangle01)
!        zprimaL(ixO^S)=-(x(ixO^S,1)-xcenter1)*sin(thetaangle01)+&
!                      (x(ixO^S,2)-ycenter1)*cos(thetaangle01)
!        lQgridleft(ixO^S)=1.d0*lQ0*dexp(-(xprimaL(ixO^S)**2/lambdamax**2+&
!                          zprimaL(ixO^S)**2/lambdamin**2))
        ! where(abs(thetaangle1(ixO^S)-thetaangle01)<dpi/2.d0)
        !     lQgridleft(ixO^S)=1.d0*lQ0*dexp(-(auxradius1(ixO^S))/lambda1(ixO^S))
        ! end where
        ! where(abs(thetaangle1(ixO^S)-thetaangle01)>dpi/2.d0)
        !   lQgridleft(ixO^S)=  1.d0*lQ0*dexp(-(auxradius1(ixO^S))/lambda11(ixO^S))
        ! end where
        ! where(x(ixO^S,2)>recta1(ixO^S))
        !   lQgridleft(ixO^S)=1.d0*lQ0*dexp(-(auxradius1(ixO^S))/lambda1(ixO^S))
        ! end where
        ! where(x(ixO^S,2)<recta1(ixO^S))
        !    lQgridleft(ixO^S)=1.d0*lQ0*dexp(-(auxradius1(ixO^S))/lambda11(ixO^S))
        ! end where
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !
    ! y11(ixO^S)=ycenter1+auxm1*(x(ixO^S,1)+3.7d0*1.d0) !volver a poner 4.1? Mira el termino derecho !pongo el soplido en 3.3*2
    ! y12(ixO^S)=ycenter1+auxm1*(x(ixO^S,1)+4.17d0*1.d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !calentamiento en el pie derecho
    !xcenter2=5.2288454d0!-0.d0!8.1d0!
    ! ycenter2=hp*0.d0
    ! auxm2=dtan(127.99343d0*dpi/180.d0)
    ! thetaangle02=127.99343d0*dpi/180.d0
    ! thetaangle2(ixO^S)=atan((x(ixO^S,2)-ycenter2)/(x(ixO^S,1)-xcenter2))
    ! lambda2(ixO^S)=lambdamin+(lambdamax-lambdamin)*cos(thetaangle2(ixO^S)-thetaangle02)**2
    ! lambda22(ixO^S)=lambdamin+(lambdamax2-lambdamin)*cos(thetaangle2(ixO^S)-thetaangle02)**2
    ! ! recta2(ixO^S)=-1.d0/auxm2*(x(ixO^S,1)-xcenter2)+ycenter2

    ! auxradius2(ixO^S)=sqrt((x(ixO^S,1)-xcenter2)**2+(x(ixO^S,2)-ycenter2)**2)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     !derecha
    !     xprimaR(ixO^S)=(x(ixO^S,1)-xcenter2)*cos(thetaangle02)+&
    !                   (x(ixO^S,2)-ycenter2)*sin(thetaangle02)
    !     zprimaR(ixO^S)=-(x(ixO^S,1)-xcenter2)*sin(thetaangle02)+&
    !                   (x(ixO^S,2)-ycenter2)*cos(thetaangle02)
    !     lQgridright(ixO^S)=1.d0*lQ0*dexp(-(xprimaR(ixO^S)**2/lambdamax**2+&
    !                       zprimaR(ixO^S)**2/lambdamin**2))
        ! where(abs(thetaangle2(ixO^S)-thetaangle02)<dpi/2.d0)
        !   lQgridright(ixO^S)=1.d0*lQ0*dexp(-(auxradius2(ixO^S))/lambda2(ixO^S))
        ! end where
        ! where(abs(thetaangle2(ixO^S)-thetaangle02)>dpi/2.d0)
        !   lQgridright(ixO^S)=1.d0*lQ0*dexp(-(auxradius2(ixO^S))/lambda22(ixO^S))
        ! end where
      ! end where
        ! where(x(ixO^S,2)>recta2(ixO^S))
        ! lQgridright(ixO^S)=1.d0*lQ0*dexp(-(auxradius2(ixO^S))/lambda2(ixO^S))
        ! end where
        ! where(x(ixO^S,2)<recta2(ixO^S))
        ! lQgridright(ixO^S)=1.d0*lQ0*dexp(-(auxradius2(ixO^S))/lambda22(ixO^S))
        ! end where
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y21(ixO^S)=ycenter2+auxm2*(x(ixO^S,1)-3.8d0*1.d0)
    ! y22(ixO^S)=ycenter2+auxm2*(x(ixO^S,1)-4.1d0*1.d0)





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calentamiento total
     ! lQgrid(ixO^S)=(lQgridleft(ixO^S)+lQgridright(ixO^S))*4.d0 !el factor 2 corrige
     !que el calentamiento esté más profundo.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!












    ! Atop=B0/kx*cos(kx*1.35d0)
    ! Alow=B0/kx*cos(kx*1.99d0)
    !MLB
    ! Atop=B0/kx*cos(kx*1.35d0)
    ! Alow=B0/kx*cos(kx*1.75d0)
    ! where(A(ixO^S)<Atop .and. A(ixO^S)>Alow)
    !   lQgrid(ixO^S)=lQ0
    ! end where
    ! where(x(ixO^S,2)>hp)
    !   lQgrid(ixO^S)=lQgrid(ixO^S)*dexp(-(x(ixO^S,2)-hp)**2/lQd)
    ! end where
    !MLB
    ! if(qt<lQt) then
    !   lQgrid(ixO^S)=qt/lQt*lQgrid(ixO^S)
    ! endif

! quito el calentamiento cuadrado en la cromosfera
    ! where(x(ixO^S,2)<hp .and. x(ixO^S,1)<=-3.8d0 .and. x(ixO^S,1)>=-4.2d0)
    !   lQgrid(ixO^S)=lQ0
    ! end where


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! where(x(ixO^S,2)>y11(ixO^S) .and. x(ixO^S,2)<y12(ixO^S) .and. x(ixO^S,2)>hp)
      ! lQgrid(ixO^S)=1.d0*lQ0*dexp(-(auxradius1(ixO^S))/lQd)
      ! lQgrid(ixO^S)=1.d0*lQ0*dexp(-(auxradius1(ixO^S))/lambda1(ixO^S))
    ! end where
    ! where(x(ixO^S,2)>y11(ixO^S) .and. x(ixO^S,2)<y12(ixO^S) .and. x(ixO^S,2)<hp)
    !   lQgrid(ixO^S)=1.d0*lQ0
    ! end where
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! where(x(ixO^S,2)<hp .and. x(ixO^S,1)<=4.3d0 .and. x(ixO^S,1)>=3.9d0)
    !   lQgrid(ixO^S)=lQ0
    ! end where


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !derecha
    ! where(x(ixO^S,2)>y21(ixO^S) .and. x(ixO^S,2)<y22(ixO^S) .and. x(ixO^S,2)>hp)
      ! lQgrid(ixO^S)=2.d0*lQ0*dexp(-(auxradius2(ixO^S))/lQd)
      ! lQgrid(ixO^S)=1.d0*lQ0*dexp(-(auxradius2(ixO^S))/lambda2(ixO^S))
    ! end where
    ! where(x(ixO^S,2)>y21(ixO^S) .and. x(ixO^S,2)<y22(ixO^S) .and. x(ixO^S,2)<hp)
    !   lQgrid(ixO^S)=2.d0*lQ0
    ! end where
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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
      impulsiveheatingduration=1.d0 ! every unit is 85 secs
      timeimpulse=75.001d0
      impulsiveheatingfactoramplitude=20.d0
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
      !auxradius1(ixO^S)=sqrt((x(ixO^S,1)-xcenter1)**2+(x(ixO^S,2)-y1)**2)

      !!CUIDADO HAY UNA ASIMETRÍA EN EL HEATING DEL 0.75 POR UNO.
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=lQ0*(impulsiveheatingfactor*dexp(&
         -(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)-xcenter1)**2/sigma_heating**2)+1.d0*dexp(-(x(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,1)-xcenter2)**2/sigma_heating**2) )
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=lQgrid(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*dexp(-(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)-y1)**2/lam1**2)

      call initvecpot_usr(ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
         ixOmax1,ixOmax2, x, A, 3)
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=lQgrid(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)*dexp(-(A(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2)-A0)**2/2.d0/(0.50d0)**2)

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
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=0.d0
    endif
    if(qt>10.d0 .and. qt<60.d0) then
      lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(qt-&
         10.d0)/50.d0*lQgrid(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    endif
    !solo se calienta un fajo de lineas
!    call initvecpot_usr(ixI^L, ixO^L, x, A, 3)
!    radiuscentral(ixO^S)=sqrt(x(ixO^S,1)**2+x(ixO^S,2)**2)
!    where(A(ixO^S)<-58.0d0 .or. A(ixO^S)>-57.0d0 .or. radiuscentral(ixO^S)<1.6) !usamos el codigo IDL para determinar estos valores
!      lQgrid(ixO^S)=0.d0
!    end where

    !!!if(qt>lQt) then
    !!!   lQgrid(ixO^S)=qt/lQt*lQgrid(ixO^S)
    !!!endif


  end subroutine getlQ

  subroutine special_refine_grid(igrid,level,ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2,qt,w,x,refine,coarsen)
  ! Enforce additional refinement or coarsening
  ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    integer, intent(in) :: igrid, level, ixImin1,ixImin2,ixImax1,ixImax2,&
        ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in) :: qt, w(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    integer, intent(inout) :: refine, coarsen

    ! fix the bottom layer to the highest level
    if (any(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)<=xprobmin2+0.05d0)) then
      refine=1
      coarsen=-1
    endif

  end subroutine special_refine_grid

  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
  ! this subroutine can be used in convert, to add auxiliary variables to the
  ! converted output file, for further analysis using tecplot, paraview, ....
  ! these auxiliary values need to be stored in the nw+1:nw+nwauxio slots
  ! the array normconv can be filled in the (nw+1:nw+nwauxio) range with
  ! corresponding normalization values (default value 1)
    use mod_radiative_cooling
    integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
       ixOmin1,ixOmin2,ixOmax1,ixOmax2
    double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim)
    double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       nw+nwauxio)
    double precision                   :: normconv(0:nw+nwauxio)

    double precision :: pth(ixImin1:ixImax1,ixImin2:ixImax2),&
       B2(ixImin1:ixImax1,ixImin2:ixImax2),tmp2(ixImin1:ixImax1,&
       ixImin2:ixImax2),dRdT(ixImin1:ixImax1,ixImin2:ixImax2)
    double precision :: ens(ixImin1:ixImax1,ixImin2:ixImax2),&
       divb(ixImin1:ixImax1,ixImin2:ixImax2),wlocal(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    double precision :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
       curlvec(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
    double precision :: Te(ixImin1:ixImax1,ixImin2:ixImax2),tco_local
    double precision, dimension(ixImin1:ixImax1,ixImin2:ixImax2,&
       1:ndim) :: gradT, bunitvec
    integer :: idirmin,idir,ix1,ix2
    logical :: lrlt(ixImin1:ixImax1,ixImin2:ixImax2)

    wlocal(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)=w(ixImin1:ixImax1,&
       ixImin2:ixImax2,1:nw)
    ! output temperature
    call mhd_get_pthermal(wlocal,x,ixImin1,ixImin2,ixImax1,ixImax2,ixImin1,&
       ixImin2,ixImax1,ixImax2,pth)
    Te(ixImin1:ixImax1,ixImin2:ixImax2)=pth(ixImin1:ixImax1,&
       ixImin2:ixImax2)/w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1)=Te(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

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
    ! B^2
    B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=sum((Btotal(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,:))**2,dim=ndim+1)

    ! output Alfven wave speed B/sqrt(rho)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2)=dsqrt(B2(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_))

    ! output divB1
    call get_divb(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,divb)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)=divb(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! output the plasma beta p*2/B**2
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4)=pth(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)*two/B2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
    ! output heating rate
    call getbQ(ens,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,global_time,wlocal,x)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)
    ! store the cooling rate
    if(mhd_radiative_cooling)call getvar_cooling(ixImin1,ixImin2,ixImax1,&
       ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,wlocal,x,ens,rc_fl)
    w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6)=ens(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2)

    ! store current
    call get_current(wlocal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,idirmin,curlvec)
    do idir=1,ndir
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6+idir)=curlvec(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,idir)
    end do

    if(nw_extra>0) w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       nw+10)=block%wextra(ixOmin1:ixOmax1,ixOmin2:ixOmax2,iw_tcoff)

    end subroutine specialvar_output

      subroutine specialvar_output_custom(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2,w,x,normconv)
      integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,ixOmax2
      double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
      double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,nw+nwauxio)
      double precision                   :: normconv(0:nw+nwauxio)
      double precision :: Btotal(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)
      integer :: idir

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

      print*, 'specialvar_output_custom'
      print*, 'B0:', Btotal(1,1,1), Btotal(1,1,2), Btotal(1,1,3)
      print*, 'B_pert:', w(1,1,mag(1)), w(1,1,mag(2)), w(1,1,mag(3))
      print*, 'B_tot:', w(1,1,bt_(1)), w(1,1,bt_(2)), w(1,1,bt_(3))
      print*,''

      end subroutine specialvar_output_custom

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

  subroutine specialvarnames_output(varnames)
  ! newly added variables need to be concatenated with the w_names/primnames string
       character(len=*) :: varnames
       varnames = 'bt1 bt2 bt3'
  end subroutine specialvarnames_output

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

  subroutine specialset_B0(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x,wB0)
 ! Here add a time-independent background magnetic field
   use mod_global_parameters
   integer, intent(in)           :: ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
      ixOmin2,ixOmax1,ixOmax2
   double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
   double precision, intent(inout) :: wB0(ixImin1:ixImax1,ixImin2:ixImax2,&
      1:ndir)
   double precision :: Bx1(ixImin1:ixImax1,ixImin2:ixImax2),&
      By1(ixImin1:ixImax1,ixImin2:ixImax2),Bx2(ixImin1:ixImax1,&
      ixImin2:ixImax2),By2(ixImin1:ixImax1,ixImin2:ixImax2),&
      Bx3(ixImin1:ixImax1,ixImin2:ixImax2),By3(ixImin1:ixImax1,&
      ixImin2:ixImax2)
   double precision :: Bradial1(ixImin1:ixImax1,ixImin2:ixImax2),&
      Btheta1(ixImin1:ixImax1,ixImin2:ixImax2),Bradial2(ixImin1:ixImax1,&
      ixImin2:ixImax2),Btheta2(ixImin1:ixImax1,ixImin2:ixImax2)
   double precision :: rho1(ixImin1:ixImax1,ixImin2:ixImax2),&
      cosphi1(ixImin1:ixImax1,ixImin2:ixImax2),sinphi1(ixImin1:ixImax1,&
      ixImin2:ixImax2),cos2phi1(ixImin1:ixImax1,ixImin2:ixImax2)
   double precision :: rho2(ixImin1:ixImax1,ixImin2:ixImax2),&
      cosphi2(ixImin1:ixImax1,ixImin2:ixImax2),sinphi2(ixImin1:ixImax1,&
      ixImin2:ixImax2),cos2phi2(ixImin1:ixImax1,ixImin2:ixImax2)
   ! double precision :: Bx1(ixI^S),By1(ixI^S)
   double precision :: xn1(ixImin1:ixImax1,ixImin2:ixImax2),&
      yn1(ixImin1:ixImax1,ixImin2:ixImax2),xn2(ixImin1:ixImax1,&
      ixImin2:ixImax2),yn2(ixImin1:ixImax1,ixImin2:ixImax2),&
      xn3(ixImin1:ixImax1,ixImin2:ixImax2),yn3(ixImin1:ixImax1,&
      ixImin2:ixImax2)

   xn1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      1)-xc1)
   yn1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      2)+hb)

   xn2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      1)-xc2)
   yn2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      2)+hb2)

   xn3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      1)-xc3)
   yn3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      2)+hb3)

   Bx1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-B0*((adipo**2-hb**2)*(adipo**2+&
      xn1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2-yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2))/((adipo**2+xn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-2.d0*adipo*yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))*(adipo**2+xn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+2.d0*adipo*yn1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
   By1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-B0*(2.d0*(adipo**2-&
      hb**2)*xn1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))/((adipo**2+xn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-2.d0*adipo*yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))*(adipo**2+xn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+2.d0*adipo*yn1(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

   Bx2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-B0*((adipo**2-hb2**2)*(adipo**2+&
      xn2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)**2-yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2))/((adipo**2+xn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-2.d0*adipo*yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))*(adipo**2+xn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+2.d0*adipo*yn2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
   By2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-B0*(2.d0*(adipo**2-&
      hb2**2)*xn2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)*yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))/((adipo**2+xn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-2.d0*adipo*yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))*(adipo**2+xn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+2.d0*adipo*yn2(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

   Bx3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-B0*factorcentralarcade*((adipo**2-&
      hb3**2)*(adipo**2+xn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2))/((adipo**2+xn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-2.d0*adipo*yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))*(adipo**2+xn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+2.d0*adipo*yn3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))
   By3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)=-&
      B0*factorcentralarcade*(2.d0*(adipo**2-hb3**2)*xn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)*yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))/((adipo**2+xn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2-2.d0*adipo*yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2))*(adipo**2+xn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+yn3(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)**2+2.d0*adipo*yn3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)))

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

   wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)=Bx1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+Bx2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+Bx3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)=By1(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+By2(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)+By3(ixOmin1:ixOmax1,ixOmin2:ixOmax2)
   wB0(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)=0.0 !*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
 end subroutine specialset_B0







  subroutine specialthreshold(wlocal,xlocal,tolerance,qt,level)
    !PURPOSE: use different tolerance in special regions for AMR to
    !reduce/increase resolution there where nothing/something interesting happens.
    use mod_global_parameters

    double precision, intent(in) :: wlocal(1:nw),xlocal(1:ndim),qt
    double precision, intent(inout) :: tolerance
    integer, intent(in) :: level

    double precision :: bczone1,bczone2,addtol,tol_add

    tol_add=0.d0
    !amplitude of additional tolerance
    addtol=0.6d0
    ! thickness of near-boundary region
    bczone1=0.8d0
    bczone2=2.d0
    ! linear changing of additional tolerance
    if(xlocal(1)-xprobmin1 < bczone1 .or. xprobmax1-xlocal(1) < bczone1) then
      tol_add=(1.d0-min(xlocal(1)-xprobmin1,&
         xprobmax1-xlocal(1))/bczone1)*addtol
    endif
    if(xprobmax2-xlocal(2) < bczone2) then
      tol_add=(1.d0-(xprobmax2-xlocal(2))/bczone2)*addtol
    endif
    tolerance=tolerance+tol_add

  end subroutine specialthreshold

end module mod_usr
