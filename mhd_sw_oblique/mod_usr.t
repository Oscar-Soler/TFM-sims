! Initial settings for MHD shocks
module mod_usr
  use mod_mhd

  implicit none

contains

  subroutine usr_init()
    usr_set_parameters => initglobaldata_usr
    usr_init_one_grid => rm1d_init_one_grid

    call set_coordinate_system("Cartesian_1.75D")
    call mhd_activate()

  end subroutine usr_init

  !Initialize gamma value
  subroutine initglobaldata_usr
    mhd_gamma = 5.0d0/3.0d0
  end subroutine initglobaldata_usr

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixG^L,ix^L,w,x)
    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision :: pi, mu_0
    double precision :: M0, beta_0, theta_0
    double precision :: rho_0, v_0, v_0x, v_0y, p_0, b_0, b_0x, b_0y, M0_ms, cs_0, csA_0, cms_0 
    double precision :: rho_1, v_1, v_1x, v_1y, p_1, b_1, b_1x, b_1y
    double precision :: u_0, u_1, v_shock
    double precision :: x_value

    pi = 3.14159265358979323
    mu_0 = 1 ! 4*pid-7

!----------------------Define parameters
    beta_0 = 1d-1
    theta_0 = pi/4 ! attack angle (below horizontal, upwards)
    theta_0 = -theta_0 + pi

    rho_0 = 0.8d0
    p_0 = 1.15d0
    v_0x = -3
    x_value = 2.596651338358

!--------------------- Below here just calculations, all parameteres set above
    cs_0 = sqrt(mhd_gamma*p_0/rho_0)    !Sound speed

    b_0 = sqrt(2*mu_0*p_0/beta_0)
    b_0x = b_0 * cos(theta_0)
    b_0y = b_0 * sin(theta_0)
    csA_0 = b_0 / sqrt(mu_0*rho_0)       !Alfven speed
    cms_0 = sqrt(cs_0**2 + csA_0**2)    !magnetosonic speed

!comoving viewer so that the fluid speed before and after shock
! is parallel to magnetic field
    v_0y = v_0x * b_0y/b_0x
    v_0 = sqrt(v_0x**2 + v_0y**2)

    rho_1 = rho_0 * x_value
    b_1x = b_0x
    v_1x = v_0x / x_value
    v_1y = v_0y * (v_0**2-csA_0**2)/(v_0**2-x_value*csA_0**2)
    v_1 = sqrt(v_1x**2+v_1y**2)
    b_1y = b_0y * (v_0**2-csA_0**2)*x_value/(v_0**2-x_value*csA_0**2)
    p_1 = p_0 * (x_value + (mhd_gamma-1)*x_value*v_0**2/(2*cs_0**2)*(1-v_1**2/v_0**2))


    where (x(ix^S,1) < 5.0d0)
        w(ix^S,rho_)   = rho_1
        w(ix^S,mom(1)) = v_1x
        w(ix^S,mom(2)) = v_1y
        w(ix^S,mom(3)) = 0
        w(ix^S,p_)     = p_1
        w(ix^S,mag(1)) = b_1x
        w(ix^S,mag(2)) = b_1y
        w(ix^S,mag(3)) = 0

    elsewhere
        w(ix^S,rho_)   = rho_0
        w(ix^S,mom(1)) = v_0x
        w(ix^S,mom(2)) = v_0y
        w(ix^S,mom(3)) = 0
        w(ix^S,p_)     = p_0
        w(ix^S,mag(1)) = b_0x
        w(ix^S,mag(2)) = b_0y
        w(ix^S,mag(3)) = 0
    end where


    call mhd_to_conserved(ixG^L,ix^L,w,x)

  end subroutine rm1d_init_one_grid

end module mod_usr
