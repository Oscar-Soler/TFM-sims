module mod_usr
  use mod_mhd

  implicit none

contains

  subroutine usr_init()

      usr_init_one_grid => rm2d_init_one_grid

      ! usr_create_particles => generate_particles

      call set_coordinate_system("Cartesian_2.5D")


      ! unit_length = 1d8  ! cm
      ! unit_numberdensity = 1d9  ! cm ^-3
      ! unit_temperature = 974204.9765757825 ! K
      ! unit_velocity = 1.5d7
      
      call mhd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm2d_init_one_grid(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,w,x)
      integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2
      double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
      double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
      double precision :: M, rho_eq, v_0, v_eq, p_eq, rho_1, v_1, v_sh, p_1, M1,&
      cs_0, gamma
    ! mhd variables
    ! for a uniform magnetic field
      double precision :: theta_b, phi_b, B0, b_0x, b_0y, b_0z, beta_0,&
      B1, b_1x, b_1y, b_1z, vA0, vA1, X_value, mu_0=1,&
      v_0x, v_0y, v_0z, v_1x, v_1y, v_1z
    ! mhd and perpendicular / oblique shock switch
      logical :: magnetic, perpendicular
    ! thread variables
      double precision :: q_th, rho_th, x_th_ini, x_th_fin, x_sh, trans_sh, y_th_ini, y_th_fin

      perpendicular = .true.
      magnetic = .true.
      M = 5.405405405405405d0
      rho_eq = 1d0
      p_eq = 1d0

      if (perpendicular .eqv. .true.) then
         theta_b = 0
         phi_b = 0
      end if
      
      beta_0 = 1
      X_value = 3.14954219493607

    q_th = 100d0
      x_th_ini = 2
      x_th_fin = 3
      y_th_ini = -0.5
      y_th_fin = 0.5
      
      x_sh = 0

      trans_sh = 0 ! thickness of shock front
      
    !------------------------------------- Above parameters, below calculations
       !if (magnetic .eqv. .false.) gamma = hd_gamma
      if (magnetic .eqv. .true.) gamma = mhd_gamma
      
    p_1 = p_eq * ((2*gamma)/(gamma+1)*M**2 - (gamma-1)/(gamma+1))
    rho_1 = rho_eq / (2/((gamma+1)*M**2) + (gamma-1)/(gamma+1))
    M1 = sqrt( (2+(gamma-1)*M**2) / (2*gamma*M**2 - (gamma-1)))
    cs_0 = sqrt(gamma * p_eq / rho_eq)
    !v_eq = sqrt(hd_gamma * p_eq / rho_eq)
    !v_0 = -M*cs_0
      v_0 = 0
      v_0x = 0
      v_sh = v_0 + M*cs_0
      v_1 = rho_eq/rho_1 * (v_0-v_sh) + v_sh
      rho_th = rho_eq * q_th

      ! mhd
      if (magnetic .eqv. .true.) then
         B0 = SQRT(2*mu_0*p_eq/beta_0)
         vA0 =  B0/ SQRT(mu_0 * rho_eq)
         if (perpendicular .eqv. .true.) then
            b_0z = B0
            b_0x = 0
            b_0y = 0
            rho_1 = rho_eq * X_value
            v_1x = (v_0x-v_sh) / X_value + v_sh
            v_1y = 0
            !v2y = v1y * (v1**2-vA1**2)/(v1**2-X_sol*vA1**2)
            v_1 = SQRT(v_1x**2 + v_1y**2)
            b_1x = b_0x
            b_1y = b_0y
            b_1z = b_0z * X_value
            p_1 = p_eq*(1+gamma*M**2*(1-1/X_value)+1/beta_0*(1-X_value**2))
         end if
      end if
      
    where (x(ixmin1:ixmax1,ixmin2:ixmax2,1) < x_sh)
        w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   = rho_1
        w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = v_1
        w(ixmin1:ixmax1,ixmin2:ixmax2,p_)     = p_1
    elsewhere ((x(ixmin1:ixmax1,ixmin2:ixmax2,1) > x_sh) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,1) < x_sh + trans_sh))
        w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   = (rho_eq-rho_1)/trans_sh * (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-x_sh) + rho_1
        w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = (v_0-v_1)/trans_sh * (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-x_sh) + v_1
        w(ixmin1:ixmax1,ixmin2:ixmax2,p_)     = (p_eq-p_1)/trans_sh * (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-x_sh) + p_1
    elsewhere ((x(ixmin1:ixmax1,ixmin2:ixmax2,1) > x_th_ini) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,1) < x_th_fin) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,2) > y_th_ini) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,2) < y_th_fin))
        w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   = rho_th
        w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = v_0
        w(ixmin1:ixmax1,ixmin2:ixmax2,p_)     = p_eq
    elsewhere
        w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   = rho_eq
        w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) = v_0
        w(ixmin1:ixmax1,ixmin2:ixmax2,p_)     = p_eq
      end where

      if (magnetic .eqv. .true.) then
         where (x(ixmin1:ixmax1,ixmin2:ixmax2,1) < x_sh)
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(1))   = b_1x
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(2))   = b_1y
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(3))   = b_1z
         elsewhere ((x(ixmin1:ixmax1,ixmin2:ixmax2,1) > x_sh) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,1) < x_sh + trans_sh))
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(1))   = (b_0x-b_1x)/trans_sh * (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-x_sh) + b_1x
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(2))   = (b_0y-b_1y)/trans_sh * (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-x_sh) + b_1y
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(3))   = (b_0z-b_1z)/trans_sh * (x(ixmin1:ixmax1,ixmin2:ixmax2,1)-x_sh) + b_1z
         elsewhere ((x(ixmin1:ixmax1,ixmin2:ixmax2,1) > x_th_ini) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,1) < x_th_fin) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,2) > y_th_ini) .and. (x(ixmin1:ixmax1,ixmin2:ixmax2,2) < y_th_fin))
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(1))   = b_th_x
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(2))   = b_th_y
            w(ixmin1:ixmax1,ixmin2:ixmax2,mag(3))   = b_th_z
         end where
      end if
    
      call hd_to_conserved(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,ixmax2,w,x)

      end subroutine rm2d_init_one_grid

subroutine generate_particles(n_particles, x, v, q, m, follow)
  integer, intent(in)           :: n_particles
  double precision, intent(out) :: x(3, n_particles)
  double precision, intent(out) :: v(3, n_particles)
  double precision, intent(out) :: q(n_particles)
  double precision, intent(out) :: m(n_particles)
  logical, intent(out)          :: follow(n_particles)
  integer :: i

      do i = 1, n_particles
         follow(i) = .true.
      end do

      do i= 1, n_particles
         x(1,i) = -4 + 2*i
      end do
end subroutine generate_particles

end module mod_usr
