module mod_usr
  use mod_mhd

  implicit none

contains

  subroutine usr_init()

      usr_init_one_grid => rm1d_init_one_grid

      !usr_create_particles => generate_particles

    call set_coordinate_system("Cartesian_1.75D")
    call mhd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
    ! hd variables for shocks  
    double precision :: M, rho_eq, v_0, v_eq, p_eq, rho_1, v_1, v_sh, p_1, M1,&
      cs_0, gamma
    double precision :: pi = ACOS(-1d0)
    ! mhd variables
    ! for a uniform magnetic field
      double precision :: theta_b, phi_b, B0, b_0x, b_0y, b_0z, beta_0,&
      B1, b_1x, b_1y, b_1z, vA0, vA1, X_value, mu_0=1,&
      v_0x, v_0y, v_0z, v_1x, v_1y, v_1z, V_y_ref = 0, u_0, u_1
    ! mhd and perpendicular / oblique shock switch
    logical :: magnetic, perpendicular, parallel, oblique  
    ! thread variables
    double precision :: q_th, rho_th, x_th_ini, x_th_fin, x_sh

      perpendicular = .false.
      parallel = .true.
      oblique = .false.
      magnetic = .true.
      M = 5.405405405405405d0
      rho_eq = 1d0
      p_eq = 1d0
      
      if (perpendicular .eqv. .true.) then
         !theta_b = pi/2
         phi_b = pi/2
      end if
      
      beta_0 = 1d10
      X_value = 3.14954219493607

      v_0x = 0
      v_0y = 0

      q_th = 100d0
      x_th_ini = 2
      x_th_fin = 3
      x_sh = 0
      
      
!------------------------------------- Above parameters, below calculations
      !if (magnetic .eqv. .false.) gamma = hd_gamma
      if (magnetic .eqv. .true.) gamma = mhd_gamma

      cs_0 = sqrt(gamma * p_eq / rho_eq)
      v_sh = v_0 + M*cs_0
      v_1y = v_0y

      if (parallel .eqv. .true.) then
         p_1 = p_eq * ((2*gamma)/(gamma+1)*M**2 - (gamma-1)/(gamma+1))
         rho_1 = rho_eq / (2/((gamma+1)*M**2) + (gamma-1)/(gamma+1))
!     M1 = sqrt( (2+(gamma-1)*M**2) / (2*gamma*M**2 - (gamma-1)))
!     v_eq = sqrt(gamma * p_eq / rho_eq)
!     v_0 = -M*cs_0
!         v_0 = 0
         v_1x = rho_eq/rho_1 * (v_0-v_sh) + v_sh
      end if

      rho_th = rho_eq * q_th

    ! mhd
      if (magnetic .eqv. .true.) then
         B0 = SQRT(2*mu_0*p_eq/beta_0)
         vA0 =  B0/ SQRT(mu_0 * rho_eq)

         if ((perpendicular .eqv. .true.) .or. (oblique .eqv. .true.)) then
            b_0z = B0* COS(theta_b)
            b_0x = B0* SIN(theta_b)*COS(phi_b)
            b_0y = B0* SIN(theta_b)*SIN(phi_b)
         end if
         

         if (parallel .eqv. .true.) then
            b_0x = B0
            b_0y = 0
            b_0z = 0
            b_1x = B0
            b_1y = 0
            b_1z = 0
         end if
         if (oblique .eqv. .true.)  V_y_ref = v_0y - (v_0x - v_sh) * b_0x/b_0y
         u_0 = SQRT( (v_0x-v_sh)**2 + (v_0y - V_y_ref)**2)
         v_1y = (v_0y-V_y_ref) * (u_0**2-vA0**2)/(u_0**2-X_value*vA0**2) + V_y_ref

         if ((perpendicular .eqv. .true.) .or. (oblique .eqv. .true.))then
            rho_1 = rho_eq * X_value
            v_1x = (v_0x-v_sh) / X_value + v_sh
            
            b_1x = b_0x
            if (perpendicular .eqv. .true.) then
               b_1y = b_0y * X_value
               b_1z = b_0z * X_value
               p_1 = p_eq*(1+gamma*M**2*(1-1/X_value)+1/beta_0*(1-X_value**2))
            end if
            if (oblique .eqv. .true.) then
               u_1 =  SQRT( (v_1x-v_sh)**2 + (v_1y - V_y_ref)**2)
               b_1y = b_0y * (u_0**2-vA0**2)*X_value/(u_0**2-X_value*vA0**2)
               p_1 = p_eq * (X_value + (gamma-1)*X_value*u_0**2/(2*cs_0**2)*(1-u_1**2/u_0**2))
            end if
            
         end if
      end if
            
      where (x(ixmin1:ixmax1,1) < x_sh)
         w(ixmin1:ixmax1,rho_)   = rho_1
         w(ixmin1:ixmax1,mom(1)) = v_1x
         w(ixmin1:ixmax1,mom(2)) = v_1y
         w(ixmin1:ixmax1,p_)     = p_1
      elsewhere ((x(ixmin1:ixmax1,1) > x_th_ini) .and. (x(ixmin1:ixmax1,1) < x_th_fin))
         w(ixmin1:ixmax1,rho_)   = rho_th
         w(ixmin1:ixmax1,mom(1)) = v_0x
         w(ixmin1:ixmax1,mom(2)) = v_0y
         w(ixmin1:ixmax1,p_)     = p_eq
      elsewhere
         w(ixmin1:ixmax1,rho_)   = rho_eq
         w(ixmin1:ixmax1,mom(1)) = v_0x
         w(ixmin1:ixmax1,mom(2)) = v_0y
         w(ixmin1:ixmax1,p_)     = p_eq
      end where
      
      if (magnetic .eqv. .true.) then
         where (x(ixmin1:ixmax1,1) < x_sh)
            w(ixmin1:ixmax1,mag(1))   = b_1x
            w(ixmin1:ixmax1,mag(2))   = b_1y
            w(ixmin1:ixmax1,mag(3))   = b_1z
         elsewhere ((x(ixmin1:ixmax1,1) > x_th_ini) .and. (x(ixmin1:ixmax1,1) < x_th_fin))
            w(ixmin1:ixmax1,mag(1))   = b_0x
            w(ixmin1:ixmax1,mag(2))   = b_0y
            w(ixmin1:ixmax1,mag(3))   = b_0z
         elsewhere
            w(ixmin1:ixmax1,mag(1))   = b_0x
            w(ixmin1:ixmax1,mag(2))   = b_0y
            w(ixmin1:ixmax1,mag(3))   = b_0z
         end where
      end if

    if (magnetic .eqv. .false.) call hd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    if (magnetic .eqv. .true.) call mhd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)

  end subroutine rm1d_init_one_grid

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
         x(1,i) = -2 + i
      end do
end subroutine generate_particles

end module mod_usr
