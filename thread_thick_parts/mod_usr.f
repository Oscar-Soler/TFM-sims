module mod_usr
  use mod_hd

  implicit none

contains

  subroutine usr_init()

      usr_init_one_grid => rm1d_init_one_grid

      usr_create_particles => generate_particles

    call set_coordinate_system("Cartesian")
    call hd_activate()

  end subroutine usr_init

  ! Initialize one grid
  subroutine rm1d_init_one_grid(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)
    integer, intent(in) :: ixGmin1,ixGmax1, ixmin1,ixmax1
    double precision, intent(in) :: x(ixGmin1:ixGmax1,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,1:nw)
    double precision :: M, rho_eq, v_0, v_eq, p_eq, rho_1, v_1, v_sh, p_1, M1,&
      cs_0
    double precision :: q_th, rho_th, x_th_ini, x_th_fin, x_sh

    M = 5.405405405405405d0
    rho_eq = 1d0
    p_eq = 1d0

    q_th = 100d0
      x_th_ini = 2
      x_th_fin = 32
      x_sh = 0
      
    !------------------------------------- Above parameters, below calculations
      
    p_1 = p_eq * ((2*hd_gamma)/(hd_gamma+1)*M**2 - (hd_gamma-1)/(hd_gamma+1))
    rho_1 = rho_eq / (2/((hd_gamma+1)*M**2) + (hd_gamma-1)/(hd_gamma+1))
    M1 = sqrt( (2+(hd_gamma-1)*M**2) / (2*hd_gamma*M**2 - (hd_gamma-1)))
    cs_0 = sqrt(hd_gamma * p_eq / rho_eq)
    !v_eq = sqrt(hd_gamma * p_eq / rho_eq)
    !v_0 = -M*cs_0
    v_0 = 0
    v_sh = v_0 + M*cs_0
    v_1 = rho_eq/rho_1 * (v_0-v_sh) + v_sh
      !print *, rho_1/rho_eq
      !print *, p_1/p_eq
    rho_th = rho_eq * q_th
      
    where (x(ixmin1:ixmax1,1) < x_sh)
        w(ixmin1:ixmax1,rho_)   = rho_1
        w(ixmin1:ixmax1,mom(1)) = v_1
        w(ixmin1:ixmax1,p_)     = p_1
    elsewhere ((x(ixmin1:ixmax1,1) > x_th_ini) .and. (x(ixmin1:ixmax1,1) < x_th_fin))
        w(ixmin1:ixmax1,rho_)   = rho_th
        w(ixmin1:ixmax1,mom(1)) = v_0
        w(ixmin1:ixmax1,p_)     = p_eq
    elsewhere
        w(ixmin1:ixmax1,rho_)   = rho_eq
        w(ixmin1:ixmax1,mom(1)) = v_0
        w(ixmin1:ixmax1,p_)     = p_eq
     
    end where
    
    call hd_to_conserved(ixGmin1,ixGmax1,ixmin1,ixmax1,w,x)

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
         x(1,i) = -10 + 5*i
      end do
end subroutine generate_particles

end module mod_usr
