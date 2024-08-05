subroutine find_orographic_precipitation (geometry,params)

  ! computes orographic precipitation

  use definitions
  use lfpm
  use interpolation
  use file_io
  use array_utils

  implicit none

  type(geom) :: geometry
  type(parm) :: params
  type(Grid) :: g

  type(netw) :: network
  type(stck) :: stack

  double precision :: x_min, x_max, y_min, y_max, dx, dy, xi, yi
  integer :: n, nx, ny, i, j, ierr
  double precision, allocatable :: x_merged(:), y_merged(:), z_merged(:), pr(:), x_trailed(:), y_trailed(:), z_trailed(:) 
  double precision, allocatable :: zpr(:,:)

  call time_in ('find_orographic_precipitation')

  ! interpolation
  nx = params%nx
  ny = params%ny

  allocate(pr(size(geometry%z)))

  ! ! non-uniform data points (x_merged, y_merged, z_merged)
  allocate(x_merged(size(geometry%x)+size(geometry%xdiv)))
  allocate(y_merged(size(geometry%y)+size(geometry%ydiv)))
  allocate(z_merged(size(geometry%z)+size(geometry%zdiv)))

  allocate(x_trailed(size(geometry%x)+size(geometry%xdiv)))
  allocate(y_trailed(size(geometry%y)+size(geometry%ydiv)))
  allocate(z_trailed(size(geometry%z)+size(geometry%zdiv)))

  x_merged = (/geometry%x, geometry%xdiv/)
  y_merged = (/geometry%y, geometry%ydiv/)
  z_merged = (/geometry%z, geometry%zdiv/)

  call remove_trailing_zeros(x_merged, x_trailed)
  call remove_trailing_zeros(y_merged, y_trailed)
  call remove_trailing_zeros(z_merged, z_trailed)

  ! n = size(geometry%x)
  n = size(x_trailed)

  allocate(zpr(nx, ny))  ! Regular grid

  ! Define the bounds of the grid
  x_min = minval(geometry%x)
  x_max = maxval(geometry%x)
  y_min = minval(geometry%y)
  y_max = maxval(geometry%y)

  ! Define grid spacing
  dx = (x_max - x_min) / (nx - 1)
  dy = (y_max - y_min) / (ny - 1)

  ! Initialize the grid to zero
  zpr = 0.0

  ! Perform bilinear interpolation
  call time_in ('bilinear_interpolation')
  !$omp parallel do private(i, j, xi, yi)
  do i = 1, nx
    do j = 1, ny
      xi = x_min + (i-1) * dx
      yi = y_min + (j-1) * dy
      call bilinear_interpolation(x_trailed, y_trailed, z_trailed, n, xi, yi, zpr(i, j))
      ! call bilinear_interpolation(geometry%x, geometry%y, geometry%z, n, xi, yi, zpr(i, j))
    end do
  end do
  !$omp end parallel do
  call time_out ('bilinear_interpolation')

  g%m = nx 
  g%n = ny
  g%per = params%per
  g%lc = params%lc
  g%lf = params%lf
  g%ll = params%ll
  g%ld = params%ld
  g%refheight = params%refheight
  g%evap = params%evap
  g%qin = params%qin

  allocate(g%diag(ny))
  allocate(g%upper(ny-1))
  allocate(g%lower(ny-1))
  allocate(g%rhs(ny))
  allocate(g%psi(ny))
  allocate(g%u(nx, ny))
  if (g%per == 1) then
  allocate(g%right(ny-2))
  allocate(g%bottom(ny-2))
  end if

  do i = 1, nx
    do j = 1, ny
    g%u(i, j)%h = zpr(i, j) ! sign the topography to it
    end do
  end do

  call time_in ('computePrecipitation')
  call g%computePrecipitation()
  call time_out ('computePrecipitation')

  ! Perform reverse interpolation
  call time_in ('reverse_bilinear_interpolation')
  call reverse_interpolation(geometry%x, geometry%y, pr, size(geometry%x), g%u%ptot, nx, ny, x_min, dx, y_min, dy)
  call time_out ('reverse_bilinear_interpolation')

  geometry%precipitation=pr*3.156d7+1.d-2

  deallocate(g%diag)
  deallocate(g%upper)
  deallocate(g%lower)
  deallocate(g%rhs)
  deallocate(g%psi)
  deallocate(g%u)
  if (params%per == 1) then
  deallocate(g%right)
  deallocate(g%bottom)
  end if

  deallocate(zpr)
  deallocate(pr)

  deallocate(x_merged)
  deallocate(y_merged)
  deallocate(z_merged)
  
  deallocate(x_trailed)
  deallocate(y_trailed)
  deallocate(z_trailed)

  call time_out ('find_orographic_precipitation')
  return

end subroutine find_orographic_precipitation
