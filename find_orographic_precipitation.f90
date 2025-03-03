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
  integer :: n, nx, ny, i, j, ierr, num_points, num_neighbours
  double precision, allocatable :: x_merged(:), y_merged(:), z_merged(:), pr(:), x_trailed(:), y_trailed(:), z_trailed(:)
  double precision, allocatable :: points(:, :), geo_points(:, :), x_axis(:), y_axis(:), zpr(:, :)

  call time_in ('find_orographic_precipitation')

  ! interpolation
  ! nx = params%nx
  ! ny = params%ny

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
  num_points = size(x_trailed)

  allocate(pr(size(geometry%x)))

  allocate(points(num_points, 2))

  points(:, 1) = x_trailed
  points(:, 2) = y_trailed

  ! Define the bounds of the grid
  x_min = minval(x_trailed)
  x_max = maxval(x_trailed)
  y_min = minval(y_trailed)
  y_max = maxval(y_trailed)

  dx = 1000.0
  dy = 1000.0

  ! Compute nx and ny, rounding down to the nearest integer
  nx = int(ceiling((x_max - x_min) / dx)) + 1
  ny = int(ceiling((y_max - y_min) / dy)) + 1

  ! Allocate arrays for x and y grids
  allocate(x_axis(nx))
  allocate(y_axis(ny))

  allocate(zpr(ny, nx))  ! Regular grid

  ! Populate x and y grid arrays
  do i = 1, nx
     x_axis(i) = x_min + (i - 1) * dx
  end do
  do j = 1, ny
     y_axis(j) = y_min + (j - 1) * dy
  end do

  ! Initialize the grid to zero
  zpr = 0.0

  num_neighbours = 6

  ! Perform inverse distance weighted interpolation (unstructured -> regular grid)
  call time_in ('inverse_distance_weighted_interpolation')
  call idw_esrg_nearest(points, num_points, z_trailed, x_axis, nx, y_axis, ny, dx, num_neighbours, zpr)
  ! idw_esrg_nearest(points, num_points, data_in, & x_axis, len_x, y_axis, len_y, grid_spac, num_neighbours, data_ip)
  call time_out ('inverse_distance_weighted_interpolation')

  call write_data_to_files(x_trailed, y_trailed, z_trailed, size(x_trailed), zpr, nx, ny, x_min, dx, y_min, dy)

  ! !$omp parallel do private(i, j, xi, yi)
  ! do i = 1, nx
  !   do j = 1, ny
  !     xi = x_min + (i-1) * dx
  !     yi = y_min + (j-1) * dy
  !     call bilinear_interpolation(x_trailed, y_trailed, z_trailed, n, xi, yi, zpr(i, j))
  !     ! call bilinear_interpolation(geometry%x, geometry%y, geometry%z, n, xi, yi, zpr(i, j))
  !   end do
  ! end do
  ! !$omp end parallel do

  g%m = ny ! number of rows
  g%n = nx ! number of columns
  g%per = params%per
  g%lc = params%lc
  g%lf = params%lf
  g%ll = params%ll
  g%ld = params%ld
  g%refheight = params%refheight
  g%evap = params%evap
  g%qin = params%qin

  allocate(g%diag(g%n))
  allocate(g%upper(g%n-1))
  allocate(g%lower(g%n-1))
  allocate(g%rhs(g%n))
  allocate(g%psi(g%n))
  allocate(g%u(g%m, g%n))
  if (g%per == 1) then
  allocate(g%right(g%n-2))
  allocate(g%bottom(g%n-2))
  end if

  do j = 1, ny
    do i = 1, nx
    g%u(j, i)%h = zpr(j, i) ! assign the topography to it
    end do
  end do

  call time_in ('compute_precipitation')
  call g%computePrecipitation()
  call time_out ('compute_precipitation')

  allocate(geo_points(size(geometry%x), 2))

  geo_points(:, 1) = geometry%x
  geo_points(:, 2) = geometry%y

  ! Perform reverse interpolation
  call time_in ('bilinear_interpolation')
  call bilinear(x_axis, nx, y_axis, ny, g%u%ptot, geo_points, size(geometry%x), pr)
  ! bilinear(x_axis, len_x, y_axis, len_y, data_in, points, num_points, data_ip)
  call time_out ('bilinear_interpolation')

  call write_data_to_files(geometry%x, geometry%y, pr, size(geometry%x), g%u%ptot, nx, ny, x_min, dx, y_min, dy)

  geometry%precipitation=pr+1.d-2  !*3.156d7+1.d-2

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
