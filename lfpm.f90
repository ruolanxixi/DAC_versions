module lfpm
  implicit none
  private
  public :: Vect, Matrix, Node, Grid

  type :: Vect
    double precision :: u1 = 0.0, u2 = 0.0
  contains
    procedure :: subtract => vect_subtract
  end type Vect

  type :: Matrix
    double precision :: a11 = 0.0, a12 = 0.0, a21 = 0.0, a22 = 0.0
  contains
    procedure :: subtract => matrix_subtract
    procedure :: multiply_vect => matrix_multiply_vect
    procedure :: multiply_matrix => matrix_multiply_matrix
    procedure :: inv => matrix_inv
  end type Matrix

  type :: Node
    double precision :: h = 0.0
    double precision :: qv = 0.0, qc = 0.0, ptot = 0.0, p = 0.0
  end type Node

  type :: Grid
    integer :: m, n, per
    double precision :: lc, lf, ll, ld, refheight, evap, qin
    type(Matrix), allocatable :: diag(:), upper(:), lower(:), right(:), bottom(:)
    type(Vect), allocatable :: rhs(:)
    double precision, allocatable :: psi(:)
    type(Node), allocatable :: u(:,:)
  contains
    procedure :: computePrecipitation => grid_compute_precipitation
  end type Grid

contains

  function vect_subtract(this, u) result(res)
    class(Vect), intent(in) :: this, u
    type(Vect) :: res
    res%u1 = this%u1 - u%u1
    res%u2 = this%u2 - u%u2
  end function vect_subtract

  function matrix_subtract(this, a) result(res)
    class(Matrix), intent(in) :: this, a
    type(Matrix) :: res
    res%a11 = this%a11 - a%a11
    res%a12 = this%a12 - a%a12
    res%a21 = this%a21 - a%a21
    res%a22 = this%a22 - a%a22
  end function matrix_subtract

  function matrix_multiply_vect(this, b) result(res)
    class(Matrix), intent(in) :: this
    class(Vect), intent(in) :: b
    type(Vect) :: res
    res%u1 = this%a11 * b%u1 + this%a12 * b%u2
    res%u2 = this%a21 * b%u1 + this%a22 * b%u2
  end function matrix_multiply_vect

  function matrix_multiply_matrix(this, b) result(res)
    class(Matrix), intent(in) :: this, b
    type(Matrix) :: res
    res%a11 = this%a11 * b%a11 + this%a12 * b%a21
    res%a12 = this%a11 * b%a12 + this%a12 * b%a22
    res%a21 = this%a21 * b%a11 + this%a22 * b%a21
    res%a22 = this%a21 * b%a12 + this%a22 * b%a22
  end function matrix_multiply_matrix

  function matrix_inv(this) result(res)
    class(Matrix), intent(in) :: this
    type(Matrix) :: res
    double precision :: det
    det = this%a11 * this%a22 - this%a12 * this%a21
    if (det /= 0.0) then
      res%a11 = this%a22 / det
      res%a12 = -this%a12 / det
      res%a21 = -this%a21 / det
      res%a22 = this%a11 / det
    end if
  end function matrix_inv

  subroutine grid_compute_precipitation(this)
    class(Grid), intent(inout) :: this
    double precision :: beta, f
    integer :: i, j
    type(Matrix) :: temp_matrix_f, inv_diag, temp_matrix_g

    beta = (1.0 - this%lc / this%ll) * (this%ll / this%lf - 1.0)

    !$OMP PARALLEL DO PRIVATE(j, f)
    do j = 1, this%n
      f = this%lf / (this%ll - this%lf) * exp(this%u(1, j)%h / this%refheight)
      this%rhs(j)%u1 = this%qin / (1.0 + f)
      this%rhs(j)%u2 = this%qin * f / (1.0 + f)
      this%u(1, j)%qv = this%rhs(j)%u1
      this%u(1, j)%qc = this%rhs(j)%u2
    end do
    !$OMP END PARALLEL DO

    ! Loop over the rows of the grid. Each rows uses the fluxes of the previous row.
    do i = 2, this%m
      !$OMP PARALLEL DO PRIVATE(j, f) SHARED(this, beta)
      do j = 1, this%n
        f = exp(-this%u(i, j)%h / this%refheight)
        this%psi(j) = 1.0 - f * this%evap
        f = f * beta / this%lc
        if (this%per == 1 .or. (j > 1 .and. j < this%n)) then
          this%diag(j)%a11 = 2.0 * this%ld + 1.0 + 1.0 / this%lc
          this%diag(j)%a12 = -f - (1.0 - this%psi(j)) / this%lf
          this%diag(j)%a21 = -1.0 / this%lc
          this%diag(j)%a22 = 2.0 * this%ld + 1.0 + f + 1.0 / this%lf
        else
          this%diag(j)%a11 = this%ld + 1.0 + 1.0 / this%lc
          this%diag(j)%a12 = -f - (1.0 - this%psi(j)) / this%lf
          this%diag(j)%a21 = -1.0 / this%lc
          this%diag(j)%a22 = this%ld + 1.0 + f + 1.0 / this%lf
        end if
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(j) SHARED(this)
      do j = 1, this%n - 1
        this%lower(j)%a11 = -this%ld
        this%lower(j)%a12 = 0.0
        this%lower(j)%a21 = 0.0
        this%lower(j)%a22 = -this%ld
        this%upper(j)%a11 = -this%ld
        this%upper(j)%a12 = 0.0
        this%upper(j)%a21 = 0.0
        this%upper(j)%a22 = -this%ld
      end do
      !$OMP END PARALLEL DO
      if (this%per == 1) then
        !$OMP PARALLEL
        this%right(1)%a11 = -this%ld
        this%right(1)%a12 = 0.0
        this%right(1)%a21 = 0.0
        this%right(1)%a22 = -this%ld
        this%bottom(1)%a11 = -this%ld
        this%bottom(1)%a12 = 0.0
        this%bottom(1)%a21 = 0.0
        this%bottom(1)%a22 = -this%ld
        !$OMP DO
        do j = 2, this%n - 2
          this%right(j)%a11 = 0.0
          this%right(j)%a12 = 0.0
          this%right(j)%a21 = 0.0
          this%right(j)%a22 = 0.0
          this%bottom(j)%a11 = 0.0
          this%bottom(j)%a12 = 0.0
          this%bottom(j)%a21 = 0.0
          this%bottom(j)%a22 = 0.0
        end do
        !$OMP END DO
        !$OMP END PARALLEL
      end if

      do j = 1, this%n - 1
        inv_diag = matrix_inv(this%diag(j))
        temp_matrix_f = matrix_multiply_matrix(inv_diag, this%lower(j))
        this%diag(j + 1) = matrix_subtract(this%diag(j + 1), matrix_multiply_matrix(temp_matrix_f, this%upper(j)))
        this%rhs(j + 1) = vect_subtract(this%rhs(j + 1), matrix_multiply_vect(temp_matrix_f, this%rhs(j)))

        if (this%per == 1 .and. j < this%n - 1) then
          temp_matrix_g = matrix_multiply_matrix(inv_diag, this%bottom(j))
          if (j < this%n - 2) then
            this%right(j + 1) = matrix_subtract(this%right(j + 1), matrix_multiply_matrix(temp_matrix_f, this%right(j)))
            this%bottom(j + 1) = matrix_subtract(this%bottom(j + 1), matrix_multiply_matrix(temp_matrix_g, this%upper(j)))
          else
            this%upper(j + 1) = matrix_subtract(this%upper(j + 1), matrix_multiply_matrix(temp_matrix_f, this%right(j)))
            this%lower(j + 1) = matrix_subtract(this%lower(j + 1), matrix_multiply_matrix(temp_matrix_g, this%upper(j)))
          end if
          this%diag(this%n) = matrix_subtract(this%diag(this%n), matrix_multiply_matrix(temp_matrix_g, this%right(j)))
          this%rhs(this%n) = vect_subtract(this%rhs(this%n), matrix_multiply_vect(temp_matrix_g, this%rhs(j)))
        end if
      end do

      this%rhs(this%n) = matrix_multiply_vect(matrix_inv(this%diag(this%n)), this%rhs(this%n))
      this%rhs(this%n - 1) = vect_subtract(this%rhs(this%n - 1), matrix_multiply_vect(this%upper(this%n - 1), this%rhs(this%n)))
      this%rhs(this%n - 1) = matrix_multiply_vect(matrix_inv(this%diag(this%n - 1)), this%rhs(this%n - 1))

      do j = this%n - 2, 1, -1
        this%rhs(j) = vect_subtract(this%rhs(j), matrix_multiply_vect(this%upper(j), this%rhs(j + 1)))
        if (this%per == 1) this%rhs(j) = vect_subtract(this%rhs(j), matrix_multiply_vect(this%right(j), this%rhs(this%n)))
        this%rhs(j) = matrix_multiply_vect(matrix_inv(this%diag(j)), this%rhs(j))
      end do

      do j = 1, this%n
        this%u(i, j)%qv = this%rhs(j)%u1
        this%u(i, j)%qc = this%rhs(j)%u2
        this%u(i, j)%ptot = this%rhs(j)%u2 / this%lf
        this%u(i, j)%p = this%psi(j) * this%rhs(j)%u2 / this%lf
      end do
    end do
  end subroutine grid_compute_precipitation
end module lfpm


module interpolation
  implicit none
  contains
  subroutine bilinear_interpolation(x, y, z, n, xi, yi, zi)
    implicit none
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: x, y, z
    double precision, intent(in) :: xi, yi
    double precision, intent(out) :: zi
    double precision :: w_sum, w, dx, dy
    integer :: k

    w_sum = 0.0
    zi = 0.0
    !$omp parallel do private(k, dx, dy, w) reduction(+:zi, w_sum)
    do k = 1, n
        dx = xi - x(k)
        dy = yi - y(k)
        w = 1.0 / (dx * dx + dy * dy + 1.0e-10)  ! Avoid division by zero
        zi = zi + w * z(k)
        w_sum = w_sum + w
    end do
    !$omp end parallel do
    if (w_sum /= 0.0) zi = zi / w_sum
  end subroutine bilinear_interpolation


  subroutine reverse_interpolation(x, y, z, n, grid, nx, ny, x_min, dx, y_min, dy)
    implicit none
    integer, intent(in) :: n, nx, ny
    double precision, dimension(n), intent(in) :: x, y
    double precision, dimension(nx, ny), intent(in) :: grid
    double precision, intent(out) :: z(n)
    double precision, intent(in) :: x_min, y_min, dx, dy
    integer :: i, j, ix, iy
    double precision :: xi, yi, x1, x2, y1, y2, f11, f12, f21, f22, w11, w12, w21, w22
    !$omp parallel do private(i, xi, yi, ix, iy, x1, x2, y1, y2, f11, f12, f21, f22, w11, w12, w21, w22)
    do i = 1, n
      xi = x(i)
      yi = y(i)
      
      ! Determine the grid cell (ix, iy) that contains (xi, yi)
      ix = int((xi - x_min) / dx) + 1
      iy = int((yi - y_min) / dy) + 1

      ! Ensure indices are within bounds
      if (ix < 1 .or. ix >= nx .or. iy < 1 .or. iy >= ny) then
        z(i) = 0.0  ! Default value to avoid NaN
      else
        x1 = x_min + (ix-1) * dx
        x2 = x1 + dx
        y1 = y_min + (iy-1) * dy
        y2 = y1 + dy
          
        f11 = grid(ix, iy)
        f12 = grid(ix, iy+1)
        f21 = grid(ix+1, iy)
        f22 = grid(ix+1, iy+1)
        
        w11 = (x2 - xi) * (y2 - yi)
        w12 = (x2 - xi) * (yi - y1)
        w21 = (xi - x1) * (y2 - yi)
        w22 = (xi - x1) * (yi - y1)
        
        z(i) = (w11 * f11 + w12 * f12 + w21 * f21 + w22 * f22) / (dx * dy)
      end if
    end do
    !$omp end parallel do
  end subroutine reverse_interpolation
end module interpolation

module nan_handling
  implicit none
  contains

  ! Function to check if a value is NaN
  logical function is_nan(value)
    implicit none
    double precision, intent(in) :: value
    is_nan = (value /= value)  ! A value is NaN if it is not equal to itself
  end function is_nan

  ! Subroutine to merge arrays and move NaN values to the end
  subroutine merge_and_move_nans(x, xdiv, x_merged)
    implicit none
    double precision, intent(in) :: x(:), xdiv(:)
    double precision, intent(out) :: x_merged(:)
    double precision, allocatable :: temp(:)
    integer :: n, m, i, non_nan_count, nan_count

    ! Get sizes of input arrays
    n = size(x)
    m = size(xdiv)

    ! Allocate the merged array
    ! allocate(x_merged(n + m))

    ! Merge the arrays
    x_merged = (/ x, xdiv /)

    ! Allocate a temporary array to store the non-NaN values
    allocate(temp(size(x_merged)))

    ! Move non-NaN values to the beginning of the array
    non_nan_count = 0
    do i = 1, size(x_merged)
      if (.not. is_nan(x_merged(i))) then
        non_nan_count = non_nan_count + 1
        temp(non_nan_count) = x_merged(i)
      end if
    end do

    ! Move NaN values to the end of the array
    nan_count = non_nan_count
    do i = 1, size(x_merged)
      if (is_nan(x_merged(i))) then
        nan_count = nan_count + 1
        temp(nan_count) = x_merged(i)
      end if
    end do

    ! Assign the rearranged values back to x_merged
    x_merged = temp

    ! Deallocate the temporary array
    deallocate(temp)
  end subroutine merge_and_move_nans

end module nan_handling

module file_io
  implicit none
contains
  subroutine write_data_to_files(x, y, z, n, grid, nx, ny, x_min, dx, y_min, dy)
    implicit none
    integer, intent(in) :: n, nx, ny
    double precision, dimension(n), intent(in) :: x, y, z
    double precision, dimension(nx, ny), intent(in) :: grid
    double precision, intent(in) :: x_min, y_min, dx, dy
    integer :: i, j

    ! Write original data to a file
    open(unit=10, file='original_data.dat', status='replace')
    do i = 1, n
      write(10, *) x(i), y(i), z(i)
    end do
    close(10)

    ! Write gridded data to a file
    open(unit=11, file='gridded_data.dat', status='replace')
    do i = 1, nx
      do j = 1, ny
        write(11, *) x_min + (i-1) * dx, y_min + (j-1) * dy, grid(i, j)
      end do
    end do
    close(11)
  end subroutine write_data_to_files
end module file_io

module array_utils
    implicit none
    private
    public :: remove_trailing_zeros

  contains

    subroutine remove_trailing_zeros(array, new_array)
        implicit none
        double precision, intent(in) :: array(:)
        double precision, allocatable, intent(out) :: new_array(:)
        integer :: i, last_nonzero

        ! Find the position of the last non-zero element
        last_nonzero = 0
        do i = size(array), 1, -1
            if (array(i) /= 0) then
                last_nonzero = i
                exit
            end if
        end do

        ! Allocate new array with size up to the last non-zero element
        allocate(new_array(last_nonzero))
        new_array = array(1:last_nonzero)
    end subroutine remove_trailing_zeros

end module array_utils
