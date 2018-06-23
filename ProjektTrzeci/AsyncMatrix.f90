   !---------------------------------------------------------------------------  
   !> @author 
   !> Maciej Klimowski
   !
   ! DESCRIPTION: 
   !> Async routine multiply two matrixes 
   !
   ! REVISION HISTORY:
   ! 23 06 2018 - Initial Version
   !
   !> @param[in] n - first dimension
   !> @param[in] m - second dimension
   !> @param[in] first - first matrix
   !> @param[in] second - second matrix
   !> @param[out] resultMatrix - the result of multiplication
   !--------------------------------------------------------------------------- 
subroutine multiplication(first, second, n, m, resultMatrix)
    implicit none
    real (kind = 8), dimension(n, m), intent(in) :: first
    real (kind = 8), dimension(m, n), intent(in) :: second
    real (kind = 8), dimension(m, m), intent(out) :: resultMatrix
    real (kind = 8), codimension[:], dimension(:,:), allocatable :: buffor
    integer, intent(in) :: n, m
    integer, codimension[:], allocatable :: minrow, maxrow
    integer :: i, j, k, rows, im
    
    rows = CEILING(real(n) / NUM_IMAGES())
    minrow = MIN(n, (THIS_IMAGE() - 1) * rows) + 1
    maxrow = MIN(n, THIS_IMAGE() * rows)
    
    allocate(buffor(rows, m)[*])
    
    do i = minrow, maxrow
        do j = 1, n
            do k = 1, m
                buffor(i - minrow + 1, j) = buffor(i - minrow + 1, j) + first(i, k) * second(k, j)
            end do
        end do
    end do
    
    sync all
    if(THIS_IMAGE() .EQ. 1) then
        do im = 1, NUM_IMAGES()
            resultMatrix(minrow[im]:maxrow[im], :) = buffor(1:(maxrow[im] - minrow[im] + 1) :)[im]
        end do
    end if
    
    deallocate(buffor)            
end subroutine multiplication
    
   !---------------------------------------------------------------------------  
   !> @author 
   !> Maciej Klimowski
   !
   ! DESCRIPTION: 
   !> Async routine doing gauss elimination
   !
   ! REVISION HISTORY:
   ! 23 06 2018 - Initial Version
   !
   !> @param[in] n - size of matrix
   !> @param[inout] first - main matrix
   !> @param[inout] second - right vector 
   !---------------------------------------------------------------------------   
subroutine gauss(first, second, n))
    implicit none
    real (kind = 8), dimension(n, n), intent(inout) :: first
    real (kind = 8), dimension(n), intent(inout) :: second
    real (kind = 8), codimension[:], allocatable :: asyncFirst(:, :)
    real (kind = 8), codimension[:], allocatable :: asyncSecond(:)
    integer, intent(in) :: n
    integer :: i, j
    real (kind = 8) :: val
    
    allocate(asyncFirst(0:n, 0:n)[*])
    allocate(asyncSecond(0:n)[*])
    
    if(THIS_IMAGE() .EQ. 1) then
        asyncFirst(:, :)[1] = A(:,:)
        asyncSecond(:)[1] = X(:)
    end if
        
    
    do i = 0, n
      if(THIS_IMAGE() .eq. 1) then
        asyncSecond(i)[1] = asyncSecond(i)[1] / asyncFirst(i, i)[1]
        asyncFirst(:, i)[1] = asyncFirst(:, i)[1] / asyncFirst(i, i)[1]
      end if
      sync all
      do j = THIS_IMAGE() - 1, n, NUM_IMAGES()
        if ((i .NE. j) .AND. (ABS(coA(i, i)[1] - 0) > 1d-6)) then
          val = asyncFirst(i, j)[1] / asyncFirst(i, i)[1]
          asyncFirst(:,j)[1] = asyncFirst(:,j)[1] - val * asyncFirst(:, i)[1]
          asyncSecond(j)[1] = asyncSecond(j)[1] - val * asyncSecond(i)[1]
        end if
      end do
    end do
    
    if (THIS_IMAGE() .EQ. 1) then
        first(:,:) = asyncFirst(:, :)[1]
        second(:) = asyncSecond(:)[1]
    end if
    
    deallocate(asyncFirst)
    deallocate(asyncSecond)
end subroutine gauss
    