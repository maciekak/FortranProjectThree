subroutine multiplication(first, second, n, m, resultMatrix)
    implicit none
    real (kind = 8), dimension(n, m), intent(in) :: first
    real (kind = 8), dimension(m, n), intent(in) :: second
    real (kind = 8), dimension(m, m), intent(out) :: resultMatrix
    integer, intent(in) :: n, m
    real (kind = 8) :: value
    integer :: i, j, k
    
    !f2py intent (out) :: resultMatrix
    do i = 1, m
        do j = 1, n
            value = 0.d
            do k = 1, n
                value = value + first(k, i) * second(j, k)
            end do
            resultMatrix(i, j) = value
        end do
    end do
end subroutine multiplication
    
    
subroutine gauss(first, second, n))
    implicit none
    real (kind = 8), dimension(n, n), intent(inout) :: first
    real (kind = 8), dimension(n), intent(inout) :: second
    integer, intent(in) :: n
    integer :: i, j
    real (kind = 8) :: val
    !f2py intent(in, out) first, second
    
    do i = 1, n
        do j = 1, n
            if (i .NE. j) then
                val = first(i, j) / first(i, i)
                first(:, j) = first(:, j) - val * first(:, i)
                second(j) = second(j) - val * second(i)
                second(i) = second(i) / first(i, i)
                first(:, i) = first(:, i) / first(i, i)
            end if
        end do
    end do
end subroutine gauss
    