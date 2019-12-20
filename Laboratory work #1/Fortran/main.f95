module GaussMethod
    integer systemSize
    real, dimension (:,:), allocatable :: coeff
    real, dimension (:), allocatable :: freeTerm
    real, dimension (:), allocatable :: solutions
    character(12) :: filename
    integer isSolutionFounded

    contains
        subroutine init()
            systemSize = 0
            filename = "matrix.txt"
            isSolutionFounded = 0
        end subroutine init

        subroutine destruct()
            deallocate (coeff)
            deallocate (freeTerm)
            deallocate (solutions)
        end subroutine destruct

        subroutine initializationFromUserInput()
            print *, 'Введите количество уравнений: '
            read *, systemSize

            allocate ( coeff(systemSize,systemSize+1) )
            allocate ( freeTerm(systemSize) )
            allocate ( solutions(systemSize) )

            print *, 'Введите коэффициенты: '
            do i = 1, systemSize
                print *, i, ' уравнение: '
                do j = 1, systemSize
                    print *, 'a[', i, '][', j, ']='
                    read *, coeff(i,j)
                end do

                print *, 'Введите свободный член: '
                print *, 'y[', i, ']= '
                read *, freeTerm(i)
            end do

        end subroutine initializationFromUserInput

        subroutine initializationFromFile()
            real tmp

            open(1,file=filename)

            do
                read(1,*, end=1) tmp
                systemSize = systemSize + 1
            end do

            1 close (1)

            allocate ( coeff(systemSize,systemSize+1) )
            allocate ( freeTerm(systemSize) )
            allocate ( solutions(systemSize) )

            open(2,file=filename)

            do i = 1,systemSize
                read(2,*)(coeff(i,j),j=1,systemSize+1)
            end do

            do i = 1,systemSize
                freeTerm(i) = coeff(i,systemSize+1)
                coeff(i,systemSize+1) = 0
            end do

            close (2)
        end subroutine

        subroutine showSystem()
            write (*,*) 'Введённая система: '
            do i = 1, systemSize
                do j = 1, systemSize
                    if (coeff(i,j) < 0) then
                        write (*,*) '(', coeff(i,j), ')*x', j
                    else
                        write (*,*) coeff(i,j), '*x', j
                    end if
                    if (j < systemSize) then
                        write (*,*) ' + '
                    end if
                end do
                write (*,*) ' = ', freeTerm(i)
            end do

        end subroutine showSystem

        integer function rankOfMatrix(typeOfMatrix)
            integer typeOfMatrix
            integer :: dimensionOfMatrix = 0
            integer :: returnValue = 0
            integer :: k, c = 0
            real :: tmp = 1

            if (typeOfMatrix == 0) then
                dimensionOfMatrix = systemSize
            else if (typeOfMatrix == 1) then
                dimensionOfMatrix = systemSize + 1
            end if

            returnValue = 0
            do i = 1,systemSize
                if (coeff(i,i+typeOfMatrix) .NE. 0) then
                    returnValue = returnValue + 1
                else
                    exit
                end if
            end do

            rankOfMatrix = returnValue
        end function rankOfMatrix

        subroutine applyGaussMethod()
            real x(systemSize)
            integer rankOfUsualMatrix, rankOfExtendedMatrix
            real c

            do i = 1, systemSize
                coeff(i,systemSize+1) = freeTerm(i)
            end do

            do i = systemSize,2,-1
                if (coeff(i-1,1) < coeff(i,1)) then
                    do j = 1,systemSize+1
                        c = coeff(i,j)
                        coeff(i,j) = coeff(i-1,j)
                        coeff(i-1,j) = c
                    end do
                end if
            end do

            do k = 1,systemSize-1
                do i = k, systemSize-1
                    c = (coeff(i+1,k) / coeff(k,k))

                    do j = 1,systemSize+1
                        coeff(i+1,j) = coeff(i+1,j) - c*coeff(k,j)
                    end do
                end do
            end do

            rankOfUsualMatrix = rankOfMatrix(0)
            rankOfExtendedMatrix = rankOfMatrix(1)

            if (rankOfUsualMatrix < rankOfExtendedMatrix) then
                print *, 'Нет решения!'
                stop
            else if ( (rankOfUsualMatrix == rankOfExtendedMatrix).and.(rankOfUsualMatrix < systemSize) ) then
                print *, 'Система имеет бесконечное множество решений'
                stop
            else if ( (rankOfUsualMatrix == rankOfExtendedMatrix).and.(rankOfUsualMatrix == systemSize) ) then
                print *, 'Система имеет единственное решение!'
            end if

            do i = systemSize, 1, -1
                !real c
                c = 0
                do j = i,systemSize
                    c = c + coeff(i,j) * x(j)
                end do

                x(i) = (coeff(i,systemSize+1) - c) / coeff(i,i)
            end do

            isSolutionFounded = 1

            do i = 1,systemSize
                solutions(i) = x(i)
            end do

        end subroutine applyGaussMethod

        subroutine showSolutions()
            if (isSolutionFounded == 1) then
                print *, 'Решение введённой системы:'
                do i = 1,systemSize
                    print *, 'x[', i, ']=', solutions(i)
                end do
            end if
        end subroutine showSolutions

end module


program main
    use gaussMethod
    call init()

    !call initializationFromUserInput()
    call initializationFromFile()
    !call showSystem()
    call applyGaussMethod()
    call showSolutions()

    call destruct()
end
