module GaussMethod
    integer systemSize
    real, dimension (:,:), allocatable :: coeff
    real, dimension (:), allocatable :: freeTerm
    !real, dimension (:), allocatable :: solutions
    character(12) :: filename
    real eps
    !integer isSolutionFounded

    contains
        subroutine init()
            systemSize = 0
            filename = "mat.txt"
            eps = 0.001
            !isSolutionFounded = 0
        end subroutine init

        subroutine destruct()
            deallocate (coeff)
            deallocate (freeTerm)
            !deallocate (solutions)
        end subroutine destruct

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
            !allocate ( solutions(systemSize) )

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

        subroutine applyJakobyMethod()
            real, dimension (:), allocatable :: TempX
            real, dimension (:), allocatable :: x
            real norm
            integer iterations

            allocate( TempX(systemSize) )
            allocate( x(systemSize) )
            norm = 1
            iterations = 0

            print *, "---------JAKOBY--------"
            do i = 1, systemsize
                x(i) = 1
            end do

            do while (norm > eps)
                do i = 1,systemSize
                    TempX(i) = freeTerm(i)
                    do g = 1,systemSize
                        if (i .NE. g) then
                            TempX(i) = TempX(i) - coeff(i,g) * x(g)
                        end if
                    end do
                    TempX(i) = TempX(i) / coeff(i,i)
                end do

                norm = abs(x(1) - TempX(1))
                do h = 1,systemSize
                    if (abs(x(h) - TempX(h)) > norm) then
                        norm = abs(x(h) - TempX(h))
                    end if
                    x(h) = TempX(h)
                end do

                iterations = iterations + 1
            end do

            print *, "Iterations:    ",iterations
            do i = 1,systemSize
                print *, "x[",i,"] = ",x(i)
            end do

            deallocate (TempX)
            deallocate (x)
        end subroutine applyJakobyMethod

        integer function converge(x,p)
            real x(systemSize)
            real p(systemSize)
            real, dimension (:), allocatable :: diff
            real maxValue

            allocate( diff(systemSize) )
            maxValue = 0

            do i = 1,systemSize
                diff(i) = abs(x(i) - p(i))
                if (diff(i) > maxValue) then
                    maxValue = diff(i)
                end if
            end do

            deallocate(diff)
            if (maxValue < eps) then
                converge = 1
            else if (maxValue .GE. eps) then
                converge = 0
            end if
        end function converge

        subroutine applyZeidelMethod()
            real, dimension (:), allocatable :: p
            real, dimension (:), allocatable :: x
            real, dimension (:), allocatable :: diff
            real var
            integer iterations, tmp, maxValue


            allocate( p(systemSize) )
            allocate( x(systemSize) )
            allocate( diff(systemSize) )
            iterations = 0

            print *, "---------ZEIDEL--------"
            do i = 1, systemsize
                x(i) = 1
            end do

            maxValue = 1
            tmp = converge(x,p)
            do while (tmp .EQ. 0)
                do i = 1,systemSize
                    p(i) = x(i)
                end do

                do i = 1, systemSize
                    var = 0
                    do j = 1,i-1
                        var = var + (coeff(i,j) * x(j))
                    end do
                    do j = i+1,systemSize
                        var = var + (coeff(i,j) * p(j))
                    end do
                    x(i) = (freeTerm(i) - var) / coeff(i,i)
                end do

                !print *, x(1)

                iterations = iterations + 1
                tmp = converge(x,p)
            end do

            print *, "Iterations:    ",iterations
            do i = 1,systemSize
                print *, "x[",i,"] = ",x(i)
            end do

            deallocate (p)
            deallocate (x)
        end subroutine applyZeidelMethod
end module


program main
    use gaussMethod
    call init()

    call initializationFromFile()
    !call showSystem()
    call applyJakobyMethod()
    call applyZeidelMethod()
    !call showSolutions()

    call destruct()
end
