real function func(point)
    real :: point;
    !func = ( 1/( EXP(point) + EXP(-point) ));
    func = ( (0.1*(point**2))/log10(point) );
end function func

module integral
    contains
    real function simpson(n, a, b)
        integer n;
        real b, a, buff;

        real, dimension (:), allocatable :: xTmp;
        real h, secondSum, firstSum
        firstSum = 0.0;
        secondSum = 0.0;
        buff=0;
        m = n/2;
        h = (b - a)/n;

        allocate ( xTmp(n) );
        do i = 0, n
            xTmp(i) = (a + h*(i));
        end do

        do i = 1, m-1
            firstSum = firstSum + func(xTmp(2*(i)));
        end do

        do i = 1, m
            secondSum = secondSum + func(xTmp(2*(i)-1));
        end do

        buff=(h/3)*(func(xTmp(0)) + 2*firstSum + 4*secondSum + func(xTmp(n)));
        deallocate (xTmp);
        simpson = buff;
    end function simpson

    real function trapeze(n,a,b)
        integer n;
        real a,b;

        real, dimension (:), allocatable :: xTmp;
        real h, msum
        msum = 0

        h = (b - a)/n;

        allocate( xTmp(n) );

        do i = 0, n
            xTmp(i) = (a + h*(i));
        end do

        do i = 1, n-1
            msum = msum + func(xTmp(i));
        end do
        trapeze = (h/2)*(func(xTmp(0)) + 2*msum + func(xTmp(n)));
        deallocate(xTmp);
    end function trapeze

    real function leftRect(n,a,b)
        integer n;
        real a,b,h,msum;

        real, dimension (:), allocatable :: xTmp;
        msum = 0.0;

        h = (b - a)/n;
        allocate( xTmp(n) );

        do i = 0, n
            xTmp(i) = (a + h*(i));
        end do

        do i = 0, n-1
            msum = msum + func(xTmp(i));
        end do

        deallocate(xTmp)
        leftRect = h * msum;
    end function leftRect

    real function rightRect(n,a,b)
        integer n;
        real a,b,h,msum;

        real, dimension (:), allocatable :: xTmp;
        msum = 0.0;

        h = (b - a)/n;
        allocate( xTmp(n) );

        do i = 0, n
            xTmp(i) = (a + h*(i));
        end do

        do i = 1, n
            msum = msum + func(xTmp(i));
        end do

        deallocate(xTmp);
        rightRect = h * msum;
    end function rightRect

    real function centralRect(n,a,b)
        integer n;
        real b, a, h, msum, buff;

        real, dimension (:), allocatable :: xTmp;
        msum = 0.0;

        h = (b - a)/n;

        allocate( xTmp(n) );

        do i = 0, n
            xTmp(i) = (a + h*(i));
        end do

        do i = 1, n-1
            msum = msum + func(xTmp(i));
        end do
	buff=h * (func(xTmp(0))/2 + msum + func(xTmp(n))/2);
	deallocate(xTmp);
        centralRect = buff
    end function centralRect

    real function checkForCondition(mtype, a, b, eps)
        integer mtype, n, tmp;
        real b, a, eps, res1, res2, diff;
	n=1
	res1=0
	res2=0
	diff=0
	tmp=1

        do while (tmp == 1)
        !do while (eps < diff)
            if (mtype == 1) then
                n=n+1
                res1 = simpson(n,a,b);
                diff = ABS(res2 - res1);
                if ((eps*15) >= diff) then
                    print *,n
                    tmp=0;
                end if
            end if

            if (mtype == 2) then
                res1 = trapeze(n,a,b);
                diff = ABS(res2 - res1);
                if (eps*3 >= diff) then
                    print *,n
                    tmp=0;
                end if
            end if

            if (mtype == 3) then
                res1 = leftRect(n,a,b);
                diff = ABS(res2 - res1);
                if (eps*3 >= diff) then
                    print *,n
                    tmp=0;
                end if
            end if

            if (mtype == 4) then
                res1 = rightRect(n,a,b);
                diff = ABS(res2 - res1);
                if (eps*3 >= diff) then
                    print *,n
                    tmp=0;
                end if
            end if

            if (mtype == 5) then
                res1 = centralRect(n,a,b);
                diff = ABS(res2 - res1);
                if (eps*3 >= diff) then
                    print *,n
                    tmp=0;
                end if
            end if

            res2 = res1;
            n = n + 1;
        end do

        checkForCondition = res1;
    end function checkForCondition

    subroutine calculateIntegral(a, b, eps)
        real a, b, eps
        print *, "Simpson method: n =", checkForCondition(1,a,b,eps);
        print *, "Trapeze method: n = ", checkForCondition(2,a,b,eps);
        print *, "Left rectangular method: n = ", checkForCondition(3,a,b,eps);
        print *, "Right rectangular method: n = ", checkForCondition(4,a,b,eps);
        print *, "Central rectangular method: n = ", checkForCondition(5,a,b,eps);

    end subroutine calculateIntegral
end module integral

program main
    use integral;
    real a,b,eps, eps1;

    a = 2.0
    b = 3.0
    eps = 0.0001
    eps1 = 0.00001

    print *, "a =", a
    print *, "b =", b
    print *, "eps =", 0.0001
    call calculateIntegral(a,b,eps);
    print *," "
    print *, "a =", a
    print *, "b =", b
    print *, "eps =", 0.00001
    call calculateIntegral(a,b,eps1);
end
