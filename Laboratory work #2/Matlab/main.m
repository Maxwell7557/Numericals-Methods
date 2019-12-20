function main()
    a = 0.0;
    b = 1.0;
    eps1 = 0.0001;
    eps2 = 0.00001;
    
    fprintf("a = %d \n", a);
    fprintf("b = %d \n", b);
    fprintf("eps = %d \n", eps1);
    calculateIntegral(a,b,eps1);
    fprintf("\n");
    fprintf("a = %d \n", a);
    fprintf("b = %d \n", b);
    fprintf("eps = %d \n", eps2);
    calculateIntegral(a,b,eps2);
end

function val = func(point)
    val = ( 1/( exp(point) + exp(-point) ));
end

function val = simpson(n,a,b)
    h = (b - a)/n;
    f0 = func(a);
    fn = func(b);
    firstSum = 0.0;
    secondSum = 0.0;
    m = n/2;

    xTmp = [n+1];
    for i = 1: n+1
        xTmp(i) = ( a + h*(i-1) );
    end
        
    for i = 2: m
        firstSum = firstSum + func(xTmp(2*(i-1)+1) );
    end
    
    for i = 2: m+1
        secondSum = secondSum + func(xTmp(2*(i-1)-1+1) );
    end
     
    val = (h/3)*(f0 + 2*firstSum + 4*secondSum + fn);
end

function val = trapeze(n,a,b)
    h = (b - a)/n;
    f0 = func(a);
    fn = func(b);
    sum = 0;

    xTmp = [n+1];
    for i = 1: n+1
        xTmp(i) = ( a + h*(i-1) );
    end
    
    for i = 2: n
        sum = sum + func(xTmp(i));
    end
    
    val = (h/2)*(f0 + 2*sum + fn);
end

function val = leftRect(n,a,b)
    h = (b - a)/n;
    sum = 0;

    xTmp = [n+1];
    for i = 1: n+1
        xTmp(i) = ( a + h*(i-1) );
    end
    
    for i = 1: n
        sum = sum + func(xTmp(i));
    end
    
    val = h * sum;
end

function val = rightRect(n,a,b)
    h = (b - a)/n;
    sum = 0;

    xTmp = [n+1];
    for i = 1: n+1
        xTmp(i) = ( a + h*(i-1) );
    end
    
    for i = 2: n+1
        sum = sum + func(xTmp(i));
    end
    
    val = h * sum;
end

function val = centralRect(n,a,b)
    f0 = func(a)/2;
    fn = func(b)/2;
    h = (b - a)/n;
    sum = 0;

    xTmp = [n+1];
    for i = 1: n+1
        xTmp(i) = ( a + h*(i-1) );
    end
    
    for i = 2: n
        sum = sum + func(xTmp(i));
    end
    
    val = h * (f0 + sum + fn);
end

function val = checkForCondition(type,a,b,eps)
    n = 2;
    res2 = 0;
    res1 = 0;
    tmp = 1;
    
    while tmp == 1
        if type == "simpson"
            res1 = simpson(n,a,b);
        elseif type == "trapeze"
            res1 = trapeze(n,a,b);
        elseif type == "leftRect"
            res1 = leftRect(n,a,b);
        elseif type == "rightRect"
            res1 = rightRect(n,a,b);
        elseif type == "centralRect"
            res1 = centralRect(n,a,b);
        end
        
        diff = abs(res1 - res2);

        res2 = res1;
        if (type == "simpson"); n = n + 2; else; n = n + 1; end

        if (eps*15 >= diff && type == "simpson")
            fprintf("n = %d", n-2);
            tmp = 0;
        elseif (eps*3 >= diff && type ~= "simpson")
            fprintf("n = %d", n-1);
            tmp = 0;
        end
    end
    val = res1;
end

function calculateIntegral(a,b,eps)
    fprintf(" Simpson method: %f \n", checkForCondition('simpson',a,b,eps));
    fprintf(" Trapeze method: %f \n", checkForCondition('trapeze',a,b,eps));
    fprintf(" Left rectangular method: %f \n", checkForCondition('leftRect',a,b,eps));
    fprintf(" Right rectangular method: %f \n", checkForCondition('rightRect',a,b,eps));
    fprintf(" Central rectangular method: %f \n", checkForCondition('centralRect',a,b,eps));
end