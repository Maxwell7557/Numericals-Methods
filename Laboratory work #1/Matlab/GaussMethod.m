classdef GaussMethod < handle
    properties
        systemSize = 0;
        coeff = [];
        freeTerm = [];
        solutions = [];
        filename = "matrix.txt";
        isSolutionFounded = 0;
    end
        
     methods
        function initializationFromUserInput(obj)
            obj.systemSize = input('Введите количество уравнений: ');
            
            obj.coeff(obj.systemSize,obj.systemSize+1) = 0;
            obj.freeTerm(obj.systemSize) = 0;
            obj.solutions(obj.systemSize) = 0;
            
            disp('Введите коэффициенты: ');
            for i = 1:obj.systemSize
               fprintf('%d уравнение: \n', i); 
               for j = 1:obj.systemSize
                   fprintf('a[%d][%d]= ', i,j);
                   obj.coeff(i,j) = input(' ');
               end
               
               disp('Введите свободный член: ');
               fprintf('y[%d]= ',i);
               obj.freeTerm(i) = input(' ');
            end
            
        end
        
        function initializationFromFile(obj)
            m = 0;
            n = 0;
            
            obj.coeff = load(obj.filename);
            [m,n] = size(obj.coeff);
            
            obj.systemSize = m;
            obj.freeTerm(obj.systemSize) = 0;
            obj.solutions(obj.systemSize) = 0;
            
            for i = 1:obj.systemSize
               obj.freeTerm(i) = obj.coeff(i,obj.systemSize+1);
               obj.coeff(i,obj.systemSize+1) = 0;
            end
        end
        
        function showSystem(obj)
            disp('Введённая система:');
            for i = 1:obj.systemSize
                for j = 1:obj.systemSize
                    if (obj.coeff(i,j) < 0)
                        fprintf('(%f)*x%d',obj.coeff(i,j), j);
                    else
                        fprintf('%f*x%d',obj.coeff(i,j), j);
                    end
                    
                    if (j < obj.systemSize)
                        fprintf( ' + ' );
                    end 
                end
                fprintf(' = %d\n', obj.freeTerm(i));
            end
        end
        
        function applyGaussMethod(obj)
            x(obj.systemSize) = 0;
            c = 0.0;

            rankOfUsualMatrix = rank(obj.coeff);
            
            for i = 1: obj.systemSize
                obj.coeff(i,obj.systemSize+1) = obj.freeTerm(i);
            end 
            
            rankOfExtendedMatrix = rank(obj.coeff);

            if (rankOfUsualMatrix < rankOfExtendedMatrix)
                disp ('Нет решения!');
                return
            elseif ( (rankOfUsualMatrix == rankOfExtendedMatrix) && (rankOfUsualMatrix < obj.systemSize) )
                disp ('Система имеет бесконечное множество решений');
                return
            elseif ( (rankOfUsualMatrix == rankOfExtendedMatrix) && (rankOfUsualMatrix == obj.systemSize) )
                disp ('Система имеет единственное решение!');
            end
            
            for i = obj.systemSize: -1:2
                if (obj.coeff(i-1,1) < obj.coeff(i,1))
                    for j = 1:obj.systemSize+1
                        c = obj.coeff(i,j);
                        obj.coeff(i,j) = obj.coeff(i-1,j);
                        obj.coeff(i-1,j) = c;
                    end 
                end 
            end
            
            for k = 1:obj.systemSize-1
                for i = k: obj.systemSize-1
                    c = (obj.coeff(i+1,k) / obj.coeff(k,k));

                    for j = 1:obj.systemSize+1
                        obj.coeff(i+1,j) = obj.coeff(i+1,j) - c*obj.coeff(k,j);
                    end
                end
            end
            
            for i = obj.systemSize: -1: 1
                c = 0;
                for j = i:obj.systemSize
                    c = c + obj.coeff(i,j) * x(j);
                end

                x(i) = (obj.coeff(i,obj.systemSize+1) - c) / obj.coeff(i,i);
            end

            obj.isSolutionFounded = 1;

            for i = 1:obj.systemSize
                obj.solutions(i) = x(i);
            end

        end

        function showSolutions(obj)
            if (obj.isSolutionFounded == 1)
                disp('Решение введённой системы:');
                for i = 1:obj.systemSize
                    fprintf ('x[%d]=%f\n',i, obj.solutions(i));
                end
            end
        end 

        
    end
end