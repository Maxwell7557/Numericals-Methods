classdef ZeidelAndJakoby < handle
    properties
        systemSize = 0;
        coeff = [];
        freeTerm = [];
        filename = "mat.txt";
        eps = 0.001
    end
        
     methods        
        function initializationFromFile(obj)
            m = 0;
            n = 0;
            
            obj.coeff = load(obj.filename);
            [m,n] = size(obj.coeff);
            
            obj.systemSize = m;

            for i = 1:obj.systemSize
               obj.freeTerm(i,1) = obj.coeff(i,obj.systemSize+1);
               obj.coeff(i,obj.systemSize+1) = 0;
            end
            obj.coeff(:,obj.systemSize+1)=[];
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
        
        function applyJakobyMethod(obj)
            disp("--------JAKOBY--------");
            diagonalElements = diag(obj.coeff);
            diagonalMatrix = diag(diagonalElements);
            inversedDiagMatrix = inv(diagonalMatrix);
            
            b = eye(obj.systemSize) - mtimes(inversedDiagMatrix,obj.coeff);
            d = mtimes(inversedDiagMatrix,obj.freeTerm);
            
            currentX = zeros(obj.systemSize,1);
            previousX = ones(obj.systemSize,1);
            
            iterations = 0;
            while 1
               currentX = mtimes(b,previousX) + d;
               diff = currentX - previousX;
               
               maxValue = abs(max(diff));
               minValue = abs(min(diff));
               
               iterations = iterations + 1; 
               if ( (maxValue > minValue) && (maxValue < obj.eps) )
                   fprintf ('Iterations:   %d\n',iterations);
                   for i = 1:obj.systemSize
                       fprintf ('x[%d] = %f\n',i, currentX(i));
                   end
                   break;
               elseif ( (maxValue < minValue) && (minValue < obj.eps) )
                   fprintf ('Iterations:   %d\n',iterations);
                   for i = 1:obj.systemSize
                       fprintf ('x[%d] = %f\n',i, currentX(i));
                   end
                   break;
               end
               
               previousX = currentX;
            end
            
        end

        function applyZeidelMethod(obj)
            disp("--------ZEIDEL--------");
            triangMatrix = tril(obj.coeff);
            inversedMatrix = inv(triangMatrix);
            
            b = eye(obj.systemSize) - mtimes(inversedMatrix,obj.coeff);
            d = mtimes(inversedMatrix,obj.freeTerm);
            
            currentX = zeros(obj.systemSize,1);
            previousX = ones(obj.systemSize,1);
            
            iterations = 0;
            while 1
               currentX = mtimes(b,previousX) + d;
               diff = currentX - previousX;
               
               maxValue = abs(max(diff));
               minValue = abs(min(diff));
               
               iterations = iterations + 1; 
               if ( (maxValue > minValue) && (maxValue < obj.eps) )
                   fprintf ('Iterations:   %d\n',iterations);
                   for i = 1:obj.systemSize
                       fprintf ('x[%d] = %f\n',i, currentX(i));
                   end
                   break;
               elseif ( (maxValue < minValue) && (minValue < obj.eps) )
                   fprintf ('Iterations:   %d\n',iterations);
                   for i = 1:obj.systemSize
                       fprintf ('x[%d] = %f\n',i, currentX(i));
                   end
                   break;
               end
               
               previousX = currentX;
            end
        end
        
    end
end