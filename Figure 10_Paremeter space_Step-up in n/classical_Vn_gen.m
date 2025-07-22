function [x, exitflag]= classical_Vn_gen(u1,u2,k,m,A2bar,A1,A2,n1,n2)

% Define the system of nonlinear equations
fun = @(x) [
    x(2) - u1 + piecewise1_Vn_gen(x,A1,m,n1);
    x(1)*x(2) - x(3)*x(4);
    0.5*(x(2)^2 - x(4)^2) + x(1).^m-x(1).^(-n1)-k*((x(3)/A2bar)^m-(x(3)/A2bar)^(-n2));
    x(4)- u2 - piecewise2_Vn_gen(x,A2,m,k,n2,A2bar)
];

% Initial guess
% for m=10
x0 = [1; 1; 0.8; 1.2]; 

% for m=0.5
%x0=[2.5901  -0.5  2.1551  -0.45];

% Solve the system of equations using fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'MaxIterations', 10000, 'FunctionTolerance', 1e-8);
[x,fval,exitflag] = fsolve(fun, x0, options);

% Display the result
disp('Solution:');
disp(x);
disp(exitflag);

% Decide the the type of solution
if exitflag==1 && x(1)>0 && x(3)>0  

    if x(1) <= A1 && x(3) <= A2
        disp('RoR')
        
    elseif x(1) <= A1 && x(3) > A2  
           disp('RoS')
                       
    elseif x(1) > A1 && x(3) <= A2
        disp('SoR')
                
    elseif x(1) > A1 && x(3) > A2
        disp('SoS')
              
    end
else
    disp('Not solved')
end


