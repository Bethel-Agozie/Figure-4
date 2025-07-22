% Define the piece-wise function for the 1-wave

function f1 = piecewise1_Vm_gen(x,A1,m1,n)

    if x(1) <= A1       % rarefaction
        fun = @(A) sqrt(m1*A.^(m1-2)+n*A.^(-n-2));
        f1 = integral(fun,A1,x(1),'ArrayValued',true);

    else                % shock
        g1=sqrt((x(1)-A1)/(A1*x(1)));
        h1=sqrt(m1/(m1+1)*(x(1).^(m1+1)-A1.^(m1+1))-n/(n-1)*(x(1).^(-n+1)-A1.^(-n+1)));
        f1 = g1*h1;
    end
end
