% Define the piece-wise function for the 1-wave

function f1 = piecewise1_Vn_gen(x,A1,m,n1)

    if x(1) <= A1       % rarefaction
        fun = @(A) sqrt(m*A.^(m-2)+n1*A.^(-n1-2));
        f1 = integral(fun,A1,x(1),'ArrayValued',true);

    else                % shock
        g1=sqrt((x(1)-A1)/(A1*x(1)));
        h1=sqrt(m/(m+1)*(x(1).^(m+1)-A1.^(m+1))-n1/(n1-1)*(x(1).^(-n1+1)-A1.^(-n1+1)));
        f1 = g1*h1;
    end
end
