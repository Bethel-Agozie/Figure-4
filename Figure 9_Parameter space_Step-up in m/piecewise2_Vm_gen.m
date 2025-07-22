% Define the piece-wise function for the 2-wave

function f2 = piecewise2_Vm_gen(x,A2,k,m2,n,A2bar)

    if x(3) <= A2       % rarefaction
        fun = @(A) sqrt(k*(m2*(A/A2bar).^(m2-2)+n*(A/A2bar).^(-n-2)));
        f2 = integral(fun,A2,x(3),'ArrayValued',true);

    else                % shock
        g=sqrt((x(3)-A2)/(x(3)*A2));
        h=sqrt(k*((m2/((m2+1)*A2bar.^m2)*(x(3).^(m2+1)-A2.^(m2+1)))- ...
            (n/((n-1)*A2bar.^(-n))*(x(3).^(-n+1)-A2.^(-n+1)))));
        f2 = g*h;
    end
end