function [exitflag,guess] = resonance_So0oR_Vn_gen(u1,u2,k,m,n1,n2,A2bar,A1,A2,guess,A22_guess)

    exitflag = 0;
    options = optimoptions('fsolve', 'FunctionTolerance', 1e-8, 'OptimalityTolerance', 1e-8);
       
    % Find A22
    [A22,~,~]=fsolve(@find_A22,A22_guess,options);
    c22=sqrt(k*(m*(A22/A2bar).^m+n2*(A22/A2bar).^(-n2)));
    u22=rRout(A22);

   % Define functions to find_A22 and rRout
    function u=rRout(A)
       f1=@(a) (k*(m*a.^(m-2)/A2bar.^m+n2*a.^(-n2-2)/A2bar.^(-n2))).^(1/2);
      int1=integral(f1,A2,A);
    u=u2+int1;
    end

    function v=find_A22(A)
    c22=sqrt(k*(m*(A/A2bar).^m+n2*(A/A2bar).^(-n2)));
    u22=rRout(A);
    v=u22+c22;
    end

    % Initial celerities
    c1=sqrt(m*A1.^m+n1*A1.^(-n1));
    c2=sqrt(k.*(m*(A2/A2bar).^m+n2*(A2/A2bar).^(-n2)));

 % Define the cases of the system of nonlinear equations and apply the initial guess
    % Case 1
    [S,~,exit1]=fsolve(@So0oR1,guess,options);
    ns=S(1);
    A11=S(2);
    A12=S(3);
    A21=S(4);

    u21=u22*A22/A21;
    c21=sqrt(k*(m*(A21/A2bar).^m+ns*(A21/A2bar).^(-ns)));
    u12a=u21*A21/A12;
    u12b=u21+sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    c12=sqrt(k*(m*A12.^m+ns*A12.^(-ns)));
    u11a=A12*u12/A11;
    u11b=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;
    c11=sqrt((m*A11.^m+n1*A11.^(-n1)));
    lS=(A1*u1-A11*u11)/(A1-A11);
    
    
    if exit1 > 0 && (u11-c11)<lS && lS<(u1-c1) && (u22+c22) < (u2+c2) && u12-c12<0 && u21+c21<0 && A1<A11 &&...
            A2>A22 && ns<=n2 && ns>=n1 && isreal(A11) && isreal(A12) && isreal(A21)&& isreal(ns)
        disp(['exit1: ' num2str(exit1)]);
        exitflag = 1;
        guess = [ns,A11,A12,A21];
    end

    % Case2
    [S,~,exit2]=fsolve(@So0oR2,guess,options);
    ns=S(1);
    A11=S(2);
    A12=S(3);
    A21=S(4);

    u21=u22*A22/A21;
    c21=sqrt(k*(m*(A21/A2bar).^m+ns*(A21/A2bar).^(-ns)));
    u12a=u21*A21/A12;
    u12b=u21+sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    c12=sqrt(k*(m*A12.^m+ns*A12.^(-ns)));
    u11a=A12*u12/A11;
    u11b=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;
    c11=sqrt((m*A11.^m+n1*A11.^(-n1)));
    lS=(A1*u1-A11*u11)/(A1-A11);
    
    if exit2 > 0 && (u11-c11)<lS && lS<(u1-c1) && (u22+c22) < (u2+c2) && u12-c12<0 && u21+c21<0 && A1<A11 &&...
            A2>A22 && ns<=n2 && ns>=n1 && isreal(A11) && isreal(A12) && isreal(A21)&& isreal(ns)
       disp(['exit2: ' num2str(exit2)]);
        exitflag = 1;
       guess = [ns,A11,A12,A21];
    end

   % Case 3
   [S,~,exit3]=fsolve(@So0oR3,guess,options);
    ns=S(1);
    A11=S(2);
    A12=S(3);
    A21=S(4);

    u21=u22*A22/A21;
    c21=sqrt(k*(m*(A21/A2bar).^m+ns*(A21/A2bar).^(-ns)));
    u12a=u21*A21/A12;
    u12b=u21-sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    c12=sqrt(k*(m*A12.^m+ns*A12.^(-ns)));
    u11a=A12*u12/A11;
    u11b=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;
    c11=sqrt((m*A11.^m+n1*A11.^(-n1)));
    lS=(A1*u1-A11*u11)/(A1-A11);
   
   if exit3 > 0 && (u11-c11)<lS && lS<(u1-c1) && (u22+c22) < (u2+c2) && u12-c12<0 && u21+c21<0 && A1<A11 &&...
           A2>A22 && ns<=n2 && ns>=n1 && isreal(A11) && isreal(A12) && isreal(A21)&& isreal(ns)
       disp(['exit3: ' num2str(exit3)]);
       exitflag = 1;
       guess = [ns,A11,A12,A21];
   end

   % Case 4
   [S,~,exit4]=fsolve(@So0oR4,guess,options);
    ns=S(1);
    A11=S(2);
    A12=S(3);
    A21=S(4);

    u21=u22*A22/A21;
    c21=sqrt(k*(m*(A21/A2bar).^m+ns*(A21/A2bar).^(-ns)));
    u12a=u21*A21/A12;
    u12b=u21-sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    c12=sqrt(k*(m*A12.^m+ns*A12.^(-ns)));
    u11a=A12*u12/A11;
    u11b=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;
    c11=sqrt((m*A11.^m+n1*A11.^(-n1)));
    lS=(A1*u1-A11*u11)/(A1-A11);
   
   if exit4 > 0 && (u11-c11)<lS && lS<(u1-c1) && (u22+c22) < (u2+c2) && u12-c12<0 && u21+c21<0 && A1<A11 &&...
           A2>A22 && ns<=n2 && ns>=n1 && isreal(A11) && isreal(A12) && isreal(A21)&& isreal(ns)
       disp(['exit4: ' num2str(exit4)]);
       exitflag = 1;
       guess = [ns,A11,A12,A21];
   end
   

   % Equations to solve
   function v=So0oR1(S)
        ns=S(1);
        A11=S(2);
        A12=S(3);
        A21=S(4);
    
    u21=u22*A22/A21;
    u12a=u21*A21/A12;
    u12b=u21+sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    u11a=A12*u12/A11;
    u11b=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;

        v(1)=1/2*u11^2+(A11.^m-A11.^(-n1))-1/2*u12^2-k*(A12.^m-A12.^(-ns));
        v(2)=u11a-u11b;
        v(3)=u12a-u12b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^m-(A21/A2bar).^(-ns))-1/2*u22^2-k*((A22/A2bar).^m-(A22/A2bar).^(-n2));
    end

    function v=So0oR2(S)
        ns=S(1);
        A11=S(2);
        A12=S(3);
        A21=S(4);
    
    u21=u22*A22/A21;
    u12a=u21*A21/A12;
    u12b=u21+sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    u11a=A12*u12/A11;
    u11b=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;

        v(1)=1/2*u11^2+(A11.^m-A11.^(-n1))-1/2*u12^2-k*(A12.^m-A12.^(-ns));
        v(2)=u11a-u11b;
        v(3)=u12a-u12b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^m-(A21/A2bar).^(-ns))-1/2*u22^2-k*((A22/A2bar).^m-(A22/A2bar).^(-n2));
    end


    function v=So0oR3(S)
        ns=S(1);
        A11=S(2);
        A12=S(3);
        A21=S(4);
    
    u21=u22*A22/A21;
    u12a=u21*A21/A12;
    u12b=u21-sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    u11a=A12*u12/A11;
    u11b=u1+sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;

        v(1)=1/2*u11^2+(A11.^m-A11.^(-n1))-1/2*u12^2-k*(A12.^m-A12.^(-ns));
        v(2)=u11a-u11b;
        v(3)=u12a-u12b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^m-(A21/A2bar).^(-ns))-1/2*u22^2-k*((A22/A2bar).^m-(A22/A2bar).^(-n2));
    end

    function v=So0oR4(S)
        ns=S(1);
        A11=S(2);
        A12=S(3);
        A21=S(4);
    
    u21=u22*A22/A21;
    u12a=u21*A21/A12;
    u12b=u21-sqrt(k*(m*(A12.^(m+1)-A21.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A21.^(-ns+1)-A12.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A12-A21)/(A12*A21));
    u12=u12a;
    u11a=A12*u12/A11;
    u11b=u1-sqrt((m*(A1.^(m+1)-A11.^(m+1))/(m+1)+...
        n1*(A11.^(-n1+1)-A1.^(-n1+1))/(n1-1))*(A1-A11)/(A1*A11));
    u11=u11a;

        v(1)=1/2*u11^2+(A11.^m-A11.^(-n1))-1/2*u12^2-k*(A12.^m-A12.^(-ns));
        v(2)=u11a-u11b;
        v(3)=u12a-u12b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^m-(A21/A2bar).^(-ns))-1/2*u22^2-k*((A22/A2bar).^m-(A22/A2bar).^(-n2));
    end


% Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('So0oR');
end
