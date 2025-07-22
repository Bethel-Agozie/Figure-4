function [exitflag,guess] = resonance_Ro0oR_Vn_gen(u1,u2,k,m,n1,n2,A2bar,A1,A2,guess,A11_guess)

    exitflag = 0;

    % Find A11
    [A11,~,~]=fsolve(@find_A11,A11_guess);
    c11=sqrt((m*A11.^m+n1*A11.^(-n1)));
    u11=lRout(A11);

     % Define functions to find_A11 and lRout 
    function u=lRout(A)
       f1=@(a) (m*a.^(m-2)+n1*a.^(-(n1+2))).^(1/2);
        int1=quad(f1,A1,A);
    u=u1-int1;
    end
  
    function v=find_A11(A)
    c11=sqrt(m*(A)^m+n1*(A)^(-n1));
    u11=lRout(A);
    v=u11-c11;
    end
    
    % Initial celerities
    c1=sqrt(m*A1.^m+n1*A1.^(-n1));
    c2=sqrt(k.*(m*(A2/A2bar).^m+n2*(A2/A2bar).^(-n2)));

    % Define the cases of the system of nonlinear equations and apply the initial guess

    % Case 1
    [S,~,exit1]=fsolve(@Ro0oR1,guess);
    ns=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    c11=sqrt(m*A11.^m+n1*A11.^(-n1));
    c12=sqrt(k*(m*A12.^m+ns*A12.^(-ns)));
    u21a=A12*u12/A21;
    u21b=u12+sqrt(k*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A12.^(-ns+1)-A21.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(k*(m*(A21/A2bar).^m+ns*(A21/A2bar).^(-ns)));
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k*(m*A.^(m-2)/A2bar.^m+n2*A.^(-n2-2)/A2bar.^(-n2)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;
    c22=sqrt(k*(m*(A22/A2bar).^m+n2*(A22/A2bar).^(-n2)));

    if exit1 > 0 && (u11-c11)>(u1-c1) && u12-c12>0 && u21+c21>0 && A1>=A11 && A2>=A22 && (u22b+c22)<(u2+c2) &&...
            ns<=n2 && ns>=n1 && A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(ns)   
        disp(['exit1: ' num2str(exit1)]);
        exitflag = 1;
        guess = [ns,A12,A21,A22];
    end

    % Case2
    [S,~,exit2]=fsolve(@Ro0oR2,guess);
    ns=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    c11=sqrt(m*A11.^m+n1*A11.^(-n1));
    c12=sqrt(k*(m*A12.^m+ns*A12.^(-ns)));
    u21a=A12*u12/A21;
    u21b=u12-sqrt(k*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A12.^(-ns+1)-A21.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(k*(m*(A21/A2bar).^m+ns*(A21/A2bar).^(-ns)));
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k*(m*A.^(m-2)/A2bar.^m+n2*A.^(-n2-2)/A2bar.^(-n2)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;
    c22=sqrt(k*(m*(A22/A2bar).^m+n2*(A22/A2bar).^(-n2)));

    if exit2 > 0 && (u11-c11)>(u1-c1) && u12-c12>0 && u21+c21>0 && A1>=A11 && A2>=A22 && (u22b+c22)<(u2+c2) &&...
            ns<=n2 && ns>=n1 && A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(ns)
       disp(['exit2: ' num2str(exit2)]);
       exitflag = 1;
       guess = [ns,A12,A21,A22];
    end
   

   % Equations to solve
   function v=Ro0oR1(S)
    ns=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12+sqrt(k*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A12.^(-ns+1)-A21.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k*(m*A.^(m-2)/A2bar.^m+n2*A.^(-n2-2)/A2bar.^(-n2)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;

        v(1)=1/2*u11^2+(A11.^m-A11.^(-n1))-1/2*u12^2-k*(A12.^m-A12.^(-ns));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^m-(A21/A2bar).^(-ns))-1/2*u22^2-k*((A22/A2bar).^m-(A22/A2bar).^(-n2));
    end

    function v=Ro0oR2(S)
    ns=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12-sqrt(k*(m*(A21.^(m+1)-A12.^(m+1))/((m+1)*A2bar.^m)+...
        ns*(A12.^(-ns+1)-A21.^(-ns+1))/((ns-1)*A2bar.^(-ns)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    fun = @(A) sqrt(k*(m*A.^(m-2)/A2bar.^m+n2*A.^(-n2-2)/A2bar.^(-n2)));
       u22b=u2+ quad(fun,A2,A22);
    u22=u22a;

        v(1)=1/2*u11^2+(A11.^m-A11.^(-n1))-1/2*u12^2-k*(A12.^m-A12.^(-ns));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^m-(A21/A2bar).^(-ns))-1/2*u22^2-k*((A22/A2bar).^m-(A22/A2bar).^(-n2));
    end

 % Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('Ro0oR');
end
