function [exitflag,guess] = resonance_Ro0oS_Vm_gen(u1,u2,k,m1,m2,n,A2bar,A1,A2,guess,A11_guess)

    exitflag = 0;

     % Find A11
    [A11,~,~]=fsolve(@find_A11,A11_guess);
    c11=sqrt((m1*A11.^m1+n*A11.^(-n)));
    u11=lRout(A11);
    
     % Define functions to find_A11 and lRout 
    function u=lRout(A)
       f1=@(a) (m1*a.^(m1-2)+n*a.^(-(n+2))).^(1/2);
        int1=quad(f1,A1,A);
    u=u1-int1;
    end
  
    function v=find_A11(A)
    c11=sqrt(m1*(A)^m1+n*(A)^(-n));
    u11=lRout(A);
    v=u11-c11;
    end

    % Initial celerities
    c1=sqrt(m1*A1.^m1+n*A1.^(-n));
    c2=sqrt(k.*(m2*(A2/A2bar).^m2+n*(A2/A2bar).^(-n)));

    % Define the cases of the system of nonlinear equations and apply the initial guess

    % Case 1
    [S,~,exit1]=fsolve(@Ro0oS1,guess);

    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
  
    u12=A11*u11/A12;
    c12=sqrt(k*(ms*A12.^ms+n*A12.^(-n)));
    u21a=A12*u12/A21;
    u21b=u12+sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(k*(ms*(A21/A2bar).^ms+n*(A21/A2bar).^(-n)));
    u22a=A21*u21/A22;
    u22b=u2+sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
    c22=sqrt(k*(m2*(A22/A2bar).^m2+n*(A22/A2bar).^(-n)));
    rS=(A2*u2-A22*u22)/(A2-A22);
    
    if exit1 > 0 && (u11-c11)>(u1-c1) &&  A1>=A11 && A22>A2 && (u2+c2)<rS && rS<(u22+c22) && ms<=m2 && ms>=m1 &&...
            A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(ms)&& u12-c12>0 && u21+c21>0
        disp(['exit1: ' num2str(exit1)]);
        exitflag = 1;
        guess = [ms,A12,A21,A22];
    end


    % Case2
    [S,~,exit2]=fsolve(@Ro0oS2,guess);
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
  
    u12=A11*u11/A12;
    c12=sqrt(k*(ms*A12.^ms+n*A12.^(-n)));
    u21a=A12*u12/A21;
    u21b=u12+sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(k*(ms*(A21/A2bar).^ms+n*(A21/A2bar).^(-n)));
    u22a=A21*u21/A22;
    u22b=u2-sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
    c22=sqrt(k*(m2*(A22/A2bar).^m2+n*(A22/A2bar).^(-n)));
    rS=(A2*u2-A22*u22)/(A2-A22);
   
    if exit2 > 0 && (u11-c11)>(u1-c1) &&  A1>=A11 && A22>A2 && (u2+c2)<rS && rS<(u22+c22) && ms<=m2 && ms>=m1 &&...
            A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(ms)&& u12-c12>0 && u21+c21>0
        disp(['exit2: ' num2str(exit2)]);
        exitflag = 1;
        guess = [ms,A12,A21,A22];
    end

   % Case 3
   [S,~,exit3]=fsolve(@Ro0oS3,guess);
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
  
    u12=A11*u11/A12;
    c12=sqrt(k*(ms*A12.^ms+n*A12.^(-n)));
    u21a=A12*u12/A21;
    u21b=u12-sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(k*(ms*(A21/A2bar).^ms+n*(A21/A2bar).^(-n)));
    u22a=A21*u21/A22;
    u22b=u2+sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
    c22=sqrt(k*(m2*(A22/A2bar).^m2+n*(A22/A2bar).^(-n)));
    rS=(A2*u2-A22*u22)/(A2-A22);
    
    if exit3 > 0 && (u11-c11)>(u1-c1) &&  A1>=A11 && A22>A2 && (u2+c2)<rS && rS<(u22+c22) && ms<=m2 && ms>=m1 &&...
            A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(ms)&& u12-c12>0 && u21+c21>0
        disp(['exit3: ' num2str(exit3)]);
        exitflag = 1;
        guess = [ms,A12,A21,A22];
    end


   % Case 4
   [S,~,exit4]=fsolve(@Ro0oS4,guess);
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
  
    u12=A11*u11/A12;
    c12=sqrt(k*(ms*A12.^ms+n*A12.^(-n)));
    u21a=A12*u12/A21;
    u21b=u12-sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    c21=sqrt(k*(ms*(A21/A2bar).^ms+n*(A21/A2bar).^(-n)));
    u22a=A21*u21/A22;
    u22b=u2-sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
    c22=sqrt(k*(m2*(A22/A2bar).^m2+n*(A22/A2bar).^(-n)));
    rS=(A2*u2-A22*u22)/(A2-A22);
   
    if exit4 > 0 && (u11-c11)>(u1-c1) &&  A1>=A11 && A22>A2 && (u2+c2)<rS && rS<(u22+c22) && ms<=m2 && ms>=m1 &&...
            A11<=A22 && isreal(A12) && isreal(A21) && isreal(A22)&& isreal(ms)&& u12-c12>0 && u21+c21>0
        disp(['exit4: ' num2str(exit4)]);
        exitflag = 1;
        guess = [ms,A12,A21,A22];
    end
         
   % Equations to solve
   function v=Ro0oS1(S)
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12+sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    u22b=u2+sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
        v(1)=1/2*u11^2+(A11.^m1-A11.^(-n))-1/2*u12^2-k*(A12.^ms-A12.^(-n));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^ms-(A21/A2bar).^(-n))-1/2*u22^2-k*((A22/A2bar).^m2-(A22/A2bar).^(-n));

    end

    function v=Ro0oS2(S)
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12+sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    u22b=u2-sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
        v(1)=1/2*u11^2+(A11.^m1-A11.^(-n))-1/2*u12^2-k*(A12.^ms-A12.^(-n));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^ms-(A21/A2bar).^(-n))-1/2*u22^2-k*((A22/A2bar).^m2-(A22/A2bar).^(-n));
    end

    function v=Ro0oS3(S)
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12-sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    u22b=u2+sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
        v(1)=1/2*u11^2+(A11.^m1-A11.^(-n))-1/2*u12^2-k*(A12.^ms-A12.^(-n));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^ms-(A21/A2bar).^(-n))-1/2*u22^2-k*((A22/A2bar).^m2-(A22/A2bar).^(-n));
    end

    function v=Ro0oS4(S)
    ms=S(1);
    A12=S(2);
    A21=S(3);
    A22=S(4);
    
    u12=A11*u11/A12;
    u21a=A12*u12/A21;
    u21b=u12-sqrt(k*(ms*(A21.^(ms+1)-A12.^(ms+1))/((ms+1)*A2bar.^ms)+...
        n*(A12.^(-n+1)-A21.^(-n+1))/((n-1)*A2bar.^(-n)))*(A21-A12)/(A21*A12));
    u21=u21a;
    u22a=A21*u21/A22;
    u22b=u2-sqrt(k*(m2*(A2.^(m2+1)-A22.^(m2+1))/((m2+1)*A2bar.^m2)+...
        n*(A22.^(-n+1)-A2.^(-n+1))/((n-1)*A2bar.^(-n)))*(A2-A22)/(A2*A22));
    u22=u22a;
        v(1)=1/2*u11^2+(A11.^m1-A11.^(-n))-1/2*u12^2-k*(A12.^ms-A12.^(-n));
        v(2)=u21a-u21b;
        v(3)=u22a-u22b;
        v(4)=1/2*u21^2+k*((A21/A2bar).^ms-(A21/A2bar).^(-n))-1/2*u22^2-k*((A22/A2bar).^m2-(A22/A2bar).^(-n));
    end

% Display the result
    disp('Resonant soln:');
    disp(guess);
    disp('Ro0oS');
end
