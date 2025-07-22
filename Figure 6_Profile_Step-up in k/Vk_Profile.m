%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A profile code for the standard non-dimensional blood flow model Riemann problem using general tube law.
%
%   Bethel fecit, AD 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vk_Profile 
% First, clear the workspace and close any open figures
clear all
clc
close all

% Set precision
format short

% Set initial conditions for the Riemann problem
u1 = 0.2;
A1 = 0.8;

u2 = 0.4;
A2 = 0.88;

A2bar= 1;
k=1.1;
m=10; % Go to classical_Vk_gen.m and choose the appropriate initial guess for the classical case
n=0;
t=0.8;

% Initial guess for the resonant solutions
A11_guess = 1.05;
A12_guess = A11_guess;
guess_Ro0oS = [1.0695    0.4126    0.7675    0.7935];
guess_Ro0oR = [1.0905    0.0993    0.6879    0.7094];
guess_RoSS = [0.4126    0.7935];
guess_RoSR = [0.4126    0.7935];

A22_guess=1.3;
A21_guess=A22_guess;
guess_So0oR= [1.0653    2.075    1.7089    1.45]; 
guess_SSoR= [0.3614    1.0723];

 % Find A11
 [A11,~,~]=fsolve(@find_A11,A11_guess);
    c11=sqrt((m*A11.^m+n*A11.^(-n)));
    u11=lRout(A11);
   
 % Find A22
    [A22,~,~]=fsolve(@find_A22,A22_guess);
    c22=sqrt(k*(m*(A22/A2bar).^m+n*(A22/A2bar).^(-n)));
    u22=rRout(A22);
 
 c2=sqrt(k*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
 
 if u2<c2 && u2>-c2  % Consider only subcritical regime 

 % solve the equations and decide what type of solution
   [exitflag_Ro0oS,guess_Ro0oS] = resonance_Ro0oS_Vk_gen(u1,u2,k,m,n,A2bar,A1,A2,guess_Ro0oS,A11_guess);
   [exitflag_Ro0oR,guess_Ro0oR] = resonance_Ro0oR_Vk_gen(u1,u2,k,m,n,A2bar,A1,A2,guess_Ro0oR,A11_guess);
   [exitflag_RoSS,guess_RoSS] = resonance_RoSS_Vk_gen(u1,u2,k,m,n,A2bar,A1,A2,guess_RoSS,A12_guess);
   [exitflag_RoSR,guess_RoSR] = resonance_RoSR_Vk_gen(u1,u2,k,m,n,A2bar,A1,A2,guess_RoSR,A12_guess);
   [exitflag_So0oR,guess_So0oR] = resonance_So0oR_Vk_gen(u1,u2,k,m,n,A2bar,A1,A2,guess_So0oR,A22_guess);
   [exitflag_SSoR,guess_SSoR] = resonance_SSoR_Vk_gen(u1,u2,k,m,n,A2bar,A1,A2,guess_SSoR,A21_guess);

   [x, exitflag] = classical_Vk_gen(u1,u2,k,m,A2bar,A1,A2,n);
 
%  Compute celerities, shock speeds, decide what type of solution and plot
 
% Plot cross-sectional area, A 
    figure('Position', [234, 150, 350, 240]) % Set the pixel size of the figure
    box on
    if exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            c1=sqrt(m*A1.^m+n*A1.^(-n));
            u12=A11*u11/guess_Ro0oS(2);
            u21=guess_Ro0oS(2)*u12/guess_Ro0oS(3);
            u22=guess_Ro0oS(3)*u21/guess_Ro0oS(4);
            ReS=(guess_Ro0oS(2)*u12-guess_Ro0oS(3)*u21)/(guess_Ro0oS(2)-guess_Ro0oS(3));
            rS=(A2*u2-guess_Ro0oS(4)*u22)/(A2-guess_Ro0oS(4));

            % distances 
            x1=(u1-c1)*t; 
            x11 = (u11-c11)*t;
            x12=x11;
            xReS=ReS*t; 
            x21=xReS;
            x22=x21;
            x2 = rS*t;

               % plots
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x11,20); a=linspace(A1,A11,20);plot(X, a, 'b',LineWidth=1);
        hold on
        X=linspace(x11,x12,20); a=linspace(A11, guess_Ro0oS(2),20); plot(X, a, 'mo',LineWidth=1.5);
        hold on
        X=linspace(x12,xReS,20); a=linspace(guess_Ro0oS(2), guess_Ro0oS(2),20); plot(X, a, 'ko',LineWidth=1.5);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_Ro0oS(2), guess_Ro0oS(3),20); plot(X, a, 'y*',LineWidth=1.5);
        hold on
        X=linspace(xReS,x21,20); a=linspace(guess_Ro0oS(3), guess_Ro0oS(3),20); plot(X, a, 'go',LineWidth=1.5);
        hold on
        X=linspace(x21,x22,20); a=linspace(guess_Ro0oS(3), guess_Ro0oS(4),20); plot(X, a, ':','color', "#7e1e9c",LineWidth=2);
        hold on
        X=linspace(x22,x2,20); a=linspace(guess_Ro0oS(4), guess_Ro0oS(4),20); plot(X, a, 'c',LineWidth=1);
        hold on
        X=linspace(x2,x2,20); a=linspace(guess_Ro0oS(4), A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1);
        legend({'$A_1$','1-$R$', '$lC$','$A_{12}$','$S_{=0}$','$A_{21}$','$rC$',...
                    '$A_{22}$','2-$S$', '$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

     elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
        c1=sqrt(m*A1.^m+n*A1.^(-n));
        u12=A11*u11/guess_Ro0oR(2);
        u21=guess_Ro0oR(2)*u12/guess_Ro0oR(3);
        u22=guess_Ro0oR(3)*u21/guess_Ro0oR(4);
        c22=sqrt(k*(m*(guess_Ro0oR(4)/A2bar).^m+n*(guess_Ro0oR(4)/A2bar).^(-n)));
        ReS=(guess_Ro0oR(2)*u12-guess_Ro0oR(3)*u21)/(guess_Ro0oR(2)-guess_Ro0oR(3));
        c2=sqrt(k.*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
          % distances 
        x1=(u1-c1)*t; 
        x11 = (u11-c11)*t;
        x12=x11;
        xReS=ReS*t;
        x21=xReS;
        x22 = (u22 + c22)*t;
        x2 = (u2 + c2)*t;

            % plots
    X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1.5);
    hold on
    X=linspace(x1,x11,20); a=linspace(A1,A11,20);plot(X, a, 'b',LineWidth=1.5);
    hold on
    X=linspace(x11,x12,20); a=linspace(A11, guess_Ro0oR(2),20); plot(X, a, 'mo',LineWidth=1.5);
    hold on
    X=linspace(x12,xReS,20); a=linspace(guess_Ro0oR(2), guess_Ro0oR(2),20); plot(X, a, 'ko',LineWidth=1.5);
    hold on
    X=linspace(xReS,xReS,20); a=linspace(guess_Ro0oR(2), guess_Ro0oR(3),20); plot(X, a, 'y*',LineWidth=1.5);
    hold on
    X=linspace(xReS,x21,20); a=linspace(guess_Ro0oR(3), guess_Ro0oR(3),20); plot(X, a, 'go',LineWidth=1.5);
    hold on
    X=linspace(x21,x21,20); a=linspace(guess_Ro0oR(3), guess_Ro0oR(4),20); plot(X, a,':','color', "#7e1e9c",LineWidth=2)
    hold on
    X=linspace(x21,x22,20); a=linspace(guess_Ro0oR(4), guess_Ro0oR(4),20); plot(X, a, 'c',LineWidth=1.5);
    hold on
    X=linspace(x22,x2,20); a=linspace(guess_Ro0oR(4), A2,20); plot(X, a, 'r',LineWidth=1.5);
    hold on
    X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1.5);

            legend({'$A_1$','1-$R$', '$lC$','$A_{12}$','$S_{=0}$','$A_{21}$','$rC$',...
                '$A_{22}$','2-$R$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

    elseif exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0  
         A12=A11;
         u12=u11;
         c12=c11;
         c1=sqrt(m*A1.^m+n*A1.^(-n));
        u21=A12*u12/guess_RoSS(1);
        ReS=(guess_RoSS(2)*guess_RoSS(3)-guess_RoSS(1)*u21)/(guess_RoSS(2)-guess_RoSS(1));
        rS=(A2*u2-guess_RoSS(2)*guess_RoSS(3))/(A2-guess_RoSS(2));

       % distances 
        x1=(u1-c1)*t; 
        x12 = (u12-c12)*t;
        x21=x12;
        xReS = ReS*t;
        x2 = rS*t;
       
       % plots
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x12,20); a=linspace(A1,A12,20);plot(X, a, 'b',LineWidth=1);
        hold on
        X=linspace(x12,x21,20); a=linspace(A12,guess_RoSS(1),20);plot(X, a, 'm:',LineWidth=1);
        hold on
        X=linspace(x21,xReS,20); a=linspace(guess_RoSS(1), guess_RoSS(1),20); plot(X, a, 'g',LineWidth=1);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_RoSS(1), guess_RoSS(2),20); plot(X, a, 'color', "#7e1e9c",LineWidth=1);
        hold on
        X=linspace(xReS,x2,20); a=linspace(guess_RoSS(2), guess_RoSS(2),20); plot(X, a, 'c',LineWidth=1);
        hold on
        X=linspace(x2,x2,20); a=linspace(guess_RoSS(2), A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1);
            legend({'$A_1$','1-$R$', '$C$','$A_{21}$','$S_{\neq 0}$',...
                '$A_{22}$','2-$S$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

     elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0
         A12=A11;
         u12=u11;
         c12=c11;
         c1=sqrt(m*A1.^m+n*A1.^(-n));
        u21=A12*u12/guess_RoSR(1);
        c22=sqrt(k*(m*(guess_RoSR(2)/A2bar).^m+n*(guess_RoSR(2)/A2bar).^(-n)));
        ReS=(guess_RoSR(2)*guess_RoSR(3)-guess_RoSR(1)*u21)/(guess_RoSR(2)-guess_RoSR(1));
        c2=sqrt(k.*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
        
        % distances 
        x1=(u1-c1)*t; 
        x12 = (u12-c12)*t;
        xReS = ReS*t;
        x22 = (guess_RoSR(3) + c22)*t; 
        x2=(u2+c2)*t;
       
        % plots
            X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
            hold on
            X=linspace(x1,x12,20); a=linspace(A1,A12,20);plot(X, a, 'b',LineWidth=1);
             hold on
            X=linspace(x12,x12,20); a=linspace(A12,guess_RoSR(1),20);plot(X, a, 'm:',LineWidth=1);
            hold on
            X=linspace(x12,xReS,20); a=linspace(guess_RoSR(1), guess_RoSR(1),20); plot(X, a, 'g',LineWidth=1);
            hold on
            X=linspace(xReS,xReS,20); a=linspace(guess_RoSR(1), guess_RoSR(2),20); plot(X, a,  'color', "#7e1e9c",LineWidth=1);
            hold on
            X=linspace(xReS,x22,20); a=linspace(guess_RoSR(2), guess_RoSR(2),20); plot(X, a, 'c',LineWidth=1);
            hold on
            X=linspace(x22,x2,20); a=linspace(guess_RoSR(2), A2,20); plot(X, a, 'r',LineWidth=1);
            hold on
            X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1); 

            legend({'$A_1$','1-$R$', '$C$','$A_{21}$','$S_{\neq 0}$',...
                '$A_{22}$','2-$R$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
            
    elseif exitflag_So0oR == 1 && guess_So0oR(1)>0 && guess_So0oR(2)>0 && guess_So0oR(3)>0 && guess_So0oR(4)>0
         u21=u22*A22/guess_So0oR(4);
         u12=u21*guess_So0oR(4)/guess_So0oR(3);
         u11=u12*guess_So0oR(3)/guess_So0oR(2);
        lS=(A1*u1-guess_So0oR(2)*u11)/(A1-guess_So0oR(2));
        ReS=(guess_So0oR(3)*u12-guess_So0oR(4)*u21)/(guess_So0oR(3)-guess_So0oR(4));
        c11=sqrt(m*guess_So0oR(2).^m+n*guess_So0oR(2).^(-n));
        % distances 
        x1=lS*t; 
        xReS=ReS*t;
        x12=xReS;
        %x11=x12;
        
        x22=(u22+c22)*t;
        x21=x22;
        x2=(u2+c2)*t;
                
          % plots
    
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x1,20); a=linspace(A1,guess_So0oR(2),20);plot(X, a, 'b',LineWidth=1);
        hold on
        X=linspace(x1,x12,20); a=linspace(guess_So0oR(2),guess_So0oR(2),20);plot(X, a, 'g',LineWidth=1);
        hold on
        X=linspace(x12,x12,20); a=linspace(guess_So0oR(2), guess_So0oR(3),20); plot(X, a, ':','color', "#7e1e9c",LineWidth=2);
        hold on
        X=linspace(x12,xReS,20); a=linspace(guess_So0oR(3), guess_So0oR(3),20); plot(X, a, 'go',LineWidth=1.5);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_So0oR(3), guess_So0oR(4),20); plot(X, a, 'y*',LineWidth=1.5);
        hold on
        X=linspace(xReS,x21,20); a=linspace(guess_So0oR(4), guess_So0oR(4),20); plot(X, a, 'ko',LineWidth=1.5);
        hold on
        X=linspace(x21,x22,20); a=linspace(guess_So0oR(4), A22,20); plot(X, a, 'm',LineWidth=1.5);
        hold on
        X=linspace(x22,x2,20); a=linspace(A22, A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1)
        legend({'$A_1$','1-$S$','$A_{11}$','$lC$','$A_{12}$','$S_{=0}$','$A_{21}$','$rC$',...
                    '2-$R$', '$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
    
  elseif exitflag_SSoR== 1 && guess_SSoR(1)>0 && guess_SSoR(2)>0  
    
    A21=A22;
    ReS=(guess_SSoR(1)*guess_SSoR(3)-guess_SSoR(2)*guess_SSoR(4))/(guess_SSoR(1)-guess_SSoR(2));
    lS=(A1*u1-guess_SSoR(1)*guess_SSoR(3))/(A1-guess_SSoR(1));
    % distances 

    x1=lS*t;
    xReS = ReS*t;
    x22=(u22+c22)*t;
    x21 = x22;
    x2=(u2+c2)*t;
  
        % Plot A against x
        X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
        hold on
        X=linspace(x1,x1,20); a=linspace(A1,guess_SSoR(1),20);plot(X, a, 'b',LineWidth=1);
         hold on
        X=linspace(x1,xReS,20); a=linspace(guess_SSoR(1),guess_SSoR(1),20);plot(X, a, 'g',LineWidth=1);
        hold on
        X=linspace(xReS,xReS,20); a=linspace(guess_SSoR(1), guess_SSoR(2),20); plot(X, a, 'color', "#7e1e9c",LineWidth=1);
        hold on
        X=linspace(xReS,x21,20); a=linspace(guess_SSoR(2), guess_SSoR(2),20); plot(X, a, 'c',LineWidth=1);   
        hold on
        X=linspace(x21,x22,20); a=linspace(guess_SSoR(2),A21,20);plot(X, a, 'm:',LineWidth=1);
        hold on
        X=linspace(x22,x2,20); a=linspace(A22, A2,20); plot(X, a, 'r',LineWidth=1);
        hold on
        X=linspace(x2,x2+1,20); a=linspace(A2, A2,20); plot(X, a, 'color', "#6e750e",LineWidth=1);
        legend({'$A_1$','1-$S$','$A_{11}$','$S_{\neq 0}$','$A_{12}$',...
                '$C$','2-$R$', '$A_2$'},...
                'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
     
    elseif exitflag==1  && x(3)>0 && x(1)>0 
          c1=sqrt(m*A1.^m+n*A1.^(-n)); 
          c11=sqrt(m*x(1).^m+n*x(1).^(-n)); 
          c22=sqrt(k.*(m*(x(3)/A2bar).^m+n*(x(3)/A2bar).^(-n)));
          c2=sqrt(k.*(m*(A2/A2bar).^m+n*(A2/A2bar).^(-n)));
          lS=(A1*u1-x(1)*x(2))/(A1-x(1));
          rS=(A2*u2-x(3)*x(4))/(A2-x(3));
          xS1=lS*t; xS2=rS*t;
          x1=(u1-c1)*t;     x2=(u2+c2)*t;
          x11=(x(2)-c11)*t;  x22=(x(4)+c22)*t;

      if x(1) <= A1 && x(3) <= A2 %RoR

                X=linspace(x1-1,x1,20); a=linspace(A1,A1,20);plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(x1,x11,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(x11,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20); plot(X,a,'m:',LineWidth=1);
                hold on
                X=linspace(0,x22,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(x22,x2,20); a=linspace(x(3),A2,20);plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(x2,x2+1,20); a=linspace(A2,A2,20);  plot(X, a, 'color', "#6e750e",LineWidth=1);
                 %title('$RoR$','interpreter','latex')
                legend({'$A_1$','1-$R$','$A_{11}$', '$C$','$A_{22}$','2-$R$','$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

     elseif x(1) <= A1 && x(3) > A2 %RoS

                X=linspace(x1-1,x1,20); a=linspace(A1,A1,20); plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(x1,x11,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(x11,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20); plot(X, a,'m:',LineWidth=1);
                hold on
                X=linspace(0,xS2,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(xS2,xS2,20); a=linspace(x(3),A2,20);plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(xS2,x2+1,20); a=linspace(A2,A2,20);plot(X, a,'color', "#6e750e",LineWidth=1);
                 %title('$RoS$','interpreter','latex')
                 legend({'$A_1$','1-$R$','$A_{11}$', '$C$','$A_{22}$','2-$S$','$A_2$'},...
                     'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');

      elseif x(1) > A1 && x(3) <= A2  %SoR

                X=linspace(x1-1,xS1,20); a=linspace(A1,A1,20);plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(xS1,xS1,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(xS1,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20);plot(X,a,'m:',LineWidth=1);
                hold on
                X=linspace(0,x22,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(x22,x2,20); a=linspace(x(3),A2,20); plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(x2,x2+1,20); a=linspace(A2,A2,20);  plot(X, a, 'color', "#6e750e",LineWidth=1);
                %title('$SoR$','interpreter','latex')
                legend({'$A_1$','1-$S$','$A_{11}$', '$C$','$A_{22}$','2-$R$','$A_2$'},...
                    'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
                %xlim([-1.2,2.6])

       elseif x(1) > A1 && x(3) > A2  %SoS

                X=linspace(x1-1,xS1,20); a=linspace(A1,A1,20);plot(X, a, 'color', "#f97306",LineWidth=1);
                hold on
                X=linspace(xS1,xS1,20); a=linspace(A1,x(1),20);plot(X, a, 'b',LineWidth=1);
                hold on
                X=linspace(xS1,0,20); a=linspace(x(1),x(1),20);plot(X, a, 'g',LineWidth=1);
                hold on
                X=linspace(0,0,20); a=linspace(x(1),x(3),20);plot(X,a,'m:',LineWidth=1);
                hold on
                X=linspace(0,xS2,20); a=linspace(x(3),x(3),20);plot(X, a, 'c',LineWidth=1);
                hold on
                X=linspace(xS2,xS2,20); a=linspace(x(3),A2,20);plot(X, a, 'r',LineWidth=1);
                hold on
                X=linspace(xS2,x2+1,20); a=linspace(A2,A2,20);plot(X, a, 'color', "#6e750e",LineWidth=1);
                 %title('$SoS$','interpreter','latex')
                 legend({'$A_1$','1-$S$','$A_{11}$', '$C$','$A_{22}$','2-$S$','$A_2$'},...
                     'interpreter','latex','FontSize',12,'Location','northeastoutside','Orientation','vertical');
        
       end
   end
   xlabel('$x$','interpreter','latex','FontSize',12)
    ylabel('$A$','interpreter','latex','FontSize',12)
   ylim([0.2,1])
   %xlim([-1.2, 2.8])
   xlim([x1-1, x2+1])
   
 end
 % Define functions to find_A11 and lRout
     function u=lRout(A)
          f1=@(a) (m*a.^(m-2)+n*a.^(-(n+2))).^(1/2);
          int1=quad(f1,A,A1);
          u=u1+int1;
        end
    
    function v=find_A11(A)
    c11=sqrt(m*(A)^m+n*(A)^(-n));
    u11=lRout(A);
    v=u11-c11;
    end

    % Define functions to find_A22 and rRout
    function u=rRout(A)
       f1=@(a) (k*(m*a.^(m-2)/A2bar.^m+n*a.^(-n-2)/A2bar.^(-n))).^(1/2);
      int1=integral(f1,A2,A);
    u=u2+int1;
    end

    function v=find_A22(A)
    c22=sqrt(k*(m*(A/A2bar).^m+n*(A/A2bar).^(-n)));
    u22=rRout(A);
    v=u22+c22;
    end
end