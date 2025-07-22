%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   A bifurcation code for the standard non-dimensional blood flow model Riemann problem using general tube law.
%
%   Bethel fecit, AD 2025.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Vk_Bifurcation
% First, clear the workspace and close any open figures
clear all
clc
close all

% Set precision
format short

% Figure
figure 
box on

% Set initial conditions for the Riemann problem
u1 = 0.2;
A1 = 0.8;
u2 = 0.4;
A2bar= 1;
k=1.1;
m=10; % Go to classical_Vk_gen.m and choose the appropriate initial guess for the classical case
n=0;
dur2=0.005;
A2range= 0.06:dur2:2;

%initial guess for the resonant solutions
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

% Define the loop for A2
for  An = 1:length(A2range)
        A2 = A2range(An);

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
 
    % plot1 A_11
    figure(1)
    box on
    subplot 211

    if exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            hold on
            plot(A2,A11, '.',"MarkerEdgeColor", "#7e1e9c", "MarkerFaceColor","#7e1e9c"); % Plots purple dots on bifurcation graph
       
        elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
            hold on
            plot(A2,A11, 'y.'); % Plots yellow dots on bifurcation graph

        elseif exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0 
            A12=A11;
            hold on
            plot(A2,A12, '.',"MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor",[0.8500 0.3250 0.0980]); % Plots orange dots on bifurcation graph
        
         elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0 
            A12=A11;
            hold on
            plot(A2,A12, 'm.'); % magenta dots on bifurcation graph

         elseif exitflag_So0oR == 1 && guess_So0oR(1)>0 && guess_So0oR(2)>0 && guess_So0oR(3)>0 && guess_So0oR(4)>0
            hold on
            plot(A2,guess_So0oR(2),'.',"MarkerEdgeColor",  "#6e750e", "MarkerFaceColor","#6e750e" ); % Plots olive dots on bifurcation graph
    
         elseif exitflag_SSoR == 1 && guess_SSoR(1)>0 && guess_SSoR(2)>0 
            hold on
            plot(A2,guess_SSoR(1), '.',"MarkerEdgeColor", "#cea2fd", "MarkerFaceColor","#cea2fd"); % Plots lilac dots on bifurcation graph
       
        elseif exitflag==1 && x(3)>0 && x(1)>0
        
        if x(1) <= A1 && x(3) <= A2 % RoR
            hold on
            plot(A2,x(1), 'b.');    % Plots blue dots on  bifurcation graph
    
        elseif x(1) <= A1 && x(3) > A2 % RoS
            hold on
            plot(A2,x(1), 'g.');    % Plots green dots on  bifurcation graph
    
        elseif x(1) > A1 && x(3) <= A2 % SoR
            hold on
            plot(A2,x(1), 'c.');   % Plots cyan dots on  bifurcation graph
    
        elseif x(1) > A1 && x(3) > A2 % SoS
            hold on
            plot(A2,x(1), 'r.');    % Plots red dots on  bifurcation graph
        end
   
     else                         % No solution
            hold on
            plot(A2,u2, 'k.');   % Plots black dots on  bifurcation graph
     end  
 
    xlabel('$A_2$','interpreter','latex','FontSize',11)
    ylabel('$A_{11}$','interpreter','latex','FontSize',11)
     ylim([0 2])

% plot2 A_22

    figure(1)% ('Position', [360, 98, 560, 177]) % Set the pixel size of the figure
    box on
    subplot 212
  if exitflag_Ro0oS == 1 && guess_Ro0oS(1)>0 && guess_Ro0oS(2)>0 && guess_Ro0oS(3)>0 && guess_Ro0oS(4)>0
            hold on
            plot(A2,guess_Ro0oS(4), '.',"MarkerEdgeColor", "#7e1e9c", "MarkerFaceColor","#7e1e9c"); % Plots purple dots on bifurcation graph
       
        elseif exitflag_Ro0oR == 1 && guess_Ro0oR(1)>0 && guess_Ro0oR(2)>0 && guess_Ro0oR(3)>0 && guess_Ro0oR(4)>0
            hold on
            plot(A2,guess_Ro0oR(4), 'y.'); % Plots yellow dots on bifurcation graph

        elseif exitflag_RoSS == 1 && guess_RoSS(1)>0 && guess_RoSS(2)>0 
             hold on
            plot(A2,guess_RoSS(2), '.',"MarkerEdgeColor", [0.8500 0.3250 0.0980], "MarkerFaceColor",[0.8500 0.3250 0.0980]); % Plots orange dots on bifurcation graph
        
        elseif exitflag_RoSR == 1 && guess_RoSR(1)>0 && guess_RoSR(2)>0 
             hold on
            plot(A2,guess_RoSR(2), 'm.'); % Pot magenta dots on  bifurcation graph

        elseif exitflag_So0oR == 1 && guess_So0oR(1)>0 && guess_So0oR(2)>0 && guess_So0oR(3)>0 && guess_So0oR(4)>0
            hold on
            plot(A2,A22,'.',"MarkerEdgeColor",  "#6e750e", "MarkerFaceColor","#6e750e" ); % Plots olive dots on bifurcation graph
    
         elseif exitflag_SSoR == 1 && guess_SSoR(1)>0 && guess_SSoR(2)>0 
            A21=A22;
            hold on
            plot(A2,A21, '.',"MarkerEdgeColor", "#cea2fd", "MarkerFaceColor","#cea2fd"); % Plots lilac dots on bifurcation graph
       elseif exitflag==1 && x(3)>0 && x(1)>0
        
            if x(1) <= A1 && x(3) <= A2
                hold on
                plot(A2,x(3), 'b.');    % Plots blue dots on  bifurcation graph
        
            elseif x(1) <= A1 && x(3) > A2
                hold on
                plot(A2,x(3), 'g.');    % Plots green dots on  bifurcation graph
        
            elseif x(1) > A1 && x(3) <= A2
                hold on
                plot(A2,x(3), 'c.');   % Plots cyan dots on  bifurcation graph
        
            elseif x(1) > A1 && x(3) > A2
                hold on
                plot(A2,x(3), 'r.');    % Plots red dots on  bifurcation graph
            end
        else                        % No solution
            hold on
            plot(A2,u2, 'k.');    % Plots black dots on  bifurcation graph
   end      
 
    xlabel('$A_2$','interpreter','latex','FontSize',11)
    ylabel('$A_{22}$','interpreter','latex','FontSize',11)
    ylim([0 2])

   end
 end %A2
 
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
    