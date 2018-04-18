% P3 Marshak's BC
clc;
clear;
close all;
seg_t=1.5;seg_s=0.9;Q=100;
seg_0=seg_t-seg_s;
L=5;                        % thickness of the whole slab
k=101;                      % number of mesh points
del=L/(k-1);

% F0 and F1 equations' coefficients
A=3*seg_t*(del)^2;
B=(7*seg_t*(del)^2)/9;
C=2/(3*del*seg_t);

% (i) is the space mesh index, and (n) is the iterations counter
% Arbitrary intial guess of F0 and F1 of the 1st iteration
for i=1:k
    F0(i,1)=5;
    F1(i,1)=10;
end

n=0;
mx_dif=1;                    % intiating the convergance criterion
while mx_dif>0.00001         % convergance criterion for the flux
    
    n=n+1;                   % next iteration
        
% +++++ Finding (F0) and (F1) for each mesh point for iteration (n) ++++++
        for i=1:k-2

        F0(i+1,n+1)=0.5*((F0(i+2,n)+F0(i,n))-A*(seg_0*F0(i,n)-2*seg_0*F1(i,n)-Q));
        F1(i+1,n+1)=0.5*((F1(i+2,n)+F1(i,n))-B*(((4*seg_0+5*seg_t)*F1(i,n))-(2*(seg_0*F0(i,n+1)+Q))));

        end
        
% Choose the two boundary values of F1 for each iteration to be equal the neighborhood values this choice
% should be done by different approaches, but I found this most convenient
        F1(k,n+1)=F1(k-1,n+1);
        F1(1,n+1)=F1(2,n+1);
        
        % implementing Marshak's BCs
        F0(1,n+1)=(1/(1+C))*((3/4)*F1(1,n+1)+C*F0(2,n+1));
        F0(k,n+1)=(1/(1+C))*((3/4)*F1(k,n+1)+C*F0(k-1,n+1));
        
        % find phi for each mesh point for this iterantion
        for i=1:k
            phi(i,n)=F0(i,n+1)-2*F1(i,n+1);

            % Specifying the difference in phi between each two successive iteration
            
            if n>1    % skip evaluating the difference for the 1st iteration       
            dphi(i,n-1)=phi(i,n)-phi(i,n-1);           % difference in phi
            else
            dphi=1;               % put the difference=1 for 1st iteration
            end
        
        end
        % find the maximum difference in flux between all the mesh points
        mx_dif=max(dphi);    
end
n               % display the number of iterations required for convergance

% Plot (phi) vs the width of the slab (x)
x=linspace(0,L,k);
plot(x,phi(:,n))
xlabel('x [Cm]')
ylabel('\phi(x)')
title('P3: Flux Distribution Inside the Slab Evaluated Using Marshak''s BC')
grid on
hold on