% P1 Marshak's BC
clc;
clear;
close all;
seg_t=1.5;seg_s=0.9;Q=100;
L=5;                        % thickness of the whole slab
k=101;                      % number of mesh points
del=L/(k-1);

% Matrices Coefficents
a=-1/(3*del*seg_t);
b=(((seg_t-seg_s)*del)-(2*a));
c=a;
s=Q*del;

% i,j are indices which are indication on the number of mesh points
% Build the S Matrix
for i=1:k
    if i==1 | i==k
        S(i,1)=0;           % Because of the Marshak's BCs
    else
    S(i,1)=s;
    end
end

% Specify the 1st and kth rows of the A Matrix from Marshak's BCs
    A(1,1)=1-2*a;
    A(1,2)=2*a;
    A(k,k-1)=2*a;
    A(k,k)=1-2*a;  
    for j=3:k
        A(1,j)=0;
    end
    for j=1:k-2
        A(k,j)=0;
    end
    
% Build the Coefficient Matrix A
for i=2:k-1
        
    for j=1:k
        
        if j==i-1
            A(i,j)=a;
        elseif j==i
            A(i,j)=b;
        elseif j==i+1
            A(i,j)=c;
        else
            A(i,j)=0;
        end
        
     end
end

% Find the Inverse of (A) and computing phi
D=det(A);
P=inv(A);
phi=P*S

% Plot (phi) vs the width of the slab (x)
x=linspace(0,L,k);
plot(x,phi(:))
    xlabel('x [Cm]')
    ylabel('\phi(x)')
    title('P1: Flux Distribution Inside the Slab Evaluated Using Marshak''s BC')
    grid on