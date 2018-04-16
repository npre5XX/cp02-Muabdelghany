% P1 Mark's BC
clc;
clear;
close all;
seg_t=1.5;seg_s=0.9;Q=100;
L=5;         % thickness of the whole slab
k=101;       % number of points
del=L/(k-1);
a=-1/(3*del*seg_t);
b=(((seg_t-seg_s)*del)-(2*a));
c=a;
s=Q*del;

for i=1:k
    if i==1 | i==k
        S(i,1)=0;
    else
    S(i,1)=s;
    end
end

    A(1,1)=1;
    A(k,k)=1;
    for j=2:k
        A(1,j)=0;
    end
    for j=1:k-1
        A(k,j)=0;
    end
    
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

D=det(A);
P=inv(A);
phi=P*S

x=linspace(0,L,k);
plot(x,phi(:))
xlabel('x [Cm]')
ylabel('\phi(x)')
title('P1: Flux Distribution Inside the Slab Evaluated Using Mark''s BC')