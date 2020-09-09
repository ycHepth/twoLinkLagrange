% test for 2-link model
% method : lagrange methods

clear all 
clc

syms q1 alpha1 l1 lc1 m1 Ixx1 Iyy1 Ixy1 Iyz1 Ixz1 Izz1 real
syms q2 alpha2 l2 lc2 m2 Ixx2 Iyy2 Ixy2 Iyz2 Ixz2 Izz2 real
syms dq1 dq2 real
syms ddq1 ddq2 real

% SDH parameters
a1 = l1;
a2 = l2;
alpha1 = 0;
alpha2 = 0;
d1 = 0;
d2 = 0;

a = [a1,a2]';
alpha = [alpha1,alpha2]';
q = [q1,q2]';
dq = [dq1,dq2]';
d = [d1,d2]';

% calculate neighbor corrinate transformed matrix
n = 2;
for k =1:n
    A(:,:,k) = SDHmatrix(a(k),alpha(k),q(k),d(k));
end

A00 = eye(4);
for i = 2:n
    T(:,:,1) = A(:,:,1); % T10
    T(:,:,i) = simplify(T(:,:,i-1)*A(:,:,i));
end

% inertia matrix of link 1
I1 = [Ixx1,Ixy1,Ixz1;
      Ixy1,Iyy1,Iyz1;
      Ixz1,Iyz1,Izz1];
% inertia matrix of link 2
I2 = [Ixx2,Ixy2,Ixz2;
      Ixy2,Iyy2,Iyz2;
      Ixz2,Iyz2,Izz2];

  
In = cat(3,I1,I2);

z0 = [0,0,1]';
for k = 1:n
    Z(:,:,k) = T(1:3,3,k);
end



  




