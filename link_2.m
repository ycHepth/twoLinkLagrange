% test for 2-link model
% method : lagrange methods

clear all 
clc
tic;

syms q1 alpha1 l1 lc1 m1 Ixx1 Iyy1 Ixy1 Iyz1 Ixz1 Izz1 real
syms q2 alpha2 l2 lc2 m2 Ixx2 Iyy2 Ixy2 Iyz2 Ixz2 Izz2 real
syms dq1 dq2 real
syms ddq1 ddq2 real

% ======================================================
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
ddq = [ddq1,ddq2]';
d = [d1,d2]';

% ======================================================
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
% why need tensor here?

z0 = [0,0,1]';
for k = 1:n
    Z(:,:,k) = T(1:3,3,k);
end

% ======================================================
o0 = [0,0,0]';
for k = 1:n
    o(:,:,k) = T(1:3,4,k);
end

% calculate mess center position in gravity coordinate
oc1x = lc1 - l1;
oc1y = 0;
oc1z = 0;
Pc1 = [oc1x,oc1y,oc1z,1]';

oc2x = lc2 - l2;
oc2y = 0;
oc2z = 0;
Pc2 = [oc2x,oc2y,oc2z,1]';

% linkage mess center coordinate matrix
Pc = [Pc1,Pc2];

for k = 1:n
    oc(:,:,k) = simplify(T(:,:,k)*Pc(1:4,k));
end

% ======================================================
% calculate jacobian matrix
Jvc1 = simplify([cross(z0,(oc(1:3,1,1)-o0)),zeros(3,1)]);
Jw1  = [z0,zeros(3,1)];

Jvc2 = simplify([cross(z0,(oc(1:3,1,2)-o0)),cross(Z(:,:,1),(oc(1:3,1,2)-o(:,:,1)))]);
Jw2  = [z0,Z(:,:,1)];

% Jacobina matrix
Jvc = cat(3,Jvc1,Jvc2);
Jwc = cat(3,Jw1,Jw2);

% ======================================================
%            calculate inertia matrix D(q)
% ======================================================
m = [m1,m2]';
D0 = zeros(2);
D = D0;
for i = 1:n
    D = D + simplify([m(i).*Jvc(:,:,1).'*Jvc(:,:,i) + Jwc(:,:,i).'*T(1:3,1:3,i) * In(:,:,i)*T(1:3,1:3,i).'*Jwc(:,:,i)]);
end

disp('inertia matrix term D: ');
D
% ======================================================
%            calculate Cq
% ======================================================

Cq = TwoLinkDirectC(2,D,q,dq);
disp('Cq: ');
Cq

% ======================================================
%            calculate Gravity term G
% ======================================================
syms g real
g = [0,-9.806,0]';
P = 0;
for k = 1:n
    P = P + m(k)*g.'*oc(1:3,:,k);
end

for i = 1:n
    gk(i) = diff(P,q(i));
end

gk = gk';
disp(' Gravity term G:');
gk



% ======================================================
%            calculate tau
% ======================================================
% dynamics: D*ddq + Cq*dq + gk = tau

tau = simplify(D*ddq+Cq*dq + gk);
disp('torque output tau:');
tau

% ======================================================
%            numerical calculation
% ======================================================
dats={alpha1,l1,lc1,m1,Ixx1,Iyy1,Izz1,...
      alpha2,l2,lc2,m2,Ixx2,Iyy2,Izz2};
datn={0,1,0.5,1,0,0,1/3,...
      0,1,0.5,1,0,0,1/3};

tau = vpa(subs(tau,dats,datn),6);

dats_angle = {q1,dq1,ddq1,...
              q2,dq2,ddq2};
          
datn_angle = {pi/6,0.5,0.5,...
              pi/3,0.25,0.45};
          
tau = vpa(subs(tau,dats_angle,datn_angle),6)
toc;


% ======================================================
%            dynamics term in numerical
% ======================================================
disp('======================');
disp('Numerical term');
D = vpa(subs(D,dats,datn),6)
C = vpa(subs(Cq,dats,datn),6)
G = vpa(subs(gk,dats,datn),6)



