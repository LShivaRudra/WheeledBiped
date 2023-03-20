thetaHip=2*pi/3;
r=0.095; % Wheel Radius;
M=5; % Mass of the bot;
L=sqrt(0.0784*sin(thetaHip/2)*sin(thetaHip/2) + (3.969*0.001));
mb=0.1; % Mass of the wheel;
Ib=M*L*L; % Moment of inertia of body abt y-axis;
Iw= 4.5125*0.0001;% Moment of inertia of the wheel;
b=0.75; % Friction coeff b/w the wheel and the rod;
g=9.82; % Acc. due to gravity;
D=0.3; % distance b/w centers of both the wheels;

It=20; % Moment of inertia of the bot about the vertical axis;

delta = (r*r)*(2*M*M*L*L + 2*mb*M*L*L - 2*mb*Ib - M*Ib) + 2*M*L*L*Iw - 2*Iw*Ib;
A22 = 2*(b*M*L*L - b*Ib + M*L*r*b)/delta;
A23 = M*M*g*L*L*r*r/delta;
A42 = (4*b*r*mb + 2*M*b*r + 4*b*Iw/r - 2*M*b*L)/delta;
A43 = (2*M*mb*g*L*r*r + M*M*g*L*r*r + 2*Iw*M*g*L)/delta;
B2 = (-M*L*L*r + Ib*r - M*L*r*r)/delta;
B4 = (-2*mb*r*r - M*r*r - 2*Iw + M*L*r)/delta;
B6 = (D*r)/(D*D*Iw + mb*D*D*r*r + 2*It*r*r);

% TL = (Tphi + Tomega)/2;
% TR = (Tphi + Tomega)/2;
% 
% T = [TL;TR];

% s = [xb;xb_dot;phi;phi_dot;alpha;alpha_dot];
% s_dot = [xb_dot;xb_ddot;phi_dot;phi_ddot;alpha_dot;alpha_ddot];

A = [0 1 0 0 0 0; 0 A22 A23 0 0 0; 0 0 0 1 0 0;
     0 A42 A43 0 0 0; 0 0 0 0 0 1; 0 0 0 0 0 0];
B = [0 0; B2 B2; 0 0; B4 B4; 0 0; B6 -B6];
C = eye(6);
D = zeros(6,2);

Q=eye(6);
R=0.8*eye(2)+0.2*ones(2,2);
[K,S,P]=lqr(A,B,Q,R);

save('LQR_vars_WB.mat','A','B','K','S','P','Q','R');
clear all;
load('LQR_vars_WB.mat');
