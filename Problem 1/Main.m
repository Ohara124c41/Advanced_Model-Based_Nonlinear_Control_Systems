clc, clear all, close all

%% Define the model parameters

L0 = 0.09; % H
R = 60; % Ohms
m = 10e-2; % kg
g = 9.81; % m/s^2

tf = 30; % simulation time (in seconds)


%% Define ref values
yr = 0.01;


%% Part 2) Constant input voltage

% This part is defined before part 1, because the reference values are used
% in the simulation

Ir = sqrt(2*g*m*yr^2 /(L0*yr));
er = R*Ir;

I0 = Ir;
v0 = 0;
y0 = yr;
u0 = er;

x10 = y0;
x20 = v0;
x30 = I0;

%% Part 1) Run Simulink model

sim('sim_model.slx');

%% Plot results

I = vars.Data(:,1);
y = vars.Data(:, 2);
v = vars.Data(:, 3);
e = vars.Data(:, 4);
t = vars.Data(:, 5);

figure(1)
% subplot(3,1,1);
subplot(3,2,1);
% plotyy(t, y, t, e), grid on
plot(t,y,'linewidth',2), grid on
xlabel('Time (s)')
ylabel('Y-position (m)');
title('Y-position vs. Time');
legend('Y-position', 'Input');

% subplot(3,1,2);
subplot(3,2,3);
% plotyy(t, v, t, e), grid on
plot(t,v,'linewidth',2), grid on
xlabel('Time (s)')
ylabel('Velocity (m/s)');
title('Velocity vs. Time');
legend('Velocity', 'Input');

% subplot(3,1,3);
subplot(3,2,5);
% plotyy(t, I, t, e), grid on
plot(t,I,'linewidth',2), grid on
xlabel('Time (s)')
ylabel('Current (A)');
title('Current vs. Time');
legend('Current', 'Input');

%% Part 3) Construct the Jacobian
K = L0*x10/(2*m);
K1 = 2*K*x30^2 /x10^3;
K2 = -2*K*x30/x10^2;

F30 = u0/L0 - R*x30/L0 + x30*x20/x10;
K3 = u0/(L0*x10) - R*x30/(L0*x10) - x30*x20/x10^2;
K4 = x30/x10;
K5 = -R/L0 + x20/x10;
K6 = 1/L0;

%% Linearized model
sim('sim_model_linearized.slx');

% Plot linearized results

x1 = vars.Data(:,1);
x2 = vars.Data(:, 2);
x3 = vars.Data(:, 3);
u = vars.Data(:, 4);
t = vars.Data(:, 5);

% figure
% subplot(3,1,1);
subplot(3,2,2);
% plotyy(t, x1, t, u), grid on
plot(t,x1,'linewidth',2), grid on
xlabel('Time (s)')
ylabel('Y-position (m)');
title('Y-position vs. Time');
legend('Y-position', 'Input');

% subplot(3,1,2);
subplot(3,2,4);
% plotyy(t, x2, t, u), grid on
plot(t,x2,'linewidth',2), grid on
xlabel('Time (s)')
ylabel('Velocity (m/s)');
title('Velocity vs. Time');
legend('Velocity', 'Input');

% subplot(3,1,3);
subplot(3,2,6);
% plotyy(t, x3, t, u), grid on
plot(t,x3,'linewidth',2), grid on
xlabel('Time (s)')
ylabel('Current (A)');
title('Current vs. Time');
legend('Current', 'Input');



% J = [0 1 0;K1 0 -K2;I0*v0/y0^2, I0/y0, (v0/y0 - R/(L0*y0))];
% J = [0, 1, 0;K1, 0, -K2; 0, I0/y0, 0];
fprintf("The Jacobian of the system is:\n");
J = [0, 1, 0;K1, 0, K2;K3, K4, K5];
disp(J);
% The Jacobian has a positive eigenvalue
% which means that the equilibrium point is unstable


% Define the model as a state-space
A = J;
B = [0;0;K6];
C1 = [1, 0, 0]; % x1 as output
C2 = [0, 1, 0]; % x2 as output
C3 = [0, 0, 1]; % x3 as output
D = 0;


%% Part 4 and 6: LQ control

% Develop an LQ full-state feedback control

% Define the penalization matrix

C = [1 0 0];

% Q = 50*C'*C;
Q = [10 0 0;0 0.01 0;0 0 0.01];
Rc = 100;

Klq = lqr(A,B,Q,Rc);
A_lq = A-B*Klq;

%% Apply the controller
sys = ss(A,B,C,D);
sysLQ = ss(A_lq, B, C, D);

%% Part 5) SDRE METHOD

u = 1; % No input
y = zeros(3,1);
Q=[100 0 0;0 1 0;0 0 1];
p=20;
Rc = 0.1;

tspan = 0:0.01:60;
N = length(tspan);

clear x1 x2 x3

x1 = x10;
x2 = x20;
x3 = x30;
% for i = 2:N
%     [G,M, eig]= lqr(A,B,Q,Rc);
%     x = [x1(i-1);x2(i-1);x3(i-1)];
%     dx = A*x+B*u;
%     x1(i) = x1(i-1)+dx(1);
%     x2(i) = x2(i-1)+dx(2);
%     x3(i) = x3(i-1)+dx(3);
%     
%     y(i) = C*x;
%     
%     uu=(B'*M*[x1(i-1);x2(i-1);x3(i-1)]);
%     u=-abs(uu)^(1/((2*p)-1))*sign(uu);
% 
%     Rc=abs(u)^(2*p-2);
% end
P = ARE_diag(A,B,Q,0.01)
Ksdre = inv(R)*B'*P;

% Construct new A
Asdre = A - B*Ksdre;
ssSDRE = ss(Asdre, B, C, D);
figure
step(ssSDRE), grid on
title('SDRE Method')

%rlocus(ssSDRE), title('Root Loci; State-Dependent Riccati Equation')
%bode(ssSDRE), title('Bode; State-Dependent Riccati Equation')
%nyquist(ssSDRE), title('Nyquist; State-Dependent Riccati Equation')

%% Part 6) State-Space model with LQ
x10 = 0.03;
sim('sim_model_ss_LQ.slx');

% Plot linearized results

x1 = vars.Data(:,1);
u = vars.Data(:, 2);
t = vars.Data(:, 3);

figure(3)
% plotyy(t, x1, t, u), grid on
plot(t, x1, 'linewidth', 2), grid on
xlabel('Time (s)');
ylabel('Position (m)');
title('Position vs. Time');

%% Part 7)

% if the input is 10v, lets calculate the current and displacement

u8 = 10;
I8 = u8/R;
y8 = fsolve(@(x)(I8 - sqrt(2*g*m*x^2/(L0*x))),10);
fprintf("For an input of 10v, the maximum displacement is %.4f m\n", y8);

%% Desired value ( Question 2)
Ir = sqrt(2*g*m*yr^2 /(L0*yr));
er = R*Ir;

%% part 8)

% Theoretical Question

%% Part 9) Full-order observer

obpoles = [-2 -2 -2];
l = acker(A', C1', obpoles);
Lc = l';

Nbar = rscale(sys,Klq);

At = [A-B*Klq, B*Klq;zeros(size(A)), A-Lc*C];

Bt = [B*Nbar; zeros(size(B))];

Ct = [C, zeros(size(C))];


sys2 = ss(At,Bt,Ct,0);
%nyquist(sys2), title('Nyquist')
%locus(sys2), title('Root Loci; Full-Order Observer')
%bode(sys2), title('Bode; Full-Order Observer')
figure(4)
step(sysLQ), hold on
step(sys2); 
legend('Full-State feedback', 'F-S feedback with observer');
grid on
title('Full-order Observer response');


%% Part 10) Reduced order observer
% Assuming only the y-position is measured

% dx1 = A11*x1 + A12*x2 + A13*x3

C = [1, 0, 0];
Q= 0.01*eye(3);
Q(1,1)=10;
Qobs = 1000;
Rc = 100;

Klq = lqr(A,B,Q,Rc);
Alq = A-B*Klq;

% plant states
x1(1) = x10;
x2(1) = x20;
x3(1) = x30;

% Observer states
x4(1) = x10;
x5(1) = x20;
x6(1) = x30;

u = u0; % No input
y = zeros(3,1);

sys = ss(Alq,B,C,D);


[kest, L, P] = kalman(sys, Qobs, Rc, []);
Aobs = Alq;
Bobs = B;
Cobs = C;

Anew = [Alq, zeros(3);
        L*C, Aobs-L*Cobs];


Bnew = [B, B;zeros(3,1), Bobs];
Cnew = [C, zeros(1,3);
        zeros(3,3), eye(3)];


for i = 2:N
    x = [x1(i-1);x2(i-1);x3(i-1);x4(i-1);x5(i-1);x6(i-1)];
    dx = Anew*x+Bnew*u;
    
    % Plant
    x1(i) = x1(i-1)+dx(1);
    x2(i) = x2(i-1)+dx(2);
    x3(i) = x3(i-1)+dx(3);
    
    % Observer
    x4(i) = x4(i-1)+dx(4);
    x5(i) = x5(i-1)+dx(5);
    x6(i) = x6(i-1)+dx(6);
    
    
%     y(i) = Cnew*x;
end

figure
subplot(3,1,1)
plot(tspan, x1, 'b', 'linewidth', 2), hold on
plot(tspan, x4, 'r--','linewidth',2)
grid on
legend('Plant', 'Observer')
xlabel('Time (s)')
ylabel('Position (m)');

subplot(3,1,2)
plot(tspan, x2, 'b', 'linewidth', 2), hold on
plot(tspan, x5, 'r--', 'linewidth', 2)
grid on
legend('Plant', 'Observer')
xlabel('Time (s)')
ylabel('Velocity (m/s)');

subplot(3,1,3)
plot(tspan, x3, 'b', 'linewidth', 2), hold on
plot(tspan, x6, 'r--', 'linewidth', 2)
grid on
legend('Plant', 'Observer')
xlabel('Time (s)')
ylabel('Current (A)');