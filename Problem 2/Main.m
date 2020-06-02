clc, clear all, close all

M = 0.15; % mass of a beam
L = 0.5;
d = 0.75*0.254; % diameter of ball, in meters

% To calculate the mass of the ball, we use the density of the steel and
% its volume

%% Part 1
Vball = (4/3)*pi*(d/2)^3;

p_steel = 8050; % kg/m3;

m = p_steel*Vball

K = 0.0275; % torque constant
R = 5; % motor resistance
6
I = M*L^2 /12; % moment of inertia of the beam,

g = 9.81;

%% Define equilibrium point

% we define a reference point, which is the position of the ball in the
% beam. We define it at the center of the beam.

x10 = 0;
% x20 = -0.5;
x20 = 0;
x30 = 0;
x40 = 0;


% u0 = (R/K)*m*g*x20*cos(x10)/(I+m*x20^2);
u0 = 0;
%% Part 1) Run simulink model

sim('BallOnBeam_model.slx');
r = vars.Data(:,1);
theta = vars.Data(:,2);
T = vars.Data(:,3); % torque
u = vars.Data(:,4); % input
t = vars.Data(:,5);

% Plot results
figure
subplot(3,1,1)
plot(t, r, 'linewidth', 2), grid on
xlabel('Time(s)');
ylabel('Ball Position (m)');

subplot(3,1,2)
plot(t, theta, 'linewidth', 2), grid on
xlabel('Time(s)');
ylabel('Angular Position of beam (rad)');

subplot(3,1,3)
plot(t, T, 'linewidth',2), grid on
xlabel('Time(s)');
ylabel('Motor Torque (Nm)');

%% Linearized system and state-space representation

F0 = (K/R)*(u0 - K*x30)/(I+m*x20^2) - m*g*x20*cos(x10)/(I+m*x20^2) - 2*m*x20*x30*x40/(I+m*x20^2);

K1 = (K/R)/(I+m*x20^2);
K2 = m*g*x20*sin(x10)/(I+m*x20^2);

K3 = -(K/R)*(u0-K*x30)*2*m*x20/(I+m*x20^2)^2;
K4 = -(m*g*cos(x10) - m*g*x20*cos(x10)*2*m*x20/(I+m*x20^2)^2);
K5 = -(2*m*x40*x30*(I+m*x20^2) - 2*m*x20*x40*x30*2*m*x20)/(I+m*x20^2)^2;

K6 = (-K^2 /R)/(I+m*x20^2) - 2*m*x20*x40/(I+m*x20^2);

K7 = -2*m*x20*x30/(I+m*x20^2);

K8 = -5*g/7 *cos(x10);
K9 = x30^2;
K10 = 2*x20*x30;

A = [0 0 1 0;
     0 0 0 1;
     K2, (K3+K4+K5), K6, K7;
     K8, K9, K10, 0];
 
B = [0;0;K1;0];

C = [0, 1, 0, 0]; % r-position and theta as output
D = 0;



sys = ss(A,B,C,D);

% step(sys);


%% Part 2) Controller

% let's define the parameters
ts = 0.2; % settling time
mp = 0.05; % % overshoot

% Find seda and wn
[x, fval, exitflag] = fsolve(@(x)[4.6/(x(1)*x(2)) - ts;exp(-pi*x(1)/sqrt(1-x(1)^2)) - mp], [1,1]) % find the values of seda and wn
seda = x(1);
wn = x(2);
% Calculate poles
p1 = -seda*wn + 1i*wn*sqrt(1-seda^2); 
p2 = -seda*wn - 1i*wn*sqrt(1-seda^2);

% State-Feedback gain
Kc = real(place(A,B,[p1,p2,-20,-40]));
% closed-loop A
Ac = A-B*Kc;

%sys2 = ss(Ac, B, C, D);

sim('BallOnBeam_Model_Controller.slx');

theta = vars.Data(:,1);
y = vars.Data(:,2);
u = vars.Data(:,3);
t = vars.Data(:,4);

figure
% step(sys2)
plotyy(t,y,t,u), grid on


%% Part 3) Full-order Observer Kalman filter method
Q= 1e4;

Rc = .01;
[kest, L, P] = kalman(sys, Q, Rc, []);
Aobs = Ac;
Bobs = B;
Cobs = C;

Anew = [Ac, zeros(4);
        L*C, Aobs-L*Cobs];
   
    
Bnew = [B, B;zeros(4,1), Bobs];
Cnew = [C, zeros(1,4);
        zeros(4,4), eye(4)];
    
sysNew = ss(Anew, Bnew, Cnew, 0);

%nyquist(sysNew)
figure(10)
t = 0:1e-3:10;
u1 = 2*cos(5*2*pi*t')+2*cos(1*2*pi*t');+2*cos(.2*2*pi*t');; % disturbance
u2 = zeros(length(t),1); % control input 
y = lsim(sysNew,[u1 u2],t);
plot(t,y(:,1),t,y(:,2)) % first state of the observer is the position of the mass
legend('plant output','observer output estimate')

%% Part 4) Reduced-Order observer

 MatrixAaa = [0];
 MatrixAab = [1 0];
 MatrixAba = [0;-10];
 MatrixAbb = [0 1;-17 -8];
 desiredObserverPoles = [-2+2.377i -2-2.377i];
 observerGain = acker(MatrixAbb.',MatrixAab.',desiredObserverPoles.').'
 newMatrixAAbb = MatrixAbb - (observerGain*MatrixAab);
 newMatrixB = eye(2);
 newMatrixC = eye(2);
 newMatrixD = eye(2);
 mysys = ss(newMatrixAAbb,newMatrixB,newMatrixC,newMatrixD);
 initialX = [1 0];
 
 timeT = 0:.01:10;
 x = initial(mysys,initialX,timeT);
 x1 = [1 0]*x';
 x2 = [0 1]*x';
 
 figure(5)
 plot(timeT,x1,'r',timeT,x2,'g');
 title('Response to Initial Condition of State Variables Reduced Observer');
 xlabel('Time -->');
 ylabel('Magnitude -->');

%% Part 5) 
% Uncertainty
% using MATLAB's ureal
p = ureal('p', 0.01);

% We construct A again
I = I*p;
F0 = (K/R)*(u0 - K*x30)/(I+m*x20^2) - m*g*x20*cos(x10)/(I+m*x20^2) - 2*m*x20*x30*x40/(I+m*x20^2);

K1 = (K/R)/(I+m*x20^2);
K2 = m*g*x20*sin(x10)/(I+m*x20^2);

K3 = -(K/R)*(u0-K*x30)*2*m*x20/(I+m*x20^2)^2;
K4 = -(m*g*cos(x10) - m*g*x20*cos(x10)*2*m*x20/(I+m*x20^2)^2);
K5 = -(2*m*x40*x30*(I+m*x20^2) - 2*m*x20*x40*x30*2*m*x20)/(I+m*x20^2)^2;

K6 = (-K^2 /R)/(I+m*x20^2) - 2*m*x20*x40/(I+m*x20^2);

K7 = -2*m*x20*x30/(I+m*x20^2);

K8 = -5*g/7 *cos(x10);
K9 = x30^2;
K10 = 2*x20*x30;

A = [0 0 1 0;
     0 0 0 1;
     K2, (K3+K4+K5), K6, K7;
     K8, K9, K10, 0];
 
B = [0;0;K1;0];

usys = ss(A,B,C,D);
figure
step(usys), grid on
%nyquist(usys), title('Inertia Uncertainty Nyquist')

%% Part 6) This is contained inside the model for Part 1)

