%% Import CASADI package

addpath('C:\Users\Samuel\OneDrive - Cranfield University\Documents\Cranfield Course\Individual Research Project (IRP)\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

%% Track Generator

% Define track parameters
OvalTrackCornerRadius = 200;
OvalTrackStraightLength = 200;
OvalTrackWidth = 7.5;
OvalTrackResApprox = 20;

Track = ovalTrackGeneration(OvalTrackStraightLength, OvalTrackCornerRadius, OvalTrackWidth, OvalTrackResApprox);

%% Optimization problem initialization

N = Track.NOPoints; %Prediction horizon defined by the number of track points

opti = casadi.Opti(); % Optimization problem

%% Definition of decision variables

X = opti.variable(5,N+1); % state trajectory defined by the states x_k = (n_k, e_k, u_k, r_k, t_k)
n_k   = X(1,:);
e_k = X(2,:);
u_k = X(3,:);
r_k = X(4,:);
t_k = X(5,:);
U = opti.variable(2,N);   % control input defined by U_k = (u_dot_k, r_dot_k)
u_dot_k = U(1,:);
r_dot_k = U(2,:);
delta_U = opti.variable(2,N+1); %input increment vector decision variable

%% Cost function definition

%Weighting matrices definition (Q,R,P)
v_Q = [10^-5 10^-8 10^-8 10^-3 10^-8];
Q = diag(v_Q);
v_R = [0.01 0.01];
R = diag(v_R);
v_P = [10^-8 10^-8 10^-8 10^-8 1];
P = diag(v_P);

%Cost function summatory
J = sum(X(:,1:N+1)'*Q*X(:,1:N+1) + delta_U(:,1:N+1)'*R*delta_U(:,1:N+1));
J = sum(J);
J = J + X(:,N+1)'*P*X(:,N+1);

%minimize
opti.minimize(J); % cost function to minimize to obtain the optimal lap time

%% ---- dynamic constraints --------
C_k = Track.Curvature';

f = @(x,u,C) [((1-x(1)*C)/(x(3)*cos(x(2))))*(x(3)*sin(x(2)));...   n_k
    ((1-x(1)*C)/(x(3)*cos(x(2))))*x(4)-C;...  e_k
    ((1-x(1)*C)/(x(3)*cos(x(2))))*u(1);...   u_k
    ((1-x(1)*C)/(x(3)*cos(x(2))))*u(2);...   r_k
    ((1-x(1)*C)/(x(3)*cos(x(2))))];         % t_k                        dx/ds = f(x,u,C)

ds = Track.Res; % length of a control interval: space-dependent
for k=1:N % loop over control intervals
   % Euler Integration Method
   y_n = f(X(:,k),U(:,k),C_k(k));
   x_next = X(:,k) + ds*y_n; 
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

for k=1:N
    %input vector definition
    u_next = U(:,k-1) + delta_U(:,k);
    opti.subject_to(U(:,k)==u_next);
end

%% ---- Initial conditions -----------
% set the initial state of the racing vehicle
opti.subject_to(n_k(1)==0);  
opti.subject_to(e_k(1)==0); 
opti.subject_to(u_k(1)==50);
opti.subject_to(r_k(1)==0);
opti.subject_to(u_dot_k(1)==0);
opti.subject_to(r_dot_k(1)==0);
%Time needs to be 0 at the initial state
opti.subject_to(t_k(1)==0);
%% ---- boundary conditions --------
width = Track.Width(1); %track boundaries for the n_k state
maxAngle = (45/180)*pi;   %max angle with the track for state e_k
a_max = 9.2;    % max acceleration (from 0 to 100km/h in 3 seconds)
maxRDot = 0.3; % max yaw acceleration (rad/s^2)
v_max = 106; % max speed (m/s)


%boundary constraints for the states
opti.subject_to(-width<=n_k<=width);
opti.subject_to(-maxAngle<=e_k<=maxAngle); 

%boundary constraints for the inputs
opti.subject_to(0<=u_dot_k<=a_max);
opti.subject_to(0<=u_k<=v_max);
opti.subject_to(-maxRDot<=r_dot_k<=maxRDot);


%% ---- Other constraints  ----------
opti.subject_to(t_k>=0); % Time must be positive

%% ---- initial values for solver ---
opti.set_initial(n_k, 0);
opti.set_initial(e_k, 0);
opti.set_initial(u_k, Track.Res);
opti.set_initial(r_k, 0);
opti.set_initial(u_dot_k, 0);
opti.set_initial(r_dot_k, 0);
opti.set_initial(t_k, 1);

%% ---- solve NLP ------

p_opts = struct('expand',true);
s_opts = struct('max_iter', 10000);
opti.solver('ipopt', p_opts, s_opts); % set numerical backend
sol = opti.solve();   % actual solve
%% ---- post-processing        ------

%Extract the solution for the defined states
sol_n = sol.value(n_k); %lateral displacement
sol_eps = sol.value(e_k); %heading angle error
sol_u = sol.value(u_k);   %longitudinal velocity
sol_r = sol.value(r_k);   %angular velocity
sol_u_dot = sol.value(u_dot_k);  %longitudinal acceleration
sol_r_dot = sol.value(r_dot_k);  %yaw acceleration
sol_t = sol.value(t_k);       %elapsed time

%Extract the solution in Cartesian coordinates
len = length(sol_n);
for i=1:len
    if i<=N
        sol_yaw(i) = Track.Heading(i)' + sol_eps(i);
        sol_x(i) = Track.CentrelineCoord(i,1)' - sol_n(i).*sin(Track.Heading(i)');
        sol_y(i) = Track.CentrelineCoord(i,2)' + sol_n(i).*cos(Track.Heading(i)');
    else
        sol_yaw(i) = Track.Heading(1)' + sol_eps(i);
        sol_x(i) = Track.CentrelineCoord(1,1)' - sol_n(i).*sin(Track.Heading(1)');
        sol_y(i) = Track.CentrelineCoord(1,2)' + sol_n(i).*cos(Track.Heading(1)');
     end 
end

%plot the track
figure
Track.CentrelineCoord(N+1,:) = [Track.CentrelineCoord(1,1), Track.CentrelineCoord(1,2)];
plot(Track.LeftEdgeCoord(:,1), Track.LeftEdgeCoord(:,2), '+',...
    Track.RightEdgeCoord(:,1), Track.RightEdgeCoord(:,2), '*',...
    sol_x, sol_y, 'k',...
    sol_x(1),sol_y(1), 'bo')
title('Shortest Time Manoeuvre')
axis equal
legend('Trajectory Line', 'Vehicle Starting Point')
%Longitudinal speed
figure
hold on
plot(sol_t,sol.value(u_k));
xlabel('Time [s]');
ylabel('speed [m/s]')
title('Speed profile u_k(m/s)')

%Lateral displacement
figure
plot(sol_t,sol.value(n_k))
xlabel('Time [s]');
ylabel('lateral displacement [m]')
title('Lateral displacement profile n_k(m)')

%Longitudinal acceleration 
figure
plot(sol_t(1:end-1),sol.value(u_dot_k))
xlabel('Time [s]');
ylabel('longitudinal acceleration [m/s^2]')
title('Longitudinal acceleration u_dot_k(m/s^2)')



