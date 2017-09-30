% Title: SPEX IREC Spring Script
% Author: James Emerson Parkus
% Date: 9/29/17
% Purpose: First order approximation to understanding the requirements for
% the balloon airbags of Hyperion. The solution is attained by modeling the
% system as a spring. The full equation is m*a + r*v + s*x = 0. The damping
% coefficient, r, is taken from the material and system properties. The
% mass of the system, m, is set at 4kg due to payload restricitions. The
% spring coefficient (a.k.a. Stiffness constant), s, is taken as a range of
% values for which the motion of the system will be analyzed. We set the
% initial conditions of the situation (impact from freefall) as the impact
% velocity at time = 0, and the time at which the velocity is 0 m/s. The
% impact velocity and the time to rest are given by the user.
%
% Assumptions:
% 1. Linear system (Not Chaotic)
% 2. m, r, & s are constant
% 3. There is no driving force

clc
clear

%% Initial Conditions
v_imp = 31.29; % m/s - Impact Velocity
t_r = 2; % s - Time to rest (after impact)

%% Constants
g_0 = 9.8066; % m/s^2 - Standard Gravity
m = 4; % kg - Mass of system
r = 100; % kg/s - Dampening coefficient

i_n = 5000; % Number of Terms - 1
s_min = 0.01; % Lower Bound
s_max = 60; % Upper Bound
n = (s_max-s_min)/i_n; % Iteration Step

s = s_min:n:s_max; % kg/s^2 - Spring Coefficient

%% Solution to Differential Eqn Constants
alpha_1 = (-r+sqrt(r^2 - 4*m*s))/(2*m); % r_1

alpha_2 = (-r-sqrt(r^2 - 4*m*s))/(2*m); % r_2

%% Solution to Differential Equation
t(1) = 0;

dt = (t_r - t(1))/i_n; % Iteration Step

for j = 1:1:5000
    
    alp_1 = alpha_1(j);
    alp_2 = alpha_2(j);
    
    for i = 1:1:5000
        
        c_1 = v_imp*(1 + exp(2*alp_1))/(alp_1*alp_2*(exp(2*alp_2) - exp(2*alp_1))); % Constant of Differentiation 1 (c_1)
        
        c_2 = -v_imp*exp(2*alp_1)/(alp_2*(exp(2*alp_2) - exp(2*alp_1))); % Constant of Differentiation 2 (c_2)
        
        x(i,j) = c_1*exp(alp_1*t(i)) + c_2*exp(alp_2*t(i)); % m - Position (Extension (+) , Compression (-))
        
        v(i,j) = alp_1*c_1*exp(alp_1*t(i)) + alp_2*c_2*exp(alp_2*t(i)); % m/s - Velocity
        
        a(i,j) = alp_1^2*c_1*exp(alp_1*t(i)) + alp_2^2*c_2*exp(alp_2*t(i)); % m/s^2 - Acceleration
        
        g(i,j) = a(i,j)/g_0;
        
        t(i+1) = t(i) + dt;
        
    end
    
end

%% Plots

figure()
plot(t(1:5000),abs(x(1:5000,5000)));
title('Position vs. Time');
ylabel('Position [m]');
xlabel('Time [s]');

figure()
plot(t(1:5000),abs(v(1:5000,5000)))
title('Velocity vs. Time');
ylabel('Velocity [m/s]');
xlabel('Time [s]');

figure()
semilogy(t(1:5000),abs(g(1:5000,5000)));
title('Acceleration vs. Time');
ylabel('Acceleration [m/s^2]');
xlabel('Time [s]');
grid on

figure()
plot(t(1:5000),abs(g(1:5000,5000)));
title('Acceleration vs. Time');
ylabel('Acceleration [m/s^2]');
xlabel('Time [s]');
grid on












































