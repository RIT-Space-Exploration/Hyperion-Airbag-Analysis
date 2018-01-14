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
v_imp = 9; % m/s - Impact Velocity

%% Constants
g_0 = 9.8066; % m/s^2 - Standard Gravity
mass = 4; % kg - Mass of system
dampening_coefficient = 1; % kg/s - Dampening coefficient
stiffness = 4; % kg/s^2
angular_frequency = sqrt(stiffness/mass);
sign = dampening_coefficient^2/(4*mass^2) - angular_frequency^2;
relax_time = 2*mass/dampening_coefficient;

if sign > 0
    disp('This creates a non-complex solution. Please fix the system constants so the sign is negative.');
end

if sign < 0
    %% Solution to Differential Eqn Roots
    alpha_1 = -dampening_coefficient/(2*mass)+1i*sqrt(angular_frequency^2-dampening_coefficient^2/(4*mass)); % r_1
    
    alpha_2 = -dampening_coefficient/(2*mass)-1i*sqrt(angular_frequency^2-dampening_coefficient^2/(4*mass)); % r_2
    
    %% Solution to Differential Equation
    %     t(1) = 0;
    j = 1;
    %     dt = 0.001; % s - Time step
    %     x(1) = 0;
    
    for t = 0:0.001:20
        c_1 = v_imp/(alpha_2-alpha_1); % Constant of Differentiation 1 (c_1)
        
        c_2 = -v_imp/(alpha_2-alpha_1); % Constant of Differentiation 2 (c_2)
        
        x(j) = c_1*exp(alpha_1*t) + c_2*exp(alpha_2*t); % m - Position (Extension (+) , Compression (-))
        
        v(j) = alpha_1*c_1*exp(alpha_1*t) + alpha_2*c_2*exp(alpha_2*t); % m/s - Velocity
        
        a(j) = alpha_1^2*c_1*exp(alpha_1*t) + alpha_2^2*c_2*exp(alpha_2*t); % m/s^2 - Acceleration
        
        f(j) = mass*a(j); % Force
        
        g(j) = a(j)/g_0; % G-Force
        
        j = j + 1;
    end
    
    t = 0:0.001:20;
    
    %% Plots
    grid on
    
    figure()
    plot(t,x);
    title('Position vs. Time');
    ylabel('Position [m]');
    xlabel('Time [s]');
    
    figure()
    plot(t,v);
    title('Velocity vs. Time');
    ylabel('Velocity [m/s]');
    xlabel('Time [s]');
    
    figure()
    plot(t,a);
    title('Acceleration vs. Time');
    ylabel('Acceleration [m/s^2]');
    xlabel('Time [s]');
    
end
