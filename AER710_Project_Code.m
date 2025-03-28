clear; clc; close all;
 
% 554 708 990
% i = 557
% j = 710
% k = 996
% Mn = 1.4949

% Akus Chhabra
% 500970974

%% Parameter Definitions

num_shock = 4; % 3 oblique shocks + 1 normal shock
M1 = 3.2; % Downstream Mach
M4 = 1.3; % Upstream Mach
gamma = 1.4;

%% Main Code
% Compute M(2) to M(n-1), B(2) to B(n-1) and θ(2) to θ(n-1), stagnation
% pressure across the individual oblique shock, and normal shock π(1) to
% π(n-1)
%
% Objective is to maximize π(d).

inc = 0.0001;
Mn_norm = 0;
beta = asind(1/M1);
Stag_ratio4 = ((((gamma+1)*M4^2)/((gamma-1)*M4^2+2))^((gamma)/(gamma-1))) * (((gamma+1)/(2*gamma*M4^2-(gamma-1)))^(1/(gamma-1)));
while Mn_norm ~= M4
    [theta, Stag_ratio, M2] = Oblique(M1, gamma, beta);
    beta2 = asind(M1*sind(beta)/M2);
    [theta2, Stag_ratio2, M3] = Oblique(M2, gamma, beta2);
    beta3 = asind(M2*sind(beta2)/M3);
    [theta3, Stag_ratio3, Mn_norm] = Oblique(M3, gamma, beta3);
    if Mn_norm > M4
        beta = beta + inc;
    else
        break
    end
end

%% Print Results
fprintf('\n\nM1 = %.3f', M1)
fprintf('\nM2 = %.3f', M2)
fprintf('\nM3 = %.3f', M3)
fprintf('\nM1n = M2n = M3n = %.3f', M1*sind(beta))
fprintf('\nM2n = %.3f', M2)
fprintf('\nM3n = %.3f', M3)
fprintf('\nBeta1 = %.2f deg', beta)
fprintf('\nBeta2 = %.2f deg', beta2)
fprintf('\nBeta3 = %.2f deg', beta3)
fprintf('\nTheta1 = %.2f deg', theta)
fprintf('\nTheta2 = %.2f deg', theta2)
fprintf('\nTheta3 = %.2f deg', theta3)
fprintf('\nPo2/Po1 = %.3f', Stag_ratio)
fprintf('\nPo3/Po2 = %.3f', Stag_ratio2)
fprintf('\nPo4/Po3 = %.3f', Stag_ratio3)
fprintf('\nPo5/Po4 = %.3f', Stag_ratio4)
fprintf('\nPressure Recovery = %.3f\n\n', Stag_ratio*Stag_ratio2*Stag_ratio3*Stag_ratio4)

%% Rendition of Inlet Geometry

beta_list = [0 beta beta2 beta3];
theta_list = [0 theta theta2 theta3];
x = [0 0.333 0.666 1]';
y = zeros(size(x,2),1);
for i = 1:length(x)
    y(i) = x(i)*tand(theta_list(i));
    y_flat = [y(1) y(1) y(1) y(1)];
end

figure(1)
plot(x,y,'color', 'black')
hold on
plot(x,y_flat,'color', 'black')
plot([x(2) x(2)],[y(1) y(2)],'color', 'black')
plot([x(3) x(3)],[y(1) y(3)],'color', 'black')
plot([x(4) x(4)],[y(1) y(4)],'color', 'black')
title('Rendition of Inlet Geometry')

delX = 0.15; % Change Cowl Edge  location

m_shock = -1/((y(4)-y(3))/(x(4)-x(3)));
b = y(4) - m_shock*x(4);
plot([(x(3)+delX) x(4)], [m_shock*(x(3)+delX)+b m_shock*x(4)+b],'color', 'red', 'linestyle', '--'); % Normal Shock
plot([x(3)+delX x(4)+delX],[m_shock*(x(3)+delX)+b ((-1/m_shock)*x(4)+b)*sind(theta_list(4))],'color', 'black'); % Cowl Edge
plot([x(1) (x(3)+delX)],[y(1) m_shock*(x(3)+delX)+b],'color', 'red', 'linestyle', '--') % Oblique Shock 1
plot([x(2) (x(3)+delX)],[y(2) m_shock*(x(3)+delX)+b],'color', 'red', 'linestyle', '--') % Oblique Shock 2
plot([x(3) (x(3)+delX)],[y(3) m_shock*(x(3)+delX)+b],'color', 'red', 'linestyle', '--') % Oblique Shock 3

hold off
ylim([-0.1 0.9])
xlim([-0.1 1.1])
grid on
grid minor

text(x(1)+0.07,y(1)+0.05,'β1')
text(x(1)+0.2,y(1)+0.025,'θ1')

text(x(1)+0.42,y(1)+0.15,'β2')
text(x(1)+0.5,y(1)+0.1,'θ2')

text(x(1)+0.73,y(1)+0.25,'β3')
text(x(1)+0.8,y(1)+0.2,'θ3')

%% Function Definition

function [theta, Stag_ratio, M1] = Oblique(M, gamma, beta)
    A = (((gamma+1)*M^2*(sind(beta))^2)/((gamma-1)*M^2*(sind(beta))^2+2))^(gamma/(gamma-1));
    B = (((gamma+1))/(2*gamma*M^2*(sind(beta))^2-(gamma-1)))^(1/(gamma-1));
    Stag_ratio = A*B;
    C = (gamma+1)*M^2;
    D = 2*(M^2*(sind(beta))^2-1);
    theta = acotd(tand(beta)*(C/D-1));
    M1_num = (gamma-1)*M^2*(sind(beta))^2 + 2;
    M1_denom = (2*gamma*M^2*(sind(beta))^2-(gamma-1)) * (sind(beta-theta))^2;
    M1 = sqrt(M1_num/M1_denom);
end
