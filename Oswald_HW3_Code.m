clear;
clc;
close all

gamma = 1.4;
P1 = 1;
T1 = 298.15;
theta_c = 5;

Mach = linspace(4, 20, 17);

for i = 1:17
    M = Mach(i);
    tau = 5.1;
    f = -999;
        while abs(f(end,1)) > 0.005
            tau = tau+0.01;
            b = theta_c/tau;
            K = M * (tau*pi/180);
            global omega;
            omega = (2*gamma*K^2-(gamma-1))/(gamma*(gamma+1)*K^2) * ((2+(gamma-1)*K^2)/((gamma+1)*K^2))^gamma;
            f = 1/2;
            fp = ((gamma+1)*K^2)/(2+(gamma-1)*K^2);
            f0 =[f fp];
            thetaspan = [1 b];
            [theta,f] = ode45(@fun,thetaspan,f0);
        end
    Cp = 2 * ((omega/b^2) * (f(end,2)/b)^gamma - 1/(gamma*M^2*(theta_c*pi/180)^2)) * (theta_c*pi/180)^2;
    P(i) = Cp*gamma*M^2/2 + 1;
end

plot(Mach,P)
hold on

%% Part 2
theta_c = theta_c * pi/180;
for i = 1:17
theta_s = theta_c;
M1 = Mach(i);
V_check = -999;
while abs(V_check) > 0.0005
theta_s = theta_s + 0.0001;
delta = atan(2*cot(theta_s) * ((M1^2*sin(theta_s)^2-1)/(M1^2*(gamma+cos(2*theta_s))+2)));
if delta < 0
continue
end
Mn1 = M1 * sin(theta_s);
Mn2 = sqrt((Mn1^2 + (2/(gamma-1)))/((2*gamma/(gamma-1)*Mn1^2)-1));
M2 = Mn2 / (sin(theta_s - delta));
Po1_p1 = (1+(((gamma-1)/2)*M1^2))^(gamma/(gamma-1));
Po2_Po1=((((gamma+1)/2*Mn1^2)/(1+(((gamma-1)/2)*Mn1^2)))^(gamma/(gamma-1)))/((((2*gamma/(gamma+1))*Mn1^2)-((gamma-1)/(gamma+1)))^(1/(gamma-1)));
Po2 = Po2_Po1 *Po1_p1;
V = (2/((gamma-1)*M2^2)+1)^-0.5;
V_theta = -V*sin(theta_s-delta);
V_r = V*cos(theta_s-delta);
v0 = [V_r V_theta];
thetaspan2 = linspace(theta_s, theta_c);
[theta,v] = ode45(@fun2, thetaspan2, v0);
V_check = v(end,2);

end

M = sqrt(2/(((v(end,1)^(-2))-1)*(gamma-1)));
Po_P_surface=(1+(((gamma-1)/2)*M^2))^(gamma/(gamma-1));
Pf_2(i) = Po2/Po_P_surface;
end


plot(Mach, Pf_2)

xlabel("Mach Number")
ylabel("Pressure at the surface (atm)")
title("Mach Number v. Pressure at the Surface")
legend("Small Disturbance Theory", "Inviscid Exact Theory", 'location', 'northwest')
hold off


%% Functions
function output = fun(theta, f)
global omega;
gamma = 1.4;
num = 2*f(1)*f(2)^2-gamma*omega*(f(2)^(gamma+2)/(theta^gamma));
den = 4*f(1)^2 - gamma*omega*(f(2)^(gamma+1))/(theta^(gamma-1));
output = [f(2);  num/den];
end

function output = fun2(theta, v)
gamma = 1.4;
vt_der = (v(2)^2*v(1) - (gamma-1)/2 * (1-v(1)^2-v(2)^2) * (2*v(1)+v(2)*cot(theta)))/((gamma-1)/2 * (1-v(1)^2-v(2)^2) -v(2)^2);
output = [v(2); vt_der];
end



