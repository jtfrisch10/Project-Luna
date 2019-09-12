clc;
clear;
%% Launch load calculations
% This code calculates whether the current Tug geometry and mass will break
% under launch loads

%% Inputs
%% Given
g_lateral = [0 0.5 0.5 1 2 2 2 2 2 2 2 1 0.5 0.5 0.5 0 -0.5 -0.5 -0.5 -1 -2 -2 -2 -2 -2 -2 -2 -1 -0.5 -0.5]; % g's in lateral direction
g_axial = [-2 -2 -1.5 -1.5 -1.5 -1 0 1 2 3 3.5 23/6 4 5 6 6 6 5 4 23/6 3.5 3 2 1 0 -1 -1.5 -1.5 -1.5 -2]; % g's in axial direction
g = 9.81; % acceleration of gravity (m/s^)
a_ax = g_axial*9.81; % axial acceleration (m/s^2)
a_lat = g_lateral*9.81; % Lateral acceleration (m/s^2)
a_c = linspace(0,30.8425,length(a_ax)); % centripetal acceleration (m/s^2)
R = 1.3369; % Radius of Tug (m)
L = 3.9837; % Length of Tug
t = 0.0026; % Thickness of skin (m)
b = 0.0011; % Width of rib web (m)
d = 0.006; % Depth of web (m)
c = 0; % Depth of flange (m)
w = 0; % Width of flange (m)
a = 0.2; % leg of triangle
h = a*sqrt(3)/2; % height of triangle
delta = d/t; % Non-dimensional parameter of isogrid
lambda = c/t; % Non-dimensional parameter of isogrid
alpha = b*d/(t*h); % Non-dimensional parameter of isogrid
mu = w*c/(t*h); % Non-dimesional parameter of isogrid
beta = sqrt((1 + alpha + mu)*(3*(1 + delta)^2 + 3*mu*(1 - lambda)^2 + 1 + alpha*delta^2 + mu*lambda^2) - 3*((1 + delta) - mu*(1 + lambda))^2);
m_Tug = 3.1971e3; % Mass of Tug (kg)
FS = 2; % Factor of safety
E = 71.7e9; % Young's Modulus of Aluminum 7075 (Pa)
tstar = t*beta/(1 + alpha + mu); % Equivalent thickness of isogrid
Estar = E*(1 + alpha + mu)^2/beta; % Equivalent Young's Modulus
sigma_yield = 503e6; % Yield stress of Aluminum 7075 (Pa)
sigma_allow = sigma_yield/FS; % Maximum allowable stress (Pa)
nu = 0.33; % Poissions ratio of aluminum
rhoAl = 2810; % Density of aluminum (kg/m^3)

%% Cross Section Calculations 
I = (pi/4)*(R^4 - (R-t)^4); % Area moment of inertia  (m^4)
A = pi*(R^2 - (R-t)^2); % Area of Tug (m^2)

%% Buckling Calculations
Pcr = pi^2*E*I/(4*L^2); % Critical buckling force (N)

%% General Instability Calculations
count = 0;
for k = 1:length(a_ax)
    P = -a_ax*m_Tug; % axial force on Tug (N)
    Mz = a_c(k)*m_Tug*L/2; % Intertial moment caused by turning (Nm)
    V = m_Tug*a_lat; % 
    sigma_xx = P/A - Mz/I; % Normal stress on Tug
    Tau_xy = V/A; % Shear stress on Tug
    sigma_max = sigma_xx./2 + sqrt((sigma_xx/2).^2 + Tau_xy.^2); % Max principle stress (Pa)
    sigma_min = sigma_xx./2 - sqrt((sigma_xx/2).^2 + Tau_xy.^2); % Min principle stress (Pa)
    
    sigma_v = sqrt(sigma_max.^2 + sigma_min.^2 - sigma_max.*sigma_min);

    %% Check Tug Against launch loads

    for i = 1:length(a_ax)
        if abs(P(i)) > Pcr
            count = count + 1;
        elseif abs(sigma_max(i)) > sigma_allow
            count = count + 1;
        elseif abs(sigma_v) > sigma_allow
            count = count + 1;
        elseif abs(sigma_min(i)) > sigma_allow
            count = count + 1;
        end
    end
end
if count >= 1
    disp('Tug will buckle during launch')
else
    disp('Tug will not break due to launch loads')
end

%% Acoustic Analysis

ma = m_Tug/L; % Mass per unit length (kg/m)
c = [1.875, 4.694, 7.855, 11.00, 14.14]; % Constant specific to first 5 modeshpaes of a cantilevered beam 
d = [0.734096, 1.018467, 0.999224, 1.000034, 0.999999]; % Constant specific to first 5 modeshpaes of a cantilevered beam 
fn = c.^2/(2*pi*L^2)*sqrt(E*I/ma);
if fn(1) <= 50
    disp('Further analysis is required to determine launch capabilities')
end






