clear;
clc;
%% Test Section
Xe_Mass = 1046; % Xenon mass (kg)
Xe_Cost = Xenon_Cost(Xe_Mass);% Xenon cost ($)
fprintf('Cost of propellant: $%.2f\n',Xe_Cost)

Tank_Vol = 0.0605; % Tank volume (m^3)
Cost = Tank_Cost(Tank_Vol);
fprintf('Cost of tug spherical propellant tank: $%.2f\n',Cost)
%% Propellant cost function
function Xe_Cost = Xenon_Cost(params_MP0)
T_STP = 288; % Temperature of air at sea level (K)
Rbar = 8314.3; % Universal gas constant (J/kmol-K)
MW = 131.2930; % Molecular weight of Xenon (kg/kmol)
RXe = Rbar/MW; % Gas constant of Xenon (J/kg-K)
p_internal = 1250*6894.76; % Pressure inside the tank (Pa)
rhoXe = p_internal/(RXe*T_STP); % Density of xenon (kg/m^3)
Xe_Vol = params_MP0/rhoXe; % Volume of xenon in tank (m^3)
c_Xe = 10/.001; % Cost density of Xenon ($/m^3)
Xe_Cost = c_Xe*Xe_Vol; % Cost of Xenon ($)
end

%% Propellant tank cost function
function Cost = Tank_Cost(Tank_V)
% tank_V is volume of tank
c_d = 20; % Cost density of titanium ($/kg)
rhoTi = 4430; % Density of Titanium (kg/m^3)
Tank_Mass = rhoTi*Tank_V; % Dry mass of tank (kg)
Cost = Tank_Mass*c_d; % Cost of tank ($)
end