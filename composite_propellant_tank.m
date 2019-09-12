
% Assume STP
%% Givens
p_i = 1250*6894.76; % Pressure inside the tank (Pa)
design_choice = 1/3;% Liner can be designed to hold between 1/2 the internal pressure and 1/3 the internal pressure
p_l = p_i*design_choice; % Pressure on the liner (Pa)
p_c = p_i*(1 - design_choice); % Pressure on the composite (Pa)
T_STP = 288; % Temperature of air at sea level (K)
Rbar = 8314.3; % Universal gas constant (J/kmol-K)
MW = 131.2930; % Molecular weight of Xenon (kg/kmol)
RXe = Rbar/MW; % Gas constant of Xenon (J/kg-K)
Xe_mass = 1000; % Mass of Xenon (kg)
rho_Ti = 4430; % Density of Titanium (kg/m^3)
FS = 2; % Factor of safety
stress_Yield_Ti = 880e6; % Yield stress of Titanium (Pa)
stress_Max_Ti = stress_Yield_Ti/FS; % Maximum allowable stress of Titanium (Pa)
stress_Yield_f = 3.53e9; % Yield stress of Carbon Fiber (Pa)
stress_Max_f =  stress_Yield_f/FS;% Maximum stress of Carbon Fiber strand (Pa)
rhoXe = 1700; % Density of xenon (kg/m^3)
Xe_Vol = Xe_mass/rhoXe; % Volume of the Xenon (m^3)
nu_f = 0.3; % poisson's ratio of Carbon Fiber (n.d.)
rho_Ca = 896; % Density of Carbon Fiber (kg/m^3)

%% Calculations
rl = (3*Xe_Vol/(4*pi))^(1/3); % Inner liner radius (m)
tl = p_l*rl/(2*stress_Max_Ti); % Liner thickness (m)
rc = rl + tl; % Inner composite radius (m)
tf = p_c*rc/stress_Max_f; % Fiber radius (m)
stress_Max_Ca = nu_f*stress_Max_f/2; % Stress of Carbon Fiber composite
tc = p_c*rc/(2*stress_Max_Ca); % Composite thickness (m)
rol = rl + tl; % Outer liner radius (m)
roc = rc + tc; % Outer composite radius (m)
Tank_Vol_Ti = ((4/3)*pi*(rol^3)) - Xe_Vol; % Volume of titanium liner (m^3)
Tank_Vol_Ca = ((4/3)*pi*(roc^3)) - (Tank_Vol_Ti + Xe_Vol); % Volume of carbon fiber (m^3)
Tank_Mass_Ti = rho_Ti*Tank_Vol_Ti;
Tank_Mass_Ca = rho_Ca*Tank_Vol_Ca;

%% Cost Calculation
cd_Ti = 20; % Cost density of titanium ($/kg)
Cost_Ti = Tank_Mass_Ti*cd_Ti; % Cost of titanium liner ($)
cd_Ca = 5*171; % Cost density of carbon fiber ($/kg)
Cost_Ca = cd_Ca*Tank_Mass_Ca; % Cost of carbon fiber ($)
Cost = Cost_Ca + Cost_Ti; % Cost of tank ($)
%% Output
fprintf('Dry mass of tug  propellant tank: %.4f kg\n', Tank_Mass_Ti + Tank_Mass_Ca)
fprintf('Diameter of tug  propellant tank: %.4f m\n', 2*roc)
fprintf('Wall Thickness of carbon fiber: %.4f m\n', tc)
fprintf('Wall Thickness of Titanium liner: %.4f m\n', tl)
fprintf('Cost of propellant tank: $%.2f\n',Cost)