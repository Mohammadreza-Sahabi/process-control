clc
clear
close all

%% Static Simiulation =====================================================
% Data & Relations---------------------------------------------------------

global R MM_pent omega_pent Tc_pent Pc_pent Zc_pent Vc_pent
global MM_but omega_but Tc_but Pc_but Zc_but Vc_but
global F_rate F_temp F_zC4 F_zC5 D_rate D_zC4 D_zC5 B_rate B_zC4 B_zC5 RR
global Col_ID Col_trays Total_stages pd Lw A_available h_weir L_t2t Col_h
global C_P C_T C_duty R_P R_T R_duty
global alpha p1 p2 p3 p4
global Y rho_func rho_nbutane rho_npentane Cp_ig_but Cp_ig_pent
global Cp_liq_func Cp_liq_but Cp_liq_pent P_sat_but P_sat_pent
global T_sat_but T_sat_pent P_mix rho_avg Mw_avg how F_L T_mix V

R = 8.3143; % J/(mo1.K)

% constants for n-Pentane:
MM_pent = 72.150e-3; % kg/mol   % molar_mass_pentane
omega_pent = 0.252; % dimensionless % acentric_factor_pentane
Tc_pent = 469.7; % K    % critical_temperature_pentane
Pc_pent = 33.70 * 1e5; % Pa  % critical_pressure_pentane (converted from bar to Pa)
Zc_pent = 0.270; % dimensionless    % critical_compressibility_factor_pentane
Vc_pent = 313.0 * 1e-6; % m³/mol  % critical_molar_volume_pentane (converted from cm³/mol to m³/mol)

% constants for n-Butane
MM_but = 58.123e-3; % kg/mol    % molar_mass_butane
omega_but = 0.200; % dimensionless  % acentric_factor_butane
Tc_but = 425.1; % K % critical_temperature_butane
Pc_but = 37.96 * 1e5; % Pa   % critical_pressure_butane (converted from bar to Pa)
Zc_but = 0.274; % dimensionless % critical_compressibility_factor_butane
Vc_but = 255.0 * 1e-6; % m³/mol % critical_molar_volume_butane (converted from cm³/mol to m³/mol)

% STEADY STATE:------------------------------------------------------------
% Feed stream properties
F_rate = 100 * (1000/3600); % mol/s (converted from Kmol/h)
F_temp = 320; % K
F_zC4 = 0.50; % Mole fraction of n-Butane
F_zC5 = 0.50; % Mole fraction of n-Pentane

% Distillate stream properties
D_rate = 50 * (1000/3600); % mol/s (converted from Kmol/h)
D_zC4 = 0.99; % Mole fraction of n-Butane
D_zC5 = 0.01; % Mole fraction of n-Pentane

% Bottoms stream properties
B_rate = 50 * (1000/3600); % mol/s (converted from Kmol/h)
B_zC4 = 0.01; % Mole fraction of n-Butane
B_zC5 = 0.99; % Mole fraction of n-Pentane

% Reflux ratio
RR = 1.323; % Reflux ratio

% Column properties
Col_ID = 0.7145; % m
Col_trays = 59; % Total number of trays
Total_stages = 61; % Total number of equilibrium stages
pd = 689.476; % Pa (0.1 psi)    % Pressure drop per tray
Lw = 0.6*Col_ID; % m    % length of weir
A_available = 0.93*pi*(Col_ID^2)/4; % m^2   
h_weir = 0.06; %m
L_t2t = 0.6096; %m  % Tray to tray distance
Col_h = 41.35;  %m  % Column height

% Condenser properties
C_P = 4.5 * 101325; % Pa (converted from 4.5 atm)
C_T = 322; % K
C_duty = 0.6322e6; % W
C_Diameter = 0.6;   % m
C_area = pi*(C_Diameter^2)/4;   % m^2

% Reboiler properties
R_P = 4.9 * 101325; % Pa (converted from 4.9 atm)
R_T = 367; % K
R_duty = 0.758e6; % W
R_Diameter = 0.9;   % m
R_area = pi*(R_Diameter^2)/4;   % m^2

% Calculation of alpha-----------------------------------------------------

% Data was obtained from Hysys:
disp('srk')
xd1 = [0.00E+00, 2.00E-02, 4.00E-02, 6.00E-02, 8.00E-02, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1];
yd1 = [0, 4.47E-02, 8.76E-02, 0.128838519, 0.16841662, 0.206412796, 0.242871606, 0.277893943, 0.311527324, 0.343820544, 0.374807497, 0.404587191, 0.433192545, 0.460684863, 0.487095714, 0.512457173, 0.536843098, 0.560302082, 0.582850308, 0.604532772, 0.62538813, 0.645471184, 0.66479037, 0.683389008, 0.701297877, 0.718546914, 0.735164694, 0.75117849, 0.766614325, 0.781497022, 0.795850261, 0.809696626, 0.823057657, 0.835953895, 0.848410304, 0.860434716, 0.872050427, 0.883274429, 0.89412293, 0.904611391, 0.91475454, 0.924566499, 0.934062989, 0.943252138, 0.952148546, 0.960763893, 0.969109331, 0.97719551, 0.985032604, 0.992630333, 1];
alphaFit(xd1, yd1)

disp('peng robinson')
xd2 = [0.000000000000000, 2.00000000000000e-002, 4.00000000000000e-002, 6.00000000000000e-002, 8.00000000000000e-002, 0.100000000000000, 0.120000000000000, 0.140000000000000, 0.160000000000000, 0.180000000000000, 0.200000000000000, 0.220000000000000, 0.240000000000000, 0.260000000000000, 0.280000000000000, 0.300000000000000, 0.320000000000000, 0.340000000000000, 0.360000000000000, 0.380000000000000, 0.400000000000000, 0.420000000000000, 0.440000000000000, 0.460000000000000, 0.480000000000000, 0.500000000000000, 0.520000000000000, 0.540000000000000, 0.560000000000000, 0.580000000000000, 0.600000000000000, 0.620000000000000, 0.640000000000000, 0.660000000000000, 0.680000000000000, 0.700000000000000, 0.720000000000000, 0.740000000000000, 0.760000000000000, 0.780000000000000, 0.800000000000000, 0.820000000000000, 0.840000000000000, 0.860000000000000, 0.880000000000000, 0.900000000000000, 0.920000000000000, 0.940000000000000, 0.960000000000000, 0.980000000000000, 1.00000000000000];
yd2 = [0.000000000000000, 4.45322292746566e-002, 8.72980656248920e-002, 0.128363557461506, 0.167806379996668, 0.205678528106102, 0.242024985263500, 0.276945317321689, 0.310487091123557, 0.342698926304944, 0.373614834557850, 0.403332357344191, 0.431884479466097, 0.459331826880687, 0.485705967603598, 0.511038727108430, 0.535403125879835, 0.558847257475426, 0.581387234445524, 0.603067636726287, 0.623926679360559, 0.644018623431191, 0.663351736044603, 0.681969062314683, 0.699900983660458, 0.717177131701022, 0.733825792643207, 0.749873962326957, 0.765347400331117, 0.780270682937779, 0.794667254793448, 0.808559479140460, 0.821968686522462, 0.834915221894026, 0.847423854667073, 0.859502264245044, 0.871173572057965, 0.882454618525231, 0.893361469047116, 0.903909426281262, 0.914113173233437, 0.923986631296062, 0.933545438292649, 0.942797635501670, 0.951757721407837, 0.960437288220692, 0.968847407135491, 0.976998653747820, 0.984901132241064, 0.992564498293839, 1.00000000000000];
alphaFit(xd2, yd2)

alpha = 2.5; % dimensionless
% y-x in saturation
Y = @(x) alpha .* x ./ (1 + (alpha - 1) .* x);

% Density:-----------------------------------------------------------------

rho_func = @(Vc, T, Tc, omega) (1/Vc) * (1 + 0.85.*(1 - T./Tc) + (1.6916 + 0.984.*omega).*(1 - T./Tc).^(1/3));
% Specific rho_func for each component
rho_nbutane = @(T) (1/Vc_but) .* (1 + 0.85.*(1 - T./Tc_but) + (1.6916 + 0.984.*omega_but).*(1 - T./Tc_but).^(1/3));
rho_npentane = @(T) (1/Vc_pent) .* (1 + 0.85.*(1 - T./Tc_pent) + (1.6916 + 0.984.*omega_pent).*(1 - T./Tc_pent).^(1/3));

% Heat capacities:---------------------------------------------------------

%These function handles calculate C_p^ig for the given temperature T in Kelvin. The valid temperature range is from 298 K to 1500 K (T_max) for both compounds.
Cp_ig_but = @(T) R.*(1.935 + 36.915e-3.*T + (-11.402e-6).*T.^2 + 0.*T.^(-2));
Cp_ig_pent = @(T) R.*(2.464 + 45.351e-3.*T + (-14.111e-6).*T.^2 + 0.*T.^(-2));

%The Rowlinson-Bondi method can be used to estimate liquid heat capacities.This method is based on the ideal gas heat capacity and the corresponding-states principle.
Cp_liq_func = @(Cp_id, R, T, Tc, w) Cp_id + 1.45*R + (0.45*R)./(1 - T/Tc) + 0.25*w*R.*(17.11 + (25.2.*(1 - T./Tc).^(1/3))/(T./Tc) + 1.742./(1 - T./Tc));
% Specific Cp_liq_func for each component
Cp_liq_but = @(T) Cp_ig_but(T) + 1.45*R + (0.45*R)./(1 - T./Tc_but) + 0.25*omega_but*R.*(17.11 + (25.2.*(1 - T./Tc_but).^(1/3))/(T./Tc_but) + 1.742./(1 - T./Tc_but));
Cp_liq_pent = @(T) Cp_ig_pent(T) + 1.45*R + (0.45*R)./(1 - T./Tc_pent) + 0.25*omega_pent*R.*(17.11 + (25.2.*(1 - T./Tc_pent).^(1/3))/(T./Tc_pent) + 1.742./(1 - T./Tc_pent));

% Antoine Eqn calculate P_sat & T_sat. The valid temperature ranges are: For n-butane: -73 to 19 °C For n-pentane: -45 to 58 °C
P_sat_but = @(T) 1e3.*exp(13.6608 - 2154.70 ./ ((T - 273.15) + 238.789));
P_sat_pent = @(T) 1e3.*exp(13.7667 - 2451.88 ./ ((T - 273.15) + 232.014));
T_sat_but = @(P) 2154.70 ./ (13.6608 - log(P*1e3)) - 238.789 + 273.15;
T_sat_pent = @(P) 2451.88 ./ (13.7667 - log(P*1e3)) - 232.014 + 273.15;

% Raoult's Law-------------------------------------------------------------
P_mix = @(x , T) x.*P_sat_but(T) + (1-x).*P_sat_pent(T);  % Pa

% simple Francis weir formula relationship---------------------------------

rho_avg = @(x, T) x.*rho_nbutane(T) + (1-x).*rho_npentane(T);
Mw_avg = @(x) x.*MM_but + (1-x).*MM_pent;
how = @(Mn, Mw_avg, rho_avg) (Mn .*Mw_avg ./A_available ./rho_avg /0.5);
F_L = @(h_ow) 1.84 .* Lw .* h_ow.^1.5;

% T-x data(from HYSYS)----------------------------------------------------

T_4_5 = [88.2113395690918; 86.9786095275229; 84.5761272486725; 83.4057796296778; ...
         82.2559599712353; 81.1269867341945; 80.0173763919771; 78.9271955687459; ...
         77.8564272215862; 76.8055113666752; 75.7726255654287; 74.7580773093107; ...
         73.7611358826471; 72.7820499120072; 71.8210319600412; 70.8766493453526; ...
         69.9483176235643; 69.0371454753988; 68.1410998525813; 66.3960039386178; ...
         65.5461335920334; 64.7108008312437; 63.8897324086249; 63.0826041625215; ...
         62.2890967145798; 61.5088957365871; 60.7416921719082; 59.9871824157029; ...
         59.2450684578368; 58.5150579921291; 57.7968644953284; 57.0902072789477; ...
         56.3945063410713; 55.7101030754815; 55.0364292038228; 54.3732274486648; ...
         53.0772400211675; 52.4439697754725; 51.8201979195605; 51.2055447145575; ...
         50.6000916437492; 50.0034677167495; 49.4154597166357; 48.8358595033393; ...
         48.2644639013294; 47.70107458432; 47.1454979592437; 46.5953002929688; ...
         ];

T_4_7 = [90.0101791381836; 88.7790821134614; 86.3786576006015; 85.2086886280696; ...
         84.0588414501664; 82.9294639289928; 81.8190811389165; 80.7277811997525; ...
         79.6555707222504; 78.6029207209017; 77.5680070837905; 76.5511631929888; ...
         75.5516633601869; 74.5697866016352; 73.6057720398301; 72.6581913574927; ...
         71.7264563428499; 70.8112140693964; 69.9119185813803; 68.1589266011244; ...
         67.3048733242145; 66.4652712601609; 65.639826964023; 64.8282257088655; ...
         64.0301570048566; 63.2453148875455; 62.4733981609749; 61.7141105996167; ...
         60.9671611128769; 60.2322638756794; 59.5091384284011; 58.7975097491973; ...
         58.0968031257565; 57.4073648837705; 56.7286313377941; 56.060349482323; ...
         54.7541561621479; 54.1157671872136; 53.4868704277988; 52.8670890157112; ...
         52.2565071193221; 51.6547562184202; 51.0616253636795; 50.4769084917288; ...
         49.9004043238882; 49.3319162613238; 48.7712522787111; 48.2177673339844; ...
         ];

T_4_9 = [91.7509399414063; 90.5229120489418; 89.3136516292892; 88.1244854964419; ...
         86.9548928631467; 85.8050324664813; 84.6752800873834; 83.5641685170795; ...
         81.3982220486274; 80.3439168104675; 79.3070642620596; 78.2880225121037; ...
         77.2860702675335; 76.3015157794924; 75.3346235311715; 74.3839694676092; ...
         73.4489613393127; 72.5302914019625; 71.627376142101; 70.7394240083692; ...
         69.8667731477453; 69.0089650654585; 68.1649620382712; 67.3352920550121; ...
         66.5193692246072; 64.927561406725; 64.1510849456819; 63.3871731748848; ...
         62.6355417137498; 61.8959110293044; 61.1680065554724; 60.451558784382; ...
         59.745998156652; 59.0516758078858; 58.3680325369969; 57.6948195141551; ...
         57.0317930936245; 56.3787147879971; 55.7353524981781; 55.1014749499055; ...
         54.4767080977051; 53.8611387527891; 52.6562856150915; 52.0665890987333; ...
         51.4851118779699; 50.9116590652536; 50.3460401906003; 49.7868286132813; ...
         ];

Xb = [0; 0.02; 0.06; 0.08; 0.1; 0.12; 0.14; 0.16; 0.18; 0.2; 0.22; 0.24; ...
      0.26; 0.28; 0.3; 0.32; 0.34; 0.36; 0.38; 0.42; 0.44; 0.46; 0.48; ...
      0.5; 0.52; 0.54; 0.56; 0.58; 0.6; 0.62; 0.64; 0.66; 0.68; 0.7; 0.72; ...
      0.74; 0.78; 0.8; 0.82; 0.84; 0.86; 0.88; 0.9; 0.92; 0.94; 0.96; 0.98; ...
      1];

T_x = [T_4_5, T_4_7, T_4_9];

createfigure(Xb, T_x)

% Linear model Poly3:  f(x) = p1*x^3 + p2*x^2 + p3*x + p4 -----------------
TXFit(Xb,T_4_7)
p1 = -6.18; 
p2 = 26.39; 
p3 = -62.01;  
p4 = 90;
T_mix = @(x) p1*x.^3 + p2*x.^2 + p3*x + p4 + 273.15; % C to K

%% Solving Steady State ===================================================

% V calculation -----------------------------------------------------------
% CONDITION OF Reboiler - Data on Saturation Curve - NIST: estimated by Saturation Properties for Pentane cause xp = 0.99
T_Rb = 367; % K
P_Rb = 4.9; % atm
Enthalpy_liq = 10.226e3;    % J/mol
Enthalpy_vap = 32.249e3;    % J/mol
Hfg = Enthalpy_vap - Enthalpy_liq;  % J/mol

V = R_duty/ Hfg;    % mol/s

% Mn & xb / x1 calculation ------------------------------------------------

x_guess = linspace(0.99,0.01,61)';
Mn_guess = linspace(18.305,47.78,61)';   % mol/s % 65.9 kmol/hr at top and 172 kmol/hr at bottom
init_guess = [Mn_guess, x_guess];

% Answers------------------------------------------------------------------
global Mn_final x_final

opt = optimoptions("fsolve",OptimalityTolerance=1e-3,MaxIterations=800, MaxFunctionEvaluations=1e8);
AnsP1 = fsolve(@stst_Equation,init_guess,opt);

x_final = AnsP1(:,2);
T_final = T_mix(AnsP1(:,2));
Mn_final = AnsP1(:,1);
y_final = Y(AnsP1(:,2));
Mw_n = Mw_avg(AnsP1(:,2));
rho_n = rho_avg(AnsP1(:,2), T_final);
how_final = how(Mn_final, Mw_n, rho_n);
Ln_final = F_L(how_final);


%% Cv valve calculation: --------------------------------------------------
% q = Cv* N* F(l) * sqrt(delta_P/ gs)  valve is linear
global q_stst N VO_stst gs_cond gs_Reb delta_stst_cond delta_stst_Reb Cv_cond Cv_Reb q_cond q_Reb

q_stst = 50 * (1000/3600); % mol/s (converted from Kmol/h)
N = 0.0865; % Correction factor
VO_stst = 0.5;   % Valve opening F(L)
gs_cond = 9.345814401727645e+03/ 54919;    % specific weight rho_substance/rho_water source: NIST
gs_Reb = 7.603376529595583e+03/ 53321;    % specific weight rho_substance/rho_water  source: NIST
delta_stst_cond = 4.5* 101325; % Pa (converted from 4.5 atm)
delta_stst_Reb = 4.9* 101325; % Pa (converted from 4.9 atm)

Cv_cond = q_stst/(N* VO_stst* sqrt(delta_stst_cond/ gs_cond));
Cv_Reb = q_stst/(N* VO_stst* sqrt(delta_stst_Reb/ gs_Reb));
Cv_cond = 0.3;
Cv_Reb = 0.2;

q_cond = @(dP) Cv_cond* N* VO_stst* sqrt(dP/ gs_cond);  % mol/s
q_Reb = @(dP) Cv_Reb* N* VO_stst* sqrt(dP/ gs_Reb); % mol/s

%% Cv valve calculation new: ----------------------------------------------
%  q = K* sqrt(delta_P)  valve is linear

g = 9.81;
K_cond = 0.2;
K_Reb = 0.2;
MQ_cond = @(M) K_cond * sqrt(M* g/ C_area);  % mol/s
MQ_Reb = @(M) K_Reb * sqrt(M* g/ R_area); % mol/s

%% Plotting stst Results Check with Rault's Law: --------------------------
P_final = P_mix(x_final, T_final);
P_real = 4.5e5*ones(1,61)';
for i = 2:60
    P_real(1) = C_P;
    P_real(i) = P_real(i-1) + pd;
    P_real(61) = R_P;
end
PMatrix = [P_final, P_real];

% plot parameters ---------------------------------------------------------
T_STSTfigure(T_final)
X_StStfigure(x_final)
P_ststfigure(PMatrix)


%% Dynamic Simiulation: Open Loop =========================================

Mn_initial = Mn_final;
MnX_initial = Mn_final.* x_final;
init_cond = [Mn_initial; MnX_initial]';
Tspan = [1 1e4];

[t, Mn_over_t] = ode45(@(t,Mn) Oploop(t,Mn), Tspan, init_cond);

Mn_OL = Mn_over_t(:,1:61);
MnX_OL = Mn_over_t(:,62:122);
X_OL = MnX_OL ./ Mn_OL;
T_OL = T_mix(X_OL);
P_OL = P_mix(X_OL, T_OL);
y_OL = Y(X_OL);
Mw_OL = Mw_avg(X_OL);
rho_OL = rho_avg(X_OL, T_OL);
how_OL = how(Mn_OL, Mw_OL, rho_OL);
Ln_OL = F_L(how_OL);

%plot(t,T_OL(:,1:10))
% plot(t,Mn_over_t(:,1:10))

% Dynamic Simiulation: Close Loop ========================================
% 2 P controllers for Distillate and Bottom products flow

% Initialize system
Mn_initial = Mn_final;
MnX_initial = Mn_final.* x_final;
init_cond = [Mn_initial; MnX_initial]';

%% Optimize controller gains
options = optimoptions('fmincon','Display','iter');
K0 = [3; 3];  % Initial guess for gains
lb = [5 ,5];      % Lower bounds
ub = [10 , 10];    % Upper bounds
[K_opt, J_opt] = fmincon(@objective, K0, [], [], [], [], lb, ub, [], options);

%% Run closed-loop simulation with optimized gains =====================
Tspan = [0 1e4];
[time, Mn_over_tP] = ode45(@(t,y) Closedloop(t, y, K_opt(1), K_opt(2)), Tspan, init_cond);

% Process results
Mn_CL = Mn_over_tP(:,1:61);
MnX_CL = Mn_over_tP(:,62:122);
X_CL = MnX_CL ./ Mn_CL;
T_CL = T_mix(X_CL);
P_CL = P_mix(X_CL, T_CL);
y_CL = Y(X_CL);
Mw_CL = Mw_avg(X_CL);
rho_CL = rho_avg(X_CL, T_CL);
how_CL = how(Mn_CL, Mw_CL, rho_CL);
Ln_CL = F_L(how_CL);

% Plot results
figure;
subplot(2,1,1);
plot(time, Mn_CL(:,1));
hold on;
plot(time, ones(size(time))*D_rate, 'r--');
title('Condenser Flow Rate');
legend('Actual', 'Setpoint');

subplot(2,1,2);
plot(time, Mn_CL(:,61));
hold on;
plot(time, ones(size(time))*B_rate, 'r--');
title('Reboiler Flow Rate');
legend('Actual', 'Setpoint');

%% Functions ==============================================================

function [fitresult, gof] = alphaFit(x, y)
%CREATEFIT(X,Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : x
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.

%% Fit: 'Equation of State fit'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( '(alpha*x)/(1+(alpha-1)*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = 0.736821294216622;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'Equation of State fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. x', 'Aplpha fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'x', 'Interpreter', 'none' );
ylabel( 'y', 'Interpreter', 'none' );
grid on
end

function createfigure(X1, YMatrix1)
%CREATEFIGURE(X1, YMatrix1)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(1),'DisplayName','4.5 atm');
set(plot1(2),'DisplayName','4.7 atm');
set(plot1(3),'DisplayName','4.9 atm');

% Create ylabel
ylabel({'T_b'});

% Create xlabel
xlabel({'x_B'});

% Create title
title({'T vs x_b of binary mixture'});

box(axes1,'on');
grid(axes1,'on');
hold(axes1,'off');
% Create legend
legend(axes1,'show');

end

function [fitresult, gof] = TXFit(Xb, T_4_7)

%CREATEFIT(XB,T_4_7)
%  Create a fit.
%
%  Data for 'T-x corelation' fit:
%      X Input : Xb
%      Y Output: T_4_7
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.

%% Fit: 'T-x correlation'.
[xData, yData] = prepareCurveData( Xb, T_4_7 );

% Set up fittype and options.
ft = fittype( 'poly3' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name', 'T-x corelation' );
h = plot( fitresult, xData, yData );
legend( h, 'T_4_7 vs. Xb', 'T-x corelation', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Xb', 'Interpreter', 'none' );
ylabel( 'T_4_7', 'Interpreter', 'none' );
grid on

end

function outputs = stst_Equation(guess)

global R MM_pent omega_pent Tc_pent Pc_pent Zc_pent Vc_pent
global MM_but omega_but Tc_but Pc_but Zc_but Vc_but
global F_rate F_temp F_zC4 F_zC5 D_rate D_zC4 D_zC5 B_rate B_zC4 B_zC5 RR
global Col_ID Col_trays Total_stages pd Lw A_available h_weir L_t2t Col_h
global C_P C_T C_duty R_P R_T R_duty
global alpha p1 p2 p3 p4
global Y rho_func rho_nbutane rho_npentane Cp_ig_but Cp_ig_pent
global Cp_liq_func Cp_liq_but Cp_liq_pent P_sat_but P_sat_pent
global T_sat_but T_sat_pent P_mix rho_avg Mw_avg how F_L T_mix V

Mn = guess(:,1);
x = guess(:,2);
yn = Y(x);
Tn = T_mix(x);
Mw_n = Mw_avg(x);
rho_n = rho_avg(x, Tn);
how_n = how(Mn, Mw_n, rho_n);
Ln = F_L(how_n);

% Condenser and Reflux Drum ----> S1
outputs(1,1) = V- (D_rate* RR) - D_rate;
outputs(1,2) = (V* yn(2)) - ((D_rate* RR) + D_rate)* D_zC4;

% Top Tray ----> S2
outputs(2,1) = Ln(1) - Ln(2);
outputs(2,2) = Ln(1)*x(1) - Ln(2)*x(2) + V*yn(3) - V*yn(2);

% nth Tray ----> S3-S22
for n = 3:22
    outputs(n,1) = Ln(n-1) - Ln(n);
    outputs(n,2) = Ln(n-1)*x(n-1) - Ln(n)*x(n) + V*yn(n+1) - V*yn(n);
end

% Feed Tray ----> S23
outputs(23,1) = Ln(22) - Ln(23) + F_rate;
outputs(23,2) = Ln(22)*x(22) - Ln(23)*x(23) + V*yn(24) - V*yn(23) + F_rate*F_zC4;

% nth Tray ----> S24-S59
for n = 24:59
    outputs(n,1) = Ln(n-1) - Ln(n);
    outputs(n,2) = Ln(n-1)*x(n-1) - Ln(n)*x(n) + V*yn(n+1) - V*yn(n);
end

% Last Tray ----> S60
outputs(60,1) = Ln(59) - Ln(60);
outputs(60,2) = Ln(59)*x(59) - Ln(60)*x(60) + V*yn(61) - V*yn(60);

% Reboiler and Column Base ----> S61
outputs(61,1) = Ln(60) - V - B_rate;
outputs(61,2) = Ln(60)*x(60) - V*yn(61) - B_rate*B_zC4;

end

function dANSdt = Oploop(t, unstst_eq)

global R MM_pent omega_pent Tc_pent Pc_pent Zc_pent Vc_pent
global MM_but omega_but Tc_but Pc_but Zc_but Vc_but
global F_rate F_temp F_zC4 F_zC5 D_rate D_zC4 D_zC5 B_rate B_zC4 B_zC5 RR
global Col_ID Col_trays Total_stages pd Lw A_available h_weir L_t2t Col_h
global C_P C_T C_duty R_P R_T R_duty
global alpha p1 p2 p3 p4
global Y rho_func rho_nbutane rho_npentane Cp_ig_but Cp_ig_pent
global Cp_liq_func Cp_liq_but Cp_liq_pent P_sat_but P_sat_pent
global T_sat_but T_sat_pent P_mix rho_avg Mw_avg how F_L T_mix V
global q_stst N VO_stst gs_cond gs_Reb delta_stst_cond delta_stst_Reb Cv_cond Cv_Reb q_cond q_Reb

dANSdt = zeros(122, 1);  % Initialize as a 122x1 column vector
Mn = ones(1,61);
MnX = ones(1,61);
x = ones(1,61);
for i = 1:61
Mn(i) = unstst_eq(i);
MnX(i) = unstst_eq(i+61);
x(i) = MnX(i) ./ Mn(i);
end


yn = Y(x);
Tn = T_mix(x);
Mw_n = Mw_avg(x);
rho_n = rho_avg(x, Tn);
how_n = how(Mn, Mw_n, rho_n);
Ln = F_L(how_n);

% Varient flow ------------------------------------------------------------
dP_cond = P_mix(x(1), Tn(1));    % Pa
dP_Reb = P_mix(x(61), Tn(61));   % Pa
Q_cond = q_cond(dP_cond);
Q_Reb = q_Reb(dP_Reb);

% Condenser and Reflux Drum ----> S1
dANSdt(1) = V- (Q_cond* RR) - Q_cond;
dANSdt(62) = (V* yn(2)) - ((Q_cond* RR) + Q_cond)* x(1);

% Top Tray ----> S2
dANSdt(2) = Ln(1) - Ln(2);
dANSdt(63) = Ln(1)*x(1) - Ln(2)*x(2) + V*yn(3) - V*yn(2);

% nth Tray ----> S3-S22
for n = 3:22
    dANSdt(n) = Ln(n-1) - Ln(n);
    dANSdt(n+61) = Ln(n-1)*x(n-1) - Ln(n)*x(n) + V*yn(n+1) - V*yn(n);
end

% Feed Tray ----> S23
dANSdt(23) = Ln(22) - Ln(23) + F_rate;
dANSdt(84) = Ln(22)*x(22) - Ln(23)*x(23) + V*yn(24) - V*yn(23) + F_rate*F_zC4;

% nth Tray ----> S24-S59
for n = 24:59
    dANSdt(n) = Ln(n-1) - Ln(n);
    dANSdt(n+61) = Ln(n-1)*x(n-1) - Ln(n)*x(n) + V*yn(n+1) - V*yn(n);
end

% Last Tray ----> S60
dANSdt(60) = Ln(59) - Ln(60);
dANSdt(121) = Ln(59)*x(59) - Ln(60)*x(60) + V*yn(61) - V*yn(60);

% Reboiler and Column Base ----> S61
dANSdt(61) = Ln(60) - V - Q_Reb;
dANSdt(122) = Ln(60)*x(60) - V*yn(61) - Q_Reb*x(61);


end

function T_STSTfigure(Y1)
%STSTFIGURE(Y1)
%  Y1:  vector of y data

% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot1 = plot(Y1,'DisplayName','T','LineWidth',1,'Color',[1 0 0]);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,1,321.628493015796);
% Create datatip
datatip(plot1,'DataIndex',1);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,2,326.005803593467);
% Create datatip
datatip(plot1,'DataIndex',2);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,16,326.700108101435);
% Create datatip
datatip(plot1,'DataIndex',16);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,23,325.49988237521);
% Create datatip
datatip(plot1,'DataIndex',23);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,24,346.045573696653);
% Create datatip
datatip(plot1,'DataIndex',24);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,38,348.820357546566);
% Create datatip
datatip(plot1,'DataIndex',38);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,51,356.611577249845);
% Create datatip
datatip(plot1,'DataIndex',51);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,61,362.698102785629);
% Create datatip
datatip(plot1,'DataIndex',61);

% Create ylabel
ylabel({'Temperature (K)'});

% Create xlabel
xlabel({'Tray number'});

% Create title
title({'Temperature distribution throughout the Tower'});

box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XGrid','on','YGrid','on');
% Create legend
legend(axes1,'show');

end

function X_StStfigure(Y1)
%CREATEFIGURE(Y1)
%  Y1:  vector of y data

% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create plot
plot1 = plot(Y1,'LineWidth',1,'LineStyle','--',...
    'Color',[0.466666666666667 0.674509803921569 0.188235294117647]);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,1,0.989999934268842);
% Create datatip
datatip(plot1,'DataIndex',1);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,2,0.840444516383016);
% Create datatip
datatip(plot1,'DataIndex',2);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,14,0.817775124631445);
% Create datatip
datatip(plot1,'DataIndex',14);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,23,0.85699420529391);
% Create datatip
datatip(plot1,'DataIndex',23);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,24,0.31492933346667);
% Create datatip
datatip(plot1,'DataIndex',24);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,36,0.275273238525078);
% Create datatip
datatip(plot1,'DataIndex',36);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,36,0.275273238525078);
% Create datatip
datatip(plot1,'DataIndex',36);

% The following line demonstrates an alternative way to create a data tip.
% datatip(plot1,61,0.00731019276267204);
% Create datatip
datatip(plot1,'DataIndex',61);

% Create ylabel
ylabel({'mole Fraction'});

% Create xlabel
xlabel({'Tray number'});

% Create title
title({'mole fraction per Tray'});

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.0470430107526842 70.0470430107527]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[-0.00132504737044381 0.998674952629556]);
box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'XGrid','on','XMinorGrid','on','YGrid','on','YMinorGrid','on');

end

function P_ststfigure(YMatrix1)
%CREATEFIGURE(YMatrix1)
%  YMATRIX1:  matrix of y data


% Create figure
figure1 = figure('WindowState','maximized');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix1);
set(plot1(1),'DisplayName','P calculated','LineWidth',1,...
    'Color',[0.494117647058824 0.184313725490196 0.556862745098039]);
set(plot1(2),'DisplayName','P real','LineStyle','--');

% Create ylabel
ylabel({'Pressure (Pa)'});

% Create xlabel
xlabel({'Tray number'});

% Create title
title({'Pressure over Tower'});

box(axes1,'on');
hold(axes1,'off');
% Create legend
legend(axes1,'show');

end

function dANSdt = Closedloop(t, unstst_eq, Kp_cond, Kp_reb)

global R MM_pent omega_pent Tc_pent Pc_pent Zc_pent Vc_pent
global MM_but omega_but Tc_but Pc_but Zc_but Vc_but
global F_rate F_temp F_zC4 F_zC5 D_rate D_zC4 D_zC5 B_rate B_zC4 B_zC5 RR
global Col_ID Col_trays Total_stages pd Lw A_available h_weir L_t2t Col_h
global C_P C_T C_duty R_P R_T R_duty
global alpha p1 p2 p3 p4
global Y rho_func rho_nbutane rho_npentane Cp_ig_but Cp_ig_pent
global Cp_liq_func Cp_liq_but Cp_liq_pent P_sat_but P_sat_pent
global T_sat_but T_sat_pent P_mix rho_avg Mw_avg how F_L T_mix V
global q_stst N VO_stst gs_cond gs_Reb delta_stst_cond delta_stst_Reb Cv_cond Cv_Reb q_cond q_Reb

dANSdt = zeros(122, 1);  % Initialize as a 122x1 column vector
Mn = ones(1,61);
MnX = ones(1,61);
x = ones(1,61);
for i = 1:61
Mn(i) = unstst_eq(i);
MnX(i) = unstst_eq(i+61);
x(i) = MnX(i) ./ Mn(i);
end


yn = Y(x);
Tn = T_mix(x);
Mw_n = Mw_avg(x);
rho_n = rho_avg(x, Tn);
how_n = how(Mn, Mw_n, rho_n);
Ln = F_L(how_n);

% Varient flow ------------------------------------------------------------
dP_cond = P_mix(x(1), Tn(1));    % Pa
dP_Reb = P_mix(x(61), Tn(61));   % Pa
Q_cond = q_cond(dP_cond);
Q_Reb = q_Reb(dP_Reb);

% P flow controllers ======================================================

% Calculate errors
error_cond = D_rate - Q_cond;
error_reb = B_rate - Q_Reb;
    
% Apply P control
Q_cond = Q_cond + Kp_cond * error_cond;
Q_Reb = Q_Reb + Kp_reb * error_reb;

% Condenser and Reflux Drum ----> S1
dANSdt(1) = V- (Q_cond* RR) - Q_cond;
dANSdt(62) = (V* yn(2)) - ((Q_cond* RR) + Q_cond)* x(1);

% Top Tray ----> S2
dANSdt(2) = Ln(1) - Ln(2);
dANSdt(63) = Ln(1)*x(1) - Ln(2)*x(2) + V*yn(3) - V*yn(2);

% nth Tray ----> S3-S22
for n = 3:22
    dANSdt(n) = Ln(n-1) - Ln(n);
    dANSdt(n+61) = Ln(n-1)*x(n-1) - Ln(n)*x(n) + V*yn(n+1) - V*yn(n);
end

% Feed Tray ----> S23
dANSdt(23) = Ln(22) - Ln(23) + F_rate;
dANSdt(84) = Ln(22)*x(22) - Ln(23)*x(23) + V*yn(24) - V*yn(23) + F_rate*F_zC4;

% nth Tray ----> S24-S59
for n = 24:59
    dANSdt(n) = Ln(n-1) - Ln(n);
    dANSdt(n+61) = Ln(n-1)*x(n-1) - Ln(n)*x(n) + V*yn(n+1) - V*yn(n);
end

% Last Tray ----> S60
dANSdt(60) = Ln(59) - Ln(60);
dANSdt(121) = Ln(59)*x(59) - Ln(60)*x(60) + V*yn(61) - V*yn(60);

% Reboiler and Column Base ----> S61
dANSdt(61) = Ln(60) - V - Q_Reb;
dANSdt(122) = Ln(60)*x(60) - V*yn(61) - Q_Reb*x(61);

end

function J = objective(K) 
% Function to optimize controller gains====================================
    global Mn_final x_final D_rate B_rate P_mix T_mix q_cond q_Reb Y
    
    Kp_cond = K(1);
    Kp_reb = K(2);
    Mn_initial = Mn_final;
    MnX_initial = Mn_final.* x_final;
    init_cond = [Mn_initial; MnX_initial]';
    
    % Simulate the system with these gains
    [t, y] = ode45(@(t,y) Closedloop(t, y, Kp_cond, Kp_reb), [0 10000], init_cond);
    
    % Results
    Mn_cl = y(:,1:61);
    MnX_cl = y(:,62:122);
    X_cl = MnX_cl ./ Mn_cl;
    T_cl = T_mix(X_cl);

    % Calculate error
    error_cond = D_rate - q_cond(P_mix(X_cl(end,1) , T_cl(end, 1)));
    error_reb = B_rate - q_Reb(P_mix(X_cl(end,61) , T_cl(end, 61)));
    
    % Objective function (sum of squared errors)
    J = error_cond^2 + error_reb^2;

end