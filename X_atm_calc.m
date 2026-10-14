clear all; clc;

% =========================================================================
% X_atm_calc.m
% 
% Calculates modified specific impact energy (Q_R_prime) and total
% atmospheric loss fraction (X_atm) for planetary collisions.
% Based on the scaling laws developed in:
%   Roche et al. (2026, ApJL) - "Giant impacts preferentially remove
%   low-mass atmosphere: a generalised scaling law for impact-driven
%   atmospheric loss"
%
% USAGE:
%   1. Edit the input parameter array `params` below (all in SI units),
%      OR define `params` in your MATLAB workspace before running:
%        params = [R_t, R_i, M_t, M_i, b, v_c, f_atm];
%   2. Run this script:
%        X_atm_calc
%
% INPUT PARAMETER ARRAY (params) - ALL IN SI SCALE:
%   params(1) : R_t   - Target planet radius [m]
%   params(2) : R_i   - Impactor planet radius [m]
%   params(3) : M_t   - Target planet mass [kg]
%   params(4) : M_i   - Impactor planet mass [kg]
%   params(5) : b     - Impact parameter [dimensionless, 0 <= b <= 1]
%   params(6) : v_c   - Impact velocity [m/s]
%   params(7) : f_atm - Target atmospheric mass fraction [dimensionless, M_atm / M_t]
% =========================================================================

% -------------------------------------------------------------------------
% 0. Input Parameter Array (SI Scale)
% -------------------------------------------------------------------------
if ~exist('params', 'var') || isempty(params)
    % Default example parameters (SI units):
    % Target: Venus-like (R_t = 6.0518e6 m, M_t = 4.86884e24 kg)
    % Impactor: R_i = 3.5e6 m, M_i = 39.584e9 kg, b = 0.7071, v_c = 20 km/s
    % Atmosphere mass fraction: f_atm = 1.0e-4
    params = [ ...
        7.5325361e6, ...          % 1. R_t   [m]
        3.342124e6, ...             % 2. R_i   [m]
        6.0154152e24, ...        % 3. M_t   [kg]
        6.630696e23, ...    % 4. M_i   [kg]
        0, ...     % 5. b     [-]
        9.53e3, ...            % 6. v_c   [m/s]
        0.01 ...             % 7. f_atm [-]
    ];
end

% Unpack parameters:
R_t   = params(1);
R_i   = params(2);
M_t   = params(3);
M_i   = params(4);
b     = params(5);
v_c   = params(6);
f_atm = params(7);

% -------------------------------------------------------------------------
% 1. Calculate Q_R_prime (Modified Specific Impact Energy)
% -------------------------------------------------------------------------
B = (R_t + R_i) * b;

% Projected length of the projectile overlapping the target:
if (B + R_i) <= R_t
    l = 2 * R_i;
    alpha = 1.0;
else
    l = R_t + R_i - B;
    alpha = (3 * R_i * l^2 - l^3) / (4 * R_i^3);
end

M_tot = M_t + M_i;
mu = (M_t * M_i) / M_tot;
mu_alpha = (alpha * M_t * M_i) / ((alpha * M_i) + M_t);

Q_R = (mu * v_c^2) / (2 * M_tot);       % Unmodified specific impact energy [J/kg]
Q_R_prime = (mu_alpha * Q_R) / mu;      % Modified specific impact energy [J/kg]
Q_R_prime_MJ = Q_R_prime / 1.0e6;       % Converted to [MJ/kg] for X_atm calculation

% -------------------------------------------------------------------------
% 2. Physical Constants & Derived Dimensionless Variables
% -------------------------------------------------------------------------
G       = 6.67408e-11;   % Gravitational constant [m^3 kg^-1 s^-2]
M_earth = 5.9724e24;     % Earth mass [kg] (Roche et al. standard)

% Masses in Earth masses:
M_t_earth = M_t / M_earth;
M_i_earth = M_i / M_earth;

% Impactor-total mass ratio:
gamma = M_i / (M_t + M_i);

% Impactor to target radius ratio:
R_ratio = R_i / R_t;

% Mutual escape velocity [m/s]:
v_esc = sqrt(2 * G * (M_t + M_i) / (R_t + R_i));

% Impact velocity in units of mutual escape velocity:
vc = v_c / v_esc;

% -------------------------------------------------------------------------
% 3. Empirical Fitting Parameters (Roche et al. 2026 ApJL)
% -------------------------------------------------------------------------
% Near-field mass fraction parameters:
q11 = -0.0011404070287601;
q12 =  0.0005853554706874;
q13 =  0.0004229310650559;
q14 =  2.0;
q15 = -0.0001681893595867;
q21 = -2.990394733907284;
q22 = -4.340509285222407;
q31 =  3.680308721661617;
q32 = -3.9046955628287696;
q33 =  1.881415631996061;
q34 =  2.0;
q35 =  2.328882600831519;
q36 = -6.797028501938182;
q38 = -0.0411818029322558;
q41 =  161.05720316915694;
q42 = -7.690384546730376;
q43 =  7.577016664697488;
q44 =  1.0308432839505288;
q45 =  1.033504078552506;
q46 = -2.3556779061398494;
q47 = -0.0012879570817664;
zeta5 = -161.02157740225127;
zeta6 =  8.28566684431704e-08;

% Near-field loss parameters:
k11 = -535.6899898220761;
k12 =  843.4184481359747;
k13 = -847.0795889608318;
k14 =  494.21158521587927;
k15 =  0.0002703723957616;
k16 = -0.0187233083484831;
k17 = -0.3845105615920718;
k21 =  0.0029091121113151;
k22 =  0.0013483201004243;
k23 = -0.0013547525798418;
k24 = -0.0029748236386274;
k25 =  0.0002703723957616;
k31 =  4811.153182686202;
k32 = -0.0364090718700348;
k33 = -0.2382117497656167;
k34 = -4020.819856042687;
k35 =  0.0002703723957616;

% Far-field loss parameters:
s11 =  1020.3918108627572;
s12 =  9.294519886102522;
s13 = -0.11166497483473946;
s14 = -1030.1005018606777;
s15 =  3.262206184294514e-05;
s16 = -0.0501934230566666;
s21 = -3.8940780439779576;
s22 = -1.1017541187112085;
s23 = -0.0004728717375858641;
s24 =  4.993174134062099;
s25 =  0.000581419775592155;
s26 =  0.003429221121849461;
s31 = -3208.3991776498683;
s32 =  3208.848548705081;
s33 =  0.00038822588810121883;
s41 = -322.486418835046;
s42 = -1.9230068604833475;
s43 =  1.415457700191505;
s44 =  322.9055347383739;
s45 =  0.0001799965131025202;
s46 = -0.1320268872697112;

% -------------------------------------------------------------------------
% 4. Near-Field & Far-Field Atmosphere Mass Fractions (f_NF, f_FF)
% -------------------------------------------------------------------------
zeta1 = q11 + (q12 * b) + (q13 * b^q14) + (q15 * f_atm);
zeta2 = q21 + (q22 * b);
zeta3 = q31 + (q32 * b) + (q33 * b^q34) + (q35 * f_atm) + (q36 * f_atm^2) + (M_t_earth)^q38;
zeta4 = q41 + (q42 * b) + (q43 * b^q44) + (q45 * f_atm) + (q46 * f_atm^2) + (q47 * M_t_earth);

f_NF = zeta4 / ((1.0 + zeta6 * exp(zeta3 * (R_ratio - zeta2)))^zeta1) + zeta5;

% Clip f_NF to [0, 1] using simple if statements:
if f_NF < 0.0
    f_NF = 0.0;
elseif f_NF > 1.0
    f_NF = 1.0;
end

f_FF = 1.0 - f_NF;

% -------------------------------------------------------------------------
% 5. Near-Field Atmospheric Loss Fraction (X_NF)
% -------------------------------------------------------------------------
log10_fatm = log10(f_atm);

xi1 = k11 + (k12 * gamma) + (k13 * (gamma + 0.05)^2) + (k14 * vc^k15) + (k16 * (M_t_earth + 1.0)) + (k17 * log10_fatm);
xi2 = k21 + (k22 * gamma) + (k23 * (gamma + 0.05)^2) + (k24 * vc^k25);
xi3 = k31 + (k32 * gamma) + (k33 * (gamma + 0.05)^2) + (k34 * vc^k35);

X_NF = f_NF * (xi1 - xi2 * (b + xi3)^2);

% Near-field loss evaluated at mutual escape velocity (vc = 1.0):
xi1_vesc = k11 + (k12 * gamma) + (k13 * (gamma + 0.05)^2) + k14 + (k16 * (M_t_earth + 1.0)) + (k17 * log10_fatm);
xi2_vesc = k21 + (k22 * gamma) + (k23 * (gamma + 0.05)^2) + k24;
xi3_vesc = k31 + (k32 * gamma) + (k33 * (gamma + 0.05)^2) + k34;

X_NF_vesc = f_NF * (xi1_vesc - xi2_vesc * (b + xi3_vesc)^2);

% Clip X_NF to [max(0, X_NF_vesc), f_NF] using simple if statements:
lower_bound_NF = max(0.0, X_NF_vesc);
if X_NF < lower_bound_NF
    X_NF = lower_bound_NF;
elseif X_NF > f_NF
    X_NF = f_NF;
end

% -------------------------------------------------------------------------
% 6. Far-Field Atmospheric Loss Fraction (X_FF)
% -------------------------------------------------------------------------
psi1 = s11 + (s12 * (gamma + 0.05)^s13) + (s14 * (M_t_earth + 0.05)^s15) + (s16 * log10_fatm);
psi2 = s21 + (s22 * (gamma + 0.05)^s23) + (s24 * (M_t_earth + 0.05)^s25) + (s26 * log10_fatm);
psi3 = s31 + (s32 * (gamma + 0.05)^s33);
psi4 = s41 + (s42 * (gamma + 0.05)^s43) + (s44 * (M_t_earth + 0.05)^s45) + (s46 * log10_fatm);

X_FF = f_FF * (psi1 * exp(-psi2 * (Q_R_prime_MJ * (1.0 + (M_i_earth / M_t_earth)) * (1.0 - b)^psi4)) + psi3);

% Clip X_FF to [0, f_FF] using simple if statements:
if X_FF < 0.0
    X_FF = 0.0;
elseif X_FF > f_FF
    X_FF = f_FF;
end

% -------------------------------------------------------------------------
% 7. Total Atmospheric Loss Fraction (X_atm)
% -------------------------------------------------------------------------
X_atm = X_NF + X_FF;

% Clip X_atm to [0, 1] using simple if statements:
if X_atm < 0.0
    X_atm = 0.0;
elseif X_atm > 1.0
    X_atm = 1.0;
end

% -------------------------------------------------------------------------
% 8. Display Results
% -------------------------------------------------------------------------
fprintf('\n=================================================================\n');
fprintf('                ATMOSPHERIC LOSS CALCULATION RESULTS             \n');
fprintf('=================================================================\n');
fprintf('INPUTS (SI scale):\n');
fprintf('  Target Radius (R_t)        : %.4e m\n', R_t);
fprintf('  Impactor Radius (R_i)      : %.4e m\n', R_i);
fprintf('  Target Mass (M_t)          : %.4e kg (%.4f M_Earth)\n', M_t, M_t_earth);
fprintf('  Impactor Mass (M_i)        : %.4e kg (%.4e M_Earth)\n', M_i, M_i_earth);
fprintf('  Impact Parameter (b)       : %.6f\n', b);
fprintf('  Impact Velocity (v_c)      : %.4e m/s (%.2f km/s)\n', v_c, v_c / 1e3);
fprintf('  Atmospheric Fraction (f_atm): %.4e\n', f_atm);
fprintf('-----------------------------------------------------------------\n');
fprintf('INTERMEDIATE DERIVED QUANTITIES:\n');
fprintf('  Mutual Escape Velocity (v_esc): %.4e m/s (%.2f km/s)\n', v_esc, v_esc / 1e3);
fprintf('  Velocity Ratio (v_c / v_esc)  : %.4f\n', vc);
fprintf('  Mass Ratio gamma (M_i / M_tot): %.4e\n', gamma);
fprintf('  Radius Ratio (R_i / R_t)      : %.4f\n', R_ratio);
fprintf('  Projected Overlap (alpha)     : %.6f\n', alpha);
fprintf('  Specific Energy Q_R_prime     : %.6e J/kg (%.6e MJ/kg)\n', Q_R_prime, Q_R_prime_MJ);
fprintf('-----------------------------------------------------------------\n');
fprintf('ATMOSPHERIC FRACTIONS:\n');
fprintf('  Near-Field Mass Fraction (f_NF): %.6f\n', f_NF);
fprintf('  Far-Field Mass Fraction (f_FF) : %.6f\n', f_FF);
fprintf('LOSS FRACTIONS:\n');
fprintf('  Near-Field Loss (X_NF)        : %.6f (%.2f%%)\n', X_NF, X_NF * 100);
fprintf('  Far-Field Loss (X_FF)         : %.6f (%.2f%%)\n', X_FF, X_FF * 100);
fprintf('  TOTAL ATMOSPHERIC LOSS (X_atm): %.6f (%.2f%%)\n', X_atm, X_atm * 100);
fprintf('=================================================================\n\n');
