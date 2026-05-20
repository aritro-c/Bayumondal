% Atmosphere loss parameterizations by
% - V. Shuvalov, Meteorit. Planet. Sci., 44, 1095-1105 (2009) and V. Shuvalov, LPSC, abstract no. 1191 (2010).
% - H. J. Melosh & A. M. Vickery, Nature, 338, 487-489 (1989) and A. M. Vickery, LPI Contributions, 746, 51-52 (1990).
% - S. T. Stewart & S. Lock, unpublished manuscript
% - J. A. Kegerreis et al., Astrophys. J. Lett. 901, L31 (2020)
%
% (last updated: 18/10/2024)
%Including various additions, such as airburst, kegerreis, volatile
%deposition, composite escapes.

% First clean up before starting the script
clear all; clc;


% Read data from N-Body impact file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
data           = readtable('VENUS.coll.250.dat', 'VariableNamingRule', 'preserve');

% Load continuous degassing data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degas_data     = readtable('tseries_outgassed_water.txt', 'VariableNamingRule', 'preserve');
% Scale normalized outgassed water to kg
M_degas_H2O_total_kg = degas_data.outgassed_water * 2.27e20;


% Choose planet: 0 - Earth, 1 - Venus, 2 - Mars
planet         = 1;

% Time shift configuration [Myr]
% (Impacts start at time zero in N-body data; this shifts them to match planetary evolution time)
time_shift_Myr = 100; 

% Plotting Limits [Myr]
% (Choose the X-axis range for the plots; simulation will run to cover this range)
plot_xlim      = [0, 350]; 

% Apply shift to impact data
data.time      = data.time + time_shift_Myr * 1e6;

% Timestep configuration [years]
% Set use_auto_dt to true for automatic selection, or false to use manual_dt
use_auto_dt    = true; 
manual_dt      = 1e3; 

if use_auto_dt
    % Calculate the smallest interval in the degassing data
    t_diffs = diff(degas_data.Time_Years);
    min_interval = min(t_diffs(t_diffs > 0)); % Ignore zeros
    
    if isempty(min_interval)
        dt = manual_dt; % Fallback if data is insufficient
        fprintf('Auto-dt: Insufficient data, using manual_dt = %d yr\n', dt);
    else
        % Round down to nearest multiple of 50
        dt = floor(min_interval / 50) * 50;
        
        % Ensure dt is at least 50 (or some reasonable floor)
        if dt < 50; dt = 50; end
        fprintf('Auto-dt: Smallest interval is %.2f yr. Setting dt = %d yr\n', min_interval, dt);
    end
else
    dt = manual_dt;
    fprintf('Manual-dt: Using dt = %d yr\n', dt);
end

t_start        = 0;
t_end          = max(max(data.time), plot_xlim(2) * 1e6);
t_vec          = t_start:dt:t_end;
num_steps      = length(t_vec);

% Interpolate cumulative degassing onto simulation time vector
[t_unique, idx] = unique(degas_data.Time_Years);
M_degas_H2O_unique = M_degas_H2O_total_kg(idx);
M_degas_H2O_interp = interp1(t_unique, M_degas_H2O_unique, t_vec, 'linear', 'extrap');

% Impactor composition: Percentage of CC (Carbonaceous Chondrite), balance is EC (Enstatite Chondrite)
cc_percent     =  2.5;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Counter value
counter        = 0;



% Define radius of interest [m]
r_int          = 3.000e5;

% Define constants
pival          = 3.141592654; % Pi [non-dim.]
G              = 6.67384e-11; % Gravitational constant [m^/(kg s^2)]
H_vap          = 13.000e6;    % Vaporization energy of silicates [J/kg]
M_sun          = 1.98847e30;  % Solar mass [kg] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Export Config

    % Figure size (inches)
    % fig_width  = 10; 
    % fig_height = 14;
    fig_width  = 14; 
    fig_height = 9;
    
    % Global Plot Style
    plot_lw    = 5;      % LineWidth for plots
    border_lw  = 1.5;    % LineWidth for axes border
    legend_fs  = 18;     % FontSize for legends
    axes_fs    = 24;     % FontSize for axes labels
    
    % Export Settings
    export_format = 'png'; % 'svg' or 'png'
    export_dpi    = 300;   % DPI for PNG export
    
    if strcmp(export_format, 'svg')
        export_ext = '.svg';
        export_driver = '-dsvg';
        export_res = '-painters'; % Use painters for SVG
        bg_color = 'none';
    else
        export_ext = '.png';
        export_driver = '-dpng';
        export_res = ['-r', num2str(export_dpi)];
        bg_color = 'w'; % White background for PNG
    end
    
    
    % Volatile delivery parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    depo_eff      = 0.36;        % Deposition efficiency
    
    % Endmember compositions (Weight fractions)
    % EC (Enstatite Chondrite)
    ec_H2O = 0.001; ec_CO2 = 0.004; ec_N2 = 0.0002;
    % CC (Carbonaceous Chondrite)
    cc_H2O = 0.08;  cc_CO2 = 0.04;  cc_N2 = 0.002;
    
    % Mix composition
    f_cc  = cc_percent / 100;
    f_ec  = 1 - f_cc;
    
    f_H2O = f_cc * cc_H2O + f_ec * ec_H2O;
    f_CO2 = f_cc * cc_CO2 + f_ec * ec_CO2;
    f_N2  = f_cc * cc_N2  + f_ec * ec_N2;
    
    % Composition label for plots
    if cc_percent == 0
        comp_label = 'EC';
    elseif cc_percent == 100
        comp_label = 'CC';
    else
        comp_label = [num2str(cc_percent), '% CC'];
    end
%%  
    
% Fit parameters for Stewart & Lock model [non-dim.]
a0             = -3.403;
a1             =  0.562;

%factor to vary composition for Kegerrreis law.
factorcomp=   1.00; %0.3 %1.00


% Physical parameters of Earth as target body
if planet==0
    M_p        = 5.9736e24;   % Earth mass [kg]
    R_p        = 6.3781e6;    % Earth radius [m]
    rho_mean_p = 5515.000;    % Mean density of planet [kg/m^3]
    rho_p      = 2900.000;    % Density of crust [kg/m^3] (assumes that basaltic crust is target material)
    M_at       = 5.100e18;    % Earth atmosphere mass [kg]
    rho_atm    = 1.217;       % Surface density of Earth's atmosphere [kg/m^3]
    height_atm = 8500.0;      % Scale height of Earth atmosphere [m]
    mu         = 29;          % molecular mass of (N2 + O2) atmosphere [non-dim.] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M_N2  = M_at * 0.7808; %%%%%%%%%%%%%%%%%
    M_H2O  = M_at * 0.01; %%%%%%%%%%%%%%%%%
    M_CO2 = M_at * 0.0004; %%%%%%%%%%%%%%%%%%


end


% Physical parameters of Venus as target body
if planet==1
    M_p        = 4.8685e24;   % Venus mass [kg]
    R_p        = 6.0518e6;    % Venus radius [m]
    rho_mean_p = 5204.000;    % Mean density of planet [kg/m^3]
    rho_p      = 2900.000;    % Density of crust [kg/m^3] (assumes that basaltic crust is target material)
    M_at       = 3.4e20;      % Venus atmosphere mass [kg]
    rho_atm    = 64.790;      % 50.000;       %64.790;      % Surface density of Venus' atmosphere [kg/m^3]
    height_atm = 15900.0;     % 18900.0;       %15900.0;     % Scale height of Venus atmosphere [m]
    mu         = 44;          % molecular mass of pure CO2 atmosphere [non-dim.] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M_N2  = M_at * 0.0000; %%%%%%%%%%%%%%%%%%%%
    M_H2O  = M_at * 0.0000; %%%%%%%%%%%%%%%%%%%%%
    M_CO2 = M_at * 1.0000; %%%%%%%%%%%%%%%%%%%


end

% Physical parameters of Mars as target body
if planet==2
    M_p        = 0.64171e24;   % Mars mass [kg]
    R_p        = 3.3962e6;    % Mars radius [m]
    rho_mean_p = 3933.000;    % Mean density of planet [kg/m^3]
    rho_p      = 2900.000;    % Density of crust [kg/m^3] (assumes that crust is target material)
    M_at       = 2.5e18;           %2.5e16     % Mars atmosphere mass [kg]
    rho_atm    = 1.2;           %0.020;      % Surface density of Mars' atmosphere [kg/m^3]
    height_atm = 14100.0;      %11100.0;     % Scale height of Mars atmosphere [m]
    mu         = 44;          % molecular mass of pure CO2 atmosphere [non-dim.] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M_N2  = M_at * 0.027; %%%%%%%%%%%%%%%%%%%%
    M_H2O  = M_at * 0.000; %%%%%%%%%%%%%%%%%%%%%
    M_CO2 = M_at * 0.9500; %%%%%%%%%%%%%%%%%%%


end


% Variables to calculate remaining atmosphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_at_rem_shu   = M_at;  % for Shuvalov    
M_at_rem_ke    = M_at;  % for Kegerreis 
M_at_rem_comp  = M_at;  % for composite

% Shuvalov Reservoir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_CO2_shu = M_CO2; 
M_N2_shu  = M_N2;  
M_H2O_shu = M_H2O; 

% Kegerreis Reservoir %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_CO2_ke = M_CO2; 
M_N2_ke  = M_N2;  
M_H2O_ke = M_H2O;

% Composite Reservoir %%%%%%%%%%%%%%%%%%%%%%%%%%%
M_CO2_comp     = M_CO2; 
M_N2_comp      = M_N2;  
M_H2O_comp     = M_H2O;

% Impactor density [kg/m^3]
rho_imp        = 2700.000;    

% Compute total number of radius steps
i_max          = height(data);  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define vector sizes
atm_loss       = zeros(i_max,1);
chi            = zeros(i_max,1);
d              = zeros(i_max,1);
Dchi           = zeros(i_max,1);
f_M            = zeros(i_max,1);
loss_frac      = zeros(i_max,1);
loss_frac2     = zeros(i_max,1);
loss_frac_mv   = zeros(i_max,1);
loss_frac2_mv  = zeros(i_max,1);
M_atm          = M_at*ones(i_max,1);
M_imp          = zeros(i_max,1);
M_loss         = zeros(i_max,1);
M_loss_ke      = zeros(i_max,1);
M_loss_mv      = zeros(i_max,1);
M_loss_Fe      = zeros(i_max,1);
M_deliv        = zeros(i_max,1);
MD             = zeros(i_max,1);
M_tot          = zeros(i_max,1);
v_esc          = zeros(i_max,1);
vel_ratio      = zeros(i_max,1);
V_cap_imp      = zeros(i_max,1);
V_cap_tar      = zeros(i_max,1);
V_imp          = zeros(i_max,1);
X              = zeros(i_max,1);
xi             = zeros(i_max,1);

R              = zeros(i_max,1);
v              = zeros(i_max,1);

M_loss_comp    = zeros(i_max,1);
M_loss_AB      = zeros(i_max,1);
Hhalf          = zeros(i_max,1);
Hfr            = zeros(i_max,1);
M_loss_S_frag_AB = zeros(i_max,1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M_imp_all      = zeros(i_max,1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vectors to store the history of the atmosphere %%%%%%%%%%%%%%%%%%%%%%%%

% Shuvalov
M_at_history_shu = zeros(num_steps, 1); 
    
M_CO2_history_shu = zeros(num_steps, 1);
M_N2_history_shu  = zeros(num_steps, 1);
M_H2O_history_shu  = zeros(num_steps, 1);
    
% Kegerreis
M_at_history_ke  = zeros(num_steps, 1); 
    
M_CO2_history_ke = zeros(num_steps, 1);
M_N2_history_ke  = zeros(num_steps, 1);
M_H2O_history_ke  = zeros(num_steps, 1);

% Composite
M_at_history_comp  = zeros(num_steps, 1); 

M_CO2_history_comp = zeros(num_steps, 1);
M_N2_history_comp  = zeros(num_steps, 1);
M_H2O_history_comp = zeros(num_steps, 1);


% Vectors to track species delivery [kg]%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_H2O_delivered = zeros(num_steps, 1);
M_CO2_delivered = zeros(num_steps, 1);
M_N2_delivered  = zeros(num_steps, 1);

% Track scale height of the atmosphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
height_atm_history = zeros(num_steps, 1);

% Track molecular mass of the atmosphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_history = zeros(num_steps, 1);

% Initial tracking for delivery
m_h2o_delivered_cum = 0;
m_co2_delivered_cum = 0;
m_n2_delivered_cum  = 0;

impact_idx = 1;

% Progress bar initialization
prog_bar_width = 50;
fprintf('Progress: [%s]', repmat(' ', 1, prog_bar_width));
next_prog_update = 0;

% Loop over discrete timesteps
for j=1:1:num_steps
    
    % Update progress bar
    if j >= ((next_prog_update + 1) / prog_bar_width) * num_steps
        % Delete the current bar content and closing bracket
        fprintf(repmat('\b', 1, prog_bar_width + 1)); 
        
        next_prog_update = next_prog_update + 1;
        
        % Reprint updated bar and closing bracket
        fprintf('%s%s]', repmat('#', 1, next_prog_update), repmat(' ', 1, prog_bar_width - next_prog_update));
    end
    
    t_curr = t_vec(j);
    t_next = t_curr + dt;
    
    % Add continuous degassing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if j == 1
        dm_degas = M_degas_H2O_interp(j);
    else
        dm_degas = M_degas_H2O_interp(j) - M_degas_H2O_interp(j-1);
    end
    
    % Update reservoirs
    M_H2O_shu   = M_H2O_shu + dm_degas;
    M_H2O_ke    = M_H2O_ke + dm_degas;
    M_H2O_comp  = M_H2O_comp + dm_degas;
    
    % Update total atmospheric mass to reflect new outgassed water
    M_at_rem_shu   = M_CO2_shu + M_N2_shu + M_H2O_shu;
    M_at_rem_ke    = M_CO2_ke + M_N2_ke + M_H2O_ke;
    M_at_rem_comp  = M_CO2_comp + M_N2_comp + M_H2O_comp;
    
    % Process all impacts in this timestep
    while impact_idx <= i_max && data.time(impact_idx) < t_next
        i = impact_idx;
        
        % A: Shuvalov (2009, 2010) model:
        % ===============================

         % Current mass of impactor [kg]
        M_imp(i)        = data.mimp(i) * M_sun; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        M_imp_all(i)    = M_imp(i); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        % Calculate Volatile Delivery for impact step %%%%%%%%%%%%%%%%%%%
        m_h2o_add = M_imp(i) * f_H2O * depo_eff;
        m_co2_add = M_imp(i) * f_CO2 * depo_eff;
        m_n2_add  = M_imp(i) * f_N2 * depo_eff;
        M_delivered_total = m_h2o_add + m_co2_add + m_n2_add;


        % Current radius of impactor (calculated from mass) [m]%%%%%%%%%%%%%%%
        R(i)        = (M_imp(i) / (4/3 * pival * rho_imp))^(1/3);


        % Mutual escape velocity [m/s]
        v_esc(i)    = sqrt(2*G*(M_p+M_imp(i))/(R_p+R(i)));

        % Define impact velocity of interest [m/s] %%%%%%%%%%%%%%%%%%%%%%%%%%%
        v(i)          = data.("vrel/vesc")(i) * v_esc(i);


        % Current atm: %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        M_at_current = M_CO2_shu + M_N2_shu + M_H2O_shu;
        
        % Molecular weight (mu) %%%%%%%%%%%%%%%%%%%%%%%%%
        moles_total = (M_H2O_shu/18) + (M_N2_shu/28) + (M_CO2_shu/44);
        mu_current_val = M_at_current / moles_total;
        
        % Scale Height and Surface Density %%%%%%%%%%%%%%%%%%%%%%
        height_atm_current_val = height_atm * (mu / mu_current_val);
        rho_atm_current = (M_at_current / M_at) * rho_atm;
        

        % Impact efficiency [non-dim.] (Shuvalov's eq. (2)) %%%%%%%%%%%%%%%%%%
        xi(i)       = (2.000*R(i))^3*rho_imp*rho_p*(v(i)^2-v_esc(i)^2)/(height_atm_current_val^3*rho_atm_current*v_esc(i)^2*(rho_p+rho_imp));
        
        xi(i)=max(xi(i),1000);
        
        % Slope of factor xi
        Dchi(i)     = (9.70992e-7)*(xi(i)^(4.239))*(0.000100423*(log(xi(i))^4)-0.00526501*(log(xi(i))^3)+0.0975584*(log(xi(i))^2)-0.800091*log(xi(i))+2.27527)*10^((log(xi(i))^2)*(0.0000200847*(log(xi(i))^3)-0.00131625*(log(xi(i))^2)+0.0325195*log(xi(i))-0.400045));
        
        
        % Polynomial fit [non-dim.] (reformulated version of Shuvalov's equation (5))
        chi(i)   = 10^(-6.375+5.239*log10(xi(i))-2.121*(log10(xi(i))^2)+0.397*(log10(xi(i))^3)-0.037*(log10(xi(i))^4)+0.0013*(log10(xi(i))^5));  
            
        % Velocity ratio [non-dim.]
        vel_ratio(i)  = ((v(i)^2-v_esc(i)^2)/v_esc(i)^2);
        
        % Atmospheric mass loss [kg] (reformulated version of Shuvalov's equation (4))
        M_loss(i)     = chi(i)*M_imp(i)*((v(i)^2-v_esc(i)^2)/v_esc(i)^2);
      


    %%
        
        
        % for AIRBURSTS for Shuvalov 2014.

        Hhalf(i)=-1.4*height_atm_current_val* log(0.5*(1/height_atm_current_val)*5.5*(2.000*R(i))*(rho_imp/rho_atm_current)^0.5);

        Hfr(i)=-1.4*height_atm_current_val* log(0.5*(1/height_atm_current_val)*3.5*(2.000*R(i))*(rho_imp/rho_atm_current)^0.5);

        Vlim= (15+0.2*height_atm_current_val/1000)*1000; % Convert km/s to m/s %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        M_loss_AB(i)= M_imp(i)*5.5*0.001*(2.000*R(i)/height_atm_current_val)*((rho_imp/1000.)^0.5)*((v(i)-Vlim)^2.);
        if(Hhalf(i)<0)
            M_loss_S_frag_AB(i)=M_loss(i);
        elseif(Hfr(i)>0)
            M_loss_S_frag_AB(i)=M_loss_AB(i);
        else
           M_loss_S_frag_AB(i)=10^((log10(M_loss_AB(i))-log10(M_loss(i)))/(Hfr(i)-Hhalf(i))*(Hfr(i))+log10(M_loss(i)));    
        end
        
        
        %%


    % Update Species Reservoirs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Fractions of the atmosphere being eroded
        f_CO2_s = M_CO2_shu / M_at_rem_shu;
        f_N2_s  = M_N2_shu  / M_at_rem_shu;
        f_H2O_s = M_H2O_shu / M_at_rem_shu;

        % Calculate erosion and delivery (with AIRBURSTS)
        M_CO2_shu = M_CO2_shu - (M_loss_S_frag_AB(i) * f_CO2_s) + m_co2_add;
        M_N2_shu  = M_N2_shu  - (M_loss_S_frag_AB(i) * f_N2_s)  + m_n2_add;
        M_H2O_shu = M_H2O_shu - (M_loss_S_frag_AB(i) * f_H2O_s) + m_h2o_add;

        % Update the total mass:
        M_at_rem_shu = M_CO2_shu + M_N2_shu + M_H2O_shu;
        


        % Atmosphere loss fraction [non-dim.]
        loss_frac(i)  = M_loss(i)/M_at;
        
        % Atmosphere loss fraction 2 [non-dim.]
        loss_frac2(i) = M_loss(i)/M_imp(i);
        
        % Read out data for radius of interest
        if(R(i)==r_int)
            M_loss_int    = round(M_loss(i));
            loss_frac_int = loss_frac(i);
        end
        
        
        % B: Melosh & Vickery (1989) + Vickery (1990) model:
        % ==================================================
        
        % Mass of atmosphere above tangent plane [kg] (based on Pham, PhD thesis, eq. (2.4), p. 22, 2012)
        m_tan         = M_at*height_atm_current_val/(2.000*R_p);
        
        % Based on Melosh & Vickery (1989)
        v_escc        = sqrt(8.000/3.000*pival*rho_mean_p*G*R_p^2.000);
        v_min_mv      = sqrt(8.000*((v_escc^2.000)/2.000+H_vap));
        
        % Based on Vickery (1990)  
        m_crit_mv     = 10.000*m_tan;
        r_crit_mv     = (3.000*m_crit_mv/(4.000*pival*rho_imp))^(1.000/3.000);
       
        if(M_imp(i)<m_crit_mv)
            M_loss_mv(i) = 0.000;
            
        elseif(M_imp(i)>=m_crit_mv)
            M_loss_mv(i) = m_tan;
            
        end
        
        % Atmosphere loss fraction [non-dim.]
        loss_frac_mv(i)  = M_loss_mv(i)/M_at;
        
        % Atmosphere loss fraction 2 [non-dim.]
        loss_frac2_mv(i) = M_loss_mv(i)/M_imp(i);
        
         % Data for radius of interest
        if(R(i)==r_int)
            M_loss_mv_int    = round(M_loss_mv(i));
            loss_frac_mv_int = loss_frac_mv(i);      
        end
        

        % D: Kegerreis et al. (2020) model:
        % =================================

        % Impact parameter [non-dim.]
        b            = data.("sin(thet)")(i); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Height of caps [m]
        d(i)         = (R_p+R(i))*(1-b);

        % Calculate volume of impactor cap [m^3]
        V_cap_imp(i) = pival/3.*d(i)^2*(3.*R(i)-d(i));

        % Calculate volume of target cap [m^3]
        V_cap_tar(i) = pival/3.*d(i)^2*(3.*R_p-d(i));

        % Calculate volume of target body [m^3]
        V_p          = 4./3.*pival*R_p^3;

        % Calculate volume of impactor body [m^3]
        V_imp(i)     = 4./3.*pival*R(i)^3;

        % Calculate fractional interacting mass [non-dim.]
        f_M(i)       = (rho_mean_p*V_cap_tar(i)+rho_imp*V_cap_imp(i))/(rho_mean_p*V_p+rho_imp*V_imp(i));

        % Total mass involved in collision [kg]
        M_tot(i)     = M_imp(i)+M_p;

        % Fraction of eroded atmosphere [non-dim.] (eq. 1)
        X(i)         = factorcomp * 0.64*( (v(i)/v_esc(i))^2. * (M_imp(i)/M_tot(i))^(0.5) * (rho_imp/rho_mean_p)^(0.5)*f_M(i))^(0.65);

        % Set upper limit
        if(X(i)>1.000)
             X(i) = 1.000;
        end

        % CUMULATIVE Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Calculate and store the mass lost in each impact step 
        M_loss_ke(i) = X(i) * M_at_rem_ke;
        
        % Update the individual species (Erosion followed by Delivery)
        M_CO2_ke = M_CO2_ke - (X(i) * M_CO2_ke) + m_co2_add;
        M_N2_ke  = M_N2_ke  - (X(i) * M_N2_ke)  + m_n2_add;
        M_H2O_ke = M_H2O_ke - (X(i) * M_H2O_ke) + m_h2o_add;
        
        % Update the master total
        M_at_rem_ke = M_CO2_ke + M_N2_ke + M_H2O_ke;
        



    % %% COPMOSITE

    %compteurs for transition in composite scheme (Shuv->Keg).
       
       iminc = 1.000e5; % Below this radius value, pure Shuvalov model
       imaxc = 10.000e5; % Above this radius value, pure Kegerreis model

    %composite lossscheme
        if (R(i)<iminc)
             M_loss_comp(i)=M_loss_S_frag_AB(i); % use AIRBIRST
         elseif (R(i)>imaxc)
             M_loss_comp(i)=M_loss_ke(i);
         else
              M_loss_comp(i)= 10^((log10(M_loss_ke(i))-log10(M_loss_S_frag_AB(i)))/(log10(imaxc)-log10(iminc))*(log10(R(i))-log10(iminc))+log10(M_loss_S_frag_AB(i)));
              % M_loss_comp(i)= exp((log(M_loss_ke(i))-log(M_loss_S_frag_AB(i)))/(log(imaxc)-log(iminc))*(logR(i)-log(iminc))+log(M_loss_S_frag_AB(i)));
        end  


        % Calculate current fractions for the composite reservoir
        f_CO2_c = M_CO2_comp / M_at_rem_comp;
        f_N2_c  = M_N2_comp  / M_at_rem_comp;
        f_H2O_c = M_H2O_comp / M_at_rem_comp;

        % Update the species (Erosion followed by Delivery)
        M_CO2_comp = M_CO2_comp - (M_loss_comp(i) * f_CO2_c) + m_co2_add;
        M_N2_comp  = M_N2_comp  - (M_loss_comp(i) * f_N2_c)  + m_n2_add;
        M_H2O_comp = M_H2O_comp - (M_loss_comp(i) * f_H2O_c) + m_h2o_add;

        % Update master total
        M_at_rem_comp = M_CO2_comp + M_N2_comp + M_H2O_comp;

        % Update cumulative delivery trackers
        m_h2o_delivered_cum = m_h2o_delivered_cum + m_h2o_add;
        m_co2_delivered_cum = m_co2_delivered_cum + m_co2_add;
        m_n2_delivered_cum  = m_n2_delivered_cum  + m_n2_add;

        impact_idx = impact_idx + 1;
        
    end % End while impacts
    
    % Record History for this timestep
    M_at_history_shu(j)  = M_at_rem_shu;
    M_CO2_history_shu(j) = M_CO2_shu;
    M_N2_history_shu(j)  = M_N2_shu;
    M_H2O_history_shu(j) = M_H2O_shu;
    
    M_at_history_ke(j)   = M_at_rem_ke;
    M_CO2_history_ke(j)  = M_CO2_ke;
    M_N2_history_ke(j)   = M_N2_ke;
    M_H2O_history_ke(j)  = M_H2O_ke;
    
    M_at_history_comp(j)  = M_at_rem_comp;
    M_CO2_history_comp(j) = M_CO2_comp;
    M_N2_history_comp(j)  = M_N2_comp;
    M_H2O_history_comp(j) = M_H2O_comp;
    
    M_H2O_delivered(j) = m_h2o_delivered_cum;
    M_CO2_delivered(j) = m_co2_delivered_cum;
    M_N2_delivered(j)  = m_n2_delivered_cum;
    
    % Track atmospheric properties
    M_at_curr_step = M_CO2_shu + M_N2_shu + M_H2O_shu;
    moles_total_step = (M_H2O_shu/18) + (M_N2_shu/28) + (M_CO2_shu/44);
    mu_history(j) = M_at_curr_step / moles_total_step;
    height_atm_history(j) = height_atm * (mu / mu_history(j));

end % End timestep loop

% Ensure we finish the bar correctly
if next_prog_update < prog_bar_width
    fprintf(repmat('\b', 1, prog_bar_width + 1));
    fprintf('%s]', repmat('#', 1, prog_bar_width));
end
fprintf(' Done!\n');


%Second run for composite

%radius loop
%for i=1:1:i_max
%composite
% if (R(i)<1.500e5);
%     M_loss_comp(i)=M_loss(i);
% elseif (R(i)>8.000e5);
%     M_loss_comp(i)=M_loss_ke(i);
% else
%     M_loss_comp(i)= 10^((log10(M_loss_ke(imaxc))-log10(M_loss(iminc)))/(imaxc-iminc)*(i-iminc)+log10(M_loss(iminc)));
    
% end  
%end


% % Include warning messages in case impactor radius and/or velocity is insufficient
% if(r_int<r_crit_mv && v >= v_min_mv)
%     disp(' ')
%     disp('Impactor radius insufficient to blow off atmosphere above tangent plane! Velocity is ok.')
% elseif(r_int>r_crit_mv && v < v_min_mv)
%     disp(' ')
%     disp('Impactor velocity insufficient to blow off atmosphere above tangent plane! Radius is ok.')    
% elseif(r_int<r_crit_mv && v < v_min_mv)
%     disp(' ')
%     disp('Impactor radius and velocity insufficient to blow off atmosphere above tangent plane!')
% end    


time_Myr_impacts = data.time / 1e6;
time_Myr_steps   = t_vec / 1e6;


% Figure: Impact Distribution
figure("Name","Impact Distribution");
stairs(time_Myr_impacts, cumsum(M_imp_all), 'k', 'LineWidth', 6);
xlabel('Time (Myr after CAI)', 'FontWeight', 'bold');
ylabel('Cumulative Impact Mass (kg)', 'FontWeight', 'bold');
legend('Scenario: 46 Impacts', 'Location', 'southeast', 'FontSize', legend_fs, 'FontWeight', 'bold', 'Color', 'none');
grid on; title('Total Impact Mass Delivered', 'FontWeight', 'bold');
xlim(plot_xlim);

% Enhanced Axis Style
set(gca, 'FontSize', axes_fs, 'FontWeight', 'bold', 'LineWidth', border_lw, 'Box', 'on');

% Transparency and Sizing
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, fig_width, fig_height]);
set(gcf, 'PaperSize', [fig_width, fig_height]);
set(gcf, 'Color', bg_color); set(gca, 'Color', bg_color);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters');

% Export Configured Format
print(gcf, ['Impact Distribution', export_ext], export_driver, export_res);


% Figure: Shuvalov Evolution
figure("Name","Shuvalov Evolution");
semilogy(time_Myr_steps, M_CO2_history_shu, 'r', 'LineWidth', plot_lw); hold on;
semilogy(time_Myr_steps, M_N2_history_shu,  'b', 'LineWidth', plot_lw);
semilogy(time_Myr_steps, M_H2O_history_shu, 'k', 'LineWidth', plot_lw);
semilogy(time_Myr_steps, M_at_history_shu,  '--', 'Color', [1, 0.5, 0], 'LineWidth', plot_lw + 1); 
xlabel('Time (Myr after CAI)', 'FontWeight', 'bold');
ylabel('Species Mass (kg)', 'FontWeight', 'bold');
title(['Species Evolution: Shuvalov Model (', comp_label, ')'], 'FontWeight', 'bold');
legend('CO_2', 'N_2', 'H_2O','Total Atmosphere', 'Location', 'southeast', 'FontSize', legend_fs, 'FontWeight', 'bold', 'Color', 'none');
grid on;
xlim(plot_xlim);

% Fix artifacts by disabling minor grid for log plots
set(gca, 'FontSize', axes_fs, 'FontWeight', 'bold', 'LineWidth', border_lw, 'Box', 'on', 'MinorGridLineStyle', 'none');
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, fig_width, fig_height]);
set(gcf, 'PaperSize', [fig_width, fig_height]);
set(gcf, 'Color', bg_color); set(gca, 'Color', bg_color);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters');

print(gcf, ['Shuvalov', export_ext], export_driver, export_res);

% Figure: Kegerreis Evolution
figure("Name","Kegerreis Evolution");
semilogy(time_Myr_steps, M_CO2_history_ke, 'r', 'LineWidth', plot_lw); hold on;
semilogy(time_Myr_steps, M_N2_history_ke,  'b', 'LineWidth', plot_lw);
semilogy(time_Myr_steps, M_H2O_history_ke, 'k', 'LineWidth', plot_lw);
semilogy(time_Myr_steps, M_at_history_ke,  '--', 'Color', [1, 0.5, 0], 'LineWidth', plot_lw + 1); 
xlabel('Time (Myr after CAI)', 'FontWeight', 'bold');
ylabel('Species Mass (kg)', 'FontWeight', 'bold');
title(['Species Evolution: Kegerreis Model (', comp_label, ')'], 'FontWeight', 'bold');
legend('CO_2', 'N_2', 'H_2O','Total Atmosphere','Location', 'southeast', 'FontSize', legend_fs, 'FontWeight', 'bold', 'Color', 'none');
grid on;
xlim(plot_xlim);

set(gca, 'FontSize', axes_fs, 'FontWeight', 'bold', 'LineWidth', border_lw, 'Box', 'on', 'MinorGridLineStyle', 'none');
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, fig_width, fig_height]);
set(gcf, 'PaperSize', [fig_width, fig_height]);
set(gcf, 'Color', bg_color); set(gca, 'Color', bg_color);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters');

print(gcf, ['Kegerreis', export_ext], export_driver, export_res);


% Figure: Composite Evolution
figure("Name","Composite Evolution");
semilogy(time_Myr_steps, M_CO2_history_comp, 'r', 'LineWidth', plot_lw); hold on;
semilogy(time_Myr_steps, M_N2_history_comp,  'b', 'LineWidth', plot_lw);
semilogy(time_Myr_steps, M_H2O_history_comp, 'k', 'LineWidth', plot_lw);
semilogy(time_Myr_steps, M_at_history_comp,  '--', 'Color', [1, 0.5, 0], 'LineWidth', plot_lw + 1); 
xlabel('Time (Myr after CAI)', 'FontWeight', 'bold');
ylabel('Species Mass (kg)', 'FontWeight', 'bold');
title(['Species Evolution: Composite Model (', comp_label, ')'], 'FontWeight', 'bold');
legend('CO_2', 'N_2', 'H_2O','Total Atmosphere', 'Location', 'southeast', 'FontSize', legend_fs, 'FontWeight', 'bold', 'Color', 'none');
grid on;
xlim(plot_xlim);

set(gca, 'FontSize', axes_fs, 'FontWeight', 'bold', 'LineWidth', border_lw, 'Box', 'on', 'MinorGridLineStyle', 'none');
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, fig_width, fig_height]);
set(gcf, 'PaperSize', [fig_width, fig_height]);
set(gcf, 'Color', bg_color); set(gca, 'Color', bg_color);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters');

print(gcf, ['Composite', export_ext], export_driver, export_res);


% Figure: Model Comparison (All)
figure("Name","Compare All Models");
plot(time_Myr_steps, M_at_history_shu, 'r', 'LineWidth', plot_lw); hold on;
plot(time_Myr_steps, M_at_history_ke, 'b', 'LineWidth', plot_lw);
plot(time_Myr_steps, M_at_history_comp, 'Color', [1, 0.5, 0], 'LineWidth', plot_lw, 'LineStyle', '-'); 
yline(4.8e20, 'k--', 'LineWidth', plot_lw - 1);
xlabel('Time (Myr after CAI)', 'FontWeight', 'bold');
ylabel('Total Atmosphere Mass (kg)', 'FontWeight', 'bold');
legend('Shuvalov Model', 'Kegerreis Model', 'Composite Model', 'Current Venus Atmosphere', 'Location', 'northwest', 'FontSize', legend_fs, 'FontWeight', 'bold', 'Color', 'none');
title(['Atmosphere Evolution: Comparison (', comp_label, ')'], 'FontWeight', 'bold');
grid on;
xlim(plot_xlim);

set(gca, 'FontSize', axes_fs, 'FontWeight', 'bold', 'LineWidth', border_lw, 'Box', 'on');
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, fig_width, fig_height]);
set(gcf, 'PaperSize', [fig_width, fig_height]);
set(gcf, 'Color', bg_color); set(gca, 'Color', bg_color);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters');

print(gcf, ['All_Models_Comparison', export_ext], export_driver, export_res);