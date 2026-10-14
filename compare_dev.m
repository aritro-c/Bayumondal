
clear all; clc;


% Choose planet: 0 - Earth, 1 - Venus, 2 - Mars

planet         = 1;           

% Impactor typical parameters
rho_imp        = 2700.000;    % Impactor density [kg/m^3]
vrel_vesc      = 2.0;         % Typical impact velocity ratio v/v_esc
b              = sin(deg2rad(45));  % Impact Angle is 45 degree

pival          = pi; % Pi [non-dim.]
G              = 6.67384e-11; % Gravitational constant [m^/(kg s^2)]


% AIRBURST for Shuvalov
use_AB         = true;     % Set to true to use Airburst for Shuvalov


% COMPOSITE Scheme Limits
iminc = 50e3; % [meters] Below this radius value, pure Shuvalov model 
imaxc = 2443e3; % [meters] Above this radius value, pure Kegerreis model 




% Physical parameters of Earth as target body
if planet==0
    M_p        = 5.9736e24;   % Earth mass [kg]
    R_p        = 6.3781e6;    % Earth radius [m]
    rho_mean_p = 5515.000;    % Mean density of planet [kg/m^3]
    rho_p      = 2900.000;    % Density of crust [kg/m^3] (assumes that basaltic crust is target material)
    M_at       = 5.100e18;    % Earth atmosphere mass [kg]
    rho_atm    = 1.217;       % Surface density of Earth's atmosphere [kg/m^3]
    height_atm = 8500.0;      % Scale height of Earth atmosphere [m]
    mu         = 29;          % molecular mass of (N2 + O2) atmosphere [non-dim.] 

    M_N2  = M_at * 0.7808; 
    M_H2O  = M_at * 0.01; 
    M_CO2 = M_at * 0.0004; 


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
    mu         = 44;          % molecular mass of pure CO2 atmosphere [non-dim.] 

    M_N2  = M_at * 0.0000; 
    M_H2O  = M_at * 0.0000; 
    M_CO2 = M_at * 1.0000; 


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
    mu         = 44;          % molecular mass of pure CO2 atmosphere [non-dim.] 
    
    M_N2  = M_at * 0.027; 
    M_H2O  = M_at * 0.000; 
    M_CO2 = M_at * 0.9500; 


end



height_atm_current_val = height_atm;
rho_atm_current = rho_atm;

%factor to vary composition for Kegerrreis law.
factorcomp=   1.00; %0.3 %1.00




% ================= Radius Ragne ================= 


% Choose the minimum and maximum impactor radius [meters]
radius_min = 10;   
radius_max = 10e6;  

% Range of Impactor Radii [m]
R = logspace(log10(radius_min), log10(radius_max), 1000)'; 
num_pts = length(R);





% ================== Plot Config ====================

    % Figure size (inches)
    fig_width  = 18; % Default is 10
    fig_height = 12;  % Default is 14
    
    % Plot Style
    plot_lw    = 5;      % LineWidth for plots
    border_lw  = 1.5;    % LineWidth for axes border
    legend_fs  = 18;     % FontSize for legends
    axes_fs    = 28;     % FontSize for axes labels
    
    % Export Settings
    export_format = 'png'; % 'svg' or 'png'
    export_dpi    = 300;   % DPI for PNG export


% =======================================================





M_loss_comp = zeros(num_pts,1);
M_loss_ke = zeros(num_pts,1);
M_loss_S_frag_AB = zeros(num_pts,1);
M_loss_S_frag_AB_extrap = zeros(num_pts,1);
M_imp = zeros(num_pts,1);


    
% Mass of atmosphere above tangent plane [kg] 
% the maximum atmospheric mass that they can eject in a single impact is given by all the mass above the tangent plane, 
% which is H/2R of the total atmosphere, where H is the scale height, and R is the radius of the planet (Schlichting+ 2015)                                                                                                    
m_tan         = M_at*height_atm_current_val/(2.000*R_p);





for i = 1:num_pts
    % 1. Impactor Mass
    M_imp(i) = 4/3 * pival * R(i)^3 * rho_imp;

    % Mutual escape velocity [m/s]
    v_esc = sqrt(2*G*(M_p+M_imp(i))/(R_p+R(i)));

    % Impact velocity [m/s]
    v = vrel_vesc * v_esc;




    % =======================================================
    % A: Shuvalov (2009, 2010) + Airburst [Shuvalov (2014)]
    % =======================================================
    


    % Impact efficiency [non-dim.] (Shuvalov's eq. (2))
    xi = (2.000*R(i))^3*rho_imp*rho_p*(v^2-v_esc^2)/(height_atm_current_val^3*rho_atm_current*v_esc^2*(rho_p+rho_imp));
    xi = max(xi, 1000); % Without this threshold the plot looks weird
    
    % Polynomial fit [non-dim.]
    chi = 10^(-6.375+5.239*log10(xi)-2.121*(log10(xi)^2)+0.397*(log10(xi)^3)-0.037*(log10(xi)^4)+0.0013*(log10(xi)^5));  
    
    % Atmospheric mass loss [kg]
    M_loss = chi * M_imp(i) * ((v^2-v_esc^2)/v_esc^2);
    
    % Airbursts
    Hhalf = -1.4*height_atm_current_val* log(0.5*(1/height_atm_current_val)*5.5*(2.000*R(i))*(rho_imp/rho_atm_current)^0.5); % [Eq (7) in Shuvalov+ 2014]
    Hfr   = -1.4*height_atm_current_val* log(0.5*(1/height_atm_current_val)*3.5*(2.000*R(i))*(rho_imp/rho_atm_current)^0.5); % [Eq (8) in Shuvalov+ 2014]
       
    Vlim  = (15+0.2*height_atm_current_val/1000); % [km] [Eq (10) in Shuvalov+ 2014]
    Vlim  = max(Vlim, 0); % Vlim shouldn't be negative

    
    M_loss_AB = M_imp(i)*5.5*0.001*(2.000*R(i)/height_atm_current_val)*((rho_imp/1000)^0.5)*(((v/1000)-Vlim)^2.); % [From eq (9) in Shuvalov+ 2014]
 
    M_loss_fr = 10^((log10(M_loss_AB)-log10(M_loss))/(Hfr-Hhalf)*(Hfr)+log10(M_loss)); % [ From eq (11) in Shuvalov+ 2014]
    
    if use_AB 

        %  If Hhalf > 0 the impact is regarded as aerial burst,
        %  and Eqs. (9) and (10) are used to work out how much atmosphere is eroded (Shuvalov+, 2014)
        if(Hhalf > 0) 
            M_loss_S_frag_AB(i) = M_loss_AB; 

        % Estimates for impacts of fragmented projectiles (Hhalf < 0; Hfr > 0) are obtained from (11) (Shuvalov+, 2014) 
        elseif(Hhalf < 0) && (Hfr > 0)

            M_loss_S_frag_AB(i) = M_loss_fr; 
        else
             
            M_loss_S_frag_AB(i) = M_loss;
        end
    
    else
        M_loss_S_frag_AB(i) = M_loss;

    end

    % Tangent Plane limit for Shuvalov
    if M_loss_S_frag_AB(i) > m_tan
        M_loss_S_frag_AB(i) = m_tan;
    end

    % EXTRAPOLATION
    M_loss_S_frag_AB_extrap(i) = M_loss_S_frag_AB(i);
    
    if i > 1 && M_loss_S_frag_AB_extrap(i) < M_loss_S_frag_AB_extrap(i-1)
        M_loss_S_frag_AB_extrap(i) = M_loss_S_frag_AB_extrap(i-1);
    end
    



    % ==========================================================
    % Kegerreis (2020)
    % ==========================================================
    
    d_cap = (R_p+R(i))*(1-b); % Height of caps [m]
    V_cap_imp = pival/3.*d_cap^2*3.*R(i)-d_cap;  % Calculate volume of impactor cap [m^3]
    V_cap_tar = pival/3.*d_cap^2*(3.*R_p-d_cap);  % Calculate volume of target cap [m^3]
    V_p       = 4./3.*pival*R_p^3; % Calculate volume of target body [m^3]
    V_imp     = 4./3.*pival*R(i)^3; % Calculate volume of impactor body [m^3]

    V_cap_imp(d_cap >= 2.*R(i)) = V_imp (d_cap >= 2.*R(i)); % When d_cap >= Impactor_diameter, Impactor_volume_cap = Impactor_volume ;
    V_cap_tar(d_cap >= R_p) = V_p; % When d_cap >= Planet_diameter, Planet_volume_cap = Planet_volume (This limit is probably unnecessary for our use case);

    f_M       = (rho_mean_p*V_cap_tar+rho_imp*V_cap_imp)/(rho_mean_p*V_p+rho_imp*V_imp);  % Calculate fractional interacting mass [non-dim.]
    M_tot     = M_imp(i)+M_p;  % Total mass involved in collision [kg]
    
    X_frac = factorcomp * 0.64*( (v/v_esc)^2. * (M_imp(i)/M_tot)^(0.5) * (rho_imp/rho_mean_p)^(0.5)*f_M)^(0.65); % Fraction of eroded atmosphere [non-dim.] (eq. 1)
    
    % Set upper limit
    if(X_frac > 1.000)
        X_frac = 1.000;
    end
    
    M_loss_ke(i) = X_frac * M_at;
    
    % ==========================================================
    % Composite Scheme
    % ==========================================================
    
    if (R(i) < iminc)
        M_loss_comp(i) = M_loss_S_frag_AB_extrap(i);
    elseif (R(i) > imaxc)
        M_loss_comp(i) = M_loss_ke(i);
    else
        M_loss_comp(i) = 10^((log10(M_loss_ke(i))-log10(M_loss_S_frag_AB_extrap(i)))/(log10(imaxc)-log10(iminc))*(log10(R(i))-log10(iminc))+log10(M_loss_S_frag_AB_extrap(i)));
    end
end



    


if use_AB 

     shuv_label = 'Shuvalov + Airburst';

     file_name = 'Impactor_Radius_vs_Erosion_AB';

else
     shuv_label = 'Shuvalov';
     file_name = 'Impactor_Radius_vs_Erosion';
end



% ==========================================================
% Plotting Export Logic
% ==========================================================
if strcmp(export_format, 'svg')
    export_ext = '.svg';
    export_driver = '-dsvg';
    export_res = '-painters'; 
    bg_color = 'none';
else
    export_ext = '.png';
    export_driver = '-dpng';
    export_res = ['-r', num2str(export_dpi)];
    bg_color = 'w'; 
end




% ==========================================================
% Plot
% ==========================================================
figure('Name', 'Impactor Radius vs Erosion');

% Convert Radius to km for X-axis
R_km = R / 1000; 

% Plot individual and composite models
loglog(R_km, M_loss_comp, '-','Color', [1 0.5 0], 'LineWidth', plot_lw);
hold on;
loglog(R_km, M_loss_S_frag_AB_extrap, 'k--', 'LineWidth', plot_lw-2);
loglog(R_km, M_loss_S_frag_AB, 'r--', 'LineWidth', plot_lw-2);
loglog(R_km, M_loss_ke, 'b--', 'LineWidth', plot_lw-2);

% Transition boundaries
xline(iminc/1000, 'k:', 'LineWidth', plot_lw-2, 'LabelVerticalAlignment', 'bottom','Label', iminc/1000 + " km",'FontSize', 30);
xline(imaxc/1000, 'k:', 'LineWidth', plot_lw-2, 'LabelVerticalAlignment', 'bottom', 'Label', imaxc/1000 + " km", 'FontSize', 30);

xlabel('Impactor Radius (km)', 'FontWeight', 'bold');
ylabel('Atmospheric Mass Loss (kg)', 'FontWeight', 'bold');
title('Impactor Radius vs Atmospheric Erosion', 'FontWeight', 'bold');
legend('Composite', 'Linear-Extrapolation for Shuvalov', shuv_label, 'Kegerreis et al.', 'Location', 'northwest', 'FontSize', legend_fs, 'FontWeight', 'bold', 'Color', 'none');
grid on;

% Axis Style
set(gca, 'FontSize', axes_fs, 'FontWeight', 'bold', 'LineWidth', border_lw, 'Box', 'on');

% Transparency and Sizing
set(gcf, 'Units', 'inches', 'Position', [1, 1, fig_width, fig_height]);
set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, fig_width, fig_height]);
set(gcf, 'PaperSize', [fig_width, fig_height]);
set(gcf, 'Color', bg_color); set(gca, 'Color', bg_color);
set(gcf, 'InvertHardcopy', 'off');
set(gcf, 'Renderer', 'painters');

% Export Configured Format
print(gcf, [file_name, export_ext], export_driver, export_res);
