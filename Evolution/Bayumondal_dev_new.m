elseif (R(i) > imaxc)
            % 3b. Pure Kegerreis evaluated on COMPOSITE atmosphere
            M_loss_comp(i) = M_loss_ke_comp;
        else
            % 3c. Sweep radii from iminc to R(i) using COMPOSITE atmospheric properties
            vrel_vesc_i = data.("vrel/vesc")(i);
            n_sweep = 30;
            R_sweep = logspace(log10(iminc), log10(R(i)), n_sweep)';
            M_shu_sweep = zeros(n_sweep, 1);

            for ks = 1:n_sweep
                Rk = R_sweep(ks);
                Mk = 4/3 * pival * Rk^3 * rho_imp;
                v_esc_k = sqrt(2*G*(M_p + Mk)/(R_p + Rk));
                vk = vrel_vesc_i * v_esc_k;

                % Shuvalov ground-impact (using height_atm_comp and rho_atm_comp)
                xi_k = (2*Rk)^3*rho_imp*rho_p*(vk^2-v_esc_k^2) / ...
                       (height_atm_comp^3*rho_atm_comp*v_esc_k^2*(rho_p+rho_imp));
                xi_k = max(xi_k, 1000);
                chi_k = 10^(-6.375+5.239*log10(xi_k)-2.121*(log10(xi_k)^2) ...
                        +0.397*(log10(xi_k)^3)-0.037*(log10(xi_k)^4)+0.0013*(log10(xi_k)^5));
                M_loss_k = chi_k * Mk * ((vk^2-v_esc_k^2)/v_esc_k^2);

                % Airburst [Shuvalov (2014)] (using height_atm_comp and rho_atm_comp)
                Hhalf_k = -1.4*height_atm_comp*log(0.5*(1/height_atm_comp)*5.5*(2*Rk)*(rho_imp/rho_atm_comp)^0.5);
                Hfr_k   = -1.4*height_atm_comp*log(0.5*(1/height_atm_comp)*3.5*(2*Rk)*(rho_imp/rho_atm_comp)^0.5);
                Vlim_k  = max(15+0.2*height_atm_comp/1000, 0);
                M_loss_AB_k = Mk*5.5*0.001*(2*Rk/height_atm_comp)*((rho_imp/1000)^0.5)*((vk/1000-Vlim_k)^2);
                M_loss_fr_k = 10^((log10(M_loss_AB_k)-log10(M_loss_k))/(Hfr_k-Hhalf_k)*(Hfr_k)+log10(M_loss_k));

                % Select regime
                if (Hhalf_k > 0)
                    val_k = M_loss_AB_k;
                elseif (Hhalf_k < 0) && (Hfr_k > 0)
                    val_k = M_loss_fr_k;
                else
                    val_k = M_loss_k;
                end

                % Tangent plane limit (using m_tan_comp)
                val_k = min(val_k, m_tan_comp);

                % Monotonic clamp
                if ks > 1 && val_k < M_shu_sweep(ks-1)
                    val_k = M_shu_sweep(ks-1);
                end
                M_shu_sweep(ks) = val_k;
            end

            % Extrapolated Shuvalov value at R(i)
            M_loss_shu_extrap = M_shu_sweep(end);

            % Composite interpolation (using M_loss_ke_comp and M_loss_shu_extrap)
            M_loss_comp(i) = 10^((log10(M_loss_ke_comp)-log10(M_loss_shu_extrap))/(log10(imaxc)-log10(iminc))*(log10(R(i))-log10(iminc))+log10(M_loss_shu_extrap));
        end
