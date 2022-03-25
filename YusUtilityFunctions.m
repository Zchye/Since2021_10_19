classdef YusUtilityFunctions < handle
    %
    methods
        function obj = YusUtilityFunctions( )
            % Constructor
        end
        
        % Methods for calculating pathloss
        function [Pathloss, sigma_SF]=calculatePathloss(obj, param, gNBPos, UEStat, LinkDir)
            %Calculate pathloss with formulas in Table 7.4.1-1 for each BS-UT link to be modelled
            %param - a structure containing the field Scenario
            %gNBPos - the Cartesian coordinates of a gNB, must be a 3-D column
            %vector
            %UEStat - a structure containing states of a UE
            %get the parameters in 3GPP TR 38.901 Figure7.4.1-1

            % Parse UEStat
            UEPos = UEStat.UEPosition;
            LOS = UEStat.LOS;
            
            %get distance and heights
            d_2D=norm(gNBPos(1:2)-UEPos(1:2));
            d_3D=norm(gNBPos(1:3)-UEPos(1:3));
            h_BS = gNBPos(3);
            h_UT=UEPos(3);
            

            %constants
            if ~LinkDir %For DL
                f_c = param.DLCarrierFreq/1e9; %f_c is carrier frequency in GHz
            else %For UL
                f_c = param.ULCarrierFreq/1e9; %f_c is carrier frequency in GHz
            end


            %select the function to calculate LOS probability
            switch param.Scenario
                case 'UMi'
                    h=@obj.UMiPathloss;
                case 'UMa'
                    h=@obj.UMaPathloss;
                case 'RMa'
                    h=@obj.RMaPathloss;
                otherwise
                    error('Scenarios other than RMa, UMi, UMa are yet to be implemented.');
            end

            % Calculate Pathloss
            % PL_b is the basic outdoor pathloss
            [PL_b, sigma_SF]=h(d_2D,d_3D,h_BS,h_UT,f_c,LOS);
            mu = 0; %band-aid
            if ~UEStat.Indoor % UE is outdoor
                PL_tw = 0; % No penetration loss for outdoor UEs
                PL_in = 0; % No inside loss for outdoor UEs
                if UEStat.InCar % UE is in a car
                    % See TR 39.901 section 7.4.3.2 for μ and σ_P. 
                    % Optionally, for metallized car windows, μ = 20 can be used.
                    mu = 9;
                    sigma_P = 5;
                end
            else % UE is indoor
                % Material penetration losses, see TR 38.901 Table 7.4.3-1
                L_glass = 2+0.2*f_c;
                L_concrete = 5+4*f_c;
                % O2I building penetration loss, see TR 38.901 Table
                % 7.4.3-2
                if UEStat.HighLoss % High penetration loss
                    PL_tw = 5 - 10*log10(0.7*10^(-L_glass/10) + 0.3*10^(-L_concrete/10));
                    sigma_P = 6.5;
                else % Low penetration loss
                    PL_tw = 5 - 10*log10(0.3*10^(-L_glass/10) + 0.7*10^(-L_concrete/10));
                    sigma_P = 4.4;
                end
                PL_in = 0.5*UEStat.d2Din; % Inside loss
            end
            
            Pathloss = PL_b + PL_tw + PL_in + UEStat.Indoor*normrnd(0,sigma_P) +...
                (~UEStat.Indoor)*UEStat.InCar*normrnd(mu, sigma_P);



        end

        function [Pathloss, sigma_SF]=UMiPathloss(~, d_2D,d_3D,h_BS,h_UT,f_c,LOS)
            %see 3GPP TR 38.901 Table 7.4.1-1

            %generate coefficients
            c = 300000000;
            h_e = 1;
            h_bs_prime = h_BS - h_e;
            h_ut_prime = h_UT - h_e;
            d_bp_prime = (4 * h_bs_prime * h_ut_prime * f_c * 1000000000) / c;

            %calculate LOS pathloss
            if (d_2D >= 10) && (d_2D <= d_bp_prime)
                PL_umi_los = 32.4 + 21*log10(d_3D) + 20*log10(f_c);
            elseif (d_2D >= d_bp_prime) && (d_2D <=5000)
                PL_umi_los = 32.4 + 40*log10(d_3D) + 20*log10(f_c) - 9.5*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
            else
        %         error('d_2D out of bound.');
                PL_umi_los = 32.4 + 40*log10(d_3D) + 20*log10(f_c) - 9.5*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
            end 

            if LOS
                Pathloss = PL_umi_los;
                sigma_SF = 4;
            else
                %calculate NLOS pathloss
                PL_umi_nlos_prime = 35.3*log10(d_3D)+22.4+21.3*log10(f_c)-0.3*(h_UT-1.5);
                Pathloss = max(PL_umi_los, PL_umi_nlos_prime);
                sigma_SF = 7.82;
            end
        end

        function [Pathloss, sigma_SF]=UMaPathloss(~, d_2D,d_3D,h_BS,h_UT,f_c,LOS)
            %see 3GPP TR 38.901 Table 7.4.1-1

            %generate coefficients
            c = 300000000;
            h_e = 1;
            h_bs_prime = h_BS - h_e;
            h_ut_prime = h_UT - h_e;
            d_bp_prime = (4 * h_bs_prime * h_ut_prime * f_c * 1000000000) / c;

            %calculate LOS pathloss
            if (d_2D >= 10) && (d_2D <= d_bp_prime)
                PL_uma_los = 28.0 + 22*log10(d_3D) + 20*log10(f_c);
            elseif (d_2D >= d_bp_prime) && (d_2D <=5000)
                PL_uma_los = 28.0 + 40*log10(d_3D) + 20*log10(f_c) - 9*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
            else
        %         error('d_2D out of bound.');
                PL_uma_los = 28.0 + 40*log10(d_3D) + 20*log10(f_c) - 9*log10((d_bp_prime)^2+(h_BS-h_UT)^2);
            end 

            if LOS
                Pathloss = PL_uma_los;
                sigma_SF = 4;
            else
                %calculate NLOS pathloss
                PL_uma_nlos_prime = 13.54 + 39.08*log10(d_3D)+20*log10(f_c)-0.6*(h_UT-1.5);
                Pathloss = max(PL_uma_los, PL_uma_nlos_prime);
                sigma_SF = 6;
            end
        end

        function [Pathloss, sigma_SF]=RMaPathloss(~, d_2D,d_3D,h_BS,h_UT,f_c,LOS)
            %see 3GPP TR 38.901 Table 7.4.1-1

            %generate coefficients
            c = 300000000;
            h = 5;
            w = 20;
            d_bp = (2 * pi * h_BS * h_UT * f_c * 1000000000) / c;

            %calculate LOS pathloss
            if (d_2D >= 10) && (d_2D <= d_bp)
                PL_rma_los = 20*log10(40*pi*d_3D*f_c/3)+min(0.03*(h^1.72),10)*log10(d_3D)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_3D;
                sigma = 4;
            elseif (d_2D >= d_bp) && (d_2D <=10000)
                PL_rma_los = 20*log10(40*pi*d_bp*f_c/3)+min(0.03*(h^1.72),10)*log10(d_bp)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_bp+40*log10(d_3D/d_bp);
                sigma = 6;
            else
        %         error('d_2D out of bound.');
                PL_rma_los = 20*log10(40*pi*d_bp*f_c/3)+min(0.03*(h^1.72),10)*log10(d_bp)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d_bp+40*log10(d_3D/d_bp);
                sigma = 6;
            end 

            if LOS
                Pathloss = PL_rma_los;
                sigma_SF = sigma;
            else
                %calculate NLOS pathloss
                PL_rma_nlos_prime = 161.04-7.1*log10(w)+7.5*log10(h)-(24.34-((h/h_BS)^2))*log10(h_BS)+(43.42-3.1*log10(h_BS))*(log10(d_3D)-3)+20*log10(f_c)-(3.2*((log10(11.75*h_UT))^2)-4.97);
                Pathloss = max(PL_rma_los, PL_rma_nlos_prime);
                sigma_SF = 8;
            end
        end
    end
end