function h = h_He_NB(current_pressure, T_GGG_surface, properties)
    % Helium heat transfer coefficient calculation
    % Includes both nucleate boiling (Z.Q. Long and P. Zhang correlation)
    % and conduction through liquid layer
    %
    % Author: yubin-kimm
    % Last modified: 2025-02-08 12:08:18 UTC
    %
    % Inputs:
    %   current_pressure - System pressure [Pa]
    %   T_GGG_surface   - GGG surface temperature (first node) [K]
    %   properties      - Table containing helium properties
    %
    % Output:
    %   h              - Heat transfer coefficient [W/(m²·K)]
    
    global L_gap  % Gap between thermosyphon and GGG [m]
    
    % 데이터 추출
    T_data = properties.Temperature;
    P_data = properties.Pressure;
    rho_l_data = properties{:, 3};   % Liquid Density
    rho_v_data = properties{:, 4};   % Vapor Density
    h_l_data = properties{:, 5};     % Liquid Enthalpy
    h_v_data = properties{:, 6};     % Vapor Enthalpy
    cp_l_data = properties{:, 9};    % Liquid Cp
    k_l_data = properties{:, 11};    % Liquid Thermal Conductivity
    mu_l_data = properties{:, 13};   % Liquid Viscosity
    sigma_data = properties{:, 15};  % Surface Tension
    
    % Pa를 kPa로 변환 (상관식이 kPa 기준으로 만들어짐)
    P_kPa = current_pressure / 1000;
    
    % 보간 함수를 통한 물성치 계산
    T_sat = interp1(P_data/1000, T_data, P_kPa, 'linear', 'extrap');  % P_data도 kPa로 변환
    k_l = interp1(T_data, k_l_data, T_sat, 'linear', 'extrap');
    
    % 온도차 계산
    del_T = T_GGG_surface - T_sat;
    
    if del_T <= 0
        % When wall temperature is lower than saturation temperature
        % Use conduction through the gap between thermosyphon and GGG
        h = k_l / L_gap;  % [W/(m²·K)]
    else
        % Nucleate boiling regime - Using Z.Q. Long and P. Zhang correlation
        rho_l = interp1(T_data, rho_l_data, T_sat, 'linear', 'extrap');
        rho_v = interp1(T_data, rho_v_data, T_sat, 'linear', 'extrap');
        cp_l = interp1(T_data, cp_l_data, T_sat, 'linear', 'extrap');
        mu_l = interp1(T_data, mu_l_data, T_sat, 'linear', 'extrap');
        sigma = interp1(T_data, sigma_data, T_sat, 'linear', 'extrap');
        h_fg = interp1(T_data, h_v_data - h_l_data, T_sat, 'linear', 'extrap');
        
        % 물리 상수
        g = 9.81;    % 중력가속도 [m/s²]
        
        % Z.Q. Long과 P. Zhang의 실험식을 통한 열전달 계수 계산
        h = ((3.25e-4 * ...
              ((del_T * cp_l * rho_l)/(h_fg * rho_v * k_l) * (sigma/(g * rho_l))^0.5)^0.6 * ...
              (g * (rho_l/mu_l)^2 * (sigma/(g * rho_l))^1.5)^0.125 * ...
              (1e3 * P_kPa/(sigma * g * rho_l)^0.5)^0.7) * ...
              k_l * (g * rho_l/sigma)^0.5)^2.5;
        
        % 임계열유속(CHF) 체크
        q = h * del_T;
        q_CHF = 0.16 * h_fg * rho_v * (sigma * g * (rho_l - rho_v))^0.25;
        
        if q > q_CHF
            warning('Heat flux (%.2f W/m²) exceeds Critical Heat Flux (%.2f W/m²)!', q, q_CHF);
        end
    end
end