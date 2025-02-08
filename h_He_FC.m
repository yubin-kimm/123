function h = h_He_FC(P, T_surface, Lcond)
    % Helium film condensation heat transfer coefficient calculation
    % Based on Nusselt theory for vertical plate condensation
    %
    % Inputs:
    %   P          - System pressure [kPa]
    %   T_surface  - Surface temperature [K]
    %   Lcond      - Characteristic length for condensation [m]
    
    % 엑셀 데이터 읽기
    opts = detectImportOptions('properties.xlsx');
    opts.VariableNamingRule = 'preserve';
    properties = readtable('properties.xlsx', opts);
    
    % 데이터 추출
    T_data = properties.Temperature;
    P_data = properties.Pressure;
    rho_l_data = properties{:, 3};  % Liquid Density
    rho_v_data = properties{:, 4};  % Vapor Density
    h_l_data = properties{:, 5};    % Liquid Enthalpy
    h_v_data = properties{:, 6};    % Vapor Enthalpy
    k_l_data = properties{:, 11};   % Liquid Thermal Conductivity
    mu_l_data = properties{:, 13};  % Liquid Viscosity
    sigma_data = properties{:, 15}; % Surface Tension
    
    % 보간 함수 생성
    T_sat = interp1(P_data, T_data, P, 'linear', 'extrap');
    rho_l = interp1(T_data, rho_l_data, T_sat, 'linear', 'extrap');
    rho_v = interp1(T_data, rho_v_data, T_sat, 'linear', 'extrap');
    k_l = interp1(T_data, k_l_data, T_sat, 'linear', 'extrap');
    mu_l = interp1(T_data, mu_l_data, T_sat, 'linear', 'extrap');
    sigma = interp1(T_data, sigma_data, T_sat, 'linear', 'extrap');
    h_fg = interp1(T_data, h_v_data - h_l_data, T_sat, 'linear', 'extrap');
    
    % 온도차 계산
    del_T = T_sat - T_surface;
    
    % 안전성 검사
    if del_T <= 0
        warning('Surface temperature is higher than or equal to saturation temperature');
        h = 0;
        return;
    end
    
    % 물리 상수
    g = 9.81;    % 중력가속도 [m/s²]
    
    % Nusselt 응축 이론을 통한 열전달 계수 계산
    h = 0.789 * (rho_l * (rho_l - rho_v) * g * h_fg * k_l^3 / (mu_l * del_T * Lcond))^(1/4);
end