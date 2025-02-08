function P = calculate_system_pressure(T_vapor, He_level, L_thermo, rho_He_interp)
    % 물리 상수
    g = 9.81;    % 중력가속도 [m/s²]
    
    % 수두압 계산
    T_evap_surface = T_vapor; % 증발면(써모사이펀 마지막 노드) 온도
    rho_liquid = rho_He_interp(T_evap_surface); % 증발면 온도에서의 액체 헬륨 밀도
    P_hydrostatic = rho_liquid * g * He_level;  % 수두압 [Pa]
    
    % 증기압 계산 - 헬륨의 증기압 곡선
    % 아래 상수들은 헬륨의 실제 데이터에 맞게 조정 필요
    A = 5.6;  % [ln(kPa)]
    B = 2.858; % [K]
    P_vapor = exp(A - B/T_vapor) * 1000;  % [Pa]
    
    % 총 압력 계산 및 kPa 단위로 변환
    P = (P_hydrostatic + P_vapor) / 1000;
end