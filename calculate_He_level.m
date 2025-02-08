function function [new_level, P_new, T_new] = calculate_He_level(T_thermo, current_level, P_old, T_old, ...
    del_t, A_thermo, V_total, m_total, rho_He_interp, rho_v_interp, h_fg_interp, ...
    P_data, rho_v_data, T_data, T_surface)


    % 1. 평균 온도 계산
    T_evap_surface = T_thermo(end); % 증발면(써모사이펀 마지막 노드) 온도

    % 2. 현재 포화 온도 계산 (현재 압력에서)
    if isempty(P_old) || isnan(P_old)
        P_old = interp1(T_data, P_data, T_old, 'linear', 'extrap');
        warning('⚠ P_old was empty. Set to interpolated value: %.3f kPa', P_old);
    end
    T_sat = interp1(P_data, T_data, P_old, 'linear', 'extrap');

     % 3. 핵 비등 열전달 계수 계산
    % 핵비등 열전달 계수 계산 (h_He_NB 적용)
    h_boiling = h_He_NB(current_pressure, T_m1(1), properties);

    % 4. GGG와의 온도 차이를 이용한 냉각량 계산
    Q_cooling = h_boiling * A_thermo * (T_evap_surface - T_sat) * del_t;

    % 5. 증발량 계산 (Q_cooling / 증발잠열)
    % 동적 증발잠열 계산 (엑셀 데이터 기반)
    L_He = interp1(T_data, H_vapor_data, T_sat) - interp1(T_data, H_liquid_data, T_sat);
    % 증발된 헬륨 질량 계산 (핵비등 열전달 적용)
    m_evap = Q_boiling / L_He;
    % 증발된 헬륨 질량을 현재 총 헬륨 질량에서 제거
    m_He = m_He - m_evap;
    % 남은 헬륨 질량을 이용하여 헬륨 레벨 업데이트
    He_level = m_He / (rho_He_interp(T_thermo(end)) * A_thermo);

    % 6. 새로운 헬륨 액체 수위 계산
    new_level = max(0, current_level - (m_evap / rho_He_interp(T_sat)));

    % 7. 새로운 액체 및 기체 부피 계산
    V_liquid_new = A_thermo * new_level;
    V_gas_new = V_total - V_liquid_new;

    % 8. 연립 방정식 풀기 (새로운 밀도 계산)
    A_matrix = [V_gas_new, V_liquid_new; 1, 1];  
    B_vector = [m_total; V_total]; 

    rho_values = A_matrix \ B_vector;  
    rho_gas_new = rho_values(1);
    rho_liquid_new = rho_values(2);

    % 9. 새로운 밀도로부터 압력 및 온도 계산
    P_new = interp1(rho_v_data, P_data, rho_gas_new, 'linear', 'extrap');
    T_new = interp1(P_data, T_data, P_new, 'linear', 'extrap');
end
