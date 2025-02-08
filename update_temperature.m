function [T_new_GGG, T_new_thermo] = update_temperature(T_old_GGG, T_old_thermo, dT, del_t, ...
    del_x_GGG, del_x_thermo, rho_GGG, rho_He_interp, k_He_interp, cp_He_interp, A_GGG, A_thermo, D_thermo, He_level, L_thermo)
    
    % GGG 부분
    N_GGG = length(T_old_GGG);
    T_new_GGG = T_old_GGG;
    alpha_GGG = k_GGG(mean(T_old_GGG)) / (rho_GGG * c_GGG(mean(T_old_GGG), 0));

    % 써모사이펀 부분
    N_thermo = length(T_old_thermo);
    T_new_thermo = T_old_thermo;
    
    % 각 격자점에서의 헬륨 물성치 계산
    T_mean_thermo = mean(T_old_thermo);
    rho_He = rho_He_interp(T_mean_thermo);
    k_He = k_He_interp(T_mean_thermo);
    cp_He = cp_He_interp(T_mean_thermo);
    alpha_He = k_He / (rho_He * cp_He);

    % 내부 격자점 업데이트 - GGG
    for i = 2:N_GGG-1
        T_new_GGG(i) = T_old_GGG(i) + alpha_GGG * del_t / (del_x_GGG^2) * ...
            (T_old_GGG(i+1) - 2*T_old_GGG(i) + T_old_GGG(i-1));
    end

    % 내부 격자점 업데이트 - 써모사이펀
    for i = 2:N_thermo-1
        T_new_thermo(i) = T_old_thermo(i) + alpha_He * del_t / (del_x_thermo^2) * ...
            (T_old_thermo(i+1) - 2*T_old_thermo(i) + T_old_thermo(i-1));
    end

    % 응축 열전달 계산
    % 증기 공간의 평균 온도 계산
    vapor_idx = ceil(He_level/del_x_thermo):N_thermo;
    if ~isempty(vapor_idx)
        T_vapor = mean(T_old_thermo(vapor_idx));
    else
        T_vapor = T_old_thermo(end);
    end
    
    % 응축 열전달 계수 계산
    h_interface = h_He_FC(calculate_system_pressure(T_vapor, He_level, L_thermo, rho_He_interp), ...
                         T_old_GGG(1), ...
                         L_thermo - He_level); % 응축 특성 길이는 증기 공간 높이
    
    interface_area = min(A_GGG, pi * D_thermo * del_x_thermo); % 접촉 면적

    % GGG와 써모사이펀 사이의 열전달
    q_interface = h_interface * interface_area * (T_old_GGG(1) - T_old_thermo(end));
    
    % 경계조건 업데이트
    % GGG 첫 번째 노드와 써모사이펀 마지막 노드의 열전달 수정
    T_new_GGG(1) = T_old_GGG(1) - q_interface * del_t / (rho_GGG * c_GGG(T_old_GGG(1), 0) * A_GGG * del_x_GGG);
    T_new_thermo(end) = T_old_thermo(end) + q_interface * del_t / (rho_He * cp_He * A_thermo * del_x_thermo);
    T_new_thermo(end) = T_old_thermo(end) + q_interface * del_t / (rho_He * cp_He * A_thermo * del_x_thermo);
    
    % 다른 경계조건
    T_new_GGG(N_GGG) = T_new_GGG(N_GGG-1);
    T_new_thermo(1) = 4.2; % 써모사이펀 상단은 액체 헬륨 온도로 고정

    % 자기열량 효과 추가
    T_new_GGG = T_new_GGG + dT;
end