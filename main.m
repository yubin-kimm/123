% All units are in SI
% Author: yubin-kimm
% Last modified: 2025-02-07 02:54:16 UTC

%% 초기 설정 및 물성치 정의
Bmax = 2.5;                % 최대 자기장 [T]
Magnetization_time = 10;   % 자화 과정에 소요되는 시간 [s]
Cooling_time = 10;         % 냉각 과정에 소요되는 시간 [s]
Demagnetization_time = 10; % 탈자화 과정에 소요되는 시간 [s]
Holding_time = 10;         % 대기 과정에 소요되는 시간 [s]
T_init = 4.2;             % 시스템의 초기 온도 [K]

%% 헬륨 물성치 데이터 읽기
opts = detectImportOptions('properties.xlsx');
opts.VariableNamingRule = 'preserve';  % 원래 열 이름 보존
properties = readtable('properties.xlsx', opts);

% 데이터 추출
T_data = properties.Temperature;
rho_He_data = properties{:, 3};  % Liquid Density 열 위치
k_He_data = properties{:, 11};   % Liquid Thermal Conductivity 열 위치
cp_He_data = properties{:, 9};   % Liquid Cp 열 위치
mu_He_data = properties{:, 13};  % Liquid Viscosity 열 위치

% 보간 함수 생성
rho_He_interp = @(T) interp1(T_data, rho_He_data, T, 'linear', 'extrap');
k_He_interp = @(T) interp1(T_data, k_He_data, T, 'linear', 'extrap');
cp_He_interp = @(T) interp1(T_data, cp_He_data, T, 'linear', 'extrap');
mu_He_interp = @(T) interp1(T_data, mu_He_data, T, 'linear', 'extrap');

%% 써모사이펀 물성치 및 구조 정의
% 구조 파라미터
L_thermo = 0.05;          % 써모사이펀 총 길이 [m] (5cm)
D_thermo = 0.01;          % 써모사이펀 직경 [m] (1cm)
A_thermo = pi * (D_thermo/2)^2; % 써모사이펀 단면적 [m²]

% 써모사이펀 격자 설정
N_thermo = 40;           % 써모사이펀 격자 수
del_x_thermo = L_thermo / N_thermo; % 격자 간격 [m]

% 초기 조건
T_thermo = T_init * ones(N_thermo, 1); % 써모사이펀 초기 온도 분포
v_thermo = zeros(N_thermo, 1);         % 초기 유체 속도 분포

% 초기 물성치 계산
rho_He = rho_He_interp(T_init);
k_He = k_He_interp(T_init);
cp_He = cp_He_interp(T_init);
mu_He = mu_He_interp(T_init);

%% 액체 헬륨 레벨 모니터링을 위한 변수
initial_He_level = 0.1 * L_thermo;  % 초기 액체 헬륨 레벨 (써모사이펀 길이의 10%)
He_level = initial_He_level;        % 현재 액체 헬륨 레벨
He_level_history_mag = zeros(1, N_t_mag);    % 자화 과정 동안의 레벨 기록
He_level_history_cool = zeros(1, N_t_cool);  % 냉각 과정 동안의 레벨 기록
He_level_history_demag = zeros(1, N_t_demag); % 탈자화 과정 동안의 레벨 기록
He_level_history_hold = zeros(1, N_t_hold);   % 대기 과정 동안의 레벨 기록

%% GGG 물성치 및 구조 정의
rho_GGG = 7080;          % GGG(가드올리늄 갈륨 가넷)의 밀도 [kg/m³]
A_GGG = 1.257e-3;        % GGG 단면적 [m²]
L_GGG = 0.07;            % GGG 길이 [m]

N = 40;                  % 격자 수, GGG 구간을 나누는 개수
del_x_GGG = L_GGG / N;   % 각 격자의 공간 간격 [m]

%% 시간 설정
del_t_mag = Magnetization_time / (Magnetization_time * 5e5); % 자화 단계의 시간 간격
N_t_mag = Magnetization_time / del_t_mag; % 자화 단계의 시간 스텝 개수

del_t_cool = Cooling_time / (Cooling_time * 5e5); % 냉각 단계의 시간 간격
N_t_cool = Cooling_time / del_t_cool; % 냉각 단계의 시간 스텝 개수

del_t_demag = Demagnetization_time / (Demagnetization_time * 5e5); % 탈자화 단계의 시간 간격
N_t_demag = Demagnetization_time / del_t_demag; % 탈자화 단계의 시간 스텝 개수

del_t_hold = Holding_time / (Holding_time * 1e4); % 대기 단계의 시간 간격
N_t_hold = Holding_time / del_t_hold; % 대기 단계의 시간 스텝 개수

%% 자기장(B-field) 설정
Bfield_mag = linspace(0, Bmax, N_t_mag);      % 자화 단계 동안 선형적으로 증가하는 자기장 값
Bfield_demag = linspace(Bmax, 0, N_t_demag);  % 탈자화 단계 동안 선형적으로 감소하는 자기장 값

%% 초기 조건 설정
T_m1 = T_init * ones(N, 1);      % GGG의 초기 온도 분포
T_avg_mag = zeros(1, N_t_mag);   % 자화 단계 온도 저장 배열
T_avg_mag(1) = mean(T_m1);       % 첫 스텝의 평균 온도 초기화

%% GGG 질량 계산
mass_GGG = rho_GGG * A_GGG * L_GGG; % GGG의 전체 질량 계산 [kg]

%% 액체 헬륨 레벨 계산 함수

function new_level = calculate_He_level(current_level, T_avg, del_t, base_evap_rate, temp_evap_coeff)

    % 기본 증발
    level_change = 0;
    
    % 온도에 따른 추가 증발 (4.2K 이상일 때)
    if T_avg > 4.2
        level_change = level_change + temp_evap_coeff * (T_avg - 4.2) * del_t;
    end
    
    % 새로운 레벨 계산
    new_level = current_level;
    
    % 레벨이 0 이하로 내려가지 않도록 제한
    if new_level < 0
        new_level = 0;
    end
end

%% 헬륨 레벨 경고 함수
function check_He_level_warning(current_level, total_length)
    warning_threshold = 0.02; % 20% 이하일 때 경고
    critical_threshold = 0.01; % 10% 이하일 때 심각
    
    level_percentage = current_level / total_length * 100;
    
    if level_percentage <= critical_threshold * 100
        warning('CRITICAL: Helium level is critically low (%.1f%%)!', level_percentage);
    elseif level_percentage <= warning_threshold * 100
        warning('WARNING: Helium level is low (%.1f%%)!', level_percentage);
    end
end


%% 자화(Magnetization) 단계 시뮬레이션
fprintf('Starting Magnetization process...\n');
for j = 2:N_t_mag
    T_avg_mag(j) = mean(T_m1);
    dT = calculate_dT(T_avg_mag(j), Bfield_mag(j), Bfield_mag(j-1), mass_GGG);
    [T_m1, T_thermo] = update_temperature(T_m1, T_thermo, dT, del_t_mag, del_x_GGG, ...
        del_x_thermo, rho_GGG, rho_He_interp, k_He_interp, cp_He_interp, A_GGG, A_thermo, D_thermo);
    
    % 액체 헬륨 레벨 업데이트
    He_level = calculate_He_level(He_level, mean(T_thermo), del_t_mag, He_evaporation_rate, He_temp_dependent_evap);
    He_level_history_mag(j) = He_level;
    
    % 레벨 경고 체크
    check_He_level_warning(He_level, L_thermo);
    
    if mod(j, 1000) == 0
        fprintf('Magnetization - Time: %.2f s, B-field: %.2f T, Avg Temp: %.4f K, He Level: %.2f%%\n', ...
            j*del_t_mag, Bfield_mag(j), T_avg_mag(j), (He_level/L_thermo)*100);
    end
end

%% 냉각(Cooling) 단계 시뮬레이션
fprintf('\nStarting Cooling process...\n');
T_avg_cool = zeros(1, N_t_cool);
T_avg_cool(1) = T_avg_mag(end);
for j = 2:N_t_cool
    [T_m1, T_thermo] = update_temperature(T_m1, T_thermo, 0, del_t_cool, del_x_GGG, ...
        del_x_thermo, rho_GGG, rho_He_interp, k_He_interp, cp_He_interp, A_GGG, A_thermo, D_thermo);
    T_avg_cool(j) = mean(T_m1);
    
    % 액체 헬륨 레벨 업데이트
    He_level = calculate_He_level(He_level, mean(T_thermo), del_t_cool, He_evaporation_rate, He_temp_dependent_evap);
    He_level_history_cool(j) = He_level;
    
    % 레벨 경고 체크
    check_He_level_warning(He_level, L_thermo);
    
    if mod(j, 1000) == 0
        fprintf('Cooling - Time: %.2f s, Avg Temp: %.4f K, He Level: %.2f%%\n', ...
            (Magnetization_time + j*del_t_cool), T_avg_cool(j), (He_level/L_thermo)*100);
    end
end

%% 탈자화(Demagnetization) 단계 시뮬레이션
fprintf('\nStarting Demagnetization process...\n');
T_avg_demag = zeros(1, N_t_demag);
T_avg_demag(1) = T_avg_cool(end);
for j = 2:N_t_demag
    T_avg_demag(j) = mean(T_m1);
    dT = calculate_dT(T_avg_demag(j), Bfield_demag(j), Bfield_demag(j-1), mass_GGG);
    [T_m1, T_thermo] = update_temperature(T_m1, T_thermo, dT, del_t_demag, del_x_GGG, ...
        del_x_thermo, rho_GGG, rho_He_interp, k_He_interp, cp_He_interp, A_GGG, A_thermo, D_thermo);
    
    % 액체 헬륨 레벨 업데이트
    He_level = calculate_He_level(He_level, mean(T_thermo), del_t_demag, He_evaporation_rate, He_temp_dependent_evap);
    He_level_history_demag(j) = He_level;
    
    % 레벨 경고 체크
    check_He_level_warning(He_level, L_thermo);
    
    if mod(j, 1000) == 0
        fprintf('Demagnetization - Time: %.2f s, B-field: %.2f T, Avg Temp: %.4f K, He Level: %.2f%%\n', ...
            (Magnetization_time + Cooling_time + j*del_t_demag), Bfield_demag(j), T_avg_demag(j), (He_level/L_thermo)*100);
    end
end

%% 대기(Holding) 단계 시뮬레이션
fprintf('\nStarting Holding process...\n');
T_avg_hold = zeros(1, N_t_hold);
T_avg_hold(1) = T_avg_demag(end);
for j = 2:N_t_hold
    [T_m1, T_thermo] = update_temperature(T_m1, T_thermo, 0, del_t_hold, del_x_GGG, ...
        del_x_thermo, rho_GGG, rho_He_interp, k_He_interp, cp_He_interp, A_GGG, A_thermo, D_thermo);
    T_avg_hold(j) = mean(T_m1);
    
    % 액체 헬륨 레벨 업데이트
    He_level = calculate_He_level(He_level, mean(T_thermo), del_t_hold, He_evaporation_rate, He_temp_dependent_evap);
    He_level_history_hold(j) = He_level;
    
    % 레벨 경고 체크
    check_He_level_warning(He_level, L_thermo);
    
    if mod(j, 1000) == 0
        fprintf('Holding - Time: %.2f s, Avg Temp: %.4f K, He Level: %.2f%%\n', ...
            (Magnetization_time + Cooling_time + Demagnetization_time + j*del_t_hold), T_avg_hold(j), (He_level/L_thermo)*100);
    end
end

% All units are in SI
% Author: yubin-kimm
% Last modified: 2025-02-07 02:57:42 UTC

%% 결과 출력
fprintf('\nSimulation completed successfully.\n');

%% 그래프 출력
% 시간 벡터 생성
time_vector_mag = linspace(0, Magnetization_time, N_t_mag);
time_vector_cool = linspace(Magnetization_time, Magnetization_time + Cooling_time, N_t_cool);
time_vector_demag = linspace(Magnetization_time + Cooling_time, ...
    Magnetization_time + Cooling_time + Demagnetization_time, N_t_demag);
time_vector_hold = linspace(Magnetization_time + Cooling_time + Demagnetization_time, ...
    Magnetization_time + Cooling_time + Demagnetization_time + Holding_time, N_t_hold);

% 그래프 생성을 위한 figure 설정
figure('Position', [100, 100, 1200, 1000]);

% 서브플롯 1: 평균 온도 변화
subplot(3,1,1);
plot(time_vector_mag, T_avg_mag, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Magnetization'); hold on;
plot(time_vector_cool, T_avg_cool, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Cooling');
plot(time_vector_demag, T_avg_demag, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Demagnetization');
plot(time_vector_hold, T_avg_hold, 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Holding');
xlabel('Time (s)');
ylabel('Average Temperature (K)');
title('Average Temperature vs. Time');
legend('show');
grid on;

% 서브플롯 2: 액체 헬륨 레벨
subplot(3,1,2);
plot(time_vector_mag, (He_level_history_mag/L_thermo)*100, 'LineWidth', 2, 'Color', 'b', 'DisplayName', 'Magnetization'); hold on;
plot(time_vector_cool, (He_level_history_cool/L_thermo)*100, 'LineWidth', 2, 'Color', 'g', 'DisplayName', 'Cooling');
plot(time_vector_demag, (He_level_history_demag/L_thermo)*100, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Demagnetization');
plot(time_vector_hold, (He_level_history_hold/L_thermo)*100, 'LineWidth', 2, 'Color', 'k', 'DisplayName', 'Holding');
xlabel('Time (s)');
ylabel('Liquid Helium Level (%)');
title('Liquid Helium Level vs. Time');
yline(20, '--r', 'Warning Level (20%)', 'LineWidth', 1.5);
yline(10, '--r', 'Critical Level (10%)', 'LineWidth', 1.5);
legend('show');
grid on;

% 서브플롯 3: 최종 온도 분포
subplot(3,1,3);
% 온도 분포 플롯
yyaxis left
plot(linspace(0, L_thermo, N_thermo), T_thermo, 'b-', 'LineWidth', 2, 'DisplayName', 'Thermosyphon Temperature');
hold on;
plot(L_thermo + linspace(0, L_GGG, N), T_m1, 'r-', 'LineWidth', 2, 'DisplayName', 'GGG Temperature');
ylabel('Temperature (K)');

% 액체 헬륨 레벨 표시
yyaxis right
area([0 L_thermo], [100 100], 'FaceColor', [0.8 0.8 1], 'FaceAlpha', 0.3, 'DisplayName', 'He Level');
ylim([0 100]);
plot([0 L_thermo], [He_level/L_thermo*100 He_level/L_thermo*100], '--b', 'LineWidth', 1.5, 'DisplayName', 'Current He Level');
ylabel('Liquid Helium Level (%)');

xlabel('Position (m)');
title('Final Temperature Distribution and Helium Level');
legend('show');
grid on;

%% 헬퍼 함수 정의
function dT = calculate_dT(T, B, B_prev, mass)
    dM_dT = vdMdt(T, B);  % vdMdt는 자기화율 함수
    dB = B - B_prev;      % 자기장의 변화량 계산
    dQ = -T * dM_dT * dB; % 열 생성 계산
    C_GGG_val = c_GGG(T, B); % 열용량 계산
    dT = dQ / (mass * C_GGG_val); % 온도 변화량 계산
end

function [T_new_GGG, T_new_thermo] = update_temperature(T_old_GGG, T_old_thermo, dT, del_t, ...
    del_x_GGG, del_x_thermo, rho_GGG, rho_He_interp, k_He_interp, cp_He_interp, A_GGG, A_thermo, D_thermo)
    
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

    % 경계면에서의 열전달
    h_interface = 500; % 계면 열전달 계수 [W/(m²·K)]
    interface_area = min(A_GGG, pi * D_thermo * del_x_thermo); % 접촉 면적

    % GGG와 써모사이펀 사이의 열전달
    q_interface = h_interface * interface_area * (T_old_GGG(1) - T_old_thermo(end));
    
    % 경계조건 업데이트
    T_new_GGG(1) = T_old_GGG(1) - q_interface * del_t / (rho_GGG * c_GGG(T_old_GGG(1), 0) * A_GGG * del_x_GGG);
    T_new_thermo(end) = T_old_thermo(end) + q_interface * del_t / (rho_He * cp_He * A_thermo * del_x_thermo);
    
    % 다른 경계조건
    T_new_GGG(N_GGG) = T_new_GGG(N_GGG-1);
    T_new_thermo(1) = 4.2; % 써모사이펀 상단은 액체 헬륨 온도로 고정

    % 자기열량 효과 추가
    T_new_GGG = T_new_GGG + dT;
end