function dT = calculate_dT(T, B, B_prev, mass)
    % Calculates temperature change due to magnetocaloric effect
    dM_dT = vdMdt(T, B);  % Get magnetic susceptibility
    dB = B - B_prev;      % Magnetic field change
    dQ = -T * dM_dT * dB; % Heat generation
    C_GGG_val = c_GGG(T, B); % Get specific heat capacity
    dT = dQ / C_GGG_val; % Temperature change
end