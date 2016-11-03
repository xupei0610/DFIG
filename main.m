%%
% This is the Matlab file for [EE5705] Project 3
%
% This file would automatically load the simulink to do simulation.

clear; clc;
j = sqrt(-1);

%% Parameters of DFIG
f_rated = 60;         % Rated Frequency; Unit: Hz
w_syn = 2*pi*f_rated; % Synchronous Electrical Rotation Speed; Unit: rad/s

V_ll_rated = 690;     % Rated Line Voltage; RMS Value; Unit: V

p = 6;                % # of Poles

s = 0.01;             % Slip at rated (full) load

J = 70;               % Moment of Inertia; Unit: kg*m^2

R_s = 2e-3;           % Subscript s For stator; Unit: Ohm
R_r = 1.5e-3;         % Subscript r For rotor;  Unit: Ohm

X_ls = 50e-3;         % Subscript l For leakage; Unit: Ohm
X_lr = 47e-3;
X_m = 860e-3;         % Subscript m For magnetizing; Unit: Ohm

L_ls = X_ls/w_syn;    % Unit: H
L_lr = X_lr/w_syn;
L_m  = X_m /w_syn;
L_s  = L_ls+L_m;
L_r  = L_lr+L_m;

tau_r = L_r/R_r;           % Time Constant referring to Rotor Winding
sigma = 1-L_m^2/(L_s*L_r); % Leakage Factor

%% Initial (Rated) Condition of DFIG

% Rotor Rotation Speed at rated (full) load
w_mech_rated = (1-s)*w_syn/(p/2); 

% Rated Stator Current, RMS value, Unit: A
I_s_rated = V_ll_rated/sqrt(3) / (R_s + j*X_ls + j*X_m*(R_r/s+j*X_lr)/(j*X_m+R_r/s+j*X_lr));

% Rated Rotor Current, RMS value, Unit: A
I_r_rated = -I_s_rated*j*X_m/(j*X_m+R_r/s+j*X_lr);

% Rated Torque, P_rated/w_mech_rated
T_em_rated = (3*abs(I_r_rated)^2*R_r*(1-s)/s) / w_mech_rated;
T_load_rated = T_em_rated;

% Voltage in dq domain
V_sd_rated = V_ll_rated;
V_rd_rated = 0;
V_rq_rated = 0;

% Stator Current in dq domain
I_sd_rated = sqrt(3)*real(I_s_rated);
I_sq_rated = sqrt(3)*imag(I_s_rated);

% Rotor Current in dq domain
I_rd_rated = sqrt(3)*real(I_r_rated);
I_rq_rated = sqrt(3)*imag(I_r_rated);

% Rated Qs
Q_s_rated = (V_sd_rated)^2/(w_syn*L_s)+(L_m/L_s)*V_sd_rated*I_rq_rated;

% Rated flux
fl_sd_rated = L_s*I_sd_rated + L_m*I_rd_rated;
fl_sq_rated = L_s*I_sq_rated + L_m*I_rq_rated;
fl_rd_rated = L_m*I_sd_rated + L_r*I_rd_rated;
fl_rq_rated = L_m*I_sq_rated + L_r*I_rq_rated;

%% Parameters for DFIG Controller
w_c_i = 10*2*pi; % crossover frequency
PM_i = pi/3;     % PM 

kpi_by_kii = (1/w_c_i)*tan(-pi/2+PM_i+atan(sigma*L_r*w_c_i/R_r));
kii = w_c_i*sqrt(R_r^2+(sigma*L_r*w_c_i)^2)/sqrt(1+(kpi_by_kii*w_c_i)^2);
kpi = kpi_by_kii * kii;

V_sd_prime_rated = V_rd_rated + s*w_syn*sigma*L_r*I_rq_rated;
V_sq_prime_rated = V_rq_rated - s*w_syn*sigma*L_r*I_rd_rated;

%% Parameters of Wind Turbine

A = 3904;            % Swept Area, Unit: m^2
R = 70.5/2;          % Rotor Radius, Unit: m
J_turb = 2.4*10^6;   % Moment of Inertia, Unit: kg * m^2
rho = 1.2;           % Density of air

v_wind = [12; 9; 6]; % Wind Speed, Unit: m/s

% Cp w.r.t lambda (Tip-Speed Radio) with beta = 0
% The expression of Cp is derived from the given Cp Model
Cp_2 = @(lambda) 0.52*(116*(1-0.035*lambda)/lambda-5)*exp(-21*(1-0.035*lambda)/lambda)+0.0068*lambda;

% Seek Cp_opt and lambda_opt
Cp_opt = -inf;
lambda_opt = -inf;
iter = 1;
while true
    lam = iter*0.01;
    res = Cp_2(lam);
    if ~isnan(res)
        if res < Cp_opt
            break;
        else
            lambda_opt = lam;
            Cp_opt  = res;
        end
    end
    iter = iter + 1;
end

% Calculate K_opt for Speed-Squared Controller
K_opt = Cp_opt * 0.5 * rho * A * R^3/lambda_opt^3;

% Calculate w_mech_opt for v_wind = 12m/s
w_mech_opt = lambda_opt*v_wind(1)/R;

% Calculate Gear Ratio for v_wind = 12m/s
GearRatio = w_mech_rated / w_mech_opt;


%% Simulation Part 1
fprintf('\nStarting Simulation Part 1.\n');
sim('Proj3_P1');
fprintf('Simulation Part 1 Completes.\n');

figure(1);
subplot(2, 1, 1);
plot(Ps(:,1), Ps(:,2));
hold on;
plot(Ps(:,1), Ps(:,3));
legend('Actual Value of P_{s}', 'Reference Signal of P_{s}');
ylabel('P_{s} & P_{s}^{*} [W]')
title('Control Effect of I_{rd} for Simulation 1');
subplot(2, 1, 2);
plot(Ird(:,1), Ird(:,2));
hold on;
plot(Ird(:,1), Ird(:,3));
legend('Actual Value of I_{rd}', 'Reference Signal of I_{rd}', 'Location', 'Southeast');
xlabel('Time [second]');
ylabel('I_{rd} & I_{rd}^{*} [A]');

figure(2);
subplot(2, 1, 1);
plot(Qs(:,1), Qs(:,2));
hold on;
plot(Qs(:,1), Qs(:,3));
legend('Actual Value of Q_{s}', 'Reference Signal of Q_{s}');
ylabel('Q_{s} & Q_{s}^{*} [Var]')
title('Control Effect of I_{rq} for Simulation 1');
subplot(2, 1, 2);
plot(Irq(:,1), Irq(:,2));
hold on;
plot(Irq(:,1), Irq(:,3));
legend('Actual Value of I_{rq}', 'Reference Signal of I_{rq}');
xlabel('Time [second]');
ylabel('I_{rq} & I_{rq}^{*} [A]');

figure(3);
plot(Tem(:,1), Tem(:,2));
hold on;
plot(Tem(:,1), Tem(:,3));
plot(Tem(:,1), Tem(:,4));
legend('Actural Value of T_{em}', 'Desired Reference Signal of T_{em}', 'Actual Reference Signal of T_{em}');
title('T_{em} and its Reference Signal generated for Simulation 1');
xlabel('Time [second]');
ylabel('T_{em} & T_{em}^{*} [N \cdot m]');

clear Ps Qs Pr Qr Ird Irq Tem;

%% Simulation Part 2
% This Part is the same to HW5

%% Simulation Part 3

ControlQs = -1;
ControlVwind = 1;
fprintf('\nStarting Simulation Part 3.\n');
fprintf('This Simulation May Take a while.\n')
sim('Proj3_P3');
fprintf('Simulation Part 3 Completes.\n');

figure(4);
subplot(3, 1, 1);
plot(Qs(:,1), Qs(:,2));
hold on;
plot(Qs(:,1), Qs(:,3));
legend('Actual Value of Q_{s}', 'Reference Signal of Q_{s}', 'Location', 'Southeast');
ylabel('Q_{s} & Q_{s}^{*} [Var]')
title('Control Effect of I_{rq} for Simulation 3');
subplot(3, 1, 2);
plot(Irq(:,1), Irq(:,2));
hold on;
plot(Irq(:,1), Irq(:,3));
legend('Actual Value of I_{rq}', 'Reference Signal of I_{rq}', 'Location', 'Southeast');
ylabel('I_{rq} & I_{rq}^{*} [A]');
subplot(3, 1, 3);
plot(Ird(:,1), Ird(:,2));
hold on;
plot(Ird(:,1), Ird(:,3));
legend('Actual Value of I_{rd}', 'Reference Signal of I_{rd}', 'Location', 'Southeast');
xlabel('Time [second]');
ylabel('I_{rd} & I_{rd}^{*} [A]');

figure(5);
subplot(3, 2, 1);
plot(Ps(:,1), Ps(:,2));
hold on;
plot(Ps(:,1), Ps(:,3));
legend('Actual Value of P_{s}', 'Reference Signal of P_{s}');
ylabel('P_{s} & P_{s}^{*} [W]');
xlim([0 90]);
title('Real Power for Simulation 3');
subplot(3, 2, 2);
plot(Qs(:,1), Qs(:,2));
hold on;
plot(Qs(:,1), Qs(:,3));
legend('Actual Value of Q_{s}', 'Reference Signal of Q_{s}', 'Location', 'Northwest');
ylabel('Q_{s} & Q_{s}^{*} [Var]');
xlim([0 90]);
title('Reactive Power for Simulation 3');
subplot(3, 2, 3);
plot(Pr(:,1), Pr(:,2));
legend('Actual Value of P_{r}', 'Location', 'Northwest');
ylabel('P_{r} [W]');
xlim([0 90]);
subplot(3, 2, 4);
plot(Qr(:,1), Qr(:,2));
legend('Actual Value of Q_{r}');
ylabel('Q_{r} [Var]');
xlim([0 90]);
subplot(3, 2, 5);
plot(Pr(:,1), Ps(:,2)+Pr(:,2));
legend('P = P_{s}+P_{r}', 'P_{out} = P_{turb}', 'Location', 'Northwest');
xlabel('Time [second]');
ylabel('P = P_{s}+P_{r} [W]');
xlim([0 90]);
subplot(3, 2, 6);
plot(Qr(:,1), Qs(:,2)+Qr(:,2));
legend('Q = Q_{s}+Q_{r}');
xlabel('Time [second]');
ylabel('Q = Q_{s}+Q_{r} [Var]');
xlim([0 90]);

figure(6);
plot(Tem(:,1), Tem(:,2));
hold on;
plot(Tem(:,1), Tem(:,3));
plot(Tem(:,1), Tem(:,4));
legend('Actural Value of T_{em}', 'Desired Reference Signal of T_{em}', 'Actual Reference Signal of T_{em}');
title('T_{em} and its Reference Signal generated for Simulation 3');
xlabel('Time [second]');
ylabel('T_{em} & T_{em}^{*} [N \cdot m]');

figure(7);
plot(Pwt(:,1), Pwt(:,2));
hold on;
plot(Pwt(:,1), Pwt(:,3));
legend('P_{out} = P_{turb}', 'P_{in} = P_{wind}');
xlabel('Time [second]');
ylabel('Power [W]');
title('P_{out} and P_{in} of the Wind Turbine for Simulation 3');

figure(8);
plot(Wgear(:,1), Wgear(:,2));
hold on;
plot(Wgear(:,1), Wgear(:,3));
legend('\omega_{mech}', '\omega_{mech}^{opt}');
xlabel('Time [second]');
ylabel('\omega_{mech} [rad/s]');
title('Mechanical Rotational Speed of the Wind Turbine for Simulation 3');

figure(9);
plot(Wmech(:,1), Wmech(:,2));
hold on;
plot(Wmech(:,1), Wmech(:,3));
legend('\omega_{mech}', '\omega_{mech}^{rated}');
xlabel('Time [second]');
ylabel('\omega_{mech} [rad/s]');
title('Mechanical Rotational Speed of the DFIG for Simulation 3');

clear Ird Irq Ps Qs Pr Qr Wmech Wgear Tem Pwt;

%% Simulation Part 4

ControlQs = 1;
ControlVwind = -1;
fprintf('\nStarting Simulation Part 4.\n');
fprintf('This Simulation May Take a while.\n')
sim('Proj3_P3');
fprintf('Simulation Part 4 Completes.\n');

figure(10);
subplot(3, 1, 1);
plot(Vwind(:, 1), Vwind(:, 2));
ylabel('v_{wind} [m/s]')
subplot(3, 1, 2);
plot(Ird(:,1), Ird(:,2));
hold on;
plot(Ird(:,1), Ird(:,3));
legend('Actual Value of I_{rd}', 'Reference Signal of I_{rd}', 'Location', 'Southeast');
ylabel('I_{rd} & I_{rd}^{*} [A]');
subplot(3, 1, 3);
plot(Irq(:,1), Irq(:,2));
hold on;
plot(Irq(:,1), Irq(:,3));
legend('Actual Value of I_{rq}', 'Reference Signal of I_{rq}', 'Location', 'Southeast');
xlabel('Time [second]');
ylabel('I_{rq} & I_{rq}^{*} [A]');

figure(11);
subplot(3, 2, 1);
plot(Ps(:,1), Ps(:,2));
hold on;
plot(Ps(:,1), Ps(:,3));
legend('Actual Value of P_{s}', 'Reference Signal of P_{s}');
ylabel('P_{s} & P_{s}^{*} [W]');
xlim([0 90]);
title('Real Power for Simulation 4');
subplot(3, 2, 2);
plot(Qs(:,1), Qs(:,2));
hold on;
plot(Qs(:,1), Qs(:,3));
legend('Actual Value of Q_{s}', 'Reference Signal of Q_{s}', 'Location', 'Southeast');
ylabel('Q_{s} & Q_{s}^{*} [Var]');
xlim([0 90]);
title('Reactive Power for Simulation 4');
subplot(3, 2, 3);
plot(Pr(:,1), Pr(:,2));
legend('Actual Value of P_{r}');
ylabel('P_{r} [W]');
xlim([0 90]);
subplot(3, 2, 4);
plot(Qr(:,1), Qr(:,2));
legend('Actual Value of Q_{r}', 'Location', 'Southeast');
ylabel('Q_{r} [Var]');
xlim([0 90]);
subplot(3, 2, 5);
plot(Pr(:,1), Ps(:,2)+Pr(:,2));
hold on;
plot(Pwt(:,1), Pwt(:,2));
legend('P = P_{s}+P_{r}', 'P_{out} = P_{turb}');
xlabel('Time [second]');
ylabel('P = P_{s}+P_{r} [W]');
xlim([0 90]);
subplot(3, 2, 6);
plot(Qr(:,1), Qs(:,2)+Qr(:,2));
legend('Q = Q_{s}+Q_{r}', 'Location', 'Southeast');
xlabel('Time [second]');
ylabel('Q = Q_{s}+Q_{r} [Var]');
xlim([0 90]);

figure(12);
plot(Tem(:,1), Tem(:,2));
hold on;
plot(Tem(:,1), Tem(:,3));
plot(Tem(:,1), Tem(:,4));
legend('Actural Value of T_{em}', 'Desired Reference Signal of T_{em}', 'Actual Reference Signal of T_{em}');
title('T_{em} and its Reference Signal generated for Simulation 4');
xlabel('Time [second]');
ylabel('T_{em} & T_{em}^{*} [N \cdot m]');

figure(13);
plot(Pwt(:,1), Pwt(:,2));
hold on;
plot(Pwt(:,1), Pwt(:,3));
legend('P_{out} = P_{turb}', 'P_{in} = P_{wind}');
xlabel('Time [second]');
ylabel('Power [W]');
title('P_{out} and P_{in} of the Wind Turbine for Simulation 4');

figure(14);
plot(Wgear(:,1), Wgear(:,2));
hold on;
plot(Wgear(:,1), Wgear(:,3));
legend('\omega_{mech}', '\omega_{mech}^{opt}');
xlabel('Time [second]');
ylabel('\omega_{mech} [rad/s]');
title('Mechanical Rotational Speed of the Wind Turbine for Simulation 4');

figure(15);
plot(Wmech(:,1), Wmech(:,2));
hold on;
plot(Wmech(:,1), Wmech(:,3));
legend('\omega_{mech}', '\omega_{mech}^{rated}');
xlabel('Time [second]');
ylabel('\omega_{mech} [rad/s]');
title('Mechanical Rotational Speed of the DFIG for Simulation 4');

clear Ird Irq Ps Qs Pr Qr Wmech Wgear Tem Pwt;