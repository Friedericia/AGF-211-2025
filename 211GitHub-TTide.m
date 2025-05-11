%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data folder and subfolders must be added to path before running
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Read data

Minilog5m = readtable('Minilog-II-T_355571_20250405_1.csv','NumHeaderLines',8);
Minilog15m = readtable('Minilog-II-T_355583_20250405_1.csv','NumHeaderLines',8);
Minilog25m = readtable('Minilog-II-T_355572_20250405_1.csv','NumHeaderLines',8);
Minilog35m = readtable('Minilog-II-T_355573_20250405_1.csv','NumHeaderLines',8);
Minilog45m = readtable('Minilog-II-T_358946_20250405_1.csv','NumHeaderLines',8);

Seabird1m = readtable('SBE37SMP-RS232_03720169_2025_04_05.csv');
Seabird20m = readtable('SBE37SMP-RS232_03727333_2025_04_05.csv');
Seabird30m = readtable('SBE37SM-RS232_03710965_2025_04_05.csv');

Seaguard10m = readtable('1711.txt');
Seaguard40m = readtable('1705.txt');

ADCP = readtable('adcp.txt');

%% Plot Current -  Seaguard
figure(1)
% Current Direction vs Time - 10m
figure(1)
subplot(4,2,1)
plot(Seaguard10m.Var2(16:end-7), Seaguard10m.Var8(16:end-7), 'LineWidth', 1.5, 'Color', "#0072BD");  
ylabel('Direction °', 'FontSize', 12)
ylim([0,360]);
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);  
xlabel('Time', 'FontSize', 12)
grid on;
title('Current Direction vs Time (10m)', 'FontSize', 14)

% Current Speed vs Time - 10m
subplot(4,2,3)
plot(Seaguard10m.Var2(16:end-7), Seaguard10m.Var7(16:end-7), 'LineWidth', 1.5, 'Color', "#7E2F8E");
ylabel('Speed cm/s', 'FontSize', 12)
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);  
xlabel('Time', 'FontSize', 12)
grid on;
title('Current Speed vs Time (10m)', 'FontSize', 14)

% North Velocity vs Time - 10m
subplot(4,2,5)
plot(Seaguard10m.Var2(16:end-7), Seaguard10m.Var9(16:end-7), 'LineWidth', 1.5, 'Color', "#77AC30");
ylabel('Velocity cm/s', 'FontSize', 12)
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);  
xlabel('Time', 'FontSize', 12)
grid on;
title('North Velocity vs Time (10m)', 'FontSize', 14)

% East Velocity vs Time - 10m
subplot(4,2,7)
plot(Seaguard10m.Var2(16:end-7), Seaguard10m.Var10(16:end-7), 'LineWidth', 1.5, 'Color', "#4DBEEE");
ylabel('Velocity cm/s', 'FontSize', 12)
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);  
xlabel('Time', 'FontSize', 12)
grid on;
title('East Veloity vs Time (10m)', 'FontSize', 14)

% 40m - Current Direction vs Time

subplot(4,2,2)
plot(Seaguard10m.Var2(16:end-7), Seaguard40m.Var14(16:end-7), 'LineWidth', 1.5, 'Color', "#0072BD");
ylabel('Direction °', 'FontSize', 12)
ylim([0,360]);
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);
xlabel('Time', 'FontSize', 12)
grid on;
title('Current Direction vs Time (40m)', 'FontSize', 14)

% Current Speed vs Time - 40m
subplot(4,2,4)
plot(Seaguard10m.Var2(16:end-7), Seaguard40m.Var13(16:end-7), 'LineWidth', 1.5, 'Color', "#7E2F8E");
ylabel('Speed cm/s', 'FontSize', 12)
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);
xlabel('Time', 'FontSize', 12)
grid on;
title('Current Speed vs Time (40m)', 'FontSize', 14)

% North Velocity vs Time - 40m
subplot(4,2,6)
plot(Seaguard10m.Var2(16:end-7), Seaguard40m.Var15(16:end-7), 'LineWidth', 1.5, 'Color', "#77AC30");
ylabel('Velocity cm/s', 'FontSize', 12)
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);
xlabel('Time', 'FontSize', 12)
grid on;
title('North Velocity vs Time (40m)', 'FontSize', 14)

% East Velocity vs Time - 40m
subplot(4,2,8)
plot(Seaguard10m.Var2(16:end-7), Seaguard40m.Var16(16:end-7), 'LineWidth', 1.5, 'Color', "#4DBEEE");
ylabel('Velocity cm/s', 'FontSize', 12)
xlim([Seaguard10m.Var2(16), Seaguard10m.Var2(end-7)]);
xlabel('Time', 'FontSize', 12)
grid on;
title('East Velocity vs Time (40m)', 'FontSize', 14)

%% New Tidal Analysis using T_Tide - Seaguard 10m


t = Seaguard10m.Var2(16:end-11);
u = Seaguard10m.Var10(16:end-11); % Eastward
v = Seaguard10m.Var9(16:end-11);  % Northward

 
t_dnum = datenum(t);


dt_hours = median(diff(t_dnum)) * 24;

 
uv = u + 1i * v;
tide_uv = t_tide(uv, 'interval', dt_hours, 'start time', t_dnum(1), 'latitude', 78.2);

% === Plot Top 6 Tidal Constituents by Major Axis ===
[sorted_major, idx_major] = sort(tide_uv.tidecon(:,3), 'descend');  % Column 3 = major axis

figure();
bar(categorical(cellstr(tide_uv.name(idx_major(1:6),:))), sorted_major(1:6));
ylabel('Major Axis Amplitude (cm/s)');
title('Top 6 Tidal Constituents (Current Ellipse Major Axis)');

% === Predict tidal current using complex output ===
uv_pred = t_predic(t_dnum, tide_uv);
 
u_pred = real(uv_pred);
v_pred = imag(uv_pred);

% === Calculate residuals ===
u_resid = u - u_pred;
v_resid = v - v_pred;

figure();

% --- Eastward component ---
subplot(2,1,1);
plot(t_dnum, u, 'k', 'DisplayName', 'Measured U');
hold on;
plot(t_dnum, u_pred, 'b', 'DisplayName', 'Tide (U)');
plot(t_dnum, u_resid, 'r', 'DisplayName', 'Residual (U)');
datetick('x', 'dd-mmm HH:MM', 'keeplimits', 'keepticks');
ylabel('U (cm/s)', 'FontSize', 18);
title('Eastward Current: Measured vs Tidal vs Residual', 'FontSize', 20);
legend('FontSize', 15);
set(gca, 'FontSize', 15);

% --- Northward component ---
subplot(2,1,2);
plot(t_dnum, v, 'k', 'DisplayName', 'Measured V');
hold on;
plot(t_dnum, v_pred, 'b', 'DisplayName', 'Tide (V)');
plot(t_dnum, v_resid, 'r', 'DisplayName', 'Residual (V)');
datetick('x', 'dd-mmm HH:MM', 'keeplimits', 'keepticks');
ylabel('V (cm/s)', 'FontSize', 18);
xlabel('Date', 'FontSize', 18);
title('Northward Current: Measured vs Tidal vs Residual', 'FontSize', 20);
legend('FontSize', 15);
set(gca, 'FontSize', 15);


% === High and low tides from zero-crossings of predicted northward current ===
isHigh = islocalmax(v_pred);
isLow  = islocalmin(v_pred);

% Zero crossings of northward predicted tide
zeroCrossings = find(diff(sign(v_pred)) ~= 0);
tideChangeTimes = t_dnum(zeroCrossings);

fprintf('Estimated High/Low Tide Times (zero crossings of v_pred):\n');
for i = 1:length(tideChangeTimes)
    fprintf('%s\n', datestr(tideChangeTimes(i), 'dd-mmm-yyyy HH:MM'));
end

maj = [0.822, 5.336, 1.073, 0.459, 0.525, 1.749, 0.354, 0.174];  
min = [-0.284, -1.540, -0.453, -0.105, -0.099, -0.253, -0.019, 0.134];  



eccentricity = sqrt(1 - (min ./ maj).^2);

 
T = table(maj', min', eccentricity', 'VariableNames', {'MajorAxis', 'MinorAxis', 'Eccentricity'});

 
disp(T);

amplitude = sqrt((maj.^2)+(min.^2));
disp(amplitude)


%% New Tidal Analysis 
%% Tidal Analysis using T_Tide - Seaguard 40m

% Extract time and velocity components
t = Seaguard40m.Var2(16:end-11);
u = Seaguard40m.Var16(16:end-11); % Eastward
v = Seaguard40m.Var15(16:end-11); % Northward

% Convert time to datenum format
t_dnum = datenum(t);

% Time interval in hours
dt_hours = median(diff(t_dnum)) * 24;

% === Use complex velocity vector ===
uv = u + 1i * v;
tide_uv = t_tide(uv, 'interval', dt_hours, 'start time', t_dnum(1), 'latitude', 78.2);

% === Plot Top 6 Tidal Constituents by Major Axis ===
[sorted_major, idx_major] = sort(tide_uv.tidecon(:,3), 'descend');  % Column 3 = major axis

figure();
bar(categorical(cellstr(tide_uv.name(idx_major(1:6),:))), sorted_major(1:6));
ylabel('Major Axis Amplitude (cm/s)');
title('Top 6 Tidal Constituents (Current Ellipse Major Axis)');

% === Predict tidal current using complex output ===
uv_pred = t_predic(t_dnum, tide_uv);

u_pred = real(uv_pred);
v_pred = imag(uv_pred);

% === Calculate residuals ===
u_resid1 = u - u_pred;
v_resid1 = v - v_pred;

% === Plot Eastward component ===
figure();
subplot(2,1,1);
plot(t_dnum, u, 'k', 'DisplayName', 'Measured U');
hold on;
plot(t_dnum, u_pred, 'b', 'DisplayName', 'Tide (U)');
plot(t_dnum, u_resid1, 'r', 'DisplayName', 'Residual (U)');
datetick('x', 'dd-mmm HH:MM', 'keeplimits', 'keepticks');
ylabel('U (cm/s)');
title('Eastward Current: Measured vs Tidal vs Residual');
legend();

% === Plot Northward component ===
subplot(2,1,2);
plot(t_dnum, v, 'k', 'DisplayName', 'Measured V');
hold on;
plot(t_dnum, v_pred, 'b', 'DisplayName', 'Tide (V)');
plot(t_dnum, v_resid1, 'r', 'DisplayName', 'Residual (V)');
datetick('x', 'dd-mmm HH:MM', 'keeplimits', 'keepticks');
ylabel('V (cm/s)');
title('Northward Current: Measured vs Tidal vs Residual');
legend();
xlabel('Date');

% === High and low tides from zero-crossings of predicted northward current ===
isHigh = islocalmax(v_pred);
isLow  = islocalmin(v_pred);

% Zero crossings of northward predicted tide
zeroCrossings = find(diff(sign(v_pred)) ~= 0);
tideChangeTimes = t_dnum(zeroCrossings);

fprintf('Estimated High/Low Tide Times (zero crossings of v_pred):\n');
for i = 1:length(tideChangeTimes)
    fprintf('%s\n', datestr(tideChangeTimes(i), 'dd-mmm-yyyy HH:MM'));
end

maj2 = [1.509, 4.381, 0.293, 0.704, 0.430, 1.422, 0.315, 0.355]; % Replace with your actual major axis data
min2 = [-0.076, 0.172, -0.088, -0.177, 0.264, 0.074, -0.119, -0.053]; % Replace with your actual minor axis data

% Calculate eccentricity for each pair of major and minor axes
eccentricity = sqrt(1 - (min2 ./ maj2).^2);

% Store the eccentricity values in a new column
% Create a table with the original major and minor axes and the calculated eccentricity
T = table(maj2', min2', eccentricity', 'VariableNames', {'MajorAxis', 'MinorAxis', 'Eccentricity'});

% Display the table with the eccentricity column
disp(T);

amplitude = sqrt((maj2.^2)+(min2.^2));
disp(amplitude)
