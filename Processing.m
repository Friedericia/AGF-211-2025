%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data folder and subfolders must be added to path before running
% A Plot folder should already exist to save png files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% Read data

Minilog05m = readtable('Minilog-II-T_355571_20250405_1.csv','NumHeaderLines',8);
Minilog15m = readtable('Minilog-II-T_355583_20250405_1.csv','NumHeaderLines',8);
Minilog25m = readtable('Minilog-II-T_355572_20250405_1.csv','NumHeaderLines',8);
Minilog35m = readtable('Minilog-II-T_355573_20250405_1.csv','NumHeaderLines',8);
Minilog45m = readtable('Minilog-II-T_358946_20250405_1.csv','NumHeaderLines',8);

Seabird01m = readtable('SBE37SMP-RS232_03720169_2025_04_05.csv');
Seabird20m = readtable('SBE37SMP-RS232_03727333_2025_04_05.csv');
Seabird30m = readtable('SBE37SM-RS232_03710965_2025_04_05.csv');

Seaguard10m = readtable('1711.txt',NumHeaderLines=5);
Seaguard40m = readtable('1705.txt',NumHeaderLines=5);

Seastar01m = readtable('1S13682DAT.CSV');
Seastar20m = readtable('1S13687DAT.CSV');
Seastar30m = readtable('1S13686DAT.CSV');

ADCP = readtable('adcp.txt',NumHeaderLines=12);

%% Extract relevant variables and crop to effective timeframe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mooring from 02/04/2025 13:00 UTC to 05/04/2025 10:40 UTC
% ADCP from 02/04/2025 15:40 UTC to 05/04/2025 10:40 UTC
% T=Temperature (°C), S=Salinity (PSU), C=Conductivity (mS/cm)
% U=Eastward current component (cm/s), V=Northward, W=Vertical
% AbsU=Magnitude of current, Dir=Direction of current (0°=North, 90°=East)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_01m = Seabird01m.Temperature_ITS_90DegC_(13:222);
S_01m = Seabird01m.PracticalSalinity_PSU_(13:222);

T_05m = Minilog05m.Var3(13:222);

T_10m = Seaguard10m.Temperature_Deg_C__1(12:221);
C_10m = Seaguard10m.Conductivity_mS_cm_(12:221);
V_10m = Seaguard10m.North_cm_s_(12:221);
U_10m = Seaguard10m.East_cm_s_(12:221);
AbsU_10m = Seaguard10m.AbsSpeed_cm_s_(12:221);
Dir_10m = Seaguard10m.Direction_Deg_M_(12:221);
P_10m = Seaguard10m.Pressure_kPa_(12:221);

T_15m = Minilog15m.Var3(13:222);

T_20m = Seabird20m.Temperature_ITS_90DegC_(13:222);
S_20m = Seabird20m.PracticalSalinity_PSU_(13:222);

T_25m = Minilog25m.Var3(13:222);

T_30m = Seabird30m.Temperature_ITS_90DegC_(13:222);
S_30m = Seabird30m.PracticalSalinity_PSU_(13:222);

T_35m = Minilog35m.Var3(13:222);

T_40m = Seaguard40m.Temperature_Deg_C_(12:221);
C_40m = Seaguard40m.Conductivity_mS_cm_(12:221);
V_40m = Seaguard40m.North_cm_s_(12:221);
U_40m = Seaguard40m.East_cm_s_(12:221);
AbsU_40m = Seaguard40m.AbsSpeed_cm_s_(12:221);
Dir_40m = Seaguard40m.Direction_Deg_M_(12:221);

T_45m = Minilog45m.Var3(13:222);

U_ADCP = NaN([210,25]);
U_ADCP(9:210,:) = table2array(ADCP(5:206,10:34)); U_ADCP = U_ADCP'./10;
V_ADCP = NaN([210,25]);
V_ADCP(9:210,:) = table2array(ADCP(5:206,35:59)); V_ADCP = V_ADCP'./10;
W_ADCP = NaN([210,25]);
W_ADCP(9:210,:) = table2array(ADCP(5:206,60:84)); W_ADCP = W_ADCP'./10;

%% Creation of arrays

time = (0:209)./3; % hours since 02/04/2025 13:00 UTC
t0 = datetime(2025,4,2,13,0,0); % 02 April 2025 13:00 UTC
date_time = hours(time)+t0; % time and date in datetime format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose a case for which combination of instruments to use
% Case 1 : All / 2 : Minis + Sbe37 / 3 : Seaguards + Sbe37 
% 4 : Minis / 5 : Sbe37 / 6 : Seaguards + Sbe37 + 45m mini
% For salinity, only 1=3 and 5 are valid (no salinity sensor on minis)
case_T = 6;
case_S = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "depth_var" is a list of the depths at which the chosen instruments appear
switch case_T
    case 1
        depths_T = [1,5,10,15,20,25,30,35,40,45];
        T = [T_01m,T_05m,T_10m,T_15m,T_20m,T_25m,T_30m,T_35m,T_40m,T_45m]';
        label_T = 'All Instruments';
    case 2
        depths_T = [1,5,15,20,25,30,35,45];
        T = [T_01m,T_05m,T_15m,T_20m,T_25m,T_30m,T_35m,T_45m]';
        label_T = 'Minilogs + Sbe37';
    case 3
        depths_T = [1,10,20,30,40];
        T = [T_01m,T_10m,T_20m,T_30m,T_40m]';
        label_T = 'Seaguards + Sbe37';
    case 4
        depths_T = [5,15,25,35,45];
        T = [T_05m,T_15m,T_25m,T_35m,T_45m]';
        label_T = 'Minilogs Only';
    case 5
        depths_T=[1,20,30];
        T = [T_01m,T_20m,T_30m]';
        label_T = 'Sbe37 Only';
    case 6
        depths_T = [1,10,20,30,40,45];
        T = [T_01m,T_10m,T_20m,T_30m,T_40m,T_45m]';
        label_T = 'Seaguards + Sbe37 + 45m Mini';
end

switch case_S
    case {1,3}
        depths_S = [1,10,20,30,40];
        S = [S_01m,S_10m,S_20m,S_30m,S_40m]';
        label_S = 'Seaguards + Sbe37';
    case 5
        depths_S = [1,20,30];
        S = [S_01m,S_20m,S_30m]';
        label_S = 'Sbe37 Only';
end

depths_ADCP = 3:2:51;

%% Plot Temperature and Salinity on same figure

figure();
tiledlayout(2,1);

% Temperature plot
nexttile;
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p1] = contourf(time,depths_T,T,17);
set(p1, 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
colormap(gca,turbo(17));
clim([-1.87,-1.7]);
cb1 = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
ylabel('Depth (m)');
title(['Water Temperature - ',label_T]);
ylabel(cb1,'Temperature (°C)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;
hold off;

% Salinity plot
nexttile;
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p2] = contourf(time,depths_S,S,10);
set(p2,'EdgeColor','none');
set(gca,'Ydir','reverse');
colormap(gca,flipud(bone(10)));
cb2 = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
box off;
xlabel('Time (h)');
ylabel('Depth (m)');
title(['Salinity (PSU) - ',label_S])
ylabel(cb2,'Salinity (PSU)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;
hold off;

saveas(gcf,strcat('Plots/1_Temp_Sal_',string(case_T),'_',string(case_S),'.png'));

%% Plot Temperature and Salinity on different figures

figure();
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p1] = contourf(time,depths_T,T,17);
set(p1, 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
colormap(gca,turbo(17));
clim([-1.87,-1.7]);
cb1 = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
ylabel('Depth (m)');
title(['Water Temperature - ',label_T]);
ylabel(cb1,'Temperature (°C)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;
hold off;

saveas(gcf,strcat('Plots/Temp_',string(case_T),'.png'));

figure;
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p2] = contourf(time,depths_S,S,10);
set(p2,'EdgeColor','none');
set(gca,'Ydir','reverse');
colormap(gca,flipud(bone(10)));
cb2 = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
box off;
xlabel('Time (h)');
ylabel('Depth (m)');
title(['Salinity (PSU) - ',label_S])
ylabel(cb2,'Salinity (PSU)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;
hold off;

saveas(gcf,strcat('Plots/Sal_',string(case_S),'.png'));

%% Plot currents

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs stick plot with color
% x=time, y=0 (horizontal line)
function quiverc(x, y, u, v, mag)

    % Set up colormap
    cmap = jet(16);
    cmin = 0;
    cmax = 15;
    cidx = round( (mag - cmin) / (cmax - cmin) * 15 ) + 1;
    cidx = max(min(cidx, 16), 1);

    hold on;
    for i = 1:length(x)
        x0 = x(i);
        y0 = y(i);
        x1 = x0 + u(i);
        y1 = y0 + v(i);

        % Color
        color = cmap(cidx(i), :);

        % Plot line (vector)
        plot([x0 x1], [y0 y1], '-', 'Color', color, 'LineWidth', 2);
    end

    % Add colorbar and format
    colormap(cmap);
    clim([cmin cmax]);
    cb = colorbar;
    ylabel(cb,'Current Speed (cm/s)')
    axis equal;

    hold off;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*

figure();
tiledlayout(2,1);

nexttile;
quiverc(time,zeros(size(time)),U_10m',V_10m',AbsU_10m');
title('Current Velocity')
xlim([-2,72]);
ylim([-15,6]);
set(gca, 'YTick', [])
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
ylabel('10m Depth');
grid on;
annotation('textarrow', ...
    [0.06 0.06], [0.8 0.85], ...  % x and y in normalized figure coords
    'String', 'N', ...
    'FontWeight', 'bold', ...
    'FontSize', 12);

nexttile;
quiverc(time,zeros(size(time)),U_40m',V_40m',AbsU_40m);
xlim([-2,72]);
ylim([-15,6]);
xt = 0:6:70;
set(gca, 'YTick', [])
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
ylabel('40m Depth');
grid on;

saveas(gcf,'Plots/2_Currents_Seaguards.png')

%% Plot ADCP data

figure();
tiledlayout(2,1);

nexttile;
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p3] = contourf(time,depths_ADCP,V_ADCP,21);
title('ADCP Currents');
set(p3, 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
colormap(gca,turbo(21));
clim([-15,6]);
cb = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
ylabel('Depth (m)');
ylabel(cb,'Northward Current Speed (cm/s)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;
hold off;

nexttile;
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p4] = contourf(time,depths_ADCP,U_ADCP,21);
set(p4, 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
colormap(gca,turbo(21));
clim([-15,6]);
cb = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
ylabel('Depth (m)');
ylabel(cb,'Eastward Current Speed (cm/s)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;
hold off;


saveas(gcf,'Plots/3_Currents_ADCP.png')

%% Comparison of temperature and tidal motion

figure();

% yyaxis left;
hold on;
fill([0,0,70,70],[0,51,51,0],[0.85,0.85,0.85],EdgeColor="none");
[~,p1] = contourf(time,depths_T,T,17);
set(p1, 'EdgeColor', 'none');
set(gca, 'YDir', 'reverse');
colormap(gca,turbo(17));
clim([-1.87,-1.7]);
cb1 = colorbar;
ylim([-5,54]);
xlim([0,69.66]);
yline(0,'-','Sea ice',Color='b',LineWidth=1);
yline(51,'-','Floor',Color='r',LineWidth=1,LabelVerticalAlignment='middle');
ylabel('Depth (m)');
title('Water Temperature, NS Currents and Tides');
ylabel(cb1,'Temperature (°C)')
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');
box off;

%yyaxis right;
scale = 1;
plot(time,-V_10m*scale+10, '-k', 'LineWidth', 1);
yline(10,'--','10m',Color='k',LineWidth=1);
plot(time,-V_40m*scale+40, '-k', 'LineWidth', 1);
yline(40,'--','40m',Color='k',LineWidth=1);

% High and low tides
x10H=[1.33,14,26.33,38.67,51.33,63.67];
x40H=[1.33,15,26,39,51,64];
x10L=[7.33,20.33,32,45,57];
x40L=[7.33,19.67,32.33,44.66,57];

plot(x10H,[10,10,10,10,10,10],'ro','MarkerFaceColor', 'r', MarkerSize=8);
plot(x40H,[40,40,40,40,40,40],'ro','MarkerFaceColor', 'r', MarkerSize=8);
plot(x10L,[10,10,10,10,10],'bo','MarkerFaceColor', 'b', MarkerSize=8);
plot(x40L,[40,40,40,40,40],'bo','MarkerFaceColor', 'b', MarkerSize=8);

saveas(gcf,strcat('Plots/4_Temp_vs_Tides.png'));

%% Demean and normalizing signals for correlation analysis

T_10m_proc = (T_10m-mean(T_10m)) / std(T_10m-mean(T_10m));
T_40m_proc = (T_40m-mean(T_40m)) / std(T_40m-mean(T_40m));
V_10m_proc = (V_10m-mean(V_10m)) / std(V_10m-mean(V_10m));
V_40m_proc = (V_40m-mean(V_40m)) / std(V_40m-mean(V_40m));

%% Correlation beetween temperatures and tides

figure();
set(gcf, 'Position', [100, 242, 1120, 420]); % [left, bottom, width, height]
tiledlayout(2,2);

% T vs V 10m
nexttile;
yyaxis right;
plot(time,V_10m_proc);
ylabel('NS current');
yyaxis left;
plot(time,T_10m_proc);
ylabel('Temperature');
title('Temperature and NS Current at 10m (De-meaned and Normalized)');
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');

% T vs V 40m
nexttile;
yyaxis right;
plot(time,V_40m_proc);
ylabel('NS current');
yyaxis left;
plot(time,T_40m_proc);
ylabel('Temperature');
title('Temperature and NS Current at 40m (De-meaned and Normalized)');
xt = 0:6:70;
xticks(xt);
xt_dt = t0 + hours(xt);
xticklabels(datestr(xt_dt, 'dd, HH:MM'));
xlabel('Day and Time in April 2025 (UTC)');

% xcorr 10m
nexttile;
[xcorr_val_T_10m, lags_T_10m] = xcorr(T_10m_proc, V_10m_proc,'coeff');
plot(-lags_T_10m/3, xcorr_val_T_10m,Color='k');
xlabel('Lag (h)');
ylabel('Cross-correlation');
title('x-corr between Temperature and NS Current at 10m');
grid on;
xlim([0,50]);
ylim([-0.4,0.4]);
%xline(-1,'--','-1h00',Color='r',LineWidth=1,LabelVerticalAlignment='bottom');

% xcorr 40m
nexttile;
[xcorr_val_T_40m, lags_T_40m] = xcorr(T_40m_proc, V_40m_proc,'coeff');
plot(-lags_T_40m/3, xcorr_val_T_40m,Color='k');
xlabel('Lag (h)');
ylabel('Cross-correlation');
title('x-corr between Temperature and NS Current at 40m');
grid on;
xlim([0,50]);
ylim([-0.4,0.4]);
%xline(-1.33,'--','-1h20',Color='r',LineWidth=1,LabelVerticalAlignment='bottom');

saveas(gcf,strcat('plots/5_xcorr_T_V.png'))

%% FFT temperature

dt = 1200;              % Time interval between samples (s)
Fs = 1 / dt;            % Sampling period (Hz)
N = length(T_10m);      % Sample number

T_10m_fft = fft(T_10m_proc);
T_40m_fft = fft(T_40m_proc);

f = (0:N-1)*(Fs/N);     % Frequencies
T_10m_fft_mag = abs(T_10m_fft) / N;  % Amplitude
T_40m_fft_mag = abs(T_40m_fft) / N;

figure;
tiledlayout(2,1);

nexttile;
plot(f, T_10m_fft_mag,Color='k');
title('Fast Fourier Transform of Temperature at 10m Depth');
ylabel('Amplitude');
xlim([0, Fs/4]);
grid on;
xlim_vals = get(gca, 'XLim');
set(gca, 'XLim', xlim_vals)
xticks_hz = get(gca, 'XTick');
tick_step = mean(diff(xticks_hz));  % espacement moyen
xticks_hz_dense = xlim_vals(1):tick_step/4:xlim_vals(2);
xticks_hours = 1 ./ xticks_hz_dense / 3600;
xticks_hours(isinf(xticks_hours)) = NaN;
set(gca, 'XTick', xticks_hz_dense, 'XTickLabel', round(xticks_hours, 2))
xlabel(gca, 'Period (h)', 'Color', 'k')

f_M2 = 1/(12.4206*3600);
f_day = 1/(24*3600);
f_6h = 1/(6*3600);
xline(f_M2,'--','M2 period',Color='r')
xline(f_day,'--','24h period',Color='r')
xline(f_6h,'--','6h period',Color='r')

nexttile;
plot(f, T_40m_fft_mag,Color='k');
title('Fast Fourier Transform of Temperature at 40m Depth');
ylabel('Amplitude');
xlim([0, Fs/4]);
grid on;
xlim_vals = get(gca, 'XLim');
set(gca, 'XLim', xlim_vals)
xticks_hz = get(gca, 'XTick');
tick_step = mean(diff(xticks_hz));  % espacement moyen
xticks_hz_dense = xlim_vals(1):tick_step/4:xlim_vals(2);
xticks_hours = 1 ./ xticks_hz_dense / 3600;
xticks_hours(isinf(xticks_hours)) = NaN;
set(gca, 'XTick', xticks_hz_dense, 'XTickLabel', round(xticks_hours, 2))
xlabel(gca, 'Period (h)', 'Color', 'k')

f_M2 = 1/(12.4206*3600);
f_day = 1/(24*3600);
f_6h = 1/(6*3600);
xline(f_M2,'--','M2 period',Color='r')
xline(f_day,'--','24h period',Color='r')
xline(f_6h,'--','6h period',Color='r')

saveas(gcf,strcat('PLots/7_FFT_T.png'));
