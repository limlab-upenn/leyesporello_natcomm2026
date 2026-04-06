% Create figures from data shared by Delware Computational collaboration

clear all

%% Import trajectories data 

UP_data = readmatrix("C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Delaware computational colab\SelectedUp.csv");
DOWN_data = readmatrix("C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Delaware computational colab\SelectedDown.csv");

% column values
% 1 - Enh Position - NaN entry in MATLAB
% 2 - Cell/nucleus index
% 3 - Time (x-value)
% 4 - MS2 Signal (raw)
% 5 - Smoothed MS2 Signal
% 6 - Model predicted signal
% 7 - Promoter state (0-inactive, 1-active)
% 8 - Threshold 

%% Sort data into cell arrays for each nucleus

% enh UP data
unique_vals_UP = unique(UP_data(:,2)); % Find unique values in column 2
UP_split= cell(length(unique_vals_UP),1); % Preallocate cell array
for i = 1:length(unique_vals_UP)
    val = unique_vals_UP(i);
    UP_split{i} = UP_data(UP_data(:,2) == val, :); % Extract rows with this value in column 2
    % Optionally, assign to a variable with dynamic name (not recommended)
    % eval(['A_' num2str(val) ' = split_matrices{i};']);
end

% enh DOWN data
unique_vals_DOWN = unique(DOWN_data(:,2)); % Find unique values in column 2
DOWN_split= cell(length(unique_vals_DOWN),1); % Preallocate cell array
for i = 1:length(unique_vals_DOWN)
    val = unique_vals_DOWN(i);
    DOWN_split{i} = DOWN_data(DOWN_data(:,2) == val, :); % Extract rows with this value in column 2
    % Optionally, assign to a variable with dynamic name (not recommended)
    % eval(['A_' num2str(val) ' = split_matrices{i};']);
end

%% Plot trajectories with predicted state overlaid

% enh UP data
figure(1);
for i = 1:length(UP_split)
    subplot(2,2,i)
    yyaxis left
    plot(UP_split{i}(:,3), UP_split{i}(:,4), 'LineWidth', 0.5, 'Color', 'b'); % plot smoothed MS2 signal
    ylim([0 max(UP_split{i}(:,5))+1])
    ylabel('MS2 Trajectory (A.U.)')
    hold on;
    yyaxis left
    plot(UP_split{i}(:,3), UP_split{i}(:,6), 'LineWidth', 0.5, 'Color', 'b', 'LineStyle','-.'); % plot smoothed MS2 signal
    ylim([0 3])
    hold on;
    yyaxis right
    plot(UP_split{i}(:,3), UP_split{i}(:,7), 'LineWidth', 1, 'Color', 'k'); % plot inferred promoter state
    ylim([0 3])
    ylabel('Inferred Activity State')
    hold off;

    xlim([0 max(UP_split{i}(:,3))])
    xlabel('Time (min)')
    
    ax = gca; % Get current axes handle
    ax.YAxis(1).Color = 'blue'; % Set the color of the right-side y-axis to green
    ax.YAxis(2).Color = 'black'; % Set the color of the right-side y-axis to green
end

% enh DOWN data
figure(2);
for i = 1:length(DOWN_split)
    subplot(2,2,i)
    yyaxis left
    plot(DOWN_split{i}(:,3), DOWN_split{i}(:,4), 'LineWidth', 0.5, 'Color', 'r'); % plot smoothed MS2 signal
    ylim([0 max(DOWN_split{i}(:,5))+1])
    ylabel('MS2 Trajectory (A.U.)')
    hold on;
    yyaxis left
    plot(DOWN_split{i}(:,3), DOWN_split{i}(:,6), 'LineWidth', 0.5, 'Color', 'r', 'LineStyle','-.'); % plot smoothed MS2 signal
    ylim([0 3])
    hold on;
    yyaxis right
    plot(DOWN_split{i}(:,3), DOWN_split{i}(:,7), 'LineWidth', 1, 'Color', 'k'); % plot inferred promoter state
    ylim([0 3])
    ylabel('Inferred Activity State')
    hold off;

    xlim([0 max(DOWN_split{i}(:,3))])
    xlabel('Time (min)')
    
    ax = gca; % Get current axes handle
    ax.YAxis(1).Color = 'red'; % Set the color of the right-side y-axis to green
    ax.YAxis(2).Color = 'black'; % Set the color of the right-side y-axis to green
end

%% select individual nuclei to plot for paper figure
% enhUP - 1, enhDOWN - 4

% enhUP
i = 1;
figure;
yyaxis left
plot(UP_split{i}(:,3), UP_split{i}(:,4), 'LineWidth', 1, 'Color', 'b'); % plot smoothed MS2 signal
ylim([0 3])
ylabel('MS2 Trajectory (A.U.)')
hold on;
yyaxis left
plot(UP_split{i}(:,3), UP_split{i}(:,6), 'LineWidth', 1.5, 'Color', [.5 .5 .5], 'LineStyle','-.'); % plot smoothed MS2 signal
ylim([0 3])
hold on;
yyaxis right
plot(UP_split{i}(:,3), UP_split{i}(:,7), 'LineWidth', 1.5, 'Color', 'k'); % plot inferred promoter state
ylim([0 5])
ylabel('Inferred Activity State')
hold off;

xlim([0 max(UP_split{i}(:,3))])
xlabel('Time (min)')

ax = gca; % Get current axes handle
ax.YAxis(1).Color = 'blue'; % Set the color of the right-side y-axis to green
ax.YAxis(2).Color = 'black'; % Set the color of the right-side y-axis to green

set(gcf,'position',[100,600,600,300])


% enhDOWN
i = 4; 
figure; 
yyaxis left
plot(DOWN_split{i}(:,3), DOWN_split{i}(:,4), 'LineWidth', 1, 'Color', 'r'); % plot smoothed MS2 signal
ylim([0 3])
ylabel('Smoothed MS2 Trajectory (A.U.)')
hold on;
yyaxis left
plot(DOWN_split{i}(:,3), DOWN_split{i}(:,6), 'LineWidth', 1.5, 'Color', [.5 .5 .5], 'LineStyle','-.'); % plot smoothed MS2 signal
ylim([0 3])
hold on;
yyaxis right
plot(DOWN_split{i}(:,3), DOWN_split{i}(:,7), 'LineWidth', 1.5, 'Color', 'k'); % plot inferred promoter state
ylim([0 5])
ylabel('Inferred Activity State')
hold off;

xlim([0 max(DOWN_split{i}(:,3))])
xlabel('Time (min)')

ax = gca; % Get current axes handle
ax.YAxis(1).Color = 'red'; % Set the color of the right-side y-axis to green
ax.YAxis(2).Color = 'black'; % Set the color of the right-side y-axis to green

set(gcf,'position',[100,200,600,300])


%% Import stat analysis of data

stat_data = readmatrix("C:\Users\Emilia Leyes Porello\OneDrive - PennO365\Penn Research\Enhancer-Promoter Spacers\Delaware computational colab\StatsEP.csv");
UP_stats = stat_data(843:end, :);
DOWN_stats = stat_data(1:843, :);

%% Plot stats - boxchart

figure; 
boxchart(ones(1,length(UP_stats)), UP_stats(:,5), 'BoxFaceColor', 'b')
hold on;
boxchart(ones(1,length(DOWN_stats)).*2, DOWN_stats(:,5), 'BoxFaceColor', 'r')
ylabel("# peaks per nucleus")
xticks([1 2])

figure; 
boxchart(ones(1,length(UP_stats)), UP_stats(:,6), 'BoxFaceColor', 'b')
hold on;
boxchart(ones(1,length(DOWN_stats)).*2, DOWN_stats(:,6), 'BoxFaceColor', 'r')
ylabel("Tau (ON)")
xticks([1 2])

figure; 
boxchart(ones(1,length(UP_stats)), UP_stats(:,7), 'BoxFaceColor', 'b')
hold on;
boxchart(ones(1,length(DOWN_stats)).*2, DOWN_stats(:,7), 'BoxFaceColor', 'r')
ylabel("Tau (OFF)")
xticks([1 2])

%% Plot stats - bar graph

figure; 
bar(1, nanmean(UP_stats(:,5)), 'b')
err = nanstd(UP_stats(:,5)) / sqrt(length(UP_stats(:,5)));
hold on; 
er = errorbar(1, nanmean(UP_stats(:,5)), err, err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold on;
bar(2, nanmean(DOWN_stats(:,5)), 'r')
err = nanstd(DOWN_stats(:,5)) / sqrt(length(DOWN_stats(:,5)));
hold on; 
er = errorbar(2, nanmean(DOWN_stats(:,5)), err, err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel("# peaks per nucleus")
xticks([1 2])

set(gcf,'position',[100,200,300,400])

figure; 
bar(1, nanmean(UP_stats(:,6)), 'b')
err = nanstd(UP_stats(:,6)) / sqrt(length(UP_stats(:,6)));
hold on; 
er = errorbar(1, nanmean(UP_stats(:,6)), err, err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold on;
bar(2, nanmean(DOWN_stats(:,6)), 'r')
err = nanstd(DOWN_stats(:,6)) / sqrt(length(DOWN_stats(:,6)));
hold on; 
er = errorbar(2, nanmean(DOWN_stats(:,6)), err, err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel("Tau (ON)")
xticks([1 2])

set(gcf,'position',[500,200,300,400])

figure; 
bar(1, nanmean(UP_stats(:,7)), 'b')
err = nanstd(UP_stats(:,7)) / sqrt(length(UP_stats(:,7)));
hold on; 
er = errorbar(1, nanmean(UP_stats(:,7)), err, err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
hold on;
bar(2, nanmean(DOWN_stats(:,7)), 'r')
err = nanstd(DOWN_stats(:,7)) / sqrt(length(DOWN_stats(:,7)));
hold on; 
er = errorbar(2, nanmean(DOWN_stats(:,7)), err, err);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
ylabel("Tau (OFF)")
xticks([1 2])

set(gcf,'position',[900,200,300,400])
