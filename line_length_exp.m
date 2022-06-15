% REANALYSIS OF SIMS, JACOBS, AND KNILL (2012) - LINE LENGTH EXPERIMENT
% Code by Dagmar Adamcova, 2022

% SIGMA - low / high variance condition
% N - set size
% DELTA - true perturbation
% DELTA_LOG - true pertubation as the difference between log transformed
% initial and perturbed stimulus length
% RESP_BIGGER - participant response

%% LOAD DATA
%data = readtable("merged_data_line_length.csv");
data = readtable("merged_data_line_length_log_bins.csv");  % this file has bins based on the log of the data, (gaussian distribution)

%% FIT PROBIT MODEL
% to estimate the psychometric function

varianceConds = [0.0748, 0.299];  % 0.0748 or 0.299 (sd in the variance conditions)
setSizes = [1,2,4,8];             % 1,2,4,8
binNums = [1:10];


% VECTORS
sig = [];
N = [];
bin = [];
bin_mean = [];
param = [];     % vector for bias (mean)
param2 = [];    % vector for standard deviation

for varianceCond = varianceConds
    for setSize = setSizes
        for binNum = binNums
        
            %x = unique(data.DELTA)';  % delta values (in cm)
            %x = unique(data.DELTA_LOG)';  % delta values (log)
            x = unique(data((data.SIGMA == varianceCond & data.N == setSize & data.BIN == binNum), :).DELTA_LOG)';
            n = [];                       % vector for total n trials corresponding to each delta 
            y = [];                       % vector for total n of positive responses corresponding to each delta
            
            for truePert = x 
            
                % Count total number of trials
                n = [n; height(data((data.SIGMA == varianceCond & data.N == setSize & data.BIN == binNum & data.DELTA_LOG == truePert), :))'];
            
                % Count RESP_BIGGER = 1
                y = [y; height(data((data.SIGMA == varianceCond & data.N == setSize & data.BIN == binNum & data.DELTA_LOG == truePert & data.RESP_BIGGER == 1), :))];
            
            end
            
            % FIT MODEL
            
            [b,dev,stats] = glmfit(x,[y n],'binomial','Link','probit');     % fit model            
            
            yfit = glmval(b,x,'probit','Size',n);                           % get model predictions
            
            %plot(x,y./n,'o',x,yfit./n,'-')                   % plot observed and predicted y
            
            % SAVE PARAMETERS
        
            sig = [sig; varianceCond];
            N = [N; setSize];
            bin = [bin; binNum];
            bin_mean = [bin_mean; unique(data((data.SIGMA == varianceCond & data.BIN == binNum), :).BIN_MEAN)];
            param = [param; -b(1)/b(2)]; % bias (mean)
            param2 = [param2; 1/b(2)];   % standard deviation
        
        end
    end
end

% SAVE TABLE
paramTable = array2table([sig, N, bin, bin_mean, param, param2], 'VariableNames', {'SIGMA', 'N', 'BIN', 'BIN_MEAN', 'MEAN', 'STD'});





%% GET MEAN OF SD & SLOPE ESTIMATES 
% for each condition (low, high), for each set size (1, 2, 4, 8)

varLow = 0.0748; 
varHigh =  0.299;
setSizes = [1,2,4,8]; 

lowN = [];      % vector for set size
lowSD = [];     % vector for computed mean SD
lowSlopes = []; % vector for slope estimates
highN = [];
highSD = [];
highSlopes = [];

% Low variance
for setSize = setSizes

    % GET MEAN OF SD
    meanSD = mean(paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).STD);
    lowSD = [lowSD; meanSD];
    lowN = [lowN; setSize];

    % GET SLOPE ESTIMATES USING POLYFIT
    x = paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).BIN_MEAN;
    y = paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).MEAN;

    p = polyfit(x,y,1);  % fit linear regression model

    lowSlopes = [lowSlopes; p(1)]; % save slope estimate

end 

% High variance
for setSize = setSizes

    % GET MEAN OF SD
    meanSD = mean(paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).STD);
    highSD = [highSD; meanSD];
    highN = [highN; setSize];

    % GET SLOPE ESTIMATES USING POLYFIT
    x = paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).BIN_MEAN;
    y = paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).MEAN;

    p = polyfit(x,y,1);  % fit linear regression model

    highSlopes = [highSlopes; p(1)]; % save slope estimate

end 

meanSDTableLow = array2table([lowN, lowSD, lowSlopes], 'VariableNames', {'SET_SIZE', 'MEAN_SD', 'SLOPE'});
meanSDTableHigh = array2table([highN, highSD, highSlopes], 'VariableNames', {'SET_SIZE', 'MEAN_SD', 'SLOPE'});

mergedStruct.Low_Variance = meanSDTableLow;
mergedStruct.High_Variance = meanSDTableHigh;

struct2table(mergedStruct)

%save('slopes_sd_table_line_lengths.mat', 'mergedStruct')


%% BAYESIAN MODEL PREDICTIONS
% sigma_s = sigma_x / (1 + bias_slope)
% for each variance condition, for each set size

% sigma_x = mergedStruct.Low_Variance.MEAN_SD
% bias_slope = mergedStruct.Low_Variance.SLOPE


% Calculate sigma_s
sigma_s_low = mergedStruct.Low_Variance.MEAN_SD / (1 + mergedStruct.Low_Variance.SLOPE);
sigma_s_high = mergedStruct.High_Variance.MEAN_SD / (1 + mergedStruct.High_Variance.SLOPE);


% Plot
plot(mergedStruct.Low_Variance.SET_SIZE, sigma_s_low(:,1), 'marker','o')
hold on
plot(mergedStruct.Low_Variance.SET_SIZE, sigma_s_high(:,1), 'marker','o')

legend(["Low variance", "High variance"])
xlabel('Set Size')
ylabel('Sigma_s')
sgtitle('Line length experiment')



%% INFORMATION-LIMIT MODEL (FIND STARTING PARAMETERS for R and SIGMA_S)
% σ_e = sqrt[ (σ_s)^2 + ( (σ_w)^2 - (σ_s)^2 ) / ( exp(2 R/N) - 1 ) ]
% where σ_s and R are unknown (free parameters) 
% σ_w is the SD of the population the samples are drawn from (differing between high/low variance conditions)

% free parameters (same for both high and low variance conditions)
%sig_s = 0.17;
%R = 3;

% use these after fminsearch
sig_s = x(1);
R = x(2);

% variables
varLow = 0.0748;            % σ_w
varHigh =  0.299;           % σ_w
setSizes = [1,2,4,8];       % N

% vectors
sigma_e_low = [];
sigma_e_high = [];

for setSize = setSizes

    sigma_e_low = [sigma_e_low; sqrt( sig_s^2 + ((varLow^2 - sig_s^2) / (exp(2*R / setSize) - 1)) )];
    sigma_e_high = [sigma_e_high; sqrt( sig_s^2 + ((varHigh^2 - sig_s^2) / (exp(2*R / setSize) - 1)) )];

end

% Plot
plot(mergedStruct.Low_Variance.SET_SIZE, sigma_s_low(:,1), 'marker','o')
hold on
plot(mergedStruct.Low_Variance.SET_SIZE, sigma_s_high(:,1), 'marker','o')

plot(mergedStruct.Low_Variance.SET_SIZE, sigma_e_low(:,1), 'marker','*')
plot(mergedStruct.Low_Variance.SET_SIZE, sigma_e_high(:,1), 'marker','*')

legend(["Low - Sigma_s", "High - Sigma_s", "Low - Sigma_e", "High - Sigma_e",])
%legend()
xlabel('Set Size')
ylabel('Sigma_e')
sgtitle('Line length experiment')


%% INFORMATION-LIMIT MODEL - FMINSEARCH
% minimise the summed squared error
% fminsearch(@(B) sum( (D - sqrt[ (σ_s)^2 + ( (σ_w)^2 - (σ_s)^2 ) / ( exp(2 R/N) - 1 ) ]).^2 ), ...
% use column vectors

options = optimset('PlotFcns',@optimplotfval); % monitor the minimization process

%f = @(x) sum(D - sqrt( x(1)^2 + ((varianceConds.^2 - x(1)^2) ./ (exp(2*x(2) ./ setSizes) - 1)) ).^2) ;
%f = @(x, varianceConds, setSizes)sqrt( x(1)^2 + ((varianceConds.^2 - x(1)^2) ./ (exp(2*x(2) ./ setSizes) - 1)) );
%varianceConds = [0.0748];  % 0.0748 or 0.299 (sd in the variance conditions)
%setSizes = [2];             % 1,2,4,8

D = [sigma_s_low(:,1)', sigma_s_high(:,1)']';  % data points from the Bayesian predictions plot
varianceConds = [0.0748, 0.0748, 0.0748, 0.0748, 0.299, 0.299, 0.299, 0.299]';  % 0.0748 or 0.299 (sd in the variance conditions)
setSizes = [1,2,4,8,1,2,4,8]';             % 1,2,4,8

fun = @(x) sum(D - sqrt( x(1)^2 + ((varianceConds.^2 - x(1)^2) ./ (exp(2*x(2) ./ setSizes) - 1)) ).^2) ;
x0 = [0.17, 3];          % starting parameters: sig_s = 0.17; R = 3;
x = fminsearch(fun, x0, options)




%% 1. PLOT FOR EACH SET SIZE - LOW VAR MEAN

varLow = 0.0748; 
varHigh =  0.299;

tiledlayout(1,4)

for setSize = setSizes

    % GET SLOPE ESTIMATES USING POLYFIT
    x = paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).BIN_MEAN;
    y = paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).MEAN;

    p = polyfit(x,y,1);  % fit linear regression model
    f = polyval(p,x);  % get predictions for y
    %plot(x,y,'o',x,f,'-') 


    nexttile
    plot(paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).BIN_MEAN, paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).MEAN, 'o')
    hold on
    plot(x,f,'-')


    %ylim([-0.6 0.4]);
    ylim([-0.20 0.15]);
    %legend(num2str(setSize))
    xlabel('Mean line length (log) in decile')
    ylabel('Bias')
    titleString = ['Low variance - Set size ', num2str(setSize)];
    title(titleString)
    sgtitle('Line length experiment')

    hold off
end

%hold off


%% 2. PLOT SEPERATE FOR EACH SET SIZE - HIGH VAR MEAN

varLow = 0.0748; 
varHigh =  0.299;

tiledlayout(1,4)

for setSize = setSizes

    % GET SLOPE ESTIMATES USING POLYFIT
    x = paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).BIN_MEAN;
    y = paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).MEAN;

    p = polyfit(x,y,1);  % fit linear regression model
    f = polyval(p,x);  % get predictions for y
    %plot(x,y,'o',x,f,'-') 


    nexttile
    plot(paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).BIN_MEAN, paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).MEAN, 'o')
    hold on
    plot(x,f,'-')
   

    %ylim([-1.7 0.8]);
    ylim([-0.5 0.5]);
    %legend(num2str(setSize))
    xlabel('Mean line length (log) in decile')
    ylabel('Bias')
    titleString = ['High variance - Set size ', num2str(setSize)];
    title(titleString)
    sgtitle('Line length experiment')


    hold off
end

%hold off


%% 3. PLOT SEPERATE FOR EACH SET SIZE - LOW VAR SD

varLow = 0.0748; 
varHigh =  0.299;    

tiledlayout(1,4)

for setSize = setSizes

    temp_table = mergedStruct.Low_Variance;
    temp_sd = zeros(height(paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).BIN_MEAN),1);
    temp_sd(:) = temp_table((temp_table.SET_SIZE == setSize),:).MEAN_SD;
    %temp_table((temp_table.SET_SIZE == setSize),:).MEAN_SD

    nexttile
    plot(paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).BIN_MEAN, paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).STD, 'o')  % 'marker'
    hold on
    plot(paramTable((paramTable.SIGMA == varLow & paramTable.N == setSize),:).BIN_MEAN, temp_sd,'-')
   

    ylim([0.1 0.23]);
    legend(["sd", "mean of sd"])
    xlabel('Mean line length (log) in decile')
    ylabel('Standard deviation')
    titleString = ['Low variance - Set size ', num2str(setSize)];
    title(titleString)
    sgtitle('Line length experiment')

    hold off
end

%hold off

%% 4. PLOT SEPERATE FOR EACH SET SIZE - HIGH VAR SD

varLow = 0.0748; 
varHigh =  0.299;

tiledlayout(1,4)

for setSize = setSizes

    temp_table = mergedStruct.High_Variance;
    temp_sd = zeros(height(paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).BIN_MEAN),1);
    temp_sd(:) = temp_table((temp_table.SET_SIZE == setSize),:).MEAN_SD;

    nexttile
    plot(paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).BIN_MEAN, paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).STD, 'o') % 'marker',
    hold on
    plot(paramTable((paramTable.SIGMA == varHigh & paramTable.N == setSize),:).BIN_MEAN, temp_sd,'-')
    

    ylim([0.1 0.6]);
    legend(["sd", "mean of sd"])
    xlabel('Mean line length (log) in decile')
    ylabel('Standard deviation')
    titleString = ['High variance - Set size ', num2str(setSize)];
    title(titleString)
    sgtitle('Line length experiment')

    hold off
end

%hold off

