%% SMIM Optimization Code Run File

%SNV to delete befre uploading
close all
clc
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
%% Part 1: User definitions

% filename - Name of the excel sheet containing experimental data
    fileName = 'TemplateExample.xlsx';

% params_upper - Upper bounds to be used in the optimization for the 6 parameters being fit, listed as [V,D,Λ,log_10⁡(T_1),log_10⁡(T_2),β]. 
     params_upper = [1.5 2.3 0.5 2 0.2 12];

% params_lower - Lower bounds to be used in the optimization for the 6 parameters being fit, listed as [V,D,Λ,log_10⁡(T_1),log_10⁡(T_2),β]. 
    params_lower = [0.05 2E-10 0.0001 0  0 0.5];

%% Part 2: Run the SMIM

%Model Outputs:
    %1 - .mat file for each excel sheet containing the complete SMIM output
    %2 - An excel file containing a summary of model outputs (fits and standard errors) for each sheet.


%Extract the site names from excel tabs.
    SiteDescriptions = sheetnames(fileName);

% p - Number of parameters being estimated (The conservative code optimizes 6, but the reactive version can optimize 8 potentially) 
    p = 6;    

%Initialize structures to save data later on.
    SMIMSum = zeros(numel(SiteDescriptions),15); %This will store model outputs and associated standard errors.
    ScaledData = cell(numel(SiteDescriptions), 2); %This will store predicted time and concentration data, with predicted data scaled to match observations.


for i=1:numel(SiteDescriptions)
     
    %Load in excel file data
        SensorLoc = xlsread(fileName, SiteDescriptions(i), 'F4:F4'); 
        Conc = xlsread(fileName, SiteDescriptions(i), 'K:K');
        ReleaseTime = xlsread(fileName, SiteDescriptions(i), 'F3:F3');
        HitTime = xlsread(fileName, SiteDescriptions(i), 'F6:F6');
        RunTime= xlsread(fileName, SiteDescriptions(i), 'B:B');
        CutTime = xlsread(fileName, SiteDescriptions(i), 'F9:F9');
        Day = xlsread(fileName, SiteDescriptions(i), 'C:C');
        DayStart = xlsread(fileName, SiteDescriptions(i), 'F7:F7');
        DayStop = xlsread(fileName, SiteDescriptions(i), 'F10:F10');
        Q = xlsread(fileName, SiteDescriptions(i),'F12:F12');
        TimeStep = xlsread(fileName, SiteDescriptions(i),'F13:F13');
        TracerMass = xlsread(fileName, SiteDescriptions(i), 'F5:F5');

    
    %Clean data

            %Remove negative values
                Conc(find(Conc < 0)) = 0;
    
    
            %Make sure vectors are the same size (empty cells at the end can occassionally
            %carry over from Excel)
                ConcSize = 1:length(Conc);
                RunTime = RunTime(ConcSize);
                Day = Day(ConcSize);
    
    
            % Set the first concentration measurement to the release start.
                IndexRelease = find((RunTime < ReleaseTime) & (Day == DayStart));
                Conc(IndexRelease) = [];
    
                % Make vectors are the same length
                    RunTime(IndexRelease) = [];
                    Day(IndexRelease) = [];
    
          % Set concentrations before the hit time to zero.          
                IndexHit = find((RunTime < HitTime)&(Day == DayStart));
                Conc(IndexHit) = 0;
                
    
          % Rewrite the time vector such that the first time is zero (corrects
          % for sensors reporting analog time)
                t = (0:TimeStep:((TimeStep*length(Conc))-1))';
    
          % Remove data past the cutoff time.
                IndexCut = find(((RunTime > CutTime)& (Day == DayStop)) | (Day > DayStop));
                Conc(IndexCut) = [];
                t(IndexCut) = [];
                Day(IndexCut) = [];


    %Running the SMIM
          % Normalize the data
            if ~isempty(Q)
                [ccNorm, massRecoveryFraction, Qdg] = cNorm(t, Conc, 'c', 'cMass', TracerMass, 'Q', Q); 
            else
                 [ccNorm, massRecoveryFraction, Qdg] = cNorm(t, Conc, 'c', 'cMass', TracerMass); 
            end

          %Initial guess for optimization      
          params_guess  = [(params_upper + params_lower)/2];

         % Run the model and save output      
            ModCon =SMIMfit(t, ccNorm, 'c', 'L', SensorLoc, 't_end', t(length(t)),'model_type', 'TPL', 'params_guess', params_guess, 'params_upper', params_upper, 'params_lower', params_lower)
            save(append(SiteDescriptions(i), '.mat'));
 

   % Calculate Covariance matrix.
     
            %lsqnonlin returns the jacobian as a sparse matrix
               Jacobian = full(ModCon.jacobian); 
                
            %Get the sum of squares weighted residuals
               SSWR = sum((ModCon.cobs - ModCon.ccfit).^2);
    
            %Number of observations
               n = numel(ModCon.cobs);
        
            %Calculated Error Variance 
               CEV = SSWR/(n-p);
    
            %Covariance Matrix
               Covariance = CEV*inv(Jacobian.'*Jacobian);


       
  % Save Model Data
            SMIMSum(i,1) = ModCon.params_fit(1); %U
            SMIMSum(i,2) = ModCon.params_fit(2); % D
            SMIMSum(i,3) = ModCon.params_fit(3); %Lambda
            SMIMSum(i,4) = ModCon.params_fit(4); %Beta
            SMIMSum(i,5) = ModCon.params_fit(5); %logT1
            SMIMSum(i,6) = ModCon.params_fit(6); %logT2
   
            SMIMSum(i, 7) = sqrt(Covariance(1,1)); %Standard Error for U
            SMIMSum(i, 8) = sqrt(Covariance(2,2)); %Standard Error for D
            SMIMSum(i, 9) = sqrt(Covariance(3,3)); %Standard Error for Lambda
            SMIMSum(i, 10) = sqrt(Covariance(4,4)); %Standard Error for Beta
            SMIMSum(i,11) = sqrt(Covariance(5,5)); %Standard Error for logT1
            SMIMSum(i,12) = sqrt(Covariance(6,6)); %Standard Error for logT2
        
            SMIMSum(i,13) = sum(ModCon.resid); %Weighted Mean Absolute Error Objective Output
            SMIMSum(i,14) = massRecoveryFraction; %Mass Recovery Fraction
            SMIMSum(i,15) = Qdg;  %Estimated flow rate
           
           save(append(SiteDescriptions(i), '.mat'));
           

            SMIMTable = table(SiteDescriptions, ...
                SMIMSum(:,1), SMIMSum(:,2), SMIMSum(:,3), SMIMSum(:,4), SMIMSum(:,5), ...
                SMIMSum(:,6), SMIMSum(:,7), SMIMSum(:,8), SMIMSum(:,9), SMIMSum(:,10), SMIMSum(:,11), ...
                SMIMSum(:,12), SMIMSum(:,13), SMIMSum(:,14), SMIMSum(:,15));          
            SMIMTable.Properties.VariableNames(1:16) = {'SiteDescription', 'U', 'D','Lambda', ...
                'Beta', 'logT1', 'logT2', 'SE_U', 'SE_D', 'SE_Lambda', 'SE_B', 'SE_logT1', 'SE_logT2', ...
                'WMAE', 'MassRecoveryFractionQAv', 'QEstimated'}; 
            
            writetable(SMIMTable, append('SMIMResults', '.xls'))    

    
 %  Plot results

        %Remove pre-hit values for plotting

            %Find first non-zero concentration value
                hitIndex = find(Conc, 1, 'first');

            %Remove pre-hit values
                concMass = Conc;
                tMass = t;
                concMass(1:(hitIndex-1)) = [];
                tMass(1:(hitIndex-1)) = [];

        %Identify the peak of the predicted curve
            [predPeak, predPeakIndex] = max(ModCon.ccfit);

        %Identify the peak of the observed data
            [obsPeak, obsPeakIndex] = max(concMass);

        %Scale the predicted peak to match the observed peak
            scalefactor = (obsPeak/predPeak);
            scaledPred = ModCon.ccfit .* scalefactor;

        %Save time and scaled prediction data
            ScaledData{i,1} = ModCon.tcfit;
            ScaledData{i,2} = scaledPred;
                
        %Plot

             NlogFig = figure(1)
            colors = lines(7);
            fontsizeticks = 14; %Size for axis tick labels
          fontsizeaxes = 16;  %Size for axis labels
            
            
            %Plot 1: BTC on a regular scale
                subplot(1,2,1)
                    plot(ModCon.tcfit, scaledPred, 'linewidth', 7, 'color', colors(2, :))
                    hold on
                    plot(tMass,concMass, 'o', 'markersize', 5, 'linewidth', 3,...
                     'color', colors(1, :))      
                    legend('\textbf{SMIM Predicted}', '\textbf{Observed}', 'Interpreter', 'latex',...
                           'FontSize', 12, 'Location', 'northeast') 
                    set(gca, 'FontSize', fontsizeticks, 'LineWidth', 1.5, 'TickDir', 'out', 'FontWeight', 'bold', 'FontAngle', 'italic', 'FontName', 'Times New Roman')
                    box off
                    xlabel('Time (s)', 'FontSize', fontsizeaxes)
                    ylabel('Tracer RWT Concentration [mg/L]', 'FontSize', fontsizeaxes)   

            %Plot 2: BTC on a log-log scale 
                subplot(1,2,2)
                    plot(ModCon.tcfit, scaledPred, 'linewidth', 7, 'color', colors(2, :))
                    hold on
                    plot(tMass,concMass, 'o', 'markersize', 5, 'linewidth', 3,...
                     'color', colors(1, :))    
                    set(gca, 'YScale', 'log', 'FontSize', fontsizeticks, 'LineWidth', 1.5, 'TickDir', 'out', 'FontWeight', 'bold', 'FontAngle', 'italic', 'FontName', 'Times New Roman')
                    box off
                    xlabel('Time (s)', 'FontSize', fontsizeaxes)
                    ylabel('Tracer RWT Concentration [mg/L]', 'FontSize', fontsizeaxes)   

             
            %Save Plots
                   saveas(NlogFig, append(SiteDescriptions(i), '.fig'));
                   close(NlogFig)
                   save(append(SiteDescriptions(i), '.mat'));


end
