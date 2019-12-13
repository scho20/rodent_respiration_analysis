%% RODENT BREATHING DATA ANALYSIS
% The objective of the code below is to aid the analysis of rodent
% breathing data. You can look over and manually edit peaks and troughs of
% rodent respiration data and electrophysiological signals of rodent
% medullary neurons, both raw and filtered.

% Initially programmed by Jo Suresh (April 12, 2016) and Tahra Eissa.
% Edited by SungJun Cho for the manual addition of undetected maxima and
% minima of the signal waves. 10/11/2019

clc; clear; close all;
%% Setting Inputs
fprintf('SETTING UP... \n')
tic;

%%%%%%%%%%%%%%%%%%%%%%%%% USER INTERACTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FileName          = input('Enter the filename to analyze (program will add .abf extension as defualt):  ','s');
sr                = input('Enter the sampling frequency in Hz:   '); % unit: Hz
tStartTime        = input('Enter START TIME in secs:   ');
tEndTime          = input('Enter END TIME in secs:   ');
minPeakSeparation = input('Enter MIN PEAK SEPARATION in secs (this will be applied to both peaks and valleys):   '); % usually 0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% May add minPeakProminence based on the need of user

% Integer operands are required for color operator when used as index.
tStartTime = ceil(tStartTime);
tEndTime = ceil(tEndTime);
%% Read the abf file
FileNameOpen = [FileName '.abf'];
[abfData,si,h] = abfload(FileNameOpen);

% Setting time array
x_raw = abfData(:,1)'; % original signal
t_raw = 1/sr:1/sr:length(x_raw)/sr; % original time array
t = tStartTime+(1/sr):1/sr:tEndTime+(1/sr);

filename = [FileName '_' num2str(tStartTime) 's_' num2str(tEndTime) 's_Results.xlsx'];

% Check whether the file already exists
if exist(filename, 'file')
    error('File with same name already exists. \n')
end

elapsed_set_up = toc;
fprintf('TOOK %.3f SECONDS. \n', elapsed_set_up)
%% Plot Original and Segmented Signals
chans = size(abfData,2); % number of channels
x_timeframed = zeros(chans,length(t));
for k = 1:chans
    if tStartTime == 0
        x_timeframed(k,:) = abfData(1:tEndTime*sr,k)';
    else
        x_timeframed(k,:) = abfData(tStartTime*sr:tEndTime*sr,k)';
    end
end

figure(1);
for k = 1:chans
    subplot(chans,1,k);
    plot(t_raw,abfData(:,k));
    title(['Entier Signal of Channel', num2str(k)])
    xlabel('Time (s)');
    ylabel('Amplitude (AU)')
end

figure(2);
for k = 1:chans
    subplot(chans,1,k);
    plot(t,x_timeframed(k,:));
    xlim([tStartTime, tEndTime]);
    xlabel('Time (s)');
    ylabel('Amplitude (AU)');
    title(['Part of signal being analyzed (' num2str(tStartTime) 's - ' num2str(tEndTime) 's)']);
end

%%%%%%%%%%%%%%%%%%%%%%%%% USER INTERACTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
good_chans = input('Which channels would you like to analyze (e.g.[1,2,5])?');
% We want to look at channel 2 and 3: ch2 - integrated breathing; ch3 - raw EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tStartTime == 0
    x_timeframed = abfData(1:tEndTime*sr, good_chans)';
else
    x_timeframed = abfData(tStartTime*sr:tEndTime*sr, good_chans)';
end

fprintf('Please check the figures. Press spacebar to continue. \n')
pause;

%% Signal Filtering
fprintf('Filtering...')
tic;

Ny = sr/2; %Nyquist frequency: half of the sampling rate

%%%%%%%%%%%%%%%%%%%%%%%%% USER INTERACTIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Band-pass filter setting in Hz
fls = input('Enter the bandpass filter setting in Hz (low):    ');
% Band-pass filter setting in Hz
fhs = input('Enter the bandpass filter setting in Hz (high):   ');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fls > fhs
    error('Wrong filter setting: low-pass frequency should be smaller than high-pass.');
end

[b,a] = butter(2,[fls/Ny, fhs/Ny],'bandpass');

% Filtering the corresponding signals
x_filt = zeros(length(good_chans),length(t));
for k = 1:length(good_chans)
    x_filt(k,:) = filtfilt(b,a,x_timeframed(k,:));
end

Breath_chan = 1;
Breath_chan_raw = 2;
EEG_chan = 3;
EEG_chan_raw = 4;

Breath_chan = find(good_chans==Breath_chan);
EEG_chan = find(good_chans==EEG_chan);
Breath_chan_raw = find(good_chans==Breath_chan_raw);
EEG_chan_raw = find(good_chans==EEG_chan_raw);

% Plotting original and filtered signals
figure(4)
if ~isempty(Breath_chan)
    subplot(4,1,1)
    hold on
    plot(t,x_timeframed(Breath_chan,:),'b');
    plot(t,x_filt(Breath_chan,:),'r');
    title('Breathing Data')
    legend('OriginalSignal',['BandPass: ', num2str(fls(1)), 'to', num2str(fhs(1)),'Hz'],'Location','northeastoutside');
    xlabel('Time in secs');
    ylabel('Amplitude AU');
end
if ~isempty(EEG_chan)
    subplot(4,1,2)
    hold on
    plot(t,x_timeframed(EEG_chan,:),'b');
    plot(t,x_filt(EEG_chan,:),'r');
    title('EEG Data')
    legend('OriginalSignal',['BandPass: ', num2str(fls(1)), 'to', num2str(fhs(1)), 'Hz'], 'Location','northeastoutside');
    xlabel('Time in secs');
    ylabel('Amplitude AU');
end
if ~isempty(Breath_chan_raw)
    subplot(4,1,3)
    hold on
    plot(t,x_timeframed(Breath_chan_raw,:),'b');
    plot(t,x_filt(Breath_chan_raw,:),'r');
    title('Rectified Breathing Data')
    legend('OriginalSignal',['BandPass: ', num2str(fls(1)),'to', num2str(fhs(1)),'Hz'], 'Location','northeastoutside');
    xlabel('Time in secs');
    ylabel('Amplitude AU');
end
if ~isempty(EEG_chan_raw)
    subplot(4,1,4)
    hold on
    plot(t,x_timeframed(EEG_chan_raw,:),'b');
    plot(t,x_filt(EEG_chan_raw,:),'r');
    title('Rectified EEG Data')
    legend('OriginalSignal',['BandPass: ', num2str(fls(1)),'to', num2str(fhs(1)),'Hz'],'Location','northeastoutside');
    xlabel('Time in secs');
    ylabel('Amplitude AU');
end

elapsed_filter = toc;
fprintf('took %.3f seconds. \n', elapsed_filter)

fprintf('Please check the figures and press spacebar to continue. \n')
pause;

%% Find Maxima/Peaks
fprintf('Peak detecting ...')
tic;

offset_time = 2; % unit: sec; 
% We need to ignore the first few secs due to the filter artefacts.
x_filt_offset = x_filt(:,offset_time*sr:length(x_filt)); % getting rid of first 2 secs due to initial filter artifact
t_offset = t(offset_time*sr:size(x_filt,2));

locs_max = cell(1, length(good_chans));
locs_min = cell(1, length(good_chans));

for k = 1:length(good_chans)    
    [~,locs_max{k}] = findpeaks(x_filt_offset(k,:),1:1:length(x_filt_offset(k,:)),...
        'Annotate','extents','MinPeakDistance',minPeakSeparation*sr,'MinPeakProminence', 0.001);
    
    x_filt_offset_inverted = -x_filt_offset(k,:);
    
    [~,locs_min{k}] = findpeaks(x_filt_offset_inverted,1:1:length(x_filt_offset(k,:)),...
        'MinPeakDistance',minPeakSeparation*sr,'MinPeakProminence',0.001);
    
    % Minimum distance between peaks is 0.2 secs
end
%% Pull out all the single waves
% We are only going to use all those maxima that are enclosed
% between the first and the last minima.

% Initialize Arrays
cellsize = cell(1,length(good_chans));

final_locs_max = cellsize;
final_locs_min = cellsize;

FINAL_locs_max = cellsize;
wave = cellsize;

wave_rt = cellsize;
wave_ft = cellsize;
FWHM = cellsize;

for k = 1:length(good_chans)
    % (1) Isolating the Single Waves
    
    % [1] Band-Pass Filtered Signal: Identifying Peaks
    % Finding the first and last maxima that occur after the first and the
    % last minima, respectively. Wave constricted by these two index points.
    first_idx_max = find(locs_max{k} > locs_min{k}(1),1);
    last_idx_max = find(locs_max{k} < locs_min{k}(length(locs_min{k})),1,'last');
    final_locs_max_template = locs_max{k}(first_idx_max:last_idx_max);
    
    % (2) Manual Screening
    % Here, the user manually decides to either accept or reject the maxima
    % and minima of the signal wave.
    BP = 'bp'; % bandpass-pass filter
    
    [finalpeaks,finaltrough,pt_max,pt_min] = manualedit_maxmin_bp(x_filt_offset(k,:), final_locs_max_template, locs_min{k},BP);
    
    % [1] Update manually edited maxima and minima (for detected peaks)
    p = 1;
    for i = 1:length(finalpeaks)
        if finalpeaks(i) == 0 % maxima accepted
            final_locs_max{k}(p) = final_locs_max_template(i);
            p = p+1;
        end
    end
    
    q = 1;
    for i = 1:length(finaltrough)
        if finaltrough(i) == 0
            final_locs_min{k}(q) = locs_min{k}(i);
            q = q+1;
        end
    end
    
    % [2] Update additionally marked maxima and minima (for undetected peaks)
    % Cartesian coordinates (x,y) on the plot where x=time(s) and y=amplitude(AU)
    
    % IMPORTANT: When marking additional points on the plot, make sure to
    % remember that the user can only mark the points that are between the
    % first and last minimum or maximum (whichever one is first or last).
    
    % NOTE: You should mark the point in the order of maxima and minima.
    % 1. Click on the maxima you want to add. If you right click your last
    % maxima, it will move on to the next part.
    % 2. Click on the minima you want to add. If you right click your last
    % minima, it will move on to the next part.
    % 3. Finally, exclude any errors by clicking red triangles and blue
    % squares.
    % 4. Click "Save selected peaks" if you are done with manual editing.
    % You will be doing this process twice per one channel.
    % Ex) One for LP and one for HP per raw EEG signal.
    
    pts_min = cat(2,pt_min{1,1},pt_min{1,2});
    pts_max = cat(2,pt_max{1,1},pt_max{1,2});
    
    % Below may need more optimization %%%%%%%%%%%%%%%%%
    pts_min = ceil(pts_min);
    pts_max = ceil(pts_max);
    
    opt_frame = diff(t_offset);
    opt_frame = ceil(10000*opt_frame(1,1)*sr);
    opt_frame_sec = opt_frame/sr; % opt_frame time window in actual seconds
    
    opt_loc_min = zeros(1,size(pts_min,1));
    opt_loc_max = zeros(1,size(pts_max,1));
    
    for i = 1:size(pts_min,1)
        opt_idx_min = min(x_filt_offset(k,pts_min(i,1)-opt_frame:pts_min(i,1)+opt_frame));
        opt_loc_min(1,i) = find(x_filt_offset(k,:) == opt_idx_min);
    end
    
    pts_min(:,1) = opt_loc_min(1,:);
    
    for i = 1:size(pts_max,1)
        opt_idx_max = max(x_filt_offset(k,pts_max(i,1)-opt_frame:pts_max(i,1)+opt_frame));
        opt_loc_max(1,i) = find(x_filt_offset(k,:) == opt_idx_max);
    end
    
    pts_max(:,1) = opt_loc_max(1,:);
    
    % Case: No addition needed
    noaddition = cell(1,2);
    
    noaddition{1,1} = input('Did you additionally marked peaks by clicking the signal plot? [Y/N]: ','s'); % for bandpass maxima
    noaddition{1,2} = input('Did you additionally marked troughs by clicking the signal plot? [Y/N]: ','s'); % for bandpass minima
    
    for i = 1:length(noaddition)
        switch noaddition{1,i}
            case 'Y'
                continue;
            case 'N'
                if i == 1
                    pts_max(:,1) = pts_max(:,1)*0;
                elseif i == 2
                    pts_min(:,1) = pts_min(:,1)*0;
                end
        end
    end
  
    % Include the points designated on the plot by the user
    % Add designated x-ccordinate indices that represent time
    final_locs_max{1,k} = sort([final_locs_max{1,k}, pts_max(:,1)']);
    final_locs_min{1,k} = sort([final_locs_min{1,k}, pts_min(:,1)']);
    
    % Case N: Need to delete zero that has no meaning
    % This method is possible, because we have offsetting 2 seconds.
    final_locs_max{1,k} = final_locs_max{1,k}(final_locs_max{1,k} ~= 0);
    final_locs_min{1,k} = final_locs_min{1,k}(final_locs_min{1,k} ~= 0);
    
    figure;
    hold on;
    plot(t_offset,x_filt_offset(k,:));
    plot((final_locs_max{k}./sr) + offset_time + tStartTime, x_filt_offset(k,final_locs_max{k}),'rv','MarkerFaceColor', 'r');
    plot((final_locs_min{k}./sr) + offset_time + tStartTime, x_filt_offset(k,final_locs_min{k}),'rs','MarkerFaceColor', 'b');
    grid on;
    legend('Filtered Signal','Maxima','Minima','Location','northeastoutside');
    xlabel('Time (s)')
    ylabel('Amplitude (AU)')
    title(['Band-Pass Filtered Signal with 2s offset: Peaks Detected for Channel', num2str(good_chans(k))]);
    
    fprintf('Check the figure and press spaceboar to proceed. \n');
    pause;
    
    % [3] Finalizing the edited waves
    % For band-pass filtered signal
    count = 1;
    for i = 1: length(final_locs_max{k})
        left_idx_min = find(final_locs_min{k} < final_locs_max{k}(i),1,'last');
        right_idx_min = find(final_locs_min{k} > final_locs_max{k}(i),1,'first');
        
        if (x_filt_offset(k,final_locs_min{k}(left_idx_min)) >= x_filt_offset(k,final_locs_max{k}(i))) || ...
            (x_filt_offset(k,final_locs_min{k}(right_idx_min)) >= x_filt_offset(k,final_locs_max{k}(i)))
            continue; % if true, remaining statements (below) in i-th for-loop is essentially skipped
        end
        
        FINAL_locs_max{k}(count) = final_locs_max{k}(i);
        wave{k}{count} = x_filt_offset(k,final_locs_min{k}(left_idx_min):final_locs_min{k}(right_idx_min));
        myWaveTemp = x_filt_offset(k,final_locs_min{k}(left_idx_min):final_locs_min{k}(right_idx_min));
        
        % Risetime and falltime set for the time it takes for the response
        % to rise and fall from 20% to 80% (or vice versa) of the
        % steady-state response, respectively.
        rt_left  = risetime(myWaveTemp,sr,'PercentReferenceLevels',[20 80],'StateLevels',[x_filt_offset(k,final_locs_min{k}(left_idx_min)), x_filt_offset(k,final_locs_max{k}(i))]);
        ft_right = falltime(myWaveTemp,sr,'PercentReferenceLevels',[20 80],'StateLevels',[x_filt_offset(k,final_locs_min{k}(right_idx_min)), x_filt_offset(k,final_locs_max{k}(i))]);
        
        if length(rt_left) ~= 1
            wave_rt{k}(count) = -99;
        else
            wave_rt{k}(count) = rt_left;
        end
        
        if length(ft_right) ~= 1
            wave_ft{k}(count) = - 99;
        else
            wave_ft{k}(count) = ft_right;
        end
        
        % State levels for full width at half maximum (FWHM) are stated as
        % below since there exist cases in which the starting point of
        % upward trace is much lower than the ending point of downward
        % trace.
        
        pw = pulsewidth(myWaveTemp,sr,'StateLevels',[max([x_filt_offset(k,final_locs_min{k}(left_idx_min)) x_filt_offset(k,final_locs_min{k}(right_idx_min))]), x_filt_offset(k,final_locs_max{k}(i))]);
        
        if length(pw) ~= 1
            FWHM{k} = - 99;
        else
            FWHM{k} = pw;
        end
        
        count = count + 1;
    end
    
    elapsed_peak = toc;
    fprintf('took %.3f seconds. \n', elapsed_peak);   
    %% Calculating the period and frequency of the detected peaks
    period = cellsize;
    freq   = cellsize;
    
    t_max_occurance = FINAL_locs_max{k}/sr; % transform unit to seconds
    
    period{1,k} = diff(t_max_occurance);    
    freq{1,k} = 1./period{1,k};
    %% Measuring scores on the irregularity index
    irreg_score = cellsize;
    
    if length(period{k}) < 3
        irreg_score{k} = 0;
        warning('Error: Need at least 3 maxima for band-pass signal for irregularity score. \n')
    else
        for i = 1:length(period{k}) - 1
            irreg_score{k}(i) = abs((period{k}(i+1)-period{k}(i))/period{k}(i));
        end
    end
    %% Writing results in the excel formation
    fprintf('Exporting data to excel...');
    
    header = {'Channel#','Time_of_peaks','Period','Freq','Irreg_Score','Risetime','Falltime','FWHM',...
        'Sample_reate(Hz)','Filter_freq','MinPeakSeparation'};
    
    filename_chan = [FileName '_' num2str(tStartTime) 's_' num2str(tEndTime) 's_Results_chan_' num2str(good_chans(k)) '.xlsx'];
    
    sheet = k;
    xlswrite(filename,header,sheet);
    
    % CAUTION: The excel file is generated only if the innermost folder containing
    % this code is located at Desktop (due to MATLAB default setting). You
    % can manually set up the path according to your interest.
    
    xlRange = 'A2';
    data = {filename_chan};
    xlswrite(filename,data,sheet,xlRange);
    
    xlRange = 'B2';
    t_correction = offset_time + tStartTime;
    xlswrite(filename,((t_max_occurance)+t_correction)',sheet,xlRange);
    
    xlRange = 'C2';
    xlswrite(filename,period{k}',sheet,xlRange);
    
    xlRange = 'D2';
    xlswrite(filename,freq{k}',sheet,xlRange);
    
    xlRange = 'E2';
    xlswrite(filename,irreg_score{k}',sheet,xlRange);
    
    xlRange = 'F2';
    xlswrite(filename,wave_rt{k}',sheet,xlRange);
    
    xlRange = 'G2';
    xlswrite(filename,wave_ft{k}',sheet,xlRange);
    
    xlRange = 'H2';
    xlswrite(filename,FWHM{k}',sheet,xlRange);
    
    xlRange = 'I2';
    data2 = {num2str(sr)};
    xlswrite(filename,data2,sheet,xlRange);
    
    xlRange = 'J2';
    data3 = {num2str(fls),num2str(fhs)};
    xlswrite(filename,data3,sheet,xlRange);
    
    xlRange = 'K2';
    data5 = {num2str(minPeakSeparation)};
    xlswrite(filename,data5,sheet,xlRange);
    
    fprintf('Completed export on excel sheet. \n');
end
