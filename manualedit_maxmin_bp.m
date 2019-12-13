function [finalpeaks,finaltrough,pt_max,pt_min] ...
    = manualedit_maxmin_bp(wave,idx_peak,idx_trough,srng)
close all

idx_peak_flags = idx_peak*0;
idx_trough_flags = idx_trough*0;

fig = figure;
hold all;
grid on;
plot(wave)
handles_peak = gobjects(1,length(idx_peak));

pt_max = cell(1,2); pt_min = cell(1,2);

for i = 1:length(idx_peak)
    handles_peak(i) = plot(idx_peak(i),wave(idx_peak(i)), 'v', ...
        'MarkerFaceColor', 'r', 'MarkerSize', 12, 'MarkerEdgeColor',...
        'k','ButtonDownFcn',@buttondownfcn_clickmarker_peak);
end

handles_trough = gobjects(1,length(idx_trough));

for j = 1:length(idx_trough)
    handles_trough(j) = plot(idx_trough(j),wave(idx_trough(j)), 's', ...
        'MarkerFaceColor','b','MarkerSize',12,'MarkerEdgeColor', ...
        'r', 'ButtonDownFcn', @buttondownfcn_clickmarker_trough);
end

finish = uicontrol('Parent',fig,'Style','pushbutton','Units','Normalized', ...
    'Position',[0.05 0.01 0.5 0.05], 'String','Save selected peaks',...
    'Callback',{@pushbutton_callback});

if srng == 'bp'
    [x_max,y_max] = getpts;
    pt_max{1,1} = x_max;
    pt_max{1,2} = y_max;
end

if srng == 'bp'
    [x_min,y_min] = getpts;
    pt_min{1,1} = x_min;
    pt_min{1,2} = y_min;
end

fprintf('Press spacebar to continue. \n')
pause;

    function buttondownfcn_clickmarker_peak(source,~)
        h_peak = find(handles_peak == source,1);
        if ~isempty(h_peak)
            if idx_peak_flags(h_peak) == 0
                set(handles_peak(h_peak),'MarkerFaceColor','k')
                idx_peak_flags(h_peak) = 1;
            elseif idx_peak_flags(h_peak) == 1
                set(handles_peak(h_peak), 'MarkerFaceColor','r')
                idx_peak_flags(h_peak) = 0;
            end
        end
    end
    
    function buttondownfcn_clickmarker_trough(source,~)
        h_trough = find(handles_trough == source,1);
        if ~isempty(h_trough)
            if idx_trough_flags(h_trough) == 0
                set(handles_trough(h_trough),'MarkerFaceColor','k')
                idx_trough_flags(h_trough) = 1;
            elseif idx_trough_flags(h_trough) == 1
                set(handles_trough(h_trough),'MarkerFaceColor','b')
                idx_trough_flags(h_trough) = 0;
            end
        end
    end
    
    function pushbutton_callback(source,~)
        switch source
            case finish
                finalpeaks = idx_peak_flags;
                finaltrough = idx_trough_flags;
                close(fig);
        end
    end
end
