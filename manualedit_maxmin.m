function [finalpeaks,finaltrough,ptl_max,ptl_min,pth_max,pth_min] ...
    = manualedit_maxmin(wave,idx_peak,idx_trough,srng)
close all

idx_peak_flags = idx_peak*0;
idx_trough_flags = idx_trough*0;

fig = figure;
hold all;
grid on;
plot(wave)
handles_peak = gobjects(1,length(idx_peak));

ptl_max = cell(1,2); ptl_min = cell(1,2);
pth_max = cell(1,2); pth_min = cell(1,2);

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

if srng == 'lp'
    [xl_max,yl_max] = getpts;
    ptl_max{1,1} = xl_max;
    ptl_max{1,2} = yl_max;
elseif srng == 'hp'
    [xh_max,yh_max] = getpts;
    pth_max{1,1} = xh_max;
    pth_max{1,2} = yh_max;
end

if srng == 'lp'
    [xl_min,yl_min] = getpts;
    ptl_min{1,1} = xl_min;
    ptl_min{1,2} = yl_min;
elseif srng == 'hp'
    [xh_min,yh_min] = getpts;
    pth_min{1,1} = xh_min;
    pth_min{1,2} = yh_min;
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
