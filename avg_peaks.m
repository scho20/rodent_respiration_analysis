function [avgpeaks_LP, avgpeaks_HP] = avg_peaks(final_locs_l_max, final_locs_h_max, timeframe, xl_filt_offset, xh_filt_offset, k)

peak_idx_LP = final_locs_l_max{1,k};
peak_idx_HP = final_locs_h_max{1,k};

avgpeaks_LP = zeros(1,2*timeframe+1);
avgpeaks_HP = zeros(1,2*timeframe+1);

for i = 1:size(peak_idx_LP)
    avgpeaks_LP = xl_filt_offset(k,peak_idx_LP(1,i)-timeframe:peak_idx_LP(1,i)+timeframe); 
end

for i = 1:size(peak_idx_HP)
    avgpeaks_HP = xh_filt_offset(k,peak_idx_HP(1,i)-timeframe:peak_idx_HP(1,i)+timeframe);
end

end