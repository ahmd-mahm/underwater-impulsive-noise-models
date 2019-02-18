function H=absorption_filter_range(range_max,bins,order,fs,d)

% Evaluates absorption filters (each of length 'order+1') for
% 'range_max/(2*bins)' up until 'range_max(2*bins-1)/(2*bins)' with
% 'range_max/bins' increments. These are returned as the columns of 'H'.

% range_max:    maximum range in meters from the source
% bins:         number of absorption filters required for ranges within [0,range_max]
% order:        order of the filter
% fs:           sampling frequency in Hz
% d:            depth of the hydrophone in meters

edges=0:range_max/bins:range_max;
range_bins=(edges(1:end-1)+edges(2:end))/2;

H=zeros(order+1,bins);
for i=1:bins
    H(:,i) = absorptionFilter(order,fs,range_bins(i),d).';
end