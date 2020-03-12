% This function takes a data series with a fine temporal scale and
% aggregates it to a more coarse temporal scale.
%
% RELEASE NOTES
%   Written by Mark Raleigh (mraleig1@uw.edu) June 2009
%   Version 1.1 by Mark Raleigh - allowed multiple aggregation options for a single station
%
% SYNTAX
%   [coarse_data, coarse_time] = aggMET(fine_DATA, fine_TIME, fine_timestep, coarse_timestep, yr_base, month_base, day_base, hr_base, tzone, agg_period, agg_option, tolerance)
%
% INPUTS
% % % fine_DATA = NxM matrix of fine temporal scale, where N = # of fine timesteps, M = # of stations
% % % fine_TIME = Nx7 matrix (built from timebuilder.m), where:
%           Column 1: Year
%           Column 2: Month
%           Column 3: Gregorian Day
%           Column 4: Hour
%           Column 5: Minute
%           Column 6: Julian Day
%           Column 7: Matlab Serial Date
% % % fine_timestep = 1x1 value, fine timestep (hrs)
% % % coarse_timestep = 1x1 value, coarse timestep (hrs)
% % % yr_base = 1x1 value, year of first coarse timestep (e.g., 2005)
% % % month_base = 1x1 value, month of first coarse timestep (e.g., 10)
% % % day_base = 1x1 value, day of first coarse timestep (e.g., 15)
% % % hr_base = 1x1 value, hour of first coarse timestep (e.g., 0)
% % % tzone = 1x1 value, time zone shift from input to output (e.g.,, enter -8 going from UTC/GMT data to PST data, +8 when going from PST to GMT)
% % % agg_period = data aggregation period at each time step, where:
%     % 1 = Period before the current time step (e.g., the hourly 11 AM data point would represent the aggregation between 10 and 11 am)
%     % 2 = Halfway before the current time step to half way after the current time step (e.g., 11 am would represent the period from 10:30 to 11:30)
%     % 3 = Snapshot at the current time step (11 am represents only the data taken at 11 am)
%     % 4 = Period after the current time step (e.g., the hourly 11 AM data point would represent the aggregation between 11 am and 12 pm.)
% % % agg_option = data aggregating option (when raw data with finer temporal resolution exists), where:
%     % 1 = Average over agg_period
%     % 2 = Maximum over agg_period
%     % 3 = Minimum over agg_period
%     % 4 = Median over agg_period
%     % 5 = Cumulative sum over agg_period
% % tolerance = fraction of fine temporal record allowed to be missing
%
% OUTPUTS
% % % coarse_data = nxM matrix (n<N) of coarse temporal scale, where n = # of coarse timesteps, and M = # of stations or M = different agg_options (if single station)
% % % coarse_time = nx7 matrix of time values, same column designations as fine_TIME
%
% SCRIPTS REQUIRED
%   time_builder.m, time_shift.m



function [coarse_data, coarse_time] = aggMET(fine_DATA, fine_TIME, fine_timestep, coarse_timestep, yr_base, month_base, day_base, hr_base, tzone, agg_period, agg_option, tolerance)


%% Checks

M= size(fine_DATA,2);

if fine_timestep<=0
    error('Timestep of fine resolution data must be greater than 0')
end

if fine_timestep > coarse_timestep
    error('fine_timestep must be less than coarse time_step')
end

if size(fine_DATA,1) ~= size(fine_TIME,1)
    error('fine_DATA and fine_TIME must have the same number of rows')
end

if size(fine_TIME,2) ~= 7
    error('fine_TIME must have 7 columns, as designated by timebuilder.m')
end


if datenum(yr_base, month_base, day_base) > nanmax(fine_TIME(:,7))
    disp('Warning: Aggregation starts after data.  No values will result')
    badt = 1;
else
    badt = 0;
end

if numel(tzone)~=1
    error('Invalid time zone shift')
end

if M==1
    if numel(agg_period)~=numel(agg_option)
        error('agg_period must have same number of inputs as agg_option')
    else
        if numel(agg_period)>1
            multagg = 1;
            M= numel(agg_period);
        else
            multagg = 0;
        end
    end
else
    if numel(agg_period)~=1 || numel(agg_option)~=1
        error('agg_period and agg_option must be single values (multiple station aggregation)')
    end
    multagg = 0;
end


%% Code

if tzone~=0
    fine_TIME = time_shift(fine_TIME, tzone);
end



[f_yr, f_mo, f_da, f_hr, ~, ~] = datevec(fine_TIME(end,7));

tolerance = 1-tolerance;

if badt==0
    coarse_time = time_builder(yr_base, month_base, day_base, hr_base, 0, f_yr, f_mo, f_da, f_hr, 0, coarse_timestep);
    
    %%%% added this in March 2019
    if agg_period==1
        ctime_f = datenum(yr_base, month_base, day_base, hr_base, 0, 0);
        if ((ctime_f-coarse_time(end,7)) - (coarse_timestep/24)) <= 0.0001
            ctime_f = coarse_time(end,7) + (coarse_timestep/24);
            coarse_time = [coarse_time; time_builder(ctime_f)];
        end
    end
    %%%%
    
    n = size(coarse_time,1);
    coarse_data = zeros(n, M)*NaN;
else
    coarse_time = [];
    coarse_data = [];
end

if badt==0
    for i=1:n
        clear time1 time2 data_sub
        
        if agg_period==1
            time1 = coarse_time(i,7) - (coarse_timestep/24);
            time2 = coarse_time(i,7);
            data_sub = find(fine_TIME(:,7) <= time2 & fine_TIME(:,7) >time1);
        elseif agg_period==2
            time1 = coarse_time(i,7) - (coarse_timestep/48);
            time2 = coarse_time(i,7) + (coarse_timestep/48);
            data_sub = find(fine_TIME(:,7) < time2 & fine_TIME(:,7) >=time1);
        elseif agg_period==3
            data_sub = find(fine_TIME(:,7) == coarse_time(i,7), 1, 'first');
        elseif agg_period==4
            time1 = coarse_time(i,7);
            time2 = coarse_time(i,7) + (coarse_timestep/24);
            data_sub = find(fine_TIME(:,7) < time2 & fine_TIME(:,7) >=time1);
        end
        
        
        for j=1:M
            
            if multagg==1
                ma = j;
                jjj = 1;
            else
                ma = 1;
                jjj = j;
            end
            
            if isempty(data_sub)==0
                
                if agg_period(ma)~=3
                    
                    Rnum = nancount(fine_DATA(data_sub,jjj));
                    
                    Rtot = time1:(fine_timestep/24):time2;
                    if Rtot(end) ~= time2
                        Rtot = horzcat(Rtot, time2);
                    end
                    if (Rnum/length(Rtot)) <= tolerance
                        lock =1;    % then there are too much missing data
                    else
                        lock =0;    % then there are enough data
                    end
                else
                    lock = 0;
                end
                
                if lock==0
                    if agg_option(ma)==1
                        coarse_data(i,j) = nanmean(fine_DATA(data_sub,jjj));
                    elseif agg_option(ma)==2
                        coarse_data(i,j) = nanmax(fine_DATA(data_sub,jjj));
                    elseif agg_option(ma)==3
                        coarse_data(i,j) = nanmin(fine_DATA(data_sub,jjj));
                    elseif agg_option(ma)==4
                        coarse_data(i,j) = nanmedian(fine_DATA(data_sub,jjj));
                    elseif agg_option(ma)==5
                        coarse_data(i,j) = nansum(fine_DATA(data_sub,jjj));
                    end
                end
            end
            
        end
    end
end
