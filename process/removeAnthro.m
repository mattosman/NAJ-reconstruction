% this function removes the anthropogenic warming secular trends from input
% data structure, which contains at minimum the variables "year" and "data" 
% (see also load_GrIS_data.m from which this function should follow)
% WRITTEN BY M OSMAN (Feb2020; mattosman@arizona.edu)
% 
% INPUTS:
% data : structure containing the variables data and year.  Note that year
%   must extend back to the year 1831 (onset of sustained industrial era 
% arctic warming via Abram et al., 2016, Nature), and beyond the year 1950.
% version : enter either 1, to remove warming trends using overlapping SW
%   Greenland temperature trends (Cappelan, 2018, DMI) or 2 to remove all
%   secular trends in proxies, regardless of magnitude/direction. 
% OUTPUTS:
% data : same structure, with Industrial Era secular trends now removed.
function [data] = removeAnthro(data,version)

home_fold = cd;

% make some checks
if ~isfield(data,'data')
    warning('Make sure the ''data'' structure input contains the variable ''data''. Try again!'); return; 
elseif ~isfield(data,'year')
    warning('Make sure the ''data'' structure input contains the variable ''year''. Try again!'); return; 
end

if nanmin(data.year) > 1831
    warning('Make sure the ''data'' structure''s oldest year is older than A.D. 1831. Try again!'); return; 
elseif nanmax(data.year) < 1950
    warning('Make sure the ''data'' structure''s youngest year is younger than A.D. 1950. Try again!'); return; 
end

if version ~= 1 && version ~= 2
    warning('Input for ''version'' must equal 1 or 2. Try again!'); return; 
end

%% Remove secular trends

if version == 2
    % standardize data relative to industrial year onset
    indexer = data.year >= 1831 & data.year <= nanmax(data.year);
    for i = 1:size(data.data,2)
        data.avgVal(1,i) = nanmean(data.data(indexer,i));
        data.stdVal(1,i) = nanstd(data.data(indexer,i));
        data.data(:,i) = (data.data(:,i) - data.avgVal(i))./data.stdVal(i);
    end
    % determine + remove trend for each record
    for i = 1:size(data.data,2)
        m = polyfit(data.year(indexer),data.data(indexer,i),1);
        fittedData = [m(1).*data.year(indexer) + m(2)];
        data.data(indexer,i) = data.data(indexer,i) - fittedData + fittedData(end);
        data.data(:,i) = data.data(:,i) - nanmean(data.data(:,i)); % adjust mean bias back to zero
    end
    % convert records back to their original mean and variance:
    for i = 1:size(data.data,2)
        data.data(:,i) = data.data(:,i).*data.stdVal(i) + data.avgVal(i);
    end
    data = rmfield(data,'stdVal'); data = rmfield(data,'avgVal');
elseif version == 1    
    %load in the SW greenland temps
    cd('/Users/matthewosman/Documents/MATLAB/Nuus_project/Nuus_temp_reconstruction/DMI_temps/');
        docdata = xlsread('SW_Greenland_temps.xlsx','SW_Greenland');
    cd(home_fold);
    series = 1:183; % back to 
	year = docdata(series,1);                  
	monthly_temp = docdata(series,2:13);   
    % index temperature data relative to overlapping period with proxies
    indexerTemp = year >= 1831 & year <= nanmax(data.year);
    monthly_temp = monthly_temp(indexerTemp,:);
    year = year(indexerTemp);
    for i = 1:size(monthly_temp,2)
        avgTemp(1,i) = nanmean(monthly_temp(:,i));
        stdTemp(1,i) = nanstd(monthly_temp(:,i));
        monthly_temp(:,i) = (monthly_temp(:,i) - avgTemp(i))./stdTemp(i);
    end
    % infill using ppca
    monthly_temp_preInfilled = monthly_temp;
    cd('/Users/matthewosman/Documents/MATLAB/ppca')
        [~, ~, ~, monthly_temp] = ppcaMO(monthly_temp, size(monthly_temp,2)-1);
    cd(home_fold);    
    % now, need to remove anomalies:
    AnomThresh = 2.5;
    for i = 1:size(monthly_temp,2)
    	monthly_temp(:,i) = filter_anoms(monthly_temp(:,i), 31, AnomThresh);
    end
    m = polyfit(year,nanmean(monthly_temp,2),1);
    fittedData = [m(1).*year + m(2)];
    % NOW, deal with proxies
    % standardize data relative to industrial year onset
    indexer = data.year >= 1831 & data.year <= nanmax(data.year);
    % standardize data
    for i = 1:size(data.data,2)
        data.avgVal(1,i) = nanmean(data.data(indexer,i));
        data.stdVal(1,i) = nanstd(data.data(indexer,i));
        data.data(:,i) = (data.data(:,i) - data.avgVal(i))./data.stdVal(i);
    end    
    % remov temperature component
    for i = 1:size(data.data,2)
        data.data(indexer,i) = data.data(indexer,i) - fittedData + fittedData(end);
        data.data(:,i) = data.data(:,i) - nanmean(data.data(:,i)); % adjust mean bias back to zero  
    end
    % convert records back to their original mean and variance:
    for i = 1:size(data.data,2)
        data.data(:,i) = data.data(:,i).*data.stdVal(i) + data.avgVal(i);
    end
    data = rmfield(data,'stdVal'); data = rmfield(data,'avgVal');    
end

%%
%     % standardize data relative to industrial year onset
%     indexer = data.year >= 1831 & data.year <= nanmax(data.year);
%     for i = 1:size(data.data,2)
%         data.avgVal(1,i) = nanmean(data.data(indexer,i));
%         data.stdVal(1,i) = nanstd(data.data(indexer,i));
%         data.data(:,i) = (data.data(:,i) - data.avgVal(i))./data.stdVal(i);
%     end
%     % determine + remove trend for each record
%     for i = 1:size(data.data,2)
%         m = polyfit([0:(sum(indexer)-1)]',data.data(indexer,i),1);
%         data.data(indexer,i) = data.data(indexer,i) - [m(1).*[0:(sum(indexer)-1)]' + m(2)] - m(2);
%         data.data(:,i) = data.data(:,i) - nanmean(data.data(:,i)); % adjust mean bias back to zero
%     end
%     % convert records back to their original mean and variance:
%     for i = 1:size(data.data,2)
%         data.data(:,i) = data.data(:,i).*data.stdVal(i) + data.avgVal(i);
%     end
%     data = rmfield(data,'stdVal'); data = rmfield(data,'avgVal');

end
