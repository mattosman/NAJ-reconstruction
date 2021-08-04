% Driver function to reconstruct the N. Atlantic Jet Stream from
% Greenlandic ice proxies.  
% Loads in and processes the observed NAJ targets (ERA or NOAA), Greeneland 
% proxies, runs the reconstruction (CCR_nested.m), and saves the output.
% Written by M. Osman (osmanm@mit.edu), May 2019
% See Osman et al. (2021, PNAS) for details. 
clear; 
addpath(genpath('ccr/')); 

%% USER INPUTS REQUIRED HERE!  User defined reconstruction params:

oldYr = 2000; % oldest year to retain for proxies
yngYr = 700; % youngest year to retain for proxies
calibYr = 1900; % youngest calibration year (oldYr assumed oldest calib year)
smoother = 1; % if 1, then no smoothing applied
target = "NOAA"; % input either "NOAA" or "ERA"
dtJet = false; % option to linearly detrend jet stream 
dtProx = false; %  option to linearly detrend proxies
proxInflVers = "eof"; % choose either "ppca", "eof", "mean", or "iwd"
crossValSig = false; % if true, computes cross validation significance tests (caution, is SLOW!)
saveOutput = true; % option to save output
plotOpt = true; % option to plot the reconstruction

%% Load in target data

% load in target data
    if target == "NOAA"
        cd target; load NOAA20C_NAJ_ann.mat; cd ../
    elseif target == "ERA"
        cd target; load ERA20C_NAJ_ann.mat; cd ../
    else 
        warning('Input either "NOAA" or "ERA" for the ''target'' variable. Try again.');
        return;
    end
% option to detrend
    if dtJet
        for i = 1:size(jet_profile,2)
            jet_profile(:,i) = detrend(jet_profile(:,i),'linear') + nanmean(jet_profile(:,i));
        end     
    end
    
%% Load in proxy data

% load in proxy data
    cd proxy    
        if strcmp(proxInflVers,"ppca")
            load data_ppca.mat
        elseif strcmp(proxInflVers,"eof")
            load data_eof.mat
        elseif strcmp(proxInflVers,"iwd")
            load data_iwd.mat
        elseif strcmp(proxInflVers,"mean")
            load data_mean.mat
        else
            warning('Input either "ppca", "eof", "iwd", or "mean" for the ''proxVers'' variable. Try again.');
            return;
        end
    cd ../
% clip data to specified yeears
    I = data.year <= oldYr & data.year >= yngYr; 
    data.data = data.data(I,:); 
    data.year = data.year(I); 
    % Option to remove remove industrial era secular trend:
    if dtProx 
        cd process
        [data] = removeAnthro(data,2); % always use version == 2
        cd ../
    end   
% remove fliers
    I = data.year >= 1775 & data.year <= oldYr; % 1775 is year of greatest overlaps
% standardize data to common overlap
	for i = 1:size(data.data,2)
        data.mu(i) = nanmean(data.data(I,i));
        data.sigma(i) = nanstd(data.data(I,i));
        data.data(:,i) = (data.data(:,i) - data.mu(i))./data.sigma(i);
    end
% Remove extreme fliers
	for i = 1:size(data.data,2)
        indexer = ~isnan(data.data(:,i));
        cd process
        data.data(indexer,i) = filterAnoms(data.data(indexer,i), 101, 6); % <--- CHANGE here for anoms, 6 is VERY conservative
        cd ../
    end
% return variance back into data
	for i = 1:size(data.data,2)
        data.data(:,i) = data.data(:,i).*data.sigma(i) + data.mu(i);
	end        

%% Define X and Y variables for reconstruction (see Osman et al., 2021, PNAS, Supporting Information)

X.data = data.data;
X.year = data.year;
Y.data = jet_profile;
Y.year = year;
Y.lat = lat;
Y.var = var(Y.data,0,1);

if smoother ~= 1
    for i = 1:size(X.data,2)
        indexer = ~isnan(X.data(:,i));
        X.data(indexer,i) = lowpass(X.data(indexer,i),1/smoother,1,1);
    end
    for i = 1:size(Y.data,2)
        Y.data(:,i) = lowpass(Y.data(:,i),1/smoother,1,1);
    end
end

% standardize relative to calibration interval - save mean + stdev for later
calib_int = X.year >= calibYr;
for i = 1:size(X.data,2)
    X.mu(i) = nanmean(X.data(calib_int,i));
    X.sigma(i) = nanstd(X.data(calib_int,i));
    X.data(:,i) = (X.data(:,i) - X.mu(i))./X.sigma(i);
end
for i = 1:size(Y.data,2)
    Y.mu(i) = nanmean(Y.data(:,i));
    Y.sigma(i) = nanstd(Y.data(:,i));
    Y.data(:,i) = (Y.data(:,i) - Y.mu(i))./Y.sigma(i);
end

% get area per lat (for a delta lat x 1deg lon profile)
earth_ellipsoid = referenceSphere('earth','km');
delta_lat = abs(Y.lat(2) - Y.lat(1));
Y.area_per_lat = zeros(1,length(Y.lat));
for j = 1:length(Y.lat)
	Y.area_per_lat(j) = areaquad((Y.lat(j)-delta_lat/2),0,(Y.lat(j)+delta_lat/2),1,earth_ellipsoid);
end

%% Run the CCR

[yrec, CCA, pctVar, CrossVal, TimeInt, pctVarp, pctVarj] = ccrNested(X,Y,1,'Aut',crossValSig); % Main reconstruction function

%% Option to save output
if saveOutput
    cd output
    file = strcat('CCA_',target,'_',num2str(yngYr),'-',num2str(oldYr),'CE_',proxInflVers,'ProxyInfill_',date,'.mat');
    save(file,'yrec', 'CCA', 'pctVar', 'CrossVal', 'TimeInt', 'pctVarp', 'pctVarj');
    cd ../
end

%% Option to plot jet lat / speed diagnostics

if plotOpt
    
    % assign colors
    cd cbrewer
    warning('off','all')
        CT1 = cbrewer('seq','Reds' ,7); 
        CT2 = cbrewer('seq','Blues' ,5);
        CT3 = cbrewer('seq','Purples' ,5);    
        CT4 = cbrewer('seq','Greys' ,5);    
    warning('on','all')
    cd ../

    cd ccr
        [JetLat_obs, JetSpeed_obs] = upscaleJetLat(yrec.obs,lat,5,0.01,true);
        [JetLat_rec, JetSpeed_rec] = upscaleJetLat(yrec.rec,lat,5,0.01,true);
    cd ../

    % Jet latitude
    h = figure; hold on;
        set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
        h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.50*Pix_SS(3),.30*Pix_SS(4)];
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');

    set(gca,'Color','none','Linewidth',1.5,'Fontsize',11);
    xlabel('Year (CE)'); 
    yearz = yrec.year;
        for i = 1:size(TimeInt,1)
            indexer = yearz >= nanmin(TimeInt(i,:)) & yearz <= nanmax(TimeInt(i,:));
            upper_lat(indexer) = JetLat_rec(indexer) + CrossVal{i}.CVstats_mu_jetlat.RMSEP;
            lower_lat(indexer) = JetLat_rec(indexer) - CrossVal{i}.CVstats_mu_jetlat.RMSEP;
        end
        h = fill([yearz(1);  yearz;  flipud([(yearz);  yearz(end)])],[upper_lat(1);  smooth(lower_lat,1);  flipud([smooth(upper_lat,1);  lower_lat(end)])],CT2(end-1,:));
        set(h,'edgecolor','none','facealpha',0.3);
        f2 = plot(yrec.year_obs, JetLat_obs,'Color',CT1(end-2,:),'linewidth',0.5); f2.Color(4) = 0.9;
        f1 = plot(yrec.year, JetLat_rec,'Color',CT2(end-1,:),'linewidth',1); f1.Color(4) = 0.9;
            ylabel('Position (^{\circ}N)');
            legend([f2,f1,h],'Observed NAJ position','Reconstructed NAJ position','Reconstruction uncertainty (\pm1\sigma)','Orientation','vertical','box','off');

    % Jet speed
    h = figure; hold on;
        set(0,'units','pixels'); Pix_SS = get(0,'screensize'); % PIX_SS - 3rd val = width of screen, 4th val = height of screen
        h.Position = [.05*Pix_SS(3),.05*Pix_SS(4),.50*Pix_SS(3),.30*Pix_SS(4)];
    set(gcf,'PaperPositionMode','auto');         
    set(gcf,'PaperOrientation','landscape');

    set(gca,'Color','none','Linewidth',1.5,'Fontsize',11);
    xlabel('Year (CE)'); 
    yearz = yrec.year;
        for i = 1:size(TimeInt,1)
            indexer = yearz >= nanmin(TimeInt(i,:)) & yearz <= nanmax(TimeInt(i,:));
            upper_speed(indexer) = JetSpeed_rec(indexer) + CrossVal{i}.CVstats_mu_jetspeed.RMSEP;
            lower_speed(indexer) = JetSpeed_rec(indexer) - CrossVal{i}.CVstats_mu_jetspeed.RMSEP;
        end
        h = fill([yearz(1);  yearz;  flipud([(yearz);  yearz(end)])],[upper_speed(1);  smooth(lower_speed,1);  flipud([smooth(upper_speed,1);  lower_speed(end)])],CT2(end-1,:));
        set(h,'edgecolor','none','facealpha',0.3);
        f2 = plot(yrec.year_obs, JetSpeed_obs,'Color',CT1(end-2,:),'linewidth',0.5); f2.Color(4) = 0.9;
        f1 = plot(yrec.year, JetSpeed_rec,'Color',CT2(end-1,:),'linewidth',1); f1.Color(4) = 0.9;
            ylabel('Intensity (m s^{-1})');
            legend([f2,f1,h],'Observed NAJ intensity','Reconstructed NAJ intensity','Reconstruction uncertainty (\pm1\sigma)','Orientation','vertical','box','off');

end