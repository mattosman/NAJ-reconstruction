% this function accompanies ccrNested.m
% written by M. Osman (osmanm@mit.edu) May 2019

cd cbrewer
	warning('off','all'); 
    CT1 = cbrewer('seq','Blues' ,17); CT1(CT1<0) = 0; CT1(CT1>1) = 1; 
    CT2 = cbrewer('seq','Reds' ,21); CT2(CT2<0) = 0; CT2(CT2>1) = 1; 
    CT3 = cbrewer('div','RdBu' ,21);  CT3 = flipud(CT3);  CT3(CT3<0) = 0; CT3(CT3>1) = 1; 
    warning('on','all')
cd ../

CVfig1 = figure(101); clf; hold on; 
set(CVfig1,'PaperPositionMode','auto','PaperOrientation','landscape','Position',[50 50 1100 600]); 

% RMSEP
    subplot(3,3,1);
        RMSEP_dp_dj = nan(dp,dj);
            for j = 1:dp
                for k = 1:dj
                    indexer = MODES(:,1) == j & MODES(:,2) == k;
                    RMSEP_dp_dj(j,k) = nanmean(RMSEP(indexer));
                end
            end        
    xlim([1.5 size(RMSEP_dp_dj,1)]);
    ylim([1.5 size(RMSEP_dp_dj,2)]);
    h = imagesc(RMSEP_dp_dj'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat)); 
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Proxy mode'); ylabel('Jet mode');
    colormap(gca,CT2);
    c = colorbar('eastoutside'); c.Label.String = ['RMSEP']; set(c,'Fontsize',12); % only need to do here 
    
   subplot(3,3,2);
        RMSEP_dj_dcca = nan(dj,dcca);
            for j = 1:dj
                for k = 1:dcca
                    indexer = MODES(:,2) == j & MODES(:,3) == k;
                    RMSEP_dj_dcca(j,k) = nanmean(RMSEP(indexer));
                end
            end
    xlim([1.5 size(RMSEP_dj_dcca,1)]);
    ylim([1.5 size(RMSEP_dj_dcca,2)]);
    h = imagesc(RMSEP_dj_dcca'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat)); 
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Jet mode'); ylabel('CCA mode');
    colormap(gca,CT2);
    c = colorbar('eastoutside'); c.Label.String = ['RMSEP']; set(c,'Fontsize',12); % only need to do here     
    
    subplot(3,3,3);
        RMSEP_dcca_dp = nan(dcca,dp);
            for j = 1:dcca
                for k = 1:dp
                    indexer = MODES(:,3) == j & MODES(:,1) == k;
                    RMSEP_dcca_dp(j,k) = nanmean(RMSEP(indexer));
                end
            end
    xlim([1.5 size(RMSEP_dcca_dp,1)]);
    ylim([1.5 size(RMSEP_dcca_dp,2)]); 
    h = imagesc(RMSEP_dcca_dp'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat));     
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('CCA mode'); ylabel('Proxy mode');
    colormap(gca,CT2);
    c = colorbar('eastoutside'); c.Label.String = ['RMSEP']; set(c,'Fontsize',12); % only need to do here 
    
% RE
    subplot(3,3,4);
        RE_dp_dj = nan(dp,dj);
            for j = 1:dp
                for k = 1:dj
                    indexer = MODES(:,1) == j & MODES(:,2) == k;
                    RE_dp_dj(j,k) = nanmean(RE(indexer));
                end
            end    
    xlim([1.5 size(RE_dp_dj,1)]);
    ylim([1.5 size(RE_dp_dj,2)]);
    h = imagesc(RE_dp_dj'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat)); 
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Proxy mode'); ylabel('Jet mode');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['RE']; set(c,'Fontsize',12); % only need to do here 
    lims = [-0.5 0.5]; caxis([lims(1) lims(2)]);
    
    subplot(3,3,5);
        RE_dj_dcca = nan(dj,dcca);
            for j = 1:dj
                for k = 1:dcca
                    indexer = MODES(:,2) == j & MODES(:,3) == k;
                    RE_dj_dcca(j,k) = nanmean(RE(indexer));
                end
            end
    xlim([1.5 size(RE_dj_dcca,1)]);
    ylim([1.5 size(RE_dj_dcca,2)]);
    h = imagesc(RE_dj_dcca'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat)); 
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Jet mode'); ylabel('CCA mode');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['RE']; set(c,'Fontsize',12); % only need to do here 
    lims = [-0.5 0.5]; caxis([lims(1) lims(2)]);

    
    subplot(3,3,6);
        RE_dcca_dp = nan(dcca,dp);
            for j = 1:dcca
                for k = 1:dp
                    indexer = MODES(:,3) == j & MODES(:,1) == k;
                    RE_dcca_dp(j,k) = nanmean(RE(indexer));
                end
            end
    xlim([1.5 size(RE_dcca_dp,1)]);
    ylim([1.5 size(RE_dcca_dp,2)]); 
    h = imagesc(RE_dcca_dp'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat));     
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('CCA mode'); ylabel('Proxy mode');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['RE']; set(c,'Fontsize',12); % only need to do here 
    lims = [-0.5 0.5]; caxis([lims(1) lims(2)]);
      
% CE
    subplot(3,3,7);
        CE_dp_dj = nan(dp,dj);
            for j = 1:dp
                for k = 1:dj
                    indexer = MODES(:,1) == j & MODES(:,2) == k;
                    CE_dp_dj(j,k) = nanmean(CE(indexer));
                end
            end    
    xlim([1.5 size(CE_dp_dj,1)]);
    ylim([1.5 size(CE_dp_dj,2)]);
    h = imagesc(CE_dp_dj'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat)); 
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Proxy mode'); ylabel('Jet mode');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['CE']; set(c,'Fontsize',12); % only need to do here 
    lims = [-0.5 0.5]; caxis([lims(1) lims(2)]);
    
    subplot(3,3,8);
        CE_dj_dcca = nan(dj,dcca);
            for j = 1:dj
                for k = 1:dcca
                    indexer = MODES(:,2) == j & MODES(:,3) == k;
                    CE_dj_dcca(j,k) = nanmean(CE(indexer));
                end
            end
    xlim([1.5 size(CE_dj_dcca,1)]);
    ylim([1.5 size(CE_dj_dcca,2)]);
    h = imagesc(CE_dj_dcca'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat)); 
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Jet mode'); ylabel('CCA mode');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['CE']; set(c,'Fontsize',12); % only need to do here 
    lims = [-0.5 0.5]; caxis([lims(1) lims(2)]);

    
    subplot(3,3,9);
        CE_dcca_dp = nan(dcca,dp);
            for j = 1:dcca
                for k = 1:dp
                    indexer = MODES(:,3) == j & MODES(:,1) == k;
                    CE_dcca_dp(j,k) = nanmean(CE(indexer));
                end
            end
    xlim([1.5 size(CE_dcca_dp,1)]);
    ylim([1.5 size(CE_dcca_dp,2)]); 
    h = imagesc(CE_dcca_dp'); shading flat;
        cdat = get(h,'CData'); 
        set(h,'alphadata',isfinite(cdat));     
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('CCA mode'); ylabel('Proxy mode');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['CE']; set(c,'Fontsize',12); % only need to do here 
    lims = [-0.5 0.5]; caxis([lims(1) lims(2)]);
    
%% RMSEP vs. RE for dp, dj, dcca

CVfig2 = figure(102); clf; hold on; 
set(CVfig2,'PaperPositionMode','auto','PaperOrientation','landscape','Position',[50 50 1100 250]); 

cd cbrewer
	warning('off','all'); 
    CT1 = cbrewer('div','RdYlBu' ,nanmax(MODES(:,1)));  CT1 = flipud(CT1); CT1(CT1<0) = 0; CT1(CT1>1) = 1; 
    CT2 = cbrewer('div','RdYlBu' ,nanmax(MODES(:,2)));  CT2 = flipud(CT2); CT2(CT2<0) = 0; CT2(CT2>1) = 1; 
    CT3 = cbrewer('div','RdYlBu' ,nanmax(MODES(:,3)));  CT3 = flipud(CT3); CT3(CT3<0) = 0; CT3(CT3>1) = 1; 
    warning('on','all')
cd ../

alph = 0.5;
scatter_size = 20;
subplot(1,3,1);% proxy mode
    s = scatter(RE,RMSEP,scatter_size,MODES(:,1),'Filled','MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph); shading flat; hold on;
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('RE'); ylabel('RMSEP');
    colormap(gca,CT1);
    c = colorbar('eastoutside'); c.Label.String = ['Proxy mode']; set(c,'Fontsize',12); % only need to do here 
   
subplot(1,3,2);% proxy mode
    s = scatter(RE,RMSEP,scatter_size,MODES(:,2),'Filled','MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph); shading flat; hold on;
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('RE'); ylabel('RMSEP');
    colormap(gca,CT2);
    c = colorbar('eastoutside'); c.Label.String = ['Jet mode']; set(c,'Fontsize',12); % only need to do here 

subplot(1,3,3);% proxy mode
    s = scatter(RE,RMSEP,scatter_size,MODES(:,3),'Filled','MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph); shading flat; hold on;
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('RE'); ylabel('RMSEP');
    colormap(gca,CT3);
    c = colorbar('eastoutside'); c.Label.String = ['CCA mode']; set(c,'Fontsize',12); % only need to do here 
  
    
%% RE vs. proxy, jet, and cca

CVfig3 = figure(103); clf; hold on; 
% set(CVfig3,'PaperPositionMode','auto','PaperOrientation','landscape','Position',[50 50 250 250]); 

cd cbrewer
	warning('off','all'); 
    CT1 = cbrewer('div','RdYlBu' ,101);  CT1 = flipud(CT1);  CT1(CT1<0) = 0; CT1(CT1>1) = 1; 
    warning('on','all')
cd ../

alph = 0.85;
scatter_size = 100;
subplot(1,1,1); grid on;
    s = scatter3(gca,MODES(:,1),MODES(:,2),MODES(:,3),scatter_size,RE,'Filled','MarkerFaceAlpha',alph,'MarkerEdgeAlpha',alph); shading flat; hold on;
    set(gca,'Linewidth',1.5,'Box','on','Fontsize',12,'Color','none')
    xlabel('Proxy mode'); ylabel('Jet mode'); zlabel('CCA mode'); 
    colormap(gca,CT1);
    c = colorbar('eastoutside'); c.Label.String = ['RE']; set(c,'Fontsize',12); % only need to do here 
    view(40,35)

%% Remove some variables

clearvars CT1 CT2 CT3 ...
    RMSEP_dp_dj RMSEP_dj_dcca RMSEP_dcca_dp ...
    RE_dp_dj RE_dj_dcca RE_dcca_dp ...
    CE_dp_dj CE_dj_dcca CE_dcca_dp ...
    h c s
