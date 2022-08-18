% 04_ConstrainCO2Flux_genTrawlSeqFracs.m
% Created June 7, 2022 by Jamie Collins, jcollins@edf.org, based on the script 
% plot_sequestration_fraction.m provided by Siegel et al. at
%  https://doi.org/10.6084/m9.figshare.15228690.v2

% Note: Assumes user has already run scripts 1-3 in this series:
% '01_ConstrainCO2Flux_genSiegelMetaData.m'
% '02_ConstrainCO2Flux_IO.R'
% '03_ConstrainCO2Flux_coordMatch.R'
% ... and has loaded into memory two files:
% 'fseq_OCIM2_48L.mat' (unmodified Siegel et al. OCIM48 model output)
% 'Sala_CO2_efflux_nonZero.csv' (generated in 03_ConstrainCO2Flux_coordMatch.R)

% variables contained in fseq_OCIM2_48L.mat are as follows
% 1. fseq = the fraction of carbon remaining sequestered at each ocean grid
% cell as a function of time since injection; dimensions (m x n) where m is
% the number of ocean grid cells and n is the number of years since
% injection
% 2. MASK = 3-dimensional land-sea mask for the model, with 1 = ocean and 0
% = land;dimensions (ny x nx x nz) where ny = 91 (number of grid points in
% latitudinal direction), nx = 180 (number of grid points in the
% longitudinal direction), and nz = 48 (number of grid points in the
% vertical direction)
% 3. LAT = 3-d array of the latitudes (degrees north) of the model grid
% cells; same dimensions as MASK
% 4. LON = 3-d array of the longitudes (degrees east) of the model grid
% cells; same dimensions as MASK
% 5. DEPTH = 3-d array of the depths (meters) of the model grid cells; same
% dimensions as MASK
% 6. VOL = 3-d array of grid box volumes (m^3) of the model grid cells; same
% dimensions as MASK
% 7. AREA = 3-d array of grid box areas (m^2) of the model grid cells; same
% dimensions as MASK
% 8. time = 1-d array of time (in years) since injection of CO2; range from
% 0-1000 years

% first load in 'Sala_CO2_efflux_nonZero.csv', a very large matrix generated in the
% previous step using '03_ConstrainCO2Flux_coordMatch.R'

Sala_CO2_efflux_nonZero_raw = csvread('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/output/Sala_CO2_efflux_nonZero.csv',1,0);

% % Example 1: Find the near-bottom ocean grid cells and plot the fraction of
% % CO2 remaining sequestered after 1 year
% 
% % make a 3-d array of the fraction sequestered after 1 year
% [ny,nx,nz] = size(MASK);
% FSEQ_1yr = 0*MASK+NaN;
% t_indx = find(time==1); % find the index of the 1st year
% FSEQ_1yr(MASK==1) = fseq(:,t_indx);
% % now find the near-bottom ocean grid cells
% fseq_bottom_1yr = 0*MASK(:,:,1)+NaN;
% TOPO = sum(MASK,3); % number of grid cells in the vertical direction
% for i = 1:ny
%   for j = 1:nx
%     if TOPO(i,j)~=0
%       fseq_bottom_1yr(i,j) = FSEQ_1yr(i,j,TOPO(i,j));
%     end
%   end
% end

% objective: make matrix of the fractions sequestered at the ~ bottom depth
% at the best-fit location for each non-zero value in the Sala et al dataset

% timescale: after 1-200 (inclusive) years, plus years 300, 400, ... 900, 1000
fseq_bottom_multyears = zeros([size(Sala_CO2_efflux_nonZero_raw,1), 208]);
years = [1:200 300 400 500 600 700 800 900 1000];
% return indices for the times we want
[tf,idx_time] = ismember(years,time);

for i = 1:size(Sala_CO2_efflux_nonZero_raw,1)
    lat = Sala_CO2_efflux_nonZero_raw(i,7); % degrees north
    lon = Sala_CO2_efflux_nonZero_raw(i,6);
    globalInd = Sala_CO2_efflux_nonZero_raw(i,5);
    % this will be in eastings & westings, so will have to check and
    % reconvert to degrees east, from 0 - 360, if necessary
    if lon < 0
        lon = lon+360;
    end
    bottomDepth = -Sala_CO2_efflux_nonZero_raw(i,2);
    % now find the right OCIM48 grid cell indices (all depths) for a 
    % given Sala data point j
    % find the equivalent 2-D index, adjusted for the mask
    % generate indices to fseq to all depths for this point location
    fseq_ind_thislocation = NaN(size(MASK,3),1);
    masks_thislocation = MASK(globalInd:91*180:globalInd+91*180*47);
    for j = 1:size(fseq_ind_thislocation)
        if masks_thislocation(j)~=1
            break
        end
        fseq_ind_thislocation(j) = sum(MASK(1:((globalInd+(j-1)*(91*180)))));
    end
    % now, find the closest depth at this point for which there is model
    % output
    [~,r]=sort(abs(bottomDepth-squeeze(DEPTH(1,1,:))),'ascend');
    ranked_fseqs = fseq_ind_thislocation(r);
    ranked_fseqs_noNaN = ranked_fseqs(~isnan(ranked_fseqs));
    fseq_thislocation = fseq(ranked_fseqs_noNaN(1),idx_time);
    fseq_bottom_multyears(i,:) = fseq_thislocation;
end

% export the fseq_bottom_multyears array as a series of .mat files; also export a
% matrix containing the years for which we retrieved data as .csv

writematrix(years,'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/trawlYears.csv')

split_ind = round([1:size(fseq_bottom_multyears,1)/7:size(fseq_bottom_multyears,1)]);
fseq_bottom_multyears_1of6 = fseq_bottom_multyears(split_ind(1):split_ind(2),:);
fseq_bottom_multyears_2of6 = fseq_bottom_multyears(split_ind(2)+1:split_ind(3),:);
fseq_bottom_multyears_3of6 = fseq_bottom_multyears(split_ind(3)+1:split_ind(4),:);
fseq_bottom_multyears_4of6 = fseq_bottom_multyears(split_ind(4)+1:split_ind(5),:);
fseq_bottom_multyears_5of6 = fseq_bottom_multyears(split_ind(5)+1:split_ind(6),:);
fseq_bottom_multyears_6of6 = fseq_bottom_multyears(split_ind(6)+1:split_ind(7),:);
save('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_1of6.mat','fseq_bottom_multyears_1of6')
save('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_2of6.mat','fseq_bottom_multyears_2of6')
save('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_3of6.mat','fseq_bottom_multyears_3of6')
save('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_4of6.mat','fseq_bottom_multyears_4of6')
save('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_5of6.mat','fseq_bottom_multyears_5of6')
save('/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_bottom_multyears_6of6.mat','fseq_bottom_multyears_6of6')

% % plot both the fraction remaining after 100 years and the bottom depth
% figure(1)
% subplot(2,1,1)
% pcolor(LON(:,:,1),LAT(:,:,1),fseq_bottom_100yr)
% set(gca,'clim',[0 1])
% colormap(parula(10))
% colorbar
% set(gca,'color',[.8 .8 .8])
% set(gcf,'color','w')
% shading('flat')
% title('Fraction of CO2 injected at seafloor remaining after 100 years')
% subplot(2,1,2)
% pcolor(LON(:,:,1),LAT(:,:,1),bottom_depth)
% set(gca,'clim',[0 5500])
% colormap(parula(11))
% colorbar
% set(gca,'color',[.8 .8 .8])
% set(gcf,'color','w')
% shading('flat')
% title('Depth of seafloor (m)')

% Example 2: Find a particular point in the ocean and plot the fraction
% remaining sequestered over time
% first specify a location by latitude, longitude, and depth
lat = 83.0769; % degrees north
lon = 333; % degrees east, from 0 - 360
depth = 1000; % meters
% the following code will find the nearest model grid cell(s)
iy = find(abs(LAT(:)-lat)==min(abs(LAT(:)-lat)));
ix = find(abs(LON(:)-lon)==min(abs(LON(:)-lon)));
iz = find(abs(DEPTH(:)-depth)==min(abs(DEPTH(:)-depth)));
indx = intersect(intersect(ix,iy),iz);
% check if that grid cell(s) is land or ocean; if in ocean, plot fseq
if sum(MASK(indx))==0
  fprintf('Location is not in the ocean. Please choose another location.\n')
else
  % if in ocean, find location in ocean grid coordinates and plot
  iy2 = find(abs(LAT(MASK==1)-lat)==min(abs(LAT(MASK==1)-lat)));
  ix2 = find(abs(LON(MASK==1)-lon)==min(abs(LON(MASK==1)-lon)));
  iz2 = find(abs(DEPTH(MASK==1)-depth)==min(abs(DEPTH(MASK==1)-depth)));
  indx2 = intersect(intersect(ix2,iy2),iz2);
  figure(2)
  plot(time,mean(fseq(indx2,:)))
  xlabel('Time since injection (years)')
  ylabel('Fraction of CO2 remaining sequestered')
  title('Sequestration fraction over time')
  set(gca,'ylim',[0 1])
end
