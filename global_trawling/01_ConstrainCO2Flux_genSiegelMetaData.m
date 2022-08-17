% 01_ConstrainCO2Flux_genSiegelMetaData.m
% Created June 7, 2022 by Jamie Collins, jcollins@edf.org, based on the script 
% plot_sequestration_fraction.m provided by Siegel et al. at
%  https://doi.org/10.6084/m9.figshare.15228690.v2

% this script collects and writes some metadata from the model output
% provided by Siegel et al., which we'll need in the R script
% '02_ConstrainCO2Flux_IO.R' and beyond

% Note: Assumes user has loaded "fseq_OCIM2_48L.mat" into memory 

% variables are as follows
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

% find the bottom depth
bottom_depth = sum(VOL.*MASK,3)./AREA(:,:,1); % depth of the ocean at each water column

% save matrices containing the land-sea mask (just surface depth), model
% time steps, lat and long,  list of depths in the model domain, and bottom
% depth for each grid cell

writematrix(MASK(:,:,1),'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/mask_surface.csv')
writematrix(time,'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/model_timesteps.csv')
writematrix(DEPTH(1,1,:),'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/OCIM_modelDepths.csv')
writematrix(LAT(:,1,1),'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/lat_degN.csv')
writematrix(LON(1,:,1),'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/long_degE.csv')
writematrix(bottom_depth,'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/bottom_depth_m.csv')

% % save .csv file containing all the seq fractions
% 
% writematrix(fseq,'/Users/jamesrco/Code/global-trawling-CO2/data/global_trawling/derived/benthic_seqfractions/fseq_all.csv')
