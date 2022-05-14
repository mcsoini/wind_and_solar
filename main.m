%%%List of used packages/toolboxes
%nctoolbox from https://code.google.com/p/nctoolbox/
%   appears to require setup for every matlab launch
%image processing toolbox
%inpoly mex function (replaces matlab's inpolygon (much faster))
%   http://www.mathworks.com/matlabcentral/fileexchange/15410-inpoly-mex-file
%distgreatcircle mex function

%DATA SETS:
%WIND SPEEDS:
%http://apps.ecmwf.int/datasets/data/interim_full_daily/?levtype=pl
%POPULATION DENSITY
%http://sedac.ciesin.columbia.edu/data/set/gpw-v3-population-count-future-estimates/data-download
%BATHYMETRY
%http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NGDC/.ETOPO5/datasetdatafiles.html

clear all

addpath('../ext/sun_position_functions/SolarAzEl/')
addpath('../ext/sun_position_functions/sun_position/')
addpath('../ext/tightfig/')
addpath('../ext/parfor_progress/')
addpath('../ext/suplabel')
addpath('../ext/csv_write_headers/')
addpath('../ext/')

%resource='solar';
resource='wind';

set(0,'defaulttextfontsize',12);
set(0,'defaultaxesfontsize',12);
addpath('inpoly')
save_data=cell(5,1)

global plotfilename
global plotlimits
global output_dir_data
global output_dir_graph

output_dir_data='../output_data/';
output_dir_graph='../output_figures/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   DEFINITION OF FILE NAMES AND DIRECTORIES   %%%%%%%%%%%%%%%%%%%%%

%Define the paths and filenames for all datasets
%CHALMERS LINUX SYSTEM
% filedir{1}='/chalmers/users/soini/wind/data/wind/';         %Directory of the wind data file
% filedir{2}='interim_UV_57of60_0p5deg_6h_2013.grib';      %Name of the wind data file
% filedir{3}='/chalmers/users/soini/wind/data/speed2power/';   %Directory of the wind turbine output characteristics file
% filedir{4}='mclean_future_lowland';                         %Name of the wind turbine output characteristics data
% filedir{5}='/chalmers/users/soini/wind/data/borders/';      %Directory border data set
% filedir{6}='ne_50m_admin_0_countries_lakes_mod.dbf';        %Filename border data set
% filedir{7}='/chalmers/users/soini/wind/data/population/';   %Directory population data set
% filedir{8}='glds00ag30.asc';                                %Filename population data set
% filedir{9}='/chalmers/users/soini/wind/data/bathymetric/';
% filedir{10}='bathymetry_data.cdf';

%Define the paths and filenames for all datasets
%TOURS SYSTEM
% filedir{1}='/mnt/data/wind/data/';                      %Directory of the wind data file
% filedir{2}='interim_UV_57of60_1p5deg_6h_2013_12.grib';  %Name of the wind data file
% filedir{3}='/mnt/data/wind/wind_2_power/';              %Directory of the wind turbine output characteristics file
% filedir{4}='mclean_future_lowland';                     %Name of the wind turbine output characteristics data
% filedir{5}='/mnt/data/wind/border/';                    %Directory border data set
% filedir{6}='ne_110m_admin_0_countries_lakes_mod.dbf';   %Filename border data set
% filedir{7}='/mnt/data/wind/population/';                %Directory population data set
% filedir{8}='glds00ag30.asc';                            %Filename population data set
% filedir{9}='/mnt/data/wind/bathymetry/';
% filedir{10}='bathymetry_data.cdf';

%Define the paths and filenames for all datasets
%LINUX RELATIVE
filedir{1}='../data/wind/';                              %Directory of the wind data file
filedir{2}='interim_UV_57of60_0p5deg_6h_2013.grib';      %Name of the wind data file
filedir{3}='../data/speed2power/';                       %Directory of the wind turbine output characteristics file
filedir{4}='mclean_future_lowland';                      %Name of the wind turbine output characteristics data
filedir{5}='../data/borders/';                           %Directory border data set
filedir{6}='ne_50m_admin_0_countries_lakes_mod.dbf';     %Filename border data set
filedir{7}='../data/population/';                        %Directory population data set
filedir{8}='glds00ag30.asc';                             %Filename population data set
filedir{9}='../data/bathymetric/';                       %Directory bathymetry data set
filedir{10}='bathymetry_data.cdf';                       %Filename bathymetry data set
filedir{11}='../data/solar/';                            %Directory of the solar data file
filedir{12}='interim_ssrd_0p5deg_3_2013_1_7.grib';       %Name of the solar data file

%%%%%%%%   DEFINITION OF FILE NAMES AND DIRECTORIES   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%              DATA EXTRACTION           %%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(resource, 'wind')
    %Wind:
    [x,y,datamat0,rawdata,h_vec,pressure]=getdata4Dm(strcat(filedir{1},filedir{2}),'U_velocity_hybrid','V_velocity_hybrid',1,'windabsolute', 'avg',0);
elseif strcmp(resource, 'solar')
    %Solar:
    [x,y,datamatrix_pv0,rawdata_pv,datamatrix_csp0,rawdata_csp,h_vec]=getdata4Dm_solar(strcat(filedir{11},filedir{12}));
end

h_vec=double(h_vec);

[X,Y]=meshgrid(x,y);
pixelarea=111*111*abs(x(end-1)-x(end))*abs(y(end-1)-y(end))*cos(Y*pi/180);

%temporal resolution of the data:
dt=double(h_vec(2)-h_vec(1)); %[hours]

if strcmp(resource, 'wind')
    %import the data mapping wind speeds to capacity factor
    pd=importdata(strcat(filedir{3},filedir{4}));
end
        
%define array of colors for the subregions
subcol = {'b',[.8 .2 .6],'g',[.5 .6 .7],'r'};

%%%%%%%%%%              DATA EXTRACTION           %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    DATA STRUCTURING AND FILTERING      %%%%%%%%%%%%%%%%%%%%%%%%%
%Copy matrices -> NECESSARY?
if strcmp(resource,'wind')
    datamat_orig=datamat0;
end

close all;

%%
%%%%%%%%%%%%%%%%     DEFINITION OF THE REGIONS       %%%%%%%%%%%%%%%%%%%%%%

%Exclude some countries for political etc reasons
exclude_countries={'Somalia', 'Central African Rep.','Somaliland','Albania'}

%Get the parameters specific to each region; possible choices for
%region_str:
% 'FSU'  -> Former Soviet Union
% 'CPA'  -> Centrally planned Asia and China
% 'PAS'  -> Other Pacific Asia
% 'LAC'  -> Latin America and the Caribbean
% 'MEA'  -> Middle East and North Africa
% 'NAM'  -> North America
% 'AFR'  -> Sub-Saharan Africa
% 'SAS'  -> South Asia
% 'PAO'  -> Pacific OECD
% 'TEU'  -> Central and Eastern Europe combined with Western Europe

region_str='TEU'
%All relevant information on the regions stored in the function get_region_parameters
[countryCellArray,plotfilename0,excludeCoords,plotlimits,initial_cluster_centroids_onsh,initial_cluster_centroids_offsh]=get_region_parameters(region_str)

countryCellArray={countryCellArray} %For compatibility with the following functions

%%%%%%%%%%%%%%%%   ACCESSING THE BORDER DATA SET     %%%%%%%%%%%%%%%%%%%%%%

% Filter the data according to the location, i.e. set it zero everyone
% outside the region of consideration. This needs an cellarray of strings
% with country names as an input argument
'Calling function readcountries'
[xbreg,ybreg,xbglob,ybglob,index_array_offshore_global, coord_country_cellarray_global,index_country_cellarray_global,coord_country_cellarray_region,index_country_cellarray_region,map_countries_lists]=readcountries(strcat(filedir{5},filedir{6}),X,Y,countryCellArray,exclude_countries);
'...readcountries done.'

%%%%%%%%%%    DATA STRUCTURING AND FILTERING      %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%   FILTERING BY POPULATION DENSITY   %%%%%%%%%%%%%%%%%%%%%%

if strcmp(resource, 'wind')
    tic
    pop_thresh_low=10;
    pop_thresh_upp=500;
    dist_thresh_offsh=500; %[km]
    dist_thresh_onsh=500; %[km]
    area_opening_threshold=3;
    area_opening_connectivity=8;
    % THE FUNCTION readcountries NEEDS TO BE CALLED FIRST ONCE (SEE BELOW)
    [pop_mask00,pop_mask_step,popdens_mat]=create_population_mask(strcat(filedir{7},filedir{8}),X,Y,xbglob,ybglob,index_array_offshore_global,vertcat(index_country_cellarray_global{:}),pop_thresh_upp,pop_thresh_low,dist_thresh_offsh,dist_thresh_onsh,area_opening_threshold,area_opening_connectivity,resource);
elseif strcmp(resource, 'solar')
    tic
    pop_thresh_low=10;
    pop_thresh_upp=500; %Lower?
    dist_thresh_offsh=0; %[km]
    dist_thresh_onsh=500; %[km]
    area_opening_threshold=3;
    area_opening_connectivity=8;
    % THE FUNCTION readcountries NEEDS TO BE CALLED FIRST ONCE (SEE BELOW)
    [pop_mask00,pop_mask_step,popdens_mat]=create_population_mask(strcat(filedir{7},filedir{8}),X,Y,xbglob,ybglob,index_array_offshore_global,vertcat(index_country_cellarray_global{:}),pop_thresh_upp,pop_thresh_low,dist_thresh_offsh,dist_thresh_onsh,area_opening_threshold,area_opening_connectivity,resource);
end
%pop_mask_step code:
%4: above population threshold
%1: populated
%3: within the onshore range


%%
%Copy binary population mask to keep the original
pop_mask0=pop_mask00;

%Define a population mask with priorities 1 and 2 to be able to allocate PV
%in populated places with a higher priority
pop_mask_priorit=zeros(size(pop_mask_step));
pop_mask_priorit(pop_mask_step==3 | pop_mask_step==1)=1;
pop_mask_priorit(pop_mask_step==4)=2;

%exclude certain regions defined above in the array excludeCoords, e.g. to
%cut the European overseas territories, esp. French Guiana
%TODO: modify the border data set to make French Guiana a part of SAmerica
%(probably negligible)
pop_mask_priorit(X>excludeCoords(1))=0;
pop_mask_priorit(X<excludeCoords(2))=0;
pop_mask_priorit(Y>excludeCoords(3))=0;
pop_mask_priorit(Y<excludeCoords(4))=0;
pop_mask0(X>excludeCoords(1))=0;
pop_mask0(X<excludeCoords(2))=0;
pop_mask0(Y>excludeCoords(3))=0;
pop_mask0(Y<excludeCoords(4))=0;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    ONSHORE        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
if strcmp(resource,'wind')
    plotfilename=[plotfilename0 'wind_onsh_']
elseif strcmp(resource,'solar')
    plotfilename=[plotfilename0 'solar_onsh_']
end

%------     GET RELEVANT ONSHORE PIXELS     -------------------------------

%value_thresh: threshold for the average data values (e.g. wind speeds) for filtering prior to
%the definition of the clusters (increases the optimality of the allocation)
if strcmp(resource,'wind')
    value_thresh=5;
elseif strcmp(resource,'solar')
    value_thresh=0;
end
%apply the population mask to the onshore pixels; returns a cellarray of
%indices for each country

if strcmp(resource,'wind')
    index_country_cellarray_region_onsh_pop_filtered=get_onshore_pixels_filtered(X,Y,xbreg,ybreg,datamat_orig,pop_mask0==1,index_country_cellarray_region,value_thresh);
elseif strcmp(resource,'solar')
    %choice of the datamatrix in the arguments is not relevant, as long as
    %value_thresh=0
    index_country_cellarray_region_onsh_pop_filtered=get_onshore_pixels_filtered(X,Y,xbreg,ybreg,datamatrix_pv0,pop_mask_priorit,index_country_cellarray_region,value_thresh);
end

%----     CALCULATE ENERGY OUTPUT FOR ALL RELEVANT ONSHORE PIXELS     -----

%Calculate the time average of the cf of all relevant pixels, stored in the
%matrix cf_time_avg_mat (for pixel ranking)
if strcmp(resource,'wind') %need to use the function speed2power since the relation power(speed) is non-linear
    ['Time averaging the cf of ' num2str(length([index_country_cellarray_region_onsh_pop_filtered{:}].')) ' pixels']
    tic; cf_time_avg_onsh_mat=zeros(size(datamat0));
    cf_time_avg_onsh_mat([index_country_cellarray_region_onsh_pop_filtered{:}].')=permute(mean(speed2power(rawdata(:,[index_country_cellarray_region_onsh_pop_filtered{:}].'),pd)), [ 2 3 1 ])/100;
    toc
end

close all

%%
%%%%%%%%%%%%%%%   CLUSTERS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Reload region parameters here (finding good centroids)
[countryCellArray,plotfilename0,excludeCoords,plotlimits,initial_cluster_centroids_onsh,initial_cluster_centroids_offsh]=get_region_parameters(region_str)

%----------     DEFINE CLUSTERS     ---------------------------------------

%create clusters of the onshore pixels
clustering_type='kmeans'
cluster_control=initial_cluster_centroids_onsh;
nf=1 %number of cluster pairs furthest apart (OBSOLETE)
[index_cluster_cellarray_onsh,~,final_cluster_centroids]=find_clusters(X,Y,xbreg,ybreg,{[index_country_cellarray_region_onsh_pop_filtered{:}].'},cluster_control,clustering_type,nf)

%----------     ALLOCATE CAPACITY TO CLUSTERS     -------------------------

%minimum fraction of pixels within a cluster to be allocated
frac_alloc_min=0.0
if strcmp(resource,'wind')
    %capacity density (only relevant for allocation according to wind speed ranking and low frac_alloc_min!)
    capdens_onsh=0.5*1e6;%W/km^2
    %maximum capacity to be allocated in the region
    Cmax_onsh=20*10^18; %[EJ]
    poppropctrl=0.3;
    alloc_strat_clusters='cl_powpopprop' %'cl_equal' or 'cl_areaprop', 'cl_popprop', 'cl_powpopprop'
    alloc_strat_pixels='px_windrank' %'px_windrank' or 'px_areaprop'
    [index_cluster_cellarray_onsh_allocated,capdens_onsh_mat,pop_cluster,Calloc_clusters,cluster_area]=allocate_capacity_clusters(alloc_strat_clusters,alloc_strat_pixels,X,Y,xbreg,ybreg,cf_time_avg_onsh_mat,index_cluster_cellarray_onsh,pd,pixelarea,Cmax_onsh,capdens_onsh,popdens_mat,index_country_cellarray_region_onsh_pop_filtered,final_cluster_centroids,value_thresh,frac_alloc_min,poppropctrl);
elseif strcmp(resource,'solar')
    %percentage of land covered with panels
    radiation_to_ac_eff_pv= 0.2*0.75; %rough guess 20% Radiation to dc efficiency and 75% dc to ac efficiency
    radiation_to_ac_eff_csp=0.3;
    fraction_pv=0.7; %fraction of pixels PV allocated outside of areas above the population threshold
    areadens_pv=0.025;
    areadens_csp=0.025;
    Cmax_total_pv =15; %[EJ] maximum pv energy per year each region
    Cmax_total_csp=5; %[EJ] maximum csp energy per year each region
    alloc_strat_clusters='cl_powpopprop' %'cl_equal' or 'cl_areaprop', 'cl_popprop', 'cl_powpopprop'
    alloc_strat_pixels='px_solarrank' %'px_solarrank' or 'px_areaprop'
    [index_cluster_cellarray_onsh_allocated,pop_cluster,Calloc_clusters_pv,Calloc_clusters_csp,cluster_area]=allocate_capacity_clusters_solar(alloc_strat_clusters, alloc_strat_pixels, X, Y, xbreg, ybreg, pop_mask_priorit, pixelarea, popdens_mat, datamatrix_pv0, datamatrix_csp0, index_cluster_cellarray_onsh, index_country_cellarray_region_onsh_pop_filtered, Cmax_total_pv, Cmax_total_csp, areadens_pv, areadens_csp, fraction_pv, final_cluster_centroids, value_thresh, frac_alloc_min, radiation_to_ac_eff_pv, radiation_to_ac_eff_csp);
end

%get a list of all possible pairs of cluster indices
cpairs = nchoosek(1:length(index_cluster_cellarray_onsh_allocated), 2)

%%
%%%%%%%
%SOLAR%
%%%%%%%
if strcmp(resource, 'solar')
    cf_clusterpair_cellarray_pv={};
    cf_clusterpair_cellarray_csp={};
    E_tot_pair_csp=zeros(size(length(cpairs),1));
    E_tot_pair_pv=zeros(size(length(cpairs),1));

    'Calculating total cf for the whole region'

    ira=vertcat(vertcat(index_cluster_cellarray_onsh_allocated{:})); %Indices Region Allocated
    [P_t_pv_region,P_t_csp_region,E_tot_pv_region,E_tot_csp_region]=cum_cf_arb_index_list_solar(ira(ira(:,2)==1 | ira(:,2)==2,1),ira(ira(:,2)==3,1),rawdata_pv,rawdata_csp,pixelarea,radiation_to_ac_eff_pv, radiation_to_ac_eff_csp,areadens_pv, areadens_csp,h_vec);
    for i=1:length(cpairs)
        %indices pv: all pixels with code "1-PV Urban" or "2-PV"
        indices_pv_cluster_1_of_pair =[index_cluster_cellarray_onsh_allocated{cpairs(i,1)}(index_cluster_cellarray_onsh_allocated{cpairs(i,1)}(:,2)==1 | index_cluster_cellarray_onsh_allocated{cpairs(i,1)}(:,2)==2,1)];
        indices_pv_cluster_2_of_pair =[index_cluster_cellarray_onsh_allocated{cpairs(i,2)}(index_cluster_cellarray_onsh_allocated{cpairs(i,2)}(:,2)==1 | index_cluster_cellarray_onsh_allocated{cpairs(i,2)}(:,2)==2,1)];
        %indices csp: all pixels with code "3-CSP"
        indices_csp_cluster_1_of_pair=[index_cluster_cellarray_onsh_allocated{cpairs(i,1)}(index_cluster_cellarray_onsh_allocated{cpairs(i,1)}(:,2)==3,1)];
        indices_csp_cluster_2_of_pair=[index_cluster_cellarray_onsh_allocated{cpairs(i,2)}(index_cluster_cellarray_onsh_allocated{cpairs(i,2)}(:,2)==3,1)];

        %'both clusters' -> cumulative power output
        [cf_clusterpair_cellarray_pv{i}(:,3),cf_clusterpair_cellarray_csp{i}(:,3),E_tot_pair_pv(i),E_tot_pair_csp(i)]=cum_cf_arb_index_list_solar([indices_pv_cluster_1_of_pair; indices_pv_cluster_2_of_pair],[indices_csp_cluster_1_of_pair; indices_csp_cluster_2_of_pair],rawdata_pv,rawdata_csp,pixelarea,radiation_to_ac_eff_pv, radiation_to_ac_eff_csp,areadens_pv, areadens_csp,h_vec)
        %'cluster 1'
        [cf_clusterpair_cellarray_pv{i}(:,1),cf_clusterpair_cellarray_csp{i}(:,1),~,~]=cum_cf_arb_index_list_solar(indices_pv_cluster_1_of_pair,indices_csp_cluster_1_of_pair,rawdata_pv,rawdata_csp,pixelarea,radiation_to_ac_eff_pv, radiation_to_ac_eff_csp,areadens_pv, areadens_csp,h_vec)    
        %'cluster 2'
        [cf_clusterpair_cellarray_pv{i}(:,2),cf_clusterpair_cellarray_csp{i}(:,2),~,~]=cum_cf_arb_index_list_solar(indices_pv_cluster_2_of_pair,indices_csp_cluster_2_of_pair,rawdata_pv,rawdata_csp,pixelarea,radiation_to_ac_eff_pv, radiation_to_ac_eff_csp,areadens_pv, areadens_csp,h_vec)    
        %'all clusters'
        cf_clusterpair_cellarray_pv{i}(:,4) =P_t_pv_region ;
        cf_clusterpair_cellarray_csp{i}(:,4)=P_t_csp_region;
    end
    % %%
    nselect=3;
    %add a column to indicate the pair number and sort by increasing power
    cum_power_vecsorted_pv =sortrows([[1:length(cf_clusterpair_cellarray_pv)].', E_tot_pair_pv.' ],2);
    cum_power_vecsorted_csp=sortrows([[1:length(cf_clusterpair_cellarray_csp)].',E_tot_pair_csp.'],2);

    %the (n_total-nselect) pairs with the lowest output
    cf_pairs_min_pv =cum_power_vecsorted_pv(1:end-nselect,1);
    cf_pairs_min_csp=cum_power_vecsorted_csp(1:end-nselect,1);

    %the nselect pairs with the highest output
    cf_pairs_max_pv=cum_power_vecsorted_pv(end-nselect+1:end,1);
    cf_pairs_max_csp=cum_power_vecsorted_csp(end-nselect+1:end,1);

    %number of pairs sorted according to output
    cf_pairs_all_pv=cum_power_vecsorted_pv(:,1);
    cf_pairs_all_csp=cum_power_vecsorted_csp(:,1);


    %%%%%%%%%%%%%%%%
    %PLOTTING SOLAR%
    %%%%%%%%%%%%%%%%

    % %% Plot the duration curves for the cluster pairs with the highest and the lowest energy output for the whole duration
    nfig=3462356;
    f=figure(nfig)
    screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );

    legend_cellarray_min_pv ={};
    legend_cellarray_max_pv ={};
    legend_cellarray_min_csp={};
    legend_cellarray_max_csp={};
    %Plot all cumulative cf curves except for those of the nselect pairs with
    %the highest cumulative output
    if length(cf_pairs_min_pv)>0
        for i=1:length(cf_pairs_min_pv)
            subplot(2,2,1)
            hold all
            h=plot(h_vec,sort(cf_clusterpair_cellarray_pv{cf_pairs_min_pv(i)}(:,3)),'Linewidth',2,'color',[0 0 0]+0.7);
            legend_cellarray_min_pv={legend_cellarray_min_pv{:}, ['pair (' num2str(cpairs(cf_pairs_min_pv(i),:)) '); Tot. energy produced=' num2str(E_tot_pair_pv(cf_pairs_min_pv(i)),'%.3f') 'EJ; Total capacity=' num2str(E_tot_pair_pv(cf_pairs_min_pv(i))/1e9,'%.3f') 'GW; Avg cf=' num2str(mean(cf_clusterpair_cellarray_pv{cf_pairs_min_pv(i)}(:,3)),'%.3f') ]}
            hAnnotation = get(h,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')

            subplot(2,2,2)
            hold all
            h=plot(h_vec,sort(cf_clusterpair_cellarray_csp{cf_pairs_min_csp(i)}(:,3)),'Linewidth',2,'color',[0 0 0]+0.7);
            legend_cellarray_min_csp={legend_cellarray_min_csp{:}, ['pair (' num2str(cpairs(cf_pairs_min_csp(i),:)) '); Tot. energy produced=' num2str(E_tot_pair_csp(cf_pairs_min_csp(i)),'%.3f') 'EJ; Total capacity=' num2str(E_tot_pair_csp(cf_pairs_min_csp(i))/1e9,'%.3f') 'GW; Avg cf=' num2str(mean(cf_clusterpair_cellarray_csp{cf_pairs_min_csp(i)}(:,3)),'%.3f') ]}
            hAnnotation = get(h,'Annotation');
            hLegendEntry = get(hAnnotation','LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off')     

            subplot(2,2,3)
            hold all
            h=plot(h_vec,(cf_clusterpair_cellarray_pv{cf_pairs_min_pv(i)}(:,3)),'Linewidth',2,'color',[0 0 0]+0.7);

            subplot(2,2,4)
            hold all
            h=plot(h_vec,(cf_clusterpair_cellarray_csp{cf_pairs_min_csp(i)}(:,3)),'Linewidth',2,'color',[0 0 0]+0.7);
        end
    end
    %plot the cf curves of the nselect pairs with the highest cumulative output
    for i=1:length(cf_pairs_max_pv)
        subplot(2,2,1)
        hold all
        plot(h_vec,sort(cf_clusterpair_cellarray_pv{cf_pairs_max_pv(i)}(:,3)),'-','Linewidth',3);
        legend_cellarray_max_pv={legend_cellarray_max_pv{:}, [ ...
            'pair (' num2str(cpairs(cf_pairs_max_pv(i),:)) '); ', ...
            'Tot.E=' num2str(E_tot_pair_pv(cf_pairs_max_pv(i)),'%.3f') 'EJ; ', ...
            'Avg cf=' num2str(mean(cf_clusterpair_cellarray_pv{cf_pairs_max_pv(i)}(:,3)),'%.3f') ]}

        subplot(2,2,2)
        hold all
        plot(h_vec,sort(cf_clusterpair_cellarray_csp{cf_pairs_max_csp(i)}(:,3)),'-','Linewidth',3);
        legend_cellarray_max_csp={legend_cellarray_max_csp{:}, [ ...
            'pair (' num2str(cpairs(cf_pairs_max_csp(i),:)) '); ', ...
            'Tot.E=' num2str(E_tot_pair_csp(cf_pairs_max_csp(i)),'%.3f') 'EJ; ', ...
            'Avg cf=' num2str(mean(cf_clusterpair_cellarray_csp{cf_pairs_max_csp(i)}(:,3)),'%.3f') ]}

        subplot(2,2,3)
        hold all
        plot(h_vec,(cf_clusterpair_cellarray_pv{cf_pairs_max_pv(i)}(:,3)),'-','Linewidth',2);

        subplot(2,2,4)
        hold all
        plot(h_vec,(cf_clusterpair_cellarray_csp{cf_pairs_max_csp(i)}(:,3)),'-','Linewidth',2);
    end

    legend_cellarray_pv={legend_cellarray_max_pv{:}};
    legend_cellarray_csp={legend_cellarray_max_csp{:}};

    subplot(2,2,1)
    legend(legend_cellarray_pv,'Location','NorthWest')
    ylim([0,inf]); xlim([0,8760]); xlabel('Duration/time [h]'); ylabel('cf')
    title('PV')
    subplot(2,2,2)
    legend(legend_cellarray_csp,'Location','NorthWest')
    ylim([0,inf]); xlim([0,8760]); xlabel('Duration/time [h]'); ylabel('cf')
    title('CSP')
    subplot(2,2,3)
    ylim([0,inf]); xlim([2200,2300]); xlabel('Duration/time [h]'); ylabel('cf')
    subplot(2,2,4)
    ylim([0,inf]); xlim([2200,2300]); xlabel('Duration/time [h]'); ylabel('cf')

    title_string=[plotfilename, '; ', ...
        'Cmaxpv=',num2str(Cmax_total_pv),'EJ; ', ...
        'Cmaxcsp=',num2str(Cmax_total_csp),'EJ; ', ...
        'Area dens.pv: ', num2str(areadens_pv), ' ', ...
        'Area dens.csp: ', num2str(areadens_csp), ' ', ...
        'PV priority: ', num2str(fraction_pv), ...
        alloc_strat_clusters,'; ', ...
        alloc_strat_pixels];
    [ax,h3]=suplabel(title_string ,'t');
    set(h3,'FontSize',15) 

    drawnow
    tightfig
    screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );

    %print(f,'-dpng',strcat(output_dir_graph, plotfilename,'_cf_curves.png'))

    %Write data to file PV
    headers_pv={'time [h]','cf region'}
    data_pv=[h_vec, P_t_pv_region];
    headers_csp={'time [h]','cf region'}
    data_csp=[h_vec, P_t_csp_region];
    for i=1:length(cpairs(:,1))
        headers_pv={headers_pv{:}, ['cf pair(' num2str(cpairs(i,1)),' ', num2str(cpairs(i,2)) ')']};
        data_pv=[data_pv, cf_clusterpair_cellarray_pv{i}(:,3)];
        headers_csp={headers_csp{:}, ['cf pair(' num2str(cpairs(i,1)),' ', num2str(cpairs(i,2)) ')']};
        data_csp=[data_csp, cf_clusterpair_cellarray_csp{i}(:,3)];
    end

    csvwrite_with_headers(strcat(output_dir_data, plotfilename,'PV_cf_data.csv'),data_pv,headers_pv)
    csvwrite_with_headers(strcat(output_dir_data, plotfilename,'CSP_cf_data.csv'),data_csp,headers_csp)

elseif strcmp(resource, 'wind')
 
    %%%%%%
    %WIND%
    %%%%%%
    %-------     CALCULATE THE CUMULATIVE OUTPUT OF THE CLUSTERS     ----------

    %for each of the cluster pairs get the aggregate cf
    cf_clusterpair_cellarray={}
    'Calculating total cf for the whole region'
    [tot_cf,~]=cum_cf_arb_index_list(vertcat(index_cluster_cellarray_onsh_allocated{:}),rawdata,pixelarea,capdens_onsh_mat,h_vec,pd);
    E_tot_pair=zeros(size(length(cpairs),1))
    for i=1:length(cpairs)
    %PARFOR USED IN SPEED2POWER-> don't use it here, since it copies
    %rawmat
        [i,length(cpairs)]
        %'both clusters' -> cumulative cf total energy output
        [cf_clusterpair_cellarray{i}(:,3),E_tot_pair(i),cap_clusterpair_final(i)]=cum_cf_arb_index_list(vertcat(index_cluster_cellarray_onsh_allocated{cpairs(i,1)},index_cluster_cellarray_onsh_allocated{cpairs(i,2)}),rawdata,pixelarea,capdens_onsh_mat,h_vec,pd);
        %'cluster 1'
        [cf_clusterpair_cellarray{i}(:,1),~]=cum_cf_arb_index_list(index_cluster_cellarray_onsh_allocated{cpairs(i,1)},rawdata,pixelarea,capdens_onsh_mat,h_vec,pd);
        %'cluster 2'
        [cf_clusterpair_cellarray{i}(:,2),~]=cum_cf_arb_index_list(index_cluster_cellarray_onsh_allocated{cpairs(i,2)},rawdata,pixelarea,capdens_onsh_mat,h_vec,pd);
        %'all clusters'
        cf_clusterpair_cellarray{i}(:,4)=tot_cf;
    end
    % %%
    nselect=3;
    %add a column to indicate the pair number and sort by increasing power
    cum_power_vecsorted=sortrows([[1:length(cf_clusterpair_cellarray)].', E_tot_pair.'],2)

    %the nselect pairs with the lowest output
    cf_pairs_min=cum_power_vecsorted(1:end-nselect,1)

    %the nselect pairs with the highest output
    cf_pairs_max=cum_power_vecsorted(end-nselect+1:end,1)

    %number of pairs sorted according to output
    cf_pairs_all=cum_power_vecsorted(:,1)

    %%%%%%%%%%%%%%%
    %PLOTTING WIND%
    %%%%%%%%%%%%%%%
    % %% Plot the duration curves for the cluster pairs with the highest and the lowest energy output for the whole duration
    nfig=3462346;
    f=figure(nfig)
    screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );
    hold all

    legend_cellarray_min={};
    legend_cellarray_max={};
    %Plot all cumulative cf curves except for those of the nselect pairs with
    %the highest cumulative output
    for i=1:length(cf_pairs_min)
        h=plot(h_vec,sort(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3)),'Linewidth',2,'color',[0 0 0]+0.7);
        legend_cellarray_min={legend_cellarray_min{:}, ['pair (' num2str(cpairs(cf_pairs_min(i),:)) '); Tot. energy produced=' num2str(E_tot_pair(cf_pairs_min(i)),'%.3f') 'EJ; Total capacity=' num2str(cap_clusterpair_final(cf_pairs_min(i))/1e9,'%.3f') 'GW; Avg cf=' num2str(mean(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3)),'%.3f') '; Median cf=' num2str(median(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3)),'%.3f')]}
        hAnnotation = get(h,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
    end
    %plot the cf curves of the nselect pairs with the highest cumulative output
    for i=1:length(cf_pairs_max)
        plot(h_vec,sort(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3)),'-','Linewidth',3);
        legend_cellarray_max={legend_cellarray_max{:}, ['pair (' num2str(cpairs(cf_pairs_max(i),:)) '); Tot. energy produced=' num2str(E_tot_pair(cf_pairs_max(i)),'%.3f') 'EJ; Total capacity=' num2str(cap_clusterpair_final(cf_pairs_max(i))/1e9,'%.3f') 'GW; Avg cf=' num2str(mean(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3)),'%.3f') '; Median cf=' num2str(median(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3)),'%.3f')]}
    end

    legend_cellarray={legend_cellarray_max{:}}
    legend(legend_cellarray,'Location','NorthWest')
    xlabel('Duration/time [h]'); ylabel('cf'),ylim([0,max(pd(:,2)/100)])
    title(strcat(plotfilename, '; Cmax=',num2str(Cmax_onsh/1e18),'EJ; ','Orig. cap. density=',num2str(capdens_onsh/1e3),'kW/km2 ; Avg wind speed threshold=',num2str(value_thresh),'m/s;',alloc_strat_clusters,'; ',alloc_strat_pixels,'; Minimum fraction pixels allocated=',num2str(frac_alloc_min)),'interpreter','none')
    drawnow
    tightfig
    screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );

    %print(f,'-dpng',strcat(output_dir_graph, plotfilename,'cf_curves.png'))
    
    %Write data to file
    headers_wind={'time [h]','cf region'}
    data_wind=[h_vec, tot_cf];
    for i=1:length(cpairs(:,1))
        headers_wind={headers_wind{:}, ['cf pair(' num2str(cpairs(i,1)),' ', num2str(cpairs(i,2)) ')']};
        data_wind=[data_wind, cf_clusterpair_cellarray{i}(:,3)];
    end

    csvwrite_with_headers(strcat(output_dir_data, plotfilename,'cf_data.csv'),data_wind,headers_wind)

    % %% Store just plotted data in variable for later use (if needed, e.g. for comparing parameters)
    % 
    % j=5
    % for i=1:nselect
    %     save_data{j}=[save_data{j} sort(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3))]
    % end
    % for i=1:nselect
    %     save_data{j}=[save_data{j} sort(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3))]
    % end
end

%%%%%%%%%%    ONSHORE        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    OFFSHORE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
plotfilename=[plotfilename0 'wind_offsh_']

%speed_thresh: threshold for the average wind speeds for filtering prior to
%the definition of the clusters (increases the optimality of the allocation)
value_thresh=5

%returns a cellarray with one cell per country (global, same order as in
%the country dataset) containing the indices of the offshore locations within
%the range of depth specified by depth_threshold and closest to a
%POPULATED pixel of the corresponding country

%maximum depth of the seafloor to accomodate wind power
depth_threshold=-50 %[m]
bathy_datafile=strcat(filedir{9},filedir{10});
[index_country_cellarray_global_offshore_pop_filtered]=get_offshore_pixels_filtered_global(bathy_datafile,X,Y,datamat_orig,depth_threshold,value_thresh,coord_country_cellarray_global,pop_mask0==1,index_array_offshore_global)

% %%
%extract those countries which are in the countryCellArray
%index_country_cellarray_region_offsh contains all the pixels for potential
%offshore capacity allocation sorted by country in the countryCellArray
index_country_cellarray_region_offsh={}
for i=1:length(countryCellArray{:})
    if sum(map_countries_lists(:,1)==i)~=0
        map_countries_lists(map_countries_lists(:,1)==i,2)
        index_country_cellarray_region_offsh{i}=index_country_cellarray_global_offshore_pop_filtered{map_countries_lists(map_countries_lists(:,1)==i,2)}
    else
        index_country_cellarray_region_offsh{i}=[]
    end
end

%----     CALCULATE ENERGY OUTPUT FOR ALL RELEVANT OFFHORE PIXELS     -----

%Calculate the time average of the cf of all relevant pixels, stored in the
%matrix cf_time_avg_mat
['Time averaging the cf of ' num2str(length(vertcat(index_country_cellarray_region_offsh{:}))) ' pixels']
tic; cf_time_avg_offsh_mat=zeros(size(datamat0));
cf_time_avg_offsh_mat(vertcat(index_country_cellarray_region_offsh{:}))=permute(mean(speed2power(rawdata(:,vertcat(index_country_cellarray_region_offsh{:})),pd)), [ 2 3 1 ])/100;
toc
close all

%%%%%%%%%%%%%%%   CLUSTERS       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------     DEFINE CLUSTERS     ---------------------------------------

%create clusters of the onshore pixels
clustering_type='kmeans'
cluster_control=initial_cluster_centroids_offsh;
nf=1 %number of cluster pairs furthest apart (OBSOLETE)
[index_cluster_cellarray_offsh,~,final_cluster_centroids]=find_clusters(X,Y,xbreg,ybreg,{vertcat(index_country_cellarray_region_offsh{:})},cluster_control,clustering_type,nf)

%----------     ALLOCATE CAPACITY TO CLUSTERS     -------------------------

%minimum fraction of pixels within a cluster to be allocated
frac_alloc_min=0.0
%capacity density (only relevant for allocation according to wind speed ranking and low frac_alloc_min!)
capdens_offsh=0.5*1e6;%W/km^2
%maximum capacity to be allocated in the region
Cmax_offsh=5*10^18; %[EJ]
poppropctrl=0.3
alloc_strat_clusters='cl_powpopprop' %'cl_equal' or 'cl_areaprop', 'cl_popprop'
alloc_strat_pixels='px_windrank' %'px_windrank' or 'px_areaprop'
[index_cluster_cellarray_offsh_allocated,capdens_offsh_mat,pop_cluster,Calloc_clusters,cluster_area]=allocate_capacity_clusters(alloc_strat_clusters,alloc_strat_pixels,X,Y,xbreg,ybreg,cf_time_avg_offsh_mat,index_cluster_cellarray_offsh,pd,pixelarea,Cmax_offsh,capdens_offsh,popdens_mat,index_country_cellarray_region,final_cluster_centroids,value_thresh,frac_alloc_min,poppropctrl);
%get a list of all possible pairs of cluster indices
cpairs = nchoosek(1:length(index_cluster_cellarray_offsh_allocated), 2)

% %%
%-------     CALCULATE THE CUMULATIVE OUTPUT OF THE CLUSTERS     ----------

%for each of the cluster pairs get the aggregate cf
cf_clusterpair_cellarray={}
'Calculating total cf for the whole region'
[tot_cf,~]=cum_cf_arb_index_list(vertcat(index_cluster_cellarray_offsh_allocated{:}),rawdata,pixelarea,capdens_offsh_mat,h_vec,pd);
E_tot_pair=zeros(size(length(cpairs),1))
for i=1:length(cpairs)
%!!!!!PARFOR USED IN SPEED2POWER!!!!!-> don't use it here, since it copies
%rawmat
    [i,length(cpairs)]
    %'both clusters' -> cumulative cf total energy output
    [cf_clusterpair_cellarray{i}(:,3),E_tot_pair(i),cap_clusterpair_final(i)]=cum_cf_arb_index_list(vertcat(index_cluster_cellarray_offsh_allocated{cpairs(i,1)},index_cluster_cellarray_offsh_allocated{cpairs(i,2)}),rawdata,pixelarea,capdens_offsh_mat,h_vec,pd);
    %'cluster 1'
    [cf_clusterpair_cellarray{i}(:,1),~]=cum_cf_arb_index_list(index_cluster_cellarray_offsh_allocated{cpairs(i,1)},rawdata,pixelarea,capdens_offsh_mat,h_vec,pd);
    %'cluster 2'
    [cf_clusterpair_cellarray{i}(:,2),~]=cum_cf_arb_index_list(index_cluster_cellarray_offsh_allocated{cpairs(i,2)},rawdata,pixelarea,capdens_offsh_mat,h_vec,pd);
    %'all clusters'
    cf_clusterpair_cellarray{i}(:,4)=tot_cf;
end

% %%
nselect=3;
%add a column to indicate the pair number and sort by increasing power
cum_power_vecsorted=sortrows([[1:length(cf_clusterpair_cellarray)].', E_tot_pair.'],2)

%the nselect pairs with the lowest output
cf_pairs_min=cum_power_vecsorted(1:end-nselect,1)

%the nselect pairs with the highest output
cf_pairs_max=cum_power_vecsorted(end-nselect+1:end,1)

%number of pairs sorted according to output
cf_pairs_all=cum_power_vecsorted(:,1)

% %% Plot the duration curves for the cluster pairs with the highest and the lowest energy output for the whole duration
nfig=3462346;
f=figure(nfig)
screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );
hold all

legend_cellarray_min={};
legend_cellarray_max={};
%Plot all cumulative cf curves except for those of the nselect pairs with
%the highest cumulative output
for i=1:length(cf_pairs_min)
    h=plot(h_vec,sort(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3)),'Linewidth',2,'color',[0 0 0]+0.7);
    legend_cellarray_min={legend_cellarray_min{:}, ['pair (' num2str(cpairs(cf_pairs_min(i),:)) '); Tot. energy produced=' num2str(E_tot_pair(cf_pairs_min(i)),'%.3f') 'EJ; Total capacity=' num2str(cap_clusterpair_final(cf_pairs_min(i))/1e9,'%.3f') 'GW; Avg cf=' num2str(mean(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3)),'%.3f') '; Median cf=' num2str(median(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3)),'%.3f')]}
    hAnnotation = get(h,'Annotation');
    hLegendEntry = get(hAnnotation','LegendInformation');
    set(hLegendEntry,'IconDisplayStyle','off')
end
%plot the cf curves of the nselect pairs with the highest cumulative output
for i=1:length(cf_pairs_max)
    plot(h_vec,sort(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3)),'-','Linewidth',3);
    legend_cellarray_max={legend_cellarray_max{:}, ['pair (' num2str(cpairs(cf_pairs_max(i),:)) '); Tot. energy produced=' num2str(E_tot_pair(cf_pairs_max(i)),'%.3f') 'EJ; Total capacity=' num2str(cap_clusterpair_final(cf_pairs_max(i))/1e9,'%.3f') 'GW; Avg cf=' num2str(mean(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3)),'%.3f') '; Median cf=' num2str(median(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3)),'%.3f')]}
end

legend_cellarray={legend_cellarray_max{:}}
legend(legend_cellarray,'Location','NorthWest')
xlabel('Duration/time [h]'); ylabel('cf'),ylim([0,max(pd(:,2)/100)])
title(strcat(plotfilename,'; Cmax=',num2str(Cmax_offsh/1e18),'EJ; ','Orig. cap. density=',num2str(capdens_offsh/1e3),'kW/km2 ; Avg wind speed threshold=',num2str(value_thresh),'m/s;',alloc_strat_clusters,'; ',alloc_strat_pixels,'; Minimum fraction pixels allocated=',num2str(frac_alloc_min)),'interpreter','none')
drawnow
tightfig
screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );

%print(f,'-dpng',strcat(output_dir_graph, plotfilename,'cf_curves.png'))

%Write data to file
headers_wind={'time [h]','cf region'}
data_wind=[h_vec, tot_cf];
for i=1:length(cpairs(:,1))
    headers_wind={headers_wind{:}, ['cf pair(' num2str(cpairs(i,1)),' ', num2str(cpairs(i,2)) ')']};
    data_wind=[data_wind, cf_clusterpair_cellarray{i}(:,3)];
end

csvwrite_with_headers(strcat(output_dir_data, plotfilename,'cf_data.csv'),data_wind,headers_wind)

%% Store just plotted data in variable for later use (if needed, e.g. for comparing parameters)

j=5
for i=1:nselect
    save_data{j}=[save_data{j} sort(cf_clusterpair_cellarray{cf_pairs_max(i)}(:,3))]
end
for i=1:nselect
    save_data{j}=[save_data{j} sort(cf_clusterpair_cellarray{cf_pairs_min(i)}(:,3))]
end

%%%%%%%%%%    OFFSHORE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
