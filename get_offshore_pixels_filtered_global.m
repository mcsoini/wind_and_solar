% GET_OFFSHORE_PIXELS_FILTERED_GLOBAL
%
% This assigns all offshore pixels to the closest country after filtering
% according to the depth threshold;
% then applies the population mask to get those pixels which are close
% enough to a populated area

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:
% bathy_datafile:                                       directory and filename of the bathymetric data
% X:                                                    latitude meshgrid
% Y:                                                    longitude meshgrid
% depth_threshold:                                      maximum sea floor depth for wind power allocation
% coord_country_cellarray_global:                       cellarray of the border coordinates of all countries (global)
% pop_mask:                                             binary population mask
% index_array_offshore_global:                          1D array containing all indices of offshore pixels
% 
% OUTPUT:
% index_country_cellarray_global_offshore_pop_filtered: cellarray of 1D arrays of datamat-indices containing the offshore pixels suited for wind power allocation
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [index_country_cellarray_global_offshore_pop_filtered]=get_offshore_pixels_filtered_global(bathy_datafile,X,Y,datamat,depth_threshold,speed_thresh,coord_country_cellarray_global,pop_mask,index_array_offshore_global)

% %% Test
% bathy_datafile
% X
% Y
% depth_threshold
% coord_country_cellarray_global
% pop_mask
% index_array_offshore_global


%load the bathymetric data
bathydataset=ncgeodataset(bathy_datafile);

xdataset=bathydataset.geovariable('X');
zdataset=bathydataset.geovariable('elev');

xoffsh=xdataset.data(:);
elev=zdataset.data(:);

clear xdataset zdataset

%shift data by 180° west to cover the same range as the other datasets
valshift=-180;
[~,indshift]=min(abs(xoffsh+valshift))
elev=[elev(:,indshift:end),elev(:,1:indshift-1)];

%filter according to sea floor depth
elevmask=elev<0 & elev > depth_threshold;

%interpolating the elevation mask to fit the wind speeds matrix
elevmask_res=imresize(elevmask, size(X));

%extract the indices where the resized elevation mask is not zero
offshore_ind=find(elevmask_res);

%limit the index list to those pixels which are outside country borders
%(i.e. in the index_array_offshore_global list)
offshore_ind=offshore_ind(ismember(offshore_ind,index_array_offshore_global));

%%%%%%%%   ASSIGN COUNTRIES TO OFFSHORE PIXELS   %%%%%%%%%%%%%%%%%%%%%%%%%%

% measure distance to pixels inside the countries
% crude approach using just the coordinates (in deg) to calculate a
% distance; probably ok, since only used for ranking, and short distances
% on the surface, and not too close to the poles/extreme longitudes

distance_pixels_countries=zeros(length(X(offshore_ind)),length(coord_country_cellarray_global));
tic
%get x and y coordinates of the offshore pixels
LV_offshore_x=X(offshore_ind);
LV_offshore_y=Y(offshore_ind);
%Loop over the global country list
parfor_progress(length(coord_country_cellarray_global));
'Calculate distances offshore pixels to countries'
parfor i=1:length(coord_country_cellarray_global)    
    %Extract x and y coordinates of the pixels of the current country (not really necessary)
    LV_country_x=coord_country_cellarray_global{i}(:,1);
    LV_country_y=coord_country_cellarray_global{i}(:,2);
    
    %Create a meshgrid to easily calculate the distances between all
    %combinations of offshore pixels and country pixels
    [LV_X_country,LV_X_offshore]=meshgrid(LV_country_x,LV_offshore_x);
    [LV_Y_country,LV_Y_offshore]=meshgrid(LV_country_y,LV_offshore_y);
    
    %Calculate the "distance" as the sum of the squares of the x and y
    %distances of country pixels and offshore pixels
    distance_pixels_country_i=(LV_X_country-LV_X_offshore).^2+(LV_Y_country-LV_Y_offshore).^2;
    
    %
    distance_pixels_countries(:,i)=min(distance_pixels_country_i,[],2)
    parfor_progress;
end
parfor_progress(0);
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%make list with length equal the number of offshore pixels and entries the
%index of the countries they are closest to
[~,pixels_closest_country_indices]=min(distance_pixels_countries,[],2);

%make a cellarray offshore_pixels_sorted_by_country with one cell per
%country containing the list of indices corresponding to the offshore
%pixels closest to that country
offshore_pixels_sorted_by_country={};
bglobvec=vertcat(coord_country_cellarray_global{:});
for i=1:length(coord_country_cellarray_global)
    offshore_pixels_sorted_by_country{i}=offshore_ind(pixels_closest_country_indices==i);
end

%%%%%%%%%%%%%%%   APPLY THE POPULATION FILTER MASK TO THE DATA %%%%%%%%%%%%
pcolorjw(X,Y,1-pop_mask)
hold on

ind_pop_mask=find(pop_mask);

for i=1:length(offshore_pixels_sorted_by_country)
    index_country_cellarray_global_offshore_pop_filtered{i}=offshore_pixels_sorted_by_country{i}(ismember(offshore_pixels_sorted_by_country{i},ind_pop_mask));
    
    index_country_cellarray_global_offshore_pop_filtered{i}=index_country_cellarray_global_offshore_pop_filtered{i}(datamat(index_country_cellarray_global_offshore_pop_filtered{i})>speed_thresh);
    
    colormap(gray)

    plot(X(index_country_cellarray_global_offshore_pop_filtered{i}),Y(index_country_cellarray_global_offshore_pop_filtered{i}),'r.')
    plot(X(offshore_pixels_sorted_by_country{i}),Y(offshore_pixels_sorted_by_country{i}),'og')
end
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')

