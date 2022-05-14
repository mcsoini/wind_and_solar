% GET_ONSHORE_PIXELS_FILTERED
%
% This looks at the pixels corresponding to the countries in the chosen
% region and checks which of them can accomodate wind power according to
% the population mask; plus: filtering according to a lower wind speed
% threshold

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:
% X:                                                  latitude meshgrid
% Y:                                                  longitude meshgrid
% xbreg:                                              1D array x-values of the region's border points
% ybreg:                                              1D array y-values of the region's border points
% datamat:                                            2D wind data array
% pop_mask:                                           2D binary matrix; population mask indicating potential pixels
% index_country_cellarray_region:                     cellarray of 1-column vectors per country REGION; contains the datamat-indices of all points within the borders of a country      
% speed_thresh:                                       wind speed threshold for pre-filtering; probably too optimizing
% 
% OUTPUT:
% index_country_cellarray_region_onsh_pop_filtered:   cellarray of 1D arrays (per country) of datamat-indices containing the onshore pixels suited for wind power allocation
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [index_country_cellarray_region_onsh_pop_filtered]=get_onshore_pixels_filtered(X,Y,xbreg,ybreg,datamat,pop_mask,index_country_cellarray_region,speed_thresh)
'... entered get_onshore_pixels_filtered ...'

% global plotlimits %declare the global variable plotlimits

%%%%%%%%%%%%%%%   APPLY THE POPULATION FILTER MASK TO THE DATA %%%%%%%%%%%%
index_country_cellarray_region_onsh_pop_filtered={}

%plot the binary population mask as a background
% pcolorjw(X,Y,1-pop_mask)
% hold all
% plot(xbreg,ybreg,'w')

%get the indices of the pixels where the binary population mask is equal 1
ind_pop_mask=find(pop_mask);

%Loop over the countries of the regions
for i=1:length(index_country_cellarray_region)    
    %for country i, get those indices for which the binary population mask
    %is non zero
    index_country_cellarray_region_onsh_pop_filtered{i}=index_country_cellarray_region{i}(ismember(index_country_cellarray_region{i},ind_pop_mask));
    %for country i, get the indices for which the wind speed is above the
    %threshold
    index_country_cellarray_region_onsh_pop_filtered{i}=index_country_cellarray_region_onsh_pop_filtered{i}(datamat(index_country_cellarray_region_onsh_pop_filtered{i})>speed_thresh);

    %plot the result
    colormap(gray)
    plot(X(index_country_cellarray_region_onsh_pop_filtered{i}),Y(index_country_cellarray_region_onsh_pop_filtered{i}),'.')
    plot(X(index_country_cellarray_region{i}),Y(index_country_cellarray_region{i}),'og')
    hold on
end
%some more plotting options
% ylim([plotlimits(2,1),plotlimits(2,2)])
% xlim([plotlimits(1,1),plotlimits(1,2)])
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% drawnow

'... leaving get_onshore_pixels_filtered ...'
