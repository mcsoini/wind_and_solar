%Reads the border dataset and returns relevant data for the subsequent
%calculations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ~~~INPUT:
% filename:                         directory and filename of the border data file
% X:                                latitude meshgrid
% Y:                                longitude meshgrid
% countryCellArray:                 cellarray of country name strings, consistent with the names in the border data file
% exclude_countries_cellArray:      cellarray of country name strings, to be excluded in any case
%
% ~~~OUTPUT:
% xbreg:                            1D array x-values of the region's border points
% ybreg:                            1D array y-values of the region's border points
% xbglob:                           1D array x-values of all countries' border points
% ybglob:                           1D array y-values of all countries' border points
% index_array_offshore_global:      array of datamat-indices -> all offshore points, i.e. outside of country borders
% coord_country_cellarray_global:   cellarray of 2-column vectors per country GLObAL; contains x and y values of the border points
% index_country_cellarray_global:   cellarray of 1-column vectors per country GLObAL; contains the datamat-indices of all points within the borders of a country      
% coord_country_cellarray_region:   cellarray of 2-column vectors per country REGION; contains x and y values of the border points              
% index_country_cellarray_region:   cellarray of 1-column vectors per country REGION; contains the datamat-indices of all points within the borders of a country      
% map_countries_lists:              2D array -> maps the indices of the countries in the countrycellarray to the indices of the same countries in the country border data set
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xbreg,ybreg,xbglob,ybglob,index_array_offshore_global,coord_country_cellarray_global,index_country_cellarray_global,coord_country_cellarray_region,index_country_cellarray_region,map_countries_lists]=readcountries(filename,X,Y,countryCellArray,exclude_countries_cellArray)
'... entered readcountries ...'


%Load the file containing the data of the country borders
%countryData = shaperead('/mnt/data/wind/border/ne_110m_admin_0_countries_lakes_mod.dbf','UseGeoCoords',true);
countryData = shaperead(filename,'UseGeoCoords',true);

%Initialize the aggregate border data arrays to accomodate the border data of all countries in
%the input country cell array
xbreg=[];
ybreg=[];
xbglob=[];
ybglob=[];
coord_country_cellarray_global={};
map_countries_lists=[];
%Loop over all countries of the country data set
for i=1:length(countryData)
    coord_country_cellarray_global{i}=[countryData(i).Lon.',countryData(i).Lat.'];
    countries_in_dataset{i}=countryData(i).name;
    strcmp(countryData(i).name,countryCellArray{:});
    
    %Check whether country i is in the input countryCellArray
    position_country_in_input_cellarray=strcmp(countryData(i).name,countryCellArray{:})
    
    %Check whether country i is in the input exclude_countries_cellArray
    position_country_in_exclude_countries_cellarray=strcmp(countryData(i).name,exclude_countries_cellArray)

    %save all border coordinates in two arrays, useful for quick plotting
    xbglob =[xbglob,countryData(i).Lon];
    ybglob =[ybglob,countryData(i).Lat];
    
    if sum(position_country_in_exclude_countries_cellarray)==0 && sum(position_country_in_input_cellarray)~=0
        %store an array which maps the index of the countries in the input
        %country list to the index of the same country in the dataset list
        map_countries_lists=[map_countries_lists;find(position_country_in_input_cellarray==1),i]
       
        %store the borders of all countries in countryCellArray individually in a
        %cellarray borderdata
        [~,index_country_in_countryCellArray]=max(strcmp(countryData(i).name,countryCellArray{:}));
        coord_country_cellarray_region{index_country_in_countryCellArray}=[countryData(i).Lon.',countryData(i).Lat.'];

        %add the border of country i to the aggregate border data arrays
        xbreg =[xbreg,countryData(i).Lon];
        ybreg =[ybreg,countryData(i).Lat];
    end
end
map_countries_lists=sortrows(map_countries_lists,1)

%Store the country names as found in the database in a separate text file;
%this is useful to have a look at how the countries are actually called
fid = fopen('countries_in_database', 'wt');
fprintf(fid, '"%s"\n', countries_in_dataset{:});
fclose(fid);
% h=figure(1)
% pcolorjw(x,y,datamat)
% hold on
% plot(xv,yv,'w','LineWidth',1.5)

%Create cellarray country_indices with one cell per country of the input region containing
%all indices of the pixels within that country
tic; index_country_cellarray_region={}; %all indices within each country (region)
for i=1:length(coord_country_cellarray_region)
    if length(coord_country_cellarray_region{i})~=0
        index_country_cellarray_region{i}=find(inpoly([X(:),Y(:)].',[coord_country_cellarray_region{i}(:,1).';coord_country_cellarray_region{i}(:,2).']));
    end
end; toc

%Create cellarray country_indices_glob with one cell per country containing
%all indices of the pixels within that country
tic; index_country_cellarray_global={}; %all indices within each country (global)
parfor i=1:length(coord_country_cellarray_global)
    if length(coord_country_cellarray_global{i})~=0
        index_country_cellarray_global{i}=find(inpoly([X(:),Y(:)].',[coord_country_cellarray_global{i}(:,1).';coord_country_cellarray_global{i}(:,2).'])).';
    end
end; toc

%create an array index_array_offshore_global with all the offshore indices
LV_all_indices_datamat=1:numel(X);
index_array_offshore_global=setdiff(LV_all_indices_datamat,vertcat(index_country_cellarray_global{:})).';

'... leaving readcountries ...'