
%DATA: ETOPO5
%http://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NGDC/.ETOPO5/datasetdatafiles.html

function [offshore_pixels_sorted_by_country_pop_filtered]=get_offshore_pixels_filtered(X,Y,datamat,depth_threshold,bglob,country_indices_glob)

bathydataset=ncgeodataset('data.cdf');

xdataset=bathydataset.geovariable('X');
zdataset=bathydataset.geovariable('elev');

xoffsh=xdataset.data(:);
elev=zdataset.data(:);

clear xdataset zdataset

%shift data by 180° west
valshift=-180;
[~,indshift]=min(abs(xoffsh+valshift))
elev=[elev(:,indshift:end),elev(:,1:indshift-1)];

%filter
elevmask=elev<0 & elev > depth_threshold;

%%
%interpolating the elevation mask to fit to the wind speeds matrix
elevmask_res=imresize(elevmask, size(datamat));

%extract the indices where the resized elevation mask is not zero
offshore_ind=find(elevmask_res);

%get the (real world) coordinates of the offshore pixels
offshore_coords=[X(offshore_ind),Y(offshore_ind)];

%calculate the distance between the offshore pixels and all countries
%parallel computing really makes a difference here
%what about -180->180... probably not relevant, since in the middle of the
%pacific ocean
distance_pixels_countries=zeros(length(offshore_coords),length(bglob));
tic
for i=1:length(bglob)
    LV_bglobi=bglob{i};
    strcat('Offshore pixels: Calculating distances to country no.',int2str(i),'/',int2str(length(bglob)))
    parfor j=1:length(offshore_coords)
        distance_pixels_countries(j,i)=p_poly_dist(offshore_coords(j,1),offshore_coords(j,2),LV_bglobi(:,1),LV_bglobi(:,2));
    end
end
toc

%remove the values of pixels that are inside country borders from both the
%distance matrix and the list of offshore pixels coordinates
[ind_pixels_inside,~]=find(distance_pixels_countries<0);
distance_pixels_countries_outside=distance_pixels_countries(setdiff(1:length(distance_pixels_countries),ind_pixels_inside),:);
offshore_ind_outside=offshore_ind(setdiff(1:length(distance_pixels_countries),ind_pixels_inside),:);

%make list with length equal the number of offshore pixels and entries the
%index of the countries they are closest to
[~,pixels_closest_country_indices]=min(distance_pixels_countries_outside,[],2);

%%
%make a cellarray offshore_pixels_sorted_by_country with one cell per
%country containing the list of indices corresponding to the offshore
%pixels closest to that country
offshore_pixels_sorted_by_country={};
bglobvec=vertcat(bglob{:});
for i=1:length(bglob)
    offshore_pixels_sorted_by_country{i}=offshore_ind_outside(pixels_closest_country_indices==i);
end

tic
for i=1:length(offshore_pixels_sorted_by_country)
    if not(isempty(bglob{i})) && not(isempty(offshore_pixels_sorted_by_country{i})) && not(isempty(country_indices_glob{i}))
        %get the coordinates of the offshore pixels associated with country i
        LV_country_offshore_coords=[X(offshore_pixels_sorted_by_country{i}),Y(offshore_pixels_sorted_by_country{i})];

        %get the index of the countries pixel closest to each offshore pixel
        ind_near=knnsearch([X(country_indices_glob{i}),Y(country_indices_glob{i})],LV_country_offshore_coords);

        %use the pop_mat_mask_res array to create a vector with entries 1 if 
        LV_country_pop_mat_mask_res=pop_mat_mask_res(country_indices_glob{i});

        %get the one-cell-per-country cellarray containing the indices of
        %the offshore pixels close to a poulated pixel
        offshore_pixels_sorted_by_country_pop_filtered{i}=offshore_pixels_sorted_by_country{i}(LV_country_pop_mat_mask_res(ind_near));

%         plot(X(offshore_pixels_sorted_by_country_pop_filtered{i}),Y(offshore_pixels_sorted_by_country_pop_filtered{i}),'go','linewidth',2)
%         hold on
%         colormap(gray)
%         plot(X(offshore_pixels_sorted_by_country{i}),Y(offshore_pixels_sorted_by_country{i}),'xr','linewidth',2)
%         pcolorjw(X,Y,pop_mat_mask_res)
%         plot(bglob{i}(:,1),bglob{i}(:,2),'linewidth',2)
%         xlabel('Longitude [deg]')
%         ylabel('Latgitude [deg]')
%         title(countryData(i).name)
%         drawnow
%         hold off
        
    end
end
toc


