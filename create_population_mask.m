% CREATE_POPULATION_MASK

% Determine the pixels suitable for allocation based on population data and
% thresholds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
% population_datafile:         directory and filename of the population data set
% X:                           meshgrid matrix of the latitudes
% Y:                           meshgrid matrix of the longitudes
% xbglob:                      global boundary x coordinates (just for plotting)
% ybglob:                      global boundary y coordinates (just for plotting)
% index_array_offshore_global: 1D array of all offshore indices (outside all country boundaries)
% index_array_onshore_global:  1D array of all onshore indices (within all country boundaries)
% pop_thresh_upp:              upper threshold of the poplation density [1/km²]
% pop_thresh_low:              lower threshold of the population density [1/km²]
% dist_thresh_offsh:           max great circle distance around populated pixels to be included (offshore) [km]
% dist_thresh_onsh:            max great circle distance around populated pixels to be included (onshore) [km]
% area_opening_threshold:      number of pixels of small populated areas to be excluded
% area_opening_connectivity:   connectivity of small populated areas to be excluded (4 or 8... 8 less stringent)
% 
% OUTPUT:
% pop_mask_out:                2D binary matrix; population mask indicating potential pixels 
% pop_mask_out_step:           2D matrix; different values for offshore pixels, onshore pixels, above upper population threshold pixels (just for illustration)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% REMARKS:
% -> the matlab function "distance" might actually be faster than
%    "distgreatcircle"! (but probably not "distgreatcircle_approx")
% -> Do some square pre-filtering before calculating the distances between
%    onshore and offshore
% -> Use of "imresize" to scale data is not optimal!! However, no change 
%    introduced if the resolution is the same from the beginning
% -> "distgreatcircle_approx" gives only a few wrong pixels for distance
%    thresholds of 500km onshore and offshore; maybe necessary to use
%    full accuracy "distgreatcircle" for larger distances
function [pop_mask_out,pop_mask_out_step,pop_mat_exp_res_output]=create_population_mask(population_datafile,X,Y,xbglob,ybglob,index_array_offshore_global,index_array_onshore_global,pop_thresh_upp,pop_thresh_low,dist_thresh_offsh,dist_thresh_onsh,area_opening_threshold,area_opening_connectivity,resource)
%%
% population_datafile=strcat(filedir{7},filedir{8});
% X;
% Y;
% xbglob;
% ybglob;
% index_array_offshore_global=index_array_offshore_global;
% index_array_onshore_global=vertcat(index_country_cellarray_global{:});
% pop_thresh_upp=pop_thresh_upp;
% pop_thresh_low=pop_thresh_low;
% dist_thresh_offsh=dist_thresh_offsh;
% dist_thresh_onsh=dist_thresh_onsh;
% area_opening_threshold=area_opening_threshold;
% area_opening_connectivity=area_opening_connectivity;
% resource=resource;

global output_dir_graph

'... entered create_population_mask ...'

%define a maximum number of pixel combinations to be calculated by the mex
%function distgreatcircle during each call. The optimal value probably depends on
%the computer characteristics
maxComb=1e7;

%import and extract data
population_dataset=importdata(population_datafile,' ',6);
pop_mat=population_dataset.data;

%The population data does not cover the whole globe;
%Specify vectors popLat and popLon manually, based on the information in the
%population_dataset.textdata information
dLon=0.5;
Lon0=-180;
LonN=720; %Number data point longitudinal direction
dLat=0.5;
Lat0=-58;
LatN=286;

%increase the geographical range of the population data to cover the whole
%globe (zero for all pixels not covered by the dataset)
popLon_exp=Lon0:dLon:180-0.0001;
popLat_exp=90:-dLat:-90;
pop_mat_exp=[zeros(sum([popLat_exp>Lat0+(LatN-1)*dLat]),LonN);pop_mat;zeros(sum([popLat_exp<Lat0]),LonN)];

%Output population density matrix for use outside of the function
pop_mat_exp_res_output=imresize(pop_mat_exp, size(X));
pop_mat_exp_res_output(pop_mat_exp_res_output<0)=0;

%create a filter mask which is zero except for where the population density
%is above a certain threshold...
pop_mat_mask=pop_mat_exp>=pop_thresh_low;
%...and another one for the upper population threshold
pop_mat_mask_upp=pop_mat_exp<=pop_thresh_upp;

%resize the population mask to fit the size of the wind data matrix
%this might be kind of dubious for low resolution data
pop_mat_mask_res=imresize(pop_mat_mask, size(X));
pop_mat_mask_upp_res=imresize(pop_mat_mask_upp, size(X));

%remove small spots of population from the population mask:
%---all spots with area_opening_threshold pixels or less
%---"area_opening_connectivity"-connected (8 is less stringent)
pop_mat_mask_res_oparea=bwareaopen(pop_mat_mask_res,area_opening_threshold,area_opening_connectivity);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    OFFSHORE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(resource,'wind')

    %Get the coordinates of all offshore pixels
    P1=X(index_array_offshore_global);
    L1=Y(index_array_offshore_global);

    %Get the coordinates of all populated pixels
    P2=X(pop_mat_mask_res_oparea==1).';
    L2=Y(pop_mat_mask_res_oparea==1).';

    %Get the parameters to slice the coordinate vectors in order to contain a
    %maximum of maxComb elements, defined above; this is necessary to limit
    %the number of pixels which are passed to the distance-calculating function
    %in one go... otherwise it might just freeze
    nlocInp=floor(maxComb/length(P2));
    nlocSlabs=floor(length(P1)/nlocInp)

    tic
    A={};
    %Loop over the coordinate slices
    'Calculating the distances between all populated pixels and all offshore pixels...'
    parfor_progress(nlocSlabs+1);
    parfor i=1:nlocSlabs+1
        [length(P2),length(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i)))),length(P2)*length(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))))];
        [nlocInp*(i-1)+1,min(length(P1),nlocInp*(i))];
        %Call the function to calculate the distances between all combinations
        %of n offshore pixels and m populated pixels -> returns a n*m vector
        A{i}=distgreatcircle_approx(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))).',L1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))).',P2,L2);
        %Reshape the n*m vector to obtain a nxm matrix
        A{i}=reshape(A{i},[length(P2),length(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))))]);
        %For each offshore pixel, only keep the distance to the closest
        %populated pixel
        A{i}=min(A{i},[],1);
        parfor_progress;
    end
    parfor_progress(0);
    toc

    %Combine all slices to obtain a vector with the length equal to the number
    %of offshore pixels and values equal to the distance to the closest
    %populated pixel
    mindist_offsh=[A{:}].';
end
%%%%%%%%%%    OFFSHORE       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    ONSHORE        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get the coordinates of all onshore pixels
P1=X(index_array_onshore_global);
L1=Y(index_array_onshore_global);

%Get the coordinates of all populated pixels
P2=X(pop_mat_mask_res_oparea==1).';
L2=Y(pop_mat_mask_res_oparea==1).';

%Get the parameters to slice the coordinate vectors in order to contain a
%maximum of maxComb elements, defined above
nlocInp=floor(maxComb/length(P2));
nlocSlabs=floor(length(P1)/nlocInp);

tic
A={};
'Calculating the distances between all populated pixels and all onshore pixels...'
parfor_progress(nlocSlabs+1);
parfor i=1:nlocSlabs+1
    [length(P2),length(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i)))),length(P2)*length(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))))];
    [nlocInp*(i-1)+1,min(length(P1),nlocInp*(i))];
    %Call the function to calculate the distances between all combinations
    %of n onshore pixels and m populated pixels -> returns a n*m vector
    A{i}=distgreatcircle_approx(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))).',L1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))).',P2,L2);
    %Reshape the n*m vector to obtain a nxm matrix
    A{i}=reshape(A{i},[length(P2),length(P1(nlocInp*(i-1)+1:min(length(P1),nlocInp*(i))))]);
    %For each onshore pixel, only keep the distance to the closest
    %populated pixel
    A{i}=min(A{i},[],1);
    parfor_progress;
end
'...done.'
parfor_progress(0);
toc

%Combine all slices to obtain a vector with the length equal to the number
%of offshore pixels and values equal to the distance to the closest
%populated pixel
mindist_onsh=[A{:}].';

%%%%%%%%%%    ONSHORE        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%create the output population mask with different values for populated
%pixels, offshore pixels, onshore pixels, and above upper threshold pixels;
%only used for visualization
pop_mask_out_step=zeros(size(pop_mat_mask_res_oparea));
if strcmp(resource,'wind')
    pop_mask_out_step(index_array_offshore_global(mindist_offsh<dist_thresh_offsh))=2;
end
pop_mask_out_step(index_array_onshore_global(mindist_onsh<dist_thresh_onsh))=3;
pop_mask_out_step(pop_mat_mask_res_oparea==1)=1;
pop_mask_out_step(pop_mat_mask_upp_res==0)=4;

%create the output mask; same as pop_mat_mask_res_new but binary; this is
%the actual population mask
pop_mask_out=zeros(size(pop_mask_out_step));
pop_mask_out(pop_mask_out_step==2)=1;
pop_mask_out(pop_mask_out_step==3)=1;
pop_mask_out(pop_mask_out_step==4)=2;
pop_mask_out(pop_mask_out_step==1)=1;

%Plot the population masks (binary and step)
nfig=1414;
f=figure(nfig)
%screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );

subplot(1,2,1)
plot(xbglob,ybglob,'w')
hold on
pcolorjw(X,Y,pop_mask_out_step)
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')

subplot(1,2,2)
plot(xbglob,ybglob,'w')
hold on
pcolorjw(X,Y,pop_mask_out)
xlim([-180,180])
ylim([-90,90])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')

%print(f,'-dpng',strcat(output_dir_graph, 'AA_population_mask_',resource,'.png'))


'... leaving create_population_mask ...'
