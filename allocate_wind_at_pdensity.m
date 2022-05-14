% ALLOCATE_WIND_AT_PDENSITY
%
% Chooses the pixels with the highest time averaged cf for wind power
% allocation and returns a cellarray with the indices of the chosen pixels
% as well as the total allocated power
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:
% X:                                                  latitude meshgrid
% Y:                                                  longitude meshgrid
% datamat_orig:                                       2D array wind data
% index_cluster_cellarray:                            cellarray of 1-column vectors per cluster; contains the datamat-indices of all points within the cluster
% pd:                                                 2D array: map speed->power
% pixelarea:                                          2D array area [km²] of all pixels
% capdens_MAT:                                        2D array capacity density
% Cmax:                                               total capacity to be allocated to that region
%
% OUTPUT:
% index_country_cellarray_allocated:                  cellarray of 1D index arrays per cluster, defining those pixels where wind power has been allocated
% Calloc:                                             2D array size(X) containing the capacity density for each pixel (constant capdens if allocated according to highest wind speeds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index_cellarray_allocated,Calloc]=allocate_wind_at_pdensity(X,Y,cf_time_avg_mat,index_cellarray,pd,pixelarea,capdens_mat,Cmax,frac_alloc_min)

%% Test
%  index_cellarray={index_cluster_cellarray{9}}
%  Cmax=Cmax_cluster_cells(8)

%Create vector with all wind speeds
%[windspeed,countryindex]
vector_pindex_cn=[];
tic;for i=1:length(index_cellarray)
    if not(isempty(index_cellarray{i}))
        vector_pindex_cn=[vector_pindex_cn;[index_cellarray{i},i*ones(length(index_cellarray{i}),1)]];
    end
end;toc

%add the power output times the pixel area as a first column for all indices in the previously
%defined array and sort the resulting array according to cf_time_avg_mat

%                                        |          1st column: total energy produced during one year given capdens_mat Wh/yr    | 2nd: Indices   |  3rd avg CF each point of the cluster |                                                                
vector_powercum_pindex_cn_sorted=flipdim(sortrows([pixelarea(vector_pindex_cn(:,1)).*capdens_mat(vector_pindex_cn(:,1))*3600*8760,vector_pindex_cn,cf_time_avg_mat(vector_pindex_cn(:,1))],4),1);

% using the matlab function cumsum create a vector containing the sum of the
% power output times the area of the corresponding pixel 
% [1:length(vector_powercum_pindex_cn_sorted)].'
% cumsum(vector_powercum_pindex_cn_sorted(:,1))
% vector_powercum_pindex_cn_sorted(:,2:end)
vector_powercum_pindex_cn_sorted=[[1:length(vector_powercum_pindex_cn_sorted)].',cumsum(vector_powercum_pindex_cn_sorted(:,1)),vector_powercum_pindex_cn_sorted(:,2:end)]

%cut the subset up to a maximum cumulative power of Cmax OR to include at
%least a fraction frac_alloc_min of the total pixels
n_cmax=sum(vector_powercum_pindex_cn_sorted(:,2)<=Cmax)
n_fam=floor(length(vector_powercum_pindex_cn_sorted(:,2))*frac_alloc_min)
vector_powercum_pindex_cn_sorted_csum=vector_powercum_pindex_cn_sorted(1:max(n_cmax,n_fam),:)

%save the actually allocated capacity
Calloc=0;
if length(vector_powercum_pindex_cn_sorted_csum)>0
    Calloc=vector_powercum_pindex_cn_sorted_csum(end,2)
end

%create a cellarray to sort the thus selected pixels by country
index_cellarray_allocated={};
tic;for i=1:length(index_cellarray)
	sum(vector_powercum_pindex_cn_sorted_csum(:,4)==i);
    index_cellarray_allocated{i}=vector_powercum_pindex_cn_sorted_csum(vector_powercum_pindex_cn_sorted_csum(:,4)==i,3);
end

index_cellarray_allocated
index_cellarray