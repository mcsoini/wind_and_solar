% FIND_CLUSTERS
%
% finds clusters according to the
% -kmeans
% -hierarchical (to be implemented)
% methods and determines the pairs of
% clusters which are furthest apart

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:
% X:                                                  latitude meshgrid
% Y:                                                  longitude meshgrid
% xbreg:                                              1D array x-values of the region's border points
% ybreg:                                              1D array y-values of the region's border points
% index_country_cellarray:                            cellarray of 1-column vectors per country; contains the datamat-indices of all points within the borders of a country
% cluster_control:                                    2xn array specifying the initial centroids of the clusters or single number n specifying the number of clusters (random initial centroids)
% clustering_type:                                    string defining the type of clustering to be used; options: 'kmeans'
% nf:                                                 number nf of pairs of clusters whos centroids are apart the furthest
%
% OUTPUT:
% index_cluster_cellarray:                            cellarray of 1D index arrays for each cluster
% cdist_f:                                            cell indices of the
% nf pairs of clusters whose centroids are furthest apart (2x10 array) -> obsolete
% final_cluster_centroids:                            2xn array specifying the final centroids of the clusters
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [index_cluster_cellarray,cdist_f,final_cluster_centroids]=find_clusters(X,Y,xbreg,ybreg,index_country_cellarray,cluster_control, clustering_type,nf)

global plotlimits


% plot(xbreg,ybreg)
% hold all
% plot(X(index_country_cellarray{:}),Y(index_country_cellarray{:}),'o')
% plot(cluster_control(:,1),cluster_control(:,2),'xk','Linewidth',5)
% hold off
% ylim([plotlimits(2,1),plotlimits(2,2)])
% xlim([plotlimits(1,1),plotlimits(1,2)])
% xlabel('Longitude [deg]')
% ylabel('Latitude [deg]')
% drawnow
% pause


% Dump the pixel indices of all countries in a single array (countries are irrelevant for the formation of the clusters)
LV_ind_list=vertcat(index_country_cellarray{:});
% Get an array with all coordinates of the relevant pixels
LV_xy=[X(LV_ind_list),Y(LV_ind_list)];

%opts = statset('Display','final');



if numel(cluster_control)<2
    nc=cluster_control
    if strcmp(clustering_type,'kmeans')
        [CL,final_cluster_centroids]=kmeans(LV_xy,nc,'Replicates',1);%,'Options',opts);
    %Hierarchical clustering -> not really useful
    elseif strcmp(clustering_type,'hierarchy')
        CL=clusterdata(LV_xy,nc);
        final_cluster_centroids=zeros(max(CL),2);
        for i=1:max(CL)
            final_cluster_centroids(i,:)=sum(LV_xy(CL==i,:))/length(LV_xy(CL==i,:))
        end
    end
else
    nc=length(cluster_control)
    if strcmp(clustering_type,'kmeans')
        cluster_control
        nc
        
        [CL,final_cluster_centroids]=kmeans(LV_xy,nc,'Replicates',1,'Start',cluster_control);%,'Options',opts);
    else
        %add other clustering strategies here
    end
end
%%%%%%%%%%%%%%%%%%%%  PLOT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(235235)

%~~~~~~~~kmeans
if strcmp(clustering_type,'kmeans')
    plot(final_cluster_centroids(:,1),final_cluster_centroids(:,2),'xk','Linewidth',20)
    hold all
    if numel(cluster_control)>=2
        plot(cluster_control(:,1),cluster_control(:,2),'xk','Linewidth',5)
    end
    plot(xbreg,ybreg)
    for i=1:length(final_cluster_centroids)
        plot(LV_xy(CL==i,1),LV_xy(CL==i,2),'o')
    end
    voronoi(final_cluster_centroids(:,1),final_cluster_centroids(:,2))
    hold off
    ylim([plotlimits(2,1),plotlimits(2,2)])
    xlim([plotlimits(1,1),plotlimits(1,2)])
    legend('Final cluster centroids','Initial cluster centroids')
    xlabel('Longitude [deg]')
    ylabel('Latitude [deg]')
    drawnow
elseif strcmp(clustering_type,'hierarchy')
    plot(final_cluster_centroids(:,1),final_cluster_centroids(:,2),'xk','Linewidth',20)
    hold all
    for i=1:max(CL)
        plot(LV_xy(CL==i,1),LV_xy(CL==i,2),'.')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create output index cellArray with one cell per cluster
index_cluster_cellarray={};
for i=1:nc
    index_cluster_cellarray{i}=LV_ind_list(CL==i);
end

%get a list of all possible pairs of cluster indices
cpairs = nchoosek(1:size(final_cluster_centroids,1), 2);

%calculate the distance between each pair of cluster centroids and sort them
cdist = flipdim(sortrows([sum((final_cluster_centroids( cpairs(:, 1), :)-final_cluster_centroids( cpairs(:, 2), :)).^2,2),cpairs],1),1);

%get the nf combinations whos centroids are furthest apart
cdist_f=cdist(1:nf,2:3);

