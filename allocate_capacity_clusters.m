% ALLOCATE_CAPACITY_CLUSTERS
%
% This allocates the total energy output to the different clusters
% according to different strategies; for each cluster it allocates the
% energy output to the pixels according to different strategies
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:
% alloc_strat_clusters:                               
% alloc_strat_pixels:                                 
% X:                                                  latitude meshgrid
% Y:                                                  longitude meshgrid
% xbreg:                                              1D array x-values of the region's border points
% ybreg:                                              1D array y-values of the region's border points
% datamat_orig:                                       2D array wind data
% index_cluster_cellarray:                            cellarray of 1-column vectors per cluster; contains the datamat-indices of all points within the cluster
% pd:                                                 2D array: map speed->power
% pixelarea:                                          2D array area [km²] of all pixels
% Cmax:                                               total capacity allocated to that region
% capdens:                                            capacity density
% final_cluster_centroids:                            cluster centroids (only needed for plotting)
% speed_thresh:                                       probably obsolete
%
% OUTPUT:
% index_cluster_cellarray_allocated:                  cellarray of 1D index arrays per cluster, defining those pixels where wind power has been allocated
% capdens_mat:                                        2D array size(X) containing the capacity density for each pixel (constant capdens if allocated according to highest wind speeds)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [index_cluster_cellarray_allocated,capdens_mat,pop_cluster,Cmax_cluster_cells,cluster_area]=allocate_capacity_clusters(alloc_strat_clusters,alloc_strat_pixels,X,Y,xbreg,ybreg,cf_time_avg_mat,index_cluster_cellarray,pd,pixelarea,Cmax,capdens,popdens_mat,index_country_cellarray_region,final_cluster_centroids,speed_thresh,frac_alloc_min,poppropctrl)

% % %%
% frac_alloc_min=0.5
% alloc_strat_clusters='cl_powpopprop'
% alloc_strat_pixels='px_windrank'
% index_cluster_cellarray=index_cluster_cellarray_offsh
% Cmax=Cmax_offsh
% capdens=capdens_offsh
% popdens_mat=popdens_mat
% %final_cluster_centroids=final_cluster_centroids_onsh
% speed_thresh=5
% cf_time_avg_mat=cf_time_avg_offsh_mat

global plotlimits
global plotfilename
global output_dir_graph

%CALCULATE THE TOTAL AREA OF THE POTENTIALLY ALLOCATED PIXELS FOR EACH
%CLUSTER
cluster_area=zeros(size(index_cluster_cellarray.'))
for i=1:length(index_cluster_cellarray)
    cluster_area(i)=sum(pixelarea(index_cluster_cellarray{i}))
end

%Allocate the total Cmax to the different clusters in proportion to their
%total area
Cmax_cluster_cells=zeros(size(index_cluster_cellarray.'));
pop_cluster=zeros(size(index_cluster_cellarray.'));
if strcmp(alloc_strat_clusters,'cl_areaprop')
    'cl_areaprop'
    Cmax_cluster_cells(:)=cluster_area(:)/sum(cluster_area)*Cmax
%OR: Allocate the total Cmax to the different clusters equally
elseif strcmp(alloc_strat_clusters,'cl_equal')
    'cl_equal'
    Cmax_cluster_cells(:)=Cmax/length(index_cluster_cellarray);
elseif strcmp(alloc_strat_clusters,'cl_popprop') || strcmp(alloc_strat_clusters,'cl_powpopprop')
    'cl_popprop or cl_powpopprop'
    
    %Get all points within each cluster and calculate the total
    %population within that cluster
        
    index_region_all=[index_country_cellarray_region{:}].'

    A=[X(index_region_all),Y(index_region_all)];
    
    %calculate the distance between all pixels of the region and all the
    %final cluster centroids
    distance_points_cluster_centroids = pdist2(A,final_cluster_centroids);

    %determine the index of the cluster centroid closest to each pixel
    [~,clindex_mindist]=min(distance_points_cluster_centroids,[],2)

    index_cluster_region_region={}    
    pop_cluster=zeros(length(index_cluster_cellarray),1)
    for i=1:length(index_cluster_cellarray)
        %assign the pixels closest to the centroid to the corresponding
        %cluster
        index_cluster_all_points_region{i}=index_region_all(clindex_mindist==i);

        %Calculate the total population in cluster i
        pop_cluster(i)=sum(popdens_mat(index_cluster_all_points_region{i}).*pixelarea(index_cluster_all_points_region{i}));
    end
 
    
    if strcmp(alloc_strat_clusters,'cl_popprop')
        %Allocate the total region capacity among the clusters in
        %proportion to the population
        Cmax_cluster_cells=Cmax*pop_cluster(:)/sum(pop_cluster);
    elseif strcmp(alloc_strat_clusters,'cl_powpopprop')
        %Calculate the total average solar resource (CSP) [kWh/yr] within cluster i
        power_cluster=zeros(length(index_cluster_cellarray),1);
        for i=1:length(index_cluster_cellarray)
            power_cluster(i)=sum(cf_time_avg_mat(index_cluster_cellarray{i}).*pixelarea(index_cluster_cellarray{i}));

%              [i power_cluster(i) sum(pixelarea(index_cluster_cellarray{i}))]
%              plot(X(index_cluster_cellarray{i}),Y(index_cluster_cellarray{i}),'o')
%              pause
%             power_cluster(i)=power_cluster(i)./sum(pixelarea(index_cluster_all_points_region{i}));
        end
       % poppropctrl=0.3

        %Allocate the total region capacity among the clusters in
        %proportion to the average fraction of population and solar
        %resource
        Cmax_cluster_cells=(poppropctrl*pop_cluster(:)/sum(pop_cluster)+(1-poppropctrl)*power_cluster(:)/sum(power_cluster));
                
        Cmax_cluster_cells=Cmax_cluster_cells*Cmax;
    end
end

frac_pop=pop_cluster(:)/sum(pop_cluster);
frac_pow=power_cluster(:)/sum(power_cluster);

Cmax_cluster_cells_frac=Cmax_cluster_cells/sum(Cmax_cluster_cells);

%Define the capacity density for the relevant pixels
capdens_mat=zeros(size(X));
% Option 1: The capacity is allocated to the pixels with the highest wind
% speeds; this means that we just use a constant
if strcmp(alloc_strat_pixels,'px_windrank')
    capdens_mat(vertcat(index_cluster_cellarray{:}))=capdens;
% Option 2: The capacity in each cluster is allocated within each region
% with constant capacity density; this requires a different capacity density for
% each region
elseif strcmp(alloc_strat_pixels,'px_areaprop') 
    for i=1:length(index_cluster_cellarray)
        %The capacity density of all pixels of cluster i is the total
        %energy produced in that cluster [J] divided by the area and the
        %duration of a year [s/yr]
        capdens_mat(index_cluster_cellarray{i})=Cmax_cluster_cells(i)/(cluster_area(i)*8760*3600);
    end
end

%Allocate capacity according to highest average wind speeds using the
%function allocate_wind_at_pdensity (if the capacity density is defined such
%that the wind power is distibuted equally among all pixels, the resulting
%list of pixel indices should be the list of pixel indices of the whole region)
index_cluster_cellarray_allocated={};
Calloc_clusters=zeros(length(index_cluster_cellarray),1);
for i=1:length(index_cluster_cellarray)
    size(index_cluster_cellarray{i})
    [LV_index_cluster_cellarray_allocated,Calloc_clusters(i)]=allocate_wind_at_pdensity(X,Y,cf_time_avg_mat,{index_cluster_cellarray{i}},pd,pixelarea,capdens_mat,Cmax_cluster_cells(i),frac_alloc_min);
    index_cluster_cellarray_allocated{i}=LV_index_cluster_cellarray_allocated{:};
end

%Recalculation of the cluster capacities necessary if the parameter
%frac_alloc_min is not equal to zero; NOTE: for these new capacity
%densities Calloc_clusters=Cmax_cluster_cells
for i=1:length(index_cluster_cellarray)
    %The capacity density of all pixels of cluster i is the total
    %energy allocated to that cluster [J] divided by the area and the
    %duration of a year [s/yr]
    capdens_mat(index_cluster_cellarray{i})=Cmax_cluster_cells(i)/(sum(pixelarea(index_cluster_cellarray_allocated{i}))*8760*3600);
    [i Calloc_clusters(i)/(sum(pixelarea(index_cluster_cellarray_allocated{i}))*8760*3600) Calloc_clusters(i) (sum(pixelarea(index_cluster_cellarray_allocated{i}))*8760*3600)]
end

%calculate the allocation without clustering, for comparison; independent
%of the allocation in the whole region -> no new capacity density required
%(same for all pixels)
[LV_index_cluster_cellarray_allocated_all_together,Calloc_all_together]=allocate_wind_at_pdensity(X,Y,cf_time_avg_mat,index_cluster_cellarray,pd,pixelarea,capdens_mat,Cmax,frac_alloc_min);
'a'
sum(Calloc_clusters)
Calloc_all_together
pop_cluster

%%%%%%%%%%%%%      PLOT      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=262467;
f=figure(nfig)
%screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );
colormap(flipud(gray))
pcolorjw(X,Y,cf_time_avg_mat)
h=colorbar
h.Label.String='Average Capacity Factor'
hold all

msize=2;

h=plot(X(vertcat(index_cluster_cellarray{:})),Y(vertcat(index_cluster_cellarray{:})),'og');
set(h,'markersize',msize)
h=plot(X(vertcat(index_cluster_cellarray_allocated{:})),Y(vertcat(index_cluster_cellarray_allocated{:})),'ro')
set(h,'MarkerEdgeColor','r','MarkerFaceColor','r','markersize',msize)

%plot(X(vertcat(LV_index_cluster_cellarray_allocated_all_together{:})),Y(vertcat(LV_index_cluster_cellarray_allocated_all_together{:})),'gx')
legend_entries={'Time Averaged Capacity Factor','Filtered Pixels','Allocated: Wind'}

title(strcat(plotfilename, '; Cmax=',num2str(Cmax/1e18),'EJ; ','Orig. cap. density=',num2str(capdens/1e3),'kW/km2 ; Avg wind speed threshold=',num2str(speed_thresh),'m/s;',alloc_strat_clusters,'; ',alloc_strat_pixels,'; Minimum fraction pixels allocated=',num2str(frac_alloc_min)),'interpreter','none')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
ylim([plotlimits(2,1),plotlimits(2,2)])
xlim([plotlimits(1,1),plotlimits(1,2)])

for i=1:length(index_cluster_cellarray)
    plot(final_cluster_centroids(i,1),final_cluster_centroids(i,2),'.w')
    text(final_cluster_centroids(i,1),final_cluster_centroids(i,2),['\fontsize{30}{\color{blue}' num2str(i) '}'])
    if isempty(index_cluster_cellarray_allocated{i})
        capdens_mat_value='-';
    else
        capdens_mat_value=num2str(capdens_mat(index_cluster_cellarray_allocated{i}(1))/1000,'%0.2f')
    end
%     legend_entries={legend_entries{:},[num2str(i) ' Population: ' num2str(pop_cluster(i)/1e6,'%0.2f') 'M; Wind energy: ' num2str(Cmax_cluster_cells(i)/1e18,'%0.2f')  'EJ; ' capdens_mat_value 'kW/km2]; ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}),'%0.2f')]};%; Number px cell: ' num2str(length(index_cluster_cellarray{i})) '; Number alloc: ' num2str(length(index_cluster_cellarray_allocated{i})) '; Fraction alloc: ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}))]};
%    legend_entries={legend_entries{:},[num2str(i) ' Population: ' num2str(pop_cluster(i)/1e6,'%0.2f') 'M; Wind energy: ' num2str(Cmax_cluster_cells(i)/1e18,'%0.2f')  'EJ']};%; Number px cell: ' num2str(length(index_cluster_cellarray{i})) '; Number alloc: ' num2str(length(index_cluster_cellarray_allocated{i})) '; Fraction alloc: ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}))]};
    legend_entries={legend_entries{:},[num2str(i) ' Population: ' num2str(frac_pop(i)*100,'%0.1f') '%' ' Power: ' num2str(frac_pow(i)*100,'%0.1f') '%' ' Allocation: ' num2str(Cmax_cluster_cells_frac(i)*100,'%0.1f') '%']};%; Number px cell: ' num2str(length(index_cluster_cellarray{i})) '; Number alloc: ' num2str(length(index_cluster_cellarray_allocated{i})) '; Fraction alloc: ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}))]};
end

h_legend=legend(legend_entries)
set(h_legend,'FontSize',9);
plot(xbreg,ybreg)

[vorx, vory] =voronoi(final_cluster_centroids(:,1),final_cluster_centroids(:,2));
plot(vorx,vory,'k')
hold off
drawnow


%tightfig
%screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );
set(h_legend,'Location','Best');
set(gcf,'Position',[128,496,1271,448])

%print(f,'-dpng',strcat(output_dir_graph, plotfilename,'_cluster_cap_alloc.png'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




