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

function [index_cluster_cellarray_allocated,pop_cluster,Calloc_clusters_pv,Calloc_clusters_csp,cluster_area]=allocate_capacity_clusters_solar(alloc_strat_clusters, alloc_strat_pixels, X, Y, xbreg, ybreg, pop_mask_priorit, pixelarea, popdens_mat, datamatrix_pv0, datamatrix_csp0, index_cluster_cellarray, index_country_cellarray_region, Cmax_pv, Cmax_csp, areadens_pv, areadens_csp, fraction_pv, final_cluster_centroids, value_thresh, frac_alloc_min, radiation_to_ac_eff_pv, radiation_to_ac_eff_csp)
%%
% alloc_strat_clusters;
% alloc_strat_pixels;
% X;
% Y;
% xbreg;
% ybreg;
% pop_mask_priorit;
% pixelarea;
% popdens_mat;
% datamatrix_pv0;
% datamatrix_csp0;
% index_cluster_cellarray=index_cluster_cellarray_onsh;
% index_country_cellarray_region;
% Cmax_pv=Cmax_total_pv;
% Cmax_csp=Cmax_total_csp;
% areadens_pv;
% areadens_csp;
% fraction_pv;
% final_cluster_centroids;
% value_thresh;
% frac_alloc_min;
% radiation_to_ac_eff_pv;
% radiation_to_ac_eff_csp;

%%
global plotlimits
global plotfilename
global output_dir_graph

%CALCULATE THE TOTAL AREA OF THE RELEVANT PIXELS FOR EACH
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
    Cmax_cluster_cells(:)=cluster_area(:)/sum(cluster_area)
%OR: Allocate the total Cmax to the different clusters equally
elseif strcmp(alloc_strat_clusters,'cl_equal')
    'cl_equal'
    Cmax_cluster_cells(:)=1/length(index_cluster_cellarray);
elseif strcmp(alloc_strat_clusters,'cl_popprop') || strcmp(alloc_strat_clusters,'cl_powpopprop')
    'cl_popprop'
    
    %Get all points within each cluster and calculate the total
    %population within that cluster
    index_region_all=[index_country_cellarray_region{:}].';

    A=[X(index_region_all),Y(index_region_all)];
    
    %calculate the distance between all pixels of the region and all the
    %final cluster centroids
    distance_points_cluster_centroids = pdist2(A,final_cluster_centroids);

    %determine the index of the cluster centroid closest to each pixel
    [~,clindex_mindist]=min(distance_points_cluster_centroids,[],2);

    index_cluster_region_region={};
    pop_cluster=zeros(length(index_cluster_cellarray),1);
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
        Cmax_cluster_cells=pop_cluster(:)/sum(pop_cluster);
    elseif strcmp(alloc_strat_clusters,'cl_powpopprop')
        %Calculate the total average solar resource (CSP) [kWh/yr] within cluster i
        power_cluster=zeros(length(index_cluster_cellarray),1);
        for i=1:length(index_cluster_cellarray)
            power_cluster(i)=sum(datamatrix_csp0(index_cluster_all_points_region{i}).*pixelarea(index_cluster_all_points_region{i}));
        end
        
        cluster_allocation_powpop=0.3
        
        %Allocate the total region capacity among the clusters in
        %proportion to the average fraction of population and solar
        %resource
        Cmax_cluster_cells=(cluster_allocation_powpop*pop_cluster(:)/sum(pop_cluster)+(1-cluster_allocation_powpop)*power_cluster(:)/sum(power_cluster));

    end
end

Cmax_cluster_cells_pv=Cmax_cluster_cells*Cmax_pv
Cmax_cluster_cells_csp=Cmax_cluster_cells*Cmax_csp

frac_pop=pop_cluster(:)/sum(pop_cluster);
frac_pow=power_cluster(:)/sum(power_cluster);

if sum(Cmax_cluster_cells_pv)>0
    frac_Cmax_cluster_cells_pv=Cmax_cluster_cells_pv/sum(Cmax_cluster_cells_pv)
end
if sum(Cmax_cluster_cells_csp)>0
    frac_Cmax_cluster_cells_csp=Cmax_cluster_cells_csp/sum(Cmax_cluster_cells_csp)
end





%%
Calloc_clusters_pv=zeros(size(index_cluster_cellarray));
Calloc_clusters_csp=zeros(size(index_cluster_cellarray));
PV_yearly_energy_clusters_cells={};
CSP_yearly_energy_clusters_cells={};

% Calculate the yearly energy generation for all relevant pixels of the
% cluster (CSP and PV)
index_cluster_cellarray_allocated={};
for i=1:length(index_cluster_cellarray)
    vector_powercum_pindex_cn_PV_above_pop_thresh=[];

    Calloc_pv=0;
    Calloc_csp=0;
    %                                   |kWh/m_panel^2/yr                            km_surf^2                             m_surf^2/km_surf^2  |m_panel^2/m_surf^2  |EJ/KWh     |kWh_ac/kWh_solar  
    PV_yearly_energy_clusters_cells{i} =datamatrix_pv0(index_cluster_cellarray{i}).* pixelarea(index_cluster_cellarray{i})*10^6*               areadens_pv*         3.6*10^-12* radiation_to_ac_eff_pv ; %[EJ/yr]
    CSP_yearly_energy_clusters_cells{i}=datamatrix_csp0(index_cluster_cellarray{i}).*pixelarea(index_cluster_cellarray{i})*10^6*               areadens_csp*        3.6*10^-12* radiation_to_ac_eff_csp; %[EJ/yr]
%                                       |kW/m_panel^2-----------------------------|
%                                       |kW/m_panel^2*km_surf^2----------------------------------------------------------|
%                                       |kW/m_panel^2*m_surf^2-----------------------------------------------------------------|
%                                       |kW--------------------------------------------------------------------------|
%                                       |kW_ac------------------------------------------------------------------------------------------------------|

    %                              |index                                    |yearly energy output pixel         |pixel matrix index         |population allocation priority     
    vector_powercum_pindex_cn_PV= [(1:size(index_cluster_cellarray{i},1)).', PV_yearly_energy_clusters_cells{i}, index_cluster_cellarray{i}, pop_mask_priorit(index_cluster_cellarray{i})]
    vector_powercum_pindex_cn_CSP=[(1:size(index_cluster_cellarray{i},1)).', CSP_yearly_energy_clusters_cells{i},index_cluster_cellarray{i}, pop_mask_priorit(index_cluster_cellarray{i})]

    %first sort by population priority, then sort according to the power production in each cell
    vector_powercum_pindex_cn_sorted_PV=flip(sortrows(sortrows(vector_powercum_pindex_cn_PV, 2), 4), 1)
    vector_powercum_pindex_cn_sorted_CSP=flip(sortrows(sortrows(vector_powercum_pindex_cn_CSP, 2), 4), 1)

    % Get pixels with population density above threshold
    vector_powercum_pindex_cn_PV_above_pop_thresh=vector_powercum_pindex_cn_sorted_PV(vector_powercum_pindex_cn_sorted_PV(:,4)==2,:)
    vector_powercum_pindex_cn_PV_below_pop_thresh=vector_powercum_pindex_cn_sorted_PV(vector_powercum_pindex_cn_sorted_PV(:,4)==1,:)
    vector_powercum_pindex_cn_CSP_below_pop_thresh=vector_powercum_pindex_cn_sorted_CSP(vector_powercum_pindex_cn_sorted_CSP(:,4)==1,:)

    %add a column containing the cumulative energy of the pixels with population density above threshold
    vector_powercum_pindex_cn_PV_above_pop_thresh=[vector_powercum_pindex_cn_PV_above_pop_thresh,cumsum(vector_powercum_pindex_cn_PV_above_pop_thresh(:,2))]
    vector_powercum_pindex_cn_PV_below_pop_thresh=[vector_powercum_pindex_cn_PV_below_pop_thresh,cumsum(vector_powercum_pindex_cn_PV_below_pop_thresh(:,2))]
    vector_powercum_pindex_cn_CSP_below_pop_thresh=[vector_powercum_pindex_cn_CSP_below_pop_thresh,cumsum(vector_powercum_pindex_cn_CSP_below_pop_thresh(:,2))]

    %for PV on pixels above the population threshold ("urban PV") cut all the pixels causing energy output above the threshold
    vector_powercum_pindex_cn_PV_above_pop_thresh_cut=vector_powercum_pindex_cn_PV_above_pop_thresh(vector_powercum_pindex_cn_PV_above_pop_thresh(:,5)<=Cmax_cluster_cells_pv(i),:)

    %Total energy allocated to pixels with above threshold population
    %density
    if isempty(vector_powercum_pindex_cn_PV_above_pop_thresh_cut)
        Calloc_pv=0;
    else
        Calloc_pv=vector_powercum_pindex_cn_PV_above_pop_thresh_cut(end,5);
    end
    
    %define the remaining pixels as allocated
    %                                     |indices of allocated pixels                            |code: 1-PV Urban, 2-PV, 3-CSP            
    index_cluster_cellarray_allocated{i}=[vector_powercum_pindex_cn_PV_above_pop_thresh_cut(:,3), 1*ones(size(vector_powercum_pindex_cn_PV_above_pop_thresh_cut(:,3))), ...
        vector_powercum_pindex_cn_PV_above_pop_thresh(vector_powercum_pindex_cn_PV_above_pop_thresh(:,5)<=Cmax_cluster_cells_pv(i),5), ...
        vector_powercum_pindex_cn_PV_above_pop_thresh(vector_powercum_pindex_cn_PV_above_pop_thresh(:,5)<=Cmax_cluster_cells_pv(i),2)]

    %continue allocating if Cmax_total has not been reached yet
    while Calloc_pv < Cmax_cluster_cells_pv(i) 
        if isempty(vector_powercum_pindex_cn_PV_below_pop_thresh(:,1))
            break
        end

        if sum(vector_powercum_pindex_cn_CSP_below_pop_thresh(:,3)==vector_powercum_pindex_cn_PV_below_pop_thresh(1,3))~=0 %next best index in pv also present in csp
            rand_fraction_pv=rand
            if rand_fraction_pv<fraction_pv %allocate to PV with probability fraction_pv and delete index in the csp list:
                Calloc_pv=Calloc_pv+vector_powercum_pindex_cn_PV_below_pop_thresh(1,2);
                pixel_index=vector_powercum_pindex_cn_PV_below_pop_thresh(1,3);                         %|code: 1-PV Urban, 2-PV, 3-CSP
                index_cluster_cellarray_allocated{i}=[index_cluster_cellarray_allocated{i}; [pixel_index,2, Calloc_pv, vector_powercum_pindex_cn_PV_below_pop_thresh(1,2)]]; %set pixel as allocated
                % Delete the newly allocated pixel from sorted power
                % vectors both pv and csp
                vector_powercum_pindex_cn_PV_below_pop_thresh=vector_powercum_pindex_cn_PV_below_pop_thresh(vector_powercum_pindex_cn_PV_below_pop_thresh(:,3)~=pixel_index,:);
                vector_powercum_pindex_cn_CSP_below_pop_thresh=vector_powercum_pindex_cn_CSP_below_pop_thresh(vector_powercum_pindex_cn_CSP_below_pop_thresh(:,3)~=pixel_index,:);
            else %delete in the pv list:
                vector_powercum_pindex_cn_PV_below_pop_thresh(1,3)
                vector_powercum_pindex_cn_PV_below_pop_thresh=vector_powercum_pindex_cn_PV_below_pop_thresh(2:end,:);                
            end
        else %allocate to pv
            Calloc_pv=Calloc_pv+vector_powercum_pindex_cn_PV_below_pop_thresh(1,2);
            pixel_index=vector_powercum_pindex_cn_PV_below_pop_thresh(1,3);                         %|code: 1-PV Urban, 2-PV, 3-CSP
            index_cluster_cellarray_allocated{i}=[index_cluster_cellarray_allocated{i}; [pixel_index,2, Calloc_pv, vector_powercum_pindex_cn_PV_below_pop_thresh(1,2)]]; %set pixel as allocated
            % Delete the newly allocated pixel from sorted power
            % vectors both pv and csp
            vector_powercum_pindex_cn_PV_below_pop_thresh=vector_powercum_pindex_cn_PV_below_pop_thresh(vector_powercum_pindex_cn_PV_below_pop_thresh(:,3)~=pixel_index,:);
        end
    end
    while Calloc_csp < Cmax_cluster_cells_csp(i)
        if isempty(vector_powercum_pindex_cn_CSP_below_pop_thresh(:,1))
            break
        end       
        Calloc_csp=Calloc_csp+vector_powercum_pindex_cn_CSP_below_pop_thresh(1,2);
        pixel_index=vector_powercum_pindex_cn_CSP_below_pop_thresh(1,3);
        index_cluster_cellarray_allocated{i}=[index_cluster_cellarray_allocated{i}; [pixel_index,3, Calloc_csp, vector_powercum_pindex_cn_CSP_below_pop_thresh(1,2)]]; %set pixel as allocated 
        vector_powercum_pindex_cn_CSP_below_pop_thresh=vector_powercum_pindex_cn_CSP_below_pop_thresh(vector_powercum_pindex_cn_CSP_below_pop_thresh(:,3)~=pixel_index,:);
    end

    [i, Calloc_pv, Cmax_cluster_cells_pv(i), Calloc_csp, Cmax_cluster_cells_csp(i)]
    Calloc_clusters_pv(i)=Calloc_pv;
    Calloc_clusters_csp(i)=Calloc_csp;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%      PLOT      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=262467;
f=figure(nfig)
%screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );
colormap(gray)
pcolorjw(X,Y,datamatrix_csp0)%max(max(datamatrix_csp0))-datamatrix_csp0)
h=colorbar;
h.Label.String='Yearly Direct Irradiation [kWh/m^2/yr]'
hold all

index_cluster_cellarray_allocated_pv={};
index_cluster_cellarray_allocated_pv_urban={};
index_cluster_cellarray_allocated_csp={};
% Sort allocated pixels by CSP/PV code: 1-PV Urban, 2-PV, 3-CSP  
for i=1:length(index_cluster_cellarray_allocated)
    index_cluster_cellarray_allocated_pv_urban{i}=index_cluster_cellarray_allocated{i}(index_cluster_cellarray_allocated{i}(:,2)==1,1)
    index_cluster_cellarray_allocated_pv{i}=index_cluster_cellarray_allocated{i}(index_cluster_cellarray_allocated{i}(:,2)==2,1)
    index_cluster_cellarray_allocated_csp{i}=index_cluster_cellarray_allocated{i}(index_cluster_cellarray_allocated{i}(:,2)==3,1)
end
msize=3;
plot(X(vertcat(index_cluster_cellarray{:})),Y(vertcat(index_cluster_cellarray{:})),'og','markersize',msize);

h=plot(X(vertcat(index_cluster_cellarray_allocated_pv_urban{:})),Y(vertcat(index_cluster_cellarray_allocated_pv_urban{:})),'oy');
set(h,'MarkerEdgeColor','y','MarkerFaceColor','y','markersize',msize)
h=plot(X(vertcat(index_cluster_cellarray_allocated_pv{:})),Y(vertcat(index_cluster_cellarray_allocated_pv{:})),'or');
set(h,'MarkerEdgeColor','r','MarkerFaceColor','r','markersize',msize)
h=plot(X(vertcat(index_cluster_cellarray_allocated_csp{:})),Y(vertcat(index_cluster_cellarray_allocated_csp{:})),'ob');
set(h,'MarkerEdgeColor','b','MarkerFaceColor','b','markersize',msize)

%plot(X(vertcat(LV_index_cluster_cellarray_allocated_all_together{:})),Y(vertcat(LV_index_cluster_cellarray_allocated_all_together{:})),'gx')
legend_entries={'CSP kWh/m2/yr','Filtered pixels','Allocated: PV Urban','Allocated: PV','Allocated: CSP'}

%title(strcat(plotfilename, '; CmaxPV=',num2str(Cmax_pv), '; CmaxCSP=',num2str(Cmax_csp),'EJ; ','kW/km2 ; Avg wind speed threshold=',num2str(value_thresh),'m/s;',alloc_strat_clusters,'; ',alloc_strat_pixels,'; Minimum fraction pixels allocated=',num2str(frac_alloc_min)),'interpreter','none')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
ylim([plotlimits(2,1),plotlimits(2,2)]) 
xlim([plotlimits(1,1),plotlimits(1,2)])

for i=1:length(index_cluster_cellarray)
    plot(final_cluster_centroids(i,1),final_cluster_centroids(i,2),'.w')
    text(final_cluster_centroids(i,1),final_cluster_centroids(i,2),['\fontsize{30}{\color{blue}' num2str(i) '}'])

%     legend_entries={legend_entries{:},[num2str(i) ...
%         'Population: ' num2str(pop_cluster(i)/1e6, '%0.2f') 'M; ' ...
%         'mPV' num2str(Cmax_cluster_cells_pv(i),'%0.2f') 'EJ; '...
%         'aPV' num2str(Calloc_clusters_pv(i),'%0.2f') 'EJ; ' ...
%         'mCSP' num2str(Cmax_cluster_cells_csp(i),'%0.2f') 'EJ; '...
%         'aCSP' num2str(Calloc_clusters_csp (i),'%0.2f') 'EJ; ' ...
%         '#CSP: ' num2str(length(index_cluster_cellarray_allocated_csp{i})) ...
%         '#PV urban: ' num2str(length(index_cluster_cellarray_allocated_pv_urban{i})) ...
%         '#PV: '  num2str(length(index_cluster_cellarray_allocated_pv{i})) ...
%         ]};%; Number px cell: ' num2str(length(index_cluster_cellarray{i})) '; Number alloc: ' num2str(length(index_cluster_cellarray_allocated{i})) '; Fraction alloc: ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}))]};

%     legend_entries={legend_entries{:},[num2str(i),' ', ...
%         'Population: ', num2str(pop_cluster(i)/1e6, '%0.2f'), 'M; ', ...
%         'PV energy', num2str(Cmax_cluster_cells_pv(i),'%0.2f'), 'EJ/yr; ',...
%         'CSP energy', num2str(Cmax_cluster_cells_csp(i),'%0.2f'), 'EJ/yr; ',...
%         ]};%; Number px cell: ' num2str(length(index_cluster_cellarray{i})) '; Number alloc: ' num2str(length(index_cluster_cellarray_allocated{i})) '; Fraction alloc: ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}))]};
    legend_entries={legend_entries{:},[num2str(i) ' Population: ' num2str(frac_pop(i)*100,'%0.1f') '%' ' Solar Resource: ' num2str(frac_pow(i)*100,'%0.1f') '%' ' Allocation: ' num2str(frac_Cmax_cluster_cells_pv(i)*100,'%0.1f') '%']};%; Number px cell: ' num2str(length(index_cluster_cellarray{i})) '; Number alloc: ' num2str(length(index_cluster_cellarray_allocated{i})) '; Fraction alloc: ' num2str(length(index_cluster_cellarray_allocated{i})/length(index_cluster_cellarray{i}))]};


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
