%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   CALCULATE THE AGGREGATE CF FOR THE PIXELS IN THE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [country_cmb_cum_cf,country_cmb_names]=combine_subregions(power_tdep_subregion,index_country_cellarray_allocated,nhighest_index,countryCellArray,nc,pixelarea,capdens,dt)

%%
% 
% 'TESTING COMBINE-SUBREGIONS'
% 
% power_tdep_subregion=power_tdep_subregion_offshore
% index_country_cellarray_allocated=index_country_cellarray_allocated_offsh
% nhighest_index
% countryCellArray
% nc
% pixelarea
% capdens
% dt
%possible combinations of a subset of m of the values nhighest_index;
%combinations returns a matrix were the each row corresponds to one
%distinct possibility of choosing m elements out of the nhighest_index
%nhighest_index_cmb=combntns(nhighest_index,nc)
nhighest_index_cmb=combntns(nhighest_index,nc)

%calculate combined capacity factor for each combination of regions
country_cmb_cum_cf=[];
country_cmb_names={};
for i=1:length(nhighest_index_cmb)
    LV_cumulative_power=[sum(cat(2,power_tdep_subregion{nhighest_index_cmb(i,:)}),2)]
    
    country_cmb_cum_cf=[country_cmb_cum_cf, LV_cumulative_power/sum(pixelarea(vertcat(index_country_cellarray_allocated{nhighest_index_cmb(i,:)})))/(capdens*dt*3600)]    
    country_cmb_names{i}=strjoin(countryCellArray{:}(nhighest_index_cmb(i,:)))
end








