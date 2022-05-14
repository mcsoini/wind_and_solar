% CUM_CF_ARB_INDEX_LIST
%
% Calculates the time dependent cumulative capacity factor for an arbitrary
% list of indices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% INPUT:        
% indices:      indices of the datamat array pixels to calculate the aggregate cf for
% rawdata:      3D time resolved wind speed array
% pixelarea:    2D array size datamat -> areas [km2] of all pixels
% capdens_mat:  2D capacity density matrix
% pd:           wind speed to power data
%
% OUTPUT:
% cum_cf:       1D array time dependent cumulative capacity factor
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [cum_cf,E_tot,cap_clusterpair_final]=cum_cf_arb_index_list(indices,rawdata,pixelarea,capdens_mat,h_vec,pd)
%%
D=double(h_vec(2)-h_vec(1));

E_t=double(speed2power(rawdata(:,indices),pd)/100*(pixelarea(indices).*capdens_mat(indices))); %[W](t)

cap_clusterpair_final=sum(pixelarea(indices).*capdens_mat(indices))

E_tot=D*sum(E_t)*3600*1e-18; %Total energy produced [EJ]

cum_cf=E_t/sum(pixelarea(indices).*capdens_mat(indices));

% %%
% cum_cf=zeros(3,length(h_vec))
% cum_cf(1,:)=(speed2power(rawdata(:,index_cluster_cellarray{1}),pd)*(pixelarea(index_cluster_cellarray{1}).*capdens_mat(index_cluster_cellarray{1})))/sum(pixelarea(index_cluster_cellarray{1}).*capdens_mat(index_cluster_cellarray{1}))/100;
% 
% cum_cf(2,:)=(speed2power(rawdata(:,index_cluster_cellarray{2}),pd)*(pixelarea(index_cluster_cellarray{2}).*capdens_mat(index_cluster_cellarray{2})))/sum(pixelarea(index_cluster_cellarray{2}).*capdens_mat(index_cluster_cellarray{2}))/100;
% 
% indices=vertcat(index_cluster_cellarray{2},index_cluster_cellarray{1})
% cum_cf(3,:)=(speed2power(rawdata(:,indices),pd)*(pixelarea(indices).*capdens_mat(indices)))/sum(pixelarea(indices).*capdens_mat(indices))/100;
% 
% plot(h_vec,sortrows(cum_cf.',3))

