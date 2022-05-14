% CUM_CF_ARB_INDEX_LIST_SOLAR
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


function [cf_t_pv,cf_t_csp,E_tot_pv,E_tot_csp]=cum_cf_arb_index_list_solar(indices_pv,indices_csp,rawdata_pv,rawdata_csp,pixelarea,radiation_to_ac_eff_pv, radiation_to_ac_eff_csp,areadens_pv, areadens_csp,h_vec)

D=double(h_vec(2)-h_vec(1));

%        |kW/m_panel^2            |km_surf^2            |m_surf^2/km_surf^2 |m_panel^2/m_surf^2 |kW_ac/kW_solar  
P_t_pv = rawdata_pv(:,indices_pv)*pixelarea(indices_pv)*10^6*               areadens_pv*        radiation_to_ac_eff_pv; %[kW_ac]
%        |kW/m_panel^2----------|
%        |kW/m_panel^2*km_surf^2----------------------|
%        |kW/m_panel^2*m_surf^2----------------------------|
%        |kW--------------------------------------------------------------------------|
%        |kW_ac------------------------------------------------------------------------------------------------------|

P_t_csp=rawdata_csp(:,indices_csp)*pixelarea(indices_csp)*10^6*areadens_csp*radiation_to_ac_eff_csp; %[kW_ac]

%Total energy (for consistency check)
E_tot_pv =sum(P_t_pv) *3.6e-12*D; %[EJ/yr]
E_tot_csp=sum(P_t_csp)*3.6e-12*D; %[EJ/yr]

%Normalize to 1
cf_t_csp=P_t_csp/(sum(pixelarea(indices_csp))*10^6*areadens_csp*radiation_to_ac_eff_csp);
cf_t_pv =P_t_pv /(sum(pixelarea(indices_pv))*10^6*areadens_pv*radiation_to_ac_eff_pv);