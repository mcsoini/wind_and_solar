function [countryCellArray,plotfilename,excludeCoords,plotlimits,initial_cluster_centroids_onsh,initial_cluster_centroids_offsh]=get_region_parameters(region_str)
'... entered get_region_parameters ...'

switch region_str
    case {'FSU'};  plotfilename='FSU_' ; excludeCoords=[nan,0,nan,nan]  ; plotlimits=[[18,180];[35,80]];
    case {'CPA'};  plotfilename='CPA_' ; excludeCoords=[nan,nan,nan,nan]; plotlimits=[[72.61,135.75];[7.58,54.56]];
    case {'PAS'};  plotfilename='PAS_' ; excludeCoords=[nan,0,nan,nan]  ; plotlimits=[[91.18,181.00];[-23.66,39.62]];
    case {'LAC'};  plotfilename='LAC_' ; excludeCoords=[nan,nan,nan,nan]; plotlimits=[[-120,-32];[-60,40]];
    case {'MEA'};  plotfilename='MEA_' ; excludeCoords=[nan,nan,nan,nan]; plotlimits=[[-18.10,64.31];[2.49,40.77]];
    case {'NAM'};  plotfilename='NAM_' ; excludeCoords=[nan,nan,nan,nan]; plotlimits=[[-175,-52];[18.5,85]];
    case {'AFR'};  plotfilename='AFR_' ; excludeCoords=[nan,nan,nan,nan]; plotlimits=[[-25,58];[-35,27]];
    case {'SAS'};  plotfilename='SAS_' ; excludeCoords=[nan,nan,nan,nan]; plotlimits=[[59.49,98.34];[2.23,39.46]];
    case {'PAO'};  plotfilename='PAO_' ; excludeCoords=[nan,0,nan,nan]  ; plotlimits=[[111.90,179.54];[-55.75,46.51]];
    case {'TEU'};  plotfilename='TEU_' ; excludeCoords=[nan,nan,nan,20] ; plotlimits=[[-40,41];[28,71.5]];
end

%Define countryCellArrays containing the country names of the regions to be
%considered (these need to be consistent with the names in the country border data set!)

%ALSO: Define initial values for the cluster centroids for each region
%First: Get a specified number n of clusters for random starting centroid
%values by setting the variable cluster_control=n; then, set
%initial_cluster_centroids_onsh (or _offsh) equal the
%final_cluster_centroids, for reproducibility; set cluster_control=initial_cluster_centroids_onsh
%-> get from % fprintf('%.2f,%.2f];[', final_cluster_centroids.')
switch region_str
    case {'FSU'}
        % FSU
        initial_cluster_centroids_onsh=[[37.54,56.49];[100.47,55.25];[139.70,53.38];[75.84,55.92];[58.11,50.98]];
        initial_cluster_centroids_offsh=[[32.97,48.13];[42.98,67.11];[50.67,42.57];[144.49,51.96]];
        countryCellArray={'Armenia','Azerbaijan','Belarus','Georgia','Kazakhstan','Kyrgyzstan','Moldova','Russia','Tajikistan','Turkmenistan','Ukraine','Uzbekistan'}
    case {'CPA'}
        % CPA
        initial_cluster_centroids_onsh=[[92.81,32.77];[107.57,40.02];[108.37,23.01];[119.60,40.32];[112.75,45.77];[84.22,34.58];[128.70,45.75];[102.12,44.80];[93.38,42.88];[122.11,48.28]];
        initial_cluster_centroids_offsh=[[121.79,35.81];[118.13,23.95];[108.62,20.27];[105.31,8.69]];
        countryCellArray={'Cambodia','China','Dem. Rep. Korea','Lao PDR','Mongolia','Vietnam'}
    case {'PAS'}
        % PAS
        initial_cluster_centroids_onsh=[[124.56,19.59];[101.57,15.29];[152.38,-11.94]];
        initial_cluster_centroids_offsh=[[121.65,12.24];[125.91,34.21];[105.51,2.95];[139.90,-8.43]];
        countryCellArray={'American Samoa','Brunei Darussalam','Fiji','French Polynesia','Gilbert-Kiribati','Indonesia','Malaysia','Myanmar','New Caledonia','Papua New Guinea','Philippines','Korea','Singapore','Solomon Islands','Taiwan (China)','Thailand','Tonga','Vanuatu','Western Samoa'}
    case {'LAC'}
        % LAC
        initial_cluster_centroids_onsh=[[-84.42,18.53];[-104.04,26.95];[-66.97,7.81];[-40.88,-8.74];[-58.92,-27.27];[-67.39,-40.71]];
        initial_cluster_centroids_offsh=[[-79.05,23.82];[-90.51,21.38];[-79.79,-8.04];[-114.50,29.21];[-57.90,7.87];[-61.17,-41.24];[-81.18,13.78];[-45.16,-25.11];[-46.70,0.22]];
        countryCellArray={'Antigua and Barbuda','Argentina','Bahamas','Barbados','Belize','Bermuda','Bolivia','Brazil','Chile','Colombia','Costa Rica','Cuba','Dominica','Dominican Republic','Ecuador','El Salvador','French Guyana','Grenada','Guadeloupe','Guatemala','Guyana','Haiti','Honduras','Jamaica','Martinique','Mexico','Netherlands Antilles','Nicaragua','Panama','Paraguay','Peru','Saint Kitts and Nevis','Santa Lucia','Saint Vincent and the Grenadines','Suriname','Trinidad and Tobago','Uruguay','Venezuela'}
    case {'MEA'}
        % MEA
        initial_cluster_centroids_onsh=[[26.63,27.84];[-6.40,28.55];[7.92,31.61];[30.32,15.81];[41.97,29.58];[51.89,23.65]];
        initial_cluster_centroids_offsh=[[-14.91,25.23];[19.07,32.50];[50.92,27.13]];
        countryCellArray={'Palestine','W. Sahara','S. Sudan','Algeria','Bahrain','Egypt','Iraq','Iran','Israel','Jordan','Kuwait','Lebanon','Libya','Morocco','Oman','Qatar','Saudi Arabia','Sudan','Syria','Tunisia','United Arab Emirates','Yemen'}
    case {'NAM'}
        % NAM
        initial_cluster_centroids_onsh=[[-64.00,49.82];[-83.08,35.17];[-105.27,37.99];[-110.31,50.85];[-88.44,45.36];[-120.27,50.90];[-95.44,34.35];[-75.58,46.43];[-152.53,62.28];[-98.06,48.66]]; % NAM 10
        initial_cluster_centroids_offsh=[[-68.78,48.34];[-107.79,49.53];[-85.62,42.61];[-98.89,36.25];[-154.23,60.35]];
        countryCellArray={'Canada','Guam','Puerto Rico','United States','Virgin Islands'}
    case {'AFR'}
        % AFR
        initial_cluster_centroids_onsh=[[40.51,6.69];[26.14,-25.67];[-10.49,17.02];[15.63,13.61];[1.56,14.96]];
        initial_cluster_centroids_offsh=[[21.30,-35.70];[47.19,-11.66];[-17.12,15.08]];
        countryCellArray={'Somaliland','Dem. Rep. Congo','Angola','Benin','Botswana','Br. Indian Ocean Ter.','Burkina Faso','Burundi','Cameroon','Cape Verde','Central African Rep.','Chad','Comoros','Cote d''Ivoire','Congo','Djibouti','Eq. Guinea','Eritrea','Ethiopia','Gabon','Gambia','Ghana','Guinea','Guinea-Bissau','Kenya','Lesotho','Liberia','Madagascar','Malawi','Mali','Mauritania','Mauritius','Mozambique','Namibia','Niger','Nigeria','Reunion','Rwanda','São Tomé and Principe','Senegal','Seychelles','Sierra Leone','Somalia','South Africa','Saint Helena','Swaziland','Tanzania','Togo','Uganda','Zambia','Zimbabwe'}
    case {'SAS'}
        % SAS
        initial_cluster_centroids_onsh=[[77.83,15.48];[80.71,9.00];[70.46,24.80];[75.74,20.61];[62.74,30.61]];
        initial_cluster_centroids_offsh=[[79.67,9.75];[68.50,22.67];[90.04,19.29];[71.57,18.82]];
        countryCellArray={'Afghanistan','Bangladesh','Bhutan','India','Maldives','Nepal','Pakistan','Sri Lanka'}
    case {'PAO'}
        % PAO
        initial_cluster_centroids_onsh=[[138.90,38.83];[173.03,-41.44];[146.03,-20.96];[117.40,-29.75];[144.20,-34.03]];
        initial_cluster_centroids_offsh=[[141.89,-16.33];[172.52,-41.23];[135.50,-35.62];[133.37,35.63]];
        countryCellArray={'Australia','Japan','New Zealand'}
    case {'TEU'}
        % TEU
        initial_cluster_centroids_onsh=[...
            [18.244,51.457];...
            [-18.449,64.968];...
            [29.95,40.672];...
            [5.4923,48.117];...
            [13.269,61.569];...
            [24.822,64.752];...
            [-4.6341,40.654];...
            [-4.0756,53.398];...
            ]
        initial_cluster_centroids_offsh=[[16.5,57];[4.5,54.5];[-3,49];[21,42]]
        countryCellArrayWEU={'Andorra','Austria','Azores','Belgium','Canary Islands','Channel Islands','Cyprus','Denmark','Faeroe Is.','Finland','France','Germany','Gibraltar','Greece','Greenland','Iceland','Ireland','Isle of Man','Italy','Liechtenstein','Luxembourg','Madeira','Malta','Monaco','Netherlands','Norway','Portugal','Spain','Sweden','Switzerland','Turkey','United Kingdom'};
        countryCellArrayEEU={'Albania','Bosnia and Herz.','Bulgaria','Croatia','Czech Rep.','Estonia','Macedonia','Latvia','Lithuania','Hungary','Poland','Romania','Slovakia','Slovenia','Montenegro','Serbia','Kosovo'};
        countryCellArray={countryCellArrayEEU{:},countryCellArrayWEU{:}}
end

'... leaving get_region_parameters ...'