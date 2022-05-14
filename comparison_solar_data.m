datafile2='../data/solar/interim_3rad_0p75deg_3_2013_01.grib';

varname_1='Surface_net_solar_radiation_clear_sky_surface';
varname_2='Surface_solar_radiation_surface';
varname_3='Surface_solar_radiation_downwards_surface';

dT=12; %According to "Select time"
dt=3;

'... entered getdata4Dm_solar ...'

% create ncgeodataset object, contains all the data from the file
nc=ncgeodataset(datafile2);

nctime=nc.geovariable('time');
timevec=nctime.data(:);
%%
% create geovariable objects
dirvar=nc.geovariable(varname_3);

% get the grid for time period 1
g=dirvar.grid_interop(1,:,:);

% extract the vectors containing the longitudes and latitudes
lon=g.lon;
lat=g.lat;   
[LON,LAT]=meshgrid(lon,lat);

rawdata=zeros([length(timevec), size(LON)]);

rawdata=dirvar.data(:,:,:);  %in [W*m^-2*s/Xh]

'convert to [kWs*m^-2/3h]'
rawdata0=rawdata;
parfor_progress(length(timevec));
for i=1:length(timevec)
    if mod(timevec(i)-dt,dT)==0 || i==1
        rawdata(i,:,:)=rawdata0(i,:,:);
    else
        rawdata(i,:,:)=rawdata0(i,:,:)-rawdata0(i-1,:,:);
    end
    parfor_progress;
end
parfor_progress(0);
clear rawdata0;
rawdata_ssrd=rawdata/dt/3600/1000; %[kW*m^2]

%%
pcolorjw(LON,LAT,rawdata_ssrs(1,:,:))


plot(lon, permute(rawdata_ssrs(1,45,:), [3 2 1]), lon, permute(rawdata_clear_sky(1,45,:), [3 2 1]), 'ro-', lon, permute(rawdata_ssrd(1,45,:), [3 2 1]), 'xg-')













