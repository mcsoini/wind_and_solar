%%Extracts data from an ECMWF datafile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Mostly suited to extract average absolute wind speed data; returns a total 3D data matrix for a given pressure level in dependence on x,y, and t
% 
% Input:
% datafile:   filename string of the dataset to be read
% varname_1:  name of the variable 1 to be extracted from the the datafile
% varname_2:  name of the variable 2 to be extracted from the the datafile
% plevel:     pressure level to be extracted from the datafile (in general 1, since pressure level is set prior to downloading)
% dataname:   data type to be returned, only "wind_absolute" absolute value of wind speeds implemented
% t:          extract a specific time step or calculate the time average 'avg' to be returned as datamatrix
% shiftLon:   shift the longitudes; not reliable/obsolete -> download the data accordingly (-180...180)
%     
% Output:
% xv:         vector length m of longitudes
% yv:         vector length n of latitudes
% datamatrix: matrix mxn for selected pressure level and time step or time average 
% rawdata:    time resolved data for selected pressure
% timevec:    vector containing the hours for each time step as integers starting from 0 (?)
% level:      level name as stated in the data set; for double  check
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lon,lat,datamatrix_pv,rawdata_pv,datamatrix_csp,rawdata_csp,timevec]=getdata4Dm_solar(datafile)

%% Test
%datafile=strcat(filedir{11},filedir{12});

global output_dir_graph

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     DIRECT NORMAL IRRADIATION           %%%%%%%%%%%%%%%%%%%%%

    varname_1='Surface_solar_radiation_downwards_surface';
    dT=12; %According to "Select time"
    dt=3;

    '... entered getdata4Dm_solar ...'

    % create ncgeodataset object, contains all the data from the file
    nc=ncgeodataset(datafile);

    nctime=nc.geovariable('time');
    timevec=nctime.data(:);

    % create geovariable objects
    dirvar=nc.geovariable(varname_1);

    % get the grid for time period 1
    g=dirvar.grid_interop(1,:,:);

    % extract the vectors containing the longitudes and latitudes
    lon=g.lon;
    lat=g.lat;   
    [LON,LAT]=meshgrid(lon,lat);

    rawdata=zeros([length(timevec), size(LON)]);
    
    rawdata=dirvar.data(:,:,:);  %in [W*m^-2*s/Xh]
    
%     'Copying data from dataset to matlab variable'
%     %Loop necessary due to java heap space limitations
%     parfor_progress(length(timevec));
%     for i=1:length(timevec)
%         rawdata(i,:,:)=dirvar.data(i,:,:);  %in [W*m^-2*s/Xh]
%     parfor_progress;
%     end
%     parfor_progress(0);
    
%     figure(3409657)
%     for i=1230:length(rawdata)
%         mod(timevec(i),24)
%         pcolorjw(lon,lat,rawdata(i,:,:));
%         hold on
%         %plot(xbglob,ybglob)
%         hold off
%         drawnow
%         pause
%     end

%%%%%%%%%%%%%%%%   CONVERTING THE TEMPORALLY SEMI-CUMULATIVE DATA TO  [kWs*m^-2/Xh] %%%%%%%    
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
rawdata=rawdata/dt/3600/1000; %[kW*m^2]

%%%%%%%%%%%%%     DIRECT NORMAL IRRADIATION           %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

'Copying rawdata to PV and CSP'
%CSP assumes 2D tracking of the sun, therefore equal to the direct normal
%irradiation
rawdata_csp=rawdata;
%PV assumes no tracking -> calculate the component of the radiation normal
%to the panel surface
rawdata_pv=zeros(size(rawdata));

clear rawdata;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%     SOLAR PV CALCULATION                %%%%%%%%%%%%%%%%%%%%%

%define the assumed angle of the solar panels as a function of latitude
%first approximation: equal the latitude
%ptilt=180*ones(size(lat)); %flat on the ground (for test purposes)
ptilt=-lat+180; %equal the latitude

%vector components of the normal vector on the panel for all latitudes in a
%coordinate system South(x)...East(y)...Zenith(z)
np_S = sin((ptilt)*pi/180);
np_Z = -cos((ptilt)*pi/180);
np_E = zeros(size(ptilt));

NP_S=zeros(size(LON));
NP_Z=zeros(size(LON));
NP_E=zeros(size(LON));
for i=1:length(lon)
    NP_S(:,i)=np_S;
    NP_Z(:,i)=np_Z;
end
NP=permute(cat(3,NP_S,NP_E,NP_Z),[3 1 2]);
size(NP)

%hour offset: alignment offset between ECMWF data and solar position
%calculation; not sure why this is necessary (all data and values should
%be in UTC)

timevec=double(timevec)

h_offs=(timevec(end-1)-timevec(end))*0.5;

timevec=timevec+h_offs

datenum_base=datenum('2013/01/01 00:00:00');
'calculating normal PV irradiation...'
parfor_progress(length(timevec));
for i=1:length(timevec)
    
    %get the point in time where to calculate the sun's position
    [h,m,s] = decimalhours(timevec(i));
    date_number=addtodate(datenum_base,h,'hour');
    date_number=addtodate(date_number,m,'minute');
    date_number=addtodate(date_number,s,'second');
    
    %Calculate the coordinates of the sun's position in the coordinate
    %system South-Zenith-East at each location
    [~,~,~,xhor,yhor,zhor]=SolarAzEl(date_number,LAT,LON,0);

    %OBSOLETE
    %vector components of the sun's position in a
    %coordinate system South(x)...East(y)...Zenith(z)
    %sun_S=-cos(Az*pi/180).*cos((El)*pi/180);
    %sun_E= sin(Az*pi/180).*cos((El)*pi/180);
    %sun_Z= sin((El)*pi/180);
    
    %r is the fraction normal component for all coordinates
    RHOR=permute(cat(3,xhor, yhor,zhor),[3 1 2]);
    r=permute(dot(RHOR,NP),[2,3,1]);
    r(r<0)=0;
    
    %for the current time step, multiply the fraction of the panel normal
    %irradiation and the direct normal irradiation
    rawdata_pv(i,:,:)=permute(rawdata_csp(i,:,:),[2 3 1]).*r;
    
% %%
% h=figure(5687)
% subplot(2,1,1)
%     pcolorjw(lon,lat,rawdata_csp(i,:,:))
%     h=colorbar;
%     title(h,'')
%     hold on
%     plot(xbglob,ybglob)
%     hold off
%     date_number=addtodate(datenum_base,timevec(i)+h_offs,'hour');
%     title(datestr(date_number))
% subplot(2,1,2)
%     pcolorjw(lon,lat,r)
%     h=colorbar;
%     title(h,'')
%     hold on
%     plot(xbglob,ybglob)
%     hold off
%     date_number=addtodate(datenum_base,timevec(i)+h_offs,'hour');
%     title(datestr(date_number))
%     drawnow
%     
%     
%     max(max(rawdata_pv(i,:,:)))
%     pause
    
    
    parfor_progress;
end
parfor_progress(0);
'... done'
%%%%%%%%%%%%%     SOLAR PV CALCULATION                %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%%%%%%%%%%%%%%%%   CALCULATING THE TIME AVERAGE OF THE DATA  %%%%%%%%%%%%%%    
%calculate the time average of the PV matrix
datamatrix_pv =sum(rawdata_pv )/length(timevec)*8760; %in [kWh/yr/m_panel^2]
datamatrix_csp=sum(rawdata_csp)/length(timevec)*8760; %in [kWh/yr/m_panel^2]

%permute the dimensions of the output matrix to get rid of the time
%dimension (which contains only one step at this point)
datamatrix_pv =permute(datamatrix_pv , [ 2 3 1 ]);
datamatrix_csp=permute(datamatrix_csp, [ 2 3 1 ]);

nfig=2355;
f=figure(nfig)
%screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );
subplot(2,1,1)
    pcolorjw(lon,lat,datamatrix_pv)
    h=colorbar;
    title(h,'kWh/yr PV')
    xlabel('Longitude')
    ylabel('Latitude')
subplot(2,1,2)
    pcolorjw(lon,lat,datamatrix_csp)
    h=colorbar;
    title(h,'kWh/yr CSP')
    xlabel('Longitude')
    ylabel('Latitude')

%tightfig
%screen_size=get(0,'ScreenSize');set(nfig, 'units','normalized','position', [0 0 screen_size(3) screen_size(4) ] );

%print(f,'-dpng',strcat(output_dir_graph, 'AA_solar_data.png'))
    
    
'... leaving getdata4Dm_solar ...'


