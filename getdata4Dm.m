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

function [xv,yv,datamatrix,rawdata,timevec,level]=getdata4Dm(datafile,varname_1,varname_2,plevel,dataname,t,shiftLon)
%%
%datafile=strcat(filedir{1},filedir{2});

global output_dir_graph

'... entered getdata4Dm ...'

% create ncgeodataset object, contains all the data from the file
nc=ncgeodataset(datafile);

nctime=nc.geovariable('time');
timevec=nctime.data(:);

if strcmp(dataname,'windabsolute')

    % create geovariable objects
    dirvarU=nc.geovariable(varname_1);
    dirvarV=nc.geovariable(varname_2);
    
    levelvec=nc.geovariable('hybrid');
    level=levelvec.data(plevel);
    
    %get the grid for time period 1
    g=dirvarU.grid_interop(1,1,:,:);
        
    %extract the vectors containing the longitudes and latitudes
    xv=g.lon;
    yv=g.lat;   
    
    %%OBSOLETE
    %get the index corresponding to the longitudinal shift
    [~,nshift]=min(abs(xv+shiftLon));

    %shift the longitudinal vector by shiftLon
    if shiftLon~=0
        xv=[xv(nshift:end,1)-360;xv(1:nshift-1,1)];
    end

    % Extract a given pressure level and shift the other dimensions accordingly
    dirvarUpconst=permute(dirvarU.data(:,plevel,:,:),[1 3 4 2]);
    dirvarVpconst=permute(dirvarV.data(:,plevel,:,:),[1 3 4 2]);
    
%%%%%%%%%%%%%%%%   CALCULATING THE ABSOLUTE WIND SPEEDS    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%   FROM THE U AND V COMPONENTS             %%%%%%%%%%%%%%%%
    tic; 'absolute wind speed...'
    %direct version... handles really large data structures eventually,
    %maybe too large for the java heap space... use the memory limited
    %version if it is not possible to use this one
    rawdata(:,:,:)=sqrt(dirvarUpconst(:,:,:).^2+dirvarVpconst(:,:,:).^2);

    %loop version: works even if the 3D matrix based approach is not
    %possible. The additional definition of the cellarray timestr (date and
    %time as string for each time step) seems to be obsolete
%     parfor i=1:tN
%         %g=dirvarU.grid_interop(i,:,:);
%         %timestr{i}=datestr(g.time);
%         rawdata(i,:,:)=sqrt(dirvarUpconst(i,:,:).^2+dirvarVpconst(i,:,:).^2);
%         %[i,tN]%,timestr{i}]
%     end
    '...absolute wind speed done.'
    toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%%%%%%%%%   CALCULATING THE TIME AVERAGE OF THE WIND SPEEDS  %%%%%%%    
    if strcmp(t,'avg')
        %calculate the time average
        'Calculating average wind speeds...'
        tic
        datamatrix=sum(rawdata)/length(rawdata(:,1,1));
        toc
        'ok'
    else %if not time average but a specific time step t
        %data for time step t
        datamatrix=rawdata(t,:,:);
        g=dirvarU.grid_interop(t,:,:);    
        date=datestr(g.time)
    end
    
    %permute the dimensions of the output matrix to get rid of the time
    %dimension (which contains only one step at this point)
    datamatrix=permute(datamatrix, [ 2 3 1 ]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
%%%%%%%%  SHIFT THE MATRIX TO MAKE THE LEFT BORDER -180W            %%%%%%%
%%%%%%%% OBSOLETE: RATHER DOWNLOAD THE DATA ACCORDINGLY  %%%%%%%%%%%%%%%%%%
    if shiftLon~=0
        datamatrix=[datamatrix(:,nshift:end),datamatrix(:,1:nshift-1)];
        for i=1:length(rawdata(:,1,1))
            rawdataipermute=permute(rawdata(i,:,:),[ 2 3 1]);
            rawdataipermute=[rawdataipermute(:,nshift:end),rawdataipermute(:,1:nshift-1)];
            rawdata(i,:,:)=rawdataipermute;
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    

    nfig=3636;
    f=figure(nfig)
    pcolorjw(xv,yv,datamatrix)
    xlabel('Longitude [deg]')
    ylabel('Latgitude [deg]')
    h = colorbar;
    ylabel(h,'Wind speed [m/s]')
    
 %   print(f,'-dpng',strcat(output_dir_graph, 'AA_wind_data.png'))

%%use elseif to access other variables depending on input 'dataname'
% elseif strcmp(dataname,'10_metre_U_wind_component_surface')
% 
%     % create geovariable object
%     dirvarU=nc.geovariable('10_metre_U_wind_component_surface');
% 
%     % get data at time step t
%     dirU=dirvarU.data(t,:,:);
% 
%     % get grid at time step t
%     g=dirvarU.grid_interop(t,:,:);    
%     
%     % define output variables
%     datamatrix=dirU;
%     xvector=g.lon;
%     yvector=g.lat;

else
    datamatrix=[];
    xv=[];
    yv=[];
end

'... leaving getdata4Dm ...'


% 
% % plot
% pcolorjw(g.lon,g.lat,dir);
% %quiver(glon,glat,permute(dirU, [2 3 1]),permute(dirV, [2 3 1]));
% title(datestr(g.time))
% drawnow
% pause(0.01)
