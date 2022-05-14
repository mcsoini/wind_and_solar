%% The following picks the windiest pixels of the region successively to 
% obtain a maximum capacity of Cmax in the region. Thereby, the capacity density in W*km^-2
% is given by capdens. The result is a list of indices indexarray_hispeeds

function [indexarray_hispeeds,datamat]=highest_wind_speeds_capacity_limited(x,y,datamat,Cmax,capdens)
% figure(2)
% pcolorjw(x,y,datamat1)
% drawnow
% hold on
%%
i=0;
Csum=0;
indexarray_hispeeds=[];
while Csum<=Cmax%add capacities until the maximum is reached
    %find value and index of the maximum average wind speed
    [max_speed,ind_max_speed]=max(datamat(:))
    
    %convert the index to 2D subscripts
    [y_ind_max_speed,~]=ind2sub(size(datamat),ind_max_speed);
    
    if length(indexarray_hispeeds)~=0 && ind_max_speed==indexarray_hispeeds(end)%something awry going on
        'error'
        break
    end
    %store the index in the index vector for access to the chosen pixels
    indexarray_hispeeds=[indexarray_hispeeds;[ind_max_speed,max_speed]];
    
    %calculate the area of the pixel with the maximum wind speed
    A_max_speed=110^2*(abs(x(end)-x(end-1)))^2*cos(pi/180*y(y_ind_max_speed));
    
    %add the pixel's contribution to the total capacity of the region
    Csum=Csum+A_max_speed*capdens*8760*3600; %km²*J/(s*km²)*h/yr*s/h
    
    %set the pixel's average wind speed zero in the copied matrix datamat1
    datamat(ind_max_speed)=0;
end