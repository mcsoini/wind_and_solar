%Calculate the power output for input wind speed arrays of arbitrary
%dimensionality

function [power]=speed2power(speed,pd)
power=zeros(size(speed));

parfor i=1:numel(speed)
    [i,numel(speed)];
    if  speed(i) <= 35
        %find next larger data point index to given speed
        i1=find(pd(:,1)>=speed(i),1,'first');
        %find next smaller data point index to given speed
        i0=find(pd(:,1)<=speed(i),1,'last');
        
        if i0==i1
            %exactly a data point
            power(i)=pd(i0,2);
        else
            %linear interpolation
            power(i)=(pd(i1,2)-pd(i0,2))/(pd(i1,1)-pd(i0,1))*(speed(i)-pd(i0,1))+pd(i0,2);
        end
    else
        power(i)=0;
    end
end

% figure(87)
% plot(speed(:),power(:),'o')
% xlabel('Wind speeds [m/s]')
% ylabel('Capacity factor [%]')
% hold on
% plot(pd(:,1),pd(:,2))
% drawnow
% hold off