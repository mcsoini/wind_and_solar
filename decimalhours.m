function [h,m,s] = decimalhours(hd)

hd_res=hd;

h=0;
m=0;
s=0;

while hd_res>0
    if hd_res>1
        hd_res=hd_res-1;
        h=h+1;
    elseif 60>1/hd_res
        hd_res=hd_res-1/60;
        m=m+1;
    elseif (60*60)>1/hd_res
        hd_res=hd_res-1/(60*60);
        s=s+1;
    else
        hd_res=-1;
    end
end