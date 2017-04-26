filename = 'mp_bmp_20170106';

for i = 1:12
    load MP_BMP_20170106.mat
    a = find(ID==i);
    angle_deg_N0=angle_deg_N0(a);
    angle_deg_E0=angle_deg_E0(a);
    angle_rad=angle_rad(a);
    dist=dist(a);
    lat=lat(a);
    leg_time=leg_time(a);
    lon=lon(a);
    range=range(a);
    speed=speed(a);
    tide_time=tide_time(a);
    time=time(a);
    Tr=Tr(a);
    u=u(a);
    v=v(a);
    WL=WL(a);
    save(strcat(filename,'_',num2str(i),'.mat'),'angle_deg_N0','angle_deg_E0','angle_rad','dist','lat','leg_time','lon','range','speed','tide_time','time','Tr','u','v','WL')
end