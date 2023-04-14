%物理参数
theta=10;
% semi-angle at half power 半功率角
m=-log10(2)/log10(cosd(theta));
%Lambertian order of emission 兰伯特发射阶数
P_total=20;
%transmitted optical power by individual LED 单个LED发射功率
Adet=1e-4;
%detector physical area of a PD 接收器中光电二极管的物理面积

%接收机参数
Ts=1;
%gain of an optical filter; ignore if no filter is used 光学滤波器的增益
index=1.5;
%refractive index of a lens at a PD; ignore if no lens is used 接收器平面上透镜折射率
FOV=60*pi/180;
%FOV of a receiver 接收器平面视场角
G_Con=(index^2)/sin(FOV);
%gain of an optical concentrator; ignore if no lens is used 光聚光器增益

%环境参数
lx=5; ly=5; lz=3;
% room dimension in metre 房间尺寸
h=2.15;
%the distance between source and receiver plane 源与接收器平面间距离

%
%[XT,YT]=meshgrid([?1.25 1.25],[?1.25 1.25]); %生成多个led（网格点坐标矩阵）
xt=[1 -1];
yt=[1 -1 ];
[XT,YT]=meshgrid(xt,yt);
%XT1=1; YT1=1;
%XT2=1; YT2=-1;
%XT3=-1; YT3=-1;
%XT4=-1; YT4=1;
% position of LED; LED位置
Nx=lx*10; Ny=ly*10;
% number of grid in the receiver plane 接收器平面上的网格数量
x=-lx/2:lx/Nx:lx/2;
y=-ly/2:ly/Ny:ly/2;
[XR,YR]=meshgrid(x,y);
% receiver plane grid 接收器平面坐标


D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2);
% distance vector from source 1 从源到接收器平面上每个点的距离
cosphi_A1=h./D1;
% angle vector
D2=sqrt((XR-XT(1,2)).^2+(YR-YT(1,2)).^2+h^2);
% distance vector from source 1 从源到接收器平面上每个点的距离
cosphi_A2=h./D2;
% angle vector
D3=sqrt((XR-XT(2,2)).^2+(YR-YT(2,2)).^2+h^2);
% distance vector from source 1 从源到接收器平面上每个点的距离
cosphi_A3=h./D3;
% angle vector
D4=sqrt((XR-XT(2,1)).^2+(YR-YT(2,1)).^2+h^2);
% distance vector from source 1 从源到接收器平面上每个点的距离
cosphi_A4=h./D4;
% angle vector


H_A1=(m+1)*Adet.*cosphi_A1.^(m+1)./(2*pi.*D1.^2);
% channel DC gain for source 1
P_rec1=P_total.*H_A1.*Ts.*G_Con;
% received power from source 1;
H_A2=(m+1)*Adet.*cosphi_A2.^(m+1)./(2*pi.*D2.^2);
% channel DC gain for source 2
P_rec2=P_total.*H_A2.*Ts.*G_Con;
% received power from source 2;
H_A3=(m+1)*Adet.*cosphi_A3.^(m+1)./(2*pi.*D3.^2);
% channel DC gain for source 3
P_rec3=P_total.*H_A3.*Ts.*G_Con;
% received power from source 3;
H_A4=(m+1)*Adet.*cosphi_A4.^(m+1)./(2*pi.*D4.^2);
% channel DC gain for source 4
P_rec4=P_total.*H_A4.*Ts.*G_Con;
% received power from source 4;
P_rec=P_rec1 + P_rec2 + P_rec3 + P_rec4;
P_rec_dBm=10*log10(P_rec);

meshc(x,y,P_rec_dBm);
xlabel('X(m)');
ylabel('Y(m)');
zlabel('Received power (dBm)');
axis([-lx/2 lx/2 -ly/2 ly/2 min(min(P_rec_dBm)) max(max(P_rec_dBm))]);
% axis([0 lx 0 ly min(min(P_rec_dBm)) max(max(P_rec_dBm))]);