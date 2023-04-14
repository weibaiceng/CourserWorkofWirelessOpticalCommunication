%�������
theta=10;
% semi-angle at half power �빦�ʽ�
m=-log10(2)/log10(cosd(theta));
%Lambertian order of emission �����ط������
P_total=20;
%transmitted optical power by individual LED ����LED���书��
Adet=1e-4;
%detector physical area of a PD �������й������ܵ��������

%���ջ�����
Ts=1;
%gain of an optical filter; ignore if no filter is used ��ѧ�˲���������
index=1.5;
%refractive index of a lens at a PD; ignore if no lens is used ������ƽ����͸��������
FOV=60*pi/180;
%FOV of a receiver ������ƽ���ӳ���
G_Con=(index^2)/sin(FOV);
%gain of an optical concentrator; ignore if no lens is used ��۹�������

%��������
lx=5; ly=5; lz=3;
% room dimension in metre ����ߴ�
h=2.15;
%the distance between source and receiver plane Դ�������ƽ������

%
%[XT,YT]=meshgrid([?1.25 1.25],[?1.25 1.25]); %���ɶ��led��������������
xt=[1 -1];
yt=[1 -1 ];
[XT,YT]=meshgrid(xt,yt);
%XT1=1; YT1=1;
%XT2=1; YT2=-1;
%XT3=-1; YT3=-1;
%XT4=-1; YT4=1;
% position of LED; LEDλ��
Nx=lx*10; Ny=ly*10;
% number of grid in the receiver plane ������ƽ���ϵ���������
x=-lx/2:lx/Nx:lx/2;
y=-ly/2:ly/Ny:ly/2;
[XR,YR]=meshgrid(x,y);
% receiver plane grid ������ƽ������


D1=sqrt((XR-XT(1,1)).^2+(YR-YT(1,1)).^2+h^2);
% distance vector from source 1 ��Դ��������ƽ����ÿ����ľ���
cosphi_A1=h./D1;
% angle vector
D2=sqrt((XR-XT(1,2)).^2+(YR-YT(1,2)).^2+h^2);
% distance vector from source 1 ��Դ��������ƽ����ÿ����ľ���
cosphi_A2=h./D2;
% angle vector
D3=sqrt((XR-XT(2,2)).^2+(YR-YT(2,2)).^2+h^2);
% distance vector from source 1 ��Դ��������ƽ����ÿ����ľ���
cosphi_A3=h./D3;
% angle vector
D4=sqrt((XR-XT(2,1)).^2+(YR-YT(2,1)).^2+h^2);
% distance vector from source 1 ��Դ��������ƽ����ÿ����ľ���
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