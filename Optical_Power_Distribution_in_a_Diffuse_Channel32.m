clear;
P_total=1;
% Total transmitted power
% rho=0.8;
rho1=0.8;
rho2=0.8;
rho3=0.8;
rho4=0.8;
%reflection coefficient
lx=5; ly=5; lz=2.15;
% room dimension in metre
FOV =60; %FOV angle
m = 1;% Lambert's order
Adet=1e-4; % active area
Nx=lx*10; Ny=ly*10; Nz=round(lz*10);
% number of grid in each surface
Ts=1;
%gain of an optical filter; ignore if no filter is used
index=1;
%refractive index of a lens at a PD; ignore if no lens is used
G_Con=(index^2)/sin(FOV)^2;
%gain of an optical concentrator; ignore if no lens is used
dA=lz*ly/(Ny*Nz);
% calculation grid area
x=linspace(-lx/2,lx/2,Nx);
y=linspace(-ly/2,ly/2,Ny);
z=linspace(-lz/2,lz/2,Nz);
[XR,YR,ZR]=meshgrid(x,y,-lz/2);
TP1=[0 0 lz/2];
% transmitter position
%%%%%%%%%%%%%%%calculation for 4 walls %%%%%%%%%%%%%%%%%%
for ii=1:Nx
    for jj=1:Ny
    RP=[x(ii) y(jj) -lz/2];
    % receiver position vector
    h1(ii,jj)=0;
    % reflection from wall 1
        for kk=1:Ny
            for ll=1:Nz
            WP1=[-lx/2 y(kk) z(ll)];
            % point of incidence in wall
            D1=sqrt(dot(TP1-WP1,TP1-WP1));
            % distance from transmitter to WP1
            cos_phi=abs(WP1(3)-TP1(3))/D1;
            cos_alpha=abs(TP1(1)-WP1(1))/D1;
            D2=sqrt(dot(WP1-RP,WP1-RP));
            % distance from WP1 to receiver
            cos_beta=abs(WP1(1)-RP(1))/D2;
            cos_psi=abs(WP1(3)-RP(3))/D2;
                if abs(acosd(cos_psi))<=FOV
                h1(ii,jj)=h1(ii,jj)+(m+1)*Adet*rho1*dA*...
                cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end    
       h2(ii,jj)=0;
    % reflection from wall 2
        for kk=1:Ny
            for ll=1:Nz
            WP1=[lx/2 y(kk) z(ll)];
            % point of incidence in wall
            D1=sqrt(dot(TP1-WP1,TP1-WP1));
            % distance from transmitter to WP1
            cos_phi=abs(WP1(3)-TP1(3))/D1;
            cos_alpha=abs(TP1(1)-WP1(1))/D1;
            D2=sqrt(dot(WP1-RP,WP1-RP));
            % distance from WP1 to receiver
            cos_beta=abs(WP1(1)-RP(1))/D2;
            cos_psi=abs(WP1(3)-RP(3))/D2;
                if abs(acosd(cos_psi))<=FOV
                h2(ii,jj)=h2(ii,jj)+(m+1)*Adet*rho2*dA*...
                cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
        h3(ii,jj)=0;
    % reflection from wall 3
        for kk=1:Nx
            for ll=1:Nz
            WP1=[x(kk) -ly/2  z(ll)];
            % point of incidence in wall
            D1=sqrt(dot(TP1-WP1,TP1-WP1));
            % distance from transmitter to WP1
            cos_phi=abs(WP1(3)-TP1(3))/D1;
            cos_alpha=abs(TP1(1)-WP1(1))/D1;
            D2=sqrt(dot(WP1-RP,WP1-RP));
            % distance from WP1 to receiver
            cos_beta=abs(WP1(1)-RP(1))/D2;
            cos_psi=abs(WP1(3)-RP(3))/D2;
                if abs(acosd(cos_psi))<=FOV
                h3(ii,jj)=h3(ii,jj)+(m+1)*Adet*rho3*dA*...
                cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
        h4(ii,jj)=0;
    % reflection from wall 4
        for kk=1:Nx
            for ll=1:Nz
            WP1=[x(kk) ly/2  z(ll)];
            % point of incidence in wall
            D1=sqrt(dot(TP1-WP1,TP1-WP1));
            % distance from transmitter to WP1
            cos_phi=abs(WP1(3)-TP1(3))/D1;
            cos_alpha=abs(TP1(1)-WP1(1))/D1;
            D2=sqrt(dot(WP1-RP,WP1-RP));
            % distance from WP1 to receiver
            cos_beta=abs(WP1(1)-RP(1))/D2;
            cos_psi=abs(WP1(3)-RP(3))/D2;
                if abs(acosd(cos_psi))<=FOV
                h4(ii,jj)=h4(ii,jj)+(m+1)*Adet*rho4*dA*...
                cos_phi^m*cos_alpha*cos_beta*cos_psi/(2*pi^2*D1^2*D2^2);
                end
            end
        end
    end
end
% calculate channel gain (h1，h2, h3 and h4) from 4 walls
P_rec_A1=(h1+h2+h3+h4)*P_total.*Ts.*G_Con;
P_rec_dBm=10*log10(P_rec_A1);
meshc(x,y,P_rec_dBm);
xlabel('X(m)');
ylabel('Y (m)');
zlabel('Received power (dBm)');