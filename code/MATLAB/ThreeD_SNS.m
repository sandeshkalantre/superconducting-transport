clear all
% Constants
q=1;
hbar=1;
%Temperature
kT = 0.00001;
%Fermi Level
mu_F =6e-3;
%Hopping Parameter in contacts
tx=0.002;
ty=0.002;
%t=tx+ty;
%Hopping Parameter in device
tdevx=tx;
tdevy=ty;
%tdev=tdevx+tdevy;
% Order Parameters
Delta1 =1e-3;
%phi = pi/2;
%Delta2 = Delta1*exp(1j*phi);
Nx=3;
Ny=3;
Nz=24; 
tz=1.7/900;
tdevz=tz;
t=tx+ty+tz;
EE=linspace(-1*Delta1,1*Delta1,4000);
%EE=9-mu_F:-0.005:3-mu_F;
%EE =linspace(0,10,50000);
zplus=(EE(2)-EE(1));
phase_count=0;
Flux_count=0;
%alpha1 = zeros(2, 2, 2, 2);
%alpha2 = zeros(2, 2, 2, 2);
% Flux
%Phi = linspace(0,8,100);
%Phi=0.000001;
Phi=0.0;
%phi=pi;
phi=1*linspace(0,2*pi,300);
%mu_F=0:0.1:6;
%phi=0.8*pi;
%T = zeros(length(EE),length(Phi));
DOS = zeros(length(EE),length(phi));
%DOS = zeros(length(EE),length(mu_F));
EXLz = zeros(length(EE),length(phi));
%I_E = zeros(length(EE),1);
for k=1:length(Nz)
    for i=1:length(phi)
        tic
            phase_count=phase_count+1;

            Flux_count=Flux_count+1;
            %Delta2=0;
             Delta2 = Delta1*exp(1j*phi(i));
            alpha1 = [2*(tx+abs(ty)+tz)+Phi^2-mu_F Delta1; Delta1' -2*(tx+abs(ty)+tz)-Phi^2+mu_F];
            alpha2 = [2*(tx+abs(ty)+tz)+Phi^2-mu_F Delta2; Delta2' -2*(tx+abs(ty)+tz)-Phi^2+mu_F];


    %% Density of States
    parfor j=1:length(EE)

               [DOS(j,phase_count,k),EXLz(j,phase_count,k)]=  Iop_3D(EE(j), tdevx, tdevy,tdevz, tx, ty,tz, zplus, Nx,Ny,Nz, mu_F, alpha1, alpha2,kT,Phi);
               %DOS(j,phase_count)=  Iop_3D(EE(j), tdevx, tdevy,tdevz, tx, ty,tz, zplus, Nx,Ny,Nz, mu_F, alpha1, alpha2,kT,Phi);
             %   100*j/length(EE)
           % j
    end
    (phi(i)/(2*pi))*100
    %k
    toc 
    %% Density of States Plots 
     phase_count=150;
            figure(10)
           % hold on
           [ax,hline1,hline2]= plotyy(EE/Delta1,real(EXLz(:,phase_count)),EE/Delta1,DOS(:,phase_count))
          % set(ax(1),'Color','b')
          %  set(ax(2),'Color','r')
           % yyaxis left
           %title('$\Phi$ = 0.1' ,'fontSize',40,'interpreter','latex')
             title(['$\Phi$ = 0.01, $\Delta \phi$ =  ' num2str(phi(1))],'fontSize',40,'interpreter','latex')
             ylabel(ax(1),'$\langle L_z \rangle$','interpreter','latex','fontsize',40)
             ylabel(ax(2),'Density Of States $(eV^{-1})$','interpreter','latex','fontsize',40)
             xlabel('$E/\Delta_0$','fontSize',40,'interpreter','latex')
             set(ax(1),'LineWidth',2,'fontSize',40)
             set(ax(2),'LineWidth',2,'fontSize',40)
             set(hline1,'LineWidth',4)
             set(hline2,'LineWidth',1)
             set(hline1,'Color','b')
             set(hline2,'Color','r')
    %         
    %        set(ax(1),'Color','r')
    %        set(ax(2),'Color','b')
            set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
            drawnow
    %         % pause
    %         figure (2)
    %         plot(EE,DOS(:,phase_count),'r', 'linewidth', 1)
    %        % yyaxis right
    %         title(['$\Delta \Phi$ =  ' num2str(phi(i))],'fontSize',40,'interpreter','latex')
    %         ylabel('Density Of States $(eV^{-1})$','interpreter','latex','fontsize',40)
    %         xlabel('$E/\Delta_0$','fontSize',40,'interpreter','latex')
    %         set(gca,'fontSize',40,'linewidth',2,'fontSize',40)

    %         figure (7)
    % %        % hold on
    %          plot(EE/Delta1,real(DOS(:,k)),'r', 'linewidth', 1)
    % %         %title('$\Phi$ = 0' ,'fontSize',40,'interpreter','latex')
    % %        % yyaxis right
    % %         title(['$\Phi$ = 0.01, \Delta \Phi$ =  ' num2str(phi(i))],'fontSize',40,'interpreter','latex')
    %          ylabel('$\langle L_z \rangle$','interpreter','latex','fontsize',40)
    %          xlabel('$E/\Delta_0$','fontSize',40,'interpreter','latex')
    %          set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
    %       
    %       drawnow
        %   pause
         % close
    %toc       
    %% CPR 
    % for j=1:length(EE)
    % %          I_E = @(E) Iop_3D(E, tdevx, tdevy,tdevz, tx, ty,tz, zplus, Nx,Ny,Nz, mu_F, alpha1, alpha2,kT,Phi);
    %           [I_E(j),EXLz(j)] = Iop_3D(EE(j), tdevx, tdevy,tdevz, tx, ty,tz, zplus, Nx,Ny,Nz, mu_F, alpha1, alpha2,kT,Phi);
    %          % I(phase_count) = integral(I_E, -5*Delta1, 0,'AbsTol',1e-6,'ArrayValued',true);
    % end
    %           phi(i)
    % % % % % 
          %  toc

     %% Transmission
    % for j=1:length(EE)
    %     
    %       T(j,i)= Iop_3D(EE(j), tdevx, tdevy,tdevz, tx, ty,tz, zplus, Nx,Ny,Nz, mu_F, alpha1, alpha2,kT,Phi(i));
    %           %  j
    % end
    % Flux_count
    %toc    
    end 
 %mu_F(k)
end
 %% DOS Plot
% 
%     DOSp = zeros(size(DOS));
%     for i=1:length(EE)
%         for j=1:length(phi)
%             if DOS(i,j) > 1.8e4
%                    DOSp(i,j)=1;
%             end
%         end
%     end    
% DOS1 = zeros(length(EE),length(phi));
% 
% for m=1:length(EE)
%     for n=1:length(phi)
%         DOS1(m,n)=DOS(2,m,n);
%     end
% end    
     figure(5)
    %surf(mu_F,EE,DOS,'edgecolor','none')
    % %hold on  
     surf((phi/(2*pi)),(EE/Delta1),real(DOS),'edgecolor','none')
     grid on
     az = 0;
     el = 90;
     view(az, el);
     title('$L$  $\approx 720$nm, $\Phi = 0.00$','fontSize',40,'interpreter','latex')
     xlabel('$\Delta \phi/2\pi$','interpreter','latex','fontsize',40)
     ylabel('$E/\Delta_0$','fontSize',40,'interpreter','latex')
     %xlabel('$\mu_F$','interpreter','latex','fontsize',40)
     %ylabel('$E$','fontSize',40,'interpreter','latex')
     set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
     map = [1,1,1; 0,0,1];
     colormap(map)
     colorbar
       caxis([0 100000])   
%% Transmission Plot
%     figure(3)
%     surf(Phi*Nx*Ny,(EE)/t,real(T),'edgecolor','none')
%     %hold on
%     %surf(phi_a/(2*pi),abs(EE/Delta1(1,1)),DOS_Eh)
%     az = 0;
%     el = 90;
%     view(az, el);
%     %title(['$\Delta \Phi$ =  ' num2str(phi)],'fontSize',40,'interpreter','latex')
%     xlabel('$\Phi/\Phi_0$','interpreter','latex','fontsize',40)
%     ylabel('$E/t$','fontSize',40,'interpreter','latex')
%     set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
%     colorbar
%     caxis([0 20])
%%  CPR Plots
% figure(3)
% hold on
% plot(EE/Delta1,real(I_E),'r','LineWidth',1) 
% %plot(phi/(2*pi),real(I)/Delta1,'g','LineWidth',2) 
% %hold on
% % plot(phi_a/(2*pi),real(I(2,:))/Delta1,'r','LineWidth',4) 
% % plot(phi_a/(2*pi),real(I(3,:))/Delta1,'k','LineWidth',4) 
% %title(['$\Delta \phi$ =  ' num2str(phi(i))],'fontSize',40,'interpreter','latex')
% title('Current Phase Relation','fontSize',40,'interpreter','latex')
% ylabel('Current $(e\Delta_0/\hbar)$','interpreter','latex','fontsize',40)
% xlabel('Energy E/$\Delta_0$','fontSize',40,'interpreter','latex')
% %xlabel('Phase Difference $\left(\Delta \phi\right)/2\pi$','fontSize',40,'interpreter','latex')
% %legend('t_y=t_x/10','t_y=t_x','interpreter','latex')
% set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
% 
%         figure(10)
%        % hold on
%        [ax,hline1,hline2]= plotyy(EE/Delta1,real(EXLz),EE/Delta1,real(I_E))
%       % set(ax(1),'Color','b')
%       %  set(ax(2),'Color','r')
%        % yyaxis left
%        %title('$\Phi$ = 0.1' ,'fontSize',40,'interpreter','latex')
%          title(['$\Phi$ = 0.0, $\Delta \phi$ =  ' num2str(phi(1))],'fontSize',40,'interpreter','latex')
%          ylabel(ax(1),'$\langle L_z \rangle$','interpreter','latex','fontsize',40)
%          ylabel(ax(2),'Density Of States $(eV^{-1})$','interpreter','latex','fontsize',40)
%          xlabel('$E/\Delta_0$','fontSize',40,'interpreter','latex')
%          set(ax(1),'LineWidth',2,'fontSize',40)
%          set(ax(2),'LineWidth',2,'fontSize',40)
%          set(hline1,'LineWidth',4)
%          set(hline2,'LineWidth',1)
%          set(hline1,'Color','b')
%          set(hline2,'Color','r')
% %         
% %        set(ax(1),'Color','r')
% %        set(ax(2),'Color','b')
%         set(gca,'fontSize',40,'linewidth',2,'fontSize',40)
%         drawnow
%%
% figure (3)
% for i=1:length(phi)
%     plot(EE,I_E(1:length(EE),i),'b')
%     title(['$\Delta \Phi$ =  ' num2str(phi(i))],'fontSize',40,'interpreter','latex')
%     pause
% end    
%%

save 24PointsSNS_30nmLC_B00_mu6meV.mat