function [DOS,EXLz] = Iop_3D(E, tdevx, tdevy,tdevz, tx, ty,tz, zplus, Nx,Ny,Nz, mu_F, alpha1, alpha2,kT,Phi)
q=1;
hbar=1;
%DOS = zeros(length(E));
%EXLz = zeros(length(E));
%% Constant Longitudinal Magnetic Field
beta_t1 = [-tx 0; 0 conj(tx)];
beta_t2 = [-ty 0; 0 conj(ty)];
beta_l = [-tz 0; 0 tz];
beta_dev_t1 = [-tdevx 0; 0 conj(tdevx)];
beta_dev_t2 = [-tdevy 0; 0 conj(tdevy)];
beta_dev_l = [-tdevz 0; 0 tdevz];

tdev=tdevx+tdevy+tdevz;
t=tx+ty+tz;

for countx = 1:Nx
    ty = ty*exp(1i*2*pi*countx*Phi/(Nx*Ny));
    beta_t2 = [-ty 0; 0 conj(ty)];
    alpha1_column_y=kron(eye(Ny),alpha1)+kron(diag(ones(1,Ny-1),+1),beta_t2)+kron(diag(ones(1,Ny-1),-1),beta_t2');
    alpha2_column_y=kron(eye(Ny),alpha2)+kron(diag(ones(1,Ny-1),+1),beta_t2)+kron(diag(ones(1,Ny-1),-1),beta_t2');
    
    if countx == 1
        alpha1_x = alpha1_column_y;
        alpha2_x = alpha2_column_y;
   
    else 
        alpha1_x = blkdiag(alpha1_x,alpha1_column_y);
        alpha2_x = blkdiag(alpha2_x,alpha2_column_y);

    end
   
end    
    beta_t1_column_x = kron(eye(Ny),beta_t1);
    alpha1_crosssection=alpha1_x+kron(diag(ones(1,Nx-1),+1),beta_t1_column_x)+kron(diag(ones(1,Nx-1),-1),beta_t1_column_x');

    alpha2_crosssection=alpha2_x+kron(diag(ones(1,Nx-1),+1),beta_t1_column_x)+kron(diag(ones(1,Nx-1),-1),beta_t1_column_x');

    beta_contact_column_y = kron(diag(exp(0*[1:1:Ny])),beta_l);
    beta_contact_crosssection=kron(diag(exp(0*[1:1:Nx])),beta_contact_column_y);



%% No Magnetic Field
% beta_t1 = [-tx 0; 0 tx];
% beta_t2 = [-ty 0; 0 ty];
% beta_l = [-tz 0; 0 tz];
% beta_dev_t1 = [-tdevx 0; 0 tdevx];
% beta_dev_t2 = [-tdevy 0; 0 tdevy];
% beta_dev_l = [-tdevz 0; 0 tdevz];
% 
% tdev=tdevx+tdevy+tdevz;
% t=tx+ty+tz;
% 
% 
%     alpha1_column_y=kron(eye(Ny),alpha1)+kron(diag(ones(1,Ny-1),+1),beta_t2)+kron(diag(ones(1,Ny-1),-1),beta_t2');
%     beta_t1_column_x = kron(eye(Ny),beta_t1);
%     alpha1_crosssection=kron(eye(Nx),alpha1_column_y)+kron(diag(ones(1,Nx-1),+1),beta_t1_column_x)+kron(diag(ones(1,Nx-1),-1),beta_t1_column_x');
% 
%     alpha2_column_y=kron(eye(Ny),alpha2)+kron(diag(ones(1,Ny-1),+1),beta_t2)+kron(diag(ones(1,Ny-1),-1),beta_t2');
%     alpha2_crosssection=kron(eye(Nx),alpha2_column_y)+kron(diag(ones(1,Nx-1),+1),beta_t1_column_x)+kron(diag(ones(1,Nx-1),-1),beta_t1_column_x');
% 
%     beta_contact_column_y = kron(diag(exp(0*[1:1:Ny])),beta_l);
%     beta_contact_crosssection=kron(diag(exp(0*[1:1:Nx])),beta_contact_column_y);

%% Alpha Device NO Magnetic Field
%     alpha_2x2 = [2*tdevx + 2*tdevy + 2*tdevz-mu_F 0; 0 -2*tdevx - 2*tdevy - 2*tdevz+mu_F];
%     alpha_column_y=kron(eye(Ny),alpha_2x2)+kron(diag(ones(1,Ny-1),+1),beta_dev_t2)+kron(diag(ones(1,Ny-1),-1),beta_dev_t2');
%     beta_dev_t1_column_x = kron(eye(Ny),beta_dev_t1);
%     alpha=kron(eye(Nx),alpha_column_y)+kron(diag(ones(1,Nx-1),+1),beta_dev_t1_column_x)+kron(diag(ones(1,Nx-1),-1),beta_dev_t1_column_x');

%% Alpha Device Constant Longitudinal Magnetic Field
    alpha_2x2 = [2*tdevx + 2*abs(tdevy) + 2*tdevz+Phi^2-mu_F 0; 0 -2*tdevx - 2*abs(tdevy) - 2*tdevz-Phi^2+mu_F];
    beta_dev_t1_column_x = kron(eye(Ny),beta_dev_t1);
   
    for countx = 1:Nx
    tdevy = tdevy*exp(1i*2*pi*countx*Phi/(Nx*Ny));
    beta_dev_t2 = [-tdevy 0; 0 conj(tdevy)];
    alpha_column_y=kron(eye(Ny),alpha_2x2)+kron(diag(ones(1,Ny-1),+1),beta_dev_t2)+kron(diag(ones(1,Ny-1),-1),beta_dev_t2');
 
    if countx == 1
        alpha_x = alpha_column_y;
   
    else 
        alpha_x = blkdiag(alpha_x,alpha_column_y);

    end
   
    end   
    
    alpha=alpha_x+kron(diag(ones(1,Nx-1),+1),beta_dev_t1_column_x)+kron(diag(ones(1,Nx-1),-1),beta_dev_t1_column_x');
%%
    beta_column_y = kron(diag(exp(0*[1:1:Ny])),beta_dev_l);
    beta=kron(diag(exp(0*[1:1:Nx])),beta_column_y);
    
    H=kron(eye(Nz),alpha);
%     H(1:2*Nx*Ny,1:2*Nx*Ny)=alpha1_crosssection;
%     H(2*Nx*Ny*Nz-2*Nx*Ny+1:2*Ny*Nx*Nz,2*Nx*Ny*Nz-2*Nx*Ny+1:2*Ny*Nx*Nz) = alpha2_crosssection;
    if Nz>1
        H=H+kron(diag(ones(1,Nz-1),+1),beta)+kron(diag(ones(1,Nz-1),-1),beta');
    end
    
    Cy=eye(Ny); Cy(1,1)=0; Cy(Ny,Ny)=0;
    Cx=eye(Nx); Cx(1,1)=0; Cx(Nx,Nx)=0;
    
%     Ax = zeros(Ny,Ny);
%     Ax(1,1) = 1; Ax(3,3)=1; Ax(Ny,Ny) = 1; Ax(Ny-2,Ny-2) = 1;
%     
    
    
    U = zeros(2*Nx*Ny*Nz,2*Nx*Ny*Nz); % confinement potential
    Uconf = [10000*mu_F 0; 0 -10000*mu_F];
    U_column_y = kron(Cy,Uconf);
    U = kron(eye(Nz),kron(Cx,U_column_y));
%% Angular Momentum Operator
%% X
x=zeros(2*Nx,1);
for xc = 1:Nx
    x(2*xc-1)=xc;
    x(2*xc)=xc;
end
X = kron(eye(Nz),kron(diag(x),eye(Ny)));
%% Y
y=zeros(2*Ny,1);
for yc = 1:Ny
    y(2*yc-1)=yc;
    y(2*yc)=yc;
end
Y = kron(eye(Nz),kron(eye(Nx),diag(y)));
%% px
px = -1i*(X*H-H*X);
% px_x = kron(eye(Nx),kron(eye(Ny),d_diag)) + kron(diag(ones(1,Nx*Ny-Ny),-Ny),d_diagm1);
% px = -1i*kron(eye(Nz),px_x);

%% py
py = -1i*(Y*H-H*Y);
% py_xy = kron(eye(Ny),d_diag) + kron(diag(ones(1,Ny-1),-1),d_diagm1);
% py_x = kron(eye(Nx),py_xy);
% py = -1i*kron(eye(Nz),py_x);



Lz = X*py-Y*px;
       

%% Surgace Green's Functions

        ig1=(E+1j*zplus).*eye(2*Ny*Nx)-alpha1_crosssection;
        ig2=(E+1j*zplus).*eye(2*Ny*Nx)-alpha2_crosssection;
       
        gs1=inv(ig1);gs2=inv(ig2);
    
        change=1;
        %Recursive calculation of surface Green's function
     %   tic
        while change >1e-6
        Gs=inv(ig1-beta_contact_crosssection*gs1*beta_contact_crosssection');
        change=sum(sum(abs(Gs-gs1)))/(sum(sum(abs(gs1)+abs(Gs))));
        gs1=0.5*Gs+0.5*gs1;
        %change
        %gs1
        end
      %  toc
        change=1;
     %   tic
        while change >1e-6
        Gs=inv(ig2-beta_contact_crosssection'*gs2*beta_contact_crosssection);
        change=sum(sum(abs(Gs-gs2)))/(sum(sum(abs(gs2)+abs(Gs))));
        gs2=0.5*Gs+0.5*gs2;
        %change
        %gs2
        end
     %  toc
    sig1=beta_contact_crosssection'*gs1*beta_contact_crosssection;gam1=1i*(sig1-sig1');
    sig2=beta_contact_crosssection*gs2*beta_contact_crosssection';gam2=1i*(sig2-sig2');

    Sigma1 = zeros(2*Nx*Ny*Nz,2*Nx*Ny*Nz); Sigma2 = zeros(2*Nx*Ny*Nz,2*Nx*Ny*Nz);
    Sigma1(1:2*Nx*Ny,1:2*Nx*Ny)=sig1;
    Sigma2(2*Nx*Ny*Nz-2*Nx*Ny+1:2*Ny*Nx*Nz,2*Nx*Ny*Nz-2*Nx*Ny+1:2*Ny*Nx*Nz)=sig2;

    Gamma1 = 1j*(Sigma1 - Sigma1');
    Gamma2 = 1j*(Sigma2 - Sigma2');


%%
    %GD =  inv((E+1j*zplus).*eye(2*Nx*Ny*Nz) - H -U);
   GD =  inv((E+1j*zplus).*eye(2*Nx*Ny*Nz) - H -U-Sigma1-Sigma2);
    %    T = real(trace(Gamma1*GD*Gamma2*GD'));
    A = 1j*(GD-GD');
    DOS = trace(A);
        
    fermi = 1.0/(1.0 + exp(E)/kT);
    
    Fermi_matrix = [fermi 0;0 fermi];
    Fermi = kron(eye(Nx*Ny*Nz),Fermi_matrix);
%     
%     Sigma_corr = (Gamma1+ Gamma2)*Fermi;
%     G_corr = GD*Sigma_corr*GD';
     EXLz = trace(A*Lz)/trace(A);
   % EXLz = trace(A*Lz)/trace(A);
%     %DOS = trace(1i*(GD-GD'));
%       I_op = (1j*q/hbar)*(H*G_corr - G_corr*H);
%       I=0;
%       for i=1:2:2*Nx*Ny-1
%           I = I+I_op(i,i)-I_op(i+1,i+1);
%       end
end