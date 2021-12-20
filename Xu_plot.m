%% clear and clc;

close all; clear; clc;
addpath('chebfun-master');

%% read and processing data
filename = 'PhipsData_20210629-0658_level_3';

PhipsData = readtable(filename);
PhipsDataMatrix = readmatrix(filename);
vbs = PhipsData.Properties.VariableNames;
Index_scattering = find(contains(vbs,'ScatteringAngle'));
Index_size = find(contains(vbs,'diameter'));
time_stamp=PhipsData.RealTimeStamp;
Start_time = PhipsData.RealTimeStamp(1);
end_time = PhipsData.RealTimeStamp(end);
ScaIn = PhipsDataMatrix(:,...
    Index_scattering);
psize_c1c2 = PhipsDataMatrix(:,Index_size);
Tl = length(time_stamp);
sizeinfo=zeros(Tl,1);

% particle size larger than 26 micrometers from two cameras;
sizelimit=15; % > 26 micrometers;
for ns = 1:Tl
    p1 = psize_c1c2(ns,1);
    p2 = psize_c1c2(ns,2);
    p1nan = isnan(p1);
    p2nan = isnan(p2);
    if p1nan && ~p2nan
        sizeinfo(ns)=p2;
    end
    if ~p1nan && p2nan
        sizeinfo(ns)=p1;
    end
    if p1nan && p2nan 
        sizeinfo(ns)=0;
    end
    if ~p1nan && ~p2nan
        sizeinfo(ns)=(p1+p2)/2;
    end
    if sizeinfo(ns)>0 && sizeinfo(ns)<sizelimit
        sizeinfo(ns)=0;
    end
end

% filtering by size;
ids = find(sizeinfo);
Sca_data=ScaIn(ids,:);
Time_data = time_stamp(ids,:);
Size_data = sizeinfo(ids);

% filtering by intensity;
sca_temp = sum(Sca_data,2);
idds = find(~isnan(sca_temp));

% results;
time_s=Time_data(idds);
size_s = Size_data(idds);
scain_s = Sca_data(idds,:);
Nsignals = length(scain_s);
scalf=zeros(size(size_s));

% coefficient for polynomial-fitting of the asymmetry factor; 
poc = [-5.92702829979099e-05,0.00130522506929561,...
    -0.0108666805589366,0.0409324338201189,0.940294146231213];


%% retreival of asymmetry factor and phase function for single particle;
maxdeg=3000;
deg_range=(1:maxdeg)-1;
legnodes= legpts(maxdeg);
nodes_indeg = rad2deg(acos(legnodes));
angs_det = linspace(18,170,20);
angs_det = flip(angs_det');
gm = zeros(Nsignals,1);
Cp = gm;

% wait bar;
f = waitbar(0,'1','Name','computing single particle', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for k = 1:Nsignals
    if getappdata(f,'canceling')
        break
    end
    P11_m = scain_s(k,:);
    P11_m = flip(P11_m');
    method = 'nearest';
    %vq = interp1(angs2_e,P11_e,legdeg,method, 'extrap') ;
    vq1 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'makima';
    vq2 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'pchip';
    vq3 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'linear';
    vq4 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    vq =(vq1 + vq2 + vq4)/3;
    legcoefs_m = legvals2legcoeffs(vq);
    scalf(k) = legcoefs_m(1);
    legcoefs_nor=legcoefs_m/legcoefs_m(1);
    g_m = legcoefs_nor(2)/3;
    % Update waitbar and message
    waitbar(k/Nsignals,f,sprintf('%12.9f',g_m))
    
    % asymmetry factors; 
    coef = legcoefs_nor./(deg_range'*2+1);
    IoC = 1/sum(abs(coef));
    Cp(k) = IoC;
    gm(k) = g_m; 
    
end
delete(f)

%% filtering out extreme values

% ray tracing asymmetry factor < 0.2, consider as fake sampling;
[r] = find( gm<0.8&gm>0.15);
gmf_temp = gm(r);
ctf_temp = time_s(r);
Cmf_temp = Cp(r);
Csf_temp = size_s(r);
scalf_temp = scalf(r);
%phasef_temp = phasef(r,:);
scain_temp = scain_s(r,:);
% to insure accuracy filtering out ray-tracing complexity <0.25;
% [r2] = find( Cmf_temp > 0.1);
% gmf = gmf_temp(r2);
% ctf = ctf_temp(r2);
% Cmf = Cmf_temp(r2);
% Csf = Csf_temp(r2);
% scal_f = scalf_temp(r2);
% sca_m = scain_temp(r2,:);
gmf = gmf_temp;
ctf = ctf_temp;
Cmf = Cmf_temp;
Csf = Csf_temp;
scal_f = scalf_temp;
sca_m = scain_temp;


%% asymmetry factor due to diffraction by spherical particles;

Sp = Csf*pi/0.532;
gdif = polyval(poc, log(Csf)); 
% total asymmetry factor;
gav = (gmf+gdif)/2;
% ray-tracing coefficients; 

%% prepare the group average  ; 
% now we can group thing together ; 
% how many particles per group;
Ng = 20; 

NTT = length(gmf); 
Ngroup = idivide(int64(NTT),int64(Ng));
remd = mod(NTT,Ng);
g_group = zeros(Ngroup,1);
C_group = g_group;
M_group = zeros(Ngroup,20);
time_group = ctf(1:Ngroup);
size_group = g_group;
% get the group avergaed intensity; 
for kg = 1:Ngroup 
    k1 = 1+(kg-1)*Ng;
    k2 = k1+Ng-1;
    % intensity; 
    sca_tmp = sca_m(k1:k2,:);
    M_group(kg,:) = sum(sca_tmp,1)/Ng;
    % size 
    size_group(kg) = sum(Csf(k1:k2))/Ng;
    dur = ctf(k2)-ctf(k1);
    time_group(kg) = ctf(k1)+dur/2;
end
scalf_group = zeros(Ngroup,1);

%% ----retreival of asymmetry factor and phase function for a group of particles;

maxdeg = 15000; 
deg_range=(1:maxdeg)-1;
legnodes= legpts(maxdeg);
nodes_indeg = rad2deg(acos(legnodes));
N_trunc = 30;
% diffraction coefficeints; 
gdif_group = polyval(poc,log(size_group)); 
coe_difg = gdif_group.^(deg_range);
% phase functions;
phasef_group= zeros(Ngroup,maxdeg);
% number of truncation for fitting data. 
% less term, more smoothing; 
RT_coeg = zeros(Ngroup, N_trunc);

% wait bar;
%Ld = zeros(Ngroup, 1);


f = waitbar(0,'1','Name','computing particle group', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);
for k = 1:Ngroup
    if getappdata(f,'canceling')
        break
    end
    P11_m = M_group(k,:);
    P11_m = flip(P11_m');
    method = 'nearest';
    %vq = interp1(angs2_e,P11_e,legdeg,method, 'extrap') ;
    vq1 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'makima';
    vq2 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
%     method = 'pchip';
%     vq3 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    method = 'linear';
    vq4 = interp1(angs_det,P11_m,nodes_indeg,method, 'extrap') ;
    vq =(vq2 +vq1+vq4 )/3;
    legcoefs_m = legvals2legcoeffs(vq);
    scalf_group(k) = legcoefs_m(1);
    legcoefs_nor=legcoefs_m/legcoefs_m(1);
    g_m = legcoefs_nor(2)/3;
    % Update waitbar and message
    waitbar(k/Ngroup,f,sprintf('%12.9f',g_m))
    
    % asymmetry factors; 
    coef = legcoefs_nor./(deg_range'*2+1);
%     [rr]=find(abs(coef)>0.5*10e-3);
%     Ld(k) = rr(end);
    
    IoC = 1/sum(abs(coef));
    
    
    C_group(k) = IoC;
    g_group(k) = (g_m+gdif_group(k))/2;
     
    % recover phase functions;
    rt_coefs=zeros(size(coef));
    rt_coefs(1:N_trunc)=coef(1:N_trunc);
    %RT_coefs(k,:)= (coef(1:N_trunc))';
    dif_nor = coe_difg(k,:);
    %tot_nor = (rt_coefs + dif_nor')/2;
    %tot_coef = tot_nor.*(deg_range'*2+1);
    phasef_temp1 =  legcoeffs2legvals(rt_coefs.*(deg_range'*2+1));
    phasef_temp2 = legcoeffs2legvals(dif_nor'.*(deg_range'*2+1));
    phasef_group(k,:) =0.5*(phasef_temp1 + phasef_temp2); 
        
end
delete(f)

%%  --- phase function plots for group of particles;
%close(f6);
close all; 
f6 = figure;
% quasi-spherical crystals index : 96-107;
indices=[56;105;118];
rt=indices;
%N= length(g_group);
% select the index of a measurement, for phase function plotting; 
rns = round(double(Ngroup)*rand(3,1));
rns = indices;
asys = rns;
cps = rns;
h = zeros(1,2*length(rns));
kks = length(rns);
for ks = 1:kks
nt = rns(ks); 
phf_se=phasef_group(nt,:);
measured_se = flip(M_group(nt,:)); 
asys(ks)= g_group(rns(ks));
cps(ks) = C_group(rns(ks));
% match the data at the 42 degree  
fitangle=42;
lc = abs(nodes_indeg-fitangle);
[m1,id1]  =  min(lc);
lc2 = abs(angs_det-fitangle);
[m2,id2]  =  min(lc2);
% ratio of measured intensity to phase function; 
rt(ks) = phf_se(id1)/measured_se(id2);
h(ks)=semilogy(nodes_indeg, phf_se', ':.');hold on; 
%semilogy(angs2, measured_select',  'o');
title('Scattering phase function recovery');
xlabel('Scattering angles (deg)');
ylabel('Scattering phase function') 
end

disp(asys); 
disp(cps);

g1=asys(1);
g2=asys(2);
g3=asys(3);

cp1 = cps(1);
cp2 = cps(2);
cp3 = cps(3);

l1 = strcat("I: ", "$g=",num2str(g1), '$  and $',  "C_{p}=", num2str(cp1), "$"); 
l2 = strcat("II: ", "$g=",num2str(g2), '$ and $',  "C_{p}=", num2str(cp2), "$"); 
l3 = strcat("III:  ", "$g=",num2str(g3), '$ and $',  "C_{p}=", num2str(cp3), "$"); 

hold on; 
for ks = 1:length(rns)
    
nt = rns(ks); 
measured_select = flip(M_group(nt,:)); 
h(kks+ks) = scatter(angs_det, (measured_select')*rt(ks),  'o');

end
hold on; 
ylim([0.05 100]);
legend(h(1:(kks+1)), l1, ...
    l2,...
    l3,'measurements','interpreter','latex')
    title('')
    
%% Load aircraft data
addpath('/Users/nikewolf/Documents/GitHub/PHIPS-Scattering-Data-Analysis-Tools/Aircraft Data Functions')
% [time,airspeed,height,lat,lon,hFig1,hFig2,T,S_ice,w] = Aircraft_data('CIRRUS-HL',...
%     './','QL-CIRRUS-HL_F05_20210629a_ADLR_BAHAMAS_v1.nc');
[time,airspeed,height,lat,lon,hFig1,hFig2,T,S_ice,w]=read_flight_speed_CIRRUSHL('./',...
    'QL-CIRRUS-HL_F05_20210629a_ADLR_BAHAMAS_v1.nc');
time = datetime(time,'ConvertFrom','datenum');

% Preparations for aircraft data plotting
alttick1 = 0:2000:round(max(height),5)+1000;
ttick = round(min(T),-1)-10:10:round(max(T),2)+10;
PlotTimeWindow = [datetime(2021,6,29,10,0,0) datetime(2021,6,29,12,30,0)];


%% Timeseries of g and Cp

fsize=19; %fontsize
singleplotheight=0.19;


% Plot aircraft data
k1 = figure;
set(k1, 'Position', [0 00 1000 940]);
subplot('Position',[0.13,0.08+(singleplotheight+0.025)*3,0.775,singleplotheight])
plot(time,T,'linewidth',2,'color','b'); grid on;
set(gca,'xlim', [PlotTimeWindow(1),PlotTimeWindow(2)],'xticklabel',{[]},'ylim',[min(ttick),max(ttick)],'ytick',ttick,'Ycolor','k','FontSize',fsize-1); 
ylabel('$T [^{\circ}C]$', 'interpreter','latex','FontSize',fsize);
set(gca,'FontSize',fsize-4,'Color','w');
yyaxis right
plot(time,height,'linewidth',2,'color','k'); grid on;
set(gca,'xlim', [PlotTimeWindow(1),PlotTimeWindow(2)],'xticklabel',{[]},'ylim',[min(alttick1),max(height)+1000],'ytick',alttick1,'Ycolor','k','FontSize',fsize-1); 
ylabel('$Altitude [m.a.s.l.]$', 'interpreter','latex','FontSize',fsize)
legend('T_{ambient}','GPS altitude','Location','NorthWest');
set(gca,'FontSize',fsize-4,'Color','w');

subplot('Position',[0.13,0.08+(singleplotheight+0.025)*2,0.775,singleplotheight]);
sz = 35;
scatter(time_group,size_group,sz,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.635 0.078 .183],...
    'LineWidth',0.2)
hold on; grid on;
ylabel('$Diameter [\mu m]$', 'interpreter','latex','FontSize',fsize);
set(gca,'xlim', [PlotTimeWindow(1),PlotTimeWindow(2)],'xticklabel',{[]},'FontSize',fsize-1); 
ylim([5 250])
set(gca,'FontSize',fsize-4,'Color','w');

subplot('Position',[0.13,0.08+(singleplotheight+0.025)*1,0.775,singleplotheight]);
scatter(time_group,C_group,sz,'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',0.2)
hold on; grid on;
yline(0.4,'r.');
ylabel('$C_{p}$', 'interpreter','latex','FontSize',fsize);
ylim([0.3 0.7])
set(gca,'xlim', [PlotTimeWindow(1),PlotTimeWindow(2)],'xticklabel',{[]},'FontSize',fsize-1); 
set(gca,'FontSize',fsize-4,'Color','w');

subplot('Position',[0.13,0.08+(singleplotheight+0.025)*0,0.775,singleplotheight]);
scatter(time_group,g_group,sz,'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',0.2)
hold on; grid on;
yline(0.7,'r.');
ylabel('$g$', 'interpreter','latex','FontSize',fsize);
ylim([0.6 0.8])
xlim([PlotTimeWindow(1),PlotTimeWindow(2)])
xlabel('time (UTC) June 29 2021');
set(gca,'FontSize',fsize-4,'Color','w');

% Create textbox
annotation(k1,'textbox',...
    [0.133 0.728255319148936 0.04675 0.0409574468085107],'String',{'A'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica Neue');

% Create textbox
annotation(k1,'textbox',...
    [0.134 0.653787234042553 0.04725 0.0409574468085107],'String',{'B'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica Neue');

% Create textbox
annotation(k1,'textbox',...
    [0.132 0.441021276595744 0.04825 0.0409574468085107],'String',{'C'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica Neue');

% Create textbox
annotation(k1,'textbox',...
    [0.132 0.226127659574468 0.04775 0.0409574468085107],'String',{'D'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica Neue');
print('./CIRRUSHL_g_RF05_group20.png','-dpng')


%%
k2 = figure;
scatter( C_group, g_group , sz,'MarkerEdgeColor',[0 .5 .5],...
    'MarkerFaceColor',[0 .7 .7],...
    'LineWidth',0.2);
ylabel('$g$', 'interpreter','latex');
xlabel('$C_{p}$', 'interpreter','latex');
xlim([0.1 0.9])

k3 = figure;
histfit(g_group, 10,'normal');
xlabel('$g$ distribution', 'interpreter','latex');
ylabel('Counts')
pfg = fitdist(g_group, 'Normal');

%close figure(4);
k4 =  figure;
histfit(C_group, 10,'normal');
xlabel('$C_{p}$ distribution', 'interpreter','latex');
ylabel('Counts')
pfc = fitdist(C_group, 'Normal');

% k5 = figure; 
% scatter(C_group, Ld); ylim([1 20]);xlim([0.1 0.7]);
% 
% k6 =  figure;
% histfit(Ld, 10,'normal');
% xlabel('$C_{p}$ distribution', 'interpreter','latex');
% pfc = fitdist(Ld, 'Normal');





