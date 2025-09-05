%% 
lat = double(ncread('uvuvuvuv.nc', 'latitude')); % 15, -55
lon = double(ncread('uvuvuvuv.nc', 'longitude')); % -85,-30
u = double(ncread('uvuvuvuv.nc', 'u')); % Componente U del viento
v = double(ncread('uvuvuvuv.nc', 'v')); % Componente V del viento
date_ndefm = double(ncread('uvuvuvuv.nc', 'date')); % Fechas nov-marzo, 1940-2023


[lat_grid, lon_grid] = meshgrid(lat, lon);
weights = sqrt(cosd(lat_grid));
u_weighted = u .* weights;
v_weighted = v .* weights;

u_weighted = squeeze(u_weighted);
v_weighted = squeeze(v_weighted);

u_weighted_reshaped = reshape(u_weighted, numel(lat)*numel(lon), size(u_weighted, 3));
nanindex_u = isnan(nanmean(u_weighted_reshaped, 2));
u_weighted_reshaped = u_weighted_reshaped(~nanindex_u, :);


% Detrend data along the temporal axis
F_u = detrend(u_weighted_reshaped, 0);

% Compute covariance matrix
C_u = (F_u * F_u') / (size(F_u, 2) - 1);

% Eigen decomposition
[EOF_u, D_u] = eig(C_u);
PC_u = EOF_u' * F_u;

% Selecciona el primer EOF
EOF1_u = EOF_u(:, end); 
PC1_u = PC_u(end, :);

%% EOF - Direct Calculation
load('ssta');
[lats,lons]=meshgrid(lat,lon);
u=u.*repmat(sqrt(cosd(lats)),1,1,504);

UUU=(reshape(u,91*31,504))';
nanindex=isnan(nanmean(ssta));
ssta=ssta(:,~nanindex);

F=detrend(ssta,0);
C=F'*F;

[EOFs,D]=eig(C);
PCs=EOFs'*F';

EOF1=EOFs(:,end);
PC1=PCs(end,:);

sEOF1=NaN(91*31,1);
sEOF1(~nanindex)=EOF1;
sEOF1=reshape(sEOF1,91,31);

sEOF1=sEOF1.*nanstd(PC1);
PC1=PC1./nanstd(PC1);

%% EOF - CDT
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);
[eof_maps,pc,expvar]=eof(ssta);
eof1=eof_maps(:,:,1);
pc1=(pc(1,:))';
eof1=eof1.*nanstd(pc1);
pc1=pc1./nanstd(pc1);

%% EOF - PCA
load('ssta');
[lats,lons]=meshgrid(lata,lona);
ssta=ssta.*repmat(sqrt(cosd(lats)),1,1,504);

ssta=(reshape(ssta,91*31,504))';
nanindex=isnan(nanmean(ssta));
ssta=ssta(:,~nanindex);

F=detrend(ssta,0);
[coef,score,latent]=pca(F);

scoef1=NaN(91*31,1);
scoef1(~nanindex)=coef(:,1);
scoef1=reshape(scoef1,91,31);

score1=score(:,1);
scoef1=scoef1.*nanstd(score1);
score1=score1./nanstd(score1);


%%

figure
subplot(2,1,1);
m_proj('miller','lon',[nanmin(lona) nanmax(lona)],'lat',[nanmin(lata) nanmax(lata)]);
m_contourf(lona,lata,sEOF1',linspace(-1.4,1.4,200),'linestyle','none');
m_coast('patch',[0.7 0.7 0.7],'linewidth',2);
m_grid('linewidth',2,'fontname','consolas');
colormap(m_colmap('diverging'));
caxis([-1.4 1.4]);
s=colorbar('fontname','consolas','fontsize',12);
title(s,'^{o}C','fontname','consolas');
set(gca,'fontsize',12)
title('EOF1: 34.92%','fontsize',16,'fontname','consolas');

subplot(2,1,2);
plot(1:504,PC1,'r','linewidth',2);
set(gca,'xtick',[6:60:504],'xticklabels',1980:5:2021,'fontname','consolas','fontsize',12);
xlabel('Year','fontname','consolas');
ylabel('PC1','fontname','consolas');
xlim([1 504]);
set(gca,'fontsize',12,'linewidth',2)
title('PC1: 34.92%','fontsize',16,'fontname','consolas');

%% 

