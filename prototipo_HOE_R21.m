function prototipo_HOE_R21(cadFile,visualScale)
% Presentacion R2.1: no cambia fisica ni recalcula el holograma.
% Usa resultados R2 del mismo CAD y logo. STL en mm, calculos en m.
root=fileparts(mfilename('fullpath'));
if nargin<1, cadFile=fullfile(root,'parabrisas.stl'); end
if nargin<2, visualScale=25; end
assert(isscalar(visualScale)&&isfinite(visualScale)&&visualScale>0);
base=fullfile(root,'resultados_HOE_R2','datos_HOE.mat');
assert(exist(base,'file')==2,'Ejecutar R2 primero.');
s=load(base);
tr=stlread(cadFile); V=tr.Points*1e-3; faces=tr.ConnectivityList;
out=fullfile(root,'resultados_HOE_R21'); if ~exist(out,'dir'), mkdir(out); end
R=s.R; C=s.center; world=@(p) p*R+C;
E=world([0 0 s.eyeDist]); H=world([0 0 0]); Q=world([0 0 -s.D]);
S=world(-.5*[sin(s.theta) 0 -cos(s.theta)]);
% Intensidad calculada, no logo original superpuesto.
ix=find(abs(s.xe)<=.004); iy=find(abs(s.xe)<=.003);
xx=s.xe(ix); yy=s.xe(iy); img=s.imgs{3}(iy,ix);
img=max(0,min(1,img)); displayI=img.^0.55; % gamma solo visual
[xxg,yyg]=meshgrid(xx,yy);
rgb=cat(3,.08*displayI,displayI,.15*displayI);
alpha=displayI; alpha(img<0.012)=0;
fig=figure('Color','w','Name','R2.1 - Sistema HOE','Position',[100 100 1250 700]);
ax=axes(fig); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
glass(ax,V,faces);
% HOE curvo aproximado, tomado del ajuste R2
[a,b]=meshgrid(s.x(1:32:end)); hh=s.h(1:32:end,1:32:end);
drawSurface(ax,world,a,b,hh,[.1 .7 .25],.9);
plot3(ax,[S(1) H(1)],[S(2) H(2)],[S(3) H(3)],'-','Color',[.9 .55 .12],'LineWidth',2);
plot3(ax,[H(1) E(1)],[H(2) E(2)],[H(3) E(3)],'-','Color',[.12 .35 .85],'LineWidth',2);
plot3(ax,[H(1) Q(1)],[H(2) Q(2)],[H(3) Q(3)],'--','Color',[.12 .35 .85],'LineWidth',1.5);
plot3(ax,E(1),E(2),E(3),'ko','MarkerFaceColor','k');
text(ax,E(1),E(2),E(3),'  Ojo');
text(ax,H(1),H(2),H(3),'  HOE 20 x 20 mm');
text(ax,S(1),S(2),S(3),'  Iluminacion (30 grados)');
text(ax,Q(1),Q(2),Q(3),sprintf('  Plano virtual: 2 m | logo x%g',visualScale));
drawLogo(ax,world,xxg,yyg,-s.D,rgb,alpha,visualScale);
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
title(ax,sprintf('Sistema optico | logo ampliado x%g SOLO para visualizacion',visualScale));
view(ax,3); camlight(ax,'headlight'); lighting(ax,'gouraud'); axis(ax,'vis3d');
rotate3d(fig,'on');
savefig(fig,fullfile(out,'sistema_3D.fig')); saveas(fig,fullfile(out,'sistema_3D.png'));
% Vista del observador. La textura virtual esta a escala fisica.
% No es una simulacion de retina: se muestra la reconstruccion inversa.
fig2=figure('Color','w','Name','R2.1 - Vista del observador','Position',[100 100 1250 600]);
ax2=subplot(1,2,1); hold(ax2,'on'); glass(ax2,V,faces);
drawLogo(ax2,world,xxg,yyg,-s.D,rgb,alpha,1);
axis(ax2,'equal'); axis(ax2,'off');
campos(ax2,E); camtarget(ax2,Q); camup(ax2,R(2,:)); camproj(ax2,'perspective');
camva(ax2,0.35); camlight(ax2,'headlight'); lighting(ax2,'gouraud');
title(ax2,{'Desde el ojo: vista angular ampliada','Logo a escala fisica, textura numerica'});
ax3=subplot(1,2,2); imagesc(ax3,xx*1e3,yy*1e3,rgb); axis(ax3,'image'); axis(ax3,'xy');
xlabel(ax3,'x virtual (mm)'); ylabel(ax3,'y virtual (mm)');
title(ax3,{'Reconstruccion desde pupila a 2 m','Gamma visual; datos fisicos sin cambios'});
savefig(fig2,fullfile(out,'vista_observador.fig')); saveas(fig2,fullfile(out,'vista_observador.png'));
% Comparacion de resultados. El caso curvo no compensado requiere muestreo.
fig3=figure('Color','w','Name','R2.1 - Resultados','Position',[100 100 1250 700]);
for j=1:4
    axj=subplot(2,3,j);
    if j==1, im=s.Te; lab='Objetivo'; else, im=s.imgs{j-1}; lab=s.labels{j-1}; end
    imagesc(axj,s.xe*1e3,s.xe*1e3,im); axis(axj,'image'); axis(axj,'xy'); xlim(axj,[-4 4]); ylim(axj,[-3 3]); caxis(axj,[0 1]);
    if j==3 && s.phaseStep>pi, lab='Curvo: submuestreado (no cuantitativo)'; end
    title(axj,lab); xlabel(axj,'mm'); ylabel(axj,'mm');
end
colormap(fig3,gray);
ax4=subplot(2,3,5); hold(ax4,'on');
plot(ax4,s.depths,s.axial(:,1),'b-','LineWidth',2);
plot(ax4,s.depths,s.axial(:,3),'--','Color',[0 .6 .25],'LineWidth',2);
if s.phaseStep<=pi
    plot(ax4,s.depths,s.axial(:,2),'Color',[.85 .4 .1]); legend(ax4,'Plano','Compensado ideal','Curvo','Location','best');
else
    legend(ax4,'Plano','Compensado ideal','Location','best');
end
grid(ax4,'on'); xlabel(ax4,'Distancia virtual (m)'); ylabel(ax4,'Similitud'); title(ax4,'Barrido axial');
ax5=subplot(2,3,6); axis(ax5,'off');
text(ax5,0,.9,sprintf('Captura plana: %.2f %%',100*s.capture(1)));
text(ax5,0,.75,sprintf('Similitud plana: %.3f',s.similarity(1)));
text(ax5,0,.6,sprintf('Pico plano: %.2f m',s.peakDepth(1)));
text(ax5,0,.45,sprintf('Paso fase: %.2f rad (limite pi)',s.phaseStep));
text(ax5,0,.25,{'Compensacion ideal por construccion','Sin refraccion multicapa'},'FontSize',10);
saveas(fig3,fullfile(out,'resultados.png'));
fprintf('R2.1: visualizacion solamente; logo global x%g; vista del ojo a escala fisica.\n',visualScale);
if s.phaseStep>pi, fprintf('ADVERTENCIA: caso curvo submuestreado; no interpretar perdida como cuantitativa.\n'); end
end
function glass(ax,V,F)
patch(ax,'Faces',F,'Vertices',V,'FaceColor',[.3 .7 .9],'FaceAlpha',.20,'EdgeColor','none');
end
function drawSurface(ax,world,x,y,z,c,a)
p=world([x(:),y(:),z(:)]);
surf(ax,reshape(p(:,1),size(x)),reshape(p(:,2),size(x)),reshape(p(:,3),size(x)),'FaceColor',c,'FaceAlpha',a,'EdgeColor','none');
end
function drawLogo(ax,world,x,y,z,rgb,alpha,scale)
p=world([scale*x(:),scale*y(:),z*ones(numel(x),1)]);
surf(ax,reshape(p(:,1),size(x)),reshape(p(:,2),size(x)),reshape(p(:,3),size(x)),rgb,'FaceColor','texturemap','EdgeColor','none','FaceAlpha','texturemap','AlphaData',alpha,'AlphaDataMapping','none','FaceLighting','none');
end
