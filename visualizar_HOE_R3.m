function visualizar_HOE_R3(logoMM,cadFile)
% Visualizacion ligera de resultados existentes. No recalcula la optica.
if nargin<1, logoMM=50; end
root=fileparts(mfilename('fullpath'));
if nargin<2, cadFile=fullfile(root,'parabrisas.stl'); end
out=fullfile(root,sprintf('resultados_HOE_R3_%gmm',logoMM));
file=fullfile(out,'datos_HOE.mat');
assert(exist(file,'file')==2,'Esperar a que R3 termine y guarde datos_HOE.mat.');
s=load(file,'R','center','W','D','eyeDist','theta','logoW','logoH','xe','Ie','capture','score1','bestDepth');
tr=stlread(cadFile); V=tr.Points*1e-3; F=tr.ConnectivityList;
% Reducir solo la malla de PRESENTACION, no la geometria de calculo.
if size(F,1)>8000
    [F,V]=reducepatch(F,V,8000);
end
world=@(p) p*s.R+s.center;
E=world([0 0 s.eyeDist]); H=world([0 0 0]); Q=world([0 0 -s.D]);
S=world(-.5*[sin(s.theta) 0 -cos(s.theta)]);
ix=find(abs(s.xe)<=.6*s.logoW); iy=find(abs(s.xe)<=.7*s.logoH);
% Remuestreo de intensidad para dibujo, sin texturas ni transparencia.
nx=min(160,numel(ix)); ny=min(80,numel(iy));
xx=linspace(s.xe(ix(1)),s.xe(ix(end)),nx);
yy=linspace(s.xe(iy(1)),s.xe(iy(end)),ny);
[X,Y]=meshgrid(xx,yy);
I=interp2(s.xe,s.xe,s.Ie,X,Y,'linear',0);
I=max(0,min(1,I)); G=I.^.55;
mask=I>.025; p=world([X(mask),Y(mask),-s.D*ones(nnz(mask),1)]);
colors=[.08*G(mask),G(mask),.15*G(mask)];
fig=figure('Color','w','Name','R3 - Visualizacion ligera','Position',[80 80 1200 650]);
ax=subplot(1,2,1); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
% Malla opaca: sin alpha, iluminacion ni texturas.
patch(ax,'Faces',F,'Vertices',V,'FaceColor',[.72 .85 .9],'EdgeColor','none','FaceLighting','none');
q=world([-s.W/2 -s.W/2 0;s.W/2 -s.W/2 0;s.W/2 s.W/2 0;-s.W/2 s.W/2 0]);
patch(ax,q(:,1),q(:,2),q(:,3),[.05 .6 .2]);
ray(ax,S,H,[.85 .5 .1],'-'); ray(ax,H,E,[.1 .3 .85],'-'); ray(ax,H,Q,[.1 .3 .85],'--');
scatter3(ax,p(:,1),p(:,2),p(:,3),5,colors,'filled');
plot3(ax,E(1),E(2),E(3),'ko','MarkerFaceColor','k');
text(ax,E(1),E(2),E(3),' Ojo'); text(ax,H(1),H(2),H(3),' HOE'); text(ax,S(1),S(2),S(3),' Iluminacion'); text(ax,Q(1),Q(2),Q(3),' Imagen virtual');
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)');
title(ax,'Sistema a escala real (CAD opaco para dibujo)'); view(ax,3);
axis(ax,'tight'); camproj(ax,'orthographic');
ax2=subplot(1,2,2);
% Vista optica numerica del observador, no camara 3D ni retina simulada.
rgb=cat(3,.08*G,G,.15*G);
imagesc(ax2,xx*1e3,yy*1e3,rgb); axis(ax2,'image'); axis(ax2,'xy');
xlabel(ax2,'x virtual (mm)'); ylabel(ax2,'y virtual (mm)');
title(ax2,{sprintf('Logo calculado desde pupila: %.0f mm',s.logoW*1e3),sprintf('Captura %.2f %% | similitud %.3f | pico %.2f m',100*s.capture,s.score1,s.bestDepth)});
drawnow;
saveas(fig,fullfile(out,'prototipo_ligero.png')); savefig(fig,fullfile(out,'prototipo_ligero.fig'));
% Vista global adicional en alambre, con logo fisico, sin oclusion opaca.
fig2=figure('Color','w','Name','R3 - Vista geometrica en alambre');
ax=axes(fig2); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
patch(ax,'Faces',F,'Vertices',V,'FaceColor','none','EdgeColor',[.7 .8 .85]);
scatter3(ax,p(:,1),p(:,2),p(:,3),5,colors,'filled');
ray(ax,S,H,[.85 .5 .1],'-'); ray(ax,H,E,[.1 .3 .85],'-'); ray(ax,H,Q,[.1 .3 .85],'--');
text(ax,E(1),E(2),E(3),' Ojo'); text(ax,Q(1),Q(2),Q(3),' Imagen virtual a 2 m');
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)'); title(ax,'Geometria y logo virtual a escala real');
view(ax,3); axis(ax,'tight'); camproj(ax,'orthographic'); drawnow;
saveas(fig2,fullfile(out,'geometria_ligera.png')); savefig(fig2,fullfile(out,'geometria_ligera.fig'));
end
function ray(ax,a,b,c,style)
plot3(ax,[a(1) b(1)],[a(2) b(2)],[a(3) b(3)],style,'Color',c,'LineWidth',1.5);
end
