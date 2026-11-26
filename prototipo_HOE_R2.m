function prototipo_HOE_R2(cadFile,logoFile)
% R2: curvatura local del vidrio y compensacion escalar de fase.
% Reutiliza el diseno R1. No modela capas ni eficiencia volumetrica.
root=fileparts(mfilename('fullpath'));
if nargin<1, cadFile=fullfile(root,'parabrisas.stl'); end
if nargin<2, logoFile=fullfile(root,'logo.png'); end
assert(exist(cadFile,'file')==2,'Se requiere el STL real.');
assert(exist(logoFile,'file')==2,'Se requiere el logo PNG utilizado en R1.');
base=fullfile(root,'resultados_HOE_R1','datos_HOE.mat');
if ~exist(base,'file'), prototipo_HOE_R1(cadFile,logoFile); end
% Si se cambia STL o PNG, borrar resultados_HOE_R1 y volver a ejecutar.
s=load(base); assert(s.cadUsed==1,'Regenerar R1 con el CAD real.');
out=fullfile(root,'resultados_HOE_R2'); if ~exist(out,'dir'), mkdir(out); end
lambda=s.lambda; k=2*pi/lambda; W=s.W; dx=s.dx;
x=s.x; N=numel(x); [X,Y]=meshgrid(x);
D=s.D; eyeDist=s.eyeDist; theta=s.theta; De=D+eyeDist;
R=s.R; center=s.center; A=double(abs(X)<W/2 & abs(Y)<W/2);
U=A.*exp(1i*s.phi);
%% 1. Medir la superficie local del CAD
tr=stlread(cadFile); V=tr.Points*1e-3; faces=tr.ConnectivityList;
% Coordenadas locales: x,y en plano tangente, z hacia el ojo.
v=(V-center)*R';
p1=v(faces(:,1),:); p2=v(faces(:,2),:); p3=v(faces(:,3),:);
normals=cross(p2-p1,p3-p1,2);
normals=normals./max(sqrt(sum(normals.^2,2)),eps);
c=(p1+p2+p3)/3;
% Seleccionar caras de la misma orientacion y cercanas en altura.
% Evitar mezclar superficies opuestas del laminado.
good=normals(:,3)>0.8;
if nnz(good)<6, good=normals(:,3)<-0.8; end
fitOK=false;
for radius=[0.025 0.04 0.06 0.09]
    use=good & hypot(c(:,1),c(:,2))<=radius & abs(c(:,3))<0.004;
    ids=unique(faces(use,:)); q=v(ids,:);
    if size(q,1)<12, continue; end
    a=q(:,1)/radius; b=q(:,2)/radius;
    G=[ones(size(a)),a,b,a.^2,a.*b,b.^2];
    if rank(G)<6 || cond(G)>1e4, continue; end
    coeff=G\q(:,3);
    res=q(:,3)-G*coeff;
    % Descartar puntos alejados del ajuste (posibles otras capas).
    keep=abs(res-median(res))<=max(5e-5,3*median(abs(res-median(res))));
    if nnz(keep)<12 || rank(G(keep,:))<6, continue; end
    coeff=G(keep,:)\q(keep,3); qfit=q(keep,:);
    if any(min(qfit(:,1:2))> -W/2) || any(max(qfit(:,1:2))<W/2), continue; end
    fitRMS=sqrt(mean((G(keep,:)*coeff-qfit(:,3)).^2));
    fitOK=true; break
end
assert(fitOK,'Malla insuficiente para ajuste local. Exportar STL mas fino o revisar cara.');
height=@(a,b) coeff(1)+coeff(2)*a/radius+coeff(3)*b/radius+ ...
    coeff(4)*(a/radius).^2+coeff(5)*a.*b/radius^2+coeff(6)*(b/radius).^2;
h=height(X,Y); h=h-height(0,0); % piston irrelevante
curv=[2*coeff(4),coeff(5);coeff(5),2*coeff(6)]/radius^2;
curvatures=eig(curv); sagPV=max(h(A>0))-min(h(A>0));
%% 2. Fase adicional por altura del vidrio (modelo reflectivo de aire)
% Incidencia k_i,z=-k*cos(theta), salida de referencia k_o,z=+k.
% Delta phi=(k_i,z-k_o,z)*h. Aproximacion de pantalla de fase.
delta=-k*(1+cos(theta))*h;
% Tres casos: plano, curvo sin compensar, curvo con compensacion.
% La cancelacion ideal es algebraica: no es validacion independiente
% ni prueba de fabricabilidad. No se incluye refraccion en el laminado.
Ucurved=U.*exp(1i*delta);
phiComp=angle(exp(1i*(s.phi-delta)));
Ucomp=A.*exp(1i*(phiComp+delta));
phaseStep=max([max(abs(diff(delta,1,1)),[],'all'),max(abs(diff(delta,1,2)),[],'all')]);
if phaseStep>pi
    warning('Fase de curvatura submuestreada: caso sin compensar no cuantitativo.');
end
%% 3. Misma propagacion de R1, misma pupila y mismo plano virtual
M=2*N; f=(-M/2:M/2-1)/(M*dx); [KX,KY]=meshgrid(2*pi*f);
KZ=sqrt(max(0,k^2-KX.^2-KY.^2)); idx=N/2+(1:N);
flim=1/(lambda*sqrt(1+(2*eyeDist/(M*dx))^2));
Hprop=exp(1i*(KZ-k)*eyeDist).*double(abs(KX)<=2*pi*flim & abs(KY)<=2*pi*flim);
clear KX KY KZ
pupil=double(X.^2+Y.^2<=0.0025^2);
xe=s.xe; Te=s.Te; target=sqrt(max(Te,0));
fields={U,Ucurved,Ucomp}; labels={'Plano R1','Curvo sin compensar','Curvo compensado ideal'};
imgs=cell(1,3); eyes=imgs; capture=zeros(1,3); similarity=capture;
for j=1:3
    p=zeros(M,'like',U); p(idx,idx)=fields{j};
    p=fftshift(ifft2(ifftshift(fftshift(fft2(ifftshift(p))).*Hprop)));
    eyes{j}=p(idx,idx).*pupil;
    capture(j)=sum(abs(eyes{j}(:)).^2)/sum(abs(U(:)).^2);
    a=reconstruct(eyes{j},X,Y,k,De,N);
    similarity(j)=overlap(a,target);
    imgs{j}=abs(a).^2; imgs{j}=imgs{j}/max(max(imgs{j}(:)),eps);
end
depths=s.depths; axial=zeros(numel(depths),3);
for j=1:numel(depths)
    z=depths(j); zz=lambda*z/W; xz=(-N/2:N/2-1)*zz;
    [XZ,YZ]=meshgrid(xz);
    t=sqrt(max(0,interp2(s.xv,s.xv,s.T,XZ,YZ,'linear',0)));
    for m=1:3, axial(j,m)=overlap(reconstruct(eyes{m},X,Y,k,z,N),t); end
end
[~,ix]=max(axial,[],1); peakDepth=depths(ix);
%% 4. Figuras faciles de interpretar
fig=figure('Color','w','Name','R2 - Comparacion optica');
subplot(2,3,1); imagesc(xe*1e3,xe*1e3,Te); axis image xy; xlim([-4 4]); ylim([-3 3]); title('Logo objetivo');
for j=1:3
    subplot(2,3,j+1); imagesc(xe*1e3,xe*1e3,imgs{j}); axis image xy;
    xlim([-4 4]); ylim([-3 3]); title(labels{j});
end
subplot(2,3,5); imagesc(x*1e3,x*1e3,angle(exp(1i*delta))); axis image xy; title('Fase extra (rad)'); colorbar;
subplot(2,3,6); plot(depths,axial,'LineWidth',1.5); grid on; xlabel('Distancia virtual (m)'); ylabel('Similitud'); legend(labels,'Location','best');
colormap(gray); saveas(fig,fullfile(out,'comparacion.png'));
fig=figure('Color','w','Name','R2 - Curvatura CAD');
subplot(1,2,1); surf(X(1:16:end,1:16:end)*1e3,Y(1:16:end,1:16:end)*1e3,h(1:16:end,1:16:end)*1e6,'EdgeColor','none'); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('Altura (um)'); title('Ajuste local del vidrio'); colorbar; view(3);
subplot(1,2,2); scatter3(qfit(:,1)*1e3,qfit(:,2)*1e3,qfit(:,3)*1e6,12,'filled'); grid on; xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (um)'); title('Puntos STL usados en ajuste');
saveas(fig,fullfile(out,'curvatura.png'));
fig=figure('Color','w','Name','R2 - Prototipo CAD');
subplot(1,2,1); hold on; axis equal; grid on; view(3);
patch('Faces',faces,'Vertices',V,'FaceColor',[.3 .7 .9],'FaceAlpha',.25,'EdgeColor','none');
world=@(p) p*R+center;
q=world([-W/2 -W/2 0; W/2 -W/2 0; W/2 W/2 0; -W/2 W/2 0]); patch(q(:,1),q(:,2),q(:,3),[.1 .8 .3]);
E=world([0 0 eyeDist]); C=world([0 0 0]); Vp=world([0 0 -D]); S=world(-.5*[sin(theta) 0 -cos(theta)]);
plot3([S(1) C(1)],[S(2) C(2)],[S(3) C(3)],'Color',[.85 .5 .1],'LineWidth',2);
plot3([C(1) E(1)],[C(2) E(2)],[C(3) E(3)],'b','LineWidth',2);
plot3([C(1) Vp(1)],[C(2) Vp(2)],[C(3) Vp(3)],'b--');
text(E(1),E(2),E(3),' Ojo'); text(C(1),C(2),C(3),' HOE'); text(S(1),S(2),S(3),' Iluminacion'); text(Vp(1),Vp(2),Vp(3),' Plano virtual: 2 m del ojo');
% Marco fisico del logo, sin marcador verde ni ampliacion geometrica.
q=world([-.003 -.002 -D;.003 -.002 -D;.003 .002 -D;-.003 .002 -D;-.003 -.002 -D]); plot3(q(:,1),q(:,2),q(:,3),'g-');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('CAD y rayos centrales de referencia'); camlight; lighting gouraud;
subplot(1,2,2); imagesc(xe*1e3,xe*1e3,imgs{3}); axis image xy; xlim([-4 4]); ylim([-3 3]); xlabel('x virtual (mm)'); ylabel('y virtual (mm)'); title('Vista ampliada desde pupila (modelo ideal)'); colormap(gca,[zeros(256,1),linspace(0,1,256)',zeros(256,1)]);
saveas(fig,fullfile(out,'prototipo.png'));
%% 5. Resultados y advertencias
fprintf('Ajuste CAD: %d puntos, radio %.1f mm, RMS %.2f um\n',size(qfit,1),radius*1e3,fitRMS*1e6);
fprintf('Variacion de altura HOE: %.2f um; paso max fase %.3f rad\n',sagPV*1e6,phaseStep);
for j=1:3, fprintf('%s: captura %.3f %%, similitud %.4f, profundidad %.2f m\n',labels{j},100*capture(j),similarity(j),peakDepth(j)); end
fprintf('Compensacion ideal, sin refraccion ni eficiencia volumetrica. RMS del ajuste no equivale a precision optica.\n');
save(fullfile(out,'datos_HOE.mat'),'lambda','W','dx','D','eyeDist','theta','x','xe','Te','R','center','coeff','radius','fitRMS','curvatures','sagPV','h','delta','phiComp','phaseStep','capture','similarity','depths','axial','peakDepth','imgs','labels');
end
function a=reconstruct(u,X,Y,k,z,N)
a=fftshift(ifft2(ifftshift(u.*exp(-1i*k*(X.^2+Y.^2)/(2*z)))))*N;
end
function s=overlap(a,t)
a=abs(a(:)); t=t(:); s=sum(a.*t)/max(norm(a)*norm(t),eps);
end
