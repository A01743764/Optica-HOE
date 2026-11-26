function HOE_Final(cadFile,logoFile)
% RFinal: simulacion escalar de un HOE reflectivo sobre parabrisas.
% Ejecutar HOE_Final con parabrisas.stl y logo.png junto a este archivo.
% HOE_Final('ruta/parabrisas.stl','ruta/logo.png') tambien es valido.
% No necesita scripts ni resultados MAT de revisiones anteriores.
% Conserva el algoritmo y parametros de R3 (50 mm).
% Todas las distancias en SI, salvo entrada STL en mm.
root=fileparts(mfilename('fullpath'));
if nargin<1, cadFile=fullfile(root,'parabrisas.stl'); end
if nargin<2, logoFile=fullfile(root,'logo.png'); end
assert(exist(cadFile,'file')==2,'Falta parabrisas.stl.');
assert(exist(logoFile,'file')==2,'Falta logo.png.');
out=fullfile(root,'resultados_Final'); if ~exist(out,'dir'), mkdir(out); end
logoMM=50;
fprintf('HOE Final: logo 50 mm, imagen virtual a 2 m.\n');
% Ajuste local antes de las FFT: detectar problemas de CAD al inicio.
tr=stlread(cadFile); Pcad=tr.Points*1e-3; faces=tr.ConnectivityList;
[R,center,coeff,radius,fitRMS,qfit]=ajusteCAD(Pcad,faces,0.020);
fprintf('CAD: %d puntos de ajuste, radio %.1f mm, RMS %.2f um.\n',size(qfit,1),radius*1e3,fitRMS*1e6);
%% 1. Parametros fisicos y malla numerica
lambda=532e-9; k=2*pi/lambda; W=0.020; N=2048;
dx=W/N; x=(-N/2:N/2-1)*dx; [X,Y]=meshgrid(x);
D=1.4; eyeDist=0.6; pupilRadius=0.0025; theta=30*pi/180;
% +z local hacia el ojo; imagen virtual z=-D. Incidencia hacia -z.
% Portadora oblicua se trata analiticamente, no se muestrea en la FFT.

% Fresnel escalado unitario desde HOE hasta plano virtual (-D).
dv=lambda*D/W; xv=(-N/2:N/2-1)*dv; [XV,YV]=meshgrid(xv);
C1=exp(-1i*k*(X.^2+Y.^2)/(2*D));
C2=exp(-1i*k*(XV.^2+YV.^2)/(2*D));
F=@(u) C2.*(fftshift(ifft2(ifftshift(C1.*u)))*N);
B=@(u) conj(C1).*(fftshift(fft2(ifftshift(conj(C2).*u)))/N);
A=double(abs(X)<W/2 & abs(Y)<W/2);
%% 2. Imagen objetivo (mismo procesamiento que R3)
[im,map,alpha]=imread(logoFile);
if ~isempty(map), im=reshape(map(double(im)+1,:),[size(im,1),size(im,2),3]); end
if isinteger(im), im=double(im)/double(intmax(class(im))); else, im=double(im); end
if size(im,3)>1, im=mean(im(:,:,1:3),3); end
if ~isempty(alpha), im=im.*(double(alpha)/max(1,max(double(alpha(:))))); end
assert(max(im(:))>0,'Logo vacio.');
% Recortar margenes negros para que ancho solicitado sea el ancho del contenido.
mask=im>0.05*max(im(:)); im=im(any(mask,2),any(mask,1));
assert(all(size(im)>=2),'Logo demasiado pequeno.');
% Mantener proporcion: ancho objetivo de 50 mm, centrado.
logoW=logoMM/1000; assert(logoW<=0.75*N*dv,'Plano virtual insuficiente.'); logoH=logoW*size(im,1)/size(im,2);
assert(logoH<=0.75*N*dv,'El logo no cabe en el plano virtual con el ancho fijado.');
% interp2 requiere imagen de al menos 2x2.
T=interp2(linspace(-logoW/2,logoW/2,size(im,2)), ...
    linspace(-logoH/2,logoH/2,size(im,1)),flipud(im),XV,YV,'nearest',0);
T=max(T,0); T=T/max(T(:)); target=sqrt(T); target=target*norm(A(:))/norm(target(:));
%% 3. Propagacion de Fresnel entre HOE, pupila y plano virtual
% Se calcula solamente el campo necesario para el ojo; el adjunto
% permite devolver la informacion de la imagen al plano del holograma.
pidx=find(abs(x)<=pupilRadius);
M=ceil((N+numel(pidx)-1)/128)*128;
dif=(pidx(1)-N):(pidx(end)-1);
assert(dx*max(abs(dif*dx))/(lambda*eyeDist)<0.5,'Muestreo Fresnel insuficiente.');
kernel=zeros(1,M);
kernel(mod(dif,M)+1)=exp(1i*k*(dif*dx).^2/(2*eyeDist));
kf=fft(kernel);
Hprop=(dx^2/(1i*lambda*eyeDist))*(kf.'*kf);
clear kernel kf
pupil=double(X.^2+Y.^2<=pupilRadius^2);
De=D+eyeDist; dve=lambda*De/W; xe=(-N/2:N/2-1)*dve; [XE,YE]=meshgrid(xe);
Ceye=exp(-1i*k*(X.^2+Y.^2)/(2*De));
Ceout=exp(-1i*k*(XE.^2+YE.^2)/(2*De));
Fe=@(u) Ceout.*(fftshift(ifft2(ifftshift(Ceye.*u)))*N);
Be=@(u) conj(Ceye).*(fftshift(fft2(ifftshift(conj(Ceout).*u)))/N);
P=@(u) pupilProp(u,Hprop,pidx,M,N,false);
Pt=@(u) pupilProp(u,Hprop,pidx,M,N,true);
L=@(u) Fe(pupil.*P(u));
Lt=@(v) Pt(pupil.*Be(v));
Te=interp2(xv,xv,T,XE,YE,'linear',0);
te=sqrt(max(Te,0)); te=te/norm(te(:));
% Comprobar adjunto antes de optimizar (identidad de producto interno).
rng(7); a=A.*exp(2i*pi*rand(N)); b=exp(2i*pi*rand(N));
la=L(a); lb=Lt(b); lhs=sum(conj(la(:)).*b(:)); rhs=sum(conj(a(:)).*lb(:));
adjErr=abs(lhs-rhs)/max([abs(lhs),abs(rhs),eps]);
assert(adjErr<1e-8,'Fallo prueba de adjunto.');
clear a b la lb
%% 4. Diseno de fase: Gerchberg-Saxton y optimizacion con pupila
rng(7); U=A.*exp(2i*pi*rand(N)); err=zeros(40,1);
for it=1:40
    if mod(it,10)==0, fprintf('Diseno convencional %d/40\n',it); end
    V=F(U); err(it)=norm(abs(V(:))-target(:))/norm(target(:));
    U=A.*exp(1i*angle(B(target.*exp(1i*angle(V)))));
end
U0=U; Ve0=L(U0); eye0=pupil.*P(U0);
capture0=sum(abs(eye0(:)).^2)/sum(abs(U0(:)).^2);
score0=scoreAmplitude(Ve0,te);
% Optimizacion por gradiente proyectado: maximizar solapamiento de
% amplitud absoluto, normalizado por potencia emitida (constante).
% Esto favorece tanto captura como fidelidad, no solo una imagen normalizada.
% Regularizacion opcional: mantener reconstruccion desde apertura completa.
weightHOE=0.15; th=target/norm(target(:));
history=zeros(40,3); best=U; bestScore=-inf;
for it=1:40
    if mod(it,10)==0, fprintf('Optimizacion con pupila %d/40\n',it); end
    v=L(U); vh=F(U);
    s=sum(te(:).*abs(v(:))); sh=sum(th(:).*abs(vh(:)));
    objective=(s^2+weightHOE*sh^2)/sum(abs(U(:)).^2);
    if objective>bestScore, bestScore=objective; best=U; end
    g=Lt(s*te.*exp(1i*angle(v)))+weightHOE*B(sh*th.*exp(1i*angle(vh)));
    U=A.*exp(1i*angle(g+0.1*U));
    history(it,:)=[objective,scoreAmplitude(v,te),sum(abs(v(:)).^2)/sum(abs(U(:)).^2)];
end
v=L(U); vh=F(U); s=sum(te(:).*abs(v(:))); sh=sum(th(:).*abs(vh(:)));
if (s^2+weightHOE*sh^2)/sum(abs(U(:)).^2)>bestScore, best=U; end
U=best; phi=angle(U); V=F(U); I=abs(V).^2; I=I/max(I(:));
Ueye=pupil.*P(U); capture=sum(abs(Ueye(:)).^2)/sum(abs(U(:)).^2);
Ve=Fe(Ueye); Ie=abs(Ve).^2; Ie=Ie/max(max(Ie(:)),eps);
Ie0=abs(Ve0).^2; Ie0=Ie0/max(max(Ie0(:)),eps);
score1=scoreAmplitude(Ve,te);
%% 5. Validacion de profundidad: comparar en coordenadas fisicas
depths=1.4:0.05:2.6; axial=zeros(size(depths)); axial0=axial;
for j=1:numel(depths)
    z=depths(j); zz=lambda*z/W; xxz=(-N/2:N/2-1)*zz;
    [XZ,YZ]=meshgrid(xxz);
    tz=sqrt(max(0,interp2(xv,xv,T,XZ,YZ,'linear',0)));
    tz=tz/max(norm(tz(:)),eps);
    chirp=exp(-1i*k*(X.^2+Y.^2)/(2*z));
    vz=fftshift(ifft2(ifftshift(Ueye.*chirp)))*N;
    vz0=fftshift(ifft2(ifftshift(eye0.*chirp)))*N;
    axial(j)=scoreAmplitude(vz,tz); axial0(j)=scoreAmplitude(vz0,tz);
end
[~,peak]=max(axial); bestDepth=depths(peak);
fprintf('Prueba adjunto: %.3g\n',adjErr);
fprintf('Captura R0 / R1: %.4f / %.4f %%\n',100*capture0,100*capture);
fprintf('Similitud amplitud pupila R0 / R1: %.4f / %.4f\n',score0,score1);
fprintf('Profundidad de maxima similitud en barrido: %.3f m\n',bestDepth);
%% 6. Curvatura y compensacion ideal de fase
h=coeff(2)*X/radius+coeff(3)*Y/radius+ ...
    coeff(4)*(X/radius).^2+coeff(5)*X.*Y/radius^2+coeff(6)*(Y/radius).^2;
delta=-k*(1+cos(theta))*h;
phiComp=angle(exp(1i*(phi-delta)));
phaseStep=max([max(abs(diff(delta,1,1)),[],'all'),max(abs(diff(delta,1,2)),[],'all')]);
% Compensacion ideal: phiComp+delta=phi (mod 2pi). No simular perdidas
% del caso sin compensar si el muestreo de fase no es suficiente.
curvedValid=phaseStep<=pi;
if curvedValid
    Uc=U.*exp(1i*delta); vc=L(Uc); Ic=abs(vc).^2; Ic=Ic/max(max(Ic(:)),eps);
    curvedScore=scoreAmplitude(vc,te);
    ec=pupil.*P(Uc); curvedCapture=sum(abs(ec(:)).^2)/sum(abs(U(:)).^2);
else
    Ic=[]; curvedScore=NaN; curvedCapture=NaN;
end
% Verificar que la aproximacion compensada coincide con el plano.
compResidual=norm(A.*exp(1i*(phiComp+delta))-U,'fro')/norm(U,'fro');
fprintf('Ancho objetivo %.1f mm; apertura HOE %.1f mm\n',logoW*1000,W*1000);
fprintf('Limite geometrico sin vineteo %.2f mm\n',(W-2*pupilRadius*D/De)*De/eyeDist*1000);
fprintf('Paso fase curva %.3f rad; caso curvo cuantitativo: %d\n',phaseStep,curvedValid);
fprintf('Residuo compensacion ideal %.3g\n',compResidual);
if ~curvedValid, warning('Curvatura sin compensar submuestreada; no se reporta perdida cuantitativa.'); end

%% 7. Guardar datos antes de dibujar las figuras
save(fullfile(out,'datos_HOE.mat'),'lambda','W','N','dx','D','eyeDist','theta','logoW','logoH','x','xv','xe','T','Te','phi','phiComp','I','Ie','Ie0','err','history','capture0','capture','score0','score1','depths','axial','axial0','bestDepth','adjErr','R','center','phaseStep','curvedValid','curvedScore','curvedCapture','compResidual','coeff','radius','fitRMS','qfit','h','Ic','pupilRadius');
% Liberar matrices de trabajo grandes antes de abrir graficos.
clear F B P Pt L Lt Fe Be Hprop C1 C2 Ceye Ceout U U0 Ueye eye0 best g v vh V Ve Ve0 XV YV XE YE XZ YZ tz chirp vz vz0 target th te pupil A delta
%% 8. Resultados opticos y visualizacion ligera
limx=logoW*600; limy=max(logoH*650,logoW*200);
fig=figure('Color','w','Name','HOE Final - Reconstruccion','Position',[50 50 1200 700]);
subplot(2,3,1); showImage(xe,Te,limx,limy,'Logo objetivo');
subplot(2,3,2); showImage(xe,Ie0,limx,limy,'Sin optimizar para pupila');
subplot(2,3,3); showImage(xe,Ie,limx,limy,'Optimizado para pupila');
subplot(2,3,4); imagesc(x*1e3,x*1e3,phi); axis image xy; title('Fase plana (rad)'); colorbar;
subplot(2,3,5); imagesc(x*1e3,x*1e3,phiComp); axis image xy; title('Fase compensada ideal (rad)'); colorbar;
subplot(2,3,6); plot(depths,axial0,'--',depths,axial,'LineWidth',1.5); grid on; xlabel('Distancia (m)'); ylabel('Similitud de amplitud'); legend('Convencional','Optimizado','Location','best'); title('Barrido axial');
colormap(fig,gray); saveas(fig,fullfile(out,'reconstruccion.png'));
fig=figure('Color','w','Name','HOE Final - Profundidad');
plot(depths,axial0,'--',depths,axial,'LineWidth',1.8); hold on; xline(D+eyeDist,':','Objetivo 2 m'); grid on;
xlabel('Distancia desde el ojo (m)'); ylabel('Similitud de amplitud'); legend('Convencional','Optimizado','Objetivo','Location','best');
title(sprintf('Maximo de similitud: %.2f m',bestDepth)); saveas(fig,fullfile(out,'profundidad.png'));
fig=figure('Color','w','Name','HOE Final - Curvatura','Position',[50 50 1000 500]);
subplot(1,2,1); surf(X(1:32:end,1:32:end)*1e3,Y(1:32:end,1:32:end)*1e3,h(1:32:end,1:32:end)*1e6,'EdgeColor','none'); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('Altura (um)'); title('Ajuste local del CAD'); colorbar;
subplot(1,2,2);
if curvedValid
    showImage(xe,Ic,limx,limy,'Curvo sin compensacion');
else
    axis off; text(0,.6,{'Caso sin compensar no cuantitativo','Muestreo insuficiente de fase','Compensacion ideal por construccion'});
end
saveas(fig,fullfile(out,'curvatura.png'));
% Usar los datos guardados: si falla el dibujo, no se pierde el calculo.
visualizarFinal(out,Pcad,faces);
fprintf('Resultados guardados en %s\n',out);
fprintf('Modelo escalar ideal: sin refraccion multicapa ni eficiencia volumetrica.\n');
end

%% Funciones auxiliares
function [R,center,coeff,radius,fitRMS,qfit]=ajusteCAD(V,faces,W)
% Misma seleccion de cara que R1, mismo ajuste cuadratico de R2.
cent=(V(faces(:,1),:)+V(faces(:,2),:)+V(faces(:,3),:))/3;
targetCAD=(min(V)+max(V))/2;
[~,j]=min(sum((cent-targetCAD).^2,2)); center=cent(j,:);
n=cross(V(faces(j,2),:)-V(faces(j,1),:),V(faces(j,3),:)-V(faces(j,1),:));
assert(norm(n)>0,'Triangulo degenerado.'); n=n/norm(n);
ref=[0 0 1]; if abs(dot(ref,n))>0.95, ref=[0 1 0]; end
u=cross(ref,n); u=u/norm(u); v=cross(n,u); R=[u;v;n];
% Verificar visualmente que la normal apunta al lado del habitaculo.
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
end
function v=pupilProp(u,H,idx,M,N,adjoint)
p=zeros(M,'like',u);
if adjoint
    p(idx,idx)=u(idx,idx);
    p=ifft2(fft2(p).*conj(H)); v=p(1:N,1:N);
else
    p(1:N,1:N)=u;
    p=ifft2(fft2(p).*H); v=zeros(N,'like',u); v(idx,idx)=p(idx,idx);
end
end
function s=scoreAmplitude(v,t)
a=abs(v(:)); b=t(:); s=sum(a.*b)/max(norm(a)*norm(b),eps);
end
function showImage(x,I,lx,ly,lab)
imagesc(x*1e3,x*1e3,I); axis image xy; xlim([-lx lx]); ylim([-ly ly]); title(lab); xlabel('mm'); ylabel('mm');
end

function visualizarFinal(out,V,F)
s=load(fullfile(out,'datos_HOE.mat'),'R','center','W','D','eyeDist','theta','logoW','logoH','xe','Ie','capture','score1','bestDepth');
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
fig=figure('Color','w','Name','HOE Final - Prototipo','Position',[80 80 1200 650]);
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
title(ax,'Sistema a escala real (CAD opaco para dibujo)'); view(ax,[1.3 0.8 0.65]*s.R);
axis(ax,'tight'); camproj(ax,'orthographic');
ax2=subplot(1,2,2);
% Vista optica numerica del observador, no camara 3D ni retina simulada.
rgb=cat(3,.08*G,G,.15*G);
imagesc(ax2,xx*1e3,yy*1e3,rgb); axis(ax2,'image'); axis(ax2,'xy');
xlabel(ax2,'x virtual (mm)'); ylabel(ax2,'y virtual (mm)');
title(ax2,{sprintf('Logo calculado desde pupila: %.0f mm',s.logoW*1e3),sprintf('Captura %.2f %% | similitud %.3f | pico %.2f m',100*s.capture,s.score1,s.bestDepth)});
drawnow;
saveas(fig,fullfile(out,'prototipo_3D.png')); savefig(fig,fullfile(out,'prototipo_3D.fig'));
figLogo=figure('Color','w','Name','HOE Final - Logo observado');
axLogo=axes(figLogo); imagesc(axLogo,xx*1e3,yy*1e3,rgb); axis(axLogo,'image'); axis(axLogo,'xy');
xlabel(axLogo,'x virtual (mm)'); ylabel(axLogo,'y virtual (mm)');
title(axLogo,'Reconstruccion desde pupila a 2 m (intensidad normalizada)');
saveas(figLogo,fullfile(out,'logo_observado.png'));
% Vista global adicional en alambre, con logo fisico, sin oclusion opaca.
fig2=figure('Color','w','Name','HOE Final - Vista en alambre');
ax=axes(fig2); hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
patch(ax,'Faces',F,'Vertices',V,'FaceColor','none','EdgeColor',[.7 .8 .85]);
scatter3(ax,p(:,1),p(:,2),p(:,3),5,colors,'filled');
ray(ax,S,H,[.85 .5 .1],'-'); ray(ax,H,E,[.1 .3 .85],'-'); ray(ax,H,Q,[.1 .3 .85],'--');
text(ax,E(1),E(2),E(3),' Ojo'); text(ax,Q(1),Q(2),Q(3),' Imagen virtual a 2 m');
xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)'); title(ax,'Geometria y logo virtual a escala real');
view(ax,[1.3 0.8 0.65]*s.R); axis(ax,'tight'); camproj(ax,'orthographic'); drawnow;
saveas(fig2,fullfile(out,'geometria_3D.png')); savefig(fig2,fullfile(out,'geometria_3D.fig'));
end
function ray(ax,a,b,c,style)
plot3(ax,[a(1) b(1)],[a(2) b(2)],[a(3) b(3)],style,'Color',c,'LineWidth',1.5);
end
