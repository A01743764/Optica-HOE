function prototipo_HOE_R3(logoMM,cadFile,logoFile)
% Prototipo HOE reflectivo escalar. MATLAB R2018b o posterior.
% prototipo_HOE('parabrisas.stl','logo.png')
% Sin argumentos: busca archivos junto al script; si faltan, usa demo.
root=fileparts(mfilename('fullpath'));
if nargin<1, logoMM=50; end
assert(isscalar(logoMM)&&isfinite(logoMM)&&logoMM>=6&&logoMM<=80,'Ancho entre 6 y 80 mm.');
if nargin<2, cadFile=fullfile(root,'parabrisas.stl'); end
if nargin<3, logoFile=fullfile(root,'logo.png'); end
out=fullfile(root,sprintf('resultados_HOE_R3_%gmm',logoMM)); if ~exist(out,'dir'), mkdir(out); end
%% PARAMETROS (SI, salvo STL en mm)
lambda=532e-9; k=2*pi/lambda; W=0.020; N=2048; if logoMM>55, N=3072; end
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
%% LOGO: imagen blanca sobre fondo negro; demo tipografica si falta PNG
if exist(logoFile,'file')
    [im,map,alpha]=imread(logoFile);
    if ~isempty(map), im=reshape(map(double(im)+1,:),[size(im,1),size(im,2),3]); end
    if isinteger(im), im=double(im)/double(intmax(class(im))); else, im=double(im); end
    if size(im,3)>1, im=mean(im(:,:,1:3),3); end
    if ~isempty(alpha), im=im.*(double(alpha)/max(1,max(double(alpha(:))))); end
else
    % Letras de prueba, no identidad grafica oficial de Vitro.
    letters={'10001;10001;10001;10001;01010;01010;00100', ...
        '111;010;010;010;010;010;111', ...
        '11111;00100;00100;00100;00100;00100;00100', ...
        '11110;10001;10001;11110;10100;10010;10001', ...
        '01110;10001;10001;10001;10001;10001;01110'};
    im=[]; for j=1:numel(letters), q=char(strsplit(letters{j},';'))-'0'; im=[im q zeros(7,1)]; end %#ok<AGROW>
end
assert(max(im(:))>0,'Logo vacio.');
% Recortar margenes negros para que ancho solicitado sea el ancho del contenido.
mask=im>0.05*max(im(:)); im=im(any(mask,2),any(mask,1));
assert(all(size(im)>=2),'Logo demasiado pequeno.');
% Mantener proporcion, ancho objetivo 6 mm, centrado en plano virtual.
logoW=logoMM/1000; assert(logoW<=0.75*N*dv,'Plano virtual insuficiente.'); logoH=logoW*size(im,1)/size(im,2);
if logoH>0.75*N*dv, logoW=logoW*(0.75*N*dv/logoH); logoH=0.75*N*dv; end
% interp2 requiere imagen de al menos 2x2.
T=interp2(linspace(-logoW/2,logoW/2,size(im,2)), ...
    linspace(-logoH/2,logoH/2,size(im,1)),flipud(im),XV,YV,'nearest',0);
T=max(T,0); T=T/max(T(:)); target=sqrt(T); target=target*norm(A(:))/norm(target(:));
%% R1: operador HOE -> pupila -> plano virtual y su adjunto
% Propagacion con padding, recorte central y pupila. Todos los operadores
% discretos usan FFT unitarias. El adjunto del recorte es zero padding.
% Fresnel de convolucion lineal exacta en la malla discreta:
% solo se calcula la region de pupila, evitando padding 2N por eje.
% Valido en aproximacion paraxial; fase global omitida.
pidx=find(abs(x)<=pupilRadius); M=2^nextpow2(N+numel(pidx)-1);
% Para N=2048, M=4096 seria innecesario: usar longitud suave menor.
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
%% Comparacion R0 reproducida con el mismo logo y propagador
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
%% Barrido axial: comparar en coordenadas fisicas, no pixeles reescalados
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
%% Geometria CAD y visualizacion
assert(exist(cadFile,'file')==2,'Se requiere el STL del parabrisas.');
tr=stlread(cadFile); Pcad=tr.Points*1e-3; faces=tr.ConnectivityList;
% Reusar la orientacion de R2, calculada con el mismo STL.
base2=fullfile(root,'resultados_HOE_R2','datos_HOE.mat');
assert(exist(base2,'file')==2,'Conservar resultados_HOE_R2/datos_HOE.mat del mismo CAD.');
c2=load(base2,'R','center','h','x','coeff','radius','fitRMS','phaseStep');
R=c2.R; center=c2.center;
h=c2.coeff(2)*X/c2.radius+c2.coeff(3)*Y/c2.radius+ ...
    c2.coeff(4)*(X/c2.radius).^2+c2.coeff(5)*X.*Y/c2.radius^2+c2.coeff(6)*(Y/c2.radius).^2;
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
if logoMM>55, warning('Logo superior al campo geometrico sin vineteo. Evaluar bordes y captura.'); end
if ~curvedValid, warning('Curvatura sin compensar submuestreada; no se reporta perdida cuantitativa.'); end
%% Figuras
limx=logoW*600; limy=max(logoH*650,logoW*200);
fig=figure('Color','w','Name','R3 - Reconstruccion','Position',[50 50 1300 700]);
subplot(2,3,1); showImage(xe,Te,limx,limy,'Objetivo');
subplot(2,3,2); showImage(xe,Ie0,limx,limy,'Convencional desde ojo');
subplot(2,3,3); showImage(xe,Ie,limx,limy,'Optimizado desde ojo');
subplot(2,3,4); imagesc(x*1e3,x*1e3,phi); axis image xy; title('Fase HOE plano'); colorbar;
subplot(2,3,5); imagesc(x*1e3,x*1e3,phiComp); axis image xy; title('Fase compensada ideal'); colorbar;
subplot(2,3,6); plot(depths,axial0,'--',depths,axial,'LineWidth',1.5); grid on; xlabel('Distancia (m)'); ylabel('Similitud'); legend('Convencional','Optimizado'); title('Barrido axial');
colormap(fig,gray); saveas(fig,fullfile(out,'reconstruccion.png'));
fig=figure('Color','w','Name','R3 - Curvatura');
subplot(1,2,1); surf(X(1:32:end,1:32:end)*1e3,Y(1:32:end,1:32:end)*1e3,h(1:32:end,1:32:end)*1e6,'EdgeColor','none'); xlabel('x (mm)'); ylabel('y (mm)'); zlabel('Altura (um)'); title('Ajuste local CAD'); colorbar;
subplot(1,2,2);
if curvedValid, showImage(xe,Ic,limx,limy,'Sin compensacion'); else, axis off; text(0,.6,{'Caso sin compensar no cuantitativo','Muestreo insuficiente de fase','Compensacion: ideal por construccion'}); end
saveas(fig,fullfile(out,'curvatura.png'));
%% Vista 3D con textura calculada a escala real
world=@(p) p*R+center; E=world([0 0 eyeDist]); H=world([0 0 0]); Q=world([0 0 -D]); S=world(-.5*[sin(theta) 0 -cos(theta)]);
ix=find(abs(xe)<=limx/1000); iy=find(abs(xe)<=limy/1000);
stride=max(1,ceil(max(numel(ix),numel(iy))/400)); ix=ix(1:stride:end); iy=iy(1:stride:end);
[gx,gy]=meshgrid(xe(ix),xe(iy)); intens=Ie(iy,ix); bright=intens.^0.55;
rgb=cat(3,.08*bright,bright,.15*bright); alpha=bright; alpha(intens<.012)=0;
fig=figure('Color','w','Name','R3 - Prototipo final','Position',[50 50 1350 650]);
ax=subplot(1,2,1); hold(ax,'on'); glass(ax,Pcad,faces); axis(ax,'equal'); grid(ax,'on'); view(ax,3);
p=world([-W/2 -W/2 0;W/2 -W/2 0;W/2 W/2 0;-W/2 W/2 0]); patch(ax,p(:,1),p(:,2),p(:,3),[.1 .8 .3]);
line3(ax,S,H,[.85 .5 .1],'-'); line3(ax,H,E,[.1 .3 .85],'-'); line3(ax,H,Q,[.1 .3 .85],'--');
text(ax,E(1),E(2),E(3),' Ojo'); text(ax,H(1),H(2),H(3),' HOE'); text(ax,S(1),S(2),S(3),' Fuente'); text(ax,Q(1),Q(2),Q(3),' 2 m del ojo');
drawLogo(ax,world,gx,gy,-D,rgb,alpha); xlabel(ax,'X (m)'); ylabel(ax,'Y (m)'); zlabel(ax,'Z (m)'); title(ax,'Sistema a escala real'); camlight(ax,'headlight');
ax=subplot(1,2,2); hold(ax,'on'); glass(ax,Pcad,faces); drawLogo(ax,world,gx,gy,-D,rgb,alpha);
axis(ax,'equal'); axis(ax,'off'); campos(ax,E); camtarget(ax,Q); camup(ax,R(2,:)); camproj(ax,'perspective'); camva(ax,max(3,2*atan(logoW/De)*180/pi)); title(ax,'Vista desde ojo (zoom angular, textura numerica)'); camlight(ax,'headlight');
savefig(fig,fullfile(out,'prototipo.fig')); saveas(fig,fullfile(out,'prototipo.png'));
fig=figure('Color','k','Name','R3 - Logo observado'); ax=axes(fig); imagesc(ax,xe(ix)*1e3,xe(iy)*1e3,rgb); axis(ax,'image'); axis(ax,'xy'); set(ax,'Color','k','XColor','w','YColor','w'); title(ax,sprintf('Reconstruccion a 2 m | ancho %.0f mm',logoW*1000),'Color','w'); xlabel(ax,'mm'); ylabel(ax,'mm'); saveas(fig,fullfile(out,'logo_observado.png'));
save(fullfile(out,'datos_HOE.mat'),'lambda','W','N','dx','D','eyeDist','theta','logoW','logoH','x','xv','xe','T','Te','phi','phiComp','I','Ie','Ie0','err','history','capture0','capture','score0','score1','depths','axial','axial0','bestDepth','adjErr','R','center','phaseStep','curvedValid','curvedScore','curvedCapture','compResidual');
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
function glass(ax,V,F)
patch(ax,'Faces',F,'Vertices',V,'FaceColor',[.3 .7 .9],'FaceAlpha',.18,'EdgeColor','none');
end
function line3(ax,a,b,c,sty)
plot3(ax,[a(1) b(1)],[a(2) b(2)],[a(3) b(3)],sty,'Color',c,'LineWidth',1.6);
end
function drawLogo(ax,world,x,y,z,rgb,alpha)
p=world([x(:),y(:),z*ones(numel(x),1)]);
surf(ax,reshape(p(:,1),size(x)),reshape(p(:,2),size(x)),reshape(p(:,3),size(x)),rgb,'FaceColor','texturemap','EdgeColor','none','FaceAlpha','texturemap','AlphaData',alpha,'AlphaDataMapping','none','FaceLighting','none');
end
