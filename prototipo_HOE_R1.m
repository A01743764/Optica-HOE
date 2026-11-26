function prototipo_HOE_R1(cadFile,logoFile)
% Prototipo HOE reflectivo escalar. MATLAB R2018b o posterior.
% prototipo_HOE('parabrisas.stl','logo.png')
% Sin argumentos: busca archivos junto al script; si faltan, usa demo.
root=fileparts(mfilename('fullpath'));
if nargin<1, cadFile=fullfile(root,'parabrisas.stl'); end
if nargin<2, logoFile=fullfile(root,'logo.png'); end
out=fullfile(root,'resultados_HOE_R1'); if ~exist(out,'dir'), mkdir(out); end
%% PARAMETROS (SI, salvo STL en mm)
lambda=532e-9; k=2*pi/lambda; N=1024; W=0.020;
dx=W/N; x=(-N/2:N/2-1)*dx; [X,Y]=meshgrid(x);
D=1.4; eyeDist=0.6; pupilRadius=0.0025; theta=30*pi/180;
% +z local hacia el ojo; imagen virtual z=-D. Incidencia hacia -z.
% Portadora oblicua se trata analiticamente, no se muestrea en la FFT.
fx=(-N/2:N/2-1)/(N*dx); [FX,FY]=meshgrid(fx);
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
% Mantener proporcion, ancho objetivo 6 mm, centrado en plano virtual.
logoW=min(0.006,0.75*N*dv); logoH=logoW*size(im,1)/size(im,2);
if logoH>0.75*N*dv, logoW=logoW*(0.75*N*dv/logoH); logoH=0.75*N*dv; end
% interp2 requiere imagen de al menos 2x2.
T=interp2(linspace(-logoW/2,logoW/2,size(im,2)), ...
    linspace(-logoH/2,logoH/2,size(im,1)),flipud(im),XV,YV,'nearest',0);
T=max(T,0); T=T/max(T(:)); target=sqrt(T); target=target*norm(A(:))/norm(target(:));
%% R1: operador HOE -> pupila -> plano virtual y su adjunto
% Propagacion con padding, recorte central y pupila. Todos los operadores
% discretos usan FFT unitarias. El adjunto del recorte es zero padding.
M=2*N; f=(-M/2:M/2-1)/(M*dx); [KX,KY]=meshgrid(2*pi*f);
KZ=sqrt(max(0,k^2-KX.^2-KY.^2)); idx=N/2+(1:N);
% Limite de muestreo del propagador (conservador, separable).
flim=1/(lambda*sqrt(1+(2*eyeDist/(M*dx))^2));
Hprop=exp(1i*(KZ-k)*eyeDist).*double(abs(KX)<=2*pi*flim & abs(KY)<=2*pi*flim);
clear KX KY KZ FX FY
pupil=double(X.^2+Y.^2<=pupilRadius^2);
De=D+eyeDist; dve=lambda*De/W; xe=(-N/2:N/2-1)*dve; [XE,YE]=meshgrid(xe);
Ceye=exp(-1i*k*(X.^2+Y.^2)/(2*De));
Ceout=exp(-1i*k*(XE.^2+YE.^2)/(2*De));
Fe=@(u) Ceout.*(fftshift(ifft2(ifftshift(Ceye.*u)))*N);
Be=@(u) conj(Ceye).*(fftshift(fft2(ifftshift(conj(Ceout).*u)))/N);
P=@(u) propCrop(u,Hprop,idx,M);
Pt=@(u) propCrop(u,conj(Hprop),idx,M);
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
%% VALIDACION INDEPENDIENTE: punto virtual, ajuste de rayos divergentes
% Se usa frente esferico ideal de referencia, no rayos extraidos del logo.
[xx,yy]=meshgrid(linspace(-0.009,0.009,9)); H=[xx(:),yy(:),zeros(numel(xx),1)];
Vpoint=[0 0 -D]; E=[0 0 eyeDist]; dirs=H-Vpoint; dirs=dirs./sqrt(sum(dirs.^2,2));
J=zeros(3); rhs=zeros(3,1);
for j=1:size(H,1), d=dirs(j,:)'; Q=eye(3)-d*d'; J=J+Q; rhs=rhs+Q*H(j,:)'; end
Vfit=(J\rhs)'; dist=norm(Vfit-E);
%% FIGURAS
fig=figure('Name','HOE R1','Color','w');
subplot(2,3,1); imagesc(xe*1e3,xe*1e3,Te); axis image xy; xlim([-4 4]); ylim([-3 3]); title('Objetivo');
subplot(2,3,2); imagesc(xe*1e3,xe*1e3,Ie0); axis image xy; xlim([-4 4]); ylim([-3 3]); title('Desde ojo R0');
subplot(2,3,3); imagesc(xe*1e3,xe*1e3,Ie); axis image xy; xlim([-4 4]); ylim([-3 3]); title('Desde ojo R1');
subplot(2,3,4); imagesc(x*1e3,x*1e3,phi); axis image xy; title('Fase R1'); colorbar;
subplot(2,3,5); imagesc(xv*1e3,xv*1e3,I); axis image xy; title('Desde HOE R1');
subplot(2,3,6); plot(depths,axial0,'--',depths,axial,'LineWidth',1.5); grid on; xlabel('Profundidad (m)'); ylabel('Similitud amplitud'); legend('R0','R1');
colormap(gray); saveas(fig,fullfile(out,'reconstruccion.png'));
fig3=figure('Name','Optimizacion R1','Color','w');
subplot(2,1,1); plot(history(:,1)); ylabel('Objetivo'); grid on;
subplot(2,1,2); plot(history(:,2:3)); legend('Similitud','Captura'); xlabel('Iteracion'); grid on;
saveas(fig3,fullfile(out,'optimizacion.png'));
%% CAD: plano tangente local. Elegir una cara del vidrio solamente en FreeCAD.
% STL no tiene unidades; aqui se asume mm. normalFlip cambia lado del ojo.
normalFlip=1; cadScale=1e-3; cadUsed=false; R=eye(3); center=[0 0 0];
fig2=figure('Name','Integracion geometrica','Color','w'); hold on; axis equal; grid on; view(3);
if exist(cadFile,'file')
    tr=stlread(cadFile); Pcad=tr.Points*cadScale; faces=tr.ConnectivityList;
    cent=(Pcad(faces(:,1),:)+Pcad(faces(:,2),:)+Pcad(faces(:,3),:))/3;
    % Cara cercana al centro de la caja; ajustar manualmente si no es vidrio.
    targetCAD=(min(Pcad)+max(Pcad))/2;
    [~,j]=min(sum((cent-targetCAD).^2,2)); center=cent(j,:);
    n=cross(Pcad(faces(j,2),:)-Pcad(faces(j,1),:),Pcad(faces(j,3),:)-Pcad(faces(j,1),:));
    assert(norm(n)>0,'Triangulo degenerado.'); n=normalFlip*n/norm(n);
    ref=[0 0 1]; if abs(dot(ref,n))>0.95, ref=[0 1 0]; end
    u=cross(ref,n); u=u/norm(u); v=cross(n,u); R=[u;v;n];
    patch('Faces',faces,'Vertices',Pcad,'FaceColor',[.3 .7 .9],'FaceAlpha',.25,'EdgeColor','none');
    cadUsed=true;
else
    [a,b]=meshgrid(linspace(-.6,.6,35),linspace(-.35,.35,25));
    surf(a,b,.08*a.^2,'FaceColor',[.3 .7 .9],'FaceAlpha',.25,'EdgeColor','none');
end
world=@(p) p*R+center;
Q=world([-W/2 -W/2 0;W/2 -W/2 0;W/2 W/2 0;-W/2 W/2 0]);
patch(Q(:,1),Q(:,2),Q(:,3),[.1 .8 .3],'FaceAlpha',.8);
% Fuente colimada: dibujar rayos paralelos de incidencia, no fuente puntual.
sdir=[sin(theta) 0 -cos(theta)]; S=H-0.5*sdir;
for j=1:8:size(H,1)
    q=world([S(j,:);H(j,:)]); plot3(q(:,1),q(:,2),q(:,3),'Color',[.85 .5 .1]);
    q=world([H(j,:);H(j,:)+dirs(j,:)*eyeDist]); plot3(q(:,1),q(:,2),q(:,3),'b');
    q=world([H(j,:);Vpoint]); plot3(q(:,1),q(:,2),q(:,3),'b--');
end
q=world(E); plot3(q(1),q(2),q(3),'ko','MarkerFaceColor','k'); text(q(1),q(2),q(3),' Ojo');
% Logo en plano virtual: visualizacion de la intensidad reconstruida,
% no una pantalla fisica ni trazado adicional.
step=4; [lx,ly]=meshgrid(xe(1:step:end)); li=Ie(1:step:end,1:step:end);
sel=li>0.15 & abs(lx)<0.004 & abs(ly)<0.003;
qlogo=world([lx(sel),ly(sel),-D*ones(nnz(sel),1)]);
scatter3(qlogo(:,1),qlogo(:,2),qlogo(:,3),8,[0 .75 .25],'filled');
q=world(Vpoint); plot3(q(1),q(2),q(3),'go','MarkerFaceColor','g'); text(q(1),q(2),q(3),' Imagen virtual');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('HOE tangente: rayos de referencia, no trazado dentro del vidrio');
camlight; lighting gouraud; saveas(fig2,fullfile(out,'integracion_CAD.png'));
%% GUARDAR RESULTADOS
fprintf('Distancia virtual de referencia: %.6f m\n',dist);
fprintf('Fraccion de potencia de envolvente captada por pupila: %.6g\n',capture);
fprintf('Periodo portadora: %.3f um\n',lambda/sin(theta)*1e6);
fprintf('CAD real cargado: %d\n',cadUsed);
save(fullfile(out,'datos_HOE.mat'),'lambda','W','dx','D','eyeDist','theta','phi','T','I','Ie','x','xv','xe','err','capture','dist','Vfit','R','center','cadUsed','Te','Ie0','capture0','score0','score1','depths','axial','axial0','bestDepth','history','adjErr');
end

function v=propCrop(u,H,idx,M)
p=zeros(M,'like',u); p(idx,idx)=u;
p=fftshift(ifft2(ifftshift(fftshift(fft2(ifftshift(p))).*H)));
v=p(idx,idx);
end
function s=scoreAmplitude(v,t)
a=abs(v(:)); b=t(:);
s=sum(a.*b)/max(norm(a)*norm(b),eps);
end
