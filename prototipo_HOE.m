function prototipo_HOE(cadFile,logoFile)
% Prototipo HOE reflectivo escalar. MATLAB R2018b o posterior.
% prototipo_HOE('parabrisas.stl','logo.png')
% Sin argumentos: busca archivos junto al script; si faltan, usa demo.
root=fileparts(mfilename('fullpath'));
if nargin<1, cadFile=fullfile(root,'parabrisas.stl'); end
if nargin<2, logoFile=fullfile(root,'logo.png'); end
out=fullfile(root,'resultados_HOE'); if ~exist(out,'dir'), mkdir(out); end
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
%% GERCHBERG-SAXTON: diseno de fase de la envolvente reflejada
rng(7); U=A.*exp(2i*pi*rand(N)); err=zeros(40,1);
for it=1:numel(err)
    V=F(U); vabs=abs(V); err(it)=norm(vabs(:)-target(:))/norm(target(:));
    U=A.*exp(1i*angle(B(target.*exp(1i*angle(V)))));
end
phi=angle(U); V=F(U); I=abs(V).^2; I=I/max(I(:));
% Fase fisica ideal en plano: phi - k*sin(theta)*x (mod 2pi).
% Periodo portadora: lambda/sin(theta). No exportar como mascara fabricable.
%% PROPAGACION AL OJO: espectro angular de envolvente, con zero padding
M=2*N; f=(-M/2:M/2-1)/(M*dx); [KX,KY]=meshgrid(2*pi*f);
KZ=sqrt(max(0,k^2-KX.^2-KY.^2)); idx=N/2+(1:N);
P=zeros(M); P(idx,idx)=U;
P=fftshift(ifft2(ifftshift(fftshift(fft2(ifftshift(P))).*exp(1i*(KZ-k)*eyeDist))));
Ueye=P(idx,idx); pupil=double(X.^2+Y.^2<=pupilRadius^2);
Ueye=Ueye.*pupil; capture=sum(abs(Ueye(:)).^2)/sum(abs(U(:)).^2);
% Reconstruccion inversa desde el ojo, no foco real en el exterior.
De=D+eyeDist; dve=lambda*De/W; xe=(-N/2:N/2-1)*dve; [XE,YE]=meshgrid(xe);
Ve=exp(-1i*k*(XE.^2+YE.^2)/(2*De)).* ...
    (fftshift(ifft2(ifftshift(Ueye.*exp(-1i*k*(X.^2+Y.^2)/(2*De)))))*N);
Ie=abs(Ve).^2; Ie=Ie/max(max(Ie(:)),eps);
%% VALIDACION INDEPENDIENTE: punto virtual, ajuste de rayos divergentes
% Se usa frente esferico ideal de referencia, no rayos extraidos del logo.
[xx,yy]=meshgrid(linspace(-0.009,0.009,9)); H=[xx(:),yy(:),zeros(numel(xx),1)];
Vpoint=[0 0 -D]; E=[0 0 eyeDist]; dirs=H-Vpoint; dirs=dirs./sqrt(sum(dirs.^2,2));
J=zeros(3); rhs=zeros(3,1);
for j=1:size(H,1), d=dirs(j,:)'; Q=eye(3)-d*d'; J=J+Q; rhs=rhs+Q*H(j,:)'; end
Vfit=(J\rhs)'; dist=norm(Vfit-E);
%% FIGURAS
fig=figure('Name','Reconstruccion holografica','Color','w');
subplot(2,3,1); imagesc(xv*1e3,xv*1e3,T); axis image xy; title('Objetivo');
subplot(2,3,2); imagesc(x*1e3,x*1e3,phi); axis image xy; title('Fase envolvente (rad)'); colorbar;
subplot(2,3,3); imagesc(xv*1e3,xv*1e3,I); axis image xy; title('Reconstruccion virtual HOE');
subplot(2,3,4); imagesc(x*1e3,x*1e3,abs(Ueye).^2); axis image xy; title('Campo en pupila');
subplot(2,3,5); imagesc(xe*1e3,xe*1e3,Ie); axis image xy; title('Reconstruccion desde ojo');
subplot(2,3,6); plot(err,'LineWidth',1.5); grid on; xlabel('Iteracion'); ylabel('Error amplitud');
colormap(gray); saveas(fig,fullfile(out,'reconstruccion.png'));
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
q=world(Vpoint); plot3(q(1),q(2),q(3),'go','MarkerFaceColor','g'); text(q(1),q(2),q(3),' Imagen virtual');
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)'); title('HOE tangente: rayos de referencia, no trazado dentro del vidrio');
camlight; lighting gouraud; saveas(fig2,fullfile(out,'integracion_CAD.png'));
%% GUARDAR RESULTADOS
fprintf('Distancia virtual de referencia: %.6f m\n',dist);
fprintf('Fraccion de potencia de envolvente captada por pupila: %.6g\n',capture);
fprintf('Periodo portadora: %.3f um\n',lambda/sin(theta)*1e6);
fprintf('CAD real cargado: %d\n',cadUsed);
save(fullfile(out,'datos_HOE.mat'),'lambda','W','dx','D','eyeDist','theta','phi','T','I','Ie','x','xv','xe','err','capture','dist','Vfit','R','center','cadUsed');
end
