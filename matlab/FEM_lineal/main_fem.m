
clear %all;clc
xnod = [-1 -1 -1;
 1 -1 -1;
 1 1 -1;
 -1 1 -1;
 -1 -1 1;
 1 -1 1;
 1 1 1 ;
 -1 1 1];

%N = ffHex8 (-1,-1,1,xnod)
[N,a,b,c] = ffHex8 (rand,rand,rand); % ojo que desp cambie las funciones sin pasar xnod, 
%ya que lo definia en ffHex8

%J = jacobiano (x,y,z,psi,eta,ceda,xnod)
 
LX = 10;
LY = 2;
LZ = 1;
NX = 20;
NY = 4;
NZ = 2;

[XYZ,CONEC,CONECEXT,CONECINT,NEIGEXT,NEIGINT] = GenPrismaticMeshHexa8(LX,LY,LZ,NX,NY,NZ);

nelem = size(CONEC,1);
nnod  = size(XYZ,1);

% vamos a calc el volumen de los elem
V = zeros(nelem,1);

E = 1e6; % Young
nu = 0.3; % Poisson
C = Hooke_3D(E,nu);
K = zeros(3*nnod);

for e = 1:nelem
  % nodos del elem
  conec = CONEC(e,:);
  % coord globales de los nodos del elem
  x = XYZ(conec,1);
  y = XYZ(conec,2);
  z = XYZ(conec,3);
  k = rigidezLinealHexa8 (x,y,z,C); % matriz de rigidez elem
  glo = zeros(1,24);

  for j = 1:3
      glo(j:3:end) = 3*(conec-1)+j;
  end
  K(glo,glo) = K(glo,glo) + k; 
  % volumen del elem
  if 1 
   v = 0;
  npg = 8;
  for i = 1:npg
      [w, psi, eta, ceda] = WeightAndCoordHexa(npg,i);
      J = jacobiano (x,y,z,psi,eta,ceda);
      f = det(J);
      v = v + f*w;
  end
  end
  if 0 
  % calc de volumen mediante funcion
    v = volumen_hexa8(x,y,z,xnod);
  end
  V(e) = v; % para verificar hacemos el volumen y tiene q ser = a LX*LY*LZ

  end

%figure();PlotMesh(XYZ,CONEC,'HEXA8');axis equal;

% cargas unitarias en dir Y en las aristas x=LX y=0
nodcar = find ((XYZ(:,1)==LX) & (XYZ(:,2)==0));
% vector de cargas nodales
F = zeros(3*nnod,1);
F(3*(nodcar-1)+2) = 1 % carga unitaria en los nodos nodcar en direccion y

% fijaciones: cara x=0 empotrada
nodfij = find(XYZ(:,1)==0)
% grados de libertad fijos
gdlfij = [3*(nodfij-1)+1;3*(nodfij-1)+2;3*(nodfij-1)+3];
ufij   = 0*gdlfij;
% grados de libertad libres
gdllib = setdiff(1:3*nnod,gdlfij); % dif entre el conjunt de todo los GDL y los fijos

% desplazamientos
F1 = F(gdllib) - K(gdllib,gdlfij)*ufij; % ver T26, q operacion, los ifjos, los conozco y los paso al t'ermino derecho
K1 = K(gdllib,gdllib);
u = zeros(3*nnod,1);
u(gdlfij) = ufij;
u(gdllib) = K1\F1;

SOL = struct('Name','U','Type','PointData','NumberOfComponents',3,'Value',u);
SOL(2) = struct('Name','F','Type','PointData','NumberOfComponents',3,'Value',F);
mtl2vtu('ejem.vtu',XYZ,struct('Connectivity',CONEC,'Type','HEXA8'),SOL)


% hay que hacer la rutina para calcular las fuerzas internas
% tmb hacer lo mismo q antes pero para no lineal
% q ser'ia B'(y??)  copiar la rutina de matriz de rigidez
