
clear %all;clc


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

K = rigidezLinealHexa8Global (XYZ,CONEC,C);

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
F1 = F(gdllib) - K(gdllib,gdlfij)*ufij; % ver T26, q operacion, los ifjos, los ocnozco y los paso al t'ermino derecho
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
