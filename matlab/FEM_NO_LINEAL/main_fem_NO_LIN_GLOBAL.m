clear %all;clc

% si quiere cambiar otro material, hay q cambiar funcTen y funcMod, ej para
% un neohookeano

% J = jacobiano (x,y,z,psi,eta,ceda,xnod)
 
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
E  = 1e6; % Young
nu = 0.3; % Poisson

Ptot = 100;

% Forma de computar la tangente: numerica, analitica o constante
% calcid = 1; % numerica
calcid = 2; % analitica
% calcid = 3; % rigidez lineal

% material elastico Kirchoff-SaintVenant;
C = Hooke_3D(E,nu);
funcMod = C; 
funcTension = @(E) C*E;  

% figure();PlotMesh(XYZ,CONEC,'HEXA8');axis equal;

% cargas unitarias en dir Y en las aristas x=LX y=0
nodcar = find ((XYZ(:,1)==LX) & (XYZ(:,2)==0));
% vector de cargas nodales
Fext = zeros(3*nnod,1);
Fext(3*nodcar-1) = 1 % carga unitaria en los nodos nodcar en direccion y, 3 por nodcar es el eje z
Fext = Fext*Ptot/sum(Fext);

% fijaciones: cara x=0 empotrada
nodfij = find(XYZ(:,1)==0)
% grados de libertad fijos
gdlfij = [3*(nodfij-1)+1;3*(nodfij-1)+2;3*(nodfij-1)+3];
ufij   = 0*gdlfij;
% grados de libertad libres
gdllib = setdiff(1:3*nnod,gdlfij); % dif entre el conjunt de todo los GDL y los fijos

% iniciacion de desplazamientos, es la semilla de la iteracion
U = zeros(3*nnod,1);
U(gdlfij) = ufij;

% loop de solucion del problema no lineal | IMP -> en la 1er iteracion (considerando el seed como u=0) vamos a tener la sol del problema lineal
maxiter = 100;
tol = 1e-4*norm(Fext); % se normaliza la tolerancia en funcion de los esfuerzos internos
tim0 = tic;
for iter = 1:maxiter
  tic
  iter
  Fint = FuerzaInternaHexa8Global(XYZ,CONEC,U,funcTension);
  Res = Fint - Fext;
  Res1 = Res(gdllib);
  normRes = norm(Res(gdllib))
  if normRes < tol
    break;
  end
  % calculo de matriz tangente
  K = TangenteHexa8Global (XYZ,CONEC,U,funcTension,funcMod,calcid);
  K1 = K(gdllib,gdllib);
  % actualizacion de desplazamientos
  DU = -K1\Res1;
  U(gdllib) = U(gdllib)+DU;
  if iter==1
    U0 = U;
  end
  toc
end
toc(tim0)
% desplazamientos


SOL = struct('Name','U','Type','PointData','NumberOfComponents',3,'Value',U);
SOL(2) = struct('Name','U0','Type','PointData','NumberOfComponents',3,'Value',U0);
SOL(3) = struct('Name','DU','Type','PointData','NumberOfComponents',3,'Value',U-U0);

mtl2vtu('ejem_NO_LIN.vtu',XYZ,struct('Connectivity',CONEC,'Type','HEXA8'),SOL)

% ploteo post-proces
% figure();PlotDefMesh(XYZ,CONEC,u,1,'HEXA8');axis equal; % plotea la malla
% deformada

% hay que hacer la rutina para calcular las fuerzas internas
% tmb hacer lo mismo q antes pero para no lineal
% q ser'ia B'(y??)  copiar la rutina de matriz de rigidez
