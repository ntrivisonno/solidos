function [XYZ,CONEC,CONECEXT,CONECINT,NEIGEXT,NEIGINT] = GenPrismaticMeshHexa8(LX,LY,LZ,NX,NY,NZ)

  %------------------------------------------------
  % [XYZ,CONEC,CONECEXT,CONECINT,NEIGEXT,NEIGINT] = GenPrismaticMeshHexa8(LX,LY,LZ,NX,NY,NZ)
  % generador de malla hexahedrica p dominios prismaticos
  % LX,LY,LZ son las dimensiones del prisma
  % NX,NY,NZ -> la cant de divisiones de los segmentos
  % NEIGEXT -> elem vecinos a esas caras externa	 
  %       y |   z
  %         | /
  %         |/_ _ _ x
  %
  %-------------------------------------------------
	 
DX = LX/NX;
X = 0:DX:LX;
DY = LY/NY;
%Y = LY:-DY:0;
Y = 0:DY:LY;
DZ = LZ/NZ;
Z = 0:DZ:LZ;
% COORDENADAS NODALES
X = repmat(X(:),1,NY+1);
X = repmat(X(:),1,NZ+1);
Y = repmat(Y,NX+1,1);
Y = repmat(Y(:),1,NZ+1);
Z = repmat(Z,(NX+1)*(NY+1),1);
XYZ = [X(:) Y(:) Z(:)];

% CONECTIVIDAD DE LA PRIMERA CARA
CONECQ = [1 2 3+NX 2+NX];
% CONECTIVIDAD DEL ELEMENTO 1
CONEC1 = [CONECQ CONECQ+(NX+1)*(NY+1)];
NELE = NX*NY*NZ;
CONEC = zeros(NELE,8);
% CONECTIVIDAD DE LOS ELEMENTOS 1 a NX (SOBRE LA PRIMERA FRANJA PARALELA AL EJE X)
CONEC(1,:) = CONEC1;
for i=2:NX
    CONEC1 = CONEC1 + 1;
    CONEC(i,:) = CONEC1;
end
% CONECTIVIDAD DE LOS ELEMENTOS 1 a NX*NY (SOBRE LA PRIMERA CAPA PARALELA AL PLANO XY)
CONEC1 = CONEC(1:NX,:);
for i=2:NY
    CONEC1 = CONEC1+NX+1;
    CONEC(NX*(i-1)+1:i*NX,:) = CONEC1;
end
% CONECTIVIDAD DE LOS ELEMENTOS SOBRE TODAS LAS RODAJAS
CONEC1 = CONEC(1:NX*NY,:);
for i=2:NZ
    CONEC1 = CONEC1+(NX+1)*(NY+1);
    CONEC(NX*NY*(i-1)+1:i*NX*NY,:) = CONEC1;
end

if nargout>2
    % ELEMENTOS SOBRE LA SUPERFICIE EXTERNA
    EZ0 = 1:NX*NY;
    EZM = NELE-NX*NY+1:NELE;
    EX0 = reshape(bsxfun(@plus,1:NX:NX*(NY-1)+1,(0:NX*NY:NELE-NX*NY)'),1,[]);
    EXM = EX0 + NX - 1;
    EY0 = reshape(bsxfun(@plus,1:NX,(0:NX*NY:NELE-NX*NY)'),1,[]);
    EYM = EY0 + NX*(NY-1);
    
    CONECX0 = CONEC(EX0,[1 4 8 5]);
    CONECXM = CONEC(EXM,[2 3 7 6]);
    CONECY0 = CONEC(EY0,[1 5 6 2]);
    CONECYM = CONEC(EYM,[4 8 7 3]);
    CONECZ0 = CONEC(EZ0,1:4);
    CONECZM = CONEC(EZM,5:8);
    CONECEXT = [CONECX0; CONECXM; CONECY0; CONECYM; CONECZ0; CONECZM];
    
    if nargout > 3
        NEIGEXT = [EX0 EXM EY0 EYM EZ0 EZM]';

        % ELEMENTOS DE INTERFACE
        CONECINT = [];
        NEIGINT = [];
        % PARALELOS AL PLANO XY
        if NZ>1
            NEIG1 = reshape(bsxfun(@plus,EZ0,(0:NX*NY:NX*NY*(NZ-2))'),[],1);
            CONECINT = [CONECINT; CONEC(NEIG1,5:8)];
            NEIG2 = NEIG1 + NX*NY; 
            NEIGINT = [NEIGINT; NEIG1 NEIG2];

        end

        % PARALELOS AL PLANO ZX
        if NY>1
            NEIG1 = reshape(bsxfun(@plus,EY0,(0:NX:NX*(NY-2))'),[],1);
            CONECINT = [CONECINT; CONEC(NEIG1,[4 8 7 3])];
            NEIG2 = NEIG1 + NX; 
            NEIGINT = [NEIGINT; NEIG1 NEIG2];
        end

        % PARALELOS AL PLANO YZ
        if NX>1
            NEIG1 = reshape(bsxfun(@plus,EX0,(0:NX-2)'),[],1);
            CONECINT = [CONECINT; CONEC(NEIG1,[2 3 7 6])];
            NEIG2 = NEIG1 + 1; 
            NEIGINT = [NEIGINT; NEIG1 NEIG2];
        end
    end
end
