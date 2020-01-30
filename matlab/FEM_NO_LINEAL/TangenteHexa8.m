function K = TangenteHexa8(x,y,z,U,funcTension,funcMod)
  % funcTension -> tension en funcion de las relaciones constitutivas
	 
  % desplazamientos de los nodos del elem
  UX = U(1:3:end);
  UY = U(2:3:end);
  UZ = U(3:3:end);
  
  npg = 8;
  nnodxelem = 8;
  ndofxelem = nnodxelem*3;
  K = zeros(ndofxelem);
  
  % I3x3 = eye(3);
  % kgeo = zeros(ndofxelem);
  
  % Modulos elasticos
  if ~isa(funcMod,'function_handle') % en mat kirchoff-saintVenait no es una funcion y es una constante
    C = funcMod;
  end
  
  for i = 1:npg
    % peso de coord naturales del punto de integracion
    [w, psi, eta, ceda] = WeightAndCoordHexa(npg,i);
    [J,dN_dpsi,dN_deta,dN_dceda] = jacobiano (x,y,z,psi,eta,ceda);
    j = inv(J);
    dNdX = zeros (nnodxelem,1);
    dNdY = zeros (nnodxelem,1);
    dNdZ = zeros (nnodxelem,1);
    for k1 = 1:nnodxelem
      dNknat = [dN_dpsi(k1);dN_deta(k1);dN_dceda(k1)];
      dNkglob = j*dNknat;
      dNdX(k1) = dNkglob(1);
      dNdY(k1) = dNkglob(2);
      dNdZ(k1) = dNkglob(3);
    end
    % Gradiente de desplazamientos
    D=zeros(3,3);
    D(1,1)=dot(dNdX,UX);
    D(1,2)=dot(dNdY,UX);
    D(1,3)=dot(dNdZ,UX);
    D(2,1)=dot(dNdX,UY);
    D(2,2)=dot(dNdY,UY);
    D(2,3)=dot(dNdZ,UY);
    D(3,1)=dot(dNdX,UZ);
    D(3,2)=dot(dNdY,UZ);
    D(3,3)=dot(dNdZ,UZ);
    % Gradiente de deformacion
    A = D +eye(3);
    % Deformacion Gree-Lagrange
    E = (A*A'-eye(3))/2;
    % en notacion de Voigt
    E = [E(1,1) E(2,2) E(3,3) 2*E(1,2) 2*E(2,3) 2*E(3,1)]';
    % segundo tensor de piola-kirchoff en notacion de Voigt
    T2  = funcTension(E);
    % Segundo tensor de piola-kirchoff en forma tensorial
    T2T = [T2(1) T2(4) T2(6)
           T2(4) T2(2) T2(5)
           T2(6) T2(5) T2(3)];
    % matriz de gradiente de funciones de forma
    Blin=zeros(6,24);
    Blin(1,1:3:end)=dNdX;
    Blin(2,2:3:end)=dNdY;
    Blin(3,3:3:end)=dNdZ;
    Blin(4,1:3:end)=dNdY;
    Blin(4,2:3:end)=dNdX;
    Blin(5,2:3:end)=dNdZ;
    Blin(5,3:3:end)=dNdY;
    Blin(6,1:3:end)=dNdZ;
    Blin(6,3:3:end)=dNdX;
    
    Bno_lin = zeros(6,24);
    Bno_lin(1,1:3:end)=dNdX*D(1,1);
    Bno_lin(1,2:3:end)=dNdX*D(2,1);
    Bno_lin(1,3:3:end)=dNdX*D(3,1);
    Bno_lin(2,1:3:end)=dNdX*D(1,2);
    Bno_lin(2,2:3:end)=dNdX*D(2,2);
    Bno_lin(3,1:3:end)=dNdX*D(1,3);
    Bno_lin(3,2:3:end)=dNdX*D(2,3);
    Bno_lin(3,3:3:end)=dNdX*D(3,3);
    Bno_lin(4,1:3:end)=dNdY*D(1,1)+dNdX*D(1,2);
    Bno_lin(4,2:3:end)=dNdY*D(2,1)+dNdX*D(2,2);
    Bno_lin(4,3:3:end)=dNdY*D(3,1)+dNdX*D(3,2);
    Bno_lin(5,1:3:end)=dNdZ*D(1,2)+dNdY*D(1,3);
    Bno_lin(5,2:3:end)=dNdZ*D(2,2)+dNdY*D(2,3);
    Bno_lin(5,3:3:end)=dNdZ*D(3,2)+dNdY*D(3,3);
    Bno_lin(6,1:3:end)=dNdX*D(1,3)+dNdZ*D(1,1);
    Bno_lin(6,2:3:end)=dNdX*D(2,3)+dNdZ*D(2,1);
    Bno_lin(6,3:3:end)=dNdX*D(3,3)+dNdZ*D(3,1);
    
    B = Blin + Bno_lin;
    
    % Modulos elasticos
    if isa(funcMod,'function_handle') % en mat kirchoff-saintVenait no es una funcion y es una constante
      C = funcMod(E);
    end
    % integrando
    kmat = B'*C*B;
    I3x3 = eye(3);
    kgeo = zeros(ndofxelem);
    
    for p=1:nnodxelem
      for q=1:nnodxelem
        kgeo_pq = [dNdX(p) dNdY(p) dNdZ(p)]*T2T*[dNdX(q) dNdY(q) dNdZ(q)]';
        kgeo(3*p-2:3*p,3*q-2:3*q) = kgeo_pq*I3x3;
      end
    end
 
    % geo = [ k(1,1)*I k(1,2)*I k(1,3)*I k(1,4)*I k(1,5)*I k(1,6)*I k(1,7)*I k(1,7)*I k(1,8)*I
    %         k(2,1)*I k(2,2)*I k(2,3)*I k(2,4)*I k(2,5)*I k(2,6)*I k(2,7)*I k(2,7)*I k(2,8)*I
    %         k(3,1)*I k(3,2)*I k(3,3)*I k(3,4)*I k(3,5)*I k(3,6)*I k(3,7)*I k(3,7)*I k(3,8)*I
    %         k(4,1)*I k(4,2)*I k(4,3)*I k(4,4)*I k(4,5)*I k(4,6)*I k(4,7)*I k(4,7)*I k(4,8)*I
    %         k(5,1)*I k(5,2)*I k(5,3)*I k(5,4)*I k(5,5)*I k(5,6)*I k(5,7)*I k(5,7)*I k(5,8)*I
    %         k(6,1)*I k(6,2)*I k(6,3)*I k(6,4)*I k(6,5)*I k(6,6)*I k(6,7)*I k(6,7)*I k(6,8)*I
    %         k(7,1)*I k(7,2)*I k(7,3)*I k(7,4)*I k(7,5)*I k(7,6)*I k(7,7)*I k(7,7)*I k(7,8)*I
    %         k(8,1)*I k(8,2)*I k(8,3)*I k(8,4)*I k(8,5)*I k(8,6)*I k(8,7)*I k(8,7)*I k(8,8)*I];
    %     
    
    % integral numerica
    K = K + (kgeo+kmat)*w*det(J);
  end
  
  
end
