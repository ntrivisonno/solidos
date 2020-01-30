function F = FuerzaInternaHexa8 (x,y,z,U,funcTension)
% funcTension -> tension en funcion de las relaciones constitutivas

  % desplazamientos de los nodos del elem
  UX = U(1:3:end);
  UY = U(2:3:end);
  UZ = U(3:3:end);
  
  npg = 8;
  nnodxelem = 8;
  ndofxelem = nnodxelem*3;
  F = zeros(ndofxelem,1);
  for i = 1:npg
    % peso de coord naturales del punto de integracion
    [w, psi, eta, ceda] = WeightAndCoordHexa(npg,i);
    [J,dN_dpsi,dN_deta,dN_dceda] = jacobiano (x,y,z,psi,eta,ceda);
    j = inv(J);
    dNdX = zeros ( nnodxelem,1);
    dNdY = zeros ( nnodxelem,1);
    dNdZ = zeros ( nnodxelem,1);
    for k1 = 1:nnodxelem
        dNknat = [dN_dpsi(k1);dN_deta(k1);dN_dceda(k1)];
        dNkglob = j*dNknat;
        dNdX(k1) = dNkglob(1);
        dNdY(k1) = dNkglob(2);
        dNdZ(k1) = dNkglob(3);
    end
    % Gradiente de desplazamientos
    D = zeros(3);
    D(1,1) = dot(dNdX,UX);
    D(1,2) = dot(dNdY,UX);
    D(1,3) = dot(dNdZ,UX);           
    D(2,1) = dot(dNdX,UY);
    D(2,2) = dot(dNdY,UY);
    D(2,3) = dot(dNdZ,UY);           
    D(3,1) = dot(dNdX,UZ);
    D(3,2) = dot(dNdY,UZ);
    D(3,3) = dot(dNdZ,UZ);
    % Gradiente de deformacion
    A = D +eye(3);
    % Deformacion Gree-Lagrange
    E = (A*A'-eye(3))/2;
    % en notacion de Voigt
    E = [E(1,1) E(2,2) E(3,3) 2*E(1,2) 2*E(2,3) 2*E(3,1)]';
    % segundo tensor de piola-kirchoff en notacion de Voigt
    T2  = funcTension(E);
    % matriz de gradiente de funciones de forma
    Blin = zeros(6,24);
    Blin(1,1:3:end) = dNdX;
    Blin(2,2:3:end) = dNdY;
    Blin(3,3:3:end) = dNdZ;
    Blin(4,1:3:end) = dNdY;
    Blin(4,2:3:end) = dNdX;
    Blin(5,2:3:end) = dNdZ;
    Blin(5,3:3:end) = dNdY;
    Blin(6,1:3:end) = dNdZ;
    Blin(6,3:3:end) = dNdX;

    Bno_lin = zeros(6,24);
    Bno_lin(1,1:3:end) = dNdX*D(1,1);
    Bno_lin(1,2:3:end) = dNdY*D(2,1);
    Bno_lin(1,3:3:end) = dNdZ*D(3,1);
    Bno_lin(2,1:3:end) = dNdX*D(1,2);
    Bno_lin(2,2:3:end) = dNdY*D(2,2);
    Bno_lin(2,3:3:end) = dNdZ*D(3,2);
    Bno_lin(3,1:3:end) = dNdX*D(1,3);
    Bno_lin(3,2:3:end) = dNdY*D(2,3);
    Bno_lin(3,3:3:end) = dNdZ*D(3,3);
    Bno_lin(4,1:3:end) = dNdY*D(1,1)+dNdX*D(1,2);
    Bno_lin(4,2:3:end) = dNdY*D(2,1)+dNdX*D(2,2);
    Bno_lin(4,3:3:end) = dNdY*D(3,1)+dNdX*D(3,2); 
    Bno_lin(5,2:3:end) = dNdZ*D(1,2)+dNdX*D(1,3); 
    Bno_lin(5,2:3:end) = dNdZ*D(2,2)+dNdX*D(2,3);
    Bno_lin(5,3:3:end) = dNdZ*D(3,2)+dNdX*D(3,3);
    Bno_lin(6,1:3:end) = dNdX*D(1,3)+dNdZ*D(1,1);
    Bno_lin(6,2:3:end) = dNdX*D(2,3)+dNdZ*D(2,1);
    Bno_lin(6,3:3:end) = dNdX*D(3,3)+dNdZ*D(3,1);
    
    B = Blin + Bno_lin;
    % integrando
    f = B'*T2;
    % integral numerica
    F = F + f*w*det(J);
end


end
