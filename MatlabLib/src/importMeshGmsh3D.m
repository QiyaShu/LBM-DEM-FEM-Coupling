function [p, e, t] = importMeshGmsh3D(fileName)

% importMeshGmsh - Read external mesh in Gmsh ASCII or binary file format.
%
% This QuickerSim CFD Toolbox function opens a Gmsh *.msh file written in
% text format and imports the mesh to MATLAB environment. It only imports
% first order mesh.
%
% [p, e, t] = importMeshGmsh(fileName)
%
% Input arguments:
% fileName - name of the *.msh mesh file to be imported.
%
% Output arguments:
% p - array containig node coordinates. It is of size 2-by-nnodes for 2-D
%     case or 3-by-nnodes for 3-D case, where nnodes denotes total number
%     of nodes in the mesh. The matrix is structured as follows: each
%     column corresponds to one node in the mesh, row 1 contains
%     x-coordinates, row 2 - y-coordinates and row 3 (if applicable) -
%     z-coordinates. For example p(2,17) is equal to the y-coordinate of
%     17th node in the mesh.
% e - array which defines all mesh boundaries. Each column of that matrix
%     corresponds to a single edge element and it has 7 rows in case of a
%     2-D linear mesh and 8 rows in case of 2-D second order mesh (second
%     order mesh is created by the convertMeshToSecondOrder function). For
%     a single edge element which is stored in j-th column of
%     array e, e(1,j) contains id of the first node along this edge, e(2,j)
%     id of the second node along the edge, elements e(3,j) and e(4,j)
%     correspond to the parameter along the curve in the original
%     description of the (p,e,t) format of Matlab PDE Toolbox - however
%     this CFD Toolbox ignores these values completely and the present
%     function assigns arbitrary values to these two elements, element
%     e(5,j) contains id of the boundary to which this particular edge
%     element belongs (when reading mesh from Gmsh mesh generator e(5,j)
%     will be equal to the Physical Group id defined in Gmsh), e(6,j) is
%     equal to the id (or Physical Group id) of the mesh subdomain which
%     is lying on the left hand side of the edge (when travelling from node
%     1 to node 2) and element e(7,j) id of the subdomain on the right hand
%     side of the given edge. However, current version of the CFD Toolbox
%     does not support multiple subdomains and hence e(6,j) is always
%     assigned with 1 (there is guarantee that the physical domain is lying
%     on the left hand side of the edge) and e(7,j) is always assigned with
%     zero in this toolbox version. For a second order 2-D mesh, there is
%     one additional row of the e array and then element e(8,j) contains
%     global id of the node which has been placed in the middle of the
%     particular edge. For 3-D cases the e array is built as follows. For
%     linear mesh rows 1-3 contain ids of nodes defining the boundary
%     triangular element, row 4 denotes physical group id and rows 5 and 6
%     contain id of the domain on the left and right hand side of the
%     boundary triangle. For a second order case, the array is built
%     correspondingly except for the fact that rows 1-6 contain node ids
%     and the rest of information is shifted to rows 7-9.
% t - array of all finite elements in the mesh. Each column corresponds to
%     a single mesh element and the matrix is of size m-by-n, where
%     n denotes number of elements in mesh and m equals:
%       m = 4, for a 2-D linear mesh
%       m = 7, for a 2-D second order mesh
%       m = 5, for a 3-D linear mesh
%       m = 11, for a 3-D second order mesh.
%     The last row always contains id of the subdomain to which this
%     particular element belongs (however, multiple subdomains are not
%     supported in the current vesion of the CFD Toolbox), and the rows
%     from 1 to m-1 contain global ids of the nodes which constitute this
%     element. The local node numbering for 2-D is as follows:
%
%     2-D first order triangle      2-D second order triangle
%
%           3                                  3
%          / \                                / \
%         /   \                              6   5
%        /     \                            /     \
%       1-------2                          1---4---2
%
% Arrays p, e, t are defined identically to the mesh format in Matlab PDE
% Toolbox, thus allowing for interchangable use of the mesh between
% toolboxes. Some exceptions arise in case of second order 2-D mesh (which
% is not handled by Matlab PDE toolbox at all) and in this case additional
% rows in the t array are defined to contain information about second order
% nodes in the middle of each edge. Additionally, array e for a second
% order 2-D mesh is supplemented with an additional row (8-th row
% containing the id of the mid-edge node). One additional exception is also
% that the value of the parameter of the node along the curve (row 3 and 4
% in the original format of e array in Matlab PDE Toolbox) is not used in
% this CFD Toolbox and is filled in with arbitrary values, which may have
% nothing in common with the original idea contained in Matlab PDE Toolbox.
% For the 3-D case e structure in the PDE Toolbox and e array of this
% toolbox are different - however, function transformBoundaryToCfdFormat,
% which is available in full version of the QuickerSim CFD Toolbox allows
% transformation of the e structure of the PDE Toolbox to the format
% understandable by QuickerSim CFD Toolbox.
%
% Examples:
%       [p, e, t] = importMeshGmsh('backwardFacingStepMesh.msh');
%
% Visit www.quickersim.com/cfd-toolbox-for-matlab/index for more info, help
% and support. Contact us by cfdtoolbox@quickersim.com
%
% See also: ASSEMBLENAVIERSTOKESMATRIX2D, ASSEMBLESTOKESMATRIX2D,
% CONVERTMESHTOSECONDORDER, EXPORTTOGMSH2D.

if ~exist(fileName)
   error('msh file with that name does not exist');
end

fd = fopen(fileName, 'r');

% txt = mgetl(fd, 4);
txt = fgets(fd);%, 4)
txt = fgets(fd);

%check if mesh is saved in ASCII or binary and set binaryflag
if str2double(txt(5)) == 0
    binaryflag = 0;
elseif str2double(txt(5)) == 1
    binaryflag = 1;
end

if binaryflag
    txt = fgets(fd);
    txt = fgets(fd);
    txt = fgets(fd);
    
    nnodes = fscanf(fd, '%d', [1, 1]);
    p = zeros(3, nnodes);
    
    txt = fgets(fd);
    
    for i = 1:nnodes
        ind = fread(fd, 1, 'int');
        coords = fread(fd, 3, 'double');
        p(1:3,i) = coords';
    end
    
    %txt = mgetl(fd, 2);
    txt = fgets(fd);
    txt = fgets(fd);
    txt = fgets(fd);
    
    %disp(p);
    
    nelements = fscanf(fd, '%d', [1,1]);
    
    txt = fgets(fd);
    
    tmp = zeros(nelements, 9);
    
    for i = 1:nelements
        %type = fread(fd, 1, 'int');
        %numb = fread(fd, 1, 'int');
        %ntags = fread(fd, 1, 'int');
        L = fread(fd, 3, 'int');
        if(L(1) == 1) % edge element
            %ind = fread(fd, 1, 'int');
            %tags = fread(fd, 2, 'int');
            %nodes = fread(fd, 2, 'int');
            K = fread(fd, 5, 'int');
            
            %tmp(i,:) = [ind type ntags tags' nodes' -1 -1];
            tmp(i,:) = [K(1) L(1) L(3) K(2:3)' K(4:5)' -1 -1];
        elseif(L(1) == 2) % triangle
            %ind = fread(fd, 1, 'int');
            %tags = fread(fd, 2, 'int');
            %nodes = fread(fd, 3, 'int');
            K = fread(fd, 6, 'int');
            
            %tmp(i,:) = [ind type ntags tags' nodes' -1];
            tmp(i,:) = [K(1) L(1) L(3) K(2:3)' K(4:6)' -1];
        elseif(L(1) == 4) % tetrahedron
            %ind = fread(fd, 1, 'int');
            %tags = fread(fd, 2, 'int');
            %nodes = fread(fd, 4, 'int');
            K = fread(fd, 7, 'int');
            
            %tmp(i,:) = [ind type ntags tags' nodes'];
            tmp(i,:) = [K(1) L(1) L(3) K(2:3)' K(4:7)'];
        end
    end
else
    txt = fgets(fd);
    txt = fgets(fd);
    
    nnodes = fscanf(fd, '%d', [1, 1]);
    p = zeros(3, nnodes);
    
    for i = 1:nnodes
        [abcd] = fscanf(fd, '%d %e %e %e', [1,4]);
        p(1:3,i) = [abcd(2:4)]';
    end
    
    %txt = mgetl(fd, 2);
    txt = fgets(fd);
    txt = fgets(fd);
    txt = fgets(fd);
    
    %disp(p);
    
    nelements = fscanf(fd, '%d', [1,1]);
    
    tmp = zeros(nelements, 9);
    
    for i = 1:nelements
        L = fscanf(fd, '%d %d', [1 2]);
        if(L(2) == 1) % edge element
            K = fscanf(fd, '%d %d %d %d %d', [1 5]);
            K = [K -1 -1];
        elseif(L(2) == 2) % triangle
            K = fscanf(fd, '%d %d %d %d %d %d', [1 6]);
            K = [K -1];
        elseif(L(2) == 4) % tetrahedron
            K = fscanf(fd, '%d %d %d %d %d %d %d', [1 7]);
        end
        tmp(i,:) = [L K];
    end
end

fclose(fd);

% Detect mesh space dimension
if(min(p(3,:))==max(p(3,:)))
    dim = 2;
else
    dim = 3;
end

if(dim == 2)
    t = tmp(tmp(:,2)==2,6:8)';
    t = [t; tmp(tmp(:,2)==2,4)'];
else
    t = tmp(tmp(:,2)==4,6:9)';
    t = [t; tmp(tmp(:,2)==4,5)'];
end

if(dim == 2)
    % Odwroc lewe elementy na druga strone
    p = p(1:2,:);
    
    for el = 1:size(t,2)
        n0 = t(1,el);
        n1 = t(2,el);
        n2 = t(3,el);
        
        J = [-p(1,n0)+p(1,n1), -p(1,n0)+p(1,n2);
            -p(2,n0)+p(2,n1), -p(2,n0)+p(2,n2)];
        
        if(det(J) < 0)
            t(2,el) = n2;
            t(3,el) = n1;
        end
    end
    
    % Zbuduj tablice e i dostosuj wlasciwe kierunki
    nedges = nelements-size(t,2);
    
    e = zeros(7,nedges);
    e(1:2,:) = tmp(tmp(:,2)==1,6:7)';
    e(5,:) = tmp(tmp(:,2)==1,4)';
    e(6,:) = ones(1,size(e,2));
    
    % Wygeneruj wlasciwe kierunki brzegu
    % 1. Stworz liste elementow wezla
    elementsOfNode = cell(nnodes,1);
    
    for el = 1:size(t,2)
        elementsOfNode{t(1,el)}(end+1) = el;
        elementsOfNode{t(2,el)}(end+1) = el;
        elementsOfNode{t(3,el)}(end+1) = el;
    end
    % 2. Znajdz wszystkie wezly
    
    for edge = 1:nedges
        n1 = e(1,edge);
        n2 = e(2,edge);
        % Znajdz elementy, ktore maja oba wezly
        el = myintersect(elementsOfNode{n1}, elementsOfNode{n2});
        % Porownaj kierunek
        if(sum(t([1 2],el)'==[n1 n2])~=2 && sum(t([2 3],el)'==[n1 n2])~=2 && sum(t([3 1],el)'==[n1 n2])~=2)
            % Jesli nie zgadza sie kierunek obiegu z zadna ze scianek, to
            % trzeba obrocic zwrot krawedzi edge
            e(2,edge) = n1;
            e(1,edge) = n2;
        end
    end
else % dim == 3
    % Check for negative volumes
    for el = 1:size(t,2)
        % Compute element level gradient matrices
        elementPressureNodes = t(1:4,el);
        x = p(:,elementPressureNodes);
        J = [x(1,2)-x(1,1) x(1,3)-x(1,1) x(1,4)-x(1,1);
            x(2,2)-x(2,1) x(2,3)-x(2,1) x(2,4)-x(2,1);
            x(3,2)-x(3,1) x(3,3)-x(3,1) x(3,4)-x(3,1)];
        
        if(det(J)<0)
            disp('Error: Mesh consists of negative volume elements.')
        end
    end
    
    % Zbuduj tablice e
    nfaces = nelements-size(t,2);
    
    e = zeros(6,nfaces);
    e(1:3,:) = tmp(tmp(:,2)==2,6:8)';
    e(4,:) = tmp(tmp(:,2)==2,5)';
    e(5,:) = ones(1,size(e,2));
    
    % Wygeneruj wlasciwe kierunki brzegu - do zrobienia odtad...
    % 1. Stworz liste elementow wezla
    elementsOfNode = cell(nnodes,1);
    
    for el = 1:size(t,2)
        elementsOfNode{t(1,el)}(end+1) = el;
        elementsOfNode{t(2,el)}(end+1) = el;
        elementsOfNode{t(3,el)}(end+1) = el;
        elementsOfNode{t(4,el)}(end+1) = el;
    end
    % 2. Znajdz wszystkie wezly
    
    for face = 1:nfaces
        n1 = e(1,face);
        n2 = e(2,face);
        n3 = e(3,face);
        % Znajdz element, ktory ma oba wezly
        el = myintersect(myintersect(elementsOfNode{n1}, elementsOfNode{n2}),elementsOfNode{n3});
        
        if size(el,2)>1
            if t(5,el(:,1))==1
                el=el(:,1);
            else
                el=el(:,2);
            end
        end
        % Wyznacz polozenie centroidu elementu i centroidu scianki
        % Policz normalna do scianki i wektor od centroidu scianki do
        % centroidu elementu - jesli iloczyn skalarny jest ujemny,
        % obroc kolejnosc wezlow scianki
        
        v1 = p(:,n2)-p(:,n1);
        v2 = p(:,n3)-p(:,n1);
        nvec = cross(v1,v2);
        
        facecent = sum(p(:,[n1 n2 n3]),2)/3;
        elcent = sum(p(:,t(1:4,el)),2)/4;
        
        v3 = elcent-facecent;
        
        if(nvec'*v3>0)
            % Odwroc kierunek
            e(1:2,face) = flipud(e(1:2,face));
        end
        
        % Kazdej sciance e przyporzadkuj odpowiednie id domeny
        e(5,face) = t(5,el);
        
        %             % Porownaj kierunek
        %             if(sum(t([1 2],el)'==[n1 n2])~=2 && sum(t([2 3],el)'==[n1 n2])~=2 && sum(t([3 1],el)'==[n1 n2])~=2)
        %                 % Jesli nie zgadza sie kierunek obiegu z zadna ze scianek, to
        %                 % trzeba obrocic zwrot krawedzi edge
        %                 e(2,edge) = n1;
        %                 e(1,edge) = n2;
        %             end
    end
    
    % Kazdej sciance e przyporzadkuj odpowiednie id domeny
    
end

% ////////////////////////////////////////////////////////////////
% Obrobka komorek
%     ncells = 0;
%     for i = 1:ncellsgmsh
%         if(tmp(i,2)==2) then
%             ncells = ncells+1;
%         end
%     end
%
%     cells = zeros(ncells,3);
%     domains = zeros(ncells,1);
%
%     k = 1;
%     for i = 1:ncellsgmsh
%         if(tmp(i,2)==2) then
%             cells(k,1:3) = tmp(i,6:8)
%             domains(k) = tmp(i,4);
%             k = k+1;
%         end
%     end
%     // Obroc "lewe" elementy
%     for el = 1:ncells
%         n0 = cells(el,1);
%         n1 = cells(el,2);
%         n2 = cells(el,3);
%
%         J = [-nodes(n0,1)+nodes(n1,1), -nodes(n0,1)+nodes(n2,1);
%             -nodes(n0,2)+nodes(n1,2), -nodes(n0,2)+nodes(n2,2)];
%
%         if(det(J) < 0) then
%             cells(el,2) = n2;
%             cells(el,3) = n1;
%         end
%     end
%
%     //bufor = cells(:,2);
%     //cells(:,2) = cells(:,3);
%     //cells(:,3) = bufor;
%
%     numb = 0;
%     for i = 1:ncellsgmsh
%         if(tmp(i,2)==1) then
%             numb = numb + 1;
%         end
%     end
%     b1 = zeros(numb,2);
%     b2 = zeros(numb,2);
%     k1 = 1;
%     k2 = 1;
%     for i = 1:ncellsgmsh
%         if(tmp(i,4)==b1id) then
%             b1(k1) = tmp(i,6);
%             b1(k1+1) = tmp(i,7);
%             k1 = k1+2;
%         end
%         if(tmp(i,4)==b2id) then
%             b2(k2) = tmp(i,6);
%             b2(k2+1) = tmp(i,7);
%             k2 = k2+2;
%         end
%     end
%
%     nwallb1 = (k1-1)/2;
%     nwallb2 = (k2-1)/2;
%
%     //disp(nwallb1);
%     //disp(nwallb2);
%
%     b1buf = b1;
%     b2buf = b2;
%
%     b1 = unique(b1);
%     b2 = unique(b2);
%
%     if(b1(1)==0) then b1 = b1(2:length(b1)); end
%     if(b2(1)==0) then b2 = b2(2:length(b2)); end
%
%     numb1 = length(b1);
%     numb2 = length(b2);
%
%     domains = domains-min(domains)+1;
%
%     B1 = zeros(nwallb1,2);
%     B2 = zeros(nwallb2,2);
%
%     //disp(B1);
%
%     // Przynaleznosc elementow do wezla
%     powt = zeros(nnodes,1);
%
%     for i = 1:ncells
%         powt(cells(i,1)) = powt(cells(i,1))+1;
%         powt(cells(i,2)) = powt(cells(i,2))+1;
%         powt(cells(i,3)) = powt(cells(i,3))+1;
%     end
%
%     ElsOfNode = zeros(nnodes,max(powt)+1);
%     last = max(powt)+1;
%     clear powt;
%
%     for i = 1:nnodes
%         ElsOfNode(i,last) = 1;
%     end
%
%     for i =1:ncells
%         k = ElsOfNode(cells(i,1),last);
%         ElsOfNode(cells(i,1),k) = i;
%         ElsOfNode(cells(i,1),last) = ElsOfNode(cells(i,1),last) + 1;
%
%         k = ElsOfNode(cells(i,2),last);
%         ElsOfNode(cells(i,2),k) = i;
%         ElsOfNode(cells(i,2),last) = ElsOfNode(cells(i,2),last) + 1;
%
%         k = ElsOfNode(cells(i,3),last);
%         ElsOfNode(cells(i,3),k) = i;
%         ElsOfNode(cells(i,3),last) = ElsOfNode(cells(i,3),last) + 1;
%     end
%
%     // disp(ElsOfNode);
%
%     for i = 1:nwallb1
%         n1 = b1buf(2*i-1);
%         n2 = b1buf(2*i);
%         przedost1 = ElsOfNode(n1, last)-1;
%         przedost2 = ElsOfNode(n2, last)-1;
%         el = intersect(ElsOfNode(n1,1:przedost1),ElsOfNode(n2,1:przedost2));
%         B1(i,1) = el;
%         if(gsort([n1 n2]) == gsort([cells(el,1) cells(el,2)]))
%             wall = 1;
%         end
%         if(gsort([n1 n2]) == gsort([cells(el,2) cells(el,3)]))
%             wall = 2;
%         end
%         if(gsort([n1 n2]) == gsort([cells(el,1) cells(el,3)]))
%             wall = 3;
%         end
%         B1(i,2) = wall;
%         i = i + 1;
%     end
%     for i = 1:nwallb2
%         n1 = b2buf(2*i-1);
%         n2 = b2buf(2*i);
%         przedost1 = ElsOfNode(n1, last)-1;
%         przedost2 = ElsOfNode(n2, last)-1;
%         el = intersect(ElsOfNode(n1,1:przedost1),ElsOfNode(n2,1:przedost2));
%         B2(i,1) = el;
%         if(gsort([n1 n2]) == gsort([cells(el,1) cells(el,2)]))
%             wall = 1;
%         end
%         if(gsort([n1 n2]) == gsort([cells(el,2) cells(el,3)]))
%             wall = 2;
%         end
%         if(gsort([n1 n2]) == gsort([cells(el,1) cells(el,3)]))
%             wall = 3;
%         end
%         B2(i,2) = wall;
%         i = i + 1;
%     end
%
%     // Dodaj dodatkowe wezly w polowie krawedzi
%     Edges = list(1:nnodes);
%     AddNodes = list(1:nnodes);
%
%     for i = 1:nnodes
%         Edges(i) = list(i);
%         AddNodes(i) = list(i);
%     end
%
%     cells = [cells zeros(ncells,3)];
%
%     el = nnodes+1;
%
%     //disp(cells)
%
%     for i = 1:ncells
%         for j = 1:3
%             n1 = cells(i,j);
%             n2 = cells(i,modulo(j,3)+1);
%             nloc = gsort([n1 n2]);
%             n1 = nloc(2);
%             n2 = nloc(1);
%             if(ZnajdzNaLiscie(Edges(n1),n2)) then
%                 ind = ZnajdzNaLiscie(Edges(n1),n2)
%                 cells(i,3+j) = AddNodes(n1)(ind);
%             else
%                 //disp([n1 n2 el]);
%                 Edges(n1)($+1)=n2;
%                 AddNodes(n1)($+1) = el;
%                 cells(i,3+j) = el;
%                 el = el+1;
%             end
%
%         end
%     end
%     el = el-1;
%
%     przyrost = el-nnodes;
%     nnodesfull = el;
%
%     // Zwieksz tablice wezlow i uzupelnij wspolrzednymi
%
%     nodes = [nodes; zeros(przyrost,2)];
%
%     for i = 1:ncells
%         for j = 1:3
%             n1 = cells(i,j);
%             n2 = cells(i,modulo(j,3)+1);
%
%             nodes(cells(i,j+3),:) = 0.5*(nodes(n1,:)+nodes(n2,:));
%         end
%     end
%
%     //disp(Edges)
%     //disp(AddNodes)

end

function c = myintersect(a, b)
    c = a(logical(sum(repmat(a, length(b), 1) == repmat(b', 1, length(a)), 1)));
end