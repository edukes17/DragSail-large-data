function [ Centroids, Normals, Areas, Faces, Tangents] = MeshData( STL_Filename, View_Geometry, Scale_Geometry, Spacecraft, View_Normals, View_Tangents, View_Centroid_Labels,outdir)
%% MeshData.m function

%Purpose: Read the input geometry file in STL format to obtain the
%centroids, normals, and areas of surfaces. 
 
%Created:  Sebastian Tamrazian 10/27/2018
%Modified: Sebastian Tamrazian 10/27/2018
%Modified: Arly Black 9/26/2019 - added tangent vectors
 
%Inputs:
%     STL_Filename:   filename of .stl geometry. .stl can be in ASCII or binary format

%Outputs:
%     Centroids: Centroid location of each geometry face
%     Normals: Normal vector of each geometry face
%     Areas: Area of each geometry face
%     Faces: list of faces of geometry
%     Tangents: Tangent vector of each geometry face

%Supporting Functions: 
%     stlread.m
%     stlGetFormat.m
%     stlReadAscii.m
%     stlReadBinary.m

%% Read STL file
[Vertices,Faces,Normals] = stlread(STL_Filename);  %import verticies, faces,and normal vectors
Vertices = Vertices/Scale_Geometry;                %scale the geometry so that all measurements are in meters
          
%% Translate Geometry so that the cm lies on [0,0,0]
Vertices = bsxfun(@plus,Vertices, Spacecraft.cm);

%% Determine area of each face
Areas  = zeros(length(Vertices)/3,1);       %initialize area vector
nn = 1;                                     %counter

for ii = 1:3:(length(Vertices)-2)           %calculate surface area of each triangular face
    
    %determine side lengths using distance formula
    s1 = sqrt((Vertices(ii  ,1)-Vertices(ii+1,1))^2+(Vertices(ii  ,2)-Vertices(ii+1,2))^2+(Vertices(ii  ,3)-Vertices(ii+1,3))^2);
    s2 = sqrt((Vertices(ii+1,1)-Vertices(ii+2,1))^2+(Vertices(ii+1,2)-Vertices(ii+2,2))^2+(Vertices(ii+1,3)-Vertices(ii+2,3))^2);
    s3 = sqrt((Vertices(ii  ,1)-Vertices(ii+2,1))^2+(Vertices(ii  ,2)-Vertices(ii+2,2))^2+(Vertices(ii  ,3)-Vertices(ii+2,3))^2);
    
    %apply Heron's formula using side lengths
    Areas(nn) = sqrt((s1+s2+s3)/2*((s1+s2+s3)/2-s1)*((s1+s2+s3)/2-s2)*((s1+s2+s3)/2-s3));
    nn = nn+1;
end

%% Determine centroid of each face
Centroids  = zeros(length(Vertices)/3,3);   %initialize centroid vector
nn = 1;                         

for ii = 1:3:(length(Vertices)-2)           %calculate centroid each triangular face
    Centroids(nn,1) = (Vertices(ii,1)+Vertices(ii+1,1)+Vertices(ii+2,1))/3;
    Centroids(nn,2) = (Vertices(ii,2)+Vertices(ii+1,2)+Vertices(ii+2,2))/3;
    Centroids(nn,3) = (Vertices(ii,3)+Vertices(ii+1,3)+Vertices(ii+2,3))/3;
    nn = nn+1;
end

%% Determine tangent vector of each face
Tangents = cross(Centroids,Normals);

%% Plot Geometry
%If the user has specified to view the geometry in the main function, this
%code will plot the stl file

if View_Geometry
    
    if ~exist(outdir,'dir') %make output directory if it doesn't exist already
    mkdir(outdir)
    end
    
    Display_Structure = struct('faces',Faces,'vertices',Vertices);
%     f(1) = figure('Visible', 'off');
    figure(1)
    hold on
    patch(Display_Structure,'FaceColor',       [0.9 0.9 0.9], ...
             'EdgeColor',       'b',        ...
             'FaceColor',    'none', 'Linewidth',1);
    if View_Normals
        quiver3(Centroids(:,1),Centroids(:,2),Centroids(:,3),Normals(:,1),Normals(:,2),Normals(:,3))
    end
    if View_Tangents
        quiver3(Centroids(:,1),Centroids(:,2),Centroids(:,3),Tangents(:,1),Tangents(:,2),Tangents(:,3))
    end
    if View_Centroid_Labels
        for i = 1:length(Centroids)
            str = sprintf('%d',i);
            text(Centroids(i,1),Centroids(i,2),Centroids(i,3),str,'FontSize',12)
        end
    end
    xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
    axis tight;
    axis equal;
    hold off
    view([-135 35]);
    
    if ispc
    savefig([pwd '\' outdir '\' 'geometry'])
elseif ismac
    savefig([pwd '/' outdir '/' 'geometry'])
end
    
end
end

%% STL Reader Functions
% BEFORE USING ELSEWHERE, READ LICENSE FILE IN GEOMETRY SUBFOLDER
% These functions determine the type of STL being input (binary or ASCII),
% and return verticies, faces, and normal vectors for the data points
% contained in the file specified

function [v, f, n, name] = stlread(fileName)
%STLREAD reads any STL file not depending on its format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

format = stlGetFormat(fileName);
if strcmp(format,'ascii')
  [v,f,n,name] = stlReadAscii(fileName);
elseif strcmp(format,'binary')
  [v,f,n,name] = stlReadBinary(fileName);
end
end

function format = stlGetFormat(fileName)
%STLGETFORMAT identifies the format of the STL file and returns 'binary' or
%'ascii'

fid = fopen(fileName);
% Check the file size first, since binary files MUST have a size of 84+(50*n)
fseek(fid,0,1);         % Go to the end of the file
fidSIZE = ftell(fid);   % Check the size of the file
if rem(fidSIZE-84,50) > 0
    format = 'ascii';
else
    % Files with a size of 84+(50*n), might be either ascii or binary...
    % Read first 80 characters of the file.
    % For an ASCII file, the data should begin immediately (give or take a few
    % blank lines or spaces) and the first word must be 'solid'.
    % For a binary file, the first 80 characters contains the header.
    % It is bad practice to begin the header of a binary file with the word
    % 'solid', so it can be used to identify whether the file is ASCII or
    % binary.
    fseek(fid,0,-1); % go to the beginning of the file
    header = strtrim(char(fread(fid,80,'uchar')')); % trim leading and trailing spaces
    isSolid = strcmp(header(1:min(5,length(header))),'solid'); % take first 5 char
    fseek(fid,-80,1); % go to the end of the file minus 80 characters
    tail = char(fread(fid,80,'uchar')');
    isEndSolid = findstr(tail,'endsolid');
    
    % Double check by reading the last 80 characters of the file.
    % For an ASCII file, the data should end (give or take a few
    % blank lines or spaces) with 'endsolid <object_name>'.
    % If the last 80 characters contains the word 'endsolid' then this
    % confirms that the file is indeed ASCII.
    if isSolid & isEndSolid
        format = 'ascii';
    else
        format = 'binary';
    end
end
fclose(fid);
end

function [v, f, n, name] = stlReadAscii(fileName)
%STLREADASCII reads a STL file written in ASCII format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

%======================
% STL ascii file format
%======================
% ASCII STL files have the following structure.  Technically each facet
% could be any 2D shape, but in practice only triangular facets tend to be
% used.  The present code ONLY works for meshes composed of triangular
% facets.
%
% solid object_name
% facet normal x y z
%   outer loop
%     vertex x y z
%     vertex x y z
%     vertex x y z
%   endloop
% endfacet
%
% <Repeat for all facets...>
%
% endsolid object_name

fid = fopen(fileName);
cellcontent = textscan(fid,'%s','delimiter','\n'); % read all the file and put content in cells
content = cellcontent{:}(logical(~strcmp(cellcontent{:},''))); % remove all blank lines
fclose(fid);

% read the STL name
line1 = char(content(1));
if (size(line1,2) >= 7)
    name = line1(7:end);
else
    name = 'Unnamed Object';
end

% read the vector normals
normals = char(content(logical(strncmp(content,'facet normal',12))));
n = str2num(normals(:,13:end));

% read the vertex coordinates (vertices)
vertices = char(content(logical(strncmp(content,'vertex',6))));
v = str2num(vertices(:,7:end));
nvert = size(vertices,1); % number of vertices
nfaces = sum(strcmp(content,'endfacet')); % number of faces
if (nvert == 3*nfaces)
    f = reshape(1:nvert,[3 nfaces])'; % create faces
end
end

function [v, f, n, name] = stlReadBinary(fileName)
%STLREADBINARY reads a STL file written in BINARY format
%V are the vertices
%F are the faces
%N are the normals
%NAME is the name of the STL object (NOT the name of the STL file)

%=======================
% STL binary file format
%=======================
% Binary STL files have an 84 byte header followed by 50-byte records, each
% describing a single facet of the mesh.  Technically each facet could be
% any 2D shape, but that would screw up the 50-byte-per-facet structure, so
% in practice only triangular facets are used.  The present code ONLY works
% for meshes composed of triangular facets.
%
% HEADER:
% 80 bytes:  Header text
% 4 bytes:   (int) The number of facets in the STL mesh
%
% DATA:
% 4 bytes:  (float) normal x
% 4 bytes:  (float) normal y
% 4 bytes:  (float) normal z
% 4 bytes:  (float) vertex1 x
% 4 bytes:  (float) vertex1 y
% 4 bytes:  (float) vertex1 z
% 4 bytes:  (float) vertex2 x
% 4 bytes:  (float) vertex2 y
% 4 bytes:  (float) vertex2 z
% 4 bytes:  (float) vertex3 x
% 4 bytes:  (float) vertex3 y
% 4 bytes:  (float) vertex3 z
% 2 bytes:  Padding to make the data for each facet 50-bytes in length
%   ...and repeat for next facet... 

fid = fopen(fileName);
header = fread(fid,80,'int8'); % reading header's 80 bytes
name = deblank(native2unicode(header,'ascii')');
if isempty(name)
    name = 'Unnamed Object'; % no object name in binary files!
end
nfaces = fread(fid,1,'int32');  % reading the number of facets in the stl file (next 4 byters)
nvert = 3*nfaces; % number of vertices
% reserve memory for vectors (increase the processing speed)
n = zeros(nfaces,3);
v = zeros(nvert,3);
f = zeros(nfaces,3);
for i = 1 : nfaces % read the data for each facet
    tmp = fread(fid,3*4,'float'); % read coordinates
    n(i,:) = tmp(1:3); % x,y,z components of the facet's normal vector
    v(3*i-2,:) = tmp(4:6); % x,y,z coordinates of vertex 1
    v(3*i-1,:) = tmp(7:9); % x,y,z coordinates of vertex 2
    v(3*i,:) = tmp(10:12); % x,y,z coordinates of vertex 3
    f(i,:) = [3*i-2 3*i-1 3*i]; % face
    fread(fid,1,'int16'); % Move to the start of the next facet (2 bytes of padding)
end
fclose(fid);
end
