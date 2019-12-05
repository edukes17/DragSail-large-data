function [a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z] = PropertiesFile(FN)

% Enter filename
fileName = sprintf(FN); 

% Open file
fid = fopen(fileName);
% hdrRows = 4;
% hdrData = textscan(fid,'%s',hdrRows, 'Delimiter','\n');     % define header rows
matData = textscan(fid,'%s%s%s%s%s%s%s', 'Delimiter',{',','=','%'}, 'CollectOutput',true);
idx1 = find(contains(matData{1}, 'STL Filename'),1);     % Search for text containing "STL Filename" and identify the row number
idx2 = find(contains(matData{1}, 'temperature'),1);
idx3 = find(contains(matData{1}, 'reflectance_coef'),1);
idx4 = find(contains(matData{1}, 'specular_coef'),1);
idx5 = find(contains(matData{1}, 'front_emiss_coef'),1);
idx6 = find(contains(matData{1}, 'back_emiss_coef'),1);
idx7 = find(contains(matData{1}, 'front_nonLamb_coef'),1);
idx8 = find(contains(matData{1}, 'back_nonLamb_coef'),1);
idx9 = find(contains(matData{1}, 'mass'),1);
idx10 = find(contains(matData{1}, 'inertia'),1);
idx11 = find(contains(matData{1}, 'cm'),1);
idx12 = find(contains(matData{1}, 'Semimajor'),1);
idx13 = find(contains(matData{1}, 'Eccentricity'),1);
idx14 = find(contains(matData{1}, 'Inclination'),1);
idx15 = find(contains(matData{1}, 'Argument'),1);
idx16 = find(contains(matData{1}, 'Ascension'),1);
idx17 = find(contains(matData{1}, 'True Anomaly'),1);
idx18 = find(contains(matData{1}, 'Output Directory'),1);
idx19 = find(contains(matData{1}, 'Roll angle'),1);
idx20 = find(contains(matData{1}, 'Pitch angle'),1);
idx21 = find(contains(matData{1}, 'Yaw angle'),1);
idx22 = find(contains(matData{1}, 'Roll rate'),1);
idx23 = find(contains(matData{1}, 'Pitch rate'),1);
idx24 = find(contains(matData{1}, 'Yaw rate'),1);
idx25 = find(contains(matData{1}, 'Epoch'),1);
idx26 = find(contains(matData{1}, 'Propagation'),1);

fclose(fid);

info = matData{1,1};

% Identify which row and column the relevant value is in. Convert from
% string to double
a = (info(idx1,2));
b = str2double(info(idx2,2));
c = str2double(info(idx3,2));
d = str2double(info(idx4,2));
e = str2double(info(idx5,2));
f = str2double(info(idx6,2));
g = str2double(info(idx7,2));
h = str2double(info(idx8,2));
i = str2double(info(idx9,2));
j = [str2double(info(idx10,2)),str2double(info(idx10,3)),str2double(info(idx10,4));str2double(info(idx10+1,2)),str2double(info(idx10+1,3)),str2double(info(idx10+1,4));str2double(info(idx10+2,2)),str2double(info(idx10+2,3)),str2double(info(idx10+2,4))];
k = [str2double(info(idx11,2)),str2double(info(idx11,3)),str2double(info(idx11,4))];
l = str2double(info(idx12,2));
m = str2double(info(idx13,2));
n = str2double(info(idx14,2));
o = str2double(info(idx15,2));
p = str2double(info(idx16,2));
q = str2double(info(idx17,2));
r = (info(idx18,2));
s = str2double(info(idx19,2));
t = str2double(info(idx20,2));
u = str2double(info(idx21,2));
v = str2double(info(idx22,2));
w = str2double(info(idx23,2));
x = str2double(info(idx24,2));
y =(info(idx25,2));
z = str2double(info(idx26,2));
end