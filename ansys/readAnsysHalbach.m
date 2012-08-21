clear all;clc;
fid=fopen('NLIST.lis');

A1 = [];
A2 = [];
A3 = [];

A_data=textscan(fid,'%f %f %f %f','HeaderLines',5);
A1 = [A1;A_data{1}];
A2 = [A2;A_data{2}];
A3 = [A3;A_data{3}];

while ~feof(fid)
    A_data=textscan(fid,'%f %f %f %f','HeaderLines',1);
    A1 = [A1;A_data{1}];
    A2 = [A2;A_data{2}];
    A3 = [A3;A_data{3}];
end

fid=fopen('PRVECT.lis');

B1 = [];
B2 = [];
B3 = [];
B4 = [];

B_data=textscan(fid,'%f %f %f %f %f','HeaderLines',9);
B1 = [B1;B_data{1}];
B2 = [B2;B_data{2}];
B3 = [B3;B_data{3}];
B4 = [B4;B_data{4}];

while ~feof(fid)
    B_data=textscan(fid,'%f %f %f %f %f','HeaderLines',1);
    B1 = [B1;B_data{1}];
    B2 = [B2;B_data{2}];
    B3 = [B3;B_data{3}];
    B4 = [B4;B_data{4}];
end

for i=1:14925
    for j=1:44545
        if A1(j)==B1(i)
           Z(j,1)=[1];
        end
        j=j+1;
    end
       i=i+1; 
end

C1=Z.*A1;
% C2=Z.*A2;
% C3=Z.*A3;

condition=C1(:,1)==0;
C1(condition,:)=[];
A2(condition,:)=[];
A3(condition,:)=[];

caxis([min(B3) max(B3)]);
scatter(A2,A3,[],B3,'filled');

figure(2);
x1=linspace(min(A2),max(A2),1000);
y1=linspace(min(A3),max(A3),1000);
[X1 Y1]=meshgrid(x1,y1);
Z1=griddata(A2,A3,B3,X1,Y1);
contour(X1,Y1,Z1);
