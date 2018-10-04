clear,clc
table = xlsread('members.xlsx');
forcetable = xlsread('forces.xlsx');
forcetable = [0;-10;0;0;0;0];
members = size(table,1); %number of members
size = max((table(:,3)));
%for each member, fills out the remainder of the table
for i = 1:members
    table(i,8) = cosd(table(i,7));
    table(i,9) = sind(table(i,7));
end
disp('Members,startpin, endpin, E, A, L, theta, costheta, sintheta');
disp(table);
KLBase = [1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
disp(KLBase);
%calculates the local stiffness matrix and trasnformation matrix for each
%member
for i = 1:members
   KL(:,:,i) = (table(i,4)*table(i,5)/table(i,6))*KLBase;
   TL(:,:,i) = [table(i,8) table(i,9) 0 0; -1*table(i,9) table(i,8) 0 0; 0 0 table(i,8) table(i,9); 0 0 -1*table(i,9) table(i,8)];
end
%calculates the stiffness matrix for each member in global coordinate
for i = 1:members
  KGL(:,:,i) = transpose(TL(:,:,i))*(KL(:,:,i))*TL(:,:,i);
end
%Assembly... this is where the problems arise,  we know pins now how do we
%format matriceis to realize that
KGG = zeros(size*2); %How do i know the size of KGG? is it number of members*2, number of pins *2?
%go through and standarize the arrays based on the sizes of thigns?? sure
%lets go with that wording
KGLupdated = zeros(size*2);

%Create an array of the values for each table to use for row and columns,
%hard to explain in text
for i=1:size
    vecvalues(:,:,i) = [[table(i,2)*2-1,table(i,2)*2],[table(i,3)*2-1,table(i,3)*2]];
end

%going through orginal KGL and updating the position for the correct size
%of KGG
for k=1:members
    for i=1:4 %rows
        for j=1:4 %cols
         KGLupdated(vecvalues(1,i,k),vecvalues(1,j,k),k) = KGL(i,j,k); 
         %fprintf('KGLupdated array location %d %d \n',vecvalues(1,i,k),vecvalues(1,j,k))
        end
    end
end
       
%add all of the matricies together?
for i=1:members
    KGG = KGG + KGLupdated(:,:,i);
end

%solving for diflection in the beams:
%invkgg = inv(KGG);
%Diflection = forcetable * inv(KGG);
%x = mldivide(KGG,forcetable);

