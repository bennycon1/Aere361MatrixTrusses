clear,clc
table = xlsread('members.xlsx');
forcetable1 = xlsread('forces.xlsx');
forcetable = forcetable1(:,1);
%forcetable = [0;0;0;0;100000*cosd(60);-100000*sind(60);0;0];
members = size(table,1); %number of members
size = max(max(table(:,2)),max(table(:,3)));
%for each member, fills out the remainder of the table
for i = 1:members
    table(i,8) = cosd(table(i,7));
    table(i,9) = sind(table(i,7));
end
disp('Members,startpin, endpin, E, A, L, theta, costheta, sintheta');
disp(table);
KLBase = [1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
%calculates the local stiffness matrix and trasnformation matrix for each
%member
for i = 1:members
    KL(:,:,i) = (table(i,4)*table(i,5)/table(i,6))*KLBase;
    TL(:,:,i) = [table(i,8) table(i,9) 0 0; -1*table(i,9) table(i,8) 0 0; 0 0 table(i,8) table(i,9); 0 0 -1*table(i,9) table(i,8)];
end
%calculates the stiffness matrix for each member in global coordinate
for i = 1:members
    KGL(:,:,i) = TL(:,:,i).'*(KL(:,:,i))*TL(:,:,i);
end
%Assembly... this is where the problems arise,  we know pins now how do we
%format matriceis to realize that
KGG = zeros(size*2); %How do i know the size of KGG? is it number of members*2, number of pins *2?
%go through and standarize the arrays based on the sizes of thigns?? sure
%lets go with that wording
KGLupdated = zeros(size*2);

%Create an array of the values for each table to use for row and columns,
%hard to explain in text
for i=1:members
    vecvalues(:,:,i) = [[table(i,2)*2-1,table(i,2)*2],[table(i,3)*2-1,table(i,3)*2]];
end

%going through orginal KGL and updating the position for the correct size
%of KGG
%
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
%I need to know the given diflections, because we know that certain joints
%are 0 we can neglect this
%invkgg = inv(KGG);
%DLETE THIS LATER AND LOAD FROM FILE!
%Enter 1 if member can move in this direction (base conditions)
%knowndiflections = [0,0,0,0,1,1,0,0];
knowndiflections = forcetable1(:,2);
%go through and basing it off the diflections, calculate the matrix we
%need, basically 0 out every other column?
%loop through the array and if the value is 0 delete that value from the
%array?
KGGTest = KGG;
for i=1:length(knowndiflections)
    if(knowndiflections(i) == 0)
        KGGTest(:,i) = 0;
        KGGTest(i,:) = 0;
    end
end
%can i take the inverse of this, NO BECAUSE DET IS STILL 0?
notsurewhattocallthis = nonzeros(KGGTest);
%we know that its a 2x2 because there are 2 nonzero values in known
%difelctions, will this always be the case though?, I think so as long
%there is not a roller, I will have to go back and update this laster

%So here I am taking a vector of all the non zero values (could cause an
%issue if the value is actually 0 in the array)
BoundaryKGMatrix = vec2mat(notsurewhattocallthis,sum(knowndiflections));

%now that we have the boundary matrix we hsould be able to solve for the
%unknown forces
%We need to find the forces that are being apllied on in the matrix, this
%should just be the forces that correspond to the 1 values of the
%knowndiflections i beleive.
counter = 1;
for i = 1:length(forcetable)
    if(knowndiflections(i)~= 0)
        BoundaryForceMatrix(counter) = forcetable(i);
        counter = counter +1;
    end
    
end
BoundaryForceMatrix = transpose(BoundaryForceMatrix);
Reactions = BoundaryForceMatrix;
DisplacemntValues = BoundaryKGMatrix\BoundaryForceMatrix;

%Now we need to solve for the reaction forces in teh strucutre, to do this
%we just multipl ythe KG by the new dispacement vales we just obtained,
%first i need to put the dispacment values I got  from before and put them
%in the correct spot in the matrix
ReactionForcesDispacementMatrix = zeros(length(knowndiflections),1);
counter = 1;
for i=1:length(knowndiflections)
    if(knowndiflections(i) ~= 0)
        ReactionForcesDispacementMatrix(i) = DisplacemntValues(counter);
        counter = counter +1;
    end
end
disp('Deflection of joints');
disp(ReactionForcesDispacementMatrix)

%now we should solve for the forces
ReactionForces = KGG*ReactionForcesDispacementMatrix;
%go through the matrix and round the values, due to some computer rounding
%error we can be left with a value so small that we can basically asumme it
%is 0
ReactionForces=round(ReactionForces,8);
disp('Reaction Forces');
disp(ReactionForces);
%Now we have to go through and find the internal forces within each
%structureal member I guess now i need ot ifgure out what kind of members i
%need to add?

knowndiflections = knowndiflections';
counter = 1;
for i = 1:members
    for j =2:3
        temp(counter) = table(i,j)*2-1;
        counter = counter +1;
        temp(counter) = table(i,j)*2;
        counter = counter + 1;
    end
    DispacmentPerMember(:,:,i) = temp';
    
    counter = 1;
end


for j = 1:members
    for i = 1:length(temp)
        DispacmentPerMemberUpdated(i,1,j) = ReactionForcesDispacementMatrix(DispacmentPerMember(i,1,j));
    end
end



for i = 1:members
    InternalForces(:,:,i) = KGL(:,:,i) * DispacmentPerMemberUpdated(:,:,i);
end
disp('Internal Forces');
disp(InternalForces);

