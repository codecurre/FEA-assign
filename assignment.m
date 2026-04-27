%%Omar Mohamed Sayed | 2200198
%%MDP411s | intro to finite elements
%%assignment submission

clc; clear all; close all; %restarting everything

%%Mesh Generation:
%%---------------------------------------------------------------

coord = csvread('nodes.csv'); %node matrix which has the coordinates of the node in space with each row index representing one node
elem = csvread('connectivity.csv'); %elements --which is a matrix that has 2 columns, each row has the two nodes that make up the element
nnod = length(coord); %this is the number of nodes that we have
nel = length(elem); %this is the number of elements I had
ndof = 3*nnod; %the number of degress of freedom in my system

%%parameter definition
%%---------------------------------------------------------------

FCR_pos = 200; % Newton
FDR = 5000; % Newton
E = 210e9; % Pascal - assumed to be steel
R = 2e-2; %meter
A = pi * (R)^2; %meter^2
%%below is the
dX = coord(elem(:,2),1) - coord(elem(:,1),1); %this is the difference between the x of the second node and the x of the first, matrix (nel,1)
dY = coord(elem(:,2),2) - coord(elem(:,1),2); %this is the difference between the y of the second node and the y of the first, matrix (nel,1)
dZ = coord(elem(:,2),3) - coord(elem(:,1),3); %this is the difference between the z of the second node and the z of the first, matrix (nel,1)
L = sqrt(dX.^2 + dY.^2 + dZ.^2); %this is the length of each element put into an (nel,1) matrix
LMN = [dX./L , dY./L , dZ./L]; %matrix that has l,m,n for each element in the truss ,corresponding to each element, matrix(nel,3)

%%Boundary Conditions:
%%------------------------------------------------------------------

%from CSV file it is clear that the first 4 nodes are grounded
fixdof = [1,2,3, 4,5,6, 7,8,9, 10,11,12]; %as the first 12 degrees of freedom(4 nodes) are suppports
freedof = setdiff(1:ndof,fixdof); %this function creates a matrix which is (1,ndof-fixdof), where the elements in the row are numbers starting from 1 to ndof ,excluding the elements fixdof
Kg = zeros(ndof,ndof); %intialization of global stiffness matrix
Fg = zeros(ndof,1); %initialization of global force matrix
Ug = zeros(ndof,1); %initialization of global displacement matrix
sigma_g = zeros(nel,1); %initialization of global stress matrix (which is sized (nel,1))

%%-------------------------------------------------------------------
%%stiffness matrix function
function k = truss_stiff(E,A,L,l,m,n)

 k = ((E*A)/L)* [ l^2,    l*m,    l*n,   -l^2,   -l*m,   -l*n;
                 l*m,    m^2,    m*n,   -l*m,   -m^2,   -m*n;
                 l*n,    m*n,    n^2,   -l*n,   -m*n,   -n^2;
                -l^2,   -l*m,   -l*n,    l^2,    l*m,    l*n;
                -l*m,   -m^2,   -m*n,    l*m,    m^2,    m*n;
                -l*n,   -m*n,   -n^2,    l*n,    m*n,    n^2 ];
endfunction



%%FORCE LOCATION
%%------------------------------------------------------------------

z_max = max(coord(:,3));
top_nodes = find(abs(coord(:,3) - z_max) < 1e-6);

while true
  choice = input("case number: ");
  switch choice
     case 1
       for i = top_nodes'
         X = coord(i,1);
         Y = coord(i,2);
         if X == 0 | Y == 0
           continue
         else
           Fg(3*i-1,1) = -FDR;
           if X > 0
             Fg(3*i,1) = FCR_pos;
           else
             Fg(3*i,1) = -FCR_pos;
           endif
         endif
       end
       break

     case 2
       for i = top_nodes'
         X = coord(i,1);
         Y = coord(i,2);
         if X == 0 | Y == 0
           continue
         else
           Fg(3*i-1,1) = -FDR;
           if X > 0
           Fg(3*i,1) = -FCR_pos;
           else
           Fg(3*i,1) = FCR_pos;
           endif
         endif
       end
       break

     case 3
       for i = top_nodes'
         X = coord(i,1);
         Y = coord(i,2);
         if X == 0 | Y == 0
           continue
         else
           Fg(3*i-1,1) = -FDR;
           if X > 0
           Fg(3*i-2,1) = FCR_pos;
           else
           Fg(3*i-2,1) = -FCR_pos;
           endif
         endif
       end
       break

     case 4
       for i = top_nodes'
         X = coord(i,1);
         Y = coord(i,2);
         if X == 0 | Y == 0
           continue
         else
           Fg(3*i-1,1) = -FDR;
           if X > 0
           Fg(3*i-2,1) = -FCR_pos;
           else
           Fg(3*i-2,1) = FCR_pos;
           endif
         endif
       end
       break

     otherwise
      printf("invalid entry please try again");
      continue;

  end
end


%%------------------------------------------------------------------------------------
%%assembly and calculations
%%------------------------------------------------------------------------------------




%%------------------------------------------------------------------------------------
%%filling global stiffness matrix

for i = 1 : nel

  dofe = [3*elem(i,1)-2 3*elem(i,1)-1 3*elem(i,1) 3*elem(i,2)-2 3*elem(i,2)-1 3*elem(i,2)];


  l = LMN(i,1);
  m = LMN(i,2);
  n = LMN(i,3);

  ke = truss_stiff(E,A,L(i),l,m,n);

  Kg(dofe,dofe) = Kg(dofe,dofe) + ke;


end

%%-------------------------------------------------------------------------------------
%%calculating displacements

Ug(freedof,1) = Kg(freedof,freedof)\Fg(freedof,1);

%%-------------------------------------------------------------------------------------
%%Calculating reaction forces at the supports

Rg = Kg * Ug - Fg;


F_reaction = Rg(1:length(fixdof),1);

%%-------------------------------------------------------------------------------------
%%stress calculation
for i = 1 : nel

  l = LMN(i,1);
  m = LMN(i,2);
  n = LMN(i,3);
  node1 = Ug(3*elem(i,1)-2,1);
  node2 = Ug(3*elem(i,1)-1,1);
  node3 = Ug(3*elem(i,1),1);
  node4 = Ug(3*elem(i,2)-2,1);
  node5 = Ug(3*elem(i,2)-1,1);
  node6 = Ug(3*elem(i,2),1);

    sigma_g(i,1) = (E/L(i)) * [-l -m -n l m n] * [node1;node2;node3;node4;node5;node6];
end

[max_tensile, elem_tensile] = max(sigma_g);
[max_compressive, elem_compressive] = min(sigma_g);

fprintf('------------------------------------------------------\n')
fprintf('stress results\n');
fprintf('max tensile stress:     %.4f Pa   element %d\n', max_tensile, elem_tensile);
fprintf('max compressive stress: %.4f Pa   element %d\n', max_compressive, elem_compressive);

%%-----------------------------------------------------------------------------------
%%visualization
scale = 10;

def_coord = zeros(size(coord));
def_coord(:,1) = coord(:,1) + scale * Ug(1:3:end);
def_coord(:,2) = coord(:,2) + scale * Ug(2:3:end);
def_coord(:,3) = coord(:,3) + scale * Ug(3:3:end);

% figure 1: original only
figure;
hold on;
for i = 1:nel
    n1 = elem(i,1);
    n2 = elem(i,2);
    plot3([coord(n1,1), coord(n2,1)], ...
          [coord(n1,2), coord(n2,2)], ...
          [coord(n1,3), coord(n2,3)], 'b');
end

grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title('Original Truss');
view(3);
hold off;

% figure 2: deformed only
figure;
hold on;
for i = 1:nel
    n1 = elem(i,1);
    n2 = elem(i,2);
    plot3([def_coord(n1,1), def_coord(n2,1)], ...
          [def_coord(n1,2), def_coord(n2,2)], ...
          [def_coord(n1,3), def_coord(n2,3)], 'r');
end
grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(['Deformed Truss (x' num2str(scale) ' scale)']);
view(3);
hold off;

% figre 3: overlay
figure;
hold on;
for i = 1:nel
    n1 = elem(i,1);
    n2 = elem(i,2);
    plot3([coord(n1,1), coord(n2,1)], ...
          [coord(n1,2), coord(n2,2)], ...
          [coord(n1,3), coord(n2,3)], 'b');
    plot3([def_coord(n1,1), def_coord(n2,1)], ...
          [def_coord(n1,2), def_coord(n2,2)], ...
          [def_coord(n1,3), def_coord(n2,3)], 'r');
end
grid on; axis equal;
xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
title(['Original (blue) vs Deformed (red) — x' num2str(scale) ' scale']);
legend('Original', 'Deformed');
view(3);
hold off;


%some code was inspired by truss13bar :)
