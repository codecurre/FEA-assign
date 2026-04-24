%%Omar Mohamed Sayed | 2200198
%%MDP411s | intro to finite elements
%%assignment submission

clc; clear all; close all;

%%Mesh Generation:
%%---------------------------------------------------------------

coord = csvread('nodes.csv'); %nodes
elem = csvread('connectivity.csv'); %elements
nnod = length(coord);
nel = length(elem);
ndof = 3*nnod;

%%parameter definition
%%---------------------------------------------------------------

FCR_pos = 200; %N
FDR = 5000; %N
E = 210e9; %Pa
R = 20e-3; %m
A = pi * (R)^2; %m^2
dX = coord(elem(:,2),1) - coord(elem(:,1),1); %m
dY = coord(elem(:,2),2) - coord(elem(:,1),2); %m
dZ = coord(elem(:,2),3) - coord(elem(:,1),3); %m
L = sqrt(dX.^2 + dY.^2 + dZ.^2); %m
LMN = [dX./L , dY./L , dZ./L]; %matrix that has l,m,n for each element in the truss

%%Boundary Conditions:
%%------------------------------------------------------------------

%from CSV file it is clear that the first 4 nodes are grounded
fixdof = [1,2,3, 4,5,6, 7,8,9, 10,11,12];
freedof = setdiff(1:ndof,fixdof);
Kg = zeros(ndof,ndof); Fg = zeros(ndof,1); Ug = zeros(ndof,1); sigma_g = zeros(nel,1);

%%------------------------------------------------------------------
%%FORCE LOCATION


while true
  choice = input("enter the case number from 1 to 4(different cases of centrifugal force direction): ");
  switch choice
     case 1
       for i = 5 : nnod
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
       for i = 5 : nnod
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
       for i = 5 : nnod
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
       for i = 5 : nnod
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
%%assembly

%TYPE HERE

%%-------------------------------------------------------------------
%%stiffness matrix
function k = truss_stiff(E,A,L,l,m,n)

 k = ((E*A)/L)* [ l^2,    l*m,    l*n,   -l^2,   -l*m,   -l*n;
                 l*m,    m^2,    m*n,   -l*m,   -m^2,   -m*n;
                 l*n,    m*n,    n^2,   -l*n,   -m*n,   -n^2;
                -l^2,   -l*m,   -l*n,    l^2,    l*m,    l*n;
                -l*m,   -m^2,   -m*n,    l*m,    m^2,    m*n;
                -l*n,   -m*n,   -n^2,    l*n,    m*n,    n^2 ];
endfunction
%calculating stress for each section
for i = 1 : nel

  dofe = [3*elem(i,1)-2 3*elem(i,1)-1 3*elem(i,1) 3*elem(i,2)-2 3*elem(i,2)-1 3*elem(i,2)];


  l = LMN(i,1);
  m = LMN(i,2);
  n = LMN(i,3);

  ke = truss_stiff(E,A,L(i),l,m,n);

  Kg(dofe,dofe) = Kg(dofe,dofe) + ke;


end

%%calculating displacements

Ug(freedof,1) = Kg(freedof,freedof)\Fg(freedof,1);

%stress calculation
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

fprintf('\n--- Stress Results ---\n');
fprintf('Max Tensile Stress:     %.4f Pa  -> Element %d\n', max_tensile, elem_tensile);
fprintf('Max Compressive Stress: %.4f Pa  -> Element %d\n', max_compressive, elem_compressive);

%visualization
scale = 100;

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

%Calculating reaction forces at the supports

Rg = Kg * Ug - Fg;


F_reaction = Rg(1:12,1);

