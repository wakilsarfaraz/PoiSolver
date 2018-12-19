
coordx = importdata('xcoordates.txt');
coordy = importdata('ycoordates.txt');
triangles = importdata('triangles.txt');
numsol = importdata('Solutions.txt');
figure(1)
trisurf(triangles,coordx,coordy,numsol)
shading interp
title('Finite element solution','fontsize',20)
xlabel('x','fontsize',16)
ylabel('y','fontsize',16)
zlabel('U','fontsize',16)

