# This script solves Poisson equation using the finite element method. See the documentation for the theory and implementation.
L = 1
N = 50
X = seq(0,L,len=N+1)
m = length(X); n=length(X);
x = matrix(rep(X,each=n),nrow=n);
y = matrix(rep(X,m),nrow=n)
x = c(x)
y = c(y)
GNodes = (N+1)^2
NTRI = 2*N^2
LocNodes = matrix(0,NTRI,3)

for (i in 1:N){
	for (j in 1:N){
		    LocNodes[i+2*(j-1)*N,1] = i+(j-1)*(N+1)
        LocNodes[i+2*(j-1)*N,2] = i+j*(N+1)
        LocNodes[i+2*(j-1)*N,3] = (i+1)+(j-1)*(N+1)
        LocNodes[i+N+2*(j-1)*N,1] = i+1+j*(N+1)
        LocNodes[i+N+2*(j-1)*N,2] = (i+1)+(j-1)*(N+1)
        LocNodes[i+N+2*(j-1)*N,3] = i+j*(N+1)
    }
}

Sparsity <- matrix(0, GNodes, GNodes)
LoadVect = matrix(0, GNodes, 1)


for (n in 1:NTRI)
	{
	r1 = matrix(c(x[LocNodes[n,1]],y[LocNodes[n,1]]),nrow=2, byrow=FALSE)
  r2 = matrix(c(x[LocNodes[n,2]],y[LocNodes[n,2]]),nrow=2, byrow=FALSE);
  r3 = matrix(c(x[LocNodes[n,3]],y[LocNodes[n,3]]),nrow=2, byrow=FALSE);
  J = matrix(c(r2[1]-r1[1],r2[2]-r1[2]
    ,r3[1]-r1[1],r3[2]-r1[2]), nrow=2, byrow=TRUE) 
  Astiff = (1/(2*det(J)))*matrix(c(sum((r2-r3)*(r2-r3))
          ,sum((r2-r3)*(r3-r1)),sum((r2-r3)*(r1-r2))
          ,sum((r2-r3)*(r3-r1)),sum((r3-r1)*(r3-r1))
          ,sum((r3-r1)*(r1-r2)),sum((r2-r3)*(r1-r2))
          ,sum((r3-r1)*(r1-r2)),sum((r1-r2)*(r1-r2))), nrow=3, byrow=TRUE);
  for (i in 1:3){
      for (j in 1:3){
        Sparsity[LocNodes[n,i],LocNodes[n,j]]=Sparsity[LocNodes[n,i],LocNodes[n,j]]+Astiff[i,j]
      }
    }
  ksi = 1/3;
  eta = 1/3;
  xx = (1-ksi-eta)*r1[1]+ksi*r2[1]+eta*r3[1];
  yy = (1-ksi-eta)*r1[2]+ksi*r2[2]+eta*r3[2];
  
  # These three lines will construct the local load vector if the right handside of the Poisson equation is 1. For example if the problem reads as $ laplace u = 1 $.
  F[1] = (1-ksi-eta)*det(J)*1/2;
  F[2] = (ksi)*det(J)*1/2;
  F[3] = (eta)*det(J)*1/2;
  
  # F[1] = (1-ksi-eta)*5*pi^2*sin(pi*xx)*sin(2*pi*yy)*det(J)*1/2;
  # F[2] = (ksi)*5*pi^2*sin(pi*xx)*sin(2*pi*yy)*det(J)*1/2;
  # F[3] = (eta)*5*pi^2*sin(pi*xx)*sin(2*pi*yy)*det(J)*1/2;
  
  # F[1] = pi^2*((1-ksi-eta)*det(J)*(1/2))*(yy*sin(2*pi*xx)+xx*sin(2*pi*yy));
  # F[2] = pi^2*((ksi)*det(J)*(1/2))*(yy*sin(2*pi*xx)+xx*sin(2*pi*yy));
  # F[3] = pi^2*((eta)*det(J)*(1/2))*(yy*sin(2*pi*xx)+xx*sin(2*pi*yy));

  
  for (i in 1:3){
    LoadVect[LocNodes[n,i]] = LoadVect[LocNodes[n,i]]+ F[i]
    }
	}
for (i in 1:GNodes){
if (x[i]==0 | y[i]==0 | x[i]==L | y[i]==L){
    LoadVect[i] = 0
    Sparsity[i,] = 0
    Sparsity[i,i] = 1
     }
}



U = solve(Sparsity,LoadVect)






# The following commands will return out put that can be processed on any external software like MATLAB for visualisation. The file 'plotsolution.m' is an example of how to obtain visualisation of the solution surface.  


write.table(x,file="./xcoordates.txt",row.names=FALSE,col.names=FALSE)
write.table(y,file="./ycoordates.txt",row.names=FALSE,col.names=FALSE)
write.table(U,file="./Solutions.txt",row.names=FALSE,col.names=FALSE)
write.table(LocNodes,file="./triangles.txt",row.names=FALSE,col.names=FALSE)



