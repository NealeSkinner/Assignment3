%%%for part 3 i will use the FD method along with monte carlo. 


x= 3;           % x dimension
y=2;            % y dimension
%%%%% ratio for L/W is 3/2

dx = .1;      
dy = .1;
nx = x/dx;
ny = y/dy;
%%create matrices for the F and G values. 
G = sparse(nx*ny,nx*ny);        %large matrix of mostly 0's
F = zeros(nx*ny,1);
%step through the each point and change the values of the matrix that are
%equivalent to n or the boundaries
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        if i == 1 
            G(n,n) = 1;
            F(n) = 1;
         elseif i == nx 
            G(n,n) = 1;
        elseif j == 1 
            G(n, n) = -3;
            G(n, n+1) = 1;
            G(n, n+ny) = 1;
            G(n, n-ny) = 1;
        elseif j == ny
            G(n, n) = -3;
            G(n, n-1) = 1;
            G(n, n+ny) = 1;
            G(n, n-ny) = 1;
        else
            G(n, n) = -4;
            G(n, n+1) = 1;
            G(n, n-1) = 1;
            G(n, n+ny) = 1;
            G(n, n-ny) = 1;   
        end
    end
end
%%%%rearrange the equation GV=F to find V,
Volt = G\F;
%%%create a matrix/ map that holds the voltages in each cell
Voltmap = zeros(nx,ny,1);
%for loop that changes cells of voltmap to the given voltage
for i = 1:nx
    for j = 1:ny
         n = j+(i-1)*ny;
         Voltmap(i,j) = Volt(n);
    end
end
% %%use surf to 3d map the linear voltage 
% figure(1)
% plot(Voltmap);
% title('Voltage Plot from FD');
% xlabel('Bounds');
% % ylabel('y');
% ylabel('Voltage');
% % view(-100,10)

%%%%%%%%%%%%%%%%%%%%%part b%%%%%
G = sparse(nx*ny,nx*ny);        %recreate the sparse matrix of G
F = zeros(nx*ny,1);             %recreate the zeros matrix of F
%need to get answer for n, nym, nyp, nxm, and nxp for each iteration
%iterate through each point and change the values of G and F for all points.  
for i = 1:nx
    for j = 1:ny
        
        n = i + (j - 1) * nx;       
        nym = i + (j - 2) * nx; 
        nyp = i + j * nx; 
        nxm = i - 1 + (j - 1) * nx; 
        nxp = i + 1 + (j - 1) * nx; 
        
        if i == 1
            G(n, n) = 1; 
            F(n) = 1;
        elseif i == nx
            G(n, n) = 1; 
            F(n) = 1;
        elseif j == 1
            G(n, n) = 1; 
        elseif j == ny
            G(n, n) = 1; 
        else
            G(n, n) = -4; 
            G(n, nxm) = 1; 
            G(n, nxp) = 1; 
            G(n, nym) = 1; 
            G(n, nyp) = 1; 
            F(n) = 0;
        end
    end
end
%create a new voltage map using matrix 
Voltmap2 = zeros(nx, ny); 
%rearrange equation to get V
V = G\F; 
%iterate through map and change the values for their proper voltage
for i = 1:nx
    for j = 1:ny
        n = i + (j - 1)*nx;
        Voltmap2(i, j) = V(n); 
    end
end
 
figure(2)
surf(Voltmap2)
title('Electrostatic Potential: Numerical')
xlabel('x')
ylabel('y')
[Ex, Ey]=gradient(Voltmap2);
figure(3)
quiver(-Ex,-Ey);
title('Potential using quiver')

% sum = 0; 
% %analytical solution to the voltage map, not using n.
% for i = 1:nx
%     for j = 1:ny
%         for n = 1:2:(nx*ny)
%             sum = sum + ((1/n) * cosh(n*pi*i/dx) * sin(n*pi*j/dx)) / cosh(n*pi*dy/dx);
%         end
%     end
%     Voltmap2(i, j) = (4 * sum ) / pi ; 
% end
% 
% figure(3)
% surf(Voltmap2)
% title('Electrostatic Potential: Analytical')
% xlabel('x')
% ylabel('y')
