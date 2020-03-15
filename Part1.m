clear
clc
%%%%Assignment 3 question 1
%%%%most of this comes from assignment 1. Note that my density graphs arent
%%%%that excited from lack of electrons. you could add more electrons
%%%%using PopE but then the electron modelling is unreadable due to the
%%%%method i have used. 
T=300;                          %Temp in K
K=1.38e-23;                     %Boltsmann constant
Tmn=0.2e-12;                    %mean time between collisions
Mo=9.11e-31;                    %rest mass
Mn=0.26*Mo;                     %effective mass of electrons
L=200e-09;                      %Length of region
W=100e-09;                      %Width of region
PopE=5000;                        %number of particles  
TimeS=15e-15;                   %time step of 15ns
q = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;           % Dirac constant
h = hb * 2 * pi;                % Planck constant
eps = 8.854187817e-12;          % vacuum permittivity
mu = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                  % speed of light
g = 9.80665;                    %metres (32.1740 ft) per s²
am = 1.66053892e-27;

iterations = 30;
AtomSpacing = 0.5430710e-9;
Vth=sqrt((K*T)/(Mn));           %Thermal velocity 
stdD=Vth/sqrt(2);               %std. Dev.
concentration = 10e-15;
Volts=0.1;
ex=Volts/L;

%%%%Electron Modeling%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MFP=Tmn * Vth;          %Mean Free Path
%fprintf('the Mean Free Path is %d\n', MFP);
    CurrentVector = zeros(iterations,1);
    TempArray = zeros(iterations,1); % temp array
    temp = (1:1:iterations);
    Angle = rand(PopE,1)*2*pi;   %random angle
    Pscat = 1 - exp(-TimeS/Tmn);    %probability of scatter
    Pos = [rand(PopE,1)*L rand(PopE,1)*W Vth*cos(Angle) Vth*sin(Angle)];        %initializing position array. column 1 holds X, 2 holds y, 3 holds velocity in X and 4 holds velocity in Y
    DVelocity = zeros(1,PopE); % Drift Velocity
    EMobility = zeros(1,PopE); % Electron Mobility
    VprobS = makedist('Normal', 'mu', 0, 'sigma', sqrt(K*T/Mn));        %Velocity Distribution
    %Calculate the Acceleration
    AccX = zeros(PopE,1);
    AccX(:,1)= ex * (q/Mn);
    initialX=Pos(:,1);                  %array of the initial X values in Position array
    initialY=Pos(:,2);                 %array of the initial Y values in Position array
    colorstring=rand(PopE,1);           %array of values corresponding to shades of color, will be used in plotting
for i= 1 : iterations
    OverV = rand(PopE,1) < Pscat ;                              %check for probability of scatter
    Pos(OverV,3:4)= random(VprobS,[sum(OverV),2]);
    Pos(:,3)=Pos(:,3) + AccX.*TimeS;        %take into account Acceleration. Acc*DT=Vel.
    NextX=initialX + Pos(:,3) * TimeS;      %the next X position will be the initial plus a time step of velocity in X direction
    NextY=initialY + Pos(:,4) * TimeS;      %the next Y position will be the initial a plus a time step of velocity in Y.
    %calculate VRMS
    Vrms = sqrt((Pos(:,3) .^ 2) + (Pos(:,4) .^ 2));   
    VrmsArray =zeros(1,PopE);
    %calculate electron mobility
    EMobility=mean(Vrms);
    %calculate Drift Velocity
    DVelocity=EMobility*ex;
    %Calculte Current from charge, area, concentration and drift velocity.
    current = concentration * L * W * (sum(DVelocity)/PopE) * q;
    CurrentVector(i) = current;
    %calculate Temperature and fill array. 
    Temp = (sqrt(2) * (mean(Vrms)) ^2 *Mn) / K; 
    TempArray(1,i) = Temp;
    
    OverX=NextX > L;                     %check for values in next position that are greater than bounds
    NextX(OverX)=NextX(OverX)-L;         %subtract from the bound to reflect
  
    UnderX=NextX < 0;                   %check for values less than bounds
    NextX(UnderX)= NextX(UnderX) + L;   %reflect bound
    
    OverY=NextY > W;                    %check for values greater then bounds
    NextY(OverY)= 2*W -NextY(OverY) ;   %subtract from bounds to reflect
    Pos(OverY,4)=- Pos(OverY,4);        %reverse the velocity
    
    UnderY= NextY < 0;                  %check for under bounds
    NextY(UnderY)=- NextY(UnderY);      %reflect
    Pos(UnderY,4)=- Pos(UnderY,4);      %reverse the velocity
 

    figure(1)   
    scatter(initialX,initialY,PopE, colorstring,'.' );      %plot the coordinates of the electron before the move
    hold on
    scatter(NextX, NextY,PopE, colorstring,'.' );           %plot the coordinates of the electron after the move
    
    axis([0 L 0 W]);
    title (['Average Temperature ', num2str(Temp)]);
    
    initialX=NextX;
    initialY=NextY;
end


    %plot of Temperature vs time

    figure(2);
    plot(temp,TempArray)
    title('Average Temperature')
    xlabel('Time (s)')
    ylabel('Temperature')
    hold on

    %Plot of Current vs Time
    figure(3);
    plot(temp, CurrentVector)
    title('Average Current')
    xlabel('Time')
    ylabel('Current')
    hold on

    %for density, divide the L and W into a grid, step through grid points
    %to see how many electrons are in each grid point
    [X Y]=meshgrid(0:L/10:L,0:W/10:W);  %grid for stepping by 1/10 of width and length
    %created an 11 x 11 vectors
    Dense=zeros(11,11);                 %will hold the density per nm
    Temp=zeros(11,11);                  %holds temperature per nm
    Tcount=0;                           %initilaize to 0
    Dcount=0;
for i=1:10
        XB1=X(1,i);            %holds the value of first bound grid coordinate
        XB2=X(1,i+1);          %holds the value of the outside grid coordinate
        for j =1:10
            YB1=Y(j,1);        %holds the value of the first Y bound
            YB2=Y(j+1,1);      %holds the value of the last Y bound
        
            %check each frame
            for k=1:PopE
                %Check to see if particle is within the bounds stored above
                if((Pos(k,1)>XB1 & Pos(k,1)<XB2) & Pos(k,2)<YB2 & Pos(k,2)>YB1)
                    Tcount=Tcount+1;                                %if it is within the bound we update out temp counter
                    Dense(i,j)=Dense(i,j)+1;                        % we also add to the density array count
                    Dcount=Dcount+sqrt(Pos(k,3)^2+Pos(k,4)^2);      % denisty count. 
                    if(Tcount >0)
                        Temp(i,j)=Mn*(Dcount^2)/(Tcount)/K/2;         %the temp is divided for each electron 
                    end
                end
            end
            Dcount=0;
            Tcount=0;
        end
end
    %density map
    figure(4)
    surf(X,Y,Dense)
    colorbar
    title 'Electron Density Map';
    zlabel 'Number of Electrons per Grid Point';
    ylabel 'Y Coordinate';
    xlabel 'X coordinate';
    %temperature map
    figure(5)
    surf(X,Y,Temp)
    colorbar
    title 'Temperature Density Map';
    zlabel 'Temperature per Grid Point';
    ylabel 'Y Coordinate';
    xlabel 'X coordinate';
