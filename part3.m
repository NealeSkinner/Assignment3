clc
clear

T=300;                  %Temp in K
K=1.38e-23;             %Boltsmann constant
Tmn=0.2e-12;            %mean time between collisions
Mo=9.11e-31;            %rest mass
Mn=0.26*Mo;             %effective mass of electrons
L=200e-09;               %Length of region
W=100e-09;               %Width of region
PopE=15;                %number of particles
Vth=sqrt((K*T)/(Mn));   %Thermal velocity    
TimeS=15e-15;           %time step of 15ns
iterations=100;
q = 1.60217653e-19;             % electron charge
AtomSpacing = 0.5430710e-9;         
stdD=Vth/sqrt(2);               %std. Dev.
concentration = 10e-15;
Volts=0.1;
ex=Volts/L;
%%%%%%%%%%%%%%%%%%%Electron Modeling%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CurrentVector = zeros(iterations,1);
        TempArray = zeros(iterations,1); % temp array
        temp = (1:1:iterations);
        Pscat = 1 - exp(-TimeS/Tmn);    %probability of scatter
        Tar=zeros(1,PopE)+T;    %Array of Temp
        MFP=Tmn * Vth;          %Mean Free Path
        Angle = rand(PopE,1)*2*pi;   %random angle
        Pos = [rand(PopE,1)*L rand(PopE,1)*W Vth*cos(Angle) Vth*sin(Angle)];
        OutBounds=((Pos(:,1)>L/4)&(Pos(:,1)<3*L/4))&((Pos(:,2)<W/4)|(Pos(:,2)>3*W/4));      %logical array for all values inside boxes.
        initialX=Pos(:,1);
        initialY=Pos(:,2);
        VprobS = makedist('Normal', 'mu', 0, 'sigma', Vth);        %Velocity Distribution
        AccX = zeros(PopE,1);
        AccX(:,1)= ex * (q/Mn);
     %the electrons initialized inside the boxes need to be randomly placed
     %again until no atoms initialize inside the box
       while (sum(OutBounds)>1)
           initialX(OutBounds)=initialX(OutBounds) .*randn();                           %randomize the electrons initialized inside the box
           initialY(OutBounds)=initialY(OutBounds) .*randn();
           OutBounds=(initialX>L/4 & initialX<3*L/4)&(initialY<W/4 |initialY>3*W/4);    %check once again for electrons initialized inside the box.
       end
colorstring= rand(1,PopE);          %array of color shades, will be used to plot so that electrons keep the same colors
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
    
    OverX=NextX > L;                 %check for values in next position that are greater than bounds
    NextX(OverX)=NextX(OverX)-L;         %subtract from the bound to reflect
  
    UnderX=NextX < 0;                   %check for values less than bounds
    NextX(UnderX)= NextX(UnderX) + L;   %reflect bound
        
    OverY=NextY > W;                    %check for values greater then bounds
    NextY(OverY)= 2*W -NextY(OverY) ;   %subtract from bounds to reflect
    Pos(OverY,4)=- Pos(OverY,4);
    
    UnderY= NextY < 0;                  %check for under bounds
    NextY(UnderY)=- NextY(UnderY);      %reflect
    Pos(UnderY,4)=- Pos(UnderY,4);
    
     Xless=((NextX <3*L/4 & NextX>L/4) &(NextY<W/3 | NextY>3*W/4)); %logical array for values that are going into box
     % below is where we check bounds of boxes
       if sum(Xless) ~=0                                            %if there are any values in the Xless array we have an electron on the boundary of the box
        if (initialX(Xless)<3*L/4 & initialX(Xless)>L/4)            %check for electron coming from the Y plane, reverse velocity in Y direction if so
            Pos(Xless,4)=-Pos(Xless,4);     
        else                                                        %if its not coming from Y plane that means its in x plane, reverse X velocity
            Pos(Xless,3)=-Pos(Xless,3); 
        end
       end
                                 
    avgV=mean(sqrt(Pos(:,3).^2 + Pos(:,4).^2));
    MFP=(avgV*i)/1000000;
    figure(1)
    scatter(initialX,initialY,PopE,colorstring, '.' );      %plot the initial positions that are determined after we know what type of reflection
    hold on
    scatter(NextX, NextY,PopE,colorstring,'.' );            %plot the next position.
    axis([0 L 0 W]);
    title (['Trajectories with MFP ', num2str(MFP),'e^6']);
    
    initialX=NextX;                     %the next time we timestep we need to save the coordinates we just timestepped as the initial values to be able to time step them
    initialY=NextY;
end

    %Plot of Current vs Time
    figure(2);
    plot(temp, CurrentVector)
    title('Average Current')
    xlabel('Bottle Neck')
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
    figure(3)
    surf(X,Y,Dense)
    colorbar
    title 'Electron Density Map';
    zlabel 'Number of Electrons per Grid Point';
    ylabel 'Y Coordinate';
    xlabel 'X coordinate';
    %temperature map
%     figure(4)
%     surf(X,Y,Temp)
%     colorbar
%     title 'Temperature Density Map';
%     zlabel 'Temperature per Grid Point';
%     ylabel 'Y Coordinate';
%     xlabel 'X coordinate';