% Load Flow Solution using N-R Method 

clc
clear all
%% Data 

% Line Data 
%            From   To    R     X      Y/2
linedata = [ 1      2     0   0.1     0.01
             1      3     0   0.1     0.01
             2      3     0   0.1     0.01];

% Bus Data
%         Bus No   Bus Type    Pg        Qg    Pd       Qd      |V|    delta    Qmin    Qmax      
busdata = [ 1       1          0          0     0        0       1       0         0       0
            2       2          0.6661     0     0        0       1.05    0         0.2     2
            3       3          0          0     2.8653   1.2244  1       0         0       0  ];
        
%  Bus Type  Slack =1, PV = 2, PQ = 3

%% Y Bus Formation
R = linedata(:,3); X=linedata(:,4); B=1i*linedata(:,5);
Z = R+1i*X; Y = 1./Z;
nline= length(linedata(:,1));
nbus = max(max(linedata(:,1),linedata(:,2)));  % Total Number of buses
ng = length(find(busdata(:,2)==2));  % Total number of generator buses

Ybus = zeros(nbus,nbus);
for k=1:nline
    
    % Off Diagaonal Elements
    Ybus(linedata(k,1),linedata(k,2)) = Ybus(linedata(k,1),linedata(k,2)) - Y(k);
    Ybus(linedata(k,2),linedata(k,1)) = Ybus(linedata(k,1),linedata(k,2));
    
    % Diagonal Elements
    Ybus(linedata(k,1),linedata(k,1)) = Ybus(linedata(k,1),linedata(k,1))+ Y(k) + B(k); 
    Ybus(linedata(k,2),linedata(k,2)) = Ybus(linedata(k,2),linedata(k,2))+ Y(k) + B(k);
end
Ymag = abs(Ybus); theta = angle(Ybus);
%% Bus Data collection
type = busdata(:,2);
Pg = busdata(:,3); Qg = busdata(:,4); Pd = busdata(:,5); Qd = busdata(:,6); Qmin=busdata(:,9); Qmax = busdata(:,10);
Vmag = busdata(:,7); delta = busdata(:,8); 
V = Vmag.*(cos(delta) + 1i*sin(delta));
P_sch = Pg-Pd; Q_sch = Qg-Qd;  accuracy = 1; 
%% Iteration
iter = 1;
while accuracy >=0.001 && iter < 10 
% Since Bus 1 is a slack bus, We need not calculate its voltage.


%% Calculation of Mismatch Vector
for i=2:nbus
% Calculation of P estimate 
    P_cal(i) = 0; 
    Q_cal(i) = 0;
        for n=1:nbus
        P_cal(i) = P_cal(i) + Vmag(i)*Vmag(n)*Ymag(i,n)*cos(theta(i,n)+delta(n) - delta(i));
        Q_cal(i) = Q_cal(i) - Vmag(i)*Vmag(n)*Ymag(i,n)*sin(theta(i,n)+delta(n) - delta(i));
        end
        
        %% Q Limit checking for PV buses
        if Qmax(i) ~=0
            if Q_cal(i) > Qmax(i)
                Q_cal(i) = Qmax(i);
                busdata(i,2) = 3; % PV will become PQ temporarily
            elseif Q_cal(i) < Qmin(i)
                Q_cal(i) = Qmin(i);
                busdata(i,2) = 3;  % PV will become PQ temporarily
            else
                busdata(i,2) = 2; % PV will be restored
                Vmag(i) = busdata(i,7);
            end
        end
        
end


DP = P_sch (2:nbus) - P_cal(2:nbus)';
DQ = Q_sch ([find(busdata(:,2)==3)]) - Q_cal([find(busdata(:,2)==3)])';


%% Calcualtion of Jacobian Matrix 

% J1 
J1 = zeros(nbus,nbus);
for i=1:nbus
    for n=1:nbus
        if n~=i
            J1(i,i) = J1(i,i) + Vmag(i)*Vmag(n)*Ymag(i,n)*sin(theta(i,n)+delta(n)-delta(i));           
            J1(i,n) = - Vmag(i)*Vmag(n)*Ymag(i,n)*sin(theta(i,n)+delta(n)-delta(i));
            J1(n,i) = J1(i,n);
        end
    end
end
J11 = J1([find(busdata(:,2)~=1)], [find(busdata(:,2)~=1)]);


% J2
J2 = zeros(nbus,nbus);
for i=1:nbus
    for n=1:nbus %[find(busdata(:,2)==3)]
        if n~=i
            J2(i,i) = J2(i,i) + Vmag(n)*Ymag(i,n)*cos(theta(i,n)+delta(n)-delta(i));           
            J2(i,n) =  Vmag(i)*Ymag(i,n)*cos(theta(i,n)+delta(n)-delta(i));
            J2(n,i) = J2(i,n);
        else
            J2(i,i) = J2(i,i) + 2*Vmag(i)*Ymag(i)*cos(theta(i,i));
        end
    end
end
J22 = J2([find(busdata(:,2)~=1)], [find(busdata(:,2)==3)]);

% J3
J3 =  zeros(nbus,nbus);
for i=1:nbus %[find(busdata(:,2)==3)]
    for n=1:nbus
        if n~=i
            J3(i,i) = J3(i,i) + Vmag(i)*Vmag(n)*Ymag(i,n)*cos(theta(i,n)+delta(n)-delta(i));
            J3(i,n) = -Vmag(i)*Vmag(n)*Ymag(i,n)*cos(theta(i,n)+delta(n)-delta(i));
            J3(n,i) = J3(i,n);
        end
    end
end
J33 =J3([find(busdata(:,2)==3)], [find(busdata(:,2)~=1)]);


% J4
J4 = zeros(nbus,nbus);
for i = 1:nbus %[find(busdata(:,2)==3)]
   
        for n=1:nbus
            if n == i
                J4(i,i) = J4(i,i) -2*Vmag(i)*Ymag(i,i)*sin(theta(i,i));
            else
                J4(i,i) = J4(i,i) - Vmag(n)*Ymag(i,n)*sin(theta(i,n)+delta(n)-delta(i));
            end
        end
   
end

 J44 = J4([find(busdata(:,2)==3)],[find(busdata(:,2)==3)]);


J = [J11 J22 ; J33 J44 ];
%%  Calcualtion of Correction Vector

DF = [DP;DQ];

DX = J\DF;

delta([find(busdata(:,2)~=1)]) = delta([find(busdata(:,2)~=1)])+ DX(1:length(find(busdata(:,2)~=1)));
Vmag([find(busdata(:,2)==3)]) = Vmag([find(busdata(:,2)==3)]) + DX(length([find(busdata(:,2)~=1)])+1:length(DX));

accuracy = norm(DF);
iter = iter+1;
end
%% To find Slack bus power
for n=1:nbus
    Pg(1) = Pg(1) + Vmag(1)*Vmag(n)*Ymag(1,n)*cos(theta(1,n)+delta(n) - delta(1));
    Qg(1) = Qg(1) - Vmag(1)*Vmag(n)*Ymag(1,n)*sin(theta(1,n)+delta(n) - delta(1));
end

%% To find Reactive power of PV buses
for i=[find(busdata(:,2)==2)]
    for n=1:nbus
        Qg(i) = Qg(i) - Vmag(i)*Vmag(n)*Ymag(i,n)*sin(theta(i,n)+delta(n) - delta(i));      
    end
end


%% To find losses
PL = sum(Pg) - sum(Pd);
QL = sum(Qg) - sum(Qd);

disp('Number of iterations')
disp(iter)
disp('Vm in p.u.') 
disp(Vmag)
disp('delta in rad')
disp(delta)
disp('Pg of slack bus in p.u')
disp(Pg(1))
disp('Qg of slack bus in p.u')
disp(Qg(1))
disp('Real Power loss in p.u')
disp(PL)
disp('Reactive power loss in p.u')
disp(QL)
