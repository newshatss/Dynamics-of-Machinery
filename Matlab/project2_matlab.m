
%kinematics --> lines 6-128
%kinetics --> lines 132-232
%plots --> lines 236-305

%kinematics
t = 0:0.002:1;
L = length(t);
alpha2 = 10;
omega2_0 = 20;
omega2 = alpha2*t + omega2_0;
theta1 = 0.271;
theta2 = (alpha2*t.^2)/2 + omega2_0*t;
theta6 = pi;
theta7 = pi/2;

%Lengths in cm
r1 = 21.45;
r2 = 12.7;
r4 = 60.45;
r5 = 90;
r7 = 36.3 + r1*sin(theta1);

a = (r1*sin(theta1)+r2*sin(theta2))./(r1*cos(theta1)+r2*cos(theta2));
theta3 = atan(a);
clear a;

r3 = (r1*cos(theta1)+r2*cos(theta2))./cos(theta3);

r3_dot = zeros(1,L);
omega3 = zeros(1,L);
for n=1:L
    th3 = theta3(n);
    th2 = theta2(n);
    o2 = omega2(n);
    
    A = [cos(th3) -r3(n)*sin(th3); sin(th3) r3(n)*cos(th3)];
    B = [-r2*o2*sin(th2); r2*o2*cos(th2)];
    
    x = A\B;
    
    r3_dot(n) = x(1);
    omega3(n) = x(2);
    
    clear th3 th2 o2 A B;
end

r3_ddot = zeros(1,L);
alpha3 = zeros(1,L);
for n=1:L
    th3 = theta3(n);
    th2 = theta2(n);
    r3d = r3_dot(n);
    o2 = omega2(n);
    o3 = omega3(n);
    
    A = [-cos(th3) r3(n)*sin(th3); sin(th3) r3(n)*cos(th3)];
    B = [r2*alpha2*sin(th2) + r2*o2^2*cos(th2) - 2*r3d*o3*sin(th3) - r3(n)*o3^2*cos(th3);
        r2*alpha2*cos(th2) - r2*o2^2*sin(th2) - 2*r3d*o3*cos(th3) + r3(n)*o3^2*sin(th3)];
    
    x = A\B;
    
    r3_ddot(n) = x(1);
    alpha3(n) = x(2);
    
    clear th3 th2 r3d o3 o2 A B;
end


theta4 = theta3;
omega4 = omega3;
alpha4 = alpha3;


r6 = zeros(1,L);
theta5 = zeros(1,L);
for n=1:L
    th4 = theta4(n);
    
    func = @(x) ([x(1)*cos(theta6) + r7*cos(theta7) - r4*cos(th4) - r5*cos(x(2));
        x(1)*sin(theta6) + r7*sin(theta7) - r4*sin(th4) - r5*sin(x(2))]);
    
    y = fsolve(func,[1,1]);
    clc;
    
    r6(n) = y(1);
    theta5(n) = y(2);
    
    clear th4 y;
end

r6_dot = zeros(1,L);
omega5 = zeros(1,L);
for n=1:L
    th4 = theta4(n);
    th5 = theta5(n);
    o4 = omega4(n);
    
    A = [cos(theta6) r5*sin(th5); -sin(theta6) r5*cos(th5)];
    B = [-r4*o4*sin(th4); -r4*o4*cos(th4)];
    
    x = A\B;
    
    r6_dot(n) = x(1);
    omega5(n) = x(2);
    
    clear th4 th5 o4 A B;
end

r6_ddot = zeros(1,L);
alpha5 = zeros(1,L);
for n=1:L
    a4 = alpha4(n);
    th4 = theta4(n);
    o4 = omega4(n);
    th5 = theta5(n);
    o5 = omega5(n);
    
    A = [-cos(theta6) -r5*sin(th5); sin(theta6) -r5*cos(th5)];
    B = [r4*a4*sin(th4) + r4*o4^2*cos(th4) + r5*o5^2*cos(th5);
        r4*a4*cos(th4) - r4*o4^2*sin(th4) - r5*o5^2*sin(th5)];
    
    x = A\B;
    
    r6_ddot(n) = x(1);
    alpha5(n) = x(2);
    
    clear a4 th4 o4 th5 o5 A B x;
end


%kinetics
m2 = 0.46;
m3 = 0.23;
m4 = 3.987;
m5 = 4.487;
m6 = 0.23;
I2 = 9.288*10^(-4);
I3 = 6.377*10^(-5);
I4 = 0.17;
I5 = 0.312;
I6 = 6.377*10^(-5);

rg2 = 4.8195;
rg4 = 32.3804;
rg5 = r5/2;

ag2x = (rg2.*alpha2.*cos(theta2+pi/2) - rg2.*omega2.^2.*cos(theta2))*10^(-2);
ag2y = (rg2.*alpha2.*sin(theta2+pi/2) - rg2.*omega2.^2.*sin(theta2))*10^(-2);

ag3x = ((r3_ddot-r3.*omega3.^2).*cos(theta3) + (2*r3_dot.*omega3+r3.*alpha3).*cos(theta3+pi/2))*10^(-2);
ag3y = ((r3_ddot-r3.*omega3.^2).*sin(theta3) + (2*r3_dot.*omega3+r3.*alpha3).*sin(theta3+pi/2))*10^(-2);

ag4x = (rg4.*alpha4.*cos(theta4+pi/2) - rg4.*omega4.^2.*cos(theta4))*10^(-2);
ag4y = (rg4.*alpha4.*sin(theta4+pi/2) - rg4.*omega4.^2.*sin(theta4))*10^(-2);

ag5x = (r4.*alpha4.*cos(theta4+pi/2) - r4.*omega4.^2.*cos(theta4) + rg5.*alpha5.*cos(theta5+pi/2) - rg5.*omega5.^2.*cos(theta5))*10^(-2);
ag5y = (r4.*alpha4.*sin(theta4+pi/2) - r4.*omega4.^2.*sin(theta4) + rg5.*alpha5.*sin(theta5+pi/2) - rg5.*omega5.^2.*sin(theta5))*10^(-2);

ag6x = -r6_ddot*10^(-2);

%link 6
F56x = m6*ag6x;

%link 5
F45x = zeros(1,L);
F45y = zeros(1,L);
F56y = zeros(1,L);
for n=1:L
    A = [1 0 0; 0 1 -1; -r5/2*sin(theta5(n)) r5/2*cos(theta5(n)) r5/2*cos(theta5(n))];
    B = [m5*ag5x(n) + F56x(n); m5*ag5y(n); I5*alpha5(n)*10^2 + F56x(n)*r5/2.*sin(theta5(n))]; 
    
    x = A\B;
    
    F45x(n) = x(1);
    F45y(n) = x(2);
    F56y(n) = x(3);
    
    clear A B x;
end

%link 4,3
F14x = zeros(1,L);
F34 = zeros(1,L);
F14y = zeros(1,L);
a = zeros(1,L);
F23x = zeros(1,L);
F23y = zeros(1,L);
for n=1:L
    th4 = theta4(n);
    func = @(x) ([x(1)+x(2)*cos(th4+pi/2)-F45x(n)-ag4x(n)*m4;
        x(3)+x(2)*sin(th4+pi/2)-F45y(n)-ag4y(n)*m4;
        F45x(n)*(r4-rg4)*sin(th4)+x(1)*rg4*sin(th4)-F45y(n)*(r4-rg4)*cos(th4)-x(2)*(rg4-r3(n)+x(4))-x(3)*rg4*cos(th4)-I4*alpha4(n)*10^2;
        x(5)+x(2)*cos(pi/2-th4)-ag3x(n)*m3;
        x(6)-x(2)*sin(th4)-ag3y(n)*m3;
        x(2)*x(4)-I3*alpha3(n)*10^2;]);

    if n==1
        y = fsolve(func,[1000,1000,1000,0,1000,1000]);
        clc;
    else
        y = fsolve(func,[y(1),y(2),y(3),y(4),y(5),y(6)]);
        clc;
    end

    F14x(n) = y(1);
    F34(n) = y(2);
    F14y(n) = y(3);
    a(n) = y(4);
    F23x(n) = y(5);
    F23y(n) = y(6);
    
    clear th4;
end

%link2
F12x = zeros(1,L);
F12y = zeros(1,L);
T12 = zeros(1,L);
for n=1:L
    A = [1 0 0; 0 1 0; rg2*sin(theta2(n)) rg2*cos(theta2(n)) 1];
    B = [m2*ag2x(n)+F23x(n); m2*ag2y(n)+F23y(n); I2*alpha2*10^2-F23x(n)*(r2-rg2)*sin(theta2(n))-F23y(n)*(r2-rg2)*cos(theta2(n))];
    
    x = A\B;
    
    F12x(n) = x(1);
    F12y(n) = x(2);
    T12(n) = x(3);
    


    clear A B x;
end


%Plots
% plot(t,omega4)
% xlabel('Time(s)')
% ylabel('Angular Velocity(rad/s)')
% title('Link 4')
% plot(t,omega5)
% xlabel('Time(s)')
% ylabel('Angular Velocity(rad/s)')
% title('Link 5')
% plot(t,alpha4)
% xlabel('Time(s)')
% ylabel('Angular Acceleration(rad/s^2)')
% title('Link 4')
% plot(t,alpha5)
% xlabel('Time(s)')
% ylabel('Angular Acceleration(rad/s^2)')
% title('Link 5')
% plot(r4*cos(theta4),r4*sin(theta4))
% xlabel('X(cm)')
% ylabel('Y(cm)')
% title('Path of Point B')
% plot(t,abs(r4*omega4))
% xlabel('Time(s)')
% ylabel('Velocity(cm/s)')
% title('Point B(velocity magnitude')
% plot(t,sqrt((r4*alpha4).^2+(r4*omega4.^2).^2)
% xlabel('Time(s)')
% ylabel('Acceleration(cm/s^2)')
% title('Point B(acceleration magnitude)')
% plot(t,r3_ddot)
% xlabel('Time(s)')
% ylabel('Acceleration')
% title('r3 double dot')
% plot(t,r6_dot)
% xlabel('Time(s)')
% ylabel('Velocity(cm/s)')
% title('r6 dot')
% plot(t,r6_ddot)
% xlabel('Time(s)')
% ylabel('Acceleration (rad/s^2)')
% title('r6 double dot')

% plot(t,sqrt(F12x.^2+F12y.^2))
% xlabel('Time(s)')
% ylabel('Force(N)')
% title('Joint Between 1 & 2')
% plot(t,sqrt(F23x.^2+F23y.^2))
% xlabel('Time(s)')
% ylabel('Force(N)')
% title('Joint Between 2 & 3')
% plot(t,abs(F34))
% xlabel('Time(s)')
% ylabel('Force(N)')
% title('Normal Force Between 3 & 4(N)')
% plot(t,sqrt(F14x.^2+F14y.^2))
% xlabel('Time(s)')
% ylabel('Force(N)')
% title('Joint Between 1 & 4')
% plot(t,sqrt(F45x.^2+F45y.^2))
% xlabel('Time(s)')
% ylabel('Force(N)')
% title('Joint Between 4 & 5')
% plot(t,sqrt(F56x.^2+F56y.^2))
% xlabel('Time(s)')
% ylabel('Force(N)')
% title('Joint Between 5 & 6')
% plot(t,abs(T12))
% xlabel('Time(s)')
% ylabel('Torque(N.m)')
% title('Input Torque')
% grid on