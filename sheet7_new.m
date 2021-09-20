%% sheet 7
clear
close all
clc
%% Assignment 9 a)
load data.mat

figure(1)
plot(pl(:,2),T_m-T30(:,2))
grid on
xlabel('Temperature drop [°C]')
ylabel('P_l [W/m]')
%% b)
R_eff = (T_m-T30(:,2))./pl(:,2);

figure(2)
plot(t30/60/60/24/365.15,R_eff)
grid on
xlabel('Time [a]')
ylabel('R_{eff} [°Cm/W]')

disp(['The maximum R_{eff} of ' num2str(max(R_eff)) '°Cm/W occures in month ' num2str(round(t30(R_eff==max(R_eff))/60/60/24/365.15*12))])
%% c) + d) + e)
T_mean = sum(pl(:,2).*(T_m-T30(:,2)))/sum(pl(:,2))

R_eff = T_mean/pl(1,1)
R_eff_e = max(T_m-T30(:,2))/pl(1,1)
% return
% close all
%% Assignment 10 a)
K = 273.15;
depth = 1.5; %m
spacing = 1; %spacing
n = 16; %tubes
size = 1.5; % size*living area
w = 15; % m
r = 0.032/2; % m - radius of tube
d = 0.003; % m - thickness of walls
T_m = 10+K; %°C
pl = [0.7 2.7 10.2 19.7 29.5 35.5 32.3 28.0 23.8 13.4 6.3 3.2];
l = 1/3;
lambda = 1.5;
rho = 1800;
c = 1100;
kappa = lambda/(rho*c);
l_PE = 0.48; % W/mK
l_s = lambda; % W/mK

l_tot = n*(size/w); % m/m^2 
pl = pl*(l/l_tot); % Wm/m^2                                    
disp(['Total length of tubes: ' num2str(l_tot) ' m/m^2']) 

%% b)
R = 1/(2*pi*l_PE)*log(r/(r-d));
disp(['Borehole resistance: ' num2str(R) ' mK/K'])
r_a = r*((r-d)/r)^(l_s/l_PE);
disp(['Apparent radius: ' num2str(r_a) ' m/m^2'])

%% c) + d) + e)
t = ((1:numel(pl))-1)*60*60*24*365.25/12;
T_h = 35+K; % °C
n_rel = 0.5;
T = [T_m T_m T_m]; % K

s = [r_a (1:15)*spacing];
d = [sqrt(s.^2+(2*depth)^2)];
num = [1 2*(1-(1:n-1)/n)];

e = DownholeHeatExchanger(lambda,kappa,r);
for j = 1:numel(T)
    for i = 1:numel(t)
        p_th(i,j) = pl(i)/(1+(1/n_rel)*((T_h-T(i,j))/T(i,j)));
        e.setpower(t(i),p_th(i,j));
        T(i+1,j) = T_m-e.temperature(i*60*60*24*365.25/12,r_a);
        if j == 2
            T(i+1,j) = T_m-(e.temperature(i*60*60*24*365.25/12,r_a)...
                -e.temperature(i*60*60*24*365.25/12,2*depth));
        elseif j == 3
            for o = 1:numel(num)
                T_s(o) = num(o)*e.temperature(i*60*60*24*365.25/12,s(o));
                T_i(o) = num(o)*e.temperature(i*60*60*24*365.25/12,d(o));
            end
            T_s = sum(T_s);
            T_i = sum(T_i);
            T(i+1,j) = T_m-(T_s-T_i);
        end
    end
end

T = T(2:end,:)-K;

%% f)
load T_soil.mat
load Ppla.mat

T_soil = repmat(T_d(2:end-1,6),30,1)+K;
t = ((1:360)-1)*60*60*24*365.25/12;
pl = repmat(Power_p_area,30,1);
T_s_30 = T_soil(1)+K;

for i = 1:numel(t)
    p_th_s(i) = pl(i)/(1+(1/n_rel)*((T_h-T_s_30(i))/T_s_30(i)));
    e.setpower(t(i),p_th_s(i));
    for o = 1:numel(num)
        T_s(o) = num(o)*e.temperature(i*60*60*24*365.25/12,s(o));
        T_i(o) = num(o)*e.temperature(i*60*60*24*365.25/12,d(o));
    end
    T_s = sum(T_s);
    T_i = sum(T_i);
    T_s_30(i+1) = T_soil(i)-(T_s-T_i);
end

T_s_30 = T_s_30(2:end)-K;

figure(3)
plot(t30/60/60/24/365.25,T_s_30)
grid on
xlabel('Time [a]')
ylabel('Temperature [°C]')
%% g)
T_min = min(T_s_30);
T_0 = T_soil(1)+K;
pl_min = pl;

for j = 1:200
    f(j) = (mean(T_soil)-0)/(mean(T_soil)-T_min(j))-1e-2;
    pl_min = pl_min.*f(j);
    for i = 1:numel(t)
        p_th_0(i) = pl_min(i)/(1+(1/n_rel)*((T_h-T_0(i))/T_0(i)));
        e.setpower(t(i),p_th_0(i));
        for o = 1:numel(num)
            T_s(o) = num(o)*e.temperature(i*60*60*24*365.25/12,s(o));
            T_i(o) = num(o)*e.temperature(i*60*60*24*365.25/12,d(o));
        end
        T_s = sum(T_s);
        T_i = sum(T_i);
        T_0(i+1) = T_soil(i)-(T_s-T_i);
    end
    T_min(j+1) = min(T_0)-K;
    if T_min(j+1) >= 0
        T_0 = T_0(2:end)-K;
        break
    end
end

f_new = prod(f);
size_new = size*1/f_new;
disp(['New size: ' num2str(size_new) ' times living area'])

figure(4)
plot(t30/60/60/24/365.25,T_0)
grid on
xlabel('Time [a]')
ylabel('Temperature [°C]')

%% h)
load E.mat
hd = 50;
P_el = pl-p_th_s;
l_tot_min = n*(size_new/w); % m/m^2
pl_0 = pl*(l/l_tot_min); % Wm/m^2
P_el_min = pl_0-p_th_0;

k = 1:12:360; 
for i = 1:numel(k)
    P_el_a(i) = sum(P_el(k(i):k(i)+11));
    P_tot_a(i) = sum(pl(k(i):k(i)+11));
    P_el_min_a(i) = sum(P_el_min(k(i):k(i)+11));
%     P_tot_min_a(i) = sum(pl_min(k(i):k(i)+11));
end

COP_s = P_tot_a./P_el_a;
COP_min = P_tot_a./P_el_min_a;
El_s = hd./COP_s;
El_min = hd./COP_min;

figure(5)
plot(1:numel(El_s),El_s)
hold on
grid on
plot(1:numel(El_min),El_min)
plot(1:numel(El_a(2,:)),El_a(2,:))
xlabel('Time [a]')
ylabel('Electric energy [kWh/m^2a]')
legend(['Heat collector with size ' num2str(size) ' times living area'],...
    ['Heat collector with size ' num2str(size_new) ' times living area'],...
    'Downhole heat exchanger with heat pump and R_b = 0.1','Location','best')

disp(['The larger heat collector saves ' num2str(sum(El_s)-sum(El_min)) ' kWh/m^2a'])
disp(['Energy consumption for size = ' num2str(size) ': ' num2str(sum(El_s)) ' kWh/m^2a'])
disp(['Energy consumption for size = ' num2str(size_new) ': ' num2str(sum(El_min)) ' kWh/m^2a'])
disp(['Energy consumption from Ass 8: ' num2str(E_tot(2)) ' kWh/m^2a'])
% return
%%  i)
n_new = 31;
spacing_new = 0.5; % m

l_tot = n_new*(size/w); % m/m^2
pl = pl*(l/l_tot); % Wm/m^2

s = [r_a (1:n_new-1)*spacing_new];
d = [sqrt(s.^2+(2*depth)^2)];
num = [1 2*(1-(1:n_new-1)/n_new)];

T_min_i = min(T_s_30);
T_0_i = T_soil(1)+K;
pl_min_i = pl;

for j = 1:200
    f_i(j) = (mean(T_soil)-0)/(mean(T_soil)-T_min_i(j))+1e-2;
    pl_min_i = pl_min_i.*f_i(j);
    for i = 1:numel(t)
        p_th_i(i) = pl_min_i(i)/(1+(1/n_rel)*((T_h-T_0_i(i))/T_0_i(i)));
        e.setpower(t(i),p_th_i(i));
        for o = 1:numel(num)
            T_s(o) = num(o)*e.temperature(i*60*60*24*365.25/12,s(o));
            T_i(o) = num(o)*e.temperature(i*60*60*24*365.25/12,d(o));
        end
        T_s = sum(T_s);
        T_i = sum(T_i);
        T_0_i(i+1) = T_soil(i)-(T_s-T_i);
    end
    T_min_i(j+1) = min(T_0_i)-K;
    if T_min_i(j+1) <= 0
        T_0_i = T_0_i(2:end)-K;
        break
    end
end

f_new_i = prod(f_i);
size_new_i = size*1/f_new_i;

figure(6)
plot(t30/60/60/24/365.25,T_0_i)
grid on
xlabel('Time [a]')
ylabel('Temperature [°C]')

disp(['New size: ' num2str(size_new_i) ' times living area'])
disp(['Saved area: ' num2str(size-size_new_i) ' times living area'])

%% Assignment 11 a)
load Ppla.mat
lambda = 2.5;
rho = 2300;
c = 900;
kappa = lambda/(rho*c);
r = 0.14/2;
r_b = 0.1;
r_a = r*exp(-2*pi*lambda*r_b);
pl = Power_p_area/l_new;
pl = repmat(pl,30,25);
T = repmat(T_m,1,25,2);
T(2:361,:,:) = NaN; 
p_th =[];
T_min = repmat(-T_m,1,25);

[X,Y] = meshgrid(-50:25:50,-50:25:50);
x = X(1:end);
y = Y(1:end);
dist = zeros(25,25);
f = NaN(numel(t30),25);

for i = 1:numel(x)
    for j = 1:numel(y)
        if i~=j
            dist(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
        else
            dist(i,j) = r_a;
        end
    end
end

e = DownholeHeatExchanger(lambda,kappa,r);
DHE = repmat(e,1,25);
c = 1;

for z = 1:2
    while min(T_min) < K
        if z == 2
        f(c,:) = (T_m-K)./(T_m-T_min);
        pl = pl.*f(c,:);
        end
        for i = 1:numel(t30)
%             if z == 2
%                 f(i,j) = (T_m-K)/(T_m-T_min(j));
%                 pl(:,j) = pl(:,j).*f(i,j);
%             end
            %     for o = 1:numel(DHE)
            for j = 1:numel(DHE)
                if T_min(j) < K
%                     if z == 2
%                         f(i,j) = (T_m-K)/(T_m-T_min(j));
%                         pl(:,j) = pl(:,j).*f(i,j);
%                     end
                    p_th(i,j) = pl(i,j)/(1+(1/n_rel)*((T_h-T(i,j,z))/T(i,j,z)));
                    DHE(j).setpower(t30(i),p_th(i,j));
                    for o = 1:numel(DHE)
                        T_temp(o) = DHE(o).temperature(i*60*60*24*365.25/12,dist(j,o));
                    end
                    %         return
                    T_temp = sum(T_temp);
                    %         end
                    T(i+1,j,z) = T_m-T_temp;
                    %         return
                    
                else
                    continue
                end
                
            end
%             if z == 2 && i == 200
%                 return
%             end
            %                 return
        end
        T_min = min(T(:,:,z));
%         if z == 2 && i == 200
%             return
%         end

%                 return
        if z == 1
            break
        end
        %         return
        c = c+1;
    end
    %         return
end
%%
h = jet(25);
b = cell(25,1);
for i = 1:25
figure(4)
plot(t30/365.25/60/60/24,T(2:end,i,1)-K,'color',h(i,:))
grid on
hold on
b{i} = strcat('DHE-',num2str(i));
end
b{end+1} = 'T from Ass 8 h)';
plot(t30/365.25/60/60/24,Temp_min,'k--')
xlabel('Time [a]')
ylabel('Temperature [°C]')
legend(b,'Location','southoutside','NumColumns',5);
%% b)
for i= 1:25
    f_temp = f(:,i);
    f_temp = f_temp(~isnan(f_temp));
    f2(i) = prod(f_temp);
end
l_grid = l*1./f2;

figure(7)
plot(1:numel(l_grid),l_grid,'+')
grid on
hold on
pq = plot(numel(l_grid)+1,l_new,'+');
xlabel('Borehole')
ylabel('Length [m/m^2]')
legend(pq,'From Ass 8','Location','best')
xlim([0 27])

for i = 1:25
figure(8)
plot(t30/365.25/60/60/24,T(2:end,i,2)-K,'color',h(i,:))
grid on
hold on
end
plot(t30/365.25/60/60/24,Temp_min,'k--')
xlabel('Time [a]')
ylabel('Temperature [°C]')
legend(b,'Location','southoutside','NumColumns',5);
%%
figure(9)
contour3(X,Y,[l_grid(1:5);l_grid(6:10);l_grid(11:15);l_grid(16:20);l_grid(21:25)],100)
hold on
plot3(x,y,l_grid,'r*')
xlabel('x_1 [m]')
ylabel('x_2 [m]')
zlabel('Borehole length [m/m^2]')
