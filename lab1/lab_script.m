clear all;
close all;

%������� ��������� ������������� �� ������� � 
%���������� �������� �� ��������� ������������� [0,10] �.
Ts=0.1;
Ns=10000;

%������� ������������� ������� (���������� ����� �������)
nqm1=8;
nqm2=3;

%������� ���������� � ���������� ��������� �������� Li (a) � Lo (b)
nf=2;
minf=[0.2 0.2];
maxf=[1.0 1.0];

%������������ �������� �������������� ����� ������������
fracfact('a b ab');
N=2^nf;
fracplan=ans;
fictfact=ones(N,1);
X=[fictfact ans]';
fraceks=zeros(N,nf);

for i=1:nf
    for j=1:N
        fraceks(j,i)=minf(i)+(fracplan(j,i)+1)*(maxf(i)-minf(i))/2;
    end
end

%% ����������� ������������ � ����������� �������������
%������� �������������� ��������� � ������ ���������� ������
dm=0.3;
alpha=0.1;
%����������� t-������������
tkr_alpha=norminv(1-alpha/2);

%���� �� ������������ ������������� ��������������� �����
for j=1:N
    
    a=fraceks(j,1); 
    b=fraceks(j,2);
    Li=a; % lambda
    Lo=b; % mu
    
    %����������� ����� �������������� ���������
    NE=1;
    e=0; l=0;
    SQ=0; D=1;
    while NE < tkr_alpha^2*D/dm^2        
        %�������� ���������������� �������
        sim('chart_model_vOld',Ts*Ns);

        %�������� ����������� ������ ����������
        u = done(end) / total(end);
        v = time(end) / done(end);

        %���������� ����������� � ������ ���������� ��������� D ����������� ���������
        e=e+u;
        l=l+v;
        SQ=SQ+u^2;
        
        if NE>1
            D=SQ/(NE-1)-(e^2)/(NE*(NE-1));
        end
        NE=NE+1;
    end
    NE=NE-1
    %������ ����������� (�������) �� NE �����������
    qm=e/NE
    tm=l/NE
    Ya(j)=qm;
    Yb(j)=tm;
end


%% ����������� ������������� ��������� ��� qm � taym
C=X*X';
b_=inv(C)*X*Ya'
b_1=inv(C)*X*Yb'
%������������ ������������ ������� ������� �� ���������
%�������� ��������
A=minf(1):0.05:maxf(1);
B=minf(2):0.05:maxf(2);
[k N1]=size(A);
[k N2]=size(B);
for i=1:N1
    for j=1:N2
        an(i)=2*(A(i)-minf(1))/(maxf(1)-minf(1))-1;
        bn(j)=2*(B(j)-minf(2))/(maxf(2)-minf(2))-1;
        %����������������� ����������� �������
        Yca(j,i)=b_(1)+an(i)*b_(2)+bn(j)*b_(3)+an(i)*bn(j)*b_(4);
        Ycb(j,i)=b_1(1)+an(i)*b_1(2)+bn(j)*b_1(3)+an(i)*bn(j)*b_1(4);
    end
end

%����������� ������������ � ���������� ������� ��� ���������� qm 
[x,y]=meshgrid(A,B);
figure;
subplot(1,2,1),plot3(x,y,Yca),
xlabel('Li'), %'Li (Intense of entry)'),
ylabel('Lo'), %'Lo (Intense of service)'),
zlabel('(done/all)'),
title('Ratio done requests/all requests'),
grid on,
subplot(1,2,2),plot3(x,y,Ycb),
xlabel('Li'), %'Li (Intense of entry)'),
ylabel('Lo'), %'Lo (Intense of service)'),
zlabel('Average time'),
title('Average Request Time in System'),
grid on;
