network0= importdata("C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\email-Eu-core.txt");
label0= importdata("C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\email-Eu-core-department-labels.txt");
label0=label0+1;
network0=network0+1;
n=size(label0,1);A=[];  
for k=1 : size(network0,1)
    edgei= network0(k,1);
    edgej= network0(k,2);
    A(edgei,edgej)= 1;
    A(edgej,edgei)= 1;
end
A=A-diag(diag(A));
exnodes= [find(label0(:,2)==9); find(label0(:,2)==10); find(label0(:,2)==12); find(label0(:,2)==22)];
%exnodes= [find(label0(:,2)==1); find(label0(:,2)==7); find(label0(:,2)==8); find(label0(:,2)==17);find(label0(:,2)==18)];
%exnodes= [find(label0(:,2)==6); find(label0(:,2)==14); find(label0(:,2)==16); find(label0(:,2)==20);find(label0(:,2)==21);find(label0(:,2)==23)];
exlabel= label0(exnodes,:);
a= A(exnodes,exnodes);

nozero= sum(a,2)~=0;
exnodes=exnodes(nozero);
a= a(nozero,nozero);
exlabel= exlabel(nozero);
n= size(a,1);
exgraph= graph(a);
plot(exgraph)
 %%
filename_g = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\truegroup.csv";
filename_cos = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\cosgroup.csv";
filename_sq = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\sqgroup.csv";
Email_truegroup= csvread(filename_g,1,1);
Email_cosgroup= csvread(filename_cos,1,1);
Email_sqgroup= csvread(filename_sq,1,1);



figure(1)
a1=repmat([0.96 0.64 0.38],sum(Email_truegroup==1),1);% É³×ØÉ«
a2=repmat([0.53 0.81 0.92],sum(Email_truegroup==2),1); % ÌìÀ¶
a3=repmat([0.74 0.99 0.79],sum(Email_truegroup==3),1); % ±¡ºÉÉ«
a4=repmat([0.87 0.63 0.87],sum(Email_truegroup==4),1); % Ã·ºìÉ«
%a5=repmat([0.98 0.5 0.45],sum(Email_truegroup==5),1); %·ÛºìÉ«
%a6=repmat([1 0.84 0],sum(Email_truegroup==6),1); %»ÆÂÌÉ«
p1=plot(exgraph)
p1.Marker = 'o'
p1.MarkerSize = 3.5;
p1.LineWidth = 0.1;
p1.EdgeColor=[0.5 0.54 0.53];
p1.NodeCData=1:n;
p1.NodeColor = [a1;a2;a3;a4];
%p1.NodeColor = [a1;a2;a3;a4;a5];
%p1.NodeColor = [a1;a2;a3;a4;a5;a6];
%colormap([c1;c2;c3;c4])

figure(2)
b1=repmat([0.96 0.64 0.38],sum(Email_cosgroup==4),1);
b2=repmat([0.53 0.81 0.92],sum(Email_cosgroup==3),1);
b3=repmat([0.74 0.99 0.79],sum(Email_cosgroup==1),1);
b4=repmat([0.87 0.63 0.87],sum(Email_cosgroup==2),1);
%b5=repmat([0.98 0.5 0.45],sum(Email_cosgroup==3),1); %·ÛºìÉ«
%b6=repmat([1 0.84 0],sum(Email_cosgroup==5),1); %½ð»ÆÉ«
p2=plot(exgraph)
p2.Marker = 'o'
p2.MarkerSize = 3.5;
p2.LineWidth = 0.1;
p2.EdgeColor=[0.5 0.54 0.53];
p2.NodeCData=1:n;
p2.NodeColor = [b1;b2;b3;b4]
%p2.NodeColor = [b1;b2;b3;b4;b5];
%p2.NodeColor = [b1;b2;b3;b4;b5;b6];
% 
% sum(Email_cosgroup(Email_truegroup==1)~=3)
% sum(Email_cosgroup(Email_truegroup==2)~=4)
% sum(Email_cosgroup(Email_truegroup==3)~=2)
% sum(Email_cosgroup(Email_truegroup==4)~=5)
% sum(Email_cosgroup(Email_truegroup==5)~=1)
% sum(Email_cosgroup(Email_truegroup==6)~=6)



figure(3)
c1=repmat([0.96 0.64 0.38],sum(Email_sqgroup==4),1);
c2=repmat([0.53 0.81 0.92],sum(Email_sqgroup==2),1);
c3=repmat([0.74 0.99 0.79],sum(Email_sqgroup==1),1);
c4=repmat([0.87 0.63 0.87],sum(Email_sqgroup==3),1);
%c5=repmat([0.98 0.5 0.45],sum(Email_sqgroup==5),1); %·ÛºìÉ«
%c6=repmat([1 0.84 0],sum(Email_sqgroup==3),1); %»ÆÂÌÉ«
p3=plot(exgraph)
p3.Marker = 'o'
p3.MarkerSize = 3.5;
p3.LineWidth = 0.1;
p3.EdgeColor=[0.5 0.54 0.53];
p3.NodeCData=1:n;
p3.NodeColor = [c1;c2;c3;c4];
%p3.NodeColor = [c1;c2;c3;c4;c5];
%p3.NodeColor = [c1;c2;c3;c4;c5;c6];

% sum(Email_sqgroup(Email_truegroup==1)~=4)
% sum(Email_sqgroup(Email_truegroup==2)~=2)
% sum(Email_sqgroup(Email_truegroup==3)~=1)
% sum(Email_sqgroup(Email_truegroup==4)~=5)
% sum(Email_sqgroup(Email_truegroup==5)~=3)
% sum(Email_sqgroup(Email_truegroup==6)~=6)