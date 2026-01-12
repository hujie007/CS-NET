
%%
filename1 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\BMP\cosembedding.csv";
filename2 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\BMP\coscenter.csv";
filename3 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\BMP\truegroup.csv";
BMP_cos_embedding = csvread(filename1,1,1);
BMP_cos_center = csvread(filename2,1,1);
BMP_true_group = csvread(filename3,1,1);
color_1 = [1 0 0];
color_2 = [0.12 0.56 1];
color_3 = [0.63 0.13 0.94];
cmap = [color_1; color_3; color_2];
INDEX_color = cmap(BMP_true_group,:);
x =  BMP_cos_embedding(:,1);y =  BMP_cos_embedding(:,2);z =  BMP_cos_embedding(:,3);
figure(1)
scatter3(x,y,z,30,INDEX_color)
view(-30,10)
% hold on
% scatter3(BMP_cos_center(:,1),BMP_cos_center(:,2),BMP_cos_center(:,3),50,'black','fill')
hold on
scatter3(0,0,0,30,'black','fill')
hold on
line([0,0.12*BMP_cos_center(1,1)],[0,0.12*BMP_cos_center(1,2)],[0,0.12*BMP_cos_center(1,3)],'Color','black','LineWidth',1.5,'LineStyle','--')
line([0,0.25*BMP_cos_center(2,1)],[0,0.25*BMP_cos_center(2,2)],[0,0.25*BMP_cos_center(2,3)],'Color','black','LineWidth',1.5,'LineStyle','--')
line([0,0.17*BMP_cos_center(3,1)],[0,0.17*BMP_cos_center(3,2)],[0,0.17*BMP_cos_center(3,3)],'Color','black','LineWidth',1.5,'LineStyle','--')
% set(gca,'xtick',[-1.5,-1,-0.5,0,0.5])
% set(gca,'ytick',[-1 -0.5 0 0.5 1 1.5])
% set(gca,'ztick',[-1,-0.5 0 0.5,1])


%%
filename4 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\BMP\sqembedding.csv";
filename5 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\BMP\sqcenter.csv";
BMP_sq_embedding = csvread(filename4,1,1);
BMP_sq_center = csvread(filename5,1,1);
a =  BMP_sq_embedding(:,1);b =  BMP_sq_embedding(:,2);c =  BMP_sq_embedding(:,3);
figure(2)
scatter3(a,b,c,30,INDEX_color)
view(-30,10)
hold on
scatter3(BMP_sq_center(:,1),BMP_sq_center(:,2),BMP_sq_center(:,3),40,'black','fill')
% hold on
% scatter3(0,0,0,30,'black','fill')
% hold on
% line([0,BMP_sq_center(1,1)],[0,BMP_sq_center(1,2)],[0,BMP_sq_center(1,3)],'Color','black','LineWidth',1.2,'LineStyle','--')
% line([0,BMP_sq_center(2,1)],[0,BMP_sq_center(2,2)],[0,BMP_sq_center(2,3)],'Color','black','LineWidth',1.2,'LineStyle','--')
% line([0,BMP_sq_center(3,1)],[0,BMP_sq_center(3,2)],[0,BMP_sq_center(3,3)],'Color','black','LineWidth',1.2,'LineStyle','--')
% % set(gca,'xtick',[-1.5 -1,-0.5,0,0.5])
% % set(gca,'ytick',[-1.5 -1 -0.5 0 0.5 1 1.5])
% % set(gca,'ztick',[-1,-0.5 0 0.5,1])
