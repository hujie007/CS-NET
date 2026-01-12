filename4 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\coscenter4.csv";
Email_cos_center4 = csvread(filename4,1,1);
coscenter4=zeros(4);
for i = 1:4
    for j = i+1:4
        coscenter4(i,j) = acos(dot(Email_cos_center4(i,:),Email_cos_center4(j,:)))*180/pi;
    end
end
coscenter4 = coscenter4'+coscenter4;
% 热图函数为heatmap;开始绘制
% 创建热图
figure(1);
hotfig = heatmap([9,10,12,22],[9,10,12,22],coscenter4,'FontSize',12);
hotfig.Colormap = jet(500); %matlab默认热图配色 help graph3d查找挑选
% hsv, hot, gray, bone, copper, pink, white, flag, lines, colorcube, vga,
% jet,prism,cool,autumn,spring,winter,summer
%hot_figure.GridVisible = 'on';
% 设置坐标区名字与图的大标题
hotfig.CellLabelFormat = '%0.2f';
%hotfig.Title = 'The angle of community center';
hotfig.XLabel = 'community label';
hotfig.YLabel = 'community label';

filename5 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\coscenter5.csv";
Email_cos_center5 = csvread(filename5,1,1);
coscenter5=zeros(5);
for i = 1:5
    for j = i+1:5
        coscenter5(i,j) = acos(dot(Email_cos_center5(i,:),Email_cos_center5(j,:)))*180/pi;
    end
end
coscenter5 = coscenter5'+coscenter5;
figure(2);
hotfig = heatmap([1,7,8,17,18],[[1,7,8,17,18]],coscenter5,'FontSize',12);
hotfig.Colormap = jet(500); %matlab默认热图配色 help graph3d查找挑选
% hsv, hot, gray, bone, copper, pink, white, flag, lines, colorcube, vga,
% jet,prism,cool,autumn,spring,winter,summer
%hot_figure.GridVisible = 'on';
% 设置坐标区名字与图的大标题
hotfig.CellLabelFormat = '%0.2f';
%hotfig.Title = 'The angle of community center';
hotfig.XLabel = 'community label';
hotfig.YLabel = 'community label';

filename6 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\coscenter6.csv";
Email_cos_center6 = csvread(filename6,1,1);
coscenter6=zeros(6);
for i = 1:6
    for j = i+1:6
        coscenter6(i,j) = acos(dot(Email_cos_center6(i,:),Email_cos_center6(j,:)))*180/pi;
    end
end
coscenter6 = coscenter6'+coscenter6;
figure(3);
hotfig = heatmap([6,14,16,20,21,23],[6,14,16,20,21,23],coscenter6,'FontSize',12);
hotfig.Colormap = jet(500); %matlab默认热图配色 help graph3d查找挑选
% hsv, hot, gray, bone, copper, pink, white, flag, lines, colorcube, vga,
% jet,prism,cool,autumn,spring,winter,summer
%hot_figure.GridVisible = 'on';
% 设置坐标区名字与图的大标题
hotfig.CellLabelFormat = '%0.2f';
%hotfig.Title = 'The angle of community center';
hotfig.XLabel = 'community label';
hotfig.YLabel = 'community label';


