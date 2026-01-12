filename = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\USblog\coscenter.csv";
USblog_cos_center = csvread(filename,1,1);
cos_center = acos(dot(USblog_cos_center(1,:),USblog_cos_center(2,:)))*180/pi;
 

