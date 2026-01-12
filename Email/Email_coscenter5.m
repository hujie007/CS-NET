filename5 = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\Email\coscenter5.csv";
Email_cos_center5 = csvread(filename5,1,1);
coscenter5=[];
for i = 1:5
    for j = i+1:5
        coscenter5(i,j) = acos(dot(Email_cos_center5(i,:),Email_cos_center5(j,:)))*180/pi;
    end
end
