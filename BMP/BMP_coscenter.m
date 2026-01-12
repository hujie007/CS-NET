filename = "C:\Users\lenovo\Desktop\cosembedding\Numerical Experiments\BMP\coscenter.csv";
BMP_cos_center = csvread(filename,1,1);
coscenter=[];
for i = 1:3
    for j = i+1:3
        coscenter(i,j) = acos(dot(BMP_cos_center(i,:),BMP_cos_center(j,:)))*180/pi;
    end
end