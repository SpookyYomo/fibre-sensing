#IQ-sampling
a=load("trial.dat"); 

#2d rotation operated on I&Q to account for frequency mismatch (no..)
theta=0;
R=[cos(theta) -sin(theta);
   sin(theta)  cos(theta)]


b=reshape(a(:,2), 4, length(a)/4);

c=transpose(R*b(1:2,:));
#b=transpose(b);

d=atan2(c(:,1), c(:,2));
e=zeros(length(d),1);
e(1)=d(1);

#dealing with points jumping between +- pi
for i=2:length(d)
  if( (d(i)-d(i-1)) > 5 )
    d(i:end)=d(i:end).-2*3.14159;
    #d(i)=e(i);
  elseif( (d(i)-d(i-1)) < -5 )
    d(i:end)=d(i:end).+2*3.14159;
    #d(i)=e(i);
  endif
endfor
#applying rotation, theta should be small
plot(d)
dlmwrite("out2.dat", d, "delimiter", " ");