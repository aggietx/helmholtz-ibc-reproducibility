function [A,F]=diri_set(A,F,diridof,diri_value)

if size(diridof,1)==1
    diridof=diridof.';
end

A(diridof,:)=sparse(1:size(diridof,1),diridof',ones(1,size(diridof,1)),size(diridof,1),size(A,2));
F(diridof,:)=diri_value;