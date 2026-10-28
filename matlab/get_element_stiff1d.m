function localstiff1d=get_element_stiff1d(p,hx)

n_node=1+p;
[x,wx]=get_gll_weight_position(p);
deri_x=derlgl(x,n_node);%first derivative GLL matrix % % saved in each column 
jacx=2/hx; 
localstiff1d=deri_x'*(diag(wx)*deri_x)*jacx;