function umsfine=smooth(Fd,Astiffd,umsfine,d)


    for ii=1:2
   umsfine = umsfine + 0.75*d.*(Fd - Astiffd*umsfine); 
    end

%    e=sqrt( diff'*Astiffd*diff )/sqrt(u'*Astiffd*u)
