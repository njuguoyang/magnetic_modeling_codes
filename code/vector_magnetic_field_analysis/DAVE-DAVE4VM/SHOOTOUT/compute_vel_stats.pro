function compute_vel_stats,XT,YT,X,Y,good=good,theta=theta

  dF1=moment(sqrt((((X-XT)^2+(Y-YT)^2)/(XT^2+YT^2))[good]),/double)
  df2=moment(((sqrt(X^2+Y^2)-sqrt(XT^2+YT^2))/sqrt(XT^2+YT^2))[good],/double)
  CVEC=total((X*XT+Y*YT)[good],/double)/sqrt(total((X^2+Y^2)[good],/double)*total((XT^2+YT^2)[good],/double))
  CS=moment(((X*XT+Y*YT)/sqrt((X^2+Y^2)*(XT^2+YT^2)))[good],/double)
  theta=180*atan(X*YT-Y*XT,X*XT+Y*YT)/!dpi


return,{fractional_error:DF1,fractional_magnitude:DF2,CVEC:CVEC,CS:CS,theta:moment(theta[good],/double)}

end

   
