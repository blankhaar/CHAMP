subroutine int_gam_w(gam,wa,wb,w0,lF1,lF2,dw,a1,b1,n,m,g_r,g_i,g_w_r,g_w_i)
implicit none 
integer n,m,nm,lF1,lF2,a1,b1
double precision x0,x1,x2,xm4,wab,dv,dvc,c0,dwm,w0
double precision dE,dw,gam,g_r,g_i,g_w_r,g_w_i
double precision fac1,fac2,fac,q_up,q_down
double precision, dimension (lF1) :: wa
double precision, dimension (lF2) :: wb

c0 = 299792458.d0

x1 = 1.d0
x2 = 2.d0

wab = w0 + wa(a1) - wb(b1)

dv = c0*dw/w0;
dvc = dv*n/c0

dwm = dw*m

q_up = wab - (dw/x2 + w0 + dwm)*(x1-dvc);
q_down = wab - (-dw/x2 + w0 + dwm)*(x1-dvc);

g_r = -x1*(datan(q_up/gam) - datan(q_down/gam)) ;

fac1 = gam**x2 + q_up**x2
fac2 = gam**x2 + q_down**x2

g_i = dlog(fac1/fac2)/x2;

g_w_r = (gam*g_i/(x1-dvc) + (wab/(x1-dvc) - (w0+dwm))*g_r);
g_w_i = (wab/(x1-dvc) - (w0+dwm))*g_i + gam/(x1-dvc)*g_r + dw
!g_w_i = 0.d0

end subroutine

subroutine int_gam_v(gam,wa,wb,w0,lF1,lF2,dv,a1,b1,n,m,g_r,g_i,g_w_r,g_w_i)
implicit none 
integer n,m,nm,lF1,lF2,a1,b1
double precision x0,x1,x2,xm4,wab,dv,dvc,c0,dwm,w0
double precision dE,dw,gam,g_r,g_i,g_w_r,g_w_i
double precision fac1,fac2,fac,q_up,q_down
double precision, dimension (lF1) :: wa
double precision, dimension (lF2) :: wb

c0 = 299792458.d0

x1 = 1.d0
x2 = 2.d0

wab = w0 + wa(a1) - wb(b1)

dw = w0*dv/c0

dvc = dv*n/c0

dwm = dw*m

q_up = wab - (dw/x2 + w0 + dwm)*(x1-dvc);
q_down = wab - (-dw/x2 + w0 + dwm)*(x1-dvc);

g_r = -x1*(datan(q_up/gam) - datan(q_down/gam)) ;

fac1 = gam**x2 + q_up**x2
fac2 = gam**x2 + q_down**x2

g_i = dlog(fac1/fac2)/x2;

g_w_r = (gam*g_i/(x1-dvc) + (wab/(x1-dvc) - (w0+dwm))*g_r)
g_w_i = (wab/(x1-dvc) - (w0+dwm))*g_i + gam/(x1-dvc)*g_r + dw

fac = c0/w0
g_r = fac*g_r
g_i = fac*g_i
g_w_r = -fac*fac*g_w_r
g_w_i = -fac*fac*g_w_i  
!g_w_i = 0.d0 

end subroutine

