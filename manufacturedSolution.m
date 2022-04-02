function a = manufacturedSolution(x,y)
a.u1    = cos((1/2).*x) + cos((1/4 ).*y );
a.du1dx = -(1/2).*sin((1/2).*x);
a.du1dy = -(1/4).*sin((1/4).*y);

a.U1 = sin((1/3).*x) + sin((1/6 ).*y );
a.dU1dx = (1/3).*cos((1/3).*x);
a.dU1dy = (1/6).*cos((1/6).*y);

a.u2 = cos((1/4).*x) + cos((1/8 ).*y );
a.du2dx = -(1/4).*sin((1/4).*x);
a.du2dy = -(1/8).*sin((1/8).*y);

a.U2 = sin((1/5).*x) + sin((1/10).*y );
a.dU2dx = (1/5).*cos((1/5).*x);
a.dU2dy = (1/10).*cos((1/10).*y);


a.du1U1dx = a.du1dx.*a.U1 + a.u1.*a.dU1dx;
a.du1U1dy = a.du1dy.*a.U1 + a.u1.*a.dU1dy;
a.du2U1dx = a.du2dx.*a.U1 + a.u2.*a.dU1dx;
a.du2U1dy = a.du2dy.*a.U1 + a.u2.*a.dU1dy;

a.du1U2dx = a.du1dx.*a.U2 + a.u1.*a.dU2dx;
a.du1U2dy = a.du1dy.*a.U2 + a.u1.*a.dU2dy;
a.du2U2dx = a.du2dx.*a.U2 + a.u2.*a.dU2dx;
a.du2U2dy = a.du2dy.*a.U2 + a.u2.*a.dU2dy;

a.divU1 = a.du1U1dx + a.du2U1dy;
a.divU2 = a.du1U2dx + a.du2U2dy;