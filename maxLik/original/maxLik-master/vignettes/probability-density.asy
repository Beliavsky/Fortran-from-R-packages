unitsize(25mm,65mm);
defaultpen(fontsize(9));

real xLeft = -2.2;
real xRight = 2.2;
real yTop = 0.5;

// normal density
real dnorm(real x) {
  return 1/sqrt(2*pi)*exp(-1/2*x^2);
}
// compute normal curve, plot later
path normalCurve;
for(real x = xLeft + 0.1; x < xRight - 0.1; x += 0.15) {
  normalCurve = normalCurve..(x, dnorm(x));
}

// Example points
real xs[] = {-1.695, 0.3};
real delta = 0.15;
int i = 1;
for(real x : xs) {
  real fx = dnorm(x);
  real xl = x - delta/2;
  real xr = x + delta/2;
  real tl = times(normalCurve, xl)[0];
  real tr = times(normalCurve, xr)[0];
  path striptop = subpath(normalCurve, tl, tr);
  path area = (xl, 0)--striptop--(xr, 0)--cycle;
  filldraw(area, lightgray, linewidth(0.2));
  draw((x, 0)--(x, dnorm(x)), dashed);
  label("$x_" + string(i) + " = " + format("%f", x) +"$", (x, 0), S + 0.2E);
  // width marks and width
  real barheight = dnorm(x) + 0.06;
  Label widthLabel = Label("width $\delta$", MidPoint, 2N);
  draw(widthLabel, (xl, barheight)--(xr, barheight), linewidth(0.4), Bars);
  arrow((xl, barheight), W, length=50delta, margin=DotMargin, linewidth(0.4));
  arrow((xr, barheight), E, length=50delta, margin=DotMargin);
  // mark the function value
  real xmarker = x + 1.5delta;
  draw((x, fx)--(xmarker,fx), dotted);
  Label valueLabel = Label("$f(x_" + string(i) + ") = " + format("%5.3f", fx) + "$",
			   position=EndPoint,  E);
  path valuePath = (xmarker, fx)--(xmarker+delta, fx);
  draw(valueLabel, valuePath, linewidth(0.4));
  pair barx = relpoint(valuePath, 0.5);
  draw((barx.x, 0)--barx, Arrow(4));
  //
  ++i;
}

// add normal curve later as filling area cuts into the curve otherwise
draw(normalCurve, linewidth(0.7));

// Add Axes after are to avoid cutting into it
path xaxis = (xLeft,0)--(xRight,0);
path yaxis = (0,0)--(0,yTop);
draw(xaxis, Arrow(TeXHead, 1));
draw(yaxis, Arrow(TeXHead, 1));
label("$x$", point(xaxis, 1), 2S);
// Axis labels
real tickLength = 0.05*yTop;
for(int x = (int)xLeft; x <= (int)xRight; ++x) {
  draw((x,0)--(x,-tickLength));
  label(string(x), (x,-tickLength), 3S);
}
