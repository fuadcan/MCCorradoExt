// This example simulates animation by updating the
// viewer window with sequential 3d graphs.

// The VWR window will close automatically after
// animation has completed."

library pgraph;
graphset;

// use the -wp argument to position viewer window
graphprt("-wp=200,200,440,380");

// setup graphics variables
pqgwin one;
_pdate=0;
_pframe=0;
_paxes=0;
_psurf= { 0, 0 };
_pzclr={ -.4 15, -.32 10, -.2 2, 0. 3, .2 11, .28 15 };
let v = 100 100 640 480 0 0 1 6 15 0 0 2 2;
_plwidth = { 2 2 0 };
_pframe = 0;
_pcolor = { 13 12 9 };
ztics(0,1,.5,0);
t = seqa(0,.0157,400);
a=10; b=1; c=.3;
x = 1.5*(b*(1-c^2*cos(a*t)^2)^.5 .* cos(t));
y = 1.5*(b*(1-c^2*cos(a*t)^2)^.5 .* sin(t));
z = c*cos(a*t) + .5;
a=10; c=.5;
x1 = c*(1+cos(a*t)) .* cos(t) * 1.2;
y1 = c*(1+cos(a*t)) .* sin(t) * 1.2;
z1 = c*(1+cos(a*t)) * .5;
i=1;

// cylce through animation frames
for i (1,10,1);
    for j (1,19,1);
        view(j,3,5);
        xyz(x~x1,y~y1,z~z1);
    endfor;

    for j (19,1,-1);
        view(j,3,5);
        xyz(x~x1,y~y1,z~z1);
    endfor;
endfor;

// draw last frame and use the -w option to close window automatically after 2 seconds
graphprt("-w=2 -wp=200,200,440,380");
view(0,3,5);
xyz(x~x1,y~y1,z~z1);

