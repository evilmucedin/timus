CONST
     EPS = 1E-9;

TYPE
    TARRAY = array[1..100] of extended;

VAR
   n           : integer;
   x, y, a     : TARRAY;
   acos, asin  : TARRAY;
   xv, yv      : extended;
   xvs, yvs    : extended;
   i           : integer;

PROCEDURE INPUTDATA;
          var
             i : integer;
          begin
               Assign(input, '');
               ReSet(input);
               Read(n);
               for i := 1 to n do
                   Read(x[i], y[i]);
               for i := 1 to n do
                   Read(a[i]);
               for i := 1 to n do
                   a[i] := a[i] / 180.0 * Pi;
               for i := 1 to n do
                   acos[i] := cos(a[i]);
               for i := 1 to n do
                   asin[i] := sin(a[i]);
               Close(input);
          end;

PROCEDURE ROTATION(x0, y0 : extended; alphai : integer; x1, y1 : extended; var xx1, yy1 : extended);
         var
            dx : extended;
            dy : extended;
         begin
              dx := x1 - x0;
              dy := y1 - y0;
              xx1 := x0 + dx*acos[alphai] - dy*asin[alphai];
              yy1 := y0 + dx*asin[alphai] + dy*acos[alphai];
         end;

PROCEDURE OUTPUTFLOAT(x : extended);
          var
             xi : longint;
          begin
              xi := Round(x*100);
              Write(xi div 100);
              if xi mod 100 <> 0 then
                 begin
                 Write('.');
                 if xi mod 10 <> 0 then
                    Write(xi mod 100)
                 else
                     Write(xi mod 100 div 10);
                 end;
          end;

BEGIN
     INPUTDATA;
     xv := 0.0;
     yv := 0.0;
     xvs := 10.0;
     yvs := 10.0;
     while Abs(xvs - xv) + Abs(yvs - yv) > EPS do
           begin
           xv := (xvs + xv) / 2.0;
           yv := (yvs + yv) / 2.0;
           xvs := xv;
           yvs := yv;
           for i := 1 to n do
               ROTATION(x[i], y[i], i, xvs, yvs, xvs, yvs);
           end;

     for i := 1 to n do
         begin
         OUTPUTFLOAT(xvs);
         Write(' ');
         OUTPUTFLOAT(yvs);
         WriteLN;
         ROTATION(x[i], y[i], i, xvs, yvs, xvs, yvs);
         end;
END.
