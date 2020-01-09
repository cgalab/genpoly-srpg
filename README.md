The code srpg generates simply-connected and multiply-connected polygons by
means of a regular grid (with quadratic cells).

# Compilation

Any standard C compiler will do. E.g.:

    gcc -O2 -Wall -o srpg srpg.c

This is also basically what the Makefile does, so you could also just run `make`.

Usage:

    srgp --Nx <X> --Ny <Y> --output <OUTPUTFILE>
         [--percent <P>] [--seed <S>] [--holes]
         [--aligned | --perturb [--smooth <M>]]
         [--hierarchy <N>] [--diagonal]

where X,Y,M,N are positive integers, S is a non-negative integer,
0.001 < P < 0.5 is a real, and OUTPUTFILE is the name of the output file.

# Generating polygonal data

The computation is based on a grid with X times Y quadratic cells. The user
has to specify X and Y. In addition, the user has to specify the name of an
output file. If no other options are specified, then srpg generates a polygon
whose edges are parallel to the coordinate axes. The number of vertices of the
polygon generated is random, but it does depend on X, Y and the percentage P:
The larger X and Y, the more vertices the polygon tends to have if P is kept
constant. Visually pleasing "random" polygons can be achieved by selecting
fairly small values for P, e.g., P:=0.1 or even P:=0.01. (However, a small
value of P will also reduce the number of vertices of the polygon.)

The option "--aligned" will cause all vertices to lie on grid points, i.e., to
have integer coordinates. The option "--diagonal" causes srpg to cut off some
corners by line segments with inclination +/-1, thus generating an octagonal
polygon. If the option "--perturb" is used then the vertices are moved away
from the grid points and (most) polygon edges will not be parallel to the
coordinate axes. The option "--hierarchy N" instructs srpg to apply N rounds
of a recursive refinement to the polygon generated. (Typically, N will be a
small positive integer.) If the option "--holes" is specified then srpg will
generate a multiply-connected polygonal area. The option "--smooth M" tells
srpg to apply M rounds of corner cutting to the polygon generated, thus
generating a polygon which resembles a polygonal approximation of a smooth
free-form curve. (Again, M will be a small positive integer.)

Please direct bug reports or suggestions to Martin Held at held@cs.sbg.ac.at.
