/*****************************************************************************/
/*                                                                           */
/*            Copyright (C)          2008-2020            M. Held            */
/*                                                                           */
/*****************************************************************************/
/*                                                                           */
/* Written by:  Martin Held                                                  */
/*                                                                           */
/* E-Mail:      held@cs.sbg.ac.at                                            */
/* Snail Mail:  Martin Held                                                  */
/*              FB Computerwissenschaften                                    */
/*              Universitaet Salzburg                                        */
/*              A-5020 Salzburg, Austria                                     */
/*                                                                           */
/*****************************************************************************/

/* get standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#ifdef BOOL_DEFINED
typedef bool boolean;
#else
#define false 0
#define true  (!false)
typedef unsigned char  boolean;
#endif

#define NIL -1

typedef struct {
   double x;         /* x-coordinate                                         */
   double y;         /* y-coordinate                                         */
} coord;             /**********  point/vector  ******************************/

int dummy_random;

#define UniformRandom(C, M) ((C) = lrand48() % (M))

#define Convert(i, j, k) \
{ \
   k = i * N_y + j; \
} 

#define Invert(k, i, j) \
 { \
   i = k / N_y; \
   j = k - i * N_y; \
 }


#define Det2D(u, v, w) \
 (((u).i1 - (v).i1) * ((v).j1 - (w).j1) + ((v).j1 - (u).j1) * ((v).i1 - (w).i1))

#define Swap(i1, i2, i) \
{(i)  = (i1);           \
 (i1) = (i2);         \
 (i2) = (i); }


#define SetTop(I, J, W) \
{\
   assert((I >= 0)  &&  (I < N_x)  &&  (J >= 0)  &&  (J < N_y)); \
   top[I][J] = W; \
}


#define SetBot(I, J, W) \
{\
   assert((I >= 0)  &&  (I < N_x)  &&  (J >= 0)  &&  (J < N_y)); \
   bot[I][J] = W; \
}

#define SetLft(I, J, W) \
{\
   assert((I >= 0)  &&  (I < N_x)  &&  (J >= 0)  &&  (J < N_y)); \
   lft[I][J] = W; \
}

#define SetRgt(I, J, W) \
{\
   assert((I >= 0)  &&  (I < N_x)  &&  (J >= 0)  &&  (J < N_y)); \
   rgt[I][J] = W; \
}

#define IsTopFull(I, J)    top[I][J]
#define IsBotFull(I, J)    bot[I][J]
#define IsLftFull(I, J)    lft[I][J]
#define IsRgtFull(I, J)    rgt[I][J]


#define GetStartI(I, J) \
(\
   assert((I >= 0)  &&  (I <= N_x)  &&  (J >= 0)  &&  (J <= N_y)), \
   edges[I][J].i1)

#define GetStartJ(I, J) \
(\
   assert((I >= 0)  &&  (I <= N_x)  &&  (J >= 0)  &&  (J <= N_y)), \
   edges[I][J].j1)

#define GetEndI(I, J) \
(\
   assert((I >= 0)  &&  (I <= N_x)  &&  (J >= 0)  &&  (J <= N_y)), \
   edges[I][J].i2)

#define GetEndJ(I, J) \
(\
   assert((I >= 0)  &&  (I <= N_x)  &&  (J >= 0)  &&  (J <= N_y)), \
   edges[I][J].j2)

#define GetVis(I, J) \
(\
   assert((I >= 0)  &&  (I <= N_x)  &&  (J >= 0)  &&  (J <= N_y)), \
   edges[I][J].vis)

#define SetVis(I, J, B) \
{\
   assert((I >= 0)  &&  (I <= N_x)  &&  (J >= 0)  &&  (J <= N_y)); \
   edges[I][J].vis = B; \
}


boolean holes = false;
boolean **full = (boolean**)NULL;
boolean **old_full = (boolean**)NULL;
boolean **keep = (boolean**)NULL;
int N_x = 0, N_y = 0, N = 0;

boolean **top = (boolean**)NULL;
boolean **bot = (boolean**)NULL;
boolean **lft = (boolean**)NULL;
boolean **rgt = (boolean**)NULL;

int *candidates = (int*)NULL;
int num_candidates = 0, max_num_candidates = 0;
int *keep_candidates = (int*)NULL;
int num_keep_candidates = 0, max_num_keep_candidates = 0;

typedef struct {
   int i1;
   int j1;
   int i2;
   int j2;
   boolean vis;
} edge_node;

edge_node **edges = (edge_node**)NULL;

typedef struct {
   int i1;
   int j1;
} vertex_node;

vertex_node *vertices = (vertex_node*)NULL;
int num_vertices = 0;
int max_num_vertices = 0;



boolean** MakeBMatrix(int Nx, int Ny, boolean value)
{
   int i, j;
   boolean **matrix = (boolean**) malloc(Nx * sizeof(boolean *));

   if (matrix == (boolean**)NULL) {
      fprintf(stderr, 
              "MakeBMatrix() - Cannot get %d elements of %ld bytes...\n", 
              Nx, sizeof(boolean *));
      exit(1);
   }

   for (i = 0;  i < Nx;  ++i) {
      matrix[i] = (boolean *) malloc(Ny * sizeof(boolean));
      if (matrix[i] == (boolean*)NULL) {
         fprintf(stderr, 
                 "MakeBMatrix() - Cannot get %d elements of %ld bytes...\n", 
                 Ny, sizeof(boolean *));
         exit(1);
      }
   }

   for (i = 0;  i < Nx;  ++i) {
      for (j = 0;  j < Ny;  ++j) {
         matrix[i][j] = value;
      }
   }

   return matrix;
}



boolean** FreeBMatrix(boolean **matrix, int Nx)
{
   int i;

   for (i = 0;  i < Nx;  ++i) {
      free(matrix[i]);
   }
   free(matrix);

   return (boolean**)NULL;
}

   

edge_node** MakeEMatrix(int Nx, int Ny)
{
   int i, j;
   edge_node **matrix = (edge_node**) malloc(Nx * sizeof(edge_node *));

   if (matrix == (edge_node**)NULL) {
      fprintf(stderr, 
              "MakeEMatrix() - Cannot get %d elements of %ld bytes...\n", 
              Nx, sizeof(edge_node *));
      exit(1);
   }

   for (i = 0;  i < Nx;  ++i) {
      matrix[i] = (edge_node*) malloc(Ny * sizeof(edge_node));
      if (matrix[i] == (edge_node*)NULL) {
         fprintf(stderr, 
                 "MakeEMatrix() - Cannot get %d elements of %ld bytes...\n", 
                 Ny, sizeof(edge_node *));
         exit(1);
      }
   }

   for (i = 0;  i < Nx;  ++i) {
      for (j = 0;  j < Ny;  ++j) {
         matrix[i][j].i1  = NIL;
         matrix[i][j].j1  = NIL;
         matrix[i][j].i2  = NIL;
         matrix[i][j].j2  = NIL;
         matrix[i][j].vis = true;
      }
   }

   return matrix;
}



edge_node** FreeEMatrix(edge_node **matrix, int Nx)
{
   int i;

   for (i = 0;  i < Nx;  ++i) {
      free(matrix[i]);
   }
   free(matrix);

   return (edge_node**)NULL;
}

   

double Perturbation(void)
{
   int c;

   UniformRandom(c, 800001);
   c -= 400000;
   
   return (((double) c) / 899000.0);
}


int e_comp(const void *e1, const void *e2)
{
   if         (((edge_node*) e1)->i1 < ((edge_node*) e2)->i1)   return -1;
   else if    (((edge_node*) e1)->i1 > ((edge_node*) e2)->i1)   return  1;
   else  {
      if      (((edge_node*) e1)->j1 < ((edge_node*) e2)->j1)   return -1;
      else if (((edge_node*) e1)->j1 > ((edge_node*) e2)->j1)   return  1;
      else                                                      return  0;
   }
}


void *ReallocateArray(void *old_ptr, int number, size_t size)
{
   void *new_ptr;

   if (old_ptr != (void*)NULL) {
      new_ptr = (void*) realloc((void*) old_ptr, number * size);
      if (new_ptr == (void*)NULL)  {
         fprintf(stderr, 
                 "ReallocateArray() - Cannot get %d elements of %ld bytes...\n", 
                 number, size);
         exit(1);
      }
   }
   else {
      new_ptr = (void*) calloc(number, size);
      if (new_ptr == (void*)NULL) {
         fprintf(stderr, 
                 "ReallocateArray() - Cannot get %d elements of %ld bytes...\n", 
                 number, size);
         exit(1);
      }
   }

   return new_ptr;
}




#define StorePnt(P) \
{\
   if (num_pnts >= max_num_pnts) { \
      max_num_pnts += 131072; \
      pnts = (coord*) ReallocateArray(pnts, max_num_pnts, sizeof(coord)); \
   } \
   pnts[num_pnts] = P; \
   \
   ++(num_pnts); \
}


void StoreEdge(int i1, int j1, int i2, int j2)
{
   if ((i1 >= 0)  &&  (j1 >= 0)  &&  (i1 <= N_x)  &&  (j1 <= N_y)  &&
       (i2 >= 0)  &&  (j2 >= 0)  &&  (i2 <= N_x)  &&  (j2 <= N_y)) {
      SetVis(i1, j1, false);
      if (edges[i1][j1].i1 == NIL) {
         assert(edges[i1][j1].j1 == NIL);
         edges[i1][j1].i1 = i2;
         edges[i1][j1].j1 = j2;
      }
      else {
         assert(edges[i1][j1].j1 != NIL);
         assert(edges[i1][j1].i2 == NIL);
         assert(edges[i1][j1].j2 == NIL);
         edges[i1][j1].i2 = i2;
         edges[i1][j1].j2 = j2;
      }
      SetVis(i2, j2, false);
      if (edges[i2][j2].i1 == NIL) {
         assert(edges[i2][j2].j1 == NIL);
         edges[i2][j2].i1 = i1;
         edges[i2][j2].j1 = j1;
      }
      else {
         assert(edges[i2][j2].j1 != NIL);
         assert(edges[i2][j2].i2 == NIL);
         assert(edges[i2][j2].j2 == NIL);
         edges[i2][j2].i2 = i1;
         edges[i2][j2].j2 = j1;
      }
   }
   else {
      printf("StoreEdge(): index out of bounds!\n");
      printf("             (%d,%d) <--> (%d,%d)\n", i1, j1, i2, j2);
      exit(1);
   }

   return;
}


void StoreVertex(int i1, int j1)
{
   if (num_vertices >= max_num_vertices) {
      max_num_vertices += 131072;
      vertices = (vertex_node*) ReallocateArray(vertices, max_num_vertices, 
                                                sizeof(vertex_node));
   }
   vertices[num_vertices].i1  = i1;
   vertices[num_vertices].j1  = j1;
   ++num_vertices;

   return;
}



void OpenFile(FILE  **output, const char *file_name)
{
   /* open file */
   if((*output = fopen(file_name, "w")) == NULL) {
      fprintf(stderr,"*** Output file not created! ***\n");
      exit(1);
   }
   
   return;
}   



boolean IsFull(int i, int j)
{
   if ((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  &&  (j < N_y)) {
      return full[i][j];
   }
   else {
      return false;
   }
}



boolean IsOldFull(int i, int j, int N_x_old, int N_y_old, boolean **old_full)
{
   if ((i >= 0)  &&  (i < N_x_old)  &&  (j >= 0)  &&  (j < N_y_old)) {
      return old_full[i][j];
   }
   else {
      return false;
   }
}





boolean IsCompletelyFull(int i, int j)
{
   if ((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  &&  (j < N_y)) {
      return (top[i][j] && bot[i][j] && lft[i][j] && rgt[i][j]);
   }
   else {
      return false;
   }
}


boolean IsCompletelyEmpty(int i, int j)
{
   if ((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  &&  (j < N_y)) {
      return (!(top[i][j] || bot[i][j] || lft[i][j] || rgt[i][j]));
   }
   else {
      return true;
   }
}



boolean ToBeKept(int i, int j)
{

   if ((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  &&  (j < N_y)) {
      return keep[i][j];
   }
   else {
      return true;
   }
}



boolean IsPossible(int i, int j)
{
   int c = 0;

   if ((i < 0)  ||  (i >= N_x)  ||  (j < 0)  ||  (j >= N_y))  
      return false;

   if (ToBeKept(i, j)) return false;
   if (IsFull(i, j))   return false;

   if (IsFull(i-1,j+1)) {
      if (!(IsFull(i-1,j)  ||  IsFull(i,j+1)))  return false;
   }
   if (IsFull(i-1,j-1)) {
      if (!(IsFull(i-1,j)  ||  IsFull(i,j-1)))  return false;
   }
   if (IsFull(i+1,j-1)) {
      if (!(IsFull(i+1,j)  ||  IsFull(i,j-1)))  return false;
   }
   if (IsFull(i+1,j+1)) {
      if (!(IsFull(i+1,j)  ||  IsFull(i,j+1)))  return false;
   }

   if (holes)  UniformRandom(c, 30);

   if (c != 0) {
      if (IsFull(i,j-1)  &&  IsFull(i,j+1)) {
         if (!(IsFull(i+1,j)  ||  IsFull(i-1,j)))  return false;
      }
      if (IsFull(i-1,j)  &&  IsFull(i+1,j)) {
         if (!(IsFull(i,j-1)  ||  IsFull(i,j+1)))  return false;
      }
   }
   
   return true;
}



void StoreCandidate(int i, int j)
{
   int k;

   if ((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  &&  (j < N_y)) {
      if (!(ToBeKept(i, j)  &&  IsFull(i, j))) {
         Convert(i, j, k);
         if (num_candidates >= max_num_candidates) {
            max_num_candidates += 131072;
            candidates = (int*) ReallocateArray(candidates, 
                                                max_num_candidates, 
                                                sizeof(int));
         }
         candidates[num_candidates] = k;
         ++num_candidates;
      }
   }

   return;
}



void StoreKeepCandidate(int i, int j)
{
   int k;

   if ((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  &&  (j < N_y)) {
      if (!(ToBeKept(i, j)  &&  IsFull(i, j))) {
         Convert(i, j, k);
         if (num_keep_candidates >= max_num_keep_candidates) {
            max_num_keep_candidates += 131072;
            keep_candidates = (int*) ReallocateArray(keep_candidates, 
                                                     max_num_keep_candidates, 
                                                     sizeof(int));
         }
         keep_candidates[num_keep_candidates] = k;
         ++num_keep_candidates;
      }
   }

   return;
}


void Mark(int i, int j, int *num_cells)
{
   assert((i >= 0)  &&  (i < N_x)  &&  (j >= 0)  && (j < N_y));
   full[i][j] = true;
   ++(*num_cells);

   if (i > 0)        StoreCandidate(i-1,j);
   if (j > 0)        StoreCandidate(i,j-1);
   if (i < (N_x-1))  StoreCandidate(i+1,j);
   if (j < (N_y-1))  StoreCandidate(i,j+1);

   return;
}



void Keep(int i, int j, boolean store_candidates, int *num_keep)
{
   if ((i >= 0)  &&  (j >= 0)  &&  (i < N_x)  &&  (j < N_y)) {
      keep[i][j] = true;
      ++(*num_keep);
      
      if (store_candidates) {
         if (i > 0)        StoreKeepCandidate(i-1,j);
         if (j > 0)        StoreKeepCandidate(i,j-1);
         if (i < (N_x-1))  StoreKeepCandidate(i+1,j);
         if (j < (N_y-1))  StoreKeepCandidate(i,j+1);
      }
   }

   return;
}



void OutputLoop(FILE *output, int i1, int j1, 
                int *total_number, boolean aligned,
                boolean perturb, int loop_cntr,
                boolean diagonal, int smooth)
{
   int number = 0, i0, j0;
   int e_comp(const void *e1, const void *e2);
   int i, j, k, i2, j2, sum = 0;
   vertex_node zero={0,0};
   coord p;
   coord *pnts = (coord*)NULL, *old_pnts = (coord*)NULL;
   int num_pnts = 0, max_num_pnts = 0, old_num_pnts;

   assert((i1 >= 0)  &&  (i1 <= N_x)  &&  (j1 >= 0)  && (j1 <= N_y));

   /*                                                                        */
   /* count the number of edges of this loop                                 */
   /*                                                                        */
   i0 = i1;
   j0 = j1;
   num_vertices = 0;

   do {
      StoreVertex(i1, j1);
      SetVis(i1, j1, true);
      i2 = GetStartI(i1, j1);
      j2 = GetStartJ(i1, j1);
      assert((i2 >= 0)  &&  (i2 <= N_x)  &&  (j2 >= 0)  && (j2 <= N_y));
      if (GetVis(i2, j2)) {
         i2 = GetEndI(i1, j1);
         j2 = GetEndJ(i1, j1);
         if (GetVis(i2, j2)) {
            i2 = i0;
            j2 = j0;
         }
         assert((i2 >= 0)  &&  (i2 <= N_x)  &&  (j2 >= 0)  && (j2 <= N_y));
      }
      i1 = i2;
      j1 = j2;
   } while (!GetVis(i1, j1));
   assert((i0 == i1)  &&  (j0 == j1));
   StoreVertex(i0, j0);
   StoreVertex(i0, j0);
   number = num_vertices - 1;      

   i = 0;
   j = 1;
   k = 2;
   while (k < number) {
      while ((k < number)  &&
             (((vertices[i].i1 == vertices[j].i1)  &&
               (vertices[i].i1 == vertices[k].i1))  ||
              ((vertices[i].j1 == vertices[j].j1)  &&
               (vertices[i].j1 == vertices[k].j1)))) {
         j = k;
         ++k;
      }
      ++i;
      vertices[i].i1 = vertices[j].i1;
      vertices[i].j1 = vertices[j].j1;
      ++j;
      ++k;
   }
   ++i;
   vertices[i].i1 = vertices[j].i1;
   vertices[i].j1 = vertices[j].j1;

   if (diagonal) {
      number = i + 1;
      i = 0;
      j = 1;
      k = 2;
      while (k < number) {
         while ((k < number)                             &&
                ((vertices[i].i1 != vertices[j].i1)  &&
                 (vertices[i].j1 != vertices[j].j1)  &&
                 (vertices[i].i1 != vertices[k].i1)  &&
                 (vertices[i].j1 != vertices[k].j1))     &&
                (((vertices[i].i1 - vertices[j].i1) *
                  (vertices[i].j1 - vertices[k].j1)) ==
                 ((vertices[i].j1 - vertices[j].j1) *
                  (vertices[i].i1 - vertices[k].i1)))) {
            j = k;
            ++k;
         }
         ++i;
         vertices[i].i1 = vertices[j].i1;
         vertices[i].j1 = vertices[j].j1;
         ++j;
         ++k;
      }
      ++i;
      vertices[i].i1 = vertices[j].i1;
      vertices[i].j1 = vertices[j].j1;
   }

   if ((vertices[i-1].i1 == vertices[i].i1)  &&
       (vertices[i-1].j1 == vertices[i].j1))        number = i;
   else if ((vertices[0].i1 == vertices[i].i1)  &&
            (vertices[0].j1 == vertices[i].j1))     number = i + 1;
   else                                             number = i;

   for (i = 1;  i < number;  ++i) {
      sum += Det2D(vertices[i-1], vertices[i], zero);
   }

   if (((loop_cntr == 0)  &&  (sum < 0))  ||
       ((loop_cntr >  0)  &&  (sum > 0))) {
      i = 1;
      j = number - 2;
      while (i < j) {
         Swap(vertices[i], vertices[j], zero);
         ++i;
         --j;
      }
   }

   if (aligned) {
      *total_number += number;
      fprintf(output, "%d\n", number);
      fprintf(output, "%d %d\n", vertices[0].i1, vertices[0].j1);
      for (i = 1;  i < number;  ++i) {
         fprintf(output, "%d %d\n", vertices[i].i1, vertices[i].j1);
      }
      fprintf(output, "\n");
   }
   else {
      max_num_pnts = (number - 1) * (smooth + 1) + 2;
      pnts = (coord*) ReallocateArray(pnts, max_num_pnts, sizeof(coord));
      if (perturb) {
         --number;
         for (i = 0;  i < number;  ++i) {
            p.x = (double) vertices[i].i1 + Perturbation();
            p.y = (double) vertices[i].j1 + Perturbation();
            StorePnt(p);
         }
         StorePnt(pnts[0]);
      }
      else {
         p.x = (double) vertices[0].i1;
         p.y = (double) vertices[0].j1;
         StorePnt(p);
         number -= 2;
         for (i = 1;  i < number;  ++i) {
            if (vertices[i].i1 == vertices[i-1].i1) {
               p.y = vertices[i].j1 + Perturbation();
            }
            else {
               p.x = vertices[i].i1 + Perturbation();
            } 
            StorePnt(p);
         }
         i = number;
         if (vertices[i].i1 == vertices[i-1].i1) {
            p.y = vertices[i].j1;
         }
         else {
            p.x = vertices[i].i1;
         }
         StorePnt(p);
         StorePnt(pnts[0]);
      }
      
      while (smooth > 0) {
         old_pnts = pnts;
         old_num_pnts = num_pnts;
         max_num_pnts *= 2;
         num_pnts = 0;
         pnts = (coord*)NULL;
         pnts = (coord*) ReallocateArray(pnts, max_num_pnts, sizeof(coord));
         for (i = 1;  i < old_num_pnts;  ++i) {
            p.x = (3.0 * old_pnts[i-1].x + old_pnts[i].x) / 4.0;
            p.y = (3.0 * old_pnts[i-1].y + old_pnts[i].y) / 4.0;
            StorePnt(p);
            p.x = (old_pnts[i-1].x + 3.0 * old_pnts[i].x) / 4.0;
            p.y = (old_pnts[i-1].y + 3.0 * old_pnts[i].y) / 4.0;
            StorePnt(p);
         }
         p.x = (3.0 * old_pnts[0].x + old_pnts[1].x) / 4.0;
         p.y = (3.0 * old_pnts[0].y + old_pnts[1].y) / 4.0;
         StorePnt(p);
         free(old_pnts);
         --smooth;
      }

      fprintf(output, "%d\n", num_pnts);
      fprintf(output, "%20.16f %20.16f\n", pnts[0].x, pnts[0].y);
      for (i = 1;  i < num_pnts;  ++i) {
         fprintf(output, "%20.16f %20.16f\n", pnts[i].x, pnts[i].y);
      }
      fprintf(output, "\n");

      *total_number += num_pnts;
   }

   free(pnts);
   pnts = (coord*)NULL;

   return;
}



void InitializeGrid() 
{
   full = MakeBMatrix(N_x, N_y, false);   
   keep = MakeBMatrix(N_x, N_y, false);   

   return;
}



void SelectRandomCells(int num_cells, int max_cells, 
                       int num_keep, int max_keep) 
{
   int k, m, mm, i, j, c;

   while ((num_cells < max_cells)  &&  (num_candidates > 0)) {
      UniformRandom(c, num_candidates);
      k = candidates[c];
      --num_candidates;
      candidates[c] = candidates[num_candidates];
      Invert(k, i, j);
      if (IsPossible(i, j)) {
         Mark(i, j, &num_cells);
      }
      if ((num_keep_candidates > 0)  &&  (num_keep < max_keep)) {
         if ((max_cells - num_cells - 1) > 0) {
            mm = 1 + 2 * (max_keep - num_keep) / (max_cells - num_cells - 1);
         }
         else {
            mm = 1 + 2 * (max_keep - num_keep);
         }
         m  = 0;
         while ((m < mm)  &&  (num_keep_candidates > 0))  {
            UniformRandom(c, num_keep_candidates);
            k = keep_candidates[c];
            --num_keep_candidates;
            keep_candidates[c] = keep_candidates[num_keep_candidates];
            Invert(k, i, j);
            if (!IsFull(i, j)) {
               Keep(i, j, true, &num_keep);
               ++m;               
            }
         }
      }
   }

   // if (num_cells < max_cells) printf("ran out of cells: %d of %d\n", num_cells, max_cells); 

   return;
}


#define R_x 10
#define R_y 10



void Compute(double mark_percent, boolean perturb, boolean aligned,
             int hierarchy, boolean diagonal, int smooth,
             const char *file_name) 
{
   static FILE *output;
   static int i = 0, j = 0, k, m, n, mm, nn;
   static int max_cells = 0, num_keep = 0, max_keep, num_cells;
   static int total_number = 0;
   static int hierarchy_cntr = 1;
   static int loop_cntr, i1, j1;
   static double keep_percent;
   static int N_x_old, N_y_old, ii, jj;
   static boolean **old_full, **old_keep;

   /**************************************************************************/
   /*                                                                        */
   /* generate polygon                                                       */
   /*                                                                        */
   /**************************************************************************/

   /*                                                                        */
   /* allocate the grid                                                      */
   /*                                                                        */
   N = N_x * N_y;

   InitializeGrid();
   
   /*                                                                        */
   /* select cells to be kept                                                */
   /*                                                                        */
   keep_percent = (1.0 - mark_percent) * 0.95;
   max_keep  = (int) (N * keep_percent);
   max_cells = (int) (N * mark_percent);
   num_keep  = 0;
   num_cells = 0;
 
   m = N_y / 10;
   for (i = 0;  i < N_x;  ++i) {
      UniformRandom(k, m);
      for (j = 0;  j <= k;  ++j)            Keep(i, j, false, &num_keep); 
      UniformRandom(k, m);
      for (j = N_y-1;  j >= N_y-k-1;  --j)  Keep(i, j, false, &num_keep); 
   }
   for (j = 0;  j < N_y;  ++j) {
      UniformRandom(k, m);
      for (i = 0;  i <= k;  ++i)            Keep(i, j, false, &num_keep); 
      UniformRandom(k, m);
      for (i = N_x-1;  i >= N_x-k-1;  --i)  Keep(i, j, false, &num_keep); 
   }

   if (num_keep < max_keep) {
      m = (max_keep - num_keep) / 4;
      n = 0;
      for (n = 0;  n < m;  ++n) {
         UniformRandom(k, N);
         Invert(k, i, j);
         Keep(i, j, true, &num_keep);
      }
   }

   UniformRandom(i, N_x);
   UniformRandom(j, N_y);
   
   if (num_keep < max_keep) {
      for (ii = i - N_x / 30;  ii < (i + N_x / 30);  ++ii) {
         for (jj = j - N_y / 30;  jj < (j + N_y / 30);  ++jj) {
            Keep(ii, jj, true, &num_keep);
            if (num_keep >= max_keep) break;
         }
         if (num_keep >= max_keep) break;
      }
   }

   /*                                                                        */
   /* select seed cell(s)                                                    */
   /*                                                                        */
   m  = 7 * N_x / 9;
   n  = 7 * N_y / 9;
   mm = N_x / 9;
   nn = N_y / 9;
   do {
      UniformRandom(i, m);
      UniformRandom(j, n);
      i += mm;
      j += nn;
   } while (!IsPossible(i, j));
   Mark(i, j, &num_cells);

   /*                                                                     */
   /* randomly add cells to already selected cells                        */
   /*                                                                     */
   SelectRandomCells(1, max_cells, num_keep, max_keep);   

   /*                                                                     */
   /* use current polygon as "seed" for a refinement                      */
   /*                                                                     */
   while (hierarchy_cntr <= hierarchy) {
      /*                                                                  */
      /* reset grid data                                                  */
      /*                                                                  */
      old_full = full;
      old_keep = keep;
      N_x_old = N_x;
      N_y_old = N_y;
      N_x *= 3;
      N_y *= 3;
      N   = N_x * N_y;

      InitializeGrid();
      num_candidates      = 0;
      num_keep_candidates = 0;
      max_keep  = N * keep_percent;
      max_cells = N * mark_percent;
      num_keep  = 0;
      num_cells = 0;
      
      /*                                                                  */
      /* copy data from coarse grid to fine grid                          */
      /*                                                                  */
      for (i = 0;  i < N_x_old;  ++i) {
         for (j = 0;  j < N_y_old;  ++j) {
            ii = i * 3;
            jj = j * 3;
            if (IsOldFull(i, j, N_x_old, N_y_old, old_full)) {
               Mark(ii+1, jj+1, &num_cells);
               if (IsOldFull(i-1, j, N_x_old, N_y_old, old_full)  &&
                   IsOldFull(i, j-1, N_x_old, N_y_old, old_full))
                  Mark(ii, jj, &num_cells);
               if (IsOldFull(i-1, j, N_x_old, N_y_old, old_full)  &&
                   IsOldFull(i, j+1, N_x_old, N_y_old, old_full))
                  Mark(ii, jj+2, &num_cells);
               if (IsOldFull(i+1, j, N_x_old, N_y_old, old_full)  &&
                   IsOldFull(i, j+1, N_x_old, N_y_old, old_full))
                  Mark(ii+2, jj+2, &num_cells);
               if (IsOldFull(i+1, j, N_x_old, N_y_old, old_full)  &&
                   IsOldFull(i, j-1, N_x_old, N_y_old, old_full))
                  Mark(ii+2, jj, &num_cells);
               if (IsOldFull(i-1, j, N_x_old, N_y_old, old_full))
                  Mark(ii, jj+1, &num_cells);
               if (IsOldFull(i+1, j, N_x_old, N_y_old, old_full))
                  Mark(ii+2, jj+1, &num_cells);
               if (IsOldFull(i, j-1, N_x_old, N_y_old, old_full))
                  Mark(ii+1, jj, &num_cells);
               if (IsOldFull(i, j+1, N_x_old, N_y_old, old_full))
                  Mark(ii+1, jj+2, &num_cells);
            }
            else {
               Keep(ii+1, jj+1, true, &num_keep);
               if (!IsOldFull(i-1, j, N_x_old, N_y_old, old_full)  &&
                   !IsOldFull(i, j-1, N_x_old, N_y_old, old_full))
                  Keep(ii, jj, true, &num_keep);
               if (!IsOldFull(i-1, j, N_x_old, N_y_old, old_full)  &&
                   !IsOldFull(i, j+1, N_x_old, N_y_old, old_full))
                  Keep(ii, jj+2, true, &num_keep);
               if (!IsOldFull(i+1, j, N_x_old, N_y_old, old_full)  &&
                   !IsOldFull(i, j+1, N_x_old, N_y_old, old_full))
                  Keep(ii+2, jj+2, true, &num_keep);
               if (!IsOldFull(i+1, j, N_x_old, N_y_old, old_full)  &&
                   !IsOldFull(i, j-1, N_x_old, N_y_old, old_full))
                  Keep(ii+2, jj, true, &num_keep);
               if (!IsOldFull(i-1, j, N_x_old, N_y_old, old_full))
                  Keep(ii, jj+1, true, &num_keep);
               if (!IsOldFull(i+1, j, N_x_old, N_y_old, old_full))
                  Keep(ii+2, jj+1, true, &num_keep);
               if (!IsOldFull(i, j-1, N_x_old, N_y_old, old_full))
                  Keep(ii+1, jj, true, &num_keep);
               if (!IsOldFull(i, j+1, N_x_old, N_y_old, old_full))
                  Keep(ii+1, jj+2, true, &num_keep);
            }
         }
      }

      old_full = FreeBMatrix(old_full, N_x_old);
      old_keep = FreeBMatrix(old_keep, N_x_old);
      
      if (hierarchy_cntr <= hierarchy) {
         //         printf("new round: num_keep = %d, max_keep = %d, num_cells = %d, max_cells = %d\n", num_keep, max_keep, num_cells, max_cells);
         SelectRandomCells(num_cells, max_cells, num_keep, max_keep);   
      }
      else {
         //         printf("thinned:   num_keep = %d, max_keep = %d, num_cells = %d, max_cells = %d\n", num_keep, max_keep, num_cells, max_cells);
      }

      hierarchy_cntr += 1;
   }

   /**************************************************************************/
   /*                                                                        */
   /* output polygon                                                         */
   /*                                                                        */
   /**************************************************************************/
   /*                                                                        */
   /* determine boundaries between marked and unmarked cells                 */
   /*                                                                        */
   edges = MakeEMatrix(N_x + 1, N_y + 1);

   if (diagonal) {
      top = MakeBMatrix(N_x, N_y, false);
      bot = MakeBMatrix(N_x, N_y, false);
      lft = MakeBMatrix(N_x, N_y, false);
      rgt = MakeBMatrix(N_x, N_y, false);
      for (i = 0;  i < N_x;  ++i) {
         for (j = 0;  j < N_y;  ++j) {
            if (IsFull(i, j))  {
               SetBot(i, j, true);
               SetTop(i, j, true);
               SetLft(i, j, true);
               SetRgt(i, j, true);
            }
         }
      }
      for (i = 0;  i <= N_x;  ++i) {
         for (j = 0;  j <= N_y;  ++j) {
            if (IsCompletelyFull(i,j)) {
               if      (IsFull(i-1,j) &&
                        IsFull(i-1,j+1) &&
                        IsFull(i,j+1) &&
                        IsCompletelyEmpty(i,j-1) &&
                        IsCompletelyEmpty(i+1,j)) {
                  SetBot(i, j, false);
                  SetRgt(i, j, false);
               }
               else if (IsFull(i+1,j) &&
                        IsFull(i+1,j+1) &&
                        IsFull(i,j+1) &&
                        IsCompletelyEmpty(i,j-1) &&
                        IsCompletelyEmpty(i-1,j)) {
                  SetBot(i, j, false);
                  SetLft(i, j, false);
               }
               else if (IsFull(i+1,j) &&
                        IsFull(i+1,j-1) &&
                        IsFull(i,j-1) &&
                        IsCompletelyEmpty(i,j+1) &&
                        IsCompletelyEmpty(i-1,j)) {
                  SetTop(i, j, false);
                  SetLft(i, j, false);
               }
               else if (IsFull(i-1,j) &&
                        IsFull(i-1,j-1) &&
                        IsFull(i,j-1) &&
                        IsCompletelyEmpty(i,j+1) &&
                        IsCompletelyEmpty(i+1,j)) {
                  SetTop(i, j, false);
                  SetRgt(i, j, false);
               }
            }
         }
      }
      for (i = 0;  i <= N_x;  ++i) {
         for (j = 0;  j <= N_y;  ++j) {
            if (IsCompletelyFull(i, j)) {
               if (IsCompletelyEmpty(i,j-1)) {
                  StoreEdge(i, j, i+1, j);
               }
               if (IsCompletelyEmpty(i-1,j)) {
                  StoreEdge(i, j, i, j+1);
               }
            }
            else if (IsCompletelyEmpty(i,j)) {
               if (IsCompletelyFull(i,j-1)) {
                  StoreEdge(i, j, i+1, j);
               }
               if (IsCompletelyFull(i-1,j)) {
                  StoreEdge(i, j+1, i, j);

               }
            }
            else {
               if (IsTopFull(i, j)) {
                  if (IsRgtFull(i, j)) {
                     StoreEdge(i, j+1, i+1, j);
                  }
                  else if (IsLftFull(i, j)) {
                     StoreEdge(i, j, i+1, j+1);
                  }
               }
               else if (IsBotFull(i, j)) {
                  if (IsLftFull(i, j)) {
                     StoreEdge(i, j+1, i+1, j);
                  }
                  else if (IsRgtFull(i, j)) {
                     StoreEdge(i, j, i+1, j+1);
                  }
               }
               else {
                  printf("cell[%d,%d]: neither top nor bottom is full!\n",
                         i, j);
               }
            }
         }
      }
   }
   else {
      for (i = 0;  i <= N_x;  ++i) {
         for (j = 0;  j <= N_y;  ++j) {
            if (IsFull(i, j)) {
               if (!IsFull(i,j-1)) {
                  StoreEdge(i, j, i+1, j);
               }
               if (!IsFull(i-1,j)) {
                  StoreEdge(i, j, i, j+1);
               }
            }
            else {
               if (IsFull(i,j-1)) {
                  StoreEdge(i, j, i+1, j);
               }
               if (IsFull(i-1,j)) {
                  StoreEdge(i, j, i, j+1);
               }
            }
         }
      }
   }

   for (i1 = 0;  i1 <= N_x;  ++i1) {
      for (j1 = 0;  j1 <= N_y;  ++j1) {
         if (edges[i1][j1].i1 != NIL) {
            assert(edges[i1][j1].j1 != NIL);
            assert(edges[i1][j1].i2 != NIL);
            assert(edges[i1][j1].j2 != NIL);
            assert(((edges[edges[i1][j1].i1][edges[i1][j1].j1].i1 == i1) &&
                    (edges[edges[i1][j1].i1][edges[i1][j1].j1].j1 == j1)) ||
                   ((edges[edges[i1][j1].i1][edges[i1][j1].j1].i2 == i1) &&
                    (edges[edges[i1][j1].i1][edges[i1][j1].j1].j2 == j1)));
         }
      }
   }

   /*                                                                        */
   /* open output file                                                       */
   /*                                                                        */
   OpenFile(&output, file_name);
   loop_cntr = 0;
   for (i1 = 0;  i1 <= N_x;  ++i1) {
      for (j1 = 0;  j1 <= N_y;  ++j1) {
         if (edges[i1][j1].i1 != NIL) {
            assert(edges[i1][j1].j1 != NIL);
            assert(edges[i1][j1].i2 != NIL);
            assert(edges[i1][j1].j2 != NIL);
            assert(((edges[edges[i1][j1].i1][edges[i1][j1].j1].i1 == i1) &&
                    (edges[edges[i1][j1].i1][edges[i1][j1].j1].j1 == j1)) ||
                   ((edges[edges[i1][j1].i1][edges[i1][j1].j1].i2 == i1) &&
                    (edges[edges[i1][j1].i1][edges[i1][j1].j1].j2 == j1)));
         }
         if (!GetVis(i1, j1)) {
            OutputLoop(output, i1, j1, &total_number, aligned, perturb, 
                       loop_cntr, diagonal, smooth);
            ++loop_cntr;
         }
      }
   }
   fclose(output);

   printf("\n%d edges generated within %d loops\n", 
          total_number - loop_cntr, loop_cntr);

   free(vertices);

   return;
}


void PrintHeader(void)
{
   printf("\n");
   printf("***********************************************************\n");
   printf("*                                                         *\n");
   printf("*                      srpg                               *\n");
   printf("*                                                         *\n");
   printf("* Generation of pseudo-random polygonal areas based on a  *\n");
   printf("* regular grid of quadratic cells.                        *\n");
   printf("*                                                         *\n");
   printf("*      Martin Held 2008-2020        held@cs.sbg.ac.at     *\n");
   printf("*                                                         *\n");
   printf("***********************************************************\n");

   return;
}


int main(int argc, char **argv)
{
   int count = 1;
   boolean success = true, perturb = false, aligned = false, diagonal = false;
   boolean name_read = false;
   int hierarchy = 0;
   int smooth = 0, seed = 0;
   double mark_percent = 0.5;
   char *file_name = "srpg.line";

   PrintHeader();

   /*                                                                        */
   /* parse command-line arguments                                           */
   /*                                                                        */
   while ((count < argc)  &&  success)    {

      if (strcmp(argv[count],"--Nx") == 0) {
         ++count;
         if ((success = (count < argc)))  N_x = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--Ny") == 0) {
         ++count;
         if ((success = (count < argc)))  N_y = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--percent") == 0) {
         ++count;
         if ((success = (count < argc)))  mark_percent = atof(argv[count]);
      }
      else if (strcmp(argv[count],"--output") == 0) {
         ++count;
         if ((success = (count < argc))) {
            file_name = argv[count];
            name_read = true;
         }
      }
      else if (strcmp(argv[count],"--perturb") == 0) {
         perturb = true;
      }
      else if (strcmp(argv[count],"--holes") == 0) {
         holes = true;
      }
      else if (strcmp(argv[count],"--aligned") == 0) {
         aligned = true;
      }
      else if (strcmp(argv[count],"--hierarchy") == 0) {
         ++count;
         if ((success = (count < argc)))  hierarchy = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--smooth") == 0) {
         ++count;
         if ((success = (count < argc)))  smooth = atoi(argv[count]);
      }
      else if (strcmp(argv[count],"--seed") == 0) {
         ++count;
         if ((success = (count < argc)))  seed = atof(argv[count]);
      }
      else if (strcmp(argv[count],"--diagonal") == 0) {
         diagonal = true;
      }
      else {
         success = false;
      }
      ++count;
   }

   if (!success  ||  !name_read) {
      if (!success) printf("\nUnrecognized command-line option: %s\n", argv[count-1]); 
         printf("\nUsage: srgp --Nx X --Ny Y --percent P --output XYZ [--seed S] [--holes]\n            [[--aligned] | [--perturb [--smooth M]]\n            [--hierarchy N] [--diagonal]\n       where X,Y,M,N are positive integers,\n             S is a non-negative integer, and\n             0.001<P<0.5 is a real.\n"); 
         exit(2);
   }

   if (aligned)         perturb = false;
   else if (perturb)    aligned = false;
   if (hierarchy < 0)   hierarchy = 0;
   if (smooth < 0)      smooth = 0;
   if (seed < 0)        seed = 0;
   if (N_x < 3)         N_x = 10;
   if (N_y < 3)         N_y = 10;
   if (mark_percent < 0.0001)       mark_percent = 0.0001;
   else if (mark_percent > 0.5)     mark_percent = 0.5;

   /*                                                                        */
   /* initialized random-number generator                                    */
   /*                                                                        */
   srand48(seed);

   if (diagonal && !perturb)  aligned = true;

   /*                                                                        */
   /* call the driver that controls the random-polygon generation            */
   /*                                                                        */
   Compute(mark_percent, perturb, aligned, hierarchy, diagonal, smooth,
           file_name);

   full  = FreeBMatrix(full, N_x);
   keep  = FreeBMatrix(keep, N_x);
   edges = FreeEMatrix(edges, N_x);
   if (diagonal) {
      top = FreeBMatrix(top, N_x);
      bot = FreeBMatrix(bot, N_x);
      lft = FreeBMatrix(lft, N_x);
      rgt = FreeBMatrix(rgt, N_x);
   }

   exit(0);
}
