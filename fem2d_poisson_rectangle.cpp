
// Homework part B
// PES- UPC (2020)
// by: Nino Guzman, Shridharan Suresh

//  Original code obtained from below reference
//  https://people.sc.fsu.edu/~jburkardt/cpp_src/cpp_src.html
//  Author: John Burkardt

//  IMPORTANT:
//    We treat this code on a user-level, to perfom some analysis and comparison
// with the scheme in part A.
//

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

int main ( );

class Node {
    public:
    
    void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] );
    
    void indx_set ( int nx, int ny, int node_num, int indx[], int *nunk );
};
class Element {
    public:
    
    void area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] );
};
class Shape {
    public:
    
    void grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] );
};
class Numerical {
    public:
    
    int dgb_fa ( int n, int ml, int mu, double a[], int pivot[] );
    
    double *dgb_sl ( int n, int ml, int mu, double a[], int pivot[],
  double b[], int job );
    
    void assemble ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int nq,
  double wq[], double xq[], double yq[], double element_area[], int indx[],
  int ib, int nunk, double a[], double f[] );
  
    int bandwidth ( int nnodes, int element_num, int element_node[],
  int node_num, int indx[] );
  
};
class BoundaryConditions {
    public:
    
    void boundary ( int nx, int ny, int node_num, double node_xy[], int indx[],
  int ib, int nunk, double a[], double f[] );
};
class Analysis: public BoundaryConditions {
    public:
    
    void errors ( double element_area[], int element_node[], int indx[],
  double node_xy[], double f[], int element_num, int nnodes,
  int nunk, int node_num, double *el2, double *eh1 );
  
    void quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] );
  
    void exact ( double x, double y, double *u, double *dudx, double *dudy );
    
    void solution_write ( double f[], int indx[], int node_num, int nunk,
  string output_filename, double node_xy[] );
};
class IntegrationPoint: public Numerical, public Analysis {
    public:
    
    void qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *bb, double *bx, double *by );
  
    void quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] );
};
class operations: public Numerical, public BoundaryConditions {
    public:
    int i4_max ( int i1, int i2 );
    
    int i4_min ( int i1, int i2 );
};

void element_write ( int nnodes, int element_num, int element_node[],
  string triangulation_txt_file_name );
  
void i4vec_print_some ( int n, int a[], int max_print, string title );

void nodes_plot ( string file_name, int node_num, double node_xy[],
  bool node_label );
  
void nodes_write ( int node_num, double node_xy[], string output_filename );

double r8_huge ( void );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
int r8_nint ( double x );
//void r8vec_print_some ( int n, double a[], int max_print, string title );
double rhs ( double x, double y );

void triangulation_order6_plot ( string file_name, int node_num, double node_xy[],
  int tri_num, int triangle_node[], int node_show, int triangle_show );

/*void output_results ( double f[], int indx[], int node_num, int nunk,
  string output_vtk, double node_xy[] );*/

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM2D_POISSON_RECTANGLE.
//
//  Discussion:
//
//    FEM2D_POISSON_RECTANGLE solves
//
//      -Laplacian U(X,Y) = F(X,Y)
//
//    in a rectangular region in the plane.  Along the boundary,
//    Dirichlet boundary conditions are imposed.
//
//      U(X,Y) = G(X,Y)
//
//    The code uses continuous piecewise quadratic basis functions on
//    triangles determined by a uniform grid of NX by NY points.
//
//  Local parameters:
//
//    Local, double A[(3*IB+1)*NUNK], the coefficient matrix.
//
//    Local, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
//
//    Local, double C[NUNK], the finite element coefficients, solution of A * C = F.
//
//    Local, double EH1, the H1 seminorm error.
//
//    Local, double EL2, the L2 error.
//
//    Local, int ELEMENT_NODE[ELEMENT_NUM*NNODES]; ELEMENT_NODE(I,J) is the
//    global node index of the local node J in element I.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, double F[NUNK], the right hand side.
//
//    Local, int IB, the half-bandwidth of the matrix.
//
//    Local, int INDX[NODE_NUM], gives the index of the unknown quantity
//    associated with the given node.
//
//    Local, int NNODES, the number of nodes used to form one element.
//
//    Local, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
//
//    Local, int NQ, the number of quadrature points used for assembly.
//
//    Local, int NUNK, the number of unknowns.
//
//    Local, int NX, the number of points in the X direction.
//
//    Local, int NY, the number of points in the Y direction.
//
//    Local, double WQ[NQ], quadrature weights.
//
//    Local, double XL, XR, YB, YT, the X coordinates of
//    the left and right sides of the rectangle, and the Y coordinates
//    of the bottom and top of the rectangle.
//
//    Local, double XQ[NQ*ELEMENT_NUM], YQ[NQ*ELEMENT_NUM], the X and Y
//    coordinates of the quadrature points in each element.
//
{
# define NNODES 6
# define NQ 3
# define NX 7
# define NY 7
# define ELEMENT_NUM ( NX - 1 ) * ( NY - 1 ) * 2
# define NODE_NUM ( 2 * NX - 1 ) * ( 2 * NY - 1 )

  double *a;
  double *c;
  double eh1;
  double el2;
  double element_area[ELEMENT_NUM];
  int element_node[NNODES*ELEMENT_NUM];
  double *f;
  int ib;
  int ierr;
  int indx[NODE_NUM];
  int job;
  string node_eps_file_name = "fem2d_poisson_rectangle_nodes.eps";
  string node_txt_file_name = "fem2d_poisson_rectangle_nodes.txt";
  bool node_label;
  int node_show;
  double node_xy[2*NODE_NUM];
  int nunk;
  int *pivot;
  string solution_txt_file_name = "fem2d_poisson_rectangle_solution.txt";
  int triangle_show;
  string triangulation_eps_file_name = "fem2d_poisson_rectangle_elements.eps";
  string triangulation_txt_file_name = "fem2d_poisson_rectangle_elements.txt";
  double wq[NQ];
  double xl = 0.0E+00;
  double xq[NQ*ELEMENT_NUM];
  double xr = 1.0E+00;
  double yb = 0.0E+00;
  double yq[NQ*ELEMENT_NUM];
  double yt = 1.0E+00;

  //timestamp ( );

  cout << "\n";
  cout << "FEM2D_POISSON_RECTANGLE:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Solution of the Poisson equation on a unit box\n";
  cout << "  in 2 dimensions.\n";
  cout << "\n";
  cout << "  - Uxx - Uyy = F(x,y) in the box\n";
  cout << "       U(x,y) = G(x,y) on the boundary.\n";
  cout << "\n";
  cout << "  The finite element method is used, with piecewise\n";
  cout << "  quadratic basis functions on 6 node triangular\n";
  cout << "  elements.\n";
  cout << "\n";
  cout << "  The corner nodes of the triangles are generated by an\n";
  cout << "  underlying grid whose dimensions are\n";
  cout << "\n";
  cout << "  NX =                 " << NX << "\n";
  cout << "  NY =                 " << NY << "\n";
  cout << "\n";
  cout << "  Number of nodes    = " << NODE_NUM << "\n";
  cout << "  Number of elements = " << ELEMENT_NUM << "\n";
//
//  Set the coordinates of the nodes.
//
Node n1;
Element e1;
Shape s1;
IntegrationPoint i1;
BoundaryConditions b1;
Analysis a1;
Numerical N1;
  n1.xy_set ( NX, NY, NODE_NUM, xl, xr, yb, yt, node_xy );
//
//  Organize the nodes into a grid of 6-node triangles.
//
  s1.grid_t6 ( NX, NY, NNODES, ELEMENT_NUM, element_node );
//
//  Set the quadrature rule for assembly.
//
  i1.quad_a ( node_xy, element_node, ELEMENT_NUM, NODE_NUM,
    NNODES, wq, xq, yq );
//
//  Determine the areas of the elements.
//
  e1.area_set ( NODE_NUM, node_xy, NNODES, ELEMENT_NUM,
    element_node, element_area );
//
//  Determine which nodes are boundary nodes and which have a
//  finite element unknown.  Then set the boundary values.
//
  n1.indx_set ( NX, NY, NODE_NUM, indx, &nunk );

  cout << "  Number of unknowns =       " << nunk << "\n";
//
//  Determine the bandwidth of the coefficient matrix.
//
  ib = N1.bandwidth ( NNODES, ELEMENT_NUM, element_node, NODE_NUM, indx );

  cout << "\n";
//
//  Make an EPS picture of the nodes.
/*
  if ( NX <= 10 && NY <= 10 )
  {
    node_label = true;
    nodes_plot ( node_eps_file_name, NODE_NUM, node_xy, node_label );

    cout << "\n";
    cout << "FEM2D_POISSON_RECTANGLE:\n";
    cout << "  Wrote an EPS file\n";
    cout << "    \"" << node_eps_file_name << "\".\n";
    cout << "  containing a picture of the nodes.\n";
  }
*/
//  Write the nodes to an ASCII file that can be read into MATLAB.
//
  nodes_write ( NODE_NUM, node_xy, node_txt_file_name );

  cout << "\n";
  cout << "FEM2D_POISSON_RECTANGLE:\n";
  cout << "  Wrote an ASCII node file\n";
  cout << "    " << node_txt_file_name << "\n";
  cout << "  of the form\n";
  cout << "    X(I), Y(I)\n";
  cout << "  which can be used for plotting.\n";
//

//  Allocate space for the coefficient matrix A and right hand side F.
//
  a = new double[(3*ib+1)*nunk];
  f = new double[nunk];
  pivot = new int[nunk];
//
//  Assemble the coefficient matrix A and the right-hand side F of the
//  finite element equations.
//
  N1.assemble ( NODE_NUM, node_xy, NNODES,
    ELEMENT_NUM, element_node, NQ,
    wq, xq, yq, element_area, indx, ib, nunk, a, f );

//  Modify the coefficient matrix and right hand side to account for
//  boundary conditions.
//
  b1.boundary ( NX, NY, NODE_NUM, node_xy, indx, ib, nunk, a, f );

//  Solve the linear system using a banded solver.
//
  ierr = N1.dgb_fa ( nunk, ib, ib, a, pivot );

  if ( ierr != 0 )
  {
    cout << "\n";
    cout << "FEM2D_POISSON_RECTANGLE - Error!\n";
    cout << "  DGB_FA returned an error condition.\n";
    cout << "\n";
    cout << "  The linear system was not factored, and the\n";
    cout << "  algorithm cannot proceed.\n";
    exit ( 1 );
  }

  job = 0;
  c = N1.dgb_sl ( nunk, ib, ib, a, pivot, f, job );

//
//  Calculate error
//
  a1.errors ( element_area, element_node, indx, node_xy, c,
    ELEMENT_NUM, NNODES, nunk, NODE_NUM, &el2, &eh1 );

//
//  Write an ASCII file that can be read into MATLAB.
//
  a1.solution_write ( c, indx, NODE_NUM, nunk, solution_txt_file_name,
    node_xy );

  cout << "\n";
  cout << "FEM2D_POISSON_RECTANGLE:\n";
  cout << "  Wrote an ASCII solution file\n";
  cout << "    " << solution_txt_file_name << "\n";
  cout << "  of the form\n";
  cout << "    U( X(I), Y(I) )\n";
  cout << "  which can be used for plotting.\n";
//
//  Deallocate memory.
//
  delete [] a;
  delete [] c;
  delete [] f;
  delete [] pivot;
/*
//  Terminate.
//
  cout << "\n";
  cout << "FEM2D_POISSON_RECTANGLE:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
*/
  return 0;
# undef NNODES
# undef NQ
# undef NX
# undef NY
# undef ELEMENT_NUM
# undef NODE_NUM
}
//****************************************************************************80

void Node::xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
  double yt, double node_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    XY_SET sets the XY coordinates of the nodes.
//
{
  int i;
  int j;

  for ( j = 1; j <= 2*ny-1; j++ )
  {
    for ( i = 1; i <= 2*nx - 1; i++ )
    {
      node_xy[0+(i-1+(j-1)*(2*nx-1))*2] =
        ( ( double ) ( 2 * nx - i - 1 ) * xl
        + ( double ) (          i - 1 ) * xr )
        / ( double ) ( 2 * nx     - 2 );

      node_xy[1+(i-1+(j-1)*(2*nx-1))*2] =
        ( ( double ) ( 2 * ny - j - 1 ) * yb
        + ( double ) (          j - 1 ) * yt )
        / ( double ) ( 2 * ny     - 2 );

    }
  }

  return;
}

//****************************************************************************80

void Element::area_set ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], double element_area[] )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_SET sets the area of each element.
//
{
  int element;
  int i1;
  int i2;
  int i3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  for ( element = 0; element < element_num; element++ )
  {
    i1 = element_node[0+element*nnodes];
    x1 = node_xy[0+(i1-1)*2];
    y1 = node_xy[1+(i1-1)*2];

    i2 = element_node[1+element*nnodes];
    x2 = node_xy[0+(i2-1)*2];
    y2 = node_xy[1+(i2-1)*2];

    i3 = element_node[2+element*nnodes];
    x3 = node_xy[0+(i3-1)*2];
    y3 = node_xy[1+(i3-1)*2];

    element_area[element] = 0.5E+00 * fabs
      ( y1 * ( x2 - x3 )
      + y2 * ( x3 - x1 )
      + y3 * ( x1 - x2 ) );
  }

  return;
}
//****************************************************************************80

void Node::indx_set ( int nx, int ny, int node_num, int indx[], int *nunk )

//****************************************************************************80
//
//  Purpose:
//
//    INDX_SET assigns a boundary value index or unknown value index at each node.
//
{
  int i;
  int in;
  int j;

  *nunk = 0;
  in = 0;

  for ( j = 1; j <= 2 * ny - 1; j++ )
  {
    for ( i = 1; i <= 2 * nx - 1; i++ )
    {
      in = in + 1;
      *nunk = *nunk + 1;
      indx[in-1] = *nunk;
    }
  }

  return;
}
//****************************************************************************80

void Numerical::assemble ( int node_num, double node_xy[], int nnodes,
  int element_num, int element_node[], int nq,
  double wq[], double xq[], double yq[], double element_area[], int indx[],
  int ib, int nunk, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE assembles the matrix and right-hand side using piecewise quadratics.
//
{
  double aij;
  int basis;
  double bi;
  double bj;
  double dbidx;
  double dbidy;
  double dbjdx;
  double dbjdy;
  int element;
  int i;
  int ip;
  int ipp;
  int j;
  int quad;
  int test;
  double w;
  double x;
  double y;
  
  IntegrationPoint i2;
//
//  Initialize the arrays to zero.
//
  for ( i = 1; i <= nunk; i++ )
  {
    f[i-1] = 0.0E+00;
  }

  for ( j = 1; j <= nunk; j++ )
  {
    for ( i = 1; i <= 3*ib + 1; i++ )
    {
      a[i-1+(j-1)*(3*ib+1)] = 0.0E+00;
    }
  }
//
//  The actual values of A and F are determined by summing up
//  contributions from all the elements.
//
  for ( element = 1; element <= element_num; element++ )
  {
    for ( quad = 1; quad <= nq; quad++ )
    {
      x = xq[quad-1+(element-1)*nq];
      y = yq[quad-1+(element-1)*nq];
      w = element_area[element-1] * wq[quad-1];

      for ( test = 1; test <= nnodes; test++ )
      {
        ip = element_node[test-1+(element-1)*nnodes];
        i = indx[ip-1];

        i2.qbf ( x, y, element, test, node_xy, element_node,
          element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        f[i-1] = f[i-1] + w * rhs ( x, y ) * bi;
//
//  We are about to compute a contribution associated with the
//  I-th test function and the J-th basis function, and add this
//  to the entry A(I,J).
//
//  Because of the compressed storage of the matrix, the element
//  will actually be stored in A(I-J+2*IB+1,J).
//
//  Two extra complications: we are storing the array as a vector,
//  and C uses 0-based indices rather than 1-based indices.
//
//  Therefore, we ACTUALLY store the entry in A[I-J+2*IB+1-1 + (J-1) * (3*IB+1)];
//
        for ( basis = 1; basis <= nnodes; basis++ )
        {
          ipp = element_node[basis-1+(element-1)*nnodes];
          j = indx[ipp-1];

          i2.qbf ( x, y, element, basis, node_xy, element_node,
            element_num, nnodes, node_num, &bj, &dbjdx, &dbjdy );

          aij = dbidx * dbjdx + dbidy * dbjdy;

          a[i-j+2*ib+(j-1)*(3*ib+1)] = a[i-j+2*ib+(j-1)*(3*ib+1)] + w * aij;
        }
      }
    }
  }

  return;
}
//****************************************************************************80

int Numerical::bandwidth ( int nnodes, int element_num, int element_node[],
  int node_num, int indx[] )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth of the coefficient matrix.
//
{
  int element;
  int i;
  int iln;
  int in;
  int j;
  int jln;
  int jn;
  int nhba;
  operations o1;
  nhba = 0;

  for ( element = 1; element <= element_num; element++ )
  {
    for ( iln = 1; iln <= nnodes; iln++ )
    {
      in = element_node[iln-1+(element-1)*nnodes];
      i = indx[in-1];
      if ( 0 < i )
      {
        for ( jln = 1; jln <= nnodes; jln++ )
        {
          jn = element_node[jln-1+(element-1)*nnodes];
          j = indx[jn-1];
          nhba = o1.i4_max ( nhba, j - i );
        }
      }
    }
  }

  return nhba;
}
//****************************************************************************80

void BoundaryConditions::boundary ( int nx, int ny, int node_num, double node_xy[], int indx[],
  int ib, int nunk, double a[], double f[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOUNDARY modifies the linear system for boundary conditions.
//
{
  int col;
  double dudx;
  double dudy;
  int i;
  int j;
  int jhi;
  int jlo;
  int node;
  int row;
  double u;
  double x;
  double y;
  operations o2;
  Analysis a2;
//
//  Consider each node.
//
  node = 0;

  for ( row = 1; row <= 2 * ny - 1; row++ )
  {
    for ( col = 1; col <= 2 * nx - 1; col++ )
    {
      node = node + 1;

      if ( row == 1 ||
           row == 2 * ny - 1 ||
           col == 1 ||
           col == 2 * nx - 1 )
      {
        i = indx[node-1];
        x = node_xy[0+(node-1)*2];
        y = node_xy[1+(node-1)*2];
        a2.exact ( x, y, &u, &dudx, &dudy );

        jlo = o2.i4_max ( i - ib, 1 );
        jhi = o2.i4_min ( i + ib, nunk );

        for ( j = jlo; j <= jhi; j++ )
        {
          a[i-j+2*ib+(j-1)*(3*ib+1)] = 0.0;
        }

        a[i-i+2*ib+(i-1)*(3*ib+1)] = 1.0;

        f[i-1] = u;
      }
    }
  }

  return;
}
//****************************************************************************80


int Numerical::dgb_fa ( int n, int ml, int mu, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_FA performs a LINPACK-style PLU factorization of an DGB matrix.
//
{
  int col = 2 * ml + mu + 1;
  int i;
  int i0;
  int j;
  int j0;
  int j1;
  int ju;
  int jz;
  int k;
  int l;
  int lm;
  int m;
  int mm;
  double t;

  operations o3;
  operations o4;
  m = ml + mu + 1;
//
//  Zero out the initial fill-in columns.
//
  j0 = mu + 2;
  j1 = o3.i4_min ( n, m ) - 1;

  for ( jz = j0; jz <= j1; jz++ )
  {
    i0 = m + 1 - jz;
    for ( i = i0; i <= ml; i++ )
    {
      a[i-1+(jz-1)*col] = 0.0E+00;
    }
  }

  jz = j1;
  ju = 0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Zero out the next fill-in column.
//
    jz = jz + 1;
    if ( jz <= n )
    {
      for ( i = 1; i <= ml; i++ )
      {
        a[i-1+(jz-1)*col] = 0.0E+00;
      }
    }
//
//  Find L = pivot index.
//
    lm = o3.i4_min ( ml, n-k );
    l = m;

    for ( j = m+1; j <= m + lm; j++ )
    {
      if ( fabs ( a[l-1+(k-1)*col] ) < fabs ( a[j-1+(k-1)*col] ) )
      {
        l = j;
      }
    }

    pivot[k-1] = l + k - m;
//
//  Zero pivot implies this column already triangularized.
//
    if ( a[l-1+(k-1)*col] == 0.0E+00 )
    {
      cout << "\n";
      cout << "DGB_FA - Fatal error!\n";
      cout << "  Zero pivot on step " << k << "\n";
      return k;
    }
//
//  Interchange if necessary.
//
    t                = a[l-1+(k-1)*col];
    a[l-1+(k-1)*col] = a[m-1+(k-1)*col];
    a[m-1+(k-1)*col] = t;
//
//  Compute multipliers.
//
    for ( i = m+1; i <= m+lm; i++ )
    {
      a[i-1+(k-1)*col] = - a[i-1+(k-1)*col] / a[m-1+(k-1)*col];
    }
//
//  Row elimination with column indexing.
//
    ju = o4.i4_max ( ju, mu + pivot[k-1] );
    ju = o4.i4_min ( ju, n );
    mm = m;

    for ( j = k+1; j <= ju; j++ )
    {
      l = l - 1;
      mm = mm - 1;

      if ( l != mm )
      {
        t                 = a[l-1+(j-1)*col];
        a[l-1+(j-1)*col]  = a[mm-1+(j-1)*col];
        a[mm-1+(j-1)*col] = t;
      }
      for ( i = 1; i <= lm; i++ )
      {
        a[mm+i-1+(j-1)*col] = a[mm+i-1+(j-1)*col]
          + a[mm-1+(j-1)*col] * a[m+i-1+(k-1)*col];
      }
    }
  }

  pivot[n-1] = n;

  if ( a[m-1+(n-1)*col] == 0.0E+00 )
  {
    cout << "\n";
    cout << "DGB_FA - Fatal error!\n";
    cout << "  Zero pivot on step " << n << "\n";
    return n;
  }

  return 0;
}
//****************************************************************************80

double* Numerical::dgb_sl ( int n, int ml, int mu, double a[], int pivot[],
  double b[], int job )

//****************************************************************************80
//
//  Purpose:
//
//    DGB_SL solves a system factored by DGB_FA.
{
  int col = 2 * ml + mu + 1;
  int i;
  int k;
  int l;
  int la;
  int lb;
  int lm;
  int m;
  double t;
  double *x;
  operations o5; operations o6; operations o7; operations o8;
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }
//
  m = mu + ml + 1;
//
//  Solve A * x = b.
//
  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    if ( 1 <= ml )
    {
      for ( k = 1; k <= n-1; k++ )
      {
        lm = o5.i4_min ( ml, n-k );
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
        for ( i = 1; i <= lm; i++ )
        {
          x[k+i-1] = x[k+i-1] + x[k-1] * a[m+i-1+(k-1)*col];
        }
      }
    }
//
//  Solve U * X = Y.
//
    for ( k = n; 1 <= k; k-- )
    {
      x[k-1] = x[k-1] / a[m-1+(k-1)*col];
      lm = o6.i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[lb+i-1] = x[lb+i-1] - x[k-1] * a[la+i-1+(k-1)*col];
      }
    }
  }
//
//  Solve A' * X = B.
//
  else
  {
//
//  Solve U' * Y = B.
//
    for ( k = 1; k <= n; k++ )
    {
      lm = o7.i4_min ( k, m ) - 1;
      la = m - lm;
      lb = k - lm;
      for ( i = 0; i <= lm-1; i++ )
      {
        x[k-1] = x[k-1] - x[lb+i-1] * a[la+i-1+(k-1)*col];
      }
      x[k-1] = x[k-1] / a[m-1+(k-1)*col];
    }
//
//  Solve L' * X = Y.
//
    if ( 1 <= ml )
    {
      for ( k = n-1; 1 <= k; k-- )
      {
        lm = o8.i4_min ( ml, n-k );
        for ( i = 1; i <= lm; i++ )
        {
          x[k-1] = x[k-1] + x[k+i-1] * a[m+i-1+(k-1)*col];
        }
        l = pivot[k-1];

        if ( l != k )
        {
          t      = x[l-1];
          x[l-1] = x[k-1];
          x[k-1] = t;
        }
      }
    }
  }

  return x;
}
//****************************************************************************80

void element_write ( int nnodes, int element_num, int element_node[],
  string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_WRITE writes the elements to a file.
//
{
  int element;
  int i;
  ofstream output;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "ELEMENT_WRITE - Warning!\n";
    cout << "  Could not write the node file.\n";
    return;
  }

  for ( element = 0; element < element_num; element++ )
  {
    for ( i = 0; i < nnodes; i++ )
    {
      output << setw(8)  << element_node[i+element*nnodes] << "  ";
    }
    output << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

void Analysis::errors ( double element_area[], int element_node[], int indx[],
  double node_xy[], double f[], int element_num, int nnodes,
  int nunk, int node_num, double *el2, double *eh1 )

//****************************************************************************80
//
//  Purpose:
//
//    ERRORS calculates the error in the L2
//
{
# define NQE 13

  double ar;
  double bi;
  double dbidx;
  double dbidy;
  double dudx;
  double dudxh;
  double dudy;
  double dudyh;
  int element;
  int i;
  int in1;
  int ip;
  int quad;
  double u;
  double uh;
  double wqe[NQE];
  double x;
  double xqe[NQE];
  double y;
  double yqe[NQE];

  *el2 = 0.0E+00;
  IntegrationPoint i3;

//
//  For each element, retrieve the nodes, area, quadrature weights,
//  and quadrature points.
//
  for ( element = 1; element <= element_num; element++ )
  {
    quad_e ( node_xy, element_node, element, element_num,
      nnodes, node_num, NQE, wqe, xqe, yqe );
//
//  For each quadrature point, evaluate the computed solution and its X and
//  Y derivatives.
//
    for ( quad = 1; quad <= NQE; quad++ )
    {
      ar = element_area[element-1] * wqe[quad-1];
      x = xqe[quad-1];
      y = yqe[quad-1];

      uh = 0.0E+00;
      dudxh = 0.0E+00;
      dudyh = 0.0E+00;

      for ( in1 = 1; in1 <= nnodes; in1++ )
      {
        ip = element_node[in1-1+(element-1)*nnodes];

        i3.qbf (x, y, element, in1, node_xy,
          element_node, element_num, nnodes, node_num, &bi, &dbidx, &dbidy );

        i = indx[ip-1];

        uh    = uh    + bi    * f[i-1];
        dudxh = dudxh + dbidx * f[i-1];
        dudyh = dudyh + dbidy * f[i-1];
      }
//
//  Evaluate the exact solution and its X and Y derivatives.
//
      exact ( x, y, &u, &dudx, &dudy );
//
//  Add the weighted value at this quadrature point to the quadrature sum.
//
      *el2 = *el2 + ar *   pow ( ( uh  - u  ), 2 );

    }
  }

  *el2 = sqrt ( *el2 );

  cout << "\n";
  cout << "   ERROR:                                  \n";
  cout << "    L2 error =          " << setw(14) << *el2 << "     \n";

  return;
# undef NQE
}
//****************************************************************************80

void Analysis::exact ( double x, double y, double *u, double *dudx, double *dudy )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT calculates the exact solution and its first derivatives.
//
{
# define PI 3.14159265358979323846264338327950288419716939937510

  *u    =      sin ( PI * x ) * sin ( PI * y ) + x;
  *dudx = PI * cos ( PI * x ) * sin ( PI * y ) + 1.0E+00;
  *dudy = PI * sin ( PI * x ) * cos ( PI * y );

  return;
# undef PI
}
//****************************************************************************80

void Shape::grid_t6 ( int nx, int ny, int nnodes, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T6 produces a grid of pairs of 6 node triangles.
//
//  Example:
//
//    Input:
//
//      NX = 4, NY = 3
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  3, 15,  2,  9,  8;
//        17, 15,  3, 16,  9, 10;
//         3,  5, 17,  4, 11, 10;
//        19, 17,  5, 18, 11, 12;
//         5,  7, 19,  6, 13, 12;
//        21, 19,  7, 20, 13, 14;
//        15, 17, 29, 16, 23, 22;
//        31, 29, 17, 30, 23, 24;
//        17, 19, 31, 18, 25, 24;
//        33, 31, 19, 32, 25, 26;
//        19, 21, 33, 20, 27, 26;
//        35, 33, 21, 34, 27, 28.
//
//  Diagram:
//
//   29-30-31-32-33-34-35
//    |\ 8  |\10  |\12  |
//    | \   | \   | \   |
//   22 23 24 25 26 27 28
//    |   \ |   \ |   \ |
//    |  7 \|  9 \| 11 \|
//   15-16-17-18-19-20-21
//    |\ 2  |\ 4  |\ 6  |
//    | \   | \   | \   |
//    8  9 10 11 12 13 14
//    |   \ |   \ |   \ |
//    |  1 \|  3 \|  5 \|
//    1--2--3--4--5--6--7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, controls the number of elements along the
//    X and Y directions.  The number of elements will be
//    2 * ( NX - 1 ) * ( NY - 1 ).
//
//    Input, int NNODES, the number of local nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Output, int ELEMENT_NODE[NNODES*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the index of the I-th node of the J-th element.
//
{
  int c;
  int e;
  int element;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element = 0;

  for ( j = 1; j <= ny - 1; j++ )
  {
    for ( i = 1; i <= nx - 1; i++ )
    {
      sw = ( j - 1 ) * 2 * ( 2 * nx - 1 ) + 2 * i - 1;
      w  = sw + 1;
      nw = sw + 2;

      s  = sw + 2 * nx - 1;
      c  = s  + 1;
      n  = s  + 2;

      se = s  + 2 * nx - 1;
      e  = se + 1;
      ne = se + 2;

      element = element + 1;
      element_node[0+(element-1)*nnodes] = sw;
      element_node[1+(element-1)*nnodes] = se;
      element_node[2+(element-1)*nnodes] = nw;
      element_node[3+(element-1)*nnodes] = s;
      element_node[4+(element-1)*nnodes] = c;
      element_node[5+(element-1)*nnodes] = w;

      element = element + 1;
      element_node[0+(element-1)*nnodes] = ne;
      element_node[1+(element-1)*nnodes] = nw;
      element_node[2+(element-1)*nnodes] = se;
      element_node[3+(element-1)*nnodes] = n;
      element_node[4+(element-1)*nnodes] = c;
      element_node[5+(element-1)*nnodes] = e;
    }
  }

  return;
}
//****************************************************************************80

int operations::i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int operations::i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

void i4vec_print_some ( int n, int a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT_SOME prints "some" of an I4VEC.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(10) << a[i] << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(10) << a[i]  << "\n";
    }
    cout << "......  ..............\n";
    i = n - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(10) << a[i]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print-1; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(10) << a[i]  << "\n";
    }
    i = max_print - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(10) << a[i]  << "...more entries...\n";
  }

  return;
}

//****************************************************************************80

void nodes_write ( int node_num, double node_xy[], string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    NODES_WRITE writes the nodes to a file.
//
{
  int node;
  ofstream output;
  double x;
  double y;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "NODES_WRITE - Warning!\n";
    cout << "  Could not write the node file.\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    output << setw(8)  << x << "  "
           << setw(8)  << y << "\n";
  }

  output.close ( );

  return;
}
//****************************************************************************80

void IntegrationPoint::qbf ( double x, double y, int element, int inode, double node_xy[],
  int element_node[], int element_num, int nnodes,
  int node_num, double *b, double *dbdx, double *dbdy )

//****************************************************************************80
//
//  Purpose:
//
//    QBF evaluates the quadratic basis functions.
//
{
  double dbdr;
  double dbds;
  double det;
  double drdx;
  double drdy;
  double dsdx;
  double dsdy;
  int i;
  double r;
  double s;
  double xn[6];
  double yn[6];

  for ( i = 0; i < 6; i++ )
  {
    xn[i] = node_xy[0+(element_node[i+(element-1)*nnodes]-1)*2];
    yn[i] = node_xy[1+(element_node[i+(element-1)*nnodes]-1)*2];
  }
//
//  Determine the (R,S) coordinates corresponding to (X,Y).
//
//  What is happening here is that we are solving the linear system:
//
//    ( X2-X1  X3-X1 ) * ( R ) = ( X - X1 )
//    ( Y2-Y1  Y3-Y1 )   ( S )   ( Y - Y1 )
//
//  by computing the inverse of the coefficient matrix and multiplying
//  it by the right hand side to get R and S.
//
//  The values of dRdX, dRdY, dSdX and dSdY are easily from the formulas
//  for R and S.
//
  det =   ( xn[1] - xn[0] ) * ( yn[2] - yn[0] )
        - ( xn[2] - xn[0] ) * ( yn[1] - yn[0] );

  r = ( ( yn[2] - yn[0] ) * ( x     - xn[0] )
      + ( xn[0] - xn[2] ) * ( y     - yn[0] ) ) / det;

  drdx = ( yn[2] - yn[0] ) / det;
  drdy = ( xn[0] - xn[2] ) / det;

  s = ( ( yn[0] - yn[1] ) * ( x     - xn[0] )
      + ( xn[1] - xn[0] ) * ( y     - yn[0] ) ) / det;

  dsdx = ( yn[0] - yn[1] ) / det;
  dsdy = ( xn[1] - xn[0] ) / det;
//
//  The basis functions can now be evaluated in terms of the
//  reference coordinates R and S.  It's also easy to determine
//  the values of the derivatives with respect to R and S.
//
  if ( inode == 1 )
  {
    *b   =   2.0E+00 *     ( 1.0E+00 - r - s ) * ( 0.5E+00 - r - s );
    dbdr = - 3.0E+00 + 4.0E+00 * r + 4.0E+00 * s;
    dbds = - 3.0E+00 + 4.0E+00 * r + 4.0E+00 * s;
  }
  else if ( inode == 2 )
  {
    *b   =   2.0E+00 * r * ( r - 0.5E+00 );
    dbdr = - 1.0E+00 + 4.0E+00 * r;
    dbds =   0.0E+00;
  }
  else if ( inode == 3 )
  {
    *b   =   2.0E+00 * s * ( s - 0.5E+00 );
    dbdr =   0.0E+00;
    dbds = - 1.0E+00               + 4.0E+00 * s;
  }
  else if ( inode == 4 )
  {
    *b   =   4.0E+00 * r * ( 1.0E+00 - r - s );
    dbdr =   4.0E+00 - 8.0E+00 * r - 4.0E+00 * s;
    dbds =           - 4.0E+00 * r;
  }
  else if ( inode == 5 )
  {
    *b   =   4.0E+00 * r * s;
    dbdr =                           4.0E+00 * s;
    dbds =             4.0E+00 * r;
  }
  else if ( inode == 6 )
  {
    *b   =   4.0E+00 * s * ( 1.0E+00 - r - s );
    dbdr =                         - 4.0E+00 * s;
    dbds =   4.0E+00 - 4.0E+00 * r - 8.0E+00 * s;
  }
  else
  {
    cout << "\n";
    cout << "QBF - Fatal error!\n";
    cout << "  Request for local basis function INODE = " << inode << "\n";
    exit ( 1 );
  }
//
//  We need to convert the derivative information from (R(X,Y),S(X,Y))
//  to (X,Y) using the chain rule.
//
  *dbdx = dbdr * drdx + dbds * dsdx;
  *dbdy = dbdr * drdy + dbds * dsdy;

  return;
}
//****************************************************************************80

void IntegrationPoint::quad_a ( double node_xy[], int element_node[],
  int element_num, int node_num, int nnodes, double wq[], double xq[],
  double yq[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_A sets the quadrature rule for assembly.
//
{
  int element;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  wq[0] = 1.0E+00 / 3.0E+00;
  wq[1] = wq[0];
  wq[2] = wq[0];

  for ( element = 1; element <= element_num; element++ )
  {
    ip1 = element_node[0+(element-1)*nnodes];
    ip2 = element_node[1+(element-1)*nnodes];
    ip3 = element_node[2+(element-1)*nnodes];

    x1 = node_xy[0+(ip1-1)*2];
    x2 = node_xy[0+(ip2-1)*2];
    x3 = node_xy[0+(ip3-1)*2];

    y1 = node_xy[1+(ip1-1)*2];
    y2 = node_xy[1+(ip2-1)*2];
    y3 = node_xy[1+(ip3-1)*2];

    xq[0+(element-1)*3] = 0.5E+00 * ( x1 + x2 );
    xq[1+(element-1)*3] = 0.5E+00 * ( x2 + x3 );
    xq[2+(element-1)*3] = 0.5E+00 * ( x1 + x3 );

    yq[0+(element-1)*3] = 0.5E+00 * ( y1 + y2 );
    yq[1+(element-1)*3] = 0.5E+00 * ( y2 + y3 );
    yq[2+(element-1)*3] = 0.5E+00 * ( y1 + y3 );
  }

  return;
}
//****************************************************************************80

void Analysis::quad_e ( double node_xy[], int element_node[],
  int element, int element_num, int nnodes, int node_num, int nqe,
  double wqe[], double xqe[], double yqe[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_E sets a quadrature rule for the error calculation.
//
{
  int i;
  int ii;
  int iii;
  int ip1;
  int ip2;
  int ip3;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;
  double z1;
  double z2;
  double z3;
  double z4;
  double z5;
  double z6;
  double z7;

  for ( i = 1; i <= 3; i++ )
  {
    wqe[i-1] = 0.175615257433204E+00;
    ii = i + 3;
    wqe[ii-1] = 0.053347235608839E+00;
    ii = i + 6;
    iii = ii + 3;
    wqe[ii-1] = 0.077113760890257E+00;
    wqe[iii-1] = wqe[ii-1];
  }

  wqe[13-1] = -0.14957004446767E+00;

  z1 = 0.479308067841923E+00;
  z2 = 0.260345966079038E+00;
  z3 = 0.869739794195568E+00;
  z4 = 0.065130102902216E+00;
  z5 = 0.638444188569809E+00;
  z6 = 0.312865496004875E+00;
  z7 = 0.048690315425316E+00;

  ip1 = element_node[0+(element-1)*nnodes];
  ip2 = element_node[1+(element-1)*nnodes];
  ip3 = element_node[2+(element-1)*nnodes];

  x1 = node_xy[0+(ip1-1)*2];
  x2 = node_xy[0+(ip2-1)*2];
  x3 = node_xy[0+(ip3-1)*2];

  y1 = node_xy[1+(ip1-1)*2];
  y2 = node_xy[1+(ip2-1)*2];
  y3 = node_xy[1+(ip3-1)*2];

  xqe[ 1-1] = z1 * x1 + z2 * x2 + z2 * x3;
  yqe[ 1-1] = z1 * y1 + z2 * y2 + z2 * y3;
  xqe[ 2-1] = z2 * x1 + z1 * x2 + z2 * x3;
  yqe[ 2-1] = z2 * y1 + z1 * y2 + z2 * y3;
  xqe[ 3-1] = z2 * x1 + z2 * x2 + z1 * x3;
  yqe[ 3-1] = z2 * y1 + z2 * y2 + z1 * y3;
  xqe[ 4-1] = z3 * x1 + z4 * x2 + z4 * x3;
  yqe[ 4-1] = z3 * y1 + z4 * y2 + z4 * y3;
  xqe[ 5-1] = z4 * x1 + z3 * x2 + z4 * x3;
  yqe[ 5-1] = z4 * y1 + z3 * y2 + z4 * y3;
  xqe[ 6-1] = z4 * x1 + z4 * x2 + z3 * x3;
  yqe[ 6-1] = z4 * y1 + z4 * y2 + z3 * y3;
  xqe[ 7-1] = z5 * x1 + z6 * x2 + z7 * x3;
  yqe[ 7-1] = z5 * y1 + z6 * y2 + z7 * y3;
  xqe[ 8-1] = z5 * x1 + z7 * x2 + z6 * x3;
  yqe[ 8-1] = z5 * y1 + z7 * y2 + z6 * y3;
  xqe[ 9-1] = z6 * x1 + z5 * x2 + z7 * x3;
  yqe[ 9-1] = z6 * y1 + z5 * y2 + z7 * y3;
  xqe[10-1] = z6 * x1 + z7 * x2 + z5 * x3;
  yqe[10-1] = z6 * y1 + z7 * y2 + z5 * y3;
  xqe[11-1] = z7 * x1 + z5 * x2 + z6 * x3;
  yqe[11-1] = z7 * y1 + z5 * y2 + z6 * y3;
  xqe[12-1] = z7 * x1 + z6 * x2 + z5 * x3;
  yqe[12-1] = z7 * y1 + z6 * y2 + z5 * y3;
  xqe[13-1] = ( x1 + x2 + x3 ) / 3.0;
  yqe[13-1] = ( y1 + y2 + y3 ) / 3.0;

  return;
}
//****************************************************************************80

double r8_huge ( void )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
{
  return ( double ) HUGE_VAL;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
//****************************************************************************80

void r8vec_print_some ( int n, double a[], int max_print, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_SOME prints "some" of an R8VEC.
//
{
  int i;

  if ( max_print <= 0 )
  {
    return;
  }

  if ( n <= 0 )
  {
    return;
  }

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  if ( n <= max_print )
  {
    for ( i = 0; i < n; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(14) << a[i]  << "\n";
    }
  }
  else if ( 3 <= max_print )
  {
    for ( i = 0; i < max_print-2; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(14) << a[i]  << "\n";
    }

    cout << "......  ..............\n";
    i = n - 1;
    cout << setw(6)  << i + 1 << "  "
         << setw(14) << a[i]  << "\n";
  }
  else
  {
    for ( i = 0; i < max_print-1; i++ )
    {
      cout << setw(6)  << i + 1 << "  "
           << setw(14) << a[i]  << "\n";
    }

    i = max_print - 1;

    cout << setw(6)  << i + 1 << "  "
         << setw(14) << a[i]  << "...more entries...\n";
  }

  return;
}
//****************************************************************************80

double rhs ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    RHS gives the right-hand side of the differential equation.
//
{
# define PI 3.14159265358979323846264338327950288419716939937510

  double value;

  value = 2.0E+00 * PI * PI * sin ( PI * x ) * sin ( PI * y );

  return value;
# undef PI
}
//****************************************************************************80

void Analysis::solution_write ( double f[], int indx[], int node_num, int nunk,
  string output_filename, double node_xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    SOLUTION_WRITE writes the solution to a file.
//
{
  double dudx;
  double dudy;
  int node;
  ofstream output;
  double u;
  double x;
  double y;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cout << "\n";
    cout << "SOLUTION_WRITE - Warning!\n";
    cout << "  Could not write the solution file.\n";
    return;
  }

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    if ( 0 < indx[node] )
    {
      u = f[indx[node]-1];
    }
    else
    {
      exact ( x, y, &u, &dudx, &dudy );
    }

    output << setw(14) << u << "\n";
  }

  output.close ( );

  return;
}

/*
void output_results ( double f[], int indx[], int node_num, int nunk,
  string output_vtk, double node_xy[] )

  double dudx;
  double dudy;
  int node;
  ofstream output;
  double u;
  double x;
  double y;

  output.open ( output_vtk.c_str ( ) );
{


  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    if ( 0 < indx[node] )
    {
      u = f[indx[node]-1];
    }

    output << setw(8)  << x << "  "
           << setw(8)  << y << "\n";
           << setw(8)  << u << "\n";
  }

  output.close ( );

  //Write results to VTK file
  std::ofstream output1("solution.vtk");
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector(D, nodal_solution_names, DataOut<dim>::type_dof_data,
               nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}
*/



