////////////////////////////////////////////
// MPI Life 0.9
// Copyright 2002, David Joiner and
//   The Shodor Education Foundation, Inc.
////////////////////////////////////////////
#include <stdlib.h>
#include <stdio.h>
#include "Life.h"
#include "mpi.h"
#include "mpidata.h"
#include <time.h>
#include <omp.h>
#include <X11/Xlib.h> // Every Xlib program must include this
#include <assert.h>   // I include this to test return values the lazy way
#include <unistd.h>   // So we got the profile for 10 seconds
#define NIL (0)       // A name for the void pointer
mpidata g_mpi;
int main(int argc, char ** argv) {
    extern void setup_mpi(int*, char***,mpidata *);
    extern void finalize_mpi();
    double wtime;
    // Set up MPI
    setup_mpi(&argc,&argv,&g_mpi);
    int rank=g_mpi.rank;
    int size=g_mpi.size;
    random_initByTime(rank);
    // defaults
    int ngrid=105;
    int ncols=ngrid;
    int nrows=ngrid;
    int max_count=10000;
    int do_display=1;
    
    // command line arguments
    if (argc > 1) {
        sscanf(argv[1],"%d",&nrows);
    }
    if (argc > 2) {
        sscanf(argv[2],"%d",&ncols);
    }
    if (argc > 3) {
        sscanf(argv[3],"%d",&max_count);
    }
    if (argc > 4) {
        sscanf(argv[4],"%d",&do_display);
    }
    if (do_display!=0) do_display=1;
    // I'm having a problem with non-square grids, there
    //  might be an nrows and ncols mixed up somewhere down
    //  the line. DAJ 03-10-02
    int ** grid;
    int ** next_grid;
    allocate_grid(ncols, nrows, &grid);
    allocate_grid(ncols, nrows, &next_grid);
    randomize_grid(ncols, nrows, grid, 0.25);
    if (do_display==1) 
        setupWindow(ncols, nrows);
    wtime = MPI_Wtime();
    bool done=false;
    int count=0;
    while(!done) {
        if (count++>max_count) done=true;
        // output
        if (count%1==0&&do_display==1) doDraw(rank,ncols,nrows,grid);
        do_step(rank,size,ncols,nrows,grid,next_grid);
        #pragma omp parallel num_threads(8)
        #pragma omp for
        for (int i=0;i<ncols+2;i++) {
            for (int j=0;j<nrows+2;j++) {
                grid[i][j]=next_grid[i][j];
            }
        }
    }
    wtime = MPI_Wtime() - wtime;
    if (rank==0) {
        printf("\n");
        printf("Number of nodes %d\n",size);
        printf("Problem size: %d by %d, for %d steps\n",
                ncols, nrows, max_count);
        printf("Elapsed wall clock time (seconds) %f\n", wtime);
        printf("Cell update rate (Mcells/second) %f\n",
                (double)ncols*nrows*max_count/(wtime*1000000));
    }
    cleanup_grid(ncols, nrows, &grid);
    cleanup_grid(ncols, nrows, &next_grid);
    finalize_mpi();
}
void do_step(int rank, int size, int ncols, int nrows, int ** grid,
        int ** next_grid) {
    
    //assume 2X2 grid? Side by side grid? Go with side by side, easiest
    // for now, arbitrary number, can program for size=1.
    
    // top and bottom we get from current cell.
    //left right and corners we get from neighboring grids.
    // start off with non blocking sends of each "row"
    // left is rank - 1 % size, right is rank + 1 % size.
    
    // copy sides
    int left_rank = (rank-1+size)%size;
    int right_rank = (rank+1)%size;
    
    if (left_rank>=rank) {
        MPI_Send(grid[1],nrows+2,MPI_INT,left_rank,MPIDATA_TOLEFT,
            MPI_COMM_WORLD);
        MPI_Recv(grid[ncols+1],nrows+2,MPI_INT,right_rank,
            MPIDATA_TOLEFT,
            MPI_COMM_WORLD, &g_mpi.status);
    } else {
        MPI_Recv(grid[ncols+1],nrows+2,MPI_INT,right_rank,
            MPIDATA_TOLEFT,MPI_COMM_WORLD, &g_mpi.status);
        MPI_Send(grid[1],nrows+2,MPI_INT,left_rank,MPIDATA_TOLEFT,
            MPI_COMM_WORLD);
    }
    
    if (right_rank>=rank) {
        MPI_Send(grid[ncols],nrows+2,MPI_INT,right_rank,MPIDATA_TOLEFT,
            MPI_COMM_WORLD);
        MPI_Recv(grid[0],nrows+2,MPI_INT,left_rank,
            MPIDATA_TOLEFT,
            MPI_COMM_WORLD, &g_mpi.status);
    } else {
        MPI_Recv(grid[0],nrows+2,MPI_INT,left_rank,
            MPIDATA_TOLEFT,
            MPI_COMM_WORLD, &g_mpi.status);
        MPI_Send(grid[ncols],nrows+2,MPI_INT,right_rank,MPIDATA_TOLEFT,
            MPI_COMM_WORLD);
    }
    
    // copy corners
        grid[0][0]=grid[0][nrows];
        grid[0][nrows+1]=grid[0][1];
        grid[ncols+1][0]=grid[ncols+1][nrows];
        grid[ncols+1][nrows+1]=grid[ncols+1][1];

        //copy top and bottom
        for (int i=1;i<=ncols;i++) {
                grid[i][0]=grid[i][nrows];
                grid[i][nrows+1]=grid[i][1];
        }
        //update
        #pragma omp parallel num_threads(8)
        #pragma omp for
        for (int i=1;i<=ncols;i++) {
                for (int j=1;j<=nrows;j++) {
                        int neighbors=0;
                        for (int k=i-1;k<=i+1; k++) {
                                for (int l=j-1;l<=j+1; l++) {
                                        if (!(k==i&&l==j)&&grid[k][l]>0) {
                                                neighbors++;
                                        }
                                }
                                if (neighbors>3) continue;
                        }
                        if (neighbors<2||neighbors>3) {
                                next_grid[i][j]=0;
                        } else if (grid[i][j]>0||neighbors==3) {
                                next_grid[i][j]=grid[i][j]+1;
                        }
                }
        }
}
void allocate_grid(int ncols, int nrows, int *** grid){
    (*grid) = new  int * [ncols+2];
    for (int i=0; i<ncols+2;i++) {
        (*grid)[i] = new int[nrows+2];
        for (int j=0;j<nrows+2;j++) {
            (*grid)[i][j]=0;
        }
    }
}
void cleanup_grid(int ncols, int nrows, int *** grid){
    for (int i=0;i<ncols+2;i++) {
        delete (*grid)[i];
    }
    delete (*grid);
}
void randomize_grid(int ncols, int nrows, int ** grid, double prob){
    #pragma omp parallel num_threads(8)
    #pragma omp for
    for (int i=1;i<=ncols;i++) {
        for (int j=1;j<=nrows;j++) {
            if (rand_double()<prob) {
                grid[i][j]=1;
            }
        }
    }
}
double rand_double() {
    return (double)rand()/(double)RAND_MAX;
}
// X information, at some point this should be cleaned up so
// that it does not use global variables
// setupWindow modified from the tutorial on
// http://tronche.com/gui/x/xlib-tutorial/
// by Christophe Tronche
Display *dpy;
int blackColor;
int whiteColor;
Window w;
GC gc;
Pixmap buffer;
Colormap theColormap;
int numXGrayscale=10;
XColor Xgrayscale[10];
int IMAGE_WIDTH=LIFE_IMAGE_WIDTH;
int IMAGE_HEIGHT=LIFE_IMAGE_HEIGHT;
void setupWindow(int ncols, int nrows) {
      // Open the display
      dpy = XOpenDisplay(NIL);
      assert(dpy);
      // Get some colors
      blackColor = BlackPixel(dpy, DefaultScreen(dpy));
      whiteColor = WhitePixel(dpy, DefaultScreen(dpy));
      // Create the window
      if (nrows>ncols) {
         printf ("Howdy\n");
         IMAGE_WIDTH = (int)((double)LIFE_IMAGE_WIDTH*(double)ncols/(double)nrows);
         IMAGE_HEIGHT = LIFE_IMAGE_HEIGHT;
      } else {
         IMAGE_HEIGHT = (int)((double)LIFE_IMAGE_HEIGHT*(double)nrows/(double)ncols);
         IMAGE_WIDTH = LIFE_IMAGE_WIDTH;
      }
      
      w = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0, 0, 
                                     IMAGE_WIDTH, IMAGE_HEIGHT, 0, blackColor,
                                     blackColor);
      buffer = XCreatePixmap(dpy,DefaultRootWindow(dpy),
          IMAGE_WIDTH,IMAGE_HEIGHT,DefaultDepth(dpy,
          DefaultScreen(dpy)));
          
      theColormap = XCreateColormap(dpy, DefaultRootWindow(dpy),
          DefaultVisual(dpy,DefaultScreen(dpy)), AllocNone);
          
      for (int i=0;i<numXGrayscale;i++) {
          int color = (int)((double)i*35535.0/(double)numXGrayscale)+30000;
          Xgrayscale[i].red=color;
          Xgrayscale[i].green=color;
          Xgrayscale[i].blue=color;
          XAllocColor(dpy,theColormap,&(Xgrayscale[i]));
      }
      // We want to get MapNotify events
      XSelectInput(dpy, w, StructureNotifyMask);
      // "Map" the window (that is, make it appear on the screen)
      XMapWindow(dpy, w);
      // Create a "Graphics Context"
      gc = XCreateGC(dpy, w, 0, NIL);
      // Tell the GC we draw using the white color
      XSetForeground(dpy, gc, whiteColor);
      // Wait for the MapNotify event
      for(;;) {
            XEvent e;
            XNextEvent(dpy, &e);
            if (e.type == MapNotify)
                  break;
      }
}
void doDraw(int rank, int ncols, int nrows, int ** grid) {
    int x1,x2,y1,y2; 
    char string[2];
    sprintf(string,"%d",rank);
    
    XSetForeground(dpy, gc, blackColor);
    XFillRectangle(dpy,buffer,gc,0,0,IMAGE_WIDTH,IMAGE_HEIGHT);
    int rect_width=(int)((double)IMAGE_WIDTH/(double)(ncols+1));
    int rect_height=(int)((double)IMAGE_HEIGHT/(double)(nrows+1));
    for (int i=1;i<=ncols;i++) {
        x1 = (int)((double)(i-1)/(double)(ncols+1)*(double)IMAGE_WIDTH);
        for (int j=1;j<=nrows;j++) {
            y1 = (int)((double)(j-1)/(double)(nrows+1)*
                (double)IMAGE_HEIGHT);
            if (grid[i][j]>0) {
                int life =grid[i][j];
                if (life>numXGrayscale-1) life=numXGrayscale-1;
                XSetForeground(dpy, gc, Xgrayscale[life].pixel);
            } else {
                XSetForeground(dpy, gc, blackColor);
            }
            XFillRectangle(dpy,buffer,gc,x1,y1,rect_width,rect_height);
         }
     }
     XSetForeground(dpy,gc,blackColor);
     XFillRectangle(dpy,buffer,gc,10,10,15,15);
     XSetForeground(dpy,gc,whiteColor);
     XDrawRectangle(dpy,buffer,gc,10,10,15,15);
     XDrawString(dpy,buffer,gc,12,23,string,2);
     
     XCopyArea(dpy, buffer, w, gc, 0, 0,
         IMAGE_WIDTH, IMAGE_HEIGHT,  0, 0);
     XFlush(dpy);
          
}
void setup_mpi(int *argc,char*** argv, mpidata* my_mpi) {
    MPI_Init(argc,argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_mpi->rank);
    MPI_Comm_size(MPI_COMM_WORLD, &my_mpi->size);
}
void finalize_mpi() {
    MPI_Finalize();
}
void random_initByTime(int rank) {
    time_t ltime;
    time(&ltime);
    srand((unsigned) ltime + 100*rank);
}
