// -----------------------------------------------------------------
// Directions, links, and macros to give their opposites
// MPI communications assume directions from 0 to 5
#ifndef _DIRS_H
#define _DIRS_H

#define NDIMS 3       // Number of dimensions
#define NODIR -1      // Not a direction
#define XUP 0
#define YUP 1
#define TUP 2
#define TDOWN 3
#define YDOWN 4
#define XDOWN 5
#define OPP_DIR(dir) (5 - (dir))  // Opposite spacetime direction

#endif
// -----------------------------------------------------------------
