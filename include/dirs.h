// -----------------------------------------------------------------
// Directions, links, and macros to give their opposites
// MPI communications assume directions from 0 to 7
#ifndef _DIRS_H
#define _DIRS_H

#define NDIMS 4       // Number of dimensions
#define NODIR -1      // Not a direction
#define XUP 0
#define YUP 1
#define ZUP 2
#define TUP 3
#define TDOWN 4
#define ZDOWN 5
#define YDOWN 6
#define XDOWN 7
#define OPP_DIR(dir) (7 - (dir))  // Opposite spacetime direction

#endif
// -----------------------------------------------------------------
