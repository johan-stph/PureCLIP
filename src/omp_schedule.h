// ======================================================================
// PureCLIP — OpenMP schedule selection (compile-time)
// ======================================================================
// Controlled by CMake option -DPURECLIP_OMP_SCHEDULE:
//   PURECLIP_OMP_SCHEDULE=guided     — faster, fewer barriers (default)
//   PURECLIP_OMP_SCHEDULE=dynamic,1  — original PureCLIP output
// ======================================================================

#ifndef PURECLIP_OMP_SCHEDULE_H_
#define PURECLIP_OMP_SCHEDULE_H_

#ifndef PURECLIP_OMP_SCHEDULE
#error "PURECLIP_OMP_SCHEDULE must be defined via CMake (e.g. -DPURECLIP_OMP_SCHEDULE=guided)"
#endif

// Build the full schedule clause from the CMake token.
// CMake passes: -DPURECLIP_OMP_SCHEDULE=guided  or  -DPURECLIP_OMP_SCHEDULE=dynamic,1
// The comma in "dynamic,1" is fine — it's inside the preprocessor token.
#define PURECLIP_SCHEDULE_ schedule(PURECLIP_OMP_SCHEDULE)

// Ready-to-use pragma clauses
#define PURECLIP_OMP_FOR              for PURECLIP_SCHEDULE_
#define PURECLIP_OMP_FOR_NOWAIT       for PURECLIP_SCHEDULE_ nowait
#define PURECLIP_OMP_PARALLEL_FOR     parallel for PURECLIP_SCHEDULE_

#endif
