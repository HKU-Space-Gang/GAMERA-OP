#ifndef ANALYSIS_IO_H
#define ANALYSIS_IO_H

/* Write one compact, ghost-free HDF5 snapshot per Yin/Yang patch.  All ranks
 * in the patch Cartesian communicator participate in a parallel hyperslab
 * write. */
int dump_analysis_hdf5(void);

/* Recover the last compact-analysis sequence from the shared patch files so
 * a restart cannot overwrite an earlier snapshot. */
int initialize_analysis_sequence(void);

#endif
