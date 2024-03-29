#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*! \file pm_periodic.c
 *  \brief routines for periodic PM-force computation
 */
/*!
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * slightly by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#include "../allvars.h"
#include "../proto.h"

#ifdef PMGRID
#ifdef BOX_PERIODIC

#ifndef USE_FFTW3
#ifdef NOTYPEPREFIX_FFTW
#include        <rfftw_mpi.h>
#else
#ifdef DOUBLEPRECISION_FFTW
#include     <drfftw_mpi.h>	/* double precision FFTW */
#else
#include     <srfftw_mpi.h>
#endif
#endif
#define cmplx_re(c) ((c).re)
#define cmplx_im(c) ((c).im)
#else 
#include "myfftw3.h"
#endif

#define  PMGRID2 (2*(PMGRID/2 + 1))

#if (PMGRID > 1024)
typedef long long large_array_offset;
#else
typedef unsigned int large_array_offset;
#endif

#define d_fftw_real fftw_real

#ifndef USE_FFTW3
static rfftwnd_mpi_plan fft_forward_plan, fft_inverse_plan;
#else 
static fftw_plan fft_forward_plan, fft_inverse_plan;
#endif


static int slab_to_task[PMGRID];
#ifndef USE_FFTW3
static int *slabs_per_task;
static int *first_slab_of_task;
#endif

#ifndef USE_FFTW3
static int slabstart_x, nslab_x, slabstart_y, nslab_y; //, smallest_slab;
static int fftsize, maxfftsize;
#else 
static ptrdiff_t *slabs_per_task;
static ptrdiff_t *first_slab_of_task;
static ptrdiff_t slabstart_x, nslab_x, slabstart_y, nslab_y; 
static ptrdiff_t fftsize, maxfftsize;
static MPI_Datatype MPI_TYPE_PTRDIFF; 
#endif

static fftw_real *rhogrid, *forcegrid, *workspace;
static d_fftw_real *d_rhogrid, *d_forcegrid, *d_workspace;

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
static fftw_real *tidal_workspace;
static fftw_real *d_tidal_workspace;
#endif


static fftw_complex *fft_of_rhogrid;


static MyFloat to_slab_fac;

void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch);
void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch);
int pm_periodic_compare_sortindex(const void *a, const void *b);

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
void pm_periodic_transposeAz(fftw_real * field, fftw_real * scratch);
void pm_periodic_transposeBz(fftw_real * field, fftw_real * scratch);
#endif

static struct part_slab_data
{
  large_array_offset globalindex;
  int partindex;
  int localindex;
} *part;

static int *part_sortindex;


/*! This routines generates the FFTW-plans to carry out the parallel FFTs
 *  later on. Some auxiliary variables are also initialized.
 */
void pm_init_periodic(void)
{
  int i, slab_to_task_local[PMGRID];
#ifdef USE_FFTW3
  double bytes_tot; bytes_tot = 0; size_t bytes;
#endif
    
  All.Asmth[0] = PM_ASMTH * All.BoxSize / PMGRID; /* note that these routines REQUIRE a uniform (BOX_LONG_X=BOX_LONG_Y=BOX_LONG_Z=1) box, so we can just use 'BoxSize' */
  All.Rcut[0] = PM_RCUT * All.Asmth[0];

#ifndef USE_FFTW3
  /* Set up the FFTW plan files. */

  fft_forward_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE | FFTW_IN_PLACE);
  fft_inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, PMGRID, PMGRID, PMGRID,
					     FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE | FFTW_IN_PLACE);

  /* Workspace out the ranges on each processor. */

  rfftwnd_mpi_local_sizes(fft_forward_plan, &nslab_x, &slabstart_x, &nslab_y, &slabstart_y, &fftsize);
#else 
  /* define MPI_TYPE_PTRDIFF */

  if (sizeof(ptrdiff_t) == sizeof(long long)) {
    MPI_TYPE_PTRDIFF = MPI_LONG_LONG; 
  } else if (sizeof(ptrdiff_t) == sizeof(long)) {
    MPI_TYPE_PTRDIFF = MPI_LONG; 
  } else if (sizeof(ptrdiff_t) == sizeof(int)) {
    MPI_TYPE_PTRDIFF = MPI_INT; 
  }

  /* get local data size and allocate */

  //fftsize = fftw_mpi_local_size_3d(PMGRID, PMGRID, PMGRID2, MPI_COMM_WORLD, &nslab_x, &slabstart_x); 
  fftsize = fftw_mpi_local_size_3d_transposed(PMGRID, PMGRID, PMGRID2, MPI_COMM_WORLD, 
	  &nslab_x, &slabstart_x, &nslab_y, &slabstart_y); 
#endif

  for(i = 0; i < PMGRID; i++)
    slab_to_task_local[i] = 0;

  for(i = 0; i < nslab_x; i++)
    slab_to_task_local[slabstart_x + i] = ThisTask;

  MPI_Allreduce(slab_to_task_local, slab_to_task, PMGRID, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

#ifndef USE_FFTW3 
  /* not used */
  /*
  MPI_Allreduce(&nslab_x, &smallest_slab, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  */

  slabs_per_task = (int *) mymalloc("slabs_per_task", NTask * sizeof(int));
  MPI_Allgather(&nslab_x, 1, MPI_INT, slabs_per_task, 1, MPI_INT, MPI_COMM_WORLD);

  first_slab_of_task = (int *) mymalloc("first_slab_of_task", NTask * sizeof(int));
  MPI_Allgather(&slabstart_x, 1, MPI_INT, first_slab_of_task, 1, MPI_INT, MPI_COMM_WORLD);

  to_slab_fac = PMGRID / All.BoxSize;

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#else 
  slabs_per_task = (ptrdiff_t *) mymalloc("slabs_per_task", NTask * sizeof(ptrdiff_t));
  MPI_Allgather(&nslab_x, 1, MPI_TYPE_PTRDIFF, slabs_per_task, 1, MPI_TYPE_PTRDIFF, MPI_COMM_WORLD);

  first_slab_of_task = (ptrdiff_t *) mymalloc("first_slab_of_task", NTask * sizeof(ptrdiff_t));
  MPI_Allgather(&slabstart_x, 1, MPI_TYPE_PTRDIFF, first_slab_of_task, 1, MPI_TYPE_PTRDIFF, MPI_COMM_WORLD);

  to_slab_fac = PMGRID / All.BoxSize;

  MPI_Allreduce(&fftsize, &maxfftsize, 1, MPI_TYPE_PTRDIFF, MPI_MAX, MPI_COMM_WORLD);

  if(!(rhogrid = (fftw_real *) mymalloc("rhogrid", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;


  if(ThisTask == 0)
    printf("Allocated %g MByte for rhogrid.\n", bytes_tot / (1024.0 * 1024.0));

  fft_of_rhogrid = (fftw_complex *) rhogrid;

  fft_forward_plan = fftw_mpi_plan_dft_r2c_3d(PMGRID, PMGRID, PMGRID, rhogrid, fft_of_rhogrid, 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_OUT); 

  fft_inverse_plan = fftw_mpi_plan_dft_c2r_3d(PMGRID, PMGRID, PMGRID, fft_of_rhogrid, rhogrid, 
	  MPI_COMM_WORLD, FFTW_ESTIMATE | FFTW_MPI_TRANSPOSED_IN); 

#endif

}


/*! This function allocates the memory neeed to compute the long-range PM
 *  force. Three fields are used, one to hold the density (and its FFT, and
 *  then the real-space potential), one to hold the force field obtained by
 *  finite differencing, and finally a workspace field, which is used both as
 *  workspace for the parallel FFT, and as buffer for the communication
 *  algorithm used in the force computation.
 */
void pm_init_periodic_allocate(void)
{
  double bytes_tot = 0;
  size_t bytes;

  /* allocate the memory to hold the FFT fields */

#ifndef USE_FFTW3
  if(!(rhogrid = (fftw_real *) mymalloc("rhogrid", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-rhogrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;
#endif

  if(!(forcegrid = (fftw_real *) mymalloc("forcegrid", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-forcegrid' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!
     (part = (struct part_slab_data *) mymalloc("part", bytes = 8 * NumPart * sizeof(struct part_slab_data))))
    {
      printf("failed to allocate memory for `part' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

  if(!(part_sortindex = (int *) mymalloc("part_sortindex", bytes = 8 * NumPart * sizeof(int))))
    {
      printf("failed to allocate memory for `part_sortindex' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
  bytes_tot += bytes;

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
  if(!(tidal_workspace = (fftw_real *) mymalloc("tidal_workspace", bytes = maxfftsize * sizeof(d_fftw_real))))
    {
      printf("failed to allocate memory for `FFT-tidal_workspace' (%g MB).\n", bytes / (1024.0 * 1024.0));
      endrun(1);
    }
#endif
  bytes_tot += bytes;


  if(ThisTask == 0)
    printf(" ..using %g MByte for periodic FFT computation. (presently allocated=%g MB)\n",
	   bytes_tot / (1024.0 * 1024.0), AllocatedBytes / (1024.0 * 1024.0));

  workspace = forcegrid;

#ifndef USE_FFTW3
  fft_of_rhogrid = (fftw_complex *) & rhogrid[0];
#endif

  d_rhogrid = (d_fftw_real *) rhogrid;
  d_forcegrid = (d_fftw_real *) forcegrid;
  d_workspace = (d_fftw_real *) workspace;
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
  d_tidal_workspace = (d_fftw_real *) tidal_workspace;
#endif
}



/*! This routine frees the space allocated for the parallel FFT algorithm.
 */
void pm_init_periodic_free(void)
{
  /* allocate the memory to hold the FFT fields */
#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
  myfree(tidal_workspace);
#endif
  myfree(part_sortindex);
  myfree(part);
  myfree(forcegrid);
#ifndef USE_FFTW3
  myfree(rhogrid);
#endif
}

#ifdef ALT_QSORT
#define KEY_TYPE int
#define KEY_BASE_TYPE large_array_offset
#define KEY_GETVAL(pk) (part[*(pk)].globalindex)
#define KEY_COPY(pk1,pk2)       \
{                               \
  *(pk2) = *(pk1);      \
}
#define QSORT qsort_pm_periodic
#include "../system/myqsort.h"
#endif

/*! Calculates the long-range periodic force given the particle positions
 *  using the PM method.  The force is Gaussian filtered with Asmth, given in
 *  mesh-cell units. We carry out a CIC charge assignment, and compute the
 *  potenial by Fourier transform methods. The potential is finite differenced
 *  using a 4-point finite differencing formula, and the forces are
 *  interpolated tri-linearly to the particle positions. The CIC kernel is
 *  deconvolved. Note that the particle distribution is not in the slab
 *  decomposition that is used for the FFT. Instead, overlapping patches
 *  between local domains and FFT slabs are communicated as needed.
 *
 *  For mode=0, normal force calculation, mode=1, only PS calculation.
 */
void pmforce_periodic(int mode, int *typelist)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, acc_dim;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, yl, zl, yr, zr, yll, zll, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  MyDouble pp[3], *pos;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

#ifdef DM_SCALARFIELD_SCREENING
  int phase;
  double kscreening2;

  kscreening2 = pow(All.BoxSize / All.ScalarScreeningLength / (2 * M_PI), 2);
#endif
 
  PRINT_STATUS("Starting periodic PM calculation.  (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */

  pm_init_periodic_allocate();

#ifdef DM_SCALARFIELD_SCREENING
  for(phase = 0; phase < 2; phase++)
    {
#endif

      /* determine the cells each particles accesses */
      for(i = 0, num_on_grid = 0; i < NumPart; i++)
	{
	  if(mode)
	    {
	      /* only power spectrum calculation */
	      if(typelist[P[i].Type] == 0)
		continue;
	    }

#ifdef DM_SCALARFIELD_SCREENING
	  if(phase == 1)
	    if(P[i].Type == 0)	/* don't bin baryonic mass in this phase */
	      continue;
#endif

        /* possible bugfix: Y.Feng:  (was if(mode)) */
	  if(mode > -1)
	    {
	      /* make sure that particles are properly box-wrapped */
	      for(j = 0; j < 3; j++)
		{
		  pp[j] = P[i].Pos[j]; 
		  pp[j] = WRAP_POSITION_UNIFORM_BOX(pp[j]);
		}
	      pos = pp;
	    }
	  else
	    pos = P[i].Pos;

	  slab_x = (int) (to_slab_fac * pos[0]);
	  slab_y = (int) (to_slab_fac * pos[1]);
	  slab_z = (int) (to_slab_fac * pos[2]);

        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

	  for(xx = 0; xx < 2; xx++)
	    for(yy = 0; yy < 2; yy++)
	      for(zz = 0; zz < 2; zz++)
		{
		  slab_xx = slab_x + xx;
		  slab_yy = slab_y + yy;
		  slab_zz = slab_z + zz;

		  if(slab_xx >= PMGRID)
		    slab_xx -= PMGRID;
		  if(slab_yy >= PMGRID)
		    slab_yy -= PMGRID;
		  if(slab_zz >= PMGRID)
		    slab_zz -= PMGRID;

		  offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

		  part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
		  part[num_on_grid].globalindex = offset;
		  part_sortindex[num_on_grid] = num_on_grid;
		  num_on_grid++;
		}
	}
      /* note: num_on_grid will be  8 times larger than the particle number,
         but num_field_points will generally be much smaller */

      /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
      mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
      qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

      /* determine the number of unique field points */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  num_field_points++;
	}

      /* allocate the local field */
      localfield_globalindex =
	(large_array_offset *) mymalloc("localfield_globalindex",
					num_field_points * sizeof(large_array_offset));
      localfield_d_data =
	(d_fftw_real *) mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
      localfield_data = (fftw_real *) localfield_d_data;
      localfield_first = (int *) mymalloc("localfield_first", NTask * sizeof(int));
      localfield_count = (int *) mymalloc("localfield_count", NTask * sizeof(int));
      localfield_offset = (int *) mymalloc("localfield_offset", NTask * sizeof(int));
      localfield_togo = (int *) mymalloc("localfield_togo", NTask * NTask * sizeof(int));

      for(i = 0; i < NTask; i++)
	{
	  localfield_first[i] = 0;
	  localfield_count[i] = 0;
	}

      /* establish the cross link between the part[] array and the local list of
         mesh points. Also, count on which CPU how many of the needed field points are stored */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	      num_field_points++;

	  part[part_sortindex[i]].localindex = num_field_points;

	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

	  slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
	  task = slab_to_task[slab];
	  if(localfield_count[task] == 0)
	    localfield_first[task] = num_field_points;
	  localfield_count[task]++;
	}
      num_field_points++;

      for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
	localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

      if(mode != 0)
	{
	  foldonitself(typelist);
	  powerspec(0, typelist);
	}

      /* now bin the local particle data onto the mesh list */

      for(i = 0; i < num_field_points; i++)
	localfield_d_data[i] = 0;

      for(i = 0; i < num_on_grid; i += 8)
	{
	  pindex = (part[i].partindex >> 3);
        if(P[pindex].Mass<=0) continue;

        /* possible bugfix: Y.Feng:  (was if(mode)) */
	  if(mode > -1)
	    {
	      /* make sure that particles are properly box-wrapped */
	      for(j = 0; j < 3; j++)
		{
		  pp[j] = P[pindex].Pos[j];
            pp[j] = WRAP_POSITION_UNIFORM_BOX(pp[j]);
		}
	      pos = pp;
	    }
	  else
	    pos = P[pindex].Pos;

	  slab_x = (int) (to_slab_fac * pos[0]);
	  slab_y = (int) (to_slab_fac * pos[1]);
	  slab_z = (int) (to_slab_fac * pos[2]);
        /* possible bugfix: Y.Feng:  */
        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

	  dx = to_slab_fac * pos[0] - slab_x;
	  dy = to_slab_fac * pos[1] - slab_y;
	  dz = to_slab_fac * pos[2] - slab_z;

	  localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
	  localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
	}

      /* clear local FFT-mesh density field */
      for(i = 0; i < fftsize; i++)
	d_rhogrid[i] = 0;

      /* exchange data and add contributions to the local mesh-path */

      MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(level > 0)
		{
		  import_d_data =
		    (d_fftw_real *) mymalloc("import_d_data", localfield_togo[recvTask * NTask + ThisTask] *
					     sizeof(d_fftw_real));
		  import_globalindex =
		    (large_array_offset *) mymalloc("import_globalindex",
						    localfield_togo[recvTask * NTask +
								    ThisTask] * sizeof(large_array_offset));

		  if(localfield_togo[sendTask * NTask + recvTask] > 0
		     || localfield_togo[recvTask * NTask + sendTask] > 0)
		    {
		      MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_A, import_d_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		    }
		}
	      else
		{
		  import_d_data = localfield_d_data + localfield_offset[ThisTask];
		  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		}

	      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		{
		  /* determine offset in local FFT slab */
		  offset =
		    import_globalindex[i] -
		    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

		  d_rhogrid[offset] += import_d_data[i];
		}

	      if(level > 0)
		{
		  myfree(import_globalindex);
		  myfree(import_d_data);
		}
	    }
	}

      /* Do the FFT of the density field */

      report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC");

#ifndef USE_FFTW3
      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
      fftw_execute(fft_forward_plan); 
#endif

      if(mode != 0)
	{
	  powerspec(1, typelist);
	}

      if(mode == 0)		/* only carry out this part for the ordinary force calculation */
	{
	  /* multiply with Green's function for the potential */

	  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
	    for(x = 0; x < PMGRID; x++)
	      for(z = 0; z < PMGRID / 2 + 1; z++)
		{
		  if(x > PMGRID / 2)
		    kx = x - PMGRID;
		  else
		    kx = x;
		  if(y > PMGRID / 2)
		    ky = y - PMGRID;
		  else
		    ky = y;
		  if(z > PMGRID / 2)
		    kz = z - PMGRID;
		  else
		    kz = z;

		  k2 = kx * kx + ky * ky + kz * kz;

		  if(k2 > 0)
		    {
#ifdef DM_SCALARFIELD_SCREENING
		      if(phase == 1)
			smth = -All.ScalarBeta * exp(-k2 * asmth2) / (k2 + kscreening2);
		      else
#endif
			smth = -exp(-k2 * asmth2) / k2;

		      /* do deconvolution */

		      fx = fy = fz = 1;
		      if(kx != 0)
			{
			  fx = (M_PI * kx) / PMGRID;
			  fx = sin(fx) / fx;
			}
		      if(ky != 0)
			{
			  fy = (M_PI * ky) / PMGRID;
			  fy = sin(fy) / fy;
			}
		      if(kz != 0)
			{
			  fz = (M_PI * kz) / PMGRID;
			  fz = sin(fz) / fz;
			}
		      ff = 1 / (fx * fy * fz);
		      smth *= ff * ff * ff * ff;

		      /* end deconvolution */

		      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
		      cmplx_re(fft_of_rhogrid[ip]) *= smth;
		      cmplx_im(fft_of_rhogrid[ip]) *= smth;
		    }
		}

	  if(slabstart_y == 0)
	    cmplx_re(fft_of_rhogrid[0]) = cmplx_im(fft_of_rhogrid[0]) = 0.0;

	  /* Do the inverse FFT to get the potential */

#ifndef USE_FFTW3
	  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
	  fftw_execute(fft_inverse_plan);  
#endif

	  /* Now rhogrid holds the potential */

#ifdef EVALPOTENTIAL		/* now read out the potential */
	  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ level;

	      if(recvTask < NTask)
		{
		  if(level > 0)
		    {
		      import_data =
			(fftw_real *) mymalloc("import_data",
					       localfield_togo[recvTask * NTask +
							       ThisTask] * sizeof(fftw_real));
		      import_globalindex =
			(large_array_offset *) mymalloc("import_globalindex",
							localfield_togo[recvTask * NTask +
									ThisTask] *
							sizeof(large_array_offset));

		      if(localfield_togo[sendTask * NTask + recvTask] > 0
			 || localfield_togo[recvTask * NTask + sendTask] > 0)
			{
			  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				       localfield_togo[sendTask * NTask +
						       recvTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_PERIODIC_C, import_globalindex,
				       localfield_togo[recvTask * NTask +
						       sendTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
			}
		    }
		  else
		    {
		      import_data = localfield_data + localfield_offset[ThisTask];
		      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		    }

		  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		    {
		      offset =
			import_globalindex[i] -
			first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
		      import_data[i] = rhogrid[offset];
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(import_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_A,
				   localfield_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		      myfree(import_globalindex);
		      myfree(import_data);
		    }
		}
	    }

	  /* read out the potential values, which all have been assembled in localfield_data */

	  double pot;

	  for(i = 0, j = 0; i < NumPart; i++)
	    {
	      while(j < num_on_grid && (part[j].partindex >> 3) != i)
              j++;

            /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
            /* make sure that particles are properly box-wrapped */
            for(xx = 0; xx < 3; xx++)
            {
                pp[xx] = P[i].Pos[xx];
                pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
            }
            slab_x = (int) (to_slab_fac * pp[0]);
            slab_y = (int) (to_slab_fac * pp[1]);
            slab_z = (int) (to_slab_fac * pp[2]);
            dx = to_slab_fac * pp[0] - slab_x;
            dy = to_slab_fac * pp[1] - slab_y;
            dz = to_slab_fac * pp[2] - slab_z;

            pot =
            +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
            + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
            + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
            + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
            + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
            + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
            + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
            + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	      P[i].PM_Potential += pot * fac * (2 * All.BoxSize / PMGRID);
	      /* compensate the finite differencing factor */ ;
	    }

#endif


	  /* get the force components by finite differencing the potential for each dimension,
	     and send back the results to the right CPUs */

	  for(dim = 2; dim >= 0; dim--)	/* Calculate each component of the force. */
	    {			/* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */
	      if(dim == 0)
		pm_periodic_transposeA(rhogrid, forcegrid);	/* compute the transpose of the potential field */

	      for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
		for(y = 0; y < PMGRID; y++)
		  for(z = 0; z < PMGRID; z++)
		    {
		      x = xx - slabstart_x;

		      yrr = yll = yr = yl = y;
		      zrr = zll = zr = zl = z;

		      switch (dim)
			{
			case 0:	/* note: for the x-direction, we difference the transposed direction (y) */
			case 1:
			  yr = y + 1;
			  yl = y - 1;
			  yrr = y + 2;
			  yll = y - 2;
			  if(yr >= PMGRID)
			    yr -= PMGRID;
			  if(yrr >= PMGRID)
			    yrr -= PMGRID;
			  if(yl < 0)
			    yl += PMGRID;
			  if(yll < 0)
			    yll += PMGRID;
			  break;
			case 2:
			  zr = z + 1;
			  zl = z - 1;
			  zrr = z + 2;
			  zll = z - 2;
			  if(zr >= PMGRID)
			    zr -= PMGRID;
			  if(zrr >= PMGRID)
			    zrr -= PMGRID;
			  if(zl < 0)
			    zl += PMGRID;
			  if(zll < 0)
			    zll += PMGRID;
			  break;
			}

		      if(dim == 0)
			{
			  forcegrid[PMGRID * (x + y * nslab_x) + z]
			    =
			    fac * ((4.0 / 3) *
				   (rhogrid[PMGRID * (x + yl * nslab_x) + zl] -
				    rhogrid[PMGRID * (x + yr * nslab_x) + zr]) -
				   (1.0 / 6) * (rhogrid[PMGRID * (x + yll * nslab_x) + zll] -
						rhogrid[PMGRID * (x + yrr * nslab_x) + zrr]));
			}
		      else
			forcegrid[PMGRID2 * (PMGRID * x + y) + z]
			  =
			  fac * ((4.0 / 3) *
				 (rhogrid[PMGRID2 * (PMGRID * x + yl) + zl] -
				  rhogrid[PMGRID2 * (PMGRID * x + yr) + zr]) -
				 (1.0 / 6) * (rhogrid[PMGRID2 * (PMGRID * x + yll) + zll] -
					      rhogrid[PMGRID2 * (PMGRID * x + yrr) + zrr]));
		    }

	      if(dim == 0)
		pm_periodic_transposeB(forcegrid, rhogrid);	/* compute the transpose of the potential field */

	      /* send the force components to the right processors */

	      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
		{
		  sendTask = ThisTask;
		  recvTask = ThisTask ^ level;

		  if(recvTask < NTask)
		    {
		      if(level > 0)
			{
			  import_data =
			    (fftw_real *) mymalloc("import_data",
						   localfield_togo[recvTask * NTask +
								   ThisTask] * sizeof(fftw_real));
			  import_globalindex =
			    (large_array_offset *) mymalloc("import_globalindex",
							    localfield_togo[recvTask * NTask +
									    ThisTask] *
							    sizeof(large_array_offset));

			  if(localfield_togo[sendTask * NTask + recvTask] > 0
			     || localfield_togo[recvTask * NTask + sendTask] > 0)
			    {
			      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
					   localfield_togo[sendTask * NTask +
							   recvTask] * sizeof(large_array_offset), MPI_BYTE,
					   recvTask, TAG_PERIODIC_C, import_globalindex,
					   localfield_togo[recvTask * NTask +
							   sendTask] * sizeof(large_array_offset), MPI_BYTE,
					   recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
			    }
			}
		      else
			{
			  import_data = localfield_data + localfield_offset[ThisTask];
			  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
			}

		      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
			{
			  /* determine offset in local FFT slab */
			  offset =
			    import_globalindex[i] -
			    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
			  import_data[i] = forcegrid[offset];
			}

		      if(level > 0)
			{
			  MPI_Sendrecv(import_data,
				       localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real),
				       MPI_BYTE, recvTask, TAG_PERIODIC_A,
				       localfield_data + localfield_offset[recvTask],
				       localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real),
				       MPI_BYTE, recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

			  myfree(import_globalindex);
			  myfree(import_data);
			}
		    }
		}

	      /* read out the forces, which all have been assembled in localfield_data */

	      for(i = 0, j = 0; i < NumPart; i++)
		{
#ifdef DM_SCALARFIELD_SCREENING
		  if(phase == 1)
		    if(P[i].Type == 0)	/* baryons don't get an extra scalar force */
		      continue;
#endif
            while(j < num_on_grid && (part[j].partindex >> 3) != i)
                j++;

            /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
            /* make sure that particles are properly box-wrapped */
            for(xx = 0; xx < 3; xx++)
            {
                pp[xx] = P[i].Pos[xx];
                pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
            }
            slab_x = (int) (to_slab_fac * pp[0]);
            slab_y = (int) (to_slab_fac * pp[1]);
            slab_z = (int) (to_slab_fac * pp[2]);
            dx = to_slab_fac * pp[0] - slab_x;
            dy = to_slab_fac * pp[1] - slab_y;
            dz = to_slab_fac * pp[2] - slab_z;

            acc_dim =
		    +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
		    + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
		    + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
		    + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
		    + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
		    + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
		    + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
		    + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;
            
		  P[i].GravPM[dim] += acc_dim;
		}

	    }			/* end of if(mode==0) block */

	}

      /* free locallist */
      myfree(localfield_togo);
      myfree(localfield_offset);
      myfree(localfield_count);
      myfree(localfield_first);
      myfree(localfield_d_data);
      myfree(localfield_globalindex);
#ifdef DM_SCALARFIELD_SCREENING
    }
#endif

  pm_init_periodic_free();
}


/*! Calculates the long-range potential using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 */
void pmpotential_periodic(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, pot;
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MyDouble pp[3];
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

  PRINT_STATUS("Starting periodic PM-potential calculation.  (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */

  pm_init_periodic_allocate();


  /* determine the cells each particles accesses */
  for(i = 0, num_on_grid = 0; i < NumPart; i++)
    {
        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

      for(xx = 0; xx < 2; xx++)
	for(yy = 0; yy < 2; yy++)
	  for(zz = 0; zz < 2; zz++)
	    {
	      slab_xx = slab_x + xx;
	      slab_yy = slab_y + yy;
	      slab_zz = slab_z + zz;

	      if(slab_xx >= PMGRID)
		slab_xx -= PMGRID;
	      if(slab_yy >= PMGRID)
		slab_yy -= PMGRID;
	      if(slab_zz >= PMGRID)
		slab_zz -= PMGRID;

	      offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

	      part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
	      part[num_on_grid].globalindex = offset;
	      part_sortindex[num_on_grid] = num_on_grid;
	      num_on_grid++;
	    }
    }

  /* note: num_on_grid will be  8 times larger than the particle number,
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
  qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

  /* determine the number of unique field points */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex =
    (large_array_offset *) mymalloc("localfield_globalindex", num_field_points * sizeof(large_array_offset));
  localfield_d_data = (d_fftw_real *) mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
  localfield_data = (fftw_real *) localfield_d_data;
  localfield_first = (int *) mymalloc("localfield_first", NTask * sizeof(int));
  localfield_count = (int *) mymalloc("localfield_count", NTask * sizeof(int));
  localfield_offset = (int *) mymalloc("localfield_offset", NTask * sizeof(int));
  localfield_togo = (int *) mymalloc("localfield_togo", NTask * NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_count[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of
     mesh points. Also, count on which CPU how many of the needed field points are stored */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	  num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

      slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
      task = slab_to_task[slab];
      if(localfield_count[task] == 0)
	localfield_first[task] = num_field_points;
      localfield_count[task]++;
    }
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

  /* now bin the local particle data onto the mesh list */

  for(i = 0; i < num_field_points; i++)
    localfield_d_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    {
      pindex = (part[i].partindex >> 3);
        if(P[pindex].Mass<=0) continue;

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[pindex].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        dx = to_slab_fac * pp[0] - slab_x;
        dy = to_slab_fac * pp[1] - slab_y;
        dz = to_slab_fac * pp[2] - slab_z;

      localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
      localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
      localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
      localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
      localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
    }

  /* clear local FFT-mesh density field */
  for(i = 0; i < fftsize; i++)
    d_rhogrid[i] = 0;

  /* exchange data and add contributions to the local mesh-path */

  MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_d_data =
		(d_fftw_real *) mymalloc("import_d_data",
					 localfield_togo[recvTask * NTask + ThisTask] * sizeof(d_fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc("import_globalindex",
						localfield_togo[recvTask * NTask +
								ThisTask] * sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A,
			       import_d_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_B, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_d_data = localfield_d_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

	      d_rhogrid[offset] += import_d_data[i];
	    }

	  if(level > 0)
	    {
	      myfree(import_globalindex);
	      myfree(import_d_data);
	    }
	}
    }

  report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC_POTENTIAL");

  /* Do the FFT of the density field */
#ifndef USE_FFTW3
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
  fftw_execute(fft_forward_plan); 
#endif

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2 * fac;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;

	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
	      cmplx_re(fft_of_rhogrid[ip]) *= smth;
	      cmplx_im(fft_of_rhogrid[ip]) *= smth;
	    }
	}

  if(slabstart_y == 0)
    cmplx_re(fft_of_rhogrid[0]) = cmplx_im(fft_of_rhogrid[0]) = 0.0;

  /* Do the inverse FFT to get the potential */

#ifndef USE_FFTW3
  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
  fftw_execute(fft_inverse_plan); 
#endif

  /* Now rhogrid holds the potential */


  /* now read out the potential */

  /* send the force components to the right processors */

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_data =
		(fftw_real *) mymalloc("import_data",
				       localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc("import_globalindex",
						localfield_togo[recvTask * NTask +
								ThisTask] * sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_C, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_data = localfield_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
	      import_data[i] = rhogrid[offset];
	    }

	  if(level > 0)
	    {
	      MPI_Sendrecv(import_data,
			   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_PERIODIC_A,
			   localfield_data + localfield_offset[recvTask],
			   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

	      myfree(import_globalindex);
	      myfree(import_data);
	    }
	}
    }

  /* read out the potential values, which all have been assembled in localfield_data */

  for(i = 0, j = 0; i < NumPart; i++)
    {
        while(j < num_on_grid && (part[j].partindex >> 3) != i)
            j++;

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        dx = to_slab_fac * pp[0] - slab_x;
        dy = to_slab_fac * pp[1] - slab_y;
        dz = to_slab_fac * pp[2] - slab_z;

        pot =
        +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
        + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
        + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
        + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
        + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
        + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
        + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
        + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

#if defined(EVALPOTENTIAL) || defined(COMPUTE_POTENTIAL_ENERGY) || defined(OUTPUT_POTENTIAL)
      P[i].Potential += pot;
#endif
    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm_init_periodic_free();
  PRINT_STATUS(" ..done PM-Potential");
}



int pm_periodic_compare_sortindex(const void *a, const void *b)
{
  if(part[*(int *) a].globalindex < part[*(int *) b].globalindex) {return -1;}
  if(part[*(int *) a].globalindex > part[*(int *) b].globalindex) {return +1;}
  return 0;
}

static void msort_pmperiodic_with_tmp(int *b, size_t n, int *t)
{
  int *tmp;
  int *b1, *b2;
  size_t n1, n2;

  if(n <= 1)
    return;

  n1 = n / 2;
  n2 = n - n1;
  b1 = b;
  b2 = b + n1;

  msort_pmperiodic_with_tmp(b1, n1, t);
  msort_pmperiodic_with_tmp(b2, n2, t);

  tmp = t;

  while(n1 > 0 && n2 > 0)
    {
      if(part[*b1].globalindex <= part[*b2].globalindex)
	{
	  --n1;
	  *tmp++ = *b1++;
	}
      else
	{
	  --n2;
	  *tmp++ = *b2++;
	}
    }

  if(n1 > 0)
    memcpy(tmp, b1, n1 * sizeof(int));

  memcpy(b, t, (n - n2) * sizeof(int));
}

void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp) (const void *, const void *))
{
  const size_t size = n * s;

  int *tmp = (int *) mymalloc("int *tmp", size);

  msort_pmperiodic_with_tmp((int *) b, n, tmp);

  myfree(tmp);
}

void pm_periodic_transposeA(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < PMGRID; z++)
	  {
	    scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
			      x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z] =
	      field[PMGRID2 * (PMGRID * x + y) + z];
	  }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif
}



void pm_periodic_transposeB(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }


  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = first_slab_of_task[task]; y < first_slab_of_task[task] + slabs_per_task[task]; y++)
	for(z = 0; z < PMGRID; z++)
	  {
	    field[PMGRID2 * (PMGRID * x + y) + z] =
	      scratch[PMGRID * (first_slab_of_task[task] * nslab_x +
				x * slabs_per_task[task] + (y - first_slab_of_task[task])) + z];
	  }

}

#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
void pm_periodic_transposeAz(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = 0; y < PMGRID; y++)
	for(z = first_slab_of_task[task]; z < first_slab_of_task[task] + slabs_per_task[task]; z++)
	  {
	    scratch[nslab_x * (first_slab_of_task[task] * PMGRID +
			       x * PMGRID + y) + (z - first_slab_of_task[task])] =
	      field[PMGRID2 * (PMGRID * x + y) + z];
	  }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }

  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);
#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif
}



void pm_periodic_transposeBz(fftw_real * field, fftw_real * scratch)
{
  int x, y, z, task;

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  MPI_Request *requests;
  int nrequests = 0;

  requests = (MPI_Request *) mymalloc("requests", 2 * NTask * sizeof(MPI_Request));

  for(task = 0; task < NTask; task++)
    {
      MPI_Isend(field + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);

      MPI_Irecv(scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, &requests[nrequests++]);
    }


  MPI_Waitall(nrequests, requests, MPI_STATUSES_IGNORE);
  myfree(requests);

#else
  int ngrp;

  for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
    {
      task = ThisTask ^ ngrp;

      if(task < NTask)
	{
	  MPI_Sendrecv(field + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY,
		       scratch + PMGRID * first_slab_of_task[task] * nslab_x,
		       PMGRID * nslab_x * slabs_per_task[task] * sizeof(fftw_real),
		       MPI_BYTE, task, TAG_KEY, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }
#endif

  for(task = 0; task < NTask; task++)
    for(x = 0; x < nslab_x; x++)
      for(y = 0; y < PMGRID; y++)
	for(z = first_slab_of_task[task]; z < first_slab_of_task[task] + slabs_per_task[task]; z++)
	  {
	    field[PMGRID2 * (PMGRID * x + y) + z] =
	      scratch[nslab_x * (first_slab_of_task[task] * PMGRID +
				 x * PMGRID + y) + (z - first_slab_of_task[task])];
	  }

}
#endif


#endif
#endif







#ifdef PMGRID
#ifdef BOX_PERIODIC



#ifdef COMPUTE_TIDAL_TENSOR_IN_GRAVTREE
/*! Calculates the long-range tidal field using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 *  Like the forces the derivatives are calculated by finite differences
 *  of the potential on the grid.
 */

void pmtidaltensor_periodic_diff(void)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, tidal_dim;
  MyDouble pp[3];
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, yl, zl, yr, zr, yll, zll, yrr, zrr, ip, dim;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;

#ifdef DM_SCALARFIELD_SCREENING
  int phase;
  double kscreening2;

  kscreening2 = pow(All.BoxSize / All.ScalarScreeningLength / (2 * M_PI), 2);
#endif
 
  PRINT_STATUS("Starting periodic PM-TIDAL calculation.  (presently allocated=%g MB)", AllocatedBytes / (1024.0 * 1024.0));
  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */
  fac *= 1 / (2 * All.BoxSize / PMGRID);	/* for finite differencing */


  pm_init_periodic_allocate();

#ifdef DM_SCALARFIELD_SCREENING
  for(phase = 0; phase < 2; phase++)
    {
#endif

      /* determine the cells each particles accesses */
      for(i = 0, num_on_grid = 0; i < NumPart; i++)
	{
#ifdef DM_SCALARFIELD_SCREENING
	  if(phase == 1)
	    if(P[i].Type == 0)	/* don't bin baryonic mass in this phase */
	      continue;
#endif

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID;

	  for(xx = 0; xx < 2; xx++)
	    for(yy = 0; yy < 2; yy++)
	      for(zz = 0; zz < 2; zz++)
		{
		  slab_xx = slab_x + xx;
		  slab_yy = slab_y + yy;
		  slab_zz = slab_z + zz;

		  if(slab_xx >= PMGRID)
		    slab_xx -= PMGRID;
		  if(slab_yy >= PMGRID)
		    slab_yy -= PMGRID;
		  if(slab_zz >= PMGRID)
		    slab_zz -= PMGRID;

		  offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz;

		  part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz;
		  part[num_on_grid].globalindex = offset;
		  part_sortindex[num_on_grid] = num_on_grid;
		  num_on_grid++;
		}
	}
      /* note: num_on_grid will be  8 times larger than the particle number,
         but num_field_points will generally be much smaller */

      /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
      mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
      qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

      /* determine the number of unique field points */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  num_field_points++;
	}

      /* allocate the local field */
      localfield_globalindex =
	(large_array_offset *) mymalloc("localfield_globalindex",
					num_field_points * sizeof(large_array_offset));
      localfield_d_data =
	(d_fftw_real *) mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
      localfield_data = (fftw_real *) localfield_d_data;
      localfield_first = (int *) mymalloc("localfield_first", NTask * sizeof(int));
      localfield_count = (int *) mymalloc("localfield_count", NTask * sizeof(int));
      localfield_offset = (int *) mymalloc("localfield_offset", NTask * sizeof(int));
      localfield_togo = (int *) mymalloc("localfield_togo", NTask * NTask * sizeof(int));

      for(i = 0; i < NTask; i++)
	{
	  localfield_first[i] = 0;
	  localfield_count[i] = 0;
	}

      /* establish the cross link between the part[] array and the local list of
         mesh points. Also, count on which CPU how many of the needed field points are stored */
      for(i = 0, num_field_points = 0; i < num_on_grid; i++)
	{
	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	      num_field_points++;

	  part[part_sortindex[i]].localindex = num_field_points;

	  if(i > 0)
	    if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	      continue;

	  localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

	  slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
	  task = slab_to_task[slab];
	  if(localfield_count[task] == 0)
	    localfield_first[task] = num_field_points;
	  localfield_count[task]++;
	}
      num_field_points++;

      for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
	localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];


      /* now bin the local particle data onto the mesh list */

      for(i = 0; i < num_field_points; i++)
	localfield_d_data[i] = 0;

      for(i = 0; i < num_on_grid; i += 8) 
	{
	  pindex = (part[i].partindex >> 3); 
	  if(P[pindex].Mass<=0) continue; 
	  
	  /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */ 
	  /* make sure that particles are properly box-wrapped */ 
	  for(xx = 0; xx < 3; xx++) 
	  { 
	      pp[xx] = P[pindex].Pos[xx]; 
	      pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]); 
	  } 
	  slab_x = (int) (to_slab_fac * pp[0]); 
	  slab_y = (int) (to_slab_fac * pp[1]); 
	  slab_z = (int) (to_slab_fac * pp[2]); 
	  dx = to_slab_fac * pp[0] - slab_x; 
	  dy = to_slab_fac * pp[1] - slab_y; 
	  dz = to_slab_fac * pp[2] - slab_z;

	  localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
	  localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	  localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
	  localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
	  localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
	} 
      
      /* clear local FFT-mesh density field */
        for(i = 0; i < fftsize; i++) {d_rhogrid[i] = 0;}

      /* exchange data and add contributions to the local mesh-path */

      MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(level > 0)
		{
		  import_d_data =
		    (d_fftw_real *) mymalloc("import_d_data", localfield_togo[recvTask * NTask + ThisTask] *
					     sizeof(d_fftw_real));
		  import_globalindex =
		    (large_array_offset *) mymalloc("import_globalindex",
						    localfield_togo[recvTask * NTask +
								    ThisTask] * sizeof(large_array_offset));

		  if(localfield_togo[sendTask * NTask + recvTask] > 0
		     || localfield_togo[recvTask * NTask + sendTask] > 0)
		    {
		      MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_A, import_d_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real),
				   MPI_BYTE, recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		    }
		}
	      else
		{
		  import_d_data = localfield_d_data + localfield_offset[ThisTask];
		  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		}

	      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		{
		  /* determine offset in local FFT slab */
		  offset =
		    import_globalindex[i] -
		    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

		  d_rhogrid[offset] += import_d_data[i];
		}

	      if(level > 0)
		{
		  myfree(import_globalindex);
		  myfree(import_d_data);
		}
	    }
	}

      /* Do the FFT of the density field */

      report_memory_usage(&HighMark_pmperiodic, "PM_PERIODIC");

#ifndef USE_FFTW3
      rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
      fftw_execute(fft_forward_plan); 
#endif

      /* multiply with Green's function for the potential */

      for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
	for(x = 0; x < PMGRID; x++)
	  for(z = 0; z < PMGRID / 2 + 1; z++)
	    {
	      if(x > PMGRID / 2)
		kx = x - PMGRID;
	      else
		kx = x;
	      if(y > PMGRID / 2)
		ky = y - PMGRID;
	      else
		ky = y;
	      if(z > PMGRID / 2)
		kz = z - PMGRID;
	      else
		kz = z;

	      k2 = kx * kx + ky * ky + kz * kz;

	      if(k2 > 0)
		{
#ifdef DM_SCALARFIELD_SCREENING
		  if(phase == 1)
		    smth = -All.ScalarBeta * exp(-k2 * asmth2) / (k2 + kscreening2);
		  else
#endif
		    smth = -exp(-k2 * asmth2) / k2;

		  /* do deconvolution */

		  fx = fy = fz = 1;
		  if(kx != 0)
		    {
		      fx = (M_PI * kx) / PMGRID;
		      fx = sin(fx) / fx;
		    }
		  if(ky != 0)
		    {
		      fy = (M_PI * ky) / PMGRID;
		      fy = sin(fy) / fy;
		    }
		  if(kz != 0)
		    {
		      fz = (M_PI * kz) / PMGRID;
		      fz = sin(fz) / fz;
		    }
		  ff = 1 / (fx * fy * fz);
		  smth *= ff * ff * ff * ff;

		  /* end deconvolution */

		  ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;
		  cmplx_re(fft_of_rhogrid[ip]) *= smth;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth;
		}
	    }

      if(slabstart_y == 0)
	cmplx_re(fft_of_rhogrid[0]) = cmplx_im(fft_of_rhogrid[0]) = 0.0;

      /* Do the inverse FFT to get the potential */

#ifndef USE_FFTW3
      rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
      fftw_execute(fft_inverse_plan); 
#endif

      /* Now rhogrid holds the potential */


#ifdef EVALPOTENTIAL		/* now read out the potential */

      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(level > 0)
		{
		  import_data =
		    (fftw_real *) mymalloc("import_data",
					   localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
		  import_globalindex =
		    (large_array_offset *) mymalloc("import_globalindex",
						    localfield_togo[recvTask * NTask +
								    ThisTask] * sizeof(large_array_offset));

		  if(localfield_togo[sendTask * NTask + recvTask] > 0
		     || localfield_togo[recvTask * NTask + sendTask] > 0)
		    {
		      MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_C, import_globalindex,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
				   MPI_BYTE, recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
		    }
		}
	      else
		{
		  import_data = localfield_data + localfield_offset[ThisTask];
		  import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		}

	      for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		{
		  offset =
		    import_globalindex[i] -
		    first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
		  import_data[i] = rhogrid[offset];
		}

	      if(level > 0)
		{
		  MPI_Sendrecv(import_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A,
			       localfield_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		  myfree(import_globalindex);
		  myfree(import_data);
		}
	    }
	}

      /* read out the potential values, which all have been assembled in localfield_data */

      double pot;

      for(i = 0, j = 0; i < NumPart; i++)
	{
	  while(j < num_on_grid && (part[j].partindex >> 3) != i)
	    j++;

	  /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */ 
	  /* make sure that particles are properly box-wrapped */ 
	  for(xx = 0; xx < 3; xx++) 
	  { 
	      pp[xx] = P[i].Pos[xx]; 
	      pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]); 
	  } 
	  slab_x = (int) (to_slab_fac * pp[0]); 
	  slab_y = (int) (to_slab_fac * pp[1]); 
	  slab_z = (int) (to_slab_fac * pp[2]); 
	  dx = to_slab_fac * pp[0] - slab_x; dy = to_slab_fac * pp[1] - slab_y; 
	  dz = to_slab_fac * pp[2] - slab_z; 
	
	  pot = 
	      +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) 
	      + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz 
	      + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz) 
	      + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz 
	      + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz) 
	      + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz 
	      + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz) 
	      + localfield_data[part[j + 7].localindex] * (dx) * dy * dz; 
	  
	  P[i].PM_Potential += pot * fac * (2 * All.BoxSize / PMGRID);
	  /* compensate the finite differencing factor */ ; 
	}
#endif




      /* get the force components by finite differencing the potential for each dimension,
         and send back the results to the right CPUs */

      for(dim = 5; dim >= 0; dim--)	/* Calculate each component of the force. */
	{			/* we do the x component last, because for differencing the potential in the x-direction, we need to contruct the transpose */

	  /* compute the transpose of the potential field */
	  if(dim == 2)
	    pm_periodic_transposeA(rhogrid, tidal_workspace);

	  if(dim == 1)
	    {
	      pm_periodic_transposeB(rhogrid, tidal_workspace);
	      pm_periodic_transposeAz(rhogrid, tidal_workspace);
	    }

	  if(dim == 0)
	    {
	      pm_periodic_transposeBz(rhogrid, tidal_workspace);
	      pm_periodic_transposeA(rhogrid, tidal_workspace);
	    }

	  for(xx = slabstart_x; xx < (slabstart_x + nslab_x); xx++)
	    for(y = 0; y < PMGRID; y++)
	      for(z = 0; z < PMGRID; z++)
		{
		  x = xx - slabstart_x;

		  yrr = yll = yr = yl = y;
		  zrr = zll = zr = zl = z;


		  switch (dim)
		    {
		    case 0:
		    case 3:
		      yr = y + 1;
		      yl = y - 1;
		      yrr = y + 2;
		      yll = y - 2;
		      if(yr >= PMGRID)
			yr -= PMGRID;
		      if(yrr >= PMGRID)
			yrr -= PMGRID;
		      if(yl < 0)
			yl += PMGRID;
		      if(yll < 0)
			yll += PMGRID;
		      break;
		    case 2:
		    case 1:
		    case 4:
		      yr = y + 1;
		      yl = y - 1;
		      yrr = y + 2;
		      yll = y - 2;
		      if(yr >= PMGRID)
			yr -= PMGRID;
		      if(yrr >= PMGRID)
			yrr -= PMGRID;
		      if(yl < 0)
			yl += PMGRID;
		      if(yll < 0)
			yll += PMGRID;
		      zr = z + 1;
		      zl = z - 1;
		      zrr = z + 2;
		      zll = z - 2;
		      if(zr >= PMGRID)
			zr -= PMGRID;
		      if(zrr >= PMGRID)
			zrr -= PMGRID;
		      if(zl < 0)
			zl += PMGRID;
		      if(zll < 0)
			zll += PMGRID;
		      break;
		    case 5:
		      zr = z + 1;
		      zl = z - 1;
		      zrr = z + 2;
		      zll = z - 2;
		      if(zr >= PMGRID)
			zr -= PMGRID;
		      if(zrr >= PMGRID)
			zrr -= PMGRID;
		      if(zl < 0)
			zl += PMGRID;
		      if(zll < 0)
			zll += PMGRID;
		      break;
		    }
		  if(dim == 0)
		    {
		      forcegrid[PMGRID * (x + y * nslab_x) + z] =
			-2.0 * PMGRID / All.BoxSize * (rhogrid[PMGRID * (x + yl * nslab_x) + z] +
						       rhogrid[PMGRID * (x + yr * nslab_x) + z] -
						       2.0 * rhogrid[PMGRID * (x + y * nslab_x) + z]) * fac;

		    }
		  if(dim == 1)
		    {
		      forcegrid[nslab_x * (y + z * PMGRID) + x] =
			-0.5 * fac * PMGRID / All.BoxSize * (rhogrid[nslab_x * (yl + zl * PMGRID) + x] +
							     rhogrid[nslab_x * (yr + zr * PMGRID) + x] -
							     rhogrid[nslab_x * (yl + zr * PMGRID) + x] -
							     rhogrid[nslab_x * (yr + zl * PMGRID) + x]);
		    }

		  if(dim == 2)
		    {
		      forcegrid[PMGRID * (x + y * nslab_x) + z] =
			-0.5 * fac * PMGRID / All.BoxSize * (rhogrid[PMGRID * (x + yl * nslab_x) + zl] +
							     rhogrid[PMGRID * (x + yr * nslab_x) + zr] -
							     rhogrid[PMGRID * (x + yl * nslab_x) + zr] -
							     rhogrid[PMGRID * (x + yr * nslab_x) + zl]);
		    }

		  if(dim == 3)
		    {
		      forcegrid[PMGRID2 * (PMGRID * x + y) + z] =
			-2.0 * PMGRID / All.BoxSize * (rhogrid[PMGRID2 * (PMGRID * x + yl) + zl] +
						       rhogrid[PMGRID2 * (PMGRID * x + yr) + zr] -
						       2.0 * rhogrid[PMGRID2 * (PMGRID * x + y) + z]) * fac;
		    }
		  if(dim == 4)
		    {
		      forcegrid[PMGRID2 * (y + x * PMGRID) + z] =
			-0.5 * fac * PMGRID / All.BoxSize * (rhogrid[PMGRID2 * (yl + x * PMGRID) + zl] +
							     rhogrid[PMGRID2 * (yr + x * PMGRID) + zr] -
							     rhogrid[PMGRID2 * (yl + x * PMGRID) + zr] -
							     rhogrid[PMGRID2 * (yr + x * PMGRID) + zl]);

		    }
		  if(dim == 5)
		    {
		      forcegrid[PMGRID2 * (PMGRID * x + y) + z] =
			-2.0 * PMGRID / All.BoxSize * (rhogrid[PMGRID2 * (PMGRID * x + yl) + zl] +
						       rhogrid[PMGRID2 * (PMGRID * x + yr) + zr] -
						       2.0 * rhogrid[PMGRID2 * (PMGRID * x + y) + z]) * fac;
		    }

		}

	  if((dim == 0) || (dim == 2))
	    pm_periodic_transposeB(forcegrid, tidal_workspace);

	  if(dim == 1)
	    pm_periodic_transposeBz(forcegrid, tidal_workspace);

	  /* send the force components to the right processors */

	  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	    {
	      sendTask = ThisTask;
	      recvTask = ThisTask ^ level;

	      if(recvTask < NTask)
		{
		  if(level > 0)
		    {
		      import_data =
			(fftw_real *) mymalloc("import_data", localfield_togo[recvTask * NTask + ThisTask] *
					       sizeof(fftw_real));
		      import_globalindex =
			(large_array_offset *) mymalloc("import_globalindex",
							localfield_togo[recvTask * NTask +
									ThisTask] *
							sizeof(large_array_offset));

		      if(localfield_togo[sendTask * NTask + recvTask] > 0
			 || localfield_togo[recvTask * NTask + sendTask] > 0)
			{
			  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
				       localfield_togo[sendTask * NTask +
						       recvTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_PERIODIC_C, import_globalindex,
				       localfield_togo[recvTask * NTask +
						       sendTask] * sizeof(large_array_offset), MPI_BYTE,
				       recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
			}
		    }
		  else
		    {
		      import_data = localfield_data + localfield_offset[ThisTask];
		      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
		    }

            /* PFH: I'm getting bugs/crashes in this step below in some runs. works fine for the
                first few dims, but when the dims gets down to dim=1, some task or tasks
                seem to segfault in the operation below. can't isolate it yet. using Fourier
                method instead for this step, as its cheap in our cosmological runs */
		  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
		    {
		      /* determine offset in local FFT slab */
		      offset =
			import_globalindex[i] -
			first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
		      import_data[i] = forcegrid[offset];
		    }

		  if(level > 0)
		    {
		      MPI_Sendrecv(import_data,
				   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_A,
				   localfield_data + localfield_offset[recvTask],
				   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
				   recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		      myfree(import_globalindex);
		      myfree(import_data);
		    }
		}
	    }


	  /* read out the forces, which all have been assembled in localfield_data */

	  for(i = 0, j = 0; i < NumPart; i++)
	    {
#ifdef DM_SCALARFIELD_SCREENING
	      if(phase == 1)
		if(P[i].Type == 0)	/* baryons don't get an extra scalar force */
		  continue;
#endif
	      while(j < num_on_grid && (part[j].partindex >> 3) != i)
		j++; 
	      
	      /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */ 
	      /* make sure that particles are properly box-wrapped */ 
	      for(xx = 0; xx < 3; xx++) 
	      {
		  pp[xx] = P[i].Pos[xx]; 
		  pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]); 
	      } 
	      slab_x = (int) (to_slab_fac * pp[0]); 
	      slab_y = (int) (to_slab_fac * pp[1]); 
	      slab_z = (int) (to_slab_fac * pp[2]); 
	      dx = to_slab_fac * pp[0] - slab_x; 
	      dy = to_slab_fac * pp[1] - slab_y; 
	      dz = to_slab_fac * pp[2] - slab_z; 
	      tidal_dim = + localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
		+ localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
		+ localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
		+ localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
		+ localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
		+ localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
		+ localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
		+ localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	      if(dim == 0)
		{
		  P[i].tidal_tensorpsPM[0][0] += tidal_dim;
		}
	      if(dim == 1)
		{
		  P[i].tidal_tensorpsPM[0][1] += tidal_dim;
		  P[i].tidal_tensorpsPM[1][0] += tidal_dim;
		}
	      if(dim == 2)
		{
		  P[i].tidal_tensorpsPM[0][2] += tidal_dim;
		  P[i].tidal_tensorpsPM[2][0] += tidal_dim;
		}
	      if(dim == 3)
		{
		  P[i].tidal_tensorpsPM[1][1] += tidal_dim;
		}
	      if(dim == 4)
		{
		  P[i].tidal_tensorpsPM[1][2] += tidal_dim;
		  P[i].tidal_tensorpsPM[2][1] += tidal_dim;
		}
	      if(dim == 5)
		{
		  P[i].tidal_tensorpsPM[2][2] += tidal_dim;
		}
	    }
	}

      /* free locallist */
      myfree(localfield_togo);
      myfree(localfield_offset);
      myfree(localfield_count);
      myfree(localfield_first);
      myfree(localfield_d_data);
      myfree(localfield_globalindex);
#ifdef DM_SCALARFIELD_SCREENING
    }
#endif

  pm_init_periodic_free();
  PRINT_STATUS(" ..done PM-TIDAL");
}



/*! Calculates the long-range tidal field using the PM method.  The potential is
 *  Gaussian filtered with Asmth, given in mesh-cell units. We carry out a CIC
 *  charge assignment, and compute the potenial by Fourier transform
 *  methods. The CIC kernel is deconvolved.
 *  Note that the k's need a pre-factor of 2 M_PI / All.BoxSize.
 *  The procedure calculates the second derivates of the gravitational potential by "pulling" down k's in fourier space.
 *  Component specifies the entry in the tidal field tensor that should be calculated:
 *  0=xx 1=xy 2=xz 3=yy 4=yz 5=zz
 */
/* NOTES:
   this procedure is not yet  optimized:
   ->the FFT is done in both directions for all components
 */
void pmtidaltensor_periodic_fourier(int component)
{
  double k2, kx, ky, kz, smth;
  double dx, dy, dz;
  double fx, fy, fz, ff;
  double asmth2, fac, tidal;
  MyDouble pp[3];
  int i, j, slab, level, sendTask, recvTask, task;
  int x, y, z, ip;
  int slab_x, slab_y, slab_z;
  int slab_xx, slab_yy, slab_zz;
  int num_on_grid, num_field_points, pindex, xx, yy, zz;
  MPI_Status status;
  int *localfield_count, *localfield_first, *localfield_offset, *localfield_togo;
  large_array_offset offset, *localfield_globalindex, *import_globalindex;
  d_fftw_real *localfield_d_data, *import_d_data;
  fftw_real *localfield_data, *import_data;
 
  PRINT_STATUS("Starting periodic PM-Tidaltensor (component=%d) calculation.  (presently allocated=%g MB)",component, AllocatedBytes / (1024.0 * 1024.0));
  asmth2 = (2 * M_PI) * All.Asmth[0] / All.BoxSize;
  asmth2 *= asmth2;

  fac = All.G / (M_PI * All.BoxSize);	/* to get potential */

  pm_init_periodic_allocate();


  /* determine the cells each particles accesses */
  for(i = 0, num_on_grid = 0; i < NumPart; i++)
    {
        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        if(slab_x >= PMGRID) slab_x -= PMGRID;
        if(slab_y >= PMGRID) slab_y -= PMGRID;
        if(slab_z >= PMGRID) slab_z -= PMGRID; 
	
	for(xx = 0; xx < 2; xx++) 
	    for(yy = 0; yy < 2; yy++) 
		for(zz = 0; zz < 2; zz++) { 
		    slab_xx = slab_x + xx; 
		    slab_yy = slab_y + yy; 
		    slab_zz = slab_z + zz; 
		    
		    if(slab_xx >= PMGRID) 
			slab_xx -= PMGRID; 
		    if(slab_yy >= PMGRID) 
			slab_yy -= PMGRID; 
		    if(slab_zz >= PMGRID) 
			slab_zz -= PMGRID; 
		    
		    offset = ((large_array_offset) PMGRID2) * (PMGRID * slab_xx + slab_yy) + slab_zz; 
		    part[num_on_grid].partindex = (i << 3) + (xx << 2) + (yy << 1) + zz; 
		    part[num_on_grid].globalindex = offset; 
		    part_sortindex[num_on_grid] = num_on_grid; 
		    num_on_grid++; 
		} 
    } 
  
  /* note: num_on_grid will be  8 times larger than the particle number, 
     but num_field_points will generally be much smaller */

  /* bring the part-field into the order of the accessed cells. This allow the removal of duplicates */
#ifdef MYSORT
  mysort_pmperiodic(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#else
  qsort(part_sortindex, num_on_grid, sizeof(int), pm_periodic_compare_sortindex);
#endif

  /* determine the number of unique field points */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      num_field_points++;
    }

  /* allocate the local field */
  localfield_globalindex =
    (large_array_offset *) mymalloc("localfield_globalindex", num_field_points * sizeof(large_array_offset));
  localfield_d_data = (d_fftw_real *) mymalloc("localfield_d_data", num_field_points * sizeof(d_fftw_real));
  localfield_data = (fftw_real *) localfield_d_data;
  localfield_first = (int *) mymalloc("localfield_first", NTask * sizeof(int));
  localfield_count = (int *) mymalloc("localfield_count", NTask * sizeof(int));
  localfield_offset = (int *) mymalloc("localfield_offset", NTask * sizeof(int));
  localfield_togo = (int *) mymalloc("localfield_togo", NTask * NTask * sizeof(int));

  for(i = 0; i < NTask; i++)
    {
      localfield_first[i] = 0;
      localfield_count[i] = 0;
    }

  /* establish the cross link between the part[] array and the local list of
     mesh points. Also, count on which CPU how many of the needed field points are stored */
  for(i = 0, num_field_points = 0; i < num_on_grid; i++)
    {
      if(i > 0)
	if(part[part_sortindex[i]].globalindex != part[part_sortindex[i - 1]].globalindex)
	  num_field_points++;

      part[part_sortindex[i]].localindex = num_field_points;

      if(i > 0)
	if(part[part_sortindex[i]].globalindex == part[part_sortindex[i - 1]].globalindex)
	  continue;

      localfield_globalindex[num_field_points] = part[part_sortindex[i]].globalindex;

      slab = part[part_sortindex[i]].globalindex / (PMGRID * PMGRID2);
      task = slab_to_task[slab];
      if(localfield_count[task] == 0)
	localfield_first[task] = num_field_points;
      localfield_count[task]++;
    } 
  num_field_points++;

  for(i = 1, localfield_offset[0] = 0; i < NTask; i++)
    localfield_offset[i] = localfield_offset[i - 1] + localfield_count[i - 1];

  /* now bin the local particle data onto the mesh list */

  for(i = 0; i < num_field_points; i++)
    localfield_d_data[i] = 0;

  for(i = 0; i < num_on_grid; i += 8)
    { 
	pindex = (part[i].partindex >> 3);
        if(P[pindex].Mass<=0) continue;

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[pindex].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        dx = to_slab_fac * pp[0] - slab_x;
        dy = to_slab_fac * pp[1] - slab_y;
        dz = to_slab_fac * pp[2] - slab_z;

	localfield_d_data[part[i + 0].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	localfield_d_data[part[i + 1].localindex] += P[pindex].Mass * (1.0 - dx) * (1.0 - dy) * dz;
	localfield_d_data[part[i + 2].localindex] += P[pindex].Mass * (1.0 - dx) * dy * (1.0 - dz);
	localfield_d_data[part[i + 3].localindex] += P[pindex].Mass * (1.0 - dx) * dy * dz;
	localfield_d_data[part[i + 4].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * (1.0 - dz);
	localfield_d_data[part[i + 5].localindex] += P[pindex].Mass * (dx) * (1.0 - dy) * dz;
	localfield_d_data[part[i + 6].localindex] += P[pindex].Mass * (dx) * dy * (1.0 - dz);
	localfield_d_data[part[i + 7].localindex] += P[pindex].Mass * (dx) * dy * dz;
    }

  /* clear local FFT-mesh density field */
  for(i = 0; i < fftsize; i++)
    d_rhogrid[i] = 0;

  /* exchange data and add contributions to the local mesh-path */

  MPI_Allgather(localfield_count, NTask, MPI_INT, localfield_togo, NTask, MPI_INT, MPI_COMM_WORLD);

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_d_data =
		(d_fftw_real *) mymalloc("import_d_data",
					 localfield_togo[recvTask * NTask + ThisTask] * sizeof(d_fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc("import_globalindex",
						localfield_togo[recvTask * NTask +
								ThisTask] * sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_d_data + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A,
			       import_d_data,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(d_fftw_real), MPI_BYTE,
			       recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_B, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_B, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_d_data = localfield_d_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);

	      d_rhogrid[offset] += import_d_data[i];
	    }

	  if(level > 0)
	    {
	      myfree(import_globalindex);
	      myfree(import_d_data);
	    }
	}
    }

  /* Do the FFT of the density field */

#ifndef USE_FFTW3
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
  fftw_execute(fft_forward_plan); 
#endif

  /* multiply with Green's function for the potential */

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID / 2 + 1; z++)
	{
	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      smth = -exp(-k2 * asmth2) / k2;

	      /* do deconvolution */

	      fx = fy = fz = 1;
	      if(kx != 0)
		{
		  fx = (M_PI * kx) / PMGRID;
		  fx = sin(fx) / fx;
		}
	      if(ky != 0)
		{
		  fy = (M_PI * ky) / PMGRID;
		  fy = sin(fy) / fy;
		}
	      if(kz != 0)
		{
		  fz = (M_PI * kz) / PMGRID;
		  fz = sin(fz) / fz;
		}
	      ff = 1 / (fx * fy * fz);
	      smth *= ff * ff * ff * ff;


	      /* end deconvolution */

	      ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + z;

	      /* modify greens function to get second derivatives of potential ("pulling" down k's) */
	      if(component == 0)
		{
		  cmplx_re(fft_of_rhogrid[ip]) *= smth * kx * kx;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth * kx * kx;
		}
	      if(component == 1)
		{
		  cmplx_re(fft_of_rhogrid[ip]) *= smth * kx * ky;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth * kx * ky;
		}
	      if(component == 2)
		{
		  cmplx_re(fft_of_rhogrid[ip]) *= smth * kx * kz;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth * kx * kz;
		}
	      if(component == 3)
		{
		  cmplx_re(fft_of_rhogrid[ip]) *= smth * ky * ky;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth * ky * ky;
		}
	      if(component == 4)
		{
		  cmplx_re(fft_of_rhogrid[ip]) *= smth * ky * kz;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth * ky * kz;
		}
	      if(component == 5)
		{
		  cmplx_re(fft_of_rhogrid[ip]) *= smth * kz * kz;
		  cmplx_im(fft_of_rhogrid[ip]) *= smth * kz * kz;
		}

	      /* prefactor = (2*M_PI) / All.BoxSize */
	      /* note: tidal tensor = - d^2 Phi/ dx_i dx_j  -- make sure the sign is correct here -- */
	      cmplx_re(fft_of_rhogrid[ip]) *= (2 * M_PI) * (2 * M_PI) / (All.BoxSize * All.BoxSize);
	      cmplx_im(fft_of_rhogrid[ip]) *= (2 * M_PI) * (2 * M_PI) / (All.BoxSize * All.BoxSize);

	    }
	}

  if(slabstart_y == 0)
    cmplx_re(fft_of_rhogrid[0]) = cmplx_im(fft_of_rhogrid[0]) = 0.0;

  /* Do the inverse FFT to get the tidal tensor component */

#ifndef USE_FFTW3
  rfftwnd_mpi(fft_inverse_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
  fftw_execute(fft_inverse_plan); 
#endif

  /* Now rhogrid holds the tidal tensor componet */


  /* now read out the tidal tensor component */

  /* send the tidal tensor component to the right processors */

  for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	{
	  if(level > 0)
	    {
	      import_data =
		(fftw_real *) mymalloc("import_data",
				       localfield_togo[recvTask * NTask + ThisTask] * sizeof(fftw_real));
	      import_globalindex =
		(large_array_offset *) mymalloc("import_globalindex",
						localfield_togo[recvTask * NTask +
								ThisTask] * sizeof(large_array_offset));

	      if(localfield_togo[sendTask * NTask + recvTask] > 0
		 || localfield_togo[recvTask * NTask + sendTask] > 0)
		{
		  MPI_Sendrecv(localfield_globalindex + localfield_offset[recvTask],
			       localfield_togo[sendTask * NTask + recvTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_C, import_globalindex,
			       localfield_togo[recvTask * NTask + sendTask] * sizeof(large_array_offset),
			       MPI_BYTE, recvTask, TAG_PERIODIC_C, MPI_COMM_WORLD, &status);
		}
	    }
	  else
	    {
	      import_data = localfield_data + localfield_offset[ThisTask];
	      import_globalindex = localfield_globalindex + localfield_offset[ThisTask];
	    }

	  for(i = 0; i < localfield_togo[recvTask * NTask + sendTask]; i++)
	    {
	      /* determine offset in local FFT slab */
	      offset =
		import_globalindex[i] -
		first_slab_of_task[ThisTask] * PMGRID * ((large_array_offset) PMGRID2);
	      import_data[i] = rhogrid[offset];
	    }

	  if(level > 0)
	    {
	      MPI_Sendrecv(import_data,
			   localfield_togo[recvTask * NTask + sendTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_PERIODIC_A,
			   localfield_data + localfield_offset[recvTask],
			   localfield_togo[sendTask * NTask + recvTask] * sizeof(fftw_real), MPI_BYTE,
			   recvTask, TAG_PERIODIC_A, MPI_COMM_WORLD, &status);

	      myfree(import_globalindex);
	      myfree(import_data);
	    }
	}
    }

  /* read out the tidal field values, which all have been assembled in localfield_data */

  for(i = 0, j = 0; i < NumPart; i++)
    {
      while(j < num_on_grid && (part[j].partindex >> 3) != i)
	j++;

        /* possible bugfix: Y.Feng:  (otherwise just pp[xx]=Pos[xx]) */
        /* make sure that particles are properly box-wrapped */
        for(xx = 0; xx < 3; xx++)
        {
            pp[xx] = P[i].Pos[xx];
            pp[xx] = WRAP_POSITION_UNIFORM_BOX(pp[xx]);
        }
        slab_x = (int) (to_slab_fac * pp[0]);
        slab_y = (int) (to_slab_fac * pp[1]);
        slab_z = (int) (to_slab_fac * pp[2]);
        dx = to_slab_fac * pp[0] - slab_x;
        dy = to_slab_fac * pp[1] - slab_y;
        dz = to_slab_fac * pp[2] - slab_z; 
	
	tidal = 
	    +localfield_data[part[j + 0].localindex] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz)
	    + localfield_data[part[j + 1].localindex] * (1.0 - dx) * (1.0 - dy) * dz
	    + localfield_data[part[j + 2].localindex] * (1.0 - dx) * dy * (1.0 - dz)
	    + localfield_data[part[j + 3].localindex] * (1.0 - dx) * dy * dz
	    + localfield_data[part[j + 4].localindex] * (dx) * (1.0 - dy) * (1.0 - dz)
	    + localfield_data[part[j + 5].localindex] * (dx) * (1.0 - dy) * dz
	    + localfield_data[part[j + 6].localindex] * (dx) * dy * (1.0 - dz)
	    + localfield_data[part[j + 7].localindex] * (dx) * dy * dz;

	tidal *= fac;

	if(component == 0) 
	{ 
	    P[i].tidal_tensorpsPM[0][0] += tidal;
	} 
	if(component == 1)
	{ 
	    P[i].tidal_tensorpsPM[0][1] += tidal; 
	    P[i].tidal_tensorpsPM[1][0] += tidal;
	} 
	if(component == 2)
	{
	    P[i].tidal_tensorpsPM[0][2] += tidal; 
	    P[i].tidal_tensorpsPM[2][0] += tidal;
	} 
	if(component == 3)
	{ 
	    P[i].tidal_tensorpsPM[1][1] += tidal;
	} 
	if(component == 4)
	{ 
	    P[i].tidal_tensorpsPM[1][2] += tidal; 
	    P[i].tidal_tensorpsPM[2][1] += tidal;
	} 
	if(component == 5) 
	{ 
	    P[i].tidal_tensorpsPM[2][2] += tidal;
	} 
    }

  /* free locallist */
  myfree(localfield_togo);
  myfree(localfield_offset);
  myfree(localfield_count);
  myfree(localfield_first);
  myfree(localfield_d_data);
  myfree(localfield_globalindex);

  pm_init_periodic_free();

  PRINT_STATUS(" ..done PM-Tidaltensor (component=%d).", component);
}

#endif /*COMPUTE_TIDAL_TENSOR_IN_GRAVTREE*/

/*           Here comes code for the power-spectrum computation.
 */
#define BINS_PS  2000		/* number of bins for power spectrum computation */
#define POWERSPEC_FOLDFAC 32
static long long CountModes[2][BINS_PS];
static double SumPower[2][BINS_PS];
static double SumPowerUncorrected[2][BINS_PS];	/* without binning correction (as for shot noise) */
static double Power[2][BINS_PS];
static double PowerUncorrected[2][BINS_PS];	/* without binning correction */
static double Delta[2][BINS_PS];
static double DeltaUncorrected[2][BINS_PS];	/* without binning correction */
static double ShotLimit[2][BINS_PS];
static double Kbin[BINS_PS];
static double K0, K1;
static double binfac;
static char power_spec_fname[500];
static long long power_spec_totnumpart;
static double power_spec_totmass;


void powerspec(int flag, int *typeflag)
{
  int i, n, x, y, z, kx, ky, kz, bin, ip, rep, zz;
  double k, k2, po, ponorm, smth, fac;
  double fx, fy, fz, ff, mass;
  double *powerbuf;
  long long *countbuf;
  double tstart, tend;

  PRINT_STATUS("begin power spectrum. (step=%d)  POWERSPEC_FOLDFAC=%d", flag, POWERSPEC_FOLDFAC);
  tstart = my_second();

  for(i = 0, mass = 0; i < NumPart; i++)
    if(typeflag[P[i].Type] && (P[i].Mass>0))
      mass += P[i].Mass;

  MPI_Allreduce(&mass, &power_spec_totmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  fac = 1.0 / power_spec_totmass;

  K0 = 2 * M_PI / All.BoxSize;	/* minimum k */
  K1 = K0 * All.BoxSize / (All.ForceSoftening[1] / 3.);	/* maximum k */
  binfac = BINS_PS / (log(K1) - log(K0));

  if(flag == 0)
    {
      for(rep = 0; rep < 2; rep++)
	for(i = 0; i < BINS_PS; i++)
	  {
	    SumPower[rep][i] = 0;
	    SumPowerUncorrected[rep][i] = 0;
	    CountModes[rep][i] = 0;
	  }
    }

  for(y = slabstart_y; y < slabstart_y + nslab_y; y++)
    for(x = 0; x < PMGRID; x++)
      for(z = 0; z < PMGRID; z++)
	{
	  zz = z;
	  if(z >= PMGRID / 2 + 1)
	    zz = PMGRID - z;

	  if(x > PMGRID / 2)
	    kx = x - PMGRID;
	  else
	    kx = x;
	  if(y > PMGRID / 2)
	    ky = y - PMGRID;
	  else
	    ky = y;
	  if(z > PMGRID / 2)
	    kz = z - PMGRID;
	  else
	    kz = z;

	  k2 = kx * kx + ky * ky + kz * kz;

	  if(k2 > 0)
	    {
	      if(k2 < (PMGRID / 2.0) * (PMGRID / 2.0))
		{
		  /* do deconvolution */

		  fx = fy = fz = 1;
		  if(kx != 0)
		    {
		      fx = (M_PI * kx) / PMGRID;
		      fx = sin(fx) / fx;
		    }
		  if(ky != 0)
		    {
		      fy = (M_PI * ky) / PMGRID;
		      fy = sin(fy) / fy;
		    }
		  if(kz != 0)
		    {
		      fz = (M_PI * kz) / PMGRID;
		      fz = sin(fz) / fz;
		    }
		  ff = 1 / (fx * fy * fz);
		  smth = ff * ff * ff * ff;

		  /* end deconvolution */

		  ip = PMGRID * (PMGRID / 2 + 1) * (y - slabstart_y) + (PMGRID / 2 + 1) * x + zz;

		  po = (cmplx_re(fft_of_rhogrid[ip]) * cmplx_re(fft_of_rhogrid[ip])
			+ cmplx_im(fft_of_rhogrid[ip]) * cmplx_im(fft_of_rhogrid[ip]));

		  po *= fac * fac * smth;

		  k = sqrt(k2) * 2 * M_PI / All.BoxSize;

		  if(flag == 0)
		    k *= POWERSPEC_FOLDFAC;

		  if(k >= K0 && k < K1)
		    {
		      bin = log(k / K0) * binfac;

		      ponorm = po / PowerSpec_Efstathiou(k);

		      SumPower[flag][bin] += ponorm;
		      SumPowerUncorrected[flag][bin] += po;


		      CountModes[flag][bin] += 1;

#ifdef RAYLEIGH_BINS
		      if(flag == 1)	/* only for course grid */
			{
			  ratio = sqrt(ponorm * RayleighFactor[bin]);

			  if(ratio > RayleighMax[bin])
			    RayleighMax[bin] = ratio;

			  if(ratio >= 10.0)
			    RayleighCountAbove[bin]++;
			  else
			    {
			      i = RAYLEIGH_BINS * ratio / 10.0;
			      RayleighCountModes[bin][i]++;
			    }
			}
#endif
		    }
		}
	    }
	}

  /* Now compute the power spectrum */

  countbuf = (long long int *) mymalloc("countbuf", NTask * BINS_PS * sizeof(long long));
  powerbuf = (double *) mymalloc("powerbuf", NTask * BINS_PS * sizeof(double));

  MPI_Allgather(CountModes[flag], BINS_PS * sizeof(long long), MPI_BYTE,
		countbuf, BINS_PS * sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      CountModes[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	CountModes[flag][i] += countbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPower[flag], BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPower[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	SumPower[flag][i] += powerbuf[n * BINS_PS + i];
    }

  MPI_Allgather(SumPowerUncorrected[flag], BINS_PS * sizeof(double), MPI_BYTE,
		powerbuf, BINS_PS * sizeof(double), MPI_BYTE, MPI_COMM_WORLD);

  for(i = 0; i < BINS_PS; i++)
    {
      SumPowerUncorrected[flag][i] = 0;
      for(n = 0; n < NTask; n++)
	SumPowerUncorrected[flag][i] += powerbuf[n * BINS_PS + i];
    }

  myfree(powerbuf);
  myfree(countbuf);


  for(i = 0; i < BINS_PS; i++)
    {
      Kbin[i] = exp((i + 0.5) / binfac + log(K0));

      if(CountModes[flag][i] > 0)
	{
	  Power[flag][i] = PowerSpec_Efstathiou(Kbin[i]) * SumPower[flag][i] / CountModes[flag][i];
	  PowerUncorrected[flag][i] = SumPowerUncorrected[flag][i] / CountModes[flag][i];
	}
      else
	{
	  Power[flag][i] = 0;
	  PowerUncorrected[flag][i] = 0;
	}

      Delta[flag][i] = 4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3) * Power[flag][i];

      DeltaUncorrected[flag][i] = 4 * M_PI * pow(Kbin[i], 3) /
	pow(2 * M_PI / All.BoxSize, 3) * PowerUncorrected[flag][i];

      ShotLimit[flag][i] = 4 * M_PI * pow(Kbin[i], 3) /
	pow(2 * M_PI / All.BoxSize, 3) * (1.0 / power_spec_totnumpart);
    }

  if(flag == 1)
    {
      powerspec_save();
#ifdef RAYLEIGH_BINS
      rayleigh_save();
#endif
    }

  tend = my_second();
  PRINT_STATUS("end power spectrum. (step=%d) took %g seconds", flag, timediff(tstart, tend));
}

double PowerSpec_Efstathiou(double k)
{
  double AA, BB, CC, nu, ShapeGamma;

  ShapeGamma = 0.21;
  AA = 6.4 / ShapeGamma * (1000. / UNIT_LENGTH_IN_KPC);
  BB = 3.0 / ShapeGamma * (1000. / UNIT_LENGTH_IN_KPC);
  CC = 1.7 / ShapeGamma * (1000. / UNIT_LENGTH_IN_KPC);
  nu = 1.13;


  return k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}



void calculate_power_spectra(int num, long long *ntot_type_all)
{
  int i, typeflag[6];

  power_spec_totnumpart = 0;

  for(i = 0; i < 6; i++)
    {
      typeflag[i] = 1;
      power_spec_totnumpart += ntot_type_all[i];
    }

  sprintf(power_spec_fname, "%s/powerspec_%03d.txt", All.OutputDir, num);

  pmforce_periodic(1, typeflag);	/* calculate power spectrum for all particle types */

#ifdef OUTPUT_POWERSPEC_EACH_TYPE
  if(ntot_type_all)
    for(i = 0; i < 6; i++)
      {
	if(ntot_type_all[i] > 0)
	  {
	    int j;
	    for(j = 0; j < 6; j++)
	      typeflag[j] = 0;

	    typeflag[i] = 1;
	    power_spec_totnumpart = ntot_type_all[i];

	    sprintf(power_spec_fname, "%s/powerspec_type%d_%03d.txt", All.OutputDir, i, num);

	    pmforce_periodic(1, typeflag);	/* calculate power spectrum for type i */
	  }
      }
#endif
}

void powerspec_save(void)
{
  FILE *fd;
  char buf[500];
  int i, flag;

  if(ThisTask == 0)
    {
      if(!(fd = fopen(power_spec_fname, "w")))
	{
	  sprintf(buf, "can't open file `%s`\n", power_spec_fname);
	  terminate(buf);
	}

      for(flag = 0; flag < 2; flag++)
	{
	  fprintf(fd, "%.16g\n", All.Time);
	  i = BINS_PS;
	  fprintf(fd, "%d\n", i);
	  fprintf(fd, "%g\n", power_spec_totmass);
	  fprintf(fd, "%d%09d\n", (int) (power_spec_totnumpart / 1000000000),
		  (int) (power_spec_totnumpart % 1000000000));

	  for(i = 0; i < BINS_PS; i++)
	    {
	      fprintf(fd, "%g %g %g %g %g %g %g %g %g %g\n", Kbin[i], Delta[flag][i], ShotLimit[flag][i],
		      Power[flag][i], (double) CountModes[flag][i], DeltaUncorrected[flag][i],
		      PowerUncorrected[flag][i], PowerSpec_Efstathiou(Kbin[i]), SumPower[flag][i],
		      4 * M_PI * pow(Kbin[i], 3) / pow(2 * M_PI / All.BoxSize, 3));
	    }
	}
      fclose(fd);
    }
}



void foldonitself(int *typelist)
{
  int i, j, level, sendTask, recvTask, istart, nbuf, n, rest, iter = 0;
  int slab_x, slab_xx, slab_y, slab_yy, slab_z, slab_zz;
  int *nsend_local, *nsend_offset, *nsend, count, buf_capacity;
  double to_slab_fac_folded, dx, dy, dz;
  double tstart0, tstart, tend, t0, t1;
  MyDouble pp[3];
  MyFloat *pos_sendbuf, *pos_recvbuf, *pos;
  MPI_Status status;

  PRINT_STATUS("begin folding for power spectrum estimation...");
  tstart0 = tstart = my_second();

  nsend_local = (int *) mymalloc("nsend_local", NTask * sizeof(int));
  nsend_offset = (int *) mymalloc("nsend_offset", NTask * sizeof(int));
  nsend = (int *) mymalloc("nsend", NTask * NTask * sizeof(int));

  buf_capacity = (maxfftsize * sizeof(d_fftw_real)) / (4 * sizeof(MyFloat));
  buf_capacity /= 2;

  pos_sendbuf = (MyFloat *) forcegrid;
  pos_recvbuf = pos_sendbuf + 4 * buf_capacity;

  to_slab_fac_folded = to_slab_fac * POWERSPEC_FOLDFAC;

  for(i = 0; i < fftsize; i++)	/* clear local density field */
    rhogrid[i] = 0;

  istart = 0;

  do
    {
      t0 = my_second();

      for(i = 0; i < NTask; i++)
	nsend_local[i] = 0;

      for(i = istart, nbuf = 0; i < NumPart; i++)
	{
	  if(typelist[P[i].Type] == 0)
	    continue;

	  if(nbuf + 1 >= buf_capacity)
	    break;


	  /* make sure that particles are properly box-wrapped */
	  pp[0] = P[i].Pos[0]; 
	  pp[0] = WRAP_POSITION_UNIFORM_BOX(pp[0]);

	  slab_x = to_slab_fac_folded * pp[0];
	  slab_xx = slab_x + 1;
	  slab_x %= PMGRID;
	  slab_xx %= PMGRID;

	  nsend_local[slab_to_task[slab_x]]++;
	  nbuf++;

	  if(slab_to_task[slab_x] != slab_to_task[slab_xx])
	    {
	      nsend_local[slab_to_task[slab_xx]]++;
	      nbuf++;
	    }
	}

      for(i = 1, nsend_offset[0] = 0; i < NTask; i++)
	nsend_offset[i] = nsend_offset[i - 1] + nsend_local[i - 1];

      for(i = 0; i < NTask; i++)
	nsend_local[i] = 0;

      for(i = istart, nbuf = 0; i < NumPart; i++)
	{
	  if(typelist[P[i].Type] == 0) continue;
        if(P[i].Mass <= 0) continue;

	  if(nbuf + 1 >= buf_capacity)
	    break;

	  /* make sure that particles are properly box-wrapped */
	  pp[0] = P[i].Pos[0]; 
	  pp[0] = WRAP_POSITION_UNIFORM_BOX(pp[0]);

	  slab_x = to_slab_fac_folded * pp[0];
	  slab_xx = slab_x + 1;
	  slab_x %= PMGRID;
	  slab_xx %= PMGRID;

	  for(j = 0; j < 3; j++)
	    pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_x]] + nsend_local[slab_to_task[slab_x]]) + j] =
	      P[i].Pos[j];

	  pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_x]] + nsend_local[slab_to_task[slab_x]]) + 3] =
	    P[i].Mass;

	  nsend_local[slab_to_task[slab_x]]++;
	  nbuf++;

	  if(slab_to_task[slab_x] != slab_to_task[slab_xx])
	    {
	      for(j = 0; j < 3; j++)
		pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_xx]] + nsend_local[slab_to_task[slab_xx]]) +
			    j] = P[i].Pos[j];

	      pos_sendbuf[4 * (nsend_offset[slab_to_task[slab_xx]] + nsend_local[slab_to_task[slab_xx]]) +
			  3] = P[i].Mass;


	      nsend_local[slab_to_task[slab_xx]]++;
	      nbuf++;
	    }
	}

      istart = i;


      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      t1 = my_second();
	  PRINT_STATUS("buffer filled (took %g sec)", timediff(t0, t1));
      t0 = my_second();
      for(level = 0; level < (1 << PTask); level++)	/* note: for level=0, target is the same task */
	{
	  sendTask = ThisTask;
	  recvTask = ThisTask ^ level;

	  if(recvTask < NTask)
	    {
	      if(recvTask != sendTask)
		{
		  MPI_Sendrecv(&pos_sendbuf[4 * nsend_offset[recvTask]],
			       4 * nsend_local[recvTask] * sizeof(MyFloat), MPI_BYTE,
			       recvTask, TAG_PM_FOLD,
			       &pos_recvbuf[0],
			       4 * nsend[recvTask * NTask + ThisTask] * sizeof(MyFloat), MPI_BYTE,
			       recvTask, TAG_PM_FOLD, MPI_COMM_WORLD, &status);

		  pos = &pos_recvbuf[0];
		  count = nsend[recvTask * NTask + ThisTask];
		}
	      else
		{
		  pos = &pos_sendbuf[4 * nsend_offset[ThisTask]];
		  count = nsend_local[ThisTask];
		}

	      for(n = 0; n < count; n++, pos += 4)
		{
		  /* make sure that particles are properly box-wrapped */
		  for(j = 0; j < 3; j++)
		    {
		      pp[j] = pos[j]; 
		      pp[j] = WRAP_POSITION_UNIFORM_BOX(pp[j]);
		    }

		  slab_x = to_slab_fac_folded * pp[0];
		  dx = to_slab_fac_folded * pp[0] - slab_x;
		  slab_xx = slab_x + 1;
		  slab_x %= PMGRID;
		  slab_xx %= PMGRID;

		  slab_y = to_slab_fac_folded * pp[1];
		  dy = to_slab_fac_folded * pp[1] - slab_y;
		  slab_yy = slab_y + 1;
		  slab_y %= PMGRID;
		  slab_yy %= PMGRID;

		  slab_z = to_slab_fac_folded * pp[2];
		  dz = to_slab_fac_folded * pp[2] - slab_z;
		  slab_zz = slab_z + 1;
		  slab_z %= PMGRID;
		  slab_zz %= PMGRID;

		  float mass = pos[3];

		  if(slab_to_task[slab_x] == ThisTask)
		    {
		      slab_x -= first_slab_of_task[ThisTask];

		      rhogrid[(slab_x * PMGRID + slab_y) * PMGRID2 + slab_z] +=
			mass * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
		      rhogrid[(slab_x * PMGRID + slab_yy) * PMGRID2 + slab_z] +=
			mass * (1.0 - dx) * dy * (1.0 - dz);
		      rhogrid[(slab_x * PMGRID + slab_y) * PMGRID2 + slab_zz] +=
			mass * (1.0 - dx) * (1.0 - dy) * dz;
		      rhogrid[(slab_x * PMGRID + slab_yy) * PMGRID2 + slab_zz] += mass * (1.0 - dx) * dy * dz;
		    }

		  if(slab_to_task[slab_xx] == ThisTask)
		    {
		      slab_xx -= first_slab_of_task[ThisTask];

		      rhogrid[(slab_xx * PMGRID + slab_y) * PMGRID2 + slab_z] +=
			mass * (dx) * (1.0 - dy) * (1.0 - dz);
		      rhogrid[(slab_xx * PMGRID + slab_yy) * PMGRID2 + slab_z] +=
			mass * (dx) * dy * (1.0 - dz);
		      rhogrid[(slab_xx * PMGRID + slab_y) * PMGRID2 + slab_zz] +=
			mass * (dx) * (1.0 - dy) * dz;
		      rhogrid[(slab_xx * PMGRID + slab_yy) * PMGRID2 + slab_zz] += mass * (dx) * dy * dz;
		    }

		}
	    }
	}

      count = NumPart - istart;	/* local remaining particles */
      MPI_Allreduce(&count, &rest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      iter++;

      t1 = my_second();
	  PRINT_STATUS("particles exchanged and binned. (took %g sec) max-rest=%d", timediff(t0, t1), rest);
    }
  while(rest > 0);

  tend = my_second();
  PRINT_STATUS("folded density field assembled (took %g seconds, iter=%d)", timediff(tstart, tend), iter);
  tstart = my_second();

  /* Do the FFT of the self-folded density field */
#ifndef USE_FFTW3
  rfftwnd_mpi(fft_forward_plan, 1, rhogrid, workspace, FFTW_TRANSPOSED_ORDER);
#else 
  fftw_execute(fft_forward_plan);
#endif

  tend = my_second();

  PRINT_STATUS("FFT for folded density done (took %g seconds)", timediff(tstart, tend));
  myfree(nsend);
  myfree(nsend_offset);
  myfree(nsend_local);

  tend = my_second();
}


#ifdef OUTPUT_LONGRANGE_POTENTIAL
void dump_potential(void)
{
  char buf[1000];
  int nprocgroup, primaryTask, groupTask, n, i, j, k;
  double asmth, fac, box, tstart, tend;
  float *potential;
  FILE *fd;

  tstart = my_second();

  PRINT_STATUS("Start dumping potential");
  sprintf(buf, "%s/snapdir_%03d/potential_%03d.%d", All.OutputDir, All.PowerSpecFlag - 1, All.PowerSpecFlag - 1, ThisTask);

  nprocgroup = NTask / All.NumFilesWrittenInParallel;

  if((NTask % All.NumFilesWrittenInParallel))
    nprocgroup++;

  primaryTask = (ThisTask / nprocgroup) * nprocgroup;

  for(groupTask = 0; groupTask < nprocgroup; groupTask++)
    {
      if(ThisTask == (primaryTask + groupTask))	/* ok, it's this processor's turn */
	{
	  if(!(fd = fopen(buf, "w")))
	    {
	      printf("Error. Can't write in file '%s'\n", buf);
	      endrun(11);
	    }

	  n = PMGRID;
	  fwrite(&n, sizeof(int), 1, fd);

	  n = sizeof(float);
	  fwrite(&n, sizeof(int), 1, fd);

	  fwrite(&slabs_per_task[ThisTask], sizeof(int), 1, fd);
	  fwrite(&first_slab_of_task[ThisTask], sizeof(int), 1, fd);

	  box = All.BoxSize;
	  asmth = All.Asmth[0];

	  fwrite(&box, sizeof(double), 1, fd);
	  fwrite(&asmth, sizeof(double), 1, fd);

	  fac = All.G * All.PartMass / (M_PI * All.BoxSize);

	  potential = (float *) forcegrid;

	  for(i = 0; i < slabs_per_task[ThisTask]; i++)
	    for(j = 0; j < PMGRID; j++)
	      for(k = 0; k < PMGRID; k++)
		*potential++ = fac * rhogrid[(i * PMGRID + j) * PMGRID2 + k];

	  potential = (float *) forcegrid;

	  fwrite(potential, sizeof(float), PMGRID * PMGRID * slabs_per_task[ThisTask], fd);

	  fclose(fd);
	}

      /* wait inside the group */
      MPI_Barrier(MPI_COMM_WORLD);
    }


  MPI_Barrier(MPI_COMM_WORLD);
  tend = my_second();
  PRINT_STATUS("finished writing potential (took=%g sec)", timediff(tstart, tend));
}
#endif



#endif
#endif
