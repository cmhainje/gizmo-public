#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../allvars.h"
#include "../proto.h"


/* This routine allocates memory for particle storage, both the collisionless and the fluid particles.
 * The memory for the ordered binary tree of the timeline is also allocated.
 */
/*
 * This file was originally part of the GADGET3 code developed by
 * Volker Springel. The code has been modified
 * in part (cleaned up a bit, dealt with some newer memory
 * structures and allocation strategies) by Phil Hopkins
 * (phopkins@caltech.edu) for GIZMO.
 */

void allocate_memory(void)
{
  size_t bytes;

  double bytes_tot = 0;

  int NTaskTimesThreads;

  NTaskTimesThreads = maxThreads * NTask;

  Exportflag = (int *) mymalloc("Exportflag", NTaskTimesThreads * sizeof(int));
  Exportindex = (int *) mymalloc("Exportindex", NTaskTimesThreads * sizeof(int));
  Exportnodecount = (int *) mymalloc("Exportnodecount", NTaskTimesThreads * sizeof(int));

  Send_count = (int *) mymalloc("Send_count", sizeof(int) * NTask);
  Send_offset = (int *) mymalloc("Send_offset", sizeof(int) * NTask);
  Recv_count = (int *) mymalloc("Recv_count", sizeof(int) * NTask);
  Recv_offset = (int *) mymalloc("Recv_offset", sizeof(int) * NTask);

  ProcessedFlag = (unsigned char *) mymalloc("ProcessedFlag", bytes = All.MaxPart * sizeof(unsigned char));
  bytes_tot += bytes;

  NextActiveParticle = (int *) mymalloc("NextActiveParticle", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  NextInTimeBin = (int *) mymalloc("NextInTimeBin", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;

  PrevInTimeBin = (int *) mymalloc("PrevInTimeBin", bytes = All.MaxPart * sizeof(int));
  bytes_tot += bytes;


  if(All.MaxPart > 0)
    {
      if(!(P = (struct particle_data *) mymalloc("ParticleData_P", bytes = All.MaxPart * sizeof(struct particle_data))))
	{
	  printf("failed to allocate memory for particle data storage structure `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0) {printf("Allocated %g MByte for particle data storage.\n", bytes_tot / (1024.0 * 1024.0));}
    }

  if(All.MaxPartGas > 0)
    {
      bytes_tot = 0;

      if(!(SphP = (struct gas_cell_data *) mymalloc("GasCellData_P", bytes = All.MaxPartGas * sizeof(struct gas_cell_data))))
    {
	  printf("failed to allocate memory for gas cell data storage structure (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1);
	}
      bytes_tot += bytes;

      if(ThisTask == 0) {printf("Allocated %g MByte for storage of hydro data.\n", bytes_tot / (1024.0 * 1024.0));}

#ifdef CHIMES 
      if (!(ChimesGasVars = (struct gasVariables *) mymalloc("gasVars", bytes = All.MaxPartGas * sizeof(struct gasVariables)))) 
	{
	  printf("failed to allocate memory for 'ChimesGasVars' (%g MB).\n", bytes / (1024.0 * 1024.0));
	  endrun(1); 
	}
      bytes_tot += bytes; 

      if(ThisTask == 0)
	printf("Allocated %g MByte for storage of ChimesGasVars data.\n", bytes_tot / (1024.0 * 1024.0));
#endif
    }




}
