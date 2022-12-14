#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../allvars.h"
#include "../../proto.h"
#include "../../kernel.h"
/*
* This file was originally part of the GADGET3 code developed by Volker Springel.
* It has been updated significantly by PFH for basic compatibility with GIZMO,
* as well as code cleanups, and accommodating new GIZMO functionality for various
* other operations. See notes in subfind.c and GIZMO User Guide for details.
*/


#ifdef SUBFIND

#include "../fof.h"
#include "subfind.h"


void subfind_distribute_groups(void)
{
  int i, nexport = 0, nimport = 0, target, ngrp, recvTask;
  group_properties *send_Group;

  /* count how many we have of each task */
  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = (Group[i].GrNr - 1) % NTask;
      if(target != ThisTask)
	Send_count[target]++;
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Recv_offset[0] = Send_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }

  send_Group = (group_properties *) mymalloc("send_Group", nexport * sizeof(group_properties));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(i = 0; i < Ngroups; i++)
    {
      target = (Group[i].GrNr - 1) % NTask;
      if(target != ThisTask)
	{
	  send_Group[Send_offset[target] + Send_count[target]] = Group[i];
	  Send_count[target]++;

	  Group[i] = Group[Ngroups - 1];
	  Ngroups--;
	  i--;
	}
    }

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      recvTask = ThisTask ^ ngrp;

      if(recvTask < NTask)
	{
	  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
	    {
	      /* get the group info */
	      MPI_Sendrecv(&send_Group[Send_offset[recvTask]],
			   Send_count[recvTask] * sizeof(group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A,
			   &Group[Ngroups + Recv_offset[recvTask]],
			   Recv_count[recvTask] * sizeof(group_properties), MPI_BYTE,
			   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }

  Ngroups += nimport;

  myfree(send_Group);
}




void subfind_distribute_particles(int mode)
{
  long nexport = 0, nimport = 0;
  int i, n, ngrp, target = 0;
  struct particle_data *partBuf;

  for(n = 0; n < NTask; n++)
    Send_count[n] = 0;

  for(n = 0; n < NumPart; n++)
    {
      if(mode == 1)
	{
	  if(P[n].targettask != ThisTask)
	    Send_count[P[n].targettask] += 1;
	}
      else if(mode == 2)
	{
	  if(P[n].origintask != ThisTask)
	    Send_count[P[n].origintask] += 1;
	}
    }

  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

  for(i = 0, Send_offset[0] = Recv_offset[0] = 0; i < NTask; i++)
    {
      nimport += Recv_count[i];
      nexport += Send_count[i];

      if(i > 0)
	{
	  Send_offset[i] = Send_offset[i - 1] + Send_count[i - 1];
	  Recv_offset[i] = Recv_offset[i - 1] + Recv_count[i - 1];
	}
    }


  if(NumPart - nexport + nimport > All.MaxPart) {PRINT_WARNING("on task=%d the maximum particle number All.MaxPart=%d is reached (NumPart=%d togo=%ld toget=%ld)",ThisTask, All.MaxPart, NumPart, nexport, nimport); endrun(8765);}

  partBuf = (struct particle_data *) mymalloc("partBuf", nexport * sizeof(struct particle_data));

  for(i = 0; i < NTask; i++)
    Send_count[i] = 0;

  for(n = 0; n < NumPart; n++)
    {
      if(mode == 0)
	{
        if(!(P[n].GrNr > Ncollective && P[n].GrNr <= TotNgroups)) {continue;}	/* particle is in small group */
	  target = (P[n].GrNr - 1) % NTask;
	}
      else if(mode == 1)
	{
	  target = P[n].targettask;
	}
      else if(mode == 2)
	{
	  target = P[n].origintask;
	}

      if(target != ThisTask)
	{
	  partBuf[Send_offset[target] + Send_count[target]] = P[n];
	  Send_count[target]++;

	  P[n] = P[NumPart - 1];
	  NumPart--;
	  n--;
	}
    }

#ifndef NO_ISEND_IRECV_IN_DOMAIN
  int n_requests = 0, max_requests = 2;

  MPI_Request *requests = (MPI_Request *) mymalloc("requests", max_requests * NTask * sizeof(MPI_Request));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(Recv_count[target] > 0)
	    MPI_Irecv(P + NumPart + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data),
		      MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	}
    }


  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(Send_count[target] > 0)
	    MPI_Isend(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data),
		      MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	}
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  if(n_requests > max_requests * NTask)
    {
      printf("Not enough memory reserved for requests: %d > %d !\n", n_requests, max_requests * NTask);
      endrun(52099);
    }

  myfree(requests);

#else

  /*
     double t0, t1;

     t0 = my_second();
   */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      /*
         if(NTask > 512)
         {
         if(ThisTask == 0)
         {
         t1 = my_second();
         printf("ngrp=%d %g sec\n", ngrp, timediff(t0, t1));
         fflush(stdout);
         t0 = my_second();
         }
         }
       */
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  MPI_Sendrecv(partBuf + Send_offset[target], Send_count[target] * sizeof(struct particle_data),
		       MPI_BYTE, target, TAG_PDATA,
		       P + NumPart + Recv_offset[target], Recv_count[target] * sizeof(struct particle_data),
		       MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
    }

#endif

  NumPart += nimport;

  myfree(partBuf);
}



void subfind_exchange(void)
{
  int count_togo = 0, count_togo_gas = 0, count_get = 0, count_get_gas = 0;
  int *count, *count_gas, *offset, *offset_gas;
  int *count_recv, *count_recv_gas, *offset_recv, *offset_recv_gas;
  int i, n, ngrp, target;
  int *local_toGet, *local_toGetGas;
  int *local_toGo, *local_toGoGas;
  struct particle_data *partBuf;
  struct gas_cell_data *gasBuf;

  count = (int *) mymalloc("count", NTask * sizeof(int));
  count_gas = (int *) mymalloc("count_gas", NTask * sizeof(int));
  offset = (int *) mymalloc("offset", NTask * sizeof(int));
  offset_gas = (int *) mymalloc("offset_gas", NTask * sizeof(int));

  count_recv = (int *) mymalloc("count_recv", NTask * sizeof(int));
  count_recv_gas = (int *) mymalloc("count_recv_gas", NTask * sizeof(int));
  offset_recv = (int *) mymalloc("offset_recv", NTask * sizeof(int));
  offset_recv_gas = (int *) mymalloc("offset_recv_gas", NTask * sizeof(int));

  local_toGet = (int *) mymalloc("local_toGet", NTask * sizeof(int));
  local_toGetGas = (int *) mymalloc("local_toGetGas", NTask * sizeof(int));
  local_toGo = (int *) mymalloc("local_toGo", NTask * sizeof(int));
  local_toGoGas = (int *) mymalloc("local_toGoGas", NTask * sizeof(int));

  for(n = 0; n < NTask; n++)
    {
      local_toGo[n] = 0;
      local_toGoGas[n] = 0;
    }


  for(n = 0; n < NumPart; n++)
    {
      if(P[n].targettask != ThisTask)
	{
	  local_toGo[P[n].targettask]++;
	  if(P[n].Type == 0)
	    local_toGoGas[P[n].targettask]++;
	}
    }

  MPI_Alltoall(local_toGo, 1, MPI_INT, local_toGet, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall(local_toGoGas, 1, MPI_INT, local_toGetGas, 1, MPI_INT, MPI_COMM_WORLD);

  int prec_offset, prec_count;
  int *decrease;

  decrease = (int *) mymalloc("decrease", NTask * sizeof(int));

  for(i = 1, offset_gas[0] = 0, decrease[0] = 0; i < NTask; i++)
    {
      offset_gas[i] = offset_gas[i - 1] + local_toGoGas[i - 1];
      decrease[i] = local_toGoGas[i - 1];
    }

  prec_offset = offset_gas[NTask - 1] + local_toGoGas[NTask - 1];

  offset[0] = prec_offset;

  for(i = 1; i < NTask; i++)
    offset[i] = offset[i - 1] + (local_toGo[i - 1] - decrease[i]);

  myfree(decrease);

  for(i = 0; i < NTask; i++)
    {
      count_togo += local_toGo[i];
      count_togo_gas += local_toGoGas[i];

      count_get += local_toGet[i];
      count_get_gas += local_toGetGas[i];

    }

  if(NumPart + count_get - count_togo > All.MaxPart)
    {
      printf("Task=%d NumPart(%d)+count_get(%d)-count_togo(%d)=%d  All.MaxPart=%d\n", ThisTask, NumPart,
	     count_get, count_togo, NumPart + count_get - count_togo, All.MaxPart);
      endrun(12787878);
    }

  if(N_gas + count_get_gas - count_togo_gas > All.MaxPartGas)
    {
      printf("Task=%d N_gas(%d)+count_get_gas(%d)-count_togo_gas(%d)=%d All.MaxPartGas=%d\n", ThisTask,
	     N_gas, count_get_gas, count_togo_gas, N_gas + count_get_gas - count_togo_gas, All.MaxPartGas);
      endrun(712187879);
    }

  partBuf = (struct particle_data *) mymalloc("partBuf", count_togo * sizeof(struct particle_data));
  gasBuf = (struct gas_cell_data *) mymalloc("gasBuf", count_togo_gas * sizeof(struct gas_cell_data));

  for(i = 0; i < NTask; i++)
    count[i] = count_gas[i] = 0;


  for(n = 0; n < NumPart; n++)
    {
      if(P[n].targettask != ThisTask)
	{
	  target = P[n].targettask;

	  if(P[n].Type == 0)
	    {
	      partBuf[offset_gas[target] + count_gas[target]] = P[n];
	      gasBuf[offset_gas[target] + count_gas[target]] = SphP[n];
	      count_gas[target]++;
	    }
	  else
	    {
	      partBuf[offset[target] + count[target]] = P[n];
	      count[target]++;
	    }


	  if(P[n].Type == 0)
	    {
	      P[n] = P[N_gas - 1];
	      SphP[n] = SphP[N_gas - 1];

	      P[N_gas - 1] = P[NumPart - 1];

	      NumPart--;
	      N_gas--;
	      n--;
	    }
	  else
	    {
	      P[n] = P[NumPart - 1];
	      NumPart--;
	      n--;
	    }
	}
    }

  int count_totget;

  count_totget = count_get_gas;

  if(count_totget)
    memmove(P + N_gas + count_totget, P + N_gas, (NumPart - N_gas) * sizeof(struct particle_data));

  for(i = 0; i < NTask; i++)
    {
      count_recv_gas[i] = local_toGetGas[i];
      count_recv[i] = local_toGet[i] - local_toGetGas[i];
    }

  for(i = 1, offset_recv_gas[0] = N_gas; i < NTask; i++)
    offset_recv_gas[i] = offset_recv_gas[i - 1] + count_recv_gas[i - 1];
  prec_count = N_gas + count_get_gas;

  offset_recv[0] = NumPart - N_gas + prec_count;

  for(i = 1; i < NTask; i++)
    offset_recv[i] = offset_recv[i - 1] + count_recv[i - 1];


#ifndef NO_ISEND_IRECV_IN_DOMAIN

  int n_requests = 0, max_requests = 6;
  MPI_Request *requests;

  requests = (MPI_Request *) mymalloc("requests", max_requests * NTask * sizeof(MPI_Request));

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_recv_gas[target] > 0)
	    {
	      MPI_Irecv(P + offset_recv_gas[target], count_recv_gas[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_GAS, MPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Irecv(SphP + offset_recv_gas[target],
			count_recv_gas[target] * sizeof(struct gas_cell_data), MPI_BYTE, target,
			TAG_GASDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }


	  if(count_recv[target] > 0)
	    {
	      MPI_Irecv(P + offset_recv[target], count_recv[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }
	}
    }


  MPI_Barrier(MPI_COMM_WORLD);	/* not really necessary, but this will guarantee that all receives are
				   posted before the sends, which helps the stability of MPI on
				   bluegene, and perhaps some mpich1-clusters */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_gas[target] > 0)
	    {
	      MPI_Isend(partBuf + offset_gas[target], count_gas[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA_GAS, MPI_COMM_WORLD, &requests[n_requests++]);

	      MPI_Isend(gasBuf + offset_gas[target], count_gas[target] * sizeof(struct gas_cell_data),
			MPI_BYTE, target, TAG_GASDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }

	  if(count[target] > 0)
	    {
	      MPI_Isend(partBuf + offset[target], count[target] * sizeof(struct particle_data),
			MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, &requests[n_requests++]);
	    }
	}
    }

  MPI_Waitall(n_requests, requests, MPI_STATUSES_IGNORE);

  if(n_requests > max_requests * NTask)
    {
      printf("Not enough memory reserved for requests: %d > %d !\n", n_requests, max_requests * NTask);
      endrun(52098);
    }

  myfree(requests);

#else /* NO_ISEND_IRECV_IN_DOMAIN */

  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
    {
      target = ThisTask ^ ngrp;

      if(target < NTask)
	{
	  if(count_gas[target] > 0 || count_recv_gas[target] > 0)
	    {
	      MPI_Sendrecv(partBuf + offset_gas[target], count_gas[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA_GAS,
			   P + offset_recv_gas[target], count_recv_gas[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA_GAS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	      MPI_Sendrecv(gasBuf + offset_gas[target], count_gas[target] * sizeof(struct gas_cell_data),
			   MPI_BYTE, target, TAG_GASDATA,
			   SphP + offset_recv_gas[target],
			   count_recv_gas[target] * sizeof(struct gas_cell_data), MPI_BYTE, target,
			   TAG_GASDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }


	  if(count[target] > 0 || count_recv[target] > 0)
	    {
	      MPI_Sendrecv(partBuf + offset[target], count[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA,
			   P + offset_recv[target], count_recv[target] * sizeof(struct particle_data),
			   MPI_BYTE, target, TAG_PDATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}
    }


#endif


  NumPart += count_get;
  if(NumPart > All.MaxPart)
    {
      printf("Task=%d NumPart=%d All.MaxPart=%d\n", ThisTask, NumPart, All.MaxPart);
      endrun(12787878);
    }

  N_gas += count_get_gas;
  if(N_gas > All.MaxPartGas)
    {
      printf("Task=%d N_gas=%d All.MaxPartGas=%d\n", ThisTask, N_gas, All.MaxPartGas);
      endrun(712187879);
    }


  myfree(gasBuf);
  myfree(partBuf);

  myfree(local_toGoGas);
  myfree(local_toGo);
  myfree(local_toGetGas);
  myfree(local_toGet);

  myfree(offset_recv_gas);
  myfree(offset_recv);
  myfree(count_recv_gas);
  myfree(count_recv);

  myfree(offset_gas);
  myfree(offset);
  myfree(count_gas);
  myfree(count);
}





















#endif
