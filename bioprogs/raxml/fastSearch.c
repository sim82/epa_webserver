/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>



#include "axml.h"

#define POSTERIOR_THRESHOLD 0.90

extern int Thorough;
extern double masterTime;


extern char permFileName[1024], resultFileName[1024], 
  logFileName[1024], checkpointFileName[1024], infoFileName[1024], run_id[128], workdir[1024], bootStrapFile[1024], bootstrapFileName[1024], 
  bipartitionsFileName[1024],bipartitionsFileNameBranchLabels[1024]; 



typedef struct 
{
  nodeptr p;
  double lh;
} scores;


typedef struct 
{
  scores *s;
  int count;
  int maxCount; 
} insertions;

static void addInsertion(nodeptr p, double lh, insertions *ins)
{ 
  if(ins->count < ins->maxCount)
    {
      ins->s[ins->count].lh = lh;
      ins->s[ins->count].p = p;
      ins->count = ins->count + 1;               
    }
  else
    {        
      ins->s = realloc(ins->s, sizeof(scores) * ins->maxCount * 2);     

      ins->maxCount *= 2;

      ins->s[ins->count].lh = lh;
      ins->s[ins->count].p = p;
      ins->count = ins->count + 1; 
    }
}


/* the two functions below just compute the subtree insertion likelihood score into 
   a branch using the thorough method.
*/

static void insertFast (tree *tr, nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r, s;
  int i;
  
  r = q->back;
  s = p->back;
      
  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];
  
  if(Thorough)
    {
      double  zqr[NUM_BRANCHES], zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;      
      double defaultArray[NUM_BRANCHES];	
      double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
      double *qz;
      
      qz = q->z;
      
      for(i = 0; i < numBranches; i++)
	defaultArray[i] = defaultz;
      
      makenewzGeneric(tr, q, r, qz, iterations, zqr, FALSE);           
      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs, FALSE);                  
      makenewzGeneric(tr, r, s, defaultArray, iterations, zrs, FALSE);
      
      
      for(i = 0; i < numBranches; i++)
	{
	  lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
	  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
	  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
	  lzsum = 0.5 * (lzqr + lzqs + lzrs);
	  
	  lzq = lzsum - lzrs;
	  lzr = lzsum - lzqs;
	  lzs = lzsum - lzqr;
	  lzmax = log(zmax);
	  
	  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
	  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
	  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
	  
	  e1[i] = exp(lzq);
	  e2[i] = exp(lzr);
	  e3[i] = exp(lzs);
	}
      hookup(p->next,       q, e1, numBranches);
      hookup(p->next->next, r, e2, numBranches);
      hookup(p,             s, e3, numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES]; 
      
      for(i = 0; i < numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  
	  if(z[i] < zmin) 
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);	                         
    }
  
  newviewGeneric(tr, p);
  
  if(Thorough)
    {          
      localSmooth(tr, p, smoothings);   
      
      for(i = 0; i < numBranches; i++)
	{
	  tr->lzq[i] = p->next->z[i];
	  tr->lzr[i] = p->next->next->z[i];
	  tr->lzs[i] = p->z[i];            
	}           
    }
}




static double testInsertFast (tree *tr, nodeptr p, nodeptr q, insertions *ins, boolean veryFast)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r, s;
  double LH;
  int i;
  
  r = q->back; 
  
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
   
  insertFast(tr, p, q, tr->numBranches);
      
  evaluateGeneric(tr, p->next->next);   

  addInsertion(q, tr->likelihood, ins);
  
  
  if(veryFast)
    if(tr->likelihood > tr->endLH)
      {			  
	tr->insertNode = q;
	tr->removeNode = p;   
	for(i = 0; i < tr->numBranches; i++)
	  tr->currentZQR[i] = tr->zqr[i];      
	tr->endLH = tr->likelihood;                            
      }  
 
  LH = tr->likelihood;                  
              
  hookup(q, r, qz, tr->numBranches);
      
  p->next->next->back = p->next->back = (nodeptr) NULL;
  
  if(Thorough)
    {
      s = p->back;
      hookup(p, s, pz, tr->numBranches);          
    }
      
  return LH;
}

static double testInsertCandidates(tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r, s;
  double LH;
  int i;
  
  r = q->back; 
  
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
   
  insertFast(tr, p, q, tr->numBranches);
      
  evaluateGeneric(tr, p->next->next);       

  if(tr->likelihood > tr->endLH)
    {			  
      tr->insertNode = q;
      tr->removeNode = p;   
      for(i = 0; i < tr->numBranches; i++)
	tr->currentZQR[i] = tr->zqr[i];      
      tr->endLH = tr->likelihood;                            
    }  
 
  LH = tr->likelihood;                  
              
  hookup(q, r, qz, tr->numBranches);
      
  p->next->next->back = p->next->back = (nodeptr) NULL;
  
  if(Thorough)
    {
      s = p->back;
      hookup(p, s, pz, tr->numBranches);          
    }
      
  return LH;
}


/* 
   function below just computes all subtree roots, 
   we go to every inner node and store the three outgoing subtree 
   roots. Assumes evidently that nodeptr p that is passed to the first 
   instance of this recursion is not a tip ;-)
*/

static void getSubtreeRoots(nodeptr p, nodeptr *ptr, int *count, int mxtips)
{
  if(isTip(p->number, mxtips))
    return;
  
  ptr[*count] = p;
  *count = *count + 1;
  
  ptr[*count] = p->next;
  *count = *count + 1;

  ptr[*count] = p->next->next;
  *count = *count + 1;
 
  getSubtreeRoots(p->next->back, ptr, count, mxtips);
  getSubtreeRoots(p->next->next->back, ptr, count, mxtips);
}
  

/*
  conducts linear SPRs, computes the approximate likelihood score 
  fo an insertion of the subtree being rearranged into the left and 
  right branch. Depending on the score it then descends into 
  either the right or the left tree. 
  
  I also added a radius that limits the number of branches away from the original pruning position into which the tree 
  will be inserted.

  We actually need to test empirically what a good setting may be !
*/
  

static void insertBeyond(tree *tr, nodeptr p, nodeptr q, int radius, insertions *ins, boolean veryFast)
{
  if(radius > 0)
    {
      int count = 0;
      
      double 
	twoScores[2];
      
      twoScores[0] = unlikely;
      twoScores[1] = unlikely;
      
      if(isTip(q->next->back->number, tr->mxtips) && isTip(q->next->next->back->number, tr->mxtips))
	return;
      
      if(!isTip(q->next->back->number, tr->mxtips))    
	{
	  twoScores[0] = testInsertFast(tr, p, q->next->back, ins, veryFast);    
	  count++;
	}
      
      if(!isTip(q->next->next->back->number, tr->mxtips))    
	{
	  twoScores[1] = testInsertFast(tr, p, q->next->next->back, ins, veryFast);
	  count++;
	}          
      
      if(count == 2 && !veryFast)
	{
	  double 
	    weight;

	  if(twoScores[0] > twoScores[1])	    
	    {
	      weight = exp(twoScores[0] - twoScores[0]) / (exp(twoScores[0] - twoScores[0]) + exp(twoScores[1] - twoScores[0]));
	      
	      if(weight >= POSTERIOR_THRESHOLD)
		insertBeyond(tr, p, q->next->back, radius - 1, ins, veryFast);
	      else
		{
		  insertBeyond(tr, p, q->next->back, radius - 1, ins, veryFast);	  
		  insertBeyond(tr, p, q->next->next->back, radius - 1, ins, veryFast);
		}
	    }
	  else
	    {
	      weight = exp(twoScores[1] - twoScores[1]) / (exp(twoScores[1] - twoScores[0]) + exp(twoScores[1] - twoScores[1]));
	      
	      if(weight >= POSTERIOR_THRESHOLD)
		insertBeyond(tr, p, q->next->next->back, radius - 1, ins, veryFast);
	      else
		{
		  insertBeyond(tr, p, q->next->back, radius - 1, ins, veryFast);	  
		  insertBeyond(tr, p, q->next->next->back, radius - 1, ins, veryFast);
		}
	    }
	}
      else
	{
	  if(twoScores[0] > twoScores[1])
	    insertBeyond(tr, p, q->next->back, radius - 1, ins, veryFast);
	  else
	    insertBeyond(tr, p, q->next->next->back, radius - 1, ins, veryFast);
	}
    }

}



typedef struct 
{
  int direction;
  double likelihood;

} fourLikelihoods;






static int fourCompare(const void *p1, const void *p2)
{
  fourLikelihoods *rc1 = (fourLikelihoods *)p1;
  fourLikelihoods *rc2 = (fourLikelihoods *)p2;

  double i = rc1->likelihood;
  double j = rc2->likelihood;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}

static int scoreCompare(const void *p1, const void *p2)
{
  scores *rc1 = (scores *)p1;
  scores *rc2 = (scores *)p2;

  double i = rc1->lh;
  double j = rc2->lh;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}




static double linearSPRs(tree *tr, int radius, boolean veryFast)
{
  int 
    numberOfSubtrees = (tr->mxtips - 2) * 3,
    count = 0,
    k,
    i;

  double 
    fourScores[4];

  nodeptr 
    *ptr = (nodeptr *)malloc(sizeof(nodeptr) * numberOfSubtrees);
  
  fourLikelihoods 
    *fourLi = (fourLikelihoods *)malloc(sizeof(fourLikelihoods) * 4);

  insertions 
    *ins = (insertions*)malloc(sizeof(insertions));

  

  ins->count = 0;
  ins->maxCount = 2048;

  ins->s = (scores *)malloc(sizeof(scores) * ins->maxCount);


  /* recursively compute the roots of all subtrees in the current tree 
     and store them in ptr */

  getSubtreeRoots(tr->start->back, ptr, &count, tr->mxtips);

  assert(count == numberOfSubtrees);
 
  tr->startLH = tr->endLH = tr->likelihood;      

  /* loop over subtrees, i.e., execute a full SPR cycle */

  for(i = 0; i < numberOfSubtrees; i++)
    {           
      nodeptr 
	p = ptr[i],
	p1 = p->next->back,
	p2 = p->next->next->back;
      
      double   
	p1z[NUM_BRANCHES], 
	p2z[NUM_BRANCHES];
	
      ins->count = 0;

      /*printf("Node %d %d\n", p->number, i);*/

      tr->bestOfNode = unlikely;  

      for(k = 0; k < 4; k++)
	fourScores[k] = unlikely;
      
      assert(!isTip(p->number, tr->rdta->numsp));                   
      
      if(!isTip(p1->number, tr->rdta->numsp) || !isTip(p2->number, tr->rdta->numsp))
	{
	  double 
	    max = unlikely;

	  int 
	    maxInt = -1;	  

	  for(k = 0; k < tr->numBranches; k++)
	    {
	      p1z[k] = p1->z[k];
	      p2z[k] = p2->z[k];	   	   
	    }
	  
	  /* remove the current subtree */

	  removeNodeBIG(tr, p,  tr->numBranches);  

	  /* pre score with fast insertions */
	  if(veryFast)
	    Thorough = 1;
	  else
	    Thorough = 0;

	  if (!isTip(p1->number, tr->rdta->numsp)) 
	    {
	      fourScores[0] = testInsertFast(tr, p, p1->next->back, ins, veryFast);
	      fourScores[1] = testInsertFast(tr, p, p1->next->next->back, ins, veryFast);		        
	    }
	  
	  if (!isTip(p2->number, tr->rdta->numsp)) 
	    {
	      fourScores[2] = testInsertFast(tr, p, p2->next->back, ins, veryFast);
	      fourScores[3] = testInsertFast(tr, p, p2->next->next->back, ins, veryFast);			          
	    }
	  
	  if(veryFast)
	    Thorough = 1;
	  else
	    Thorough = 0;

	  /* find the most promising direction */
	  

	  if(!veryFast)
	    {
	      int 
		j = 0,
		validEntries = 0;

	      double 
		lmax = unlikely,
		posterior = 0.0;	  

	      for(k = 0; k < 4; k++)
		{
		  fourLi[k].direction = k;
		  fourLi[k].likelihood = fourScores[k];
		}

	      qsort(fourLi, 4, sizeof(fourLikelihoods), fourCompare);
	      
	      for(k = 0; k < 4; k++)
		if(fourLi[k].likelihood > unlikely)
		  validEntries++;
	      
	      lmax = fourLi[0].likelihood;		  	     

	      while(posterior <= POSTERIOR_THRESHOLD && j < validEntries)	  
		{ 	       	      
		  double 
		    all = 0.0,
		    prob = 0.0;

		  for(k =  0; k < validEntries; k++) 	   
		    all += exp(fourLi[k].likelihood - lmax);	     
	      
		  posterior += (prob = (exp(fourLi[j].likelihood - lmax) / all));		  

		  switch(fourLi[j].direction)
		    { 
		    case 0:
		      insertBeyond(tr, p, p1->next->back, radius, ins, veryFast);
		      break;
		    case 1:
		      insertBeyond(tr, p, p1->next->next->back, radius, ins, veryFast);
		      break;
		    case 2:
		      insertBeyond(tr, p, p2->next->back, radius, ins, veryFast);
		      break;
		    case 3:
		      insertBeyond(tr, p, p2->next->next->back, radius, ins, veryFast);
		      break;
		    default:
		      assert(0);
		    }
	      		  	      
		  j++;
		}	 	    

	      qsort(ins->s, ins->count, sizeof(scores), scoreCompare);	     

	      Thorough = 1;

	      for(k = 0; k < MIN(ins->count, 20); k++)
		testInsertCandidates(tr, p, ins->s[k].p);	      	      
	    }
	  else
	    {
	      Thorough = 1;
	      

	      for(k = 0; k < 4; k++)
		{
		  if(max < fourScores[k])
		    {
		      max = fourScores[k];
		      maxInt = k;
		    }
		}
	      
	      /* descend into this direction and re-insert subtree there */
	      
	      if(maxInt >= 0)
		{
		  switch(maxInt)
		    { 
		    case 0:
		      insertBeyond(tr, p, p1->next->back, radius, ins, veryFast);
		      break;
		    case 1:
		      insertBeyond(tr, p, p1->next->next->back, radius, ins, veryFast);
		      break;
		    case 2:
		      insertBeyond(tr, p, p2->next->back, radius, ins, veryFast);
		      break;
		    case 3:
		      insertBeyond(tr, p, p2->next->next->back, radius, ins, veryFast);
		      break;
		    default:
		      assert(0);
		    }
		}
	    }	      
	  
	  /* repair branch and reconnect subtree to its original position from which it was pruned */

	  hookup(p->next,       p1, p1z, tr->numBranches); 
	  hookup(p->next->next, p2, p2z, tr->numBranches);	  	 
	  
	  /* repair likelihood vectors */

	  newviewGeneric(tr, p);

	  /* if the rearrangement of subtree rooted at p yielded a better likelihood score 
	     restore the altered topology and use it from now on 
	  */

	  if(tr->endLH > tr->startLH)                 	
	    {			   	     
	      restoreTreeFast(tr);	 	 
	      tr->startLH = tr->endLH = tr->likelihood;	 	       
	    }
	   	    	
	}             
    }
 
  return tr->startLH;     
}




static boolean allSmoothed(tree *tr)
{
  int i;
  boolean result = TRUE;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      if(tr->partitionSmoothed[i] == FALSE)
	result = FALSE;
      else
	tr->partitionConverged[i] = TRUE;
    }

  return result;
}

static void nniSmooth(tree *tr, nodeptr p, int maxtimes)
{
  int
    i;

  for(i = 0; i < tr->numBranches; i++)	
    tr->partitionConverged[i] = FALSE;	

 

  while (--maxtimes >= 0) 
    {     
      

      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;
      
      

      assert(!isTip(p->number, tr->mxtips)); 	

     

      assert(!isTip(p->back->number, tr->mxtips));  
      
      update(tr, p);
     
      update(tr, p->next);
     
      update(tr, p->next->next);
      
      update(tr, p->back->next);
      
      update(tr, p->back->next->next);           
     
      if (allSmoothed(tr)) 
	break;
      
    }

  

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }

  
}


static void storeBranches(tree *tr, nodeptr p, double *pqz, double *pz1, double *pz2, double *qz1, double *qz2)
{
  int 
    i;
  
  nodeptr 
    q = p->back;

  for(i = 0; i < tr->numBranches; i++)
    {
      pqz[i] = p->z[i];
      pz1[i] = p->next->z[i];
      pz2[i] = p->next->next->z[i];
      qz1[i] = q->next->z[i];
      qz2[i] = q->next->next->z[i];
    }  
}


static int SHSupport(int nPos, int nBootstrap, int *col, double loglk[3], double *siteloglk[3]) 
{
  double
    resampleDelta,
    resample1,
    resample2,  
    support = 0.0,
    delta1 = loglk[0] - loglk[1],
    delta2 = loglk[0] - loglk[2],
    delta = delta1 < delta2 ? delta1 : delta2;
  
  int
    iBest,
    i,
    j,
    nSupport = 0,
    iBoot;
  
  assert(loglk[0] >= loglk[1] && loglk[0] >= loglk[2]);

  /*printf("%f %f %f\n", loglk[0], loglk[1], loglk[2]);*/

  for(iBoot = 0; iBoot < nBootstrap; iBoot++) 
    {
      double resampled[3];
      
      for (i = 0; i < 3; i++)
	resampled[i] = -loglk[i];
      
      for (j = 0; j < nPos; j++) 
	{
	  int pos = col[iBoot * nPos + j];
	  for (i = 0; i < 3; i++)
	    resampled[i] += pos * siteloglk[i][j];
	}

      iBest = 0;
      
      /*printf("%d %f %f %f\n", iBoot, resampled[0], resampled[1], resampled[2]);*/

      for (i = 1; i < 3; i++)
	if (resampled[i] > resampled[iBest])
	  iBest = i;
      
      resample1 = resampled[iBest] - resampled[(iBest+1)%3];
      resample2 = resampled[iBest] - resampled[(iBest+2)%3];
      resampleDelta = resample1 < resample2 ? resample1 : resample2;
      
      if(resampleDelta < delta)
	nSupport++;
    }
  
  support = (nSupport/(double)nBootstrap);
    
  return ((int)((support * 100.0) + 0.5));
}

static void setupBranchInfo(nodeptr p, tree *tr, int *counter)
{
  if(!isTip(p->number, tr->mxtips))       
    {
      nodeptr q;          

      if(!(isTip(p->back->number, tr->mxtips)))
	{	 	  
	  p->bInf = p->back->bInf = &(tr->bInf[*counter]);	  	  	  	 	       	  

	  p->bInf->oP = p;
	  p->bInf->oQ = p->back;
	  
	  *counter = *counter + 1;
	}
      
      q = p->next;

      while(q != p)
	{
	  setupBranchInfo(q->back, tr, counter);	
	  q = q->next;
	}
        
      return;
    }
}

static void doNNIs(tree *tr, nodeptr p, double *lhVectors[3], boolean shSupport, int *interchanges, boolean brOpt)
{  
  if(isTip(p->number, tr->mxtips))
    return;

  {
    nodeptr q = p->back;

    assert(!isTip(p->number, tr->mxtips));

    if(!isTip(p->number, tr->mxtips) && !isTip(q->number, tr->mxtips))
      {	
	int 
	  whichNNI;
	       
	double 		 
	  lh[3],	 
	  pqz_0[NUM_BRANCHES],
	  pz1_0[NUM_BRANCHES],
	  pz2_0[NUM_BRANCHES],
	  qz1_0[NUM_BRANCHES],
	  qz2_0[NUM_BRANCHES],
	  pqz_1[NUM_BRANCHES],
	  pz1_1[NUM_BRANCHES],
	  pz2_1[NUM_BRANCHES],
	  qz1_1[NUM_BRANCHES],
	  qz2_1[NUM_BRANCHES],
	  pqz_2[NUM_BRANCHES],
	  pz1_2[NUM_BRANCHES],
	  pz2_2[NUM_BRANCHES],
	  qz1_2[NUM_BRANCHES],
	  qz2_2[NUM_BRANCHES];

	nodeptr
	  pb1 = p->next->back,
	  pb2 = p->next->next->back,
	  qb1 = q->next->back,
	  qb2 = q->next->next->back;
	
	if(brOpt)
	  nniSmooth(tr, p, 16);	

	if(shSupport)
	  {	   
	    evaluateGenericVector(tr, p);
	    memcpy(lhVectors[0], tr->perSiteLL, sizeof(double) * tr->cdta->endsite);
	  }
	else
	  evaluateGeneric(tr, p);
	
	

	lh[0] =tr->likelihood;
	
	storeBranches(tr, p, pqz_0, pz1_0, pz2_0, qz1_0, qz2_0);
	

	
	whichNNI = 0;

	/*******************************************/

	hookup(p, q, pqz_0, tr->numBranches); 

	hookup(p->next,       qb1, qz1_0, tr->numBranches); 
	hookup(p->next->next, pb2, pz2_0, tr->numBranches); 

	hookup(q->next,       pb1, pz1_0, tr->numBranches); 	
	hookup(q->next->next, qb2, qz2_0, tr->numBranches); 

	newviewGeneric(tr, p);
	newviewGeneric(tr, p->back);
	
	if(brOpt)
	  nniSmooth(tr, p, 16);
     
	if(shSupport)
	  {
	    evaluateGenericVector(tr, p);
	    memcpy(lhVectors[1], tr->perSiteLL, sizeof(double) * tr->cdta->endsite);
	  }
	else
	  evaluateGeneric(tr, p);
	
	lh[1] = tr->likelihood;		
	
	storeBranches(tr, p, pqz_1, pz1_1, pz2_1, qz1_1, qz2_1);

	if(lh[1] > lh[0])
	  whichNNI = 1;
	
	/*******************************************/

	hookup(p, q, pqz_0, tr->numBranches); 

	hookup(p->next,       qb1, qz1_0, tr->numBranches); 
	hookup(p->next->next, pb1, pz1_0, tr->numBranches); 

	hookup(q->next,       pb2, pz2_0, tr->numBranches); 		
	hookup(q->next->next, qb2, qz2_0, tr->numBranches); 
	
	newviewGeneric(tr, p);
	newviewGeneric(tr, p->back);

	if(brOpt)
	  nniSmooth(tr, p, 16);
	
	if(shSupport)
	  {
	    evaluateGenericVector(tr, p);
	    memcpy(lhVectors[2], tr->perSiteLL, sizeof(double) * tr->cdta->endsite);
	  }
	else
	  evaluateGeneric(tr, p);
	
	lh[2] =tr->likelihood;

	storeBranches(tr, p, pqz_2, pz1_2, pz2_2, qz1_2, qz2_2);	 
	     
	if(lh[2] > lh[0] && lh[2] > lh[1])
	  whichNNI = 2;

	/*******************************************/

	/*printf("%f %f %f NNI: %d\n", lh[0], lh[1], lh[2], whichNNI);*/

	switch(whichNNI)
	  {
	  case 0:
	    hookup(p, q, pqz_0, tr->numBranches); 	  
	    
	    hookup(p->next,       pb1, pz1_0, tr->numBranches); 
	    hookup(p->next->next, pb2, pz2_0, tr->numBranches); 
	    
	    hookup(q->next,       qb1, qz1_0, tr->numBranches); 	
	    hookup(q->next->next, qb2, qz2_0, tr->numBranches);
	    break;
	  case 1:
	    hookup(p, q, pqz_1, tr->numBranches); 	    
	    
	    hookup(p->next,       qb1, pz1_1, tr->numBranches); 
	    hookup(p->next->next, pb2, pz2_1, tr->numBranches); 
	    
	    hookup(q->next,       pb1, qz1_1, tr->numBranches); 	
	    hookup(q->next->next, qb2, qz2_1, tr->numBranches); 
	    break;
	  case 2:	   
	    hookup(p, q, pqz_2, tr->numBranches); 
	    
	    hookup(p->next,       qb1, pz1_2, tr->numBranches); 
	    hookup(p->next->next, pb1, pz2_2, tr->numBranches); 
	    
	    hookup(q->next,       pb2, qz1_2, tr->numBranches); 		
	    hookup(q->next->next, qb2, qz2_2, tr->numBranches); 
	    break;
	  default:
	    assert(0);
	  }       

	newviewGeneric(tr, p);
	newviewGeneric(tr, p->back);
	evaluateGeneric(tr, p);

	if(whichNNI > 0)
	  *interchanges = *interchanges + 1;

	/*
	  if(whichNNI > 0)
	  printf("PAPA %f\n", tr->likelihood);
	else
	printf("Branch %f \n", tr->likelihood);
	*/
	
	if(shSupport)
	  {
	    int support = SHSupport(tr->cdta->endsite, 1000, tr->resample, lh, lhVectors);
	    p->bInf->support = support;
	    /*printf("support: %d\n", support);*/
	  }

	
	doNNIs(tr, pb1, lhVectors, shSupport, interchanges, brOpt);		
	doNNIs(tr, pb2, lhVectors, shSupport, interchanges, brOpt);     	
      }
    else
      {   	
	doNNIs(tr, p->next->back, lhVectors, shSupport, interchanges, brOpt);		
	doNNIs(tr, p->next->next->back, lhVectors, shSupport, interchanges, brOpt);     		
      }
  }
}

static int *PermutationSH(tree *tr, int nBootstrap) 
{
  int 
    *buffer,
    *col = (int*)calloc(tr->cdta->endsite * nBootstrap, sizeof(int)),    
    *nonzero = (int*)calloc(tr->NumberOfModels, sizeof(int)),
    maxNonZero = 0,
    model,
    i,
    replicate;

  size_t 
    bufferSize;

  long 
    randomSeed = 12345;
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;

      for(i = 0; i < width; i++)
	nonzero[model] += tr->partitionData[model].wgt[i];

      if(nonzero[model] > maxNonZero)
	maxNonZero = nonzero[model];
    }

  bufferSize = ((size_t)maxNonZero) * sizeof(int);
  buffer = (int*)malloc(bufferSize);
   
  for(replicate = 0; replicate < nBootstrap; replicate++)
    {      
      int offset = 0;

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int
	    w, 
	    pos,
	    width = tr->partitionData[model].upper - tr->partitionData[model].lower,
	    *wgt = &col[tr->cdta->endsite * replicate + offset];
	  
	  memset(buffer, 0, bufferSize);

	  for(i = 0; i < nonzero[model]; i++)
	    buffer[(int) (nonzero[model] * randum(&randomSeed))]++; 

	  for(i = 0, pos = 0; i < width; i++)
	    for(w = 0; w < tr->partitionData[model].wgt[i]; w++)	  	  	 
	      {
		wgt[i] += buffer[pos];
		pos++;		      
	      }

	  offset += width;
	}
    }
  
  free(buffer);
  free(nonzero);

  return col;
}




void fastSearch(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  double    
    likelihood, 
    startLikelihood,
    *lhVectors[3];

  char 
    bestTreeFileName[1024],
    shSupportFileName[1024];
  
  FILE 
    *f;

  int
    interchanges;

 
  if(adef->shSupports)
    {     
      tr->resample = PermutationSH(tr, 1000);
      
      lhVectors[0] = (double *)malloc(sizeof(double) * tr->cdta->endsite);
      lhVectors[1] = (double *)malloc(sizeof(double) * tr->cdta->endsite);
      lhVectors[2] = (double *)malloc(sizeof(double) * tr->cdta->endsite);
      tr->bInf = (branchInfo*)malloc(sizeof(branchInfo) * (tr->mxtips - 3));      
    }
  else
    {
      lhVectors[0] = (double *)NULL;
      lhVectors[1] = (double *)NULL;
      lhVectors[2] = (double *)NULL;
    }
      
  /* initialize model parameters with standard starting values */

  initModel(tr, rdta, cdta, adef);      

  printBothOpen("Time after init : %f\n", gettime() - masterTime);

  /* 
     compute starting tree, either by reading in a tree specified via -t 
     or by building one 
  */

  getStartingTree(tr, adef);
 
  printBothOpen("Time after init and starting tree: %f\n", gettime() - masterTime);
  
  /* 
     rough model parameter optimization, the log likelihood epsilon should 
     actually be determined based on the initial tree score and not be hard-coded 
  */

  modOpt(tr, adef, FALSE, 10.0);
  
  printBothOpen("Time after init, starting tree, mod opt: %f\n", gettime() - masterTime);

  /* print out the number of rate categories used for the CAT model, one should 
     use less then the default, e.g., -c 16 works quite well */

  printBothOpen("Cats: %d\n", tr->NumberOfCategories);

  /* 
     means that we are going to do thorough insertions 
     with real newton-raphson based br-len opt at the three branches 
     adjactent to every insertion point 
  */

  Thorough = 1;

  
  /*
    loop over SPR cycles until the likelihood difference 
     before and after the SPR cycle is <= 0.5 log likelihood units.
     Rather than being hard-coded this should also be determined based on the 
     actual likelihood of the tree 
  */
 
  do
    {      
      startLikelihood = tr->likelihood;
   
      /* conduct a cycle of linear SPRs */

    

      likelihood = linearSPRs(tr, 20, adef->veryFast);          
      
     

      /* optimize br-lens of resulting topology a bit */
      evaluateGeneric(tr, tr->start); 
               
      interchanges = 0;
      doNNIs(tr, tr->start->back, lhVectors, FALSE, &interchanges, TRUE);
	           
      /*treeEvaluate(tr, 1);             */

      printBothOpen("%f\n", tr->likelihood);
    }
  while(ABS(likelihood - startLikelihood) > 0.5);

  if(adef->shSupports)
    {
      int 
	counter = 0;

      do
	{
	  interchanges = 0;
	  doNNIs(tr, tr->start->back, lhVectors, FALSE, &interchanges, FALSE);
	  
	  evaluateGeneric(tr, tr->start); 
	  
	  /*printf("Inter %d %f\n", interchanges, tr->likelihood);*/
	}
      while(interchanges > 0);

      printBothOpen("Final Likelihood SH-Supports: %f\n", tr->likelihood);

      setupBranchInfo(tr->start->back, tr, &counter);
      assert(counter == tr->mxtips - 3);
      interchanges = 0;
      doNNIs(tr, tr->start->back, lhVectors, TRUE, &interchanges, FALSE);
    }
  
 

  
  /* print out the resulting tree to the RAxML_bestTree. file. 
     note that boosttrapping or doing multiple inferences won't work.
     This thing computes a single tree and that's it */
      
  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_fastTree.");
  strcat(bestTreeFileName,         run_id);

 

  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, adef, SUMMARIZE_LH, FALSE);
    
  f = myfopen(bestTreeFileName, "wb");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);  

  if(adef->shSupports)
    {
      strcpy(shSupportFileName, workdir); 
      strcat(shSupportFileName, "RAxML_fastTreeSH_Support.");
      strcat(shSupportFileName,         run_id);
    
      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, adef, SUMMARIZE_LH, TRUE);
    
      f = myfopen(shSupportFileName, "wb");
      fprintf(f, "%s", tr->tree_string);
      fclose(f);  

    }
    
  printBothOpen("RAxML fast tree written to file: %s\n", bestTreeFileName);
  if(adef->shSupports)
    printBothOpen("RAxML fast tree with SH-like supports written to file: %s\n", shSupportFileName);
  
  printBothOpen("Total execution time: %f\n", gettime() - masterTime);

  printBothOpen("Good bye ... \n");
}
