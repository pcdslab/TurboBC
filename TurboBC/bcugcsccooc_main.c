/*
 *  Copyright (C) 2020 Oswaldo Artiles and Fahad Saeed 
 *  Florida International University, Florida, USA.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 3
 * of the License, or (at your option) any later version.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE   
 * Please refer to the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/* 
*   Betweenness centrality (BC) 
*   Single precision (float data type)
*   TurboBC:bcugcsccooc_main.c
*
*   This program:
*   1) Reads a real (non-complex) sparse matrix from a Matrix Market (v. 2.0) file.
*   2) computes a sequential BC for unweighted, and weighted
*      graphs with sparse matrices in the CSC  format.            
*   3) computes a single precision GPU-based parallel BC 
*      for undirected and directed graphs with sparse 
*      matrices in the CSC and COOC formats.
*   
*   Compile: make with Makefile
*
*   Run: ./exec /lclhome/oarti001/documents/SpMVcuda/MMdataFiles/[filename]  --w --ug --format --p --nr --rs --repet --seq
*   
*   --w is an integer equal to 1 if there is a weighted graph, and 0 otherwise,
*   --ug is an integer equal to 1 if there is an undirected graph (matrix A is symmetric),
*   and equal to 0 otherwise,
*   --format: integer equal to: 0 (CSC_sc), 1 (CSC_wa), 2 (COOC_sc)
*   --p: integer equal to 1 to print results, 0 otherwise,
*   --nr: number root vertices to compute the BC (also source vertex)
*   --rs: if nr =1, rs is the root vertex.
*   --repet: the number of computations of BC to obtain average time.
*   --seq: 1 when the sequential BC is computed, 0 otherwise
*
*   run with the: run_bcugcsccooc.sh script
 *   
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include "timer.h"
#include "mmiof.h"
#include "bcug.h"
#include "sparseformatransf.h"
#include "readMMfilef.h"
#include "degree.h"
#include "utils_bc.h"

int main(int argc, char *argv[]){
   
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M,N;
    long nz;
    int i,j;
    int *RP,*S,*S_hgpu;
    int *out_degree,*in_degree,*degree;
    int *degreeDist,*degreeDistAc;
    double *degreeDistRel,*degreeDistRelAc;
    int *ICLT,*JCLT,*CPLT,*RPLT;
    struct Csc CscA;
    struct Cooc CoocA;
    float *valCLT,*valC,*sigma,*sigma_hgpu,*bc,*bc_hgpu;    
    int maxDegree = 0,idxMaxDegree;
    double meanDegree,sdDegree;
    int maxInDegree = 0,idxMaxInDegree,maxOutDegree = 0,idxMaxOutDegree;
    double meanInDegree,sdInDegree,meanOutDegree,sdOutDegree;
    int w = atoi(argv[2]);
    int ug = atoi(argv[3]);
    int format = atoi(argv[4]);
    int p = atoi(argv[5]);
    int nr = atoi(argv[6]);
    int rs = atoi(argv[7]);
    int repet = atoi(argv[8]);
    int seq = atoi(argv[9]);

    if (argc < 9){
      printf("number of arguments not valid, should be 9.\n");
      fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
      exit(1);	
    }else{        
       if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }
    
    if (mm_read_banner(f, &matcode) != 0){
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) ){
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* compute size of sparse matrix .... */
    printf("compute size of sparse adjacency matrix starting....\n ");
    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0) exit(1);        
    if (ug){nz = 2*nz;}
    printf("\nN = %d,M = %d, nz = %ld\n",N,M,nz);

    /* allocate memory for vectors */   
    CscA.IC = (int *) calloc(nz,sizeof(*CscA.IC));
    CscA.CP = (int *) calloc((N+1),sizeof(*CscA.CP));
    CoocA.ICOOC = (int *) calloc(nz,sizeof(*CoocA.ICOOC));
    CoocA.JCOOC = (int *) calloc(nz,sizeof(*CoocA.JCOOC));    
    if (w) valC = (float *) calloc(nz,sizeof(*valC));
    out_degree = (int *) calloc(N,sizeof(*out_degree));
    in_degree = (int *) calloc(N,sizeof(*in_degree));
    RP = (int *) calloc((N+1),sizeof(*RP));
    sigma = (float *) calloc(N,sizeof(*sigma));
    sigma_hgpu = (float *) calloc(N,sizeof(*sigma_hgpu));
    bc = (float *) calloc(N,sizeof(*bc));
    bc_hgpu = (float *) calloc(N,sizeof(*bc_hgpu));
    S = (int *) calloc(N,sizeof(*S));
    S_hgpu = (int *) calloc(N,sizeof(*S_hgpu));
    if (ug){
      CPLT = (int *) calloc((N+1),sizeof(*CPLT));
      RPLT = (int *) calloc((N+1),sizeof(*RPLT)); 
      ICLT = (int *) calloc(nz,sizeof(*ICLT));
      JCLT = (int *) calloc(nz,sizeof(*JCLT));   
      if (w) valCLT = (float *) calloc(nz,sizeof(*valCLT));
      degree = (int *) calloc(N,sizeof(*degree));
    }

    /***************************************************************************
    read Matrix Market file and compute the in-degree and out-degree vectors            
    ***************************************************************************/      
    //printf("reading MM file starting....\n ");
    if(ug){
      if(w){
      readMMfile_w (f,ICLT,JCLT,valCLT,out_degree,in_degree,nz/2);
      }else{
	readMMfile_uw (f,ICLT,JCLT,out_degree,in_degree,nz/2);
      }
    }else{
      if(w){
      readMMfile_w (f,CoocA.ICOOC,CoocA.JCOOC,valC,out_degree,in_degree,nz);
      for (i=0; i<nz; i++) CscA.IC[i] = CoocA.ICOOC[i];
      }else{
	readMMfile_uw (f,CoocA.ICOOC,CoocA.JCOOC,out_degree,in_degree,nz);
	for (i=0; i<nz; i++) CscA.IC[i] = CoocA.ICOOC[i];
      } 
    }

    /**************************************************************************
     close Matrix Market file       
    **************************************************************************/
    if (f !=stdin) fclose(f);
       
    if (ug){
      /************************************************************************
       compute degree vector and parameters (maximum, mean and standard 
       deviation) for undirected graphs 
      ************************************************************************/
      degree_ug (out_degree,in_degree,degree,N);
      degree_max (degree,&maxDegree,&idxMaxDegree,N);
      degree_sta (degree,&meanDegree,&sdDegree,N);
      
      printf("maximum degree = %d of node = %d\n", maxDegree,idxMaxDegree);
      printf("degree statistics:meanDegree=%6.2f,sdDegree = %6.2f,\n",meanDegree,sdDegree);

      /************************************************************************
       compute degree distribution and relative degree distribution vectors for 
        undirected graphs 
      ************************************************************************/
      degreeDist = (int *) calloc((maxDegree+1),sizeof(*degreeDist));
      degreeDistRel = (double *) calloc((maxDegree+1),sizeof(*degreeDistRel));
      degreeDistAc = (int *) calloc(maxDegree,sizeof(*degreeDistAc));
      degreeDistRelAc = (double *) calloc(maxDegree,sizeof(*degreeDistRelAc));
     
      degree_dist (degree,degreeDist,N);
      degree_dist_rel (degreeDist,degreeDistRel,maxDegree,N);
      degree_dist_ac (degreeDist,degreeDistRel,degreeDistAc,degreeDistRelAc,maxDegree);
      
      
      j = 0;
      printf("Degree %d: %d (%1.3f %%)\n",j,degreeDistAc[j],degreeDistRelAc[j]);
      for (i=1; i<=maxDegree; i=i*2){
	printf("Degree 2^%d: %d (%1.3f %%))\n",j,degreeDistAc[i],degreeDistRelAc[i]);j++;}
      /************************************************************************
       compute CPLT, RPLT, CscA.CP and RP arrays for undirected graphs
      ************************************************************************/
      CP_RP_ug (out_degree,in_degree,degree,CPLT,RPLT,CscA.CP,RP,N);
       	
      /************************************************************************
       Transform Symmetric Sparse Column (SSC) to Compressed Sparse Column (CSC) 
       format for symmetric sparse adjacency  matrices representing  
       undirected graphs
      ************************************************************************/
      SSCCSC_uw (ICLT,CPLT,RPLT,RP,CscA.CP,CscA.IC,N);

      /************************************************************************
       Transform Symmetric Sparse Column (SSC) to Coordinated Ordered by Column 
       (COOC) format for symmetric sparse adjacency  matrices representing  
       undirected graphs                   
      ************************************************************************/
      SSCCOOC_uw (ICLT,CPLT,RPLT,RP,CscA.CP,CoocA.ICOOC,CoocA.JCOOC,N);
      
    }
    if (!ug){
       /************************************************************************
        compute degree distribution parameters (maximum, mean and standard 
        deviation) for directed graphs 
      ************************************************************************/
      degree_max (in_degree,&maxInDegree,&idxMaxInDegree,N);
      degree_max (out_degree,&maxOutDegree,&idxMaxOutDegree,N);
      degree_sta (in_degree,&meanInDegree,&sdInDegree,N);
      degree_sta (out_degree,&meanOutDegree,&sdOutDegree,N); 
      
      printf("maximum in_degree = %d of node = %d\n", maxInDegree,idxMaxInDegree);
      printf("maximum out_degree = %d of node = %d\n", maxOutDegree,idxMaxOutDegree);
      printf("degree statistics:meanInDegree=%6.2f,sdInDegree = %6.2f,\n",
	     meanInDegree,sdInDegree);
      printf("degree statistics:meanOutDegree=%6.2f,sdOutDegree = %6.2f,\n",
	     meanOutDegree,sdOutDegree);
      /************************************************************************
       compute degree distribution and relative degree distribution vectors for 
       directed graphs 
      ************************************************************************/      
      int *inDegreeDist = (int *) calloc((maxInDegree+1),sizeof(*inDegreeDist));
      double *inDegreeDistRel = (double *) malloc((maxInDegree+1)*sizeof(*inDegreeDistRel));
      int *outDegreeDist = (int *) calloc((maxOutDegree+1),sizeof(*outDegreeDist));
      double *outDegreeDistRel = (double *) malloc((maxOutDegree+1)*sizeof(*outDegreeDistRel));
      int *inDegreeDistAc = (int *) calloc(maxInDegree,sizeof(*inDegreeDistAc));
      double *inDegreeDistRelAc = (double *) calloc(maxInDegree,sizeof(*inDegreeDistRelAc));
      int *outDegreeDistAc = (int *) calloc(maxOutDegree,sizeof(*outDegreeDistAc));
      double *outDegreeDistRelAc = (double *) calloc(maxOutDegree,sizeof(*outDegreeDistRelAc));
      degree_dist (in_degree,inDegreeDist,N);
      degree_dist (out_degree,outDegreeDist,N);
      degree_dist_rel (inDegreeDist,inDegreeDistRel,maxInDegree,N);
      degree_dist_rel (outDegreeDist,outDegreeDistRel,maxOutDegree,N);
      degree_dist_ac (inDegreeDist,inDegreeDistRel,inDegreeDistAc,inDegreeDistRelAc,maxInDegree);
      degree_dist_ac (outDegreeDist,outDegreeDistRel,outDegreeDistAc,outDegreeDistRelAc,maxOutDegree);
      
      printf("accumulated in-degree distribution\n");
      printf("Degree %d: %d (%1.3f %%)\n",0,inDegreeDistAc[0],inDegreeDistRelAc[0]);
      j = 0;
      for (i=1;i<=maxInDegree; i=i*2){
	printf("Degree 2^%d: %d (%1.3f %%))\n",j,inDegreeDistAc[i],inDegreeDistRelAc[i]);j++;}
      printf("\naccumulated out-degree distribution\n");
      printf("Degree %d: %d (%1.3f %%)\n",0,outDegreeDistAc[0],outDegreeDistRelAc[0]);
      j = 0;
      for (i=1;i<=maxOutDegree; i=i*2){
	printf("Degree 2^%d: %d (%1.3f %%))\n",j,outDegreeDistAc[i],outDegreeDistRelAc[i]);j++;}
      /**************************************************************************
       compute CscA.CP and RP arrays for directed graphs
      **************************************************************************/
      CP_RP_dg (out_degree,in_degree,CscA.CP,RP,N);
      
    }
    
    /**************************************************************************
       compute sequential BC for unweighted graphs represented by sparse adjacency
       matrices in the CSC format
    **************************************************************************/
    printf("computing sequential BC starting");
    if (seq){// CSC format	 
      bc_seq_ug_csc (CscA.IC,CscA.CP,S,sigma,bc,nr,rs,nz,N);
    }          
         
    /**************************************************************************
      compute GPU-based parallel BC for unweighted graphs represented by sparse
      adjacency matrices in the CSC and COOC formats.         
    **************************************************************************/        
    printf("computing GPU-based parallel BC starting:nr=%d\n",nr);
    if (format == 0){// CSC(sc) format
      bc_gpu_ug_csc_sc (CscA.IC,CscA.CP,S_hgpu,sigma_hgpu,bc_hgpu,nr,rs,nz,N,repet);
      if (seq){
	printf("\ncheck sigma in CSC(sc) format: ");
        bfs_check(sigma,sigma_hgpu,N);
	printf("\ncheck values of S_hgpu: ");
	S_check(S,S_hgpu,N);
	printf("\ncheck bc in CSC(sc) format: ");
	bc_check(bc,bc_hgpu,N);	   
      }
    }else if (format == 1){ // CSC(wa) format  
      bc_gpu_ug_csc_wa (CscA.IC,CscA.CP,S_hgpu,sigma_hgpu,bc_hgpu,nr,rs,nz,N,repet);
      if (seq){
	  printf("\ncheck sigma in CSC(wa) format: ");
	  bfs_check(sigma,sigma_hgpu,N);
	  printf("\ncheck values of S_hgpu: ");
	  S_check(S,S_hgpu,N);
	  printf("\ncheck bc in CSC(wa) format: ");
	  bc_check(bc,bc_hgpu,N);
      }
    }else if (format == 2){ // COOC(sc) format  
      bc_gpu_ug_cooc_sc (CoocA.ICOOC,CoocA.JCOOC,S_hgpu,sigma_hgpu,bc_hgpu,nr,rs,nz,N,repet);
      if (seq){
	  printf("\ncheck sigma in COOC(sc) format: ");
	  bfs_check(sigma,sigma_hgpu,N);
	  printf("\ncheck values of S_hgpu: ");
	  S_check(S,S_hgpu,N);
	  printf("\ncheck bc in COOC(sc) format: ");
	  bc_check(bc,bc_hgpu,N);
      }
    }
    
    /**************************************************************************
      print out sparse formats and results
    **************************************************************************/
    
    if (p) {printBC(CscA.IC,CscA.CP,CoocA.ICOOC,CoocA.JCOOC,S,S_hgpu,sigma,
		    sigma_hgpu,bc,bc_hgpu,nz,N);}
    
    
    return 0;
}
