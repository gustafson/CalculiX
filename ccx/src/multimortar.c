/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2011 Guido Dhondt                     */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h> 
#include "CalculiX.h"

/*
 Multiplication of (Bd)^T*A*Bd
*/

void multimortar(double *au, double *ad, int *irow, int *jq, int *nzs,
	   double *aubd, double *bdd, int *irowbd, int *jqbd, int *nzsbd,
	   double **aucp, double *adc, int **irowcp, int *jqc, int *nzsc,
           double *auqdt,int *irowqdt,int *jqqdt,int *nzsqdt,
	   int *neq,double *b, double *bhat,int* islavnode, int*imastnode,int*nactdof,
           int nslavs,int nmasts,int * mi){ 

  /*compteurs*/

  int i,j,k,l,m,icol,mt=mi[1]+1,nodes,nodem,row_ln,row_lm,row_ls,kflag,flag_diag=0;
  /* Different matrices */	
  int *irow_nn=NULL,*jq_nn=NULL,*irow_sn=NULL,*jq_sn=NULL,numb,
      *irow_mn=NULL,*jq_mn=NULL,*irow_mm=NULL,*jq_mm=NULL,
      *irow_sm=NULL,*jq_sm=NULL,*irow_ss=NULL,*jq_ss=NULL,
      *irowc=NULL,*irow_bdtil=NULL,*jq_bdtil=NULL,
      *irow_ssd=NULL,*irow_mmd=NULL,*jq_ssd=NULL,*jq_mmd=NULL;

  double *au_nn=NULL,*bd_nn=NULL,*au_sn=NULL,
         *au_mn=NULL,*bd_mn=NULL,*au_mm=NULL,*au_ms=NULL,
 	 *au_sm=NULL,*au_ss=NULL,*auc=NULL,*au_bdtil=NULL,
	 *au_ssd=NULL,*au_mmd=NULL;
 	
	clock_t debut;
	clock_t fin; 
  irowc = *irowcp; auc=*aucp;

 /* Flag to produce the bijection between local and global dof*/
 
	/*Au is symmetric compute the non_symmeric whole au_matrix*/
	double *au_w=NULL,*au_t=NULL;
	int *irow_w=NULL,*irow_t=NULL,*jq_w=NULL,*jq_t=NULL,nzs_t=*nzs,*mast1=NULL;
	/*transpose*/
	au_t=NNEW(double,nzs_t);
	irow_t=NNEW(int,nzs_t);
	jq_t=NNEW(int,neq[1]+1);
	jq_w=NNEW(int,neq[1]+1);
	mast1=NNEW(int,nzs_t);


for(j=0;j<nzs_t;j++){au_t[j]=au[j];}

  /* mast1 contains the rows of the original matrix,
     irowbt the columns */
  
  for(j=0;j<neq[1];j++){
    for(i=jq[j]-1;i<jq[j+1]-1;i++){
      mast1[i]=irow[i];
      irow_t[i]=j+1;
    }
  }
  
  kflag = 2;
  FORTRAN(isortiid, (mast1, irow_t, au_t, &nzs_t, &kflag));
   
  j=0;
  for(i=0; i<neq[1]; i++){
    
    if(j == nzs_t){
      for(k=i; k<neq[1]; k++) //Guido fiche <= 21.04.09
	jq_t[k] = nzs_t+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jq_t[i] = j+1;
      continue;
    }
    
    jq_t[i] = j+1;
    
    while(1){
      j++;
      if(j == nzs_t) break;
      if(mast1[j] != i+1) break;
    }
  }
  jq_t[neq[1]] = nzs_t + 1;	
  
  free(mast1);

  for (i=0;i<neq[1];i++){
    if(jq_t[i+1]-jq_t[i]>0){
   numb=jq_t[i+1]-jq_t[i]; 
   FORTRAN(isortid,(&irow_t[jq_t[i]-1],&au_t[jq_t[i]-1],&numb,&kflag));
    }
  }

  for (i=0;i<neq[1];i++){
    if(jq[i+1]-jq[i]>0){
   numb=jq[i+1]-jq[i]; 
   FORTRAN(isortid,(&irow[jq[i]-1],&au[jq[i]-1],&numb,&kflag));
    }
  }

        irow_w=NNEW(int,*nzs*2);
        au_w=NNEW(double,*nzs*2);


  	add_rect(au,irow,jq,neq[1],neq[1],
		au_t,irow_t,jq_t,neq[1],neq[1],
               &au_w,&irow_w,jq_w);
	

 int *l_flag=NULL,*n_flag=NULL,*m_flag=NULL,*s_flag=NULL,number=1;
     /*
	l_flag[i]=1 for M Master dof
	l_flag[i]=2 for S Slave dof
 	l_flag[i]=0 for N rest of the dof

	n_flag contains local N_row number
	m_flag contains local M_row number
	s_flag contains local S_row number
	*/
 	l_flag=NNEW(int,neq[1]);
	n_flag=NNEW(int,neq[1]);
	m_flag=NNEW(int,neq[1]);
	s_flag=NNEW(int,neq[1]);

   /* Fill l_flag*/
        //Master
	for (i=0;i<nmasts;i++){
         nodem=imastnode[i];
	 for (l=0;l<3;l++){ //test of dof
	    k = nactdof[mt*(nodem-1)+l+1];
            if(k>0) l_flag[k-1]=1;
         }
	}

        //Slave
	for (i=0;i<nslavs;i++){
         nodes=islavnode[i];
	 for (l=0;l<3;l++){ //test of dof
	    k = nactdof[mt*(nodes-1)+l+1];
            if(k>0) l_flag[k-1]=2;
         }
	}

	/*** Fill the local row ***/
        row_ln=0;
        row_lm=0;
        row_ls=0;
	
	/* Stock of the diagonale */	
	bd_nn=NNEW(double,neq[1]);


	au_mmd=NNEW(double,neq[1]);
	irow_mmd=NNEW(int,neq[1]);
	jq_mmd=NNEW(int,neq[1]);

	au_ssd=NNEW(double,neq[1]);
	irow_ssd=NNEW(int,neq[1]);
	jq_ssd=NNEW(int,neq[1]);

	/**** For the construction of Bhat, the new rhs***/
	double *f_s=NULL,*f_m=NULL,*f_fm=NULL;
//	bhat=NNEW(double,neq[1]);
	f_s=NNEW(double,neq[1]);
	f_m=NNEW(double,neq[1]);
	

        for (i=0;i<neq[1];i++){
         switch(l_flag[i]){
           case 0 : 	{//N			
                        bd_nn[row_ln]=ad[i];
                        bhat[i]=b[i];
			n_flag[i]=(++row_ln);
                        break;}
           case 1 : 	{//M
//                        bd_mm[row_lm]=ad[i];
			f_m[row_lm]=b[i];
			m_flag[i]=(++row_lm);	
			au_mmd[row_lm-1]=ad[i];
			jq_mmd[row_lm-1]=row_lm;
			irow_mmd[row_lm-1]=row_lm;
                        break;}
           case 2 : 	{//S
			f_s[row_ls]=b[i];
			s_flag[i]=(++row_ls);
                        au_ssd[row_ls-1]=ad[i];
			jq_ssd[row_ls-1]=row_ls;
			irow_ssd[row_ls-1]=row_ls;
			bhat[i]=b[i];
                        break;}
           default :    {printf("error l_flag\n");
                        break;} 	
         }
	
	}


//	for (i=0;i<neq[1];i++){
//		printf("ad[%d]=%f\n",i,ad[i]);
//	}
	
//	for (i=0;i<neq[1];i++) {
//	 printf("dof :%d, L_flag = %d, N_flag = %d, M_flag = %d, S_flag = %d\n",i+1,l_flag[i],n_flag[i],m_flag[i],s_flag[i]);
//	}

	jq_ssd[row_ls]=row_ls+1;
	jq_mmd[row_lm]=row_lm+1;
	
	RENEW(f_s,double,row_ls);
	RENEW(f_m,double,row_lm);
	RENEW(bd_nn,double,row_ln);

	RENEW(au_mmd,double,row_lm);
	RENEW(irow_mmd,int,row_lm);
	RENEW(jq_mmd,int,row_lm+1);


	RENEW(au_ssd,double,row_ls);
	RENEW(irow_ssd,int,row_ls);
	RENEW(jq_ssd,int,row_ls+1);

/*	for (i=0;i<row_ln;i++){
		printf("bd_nn[%d]=%f\n",i,bd_nn[i]);
	}	
	for (i=0;i<row_lm;i++){
		printf("bd_mm[%d]=%f\n",i,bd_mm[i]);
	}	
	for (i=0;i<row_ls;i++){
		printf("bd_ss[%d]=%f\n",i,au_ssd[i]);
	}
*/	
  // The maximal size of the local matrices is now defined with row_ln,row_lm and row_ls
     
	/*initialization of the differents fields to extract the 6 sub matrices - au, symmetric*/
        au_nn=NNEW(double,*nzs*2);
        au_mn=NNEW(double,*nzs*2);
        au_sn=NNEW(double,*nzs*2);
        au_mm=NNEW(double,*nzs*2);
        au_sm=NNEW(double,*nzs*2);
        au_ss=NNEW(double,*nzs*2);
	
        irow_nn=NNEW(int,*nzs*2);
        irow_mn=NNEW(int,*nzs*2);
        irow_sn=NNEW(int,*nzs*2);
        irow_mm=NNEW(int,*nzs*2);
        irow_sm=NNEW(int,*nzs*2);
        irow_ss=NNEW(int,*nzs*2);

        jq_nn=NNEW(int,row_ln+1);
        jq_mn=NNEW(int,row_ln+1);
        jq_sn=NNEW(int,row_ln+1);
        jq_mm=NNEW(int,row_lm+1);
        jq_sm=NNEW(int,row_lm+1);
        jq_ss=NNEW(int,row_ls+1);

        jq_nn[0]=1;
        jq_mn[0]=1;
        jq_sn[0]=1;
        jq_mm[0]=1;
        jq_sm[0]=1;
        jq_ss[0]=1;
        int nzs_nn=0,nzs_mn=0,nzs_sn=0,nzs_mm=0,nzs_sm=0,nzs_ss=0;
	for(j=0;j<neq[1];j++){
	   for(i=jq_w[j]-1;i<jq_w[j+1]-1;i++){
             switch(l_flag[j]){
               case 0 : {// matrices XN
                        switch(l_flag[irow_w[i]-1]){
                           case 0 : {//NN
                                    	irow_nn[nzs_nn]=n_flag[irow_w[i]-1];
                                    	au_nn[nzs_nn]=au_w[i];
					nzs_nn++;
					break;}

                           case 1 : //MN
                                    	irow_mn[nzs_mn]=m_flag[irow_w[i]-1];
                                    	au_mn[nzs_mn]=au_w[i];
					nzs_mn++;
                                    	break;

                           case 2 : //SN
                                    	irow_sn[nzs_sn]=s_flag[irow_w[i]-1];
                                    	au_sn[nzs_sn]=au_w[i];
					nzs_sn++;
                                    	break;

                           default : 	break;
                        }
                        break;}
 
              case 1 : {// matrices XM
                        switch(l_flag[irow_w[i]-1]){
 
                           case 1 : {//MM
                                    	irow_mm[nzs_mm]=m_flag[irow_w[i]-1];
                                    	au_mm[nzs_mm]=au_w[i];
					nzs_mm++;
                                    	break;}

                           case 2 : //SM
                                    	irow_sm[nzs_sm]=s_flag[irow_w[i]-1];
                                    	au_sm[nzs_sm]=au_w[i];
					nzs_sm++;
                                    	break;

                           default : 	break;
                        }
                        break;}

               case 2 : {// matrices XS
                        switch(l_flag[irow_w[i]-1]){

                           case 2 : {//SS
                                    	irow_ss[nzs_ss]=s_flag[irow_w[i]-1];
                                    	au_ss[nzs_ss]=au_w[i];
					nzs_ss++;
                                    	break;}

                           default : 	break;
                        }
                        break;}

		default : break;
             } // end switch column lflag
           } // end loop nzs column
		/* actualisation of the jq*****/
		if (n_flag[j]!=0) jq_nn[n_flag[j]]=nzs_nn+1;
		if (n_flag[j]!=0)jq_mn[n_flag[j]]=nzs_mn+1;
		if (n_flag[j]!=0)jq_sn[n_flag[j]]=nzs_sn+1;
		if (m_flag[j]!=0)jq_mm[m_flag[j]]=nzs_mm+1;
		if (m_flag[j]!=0)jq_sm[m_flag[j]]=nzs_sm+1;
		if (s_flag[j]!=0)jq_ss[s_flag[j]]=nzs_ss+1;
	} // end loop over the global column 


	/*****************************/
        /* Adapt the size of the differents fields*/
	free(au_t);free(irow_t);free(jq_t); 
	free(au_w);free(irow_w);free(jq_w);
/****************************************************************/
	   if (nzs_nn!=0){ 
           RENEW(au_nn,double,nzs_nn);
           RENEW(irow_nn,int,nzs_nn);
	   }else{
		printf("A_NN null\n");
           RENEW(au_nn,double,1);
           RENEW(irow_nn,int,1);
		}		
	   if (nzs_mn!=0){ 
           RENEW(au_mn,double,nzs_mn);
           RENEW(irow_mn,int,nzs_mn);
	   }else{
		printf("A_MN null\n");
           RENEW(au_mn,double,1);
           RENEW(irow_mn,int,1);
		}
	   if (nzs_sn!=0){ 
           RENEW(au_sn,double,nzs_sn);
           RENEW(irow_sn,int,nzs_sn);
	   }else{
		printf("A_SN null\n");
           RENEW(au_sn,double,1);
           RENEW(irow_sn,int,1);
		}
	   if (nzs_mm!=0){ 
           RENEW(au_mm,double,nzs_mm);
           RENEW(irow_mm,int,nzs_mm);
	   }else{
		printf("A_MM null\n");
           RENEW(au_mm,double,1);
           RENEW(irow_mm,int,1);
		}
	   if (nzs_sm!=0){ 
           RENEW(au_sm,double,nzs_sm);
           RENEW(irow_sm,int,nzs_sm);
	   }else{
		printf("A_SM null\n");
           RENEW(au_sm,double,1);
           RENEW(irow_sm,int,1);
		}
	   if (nzs_ss!=0){ 
           RENEW(au_ss,double,nzs_ss);
           RENEW(irow_ss,int,nzs_ss);
	   }else{
		printf("A_SS null\n");
           RENEW(au_ss,double,1);
           RENEW(irow_ss,int,1);
		}

	   /* Extraction of Bd_tild of Qd */
	   
	  for (j=0;j<*nzsbd;j++){
	    auqdt[j]=aubd[j];
	    irowqdt[j]=irowbd[j];
	  }
	  for(j=0;j<neq[1]+1;j++){
	  jqqdt[j]=jqbd[j];
	  } 
	  *nzsqdt=*nzsbd;
    
	/*	Premultiply Bd by -Dd^(-1)
		before multiplication: Bd is stored in aubd,
                                  Dd is stored in bdd
		after multiplication: -Dd^(-1).Bd is stored in auqdt
                                The diagonal contains only 1's                      
	*/

	  for(j=0; j<neq[1]; j++){
	   for(i=jqqdt[j]-1; i<jqqdt[j+1]-1; i++){
	    icol=irowqdt[i]-1;
	    auqdt[i] /= -bdd[icol];
	    //printf("aubd[%d]=%e\n",i,aubd[i]);
	    }
	  }

	

	  au_bdtil=NNEW(double,row_ls*row_lm);
	  irow_bdtil=NNEW(int,row_ls*row_lm);
	  jq_bdtil=NNEW(int,row_lm+1);
	  int nzs_bdtil=0;
	  jq_bdtil[0]=1;		
	  for(j=0;j<neq[1];j++){
	    for (i=jqqdt[j]-1;i<jqqdt[j+1]-1;i++){
	      switch(l_flag[j]){
		case 1 : //Matrix XM
			switch(l_flag[irowqdt[i]-1]){
			  case 2 : //SM here Bdtild
				    irow_bdtil[nzs_bdtil]=s_flag[irowqdt[i]-1];
				    au_bdtil[nzs_bdtil]=auqdt[i];
				    nzs_bdtil++;
				    break;
			  default : break;
				    
			}
			  default : break;
	      }
	    }
	    if (m_flag[j]!=0) jq_bdtil[m_flag[j]]=nzs_bdtil+1;
	  }
	 
	 RENEW(irow_bdtil,int,nzs_bdtil);
	 RENEW(au_bdtil,double,nzs_bdtil);


  for (i=0;i<row_lm;i++){
    if(jq_bdtil[i+1]-jq_bdtil[i]>0){
   numb=jq_bdtil[i+1]-jq_bdtil[i]; 
   FORTRAN(isortid,(&irow_bdtil[jq_bdtil[i]-1],&au_bdtil[jq_bdtil[i]-1],&numb,&kflag));
    }
  }
/*************************************** ALL THE SUB MATRICES have been yet computed ****************************
***************************************** Calculation of the submuliplication ***********************************
******************************************************************************************************************/

        
	/***************************************************************************
	************************ TEST - Multiplication *****************************
	****************************************************************************

	double a_t1[8]={1,3,4,2,4,1,1,1};
	int irow_t1[8]={1,3,4,2,3,1,4,1};
	int jq_at1[5]={1,4,6,8,9};
	double a_t2[8]={1,3,1,2,1,2,1,1},*a_t3=NULL;
	int irow_t2[8]={1,3,2,3,2,3,1,4},*irow_t3=NULL;
	int jq_at2[5]={1,3,5,7,9},*jq_at3=NULL;
	double *dum2=NULL,*dum1=NULL,*dum3=NULL;
	a_t3=NNEW(double,16);
	irow_t3=NNEW(int,16);
	jq_at3=NNEW(int,5);
	dum3=NNEW(double,3);
	flag_diag=1;
	multi_rect(a_t1,irow_t1,jq_at1,4,4,
		   a_t2,irow_t2,jq_at2,4,4,
                   &a_t3,&irow_t3,jq_at3);

	add_rect(a_t1,irow_t1,jq_at1,4,4,
		   a_t2,irow_t2,jq_at2,4,4,
                   &a_t3,&irow_t3,jq_at3);
	
	flag_diag=0;
//	delete(a_t1);delete(irow_t1);delete(jq_at1);
//	delete(a_t2);delete(irow_t2);delete(jq_at2);
	free(a_t3);free(irow_t3);free(jq_at3);
	/***************************************************************************
	****************************************************************************/


       /************* Calculation of A_MN *****************/

		/* Bt_T * A_SN*/
	 int *irow_intmn=NULL,*jq_intmn=NULL;	
	 double *au_intmn=NULL,*bdummy1=NULL,*bdummy2=NULL,*bdummy3=NULL;
	
	jq_intmn=NNEW(int,row_ln+1);
        au_intmn=NNEW(double,row_lm*row_ln);
	irow_intmn=NNEW(int,row_lm*row_ln);
        
	multi_rect(au_bdtil,irow_bdtil,jq_bdtil,row_ls,row_lm,
		au_sn,irow_sn,jq_sn,row_ls,row_ln,
               &au_intmn,&irow_intmn,jq_intmn);


       /************** Calculation of new A_MN=A_mn+BdT*A_SN **********************/

	 int *irow_fmn=NULL,*jq_fmn=NULL;	
	 double *au_fmn=NULL;
	jq_fmn=NNEW(int,row_ln+1);
	au_fmn=NNEW(double,row_lm*row_ln);
	irow_fmn=NNEW(int,row_lm*row_ln);
	
	add_rect(au_intmn,irow_intmn,jq_intmn,row_lm,row_ln,
	       au_mn,irow_mn,jq_mn,row_lm,row_ln,
               &au_fmn,&irow_fmn,jq_fmn);

   //    FORTRAN(writematrix,(au_fmn,ad,irow_fmn,jq_fmn,&row_ln,&number));
       /**************************************************
	**************************************************/

	/**************** Calculation of A_SM *******************/

	 int *irow_intsm1=NULL,*jq_intsm1=NULL;	
	 double *au_intsm1=NULL;
	 int *irow_ssf=NULL,*jq_ssf=NULL;	
	 double *au_ssf=NULL;
	
	jq_intsm1=NNEW(int,row_lm+1);
	au_intsm1=NNEW(double,row_ls*row_lm);
	irow_intsm1=NNEW(int,row_ls*row_lm);

	jq_ssf=NNEW(int,row_ls+1);
	au_ssf=NNEW(double,row_ls*row_ls);
	irow_ssf=NNEW(int,row_ls*row_ls);

	add_rect(au_ssd,irow_ssd,jq_ssd,row_ls,row_ls,
	       au_ss,irow_ss,jq_ss,row_ls,row_ls,
               &au_ssf,&irow_ssf,jq_ssf);


	multi_rect(au_ssf,irow_ssf,jq_ssf,row_ls,row_ls,
	       au_bdtil,irow_bdtil,jq_bdtil,row_ls,row_lm,
               &au_intsm1,&irow_intsm1,jq_intsm1);

//       FORTRAN(writematrix,(au_intsm1,ad,irow_intsm1,jq_intsm1,&row_lm,&number));
//       FORTRAN(writematrix,(au_ss,ad,irow_ss,jq_ss,&row_ls,&number));
	/********************* the matrix intsm1 needed for new_AMM ****************/


        /*************** Calculation of new A_SM *********************************/



	 int *irow_fsm=NULL,*jq_fsm=NULL;	
	 double *au_fsm=NULL;
	
	jq_fsm=NNEW(int,row_lm+1);
	au_fsm=NNEW(double,row_ls*row_lm);
	irow_fsm=NNEW(int,row_ls*row_lm);


	add_rect(au_sm,irow_sm,jq_sm,row_ls,row_lm,
	       au_intsm1,irow_intsm1,jq_intsm1,row_ls,row_lm,
               &au_fsm,&irow_fsm,jq_fsm);

//       FORTRAN(writematrix,(au_fsm,ad,irow_fsm,jq_fsm,&row_lm,&number));

	/********************************************************
	*********************************************************/

	/******************** Calculation of  A_MM ****************/

		/* Calculation of Bd_T*A_SM */

	 int *irow_intmm1=NULL,*jq_intmm1=NULL;	
	 double *au_intmm1=NULL,*bd_intmm1=NULL;
	 	
	jq_intmm1=NNEW(int,row_lm+1);
	bd_intmm1=NNEW(double,row_lm);
	au_intmm1=NNEW(double,row_lm*row_lm);
	irow_intmm1=NNEW(int,row_lm*row_lm);


	multi_rect(au_bdtil,irow_bdtil,jq_bdtil,row_ls,row_lm,
	       au_sm,irow_sm,jq_sm,row_ls,row_lm,
               &au_intmm1,&irow_intmm1,jq_intmm1);

//	 FORTRAN(writematrix,(au_intmm1,ad,irow_intmm1,jq_intmm1,&row_lm,&number));

		/* Calcul of A_MS*Bd = (Bd_T*A_SM)_T *************/

	 int *irow_tmm=NULL,*jq_tmm=NULL;	
	 double *au_tmm=NULL;

	int nzs_tmm=jq_intmm1[row_lm]-1;
	if (nzs_tmm!=0){
	au_tmm=NNEW(double,nzs_tmm);
	}else{
	au_tmm=NNEW(double,1);
	}
	if (nzs_tmm!=0){
	irow_tmm=NNEW(int,nzs_tmm);
	}else{
	irow_tmm=NNEW(int,1);
	}
	if (nzs_tmm!=0){
	mast1=NNEW(int,nzs_tmm);
	}else{
	mast1=NNEW(int,1);
	}
	jq_tmm=NNEW(int,row_lm+1);
  
  for(j=0;j<nzs_tmm;j++){au_tmm[j]=au_intmm1[j];}

  /* mast1 contains the rows of the original matrix,
     irowbt the columns */
  
  for(j=0;j<row_lm;j++){
    for(i=jq_intmm1[j]-1;i<jq_intmm1[j+1]-1;i++){
      mast1[i]=irow_intmm1[i];
      irow_tmm[i]=j+1;
    }
  }
  
  kflag = 2;
  FORTRAN(isortiid, (mast1, irow_tmm, au_tmm, &nzs_tmm, &kflag));
   
  j=0;
  for(i=0; i<row_lm; i++){
    
    if(j == nzs_tmm){
      for(k=i; k<row_lm; k++) //Guido fiche <= 21.04.09
	jq_tmm[k] = nzs_tmm+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jq_tmm[i] = j+1;
      continue;
    }
    
    jq_tmm[i] = j+1;
    
    while(1){
      j++;
      if(j == nzs_tmm) break;
      if(mast1[j] != i+1) break;
    }
  }
  jq_tmm[row_lm] = nzs_tmm + 1;	
 
  free(mast1);

   for (i=0;i<row_lm;i++){
    if(jq_tmm[i+1]-jq_tmm[i]>0){
   numb=jq_tmm[i+1]-jq_tmm[i]; 
   FORTRAN(isortid,(&irow_tmm[jq_tmm[i]-1],&au_tmm[jq_tmm[i]-1],&numb,&kflag));
    }
  }

//	FORTRAN(writematrix,(au_tmm,ad,irow_tmm,jq_tmm,&row_lm,&number));

		/*Calculation of  Bd_T * Ass_Bd ***/
	/******** Remark Ass_Bd : intsm1 ****/

	 int *irow_intmm3=NULL,*jq_intmm3=NULL;	
	 double *au_intmm3=NULL;
	 	

	jq_intmm3=NNEW(int,row_lm+1);
	au_intmm3=NNEW(double,row_lm*row_lm);
	irow_intmm3=NNEW(int,row_lm*row_lm);


	multi_rect(au_bdtil,irow_bdtil,jq_bdtil,row_ls,row_lm,
	       au_intsm1,irow_intsm1,jq_intsm1,row_ls,row_lm,
               &au_intmm3,&irow_intmm3,jq_intmm3);

		/*** Calculation of the new A_MM ***/

	 int *irow_fmmi=NULL,*jq_fmmi=NULL;	
	 double *au_fmmi=NULL;
	 	

	jq_fmmi=NNEW(int,row_lm+1);
	au_fmmi=NNEW(double,row_lm*row_lm);
	irow_fmmi=NNEW(int,row_lm*row_lm);

	add_rect(au_intmm1,irow_intmm1,jq_intmm1,row_lm,row_lm,
	       au_tmm,irow_tmm,jq_tmm,row_lm,row_lm,
               &au_fmmi,&irow_fmmi,jq_fmmi);	

	 int *irow_fmmi2=NULL,*jq_fmmi2=NULL;	
	 double *au_fmmi2=NULL;
	 	

	jq_fmmi2=NNEW(int,row_lm+1);
	au_fmmi2=NNEW(double,row_lm*row_lm);
	irow_fmmi2=NNEW(int,row_lm*row_lm);

	add_rect(au_fmmi,irow_fmmi,jq_fmmi,row_lm,row_lm,
	       au_intmm3,irow_intmm3,jq_intmm3,row_lm,row_lm,
               &au_fmmi2,&irow_fmmi2,jq_fmmi2);


	 int *irow_fmm=NULL,*jq_fmm=NULL;	
	 double *au_fmm=NULL;
	 int *irow_mmf=NULL,*jq_mmf=NULL;	
	 double *au_mmf=NULL;	 	

	jq_fmm=NNEW(int,row_lm+1);
	au_fmm=NNEW(double,row_lm*row_lm);
	irow_fmm=NNEW(int,row_lm*row_lm);

	jq_mmf=NNEW(int,row_lm+1);
	au_mmf=NNEW(double,row_lm*row_lm);
	irow_mmf=NNEW(int,row_lm*row_lm);


	add_rect(au_mmd,irow_mmd,jq_mmd,row_lm,row_lm,
	       au_mm,irow_mm,jq_mm,row_lm,row_lm,
               &au_mmf,&irow_mmf,jq_mmf);

	add_rect(au_fmmi2,irow_fmmi2,jq_fmmi2,row_lm,row_lm,
	       au_mmf,irow_mmf,jq_mmf,row_lm,row_lm,
               &au_fmm,&irow_fmm,jq_fmm);


//	 FORTRAN(writematrix,(au_fmm,ad,irow_fmm,jq_fmm,&row_lm,&number));

	/*******************************************************
 	*************** Free intermediate fields  **************	
	********************************************************/
	free(au_fmmi2);free(irow_fmmi2);free(jq_fmmi2);
	free(au_fmmi);free(irow_fmmi);free(jq_fmmi);
	free(au_intmm3);free(irow_intmm3);free(jq_intmm3);
	free(au_tmm);free(irow_tmm);free(jq_tmm);
	free(au_intmm1);free(irow_intmm1);free(jq_intmm1);
	free(au_intsm1);free(irow_intsm1);free(jq_intsm1);
	free(au_intmn);free(irow_intmn);free(jq_intmn);
	/**********************************************************
	***********************************************************/


/******************************************** Submultiplications done **********************************************
************************************************ ASSEMBLEE *********************************************************/


	/************** Construction of the Bhat = Qd*b 0 <=> f_mh = f_m + Bd_T*fs ************/

	/*** Local => Global topology ***/
	int *n_flagr=NULL,*m_flagr=NULL,*s_flagr=NULL;
	n_flagr=NNEW(int,row_ln);
	m_flagr=NNEW(int,row_lm);
	s_flagr=NNEW(int,row_ls);

	for(j=0;j<neq[1];j++){
		switch(l_flag[j]){
			case 0 : {n_flagr[n_flag[j]-1]=j+1;break;} //N
			case 1 : {m_flagr[m_flag[j]-1]=j+1;break;} //M
			case 2 : {s_flagr[s_flag[j]-1]=j+1;break;} //S
			default : break;
		}
	}

	/**** Calculation of f_mhat = f_m + Bd_T*f_s ****/

	double *v_r=NULL;

	multi_rectv(au_bdtil,irow_bdtil,jq_bdtil,row_ls,row_lm,f_s,&v_r); //Bd_T*f_s

	for (i=0;i<row_lm;i++) bhat[m_flagr[i]-1]=f_m[i]+v_r[i]; // Just Master part is concerned	

//	for (i=0;i<neq[1];i++) printf("multimortar bhat[%d]=%e\n",i,bhat[i]);
	/* ASSEMBLE */

//	*nzsc=neq[1]*(neq[1]+1)/2-neq[1];
	*nzsc=jq_nn[row_ln]+jq_fmn[row_ln]*2+jq_sn[row_ln]*2+jq_fmm[row_lm]+jq_fsm[row_lm]*2+jq_ss[row_ls]-9;
	int ifree=1;
	int comptn,compts,comptm;
	auc=NNEW(double,*nzsc);
	irowc=NNEW(int,*nzsc);
	mast1=NNEW(int,*nzsc);
	int loc_n,loc_s,loc_m;
	jqc[0]=1;
	for (j=0;j<neq[1];j++){
		m=j+1;
		comptn=0;compts=0;comptm=0;
		switch(l_flag[j]){
			case 0 : // insert Matrices NN, MN,SN
				    	loc_n=n_flag[j];
				    	for(l=jq_nn[loc_n-1]-1;l<jq_nn[loc_n]-1;l++){ //NN
					 insertas(&irowc,&mast1,&n_flagr[irow_nn[l]-1],&m,&ifree,nzsc,&au_nn[l],&auc);
					 comptn++;
			  		}
				    	for(l=jq_fmn[loc_n-1]-1;l<jq_fmn[loc_n]-1;l++){ //MN
					 insertas(&irowc,&mast1,&m_flagr[irow_fmn[l]-1],&m,&ifree,nzsc,&au_fmn[l],&auc);
					 insertas(&irowc,&mast1,&m,&m_flagr[irow_fmn[l]-1],&ifree,nzsc,&au_fmn[l],&auc);
					 comptn++;
			  		}		
				    	for(l=jq_sn[loc_n-1]-1;l<jq_sn[loc_n]-1;l++){ //SN
					 insertas(&irowc,&mast1,&s_flagr[irow_sn[l]-1],&m,&ifree,nzsc,&au_sn[l],&auc);
					 insertas(&irowc,&mast1,&m,&s_flagr[irow_sn[l]-1],&ifree,nzsc,&au_sn[l],&auc);
					 comptn++;
			  		}			
					adc[j]=bd_nn[loc_n-1];
				    break;
			case 1 : // insert Matrices MM and SM
				    	loc_m=m_flag[j];
				    	for(l=jq_fmm[loc_m-1]-1;l<jq_fmm[loc_m]-1;l++){ //MM
					 if ((m==m_flagr[irow_fmm[l]-1])){ //diagonal of mm 
						adc[j]=au_fmm[l];
					}else{
					 insertas(&irowc,&mast1,&m_flagr[irow_fmm[l]-1],&m,&ifree,nzsc,&au_fmm[l],&auc);
					 comptm++; }
			  		}
				    	for(l=jq_fsm[loc_m-1]-1;l<jq_fsm[loc_m]-1;l++){ //SM
					 insertas(&irowc,&mast1,&s_flagr[irow_fsm[l]-1],&m,&ifree,nzsc,&au_fsm[l],&auc);
					 insertas(&irowc,&mast1,&m,&s_flagr[irow_fsm[l]-1],&ifree,nzsc,&au_fsm[l],&auc);
					 comptm++; 
			  		}
				    break;	
			case 2 : // insert Matrix SS
				    	loc_s=s_flag[j];
				    	for(l=jq_ssf[loc_s-1]-1;l<jq_ssf[loc_s]-1;l++){ //SS
						if (m==s_flagr[irow_ssf[l]-1]){
							adc[j]=au_ssf[l];
					}else{
					 insertas(&irowc,&mast1,&s_flagr[irow_ssf[l]-1],&m,&ifree,nzsc,&au_ssf[l],&auc);
					 compts++; }
			  		}
				    break;
			default :   break;
		}	
	}
	// sort pro column 
	*nzsc=ifree-1;

/*****************************************************/
 kflag = 2;
  FORTRAN(isortiid, (mast1, irowc, auc,nzsc, &kflag));
  
  j=0;
  for(i=0; i<neq[1]; i++){
    if(j == *nzsc){
      for(k=i; k<neq[1]; k++) //Guido fiche <= 21.04.09
	jqc[k] = *nzsc+1;
      break;
    }
    
    if(mast1[j] != i+1){
      jqc[i] = j+1;
      continue;
    }
    
    jqc[i] = j+1;
    
    while(1){
      j++;
      if(j == *nzsc) break;
      if(mast1[j] != i+1) break;
    }
  }
  jqc[neq[1]] = *nzsc + 1;

	debut=clock();
  for (i=0;i<neq[1];i++){
    if(jqc[i+1]-jqc[i]>0){
   numb=jqc[i+1]-jqc[i]; 
   FORTRAN(isortid,(&irowc[jqc[i]-1],&auc[jqc[i]-1],&numb,&kflag));
    }
  }
	fin=clock();
	printf("multimortar tri fin : %f s\n",((double)(fin-debut))/CLOCKS_PER_SEC);


/******************************************************/

	
	free(mast1);
	
	
     number=6;

//   FORTRAN(writematrix,(auc,adc,irowc,jqc,&neq[1],&number));
 
  
	/*********** Free the intermediate matrices ********/
	free(au_nn);free(irow_nn);free(jq_nn);
	free(au_mn);free(irow_mn);free(jq_mn);
	free(au_sn);free(irow_sn);free(jq_sn);
	free(au_mm);free(irow_mm);free(jq_mm);
	free(au_sm);free(irow_sm);free(jq_sm);
	free(au_ss);free(irow_ss);free(jq_ss);
	free(au_fmm);free(irow_fmm);free(jq_fmm);
	free(au_fmn);free(irow_fmn);free(jq_fmn);
	free(au_fsm);free(irow_fsm);free(jq_fsm);
	free(au_bdtil);free(irow_bdtil);free(jq_bdtil);
	free(au_ssd);free(irow_ssd);free(jq_ssd);
	free(au_mmd);free(irow_mmd);free(jq_mmd);
	free(au_ssf);free(irow_ssf);free(jq_ssf);
	free(au_mmf);free(irow_mmf);free(jq_mmf);
	free(f_s);
	free(f_m);
	free(bd_nn);
	/*************************/

	/*END transmit the new stiffness matrix*/
	RENEW(auc,double,*nzsc);
	RENEW(irowc,int,*nzsc);
	
        *irowcp = irowc; *aucp=auc;

  return;
}
