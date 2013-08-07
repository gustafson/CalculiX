/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2013 Guido Dhondt                     */

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
#include <string.h>
#include <time.h>
#include "CalculiX.h"

/**
 *Calculate the entries of Bd and Dd, and insert them into the data structure
 * 
 * @param [out] irowbdp		field containing row numbers of aubd
 * @param [out] jqbd		pointer into field irowbd
 * @param [out] aubdp		pointer to matrix Bd
 * @param [out] bdd		matrix Dd
 * @param [out] nzsbd		size of aubd
 * @param [in] ntie		number of contraints
 * @param [in] ipkon		pointer into field kon...
 * @param [in] kon 		.. for element i storing the connectivity list of elem. in succ. order
 * @param [in] lakon		(i) label for element i
 * @param [in] nslavnode	(i)pointer into field isalvnode for contact tie i 
 * @param [in] nmastnode	(i)pointer into field imastnode for contact tie i
 * @param [in] imastnode	field storing the nodes of the master surfaces
 * @param [in] islavnode	field storing the nodes of the slave surface
 * @param [in] islavsurf	islavsurf(1,i) slaveface i islavsurf(2,i) # integration points generated before looking at face i 
 * @param [in] imastsurf	index of masterface corresponding to integration point i
 * @param [in] pmastsurf 	field storing position and etal for integration points on master side
 * @param [in] itiefac		pointer into field islavsurf: (1,i) beginning slave_i (2,i) end of slave_i
 * @param [in] neq		(0) # of mechanical equations (1) sum of mechanical and thermal equations (2) neq(1+ # of single point contraints)
 * @param [in] nactdof		(i,j) actual degree of freedom of DOF i if node j
 * @param [in] co		coordinates of nodesf
 * @param [in] vold		displacement of nodesf
 * @param [in] iponoels		(i) pointer into field inoels...
 * @param [in] inoels		...which stores 1D&2D elements belonging to node 
 *                               (1,i)el. number (2,i) # nodes (3,i) pointer to next entry 
 * @param [in] mi
 * @param [in] gapmints		(i) gap between slave surface and master surface in integration point i
 * @param [out] gap		(i) \f$ g_i= \frac{1}{D_i} <g, \Psi_i> \f$ for node i on slave surface
 * @param [in] pslavsurf	field storing  position xil, etal and weight for integration point on slave side
 * @param [in] pslavdual	(1:4,i) dual shape functions for face i
 * @param [in] nintpoint	number of integration points
*/

void bdfill(int **irowbdp, int *jqbd,
        double **aubdp, double *bdd,int *nzsbd, int *ntie, int *ipkon, int *kon, 
        char *lakon, int *nslavnode, int *nmastnode, int *imastnode,
        int *islavnode, int *islavsurf, int *imastsurf, double *pmastsurf, 
        int *itiefac,char *tieset, int *neq, int *nactdof, double *co, double *vold,
	int *iponoels, int *inoels, int *mi, double *gapmints, double *gap,
        double* pslavsurf,double* pslavdual, int *nintpoint,double *slavnor,int *nk,
        int *nboun,int *ndirboun,int *nodeboun,double *xboun,
        int *nmpc,int *ipompc,int *nodempc,double *coefmpc,
        int *ikboun,int *ilboun,int *ikmpc,int *ilmpc,
        int *nslavspc,int *islavspc,int *nsspc,int *nslavmpc,int *islavmpc,int *nsmpc,
        int *nmastspc,int *imastspc,int *nmspc,int *nmastmpc,int *imastmpc,int *nmmpc,
	double **Bdp,double *Dd,int *jqb,int **irowbp, int *nzsbd2, double *dhinv){
		
  int i, j,jj, k,kk,l,m, ll,icounter,icounter2,idof1,idofs,idofm, nodesf, nodem, kflag,numb,
      *mast1=NULL,number, *irowbd=NULL,ifree,mt=mi[1]+1,istart, *iscontr=NULL, *imcontr=NULL,
      *idcontr1=NULL, *idcontr2=NULL,*igcontr=NULL, debug,intpointl, *irowb=NULL, *mast2=NULL,ifree2,
      ist,ist2,dir,node1,jslav,dim,id,dof,index,dirdep,dirind,index2,id2,dirdep2,dirind2,node2,idof2,
      idof3,islavk2;
  
  double contribution=0.0, *aubd=NULL, *contr=NULL, *dcontr=NULL, *gcontr=NULL,*anull=NULL, *Bd=NULL,
         gap2,*xs=NULL,*xm=NULL,*help=NULL,dm[3],n[3],c2,coefdep,coefind,coefdep2,coefind2,c3,detdh,dh[9]; 

       clock_t debut;
       clock_t fin;
  
  irowbd = *irowbdp; aubd=*aubdp; irowb= *irowbp; Bd=*Bdp;
  ifree = 1; // position in the fieds FORTRAN condition
  ifree2=1;
//  xs=NNEW(double,3);
//  xm=NNEW(double,3);
  help=NNEW(double,3);
  mast1=NNEW(int,*nzsbd);
  mast2=NNEW(int,*nk);
//  anull=NNEW(double,1); 
  
  debug=0;
  /* calculating the off-diagonal terms and storing them in aubd */

  /* meaning of the fields in FORTRAN notation:
     ipointer(i): points to an element in field aubd belonging to column i
                  aubd(ipointer(i)): value of that element
                  irowbd(ipointer(i)): row to which that element belongs
     mast1(ipointer(i)): points to another element in field aubd belonging
                         to column i, unless zero.
  */

  debut=clock();

  for (i=0;i<neq[1];i++){bdd[i]=0.0;}
  for (i=0;i<*nk;i++){Dd[i]=0.0;}
    

  for( i=0; i<*ntie; i++){    
    if(tieset[i*(81*3)+80]=='C'){ 
        for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
      	         gap[j]=0.0;
	}
	for (l= itiefac[2*i]; l<=itiefac[2*i+1];l++){
               	if(debug==1)printf("bdfill face %d intpoints %d %d \n",l, islavsurf[2*(l-1)+1],islavsurf[2*l+1]);
                intpointl=islavsurf[2*l+1]-islavsurf[2*(l-1)+1];
                if(intpointl>0){
		  contr=NNEW(double,9*9*intpointl/7+1);
		  iscontr=NNEW(int, 9*9*intpointl/7+1);
		  imcontr=NNEW(int, 9*9*intpointl/7+1);
		  dcontr=NNEW(double,9*9*intpointl/7+1);
		  idcontr1=NNEW(int, 9*9*intpointl/7+1);
		  idcontr2=NNEW(int, 9*9*intpointl/7+1);
		  gcontr=NNEW(double, 9*intpointl/7+1);
		  igcontr=NNEW(int, 9*intpointl/7+1);
		  FORTRAN(createbd,(&i,&l,ipkon,kon,lakon,co,vold,gapmints,
     			 islavsurf,imastsurf,pmastsurf,itiefac,contr,iscontr,imcontr,
  			dcontr,idcontr1,idcontr2,gcontr,igcontr,iponoels,inoels,mi,pslavsurf,pslavdual,
			nslavnode,islavnode,nmastnode,imastnode,&icounter,&icounter2));

		  for(j=0; j<icounter;j++){				
		    contribution=-contr[j];				
		    nodesf=islavnode[iscontr[j]-1];				
		    nodem=imastnode[imcontr[j]-1];				
		    if(debug==1)printf("\tbdfill: B_d nodesf %d nodem %d c %e \n",nodesf,nodem,contribution);                                
		    if ((contribution>1e-14 ||contribution<-1e-14)){ 				
		      insertas(&irowb, &mast2,&nodesf, &nodem,  &ifree2, nzsbd2,
		       			&contribution, &Bd);				
		    }				        
		    dof=0;				
		    for(ll=0; ll<3; ll++){	  				
		      idofs = nactdof[mt*(nodesf-1)+ll+1];	    				
		      idofm = nactdof[mt*(nodem-1)+ll+1];					
		      if(debug==1)printf("\t idofs %d idofm %d \n",idofs,idofm);	    				
		      if ((idofs>0)&&(idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){ ///insertion for active dofs  	      				
			insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			&contribution, &aubd);				        
		      }else if((idofs>0)&&(contribution>1e-18 ||contribution<-1e-18)){/// mpc on master node			                  
			for(jj=nmastmpc[2*(imcontr[j]-1)];jj<nmastmpc[2*(imcontr[j]-1)+1];jj++){                                           
			  ist=imastmpc[2*jj];                                           
			  dirdep=nodempc[3*(ist-1)+1];                                           
			  coefdep=coefmpc[ist-1];                                           
			  index=nodempc[3*(ist-1)+2];					   
			  if(ll==(dirdep-1)){                                           
			    while(index!=0){				   				            
			      node1=nodempc[3*(index-1)];                                             
			      dirind=nodempc[3*(index-1)+1];	                                    
			      c2=-coefmpc[index-1]*contribution/coefdep;	    			            
			      idofm = nactdof[mt*(node1-1)+(dirind-1)+1];                                            
			      if(idofm>0){	      				    
				insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			     &c2, &aubd);					    
			      }                                            
                              index=nodempc[3*(index-1)+2];	                                    
			    }	                                   
			  }					  
			}					
		      }else if((idofm>0)&&(contribution>1e-18 ||contribution<-1e-18)){///mpc on slave node			                  
			for(jj=nslavmpc[2*(iscontr[j]-1)];jj<nslavmpc[2*(iscontr[j]-1)+1];jj++){                                           
			  ist=islavmpc[2*(jj)];                                           
			  dirdep=nodempc[3*(ist-1)+1];                                           
			  coefdep=coefmpc[ist-1];                                           
			  index=nodempc[3*(ist-1)+2];					   
			  if(ll==(dirdep-1)){                                           
			    while(index!=0){				   				            
			      node1=nodempc[3*(index-1)];                                             
			      dirind=nodempc[3*(index-1)+1];	                                    
			      c2=-coefmpc[index-1]*contribution/coefdep;	    			            
			      idofs = nactdof[mt*(node1-1)+(dirind-1)+1];                                             
			      if(idofs>0){
				insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			     &c2, &aubd);
					    
			      }
                              index=nodempc[3*(index-1)+2];	                                    
			    }	                                   
			  }					  
			}										
		      }else if((idofs==0)&&(idofm==0)&&(contribution>1e-18 ||contribution<-1e-18)){///mpc on master and slave node			                  
			for(jj=nmastmpc[2*(imcontr[j]-1)];jj<nmastmpc[2*(imcontr[j]-1)+1];jj++){                                           
			  ist=imastmpc[2*jj];                                           
			  dirdep=nodempc[3*(ist-1)+1];                                           
			  coefdep=coefmpc[ist-1];                                           
			  index=nodempc[3*(ist-1)+2];					   
			  if(ll==(dirdep-1)){                                          
			    while(index!=0){				   				            
			      node1=nodempc[3*(index-1)];                                            
			      dirind=nodempc[3*(index-1)+1];	                                    
			      c2=-coefmpc[index-1]*contribution/coefdep;	    			            
			      idofm = nactdof[mt*(node1-1)+(dirind-1)+1];			    			                    
			      for(kk=nslavmpc[2*(iscontr[j]-1)];kk<nslavmpc[2*(iscontr[j]-1)+1];kk++){                                             
				ist2=islavmpc[2*kk];                                             
				dirdep2=nodempc[3*(ist2-1)+1];                                             
				coefdep2=coefmpc[ist2-1];                                             
				index2=nodempc[3*(ist2-1)+2];					     
				if(ll==(dirdep2-1)){                                              
				  while(index2!=0){				   				               
				    node2=nodempc[3*(index2-1)];                                               
				    dirind2=nodempc[3*(index2-1)+1];	                                       
				    c3=-coefmpc[index2-1]*c2/coefdep2;	    			               
				    idofs = nactdof[mt*(node2-1)+(dirind2-1)+1];                                               
				    if(idofs>0&&idofm>0){	      				       
				      insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			       &c3, &aubd);					       
				    }                                              
                                    index2=nodempc[3*(index2-1)+2];	                                      
				  }	                                     
				}					    
			      }					    				                                              
                              index=nodempc[3*(index-1)+2];	                                    
			    }	                                   
			  }					  
			}										
		      }  				
		    }        		  
		  }
		  for(j=0; j<icounter2;j++){			
		    contribution=dcontr[j];			
		    nodesf=islavnode[idcontr1[j]-1];		        
		    nodem=islavnode[idcontr2[j]-1];						
		    if(debug==1)printf("\tbdfill: face %d node %d %d dbb %e\n",l,nodesf,nodem,contribution);                       
		    if(nodesf==nodem){                          
		      Dd[nodesf-1]+=contribution;      			  
		      for(ll=0; ll<3; ll++){				
			idof1 = nactdof[mt*(nodesf-1)+ll+1];				
			if(debug==1) printf("\t idofs %d idofm %d c %e \n",idof1,idof1,contribution);				
			if (idof1>0){	  			          
			  bdd[idof1-1]+=contribution;				
			}else{			                  
			  for(jj=nslavmpc[2*(idcontr1[j]-1)];jj<nslavmpc[2*(idcontr1[j]-1)+1];jj++){                                           
			    ist=islavmpc[2*jj];                                           
			    dirdep=nodempc[3*(ist-1)+1];                                           
			    coefdep=coefmpc[ist-1];                                           
			    index=nodempc[3*(ist-1)+2];					   
			    if(ll==(dirdep-1)){                                           
			      while(index!=0){				   				            
				node1=nodempc[3*(index-1)];                                            
				dirind=nodempc[3*(index-1)+1];	                                    
				c2=-coefmpc[index-1]*contribution/coefdep;	    			            
				idofm = nactdof[mt*(node1-1)+(dirind-1)+1];	                                             
				index2=nodempc[3*(ist-1)+2];                                              
				while(index2!=0){					               
				  node1=nodempc[3*(index2-1)];                                                
				  dirind=nodempc[3*(index2-1)+1];	                                       
				  c3=-coefmpc[index2-1]*c2/coefdep;	    			               
				  idofs = nactdof[mt*(node1-1)+(dirind-1)+1];					       
				  if(idofs==idofm && idofs>0){	      				         
				    bdd[idofs-1]+=c3;					       
				  }else if(idofs>0 && idofm>0){						 
				    if(islavmpc[2*jj+1]==1){	      				           
				      insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			           &c3, &aubd);						       						 
				    }					       
				  } 
                                  index2=nodempc[3*(index2-1)+2];	                                      
				}					    				  
                                index=nodempc[3*(index-1)+2];					    
			      }	                                    
			    }					  
			  }				  				
			}        			  
		      }			
		    }else{  			        
		      if(contribution>1e-14 ||contribution<-1e-14){				  
			insertas(&irowb, &mast2,&nodesf, &nodem,  &ifree2, nzsbd2,
		       			&contribution, &Bd);
		      }			
		      for(ll=0; ll<3; ll++){	  							  
			idofs = nactdof[mt*(nodesf-1)+ll+1];	    							  
			idofm = nactdof[mt*(nodem-1)+ll+1];	    							  
			if ((idofs>0)&&(idofm>0)){ ///insertion for active dofs  	      							    
			  insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,		       			
				     &contribution, &aubd);			  
			}else if((idofs>0)&&(idofm==0)){/// mpc on slavenode1 node			                  			    
			  for(jj=nslavmpc[2*(idcontr2[j]-1)];jj<nslavmpc[2*(idcontr2[j]-1)+1];jj++){                                           			      
			    ist=islavmpc[2*(jj)];                                           			      
			    dirdep=nodempc[3*(ist-1)+1];                                           			      
			    coefdep=coefmpc[ist-1];                                           			      
			    index=nodempc[3*(ist-1)+2];					   			      
			    if(ll==(dirdep-1)){                                           				
			      while(index!=0){				   				            				  
				node1=nodempc[3*(index-1)];                                             				  
				dirind=nodempc[3*(index-1)+1];	                                    				  
				c2=-coefmpc[index-1]*contribution/coefdep;	    			            				  
				idofm = nactdof[mt*(node1-1)+(dirind-1)+1];                                            				  
				if(idofm>0){	      				    				    
				  insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			     &c2, &aubd);					    				  
				}                                            
                                 index=nodempc[3*(index-1)+2];	                                    				
			      }	                                   			      
			    }					  			    
			  }								  
			}else if((idofm>0)&&(idofs==0)){///mpc on slavenode2 node			                  
			  for(jj=nslavmpc[2*(idcontr1[j]-1)];jj<nslavmpc[2*(idcontr1[j]-1)+1];jj++){                                           
			    ist=islavmpc[2*(jj)];                                          
			    dirdep=nodempc[3*(ist-1)+1];                                           
			    coefdep=coefmpc[ist-1];                                          
			    index=nodempc[3*(ist-1)+2];					   
			    if(ll==(dirdep-1)){                                           
			      while(index!=0){				   				            
				node1=nodempc[3*(index-1)];                                             
				dirind=nodempc[3*(index-1)+1];	                                    
				c2=-coefmpc[index-1]*contribution/coefdep;	    			            
				idofs = nactdof[mt*(node1-1)+(dirind-1)+1];                                            
				if(idofs>0){	      				    
				  insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			     &c2, &aubd);					    
				}                                            
                                index=nodempc[3*(index-1)+2];	                                    
			      }	                                   
			    }					  
			  }										
			}else if((idofs==0)&&(idofm==0)){///mpc on both slave node (possible?)			                  
			  for(jj=nslavmpc[2*(idcontr2[j]-1)];jj<nslavmpc[2*(idcontr2[j]-1)+1];jj++){                                           
			    ist=islavmpc[2*(jj)];                                           
			    dirdep=nodempc[3*(ist-1)+1];                                           
			    coefdep=coefmpc[ist-1];                                           
			    index=nodempc[3*(ist-1)+2];					   
			    if(ll==(dirdep-1)){                                           
			      while(index!=0){				   				            
				node1=nodempc[3*(index-1)];                                             
				dirind=nodempc[3*(index-1)+1];	                                    
				c2=-coefmpc[index-1]*contribution/coefdep;	    			            
				idofm = nactdof[mt*(node1-1)+(dirind-1)+1];				    			                    
				for(kk=nslavmpc[2*(idcontr1[j]-1)];kk<nslavmpc[2*(idcontr1[j]-1)+1];kk++){                                             
				  ist2=islavmpc[2*kk];                                             
				  dirdep2=nodempc[3*(ist2-1)+1];                                             
				  coefdep2=coefmpc[ist2-1];                                             
				  index2=nodempc[3*(ist2-1)+2];					     
				  if(ll==(dirdep2-1)){                                              
				    while(index2!=0){				   				               
				      node2=nodempc[3*(index2-1)];                                                
				      dirind2=nodempc[3*(index2-1)+1];	                                       
				      c3=-coefmpc[index2-1]*c2/coefdep2;	    			               
				      idofs = nactdof[mt*(node2-1)+(dirind2-1)+1];                                               
				      if(idofs>0&&idofm>0){	      				       
					insertas(&irowbd, &mast1, &idofs, &idofm, &ifree, nzsbd,
		       			       &c3, &aubd);					       
				      }                                              
                                      index2=nodempc[3*(index2-1)+2];	                                      
				    }	                                     
				  }					    
				}					    				  
                                index=nodempc[3*(index-1)+2];	                                    
			      }	                                   
			    }					  
			  }										
			}  				  
		      }				
		    }		        
		  }  		  			
		  for(j=0; j<sqrt(icounter2);j++){					
		    contribution=gcontr[j];			
		    ll=igcontr[j]-1;			
		    gap[ll]+=contribution;			
		    if(debug==1)printf("nodes %d c %e gap %e  \n",islavnode[ll],contribution,gap[ll]);						
		  }
		  free(contr);
		  free(iscontr);
		  free(imcontr);
		  free(dcontr);
		  free(idcontr1);
		  free(idcontr2);
		  free(gcontr);
		  free(igcontr);              
		}
	}
	for(j=nslavnode[i]; j<nslavnode[i+1]; j++){
	  nodesf = islavnode[j];
	  /* calculating the gap at the slave nodesf */
	  contribution=Dd[nodesf-1];
	  gap[j]=gap[j]/contribution; 	  
	}    
    }
  }
  fin=clock();
  printf("bdfill_createbd: %f s \n",((double) (fin-debut))/CLOCKS_PER_SEC);
  
/** Sort aubd **/
  *nzsbd=ifree-1;
  *nzsbd2=ifree2-1;
   printf("bdfill: sumicounter %i \n",(ifree-1));

  /* Sort mast1, irowbd and aubd; 
     Outcome: the values in field aubd are sorted, column by
     column; no sorting is done within the columns */

  debut=clock();
  
  dim=neq[1];
  matrixsort(aubd,mast1,irowbd,jqbd,nzsbd,&dim); 
  /*Calulation ot the real contribution*/
  
 icounter=0;
 
 for (i=0;i<neq[1];i++){
   if(jqbd[i]!=jqbd[i+1]){
     irowbd[icounter]=irowbd[jqbd[i]-1];
     aubd[icounter]=aubd[jqbd[i]-1];
     icounter++;
     istart=icounter;
     for (j=jqbd[i];j<jqbd[i+1]-1;j++){
       if (irowbd[j]==irowbd[icounter-1]){
	 aubd[icounter-1]+=aubd[j];   
       }else{
	 irowbd[icounter]=irowbd[j];
	 aubd[icounter]=aubd[j];
	 icounter++;
       }
   }
   }else{ istart=icounter+1;}
  
  jqbd[i]=istart;
 }
 jqbd[neq[1]]=icounter+1; 
   printf("bdfill: size aubd %i \n",icounter);
 *nzsbd=icounter;
 RENEW(irowbd,int,*nzsbd);
 RENEW(aubd,double,*nzsbd);
 
   number=6;
 free(mast1);
  
 /** sort Bd **/
 debut=clock();
 dim=*nk;
 matrixsort(Bd,mast2,irowb,jqb,nzsbd2,&dim);
     
 icounter=0;
 for (i=0;i<*nk;i++){
   if(jqb[i]!=jqb[i+1]){
     irowb[icounter]=irowb[jqb[i]-1];
     Bd[icounter]=Bd[jqb[i]-1];
     icounter++;
     istart=icounter;
     for (j=jqb[i];j<jqb[i+1]-1;j++){
       if (irowb[j]==irowb[icounter-1]){
	 Bd[icounter-1]+=Bd[j];   
       }else{
	 irowb[icounter]=irowb[j];
	 Bd[icounter]=Bd[j];
	 icounter++;
       }
   }
   }else{ istart=icounter+1;}
  
  jqb[i]=istart;
 }
 jqb[*nk]=icounter+1; 
   printf("bdfill: size Bd %i \n",icounter);
 *nzsbd2=icounter;
 RENEW(irowb,int,*nzsbd2);
 RENEW(Bd,double,*nzsbd2);

free(mast2);


/** handle SPC's on master nodes **/
  double ndm=0.0;
  for( i=0; i<*ntie; i++){            
    if(tieset[i*(81*3)+80]=='C'){             
      for(j=nmastnode[i]; j<nmastnode[i+1]; j++){	   
	nodem = imastnode[j];	   
	dm[0]=0.0;dm[1]=0.0;dm[2]=0.0;	   
	for(jj=nmastspc[2*(j)];jj<nmastspc[2*(j)+1];jj++){	     
	  ist=imastspc[2*jj];             
	  dir=ndirboun[ist-1];	     
	  node1=nodeboun[ist-1];	     
	  dm[dir-1]=xboun[ist-1];	     
	  if(debug==1){printf("jj %d ist %d dir %d node %d\n",jj,ist,dir,node1);}	   
	}
	for(jj=nmastmpc[2*(j)];jj<nmastmpc[2*(j)+1];jj++){ 	    
	  if(imastmpc[2*jj+1]==3){                   	     
	    ist=imastmpc[2*jj];                                           	     
	    dirdep=nodempc[3*(ist-1)+1];                                           	     
	    coefdep=coefmpc[ist-1];                                           	     
	    index=nodempc[3*(ist-1)+2];                             	     
	    while(index!=0){				   				            		 
	      node1=nodempc[3*(index-1)];                                             		 
	      dirind=nodempc[3*(index-1)+1];	                                    		 
	      c2=-coefmpc[index-1]*contribution/coefdep;	    			            		 
	      idofm = nactdof[mt*(node1-1)+(dirind-1)+1];		 
	      if(idofm==0){		   
		idofm=8*(node1-1)+dirind;		   
		FORTRAN(nident,(ikboun,&idofm,nboun,&id));		   
		if(id>0 &&ikboun[id-1]==idofm){		     
		  dm[dirdep-1]=dm[dirdep-1]-coefmpc[index-1]*xboun[ilboun[id-1]]/coefdep;		   
		}		 
	      }
	      index=nodempc[3*(index-1)+2];	                                    	       
	    }	     
	  }	    	   
	}
        ndm=sqrt(dm[0]*dm[0]+dm[1]*dm[1]+dm[2]*dm[2]);	    
	if(ndm>1.e-16){            
	  for (jj=jqb[nodem-1]-1;jj<jqb[nodem]-1;jj++){	      
	    nodesf=irowb[jj];             
	    for(k=0;k<3;k++){              
	      help[k]=1/(Dd[nodesf-1])*Bd[jj]*dm[k];             
	    }        	      
	    dim=nslavnode[i+1]-nslavnode[i];	          	      
	    FORTRAN(nident,(&islavnode[nslavnode[i]], &nodesf,&dim, &id));	 	      
	    if(id>0 && islavnode[nslavnode[i]+id-1]==nodesf){	   		
	      jslav=nslavnode[i]+id-1;	 	      
	    }else{
	      printf("\n bdfill: handle SPC's on master side. something went wrong...\n");
	      FORTRAN(stop,());	 	      
	    }                   
            gap2=help[0]*slavnor[3*jslav]
	          +help[1]*slavnor[3*jslav+1]
	          +help[2]*slavnor[3*jslav+2];
	    gap[jslav]=gap[jslav]-gap2;  	    
	  }	    
	}	
      }	
    }
  }     
  printf("bdfill: SPC's on master side handled \n");	   
  /** calculate inverse of D_d in general setting with SPC's **/
  
  for (i=0;i<*ntie;i++){      
    if(tieset[i*(81*3)+80]=='C'){	
      for(j=nslavnode[i];j<nslavnode[i+1];j++){	    
	nodesf=islavnode[j];	    
	idof1=nactdof[mt*nodesf-3]-1;	    
	idof2=nactdof[mt*nodesf-2]-1;	    
	idof3=nactdof[mt*nodesf-1]-1;	    	    
	if(idof1>-1) {dh[0]=bdd[idof1];}else{dh[0]=0.0;}	    
	if(idof2>-1) {dh[4]=bdd[idof2];}else{dh[4]=0.0;}	    
	if(idof3>-1) {dh[8]=bdd[idof3];}else{dh[8]=0.0;}	    
	dh[1]=0.0;dh[2]=0.0;dh[3]=0.0;dh[5]=0.0;dh[6]=0.0;dh[7]=0.0;
	if(idof1>-1 && bdd[idof1]>1.e-15) {dhinv[j*9+0]=1.0/dh[0];}else{dhinv[j*9+0]=0.0;}
	if(idof2>-1 && bdd[idof2]>1.e-15) {dhinv[j*9+4]=1.0/dh[4];}else{dhinv[j*9+4]=0.0;}
	if(idof3>-1 && bdd[idof3]>1.e-15) {dhinv[j*9+8]=1.0/dh[8];}else{dhinv[j*9+8]=0.0;}	    
	dhinv[j*9+1]=0.0;dhinv[j*9+2]=0.0;dhinv[j*9+3]=0.0;dhinv[j*9+5]=0.0;dhinv[j*9+6]=0.0;dhinv[j*9+7]=0.0;
	for(jj=0;jj<3;jj++){	     
	  idofs=nactdof[mt*nodesf-3+jj]-1;	     
	  if(idofs>-1){	      
	    if(jqbd[idofs+1]-jqbd[idofs]>0){	       
	      for(kk=jqbd[idofs]-1;kk<jqbd[idofs+1]-1;kk++){	        
		if(irowbd[kk]-1==idof1)dirind = 1;
	        if(irowbd[kk]-1==idof2)dirind = 2;
	        if(irowbd[kk]-1==idof3)dirind = 3;
		if(irowbd[kk]-1==idof1 ||irowbd[kk]-1==idof2||irowbd[kk]-1==idof3){
	         dh[3*jj+dirind-1]=aubd[kk];
		}	       
	      }	      
	    } 	     
	  }	    
	}	    
	if(idof1==-1 && idof2>-1 && idof3>-1){ 	     
	  detdh=dh[4]*dh[8]-dh[5]*dh[7];	     
	  dhinv[j*9+4]=1/detdh*dh[8];	     
	  dhinv[j*9+8]=1/detdh*dh[4]; 	     
	  dhinv[j*9+5]=-1/detdh*dh[5];	     
	  dhinv[j*9+7]=-1/detdh*dh[7]; 	     	    	    
	}else if(idof2==-1 && idof1>-1 && idof3>-1){ 	     
	  detdh=dh[0]*dh[8]-dh[2]*dh[6];	     
	  dhinv[j*9+0]=1/detdh*dh[8];	     
	  dhinv[j*9+8]=1/detdh*dh[0]; 	     
	  dhinv[j*9+2]=-1/detdh*dh[2];	     
	  dhinv[j*9+6]=-1/detdh*dh[6];	      	    
	}else if(idof3==-1 && idof2>-1 && idof1>-1){ 	      	     
	  detdh=dh[4]*dh[0]-dh[1]*dh[3];	     
	  dhinv[j*9+4]=1/detdh*dh[0];	     
	  dhinv[j*9+0]=1/detdh*dh[4]; 	     
	  dhinv[j*9+1]=-1/detdh*dh[1];	     
	  dhinv[j*9+3]=-1/detdh*dh[3];	    	    
	}	    
	if(debug==1){	    
	  printf("bdfill: nodes %d\n",nodesf);	    
	  printf("\t %e %e %e \t %e %e %e\n",dh[0],dh[1],dh[2],dhinv[j*9+0],dhinv[j*9+1],dhinv[j*9+2]);	    
	  printf(" D=\t %e %e %e Di=\t %e %e %e\n",dh[3],dh[4],dh[5],dhinv[j*9+3],dhinv[j*9+4],dhinv[j*9+5]);	    
	  printf("\t %e %e %e \t %e %e %e\n",dh[6],dh[7],dh[8],dhinv[j*9+6],dhinv[j*9+7],dhinv[j*9+8]);	    	    
	}	
      }      
    }    
  }
  printf("bdfill: inverse of D calculated\n");	
  

  *irowbdp = irowbd; *aubdp=aubd; *irowbp = irowb; *Bdp=Bd;
  free(help);
  return;
}


