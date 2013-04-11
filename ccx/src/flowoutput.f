!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!     
!     construction of the B matrix      
!     
      subroutine flowoutput(itg,ieg,ntg,ntm,bc,lakon,ntmat_,
     &     v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
     &     ielmat,prop,ielprop,nactdog,nacteq,iin,physcon,
     &     camt,camf,camp,uamt,uamf,uamp,rhcon,nrhcon,
     &     vold,jobnamef,set,istartset,iendset,ialset,nset,mi)
!     
      implicit none
!     
      logical identity
      character*10 tim
      character*8 lakon(*)
      character*10 date
      character*132 jobnamef(*)
      character*30 jobname(30)
      character*81 set(*)
!     
      integer itg(*),ieg(*),ntg,nflow,ielmat(*),
     &     i,j,k,nrhcon(*),node,imat,ntmat_,ntm,numf,
     &     node1,node2,nshcon(*),nelem,index,ipkon(*),
     &     kon(*),idof,nodem,idirf(5),ieq,nactdog(0:3,*),
     &     nacteq(0:3,*),ielprop(*),nodef(5),iin,kflag,
     &     istartset(*),iendset(*),ialset(*),nset,count_set,
     &     istart,iend,mi(2)
!     
      real*8 bc(ntm),cp,physcon(3),r,dvi,rho,
     &     gastemp,v(0:mi(2),*),shcon(0:3,ntmat_,*),co(3,*),
     &     prop(*),tg1,tg2,dtime,ttime,time,g(3),
     &     xflow,tvar(2),f,df(5),camt(*),camf(*),camp(*),
     &     rhcon(0:1,ntmat_,*),ts1,ts2,pi,vold(0:mi(2),*),
     &     tim2,uamt,uamf,uamp,maxrt,maxrf,maxrp
!          
!     After convergence, this subroutine runs once again the converged 
!     results through the element dedicated subroutines and creates a 
!     formatted file containing for each elements miscelaneous information 
!     such as mach and reynolds numbers, cds...
!
!      write(*,*) 'flowoutput'
      kflag=3
!     
      tvar(1)=time
      tvar(2)=ttime+dtime
!     
      pi=4.d0*datan(1.d0)
!
      open(UNIT=1,FILE='fort.10_ccx',STATUS='REPLACE',RECL=300)
!     
      write(1,*) 'Air and Oil System calculation with CalculiX:'
      
      call  DATE_AND_TIME(date,tim)
      write(1,*) 'created the ',date(7:8),'/',date(5:6),'/',date(1:4),
     &' at ',tim(1:2),':',tim(3:4)
!
      i=1
      if(jobnamef(i).ne.'''') then 
         jobname(i)=jobnamef(i)
         i=1+1
      endif
!
      write(1,*) 'Input file name: ',jobname(1:i-1)
!     
      Write(1,*) ''
      Write(1,*) 'System output and convergence information'
      write(1,*) '========================================='
      Write(1,*) ''
!
      call CPU_TIME(tim2) 
      write(1,1) 'CPU time: ',tim2,' s'
 1    FORMAT(1X,A,F7.5,A)
 2    FORMAT(1X,A,G12.6)
 3    FORMAT(1X,A,G12.6,A,I6)
      write(1,*) 'Iteration number: ',iin
      write(1,2) '      largest increment of gas temperature  = ',uamt
      write(1,3) '      largest correction to gas temperature = ',
     &     camt(1),' in node',int(camt(2))
      write(1,2) '      largest increment of gas massflow     = ',uamf
      write(1,3) '      largest correction to gas massflow    = ',
     &     camf(1),' in node',int(camf(2))
      write(1,2) '      largest increment of gas pressure     = ',uamp
      write(1,3) '      largest correction to gas pressure    = ',
     &     camp(1),' in node',int(camp(2))
!
      maxrt=0
      maxrf=0
      maxrp=0
!
!
!     set names in set(*) field
!
      i=1
      numf=0
!
C      do
C         if(set(i)(1:5).eq.'NLINE') then          
C            numf=numf+1
C           write(*,*) 'numf',numf
C            i=i+1
C         else
C            exit
C         endif
C      enddo 
!
      do i=1,ntg
         node=itg(i)
         do j=0,2
            if(nactdog(j,node).eq.0) cycle
            idof=nactdog(j,node)
!     temperature residual
            if(j.eq.0) then
               maxrt=max(maxrt,bc(idof))
!
!     massflow residual
            elseif(j.eq.1) then
               maxrf=max(maxrf,bc(idof))
!
!     pressure residual
            else
               maxrp=max(maxrp,bc(idof))
!
            endif
         enddo
      enddo

!
      write(1,*) ''  
      write(1,2) '      largest energy equation residual      = ',maxrt
      write(1,2) '      largest continuity equation residual  = ',maxrf
      write(1,2) '      largest momentum equation residual    = ',maxrp
!
      write(1,*) ''
      write(1,*) ''
      write(1,*)'Steady state solution detailled for each lines'
      write(1,*) '============================================='
      write(1,*) ''
      write(1,*) ''  
!    
      count_set=0
      do i=1, nset
C         write(*,*) 'set',i,set(i)(1:20)
C         write(*,*) '    istartset',istartset(i)
C         write(*,*) '    iendset',iendset(i)
C         write(*,*) '    ialset(istartset)',ialset(istartset(i))
C         write(*,*) '    ialset(iendset)',ialset(iendset(i))
         if((set(i)(1:4).eq.'NLIN')) then
            count_set=count_set+1
C         else if (set(i)(1:4).eq.'NALL') then
         endif

      enddo
C      write(*,*) 'count_set',count_set
      do i=1,nflow
C         write(*,*)'nflow',nflow
         nelem=ieg(i)
C         write(*,*) 'nelem', nelem, lakon(nelem)(1:6)
!         
         index=ipkon(nelem)
         node1=kon(index+1)
         nodem=kon(index+2)
         node2=kon(index+3)
!
         if(node1.eq.0) then
            tg1=v(0,node2)
            tg2=tg1
            ts1=v(3,node2)
            ts2=Ts1
         elseif(node2.eq.0) then
            tg1=v(0,node1)
            tg2=tg1
            ts1=v(3,node1)
            ts2=ts1
         else
            tg1=v(0,node1)
            tg2=v(0,node2)
            ts1=v(3,node1)
            ts2=v(3,node2)
         endif
         gastemp=(ts1+ts2)/2.d0
!     
!         write(*,*) 'nelem', nelem, set(nelem)(1:6)
         imat=ielmat(nelem)
!     
         call materialdata_tg(imat,ntmat_,gastemp,shcon,nshcon,cp,r,dvi,
     &        rhcon,nrhcon,rho)
!     
         if (nacteq(2,nodem).ne.0) then
            ieq=nacteq (2,nodem)
!
            do j = count_set+1, nset
               if(set(j)(1:4).ne.'NALL') then
                  if(ialset(istartset(j)).eq.ialset(iendset(j))) then
                     if(ialset(istartset(j)).eq.nelem) then
                        numf=j
C                        write(*,*) nelem, set(numf)(1:20)
                     endif
                  elseif(ialset(istartset(j)).ne.ialset(iendset(j)))then
                     istart=istartset(j)
C     write(*,*) 'istart',istart
                     iend=iendset(j)                  
C     write(*,*) 'iend',iend
                     do k=istart, iend
C                     write(*,*) 'nelem',nelem
C     write(*,*) 'k',k
                        if(ialset(k).eq.nelem) then
                           numf=j
C                           write(*,*) nelem, set(numf)(1:20)
                       endif
                     enddo
                     
                  endif
               endif
            enddo
            xflow=v(1,nodem)
!            write(*,*) 'set(nelem+numf)(1:20)',set(nelem+numf)(1:20)
            call flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &           nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           cp,r,rho,physcon,g,co,dvi,numf,vold,set,shcon,
     &           nshcon,rhcon,nrhcon,ntmat_,mi)
         endif
      enddo
!     
      close(1)

      end
