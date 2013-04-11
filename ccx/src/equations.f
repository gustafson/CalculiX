!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2011 Guido Dhondt
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
      subroutine equations(inpc,textpart,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,nk,co,trab,inotr,ntrans,ikmpc,ilmpc,
     &  labmpc,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc)
!
!     reading the input deck: *EQUATION
!
      implicit none
!
      character*1 inpc(*)
      character*20 labmpc(*)
      character*132 textpart(16)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,istep,istat,
     &  n,i,j,ii,key,nterm,number,nk,inotr(2,*),ntrans,node,ndir,
     &  mpcfreeold,ikmpc(*),ilmpc(*),id,idof,itr,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*)
!
      real*8 coefmpc(*),co(3,*),trab(7,*),a(3,3),x
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in equations: *EQUATION should be placed'
         write(*,*) '  before all step definitions'
         stop
      endif
!
      do i=2,n
         write(*,*) 
     &        '*WARNING in equations: parameter not recognized:'
         write(*,*) '         ',
     &        textpart(i)(1:index(textpart(i),' ')-1)
         call inputwarning(inpc,ipoinpc,iline)
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
         read(textpart(1)(1:10),'(i10)',iostat=istat) nterm
!
         nmpc=nmpc+1
         if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR in equations: increase nmpc_'
            stop
         endif
!
         labmpc(nmpc)='                    '
         ipompc(nmpc)=mpcfree
         ii=0
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) then
               write(*,*) '*ERROR in equations: mpc definition ',nmpc
               write(*,*) '  is not complete. '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
!
            do i=1,n/3
!
               read(textpart((i-1)*3+1)(1:10),'(i10)',iostat=istat) node
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               if((node.gt.nk).or.(node.le.0)) then
                  write(*,*) '*ERROR in equations:'
                  write(*,*) '       node ',node,' is not defined'
                  stop
               endif
!
               read(textpart((i-1)*3+2)(1:10),'(i10)',iostat=istat) ndir 
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               if(ndir.le.3) then
               elseif(ndir.eq.4) then
                  ndir=5
               elseif(ndir.eq.5) then
                  ndir=6
               elseif(ndir.eq.6) then
                  ndir=7
               elseif(ndir.eq.8) then
                  ndir=4
               elseif(ndir.eq.11) then
                  ndir=0
               else
                  write(*,*) '*ERROR in equations:'
                  write(*,*) '       direction',ndir,' is not defined'
                  stop
               endif
c               if(ndir.eq.11) ndir=0
c               if(ndir.gt.3) then
c                  write(*,*) '*ERROR in equations:'
c                  write(*,*) '       direction',ndir,' is not defined'
c                  stop
c               endif
!
               read(textpart((i-1)*3+3)(1:20),'(f20.0)',iostat=istat) x
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!
!              check whether the node is transformed
!
               if(ntrans.le.0) then
                  itr=0
               elseif(inotr(1,node).eq.0) then
                  itr=0
               else
                  itr=inotr(1,node)
               endif
!
               if((itr.eq.0).or.(ndir.eq.0).or.(ndir.eq.4)) then
                  nodempc(1,mpcfree)=node
                  nodempc(2,mpcfree)=ndir
                  coefmpc(mpcfree)=x
!
!                 updating ikmpc and ilmpc
!
                  if(ii.eq.0) then
                     idof=8*(node-1)+ndir
                     call nident(ikmpc,idof,nmpc-1,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           write(*,100)
     &                   (ikmpc(id))/8+1,ikmpc(id)-8*((ikmpc(id))/8)
                           stop
                        endif
                     endif
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                  endif
!
                  mpcfreeold=mpcfree
                  mpcfree=nodempc(3,mpcfree)
                  if(mpcfree.eq.0) then
                     write(*,*) '*ERROR in equations: increase nmpc_'
                     stop
                  endif
               else
                  call transformatrix(trab(1,inotr(1,node)),
     &                 co(1,node),a)
!
                  number=ndir-1
                  if(ii.eq.0) then
!
!                    determining which direction to use for the
!                    dependent side: should not occur on the dependent
!                    side in another MPC and should have a nonzero
!                    coefficient
!
                     do j=1,3
                        number=number+1
                        if(number.gt.3) number=1
                        idof=8*(node-1)+number
                        call nident(ikmpc,idof,nmpc-1,id)
                        if(id.gt.0) then
                           if(ikmpc(id).eq.idof) then
                              cycle
                           endif
                        endif
                        if(dabs(a(number,ndir)).lt.1.d-5) cycle
                        exit
                     enddo
                     if(j.gt.3) then
                        write(*,*) '*ERROR in equations: SPC in node'
                        write(*,*) node,' in transformed coordinates'
                        write(*,*) ' cannot be converted in MPC: all'
                        write(*,*) ' DOFs in the node are used as'
                        write(*,*) ' dependent nodes in other MPCs'
                        stop
                     endif
                     number=number-1
!
!                    updating ikmpc and ilmpc
!
                     do j=nmpc,id+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(id+1)=idof
                     ilmpc(id+1)=nmpc
                  endif
!
                  do j=1,3
                     number=number+1
                     if(number.gt.3) number=1
                     if(dabs(a(number,ndir)).lt.1.d-5) cycle
                     nodempc(1,mpcfree)=node
                     nodempc(2,mpcfree)=number
                     coefmpc(mpcfree)=x*a(number,ndir)
                     mpcfreeold=mpcfree
                     mpcfree=nodempc(3,mpcfree)
                     if(mpcfree.eq.0) then
                        write(*,*) '*ERROR in equations: increase nmpc_'
                        stop
                     endif
                  enddo
               endif
!
               ii=ii+1
            enddo
!
            if(ii.eq.nterm) then
               nodempc(3,mpcfreeold)=0
               exit
            endif
         enddo
      enddo
!
 100  format(/,'*ERROR in equations: the DOF corresponding to',
     &           /,'node ',i7,' in direction',i5,' is detected on',
     &           /,'the dependent side of two different MPC''s') 
      return
      end

