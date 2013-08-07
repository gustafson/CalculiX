!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2013 Guido Dhondt
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
      subroutine equationremove(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nk,ipompc,nodempc,nmpc,mpcfree,ikmpc,ilmpc,
     &  labmpc,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc,itextpart)
!
!     removing equations
!
      implicit none
!
      character*1 inpc(*)
      character*20 labmpc(*)
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),itextpart,nset,istat,n,
     &  i,j,k,l,ibounstart,ibounend,key,nk,ipompc(*),nodempc(3,*),
     &  nmpc,mpcfree,ikmpc(*),ilmpc(*),nmpcold,id,idof,index1,ipos,m,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*)
!
!     all MPC's are removed
!
      if(textpart(itextpart)(1:10).eq.'REMOVE=ALL') then
         nmpc=0
         return
      endif
!
!     only some MPC's are removed
!
      nmpcold=nmpc
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(2)(1:10),'(i10)',iostat=istat) ibounstart
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
!     
         if(textpart(3)(1:1).eq.' ') then
            ibounend=ibounstart
         else
            read(textpart(3)(1:10),'(i10)',iostat=istat) ibounend
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         endif
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) k
         if(istat.eq.0) then
            if((k.gt.nk).or.(k.le.0)) then
               write(*,*) '*ERROR reading *BOUNDARY:'
               write(*,*) '       node ',k,' is not defined'
               stop
            endif
!
!           removing the MPC from the fields ikmpc and ilmpc
!
            do l=1,3
               idof=8*(k-1)+l
               call nident(ikmpc,idof,nmpc,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                     index1=ipompc(ilmpc(id))
                     do
                        if(index1.eq.0) then
                           nodempc(3,index1)=mpcfree
                           mpcfree=ipompc(ilmpc(id))
                           ipompc(ilmpc(id))=0
                           do m=id,nmpc-1
                              ikmpc(m)=ikmpc(m+1)
                              ilmpc(m)=ilmpc(m+1)
                           enddo
                           ikmpc(nmpc)=0
                           ilmpc(nmpc)=0
                           nmpc=nmpc-1
                           exit
                        endif
                        index1=nodempc(3,index1)
                     enddo
                  endif
               endif
            enddo
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            do i=1,nset
               if(set(i).eq.noset) exit
            enddo
            if(i.gt.nset) then
               noset(ipos:ipos)=' '
               write(*,*) '*ERROR reading *BOUNDARY: node set ',noset
               write(*,*) '  has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  k=ialset(j)
!
!                 removing the MPC from the fields ikmpc and ilmpc
!
                  do l=1,3
                     idof=8*(k-1)+l
                     call nident(ikmpc,idof,nmpc,id)
                     if(id.gt.0) then
                        if(ikmpc(id).eq.idof) then
                           index1=ipompc(ilmpc(id))
                           do
                              if(index1.eq.0) then
                                 nodempc(3,index1)=mpcfree
                                 mpcfree=ipompc(ilmpc(id))
                                 ipompc(ilmpc(id))=0
                                 do m=id,nmpc-1
                                    ikmpc(m)=ikmpc(m+1)
                                    ilmpc(m)=ilmpc(m+1)
                                 enddo
                                 ikmpc(nmpc)=0
                                 ilmpc(nmpc)=0
                                 nmpc=nmpc-1
                                 exit
                              endif
                              index1=nodempc(3,index1)
                           enddo
                        endif
                     endif
                  enddo
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
!
!                    removing the MPC from the fields ikmpc and ilmpc
!
                     do l=1,3
                        idof=8*(k-1)+l
                        call nident(ikmpc,idof,nmpc,id)
                        if(id.gt.0) then
                           if(ikmpc(id).eq.idof) then
                              index1=ipompc(ilmpc(id))
                              do
                                 if(index1.eq.0) then
                                    nodempc(3,index1)=mpcfree
                                    mpcfree=ipompc(ilmpc(id))
                                    ipompc(ilmpc(id))=0
                                    do m=id,nmpc-1
                                       ikmpc(m)=ikmpc(m+1)
                                       ilmpc(m)=ilmpc(m+1)
                                    enddo
                                    ikmpc(nmpc)=0
                                    ilmpc(nmpc)=0
                                    nmpc=nmpc-1
                                    exit
                                 endif
                                 index1=nodempc(3,index1)
                              enddo
                           endif
                        endif
                     enddo
                  enddo
               endif
            enddo
         endif
      enddo
!
!     removing the deleted MPC's from ipompc and labmpc
!
      k=0
      do j=1,nmpcold
         if(ipompc(j).ne.0) then
            k=k+1
            ipompc(k)=ipompc(j)
            labmpc(k)=labmpc(j)
            index1=ipompc(j)
            idof=8*(nodempc(1,index1)-1)+nodempc(2,index1)
            call nident(ikmpc,idof,nmpc,id)
            if(id.eq.0) then
               write(*,*) '*ERROR reading *BOUNDARY'
               stop
            elseif(ikmpc(id).ne.idof) then
               write(*,*) '*ERROR reading *BOUNDARY'
               stop
            endif
            ilmpc(id)=k
         endif
      enddo
!
      return
      end

