!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine surfaces(inpc,textpart,set,istartset,iendset,ialset,
     &  nset,nset_,nalset,nalset_,nk,ne,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,lakon,ipoinpc)
!
!     reading the input deck: *SURFACE
!
      implicit none
!
      character*1 type,inpc(*)
      character*8 lakon(*)
      character*20 label
      character*81 set(*),noset,elset
      character*132 textpart(16)
!
      integer nset,nset_,nalset,nalset_,istep,istat,n,key,i,nk,ne,
     &  j,istartset(*),iendset(*),ialset(*),ipos,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),iside,l,k,kstart,kend,ipoinpc(0:*)
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in surfaces: *SURFACE should be placed'
         write(*,*) '  before all step definitions'
         stop
      endif
!
      nset=nset+1
      if(nset.gt.nset_) then
         write(*,*) '*ERROR in surfaces: increase nset_'
         stop
      endif
!
      kstart=0
      kend=0
!
!     reading the name of the set and the type of set
!
      do i=1,81
         set(nset)(i:i)=' '
      enddo
!
      type='T'
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=')
     &        then
            set(nset)(1:80)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) '*ERROR in surfaces: surface name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       surface name:',textpart(i)(1:132)
               stop
            endif
         elseif(textpart(i)(1:5).eq.'TYPE=') then
            if(textpart(i)(6:12).eq.'ELEMENT') then
               type='T'
            elseif(textpart(i)(6:9).eq.'NODE') then
               type='S'
            else
               write(*,*) '*ERROR in surfaces: unknown surface type'
               stop
            endif
         endif
      enddo
!
      ipos=index(set(nset),' ')
      if(ipos.eq.1) then
         write(*,*) '*ERROR in surfaces: no name specified'
         stop
      endif
      set(nset)(ipos:ipos)=type
!
      istartset(nset)=nalset+1
      iendset(nset)=nalset
!
      if(type.eq.'S') then
!
!        node surface
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            if(nalset+1.gt.nalset_) then
               write(*,*) '*ERROR in surfaces: increase nalset_'
               stop
            endif
!     
            read(textpart(1)(1:10),'(i10)',iostat=istat)ialset(nalset+1)
            if(istat.gt.0) then
               noset=textpart(1)(1:80)
               noset(81:81)=' '
               ipos=index(noset,' ')
               noset(ipos:ipos)='N'
               do i=1,nset
                  if(set(i).eq.noset) then
                     do j=istartset(i),iendset(i)
                        if(ialset(j).gt.0) then
                           nalset=nalset+1
                           if(nalset.gt.nalset_) then
                              write(*,*) 
     &                           '*ERROR in surfaces: increase nalset_'
                              stop
                           endif
                           ialset(nalset)=ialset(j)
                        else
                           kstart=ialset(nalset-1)+1
                           kend=ialset(nalset)
                           nalset=nalset-1
                           do k=kstart,kend
                              nalset=nalset+1
                              if(nalset.gt.nalset_) then
                                 write(*,*) 
     &                            '*ERROR in surfaces: increase nalset_'
                                 stop
                              endif
                              ialset(nalset)=k
                           enddo
                        endif
                     enddo
                     iendset(nset)=nalset
                     exit
                  endif
               enddo
               if(i.gt.nset) then
                  noset(ipos:ipos)=' '
                  write(*,*) '*ERROR in surfaces: node set ',noset
                  write(*,*) '       does not exist'
                  stop
               endif
            else
               if(ialset(nalset+1).gt.nk) then
                  write(*,*) '*WARNING in surfaces: value ',
     &                 ialset(nalset+1)
                  write(*,*) '         in set ',set(nset),' > nk'
               else
                  nalset=nalset+1
                  iendset(nset)=nalset
               endif
            endif
         enddo
!
      else
!
!        element surface
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) return
            if(nalset+1.gt.nalset_) then
               write(*,*) '*ERROR in surfaces: increase nalset_'
               stop
            endif
!     
c            read(textpart(2)(2:2),'(i1)',iostat=istat) iside
c            if(istat.gt.0) then
c               write(*,*) '*ERROR in surfaces: invalid element side'
c               stop
c            endif
c            if((iside.lt.1).or.(iside.gt.6)) then
c               write(*,*) '*ERROR in surfaces: invalid element side'
c               stop
c            endif
!
            read(textpart(2)(1:20),'(a20)',iostat=istat) label
!     
            if(label(2:4).eq.'NEG') then
               label(2:4)='1  '
            elseif(label(2:4).eq.'POS') then
               label(2:4)='2  '
            endif
            if(label(2:2).eq.'N') then
               label(2:2)='5'
            elseif(label(2:2).eq.'P') then
               label(2:2)='6'
            endif
!
            if((label(1:2).ne.'S1').and.(label(1:2).ne.'S2').and.
     &         (label(1:2).ne.'S3').and.(label(1:2).ne.'S4').and.
     &         (label(1:2).ne.'S5').and.(label(1:2).ne.'S6')) then
               call inputerror(inpc,ipoinpc,iline)
            endif
!            
            read(textpart(1)(1:10),'(i10)',iostat=istat)l
            if(istat.gt.0) then
               elset=textpart(1)(1:80)
               elset(81:81)=' '
               ipos=index(elset,' ')
               elset(ipos:ipos)='E'
               do i=1,nset
                  if(set(i).eq.elset) then
                     do j=istartset(i),iendset(i)
                        l=ialset(j)
                        if(l.gt.0) then
                           kstart=kend
                           kend=l
                           nalset=nalset+1
                           if(nalset.gt.nalset_) then
                              write(*,*) 
     &                           '*ERROR in surfaces: increase nalset_'
                              stop
                           endif
                           if((lakon(l)(1:2).eq.'CP').or.
     &                          (lakon(l)(2:2).eq.'A')) then
                              if(label(1:2).eq.'S1') then
                                 label(1:2)='S3'
                              elseif(label(1:2).eq.'S2') then
                                 label(1:2)='S4'
                              elseif(label(1:2).eq.'S3') then
                                 label(1:2)='S5'
                              elseif(label(1:2).eq.'S4') then
                                 label(1:2)='S6'
                              elseif(label(1:2).eq.'S5') then
                                 label(1:2)='S1'
                              elseif(label(1:2).eq.'S6') then
                                 label(1:2)='S2'
                              endif
                           endif
                           read(label(2:2),'(i1)',iostat=istat) iside
                           ialset(nalset)=iside+10*l
                        else
                           kstart=kstart+1
                           nalset=nalset-1
                           do l=kstart,kend
                              nalset=nalset+1
                              if(nalset.gt.nalset_) then
                                 write(*,*) 
     &                            '*ERROR in surfaces: increase nalset_'
                                 stop
                              endif
                              if((lakon(l)(1:2).eq.'CP').or.
     &                             (lakon(l)(2:2).eq.'A')) then
                                 if(label(1:2).eq.'S1') then
                                    label(1:2)='S3'
                                 elseif(label(1:2).eq.'S2') then
                                    label(1:2)='S4'
                                 elseif(label(1:2).eq.'S3') then
                                    label(1:2)='S5'
                                 elseif(label(1:2).eq.'S4') then
                                    label(1:2)='S6'
                                 elseif(label(1:2).eq.'S5') then
                                    label(1:2)='S1'
                                 elseif(label(1:2).eq.'S6') then
                                    label(1:2)='S2'
                                 endif
                              endif
                              read(label(2:2),'(i1)',iostat=istat) iside
                              ialset(nalset)=iside+10*l
                           enddo
                        endif
                     enddo
                     iendset(nset)=nalset
                     exit
                  endif
               enddo
               if(i.gt.nset) then
                  elset(ipos:ipos)=' '
                  write(*,*) '*ERROR in surfaces: element set ',elset
                  write(*,*) '       does not exist'
                  stop
               endif
            else
               if(l.gt.ne) then
                  write(*,*) '*WARNING in surfaces: value ',
     &                 ialset(nalset+1)
                  write(*,*) '         in set ',set(nset),' > ne'
               else
                  if((lakon(l)(1:2).eq.'CP').or.
     &                 (lakon(l)(2:2).eq.'A')) then
                     if(label(1:2).eq.'S1') then
                        label(1:2)='S3'
                     elseif(label(1:2).eq.'S2') then
                        label(1:2)='S4'
                     elseif(label(1:2).eq.'S3') then
                        label(1:2)='S5'
                     elseif(label(1:2).eq.'S4') then
                        label(1:2)='S6'
                     elseif(label(1:2).eq.'S5') then
                        label(1:2)='S1'
                     elseif(label(1:2).eq.'S6') then
                        label(1:2)='S2'
                     endif
                  endif
                  read(label(2:2),'(i1)',iostat=istat) iside
                  nalset=nalset+1
                  ialset(nalset)=iside+10*l
                  iendset(nset)=nalset
               endif
            endif
         enddo
      endif
!     
      return
      end

