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
      subroutine tempload(xforcold,xforc,xforcact,iamforc,nforc,
     &  xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
     &  xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,
     &  amta,namta,nam,ampli,time,reltime,ttime,dtime,ithermal,nmethod,
     &  xbounold,xboun,xbounact,iamboun,nboun,
     &  nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
     &  co,vold,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi)
!
!     calculates the loading at a given time
!
      implicit none
!
      logical gasnode
!
      character*20 sideload(*)
      character*80 amname(*)
!
      integer iamforc(*),iamload(2,*),iamt1(*),nelemload(2,*),
     &  nam,i,istart,iend,id,nforc,nload,nk,namta(3,*),ithermal,
     &  nmethod,iamt1i,iamboun(*),nboun,iamforci,iambouni,
     &  iamloadi1,iamloadi2,ibody(3,*),itg(*),ntg,idof,
     &  nbody,iambodyi,nodeboun(*),ndirboun(*),nodeforc(2,*),
     &  ndirforc(*),istep,iinc,msecpt,node,j,ikboun(*),ilboun(*),
     &  ipresboun,mi(2)
!
      real*8 xforc(*),xforcact(*),xload(2,*),xloadact(2,*),
     &  t1(*),t1act(*),amta(2,*),ampli(*),time,
     &  xforcold(*),xloadold(2,*),t1old(*),reltime,
     &  xbounold(*),xboun(*),xbounact(*),ttime,dtime,reftime,
     &  xbody(7,*),xbodyold(7,*),
     &  xbodyact(7,*),co(3,*),vold(0:mi(2),*),abqtime(2),coords(3)
!
      data msecpt /1/
!
!     if an amplitude is active, the loading is scaled according to
!     the actual time. If no amplitude is active, then the load is
!     - scaled according to the relative time for a static step
!     - applied as a step loading for a dynamic step
!
!     calculating all amplitude values for the current time
!
      do i=1,nam
         if(namta(3,i).lt.0) then
            reftime=ttime+dtime
         else
            reftime=time
         endif
         if(abs(namta(3,i)).ne.i) then
            reftime=reftime-amta(1,namta(1,i))
            istart=namta(1,abs(namta(3,i)))
            iend=namta(2,abs(namta(3,i)))
            if(istart.eq.0) then
               call uamplitude(reftime,amname(namta(3,i)),ampli(i))
               cycle
            endif
         else
            istart=namta(1,i)
            iend=namta(2,i)
            if(istart.eq.0) then
               call uamplitude(reftime,amname(i),ampli(i))
               cycle
            endif
         endif
         call identamta(amta,reftime,istart,iend,id)
         if(id.lt.istart) then
            ampli(i)=amta(2,istart)
         elseif(id.eq.iend) then
            ampli(i)=amta(2,iend)
         else
            ampli(i)=amta(2,id)+(amta(2,id+1)-amta(2,id))
     &           *(reftime-amta(1,id))/(amta(1,id+1)-amta(1,id))
         endif
      enddo
!
!     scaling the boundary conditions
!
      do i=1,nboun
c         if((ithermal.le.1).and.(ndirboun(i).eq.0)) cycle
c         if((ithermal.eq.2).and.(ndirboun(i).gt.0).and.
c     &      (ndirboun(i).le.3)) cycle
         if((xboun(i).lt.1.2357111318d0).and.
     &        (xboun(i).gt.1.2357111316d0)) then
!     
!           user subroutine for boundary conditions
!     
            node=nodeboun(i)
!
!           check whether node is a gasnode
!
            gasnode=.false.
            call nident(itg,node,ntg,id)
            if(id.gt.0) then
               if(itg(id).eq.node) then
                  gasnode=.true.
               endif
            endif
!
            abqtime(1)=time
            abqtime(2)=ttime+dtime
!
!           a gasnode cannot move (displacement DOFs are used
!           for other purposes, e.g. mass flow and pressure)
!
            if(gasnode) then
               do j=1,3
                  coords(j)=co(j,node)
               enddo
            else
               do j=1,3
                  coords(j)=co(j,node)+vold(j,node)
               enddo
            endif
!
            if(ndirboun(i).eq.0) then
               call utemp(xbounact(i),msecpt,istep,iinc,abqtime,node,
     &              coords,vold,mi)
            else
               call uboun(xbounact(i),istep,iinc,abqtime,node,
     &              ndirboun(i),coords,vold,mi)
            endif
            cycle
         endif
!
         if(nam.gt.0) then
            iambouni=iamboun(i)
         else
            iambouni=0
         endif
         if(iambouni.gt.0) then
            xbounact(i)=xboun(i)*ampli(iambouni)
         elseif(nmethod.eq.1) then
            xbounact(i)=xbounold(i)+
     &         (xboun(i)-xbounold(i))*reltime
         else
            xbounact(i)=xboun(i)
         endif
      enddo
!
!     scaling the loading
!
      do i=1,nforc
         if(ndirforc(i).eq.0) then
            if((xforc(i).lt.1.2357111318d0).and.
     &         (xforc(i).gt.1.2357111316d0)) then
!
!              user subroutine for the concentrated heat flux
!
               node=nodeforc(1,i)
!
!              check whether node is a gasnode
!
               gasnode=.false.
               call nident(itg,node,ntg,id)
               if(id.gt.0) then
                  if(itg(id).eq.node) then
                     gasnode=.true.
                  endif
               endif
!
               abqtime(1)=time
               abqtime(2)=ttime+dtime
!
!              a gasnode cannot move (displacement DOFs are used
!              for other purposes, e.g. mass flow and pressure)
!
               if(gasnode) then
                  do j=1,3
                     coords(j)=co(j,node)
                  enddo
               else
                  do j=1,3
                     coords(j)=co(j,node)+vold(j,node)
                  enddo
               endif
!
               call cflux(xforcact(i),msecpt,istep,iinc,abqtime,node,
     &              coords,vold,mi)
               cycle
            endif
         endif
         if(nam.gt.0) then
            iamforci=iamforc(i)
         else
            iamforci=0
         endif
         if(iamforci.gt.0) then
            xforcact(i)=xforc(i)*ampli(iamforci)
         elseif(nmethod.eq.1) then
            xforcact(i)=xforcold(i)+
     &         (xforc(i)-xforcold(i))*reltime
         else
            xforcact(i)=xforc(i)
         endif
      enddo
!
      do i=1,nload
         ipresboun=0
!
!        check for pressure boundary conditions
!
         if(sideload(i)(3:4).eq.'NP') then
            node=nelemload(2,i)
            idof=8*(node-1)+2
            call nident(ikboun,idof,nboun,id)
            if(id.gt.0) then
               if(ikboun(id).eq.idof) then
                  ipresboun=1
                  xloadact(1,i)=xbounact(ilboun(id))
               endif
            endif
         endif
!
         if(ipresboun.eq.0) then
            if(nam.gt.0) then
               iamloadi1=iamload(1,i)
               iamloadi2=iamload(2,i)
            else
               iamloadi1=0
               iamloadi2=0
            endif
            if(iamloadi1.gt.0) then
               xloadact(1,i)=xload(1,i)*ampli(iamloadi1)
            elseif(nmethod.eq.1) then
               xloadact(1,i)=xloadold(1,i)+
     &              (xload(1,i)-xloadold(1,i))*reltime
            else
               xloadact(1,i)=xload(1,i)
            endif
            if(iamloadi2.gt.0) then
               xloadact(2,i)=xload(2,i)*ampli(iamloadi2)
            elseif(nmethod.eq.1) then
               xloadact(2,i)=xload(2,i)
            else
               xloadact(2,i)=xload(2,i)
            endif
         endif
      enddo
!
      do i=1,nbody
         if(nam.gt.0) then
            iambodyi=ibody(2,i)
         else
            iambodyi=0
         endif
         if(iambodyi.gt.0) then
            xbodyact(1,i)=xbody(1,i)*ampli(iambodyi)
         elseif(nmethod.eq.1) then
            xbodyact(1,i)=xbodyold(1,i)+
     &           (xbody(1,i)-xbodyold(1,i))*reltime
         else
            xbodyact(1,i)=xbody(1,i)
         endif
      enddo
!
!     scaling the temperatures
!
      if(ithermal.eq.1) then
         do i=1,nk
            if((t1(i).lt.1.2357111318d0).and.
     &           (t1(i).gt.1.2357111316d0)) then
!
               abqtime(1)=time
               abqtime(2)=ttime+dtime
!
               do j=1,3
                  coords(j)=co(j,i)+vold(j,i)
               enddo
               call utemp(t1act(i),msecpt,istep,iinc,abqtime,i,
     &              coords,vold,mi)
               cycle
            endif
            if(nam.gt.0) then
               iamt1i=iamt1(i)
            else
               iamt1i=0
            endif
            if(iamt1i.gt.0) then
               t1act(i)=t1(i)*ampli(iamt1i)
            elseif(nmethod.eq.1) then
               t1act(i)=t1old(i)+(t1(i)-t1old(i))*reltime
            else
               t1act(i)=t1(i)
            endif
         enddo
      endif
c      write(*,*) 'nboun'
c      do i=1,nboun
c         write(*,'(i7,1x,e11.4,1x,e11.4)') i,xbounact(i),xboun(i)
c      enddo
c      write(*,*) 'nforc'
c      do i=1,nforc
c         write(*,'(i7,1x,e11.4,1x,e11.4)') i,xforcact(i),xforc(i)
c      enddo
c      write(*,*) 'nload'
c      do i=1,nload
c         write(*,'(i7,1x,e11.4,1x,e11.4)') i,xloadact(1,i),xload(1,i)
c      enddo
!
      return
      end
