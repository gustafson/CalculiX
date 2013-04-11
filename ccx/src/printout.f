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
      subroutine printout(set,nset,istartset,iendset,ialset,nprint,
     &  prlab,prset,v,t1,fn,ipkon,lakon,stx,eme,xstate,ener,
     &  mint_,nstate_,ithermal,co,kon,qfx,ttime,trab,inotr,ntrans,
     &  orab,ielorien,norien,nk,ne,inum,filab)
!
!     stores results in the .dat file
!
      implicit none
!
      character*1 cflag
      character*6 prlab(*),filab(*)
      character*8 lakon(*)
      character*80 noset,elset
      character*81 set(*),prset(*)
!
      integer nset,istartset(*),iendset(*),ialset(*),nprint,ipkon(*),
     &  mint_,nstate_,k,ii,jj,iset,l,lb,limit,node,ipos,ithermal,
     &  nelem,kon(*),inotr(2,*),ntrans,ielorien(*),norien,nk,ne,
     &  inum(*),nfield
!
      real*8 v(0:3,*),t1(*),fn(0:3,*),stx(6,mint_,*),
     &  eme(6,mint_,*),xstate(nstate_,mint_,*),ener(mint_,*),energytot,
     &  volumetot,co(3,*),qfx(3,mint_,*),rftot(0:3),ttime,
     &  trab(7,*),orab(7,*)
!
!     interpolation in the original nodes of 1d and 2d elements
!
      do ii=1,nprint
         if((prlab(ii)(1:4).eq.'U   ').or.
     &      ((prlab(ii)(1:4).eq.'NT  ').and.(ithermal.gt.1))) then
            if(filab(1)(5:5).ne.' ') then
               nfield=4
               cflag=filab(1)(5:5)
               call map3dto1d2d(v,ipkon,inum,kon,lakon,nfield,nk,
     &              ne,cflag,co)
            endif
            exit
          endif
      enddo
      do ii=1,nprint
         if((prlab(ii)(1:4).eq.'NT  ').and.(ithermal.le.1)) then
            if(filab(2)(5:5).ne.' ') then
               nfield=1
               cflag=filab(2)(5:5)
               call map3dto1d2d(t1,ipkon,inum,kon,lakon,nfield,nk,
     &              ne,cflag,co)
            endif
            exit
          endif
      enddo
c      do ii=1,nprint
c         if(prlab(ii)(1:2).eq.'RF') then
c            if(filab(1)(5:5).ne.' ') then
c               nfield=4
c               call map3dto1d2d(fn,ipkon,inum,kon,lakon,nfield,nk,
c     &              ne,cflag,co)
c            endif
c            exit
c          endif
c      enddo
!
      do ii=1,nprint
!
!        nodal values
!
         if((prlab(ii)(1:4).eq.'U   ').or.(prlab(ii)(1:4).eq.'NT  ').or.
     &      (prlab(ii)(1:4).eq.'RF  ').or.(prlab(ii)(1:4).eq.'RFL ')) 
     &      then
!
            ipos=index(prset(ii),' ')
            noset='                    '
            noset(1:ipos-1)=prset(ii)(1:ipos-1)
!
!           printing the header
!
            if(prlab(ii)(1:4).eq.'U   ') then
               write(5,*)
c               write(5,*) 'displacements (vx,vy,vz) for set ',
c     &            noset(1:ipos-2),' and time ',ttime
               write(5,100) noset(1:ipos-2),ttime
 100           format(' displacements (vx,vy,vz) for set ',A,
     &             ' and time ',e11.4)
               write(5,*)
            elseif(prlab(ii)(1:4).eq.'NT  ') then
               write(5,*)
c               write(5,*) 'temperatures for set ',noset(1:ipos-2),
c     &               ' and time ',ttime
               write(5,101) noset(1:ipos-2),ttime
 101           format(' temperatures for set ',A,' and time ',e11.4)
               write(5,*)
            elseif((prlab(ii)(1:5).eq.'RF   ').or.
     &             (prlab(ii)(1:5).eq.'RF  T')) then
               write(5,*)
c               write(5,*) 'forces (fx,fy,fz) for set ',noset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,102) noset(1:ipos-2),ttime
 102           format(' forces (fx,fy,fz) for set ',A,
     &                ' and time ',e11.4)
               write(5,*)
            elseif((prlab(ii)(1:5).eq.'RFL ').or.
     &             (prlab(ii)(1:5).eq.'RFL T')) then
               write(5,*)
c               write(5,*) 'heat generation for set ',noset(1:ipos-2),
c     &           ' and time ',ttime
               write(5,103) noset(1:ipos-2),ttime
 103           format(' heat generation for set ',A,' and time ',e11.4)
               write(5,*)
            endif
!
!           printing the data
!
            do iset=1,nset
               if(set(iset).eq.prset(ii)) exit
            enddo
            do jj=0,3
               rftot(jj)=0.d0
            enddo
            do jj=istartset(iset),iendset(iset)
               if(ialset(jj).lt.0) cycle
               if(jj.eq.iendset(iset)) then
                  node=ialset(jj)
                  call printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &              rftot,trab,inotr,ntrans,co)
               elseif(ialset(jj+1).gt.0) then
                  node=ialset(jj)
                  call printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &              rftot,trab,inotr,ntrans,co)
               else
                  do node=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &                 -ialset(jj+1)
                  call printoutnode(prlab,v,t1,fn,ithermal,ii,node,
     &              rftot,trab,inotr,ntrans,co)
                  enddo
               endif
            enddo
!
!           writing total values to file
!
            if((prlab(ii)(1:5).eq.'RF  O').or.
     &           (prlab(ii)(1:5).eq.'RF  T')) then
               write(5,*)
c               write(5,*) 
c     &              'total force (fx,fy,fz) for set ',noset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,104) noset(1:ipos-2),ttime
 104           format(' total force (fx,fy,fz) for set ',A,
     &                 ' and time ',e11.4)
               write(5,*)
               write(5,'(6x,1p,3(1x,e11.4))') rftot(1),rftot(2),rftot(3)
            elseif((prlab(ii)(1:5).eq.'RFL O').or.
     &              (prlab(ii)(1:5).eq.'RFL T')) then
               write(5,*)
c               write(5,*) 
c     &              'total heat generation for set ',noset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,105)noset(1:ipos-2),ttime
 105           format(' total heat generation for set ',A,
     &                ' and time ',e11.4)
               write(5,*)
               write(5,'(6x,1p,1x,e11.4)') rftot(0)
            endif
!
!        integration point values
!
         elseif((prlab(ii)(1:4).eq.'S   ').or.
     &          (prlab(ii)(1:4).eq.'E   ').or.
     &          (prlab(ii)(1:4).eq.'PE  ').or.
     &          (prlab(ii)(1:4).eq.'ENER').or.
     &          (prlab(ii)(1:4).eq.'SDV ').or.
     &          (prlab(ii)(1:4).eq.'HFL ')) then
!
            ipos=index(prset(ii),' ')
            elset='                    '
            elset(1:ipos-1)=prset(ii)(1:ipos-1)
!
!           for SDV more than six columns result
!
            if(prlab(ii)(1:4).eq.'SDV ') then
               limit=(nstate_+5)/6
            else
               limit=1
            endif
!
            do l=1,limit
!
!              printing the header
!
               if(prlab(ii)(1:4).eq.'S   ') then
                  write(5,*)
c                  write(5,*) 
c     &                 'stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sxz,s
c     &yz) for set ',elset(1:ipos-2),' and time ',ttime
                  write(5,106) elset(1:ipos-2),ttime
 106              format(' stresses (elem, integ.pnt.,sxx,syy,szz,sxy,sx
     &z,syz) for set ',A,' and time ',e11.4)
                  write(5,*)
               elseif(prlab(ii)(1:4).eq.'E   ') then
                  write(5,*)
c                  write(5,*) 
c     &                 'strains (elem, integ.pnt.,exx,eyy,ezz,exy,exz,ey
c     &z) forset ',elset(1:ipos-2),' and time ',ttime
                  write(5,107) elset(1:ipos-2),ttime
 107              format(' strains (elem, integ.pnt.,exx,eyy,ezz,exy,exz
     &,eyz) forset ',A,' and time ',e11.4)
                  write(5,*)
               elseif(prlab(ii)(1:4).eq.'PE  ') then
                  write(5,*)
c                  write(5,*) 
c     &                 'equivalent plastic strain (elem, integ.pnt.,pe) 
c     &for set ',elset(1:ipos-2),' and time ',ttime
                  write(5,108) elset(1:ipos-2),ttime
 108              format(' equivalent plastic strain (elem, integ.pnt.,p 
     &e)for set ',A,' and time ',e11.4)
                  write(5,*)
               elseif(prlab(ii)(1:4).eq.'ENER') then
                  write(5,*)
c                  write(5,*) 
c     &                 'energy density (elem, integ.pnt.,energy) for set
c     & ',elset(1:ipos-2),' and time ',ttime
                  write(5,109) elset(1:ipos-2),ttime
 109              format(' energy density (elem, integ.pnt.,energy) for 
     &set ',A,' and time ',e11.4)
                  write(5,*)
               elseif(prlab(ii)(1:4).eq.'SDV ') then
                  lb=(l-1)*6
                  write(5,*)
                  if(l.eq.(nstate_+5)/6) then
                     write(5,300) (lb+k,k=1,nstate_-lb)
                  else
                     write(5,300) (lb+k,k=1,6)
                  endif
 300              format(' internal state variables (elem, integ.pnt.,',
     &                    6(i2),')')
c                  write(5,*) 'for set ',elset,' and time ',ttime
                  write(5,110) elset,ttime
 110              format(' for set ',A,' and time ',e11.4)
                  write(5,*)
               elseif(prlab(ii)(1:4).eq.'HFL ') then
                  write(5,*)
c                  write(5,*) 
c     &                 'heat flux (elem, integ.pnt.,qx,qy,qz) for set ',
c     &                    elset(1:ipos-2),' and time ',ttime
                  write(5,111) elset(1:ipos-2),ttime
 111              format(' heat flux (elem, integ.pnt.,qx,qy,qz) for set 
     & ',A,' and time ',e11.4)
                  write(5,*)
               endif
!
!           printing the data
!
               do iset=1,nset
                  if(set(iset).eq.prset(ii)) exit
               enddo
               do jj=istartset(iset),iendset(iset)
                  if(ialset(jj).lt.0) cycle
                  if(jj.eq.iendset(iset)) then
                     nelem=ialset(jj)
                     call printoutint(prlab,ipkon,lakon,stx,eme,xstate,
     &                    ener,mint_,nstate_,l,lb,ii,nelem,qfx,
     &                    orab,ielorien,norien,co,kon)
                  elseif(ialset(jj+1).gt.0) then
                     nelem=ialset(jj)
                     call printoutint(prlab,ipkon,lakon,stx,eme,xstate,
     &                    ener,mint_,nstate_,l,lb,ii,nelem,qfx,orab,
     &                    ielorien,norien,co,kon)
                  else
                     do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &                    -ialset(jj+1)
                        call printoutint(prlab,ipkon,lakon,stx,eme,
     &                       xstate,ener,mint_,nstate_,l,lb,ii,nelem,
     &                       qfx,orab,ielorien,norien,co,kon)
                     enddo
                  endif
               enddo
!
            enddo
!
!        whole element values
!
         elseif((prlab(ii)(1:4).eq.'ELSE').or.
     &          (prlab(ii)(1:4).eq.'EVOL')) then
!
            ipos=index(prset(ii),' ')
            elset='                    '
            elset(1:ipos-1)=prset(ii)(1:ipos-1)
!
!           printing the header
!
            if((prlab(ii)(1:5).eq.'ELSE ').or.
     &         (prlab(ii)(1:5).eq.'ELSET')) then
               write(5,*)
c               write(5,*) 
c     &              'energy (element, energy) for set ',elset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,112) elset(1:ipos-2),ttime
 112           format(' energy (element, energy) for set ',A,
     &                ' and time ',e11.4)
               write(5,*)
            elseif((prlab(ii)(1:5).eq.'EVOL ').or.
     &             (prlab(ii)(1:5).eq.'EVOLT')) then
               write(5,*)
c               write(5,*) 
c     &              'volume (element, volume) for set ',elset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,113) elset(1:ipos-2),ttime
 113           format(' volume (element, volume) for set ',A,
     &                ' and time ',e11.4)
               write(5,*)
            endif
!
!           printing the data
!
            do iset=1,nset
               if(set(iset).eq.prset(ii)) exit
            enddo
            volumetot=0.d0
            energytot=0.d0
            do jj=istartset(iset),iendset(iset)
               if(ialset(jj).lt.0) cycle
               if(jj.eq.iendset(iset)) then
                  nelem=ialset(jj)
                  call printoutelem(prlab,ipkon,lakon,kon,co,
     &                 ener,mint_,ii,nelem,energytot,volumetot)
               elseif(ialset(jj+1).gt.0) then
                  nelem=ialset(jj)
                  call printoutelem(prlab,ipkon,lakon,kon,co,
     &                 ener,mint_,ii,nelem,energytot,volumetot)
               else
                  do nelem=ialset(jj-1)-ialset(jj+1),ialset(jj),
     &                 -ialset(jj+1)
                     call printoutelem(prlab,ipkon,lakon,kon,co,
     &                    ener,mint_,ii,nelem,energytot,volumetot)
                  enddo
               endif
            enddo
!
!     writing total values to file
!
            if((prlab(ii)(1:5).eq.'ELSEO').or.
     &           (prlab(ii)(1:5).eq.'ELSET')) then
               write(5,*)
c               write(5,*) 
c     &              'total energy for set ',elset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,*) elset(1:ipos-2),ttime
 114           format(' total energy for set ',A,' and time ',e11.4)
               write(5,*)
               write(5,'(6x,1p,1x,e11.4)') energytot
            elseif((prlab(ii)(1:5).eq.'EVOLO').or.
     &              (prlab(ii)(1:5).eq.'EVOLT')) then
               write(5,*)
c               write(5,*) 
c     &              'total volume for set ',elset(1:ipos-2),
c     &              ' and time ',ttime
               write(5,*) elset(1:ipos-2),ttime
 116           format(' total volume for set ',A,' and time ',e11.4)
               write(5,*)
               write(5,'(6x,1p,1x,e11.4)') volumetot
            endif
!
         endif
      enddo
!                     
      return
      end
