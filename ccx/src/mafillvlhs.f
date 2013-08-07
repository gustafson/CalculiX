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
      subroutine mafillvlhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,
     &  nactdoh,icolv,jqv,irowv,neqv,nzlv,
     &  ikmpc,ilmpc,ikboun,ilboun,nzsv,adbv,aubv,ipvar,var)
!
!     filling the stiffness matrix in spare matrix format (sm)
!
      implicit none
!
      character*8 lakon(*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  icolv(*),jqv(*),ikmpc(*),nzsv,nmethod,
     &  ilmpc(*),ikboun(*),ilboun(*),nactdoh(0:4,*),konl(20),irowv(*),
     &  ipkon(*),ipvar(*)
!
      integer nk,ne,nboun,nmpc,neqv,nzlv,i,j,k,l,m,jj,
     &  ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &  mpc1,mpc2,index1,index2,node1,node2,
     &  indexe,nope,i0
!
      real*8 co(3,*),xboun(*),coefmpc(*),sm(78,78),adbv(*),aubv(*),
     &  var(*)
!
      real*8 value
!
c      write(*,*) 'print nactdoh'
c      do i=1,nk
c         write(*,*) i,(nactdoh(j,i),j=0,4)
c      enddo
!
      i0=0
!
!     determining nzlv
!
      nzlv=0
      do i=neqv,1,-1
         if(icolv(i).gt.0) then
            nzlv=i
            exit
         endif
      enddo
!
      do i=1,neqv
         adbv(i)=0.d0
      enddo
      do i=1,nzsv
         aubv(i)=0.d0
      enddo
!
!     loop over all fluid elements
!
      do i=1,ne
!
        if(ipkon(i).lt.0) cycle
        if(lakon(i)(1:1).ne.'F') cycle
        indexe=ipkon(i)
        if(lakon(i)(4:4).eq.'2') then
           nope=20
        elseif(lakon(i)(4:4).eq.'8') then
           nope=8
        elseif(lakon(i)(4:5).eq.'10') then
           nope=10
        elseif(lakon(i)(4:4).eq.'4') then
           nope=4
        elseif(lakon(i)(4:5).eq.'15') then
           nope=15
        elseif(lakon(i)(4:4).eq.'6') then
           nope=6
        else
           cycle
        endif
!
        do j=1,nope
          konl(j)=kon(indexe+j) 
        enddo
!
        call e_c3d_vlhs(co,nk,konl,lakon(i),sm,i,ipvar,var)
!
        do jj=1,3*nope
!
          j=(jj-1)/3+1
          k=jj-3*(j-1)
!
          node1=kon(indexe+j)
          jdof1=nactdoh(k,node1)
!
          do ll=jj,3*nope
!
            l=(ll-1)/3+1
            m=ll-3*(l-1)
!
            node2=kon(indexe+l)
            jdof2=nactdoh(m,node2)
!
!           check whether one of the DOF belongs to a SPC or MPC
!
            if((jdof1.ne.0).and.(jdof2.ne.0)) then
                  call add_sm_fl(aubv,adbv,jqv,irowv,jdof1,jdof2,
     &                 sm(jj,ll),jj,ll)
            elseif((jdof1.ne.0).or.(jdof2.ne.0)) then
!
!              idof1: genuine DOF
!              idof2: nominal DOF of the SPC/MPC
!
               if(jdof1.eq.0) then
                  idof1=jdof2
                  idof2=(node1-1)*8+k
               else
                  idof1=jdof1
                  idof2=(node2-1)*8+m
               endif
               if(nmpc.gt.0) then
                  call nident(ikmpc,idof2,nmpc,id)
                  if((id.gt.0).and.(ikmpc(id).eq.idof2)) then
!
!                    regular DOF / MPC
!
                     id=ilmpc(id)
                     ist=ipompc(id)
                     index=nodempc(3,ist)
                     if(index.eq.0) cycle
                     do
                        idof2=nactdoh(nodempc(2,index),nodempc(1,index))
                        if(idof2.ne.0) then
                              value=-coefmpc(index)*sm(jj,ll)/
     &                               coefmpc(ist)
                              if(idof1.eq.idof2) value=2.d0*value
                              call add_sm_fl(aubv,adbv,jqv,irowv,
     &                             idof1,idof2,value,i0,i0)
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
            else
               idof1=(node1-1)*8+k
               idof2=(node2-1)*8+m
               mpc1=0
               mpc2=0
               if(nmpc.gt.0) then
                  call nident(ikmpc,idof1,nmpc,id1)
                  if((id1.gt.0).and.(ikmpc(id1).eq.idof1)) mpc1=1
                  call nident(ikmpc,idof2,nmpc,id2)
                  if((id2.gt.0).and.(ikmpc(id2).eq.idof2)) mpc2=1
               endif
               if((mpc1.eq.1).and.(mpc2.eq.1)) then
                  id1=ilmpc(id1)
                  id2=ilmpc(id2)
                  if(id1.eq.id2) then
!
!                    MPC id1 / MPC id1
!
                     ist=ipompc(id1)
                     index1=nodempc(3,ist)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdoh(nodempc(2,index1),
     &                                nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdoh(nodempc(2,index2),
     &                                   nodempc(1,index2))
                           if((idof1.ne.0).and.(idof2.ne.0)) then
                                 value=coefmpc(index1)*coefmpc(index2)*
     &                             sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                                 call add_sm_fl(aubv,adbv,jqv,
     &                             irowv,idof1,idof2,value,i0,i0)
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                   else
!
!                    MPC id1 / MPC id2
!
                     ist1=ipompc(id1)
                     index1=nodempc(3,ist1)
                     if(index1.eq.0) cycle
                     do
                        idof1=nactdoh(nodempc(2,index1),
     &                                nodempc(1,index1))
                        ist2=ipompc(id2)
                        index2=nodempc(3,ist2)
                        if(index2.eq.0) then
                           index1=nodempc(3,index1)
                           if(index1.eq.0) then
                              exit
                           else
                              cycle
                           endif
                        endif
                        do
                           idof2=nactdoh(nodempc(2,index2),
     &                                   nodempc(1,index2))
                           if((idof1.ne.0).and.(idof2.ne.0)) then
                                 value=coefmpc(index1)*coefmpc(index2)*
     &                             sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                                 if(idof1.eq.idof2) value=2.d0*value
                                 call add_sm_fl(aubv,adbv,jqv,
     &                             irowv,idof1,idof2,value,i0,i0)
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
               endif
            endif
          enddo
        enddo
      enddo
!
      return
      end
