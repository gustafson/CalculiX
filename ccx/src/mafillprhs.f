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
      subroutine mafillprhs(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,
     &  xboun,nboun,ipompc,nodempc,coefmpc,nmpc,nelemface,sideface,
     &  nface,b,nactdoh,icolp,jqp,irowp,neqp,nzlp,nmethod,ikmpc,ilmpc,
     &  ikboun,ilboun,rhcon,nrhcon,ielmat,ntmat_,vold,voldcon,nzsp,
     &  dtl,matname,mi,ncmat_,shcon,nshcon,v,theta1,
     &  iexplicit,physcon,nea,neb,dtimef,ipvar,var,ipvarf,varf)
!
!     filling the rhs b of the pressure equations (step 2)
!
      implicit none
!
      character*1 sideface(*)
      character*8 lakon(*)
      character*80 matname(*)
!
      integer kon(*),nodeboun(*),ndirboun(*),ipompc(*),nodempc(3,*),
     &  nelemface(*),icolp(*),jqp(*),ikmpc(*),ilmpc(*),ikboun(*),
     &  ilboun(*),nactdoh(0:4,*),konl(20),irowp(*),nrhcon(*),
     &  mi(*),ielmat(mi(3),*),
     &  ipkon(*),nshcon(*),iexplicit,nea,neb,ipvar(*),ipvarf(*)
!
      integer nk,ne,nboun,nmpc,nface,neqp,nzlp,nmethod,nzsp,i,j,k,l,jj,
     &  ll,id,id1,id2,ist,ist1,ist2,index,jdof1,jdof2,idof1,idof2,
     &  mpc1,mpc2,index1,index2,node1,node2,kflag,ntmat_,indexe,nope,
     &  i0,ncmat_,idof3
!
      real*8 co(3,*),xboun(*),coefmpc(*),b(*),v(0:mi(2),*),
     &  vold(0:mi(2),*),
     &  voldcon(0:4,*),ff(60),sm(60,60),rhcon(0:1,ntmat_,*),
     &  shcon(0:3,ntmat_,*),theta1,physcon(*),var(*),varf(*)
!
      real*8 value,dtl(*),dtimef
!
      kflag=2
      i0=0
!
!     determining nzlp
!
c      if(iexplicit) then
c         nzlp=0
c         do i=neqp,1,-1
c            if(icolp(i).gt.0) then
c               nzlp=i
c               exit
c            endif
c         enddo
c      endif
!
      do i=1,neqp
         b(i)=0.d0
      enddo
!
      do i=nea,neb
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
c        do j=1,nope
c          konl(j)=kon(indexe+j) 
c        enddo
!
        call e_c3d_prhs(co,nk,kon(indexe+1),lakon(i),sm,ff,i,nmethod,
     &       rhcon,
     &       nrhcon,ielmat,ntmat_,v,vold,voldcon,nelemface,sideface,
     &       nface,dtimef,matname,mi(1),shcon,nshcon,theta1,physcon,
     &       iexplicit,ipvar,var,ipvarf,varf)
!
        do jj=1,nope
!
          j=jj
          k=jj-3*(j-1)
!
          node1=kon(indexe+j)
c          ff(jj)=ff(jj)*dtl(node1)/dtimef
          jdof1=nactdoh(4,node1)
!
          do ll=jj,nope
!
            l=ll
!
            node2=kon(indexe+l)
            jdof2=nactdoh(4,node2)
!
!           check whether one of the DOF belongs to a SPC or MPC
!
            if((jdof1.ne.0).and.(jdof2.ne.0)) then
c                  call add_sm_fl(aubp,adbp,jqp,irowp,jdof1,jdof2,
c     &                sm(jj,ll),jj,ll)
            elseif((jdof1.ne.0).or.(jdof2.ne.0)) then
!
!              idof1: genuine DOF
!              idof2: nominal DOF of the SPC/MPC
!
               if(jdof1.eq.0) then
                  idof1=jdof2
                  idof2=(node1-1)*8+4
               else
                  idof1=jdof1
                  idof2=(node2-1)*8+4
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
                        idof2=nactdoh(4,nodempc(1,index))
                        value=-coefmpc(index)*sm(jj,ll)/coefmpc(ist)
                        if(idof1.eq.idof2) value=2.d0*value
                        if(idof2.ne.0) then
c                              call add_sm_fl(aubp,adbp,jqp,irowp,
cd     &                             idof1,idof2,value,i0,i0)
c                        else
c                           if(iexplicit.eq.0) then
c                              idof2=8*(nodempc(1,index)-1)+
c     &                            nodempc(2,index)
c                              call nident(ikboun,idof2,nboun,id)
c                              b(idof1)=b(idof1)
c     &                             -xboun(ilboun(id))*value
c                           endif
cd
                        endif
                        index=nodempc(3,index)
                        if(index.eq.0) exit
                     enddo
                     cycle
                  endif
               endif
c!
c!              regular DOF / SPC
c!
cd
c               if(iexplicit.eq.0) then
cc                  idof2=idof2+4
c                  call nident(ikboun,idof2,nboun,id)
c                  b(idof1)=b(idof1)-xboun(ilboun(id))*sm(jj,ll)
c               endif
cd
            else
               idof1=(node1-1)*8+4
               idof2=(node2-1)*8+4
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
                        idof1=nactdoh(4,nodempc(1,index1))
                        index2=index1
                        do
                           idof2=nactdoh(4,nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*
     &                          sm(jj,ll)/coefmpc(ist)/coefmpc(ist)
                           if((idof1.ne.0).and.(idof2.ne.0)) then
c                                 call add_sm_fl(aubp,adbp,jqp,
c     &                             irowp,idof1,idof2,value,i0,i0)
cd
c                           elseif((iexplicit.eq.0).and.
c     &                       ((idof1.ne.0).or.(idof2.ne.0))) then
c                              if(idof2.ne.0) then
c                                 idof3=idof2
c                                 idof2=8*(nodempc(1,index1)-1)+
c     &                                 nodempc(2,index1)
c                              else
c                                 idof3=idof1
c                                 idof2=8*(nodempc(1,index2)-1)+
c     &                                 nodempc(2,index2)
c                              endif
c                              call nident(ikboun,idof2,nboun,id)
c                              b(idof3)=b(idof3)
c     &                             -value*xboun(ilboun(id))
cd
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
                        idof1=nactdoh(4,nodempc(1,index1))
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
                           idof2=nactdoh(4,nodempc(1,index2))
                           value=coefmpc(index1)*coefmpc(index2)*
     &                          sm(jj,ll)/coefmpc(ist1)/coefmpc(ist2)
                           if(idof1.eq.idof2) value=2.d0*value
                           if((idof1.ne.0).and.(idof2.ne.0)) then
c                                 call add_sm_fl(aubp,adbp,jqp,
c     &                             irowp,idof1,idof2,value,i0,i0)
cd
c                           elseif((iexplicit.eq.0).and.
c     &                       ((idof1.ne.0).or.(idof2.ne.0))) then
c                              if(idof2.ne.0) then
c                                 idof3=idof2
c                                 idof2=8*(nodempc(1,index1)-1)+
c     &                                 nodempc(2,index1)
c                              else
c                                 idof3=idof1
c                                 idof2=8*(nodempc(1,index2)-1)+
c     &                                 nodempc(2,index2)
c                              endif
c                              call nident(ikboun,idof2,nboun,id)
c                              b(idof3)=b(idof3)
c     &                               -value*xboun(ilboun(id))
cd
                           endif
!
                           index2=nodempc(3,index2)
                           if(index2.eq.0) exit
                        enddo
                        index1=nodempc(3,index1)
                        if(index1.eq.0) exit
                     enddo
                  endif
cd
c               elseif(((mpc1.eq.1).or.(mpc2.eq.1)).and.(iexplicit.eq.0))
c     &           then
c                  if(mpc1.eq.1) then
c!
c!                    MPC id1 / SPC
c!
cc                     idof2=idof2+4
c                     call nident(ikboun,idof2,nboun,id2)
c                     idof2=ilboun(id2)
c                     ist1=ipompc(id1)
c                     index1=nodempc(3,ist1)
c                     if(index1.eq.0) cycle
c                     do
c                        idof1=nactdoh(nodempc(2,index1),
c     &                                nodempc(1,index1))
c                        if(idof1.ne.0) then
c                           b(idof1)=b(idof1)+xboun(idof2)*
c     &                          coefmpc(index1)*sm(jj,ll)/coefmpc(ist1)
c                        endif
c                        index1=nodempc(3,index1)
c                        if(index1.eq.0) exit
c                     enddo
c                  elseif(mpc2.eq.1) then
c!
c!                    MPC id2 / SPC
c!
cc                     idof1=idof1+4
c                     call nident(ikboun,idof1,nboun,id1)
c                     idof1=ilboun(id1)
c                     ist2=ipompc(id2)
c                     index2=nodempc(3,ist2)
c                     if(index2.eq.0) cycle
c                     do
c                        idof2=nactdoh(nodempc(2,index2),
c     &                                nodempc(1,index2))
c                        if(idof2.ne.0) then
c                           b(idof2)=b(idof2)+xboun(idof1)*
c     &                          coefmpc(index2)*sm(jj,ll)/coefmpc(ist2)
c                        endif
c                        index2=nodempc(3,index2)
c                        if(index2.eq.0) exit
c                     enddo
c                  endif
cd
               endif
            endif
          enddo
!
!            inclusion of ff
!
                if(jdof1.eq.0) then
                   if(nmpc.ne.0) then
                      idof1=(node1-1)*8+4
                      call nident(ikmpc,idof1,nmpc,id)
                      if((id.gt.0).and.(ikmpc(id).eq.idof1)) then
                         id=ilmpc(id)
                         ist=ipompc(id)
                         index=nodempc(3,ist)
                         if(index.eq.0) cycle
                         do
                            jdof1=nactdoh(4,nodempc(1,index))
                            if(jdof1.ne.0) then
                               b(jdof1)=b(jdof1)
     &                              -coefmpc(index)*ff(jj)
     &                              /coefmpc(ist)
                            endif
                            index=nodempc(3,index)
                            if(index.eq.0) exit
                         enddo
                      endif
                   endif
                   cycle
                endif
                b(jdof1)=b(jdof1)+ff(jj)
!
        enddo
      enddo
!
!     nonlocal time stepping for compressible steady state calculations
!     
c      if((iexplicit.eq.1).and.(nmethod.eq.1)) then
c         do i=1,nk
c            if(nactdoh(4,i).gt.0) then
c               b(nactdoh(4,i))=b(nactdoh(4,i))*dtl(i)/dtimef
c            endif
c         enddo
c      endif
!
      return
      end
