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
      subroutine nonlinmpc(co,vold,ipompc,nodempc,coefmpc,labmpc,
     &  nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,maxlenmpc,ikmpc,
     &  ilmpc,icascade,kon,ipkon,lakon,ne,reltime,newstep,xboun,fmpc,
     &  iit,idiscon,ncont,trab,ntrans,ithermal)
!
!     updates the coefficients in nonlinear MPC's
!
      implicit none
!
      logical isochoric
!
      character*8 lakon(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),irefnode,irotnode,idir,
     &  nmpc,index,ii,inode,node,id,ikboun(*),ilboun(*),nboun,
     &  i,j,k,idof,na,nb,nc,np,i1,i2,i3,iaux(*),maxlenmpc,n,
     &  l,m,lmax,mmax,ikmpc(*),ilmpc(*),icascade,neigh(7,8),
     &  mpc,kon(*),ipkon(*),indexe,ne,idofrem,idofins,nmpc0,nmpc01,
     &  newstep,iit,idiscon,ncont,iexpnode,indexexp,nmpcdif,ntrans,
     &  nodei,noded,lathyp(3,6),inum,ndir,number,ithermal
!
      real*8 co(3,*),coefmpc(*),vold(0:4,*),c(3,3),dc(3,3,3),ww,
     &  e(3,3,3),d(3,3),w(3),f(3,3),c1,c2,c3,c4,c5,c6,xbounact(*),
     &  xboun(*),fmpc(*),expan
      real*8 dd,a11,a12,a13,a21,a22,a23,a31,a32,a33,
     &       b11,b12,b13,b21,b22,b23,b31,b32,b33,aux(*),const,
     &       ddmax,a(3,3),b(3,3),xj,xi,et,ze,xlag(3,20),xeul(3,20),
     &       coloc(3,8),reltime,csab(7),trab(7,*),pd(3),pi(3),
     &       ad(3,3),ai(3,3)
!
      data d /1.,0.,0.,0.,1.,0.,0.,0.,1./
      data e /0.,0.,0.,0.,0.,-1.,0.,1.,0.,
     &        0.,0.,1.,0.,0.,0.,-1.,0.,0.,
     &        0.,-1.,0.,1.,0.,0.,0.,0.,0./
      data neigh /1,9,2,12,4,17,5,2,9,1,10,3,18,6,
     &            3,11,4,10,2,19,7,4,11,3,12,1,20,8,
     &            5,13,6,16,8,17,1,6,13,5,14,7,18,2,
     &            7,15,8,14,6,19,3,8,15,7,16,5,20,4/
      data coloc /-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,
     &            -1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1./
!     
!     latin hypercube positions in a 3 x 3 matrix
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
!
      irotnode=0
      if((icascade.eq.1).and.(newstep.ne.1).and.(ncont.eq.0)) icascade=0
      isochoric=.false.
!
      ii=0
      loop: do
         ii=ii+1
         if(ii.gt.nmpc) exit
         if(labmpc(ii)(1:5).eq.'RIGID') then
!
            index=ipompc(ii)
            inode=nodempc(1,index)
            idir=nodempc(2,index)
            coefmpc(index)=1.d0
!
            index=nodempc(3,index)
            irefnode=nodempc(1,index)
            coefmpc(index)=-1.d0
!
            index=nodempc(3,index)
            node=nodempc(1,index)
!
!           check whether the rotational node is the same as in
!           the last rigid body MPC
!
            if(node.ne.irotnode) then
               irotnode=node
               w(1)=vold(1,node)
               w(2)=vold(2,node)
               w(3)=vold(3,node)
c               write(*,*) 'w   ',w(1),w(2),w(3)
               ww=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
!
               c1=dcos(ww)
               if(ww.gt.1.d-10) then
                  c2=dsin(ww)/ww
               else
                  c2=1.d0
               endif
               if(ww.gt.1.d-5) then
                  c3=(1.d0-c1)/ww**2
               else
                  c3=0.5d0
               endif
!
!              rotation matrix c
!
               do i=1,3
                  do j=1,3
                     c(i,j)=c1*d(i,j)+
     &                  c2*(e(i,1,j)*w(1)+e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                  c3*w(i)*w(j)
                  enddo
               enddo
!
               c4=-c2
               if(ww.gt.0.00464159) then
                  c5=(ww*dcos(ww)-dsin(ww))/ww**3
               else
                  c5=-1.d0/3.d0
               endif
               if(ww.gt.0.0031623) then
                  c6=(ww*dsin(ww)-2.d0+2.d0*dcos(ww))/ww**4
               else
                  c6=-1.d0/12.d0
               endif
!
!              derivative of the rotation matrix c with respect to
!              the rotation vector w
!
               do i=1,3
                  do j=1,3
                     do k=1,3
                        dc(i,j,k)=c4*w(k)*d(i,j)+
     &                            c5*w(k)*(e(i,1,j)*w(1)+
     &                                     e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                            c2*e(i,k,j)+
     &                            c6*w(k)*w(i)*w(j)+
     &                            c3*(d(i,k)*w(j)+d(j,k)*w(i))
                     enddo
                  enddo
               enddo
!
!              dummy variable
!
               do i=1,3
                  do j=1,3
c                     f(i,j)=c(i,j)-d(i,j)-dc(i,j,1)*w(1)-dc(i,j,2)*w(2)-
c     &                                    dc(i,j,3)*w(3)
                     f(i,j)=c(i,j)-d(i,j)
                  enddo
               enddo
            endif
!
!           determining the coefficients of the rotational degrees
!           of freedom
!
            coefmpc(index)=dc(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,1)*(co(3,irefnode)-co(3,inode))
!
            index=nodempc(3,index)
            coefmpc(index)=dc(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,2)*(co(3,irefnode)-co(3,inode))
!
            index=nodempc(3,index)
            coefmpc(index)=dc(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,3)*(co(3,irefnode)-co(3,inode))
!
!           determining the nonhomogeneous part
!
            index=nodempc(3,index)
            coefmpc(index)=1.d0
!
!           old value of the nonhomogeneous term must be zero
!
            vold(nodempc(2,index),nodempc(1,index))=0.d0
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=f(idir,1)*(co(1,irefnode)-co(1,inode))+
     &           f(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           f(idir,3)*(co(3,irefnode)-co(3,inode))-
     &           vold(idir,irefnode)+vold(idir,inode)
!
         elseif(labmpc(ii)(1:4).eq.'KNOT') then
!
!           dependent node
!
            index=ipompc(ii)
            inode=nodempc(1,index)
            idir=nodempc(2,index)
            coefmpc(index)=1.d0
!
!           translation node
!
            index=nodempc(3,index)
            irefnode=nodempc(1,index)
            coefmpc(index)=-1.d0
!
!           expansion node
!
            index=nodempc(3,index)
            iexpnode=nodempc(1,index)
            expan=1.d0+vold(1,iexpnode)
            indexexp=index
!
!           rotation node
!
            index=nodempc(3,index)
            node=nodempc(1,index)
!
!           check whether the rotational node is the same as in
!           the last rigid body MPC
!
            if(node.ne.irotnode) then
               irotnode=node
               w(1)=vold(1,node)
               w(2)=vold(2,node)
               w(3)=vold(3,node)
               ww=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
!
               c1=dcos(ww)
               if(ww.gt.1.d-10) then
                  c2=dsin(ww)/ww
               else
                  c2=1.d0
               endif
               if(ww.gt.1.d-5) then
                  c3=(1.d0-c1)/ww**2
               else
                  c3=0.5d0
               endif
!
!              rotation matrix c
!
               do i=1,3
                  do j=1,3
                     c(i,j)=c1*d(i,j)+
     &                  c2*(e(i,1,j)*w(1)+e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                  c3*w(i)*w(j)
                  enddo
               enddo
!
               c4=-c2
               if(ww.gt.0.00464159) then
                  c5=(ww*dcos(ww)-dsin(ww))/ww**3
               else
                  c5=-1.d0/3.d0
               endif
               if(ww.gt.0.0031623) then
                  c6=(ww*dsin(ww)-2.d0+2.d0*dcos(ww))/ww**4
               else
                  c6=-1.d0/12.d0
               endif
!
!              derivative of the rotation matrix c with respect to
!              the rotation vector w
!
               do i=1,3
                  do j=1,3
                     do k=1,3
                        dc(i,j,k)=c4*w(k)*d(i,j)+
     &                            c5*w(k)*(e(i,1,j)*w(1)+
     &                                     e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                            c2*e(i,k,j)+
     &                            c6*w(k)*w(i)*w(j)+
     &                            c3*(d(i,k)*w(j)+d(j,k)*w(i))
                     enddo
                  enddo
               enddo
!
!              dummy variable
!
               do i=1,3
                  do j=1,3
c                     f(i,j)=c(i,j)-d(i,j)-dc(i,j,1)*w(1)-dc(i,j,2)*w(2)-
c     &                                    dc(i,j,3)*w(3)
                     f(i,j)=expan*c(i,j)-d(i,j)
                  enddo
               enddo
            endif
!
            coefmpc(indexexp)=c(idir,1)*(co(1,irefnode)-co(1,inode))+
     &           c(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           c(idir,3)*(co(3,irefnode)-co(3,inode))
!
!           determining the coefficients of the rotational degrees
!           of freedom
!
            coefmpc(index)=(dc(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,1)*(co(3,irefnode)-co(3,inode)))*expan
!
            index=nodempc(3,index)
            coefmpc(index)=(dc(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,2)*(co(3,irefnode)-co(3,inode)))*expan
!
            index=nodempc(3,index)
            coefmpc(index)=(dc(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,3)*(co(3,irefnode)-co(3,inode)))*expan
!
!           determining the nonhomogeneous part
!
            index=nodempc(3,index)
            coefmpc(index)=1.d0
!
!           old value of the nonhomogeneous term must be zero
!
            vold(nodempc(2,index),nodempc(1,index))=0.d0
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=f(idir,1)*(co(1,irefnode)-co(1,inode))+
     &           f(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           f(idir,3)*(co(3,irefnode)-co(3,inode))-
     &           vold(idir,irefnode)+vold(idir,inode)
!
         elseif(labmpc(ii)(1:8).eq.'STRAIGHT') then
!
!           determining nodes and directions involved in MPC
!            
            index=ipompc(ii)
            np=nodempc(1,index)
            j=nodempc(2,index)
            index=nodempc(3,index)
            i=nodempc(2,index)
            index=nodempc(3,index)
            na=nodempc(1,index)
            index=nodempc(3,nodempc(3,index))
            nb=nodempc(1,index)
!
!           determining the coefficients
!
            index=ipompc(ii)
            c2=co(i,na)+vold(i,na)-co(i,nb)-vold(i,nb)
            if(dabs(c2).lt.1.d-5) then
               write(*,*) '*WARNING in nonlinmpc: coefficient of'
               write(*,*)
     &     '          dependent node in STRAIGHT MPC is zero'
               idofrem=8*(np-1)+j
!
!              determining a new dependent term       
!               
               ddmax=abs(c2)
               l=i
               m=j
               do k=1,2
                  l=l+1
                  m=m+1
                  if(l.gt.3) l=l-3
                  if(m.gt.3) m=m-3
                  dd=dabs(co(l,na)+vold(l,na)-co(l,nb)-vold(l,nb))
                  if(dd.gt.ddmax) then
                     ddmax=dd
                     lmax=l
                     mmax=m
                  endif
               enddo
               i=lmax
               j=mmax
               idofins=8*(np-1)+j
!
               call changedepterm(ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
!
               index=ipompc(ii)
               nodempc(2,index)=j
               index=nodempc(3,index)
               nodempc(2,index)=i
               index=nodempc(3,index)
               nodempc(2,index)=j
               index=nodempc(3,index)
               nodempc(2,index)=i
               index=nodempc(3,index)
               nodempc(2,index)=j
               index=nodempc(3,index)
               nodempc(2,index)=i
               index=nodempc(3,index)
               nodempc(2,index)=j
               if(icascade.eq.0) icascade=1
               c2=co(i,na)+vold(i,na)-co(i,nb)-vold(i,nb)
            endif
            coefmpc(index)=c2
            index=nodempc(3,index)
            c3=co(j,nb)+vold(j,nb)-co(j,na)-vold(j,na)
            coefmpc(index)=c3
            index=nodempc(3,index)
            c5=co(i,nb)+vold(i,nb)-co(i,np)-vold(i,np)
            coefmpc(index)=c5
            index=nodempc(3,index)
            c6=co(j,np)+vold(j,np)-co(j,nb)-vold(j,nb)
            coefmpc(index)=c6
            index=nodempc(3,index)
            c4=co(i,np)+vold(i,np)-co(i,na)-vold(i,na)
            coefmpc(index)=c4
            index=nodempc(3,index)
            c1=co(j,na)+vold(j,na)-co(j,np)-vold(j,np)
            coefmpc(index)=c1
            index=nodempc(3,index)
!
!           nonhomogeneous term
!
            coefmpc(index)=1.d0
!
!           old value of the nonhomogeneous term must be zero
!
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=-c1*c2+c3*c4
            if(newstep.eq.1) xboun(ilboun(id))=xbounact(ilboun(id))
            vold(nodempc(2,index),nodempc(1,index))=
     &          (1.d0-reltime)*xboun(ilboun(id))
         elseif(labmpc(ii)(1:5).eq.'PLANE') then
!
!           determining nodes and directions involved in MPC
!            
            index=ipompc(ii)
            np=nodempc(1,index)
            i1=nodempc(2,index)
            index=nodempc(3,index)
            i2=nodempc(2,index)
            index=nodempc(3,index)
            i3=nodempc(2,index)
            index=nodempc(3,index)
            na=nodempc(1,index)
            index=nodempc(3,nodempc(3,nodempc(3,index)))
            nb=nodempc(1,index)
            index=nodempc(3,nodempc(3,nodempc(3,index)))
            nc=nodempc(1,index)
!
!           determining the coefficients
!
            a11=co(i1,np)+vold(i1,np)-co(i1,nc)-vold(i1,nc)
            a12=co(i2,np)+vold(i2,np)-co(i2,nc)-vold(i2,nc)
            a13=co(i3,np)+vold(i3,np)-co(i3,nc)-vold(i3,nc)
            a21=co(i1,na)+vold(i1,na)-co(i1,nc)-vold(i1,nc)
            a22=co(i2,na)+vold(i2,na)-co(i2,nc)-vold(i2,nc)
            a23=co(i3,na)+vold(i3,na)-co(i3,nc)-vold(i3,nc)
            a31=co(i1,nb)+vold(i1,nb)-co(i1,nc)-vold(i1,nc)
            a32=co(i2,nb)+vold(i2,nb)-co(i2,nc)-vold(i2,nc)
            a33=co(i3,nb)+vold(i3,nb)-co(i3,nc)-vold(i3,nc)
!
            b11=a22*a33-a23*a32
            b12=a31*a23-a21*a33
            b13=a21*a32-a31*a22
            b21=a32*a13-a12*a33
            b22=a11*a33-a31*a13
            b23=a31*a12-a11*a32
            b31=a12*a23-a22*a13
            b32=a21*a13-a11*a23
            b33=a11*a22-a12*a21
!
            index=ipompc(ii)
            if(dabs(b11).lt.1.d-5) then
               write(*,*) '*WARNING in nonlinmpc: coefficient of'
               write(*,*) '         dependent node in PLANE MPC is zero'
!
               idofrem=8*(nodempc(1,index)-1)+i1
!
               if(dabs(b12).gt.dabs(b13)) then
                  idofins=8*(nodempc(1,index)-1)+i2
                  call changedepterm
     &                   (ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
                  coefmpc(index)=b12
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=b11
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b13
                  index=nodempc(3,index)
                  coefmpc(index)=b22
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=b21
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b23
                  index=nodempc(3,index)
                  coefmpc(index)=b32
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=b31
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b33
                  index=nodempc(3,index)
                  coefmpc(index)=-b12-b22-b32
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=-b11-b21-b31
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=-b13-b23-b33
                  if(icascade.eq.0) icascade=1
               else
                  idofins=8*(nodempc(1,index)-1)+i3
                  call changedepterm
     &                  (ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
                  coefmpc(index)=b13
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=b12
                  index=nodempc(3,index)
                  coefmpc(index)=b11
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b23
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=b22
                  index=nodempc(3,index)
                  coefmpc(index)=b21
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b33
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=b32
                  index=nodempc(3,index)
                  coefmpc(index)=b31
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=-b13-b23-b33
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=-b12-b22-b32
                  index=nodempc(3,index)
                  coefmpc(index)=-b11-b21-b31
                  nodempc(2,index)=i1
                  if(icascade.eq.0) icascade=1
               endif
            else
               coefmpc(index)=b11
               index=nodempc(3,index)
               coefmpc(index)=b12
               index=nodempc(3,index)
               coefmpc(index)=b13
               index=nodempc(3,index)
               coefmpc(index)=b21
               index=nodempc(3,index)
               coefmpc(index)=b22
               index=nodempc(3,index)
               coefmpc(index)=b23
               index=nodempc(3,index)
               coefmpc(index)=b31
               index=nodempc(3,index)
               coefmpc(index)=b32
               index=nodempc(3,index)
               coefmpc(index)=b33
               index=nodempc(3,index)
               coefmpc(index)=-b11-b21-b31
               index=nodempc(3,index)
               coefmpc(index)=-b12-b22-b32
               index=nodempc(3,index)
               coefmpc(index)=-b13-b23-b33
            endif
            index=nodempc(3,index)
            coefmpc(index)=1.d0
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
!
!           old value of the nonhomogeneous term must be zero
!
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=a11*b11+a12*b12+a13*b13
            if(newstep.eq.1) xboun(ilboun(id))=xbounact(ilboun(id))
            vold(nodempc(2,index),nodempc(1,index))=0.d0
         elseif(labmpc(ii)(1:9).eq.'ISOCHORIC') then
            isochoric=.true.
!
!           next segment is deactivated (CYCLID instead of CYCLIC):
!           cylic MPC's are considered to be linear
!
         elseif((labmpc(ii)(1:6).eq.'CYCLID').and.(ithermal.ne.2)) then
            index=ipompc(ii)
            noded=nodempc(1,index)
!
!           check for thermal MPC
!
            if(nodempc(2,index).eq.0) cycle loop
!
!           check whether the next two MPC's are cyclic MPC's
!           applied to the same dependent node
!            
            if((nodempc(1,ipompc(ii+1)).ne.noded).or.
     &         (labmpc(ii+1)(1:6).ne.'CYCLIC').or.
     &         (nodempc(1,ipompc(ii+2)).ne.noded).or.
     &         (labmpc(ii+2)(1:6).ne.'CYCLIC')) then
               write(*,*) '*WARNING in nonlinmpc: no three'
               write(*,*) '         cyclic MPCs pertaining'
               write(*,*) '         to the same dependent node;'
               write(*,*) '         no update'
               cycle loop
            endif
!
!           finding the cyclic symmetry axis
!
            do i=1,ntrans
               if(trab(7,i).eq.2) exit
            enddo
            if(i.gt.ntrans) then
               write(*,*) '*ERROR in nonlinmpc: cyclic symmetry'
               write(*,*) '       axis not found'
               stop
            endif
            do j=1,6
               csab(j)=trab(j,i)
            enddo
            csab(7)=-1
!
!           determining the independent node
!                     
            nodei=0
            do
               if(nodempc(1,index).ne.noded) then
                  if(nodei.eq.0) then
                     nodei=nodempc(1,index)
                  elseif(nodei.ne.nodempc(1,index)) then
                     write(*,*) '*WARNING in nonlinmpc:'
                     write(*,*) '          cyclic symmetry conditions'
                     write(*,*) '          between unequal meshes'
                     write(*,*) '          no update'
                     cycle loop
                  endif
               endif
               index=nodempc(3,index)
               if(index.eq.0) then
                  if(nodei.eq.0) then
                     write(*,*) '*ERROR in nonlinmpc:'
                     write(*,*) '       no independent node found'
                     stop
                  else
                     exit
                  endif
               endif
            enddo
!
!           actual location of dependent and independent node
!
            do i=1,3
               pd(i)=co(i,noded)+vold(i,noded)
               pi(i)=co(i,nodei)+vold(i,nodei)
            enddo
!
!           update transformation matrix
!
            call transformatrix(csab,pd,ad)
            call transformatrix(csab,pi,ai)
!     
!     checking for latin hypercube positions in matrix al none of
!     which are zero
!     
            do inum=1,6
               if((dabs(ad(lathyp(1,inum),1)).gt.1.d-3).and.
     &            (dabs(ad(lathyp(2,inum),2)).gt.1.d-3).and.
     &            (dabs(ad(lathyp(3,inum),3)).gt.1.d-3)) exit
            enddo
!
!           remove old DOFs
!
            do j=1,3
               idof=8*(noded-1)+j
               call nident(ikmpc,idof,nmpc,id)
               if(id.lt.0) then
                  write(*,*) '*ERROR in nonlinmpc: error in'
                  write(*,*) '       MPC database'
                  stop
               elseif(ikmpc(id).ne.idof) then
                  write(*,*) '*ERROR in nonlinmpc: error in'
                  write(*,*) '       MPC database'
                  stop
               endif
!
               do k=id,nmpc-1
                  ikmpc(k)=ikmpc(k+1)
                  ilmpc(k)=ilmpc(k+1)
               enddo
            enddo
!
!           add new MPCs
!
            ii=ii-1
            do ndir=1,3
               ii=ii+1
               number=lathyp(ndir,inum)
               idof=8*(noded-1)+number
               call nident(ikmpc,idof,nmpc-1,id)
               if(id.gt.0) then
                  if(ikmpc(id).eq.idof) then
                     write(*,*) '*WARNING in generatecycmpcs: cyclic MPC
     & in node'
                     write(*,*) '         ',noded,' and direction ',ndir
                     write(*,*) '         cannot be created: the'
                     write(*,*) '         DOF in this node is already us
     &ed'
                     cycle
                  endif
               endif
               number=number-1
!     
!     updating ikmpc and ilmpc
!     
               do j=nmpc,id+2,-1
                  ikmpc(j)=ikmpc(j-1)
                  ilmpc(j)=ilmpc(j-1)
               enddo
               ikmpc(id+1)=idof
               ilmpc(id+1)=nmpc
!
!              update the MPC coefficients
!
               index=ipompc(ii)
               do j=1,3
                  number=number+1
                  if(number.gt.3) number=1
                  if(dabs(ad(number,ndir)).lt.1.d-5) cycle
                  if(index.eq.0) then
                     write(*,*)'*ERROR in nonlinmpc: index=0'
                     stop
                  endif
                  nodempc(1,index)=noded
                  nodempc(2,index)=number
                  coefmpc(index)=ad(number,ndir)
                  index=nodempc(3,index)
               enddo
               do j=1,3
                  number=number+1
                  if(number.gt.3) number=1
                  if(dabs(ai(number,ndir)).lt.1.d-5) cycle
                  if(index.eq.0) then
                     write(*,*)'*ERROR in nonlinmpc: index=0'
                     stop
                  endif
                  nodempc(1,index)=nodei
                  nodempc(2,index)=number
                  coefmpc(index)=-ai(number,ndir)
                  index=nodempc(3,index)
               enddo
            enddo
         elseif((labmpc(ii)(1:20).ne.'                    ').and.
     &          (labmpc(ii)(1:7).ne.'CONTACT').and.
     &          (labmpc(ii)(1:6).ne.'CYCLIC').and.
     &          (labmpc(ii)(1:9).ne.'SUBCYCLIC')) then
            index=ipompc(ii)
            i=0
            do
               if(index.eq.0) exit
               node=nodempc(1,index)
               i=i+1
               iaux(i)=nodempc(2,index)
               aux(6*maxlenmpc+i)=coefmpc(index)
               do j=1,3
                  aux(3*(i-1)+j)=co(j,node)
                  aux(3*(maxlenmpc+i-1)+j)=vold(j,node)
               enddo
               index=nodempc(3,index)
            enddo
            n=i-1
            if((labmpc(ii)(1:7).eq.'MEANROT').or.
     &         (labmpc(ii)(1:1).eq.'1')) then
               call umpc_mean_rot(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            elseif(labmpc(ii)(1:4).eq.'DIST') then
               call umpc_dist(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            elseif(labmpc(ii)(1:3).eq.'GAP') then
               call umpc_gap(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            elseif(labmpc(ii)(1:4).eq.'USER') then
               call umpc_user(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            else
               write(*,*) '*ERROR in nonlinmpc: mpc of type',labmpc(ii)
               write(*,*) '       is unknown'
               stop
            endif
            index=ipompc(ii)
!
            if(iaux(1).ne.nodempc(2,index)) then
!
!              dependent MPC has changed
!  
               idofrem=8*(nodempc(1,index)-1)+nodempc(2,index)
               idofins=8*(nodempc(1,index)-1)+iaux(1)
               call changedepterm(ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
               if(icascade.eq.0) icascade=1
            endif
!
            i=0
            do
               if(index.eq.0) exit
               i=i+1
               if(i.le.n) then
!
!                 check whether any directions have changed:
!                 necessitates calling of remastruct
!
                  if(iaux(i).ne.nodempc(2,index)) then
                     if(icascade.eq.0) icascade=1
                  endif
                  nodempc(2,index)=iaux(i)
                  coefmpc(index)=aux(6*maxlenmpc+i)
               else
                  coefmpc(index)=1.d0
!
!                 old value of the nonhomogeneous term must be zero
!
                  vold(nodempc(2,index),nodempc(1,index))=0.d0
                  idof=8*(nodempc(1,index)-1)+nodempc(2,index)
                  call nident(ikboun,idof,nboun,id)
                  xbounact(ilboun(id))=const
               endif
               index=nodempc(3,index)
            enddo
         endif
      enddo loop
!
!     incompressible material
!
      if(.not.isochoric) return
!
!     initialization of the mpc's
!
      nmpc01=0
      nmpcdif=0
      do i=1,nmpc
         if(labmpc(i)(1:9).eq.'ISOCHORIC') then
            if(nmpc01.eq.0) nmpc01=i
            nmpcdif=i
            index=ipompc(i)
            do
               if(nodempc(3,index).eq.0) then
                  idof=8*(nodempc(1,index)-1)+nodempc(2,index)
                  call nident(ikboun,idof,nboun,id)
                  xbounact(ilboun(id))=0.d0
                  exit
               endif
               coefmpc(index)=0.d0
               index=nodempc(3,index)
            enddo
         endif
      enddo
      nmpc0=nmpc01-1
      nmpcdif=nmpcdif-nmpc0
!
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:7).eq.'C3D20RI') then
            indexe=ipkon(i)
!
            do j=1,20
               node=kon(indexe+j)
               do k=1,3
                  xlag(k,j)=co(k,node)
                  xeul(k,j)=xlag(k,j)+vold(k,node)
               enddo
            enddo
!
            do j=1,8
               mpc=0
               node=kon(indexe+j)
               label(1:9)='ISOCHORIC'
               write(label(10:20),'(i11)') node
c               write(*,*) 'nonlinmpclab ',label
               call cident20(labmpc(nmpc01),label,nmpcdif,id)
               id=id+nmpc0
c               write(*,*) 'nonlinmpclab ',id,label,labmpc(id)
               if(id.gt.0) then
                  if(labmpc(id).eq.label) then
                     mpc=id
                  endif
               endif
               if(mpc.eq.0) cycle
!
               xi=coloc(1,j)
               et=coloc(2,j)
               ze=coloc(3,j)
!
               call deuldlag(xi,et,ze,xlag,xeul,xj,a)
!
               b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
               b(1,2)=a(3,1)*a(2,3)-a(2,1)*a(3,3)
               b(1,3)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
               b(2,1)=a(3,2)*a(1,3)-a(1,2)*a(3,3)
               b(2,2)=a(1,1)*a(3,3)-a(3,1)*a(1,3)
               b(2,3)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
               b(3,1)=a(1,2)*a(2,3)-a(2,2)*a(1,3)
               b(3,2)=a(2,1)*a(1,3)-a(1,1)*a(2,3)
               b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
!
               index=ipompc(mpc)
               do
                  if(nodempc(3,index).eq.0) then
                     coefmpc(index)=1.d0
                     idof=8*(nodempc(1,index)-1)+nodempc(2,index)
                     call nident(ikboun,idof,nboun,id)
                     xbounact(ilboun(id))=xbounact(ilboun(id))+
     &                    a(1,1)*b(1,1)+a(1,2)*b(1,2)+a(1,3)*b(1,3)
     &                    -1.d0/xj
c                     write(*,*) 'nonlinmpcboun ',nodempc(1,index),
c     &                        nodempc(2,index),ilboun(id),
c     &                       xbounact(ilboun(id))
                     exit
                  else
                     node=nodempc(1,index)
                     idir=nodempc(2,index)
                     do k=1,7
                        if(kon(indexe+neigh(k,j)).eq.node) then
                           if(k.eq.1) then
                              if(idir.eq.1) then
                                 coefmpc(index)=coefmpc(index)+
     &                             1.5d0*(xi*b(1,1)+et*b(1,2)+ze*b(1,3))
                              elseif(idir.eq.2) then
                                 coefmpc(index)=coefmpc(index)+
     &                             1.5d0*(xi*b(2,1)+et*b(2,2)+ze*b(2,3))
                              elseif(idir.eq.3) then
                                 coefmpc(index)=coefmpc(index)+
     &                             1.5d0*(xi*b(3,1)+et*b(3,2)+ze*b(3,3))
                              endif
                           elseif(k.eq.2) then
                              if(idir.eq.1) then
                               coefmpc(index)=coefmpc(index)-
     &                                        2.d0*xi*b(1,1)
                              elseif(idir.eq.2) then
                               coefmpc(index)=coefmpc(index)-
     &                                        2.d0*xi*b(2,1)
                              elseif(idir.eq.3) then
                               coefmpc(index)=coefmpc(index)-
     &                                        2.d0*xi*b(3,1)
                              endif
                           elseif(k.eq.3) then
                              if(idir.eq.1) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*xi*b(1,1)
                              elseif(idir.eq.2) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*xi*b(2,1)
                              elseif(idir.eq.3) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*xi*b(3,1)
                              endif
                           elseif(k.eq.4) then
                              if(idir.eq.1) then
                               coefmpc(index)=coefmpc(index)-
     &                                        2.d0*et*b(1,2)
                              elseif(idir.eq.2) then
                               coefmpc(index)=coefmpc(index)-
     &                                        2.d0*et*b(2,2)
                              elseif(idir.eq.3) then
                               coefmpc(index)=coefmpc(index)-
     &                                        2.d0*et*b(3,2)
                              endif
                           elseif(k.eq.5) then
                              if(idir.eq.1) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*et*b(1,2)
                              elseif(idir.eq.2) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*et*b(2,2)
                              elseif(idir.eq.3) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*et*b(3,2)
                              endif
                           elseif(k.eq.6) then
                              if(idir.eq.1) then
                              coefmpc(index)=coefmpc(index)-
     &                                       2.d0*ze*b(1,3)
                              elseif(idir.eq.2) then
                              coefmpc(index)=coefmpc(index)-
     &                                       2.d0*ze*b(2,3)
                              elseif(idir.eq.3) then
                              coefmpc(index)=coefmpc(index)-
     &                                       2.d0*ze*b(3,3)
                              endif
                           elseif(k.eq.7) then
                              if(idir.eq.1) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*ze*b(1,3)
                              elseif(idir.eq.2) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*ze*b(2,3)
                              elseif(idir.eq.3) then
                              coefmpc(index)=coefmpc(index)+
     &                                       0.5d0*ze*b(3,3)
                              endif
                           endif
                           exit
                        endif
                     enddo
                  endif
                  index=nodempc(3,index)
               enddo
!
            enddo
         endif
      enddo
!
      return
      end






