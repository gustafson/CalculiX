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
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine radmatrix(ntr,
     &     acr,bcr,sideload,nelemload,xloadact,lakon,vold,
     &     ipkon,kon,co,pmid,e1,e2,e3,iptri,
     &     kontri,ntri,nloadtr,tarea,tenv,physcon,erad,f,
     &     dist,idist,area,ithermal,iinc,iit,
     &     cs,mcs,inocs,ntrit,nk,fenv,istep,dtime,ttime,
     &     time,iviewfile,jobnamef,xloadold,reltime,nmethod,mi,
     &     iemchange,nam,iamload)
!     
      implicit none
!     
      logical covered(160,160),exi
!
!     change following line if nlabel is increased
!     
      character*87 label(30)
      character*8 lakonl,lakon(*)
      character*20 sideload(*)
      character*132 jobnamef(*),fnvw
!     
      integer ntr,nelemload(2,*),nope,nopes,mint2d,i,j,k,l,
     &     node,ifaceq(8,6),ifacet(6,4),iviewfile,mi(2),
     &     ifacew(8,5),nelem,ig,index,konl(20),iflag,
     &     ipkon(*),kon(*),ncovered,kontri(3,*),iptri(*),nloadtr(*),
     &     i1,j1,istart,iend,jstart,jend,imin,imid,imax,mcs,inocs(*),
     &     k1,kflag,idist(*),ndist,i2,i3,ng,idi,idj,ntri,
     &     ithermal,iinc,iit,ix,iy,ntrit,jj,is,m,jmod,nkt,
     &     icntrl,imag,nk,istep,jltyp,nfield,nonzero,nmethod,
     &     limev,ier,nw,idata(1),ncalls,nlabel,iemchange,nam,
     &     iamload(2,*)
!     
      real*8 acr(ntr,*),bcr(ntr,1),xloadact(2,*),h(2),w(239),
     &     xl2(3,8),coords(3),dxsj2,temp,xi,et,weight,xsj2(3),
     &     vold(0:mi(2),*),co(3,*),shp2(7,8),xs2(3,7),xn(3),xxn,
     &     pmid(3,*),e3(4,*),e1(3,*),e2(3,*),p1(3),p2(3),p3(3),
     &     areamean,tarea(*),tenv(*),x,y,cs(17,*),porigin(3),
     &     erad(*),fenv(*),e,ec,physcon(*),yymin,yymax,xxmin,
     &     xxmid,xxmax,dummy,a(3,3),b(3,3),c(3,3),ddd(3),p31(3),
     &     xx(3),yy(3),ftij,f(ntr,*),dint,dir(3),tl2(8),
     &     dirloc(3),dist(*),area(*),dd,p21(3),p32(3),pi,
     &     totarea,fn,stn,qfn,een,t(3),sidemean,tvar(2),field,
     &     dtime,ttime,time,areaj,xloadold(2,*),reltime,p(3,3),
     &     fform,ver(2,3),epsabs,epsrel,abserr,vj(3,3),unitvec(3,3),
     &     rdata(1),vertex(3,3),vertexl(2,3),factor,argument
!     
      include "gauss.f"
!     
      data ifaceq /4,3,2,1,11,10,9,12,
     &     5,6,7,8,13,14,15,16,
     &     1,2,6,5,9,18,13,17,
     &     2,3,7,6,10,19,14,18,
     &     3,4,8,7,11,20,15,19,
     &     4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &     1,2,4,5,9,8,
     &     2,3,4,6,10,9,
     &     1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &     4,5,6,10,11,12,0,0,
     &     1,2,5,4,7,14,10,13,
     &     2,3,6,5,8,15,11,14,
     &     4,6,3,1,12,15,9,13/
      data iflag /2/
!     
      common /formfactor/ vj,unitvec,porigin
!
      external fform
!
!     change following line if nlabel is increased and the dimension
!     of field label above!
!     
      nlabel=30
!
!     factor determines when the numerical integration using cubtri
!     is replaced by a simplified formula using only the center
!     of gravity of one of the triangles. The integration over the
!     other triangle is exact (analytical formula, see
!     "Radiosity: a Programmer's Perspective", by Ian Ashdown, Wiley, 1994)
!     If the distance between the center of gravity of the triangles
!     is larger then factor*the projected sqrt(area) of the triangle on the 
!     hemisphere, the simplified formula is taken
!
      factor=0.d0
!
      pi=4.d0*datan(1.d0)
!     
      tvar(1)=time
      tvar(2)=ttime+dtime
!     
!     cavity radiation!            
!     
!     the default sink temperature is updated at the start of each
!     increment
!     
      do i=1,ntr
         node=nelemload(2,nloadtr(i))
         if(node.ne.0) then
            tenv(i)=vold(0,node)-physcon(1)
         elseif(iit.le.0) then
            tenv(i)=xloadact(2,nloadtr(i))-physcon(1)
         endif
      enddo
!     
!     for pure thermal steps the viewfactors have to be
!     calculated only once, for thermo-mechanical steps
!     (ithermal=3) they are recalculated in each iteration
!     unless they are read from file
!     
      if(((ithermal.eq.3).and.(iviewfile.ge.0)).or.(iit.eq.-1)) then
         if(iviewfile.lt.0) then
            if(ithermal.eq.3) then
               write(*,*) '*WARNING in radmatrix: viewfactors are being'
               write(*,*) '         read from file for a thermomechani-'
               write(*,*) '         cal calculation: they will not be '
               write(*,*) '         recalculated in every iteration.'
            endif
!     
            write(*,*) 'Reading the viewfactors from file'
            write(*,*)
!
            if(jobnamef(2)(1:1).eq.' ') then
               do i=1,132
                  if(jobnamef(1)(i:i).eq.' ') exit
               enddo
               i=i-1
               fnvw=jobnamef(1)(1:i)//'.vwf'
            else
               fnvw=jobnamef(2)
            endif
            inquire(file=fnvw,exist=exi)
            if(exi) then
               open(10,file=fnvw,status='old',form='unformatted',
     &              access='sequential',err=10)
            else
               write(*,*) '*ERROR in radmatrix: viewfactor file ',fnvw
               write(*,*) 'does not exist'
               stop
            endif
!
            read(10) nonzero
            do k=1,nonzero
               read(10) i,j,f(i,j)
            enddo
            read(10)(fenv(i),i=1,ntr)
!
            close(10)
         else
!     
         write(*,*) 'Calculating the viewfactors'
         write(*,*)
!     
         ng=160
         dint=2.d0/ng
!     
!     updating the displacements for cyclic symmetric structures
!     
         if(mcs.gt.0) then
            nkt=0
            do i=1,mcs
               if(int(cs(1,i)).gt.nkt) nkt=int(cs(1,i))
            enddo
            nkt=nk*nkt
            do i=1,nlabel
               do l=1,87
                  label(i)(l:l)=' '
               enddo
            enddo
            label(1)(1:1)='U'
            imag=0
            icntrl=2
            call rectcyl(co,vold,fn,stn,qfn,een,cs,nk,icntrl,t,
     &           label,imag,mi)
            
            do jj=0,mcs-1
               is=cs(1,jj+1)
!               
               do i=1,is-1
                  do l=1,nk
                     if(inocs(l).ne.jj) cycle
                     do m=1,mi(2)
                        vold(m,l+nk*i)=vold(m,l)
                     enddo
                  enddo
               enddo
            enddo
            icntrl=-2
            call rectcyl(co,vold,fn,stn,qfn,een,cs,nkt,icntrl,t,
     &           label,imag,mi)
         endif
!     
!     calculating the momentaneous center of the triangles,
!     area of the triangles and normal to the triangles
!     
         sidemean=0.d0
         do i=1,ntrit
            i1=kontri(1,i)
            if(i1.eq.0) cycle
            i2=kontri(2,i)
            i3=kontri(3,i)
            do j=1,3
               p1(j)=co(j,i1)+vold(j,i1)
               p2(j)=co(j,i2)+vold(j,i2)
               p3(j)=co(j,i3)+vold(j,i3)
               pmid(j,i)=(p1(j)+p2(j)+p3(j))/3.d0
               p21(j)=p2(j)-p1(j)
               p32(j)=p3(j)-p2(j)
            enddo
!     
!           normal to the triangle
!     
            e3(1,i)=p21(2)*p32(3)-p32(2)*p21(3)
            e3(2,i)=p21(3)*p32(1)-p32(3)*p21(1)
            e3(3,i)=p21(1)*p32(2)-p32(1)*p21(2)
!     
            dd=dsqrt(e3(1,i)*e3(1,i)+e3(2,i)*e3(2,i)+
     &           e3(3,i)*e3(3,i))
!
!           check for degenerated triangles
!
            if(dd.lt.1.d-20) then
               area(i)=0.d0
               cycle
            endif
!     
            do j=1,3
               e3(j,i)=e3(j,i)/dd
            enddo
!     
!           area of the triangle
!     
            area(i)=dd/2.d0
!     
!           unit vector parallel to side 1-2
!     
            dd=dsqrt(p21(1)*p21(1)+p21(2)*p21(2)+p21(3)*p21(3))
            sidemean=sidemean+dd
            do j=1,3
               e1(j,i)=p21(j)/dd
            enddo
!     
!           unit vector orthogonal to e1 and e3
!     
            e2(1,i)=e3(2,i)*e1(3,i)-e3(3,i)*e1(2,i)
            e2(2,i)=e3(3,i)*e1(1,i)-e3(1,i)*e1(3,i)
            e2(3,i)=e3(1,i)*e1(2,i)-e3(2,i)*e1(1,i)
!     
!           the fourth component in e3 is the constant term in the
!           equation of the triangle plane in the form
!           e3(1)*x+e3(2)*y+e3(3)*z+e3(4)=0
!     
            e3(4,i)=-(e3(1,i)*p1(1)+e3(2,i)*p1(2)
     &           +e3(3,i)*p1(3))
         enddo
         sidemean=sidemean/ntrit
!     
!        determine the geometrical factors
!     
!        initialization of the fields
!     
         do i=1,ntr
            do j=1,ntr
               f(i,j)=0.d0
            enddo
         enddo
!     
         do i=1,ntri
            if(area(i).lt.1.d-20) cycle
!
!     vertices of triangle i in local coordinates
!     
            i1=kontri(1,i)
            if(i1.eq.0) cycle
            i2=kontri(2,i)
            i3=kontri(3,i)
            do j=1,3
               porigin(j)=co(j,i1)+vold(j,i1)
               p2(j)=co(j,i2)+vold(j,i2)
               p3(j)=co(j,i3)+vold(j,i3)
               p21(j)=p2(j)-porigin(j)
               p31(j)=p3(j)-porigin(j)
            enddo
            ver(1,1)=0.d0
            ver(2,1)=0.d0
            ver(1,2)=dsqrt(p21(1)**2+p21(2)**2+p21(3)**2)
            ver(2,2)=0.d0
            ver(1,3)=p31(1)*e1(1,i)+p31(2)*e1(2,i)+p31(3)*e1(3,i)
            ver(2,3)=p31(1)*e2(1,i)+p31(2)*e2(2,i)+p31(3)*e2(3,i)
!
            do k=1,3
               unitvec(k,1)=e1(k,i)
               unitvec(k,2)=e2(k,i)
               unitvec(k,3)=e3(k,i)
            enddo
!     
!           checking which triangles face triangle i
!     
            ndist=0
            call nident(iptri,i,ntr,idi)
            do j=1,ntrit
               if((kontri(1,j).eq.0).or.(area(j).lt.1.d-20)) cycle
               if(pmid(1,j)*e3(1,i)+pmid(2,j)*e3(2,i)+
     &              pmid(3,j)*e3(3,i)+e3(4,i).le.sidemean/800.d0) cycle
               if(pmid(1,i)*e3(1,j)+pmid(2,i)*e3(2,j)+
     &              pmid(3,i)*e3(3,j)+e3(4,j).le.sidemean/800.d0) cycle
!
               if(j.gt.ntri) then
                  jmod=mod(j,ntri)
                  if(jmod.eq.0) jmod=ntri
               else
                  jmod=j
               endif
!     
c               call nident(iptri,i,ntr,idi)
               call nident(iptri,jmod,ntr,idj)
               if(sideload(nloadtr(idi))(18:20).ne.
     &            sideload(nloadtr(idj))(18:20)) cycle
!
               ndist=ndist+1
               dist(ndist)=dsqrt((pmid(1,j)-pmid(1,i))**2+
     &              (pmid(2,j)-pmid(2,i))**2+
     &              (pmid(3,j)-pmid(3,i))**2)
               idist(ndist)=j
            enddo
            if(ndist.eq.0) cycle
!     
!           ordering the triangles
!     
            kflag=2
            call dsort(dist,idist,ndist,kflag)
!     
!           initializing the coverage matrix
!
c            write(*,*) i,(idist(i1),i1=1,ndist)
            ncovered=0
            do i1=1,ng
               x=((i1-0.5d0)*dint-1.d0)**2
               do j1=1,ng
                  y=((j1-0.5d0)*dint-1.d0)**2
                  if(x+y.gt.1.d0) then
                     covered(i1,j1)=.true.
                     ncovered=ncovered+1
                  else
                     covered(i1,j1)=.false.
                  endif
               enddo
            enddo
!     
            do k1=1,ndist
               j=idist(k1)
!
!              determining the 2-D projection of the vertices
!              of triangle j
!
               do l=1,3
                  do k=1,3
                     vertex(k,l)=co(k,kontri(l,j))-pmid(k,i)
                  enddo
                  dd=dsqrt(vertex(1,l)**2+vertex(2,l)**2+
     &                     vertex(3,l)**2)
                  do k=1,3
                     vertex(k,l)=vertex(k,l)/dd
                  enddo
                  vertexl(1,l)=vertex(1,l)*e1(1,i)+
     &                         vertex(2,l)*e1(2,i)+
     &                         vertex(3,l)*e1(3,i)
                  vertexl(2,l)=vertex(1,l)*e2(1,i)+
     &                         vertex(2,l)*e2(2,i)+
     &                         vertex(3,l)*e2(3,i)
               enddo
!
!              determining the center of gravity of the projected
!              triangle
!
               do k=1,2
                  dirloc(k)=(vertexl(k,1)+vertexl(k,2)+
     &                       vertexl(k,3))/3.d0
               enddo
!     
!              determine the direction vector in global coordinates
!     
               do k=1,3
                  dir(k)=(pmid(k,j)-pmid(k,i))/dist(k1)
               enddo
!     
!              direction vector in local coordinates of triangle i
!     
               dirloc(3)=dir(1)*e3(1,i)+dir(2)*e3(2,i)+dir(3)*e3(3,i)
!     
!     check whether this direction was already covered
!     
               ix=int((dirloc(1)+1.d0)/dint)+1
               iy=int((dirloc(2)+1.d0)/dint)+1
               if(covered(ix,iy)) then
c                  write(*,*) 'triangle ',j,' was already covered'
                  cycle
               endif
!
!              if surfaces are close, numerical integration with
!              cubtri is performed
!
               if(dist(k1).le.factor*dsqrt(area(i))*dirloc(3)) then
!     
!     vertices of triangle j
!     
                  do k=1,3
                     do l=1,3
                        vj(l,k)=co(l,kontri(k,j))+vold(l,kontri(k,j))
                     enddo
                  enddo
!     
!     formfactor contribution
!     
                  epsrel=0.01d0
                  epsabs=0.d0
                  limev=100
                  nw=239
                  ncalls=0
!     
!     max 1000 evaluations for nw=239
!     
                  call cubtri(fform,ver,epsrel,limev,ftij,abserr,ncalls,
     &                 w,nw,idata,rdata,ier)
                  ftij=ftij/2.d0
c                  write(*,*) 'formfactor contri ',i,j,ftij/area(i),ier,
c     &              abserr,ncalls
               endif
!     
!              updating the coverage matrix
!     
               do k=1,3
                  p(k,1)=co(k,kontri(1,j))+vold(k,kontri(1,j))-pmid(k,i)
               enddo
               ddd(1)=dsqrt(p(1,1)*p(1,1)+p(2,1)*p(2,1)+p(3,1)*p(3,1))
               do k=1,3
                  p1(k)=p(k,1)/ddd(1)
               enddo
               xx(1)=p1(1)*e1(1,i)+p1(2)*e1(2,i)+p1(3)*e1(3,i)
               yy(1)=p1(1)*e2(1,i)+p1(2)*e2(2,i)+p1(3)*e2(3,i)
!     
               do k=1,3
                  p(k,2)=co(k,kontri(2,j))+vold(k,kontri(2,j))-pmid(k,i)
               enddo
               ddd(2)=dsqrt(p(1,2)*p(1,2)+p(2,2)*p(2,2)+p(3,2)*p(3,2))
               do k=1,3
                  p2(k)=p(k,2)/ddd(2)
               enddo
               xx(2)=p2(1)*e1(1,i)+p2(2)*e1(2,i)+p2(3)*e1(3,i)
               yy(2)=p2(1)*e2(1,i)+p2(2)*e2(2,i)+p2(3)*e2(3,i)
!     
               do k=1,3
                  p(k,3)=co(k,kontri(3,j))+vold(k,kontri(3,j))-pmid(k,i)
               enddo
               ddd(3)=dsqrt(p(1,3)*p(1,3)+p(2,3)*p(2,3)+p(3,3)*p(3,3))
               do k=1,3
                  p3(k)=p(k,3)/ddd(3)
               enddo
               xx(3)=p3(1)*e1(1,i)+p3(2)*e1(2,i)+p3(3)*e1(3,i)
               yy(3)=p3(1)*e2(1,i)+p3(2)*e2(2,i)+p3(3)*e2(3,i)
!     
               if(dabs(xx(2)-xx(1)).lt.1.d-5) xx(2)=xx(1)+1.d-5
               if(dabs(xx(2)-xx(1)).lt.1.d-5) xx(2)=xx(1)+1.d-5
!
!              if the surfaces are far enough away, one-point
!              integration is used
!
               if(dist(k1).gt.factor*dsqrt(area(i))*dirloc(3)) then
                  ftij=0.d0
                  do k=1,3
                     l=k-1
                     if(l.lt.1) l=3
                     xn(1)=p(2,k)*p(3,l)-p(2,l)*p(3,k)
                     xn(2)=p(3,k)*p(1,l)-p(3,l)*p(1,k)
                     xn(3)=p(1,k)*p(2,l)-p(1,l)*p(2,k)
                     xxn=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
!
!                    argument of dacos must have an absolute value
!                    smaller than or equal to 1.d0; due to
!                    round-off the value can slightly exceed one
!                    and has to be cut-off.
!
                     argument=
     &                    (p(1,k)*p(1,l)+p(2,k)*p(2,l)+p(3,k)*p(3,l))/
     &                    (ddd(k)*ddd(l))
                     if(dabs(argument).gt.1.d0) then
                        if(argument.gt.0.d0) then
                           argument=1.d0
                        else
                           argument=-1.d0
                        endif
                     endif
                     ftij=ftij+
     &                    (e3(1,i)*xn(1)
     &                    +e3(2,i)*xn(2)
     &                    +e3(3,i)*xn(3))/xxn
     &                    *dacos(argument)
c     &                    (p(1,k)*p(1,l)+p(2,k)*p(2,l)+p(3,k)*p(3,l))/
c     &                    (ddd(k)*ddd(l)))
                  enddo
                  ftij=ftij*area(i)/2.d0
c                  write(*,*) 'formfactor contri: one-point ',
c     &                  i,j,ftij/area(i)
               endif
!     
!     localizing which surface interaction the
!     triangle interaction is part of (the modulus is
!     necessary for cyclic structures)
!     
               if(j.gt.ntri) then
                  jmod=mod(j,ntri)
                  if(jmod.eq.0) jmod=ntri
               else
                  jmod=j
               endif
!     
c               call nident(iptri,i,ntr,idi)
               call nident(iptri,jmod,ntr,idj)
               f(idi,idj)=f(idi,idj)+ftij
!     
!     determining maxima and minima
!     
               xxmin=2.d0
               xxmax=-2.d0
               do k=1,3
                  if(xx(k).lt.xxmin) then
                     xxmin=xx(k)
                     imin=k
                  endif
                  if(xx(k).gt.xxmax) then
                     xxmax=xx(k)
                     imax=k
                  endif
               enddo
!     
               if(((imin.eq.1).and.(imax.eq.2)).or.
     &              ((imin.eq.2).and.(imax.eq.1))) then
                  imid=3
                  xxmid=xx(3)
               elseif(((imin.eq.2).and.(imax.eq.3)).or.
     &                 ((imin.eq.3).and.(imax.eq.2))) then
                  imid=1
                  xxmid=xx(1)
               else
                  imid=2
                  xxmid=xx(2)
               endif
!     
!     check for equal x-values
!     
               if(xxmid-xxmin.lt.1.d-5) then
                  xxmin=xxmin-1.d-5
                  xx(imin)=xxmin
               endif
               if(xxmax-xxmid.lt.1.d-5) then
                  xxmax=xxmax+1.d-5
                  xx(imax)=xxmax
               endif
!     
!     equation of the straight lines connecting the
!     triangle vertices in the local x-y plane
!     
               a(1,2)=yy(2)-yy(1)
               b(1,2)=xx(1)-xx(2)
               c(1,2)=yy(1)*xx(2)-xx(1)*yy(2)
!     
               a(2,3)=yy(3)-yy(2)
               b(2,3)=xx(2)-xx(3)
               c(2,3)=yy(2)*xx(3)-xx(2)*yy(3)
!     
               a(3,1)=yy(1)-yy(3)
               b(3,1)=xx(3)-xx(1)
               c(3,1)=yy(3)*xx(1)-xx(3)*yy(1)
!     
               a(2,1)=a(1,2)
               b(2,1)=b(1,2)
               c(2,1)=c(1,2)
               a(3,2)=a(2,3)
               b(3,2)=b(2,3)
               c(3,2)=c(2,3)
               a(1,3)=a(3,1)
               b(1,3)=b(3,1)
               c(1,3)=c(3,1)
!     
               istart=int((xxmin+1.d0+dint/2.d0)/dint)+1
               iend=int((xxmid+1.d0+dint/2.d0)/dint)
               do i1=istart,iend
                  x=dint*(i1-0.5d0)-1.d0
                  yymin=-(a(imin,imid)*x+c(imin,imid))/b(imin,imid)
                  yymax=-(a(imin,imax)*x+c(imin,imax))/b(imin,imax)
                  if(yymin.gt.yymax) then
                     dummy=yymin
                     yymin=yymax
                     yymax=dummy
                  endif
                  jstart=int((yymin+1.d0+dint/2.d0)/dint)+1
                  jend=int((yymax+1.d0+dint/2.d0)/dint)
                  do j1=jstart,jend
                     covered(i1,j1)=.true.
                  enddo
                  ncovered=ncovered+jend-jstart+1
               enddo
!     
               istart=int((xxmid+1.d0+dint/2.d0)/dint)+1
               iend=int((xxmax+1.d0+dint/2.d0)/dint)
               do i1=istart,iend
                  x=dint*(i1-0.5d0)-1.d0
                  yymin=-(a(imid,imax)*x+c(imid,imax))/b(imid,imax)
                  yymax=-(a(imin,imax)*x+c(imin,imax))/b(imin,imax)
                  if(yymin.gt.yymax) then
                     dummy=yymin
                     yymin=yymax
                     yymax=dummy
                  endif
                  jstart=int((yymin+1.d0+dint/2.d0)/dint)+1
                  jend=int((yymax+1.d0+dint/2.d0)/dint)
                  do j1=jstart,jend
                     covered(i1,j1)=.true.
                  enddo
                  ncovered=ncovered+jend-jstart+1
               enddo
               if(ncovered.eq.ng*ng)exit
!     
            enddo
         enddo
!     
!     division through total area and through pi
!     
         do i=1,ntr
            totarea=0.d0
            if(i.lt.ntr) then
               do j=iptri(i),iptri(i+1)-1
                  totarea=totarea+area(j)
               enddo
            else
               do j=iptri(i),ntri
                  totarea=totarea+area(j)
               enddo
            endif
            totarea=totarea*4.d0*datan(1.d0)
            do j=1,ntr
               f(i,j)=f(i,j)/totarea
            enddo
         enddo
!     
!     checking whether the sum of the viewfactors does not
!        exceed 1
!     
         do i=1,ntr
            fenv(i)=0.d0
            do j=1,ntr
               fenv(i)=fenv(i)+f(i,j)
            enddo
c      write(*,*) nelemload(1,i),',',sideload(i),',',fenv(i),
c     &          ',',1.d0-fenv(i)
            if((fenv(i).gt.1.d0).or.(tenv(i).lt.0)) then
               if(fenv(i).gt.0.d0) then
                  do j=1,ntr
                     f(i,j)=f(i,j)/fenv(i)
                  enddo
                  fenv(i)=1.d0
               else
                  write(*,*) '*WARNING in radmatrix: viewfactors'
                  write(*,*) '         for 3D-face''',
     &                    sideload(nloadtr(i)),''''
                  write(*,*) '         of element',
     &                    nelemload(1,nloadtr(i))
                  write(*,*) '         cannot be scaled since they are'
                  write(*,*) '         all zero'
                  write(*,*)
               endif
            endif
            fenv(i)=1.d0-fenv(i)
         enddo
!
         endif
!
         nonzero=0
         do i=1,ntr
            do j=1,ntr
               if(dabs(f(i,j)).gt.1.d-20) nonzero=nonzero+1
            enddo
         enddo
!
         if(abs(iviewfile).eq.2) then
!     
            write(*,*) 'Writing the viewfactors to file'
            write(*,*)
!     
            if(jobnamef(3)(1:1).eq.' ') then
               do i=1,132
                  if(jobnamef(1)(i:i).eq.' ') exit
               enddo
               i=i-1
               fnvw=jobnamef(1)(1:i)//'.vwf'
            else
               fnvw=jobnamef(3)
            endif
            open(10,file=fnvw,status='unknown',form='unformatted',
     &           access='sequential',err=10)
!
            write(10) nonzero
            do i=1,ntr
               do j=1,ntr
                  if(dabs(f(i,j)).gt.1.d-20) write(10) i,j,f(i,j)
               enddo
            enddo
            write(10)(fenv(i),i=1,ntr)
            close(10)
         endif
!     
      endif
!     
!        initialization of acr and bcr
!     
c      do i=1,ntr
c!     
c!        acr is (re)initialized only if the viewfactors changed
c!        or the emissivity
c!
c         if(((ithermal.eq.3).and.(iviewfile.ge.0)).or.
c     &      (iit.eq.-1).or.(iemchange.eq.1).or.
c     &      ((iit.eq.0).and.(nmethod.eq.1))) then
c            do j=1,ntr
c               acr(i,j)=0.d0
c            enddo
c         endif
c         bcr(i,1)=0.d0
c      enddo
!     
!     filling acr and bcr
!     
      do i1=1,ntr
c         if(((ithermal.eq.3).and.(iviewfile.ge.0)).or.
c     &      (iit.eq.-1).or.(iemchange.eq.1).or.
c     &      ((iit.eq.0).and.(nmethod.eq.1))) then
c            acr(i1,i1)=1.d0
c         endif
         i=nloadtr(i1)
         nelem=nelemload(1,i)
         lakonl=lakon(nelem)
!     
!     calculate the mean temperature of the face
!     
         read(sideload(i)(2:2),'(i1)') ig
!     
!     number of nodes and integration points in the face
!     
         if(lakonl(4:4).eq.'2') then
            nope=20
            nopes=8
         elseif(lakonl(4:4).eq.'8') then
            nope=8
            nopes=4
         elseif(lakonl(4:5).eq.'10') then
            nope=10
            nopes=6
         elseif(lakonl(4:4).eq.'4') then
            nope=4
            nopes=3
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         else
            nope=6
         endif
!     
         if(lakonl(4:5).eq.'8R') then
            mint2d=1
         elseif((lakonl(4:4).eq.'8').or.(lakonl(4:6).eq.'20R'))
     &           then
            if(lakonl(7:7).eq.'A') then
               mint2d=2
            else
               mint2d=4
            endif
         elseif(lakonl(4:4).eq.'2') then
            mint2d=9
         elseif(lakonl(4:5).eq.'10') then
            mint2d=3
         elseif(lakonl(4:4).eq.'4') then
            mint2d=1
         endif
!     
         if(lakonl(4:4).eq.'6') then
            mint2d=1
            if(ig.le.2) then
               nopes=3
            else
               nopes=4
            endif
         endif
         if(lakonl(4:5).eq.'15') then
            if(ig.le.2) then
               mint2d=3
               nopes=6
            else
               mint2d=4
               nopes=8
            endif
         endif
!     
!     connectivity of the element
!     
         index=ipkon(nelem)
         if(index.lt.0) then
            write(*,*) '*ERROR in radmatrix: element ',nelem
            write(*,*) '       is not defined'
            stop
         endif
         do k=1,nope
            konl(k)=kon(index+k)
         enddo
!     
!     coordinates of the nodes belonging to the face
!     
         if((nope.eq.20).or.(nope.eq.8)) then
            do k=1,nopes
               tl2(k)=vold(0,konl(ifaceq(k,ig)))
!     
               do j=1,3
                  xl2(j,k)=co(j,konl(ifaceq(k,ig)))+
     &                 vold(j,konl(ifaceq(k,ig)))
               enddo
            enddo
         elseif((nope.eq.10).or.(nope.eq.4)) then
            do k=1,nopes
               tl2(k)=vold(0,konl(ifacet(k,ig)))
               do j=1,3
                  xl2(j,k)=co(j,konl(ifacet(k,ig)))+
     &                 vold(j,konl(ifacet(k,ig)))
               enddo
            enddo
         else
            do k=1,nopes
               tl2(k)=vold(0,konl(ifacew(k,ig)))
               do j=1,3
                  xl2(j,k)=co(j,konl(ifacew(k,ig)))+
     &                 vold(j,konl(ifacew(k,ig)))
               enddo
            enddo
         endif
!     
!     integration to obtain the center of gravity and the mean
!     temperature; radiation coefficient
!     
         areamean=0.d0
         tarea(i1)=0.d0
!     
         read(sideload(i)(2:2),'(i1)') jltyp
         jltyp=jltyp+10
         if(sideload(i)(5:6).ne.'NU') then
            erad(i1)=xloadact(1,i)
!
!           if an amplitude was defined for the emissivity it is
!           assumed that the emissivity changes with the step, so
!           acr has to be calculated anew in every iteration
!
            if(nam.gt.0) then
               if(iamload(1,i).ne.0) then
                  iemchange=1
               endif
            endif
         else
            erad(i1)=0.d0
         endif
!     
         do l=1,mint2d
            if((lakonl(4:5).eq.'8R').or.
     &           ((lakonl(4:4).eq.'6').and.(nopes.eq.4))) then
               xi=gauss2d1(1,l)
               et=gauss2d1(2,l)
               weight=weight2d1(l)
            elseif((lakonl(4:4).eq.'8').or.
     &              (lakonl(4:6).eq.'20R').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.8))) then
               xi=gauss2d2(1,l)
               et=gauss2d2(2,l)
               weight=weight2d2(l)
            elseif(lakonl(4:4).eq.'2') then
               xi=gauss2d3(1,l)
               et=gauss2d3(2,l)
               weight=weight2d3(l)
            elseif((lakonl(4:5).eq.'10').or.
     &              ((lakonl(4:5).eq.'15').and.(nopes.eq.6))) then
               xi=gauss2d5(1,l)
               et=gauss2d5(2,l)
               weight=weight2d5(l)
            elseif((lakonl(4:4).eq.'4').or.
     &              ((lakonl(4:4).eq.'6').and.(nopes.eq.3))) then
               xi=gauss2d4(1,l)
               et=gauss2d4(2,l)
               weight=weight2d4(l)
            endif
!     
            if(nopes.eq.8) then
               call shape8q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.4) then
               call shape4q(xi,et,xl2,xsj2,xs2,shp2,iflag)
            elseif(nopes.eq.6) then
               call shape6tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            else
               call shape3tri(xi,et,xl2,xsj2,xs2,shp2,iflag)
            endif
!     
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
!     
            temp=0.d0
            do j=1,nopes
               temp=temp+tl2(j)*shp2(4,j)
            enddo
!     
            tarea(i1)=tarea(i1)+temp*dxsj2*weight
            areamean=areamean+dxsj2*weight
!     
            if(sideload(i)(5:6).eq.'NU') then
               areaj=dxsj2*weight
               do k=1,3
                  coords(k)=0.d0
               enddo
               do j=1,nopes
                  do k=1,3
                     coords(k)=coords(k)+xl2(k,j)*shp2(4,j)
                  enddo
               enddo
               call radiate(h(1),tenv(i1),temp,istep,
     &              iinc,tvar,nelem,l,coords,jltyp,field,nfield,
     &              sideload(i),node,areaj,vold,mi,iemchange)
               if(nmethod.eq.1) h(1)=xloadold(1,i)+
     &              (h(1)-xloadold(1,i))*reltime
               erad(i1)=erad(i1)+h(1)
            endif
!     
         enddo
         tarea(i1)=tarea(i1)/areamean-physcon(1)
         if(sideload(i)(5:6).eq.'NU') then
            erad(i1)=erad(i1)/mint2d
         endif
!     
!     radiation coefficient
!     
!     
         e=erad(i1)
         ec=1.d0-e
!     
!        acr is recalculated only if the viewfactors changed
!        or the emissivity
!
         if(((ithermal.eq.3).and.(iviewfile.ge.0)).or.
     &      (iit.eq.-1).or.(iemchange.eq.1).or.
     &      ((iit.eq.0).and.(nmethod.eq.1))) then
            do j1=1,ntr
c               acr(i1,j1)=acr(i1,j1)-ec*f(i1,j1)
               acr(i1,j1)=-ec*f(i1,j1)
            enddo
            acr(i1,i1)=1.d0+acr(i1,i1)
         endif
         bcr(i1,1)=physcon(2)*(e*tarea(i1)**4+
     &        ec*fenv(i1)*tenv(i1)**4)
!     
      enddo
!     
      return
!
 10   write(*,*) '*ERROR in radmatrix: could not open file ',fnvw
      stop
      end
!     
!     function to be integrated
!     
      real*8 function fform(x,y,idata,rdata)
!
      implicit none
!     
      integer k,l,number,idata(1)
!     
      real*8 pint(3),ddd(3),xn(3),vj(3,3),
     &   unitvec(3,3),p(3,3),xxn,x,y,porigin(3),rdata(1)
!     
      common /formfactor/ vj,unitvec,porigin
!
      data number /0/
      save number
!
      number=number+1
      do k=1,3
         pint(k)=porigin(k)+x*unitvec(k,1)+y*unitvec(k,2)
      enddo
!
      do k=1,3
         p(k,1)=vj(k,1)-pint(k)
      enddo
      ddd(1)=dsqrt(p(1,1)*p(1,1)+p(2,1)*p(2,1)+p(3,1)*p(3,1))
!
      do k=1,3
         p(k,2)=vj(k,2)-pint(k)
      enddo
      ddd(2)=dsqrt(p(1,2)*p(1,2)+p(2,2)*p(2,2)+p(3,2)*p(3,2))
!
      do k=1,3
         p(k,3)=vj(k,3)-pint(k)
      enddo
      ddd(3)=dsqrt(p(1,3)*p(1,3)+p(2,3)*p(2,3)+p(3,3)*p(3,3))
!     
!     calculating the contribution
!     
      fform=0.d0
      do k=1,3
         l=k-1
         if(l.lt.1) l=3
         xn(1)=p(2,k)*p(3,l)-p(2,l)*p(3,k)
         xn(2)=p(3,k)*p(1,l)-p(3,l)*p(1,k)
         xn(3)=p(1,k)*p(2,l)-p(1,l)*p(2,k)
         xxn=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
         fform=fform+
     &        (unitvec(1,3)*xn(1)
     &        +unitvec(2,3)*xn(2)
     &        +unitvec(3,3)*xn(3))/xxn
     &        *dacos(
     &        (p(1,k)*p(1,l)+p(2,k)*p(2,l)+p(3,k)*p(3,l))/
     &        (ddd(k)*ddd(l)))
      enddo
!
      return
      end
      
