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
c>
c>    cuts a triangle of the master surface with a slave surface
c>    inserts new active edges into iactiveline for current triangle
c> @see checktriaedge
c> @param [in,out]   inodesin,nnodesin !
c> @param [in,out]   inodesout, nnodesout ! 
c> @param [in]       nopes
c> @param [in]       slavstraight
c> @param [in]       xn
c> @param [in]       co
c> @param [in]       xl2s
c> @param [in]       ipe
c> @param [in]       ime
c> @param [in,out]   iactiveline !
c> @param [in,out]   nactiveline !
c> @param [in,out]   intersec !
c> @param [in,out]   xntersec !
c> @param [in,out]   ifreeintersec !
c> @param [in]       itri
c> @param [in]       koncont
c> @param [in]       itriacornerl
c> @param [in,out]   nintpoint
c> @param [out]      pslavsurf
c> @param [in]       ncont
c> @param [out]      imastsurf
c> @param [out]      pmastsurf
c> @param [in]       xl2m
c> @param [in]       nnodelem
c> @param [in]       vold
c> @param [in]       mi
c> @param [in]       pnodesin
c> @param [in]       straight
c> @param [out]      gapmints
c> @param [in]       ssurf
!
      subroutine treattriangleS(
     &  nopes,slavstraight,xn,xns,co,xl2s,xl2sp,
     &  ipe,ime,iactiveline,nactiveline,
     &  ifreeintersec,nelemm,itriacornerl,nintpoint,pslavsurf,
     &  ncont,imastsurf,pmastsurf,xl2m,nnodelem,nodem,mi,pnodesin,
     &  straight,gapmints,issurf,areaslav,debug)
!     
!    Author: Saskia Sitzmann     
!     
      implicit none
!
      logical debug
!
      integer nvertex,
     &  nnodesout,nopes,ipe(*),ime(4,*),iactiveline(3,*),nactiveline,
     &  ifreeintersec,itriacornerl(4),
     &  i,j,k,nintpoint,ncont,imastsurf(*),nnodelem,mi(*),ijk,
     &  nodem(*),modf,nelemm,itri,k_max,issurf,
     &  ii,jj,nipold,il
!
      real*8 pvertex(3,13),pnodesin(3,*),slavstraight(20),xn(3),
     &  co(3,*),xilm,etlm,straight(16,*),xnl(3),
     &  xl2s(3,*),p(3,7),p1(3),p2(3),pslavsurf(3,*),
     &  ratio(8),dist,xil,etl,area,areax,areay,areaz,pmastsurf(2,*),
     &  xl2m(3,8),al,gapmints(*),err,xns(3,8),
     &  xl2sp(3,8),xl2mp(3,8),cgp(3),ph2(3),spm,
     &  xm(3),pm(3),xs2m(3,2),xsj2m(3),
     &  shp2m(7,8),ph(3),xilm2,
     &  shp2s(7,8),ps(3),xit(3),etat(3),
     &  areaslav
!     
      data ijk /0/
      save ijk
!     
      include "gauss.f"
!     
      err=1.d-6
      nvertex=0
      nipold=nintpoint
!     
      if(debug) write(20,*) 'TT:slavsurf',issurf,' melem',nelemm
!     
!     Project master nodes to meanplane, needed for Sutherland-Hodgman
!     
      do j=1, nnodelem
         al=-xn(1)*xl2m(1,j)-xn(2)*
     &        xl2m(2,j)-xn(3)*
     &        xl2m(3,j)-slavstraight(nopes*4+4)
         do k=1,3
            xl2mp(k,j)= xl2m(k,j)+al*xn(k)    
         enddo
      enddo 
!     
      if(debug) then
         write(20,*) 'TT: xl2sp'
         do j=1,nopes
            write(20,*)(xl2sp(k,j),k=1,3)
         enddo       
         write(20,*)'TT:xl2mp'
         do j=1,nnodelem
            write(20,*)(xl2mp(k,j),k=1,3)
         enddo  
      endif
      
!     
 111  format(3(e27.20))
!     
!     call Sutherland-Hodgman Algo
!     
      call sutherland_hodgman(nopes,xn,xl2sp,xl2mp,nodem,
     &     ipe,ime,iactiveline,nactiveline,
     &     ifreeintersec,nelemm,nnodelem,
     &     nvertex,pvertex) 
!     
!     do we have a degenerated triangle?
!     
      if(debug)then
         write(20,*) 'nelemm',nelemm ,'p_new'
         do k=1,nvertex
            write(20,111)(pvertex(i,k),i=1,3)
         enddo
      endif
!     
      do k=1,3
         cgp(k)=0.0
      enddo
      if(nvertex.lt.3) return       
!     
      if(nvertex==3)then
         do k=1,3
            cgp(k)=pvertex(k,nvertex)
         enddo
         nvertex=nvertex-1
         k_max=0
      else
         do i=1,nvertex
            do k=1,3
               cgp(k)=cgp(k)+pvertex(k,i)/nvertex
            enddo
         enddo
         k_max=nvertex-1
      endif 
!     
      if(debug)then
         write(20,*) 'TT: nactiveline',nactiveline,'iactiveline'
         write(20,*) (iactiveline(1,k),k=1, nactiveline)
         write(20,*) 'cg' ,(cgp(k),k=1,3)
         write(20,*)'********************************************' 
      endif    
!     
!     Project center point back on slave face
!     
      call attachlineS(xl2s,cgp,nopes,ratio,dist,xit(3),etat(3),xn)
!     
!     generating integration points on the slave surface S
!     
      do k=0,k_max
!     
!     storing the triangulation of the slave surfaces
!     
c         write(40,*) '# islavsurf',issurf,' melem',nelemm
         ijk=ijk+1
         write(40,100) ijk,(cgp(i),i=1,3)
         ijk=ijk+1
         write(40,100) ,ijk,(pvertex(i,modf(nvertex,k)),i=1,3)
         ijk=ijk+1
         write(40,100) ijk,(pvertex(i,modf(nvertex,k+1)),i=1,3)
         write(40,101) ijk-2,ijk-2,ijk-1
         write(40,101) ijk-1,ijk-1,ijk
         write(40,101) ijk,ijk,ijk-2
 100     format('PNT ',i10,'P',3(1x,e21.14))
 101     format('LINE ',i10,'L',i10,'P ',i10,'P')
!     
!     Project back on slave surface
!     
         call attachlineS(xl2s,pvertex(1:3,modf(nvertex,k)),
     &        nopes,ratio,dist,xit(1),etat(1),xn)
         call attachlineS(xl2s,pvertex(1:3,modf(nvertex,k+1)),
     &        nopes,ratio,dist,xit(2),etat(2),xn)
         p1(1)=xit(1)-xit(3)
         p1(2)= etat(1)-etat(3)
         p1(3)=0
         p2(1)=xit(2)-xit(3)
         p2(2)=etat(2)-etat(3)
         p2(3)=0
         areax=((p1(2)*p2(3))-(p2(2)*p1(3)))**2
         areay=(-(p1(1)*p2(3))+(p2(1)*p1(3)))**2
         areaz=((p1(1)*p2(2))-(p2(1)*p1(2)))**2
         area=dsqrt(areax+areay+areaz)/2.
c         if(area.lt.1.e-4) cycle
         if(nopes.eq.4.and.areaslav+area-4.0.gt.1.e-3
     &        .and.nactiveline.gt.0)then
           write(*,*)'TT: face',issurf,'loop in slavintmortar'
           write(*,*)'area',areaslav,'+',area,'.gt.4!',nactiveline
           nactiveline=0
           return
         endif
         if(nopes.eq.3.and.areaslav+area-0.5.gt.1.e-3
     &       .and.nactiveline.gt.0)then
           write(*,*)'TT: face',issurf,'loop in slavintmortar'
           write(*,*)'area',areaslav,'+',area,'.gt.0.5!'
           nactiveline=0
           return
         endif
         areaslav=areaslav+area
         if(debug)then
            write(20,106) k, area, areaslav
 106        format('tri',i10,' area',e15.8,' atot',e15.8)
            write(20,*) '# für itri',k, 'werden 7intp gen.'
         endif
!     
!     7 points scheme
!     
         do i=1,7
            xil= xit(3)*gauss2d6(1,i)+
     &           xit(1)*gauss2d6(2,i)+
     &           xit(2)*(1-gauss2d6(1,i)-gauss2d6(2,i))
            
            etl= etat(3)*gauss2d6(1,i)+
     &           etat(1)*gauss2d6(2,i)+
     &           etat(2)*(1-gauss2d6(1,i)-gauss2d6(2,i))
!     
            call evalshapefunc(xil,etl,xns,nopes,xnl)
            call evalshapefunc(xil,etl,xl2s,nopes,ps)
!     
            nintpoint=nintpoint+1
!     
!     projection of the integration point in the mean
!     slave plane onto the slave surface
!     
!     projection of the master integration point onto the
!     master surface in order to get the local coordinates
!     own xn for every integration point?
!     
c            call attachlineS(xl2m,ps,nnodelem,ratio,dist,xilm,etlm,xnl)
            call attachlineS(xl2m,ps,nnodelem,ratio,dist,xilm,etlm,xn)
            call evalshapefunc(xilm,etlm,xl2m,nnodelem,pm)   
!
            if(debug)then     
            ijk=ijk+1
            write(40,100) ijk,(ps(jj),jj=1,3)
            ijk=ijk+1
            write(40,100) ijk,(pm(jj),jj=1,3)
            write(40,101) ijk-1,ijk-1,ijk 
            endif
!     
!     Calculation of the gap function at the integration point
!     
            do ii=1,3
               ph(ii)=pm(ii)-ps(ii)
            enddo

            al= sqrt(ph(1)**2+ph(2)**2+ph(3)**2)
c            spm= ph(1)*xnl(1)
c     &           +ph(2)*xnl(2)
c     &           +ph(3)*xnl(3)
            spm= ph(1)*xn(1)
     &           +ph(2)*xn(2)
     &           +ph(3)*xn(3)
c            do ii=1,3
c               ph2(ii)=ps(ii)+spm*xnl(ii)
c            enddo
            do ii=1,3
               ph2(ii)=ps(ii)+spm*xn(ii)
            enddo     
c            ijk=ijk+1
c            write(40,100) ijk,(ph2(jj),jj=1,3)
c            write(40,101) ijk-1,ijk-2,ijk     
            gapmints(nintpoint)=spm
            pslavsurf(1,nintpoint)=xil
            pslavsurf(2,nintpoint)=etl
            pslavsurf(3,nintpoint)=area*weight2d6(i)
            pmastsurf(1,nintpoint)=xilm
            pmastsurf(2,nintpoint)=etlm
            imastsurf(nintpoint)=nelemm
            if(debug)then
               write(30,201) xil,etl,pslavsurf(3,nintpoint),spm,nelemm 
            endif 
 201        format('xil ',1x,e15.8,2x,'etl',2x,e15.8,2x,'w ',e15.8,2x,
     &           e15.8,2x,'M',i10)
         enddo
!     
      enddo
c      write(20,*)'********************'
      return
      end
!*************************************************************************
      subroutine evalshapefunc(xil,etl,xl2,nopes,p)
!     
      implicit none
!     
      integer i,j,k,nopes,iflag
!     
      real*8  xl2(3,8),maststraight(20),xm(3),dd,xs2(3,2),
     &     shp2(7,8),p(3),xsj2s(3),spm,xsj2(3),xil,etl
!     
      iflag=1
      do j=1,3
         p(j)=0.0
      enddo 
      if(nopes.eq.8)then
         call shape8q(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.4)then
         call shape4q(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.6)then
         call shape6tri(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xil,etl,xl2,xsj2,xs2,shp2,iflag)
      endif 
      do i=1,nopes
         do j=1,3
            p(j)=p(j)+xl2(j,i)*shp2(4,i)
         enddo
      enddo 
!     
      return
      end
