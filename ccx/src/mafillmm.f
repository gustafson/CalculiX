!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2024 Guido Dhondt
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
      subroutine mafillmm(co,nodedesiinv,iregion,au,ad,
     &     aub,adb,irow,jq,ipkonfa,konfa,lakonfa,nodedesipos,
     &     idesiface,nsurfa,nsurfb,area)
!     
!     calculation of the entries of the mass matrix and it's derivative
!     
      implicit none
!     
      character*8 lakonfa(*)
!     
      integer jj,ll,kk,l,jdof1,jdof2,n,m,
     &     nopes,ifour,ithree,indexel,mint2d,iflag,indexs,
     &     i,node1,node2,nodedesiinv(*),iregion,irow(*),
     &     jq(*),nsurfa,nsurfb,isurf,ipkonfa(*),konfa(*),
     &     konfal(8),nodedesipos(*),id,ipointer,idesiface(*)
!     
      real*8 xi,et,xl(3,9),xs(3,2),xsj(3),shp(7,9),co(3,*),au(*),ad(*),
     &     aub(*),adb(*),xsjj,weight,val,area(*),xn(3),
     &     shpprj1(3),shpprj2(3),sclprd1,sclprd2
!     
!     flag for shape functions
!     
      data iflag /3/
      data indexel /0/
      save indexel
      include "gauss.f"
!     
!     flag for shape functions
!     
      ifour=4
      ithree=3
!     
!---------------------------------------------------------------------
!     Calculating the entries of the matrices
!---------------------------------------------------------------------
!     
!     Loop over all external surfaces
!     
      do i=nsurfa,nsurfb
        isurf=idesiface(i)
        indexs=ipkonfa(isurf)
!     
!     mint2d: # of integration points on the surface
!     nopes:  # of nodes in the surface
!     nope:   # of nodes in the element   
!     
        if(lakonfa(isurf)(1:2).eq.'S3') then
          mint2d=3
          nopes=3
        elseif(lakonfa(isurf)(1:2).eq.'S4') then
          mint2d=4
          nopes=4
        elseif(lakonfa(isurf)(1:2).eq.'S6') then
          mint2d=7
          nopes=6
        elseif(lakonfa(isurf)(1:2).eq.'S8') then
          mint2d=9
          nopes=8
        else
          exit
        endif
!     
!     storing node numbers and coordinates of nodes of the surface   
!     
        do l=1,nopes
          konfal(l)=konfa(indexs+l)
          do n=1,3
            xl(n,l)=co(n,konfal(l))
          enddo
        enddo 
!     
!     Evaluate Shape functions and their derivatives 
!     
        do kk=1,mint2d   
          if(lakonfa(isurf)(1:2).eq.'S3') then
            xi=gauss2d5(1,kk)
            et=gauss2d5(2,kk)
            weight=weight2d5(kk)
            call shape3tri(xi,et,xl,xsj,xs,shp,iflag)
          elseif(lakonfa(isurf)(1:2).eq.'S4') then
            xi=gauss2d2(1,kk)
            et=gauss2d2(2,kk)
            weight=weight2d2(kk)
            call shape4q(xi,et,xl,xsj,xs,shp,iflag)
          elseif(lakonfa(isurf)(1:2).eq.'S6') then
            xi=gauss2d6(1,kk)
            et=gauss2d6(2,kk)
            weight=weight2d6(kk)
            call shape6tri(xi,et,xl,xsj,xs,shp,iflag)
          elseif(lakonfa(isurf)(1:2).eq.'S8') then
            xi=gauss2d3(1,kk)
            et=gauss2d3(2,kk)
            weight=weight2d3(kk)
            call shape8q(xi,et,xl,xsj,xs,shp,iflag)       
          endif
!     
!     Loop over all nodes on the surface element
!     
          do jj=1,nopes
            node1=konfal(jj)
            jdof1=nodedesipos(node1)
!     
!     Loop over all nodes on the surface element
!     
            do ll=1,nopes
              node2=konfal(ll)
              jdof2=nodedesipos(node2)
!     
!     check if both nodes are designvariables
!     
              if((nodedesiinv(node1).eq.1).and.
     &             (nodedesiinv(node2).eq.1)) then      
!     
!     Calculate Jacobian determinant     
!     
                xsjj=dsqrt(xsj(1)**2+xsj(2)**2+xsj(3)**2)
                do m=1,3
                  xn(m)=xsj(m)/xsjj
                enddo
!     
!     entry on the main diagonal
!     
                if(node1.eq.node2) then
                  adb(jdof1)=adb(jdof1)
     &                 +shp(4,jj)**2*weight*xsjj
!     
                  area(jdof1)=area(jdof1)
     &                 +shp(4,jj)*weight*xsjj    
!     
                  sclprd1=shp(1,jj)*xn(1)
     &                 +shp(2,jj)*xn(2)
     &                 +shp(3,jj)*xn(3)
!
!                 projection of the gradient of the shape function
!                 on the face
!
                  shpprj1(1)=shp(1,jj)-xn(1)*sclprd1
                  shpprj1(2)=shp(2,jj)-xn(2)*sclprd1
                  shpprj1(3)=shp(3,jj)-xn(3)*sclprd1
!     
                  val=(shpprj1(1)**2
     &                 +shpprj1(2)**2
     &                 +shpprj1(3)**2)*weight*xsjj
                  ad(jdof1)=ad(jdof1)+val
!     
!     entry on sub diagonal
!     
                elseif(jdof2.gt.jdof1) then
!     
                  call nident(irow(jq(jdof1)),jdof2,
     &                 jq(jdof1+1)-jq(jdof1),id)
!     
                  ipointer=jq(jdof1)+id-1
!     
                  if(irow(ipointer).ne.jdof2) then
                    write(*,*) '*ERROR in mafillmm: 
     &coefficient should be 0'
                    call exit(201)
                  else
                    aub(ipointer)=aub(ipointer)
     &                   +shp(4,jj)*shp(4,ll)*weight*xsjj
!     
                    sclprd1=shp(1,jj)*xn(1)
     &                   +shp(2,jj)*xn(2)
     &                   +shp(3,jj)*xn(3)
                    shpprj1(1)=shp(1,jj)-xn(1)*sclprd1
                    shpprj1(2)=shp(2,jj)-xn(2)*sclprd1
                    shpprj1(3)=shp(3,jj)-xn(3)*sclprd1
!     
                    sclprd2=shp(1,ll)*xn(1)
     &                   +shp(2,ll)*xn(2)
     &                   +shp(3,ll)*xn(3)
                    shpprj2(1)=shp(1,ll)-xn(1)*sclprd2
                    shpprj2(2)=shp(2,ll)-xn(2)*sclprd2
                    shpprj2(3)=shp(3,ll)-xn(3)*sclprd2
!     
                    val=(shpprj1(1)*shpprj2(1)
     &                   +shpprj1(2)*shpprj2(2)
     &                   +shpprj1(3)*shpprj2(3))*weight*xsjj
                    au(ipointer)=au(ipointer)+val
!     
                  endif            
                endif
              endif
            enddo
          enddo
        enddo
      enddo   
!     
      return
      end
