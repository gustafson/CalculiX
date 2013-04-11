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
      subroutine updatecontpen(koncont,ncont,co,vold,cg,straight,mi,
     &	             	    imastnode,nmastnode,xmastnor,ntie,tieset,
     &	             	    nset,set,istartset,
     &	             	    iendset,ialset,ipkon,lakon,kon)
!
!     update geometric date of the contact master surface triangulation
!
      implicit none
!
      character*1 kind1,kind2
      character*8 lakon(*)
      character*81 tieset(3,*),mastset,set(*)
!
      integer koncont(4,*),ncont,i,j,k,node,mi(*),imastnode(*),
     &	      nmastnode(*),ntie,imast,
     &	      istartset(*),iendset(*),ialset(*),ifacem,nelemm,
     &	      jfacem,indexe,ipkon(*),nopem,nope,konl(20),kon(*),m,
     &	      ifaceq(8,6),ifacet(6,4),ifacew1(4,5),xquad(2,8),
     &	      ifacew2(8,5),id,index1,indexnode(8),l,iflag,nset,itriact,
     &	      ipos,ntrifaces
!
      real*8 co(3,*),vold(0:mi(2),*),cg(3,*),straight(16,*),col(3,3),
     &	     xmastnor(3,*),xl2m(3,8),xi,et,dd,xsj2(3),shp2(7,8),
     &	     xs2(3,2),xnor(3,3)
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             4,6,3,1/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
!
      data iflag /2/
!     
!     new added data for the local coodinates for nodes
!
      data xquad /-1, -1,
     &           1, -1,
     &           1, 1,
     &           -1, 1,
     &           0, -1,
     &          1, 0,
     &           0, 1,
     &           -1, 0/
!
      kind1='C'
      kind2='-'
!
!     calculation of the normals on the master surface 
!     at the vertex nodes of the master surface -> xmastnor
!     Iteration over all tiesets
!     
      do i=1,ntie
!
!        check for contact conditions
!
         if((tieset(1,i)(81:81).eq.'C').or.
     &      (tieset(1,i)(81:81).eq.'-')) then

	    mastset=tieset(3,i)
            do j=1,nset
               if(set(j).eq.mastset) exit
            enddo
            if(j.gt.nset) then
               write(*,*) '*ERROR in tiefaccont: master surface'
               write(*,*) '       does not exist'
               stop
            endif
            imast=j
!
            do j=istartset(imast),iendset(imast)
               if(ialset(j).gt.0) then
!               
!           Decide imastnode, and nmastnode
!
                  ifacem = ialset(j)
                  nelemm = int(ifacem/10)
                  jfacem = ifacem - nelemm*10
                  indexe = ipkon(nelemm)
!     
!     Decide on the max integration point number, just consider 2D situation 
!     
                  if(lakon(nelemm)(4:5).eq.'8R') then
                     nopem=4
                     nope=8
                  elseif(lakon(nelemm)(4:4).eq.'8') then
                     nopem=4
                     nope=8
                  elseif(lakon(nelemm)(4:6).eq.'20R') then
                     nopem=8
                     nope=20
                  elseif(lakon(nelemm)(4:4).eq.'2') then
                     nopem=8
                     nope=20
                  elseif(lakon(nelemm)(4:5).eq.'10') then
                     nopem=6
                     nope=10
                  elseif(lakon(nelemm)(4:4).eq.'4') then
                     nopem=3
                     nope=4
!     
!     treatment of wedge faces
!     
                  elseif(lakon(nelemm)(4:4).eq.'6') then
                     nope=6
                     if(jfacem.le.2) then
                     	nopem=3
                     else
                     	nopem=4
                     endif
                  elseif(lakon(nelemm)(4:5).eq.'15') then
                     nope=15
                     if(jfacem.le.2) then
                     	nopem=6
                     else
                     	nopem=8
                     endif
                  endif
!     
!     actual position of the nodes belonging to the
!     master surface
!     
                  do k=1,nope
                     konl(k)=kon(ipkon(nelemm)+k)
                  enddo
!     
                  if((nope.eq.20).or.(nope.eq.8)) then
                     do m=1,nopem
                     	do k=1,3
                     	   xl2m(k,m)=co(k,konl(ifaceq(m,jfacem)))+
     &                        vold(k,konl(ifaceq(m,jfacem)))
                     	enddo
                     enddo
                  elseif((nope.eq.10).or.(nope.eq.4)) then
                     do m=1,nopem
                     	do k=1,3
                     	   xl2m(k,m)=co(k,konl(ifacet(m,jfacem)))+
     &                        vold(k,konl(ifacet(m,jfacem)))
                     	enddo
                     enddo
                  else
                     do m=1,nopem
                     	do k=1,3
                     	   xl2m(k,m)=co(k,konl(ifacew1(m,jfacem)))+
     &                        vold(k,konl(ifacew1(m,jfacem)))
                     	enddo
                     enddo
                  endif
            
!     calculate the normal vector in the nodes belonging to the master surface
!     
                  if(nopem.eq.8) then
                     do m = 1, nopem
                     	xi = xquad(1,m)
                     	et = xquad(2,m)
                     	call shape8q(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                     	dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)
     &                     + xsj2(3)*xsj2(3))
                     	xsj2(1) = xsj2(1)/dd
                     	xsj2(2) = xsj2(2)/dd
                     	xsj2(3) = xsj2(3)/dd
!     
                     	if(nope.eq.20) then
                     	   node = konl(ifaceq(m,jfacem))
                     	elseif(nope.eq.15) then
                     	   node=konl(ifacew2(m,jfacem))
                     	endif
!     
                     	call nident(imastnode(nmastnode(i)+1), node, 
     &                     nmastnode(i+1)-nmastnode(i), id)
                     	index1=nmastnode(i)+id
                     	indexnode(m)=index1
                     	xmastnor(1,index1) = xmastnor(1,index1)
     &                     +xsj2(1)
                     	xmastnor(2,index1) = xmastnor(2,index1)
     &                     +xsj2(2)
                     	xmastnor(3,index1) = xmastnor(3,index1)
     &                     +xsj2(3)
                     enddo
                  elseif(nopem.eq.4) then
                     do m = 1, nopem
                     	xi = xquad(1,m)
                     	et = xquad(2,m)
                     	call shape4q(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                     	dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &                     + xsj2(3)*xsj2(3))
                     	xsj2(1) = xsj2(1)/dd
                     	xsj2(2) = xsj2(2)/dd
                     	xsj2(3) = xsj2(3)/dd
!     
                     	if(nope.eq.8) then
                     	   node = konl(ifaceq(m,jfacem))
                     	elseif(nope.eq.6) then
                     	   node=konl(ifacew1(m,jfacem))
                     	endif
!     
                     	call nident(imastnode(nmastnode(i)+1), node, 
     &                     nmastnode(i+1)-nmastnode(i), id)
!     
                     	index1=nmastnode(i)+id
                     	indexnode(m)=index1
                     	xmastnor(1,index1) = xmastnor(1,index1)
     &                     +xsj2(1)
                     	xmastnor(2,index1) = xmastnor(2,index1)
     &                     +xsj2(2)
                     	xmastnor(3,index1) = xmastnor(3,index1)
     &                     +xsj2(3)
                     enddo
                  elseif(nopem.eq.6) then
                     do m = 1, nopem
                     	xi = xquad(1,m)
                     	et = xquad(2,m)
                     	call shape6tri(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                     	dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &                 + xsj2(3)*xsj2(3))
                     	xsj2(1) = xsj2(1)/dd
                     	xsj2(2) = xsj2(2)/dd
                     	xsj2(3) = xsj2(3)/dd
!     
                     	if(nope.eq.10) then
                     	   node = konl(ifacet(m,jfacem))
                     	elseif(nope.eq.15) then
                     	   node = konl(ifacew2(m,jfacem))
                     	endif
!     
                     	call nident(imastnode(nmastnode(i)+1), node, 
     &                     nmastnode(i+1)-nmastnode(i), id)
                     	index1=nmastnode(i)+id
                     	indexnode(m)=index1
                     	xmastnor(1,index1) = xmastnor(1,index1)
     &                     +xsj2(1)
                     	xmastnor(2,index1) = xmastnor(2,index1)
     &                     +xsj2(2)
                     	xmastnor(3,index1) = xmastnor(3,index1)
     &                     +xsj2(3)
                     enddo
                  else
                     do m = 1, nopem
                     	xi = xquad(1,m)
                     	et = xquad(2,m)
                     	call shape3tri(xi,et,xl2m,xsj2,xs2,shp2,iflag)
                     	dd = dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
     &                     + xsj2(3)*xsj2(3))
                     	xsj2(1) = xsj2(1)/dd
                     	xsj2(2) = xsj2(2)/dd
                     	xsj2(3) = xsj2(3)/dd
!     
                     	if(nope.eq.6) then
                     	   node = konl(ifacew1(m,jfacem))
                     	elseif(nope.eq.4) then
                     	   node = konl(ifacet(m,jfacem))
                     	endif
!     
                     	call nident(imastnode(nmastnode(i)+1), node, 
     &                     nmastnode(i+1)-nmastnode(i), id)
                     	index1=nmastnode(i)+id
                     	indexnode(m)=index1
                     	xmastnor(1,nmastnode(i)+id) = xmastnor(1,index1)
     &                     +xsj2(1)
                     	xmastnor(2,nmastnode(i)+id) = xmastnor(2,index1)
     &                     +xsj2(2)
                     	xmastnor(3,nmastnode(i)+id) = xmastnor(3,index1)
     &                     +xsj2(3)
                     enddo
                  endif
               endif
            enddo
!     
!     normalizing the normals
!     
            do l=nmastnode(i)+1,nmastnode(i+1)
               dd=dsqrt(xmastnor(1,l)**2+xmastnor(2,l)**2+
     &            xmastnor(3,l)**2)
               do m=1,3
                  xmastnor(m,l)=xmastnor(m,l)/dd
               enddo
	    enddo	    
         endif 
      enddo
!
!     Loop over all tiesets
!
      itriact=0
c      do i=1,ntie
c         write(*,*) 'tie ',i,nmastnode(i),nmastnode(i+1)
c      enddo
      do i=1,ntie
c         do l=nmastnode(i)+1,nmastnode(i+1)
c            write(*,*) 'master ',imastnode(l),(xmastnor(m,l),m=1,3)
c         enddo
!  
!        check for contact conditions
!
         if((tieset(1,i)(81:81).eq.kind1).or.
     &      (tieset(1,i)(81:81).eq.kind2)) then
            mastset=tieset(3,i)
!
!           determining the master surface
!
            do j=1,nset
               if(set(j).eq.mastset) exit
            enddo
            if(j.gt.nset) then
               ipos=index(mastset,' ')
               write(*,*) '*ERROR in updatecont: master surface',
     &               mastset(1:ipos-2)
               write(*,*) '       does not exist or does not contain'
               write(*,*) '       element faces'
               stop
            endif
            imast=j
	    
      	    do j=istartset(imast),iendset(imast)
                  nelemm=int(ialset(j)/10.d0)
                  jfacem=ialset(j)-10*nelemm
!
                  if(lakon(nelemm)(4:4).eq.'2') then
                     ntrifaces=6
                  elseif(lakon(nelemm)(4:4).eq.'8') then
                     ntrifaces=2
                  elseif(lakon(nelemm)(4:5).eq.'10') then
                     ntrifaces=4
                  elseif(lakon(nelemm)(4:4).eq.'4') then
                     ntrifaces=1
                  elseif(lakon(nelemm)(4:5).eq.'15') then
                     if(jfacem.le.2) then
                        ntrifaces=4
                     else
                        ntrifaces=6
                     endif
                  elseif(lakon(nelemm)(4:4).eq.'6') then
                     if(jfacem.le.2) then
                        ntrifaces=1
                     else
                        ntrifaces=2
                     endif
                  endif

!
!     	 loop over the master triangles
!
		  do k=1,ntrifaces
		     itriact=itriact+1
		     do l=1,3
		     	node=koncont(l,itriact)
c                        write(*,*) 'node ',node
		     	do m=1,3
			   col(m,l)=co(m,node)+vold(m,node)
			enddo
			call nident(imastnode(nmastnode(i)+1),node,
     &	             	      	    nmastnode(i+1)-nmastnode(i),id)
      	             	index1=nmastnode(i)+id
			do m=1,3
			   xnor(m,l)=xmastnor(m,index1)		   
			enddo
c			if(itriact.eq.69)write(*,*)'xnor=',(xnor(m,l),m=1,3)		    		     
		     enddo
!
!        center of gravity of the triangles
!
      	             do l=1,3
		     	cg(l,itriact)=col(l,1)
		     enddo
		     do m=2,3
		     	do l=1,3
			   cg(l,itriact)=cg(l,itriact)+col(l,m)
			enddo
		     enddo
		     do l=1,3
		     	cg(l,itriact)=cg(l,itriact)/3.d0
		     enddo
!
!        calculating the equation of the triangle plane and the planes
!        through the edges of a triangle and perpendicular to the
!     	 representative normal of each triangle edge
!
                     call straighteq3dpen(col,straight(1,itriact),xnor)
c                     if(itriact.eq.69) then
c             write(*,*) 'updatecont 69 ',(koncont(l,itriact),l=1,3)
c             write(*,*) 'updatecont xnor',(xnor(l,1),l=1,3)
c             write(*,*) 'updatecont xnor',(xnor(l,2),l=1,3)
c             write(*,*) 'updatecont xnor',(xnor(l,3),l=1,3)
c             write(*,*) 'updatecont strai',(straight(l,69),l=1,4)
c             write(*,*) 'updatecont strai',(straight(l,69),l=5,8)
c             write(*,*) 'updatecont strai',(straight(l,69),l=9,12)
c             write(*,*) 'updatecont strai',(straight(l,69),l=13,16)
c                     endif
		  enddo
	    enddo
      	 endif	    
      enddo
!
      return
      end
