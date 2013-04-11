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
      subroutine flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,rho,physcon,g,co,dvi,numf,
     &     vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_)
!
!     determine whether the flux in the element is an unknown 
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!      
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(5),idirf(5),kflag,ipkon(*),kon(*),
     &     nshcon(*), nrhcon(*),ntmat_
!      
      real*8 prop(*),v(0:4,*),xflow,f,df(5),R,cp,physcon(*),rho,
     &     g(3),co(3,*),dvi,vold(0:4,*),shcon(0:3,ntmat_,*),
     &     rhcon(0:1,ntmat_,*)
!
      if(lakon(nelem)(2:6).eq.'CARBS') then  
!
         call carbon_seal(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,R,physcon,dvi,numf,set)
!     
      elseif(lakon(nelem)(2:5).eq.'CHAR') then 
!     
         call characteristic(node1,node2,nodem,nelem,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,physcon,numf,set)
!
      elseif(lakon(nelem)(2:4).eq.'LAB') then 
!         
         call labyrinth(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,co,dvi,numf,vold,set,kon,ipkon)
! 
      elseif(lakon(nelem)(2:5).eq.'GAPI') then 
!         
         call gaspipe(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_)
!
      elseif(lakon(nelem)(2:5).eq.'GAPX') then 

         call gaspipe_fanno(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold)
!
      elseif(lakon(nelem)(2:5).eq.'LIPI') then
!         
         call liquidpipe(node1,node2,nodem,nelem,lakon,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,dvi,numf,vold)
!
      elseif(lakon(nelem)(2:5).eq.'LIPU') then
!         
         call liquidpump(node1,node2,nodem,nelem,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,numf)
!      
      elseif(lakon(nelem)(2:4).eq.'MRG') then 
!     
         call moehring(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,dvi,numf,set)
!
      elseif(lakon(nelem)(2:3).eq.'OR') then 
!         
         call orifice(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,co,vold)
!
      elseif(lakon(nelem)(2:3).eq.'RE') then 
!         
         call restrictor(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_)
!
      elseif(lakon(nelem)(2:5).eq.'RIMS') then   
!
         call rimseal(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set)
!
      elseif(lakon(nelem)(2:6).eq.'SPUMP') then 
! 
        call scavenge_pump(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_)   
!
      elseif(lakon(nelem)(2:3).eq.'VO') then 
!     
         call vortex(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,numf,set) 
!     
      elseif((lakon(nelem)(2:4).eq.'ATR')
     &        .or.(lakon(nelem)(2:4).eq.'RTA')) then
!     
         call absolute_relative(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,numf,set) 
!
      else
!     
c         write(*,*) '*WARNING in flux: no fluid section type given'
c         write(*,*) '         to fluid element ',nelem
         identity=.true.
!          
      endif
!    
      return
      end
      
