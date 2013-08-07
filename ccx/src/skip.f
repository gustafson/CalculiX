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
      subroutine skip(nset,nalset,nload,nbody,
     &  nforc,nboun,nflow,nk,ne,nkon,
     &  mi,nmpc,memmpc_,nmat,ntmat_,npmat_,ncmat_,norien,ntrans,nam,
     &  nprint,nlabel,ncs_,ne1d,ne2d,infree,nmethod,
     &  iperturb,nener,iplas,ithermal,nstate_,iprestr,mcs,ntie,
     &  nslavs)
!
      implicit none
!
      integer nset,nalset,nload,nforc,nboun,nflow,nk,ne,nkon,mi(*),
     &  nmpc,memmpc_,nmat,ntmat_,npmat_,ncmat_,norien,ntrans,nam,
     &  nprint,nlabel,ncs_,ne1d,ne2d,infree(4),i,mt,
     &  nmethod,iperturb(*),nener,iplas,ithermal,nstate_,iprestr,i4,
     &  maxamta,mcs,ntie,nbody,nslavs
!
      character*1 c1
      character*3 c3
      character*4 c4
      character*5 c5
      character*8 c8
      character*20 c20
      character*80 c80
      character*81 c81
      character*87 c87
!
      real*8 r8
!
      mt=mi(2)+1
!
!        skipping the next entries
!     
      read(15)(c81,i=1,nset)
      read(15)(i4,i=1,nset)
      read(15)(i4,i=1,nset)
      do i=1,nalset
         read(15)i4
      enddo
      read(15)(r8,i=1,3*nk)
      read(15)(i4,i=1,nkon)
      read(15)(i4,i=1,ne)
      read(15)(c8,i=1,ne)
      read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(c1,i=1,nboun)
      read(15)(r8,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      if(nam.gt.0) read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(i4,i=1,nboun)
      read(15)(r8,i=1,nboun)
      read(15)(i4,i=1,nmpc)
      read(15)(c20,i=1,nmpc)
      read(15)(i4,i=1,nmpc)
      read(15)(i4,i=1,nmpc)
      read(15)(r8,i=1,nmpc)
      read(15)(i4,i=1,3*memmpc_)
      read(15)(r8,i=1,memmpc_)
      read(15)(i4,i=1,nforc)
      read(15)(i4,i=1,nforc)
      read(15)(r8,i=1,nforc)
      read(15)(i4,i=1,nforc)
      read(15)(i4,i=1,nforc)
      if(nam.gt.0) read(15)(i4,i=1,nforc)
      read(15)(r8,i=1,nforc)
      read(15)(i4,i=1,2*nload)
      read(15)(c5,i=1,nload)
      read(15)(r8,i=1,2*nload)
      if(nam.gt.0) read(15)(i4,i=1,2*nload)
      read(15)(r8,i=1,2*nload)
      read(15)(c81,i=1,nbody)
      read(15)(i4,i=1,2*nbody)
      read(15)(r8,i=1,7*nbody)
      read(15)(r8,i=1,7*nbody)
      if(iprestr.gt.0) read(15) (r8,i=1,6*mi(1)*ne)
c      read(15)(i4,i=1,2*nflow)
c      read(15)(r8,i=1,nflow)
c      if(nam.gt.0) read(15)(i4,i=1,nflow)
c      read(15)(r8,i=1,nflow)
      read(15)(c5,i=1,nprint)
      read(15)(c81,i=1,nprint)
      read(15)(c87,i=1,nlabel)
      read(15)(r8,i=1,(ncmat_+1)*ntmat_*nmat)
      read(15)(i4,i=1,2*nmat)
      read(15)(r8,i=1,2*ntmat_*nmat)
      read(15)(i4,i=1,nmat)
      read(15)(r8,i=1,4*ntmat_*nmat)
      read(15)(i4,i=1,nmat)
      read(15)(r8,i=1,7*ntmat_*nmat)
      read(15)(i4,i=1,2*nmat)
      read(15)(r8,i=1,7*ntmat_*nmat)
      read(15)(i4,i=1,2*nmat)
      read(15)(r8,i=1,nmat)
      read(15)(r8,i=1,3)
      if(npmat_.ne.0)then
         read(15)(r8,i=1,(2*npmat_+1)*ntmat_*nmat)
         read(15)(i4,i=1,(ntmat_+1)*nmat)
         read(15)(r8,i=1,(2*npmat_+1)*ntmat_*nmat)
         read(15)(i4,i=1,(ntmat_+1)*nmat)
      endif
      if(norien.ne.0)then
         read(15)(c80,i=1,norien)
         read(15)(r8,i=1,7*norien)
         read(15)(i4,i=1,mi(3)*ne)
      endif
      if(ntrans.ne.0)then
         read(15)(r8,i=1,7*ntrans)
         read(15)(i4,i=1,2*nk)
      endif
      if(nam.gt.0)then
         read(15)(c80,i=1,nam)
         read(15)(i4,i=1,3*nam-1)
         maxamta=2*i4
         read(15)i4
         read(15)(r8,i=1,maxamta)
      endif
      if(ithermal.gt.0)then
         read(15)(r8,i=1,nk)
         read(15)(r8,i=1,nk)
         if((ne1d.gt.0).or.(ne2d.gt.0))then
            read(15)(r8,i=1,2*nk)
            read(15)(r8,i=1,2*nk)
         endif
         if(nam.gt.0) read(15)(i4,i=1,nk)
         read(15)(r8,i=1,nk)
      endif
      read(15)(c80,i=1,nmat)
      read(15)(i4,i=1,mi(3)*ne)
      read(15)(r8,i=1,mt*nk)
      if((nmethod.eq.4).or.((nmethod.eq.1).and.(iperturb(1).ge.2))) 
     &     then
         read(15)(r8,i=1,mt*nk)
      endif
      read(15)(i4,i=1,nk)
      if((ne1d.gt.0).or.(ne2d.gt.0))then
         read(15)(i4,i=1,2*nkon)
         read(15)(r8,i=1,infree(1)-1)
         read(15)(i4,i=1,infree(2)-1)
         read(15)(r8,i=1,mi(3)*nkon)
         read(15)(r8,i=1,2*ne)
         read(15)(i4,i=1,infree(4))
         read(15)(i4,i=1,3*(infree(3)-1))
         read(15)(i4,i=1,infree(4))
      endif
      if(ntie.gt.0) then
         read(15)(c81,i=1,3*ntie)
         read(15)(r8,i=1,2*ntie)
      endif
      if(ncs_.gt.0)then
         read(15)(i4,i=1,ncs_)
      endif
      if(mcs.gt.0) then
         read(15)(r8,i=1,17*mcs)
      endif
      read(15)(r8,i=1,6*mi(1)*ne)
      read(15)(r8,i=1,6*mi(1)*ne)
      if(nener.eq.1) read(15)(r8,i=1,mi(1)*ne)
      if(nstate_.gt.0) read(15)(r8,i=1,nstate_*mi(1)*(ne+nslavs))
      read(15) (r8,i=1,27)
      read(15) (r8,i=1,2)
      read(15) c3
      read(15) r8
!
      return
      end
