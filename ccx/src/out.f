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
      subroutine out(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
     &  kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,xstaten,
     &  nstate_,istep,iinc,iperturb,ener,mint_,output,ithermal,qfn,
     &  mode,noddiam,trab,inotr,ntrans,orab,ielorien,norien,description,
     &  ipneigh,neigh,sti)
!
!     stores the results in frd format
!
      implicit none
!
      character*3 output
      character*6 filab(*)
      character*8 lakon(*)
      character*12 description
      character*80 matname(*)
!
      integer kon(*),inum(*),nk,ne,nmethod,kode,ipkon(*),mode,noddiam,
     &  ielmat(*),nstate_,istep,iinc,iperturb,mint_,ithermal,inotr(2,*),
     &  ntrans,ielorien(*),norien
!
      real*8 co(3,*),v(0:3,*),stn(6,*),een(6,*),t1(*),fn(0:3,*),time,
     &  epn(*),enern(*),xstaten(nstate_,*),ener(mint_,*),qfn(3,*),
     &  trab(7,*),orab(7,*)
!
      integer ipneigh(*),neigh(2,*)
      real*8 sti(6,mint_,*)
!
      if((output.eq.'frd').or.(output.eq.'FRD')) then
         call frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
     &        kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,
     &        xstaten,nstate_,istep,iinc,ithermal,qfn,mode,noddiam,
     &        trab,inotr,ntrans,orab,ielorien,norien,description,
     &        ipneigh,neigh,mint_,sti)
      else
         if(nmethod.ne.0) then
            call onf(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
     &        kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,
     &        xstaten,nstate_,istep,iinc,iperturb,ener,mint_)
         endif
      endif
!
      return
      end


