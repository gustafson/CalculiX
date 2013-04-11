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
      subroutine out(co,nk,kon,ipkon,lakon,ne0,v,stn,inum,nmethod,
     &  kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,xstaten,
     &  nstate_,istep,iinc,iperturb,ener,mi,output,ithermal,qfn,
     &  mode,noddiam,trab,inotr,ntrans,orab,ielorien,norien,description,
     &  ipneigh,neigh,stx,vr,vi,stnr,stni,vmax,stnmax,ngraph,veold,ne,
     &  cs,set,nset,istartset,iendset,ialset)
!
!     stores the results in frd format
!
      implicit none
!
      character*3 output
      character*8 lakon(*)
      character*12 description
      character*80 matname(*)
      character*81 set(*)
      character*87 filab(*)
!
      integer kon(*),inum(*),nk,ne0,nmethod,kode,ipkon(*),mode,noddiam,
     &  ielmat(*),nstate_,istep,iinc,iperturb,mi(2),ithermal,inotr(2,*),
     &  ntrans,ielorien(*),norien,ngraph,ne,nset,istartset(*),
     &  iendset(*),ialset(*),ipneigh(*),neigh(2,*)
!
      real*8 co(3,*),v(0:mi(2),*),stn(6,*),een(6,*),t1(*),fn(0:mi(2),*),
     &  time,epn(*),enern(*),xstaten(nstate_,*),ener(mi(1),*),qfn(3,*),
     &  trab(7,*),orab(7,*),vr(0:mi(2),*),vi(0:mi(2),*),stnr(6,*),
     &  cs(17,*),stni(6,*),pi,vmax(0:3,*),stnmax(0:6,*),
     &  veold(0:mi(2),*),stx(6,mi(1),*)
!
      if((output.eq.'frd').or.(output.eq.'FRD')) then
         call frd(co,nk,kon,ipkon,lakon,ne0,v,stn,inum,nmethod,
     &        kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,
     &        xstaten,nstate_,istep,iinc,ithermal,qfn,mode,noddiam,
     &        trab,inotr,ntrans,orab,ielorien,norien,description,
     &        ipneigh,neigh,mi(1),stx,vr,vi,stnr,stni,vmax,
     &        stnmax,ngraph,veold,ener,ne,cs,set,nset,istartset,
     &        iendset,ialset)
      else
         if(nmethod.ne.0) then
            call onf(co,nk,kon,ipkon,lakon,ne0,v,stn,inum,nmethod,
     &        kode,filab,een,t1,fn,time,epn,ielmat,matname,enern,
     &        xstaten,nstate_,istep,iinc,iperturb,ener,mi(1))
         endif
      endif
!
      return
      end


