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
      subroutine restarts(istep,nset,nload,nforc, nboun,nk,ne,
     &  nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,mi,
     &  ntrans,ncs_,namtot_,ncmat_,mpcfree,maxlenmpc,
     &  ne1d,ne2d,nflow,nlabel,iplas,
     &  nkon,ithermal,nmethod,iperturb,nstate_,nener,set,istartset,
     &  iendset,ialset,co,kon,ipkon,lakon,nodeboun,ndirboun,iamboun,
     &  xboun,ikboun,ilboun,ipompc,nodempc,coefmpc,labmpc,ikmpc,ilmpc,
     &  nodeforc,ndirforc,iamforc,xforc,ikforc,ilforc,nelemload,iamload,
     &  sideload,xload,elcon,nelcon,rhcon,nrhcon,
     &  alcon,nalcon,alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,
     &  ielorien,trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
     &  ielmat,matname,prlab,prset,filab,vold,nodebounold,
     &  ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
     &  iponor,xnor,knor,thickn,thicke,offset,iponoel,inoel,rig,
     &  shcon,nshcon,cocon,ncocon,ics,sti,
     &  ener,xstate,jobnamec,infree,nnn,irstrt,inpc,textpart,istat,n,
     &  key,prestr,iprestr,cbody,ibody,xbody,nbody,xbodyold,
     &  ttime,qaold,cs,mcs,output,physcon,ctrl,typeboun,iline,ipol,inl,
     &  ipoinp,inp,fmpc,tieset,ntie,tietol,ipoinpc,nslavs,t0g,t1g)
!
      implicit none
!
      character*1 typeboun(*),inpc(*)
      character*3 output
      character*6 prlab(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
      character*80 orname(*),amname(*),matname(*)
      character*81 set(*),prset(*),tieset(3,*),cbody(*)
      character*87 filab(*)
      character*132 jobnamec(*),textpart(16)
!
      integer istep,nset,nload,nforc,nboun,nk,ne,nmpc,nalset,nmat,
     &  ntmat_,npmat_,norien,nam,nprint,mi(*),ntrans,ncs_,
     &  namtot_,ncmat_,mpcfree,ne1d,ne2d,nflow,nlabel,iplas,nkon,
     &  ithermal,nmethod,iperturb(*),nstate_,istartset(*),iendset(*),
     &  ialset(*),kon(*),ipkon(*),nodeboun(*),ndirboun(*),iamboun(*),
     &  ikboun(*),ilboun(*),ipompc(*),nodempc(*),ikmpc(*),ilmpc(*),
     &  nodeforc(*),ndirforc(*),iamforc(*),ikforc(*),ilforc(*),
     &  nelemload(*),iamload(*),nelcon(*),ipoinpc(0:*),
     &  nrhcon(*),nalcon(*),nplicon(*),nplkcon(*),ielorien(mi(3),*),
     &  inotr(*),
     &  namta(*),iamt1(*),ielmat(mi(3),*),nodebounold(*),ndirbounold(*),
     &  iponor(*),knor(*),iponoel(*),inoel(*),rig(*),
     &  nshcon(*),ncocon(*),ics(*),infree(*),nnn(*),
     &  nener,irestartstep,irestartread,irstrt,istat,n,i,key,
     &  iprestr,mcs,maxlenmpc,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ntie,ibody(*),nbody,nslavs
!
      real*8 co(*),xboun(*),coefmpc(*),xforc(*),xload(*),elcon(*),
     &  rhcon(*),alcon(*),alzero(*),plicon(*),plkcon(*),orab(*),
     &  trab(*),amta(*),t0(*),t1(*),prestr(*),veold(*),tietol(2,*),
     &  vold(*),xbounold(*),xforcold(*),xloadold(*),t1old(*),eme(*),
     &  xnor(*),thickn(*),thicke(*),offset(*),t0g(2,*),t1g(2,*),
     &  shcon(*),cocon(*),sti(*),ener(*),xstate(*),
     &  ttime,qaold(2),cs(17,*),physcon(*),
     &  ctrl(*),fmpc(*),xbody(*),xbodyold(*)
!
      irestartread=0
      irestartstep=0
!
      do i=2,n
         if(textpart(i)(1:4).eq.'READ') then
            irestartread=1
c            if(irestartstep.eq.0) irestartstep=1
         elseif(textpart(i)(1:5).eq.'STEP=') then
            read(textpart(i)(6:15),'(i10)',iostat=istat) irestartstep
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         elseif(textpart(i)(1:5).eq.'WRITE') then
            irstrt=1
         elseif(textpart(i)(1:10).eq.'FREQUENCY=') then
            read(textpart(i)(11:20),'(i10)',iostat=istat) irstrt
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         else
            write(*,*) 
     &        '*WARNING in restarts: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
!
      if(irestartread.eq.1) then
        call restartread(istep,nset,nload,nforc, nboun,nk,ne,
     &  nmpc,nalset,nmat,ntmat_,npmat_,norien,nam,nprint,mi,
     &  ntrans,ncs_,namtot_,ncmat_,mpcfree,maxlenmpc,
     &  ne1d,ne2d,nflow,nlabel,iplas,
     &  nkon,ithermal,nmethod,iperturb,nstate_,nener,set,istartset,
     &  iendset,ialset,co,kon,ipkon,lakon,nodeboun,ndirboun,iamboun,
     &  xboun,ikboun,ilboun,ipompc,nodempc,coefmpc,labmpc,ikmpc,ilmpc,
     &  nodeforc,ndirforc,iamforc,xforc,ikforc,ilforc,nelemload,iamload,
     &  sideload,xload,elcon,nelcon,rhcon,nrhcon,
     &  alcon,nalcon,alzero,plicon,nplicon,plkcon,nplkcon,orname,orab,
     &  ielorien,trab,inotr,amname,amta,namta,t0,t1,iamt1,veold,
     &  ielmat,matname,prlab,prset,filab,vold,nodebounold,
     &  ndirbounold,xbounold,xforcold,xloadold,t1old,eme,
     &  iponor,xnor,knor,thickn,thicke,offset,iponoel,inoel,rig,
     &  shcon,nshcon,cocon,ncocon,ics,sti,
     &  ener,xstate,jobnamec,infree,nnn,irestartstep,prestr,iprestr,
     &  cbody,ibody,xbody,nbody,xbodyold,ttime,qaold,cs,mcs,
     &  output,physcon,ctrl,typeboun,fmpc,tieset,ntie,tietol,nslavs,
     &  t0g,t1g)
      endif
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end


