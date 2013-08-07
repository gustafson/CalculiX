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
      subroutine cloads(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nodeforc,ndirforc,xforc,nforc,nforc_,iamforc,
     &  amname,nam,ntrans,trab,inotr,co,ikforc,ilforc,nk,
     &  cload_flag,istep,istat,n,iline,ipol,inl,ipoinp,inp,nam_,
     &  namtot_,namta,amta,nmethod,iaxial,iperturb,ipoinpc,
     &  maxsectors,idefforc)
!
!     reading the input deck: *CLOADS
!
      implicit none
!
      logical cload_flag,add,user
!
      character*1 inpc(*)
      character*80 amplitude,amname(*)
      character*81 set(*),noset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),nodeforc(2,*),
     &  nset,nforc,nforc_,istep,istat,n,i,j,k,l,iforcdir,key,
     &  iamforc(*),nam,iamplitude,ntrans,inotr(2,*),ipos,ikforc(*),
     &  ilforc(*),nk,iline,ipol,inl,ipoinp(2,*),inp(3,*),nam_,namtot,
     &  namtot_,namta(3,*),idelay,lc,nmethod,ndirforc(*),isector,
     &  iperturb,iaxial,ipoinpc(0:*),maxsectors,jsector,idefforc(*)
!
      real*8 xforc(*),forcval,co(3,*),trab(7,*),amta(2,*)
!
      iamplitude=0
      idelay=0
      lc=1
      isector=0
      user=.false.
      add=.false.
!
      if(istep.lt.1) then
         write(*,*) '*ERROR in cloads: *CLOAD should only be used'
         write(*,*) '  within a STEP'
         stop
      endif
!
      do i=2,n
         if((textpart(i)(1:6).eq.'OP=NEW').and.(.not.cload_flag)) then
            do j=1,nforc
               xforc(j)=0.d0
            enddo
         elseif(textpart(i)(1:10).eq.'AMPLITUDE=') then
            read(textpart(i)(11:90),'(a80)') amplitude
            do j=1,nam
               if(amname(j).eq.amplitude) then
                  iamplitude=j
                  exit
               endif
            enddo
            if(j.gt.nam) then
               write(*,*)'*ERROR in cloads: nonexistent amplitude'
               write(*,*) '  '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
            iamplitude=j
         elseif(textpart(i)(1:10).eq.'TIMEDELAY=') THEN
            if(idelay.ne.0) then
               write(*,*) '*ERROR in cloads: the parameter TIME DELAY'
               write(*,*) '       is used twice in the same keyword'
               write(*,*) '       '
               call inputerror(inpc,ipoinpc,iline)
               stop
            else
               idelay=1
            endif
            nam=nam+1
            if(nam.gt.nam_) then
               write(*,*) '*ERROR in cloads: increase nam_'
               stop
            endif
            amname(nam)='
     &                                 '
            if(iamplitude.eq.0) then
               write(*,*) '*ERROR in cloads: time delay must be'
               write(*,*) '       preceded by the amplitude parameter'
               stop
            endif
            namta(3,nam)=sign(iamplitude,namta(3,iamplitude))
            iamplitude=nam
            if(nam.eq.1) then
               namtot=0
            else
               namtot=namta(2,nam-1)
            endif
            namtot=namtot+1
            if(namtot.gt.namtot_) then
               write(*,*) '*ERROR cloads: increase namtot_'
               stop
            endif
            namta(1,nam)=namtot
            namta(2,nam)=namtot
            read(textpart(i)(11:30),'(f20.0)',iostat=istat) 
     &           amta(1,namtot)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         elseif(textpart(i)(1:9).eq.'LOADCASE=') then
            read(textpart(i)(10:19),'(i10)',iostat=istat) lc
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            if(nmethod.ne.5) then
               write(*,*) '*ERROR in cloads: the parameter LOAD CASE'
               write(*,*) '       is only allowed in STEADY STATE'
               write(*,*) '       DYNAMICS calculations'
               stop
            endif
         elseif(textpart(i)(1:7).eq.'SECTOR=') then
            read(textpart(i)(8:17),'(i10)',iostat=istat) isector
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            if((nmethod.le.3).or.(iperturb.gt.1)) then
               write(*,*) '*ERROR in cloads: the parameter SECTOR'
               write(*,*) '       is only allowed in MODAL DYNAMICS or'
               write(*,*) '       STEADY STATE DYNAMICS calculations'
               stop
            endif
            if(isector.gt.maxsectors) then
               write(*,*) '*ERROR in cloads: sector ',isector
               write(*,*) '       exceeds number of sectors'
               stop
            endif
            isector=isector-1
         elseif(textpart(i)(1:4).eq.'USER') then
            user=.true.
         else
            write(*,*) 
     &        '*WARNING in cloads: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
!
      if(user.and.(iamplitude.ne.0)) then
         write(*,*) '*WARNING: no amplitude definition is allowed'
         write(*,*) '          for concentrated loads defined by a'
         write(*,*) '          user routine'
         iamplitude=0
      endif
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) return
!
         read(textpart(2)(1:10),'(i10)',iostat=istat) iforcdir
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         if((iforcdir.lt.1).or.(iforcdir.gt.6)) then
            write(*,*) '*ERROR in cloads: nonexistent degree of freedom'
            write(*,*) '       '
            call inputerror(inpc,ipoinpc,iline)
            stop
         endif
         if(iforcdir.gt.3) iforcdir=iforcdir+1
!
         if(textpart(3)(1:1).eq.' ') then
            forcval=0.d0
         else
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) forcval
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            if(iaxial.ne.0) forcval=forcval/iaxial
         endif
!
!        dummy flux consisting of the first primes
!
         if(user) forcval=1.2357111317d0
!
         read(textpart(1)(1:10),'(i10)',iostat=istat) l
         if(istat.eq.0) then
            if(l.gt.nk) then
               write(*,*) '*ERROR in cloads: node ',l
               write(*,*) '       is not defined'
               stop
            endif
            if(lc.ne.1) then
               jsector=isector+maxsectors
            else
               jsector=isector
            endif
c            if(isector.ne.0) then
c               if(ntrans.ne.0) then
c                  if(inotr(1,l).ne.0) then
c                     write(*,*) '*ERROR in cloads: in node ',l
c                     write(*,*) '       a force is applied in a local'
c                     write(*,*) '       coordinate system in a sector'
c                     write(*,*) '       different from the basis sector'
c                     write(*,*) '       this is not allowed'
c                     stop
c                  endif
c               endif
c            endif
            call forcadd(l,iforcdir,forcval,nodeforc,ndirforc,xforc,
     &        nforc,nforc_,iamforc,iamplitude,nam,ntrans,trab,inotr,co,
     &        ikforc,ilforc,jsector,add,user,idefforc)
         else
            read(textpart(1)(1:80),'(a80)',iostat=istat) noset
            noset(81:81)=' '
            ipos=index(noset,' ')
            noset(ipos:ipos)='N'
            do i=1,nset
               if(set(i).eq.noset) exit
            enddo
            if(i.gt.nset) then
               noset(ipos:ipos)=' '
               write(*,*) '*ERROR in cloads: node set ',noset
               write(*,*) '  has not yet been defined. '
               call inputerror(inpc,ipoinpc,iline)
               stop
            endif
            do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  k=ialset(j)
                  if(lc.ne.1) then
                     jsector=isector+maxsectors
                  else
                     jsector=isector
                  endif
c                  if(isector.ne.0) then
c                     if(ntrans.ne.0) then
c                        if(inotr(1,k).ne.0) then
c                           write(*,*) '*ERROR in cloads: in node ',k
c                           write(*,*) 
c     &                      '       a force is applied in a local'
c                           write(*,*) 
c     &                      '       coordinate system in a sector'
c                           write(*,*)
c     &                      '       different from the basis sector'
c                           write(*,*) '       this is not allowed'
c                           stop
c                        endif
c                     endif
c                  endif
                  call forcadd(k,iforcdir,forcval,
     &               nodeforc,ndirforc,xforc,nforc,nforc_,iamforc,
     &               iamplitude,nam,ntrans,trab,inotr,co,ikforc,ilforc,
     &               jsector,add,user,idefforc)
               else
                  k=ialset(j-2)
                  do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
                     if(lc.ne.1) then
                        jsector=isector+maxsectors
                     else
                        jsector=isector
                     endif
c                     if(isector.ne.0) then
c                        if(ntrans.ne.0) then
c                           if(inotr(1,k).ne.0) then
c                              write(*,*) '*ERROR in cloads: in node ',k
c                              write(*,*) 
c     &                         '       a force is applied in a local'
c                              write(*,*) 
c     &                         '       coordinate system in a sector'
c                              write(*,*) 
c     &                         '       different from the basis sector'
c                              write(*,*) '       this is not allowed'
c                              stop
c                           endif
c                        endif
c                     endif
                     call forcadd(k,iforcdir,forcval,
     &                 nodeforc,ndirforc,xforc,nforc,nforc_,
     &                 iamforc,iamplitude,nam,ntrans,trab,inotr,co,
     &                 ikforc,ilforc,jsector,add,user,idefforc)
                  enddo
               endif
            enddo
         endif
      enddo
!
      return
      end

