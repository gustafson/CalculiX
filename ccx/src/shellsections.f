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
      subroutine shellsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,ielorien,orname,norien,
     &  thicke,kon,ipkon,offset,irstrt,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,lakon,iaxial,ipoinpc)
!
!     reading the input deck: *SHELL SECTION
!
      implicit none
!
      logical nodalthickness
!
      character*1 inpc(*)
      character*8 lakon(*)
      character*80 matname(*),orname(*),material,orientation
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),ielmat(*),
     &  ielorien(*),kon(*),ipkon(*),indexe,irstrt,nset,nmat,norien,
     &  istep,istat,n,key,i,j,k,l,imaterial,iorientation,ipos,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),iaxial,ipoinpc(0:*)
!
      real*8 thicke(2,*),thickness,offset(2,*),offset1
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in shellsections: *SHELL SECTION should'
         write(*,*) '  be placed before all step definitions'
         stop
      endif
!
      nodalthickness=.false.
      offset1=0.d0
      orientation='                    '
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:12).eq.'ORIENTATION=') then
            orientation=textpart(i)(13:92)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         elseif(textpart(i)(1:14).eq.'NODALTHICKNESS') then
            nodalthickness=.true.
         elseif(textpart(i)(1:7).eq.'OFFSET=') then
            read(textpart(i)(8:27),'(f20.0)',iostat=istat) offset1
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         else
            write(*,*) 
     &        '*WARNING in shellsections: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline)
         endif
      enddo
!
!     check for the existence of the set,the material and orientation
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR in shellsections: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
      imaterial=i
!
      if(orientation.eq.'                    ') then
         iorientation=0
      else
         do i=1,norien
            if(orname(i).eq.orientation) exit
         enddo
         if(i.gt.norien) then
            write(*,*)'*ERROR in shellsections: nonexistent orientation'
            write(*,*) '  '
            call inputerror(inpc,ipoinpc,iline)
            stop
         endif
         iorientation=i
      endif
!
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR in shellsections: element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
!
!     assigning the elements of the set the appropriate material,
!     orientation number and offset
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:1).ne.'S') then
               write(*,*) '*ERROR in shellsections: *SHELL SECTION can'
               write(*,*) '       only be used for shell elements.'
               write(*,*) '       Element ',ialset(j),' is not a shell e
     &lement.'
               stop
            endif
            ielmat(ialset(j))=imaterial
            ielorien(ialset(j))=iorientation
            offset(1,ialset(j))=offset1
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:1).ne.'S') then
                  write(*,*) '*ERROR in shellsections: *SHELL SECTION ca
     &n'
                  write(*,*) '       only be used for shell elements.'
                  write(*,*) '       Element ',k,' is not a shell elemen
     &t.'
                  stop
               endif
               ielmat(k)=imaterial
               ielorien(k)=iorientation
               offset(1,k)=offset1
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
!     assigning a thickness to the elements
!
c      read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness
c      if(istat.gt.0) then
c         write(*,*) 
c     &        '*ERROR in shellsections: shell thickness is lacking'
c         call inputerror(inpc,ipoinpc,iline)
c      endif
!
      if(.not.nodalthickness) then
         read(textpart(1)(1:20),'(f20.0)',iostat=istat) thickness
         if(istat.gt.0) then
            write(*,*) 
     &           '*ERROR in shellsections: shell thickness is lacking'
            call inputerror(inpc,ipoinpc,iline)
         endif
         if(iaxial.ne.0) thickness=thickness/iaxial
         do j=istartset(i),iendset(i)
            if(ialset(j).gt.0) then
               indexe=ipkon(ialset(j))
               do l=1,8
                  thicke(1,indexe+l)=thickness
               enddo
            else
               k=ialset(j-2)
               do
                  k=k-ialset(j)
                  if(k.ge.ialset(j-1)) exit
                  indexe=ipkon(k)
                  do l=1,8
                     thicke(1,indexe+l)=thickness
                  enddo
               enddo
            endif
         enddo
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!
      return
      end

