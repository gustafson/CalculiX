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
      subroutine noelfiles(inpc,textpart,jout,filab,nmethod,
     &  nodefile_flag,elfile_flag,node_flag,nener,ithermal,
     &  istep,istat,n,iline,ipol,inl,ipoinp,inp,out3d,nlabel,
     &  amname,nam,itpamp,idrct,ipoinpc)
!
!     reading the *NODE FILE and *EL FILE cards in the input deck
!
      implicit none
!
      logical node_flag,nodefile_flag,elfile_flag,out3d,sectionforces
!
      character*1 nodesys,elemsys,inpc(*)
      character*6 filab(*)
      character*80 amname(*),timepointsname
      character*132 textpart(16)
!
      integer istep,istat,n,key,ii,jout,joutl,nmethod,nener,ithermal,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),j,nlabel,nam,itpamp,i,
     &  idrct,ipoinpc(0:*)
!
      save sectionforces
!
      nodesys='G'
      elemsys='L'
!
      if(istep.lt.1) then
         write(*,*) '*ERROR in noelfiles: *NODE FILE and *EL FILE'
         write(*,*) '       should only be used within a *STEP' 
         write(*,*) '       definition'
         stop
      endif
!
      if(node_flag) then
!
!        reset the nodal print requests
!
         if(.not.nodefile_flag) then
            filab(1)(1:4)='    '
            filab(2)(1:4)='    '
            filab(5)(1:4)='    '
            do j=10,12
               filab(j)(1:4)='    '
            enddo
            do j=14,17
               filab(j)(1:4)='    '
            enddo
!
            filab(1)(6:6)=' '
            filab(2)(6:6)=' '
            filab(5)(6:6)=' '
            do j=10,12
               filab(j)(6:6)=' '
            enddo
            do j=14,17
               filab(j)(6:6)=' '
            enddo
         endif
      else
!
!        reset the element print requests
!
         if(.not.elfile_flag) then
            filab(3)(1:4)='    '
            filab(4)(1:4)='    '
            do j=6,9
               filab(j)(1:4)='    '
            enddo
            filab(13)(1:4)='    '
!
            filab(3)(6:6)=' '
            filab(4)(6:6)=' '
            do j=6,9
               filab(j)(6:6)=' '
            enddo
            filab(13)(6:6)=' '
!
            sectionforces=.false.
         endif
      endif
!
      do ii=2,n
        if(textpart(ii)(1:10).eq.'FREQUENCY=') then
           read(textpart(ii)(11:20),'(i10)',iostat=istat) joutl
           if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
           if(joutl.eq.0) then
              do
                 call getnewline(inpc,textpart,istat,n,key,iline,ipol,
     &                inl,ipoinp,inp,ipoinpc)
                 if((key.eq.1).or.(istat.lt.0)) return
              enddo
           endif
           if(joutl.gt.0) jout=joutl
        elseif(textpart(ii)(1:10).eq.'GLOBAL=YES') then
           nodesys='G'
           elemsys='G'
        elseif(textpart(ii)(1:9).eq.'GLOBAL=NO') then
           nodesys='L'
           elemsys='L'
        elseif(textpart(ii)(1:9).eq.'OUTPUT=3D') then
           if(istep.eq.1) then
              out3d=.true.
              do j=1,nlabel
                 filab(j)(5:5)='E'
              enddo
           elseif(.not.out3d) then
              write(*,*) '*WARNING in noelfiles: OUTPUT=3D has no'
              write(*,*) '         effect in all but the first step'
           endif
        elseif(textpart(ii)(1:13).eq.'SECTIONFORCES') then
           if(out3d) then
              write(*,*) '*WARNING in noelfiles: SECTION FORCES cannot'
              write(*,*) '         be selected for 3D output'
           else
              filab(3)(5:5)='M'
           endif
        elseif(textpart(ii)(1:11).eq.'TIMEPOINTS=') then
           timepointsname=textpart(ii)(12:91)
           do i=1,nam
              if(amname(i).eq.timepointsname) then
                 itpamp=i
                 exit
              endif
           enddo
           if(i.gt.nam) then
              write(*,*) '*ERROR in noelfiles: time'
              write(*,*) '       points definition',
     &               timepointsname,' is unknown'
              stop
           endif
           if(idrct.eq.1) then
              write(*,*) '*ERROR in noelfiles: the DIRECT option'
              write(*,*) '       collides with a TIME POINTS '
              write(*,*) '       specification'
              stop
           endif
        endif
      enddo
!
      do
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((key.eq.1).or.(istat.lt.0)) return
         do ii=1,n
            if(textpart(ii)(1:4).eq.'U   ') then
               filab(1)(1:4)='U   '
               filab(1)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'NT  ') then
               filab(2)(1:4)='NT  '
               filab(2)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'S   ') then
               filab(3)(1:4)='S   '
               filab(3)(6:6)=elemsys
            elseif(textpart(ii)(1:4).eq.'E   ') then
               filab(4)(1:4)='E   '
               filab(4)(6:6)=elemsys
            elseif(textpart(ii)(1:4).eq.'RF  ') then
               filab(5)(1:4)='RF  '
               filab(5)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'PE  ') then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) '*WARNING in noelfiles: selection of PE'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
               else
                  filab(6)(1:4)='PE  '
                  filab(6)(6:6)=elemsys
               endif
            elseif(textpart(ii)(1:4).eq.'CE  ') then
               textpart(ii)(1:2)='PE'
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) '*WARNING in noelfiles: selection of CE'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
               else
                  write(*,*) '*WARNING in elprints: selection of CE'
                  write(*,*) '         is converted into PE; no distinct
     &ion'
                  write(*,*) '        is made between PE and CE'
                  filab(6)(1:4)='PE  '
                  filab(6)(6:6)=elemsys
               endif
            elseif(textpart(ii)(1:4).eq.'ENER') then
               filab(7)(1:4)='ENER'
               filab(7)(6:6)=elemsys
              nener=1
            elseif(textpart(ii)(1:4).eq.'SDV ') then
               if((nmethod.eq.2).or.(nmethod.eq.3)) then
                  write(*,*) '*WARNING in noelfiles: selection of SDV'
                  write(*,*) '         does not make sense for a'
                  write(*,*) '         frequency or bucking calculation'
               else
                  filab(8)(1:4)='SDV '
                  filab(8)(6:6)=elemsys
               endif
            elseif(textpart(ii)(1:4).eq.'HFL ') then
               if(ithermal.le.1) then
                  write(*,*) '*WARNING in nodeprints: HFL only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '          calculations'
               else
                  filab(9)(1:4)='HFL '
c                  if(.not.out3d) filab(9)(5:5)='I'
                  filab(9)(6:6)=elemsys
               endif
            elseif(textpart(ii)(1:4).eq.'RFL ') then
               if(ithermal.le.1) then
                  write(*,*) '*WARNING in nodeprints: RFL only makes '
                  write(*,*) '         sense for heat transfer '
                  write(*,*) '          calculations'
               else
                  filab(10)(1:4)='RFL '
                  filab(10)(6:6)=nodesys
               endif
            elseif(textpart(ii)(1:4).eq.'PU  ') then
               filab(11)(1:4)='PU  '
            elseif(textpart(ii)(1:4).eq.'PNT ') then
               filab(12)(1:4)='PNT '
            elseif(textpart(ii)(1:3).eq.'ZZS') then
               filab(13)(1:4)='ZZS '
               filab(13)(6:6)=elemsys
            elseif(textpart(ii)(1:4).eq.'TT  ') then
               filab(14)(1:4)='TT  '
               filab(14)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'MF  ') then
               filab(15)(1:4)='MF  '
               filab(15)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'PT  ') then
               filab(16)(1:4)='PT  '
               filab(16)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'TS  ') then
               filab(17)(1:4)='TS  '
               filab(17)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'PHS ') then
               filab(18)(1:4)='PHS '
               filab(18)(6:6)=elemsys
            elseif(textpart(ii)(1:4).eq.'MAXU') then
               filab(19)(1:4)='MAXU'
               filab(19)(6:6)=nodesys
            elseif(textpart(ii)(1:4).eq.'MAXS') then
               filab(20)(1:4)='MAXS'
               filab(20)(6:6)=elemsys
            else
               write(*,*) '*WARNING in noelfiles: label not applicable'
               write(*,*) '         or unknown; '
               call inputwarning(inpc,ipoinpc,iline)
            endif
         enddo
      enddo
!
      return
      end






