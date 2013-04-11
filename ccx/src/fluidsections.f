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
      subroutine fluidsections(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,ielmat,matname,nmat,
     &  kon,ipkon,irstrt,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,lakon,ielprop,nprop,nprop_,prop,iaxial,
     &  ipoinpc)
!
!     reading the input deck: *FLUID SECTION
!
      implicit none
!
      character*1 inpc(*)
      character*7 elname
      character*8 lakon(*)
      character*80 matname(*),material,typename
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),ielmat(*),iaxial,
     &  kon(*),ipkon(*),irstrt,nset,nmat,ndprop,npropstart,
     &  istep,istat,n,key,i,j,k,imaterial,ipos,lprop,ipoinpc(0:*),
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ielprop(*),nprop,nprop_,
     &  nodea,nodeb   
!
      real*8 prop(*)
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR in fluidsections: *FLUID SECTION should'
         write(*,*) '  be placed before all step definitions'
         stop
      endif
!
      typename='
     &                        '
!
      do i=2,n
         if(textpart(i)(1:9).eq.'MATERIAL=') then
            material=textpart(i)(10:89)
         elseif(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         elseif(textpart(i)(1:5).eq.'TYPE=') then
            read(textpart(i)(6:85),'(a80)',iostat=istat) typename
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
         endif
      enddo
!
!     types of sections
!
      if(typename(1:18).eq.'ABSOLUTETORELATIVE') then
         elname='ATR'
         ndprop=3
!
      elseif(typename(1:18).eq.'RELATIVETOABSOLUTE') then
         elname='RTA'
         ndprop=3
!
      elseif(typename(1:12).eq.'BLEEDTAPPING') then
         elname='ORBT'
         ndprop=3
!
      elseif(typename(1:10).eq.'CARBONSEAL') then
         ndprop=3
         elname='CARBS   '
!
      elseif(typename(1:14).eq.'CHARACTERISTIC') then
         ndprop=20
         elname='CHAR   '
!
      elseif(typename(1:7).eq.'GASPIPE') then
         if(typename(8:16).eq.'ADIABATIC') then
            elname='GAPIA   '
            ndprop=6
         elseif(typename(8:17).eq.'ISOTHERMAL') then
            elname='GAPII   '
            ndprop=6
         endif
!
      elseif(typename(1:9).eq.'LABYRINTH') then
         elname='LAB     '
         ndprop=10 
!
      elseif(typename(1:10).eq.'LIQUIDPUMP') then
         ndprop=20
         elname='LIPU   '
!
      elseif(typename(1:8).eq.'MOEHRING') then
         if(typename(9:19).eq.'CENTRIFUGAL') then
            elname='MRGF'
         elseif(typename(9:19).eq.'CENTRIPETAL') then
            elname='MRGP'
         endif
         ndprop=8
!
      elseif(typename(1:7).eq.'ORIFICE') then
         ndprop=8
         if(typename(9:14).eq.'OWN_AL') then
            elname='ORMA   '
         elseif(typename(9:13).eq.'PK_AL') then
            elname='ORPA   '
         elseif(typename(9:13).eq.'MS_MS') then
            elname='ORMM   '
         elseif(typename(9:13).eq.'PK_MS') then
            elname='ORPM   '
         else 
            elname='ORC1   '
         endif
!
      elseif(typename(1:4).eq.'PIPE') then
         if(typename(5:8).eq.'BEND') then
            elname='LIPIBE '
            ndprop=4
         elseif(typename(5:10).eq.'BRANCH') then
            elname='LIPIBR '
            ndprop=6
         elseif(typename(5:15).eq.'CONTRACTION') then
            elname='LIPICO '
            ndprop=2
         elseif(typename(5:13).eq.'DIAPHRAGM') then
            elname='LIPIDI '
            ndprop=2
         elseif(typename(5:15).eq.'ENLARGEMENT') then
            elname='LIPIEL '
            ndprop=2
         elseif(typename(5:12).eq.'ENTRANCE') then
            elname='LIPIEN '
            ndprop=2
         elseif(typename(5:13).eq.'GATEVALVE') then
            elname='LIPIGV '
            ndprop=2
         elseif(typename(5:11).eq.'MANNING') then
            if(typename(12:19).eq.'FLEXIBLE') then
               elname='LIPIMAF'
               ndprop=4
            else
               elname='LIPIMA '
               ndprop=3
            endif
         elseif(typename(5:19).eq.'WHITE-COLEBROOK') then
            if(typename(20:27).eq.'FLEXIBLE') then
               elname='LIPIWCF'
               ndprop=4
            else
               elname='LIPIWC '
               ndprop=3
            endif
         endif
!
      elseif(typename(1:14).eq.'PRESWIRLNOZZLE') then
         elname='ORPN'
         ndprop=5
!      
      elseif(typename(1:10).eq.'RESTRICTOR') then
         if(typename(11:15).eq.'USER') then
            elname='REUS'
            ndprop=4
         elseif(typename(11:29).eq.'LONGORIFICEIDELCHIK') then
            elname='RELOID'
            ndprop=4
         elseif(typename(11:21).eq.'WALLORIFICE') then
            elname='REWAOR'
            ndprop=3
         elseif(typename(11:18).eq.'ENTRANCE') then
            elname='REEN'
            ndprop=4 
         elseif(typename(11:21).eq.'ENLARGEMENT') then
            elname='REEL'
            ndprop=4
         elseif(typename(11:21).eq.'CONTRACTION') then
            elname='RECO'
            ndprop=5
         elseif(typename(11:26).eq.'BENDIDELCHIKCIRC') then
            elname='REBEIDC'
            ndprop=7
         elseif(typename(11:26).eq.'BENDIDELCHIKRECT') then
            elname='REBEIDR'
            ndprop=7
         elseif(typename(11:20).eq.'BENDMILLER') then
            elname='REBEMI'
            ndprop=5
         elseif(typename(11:17).eq.'BENDMAN') then
            elname='REBEMA'
            ndprop=5
         elseif(typename(11:14).eq.'EXIT') then
            elname='REEX'
            ndprop=4
         elseif(typename(11:33).eq.'LONGORIFICELICHTAROWICZ') then
            elname='RELOLI'
            ndprop=5
         endif
!     
      elseif(typename(1:6).eq.'BRANCH') then
         if(typename(7:13).eq.'JOINTGE')then
            elname='REBRJG'
            ndprop=8
         elseif(typename(7:20).eq.'JOINTIDELCHIK1')then
            elname='REBRJI1'
            ndprop=8
         elseif(typename(7:20).eq.'JOINTIDELCHIK2')then
            elname='REBRJI2'
            ndprop=8  
         elseif(typename(7:13).eq.'SPLITGE')then
            elname='REBRSG'
            ndprop=8  
         elseif(typename(7:20).eq.'SPLITIDELCHIK1')then
            elname='REBRSI1'
            ndprop=8
         elseif(typename(7:20).eq.'SPLITIDELCHIK2')then
            elname='REBRSI2'
            ndprop=8
         endif
!
      elseif(typename(1:6).eq.'VORTEX') then
         if(typename(7:10).eq.'FREE') then
            elname='VOFR'
            ndprop=9
         elseif(typename(7:12).eq.'FORCED') then
            elname='VOFO'
            ndprop=7
         endif
!     
      elseif(typename(1:1).eq.' ') then
         elname='       '
         ndprop=0
      else
         write(*,*) '*ERROR in fluidsections: ',typename
         write(*,*) '       is an unknown fluid section type'
         stop
      endif
!
!     check for the existence of the set and the material
!
      do i=1,nmat
         if(matname(i).eq.material) exit
      enddo
      if(i.gt.nmat) then
         write(*,*) '*ERROR in fluidsections: nonexistent material'
         write(*,*) '  '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
      imaterial=i
!
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR in fluidsections: element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline)
         stop
      endif
!
!     assigning the elements of the set the appropriate material
!     and property pointer
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:1).ne.'D') then
               write(*,*) '*ERROR in fluidsections: element ',
     &            ialset(j),' is no fluid element.'
               stop
            endif
            lakon(ialset(j))(2:8)=elname(1:7)
            ielmat(ialset(j))=imaterial
            ielprop(ialset(j))=nprop
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:1).ne.'D') then
                  write(*,*) '*ERROR in fluidsections: element ',
     &                 k,' is no fluid element.'
                  stop
               endif
               lakon(k)(2:8)=elname(1:7)
               ielmat(k)=imaterial
               ielprop(k)=nprop
            enddo
         endif
      enddo
!
      npropstart=nprop
!
!     storing the element properties
!
      if((elname(1:4).eq.'LIPU').or.(elname(1:4).eq.'CHAR')) then
         ndprop=0
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            do j=1,n
               ndprop=ndprop+1
               if(ndprop.gt.20) then
                  write(*,*) '*ERROR in fluidsections: no more than'
                  write(*,*) '       10 data points are allowed for a'
                  write(*,*) '       liquid pump or 9 data points for'
                  write(*,*) '       a gas characteristic'
                  stop
               endif
               read(textpart(j),'(f40.0)',iostat=istat) 
     &              prop(nprop+ndprop+1)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               if(ndprop.gt.2) then
                  if(2*(ndprop/2).ne.ndprop) then
                     if((prop(nprop+ndprop+1).le.prop(nprop+ndprop-1))             
     &                    .and.(elname(1:4).ne.'CHAR')) then
                        write(*,*) '*ERROR in fluidsections: volumetric'
                        write(*,*) '       flow data points must be in'
                        write(*,*) '       strict increasing order'
                        stop
                     endif
                  else
                     if((prop(nprop+ndprop+1).gt.prop(nprop+ndprop-1))
     &                 .and.(elname(1:4).ne.'CHAR')) then
                        write(*,*) '*ERROR in fluidsections: total'
                        write(*,*) '       head data points must be '
                        write(*,*) '       decreasing'
                        stop
                     endif
                  endif
               endif
            enddo
         enddo
         if(2*(ndprop/2).ne.ndprop) then
            if(elname(1:4).eq.'LIPU') then
               write(*,*) '*ERROR in fluidsections: less head data'
               write(*,*) '       points than volumetric flow'
               write(*,*) '       data points'
               stop
            else
               write(*,*) '*ERROR in fluidsections: less reduced flow'
               write(*,*) '       data points than pressure ratio data'
               write(*,*) '       points'
               stop
            endif
         endif
         prop(nprop+1)=ndprop/2
         if(iaxial.ne.0) then
            do j=1,ndprop/2,2
               prop(nprop+j+1)=prop(nprop+j+1)/iaxial
            enddo
            if(elname(1:4).eq.'CHAR') prop(1)=prop(1)*iaxial
         endif
         nprop=nprop+ndprop+1
      elseif(ndprop.gt.0) then
         if((elname.eq.'LIPIMAF').or.(elname.eq.'LIPIWCF')) then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
            read(textpart(1),'(i10)',iostat=istat) nodea
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            prop(nprop+1)=nodea+0.5d0
            read(textpart(2),'(i10)',iostat=istat) nodeb
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            prop(nprop+2)=nodeb+0.5d0
            read(textpart(3),'(f40.0)',iostat=istat)
     &            prop(nprop+3)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
            prop(nprop+4)=iaxial+0.5d0
            nprop=nprop+ndprop
         else
            lprop=0
            do j=1,(ndprop-1)/8+1
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               do k=1,8
                  lprop=lprop+1
                  if(lprop.gt.ndprop) exit
                  read(textpart(k),'(f40.0)',iostat=istat)
     &                 prop(nprop+lprop)
                  if(istat.gt.0) call inputerror(inpc,ipoinpc,iline)
               enddo
            enddo
!     
!           reducing the area in case of an axisymmetric model
!
            if(iaxial.ne.0) prop(nprop+1)=prop(nprop+1)/iaxial
            nprop=nprop+ndprop
         endif
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      else
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
      endif
!
      if(nprop.gt.nprop_) then
         write(*,*) '*ERROR in fluidsections: increase nprop_'
         stop
      endif
!
      if((elname(1:7).eq.'REBEIDC').or.
     &     (elname(1:7).eq.'REBEMI ').or.
     &     (elname(1:7).eq.'REBEMA ')) then
         if(elname(1:7).eq.'REBEIDC') then
            elname='REBEID '
         endif
         prop(npropstart+5)=prop(npropstart+3)
         prop(npropstart+4)=prop(npropstart+2)
         prop(npropstart+2)=prop(npropstart+1)
         prop(npropstart+3)=dsqrt(4*prop(npropstart+1)
     &        /(4.d0*datan(1.d0)))
!
      elseif(elname(1:7).eq.'REBEIDR') then
         elname='REBEID '
         prop(npropstart+7)=prop(npropstart+2)
         prop(npropstart+6)=prop(npropstart+1)
         prop(npropstart+5)=prop(npropstart+4)
         prop(npropstart+4)=prop(npropstart+3)
         prop(npropstart+1)=prop(npropstart+7)*prop(npropstart+6)
         prop(npropstart+2)=prop(npropstart+1)
         prop(npropstart+3)=2*prop(npropstart+1)/(prop(npropstart+6)+
     &        prop(npropstart+7))
!     
      elseif(elname(1:4).eq.'REEX') then
            prop(npropstart+4)=prop(npropstart+3)
            prop(npropstart+3)=prop(npropstart+2)
            prop(npropstart+2)=100000*prop(npropstart+1)
!     
      elseif(elname(1:4).eq.'REEN') then
         prop(npropstart+4)=0.5d0
         prop(npropstart+3)=prop(npropstart+2)
         prop(npropstart+2)=prop(npropstart+1)
         prop(npropstart+1)=100000*prop(npropstart+2)
!
      elseif(elname(1:7).eq.'REBRJI1') then
         prop(npropstart+8)= prop(npropstart+6)
         prop(npropstart+7)= 0d0
         prop(npropstart+6)=prop(npropstart+5)
         prop(npropstart+5)=prop(npropstart+4)
!
      elseif(elname(1:7).eq.'REBRJI2') then
         prop(npropstart+8)=prop(npropstart+7)
         prop(npropstart+7)=0.d0
         if(prop(npropstart+5)+prop(npropstart+6).ne.
     &        prop(npropstart+4))then
            write(*,*) '*ERROR: in fluidsections:'
            write(*,*) '        in element type RESTRICTOR 
     &                           BRANCH JOINT IDELCHIK2'
            write(*,*) '        A0 ist not equal to A1+A2'
         endif
!
      elseif(elname(1:7).eq.'REBRSI1') then
         prop(npropstart+8)=prop(npropstart+6)
         prop(npropstart+7)=0.d0
         prop(npropstart+6)=prop(npropstart+5)
         prop(npropstart+5)=prop(npropstart+4)
!     
      elseif(elname(1:7).eq.'REBRSI2') then
         prop(npropstart+8)=90.d0
         prop(npropstart+7)=90.d0
!     
      endif
!
!     "transforming" real inputs into integer inputs
!
      if(elname(1:4).eq.'REBR') then
!
!     for branches element the three first inputs should be integers 
!     since they are respectively the label number of the element
!     modeling branch 0 , 1 and 2.
!     Since all inputs are real the "integer input" are augmented by 0.1
!     and will be converted in integer using fucntion INT
         prop(npropstart+3)=prop(npropstart+3)+0.1d0
         prop(npropstart+2)=prop(npropstart+2)+0.1d0
         prop(npropstart+1)=prop(npropstart+1)+0.1d0
!
!     number of the table to be used for Preswirlnozzles
!
      elseif(elname(1:4).eq.'ORPN') then
         prop(npropstart+4)=prop(npropstart+4)+0.1d0
!
!     number of the table to be used for bleed tappings
!
      elseif(elname(1:4).eq.'ORBT') then
         prop(npropstart+3)=prop(npropstart+3)+0.1d0
!
!     label number of the swirl generating element for free vortices
!
      elseif(elname(1:4).eq.'VOFR') then
         if(prop(npropstart+6).ne.0d0) then
            prop(npropstart+6)=prop(npropstart+6)+0.1d0
         endif
!
!     label number of the pipe element preceding the exit loss element
!
      elseif(elname(1:4).eq.'REEX') then
         prop(npropstart+3)= prop(npropstart+3)+0.1d0
!
!     label number of the upstream and downstream nodes for Moehrings
!
      elseif(elname(1:2).eq.'MR') then
         prop(npropstart+5)= prop(npropstart+5)+0.1d0
         prop(npropstart+6)= prop(npropstart+6)+0.1d0
      endif         
!     
!     check the range of the parameters
!
      if((elname(1:4).ne.'LIPU').and.(elname.ne.'       ').and.
     &   (elname(1:3).ne.'LAB').and.(elname(1:4).ne.'CHAR').and.
     &   (elname(1:5).ne.'CARBS')) then
         if(prop(npropstart+1).le.0.d0) then
            write(*,*) '*ERROR in fluidsections: section area'
            write(*,*) '       is not positive'
            stop
         endif
      endif
!
      if((elname(1:2).eq.'OR').and.(elname(1:4).ne.'ORC1')) then
         if(prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: diameter of the'
            write(*,*) '       orifice is not positive'
            stop
         endif
         if(prop(npropstart+3).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: length of the'
            write(*,*) '       orifice is not positive'
            stop
         endif
         if((prop(npropstart+4).gt.1.d-20).and.
     &      (prop(npropstart+5).gt.1.d-20)) then
            write(*,*) '*ERROR in fluidsections: either the radius'
            write(*,*) '       of the orifice must be zero or the'
            write(*,*) '       chamfer angle'
            stop
         endif
         if(prop(npropstart+4)/prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: r/d of an orifice'
            write(*,*) '       must not be negative'
            stop
         endif
         if(prop(npropstart+4)/prop(npropstart+2).gt.0.82) then
            write(*,*) '*ERROR in fluidsections: r/d of an orifice'
            write(*,*) '       must not not exceed 0.82'
            stop
         endif
         if(prop(npropstart+5).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: the chamfer angle'
            write(*,*) '       of an orifice must not be negative'
            stop
         endif
         if(prop(npropstart+5).gt.90.d0) then
            write(*,*) '*ERROR in fluidsections: the chamfer angle'
            write(*,*) '       of an orifice must not exceed 90°'
            stop
         endif
         if(prop(npropstart+6).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: d/D (orifice)'
            write(*,*) '       must not be negative'
            stop
         endif
         if(prop(npropstart+6).gt.0.5d0) then
            write(*,*) '*ERROR in fluidsections: d/D (orifice)'
            write(*,*) '       must not exceed 0.5'
            stop
         endif
         if(prop(npropstart+3)/prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: L/d of an orifice'
            write(*,*) '       must not be negative'
            stop
         endif
         if(elname(3:3).eq.'P') then
            if(prop(npropstart+3)/prop(npropstart+2).gt.2.d0) then
               write(*,*) '*ERROR in fluidsections: L/d of an orifice'
               write(*,*) '       with Parker must not exceed 2'
               stop
            endif
         else
            if(prop(npropstart+3)/prop(npropstart+2).gt.2.d0) then
               write(*,*) '*ERROR in fluidsections: L/d of an orifice'
               write(*,*) '        must not exceed 10'
               stop
            endif
         endif
      elseif(elname(1:4).eq.'ORBT') then
         if(prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: ps1/pt1 (bleed'
            write(*,*) '       tapping) must not be negative'
            stop
         endif
         if(prop(npropstart+2).gt.1.d0) then
            write(*,*) '*ERROR in fluidsections: ps1/pt1 (bleed'
            write(*,*) '       tapping) must not exceed 1'
            stop
         endif
      elseif(elname(1:4).eq.'ORPN') then
         if(prop(npropstart+2).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: theta (preswirl'
            write(*,*) '       nozzle) must not be negative'
            stop
         endif
         if(prop(npropstart+2).gt.90.d0) then
            write(*,*) '*ERROR in fluidsections: theta (preswirl'
            write(*,*) '       nozzle) must not exceed 90°'
            stop
         endif
         if(prop(npropstart+3).lt.0.d0) then
            write(*,*) '*ERROR in fluidsections: k_phi (preswirl'
            write(*,*) '       nozzle) must not be negative'
            stop
         endif
         if(prop(npropstart+3).gt.1.05d0) then
            write(*,*) '*ERROR in fluidsections: k_phi (preswirl'
            write(*,*) '       nozzle) must not exceed 1.05'
            stop
         endif
      endif
!
      if(elname(1:3).eq.'LAB') then
         if((prop(npropstart+1).gt.1000.d0)
     &        .or.(prop(npropstart+1).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections: the selected pitch t'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=t<=1000 mm'
            stop
         endif
         if((prop(npropstart+2).gt.100d0)
     &        .or.(prop(npropstart+2).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections: the selected gap s'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=s<=100 mm'
            stop
         endif
         if((prop(npropstart+3).gt.5000.d0)
     &        .or.(prop(npropstart+3).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected diameter d'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=d<=5000'
            stop
         endif
         if((prop(npropstart+4).gt.9.d0)
     &        .or.(prop(npropstart+4).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*)'the selected spike number n'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=n<=9'
           stop
        endif
         if((prop(npropstart+5).gt.100.d0)
     &       .or.(prop(npropstart+5).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected spike breadth'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=b<=100 mm'
            stop
         endif
         if((prop(npropstart+6).gt.9.d0)
     &        .or.(prop(npropstart+6).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected spike height'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=b<=20 mm'
            stop
         endif
         if((prop(npropstart+7).gt.4.d0)
     &        .or.(prop(npropstart+7).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected Honeycomb cell width'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=L<=4 mm'
            stop
         endif
         if((prop(npropstart+8).gt.5.d0)
     &        .or.(prop(npropstart+8).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected edge radius'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=r<=5 mm'
            stop
         endif
         if((prop(npropstart+9).gt.100.d0)
     &        .or.(prop(npropstart+9).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected position of the spike'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=X<=0.1 mm'
            stop
         endif
         if((prop(npropstart+10).gt.100.d0)
     &        .or.(prop(npropstart+10).lt.0.d0)) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected height of the step'
            write(*,*) 'for the labyrinth is not correct'
            write(*,*) '0<=Hst<=0.1 mm'
            stop
         endif
      endif
!
      if(elname(1:5).eq.'CARBS') then
         if((prop(npropstart+1).le.0.d0)
     &        .or.(prop(npropstart+2).le.0.d0)
     &        .or.(prop(npropstart+3).le.0.d0) ) then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'the selected height of the diameter'
            write(*,*) 'or the selected lenght' 
            write(*,*) 'or the selected gap '
            write(*,*) 'has been defined as less or equal to 0'
            stop
         endif
      endif
!
      if(elname(1:2).eq.'RE') then
         if((prop(npropstart+1).le.0)
     &        .or.(prop(npropstart+2).le.0) 
     &        .or.(prop(npropstart+3).le.0))then
            write(*,*) '*ERROR in fluidsections:'
            write(*,*) 'A1,A2 or Dh less or equal 0'
            stop
!     
         elseif(elname(1:4).eq.'REEL') then
            if(prop(npropstart+1).ge.(prop(npropstart+2))) then
               write(*,*) '*ERROR in fluidsections:'
               write(*,*) 'Section A1 is greater than section A2'
               stop
            endif
!
         elseif(elname(1:4).eq.'RECO') then
            if(prop(npropstart+1).lt.(prop(npropstart+2))) then
               write(*,*) '*ERROR in fluidsections:'
               write(*,*) 'Section A2 is greater than section A1'
               stop
            endif
         endif
!     
         if(elname(1:4).eq.'REBR') then
            if((prop(npropstart+1).le.0)
     &           .or.(prop(npropstart+2).le.0) 
     &           .or.(prop(npropstart+3).le.0))then
               write(*,*) '*ERROR in fluidsections:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'all three elements must be different from 0'
               stop
!     
            elseif((prop(npropstart+4).le.0)
     &              .or.(prop(npropstart+5).le.0) 
     &              .or.(prop(npropstart+6).le.0))then 
               write(*,*) '*ERROR in fluidsections:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'all sections must be positive'
               stop
!
            elseif((prop(npropstart+7).lt.0)
     &              .or.(prop(npropstart+8).lt.0))then
               write(*,*) '*ERROR in fluidsections:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'alpha1 & 2 cannot be negative'
               stop
!     
            elseif((prop(npropstart+7).gt.90)
     &              .or.(prop(npropstart+8).gt.90)) then
               write(*,*) '*ERROR in fluidsections:'
               write(*,*) 'trying to define a branch '
               write(*,*) 'alpha1 & 2 cannot greater than 90 gegrees'
               stop
!     
            elseif((elname(5:7).eq.'SI1')
     &              .or.(elname(5:7).eq.'JI1')
     &              .or.(elname(5:7).eq.'JI2')) then
               if(prop(npropstart+7).gt.0) then
                  write(*,*) '*ERROR in fluidsections:'
                  write(*,*) 'trying to define a branch '
                  write(*,*) 'Type IDELCHIK JOINT1 or SPLIT1&2'
                  write(*,*) 'alpha1 must be 0 degrees'
                  stop
               endif
!
            elseif((elname(5:7).eq.'SI1').or.
     &              (elname(5:7).eq.'JI1'))then
               if(prop(npropstart+4).ne.(prop(npropstart+5)))then
                  write(*,*) '*ERROR in fluidsections:'
                  write(*,*) 'trying to define a branch '
                  write(*,*) 'Type IDELCHIK SPLIT1 or JOINT1 '
                  write(*,*) 'A1=A0'
                  stop
               endif
! 
            elseif(elname(5:7).eq.'JI2') then
               if((prop(npropstart+5)+(prop(npropstart+6)))
     &              .ne.prop(npropstart+4))then
                  write(*,*) '*ERROR in fluidsections:'
                  write(*,*) 'trying to define a branch '
                  write(*,*) 'Type IDELCHIK JOINT2 '
                  write(*,*) 'A1+A2 must be equal to A0'
                  stop
               endif
            endif
!
!     General Vortex
!     
         elseif(elname(1:2).eq.'VO') then
!
!     inner and outer radius less or equal to 0
            if((prop(npropstart+1).le.0) .or.
     &           (prop(npropstart+2).le.0)) then
               write(*,*)'*ERROR in fluidsections:'
               write(*,*)'trying to define a VORTEX'
               write(*,*)'R1 and R2 must be positive'
               stop
!     
            elseif(prop(npropstart+3).le.0) then
               write(*,*)'*ERROR in fluidsections:'
               write(*,*)'trying to define a VORTEX'
               write(*,*)'eta must be different positive'
               stop
            endif
!
!     FREE VORTEX
!
         elseif(elname(1:4).eq.'VOFR') then
!     
!     previous element must be refered to
!            if(prop(npropstart+9).le.0) then
!               write(*,*)'*ERROR in fluidsections:'
!               write(*,*)'trying to define a FREE VORTEX'
!               write(*,*)'no reference to previous upstream element'
!               stop
!            endif
!     
!     the swirl comes from another upstream  element
            if(prop(npropstart+6).ne.0d0) then
!     
!     the rotation speed must be 0
               if(prop(npropstart+7).ne.0) then
                  write(*,*)'*ERROR in fluidsections:'
                  write(*,*)'trying to define a FREE VORTEX'
                  write(*,*)'rotation speed and upstream element'
                  write(*,*)'cannot be simultaneously used '
                  stop
               endif
            endif
!
!     FORCED VORTEX
!
         elseif(elname(1:4).eq.'VOFO') then   
!
!     Core swirl ratio must be defined and positive
            if(prop(npropstart+4).le.0) then
               write(*,*)'*ERROR in fluidsections:'
               write(*,*)'trying to define a FORCED VORTEX'
               write(*,*)'Core swirl ratio Kr is strictly positive'
               stop
            endif
            if(prop(npropstart+4).gt.1) then
               write(*,*)'*ERROR in fluidsections:'
               write(*,*)'trying to define a FORCED VORTEX'
               write(*,*)'Core swirl ratio Kr cannot be greaer than 1'
               stop
            endif
!
!     Rotation speed must be defined and positive
            if(prop(npropstart+5).le.0) then
               write(*,*)'*ERROR in fluidsections:'
               write(*,*)'trying to define a FORCED VORTEX'
               write(*,*)'Rotation speed n is strictly positive'
            endif
         endif
!
!     Absolute/relative system
!
         if((elname(1:3).eq.'ATR').or.
     &        (elname(1:3).eq.'RTA')) then
            if(prop(npropstart+1).eq.0) then
               write(*,*)'*ERROR in fluidsections:'
               write(*,*)'trying to define an element' 
               write(*,*)'TYPE= ABSOLUTE TO RELATIVE or'
               write(*,*)'TYPE= RELATIVE TO ABSOLUTE'
               write(*,*)'Rotation speed is strictly positive'
               stop
            elseif(prop(npropstart+2).eq.0) then
               if(prop(npropstart+3).eq.0) then
                  write(*,*)'*ERROR in fluidsections:'
                  write(*,*)'trying to define an element' 
                  write(*,*)'TYPE= ABSOLUTE TO RELATIVE or'
                  write(*,*)'TYPE= RELATIVE TO ABSOLUTE'
                  write(*,*)'Tengential velocity is 0 but'
                  write(*,*)'no reference element has been provided'
                  stop
               endif
            elseif(prop(npropstart+3).ne.0) then
               if(prop(npropstart+2).ne.0) then
                  write(*,*)'*ERROR in fluidsections:'
                  write(*,*)'trying to define an element' 
                  write(*,*)'TYPE= ABSOLUTE TO RELATIVE or'
                  write(*,*)'TYPE= RELATIVE TO ABSOLUTE'
                  write(*,*)'reference element has been provided'
                  write(*,*)'but tengential velocity is already defined'
                  stop
               endif 
            endif
         endif
      endif
      do i=1,9
         write(*,*) prop(npropstart+i)
      enddo
      return
      end


