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
      subroutine umat_main(amat,iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,
     &        orab,pnewdt,istep,iinc)
!
!     calculates stiffness and stresses for a user defined material
!     law
!
      implicit none
!
      character*80 amat
!
      integer ithermal,icmd,kode,ielas,iel,iint,nstate_,mint_,iorien,
     &  istep,iinc
!
      real*8 elconloc(21),stiff(21),emec(6),emec0(6),beta(6),stre(6),
     &  vj,t1l,dtime,xkl(3,3),xikl(3,3),vij,pgauss(3),orab(7,*),
     &  time,ttime,pnewdt
!
      real*8 xstate(nstate_,mint_,*),xstateini(nstate_,mint_,*)
!
      if(amat(1:8).eq.'ABAQUSNL') then
!
         call umat_abaqusnl(amat(9:80),iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc)
!
      elseif(amat(1:6).eq.'ABAQUS') then
!
         call umat_abaqus(amat(7:80),iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,istep,iinc)
!
      elseif(amat(1:10).eq.'ANISO_PLAS') then
!
         call umat_aniso_plas(amat(11:80),
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:11).eq.'ANISO_CREEP') then
!
         call umat_aniso_creep(amat(12:80),
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:13).eq.'ELASTIC_FIBER') then
!
         call umat_elastic_fiber(amat(14:80),
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:10).eq.'LIN_ISO_EL') then
!
         call umat_lin_iso_el(amat(11:80),
     &        iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:14).eq.'SINGLE_CRYSTAL') then
!
         call umat_single_crystal(amat(15:80),
     &        iel,iint,kode,elconloc,emec,
     &        emec0,beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,
     &        nstate_,xstateini,xstate,stre,stiff,iorien,pgauss,orab)
!
      elseif(amat(1:4).eq.'USER') then
!
         call umat_user(amat(5:80),iel,iint,kode,elconloc,emec,emec0,
     &        beta,xikl,vij,xkl,vj,ithermal,t1l,dtime,time,ttime,
     &        icmd,ielas,mint_,nstate_,xstateini,xstate,stre,stiff,
     &        iorien,pgauss,orab,pnewdt)
      else
         write(*,*) '*ERROR in umat: no user material subroutine'
         write(*,*) '       defined for material ',amat
         stop
      endif
!
      return
      end
