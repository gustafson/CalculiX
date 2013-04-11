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
      subroutine subspace(d,aa,bb,cc,tinc,alpham,betam,nev,xini,
     &  cd,cv,time,rwork,lrw,k,jout,rpar,bj,iwork,liw,iddebdf)
!
!     solves the linear dynamic equations mapped on the subspace
!     of the eigenvectors (only if there are dashpots in the
!     model)
!
      implicit none
!
      integer nev,nev2,info(15),idid,lrw,iwork(*),liw,jout,id,
     &  iab,iaa,ibb,i,j,k,iddebdf
!
      real*8 d(*),aa(nev,*),bb(nev,*),cc(nev,*),tinc,alpham,betam,
     &  xini(*),cd(*),cv(*),time0,time,rtol,atol,rwork(*),rpar(*),
     &  bj(*)
!
      external df,djac
!
      save time0
!
      nev2=2*nev
!
!     transferring fields into global field rpar
!     (needed for subroutine fd)
!     rpar contains (field, size): tinc, 1
!                                  alpham, 1
!                                  betam, 1
!                                  cc, nev**2
!                                  d, nev
!                                  aa(1,1)..aa(nev,1), nev
!                                  bb(1,1)..bb(nev,1), nev
!                                  aa(1,2)..aa(nev,2), nev
!                                  ...
!
      if(k.eq.0) then
         rpar(1)=tinc
         rpar(2)=alpham
         rpar(3)=betam
         do j=1,nev
            do i=1,nev
               rpar(3+(j-1)*nev+i)=cc(i,j)
            enddo
         enddo
         id=3+nev*nev
         do i=1,nev
            rpar(id+i)=d(i)
         enddo
         time0=0.d0
!
!        copying the initial conditions for the system of first order
!        differential equations
!
         do i=1,nev
            xini(i)=cd(i)
            xini(nev+i)=cv(i)
         enddo
      endif
      iab=3+nev*(1+nev)
      do j=k+1,k+jout
         iaa=iab+(j-1)*nev2
         ibb=iaa+nev
         do i=1,nev
            rpar(iaa+i)=aa(i,j)
            rpar(ibb+i)=bb(i,j)
         enddo
      enddo
!
      do i=1,6
         info(i)=0
      enddo
      info(5)=1
!
!     absolute and relative tolerance for dderkf
!
      rtol=1.d-5
      atol=1.d-3
!
      if(iddebdf.eq.0) then
         call ddeabm(df,nev2,time0,xini,time,info,rtol,atol,idid,rwork,
     &        lrw,iwork,liw,rpar,nev)
!
         if((idid.ne.2).and.(idid.ne.3)) then
            write(*,*) 
     &         '*WARNING in subspace: ddeabm did not converge properly'
            write(*,*) '         idid= ',idid
            write(*,*) '         switch to routine ddebdf'
            iddebdf=2
            time0=0.d0
!
!     copying the initial conditions for the system of first order
!     differential equations
!
            do i=1,nev
               xini(i)=cd(i)
               xini(nev+i)=cv(i)
            enddo
!     
            return
         endif
      else
         call ddebdf(df,nev2,time0,xini,time,info,rtol,atol,idid,rwork,
     &        lrw,iwork,liw,rpar,nev,djac)
         if((idid.ne.2).and.(idid.ne.3)) then
            write(*,*) 
     &           '*ERROR in subspace: ddebdf did not converge properly'
            write(*,*) '       idid= ',idid
            stop
         endif
      endif
!
!     copying the solution into field bj
!
      do i=1,nev
         bj(i)=xini(i)
      enddo
!
      return
      end
!
!     subroutine df expressing the first order derivative as a function
!     of time and the function itself
!
      subroutine df(x,u,uprime,rpar,nev)
!
      implicit none
!
      integer nev,i,j,k,id,iab,iaa,ibb
!
      real*8 rpar(*),x,u(*),uprime(*)
!
      k=int(x/rpar(1))+1
      id=3+nev*nev
      iab=id+nev
      iaa=iab+(k-1)*2*nev
      ibb=iaa+nev
!
      do i=1,nev
         uprime(i)=u(nev+i)
         uprime(nev+i)=rpar(iaa+i)+x*rpar(ibb+i)
     &             -rpar(id+i)*rpar(id+i)*u(i)
     &             -(rpar(2)+rpar(3)*rpar(id+i)*rpar(id+i))*u(nev+i)
c         uprime(nev+i)=aa(i,k)+x*bb(i,k)
c     &             -d(i)*d(i)*u(i)
c     &             -(apham+betam*d(i)*d(i))*u(nev+i)
!
!        contribution of the dashpots
!   
         do j=1,nev
            uprime(nev+i)=uprime(nev+i)-rpar(3+(j-1)*nev+i)*u(nev+j)
c            uprime(nev+i)=uprime(nev+i)-cc(i,j)*u(nev+j)
         enddo
      enddo
!
      return
      end
!
!     subroutine djac 
!
      subroutine djac(x,u,pd,nrowpd,rpar,nev)
!
      implicit none
!
      integer nrowpd,nev,id,i,j,k
!
      real*8 rpar(*),x,u(*),pd(nrowpd,*)
!
      k=int(x/rpar(1))+1
      id=3+nev*nev
!
      do i=1,nev
         pd(i,nev+i)=1.d0
         pd(nev+i,i)=-rpar(id+i)*rpar(id+i)
         pd(nev+i,nev+i)=-(rpar(2)+rpar(3)*rpar(id+i)*rpar(id+i))
!
!        contribution of the dashpots
!   
         do j=1,nev
            pd(nev+i,nev+j)=pd(nev+i,nev+j)-rpar(3+(j-1)*nev+i)
         enddo
      enddo
!
      return
      end
