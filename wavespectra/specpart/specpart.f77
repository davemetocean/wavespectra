! -*- f90 -*
      real function minval(arr,n)
      integer i,n
      real arr(n)
      minval=arr(1)
      do i=1,n
            if (arr(i).lt.minval) minval=arr(i)
      enddo
      return 
      end

      real function maxval(arr,n)
      integer i,n
      real arr(n)
      maxval=arr(1)
      do i=1,n
            if (arr(i).gt.maxval) maxval=arr(i)
      enddo
      return 
      end

      integer function iminval(arr,n)
      integer i,n
      integer arr(n)
      iminval=arr(1)
      do i=1,n
            if (arr(i).lt.iminval) iminval=arr(i)
      enddo
      return 
      end


      
      subroutine partition(spec,ipart,nk,nth)

      integer ihmax
      parameter(ihmax=200)
      integer nspec,nnspec
      parameter(nnspec=600)
      integer nk,nth
      integer npart
!     ----------------------------------------------------------------
!       imi     i.a.   i   input discretized spectrum.
!       ind     i.a.   i   sorted addresses.
!       imo     i.a.   o   output partitioned spectrum.
!       zp      r.a.   i   spectral array.
!       npart   int.   o   number of partitions found.
!     ----------------------------------------------------------------
      
      real spec(nk,nth)
      real zp(nnspec)
      integer imi(nnspec),ind(nnspec),imo(nnspec)
      integer ipart(nk,nth), neigh(9,nnspec)
      integer iang,ik,nk,nth,nspec
      real zmin,zmax

      real maxval,minval

      npart=0
      nspec=nk*nth
      
      if (nspec.gt.nnspec) then 
            write(*,*) 'Error: Spectrum size exceeds maximum of ',nnspec
            stop
      endif

      call ptnghb(neigh,nspec,nk,nth)
      
      do iang=1, mth
        do ik=1, mk
            zp((iang-1)*mk+ik) = spec(ik,iang)
        enddo
      enddo
      zmin=minval(zp,nspec)
      zmax=maxval(zp,nspec)
      if (zmax-zmin.lt.1.e-9) then
        do iang=1, mth
        do ik=1, nk
        ipart(ik,iang)=0
        enddo
        enddo
        npart=0
        return
      endif

      fact   = real(ihmax-1) / ( zmax - zmin )

      do ik=1,nspec
      zp(ik)=zmax-zp(ik)
      imi(ik)    = max ( 1 , min ( ihmax , nint ( 1. + zp(ik)*fact ) ) )
      enddo
      

      call ptsort (ind, imi, nspec, ihmax)
      call pt_fld (neigh,zp,imo,imi,ind,nspec,npart,ihmax)
      
      do iang=1, mth
            do ik=1,nk
            ipart(ik,iang)=imo(1+(iang-1)*mk+ik)
            enddo
      enddo
      
      end

      subroutine ptsort(ind, imi, nspec, ihmax)
!
      implicit none
      integer                 i, in, iv, ihmax, iihmax, nspec, nnspec
      parameter(iihmax=200,nnspec=600)
      integer                 numv(iihmax), iaddr(iihmax),iorder(nnspec)
      integer                 imi(nspec), ind(nspec)
! -------------------------------------------------------------------- /
! 1.  occurences per height
!     
      do i=1, ihmax
      numv(i)   = 0
      enddo
      do i=1, nspec
        numv(imi(i)) = numv(imi(i)) + 1
      end do
!
! -------------------------------------------------------------------- /
! 2.  starting address per height
!
      iaddr(1) = 1
      do i=1, ihmax-1
        iaddr(i+1) = iaddr(i) + numv(i)
      end do
!
! -------------------------------------------------------------------- /
! 3.  order points
!
      do i=1, nspec
        iv        = imi(i)
        in        = iaddr(iv)
        iorder(i) = in
        iaddr(iv) = in + 1
      end do
!
! -------------------------------------------------------------------- /
! 4.  sort points
!
      do i=1, nspec
        ind(iorder(i)) = i
      end do
!
      return
!/
!/ end of ptsort ----------------------------------------------------- /
!/
      end

!/ ------------------------------------------------------------------- /
      subroutine ptnghb(neigh,nspec,mk,mth)

      integer  n, j, i, k, mk, mth, nspec
      integer neigh(9,nspec)

! -------------------------------------------------------------------- /
! 2.  build map
!
      do i=1,nspec
      do j=1,9
      neigh(j,i)  = 0
      enddo
      enddo
!
! ... base loop
!
      do n = 1, nspec
!
        j      = (n-1) / mk + 1
        i      = n - (j-1) * mk
        k      = 0
!
! ... point at the left(1)
!
        if ( i .ne. 1 ) then
            k           = k + 1
            neigh(k, n) = n - 1
          end if
!
! ... point at the right (2)
!
        if ( i .ne. mk ) then 
            k           = k + 1
            neigh(k, n) = n + 1
          end if
!
! ... point at the bottom(3)
!
        if ( j .ne. 1 ) then
            k           = k + 1
            neigh(k, n) = n - mk
          end if
!
! ... add point at bottom_wrap to top
!
        if ( j .eq. 1 ) then
            k          = k + 1
            neigh(k,n) = nspec - (mk-i)
          end if
!
! ... point at the top(4)
!
        if ( j .ne. mth ) then
            k           = k + 1
            neigh(k, n) = n + mk
          end if
!
! ... add point to top_wrap to bottom
!
         if ( j .eq. mth ) then
             k          = k + 1
             neigh(k,n) = n - (mth-1) * mk
            end if
!
! ... point at the bottom, left(5)
!
        if ( (i.ne.1) .and. (j.ne.1) ) then
            k           = k + 1
            neigh(k, n) = n - mk - 1
          end if
!
! ... point at the bottom, left with wrap.
!
         if ( (i.ne.1) .and. (j.eq.1) ) then
             k          = k + 1
             neigh(k,n) = n - 1 + mk * (mth-1)
           end if
!
! ... point at the bottom, right(6)
!
        if ( (i.ne.mk) .and. (j.ne.1) ) then
            k           = k + 1
            neigh(k, n) = n - mk + 1
          end if
!
! ... point at the bottom, right with wrap
!
        if ( (i.ne.mk) .and. (j.eq.1) ) then
            k           = k + 1
            neigh(k,n) = n + 1 + mk * (mth - 1)
          end  if
!
! ... point at the top, left(7)
!
        if ( (i.ne.1) .and. (j.ne.mth) ) then
            k           = k + 1
            neigh(k, n) = n + mk - 1
          end if
!
! ... point at the top, left with wrap
!
         if ( (i.ne.1) .and. (j.eq.mth) ) then
             k           = k + 1
             neigh(k,n) = n - 1 - (mk) * (mth-1)
           end if
!
! ... point at the top, right(8)
!
        if ( (i.ne.mk) .and. (j.ne.mth) ) then
            k           = k + 1
            neigh(k, n) = n + mk + 1
          end if
!
! ... point at top, right with wrap
!
!
        if ( (i.ne.mk) .and. (j.eq.mth) ) then
            k           = k + 1
            neigh(k,n) = n + 1 - (mk) * (mth-1)
          end if
!
        neigh(9,n) = k
!
        end do
!
      return
!/
!/ end of ptnghb ----------------------------------------------------- /
!/
      end
!/ ------------------------------------------------------------------- /
      subroutine pt_fld (neigh,zp,imo,imi,ind,nspec,npart,ihmax)


      integer neigh(9,nspec)
      real zp(nspec)
      integer imo(nspec), imi(nspec), ind(nspec)
      integer nspec
      integer npart,ihmax
!/
!/ ------------------------------------------------------------------- /
!/ local parameters
!/
      integer nnspec
      parameter(nnspec=300)
      integer                 mask, init, iwshed, imd(nnspec)
      integer                 ic_label, ifict_pixel, m, ih, msave
      integer                 ip, i, ipp, ic_dist, iempty, ippp
      integer                 jl, jn, ipt, j, ispec
      integer                 iq(nnspec), iq_start, iq_end
      real                    zpmax, ep1, diff

      real maxval
      integer iminval
      
! -------------------------------------------------------------------- /
! 0.  initializations
!
      mask        = -2
      init        = -1
      iwshed      =  0
      ic_label    =  0
      ifict_pixel = -100
!
      iq_start    =  1
      iq_end      =  1
!
      zpmax       = maxval ( zp, nspec )

      do i=1,nspec
            imo(i)         = init
            imd(i)          = 0
      enddo
!
! -------------------------------------------------------------------- /
! 1.  loop over levels
!
      m      =  1
!
      do ih=1, ihmax
        msave  = m
!
! 1.a pixels at level ih
!
10        continue
          ip     = ind(m)
          if ( imi(ip) .ne. ih ) goto 11
!
!     flag the point, if it stays flagge, it is a separate minimum.
!
          imo(ip) = mask
!
!     consider neighbors. if there is neighbor, set distance and add
!     to queue.
!
            do i=1, neigh(9,ip)
              ipp    = neigh(i,ip)
              if ( (imo(ipp).gt.0) .or. (imo(ipp).eq.iwshed) ) then
                imd(ip) = 1
                call fifo_add (ip, iq, iq_end, nspec)
                goto 11
              end if
            end do
!
            if ( m+1 .gt. nspec ) then
              goto 11
            else
              m = m + 1
            end if
!
      goto 10
!
11    continue
! 1.b process the queue
!
        ic_dist = 1
        call fifo_add (ifict_pixel, iq, iq_end, nspec)
!
20      continue
          call fifo_first (ip, iq, iq_start, nspec)
!
!     check for end of processing
!
          if ( ip .eq. ifict_pixel ) then
              call fifo_empty (iempty, iq_start, iq_end)
              if ( iempty .eq. 1 ) then
                  goto 21
                else
                  call fifo_add (ifict_pixel, iq, iq_end, nspec)
                  ic_dist = ic_dist + 1
                  call fifo_first (ip, iq, iq_start, nspec)
                end if
            end if
!
!     process queue
!
          do i=1, neigh(9,ip)
            ipp = neigh(i,ip)
!
!     check for labeled watersheds or basins
!
            if ( (imd(ipp).lt.ic_dist) .and. ( (imo(ipp).gt.0) .or. (imo(ipp).eq.iwshed) ) ) then
!
                if ( imo(ipp) .gt. 0 ) then
!
                  if ((imo(ip) .eq. mask) .or. (imo(ip) .eq. iwshed)) then
                        imo(ip) = imo(ipp)
                      else if (imo(ip) .ne. imo(ipp)) then
                        imo(ip) = iwshed
                      end if
!
                  else if (imo(ip) .eq. mask) then
!
                    imo(ip) = iwshed
!
                  end if
!
              else if ( (imo(ipp).eq.mask) .and. (imd(ipp).eq.0) ) then
!
                 imd(ipp) = ic_dist + 1
                 call fifo_add (ipp, iq, iq_end, nspec)
!
              end if
!
            end do
!
          goto 20
21         continue
!
! 1.c check for mask values in imo to identify new basins
!
        m = msave
!
30          continue
          ip     = ind(m)
          if ( imi(ip) .ne. ih ) goto 31
          imd(ip) = 0
!
          if (imo(ip) .eq. mask) then
!
! ... new label for pixel
!
              ic_label = ic_label + 1
              call fifo_add (ip, iq, iq_end, nspec)
              imo(ip) = ic_label
!
! ... and all connected to it ...
!
40          continue
                call fifo_empty (iempty, iq_start, iq_end)
                if ( iempty .eq. 1 ) goto 41
                call fifo_first (ipp, iq, iq_start, nspec )
!
                do i=1, neigh(9,ipp)
                  ippp   = neigh(i,ipp)
                  if ( imo(ippp) .eq. mask ) then
                      call fifo_add (ippp, iq, iq_end, nspec)
                      imo(ippp) = ic_label
                    end if
                  end do
!
            goto 40
41          continue
!
            end if
!
          if ( m + 1 .gt. nspec ) then
              goto 31
            else
              m = m + 1
            end if
!
          goto 30
31          continue
!
        end do
!
! -------------------------------------------------------------------- /
! 2.  find nearest neighbor of 0 watershed points and replace
!     use original input to check which group to affiliate with 0
!     soring changes first in imd to assure symetry in adjustment.
!
      do j=1, 5
        do ispec=1,nspec
        imd(ispec)    = imo(ispec)
        enddo
        do jl=1 , nspec
          ipt    = -1
          if ( imo(jl) .eq. 0 ) then
              ep1    = zpmax
              do jn=1, neigh (9,jl)
                diff   = abs ( zp(jl) - zp(neigh(jn,jl)))
                if ( (diff.le.ep1) .and. (imo(neigh(jn,jl)).ne.0) ) then
                    ep1    = diff
                    ipt    = jn
                  end if
                end do
              if ( ipt .gt. 0 ) imd(jl) = imo(neigh(ipt,jl))
            end if
          end do
          do ispec=1,nspec
            imo(ispec)    = imd(ispec)
          enddo
        if ( iminval(imo, nspec) .gt. 0 ) goto 60
        end do
60      continue
!
      npart = ic_label
!
      return
!
      end

!/ ------------------------------------------------------------------- /
      subroutine fifo_add ( iv, iq, iq_end, nspec )
!
!     add point to fifo queue.
!
      integer iv,iq_end,nspec
      integer iq(nspec)
!
      iq(iq_end) = iv
!
      iq_end = iq_end + 1
      if ( iq_end .gt. nspec ) iq_end = 1
!
      return
      end
!/ ------------------------------------------------------------------- /
      subroutine fifo_empty ( iempty, iq_start, iq_end )
!
!     check if queue is empty.
!
      integer iempty
      
!
      if ( iq_start .ne. iq_end ) then
        iempty = 0
      else
        iempty = 1
      end if
!
      return
      end
!/ ------------------------------------------------------------------- /
      subroutine fifo_first ( iv, iq, iq_start, nspec )
!
!     get point out of queue.
!
      integer iv
      integer iq(nspec)
!
      iv = iq(iq_start)
!
      iq_start = iq_start + 1
      if ( iq_start .gt. nspec ) iq_start = 1
!
      return
      end
!/
!/ end of pt_fld ----------------------------------------------------- /
!/
     

      
