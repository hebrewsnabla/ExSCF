

      subroutine indgrd (iout, iprint, indarr, ngrdt, ngrda, ngrdb,
     $     ngrdg, itype)

      implicit none


C +------------------------------------------------------------+
C |  indgrd   --   CAJH, 07.2011                               |
C |                                                            |
C |                                    (based on PHFB indgrd)  |
C |                                                            |
C |                                                            |
C |  Prepare the mapping array igrdt => igrdX. Here, igrdt     |
C |  should be interpreted as an indexing array that loops     |
C |  over the total number of grid points (ngrdt).             |
C |                                                            |
C |  Then, the array indarr (dimension ngrdt*3), stores the    |
C |  indices igrda, igrdb, igrdg for each and every value of   |
C |  igrdt.                                                    |
C |                                                            |
C |  We currently prepare the indexing array according to      |
C |    indarr(igrdt,1) = igrda                                 |
C |          (igrdt,2) = igrdb                                 |
C |          (igrdt,3) = igrdg                                 |
C |                                                            |
C +------------------------------------------------------------+
C |                                                            |
C |  CAJH, 08.2011                                             |
C |                                                            |
C |  Added the input argument itype, such that:                |
C |                                                            |
C |    itype = 1,  single integration is performed             |
C |          = 2,  double integration is performed             |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       indarr - mapping array igrdt => igrdX
C       ngrdt  - total number of grid points
C       ngrda  - number of grid points for alpha integration
C       ngrdb  - number of grid points for beta integration
C       ngrdg  - number of grid points for gamma integration
C       itype  - whether single / double integration is desired

      integer     iout, iprint, ngrdt, indarr(ngrdt,*)
      integer     ngrda, ngrdb, ngrdg, itype

C     other variables

      integer     ntst, cnt, i1a, i1b, i1g, i2a, i2b, i2g


C     Form indexing array.


      if ( itype .eq. 1 ) then

        ntst = ngrda * ngrdb * ngrdg

        if ( ngrdt .lt. ntst ) then
          call gauerr ('Call to indgrd is broken.')
        endif

        cnt = 1

        do 10 i1a = 1, ngrda
          do 20 i1b = 1, ngrdb
            do 30 i1g = 1, ngrdg
              indarr(cnt,1) = i1a
              indarr(cnt,2) = i1b
              indarr(cnt,3) = i1g

              cnt = cnt + 1
 30         continue
 20       continue
 10     continue


      elseif ( itype .eq. 2 ) then

        ntst = ngrda**2 * ngrdb**2 * ngrdg**2

        if ( ngrdt .lt. ntst ) then
          call gauerr ('Call to indgrd is broken.')
        endif

        cnt = 1

        do 110 i1a = 1, ngrda
          do 120 i1b = 1, ngrdb
            do 130 i1g = 1, ngrdg
              do 140 i2a = 1, ngrda
                do 150 i2b = 1, ngrdb
                  do 160 i2g = 1, ngrdg
                    indarr(cnt,1) = i1a
                    indarr(cnt,2) = i1b
                    indarr(cnt,3) = i1g
                    indarr(cnt,4) = i2a
                    indarr(cnt,5) = i2b
                    indarr(cnt,6) = i2g

                    cnt = cnt + 1
 160              continue
 150            continue
 140          continue
 130        continue
 120      continue
 110    continue

      endif


      return
      end


