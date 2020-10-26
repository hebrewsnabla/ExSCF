

      subroutine prtdmt (dotwom, norbs, itype, onem, twom)

      implicit none


C +------------------------------------------------------------+
C |  prtdmt   ---   CAJH, 08.2011                              |
C |                                                            |
C |                                                            |
C |  Print the one- and the two-particle density matrices      |
C |  associated with the PHF state to external files:          |
C |                                                            |
C |    1pdm.dat,  2pdm.dat                                     |
C |                                                            |
C |  in a convenient format.                                   |
C |                                                            |
C |  Note that if the files exist in the current directory,    |
C |  they will become overwritten.                             |
C |                                                            |
C +------------------------------------------------------------+


C     input / output variables

C       dotwom - whether to print the two-particle density matrix
C       norbs  - number of basis functions
C       itype  - type of density matrix
C                ( 1 - RHF, 2 - UHF, 3 - GHF )
C       onem   - one-particle density matrix
C       twom   - two-particle density matrix

      real*8      onem(*), twom(*)
      integer     norbs, itype
      logical     dotwom

C     other variables

      real*8      thresh
      integer     i, j, k, l
      integer     indi, indj, ind
      integer     nosq, nocb
      logical     tst
      character   stri*1, strj*1

      parameter   ( thresh = 1.0d-12 )


C     format statements

 101  format (2(2X, I4, A1), 2(E30.15))
 102  format (2X, I4, 2X, I2)
 201  format (4(2X, I4, A1), 2(E30.15))
 202  format (2X, I4, 2X, I2)



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Print the one-particle density matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      inquire (file = './1pdm.dat', exist = tst)

      if ( tst ) then
        open (unit = 100, file = './1pdm.dat', status = 'old')
        close (100, status = 'delete')
      endif

      open (unit = 100, file = './1pdm.dat', status = 'new')


      write (100, *) '!'
      write (100, *) '! 1pdm.dat - prepared by PHF-SCF '
      write (100, *) '!'
      write (100, *) '! This file contains the one-particle density '
      write (100, *) '! matrix associated with the PHF state of '
      write (100, *) '! interest. '
      write (100, *) '!'
      write (100, *) '! The first line contains the number of basis '
      write (100, *) '! functions (in the basis selected) and the '
      write (100, *) '! type of the density matrix printed: '
      write (100, *) '!'
      write (100, *) '!    type = 1,  RHF '
      write (100, *) '!         = 2,  UHF '
      write (100, *) '!         = 3,  GHF '
      write (100, *) '!'
      write (100, *) '! The one-particle density matrix is then '
      write (100, *) '! printed with the format: '
      write (100, *) '!'
      write (100, *) '!   i, j, real( ), imag ( ) '
      write (100, *) '!'
      write (100, *) '! which stands for the matrix element: '
      write (100, *) '!'
      write (100, *) '!   < i! j > '
      write (100, *) '!'

      write (100, 102) norbs, itype

      write (100, *) '!'

      if ( itype .eq. 1 ) then

        do 11 j = 1, norbs
          do 12 i = 1, norbs
            ind = (i-1)*2*norbs + (j-1)*2 + 1

            if ( abs (onem(ind)) .gt. thresh .or.
     $           abs (onem(ind+1)) .gt. thresh ) then
              write (100, 101) j, 'a', i, 'a', onem(ind), onem(ind+1)
            endif
 12       continue
 11     continue

      elseif ( itype .eq. 2 ) then

        do 21 j = 1, norbs
          do 22 i = 1, norbs
            ind = (i-1)*2*norbs + (j-1)*2 + 1

            if ( abs (onem(ind)) .gt. thresh .or.
     $           abs (onem(ind+1)) .gt. thresh ) then
              write (100, 101) j, 'a', i, 'a', onem(ind), onem(ind+1)
            endif
 22       continue
 21     continue

        do 23 j = 1, norbs
          do 24 i = 1, norbs
            ind = (i-1)*2*norbs + (j-1)*2 + 1 + 2*norbs*norbs

            if ( abs (onem(ind)) .gt. thresh .or.
     $           abs (onem(ind+1)) .gt. thresh ) then
              write (100, 101) j, 'b', i, 'b', onem(ind), onem(ind+1)
            endif
 24       continue
 23    continue

      elseif ( itype .eq. 3 ) then

        do 31 j = 1, 2*norbs
          do 32 i = 1, 2*norbs

            indi = i
            indj = j

            stri = 'a'
            strj = 'a'

            if ( i .gt. norbs ) indi = i-norbs
            if ( j .gt. norbs ) indj = j-norbs

            if ( i .gt. norbs ) stri = 'b'
            if ( j .gt. norbs ) strj = 'b'

            ind = (i-1)*4*norbs + (j-1)*2 + 1

            if ( abs (onem(ind)) .gt. thresh .or.
     $           abs (onem(ind+1)) .gt. thresh ) then
              write (100, 101) indj, strj, indi, stri,
     $                         onem(ind), onem(ind+1)
            endif
 32       continue
 31     continue

      endif

      write (100, *) '!'

      close (100)



C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C     %  Print the two-particle density matrix  %
C     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


      if ( .not. dotwom ) then
        return
      endif

      nosq = norbs * norbs
      nocb = nosq * norbs

      inquire (file = './2pdm.dat', exist = tst)

      if ( tst ) then
        open (unit = 100, file = './2pdm.dat', status = 'old')
        close (100, status = 'delete')
      endif

      open (unit = 100, file = './2pdm.dat', status = 'new')


      write (100, *) '!'
      write (100, *) '! 2pdm.dat - prepared by PHF-SCF '
      write (100, *) '!'
      write (100, *) '! This file contains the two-particle density '
      write (100, *) '! matrix associated with the PHF state of '
      write (100, *) '! interest. '
      write (100, *) '!'
      write (100, *) '! The first line contains the number of basis '
      write (100, *) '! functions (in the basis selected) and the '
      write (100, *) '! type of the density matrix printed: '
      write (100, *) '!'
      write (100, *) '!    type = 1,  RHF '
      write (100, *) '!         = 2,  UHF '
      write (100, *) '!         = 3,  GHF '
      write (100, *) '!'
      write (100, *) '! The two-particle density matrix is then '
      write (100, *) '! printed with the format: '
      write (100, *) '!'
      write (100, *) '!   i, j, k, l, real( ), imag ( ) '
      write (100, *) '!'
      write (100, *) '! which stands for the matrix element: '
      write (100, *) '!'
      write (100, *) '!   1/2 * < i! j! k l > '
      write (100, *) '!'
      write (100, *) '!'
      write (100, *) '! NOTE: Not all elements are printed. We use '
      write (100, *) '! the fact that '
      write (100, *) '!   < i! j! k l >  =  - < j! i! k l > '
      write (100, *) '!                  =  - < i! j! l k > '
      write (100, *) '!                  =    < j! i! l k > '
      write (100, *) '!'

      write (100, 202) norbs, itype

      write (100, *) '!'

      if ( itype .eq. 1 ) then

        do 111 i = 1, norbs
          do 112 j = i+1, norbs
            do 113 k = 1, norbs
              do 114 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'a', k, 'a', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 114          continue
 113        continue
 112      continue
 111    continue

        do 115 i = 1, norbs
          do 116 j = 1, norbs
            do 117 k = 1, norbs
              do 118 l = 1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 2*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'b', k, 'b', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 118          continue
 117        continue
 116      continue
 115    continue

      elseif ( itype .eq. 2 ) then

        do 121 i = 1, norbs
          do 122 j = i+1, norbs
            do 123 k = 1, norbs
              do 124 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'a', k, 'a', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 124          continue
 123        continue
 122      continue
 121    continue

        do 125 i = 1, norbs
          do 126 j = i+1, norbs
            do 127 k = 1, norbs
              do 128 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 2*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'b', j, 'b', k, 'b', l, 'b',
     $                             twom(ind), twom(ind+1)
                endif
 128          continue
 127        continue
 126      continue
 125    continue

        do 129 i = 1, norbs
          do 130 j = 1, norbs
            do 131 k = 1, norbs
              do 132 l = 1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 4*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'b', k, 'b', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 132          continue
 131        continue
 130      continue
 129    continue


      elseif ( itype .eq. 3 ) then

        do 141 i = 1, norbs
          do 142 j = i+1, norbs
            do 143 k = 1, norbs
              do 144 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'a', k, 'a', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 144          continue
 143        continue
 142      continue
 141    continue

        do 145 i = 1, norbs
          do 146 j = i+1, norbs
            do 147 k = 1, norbs
              do 148 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 2*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'b', j, 'b', k, 'b', l, 'b',
     $                             twom(ind), twom(ind+1)
                endif
 148          continue
 147        continue
 146      continue
 145    continue

        do 149 i = 1, norbs
          do 150 j = 1, norbs
            do 151 k = 1, norbs
              do 152 l = 1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 4*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'b', k, 'b', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 152          continue
 151        continue
 150      continue
 149    continue

        do 153 i = 1, norbs
          do 154 j = i+1, norbs
            do 155 k = 1, norbs
              do 156 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 6*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'a', k, 'b', l, 'b',
     $                             twom(ind), twom(ind+1)
                endif
 156          continue
 155        continue
 154      continue
 153    continue

        do 157 i = 1, norbs
          do 158 j = i+1, norbs
            do 159 k = 1, norbs
              do 160 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 8*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'b', j, 'b', k, 'a', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 160          continue
 159        continue
 158      continue
 157    continue

        do 161 i = 1, norbs
          do 162 j = i+1, norbs
            do 163 k = 1, norbs
              do 164 l = 1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 10*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'a', k, 'a', l, 'b',
     $                             twom(ind), twom(ind+1)
                endif
 164          continue
 163        continue
 162      continue
 161    continue

        do 165 i = 1, norbs
          do 166 j = 1, norbs
            do 167 k = 1, norbs
              do 168 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 12*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'a', j, 'b', k, 'a', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 168          continue
 167        continue
 166      continue
 165    continue

        do 169 i = 1, norbs
          do 170 j = i+1, norbs
            do 171 k = 1, norbs
              do 172 l = 1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 14*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'b', j, 'b', k, 'b', l, 'a',
     $                             twom(ind), twom(ind+1)
                endif
 172          continue
 171        continue
 170      continue
 169    continue

        do 173 i = 1, norbs
          do 174 j = 1, norbs
            do 175 k = 1, norbs
              do 176 l = k+1, norbs

                ind = (i-1)*2*nocb + (j-1)*2*nosq
     $              + (k-1)*2*norbs + (l-1)*2 + 1 + 16*nosq*nosq

                if ( abs (twom(ind)) .gt. thresh .or.
     $               abs (twom(ind+1)) .gt. thresh ) then
                  write (100, 201) i, 'b', j, 'a', k, 'b', l, 'b',
     $                             twom(ind), twom(ind+1)
                endif
 176          continue
 175        continue
 174      continue
 173    continue

      endif

      write (100, *) '!'

      close (100)


      return
      end


