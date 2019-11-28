C     Anna's Filter - a three-component, frequency-dependent, 
C     data-adaptive polarisation filter.

C     Anna Horleston, December 2005.
C     Hacked to spit out P values. June 2019

C     This is a program to perform a frequency dependent polarisation
C     filter on three component data. The filter is not spatially
C     dependent and can therefore be used to remove (or greatly
C     diminish) un-polarised noise and also to help pick out 
C     polarised phases. The method is based on that of Samson &
C     Olson (1981).

C     Requires the subroutine 'four1' from Numerical Recipes.

C     Basic outline of method:
C     The three components (three sac files) are read into three arrays.
C     A window length of N points is defined (centered on the value 'l')
C     and the data for the first sliding window is fourier transformed.
C     through the time window and for each frequency band
C     a spectral matrix, S, is formed. From the spectral matrix for each
C     frequency window (within the time window) a polarisation parameter
C     P is calculated, for the central frequency in that frequency band,
C      and stored in the array P(N).
C     In this way a polarisation parameter is obtained for each frequency
C     in each time window. The initial data array for that time window is
C     then multiplied by the P(N) (to the power 4) and then the reverse 
C     transfrom is performed and the central value for that time window
C     written out to the sac output file. The time window is then slid
C     through the entire time frame of the three input files.

C     N must be a power of two (so that the fourier transform subroutine
C     works) and should be chosen so that it represents a time longer
C     than the key arrival to be studied.

C     max is the maximum number of points in each data file.

      program fdpfil 
      parameter (max=10000, N=2**7)
      integer g, a
      real arr1(max),arr2(max),arr3(max),SS(6,6), P(N)
      real sac1(max),sac2(max),sac3(max), d(6), v(6,6)
      complex temp1(N),temp2(N),temp3(N),temp(N,3),sumS(3,3)
      complex ptemp1(N),ptemp2(N),ptemp3(N)
      equivalence (temp(1,1),temp1),(temp(1,2),temp2),(temp(1,3),temp3)

C     The following lines read in the sac filnemanes
C     May want to adjust the length of the following depending on 
C     file names to be used.

      character*10 kname
      character*10 knome
      character*10 knime
      character*8 knameo
      character*8 knomeo
      character*8 knimeo
      write(*,*) 'enter sac filenames'
      read(*,*) kname, knome, knime
      write(*,*) 'enter output sac filenames (3)'
      read(*,*) knameo, knomeo, knimeo

C     rsac1 is the subroutine that reads in the sac file.

       call rsac1(kname,arr1,nleni,beg,del,max,nerr)
       if (nerr.ne.0) stop 'invalid filename' 
       if (nleni .gt. max) stop 'too many data points'
       call rsac1(knome,arr2,nlen,beg,del,max,nerr)
       if (nerr.ne.0) stop 'invalid filename' 
       if (nlen .gt. max) stop 'too many data points'
       if (nlen .ne. nleni) stop 'different file size'
       call rsac1(knime,arr3,nlen,beg,del,max,nerr)
       if (nerr.ne.0) stop 'invalid filename' 
       if (nlen .gt. max) stop 'too many data points'
       if (nlen .ne. nleni) stop 'different file size'
       write(0,*) 'Filtering ',nlen,' data points.'


C     Next a temporary array is formed from each original 
C     array to include N points centered around each l value.
C     This sets up the sliding time window.

      do l=(N/2),(nlen-(N/2))
         do i=-(N/2),(N/2)-1
            it =1+(N/2)+i
            ia = l+i+1
            temp1(it) = arr1(ia)
            temp2(it) = arr2(ia)
            temp3(it) = arr3(ia) 
         enddo
         temp1(N) = arr1(ia+1) 
         temp2(N) = arr2(ia+1)
         temp3(N) = arr3(ia+1) 

C     Here is the call to the subroutine four1 which replaces the
C     data in the input array with its fourier transform.
C     four1(data,nn,isign) is the basic fft subroutine. It replaces the 
C     data(1:2*nn) by its discrete fourier transform if isign is +1 or 
C     the inverse transform if isign is -1. data is either a complex 
C     array with nn elements or a real array of length 2nn.
C     nn must be an integer power of 2.

C         write(0,*) 'calling four1'
         call four1(temp1, N, 1)
         call four1(temp2, N, 1)
         call four1(temp3, N, 1)

C     So now temp1, temp2, and temp3 are the fourier transforms of 
C     themselves and have been stored in temp.    
C     The following lines of code set up the sliding frequency window
C     and form the spectral matrix S - or in this case, sumS since S
C     is generated for the 7 frequencies in the window.
      
         do a=1,N-6
            do i=1,3
               do j=1,3
                  sumS(i,j)=0.0
               enddo
            enddo
            do k=a,a+6
               do i=1,3
                   sumS(i,i) = sumS(i,i)+(temp(k,i)*conjg(temp(k,i)))
                   do j=i+1,3
                       sumS(i,j) =sumS(i,j)+(temp(k,i)*conjg(temp(k,j)))
                       sumS(j,i) =conjg(sumS(i,j))
                   enddo
               enddo
            enddo

C     The following code calculates P (as defined in Samson & Olson
C     (1981)).
C     trS is the trace of the spectral matrix, sumS.
C     trS2 is the trace of the spectral matrix squared, sumS^2.

          trS = real(sumS(1,1)+sumS(2,2)+sumS(3,3))
          trS2 = 0.0
              do i=1,3
                 sum = 0.0
                    do j=1,3
                       sum = sum + real(sumS(i,j)*sumS(j,i))
                    enddo
                 trS2 = trS2 + sum
              enddo
          P(a)=(3*trS2-(trS**2))/(2*(trS**2))
      
C     Here we multiply the original temporary array (for this time window)
C     by the P array raised to the power g (g=4 is used here, since it gives
C     best results, but g may range between 2 and 6)
    
         g=4
           ptemp1(a+3)=(temp1(a+3))*(P(a)**g)
            ptemp2(a+3)=(temp2(a+3))*(P(a)**g)
            ptemp3(a+3)=(temp3(a+3))*(P(a)**g)

C     Having calculated all the P(N) we enddo the loop of a values

         enddo

C     The following code sorts out the beginning and ends of the time
C     window since the P values are calculated for the central point of 
C     each frequency band.
C     Here the first 3 points of the time window are multiplied by P(1)
C     and the last 3 by P(N-6) (the P value for the last frequency band.
 
         do i=1,3
            ptemp1(i)=(temp1(i))*(P(1)**g)
            ptemp2(i)=(temp2(i))*(P(1)**g)
            ptemp3(i)=(temp3(i))*(P(1)**g)
         enddo
         do i=N-2,N
            ptemp1(i)=(temp1(i))*(P(N-6)**g)
            ptemp2(i)=(temp2(i))*(P(N-6)**g)
            ptemp3(i)=(temp3(i))*(P(N-6)**g)
         enddo

C     And now we do the inverse transform by calling four1 with isign=-1.
         call four1(ptemp1, N, -1)
         call four1(ptemp2, N, -1)
         call four1(ptemp3, N, -1)

C     And now we write out the value of each ptemp array for element
C     'l' - the central value for the time window - to the arrays sac1,
C     sac2, and sac3.

         sac1(l)=ptemp1(N/2)
         sac2(l)=ptemp2(N/2)
         sac3(l)=ptemp3(N/2)

C     After the time window has been slid through the entire data file
C     we enddo the loop of l's

      enddo

C     And finally we write out the sac1, sac2 and sac3 arrays to the
C     sac output files specified at the beginning of the program.

      write(*,*) 'calling wsac1'
      call wsac1(knameo,sac1,nlen,beg,del,nerr)
      call wsac1(knomeo,sac2,nlen,beg,del,nerr)
      call wsac1(knimeo,sac3,nlen,beg,del,nerr)

      end
