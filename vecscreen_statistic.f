	program statistic
	implicit none
	integer var
	parameter (var=100000)
	character line*10000, char_n_pid*10,arg(4)*100,name_file*100
	character n_char_bins*2,lline*100,Chw*4,Chm*8,Chs*6,Chso*7
	character time*100,charline(var)*1000,lineb*1000,wms(5)*1
	character wmscharline(var)*1
c	data charline /var*'  '/
	integer i,ir,ik,k,j,l,dim_bar,ip,ipl
	integer n_pid,Iw,Im,Is,Iso,pIlen,ipp
	real, allocatable :: Nw(:,:),Nm(:,:),Ns(:,:),Nso(:,:)
	real range_bar,bins,sT,so,sf,a,b,Icharline(var),TIcharline
	real wIcharline(var),mIcharline(var),sIcharline(var)
	real wTIcharline,mTIcharline,sTIcharline,Ic,Iv,Ii
c	data Icharline /var*0/ 
	real Vw(3),Vm(3),Vs(3),Vso(3),pp,pIw,pIm,pIs,pIso
	real :: start, finish

c	Icharline(1:var)=(/(0,i=1,var)/)
	
	do i=1,var
	 charline(i)(1:2)='  '
	 Icharline(i)=0.0 
	 wmscharline(i)=' ' 
	enddo
	
	do i=1,5
	 wms(i)=' '
	enddo

	 if(iargc().lt.2.or.iargc().eq.3)then
         call getarg(1,arg(1))
         if(arg(1).eq.'-H'.or.arg(1).eq.'-h'.or.arg(1).eq.'--help')then
         write(*,*)
	 write(*,'(a53)')'Usage: statistic_vecscreen.x -ivf2 file.txt -bin
     + i'
	 write(*,*)
	 write(*,'(T5,a54,/,T8,a33,/)')'-vf2 Output file obtained from the
     + vecscreen execution',' with the option -f defined as 2'
	 write(*,'(T5,a42,/)')'-bin Intervals of histogram (defaults=0.5)'
         stop
         elseif(arg(1).eq.'--version')then	
	 write(*,'(/,T5,a27,2/,T5,a30,/,T3,a27,/)')'Version: 0.1 (Nov 19,
     +2015)','Author: J. Rodriguez-Salarichs','e-mail: javierr@cib.csic.e
     +s'
         stop
	 else
	 write(*,'(/,T5,a24)')'ERROR: Missing arguments'
         write(*,*)
	 write(*,'(a53)')'Usage: statistic_vecscreen.x -ivf2 file.txt -bin
     + i'
	 write(*,*)
	 write(*,'(T5,a54,/,T8,a33,/)')'-vf2 Output file obtained from the
     + vecscreen execution',' with the option -f defined as 2'
	 write(*,'(T5,a42,/)')'-bin Intervals of histogram (defaults=0.5)'
         stop

	 endif
	 elseif(iargc().eq.2)then
	  call getarg(2,arg(2))
	  name_file=arg(2)
	  bins=0.5
	 elseif(iargc().ge.4)then
	  call getarg(2,arg(2))
	  call getarg(4,arg(4))
	  name_file=arg(2)
	  read(arg(4),*)bins
         endif

        call cpu_time(start)
	
	dim_bar=100.0/bins


	n_pid = getpid()
	write(char_n_pid,'(i10)')n_pid
	
	char_n_pid=ADJUSTL(char_n_pid)
	name_file=ADJUSTL(name_file)

	allocate(Nw(dim_bar,dim_bar),Nm(dim_bar,dim_bar))
	allocate(Ns(dim_bar,dim_bar),Nso(dim_bar,dim_bar))
	
	Nw(1,1)=bins
	Nw(2,1)=0	
	Nm(1,1)=bins
	Nm(2,1)=0	
	Ns(1,1)=bins
	Ns(2,1)=0	
	Nso(1,1)=bins
	Nso(2,1)=0	

	do ik=2,dim_bar
	 Nw(1,ik)=Nw(1,ik-1)+bins
	 Nw(2,ik)=0
	 Nm(1,ik)=Nm(1,ik-1)+bins
	 Nm(2,ik)=0
	 Ns(1,ik)=Ns(1,ik-1)+bins
	 Ns(2,ik)=0
	 Nso(1,ik)=Nso(1,ik-1)+bins
	 Nso(2,ik)=0
	enddo



        open(10,file='statistics_vecscreen_'//char_n_pid(1:len_trim(char
     +_n_pid))//'.log')
	
	write(10,*)'==Program extracts a valuable information from an outp
     +ut vecscreen file=='

	write(10,'(/)')

	write(10,*)'The sequence regions are cataloged in weak, moderate,'
	write(10,*)'strong and suspect origin according to their likelihoo
     +d'
	write(10,*)'to be a vector'

	write(10,'(/)')

	write(10,'(a57)')'================================================
     +========='
	write(10,*)
	write(10,'(a55)')'Strong Match to Vector: Strong matches usually i 
     +ndicate' 
	write(10,'(a53)')'that the segment originated from foreign DNA (ve
     +ctor,'
	write(10,'(a52)')'adapter, linker, or primer) that was attached to 
     + the'
	write(10,'(a51)')'source DNA/RNA during the cloning process.(Expec
     +t 1'
	write(10,'(a52,/)')'random match in 1,000,000 queries of length 35 
     +0 kb.)'

	write(10,'(a57)')'------------------------------------------------
     +---------'
	write(10,'(a57)')'Moderate Match to Vector: Strong matches usually 
     + indicate' 
	write(10,'(a53)')'that the segment originated from foreign DNA (ve
     +ctor,'
	write(10,'(a52)')'adapter, linker, or primer) that was attached to 
     + the'
	write(10,'(a51)')'source DNA/RNA during the cloning process.(Expec
     +t 1'
	write(10,'(a48,/)')'random match in 1,000 queries of length 350 kb 
     +.)'

	write(10,'(a57)')'------------------------------------------------
     +---------'
	write(10,'(a55)')'Weak Match to Vector: Although weak matches ofte
     +n occur'
	write(10,'(a56)')'by chance, they indicate foreign sequence whenev
     +er there'
	write(10,'(a52)')'is corroborating evidence of contamination.(Expe
     +ct 1'
	write(10,'(a45,/)')'random match in 40 queries of length 350 kb.)'

	write(10,'(a57)')'------------------------------------------------
     +---------'

	write(10,'(a52)')'Segment of Suspect Origin: Any segment of fewer 
     +than'
	write(10,'(a54)')'50 bases between two vector matches or between a
     + match'
	write(10,'(a11,/)')'and an end.'

	write(10,'(a57)')'================================================
     +========='

	write(10,'(2/)')
	write(10,*)'Input file: '//name_file(1:len_trim(name_file))

	write(10,'(2/,a55,2/)')'The program is running. Please, wait few m
     +inutes.'

	open(1,file=name_file(1:len_trim(name_file)),status='old')
        open(2,file='.vectors_from_text_'//char_n_pid(1:len_trim(char_n_
     +pid))//'.txt',status='unknown')
        i=1
	ir=0
	

	write(10,*)'=> First step: The vectors are being identified...'

	do while(ir.eq.0)
	 read(1,'(a1000)',iostat=ir)line   
         if(line(1:7).eq.">Vector")then
          if(i.ge.2)then
	   if(i-j.gt.41)then
c	   if(i-j.gt.3)then !en el caso de Crugosa_S1_L001_R1_001_02_vecscreen.txt
c	    write(*,*)line(:len_trim(line))	    
	    write(2,*)j,i-1
	   endif
	  endif 
          j=i 
         endif 
         i=i+1
	enddo

	rewind (2)
	rewind (1)

c	open(3,file='.w_vector_'//char_n_pid(1:len_trim(char_n_pid))//
c     +'.txt',status='new')
c        open(4,file='.m_vector_'//char_n_pid(1:len_trim(char_n_pid))//
c     +'.txt',status='new')
c        open(5,file='.s_vector_'//char_n_pid(1:len_trim(char_n_pid))//
c     +'.txt',status='new')
c        open(6,file='.so_vector_'//char_n_pid(1:len_trim(char_n_pid))//
c     +'.txt',status='new')
c	open(7,file='.ls_vector_'//char_n_pid(1:len_trim(char_n_pid))//
c     +'.txt',status='new')
c        open(8,file='.t_vector_'//char_n_pid(1:len_trim(char_n_pid))//
c     +'.txt',status='new')

c	close(3,STATUS='delete')
c	close(4,STATUS='delete')	
c	close(5,STATUS='delete')
c	close(6,STATUS='delete')
c	close(7,STATUS='delete')



	ir=0
	l=0
	write(10,*)
	write(10,*)'=> Second step: The regions of sequences are being cat
     +aloged'
	write(10,*)'in weak, moderate, strong and suspect ...'

	pp=0
	pIw=0
	pIm=0
	pIs=0
	pIso=0
	pIlen=0
	
	do while(ir.eq.0)
	 
	 read(2,*,iostat=ir)j,i
	 do k=1+l,j-1
	  read(1,*)
	 enddo
	  Iw=0
	  Im=0
	  Is=0
	  Iso=0
	  ipp=1
	  do ip=1,5
	   wms(ip)=' '
	  enddo
	 do k=j,i
	  read(1,'(a1000)')line

	  if(line(1:4).eq."Weak")then
	   read(1,*)a,b
	   Vw(1)=a
	   Vw(2)=b
	   Iw=1
	   pIw=pIw+(Vw(2)-Vw(1)+1)
	   backspace (1)	
	   do ip=1,5
	    if(wms(ip).eq.' ')then
	     wms(ip)='w'
	     goto 11
	    endif
	   enddo
  11       continue
   
	  elseif(line(1:8).eq."Moderate")then
	   read(1,*)a,b
	   Vm(1)=a
	   Vm(2)=b
	   Im=1
  	   pIm=pIm+(Vm(2)-Vm(1)+1)
	   backspace (1)
	   do ip=1,5
	    if(wms(ip).eq.' ')then
	     wms(ip)='m'
	     goto 12
	    endif
	   enddo
  12       continue   
	  elseif(line(1:6).eq."Strong")then
	   read(1,*)a,b
	   Vs(1)=a
	   Vs(2)=b
	   Is=1
	   pIs=pIs+(Vs(2)-Vs(1)+1)

	   backspace (1)
	   do ip=1,5
	    if(wms(ip).eq.' ')then
	     wms(ip)='s'
	     goto 13
	    endif
	   enddo
  13       continue   
	  elseif(line(1:7).eq."Suspect")then
	   read(1,*)a,b
	   Vso(1)=a
	   Vso(2)=b
	   Iso=1
	   pIso=pIso+(Vso(2)-Vso(1)+1)

	   backspace (1)
	  elseif(line(1:16).eq."length of query:")then
	   line=ADJUSTL(line(17:))
	   read(line(:len_trim(line)),*)sT
	   pIlen=pIlen+sT
	   pp=pp+1
	   if(Iw.eq.1)then
	    do ik=1,dim_bar
	     if(Nw(1,ik).ge.(Vw(1)/sT)*100.and.Nw(1,ik).le.(Vw(2)/sT)*100)
     +then	      
	      Nw(2,ik)=Nw(2,ik)+1.0
	     endif
	    enddo
	    Iw=0
	   endif

	   if(Im.eq.1)then
	    do ik=1,dim_bar
	     if(Nm(1,ik).ge.(Vm(1)/sT)*100.and.Nm(1,ik).le.(Vm(2)/sT)*100)
     +then	      
	      Nm(2,ik)=Nm(2,ik)+1.0
	     endif
	    enddo
	    Im=0
	   endif

	   if(Is.eq.1)then
	    do ik=1,dim_bar
	     if(Ns(1,ik).ge.(Vs(1)/sT)*100.and.Ns(1,ik).le.(Vs(2)/sT)*100)
     +then	      
	      Ns(2,ik)=Ns(2,ik)+1.0
	     endif
	    enddo
	    Is=0
	   endif

	   if(Iso.eq.1)then
	    do ik=1,dim_bar
	     if(Nso(1,ik).ge.(Vso(1)/sT)*100.and.Nso(1,ik).le.(Vso(2)/sT)*
     +100)then	      
	      Nso(2,ik)=Nso(2,ik)+1.0
	     endif
	    enddo
	    Iso=0	    
	   endif



	  elseif(line(1:7).eq.">gnl|uv")then
	   line=ADJUSTL(line(1:))
   

	   do ik=1,10
	    read(1,'(a1000)')lineb
	    lineb=ADJUSTL(lineb(1:))
	    if(lineb(1:8).eq.'Length =')then
	     do ip=1,ik
	      backspace (1)
	     enddo
	     goto 30
	    else
	     line=line(1:len_trim(line))//' '//lineb(1:len_trim(lineb))
	     line=ADJUSTL(line(1:))
	    endif
	   enddo 
  30       continue
 
	   ipl=50
	   do ik=1,var
	    do ip=1,50
	     if(line(ip:ip).eq.':')then
	      ipl=ip-1	      
	      goto 20
	     endif
	    enddo
  20	    continue
	    if(charline(ik)(1:ipl).eq.line(1:ipl))then
	     Icharline(ik)=Icharline(ik)+1.0
	     do ip=ipp,5
	      if(wms(ip).eq.' ')then
	       goto 14
	      else
	       wmscharline(ik)=wms(ip)
	       ipp=ip+1
	       if(wms(ip).eq.'w')then
	        wIcharline(ik)=wIcharline(ik)+1.0
	       elseif(wms(ip).eq.'m')then
	        mIcharline(ik)=mIcharline(ik)+1.0
	       elseif(wms(ip).eq.'s')then
	        sIcharline(ik)=sIcharline(ik)+1.0
	       endif 
	       goto 14
	      endif
	     enddo
  14         continue

	     goto 10
	    elseif(charline(ik)(1:2).eq.'  ')then
	     charline(ik)=line(:len_trim(line))
	     charline(ik)=ADJUSTL(charline(ik)(1:))
	     goto 10
	    endif
	   enddo
  10       continue
c	   write(8,*)j,line(:len_trim(line))



	  endif
	 enddo
	 l=i
	enddo


	close(1)
	close(2,STATUS='delete')
c	close(8,STATUS='delete')

	write(10,*)
	write(10,*)'=> Third step: Results ...'
	write(10,*)

	open(11,file='bar_plot_vector_'//char_n_pid(1:len_trim(char_n_pid)
     +)//'.txt',status='new')
	do i=1,dim_bar
	 write(11,*)Nw(1,i),(Nw(2,i)/pp)*100,(Nm(2,i)/pp)*100,(Ns(2,i)/pp
     +)*100,(Nso(2,i)/pp)*100,((Nm(2,i)+Ns(2,i))/pp)*100
	enddo
	close(11)

	write(10,*)'Result 1: Relative abundance of contaminant along the 
     +DNA sequences'
	write(10,*)'and relative postion in the sequences was saved in ...
     +'

	write(10,*)
	write(10,*)'File: bar_plot_vector_'//char_n_pid(1:len_trim(char_n_
     +pid))//'.txt'
	write(10,*)
	write(10,*)'column 1: relative total presence of contaminant'
	write(10,*)'column 2-6: relative postion in the sequences'
	write(10,*)'column 2: Weak Matches'
	write(10,*)'column 3: Moderate Matches'	
	write(10,*)'column 4: Strong Matches'
	write(10,*)'column 5: Suspect Origin'
	write(10,*)'column 6: Strong and Moderate Matches'


	write(10,*)
	write(10,*)
	write(10,*)'Result 2: Contaminant (%)'
	write(10,*)
	write(10,'(2x,a12,3x,a16,2x,a14,3x,a14,5x,a10)')'Weak Matches','Mo
     +derate Matches','Strong Matches','Suspect Origin','No Matches'
	write(10,*)(pIw/pIlen)*100,(pIm/pIlen)*100,(pIs/pIlen)*100,(pIso/p
     +Ilen)*100,((pIlen-pIw-pIm-pIs-pIso)/pIlen)*100

	
	write(10,*)
	write(10,*)
	write(10,*)'Result 3: Types of contaminants'
	write(10,*)
	write(10,*)'The contaminant types were saved in the w_vector.txt,'
	write(10,*)'m_vector and s_vector following their matches.'
	write(10,*)'The T_vector.txt contains all contaminant types'
	write(10,*)
	write(10,*)'     number match  percent          type'
	write(10,*)'----------------------------------------'
	write(10,*)

	TIcharline=0
	wTIcharline=0
	mTIcharline=0
	sTIcharline=0
	do ik=1,var
	 TIcharline=TIcharline+Icharline(ik)
	 wTIcharline=wTIcharline+wIcharline(ik)
	 mTIcharline=mTIcharline+mIcharline(ik)
	 sTIcharline=sTIcharline+sIcharline(ik)
	enddo

	

	open(3,file='w_vector_'//char_n_pid(1:len_trim(char_n_pid))//
     +'.txt',status='new')
        open(4,file='m_vector_'//char_n_pid(1:len_trim(char_n_pid))//
     +'.txt',status='new')
        open(5,file='s_vector_'//char_n_pid(1:len_trim(char_n_pid))//
     +'.txt',status='new')
        open(6,file='T_vector_'//char_n_pid(1:len_trim(char_n_pid))//
     +'.txt',status='new')





	do ik=1,var
	 if(Icharline(ik).ne.0)then
	    do ip=1,1000
	     if(charline(ik)(ip:ip).eq.' ')then
	      ipl=ip+1	      
	      goto 50
	     endif
	    enddo
  50	    continue

	 if(wmscharline(ik).eq.'w')then
	  write(3,*)(wIcharline(ik)/wTIcharline)*100,'"',charline(ik)(ipl:
     +len_trim(charline(ik))),'"'
	 elseif(wmscharline(ik).eq.'m')then
	  write(4,*)(mIcharline(ik)/mTIcharline)*100,'"',charline(ik)(ipl:
     +len_trim(charline(ik))),'"'
	 elseif(wmscharline(ik).eq.'s')then
	  write(5,*)(sIcharline(ik)/sTIcharline)*100,'"',charline(ik)(ipl:
     +len_trim(charline(ik))),'"'

	 endif
	  write(6,*)(Icharline(ik)/TIcharline)*100,'"',charline(ik)(ipl:le
     +n_trim(charline(ik))),'"'
	  write(10,*)ik,wmscharline(ik),' ',(Icharline(ik)/TIcharline)*100
     +,charline(ik)(ipl:len_trim(charline(ik)))
	 endif
	enddo	

	ir=0
	Ic=0
	Ii=0
	Iv=0
	rewind (3)
	do while(ir.eq.0)
	 read(3,*,iostat=ir)a,line
	 line=adjustl(line)
	 do j=1,len_trim(line)-12
	  if(line(j:j+12).eq.'loning vector')then
	   Ic=Ic+a
	   goto 45
	  endif
	 enddo
	 do j=1,len_trim(line)-6
	  if(line(j:j+6).eq.'llumina')then
	   Ii=Ii+a
	   goto 45
	  endif
	 enddo
	 do j=1,len_trim(line)-5
	  if(line(j:j+5).eq.'vector')then
	   Iv=Iv+a
	   goto 45
	  endif
	 enddo
  45     continue 

	enddo
	
	write(10,*)
	write(10,*)'Illumina''s_sequences Cloning_Vectors   Other_vectors
     +Other_contaminants'
	write(10,*)'------------------------------------------------------
     +-------------------'
	write(10,*)
	write(10,*)'Weak Matches ...'
	write(10,*)'  ',Ii,Ic,Iv,100-(Ic+Ii+Iv)

	ir=0
	Ic=0
	Ii=0
	Iv=0
	rewind (4)
	do while(ir.eq.0)
	 read(4,*,iostat=ir)a,line
	 line=adjustl(line)
	 do j=1,len_trim(line)-12
	  if(line(j:j+12).eq.'loning vector')then
	   Ic=Ic+a
	   goto 46
	  endif
	 enddo
	 do j=1,len_trim(line)-6
	  if(line(j:j+6).eq.'llumina')then
	   Ii=Ii+a
	   goto 46
	  endif
	 enddo
	 do j=1,len_trim(line)-5
	  if(line(j:j+5).eq.'vector')then
	   Iv=Iv+a
	   goto 46
	  endif
	 enddo
  46     continue 

	enddo

	write(10,*)
	write(10,*)'Moderate Matches ...'
	write(10,*)'  ',Ii,Ic,Iv,100-(Ic+Ii+Iv)


	ir=0
	Ic=0
	Ii=0
	Iv=0
	rewind (5)
	do while(ir.eq.0)
	 read(5,*,iostat=ir)a,line
	 line=adjustl(line)
	 do j=1,len_trim(line)-12
	  if(line(j:j+12).eq.'loning vector')then
	   Ic=Ic+a
	   goto 47
	  endif
	 enddo
	 do j=1,len_trim(line)-6
	  if(line(j:j+6).eq.'llumina')then
	   Ii=Ii+a
	   goto 47
	  endif
	 enddo
	 do j=1,len_trim(line)-5
	  if(line(j:j+5).eq.'vector')then
	   Iv=Iv+a
	   goto 47
	  endif
	 enddo
  47     continue 

	enddo

	write(10,*)
	write(10,*)'Strong Matches ...'
	write(10,*)'  ',Ii,Ic,Iv,100-(Ic+Ii+Iv)


	close(3)
	close(4)	
	close(5)
	close(6)

        call cpu_time(finish)
	write(time,*)finish-start
	time=ADJUSTL(time)
	write(10,'(4/)')
	write(10,*)'(Time elapsed = ',time(1:len_trim(time)),' seconds.)'
	write(10,'(2/,a9,/)')' ==Done=='


	close(10)
	end
