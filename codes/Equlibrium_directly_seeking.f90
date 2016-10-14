!------------------------------------------
!             ƽ��̬��������
!------------------------------------------
!   Version      : V1.0
!   Coded by     : Xiashaopeng
!   Date         : 2016-8-15
!------------------------------------------
!   Description  : ֱ���������ζѵ�ƽ��̬
!                  ����������ϡ�����Ҫ��
!------------------------------------------
!   Code Flow    :
!          1.   ��ʼ�ٽ���㣬�ṩȼ�ļ���
!               ���跴Ӧ����Ϣ
!          2.   ȼ�ļ��㣬ֱ������Ũ�Ȳ���
!          3.   ��������ʴ��������飬����
!               �������Լ���Ҫ����Ũ��
!          4.   �ж����ε������У���Ҫ����
!               ��Ũ���Ƿ���������Ҫ����
!               �����㣬����벽�� 5�����
!               ������Ҫ���򽫼���õ���
!               �µ������ʴ��벽�� 2�У���
!               �½���ȼ�ļ��㡣
!          5.   ���µĺ���Ũ�Ƚ����ٽ���㣬
!               ��������������ٽ�Ҫ��
!               ��ɹ���������ƽ��̬�����
!               �������ٽ�Ҫ���򽫸ò���
!               �Ľ��桢ͨ������Ϣ���벽��
!               2�����½���ȼ�ļ��㡣
!-----------------------------------------


module INITIAL_CSAS
!------------------------------------------module coment
!   Version      : V1.0
!   Coded by     : Xiashaopeng
!   Date         : 2016-8-15
!------------------------------------------
!   Description  : ���������еĵ�һ��
!            Ϊȼ�ļ����ṩ��������
!------------------------------------------
!   Contains     :
!       1.  Initial_Keff ��ʼ�ٽ���� ���һ��Ǹ�module�Ľӿ�
!       2.  sumhm10      ��������ͳ��,
!       3.  molten10     �ṩscale�����������Ũ��
!------------------------------------------

contains

	subroutine Initial_Keff
!------------------------------------------
!
	character(50)rd
	character(7)nucly

		call sumhm10
!		call sumhm20

		sum1=0.0
		sum2=0.0
		open(1,file='M_sum_HM1.dat')
		read(1,*)
		read(1,*)sum1
		close(1)

!		open(1,file='M_sum_HM2.dat')
!		read(1,*)
!		read(1,*)sum2
!		close(1)

		sumhm0=sum1+sum2

		open(1,file='M_sum_HM.dat')
		write(1,*)"�ؽ���������"
		write(1,'(2x,E15.8)')sumhm0
		close(1)

		call molten_10
!		call molten_20

	call new_CSAS 
!	stop

!	call system("E:\scale6.1\cmds\runscale CSAS.inp")
	call read_keff(ncheck)
	if(ncheck.eq.0)then
		call system("E:\scale6.1\cmds\runscale CSAS.inp")
		call read_keff(ncheck)
		if(ncheck.eq.0)stop
    endif

!    call read_kmt
    
	open(1,file='input.inp')
	read(1,*)  !
	read(1,*)v1  !���
	read(1,*)  !
	read(1,*)v2  !���
	read(1,*)  
	read(1,*)fkeffdown,fkeffup    !keff_down and keff_up
	close(1)

	open(2,file='keff0.dat')
	read(2,*)fkeff
	close(2)

	if(fkeff.gt.fkeffup .or. fkeff.lt.fkeffdown)then
		write(*,*)"change the molten ",fkeff
!		stop
	endif

 100	format(8a)
    end subroutine Initial_Keff
   subroutine read_kmt
        character(200)rd
		character(5)nuclx
        open(1,file='TRITON.out')
        open(2,file='KMT_act_temp.dat')
        do while(.true.)
            read(1,100,iostat=nn)rd
            if(nn/=0)exit
            if(rd(1:9).eq." Activity")then               
				nuclx=rd(28:32)
                    do while(.true.)
                        read(1,100)rd
                        if(rd(14:18).eq."total")then
                        write(2,*)nuclx,rd(23:35)
                        exit
                        endif
                    enddo
            endif
        enddo
        close(1)
        close(2)
		
		open(1,file='KMT_act_temp.dat')
		open(2,file='KMT_act.dat')
		write(2,*)"��������     ��Ӧ��"
		do while(.true.)
			read(1,*,iostat=nn)nza,fac
            if(nn/=0)exit
			if(nza.eq.52601)then
				write(2,*)521271,fac
			else if(nza.eq.52611)then
				write(2,*)521291,fac
			else if(nza.eq.27601)then
				write(2,*)270581,fac
			else if(nza.eq.47601)then
				write(2,*)471101,fac
			else if(nza.eq.61601)then
				write(2,*)611481,fac
			else if(nza.eq.66601)then
				write(2,*)671661,fac
			else if(nza.eq.95601)then
				write(2,*)952421,fac
			else if(nza.eq.95611)then
				write(2,*)952441,fac
			else if(nza.eq.48601)then
				write(2,*)481151,fac
			else
				write(2,*)nza*10,fac
			endif
		end do
		close(1)
		close(2)
100 format(a)
    end subroutine read_kmt
    
    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine molten_10
	character(60) rd
	character(8) nucl,nucl2
	character(7) nuclx,nucly
	character(2) nucla
	character(1) tab


	dimension nac(50),fac1(50),fac2(50),fmac(50)
	dimension fwac(50),fwnac(50),wmac(50)

	nac=0
	fac1=0.0
	fac2=0.0
	fmac=0.0
	fwac=0.0
	fwnac=0.0
	wmac=0.0
      xu = 0.0
      xpu = 0.0
      xma = 0.0


	open(1,file='molten-10.inp')
	read(1,*)
	read(1,*)v0  !���������cm^3����

	read(1,*)
	read(1,*)fp0  !�����ܶ�(g/cm^3):

	read(1,*)
	read(1,*)ntemp  

	read(1,*)
	read(1,*)wpower

!c	!LiF-BeF2-ThF4-UF4-PuF4-MAF4��Ħ��Ũ�������
	i=1
	sumac=0.0
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit
			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy
			fmac(i)=fa*fy
			if(nza.eq.92000)xu=fy
			if(nza.eq.94000)xpu=fy
			if(nza.eq.95000)xma=fy
			sumac=sumac+fmac(i)
!c	write(*,*)nza,fa,fy,fmac(i),sumac

		i=i+1
		end do
		close(1)
	number=i-1

	close(1)

	open(1,file='power.dat')
	write(1,*)" ���ʣ�"
	write(1,*)wpower
	close(1)


!c	U��ʼ����
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c	i=1	

	numu233=1
	number=i-1
	sumacx=0.0

	if(xu.ne.0.0)then
		open(1,file='Proportion_U.inp')
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit

			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy*xu
			fmac(i)=fa*fy*xu
			sumac=sumac+fmac(i)
			sumacx=sumacx+fmac(i)
			if(nza.eq.92233)numu233=i
!c	write(*,*)nza,fa,fy,fmac(i),xu,sumac

		i=i+1
		end do
		close(1)
	number=i-1
	endif
	
!c	write(*,*)numu233
!c	Pu��ʼ����
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	if(xpu.ne.0.0)then
		open(1,file='Proportion_Pu.inp')
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit

			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy*xpu
			fmac(i)=fa*fy*xpu
			sumac=sumac+fmac(i)
			sumacx=sumacx+fmac(i)
!c	write(*,*)nza,fa,fy,fmac(i),xpu,sumac
		i=i+1
		end do
		close(1)
	number=i-1
	endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!c	MA��ʼ����
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	if(xma.ne.0.0)then
		open(1,file='Proportion_MA.inp')
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit

			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy*xma
			fmac(i)=fa*fy*xma
			sumac=sumac+fmac(i)
			sumacx=sumacx+fmac(i)
!c	write(*,*)nza,fa,fy,fmac(i),xma,sumac
		i=i+1
		end do
		close(1)
	number=i-1
	endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c	write(*,*)numu233,number

!c	����������
	summ0=v0*fp0

	sum=sumac

!c	��������ϵ��
	rat=summ0/sumac
!c	write(*,*)summ0,sum,rat,sumacx


	ts=0.0
	do j=1,number
		fwac(j)=fmac(j)*rat
		if(fac1(j).ge.220)then
			ts=ts+fwac(j)
!c			write(*,*)fac1(j),fwac(j),ts
		endif
	enddo




	open(1,file='M_sum_HM.dat')
	read(1,*)
	read(1,*)sumhm0
	close(1)

!c	�ؽ�����һ��������ص�����


	st=0.0
	do j=1,number
	wmac(j)=fwac(j)/sumhm0*1.0e6
	st=st+wmac(j)
!c	write(*,*)nac(j),	wmac(j)
	enddo	

!c	�����ص�ԭ���ܶ�
	
	st=0.0
	do j=1,number
	fwnac(j)=fwac(j)/v0/fac1(j)*0.602214199
	enddo
!c	stop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	open(2,file='com.inp')


	do i=1,number
		nz=int(nac(i)/1000)
		na=nac(i)-nz*1000

		open(1,file='elment.inp')
		do while(.true.)
			read(1,1001,iostat=nn)nz0,nulca
			if(nn/=0)exit
			if(nz0.eq.nz)exit
		end do
		close(1)

		fw=fwnac(i)
		if(fw.ne.0.0 .and. na.ne.0)then

			if(na.le.9)then
				write(2,1002)nulca,"-",na,fw
			elseif(na.le.99)then
				write(2,1003)nulca,"-",na,fw
			else
				write(2,1004)nulca,"-",na,fw
			endif
		endif

	enddo
	close(2)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	open(8,file='com_1_old.dat')
	write(8,100)" nucl nza  na con(atoms) con_real(grams)  con(grams) "
	write(8,*)

    open(9,file='com_FLiBe.dat')
	open(3,file='posNC1.inp')
	open(4,file='nucl_ID.inp')
	open(5,file='com_1.inp')

	j=0
	read(4,*)
	read(4,*)
    jj=0!--���ڽض�com_1.inp�к�������
	do while(.true.)
		read(4,*,iostat=nn)nuclx,nza,fna,ntype
		if(nn/=0)exit

        if(nza.eq.892240)jj=1
		fatoms=0.0

		do i=1,number
			nzax=nac(i)*10
			if(nzax.eq.nza)then
			fatoms=fwnac(i)
!c			write(*,*)fw
			exit
			endif		
		enddo	
!c	write(*,*)sumhm0

		fw=fatoms*v0*fna/0.602214199
		fq=fw/sumhm0*(1.0e6)
		fs=fq/fna
!c		fatoms=fq/v0/na*0.602214199

        if(jj.eq.0)then
		if(fatoms.ge.1e-20)then
			if(nza.eq.60120)then
				write(5,101)"c      ",' 1 0 ',fatoms,ntemp," end   "
			else
				write(5,101)nuclx,' 1 0 ',fatoms,ntemp," end     "
			endif
        endif
        endif

	write(8,"(a,(2x,i6),(2x,f9.5),3(2x,E15.8))")nuclx,nza,fna,fatoms,fw,fq
    if((nza.eq.30060).or.(nza.eq.30070).or.(nza.eq.40090).or.(nza.eq.90190))then
        write(9,"(a,(2x,i6),(2x,f9.5),3(2x,E15.8))")nuclx,nza,fna,fatoms,fw,fq
    endif
    
		
		if(fw.le.1e-20)fw=1e-20
	
		if(j.eq.0)then
		if(fw.le.1e-20)fq=1.01e-20
			write(3,*)" 74**   ",fq
		else
		if(fw.le.1e-20)fq=1.01e-20
			write(3,*)fq
		endif
	

	j=1
	enddo	
	close(3)
	close(4)
	close(8)
	close(5)
    close(9)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  ͳ��ȼ��ǰ��TRU��mol%   cccccccccccccccccccccccccccccc
      open(8,file='com_1_old.dat')
          read(8,*)	
          read(8,*)
          sumtru0=0.0
      	do while(.true.)
			read(8,*,iostat=nn)nucly,nza,fna,fatoms,fw,fq,fs
              if(nn.ne.0)exit
              if(nza.ge.930000)sumtru0=sumtru0+fatoms !!!!!!!ͳ��TRU���ص�mol�ٷֱ�
          enddo
	close(8)     
      
	open(1,file='M_sum_TRU1.dat')
	write(1,*)"TRU��mol%��"
	write(1,'(2x,E15.8)')sumtru0*v0
	close(1)
!
!
	open(3,file='com_31.inp')
	open(2,file='com_10.inp')
	open(1,file='com_1.inp')
	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
		write(2,100)rd
		write(3,100)rd
	enddo
	close(1)
	close(2)
	close(3)

!	stop	

101	format(2a,2x,E15.8,i6,a)
100	format(a)
200	format(4a,i6)
300	format(a,E15.8,2x,E15.8)
400	format(4x,i3,a)
500	format(4x,i3,3a)
600	format(a,2x,i3,a)
700	format(a,3(2x,i3,a))
800	format(a,2x,i3,5(2x,E15.8))
900	format(a,2x,a,4(2x,E15.8))
1000	format(2x,i8,4(2x,E15.8))
1001	format(i7,2x,a)
1002	format(a2,a,i1,2x,E15.8)
1003	format(a2,a,i2,2x,E15.8)
1004	format(a2,a,i3,2x,E15.8)
    end subroutine molten_10

    
    
	subroutine sumhm10
	character(60) rd
	character(8) nucl,nucl2
	character(12) nuclx,nucly

	dimension nac(50),fac1(50),fac2(50),fmac(50)
	dimension fwac(50),fwnac(50),wmac(50)

	nac=0
	fac1=0.0
	fac2=0.0
	fmac=0.0
	fwac=0.0
	fwnac=0.0
	wmac=0.0
      xu=0.0
      xpu=0.0
      xma=0.0


	open(1,file='molten-10.inp')
	read(1,*)
	read(1,*)v0  !���������cm^3����

	read(1,*)
	read(1,*)fp0  !�����ܶ�(g/cm^3):

	read(1,*)
	read(1,*)ntemp  

	read(1,*)
	read(1,*)wpower

!c	!LiF-BeF2-ThF4-UF4-PuF4-MAF4��Ħ��Ũ�������
	i=1
	sumac=0.0
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit
			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy
			fmac(i)=fa*fy
			if(nza.eq.92000)xu=fy
			if(nza.eq.94000)xpu=fy
			if(nza.eq.95000)xma=fy
			sumac=sumac+fmac(i)
!c	write(*,*)nza,fa,fy,fmac(i),sumac

		i=i+1
		end do
		close(1)
	number=i-1

	close(1)

	open(1,file='power.dat')
	write(1,*)" ���ʣ�"
	write(1,*)wpower
	close(1)


!c	U��ʼ����
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c	i=1	

	numu233=1
	number=i-1
	sumacx=0.0

	if(xu.ne.0.0)then
		open(1,file='Proportion_U.inp')
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit

			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy*xu
			fmac(i)=fa*fy*xu
			sumac=sumac+fmac(i)
			sumacx=sumacx+fmac(i)
			if(nza.eq.92233)numu233=i
!c	write(*,*)nza,fa,fy,fmac(i),xu,sumac

		i=i+1
		end do
		close(1)
	number=i-1
	endif
	
!c	write(*,*)numu233
!c	Pu��ʼ����
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	if(xpu.ne.0.0)then
		open(1,file='Proportion_Pu.inp')
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit

			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy*xpu
			fmac(i)=fa*fy*xpu
			sumac=sumac+fmac(i)
			sumacx=sumacx+fmac(i)
!c	write(*,*)nza,fa,fy,fmac(i),xpu,sumac
		i=i+1
		end do
		close(1)
	number=i-1
	endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!c	MA��ʼ����
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	if(xma.ne.0.0)then
		open(1,file='Proportion_MA.inp')
		read(1,*)
		do while(.true.)
			read(1,*,iostat=nn)nza,fa,fy
			if(nn/=0)exit

			nac(i)=nza
			fac1(i)=fa
			fac2(i)=fy*xma
			fmac(i)=fa*fy*xma
			sumac=sumac+fmac(i)
			sumacx=sumacx+fmac(i)
!c	write(*,*)nza,fa,fy,fmac(i),xma,sumac
		i=i+1
		end do
		close(1)
	number=i-1
	endif
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c	write(*,*)numu233,number

!c	����������
	summ0=v0*fp0

	sum=sumac

!c	��������ϵ��
	rat=summ0/sumac
!c	write(*,*)summ0,sum,rat,sumacx


	ts=0.0
      mth=0.0
	do j=1,number
		fwac(j)=fmac(j)*rat
		if(fac1(j).ge.220)then
			ts=ts+fwac(j)
              if((nac(j)/1000).eq.90)mth=mth+fwac(j)
              if((nac(j)/1000).eq.92)mu=mu+fwac(j)              
!c			write(*,*)fac1(j),fwac(j),ts
		endif
	enddo

!c	�غ�����	
	sumac1=ts*1.0e-6
	sumhm0=ts
	open(1,file='M_sum_HM1.dat')
	write(1,*)"�ؽ���������"
	write(1,'(2x,E15.8)')sumhm0
	close(1)
      
      open(1,file='M_sum_Th1.dat')
      write(1,*)"Th������"
      write(1,'(2x,E15.8)')mth
	close(1)
      
      open(1,file='M_sum_U31.dat')
      write(1,*)"U3������"
      write(1,'(2x,E15.8)')mu
	close(1)
!c	stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
100	format(a)
200	format(4a)
300	format(a,E15.8,2x,E15.8)
400	format(4x,i3,a)
500	format(4x,i3,3a)
600	format(a,2x,i3,a)
700	format(a,3(2x,i3,a))
800	format(a,2x,i3,5(2x,E15.8))
900	format(a,2x,a,4(2x,E15.8))

	end subroutine sumhm10
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	subroutine new_CSAS
	character(400)rd

	open(5,file='CSAS.inp')

	open(1,file='scale_CSAS.inp')

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit

		if(rd(1:1).eq.'<')then
			open(2,file=''//rd(2:7)//'.inp')	
			do while(.true.)
				read(2,100,iostat=nn)rd
				if(nn/=0)exit
				write(5,100)rd
			enddo
			close(2)
		else
			write(5,100)rd
		endif
	enddo
	close(1)
	close(5)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

100	format(a)	
    end subroutine new_CSAS
    
    
!    ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine read_keff(ncheck)
	character(120)rd
	ncheck=0
	open(1,file='CSAS.out')
	open(2,file='keff0.dat')

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	
	if(rd(2:46).eq."       ***        best estimate system k-eff")then
		write(2,200)rd(75:101)
		ncheck=1
	endif
	enddo
	close(1)
	close(2)
!c	stop
100	format(a)
200	format(2a)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	endsubroutine read_keff
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    end module INITIAL_CSAS
    
module DEPLETION

    contains
    subroutine Seek_Equilibrum(mm,wpower)!m��ʾ��ѭ������
        character(200)rdd
        character(5)rd
        dimension num_nucl(20),num_1(20),f_nucl(20)
        dimension FN_NUCL(13)
        
        call nucl_id                     
        
        open(1,file='neutron_flux.dat')
        read(1,*)
        read(1,*) flux
        close(1)    
        
     !   open(1,file="powe_1.inp")
	    !write(1,*)" 59** ",flux,flux,flux,flux,flux
     !   write(1,*)flux,flux,flux,flux,flux
	    !close(1)
        
        if(mm.gt.1)then
            open(1,file='con_fe_old.inp')
            read(1,*)feedth,feedu3
            close(1)
        endif
        
            open(1,file='powe_1.inp')
            write(1,*)" 59**   f",flux
            close(1)
        
            open(2,file='56$$a3.inp')
            write(2,1000)" 56$$ a3 ",1," e"
            close(2)
        
        nth=2
        num_1=2
        num_nucl(1:2)=(/902320,922330/)
        mmm=1
        do while(.true.)
        !if((mm.eq.1).and.(mmm.le.2))then
        !    open(1,file='powe_1.inp')
        !    write(1,*)" 59** ",flux,flux,flux,flux,flux
        !    write(1,*)flux,flux,flux,flux,flux
        !    close(1)
        !
        !    open(2,file='56$$a3.inp')
        !    write(2,1000)" 56$$ a3 ",1," e"
        !    close(2)
        !else
        !    open(1,file='powe_1.inp')
        !    write(1,*)" 58** ",wpower,wpower,wpower,wpower,wpower
        !    write(1,*)wpower,wpower,wpower,wpower,wpower
        !    close(1)
        !
        !    open(2,file='56$$a3.inp')
        !    write(2,1000)" 56$$ a3 ",0," e"
        !    close(2)
        !endif
        if((mm.eq.1).and.(mmm.eq.1))then
            f_nucl(1)=0.2e-5
            f_nucl(2)=0.2e-5
        else
            f_nucl(1)=feedth
            f_nucl(2)=feedu3
        endif   
        open(5,file='con_fe.inp')
            write(5,1001)"76$$ ",(num_nucl(i),i=1,nth)
            write(5,1002)"77** ",(f_nucl(i),i=1,nth)
            write(5,1003)"78$$ ",(num_1(i),i=1,nth)
            close(5)
        open(2,file='fe_561.inp')
          write(2,1000)"56$$ a9 ",nth,' e'
	    close(2)
                  
        ii=0
        NNN=1!---------�ж��Ƿ�ﵽ����Ҫ���һ��ָ��֮һ��NNN=9�������ı�Ҫ����
        nnnn=1
        
       call pow_60(ii,1000.0)  
       do while(.true.)
                     
            call new_origens(ii)
            call system('E:\scale6.1\cmds\runscale origens.inp')
            call out_M_stat1(ii,ncheck,NNN)
!            ncheck=1
            if(ncheck.eq.0)then
                call save_N(mm,mmm)
                exit
            endif
            
            call change(ii)
         !   call resumhm
         !   
         !   open(1,file='M_sum_HM.dat')
	        !read(1,*)
	        !read(1,*)sumhm0
	        !close(1)
         !
	        !fx=wpower/sumhm0*1.0e6
         !   open(1,file="powe_1.inp")
	        !write(1,*)" 58** ",fx,fx,fx,fx
	        !close(1)
            
            !if((mod(ii,50).eq.0).and.(ii.ne.50))then
            if(ii.gt.1)then
            open(1,file='criteria.dat')
            read(1,*)
            read(1,*)(FN_NUCL(i),i=1,2)
            read(1,*)
            !read(1,*)
            !read(1,*)(FN_NUCL(i),i=3,6)
            !read(1,*)
            !read(1,*)(FN_NUCL(i),i=7,10)
            !read(1,*)
            !read(1,*)(FN_NUCL(i),i=11,13)
            !close(1)
            njudge=0
            do while(.true.)
                read(1,*,iostat=nn)nza,ferror
                if(nn/=0)exit
                if(ferror.gt.0.0)then
                    njudge=1
                    exit
                endif
            enddo
            
            !do III=3,4
            !    if(FN_NUCL(III).gt.0.0)then 
            !        njudge=1
            !        exit
            !    endif
            !enddo
            
            
            
            if(njudge.eq.0)then
                if(NNN.eq.9)then
                    if(nnnn.eq.1)then
                        call pow_60(ii,100.0)
                        nnnn=nnnn+1
                        NNN=1
                    elseif(nnnn.eq.2)then
                        call pow_60(ii,10.0)
                        nnnn=nnnn+1
                        NNN=1
                    else
                        call save_N(mm,mmm)
                        exit
                    endif                   
                else
                    NNN=NNN+1
                endif
            endif
            
            endif
            
            !if(ii.ge.100)then
            !    call save_N(mmm)
            !    exit
            !endif
            
            
            ii=ii+1
        enddo
        
        
        !if(mmm.gt.1)then
        !call check_N(mmm)!---�ж���ѭ���Ƿ�����       
        !open(1,file='criteria_N.dat')
        !read(1,*)
        !read(1,*)(FN_NUCL(i),i=1,2)
        !read(1,*)
        !read(1,*)
        !read(1,*)(FN_NUCL(i),i=3,6)
        !read(1,*)
        !read(1,*)(FN_NUCL(i),i=7,10)
        !read(1,*)
        !read(1,*)(FN_NUCL(i),i=11,13)
        !close(1)
        !njudge=0
        !do III=1,12
        !    if(FN_NUCL(III).gt.0.003)then 
        !        njudge=1
        !        exit
        !        endif
        !enddo     
        !if(njudge.eq.0)exit
        !endif
!        if(mmm.gt.100)exit
        
        call KEY(feedth,feedu3)
        if(mmm.eq.1)then
            ferror_th_old=0.0
            ferror_u3_old=0.0
        endif
        
        if(mmm.gt.1)then
            open(1,file='criteria_NNN.dat')
            read(1,*)
            read(1,*)ferror_th,ferror_u3
            close(1)
            !if((ferror_th.lt.0.001).and.(ferror_u3.lt.0.001))exit
            !if(mmm.gt.25)then
            !    if((ferror_th.lt.0.003).and.(ferror_u3.lt.0.003))then
            !        exit
            !    else
            !        if(mmm.gt.50)exit
            !    endif
            !endif
            if((ferror_th.eq.ferror_th_old).and.(ferror_u3.eq.ferror_u3_old))exit
            
            ferror_th_old=ferror_th
            ferror_u3_old=ferror_u3
            
            if(mmm.gt.50)exit
            
        endif       
        mmm=mmm+1
        end do
        
        call resumhm
        call new_com
        
        open(1,file='con_fe.inp')
        open(2,file='con_fe_old.inp')
        read(1,*)
        read(1,100)rdd
        write(2,100)rdd(6:60)
        close(1)
        close(2)
        
        
        
!        call sumhm !---����ͳ���ؽ������������ڸ��¹����ܶ�
100     format(a)
1000	format(a,x,i2,a)
1001	format(a,6(x,i8))
1003	format(a,15(x,i3))
1002	format(a,4(x,E15.8))
200	    format(18(2x,E15.8))
    end subroutine Seek_Equilibrum

      
    subroutine new_com
    character(7)nuclx
    character(200)rd
    
    open(1,file='molten-10.inp')
	read(1,*)
	read(1,*)v0  !���������cm^3����

	read(1,*)
	read(1,*)fp0  !�����ܶ�(g/cm^3):

	read(1,*)
	read(1,*)ntemp
    close(1)
    
	
	open(3,file='com_1_x.dat')
	open(4,file='com_1_old.dat')
	do while(.true.)
		read(3,100,iostat=nn)rd
		if(nn/=0)exit
		write(4,*)rd
    end do
    
	close(3)
	close(4)
	
    open(3,file='com_1_x.dat')
	open(5,file='com_10.inp')

	j=0
	read(3,*)
	read(3,*)
	do while(.true.)
		read(3,*,iostat=nn)nuclx,nza,fna,fatoms
		if(nn/=0)exit
        if(nza.eq.892240)exit
        !if((nza.eq.30060).or.(nza.eq.30070).or.(nza.eq.90190).or.(nza.eq.40090))then
        !    open(1,file='com_FLiBe.dat')
        !    do while(.true.)
        !        read(1,*,iostat=mm)nuclx,nza1,fna,fatoms
        !        if(mm/=0)exit
        !        if(nza1.eq.nza)exit
        !    end do
        !    close(1)
        !endif
        
		if(fatoms.ge.1e-20)then
			if(nza.eq.60120)then
				write(5,101)"c      ",' 1 0 ',fatoms,ntemp," end   "
			else
				write(5,101)nuclx,' 1 0 ',fatoms,ntemp," end     "
			endif
        endif
    enddo
    close(3)
    close(5)
    
	open(3,file='com_31.inp')
	open(2,file='com_10.inp')
	do while(.true.)
		read(2,100,iostat=nn)rd
		if(nn/=0)exit
		write(3,100)rd
	enddo
	close(2)
	close(3)
    
101 format(2a,2x,E15.8,i6,a)    
100 format(a)
    end subroutine new_com

 !   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine out_M_stat1(ncycle,ncheck,NNN)
	character(160)rd,rd1
	character(7)nuclx,nucly
	character(2)elment

	ncheck=1

	kn1=mod(ncycle,10000)/1000+48	!������ת�����ַ�,���ļ�����
	kn2=mod(ncycle,1000)/100+48	!������ת�����ַ�,���ļ�����
	kn3=mod(ncycle,100)/10+48	!������ת�����ַ�,���ļ�����
	kn4=mod(ncycle,10)+48		!������ת�����ַ�,���ļ�����


	open(1,file='origens.out')
	open(2,file='out_M_stat.dat')
	ntitle=1
	rd1="    "

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	if(rd(1:10).eq." stop code")then
        ncheck=0
        PAUSE
        close(1)
        close(2)
        return
    endif
    
	if(rd(1:63).eq."      Pass no.   0 First depletion calculation, mixture no.   1")then
		read(1,100)rd1
	endif
		if((rd1(1:11).eq."     power=").or.(rd1(2:12).eq."     power=").or.(rd1(3:13).eq."     power="))then
!c	write(*,100)rd
			if(ntitle.eq.1)then
				write(2,100)rd
				write(2,100)rd1
				ncheck=1
				ntitle=ntitle+1
			endif
			
			if(rd(1:23).eq."                 charge")then
			do while(.true.)
				read(1,100,iostat=nn)rd
					if(rd(1:9).eq."   totals")exit
					if(rd(1:9).eq."1        ")exit
					write(2,100)rd
			enddo
			endif
	if(rd(1:43).eq."        *****              program:  origen")then
	write(2,*)
	endif
	
		endif
	enddo

	close(1)
	close(2)
!c	stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ntype=1
	i=1
	open(2,file='out_M_stat.dat')
	open(1,file='out_M_stat1.dat')
!    open(3,file='out_M_stat_LAST.dat')
	read(2,100)
	read(2,100)
	do while(.true.)
		read(2,100,iostat=nn)rd
		if(nn/=0)exit
		if(rd(1:11).eq."           ")exit

        
  !      if(rd(9:10).ne."  ")then
		!if(rd(8:8).eq." ")then
		!	if(i.eq.2)ntype=ntype+1
		!	i=1
		!	if(rd(9:9).eq." ")then
	 !           write(3,100)rd(1:7),"-",rd(10:11),"  ",rd(143:155),ntype,rd(10:10),rd(11:11)
		!	else
  !              write(3,100)rd(1:7),"-",rd(9:11)," ",rd(143:155),ntype,rd(9:10),rd(11:11)
		!	endif
		!else
		!	if(i.eq.2)ntype=ntype+1
		!	    i=1
	 !           write(3,100)rd(1:7),"-",rd(8:11),rd(143:155)," ",ntype,rd(8:10),rd(11:11)
		!    endif
		!else
		!	i=2
  !      endif
  !      
        
        if(NNN.eq.1)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(143:155),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(143:155),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(143:155)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
        endif
        
        if(NNN.eq.2)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(130:142),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(130:142),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(130:142)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
        endif
        
                if(NNN.eq.3)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(117:129),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(117:129),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(117:129)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
                endif
                
                        if(NNN.eq.4)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(104:116),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(104:116),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(104:116)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
                        endif
                        
                                if(NNN.eq.5)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(91:103),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(91:103),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(91:103)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
                                endif
                                
                                        if(NNN.eq.6)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(78:90),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(78:90),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(78:90)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
                                        endif
                                        
                                                if(NNN.eq.7)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(65:77),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(65:77),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(65:77)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
                                                endif
                                                
                                                        if(NNN.eq.8)then
		if(rd(9:10).ne."  ")then
		if(rd(8:8).eq." ")then
			if(i.eq.2)ntype=ntype+1
			i=1
			if(rd(9:9).eq." ")then
	            write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(52:64),ntype,rd(10:10),rd(11:11)
			else
                write(1,100)rd(1:7),"-",rd(9:11)," ",rd(52:64),ntype,rd(9:10),rd(11:11)
			endif
		else
			if(i.eq.2)ntype=ntype+1
			    i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(52:64)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
                                                        endif
                                                        
    if(NNN.eq.9)then
		if(rd(9:10).ne."  ")then
		    if(rd(8:8).eq." ")then
			    if(i.eq.2)ntype=ntype+1
			    i=1
			    if(rd(9:9).eq." ")then
	                write(1,100)rd(1:7),"-",rd(10:11),"  ",rd(39:51),ntype,rd(10:10),rd(11:11)
			    else
                    write(1,100)rd(1:7),"-",rd(9:11)," ",rd(39:51),ntype,rd(9:10),rd(11:11)
			    endif
		    else
			    if(i.eq.2)ntype=ntype+1
                i=1
	            write(1,100)rd(1:7),"-",rd(8:11),rd(39:51)," ",ntype,rd(8:10),rd(11:11)
		    endif
		else
			i=2
        endif
    endif
    
	enddo
	close(1)

	close(2)
!    close(3)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	open(1,file='M_sum_HM.dat')
	read(1,*)
	read(1,*)sumhm
	close(1)

	open(1,file='input.inp')
	read(1,*)  !
	read(1,*)v1  !���
	read(1,*)  !
	read(1,*)v2  !���
	close(1)

	ntile=1

	open(7,file='com_1_x.dat')
	write(7,100)" nucl nza  na con(atoms) con_real(grams)  con(grams)"
	write(7,100)
	open(4,file='nucl_ID.inp')
	j=1
	read(4,*)
	read(4,*)
	do while(.true.)
		read(4,*,iostat=nn)nuclx,nza,fna,ntype
		if(nn/=0)exit
		ki=0

		open(2,file='out_M_stat1.dat')
		do while(.true.)
			read(2,*,iostat=nn)nucly,fw
			if(nn/=0)exit


			if(nuclx.eq.nucly)then
			fq=fw/1.0e6*sumhm
			fatoms=fq/v1/fna*0.602214199
			write(7,400)nuclx,nza,fna,fatoms,fq,fw

			ki=1
			exit
			endif
		enddo
		close(2)

		if(ki.eq.0)then
			fw=1.0e-14
			fq=fw/1.0e6*sumhm
			fatoms=fq/v1/fna*0.602214199
			write(7,200)nuclx,nza,fna,fatoms,fq,fw
		endif
	enddo
	close(4)
	close(7)

!!----------------------------------------------------------------
!    open(7,file='com_1_xx.dat')
!	write(7,100)" nucl nza  na con(atoms) con_real(grams)  con(grams)"
!	write(7,100)
!	open(4,file='nucl_ID.inp')
!	j=1
!	read(4,*)
!	read(4,*)
!	do while(.true.)
!		read(4,*,iostat=nn)nuclx,nza,fna,ntype
!		if(nn/=0)exit
!		ki=0
!
!		open(2,file='out_M_stat_LAST.dat')
!		do while(.true.)
!			read(2,*,iostat=nn)nucly,fw
!			if(nn/=0)exit
!
!
!			if(nuclx.eq.nucly)then
!			fq=fw/1.0e6*sumhm
!			fatoms=fq/v1/fna*0.602214199
!			write(7,400)nuclx,nza,fna,fatoms,fq,fw
!
!			ki=1
!			exit
!			endif
!		enddo
!		close(2)
!
!		if(ki.eq.0)then
!			fw=1.0e-14
!			fq=fw/1.0e6*sumhm
!			fatoms=fq/v1/fna*0.602214199
!			write(7,200)nuclx,nza,fna,fatoms,fq,fw
!		endif
!	enddo
!	close(4)
!	close(7)
!!------------------------------------------------------------
    
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	open(1,file='com_1_o.dat')
	write(1,100)" nucl nza  na type con(atoms) con_real(grams) con(grams)"
	write(1,100)

	open(2,file='out_M_stat1.dat')
	do while(.true.)
		read(2,*,iostat=nn)nucly,fw,ntype,fna
		if(nn/=0)exit
	
			fq=fw/1.0e6*sumhm
			fatoms=fq/v1/fna*0.602214199
	
		open(3,file='elment.inp')
		do while(.true.)
			read(3,*,iostat=nn)nz,elment
			if(nn/=0)exit

!c			write(*,*)nz,elment
	mx=0
	if(nucly(4:4).eq."m" .or. nucly(5:5).eq."m" .or. nucly(6:6).eq."m" .or. nucly(7:7).eq."m")mx=1

		if(nucly(2:2).eq."-" .and. elment(2:2).eq." ")then
			if(nucly(1:1).eq.elment(1:1))then
				nza=nz*10000+ANINT(fna)*10+mx
!c				write(*,*)nucly,fw,ntype,nza,na,mx
				write(1,300)nucly,nza,na,ntype,fatoms,fq,fw

			endif
			else
			if(nucly(1:2).eq.elment(1:2))then
				nza=nz*10000+ANINT(fna)*10+mx
				write(1,300)nucly,nza,fna,ntype,fatoms,fq,fw
			endif
		endif

		enddo
	close(3)

	enddo
	close(2)
	close(1)


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

100	format(5a,2x,i4,5(2x,a))
200	format(a,(2x,i6),(2x,f9.5),8(2x,E15.8))
300	format(a,3(2x,i6),8(2x,E15.8))
400	format(a,(2x,i6),(2x,f9.5),3(2x,E15.8),2(2x,i6))

	endsubroutine out_M_stat1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    
!    subroutine change(ncycle)
!    character(7)nucly,nuclx
!    character(80)rd
!    
!    ncheck=mod(ncycle,50)
!    
!    if((ncheck.eq.0).and.(ncycle.eq.50))then
!        open(1,file='com_o.dat')
!        open(2,file='com_1_x.dat')
!        do while(.true.)
!            read(2,100,iostat=nn)rd
!            if(nn/=0)exit
!            write(1,100)rd
!        enddo
!        close(1)
!        close(2)
!    else if(ncheck.eq.0)then
!        nr=mod(ncycle/50,2)
!        if(nr.eq.1)then
!            open(1,file='com_o.dat')
!            open(2,file='com_1_x.dat')
!            do while(.true.)
!                read(2,100,iostat=nn)rd
!                if(nn/=0)exit
!                write(1,100)rd
!            enddo
!            close(1)
!            close(2)
!        endif
!        if(nr.eq.0)then
!            open(1,file='com_x.dat')
!            open(2,file='com_1_x.dat')
!            do while(.true.)
!                read(2,100,iostat=nn)rd
!                if(nn/=0)exit
!                write(1,100)rd
!            enddo
!            close(1)
!            close(2)
!        endif
!        call check_converge(nr,ncycle)
!    endif
!    
!    
!    open(1,file='posNC1.inp')
!    open(2,file='com_1_x.dat')
!    read(2,*)
!    read(2,*)
!    j=0
!    do while(.true.)
!            read(2,*,iostat=nn)nucly,nza,fna,ntype,fatoms,fw
!            if(nn/=0)exit
!            if(j.eq.0)then
!                write(1,*)"74**   ",fw
!            else
!                write(1,*)fw
!            endif          
!        j=1
!    enddo
!    close(1)
!    close(2)
!100 format(a)
!    end subroutine change

        
    subroutine change(ncycle)
    character(7)nucly,nuclx
    character(80)rd
    
!    ncheck=ncycle
    
    if(ncycle.eq.1)then
        open(1,file='com_o.dat')
        open(2,file='com_1_x.dat')
        do while(.true.)
            read(2,100,iostat=nn)rd
            if(nn/=0)exit
            write(1,100)rd
        enddo
        close(1)
        close(2)
    else 
        nr=mod(ncycle,2)
        if(nr.eq.1)then
            open(1,file='com_o.dat')
            open(2,file='com_1_x.dat')
            do while(.true.)
                read(2,100,iostat=nn)rd
                if(nn/=0)exit
                write(1,100)rd
            enddo
            close(1)
            close(2)
        endif
        if(nr.eq.0)then
            open(1,file='com_x.dat')
            open(2,file='com_1_x.dat')
            do while(.true.)
                read(2,100,iostat=nn)rd
                if(nn/=0)exit
                write(1,100)rd
            enddo
            close(1)
            close(2)
        endif
        call check_converge(nr,ncycle)
    endif
    
    
    open(1,file='posNC1.inp')
    open(2,file='com_1_x.dat')
    read(2,*)
    read(2,*)
    j=0
    do while(.true.)
            read(2,*,iostat=nn)nucly,nza,fna,fatoms,fq,fw
            if(nn/=0)exit
        if((nza.eq.30070).or.(nza.eq.90190).or.(nza.eq.40090))then
            open(5,file='com_FLiBe.dat')
            do while(.true.)
                read(5,*,iostat=mm)nucly,nza1,fna,fatoms,fq,fw
                if(mm/=0)exit
                if(nza1.eq.nza)exit
            end do
            close(5)
        endif
            if(j.eq.0)then
                write(1,*)"74**   ",fw
            else
                write(1,*)fw
            endif          
        j=1
    enddo
    close(1)
    close(2)
100 format(a)
    end subroutine change
  
    
    subroutine check_converge(ncheck,ncycle)
    character(7) nucly1,nucly2
    dimension FNX(106),FNO(106),N(106)
    
    open(1,file='com_x.dat')
    read(1,*)
    read(1,*)
    fallnucl1=0.0
    sumhm1=0.0    
    fcs71=0.0
    i=1
    do while(.true.)
        read(1,*,iostat=nn)nucly1,nza1,fna1,fatoms1,fq1,fw1
        if(nn/=0)exit
        fallnucl1=fallnucl1+fatoms1
        if(ANINT(fna1).ge.220)then
            sumhm1=sumhm1+fatoms1
        endif
        if(nza1.ge.890000)then
            FNX(i)=fatoms1
            N(i)=nza1
            i=i+1
        endif
        !if(nza1.eq.902320)then
        !    fth21=fatoms1
        !endif
        !if(nza1.eq.922330)then
        !    fu31=fatoms1
        !endif
        !if(nza1.eq.932370)then
        !    fnp1=fatoms1
        !endif
        !if(nza1.eq.942380)then
        !    fpu81=fatoms1
        !endif
        !if(nza1.eq.942390)then
        !    fpu91=fatoms1
        !endif
        !if(nza1.eq.952410)then
        !    fam11=fatoms1
        !endif
        !if(nza1.eq.952430)then
        !    fam31=fatoms1
        !endif
        !if(nza1.eq.962440)then
        !    fcm41=fatoms1
        !endif
        !if(nza1.eq.962450)then
        !    fcm51=fatoms1
        !endif
        !if(nza1.eq.962460)then
        !    fcm61=fatoms1
        !endif
        !if(nza1.eq.551370)then
        !    fcs71=fcs71+fatoms1
        !endif
    enddo
    close(1)
    
    open(2,file='com_o.dat')
    read(2,*)
    read(2,*)
    fallnucl2=0.0
    sumhm2=0.0
    fcs72=0.0
    i=1
    do while(.true.)
        read(2,*,iostat=mm)nucly2,nza2,fna2,fatoms2,fq2,fw2
        if(mm/=0)exit
        fallnucl2=fallnucl2+fatoms2
        if(ANINT(fna2).ge.220)then
            sumhm2=sumhm2+fatoms2
        endif
        if(nza2.ge.890000)then
            FNO(i)=fatoms2
            i=i+1
        endif
        
        !if(nza2.eq.902320)then
        !    fth22=fatoms2
        !endif
        !if(nza2.eq.922330)then
        !    fu32=fatoms2
        !endif
        !if(nza2.eq.932370)then
        !    fnp2=fatoms2
        !endif
        !if(nza2.eq.942380)then
        !    fpu82=fatoms2
        !endif
        !if(nza2.eq.942390)then
        !    fpu92=fatoms2
        !endif
        !if(nza2.eq.952410)then
        !    fam12=fatoms2
        !endif
        !if(nza2.eq.952430)then
        !    fam32=fatoms2
        !endif
        !if(nza2.eq.962440)then
        !    fcm42=fatoms2
        !endif
        !if(nza2.eq.962450)then
        !    fcm52=fatoms2
        !endif
        !if(nza2.eq.962460)then
        !    fcm62=fatoms2
        !endif
        !if(nza2.eq.551370)then
        !    fcs72=fcs72+fatoms2
        !endif
    enddo
    close(2)
    
    !if(ncheck.eq.0)then
        ferror_all=abs(fallnucl1-fallnucl2)/fallnucl2
        ferror_hm=abs(sumhm1-sumhm2)/sumhm2
    !    fth2=abs(fth21-fth22)/fth22
    !    fu3=abs(fu31-fu32)/fu32
    !    fnp=abs(fnp1-fnp2)/fnp2
    !    fpu8=abs(fpu81-fpu82)/fpu82
    !    fpu9=abs(fpu91-fpu92)/fpu92
    !    fam1=abs(fam11-fam12)/fam12
    !    fam3=abs(fam31-fam32)/fam32
    !    fcm4=abs(fcm41-fcm42)/fcm42
    !    fcm5=abs(fcm51-fcm52)/fcm52
    !    fcm6=abs(fcm61-fcm62)/fcm62
    !    fcs7=abs(fcs71-fcs72)/fcs72
    !endif
    !if(ncheck.eq.1)then
    !    ferror_all=abs(fallnucl1-fallnucl2)/fallnucl1
    !    ferror_hm=abs(sumhm1-sumhm2)/sumhm1
    !    fth2=abs(fth21-fth22)/fth21
    !    fu3=abs(fu31-fu32)/fu31
    !    fnp=abs(fnp1-fnp2)/fnp1
    !    fpu8=abs(fpu81-fpu82)/fpu81
    !    fpu9=abs(fpu91-fpu92)/fpu91
    !    fam1=abs(fam11-fam12)/fam11
    !    fam3=abs(fam31-fam32)/fam31
    !    fcm4=abs(fcm41-fcm42)/fcm41
    !    fcm5=abs(fcm51-fcm52)/fcm51
    !    fcm6=abs(fcm61-fcm62)/fcm61
    !    fcs7=abs(fcs71-fcs72)/fcs71
    !endif
    !open(1,file='criteria.dat')
    !write(1,*)"�ܵĺ���Ũ�Ȳв�  �ؽ���Ũ�Ȳв�"
    !write(1,*)ferror_all,ferror_hm
    !write(1,*)"����Ҫ���زв�"
    !write(1,*)"th u3 np pu8"
    !write(1,*)fth2,fu3,fnp,fpu8
    !write(1,*)"pu9 am1 am3 cm4"
    !write(1,*)fpu9,fam1,fam3,fcm4
    !write(1,*)"cm5 cm6 cs137"
    !write(1,*)fcm5,fcm6,fcs7
    !write(1,*)"ȼ���������"
    !write(1,*)ncycle
    !close(1)
    open(1,file='criteria.dat')
    write(1,*)"�ܵĺ���Ũ�Ȳв�  �ؽ���Ũ�Ȳв�"
    write(1,*)ferror_all,ferror_hm
    write(1,*)"���ϵ��������    ������"
    i=1
    do while(.true.)
        if(i.gt.106)exit
        write(1,*)N(i),abs(FNX(i)-FNO(i))/FNO(i)
        i=i+1
    end do
    close(1)
    end subroutine check_converge
    
    subroutine resumhm
    character(7) nucly
    open(1,file='com_1_x.dat')
    read(1,*)
    read(1,*)
    sumhm0=0.0
    do while(.true.)
        read(1,*,iostat=nn)nucly,nza,fna,fatoms,fq,fw
        if(nn/=0)exit
        if(ANINT(fna).ge.220)then
            sumhm0=sumhm0+fq
        endif
    enddo
    close(1)
    open(1,file='M_sum_HM.dat')
    write(1,*)"�ؽ���������"
    write(1,'(2x,E15.8)')sumhm0
    close(1)
    end subroutine resumhm
    

	subroutine pow_60(ncy,time)

	!if(ncy.le.36)then
	!time=10.0
	!elseif(ncy.gt.36 .and. ncy.le.54)then
	!time=20.0
	!else
	!time=30.0
 !   endif 
 !   
 !   time=1000.0

	nt=int(time)
	if(ncy.eq.1)then
		open(1,file="time.dat")
		write(1,*)nt,time
		close(1)
	else
		open(1,file="time.dat")
		read(1,*)ntx
		close(1)

		ntx=ntx+nt
		open(1,file="time.dat")
		write(1,*)ntx,time
		close(1)
	endif

	open(1,file="pow_60.inp")
	write(1,*)" 60**  8i ",1e-10,time
	close(1)

	end subroutine  pow_60

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine nucl_id

	character(4)nuclx
	character(2)nucly
	character(7)nuclz

	character(100)rd

	dimension nza(500),nt(500)

	open(2,file='postID.dat')
	open(1,file='postIDo.inp')
	i=1
	read(1,*)nuclx,n1,n2,n3,n4,n5,n6,n7,n8,n9
	write(2,"((2x,i6))")n1,n2,n3,n4,n5,n6,n7,n8,n9

	do while(.true.)
		read(1,*,iostat=nn)n1,n2,n3,n4,n5,n6,n7,n8,n9
		if(nn/=0)exit
		write(2,"((2x,i6))")n1,n2,n3,n4,n5,n6,n7,n8,n9

	enddo
	write(2,"((2x,i6))")n1

	close(1)
	close(2)


	open(2,file='posLIB.dat')
	open(1,file='posLIBo.inp')
	i=6
	read(1,100)rd
	do j=1,60/2
	write(2,*)rd(i:i+1)
	i=i+2
	enddo

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
!c	write(*,*)rd
	i=2
	do j=1,60/2
	write(2,*)rd(i:i+1)
	i=i+2
	enddo
	enddo

	close(1)
	close(2)
	
	i=1
	open(2,file='posLIB.inp')
	do while(.true.)
		read(2,*,iostat=nn)ntype
		if(nn/=0)exit
	nt(i)=ntype
	i=i+1
	enddo
	close(2)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc


	i=1

	open(1,file='nucl_ID1.dat')
	open(2,file='postID.dat')
	do while(.true.)
		read(2,*,iostat=nn)n1
		if(nn/=0)exit
		nza(i)=n1
		nz=nza(i)/10000
		na=(nza(i)-nz*10000)/10
		m=nza(i)-nz*10000-na*10
		open(3,file="elment.inp")
		do while(.true.)
			read(3,*,iostat=nn)nz2,nucly
			if(nn/=0)exit
		
		if(nz.eq.nz2)then
			if(nucly(2:2).ne." ")then

	if(na.le.9)then
	if(m.eq.1)write(1,"(2a,i1,a,5i9)")nucly,"-",na,"m",nza(i),na,nt(i)
	if(m.ne.1)write(1,"(2a,i1,5i9)")nucly,"-",na,nza(i),na,nt(i)
	endif

	if(na.ge.10 .and. na.le.99)then
	if(m.eq.1)write(1,"(2a,i2,a,5i9)")nucly,"-",na,"m",nza(i),na,nt(i)
	if(m.ne.1)write(1,"(2a,i2,5i9)")nucly,"-",na,nza(i),na,nt(i)
	endif

	if(na.ge.100)then
	if(m.eq.1)write(1,"(2a,i3,a,5i9)")nucly,"-",na,"m",nza(i),na,nt(i)
	if(m.ne.1)write(1,"(2a,i3,5i9)")nucly,"-",na,nza(i),na,nt(i)
	endif

				else


	if(na.le.9)then
	if(m.eq.1)write(1,"(2a,i1,a,5i9)")nucly(1:1),"-",na,"m",nza(i),na,nt(i)
	if(m.ne.1)write(1,"(2a,i1,5i9)")nucly(1:1),"-",na,nza(i),na,nt(i)
	endif

	if(na.ge.10 .and. na.le.99)then
	if(m.eq.1)write(1,"(2a,i2,a,5i9)")nucly(1:1),"-",na,"m",nza(i),na,nt(i)
	if(m.ne.1)write(1,"(2a,i2,5i9)")nucly(1:1),"-",na,nza(i),na,nt(i)
	endif

	if(na.ge.100)then
	if(m.eq.1)write(1,"(2a,i3,a,5i9)")nucly(1:1),"-",na,"m",nza(i),na,nt(i)
	if(m.ne.1)write(1,"(2a,i3,5i9)")nucly(1:1),"-",na,nza(i),na,nt(i)
	endif
				endif
!c			else
!c			write(*,*)i,nza(i)
			exit
			endif
		enddo
		close(3)
	i=i+1
	enddo
	close(2)
	close(1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

	open(1,file='postIDo.inp')
	open(2,file='postID.inp')
	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	write(2,100)rd
	enddo
	close(2)
	close(1)

	open(1,file='posLIBo.inp')
	open(2,file='posLIB.inp')
	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	write(2,100)rd
	enddo
	close(2)
	close(1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc

	!open(1,file='nucl_ID1.dat')
	!open(2,file='nucl_ID.inp')
	!write(2,"(a)")"nuclear    nza       na     ntype"
	!write(2,"(a)")
	!do while(.true.)
	!	read(1,*,iostat=nn)nuclz,nza1,na1,ntype
	!	if(nn/=0)exit
	!write(2,"(a,5(2x,i6))")nuclz,nza1,na1,ntype
	!enddo
	!close(2)
	!close(1)

	open(2,file='56$$13.inp')
	write(2,"(a,i6,a)")" 56$$ a13  ",429," e"
	close(2)

!c	stop

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
100	format(9a)
200	format(2x,e10.4)
300	format(a,2x,i3,3(2x,e10.4))

	end subroutine nucl_id
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine new_origens(ncycle)
	character(400)rd
	character(7)nuclx
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	open(5,file='origens.inp')

	open(1,file='tmsr-origens.inp')

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit

		if(rd(1:1).eq.'<')then
			open(2,file=''//rd(2:7)//'.inp')	
			do while(.true.)
				read(2,100,iostat=nn)rd
				if(nn/=0)exit
				write(5,100)rd
			enddo
			close(2)
		else
			write(5,100)rd
		endif
	enddo
	close(1)
	close(5)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

100	format(a)	
200	format(2(x,i6))	
	end subroutine new_origens
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine new_TRITON
	character(400)rd

	open(5,file='TRITON.inp')

	open(1,file='scale_TRITON.inp')

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit

		if(rd(1:1).eq.'<')then
			open(2,file=''//rd(2:7)//'.inp')	
			do while(.true.)
				read(2,100,iostat=nn)rd
				if(nn/=0)exit
				write(5,100)rd
			enddo
			close(2)
		else
			write(5,100)rd
		endif
	enddo
	close(1)
	close(5)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

100	format(a)	
    end subroutine new_TRITON
    
    subroutine save_N(mm,mmm)
    character(200)rd
    km1=mod(mm,10000)/1000+48	!������ת�����ַ�,���ļ�����
	km2=mod(mm,1000)/100+48	!������ת�����ַ�,���ļ�����
	km3=mod(mm,100)/10+48	!������ת�����ַ�,���ļ�����
	km4=mod(mm,10)+48		!������ת�����ַ�,���ļ�����
    
    kn1=mod(mmm,10000)/1000+48	!������ת�����ַ�,���ļ�����
	kn2=mod(mmm,1000)/100+48	!������ת�����ַ�,���ļ�����
	kn3=mod(mmm,100)/10+48	!������ת�����ַ�,���ļ�����
	kn4=mod(mmm,10)+48		!������ת�����ַ�,���ļ�����
    open(1,file='com_1_x.dat')
    open(2,file='com'//char(km1)//char(km2)//char(km3)//char(km4)//'_'//char(kn1)//char(kn2)//char(kn3)//char(kn4)//'.dat')
    do while(.true.)
        read(1,100,iostat=nn)rd
        if(nn/=0)exit
        write(2,100)rd
    end do
    close(1)
    close(2)
100 format(a)
    end subroutine save_N
    
    subroutine check_N(mmm)
    character(7)nucly1,nucly2
    nnn=mmm-1
    km1=mod(mmm,10000)/1000+48	!������ת�����ַ�,���ļ�����
	km2=mod(mmm,1000)/100+48	!������ת�����ַ�,���ļ�����
	km3=mod(mmm,100)/10+48	!������ת�����ַ�,���ļ�����
	km4=mod(mmm,10)+48		!������ת�����ַ�,���ļ�����
    
    kn1=mod(nnn,10000)/1000+48	!������ת�����ַ�,���ļ�����
	kn2=mod(nnn,1000)/100+48	!������ת�����ַ�,���ļ�����
	kn3=mod(nnn,100)/10+48	!������ת�����ַ�,���ļ�����
	kn4=mod(nnn,10)+48		!������ת�����ַ�,���ļ�����
    
    open(1,file='com_'//char(km1)//char(km2)//char(km3)//char(km4)//'.dat')
    read(1,*)
    read(1,*)
    fallnucl1=0.0
    sumhm1=0.0    
    fcs71=0.0
    do while(.true.)
        read(1,*,iostat=nn)nucly1,nza1,fna1,fatoms1,fq1,fw1
        if(nn/=0)exit
        fallnucl1=fallnucl1+fatoms1
        if(ANINT(fna1).ge.220)then
            sumhm1=sumhm1+fatoms1
        endif
        if(nza1.eq.902320)then
            fth21=fatoms1
        endif
        if(nza1.eq.922330)then
            fu31=fatoms1
        endif
        if(nza1.eq.932370)then
            fnp1=fatoms1
        endif
        if(nza1.eq.942380)then
            fpu81=fatoms1
        endif
        if(nza1.eq.942390)then
            fpu91=fatoms1
        endif
        if(nza1.eq.952410)then
            fam11=fatoms1
        endif
        if(nza1.eq.952430)then
            fam31=fatoms1
        endif
        if(nza1.eq.962440)then
            fcm41=fatoms1
        endif
        if(nza1.eq.962450)then
            fcm51=fatoms1
        endif
        if(nza1.eq.962460)then
            fcm61=fatoms1
        endif
        if(nza1.eq.551370)then
            fcs71=fcs71+fatoms1
        endif
    enddo
    close(1)

    open(2,file='com_'//char(kn1)//char(kn2)//char(kn3)//char(kn4)//'.dat')
    read(2,*)
    read(2,*)
    fallnucl2=0.0
    sumhm2=0.0
    do while(.true.)
        read(2,*,iostat=mm)nucly2,nza2,fna2,fatoms2,fq2,fw2
        if(mm/=0)exit
        fallnucl2=fallnucl2+fatoms2
        if(ANINT(fna2).ge.220)then
            sumhm2=sumhm2+fatoms2
        endif
        if(nza2.eq.902320)then
            fth22=fatoms2
        endif
        if(nza2.eq.922330)then
            fu32=fatoms2
        endif
        if(nza2.eq.932370)then
            fnp2=fatoms2
        endif
        if(nza2.eq.942380)then
            fpu82=fatoms2
        endif
        if(nza2.eq.942390)then
            fpu92=fatoms2
        endif
        if(nza2.eq.952410)then
            fam12=fatoms2
        endif
        if(nza2.eq.952430)then
            fam32=fatoms2
        endif
        if(nza2.eq.962440)then
            fcm42=fatoms2
        endif
        if(nza2.eq.962450)then
            fcm52=fatoms2
        endif
        if(nza2.eq.962460)then
            fcm62=fatoms2
        endif
        if(nza2.eq.551370)then
            fcs72=fatoms2
        endif
    enddo
    close(2)
    
        ferror_all=abs(fallnucl1-fallnucl2)/fallnucl2
        ferror_hm=abs(sumhm1-sumhm2)/sumhm2
        fth2=abs(fth21-fth22)/fth22
        fu3=abs(fu31-fu32)/fu32
        fnp=abs(fnp1-fnp2)/fnp2
        fpu8=abs(fpu81-fpu82)/fpu82
        fpu9=abs(fpu91-fpu92)/fpu92
        fam1=abs(fam11-fam12)/fam12
        fam3=abs(fam31-fam32)/fam32
        fcm4=abs(fcm41-fcm42)/fcm42
        fcm5=abs(fcm51-fcm52)/fcm52
        fcm6=abs(fcm61-fcm62)/fcm62
        fcs7=abs(fcs71-fcs72)/fcs72 
        
    open(1,file='criteria_N.dat')
    write(1,*)"�ܵĺ���Ũ�Ȳв�  �ؽ���Ũ�Ȳв�"
    write(1,*)ferror_all,ferror_hm
    write(1,*)"����Ҫ���زв�"
    write(1,*)"th u3 np pu8"
    write(1,*)fth2,fu3,fnp,fpu8
    write(1,*)"pu9 am1 am3 cm4"
    write(1,*)fpu9,fam1,fam3,fcm4
    write(1,*)"cm5 cm6 cs137"
    write(1,*)fcm5,fcm6,fcs7
    write(1,*)"ƽ��̬��ѭ���������"
    write(1,*)mmm
    close(1)        
        end subroutine check_N
    
 !   include 'Source1.f90'  
 subroutine KEY(feedth,feedu3)
    character rd
    character(7) rd2
    character(9) rd1
    character(200)rdd
    dimension N(429),FN(429)  !!!����NΪ���صı�ţ�FNΪ������������ʵ�����
	dimension NO(429),FNO(429)  !---����NO��FNO�ֱ���������ٽ����ĺ��ر�ź����ʵ���
    
    fkeff=1.00250
    
    
    open(1,file='non-leakge.dat')
    read(1,*)fabs
    read(1,*)fleak
    close(1)
        
    fnonleak=fabs/(fabs+fleak)
        
        
    open(1,file='com_1_x.dat')
    open(2,file='com_Ac_old.dat')
    read(1,*)
    read(1,*)
    i=1
    do while(.true.)
        read(1,*,iostat=nn)rd,nza,fna,fatoms,frmass,fmass
        if(nn/=0)exit
        N(i)=nza
        FN(i)=frmass/(fna*1.0)*0.602214199 !!!!���Ӹ������Ѿ����������
        if(nza.eq.902320)write(2,*)nza,FN(i)
        if(nza.eq.922330)write(2,*)nza,FN(i)
        if(nza.eq.932370)write(2,*)nza,FN(i)
        if(nza.eq.952410)write(2,*)nza,FN(i)
        if(nza.eq.952430)write(2,*)nza,FN(i)
        if(nza.eq.962440)write(2,*)nza,FN(i)
        if(nza.eq.962450)write(2,*)nza,FN(i)
        i=i+1
    end do
    close(1)
    close(2)
    
	
	open(1,file='com_1_old.dat')
    read(1,*)
    read(1,*)
    i=1
    do while(.true.)
        read(1,*,iostat=nn)rd,nza,fna,fatoms,frmass,fmass
        if(nn/=0)exit
        NO(i)=nza
        FNO(i)=frmass/(fna*1.0)*0.602214199 !!!!���Ӹ������Ѿ����������
        i=i+1
    end do
    close(1)

!------------------ ���¿�����Ϊһ���ӳ��� --------------------------
!------------------ ��ȡ����ͨ���ܶ� --------------------------------
    open(1,file='neutron_flux.dat')
    read(1,*)
    read(1,*) flux
    close(1)
    

!---------------
    facp=0.0
    open(1,file='KMT_act.dat')
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn)nza,fac
        if(nn/=0)exit
		if(fac.gt.0.0)then
		if(nza.eq.902320)then
			i=1
			do while(.true.)
                if(i.gt.429)exit
                if(NO(i).eq.nza)then
                    fth=fac/FNO(i)
                exit                
                endif
                i=i+1
            end do
		else if(nza.eq.922330)then
			i=1
			do while(.true.)
                if(i.gt.429)exit
                if(NO(i).eq.nza)then
                    fu3=fac/FNO(i)
                exit                
                endif
                i=i+1
            end do
		else
			i=1
			do while(.true.)
                if(i.gt.429)exit
                if(NO(i).eq.nza)then
                    j=1
					do while(.true.)
						if(j.gt.429)exit
						if(N(j).eq.nza)then
                            if(FNO(i).gt.0.0)then
							    facp=facp+fac/FNO(i)*FN(j)
                            endif
							i=429
							exit                
						endif
						j=j+1
					end do
                exit                
                endif
                i=i+1
            end do
		endif
		endif
    enddo
    close(1)
    
!---------------
    
!!!!!!!!!!!!!!!!!!!!!!!!!! ���Ͽ�����Ϊһ���ӳ��� !!!!!!!!!!!!!!!!!!!!!!!!!!
    
    open(1,file='tot_cs.dat')
    open(2,file='loss_cs_feed.dat')!!!!!!!!���������ʷ����к�����ʧ�������Ϣ
    do while(.true.)
        read(1,*,iostat=nn)nza,rd1,fcs
        if(nn/=0)exit
		if((nza.eq.902320).and.(rd1.eq."tot-cap"))write(2,*)nza,fcs
		if((nza.eq.922330).and.(rd1.eq."tot-cap"))write(2,*)nza,fcs
    end do
    close(1)
    close(2)
    
!!!!!!!!!!!!!!!!!!!!! Ҫ�����ӵ�ʧ�غ㷽�̵�ϵ�����浽�ļ��У�matrix_NCEs.dat(neutron conservation equations) !!!!!!!!!!!!!!!!!!!!!!!!
    fnp=0.0
	fam1=0.0
	fam3=0.0
	fcm4=0.0
	fcm5=0.0
	open(1,file='matrix_NCEs.dat')
    write(1,*)"co_th2  co_u3  co_np7  co_am1  "
    write(1,*)fth,fu3,fnp,fam1
    write(1,*)"co_am3  co_cm4  co_cm5  " !!!!!!co��ʾ��coefficient����������ϵ��
    write(1,*)fam3,fcm4,fcm5
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,0.0
    write(1,*)"co_const  "
    write(1,*)fkeff-facp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!! ���㷴�������ʷ������ϵ�� !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------�����������ǰ��ϵ��-------------   
    feed_Th=-1.0E+16
    feed_U3=-1.0E+16
    feed_MA=0.0
    !feed_MA=-1.0E+16
    !open(1,file='Proportion_MA.inp')
    !read(1,*)
    !do while(.true.)
    !    read(1,*,iostat=nn)nza,fna,fpro
    !    if(nn/=0)exit
    !    if(nza.eq.93237)feed_Np=feed_MA*fpro
    !    if(nza.eq.95241)feed_Am1=feed_MA*fpro
    !    if(nza.eq.95243)feed_Am3=feed_MA*fpro
    !    if(nza.eq.96244)feed_Cm4=feed_MA*fpro
    !    if(nza.eq.96245)feed_Cm5=feed_MA*fpro
    !enddo
    !close(1)
!------------------------------------------------------    
!!!!!!!!!Th232���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Th.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
        if(nn/=0)exit
        i=1
        do while(.true.)
            if(i.gt.429)exit
            if(N(i).eq.nza)then
                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
                exit
            endif
            i=i+1
        end do
    end do
    close(1)
    
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.902310).and.(mt.eq.102))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        
        if((nza.eq.902340).and.(mt.eq.17))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.912320).and.(mt.eq.103))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.912330).and.(mt.eq.28))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.922350).and.(mt.eq.107))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.902330).and.(mt.eq.16))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Th232��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
    open(2,file='DECAY_Th.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.902320)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! ��Th232���������ʷ��̵ĸ���浽�ļ�feed_equation_Th.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_Th.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)fcl,0.0,0.0,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,0.0,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)feed_Th,0,0
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!!!!!!!!!U233���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_U3.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
        if(nn/=0)exit
        i=1
        do while(.true.)
            if(i.gt.429)exit
            if(N(i).eq.nza)then
                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
                exit
            endif
            i=i+1
        end do
    end do
    close(1)
    
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.922320).and.(mt.eq.102))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        
        if((nza.eq.922340).and.(mt.eq.16))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.922350).and.(mt.eq.17))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.922360).and.(mt.eq.37))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.942360).and.(mt.eq.107))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp=fp+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!U233��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
    open(2,file='DECAY_U3.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.922330)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! ��U233���������ʷ��̵ĸ���浽�ļ�feed_equation_U3.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_U3.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)0.0,fcl,0.0,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,0.0,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0,feed_U3,0
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!Np237���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fp=0.0
!    open(1,file='DECAY_Np.dat')
!    read(1,*)
!    read(1,*)
!    do while(.true.)
!        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
!        if(nn/=0)exit
!        i=1
!        do while(.true.)
!            if(i.gt.429)exit
!            if(N(i).eq.nza)then
!                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
!                exit
!            endif
!            i=i+1
!        end do
!    end do
!    close(1)
!    
!    open(1,file='prod_cs_feed.dat')
!    do while(.true.)
!        read(1,*,iostat=nn)nza,fpcs,mt
!        if(nn/=0)exit
!        if((nza.eq.932360).and.(mt.eq.102))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        
!        if((nza.eq.932380).and.(mt.eq.16))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.932390).and.(mt.eq.17))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.942370).and.(mt.eq.103))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952400).and.(mt.eq.107))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!    end do
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Np237��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
!    open(2,file='DECAY_Np.dat')
!    read(2,*)
!    read(2,*)nza,fthalf,fprob
!    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
!    close(2)
!    
!    open(2,file='loss_cs_feed.dat')
!    do while(.true.)
!        read(2,*,iostat=nn)nza,flcs
!        if(nn/=0)exit
!        if(nza.eq.932370)then
!            fcl=fcl+flux*flcs
!            exit
!        endif
!    end do
!    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! ��Np237���������ʷ��̵ĸ���浽�ļ�feed_equation_Np.dat!!!!!!!!!!!!!
!    open(1,file='feed_equation_Np.dat') 
!    write(1,*)"th2 u3 np7 am1"
!    write(1,*)0.0,0.0,fcl,0.0
!    write(1,*)"am3 cm4 cm5"
!    write(1,*)0.0,0.0,0.0
!    write(1,*)"feedTh feedU3 feedMA"
!    write(1,*)0.0,0.0,feed_Np
!    write(1,*)"const"
!    write(1,*)fp
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Am241���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fp=0.0
!    open(1,file='DECAY_Am1.dat')
!    read(1,*)
!    read(1,*)
!    do while(.true.)
!        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
!        if(nn/=0)exit
!        i=1
!        do while(.true.)
!            if(i.gt.429)exit
!            if(N(i).eq.nza)then
!                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
!                exit
!            endif
!            i=i+1
!        end do
!    end do
!    close(1)
!    
!    open(1,file='prod_cs_feed.dat')
!    do while(.true.)
!        read(1,*,iostat=nn)nza,fpcs,mt
!        if(nn/=0)exit
!        if((nza.eq.962410).and.(mt.eq.103))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952400).and.(mt.eq.102))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        
!        if((nza.eq.952420).and.(mt.eq.16))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952430).and.(mt.eq.17))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952421).and.(mt.eq.16))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952440).and.(mt.eq.37))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!    end do
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Am241��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
!    open(2,file='DECAY_Am1.dat')
!    read(2,*)
!    read(2,*)nza,fthalf,fprob
!    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
!    close(2)
!    
!    open(2,file='loss_cs_feed.dat')
!    do while(.true.)
!        read(2,*,iostat=nn)nza,flcs
!        if(nn/=0)exit
!        if(nza.eq.952410)then
!            fcl=fcl+flux*flcs
!            exit
!        endif
!    end do
!    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! ��Am241���������ʷ��̵ĸ���浽�ļ�feed_equation_Am1.dat!!!!!!!!!!!!!
!    open(1,file='feed_equation_Am1.dat') 
!    write(1,*)"th2 u3 np7 am1"
!    write(1,*)0.0,0.0,0.0,fcl
!    write(1,*)"am3 cm4 cm5"
!    write(1,*)0.0,0.0,0.0
!    write(1,*)"feedTh feedU3 feedMA"
!    write(1,*)0.0,0.0,feed_Am1
!    write(1,*)"const"
!    write(1,*)fp
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Am243���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fp=0.0
!    open(1,file='DECAY_Am3.dat')
!    read(1,*)
!    read(1,*)
!    do while(.true.)
!        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
!        if(nn/=0)exit
!        i=1
!        do while(.true.)
!            if(i.gt.429)exit
!            if(N(i).eq.nza)then
!                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
!                exit
!            endif
!            i=i+1
!        end do
!    end do
!    close(1)
!    
!    open(1,file='prod_cs_feed.dat')
!    do while(.true.)
!        read(1,*,iostat=nn)nza,fpcs,mt
!        if(nn/=0)exit
!        if((nza.eq.962430).and.(mt.eq.103))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952420).and.(mt.eq.102))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        
!        if((nza.eq.952421).and.(mt.eq.102))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.972460).and.(mt.eq.107))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.952440).and.(mt.eq.16))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!    end do
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Am243��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
!    open(2,file='DECAY_Am3.dat')
!    read(2,*)
!    read(2,*)nza,fthalf,fprob
!    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
!    close(2)
!    
!    open(2,file='loss_cs_feed.dat')
!    do while(.true.)
!        read(2,*,iostat=nn)nza,flcs
!        if(nn/=0)exit
!        if(nza.eq.952430)then
!            fcl=fcl+flux*flcs
!            exit
!        endif
!    end do
!    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! ��Am243���������ʷ��̵ĸ���浽�ļ�feed_equation_Am3.dat!!!!!!!!!!!!!
!    open(1,file='feed_equation_Am3.dat') 
!    write(1,*)"th2 u3 np7 am1"
!    write(1,*)0.0,0.0,0.0,0.0
!    write(1,*)"am3 cm4 cm5"
!    write(1,*)fc1,0.0,0.0
!    write(1,*)"feedTh feedU3 feedMA"
!    write(1,*)0.0,0.0,feed_Am3
!    write(1,*)"const"
!    write(1,*)fp
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!Cm244���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fp=0.0
!    open(1,file='DECAY_Cm4.dat')
!    read(1,*)
!    read(1,*)
!    do while(.true.)
!        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
!        if(nn/=0)exit
!        i=1
!        do while(.true.)
!            if(i.gt.429)exit
!            if(N(i).eq.nza)then
!                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
!                exit
!            endif
!            i=i+1
!        end do
!    end do
!    close(1)
!    
!    open(1,file='prod_cs_feed.dat')
!    do while(.true.)
!        read(1,*,iostat=nn)nza,fpcs,mt
!        if(nn/=0)exit
!        if((nza.eq.962430).and.(mt.eq.102))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.962450).and.(mt.eq.16))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        
!        if((nza.eq.962460).and.(mt.eq.17))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.962470).and.(mt.eq.37))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!    end do
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Cm244��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
!    open(2,file='DECAY_Cm4.dat')
!    read(2,*)
!    read(2,*)nza,fthalf,fprob
!    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
!    close(2)
!    
!    open(2,file='loss_cs_feed.dat')
!    do while(.true.)
!        read(2,*,iostat=nn)nza,flcs
!        if(nn/=0)exit
!        if(nza.eq.962440)then
!            fcl=fcl+flux*flcs
!            exit
!        endif
!    end do
!    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! ��Cm244���������ʷ��̵ĸ���浽�ļ�feed_equation_Cm4.dat!!!!!!!!!!!!!
!    open(1,file='feed_equation_Cm4.dat') 
!    write(1,*)"th2 u3 np7 am1"
!    write(1,*)0.0,0.0,0.0,0.0
!    write(1,*)"am3 cm4 cm5"
!    write(1,*)0.0,fcl,0.0
!    write(1,*)"feedTh feedU3 feedMA"
!    write(1,*)0.0,0.0,feed_Cm4
!    write(1,*)"const"
!    write(1,*)fp
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!Cm245���������ʷ���ϵ��!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!������fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fp=0.0
!    open(1,file='DECAY_Cm5.dat')
!    read(1,*)
!    read(1,*)
!    do while(.true.)
!        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalfΪ���ذ�˥�ڣ�fprobΪ����˥���֧��
!        if(nn/=0)exit
!        i=1
!        do while(.true.)
!            if(i.gt.429)exit
!            if(N(i).eq.nza)then
!                fp=fp+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
!                exit
!            endif
!            i=i+1
!        end do
!    end do
!    close(1)
!    
!    open(1,file='prod_cs_feed.dat')
!    do while(.true.)
!        read(1,*,iostat=nn)nza,fpcs,mt
!        if(nn/=0)exit
!        if((nza.eq.962440).and.(mt.eq.102))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.962460).and.(mt.eq.16))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        
!        if((nza.eq.962470).and.(mt.eq.17))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.962480).and.(mt.eq.37))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.972450).and.(mt.eq.103))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!        if((nza.eq.982480).and.(mt.eq.107))then
!           i=1
!            do while(.true.)
!                if(i.gt.429)exit
!                if(N(i).eq.nza)then
!                    fp=fp+FN(i)*flux*fpcs
!                    exit
!                endif
!                i=i+1
!            end do 
!        endif
!    end do
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!Cm245��ʧ��ϵ��!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    fcl=0.0 !!!!!!coefficient loss��ʧ��ϵ��!!!!!!!!!
!    open(2,file='DECAY_Cm5.dat')
!    read(2,*)
!    read(2,*)nza,fthalf,fprob
!    fc1=fc1+log(2.0)/fthalf*fprob*(1.0E+24)
!    close(2)
!    
!    open(2,file='loss_cs_feed.dat')
!    do while(.true.)
!        read(2,*,iostat=nn)nza,flcs
!        if(nn/=0)exit
!        if(nza.eq.962450)then
!            fcl=fcl+flux*flcs
!            exit
!        endif
!    end do
!    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! ��Cm245���������ʷ��̵ĸ���浽�ļ�feed_equation_Cm5.dat!!!!!!!!!!!!!
!    open(1,file='feed_equation_Cm5.dat') 
!    write(1,*)"th2 u3 np7 am1"
!    write(1,*)0.0,0.0,0.0,0.0
!    write(1,*)"am3 cm4 cm5"
!    write(1,*)0.0,0.0,fcl
!    write(1,*)"feedTh feedU3 feedMA"
!    write(1,*)0.0,0.0,feed_Cm5
!    write(1,*)"const"
!    write(1,*)fp
!    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call extra_condition
    
    call solve(ff_th,ff_u3)
    
    open(1,file='M_sum_HM.dat')
	read(1,*)
	read(1,*)sumhm
	close(1)
    
    feedth=ff_th*(1.0e-8)*(1.0e6/sumhm)/0.602214199
    feedu3=ff_u3*(1.0e-8)*(1.0e6/sumhm)/0.602214199
    
100 format(a)
    end subroutine KEY
    
    subroutine extra_condition
    !---------------TRU��molŨ�Ȳ���------------------
    !---------------���ؽ�����������/Th��������-------
    character rd
    
    open(1,file='com_1_x.dat')
    read(1,*)
    read(1,*)
    sumhm1=0.0
    sumtru1=0.0
    ftot=0.0
    do while(.true.)
        read(1,*,iostat=nn)rd,nza,fna,fatoms,frmass,fmass
        if(nn/=0)exit
        if(ANINT(fna).ge.220)then
            if(nza.eq.902320)then
                fco_th=(fna*1.0)/0.602214199
            else if(nza.eq.922330)then
                fco_u3=(fna*1.0)/0.602214199
            else
                sumhm1=sumhm1+frmass
            endif
        endif
        
    end do    
    close(1)
    
    open(1,file='M_sum_HM.dat')
    read(1,*)
    read(1,*)sumhm0
    close(1)
    
    open(1,file='Equation_HMconst.dat')
    write(1,*)"th u3 np am1"
    write(1,*)fco_th,fco_u3,0.0,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,0.0,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,0.0
    write(1,*)"const"
    write(1,*)sumhm0-sumhm1
    close(1)
    end subroutine extra_condition
       
!-------------------������ⷽ������ӳ���-----------------------
!-------------------���ı�������3�������ʣ�7������Ũ��
!   ������˳��    ���ӵ�ʧ�غ㷽��       matrix_NCEs.dat
!                   �ؽ��������غ㷽��     Equation_HMconst.dat
!                   TRUĦ��Ũ���غ㷽��    Equation_TRUconst.dat
!                   Th���������ʷ���       feed_equation_Th.dat
!                   U3���������ʷ���       feed_equation_U3.dat
!                   Np���������ʷ���       feed_equation_Np.dat
!                   Am1���������ʷ���      feed_equation_Am1.dat
!                   Am3���������ʷ���      feed_equation_Am3.dat
!                   Cm4���������ʷ���      feed_equation_Cm4.dat
!                   Cm5���������ʷ���      feed_equation_Cm5.dat
!
!
 

 subroutine solve(ff_th,ff_u3)
    integer::i,k,j
    real*8::A(4,4),b(4),x(4)
    real*8::Aup(4,4),bup(4)
    !AbΪ�������  [Ab]
    real*8::Ab(4,5)
    
    open(1,file='matrix_NCEs.dat')
    read(1,*)
    read(1,*)A(1,1),A(1,2)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)A(1,3),A(1,4)
    read(1,*)
    read(1,*)b(1)
    close(1)
    
    open(1,file='Equation_HMconst.dat')
    read(1,*)
    read(1,*)A(2,1),A(2,2)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)A(2,3),A(2,4)
    read(1,*)
    read(1,*)b(2)
    close(1)
    
    open(1,file='feed_equation_Th.dat')
    read(1,*)
    read(1,*)A(3,1),A(3,2)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)A(3,3),A(3,4)
    read(1,*)
    read(1,*)b(3)
    close(1)

    open(1,file='feed_equation_U3.dat')
    read(1,*)
    read(1,*)A(4,1),A(4,2)
    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*)A(4,3),A(4,4)
    read(1,*)
    read(1,*)b(4)
    close(1)
    
    Ab(1:4,1:4)=A

    Ab(:,4+1)=b

    fn_u3=(b(1)*A(2,1)/A(1,1)-b(2))/(A(2,1)/A(1,1)*A(1,2)-A(2,2))
    fn_th=(b(1)-fn_u3*A(1,2))/A(1,1)
    
    open(1,file='com_Ac_old.dat')
    read(1,*)nza,fatomsth
    read(1,*)nza,fatomsu3
    close(1)
    
    ferror_th=abs(fn_th-fatomsth)/fatomsth
    ferror_u3=abs(fn_u3-fatomsu3)/fatomsu3
    
    open(1,file='criteria_NNN.dat')
    write(1,*)"Th-232  U-233"
    write(1,*)ferror_th,ferror_u3
    close(1)
    
    ff_th=(b(3)-A(3,1)*fn_th)/A(3,3)
    ff_u3=(b(4)-A(4,2)*fn_u3)/A(4,4)
end subroutine solve
!    subroutine solve
!    integer::i,k,j
!    real*8::A(10,10),b(10),x(10)
!    real*8::Aup(10,10),bup(10)
!    !AbΪ�������  [Ab]
!    real*8::Ab(10,11)
!    
!    open(1,file='matrix_NCEs.dat')
!    read(1,*)
!    read(1,*)(A(1,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(1,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(1,j),j=8,10)
!    read(1,*)
!    read(1,*)b(1)
!    close(1)
!    
!    open(1,file='feed_equation_Th.dat')
!    read(1,*)
!    read(1,*)(A(2,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(2,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(2,j),j=8,10)
!    read(1,*)
!    read(1,*)b(2)
!    close(1)
!    
!    open(1,file='feed_equation_U3.dat')
!    read(1,*)
!    read(1,*)(A(3,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(3,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(3,j),j=8,10)
!    read(1,*)
!    read(1,*)b(3)
!    close(1)
!    
!    open(1,file='feed_equation_Np.dat')
!    read(1,*)
!    read(1,*)(A(4,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(4,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(4,j),j=8,10)
!    read(1,*)
!    read(1,*)b(4)
!    close(1)
!    
!    open(1,file='feed_equation_Am1.dat')
!    read(1,*)
!    read(1,*)(A(5,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(5,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(5,j),j=8,10)
!    read(1,*)
!    read(1,*)b(5)
!    close(1)
!    
!    open(1,file='feed_equation_Am3.dat')
!    read(1,*)
!    read(1,*)(A(6,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(6,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(6,j),j=8,10)
!    read(1,*)
!    read(1,*)b(6)
!    close(1)
!    
!    open(1,file='feed_equation_Cm4.dat')
!    read(1,*)
!    read(1,*)(A(7,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(7,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(7,j),j=8,10)
!    read(1,*)
!    read(1,*)b(7)
!    close(1)
!    
!    open(1,file='feed_equation_Cm5.dat')
!    read(1,*)
!    read(1,*)(A(8,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(8,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(8,j),j=8,10)
!    read(1,*)
!    read(1,*)b(8)
!    close(1)
!    
!    open(1,file='Equation_TRUconst.dat')
!    read(1,*)
!    read(1,*)(A(9,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(9,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(9,j),j=8,10)
!    read(1,*)
!    read(1,*)b(9)
!    close(1)
!    
!    open(1,file='Equation_HMconst.dat')
!    read(1,*)
!    read(1,*)(A(10,j),j=1,4)
!    read(1,*)
!    read(1,*)(A(10,j),j=5,7)
!    read(1,*)
!    read(1,*)(A(10,j),j=8,10)
!    read(1,*)
!    read(1,*)b(10)
!    close(1)
!    
!    Ab(1:10,1:10)=A
!
!    Ab(:,10+1)=b
!
!
!!-------------------------------
!!  ����� ��˹��ȥ���ĺ��Ĳ���
!    do k=1,10-1
!
!        do i=k+1,10
!   
!            temp=Ab(i,k)/Ab(k,k)
!     
!            Ab(i,:)=Ab(i,:)-temp*Ab(k,:)
!   
!        end do
!
!    end do
!
!!-----------------------------
!! ������һ����Ab�Ѿ���Ϊ������ʽ�ľ���
!!            | *  *  *  *  # |
!!     [A b]= | 0  *  *  *  # |
!!            | 0  0  *  *  # |
!!            | 0  0  0  *  # |
!!
!    Aup(:,:)=Ab(1:10,1:10)
!
!    bup(:)=Ab(:,10+1)
!
!!�����������Ƿ�����Ļش�����
!    call uptri(Aup,bup,x,10)
!
!    open(1,file='Equations_Result.dat')
!    write(1,*)x
!    close(1)
!    end subroutine solve
!    
!    
!    subroutine uptri(A,b,x,NN)
!!---------------------------------subroutine  comment
!!  Version   :  V1.0    
!!  Coded by  :  syz 
!!  Date      :  2010-4-8
!!-----------------------------------------------------
!!  Purpose   :  �����Ƿ�����Ļش�����
!!                 Ax=b
!!-----------------------------------------------------
!!  Input  parameters  :
!!       1.   A(NN,NN)ϵ������
!!       2.   b(NN)������
!!       3.   NN����ά��
!!  Output parameters  :
!!       1.  x  ���̵ĸ�
!!       2.
!!  Common parameters  :
!!
!!----------------------------------------------------
!
!    implicit real*8(a-z)
!
!    integer::i,j,NN
!
!    real*8::A(NN,NN),b(NN),x(NN)
!
!    x(NN)=b(NN)/A(NN,NN)
!
!!�ش�����
!    do i=n-1,1,-1
!   
!        x(i)=b(i)
!        do j=i+1,NN
!            x(i)=x(i)-a(i,j)*x(j)
!        end do
!        x(i)=x(i)/A(i,i)
!
!    end do
!
!    end subroutine uptri
   
    end module DEPLETION
   
    
  
    
    program main
    use INITIAL_CSAS
    use DEPLETION

    call Initial_Keff
    
    	open(1,file='power.dat')
	    read(1,*)
	    read(1,*)wpower
	    close(1)
    
	    open(1,file='M_sum_HM.dat')
	    read(1,*)
	    read(1,*)sumhm0
	    close(1)
    
	    fx=wpower/sumhm0*1.0e6
    
	    open(2,file='burnup.inp')
	    write(2,*)" power=",fx," burn=1.0e-8  end"
	    close(2)
	   
    
        call new_TRITON
        !call system('E:\scale6.1\cmds\runscale TRITON.inp')
        call read_kmt
        call read_keffd(0)
        call cs_file(0)

        mmn=1
        NFINAL=0
        do while(.true.)             

            
            call Seek_Equilibrum(mmn,fx)
            
            !call new_CSAS
            !call system('E:\scale6.1\cmds\runscale CSAS.inp')
            !call read_kmt
            
            
            call new_TRITON
            call system('E:\scale6.1\cmds\runscale TRITON.inp')
            call read_kmt
            call read_keffd(mmn)
            call cs_file(mmn)
                        
            if(mmn.gt.1)then
                call check_keff(mmn)
                open(1,file='keff0_error.dat')
                read(1,*)ferror
                close(1)
                if(ferror.lt.0.001)then
                    NFINAL=NFINAL+1
                else
                    NFINAL=0
                endif
                
 !               if(NFINAL.gt.3)exit
            endif
            
            mmn=mmn+1
        enddo
        
 
    end program main
    
    
    subroutine read_keffd(ncheck)
	character(120)rd
    
    kn1=mod(ncheck,10000)/1000+48	!������ת�����ַ�,���ļ�����
	kn2=mod(ncheck,1000)/100+48	!������ת�����ַ�,���ļ�����
	kn3=mod(ncheck,100)/10+48	!������ת�����ַ�,���ļ�����
	kn4=mod(ncheck,10)+48		!������ת�����ַ�,���ļ�����
    
!	ncheck=0
	open(1,file='TRITON.out')
	open(2,file='keff0.dat')
    
    do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	
	if(rd(2:46).eq."       ***        best estimate system k-eff")exit
    enddo
    
	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	
	if(rd(2:46).eq."       ***        best estimate system k-eff")then
		write(2,200)rd(75:101)
!		ncheck=1
	endif
	enddo
	close(1)
	close(2)
    
    open(2,file='keff0.dat')
    open(3,file='keff0_'//char(kn1)//char(kn2)//char(kn3)//char(kn4)//'.dat')
    read(2,100)rd
    write(3,100)rd
    close(2)
    close(3)
    
!c	stop
100	format(a)
200	format(2a)
    end subroutine read_keffd  

    
    subroutine check_keff(mmm)
    km1=mod(mmm,10000)/1000+48	!������ת�����ַ�,���ļ�����
	km2=mod(mmm,1000)/100+48	!������ת�����ַ�,���ļ�����
	km3=mod(mmm,100)/10+48	!������ת�����ַ�,���ļ�����
	km4=mod(mmm,10)+48		!������ת�����ַ�,���ļ�����
    
    nnn=mmm-1
    
    kn1=mod(nnn,10000)/1000+48	!������ת�����ַ�,���ļ�����
	kn2=mod(nnn,1000)/100+48	!������ת�����ַ�,���ļ�����
	kn3=mod(nnn,100)/10+48	!������ת�����ַ�,���ļ�����
	kn4=mod(nnn,10)+48		!������ת�����ַ�,���ļ�����
 
    open(2,file='keff0_'//char(km1)//char(km2)//char(km3)//char(km4)//'.dat')
    read(2,*)fkeff1
    close(2)
    
    open(3,file='keff0_'//char(kn1)//char(kn2)//char(kn3)//char(kn4)//'.dat')    
    read(3,*)fkeff2
    close(3)
    
    ferror=abs(fkeff2-fkeff1)/fkeff2
    
    open(1,file='keff0_error.dat')
    write(1,*)ferror
    close(1)
    
    end subroutine check_keff
    
    
    
    
    subroutine cs_file(ncycle)    
    character(200)rdd
    
	kn1=mod(ncycle,10000)/1000+48	!������ת�����ַ�,���ļ�����
	kn2=mod(ncycle,1000)/100+48	!������ת�����ַ�,���ļ�����
	kn3=mod(ncycle,100)/10+48	!������ת�����ַ�,���ļ�����
	kn4=mod(ncycle,10)+48		!������ת�����ַ�,���ļ�����    
    
        open(2,file='TRITON.out') 
        open(5,file='neutron_flux.dat')!!!!����ͨ���ܶ��ļ�
        write(5,*)"neutron flux:"
        do while(.true.)
            read(2,100,iostat=nn)rdd
            if(nn/=0)exit
            if(rdd(1:27).eq."          Number (MW/MTIHM)")then
                read(2,100)rdd
                write(5,100)rdd(67:76)
            exit
            endif
            end do
        close(5) 
        close(2)
 
        open(2,file='TRITON.out') 
        open(3,file='non-leakge.dat')!-----��Ӧ�Ѳ�й¶����
        do while(.true.)
            read(2,100,iostat=nn)rdd
            if(nn/=0)exit
            if(rdd(1:13).eq." system total")then
                write(3,100)rdd(61:73),rdd(90:102)
                exit
            endif
        end do
        close(2)
        close(3)
    
    open(2,file='TRITON.out')
    open(3,file='tot_cs.dat')!!!!!�������ӵ�ʧ�غ㷽�̵Ľ�����Ϣ��tot-cap��nu-sigf
    !open(6,file='tot_cs_b.dat')
    !open(7,file='tot_cs_ni.dat')
    
    open(5,file='NN_cs.dat')!-----�������ӵ�ʧ�غ㷽�̵�(n,2n),(n,3n)��Ӧ��
    
    open(4,file='prod_cs_feed.dat')!!!!���������ʷ��̱���Ľ�����Ϣ  

    do while(.true.)
        read(2,100,iostat=nn)rdd
        if(nn/=0)exit
        if(rdd(1:13).eq." the reaction")exit
    end do
    
        
        
    do while(.true.)
        read(2,100,iostat=nn)rdd
        if(nn/=0)exit
        if(rdd(1:6).eq."1cross")exit
    end do
    read(2,*)
    do while(.true.)
        read(2,100,iostat=nn)rdd
        if(nn/=0)exit
        if(rdd(1:13).eq." the reaction")exit
        if((rdd(29:35).eq."tot-cap").or.(rdd(29:35).eq."nu-sigf"))then
            write(3,100)rdd
        endif
        if(rdd(66:74).ne."byproduct")then
            if((rdd(57:63).eq."mt=  17").or.(rdd(57:63).eq."mt=  16"))then
                write(5,*)rdd(22:28),rdd(40:50),rdd(60:63)
            endif
            
!-------------------��ȡTh232��mt=102������Ϣ������ʯī��Ӧ�ʵļ���-------------
            if((rdd(22:27).eq."902320").and.(rdd(57:63).eq."mt= 102"))then
                write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
            endif
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!��ȡ���з��������ʱ���Ľ�����Ϣ��������������!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Th232�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."902310").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."902340").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif        
        if((rdd(22:27).eq."912320").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif        
        if((rdd(22:27).eq."912330").and.(rdd(57:63).eq."mt=  28"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922350").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."902330").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! U233�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."922320").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922340").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922350").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922360").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942360").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Np237�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."932360").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."932380").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."932390").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942370").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952400").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Am241�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."952400").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952420").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952421").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952430").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952440").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962410").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Am243�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."952420").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952421").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952440").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."972460").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962430").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Cm244�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."962430").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962450").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962460").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962470").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Cm245�������������������Ϣ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."962440").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962460").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962470").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962480").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."972450").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."982480").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
    end do
    close(3)
    close(4)
    close(5)
!-----------------------------------------------------------------------------------
!------B4C�������ӵ�ʧ�غ㷽�̵ĺ��ؽ�����Ϣ
!    do while(.true.)
!        read(2,100,iostat=nn)rdd
!        if(nn/=0)exit
!        if(rdd(1:6).eq."1cross")exit
!    end do
!    read(2,*)
!    do while(.true.)
!        read(2,100,iostat=nn)rdd
!        if(nn/=0)exit
!        if(rdd(1:13).eq." the reaction")exit
!        if(rdd(29:35).eq."tot-cap")then
!            write(6,*)rdd(22:28),rdd(40:50)
!        endif
!    enddo
!!-----------------------------------------------------------------------------------
!!------Ni�Ͻ��������ӵ�ʧ�غ㷽�̵ĺ��ؽ�����Ϣ    
!    do while(.true.)
!        read(2,100,iostat=nn)rdd
!        if(nn/=0)exit
!        if(rdd(1:6).eq."1cross")exit
!    end do
!    read(2,*)
!    do while(.true.)
!        read(2,100,iostat=nn)rdd
!        if(nn/=0)exit
!        if(rdd(1:13).eq." the reaction")exit
!        if(rdd(29:35).eq."tot-cap")then
!            write(7,*)rdd(22:28),rdd(40:50)
!        endif
!    enddo    
!    
!    close(6)
!    close(7)
    close(2)


    
    
100 format(a)    
    end subroutine cs_file