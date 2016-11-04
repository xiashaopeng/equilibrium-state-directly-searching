!------------------------------------------
!             平衡态搜索程序
!------------------------------------------
!   Version      : V1.0
!   Coded by     : Xiashaopeng
!   Date         : 2016-8-15
!------------------------------------------
!   Description  : 直接搜索熔盐堆的平衡态
!                  满足各种添料、后处理要求
!------------------------------------------
!   Code Flow    :
!          1.   初始临界计算，提供燃耗计算
!               所需反应率信息
!          2.   燃耗计算，直至核素浓度不变
!          3.   求解添料率代数方程组，更新
!               添料率以及重要核素浓度
!          4.   判断两次迭代步中，重要核素
!               的浓度是否满足收敛要求。如
!               果满足，则进入步骤 5；如果
!               不满足要求，则将计算得到的
!               新的添料率带入步骤 2中，重
!               新进行燃耗计算。
!          5.   用新的核素浓度进行临界计算，
!               如果计算结果满足临界要求，
!               则成功搜索到了平衡态；如果
!               不满足临界要求，则将该步骤
!               的截面、通量等信息带入步骤
!               2中重新进行燃耗计算。
!-----------------------------------------


module INITIAL_CSAS
!------------------------------------------module coment
!   Version      : V1.0
!   Coded by     : Xiashaopeng
!   Date         : 2016-8-15
!------------------------------------------
!   Description  : 程序流程中的第一步
!            为燃耗计算提供所需数据
!------------------------------------------
!   Contains     :
!       1.  Initial_Keff 初始临界计算 而且还是该module的接口
!       2.  sumhm10      核素质量统计
!       3.  molten10     提供scale计算所需核素浓度
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
		write(1,*)"重金属总质量"
		write(1,'(2x,e12.6)')sumhm0
		close(1)

		call molten_10
!		call molten_20

!	call new_CSAS 
!!	stop
!
!	call system("E:\scale6.1\cmds\runscale CSAS.inp")
!	call read_keff(ncheck)
!	if(ncheck.eq.0)then
!		call system("E:\scale6.1\cmds\runscale CSAS.inp")
!		call read_keff(ncheck)
!		if(ncheck.eq.0)stop
!    endif
!
!    call read_kmt
!    
!	open(1,file='input.inp')
!	read(1,*)  !
!	read(1,*)v1  !体积
!	read(1,*)  !
!	read(1,*)v2  !体积
!	read(1,*)  
!	read(1,*)fkeffdown,fkeffup    !keff_down and keff_up
!	close(1)
!
!	open(2,file='keff0.dat')
!	read(2,*)fkeff
!	close(2)
!
!	if(fkeff.gt.fkeffup .or. fkeff.lt.fkeffdown)then
!		write(*,*)"change the molten ",fkeff
!!		stop
!	endif

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
		write(2,*)"核素种类     反应率"
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
    end subroutine
        
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
	read(1,*)v0  !熔盐体积（cm^3）：

	read(1,*)
	read(1,*)fp0  !熔盐密度(g/cm^3):

	read(1,*)
	read(1,*)ntemp  

	read(1,*)
	read(1,*)wpower

!c	!LiF-BeF2-ThF4-UF4-PuF4-MAF4的摩尔浓度配比情
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
	write(1,*)" 功率："
	write(1,*)wpower
	close(1)


!c	U初始比例
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
!c	Pu初始比例
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


!c	MA初始比例
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

!c	熔盐总质量
	summ0=v0*fp0

	sum=sumac

!c	质量比例系数
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

!c	重金属归一化后各核素的质量


	st=0.0
	do j=1,number
	wmac(j)=fwac(j)/sumhm0*1.0e6
	st=st+wmac(j)
!c	write(*,*)nac(j),	wmac(j)
	enddo	

!c	各核素的原子密度
	
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
    jj=0!--用于截断com_1.inp中核素种类
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

	write(8,"(a,2(2x,i6),8(2x,e12.6))")nuclx,nza,na,fatoms,fw,fq,fs
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  统计燃耗前的TRU总mol%   cccccccccccccccccccccccccccccc
      open(8,file='com_1_old.dat')
          read(8,*)	
          read(8,*)
          sumtru0=0.0
      	do while(.true.)
			read(8,*,iostat=nn)nucly,nza,na,fatoms,fw,fq,fs
              if(nn.ne.0)exit
              if(nza.ge.930000)sumtru0=sumtru0+fatoms !!!!!!!统计TRU核素的mol百分比
          enddo
	close(8)     
      
	open(1,file='M_sum_TRU1.dat')
	write(1,*)"TRU总mol%量"
	write(1,'(2x,e12.6)')sumtru0*v0
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

101	format(2a,2x,e12.6,i6,a)
100	format(a)
200	format(4a,i6)
300	format(a,e12.6,2x,e12.6)
400	format(4x,i3,a)
500	format(4x,i3,3a)
600	format(a,2x,i3,a)
700	format(a,3(2x,i3,a))
800	format(a,2x,i3,5(2x,e12.6))
900	format(a,2x,a,4(2x,e12.6))
1000	format(2x,i8,4(2x,e12.6))
1001	format(i7,2x,a)
1002	format(a2,a,i1,2x,e12.6)
1003	format(a2,a,i2,2x,e12.6)
1004	format(a2,a,i3,2x,e12.6)
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
	read(1,*)v0  !熔盐体积（cm^3）：

	read(1,*)
	read(1,*)fp0  !熔盐密度(g/cm^3):

	read(1,*)
	read(1,*)ntemp  

	read(1,*)
	read(1,*)wpower

!c	!LiF-BeF2-ThF4-UF4-PuF4-MAF4的摩尔浓度配比情
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
	write(1,*)" 功率："
	write(1,*)wpower
	close(1)


!c	U初始比例
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
!c	Pu初始比例
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


!c	MA初始比例
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

!c	熔盐总质量
	summ0=v0*fp0

	sum=sumac

!c	质量比例系数
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

!c	重核质量	
	sumac1=ts*1.0e-6
	sumhm0=ts
	open(1,file='M_sum_HM1.dat')
	write(1,*)"重金属总质量"
	write(1,'(2x,e12.6)')sumhm0
	close(1)
      
      open(1,file='M_sum_Th1.dat')
      write(1,*)"Th总质量"
      write(1,'(2x,e12.6)')mth
	close(1)
      
      open(1,file='M_sum_U31.dat')
      write(1,*)"U3总质量"
      write(1,'(2x,e12.6)')mu
	close(1)
!c	stop
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
100	format(a)
200	format(4a)
300	format(a,e12.6,2x,e12.6)
400	format(4x,i3,a)
500	format(4x,i3,3a)
600	format(a,2x,i3,a)
700	format(a,3(2x,i3,a))
800	format(a,2x,i3,5(2x,e12.6))
900	format(a,2x,a,4(2x,e12.6))

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
    subroutine Seek_Equilibrum(mm)!m表示大循环次数
        character(200)rdd
        dimension num_nucl(20),num_1(20),f_nucl(20)
        dimension FN_NUCL(13),FERROR(7)
        dimension FEED(20)
        
        call nucl_id
        
        open(1,file='neutron_flux.dat')
        read(1,*)
        read(1,*)flux
        close(1)
        
        open(1,file='powe_1.inp')
        write(1,*)" 59**   f",flux
        close(1)
        
        open(2,file='56$$a3.inp')
        write(2,1000)" 56$$ a3 ",1," e"
        close(2)
        
        if(mm.gt.2)then
            open(1,file='con_fe_old.inp')
            read(1,*)(FEED(j),j=1,4)
            read(1,*)(FEED(j),j=5,7)
            close(1)
        endif  
        !if(mm.gt.1)then
        !     open(1,file='con_fe_old.inp')
        !    read(1,*)(FEED(j),j=1,4)
        !    read(1,*)(FEED(j),j=5,7)
        !    close(1)
        !endif
        
        nth=7
        num_1=2
        num_nucl(1:7)=(/902320,922330,932370,952410,952430,962440,962450/)
        
        mmm=1
        do while(.true.)
            !if((mm.eq.1).and.(mmm.eq.1))then
            if((mm.le.2).and.(mmm.eq.1))then
                    feed_Th=0.1E-05!----0.2E-05
                    feed_U3=0.1E-05
                    feed_MA=0.1E-06
                    
                    open(1,file='feed_ma.dat')
                    write(1,*)feed_MA
                    close(1)
                    
                    open(1,file='Proportion_MA.inp')
                    read(1,*)
                    do while(.true.)
                        read(1,*,iostat=nn)nza,fna,fpro
                        if(nn/=0)exit
                        if(nza.eq.93237)feed_Np=feed_MA*fpro
                        if(nza.eq.95241)feed_Am1=feed_MA*fpro
                        if(nza.eq.95243)feed_Am3=feed_MA*fpro
                        if(nza.eq.96244)feed_Cm4=feed_MA*fpro
                        if(nza.eq.96245)feed_Cm5=feed_MA*fpro
                    enddo
                    close(1)
                    f_nucl(1) = feed_Th
                    f_nucl(2) = feed_U3
                    f_nucl(3) = feed_Np
                    f_nucl(4) = feed_Am1
                    f_nucl(5) = feed_Am3
                    f_nucl(6) = feed_Cm4
                    f_nucl(7) = feed_Cm5
            else
                f_nucl=FEED
            endif
        
            open(5,file='con_fe.inp')
            write(5,1001)"76$$ ",(num_nucl(i),i=1,nth)
            write(5,1002)"77** ",(f_nucl(i),i=1,nth)
            write(5,1003)"78$$ ",(num_1(i),i=1,nth)
            close(5)
        !endif    
        open(2,file='fe_561.inp')
          write(2,1000)"56$$ a9 ",nth,' e'
	    close(2)
                  
        ii=1
        NNN=1!---------判断是否达到收敛要求的一个指标之一，NNN=9是收敛的必要条件
        nnnn=1
        
        call pow_60(ii,1000.0)
        do while(.true.)
            write(*,*)"临界搜索周期：",mm
            write(*,*)"添料率搜索周期：",mmm
            write(*,*)"燃耗循环周期：",ii
           ! call pow_60(ii,time)           
            call new_origens(ii)
            call system('E:\scale6.1\cmds\runscale origens.inp')
            call out_M_stat1(ii,ncheck,NNN)
            if(ncheck.eq.0)then
                call save_N(mm,mmm)
                exit
            endif
            call change(ii)

            if(ii.gt.1)then
            open(1,file='criteria.dat')
            read(1,*)
            read(1,*)(FN_NUCL(i),i=1,2)
            read(1,*)
            njudge=0
            do while(.true.)
                read(1,*,iostat=nn)nza,ferror1
                if(nn/=0)exit
                if((nza.eq.902320).or.(nza.eq.922330).or.(nza.eq.922320).or.(nza.eq.922340).or.(nza.eq.922350).or.(nza.eq.922360).or.(nza.eq.932370).or.(nza.eq.942380).or.(nza.eq.952410).or.(nza.eq.952430).or.(nza.eq.962440).or.(nza.eq.962450).or.(nza.eq.962460).or.(nza.eq.942390).or.(nza.eq.942400))then
                if(ferror1.gt.0.0)then
                    njudge=1
                    exit
                endif
                endif
            enddo
                  
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
            
            ii=ii+1
        enddo
                
        call KEY(FEED)
            
        
        if(mmm.gt.1)then
!        call check_N(mmm)!---判断子循环是否收敛       
            open(1,file='criteria_NNN.dat')
            
            do j=1,7
                read(1,*,iostat=nn)FERROR(j)
            enddo
            close(1)
            njudge=0
            do III=1,7
                if(FERROR(III).gt.0.003)then 
                    njudge=1
                    exit
                endif
            enddo     
            if(njudge.eq.0)exit
            if(mmm.gt.25)then
                njudge=0
                do III=1,7
                    if(FERROR(III).gt.0.005)then 
                        njudge=1
                        exit
                    endif
                enddo
                if(njudge.eq.0)exit
            endif
            if(mmm.gt.50)exit
        endif
        mmm=mmm+1
        end do
        
        call resumhm
        call new_com
        
        open(1,file='con_fe.inp')
        open(2,file='con_fe_old.inp')
        read(1,*)
        read(1,*)
        read(1,100)rdd
        write(2,100)rdd(6:60)
        read(1,100)rdd
        write(2,100)rdd
        close(1)
        close(2)
        
!        call sumhm !---重新统计重金属质量，用于更新功率密度
1000	format(a,x,i2,a)
1001	format(a,6(x,i8))
1003	format(a,15(x,i3))
1002	format(a,4(x,e12.6))
200     format(18(2x,e12.6))
100     format(a)
    end subroutine Seek_Equilibrum

    
    subroutine new_com
    character(7)nuclx
    character(100)rd
	open(1,file='molten-10.inp')
	read(1,*)
	read(1,*)v0  !熔盐体积（cm^3）：

	read(1,*)
	read(1,*)fp0  !熔盐密度(g/cm^3):

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
    
    read(3,*)
    read(3,*)
    do while(.true.)
        read(3,*,iostat=nn)nuclx,nza,fna,fatoms
        if(nn/=0)exit
        if(nza.eq.892240)exit
        if(fatoms.ge.1e-20)then
            if(nza.eq.60120)then
                write(5,101)"c      ",' 1 0 ',fatoms,ntemp," end   "
            else
                write(5,101)nuclx,' 1 0 ',fatoms,ntemp," end   "
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
    
100 format(a)
101 format(2a,2x,e12.6,i6,a)
    end subroutine new_com
    
 !   ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine out_M_stat1(ncycle,ncheck,NNN)
	character(160)rd,rd1
	character(7)nuclx,nucly
	character(2)elment

	ncheck=1

	kn1=mod(ncycle,10000)/1000+48	!将数字转换成字符,做文件名用
	kn2=mod(ncycle,1000)/100+48	!将数字转换成字符,做文件名用
	kn3=mod(ncycle,100)/10+48	!将数字转换成字符,做文件名用
	kn4=mod(ncycle,10)+48		!将数字转换成字符,做文件名用


	open(1,file='origens.out')
	open(2,file='out_M_stat.dat')
	ntitle=1
	rd1="    "

	do while(.true.)
		read(1,100,iostat=nn)rd
		if(nn/=0)exit
	if(rd(1:10).eq." stop code")then
        ncheck=0
 !       PAUSE
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
	read(1,*)v1  !体积
	read(1,*)  !
	read(1,*)v2  !体积
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
    
    subroutine change(ncycle)
    character(7)nucly,nuclx
    character(80)rd
    
    !ncheck=mod(ncycle,50)
    
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
            read(2,*,iostat=nn)nucly,nza,na,fatoms,fg,fw
            if(nn/=0)exit
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
        endif
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
        endif
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
    open(1,file='criteria.dat')
    write(1,*)"总的核素浓度残差  重金属浓度残差"
    write(1,*)ferror_all,ferror_hm
    write(1,*)"各重要核素残差"
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
    write(1,*)"重金属总质量"
    write(1,'(2x,e12.6)')sumhm0
    close(1)
    end subroutine resumhm
    

	subroutine pow_60(ncy,time)

	!if(ncy.le.36)then
	!time=10.0
	!elseif(ncy.gt.36 .and. ncy.le.54)then
	!time=20.0
	!else
	!time=30.0
	!endif 

    !time=3650.0
    
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
	write(1,*)" 60** ",time/10.0,2.0*time/10.0,3.0*time/10.0,4.0*time/10.0,5.0*time/10.0
	write(1,*)6.0*time/10.0,7.0*time/10.0,8.0*time/10.0,9.0*time/10.0,time
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
    kn1=mod(mmm,10000)/1000+48	!将数字转换成字符,做文件名用
	kn2=mod(mmm,1000)/100+48	!将数字转换成字符,做文件名用
	kn3=mod(mmm,100)/10+48	!将数字转换成字符,做文件名用
	kn4=mod(mmm,10)+48		!将数字转换成字符,做文件名用
    open(1,file='com_1_x.dat')
    open(2,file='com_'//char(kn1)//char(kn2)//char(kn3)//char(kn4)//'.dat')
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
    km1=mod(mmm,10000)/1000+48	!将数字转换成字符,做文件名用
	km2=mod(mmm,1000)/100+48	!将数字转换成字符,做文件名用
	km3=mod(mmm,100)/10+48	!将数字转换成字符,做文件名用
	km4=mod(mmm,10)+48		!将数字转换成字符,做文件名用
    
    kn1=mod(nnn,10000)/1000+48	!将数字转换成字符,做文件名用
	kn2=mod(nnn,1000)/100+48	!将数字转换成字符,做文件名用
	kn3=mod(nnn,100)/10+48	!将数字转换成字符,做文件名用
	kn4=mod(nnn,10)+48		!将数字转换成字符,做文件名用
    
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
    fcs72=0.0
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
            fcs72=fcs72+fatoms2
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
    write(1,*)"总的核素浓度残差  重金属浓度残差"
    write(1,*)ferror_all,ferror_hm
    write(1,*)"各重要核素残差"
    write(1,*)"th u3 np pu8"
    write(1,*)fth2,fu3,fnp,fpu8
    write(1,*)"pu9 am1 am3 cm4"
    write(1,*)fpu9,fam1,fam3,fcm4
    write(1,*)"cm5 cm6 cs137"
    write(1,*)fcm5,fcm6,fcs7
    close(1)        
        end subroutine check_N
    
!    include 'Source1.f90'  
subroutine KEY(FEED)
    character rd
    character(7) rd2
    character(9) rd1
    character(200)rdd
    dimension FEED(20) !----存储各个核素的添料率信息
    dimension N(429),FN(429)  !!!数组N为核素的编号，FN为各核素相对物质的量。
    
    fnonleak=1.0
    fkeff=1.0725
    
    open(1,file='com_1_x.dat')
    open(2,file='com_Ac_old.dat')
    read(1,*)
    read(1,*)
    i=1
    do while(.true.)
        read(1,*,iostat=nn)rd,nza,fna,fatoms,frmass,fmass
        if(nn/=0)exit
        N(i)=nza
        FN(i)=frmass/(fna*1.0)*0.602214199 !!!!粒子个数，已经包含体积了
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

    !------------提取中子通量密度---------------
    open(1,file='neutron_flux.dat')
    read(1,*)
    read(1,*) flux
    close(1)
    !-----------------------------------------------------
    !---------------
    fac_el=0.0
    open(1,file='KMT_act.dat')
    read(1,*)
    read(1,*)facth
    do while(.true.)
        read(1,*,iostat=nn)fac
        if(nn/=0)exit
        fac_el=fac_el+fac
    enddo
    close(1)
    
    open(1,file='com_1_old.dat')
    read(1,*)
    read(1,*)
    i=1
    do while(.true.)
        read(1,*,iostat=nn)rd,nza,fna,fatoms,frmass,fmass
        if(nn/=0)exit
        if(nza.eq.902320)then
            fnth=frmass/(fna*1.0)*0.602214199
            exit
        endif
        
    end do
    close(1)
!---------------
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fcs,mt
        if(nn/=0)exit
        if((nza.eq.902320).and.(mt.eq.102))then
            i=1
            fag=(fac_el/facth)*fnth*fcs
            exit                
        endif
    end do
    close(1)
    
    open(1,file='tot_cs.dat')
    open(2,file='loss_cs_feed.dat')!!!!!!!!反推添料率方程中核素消失项截面信息
    fa = fag !!!!!!!中子得失守恒方程的中子消失项
    ff = 0.0 !!!!!!!中子得失守恒方程的中子产生项
    fu3=0.0 !!!!!!!中子得失守恒方程的U3前系数
    fth=0.0 !!!!!!!中子得失守恒方程的Th前系数
    fnp=0.0 !!!!!!!中子得失守恒方程的Np237前系数
    fam1=0.0 !!!!!!!中子得失守恒方程Am241前系数
    fam3=0.0 !!!!!!!中子得失守恒方程Am243前系数
    fcm4=0.0 !!!!!!!中子得失守恒方程Cm244前系数
    fcm5=0.0 !!!!!!!中子得失守恒方程Cm245前系数
    do while(.true.)
        read(1,*,iostat=nn)nza,rd1,fcs
        if(nn/=0)exit
        
        i=1
        do while(.true.)
            if(i.gt.429)exit
            if(N(i).eq.nza)then
                if(nza.eq.922330) then
                    if(rd1.eq."tot-cap") then
                        fu3=fu3-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fu3=fu3+fnonleak/fkeff*fcs
                elseif(nza.eq.902320)then
                    if(rd1.eq."tot-cap") then
                        fth=fth-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fth=fth+fnonleak/fkeff*fcs
                elseif(nza.eq.932370)then
                    if(rd1.eq."tot-cap") then
                        fnp=fnp-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fnp=fnp+fnonleak/fkeff*fcs
                elseif(nza.eq.952410)then
                    if(rd1.eq."tot-cap") then
                        fam1=fam1-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fam1=fam1+fnonleak/fkeff*fcs 
                elseif(nza.eq.952430)then
                    if(rd1.eq."tot-cap") then
                        fam3=fam3-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fam3=fam3+fnonleak/fkeff*fcs 
                elseif(nza.eq.962440)then
                    if(rd1.eq."tot-cap") then
                        fcm4=fcm4-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fcm4=fcm4+fnonleak/fkeff*fcs 
                elseif(nza.eq.962450)then
                    if(rd1.eq."tot-cap") then
                        fcm5=fcm5-fcs
                        write(2,*)nza,fcs
                    endif
                    if(rd1.eq."nu-sigf") fcm5=fcm5+fnonleak/fkeff*fcs 
                elseif(nza.eq.942380)then
                    if(rd1.eq."tot-cap") then
                        write(2,*)nza,fcs
                    endif
                elseif(nza.eq.932380)then
                    if(rd1.eq."tot-cap") then
                        write(2,*)nza,fcs
                    endif
                elseif(nza.eq.962420)then
                    if(rd1.eq."tot-cap") then
                        write(2,*)nza,fcs
                    endif
                elseif(nza.eq.952420)then
                    if(rd1.eq."tot-cap") then
                        write(2,*)nza,fcs
                    endif
                elseif(nza.eq.952421)then
                    if(rd1.eq."tot-cap") then
                        write(2,*)nza,fcs
                    endif
                else
                    if(rd1.eq."tot-cap") fa=fa+FN(i)*fcs
                    if(rd1.eq."nu-sigf") ff=ff+fnonleak/fkeff*FN(i)*fcs 
                endif
                exit                
            endif
            i=i+1
        end do       
    end do
    close(1)
    close(2)

!----------------------考虑(n,2n),(n,3n)反应道对中子产生率的贡献--------------------------
!----------------------对每一个添料核素，系数上增加两项，对于其他核素，则在常数项上增加---------------------
    open(3,file='NN_cs.dat')
        do while(.true.)
        read(3,*,iostat=nn)nza,fcs,mt
        if(nn/=0)exit
        i=1
        do while(.true.)
            if(i.gt.429)exit
            if(N(i).eq.nza)then
                if(nza.eq.922330) then
                    if(mt.eq.16) fu3=fu3+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fu3=fu3+fnonleak/fkeff*fcs*3
                elseif(nza.eq.902320)then
                    if(mt.eq.16) fth=fth+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fth=fth+fnonleak/fkeff*fcs*3
                elseif(nza.eq.932370)then
                    if(mt.eq.16) fnp=fnp+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fnp=fnp+fnonleak/fkeff*fcs*3
                elseif(nza.eq.952410)then
                    if(mt.eq.16) fam1=fam1+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fam1=fam1+fnonleak/fkeff*fcs*3 
                elseif(nza.eq.952430)then
                    if(mt.eq.16) fam3=fam3+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fam3=fam3+fnonleak/fkeff*fcs*3  
                elseif(nza.eq.962440)then
                    if(mt.eq.16) fcm4=fcm4+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fcm4=fcm4+fnonleak/fkeff*fcs*3 
                elseif(nza.eq.962450)then
                    if(mt.eq.16) fcm5=fcm5+fnonleak/fkeff*fcs*2
                    if(mt.eq.17) fcm5=fcm5+fnonleak/fkeff*fcs*3
                else
                    if(mt.eq.16) ff=ff+fnonleak/fkeff*FN(i)*fcs*2
                    if(mt.eq.17) ff=ff+fnonleak/fkeff*FN(i)*fcs*3 
                endif
                exit                
            endif
            i=i+1
        end do       
    end do
    close(3)
!-----------------------------------------------------------------------------------------
    
    fcon=ff-fa
    
!!!!!!!!!!!!!!!!!!!!! 要把中子得失守恒方程的系数保存到文件中，matrix_NCEs.dat(neutron conservation equations) !!!!!!!!!!!!!!!!!!!!!!!!
    open(1,file='matrix_NCEs.dat')
    write(1,*)"co_th2  co_u3  co_np7  co_am1  "
    write(1,*)fth,fu3,fnp,fam1
    write(1,*)"co_am3  co_cm4  co_cm5  " !!!!!!co表示“coefficient”，方程组系数
    write(1,*)fam3,fcm4,fcm5
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,0.0
    write(1,*)"co_const  "
    write(1,*)-fcon
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!! 计算反推添料率方程相关系数 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------计算各添料率前面系数-------------   
    feed_Th=-1.0E+16
    feed_U3=-1.0E+16
    feed_MA=-1.0E+16
    open(1,file='Proportion_MA.inp')
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn)nza,fna,fpro
        if(nn/=0)exit
        if(nza.eq.93237)feed_Np=feed_MA*fpro
        if(nza.eq.95241)feed_Am1=feed_MA*fpro
        if(nza.eq.95243)feed_Am3=feed_MA*fpro
        if(nza.eq.96244)feed_Cm4=feed_MA*fpro
        if(nza.eq.96245)feed_Cm5=feed_MA*fpro
    enddo
    close(1)
!------------------------------------------------------    
!!!!!!!!!Th232反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Th.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
!!!!!!!!!Th232消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Th.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
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
!!!!!!! 将Th232反推添料率方程的各项保存到文件feed_equation_Th.dat!!!!!!!!!!!!!
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
    
!!!!!!!!!U233反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_U3.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
!!!!!!!!!U233消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_U3.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
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
!!!!!!! 将U233反推添料率方程的各项保存到文件feed_equation_U3.dat!!!!!!!!!!!!!
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
    
    
!-----------------------------------------------------------Pu238 equilibrium burnup equation-----------------------------------------------------
!-------------------------------------Pu-238's two main precusors are Np-238 and Cm-242---------------------------------------------------
!-------------------------------------Np-238 is produced mainly by captruing a neutron by Np-237-----------------------------------------
!-------------------------------------Cm-242 is 
!---------------Am2421消失项系数-----------------------------
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Am21.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.952421)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!--------------------------------------------------------------------
!--------------Am2421产生率----------------------------------
    fp1=0.0
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.952410).and.(mt.eq.1021))then
           fam11=flux*fpcs
        endif
        
        if((nza.eq.952430).and.(mt.eq.161))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp1=fp1+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.952420).and.(mt.eq.511))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp1=fp1+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.962420).and.(mt.eq.1031))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp1=fp1+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.972450).and.(mt.eq.1071))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp1=fp1+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)    
    
    fp1=fp1/fcl
    fam11=fam11/fcl !---fam11代表核素Am241前面的系数
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    
!---------------Am2420消失项系数-----------------------------
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Am20.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.952420)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!--------------------------------------------------------------------
!--------------Am2420产生率----------------------------------
    fp2=0.0
    open(1,file='DECAY_Am20.dat')
    read(1,*)
    read(1,*)
    read(1,*)nza,fthalf,fprob
    fam21=log(2.0)/fthalf*fprob*(1.0E+24)
    close(1)
    
    
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.952410).and.(mt.eq.1020))then
           fam12=flux*fpcs
        endif
        
        if((nza.eq.952430).and.(mt.eq.160))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp2=fp2+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.952440).and.(mt.eq.170))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp2=fp2+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.952441).and.(mt.eq.170))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp2=fp2+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.962420).and.(mt.eq.1030))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp2=fp2+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.972450).and.(mt.eq.1070))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp2=fp2+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)    
    
    fp2=fp2/fcl
    fam21=fam21/fcl
    fam12=fam12/fcl !---fam1代表核素Am241前面的系数
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    fam1=fam12+fam21*fam11
    fp=fp2+fam21*fp1
!------------------------------------------------------------------------------------------
    
!---------------Cm242消失项系数-----------------------------
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Cm2.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.962420)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!--------------------------------------------------------------------
!--------------Cm242产生率----------------------------------
    fp3=0.0
    open(1,file='DECAY_Cm2.dat')
    read(1,*)
    read(1,*)
    read(1,*)nza,fthalf,fprob
    fam20=log(2.0)/fthalf*fprob*(1.0E+24)
    close(1)
    
    
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.962410).and.(mt.eq.102))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp3=fp3+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.962430).and.(mt.eq.16))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp3=fp3+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.962440).and.(mt.eq.17))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp3=fp3+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)    
    
    fp3=fp3/fcl
    fam20=fam20/fcl
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    fp=fp*fam20+fp3
    fam1=fam1*fam20
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    
    
!---------------Np238消失项系数-----------------------------
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Np8.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.932380)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!--------------------------------------------------------------------
!--------------Np238产生率----------------------------------
    fp4=0.0
    open(1,file='DECAY_Np8.dat')
    read(1,*)
    read(1,*)
    read(1,*)nza,fthalf,fprob
    close(1)
    i=1
    do while(.true.)
        if(i.gt.429)exit
        if(N(i).eq.nza)then
            fp4=fp4+FN(i)*(1.0E+24)*(log(2.0)/fthalf)*fprob
            exit
        endif
        i=i+1
    end do
    
    
    
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.932370).and.(mt.eq.102))then
           fnp7=flux*fpcs 
        endif
        if((nza.eq.932390).and.(mt.eq.16))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp4=fp4+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.942380).and.(mt.eq.103))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp4=fp4+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.952410).and.(mt.eq.107))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fp4=fp4+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)    
    
    fp4=fp4/fcl
    fnp7=fnp7/fcl
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------    
    
!---------------Pu238消失项系数-----------------------------
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Pu8.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.942380)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!--------------------------------------------------------------------
!--------------Pu238产生率----------------------------------
    fpr=0.0
    open(1,file='DECAY_Pu8.dat')
    read(1,*)
    read(1,*)
    read(1,*)nza,fthalf,fprob
    fnp8=log(2.0)/fthalf*fprob*(1.0E+24)
    read(1,*)nza,fthalf,fprob
    fcm2=log(2.0)/fthalf*fprob*(1.0E+24)    
    close(1)
    
    
    open(1,file='prod_cs_feed.dat')
    do while(.true.)
        read(1,*,iostat=nn)nza,fpcs,mt
        if(nn/=0)exit
        if((nza.eq.942370).and.(mt.eq.102))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fpr=fpr+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do  
        endif
        if((nza.eq.942390).and.(mt.eq.16))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fpr=fpr+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.942400).and.(mt.eq.17))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fpr=fpr+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.942410).and.(mt.eq.37))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fpr=fpr+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
        if((nza.eq.962410).and.(mt.eq.107))then
           i=1
            do while(.true.)
                if(i.gt.429)exit
                if(N(i).eq.nza)then
                    fpr=fpr+FN(i)*flux*fpcs
                    exit
                endif
                i=i+1
            end do 
        endif
    end do
    close(1)    
    
    fpr=fpr/fcl
    fnp8=fnp8/fcl
    fcm2=fcm2/fcl
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    fam1=fam1*fcm2
    fnp7=fnp7*fnp8
    fpr=fpr+fp*fcm2+fp4*fnp8
!-----------------------------------------------------------------------------------------
!-------------------- 将Pu238核素浓度保存到文件N_Pu238.dat--------------
    open(1,file='N_Pu238.dat') 
    write(1,*)"np7 am1"
    write(1,*)fnp7,fam1
    write(1,*)"const"
    write(1,*)fpr
    close(1)
!----------------------------------------------------------------------------------------  
    
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!Np237反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Np.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
        if((nza.eq.932360).and.(mt.eq.102))then
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
        
        if((nza.eq.932380).and.(mt.eq.16))then
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
        if((nza.eq.932390).and.(mt.eq.17))then
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
        if((nza.eq.942370).and.(mt.eq.103))then
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
        if((nza.eq.952400).and.(mt.eq.107))then
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
!!!!!!!!!Np237消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Np.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.932370)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 将Np237反推添料率方程的各项保存到文件feed_equation_Np.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_Np.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)0.0,0.0,fcl,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,0.0,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,feed_Np
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Am241反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Am1.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
        if((nza.eq.962410).and.(mt.eq.103))then
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
        if((nza.eq.952400).and.(mt.eq.102))then
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
        
        if((nza.eq.952420).and.(mt.eq.16))then
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
        if((nza.eq.952430).and.(mt.eq.17))then
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
        if((nza.eq.952421).and.(mt.eq.16))then
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
        if((nza.eq.952440).and.(mt.eq.37))then
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
!!!!!!!!!Am241消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Am1.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.952410)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 将Am241反推添料率方程的各项保存到文件feed_equation_Am1.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_Am1.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)0.0,0.0,0.0,fcl
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,0.0,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,feed_Am1
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Am243反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Am3.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
        if((nza.eq.962430).and.(mt.eq.103))then
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
        if((nza.eq.952420).and.(mt.eq.102))then
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
        
        if((nza.eq.952421).and.(mt.eq.102))then
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
        if((nza.eq.972460).and.(mt.eq.107))then
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
        if((nza.eq.952440).and.(mt.eq.16))then
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
!!!!!!!!!Am243消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Am3.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.952430)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 将Am243反推添料率方程的各项保存到文件feed_equation_Am3.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_Am3.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)0.0,0.0,0.0,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)fcl,0.0,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,feed_Am3
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!Cm244反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Cm4.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
        if((nza.eq.962430).and.(mt.eq.102))then
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
        if((nza.eq.962450).and.(mt.eq.16))then
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
        
        if((nza.eq.962460).and.(mt.eq.17))then
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
        if((nza.eq.962470).and.(mt.eq.37))then
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
!!!!!!!!!Cm244消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Cm4.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.962440)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 将Cm244反推添料率方程的各项保存到文件feed_equation_Cm4.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_Cm4.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)0.0,0.0,0.0,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,fcl,0.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,feed_Cm4
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!Cm245反推添料率方程系数!!!!!!!!!!!!!!!!!!!
!!!!!!!!!产生率fp!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fp=0.0
    open(1,file='DECAY_Cm5.dat')
    read(1,*)
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn) nza,fthalf,fprob !!!!!!!!!fthalf为核素半衰期，fprob为核素衰变分支比
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
        if((nza.eq.962440).and.(mt.eq.102))then
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
        if((nza.eq.962460).and.(mt.eq.16))then
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
        
        if((nza.eq.962470).and.(mt.eq.17))then
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
        if((nza.eq.962480).and.(mt.eq.37))then
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
        if((nza.eq.972450).and.(mt.eq.103))then
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
        if((nza.eq.982480).and.(mt.eq.107))then
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
!!!!!!!!!Cm245消失项系数!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fcl=0.0 !!!!!!coefficient loss消失项系数!!!!!!!!!
    open(2,file='DECAY_Cm5.dat')
    read(2,*)
    read(2,*)nza,fthalf,fprob
    fcl=fcl+log(2.0)/fthalf*fprob*(1.0E+24)
    close(2)
    
    open(2,file='loss_cs_feed.dat')
    do while(.true.)
        read(2,*,iostat=nn)nza,flcs
        if(nn/=0)exit
        if(nza.eq.962450)then
            fcl=fcl+flux*flcs
            exit
        endif
    end do
    close(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! 将Cm245反推添料率方程的各项保存到文件feed_equation_Cm5.dat!!!!!!!!!!!!!
    open(1,file='feed_equation_Cm5.dat') 
    write(1,*)"th2 u3 np7 am1"
    write(1,*)0.0,0.0,0.0,0.0
    write(1,*)"am3 cm4 cm5"
    write(1,*)0.0,0.0,fcl
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,feed_Cm5
    write(1,*)"const"
    write(1,*)fp
    close(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call extra_condition
    
    call solve(ff_th,ff_u3,ff_ma)
    
    open(1,file='M_sum_HM.dat')
	read(1,*)
	read(1,*)sumhm
	close(1)
    
    FEED(1)=ff_th*(1.0e-8)*(1.0e6/sumhm)/0.602214199
    FEED(2)=ff_u3*(1.0e-8)*(1.0e6/sumhm)/0.602214199
    feedma=ff_ma*(1.0e-8)*(1.0e6/sumhm)/0.602214199
    
    if(feedma.lt.0.0)then
        open(1,file='feed_ma.dat')
        read(1,*)feedma2
        close(1)
        feedma=feedma2/2.0
    endif
    
        open(1,file='feed_ma.dat')
        write(1,*)feedma
        close(1)
    
    
    open(1,file='Proportion_MA.inp')
    read(1,*)
    do while(.true.)
        read(1,*,iostat=nn)nza,fna,fpro
        if(nn/=0)exit
        if(nza.eq.93237)FEED(3)=feedma*fpro
        if(nza.eq.95241)FEED(4)=feedma*fpro
        if(nza.eq.95243)FEED(5)=feedma*fpro
        if(nza.eq.96244)FEED(6)=feedma*fpro
        if(nza.eq.96245)FEED(7)=feedma*fpro
    enddo
    close(1)    
100 format(a)
end subroutine KEY

!---------------TRU总mol浓度不变------------------
!---------------总重金属质量不变/Th质量不变-------
!---------------注意包含Pu238，用Np237和Am241来代表--------
    subroutine extra_condition

    character rd
    
    open(1,file='com_1_x.dat')
    read(1,*)
    read(1,*)
    sumhm1=0.0
    sumtru1=0.0
    do while(.true.)
        read(1,*,iostat=nn)rd,nza,fna,fatoms,frmass,fmass
        if(nn/=0)exit
        if(ANINT(fna).ge.220)then
            if(nza.eq.932370)then
                fco_np=(fna*1.0)/0.602214199
            else if(nza.eq.942380)then
                fco_pu8=fna/0.602214199
            else if(nza.eq.952410)then
                fco_am1=(fna*1.0)/0.602214199
            else if(nza.eq.952430)then
                fco_am3=(fna*1.0)/0.602214199
            else if(nza.eq.962440)then
                fco_cm4=(fna*1.0)/0.602214199
            else if(nza.eq.962450)then
                fco_cm5=(fna*1.0)/0.602214199
            else if(nza.eq.902320)then
                fco_th=(fna*1.0)/0.602214199
            else if(nza.eq.922330)then
                fco_u3=(fna*1.0)/0.602214199
            else
                sumhm1=sumhm1+frmass
            endif
        endif
        
        if(nza.ge.930000)then
            if(nza.eq.932370)then
            else if(nza.eq.942380)then
            else if(nza.eq.952410)then
            else if(nza.eq.952430)then
            else if(nza.eq.962440)then
            else if(nza.eq.962450)then
            else
                sumtru1=sumtru1+frmass/(fna*1.0)*0.602214199
            endif
        endif
    end do    
    close(1)
    
    open(1,file='N_Pu238.dat')
    read(1,*)
    read(1,*)fnp7,fam1
    read(1,*)
    read(1,*)fp
    close(1)
    
    open(1,file='M_sum_TRU1.dat')
    read(1,*)
    read(1,*)sumtru0
    close(1)
    
    open(1,file='M_sum_HM.dat')
    read(1,*)
    read(1,*)sumhm0
    close(1)
    
    open(1,file='Equation_TRUconst.dat')
    write(1,*)"th u3 np am1"
    write(1,*)0.0,0.0,1.0+fnp7,1.0+fam1
    write(1,*)"am3 cm4 cm5"
    write(1,*)1.0,1.0,1.0
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,0.0
    write(1,*)"const"
    write(1,*)sumtru0-sumtru1-fp
    close(1)
    
    open(1,file='Equation_HMconst.dat')
    write(1,*)"th u3 np am1"
    write(1,*)fco_th,fco_u3,fco_np+fco_pu8*fnp7,fco_am1+fco_pu8*fam1
    write(1,*)"am3 cm4 cm5"
    write(1,*)fco_am3,fco_cm4,fco_cm5
    write(1,*)"feedTh feedU3 feedMA"
    write(1,*)0.0,0.0,0.0
    write(1,*)"const"
    write(1,*)sumhm0-sumhm1-fco_pu8*fp
    close(1)
    end subroutine extra_condition
       
!-------------------用来求解方程组的子程序-----------------------
!-------------------求解的变量数：3个添料率，7个核素浓度
!   方程组顺序：    中子得失守恒方程       matrix_NCEs.dat
!                   重金属质量守恒方程     Equation_HMconst.dat
!                   TRU摩尔浓度守恒方程    Equation_TRUconst.dat
!                   Th反推添料率方程       feed_equation_Th.dat
!                   U3反推添料率方程       feed_equation_U3.dat
!                   Np反推添料率方程       feed_equation_Np.dat
!                   Am1反推添料率方程      feed_equation_Am1.dat
!                   Am3反推添料率方程      feed_equation_Am3.dat
!                   Cm4反推添料率方程      feed_equation_Cm4.dat
!                   Cm5反推添料率方程      feed_equation_Cm5.dat
!----------------------------------------------------------------
    
    subroutine solve(ff_th,ff_u3,ff_ma)
    integer::i,k,j
    real*8::A(10,10),FN(3,3),b(10),Fb(3),FNEW(7),FOLD(7)    
    open(1,file='matrix_NCEs.dat')
    read(1,*)
    read(1,*)(A(1,j),j=1,4)
    read(1,*)
    read(1,*)(A(1,j),j=5,7)
    read(1,*)
    read(1,*)(A(1,j),j=8,10)
    read(1,*)
    read(1,*)b(1)
    close(1)
    
    open(1,file='Equation_HMconst.dat')
    read(1,*)
    read(1,*)(A(2,j),j=1,4)
    read(1,*)
    read(1,*)(A(2,j),j=5,7)
    read(1,*)
    read(1,*)(A(2,j),j=8,10)
    read(1,*)
    read(1,*)b(2)
    close(1)
    
    open(1,file='Equation_TRUconst.dat')
    read(1,*)
    read(1,*)(A(3,j),j=1,4)
    read(1,*)
    read(1,*)(A(3,j),j=5,7)
    read(1,*)
    read(1,*)(A(3,j),j=8,10)
    read(1,*)
    read(1,*)b(3)
    close(1)
    
    open(1,file='feed_equation_Th.dat')
    read(1,*)
    read(1,*)(A(4,j),j=1,4)
    read(1,*)
    read(1,*)(A(4,j),j=5,7)
    read(1,*)
    read(1,*)(A(4,j),j=8,10)
    read(1,*)
    read(1,*)b(4)
    close(1)
    
    open(1,file='feed_equation_U3.dat')
    read(1,*)
    read(1,*)(A(5,j),j=1,4)
    read(1,*)
    read(1,*)(A(5,j),j=5,7)
    read(1,*)
    read(1,*)(A(5,j),j=8,10)
    read(1,*)
    read(1,*)b(5)
    close(1)
    
    open(1,file='feed_equation_Np.dat')
    read(1,*)
    read(1,*)(A(6,j),j=1,4)
    read(1,*)
    read(1,*)(A(6,j),j=5,7)
    read(1,*)
    read(1,*)(A(6,j),j=8,10)
    read(1,*)
    read(1,*)b(6)
    close(1)
    
    open(1,file='feed_equation_Am1.dat')
    read(1,*)
    read(1,*)(A(7,j),j=1,4)
    read(1,*)
    read(1,*)(A(7,j),j=5,7)
    read(1,*)
    read(1,*)(A(7,j),j=8,10)
    read(1,*)
    read(1,*)b(7)
    close(1)
    
    open(1,file='feed_equation_Am3.dat')
    read(1,*)
    read(1,*)(A(8,j),j=1,4)
    read(1,*)
    read(1,*)(A(8,j),j=5,7)
    read(1,*)
    read(1,*)(A(8,j),j=8,10)
    read(1,*)
    read(1,*)b(8)
    close(1)
    
    open(1,file='feed_equation_Cm4.dat')
    read(1,*)
    read(1,*)(A(9,j),j=1,4)
    read(1,*)
    read(1,*)(A(9,j),j=5,7)
    read(1,*)
    read(1,*)(A(9,j),j=8,10)
    read(1,*)
    read(1,*)b(9)
    close(1)
    
    open(1,file='feed_equation_Cm5.dat')
    read(1,*)
    read(1,*)(A(10,j),j=1,4)
    read(1,*)
    read(1,*)(A(10,j),j=5,7)
    read(1,*)
    read(1,*)(A(10,j),j=8,10)
    read(1,*)
    read(1,*)b(10)
    close(1)
    
    FN=0
    Fb=0
    
    FN(1,1)=A(1,1)*A(4,8)/A(4,1)
    FN(1,2)=A(1,2)*A(5,9)/A(5,2)
    FN(2,1)=A(2,1)*A(4,8)/A(4,1)
    FN(2,2)=A(2,2)*A(5,9)/A(5,2)
    do mm=1,3
        do nn=3,7
            FN(mm,3)=FN(mm,3)+A(mm,nn)*A(nn+3,10)/A(nn+3,nn)
            Fb(mm)=Fb(mm)+A(mm,nn)*b(nn+3)/A(nn+3,nn)
        enddo
        Fb(mm)=Fb(mm)-b(mm)
    enddo
    Fb(1)=Fb(1)+A(1,1)*b(4)/A(4,1)+A(1,2)*b(5)/A(5,2)
    Fb(2)=Fb(2)+A(2,1)*b(4)/A(4,1)+A(2,2)*b(5)/A(5,2)
    ff_ma=Fb(3)/FN(3,3)
    
    ff_u3up=((Fb(1)-FN(1,3)*ff_ma)*FN(2,1)-(Fb(2)-FN(2,3)*ff_ma)*FN(1,1))
    ff_u3down=FN(2,1)*FN(1,2)-FN(1,1)*FN(2,2)
    ff_u3=ff_u3up/ff_u3down
    
    ff_th=(Fb(1)-FN(1,3)*ff_ma-FN(1,2)*ff_u3)/FN(1,1)

    FNEW(1)=b(4)/A(4,1)-A(4,8)/A(4,1)*ff_th
    FNEW(2)=b(5)/A(5,2)-A(5,9)/A(5,2)*ff_u3
    
    do ii=3,7
        FNEW(ii)=b(ii+3)/A(ii+3,ii)-A(ii+3,10)/A(ii+3,ii)*ff_ma
    enddo
    
    open(1,file='com_Ac_old.dat')
    i=1
    do while(.true.)
        read(1,*,iostat=nn)nza,FOLD(i)
        if(nn/=0)exit
        i=i+1
    end do
    close(1)
    
    open(1,file='criteria_NNN.dat')
    do j=1,7
        write(1,*)ABS(FNEW(j)-FOLD(j))/FOLD(j)
    end do
    close(1)
        
    end subroutine solve
    
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
	    write(2,*)" power=",fx," burn=1.0e-16  end"
	    close(2)
	 

        call new_TRITON
        call system('E:\scale6.1\cmds\runscale TRITON.inp')
        call read_kmt
        call read_keffd(0)
        call cs_file(0)

        mmn=1
        NFINAL=0
        do while(.true.)             

            
            call Seek_Equilibrum(mmn)
            
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
                
                !if(NFINAL.gt.3)exit
            endif
            
            mmn=mmn+1
        enddo
        
 
    end program main
    
    
    subroutine read_keffd(ncheck)
	character(120)rd
    
    kn1=mod(ncheck,10000)/1000+48	!将数字转换成字符,做文件名用
	kn2=mod(ncheck,1000)/100+48	!将数字转换成字符,做文件名用
	kn3=mod(ncheck,100)/10+48	!将数字转换成字符,做文件名用
	kn4=mod(ncheck,10)+48		!将数字转换成字符,做文件名用
    
!	ncheck=0
	open(1,file='TRITON.out')
	open(2,file='keff0.dat')

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
    km1=mod(mmm,10000)/1000+48	!将数字转换成字符,做文件名用
	km2=mod(mmm,1000)/100+48	!将数字转换成字符,做文件名用
	km3=mod(mmm,100)/10+48	!将数字转换成字符,做文件名用
	km4=mod(mmm,10)+48		!将数字转换成字符,做文件名用
    
    nnn=mmm-1
    
    kn1=mod(nnn,10000)/1000+48	!将数字转换成字符,做文件名用
	kn2=mod(nnn,1000)/100+48	!将数字转换成字符,做文件名用
	kn3=mod(nnn,100)/10+48	!将数字转换成字符,做文件名用
	kn4=mod(nnn,10)+48		!将数字转换成字符,做文件名用
 
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
    
	kn1=mod(ncycle,10000)/1000+48	!将数字转换成字符,做文件名用
	kn2=mod(ncycle,1000)/100+48	!将数字转换成字符,做文件名用
	kn3=mod(ncycle,100)/10+48	!将数字转换成字符,做文件名用
	kn4=mod(ncycle,10)+48		!将数字转换成字符,做文件名用    
    
        open(2,file='TRITON.out') 
        open(5,file='neutron_flux.dat')!!!!中子通量密度文件
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
        open(3,file='non-leakge.dat')!-----反应堆不泄露概率
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
    open(3,file='tot_cs.dat')!!!!!用于中子得失守恒方程的截面信息：tot-cap、nu-sigf
    !open(6,file='tot_cs_b.dat')
    !open(7,file='tot_cs_ni.dat')
    
    open(5,file='NN_cs.dat')!-----用于中子得失守恒方程的(n,2n),(n,3n)反应道
    
    open(4,file='prod_cs_feed.dat')!!!!反推添料率方程必需的截面信息  

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
                write(5,*)rdd(20:28),rdd(40:50),rdd(60:63)
            endif
            
!-------------------提取Th232的mt=102截面信息，用于石墨反应率的计算-------------
            if((rdd(22:27).eq."902320").and.(rdd(57:63).eq."mt= 102"))then
                write(4,*)rdd(22:28),rdd(40:50),rdd(60:63)
            endif            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!提取所有反推添料率必需的截面信息，无论添料条件!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Th232反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."902310").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."902340").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif        
        if((rdd(22:27).eq."912320").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif        
        if((rdd(22:27).eq."912330").and.(rdd(57:63).eq."mt=  28"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922350").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."902330").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! U233反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."922320").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922340").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922350").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."922360").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942360").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Np237反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."932360").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."932380").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."932390").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942370").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952400").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Am241反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."952400").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952420").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952421").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952430").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952440").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962410").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Am243反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."952420").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952421").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952440").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."972460").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962430").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Cm244反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."962430").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962450").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962460").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962470").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Cm245反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."962440").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962460").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962470").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962480").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."972450").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."982480").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Pu238反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."942370").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942390").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942400").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942410").and.(rdd(57:63).eq."mt=  37"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962410").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Np238反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."932370").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."932390").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."942380").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."952410").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Cm242反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."9623410").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962430").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
        if((rdd(22:27).eq."962440").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),rdd(60:63)
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Am242反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."952410").and.(rdd(32:37).eq."952420").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),1020
        endif
        if((rdd(22:27).eq."952430").and.(rdd(32:37).eq."952420").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),160
        endif
        if((rdd(22:27).eq."952440").and.(rdd(32:37).eq."952420").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),170
        endif
        if((rdd(22:27).eq."952441").and.(rdd(32:37).eq."952420").and.(rdd(57:63).eq."mt=  17"))then
            write(4,*)rdd(20:28),rdd(40:50),370
        endif
        if((rdd(22:27).eq."962420").and.(rdd(32:37).eq."952420").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),1030
        endif
        if((rdd(22:27).eq."972450").and.(rdd(32:37).eq."952420").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),1070
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Am242m反推添料率所需截面信息!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if((rdd(22:27).eq."952410").and.(rdd(32:37).eq."952421").and.(rdd(57:63).eq."mt= 102"))then
            write(4,*)rdd(20:28),rdd(40:50),1021
        endif
        if((rdd(22:27).eq."952430").and.(rdd(32:37).eq."952421").and.(rdd(57:63).eq."mt=  16"))then
            write(4,*)rdd(20:28),rdd(40:50),161
        endif
        if((rdd(22:27).eq."952420").and.(rdd(32:37).eq."952421").and.(rdd(57:63).eq."mt=  51"))then
            write(4,*)rdd(20:28),rdd(40:50),511
        endif
        if((rdd(22:27).eq."962420").and.(rdd(32:37).eq."952421").and.(rdd(57:63).eq."mt= 103"))then
            write(4,*)rdd(20:28),rdd(40:50),1031
        endif
        if((rdd(22:27).eq."972450").and.(rdd(32:37).eq."952421").and.(rdd(57:63).eq."mt= 107"))then
            write(4,*)rdd(20:28),rdd(40:50),1071
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        endif
    end do
    close(3)
    close(4)
    close(5)
!!-----------------------------------------------------------------------------------
!!------B4C用于中子得失守恒方程的核素截面信息
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
!!------Ni合金用于中子得失守恒方程的核素截面信息    
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