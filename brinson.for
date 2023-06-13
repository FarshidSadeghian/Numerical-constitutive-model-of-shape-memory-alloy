	SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPC, SCD, RPL,DDSDDT,
     *                DRPLDE ,DRPLDT , STRAN, DSTRAN, TIME, DTIME, TEMP,
     *                DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS,
     *                NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, 
     *                CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, 
     *                KSTEP, KINC)
        INCLUDE 'ABA_PARAM.INC'
	  CHARACTER*80 CMNAME
	  DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),
     *     	    DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), 
     *          DSTRAN(NTENS), TIME(2), PREDEF(1), DPRED(1),
     *          PROPS(NPROPS), COORDS(3), DROT(3,3), DFGRD0(3,3), 
     *          DFGRD1(3,3)

        real*8 :: e(2),T(2),sigma0,kisi0(2,4),e0,T0,kisi_in(2,4),
     *          property(23),D(4),theta,C(4),Tmat(4),eL(2),sigmacrit(4),
     *          Tol(4),sigma(2),tang_matrix_in(2)
        integer :: M,N,regime(2),i,j,Start_trans(2)
        property(1:23)=PROPS(1:23)
        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)

        kisi_in(1,1:4)=STATEV(1:4) 
        kisi0(1,1:4)=STATEV(5:8)
        T0=STATEV(9)
        e0=STATEV(10)
        sigma0=STATEV(11)
        regime(1)=STATEV(12)
        M=STATEV(13)
        N=STATEV(14)
        sigma(1)=STATEV(15)
        tang_matrix_in(1)=STATEV(16)
        
        T=(/TEMP,TEMP+DTEMP/)
        e=(/STRAN,STRAN+DSTRAN/)
        write(6,*) 'ff',NOEL, NPT
        if (KINC==1 .and. KSTEP==1) then
            tang_matrix_in(1)=D(2)
            kisi0(1,1:4)=PROPS(24:27)
            kisi_in(1,1:4)=PROPS(24:27)        
            if (T0==0) then
                T0=TEMP
            end if
        end if
        call knewtonraphson_asym21(sigma(2),kisi_in(2,1:4),
     *         kisi0(2,1:4),T0,e0,sigma0,tang_matrix_in(2),regime(2),
     *         Start_trans(2),kisi_in(1,1:4),kisi0(1,1:4),T0,
     *         (/T(1),T(2)/),sigma(1),(/e(1),e(2)/),e0,sigma0,
     *         regime(1),M,N,property) 
        if ((Start_trans(2)==1) ) then
            M=1    
            N=regime(2) 
        end if
        if ((regime(2)==N .and. M==1)) then
        else
            M=0 
        end if

        STATEV(1:4)=kisi_in(2,1:4) 
        STATEV(5:8)=kisi0(2,1:4)
        STATEV(9)=T0
        STATEV(10)=e0
        STATEV(11)=sigma0
        STATEV(12)=regime(2)
        STATEV(13)=M
        STATEV(14)=N
        STATEV(15)=sigma(2)
        STATEV(16)=tang_matrix_in(2)    
        do i=1,NTENS
		    STRESS(i)=0.0d0
		    do j=1,NTENS
		        DDSDDE(i,j)=0.0d0
		    end do
	  end do 
        STRESS(1)=sigma(2)*1000000.d0         
        DDSDDE(1,1)=tang_matrix_in(2)*1000000.d0
      end 

      subroutine knewtonraphson_asym21(sigma_corr,kisi_new,kisi0_new,
     *         T0_new,e0_new,sigma0_new,tang_matrix_out,regime_corr,
     *         Start_trans,kisi_en,kisi0_en,T0_en,T,sigma1,e,e0_en,
     *         sigma0_en,regime0,M,N,property_en)   
        implicit none

        real*8, intent (in) :: kisi_en(4),kisi0_en(4),T0_en,T(2),sigma1,
     *                         sigma0_en,e0_en,e(2),property_en(23)
        integer, intent (in) :: regime0,M,N
        real*8, intent (out) :: sigma_corr,kisi_new(4),kisi0_new(4),
     *                          T0_new,sigma0_new,e0_new,tang_matrix_out
        integer,intent (out) :: regime_corr,Start_trans

        real*8 :: kisi(2,4),sigmanew,deltaD,Dm,D0_m,R,property(23),
     *            R0,f_enew,nn,delT,Rp,e0,T0,sigma0,D(4),kisi0(4),
     *           theta,C(4),Tmat(4),eL(2),sigmacrit(4),Tol(4),kisiT(2,4)
        real*8 :: kisi_new1(4),kisip(4),fp,dDm,
     *            x_Newt1,f_Newt1,init_newt,x_Newt2,f_Newt2,x_correct,
     *            x_Bise1,g_Bise1,init_Bise(2),x_Bise2,g_Bise2   
        integer :: M1,N1,regime,i,rr,j,iii,error_Newt1,error_Newt2,
     *             error_Bise1,error_Bise2
 
        property=property_en
        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)

        kisi0=kisi0_en
        kisi(1,1:4)=kisi_en
        e0=e0_en
        T0=T0_en
        sigma0=sigma0_en
        
        sigmanew=sigma1 
        delT=0.1d-6 

        call ksma21_asym1(regime,kisi0_new,T0_new,sigma0_new,e0_new,
     *              kisi_new1,Start_trans,kisi(1,1:4),kisi0,T0,T,
     *          (/sigma1,sigmanew/),sigma0,e0,e(1),regime0,M,N,property)
        if (Start_trans==1 ) then
            M1=1 
            N1=regime 
        end if
        if ((regime==N1 .and. M1==1) .or. (regime==N .and. M==1)) then
            M1=1 
        else
            M1=0 
        end if
        call kisisma21_asym(kisi(2,1:4),regime,sigmanew,T(2),kisi0_new,
     *                 kisi_new1,T0_new,sigma0_new,M1,property) 
        Dm=D(2)+(kisi(2,3)+kisi(2,2))*(D(1)-D(2))+
     *   			kisi(2,4)*(D(3)-D(2))                       
        D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *    			kisi0_new(4)*(D(3)-D(2))        
        R=-eL(1)*Dm*kisi(2,3)+eL(2)*Dm*kisi(2,4) 
        R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4) 
        f_enew=-(-sigmanew+sigma0_new-D0_m*e0_new+R-R0
     *             -theta*(T(2)-T0_new))/Dm-e(2) 
c        Tol(4)=min(0.0000001d0,0.000001d0*abs(f_enew))
c        write (6,*) ,'ccc',Tol(4)
c        property(23)=Tol(4)            
        if (abs(f_enew)>Tol(4)) then
            nn=(e(2)-e(1))*Dm 
        else
            nn=0.0d0
            call kisisma21_asym(kisiT(2,1:4),regime,sigmanew+delT,T(2),
     *                kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)
            kisip=(-kisi(2,1:4)+kisiT(2,1:4))/delT 
            dDm=((kisip(3)+kisip(2))*(D(1)-D(2))+kisip(4)*(D(3)-D(2)))
            Rp=(dDm)*(-eL(1)*kisi(2,3)+eL(2)*kisi(2,4))
     *               -eL(1)*Dm*kisip(3)+eL(2)*Dm*kisip(4)
        end if    
        if (nn/=0.0d0) then
            init_newt=sigma1+nn
            call knewtonraphson(error_Newt1,x_Newt1,f_Newt1,kisi(1,1:4),
     *                kisi0,T0,T,sigma1,sigma0,e0,e,
     *                property,regime0,M,N,init_newt,15)
            if (error_Newt1==0) then
                x_correct=x_Newt1
            else
                init_Bise=(/sigma1,sigma1+nn/)
                call kBisection(error_Bise1,x_Bise1,g_Bise1,kisi(1,1:4),
     *                kisi0,T0,T,sigma1,sigma0,e0,e,
     *                property,regime0,M,N,init_Bise,25)
                if (error_Bise1==0) then
                    x_correct=x_Bise1
                else
                    init_newt=x_Bise1
                    call knewtonraphson(error_Newt2,x_Newt2,f_Newt2,
     *                     kisi(1,1:4),kisi0,T0,T,sigma1,sigma0,e0,e,
     *                     property,regime0,M,N,init_newt,15)
                    if (error_Newt2==0 .or. abs(f_Newt2)<abs(g_Bise1))
     *                         then 
                        x_correct=x_Newt2
                    else
                        x_correct=x_Bise1
                    end if    
                end if
            end if
            if (error_Newt2==1 .and. error_Newt1==1 .and. error_Bise1==1
     *                  ) then
                init_Bise=(/sigma1,sigma1-nn/)
                call kBisection(error_Bise2,x_Bise2,g_Bise2,kisi(1,1:4),
     *                kisi0,T0,T,sigma1,sigma0,e0,e,
     *                property,regime0,M,N,init_Bise,25) 
                if (error_Bise2==0) then
                    x_correct=x_Bise2 
                else
                    if (abs(f_Newt1)<abs(f_Newt2) .and.
     *                      abs(f_Newt1)<abs(g_Bise1) .and.
     *                      abs(f_Newt1)<abs(g_Bise2)) then
                        x_correct=x_Newt1
                    else if (abs(f_Newt2)<abs(g_Bise1) .and.
     *                       abs(f_Newt2)<abs(g_Bise2)) then
                        x_correct=x_Newt2
                    else if (abs(g_Bise1)<abs(g_Bise2)) then
                        x_correct=x_Bise1
                    else
                        x_correct=x_Bise2
                    end if
                end if                          
            end if   
            
            call ksma21_asym1(regime,kisi0_new,T0_new,sigma0_new,
     *                e0_new,kisi_new1,Start_trans,kisi(1,1:4),kisi0,
     *                T0,T,(/sigma1,x_correct/),sigma0,e0,e(1),regime0,
     *                M,N,property) 
            if (Start_trans==1) then
                M1=1 
                N1=regime 
            end if
            if ((regime==N1 .and. M1==1) .or. 
     *              (regime==N .and. M==1))  then
                M1=1 
            else
                M1=0 
            end if
            call kisisma21_asym(kisi(2,1:4),regime,x_correct,T(2),
     *                kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)
            Dm=D(2)+(kisi(2,3)+kisi(2,2))*(D(1)-D(2))+
     *   			  kisi(2,4)*(D(3)-D(2))                       
            D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *      		  kisi0_new(4)*(D(3)-D(2))      
            R=-eL(1)*Dm*kisi(2,3)+eL(2)*Dm*kisi(2,4) 
            R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4)       
            
            regime_corr=regime 
            kisi_new=kisi(2,1:4) 
            sigma_corr=x_correct
            call kisisma21_asym(kisiT(2,1:4),regime,x_correct+delT,T(2),
     *                kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)
            kisip=(-kisi(2,1:4)+kisiT(2,1:4))/delT 
            dDm=((kisip(3)+kisip(2))*(D(1)-D(2))+kisip(4)*(D(3)-D(2)))
            Rp=(dDm)*(-eL(1)*kisi(2,3)+eL(2)*kisi(2,4))
     *               -eL(1)*Dm*kisip(3)+eL(2)*Dm*kisip(4)  
                  
                tang_matrix_out=Dm/(1.0d0-dDm*e(2)-Rp)  
c                tang_matrix_out=Dm                
        else
            regime_corr=regime 
            kisi_new=kisi(2,1:4) 
            sigma_corr=sigma1 
            tang_matrix_out=Dm/(1.0d0-dDm*e(2)-Rp)
        end if
        return
       end subroutine


      subroutine ksma21_asym1(regime_new,kisi0_new,T0_new,sigma0_new,
     *          e0_new,kisi_new,Start_trans,kisi_en,kisi0_1,T0_en,
     *          T_en,sigma_en,sigma0_en,e0_en,e,regime1_en,M,N,property)
        implicit none
        
        real*8,intent (in) :: kisi_en(4),kisi0_1(4),T0_en,T_en(2),
     *                        sigma0_en,e0_en,e,property(23),sigma_en(2)
        integer, intent (in) :: regime1_en,M,N
        real*8, intent (out) :: kisi0_new(4),T0_new,sigma0_new,e0_new,
     *                          kisi_new(4)
        integer,intent (out) :: regime_new,Start_trans       

        real*8 :: kisi0(2,4),Slopex,Slopey,h,sigma_C,T_C,kisi_new1(4),
     *            sigma0,e0,T0,D(4),theta,C(4),Tmat(4),eL(2),sigma(2),
     *            sigmacrit(4),Tol(4),sigma_new2(2),e_new2,kisi(4),
     *            bound(2,4),bound_comp(2,4),kisi_new2(4),T_new2(2),T(2)
        integer :: regime(2),Start_trans1,b,Line_zero,regime1,Line(10)
        
        sigma0=sigma0_en
        T0=T0_en
        e0=e0_en
        regime1=regime1_en
        kisi=kisi_en
        sigma=sigma_en
        T=T_en

        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)

        Line=(/0,0,0,0,0,0,0,0,0,0/)
        kisi0(1,1:4)=kisi0_1 
        kisi0(2,1:4)=kisi0_1
        Start_trans1=0
        Start_trans=0
        
        bound(2,1:4)=(/C(2)*(T(2)-Tmat(4)),C(2)*(T(2)-Tmat(3)),
     *              C(1)*(T(2)-Tmat(2))+sigmacrit(1),
     *              C(1)*(T(2)-Tmat(2))+sigmacrit(2)/) 
        bound(1,1:4)=(/C(2)*(T(1)-Tmat(4)),C(2)*(T(1)-Tmat(3)),
     *              C(1)*(T(1)-Tmat(2))+sigmacrit(1),
     *              C(1)*(T(1)-Tmat(2))+sigmacrit(2)/) 
        bound_comp(2,1:4)=(/C(4)*(T(2)-Tmat(4)),C(4)*(T(2)-Tmat(3)),
     *              C(3)*(T(2)-Tmat(2))+sigmacrit(3),
     *              C(3)*(T(2)-Tmat(2))+sigmacrit(4)/) 
        bound_comp(1,1:4)=(/C(4)*(T(1)-Tmat(4)),C(4)*(T(1)-Tmat(3)),
     *              C(3)*(T(1)-Tmat(2))+sigmacrit(3),
     *              C(3)*(T(1)-Tmat(2))+sigmacrit(4)/)
c        if (sigma(2)==sigma(1)) then
c        else
c        Tol(2)=0.00001d0*(sigma(2)-sigma(1))/abs(sigma(2)-sigma(1))
c        sigma(2)=sigma(2)-Tol(2)
       
c        T(2)=T(2)-0.00001d0*(T(2)-T(1))/(sigma(2)-sigma(1))
c        end if
        Slopex=-(sigma(2)-sigma(1)) 
        Slopey=T(2)-T(1) 
        h=Slopex*T(2)+sigma(2)*Slopey 
        if (bound(2,3)+Tol(2)<sigma(2) .and. 
     *         bound(1,3)+Tol(2)>=sigma(1)) then 
            sigma_C=(-h*C(1)-Slopex*(-C(1)*Tmat(2)+
     *                 sigmacrit(1)))/(-Slopey*C(1)-Slopex) 
            if (sigma_C>sigmacrit(1)) then
                Line(1)=1 
            end if
        end if
        if ((sigmacrit(1)<sigma(2) .and. sigmacrit(1)>=sigma(1)) .or. 
     *        (sigmacrit(1)>sigma(2) .and. sigmacrit(1)<=sigma(1))) then
            T_C=(h-Slopey*sigmacrit(1))/Slopex 
            if (T_C<Tmat(2)) then
                if (T_C>Tmat(1)) then
                    Start_trans1=1
                end if
                Line(2)=1 
            end if
        end if
        if (T(2)<Tmat(2) .and. T(1)>=Tmat(2)) then
            sigma_C=(h-Slopex*Tmat(2))/Slopey 
            if (sigma_C<sigmacrit(1) .and. -sigmacrit(3)<sigma_C) then
                Line(3)=1 
            end if
        end if
        if (T(2)<Tmat(1) .and. T(1)>=Tmat(1)) then
            sigma_C=(h-Slopex*Tmat(1))/Slopey 
            if (sigma_C>sigmacrit(1) .and. sigma_C<sigmacrit(2)) then
                Start_trans1=1 
                Line(4)=1 
            end if
        end if
        if (T(2)>Tmat(3) .and. bound(2,2)>sigma(2) .and.
     *             bound(1,2)<=sigma(1) .and. sigma(1)>=0.0d0)  then
            sigma_C=(-h*C(2)-Slopex*(-C(2)*Tmat(3)))/(-Slopey*C(2)-
     *             Slopex)
            if (sigma_C>=0.0d0) then
                Line(5)=1 
            end if
        end if
        if ((sigma(1)>0.0d0 .and. sigma(2)<=0.0d0 ) .or.
     *         (sigma(1)<0.0d0 .and. sigma(2)>=0.0d0 ) .or.
     *         (sigma(1)==0.0d0 .and. sigma(2)<0.0d0 .and. 
     *         sigma0==0.d0)) then
            Line(6)=1
        end if
        if (-bound_comp(2,3)-Tol(2)>=sigma(2) .and. 
     *              -bound_comp(1,3)-Tol(2)<sigma(1)) then
            sigma_C=(h*C(3)-Slopex*(C(3)*Tmat(2)
     *                  -sigmacrit(3)))/(Slopey*C(3)-Slopex) 
            if (sigma_C<-sigmacrit(3)) then
                Line(7)=1 
            end if
        end if
        if ((-sigmacrit(3)>=sigma(2) .and. -sigmacrit(3)<sigma(1)) .or.
     *          (-sigmacrit(3)<=sigma(2) .and. -sigmacrit(3)>sigma(1)))
     *            then
            T_C=(h+Slopey*sigmacrit(3))/Slopex 
            if (T_C<Tmat(2)) then
                if (T_C>Tmat(1)) then
                    Start_trans1=1 
                end if
                Line(8)=1 
            end if
        end if
        if (T(2)<Tmat(1) .and. T(1)>=Tmat(1)) then
            sigma_C=(h-Slopex*Tmat(1))/Slopey 
            if (sigma_C<-sigmacrit(3) .and. sigma_C>-sigmacrit(4)) then
                Start_trans1=1 
                Line(9)=1 
            end if
        end if
        if (-bound_comp(2,2)<=sigma(2) .and.
     *          -bound_comp(1,2)>sigma(1))  then
            sigma_C=(h*C(4)-Slopex*C(4)*Tmat(3))/(Slopey*C(4)-Slopex) 
            if (sigma_C<0) then
                Line(10)=1 
            end if
        end if

        kisi_new2=kisi
        T_new2=T
        sigma_new2=sigma
        e_new2=e

        b=max(Line(1),Line(2),Line(3),Line(4),Line(5),Line(6),Line(7),
     *              Line(8),Line(9),Line(10))   
        if (b==0) then
        else
            call kfirstpoint_sma21(kisi0(2,1:4),kisi_new1,Line_zero,
     *                 Line,sigma,T,kisi,sigma0,T0,kisi0(1,1:4),e0,
     *                 regime1,M,N,property) 
            T_new2(1)=T0 
            sigma_new2(1)=sigma0 
            kisi_new2=kisi0(2,1:4)          
            e_new2=e0 
        end if

        call ksma21_asym(regime_new,Start_trans,kisi_new2,T0,sigma0,
     *                   T_new2,sigma_new2,regime1,property) 
        T0_new=T0  
        sigma0_new=sigma0 
        e0_new=e0 

        if (Start_trans==1 .or. Start_trans1==1) then
            Start_trans=1 
            kisi0(2,1:4)=kisi_new2
            sigma0_new=sigma_new2(1) 
            e0_new=e_new2
            T0_new=T_new2(1)  
        end if
        kisi0_new=kisi0(2,1:4) 
        kisi_new=kisi_new2
        return
      end subroutine

      subroutine kfirstpoint_sma21(kisi0_new_out,kisi_new1,Line_zero,
     *                    Line_en,sigma,T,kisi_en,sigma0_out,T0_out,
     *                    kisi0_en,e0_out,regime1_out,M,N,property)
        implicit none

        real*8, intent (in) :: sigma(2),T(2),kisi_en(4),kisi0_en(4),
     *                         property(23)
        integer, intent (in) :: Line_en(10),M,N
        real*8,intent (out) :: kisi_new1(4),kisi0_new_out(4)
        integer,intent (out) :: Line_zero
        integer, intent (inout) :: regime1_out
        real*8,intent (inout) :: sigma0_out,T0_out,e0_out        

        real*8 :: Slopex,Slopey,h,S,K(10,4),T1,sigma1,deltaKK(2),Dm,
     *            Dm_C,D0_m,R,R0,e,kisi(10,4),sigma0_new,T0_new,e0_new,
     *        Line(10),D(4),theta,C(4),Tmat(4),eL(2),sigmacrit(4),Tol(4)
        integer :: j,M1,N1,Start_trans,i,regime1_new
        integer :: l,jj,nnn,ii,nn,iiii
        real*8 :: zeromatrix(10,4),b,HH(10,4),FF(10,4),KK(10,4),
     *           kisi0_new(4),del_KK_S,del_KK_T
        
        do iiii=1,10
            zeromatrix(iiii,1:4)=(/0.d0,0.d0,0.d0,0.d0/)            
        end do   
        KK=zeromatrix
        FF=zeromatrix
        HH=zeromatrix
        K=zeromatrix
        kisi=zeromatrix
          
        Line=Line_en
        regime1_new=regime1_out
        sigma0_new=sigma0_out
        T0_new=T0_out
        e0_new=e0_out
        kisi0_new=kisi0_en
        kisi(1,1:4)=kisi_en

        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)

        Slopex=-(sigma(2)-sigma(1))
        Slopey=T(2)-T(1)
        h=Slopex*T(2)+sigma(2)*Slopey
        j=1
        S=(sigma(2)-sigma(1))**2+(T(2)-T(1))**2
        K(1,1:4)=(/sigma(1),T(1),S,0.0d0/)

        if (Line(1)==1) then
            j=j+1
            T1=(Slopey*(-C(1)*Tmat(2)+sigmacrit(1))
     *                  -h)/(-Slopey*C(1)-Slopex)
            sigma1=(-h*C(1)-Slopex*(-C(1)*Tmat(2)
     *                 +sigmacrit(1)))/(-Slopey*C(1)-Slopex)
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,1.d0/)
        end if            
        if (Line(2)==1) then
            j=j+1
            sigma1=sigmacrit(1)
            T1=(h-Slopey*sigma1)/Slopex
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,2.d0/)
        end if
        if (Line(3)==1) then
            j=j+1
            T1=Tmat(2)
            sigma1=(h-Slopex*T1)/Slopey
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,3.d0/)
        end if
        if (Line(4)==1) then
            j=j+1
            T1=Tmat(1)
            sigma1=(h-Slopex*T1)/Slopey
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,4.d0/)
        end if
        if (Line(5)==1) then
            j=j+1
            T1=(Slopey*(-C(2)*Tmat(3))-h)/(-Slopey*C(2)-Slopex)
            sigma1=(-h*C(2)-Slopex*(-C(2)*Tmat(3)))/
     *                 (-Slopey*C(2)-Slopex)
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,5.d0/)
        end if
        if (Line(6)==1) then
            j=j+1
            T1=h/Slopex
            sigma1=0
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,6.d0/)
        end if
        if (Line(7)==1) then
            j=j+1
            T1=(Slopey*(C(3)*Tmat(2)-sigmacrit(3))
     *                    -h)/(Slopey*C(3)-Slopex)
            sigma1=(h*C(3)-Slopex*(C(3)*Tmat(2)
     *                    -sigmacrit(3)))/(Slopey*C(3)-Slopex)
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,7.d0/)
        end if
        if (Line(8)==1) then
            j=j+1
            sigma1=-sigmacrit(3)
            T1=(h-Slopey*sigma1)/Slopex
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,8.d0/)
        end if
        if (Line(9)==1) then
            j=j+1
            T1=Tmat(1)
            sigma1=(h-Slopex*T1)/Slopey
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,9.d0/)
        end if
        if (Line(10)==1) then
            j=j+1
            T1=(Slopey*C(4)*Tmat(3)-h)/(Slopey*C(4)-Slopex)
            sigma1=(h*C(4)-Slopex*C(4)*Tmat(3))/(Slopey*C(4)-Slopex)
            S=(sigma(2)-sigma1)**2+(T(2)-T1)**2
            K(j,1:4)=(/sigma1,T1,S,10.d0/)
        end if
        l=j 
        do jj=1,j
            b=max(K(1,3),K(2,3),K(3,3),K(4,3),K(5,3),K(6,3),K(7,3),
     *                      K(8,3),K(9,3),K(10,3))
            HH=zeromatrix
            nnn=1
            do ii=1,l     
                if (K(ii,3)==b) then
                    FF(jj,1:4)=K(ii,1:4)
                    HH(nnn:l-1,1:4)=K(ii+1:l,1:4)
                    exit            
                else
                    HH(nnn,1:4)=K(ii,1:4)
                    nnn=nnn+1
                end if
            end do
            l=l-1
            K=HH
        end do
        KK=FF

        M1=0
        N1=0
        Line_zero=0
        do i=2,j
            if ((abs(KK(i-1,1)-KK(i,1))<Tol(3)) .and. 
     *                  (abs(KK(i-1,2)-KK(i,2))<Tol(3))) then
            kisi(i,1:4)=kisi(i-1,1:4)  
            else
                deltaKK=(KK(i-1,1:2)-KK(i,1:2))/10.0d0**7
                del_KK_S=KK(i,2)+deltaKK(2)
                del_KK_T=KK(i,1)+deltaKK(1)
                call ksma21_asym(regime1_new,Start_trans,kisi(i-1,1:4),
     *               T0_new,sigma0_new,(/KK(i-1,2),del_KK_S/),
     *               (/KK(i-1,1),del_KK_T/),regime1_new,property)
                if ((Start_trans==1)) then
                    M1=1
                    N1=regime1_new
                end if
                if ((regime1_new==N1 .and. M1==1) .or.
     *                      (regime1_new==N .and. M==1)) then
                    M1=1
                else
                    M1=0
                end if
                call kisisma21_asym(kisi(i,1:4),regime1_new,KK(i,1),
     *                      KK(i,2),kisi0_new,kisi(i-1,1:4),T0_new,
     *                      sigma0_new,M1,property)
            end if
            
            Dm=D(2)+(kisi(2,3)+kisi(2,2))*(D(1)-D(2))+
     *   			kisi(2,4)*(D(3)-D(2))                       
            D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *      			kisi0_new(4)*(D(3)-D(2))          
            R=-eL(1)*Dm*kisi(i,3)+eL(2)*Dm*kisi(i,4)
            R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4)
            e=((KK(i,1)-sigma0_new)+D0_m*e0_new-R+R0
     *                +theta*(KK(i,2)-T0_new))/Dm
                  
            T0_new=KK(i,2)
            kisi0_new=kisi(i,1:4)
            sigma0_new=KK(i,1)
            e0_new=e
            kisi_new1=kisi0_new
        end do
        T0_out=T0_new
        sigma0_out=sigma0_new
        e0_out=e0_new
        regime1_out=regime1_new
        kisi0_new_out=kisi0_new
        return
      end subroutine
 
      subroutine ksma21_asym(regime_new,Start_trans,kisi_en,T0,sigma0,
     *                         T,sigma,regime0,property)
        implicit none
        real*8, intent (in) :: sigma(2),T(2),kisi_en(4),T0,sigma0,
     *                         property(23)
        integer, intent (in) :: regime0
        integer,intent (out) :: regime_new,Start_trans

        real*8 :: vector(2),vector1(2),vector2(2),vector3(2),vector4(2),
     *            vector1_comp(2),vector2_comp(2),vector3_comp(2),
     *            vector4_comp(2),deltaT,deltaS,
     *            direct(4,2),direct_comp(4,2),bound(2,4),
     *            bound_comp(2,4),kisi(4),D(4),theta,C(4),Tmat(4),
     *            eL(2),sigmacrit(4),Tol(4)
        integer :: regime(2),hhh

        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)

        kisi=kisi_en
        regime=(/regime0,0/)
        Start_trans=0
        deltaT=T(2)-T(1)
        deltaS=sigma(2)-sigma(1)
        vector=(/deltaT,deltaS/) 
        vector1=(/-1.0d0,1.0d0/C(1)/) 
        vector2=(/0.0d0,1.0d0/) 
        vector3=(/-1.0d0,0.0d0/) 
        vector4=(/1.0d0,-1.0d0/C(2)/) 
        direct(1:4,2)=(/dot_product(vector,vector1),
     *          dot_product(vector,vector2),dot_product(vector,vector3),
     *          dot_product(vector,vector4)/) 
        vector1_comp=(/-1.0d0,-1.0d0/C(3)/) 
        vector2_comp=(/0.0d0,-1.0d0/) 
        vector4_comp=(/1.0d0,1.0d0/C(4)/) 
        direct_comp(1:4,2)=(/dot_product(vector,vector1_comp),
     *     dot_product(vector,vector2_comp),dot_product(vector,vector3),
     *     dot_product(vector,vector4_comp)/) 

        bound(2,1:4)=(/C(2)*(T(2)-Tmat(4)),C(2)*(T(2)-Tmat(3)),
     *                  C(1)*(T(2)-Tmat(2))+sigmacrit(1),
     *                  C(1)*(T(2)-Tmat(2))+sigmacrit(2)/) 
        bound(1,1:4)=(/C(2)*(T(1)-Tmat(4)),C(2)*(T(1)-Tmat(3)),
     *                  C(1)*(T(1)-Tmat(2))+sigmacrit(1),
     *                  C(1)*(T(1)-Tmat(2))+sigmacrit(2)/) 
        bound_comp(2,1:4)=(/C(4)*(T(2)-Tmat(4)),C(4)*(T(2)-Tmat(3)),
     *                  C(3)*(T(2)-Tmat(2))+sigmacrit(3),
     *                  C(3)*(T(2)-Tmat(2))+sigmacrit(4)/) 
        bound_comp(1,1:4)=(/C(4)*(T(1)-Tmat(4)),C(4)*(T(1)-Tmat(3)),
     *                  C(3)*(T(1)-Tmat(2))+sigmacrit(3),
     *                  C(3)*(T(1)-Tmat(2))+sigmacrit(4)/)

        if ((T(2)>Tmat(2) .and. bound(2,4)<sigma(2) .and. 
     *          bound(1,4)>=sigma(1) .and. sigma(2)>=0.0d0) .or.      
     *         (T(2)<Tmat(2) .and. sigmacrit(2)<sigma(2) .and. 
     *              sigmacrit(2)>=sigma(1)))  then           
            regime(2)=5     
        end if
        if ((T(2)>Tmat(2) .and. -bound_comp(2,4)>sigma(2) .and. 
     *         -bound_comp(1,4)<=sigma(1) .and. sigma(2)<0.0d0) .or.    
     *           (T(2)<Tmat(2) .and. -sigmacrit(4)>sigma(2) .and. 
     *          -sigmacrit(4)<=sigma(1) .and. sigma(2)<0.0d0)) then    
            regime(2)=51     
        end if

        if (T(2)<Tmat(1) .and. T(1)>=Tmat(1) .and. 
     *        sigmacrit(1)>=sigma(2) .and. -sigmacrit(3)<=sigma(2) .and.
     *               kisi(3)/=1.0d0 )  then                
            regime(2)=6     
        end if
        if ((T(2)>Tmat(4) .and. bound(2,1)>sigma(2) .and. 
     *         bound(1,1)<=sigma(1) .and. sigma(2)>=0.0d0) .or.   
     *           (T(2)>Tmat(4) .and. -bound_comp(2,1)<sigma(2) .and. 
     *          -bound_comp(1,1)>=sigma(1) .and. sigma(2)<0.0d0))   then
            regime(2)=7   
        end if

        if ((sigma(2)>0.0d0) .or. (sigma(2)==0.0d0 .and. 
     *              sigma0>=0.0d0)) then
            if (kisi(3)/=1.0d0) then      
                if(T(2)>=Tmat(2) .and. bound(2,3)+Tol(2)<sigma(2) .and.
     *               bound(2,4)+Tol(2)>sigma(2)) then
                    regime(2)=10 
                    if (direct(1,2)>0.0d0) then
                        if (regime(1)==10 .or. regime(1)==41 .or. 
     *                          regime(1)==40)  then                   
                            Start_trans=1   
                        end if
                        regime(2)=11 
                    end if
                end if
                if (T(2)<Tmat(2) ) then
                    if (T(2)>Tmat(1) .and. T(2)<=T0 .and. 
     *                     sigmacrit(2)>=sigma(2))  then
                        regime(2)=30 
                        if (direct(3,2)>0.0d0) then
                            if (regime(1)==30 .or. (regime(1)==30 .and.
     *                           abs(sigma(1))==0.d0)) then
                                Start_trans=1                 
                            end if
                            regime(2)=31 
                        end if
                    end if
                    if (sigmacrit(1)<=sigma(2) .and. 
     *                        sigmacrit(2)>=sigma(2)) then
                        regime(2)=20 
                        if (direct(2,2)>0.0d0) then
                            if (regime(1)==20) then   
                                Start_trans=1                  
                            end if
                            regime(2)=21 
                        end if
                    end if            
                end if
            end if
            hhh=0
            if (kisi(3)/=1.0d0 .and. T(2)>=Tmat(2) .and. 
     *         bound(2,3)<sigma(2) .and. bound(2,4)>sigma(2) .and.
     *      kisi(1)/=0.0d0 .and. T(2)>Tmat(3) .and. bound(2,1)<sigma(2)
     *      .and. bound(2,2)>sigma(2) .and. direct(4,2)<0.d0 ) then
                hhh=1
            end if
            if (kisi(1)/=0.0d0 .and. T(2)>Tmat(3) .and. hhh==0 .and.
     *              bound(2,1)<sigma(2) .and. bound(2,2)>sigma(2)) then 
                regime(2)=40 
                if (direct(4,2)>0.0d0) then
                    if (regime(1)==40 .or. regime(1)==11 .or.
     *                         regime(1)==10 .or. regime(1)==81) then   
                        Start_trans=1                      
                    end if
                    regime(2)=41 
                end if
            end if   
        else    
            if (kisi(4)/=1.0d0)  then   
            if(T(2)>=Tmat(2) .and. bound_comp(2,3)+Tol(2)<abs(sigma(2)) 
     *                .and. bound_comp(2,4)+Tol(2)>abs(sigma(2))) then  
                    regime(2)=80 
                    if (direct_comp(1,2)>0.0d0) then
                        if (regime(1)==80 .or. regime(1)==100 .or. 
     *                          regime(1)==101) then             
                            Start_trans=1                              
                        end if               
                        regime(2)=81 
                    end if
                end if
                if (T(2)<Tmat(2) ) then
                    if (T(2)>Tmat(1) .and. T(2)<=T0 .and.
     *                          sigmacrit(4)>=abs(sigma(2))) then
                        regime(2)=30 
                        if (direct_comp(3,2)>0.0d0) then
                            if (regime(1)==30 .or. (regime(1)==30 .and.
     *                           abs(sigma(1))==0.d0)) then            
                                Start_trans=1                           
                            end if
                            regime(2)=31  
                        end if
                    end if
                    if (sigmacrit(3)<=abs(sigma(2)) .and. 
     *                        sigmacrit(4)>=abs(sigma(2))) then
                        regime(2)=90 
                        if (direct_comp(2,2)>0.0d0) then
                            if (regime(1)==90) then             
                                Start_trans=1            
                            end if
                            regime(2)=91 
                        end if
                    end if     
                end if
            end if
            hhh=0
            if (kisi(4)/=1.0d0 .and. T(2)>=Tmat(2) .and. 
     *          bound_comp(2,3)<abs(sigma(2)) .and. 
     *          bound_comp(2,4)>abs(sigma(2)) .and. kisi(1)/=0.0d0 .and.
     *          T(2)>Tmat(3) .and. bound_comp(2,1)<abs(sigma(2)) .and. 
     *          bound_comp(2,2)>abs(sigma(2)) .and. 
     *          direct_comp(4,2)<0.0d0 ) then
                  hhh=1
            end if
            if (kisi(1)/=0.0d0 .and. T(2)>Tmat(3) .and. hhh==0 .and.
     *                bound_comp(2,1)<abs(sigma(2)) .and. 
     *                bound_comp(2,2)>abs(sigma(2))) then
                regime(2)=100 
                if (direct_comp(4,2)>0.0d0) then
                    if (regime(1)==100 .or. regime(1)==81 .or. 
     *                      regime(1)==80 .or. regime(1)==41) then 
                    Start_trans=1         
                    end if
                    regime(2)=101 
                end if
            end if
        end if
        regime_new=regime(2) 
        return
      end subroutine


      subroutine kisisma21_asym(kisi,m,sigma,T,kisi0,kisi_old,T0,
     *                     sigma0,MM,property)
        implicit none

        real*8, intent (in) :: sigma,T,kisi0(4),kisi_old(4),T0,sigma0,
     *                         property(23)
        integer, intent (in) :: m,MM
        real*8,intent (out) :: kisi(4)

        real*8 :: pi,aA,Am,delta,sigmacrit_new,T0_norm,D(4),theta,C(4),
     *            Tmat(4),eL(2),sigmacrit(4),Tol(4)

        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)

        pi=3.14159265358979323846d0
        aA=pi/(Tmat(4)-Tmat(3))
        aM=pi/(Tmat(2)-Tmat(1))
        kisi=kisi_old
        if (MM==0) then
            if (T>=Tmat(1) .and. T<=T0 .and. sigmacrit(2)>=sigma .and.
     *                   -sigmacrit(4)<=sigma) then  
                delta=(1.0d0-kisi0(2)-kisi0(3)-kisi0(4))
     *                 *(cos(aM*(T-Tmat(1))))/2.0d0+(1.0d0+kisi0(2)
     *                  -kisi0(3)-kisi0(4))/2.0d0    
            else
                delta=kisi0(2)
            end if
            select case (m)
                case (11)                 
                    kisi(3)=(1.0d0-kisi0(3))*cos(pi*(sigma
     *                   -(C(1)*(T-Tmat(2))+sigmacrit(2)))/(sigmacrit(1)
     *                      -sigmacrit(2)))/2.0d0+(1.0d0+kisi0(3))/2.0d0
                    kisi(2)=kisi0(2)*(1.0d0-kisi(3))/(1.0d0-kisi0(3))
                    kisi(4)=kisi0(4)*(1.0d0-kisi(3))/(1.0d0-kisi0(3))
                case (21)              
                    kisi(3)=(1.0d0-kisi0(3))*cos(pi*(sigma
     *                 -sigmacrit(2))/(sigmacrit(1)-sigmacrit(2)))/2.0d0
     *                  +(1.0d0+kisi0(3))/2.0d0
                 kisi(2)=delta-delta*(kisi(3)-kisi0(3))/(1.0d0-kisi0(3))
                    kisi(4)=kisi0(4)*(1.0d0-kisi(3))/(1.0d0-kisi0(3))
                case (31)  
                    kisi(2)=delta
                case (41)     
                    kisi(1)=kisi0(1)*(cos(aA*(T-Tmat(3)-sigma/C(2)))
     *                          +1.0d0)/2.0d0
                    kisi(3)=kisi0(3)*kisi(1)/kisi0(1)
                    kisi(2)=kisi0(2)*kisi(1)/kisi0(1)
                    kisi(4)=kisi0(4)*kisi(1)/kisi0(1)
                case (5)                
                    kisi=(/1.0d0,0.0d0,1.0d0,0.0d0/)    
                case (51)
                    kisi=(/1.0d0,0.0d0,0.0d0,1.0d0/)
                case (6)
                    kisi=(/1.0d0,1.0d0-kisi_old(3)-kisi_old(4),
     *                       kisi_old(3),kisi_old(4)/)
                case (7)
                    kisi=(/0.0d0,0.0d0,0.0d0,0.0d0/)
                case (81)
                    kisi(4)=(1.0d0-kisi0(4))*cos(pi*(abs(sigma)-(C(3)
     *                        *(T-Tmat(2))+sigmacrit(4)))/(sigmacrit(3)
     *                      -sigmacrit(4)))/2.0d0+(1.0d0+kisi0(4))/2.0d0
                    kisi(2)=kisi0(2)*(1.0d0-kisi(4))/(1.0d0-kisi0(4))  
                    kisi(3)=kisi0(3)*(1.0d0-kisi(4))/(1.0d0-kisi0(4))
                case (91)
                    kisi(4)=(1.0d0-kisi0(4))*cos(pi*(abs(sigma)
     *                      -sigmacrit(4))/(sigmacrit(3)
     *                      -sigmacrit(4)))/2.0d0+(1+kisi0(4))/2.0d0
                 kisi(2)=delta-delta*(kisi(4)-kisi0(4))/(1.0d0-kisi0(4))
                    kisi(3)=kisi0(3)*(1.0d0-kisi(4))/(1.0d0-kisi0(4))
                case (101)
                   kisi(1)=kisi0(1)*(cos(aA*(T-Tmat(3)-abs(sigma)/C(4)))
     *                       +1.0d0)/2.0d0
                    kisi(3)=kisi0(3)*kisi(1)/kisi0(1)
                    kisi(2)=kisi0(2)*kisi(1)/kisi0(1)  
                    kisi(4)=kisi0(4)*kisi(1)/kisi0(1)
            end select
        else
            if (T>=Tmat(1) .and. T<=T0 .and. sigmacrit(2)>=sigma .and.
     *                   -sigmacrit(4)<=sigma) then
                delta=(1.0d0-kisi0(2)-kisi0(3)-kisi0(4))*(cos(pi*(T
     *                 -Tmat(1))/(T0-Tmat(1))))/2.0d0+(1.0d0+kisi0(2)
     *                  -kisi0(3)-kisi0(4))/2.0d0    
            else
                delta=kisi0(2)
            end if
            select case (m)
                case (11)  
                    sigmacrit_new=sigma0-C(1)*(T0-Tmat(2))
                    kisi(3)=(1.0d0-kisi0(3))*cos(pi*(sigma-(C(1)*(T
     *                         -Tmat(2))+sigmacrit(2)))/(sigmacrit_new
     *                      -sigmacrit(2)))/2.0d0+(1.0d0+kisi0(3))/2.0d0
                    kisi(2)=kisi0(2)*(1.0d0-kisi(3))/(1.0d0-kisi0(3))
                    kisi(4)=kisi0(4)*(1.0d0-kisi(3))/(1.0d0-kisi0(3))
                case (21)              
                    kisi(3)=(1.0d0-kisi0(3))*cos(pi*(sigma
     *                    -sigmacrit(2))/(sigma0-sigmacrit(2)))/2.0d0
     *                    +(1.0d0+kisi0(3))/2.0d0
                 kisi(2)=delta-delta*(kisi(3)-kisi0(3))/(1.0d0-kisi0(3))
                    kisi(4)=kisi0(4)*(1.0d0-kisi(3))/(1.0d0-kisi0(3))
                case (31)  
                 kisi(2)=delta-delta*(kisi(3)-kisi0(3))/(1.0d0-kisi0(3))
                case (41)        
                    T0_norm=T0-sigma0/C(2)
                    kisi(1)=kisi0(1)*(cos(pi*(T-T0_norm
     *                     -sigma/C(2))/(Tmat(4)-T0_norm))+1.0d0)/2.0d0
                    kisi(3)=kisi0(3)*kisi(1)/kisi0(1)
                    kisi(2)=kisi0(2)*kisi(1)/kisi0(1)
                    kisi(4)=kisi0(4)*kisi(1)/kisi0(1)
                case (5)                
                    kisi=(/1.0d0,0.0d0,1.0d0,0.0d0/)
                case (51)
                    kisi=(/1.0d0,0.0d0,0.0d0,1.0d0/)
                case (6)
                    kisi=(/1.0d0,1.0d0-kisi_old(3)-kisi_old(4),
     *                         kisi_old(3),kisi_old(4)/)
                case (7)
                    kisi=(/0.0d0,0.0d0,0.0d0,0.0d0/)
                case (81)
                    sigmacrit_new=-sigma0-C(1)*(T0-Tmat(2))
                    kisi(4)=(1.0d0-kisi0(4))*cos(pi*(abs(sigma)-(C(3)
     *                     *(T-Tmat(2))+sigmacrit(4)))/(sigmacrit_new
     *                  -sigmacrit(4)))/2.0d0+(1.0d0+kisi0(4))/2.0d0  
                    kisi(2)=kisi0(2)*(1.0d0-kisi(4))/(1.0d0-kisi0(4))
                    kisi(3)=kisi0(3)*(1.0d0-kisi(4))/(1.0d0-kisi0(4))
                case (91)
                    kisi(4)=(1.0d0-kisi0(4))*cos(pi*(abs(sigma)
     *                  -sigmacrit(4))/(abs(sigma0)-sigmacrit(4)))/2.0d0
     *                  +(1.0d0+kisi0(4))/2.0d0
                 kisi(2)=delta-delta*(kisi(4)-kisi0(4))/(1.0d0-kisi0(4))
                    kisi(3)=kisi0(3)*(1.0d0-kisi(4))/(1.0d0-kisi0(4))
                case (101)
                    T0_norm=T0-abs(sigma0)/C(4)
                    kisi(1)=kisi0(1)*(cos(pi*(T-T0_norm-
     *                  abs(sigma)/C(4))/(Tmat(4)-T0_norm))+1.0d0)/2.0d0
                    kisi(3)=kisi0(3)*kisi(1)/kisi0(1)
                    kisi(2)=kisi0(2)*kisi(1)/kisi0(1)
                    kisi(4)=kisi0(4)*kisi(1)/kisi0(1)
            end select
        end if
        kisi(1)=kisi(3)+kisi(2)+kisi(4)
        return
      end subroutine



      subroutine knewtonraphson(error_N,x_final,f_final,kisi1_N,
     *                kisi0_N,T0_N,T_N,sigma1_N,sigma0_N,e0_N,e_N,
     *                property,regime0,M,N,init_x,i_max)
        implicit none

        real*8, intent (in) :: kisi1_N(4),kisi0_N(4),T0_N,T_N(2),init_x,
     *                      sigma1_N,sigma0_N,e0_N,e_N(2),property(23) 
        integer, intent (in) :: regime0,M,N,i_max
        real*8,intent (out) :: x_final,f_final
        integer, intent (out) :: error_N
      
        real*8 :: x(i_max+1),f(i_max+1),kisip(4),Dm,D0_m,dDm,R,R0,Rp,fp,
     *            D(4),theta,C(4),Tmat(4),eL(2),sigmacrit(4),Tol(4),
     *            kisi0_new(4),T0_new,sigma0_new,e0_new,kisi_new1(4),
     *            kisi2_N(4),kisiT_N(1:4),delT
        integer :: i,regime,Start_trans,M1,N1  
        
        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)
     
        M1=-1
        N1=-1
        error_N=1
        x(1)=init_x     
        delT=.1d-7
        
        do i=1,i_max
            f(i)=0.d0
            call ksma21_asym1(regime,kisi0_new,T0_new,sigma0_new,
     *               e0_new,kisi_new1,Start_trans,kisi1_N,kisi0_N,
     *               T0_N,T_N,(/sigma1_N,x(i)/),sigma0_N,e0_N,e_N(1),
     *               regime0,M,N,property) 
            if (Start_trans==1) then
                M1=1 
                N1=regime 
            end if
            if ((regime==N1 .and. M1==1) .or. (regime==N .and. M==1))
     *             then
                M1=1 
            else
                M1=0 
            end if
            call kisisma21_asym(kisi2_N(1:4),regime,x(i),T_N(2),
     *                kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)
            call kisisma21_asym(kisiT_N(1:4),regime,x(i)+delT,T_N(2),
     *                kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)
            kisip=(-kisi2_N(1:4)+kisiT_N(1:4))/delT  
            Dm=D(2)+(kisi2_N(3)+kisi2_N(2))*(D(1)-D(2))+
     *   			kisi2_N(4)*(D(3)-D(2))                       
            D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *      			kisi0_new(4)*(D(3)-D(2))
       		dDm=((kisip(3)+kisip(2))*(D(1)-D(2))+kisip(4)*(D(3)-D(2)))     
            R=-eL(1)*Dm*kisi2_N(3)+eL(2)*Dm*kisi2_N(4) 
            R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4) 
            Rp=(dDm)*(-eL(1)*kisi2_N(3)+eL(2)*kisi2_N(4))
     *                 -eL(1)*Dm*kisip(3)+eL(2)*Dm*kisip(4)
            f(i)=-x(i)+sigma0_new+Dm*e_N(2)-D0_m*e0_new+R-R0
     *                 -theta*(T_N(2)-T0_new)
            fp=-1.0d0+dDm*e_N(2)+Rp             
            if (abs(f(i))<Tol(4)) then
                x_final=x(i)
                f_final=f(i)
                error_N=0
                exit 
            end if   
            x(i+1)=0.d0        
            x(i+1)=x(i)-f(i)/fp 
        end do
        if (i==i_max+1 .and. abs(f(i-1))>Tol(4)) then
            x_final=x(i-1)
            f_final=f(i-1)
            error_N=1
        end if                                        
        return
       end subroutine
         
      subroutine kBisection(error_B,xx_final,g_final,kisi1_N,
     *                kisi0_N,T0_N,T_N,sigma1_N,sigma0_N,e0_N,e_N,
     *                property,regime0,M,N,init_x,j_max)
        implicit none

        real*8, intent (in) :: kisi1_N(4),kisi0_N(4),T0_N,T_N(2),
     *                      sigma1_N,sigma0_N,e0_N,e_N(2),property(23),
     *                      init_x(2)
        integer, intent (in) :: regime0,M,N,j_max
        real*8,intent (out) :: xx_final,g_final
        integer, intent (out) :: error_B
      
        real*8 :: xx(j_max+1),g(j_max+1),kisip(4),Dm,D0_m,dDm,R,R0,Rp,
     *            fp,D(4),theta,C(4),Tmat(4),eL(2),sigmacrit(4),Tol(4),
     *            kisi0_new(4),T0_new,sigma0_new,e0_new,kisi_new1(4),
     *            kisi2_N(4),kisiT_N(1:4),delT,a(2),b(2)
        integer :: j,regime,Start_trans,M1,N1,rr  
        
        D=property(1:4)
        theta=property(5)
        C=property(6:9)
        Tmat=property(10:13)
        eL=property(14:15)
        sigmacrit=property(16:19)
        Tol=property(20:23)
     
        M1=-1
        N1=-1
        error_B=1                
        rr=0  
        xx(1)=init_x(1)
        xx(2)=init_x(2)
                
        call ksma21_asym1(regime,kisi0_new,T0_new,sigma0_new,
     *           e0_new,kisi_new1,Start_trans,kisi1_N,kisi0_N,
     *           T0_N,T_N,(/sigma1_N,xx(1)/),sigma0_N,e0_N,e_N(1),
     *           regime0,M,N,property) 
        if (Start_trans==1) then
            M1=1 
            N1=regime 
        end if
        if ((regime==N1 .and. M1==1) .or. (regime==N .and. M==1))
     *            then
            M1=1 
        else
            M1=0 
        end if
        call kisisma21_asym(kisi2_N(1:4),regime,xx(1),T_N(2),
     *            kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)
        Dm=D(2)+(kisi2_N(3)+kisi2_N(2))*(D(1)-D(2))+
     *            kisi2_N(4)*(D(3)-D(2))                       
        D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *      	  kisi0_new(4)*(D(3)-D(2))   
        R=-eL(1)*Dm*kisi2_N(3)+eL(2)*Dm*kisi2_N(4) 
        R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4) 
        g(1)=-xx(1)+sigma0_new+Dm*e_N(2)-D0_m*e0_new+R-R0
     *           -theta*(T_N(2)-T0_new) 
        
        call ksma21_asym1(regime,kisi0_new,T0_new,sigma0_new,
     *           e0_new,kisi_new1,Start_trans,kisi1_N,kisi0_N,
     *           T0_N,T_N,(/sigma1_N,xx(2)/),sigma0_N,e0_N,e_N(1),
     *           regime0,M,N,property) 
        if (Start_trans==1) then
            M1=1 
            N1=regime 
        end if
        if ((regime==N1 .and. M1==1) .or. (regime==N .and. M==1))
     *             then
            M1=1 
        else
            M1=0 
        end if
        call kisisma21_asym(kisi2_N(1:4),regime,xx(2),T_N(2),
     *            kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)

        Dm=D(2)+(kisi2_N(3)+kisi2_N(2))*(D(1)-D(2))+
     *            kisi2_N(4)*(D(3)-D(2))                       
        D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *  		  kisi0_new(4)*(D(3)-D(2))   
        R=-eL(1)*Dm*kisi2_N(3)+eL(2)*Dm*kisi2_N(4) 
        R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4) 
        g(2)=-xx(2)+sigma0_new+Dm*e_N(2)-D0_m*e0_new+R-R0
     *           -theta*(T_N(2)-T0_new)  
 
        if (g(2)*g(1)<0.0d0) then
            xx(3)=-(xx(2)-xx(1))/2.0d0+xx(2) 
            rr=1 
            a=(/xx(1),g(1)/) 
            b=(/xx(2),g(2)/)
        else
            xx(3)=(xx(2)-xx(1))*2.0d0+xx(2) 
            end if
        do j=3,j_max
            g(j)=0.d0
            xx(j+1)=0.d0
            call ksma21_asym1(regime,kisi0_new,T0_new,sigma0_new,
     *               e0_new,kisi_new1,Start_trans,kisi1_N,kisi0_N,
     *               T0_N,T_N,(/sigma1_N,xx(j)/),sigma0_N,e0_N,e_N(1),
     *               regime0,M,N,property) 
            if (Start_trans==1) then
                M1=1 
                N1=regime 
            end if
            if ((regime==N1 .and. M1==1) .or. (regime==N .and. M==1))
     *             then
                M1=1 
            else
                M1=0 
            end if
            call kisisma21_asym(kisi2_N(1:4),regime,xx(j),T_N(2),
     *                kisi0_new,kisi_new1,T0_new,sigma0_new,M1,property)

            Dm=D(2)+(kisi2_N(3)+kisi2_N(2))*(D(1)-D(2))+
     *   			kisi2_N(4)*(D(3)-D(2))                       
            D0_m=D(2)+(kisi0_new(3)+kisi0_new(2))*(D(1)-D(2))+
     *      			kisi0_new(4)*(D(3)-D(2))    
            R=-eL(1)*Dm*kisi2_N(3)+eL(2)*Dm*kisi2_N(4) 
            R0=-eL(1)*D0_m*kisi0_new(3)+eL(2)*D0_m*kisi0_new(4) 
            g(j)=-xx(j)+sigma0_new+Dm*e_N(2)-D0_m*e0_new+R-R0
     *                -theta*(T_N(2)-T0_new)                    
            if (abs(g(j))<Tol(4)) then
                xx_final=xx(j)
                g_final=g(j)
                error_B=0
                exit           
            end if
            if (rr==1) then
                if (a(2)*g(j)<0.0d0) then
                    xx(j+1)=xx(j)-(xx(j)-a(1))/2.0d0
                    b=(/xx(j),g(j)/)                     
                else
                    xx(j+1)=xx(j)-(xx(j)-b(1))/2.0d0 
                    a=(/xx(j),g(j)/)                      
                end if
            else if (g(j)*g(j-1)<0.0d0) then
                xx(j+1)=-(xx(j)-xx(j-1))/2.0d0+xx(j) 
                rr=1 
                a=(/xx(j-1),g(j-1)/) 
                b=(/xx(j),g(j)/)
            else
                xx(j+1)=(xx(j)-xx(j-1))*2.0d0+xx(j) 
            end if
        end do
        if (j==j_max+1 .and. abs(g(j-1))>Tol(4)) then
            xx_final=xx(j-1)
            g_final=g(j-1)
            error_B=1
        end if                
        return
      end subroutine