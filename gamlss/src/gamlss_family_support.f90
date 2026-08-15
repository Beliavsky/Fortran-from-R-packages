! Family-level CDF and support helpers used by censoring and residual code.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_family_support
   use gamlss_kinds, only : dp
   use gamlss_fit
   use gamlss_continuous
   use gamlss_discrete
   use gamlss_boxcox
   use gamlss_continuous_v02
   use gamlss_discrete_v02
   use gamlss_continuous_v03
   use gamlss_discrete_v03
   use gamlss_flexible_v03
   implicit none
   private
   public :: family_cdf, family_cdf_left, family_is_discrete, family_observation_is_atom
contains

   real(dp) function family_cdf(family,q,a,b,c,d) result(p)
      integer,intent(in)::family
      real(dp),intent(in)::q,a,b,c,d
      select case(family)
      case(GAMLSS_NO);       p=pNO(q,a,b)
      case(GAMLSS_GA);       p=pGA(q,a,b)
      case(GAMLSS_BE);       p=pBE(q,a,b)
      case(GAMLSS_NBI);      p=pNBI(q,a,b)
      case(GAMLSS_NBII);     p=pNBII(q,a,b)
      case(GAMLSS_ZIP);      p=pZIP(q,a,b)
      case(GAMLSS_GG);       p=pGG(q,a,b,c)
      case(GAMLSS_EGB2);     p=pEGB2(q,a,b,c,d)
      case(GAMLSS_GB2);      p=pGB2(q,a,b,c,d)
      case(GAMLSS_JSUO);     p=pJSUo(q,a,b,c,d)
      case(GAMLSS_TF);       p=pTF(q,a,b,c)
      case(GAMLSS_PIG);      p=pPIG(q,a,b)
      case(GAMLSS_BEINF);    p=pBEINF(q,a,b,c,d)
      case(GAMLSS_WEI);      p=pWEI(q,a,b)
      case(GAMLSS_LNO);      p=pLNO(q,a,b,c)
      case(GAMLSS_BCCG);     p=pBCCG(q,a,b,c)
      case(GAMLSS_BCT);      p=pBCT(q,a,b,c,d)
      case(GAMLSS_BCPE);     p=pBCPE(q,a,b,c,d)
      case(GAMLSS_GIG);      p=pGIG(q,a,b,c)
      case(GAMLSS_SHASHO);   p=pSHASHo(q,a,b,c,d)
      case(GAMLSS_SHASH);    p=pSHASH(q,a,b,c,d)
      case(GAMLSS_SIMPLEX);  p=pSIMPLEX(q,a,b)
      case(GAMLSS_SEP);      p=pSEP(q,a,b,c,d)
      case(GAMLSS_SEP1);     p=pSEP1(q,a,b,c,d)
      case(GAMLSS_SEP2);     p=pSEP2(q,a,b,c,d)
      case(GAMLSS_ST1);      p=pST1(q,a,b,c,d)
      case(GAMLSS_ST2);      p=pST2(q,a,b,c,d)
      case(GAMLSS_ST3);      p=pST3(q,a,b,c,d)
      case(GAMLSS_ST4);      p=pST4(q,a,b,c,d)
      case(GAMLSS_ST5);      p=pST5(q,a,b,c,d)
      case(GAMLSS_SEP3);     p=pSEP3(q,a,b,c,d)
      case(GAMLSS_SEP4);     p=pSEP4(q,a,b,c,d)
      case(GAMLSS_NET);      p=pNET(q,a,b,c,d)
      case(GAMLSS_GPO);      p=pGPO(q,a,b)
      case(GAMLSS_DPO);      p=pDPO(q,a,b)
      case(GAMLSS_DEL);      p=pDEL(q,a,b,c)
      case(GAMLSS_SI);       p=pSI(q,a,b,c)
      case(GAMLSS_SICHEL);   p=pSICHEL(q,a,b,c)
      case(GAMLSS_YULE);     p=pYULE(q,a)
      case(GAMLSS_WARING);   p=pWARING(q,a,b)
      case(GAMLSS_ZIPF);     p=pZIPF(q,a)
      case(GAMLSS_ST3C);     p=pST3C(q,a,b,c,d)
      case(GAMLSS_SN1);      p=pSN1(q,a,b,c)
      case(GAMLSS_SN2);      p=pSN2(q,a,b,c)
      case(GAMLSS_SST);      p=pSST(q,a,b,c,d)
      case(GAMLSS_GT);       p=pGT(q,a,b,c,d)
      case(GAMLSS_EXGAUS);   p=pexGAUS(q,a,b,c)
      case(GAMLSS_PARETO);   p=pPARETO(q,a)
      case(GAMLSS_PARETO1);  p=pPARETO1(q,a)
      case(GAMLSS_PARETO2);  p=pPARETO2(q,a,b)
      case(GAMLSS_PARETO2O); p=pPARETO2o(q,a,b)
      case(GAMLSS_PIG2);     p=pPIG2(q,a,b)
      case(GAMLSS_ZIPIG);    p=pZIPIG(q,a,b,c)
      case(GAMLSS_ZAPIG);    p=pZAPIG(q,a,b,c)
      case(GAMLSS_ZISICHEL); p=pZISICHEL(q,a,b,c,d)
      case(GAMLSS_ZASICHEL); p=pZASICHEL(q,a,b,c,d)
      case(GAMLSS_ZIBNB);    p=pZIBNB(q,a,b,c,d)
      case(GAMLSS_ZABNB);    p=pZABNB(q,a,b,c,d)
      case(GAMLSS_ZAZIPF);   p=pZAZIPF(q,a,b)
      case(GAMLSS_GAF);      p=pGAF(q,a,b,c)
      case(GAMLSS_NBF);      p=pNBF(q,a,b,c)
      case(GAMLSS_ZINBF);    p=pZINBF(q,a,b,c,d)
      case default;           p=-1.0_dp
      end select
      if(p>=0.0_dp)p=min(1.0_dp,max(0.0_dp,p))
   end function family_cdf

   real(dp) function family_cdf_left(family,y,a,b,c,d) result(p)
      integer,intent(in) :: family
      real(dp),intent(in) :: y,a,b,c,d
      if(family_is_discrete(family))then
         if(abs(y-real(nint(y),dp))>1.0e-9_dp)then
            p=family_cdf(family,y,a,b,c,d)
         else
            p=family_cdf(family,y-1.0_dp,a,b,c,d)
         end if
      else if(family==GAMLSS_BEINF.and.abs(y)<=1.0e-12_dp)then
         p=0.0_dp
      else if(family==GAMLSS_BEINF.and.abs(y-1.0_dp)<=1.0e-12_dp)then
         p=(c+1.0_dp)/(1.0_dp+c+d)
      else
         p=family_cdf(family,y,a,b,c,d)
      end if
      if(p>=0.0_dp)p=min(1.0_dp,max(0.0_dp,p))
   end function family_cdf_left

   pure logical function family_observation_is_atom(family,y) result(atom)
      integer,intent(in) :: family
      real(dp),intent(in) :: y
      atom=family_is_discrete(family)
      if(family==GAMLSS_BEINF)then
         atom=abs(y)<=1.0e-12_dp.or.abs(y-1.0_dp)<=1.0e-12_dp
      end if
   end function family_observation_is_atom

   pure logical function family_is_discrete(family) result(discrete)
      integer,intent(in)::family
      select case(family)
      case(GAMLSS_NBI,GAMLSS_NBII,GAMLSS_ZIP,GAMLSS_PIG,GAMLSS_GPO,GAMLSS_DPO, &
           GAMLSS_DEL,GAMLSS_SI,GAMLSS_SICHEL,GAMLSS_YULE,GAMLSS_WARING,GAMLSS_ZIPF, &
           GAMLSS_PIG2,GAMLSS_ZIPIG,GAMLSS_ZAPIG,GAMLSS_ZISICHEL,GAMLSS_ZASICHEL, &
           GAMLSS_ZIBNB,GAMLSS_ZABNB,GAMLSS_ZAZIPF,GAMLSS_NBF,GAMLSS_ZINBF)
         discrete=.true.
      case default
         discrete=.false.
      end select
   end function family_is_discrete

end module gamlss_family_support
