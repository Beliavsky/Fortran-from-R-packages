module rpart_control_api
   use rpart_kinds, only : dp
   use rpart_types, only : rpart_control
   use rpart_utils, only : r_round_even, validate_control
   implicit none
   private
   public :: rpart_make_control

contains

   function rpart_make_control(minsplit,minbucket,cp,maxcompete,maxsurrogate,usesurrogate,xval, &
                               surrogatestyle,maxdepth,stat) result(control)
      integer,intent(in),optional::minsplit,minbucket,maxcompete,maxsurrogate,usesurrogate,xval
      integer,intent(in),optional::surrogatestyle,maxdepth
      real(dp),intent(in),optional::cp
      integer,intent(out),optional::stat
      type(rpart_control)::control
      integer::s

      if(present(minsplit).and.present(minbucket))then
         control%minsplit=minsplit
         control%minbucket=minbucket
      else if(present(minsplit))then
         control%minsplit=minsplit
         control%minbucket=r_round_even(real(minsplit,dp)/3.0_dp)
      else if(present(minbucket))then
         control%minbucket=minbucket
         control%minsplit=3*minbucket
      end if
      if(present(cp))control%cp=cp
      if(present(maxcompete))control%maxcompete=maxcompete
      if(present(maxsurrogate))control%maxsurrogate=maxsurrogate
      if(present(usesurrogate))control%usesurrogate=usesurrogate
      if(present(xval))control%xval=xval
      if(present(surrogatestyle))control%surrogatestyle=surrogatestyle
      if(present(maxdepth))control%maxdepth=maxdepth
      call validate_control(control,s)
      if(present(stat))stat=s
   end function rpart_make_control

end module rpart_control_api
