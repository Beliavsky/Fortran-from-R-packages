module spatialextremes_partitions
   use spatialextremes_utils, only: bell_number
   implicit none
   private
   public :: list_set_partitions,canonicalize_partition,partition_block_count
contains
   subroutine list_set_partitions(n,part,block_count,info)
      ! Restricted-growth-string representation used by condsimMaxStab.c:
      ! {{1,3},{2,5},{4,6}} -> [0,1,0,2,1,2].
      integer,intent(in)::n
      integer,allocatable,intent(out)::part(:,:),block_count(:)
      integer,intent(out),optional::info
      integer::nb,current
      integer,allocatable::work(:)
      if(n<0)then
      allocate(part(0,0),block_count(0))
      if(present(info))info=1
      return
      end if
      nb=bell_number(n)
      if(nb<=0)then
      allocate(part(0,0),block_count(0))
      if(present(info))info=2
      return
      end if
      allocate(part(n,nb),block_count(nb),work(max(1,n)))
      current=0
      if(n==0)then
         block_count(1)=0
      else
         work=0
         work(1)=0
         call generate(2,0)
      end if
      if(present(info))info=0
   contains
      recursive subroutine generate(pos,mx)
         integer,intent(in)::pos,mx
         integer::label
         if(pos>n)then
            current=current+1
            part(:,current)=work(1:n)
            block_count(current)=mx+1
            return
         end if
         do label=0,mx+1
            work(pos)=label
            call generate(pos+1,max(mx,label))
         end do
      end subroutine generate
   end subroutine list_set_partitions

   pure function canonicalize_partition(partition) result(out)
      integer,intent(in)::partition(:)
      integer::out(size(partition)),seen(size(partition)),i,j,nseen
      nseen=0
      seen=0
      out=0
      do i=1,size(partition)
         j=1
         do while(j<=nseen.and.seen(j)/=partition(i))
         j=j+1
         end do
         if(j>nseen)then
         nseen=nseen+1
         seen(nseen)=partition(i)
         j=nseen
         end if
         out(i)=j-1
      end do
   end function canonicalize_partition

   pure integer function partition_block_count(partition) result(n)
      integer,intent(in)::partition(:)
      integer::c(size(partition))
      if(size(partition)==0)then
      n=0
      return
      end if
      c=canonicalize_partition(partition)
      n=maxval(c)+1
   end function partition_block_count
end module spatialextremes_partitions
