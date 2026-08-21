module relsurv_parsers
  use relsurv_kinds, only : dp
  use relsurv_ratetable, only : ratetable_type, make_ratetable
  implicit none
  private
  public :: transrate_hmd, transrate_hld

  type :: hld_record
    integer :: year1=0, year2=0, typelt=0, sex=0, age=0
    real(dp) :: ageint=0.0_dp, qx=0.0_dp
    logical :: ageint_ok=.false., qx_ok=.false.
    integer :: file_index=0
  end type hld_record
contains

  function transrate_hmd(male_file,female_file) result(tab)
    character(len=*),intent(in)::male_file,female_file
    type(ratetable_type)::tab
    integer,allocatable::ym(:),yf(:)
    real(dp),allocatable::qm(:),qf(:),men(:,:),women(:,:),cuts(:,:),rates(:)
    integer::ny,na,i,y,a,idx,maxc
    integer,allocatable::dims(:),factor(:),ncuts(:)
    call read_hmd_file(male_file,ym,qm)
    call read_hmd_file(female_file,yf,qf)
    if(size(qm)/=size(qf))error stop 'transrate_hmd: male/female row count differs'
    if(size(ym)/=size(yf).or.any(ym/=yf))error stop 'transrate_hmd: years differ by sex'
    call hmd_shape(ym,ny,na)
    if(na/=111)error stop 'transrate_hmd: each year must contain ages 0..110'
    allocate(men(na,ny),women(na,ny));men=reshape(qm,[na,ny]);women=reshape(qf,[na,ny])
    do y=1,ny
      call hmd_fix_tail(men(:,y));call hmd_fix_tail(women(:,y))
    end do
    allocate(dims(3),factor(3),ncuts(3));dims=[na,ny,2];factor=[0,0,1];ncuts=[na,ny,0]
    maxc=max(na,ny);allocate(cuts(maxc,3));cuts=0.0_dp
    do a=1,na;cuts(a,1)=real(a-1,dp)*365.241_dp;end do
    do y=1,ny;cuts(y,2)=real(days_from_civil(ym((y-1)*na+1),1,1),dp);end do
    allocate(rates(na*ny*2))
    do i=1,2
      do y=1,ny
        do a=1,na
          idx=a+(y-1)*na+(i-1)*na*ny
          if(i==1)rates(idx)=-log(1.0_dp-men(a,y))/365.241_dp
          if(i==2)rates(idx)=-log(1.0_dp-women(a,y))/365.241_dp
        end do
      end do
    end do
    tab=make_ratetable(dims,factor,cuts,ncuts,rates)
  end function transrate_hmd

  function transrate_hld(files,cut_year,race) result(tab)
    character(len=*),intent(in)::files(:)
    integer,intent(in),optional::cut_year(:)
    character(len=*),intent(in),optional::race(:)
    type(ratetable_type)::tab
    type(hld_record),allocatable::rec(:),tmp(:)
    integer,allocatable::years(:),dims(:),factor(:),ncuts(:)
    real(dp),allocatable::cuts(:,:),rates(:),qx(:)
    integer::f,i,k,na,ny,nrace,y,s,a,idx,maxc,nr,yr
    nr=0;allocate(rec(0))
    do f=1,size(files)
      call read_hld_file(files(f),f,tmp)
      call append_records(rec,tmp)
    end do
    if(size(rec)==0)error stop 'transrate_hld: no TypeLT=1 finite-age rows'
    na=maxval(rec%age)+1
    call unique_years(rec,years)
    ny=size(years)
    if(present(cut_year))then
      if(size(cut_year)/=ny)error stop "transrate_hld: length cut_year must match unique Year1"
    end if
    nrace=1
    if(present(race))then
      if(size(race)/=size(files))error stop 'transrate_hld: race/file length mismatch'
      call count_unique_strings(race,nrace)
    end if
    if(nrace==1)then
      allocate(dims(3),factor(3),ncuts(3));dims=[na,ny,2];factor=[0,0,1];ncuts=[na,ny,0]
      maxc=max(na,ny);allocate(cuts(maxc,3));cuts=0.0_dp
    else
      allocate(dims(4),factor(4),ncuts(4));dims=[na,ny,2,nrace];factor=[0,0,1,1];ncuts=[na,ny,0,0]
      maxc=max(na,ny);allocate(cuts(maxc,4));cuts=0.0_dp
    end if
    do a=1,na;cuts(a,1)=real(a-1,dp)*365.241_dp;end do
    do y=1,ny
      yr=years(y);if(present(cut_year))yr=cut_year(y)
      cuts(y,2)=real(days_from_civil(yr,1,1),dp)
    end do
    allocate(rates(product(dims)));rates=0.0_dp
    do k=1,nrace
      do s=1,2
        do y=1,ny
          allocate(qx(na));qx=-1.0_dp
          do i=1,size(rec)
            if(rec(i)%sex/=s.or.rec(i)%year1/=years(y))cycle
            if(present(race))then
              if(race_index(race,rec(i)%file_index)/=k)cycle
            end if
            if(rec(i)%age>=0.and.rec(i)%age<na.and.rec(i)%qx_ok)qx(rec(i)%age+1)=rec(i)%qx
          end do
          call carry_hld_tail(qx)
          do a=1,na
            idx=a+(y-1)*na+(s-1)*na*ny+(k-1)*na*ny*2
            rates(idx)=-log(1.0_dp-qx(a))/365.241_dp
          end do
          deallocate(qx)
        end do
      end do
    end do
    tab=make_ratetable(dims,factor,cuts,ncuts,rates)
  end function transrate_hld

  subroutine read_hmd_file(file,year,qx)
    character(len=*),intent(in)::file
    integer,allocatable,intent(out)::year(:)
    real(dp),allocatable,intent(out)::qx(:)
    character(len=2048)::line
    character(len=128),allocatable::tok(:)
    integer::u,ios,iy,iq,n,yy
    real(dp)::qq
    logical::header
    n=0;allocate(year(0),qx(0));header=.false.;iy=0;iq=0
    open(newunit=u,file=trim(file),status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'transrate_hmd: cannot open file'
    do
      read(u,'(A)',iostat=ios)line;if(ios/=0)exit
      if(len_trim(line)==0)cycle
      call split_ws(line,tok)
      if(.not.header)then
        call header_positions(tok,'Year','qx',iy,iq)
        if(iy>0.and.iq>0)header=.true.
        cycle
      end if
      if(size(tok)<max(iy,iq))cycle
      read(tok(iy),*,iostat=ios)yy;if(ios/=0)cycle
      if(trim(tok(iq))=='.'.or.trim(tok(iq))=='NA')then
        qq=huge(1.0_dp)
      else
        read(tok(iq),*,iostat=ios)qq;if(ios/=0)qq=huge(1.0_dp)
      end if
      call append_int(year,yy);call append_real(qx,qq)
    end do
    close(u)
    if(.not.header.or.size(year)==0)error stop 'transrate_hmd: no readable data'
  end subroutine read_hmd_file

  subroutine read_hld_file(file,file_index,rec)
    character(len=*),intent(in)::file
    integer,intent(in)::file_index
    type(hld_record),allocatable,intent(out)::rec(:)
    character(len=4096)::line
    character(len=256),allocatable::tok(:),head(:)
    integer::u,ios,ic,iy1,iy2,it,is,ia,iai,iq
    type(hld_record)::r
    allocate(rec(0))
    open(newunit=u,file=trim(file),status='old',action='read',iostat=ios)
    if(ios/=0)error stop 'transrate_hld: cannot open file'
    read(u,'(A)',iostat=ios)line;if(ios/=0)error stop 'transrate_hld: empty file'
    call split_csv(line,head);call strip_dots(head)
    ic=find_token(head,'Country');iy1=find_token(head,'Year1');iy2=find_token(head,'Year2')
    it=find_token(head,'TypeLT');is=find_token(head,'Sex');ia=find_token(head,'Age')
    iai=find_token(head,'AgeInt');iq=find_token(head,'qx')
    if(min(iy1,iy2,it,is,ia,iai,iq)<=0)error stop 'transrate_hld: missing required columns'
    do
      read(u,'(A)',iostat=ios)line;if(ios/=0)exit
      call split_csv(line,tok);if(size(tok)<max(iq,max(iai,max(ia,max(is,max(it,max(iy1,iy2)))))))cycle
      read(tok(it),*,iostat=ios)r%typelt;if(ios/=0.or.r%typelt/=1)cycle
      read(tok(iai),*,iostat=ios)r%ageint;r%ageint_ok=(ios==0);if(.not.r%ageint_ok)cycle
      read(tok(iy1),*,iostat=ios)r%year1;if(ios/=0)cycle
      read(tok(iy2),*,iostat=ios)r%year2;if(ios/=0)cycle
      read(tok(is),*,iostat=ios)r%sex;if(ios/=0)cycle
      read(tok(ia),*,iostat=ios)r%age;if(ios/=0)cycle
      r%qx_ok=.false.
      if(trim(tok(iq))/='.'.and.trim(tok(iq))/='NA'.and.len_trim(tok(iq))>0)then
        read(tok(iq),*,iostat=ios)r%qx;r%qx_ok=(ios==0)
      end if
      r%file_index=file_index
      call append_record(rec,r)
    end do
    close(u)
  end subroutine read_hld_file

  subroutine hmd_shape(year,ny,na)
    integer,intent(in)::year(:)
    integer,intent(out)::ny,na
    integer::i
    ny=1;do i=2,size(year);if(year(i)/=year(i-1))ny=ny+1;end do
    if(mod(size(year),ny)/=0)error stop 'transrate_hmd: inconsistent year blocks'
    na=size(year)/ny
  end subroutine hmd_shape

  subroutine hmd_fix_tail(q)
    real(dp),intent(inout)::q(:)
    integer::i,first
    first=0
    do i=1,size(q)
      if(q(i)>=1.0_dp.or.q(i)>1.0e100_dp)then;first=i;exit;end if
    end do
    if(first>0)q(first:)=0.999_dp
  end subroutine hmd_fix_tail

  subroutine carry_hld_tail(q)
    real(dp),intent(inout)::q(:)
    integer::i,last
    last=0
    do i=1,size(q);if(q(i)>=0.0_dp)last=i;end do
    if(last==0)error stop 'transrate_hld: no qx for sex/year block'
    if(last<size(q))q(last+1:)=q(last)
    if(q(1)<0.0_dp)error stop 'transrate_hld: age must start at zero'
    do i=2,last
      if(q(i)<0.0_dp)q(i)=q(i-1)
    end do
  end subroutine carry_hld_tail

  subroutine unique_years(rec,years)
    type(hld_record),intent(in)::rec(:)
    integer,allocatable,intent(out)::years(:)
    integer::i
    allocate(years(0))
    do i=1,size(rec)
      if(.not.any(years==rec(i)%year1))call append_int(years,rec(i)%year1)
    end do
    call sort_int(years)
  end subroutine unique_years

  integer function race_index(race,file_index) result(ix)
    character(len=*),intent(in)::race(:)
    integer,intent(in)::file_index
    integer::i,j,n,first
    logical::seen
    first=file_index
    do i=1,file_index-1
      if(trim(race(i))==trim(race(file_index)))then;first=i;exit;end if
    end do
    n=0;ix=0
    do i=1,first
      seen=.false.
      do j=1,i-1
        if(trim(race(j))==trim(race(i)))seen=.true.
      end do
      if(.not.seen)then;n=n+1;if(i==first)ix=n;end if
    end do
  end function race_index

  subroutine count_unique_strings(x,n)
    character(len=*),intent(in)::x(:)
    integer,intent(out)::n
    integer::i,j
    logical::seen
    n=0
    do i=1,size(x)
      seen=.false.;do j=1,i-1;if(trim(x(j))==trim(x(i)))seen=.true.;end do
      if(.not.seen)n=n+1
    end do
  end subroutine count_unique_strings

  integer function days_from_civil(y,m,d) result(z)
    integer,intent(in)::y,m,d
    integer::yy,era,yoe,mp,doy,doe
    yy=y;if(m<=2)yy=yy-1
    if(yy>=0)then;era=yy/400;else;era=(yy-399)/400;end if
    yoe=yy-era*400
    if(m>2)then;mp=m-3;else;mp=m+9;end if
    doy=(153*mp+2)/5+d-1;doe=yoe*365+yoe/4-yoe/100+doy
    z=era*146097+doe-719468
  end function days_from_civil

  subroutine header_positions(tok,a,b,ia,ib)
    character(len=*),intent(in)::tok(:),a,b
    integer,intent(out)::ia,ib
    ia=find_token(tok,a);ib=find_token(tok,b)
  end subroutine header_positions

  integer function find_token(tok,key) result(ix)
    character(len=*),intent(in)::tok(:),key
    integer::i
    ix=0;do i=1,size(tok);if(trim(tok(i))==trim(key))then;ix=i;return;end if;end do
  end function find_token

  subroutine strip_dots(tok)
    character(len=*),intent(inout)::tok(:)
    integer::i,j,k
    character(len=len(tok))::tmp
    do i=1,size(tok)
      tmp='';k=0
      do j=1,len_trim(tok(i))
        if(tok(i)(j:j)/='.')then;k=k+1;tmp(k:k)=tok(i)(j:j);end if
      end do
      tok(i)=adjustl(tmp)
    end do
  end subroutine strip_dots

  subroutine split_csv(line,tok)
    character(len=*),intent(in)::line
    character(len=256),allocatable,intent(out)::tok(:)
    integer::i,n,k,start
    n=1;do i=1,len_trim(line);if(line(i:i)==',')n=n+1;end do
    allocate(tok(n));tok='';k=1;start=1
    do i=1,len_trim(line)+1
      if(i>len_trim(line).or.line(i:i)==',')then
        if(i>start)tok(k)=adjustl(line(start:i-1));k=k+1;start=i+1
      end if
    end do
  end subroutine split_csv

  subroutine split_ws(line,tok)
    character(len=*),intent(in)::line
    character(len=128),allocatable,intent(out)::tok(:)
    integer::i,n,k,start,l
    logical::inword
    l=len_trim(line);n=0;inword=.false.
    do i=1,l
      if(line(i:i)/=' '.and.line(i:i)/=char(9))then
        if(.not.inword)n=n+1;inword=.true.
      else;inword=.false.;end if
    end do
    allocate(tok(n));tok='';k=0;i=1
    do while(i<=l)
      do while(i<=l.and.(line(i:i)==' '.or.line(i:i)==char(9)));i=i+1;end do
      if(i>l)exit;start=i
      do while(i<=l.and.line(i:i)/=' '.and.line(i:i)/=char(9));i=i+1;end do
      k=k+1;tok(k)=line(start:i-1)
    end do
  end subroutine split_ws

  subroutine append_int(a,v)
    integer,allocatable,intent(inout)::a(:);integer,intent(in)::v
    integer,allocatable::b(:);allocate(b(size(a)+1));if(size(a)>0)b(:size(a))=a;b(size(b))=v;call move_alloc(b,a)
  end subroutine append_int
  subroutine append_real(a,v)
    real(dp),allocatable,intent(inout)::a(:);real(dp),intent(in)::v
    real(dp),allocatable::b(:);allocate(b(size(a)+1));if(size(a)>0)b(:size(a))=a;b(size(b))=v;call move_alloc(b,a)
  end subroutine append_real
  subroutine append_record(a,v)
    type(hld_record),allocatable,intent(inout)::a(:);type(hld_record),intent(in)::v
    type(hld_record),allocatable::b(:);allocate(b(size(a)+1));if(size(a)>0)b(:size(a))=a;b(size(b))=v;call move_alloc(b,a)
  end subroutine append_record
  subroutine append_records(a,b)
    type(hld_record),allocatable,intent(inout)::a(:)
    type(hld_record),intent(in)::b(:)
    type(hld_record),allocatable::c(:)
    allocate(c(size(a)+size(b)))
    if(size(a)>0)c(:size(a))=a
    if(size(b)>0)c(size(a)+1:)=b
    call move_alloc(c,a)
  end subroutine append_records
  subroutine sort_int(x)
    integer,intent(inout)::x(:)
    integer::i,j,v
    do i=2,size(x)
      v=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=v)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=v
    end do
  end subroutine sort_int
end module relsurv_parsers
