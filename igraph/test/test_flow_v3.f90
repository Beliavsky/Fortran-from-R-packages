program test_flow_v3
  use igraph
  implicit none
  type(graph_t) :: g
  type(mincost_result_t) :: mc
  type(gomory_hu_result_t) :: gh
  type(flow_result_t) :: f
  integer :: e(2,5),s,t
  real(dp) :: cap(5),cost(5),treecut

  e(:,1)=[1,2];e(:,2)=[1,3];e(:,3)=[2,3];e(:,4)=[2,4];e(:,5)=[3,4]
  cap=[2.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp]
  cost=[1.0_dp,0.0_dp,0.0_dp,3.0_dp,1.0_dp]
  g=make_graph(4,e,.true.)
  mc=min_cost_flow(g,1,4,3.0_dp,cap,cost)
  if(.not.mc%feasible)error stop 'min cost flow feasibility'
  if(abs(mc%flow-3.0_dp)>1.0e-12_dp .or. abs(mc%cost-7.0_dp)>1.0e-12_dp)then
    error stop 'min cost flow value/cost'
  end if
  if(any(abs(mc%edge_flow-[2.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp])>1.0e-12_dp))then
    error stop 'min cost edge flows'
  end if

  e(:,1)=[1,2];e(:,2)=[2,3];e(:,3)=[3,4];e(:,4)=[4,1];e(:,5)=[1,3]
  g=make_graph(4,e,.false.)
  gh=gomory_hu_tree(g)
  if(gh%tree%m/=3)error stop 'gomory hu tree size'
  do s=1,4
    do t=s+1,4
      f=maxflow(g,s,t)
      treecut=tree_path_min(gh%tree,s,t)
      if(abs(f%value-treecut)>1.0e-12_dp)error stop 'gomory hu pairwise cut mismatch'
    end do
  end do
  print *, 'test_flow_v3: PASS'
contains
  real(dp) function tree_path_min(tree,s,t) result(val)
    type(graph_t),intent(in)::tree
    integer,intent(in)::s,t
    integer::q(tree%n),par(tree%n),pe(tree%n),head,tail,u,v,p
    logical::seen(tree%n)
    val=huge(1.0_dp);seen=.false.;par=0;pe=0;head=1;tail=1;q(1)=s;seen(s)=.true.
    do while(head<=tail)
      u=q(head);head=head+1
      if(u==t)exit
      do p=tree%out_ptr(u),tree%out_ptr(u+1)-1
        v=tree%out_nei(p)
        if(seen(v))cycle
        seen(v)=.true.;par(v)=u;pe(v)=tree%out_eid(p);tail=tail+1;q(tail)=v
      end do
    end do
    u=t
    do while(u/=s)
      if(par(u)==0)error stop 'tree_path_min disconnected'
      val=min(val,tree%weight(pe(u)));u=par(u)
    end do
  end function tree_path_min
end program test_flow_v3
