program test_canonical_v3
  use igraph
  implicit none
  type(graph_t) :: g1,g2,c1,c2
  type(canonical_result_t) :: r1,r2
  type(automorphism_result_t) :: a
  integer :: e1(2,5), e2(2,5), i
  integer, allocatable :: code1(:),code2(:)
  integer :: p(5),invp(5)

  g1=ring_graph(4)
  if(automorphism_count(g1)/=8_i8)error stop 'C4 automorphism count'
  a=automorphisms(g1,max_maps=20)
  if(a%count/=8_i8 .or. a%stored/=8 .or. a%truncated)error stop 'C4 automorphisms'
  if(automorphism_count(full_graph(4))/=24_i8)error stop 'K4 automorphism count'
  if(automorphism_count(g1,[1,2,1,2])/=4_i8)error stop 'colored automorphism count'

  e1(:,1)=[1,2];e1(:,2)=[1,3];e1(:,3)=[2,4];e1(:,4)=[3,4];e1(:,5)=[4,5]
  g1=make_graph(5,e1,.false.)
  p=[3,5,1,4,2]
  do i=1,5;invp(p(i))=i;end do
  do i=1,5
    e2(1,i)=p(e1(1,i));e2(2,i)=p(e1(2,i))
  end do
  g2=make_graph(5,e2,.false.)
  r1=canonical_labeling(g1);r2=canonical_labeling(g2)
  code1=canonical_code(g1);code2=canonical_code(g2)
  if(any(code1/=code2))error stop 'canonical code relabel invariance'
  c1=canonical_form(g1);c2=canonical_form(g2)
  if(any(canonical_code(c1)/=canonical_code(c2)))error stop 'canonical form invariance'
  if(any(r1%inverse(r1%permutation)/=[(i,i=1,5)]))error stop 'canonical inverse'
  if(any(r2%inverse(r2%permutation)/=[(i,i=1,5)]))error stop 'canonical inverse relabeled'

  print *, 'test_canonical_v3: PASS'
end program test_canonical_v3
