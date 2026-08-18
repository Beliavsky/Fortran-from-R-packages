! SPDX-License-Identifier: GPL-3.0-or-later
program test_graph
  use degreenet_graph, only : edge_list, reed_molloy
  implicit none
  type(edge_list)::g
  logical::ok
  integer::i,d(4)
  call reed_molloy([2,2,2,2],g,ok)
  if(.not.ok.or.g%nedge/=4)then;print *,'FAIL graph';error stop 1;end if
  d=0
  do i=1,g%nedge;d(g%tail(i))=d(g%tail(i))+1;d(g%head(i))=d(g%head(i))+1;end do
  if(any(d/=[2,2,2,2]))then;print *,'FAIL degrees',d;error stop 1;end if
  print *, 'test_graph: PASS'
end program
