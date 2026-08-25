* Subroutine rtmnG runs a Gibbs recursion for sampling from a multivariate
* normal variable X=(X[1],...,X[d])^T ~ N_d(mean, varcov) truncated by
* the limits (lower < X < upper), where the inequalities hold component-wise.
* Is it assumed that some preliminary quantities are supplied in arrays regr
* and sdc, which defined as  
*   regr: (d, d-1) matrix whose j-th row is varcov[j,] (varcov[-j,-j])^{-1},
*   sdc:  d-vector with the std.dev's the conditional variables (X[j]|X[-j]),
*           sqrt(varcov[j,j] -  varcov[j,] (varcov[-j,-j])^{-1} varcov[,j]) 
* where varcov[j,] denotes the j-th row of varcov,  varcov[j,j] the [j,j]th
* entry of varcov, and varcov[-j,-j] the matrix obtained removing by the j-th 
* row and column. 
*
* Arguments:
*   name    type     size      description
*   n       integer  1         the number of random vectors to be generated
*   d       integer  1         the dimension of X (condition: d>1) 
*   mean    double   d         the mean vector of X
*   regr    double   (d, d-1)  see above
*   sdc     double   d         see above
*   lower   double   d         lower truncation limits
*   upper   double   d         upper truncation limits
*   x       double   (n,d)     (output) matrix of sampled values
*   start   double   d         initial value of the Gibbs iteration
*
* Author: Adelchi Azzalini, Università degli Studi di Padova, 2022
*
******************************************************************************
      subroutine rtmng(n, d, mean, regr, sdc, lower, upper, x, start)
      integer n, d, i, j, k
      double precision x(n, d), mean(d), regr(d,d-1), sdc(d)
      double precision lower(d), upper(d), start(d), meanc, work(d-1) 
      double precision  p1, p2, u, z, qnormr, pnormr, unifrnd
      if(d .lt. 2) return
      call rndstart()
      do i = 1, n                         ! DO loop over the row index
        if(i .eq. 1) then
          do j = 1, d
            x(i,j) = start(j)
          end do
        else
          do j = 1, d
            x(i,j) = x(i-1, j)
          end do    
        endif
        do j = 1, d                       ! DO loop over column index
          do k = 1, j-1                   ! This and next DO loop build 
            work(k) = x(i, k) - mean(k)   ! (x[i,-j] - mean[-j])
          end do
          do k = j+1, d
            work(k-1) = x(i, k) - mean(k)
          end do  
          ! call dblepr("work:", -1, work, d-1)
          meanc = mean(j)                 ! initialize meanc
          do k = 1, d-1                   ! DO loop building meanc 
            meanc = meanc + regr(j, k) * work(k)
          end do
          p1 = pnormr(lower(j), meanc, sdc(j), 1, 0)
          p2 = pnormr(upper(j), meanc, sdc(j), 1, 0)
          u = unifrnd()  
          ! call dblepr1("u", -1, u)
          z = qnormr(p1 + u*(p2-p1), 0.0d0, 1.0d0, 1, 0)
          x(i, j) = meanc + sdc(j) *z
        end do
      end do
      call rndend()
      return
      end 
          
