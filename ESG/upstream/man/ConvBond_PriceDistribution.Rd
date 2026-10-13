\name{ConvBond_PriceDistribution}
\alias{ConvBond_PriceDistribution}
\alias{ConvBond_PriceDistribution,Scenarios-method}
\title{ConvBond_PriceDistribution method}
\arguments{
  \item{type}{The name of the asset for which a projection
  has to be proceeded. Can be 'shortRate', 'stock',
  'realEstate', 'liquiditySpread' or 'defaultSpread'. If
  NULL, all assets will be projected.}

  \item{t}{Date of pricing (has to be an integer)}

  \item{T}{Date of maturity for the option}

  \item{nCoupons}{Number of coupons}

  \item{couponsRate}{Rate of coupons}
}
\description{
  Proceed to the projection using the parameters that were
  previously set into the Scenarios objet.
}

