\name{rAssetDistribution}
\alias{rAssetDistribution}
\title{rAssetDistribution}
\usage{
  rAssetDistribution(type, t, T, vol, k, ZC,
    nScenarios = NULL, volStock = NULL, stock0 = NULL,
    rho = NULL, volRealEstate = NULL, realEstate0 = NULL,
    eta = NULL, liquiditySpread0 = NULL,
    defaultSpread0 = NULL, volDefault = NULL, alpha = NULL,
    beta = NULL, nCoupons = NULL, couponsRate = NULL,
    omega = NULL, s = NULL, Strike = NULL)
}
\arguments{
  \item{type}{Type of asset. Can be : Zero-Coupon, Bond,
  CBond, ConvBond, EuroCall_S, EuroPut_Stock, EuroCall_ZC,
  EuroPut_ZC or CDS.}

  \item{t}{Date of pricing (has to be an integer)}

  \item{T}{Date of maturity for the option}

  \item{vol}{Volatility for short rates}

  \item{k}{k for rates in vasicek model}

  \item{ZC}{ZC rate input}

  \item{nScenarios}{Number of scenarios}

  \item{volStock}{Volatility for stock}

  \item{stock0}{Initial value for stock}

  \item{rho}{Correlation between stock and short rates}

  \item{volRealEstate}{Volatility for real estate}

  \item{realEstate0}{Initial value for real estate}

  \item{eta}{eta Volatility for LMN model}

  \item{liquiditySpread0}{Initial value for liquidity
  spread}

  \item{defaultSpread0}{Initial value for default spread}

  \item{volDefault}{Volatilty for default spread}

  \item{alpha}{alpha for LMN model}

  \item{beta}{beta Volatility for LMN model}

  \item{nCoupons}{Number of coupons}

  \item{couponsRate}{Rate of coupons}

  \item{omega}{Recoverables in case of default}

  \item{s}{Date of maturity for the underlying}

  \item{Strike}{Strike for options}
}
\description{
  Direct generation for all assets values. Object creation
  is managed internally.
}
\examples{
data(ZC)
rAssetDistribution(type="Zero-Coupon",t=2,T=3,vol=.1, k=2, ZC=ZC, nScenarios=100)
rAssetDistribution(type="Bond",t=3,T=35,nCoupons=20, couponsRate=0.3,vol=.1, k=2,
ZC=ZC, nScenarios=10)
rAssetDistribution(type="CBond",t=5,T=35,nCoupons=5, couponsRate=0.3, omega=5,vol=.1, k=2, ZC=ZC,
nScenarios=10,eta=.05, liquiditySpread0=.01, defaultSpread0=.01, volDefault=.2, alpha=.1, beta=1)
rAssetDistribution(type="EuroPut_Stock",5,25,Strike=98.5,vol=.1,k=2,ZC=ZC,volStock=.2,
stock0=100, rho=.5,nScenarios=10)
rAssetDistribution(type="EuroCall_ZC",4,4.5,s=5, Strike=.985,vol=.1, k=2, ZC=ZC,nScenarios=10)
rAssetDistribution(type="EuroPut_ZC",4,4.5,s=5, Strike=.9385,vol=.1, k=2, ZC=ZC,nScenarios=10)
}

