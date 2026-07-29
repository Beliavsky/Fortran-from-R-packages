#'Calculates the estimated value of the firm as the present value of given amount of future Free Cash Flow to the Firm (FCFF) that is discounted at WACC.
#'@description
#'Discounted cash flow (DCF) valuation views the intrinsic value of a security as the present value of its expected future cash flows. When applied to dividends, the DCF model is the discounted dividend approach or dividend discount model (DDM). Free Cash Flow Approach extends DCF analysis to value a firm (company) and its equity securities by valuing free cash flow to the firm (FCFF) and free cash flow to equity (FCFE). Whereas, dividends are the cash flows actually paid to stockholders; however, free cash flows are the cash flows available for distribution to shareholders. Common equity can be valued directly by using FCFE or indirectly by first using a FCFF model to estimate the value of the firm and then subtracting the value of non-common-stock capital (usually the debt) from FCFF to arrive at an estimate of the value of equity. Free cash flow to the firm is the cash flow available to the company’s suppliers of capital after all operating expenses (including taxes) have been paid and necessary investments in working capital (e.g., inventory) and fixed capital (e.g., equipment) have been made. FCFF is the cash flow from operations minus capital expenditures. A  suppliers of capital include common stockholders, bondholders, and sometimes, preferred stockholders. Unlike dividends, FCFF and FCFE are not readily available data. The equations analysts use to calculate FCFF depend on the accounting information available. Analysts need to compute these quantities from available financial information which requires a clear understanding of free cash flows and the ability to interpret and use the information correctly. Forecasting future free cash flows is also a rich and demanding exercise and requires understanding of a corporate financial statements, its operations, investments, and financing (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{firmValueUsingDiscFCFF} is developed to compute the estimated value of the firm as the present value of given amount of FCFF that is discounted at WACC for the values passed to its three arguments. Here, \code{FCFF} is given amount of future Free Cash Flow to the Firm (FCFF) in millions of dollars at time \code{t}. For example a value of 0.04 means 0.4 millions or 400,000 dollars, \code{times} is a vector of number of years ranging from 1 to any specified number of years for which FCFF is to be discounted, and  \code{WACC} is Weighted Average Cost of Capital. Values used for FCFF and the output obtained are in millions of dollars. An output of 1.21494 means 1,214,940 dollars.
#'@param FCFF A vector.
#'@param times A vector.
#'@param WACC A number.
#'@return Input values to three arguments  \code{FCFF}, \code{times}, and  \code{WACC}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'firmValueUsingDiscFCFF(FCFF=c(0.4,0.4,0.4,0.4),times=c(1,2,3,4),WACC=0.12)
#'@export
firmValueUsingDiscFCFF <-function (FCFF,times,WACC)
{
  firmValue <- sum(FCFF/(1 + WACC)^times)
  (firmValue = round(firmValue, digits=6))
}
